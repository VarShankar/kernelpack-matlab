function [X, info] = generatePoissonNodesInBox(radiusOrFunc, x_min, x_max, varargin)
%GENERATEPOISSONNODESINBOX KernelPack-style Poisson disk sampling on a box.
%
%   X = kp.nodes.generatePoissonNodesInBox(radius, x_min, x_max)
%   uses a fixed-radius sampler.
%
%   X = kp.nodes.generatePoissonNodesInBox(radFunc, x_min, x_max, 'MinRadius', hmin)
%   uses a variable-density sampler driven by radFunc(p, hmin).
%
%   Both modes support KernelPack-style outer refinement with:
%     'BoundaryPoints', Xb
%     'BoundaryRefinementFraction', frac
%     'BoundaryDistance', dist

    % Parse once up front so every strip sees the same deterministic
    % sampling rules and refinement settings.
    opts = parseInputs(radiusOrFunc, x_min, x_max, varargin{:});

    if any(opts.x_max <= opts.x_min)
        X = zeros(0, numel(opts.x_min));
        info = emptyInfo(opts);
        return;
    end

    % KernelPack's current deterministic path collapses to one canonical
    % strip, while nondeterministic runs may split across strips for
    % parallel speed.
    stripBoxes = buildStripBoxes(opts.x_min, opts.x_max, opts.split_tol, opts.strip_count);
    localClouds = cell(opts.strip_count, 1);
    stripSeeds = uint32(mod(double(opts.base_seed) + 104729 * (0:opts.strip_count-1), 2^32 - 1));

    if opts.use_parallel
        parfor k = 1:opts.strip_count
            localClouds{k} = poissonStripSample(stripBoxes{k}.sample_min, stripBoxes{k}.sample_max, ...
                opts, stripSeeds(k));
        end
    else
        for k = 1:opts.strip_count
            localClouds{k} = poissonStripSample(stripBoxes{k}.sample_min, stripBoxes{k}.sample_max, ...
                opts, stripSeeds(k));
        end
    end

    % Concatenate the strip clouds in canonical order, then clip them back
    % to the requested axis-aligned box.
    X = flattenStripClouds(localClouds, opts.x_min, opts.x_max);
    info = struct( ...
        'dimension', numel(opts.x_min), ...
        'mode', opts.mode, ...
        'radius', opts.radius, ...
        'min_radius', opts.min_radius, ...
        'attempts', opts.attempts, ...
        'seed', opts.seed, ...
        'deterministic', opts.deterministic, ...
        'strip_count', opts.strip_count, ...
        'used_parallel', opts.use_parallel, ...
        'boundary_refinement_fraction', opts.boundary_refinement_fraction, ...
        'boundary_distance', opts.boundary_distance, ...
        'num_points', size(X, 1));
end

function opts = parseInputs(radiusOrFunc, x_min, x_max, varargin)
    parser = inputParser();
    parser.addRequired('radiusOrFunc');
    parser.addRequired('x_min', @(x) validateattributes(x, {'numeric'}, {'vector', 'real', 'finite'}));
    parser.addRequired('x_max', @(x) validateattributes(x, {'numeric'}, {'vector', 'real', 'finite'}));
    parser.addParameter('Attempts', 30, @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', '>=', 1}));
    parser.addParameter('Seed', [], @(x) isempty(x) || (isscalar(x) && isnumeric(x) && isfinite(x)));
    parser.addParameter('Deterministic', [], @(x) isempty(x) || islogical(x) || isnumeric(x));
    parser.addParameter('StripCount', [], @(x) isempty(x) || (isscalar(x) && isnumeric(x) && isfinite(x) && x >= 1));
    parser.addParameter('UseParallel', true, @(x) islogical(x) || isnumeric(x));
    parser.addParameter('MinRadius', [], @(x) isempty(x) || (isscalar(x) && isnumeric(x) && isfinite(x) && x > 0));
    parser.addParameter('BoundaryPoints', zeros(0, 0), @(x) isnumeric(x) && ismatrix(x));
    parser.addParameter('BoundaryRefinementFraction', 1.0, @(x) validateattributes(x, {'numeric'}, {'scalar', 'real', 'finite', '>', 0, '<=', 1}));
    parser.addParameter('BoundaryDistance', 0, @(x) validateattributes(x, {'numeric'}, {'scalar', 'real', 'finite', '>=', 0}));
    parser.parse(radiusOrFunc, x_min, x_max, varargin{:});

    opts = parser.Results;
    opts.x_min = opts.x_min(:).';
    opts.x_max = opts.x_max(:).';
    if numel(opts.x_min) ~= numel(opts.x_max)
        error('kp:nodes:BoxSizeMismatch', 'x_min and x_max must have the same length.');
    end

    if isempty(opts.Seed)
        opts.seed = double(randi(2^31 - 1));
    else
        opts.seed = double(opts.Seed);
    end
    opts.base_seed = uint32(mod(opts.seed, 2^32));

    if isempty(opts.Deterministic)
        opts.deterministic = ~isempty(opts.Seed);
    else
        opts.deterministic = logical(opts.Deterministic);
    end

    requestedStripCount = [];
    if ~isempty(opts.StripCount)
        requestedStripCount = max(1, floor(opts.StripCount));
    end
    opts.strip_count = defaultStripCount(logical(opts.UseParallel), opts.deterministic, requestedStripCount);
    opts.use_parallel = logical(opts.UseParallel) && opts.strip_count > 1 && canUseParfor();
    opts.attempts = opts.Attempts;

    opts.boundary_points = opts.BoundaryPoints;
    if ~isempty(opts.boundary_points) && size(opts.boundary_points, 2) ~= numel(opts.x_min)
        error('kp:nodes:BoundaryPointDimensionMismatch', ...
            'BoundaryPoints must have the same column count as the box dimension.');
    end
    opts.boundary_refinement_fraction = opts.BoundaryRefinementFraction;
    opts.boundary_distance = opts.BoundaryDistance;
    opts.has_boundary_refinement = ~isempty(opts.boundary_points) && size(opts.boundary_points, 1) > 0 && ...
        opts.boundary_refinement_fraction < 1 && opts.boundary_distance > 0;

    % Normalize both the fixed-radius and variable-density interfaces into
    % one internal option struct that the sampler uses everywhere.
    if isa(radiusOrFunc, 'function_handle')
        if isempty(opts.MinRadius)
            error('kp:nodes:MissingMinRadius', ...
                'Variable-density sampling requires a ''MinRadius'' value.');
        end
        opts.rad_func = radiusOrFunc;
        opts.radius = NaN;
        opts.min_radius = opts.MinRadius;
        if opts.has_boundary_refinement
            opts.mode = 'variable_radius_with_boundary_refinement';
            opts.grid_radius = opts.boundary_refinement_fraction * opts.min_radius;
        else
            opts.mode = 'variable_radius';
            opts.grid_radius = opts.min_radius;
        end
    else
        validateattributes(radiusOrFunc, {'numeric'}, {'scalar', 'real', 'finite', 'positive'});
        opts.rad_func = [];
        opts.radius = double(radiusOrFunc);
        opts.min_radius = opts.radius;
        if opts.has_boundary_refinement
            opts.mode = 'fixed_radius_with_boundary_refinement';
            opts.grid_radius = opts.boundary_refinement_fraction * opts.radius;
        else
            opts.mode = 'fixed_radius';
            opts.grid_radius = opts.radius;
        end
    end

    opts.split_tol = opts.min_radius;
end

function info = emptyInfo(opts)
    info = struct( ...
        'dimension', numel(opts.x_min), ...
        'mode', opts.mode, ...
        'radius', opts.radius, ...
        'min_radius', opts.min_radius, ...
        'attempts', opts.attempts, ...
        'seed', opts.seed, ...
        'deterministic', opts.deterministic, ...
        'strip_count', opts.strip_count, ...
        'used_parallel', false, ...
        'boundary_refinement_fraction', opts.boundary_refinement_fraction, ...
        'boundary_distance', opts.boundary_distance, ...
        'num_points', 0);
end

function stripCount = defaultStripCount(useParallel, deterministic, requestedStripCount)
    if deterministic
        stripCount = 1;
        return;
    end
    if ~isempty(requestedStripCount)
        stripCount = requestedStripCount;
        return;
    end
    if ~useParallel || ~canUseParfor()
        stripCount = 1;
        return;
    end
    pool = gcp('nocreate');
    if isempty(pool)
        stripCount = 1;
    else
        stripCount = max(1, pool.NumWorkers);
    end
end

function tf = canUseParfor()
    tf = license('test', 'Distrib_Computing_Toolbox');
end

function stripBoxes = buildStripBoxes(x_min, x_max, split_tol, stripCount)
    % The box is split only along the first coordinate, following the
    % strip-decomposition shape used in the C++ node generator.
    dim0 = (x_max(1) - x_min(1)) / stripCount;
    stripBoxes = cell(stripCount, 1);
    for k = 1:stripCount
        sample_min = x_min;
        sample_max = x_max;
        sample_min(1) = x_min(1) + (k - 1) * dim0;
        sample_max(1) = sample_min(1) + dim0 - 0.33 * split_tol;
        sample_max(1) = min(x_max(1), sample_max(1));
        stripBoxes{k} = struct('sample_min', sample_min, 'sample_max', sample_max);
    end
end

function X = poissonStripSample(sample_min, sample_max, opts, seed)
    dim = numel(sample_min);
    if any(sample_max <= sample_min)
        X = zeros(0, dim);
        return;
    end

    cellSize = opts.grid_radius / sqrt(dim);
    gridSize = max(ones(1, dim), ceil((sample_max - sample_min) / cellSize));
    grid = containers.Map('KeyType', 'char', 'ValueType', 'uint32');
    points = zeros(128, dim);
    active = zeros(128, 1, 'uint32');

    stream = RandStream('mt19937ar', 'Seed', double(seed));
    x0 = sample_min + rand(stream, 1, dim) .* (sample_max - sample_min);
    [points, active, grid, nPoints, nActive] = addPoint(points, active, grid, x0, sample_min, cellSize, 0, 0);

    % This is Bridson-style active-list sampling, but with the active
    % radius allowed to vary from point to point.
    while nActive > 0
        pick = randi(stream, nActive);
        activeIdx = active(pick);
        base = points(activeIdx, :);
        activeRadius = localRadius(base, opts);
        accepted = false;
        for attempt = 1:opts.attempts
            candidate = proposeCandidate(base, activeRadius, dim, stream);
            if any(candidate < sample_min) || any(candidate > sample_max)
                continue;
            end
            if hasConflictingNeighbor(candidate, activeIdx, points, nPoints, grid, sample_min, cellSize, gridSize, opts)
                continue;
            end
            [points, active, grid, nPoints, nActive] = addPoint(points, active, grid, candidate, sample_min, cellSize, nPoints, nActive);
            accepted = true;
            break;
        end
        if ~accepted
            active(pick) = active(nActive);
            nActive = nActive - 1;
        end
    end

    X = points(1:nPoints, :);
end

function radius = localRadius(point, opts)
    % The active radius may come from a fixed h, a user-supplied radius
    % function, or the near-boundary refinement rule.
    switch opts.mode
        case 'fixed_radius'
            radius = opts.radius;
        case 'fixed_radius_with_boundary_refinement'
            radius = boundaryRefinedRadius(point, opts);
        case 'variable_radius'
            radius = opts.rad_func(point, opts.min_radius);
            if radius <= 1e-10
                radius = opts.min_radius;
            end
        case 'variable_radius_with_boundary_refinement'
            baseRadius = opts.rad_func(point, opts.min_radius);
            if baseRadius <= 1e-10
                baseRadius = opts.min_radius;
            end
            if boundaryRadFrac(point, opts) < 1
                radius = opts.boundary_refinement_fraction * opts.min_radius;
            else
                radius = baseRadius;
            end
        otherwise
            error('kp:nodes:UnknownMode', 'Unknown Poisson sampling mode.');
    end
end

function frac = boundaryRadFrac(point, opts)
    if ~opts.has_boundary_refinement
        frac = 1.0;
        return;
    end
    dist = nearestBoundaryDistance(point, opts.boundary_points);
    if dist <= opts.boundary_distance
        frac = opts.boundary_refinement_fraction;
    else
        frac = 1.0;
    end
end

function radius = boundaryRefinedRadius(point, opts)
    radius = boundaryRadFrac(point, opts) * opts.radius;
end

function tf = hasConflictingNeighbor(point, activeIdx, points, nPoints, grid, x_min, cellSize, gridSize, opts)
    candidateRadius = localRadius(point, opts);
    switch opts.mode
        case 'fixed_radius_with_boundary_refinement'
            radiusForReach = max(candidateRadius, opts.radius);
            excludeActive = false;
            pairwise = true;
        case {'fixed_radius', 'variable_radius', 'variable_radius_with_boundary_refinement'}
            radiusForReach = candidateRadius;
            excludeActive = true;
            pairwise = false;
        otherwise
            error('kp:nodes:UnknownMode', 'Unknown Poisson sampling mode.');
    end

    % The search radius is conservative enough to catch any neighbor that
    % could violate the active-radius acceptance rule.
    idx = pointToCell(point, x_min, cellSize);
    reach = max(1, ceil(radiusForReach / cellSize));
    ranges = cell(1, numel(idx));
    for d = 1:numel(idx)
        ranges{d} = max(1, idx(d) - reach):min(gridSize(d), idx(d) + reach);
    end

    tf = false;
    neighborCells = enumerateNeighborCells(ranges);
    for k = 1:size(neighborCells, 1)
        key = cellKey(neighborCells(k, :));
        if ~isKey(grid, key)
            continue;
        end
        j = double(grid(key));
        if j < 1 || j > nPoints
            continue;
        end
        if excludeActive && j == activeIdx
            continue;
        end
        if pairwise
            threshold = max(candidateRadius, localRadius(points(j, :), opts));
        else
            threshold = candidateRadius;
        end
        if norm(point - points(j, :), 2) < threshold
            tf = true;
            return;
        end
    end
end

function candidate = proposeCandidate(base, radius, dim, stream)
    direction = randn(stream, 1, dim);
    direction = direction ./ max(norm(direction, 2), eps);
    shellRadius = radius * (1 + rand(stream, 1) * (2^dim - 1))^(1 / dim);
    candidate = base + shellRadius * direction;
end

function [points, active, grid, nPoints, nActive] = addPoint(points, active, grid, point, x_min, cellSize, nPoints, nActive)
    if nPoints == 0
        nPoints = 1;
    else
        nPoints = nPoints + 1;
    end
    if nPoints > size(points, 1)
        points = [points; zeros(size(points, 1), size(points, 2))];
    end
    points(nPoints, :) = point;

    nActive = nActive + 1;
    if nActive > numel(active)
        active = [active; zeros(numel(active), 1, 'uint32')];
    end
    active(nActive) = uint32(nPoints);

    idx = pointToCell(point, x_min, cellSize);
    grid(cellKey(idx)) = uint32(nPoints);
end

function idx = pointToCell(point, x_min, cellSize)
    idx = floor((point - x_min) ./ cellSize) + 1;
    idx = max(idx, 1);
end

function key = cellKey(idx)
    key = sprintf('%d_', idx);
end

function cells = enumerateNeighborCells(ranges)
    dim = numel(ranges);
    grids = cell(1, dim);
    [grids{:}] = ndgrid(ranges{:});
    n = numel(grids{1});
    cells = zeros(n, dim);
    for d = 1:dim
        cells(:, d) = grids{d}(:);
    end
end

function X = flattenStripClouds(localClouds, x_min, x_max)
    dim = numel(x_min);
    totalCount = sum(cellfun(@(x) size(x, 1), localClouds));
    X = zeros(totalCount, dim);
    row0 = 1;
    for k = 1:numel(localClouds)
        pts = localClouds{k};
        m = size(pts, 1);
        if m == 0
            continue;
        end
        rows = row0:(row0 + m - 1);
        X(rows, :) = pts;
        row0 = row0 + m;
    end
    X = X(1:row0-1, :);
    if ~isempty(X)
        inBox = all(X >= x_min, 2) & all(X <= x_max, 2);
        X = X(inBox, :);
    end
end

function dist = nearestBoundaryDistance(point, boundaryPoints)
    if isempty(boundaryPoints)
        dist = inf;
        return;
    end
    diffs = boundaryPoints - point;
    dist = sqrt(min(sum(diffs.^2, 2), [], 'omitnan'));
end
