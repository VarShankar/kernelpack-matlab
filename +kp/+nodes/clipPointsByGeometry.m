function [Xkeep, keepMask, phi] = clipPointsByGeometry(X, geometry, varargin)
%CLIPPOINTSBYGEOMETRY Keep points inside or outside a geometry level set.

    parser = inputParser();
    parser.addRequired('X', @(x) validateattributes(x, {'numeric'}, {'2d', 'real', 'finite'}));
    parser.addRequired('geometry');
    parser.addParameter('Keep', 'inside', @(x) any(strcmpi(x, {'inside', 'outside'})));
    parser.addParameter('Tolerance', 0, @(x) validateattributes(x, {'numeric'}, {'scalar', 'real', 'finite'}));
    parser.addParameter('BoundaryClearance', 0, @(x) validateattributes(x, {'numeric'}, {'scalar', 'real', 'finite', 'nonnegative'}));
    parser.addParameter('MinSignedDistance', -inf, @(x) validateattributes(x, {'numeric'}, {'scalar', 'real', 'finite'}));
    parser.addParameter('MaxSignedDistance', inf, @(x) validateattributes(x, {'numeric'}, {'scalar', 'real', 'finite'}));
    parser.addParameter('AutoBuildLevelSet', true, @(x) islogical(x) || isnumeric(x));
    parser.addParameter('UseParallel', true, @(x) islogical(x) || isnumeric(x));
    parser.addParameter('ChunkSize', 5000, @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', '>=', 1}));
    parser.addParameter('MinParallelPoints', 20000, @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', '>=', 1}));
    parser.parse(X, geometry, varargin{:});
    opts = parser.Results;

    levelSet = ensureGeometryLevelSet(geometry, logical(opts.AutoBuildLevelSet));
    phi = evaluateLevelSetValues(levelSet, X, opts);

    switch lower(opts.Keep)
        case 'inside'
            keepMask = phi >= (opts.BoundaryClearance - opts.Tolerance);
        case 'outside'
            keepMask = phi <= (opts.Tolerance - opts.BoundaryClearance);
        otherwise
            error('kp:nodes:BadKeepMode', 'Keep must be ''inside'' or ''outside''.');
    end

    keepMask = keepMask & (phi >= opts.MinSignedDistance) & (phi <= opts.MaxSignedDistance);

    Xkeep = X(keepMask, :);
end

function phi = evaluateLevelSetValues(levelSet, X, opts)
    nPts = size(X, 1);
    useParallel = logical(opts.UseParallel) && canUseParfor() && nPts >= opts.MinParallelPoints;
    if ~useParallel
        phi = levelSet.Evaluate(X);
        return;
    end

    nChunks = ceil(nPts / opts.ChunkSize);
    chunks = cell(nChunks, 1);
    for k = 1:nChunks
        i0 = (k - 1) * opts.ChunkSize + 1;
        i1 = min(k * opts.ChunkSize, nPts);
        chunks{k} = X(i0:i1, :);
    end

    phiChunks = cell(nChunks, 1);
    model = levelSet.getEvaluationModel();
    parfor k = 1:nChunks
        phiChunks{k} = kp.geometry.RBFLevelSet.evaluateModel(model, chunks{k});
    end
    phi = vertcat(phiChunks{:});
end

function tf = canUseParfor()
    tf = license('test', 'Distrib_Computing_Toolbox');
end

function levelSet = ensureGeometryLevelSet(geometry, autoBuild)
    if ~ismethod(geometry, 'getLevelSet')
        error('kp:nodes:MissingLevelSet', 'Geometry object must provide getLevelSet().');
    end

    levelSet = geometry.getLevelSet();
    needsBuild = isempty(levelSet) || ~isa(levelSet, 'kp.geometry.RBFLevelSet') || levelSet.n == 0;
    if needsBuild
        if ~autoBuild
            error('kp:nodes:MissingBuiltLevelSet', 'Geometry level set is not built.');
        end
        if ismethod(geometry, 'buildLevelSetFromGeometricModel')
            geometry.buildLevelSetFromGeometricModel([]);
        elseif ismethod(geometry, 'buildLevelSet')
            geometry.buildLevelSet();
        else
            error('kp:nodes:CannotBuildLevelSet', 'Geometry object does not expose a recognized level-set builder.');
        end
        levelSet = geometry.getLevelSet();
    end
end
