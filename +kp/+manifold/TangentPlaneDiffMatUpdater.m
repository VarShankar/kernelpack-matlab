classdef TangentPlaneDiffMatUpdater < handle
    %TANGENTPLANEDIFFMATUPDATER Update tangent-plane RBF-FD matrices.
    %   This updater follows the moving-domain differentiation-matrix idea:
    %   unchanged stencils avoid fresh local factorizations. On manifolds,
    %   unchanged stencils still have changed local tangent-plane matrices,
    %   so the updater computes the new weights by defect correction using
    %   the last direct local factorization as the approximate inverse.

    properties
        xi (1,1) double = 4
        theta (1,1) double = 2
        StencilSize (1,1) double = NaN
        DefectTolerance (1,1) double = 1.0e-4
        MaxDefectIterations (1,1) double = 4
        UseParallel (1,1) logical = true
        NeighborUpdateMode (1,1) string = "periodic"
        NeighborSearchInterval (1,1) double = 5
        AssembliesSinceNeighborSearch (1,1) double = 0
        op kp.manifold.rbffdop = kp.manifold.rbffdop()
        NeighborIds double = zeros(0, 0)
        NeighborMargins double = zeros(0, 1)
        PreviousX double = zeros(0, 3)
        FactorKind uint8 = zeros(0, 1)
        FactorInv double = zeros(0, 0, 0)
        QrFactors cell = {}
        LastStats struct = struct( ...
            'numStencils', 0, ...
            'direct', 0, ...
            'defectCorrected', 0, ...
            'defectFailedRefactored', 0, ...
            'meanDefectIterations', NaN, ...
            'maxDefectIterations', 0, ...
            'maxRelativeResidual', 0, ...
            'neighborSearches', 0, ...
            'neighborReuses', 0, ...
            'neighborReuseRejected', 0)
    end

    methods
        function obj = TangentPlaneDiffMatUpdater(xi, varargin)
            if nargin >= 1 && ~isempty(xi)
                obj.xi = xi;
            end

            parser = inputParser();
            parser.addParameter('Theta', obj.theta);
            parser.addParameter('StencilSize', obj.StencilSize);
            parser.addParameter('DefectTolerance', obj.DefectTolerance);
            parser.addParameter('MaxDefectIterations', obj.MaxDefectIterations);
            parser.addParameter('UseParallel', obj.UseParallel);
            parser.addParameter('NeighborUpdateMode', obj.NeighborUpdateMode);
            parser.addParameter('NeighborSearchInterval', obj.NeighborSearchInterval);
            parser.parse(varargin{:});

            obj.theta = parser.Results.Theta;
            obj.StencilSize = parser.Results.StencilSize;
            obj.DefectTolerance = parser.Results.DefectTolerance;
            obj.MaxDefectIterations = parser.Results.MaxDefectIterations;
            obj.UseParallel = parser.Results.UseParallel;
            obj.NeighborUpdateMode = lower(string(parser.Results.NeighborUpdateMode));
            obj.NeighborSearchInterval = parser.Results.NeighborSearchInterval;
            obj.op = kp.manifold.rbffdop(2, obj.xi, obj.theta, 0);
            if isfinite(obj.StencilSize) && obj.StencilSize > 0
                obj.op.stencilSize = max(obj.op.stencilSize, round(obj.StencilSize));
            end
        end

        function reset(obj)
            obj.NeighborIds = zeros(0, 0);
            obj.NeighborMargins = zeros(0, 1);
            obj.PreviousX = zeros(0, 3);
            obj.AssembliesSinceNeighborSearch = 0;
            obj.FactorKind = zeros(0, 1, 'uint8');
            obj.FactorInv = zeros(0, 0, 0);
            obj.QrFactors = {};
            obj.LastStats = defaultStats();
        end

        function [L, Gx, Gy, Gz, stats] = assemble(obj, X, normals, neighborCoordinates)
            validateInputs(X, normals);
            if nargin < 4 || isempty(neighborCoordinates)
                neighborCoordinates = X;
            end
            if size(neighborCoordinates, 1) ~= size(X, 1)
                error('kp:manifold:BadNeighborCoordinates', ...
                    'Neighbor coordinates must contain one row per surface node.');
            end
            Np = size(X, 1);
            n = obj.op.stencilSize;
            m = n + obj.op.polyM;

            oldNeighbors = obj.NeighborIds;
            oldKind = obj.FactorKind;
            oldInv = obj.FactorInv;
            oldQrFactors = obj.QrFactors;
            hasCache = isequal(size(oldNeighbors), [Np, n]) && numel(oldKind) == Np && ...
                isequal(size(oldInv), [m, m, Np]);
            [IDX, neighborMargins, neighborStats] = neighborIndices( ...
                obj, neighborCoordinates, n, oldNeighbors, hasCache);

            row_index = zeros(n, Np);
            col_index = zeros(n, Np);
            wghts_lap = zeros(n, Np);
            wghts_gx = zeros(n, Np);
            wghts_gy = zeros(n, Np);
            wghts_gz = zeros(n, Np);
            tarray = ones(1, n);
            factorUpdates = cell(Np, 1);
            modes = zeros(Np, 1);
            iters = zeros(Np, 1);
            relres = zeros(Np, 1);

            rbf = obj.op.rbf;
            drbfor = obj.op.drbfor;
            d2rbf = obj.op.d2rbf;
            ell = obj.op.ell;
            defectTol = obj.DefectTolerance;
            maxIters = obj.MaxDefectIterations;
            recurrence = @(N) kp.poly.jacobi_recurrence(N, 0, 0);
            polyIndices = kp.poly.total_degree_indices(2, ell);

            if obj.UseParallel && hasCache
                parfor i = 1:Np
                    oldNeighbor = oldNeighbors(i, :);
                    oldKindI = oldKind(i);
                    oldInvI = [];
                    if oldKindI > 0
                        oldInvI = oldInv(:, :, i);
                    end
                    [wlap, wgrad, factorUpdate, mode, iter, rr, idx] = assembleOne( ...
                        i, X, normals, IDX, oldNeighbor, oldKindI, oldInvI, true, ...
                        rbf, drbfor, d2rbf, polyIndices, recurrence, n, defectTol, maxIters);
                    wghts_lap(:, i) = wlap;
                    wghts_gx(:, i) = wgrad(:, 1);
                    wghts_gy(:, i) = wgrad(:, 2);
                    wghts_gz(:, i) = wgrad(:, 3);
                    row_index(:, i) = i .* tarray;
                    col_index(:, i) = idx;
                    factorUpdates{i} = factorUpdate;
                    modes(i) = mode;
                    iters(i) = iter;
                    relres(i) = rr;
                end
            elseif obj.UseParallel
                parfor i = 1:Np
                    [wlap, wgrad, factorUpdate, mode, iter, rr, idx] = assembleOne( ...
                        i, X, normals, IDX, [], uint8(0), [], false, ...
                        rbf, drbfor, d2rbf, polyIndices, recurrence, n, defectTol, maxIters);
                    wghts_lap(:, i) = wlap;
                    wghts_gx(:, i) = wgrad(:, 1);
                    wghts_gy(:, i) = wgrad(:, 2);
                    wghts_gz(:, i) = wgrad(:, 3);
                    row_index(:, i) = i .* tarray;
                    col_index(:, i) = idx;
                    factorUpdates{i} = factorUpdate;
                    modes(i) = mode;
                    iters(i) = iter;
                    relres(i) = rr;
                end
            elseif hasCache
                for i = 1:Np
                    oldNeighbor = oldNeighbors(i, :);
                    oldKindI = oldKind(i);
                    oldInvI = [];
                    if oldKindI > 0
                        oldInvI = oldInv(:, :, i);
                    end
                    [wlap, wgrad, factorUpdate, mode, iter, rr, idx] = assembleOne( ...
                        i, X, normals, IDX, oldNeighbor, oldKindI, oldInvI, true, ...
                        rbf, drbfor, d2rbf, polyIndices, recurrence, n, defectTol, maxIters);
                    wghts_lap(:, i) = wlap;
                    wghts_gx(:, i) = wgrad(:, 1);
                    wghts_gy(:, i) = wgrad(:, 2);
                    wghts_gz(:, i) = wgrad(:, 3);
                    row_index(:, i) = i .* tarray;
                    col_index(:, i) = idx;
                    factorUpdates{i} = factorUpdate;
                    modes(i) = mode;
                    iters(i) = iter;
                    relres(i) = rr;
                end
            else
                for i = 1:Np
                    [wlap, wgrad, factorUpdate, mode, iter, rr, idx] = assembleOne( ...
                        i, X, normals, IDX, [], uint8(0), [], false, ...
                        rbf, drbfor, d2rbf, polyIndices, recurrence, n, defectTol, maxIters);
                    wghts_lap(:, i) = wlap;
                    wghts_gx(:, i) = wgrad(:, 1);
                    wghts_gy(:, i) = wgrad(:, 2);
                    wghts_gz(:, i) = wgrad(:, 3);
                    row_index(:, i) = i .* tarray;
                    col_index(:, i) = idx;
                    factorUpdates{i} = factorUpdate;
                    modes(i) = mode;
                    iters(i) = iter;
                    relres(i) = rr;
                end
            end

            L = sparse(row_index(:), col_index(:), wghts_lap(:), Np, Np);
            Gx = sparse(row_index(:), col_index(:), wghts_gx(:), Np, Np);
            Gy = sparse(row_index(:), col_index(:), wghts_gy(:), Np, Np);
            Gz = sparse(row_index(:), col_index(:), wghts_gz(:), Np, Np);

            obj.NeighborIds = IDX;
            obj.NeighborMargins = neighborMargins;
            obj.PreviousX = neighborCoordinates;
            if neighborStats.searches > 0
                obj.AssembliesSinceNeighborSearch = 0;
            else
                obj.AssembliesSinceNeighborSearch = obj.AssembliesSinceNeighborSearch + 1;
            end
            [obj.FactorKind, obj.FactorInv, obj.QrFactors] = ...
                mergeFactorUpdates(oldKind, oldInv, oldQrFactors, ...
                factorUpdates, modes, hasCache, m, Np);
            stats = summarizeStats(modes, iters, relres);
            stats.neighborSearches = neighborStats.searches;
            stats.neighborReuses = neighborStats.reuses;
            stats.neighborReuseRejected = neighborStats.rejected;
            obj.LastStats = stats;
        end
    end
end

function [IDX, margins, stats] = neighborIndices(obj, neighborCoordinates, n, oldNeighbors, hasCache)
stats = struct('searches', 0, 'reuses', 0, 'rejected', 0);
Np = size(neighborCoordinates, 1);

if hasCache && shouldReuseNeighbors(obj, neighborCoordinates)
    IDX = oldNeighbors;
    margins = obj.NeighborMargins;
    stats.reuses = 1;
    return;
end

if hasCache
    stats.rejected = 1;
end
stats.searches = 1;
kSearch = min(n + 1, Np);
tree = KDTreeSearcher(neighborCoordinates);
[idxSearch, distSearch] = knnsearch(tree, neighborCoordinates, 'k', kSearch);
IDX = idxSearch(:, 1:n);
if kSearch >= n + 1
    margins = distSearch(:, n + 1) - distSearch(:, n);
else
    margins = inf(Np, 1);
end
end

function reuse = shouldReuseNeighbors(obj, neighborCoordinates)
mode = lower(string(obj.NeighborUpdateMode));
reuse = false;
switch mode
    case "always"
        return;
    case "reuse"
        reuse = ~isempty(obj.NeighborIds);
    case "periodic"
        reuse = ~isempty(obj.NeighborIds) && ...
            obj.AssembliesSinceNeighborSearch < obj.NeighborSearchInterval;
    case "motionbound"
        if isempty(obj.PreviousX) || isempty(obj.NeighborMargins) || ...
                size(obj.PreviousX, 1) ~= size(neighborCoordinates, 1)
            return;
        end
        displacement = vecnorm(neighborCoordinates - obj.PreviousX, 2, 2);
        maxDisplacement = max(displacement);
        % Any pairwise distance can change by at most twice the maximum
        % nodal displacement. Preserve neighbor ordering with a conservative
        % extra factor of two on the previous n/(n+1) distance margin.
        reuse = 4 * maxDisplacement < min(obj.NeighborMargins);
    otherwise
        error('kp:manifold:BadNeighborUpdateMode', ...
            'Unknown neighbor update mode "%s".', obj.NeighborUpdateMode);
end
end

function [wlap, wgrad, factor, mode, iter, rr, idx] = assembleOne( ...
    i, X, normals, IDX, oldNeighbor, oldKindI, oldInvI, hasCache, ...
    rbf, drbfor, d2rbf, polyIndices, recurrence, n, defectTol, maxIters)

idx = IDX(i, :);
centerNormal = normals(idx(1), :).';
[A, B, geom] = kp.manifold.detail.tangentPlaneLocalSystem(X(idx, :), centerNormal, ...
    rbf, drbfor, d2rbf, polyIndices, recurrence);

canDefect = hasCache && oldKindI > 0 && isequal(oldNeighbor, idx);
if canDefect
    [W, ok, iter, rr] = defectCorrectInverse(A, B, oldInvI, defectTol, maxIters);
    if ok
        factor = [];
        mode = 2;
    else
        [W, factor, rr] = kp.manifold.detail.solveAndFactorLocalAugmentedSystem(A, B);
        mode = 3;
    end
else
    [W, factor, rr] = kp.manifold.detail.solveAndFactorLocalAugmentedSystem(A, B);
    mode = 1;
    iter = 0;
end

wlap = W(1:n, 1);
gradPlane = W(1:n, 2:3);
wgrad = [gradPlane, zeros(n, 1)] * [geom.R, centerNormal].';
end

function [factorKind, factorInv, qrFactors] = mergeFactorUpdates( ...
    oldKind, oldInv, oldQrFactors, factorUpdates, modes, hasCache, m, Np)
if hasCache
    factorKind = oldKind(:);
    factorInv = oldInv;
    qrFactors = oldQrFactors(:);
else
    factorKind = zeros(Np, 1, 'uint8');
    factorInv = zeros(m, m, Np);
    qrFactors = cell(Np, 1);
end

if all(modes == 2)
    return;
end

for i = 1:numel(factorUpdates)
    if modes(i) ~= 2
        factor = factorUpdates{i};
        if isempty(factor)
            continue;
        end
        switch string(factor.kind)
            case "lu"
                factorKind(i) = 1;
                factorInv(:, :, i) = inverseFromFactor(factor, m);
                qrFactors{i} = [];
            case "qr"
                factorKind(i) = 2;
                factorInv(:, :, i) = inverseFromFactor(factor, m);
                qrFactors{i} = factor;
            otherwise
                error('kp:manifold:BadLocalFactor', ...
                    'Unknown local factorization kind "%s".', factor.kind);
        end
    end
end
end

function invA = inverseFromFactor(factor, m)
I = eye(m);
switch string(factor.kind)
    case "lu"
        invA = factor.U \ (factor.L \ I(factor.P, :));
    case "qr"
        invA = kp.manifold.detail.solveLocalFactor(factor, I);
    otherwise
        error('kp:manifold:BadLocalFactor', ...
            'Unknown local factorization kind "%s".', factor.kind);
end
end

function [W, ok, iter, rr] = defectCorrectInverse(A, B, invA, defectTol, maxIters)
W = invA * B;
rr = norm(B - A * W, 'fro') / max(norm(B, 'fro'), eps);
ok = rr <= defectTol;
iter = 0;

while ~ok && iter < maxIters
    residual = B - A * W;
    W = W + invA * residual;
    iter = iter + 1;
    rr = norm(B - A * W, 'fro') / max(norm(B, 'fro'), eps);
    ok = rr <= defectTol;
end
end

function stats = summarizeStats(modes, iters, relres)
stats = defaultStats();
stats.numStencils = numel(modes);
stats.direct = nnz(modes == 1);
stats.defectCorrected = nnz(modes == 2);
stats.defectFailedRefactored = nnz(modes == 3);
stats.meanDefectIterations = mean(iters(modes == 2), 'omitnan');
stats.maxDefectIterations = max([0; iters(:)]);
stats.maxRelativeResidual = max([0; relres(:)]);
end

function stats = defaultStats()
stats = struct( ...
    'numStencils', 0, ...
    'direct', 0, ...
    'defectCorrected', 0, ...
    'defectFailedRefactored', 0, ...
        'meanDefectIterations', NaN, ...
        'maxDefectIterations', 0, ...
        'maxRelativeResidual', 0, ...
        'neighborSearches', 0, ...
        'neighborReuses', 0, ...
        'neighborReuseRejected', 0);
end

function validateInputs(X, normals)
if size(X, 2) ~= 3 || size(normals, 2) ~= 3 || size(X, 1) ~= size(normals, 1)
    error('kp:manifold:BadUpdaterInputs', ...
        'Expected matching N-by-3 node and normal arrays.');
end
end
