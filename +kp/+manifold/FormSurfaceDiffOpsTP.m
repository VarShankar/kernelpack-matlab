function [L, Gx, Gy, Gz] = FormSurfaceDiffOpsTP( ...
        Xd, rbf, drbfor, d2rbf, nr, tree, n, ell, neighborIds)
%FORMSURFACEDIFFOPSTP Tangent-plane Laplace-Beltrami and gradient assembly.

Np = size(Xd, 1);
row_index = zeros(n, Np);
col_index = zeros(n, Np);
wghts_lap = zeros(n, Np);
wghts_gx = zeros(n, Np);
wghts_gy = zeros(n, Np);
wghts_gz = zeros(n, Np);
tarray = ones(1, n);

recurrence = @(N) kp.poly.jacobi_recurrence(N, 0, 0);
polyIndices = kp.poly.total_degree_indices(2, ell);
if nargin < 9 || isempty(neighborIds)
    IDX = knnsearch(tree, Xd, 'k', n);
else
    IDX = neighborIds;
    if ~isequal(size(IDX), [Np, n])
        error('kp:manifold:BadTangentPlaneNeighbors', ...
            'NeighborIds must contain one %d-node stencil per surface node.', n);
    end
end

parfor i = 1:Np
    idx = IDX(i, :);
    x = Xd(idx, :); %#ok<PFBNS>
    centerNormal = nr(idx(1), :).'; %#ok<PFBNS>
    [A, B, geom] = kp.manifold.detail.tangentPlaneLocalSystem( ...
        x, centerNormal, rbf, drbfor, d2rbf, polyIndices, recurrence);
    W = kp.manifold.detail.solveLocalAugmentedSystem(A, B);
    grad = [W(1:n, 2:3), zeros(n, 1)] * [geom.R, centerNormal].';

    wghts_lap(:, i) = W(1:n, 1);
    wghts_gx(:, i) = grad(1:n, 1);
    wghts_gy(:, i) = grad(1:n, 2);
    wghts_gz(:, i) = grad(1:n, 3);
    row_index(:, i) = i .* tarray;
    col_index(:, i) = idx;
end

L = sparse(row_index(:), col_index(:), wghts_lap(:), Np, Np);
Gx = sparse(row_index(:), col_index(:), wghts_gx(:), Np, Np);
Gy = sparse(row_index(:), col_index(:), wghts_gy(:), Np, Np);
Gz = sparse(row_index(:), col_index(:), wghts_gz(:), Np, Np);
end
