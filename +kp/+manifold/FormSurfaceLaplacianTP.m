function [L_rbf, Np, interp_struct] = FormSurfaceLaplacianTP(Xd, rbf, drbfor, d2rbf, nr, tree, n, ell)
%FORMSURFACELAPLACIANTP Tangent-plane Laplace-Beltrami RBF-FD assembly.

Np = size(Xd, 1);
row_index = zeros(n, Np);
col_index = zeros(n, Np);
wghts_lap = zeros(n, Np);
tarray = ones(1, n);

alph = 0;
bet = 0;
recurrence = @(N) kp.poly.jacobi_recurrence(N, alph, bet);
a = kp.poly.total_degree_indices(2, ell);
interp_struct = cell(Np + 1, 6);
IDX = knnsearch(tree, Xd, 'k', n);

parfor i = 1:Np
    idx = IDX(i, :);
    x = Xd(idx, :);
    geom = kp.manifold.detail.tangentPlaneLocalGeometry(x, nr(idx(1), :).');

    P = kp.poly.mpoly_eval(geom.xc, a, recurrence);
    Lp = (kp.poly.mpoly_eval(geom.xc(1, :), a, recurrence, [2, 0]) + ...
          kp.poly.mpoly_eval(geom.xc(1, :), a, recurrence, [0, 2])) ./ (geom.w2^2);

    Ar = rbf(1, geom.rd);
    Lr = d2rbf(1, geom.rd(:, 1)) + drbfor(1, geom.rd(:, 1));
    A = [[Ar, P]; [P.', zeros(size(P, 2))]];
    dA = decomposition(A, 'auto');
    lit = dA \ [Lr; Lp.'];

    wghts_lap(:, i) = lit(1:n, 1);
    row_index(:, i) = i .* tarray;
    col_index(:, i) = idx;

    interp_struct(i, :) = {x(1, :), idx, geom.xw, dA, geom.w2, geom.R};
end

L_rbf = sparse(row_index(:), col_index(:), wghts_lap(:), Np, Np);
X_stencil = cell2mat(interp_struct(1:Np, 1));
interp_struct{Np + 1, 1} = KDTreeSearcher(X_stencil);
end
