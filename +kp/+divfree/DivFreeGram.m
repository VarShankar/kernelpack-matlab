function A = DivFreeGram(dfrbf, X, Y, r)
%DIVFREEGRAM Build the divergence-free PHS Gram matrix.

if size(X, 2) ~= size(Y, 2)
    error('kp:divfree:DimensionMismatch', 'X and Y must have the same spatial dimension.');
end

N = size(X, 1);
M = size(Y, 1);
dim = size(X, 2);

switch dim
    case 2
        [xj, xk] = ndgrid(X(:, 1), Y(:, 1));
        [yj, yk] = ndgrid(X(:, 2), Y(:, 2));
        A = zeros(2 * N, 2 * M);
        A(1:N, 1:M) = dfrbf{1, 1}(r, xj, xk);
        A(1:N, M + 1:2 * M) = dfrbf{1, 2}(r, xj, xk, yj, yk);
        A(N + 1:2 * N, 1:M) = dfrbf{2, 1}(r, xj, xk, yj, yk);
        A(N + 1:2 * N, M + 1:2 * M) = dfrbf{2, 2}(r, yj, yk);
    case 3
        [xj, xk] = ndgrid(X(:, 1), Y(:, 1));
        [yj, yk] = ndgrid(X(:, 2), Y(:, 2));
        [zj, zk] = ndgrid(X(:, 3), Y(:, 3));
        A = zeros(3 * N, 3 * M);

        A(1:N, 1:M) = dfrbf{1, 1}(r, xj, xk);
        A(1:N, M + 1:2 * M) = dfrbf{1, 2}(r, xj, xk, yj, yk);
        A(1:N, 2 * M + 1:3 * M) = dfrbf{1, 3}(r, xj, xk, zj, zk);

        A(N + 1:2 * N, 1:M) = dfrbf{2, 1}(r, xj, xk, yj, yk);
        A(N + 1:2 * N, M + 1:2 * M) = dfrbf{2, 2}(r, yj, yk);
        A(N + 1:2 * N, 2 * M + 1:3 * M) = dfrbf{2, 3}(r, yj, yk, zj, zk);

        A(2 * N + 1:3 * N, 1:M) = dfrbf{3, 1}(r, xj, xk, zj, zk);
        A(2 * N + 1:3 * N, M + 1:2 * M) = dfrbf{3, 2}(r, yj, yk, zj, zk);
        A(2 * N + 1:3 * N, 2 * M + 1:3 * M) = dfrbf{3, 3}(r, zj, zk);
    otherwise
        error('kp:divfree:BadDimension', 'Only 2D and 3D are supported.');
end
end
