function [P, poly] = df_poly_basis_from_jacobi(p, p_out, recurrence, opts)
%DF_POLY_BASIS_FROM_JACOBI Divergence-free polynomial basis on a stencil.
%   This ports the older Jacobi/Legendre-based divergence-free polynomial
%   construction into the kernelpack-matlab namespace.

if nargin < 3 || isempty(recurrence)
    recurrence = @(N) kp.poly.jacobi_recurrence(N, 0, 0);
end
if nargin < 4
    opts = struct();
end
if ~isfield(opts, 'center') || isempty(opts.center)
    opts.center = mean(p, 1);
end
if ~isfield(opts, 'scale') || isempty(opts.scale)
    X0 = p - opts.center;
    opts.scale = max(vecnorm(X0.', 2, 1));
    if ~(opts.scale > 0)
        opts.scale = 1;
    end
end
if ~isfield(opts, 'use_pca_frame') || isempty(opts.use_pca_frame)
    opts.use_pca_frame = false;
end
if ~isfield(opts, 'rank_tol') || isempty(opts.rank_tol)
    opts.rank_tol = 1e-12;
end

[m, d] = size(p);
if ~(d == 2 || d == 3)
    error('kp:divfree:BadDimension', 'p must be m x 2 or m x 3.');
end

X = p - opts.center;
R = eye(d);
if opts.use_pca_frame
    C = (X.' * X) / max(m, 1);
    [U, ~, ~] = svd(C);
    R = U;
    for j = 1:d
        [~, k] = max(abs(R(:, j)));
        if R(k, j) < 0
            R(:, j) = -R(:, j);
        end
    end
    if det(R) < 0
        R(:, end) = -R(:, end);
    end
end
xi = (X * R) / opts.scale;

pA = p_out + 1;
aA = kp.poly.total_degree_indices(d, pA);
K = size(aA, 1);

if d == 2
    dPhidx = kp.poly.mpoly_eval(xi, aA, recurrence, [1, 0]);
    dPhidy = kp.poly.mpoly_eval(xi, aA, recurrence, [0, 1]);
    V = zeros(2 * m, K);
    V(1:m, :) = dPhidy;
    V(m + 1:2 * m, :) = -dPhidx;
else
    dPhidx = kp.poly.mpoly_eval(xi, aA, recurrence, [1, 0, 0]);
    dPhidy = kp.poly.mpoly_eval(xi, aA, recurrence, [0, 1, 0]);
    dPhidz = kp.poly.mpoly_eval(xi, aA, recurrence, [0, 0, 1]);
    V = zeros(3 * m, 3 * K);
    V(1:m, 1:K) = 0;
    V(m + 1:2 * m, 1:K) = dPhidz;
    V(2 * m + 1:3 * m, 1:K) = -dPhidy;

    V(1:m, K + 1:2 * K) = -dPhidz;
    V(m + 1:2 * m, K + 1:2 * K) = 0;
    V(2 * m + 1:3 * m, K + 1:2 * K) = dPhidx;

    V(1:m, 2 * K + 1:3 * K) = dPhidy;
    V(m + 1:2 * m, 2 * K + 1:3 * K) = -dPhidx;
    V(2 * m + 1:3 * m, 2 * K + 1:3 * K) = 0;
end

[Qref, Rq, piv] = qr(V, 0);
dd = abs(diag(Rq));
if isempty(dd) || dd(1) == 0
    r = 0;
else
    r = find(dd > opts.rank_tol * dd(1), 1, 'last');
    if isempty(r)
        r = 0;
    end
end
Qref = Qref(:, 1:r);
R11 = Rq(1:r, 1:r);
piv = piv(1:r);

P = apply_vec_rotation(Qref, R);
[Qfull, ~] = qr(P, 0);
Nnull = Qfull(:, r + 1:end);

poly = struct();
poly.d = d;
poly.p_out = p_out;
poly.pA = pA;
poly.aA = aA;
poly.center = opts.center;
poly.scale = opts.scale;
poly.R = R;
poly.piv = piv;
poly.R11 = R11;
poly.Nnull = Nnull;
poly.eval = @(xq) df_poly_eval_same_basis(xq, recurrence, poly);
end

function Pq = df_poly_eval_same_basis(xq, recurrence, poly)
[d, aA, R, scale, center, piv, R11] = deal(poly.d, poly.aA, poly.R, poly.scale, poly.center, poly.piv, poly.R11);
nq = size(xq, 1);
if size(xq, 2) ~= d
    error('kp:divfree:BadEvalDimension', 'xq must have %d columns.', d);
end

xiq = ((xq - center) * R) / scale;
K = size(aA, 1);

if d == 2
    dPhidx = kp.poly.mpoly_eval(xiq, aA, recurrence, [1, 0]);
    dPhidy = kp.poly.mpoly_eval(xiq, aA, recurrence, [0, 1]);
    Vq = zeros(2 * nq, K);
    Vq(1:nq, :) = dPhidy;
    Vq(nq + 1:2 * nq, :) = -dPhidx;
else
    dPhidx = kp.poly.mpoly_eval(xiq, aA, recurrence, [1, 0, 0]);
    dPhidy = kp.poly.mpoly_eval(xiq, aA, recurrence, [0, 1, 0]);
    dPhidz = kp.poly.mpoly_eval(xiq, aA, recurrence, [0, 0, 1]);
    Vq = zeros(3 * nq, 3 * K);

    Vq(1:nq, 1:K) = 0;
    Vq(nq + 1:2 * nq, 1:K) = dPhidz;
    Vq(2 * nq + 1:3 * nq, 1:K) = -dPhidy;

    Vq(1:nq, K + 1:2 * K) = -dPhidz;
    Vq(nq + 1:2 * nq, K + 1:2 * K) = 0;
    Vq(2 * nq + 1:3 * nq, K + 1:2 * K) = dPhidx;

    Vq(1:nq, 2 * K + 1:3 * K) = dPhidy;
    Vq(nq + 1:2 * nq, 2 * K + 1:3 * K) = -dPhidx;
    Vq(2 * nq + 1:3 * nq, 2 * K + 1:3 * K) = 0;
end

Qref_q = Vq(:, piv) / R11;
Pq = apply_vec_rotation(Qref_q, R);
end

function U = apply_vec_rotation(Uref, R)
d = size(R, 1);
n3 = size(Uref, 1);
if mod(n3, d) ~= 0
    error('kp:divfree:BadVectorStack', 'Stacked vector field length must be divisible by the dimension.');
end
n = n3 / d;

U = zeros(size(Uref));
blocks = cell(1, d);
for j = 1:d
    blocks{j} = Uref((j - 1) * n + 1:j * n, :);
end
for i = 1:d
    Ui = 0;
    for j = 1:d
        Ui = Ui + R(i, j) * blocks{j};
    end
    U((i - 1) * n + 1:i * n, :) = Ui;
end
end
