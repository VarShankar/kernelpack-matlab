function p = chebyshev_eval(x, alpha, d)
%CHEBYSHEV_EVAL Evaluate tensor-product orthonormal Chebyshev polynomials.

[m, dim] = size(x);
if nargin < 3
    d = zeros(1, dim);
else
    assert(size(d, 2) == dim, 'Inputs x and d must have the same number of columns');
end

n = size(alpha, 1);
pcount = size(d, 1);
assert(size(alpha, 2) == dim, 'Inputs x and alpha must have the same number of columns');

p = ones(m, n, pcount);
[a, b] = kp.poly.chebyshev_recurrence(max(alpha(:)) + 1);

for qd = 1:pcount
    for qdim = 1:dim
        temp = kp.poly.poly_eval(a, b, x(:, qdim), max(alpha(:, qdim)), d(qd, qdim));
        p(:, :, qd) = p(:, :, qd) .* temp(:, alpha(:, qdim) + 1);
    end
end
end
