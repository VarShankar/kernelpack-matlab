function poly_checks()
%POLY_CHECKS Lightweight checks for the shared Jacobi polynomial module.

[a, b] = kp.poly.jacobi_recurrence(5, 0, 0);
x = [-1; -0.5; 0; 0.5; 1];
p = kp.poly.poly_eval(a, b, x, 3);
dp = kp.poly.poly_eval(a, b, x, 3, 1);

assert(all(abs(p(:, 1) - 1 / sqrt(2)) < 1e-12), 'Legendre p0 should be constant 1/sqrt(2).');
assert(all(abs(p(:, 2) - sqrt(3 / 2) * x) < 1e-12), 'Legendre p1 should match the orthonormal linear mode.');
assert(all(abs(dp(:, 2) - sqrt(3 / 2)) < 1e-12), 'Derivative of orthonormal p1 should be constant.');

alpha = [0 0; 1 0; 0 1; 1 1];
recurrence = @(N) kp.poly.jacobi_recurrence(N, 0, 0);
mp = kp.poly.mpoly_eval([0 0; 0.5 -0.25], alpha, recurrence);

assert(size(mp, 1) == 2 && size(mp, 2) == size(alpha, 1), 'Tensor evaluation should return M x N output for one derivative multi-index.');
assert(abs(mp(1, 1) - 1 / 2) < 1e-12, 'Tensor-product constant mode should equal the product of the univariate constants.');

disp('poly checks passed');
end
