function poly_checks()
%POLY_CHECKS Lightweight checks for the shared Jacobi polynomial module.

[a, b] = kp.poly.jacobi_recurrence(5, 0, 0);
[al, bl] = kp.poly.legendre_recurrence(5);
x = [-1; -0.5; 0; 0.5; 1];
p = kp.poly.poly_eval(a, b, x, 3);
dp = kp.poly.poly_eval(a, b, x, 3, 1);

assert(isequal(a, al) && isequal(b, bl), 'Legendre recurrence should be the Jacobi(0,0) recurrence.');

assert(all(abs(p(:, 1) - 1 / sqrt(2)) < 1e-12), 'Legendre p0 should be constant 1/sqrt(2).');
assert(all(abs(p(:, 2) - sqrt(3 / 2) * x) < 1e-12), 'Legendre p1 should match the orthonormal linear mode.');
assert(all(abs(dp(:, 2) - sqrt(3 / 2)) < 1e-12), 'Derivative of orthonormal p1 should be constant.');

alpha = [0 0; 1 0; 0 1; 1 1];
recurrence = @(N) kp.poly.jacobi_recurrence(N, 0, 0);
mp = kp.poly.mpoly_eval([0 0; 0.5 -0.25], alpha, recurrence);

assert(size(mp, 1) == 2 && size(mp, 2) == size(alpha, 1), 'Tensor evaluation should return M x N output for one derivative multi-index.');
assert(abs(mp(1, 1) - 1 / 2) < 1e-12, 'Tensor-product constant mode should equal the product of the univariate constants.');

td = kp.poly.total_degree_indices(2, 2);
assert(isequal(td, [0 0; 1 0; 0 1; 2 0; 1 1; 0 2]), 'Total-degree indices should match the reference ordering.');

hc = kp.poly.hyperbolic_cross_indices(2, 3);
assert(ismember([0 0], hc, 'rows') && ismember([1 1], hc, 'rows'), 'Hyperbolic-cross indices should include admissible pairs.');
assert(~ismember([3 3], hc, 'rows'), 'Hyperbolic-cross indices should exclude inadmissible pairs.');

[ac, bc] = kp.poly.chebyshev_recurrence(5);
tc = kp.poly.chebyshev_eval([0; 0.5], [0; 1]);
assert(abs(ac(1)) < 1e-12 && abs(bc(2) - 0.5) < 1e-12, 'Chebyshev recurrence should match the arcsine-measure recurrence.');
assert(size(tc, 1) == 2 && size(tc, 2) == 2, 'Chebyshev tensor evaluation should produce the expected shape.');

r = kp.poly.ratio_eval(a, b, x, 3);
assert(size(r, 1) == numel(x) && size(r, 2) == 3, 'Ratio evaluation should return numel(x) by N values.');

basis = kp.poly.PolynomialBasis.fromTotalDegree(2, 2);
assert(strcmp(basis.Family, 'legendre'), 'PolynomialBasis should default to Legendre.');
basis.fitNormalizationFromPoints([1 0; 0 1; -1 0; 0 -1]);
assert(norm(basis.Center, 2) < 1e-12, 'Symmetric stencil should center at the origin.');
assert(abs(basis.Scale - 1) < 1e-12, 'Unit circle stencil should scale to radius 1.');
vals = basis.evaluate([0 0; 0.5 0]);
assert(size(vals, 1) == 2 && size(vals, 2) == size(basis.IndexSet, 1), 'Basis evaluation should return one column per multi-index.');

basis2 = kp.poly.PolynomialBasis.fromTotalDegree(2, 1, 'Center', [2 -1], 'Scale', 2);
xn = basis2.normalizePoints([4 -1; 2 1]);
assert(isequal(xn, [1 0; 0 1]), 'Stored center and scale should normalize points consistently.');
dv = basis2.evaluate([4 -1], [1 0]);
expectedDx = (sqrt(3 / 2) * (1 / sqrt(2))) / 2;
assert(abs(dv(1, 2) - expectedDx) < 1e-12, 'Derivative scaling should account for the stored radius factor and the untouched tensor-product factor.');

disp('poly checks passed');
end
