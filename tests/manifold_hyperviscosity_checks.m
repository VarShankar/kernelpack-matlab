function manifold_hyperviscosity_checks()
%MANIFOLD_HYPERVISCOSITY_CHECKS Verify the surface hyperviscosity scaling formula.

N = 256;
k = 3;
xi = 4;
X = kp.geometry.fibonacciSphere(N);
X = X ./ vecnorm(X, 2, 2);
nr = X;
op = kp.manifold.rbffdop(2, xi, 2, 0);
tree = KDTreeSearcher(X);
[L, Gx, Gy, Gz] = kp.manifold.FormSurfaceDiffOpsTP( ...
    X, op.rbf, op.drbfor, op.d2rbf, nr, tree, op.stencilSize, op.ell);
h = sqrt(1 / N);

[gamma, info] = kp.manifold.hyperviscosityCoefficient(L, Gx, Gy, Gz, X, nr, k, h);
gammaExpected = 3^(-k) * (-1)^(1 - k) * info.hypterm / info.etaMean;

assert(isfinite(gamma), 'Hyperviscosity coefficient is not finite.');
assert(abs(gamma - gammaExpected) <= 100 * eps(max(1, abs(gammaExpected))), ...
    'Hyperviscosity coefficient does not match the paper/source formula.');
assert(numel(info.gammaComponents) == 3 && ...
        abs(sum(info.gammaComponents) - gamma) <= 100 * eps(max(1, abs(gamma))), ...
    'Coordinate hyperviscosity coefficients do not recover the scalar aggregate.');
assert(info.k == k && abs(info.h - h) <= eps(h), ...
    'Hyperviscosity diagnostic metadata does not match the requested parameters.');

[gammaAuto, infoAuto] = kp.manifold.hyperviscosityCoefficient( ...
    L, Gx, Gy, Gz, X, nr, NaN, h, xi);
expectedPower = ceil((xi + max(infoAuto.q)) / 2);
assert(infoAuto.automaticPower && infoAuto.k == expectedPower, ...
    'Automatic hyperviscosity power does not match the q-based rule.');
assert(2 * infoAuto.k - max(infoAuto.q) >= xi, ...
    'Automatic hyperviscosity power does not preserve the target order.');
assert(isfinite(gammaAuto), ...
    'Automatic hyperviscosity coefficient is not finite.');

disp('manifold hyperviscosity checks passed');
fprintf('  gamma: %.6e\n', gamma);
fprintf('  q: [%.3f, %.3f, %.3f]\n', info.q(1), info.q(2), info.q(3));
fprintf('  tau: [%.3e, %.3e, %.3e]\n', info.tau(1), info.tau(2), info.tau(3));
fprintf('  automatic k for xi=%d: %d\n', xi, infoAuto.k);
end
