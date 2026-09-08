function results = moving_surface_adr_tp_convergence_breathing_sphere(xiVals, Nvals, dtScale)
%MOVING_SURFACE_ADR_TP_CONVERGENCE_BREATHING_SPHERE Cleaner moving-sphere ADR study.

if nargin < 1 || isempty(xiVals)
    xiVals = [2, 4, 6];
end
if nargin < 2 || isempty(Nvals)
    Nvals = [256, 576, 1024, 1600];
end
if nargin < 3 || isempty(dtScale)
    dtScale = 0.5;
end

mu = 0.1;
T = 0.12;
hyppow = 3;

Nvals = Nvals(:).';
xiVals = xiVals(:).';
numLevels = numel(Nvals);
numXi = numel(xiVals);
results = repmat(struct( ...
    'xi', 0, ...
    'N', zeros(1, numLevels), ...
    'h', zeros(1, numLevels), ...
    'dt', zeros(1, numLevels), ...
    'nsteps', zeros(1, numLevels), ...
    'relerr', zeros(1, numLevels), ...
    'rate', nan(1, numLevels)), 1, numXi);

for ix = 1:numXi
    xi = xiVals(ix);
    results(ix).xi = xi;
    for level = 1:numLevels
        N = Nvals(level);
        h = sqrt(4 * pi / N);
        dtTarget = dtScale * h^(xi / 3);
        nsteps = max(3, ceil(T / dtTarget));
        dt = T / nsteps;

        fprintf('Breathing sphere ADR TP convergence\n');
        fprintf('  xi: %d\n', xi);
        fprintf('  level %d / %d\n', level, numLevels);
        fprintf('  N: %d\n', N);
        fprintf('  h: %.6e\n', h);
        fprintf('  dt: %.6e\n', dt);
        fprintf('  steps: %d\n', nsteps);

        err = runOneLevel(N, xi, mu, T, dt, nsteps, hyppow);
        results(ix).N(level) = N;
        results(ix).h(level) = h;
        results(ix).dt(level) = dt;
        results(ix).nsteps(level) = nsteps;
        results(ix).relerr(level) = err;
        fprintf('  relative L2 error: %.6e\n', err);
    end

    for k = 2:numLevels
        results(ix).rate(k) = log(results(ix).relerr(k - 1) / results(ix).relerr(k)) / ...
            log(results(ix).h(k - 1) / results(ix).h(k));
    end
end

fig = figure('Color', 'w', 'Position', [100, 100, 760, 560]);
hold on;
colors = lines(numXi);
for ix = 1:numXi
    plot(results(ix).h, results(ix).relerr, '-o', 'LineWidth', 1.6, ...
        'MarkerSize', 7, 'Color', colors(ix, :), ...
        'DisplayName', sprintf('\\xi = %d', results(ix).xi));
end
hold off;
set(gca, 'XScale', 'log', 'YScale', 'log');
grid on;
xlabel('h');
ylabel('Relative L2 error');
title('Breathing sphere ADR tangent-plane convergence');
legend('Location', 'best');
set(gca, 'XDir', 'reverse');

imgPath = fullfile(pwd, 'docs', 'figures', 'moving_surface_adr_tp_convergence_breathing_sphere.png');
exportgraphics(fig, imgPath, 'Resolution', 180);
save(fullfile(pwd, 'moving_surface_adr_tp_convergence_breathing_sphere_results.mat'), 'results');
end

function relerr = runOneLevel(N, xi, mu, T, dt, nsteps, hyppow)
op = kp.manifold.rbffdop(2, xi, 2, 0);

[x0, ~, U] = breathingSphereGeometry(N, 0.0);
c0 = exactField(0.0, U);

[x1, nr1] = breathingSphereGeometry(N, dt);
v1 = (x1 - x0) / dt;
[L, Gx, Gy, Gz, H, hyp_gamma] = surfaceOpsAndHyperviscosity(x1, nr1, op, hyppow);
I = speye(N, N);
lhs = I - mu * dt * L;
rhs = c0 + dt * exactForcing(dt, U, mu) ...
    - dt * c0 .* (Gx * v1(:, 1) + Gy * v1(:, 2) + Gz * v1(:, 3)) ...
    + dt * hyp_gamma * c0 .* (H * (v1(:, 1) + v1(:, 2) + v1(:, 3)));
c1 = lhs \ rhs;

[x2, nr2] = breathingSphereGeometry(N, 2 * dt);
v2 = (3 * x2 - 4 * x1 + x0) / (2 * dt);
[L, Gx, Gy, Gz, H, hyp_gamma] = surfaceOpsAndHyperviscosity(x2, nr2, op, hyppow);
lhs = I - (2 / 3) * mu * dt * L;
c2a = 2 * c1 - c0;
rhs = (4 / 3) * c1 - (1 / 3) * c0 + (2 / 3) * dt * ...
    exactForcing(2 * dt, U, mu) ...
    - (2 / 3) * dt * c2a .* (Gx * v2(:, 1) + Gy * v2(:, 2) + Gz * v2(:, 3)) ...
    + (2 / 3) * dt * hyp_gamma * c2a .* (H * (v2(:, 1) + v2(:, 2) + v2(:, 3)));
c2 = lhs \ rhs;
c3 = c2;

for step = 3:nsteps
    tnow = step * dt;
    [x3, nr3] = breathingSphereGeometry(N, tnow);
    v3 = (11 * x3 - 18 * x2 + 9 * x1 - 2 * x0) / (6 * dt);
    [L, Gx, Gy, Gz, H, hyp_gamma] = surfaceOpsAndHyperviscosity(x3, nr3, op, hyppow);
    lhs = I - (6 / 11) * mu * dt * L;

    c3a = 3 * c2 - 3 * c1 + c0;
    rhs = (18 / 11) * c2 - (9 / 11) * c1 + (2 / 11) * c0 ...
        + (6 / 11) * dt * exactForcing(step * dt, U, mu) ...
        - (6 / 11) * dt * c3a .* (Gx * v3(:, 1) + Gy * v3(:, 2) + Gz * v3(:, 3)) ...
        + (6 / 11) * dt * hyp_gamma * c3a .* (H * (v3(:, 1) + v3(:, 2) + v3(:, 3)));
    c3 = lhs \ rhs;

    c0 = c1; c1 = c2; c2 = c3;
    x0 = x1; x1 = x2; x2 = x3;
end

cex = exactField(T, U);
relerr = norm(c3 - cex) / max(norm(cex), 1.0e-14);
end

function c = exactField(t, U)
b = baseField(U);
c = exp(-t) .* b;
end

function f = exactForcing(t, U, mu)
b = baseField(U);
lapS = surfaceLaplacianBase(U);
R = sphereRadius(t);
Rt = sphereRadiusPrime(t);
divu = 2 * Rt / R;
c = exp(-t) .* b;
ct = -c;
lapc = exp(-t) .* lapS ./ (R^2);
f = ct + c .* divu - mu .* lapc;
end

function b = baseField(U)
ux = U(:, 1);
uy = U(:, 2);
b = 0.8 + 0.15 * ux + 0.05 * (ux.^2 - uy.^2);
end

function lapb = surfaceLaplacianBase(U)
ux = U(:, 1);
uy = U(:, 2);
lapb = -0.30 * ux - 0.30 * (ux.^2 - uy.^2);
end

function [X, nr, U] = breathingSphereGeometry(N, t)
U = kp.geometry.fibonacciSphere(N);
U = U ./ vecnorm(U, 2, 2);
R = sphereRadius(t);
X = R .* U;
nr = U;
end

function R = sphereRadius(t)
R = 1 + 0.2 * sin(t);
end

function Rp = sphereRadiusPrime(t)
Rp = 0.2 * cos(t);
end

function [L, Gx, Gy, Gz, H, hyp_gamma] = surfaceOpsAndHyperviscosity(X, nr, op, hyppow)
tree = KDTreeSearcher(X);
[L, Gx, Gy, Gz] = kp.manifold.FormSurfaceDiffOpsTP( ...
    X, op.rbf, op.drbfor, op.d2rbf, nr, tree, op.stencilSize, op.ell);
assertNegativeLargestReal(L);

H = L^hyppow;
hx = sqrt(1 / size(X, 1));
hyp_gamma = kp.manifold.hyperviscosityCoefficient(H, Gx, Gy, Gz, X, nr, hyppow, hx);
end

function assertNegativeLargestReal(L)
tol = 1.0e-2;
[~, sig] = tryLargestAlgebraic(L, tol);
if sig > tol
    error('kp:manifold:BadSpectrum', ...
        ['+ve real eigenvalue detected, smooth normals with more nearest ' ...
         'neighbors or check points.']);
end
end

function [V, D] = tryLargestAlgebraic(A, tol)
warnState = warning();
cleanup = onCleanup(@() warning(warnState));
warning('off', 'all');
[V, D] = eigs(A, 1, 'largestreal', 'Tolerance', tol);
end
