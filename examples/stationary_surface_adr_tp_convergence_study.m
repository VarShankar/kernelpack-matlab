function results = stationary_surface_adr_tp_convergence_study(xiVals, Nvals, dtScale)
%STATIONARY_SURFACE_ADR_TP_CONVERGENCE_STUDY Manufactured ADR convergence on the unit sphere.

if nargin < 1 || isempty(xiVals)
    xiVals = [2, 4, 6];
end
if nargin < 2 || isempty(Nvals)
    Nvals = [256, 576, 1024, 1600];
end
if nargin < 3 || isempty(dtScale)
    dtScale = 0.1;
end

mu = 0.1;
T = 0.12;

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

        fprintf('Stationary surface ADR TP convergence\n');
        fprintf('  xi: %d\n', xi);
        fprintf('  level %d / %d\n', level, numLevels);
        fprintf('  N: %d\n', N);
        fprintf('  h: %.6e\n', h);
        fprintf('  dt: %.6e\n', dt);
        fprintf('  steps: %d\n', nsteps);

        err = runOneLevel(N, xi, mu, T, dt, nsteps);
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
title('Stationary surface ADR tangent-plane convergence');
legend('Location', 'best');
set(gca, 'XDir', 'reverse');

imgPath = fullfile(pwd, 'docs', 'figures', 'stationary_surface_adr_tp_convergence.png');
exportgraphics(fig, imgPath, 'Resolution', 180);
save(fullfile(pwd, 'stationary_surface_adr_tp_convergence_results.mat'), 'results');
end

function relerr = runOneLevel(N, xi, mu, T, dt, nsteps)
op = kp.manifold.rbffdop(2, xi, 2, 0);
X = kp.geometry.fibonacciSphere(N);
X = X ./ vecnorm(X, 2, 2);
nr = X;
tree = KDTreeSearcher(X);
[L, ~, ~, ~] = kp.manifold.FormSurfaceDiffOpsTP( ...
    X, op.rbf, op.drbfor, op.d2rbf, nr, tree, op.stencilSize, op.ell);
assertNegativeLargestReal(L);

u0 = exactField(0.0, X);
I = speye(N, N);
lhs = I - mu * dt * L;
u1 = lhs \ (u0 + dt * exactForcing(dt, X, mu));

lhs = I - (2 / 3) * mu * dt * L;
u2 = lhs \ ((4 / 3) * u1 - (1 / 3) * u0 + (2 / 3) * dt * exactForcing(2 * dt, X, mu));
u3 = u2;

for step = 3:nsteps
    tnow = step * dt;
    lhs = I - (6 / 11) * mu * dt * L;
    rhs = (18 / 11) * u2 - (9 / 11) * u1 + (2 / 11) * u0 + ...
        (6 / 11) * dt * exactForcing(tnow, X, mu);
    u3 = lhs \ rhs;
    u0 = u1; u1 = u2; u2 = u3;
end

uex = exactField(T, X);
relerr = norm(u3 - uex) / max(norm(uex), 1.0e-14);
end

function u = exactField(t, X)
u = exp(-t) .* baseField(X);
end

function f = exactForcing(t, X, mu)
b = baseField(X);
lapb = surfaceLaplacianBase(X);
u = exp(-t) .* b;
ut = -u;
f = ut - mu * exp(-t) .* lapb;
end

function b = baseField(X)
x = X(:, 1);
y = X(:, 2);
b = 0.8 + 0.15 * x + 0.05 * (x.^2 - y.^2);
end

function lapb = surfaceLaplacianBase(X)
x = X(:, 1);
y = X(:, 2);
lapb = -0.30 * x - 0.30 * (x.^2 - y.^2);
end

function assertNegativeLargestReal(L)
tol = 1.0e-2;
warnState = warning();
cleanup = onCleanup(@() warning(warnState));
warning('off', 'all');
[~, sig] = eigs(L, 1, 'largestreal', 'Tolerance', tol);
if sig > tol
    error('kp:manifold:BadSpectrum', ...
        ['+ve real eigenvalue detected, smooth normals with more nearest ' ...
         'neighbors or check points.']);
end
end
