function results = mean_curvature_flow_sphere_example(xiVals, Nvals)
%MEAN_CURVATURE_FLOW_SPHERE_EXAMPLE Evolve a sphere by X_t = Delta_Gamma X.

if nargin < 1 || isempty(xiVals)
    xiVals = [2, 4, 6];
end
if nargin < 2 || isempty(Nvals)
    Nvals = [256, 576, 1024];
end

R0 = 1.4;
T = 0.04;
dtScale = 0.02;

xiVals = xiVals(:).';
Nvals = Nvals(:).';
numXi = numel(xiVals);
numLevels = numel(Nvals);
results = repmat(struct( ...
    'xi', 0, ...
    'N', zeros(1, numLevels), ...
    'h', zeros(1, numLevels), ...
    'dt', zeros(1, numLevels), ...
    'nsteps', zeros(1, numLevels), ...
    'radiusError', zeros(1, numLevels), ...
    'shapeError', zeros(1, numLevels), ...
    'rate', nan(1, numLevels)), 1, numXi);

for ix = 1:numXi
    xi = xiVals(ix);
    results(ix).xi = xi;
    for level = 1:numLevels
        N = Nvals(level);
        U = kp.geometry.fibonacciSphere(N);
        U = U ./ vecnorm(U, 2, 2);
        X0 = R0 .* U;
        h = R0 * sqrt(4 * pi / N);
        dtTarget = dtScale * h^2;
        nsteps = max(1, ceil(T / dtTarget));
        dt = T / nsteps;

        fprintf('Mean-curvature-flow sphere example\n');
        fprintf('  xi: %d\n', xi);
        fprintf('  N: %d\n', N);
        fprintf('  h: %.6e\n', h);
        fprintf('  dt: %.6e\n', dt);
        fprintf('  steps: %d\n', nsteps);

        solver = kp.manifold.MeanCurvatureFlowSolver( ...
            'Mode', "semiimplicit", ...
            'ProjectNormal', true, ...
            'NormalProvider', @(X, ~, ~) kp.geometry.normalizeRows(X));
        solver.init(X0, U, xi, dt);

        for step = 1:nsteps
            solver.step();
        end

        X = solver.nodes();
        Rex = sqrt(R0^2 - 4 * T);
        r = vecnorm(X, 2, 2);
        radiusError = abs(mean(r) - Rex) / Rex;
        shapeError = norm(r - Rex) / max(norm(Rex * ones(size(r))), 1.0e-14);

        results(ix).N(level) = N;
        results(ix).h(level) = h;
        results(ix).dt(level) = dt;
        results(ix).nsteps(level) = nsteps;
        results(ix).radiusError(level) = radiusError;
        results(ix).shapeError(level) = shapeError;
        fprintf('  relative radius error: %.6e\n', radiusError);
        fprintf('  relative shape error: %.6e\n', shapeError);
    end

    for k = 2:numLevels
        results(ix).rate(k) = log(results(ix).radiusError(k - 1) / ...
            results(ix).radiusError(k)) / log(results(ix).h(k - 1) / results(ix).h(k));
    end
end

fig = figure('Color', 'w', 'Position', [100, 100, 760, 560]);
hold on;
colors = lines(numXi);
for ix = 1:numXi
    plot(results(ix).h, results(ix).radiusError, '-o', ...
        'LineWidth', 1.6, 'MarkerSize', 7, 'Color', colors(ix, :), ...
        'DisplayName', sprintf('\\xi = %d', results(ix).xi));
end
hold off;
set(gca, 'XScale', 'log', 'YScale', 'log', 'XDir', 'reverse');
grid on;
xlabel('h');
ylabel('Relative radius error');
title('Mean-curvature-flow sphere evolution');
legend('Location', 'best');

exportgraphics(fig, fullfile(pwd, 'docs', 'figures', ...
    'mean_curvature_flow_sphere_convergence.png'), 'Resolution', 180);
save(fullfile(pwd, 'mean_curvature_flow_sphere_results.mat'), 'results');
end
