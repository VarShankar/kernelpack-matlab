function metrics = mean_curvature_flow_ellipsoid_example()
%MEAN_CURVATURE_FLOW_ELLIPSOID_EXAMPLE Short-time ellipsoid smoothing demo.

N = 1024;
xi = 4;
dt = 2.0e-4;
nsteps = 30;
axesLengths = [1.35, 0.85, 0.65];

geom = kp.geometry.ellipsoidSurface(N, axesLengths);
solver = kp.manifold.MeanCurvatureFlowSolver( ...
    'Mode', "semiimplicit", ...
    'ProjectNormal', true, ...
    'NormalStencilSize', 28);
solver.init(geom.X, geom.normals, xi, dt);

initialR = vecnorm(geom.X - mean(geom.X, 1), 2, 2);
for step = 1:nsteps
    solver.step();
end

Xf = solver.nodes();
finalR = vecnorm(Xf - mean(Xf, 1), 2, 2);
metrics = struct( ...
    'initialRadiusStd', std(initialR), ...
    'finalRadiusStd', std(finalR), ...
    'initialRadiusRange', max(initialR) - min(initialR), ...
    'finalRadiusRange', max(finalR) - min(finalR));

fig = figure('Color', 'w', 'Position', [100, 100, 1100, 460]);
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
scatter3(geom.X(:, 1), geom.X(:, 2), geom.X(:, 3), 14, initialR, 'filled');
axis equal tight;
grid on;
title('Initial Ellipsoid');
colorbar;
view(35, 25);

nexttile;
scatter3(Xf(:, 1), Xf(:, 2), Xf(:, 3), 14, finalR, 'filled');
axis equal tight;
grid on;
title('After Mean-Curvature Flow');
colorbar;
view(35, 25);

arrayfun(@disableDefaultInteractivity, findall(fig, 'Type', 'axes'));
exportgraphics(fig, fullfile(pwd, 'docs', 'figures', ...
    'mean_curvature_flow_ellipsoid_smoothing.png'), 'Resolution', 180);
save(fullfile(pwd, 'mean_curvature_flow_ellipsoid_results.mat'), 'metrics');

fprintf('Mean-curvature-flow ellipsoid smoothing\n');
fprintf('  initial radius std/range: %.6e / %.6e\n', ...
    metrics.initialRadiusStd, metrics.initialRadiusRange);
fprintf('  final radius std/range: %.6e / %.6e\n', ...
    metrics.finalRadiusStd, metrics.finalRadiusRange);
end
