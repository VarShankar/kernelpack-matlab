function poisson_solver_example_3d()
%POISSON_SOLVER_EXAMPLE_3D Solve a smooth-domain 3D Poisson problem with RBF-FD.

pts = makeSmoothSurfaceSites(180);

surface = kp.geometry.EmbeddedSurface();
surface.setDataSites(pts);
surface.buildClosedGeometricModelPS(3, 0.22, size(pts, 1));
surface.buildLevelSetFromGeometricModel([]);

generator = kp.nodes.DomainNodeGenerator();
domain = generator.buildDomainDescriptorFromGeometry(surface, 0.18, ...
    'Seed', 17, ...
    'StripCount', 6, ...
    'DoOuterRefinement', true, ...
    'OuterFractionOfh', 0.5, ...
    'OuterRefinementZoneSizeAsMultipleOfh', 2.0);

solver = kp.solvers.PoissonSolver( ...
    'LapAssembler', 'fd', ...
    'BCAssembler', 'fd', ...
    'LapStencil', 'rbf', ...
    'BCStencil', 'rbf');
solver.init(domain, 4);

uExact = @(X) exp(0.35 * X(:, 1) - 0.2 * X(:, 2) + 0.25 * X(:, 3)) + ...
    0.15 * sin(1.1 * X(:, 1) - 0.7 * X(:, 2) + 0.9 * X(:, 3));
forcing = @(Xeq) -( ...
    (0.35^2 + (-0.2)^2 + 0.25^2) .* exp(0.35 * Xeq(:, 1) - 0.2 * Xeq(:, 2) + 0.25 * Xeq(:, 3)) + ...
    0.15 * (-(1.1^2 + (-0.7)^2 + 0.9^2)) .* sin(1.1 * Xeq(:, 1) - 0.7 * Xeq(:, 2) + 0.9 * Xeq(:, 3)));
neuCoeff = @(Xb) zeros(size(Xb, 1), 1);
dirCoeff = @(Xb) ones(size(Xb, 1), 1);
bc = @(NeuCoeffs, DirCoeffs, nr, Xb) uExact(Xb); %#ok<INUSD>

result = solver.solve(forcing, neuCoeff, dirCoeff, bc);
Xphys = domain.getIntBdryNodes();
uTrue = uExact(Xphys);
absErr = abs(result.u - uTrue);

fprintf('3D smooth-domain Poisson example\n');
fprintf('  physical nodes: %d\n', size(Xphys, 1));
fprintf('  total nodes:    %d\n', domain.getNumTotalNodes());
fprintf('  max error:      %.6e\n', max(absErr));
fprintf('  rms error:      %.6e\n', norm(absErr) / sqrt(numel(absErr)));

plotResults(domain, result.u, absErr);
end

function pts = makeSmoothSurfaceSites(n)
X = kp.geometry.fibonacciSphere(n);
uv = kp.geometry.cart2sphRows(X);
r = 1 + 0.12 * cos(3 * uv(:, 1)) .* cos(2 * uv(:, 2));
pts = X .* r;
end

function plotResults(domain, u, absErr)
Xphys = domain.getIntBdryNodes();
Xb = domain.getBdryNodes();

figure('Name', 'Poisson Solver 3D', 'Color', 'w');
tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
plot3(Xphys(:, 1), Xphys(:, 2), Xphys(:, 3), 'k.', 'MarkerSize', 8);
hold on;
plot3(Xb(:, 1), Xb(:, 2), Xb(:, 3), 'r.', 'MarkerSize', 10);
axis equal;
grid on;
view(3);
title('Interior + boundary nodes');
xlabel('x');
ylabel('y');
zlabel('z');

nexttile;
scatter3(Xphys(:, 1), Xphys(:, 2), Xphys(:, 3), 18, u, 'filled');
axis equal;
grid on;
view(3);
title('Numerical solution');
xlabel('x');
ylabel('y');
zlabel('z');
colorbar;

nexttile;
tri = kp.geometry.triangulateClosedSurface(Xb);
trisurf(tri, Xb(:, 1), Xb(:, 2), Xb(:, 3), ...
    'FaceColor', [0.2 0.6 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
hold on;
scatter3(Xphys(:, 1), Xphys(:, 2), Xphys(:, 3), 18, absErr, 'filled');
axis equal;
grid on;
view(3);
camlight headlight;
lighting gouraud;
title('Absolute error');
xlabel('x');
ylabel('y');
zlabel('z');
colorbar;
end
