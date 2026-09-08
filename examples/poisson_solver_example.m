function poisson_solver_example()
%POISSON_SOLVER_EXAMPLE Solve a fixed-domain Poisson problem.

t = linspace(0, 2*pi, 80).';
t(end) = [];
curve = [cos(t), 0.8 * sin(t)];

% Build a smooth closed boundary and its stacked implicit representation.
surface = kp.geometry.EmbeddedSurface();
surface.setDataSites(curve);
surface.buildClosedGeometricModelPS(2, 0.06, size(curve, 1));
surface.buildLevelSetFromGeometricModel([]);

% Generate interior, boundary, and ghost nodes from that geometry.
generator = kp.nodes.DomainNodeGenerator();
domain = generator.buildDomainDescriptorFromGeometry(surface, 0.1, ...
    'Seed', 17, ...
    'StripCount', 5, ...
    'DoOuterRefinement', true, ...
    'OuterFractionOfh', 0.5, ...
    'OuterRefinementZoneSizeAsMultipleOfh', 2.0);

% Assemble and solve the Poisson problem with the requested stencil
% backend.
solver = kp.solvers.PoissonSolver( ...
    'LapAssembler', 'fd', ...
    'BCAssembler', 'fd', ...
    'LapStencil', 'wls', ...
    'BCStencil', 'wls');
solver.init(domain, 3);

uExact = @(X) X(:, 1).^2 + X(:, 2).^2;
forcing = @(Xeq) -4 * ones(size(Xeq, 1), 1);
neuCoeff = @(Xb) zeros(size(Xb, 1), 1);
dirCoeff = @(Xb) ones(size(Xb, 1), 1);
bc = @(NeuCoeffs, DirCoeffs, nr, Xb) uExact(Xb); %#ok<INUSD>

result = solver.solve(forcing, neuCoeff, dirCoeff, bc);
u = result.u;
Xphys = domain.getIntBdryNodes();
uTrue = uExact(Xphys);

figure('Color', 'w');
tiledlayout(1, 2);

nexttile;
scatter(Xphys(:, 1), Xphys(:, 2), 22, u, 'filled');
axis equal;
grid on;
title('Numerical solution');
colorbar;

nexttile;
scatter(Xphys(:, 1), Xphys(:, 2), 22, abs(u - uTrue), 'filled');
axis equal;
grid on;
title('Absolute error');
colorbar;
end
