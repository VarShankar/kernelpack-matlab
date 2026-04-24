function diffusion_solver_example()
%DIFFUSION_SOLVER_EXAMPLE Advance a fixed-domain diffusion problem with BDF steps.

t = linspace(0, 2*pi, 80).';
t(end) = [];
curve = [cos(t), 0.8 * sin(t)];

% Build a smooth closed boundary and its implicit model.
surface = kp.geometry.EmbeddedSurface();
surface.setDataSites(curve);
surface.buildClosedGeometricModelPS(2, 0.06, size(curve, 1));
surface.buildLevelSetFromGeometricModel([]);

% Generate the geometry-backed domain used by the diffusion stepper.
generator = kp.nodes.DomainNodeGenerator();
domain = generator.buildDomainDescriptorFromGeometry(surface, 0.1, ...
    'Seed', 23, ...
    'StripCount', 5, ...
    'DoOuterRefinement', true, ...
    'OuterFractionOfh', 0.5, ...
    'OuterRefinementZoneSizeAsMultipleOfh', 2.0);

% Initialize the diffusion solver once, then step it with a BDF update.
solver = kp.solvers.DiffusionSolver( ...
    'LapAssembler', 'fd', ...
    'BCAssembler', 'fd', ...
    'LapStencil', 'wls', ...
    'BCStencil', 'wls');
nu = 0.25;
dt = 0.02;
solver.init(domain, 3, dt, nu);

uExact = @(time, X) exp(-time) .* (X(:, 1).^2 + X(:, 2).^2);
forcing = @(nuValue, time, X) -exp(-time) .* (X(:, 1).^2 + X(:, 2).^2) - 4 * nuValue * exp(-time);
neuCoeff = @(Xb) zeros(size(Xb, 1), 1);
dirCoeff = @(Xb) ones(size(Xb, 1), 1);
bc = @(NeuCoeffs, DirCoeffs, nr, time, Xb) uExact(time, Xb); %#ok<INUSD>

solver.setInitialState(uExact(0, domain.getIntBdryNodes()));
u1 = solver.bdf1Step(dt, forcing, neuCoeff, dirCoeff, bc);
u2 = solver.bdf2Step(2 * dt, forcing, neuCoeff, dirCoeff, bc);
u3 = solver.bdf3Step(3 * dt, forcing, neuCoeff, dirCoeff, bc);

Xphys = domain.getIntBdryNodes();

figure('Color', 'w');
tiledlayout(1, 3);

nexttile;
scatter(Xphys(:, 1), Xphys(:, 2), 22, u1, 'filled');
axis equal;
grid on;
title('BDF1 state');
colorbar;

nexttile;
scatter(Xphys(:, 1), Xphys(:, 2), 22, u2, 'filled');
axis equal;
grid on;
title('BDF2 state');
colorbar;

nexttile;
scatter(Xphys(:, 1), Xphys(:, 2), 22, abs(u3 - uExact(3 * dt, Xphys)), 'filled');
axis equal;
grid on;
title('BDF3 absolute error');
colorbar;
end
