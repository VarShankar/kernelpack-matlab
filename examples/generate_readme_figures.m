function generate_readme_figures()
%GENERATE_README_FIGURES Render the figures shown in the top-level README.

rootDir = fileparts(fileparts(mfilename('fullpath')));
imageDir = fullfile(rootDir, 'docs', 'images');
if ~exist(imageDir, 'dir')
    mkdir(imageDir);
end

renderSmooth2DGeometry(fullfile(imageDir, 'readme_smooth_2d_geometry.png'));
renderGeometryClippedNodes(fullfile(imageDir, 'readme_geometry_clipped_nodes.png'));
renderRbfFdAssembly(fullfile(imageDir, 'readme_rbffd_operator.png'));
renderPoissonNeumann(fullfile(imageDir, 'readme_poisson_neumann.png'));
renderDiffusionStepping(fullfile(imageDir, 'readme_diffusion_stepping.png'));
end

function renderSmooth2DGeometry(outFile)
t = linspace(0, 2*pi, 50).';
t(end) = [];
curve = [cos(t), 0.7*sin(t)];

surface = kp.geometry.EmbeddedSurface();
surface.setDataSites(curve);
surface.buildClosedGeometricModelPS(2, 0.05, size(curve, 1));
surface.buildLevelSetFromGeometricModel([]);

xb = surface.getSampleSites();
nrmls = surface.getNrmls();

fig = figure('Color', 'w', 'Position', [100 100 1100 360]);
tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
plot(curve(:, 1), curve(:, 2), 'ko', 'MarkerFaceColor', [0.15 0.15 0.15], 'MarkerSize', 5);
axis equal;
grid on;
title('Data Sites');

nexttile;
plot(xb(:, 1), xb(:, 2), 'b.', 'MarkerSize', 12);
axis equal;
grid on;
title('Boundary Samples');

nexttile;
plot(xb(:, 1), xb(:, 2), 'b.', 'MarkerSize', 12);
hold on;
quiver(xb(:, 1), xb(:, 2), nrmls(:, 1), nrmls(:, 2), 0.05, 'Color', [0.82 0.18 0.18]);
axis equal;
grid on;
title('Boundary Normals');

finalizeReadmeFigure(fig);
exportgraphics(fig, outFile, 'Resolution', 180);
close(fig);
end

function renderGeometryClippedNodes(outFile)
t = linspace(0, 2*pi, 50).';
t(end) = [];
curve = [cos(t), 0.7*sin(t)];

surface = kp.geometry.EmbeddedSurface();
surface.setDataSites(curve);
surface.buildClosedGeometricModelPS(2, 0.05, size(curve, 1));
surface.buildLevelSetFromGeometricModel([]);

generator = kp.nodes.DomainNodeGenerator();
domain = generator.buildDomainDescriptorFromGeometry(surface, 0.08, ...
    'Seed', 17, ...
    'StripCount', 5, ...
    'DoOuterRefinement', true, ...
    'OuterFractionOfh', 0.5, ...
    'OuterRefinementZoneSizeAsMultipleOfh', 2.0);

Xi = domain.getInteriorNodes();
Xb = domain.getBdryNodes();
Xg = domain.getGhostNodes();

fig = figure('Color', 'w', 'Position', [100 100 720 560]);
plot(Xi(:, 1), Xi(:, 2), '.', 'Color', [0.15 0.15 0.15], 'MarkerSize', 10);
hold on;
plot(Xb(:, 1), Xb(:, 2), '.', 'Color', [0.85 0.15 0.15], 'MarkerSize', 12);
plot(Xg(:, 1), Xg(:, 2), '.', 'Color', [0.2 0.45 0.9], 'MarkerSize', 10);
axis equal;
grid on;
legend({'Interior', 'Boundary', 'Ghost'}, 'Location', 'best');
title('Geometry-Clipped Domain Nodes');
finalizeReadmeFigure(fig);
exportgraphics(fig, outFile, 'Resolution', 180);
close(fig);
end

function renderRbfFdAssembly(outFile)
t = linspace(0, 2*pi, 50).';
t(end) = [];
curve = [cos(t), 0.7*sin(t)];

surface = kp.geometry.EmbeddedSurface();
surface.setDataSites(curve);
surface.buildClosedGeometricModelPS(2, 0.05, size(curve, 1));
surface.buildLevelSetFromGeometricModel([]);

generator = kp.nodes.DomainNodeGenerator();
domain = generator.buildDomainDescriptorFromGeometry(surface, 0.08, ...
    'Seed', 17, ...
    'StripCount', 5, ...
    'DoOuterRefinement', true, ...
    'OuterFractionOfh', 0.5, ...
    'OuterRefinementZoneSizeAsMultipleOfh', 2.0);

sp = kp.rbffd.StencilProperties.fromAccuracy( ...
    'Operator', 'lap', ...
    'ConvergenceOrder', 4, ...
    'Dimension', 2, ...
    'Approximation', 'rbf', ...
    'treeMode', 'all', ...
    'pointSet', 'interior_boundary');
op = kp.rbffd.OpProperties('recordStencils', true);

assembler = kp.rbffd.FDDiffOp(@() kp.rbffd.RBFStencil());
assembler.AssembleOp(domain, 'lap', sp, op);
L = assembler.getOp();

fig = figure('Color', 'w', 'Position', [100 100 1000 420]);
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
Xi = domain.getInteriorNodes();
Xb = domain.getBdryNodes();
plot(Xi(:, 1), Xi(:, 2), '.', 'Color', [0.15 0.15 0.15], 'MarkerSize', 10);
hold on;
plot(Xb(:, 1), Xb(:, 2), '.', 'Color', [0.85 0.15 0.15], 'MarkerSize', 12);
axis equal;
grid on;
title('Domain Nodes');

nexttile;
spy(L);
title(sprintf('Laplacian Sparsity (%d x %d)', size(L, 1), size(L, 2)));

finalizeReadmeFigure(fig);
exportgraphics(fig, outFile, 'Resolution', 180);
close(fig);
end

function renderPoissonNeumann(outFile)
t = linspace(0, 2*pi, 120).';
t(end) = [];
curve = [cos(t), sin(t)];

surface = kp.geometry.EmbeddedSurface();
surface.setDataSites(curve);
surface.buildClosedGeometricModelPS(2, 0.06, size(curve, 1));
surface.buildLevelSetFromGeometricModel([]);

generator = kp.nodes.DomainNodeGenerator();
domain = generator.buildDomainDescriptorFromGeometry(surface, 0.08, ...
    'Seed', 17, ...
    'StripCount', 5);

solver = kp.solvers.PoissonSolver( ...
    'LapAssembler', 'fd', ...
    'BCAssembler', 'fd', ...
    'LapStencil', 'rbf', ...
    'BCStencil', 'rbf');
targetOrder = 4;
solver.init(domain, targetOrder);

uExact = @(X) (X(:,1).^2 + X(:,2).^2).^2 - (X(:,1).^2 + X(:,2).^2) + 1/6;
forcing = @(X) 4 - 16*(X(:,1).^2 + X(:,2).^2);
neuCoeff = @(Xb) ones(size(Xb,1), 1);
dirCoeff = @(Xb) zeros(size(Xb,1), 1);
bc = @(NeuCoeffs, DirCoeffs, nr, Xb) ...
    sum(([4*Xb(:,1).*(Xb(:,1).^2 + Xb(:,2).^2) - 2*Xb(:,1), ...
          4*Xb(:,2).*(Xb(:,1).^2 + Xb(:,2).^2) - 2*Xb(:,2)]).*nr, 2);

result = solver.solve(forcing, neuCoeff, dirCoeff, bc);
Xphys = domain.getIntBdryNodes();
tri = delaunay(Xphys(:, 1), Xphys(:, 2));
u = result.u;
uTrue = uExact(Xphys);
u = u - mean(u - uTrue);
absErr = abs(u - uTrue);

fig = figure('Color', 'w', 'Position', [100 100 1000 420]);
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
trisurf(tri, Xphys(:, 1), Xphys(:, 2), u, u, 'EdgeColor', 'none');
shading interp;
view(2);
axis equal tight;
grid on;
title('Pure-Neumann Poisson Solution');
colorbar;

nexttile;
trisurf(tri, Xphys(:, 1), Xphys(:, 2), absErr, absErr, 'EdgeColor', 'none');
shading interp;
view(2);
axis equal tight;
grid on;
title('Absolute Error');
colorbar;

finalizeReadmeFigure(fig);
exportgraphics(fig, outFile, 'Resolution', 180);
close(fig);
end

function renderDiffusionStepping(outFile)
t = linspace(0, 2*pi, 80).';
t(end) = [];
curve = [cos(t), 0.8*sin(t)];

surface = kp.geometry.EmbeddedSurface();
surface.setDataSites(curve);
surface.buildClosedGeometricModelPS(2, 0.06, size(curve, 1));
surface.buildLevelSetFromGeometricModel([]);

generator = kp.nodes.DomainNodeGenerator();
domain = generator.buildDomainDescriptorFromGeometry(surface, 0.1, ...
    'Seed', 17, ...
    'StripCount', 5, ...
    'DoOuterRefinement', true, ...
    'OuterFractionOfh', 0.5, ...
    'OuterRefinementZoneSizeAsMultipleOfh', 2.0);

solver = kp.solvers.DiffusionSolver( ...
    'LapAssembler', 'fd', ...
    'BCAssembler', 'fd', ...
    'LapStencil', 'rbf', ...
    'BCStencil', 'rbf');
nu = 0.25;
dt = 0.02;
targetOrder = 4;
solver.init(domain, targetOrder, dt, nu);

uExact = @(time, X) exp(-time) .* (X(:,1).^2 + X(:,2).^2);
forcing = @(nuValue, time, X) -exp(-time) .* (X(:,1).^2 + X(:,2).^2) ...
    - 4 * nuValue * exp(-time);
neuCoeff = @(Xb) zeros(size(Xb, 1), 1);
dirCoeff = @(Xb) ones(size(Xb, 1), 1);
bc = @(NeuCoeffs, DirCoeffs, nr, time, Xb) uExact(time, Xb);

tFinal = 0.5;
nSteps = round(tFinal / dt);
Xphys = domain.getIntBdryNodes();
tri = delaunay(Xphys(:, 1), Xphys(:, 2));

solver.setInitialState(uExact(0, Xphys));
states = {solver.currentPhysicalState()};
for step = 1:nSteps
    time = step * dt;
    if step == 1
        uNext = solver.bdf1Step(time, forcing, neuCoeff, dirCoeff, bc);
    elseif step == 2
        uNext = solver.bdf2Step(time, forcing, neuCoeff, dirCoeff, bc);
    else
        uNext = solver.bdf3Step(time, forcing, neuCoeff, dirCoeff, bc);
    end
    states{end + 1, 1} = uNext; %#ok<AGROW>
end

u1 = states{2};
u2 = states{3};
uFinal = states{end};
uTrueFinal = uExact(tFinal, Xphys);
absErr = abs(uFinal - uTrueFinal);

fig = figure('Color', 'w', 'Position', [100 100 1200 380]);
tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
trisurf(tri, Xphys(:, 1), Xphys(:, 2), u1, u1, 'EdgeColor', 'none');
shading interp;
view(2);
axis equal tight;
grid on;
title('BDF1 State');
colorbar;

nexttile;
trisurf(tri, Xphys(:, 1), Xphys(:, 2), u2, u2, 'EdgeColor', 'none');
shading interp;
view(2);
axis equal tight;
grid on;
title('BDF2 State');
colorbar;

nexttile;
trisurf(tri, Xphys(:, 1), Xphys(:, 2), absErr, absErr, 'EdgeColor', 'none');
shading interp;
view(2);
axis equal tight;
grid on;
title('Final Absolute Error');
colorbar;

finalizeReadmeFigure(fig);
exportgraphics(fig, outFile, 'Resolution', 180);
close(fig);
end

function finalizeReadmeFigure(fig)
axs = findall(fig, 'Type', 'axes');
for k = 1:numel(axs)
    disableDefaultInteractivity(axs(k));
end
end
