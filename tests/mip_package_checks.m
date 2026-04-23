function mip_package_checks()
%MIP_PACKAGE_CHECKS Fast package check for MIP installs.

domain = create_disk_domain(0.18);

sp = kp.rbffd.StencilProperties.fromAccuracy( ...
    'Operator', 'lap', ...
    'ConvergenceOrder', 4, ...
    'Dimension', 2, ...
    'Approximation', 'rbf', ...
    'treeMode', 'all', ...
    'pointSet', 'interior_boundary');

assembler = kp.rbffd.FDDiffOp(@() kp.rbffd.RBFStencil());
assembler.AssembleOp(domain, 'lap', sp, kp.rbffd.OpProperties());
L = assembler.getOp();
assert(~isempty(L) && size(L, 2) == domain.getNumTotalNodes(), ...
    'RBF-FD assembly failed in the MIP package check.');

solver = kp.solvers.PoissonSolver( ...
    'LapAssembler', 'fd', ...
    'BCAssembler', 'fd', ...
    'LapStencil', 'rbf', ...
    'BCStencil', 'rbf');
targetOrder = 4;
solver.init(domain, targetOrder);
X = domain.getIntBdryNodes();
uExact = @(Xq) Xq(:, 1).^2 - 0.5 * Xq(:, 2).^2;
forcing = @(Xq) -ones(size(Xq, 1), 1);

result = solver.solve( ...
    forcing, ...
    @(Xb) zeros(size(Xb, 1), 1), ...
    @(Xb) ones(size(Xb, 1), 1), ...
    @(NeuCoeffs, DirCoeffs, nr, Xb) uExact(Xb));

relerr = norm(result.u - uExact(X)) / norm(uExact(X));
assert(relerr < 2e-1, 'PoissonSolver package check solve failed.');

disp('mip package checks passed');
end

function domain = create_disk_domain(h)
t = linspace(0, 2 * pi, 101).';
t(end) = [];
data_sites = [cos(t), sin(t)];

surface = kp.geometry.EmbeddedSurface();
surface.setDataSites(data_sites);
surface.buildClosedGeometricModelPS(2, h, size(data_sites, 1));
surface.buildLevelSetFromGeometricModel([]);

generator = kp.nodes.DomainNodeGenerator();
domain = generator.buildDomainDescriptorFromGeometry(surface, h, ...
    'Seed', 17, ...
    'StripCount', 5);
end
