function variable_poisson_solver_example()
%VARIABLE_POISSON_SOLVER_EXAMPLE End-to-end variable-coefficient Poisson solve.

t = linspace(0, 2*pi, 120).';
t(end) = [];
curve = [cos(t), 0.85 * sin(t)];

surface = kp.geometry.EmbeddedSurface();
surface.setDataSites(curve);
surface.buildClosedGeometricModelPS(2, 0.05, size(curve, 1));
surface.buildLevelSetFromGeometricModel([]);

generator = kp.nodes.DomainNodeGenerator();
domain = generator.buildDomainDescriptorFromGeometry(surface, 0.085, ...
    'Seed', 21, ...
    'StripCount', 5, ...
    'DoOuterRefinement', true, ...
    'OuterFractionOfh', 0.5, ...
    'OuterRefinementZoneSizeAsMultipleOfh', 2.0);

uExact = @(X) sin(X(:, 1)) + cos(0.75 * X(:, 2)) + 0.2 * X(:, 1) .* X(:, 2);
aCoeff = @(X) 1.5 + 0.3 * X(:, 1) - 0.2 * X(:, 2);
forcing = @(X) manufacturedForcing(X, aCoeff);

solver = kp.solvers.VariablePoissonSolver( ...
    'LapAssembler', 'fd', ...
    'BCAssembler', 'fd', ...
    'LapStencil', 'rbf', ...
    'BCStencil', 'rbf');
targetOrder = 4;
solver.init(domain, targetOrder);

result = solver.solve( ...
    forcing, ...
    aCoeff, ...
    @(Xb) zeros(size(Xb, 1), 1), ...
    @(Xb) ones(size(Xb, 1), 1), ...
    @(NeuCoeffs, DirCoeffs, nr, Xb) uExact(Xb)); %#ok<INUSD>

uTrue = uExact(domain.getIntBdryNodes());
err = result.u - uTrue;
fprintf('Variable-coefficient Poisson example\\n');
fprintf('  physical nodes: %d\\n', numel(result.u));
fprintf('  max error:      %.6e\\n', max(abs(err)));
fprintf('  rms error:      %.6e\\n', norm(err) / sqrt(numel(err)));
end

function f = manufacturedForcing(X, aCoeff)
x = X(:, 1);
y = X(:, 2);
a = aCoeff(X);
ux = cos(x) + 0.2 * y;
uy = -0.75 * sin(0.75 * y) + 0.2 * x;
lapU = -sin(x) - 0.5625 * cos(0.75 * y);
gradA_dot_gradU = 0.3 * ux - 0.2 * uy;
f = -(a .* lapU + gradA_dot_gradU);
end
