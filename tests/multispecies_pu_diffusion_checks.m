function multispecies_pu_diffusion_checks()
%MULTISPECIES_PU_DIFFUSION_CHECKS Basic checks for the multispecies PU diffusion wrapper.

domain = create_disk_domain(0.18);
solver = kp.solvers.MultiSpeciesPUDiffusionSolver();
solver.init(domain, 4, 0.01, 0.05);

X = solver.getOutputNodes();
u0 = [exact_solution(0.0, X), 0.5 * exact_solution(0.0, X)];
solver.setInitialState(u0);

u1 = solver.bdf1Step(0.01, ...
    @(nuValue, t, Xq) [forcing(nuValue, t, Xq), 0.5 * forcing(nuValue, t, Xq)], ...
    @(t, Xb) ones(size(Xb, 1), 1), ...
    @(t, Xb) zeros(size(Xb, 1), 1), ...
    @(NeuCoeffs, DirCoeffs, nr, t, Xb) zeros(size(Xb, 1), 2));

uTrue = [exact_solution(0.01, X), 0.5 * exact_solution(0.01, X)];
relerr = norm(u1 - uTrue, 'fro') / norm(uTrue, 'fro');
assert(relerr < 3e-1, 'MultiSpeciesPUDiffusionSolver should track the manufactured solution to modest accuracy.');

disp('multispecies pu diffusion checks passed');
end

function values = phi(X)
r2 = sum(X.^2, 2);
values = X(:, 1) .* (1 - r2).^2;
end

function values = lap_phi(X)
r2 = sum(X.^2, 2);
values = 4 * X(:, 1) .* (6 * r2 - 4);
end

function values = exact_solution(t, X)
values = ones(size(X, 1), 1) + exp(-t) .* phi(X);
end

function values = forcing(nu, t, X)
values = exp(-t) .* (-phi(X) - nu .* lap_phi(X));
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
