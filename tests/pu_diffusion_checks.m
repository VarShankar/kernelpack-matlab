function pu_diffusion_checks()
%PU_DIFFUSION_CHECKS Basic checks for the MATLAB PU diffusion solver.

domain = create_disk_domain(0.18);
solver = kp.solvers.PUDiffusionSolver();
solver.init(domain, 4, 0.01, 0.05);

X = solver.getOutputNodes();
u0 = exact_solution(0.0, X);
solver.setInitialState(u0);

u1 = solver.bdf1Step(0.01, ...
    @(nuValue, t, Xq) forcing(nuValue, t, Xq), ...
    @(t, Xb) ones(size(Xb, 1), 1), ...
    @(t, Xb) zeros(size(Xb, 1), 1), ...
    @(NeuCoeffs, DirCoeffs, nr, t, Xb) zeros(size(Xb, 1), 1));

uTrue = exact_solution(0.01, X);
relerr = norm(u1 - uTrue) / norm(uTrue);
assert(relerr < 3e-1, 'PU diffusion BDF1 should track the manufactured solution to modest accuracy.');

disp('pu diffusion checks passed');
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
