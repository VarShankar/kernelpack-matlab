function results = multispecies_pu_diffusion_convergence_2d()
%MULTISPECIES_PU_DIFFUSION_CONVERGENCE_2D Convergence study for the multispecies wrapper.

orders = [2, 4];
hvals = [0.22, 0.16, 0.12];
dt_factor = 0.25;
nu = 0.05;
tfinal = 0.03;

fprintf('2D multispecies PU diffusion convergence study\n\n');

results = struct();
for oi = 1:numel(orders)
    xi = orders(oi);
    fprintf('Order %d\n', xi);
    fprintf('  h        dt         Nphys      Fro error\n');
    err = zeros(numel(hvals), 1);
    counts = zeros(numel(hvals), 1);
    for hi = 1:numel(hvals)
        h = hvals(hi);
        dt = dt_factor * h^2;
        steps = round(tfinal / dt);
        dt = tfinal / steps;

        domain = create_disk_domain(h);
        solver = kp.solvers.MultiSpeciesPUDiffusionSolver();
        solver.init(domain, xi, dt, nu);
        X = solver.getOutputNodes();
        counts(hi) = size(X, 1);

        U2 = exact_solution(0, X);
        U1 = exact_solution(dt, X);
        U0 = exact_solution(2 * dt, X);
        solver.setStateHistory(U2, U1, U0);
        U = U0;
        t = 2 * dt;
        while t < tfinal - 1.0e-12
            t = t + dt;
            U = solver.bdf3Step(t, ...
                @(nuValue, time, Xq) forcing(nuValue, time, Xq), ...
                @(time, Xb) ones(size(Xb, 1), 1), ...
                @(time, Xb) zeros(size(Xb, 1), 1), ...
                @(NeuCoeffs, DirCoeffs, nr, time, Xb) zeros(size(Xb, 1), 2));
        end

        U_true = exact_solution(tfinal, X);
        err(hi) = norm(U - U_true, 'fro') / norm(U_true, 'fro');
        fprintf('  %.3f    %.4f    %-9d%.6e\n', h, dt, counts(hi), err(hi));
    end
    fprintf('\n');
    results(oi).order = xi;
    results(oi).h = hvals;
    results(oi).count = counts;
    results(oi).fro = err;
end

save('multispecies_pu_diffusion_convergence_2d_results.mat', 'results');
end

function U = exact_solution(t, X)
u = ones(size(X, 1), 1) + exp(-t) .* phi(X);
U = [u, 0.5 * u];
end

function F = forcing(nu, t, X)
f = exp(-t) .* (-phi(X) - nu .* lap_phi(X));
F = [f, 0.5 * f];
end

function values = phi(X)
r2 = sum(X.^2, 2);
values = X(:, 1) .* (1 - r2).^2;
end

function values = lap_phi(X)
r2 = sum(X.^2, 2);
values = 4 * X(:, 1) .* (6 * r2 - 4);
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
