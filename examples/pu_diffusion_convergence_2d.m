function results = pu_diffusion_convergence_2d()
%PU_DIFFUSION_CONVERGENCE_2D Convergence study for the MATLAB PU diffusion solver.

orders = [2, 4, 6];
hvals = [0.22, 0.16, 0.12];
dt_factor = 0.25;
nu = 0.05;
tfinal = 0.03;

fprintf('2D PU diffusion convergence study\n');
fprintf('u(t,x,y) = 1 + exp(-t) x (1-r^2)^2\n\n');

results = struct();
for oi = 1:numel(orders)
    xi = orders(oi);
    fprintf('Order %d\n', xi);
    fprintf('  h        dt         Nphys      Linf error      L2 error\n');
    errs_inf = zeros(numel(hvals), 1);
    errs_l2 = zeros(numel(hvals), 1);
    counts = zeros(numel(hvals), 1);
    for hi = 1:numel(hvals)
        h = hvals(hi);
        dt = dt_factor * h^2;
        steps = round(tfinal / dt);
        dt = tfinal / steps;

        domain = create_disk_domain(h);
        solver = kp.solvers.PUDiffusionSolver();
        solver.init(domain, xi, dt, nu);
        X = solver.getOutputNodes();
        counts(hi) = size(X, 1);

        t_hist = [0, dt, 2 * dt];
        U_hist = [exact_solution(t_hist(1), X), exact_solution(t_hist(2), X), exact_solution(t_hist(3), X)];
        solver.setStateHistory(U_hist(:, 1), U_hist(:, 2), U_hist(:, 3));
        u = U_hist(:, 3);
        t = 2 * dt;
        while t < tfinal - 1.0e-12
            t = t + dt;
            u = solver.bdf3Step(t, ...
                @(nuValue, time, Xq) forcing(nuValue, time, Xq), ...
                @(time, Xb) ones(size(Xb, 1), 1), ...
                @(time, Xb) zeros(size(Xb, 1), 1), ...
                @(NeuCoeffs, DirCoeffs, nr, time, Xb) zeros(size(Xb, 1), 1));
        end

        u_true = exact_solution(tfinal, X);
        diffv = u - u_true;
        errs_inf(hi) = norm(diffv, inf);
        errs_l2(hi) = norm(diffv) / sqrt(numel(diffv));
        fprintf('  %.3f    %.4f    %-9d%.6e    %.6e\n', h, dt, counts(hi), errs_inf(hi), errs_l2(hi));
    end
    fprintf('\n');
    results(oi).order = xi;
    results(oi).h = hvals;
    results(oi).count = counts;
    results(oi).linf = errs_inf;
    results(oi).l2 = errs_l2;
end

save('pu_diffusion_convergence_2d_results.mat', 'results');
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
