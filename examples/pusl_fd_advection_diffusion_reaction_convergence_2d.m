function results = pusl_fd_advection_diffusion_reaction_convergence_2d()
%PUSL_FD_ADVECTION_DIFFUSION_REACTION_CONVERGENCE_2D Convergence study for PU-SL + FD diffusion + reaction.

common = pusl_advection_diffusion_reaction_common();
orders = [2, 4, 6];
hValues = [0.14, 0.10, 0.08];
T = 0.25;
nu = 0.05;

results = struct();
results.orders = orders;
results.h = hValues;
results.rows = repmat(struct( ...
    'order', 0, ...
    'h', 0, ...
    'dt', 0, ...
    'nPhysical', 0, ...
    'rel', 0, ...
    'inf', 0, ...
    'mass', 0), numel(orders), numel(hValues));

for io = 1:numel(orders)
    xi = orders(io);
    for ih = 1:numel(hValues)
        h = hValues(ih);
        dt = 0.2 * h^2;
        domain = common.create_disk_domain(h);

        solver = kp.solvers.PUSLFDAdvectionDiffusionReactionSolver();
        solver.init(domain, xi, xi, dt, nu, "backward");
        solver.enableHomogeneousNeumannMassConservation(true);

        Xout = solver.getOutputNodes();
        u0 = common.exact_solution(0.0, Xout);
        m0 = solver.totalMass(u0);

        nsteps = max(ceil(T / max(dt, 1.0e-14)), 1);
        local_dt = T / nsteps;
        solver.setStepSize(local_dt);
        solver.setInitialState(u0);

        u = solver.bdf1Step(local_dt, ...
            @(t, X) common.boundary_zero_swirl(domain, t, X), [], ...
            @(t, state, X) common.reaction(t, state, X), ...
            @(nuValue, t, X) common.forcing(domain, nuValue, t, X), ...
            @(t, X) common.unit_coeff(t, X), ...
            @(t, X) common.zero_coeff(t, X), ...
            @(NeuCoeffs, DirCoeffs, nr, t, X) common.homogeneous_neumann(NeuCoeffs, DirCoeffs, nr, t, X));

        if nsteps >= 2
            u = solver.bdf2Step(2 * local_dt, ...
                @(t, X) common.boundary_zero_swirl(domain, t, X), [], ...
                @(t, state, X) common.reaction(t, state, X), ...
                @(nuValue, t, X) common.forcing(domain, nuValue, t, X), ...
                @(t, X) common.unit_coeff(t, X), ...
                @(t, X) common.zero_coeff(t, X), ...
                @(NeuCoeffs, DirCoeffs, nr, t, X) common.homogeneous_neumann(NeuCoeffs, DirCoeffs, nr, t, X));
        end

        if nsteps >= 3
            for step = 3:nsteps
                u = solver.bdf3Step(step * local_dt, ...
                    @(t, X) common.boundary_zero_swirl(domain, t, X), [], ...
                    @(t, state, X) common.reaction(t, state, X), ...
                    @(nuValue, t, X) common.forcing(domain, nuValue, t, X), ...
                    @(t, X) common.unit_coeff(t, X), ...
                    @(t, X) common.zero_coeff(t, X), ...
                    @(NeuCoeffs, DirCoeffs, nr, t, X) common.homogeneous_neumann(NeuCoeffs, DirCoeffs, nr, t, X));
            end
        end

        qex = common.exact_solution(T, Xout);
        row = struct();
        row.order = xi;
        row.h = h;
        row.dt = local_dt;
        row.nPhysical = size(Xout, 1);
        row.rel = norm(u - qex) / norm(qex);
        row.inf = norm(u - qex, inf);
        row.mass = abs(solver.totalMass(u) - m0);
        results.rows(io, ih) = row;
    end
end

results.rates = estimateRates(results.rows);
printResults(results);
save('pusl_fd_advection_diffusion_reaction_convergence_2d_results.mat', 'results');
end

function rates = estimateRates(rows)
rates = struct('order', [], 'rel', [], 'inf', []);
for io = 1:size(rows, 1)
    h = [rows(io, :).h];
    rel = [rows(io, :).rel];
    infv = [rows(io, :).inf];
    rateRow = struct();
    rateRow.order = rows(io, 1).order;
    rateRow.rel = NaN(1, numel(h) - 1);
    rateRow.inf = NaN(1, numel(h) - 1);
    for k = 1:numel(h) - 1
        rateRow.rel(k) = log(rel(k) / rel(k + 1)) / log(h(k) / h(k + 1));
        rateRow.inf(k) = log(infv(k) / infv(k + 1)) / log(h(k) / h(k + 1));
    end
    rates(io) = rateRow; %#ok<AGROW>
end
end

function printResults(results)
fprintf('\n2D PU-SL + FD diffusion + reaction convergence study\n');
for io = 1:size(results.rows, 1)
    fprintf('Order %d\n', results.rows(io, 1).order);
    fprintf('  h        dt        Nphys      rel error       inf error       rel rate    inf rate    mass drift\n');
    for ih = 1:size(results.rows, 2)
        row = results.rows(io, ih);
        if ih == 1
            fprintf('  %-7.3f  %-8.5f  %-9d  %-14.6e  %-14.6e  %-10s  %-10s  %-10.3e\n', ...
                row.h, row.dt, row.nPhysical, row.rel, row.inf, '-', '-', row.mass);
        else
            fprintf('  %-7.3f  %-8.5f  %-9d  %-14.6e  %-14.6e  %-10.4f  %-10.4f  %-10.3e\n', ...
                row.h, row.dt, row.nPhysical, row.rel, row.inf, ...
                results.rates(io).rel(ih - 1), results.rates(io).inf(ih - 1), row.mass);
        end
    end
    fprintf('\n');
end
end
