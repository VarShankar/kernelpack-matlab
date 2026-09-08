function results = diffusion_convergence_3d(backend)
%DIFFUSION_CONVERGENCE_3D Run a 3D smooth-domain diffusion convergence study.

if nargin < 1 || isempty(backend)
    backend = "rbf";
end
backend = lower(string(backend));

orders = [2, 4, 6];
hValues = [0.28, 0.22, 0.18];
nu = 0.2;
tFinal = 0.5; % Temporal peak of 1 + 0.3 sin(pi t)

results = struct();
results.orders = orders;
results.h = hValues;
results.nu = nu;
results.tFinal = tFinal;
results.backend = backend;
results.assembler = "fd";
results.rows = repmat(struct( ...
    'order', 0, ...
    'h', 0, ...
    'dt', 0, ...
    'nPhysical', 0, ...
    'linf', 0, ...
    'l2', 0), numel(orders), numel(hValues));

surfaceSites = makeSmoothSurfaceSites(220);

phi = @(X) exp(0.35 * X(:, 1) - 0.2 * X(:, 2) + 0.25 * X(:, 3)) + ...
    0.15 * sin(1.1 * X(:, 1) - 0.7 * X(:, 2) + 0.9 * X(:, 3));
lapPhi = @(X) ...
    0.2250 * exp(0.35 * X(:, 1) - 0.2 * X(:, 2) + 0.25 * X(:, 3)) - ...
    0.3765 * sin(1.1 * X(:, 1) - 0.7 * X(:, 2) + 0.9 * X(:, 3));
amp = @(t) 1 + 0.3 * sin(pi * t);
ampDt = @(t) 0.3 * pi * cos(pi * t);

uExact = @(t, X) 2 + amp(t) .* phi(X);
forcing = @(nuValue, t, X) ampDt(t) .* phi(X) - nuValue * amp(t) .* lapPhi(X);
neuCoeff = @(Xb) zeros(size(Xb, 1), 1);
dirCoeff = @(Xb) ones(size(Xb, 1), 1);
bc = @(NeuCoeffs, DirCoeffs, nr, t, Xb) uExact(t, Xb); %#ok<INUSD>

for io = 1:numel(orders)
    xi = orders(io);
    for ih = 1:numel(hValues)
        h = hValues(ih);
        dt = 0.04 * h^2;
        domain = buildDomain(surfaceSites, h);

        solver = kp.solvers.DiffusionSolver( ...
            'LapAssembler', 'fd', ...
            'BCAssembler', 'fd', ...
            'LapStencil', backend, ...
            'BCStencil', backend);
        solver.init(domain, xi, dt, nu);

        Xphys = domain.getIntBdryNodes();
        solver.setStateHistory( ...
            uExact(tFinal - 3 * dt, Xphys), ...
            uExact(tFinal - 2 * dt, Xphys), ...
            uExact(tFinal - dt, Xphys));
        uNum = solver.bdf3Step(tFinal, forcing, neuCoeff, dirCoeff, bc);

        uTrue = uExact(tFinal, Xphys);
        err = uNum - uTrue;

        row = struct();
        row.order = xi;
        row.h = h;
        row.dt = dt;
        row.nPhysical = size(Xphys, 1);
        row.linf = max(abs(err));
        row.l2 = norm(err) / sqrt(numel(err));
        results.rows(io, ih) = row;
    end
end

results.rates = estimateRates(results.rows);
printResults(results);
end

function pts = makeSmoothSurfaceSites(n)
X = kp.geometry.fibonacciSphere(n);
uv = kp.geometry.cart2sphRows(X);
r = 1 + 0.12 * cos(3 * uv(:, 1)) .* cos(2 * uv(:, 2));
pts = X .* r;
end

function domain = buildDomain(surfaceSites, h)
surface = kp.geometry.EmbeddedSurface();
surface.setDataSites(surfaceSites);
surface.buildClosedGeometricModelPS(3, h, size(surfaceSites, 1));
surface.buildLevelSetFromGeometricModel([]);

generator = kp.nodes.DomainNodeGenerator();
domain = generator.buildDomainDescriptorFromGeometry(surface, h, ...
    'Seed', 17, ...
    'StripCount', 6, ...
    'DoOuterRefinement', true, ...
    'OuterFractionOfh', 0.5, ...
    'OuterRefinementZoneSizeAsMultipleOfh', 2.0);
end

function rates = estimateRates(rows)
rates = struct('order', [], 'linf', [], 'l2', []);
for io = 1:size(rows, 1)
    h = [rows(io, :).h];
    linf = [rows(io, :).linf];
    l2 = [rows(io, :).l2];
    rateRow = struct();
    rateRow.order = rows(io, 1).order;
    rateRow.linf = NaN(1, numel(h) - 1);
    rateRow.l2 = NaN(1, numel(h) - 1);
    for k = 1:numel(h) - 1
        rateRow.linf(k) = log(linf(k) / linf(k + 1)) / log(h(k) / h(k + 1));
        rateRow.l2(k) = log(l2(k) / l2(k + 1)) / log(h(k) / h(k + 1));
    end
    rates(io) = rateRow; %#ok<AGROW>
end
end

function printResults(results)
fprintf('\n3D diffusion convergence study (%s, %s)\n', upper(results.backend), upper(results.assembler));
fprintf('u(t,x,y,z) = 2 + (1 + 0.3 sin(pi t)) * phi(x,y,z)\n');
fprintf('Final time: t = %.3f (temporal peak)\n\n', results.tFinal);
for io = 1:size(results.rows, 1)
    fprintf('Order %d\n', results.rows(io, 1).order);
    fprintf('  h        dt        Nphys      Linf error      L2 error        Linf rate   L2 rate\n');
    for ih = 1:size(results.rows, 2)
        row = results.rows(io, ih);
        if ih == 1
            fprintf('  %-7.3f  %-8.5f  %-9d  %-14.6e  %-14.6e  %-10s  %-10s\n', ...
                row.h, row.dt, row.nPhysical, row.linf, row.l2, '-', '-');
        else
            fprintf('  %-7.3f  %-8.5f  %-9d  %-14.6e  %-14.6e  %-10.4f  %-10.4f\n', ...
                row.h, row.dt, row.nPhysical, row.linf, row.l2, ...
                results.rates(io).linf(ih - 1), results.rates(io).l2(ih - 1));
        end
    end
    fprintf('\n');
end
end
