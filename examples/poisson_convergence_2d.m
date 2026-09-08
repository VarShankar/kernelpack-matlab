function results = poisson_convergence_2d(backend)
%POISSON_CONVERGENCE_2D Run a 2D Poisson convergence study.

if nargin < 1 || isempty(backend)
    backend = "wls";
end
backend = lower(string(backend));

orders = [2, 4, 6];
hValues = [0.14, 0.10, 0.07, 0.05];

results = struct();
results.orders = orders;
results.h = hValues;
results.backend = backend;
results.assembler = "fd";
results.rows = repmat(struct( ...
    'order', 0, ...
    'h', 0, ...
    'nPhysical', 0, ...
    'linf', 0, ...
    'l2', 0), numel(orders), numel(hValues));

curveSites = makeCircleDataSites(240);
uExact = @(X) exp(X(:, 1) + 0.5 * X(:, 2));
forcing = @(Xeq) -1.25 * exp(Xeq(:, 1) + 0.5 * Xeq(:, 2));
neuCoeff = @(Xb) zeros(size(Xb, 1), 1);
dirCoeff = @(Xb) ones(size(Xb, 1), 1);
bc = @(NeuCoeffs, DirCoeffs, nr, Xb) uExact(Xb); %#ok<INUSD>

for io = 1:numel(orders)
    xi = orders(io);
    for ih = 1:numel(hValues)
        h = hValues(ih);
        domain = buildDomain(curveSites, h);

        solver = kp.solvers.PoissonSolver( ...
            'LapAssembler', 'fd', ...
            'BCAssembler', 'fd', ...
            'LapStencil', backend, ...
            'BCStencil', backend);
        solver.init(domain, xi);
        solveResult = solver.solve(forcing, neuCoeff, dirCoeff, bc);

        Xphys = domain.getIntBdryNodes();
        uTrue = uExact(Xphys);
        err = solveResult.u - uTrue;

        row = struct();
        row.order = xi;
        row.h = h;
        row.nPhysical = size(Xphys, 1);
        row.linf = max(abs(err));
        row.l2 = norm(err) / sqrt(numel(err));
        results.rows(io, ih) = row;
    end
end

results.rates = estimateRates(results.rows);
printResults(results);
end

function domain = buildDomain(curveSites, h)
surface = kp.geometry.EmbeddedSurface();
surface.setDataSites(curveSites);
surface.buildClosedGeometricModelPS(2, h, size(curveSites, 1));
surface.buildLevelSetFromGeometricModel([]);

generator = kp.nodes.DomainNodeGenerator();
domain = generator.buildDomainDescriptorFromGeometry(surface, h, ...
    'Seed', 17, ...
    'StripCount', 5);
end

function X = makeCircleDataSites(n)
t = linspace(0, 2 * pi, n + 1).';
t(end) = [];
X = [cos(t), sin(t)];
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
fprintf('\n2D Poisson convergence study (%s, %s)\n', upper(results.backend), upper(results.assembler));
fprintf('Exact solution: u(x,y) = exp(x + 0.5 y)\n\n');
for io = 1:size(results.rows, 1)
    fprintf('Order %d\n', results.rows(io, 1).order);
    fprintf('  h        Nphys      Linf error      L2 error        Linf rate   L2 rate\n');
    for ih = 1:size(results.rows, 2)
        row = results.rows(io, ih);
        if ih == 1
            fprintf('  %-7.3f  %-9d  %-14.6e  %-14.6e  %-10s  %-10s\n', ...
                row.h, row.nPhysical, row.linf, row.l2, '-', '-');
        else
            fprintf('  %-7.3f  %-9d  %-14.6e  %-14.6e  %-10.4f  %-10.4f\n', ...
                row.h, row.nPhysical, row.linf, row.l2, ...
                results.rates(io).linf(ih - 1), results.rates(io).l2(ih - 1));
        end
    end
    fprintf('\n');
end
end
