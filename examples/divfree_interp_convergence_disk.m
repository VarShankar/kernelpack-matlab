function results = divfree_interp_convergence_disk()
%DIVFREE_INTERP_CONVERGENCE_DISK Serial convergence study for local
% divergence-free PHS+poly interpolation on the unit disk.

orders = [2, 4, 6];
hVals = [0.22, 0.18, 0.14, 0.11];
curveSites = makeCircleDataSites(320);
Xeval = makeEvaluationCloud();
[Ueval, target] = targetField(Xeval);

results = struct();
results.orders = orders;
results.h = hVals;
results.rows = repmat(struct('order', 0, 'h', 0, 'n', 0, 'linf', 0, 'l2', 0), numel(orders), numel(hVals));

for io = 1:numel(orders)
    xi = orders(io);
    for ih = 1:numel(hVals)
        h = hVals(ih);
        domain = buildDomain(curveSites, h);
        X = domain.getIntBdryNodes();
        U = targetField(X);
        [linfErr, l2Err] = runLocalInterpolation(X, U, Xeval, Ueval, xi);
        results.rows(io, ih) = struct( ...
            'order', xi, ...
            'h', h, ...
            'n', size(X, 1), ...
            'linf', linfErr, ...
            'l2', l2Err);
        fprintf('xi=%d, h=%.3f, N=%d, Linf=%.6e, L2=%.6e\\n', xi, h, size(X, 1), linfErr, l2Err);
    end
end

results.rates = estimateRates(results.rows);
results.figurePath = fullfile(pwd, 'docs', 'figures', 'divfree_interp_convergence_disk.png');
plotResults(results, target, results.figurePath);
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

function X = makeEvaluationCloud()
grid = linspace(-0.92, 0.92, 81);
[xx, yy] = meshgrid(grid, grid);
mask = (xx.^2 + yy.^2) <= 0.92^2;
X = [xx(mask), yy(mask)];
end

function [U, target] = targetField(X)
x = X(:, 1);
y = X(:, 2);
target.psi = @(xv, yv) sin(pi * xv) .* sin(pi * yv) + 0.2 * xv.^3 .* yv;
target.u = @(xv, yv) pi * sin(pi * xv) .* cos(pi * yv) + 0.2 * xv.^3;
target.v = @(xv, yv) -pi * cos(pi * xv) .* sin(pi * yv) - 0.6 * xv.^2 .* yv;
U = [target.u(x, y), target.v(x, y)];
end

function [linfErr, l2Err] = runLocalInterpolation(X, U, Xeval, Ueval, xi)
ell = xi - 1;
dim = size(X, 2);
npoly = dim * nchoosek(ell + dim, dim) - nchoosek(ell + dim - 1, dim);
stencilSize = 2 * ceil(npoly / dim) + 1;

% The r^1 divergence-free kernel is unstable on this disk test. Keep the
% polynomial order tied to xi, but use at least the cubic PHS kernel. The
% lowest-order run also benefits from a slightly richer local stencil.
phsDegree = max(3, min(ell - mod(ell + 1, 2), 7));
if xi == 2
    phsDegree = max(phsDegree, 5);
    stencilSize = max(stencilSize, 13);
end

recurrence = @(N) kp.poly.jacobi_recurrence(N, 0, 0);
polyOpts = struct('use_pca_frame', false, 'rank_tol', 1e-12);
tree = KDTreeSearcher(X);
centerIdx = knnsearch(tree, Xeval, 'K', 1);
interp = kp.divfree.LocalDivFreeInterpolator.fit(X, U, ell, phsDegree, stencilSize, ...
    'Recurrence', recurrence, ...
    'PolyOptions', polyOpts, ...
    'ActiveCenters', unique(centerIdx, 'stable'));
Uapprox = interp.evaluate(Xeval, 'CenterIndices', centerIdx);

err = Uapprox - Ueval;
rowNorm = sqrt(sum(err.^2, 2));
trueNorm = sqrt(sum(Ueval.^2, 2));
linfErr = max(rowNorm);
l2Err = norm(rowNorm) / norm(trueNorm);
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

function plotResults(results, target, figurePath)
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 700 520]);
hold on;
markers = {'o-', 's-', 'd-'};
for io = 1:numel(results.orders)
    plot([results.rows(io, :).h], [results.rows(io, :).l2], markers{io}, ...
        'LineWidth', 1.6, 'MarkerSize', 8);
end
set(gca, 'XScale', 'log', 'YScale', 'log');
grid on;
xlabel('h');
ylabel('Relative L2 error');
legend(compose('\\xi = %d', results.orders), 'Location', 'northwest');
title('Div-free PHS+poly interpolation on the disk');
drawnow;
exportgraphics(fig, figurePath, 'Resolution', 180);
close(fig);
end

function printResults(results)
fprintf('\\nDiv-free disk interpolation convergence\\n');
for io = 1:numel(results.orders)
    fprintf('Order %d\\n', results.orders(io));
    fprintf('  h        N        Linf error      L2 error        Linf rate   L2 rate\\n');
    for ih = 1:numel(results.h)
        row = results.rows(io, ih);
        if ih == 1
            linfRate = '-';
            l2Rate = '-';
        else
            linfRate = sprintf('%8.4f', results.rates(io).linf(ih - 1));
            l2Rate = sprintf('%8.4f', results.rates(io).l2(ih - 1));
        end
        fprintf('  %.3f    %-7d  %.6e    %.6e    %s   %s\\n', ...
            row.h, row.n, row.linf, row.l2, linfRate, l2Rate);
    end
end
end
