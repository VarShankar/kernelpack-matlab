function figures = moving_surface_adr_tp_spheroid_rearrangement_figures(varargin)
%MOVING_SURFACE_ADR_TP_SPHEROID_REARRANGEMENT_FIGURES Plot saved study tables.
%   Regenerates convergence, timing, quality, anisotropy, and comparison
%   figures from a saved spheroid rearrangement CSV without rerunning solves.

parser = inputParser();
parser.addParameter('InputFile', fullfile(pwd, ...
    'moving_surface_adr_tp_spheroid_rearrangement_transfer_all.csv'));
parser.addParameter('ImageDir', fullfile(pwd, 'docs', 'figures'));
parser.parse(varargin{:});
opts = parser.Results;

T = readtable(opts.InputFile, 'TextType', 'string');
T.elapsedPerStep = T.elapsedSeconds ./ max(T.nsteps, 1);
figures = writeStudyFigures(T, opts.ImageDir);
end

function figures = writeStudyFigures(T, imageDir)
if ~exist(imageDir, 'dir')
    mkdir(imageDir);
end
figures = struct();
figures.convergence = writeMetricFigures(T, imageDir, "relerr", ...
    '\textbf{relative $\ell_2$ error}', 'convergence');
figures.timing = writeMetricFigures(T, imageDir, "elapsedPerStep", ...
    'average runtime per timestep (s)', 'timing');
figures.quality = writeMetricFigures(T, imageDir, "maxQuality", ...
    'maximum nearest-neighbor quality', 'quality');
figures.anisotropy = writeMetricFigures(T, imageDir, "maxAnisotropy", ...
    'maximum tangent-stencil anisotropy', 'anisotropy');
figures.comparison = writeComparisonFigures(T, imageDir);
end

function paths = writeMetricFigures(T, imageDir, metricName, yLabelText, tag)
xiVals = unique(T.xi(T.status == "ok" & T.mode == "rearranged")).';
paths = strings(numel(xiVals), 1);
for ixi = 1:numel(xiVals)
    xi = xiVals(ixi);
    fig = figure('Color', 'w', 'Position', [100, 100, 760, 560]);
    plotMetricByM(T, xi, metricName, yLabelText);
    if metricName == "relerr"
        addReferenceSlope(T, xi);
    end
    title(sprintf('Spheroid rearrangement %s, \\xi = %d', ...
        readableMetricName(metricName), xi), 'Interpreter', 'tex');
    paths(ixi) = fullfile(imageDir, ...
        sprintf('moving_surface_adr_tp_spheroid_rearrangement_%s_xi%d.png', tag, xi));
    removeAxesToolbars(fig);
    kp.plot.exportPaperFigure(fig, paths(ixi));
    close(fig);
end
end

function paths = writeComparisonFigures(T, imageDir)
xiVals = unique(T.xi(T.status == "ok")).';
paths = strings(numel(xiVals), 1);
for ixi = 1:numel(xiVals)
    xi = xiVals(ixi);
    fig = figure('Color', 'w', 'Position', [100, 100, 760, 560]);
    hold on;
    plotRows(T, T.xi == xi & T.status == "ok" & T.mode == "lagrangian", ...
        "relerr", '-k', 'pure Lagrangian');
    plotRows(T, T.xi == xi & T.status == "ok" & T.mode == "rearranged" & ...
        T.transferMode == "sbf" & abs(T.MFactor - 0.5) < 10 * eps, ...
        "relerr", '-o', 'reduced SBF');
    plotRows(T, T.xi == xi & T.status == "ok" & T.mode == "rearranged" & ...
        T.transferMode == "localTp", "relerr", '-s', 'local tangent-plane');
    hold off;
    set(gca, 'XScale', 'log', 'YScale', 'log');
    grid on;
    xlabel('\textbf{$\sqrt{N}$}', 'Interpreter', 'latex');
    ylabel('\textbf{relative $\ell_2$ error}', 'Interpreter', 'latex');
    title(sprintf('Spheroid transfer comparison, \\xi = %d', xi), ...
        'Interpreter', 'tex');
    legend('Location', 'best');
    paths(ixi) = fullfile(imageDir, ...
        sprintf('moving_surface_adr_tp_spheroid_rearrangement_comparison_xi%d.png', xi));
    removeAxesToolbars(fig);
    kp.plot.exportPaperFigure(fig, paths(ixi));
    close(fig);
end
end

function plotMetricByM(T, xi, metricName, yLabelText)
hold on;
base = T.xi == xi & T.status == "ok" & T.mode == "rearranged";
mFactors = unique(T.MFactor(base & T.transferMode == "sbf")).';
numSeries = numel(mFactors) + any(base & T.transferMode == "localTp");
colors = lines(max(numSeries, 1));
colorIndex = 0;
for k = 1:numel(mFactors)
    colorIndex = colorIndex + 1;
    mFactor = mFactors(k);
    mask = base & T.transferMode == "sbf" & ...
        abs(T.MFactor - mFactor) < 10 * eps;
    plotRows(T, mask, metricName, '-o', sbfTransferLabel(mFactor), colors(colorIndex, :));
end
if any(base & T.transferMode == "localTp")
    colorIndex = colorIndex + 1;
    mask = base & T.transferMode == "localTp";
    plotRows(T, mask, metricName, '-s', 'local tangent-plane', colors(colorIndex, :));
end
hold off;
set(gca, 'XScale', 'log', 'YScale', 'log');
grid on;
xlabel('\textbf{$\sqrt{N}$}', 'Interpreter', 'latex');
if contains(string(yLabelText), "\ell")
    ylabel(yLabelText, 'Interpreter', 'latex');
else
    ylabel(yLabelText, 'Interpreter', 'none');
end
legend('Location', 'best');
end

function plotRows(T, mask, metricName, lineSpec, displayName, color)
if nargin < 6
    color = [];
end
S = T(mask, :);
if isempty(S)
    return;
end
[~, order] = sort(S.N);
S = S(order, :);
x = sqrt(S.N);
y = S.(metricName);
if isempty(color)
    loglog(x, y, lineSpec, 'LineWidth', 1.5, 'MarkerSize', 6, ...
        'DisplayName', displayName);
else
    loglog(x, y, lineSpec, 'LineWidth', 1.5, 'MarkerSize', 6, ...
        'Color', color, 'DisplayName', displayName);
end
end

function addReferenceSlope(T, xi)
hold on;
mFactors = T.MFactor(T.xi == xi & T.status == "ok" & T.mode == "rearranged" & ...
    T.transferMode == "sbf");
if isempty(mFactors)
    mask = T.xi == xi & T.status == "ok" & T.mode == "rearranged" & ...
        T.transferMode == "localTp";
else
    bestM = max(mFactors);
    mask = T.xi == xi & T.status == "ok" & T.mode == "rearranged" & ...
        T.transferMode == "sbf" & abs(T.MFactor - bestM) < 10 * eps;
end
S = T(mask, :);
if height(S) < 2
    return;
end
[~, order] = sort(S.N);
S = S(order, :);
x = sqrt(S.N);
y = S.relerr;
valid = isfinite(y) & y > 0;
if nnz(valid) < 2
    return;
end
x = x(valid);
y = y(valid);
x0 = x(end);
y0 = y(end);
xref = [x(1), x(end)];
yref = y0 * (xref / x0) .^ (-xi);
loglog(xref, yref, '--', 'Color', [0.25, 0.25, 0.25], ...
    'LineWidth', 1.0, 'DisplayName', sprintf('reference slope -%d', xi));
end

function label = sbfTransferLabel(mFactor)
if abs(mFactor - 1) < 10 * eps
    label = 'full SBF';
elseif mFactor < 1
    label = 'reduced SBF';
else
    label = 'SBF';
end
end

function removeAxesToolbars(fig)
axesList = findall(fig, 'Type', 'axes');
for iax = 1:numel(axesList)
    try
        axtoolbar(axesList(iax), {});
    catch
    end
end
end

function label = readableMetricName(metricName)
switch string(metricName)
    case "relerr"
        label = "convergence";
    case "elapsedPerStep"
        label = "timing";
    case "maxQuality"
        label = "marker quality";
    case "maxAnisotropy"
        label = "stencil anisotropy";
    otherwise
        label = char(metricName);
end
end
