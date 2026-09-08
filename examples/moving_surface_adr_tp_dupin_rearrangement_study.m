function study = moving_surface_adr_tp_dupin_rearrangement_study(varargin)
%MOVING_SURFACE_ADR_TP_DUPIN_REARRANGEMENT_STUDY Dupin rearrangement sweep.
%   Sweeps RBF-FD order xi, marker count N, and SBF control count M for the
%   fixed-shape Dupin cyclide with sliding labels.  The default PDE is forced
%   conservative advection with no diffusion or reaction; the numerical solve
%   uses RBF-FD surface divergence while the forcing is evaluated exactly.
%   Each run is sequential; optional parallelism is only used inside a run.

parser = inputParser();
parser.addParameter('XiVals', [2, 4, 6], @(x) isnumeric(x) && isvector(x));
parser.addParameter('NVals', [256, 512, 1024, 1600], @(x) isnumeric(x) && isvector(x));
parser.addParameter('MFactors', [0.25, 0.5, 1.0], @(x) isnumeric(x) && isvector(x));
parser.addParameter('FinalTime', 0.12, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('DtScale', 0.004, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('Nu', 0.0, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('FlowScale', 0.35, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('QualityThreshold', 1.75, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('MinStepsBetweenRearrangements', 20, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('DivergenceMode', "rbffd", @(x) isstring(x) || ischar(x));
parser.addParameter('UseExactStartup', true, @(x) islogical(x) && isscalar(x));
parser.addParameter('IncludeLagrangianComparison', true, @(x) islogical(x) && isscalar(x));
parser.addParameter('UseParallel', true, @(x) islogical(x) && isscalar(x));
parser.addParameter('ContinueOnError', true, @(x) islogical(x) && isscalar(x));
parser.addParameter('SaveAfterEach', true, @(x) islogical(x) && isscalar(x));
parser.addParameter('WriteFigures', true, @(x) islogical(x) && isscalar(x));
parser.addParameter('OutputPrefix', fullfile(pwd, ...
    'moving_surface_adr_tp_dupin_rearrangement_study'));
parser.addParameter('ImageDir', fullfile(pwd, 'docs', 'figures'));
parser.parse(varargin{:});
opts = parser.Results;

xiVals = unique(round(opts.XiVals(:).'), 'stable');
nVals = unique(round(opts.NVals(:).'), 'stable');
mFactors = unique(opts.MFactors(:).', 'stable');

rows = repmat(emptyRow(), 0, 1);
rowIndex = 0;

for xi = xiVals
    for N = nVals
        if opts.IncludeLagrangianComparison
            rowIndex = rowIndex + 1;
            rows(rowIndex, 1) = runStudyCase(opts, xi, N, NaN, false);
            saveProgress(rows, opts, false, false);
        end

        for mFactor = mFactors
            M = max(1, min(N, round(mFactor * N)));
            rowIndex = rowIndex + 1;
            rows(rowIndex, 1) = runStudyCase(opts, xi, N, M, true);
            saveProgress(rows, opts, false, false);
        end
    end
end

study = saveProgress(rows, opts, opts.WriteFigures, true);
end

function row = runStudyCase(opts, xi, N, M, useRearrangement)
row = emptyRow();
row.xi = xi;
row.N = N;
row.M = M;
row.useRearrangement = useRearrangement;
if useRearrangement
    row.mode = "rearranged";
    row.MFactor = M / N;
    mForRun = M;
else
    row.mode = "lagrangian";
    row.MFactor = NaN;
    mForRun = Inf;
end

fprintf('\nDupin rearrangement study xi=%d N=%d mode=%s M=%g\n', ...
    xi, N, row.mode, row.M);

try
    tic;
    single = moving_surface_adr_tp_dupin_rearrangement( ...
        'N', N, ...
        'Xi', xi, ...
        'SBFControlPointCount', mForRun, ...
        'UseRearrangement', useRearrangement, ...
        'FinalTime', opts.FinalTime, ...
        'DtScale', opts.DtScale, ...
        'Nu', opts.Nu, ...
        'FlowScale', opts.FlowScale, ...
        'QualityThreshold', opts.QualityThreshold, ...
        'MinStepsBetweenRearrangements', opts.MinStepsBetweenRearrangements, ...
        'DivergenceMode', opts.DivergenceMode, ...
        'UseExactStartup', opts.UseExactStartup, ...
        'UseParallel', opts.UseParallel, ...
        'ImagePath', '', ...
        'SaveOutputs', false);
    elapsed = toc;

    result = single.results(1);
    row.status = "ok";
    row.message = "";
    row.h = single.h;
    row.dt = single.dt;
    row.nsteps = single.nsteps;
    row.relerr = result.relerr;
    row.numRearrangements = result.numRearrangements;
    row.initialQuality = result.initialQuality;
    row.finalQuality = result.finalQuality;
    row.maxQuality = result.maxQuality;
    row.sbfControlPointCount = result.sbfControlPointCount;
    row.sbfControlPointFraction = result.sbfControlPointFraction;
    row.sbfFillDistance = result.sbfFillDistance;
    row.elapsedSeconds = elapsed;
    fprintf('  relerr=%.6e qmax=%.3f rearr=%d wall=%.2fs\n', ...
        row.relerr, row.maxQuality, row.numRearrangements, row.elapsedSeconds);
catch err
    row.status = "failed";
    row.message = string(err.message);
    fprintf('  FAILED: %s\n', err.message);
    if ~opts.ContinueOnError
        rethrow(err);
    end
end
end

function study = saveProgress(rows, opts, writeFigures, forceSave)
study = struct();
study.rows = rows;
study.table = struct2table(rows);
if opts.SaveAfterEach || writeFigures || forceSave
    save(string(opts.OutputPrefix) + ".mat", 'study');
    writetable(study.table, string(opts.OutputPrefix) + ".csv");
end
if writeFigures
    study.figures = writeStudyFigures(study.table, opts.ImageDir);
    save(string(opts.OutputPrefix) + ".mat", 'study');
end
end

function figures = writeStudyFigures(T, imageDir)
if ~exist(imageDir, 'dir')
    mkdir(imageDir);
end

figures = struct();
figures.convergence = writeConvergenceFigures(T, imageDir);
figures.timing = writeTimingFigures(T, imageDir);
figures.quality = writeQualityFigures(T, imageDir);
figures.lagrangianComparison = writeLagrangianComparisonFigures(T, imageDir);
end

function paths = writeConvergenceFigures(T, imageDir)
xiVals = unique(T.xi(T.status == "ok" & T.mode == "rearranged")).';
paths = strings(numel(xiVals), 1);
for ixi = 1:numel(xiVals)
    xi = xiVals(ixi);
    fig = figure('Color', 'w', 'Position', [100, 100, 760, 560]);
    plotMetricByM(T, xi, "relerr", '$\|e\|_2/\|c\|_2$');
    addReferenceSlope(T, xi);
    title(sprintf('Dupin forced-advection rearrangement, \\xi = %d', xi), 'Interpreter', 'tex');
    paths(ixi) = fullfile(imageDir, ...
        sprintf('moving_surface_adr_tp_dupin_rearrangement_convergence_xi%d.png', xi));
    exportgraphics(fig, paths(ixi), 'Resolution', 220);
    close(fig);
end
end

function paths = writeTimingFigures(T, imageDir)
xiVals = unique(T.xi(T.status == "ok" & T.mode == "rearranged")).';
paths = strings(numel(xiVals), 1);
for ixi = 1:numel(xiVals)
    xi = xiVals(ixi);
    fig = figure('Color', 'w', 'Position', [100, 100, 760, 560]);
    plotMetricByM(T, xi, "elapsedSeconds", 'wall time (s)');
    title(sprintf('Dupin forced-advection timings, \\xi = %d', xi), 'Interpreter', 'tex');
    paths(ixi) = fullfile(imageDir, ...
        sprintf('moving_surface_adr_tp_dupin_rearrangement_timing_xi%d.png', xi));
    exportgraphics(fig, paths(ixi), 'Resolution', 220);
    close(fig);
end
end

function paths = writeQualityFigures(T, imageDir)
xiVals = unique(T.xi(T.status == "ok" & T.mode == "rearranged")).';
paths = strings(numel(xiVals), 1);
for ixi = 1:numel(xiVals)
    xi = xiVals(ixi);
    fig = figure('Color', 'w', 'Position', [100, 100, 760, 560]);
    plotMetricByM(T, xi, "maxQuality", 'maximum nearest-neighbor quality');
    title(sprintf('Dupin forced-advection marker quality, \\xi = %d', xi), 'Interpreter', 'tex');
    paths(ixi) = fullfile(imageDir, ...
        sprintf('moving_surface_adr_tp_dupin_rearrangement_quality_xi%d.png', xi));
    exportgraphics(fig, paths(ixi), 'Resolution', 220);
    close(fig);
end
end

function paths = writeLagrangianComparisonFigures(T, imageDir)
xiVals = unique(T.xi(T.status == "ok")).';
paths = strings(numel(xiVals), 1);
for ixi = 1:numel(xiVals)
    xi = xiVals(ixi);
    fig = figure('Color', 'w', 'Position', [100, 100, 760, 560]);
    hold on;
    plotRows(T, T.xi == xi & T.status == "ok" & T.mode == "lagrangian", ...
        "relerr", '-k', 'pure Lagrangian');
    plotRows(T, T.xi == xi & T.status == "ok" & T.mode == "rearranged" & ...
        abs(T.MFactor - 1) < 10 * eps, "relerr", '-o', 'rearranged, M=N');
    hold off;
    set(gca, 'XScale', 'log', 'YScale', 'log');
    grid on;
    xlabel('$\sqrt{N}$', 'Interpreter', 'latex');
    ylabel('$\|e\|_2/\|c\|_2$', 'Interpreter', 'latex');
    legend('Location', 'best');
    title(sprintf('Pure Lagrangian versus rearranged, \\xi = %d', xi), 'Interpreter', 'tex');
    paths(ixi) = fullfile(imageDir, ...
        sprintf('moving_surface_adr_tp_dupin_rearrangement_lagrangian_comparison_xi%d.png', xi));
    exportgraphics(fig, paths(ixi), 'Resolution', 220);
    close(fig);
end
end

function plotMetricByM(T, xi, metricName, yLabelText)
hold on;
mFactors = unique(T.MFactor(T.xi == xi & T.status == "ok" & T.mode == "rearranged")).';
colors = lines(max(numel(mFactors), 1));
for k = 1:numel(mFactors)
    mFactor = mFactors(k);
    mask = T.xi == xi & T.status == "ok" & T.mode == "rearranged" & ...
        abs(T.MFactor - mFactor) < 10 * eps;
    plotRows(T, mask, metricName, '-o', sprintf('M/N = %.2g', mFactor), colors(k, :));
end
hold off;
set(gca, 'XScale', 'log', 'YScale', 'log');
grid on;
xlabel('$\sqrt{N}$', 'Interpreter', 'latex');
ylabel(yLabelText, 'Interpreter', 'latex');
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
mask = T.xi == xi & T.status == "ok" & T.mode == "rearranged" & ...
    abs(T.MFactor - max(T.MFactor(T.xi == xi & T.status == "ok" & ...
    T.mode == "rearranged"))) < 10 * eps;
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
    'LineWidth', 1.0, 'DisplayName', sprintf('slope -%d', xi));
end

function row = emptyRow()
row = struct( ...
    'status', "", ...
    'message', "", ...
    'xi', 0, ...
    'N', 0, ...
    'M', NaN, ...
    'MFactor', NaN, ...
    'mode', "", ...
    'useRearrangement', false, ...
    'h', NaN, ...
    'dt', NaN, ...
    'nsteps', 0, ...
    'relerr', NaN, ...
    'numRearrangements', 0, ...
    'initialQuality', NaN, ...
    'finalQuality', NaN, ...
    'maxQuality', NaN, ...
    'sbfControlPointCount', NaN, ...
    'sbfControlPointFraction', NaN, ...
    'sbfFillDistance', NaN, ...
    'elapsedSeconds', NaN);
end
