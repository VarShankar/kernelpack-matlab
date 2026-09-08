function study = moving_surface_adr_tp_spheroid_rearrangement_study(varargin)
%MOVING_SURFACE_ADR_TP_SPHEROID_REARRANGEMENT_STUDY Spheroid rearrangement sweep.
%   Runs a forced-advection convergence/timing study for Lagrangian marker
%   rearrangement with semi-Lagrangian BDF history reconstruction.  The
%   numerical update uses tangent-plane RBF-FD divergence.  Exact calculus is
%   used only for forcing and final reference evaluation.

parser = inputParser();
parser.addParameter('XiVals', [2, 4, 6], @(x) isnumeric(x) && isvector(x));
parser.addParameter('NVals', [512, 768, 1024, 1152], @(x) isnumeric(x) && isvector(x));
parser.addParameter('MFactors', 1.0, @(x) isnumeric(x) && isvector(x));
parser.addParameter('TransferModes', "localTp");
parser.addParameter('FinalTime', 0.35, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('DtScale', 0.08, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('FlowScale', 0.8, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('RearrangementTimes', 0.175, @(x) isnumeric(x));
parser.addParameter('ExactRearrangementBootstrap', false, @(x) islogical(x) && isscalar(x));
parser.addParameter('MassCorrectionMode', "balance", @(x) isstring(x) || ischar(x));
parser.addParameter('MassCorrectionMaxRelativeCorrection', Inf, ...
    @(x) isnumeric(x) && isscalar(x) && x > 0);
parser.addParameter('MaterialSampler', "fps", @(x) isstring(x) || ischar(x));
parser.addParameter('SamplerCandidateFactor', 8, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('QualityThreshold', 1.55, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('PredictiveLookaheadSteps', 3, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('MinStepsBetweenRearrangements', 8, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('IncludeLagrangianComparison', true, @(x) islogical(x) && isscalar(x));
parser.addParameter('UseParallel', false, @(x) islogical(x) && isscalar(x));
parser.addParameter('ContinueOnError', true, @(x) islogical(x) && isscalar(x));
parser.addParameter('SaveAfterEach', true, @(x) islogical(x) && isscalar(x));
parser.addParameter('WriteFigures', true, @(x) islogical(x) && isscalar(x));
parser.addParameter('OutputPrefix', fullfile(pwd, ...
    'moving_surface_adr_tp_spheroid_rearrangement_study'));
parser.addParameter('ImageDir', fullfile(pwd, 'docs', 'figures'));
parser.parse(varargin{:});
opts = parser.Results;

xiVals = unique(round(opts.XiVals(:).'), 'stable');
nVals = unique(round(opts.NVals(:).'), 'stable');
mFactors = unique(opts.MFactors(:).', 'stable');
transferModes = string(opts.TransferModes(:).');

rows = repmat(emptyRow(), 0, 1);
rowIndex = 0;
for xi = xiVals
    for N = nVals
        if opts.IncludeLagrangianComparison
            rowIndex = rowIndex + 1;
            rows(rowIndex, 1) = runStudyCase(opts, xi, N, NaN, "none", false);
            saveProgress(rows, opts, false, false);
        end
        for transferMode = transferModes
            switch lower(transferMode)
                case "sbf"
                    for mFactor = mFactors
                        M = max(1, min(N, round(mFactor * N)));
                        rowIndex = rowIndex + 1;
                        rows(rowIndex, 1) = runStudyCase(opts, xi, N, M, transferMode, true);
                        saveProgress(rows, opts, false, false);
                    end
                case {"localtp", "localtangentplane"}
                    rowIndex = rowIndex + 1;
                    rows(rowIndex, 1) = runStudyCase(opts, xi, N, NaN, "localTp", true);
                    saveProgress(rows, opts, false, false);
                otherwise
                    error('kp:examples:BadTransferMode', ...
                        'Unknown transfer mode "%s".', transferMode);
            end
        end
    end
end

study = saveProgress(rows, opts, opts.WriteFigures, true);
end

function row = runStudyCase(opts, xi, N, M, transferMode, useRearrangement)
row = emptyRow();
row.xi = xi;
row.N = N;
row.M = M;
row.transferMode = string(transferMode);
row.useRearrangement = useRearrangement;
if useRearrangement
    row.mode = "rearranged";
    if lower(string(transferMode)) == "sbf"
        row.MFactor = M / N;
        mForRun = M;
    else
        row.MFactor = NaN;
        mForRun = NaN;
    end
else
    row.mode = "lagrangian";
    row.transferMode = "none";
    row.MFactor = NaN;
    mForRun = Inf;
end

fprintf('\nSpheroid rearrangement study xi=%d N=%d mode=%s transfer=%s M=%g\n', ...
    xi, N, row.mode, row.transferMode, row.M);
try
    tic;
    single = moving_surface_adr_tp_spheroid_rearrangement( ...
        'N', N, ...
        'Xi', xi, ...
        'SBFControlPointCount', mForRun, ...
        'TransferMode', transferMode, ...
        'UseRearrangement', useRearrangement, ...
        'FinalTime', opts.FinalTime, ...
        'DtScale', opts.DtScale, ...
        'FlowScale', opts.FlowScale, ...
        'RearrangementTimes', opts.RearrangementTimes, ...
        'ExactRearrangementBootstrap', opts.ExactRearrangementBootstrap, ...
        'MassCorrectionMode', opts.MassCorrectionMode, ...
        'MassCorrectionMaxRelativeCorrection', opts.MassCorrectionMaxRelativeCorrection, ...
        'MaterialSampler', opts.MaterialSampler, ...
        'SamplerCandidateFactor', opts.SamplerCandidateFactor, ...
        'QualityThreshold', opts.QualityThreshold, ...
        'PredictiveLookaheadSteps', opts.PredictiveLookaheadSteps, ...
        'MinStepsBetweenRearrangements', opts.MinStepsBetweenRearrangements, ...
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
    row.massRelError = result.massRelError;
    row.balanceRelResidual = result.balanceRelResidual;
    row.massCorrectionSteps = result.massCorrectionSteps;
    row.massCorrectionMaxRelative = result.massCorrectionMaxRelative;
    row.massCorrectionMaxPointShift = result.massCorrectionMaxPointShift;
    row.numRearrangements = result.numRearrangements;
    row.numPredictiveRearrangements = result.numPredictiveRearrangements;
    row.numThresholdRearrangements = result.numThresholdRearrangements;
    row.numFixedRearrangements = result.numFixedRearrangements;
    row.numLateRearrangements = result.numLateRearrangements;
    row.initialQuality = result.initialQuality;
    row.finalQuality = result.finalQuality;
    row.maxQuality = result.maxQuality;
    row.initialAnisotropy = result.initialAnisotropy;
    row.finalAnisotropy = result.finalAnisotropy;
    row.maxAnisotropy = result.maxAnisotropy;
    row.sbfControlPointCount = result.sbfControlPointCount;
    row.sbfControlPointFraction = result.sbfControlPointFraction;
    row.sbfFillDistance = result.sbfFillDistance;
    row.transferStencilSize = result.transferStencilSize;
    row.predictedQualityCrossingTime = result.predictedQualityCrossingTime;
    row.actualQualityCrossingTime = result.actualQualityCrossingTime;
    row.predictorLeadSteps = result.predictorLeadSteps;
    row.elapsedSeconds = elapsed;
    fprintf('  relerr=%.6e qmax=%.3f amax=%.3f rearr=%d wall=%.2fs\n', ...
        row.relerr, row.maxQuality, row.maxAnisotropy, ...
        row.numRearrangements, row.elapsedSeconds);
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
figures.convergence = writeMetricFigures(T, imageDir, "relerr", ...
    '$\|e\|_2/\|q^\star\|_2$', 'convergence');
figures.timing = writeMetricFigures(T, imageDir, "elapsedSeconds", ...
    'wall time (s)', 'timing');
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
    title(sprintf('Spheroid SL rearrangement %s, \\xi = %d', tag, xi), ...
        'Interpreter', 'tex');
    paths(ixi) = fullfile(imageDir, ...
        sprintf('moving_surface_adr_tp_spheroid_rearrangement_%s_xi%d.png', tag, xi));
    exportgraphics(fig, paths(ixi), 'Resolution', 220);
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
        "relerr", '-o', 'SBF SL, M/N=0.5');
    plotRows(T, T.xi == xi & T.status == "ok" & T.mode == "rearranged" & ...
        T.transferMode == "localTp", "relerr", '-s', 'local TP SL');
    hold off;
    set(gca, 'XScale', 'log', 'YScale', 'log');
    grid on;
    xlabel('$\sqrt{N}$', 'Interpreter', 'latex');
    ylabel('$\|e\|_2/\|q^\star\|_2$', 'Interpreter', 'latex');
    legend('Location', 'best');
    title(sprintf('Pure Lagrangian vs SL rearranged, \\xi = %d', xi), ...
        'Interpreter', 'tex');
    paths(ixi) = fullfile(imageDir, ...
        sprintf('moving_surface_adr_tp_spheroid_rearrangement_comparison_xi%d.png', xi));
    exportgraphics(fig, paths(ixi), 'Resolution', 220);
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
    plotRows(T, mask, metricName, '-o', sprintf('SBF M/N = %.2g', mFactor), ...
        colors(colorIndex, :));
end
if any(base & T.transferMode == "localTp")
    colorIndex = colorIndex + 1;
    mask = base & T.transferMode == "localTp";
    plotRows(T, mask, metricName, '-s', 'local TP', colors(colorIndex, :));
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
hold on;
mFactors = T.MFactor(T.xi == xi & T.status == "ok" & T.mode == "rearranged" & ...
    T.transferMode == "sbf");
if isempty(mFactors)
    return;
end
bestM = max(mFactors);
mask = T.xi == xi & T.status == "ok" & T.mode == "rearranged" & ...
    T.transferMode == "sbf" & abs(T.MFactor - bestM) < 10 * eps;
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
    'transferMode', "", ...
    'useRearrangement', false, ...
    'h', NaN, ...
    'dt', NaN, ...
    'nsteps', 0, ...
    'relerr', NaN, ...
    'massRelError', NaN, ...
    'balanceRelResidual', NaN, ...
    'massCorrectionSteps', NaN, ...
    'massCorrectionMaxRelative', NaN, ...
    'massCorrectionMaxPointShift', NaN, ...
    'numRearrangements', 0, ...
    'numPredictiveRearrangements', 0, ...
    'numThresholdRearrangements', 0, ...
    'numFixedRearrangements', 0, ...
    'numLateRearrangements', 0, ...
    'initialQuality', NaN, ...
    'finalQuality', NaN, ...
    'maxQuality', NaN, ...
    'initialAnisotropy', NaN, ...
    'finalAnisotropy', NaN, ...
    'maxAnisotropy', NaN, ...
    'sbfControlPointCount', NaN, ...
    'sbfControlPointFraction', NaN, ...
    'sbfFillDistance', NaN, ...
    'transferStencilSize', NaN, ...
    'predictedQualityCrossingTime', NaN, ...
    'actualQualityCrossingTime', NaN, ...
    'predictorLeadSteps', NaN, ...
    'elapsedSeconds', NaN);
end
