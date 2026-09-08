function figures = moving_surface_adr_tp_results_figures(varargin)
%MOVING_SURFACE_ADR_TP_RESULTS_FIGURES Generate paper figures for moving ADR.
%   The convergence plots use archived result files from the moving-surface
%   geometric-flow suite.  The timing plots use the local-defect,
%   hyperviscosity-update, global-solver, and combined-acceleration benchmark
%   tables.

parser = inputParser();
parser.addParameter('RunTimingBenchmark', false, @(x) islogical(x) && isscalar(x));
parser.addParameter('TimingXiVals', 4, @(x) isnumeric(x) && isvector(x));
parser.addParameter('TimingNVals', paperNVals(), @(x) isnumeric(x) && isvector(x));
parser.addParameter('TimingDtScale', 0.05, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('TimingHyperviscosityUpdateMode', "adaptiveGeometry");
parser.addParameter('TimingPrefix', fullfile(pwd, 'moving_surface_adr_tp_defect_benchmark_results_xi4'));
parser.addParameter('HyperviscosityTimingPrefix', string.empty(1, 0));
parser.addParameter('GlobalSolverTimingPrefix', string.empty(1, 0));
parser.addParameter('CombinedTimingPrefix', string.empty(1, 0));
parser.addParameter('ImageDir', fullfile(pwd, 'docs', 'figures'));
parser.parse(varargin{:});
opts = parser.Results;

if ~exist(opts.ImageDir, 'dir')
    mkdir(opts.ImageDir);
end

function Nvals = paperNVals()
Nvals = [576, 1024, 1600, 2500, 3600, 4900];
end

cases = movingSurfaceCases();
figures = struct();
figures.convergence = writeConvergenceSqrtNFigures(cases, opts.ImageDir);

if opts.RunTimingBenchmark
    moving_surface_adr_tp_defect_benchmark( ...
        [cases.caseName], opts.TimingXiVals, opts.TimingNVals, opts.TimingDtScale, ...
        'OutputPrefix', opts.TimingPrefix, ...
        'HyperviscosityUpdateMode', opts.TimingHyperviscosityUpdateMode, ...
        'SpectrumCheck', false, ...
        'WarmupParallelPool', true, ...
        'ContinueOnError', true, ...
        'SaveAfterEach', true);
end

timingPrefixes = string(opts.TimingPrefix);
timingMats = timingPrefixes + ".mat";
timingMats = timingMats(isfile(timingMats));
if ~isempty(timingMats)
    figures.timing = strings(0, 0);
    figures.speedup = strings(0, 0);
    for imat = 1:numel(timingMats)
        [timingPaths, speedupPaths] = writeTimingAndSpeedupForBenchmark(timingMats(imat), cases, opts.ImageDir);
        figures.timing = [figures.timing; timingPaths(:)];
        figures.speedup = [figures.speedup; speedupPaths(:)];
    end
else
    figures.timing = "";
    figures.speedup = "";
    if ~isempty(timingPrefixes)
        warning('kp:examples:MissingTimingResults', ...
            'Timing benchmark file not found for prefix: %s', strjoin(timingPrefixes, ", "));
    end
end

hyperPrefixes = string(opts.HyperviscosityTimingPrefix);
hyperMats = hyperPrefixes + ".mat";
hyperMats = hyperMats(isfile(hyperMats));
if ~isempty(hyperMats)
    figures.hyperviscosityTiming = strings(0, 0);
    figures.hyperviscositySpeedup = strings(0, 0);
    for imat = 1:numel(hyperMats)
        [timingPaths, speedupPaths] = writeHyperviscosityTimingAndSpeedupForBenchmark( ...
            hyperMats(imat), cases, opts.ImageDir);
        figures.hyperviscosityTiming = [figures.hyperviscosityTiming; timingPaths(:)];
        figures.hyperviscositySpeedup = [figures.hyperviscositySpeedup; speedupPaths(:)];
    end
else
    figures.hyperviscosityTiming = "";
    figures.hyperviscositySpeedup = "";
end

globalSolverPrefixes = string(opts.GlobalSolverTimingPrefix);
globalSolverMats = globalSolverPrefixes + ".mat";
globalSolverMats = globalSolverMats(isfile(globalSolverMats));
if ~isempty(globalSolverMats)
    figures.globalSolverTiming = strings(0, 0);
    figures.globalSolverSpeedup = strings(0, 0);
    for imat = 1:numel(globalSolverMats)
        [timingPaths, speedupPaths] = writeGlobalSolverTimingAndSpeedupForBenchmark( ...
            globalSolverMats(imat), cases, opts.ImageDir);
        figures.globalSolverTiming = [figures.globalSolverTiming; timingPaths(:)];
        figures.globalSolverSpeedup = [figures.globalSolverSpeedup; speedupPaths(:)];
    end
else
    figures.globalSolverTiming = "";
    figures.globalSolverSpeedup = "";
end

combinedPrefixes = string(opts.CombinedTimingPrefix);
combinedMats = combinedPrefixes + ".mat";
combinedMats = combinedMats(isfile(combinedMats));
if ~isempty(combinedMats)
    figures.combinedTiming = strings(0, 0);
    figures.combinedSpeedup = strings(0, 0);
    for imat = 1:numel(combinedMats)
        [timingPaths, speedupPaths] = writeCombinedTimingAndSpeedupForBenchmark( ...
            combinedMats(imat), cases, opts.ImageDir);
        figures.combinedTiming = [figures.combinedTiming; timingPaths(:)];
        figures.combinedSpeedup = [figures.combinedSpeedup; speedupPaths(:)];
    end
else
    figures.combinedTiming = "";
    figures.combinedSpeedup = "";
end
end

function cases = movingSurfaceCases()
cases = struct( ...
    'caseName', { ...
        "mcf_sphere", ...
        "imcf_sphere", ...
        "rotating_breathing_sphere", ...
        "anisotropic_ellipsoid", ...
        "breathing_torus"}, ...
    'title', { ...
        "MCF sphere", ...
        "IMCF sphere", ...
        "Rotating breathing sphere", ...
        "Anisotropic ellipsoid", ...
        "Breathing torus"}, ...
    'resultsFile', { ...
        "moving_surface_adr_tp_convergence_mcf_sphere_results.csv", ...
        "moving_surface_adr_tp_convergence_imcf_sphere_results.csv", ...
        "moving_surface_adr_tp_convergence_rotating_breathing_sphere_results.csv", ...
        "moving_surface_adr_tp_convergence_anisotropic_ellipsoid_results.csv", ...
        "moving_surface_adr_tp_convergence_breathing_torus_results.csv"});
end

function imagePaths = writeConvergenceSqrtNFigures(cases, imageDir)
imagePaths = strings(1, numel(cases));
colors = lines(6);

for ic = 1:numel(cases)
    fig = figure('Color', 'w', 'Position', [100, 100, 760, 560]);
    results = readConvergenceResults(cases(ic).resultsFile);
    hold on;
    for ir = 1:numel(results)
        mask = convergencePlotMask(cases(ic).caseName, results(ir));
        x = sqrt(results(ir).N(mask));
        y = results(ir).relerr(mask);
        loglog(x, y, '-o', ...
            'LineWidth', 1.5, ...
            'MarkerSize', 6, ...
            'Color', colors(ir, :), ...
            'DisplayName', sprintf('\\xi = %d', results(ir).xi));
    end
    hold off;
    set(gca, 'XScale', 'log', 'YScale', 'log');
    grid on;
    xlabel('\textbf{$\sqrt{N}$}', 'Interpreter', 'latex');
    ylabel('\textbf{relative $\ell_2$ error}', 'Interpreter', 'latex');
    title(sprintf('%s convergence', cases(ic).title), 'Interpreter', 'none');
    legend('Location', 'best');
    imagePaths(ic) = fullfile(imageDir, ...
        sprintf('moving_surface_adr_tp_convergence_sqrtN_%s.png', cases(ic).caseName));
    kp.plot.exportPaperFigure(fig, imagePaths(ic));
    close(fig);
end
end

function results = readConvergenceResults(resultsFile)
T = readtable(resultsFile);
xiVals = unique(T.xi, 'stable');
results = repmat(struct('xi', 0, 'N', [], 'relerr', []), 1, numel(xiVals));
for k = 1:numel(xiVals)
    rows = T.xi == xiVals(k);
    results(k).xi = xiVals(k);
    results(k).N = T.N(rows).';
    results(k).relerr = T.relerr(rows).';
end
end

function mask = convergencePlotMask(caseName, result)
mask = true(size(result.N));

% The coarse high-order torus runs are intentionally retained in the raw
% output files, but they are outside the asymptotic/stable regime and make the
% publication plot unreadable.
if string(caseName) == "breathing_torus"
    if result.xi == 4
        mask = result.N >= 1024;
    elseif result.xi == 6
        mask = result.N >= 1600;
    end
end
end

function [timingPaths, speedupPaths] = writeTimingAndSpeedupForBenchmark(timingMat, cases, imageDir)
S = load(timingMat);
T = S.bench.table;
presentCases = unique(string(T.caseName(:))).';
caseNames = [cases.caseName];
cases = cases(ismember(caseNames, presentCases));
timingPaths = writeTimingSqrtNFiguresFromTable(T, cases, imageDir);
speedupPaths = writeSpeedupSqrtNFiguresFromTable(T, cases, imageDir);
end

function imagePaths = writeTimingSqrtNFiguresFromTable(T, cases, imageDir)
imagePaths = strings(1, numel(cases));

for ic = 1:numel(cases)
    fig = figure('Color', 'w', 'Position', [100, 100, 820, 600]);
    idx = T.caseName == cases(ic).caseName & T.status == "ok";
    plotTimingComparisonByXi(T(idx, :), ...
        sprintf('%s local stencil-update timing', cases(ic).title), ...
        "directTime", "defectTime", ...
        "direct update", "defect-corrected update", ...
        'Average update time per timestep (s)');
    imagePaths(ic) = fullfile(imageDir, ...
        sprintf('moving_surface_adr_tp_timing_sqrtN_%s.png', cases(ic).caseName));
    kp.plot.exportPaperFigure(fig, imagePaths(ic));
    close(fig);
end
end

function imagePaths = writeSpeedupSqrtNFiguresFromTable(T, cases, imageDir)
imagePaths = strings(1, numel(cases));

for ic = 1:numel(cases)
    fig = figure('Color', 'w', 'Position', [100, 100, 820, 600]);
    idx = T.caseName == cases(ic).caseName & T.status == "ok";
    plotSpeedupByXi(T(idx, :), ...
        sprintf('%s local stencil-update speedup', cases(ic).title), ...
        "directTime", "defectTime", ...
        'Direct/defect-corrected per-step time');
    imagePaths(ic) = fullfile(imageDir, ...
        sprintf('moving_surface_adr_tp_speedup_sqrtN_%s.png', cases(ic).caseName));
    kp.plot.exportPaperFigure(fig, imagePaths(ic));
    close(fig);
end
end

function [timingPaths, speedupPaths] = writeHyperviscosityTimingAndSpeedupForBenchmark(timingMat, cases, imageDir)
S = load(timingMat);
T = S.bench.table;
presentCases = unique(string(T.caseName(:))).';
caseNames = [cases.caseName];
cases = cases(ismember(caseNames, presentCases));
timingPaths = writeHyperviscosityTimingSqrtNFiguresFromTable(T, cases, imageDir);
speedupPaths = writeHyperviscositySpeedupSqrtNFiguresFromTable(T, cases, imageDir);
end

function imagePaths = writeHyperviscosityTimingSqrtNFiguresFromTable(T, cases, imageDir)
imagePaths = strings(1, numel(cases));

for ic = 1:numel(cases)
    fig = figure('Color', 'w', 'Position', [100, 100, 820, 600]);
    idx = T.caseName == cases(ic).caseName & T.status == "ok";
    plotTimingComparisonByXi(T(idx, :), ...
        sprintf('%s hyperviscosity-update timing', cases(ic).title), ...
        "everyStepTime", "adaptiveTime", ...
        "recompute every step", "spectral predictor", ...
        'Average hyperviscosity-update time per timestep (s)');
    imagePaths(ic) = fullfile(imageDir, ...
        sprintf('moving_surface_adr_tp_hyperviscosity_timing_sqrtN_%s.png', ...
        cases(ic).caseName));
    kp.plot.exportPaperFigure(fig, imagePaths(ic));
    close(fig);
end
end

function imagePaths = writeHyperviscositySpeedupSqrtNFiguresFromTable(T, cases, imageDir)
imagePaths = strings(1, numel(cases));

for ic = 1:numel(cases)
    fig = figure('Color', 'w', 'Position', [100, 100, 820, 600]);
    idx = T.caseName == cases(ic).caseName & T.status == "ok";
    plotSpeedupByXi(T(idx, :), ...
        sprintf('%s hyperviscosity-update speedup', cases(ic).title), ...
        "everyStepTime", "adaptiveTime", ...
        'Every-step/predictor per-step time');
    imagePaths(ic) = fullfile(imageDir, ...
        sprintf('moving_surface_adr_tp_hyperviscosity_speedup_sqrtN_%s.png', ...
        cases(ic).caseName));
    kp.plot.exportPaperFigure(fig, imagePaths(ic));
    close(fig);
end
end

function [timingPaths, speedupPaths] = writeGlobalSolverTimingAndSpeedupForBenchmark(timingMat, cases, imageDir)
S = load(timingMat);
T = S.bench.table;
presentCases = unique(string(T.caseName(:))).';
caseNames = [cases.caseName];
cases = cases(ismember(caseNames, presentCases));
timingPaths = writeGlobalSolverTimingSqrtNFiguresFromTable(T, cases, imageDir);
speedupPaths = writeGlobalSolverSpeedupSqrtNFiguresFromTable(T, cases, imageDir);
end

function imagePaths = writeGlobalSolverTimingSqrtNFiguresFromTable(T, cases, imageDir)
imagePaths = strings(1, numel(cases));

for ic = 1:numel(cases)
    fig = figure('Color', 'w', 'Position', [100, 100, 820, 600]);
    idx = T.caseName == cases(ic).caseName & T.status == "ok";
    plotTimingComparisonByXi(T(idx, :), ...
        sprintf('%s global linear-solve timing', cases(ic).title), ...
        "noDefectSolveTime", "defectSolveTime", ...
        "GMRES+ILU", "global defect correction", ...
        'Average global-solve time per timestep (s)');
    imagePaths(ic) = fullfile(imageDir, ...
        sprintf('moving_surface_adr_tp_global_solver_timing_sqrtN_%s.png', ...
        cases(ic).caseName));
    kp.plot.exportPaperFigure(fig, imagePaths(ic));
    close(fig);
end
end

function imagePaths = writeGlobalSolverSpeedupSqrtNFiguresFromTable(T, cases, imageDir)
imagePaths = strings(1, numel(cases));

for ic = 1:numel(cases)
    fig = figure('Color', 'w', 'Position', [100, 100, 820, 600]);
    idx = T.caseName == cases(ic).caseName & T.status == "ok";
    plotSpeedupByXi(T(idx, :), ...
        sprintf('%s global linear-solve speedup', cases(ic).title), ...
        "noDefectSolveTime", "defectSolveTime", ...
        'GMRES+ILU/global-defect per-step solve time');
    imagePaths(ic) = fullfile(imageDir, ...
        sprintf('moving_surface_adr_tp_global_solver_speedup_sqrtN_%s.png', ...
        cases(ic).caseName));
    kp.plot.exportPaperFigure(fig, imagePaths(ic));
    close(fig);
end
end

function [timingPaths, speedupPaths] = writeCombinedTimingAndSpeedupForBenchmark(timingMat, cases, imageDir)
S = load(timingMat);
T = S.bench.table;
presentCases = unique(string(T.caseName(:))).';
caseNames = [cases.caseName];
cases = cases(ismember(caseNames, presentCases));
timingPaths = writeCombinedTimingSqrtNFiguresFromTable(T, cases, imageDir);
speedupPaths = writeCombinedSpeedupSqrtNFiguresFromTable(T, cases, imageDir);
end

function imagePaths = writeCombinedTimingSqrtNFiguresFromTable(T, cases, imageDir)
imagePaths = strings(1, numel(cases));

for ic = 1:numel(cases)
    fig = figure('Color', 'w', 'Position', [100, 100, 820, 600]);
    idx = T.caseName == cases(ic).caseName & T.status == "ok";
    plotTimingComparisonByXi(T(idx, :), ...
        sprintf('%s aggregate timestep timing', cases(ic).title), ...
        "baselineTime", "acceleratedTime", ...
        "direct rebuilds", "update strategy", ...
        'Average runtime per timestep (s)');
    imagePaths(ic) = fullfile(imageDir, ...
        sprintf('moving_surface_adr_tp_combined_timing_sqrtN_%s.png', ...
        cases(ic).caseName));
    kp.plot.exportPaperFigure(fig, imagePaths(ic));
    close(fig);
end
end

function imagePaths = writeCombinedSpeedupSqrtNFiguresFromTable(T, cases, imageDir)
imagePaths = strings(1, numel(cases));

for ic = 1:numel(cases)
    fig = figure('Color', 'w', 'Position', [100, 100, 820, 600]);
    idx = T.caseName == cases(ic).caseName & T.status == "ok";
    plotSpeedupByXi(T(idx, :), ...
        sprintf('%s aggregate timestep speedup', cases(ic).title), ...
        "baselineTime", "acceleratedTime", ...
        'Direct-rebuild/update-strategy per-step time');
    imagePaths(ic) = fullfile(imageDir, ...
        sprintf('moving_surface_adr_tp_combined_speedup_sqrtN_%s.png', ...
        cases(ic).caseName));
    kp.plot.exportPaperFigure(fig, imagePaths(ic));
    close(fig);
end
end

function plotTimingComparisonByXi(T, titleText, referenceColumn, updateColumn, ...
    referenceLabel, updateLabel, yLabelText)
if isempty(T)
    text(0.5, 0.5, 'no successful runs', 'HorizontalAlignment', 'center');
    axis off;
    return;
end

xiVals = unique(T.xi(:).');
colors = lines(max(numel(xiVals), 1));
hold on;
for ixi = 1:numel(xiVals)
    xi = xiVals(ixi);
    S = sortedSuccessfulRows(T, xi, referenceColumn, updateColumn);
    if isempty(S)
        continue;
    end
    x = sqrt(S.N);
    yReference = perStepValue(S.(referenceColumn), S.nsteps);
    yUpdate = perStepValue(S.(updateColumn), S.nsteps);
    loglog(x, yReference, '-o', ...
        'LineWidth', 1.6, ...
        'MarkerSize', 6, ...
        'Color', colors(ixi, :), ...
        'DisplayName', sprintf('\\xi=%d %s', xi, referenceLabel));
    loglog(x, yUpdate, '--s', ...
        'LineWidth', 1.6, ...
        'MarkerSize', 6, ...
        'Color', colors(ixi, :), ...
        'DisplayName', sprintf('\\xi=%d %s', xi, updateLabel));
end
hold off;
set(gca, 'XScale', 'log', 'YScale', 'log');
grid on;
xlabel('\textbf{$\sqrt{N}$}', 'Interpreter', 'latex');
ylabel(yLabelText, 'Interpreter', 'none');
title(titleText, 'Interpreter', 'none');
legend('Location', 'best');
end

function plotSpeedupByXi(T, titleText, referenceColumn, updateColumn, yLabelText)
if isempty(T)
    text(0.5, 0.5, 'no successful runs', 'HorizontalAlignment', 'center');
    axis off;
    return;
end

xiVals = unique(T.xi(:).');
colors = lines(max(numel(xiVals), 1));
hold on;
for ixi = 1:numel(xiVals)
    xi = xiVals(ixi);
    S = sortedSuccessfulRows(T, xi, referenceColumn, updateColumn);
    if isempty(S)
        continue;
    end
    yReference = perStepValue(S.(referenceColumn), S.nsteps);
    yUpdate = perStepValue(S.(updateColumn), S.nsteps);
    loglog(sqrt(S.N), yReference ./ max(yUpdate, eps), '-o', ...
        'LineWidth', 1.6, ...
        'MarkerSize', 6, ...
        'Color', colors(ixi, :), ...
        'DisplayName', sprintf('\\xi=%d', xi));
end
hold off;
set(gca, 'XScale', 'log', 'YScale', 'log');
grid on;
yline(1, '--k', 'LineWidth', 1.0, 'HandleVisibility', 'off');
xlabel('\textbf{$\sqrt{N}$}', 'Interpreter', 'latex');
ylabel(yLabelText, 'Interpreter', 'none');
title(titleText, 'Interpreter', 'none');
legend('Location', 'best');
end

function S = sortedSuccessfulRows(T, xi, referenceColumn, updateColumn)
S = T(T.xi == xi, :);
if isempty(S)
    return;
end
valid = isfinite(S.(referenceColumn)) & isfinite(S.(updateColumn)) & ...
    isfinite(S.nsteps) & S.nsteps > 0;
S = S(valid, :);
if isempty(S)
    return;
end
S = sortrows(S, "N");
end

function y = perStepValue(totalTime, nsteps)
y = totalTime ./ max(nsteps, 1);
end
