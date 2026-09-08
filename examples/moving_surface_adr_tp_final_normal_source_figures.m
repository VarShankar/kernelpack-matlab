function figures = moving_surface_adr_tp_final_normal_source_figures(varargin)
%MOVING_SURFACE_ADR_TP_FINAL_NORMAL_SOURCE_FIGURES Plot normal-source study.
%   Reads the saved normal-source benchmark tables and writes
%   focused paper plots for accuracy, per-step runtime, and computed-normal error.

parser = inputParser();
defaultDataFiles = [ ...
    string(fullfile(pwd, 'moving_surface_adr_tp_final_normal_source_benchmark_anisotropic_ellipsoid.csv')), ...
    string(fullfile(pwd, 'moving_surface_adr_tp_final_normal_source_benchmark_breathing_torus.csv'))];
parser.addParameter('DataFiles', defaultDataFiles);
parser.addParameter('ImageDir', fullfile(pwd, 'docs', 'figures'));
parser.addParameter('CombinedCsv', fullfile(pwd, ...
    'moving_surface_adr_tp_final_normal_source_benchmark.csv'));
parser.addParameter('WriteAuxiliaryFigures', false, @(x) islogical(x) && isscalar(x));
parser.parse(varargin{:});
opts = parser.Results;

if ~exist(opts.ImageDir, 'dir')
    mkdir(opts.ImageDir);
end

T = readBenchmarkTables(string(opts.DataFiles));
T = ensureTimeStepMetadata(T, 0.05);
writetable(T, opts.CombinedCsv);

figures = struct();
figures.combinedCsv = string(opts.CombinedCsv);
figures.error = writeMetricFigures(T, opts.ImageDir, "relerr", ...
    '\textbf{relative $\ell_2$ error}', "error");
figures.timing = writeMetricFigures(T, opts.ImageDir, "elapsedPerStep", ...
    "Average runtime per timestep (s)", "timing");
figures.normalAngle = writeNormalAngleFigures(T, opts.ImageDir);
if opts.WriteAuxiliaryFigures
    figures.sbfFraction = writeSbfFractionFigures(T, opts.ImageDir);
else
    figures.sbfFraction = strings(0, 0);
end
end

function T = readBenchmarkTables(dataFiles)
tables = cell(numel(dataFiles), 1);
for k = 1:numel(dataFiles)
    tables{k} = readtable(dataFiles(k), 'TextType', 'string');
end
T = vertcat(tables{:});
T = T(T.status == "ok", :);
T.caseName = string(T.caseName);
T.normalMode = string(T.normalMode);
T.title = string(T.title);
end

function T = ensureTimeStepMetadata(T, defaultDtScale)
if ~ismember("finalTime", string(T.Properties.VariableNames))
    T.finalTime = nan(height(T), 1);
end
if ~ismember("dt", string(T.Properties.VariableNames))
    T.dt = nan(height(T), 1);
end
if ~ismember("nsteps", string(T.Properties.VariableNames))
    T.nsteps = nan(height(T), 1);
end

for i = 1:height(T)
    if ~isfinite(T.finalTime(i)) || T.finalTime(i) <= 0 || ...
            ~isfinite(T.dt(i)) || T.dt(i) <= 0 || ...
            ~isfinite(T.nsteps(i)) || T.nsteps(i) <= 0
        problem = moving_surface_adr_tp_problem(T.caseName(i));
        geom0 = problem.geometry(T.N(i), 0.0);
        dtTarget = defaultDtScale * geom0.h^(T.xi(i) / 3);
        nsteps = max(3, ceil(problem.finalTime / dtTarget));
        T.finalTime(i) = problem.finalTime;
        T.nsteps(i) = nsteps;
        T.dt(i) = problem.finalTime / nsteps;
    end
end
T.elapsedPerStep = T.elapsedSeconds ./ max(T.nsteps, 1);
end

function imagePaths = writeMetricFigures(T, imageDir, metricName, yLabelText, suffix)
cases = unique(T.caseName(:), 'stable').';
xiVals = unique(T.xi(:).');
imagePaths = strings(numel(cases), numel(xiVals));

for ic = 1:numel(cases)
    caseName = cases(ic);
    for ix = 1:numel(xiVals)
        xi = xiVals(ix);
        fig = figure('Color', 'w', 'Position', [100, 100, 760, 560]);
        idx = T.caseName == caseName & T.xi == xi;
        plotNormalModeMetric(T(idx, :), metricName);
        title('Impact of normals', 'Interpreter', 'none');
        xlabel('\textbf{$\sqrt{N}$}', 'Interpreter', 'latex');
        if contains(string(yLabelText), "\textbf")
            ylabel(char(yLabelText), 'Interpreter', 'latex');
        else
            ylabel(char(yLabelText), 'Interpreter', 'none');
        end
        imagePaths(ic, ix) = fullfile(imageDir, sprintf( ...
            'moving_surface_adr_tp_final_normal_%s_%s_xi%d.png', ...
            suffix, caseName, xi));
        kp.plot.exportPaperFigure(fig, imagePaths(ic, ix));
        close(fig);
    end
end
end

function plotNormalModeMetric(T, metricName)
hold on;
styles = normalModeStyles();
for mode = ["exact", "pca", "calibratedSbf"]
    idx = T.normalMode == mode;
    if ~any(idx)
        continue;
    end
    S = sortrows(T(idx, :), "N");
    style = styles.(mode);
    loglog(sqrt(S.N), S.(metricName), style.line, ...
        'LineWidth', 1.6, ...
        'MarkerSize', 6, ...
        'DisplayName', style.label);
end
hold off;
set(gca, 'XScale', 'log', 'YScale', 'log');
grid on;
legend('Location', 'best');
end

function imagePaths = writeNormalAngleFigures(T, imageDir)
cases = unique(T.caseName(:), 'stable').';
xiVals = unique(T.xi(:).');
imagePaths = strings(numel(cases), numel(xiVals));

for ic = 1:numel(cases)
    caseName = cases(ic);
    for ix = 1:numel(xiVals)
        xi = xiVals(ix);
        fig = figure('Color', 'w', 'Position', [100, 100, 760, 560]);
        idx = T.caseName == caseName & T.xi == xi & T.normalMode ~= "exact";
        plotNormalModeMetric(T(idx, :), "normalRMSAngle");
        title('Impact of normals', 'Interpreter', 'none');
        xlabel('\textbf{$\sqrt{N}$}', 'Interpreter', 'latex');
        ylabel('RMS normal angle (degrees)');
        imagePaths(ic, ix) = fullfile(imageDir, sprintf( ...
            'moving_surface_adr_tp_final_normal_angle_%s_xi%d.png', caseName, xi));
        kp.plot.exportPaperFigure(fig, imagePaths(ic, ix));
        close(fig);
    end
end
end

function imagePaths = writeSbfFractionFigures(T, imageDir)
cases = unique(T.caseName(:), 'stable').';
xiVals = unique(T.xi(:).');
imagePaths = strings(numel(cases), numel(xiVals));

for ic = 1:numel(cases)
    caseName = cases(ic);
    for ix = 1:numel(xiVals)
        xi = xiVals(ix);
        S = T(T.caseName == caseName & T.xi == xi & T.normalMode == "calibratedSbf", :);
        S = sortrows(S, "N");
        fig = figure('Color', 'w', 'Position', [100, 100, 760, 560]);
        yyaxis left;
        loglog(sqrt(S.N), S.normalControlPointCount, '-o', ...
            'LineWidth', 1.6, 'MarkerSize', 6, 'DisplayName', '$M$');
        ylabel('SBF control sites', 'Interpreter', 'none');
        yyaxis right;
        loglog(sqrt(S.N), S.normalControlPointFraction, '-s', ...
            'LineWidth', 1.6, 'MarkerSize', 6, 'DisplayName', 'control-site fraction');
        ylabel('SBF control-site fraction', 'Interpreter', 'none');
        set(gca, 'XScale', 'log', 'YScale', 'log');
        grid on;
        xlabel('\textbf{$\sqrt{N}$}', 'Interpreter', 'latex');
        legend('Location', 'best', 'Interpreter', 'latex');
        imagePaths(ic, ix) = fullfile(imageDir, sprintf( ...
            'moving_surface_adr_tp_final_normal_sbf_fraction_%s_xi%d.png', caseName, xi));
        kp.plot.exportPaperFigure(fig, imagePaths(ic, ix));
        close(fig);
    end
end
end

function styles = normalModeStyles()
styles = struct();
styles.exact = struct('line', '-o', 'label', 'exact');
styles.pca = struct('line', '-s', 'label', 'PCA');
styles.calibratedSbf = struct('line', '-^', 'label', 'SBF');
end
