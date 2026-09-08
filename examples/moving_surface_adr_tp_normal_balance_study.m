function study = moving_surface_adr_tp_normal_balance_study(varargin)
%MOVING_SURFACE_ADR_TP_NORMAL_BALANCE_STUDY Exact vs computed normals.
%   Runs the manufactured moving-surface ADR tests with exact analytic normals,
%   PCA point-cloud normals, and global parametric SBF normals.  The output
%   records PDE error, global mass/balance diagnostics, and normal-angle errors.

parser = inputParser();
parser.addParameter('Xi', 4, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('DtScale', 0.05, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('NormalNeighborCount', 32, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('GlobalSBFNormalDegree', 7, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('GlobalSBFControlPointCount', Inf, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('NormalModes', ["exact", "pca", "globalSbf"]);
parser.addParameter('OutputPrefix', fullfile(pwd, 'moving_surface_adr_tp_normal_balance_study_xi4'));
parser.addParameter('ImageDir', fullfile(pwd, 'docs', 'figures'));
parser.parse(varargin{:});
opts = parser.Results;
normalModes = string(opts.NormalModes);

if ~exist(opts.ImageDir, 'dir')
    mkdir(opts.ImageDir);
end

cases = localCases();
rows = repmat(emptyRow(), 0, 1);
rowIndex = 0;

for ic = 1:numel(cases)
    for mode = normalModes
        fprintf('\nNormal/balance study case=%s mode=%s\n', cases(ic).caseName, mode);
        run = moving_surface_adr_tp_geometric_flow_suite( ...
            cases(ic).caseName, opts.Xi, cases(ic).Nvals, opts.DtScale, ...
            'NormalMode', mode, ...
            'NormalNeighborCount', opts.NormalNeighborCount, ...
            'GlobalSBFNormalDegree', opts.GlobalSBFNormalDegree, ...
            'GlobalSBFControlPointCount', opts.GlobalSBFControlPointCount, ...
            'SpectrumCheck', false, ...
            'WriteOutputs', false);
        R = run.results;
        for j = 1:numel(R.N)
            rowIndex = rowIndex + 1;
            rows(rowIndex, 1) = makeRow(cases(ic).caseName, cases(ic).title, ...
                mode, R, j);
        end
    end
end

study = struct();
study.xi = opts.Xi;
study.normalNeighborCount = opts.NormalNeighborCount;
study.globalSBFNormalDegree = opts.GlobalSBFNormalDegree;
study.globalSBFControlPointCount = opts.GlobalSBFControlPointCount;
study.normalModes = normalModes;
study.rows = rows;
study.table = struct2table(rows);
study.figures = writeFigures(study.table, cases, opts.ImageDir, opts.Xi);

save(string(opts.OutputPrefix) + ".mat", 'study');
writetable(study.table, string(opts.OutputPrefix) + ".csv");
end

function cases = localCases()
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
    'Nvals', { ...
        [576, 1024, 1600], ...
        [576, 1024, 1600], ...
        [576, 1024, 1600], ...
        [576, 1024, 1600], ...
        [1600, 2304, 3136]});
end

function row = makeRow(caseName, titleText, mode, R, j)
row = emptyRow();
row.caseName = string(caseName);
row.title = string(titleText);
row.normalMode = string(mode);
row.xi = R.xi;
row.N = R.N(j);
row.sqrtN = sqrt(R.N(j));
row.h = R.h(j);
row.dt = R.dt(j);
row.nsteps = R.nsteps(j);
row.relerr = R.relerr(j);
row.massRelError = R.massRelError(j);
row.balanceRelResidual = R.balanceRelResidual(j);
row.normalRMSAngle = R.normalRMSAngle(j);
row.normalMaxAngle = R.normalMaxAngle(j);
end

function figures = writeFigures(T, cases, imageDir, xi)
figures = struct('error', strings(1, numel(cases)), 'normal', strings(1, numel(cases)));
for ic = 1:numel(cases)
    idx = T.caseName == cases(ic).caseName;
    Tc = T(idx, :);
    figures.error(ic) = writeErrorFigure(Tc, cases(ic), imageDir, xi);
    figures.normal(ic) = writeNormalFigure(Tc, cases(ic), imageDir, xi);
end
end

function imagePath = writeErrorFigure(T, caseInfo, imageDir, xi)
fig = figure('Color', 'w', 'Position', [100, 100, 760, 560]);
hold on;
plotMode(T, "exact", 'o-', 'exact normals');
plotMode(T, "pca", 's-', 'PCA normals');
plotMode(T, "globalSbf", '^-', 'global SBF normals');
hold off;
set(gca, 'XScale', 'log', 'YScale', 'log');
grid on;
xlabel('$\sqrt{N}$', 'Interpreter', 'latex');
ylabel('$\|e\|_2/\|c\|_2$', 'Interpreter', 'latex');
title(sprintf('%s, normal comparison, \\xi = %d', caseInfo.title, xi), ...
    'Interpreter', 'tex');
legend('Location', 'southwest');
imagePath = fullfile(imageDir, ...
    sprintf('moving_surface_adr_tp_normal_error_sqrtN_%s_xi%d.png', caseInfo.caseName, xi));
exportgraphics(fig, imagePath, 'Resolution', 220);
close(fig);
end

function imagePath = writeNormalFigure(T, caseInfo, imageDir, xi)
fig = figure('Color', 'w', 'Position', [100, 100, 760, 560]);
hold on;
plotNormalAngles(T, "pca", 'o-', 's-', 'PCA');
plotNormalAngles(T, "globalSbf", '^-', 'd-', 'global SBF');
hold off;
set(gca, 'XScale', 'log', 'YScale', 'log');
grid on;
xlabel('$\sqrt{N}$', 'Interpreter', 'latex');
ylabel('Normal error (degrees)');
title(sprintf('%s, computed normal error, \\xi = %d', caseInfo.title, xi), ...
    'Interpreter', 'tex');
legend('Location', 'northeast');
imagePath = fullfile(imageDir, ...
    sprintf('moving_surface_adr_tp_normal_angle_sqrtN_%s_xi%d.png', caseInfo.caseName, xi));
exportgraphics(fig, imagePath, 'Resolution', 220);
close(fig);
end

function plotNormalAngles(T, mode, rmsStyle, maxStyle, labelText)
Tm = T(T.normalMode == mode, :);
if isempty(Tm)
    return;
end
[~, order] = sort(Tm.N);
Tm = Tm(order, :);
loglog(sqrt(Tm.N), Tm.normalRMSAngle, rmsStyle, ...
    'LineWidth', 1.5, 'MarkerSize', 6, ...
    'DisplayName', sprintf('%s RMS angle', labelText));
loglog(sqrt(Tm.N), Tm.normalMaxAngle, maxStyle, ...
    'LineWidth', 1.5, 'MarkerSize', 6, ...
    'DisplayName', sprintf('%s max angle', labelText));
end

function plotMode(T, mode, style, labelText)
Tm = T(T.normalMode == mode, :);
[~, order] = sort(Tm.N);
Tm = Tm(order, :);
loglog(sqrt(Tm.N), Tm.relerr, style, ...
    'LineWidth', 1.5, ...
    'MarkerSize', 6, ...
    'DisplayName', labelText);
end

function row = emptyRow()
row = struct( ...
    'caseName', "", ...
    'title', "", ...
    'normalMode', "", ...
    'xi', 0, ...
    'N', 0, ...
    'sqrtN', 0, ...
    'h', NaN, ...
    'dt', NaN, ...
    'nsteps', 0, ...
    'relerr', NaN, ...
    'massRelError', NaN, ...
    'balanceRelResidual', NaN, ...
    'normalRMSAngle', NaN, ...
    'normalMaxAngle', NaN);
end
