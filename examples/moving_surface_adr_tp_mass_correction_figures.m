function moving_surface_adr_tp_mass_correction_figures(varargin)
%MOVING_SURFACE_ADR_TP_MASS_CORRECTION_FIGURES Plot mass-fixer diagnostics.
%   Reads the CSV produced by moving_surface_adr_tp_mass_correction_benchmark
%   and writes one publication figure per challenging moving-surface problem.

parser = inputParser();
parser.addParameter('ResultsPath', fullfile(pwd, ...
    'moving_surface_adr_tp_mass_correction_benchmark.csv'));
parser.addParameter('FigureDir', fullfile(pwd, 'docs', 'figures'));
parser.parse(varargin{:});
opts = parser.Results;

T = readtable(opts.ResultsPath, 'TextType', 'string');
if ~exist(opts.FigureDir, 'dir')
    mkdir(opts.FigureDir);
end

caseNames = ["anisotropic_ellipsoid", "breathing_torus"];
for k = 1:numel(caseNames)
    plotCase(T(T.caseName == caseNames(k), :), ...
        fullfile(opts.FigureDir, "moving_surface_adr_tp_mass_correction_" + ...
        caseNames(k) + ".png"));
end
end

function plotCase(T, outputPath)
xiVals = unique(T.xi).';
colors = lines(numel(xiVals));
fig = figure('Color', 'w', 'Position', [100, 100, 1050, 420]);
tiledlayout(fig, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
sgtitle(fig, sprintf('%s mass correction', readableCaseName(T.caseName(1))), ...
    'FontSize', 16, 'FontWeight', 'bold');

nexttile;
hold on;
for k = 1:numel(xiVals)
    xi = xiVals(k);
    plotMode(T, xi, "off", "relerr", colors(k, :), '--o');
    plotMode(T, xi, "balance", "relerr", colors(k, :), '-s');
end
grid on;
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('\textbf{$\sqrt{N}$}', 'Interpreter', 'latex');
ylabel('\textbf{relative $\ell_2$ error}', 'Interpreter', 'latex');
title('Solution error', 'Interpreter', 'none');

nexttile;
hold on;
for k = 1:numel(xiVals)
    xi = xiVals(k);
    plotMode(T, xi, "off", "balanceRelResidual", colors(k, :), '--o');
    plotMode(T, xi, "balance", "balanceRelResidual", colors(k, :), '-s');
end
grid on;
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('\textbf{$\sqrt{N}$}', 'Interpreter', 'latex');
ylabel('relative residual', 'Interpreter', 'none');
title('Mass residual', 'Interpreter', 'none');

legendStrings = strings(1, 2 * numel(xiVals));
for k = 1:numel(xiVals)
    legendStrings(2 * k - 1) = sprintf('$\\xi=%d$, off', xiVals(k));
    legendStrings(2 * k) = sprintf('$\\xi=%d$, corrected', xiVals(k));
end
legend(legendStrings, 'Interpreter', 'latex', 'Location', 'best');

kp.plot.exportPaperFigure(fig, outputPath);
close(fig);
end

function label = readableCaseName(caseName)
label = char(strrep(string(caseName), "_", " "));
label = regexprep(label, '(^|\s).', '${upper($0)}');
end

function plotMode(T, xi, modeName, valueName, color, style)
S = T(T.xi == xi & T.massCorrectionMode == modeName, :);
S = sortrows(S, 'sqrtN');
x = S.sqrtN;
y = max(S.(valueName), eps);
loglog(x, y, style, 'Color', color, 'LineWidth', 1.4, ...
    'MarkerSize', 6, 'MarkerFaceColor', 'w');
end
