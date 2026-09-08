function study = moving_surface_adr_tp_long_time_study(varargin)
%MOVING_SURFACE_ADR_TP_LONG_TIME_STUDY Longer moving-surface ADR run.
%   Uses the rotating breathing sphere because it combines changing geometry
%   with rigid ambient motion while retaining an analytic manufactured
%   solution.

parser = inputParser();
parser.addParameter('CaseName', "rotating_breathing_sphere");
parser.addParameter('XiVals', [2, 4, 6], @(x) isnumeric(x) && isvector(x));
parser.addParameter('NVals', [576, 1024, 1600], @(x) isnumeric(x) && isvector(x));
parser.addParameter('FinalTime', 0.40, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('DtScale', 0.05, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('OutputPrefix', fullfile(pwd, 'moving_surface_adr_tp_long_time_rotating_breathing_sphere'));
parser.addParameter('ImageDir', fullfile(pwd, 'docs', 'figures'));
parser.addParameter('RedrawOnly', false, @(x) islogical(x) && isscalar(x));
parser.addParameter('WriteFigure', true, @(x) islogical(x) && isscalar(x));
parser.addParameter('ExactStartup', true, @(x) islogical(x) && isscalar(x));
parser.addParameter('GlobalGMRESTolerance', 1.0e-12, ...
    @(x) isnumeric(x) && isscalar(x) && x > 0);
parser.parse(varargin{:});
opts = parser.Results;

if ~exist(opts.ImageDir, 'dir')
    mkdir(opts.ImageDir);
end

if opts.RedrawOnly
    savedPath = string(opts.OutputPrefix) + ".mat";
    if ~isfile(savedPath)
        error('kp:examples:MissingSavedStudy', ...
            'Cannot redraw because saved study does not exist: %s', savedPath);
    end
    S = load(savedPath, 'study');
    study = S.study;
    study.figurePath = writeFigure(study.results, study.finalTime, opts.ImageDir);
    save(savedPath, 'study');
    return;
end

run = moving_surface_adr_tp_geometric_flow_suite(opts.CaseName, opts.XiVals, opts.NVals, opts.DtScale, ...
    'FinalTime', opts.FinalTime, ...
    'SpectrumCheck', false, ...
    'ExactStartup', opts.ExactStartup, ...
    'GlobalGMRESTolerance', opts.GlobalGMRESTolerance, ...
    'WriteOutputs', false);

study = struct();
study.caseName = string(opts.CaseName);
study.finalTime = opts.FinalTime;
study.results = run.results;
study.table = resultsTable(run.results);
study.figurePath = "";
if opts.WriteFigure
    study.figurePath = writeFigure(run.results, opts.FinalTime, opts.ImageDir);
end

save(string(opts.OutputPrefix) + ".mat", 'study');
writetable(study.table, string(opts.OutputPrefix) + ".csv");
end

function T = resultsTable(results)
rows = repmat(emptyRow(), 0, 1);
idx = 0;
for ir = 1:numel(results)
    for j = 1:numel(results(ir).N)
        idx = idx + 1;
        rows(idx, 1) = struct( ...
            'xi', results(ir).xi, ...
            'N', results(ir).N(j), ...
            'sqrtN', sqrt(results(ir).N(j)), ...
            'h', results(ir).h(j), ...
            'dt', results(ir).dt(j), ...
            'nsteps', results(ir).nsteps(j), ...
            'relerr', results(ir).relerr(j), ...
            'rate', results(ir).rate(j), ...
            'massRelError', results(ir).massRelError(j), ...
            'balanceRelResidual', results(ir).balanceRelResidual(j));
    end
end
T = struct2table(rows);
end

function imagePath = writeFigure(results, ~, imageDir)
fig = figure('Color', 'w', 'Position', [100, 100, 760, 560]);
colors = lines(numel(results));
hold on;
for ir = 1:numel(results)
    loglog(sqrt(results(ir).N), results(ir).relerr, '-o', ...
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
title('Long-time rotating breathing sphere', 'Interpreter', 'none');
legend('Location', 'best');
imagePath = fullfile(imageDir, 'moving_surface_adr_tp_long_time_rotating_breathing_sphere.png');
kp.plot.exportPaperFigure(fig, imagePath);
close(fig);
end

function row = emptyRow()
row = struct( ...
    'xi', 0, ...
    'N', 0, ...
    'sqrtN', 0, ...
    'h', NaN, ...
    'dt', NaN, ...
    'nsteps', 0, ...
    'relerr', NaN, ...
    'rate', NaN, ...
    'massRelError', NaN, ...
    'balanceRelResidual', NaN);
end
