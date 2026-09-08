function study = moving_surface_adr_tp_reduced_sbf_normal_study(varargin)
%MOVING_SURFACE_ADR_TP_REDUCED_SBF_NORMAL_STUDY Reduced SBF normal accuracy.
%   Selects deterministic control-point subsets for the global parametric SBF
%   normal model using the rate-balance rule H_M^q <= TargetSafety*h^xi,
%   then uses ControlCountScale times that count and validates the selected
%   geometry error against exact normals.

parser = inputParser();
parser.addParameter('XiVals', [2, 4, 6], @(x) isnumeric(x) && isvector(x));
parser.addParameter('Cases', defaultCases());
parser.addParameter('Degree', 7, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('TargetSafety', 0.1, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('NormalErrorOrder', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
parser.addParameter('SelectionRule', "rateBalance");
parser.addParameter('ControlCountScale', 1 / 3, @(x) isnumeric(x) && isscalar(x) && x > 0);
parser.addParameter('CandidateCounts', [], @(x) isempty(x) || isnumeric(x));
parser.addParameter('MaxControlPointCount', 2048, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('OutputPrefix', fullfile(pwd, 'moving_surface_adr_tp_reduced_sbf_normal_study'));
parser.addParameter('ImageDir', fullfile(pwd, 'docs', 'figures'));
parser.parse(varargin{:});
opts = parser.Results;

if ~exist(opts.ImageDir, 'dir')
    mkdir(opts.ImageDir);
end

cases = string(opts.Cases);
xiVals = opts.XiVals(:).';
normalErrorOrder = opts.NormalErrorOrder;
if isempty(normalErrorOrder)
    normalErrorOrder = opts.Degree + 1;
end
selectionRule = lower(string(opts.SelectionRule));
candidateRows = repmat(emptyCandidateRow(), 0, 1);
selectedRows = repmat(emptySelectedRow(), 0, 1);
candidateIndex = 0;
selectedIndex = 0;

for caseName = cases
    problem = moving_surface_adr_tp_problem(caseName);
    nvals = diagnosticNVals(caseName);
    sampleTimes = unique([0, 0.5 * problem.finalTime, problem.finalTime]);
    for xi = xiVals
        for N = nvals
            geom0 = problem.geometry(N, 0);
            h = geom0.h;
            targetRad = opts.TargetSafety * h^xi;
            targetDeg = rad2deg(targetRad);
            targetFillDistance = targetRad^(1 / normalErrorOrder);
            counts = controlCandidates(N, opts.CandidateCounts, opts.MaxControlPointCount);
            best = [];
            fprintf('\nReduced SBF normals case=%s xi=%d N=%d rule=%s target=%.3e deg\n', ...
                caseName, xi, N, selectionRule, targetDeg);

            if selectionRule == "ratebalance"
                [~, controlInfo] = kp.manifold.selectSBFControlPoints(geom0, ...
                    'FillDistanceTarget', targetFillDistance, ...
                    'MaxCount', opts.MaxControlPointCount);
                theoryCount = controlInfo.controlPointCount;
                count = min(N, max(1, ceil(opts.ControlCountScale * theoryCount)));
                [controlIds, controlInfo] = kp.manifold.selectSBFControlPoints(geom0, ...
                    'Count', count);
                [rmsDeg, maxDeg, elapsed, fillDistance] = evaluateControlIds( ...
                    problem, N, sampleTimes, opts.Degree, controlIds);
                predictedDeg = rad2deg(fillDistance^normalErrorOrder);
                meetsTarget = rmsDeg <= targetDeg;
                candidateIndex = candidateIndex + 1;
                candidateRows(candidateIndex, 1) = makeCandidateRow(problem, caseName, ...
                    xi, N, h, targetDeg, targetFillDistance, fillDistance, ...
                    predictedDeg, count, theoryCount, opts.ControlCountScale, ...
                    rmsDeg, maxDeg, elapsed, meetsTarget);
                fprintf('  M=%d theoryM=%d fraction=%.3f H=%.3e pred=%.3e deg rms=%.3e deg max=%.3e deg %s\n', ...
                    count, theoryCount, count / N, controlInfo.fillDistance, ...
                    predictedDeg, rmsDeg, maxDeg, passFail(meetsTarget));
                best = candidateRows(candidateIndex);
            elseif selectionRule == "normalerror"
                for count = counts
                    [rmsDeg, maxDeg, elapsed, fillDistance] = evaluateCount( ...
                        problem, N, sampleTimes, opts.Degree, count);
                    predictedDeg = rad2deg(fillDistance^normalErrorOrder);
                    meetsTarget = rmsDeg <= targetDeg;
                    candidateIndex = candidateIndex + 1;
                    candidateRows(candidateIndex, 1) = makeCandidateRow(problem, caseName, ...
                        xi, N, h, targetDeg, targetFillDistance, fillDistance, ...
                        predictedDeg, count, count, 1.0, rmsDeg, maxDeg, ...
                        elapsed, meetsTarget);
                    fprintf('  M=%d fraction=%.3f H=%.3e pred=%.3e deg rms=%.3e deg max=%.3e deg %s\n', ...
                        count, count / N, fillDistance, predictedDeg, rmsDeg, ...
                        maxDeg, passFail(meetsTarget));

                    if meetsTarget && isempty(best)
                        best = candidateRows(candidateIndex);
                        break;
                    end
                end
            else
                error('kp:examples:BadSBFSelectionRule', ...
                    'SelectionRule must be "rateBalance" or "normalError".');
            end

            if isempty(best)
                best = candidateRows(candidateIndex);
            end
            selectedIndex = selectedIndex + 1;
            selectedRows(selectedIndex, 1) = makeSelectedRow(best);
        end
    end
end

study = struct();
study.degree = opts.Degree;
study.targetSafety = opts.TargetSafety;
study.normalErrorOrder = normalErrorOrder;
study.selectionRule = selectionRule;
study.controlCountScale = opts.ControlCountScale;
study.xiVals = xiVals;
study.candidateTable = struct2table(candidateRows);
study.selectedTable = struct2table(selectedRows);
study.figurePaths = writeFigures(study.selectedTable, opts.ImageDir);

save(string(opts.OutputPrefix) + ".mat", 'study');
writetable(study.candidateTable, string(opts.OutputPrefix) + "_candidates.csv");
writetable(study.selectedTable, string(opts.OutputPrefix) + "_selected.csv");
end

function names = defaultCases()
names = ["mcf_sphere", "imcf_sphere", "rotating_breathing_sphere", ...
    "anisotropic_ellipsoid", "breathing_torus"];
end

function nvals = diagnosticNVals(caseName)
if string(caseName) == "breathing_torus"
    nvals = [1600, 2304, 3136];
else
    nvals = [576, 1024, 1600];
end
end

function counts = controlCandidates(N, requested, maxCount)
if isempty(requested)
    base = [32, 48, 64, 96, 128, 192, 256, 384, 512, 768, 1024, 1536, 2048];
else
    base = requested(:).';
end
maxUsable = min(N, max(1, round(maxCount)));
counts = unique(round(base(base <= maxUsable & base < N)));
if isempty(counts) || counts(end) ~= maxUsable
    counts = unique([counts, maxUsable]);
end
counts = counts(counts > 0);
end

function [rmsDeg, maxDeg, elapsed, fillDistance] = evaluateCount(problem, N, sampleTimes, degree, count)
state = [];
rmsDeg = 0;
maxDeg = 0;
elapsed = 0;
fillDistance = NaN;
for t = sampleTimes
    geom = problem.geometry(N, t);
    tic;
    [~, info, state] = kp.manifold.estimateNormalsGlobalSBF(geom, state, ...
        'Degree', degree, ...
        'ControlPointCount', count);
    elapsed = elapsed + toc;
    rmsDeg = max(rmsDeg, info.rmsAngleDegrees);
    maxDeg = max(maxDeg, info.maxAngleDegrees);
    fillDistance = info.controlFillDistance;
end
end

function [rmsDeg, maxDeg, elapsed, fillDistance] = evaluateControlIds( ...
    problem, N, sampleTimes, degree, controlIds)
state = [];
rmsDeg = 0;
maxDeg = 0;
elapsed = 0;
fillDistance = NaN;
for t = sampleTimes
    geom = problem.geometry(N, t);
    tic;
    [~, info, state] = kp.manifold.estimateNormalsGlobalSBF(geom, state, ...
        'Degree', degree, ...
        'ControlPointIds', controlIds);
    elapsed = elapsed + toc;
    rmsDeg = max(rmsDeg, info.rmsAngleDegrees);
    maxDeg = max(maxDeg, info.maxAngleDegrees);
    fillDistance = info.controlFillDistance;
end
end

function row = makeCandidateRow(problem, caseName, xi, N, h, targetDeg, ...
    targetFillDistance, fillDistance, predictedDeg, count, theoryCount, ...
    controlCountScale, rmsDeg, maxDeg, elapsed, meetsTarget)
row = emptyCandidateRow();
row.caseName = string(caseName);
row.title = string(problem.title);
row.xi = xi;
row.N = N;
row.sqrtN = sqrt(N);
row.h = h;
row.targetRMSAngleDegrees = targetDeg;
row.targetFillDistance = targetFillDistance;
row.controlFillDistance = fillDistance;
row.predictedRMSAngleDegrees = predictedDeg;
row.controlPointCount = count;
row.theoryControlPointCount = theoryCount;
row.controlCountScale = controlCountScale;
row.controlPointFraction = count / N;
row.rmsAngleDegrees = rmsDeg;
row.maxAngleDegrees = maxDeg;
row.elapsedSeconds = elapsed;
row.meetsTarget = meetsTarget;
end

function row = makeSelectedRow(candidate)
row = emptySelectedRow();
row.caseName = candidate.caseName;
row.title = candidate.title;
row.xi = candidate.xi;
row.N = candidate.N;
row.sqrtN = candidate.sqrtN;
row.h = candidate.h;
row.targetRMSAngleDegrees = candidate.targetRMSAngleDegrees;
row.targetFillDistance = candidate.targetFillDistance;
row.controlFillDistance = candidate.controlFillDistance;
row.predictedRMSAngleDegrees = candidate.predictedRMSAngleDegrees;
row.selectedControlPointCount = candidate.controlPointCount;
row.theoryControlPointCount = candidate.theoryControlPointCount;
row.controlCountScale = candidate.controlCountScale;
row.selectedControlPointFraction = candidate.controlPointFraction;
row.rmsAngleDegrees = candidate.rmsAngleDegrees;
row.maxAngleDegrees = candidate.maxAngleDegrees;
row.elapsedSeconds = candidate.elapsedSeconds;
row.meetsTarget = candidate.meetsTarget;
end

function figures = writeFigures(T, imageDir)
cases = unique(T.caseName, 'stable').';
figures = strings(2, numel(cases));
for ic = 1:numel(cases)
    Tc = T(T.caseName == cases(ic), :);
    figures(1, ic) = writeFractionFigure(Tc, imageDir, cases(ic));
    figures(2, ic) = writeAngleFigure(Tc, imageDir, cases(ic));
end
end

function imagePath = writeFractionFigure(T, imageDir, caseName)
fig = figure('Color', 'w', 'Position', [100, 100, 760, 560]);
hold on;
plotSelected(T, "selectedControlPointFraction");
hold off;
set(gca, 'XScale', 'log', 'YScale', 'log');
grid on;
xlabel('$\sqrt{N}$', 'Interpreter', 'latex');
ylabel('$M/N$', 'Interpreter', 'latex');
title(sprintf('%s reduced SBF control fraction', caseLabel(T)), 'Interpreter', 'none');
legend('Location', 'best');
imagePath = fullfile(imageDir, ...
    sprintf('moving_surface_adr_tp_reduced_sbf_fraction_%s.png', caseName));
exportgraphics(fig, imagePath, 'Resolution', 220);
close(fig);
end

function imagePath = writeAngleFigure(T, imageDir, caseName)
fig = figure('Color', 'w', 'Position', [100, 100, 760, 560]);
xiVals = unique(T.xi(:).');
colors = lines(numel(xiVals));
for i = 1:numel(xiVals)
    xi = xiVals(i);
    idx = T.xi == xi;
    [~, order] = sort(T.sqrtN(idx));
    Tc = T(idx, :);
    Tc = Tc(order, :);
    loglog(Tc.sqrtN, Tc.rmsAngleDegrees, '-o', ...
        'LineWidth', 1.5, 'MarkerSize', 6, 'Color', colors(i, :), ...
        'DisplayName', sprintf('\\xi = %d, observed', xi));
    hold on;
    loglog(Tc.sqrtN, Tc.targetRMSAngleDegrees, '--', ...
        'LineWidth', 1.2, 'Color', colors(i, :), ...
        'DisplayName', sprintf('\\xi = %d, target', xi));
end
hold off;
set(gca, 'XScale', 'log', 'YScale', 'log');
grid on;
xlabel('$\sqrt{N}$', 'Interpreter', 'latex');
ylabel('RMS normal angle (degrees)');
title(sprintf('%s reduced SBF normal target', caseLabel(T)), 'Interpreter', 'none');
legend('Location', 'best');
imagePath = fullfile(imageDir, ...
    sprintf('moving_surface_adr_tp_reduced_sbf_angle_%s.png', caseName));
exportgraphics(fig, imagePath, 'Resolution', 220);
close(fig);
end

function plotSelected(T, variableName)
xiVals = unique(T.xi(:).');
for xi = xiVals
    idx = T.xi == xi;
    [~, order] = sort(T.sqrtN(idx));
    Tc = T(idx, :);
    Tc = Tc(order, :);
    loglog(Tc.sqrtN, Tc.(variableName), '-o', ...
        'LineWidth', 1.5, 'MarkerSize', 6, ...
        'DisplayName', sprintf('\\xi = %d', xi));
end
end

function text = caseLabel(T)
text = char(T.title(1));
end

function text = passFail(tf)
if tf
    text = 'ok';
else
    text = 'miss';
end
end

function row = emptyCandidateRow()
row = struct( ...
    'caseName', "", ...
    'title', "", ...
    'xi', 0, ...
    'N', 0, ...
    'sqrtN', NaN, ...
    'h', NaN, ...
    'targetRMSAngleDegrees', NaN, ...
    'targetFillDistance', NaN, ...
    'controlFillDistance', NaN, ...
    'predictedRMSAngleDegrees', NaN, ...
    'controlPointCount', 0, ...
    'theoryControlPointCount', 0, ...
    'controlCountScale', NaN, ...
    'controlPointFraction', NaN, ...
    'rmsAngleDegrees', NaN, ...
    'maxAngleDegrees', NaN, ...
    'elapsedSeconds', NaN, ...
    'meetsTarget', false);
end

function row = emptySelectedRow()
row = struct( ...
    'caseName', "", ...
    'title', "", ...
    'xi', 0, ...
    'N', 0, ...
    'sqrtN', NaN, ...
    'h', NaN, ...
    'targetRMSAngleDegrees', NaN, ...
    'targetFillDistance', NaN, ...
    'controlFillDistance', NaN, ...
    'predictedRMSAngleDegrees', NaN, ...
    'selectedControlPointCount', 0, ...
    'theoryControlPointCount', 0, ...
    'controlCountScale', NaN, ...
    'selectedControlPointFraction', NaN, ...
    'rmsAngleDegrees', NaN, ...
    'maxAngleDegrees', NaN, ...
    'elapsedSeconds', NaN, ...
    'meetsTarget', false);
end
