function study = moving_surface_adr_tp_calibrated_sbf_pde_study(varargin)
%MOVING_SURFACE_ADR_TP_CALIBRATED_SBF_PDE_STUDY PDE validation for calibrated SBF normals.
%   Uses a precomputed calibrated reduced-SBF normal table, runs the moving
%   surface ADR solve with those normals, and compares against archived
%   exact-normal convergence results.

parser = inputParser();
parser.addParameter('Cases', defaultCases());
parser.addParameter('XiVals', 6, @(x) isnumeric(x) && isvector(x));
parser.addParameter('NVals', 1600, @(x) isnumeric(x) && isvector(x));
parser.addParameter('DtScale', 0.05, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('UseRearrangement', true, @(x) islogical(x) && isscalar(x));
parser.addParameter('QualityThreshold', 1.75, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('PredictiveLookaheadSteps', 3, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('MinStepsBetweenRearrangements', 8, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('RearrangementTransferMode', "localTp", @(x) isstring(x) || ischar(x));
parser.addParameter('MassCorrectionMode', "balance", @(x) isstring(x) || ischar(x));
parser.addParameter('SelectionTable', fullfile(pwd, ...
    'moving_surface_adr_tp_reduced_sbf_normal_study_calibrated_selected.csv'));
parser.addParameter('OutputPrefix', fullfile(pwd, ...
    'moving_surface_adr_tp_calibrated_sbf_pde_study'));
parser.parse(varargin{:});
opts = parser.Results;

selectionTable = readtable(opts.SelectionTable, 'TextType', 'string');
cases = string(opts.Cases);
xiVals = opts.XiVals(:).';
nvals = opts.NVals(:).';
rows = repmat(emptyRow(), 0, 1);
rowIndex = 0;

for caseName = cases
    problem = moving_surface_adr_tp_problem(caseName);
    for xi = xiVals
        for N = nvals
            selection = selectedRow(selectionTable, caseName, xi, N);
            exact = archivedExactResult(problem, xi, N);
            fprintf('\nCalibrated SBF PDE case=%s xi=%d N=%d M=%d theoryM=%d\n', ...
                caseName, xi, N, selection.selectedControlPointCount, ...
                selection.theoryControlPointCount);
            tic;
            result = kp.manifold.runLagrangianMovingADRConvergence( ...
                problem, xi, N, opts.DtScale, ...
                'FinalTime', problem.finalTime, ...
                'SpectrumCheck', false, ...
                'NormalMode', "globalSbf", ...
                'GlobalSBFControlPointCount', selection.selectedControlPointCount, ...
                'UseRearrangement', opts.UseRearrangement, ...
                'QualityThreshold', opts.QualityThreshold, ...
                'PredictiveLookaheadSteps', opts.PredictiveLookaheadSteps, ...
                'MinStepsBetweenRearrangements', opts.MinStepsBetweenRearrangements, ...
                'RearrangementTransferMode', opts.RearrangementTransferMode, ...
                'MassCorrectionMode', opts.MassCorrectionMode, ...
                'RecordMassDiagnostics', true);
            elapsed = toc;
            rowIndex = rowIndex + 1;
            rows(rowIndex, 1) = makeRow(caseName, problem.title, xi, N, ...
                selection, exact, result, elapsed);
            fprintf('  reduced relerr=%.6e exact relerr=%.6e ratio=%.6f normalRMS=%.3e wall=%.2fs\n', ...
                result.relerr, exact.relerr, result.relerr / exact.relerr, ...
                result.normalRMSAngle, elapsed);
        end
    end
end

study = struct();
study.rows = rows;
study.table = struct2table(rows);
save(string(opts.OutputPrefix) + ".mat", 'study');
writetable(study.table, string(opts.OutputPrefix) + ".csv");
end

function names = defaultCases()
names = ["mcf_sphere", "imcf_sphere", "rotating_breathing_sphere", ...
    "anisotropic_ellipsoid", "breathing_torus"];
end

function row = selectedRow(T, caseName, xi, N)
idx = T.caseName == caseName & T.xi == xi & T.N == N;
if ~any(idx)
    error('kp:examples:MissingCalibratedSelection', ...
        'No calibrated reduced-SBF row for case=%s xi=%d N=%d.', caseName, xi, N);
end
row = T(find(idx, 1), :);
end

function exact = archivedExactResult(problem, xi, N)
S = load(problem.resultsName, 'results');
hits = find([S.results.xi] == xi, 1);
if isempty(hits)
    error('kp:examples:MissingExactXi', ...
        'No archived exact-normal result for xi=%d in %s.', xi, problem.resultsName);
end
result = S.results(hits);
level = find(result.N == N, 1);
if isempty(level)
    error('kp:examples:MissingExactN', ...
        'No archived exact-normal result for N=%d in %s.', N, problem.resultsName);
end
exact = struct('relerr', result.relerr(level), ...
    'massRelError', result.massRelError(level), ...
    'balanceRelResidual', result.balanceRelResidual(level));
end

function row = makeRow(caseName, titleText, xi, N, selection, exact, result, elapsed)
row = emptyRow();
row.caseName = string(caseName);
row.title = string(titleText);
row.xi = xi;
row.N = N;
row.selectedControlPointCount = selection.selectedControlPointCount;
row.theoryControlPointCount = selection.theoryControlPointCount;
row.controlPointFraction = selection.selectedControlPointFraction;
row.normalRMSAngle = result.normalRMSAngle;
row.normalMaxAngle = result.normalMaxAngle;
row.reducedRelErr = result.relerr;
row.exactRelErr = exact.relerr;
row.relErrRatio = result.relerr / exact.relerr;
row.relErrDelta = result.relerr - exact.relerr;
row.reducedMassRelError = result.massRelError;
row.exactMassRelError = exact.massRelError;
row.reducedBalanceRelResidual = result.balanceRelResidual;
row.exactBalanceRelResidual = exact.balanceRelResidual;
row.elapsedSeconds = elapsed;
end

function row = emptyRow()
row = struct( ...
    'caseName', "", ...
    'title', "", ...
    'xi', 0, ...
    'N', 0, ...
    'selectedControlPointCount', 0, ...
    'theoryControlPointCount', 0, ...
    'controlPointFraction', NaN, ...
    'normalRMSAngle', NaN, ...
    'normalMaxAngle', NaN, ...
    'reducedRelErr', NaN, ...
    'exactRelErr', NaN, ...
    'relErrRatio', NaN, ...
    'relErrDelta', NaN, ...
    'reducedMassRelError', NaN, ...
    'exactMassRelError', NaN, ...
    'reducedBalanceRelResidual', NaN, ...
    'exactBalanceRelResidual', NaN, ...
    'elapsedSeconds', NaN);
end
