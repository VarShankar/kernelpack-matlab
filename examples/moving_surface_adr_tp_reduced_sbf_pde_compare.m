function comparison = moving_surface_adr_tp_reduced_sbf_pde_compare(varargin)
%MOVING_SURFACE_ADR_TP_REDUCED_SBF_PDE_COMPARE Compare reduced SBF choices in PDE solves.
%   Runs exact normals, empirically selected reduced-SBF normals, and
%   rate-balanced reduced-SBF normals for selected moving-surface ADR cases.

parser = inputParser();
parser.addParameter('Cases', ["rotating_breathing_sphere", "breathing_torus"]);
parser.addParameter('Xi', 6, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('NByCase', struct(), @(x) isstruct(x));
parser.addParameter('DtScale', 0.05, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('UseRearrangement', true, @(x) islogical(x) && isscalar(x));
parser.addParameter('QualityThreshold', 1.75, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('PredictiveLookaheadSteps', 3, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('MinStepsBetweenRearrangements', 8, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('RearrangementTransferMode', "localTp", @(x) isstring(x) || ischar(x));
parser.addParameter('MassCorrectionMode', "balance", @(x) isstring(x) || ischar(x));
parser.addParameter('EmpiricalTable', fullfile(pwd, 'moving_surface_adr_tp_reduced_sbf_normal_study_selected.csv'));
parser.addParameter('RateTable', fullfile(pwd, 'moving_surface_adr_tp_reduced_sbf_normal_study_rate_selected.csv'));
parser.addParameter('OutputPrefix', fullfile(pwd, 'moving_surface_adr_tp_reduced_sbf_pde_compare_xi6'));
parser.parse(varargin{:});
opts = parser.Results;

empiricalTable = readtable(opts.EmpiricalTable, 'TextType', 'string');
rateTable = readtable(opts.RateTable, 'TextType', 'string');
cases = string(opts.Cases);
rows = repmat(emptyRow(), 0, 1);
rowIndex = 0;

for caseName = cases
    problem = moving_surface_adr_tp_problem(caseName);
    N = selectedN(caseName, opts.NByCase);
    selections = [
        struct('label', "exact", 'normalMode', "exact", 'M', Inf)
        selectedEntry("empirical", empiricalTable, caseName, opts.Xi, N)
        selectedEntry("rate", rateTable, caseName, opts.Xi, N)
        ];

    for k = 1:numel(selections)
        sel = selections(k);
        fprintf('\nPDE normal comparison case=%s xi=%d N=%d selector=%s M=%g\n', ...
            caseName, opts.Xi, N, sel.label, sel.M);
        tic;
        result = kp.manifold.runLagrangianMovingADRConvergence( ...
            problem, opts.Xi, N, opts.DtScale, ...
            'FinalTime', problem.finalTime, ...
            'SpectrumCheck', false, ...
            'NormalMode', sel.normalMode, ...
            'GlobalSBFControlPointCount', sel.M, ...
            'UseRearrangement', opts.UseRearrangement, ...
            'QualityThreshold', opts.QualityThreshold, ...
            'PredictiveLookaheadSteps', opts.PredictiveLookaheadSteps, ...
            'MinStepsBetweenRearrangements', opts.MinStepsBetweenRearrangements, ...
            'RearrangementTransferMode', opts.RearrangementTransferMode, ...
            'MassCorrectionMode', opts.MassCorrectionMode, ...
            'RecordMassDiagnostics', true);
        elapsed = toc;
        rowIndex = rowIndex + 1;
        rows(rowIndex, 1) = makeRow(caseName, problem.title, opts.Xi, N, sel, ...
            result, elapsed);
        fprintf('  relerr=%.6e normalRMS=%.3e wall=%.2fs\n', ...
            result.relerr, result.normalRMSAngle, elapsed);
    end
end

comparison = struct();
comparison.rows = rows;
comparison.table = struct2table(rows);
save(string(opts.OutputPrefix) + ".mat", 'comparison');
writetable(comparison.table, string(opts.OutputPrefix) + ".csv");
end

function N = selectedN(caseName, nByCase)
field = matlab.lang.makeValidName(char(caseName));
if isfield(nByCase, field)
    N = nByCase.(field);
elseif string(caseName) == "breathing_torus"
    N = 3136;
else
    N = 1600;
end
end

function sel = selectedEntry(label, tableData, caseName, xi, N)
idx = tableData.caseName == caseName & tableData.xi == xi & tableData.N == N;
if ~any(idx)
    error('kp:examples:MissingSBFSelection', ...
        'No %s reduced-SBF selection for case=%s xi=%d N=%d.', label, caseName, xi, N);
end
row = tableData(find(idx, 1), :);
sel = struct('label', string(label), ...
    'normalMode', "globalSbf", ...
    'M', row.selectedControlPointCount);
end

function row = makeRow(caseName, titleText, xi, N, selection, result, elapsed)
row = emptyRow();
row.caseName = string(caseName);
row.title = string(titleText);
row.xi = xi;
row.N = N;
row.selector = selection.label;
row.controlPointCount = selection.M;
if isfinite(selection.M)
    row.controlPointFraction = selection.M / N;
else
    row.controlPointFraction = NaN;
end
row.relerr = result.relerr;
row.massRelError = result.massRelError;
row.balanceRelResidual = result.balanceRelResidual;
row.normalRMSAngle = result.normalRMSAngle;
row.normalMaxAngle = result.normalMaxAngle;
row.elapsedSeconds = elapsed;
end

function row = emptyRow()
row = struct( ...
    'caseName', "", ...
    'title', "", ...
    'xi', 0, ...
    'N', 0, ...
    'selector', "", ...
    'controlPointCount', NaN, ...
    'controlPointFraction', NaN, ...
    'relerr', NaN, ...
    'massRelError', NaN, ...
    'balanceRelResidual', NaN, ...
    'normalRMSAngle', NaN, ...
    'normalMaxAngle', NaN, ...
    'elapsedSeconds', NaN);
end
