function study = moving_surface_adr_tp_exact_pde_study(varargin)
%MOVING_SURFACE_ADR_TP_EXACT_PDE_STUDY Current exact-normal PDE reference runs.

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
parser.addParameter('OutputPrefix', fullfile(pwd, ...
    'moving_surface_adr_tp_exact_pde_study'));
parser.parse(varargin{:});
opts = parser.Results;

cases = string(opts.Cases);
xiVals = opts.XiVals(:).';
nvals = opts.NVals(:).';
rows = repmat(emptyRow(), 0, 1);
rowIndex = 0;

for caseName = cases
    problem = moving_surface_adr_tp_problem(caseName);
    for xi = xiVals
        for N = nvals
            fprintf('\nExact-normal PDE case=%s xi=%d N=%d\n', caseName, xi, N);
            tic;
            result = kp.manifold.runLagrangianMovingADRConvergence( ...
                problem, xi, N, opts.DtScale, ...
                'FinalTime', problem.finalTime, ...
                'SpectrumCheck', false, ...
                'NormalMode', "exact", ...
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
                result, elapsed);
            fprintf('  exact relerr=%.6e wall=%.2fs\n', result.relerr, elapsed);
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

function row = makeRow(caseName, titleText, xi, N, result, elapsed)
row = emptyRow();
row.caseName = string(caseName);
row.title = string(titleText);
row.xi = xi;
row.N = N;
row.relerr = result.relerr;
row.massRelError = result.massRelError;
row.balanceRelResidual = result.balanceRelResidual;
row.elapsedSeconds = elapsed;
end

function row = emptyRow()
row = struct( ...
    'caseName', "", ...
    'title', "", ...
    'xi', 0, ...
    'N', 0, ...
    'relerr', NaN, ...
    'massRelError', NaN, ...
    'balanceRelResidual', NaN, ...
    'elapsedSeconds', NaN);
end
