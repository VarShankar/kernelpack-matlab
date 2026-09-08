function bench = moving_surface_adr_tp_global_solver_benchmark(caseNames, xiVals, Nvals, dtScale, varargin)
%MOVING_SURFACE_ADR_TP_GLOBAL_SOLVER_BENCHMARK Isolate global solve defect sweeps.
%   Compares frozen-ILU GMRES with no global defect sweeps against the same
%   solver with a small number of frozen-ILU defect-correction sweeps before
%   GMRES.  Differentiation matrices are direct-reassembled by default so the
%   timing isolates the global BDF linear solve.

if nargin < 1 || isempty(caseNames)
    caseNames = ["mcf_sphere", "imcf_sphere", "rotating_breathing_sphere", ...
        "anisotropic_ellipsoid", "breathing_torus"];
end
if nargin < 2 || isempty(xiVals)
    xiVals = 4;
end
if nargin < 3 || isempty(Nvals)
    Nvals = paperNVals();
end
if nargin < 4 || isempty(dtScale)
    dtScale = 0.05;
end

parser = inputParser();
parser.addParameter('SaveResults', true);
parser.addParameter('OutputPrefix', fullfile(pwd, 'moving_surface_adr_tp_global_solver_benchmark_results_xi4'));
parser.addParameter('WarmupParallelPool', true);
parser.addParameter('ContinueOnError', true);
parser.addParameter('SaveAfterEach', true);
parser.addParameter('DiffMatUpdateMethod', "direct");
parser.addParameter('DefectTolerance', 1.0e-6);
parser.addParameter('MaxDefectIterations', 4);
parser.addParameter('NeighborUpdateMode', "periodic");
parser.addParameter('NeighborSearchInterval', 5);
parser.addParameter('SpectrumCheck', false);
parser.addParameter('HyperviscosityUpdateMode', "adaptiveGeometry");
parser.addParameter('HyperviscosityDriftTolerance', 0.05);
parser.addParameter('MaxHyperviscositySkippedSteps', 5);
parser.addParameter('GlobalSolveMethod', "gmresIluDefect");
parser.addParameter('GlobalGMRESTolerance', 1.0e-6);
parser.addParameter('ExactStartup', false, @(x) islogical(x) && isscalar(x));
parser.addParameter('GlobalGMRESRestart', 40);
parser.addParameter('GlobalGMRESMaxIterations', 30);
parser.addParameter('GlobalDefectSweepsOff', 0);
parser.addParameter('GlobalDefectSweepsOn', 4);
parser.addParameter('GlobalILUDropTolerance', 1.0e-4);
parser.addParameter('GlobalILURefreshInterval', Inf);
parser.addParameter('GlobalILURefreshIterationThreshold', 30);
parser.addParameter('MassCorrectionMode', "balance", @(x) isstring(x) || ischar(x));
parser.parse(varargin{:});

caseNames = string(caseNames);
xiVals = xiVals(:).';
Nvals = Nvals(:).';

if parser.Results.WarmupParallelPool
    warmupParallelPool();
end

rows = repmat(emptyRow(), 0, 1);
rowIndex = 0;

for icase = 1:numel(caseNames)
    caseName = caseNames(icase);
    for ixi = 1:numel(xiVals)
        xi = xiVals(ixi);
        for iN = 1:numel(Nvals)
            N = Nvals(iN);
            fprintf('\nGlobal solver benchmark case=%s xi=%d N=%d\n', caseName, xi, N);

            noDefect = runTimed(caseName, xi, N, dtScale, ...
                parser.Results.GlobalDefectSweepsOff, parser.Results);
            defect = runTimed(caseName, xi, N, dtScale, ...
                parser.Results.GlobalDefectSweepsOn, parser.Results);

            rowIndex = rowIndex + 1;
            rows(rowIndex, 1) = makeRow(caseName, xi, N, noDefect, defect);
            if rows(rowIndex).status == "ok"
                fprintf('  global solve speedup: %.3f\n', rows(rowIndex).globalSolveSpeedup);
                fprintf('  GMRES iterations off/on: %d / %d\n', ...
                    rows(rowIndex).noDefectGmresIterations, rows(rowIndex).defectGmresIterations);
            else
                fprintf('  status: %s\n', rows(rowIndex).status);
            end

            if parser.Results.SaveResults && parser.Results.SaveAfterEach
                bench = buildBench(caseNames, xiVals, Nvals, dtScale, rows, parser.Results);
                saveBenchmark(bench, parser.Results.OutputPrefix);
            end
        end
    end
end

bench = buildBench(caseNames, xiVals, Nvals, dtScale, rows, parser.Results);

if parser.Results.SaveResults
    saveBenchmark(bench, parser.Results.OutputPrefix);
end
dispBenchmarkSummary(bench);
end

function Nvals = paperNVals()
Nvals = [576, 1024, 1600, 2500, 3600, 4900];
end

function out = runTimed(caseName, xi, N, dtScale, globalDefectSweeps, opts)
tic;
try
    result = moving_surface_adr_tp_geometric_flow_suite(caseName, xi, N, dtScale, ...
        'DiffMatUpdateMethod', opts.DiffMatUpdateMethod, ...
        'DefectTolerance', opts.DefectTolerance, ...
        'MaxDefectIterations', opts.MaxDefectIterations, ...
        'NeighborUpdateMode', opts.NeighborUpdateMode, ...
        'NeighborSearchInterval', opts.NeighborSearchInterval, ...
        'SpectrumCheck', opts.SpectrumCheck, ...
        'HyperviscosityUpdateMode', opts.HyperviscosityUpdateMode, ...
        'HyperviscosityDriftTolerance', opts.HyperviscosityDriftTolerance, ...
        'MaxHyperviscositySkippedSteps', opts.MaxHyperviscositySkippedSteps, ...
        'GlobalSolveMethod', opts.GlobalSolveMethod, ...
        'GlobalGMRESTolerance', opts.GlobalGMRESTolerance, ...
        'ExactStartup', opts.ExactStartup, ...
        'GlobalGMRESRestart', opts.GlobalGMRESRestart, ...
        'GlobalGMRESMaxIterations', opts.GlobalGMRESMaxIterations, ...
        'GlobalDefectSweeps', globalDefectSweeps, ...
        'GlobalILUDropTolerance', opts.GlobalILUDropTolerance, ...
        'GlobalILURefreshInterval', opts.GlobalILURefreshInterval, ...
        'GlobalILURefreshIterationThreshold', opts.GlobalILURefreshIterationThreshold, ...
        'MassCorrectionMode', opts.MassCorrectionMode, ...
        'WriteOutputs', false);
    elapsed = toc;

    R = result.results;
    out = struct();
    out.ok = true;
    out.elapsed = elapsed;
    out.err = R.relerr;
    out.h = R.h;
    out.dt = R.dt;
    out.nsteps = R.nsteps;
    out.massRelError = R.massRelError;
    out.balanceRelResidual = R.balanceRelResidual;
    out.solveStats = R.solveStats{1};
    out.errorIdentifier = "";
    out.errorMessage = "";
catch ME
    elapsed = toc;
    if ~opts.ContinueOnError
        rethrow(ME);
    end
    out = failedRun(elapsed, ME);
    fprintf('  sweeps=%d failed: %s\n', globalDefectSweeps, ME.message);
end
end

function out = failedRun(elapsed, ME)
out = struct();
out.ok = false;
out.elapsed = elapsed;
out.err = NaN;
out.h = NaN;
out.dt = NaN;
out.nsteps = 0;
out.massRelError = NaN;
out.balanceRelResidual = NaN;
out.solveStats = [];
out.errorIdentifier = string(ME.identifier);
out.errorMessage = string(ME.message);
end

function row = makeRow(caseName, xi, N, noDefect, defect)
row = emptyRow();
row.caseName = string(caseName);
row.xi = xi;
row.N = N;
if noDefect.ok
    row.h = noDefect.h;
    row.dt = noDefect.dt;
    row.nsteps = noDefect.nsteps;
elseif defect.ok
    row.h = defect.h;
    row.dt = defect.dt;
    row.nsteps = defect.nsteps;
end

row.noDefectTime = noDefect.elapsed;
row.defectTime = defect.elapsed;
row.noDefectError = noDefect.err;
row.defectError = defect.err;
row.errorDifference = abs(defect.err - noDefect.err);
row.noDefectMassRelError = noDefect.massRelError;
row.defectMassRelError = defect.massRelError;
row.noDefectBalanceRelResidual = noDefect.balanceRelResidual;
row.defectBalanceRelResidual = defect.balanceRelResidual;
row.noDefectErrorIdentifier = noDefect.errorIdentifier;
row.noDefectErrorMessage = noDefect.errorMessage;
row.defectErrorIdentifier = defect.errorIdentifier;
row.defectErrorMessage = defect.errorMessage;

if noDefect.ok && defect.ok
    row.overallSpeedup = noDefect.elapsed / max(defect.elapsed, eps);
    row.status = "ok";
elseif noDefect.ok
    row.status = "defect_failed";
elseif defect.ok
    row.status = "no_defect_failed";
else
    row.status = "both_failed";
end

if ~isempty(noDefect.solveStats)
    row = copySolveStats(row, "noDefect", noDefect.solveStats);
end
if ~isempty(defect.solveStats)
    row = copySolveStats(row, "defect", defect.solveStats);
end
if noDefect.ok && defect.ok
    row.globalSolveSpeedup = row.noDefectSolveTime / max(row.defectSolveTime, eps);
    row.gmresIterationReduction = row.noDefectGmresIterations - row.defectGmresIterations;
end
end

function row = copySolveStats(row, prefix, stats)
row.(prefix + "SolveTime") = stats.solveTime;
row.(prefix + "IluRefreshTime") = stats.iluRefreshTime;
row.(prefix + "DefectSweepTime") = stats.defectSweepTime;
row.(prefix + "GmresTime") = stats.gmresTime;
row.(prefix + "IluRefreshes") = stats.iluRefreshes;
row.(prefix + "IluReuses") = stats.iluReuses;
row.(prefix + "GmresSolves") = stats.gmresSolves;
row.(prefix + "GmresInvocations") = stats.gmresInvocations;
row.(prefix + "GmresIterations") = stats.gmresIterations;
row.(prefix + "MaxGmresIterations") = stats.maxGMRESIterations;
row.(prefix + "GmresFailures") = stats.gmresFailures;
row.(prefix + "GlobalDefectSweeps") = stats.globalDefectSweeps;
row.(prefix + "InitialGuessAcceptedSolves") = stats.initialGuessAcceptedSolves;
row.(prefix + "DefectAcceptedSolves") = stats.defectAcceptedSolves;
row.(prefix + "MaxInitialResidual") = stats.maxInitialRelativeResidual;
row.(prefix + "MaxPostDefectResidual") = stats.maxPostDefectRelativeResidual;
row.(prefix + "MaxSolveResidual") = stats.maxRelativeResidual;
end

function bench = buildBench(caseNames, xiVals, Nvals, dtScale, rows, opts)
bench = struct();
bench.caseNames = caseNames;
bench.xiVals = xiVals;
bench.Nvals = Nvals;
bench.dtScale = dtScale;
bench.diffMatUpdateMethod = string(opts.DiffMatUpdateMethod);
bench.hyperviscosityUpdateMode = string(opts.HyperviscosityUpdateMode);
bench.globalSolveMethod = string(opts.GlobalSolveMethod);
bench.globalGMRESTolerance = opts.GlobalGMRESTolerance;
bench.globalDefectSweepsOff = opts.GlobalDefectSweepsOff;
bench.globalDefectSweepsOn = opts.GlobalDefectSweepsOn;
bench.globalILUDropTolerance = opts.GlobalILUDropTolerance;
bench.globalILURefreshInterval = opts.GlobalILURefreshInterval;
bench.massCorrectionMode = string(opts.MassCorrectionMode);
bench.rows = rows;
bench.table = benchmarkTable(rows);
bench.summary = summarizeRows(rows);
end

function saveBenchmark(bench, outputPrefix)
outputPrefix = string(outputPrefix);
save(outputPrefix + ".mat", 'bench');
writetable(bench.table, outputPrefix + ".csv");
end

function T = benchmarkTable(rows)
if isempty(rows)
    T = table();
    return;
end
T = struct2table(rows);
T.caseName = string(T.caseName);
end

function summary = summarizeRows(rows)
summary = struct();
summary.geomeanOverallSpeedup = NaN;
summary.minOverallSpeedup = NaN;
summary.maxOverallSpeedup = NaN;
summary.geomeanGlobalSolveSpeedup = NaN;
summary.minGlobalSolveSpeedup = NaN;
summary.maxGlobalSolveSpeedup = NaN;
summary.gmresIterationReduction = 0;
if isempty(rows)
    return;
end

overallSpeedups = [rows.overallSpeedup];
overallSpeedups = overallSpeedups(isfinite(overallSpeedups) & overallSpeedups > 0);
if ~isempty(overallSpeedups)
    summary.geomeanOverallSpeedup = exp(mean(log(overallSpeedups)));
    summary.minOverallSpeedup = min(overallSpeedups);
    summary.maxOverallSpeedup = max(overallSpeedups);
end

globalSpeedups = [rows.globalSolveSpeedup];
globalSpeedups = globalSpeedups(isfinite(globalSpeedups) & globalSpeedups > 0);
if ~isempty(globalSpeedups)
    summary.geomeanGlobalSolveSpeedup = exp(mean(log(globalSpeedups)));
    summary.minGlobalSolveSpeedup = min(globalSpeedups);
    summary.maxGlobalSolveSpeedup = max(globalSpeedups);
end

summary.gmresIterationReduction = sum([rows.gmresIterationReduction]);
end

function dispBenchmarkSummary(bench)
fprintf('\nMoving-surface global-solver benchmark summary\n');
fprintf('  rows: %d\n', numel(bench.rows));
fprintf('  global solve speedup geomean/min/max: %.3f / %.3f / %.3f\n', ...
    bench.summary.geomeanGlobalSolveSpeedup, ...
    bench.summary.minGlobalSolveSpeedup, ...
    bench.summary.maxGlobalSolveSpeedup);
fprintf('  overall speedup geomean/min/max: %.3f / %.3f / %.3f\n', ...
    bench.summary.geomeanOverallSpeedup, ...
    bench.summary.minOverallSpeedup, ...
    bench.summary.maxOverallSpeedup);
fprintf('  total GMRES iteration reduction: %d\n', ...
    bench.summary.gmresIterationReduction);
end

function warmupParallelPool()
try
    pool = gcp('nocreate');
    if isempty(pool)
        parpool;
    end
catch ME
    warning('kp:examples:ParallelWarmupSkipped', ...
        'Could not warm up the parallel pool before timing: %s', ME.message);
end
end

function row = emptyRow()
row = struct( ...
    'caseName', "", ...
    'status', "", ...
    'xi', 0, ...
    'N', 0, ...
    'h', NaN, ...
    'dt', NaN, ...
    'nsteps', 0, ...
    'noDefectTime', NaN, ...
    'defectTime', NaN, ...
    'overallSpeedup', NaN, ...
    'noDefectError', NaN, ...
    'defectError', NaN, ...
    'errorDifference', NaN, ...
    'noDefectMassRelError', NaN, ...
    'defectMassRelError', NaN, ...
    'noDefectBalanceRelResidual', NaN, ...
    'defectBalanceRelResidual', NaN, ...
    'noDefectErrorIdentifier', "", ...
    'noDefectErrorMessage', "", ...
    'defectErrorIdentifier', "", ...
    'defectErrorMessage', "", ...
    'noDefectSolveTime', NaN, ...
    'defectSolveTime', NaN, ...
    'globalSolveSpeedup', NaN, ...
    'noDefectIluRefreshTime', NaN, ...
    'defectIluRefreshTime', NaN, ...
    'noDefectDefectSweepTime', NaN, ...
    'defectDefectSweepTime', NaN, ...
    'noDefectGmresTime', NaN, ...
    'defectGmresTime', NaN, ...
    'noDefectIluRefreshes', 0, ...
    'defectIluRefreshes', 0, ...
    'noDefectIluReuses', 0, ...
    'defectIluReuses', 0, ...
    'noDefectGmresSolves', 0, ...
    'defectGmresSolves', 0, ...
    'noDefectGmresInvocations', 0, ...
    'defectGmresInvocations', 0, ...
    'noDefectGmresIterations', 0, ...
    'defectGmresIterations', 0, ...
    'gmresIterationReduction', 0, ...
    'noDefectMaxGmresIterations', 0, ...
    'defectMaxGmresIterations', 0, ...
    'noDefectGmresFailures', 0, ...
    'defectGmresFailures', 0, ...
    'noDefectGlobalDefectSweeps', 0, ...
    'defectGlobalDefectSweeps', 0, ...
    'noDefectInitialGuessAcceptedSolves', 0, ...
    'defectInitialGuessAcceptedSolves', 0, ...
    'noDefectDefectAcceptedSolves', 0, ...
    'defectDefectAcceptedSolves', 0, ...
    'noDefectMaxInitialResidual', NaN, ...
    'defectMaxInitialResidual', NaN, ...
    'noDefectMaxPostDefectResidual', NaN, ...
    'defectMaxPostDefectResidual', NaN, ...
    'noDefectMaxSolveResidual', NaN, ...
    'defectMaxSolveResidual', NaN);
end
