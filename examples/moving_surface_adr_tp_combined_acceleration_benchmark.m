function bench = moving_surface_adr_tp_combined_acceleration_benchmark(caseNames, xiVals, Nvals, dtScale, varargin)
%MOVING_SURFACE_ADR_TP_COMBINED_ACCELERATION_BENCHMARK End-to-end timing.
%   Baseline: direct differentiation matrices, every-step hyperviscosity, and
%   frozen-ILU GMRES without global defect sweeps.  Accelerated: defect-corrected
%   matrices, adaptive hyperviscosity, and global defect sweeps.

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
parser.addParameter('OutputPrefix', fullfile(pwd, 'moving_surface_adr_tp_combined_acceleration_results'));
parser.addParameter('WarmupParallelPool', true);
parser.addParameter('ContinueOnError', true);
parser.addParameter('SaveAfterEach', true);
parser.addParameter('DefectTolerance', 1.0e-6);
parser.addParameter('MaxDefectIterations', 4);
parser.addParameter('NeighborUpdateMode', "periodic");
parser.addParameter('NeighborSearchInterval', 5);
parser.addParameter('SpectrumCheck', false);
parser.addParameter('HyperviscosityDriftTolerance', 0.05);
parser.addParameter('MaxHyperviscositySkippedSteps', 5);
parser.addParameter('GlobalSolveMethod', "gmresIluDefect");
parser.addParameter('GlobalGMRESTolerance', 1.0e-6);
parser.addParameter('ExactStartup', false, @(x) islogical(x) && isscalar(x));
parser.addParameter('GlobalGMRESRestart', 40);
parser.addParameter('GlobalGMRESMaxIterations', 30);
parser.addParameter('BaselineGlobalDefectSweeps', 0);
parser.addParameter('AcceleratedGlobalDefectSweeps', 4);
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
            fprintf('\nCombined acceleration case=%s xi=%d N=%d\n', caseName, xi, N);

            baseline = runTimed(caseName, xi, N, dtScale, "baseline", parser.Results);
            accelerated = runTimed(caseName, xi, N, dtScale, "accelerated", parser.Results);

            rowIndex = rowIndex + 1;
            rows(rowIndex, 1) = makeRow(caseName, xi, N, baseline, accelerated);
            if rows(rowIndex).status == "ok"
                fprintf('  speedup: %.3f\n', rows(rowIndex).speedup);
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

function out = runTimed(caseName, xi, N, dtScale, mode, opts)
tic;
try
    if mode == "baseline"
        diffMethod = "direct";
        hypMode = "everyStep";
        globalDefectSweeps = opts.BaselineGlobalDefectSweeps;
    else
        diffMethod = "defect";
        hypMode = "adaptiveGeometry";
        globalDefectSweeps = opts.AcceleratedGlobalDefectSweeps;
    end
    globalArgs = globalSolveArgs(opts);
    result = moving_surface_adr_tp_geometric_flow_suite(caseName, xi, N, dtScale, ...
        'DiffMatUpdateMethod', diffMethod, ...
        'DefectTolerance', opts.DefectTolerance, ...
        'MaxDefectIterations', opts.MaxDefectIterations, ...
        'NeighborUpdateMode', opts.NeighborUpdateMode, ...
        'NeighborSearchInterval', opts.NeighborSearchInterval, ...
        'SpectrumCheck', opts.SpectrumCheck, ...
        'HyperviscosityUpdateMode', hypMode, ...
        'HyperviscosityDriftTolerance', opts.HyperviscosityDriftTolerance, ...
        'MaxHyperviscositySkippedSteps', opts.MaxHyperviscositySkippedSteps, ...
        globalArgs{:}, ...
        'GlobalDefectSweeps', globalDefectSweeps, ...
        'MassCorrectionMode', opts.MassCorrectionMode, ...
        'ExactStartup', opts.ExactStartup, ...
        'WriteOutputs', false);
    elapsed = toc;

    R = result.results;
    stats = R.updateStats{1};
    out = struct();
    out.ok = true;
    out.elapsed = elapsed;
    out.err = R.relerr;
    out.massRelError = R.massRelError;
    out.balanceRelResidual = R.balanceRelResidual;
    out.h = R.h;
    out.dt = R.dt;
    out.nsteps = R.nsteps;
    out.stats = stats;
    out.solveStats = R.solveStats{1};
    out.errorIdentifier = "";
    out.errorMessage = "";
catch ME
    elapsed = toc;
    if ~opts.ContinueOnError
        rethrow(ME);
    end
    out = failedRun(elapsed, ME);
    fprintf('  %s failed: %s\n', mode, ME.message);
end
end

function out = failedRun(elapsed, ME)
out = struct();
out.ok = false;
out.elapsed = elapsed;
out.err = NaN;
out.massRelError = NaN;
out.balanceRelResidual = NaN;
out.h = NaN;
out.dt = NaN;
out.nsteps = 0;
out.stats = [];
out.solveStats = [];
out.errorIdentifier = string(ME.identifier);
out.errorMessage = string(ME.message);
end

function args = globalSolveArgs(opts)
args = { ...
    'GlobalSolveMethod', opts.GlobalSolveMethod, ...
    'GlobalGMRESTolerance', opts.GlobalGMRESTolerance, ...
    'GlobalGMRESRestart', opts.GlobalGMRESRestart, ...
    'GlobalGMRESMaxIterations', opts.GlobalGMRESMaxIterations, ...
    'GlobalILUDropTolerance', opts.GlobalILUDropTolerance, ...
    'GlobalILURefreshInterval', opts.GlobalILURefreshInterval, ...
    'GlobalILURefreshIterationThreshold', opts.GlobalILURefreshIterationThreshold};
end

function row = makeRow(caseName, xi, N, baseline, accelerated)
row = emptyRow();
row.caseName = string(caseName);
row.xi = xi;
row.N = N;
if baseline.ok
    row.h = baseline.h;
    row.dt = baseline.dt;
    row.nsteps = baseline.nsteps;
elseif accelerated.ok
    row.h = accelerated.h;
    row.dt = accelerated.dt;
    row.nsteps = accelerated.nsteps;
end
row.baselineTime = baseline.elapsed;
row.acceleratedTime = accelerated.elapsed;
if baseline.ok && accelerated.ok
    row.speedup = baseline.elapsed / max(accelerated.elapsed, eps);
end
row.baselineError = baseline.err;
row.acceleratedError = accelerated.err;
row.errorDifference = abs(accelerated.err - baseline.err);
row.baselineMassRelError = baseline.massRelError;
row.acceleratedMassRelError = accelerated.massRelError;
row.baselineBalanceRelResidual = baseline.balanceRelResidual;
row.acceleratedBalanceRelResidual = accelerated.balanceRelResidual;
if ~isempty(baseline.solveStats)
    row.baselineIluRefreshes = baseline.solveStats.iluRefreshes;
    row.baselineIluReuses = baseline.solveStats.iluReuses;
    row.baselineGmresIterations = baseline.solveStats.gmresIterations;
    row.baselineGmresFailures = baseline.solveStats.gmresFailures;
    row.baselineMaxSolveResidual = baseline.solveStats.maxRelativeResidual;
end
if ~isempty(accelerated.solveStats)
    row.acceleratedIluRefreshes = accelerated.solveStats.iluRefreshes;
    row.acceleratedIluReuses = accelerated.solveStats.iluReuses;
    row.acceleratedGmresIterations = accelerated.solveStats.gmresIterations;
    row.acceleratedGmresFailures = accelerated.solveStats.gmresFailures;
    row.acceleratedMaxSolveResidual = accelerated.solveStats.maxRelativeResidual;
end
if ~isempty(accelerated.stats)
    row.updateDirect = accelerated.stats.direct;
    row.updateDefect = accelerated.stats.defectCorrected;
    row.updateFallback = accelerated.stats.defectFailedRefactored;
    row.hyperviscosityFullUpdates = accelerated.stats.hyperviscosityFullUpdates;
    row.hyperviscosityPredictedUpdates = accelerated.stats.hyperviscosityPredictedUpdates;
    row.maxUpdateResidual = accelerated.stats.maxRelativeResidual;
end
row.baselineErrorIdentifier = baseline.errorIdentifier;
row.baselineErrorMessage = baseline.errorMessage;
row.acceleratedErrorIdentifier = accelerated.errorIdentifier;
row.acceleratedErrorMessage = accelerated.errorMessage;

if baseline.ok && accelerated.ok
    row.status = "ok";
elseif baseline.ok
    row.status = "accelerated_failed";
elseif accelerated.ok
    row.status = "baseline_failed";
else
    row.status = "both_failed";
end
end

function bench = buildBench(caseNames, xiVals, Nvals, dtScale, rows, opts)
bench = struct();
bench.caseNames = caseNames;
bench.xiVals = xiVals;
bench.Nvals = Nvals;
bench.dtScale = dtScale;
bench.defectTolerance = opts.DefectTolerance;
bench.maxDefectIterations = opts.MaxDefectIterations;
bench.hyperviscosityDriftTolerance = opts.HyperviscosityDriftTolerance;
bench.maxHyperviscositySkippedSteps = opts.MaxHyperviscositySkippedSteps;
bench.globalSolveMethod = string(opts.GlobalSolveMethod);
bench.globalGMRESTolerance = opts.GlobalGMRESTolerance;
bench.baselineGlobalDefectSweeps = opts.BaselineGlobalDefectSweeps;
bench.acceleratedGlobalDefectSweeps = opts.AcceleratedGlobalDefectSweeps;
bench.globalILUDropTolerance = opts.GlobalILUDropTolerance;
bench.globalILURefreshInterval = opts.GlobalILURefreshInterval;
bench.massCorrectionMode = string(opts.MassCorrectionMode);
bench.rows = rows;
bench.table = struct2table(rows);
bench.summary = summarizeRows(rows);
end

function saveBenchmark(bench, outputPrefix)
outputPrefix = string(outputPrefix);
save(outputPrefix + ".mat", 'bench');
writetable(bench.table, outputPrefix + ".csv");
end

function summary = summarizeRows(rows)
speedups = [rows.speedup];
speedups = speedups(isfinite(speedups) & speedups > 0);
if isempty(speedups)
    summary = struct('geomeanSpeedup', NaN, 'minSpeedup', NaN, 'maxSpeedup', NaN);
else
    summary = struct( ...
        'geomeanSpeedup', exp(mean(log(speedups))), ...
        'minSpeedup', min(speedups), ...
        'maxSpeedup', max(speedups));
end
end

function dispBenchmarkSummary(bench)
fprintf('\nMoving-surface combined-acceleration benchmark summary\n');
fprintf('  rows: %d\n', numel(bench.rows));
fprintf('  speedup geomean/min/max: %.3f / %.3f / %.3f\n', ...
    bench.summary.geomeanSpeedup, bench.summary.minSpeedup, bench.summary.maxSpeedup);
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
    'baselineTime', NaN, ...
    'acceleratedTime', NaN, ...
    'speedup', NaN, ...
    'baselineError', NaN, ...
    'acceleratedError', NaN, ...
    'errorDifference', NaN, ...
    'baselineMassRelError', NaN, ...
    'acceleratedMassRelError', NaN, ...
    'baselineBalanceRelResidual', NaN, ...
    'acceleratedBalanceRelResidual', NaN, ...
    'baselineIluRefreshes', 0, ...
    'baselineIluReuses', 0, ...
    'baselineGmresIterations', 0, ...
    'baselineGmresFailures', 0, ...
    'baselineMaxSolveResidual', NaN, ...
    'acceleratedIluRefreshes', 0, ...
    'acceleratedIluReuses', 0, ...
    'acceleratedGmresIterations', 0, ...
    'acceleratedGmresFailures', 0, ...
    'acceleratedMaxSolveResidual', NaN, ...
    'updateDirect', 0, ...
    'updateDefect', 0, ...
    'updateFallback', 0, ...
    'hyperviscosityFullUpdates', 0, ...
    'hyperviscosityPredictedUpdates', 0, ...
    'maxUpdateResidual', NaN, ...
    'baselineErrorIdentifier', "", ...
    'baselineErrorMessage', "", ...
    'acceleratedErrorIdentifier', "", ...
    'acceleratedErrorMessage', "");
end
