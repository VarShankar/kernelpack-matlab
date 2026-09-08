function bench = moving_surface_adr_tp_hyperviscosity_benchmark(caseNames, xiVals, Nvals, dtScale, varargin)
%MOVING_SURFACE_ADR_TP_HYPERVISCOSITY_BENCHMARK Every-step vs adaptive timing.
%   Compares full hyperviscosity coefficient recomputation at every time step
%   against the geometry-aware adaptive update, while keeping the
%   differentiation-matrix update path fixed.

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
parser.addParameter('OutputPrefix', fullfile(pwd, 'moving_surface_adr_tp_hyperviscosity_benchmark_results_xi4'));
parser.addParameter('WarmupParallelPool', true);
parser.addParameter('ContinueOnError', true);
parser.addParameter('SaveAfterEach', true);
parser.addParameter('DiffMatUpdateMethod', "direct");
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
parser.addParameter('GlobalDefectSweeps', 4);
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
            fprintf('\nHyperviscosity benchmark case=%s xi=%d N=%d\n', caseName, xi, N);

            everyStep = runTimed(caseName, xi, N, dtScale, "everyStep", parser.Results);
            adaptive = runTimed(caseName, xi, N, dtScale, "adaptiveGeometry", parser.Results);

            rowIndex = rowIndex + 1;
            rows(rowIndex, 1) = makeRow(caseName, xi, N, everyStep, adaptive);
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
    globalArgs = globalSolveArgs(opts);
    result = moving_surface_adr_tp_geometric_flow_suite(caseName, xi, N, dtScale, ...
        'DiffMatUpdateMethod', opts.DiffMatUpdateMethod, ...
        'DefectTolerance', opts.DefectTolerance, ...
        'MaxDefectIterations', opts.MaxDefectIterations, ...
        'NeighborUpdateMode', opts.NeighborUpdateMode, ...
        'NeighborSearchInterval', opts.NeighborSearchInterval, ...
        'SpectrumCheck', opts.SpectrumCheck, ...
        'HyperviscosityUpdateMode', mode, ...
        'HyperviscosityDriftTolerance', opts.HyperviscosityDriftTolerance, ...
        'MaxHyperviscositySkippedSteps', opts.MaxHyperviscositySkippedSteps, ...
        globalArgs{:}, ...
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
    out.h = R.h;
    out.dt = R.dt;
    out.nsteps = R.nsteps;
    out.massRelError = R.massRelError;
    out.balanceRelResidual = R.balanceRelResidual;
    out.fullUpdates = stats.hyperviscosityFullUpdates;
    out.predictedUpdates = stats.hyperviscosityPredictedUpdates;
    out.maxRelativeDrift = stats.maxHyperviscosityRelativeDrift;
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
out.h = NaN;
out.dt = NaN;
out.nsteps = 0;
out.massRelError = NaN;
out.balanceRelResidual = NaN;
out.fullUpdates = 0;
out.predictedUpdates = 0;
out.maxRelativeDrift = NaN;
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
    'GlobalDefectSweeps', opts.GlobalDefectSweeps, ...
    'GlobalILUDropTolerance', opts.GlobalILUDropTolerance, ...
    'GlobalILURefreshInterval', opts.GlobalILURefreshInterval, ...
    'GlobalILURefreshIterationThreshold', opts.GlobalILURefreshIterationThreshold};
end

function row = makeRow(caseName, xi, N, everyStep, adaptive)
row = emptyRow();
row.caseName = string(caseName);
row.xi = xi;
row.N = N;
if everyStep.ok
    row.h = everyStep.h;
    row.dt = everyStep.dt;
    row.nsteps = everyStep.nsteps;
elseif adaptive.ok
    row.h = adaptive.h;
    row.dt = adaptive.dt;
    row.nsteps = adaptive.nsteps;
end
row.everyStepTime = everyStep.elapsed;
row.adaptiveTime = adaptive.elapsed;
if everyStep.ok && adaptive.ok
    row.speedup = everyStep.elapsed / max(adaptive.elapsed, eps);
end
row.everyStepError = everyStep.err;
row.adaptiveError = adaptive.err;
row.errorDifference = abs(adaptive.err - everyStep.err);
row.everyStepMassRelError = everyStep.massRelError;
row.adaptiveMassRelError = adaptive.massRelError;
row.everyStepBalanceRelResidual = everyStep.balanceRelResidual;
row.adaptiveBalanceRelResidual = adaptive.balanceRelResidual;
row.everyStepFullUpdates = everyStep.fullUpdates;
row.everyStepPredictedUpdates = everyStep.predictedUpdates;
row.adaptiveFullUpdates = adaptive.fullUpdates;
row.adaptivePredictedUpdates = adaptive.predictedUpdates;
row.adaptiveMaxRelativeDrift = adaptive.maxRelativeDrift;
if ~isempty(everyStep.solveStats)
    row.everyStepIluRefreshes = everyStep.solveStats.iluRefreshes;
    row.everyStepIluReuses = everyStep.solveStats.iluReuses;
    row.everyStepGmresIterations = everyStep.solveStats.gmresIterations;
    row.everyStepGmresFailures = everyStep.solveStats.gmresFailures;
    row.everyStepMaxSolveResidual = everyStep.solveStats.maxRelativeResidual;
end
if ~isempty(adaptive.solveStats)
    row.adaptiveIluRefreshes = adaptive.solveStats.iluRefreshes;
    row.adaptiveIluReuses = adaptive.solveStats.iluReuses;
    row.adaptiveGmresIterations = adaptive.solveStats.gmresIterations;
    row.adaptiveGmresFailures = adaptive.solveStats.gmresFailures;
    row.adaptiveMaxSolveResidual = adaptive.solveStats.maxRelativeResidual;
end
row.everyStepErrorIdentifier = everyStep.errorIdentifier;
row.everyStepErrorMessage = everyStep.errorMessage;
row.adaptiveErrorIdentifier = adaptive.errorIdentifier;
row.adaptiveErrorMessage = adaptive.errorMessage;

if everyStep.ok && adaptive.ok
    row.status = "ok";
elseif everyStep.ok
    row.status = "adaptive_failed";
elseif adaptive.ok
    row.status = "every_step_failed";
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
bench.diffMatUpdateMethod = string(opts.DiffMatUpdateMethod);
bench.hyperviscosityDriftTolerance = opts.HyperviscosityDriftTolerance;
bench.maxHyperviscositySkippedSteps = opts.MaxHyperviscositySkippedSteps;
bench.globalSolveMethod = string(opts.GlobalSolveMethod);
bench.globalGMRESTolerance = opts.GlobalGMRESTolerance;
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
if isempty(rows)
    summary.geomeanSpeedup = NaN;
    summary.minSpeedup = NaN;
    summary.maxSpeedup = NaN;
    return;
end
speedups = [rows.speedup];
speedups = speedups(isfinite(speedups) & speedups > 0);
if isempty(speedups)
    summary.geomeanSpeedup = NaN;
    summary.minSpeedup = NaN;
    summary.maxSpeedup = NaN;
    return;
end
summary.geomeanSpeedup = exp(mean(log(speedups)));
summary.minSpeedup = min(speedups);
summary.maxSpeedup = max(speedups);
end

function dispBenchmarkSummary(bench)
fprintf('\nMoving-surface hyperviscosity-update benchmark summary\n');
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
    'everyStepTime', NaN, ...
    'adaptiveTime', NaN, ...
    'speedup', NaN, ...
    'everyStepError', NaN, ...
    'adaptiveError', NaN, ...
    'errorDifference', NaN, ...
    'everyStepMassRelError', NaN, ...
    'adaptiveMassRelError', NaN, ...
    'everyStepBalanceRelResidual', NaN, ...
    'adaptiveBalanceRelResidual', NaN, ...
    'everyStepFullUpdates', 0, ...
    'everyStepPredictedUpdates', 0, ...
    'adaptiveFullUpdates', 0, ...
    'adaptivePredictedUpdates', 0, ...
    'adaptiveMaxRelativeDrift', NaN, ...
    'everyStepIluRefreshes', 0, ...
    'everyStepIluReuses', 0, ...
    'everyStepGmresIterations', 0, ...
    'everyStepGmresFailures', 0, ...
    'everyStepMaxSolveResidual', NaN, ...
    'adaptiveIluRefreshes', 0, ...
    'adaptiveIluReuses', 0, ...
    'adaptiveGmresIterations', 0, ...
    'adaptiveGmresFailures', 0, ...
    'adaptiveMaxSolveResidual', NaN, ...
    'everyStepErrorIdentifier', "", ...
    'everyStepErrorMessage', "", ...
    'adaptiveErrorIdentifier', "", ...
    'adaptiveErrorMessage', "");
end
