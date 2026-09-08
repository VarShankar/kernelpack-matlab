function bench = moving_surface_adr_tp_defect_benchmark(caseNames, xiVals, Nvals, dtScale, varargin)
%MOVING_SURFACE_ADR_TP_DEFECT_BENCHMARK Direct vs defect update timings.
%   Runs prescribed moving-surface ADR cases with both full direct
%   reassembly and the defect-corrected differentiation-matrix updater.
%   Results are saved to moving_surface_adr_tp_defect_benchmark_results.mat.

if nargin < 1 || isempty(caseNames)
    caseNames = ["mcf_sphere", "imcf_sphere", "rotating_breathing_sphere", ...
        "anisotropic_ellipsoid", "breathing_torus"];
end
if nargin < 2 || isempty(xiVals)
    xiVals = [2, 4, 6];
end
if nargin < 3 || isempty(Nvals)
    Nvals = paperNVals();
end
if nargin < 4 || isempty(dtScale)
    dtScale = 0.05;
end

parser = inputParser();
parser.addParameter('SaveResults', true);
parser.addParameter('OutputPrefix', fullfile(pwd, 'moving_surface_adr_tp_defect_benchmark_results'));
parser.addParameter('WarmupParallelPool', true);
parser.addParameter('ContinueOnError', true);
parser.addParameter('SaveAfterEach', true);
parser.addParameter('DefectTolerance', 1.0e-6);
parser.addParameter('MaxDefectIterations', 4);
parser.addParameter('NeighborUpdateMode', "periodic");
parser.addParameter('NeighborSearchInterval', 5);
parser.addParameter('SpectrumCheck', false);
parser.addParameter('HyperviscosityUpdateMode', "everyStep");
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
            fprintf('\nBenchmark case=%s xi=%d N=%d\n', caseName, xi, N);

            direct = runTimed(caseName, xi, N, dtScale, "direct", parser.Results);
            defect = runTimed(caseName, xi, N, dtScale, "defect", parser.Results);

            rowIndex = rowIndex + 1;
            rows(rowIndex, 1) = makeRow(caseName, xi, N, direct, defect);
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

function bench = buildBench(caseNames, xiVals, Nvals, dtScale, rows, opts)
bench = struct();
bench.caseNames = caseNames;
bench.xiVals = xiVals;
bench.Nvals = Nvals;
bench.dtScale = dtScale;
bench.defectTolerance = opts.DefectTolerance;
bench.maxDefectIterations = opts.MaxDefectIterations;
bench.neighborUpdateMode = string(opts.NeighborUpdateMode);
bench.neighborSearchInterval = opts.NeighborSearchInterval;
bench.spectrumCheck = opts.SpectrumCheck;
bench.hyperviscosityUpdateMode = string(opts.HyperviscosityUpdateMode);
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

function out = runTimed(caseName, xi, N, dtScale, method, opts)
tic;
try
    globalArgs = globalSolveArgs(opts);
    if method == "direct"
        result = moving_surface_adr_tp_geometric_flow_suite(caseName, xi, N, dtScale, ...
            'DiffMatUpdateMethod', 'direct', ...
            'SpectrumCheck', opts.SpectrumCheck, ...
            'HyperviscosityUpdateMode', opts.HyperviscosityUpdateMode, ...
            'HyperviscosityDriftTolerance', opts.HyperviscosityDriftTolerance, ...
            'MaxHyperviscositySkippedSteps', opts.MaxHyperviscositySkippedSteps, ...
            globalArgs{:}, ...
            'MassCorrectionMode', opts.MassCorrectionMode, ...
            'ExactStartup', opts.ExactStartup, ...
            'WriteOutputs', false);
    else
        result = moving_surface_adr_tp_geometric_flow_suite(caseName, xi, N, dtScale, ...
            'DiffMatUpdateMethod', 'defect', ...
            'DefectTolerance', opts.DefectTolerance, ...
            'MaxDefectIterations', opts.MaxDefectIterations, ...
            'NeighborUpdateMode', opts.NeighborUpdateMode, ...
            'NeighborSearchInterval', opts.NeighborSearchInterval, ...
            'SpectrumCheck', opts.SpectrumCheck, ...
            'HyperviscosityUpdateMode', opts.HyperviscosityUpdateMode, ...
            'HyperviscosityDriftTolerance', opts.HyperviscosityDriftTolerance, ...
            'MaxHyperviscositySkippedSteps', opts.MaxHyperviscositySkippedSteps, ...
            globalArgs{:}, ...
            'MassCorrectionMode', opts.MassCorrectionMode, ...
            'ExactStartup', opts.ExactStartup, ...
            'WriteOutputs', false);
    end
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
    out.errorIdentifier = "";
    out.errorMessage = "";
    out.solveStats = R.solveStats{1};
    if method == "defect"
        out.updateStats = R.updateStats{1};
    else
        out.updateStats = [];
    end
catch ME
    elapsed = toc;
    if ~opts.ContinueOnError
        rethrow(ME);
    end
    out = failedRun(elapsed, ME);
    fprintf('  %s failed: %s\n', method, ME.message);
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
out.updateStats = [];
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

function row = makeRow(caseName, xi, N, direct, defect)
row = emptyRow();
row.caseName = string(caseName);
row.xi = xi;
row.N = N;
if direct.ok
    row.h = direct.h;
    row.dt = direct.dt;
    row.nsteps = direct.nsteps;
elseif defect.ok
    row.h = defect.h;
    row.dt = defect.dt;
    row.nsteps = defect.nsteps;
end
row.directTime = direct.elapsed;
row.defectTime = defect.elapsed;
if direct.ok && defect.ok
    row.speedup = direct.elapsed / max(defect.elapsed, eps);
end
row.directError = direct.err;
row.defectError = defect.err;
row.errorDifference = abs(defect.err - direct.err);
row.directMassRelError = direct.massRelError;
row.defectMassRelError = defect.massRelError;
row.directBalanceRelResidual = direct.balanceRelResidual;
row.defectBalanceRelResidual = defect.balanceRelResidual;
row.directErrorIdentifier = direct.errorIdentifier;
row.directErrorMessage = direct.errorMessage;
row.defectErrorIdentifier = defect.errorIdentifier;
row.defectErrorMessage = defect.errorMessage;
if ~isempty(direct.solveStats)
    row.directIluRefreshes = direct.solveStats.iluRefreshes;
    row.directIluReuses = direct.solveStats.iluReuses;
    row.directGmresIterations = direct.solveStats.gmresIterations;
    row.directGmresFailures = direct.solveStats.gmresFailures;
    row.directMaxSolveResidual = direct.solveStats.maxRelativeResidual;
end
if ~isempty(defect.solveStats)
    row.defectIluRefreshes = defect.solveStats.iluRefreshes;
    row.defectIluReuses = defect.solveStats.iluReuses;
    row.defectGmresIterations = defect.solveStats.gmresIterations;
    row.defectGmresFailures = defect.solveStats.gmresFailures;
    row.defectMaxSolveResidual = defect.solveStats.maxRelativeResidual;
end

stats = defect.updateStats;
if direct.ok && defect.ok
    row.status = "ok";
elseif direct.ok
    row.status = "defect_failed";
elseif defect.ok
    row.status = "direct_failed";
else
    row.status = "both_failed";
end

if ~isempty(stats)
    row.updateDirect = stats.direct;
    row.updateDefect = stats.defectCorrected;
    row.updateFallback = stats.defectFailedRefactored;
    row.neighborSearches = stats.neighborSearches;
    row.neighborReuses = stats.neighborReuses;
    row.neighborReuseRejected = stats.neighborReuseRejected;
    row.maxUpdateResidual = stats.maxRelativeResidual;
    row.meanDefectIterations = stats.meanDefectIterations;
end
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

function T = benchmarkTable(rows)
if isempty(rows)
    T = table();
    return;
end

T = struct2table(rows);
T.caseName = string(T.caseName);
T = addRates(T, "directError", "directRate");
T = addRates(T, "defectError", "defectRate");
end

function T = addRates(T, errName, rateName)
T.(rateName) = nan(height(T), 1);
groups = findgroups(T.caseName, T.xi);
for g = 1:max(groups)
    idx = find(groups == g);
    idx = idx(isfinite(T.h(idx)) & T.h(idx) > 0 & ...
        isfinite(T.(errName)(idx)) & T.(errName)(idx) > 0);
    [~, order] = sort(T.h(idx), 'descend');
    idx = idx(order);
    for k = 2:numel(idx)
        T.(rateName)(idx(k)) = log(T.(errName)(idx(k - 1)) / T.(errName)(idx(k))) / ...
            log(T.h(idx(k - 1)) / T.h(idx(k)));
    end
end
end

function dispBenchmarkSummary(bench)
fprintf('\nMoving-surface defect-update benchmark summary\n');
fprintf('  rows: %d\n', numel(bench.rows));
fprintf('  speedup geomean/min/max: %.3f / %.3f / %.3f\n', ...
    bench.summary.geomeanSpeedup, bench.summary.minSpeedup, bench.summary.maxSpeedup);
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
    'directTime', NaN, ...
    'defectTime', NaN, ...
    'speedup', NaN, ...
    'directError', NaN, ...
    'defectError', NaN, ...
    'errorDifference', NaN, ...
    'directMassRelError', NaN, ...
    'defectMassRelError', NaN, ...
    'directBalanceRelResidual', NaN, ...
    'defectBalanceRelResidual', NaN, ...
    'directErrorIdentifier', "", ...
    'directErrorMessage', "", ...
    'defectErrorIdentifier', "", ...
    'defectErrorMessage', "", ...
    'directIluRefreshes', 0, ...
    'directIluReuses', 0, ...
    'directGmresIterations', 0, ...
    'directGmresFailures', 0, ...
    'directMaxSolveResidual', NaN, ...
    'defectIluRefreshes', 0, ...
    'defectIluReuses', 0, ...
    'defectGmresIterations', 0, ...
    'defectGmresFailures', 0, ...
    'defectMaxSolveResidual', NaN, ...
    'updateDirect', 0, ...
    'updateDefect', 0, ...
    'updateFallback', 0, ...
    'neighborSearches', 0, ...
    'neighborReuses', 0, ...
    'neighborReuseRejected', 0, ...
    'maxUpdateResidual', NaN, ...
    'meanDefectIterations', NaN);
end
