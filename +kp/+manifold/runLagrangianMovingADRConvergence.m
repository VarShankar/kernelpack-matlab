function results = runLagrangianMovingADRConvergence(problem, xiVals, Nvals, dtScale, opts)
%RUNLAGRANGIANMOVINGADRCONVERGENCE Prescribed moving-surface ADR study.
%   The geometry callback supplies material nodes X(t), outward normals, a
%   nominal spacing h, and any material coordinates needed by the exact
%   solution/forcing callbacks. The time step follows the source Lagrangian
%   BDF pattern: node labels move with the surface, while tangent-plane
%   differential operators are rebuilt on each new surface.

arguments
    problem (1,1) struct
    xiVals (1,:) double = [2, 4, 6]
    Nvals (1,:) double = [256, 576, 1024, 1600]
    dtScale (1,1) double = 0.5
    opts.Mu (1,1) double = 0.1
    opts.FinalTime (1,1) double = 0.12
    opts.HyperviscosityPower (1,1) double = NaN
    opts.StencilSize (1,1) double = NaN
    opts.StencilSizeFactor (1,1) double = 1
    opts.DiffMatUpdateMethod (1,1) string = "direct"
    opts.DefectTolerance (1,1) double = 1.0e-4
    opts.MaxDefectIterations (1,1) double = 4
    opts.NeighborUpdateMode (1,1) string = "periodic"
    opts.NeighborSearchInterval (1,1) double = 5
    opts.SpectrumCheck (1,1) logical = true
    opts.HyperviscosityUpdateMode (1,1) string = "everyStep"
    opts.HyperviscosityDriftTolerance (1,1) double = 0.05
    opts.MaxHyperviscositySkippedSteps (1,1) double = 5
    opts.NormalMode (1,1) string = "exact"
    opts.NormalNeighborCount (1,1) double = 32
    opts.GlobalSBFNormalDegree (1,1) double = 7
    opts.GlobalSBFControlPointCount (1,1) double = Inf
    opts.GlobalSBFNormalOrder (1,1) double = NaN
    opts.GlobalSBFBalanceSafety (1,1) double = 0.1
    opts.GlobalSBFControlPointScale (1,1) double = 1 / 3
    opts.GlobalSBFMinControlPointCount (1,1) double = 48
    opts.RecordMassDiagnostics (1,1) logical = true
    opts.ReturnSolution (1,1) logical = false
    opts.GlobalSolveMethod (1,1) string = "gmresIluDefect"
    opts.GlobalGMRESTolerance (1,1) double = NaN
    opts.GlobalGMRESRestart (1,1) double = 40
    opts.GlobalGMRESMaxIterations (1,1) double = 30
    opts.GlobalDefectSweeps (1,1) double = 4
    opts.GlobalILUDropTolerance (1,1) double = 1.0e-4
    opts.GlobalILURefreshInterval (1,1) double = Inf
    opts.GlobalILURefreshIterationThreshold (1,1) double = 30
    opts.FixedTimeStep (1,1) double = NaN
    opts.TrackTimeErrors (1,1) logical = false
    opts.ExactStartup (1,1) logical = false
    opts.FixedMassProjection (1,1) logical = false
    opts.FixedMassTarget (1,1) double = NaN
    opts.MassCorrectionMode (1,1) string = "off"
    opts.MassCorrectionMaxRelativeCorrection (1,1) double = Inf
    opts.UseRearrangement (1,1) logical = true
    opts.QualityThreshold (1,1) double = 1.75
    opts.PredictiveLookaheadSteps (1,1) double = 3
    opts.MinStepsBetweenRearrangements (1,1) double = 8
    opts.RearrangementTimes (1,:) double = []
    opts.RearrangementTransferMode (1,1) string = "localTp"
    opts.RearrangementSBFControlPointCount (1,1) double = Inf
    opts.ImagePath char = ''
    opts.ResultsPath char = ''
end

validateProblem(problem);

Nvals = Nvals(:).';
xiVals = xiVals(:).';
numLevels = numel(Nvals);
numXi = numel(xiVals);
results = repmat(struct( ...
    'xi', 0, ...
    'exactStartup', opts.ExactStartup, ...
    'N', zeros(1, numLevels), ...
    'h', zeros(1, numLevels), ...
    'dt', zeros(1, numLevels), ...
    'nsteps', zeros(1, numLevels), ...
    'relerr', zeros(1, numLevels), ...
    'rate', nan(1, numLevels), ...
    'initialMass', nan(1, numLevels), ...
    'finalMass', nan(1, numLevels), ...
        'exactFinalMass', nan(1, numLevels), ...
        'sourceIntegral', nan(1, numLevels), ...
        'massRelError', nan(1, numLevels), ...
        'balanceRelResidual', nan(1, numLevels), ...
        'massCorrectionSteps', nan(1, numLevels), ...
        'massCorrectionMaxAbs', nan(1, numLevels), ...
        'massCorrectionMaxRelative', nan(1, numLevels), ...
        'massCorrectionTotalAbs', nan(1, numLevels), ...
        'massCorrectionFinalAbs', nan(1, numLevels), ...
        'massCorrectionMaxPointShift', nan(1, numLevels), ...
        'normalRMSAngle', nan(1, numLevels), ...
    'normalMaxAngle', nan(1, numLevels), ...
    'normalControlPointCount', nan(1, numLevels), ...
    'normalControlPointFraction', nan(1, numLevels), ...
    'normalControlFillDistance', nan(1, numLevels), ...
    'finalAbsL2Error', nan(1, numLevels), ...
    'finalRelL2Error', nan(1, numLevels), ...
    'finalAbsLInfError', nan(1, numLevels), ...
    'finalRelLInfError', nan(1, numLevels), ...
    'maxAbsH1SemiError', nan(1, numLevels), ...
    'maxRelH1SemiError', nan(1, numLevels), ...
    'maxAbsH2SemiError', nan(1, numLevels), ...
    'maxRelH2SemiError', nan(1, numLevels), ...
    'maxAbsL2Error', nan(1, numLevels), ...
    'maxRelL2Error', nan(1, numLevels), ...
    'timeAbsL2L2Error', nan(1, numLevels), ...
    'timeRelL2L2Error', nan(1, numLevels), ...
    'numRearrangements', zeros(1, numLevels), ...
    'numPredictiveRearrangements', zeros(1, numLevels), ...
    'numThresholdRearrangements', zeros(1, numLevels), ...
    'numFixedRearrangements', zeros(1, numLevels), ...
    'initialQuality', nan(1, numLevels), ...
    'finalQuality', nan(1, numLevels), ...
    'maxQuality', nan(1, numLevels), ...
    'rearrangementTimes', {cell(1, numLevels)}, ...
    'updateStats', {cell(1, numLevels)}, ...
    'solveStats', {cell(1, numLevels)}), 1, numXi);

globalSolveOptions = struct( ...
    'method', opts.GlobalSolveMethod, ...
    'gmresTolerance', opts.GlobalGMRESTolerance, ...
    'gmresRestart', opts.GlobalGMRESRestart, ...
    'gmresMaxIterations', opts.GlobalGMRESMaxIterations, ...
    'defectSweeps', opts.GlobalDefectSweeps, ...
    'iluDropTolerance', opts.GlobalILUDropTolerance, ...
    'iluRefreshInterval', opts.GlobalILURefreshInterval, ...
    'iluRefreshIterationThreshold', opts.GlobalILURefreshIterationThreshold);

rearrangementOptions = struct( ...
    'enabled', opts.UseRearrangement, ...
    'qualityThreshold', opts.QualityThreshold, ...
    'predictiveLookaheadSteps', round(opts.PredictiveLookaheadSteps), ...
    'minStepsBetweenRearrangements', round(opts.MinStepsBetweenRearrangements), ...
    'fixedTimes', sort(opts.RearrangementTimes(:)), ...
    'transferMode', opts.RearrangementTransferMode, ...
    'sbfControlPointCount', opts.RearrangementSBFControlPointCount);

for ix = 1:numXi
    xi = xiVals(ix);
    results(ix).xi = xi;
    for level = 1:numLevels
        N = Nvals(level);
        geom0 = problem.geometry(N, 0.0);
        h = geom0.h;
        if isfinite(opts.FixedTimeStep) && opts.FixedTimeStep > 0
            nsteps = max(3, round(opts.FinalTime / opts.FixedTimeStep));
        else
            dtTarget = dtScale * h^(xi / 3);
            nsteps = max(3, ceil(opts.FinalTime / dtTarget));
        end
        dt = opts.FinalTime / nsteps;

        fprintf('%s\n', problem.label);
        fprintf('  xi: %d\n', xi);
        fprintf('  level %d / %d\n', level, numLevels);
        fprintf('  N: %d\n', N);
        fprintf('  h: %.6e\n', h);
        fprintf('  dt: %.6e\n', dt);
        fprintf('  steps: %d\n', nsteps);

        levelGlobalSolveOptions = globalSolveOptions;
        levelGlobalSolveOptions.gmresTolerance = selectGlobalGMRESTolerance( ...
            opts.GlobalGMRESTolerance, h, xi);
        fprintf('  GMRES tolerance: %.3e\n', levelGlobalSolveOptions.gmresTolerance);

        [err, updateStats, solveStats, diagnostics] = solveOneLevel(problem, geom0, xi, opts.Mu, opts.FinalTime, ...
            dt, nsteps, opts.HyperviscosityPower, ...
            opts.StencilSize, opts.StencilSizeFactor, ...
            opts.DiffMatUpdateMethod, ...
            opts.DefectTolerance, opts.MaxDefectIterations, opts.NeighborUpdateMode, ...
            opts.NeighborSearchInterval, opts.SpectrumCheck, opts.HyperviscosityUpdateMode, ...
            opts.HyperviscosityDriftTolerance, opts.MaxHyperviscositySkippedSteps, ...
            opts.NormalMode, opts.NormalNeighborCount, opts.GlobalSBFNormalDegree, ...
            opts.GlobalSBFControlPointCount, opts.GlobalSBFNormalOrder, ...
            opts.GlobalSBFBalanceSafety, opts.GlobalSBFControlPointScale, ...
            opts.GlobalSBFMinControlPointCount, ...
        opts.RecordMassDiagnostics, ...
        opts.ReturnSolution, opts.TrackTimeErrors, opts.ExactStartup, opts.FixedMassProjection, ...
            opts.FixedMassTarget, opts.MassCorrectionMode, ...
            opts.MassCorrectionMaxRelativeCorrection, ...
            levelGlobalSolveOptions, rearrangementOptions);
        results(ix).N(level) = N;
        results(ix).h(level) = h;
        results(ix).dt(level) = dt;
        results(ix).nsteps(level) = nsteps;
        results(ix).relerr(level) = err;
        results(ix).initialMass(level) = diagnostics.initialMass;
        results(ix).finalMass(level) = diagnostics.finalMass;
        results(ix).exactFinalMass(level) = diagnostics.exactFinalMass;
        results(ix).sourceIntegral(level) = diagnostics.sourceIntegral;
        results(ix).massRelError(level) = diagnostics.massRelError;
        results(ix).balanceRelResidual(level) = diagnostics.balanceRelResidual;
        results(ix).massCorrectionSteps(level) = diagnostics.massCorrectionSteps;
        results(ix).massCorrectionMaxAbs(level) = diagnostics.massCorrectionMaxAbs;
        results(ix).massCorrectionMaxRelative(level) = diagnostics.massCorrectionMaxRelative;
        results(ix).massCorrectionTotalAbs(level) = diagnostics.massCorrectionTotalAbs;
        results(ix).massCorrectionFinalAbs(level) = diagnostics.massCorrectionFinalAbs;
        results(ix).massCorrectionMaxPointShift(level) = diagnostics.massCorrectionMaxPointShift;
        results(ix).normalRMSAngle(level) = diagnostics.normalRMSAngle;
        results(ix).normalMaxAngle(level) = diagnostics.normalMaxAngle;
        results(ix).normalControlPointCount(level) = diagnostics.normalControlPointCount;
        results(ix).normalControlPointFraction(level) = diagnostics.normalControlPointFraction;
        results(ix).normalControlFillDistance(level) = diagnostics.normalControlFillDistance;
        results(ix).finalAbsL2Error(level) = diagnostics.finalAbsL2Error;
        results(ix).finalRelL2Error(level) = diagnostics.finalRelL2Error;
        results(ix).finalAbsLInfError(level) = diagnostics.finalAbsLInfError;
        results(ix).finalRelLInfError(level) = diagnostics.finalRelLInfError;
        results(ix).maxAbsH1SemiError(level) = diagnostics.maxAbsH1SemiError;
        results(ix).maxRelH1SemiError(level) = diagnostics.maxRelH1SemiError;
        results(ix).maxAbsH2SemiError(level) = diagnostics.maxAbsH2SemiError;
        results(ix).maxRelH2SemiError(level) = diagnostics.maxRelH2SemiError;
        results(ix).maxAbsL2Error(level) = diagnostics.maxAbsL2Error;
        results(ix).maxRelL2Error(level) = diagnostics.maxRelL2Error;
        results(ix).timeAbsL2L2Error(level) = diagnostics.timeAbsL2L2Error;
        results(ix).timeRelL2L2Error(level) = diagnostics.timeRelL2L2Error;
        results(ix).numRearrangements(level) = diagnostics.numRearrangements;
        results(ix).numPredictiveRearrangements(level) = diagnostics.numPredictiveRearrangements;
        results(ix).numThresholdRearrangements(level) = diagnostics.numThresholdRearrangements;
        results(ix).numFixedRearrangements(level) = diagnostics.numFixedRearrangements;
        results(ix).initialQuality(level) = diagnostics.initialQuality;
        results(ix).finalQuality(level) = diagnostics.finalQuality;
        results(ix).maxQuality(level) = diagnostics.maxQuality;
        results(ix).rearrangementTimes{level} = diagnostics.rearrangementTimes;
        if opts.ReturnSolution
            results(ix).solution{level} = diagnostics.solution;
            results(ix).finalGeometry{level} = diagnostics.finalGeometry;
            results(ix).finalMaterial{level} = diagnostics.finalMaterial;
        end
        results(ix).updateStats{level} = updateStats;
        results(ix).solveStats{level} = solveStats;
        fprintf('  relative L2 error: %.6e\n', err);
        fprintf('  global solves ilu/gmres/failures: %d / %d / %d\n', ...
            solveStats.iluRefreshes, solveStats.gmresSolves, solveStats.gmresFailures);
        fprintf('  global gmres iterations/max residual: %d / %.3e\n', ...
            solveStats.gmresIterations, solveStats.maxRelativeResidual);
        if updateStats.maxHyperviscosityPower > 0
            fprintf('  hyperviscosity power min/max: %d / %d\n', ...
                updateStats.minHyperviscosityPower, updateStats.maxHyperviscosityPower);
        end
        if strcmpi(opts.DiffMatUpdateMethod, "defect")
            fprintf('  update direct/defect/fallback: %d / %d / %d\n', ...
                updateStats.direct, updateStats.defectCorrected, ...
                updateStats.defectFailedRefactored);
            fprintf('  neighbor searches/reuses/rejections: %d / %d / %d\n', ...
                updateStats.neighborSearches, updateStats.neighborReuses, ...
                updateStats.neighborReuseRejected);
            fprintf('  update max residual: %.3e\n', updateStats.maxRelativeResidual);
        end
    end

    for k = 2:numLevels
        results(ix).rate(k) = log(results(ix).relerr(k - 1) / results(ix).relerr(k)) / ...
            log(results(ix).h(k - 1) / results(ix).h(k));
    end
end

if ~isempty(opts.ImagePath)
    writeConvergenceFigure(results, problem.title, opts.ImagePath);
end
if ~isempty(opts.ResultsPath)
    save(opts.ResultsPath, 'results', 'problem');
    [resultsFolder, resultsName] = fileparts(opts.ResultsPath);
    csvPath = fullfile(resultsFolder, [resultsName, '.csv']);
    writetable(kp.manifold.convergenceResultsTable(results), csvPath);
end
end

function op = configureSurfaceOperator(xi, stencilSize, stencilSizeFactor)
op = kp.manifold.rbffdop(2, xi, 2, 0);
defaultSize = op.stencilSize;
if isfinite(stencilSize) && stencilSize > 0
    requested = round(stencilSize);
else
    requested = round(stencilSizeFactor * defaultSize);
end
op.stencilSize = max(defaultSize, requested);
end

function [relerr, updateStats, solveStats, diagnostics] = solveOneLevel(problem, geom0, xi, mu, T, dt, nsteps, hyppow, ...
    stencilSize, stencilSizeFactor, updateMethod, defectTolerance, maxDefectIterations, neighborUpdateMode, neighborSearchInterval, ...
    spectrumCheck, hyperviscosityUpdateMode, hyperviscosityDriftTolerance, maxHyperviscositySkippedSteps, ...
    normalMode, normalNeighborCount, globalSBFNormalDegree, globalSBFControlPointCount, ...
    globalSBFNormalOrder, globalSBFBalanceSafety, globalSBFControlPointScale, ...
    globalSBFMinControlPointCount, recordMassDiagnostics, returnSolution, ...
    trackTimeErrors, exactStartup, fixedMassProjection, fixedMassTarget, massCorrectionMode, ...
    massCorrectionMaxRelativeCorrection, globalSolveOptions, rearrangementOptions)
op = configureSurfaceOperator(xi, stencilSize, stencilSizeFactor);
updater = [];
updateStats = defaultAccumulatedUpdateStats();
solveStats = defaultGlobalSolveStats();
if strcmpi(updateMethod, "defect")
    updater = kp.manifold.TangentPlaneDiffMatUpdater(xi, ...
        'Theta', 2, ...
        'StencilSize', op.stencilSize, ...
        'DefectTolerance', defectTolerance, ...
        'MaxDefectIterations', maxDefectIterations, ...
        'NeighborUpdateMode', neighborUpdateMode, ...
        'NeighborSearchInterval', neighborSearchInterval);
elseif ~strcmpi(updateMethod, "direct")
    error('kp:manifold:BadDiffMatUpdateMethod', ...
        'Unknown differentiation-matrix update method "%s".', updateMethod);
end

diagnostics = defaultRunDiagnostics();
normalState = [];
[geom0, normalStats, normalState] = prepareGeometryForOperators( ...
    geom0, normalMode, normalNeighborCount, globalSBFNormalDegree, ...
    globalSBFControlPointCount, globalSBFNormalOrder, globalSBFBalanceSafety, ...
    globalSBFControlPointScale, globalSBFMinControlPointCount, xi, normalState);
diagnostics = accumulateNormalDiagnostics(diagnostics, normalStats);
x0 = geom0.X;
material0 = geom0.material;
materialInitial = material0;
geomInitial = geom0;
if isfield(problem, 'initial')
    c0 = problem.initial(material0);
else
    c0 = problem.exact(0.0, material0);
end
cInitial = c0;
massState = initMassBalanceState(problem, geom0, material0, cInitial, mu, ...
    nsteps, fixedMassProjection, fixedMassTarget, massCorrectionMode, ...
    massCorrectionMaxRelativeCorrection);
errorState = initTimeErrorState(trackTimeErrors, nsteps);
errorState = recordTimeError(errorState, problem, material0, geom0, c0, 0.0);
I = speye(size(x0, 1), size(x0, 1));
hyperviscosityState = [];
globalSolveState = [];

qualityTimes = nan(nsteps + 1, 1);
qualityValues = nan(nsteps + 1, 1);
qualityTimes(1) = 0;
qualityValues(1) = pointQuality(x0);
maxQuality = qualityValues(1);
rearrangementTimes = zeros(0, 1);
rearrangementTriggerValues = zeros(0, 1);
rearrangementWasPredicted = false(0, 1);
rearrangementWasThreshold = false(0, 1);
rearrangementWasFixed = false(0, 1);
lastRearrangementStep = -inf;
fixedRearrangementTimes = rearrangementOptions.fixedTimes(:);
nextFixedRearrangement = 1;
useFixedSchedule = ~isempty(fixedRearrangementTimes);

material1 = advanceMaterial(problem, material0, 0.0, dt, size(x0, 1));
geom1Raw = geometryAt(problem, size(x0, 1), dt, material1);
material1 = geom1Raw.material;
x1 = geom1Raw.X;
v1 = (x1 - x0) / dt;
[geom1, normalStats, normalState] = prepareGeometryForOperators( ...
    geom1Raw, normalMode, normalNeighborCount, globalSBFNormalDegree, ...
    globalSBFControlPointCount, globalSBFNormalOrder, globalSBFBalanceSafety, ...
    globalSBFControlPointScale, globalSBFMinControlPointCount, xi, normalState);
diagnostics = accumulateNormalDiagnostics(diagnostics, normalStats);
[L, Gx, Gy, Gz, hypGamma, hypPower, hyperviscosityState, stats] = surfaceOpsAndHyperviscosity( ...
    geom1, op, hyppow, xi, updater, spectrumCheck, hyperviscosityUpdateMode, ...
    hyperviscosityState, hyperviscosityDriftTolerance, maxHyperviscositySkippedSteps);
updateStats = accumulateUpdateStats(updateStats, stats);
diffusionOp = mu * L;
lhs = I - dt * diffusionOp;
rhs = c0 + dt * forcingAt(problem, dt, material1, mu) ...
    - dt * c0 .* surfaceDivergence(Gx, Gy, Gz, v1) ...
    - dt * c0 .* hyperviscosityDivergenceCorrection( ...
        L, v1, hypGamma, hypPower);
if exactStartup && isfield(problem, 'exact')
    c1 = problem.exact(dt, material1);
else
    [c1, globalSolveState, solveStats] = solveGlobalBDFSystem( ...
        lhs, rhs, c0, 1.0, globalSolveState, solveStats, globalSolveOptions);
end
massState = advanceMassBalanceState(massState, problem, geom1, material1, mu, dt, 1);
[c1, massState] = applyMassCorrectionAtStep(c1, geom1, massState, 1, true);
errorState = recordTimeError(errorState, problem, material1, geom1, c1, dt);
errorState = recordOperatorTimeError(errorState, problem, material1, geom1, c1, dt, L, Gx, Gy, Gz);
qualityTimes(2) = dt;
qualityValues(2) = pointQuality(x1);
maxQuality = max(maxQuality, qualityValues(2));

material2 = advanceMaterial(problem, material1, dt, 2 * dt, size(x0, 1));
geom2Raw = geometryAt(problem, size(x0, 1), 2 * dt, material2);
material2 = geom2Raw.material;
x2 = geom2Raw.X;
v2 = (3 * x2 - 4 * x1 + x0) / (2 * dt);
[geom2, normalStats, normalState] = prepareGeometryForOperators( ...
    geom2Raw, normalMode, normalNeighborCount, globalSBFNormalDegree, ...
    globalSBFControlPointCount, globalSBFNormalOrder, globalSBFBalanceSafety, ...
    globalSBFControlPointScale, globalSBFMinControlPointCount, xi, normalState);
diagnostics = accumulateNormalDiagnostics(diagnostics, normalStats);
[L, Gx, Gy, Gz, hypGamma, hypPower, hyperviscosityState, stats] = surfaceOpsAndHyperviscosity( ...
    geom2, op, hyppow, xi, updater, spectrumCheck, hyperviscosityUpdateMode, ...
    hyperviscosityState, hyperviscosityDriftTolerance, maxHyperviscositySkippedSteps);
updateStats = accumulateUpdateStats(updateStats, stats);
diffusionOp = mu * L;
lhs = I - (2 / 3) * dt * diffusionOp;
c2a = 2 * c1 - c0;
rhs = (4 / 3) * c1 - (1 / 3) * c0 ...
    + (2 / 3) * dt * forcingAt(problem, 2 * dt, material2, mu) ...
    - (2 / 3) * dt * c2a .* surfaceDivergence(Gx, Gy, Gz, v2) ...
    - (2 / 3) * dt * c2a .* hyperviscosityDivergenceCorrection( ...
        L, v2, hypGamma, hypPower);
if exactStartup && isfield(problem, 'exact')
    c2 = problem.exact(2 * dt, material2);
else
    [c2, globalSolveState, solveStats] = solveGlobalBDFSystem( ...
        lhs, rhs, c2a, 2 / 3, globalSolveState, solveStats, globalSolveOptions);
end
massState = advanceMassBalanceState(massState, problem, geom2, material2, mu, dt, 2);
[c2, massState] = applyMassCorrectionAtStep(c2, geom2, massState, 2, true);
errorState = recordTimeError(errorState, problem, material2, geom2, c2, 2 * dt);
errorState = recordOperatorTimeError(errorState, problem, material2, geom2, c2, 2 * dt, L, Gx, Gy, Gz);
qualityTimes(3) = 2 * dt;
qualityValues(3) = pointQuality(x2);
maxQuality = max(maxQuality, qualityValues(3));

for step = 3:nsteps
    tnow = step * dt;
    material3 = advanceMaterial(problem, material2, tnow - dt, tnow, size(x0, 1));
    geom3Raw = geometryAt(problem, size(x0, 1), tnow, material3);
    material3 = geom3Raw.material;
    x3 = geom3Raw.X;
    v3 = (11 * x3 - 18 * x2 + 9 * x1 - 2 * x0) / (6 * dt);
    [geom3, normalStats, normalState] = prepareGeometryForOperators( ...
        geom3Raw, normalMode, normalNeighborCount, globalSBFNormalDegree, ...
        globalSBFControlPointCount, globalSBFNormalOrder, globalSBFBalanceSafety, ...
        globalSBFControlPointScale, globalSBFMinControlPointCount, xi, normalState);
    diagnostics = accumulateNormalDiagnostics(diagnostics, normalStats);
    [L, Gx, Gy, Gz, hypGamma, hypPower, hyperviscosityState, stats] = surfaceOpsAndHyperviscosity( ...
        geom3, op, hyppow, xi, updater, spectrumCheck, hyperviscosityUpdateMode, ...
        hyperviscosityState, hyperviscosityDriftTolerance, maxHyperviscositySkippedSteps);
    updateStats = accumulateUpdateStats(updateStats, stats);

    diffusionOp = mu * L;
    lhs = I - (6 / 11) * dt * diffusionOp;
    c3a = 3 * c2 - 3 * c1 + c0;
    rhs = (18 / 11) * c2 - (9 / 11) * c1 + (2 / 11) * c0 ...
        + (6 / 11) * dt * forcingAt(problem, tnow, material3, mu) ...
        - (6 / 11) * dt * c3a .* surfaceDivergence(Gx, Gy, Gz, v3) ...
        - (6 / 11) * dt * c3a .* hyperviscosityDivergenceCorrection( ...
            L, v3, hypGamma, hypPower);
    [c3, globalSolveState, solveStats] = solveGlobalBDFSystem( ...
        lhs, rhs, c3a, 6 / 11, globalSolveState, solveStats, globalSolveOptions);
    massState = advanceMassBalanceState(massState, problem, geom3, material3, mu, dt, step);
    [c3, massState] = applyMassCorrectionAtStep(c3, geom3, massState, step, true);
    errorState = recordTimeError(errorState, problem, material3, geom3, c3, tnow);
    errorState = recordOperatorTimeError(errorState, problem, material3, geom3, c3, tnow, L, Gx, Gy, Gz);

    qBeforeRearrangement = pointQuality(x3);
    maxQuality = max(maxQuality, qBeforeRearrangement);
    [willCross, triggerMetric] = predictQualityCrossing( ...
        qualityValues, step, qBeforeRearrangement, ...
        rearrangementOptions.qualityThreshold, ...
        rearrangementOptions.predictiveLookaheadSteps);
    fixedTrigger = useFixedSchedule && nextFixedRearrangement <= numel(fixedRearrangementTimes) && ...
        tnow >= fixedRearrangementTimes(nextFixedRearrangement) - 0.5 * dt;
    thresholdTrigger = ~useFixedSchedule && qBeforeRearrangement > rearrangementOptions.qualityThreshold;
    predictiveTrigger = ~useFixedSchedule && ~thresholdTrigger && willCross;
    canRearrange = rearrangementOptions.enabled && step >= 3 && ...
        step - lastRearrangementStep >= rearrangementOptions.minStepsBetweenRearrangements && ...
        (fixedTrigger || thresholdTrigger || predictiveTrigger);

    if canRearrange
        [geom0, c0, geom1, c1, geom2, c2] = rearrangeBDFHistory( ...
            problem, size(x0, 1), xi, rearrangementOptions, ...
            geom1, c1, geom2, c2, geom3, c3, tnow, dt);
        [c0, massState] = applyMassCorrectionAtStep(c0, geom0, massState, step - 2, true);
        [c1, massState] = applyMassCorrectionAtStep(c1, geom1, massState, step - 1, true);
        [c2, massState] = applyMassCorrectionAtStep(c2, geom2, massState, step, true);
        x0 = geom0.X; x1 = geom1.X; x2 = geom2.X;
        material2 = geom2.material;
        if ~isempty(updater)
            updater.reset();
        end
        hyperviscosityState = [];
        globalSolveState = [];
        normalState = resetNormalStateForRearrangement(normalState);
        lastRearrangementStep = step;
        if fixedTrigger
            nextFixedRearrangement = nextFixedRearrangement + 1;
        end
        rearrangementTimes(end + 1, 1) = tnow; %#ok<AGROW>
        rearrangementTriggerValues(end + 1, 1) = triggerMetric; %#ok<AGROW>
        rearrangementWasPredicted(end + 1, 1) = predictiveTrigger; %#ok<AGROW>
        rearrangementWasThreshold(end + 1, 1) = thresholdTrigger; %#ok<AGROW>
        rearrangementWasFixed(end + 1, 1) = fixedTrigger; %#ok<AGROW>
        qualityValues(step + 1) = pointQuality(x2);
    else
        geom1 = geom2; geom2 = geom3;
        material2 = material3;
        c0 = c1; c1 = c2; c2 = c3;
        x0 = x1; x1 = x2; x2 = x3;
        qualityValues(step + 1) = qBeforeRearrangement;
    end
    qualityTimes(step + 1) = tnow;
    maxQuality = max(maxQuality, qualityValues(step + 1));
end

if isfield(problem, 'exact')
    cex = problem.exact(T, material2);
    relerr = norm(c2 - cex) / max(norm(cex), 1.0e-14);
else
    cex = NaN(size(c2));
    relerr = NaN;
end
if recordMassDiagnostics
    diagnostics = computeMassDiagnostics(diagnostics, problem, materialInitial, mu, ...
        geomInitial, geom2, cInitial, c2, cex, T, dt, nsteps, massState);
end
diagnostics = finalizeTimeErrorDiagnostics(diagnostics, errorState, c2, cex, geom2);
diagnostics.numRearrangements = numel(rearrangementTimes);
diagnostics.numPredictiveRearrangements = nnz(rearrangementWasPredicted);
diagnostics.numThresholdRearrangements = nnz(rearrangementWasThreshold);
diagnostics.numFixedRearrangements = nnz(rearrangementWasFixed);
diagnostics.initialQuality = qualityValues(1);
diagnostics.finalQuality = qualityValues(nsteps + 1);
diagnostics.maxQuality = maxQuality;
diagnostics.qualityTimes = qualityTimes;
diagnostics.qualityValues = qualityValues;
diagnostics.rearrangementTimes = rearrangementTimes;
diagnostics.rearrangementTriggerValues = rearrangementTriggerValues;
if returnSolution
    diagnostics.solution = c2;
    diagnostics.finalGeometry = geom2;
    diagnostics.finalMaterial = material2;
end
end

function geom = geometryAt(problem, N, t, material)
if nargin >= 4 && ~isempty(material) && isfield(problem, 'geometryFromMaterial')
    geom = problem.geometryFromMaterial(material, t);
else
    geom = problem.geometry(N, t);
end
if ~isfield(geom, 'material') || isempty(geom.material)
    geom.material = material;
end
if isfield(geom, 'material') && isstruct(geom.material)
    geom.material.X = geom.X;
end
end

function materialNext = advanceMaterial(problem, material, t0, t1, N)
if isfield(problem, 'advanceMaterial')
    materialNext = problem.advanceMaterial(material, t0, t1);
elseif isfield(problem, 'geometryFromMaterial')
    materialNext = material;
else
    geom = problem.geometry(N, t1);
    materialNext = geom.material;
end
if isstruct(materialNext) && isstruct(material) && ~isfield(materialNext, 'X') && isfield(material, 'X')
    materialNext.X = material.X;
end
end

function materialPast = backtraceMaterial(problem, materialNow, tNow, tPast)
if isfield(problem, 'backtraceMaterial')
    materialPast = problem.backtraceMaterial(materialNow, tNow, tPast);
else
    materialPast = materialNow;
end
end

function f = forcingAt(problem, t, material, mu)
f = problem.forcing(t, material, mu);
end

function [geom0New, c0New, geom1New, c1New, geom2New, c2New] = rearrangeBDFHistory( ...
    problem, N, xi, rearrangementOptions, geom0Old, c0Old, geom1Old, c1Old, ...
    geom2Old, c2Old, tCurrent, dt)
if ~isfield(problem, 'sampleMaterial')
    error('kp:manifold:MissingRearrangementSampler', ...
        'Rearrangement is enabled, but the problem does not provide sampleMaterial.');
end
if ~isfield(problem, 'geometryFromMaterial')
    error('kp:manifold:MissingRearrangementGeometry', ...
        'Rearrangement is enabled, but the problem does not provide geometryFromMaterial.');
end

material2New = problem.sampleMaterial(N, tCurrent);
material1New = backtraceMaterial(problem, material2New, tCurrent, tCurrent - dt);
material0New = backtraceMaterial(problem, material2New, tCurrent, tCurrent - 2 * dt);
geom2New = geometryAt(problem, N, tCurrent, material2New);
geom1New = geometryAt(problem, N, tCurrent - dt, material1New);
geom0New = geometryAt(problem, N, tCurrent - 2 * dt, material0New);

switch lower(string(rearrangementOptions.transferMode))
    case {"localtp", "localtangentplane"}
        c2New = localTangentPlaneScalarTransfer(geom2Old, c2Old, geom2New, xi);
        c1New = localTangentPlaneScalarTransfer(geom1Old, c1Old, geom1New, xi);
        c0New = localTangentPlaneScalarTransfer(geom0Old, c0Old, geom0New, xi);
    case {"sbf", "globalsbf"}
        c2New = scalarSBFTransfer(geom2Old, c2Old, geom2New, rearrangementOptions.sbfControlPointCount);
        c1New = scalarSBFTransfer(geom1Old, c1Old, geom1New, rearrangementOptions.sbfControlPointCount);
        c0New = scalarSBFTransfer(geom0Old, c0Old, geom0New, rearrangementOptions.sbfControlPointCount);
    otherwise
        error('kp:manifold:BadRearrangementTransferMode', ...
            'Unknown rearrangement transfer mode "%s".', rearrangementOptions.transferMode);
end

end

function values = localTangentPlaneScalarTransfer(sourceGeom, sourceValues, targetGeom, xi)
op = kp.manifold.rbffdop(2, xi, 2, 0);
stencilSize = min(op.stencilSize, size(sourceGeom.X, 1));
tree = KDTreeSearcher(sourceGeom.X);
idx = knnsearch(tree, targetGeom.X, 'K', stencilSize);
values = zeros(size(targetGeom.X, 1), 1);
polyBasis = kp.poly.PolynomialBasis.fromTotalDegree(2, op.ell, 'Family', 'legendre');
targetNormals = transferNormals(targetGeom);
for i = 1:size(targetGeom.X, 1)
    stencil = idx(i, :);
    normal = targetNormals(i, :).';
    R = tangentBasis(normal(:));
    local = (sourceGeom.X(stencil, :) - targetGeom.X(i, :)) * R;
    width = max(abs(local), [], 'all');
    if width <= eps
        values(i) = sourceValues(stencil(1));
        continue;
    end
    localScaled = local ./ width;
    r = kp.geometry.distanceMatrix(localScaled, localScaled);
    phi = (r + eps) .^ op.rbfexp;
    P = polyBasis.evaluate(localScaled, zeros(1, 2), true);
    A = [[phi, P]; [P.', zeros(size(P, 2))]];
    rhs = [sourceValues(stencil); zeros(size(P, 2), 1)];
    coeff = kp.manifold.detail.solveLocalAugmentedSystem(A, rhs);
    rq = sqrt(sum(localScaled .^ 2, 2)).';
    phiq = (rq + eps) .^ op.rbfexp;
    pq = polyBasis.evaluate([0, 0], zeros(1, 2), true);
    values(i) = [phiq, pq] * coeff;
end
end

function values = scalarSBFTransfer(sourceGeom, sourceValues, targetGeom, controlPointCount)
if ~isfield(sourceGeom, 'material') || ~isfield(targetGeom, 'material') || ...
        ~isfield(sourceGeom.material, 'U') || ~isfield(targetGeom.material, 'U')
    error('kp:manifold:SBFTransferNeedsSphereMaterial', ...
        'SBF rearrangement transfer currently needs material.U on source and target geometries.');
end
U = sourceGeom.material.U;
Uq = targetGeom.material.U;
n = size(U, 1);
if ~isfinite(controlPointCount)
    controlPointCount = n;
else
    controlPointCount = min(n, max(1, round(controlPointCount)));
end
ids = farthestPointSubset(U ./ max(vecnorm(U, 2, 2), eps), controlPointCount);
degree = 7;
centers = U(ids, :) ./ max(vecnorm(U(ids, :), 2, 2), eps);
targets = Uq ./ max(vecnorm(Uq, 2, 2), eps);
[r, ~] = kp.geometry.sphereChordDistance(centers, centers);
kernel = kp.geometry.phsKernel(r, degree);
reg = 1.0e-12 * max(1.0, max(abs(kernel), [], 'all'));
coeff = (kernel + reg * eye(numel(ids))) \ sourceValues(ids);
[rq, ~] = kp.geometry.sphereChordDistance(targets, centers);
values = kp.geometry.phsKernel(rq, degree) * coeff;
end

function normals = transferNormals(geom)
if isfield(geom, 'normals') && ~isempty(geom.normals)
    normals = geom.normals;
elseif isfield(geom, 'material') && ...
        ((isfield(geom.material, 'U') && size(geom.material.U, 2) == 3) || ...
        (isfield(geom.material, 'theta') && isfield(geom.material, 'phi')))
    count = min(size(geom.X, 1), max(48, ceil(size(geom.X, 1) / 3)));
    normals = kp.manifold.estimateNormalsGlobalSBF(geom, [], ...
        'Degree', 7, ...
        'ControlPointCount', count);
else
    normals = kp.manifold.estimateNormalsPCA(geom.X, 'NumNeighbors', 32);
end
normals = kp.geometry.normalizeRows(normals);
end

function R = tangentBasis(normal)
normal = normal / max(norm(normal), eps);
if abs(normal(3)) < 0.9
    seed = [0; 0; 1];
else
    seed = [1; 0; 0];
end
t1 = cross(normal, seed);
t1 = t1 / max(norm(t1), eps);
t2 = cross(normal, t1);
R = [t1, t2];
end

function ids = farthestPointSubset(Y, count)
n = size(Y, 1);
count = min(n, max(1, round(count)));
ids = zeros(count, 1);
[~, ids(1)] = max(sum((Y - mean(Y, 1)) .^ 2, 2));
dist2 = inf(n, 1);
for k = 2:count
    d2 = sum((Y - Y(ids(k - 1), :)) .^ 2, 2);
    dist2 = min(dist2, d2);
    [~, ids(k)] = max(dist2);
end
ids = sort(ids);
end

function q = pointQuality(X)
tree = KDTreeSearcher(X);
[~, d] = knnsearch(tree, X, 'K', 2);
nearest = d(:, 2);
q = max(nearest) / max(min(nearest), eps);
end

function [willCross, triggerMetric] = predictQualityCrossing( ...
    qualityValues, currentStep, currentQuality, threshold, lookaheadSteps)
previous = qualityValues(max(1, currentStep));
if ~isfinite(previous)
    previous = currentQuality;
end
slope = max(0, currentQuality - previous);
triggerMetric = currentQuality + max(0, lookaheadSteps) * slope;
willCross = triggerMetric >= threshold;
end

function normalState = resetNormalStateForRearrangement(normalState)
if isstruct(normalState) && isfield(normalState, 'calibratedControlPointCount')
    count = normalState.calibratedControlPointCount;
    normalState = struct('calibratedControlPointCount', count);
else
    normalState = [];
end
end

function [geom, normalStats, normalState] = prepareGeometryForOperators( ...
    geom, normalMode, normalNeighborCount, globalSBFNormalDegree, ...
    globalSBFControlPointCount, globalSBFNormalOrder, globalSBFBalanceSafety, ...
    globalSBFControlPointScale, globalSBFMinControlPointCount, xi, normalState)
normalStats = struct( ...
    'rmsAngleDegrees', NaN, ...
    'maxAngleDegrees', NaN, ...
    'controlPointCount', NaN, ...
    'controlPointFraction', NaN, ...
    'controlFillDistance', NaN);
mode = lower(string(normalMode));
if mode == "exact" || mode == "provided" || mode == "geometry"
    if isfield(geom, 'normals') && ~isempty(geom.normals)
        normalStats.rmsAngleDegrees = 0;
        normalStats.maxAngleDegrees = 0;
    end
    return;
elseif mode == "pca" || mode == "estimated"
    ref = [];
    if isfield(geom, 'normals') && ~isempty(geom.normals)
        ref = geom.normals;
    end
    [geom.normals, info] = kp.manifold.estimateNormalsPCA(geom.X, ...
        'NumNeighbors', normalNeighborCount, ...
        'ReferenceNormals', ref);
    if isfield(info, 'rmsAngleDegrees')
        normalStats.rmsAngleDegrees = info.rmsAngleDegrees;
    end
    if isfield(info, 'maxAngleDegrees')
        normalStats.maxAngleDegrees = info.maxAngleDegrees;
    end
elseif mode == "globalsbf" || mode == "sbf" || mode == "parametricsbf" || ...
        mode == "calibratedsbf" || mode == "calibratedglobalsbf"
    if mode == "calibratedsbf" || mode == "calibratedglobalsbf"
        if isstruct(normalState) && isfield(normalState, 'calibratedControlPointCount') && ...
                ~isempty(normalState.calibratedControlPointCount)
            globalSBFControlPointCount = normalState.calibratedControlPointCount;
        else
            globalSBFControlPointCount = calibratedSBFControlPointCount(geom, xi, ...
                globalSBFNormalDegree, globalSBFNormalOrder, globalSBFBalanceSafety, ...
                globalSBFControlPointScale, globalSBFMinControlPointCount);
        end
    end
    [geom.normals, info, normalState] = kp.manifold.estimateNormalsGlobalSBF( ...
        geom, normalState, ...
        'Degree', globalSBFNormalDegree, ...
        'ControlPointCount', globalSBFControlPointCount);
    if isfield(info, 'rmsAngleDegrees')
        normalStats.rmsAngleDegrees = info.rmsAngleDegrees;
    end
    if isfield(info, 'maxAngleDegrees')
        normalStats.maxAngleDegrees = info.maxAngleDegrees;
    end
    if isfield(info, 'controlPointCount')
        normalStats.controlPointCount = info.controlPointCount;
    end
    if isfield(info, 'controlPointFraction')
        normalStats.controlPointFraction = info.controlPointFraction;
    end
    if isfield(info, 'controlFillDistance')
        normalStats.controlFillDistance = info.controlFillDistance;
    end
    if mode == "calibratedsbf" || mode == "calibratedglobalsbf"
        normalState.calibratedControlPointCount = globalSBFControlPointCount;
    end
else
    error('kp:manifold:BadNormalMode', ...
        ['Unknown normal mode "%s". Expected "provided", "exact", "pca", ' ...
        '"globalSbf", or "calibratedSbf".'], normalMode);
end
end

function count = calibratedSBFControlPointCount(geom, xi, degree, normalOrder, safety, scale, minCount)
count = kp.manifold.calibratedSBFControlPointCount(geom, xi, ...
    'Degree', degree, ...
    'NormalOrder', normalOrder, ...
    'BalanceSafety', safety, ...
    'ControlPointScale', scale, ...
    'MinControlPointCount', minCount);
end

function diagnostics = defaultRunDiagnostics()
diagnostics = struct( ...
    'initialMass', NaN, ...
    'finalMass', NaN, ...
    'exactFinalMass', NaN, ...
    'sourceIntegral', NaN, ...
    'massRelError', NaN, ...
    'balanceRelResidual', NaN, ...
    'massCorrectionSteps', NaN, ...
    'massCorrectionMaxAbs', NaN, ...
    'massCorrectionMaxRelative', NaN, ...
    'massCorrectionTotalAbs', NaN, ...
    'massCorrectionFinalAbs', NaN, ...
    'massCorrectionMaxPointShift', NaN, ...
    'normalRMSAngle', NaN, ...
    'normalMaxAngle', NaN, ...
    'normalControlPointCount', NaN, ...
    'normalControlPointFraction', NaN, ...
    'normalControlFillDistance', NaN, ...
    'finalAbsL2Error', NaN, ...
    'finalRelL2Error', NaN, ...
    'finalAbsLInfError', NaN, ...
    'finalRelLInfError', NaN, ...
    'maxAbsH1SemiError', NaN, ...
    'maxRelH1SemiError', NaN, ...
    'maxAbsH2SemiError', NaN, ...
    'maxRelH2SemiError', NaN, ...
    'maxAbsL2Error', NaN, ...
    'maxRelL2Error', NaN, ...
    'timeAbsL2L2Error', NaN, ...
    'timeRelL2L2Error', NaN, ...
    'numRearrangements', 0, ...
    'numPredictiveRearrangements', 0, ...
    'numThresholdRearrangements', 0, ...
    'numFixedRearrangements', 0, ...
    'initialQuality', NaN, ...
    'finalQuality', NaN, ...
    'maxQuality', NaN, ...
    'qualityTimes', [], ...
    'qualityValues', [], ...
    'rearrangementTimes', [], ...
    'rearrangementTriggerValues', [], ...
    'solution', []);
end

function state = initTimeErrorState(enabled, nsteps)
state = struct( ...
    'enabled', enabled, ...
    'times', nan(nsteps + 1, 1), ...
    'absL2', nan(nsteps + 1, 1), ...
    'relL2', nan(nsteps + 1, 1), ...
    'absH1Semi', nan(nsteps + 1, 1), ...
    'relH1Semi', nan(nsteps + 1, 1), ...
    'absH2Semi', nan(nsteps + 1, 1), ...
    'relH2Semi', nan(nsteps + 1, 1), ...
    'count', 0);
end

function state = recordTimeError(state, problem, material, geom, c, t)
if ~state.enabled || ~isfield(problem, 'exact')
    return;
end
cex = problem.exact(t, material);
if ~all(isfinite(cex))
    return;
end
w = surfaceWeights(geom);
err = c(:) - cex(:);
absL2 = sqrt(max(sum(w .* abs(err).^2), 0));
exactL2 = sqrt(max(sum(w .* abs(cex(:)).^2), 0));
state.count = state.count + 1;
state.times(state.count) = t;
state.absL2(state.count) = absL2;
state.relL2(state.count) = absL2 / max(exactL2, 1.0e-14);
end

function state = recordOperatorTimeError(state, problem, material, geom, c, t, L, Gx, Gy, Gz)
if ~state.enabled || state.count == 0
    return;
end
idx = state.count;
w = surfaceWeights(geom);
if isfield(problem, 'exactGradient')
    exactGrad = problem.exactGradient(t, material);
    if all(isfinite(exactGrad), 'all')
        numericalGrad = [Gx * c(:), Gy * c(:), Gz * c(:)];
        gradErr = numericalGrad - exactGrad;
        absH1 = sqrt(max(sum(w .* sum(abs(gradErr).^2, 2)), 0));
        exactH1 = sqrt(max(sum(w .* sum(abs(exactGrad).^2, 2)), 0));
        state.absH1Semi(idx) = absH1;
        state.relH1Semi(idx) = absH1 / max(exactH1, 1.0e-14);
    end
end
if isfield(problem, 'exactLaplacian')
    exactLap = problem.exactLaplacian(t, material);
    if all(isfinite(exactLap))
        numericalLap = L * c(:);
        lapErr = numericalLap - exactLap(:);
        absH2 = sqrt(max(sum(w .* abs(lapErr).^2), 0));
        exactH2 = sqrt(max(sum(w .* abs(exactLap(:)).^2), 0));
        state.absH2Semi(idx) = absH2;
        state.relH2Semi(idx) = absH2 / max(exactH2, 1.0e-14);
    end
end
end

function diagnostics = finalizeTimeErrorDiagnostics(diagnostics, state, cFinal, cExactFinal, geomFinal)
if all(isfinite(cExactFinal))
    w = surfaceWeights(geomFinal);
    err = cFinal(:) - cExactFinal(:);
    diagnostics.finalAbsL2Error = sqrt(max(sum(w .* abs(err).^2), 0));
    exactL2 = sqrt(max(sum(w .* abs(cExactFinal(:)).^2), 0));
    diagnostics.finalRelL2Error = diagnostics.finalAbsL2Error / max(exactL2, 1.0e-14);
    diagnostics.finalAbsLInfError = max(abs(err));
    diagnostics.finalRelLInfError = diagnostics.finalAbsLInfError / max(max(abs(cExactFinal(:))), 1.0e-14);
end
if ~state.enabled || state.count == 0
    return;
end
t = state.times(1:state.count);
absL2 = state.absL2(1:state.count);
relL2 = state.relL2(1:state.count);
valid = isfinite(t) & isfinite(absL2) & isfinite(relL2);
if ~any(valid)
    return;
end
t = t(valid);
absL2 = absL2(valid);
relL2 = relL2(valid);
diagnostics.maxAbsL2Error = max(absL2);
diagnostics.maxRelL2Error = max(relL2);
absH1 = state.absH1Semi(1:state.count);
relH1 = state.relH1Semi(1:state.count);
validH1 = isfinite(absH1) & isfinite(relH1);
if any(validH1)
    diagnostics.maxAbsH1SemiError = max(absH1(validH1));
    diagnostics.maxRelH1SemiError = max(relH1(validH1));
end
absH2 = state.absH2Semi(1:state.count);
relH2 = state.relH2Semi(1:state.count);
validH2 = isfinite(absH2) & isfinite(relH2);
if any(validH2)
    diagnostics.maxAbsH2SemiError = max(absH2(validH2));
    diagnostics.maxRelH2SemiError = max(relH2(validH2));
end
if numel(t) >= 2
    diagnostics.timeAbsL2L2Error = sqrt(max(trapz(t, absL2.^2), 0));
    diagnostics.timeRelL2L2Error = sqrt(max(trapz(t, relL2.^2), 0));
end
end

function diagnostics = accumulateNormalDiagnostics(diagnostics, normalStats)
if isfinite(normalStats.rmsAngleDegrees)
    if ~isfinite(diagnostics.normalRMSAngle)
        diagnostics.normalRMSAngle = normalStats.rmsAngleDegrees;
    else
        diagnostics.normalRMSAngle = max(diagnostics.normalRMSAngle, ...
            normalStats.rmsAngleDegrees);
    end
end
if isfinite(normalStats.maxAngleDegrees)
    if ~isfinite(diagnostics.normalMaxAngle)
        diagnostics.normalMaxAngle = normalStats.maxAngleDegrees;
    else
        diagnostics.normalMaxAngle = max(diagnostics.normalMaxAngle, ...
            normalStats.maxAngleDegrees);
    end
end
if isfinite(normalStats.controlPointCount)
    if ~isfinite(diagnostics.normalControlPointCount)
        diagnostics.normalControlPointCount = normalStats.controlPointCount;
    else
        diagnostics.normalControlPointCount = max(diagnostics.normalControlPointCount, ...
            normalStats.controlPointCount);
    end
end
if isfinite(normalStats.controlPointFraction)
    if ~isfinite(diagnostics.normalControlPointFraction)
        diagnostics.normalControlPointFraction = normalStats.controlPointFraction;
    else
        diagnostics.normalControlPointFraction = max(diagnostics.normalControlPointFraction, ...
            normalStats.controlPointFraction);
    end
end
if isfinite(normalStats.controlFillDistance)
    if ~isfinite(diagnostics.normalControlFillDistance)
        diagnostics.normalControlFillDistance = normalStats.controlFillDistance;
    else
        diagnostics.normalControlFillDistance = max(diagnostics.normalControlFillDistance, ...
            normalStats.controlFillDistance);
    end
end
end

function state = initMassBalanceState(problem, geom0, material0, cInitial, mu, ...
    nsteps, fixedMassProjection, fixedMassTarget, requestedMode, maxRelativeCorrection)
mode = lower(string(requestedMode));
if fixedMassProjection && (mode == "off" || mode == "none")
    mode = "constant";
end
if mode == "fixed"
    mode = "constant";
elseif mode == "source" || mode == "conservative"
    mode = "balance";
end
if ~(mode == "off" || mode == "none" || mode == "constant" || mode == "balance")
    error('kp:manifold:BadMassCorrectionMode', ...
        'Unknown mass correction mode "%s".', requestedMode);
end

enabled = mode == "constant" || mode == "balance";
initialMass = sum(surfaceWeights(geom0, enabled) .* cInitial(:));
if isfinite(fixedMassTarget)
    constantTarget = fixedMassTarget;
else
    constantTarget = initialMass;
end
sourceMass = nan(nsteps + 1, 1);
sourceIntegral = nan(nsteps + 1, 1);
balanceTarget = nan(nsteps + 1, 1);
correctionTarget = nan(nsteps + 1, 1);
massModel = massBalanceModel(problem);
if mode == "balance" && massModel == "none"
    error('kp:manifold:MissingMassBalanceModel', ...
        ['MassCorrectionMode="balance" requires problem.exactMass, ' ...
        'problem.massRate, problem.massSource, or problem.balanceSource. ' ...
        'Use MassCorrectionMode="constant" for a source-free fixed-mass correction.']);
end
exactMass0 = NaN;
if massModel == "exactMass"
    exactMass0 = massBalanceExactMass(problem, geom0, material0, 0.0, mu);
    sourceIntegral(1) = 0;
    balanceTarget(1) = initialMass;
elseif massModel == "none"
    sourceIntegral(1) = NaN;
    balanceTarget(1) = NaN;
else
    sourceMass(1) = massBalanceRate(problem, massModel, geom0, material0, 0.0, mu, enabled);
    sourceIntegral(1) = 0;
    balanceTarget(1) = initialMass;
end
if mode == "constant"
    correctionTarget(1) = constantTarget;
else
    correctionTarget(1) = balanceTarget(1);
end

state = struct( ...
    'mode', mode, ...
    'enabled', enabled, ...
    'requireExplicitWeights', enabled, ...
    'initialMass', initialMass, ...
    'constantTarget', constantTarget, ...
    'massModel', massModel, ...
    'exactMass0', exactMass0, ...
    'maxRelativeCorrectionAllowed', maxRelativeCorrection, ...
    'sourceMass', sourceMass, ...
    'sourceIntegral', sourceIntegral, ...
    'balanceTarget', balanceTarget, ...
    'correctionTarget', correctionTarget, ...
    'correctionSteps', 0, ...
    'maxAbsCorrection', 0, ...
    'maxRelativeCorrection', 0, ...
    'totalAbsCorrection', 0, ...
    'finalAbsCorrection', 0, ...
    'maxPointShift', 0);
end

function state = advanceMassBalanceState(state, problem, geom, material, mu, dt, step)
timeIndex = step + 1;
t = step * dt;
if state.massModel == "exactMass"
    exactMass = massBalanceExactMass(problem, geom, material, t, mu);
    state.sourceIntegral(timeIndex) = exactMass - state.exactMass0;
    state.balanceTarget(timeIndex) = state.initialMass + state.sourceIntegral(timeIndex);
elseif state.massModel == "none"
    state.sourceMass(timeIndex) = NaN;
    state.sourceIntegral(timeIndex) = NaN;
    state.balanceTarget(timeIndex) = NaN;
else
    state.sourceMass(timeIndex) = massBalanceRate(problem, state.massModel, ...
        geom, material, t, mu, state.requireExplicitWeights);
    if step == 1
        increment = 0.5 * dt * (state.sourceMass(1) + state.sourceMass(2));
        state.sourceIntegral(timeIndex) = increment;
    elseif step == 2
        increment = (dt / 3) * (state.sourceMass(1) + ...
            4 * state.sourceMass(2) + state.sourceMass(3));
        state.sourceIntegral(timeIndex) = increment;
    else
        increment = (dt / 12) * (5 * state.sourceMass(timeIndex) + ...
            8 * state.sourceMass(timeIndex - 1) - state.sourceMass(timeIndex - 2));
        state.sourceIntegral(timeIndex) = state.sourceIntegral(timeIndex - 1) + increment;
    end
    state.balanceTarget(timeIndex) = state.initialMass + state.sourceIntegral(timeIndex);
end
if state.mode == "constant"
    state.correctionTarget(timeIndex) = state.constantTarget;
else
    state.correctionTarget(timeIndex) = state.balanceTarget(timeIndex);
end
end

function kind = massBalanceModel(problem)
if isfield(problem, 'exactMass') && ~isempty(problem.exactMass)
    kind = "exactMass";
elseif isfield(problem, 'massRate') && ~isempty(problem.massRate)
    kind = "massRate";
elseif isfield(problem, 'massSource') && ~isempty(problem.massSource)
    kind = "massSource";
elseif isfield(problem, 'balanceSource') && ~isempty(problem.balanceSource)
    kind = "balanceSource";
else
    kind = "none";
end
end

function val = massBalanceExactMass(problem, geom, material, t, mu)
val = callMassScalarCallback(problem.exactMass, t, geom, material, mu);
end

function val = massBalanceRate(problem, kind, geom, material, t, mu, requireExplicitWeights)
if nargin < 7
    requireExplicitWeights = false;
end
if kind == "massRate"
    val = callMassScalarCallback(problem.massRate, t, geom, material, mu);
    return;
end
w = surfaceWeights(geom, requireExplicitWeights);
if kind == "massSource"
    src = callMassPointCallback(problem.massSource, t, geom, material, mu);
elseif kind == "balanceSource"
    src = callMassPointCallback(problem.balanceSource, t, geom, material, mu);
else
    val = NaN;
    return;
end
val = sum(w .* src(:));
end

function val = callMassScalarCallback(callback, t, geom, material, mu)
nargs = nargin(callback);
if nargs < 0 || nargs >= 4
    val = callback(t, material, geom, mu);
elseif nargs == 3
    val = callback(t, material, mu);
elseif nargs == 2
    val = callback(t, material);
elseif nargs == 1
    val = callback(t);
else
    val = callback();
end
if ~isscalar(val) || ~isfinite(val)
    error('kp:manifold:BadMassScalar', ...
        'Mass scalar callbacks must return a finite scalar.');
end
end

function src = callMassPointCallback(callback, t, geom, material, mu)
nargs = nargin(callback);
if nargs < 0 || nargs >= 4
    src = callback(t, material, geom, mu);
elseif nargs == 3
    src = callback(t, material, mu);
elseif nargs == 2
    src = callback(t, material);
elseif nargs == 1
    src = callback(t);
else
    src = callback();
end
end

function [c, state] = applyMassCorrectionAtStep(c, geom, state, step, countStats)
if ~state.enabled
    return;
end
timeIndex = step + 1;
if timeIndex < 1 || timeIndex > numel(state.correctionTarget)
    return;
end
targetMass = state.correctionTarget(timeIndex);
if ~isfinite(targetMass)
    return;
end
w = surfaceWeights(geom, state.requireExplicitWeights);
denom = sum(w);
if abs(denom) < 1.0e-14
    return;
end
currentMass = sum(w .* c(:));
massDelta = targetMass - currentMass;
relativeCorrection = abs(massDelta) / max([abs(targetMass), abs(currentMass), 1.0e-14]);
if relativeCorrection > state.maxRelativeCorrectionAllowed
    error('kp:manifold:MassCorrectionTooLarge', ...
        ['Mass correction relative size %.3e exceeds the configured ' ...
        'maximum %.3e at step %d.'], ...
        relativeCorrection, state.maxRelativeCorrectionAllowed, step);
end
pointShift = massDelta / denom;
c = c + pointShift;
if countStats
    state.correctionSteps = state.correctionSteps + double(abs(massDelta) > 0);
    state.maxAbsCorrection = max(state.maxAbsCorrection, abs(massDelta));
    state.maxRelativeCorrection = max(state.maxRelativeCorrection, relativeCorrection);
    state.totalAbsCorrection = state.totalAbsCorrection + abs(massDelta);
    state.finalAbsCorrection = abs(massDelta);
    state.maxPointShift = max(state.maxPointShift, abs(pointShift));
end
end

function diagnostics = computeMassDiagnostics(diagnostics, problem, material, mu, ...
    geom0, geomFinal, cInitial, cFinal, cExactFinal, ~, dt, nsteps, massState)
N = size(geom0.X, 1);
requireExplicitWeights = nargin >= 13 && isstruct(massState) && ...
    isfield(massState, 'requireExplicitWeights') && massState.requireExplicitWeights;
w0 = surfaceWeights(geom0, requireExplicitWeights);
wT = surfaceWeights(geomFinal, requireExplicitWeights);
diagnostics.initialMass = sum(w0 .* cInitial);
diagnostics.finalMass = sum(wT .* cFinal);
diagnostics.exactFinalMass = sum(wT .* cExactFinal);
if nargin >= 13 && isstruct(massState) && ...
        numel(massState.sourceIntegral) >= nsteps + 1 && ...
        isfinite(massState.sourceIntegral(nsteps + 1))
    diagnostics.sourceIntegral = massState.sourceIntegral(nsteps + 1);
    diagnostics.massCorrectionSteps = massState.correctionSteps;
    diagnostics.massCorrectionMaxAbs = massState.maxAbsCorrection;
    diagnostics.massCorrectionMaxRelative = massState.maxRelativeCorrection;
    diagnostics.massCorrectionTotalAbs = massState.totalAbsCorrection;
    diagnostics.massCorrectionFinalAbs = massState.finalAbsCorrection;
    diagnostics.massCorrectionMaxPointShift = massState.maxPointShift;
else
    diagnostics.sourceIntegral = integrateBalanceSource(problem, material, mu, dt, nsteps, N);
end
if all(isfinite(cExactFinal))
    diagnostics.massRelError = abs(diagnostics.finalMass - diagnostics.exactFinalMass) / ...
        max(abs(diagnostics.exactFinalMass), 1.0e-14);
else
    diagnostics.massRelError = NaN;
end
balanceResidual = (diagnostics.finalMass - diagnostics.initialMass) - diagnostics.sourceIntegral;
balanceScale = max(abs(diagnostics.initialMass) + abs(diagnostics.sourceIntegral), 1.0e-14);
diagnostics.balanceRelResidual = abs(balanceResidual) / balanceScale;
end

function sourceIntegral = integrateBalanceSource(problem, material, mu, dt, nsteps, N)
kind = massBalanceModel(problem);
if kind == "none"
    sourceIntegral = NaN;
    return;
elseif kind == "exactMass"
    geom0 = problem.geometry(N, 0.0);
    material0 = materialFromGeometry(geom0, material);
    exactMass0 = massBalanceExactMass(problem, geom0, material0, 0.0, mu);
    tFinal = nsteps * dt;
    geomFinal = problem.geometry(N, tFinal);
    materialFinal = materialFromGeometry(geomFinal, material);
    exactMassFinal = massBalanceExactMass( ...
        problem, geomFinal, materialFinal, tFinal, mu);
    sourceIntegral = exactMassFinal - exactMass0;
    return;
end
vals = zeros(nsteps + 1, 1);
for step = 0:nsteps
    t = step * dt;
    geom = problem.geometry(N, t);
    stepMaterial = materialFromGeometry(geom, material);
    w = surfaceWeights(geom);
    if kind == "massRate"
        vals(step + 1) = massBalanceRate(problem, kind, geom, stepMaterial, t, mu);
    else
        src = massBalancePointSource(problem, kind, geom, stepMaterial, t, mu);
        vals(step + 1) = sum(w .* src(:));
    end
end
sourceIntegral = dt * (0.5 * vals(1) + sum(vals(2:end-1)) + 0.5 * vals(end));
end

function material = materialFromGeometry(geom, fallback)
if isfield(geom, 'material') && ~isempty(geom.material)
    material = geom.material;
else
    material = fallback;
end
end

function src = massBalancePointSource(problem, kind, geom, material, t, mu)
if kind == "massSource"
    src = callMassPointCallback(problem.massSource, t, geom, material, mu);
elseif kind == "balanceSource"
    src = callMassPointCallback(problem.balanceSource, t, geom, material, mu);
else
    src = NaN(size(geom.X, 1), 1);
end
end

function w = surfaceWeights(geom, requireExplicitWeights)
if nargin < 2
    requireExplicitWeights = false;
end
N = size(geom.X, 1);
if isfield(geom, 'weights') && ~isempty(geom.weights)
    w = geom.weights(:);
    validateSurfaceWeights(w, N);
elseif requireExplicitWeights
    error('kp:manifold:MissingQuadratureWeights', ...
        ['Mass correction requires explicit geom.weights quadrature values. ', ...
        'Add positive, finite weights to the geometry callback instead of ', ...
        'using area, spacing, or uniform fallback weights.']);
elseif isfield(geom, 'area') && ~isempty(geom.area)
    w = repmat(geom.area / N, N, 1);
elseif isfield(geom, 'h') && ~isempty(geom.h)
    w = repmat(geom.h^2, N, 1);
else
    w = ones(N, 1) / N;
end
end

function validateSurfaceWeights(w, N)
if numel(w) ~= N || any(~isfinite(w)) || any(w <= 0)
    error('kp:manifold:BadQuadratureWeights', ...
        'geom.weights must contain one positive finite quadrature weight per node.');
end
end

function divv = surfaceDivergence(Gx, Gy, Gz, v)
divv = Gx * v(:, 1) + Gy * v(:, 2) + Gz * v(:, 3);
end

function [x, state, stats] = solveGlobalBDFSystem(A, b, xGuess, bdfGamma, state, stats, opts)
solveTimer = tic;
method = lower(string(opts.method));
if method == "direct"
    directTimer = tic;
    x = A \ b;
    stats.directSolveTime = stats.directSolveTime + toc(directTimer);
    stats.directSolves = stats.directSolves + 1;
    stats.maxRelativeResidual = max(stats.maxRelativeResidual, relativeResidual(A, x, b));
    stats.solveTime = stats.solveTime + toc(solveTimer);
    return;
elseif method ~= "gmresiludefect"
    error('kp:manifold:BadGlobalSolveMethod', ...
        'Unknown global solve method "%s".', opts.method);
end

if isempty(xGuess) || numel(xGuess) ~= numel(b)
    xGuess = zeros(size(b));
end

[state, stats] = ensureGlobalILU(A, bdfGamma, state, stats, opts, "auto");
[x0, defectInfo, stats] = applyFrozenILUDefectSweeps( ...
    A, b, xGuess(:), state, stats, opts);
if defectInfo.converged
    x = x0;
    if defectInfo.sweeps == 0
        stats.initialGuessAcceptedSolves = stats.initialGuessAcceptedSolves + 1;
    else
        stats.defectAcceptedSolves = stats.defectAcceptedSolves + 1;
    end
    stats.maxRelativeResidual = max( ...
        stats.maxRelativeResidual, defectInfo.relativeResidual);
    state.solvesSinceRefresh = state.solvesSinceRefresh + 1;
    state.forceRefreshNext = false;
    stats.solveTime = stats.solveTime + toc(solveTimer);
    return;
end

[x, info, stats] = runFrozenILUGMRES(A, b, x0, state, opts, stats);
stats.gmresInvocations = stats.gmresInvocations + 1;

if info.flag ~= 0
    [state, stats] = ensureGlobalILU(A, bdfGamma, state, stats, opts, "failure");
    [x0, defectInfo, stats] = applyFrozenILUDefectSweeps( ...
        A, b, xGuess(:), state, stats, opts);
    if defectInfo.converged
        x = x0;
        stats.defectAcceptedSolves = stats.defectAcceptedSolves + 1;
        info = struct( ...
            'flag', 0, ...
            'reportedRelativeResidual', defectInfo.relativeResidual, ...
            'trueRelativeResidual', defectInfo.relativeResidual, ...
            'iterations', 0, ...
            'residualHistoryLength', 1);
    else
        [x, info, stats] = runFrozenILUGMRES(A, b, x0, state, opts, stats);
        stats.gmresInvocations = stats.gmresInvocations + 1;
    end
end

stats.gmresSolves = stats.gmresSolves + 1;
stats.gmresIterations = stats.gmresIterations + info.iterations;
stats.maxGMRESIterations = max(stats.maxGMRESIterations, info.iterations);
stats.maxRelativeResidual = max(stats.maxRelativeResidual, info.trueRelativeResidual);
stats.maxGMRESFlag = max(stats.maxGMRESFlag, info.flag);
if info.flag ~= 0
    stats.gmresFailures = stats.gmresFailures + 1;
end

state.solvesSinceRefresh = state.solvesSinceRefresh + 1;
state.forceRefreshNext = info.flag ~= 0 || ...
    info.iterations >= opts.iluRefreshIterationThreshold;
stats.solveTime = stats.solveTime + toc(solveTimer);
end

function tol = selectGlobalGMRESTolerance(userTol, h, xi)
if isfinite(userTol) && userTol > 0
    tol = userTol;
    return;
end
tol = max(1.0e-10, min(1.0e-8, 0.05 * h^xi));
end

function [state, stats] = ensureGlobalILU(A, bdfGamma, state, stats, opts, reason)
if isempty(state)
    needsRefresh = true;
    refreshReason = "startup";
elseif state.n ~= size(A, 1)
    needsRefresh = true;
    refreshReason = "size";
elseif abs(state.bdfGamma - bdfGamma) > 10 * eps(max(1, abs(bdfGamma)))
    needsRefresh = true;
    refreshReason = "coefficient";
elseif strcmp(reason, "failure")
    needsRefresh = true;
    refreshReason = "failure";
elseif state.forceRefreshNext
    needsRefresh = true;
    refreshReason = "iteration";
elseif isfinite(opts.iluRefreshInterval) && ...
        state.solvesSinceRefresh >= opts.iluRefreshInterval
    needsRefresh = true;
    refreshReason = "interval";
else
    needsRefresh = false;
    refreshReason = "";
end

if needsRefresh
    refreshTimer = tic;
    state = refreshGlobalILU(A, bdfGamma, opts);
    stats.iluRefreshTime = stats.iluRefreshTime + toc(refreshTimer);
    stats.iluRefreshes = stats.iluRefreshes + 1;
    switch refreshReason
        case "startup"
            stats.iluRefreshesByStartup = stats.iluRefreshesByStartup + 1;
        case "coefficient"
            stats.iluRefreshesByCoefficient = stats.iluRefreshesByCoefficient + 1;
        case "failure"
            stats.iluRefreshesByFailure = stats.iluRefreshesByFailure + 1;
        case "iteration"
            stats.iluRefreshesByIteration = stats.iluRefreshesByIteration + 1;
        case "interval"
            stats.iluRefreshesByInterval = stats.iluRefreshesByInterval + 1;
        otherwise
            stats.iluRefreshesByOther = stats.iluRefreshesByOther + 1;
    end
else
    stats.iluReuses = stats.iluReuses + 1;
end
end

function state = refreshGlobalILU(A, bdfGamma, opts)
setup = struct();
setup.type = 'ilutp';
setup.droptol = opts.iluDropTolerance;
[L, U, P] = ilu(sparse(A), setup);
state = struct( ...
    'n', size(A, 1), ...
    'bdfGamma', bdfGamma, ...
    'L', L, ...
    'U', U, ...
    'P', P, ...
    'solvesSinceRefresh', 0, ...
    'forceRefreshNext', false);
end

function [x, info, stats] = applyFrozenILUDefectSweeps(A, b, x, state, stats, opts)
defectTimer = tic;
residual = b - A * x;
initialResidual = residualNorm(residual, b);
stats.maxInitialRelativeResidual = max(stats.maxInitialRelativeResidual, initialResidual);
numSweeps = 0;
for sweep = 1:max(0, round(opts.defectSweeps))
    relativeResidualNow = residualNorm(residual, b);
    if relativeResidualNow <= opts.gmresTolerance
        break;
    end
    x = x + applyFrozenILU(state, residual);
    residual = b - A * x;
    numSweeps = numSweeps + 1;
    stats.globalDefectSweeps = stats.globalDefectSweeps + 1;
end
postDefectResidual = residualNorm(residual, b);
stats.maxPostDefectRelativeResidual = max( ...
    stats.maxPostDefectRelativeResidual, postDefectResidual);
stats.defectSweepTime = stats.defectSweepTime + toc(defectTimer);
info = struct( ...
    'initialRelativeResidual', initialResidual, ...
    'relativeResidual', postDefectResidual, ...
    'sweeps', numSweeps, ...
    'converged', postDefectResidual <= opts.gmresTolerance);
end

function [x, info, stats] = runFrozenILUGMRES(A, b, x0, state, opts, stats)
restart = max(1, min(size(A, 1), round(opts.gmresRestart)));
maxit = max(1, round(opts.gmresMaxIterations));
tol = max(opts.gmresTolerance, eps);
preconditioner = @(r) applyFrozenILU(state, r);
gmresTimer = tic;
[x, flag, relres, iter, resvec] = gmres(A, b, restart, tol, maxit, preconditioner, [], x0);
stats.gmresTime = stats.gmresTime + toc(gmresTimer);
info = struct( ...
    'flag', flag, ...
    'reportedRelativeResidual', relres, ...
    'trueRelativeResidual', relativeResidual(A, x, b), ...
    'iterations', gmresIterationCount(iter, restart), ...
    'residualHistoryLength', numel(resvec));
end

function y = applyFrozenILU(state, x)
y = state.U \ (state.L \ (state.P * x));
end

function rr = relativeResidual(A, x, b)
rr = residualNorm(b - A * x, b);
end

function rr = residualNorm(r, b)
rr = norm(r) / max(norm(b), 1.0e-14);
end

function niter = gmresIterationCount(iter, restart)
if isempty(iter)
    niter = 0;
elseif numel(iter) == 2
    niter = max(0, (iter(1) - 1) * restart + iter(2));
else
    niter = iter(1);
end
end

function [L, Gx, Gy, Gz, hypGamma, hypPower, hyperviscosityState, stats] = surfaceOpsAndHyperviscosity( ...
    geom, op, requestedHypPower, targetOrder, updater, spectrumCheck, hyperviscosityUpdateMode, ...
    hyperviscosityState, hyperviscosityDriftTolerance, maxHyperviscositySkippedSteps)
neighborCoordinates = surfaceNeighborCoordinates(geom);
if isempty(updater)
    tree = KDTreeSearcher(neighborCoordinates);
    neighborIds = knnsearch(tree, neighborCoordinates, 'K', op.stencilSize);
    [L, Gx, Gy, Gz] = kp.manifold.FormSurfaceDiffOpsTP( ...
        geom.X, op.rbf, op.drbfor, op.d2rbf, geom.normals, ...
        tree, op.stencilSize, op.ell, neighborIds);
    stats = defaultUpdateStats();
    stats.direct = size(geom.X, 1);
    stats.numStencils = size(geom.X, 1);
else
    [L, Gx, Gy, Gz, stats] = updater.assemble( ...
        geom.X, geom.normals, neighborCoordinates);
end
stats.hyperviscosityFullUpdates = 0;
stats.hyperviscosityPredictedUpdates = 0;
stats.maxHyperviscosityRelativeDrift = 0;
stats.minHyperviscosityPower = 0;
stats.maxHyperviscosityPower = 0;

if spectrumCheck
    assertNegativeLargestReal(L);
end

mode = lower(string(hyperviscosityUpdateMode));
if mode == "none" || mode == "off"
    hypGamma = zeros(1, 3);
    hypPower = 1;
elseif mode == "everystep"
    [hypGamma, hyperviscosityState] = refreshHyperviscosityState( ...
        L, Gx, Gy, Gz, geom, requestedHypPower, targetOrder, ...
        acceptedHyperviscosityPower(hyperviscosityState), 0);
    stats.hyperviscosityFullUpdates = 1;
elseif mode == "once"
    if isempty(hyperviscosityState)
        [hypGamma, hyperviscosityState] = refreshHyperviscosityState( ...
            L, Gx, Gy, Gz, geom, requestedHypPower, targetOrder, 1, 0);
        stats.hyperviscosityFullUpdates = 1;
    else
        hypGamma = hyperviscosityState.gamma;
        stats.hyperviscosityPredictedUpdates = 1;
    end
elseif mode == "adaptivegeometry"
    [hypGamma, hyperviscosityState, fullUpdate, relativeDrift] = adaptiveGeometryHyperviscosity( ...
        L, Gx, Gy, Gz, geom, requestedHypPower, targetOrder, hyperviscosityState, ...
        hyperviscosityDriftTolerance, maxHyperviscositySkippedSteps);
    stats.hyperviscosityFullUpdates = double(fullUpdate);
    stats.hyperviscosityPredictedUpdates = double(~fullUpdate);
    stats.maxHyperviscosityRelativeDrift = relativeDrift;
else
    error('kp:manifold:BadHyperviscosityUpdateMode', ...
        'Unknown hyperviscosity update mode "%s".', hyperviscosityUpdateMode);
end
if mode ~= "none" && mode ~= "off"
    hypPower = hyperviscosityState.power;
    stats.minHyperviscosityPower = hypPower;
    stats.maxHyperviscosityPower = hypPower;
end
end

function coordinates = surfaceNeighborCoordinates(geom)
if isfield(geom, 'neighborCoordinates') && ~isempty(geom.neighborCoordinates)
    coordinates = geom.neighborCoordinates;
else
    coordinates = geom.X;
end
end

function [gamma, state, fullUpdate, relativeDrift] = adaptiveGeometryHyperviscosity( ...
    L, Gx, Gy, Gz, geom, requestedHypPower, targetOrder, state, ...
    driftTolerance, maxSkippedSteps)
fullUpdate = true;
relativeDrift = 0;
if isempty(state)
    [gamma, state] = refreshHyperviscosityState( ...
        L, Gx, Gy, Gz, geom, requestedHypPower, targetOrder, 1, 0);
    return;
end

curvatureComponents = surfaceCurvatureComponents(L, geom.X);
if any(~isfinite(curvatureComponents)) || any(curvatureComponents <= 0) || ...
        any(~isfinite(state.curvatureComponents)) || any(state.curvatureComponents <= 0)
    [gamma, state] = refreshHyperviscosityState( ...
        L, Gx, Gy, Gz, geom, requestedHypPower, targetOrder, state.power, 0);
    return;
end

componentScale = state.curvatureComponents ./ curvatureComponents;
componentExponents = 2 * state.power - state.q - 1;
predictedRawComponents = state.rawComponents .* componentScale.^componentExponents;
[eta, ~] = kp.manifold.GetGrowthFunction(L, geom.X, geom.h, state.power);
etaMean = mean(eta);
prefactor = 3^(-state.power) * (-1)^(1 - state.power);
predictedGamma = prefactor * predictedRawComponents / etaMean;
relativeDrift = norm(predictedGamma - state.gamma) / ...
    max(norm(state.gamma), eps);
mustRefresh = relativeDrift > driftTolerance || state.skippedSteps >= maxSkippedSteps;
if mustRefresh
    [gamma, state] = refreshHyperviscosityState( ...
        L, Gx, Gy, Gz, geom, requestedHypPower, targetOrder, state.power, 0);
else
    fullUpdate = false;
    gamma = predictedGamma;
    state.gamma = gamma;
    state.rawComponents = predictedRawComponents;
    state.etaMean = etaMean;
    state.curvatureComponents = curvatureComponents;
    state.skippedSteps = state.skippedSteps + 1;
end
end

function [gamma, state] = refreshHyperviscosityState( ...
    L, Gx, Gy, Gz, geom, requestedHypPower, targetOrder, minimumPower, skippedSteps)
[~, info] = kp.manifold.hyperviscosityCoefficient( ...
    L, Gx, Gy, Gz, geom.X, geom.normals, requestedHypPower, geom.h, ...
    targetOrder, minimumPower);
gamma = info.gammaComponents;
state = struct( ...
    'gamma', gamma, ...
    'power', info.k, ...
    'info', info, ...
    'q', info.q, ...
    'rawComponents', info.hyptermComponents, ...
    'etaMean', info.etaMean, ...
    'curvatureComponents', surfaceCurvatureComponents(L, geom.X), ...
    'skippedSteps', skippedSteps);
end

function correction = hyperviscosityDivergenceCorrection(L, velocity, gamma, power)
poweredVelocity = kp.manifold.applyMatrixPower(L, velocity, power);
if isscalar(gamma)
    gamma = repmat(gamma, 1, size(velocity, 2));
end
correction = poweredVelocity * gamma(:);
end

function power = acceptedHyperviscosityPower(state)
if isempty(state)
    power = 1;
else
    power = state.power;
end
end

function components = surfaceCurvatureComponents(L, X)
meanCurvatureCoordinates = L * X;
components = sqrt(mean(meanCurvatureCoordinates.^2, 1));
end

function stats = defaultUpdateStats()
stats = struct( ...
    'numStencils', 0, ...
    'direct', 0, ...
    'defectCorrected', 0, ...
    'defectFailedRefactored', 0, ...
    'meanDefectIterations', NaN, ...
    'maxDefectIterations', 0, ...
    'maxRelativeResidual', 0, ...
    'neighborSearches', 0, ...
    'neighborReuses', 0, ...
    'neighborReuseRejected', 0, ...
    'hyperviscosityFullUpdates', 0, ...
    'hyperviscosityPredictedUpdates', 0, ...
    'maxHyperviscosityRelativeDrift', 0, ...
    'minHyperviscosityPower', Inf, ...
    'maxHyperviscosityPower', 0);
end

function stats = defaultGlobalSolveStats()
stats = struct( ...
    'directSolves', 0, ...
    'gmresSolves', 0, ...
    'gmresInvocations', 0, ...
    'gmresFailures', 0, ...
    'gmresIterations', 0, ...
    'maxGMRESIterations', 0, ...
    'maxGMRESFlag', 0, ...
    'iluRefreshes', 0, ...
    'iluReuses', 0, ...
    'iluRefreshesByStartup', 0, ...
    'iluRefreshesByCoefficient', 0, ...
    'iluRefreshesByFailure', 0, ...
    'iluRefreshesByIteration', 0, ...
    'iluRefreshesByInterval', 0, ...
    'iluRefreshesByOther', 0, ...
    'globalDefectSweeps', 0, ...
    'initialGuessAcceptedSolves', 0, ...
    'defectAcceptedSolves', 0, ...
    'maxInitialRelativeResidual', 0, ...
    'maxPostDefectRelativeResidual', 0, ...
    'maxRelativeResidual', 0, ...
    'solveTime', 0, ...
    'directSolveTime', 0, ...
    'iluRefreshTime', 0, ...
    'defectSweepTime', 0, ...
    'gmresTime', 0);
end

function stats = defaultAccumulatedUpdateStats()
stats = defaultUpdateStats();
stats.numAssemblies = 0;
stats.totalDefectIterations = 0;
end

function total = accumulateUpdateStats(total, stats)
if isempty(stats)
    return;
end
total.numAssemblies = total.numAssemblies + 1;
total.numStencils = total.numStencils + stats.numStencils;
total.direct = total.direct + stats.direct;
total.defectCorrected = total.defectCorrected + stats.defectCorrected;
total.defectFailedRefactored = total.defectFailedRefactored + stats.defectFailedRefactored;
total.maxDefectIterations = max(total.maxDefectIterations, stats.maxDefectIterations);
total.maxRelativeResidual = max(total.maxRelativeResidual, stats.maxRelativeResidual);
total.neighborSearches = total.neighborSearches + stats.neighborSearches;
total.neighborReuses = total.neighborReuses + stats.neighborReuses;
total.neighborReuseRejected = total.neighborReuseRejected + stats.neighborReuseRejected;
total.hyperviscosityFullUpdates = total.hyperviscosityFullUpdates + stats.hyperviscosityFullUpdates;
total.hyperviscosityPredictedUpdates = total.hyperviscosityPredictedUpdates + stats.hyperviscosityPredictedUpdates;
total.maxHyperviscosityRelativeDrift = max(total.maxHyperviscosityRelativeDrift, ...
    stats.maxHyperviscosityRelativeDrift);
total.minHyperviscosityPower = min(total.minHyperviscosityPower, ...
    stats.minHyperviscosityPower);
total.maxHyperviscosityPower = max(total.maxHyperviscosityPower, ...
    stats.maxHyperviscosityPower);
if stats.defectCorrected > 0
    total.totalDefectIterations = total.totalDefectIterations + ...
        stats.meanDefectIterations * stats.defectCorrected;
    total.meanDefectIterations = total.totalDefectIterations / total.defectCorrected;
end
end

function assertNegativeLargestReal(L)
tol = 1.0e-2;
warnState = warning();
cleanup = onCleanup(@() warning(warnState));
warning('off', 'all');
[~, sig] = eigs(L, 1, 'largestreal', 'Tolerance', tol);
if sig > tol
    error('kp:manifold:BadSpectrum', ...
        ['+ve real eigenvalue detected, smooth normals with more nearest ' ...
         'neighbors or check points.']);
end
end

function writeConvergenceFigure(results, titleText, imgPath)
fig = figure('Color', 'w', 'Position', [100, 100, 760, 560]);
hold on;
colors = lines(numel(results));
for ix = 1:numel(results)
    plot(results(ix).h, results(ix).relerr, '-o', 'LineWidth', 1.6, ...
        'MarkerSize', 7, 'Color', colors(ix, :), ...
        'DisplayName', sprintf('\\xi = %d', results(ix).xi));
end
hold off;
set(gca, 'XScale', 'log', 'YScale', 'log', 'XDir', 'reverse');
grid on;
xlabel('h');
ylabel('Relative L2 error');
title(titleText);
legend('Location', 'best');
exportgraphics(fig, imgPath, 'Resolution', 180);
end

function validateProblem(problem)
required = {'label', 'title', 'geometry', 'forcing'};
for k = 1:numel(required)
    if ~isfield(problem, required{k})
        error('kp:manifold:InvalidProblem', ...
            'Moving ADR problem is missing required field "%s".', required{k});
    end
end
if ~isfield(problem, 'exact') && ~isfield(problem, 'initial')
    error('kp:manifold:InvalidProblem', ...
        'Moving ADR problem must provide either an exact or initial callback.');
end
end
