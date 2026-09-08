function study = moving_surface_adr_tp_spheroid_rearrangement(varargin)
%MOVING_SURFACE_ADR_TP_SPHEROID_REARRANGEMENT Forced-advection rearrangement test.
%   Solves a manufactured conservative advection problem on a fixed spheroid.
%   The numerical update uses tangent-plane RBF-FD surface divergence only.
%   Exact parameter calculus is used only to evaluate the forcing and the
%   final reference solution.  Rearrangement uses a semi-Lagrangian bootstrap
%   to rebuild the BDF history on a quasi-uniform marker cloud.

parser = inputParser();
parser.addParameter('N', 512, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('Xi', 2, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('FinalTime', 0.14, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('DtScale', 0.08, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('FlowScale', 0.8, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('QualityThreshold', 1.55, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('PredictiveLookaheadSteps', 3, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('MinStepsBetweenRearrangements', 8, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('RearrangementTimes', [], @(x) isnumeric(x));
parser.addParameter('ExactRearrangementBootstrap', false, @(x) islogical(x) && isscalar(x));
parser.addParameter('MaterialSampler', "fps", @(x) isstring(x) || ischar(x));
parser.addParameter('SamplerCandidateFactor', 8, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('SBFControlPointCount', NaN, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('SBFControlPointFraction', 0.5, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('TransferMode', "sbf", @(x) isstring(x) || ischar(x));
parser.addParameter('DivergenceMode', "discrete", @(x) isstring(x) || ischar(x));
parser.addParameter('MassCorrectionMode', "balance", @(x) isstring(x) || ischar(x));
parser.addParameter('MassCorrectionMaxRelativeCorrection', Inf, ...
    @(x) isnumeric(x) && isscalar(x) && x > 0);
parser.addParameter('UseRearrangement', [false, true], @(x) islogical(x) || isnumeric(x));
parser.addParameter('UseParallel', false, @(x) islogical(x) && isscalar(x));
parser.addParameter('OutputPrefix', fullfile(pwd, ...
    'moving_surface_adr_tp_spheroid_rearrangement'));
parser.addParameter('ImagePath', fullfile(pwd, 'docs', 'figures', ...
    'moving_surface_adr_tp_spheroid_rearrangement.png'));
parser.addParameter('SaveOutputs', true, @(x) islogical(x) && isscalar(x));
parser.parse(varargin{:});
opts = parser.Results;

N = round(opts.N);
xi = round(opts.Xi);
theta = 2;
if isnan(opts.SBFControlPointCount)
    sbfControlPointCount = max(1, round(opts.SBFControlPointFraction * N));
else
    sbfControlPointCount = opts.SBFControlPointCount;
end
h = sqrt(spheroidAreaApprox() / N);
dtTarget = opts.DtScale * h^(xi / 3);
nsteps = max(3, ceil(opts.FinalTime / dtTarget));
dt = opts.FinalTime / nsteps;

modes = logical(opts.UseRearrangement(:).');
results = repmat(emptyRunResult(), 1, numel(modes));
for k = 1:numel(modes)
    results(k) = runOneMode(N, xi, theta, opts.FinalTime, dt, nsteps, ...
        modes(k), opts.QualityThreshold, round(opts.PredictiveLookaheadSteps), ...
        round(opts.MinStepsBetweenRearrangements), opts.FlowScale, ...
        sbfControlPointCount, string(opts.TransferMode), string(opts.DivergenceMode), ...
        opts.RearrangementTimes, opts.ExactRearrangementBootstrap, ...
        string(opts.MaterialSampler), opts.SamplerCandidateFactor, ...
        string(opts.MassCorrectionMode), opts.MassCorrectionMaxRelativeCorrection, ...
        opts.UseParallel);
end

study = struct();
study.N = N;
study.xi = xi;
study.h = h;
study.dt = dt;
study.nsteps = nsteps;
study.sbfControlPointCount = sbfControlPointCount;
study.transferMode = string(opts.TransferMode);
study.materialSampler = string(opts.MaterialSampler);
study.samplerCandidateFactor = opts.SamplerCandidateFactor;
study.results = results;
study.table = struct2table(results);

if opts.SaveOutputs
    save(string(opts.OutputPrefix) + ".mat", 'study');
    writetable(study.table, string(opts.OutputPrefix) + ".csv");
end
if ~isempty(opts.ImagePath)
    writeStudyFigure(study, opts.ImagePath);
end
end

function result = runOneMode(N, xi, theta, T, dt, nsteps, useRearrangement, ...
    qualityThreshold, predictiveLookaheadSteps, minStepsBetweenRearrangements, ...
    flowScale, sbfControlPointCount, transferMode, divergenceMode, rearrangementTimesTarget, ...
    exactRearrangementBootstrap, materialSampler, samplerCandidateFactor, ...
    massCorrectionMode, massCorrectionMaxRelativeCorrection, useParallel)
flow = spheroidFlowParameters(flowScale);
updater = kp.manifold.TangentPlaneDiffMatUpdater(xi, ...
    'Theta', theta, ...
    'DefectTolerance', 1.0e-6, ...
    'MaxDefectIterations', 4, ...
    'NeighborUpdateMode', "motionbound", ...
    'UseParallel', useParallel);

[lambda0, eta0] = spheroidMaterialSites(N, materialSampler, samplerCandidateFactor);
t0 = 0;
geom0 = spheroidGeometry(lambda0, eta0);
x0 = geom0.X;
c0 = exactConcentration(t0, lambda0, eta0);
massState = initSpheroidMassCorrection(massCorrectionMode, ...
    massCorrectionMaxRelativeCorrection, geom0, c0, t0, lambda0, eta0);
[c0, massState] = applySpheroidMassCorrection(c0, geom0, t0, ...
    lambda0, eta0, massState, false);

[lambda1, eta1] = advanceLabels(lambda0, eta0, dt, flow);
t1 = dt;
geom1 = spheroidGeometry(lambda1, eta1);
x1 = geom1.X;
c1 = exactConcentration(t1, lambda1, eta1);
[c1, massState] = applySpheroidMassCorrection(c1, geom1, t1, ...
    lambda1, eta1, massState, false);

[lambda2, eta2] = advanceLabels(lambda1, eta1, dt, flow);
t2 = 2 * dt;
geom2 = spheroidGeometry(lambda2, eta2);
x2 = geom2.X;
c2 = exactConcentration(t2, lambda2, eta2);
[c2, massState] = applySpheroidMassCorrection(c2, geom2, t2, ...
    lambda2, eta2, massState, false);

qualityTimes = zeros(nsteps + 1, 1);
qualityValues = zeros(nsteps + 1, 1);
anisotropyValues = zeros(nsteps + 1, 1);
qualityTimes(1:3) = [t0; t1; t2];
qualityValues(1:3) = [pointQuality(x0); pointQuality(x1); pointQuality(x2)];
anisotropyValues(1:3) = [stencilAnisotropy(x0, geom0.normals); ...
    stencilAnisotropy(x1, geom1.normals); stencilAnisotropy(x2, geom2.normals)];
rearrangementTimes = zeros(0, 1);
triggerValues = zeros(0, 1);
triggerQualities = zeros(0, 1);
triggeredByPrediction = false(0, 1);
triggeredByThreshold = false(0, 1);
triggeredByFixedSchedule = false(0, 1);
sbfControlPointCounts = zeros(0, 1);
sbfFillDistances = zeros(0, 1);
transferStencilSizes = zeros(0, 1);
operatorStats = repmat(emptyOperatorStats(), 0, 1);
lastRearrangementStep = -inf;
fixedRearrangementTimes = sort(rearrangementTimesTarget(:));
nextFixedRearrangement = 1;
useFixedSchedule = ~isempty(fixedRearrangementTimes);

for step = 3:nsteps
    t3 = step * dt;
    [lambda3, eta3] = advanceLabels(lambda2, eta2, dt, flow);
    geom3 = spheroidGeometry(lambda3, eta3);
    x3 = geom3.X;
    velocity3 = bdfVelocity(x0, x1, x2, x3, dt);
    [c3, stats3] = solveBDF3Step(updater, geom3, velocity3, t3, ...
        lambda3, eta3, c0, c1, c2, dt, flow, divergenceMode);
    [c3, massState] = applySpheroidMassCorrection(c3, geom3, t3, ...
        lambda3, eta3, massState, true);
    operatorStats = [operatorStats; stats3]; %#ok<AGROW>

    lambda0 = lambda1; eta0 = eta1; x0 = x1; c0 = c1;
    lambda1 = lambda2; eta1 = eta2; x1 = x2; c1 = c2;
    lambda2 = lambda3; eta2 = eta3; x2 = x3; c2 = c3;

    q = pointQuality(x2);
    amax = stencilAnisotropy(x2, geom3.normals);
    qualityTimes(step + 1) = t3;
    qualityValues(step + 1) = q;
    anisotropyValues(step + 1) = amax;

    [willCross, triggerMetric] = predictQualityCrossing( ...
        qualityValues, step + 1, qualityThreshold, predictiveLookaheadSteps);
    fixedTrigger = useFixedSchedule && nextFixedRearrangement <= numel(fixedRearrangementTimes) && ...
        t3 >= fixedRearrangementTimes(nextFixedRearrangement) - 0.5 * dt;
    adaptiveTrigger = ~useFixedSchedule && (q > qualityThreshold || willCross);
    canRearrange = useRearrangement && step >= 3 && ...
        step - lastRearrangementStep >= minStepsBetweenRearrangements && ...
        (fixedTrigger || adaptiveTrigger);
    if canRearrange
        thresholdTrigger = ~useFixedSchedule && q > qualityThreshold;
        predictiveTrigger = ~useFixedSchedule && ~thresholdTrigger && willCross;
        [lambda0, eta0, x0, c0, lambda1, eta1, x1, c1, lambda2, eta2, x2, c2, transferInfo] = ...
            rearrangeHistory(lambda0, eta0, c0, lambda1, eta1, c1, ...
            lambda2, eta2, c2, dt, flow, sbfControlPointCount, transferMode, ...
            xi, theta, t3, exactRearrangementBootstrap, materialSampler, ...
            samplerCandidateFactor);
        geom0 = spheroidGeometry(lambda0, eta0);
        geom1 = spheroidGeometry(lambda1, eta1);
        geom2 = spheroidGeometry(lambda2, eta2);
        [c0, massState] = applySpheroidMassCorrection(c0, geom0, t3 - 2 * dt, ...
            lambda0, eta0, massState, true);
        [c1, massState] = applySpheroidMassCorrection(c1, geom1, t3 - dt, ...
            lambda1, eta1, massState, true);
        [c2, massState] = applySpheroidMassCorrection(c2, geom2, t3, ...
            lambda2, eta2, massState, true);
        updater.reset();
        lastRearrangementStep = step;
        if fixedTrigger
            nextFixedRearrangement = nextFixedRearrangement + 1;
        end
        rearrangementTimes(end + 1, 1) = t3; %#ok<AGROW>
        triggerValues(end + 1, 1) = triggerMetric; %#ok<AGROW>
        triggerQualities(end + 1, 1) = q; %#ok<AGROW>
        triggeredByPrediction(end + 1, 1) = predictiveTrigger; %#ok<AGROW>
        triggeredByThreshold(end + 1, 1) = thresholdTrigger; %#ok<AGROW>
        triggeredByFixedSchedule(end + 1, 1) = fixedTrigger; %#ok<AGROW>
        sbfControlPointCounts(end + 1, 1) = transferInfo.controlPointCount; %#ok<AGROW>
        sbfFillDistances(end + 1, 1) = transferInfo.fillDistance; %#ok<AGROW>
        transferStencilSizes(end + 1, 1) = transferInfo.localStencilSize; %#ok<AGROW>
        qualityValues(step + 1) = pointQuality(x2);
        anisotropyValues(step + 1) = stencilAnisotropy(x2, spheroidNormals(x2));
    end
end

cex = exactConcentration(T, lambda2, eta2);
geomFinal = spheroidGeometry(lambda2, eta2);
[initialMass, finalMass, exactFinalMass, massRelError, balanceRelResidual] = ...
    spheroidMassDiagnostics(geom0, geomFinal, c0, c2, cex, massState);
[predictedCrossTime, actualCrossTime, predictorLeadSteps] = predictorSummary( ...
    qualityTimes, qualityValues, qualityThreshold, predictiveLookaheadSteps);
result = emptyRunResult();
if useRearrangement
    result.mode = "rearranged";
    result.transferMode = string(transferMode);
else
    result.mode = "lagrangian";
    result.transferMode = "none";
end
result.relerr = norm(c2 - cex) / max(norm(cex), 1.0e-14);
result.initialMass = initialMass;
result.finalMass = finalMass;
result.exactFinalMass = exactFinalMass;
result.massRelError = massRelError;
result.balanceRelResidual = balanceRelResidual;
result.massCorrectionSteps = massState.correctionSteps;
result.massCorrectionMaxAbs = massState.maxAbsCorrection;
result.massCorrectionMaxRelative = massState.maxRelativeCorrection;
result.massCorrectionTotalAbs = massState.totalAbsCorrection;
result.massCorrectionFinalAbs = massState.finalAbsCorrection;
result.massCorrectionMaxPointShift = massState.maxPointShift;
result.numRearrangements = numel(rearrangementTimes);
result.numPredictiveRearrangements = nnz(triggeredByPrediction);
result.numThresholdRearrangements = nnz(triggeredByThreshold);
result.numFixedRearrangements = nnz(triggeredByFixedSchedule);
result.numLateRearrangements = result.numThresholdRearrangements;
result.initialQuality = qualityValues(1);
result.finalQuality = qualityValues(nsteps + 1);
result.maxQuality = max(qualityValues);
result.initialAnisotropy = anisotropyValues(1);
result.finalAnisotropy = anisotropyValues(nsteps + 1);
result.maxAnisotropy = max(anisotropyValues);
result.rearrangementTimes = {rearrangementTimes};
result.rearrangementTriggerValues = {triggerValues};
result.rearrangementQualityValues = {triggerQualities};
result.rearrangementWasPredicted = {triggeredByPrediction};
result.rearrangementWasThresholdTriggered = {triggeredByThreshold};
result.rearrangementWasFixedScheduled = {triggeredByFixedSchedule};
if isempty(sbfControlPointCounts)
    result.sbfControlPointCount = NaN;
    result.sbfControlPointFraction = NaN;
    result.sbfFillDistance = NaN;
else
    result.sbfControlPointCount = max(sbfControlPointCounts);
    result.sbfControlPointFraction = result.sbfControlPointCount / N;
    result.sbfFillDistance = max(sbfFillDistances);
end
if isempty(transferStencilSizes)
    result.transferStencilSize = NaN;
else
    result.transferStencilSize = max(transferStencilSizes);
end
result.predictedQualityCrossingTime = predictedCrossTime;
result.actualQualityCrossingTime = actualCrossTime;
result.predictorLeadSteps = predictorLeadSteps;
result.qualityTimes = {qualityTimes};
result.qualityValues = {qualityValues};
result.anisotropyValues = {anisotropyValues};
result.finalLambda = {lambda2};
result.finalEta = {eta2};
result.finalX = {x2};
result.finalNumerical = {c2};
result.finalExact = {cex};
result.operatorStats = {operatorStats};
end

function [c3, stats] = solveBDF3Step(updater, geom, velocity, t, lambda, eta, ...
    c0, c1, c2, dt, flow, divergenceMode)
switch lower(string(divergenceMode))
    case "discrete"
        [~, Gx, Gy, Gz, stats] = updater.assemble(geom.X, geom.normals);
        divv = Gx * velocity(:, 1) + Gy * velocity(:, 2) + Gz * velocity(:, 3);
    case "exact"
        stats = emptyOperatorStats();
        divv = exactSurfaceDivergence(lambda, eta, flow);
    otherwise
        error('kp:examples:BadDivergenceMode', ...
            'Unknown divergence mode "%s".', divergenceMode);
end
forcing = manufacturedForcing(t, lambda, eta, flow);
beta = 6 / 11;
extrap = 3 * c2 - 3 * c1 + c0;
c3 = (18 / 11) * c2 - (9 / 11) * c1 + (2 / 11) * c0 + ...
    beta * dt * forcing - beta * dt * extrap .* divv;
end

function velocity = bdfVelocity(x0, x1, x2, x3, dt)
velocity = (11 * x3 - 18 * x2 + 9 * x1 - 2 * x0) / (6 * dt);
end

function [lambda0New, eta0New, x0New, c0New, lambda1New, eta1New, x1New, c1New, ...
    lambda2New, eta2New, x2New, c2New, transferInfo] = rearrangeHistory( ...
    lambda0, eta0, c0, lambda1, eta1, c1, lambda2, eta2, c2, dt, flow, ...
    sbfControlPointCount, transferMode, xi, theta, tCurrent, exactBootstrap, ...
    materialSampler, samplerCandidateFactor)
[lambda2New, eta2New] = spheroidMaterialSites(numel(lambda2), materialSampler, ...
    samplerCandidateFactor);
[lambda1New, eta1New] = advanceLabels(lambda2New, eta2New, -dt, flow);
[lambda0New, eta0New] = advanceLabels(lambda2New, eta2New, -2 * dt, flow);

x2New = spheroidPosition(lambda2New, eta2New);
x1New = spheroidPosition(lambda1New, eta1New);
x0New = spheroidPosition(lambda0New, eta0New);

if exactBootstrap
    transferInfo = struct( ...
        'mode', "exact", ...
        'controlPointCount', NaN, ...
        'controlPointFraction', NaN, ...
        'fillDistance', NaN, ...
        'localStencilSize', NaN);
    c2New = exactConcentration(tCurrent, lambda2New, eta2New);
    c1New = exactConcentration(tCurrent - dt, lambda1New, eta1New);
    c0New = exactConcentration(tCurrent - 2 * dt, lambda0New, eta0New);
    return;
end

switch lower(string(transferMode))
    case "sbf"
        controlIds = scalarSBFControlIds(lambda2, eta2, sbfControlPointCount);
        transferInfo = scalarSBFControlInfo(lambda2, eta2, controlIds);
        transferInfo.mode = "sbf";
        transferInfo.localStencilSize = NaN;
        c2New = sphericalScalarSBFInterpolate(lambda2, eta2, c2, lambda2New, eta2New, controlIds);
        c1New = sphericalScalarSBFInterpolate(lambda1, eta1, c1, lambda1New, eta1New, controlIds);
        c0New = sphericalScalarSBFInterpolate(lambda0, eta0, c0, lambda0New, eta0New, controlIds);
    case {"localtp", "localtangentplane"}
        op = kp.manifold.rbffdop(2, xi, theta, 0);
        transferInfo = struct();
        transferInfo.mode = "localTp";
        transferInfo.controlPointCount = NaN;
        transferInfo.controlPointFraction = NaN;
        transferInfo.fillDistance = NaN;
        transferInfo.localStencilSize = op.stencilSize;
        c2New = localTangentPlaneScalarInterpolate( ...
            spheroidPosition(lambda2, eta2), c2, x2New, xi, theta);
        c1New = localTangentPlaneScalarInterpolate( ...
            spheroidPosition(lambda1, eta1), c1, x1New, xi, theta);
        c0New = localTangentPlaneScalarInterpolate( ...
            spheroidPosition(lambda0, eta0), c0, x0New, xi, theta);
    otherwise
        error('kp:examples:BadTransferMode', ...
            'Unknown transfer mode "%s".', transferMode);
end
end

function c = exactConcentration(t, lambda, eta)
decay = exp(-t);
base = smoothConcentrationBase(lambda, eta);
c = decay .* base;
end

function state = initSpheroidMassCorrection(mode, maxRelativeCorrection, ...
    geom0, c0, t0, lambda0, eta0)
mode = lower(string(mode));
if mode == "off" || mode == "none"
    enabled = false;
elseif mode == "fixed" || mode == "constant"
    mode = "constant";
    enabled = true;
elseif mode == "source" || mode == "conservative" || mode == "exactmass"
    mode = "balance";
    enabled = true;
elseif mode == "balance"
    enabled = true;
else
    error('kp:examples:BadMassCorrectionMode', ...
        'Unknown mass correction mode "%s".', mode);
end

initialMass = sum(geom0.weights(:) .* c0(:));
exactInitialMass = spheroidExactMass(geom0, t0, lambda0, eta0);
state = struct( ...
    'mode', mode, ...
    'enabled', enabled, ...
    'initialMass', initialMass, ...
    'exactInitialMass', exactInitialMass, ...
    'maxRelativeCorrectionAllowed', maxRelativeCorrection, ...
    'correctionSteps', 0, ...
    'maxAbsCorrection', 0, ...
    'maxRelativeCorrection', 0, ...
    'totalAbsCorrection', 0, ...
    'finalAbsCorrection', 0, ...
    'maxPointShift', 0);
end

function [c, state] = applySpheroidMassCorrection(c, geom, t, lambda, eta, ...
    state, countStats)
if ~state.enabled
    return;
end
if state.mode == "constant"
    targetMass = state.initialMass;
else
    targetMass = state.initialMass + ...
        (spheroidExactMass(geom, t, lambda, eta) - state.exactInitialMass);
end
w = geom.weights(:);
denom = sum(w);
currentMass = sum(w .* c(:));
massDelta = targetMass - currentMass;
relativeCorrection = abs(massDelta) / max([abs(targetMass), abs(currentMass), 1.0e-14]);
if relativeCorrection > state.maxRelativeCorrectionAllowed
    error('kp:examples:MassCorrectionTooLarge', ...
        ['Mass correction relative size %.3e exceeds the configured ', ...
        'maximum %.3e.'], relativeCorrection, state.maxRelativeCorrectionAllowed);
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

function exactMass = spheroidExactMass(geom, t, lambda, eta)
exactMass = sum(geom.weights(:) .* exactConcentration(t, lambda, eta));
end

function [initialMass, finalMass, exactFinalMass, massRelError, balanceRelResidual] = ...
    spheroidMassDiagnostics(~, geomFinal, ~, cFinal, cExactFinal, massState)
wT = geomFinal.weights(:);
initialMass = massState.initialMass;
finalMass = sum(wT .* cFinal(:));
exactFinalMass = sum(wT .* cExactFinal(:));
sourceIntegral = exactFinalMass - massState.exactInitialMass;
massRelError = abs(finalMass - exactFinalMass) / max(abs(exactFinalMass), 1.0e-14);
balanceResidual = (finalMass - initialMass) - sourceIntegral;
balanceScale = max(abs(initialMass) + abs(sourceIntegral), 1.0e-14);
balanceRelResidual = abs(balanceResidual) / balanceScale;
end

function f = manufacturedForcing(t, lambda, eta, flow)
decay = exp(-t);
base = smoothConcentrationBase(lambda, eta);
gradBase = smoothConcentrationGradient(lambda, eta);
velocity = analyticVelocity(lambda, eta, flow);
c = decay .* base;
ct = -c;
divv = exactSurfaceDivergence(lambda, eta, flow);
f = ct + decay .* sum(gradBase .* velocity, 2) + c .* divv;
end

function base = smoothConcentrationBase(lambda, eta)
X = spheroidPosition(lambda, eta);
x = X(:, 1);
y = X(:, 2);
z = X(:, 3);
base = 1 + 0.12 * x + 0.08 * y .* z + ...
    0.06 * (x .^ 2 - y .^ 2) + 0.05 * x .* y .* z;
end

function grad = smoothConcentrationGradient(lambda, eta)
X = spheroidPosition(lambda, eta);
x = X(:, 1);
y = X(:, 2);
z = X(:, 3);
grad = [ ...
    0.12 + 0.12 * x + 0.05 * y .* z, ...
    0.08 * z - 0.12 * y + 0.05 * x .* z, ...
    0.08 * y + 0.05 * x .* y];
end

function divv = exactSurfaceDivergence(lambda, eta, flow)
[~, etaDot] = labelVelocity(lambda, eta, flow);
lambdaDotLambda = zeros(size(lambda));
etaDotEta = -2 * flow.B * eta;
p = spheroidParameters();
metric = p.c^2 + (p.a^2 - p.c^2) * eta .^ 2;
logSqrtgEta = (p.a^2 - p.c^2) * eta ./ metric;
divv = lambdaDotLambda + etaDotEta + logSqrtgEta .* etaDot;
end

function [lambdaDot, etaDot] = labelVelocity(~, eta, flow)
lambdaDot = flow.Omega + flow.S * eta;
etaDot = flow.B * (1 - eta .^ 2);
end

function velocity = analyticVelocity(lambda, eta, flow)
[lambdaDot, etaDot] = labelVelocity(lambda, eta, flow);
p = spheroidParameters();
r = sqrt(max(1 - eta(:) .^ 2, 0));
Xlambda = [-p.a * r .* sin(lambda(:)), p.a * r .* cos(lambda(:)), zeros(numel(eta), 1)];
safeR = max(r, sqrt(eps));
Xeta = [-p.a * eta(:) .* cos(lambda(:)) ./ safeR, ...
    -p.a * eta(:) .* sin(lambda(:)) ./ safeR, ...
    p.c * ones(numel(eta), 1)];
velocity = lambdaDot(:) .* Xlambda + etaDot(:) .* Xeta;
end

function [lambda, eta] = advanceLabels(lambda, eta, dt, flow)
maxStep = 2.5e-3;
nsub = max(1, ceil(abs(dt) / maxStep));
subdt = dt / nsub;
for k = 1:nsub
    [k1l, k1e] = labelVelocity(lambda, eta, flow);
    [k2l, k2e] = labelVelocity(lambda + 0.5 * subdt * k1l, ...
        clampEta(eta + 0.5 * subdt * k1e), flow);
    [k3l, k3e] = labelVelocity(lambda + 0.5 * subdt * k2l, ...
        clampEta(eta + 0.5 * subdt * k2e), flow);
    [k4l, k4e] = labelVelocity(lambda + subdt * k3l, ...
        clampEta(eta + subdt * k3e), flow);
    lambda = lambda + (subdt / 6) * (k1l + 2 * k2l + 2 * k3l + k4l);
    eta = eta + (subdt / 6) * (k1e + 2 * k2e + 2 * k3e + k4e);
    eta = clampEta(eta);
end
lambda = mod(lambda, 2 * pi);
end

function eta = clampEta(eta)
eta = max(min(eta, 1 - 1.0e-10), -1 + 1.0e-10);
end

function geom = spheroidGeometry(lambda, eta)
geom = struct();
geom.X = spheroidPosition(lambda, eta);
geom.normals = spheroidNormals(geom.X);
geom.weights = spheroidVertexAreaWeights(geom.X);
geom.area = sum(geom.weights);
geom.h = sqrt(geom.area / numel(lambda));
geom.material = struct('lambda', lambda(:), 'eta', eta(:));
end

function weights = spheroidVertexAreaWeights(X)
faces = convhull(X(:, 1), X(:, 2), X(:, 3));
e1 = X(faces(:, 2), :) - X(faces(:, 1), :);
e2 = X(faces(:, 3), :) - X(faces(:, 1), :);
area = 0.5 * vecnorm(cross(e1, e2, 2), 2, 2);
weights = accumarray(faces(:), repmat(area / 3, 3, 1), ...
    [size(X, 1), 1], @sum, 0);
if numel(weights) ~= size(X, 1) || any(~isfinite(weights)) || any(weights <= 0)
    error('kp:examples:BadSpheroidQuadrature', ...
        'Spheroid vertex-area quadrature failed to produce positive node weights.');
end
end

function X = spheroidPosition(lambda, eta)
p = spheroidParameters();
lambda = lambda(:);
eta = eta(:);
r = sqrt(max(1 - eta .^ 2, 0));
X = [p.a * r .* cos(lambda), p.a * r .* sin(lambda), p.c * eta];
end

function normals = spheroidNormals(X)
p = spheroidParameters();
raw = [X(:, 1) / p.a^2, X(:, 2) / p.a^2, X(:, 3) / p.c^2];
normals = kp.geometry.normalizeRows(raw);
end

function [lambda, eta] = spheroidMaterialSites(N, sampler, candidateFactor)
switch lower(string(sampler))
    case {"fps", "farthest", "farthestpoint", "farthestpointsampling"}
        [lambda, eta] = spheroidFarthestPointMaterial(N, candidateFactor);
    case {"fibonacci", "spiral"}
        [lambda, eta] = spheroidFibonacciMaterial(N);
    otherwise
        error('kp:examples:BadMaterialSampler', ...
            'Unknown material sampler "%s".', sampler);
end
end

function [lambda, eta] = spheroidFibonacciMaterial(N)
i = (1:N).';
eta = 1 - 2 * (i - 0.5) / N;
golden = pi * (3 - sqrt(5));
lambda = mod(golden * (i - 1), 2 * pi);
end

function [lambda, eta] = spheroidFarthestPointMaterial(N, candidateFactor)
candidateFactor = max(1, candidateFactor);
numCandidates = max(N, ceil(candidateFactor * N));
[lambdaCandidates, etaCandidates] = spheroidSurfaceAreaCandidates(numCandidates);
candidateX = spheroidPosition(lambdaCandidates, etaCandidates);
ids = farthestPointSubset(candidateX, N);
lambda = lambdaCandidates(ids);
eta = etaCandidates(ids);
end

function [lambda, eta] = spheroidSurfaceAreaCandidates(numCandidates)
k = (0:numCandidates - 1).';
u = radicalInverseSequence(numCandidates, 2);
margin = 0.5 / numCandidates;
u = min(max(u, margin), 1 - margin);
eta = inverseSpheroidAreaCoordinate(u);
goldenRatioConjugate = (sqrt(5) - 1) / 2;
lambda = mod(2 * pi * (k * goldenRatioConjugate), 2 * pi);
end

function eta = inverseSpheroidAreaCoordinate(u)
persistent cdfGrid etaGrid
if isempty(cdfGrid)
    p = spheroidParameters();
    etaGrid = linspace(-1, 1, 20001).';
    density = sqrt(p.c^2 + (p.a^2 - p.c^2) * etaGrid .^ 2);
    cdfGrid = cumtrapz(etaGrid, density);
    cdfGrid = (cdfGrid - cdfGrid(1)) ./ (cdfGrid(end) - cdfGrid(1));
end
eta = interp1(cdfGrid, etaGrid, u(:), 'linear');
eta = clampEta(eta);
end

function u = radicalInverseSequence(n, base)
indices = (0:n - 1).';
u = zeros(n, 1);
factor = 1 / base;
while any(indices > 0)
    digit = mod(indices, base);
    u = u + digit * factor;
    indices = floor(indices / base);
    factor = factor / base;
end
end

function ids = scalarSBFControlIds(lambda, eta, count)
n = numel(lambda);
count = normalizeControlPointCount(count, n);
embedding = materialEmbedding(lambda, eta);
ids = farthestPointSubset(embedding, count);
end

function info = scalarSBFControlInfo(lambda, eta, ids)
embedding = materialEmbedding(lambda, eta);
minDist2 = inf(numel(lambda), 1);
for k = 1:numel(ids)
    d2 = sum((embedding - embedding(ids(k), :)) .^ 2, 2);
    minDist2 = min(minDist2, d2);
end
info = struct();
info.controlPointCount = numel(ids);
info.controlPointFraction = numel(ids) / numel(lambda);
info.fillDistance = sqrt(max(minDist2));
end

function values = sphericalScalarSBFInterpolate(lambda, eta, data, lambdaq, etaq, controlIds)
degree = 7;
embedding = materialEmbedding(lambda, eta);
queryEmbedding = materialEmbedding(lambdaq, etaq);
controlEmbedding = embedding(controlIds, :);
dataControl = data(controlIds);
r = kp.geometry.distanceMatrix(controlEmbedding, controlEmbedding);
kernel = kp.geometry.phsKernel(r, degree);
reg = 1.0e-12 * max(1.0, max(abs(kernel), [], 'all'));
coefficients = (kernel + reg * eye(numel(controlIds))) \ dataControl(:);
rq = kp.geometry.distanceMatrix(queryEmbedding, controlEmbedding);
values = kp.geometry.phsKernel(rq, degree) * coefficients;
end

function values = localTangentPlaneScalarInterpolate(Xsource, data, Xquery, xi, theta)
op = kp.manifold.rbffdop(2, xi, theta, 0);
stencilSize = min(op.stencilSize, size(Xsource, 1));
tree = KDTreeSearcher(Xsource);
idx = knnsearch(tree, Xquery, 'K', stencilSize);
values = zeros(size(Xquery, 1), 1);
polyIndices = kp.poly.total_degree_indices(2, op.ell);
recurrence = @(N) kp.poly.jacobi_recurrence(N, 0, 0);

for i = 1:size(Xquery, 1)
    stencil = idx(i, :);
    normal = spheroidNormals(Xquery(i, :));
    R = tangentBasis(normal(:));
    xLocal = (Xsource(stencil, :) - Xquery(i, :)) * R;
    width = max(abs(xLocal), [], 'all');
    if width <= eps
        values(i) = data(stencil(1));
        continue;
    end
    xScaled = xLocal ./ width;
    r = kp.geometry.distanceMatrix(xScaled, xScaled);
    phi = (r + eps) .^ op.rbfexp;
    P = kp.poly.mpoly_eval(xScaled, polyIndices, recurrence);
    A = [[phi, P]; [P.', zeros(size(P, 2))]];
    rhs = [data(stencil); zeros(size(P, 2), 1)];
    coeff = solveInterpolationSystem(A, rhs);
    rq = sqrt(sum(xScaled .^ 2, 2)).';
    phiQuery = (rq + eps) .^ op.rbfexp;
    pQuery = kp.poly.mpoly_eval([0, 0], polyIndices, recurrence);
    values(i) = [phiQuery, pQuery] * coeff;
end
end

function coeff = solveInterpolationSystem(A, rhs)
if rcond(A) > 1.0e-13
    coeff = A \ rhs;
else
    coeff = lsqminnorm(A, rhs);
end
end

function embedding = materialEmbedding(lambda, eta)
r = sqrt(max(1 - eta(:) .^ 2, 0));
embedding = [r .* cos(lambda(:)), r .* sin(lambda(:)), eta(:)];
end

function ids = farthestPointSubset(Y, count)
n = size(Y, 1);
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

function count = normalizeControlPointCount(count, n)
if ~isfinite(count)
    count = n;
else
    count = min(n, max(1, round(count)));
end
end

function q = pointQuality(X)
tree = KDTreeSearcher(X);
[~, d] = knnsearch(tree, X, 'K', 2);
nearest = d(:, 2);
q = max(nearest) / max(min(nearest), eps);
end

function maxAnisotropy = stencilAnisotropy(X, normals)
k = min(16, size(X, 1));
tree = KDTreeSearcher(X);
idx = knnsearch(tree, X, 'K', k);
anisotropy = ones(size(X, 1), 1);
for i = 1:size(X, 1)
    normal = normals(i, :).';
    R = tangentBasis(normal);
    local = (X(idx(i, :), :) - X(i, :)) * R;
    C = local.' * local / max(size(local, 1), 1);
    ev = eig((C + C.') / 2);
    anisotropy(i) = sqrt(max(ev) / max(min(ev), eps));
end
maxAnisotropy = max(anisotropy);
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

function [willCross, triggerMetric] = predictQualityCrossing( ...
    qualityValues, currentIndex, threshold, lookaheadSteps)
current = qualityValues(currentIndex);
previous = qualityValues(max(1, currentIndex - 1));
slope = max(0, current - previous);
triggerMetric = current + lookaheadSteps * slope;
willCross = triggerMetric >= threshold;
end

function [predictedCrossTime, actualCrossTime, leadSteps] = predictorSummary( ...
    qualityTimes, qualityValues, threshold, lookaheadSteps)
numValues = numel(qualityValues);
predictedIndex = NaN;
actualIndex = NaN;
for k = 2:numValues
    [willCross, ~] = predictQualityCrossing(qualityValues, k, threshold, lookaheadSteps);
    if isnan(predictedIndex) && willCross
        predictedIndex = k;
    end
    if isnan(actualIndex) && qualityValues(k) >= threshold
        actualIndex = k;
    end
end

if isnan(predictedIndex)
    predictedCrossTime = NaN;
else
    predictedCrossTime = qualityTimes(predictedIndex);
end
if isnan(actualIndex)
    actualCrossTime = NaN;
else
    actualCrossTime = qualityTimes(actualIndex);
end
if isnan(predictedIndex) || isnan(actualIndex)
    leadSteps = NaN;
else
    leadSteps = actualIndex - predictedIndex;
end
end

function area = spheroidAreaApprox()
persistent cachedArea
if isempty(cachedArea)
    p = spheroidParameters();
    m = 400;
    eta = linspace(-1 + 1 / m, 1 - 1 / m, m).';
    metric = p.c^2 + (p.a^2 - p.c^2) * eta .^ 2;
    sqrtg = p.a * sqrt(metric);
    cachedArea = 2 * pi * trapz(eta, sqrtg);
end
area = cachedArea;
end

function p = spheroidParameters()
p = struct('a', 1.0, 'c', 0.62);
end

function flow = spheroidFlowParameters(scale)
flow = struct();
flow.Omega = 0.75 * scale;
flow.S = 2.40 * scale;
flow.A = 0.0;
flow.B = 0.35 * scale;
end

function result = emptyRunResult()
result = struct( ...
    'mode', "", ...
    'transferMode', "", ...
    'relerr', NaN, ...
    'initialMass', NaN, ...
    'finalMass', NaN, ...
    'exactFinalMass', NaN, ...
    'massRelError', NaN, ...
    'balanceRelResidual', NaN, ...
    'massCorrectionSteps', NaN, ...
    'massCorrectionMaxAbs', NaN, ...
    'massCorrectionMaxRelative', NaN, ...
    'massCorrectionTotalAbs', NaN, ...
    'massCorrectionFinalAbs', NaN, ...
    'massCorrectionMaxPointShift', NaN, ...
    'numRearrangements', 0, ...
    'numPredictiveRearrangements', 0, ...
    'numThresholdRearrangements', 0, ...
    'numFixedRearrangements', 0, ...
    'numLateRearrangements', 0, ...
    'initialQuality', NaN, ...
    'finalQuality', NaN, ...
    'maxQuality', NaN, ...
    'initialAnisotropy', NaN, ...
    'finalAnisotropy', NaN, ...
    'maxAnisotropy', NaN, ...
    'sbfControlPointCount', NaN, ...
    'sbfControlPointFraction', NaN, ...
    'sbfFillDistance', NaN, ...
    'transferStencilSize', NaN, ...
    'predictedQualityCrossingTime', NaN, ...
    'actualQualityCrossingTime', NaN, ...
    'predictorLeadSteps', NaN, ...
    'rearrangementTimes', {cell(1, 1)}, ...
    'rearrangementTriggerValues', {cell(1, 1)}, ...
    'rearrangementQualityValues', {cell(1, 1)}, ...
    'rearrangementWasPredicted', {cell(1, 1)}, ...
    'rearrangementWasThresholdTriggered', {cell(1, 1)}, ...
    'rearrangementWasFixedScheduled', {cell(1, 1)}, ...
    'qualityTimes', {cell(1, 1)}, ...
    'qualityValues', {cell(1, 1)}, ...
    'anisotropyValues', {cell(1, 1)}, ...
    'finalLambda', {cell(1, 1)}, ...
    'finalEta', {cell(1, 1)}, ...
    'finalX', {cell(1, 1)}, ...
    'finalNumerical', {cell(1, 1)}, ...
    'finalExact', {cell(1, 1)}, ...
    'operatorStats', {cell(1, 1)});
end

function stats = emptyOperatorStats()
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
    'neighborReuseRejected', 0);
end

function writeStudyFigure(study, imagePath)
if ~exist(fileparts(imagePath), 'dir')
    mkdir(fileparts(imagePath));
end
fig = figure('Color', 'w', 'Position', [100, 100, 1380, 540]);
tiledlayout(fig, 1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
sgtitle(fig, sprintf('Forced advection with SL rearrangement on a spheroid (N=%d, \\xi=%d)', ...
    study.N, study.xi));

nexttile;
hold on;
for k = 1:numel(study.results)
    plot(study.results(k).qualityTimes{1}, study.results(k).qualityValues{1}, ...
        'LineWidth', 1.6, 'DisplayName', study.results(k).mode);
    rt = study.results(k).rearrangementTimes{1};
    yl = ylim;
    for j = 1:numel(rt)
        plot([rt(j), rt(j)], yl, ':', 'Color', [0.3, 0.3, 0.3], ...
            'HandleVisibility', 'off');
    end
end
hold off;
grid on;
xlabel('time');
ylabel('nearest-neighbor quality');
legend('Location', 'northwest');
title('Predictive rearrangement trigger');

allExactCells = vertcat(study.results.finalExact);
allExact = vertcat(allExactCells{:});
colorLimits = [min(allExact(:)), max(allExact(:))];

for k = 1:numel(study.results)
    nexttile;
    X = study.results(k).finalX{1};
    scatter3(X(:, 1), X(:, 2), X(:, 3), 18, study.results(k).finalNumerical{1}, ...
        'filled');
    axis equal off;
    view(40, 24);
    title(sprintf('%s: E=%.2e, q=%.2f', study.results(k).mode, ...
        study.results(k).relerr, study.results(k).finalQuality), ...
        'Interpreter', 'none');
    clim(colorLimits);
    colorbar;
end
axesList = findall(fig, 'Type', 'axes');
for iax = 1:numel(axesList)
    try
        axtoolbar(axesList(iax), {});
    catch
    end
end
exportgraphics(fig, imagePath, 'Resolution', 220);
close(fig);
end
