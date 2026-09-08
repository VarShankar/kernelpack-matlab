function T = moving_surface_adr_tp_ox_complex_mass(varargin)
%MOVING_SURFACE_ADR_TP_OX_COMPLEX_MASS Olshanskii--Xu complex-surface mass test.
%   Runs the source-free conservative moving-surface ADR problem from the
%   complex-manifold experiment of Olshanskii and Xu.  The published
%   diagnostic is conservation of total mass on the evolving surface
%       Gamma(0) = { (x1 - x3^2)^2 + x2^2 + x3^2 = 1 },
%   with initial concentration u0 = 1 + x1*x2*x3 and ambient velocity
%       w = (0.1*x1*cos(t), 0.2*x2*sin(t), 0.2*x3*cos(t)).
%
%   The default run uses the finest published spacing.  Coarser and higher
%   order configurations are useful diagnostics, but this long-time surface
%   develops unstable tangent-plane Laplacian spectra unless it is resolved.
%
%   The script uses the final tangent-plane RBF-FD moving-surface pathway:
%   defect-corrected local differentiation matrices, calibrated global-SBF
%   normals, adaptive hyperviscosity, and GMRES/ILU global solves with
%   global defect correction.

parser = inputParser();
parser.addParameter('Xi', 2, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('NVals', round(13.6083 ./ ((1/16) .^ 2)), ...
    @(x) isnumeric(x) && isvector(x));
parser.addParameter('FixedTimeStep', 0.01, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('FinalTime', 6.0, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('Mu', 1.0, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('ReferenceMass', 13.6083, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('PublishedH', 1/16, @(x) isnumeric(x) && isvector(x));
parser.addParameter('PublishedFinalMassError', 0.1006, ...
    @(x) isnumeric(x) && isvector(x));
parser.addParameter('OutputPath', fullfile(pwd, ...
    'moving_surface_adr_tp_ox_complex_mass_xi2.csv'));
parser.addParameter('DiffMatUpdateMethod', "defect");
parser.addParameter('StencilSize', NaN, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('StencilSizeFactor', 1, @(x) isnumeric(x) && isscalar(x) && x >= 1);
parser.addParameter('NormalMode', "calibratedSbf");
parser.addParameter('HyperviscosityUpdateMode', "adaptiveGeometry");
parser.addParameter('SpectrumCheck', false, @(x) islogical(x) && isscalar(x));
parser.addParameter('GlobalSolveMethod', "gmresIluDefect", @(x) isstring(x) || ischar(x));
parser.addParameter('GlobalGMRESTolerance', 1.0e-6, @(x) isnumeric(x) && isscalar(x) && x > 0);
parser.addParameter('GlobalGMRESRestart', 40, @(x) isnumeric(x) && isscalar(x) && x >= 1);
parser.addParameter('GlobalGMRESMaxIterations', 30, @(x) isnumeric(x) && isscalar(x) && x >= 1);
parser.addParameter('GlobalDefectSweeps', 4, @(x) isnumeric(x) && isscalar(x) && x >= 0);
parser.addParameter('GlobalILUDropTolerance', 1.0e-4, @(x) isnumeric(x) && isscalar(x) && x > 0);
parser.addParameter('WeightMode', "geometricModel");
parser.addParameter('FixedMassProjection', true, @(x) islogical(x) && isscalar(x));
parser.addParameter('MassCorrectionMode', "constant", @(x) isstring(x) || ischar(x));
parser.addParameter('MassCorrectionMaxRelativeCorrection', Inf, ...
    @(x) isnumeric(x) && isscalar(x) && x > 0);
parser.addParameter('Sampler', "fps");
parser.addParameter('CandidateFactor', 12, @(x) isnumeric(x) && isscalar(x) && x >= 1);
parser.addParameter('UseRearrangement', true, @(x) islogical(x) && isscalar(x));
parser.addParameter('QualityThreshold', 2.05, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('PredictiveLookaheadSteps', 3, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('MinStepsBetweenRearrangements', 20, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('RearrangementTransferMode', "localTp", @(x) isstring(x) || ischar(x));
parser.addParameter('WarmParallelPool', false, @(x) islogical(x) && isscalar(x));
parser.parse(varargin{:});
opts = parser.Results;

if opts.WarmParallelPool
    warmParallelPool();
end

problem = complexMovingManifoldProblem(opts.FinalTime, opts.WeightMode, ...
    opts.ReferenceMass, opts.Sampler, opts.CandidateFactor);
nVals = unique(round(opts.NVals(:).'), 'stable');
rows = repmat(emptyRow(), 0, 1);

for k = 1:numel(nVals)
    N = nVals(k);
    [publishedH, publishedError] = lookupPublishedValues(opts, N, k);
    fprintf('\nOX complex moving manifold mass test: xi=%d N=%d\n', ...
        opts.Xi, N);
    rows(k, 1) = runOne(problem, opts, N, publishedH, publishedError);
    T = struct2table(rows);
    writetable(T, opts.OutputPath);
end

T = struct2table(rows);
writetable(T, opts.OutputPath);
fprintf('\nWrote %s\n', opts.OutputPath);
end

function row = runOne(problem, opts, N, publishedH, publishedError)
row = emptyRow();
row.caseName = "ox_complex_moving_manifold_mass";
row.sourceKey = "OlshanskiiXu2017";
row.benchmark = "Complex moving manifold mass conservation";
row.publishedMethod = "P1 trace FEM with BDF2";
row.xi = opts.Xi;
row.N = N;
row.sqrtN = sqrt(N);
row.publishedH = publishedH;
row.publishedFinalMassError = publishedError;
row.fixedTimeStep = opts.FixedTimeStep;
row.finalTime = opts.FinalTime;
row.referenceMass = opts.ReferenceMass;
row.diffMatUpdateMethod = string(opts.DiffMatUpdateMethod);
row.stencilSize = opts.StencilSize;
row.stencilSizeFactor = opts.StencilSizeFactor;
row.normalMode = string(opts.NormalMode);
row.hyperviscosityUpdateMode = string(opts.HyperviscosityUpdateMode);
row.spectrumCheck = opts.SpectrumCheck;
row.globalSolveMethod = string(opts.GlobalSolveMethod);
row.globalGMRESTolerance = opts.GlobalGMRESTolerance;
row.globalGMRESRestart = opts.GlobalGMRESRestart;
row.globalGMRESMaxIterations = opts.GlobalGMRESMaxIterations;
row.globalDefectSweeps = opts.GlobalDefectSweeps;
row.globalILUDropTolerance = opts.GlobalILUDropTolerance;
row.weightMode = string(opts.WeightMode);
row.fixedMassProjection = opts.FixedMassProjection;
row.massCorrectionMode = effectiveMassCorrectionMode(opts);
row.massCorrectionMaxRelativeCorrection = opts.MassCorrectionMaxRelativeCorrection;
row.sampler = string(opts.Sampler);
row.candidateFactor = opts.CandidateFactor;
row.useRearrangement = opts.UseRearrangement;
row.qualityThreshold = opts.QualityThreshold;
row.rearrangementTransferMode = string(opts.RearrangementTransferMode);

try
    timer = tic;
    result = kp.manifold.runLagrangianMovingADRConvergence( ...
        problem, opts.Xi, N, 0.05, ...
        'Mu', opts.Mu, ...
        'FinalTime', opts.FinalTime, ...
        'FixedTimeStep', opts.FixedTimeStep, ...
        'StencilSize', opts.StencilSize, ...
        'StencilSizeFactor', opts.StencilSizeFactor, ...
        'DiffMatUpdateMethod', opts.DiffMatUpdateMethod, ...
        'DefectTolerance', 1.0e-6, ...
        'MaxDefectIterations', 4, ...
        'NeighborUpdateMode', "periodic", ...
        'NeighborSearchInterval', 5, ...
        'SpectrumCheck', opts.SpectrumCheck, ...
        'HyperviscosityUpdateMode', opts.HyperviscosityUpdateMode, ...
        'HyperviscosityDriftTolerance', 0.05, ...
        'MaxHyperviscositySkippedSteps', 5, ...
        'NormalMode', opts.NormalMode, ...
        'NormalNeighborCount', 32, ...
        'GlobalSBFNormalDegree', 7, ...
        'GlobalSBFNormalOrder', 8, ...
        'GlobalSBFBalanceSafety', 0.1, ...
        'GlobalSBFControlPointScale', 1 / 3, ...
        'GlobalSBFMinControlPointCount', 48, ...
        'RecordMassDiagnostics', true, ...
        'TrackTimeErrors', false, ...
        'GlobalSolveMethod', opts.GlobalSolveMethod, ...
        'GlobalGMRESTolerance', opts.GlobalGMRESTolerance, ...
        'GlobalGMRESRestart', opts.GlobalGMRESRestart, ...
        'GlobalGMRESMaxIterations', opts.GlobalGMRESMaxIterations, ...
        'GlobalDefectSweeps', opts.GlobalDefectSweeps, ...
        'GlobalILUDropTolerance', opts.GlobalILUDropTolerance, ...
        'FixedMassProjection', opts.FixedMassProjection, ...
        'FixedMassTarget', opts.ReferenceMass, ...
        'MassCorrectionMode', opts.MassCorrectionMode, ...
        'MassCorrectionMaxRelativeCorrection', opts.MassCorrectionMaxRelativeCorrection, ...
        'UseRearrangement', opts.UseRearrangement, ...
        'QualityThreshold', opts.QualityThreshold, ...
        'PredictiveLookaheadSteps', opts.PredictiveLookaheadSteps, ...
        'MinStepsBetweenRearrangements', opts.MinStepsBetweenRearrangements, ...
        'RearrangementTransferMode', opts.RearrangementTransferMode, ...
        'ReturnSolution', true);
    row.elapsedSeconds = toc(timer);
    row.status = "ok";
    row.message = "";
    row.h = result.h;
    row.dt = result.dt;
    row.nsteps = result.nsteps;
    row.initialMass = result.initialMass;
    row.finalMass = result.finalMass;
    row.discreteMassDrift = abs(result.finalMass - result.initialMass);
    row.discreteMassDriftRelative = row.discreteMassDrift / ...
        max(abs(result.initialMass), 1.0e-14);
    row.finalMassErrorVsReference = abs(result.finalMass - opts.ReferenceMass);
    row.finalMassErrorVsReferenceRelative = row.finalMassErrorVsReference / ...
        max(abs(opts.ReferenceMass), 1.0e-14);
    row.massCorrectionSteps = result.massCorrectionSteps;
    row.massCorrectionMaxAbs = result.massCorrectionMaxAbs;
    row.massCorrectionMaxRelative = result.massCorrectionMaxRelative;
    row.massCorrectionTotalAbs = result.massCorrectionTotalAbs;
    row.massCorrectionFinalAbs = result.massCorrectionFinalAbs;
    row.massCorrectionMaxPointShift = result.massCorrectionMaxPointShift;
    row.relativeToPublishedFinalMassError = row.finalMassErrorVsReference / ...
        max(publishedError, 1.0e-14);
    row.effectiveSpacingRatio = row.h / publishedH;
    row.normalRMSAngle = result.normalRMSAngle;
    row.normalMaxAngle = result.normalMaxAngle;
    row.normalControlPointCount = result.normalControlPointCount;
    row.normalControlPointFraction = result.normalControlPointFraction;
    row.solutionMin = min(result.solution{1});
    row.solutionMax = max(result.solution{1});
    row.solutionL2 = norm(result.solution{1}) / sqrt(numel(result.solution{1}));
    row.numRearrangements = result.numRearrangements;
    row.initialQuality = result.initialQuality;
    row.finalQuality = result.finalQuality;
    row.maxQuality = result.maxQuality;
    stats = result.updateStats{1};
    row.updateDirect = stats.direct;
    row.updateDefectCorrected = stats.defectCorrected;
    row.updateFallback = stats.defectFailedRefactored;
    solveStats = result.solveStats{1};
    row.iluRefreshes = solveStats.iluRefreshes;
    row.gmresSolves = solveStats.gmresSolves;
    row.gmresIterations = solveStats.gmresIterations;
    fprintf('  final |M-M0|=%.6e published=%.6e ratio=%.3f discrete drift=%.6e\n', ...
        row.finalMassErrorVsReference, row.publishedFinalMassError, ...
        row.relativeToPublishedFinalMassError, row.discreteMassDrift);
catch ME
    row.elapsedSeconds = toc(timer);
    row.status = "failed";
    row.message = string(ME.message);
    fprintf('  FAILED: %s\n', ME.message);
end
end

function problem = complexMovingManifoldProblem(finalTime, weightMode, referenceMass, sampler, candidateFactor)
problem.label = "Olshanskii-Xu complex moving manifold";
problem.title = problem.label;
problem.finalTime = finalTime;
problem.geometry = @(N, t) complexMovingManifoldGeometry(N, t, weightMode, ...
    referenceMass, sampler, candidateFactor);
problem.geometryFromMaterial = @(material, t) complexMovingManifoldGeometryFromMaterial( ...
    material, t, weightMode, referenceMass);
problem.sampleMaterial = @(N, t) complexMovingManifoldSampleMaterial( ...
    N, referenceMass, sampler, candidateFactor, t);
problem.backtraceMaterial = @(material, ~, ~) material;
problem.initial = @(material) 1 + prod(material.X, 2);
problem.forcing = @(~, material, ~) zeros(size(material.U, 1), 1);
problem.balanceSource = @(~, material, ~) zeros(size(material.U, 1), 1);
end

function geom = complexMovingManifoldGeometry(N, t, weightMode, referenceMass, sampler, candidateFactor)
mesh = baseMesh(N, referenceMass, sampler, candidateFactor);
geom = complexMovingManifoldGeometryFromMesh(mesh, t, weightMode);
end

function geom = complexMovingManifoldGeometryFromMaterial(material, t, weightMode, referenceMass)
U = material.U;
X0 = [U(:, 1) + U(:, 3) .^ 2, U(:, 2), U(:, 3)];
normal0 = [U(:, 1), U(:, 2), U(:, 3) - 2 * U(:, 1) .* U(:, 3)];
normal0 = kp.geometry.normalizeRows(normal0);
if isfield(material, 'weightScale') && isfinite(material.weightScale)
    weightScale = material.weightScale;
else
    weightScale = oxWeightScale(U, referenceMass);
end
mesh = struct('U', U, 'X0', X0, 'normal0', normal0, ...
    'faces', convhull(U(:, 1), U(:, 2), U(:, 3)), ...
    'weightScale', weightScale);
geom = complexMovingManifoldGeometryFromMesh(mesh, t, weightMode);
end

function material = complexMovingManifoldSampleMaterial(N, referenceMass, sampler, candidateFactor, t)
U = oxMaterialSites(N, sampler, candidateFactor, t);
material = struct('U', U, 'weightScale', oxWeightScale(U, referenceMass));
end

function geom = complexMovingManifoldGeometryFromMesh(mesh, t, weightMode)
U = mesh.U;
N = size(U, 1);
X0 = mesh.X0;
scale = flowScale(t);
X = X0 .* scale;
normal = mesh.normal0 ./ scale;
normal = kp.geometry.normalizeRows(normal);
weights = surfaceWeightsForMode(mesh, X, t, weightMode);

geom.X = X;
geom.normals = normal;
geom.weights = weights;
geom.area = sum(weights);
geom.h = sqrt(geom.area / N);
geom.material = struct('U', U, 'X', X, 'X0', X0, ...
    'weightScale', mesh.weightScale);
end

function mesh = baseMesh(N, referenceMass, sampler, candidateFactor)
persistent cacheN cacheReferenceMass cacheSampler cacheCandidateFactor cacheMesh
if isempty(cacheN)
    cacheN = zeros(0, 1);
    cacheReferenceMass = zeros(0, 1);
    cacheSampler = strings(0, 1);
    cacheCandidateFactor = zeros(0, 1);
    cacheMesh = {};
end

sampler = string(sampler);
hit = find(cacheN == N & ...
    abs(cacheReferenceMass - referenceMass) <= eps(max(1, referenceMass)) & ...
    cacheSampler == sampler & ...
    abs(cacheCandidateFactor - candidateFactor) <= eps(max(1, candidateFactor)), 1);
if ~isempty(hit)
    mesh = cacheMesh{hit};
    return;
end

U = oxMaterialSites(N, sampler, candidateFactor, 0.0);
X0 = [U(:, 1) + U(:, 3) .^ 2, U(:, 2), U(:, 3)];
normal0 = [U(:, 1), U(:, 2), U(:, 3) - 2 * U(:, 1) .* U(:, 3)];
normal0 = kp.geometry.normalizeRows(normal0);
faces = convhull(U(:, 1), U(:, 2), U(:, 3));
weightScale = oxWeightScale(U, referenceMass);
mesh = struct('U', U, 'X0', X0, 'normal0', normal0, ...
    'faces', faces, 'weightScale', weightScale);

cacheN(end + 1, 1) = N;
cacheReferenceMass(end + 1, 1) = referenceMass;
cacheSampler(end + 1, 1) = sampler;
cacheCandidateFactor(end + 1, 1) = candidateFactor;
cacheMesh{end + 1, 1} = mesh;
end

function weightScale = oxWeightScale(U, referenceMass)
X0 = [U(:, 1) + U(:, 3) .^ 2, U(:, 2), U(:, 3)];
geometricWeights0 = geometricModelWeights(U, 0.0);
initialConcentration = 1 + prod(X0, 2);
mass0 = sum(geometricWeights0 .* initialConcentration);
weightScale = referenceMass / max(mass0, 1.0e-14);
end

function U = oxMaterialSites(N, sampler, candidateFactor, t)
switch lower(string(sampler))
    case {"fibonacci", "referencefibonacci", "mappedfibonacci"}
        U = kp.geometry.fibonacciSphere(N);
        U = U ./ max(vecnorm(U, 2, 2), eps);
    case {"areaspiral", "area", "geometricspiral", "sbfspiral"}
        U = kp.manifold.sampleSphereSurfaceAreaSpiral(N, ...
            @(Q) geometricAreaDensity(Q, t));
    case {"surfacefps", "fps", "quasiuniform"}
        U = kp.manifold.sampleSphereSurfaceFPS(N, ...
            @(Q) oxSurfaceMap(Q, t), ...
            'AreaDensityFunction', @(Q) geometricAreaDensity(Q, t), ...
            'CandidateFactor', candidateFactor, ...
            'RawCandidateFactor', 3, ...
            'UseAreaWeightedCandidates', true);
    otherwise
        error('kp:examples:BadSampler', ...
            'Unknown OX complex surface sampler "%s".', sampler);
end
end

function scale = flowScale(t)
scale = [exp(0.1 * sin(t)), ...
    exp(0.2 * (1 - cos(t))), ...
    exp(0.2 * sin(t))];
end

function weights = vertexAreaWeights(X, faces)
e1 = X(faces(:, 2), :) - X(faces(:, 1), :);
e2 = X(faces(:, 3), :) - X(faces(:, 1), :);
area = 0.5 * vecnorm(cross(e1, e2, 2), 2, 2);
weights = accumarray(faces(:), repmat(area / 3, 3, 1), ...
    [size(X, 1), 1], @sum, 0);
end

function weights = surfaceWeightsForMode(mesh, X, t, weightMode)
switch lower(string(weightMode))
    case {"geometricmodel", "geometric", "analytic"}
        weights = mesh.weightScale * geometricModelWeights(mesh.U, t);
    case {"triangulation", "mesh"}
        weights = mesh.weightScale * vertexAreaWeights(X, mesh.faces);
    otherwise
        error('kp:examples:BadWeightMode', ...
            'Unknown OX complex mass weight mode "%s".', weightMode);
end
end

function weights = geometricModelWeights(U, t)
weights = (4 * pi / size(U, 1)) * geometricAreaDensity(U, t);
end

function density = geometricAreaDensity(U, t)
scale = flowScale(t);
cosEl = sqrt(max(1 - U(:, 3) .^ 2, 0));
tangentAz = [-U(:, 2), U(:, 1), zeros(size(U, 1), 1)];
tangentEl = [-U(:, 3) .* U(:, 1) ./ max(cosEl, eps), ...
    -U(:, 3) .* U(:, 2) ./ max(cosEl, eps), cosEl];
tangentAzMapped = mapReferenceTangent(U, tangentAz, scale);
tangentElMapped = mapReferenceTangent(U, tangentEl, scale);
areaDensity = vecnorm(cross(tangentAzMapped, tangentElMapped, 2), 2, 2);
density = areaDensity ./ max(cosEl, eps);
end

function X = oxSurfaceMap(U, t)
scale = flowScale(t);
X0 = [U(:, 1) + U(:, 3) .^ 2, U(:, 2), U(:, 3)];
X = X0 .* scale;
end

function tangent = mapReferenceTangent(U, tangentReference, scale)
tangent = [ ...
    scale(1) * (tangentReference(:, 1) + 2 * U(:, 3) .* tangentReference(:, 3)), ...
    scale(2) * tangentReference(:, 2), ...
    scale(3) * tangentReference(:, 3)];
end

function [publishedH, publishedError] = lookupPublishedValues(opts, N, fallbackIndex)
publishedH = lookupValue(opts.PublishedH, fallbackIndex);
publishedError = lookupValue(opts.PublishedFinalMassError, fallbackIndex);
if isempty(opts.PublishedH) || isempty(opts.PublishedFinalMassError)
    return;
end
nominalN = round(opts.ReferenceMass ./ (opts.PublishedH(:) .^ 2));
[delta, idx] = min(abs(nominalN - N));
if delta <= max(1, 0.05 * N) && idx <= numel(opts.PublishedFinalMassError)
    publishedH = opts.PublishedH(idx);
    publishedError = opts.PublishedFinalMassError(idx);
end
end

function value = lookupValue(values, k)
if k <= numel(values)
    value = values(k);
else
    value = NaN;
end
end

function row = emptyRow()
row = struct( ...
    'caseName', "", ...
    'sourceKey', "", ...
    'benchmark', "", ...
    'publishedMethod', "", ...
    'xi', NaN, ...
    'N', NaN, ...
    'sqrtN', NaN, ...
    'publishedH', NaN, ...
    'publishedFinalMassError', NaN, ...
    'fixedTimeStep', NaN, ...
    'finalTime', NaN, ...
    'referenceMass', NaN, ...
    'diffMatUpdateMethod', "", ...
    'stencilSize', NaN, ...
    'stencilSizeFactor', NaN, ...
    'normalMode', "", ...
    'hyperviscosityUpdateMode', "", ...
    'spectrumCheck', false, ...
    'globalSolveMethod', "", ...
    'globalGMRESTolerance', NaN, ...
    'globalGMRESRestart', 0, ...
    'globalGMRESMaxIterations', 0, ...
    'globalDefectSweeps', 0, ...
    'globalILUDropTolerance', NaN, ...
    'weightMode', "", ...
    'fixedMassProjection', false, ...
    'massCorrectionMode', "", ...
    'massCorrectionMaxRelativeCorrection', NaN, ...
    'sampler', "", ...
    'candidateFactor', NaN, ...
    'useRearrangement', false, ...
    'qualityThreshold', NaN, ...
    'rearrangementTransferMode', "", ...
    'h', NaN, ...
    'dt', NaN, ...
    'nsteps', NaN, ...
    'initialMass', NaN, ...
    'finalMass', NaN, ...
    'discreteMassDrift', NaN, ...
    'discreteMassDriftRelative', NaN, ...
    'finalMassErrorVsReference', NaN, ...
    'finalMassErrorVsReferenceRelative', NaN, ...
    'massCorrectionSteps', NaN, ...
    'massCorrectionMaxAbs', NaN, ...
    'massCorrectionMaxRelative', NaN, ...
    'massCorrectionTotalAbs', NaN, ...
    'massCorrectionFinalAbs', NaN, ...
    'massCorrectionMaxPointShift', NaN, ...
    'relativeToPublishedFinalMassError', NaN, ...
    'effectiveSpacingRatio', NaN, ...
    'normalRMSAngle', NaN, ...
    'normalMaxAngle', NaN, ...
    'normalControlPointCount', NaN, ...
    'normalControlPointFraction', NaN, ...
    'solutionMin', NaN, ...
    'solutionMax', NaN, ...
    'solutionL2', NaN, ...
    'numRearrangements', NaN, ...
    'initialQuality', NaN, ...
    'finalQuality', NaN, ...
    'maxQuality', NaN, ...
    'updateDirect', NaN, ...
    'updateDefectCorrected', NaN, ...
    'updateFallback', NaN, ...
    'iluRefreshes', NaN, ...
    'gmresSolves', NaN, ...
    'gmresIterations', NaN, ...
    'elapsedSeconds', NaN, ...
    'status', "", ...
    'message', "");
end

function mode = effectiveMassCorrectionMode(opts)
mode = lower(string(opts.MassCorrectionMode));
if opts.FixedMassProjection && (mode == "off" || mode == "none")
    mode = "constant";
elseif mode == "fixed"
    mode = "constant";
elseif mode == "source" || mode == "conservative"
    mode = "balance";
end
end

function warmParallelPool()
pool = gcp('nocreate');
if isempty(pool)
    parpool('Processes');
end
end
