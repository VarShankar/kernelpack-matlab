function study = moving_surface_adr_tp_ling_comparison(varargin)
%MOVING_SURFACE_ADR_TP_LING_COMPARISON Compare against Ling surface-PDE papers.
%   The direct numerical comparison is the oscillating-ellipsoid benchmark
%   from Petras, Ling, Piret, and Ruuth (JCP 2019, Table 2).  The remaining
%   Ling benchmarks are recorded in the inventory with the reason they are or
%   are not directly comparable to the scalar moving-surface ADR solver.

parser = inputParser();
parser.addParameter('DeltaX', [0.2, 0.1, 0.05], @(x) isnumeric(x) && (isempty(x) || isvector(x)));
parser.addParameter('Times', [0.08, 0.16], @(x) isnumeric(x) && (isempty(x) || isvector(x)));
parser.addParameter('ChenLingH', [3/4, 1/2, 1/4, 1/8, 1/16], ...
    @(x) isnumeric(x) && (isempty(x) || isvector(x)));
parser.addParameter('Xi', 2, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('DtScale', 0.05, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('ContinueOnError', true, @(x) islogical(x) && isscalar(x));
parser.addParameter('SaveAfterEach', true, @(x) islogical(x) && isscalar(x));
parser.addParameter('SolverProfile', "final", ...
    @(x) any(strcmpi(string(x), ["final", "minimal"])));
parser.addParameter('ExactStartup', true, @(x) islogical(x) && isscalar(x));
parser.addParameter('GlobalGMRESTolerance', 1.0e-12, ...
    @(x) isnumeric(x) && isscalar(x) && x > 0);
parser.addParameter('OutputPrefix', fullfile(pwd, ...
    'moving_surface_adr_tp_ling_petras_comparison'));
parser.parse(varargin{:});
opts = parser.Results;

published = petrasLingOscillatingEllipsoidTable();
chenLingPublished = chenLingMovingEllipsoidTable();
inventory = lingBenchmarkInventory();
rows = repmat(emptyRow(), 0, 1);

rowIndex = 0;
for dx = opts.DeltaX(:).'
    for finalTime = opts.Times(:).'
        rowIndex = rowIndex + 1;
        rows(rowIndex, 1) = runPetrasLingCase(dx, finalTime, opts, published);
        if opts.SaveAfterEach
            study = buildStudy(rows, published, chenLingPublished, inventory, opts);
            saveStudy(study, opts.OutputPrefix);
        end
    end
end

for h = opts.ChenLingH(:).'
    rowIndex = rowIndex + 1;
    rows(rowIndex, 1) = runChenLingCase(h, opts, chenLingPublished);
    if opts.SaveAfterEach
        study = buildStudy(rows, published, chenLingPublished, inventory, opts);
        saveStudy(study, opts.OutputPrefix);
    end
end

study = buildStudy(rows, published, chenLingPublished, inventory, opts);
saveStudy(study, opts.OutputPrefix);
end

function row = runPetrasLingCase(dx, finalTime, opts, published)
row = emptyRow();
row.sourceKey = "PetrasLingPiretRuuth2019";
row.comparisonFamily = "least-squares implicit RBF-CPM";
row.benchmark = "Oscillating ellipsoid";
row.metric = "relative final ell_infty error";
row.publishedDeltaX = dx;
row.finalTime = finalTime;
row.publishedDt = 2 * dx^2;
row.xi = opts.Xi;

match = published(published.deltaX == dx & published.finalTime == finalTime, :);
if ~isempty(match)
    row.publishedError = match.relativeLInfError(1);
end

problem = petrasLingOscillatingEllipsoidProblem(finalTime);
geom0 = problem.geometryFromMaterial(problem.sampleMaterial(384, 0), 0);
area0 = geom0.area;
N = max(48, round(area0 / dx^2));
row.N = N;
row.sqrtN = sqrt(N);

fprintf('\nPetras-Ling oscillating ellipsoid: dx=%.5g T=%.5g xi=%d N=%d\n', ...
    dx, finalTime, opts.Xi, N);
try
    tic;
    result = kp.manifold.runLagrangianMovingADRConvergence( ...
        problem, opts.Xi, N, opts.DtScale, ...
        'Mu', 1.0, ...
        'FinalTime', finalTime, ...
        'FixedTimeStep', row.publishedDt, ...
        'DiffMatUpdateMethod', 'defect', ...
        'DefectTolerance', 1.0e-6, ...
        'MaxDefectIterations', 4, ...
        'NeighborUpdateMode', "periodic", ...
        'NeighborSearchInterval', 5, ...
        'SpectrumCheck', false, ...
        'HyperviscosityUpdateMode', 'adaptiveGeometry', ...
        'HyperviscosityDriftTolerance', 0.05, ...
        'MaxHyperviscositySkippedSteps', 5, ...
        'NormalMode', 'calibratedSbf', ...
        'NormalNeighborCount', 32, ...
        'GlobalSBFNormalDegree', 7, ...
        'GlobalSBFNormalOrder', 8, ...
        'GlobalSBFBalanceSafety', 0.1, ...
        'GlobalSBFControlPointScale', 1 / 3, ...
        'GlobalSBFMinControlPointCount', 48, ...
        'RecordMassDiagnostics', true, ...
        'TrackTimeErrors', true, ...
        'ExactStartup', opts.ExactStartup, ...
        'GlobalSolveMethod', "gmresIluDefect", ...
        'GlobalGMRESTolerance', opts.GlobalGMRESTolerance, ...
        'GlobalGMRESRestart', 40, ...
        'GlobalGMRESMaxIterations', 30, ...
        'GlobalDefectSweeps', 4, ...
        'GlobalILUDropTolerance', 1.0e-4, ...
        'UseRearrangement', true, ...
        'QualityThreshold', 1.55, ...
        'PredictiveLookaheadSteps', 3, ...
        'MinStepsBetweenRearrangements', 8, ...
        'RearrangementTransferMode', "localTp", ...
        'MassCorrectionMode', "balance");
    row.elapsedSeconds = toc;
    row.status = "ok";
    row.h = result.h;
    row.dt = result.dt;
    row.nsteps = result.nsteps;
    row.ourRelativeLInfError = result.finalRelLInfError;
    row.ourRelativeL2Error = result.finalRelL2Error;
    row.maxRelativeL2Error = result.maxRelL2Error;
    row.timeRelativeL2L2Error = result.timeRelL2L2Error;
    row.massRelError = result.massRelError;
    row.balanceRelResidual = result.balanceRelResidual;
    row.normalRMSAngle = result.normalRMSAngle;
    row.normalControlPointCount = result.normalControlPointCount;
    row.numRearrangements = result.numRearrangements;
    row.relativeToPublished = row.ourRelativeLInfError / row.publishedError;
    row.hOverDeltaX = row.h / dx;
    row.dtOverPublishedDt = row.dt / row.publishedDt;
    fprintf('  rel Linf=%.6e published=%.6e ratio=%.3f h/dx=%.3f dt/dt_lit=%.3f\n', ...
        row.ourRelativeLInfError, row.publishedError, row.relativeToPublished, ...
        row.hOverDeltaX, row.dtOverPublishedDt);
catch ME
    row.elapsedSeconds = toc;
    row.status = "failed";
    row.message = string(ME.message);
    fprintf('  FAILED: %s\n', ME.message);
    if ~opts.ContinueOnError
        rethrow(ME);
    end
end
end

function row = runChenLingCase(h, opts, published)
row = emptyRow();
row.sourceKey = "ChenLing2020";
row.comparisonFamily = "extrinsic kernel collocation";
row.benchmark = "Moving ellipsoid";
row.metric = "Linf_t H1/H2 seminorm errors";
row.publishedH = h;
row.publishedDt = h;
row.finalTime = 4.0;
row.xi = opts.Xi;

match = published(published.h == h, :);
if isempty(match)
    row.status = "skipped";
    row.message = "No published Chen-Ling row for requested h.";
    return;
end
row.N = match.N(1);
row.sqrtN = sqrt(row.N);
row.publishedH1Error = match.H1Error(1);
row.publishedH2Error = match.H2Error(1);

problem = chenLingMovingEllipsoidProblem(row.finalTime);
solverProfile = lower(string(opts.SolverProfile));
fprintf('\nChen-Ling moving ellipsoid: h=%.5g T=%.5g xi=%d N=%d\n', ...
    h, row.finalTime, opts.Xi, row.N);
try
    tic;
    if solverProfile == "minimal"
        result = kp.manifold.runLagrangianMovingADRConvergence( ...
            problem, opts.Xi, row.N, opts.DtScale, ...
            'Mu', 1.0, ...
            'FinalTime', row.finalTime, ...
            'FixedTimeStep', row.publishedDt, ...
            'DiffMatUpdateMethod', 'direct', ...
            'SpectrumCheck', false, ...
            'HyperviscosityUpdateMode', 'none', ...
            'NormalMode', 'exact', ...
            'RecordMassDiagnostics', true, ...
            'TrackTimeErrors', true, ...
            'ExactStartup', opts.ExactStartup, ...
            'GlobalSolveMethod', "direct", ...
            'UseRearrangement', false, ...
            'MassCorrectionMode', "off");
    else
        result = kp.manifold.runLagrangianMovingADRConvergence( ...
            problem, opts.Xi, row.N, opts.DtScale, ...
            'Mu', 1.0, ...
            'FinalTime', row.finalTime, ...
            'FixedTimeStep', row.publishedDt, ...
            'DiffMatUpdateMethod', 'defect', ...
            'DefectTolerance', 1.0e-6, ...
            'MaxDefectIterations', 4, ...
            'NeighborUpdateMode', "periodic", ...
            'NeighborSearchInterval', 5, ...
            'SpectrumCheck', false, ...
            'HyperviscosityUpdateMode', 'adaptiveGeometry', ...
            'HyperviscosityDriftTolerance', 0.05, ...
            'MaxHyperviscositySkippedSteps', 5, ...
            'NormalMode', 'calibratedSbf', ...
            'NormalNeighborCount', 32, ...
            'GlobalSBFNormalDegree', 7, ...
            'GlobalSBFNormalOrder', 8, ...
            'GlobalSBFBalanceSafety', 0.1, ...
            'GlobalSBFControlPointScale', 1 / 3, ...
            'GlobalSBFMinControlPointCount', 48, ...
            'RecordMassDiagnostics', true, ...
            'TrackTimeErrors', true, ...
            'ExactStartup', opts.ExactStartup, ...
            'GlobalSolveMethod', "gmresIluDefect", ...
            'GlobalGMRESTolerance', opts.GlobalGMRESTolerance, ...
            'GlobalGMRESRestart', 40, ...
            'GlobalGMRESMaxIterations', 30, ...
            'GlobalDefectSweeps', 4, ...
            'GlobalILUDropTolerance', 1.0e-4, ...
            'UseRearrangement', true, ...
            'QualityThreshold', 1.55, ...
            'PredictiveLookaheadSteps', 3, ...
            'MinStepsBetweenRearrangements', 8, ...
            'RearrangementTransferMode', "localTp", ...
            'MassCorrectionMode', "balance");
    end
    row.elapsedSeconds = toc;
    row.status = "ok";
    row.h = result.h;
    row.dt = result.dt;
    row.nsteps = result.nsteps;
    row.ourMaxAbsH1SemiError = result.maxAbsH1SemiError;
    row.ourMaxAbsH2SemiError = result.maxAbsH2SemiError;
    row.ourRelativeLInfError = result.finalRelLInfError;
    row.ourRelativeL2Error = result.finalRelL2Error;
    row.maxRelativeL2Error = result.maxRelL2Error;
    row.timeRelativeL2L2Error = result.timeRelL2L2Error;
    row.massRelError = result.massRelError;
    row.balanceRelResidual = result.balanceRelResidual;
    row.normalRMSAngle = result.normalRMSAngle;
    row.normalControlPointCount = result.normalControlPointCount;
    row.numRearrangements = result.numRearrangements;
    row.relativeH1ToPublished = row.ourMaxAbsH1SemiError / row.publishedH1Error;
    row.relativeH2ToPublished = row.ourMaxAbsH2SemiError / row.publishedH2Error;
    row.hOverPublishedH = row.h / h;
    row.dtOverPublishedDt = row.dt / row.publishedDt;
    fprintf('  H1=%.6e published=%.6e ratio=%.3f; H2=%.6e published=%.6e ratio=%.3f\n', ...
        row.ourMaxAbsH1SemiError, row.publishedH1Error, row.relativeH1ToPublished, ...
        row.ourMaxAbsH2SemiError, row.publishedH2Error, row.relativeH2ToPublished);
catch ME
    row.elapsedSeconds = toc;
    row.status = "failed";
    row.message = string(ME.message);
    fprintf('  FAILED: %s\n', ME.message);
    if ~opts.ContinueOnError
        rethrow(ME);
    end
end
end

function problem = petrasLingOscillatingEllipsoidProblem(finalTime)
problem.label = "Petras-Ling oscillating ellipsoid";
problem.title = problem.label;
problem.finalTime = finalTime;
problem.sampleMaterial = @(N, t) sphereMaterial(N);
problem.geometry = @(N, t) geometryFromMaterial(sphereMaterial(N), t);
problem.geometryFromMaterial = @geometryFromMaterial;
problem.backtraceMaterial = @(material, tNow, tPast) material;
problem.exact = @exactSolution;
problem.exactGradient = @(t, material) exactSurfaceGradient(t, material, @axesAt, 6.0);
problem.exactLaplacian = @(t, material) exactSurfaceLaplacian(t, material, @axesAt, 6.0);
problem.forcing = @forcing;
problem.exactMass = @(t, material, geom, mu) sum(geom.weights(:) .* exactSolution(t, material));
end

function problem = chenLingMovingEllipsoidProblem(finalTime)
problem.label = "Chen-Ling moving ellipsoid";
problem.title = problem.label;
problem.finalTime = finalTime;
problem.sampleMaterial = @(N, t) sphereMaterial(N);
problem.geometry = @(N, t) chenLingGeometryFromMaterial(sphereMaterial(N), t);
problem.geometryFromMaterial = @chenLingGeometryFromMaterial;
problem.backtraceMaterial = @(material, tNow, tPast) material;
problem.exact = @chenLingExactSolution;
problem.exactGradient = @(t, material) exactSurfaceGradient(t, material, @chenLingAxesAt, 1.0);
problem.exactLaplacian = @(t, material) exactSurfaceLaplacian(t, material, @chenLingAxesAt, 1.0);
problem.forcing = @chenLingForcing;
problem.exactMass = @(t, material, geom, mu) sum(geom.weights(:) .* chenLingExactSolution(t, material));
end

function material = sphereMaterial(N)
persistent cachedN cachedU
if isempty(cachedN) || cachedN ~= N
    cachedN = N;
    cachedU = kp.manifold.sampleSphereSurfaceFPS(N, @(U) U, ...
        'CandidateFactor', 4, ...
        'RawCandidateFactor', 2, ...
        'UseAreaWeightedCandidates', false);
end
material = struct('U', cachedU);
end

function geom = geometryFromMaterial(material, t)
U = material.U;
axesLengths = axesAt(t);
X = U .* axesLengths;
g = X ./ (axesLengths .^ 2);
normals = g ./ max(vecnorm(g, 2, 2), eps);
geom.X = X;
geom.normals = normals;
geom.weights = ellipsoidSurfaceWeights(U, axesLengths);
geom.area = sum(geom.weights);
geom.h = sqrt(geom.area / size(U, 1));
geom.material = struct('U', U, 'X', X);
end

function geom = chenLingGeometryFromMaterial(material, t)
U = material.U;
axesLengths = chenLingAxesAt(t);
X = U .* axesLengths;
g = X ./ (axesLengths .^ 2);
normals = g ./ max(vecnorm(g, 2, 2), eps);
geom.X = X;
geom.normals = normals;
geom.weights = ellipsoidSurfaceWeights(U, axesLengths);
geom.area = sum(geom.weights);
geom.h = sqrt(geom.area / size(U, 1));
geom.material = struct('U', U, 'X', X);
end

function c = exactSolution(t, material)
U = material.U;
X = U .* axesAt(t);
c = exp(-6 * t) .* X(:, 1) .* X(:, 2);
end

function c = chenLingExactSolution(t, material)
U = material.U;
X = U .* chenLingAxesAt(t);
c = exp(-t) .* X(:, 1) .* X(:, 2);
end

function grad = exactSurfaceGradient(t, material, axesFunction, decayRate)
U = material.U;
axesLengths = axesFunction(t);
X = U .* axesLengths;
g = X ./ (axesLengths .^ 2);
normals = g ./ max(vecnorm(g, 2, 2), eps);
ambient = exp(-decayRate * t) .* [X(:, 2), X(:, 1), zeros(size(X, 1), 1)];
normalPart = sum(ambient .* normals, 2);
grad = ambient - normals .* normalPart;
end

function lap = exactSurfaceLaplacian(t, material, axesFunction, decayRate)
U = material.U;
axesLengths = axesFunction(t);
X = U .* axesLengths;
g = X ./ (axesLengths .^ 2);
normal = g ./ max(vecnorm(g, 2, 2), eps);
H = ellipsoidMeanCurvature(X, axesLengths);
lapX = -2 * H .* normal;
gradDot = -normal(:, 1) .* normal(:, 2);
lapBase = X(:, 2) .* lapX(:, 1) + X(:, 1) .* lapX(:, 2) + 2 * gradDot;
lap = exp(-decayRate * t) .* lapBase;
end

function f = forcing(t, material, mu)
U = material.U;
axesLengths = axesAt(t);
axesPrime = axesPrimeAt(t);
X = U .* axesLengths;
V = U .* axesPrime;
geom = geometryFromMaterial(material, t);
normal = geom.normals;
H = ellipsoidMeanCurvature(X, axesLengths);
lapX = -2 * H .* normal;
base = X(:, 1) .* X(:, 2);
c = exp(-6 * t) .* base;
materialDerivative = exp(-6 * t) .* ( ...
    -6 * base + V(:, 1) .* X(:, 2) + X(:, 1) .* V(:, 2));
divv = linearMapSurfaceDivergence(U, axesLengths, axesPrime);
gradDot = -normal(:, 1) .* normal(:, 2);
lapBase = X(:, 2) .* lapX(:, 1) + X(:, 1) .* lapX(:, 2) + 2 * gradDot;
lapc = exp(-6 * t) .* lapBase;
f = materialDerivative + c .* divv - mu .* lapc;
end

function f = chenLingForcing(t, material, mu)
U = material.U;
axesLengths = chenLingAxesAt(t);
axesPrime = chenLingAxesPrimeAt(t);
X = U .* axesLengths;
V = U .* axesPrime;
geom = chenLingGeometryFromMaterial(material, t);
normal = geom.normals;
H = ellipsoidMeanCurvature(X, axesLengths);
lapX = -2 * H .* normal;
base = X(:, 1) .* X(:, 2);
c = exp(-t) .* base;
materialDerivative = exp(-t) .* ( ...
    -base + V(:, 1) .* X(:, 2) + X(:, 1) .* V(:, 2));
divv = linearMapSurfaceDivergence(U, axesLengths, axesPrime);
gradDot = -normal(:, 1) .* normal(:, 2);
lapBase = X(:, 2) .* lapX(:, 1) + X(:, 1) .* lapX(:, 2) + 2 * gradDot;
lapc = exp(-t) .* lapBase;
f = materialDerivative + c .* divv - mu .* lapc;
end

function axesLengths = axesAt(t)
a = 1 + sin(2 * t);
axesLengths = [sqrt(a), 1, 1];
end

function axesPrime = axesPrimeAt(t)
a = 1 + sin(2 * t);
ap = 2 * cos(2 * t);
axesPrime = [0.5 * ap / sqrt(a), 0, 0];
end

function axesLengths = chenLingAxesAt(t)
axesLengths = [1 + 0.25 * sin(t), 1, 1];
end

function axesPrime = chenLingAxesPrimeAt(t)
axesPrime = [0.25 * cos(t), 0, 0];
end

function divv = linearMapSurfaceDivergence(U, axesLengths, axesPrime)
logDetRate = sum(axesPrime ./ axesLengths);
s2 = sum((U .^ 2) ./ (axesLengths .^ 2), 2);
logStretchRate = -sum((U .^ 2) .* (axesPrime ./ (axesLengths .^ 3)), 2) ./ s2;
divv = logDetRate + logStretchRate;
end

function H = ellipsoidMeanCurvature(X, axesLengths)
A = 1 ./ (axesLengths .^ 2);
g = X .* A;
s = vecnorm(g, 2, 2);
trA = sum(A);
gAg = sum((g .^ 2) .* A, 2);
divn = trA ./ s - gAg ./ (s .^ 3);
H = 0.5 * divn;
end

function weights = ellipsoidSurfaceWeights(U, axesLengths)
baseWeight = 4 * pi / size(U, 1);
jacobian = prod(axesLengths) .* vecnorm(U ./ axesLengths, 2, 2);
weights = baseWeight .* jacobian;
end

function T = petrasLingOscillatingEllipsoidTable()
deltaX = [0.2; 0.2; 0.1; 0.1; 0.05; 0.05; 0.025; 0.025];
finalTime = [0.08; 0.16; 0.08; 0.16; 0.08; 0.16; 0.08; 0.16];
relativeLInfError = [8.89e-2; 1.76e-1; 2.60e-2; 4.91e-2; ...
    6.80e-3; 1.26e-2; 1.70e-3; 3.20e-3];
T = table(deltaX, finalTime, relativeLInfError);
end

function T = chenLingMovingEllipsoidTable()
h = [3/4; 1/2; 1/4; 1/8; 1/16];
N = [36; 78; 312; 1206; 4836];
H2Error = [5.5926e-3; 1.9828e-3; 4.3551e-4; 9.6060e-5; 2.3407e-5];
H1Error = [3.7442e-3; 9.2059e-4; 1.8276e-4; 4.0135e-5; 9.7715e-6];
T = table(h, N, H1Error, H2Error);
end

function T = lingBenchmarkInventory()
sourceKey = [
    "ChenLing2020";
    "ChenLing2020";
    "ChenLing2020";
    "ChenLing2020";
    "PetrasLingPiretRuuth2019";
    "PetrasLingPiretRuuth2019";
    "PetrasLingRuuth2022"];
benchmark = [
    "Evolving curve";
    "Evolving ellipsoid";
    "Expanding sphere mass";
    "Kissing spheres under mean-curvature motion";
    "Expanding circle";
    "Oscillating ellipsoid";
    "Meshfree semi-Lagrangian surface advection"];
publishedMetric = [
    "Linf_t L2 and H2 on a moving curve";
    "Linf_t H1 and H2 on a moving surface";
    "Mass conservation plot";
    "Qualitative high-curvature snapshots and mass curves";
    "Relative Linf on a moving curve";
    "Relative Linf on a moving ellipsoid";
    "Advection/conservation examples on static surfaces"];
useInThisStudy = [
    "not direct: one-dimensional surface";
    "inventory: requires H1/H2 diagnostic for exact metric";
    "inventory: no tabulated error values";
    "inventory: no manufactured scalar-error table";
    "not direct: one-dimensional surface";
    "direct numerical comparison";
    "context for SL backfill, not a moving-surface ADR table"];
T = table(sourceKey, benchmark, publishedMetric, useInThisStudy);
end

function study = buildStudy(rows, published, chenLingPublished, inventory, opts)
study = struct();
study.table = struct2table(rows);
study.petrasLingPublished = published;
study.chenLingPublished = chenLingPublished;
study.inventory = inventory;
study.xi = opts.Xi;
study.comparisonType = "LingPetrasMovingSurface";
study.metric = "relative ell_infty plus H1/H2 seminorm comparisons";
end

function saveStudy(study, outputPrefix)
save(string(outputPrefix) + ".mat", 'study');
writetable(study.table, string(outputPrefix) + ".csv");
writetable(study.petrasLingPublished, string(outputPrefix) + "_petras_ling_published.csv");
writetable(study.chenLingPublished, string(outputPrefix) + "_chen_ling_published.csv");
writetable(study.inventory, string(outputPrefix) + "_inventory.csv");
end

function row = emptyRow()
row = struct( ...
    'sourceKey', "", ...
    'comparisonFamily', "", ...
    'benchmark', "", ...
    'metric', "", ...
    'publishedDeltaX', NaN, ...
    'publishedH', NaN, ...
    'publishedDt', NaN, ...
    'finalTime', NaN, ...
    'publishedError', NaN, ...
    'publishedH1Error', NaN, ...
    'publishedH2Error', NaN, ...
    'xi', NaN, ...
    'N', NaN, ...
    'sqrtN', NaN, ...
    'h', NaN, ...
    'dt', NaN, ...
    'nsteps', NaN, ...
    'ourRelativeLInfError', NaN, ...
    'ourRelativeL2Error', NaN, ...
    'ourMaxAbsH1SemiError', NaN, ...
    'ourMaxAbsH2SemiError', NaN, ...
    'maxRelativeL2Error', NaN, ...
    'timeRelativeL2L2Error', NaN, ...
    'relativeToPublished', NaN, ...
    'relativeH1ToPublished', NaN, ...
    'relativeH2ToPublished', NaN, ...
    'hOverDeltaX', NaN, ...
    'hOverPublishedH', NaN, ...
    'dtOverPublishedDt', NaN, ...
    'massRelError', NaN, ...
    'balanceRelResidual', NaN, ...
    'normalRMSAngle', NaN, ...
    'normalControlPointCount', NaN, ...
    'numRearrangements', NaN, ...
    'elapsedSeconds', NaN, ...
    'status', "", ...
    'message', "");
end
