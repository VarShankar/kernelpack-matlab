function study = moving_surface_adr_tp_dupin_rearrangement(varargin)
%MOVING_SURFACE_ADR_TP_DUPIN_REARRANGEMENT Dupin cyclide rearrangement test.
%   Compares a pure Lagrangian marker evolution against occasional
%   SBF-based rearrangement on a fixed Dupin cyclide with sliding labels.
%   The default problem is forced conservative advection,
%       D_t c + c div_Gamma(v) = f,
%   with no diffusion or reaction.  The prescribed tangential flow is
%   analytic, but the default numerical solve uses RBF-FD surface divergence.
%   Exact parameter derivatives are used only to evaluate the forcing and the
%   reference solution.

parser = inputParser();
parser.addParameter('N', 512, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('Xi', 2, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('FinalTime', 0.12, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('DtScale', 0.004, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('Nu', 0.0, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('QualityThreshold', 1.75, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('MinStepsBetweenRearrangements', 20, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('FlowScale', 0.35, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('DivergenceMode', "rbffd", @(x) isstring(x) || ischar(x));
parser.addParameter('SBFControlPointCount', Inf, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('UseExactStartup', true, @(x) islogical(x) && isscalar(x));
parser.addParameter('RecordOperatorDiagnostics', false, @(x) islogical(x) && isscalar(x));
parser.addParameter('UseRearrangement', [false, true], @(x) islogical(x) || isnumeric(x));
parser.addParameter('UseParallel', false, @(x) islogical(x) && isscalar(x));
parser.addParameter('OutputPrefix', fullfile(pwd, ...
    'moving_surface_adr_tp_dupin_rearrangement'));
parser.addParameter('ImagePath', fullfile(pwd, 'docs', 'figures', ...
    'moving_surface_adr_tp_dupin_rearrangement.png'));
parser.addParameter('SaveOutputs', true, @(x) islogical(x) && isscalar(x));
parser.parse(varargin{:});
opts = parser.Results;

N = round(opts.N);
xi = round(opts.Xi);
theta = 2;
h = sqrt(dupinAreaApprox() / N);
dtTarget = opts.DtScale * h^(xi / 3);
nsteps = max(3, ceil(opts.FinalTime / dtTarget));
dt = opts.FinalTime / nsteps;

modes = logical(opts.UseRearrangement(:).');
results = repmat(emptyRunResult(), 1, numel(modes));
for k = 1:numel(modes)
    results(k) = runOneMode(N, xi, theta, opts.Nu, opts.FinalTime, dt, nsteps, ...
        modes(k), opts.QualityThreshold, round(opts.MinStepsBetweenRearrangements), ...
        opts.FlowScale, string(opts.DivergenceMode), opts.SBFControlPointCount, ...
        opts.UseExactStartup, opts.RecordOperatorDiagnostics, opts.UseParallel);
end

study = struct();
study.N = N;
study.xi = xi;
study.h = h;
study.dt = dt;
study.nsteps = nsteps;
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

function result = runOneMode(N, xi, theta, nu, T, dt, nsteps, useRearrangement, ...
    qualityThreshold, minStepsBetweenRearrangements, flowScale, divergenceMode, ...
    sbfControlPointCount, useExactStartup, recordOperatorDiagnostics, useParallel)
flow = dupinFlowParameters(flowScale);
updater = kp.manifold.TangentPlaneDiffMatUpdater(xi, ...
    'Theta', theta, ...
    'DefectTolerance', 1.0e-6, ...
    'MaxDefectIterations', 4, ...
    'NeighborUpdateMode', "motionbound", ...
    'UseParallel', useParallel);

[theta0, phi0] = dupinQuasiUniformMaterial(N, []);
t0 = 0.0;
geom0 = dupinGeometry(theta0, phi0);
c0 = exactConcentration(t0, theta0, phi0);
x0 = geom0.X;
operatorDiagnostics = computeOperatorDiagnostics( ...
    recordOperatorDiagnostics, xi, theta, useParallel, geom0, c0, t0, theta0, phi0, flow);

[theta1, phi1] = advanceLabels(theta0, phi0, dt, flow);
t1 = dt;
geom1 = dupinGeometry(theta1, phi1);
x1 = geom1.X;
if useExactStartup
    c1 = exactConcentration(t1, theta1, phi1);
    stats1 = emptyOperatorStats();
else
    v1 = (x1 - x0) / dt;
    [c1, stats1] = solveBDFStep(1, updater, geom1, v1, t1, theta1, phi1, ...
        {c0}, nu, dt, flow, divergenceMode);
end

[theta2, phi2] = advanceLabels(theta1, phi1, dt, flow);
t2 = 2 * dt;
geom2 = dupinGeometry(theta2, phi2);
x2 = geom2.X;
if useExactStartup
    c2 = exactConcentration(t2, theta2, phi2);
    stats2 = emptyOperatorStats();
else
    v2 = (3 * x2 - 4 * x1 + x0) / (2 * dt);
    [c2, stats2] = solveBDFStep(2, updater, geom2, v2, t2, theta2, phi2, ...
        {c0, c1}, nu, dt, flow, divergenceMode);
end

qualityTimes = zeros(nsteps + 1, 1);
qualityValues = zeros(nsteps + 1, 1);
qualityTimes(1:3) = [t0; t1; t2];
qualityValues(1:3) = [pointQuality(x0); pointQuality(x1); pointQuality(x2)];
rearrangementTimes = zeros(0, 1);
sbfControlPointCounts = zeros(0, 1);
sbfFillDistances = zeros(0, 1);
operatorStats = [stats1; stats2];
lastRearrangementStep = -inf;

for step = 3:nsteps
    t3 = step * dt;
    [theta3, phi3] = advanceLabels(theta2, phi2, dt, flow);
    geom3 = dupinGeometry(theta3, phi3);
    x3 = geom3.X;
    v3 = (11 * x3 - 18 * x2 + 9 * x1 - 2 * x0) / (6 * dt);
    [c3, stats3] = solveBDFStep(3, updater, geom3, v3, t3, theta3, phi3, ...
        {c0, c1, c2}, nu, dt, flow, divergenceMode);
    operatorStats = [operatorStats; stats3]; %#ok<AGROW>

    theta0 = theta1; phi0 = phi1; x0 = x1; c0 = c1;
    theta1 = theta2; phi1 = phi2; x1 = x2; c1 = c2;
    theta2 = theta3; phi2 = phi3; x2 = x3; c2 = c3;

    q = pointQuality(x2);
    qualityTimes(step + 1) = t3;
    qualityValues(step + 1) = q;
    canRearrange = useRearrangement && step >= 3 && ...
        step - lastRearrangementStep >= minStepsBetweenRearrangements && ...
        q > qualityThreshold;
    if canRearrange
        [theta0, phi0, x0, c0, theta1, phi1, x1, c1, theta2, phi2, x2, c2] = ...
            rearrangeHistory(theta0, phi0, c0, theta1, phi1, c1, ...
            theta2, phi2, x2, c2, dt, flow, sbfControlPointCount); %#ok<ASGLU>
        [~, sbfInfo] = kp.manifold.selectSBFControlPoints(dupinGeometry(theta2, phi2), ...
            'Count', normalizeControlPointCount(sbfControlPointCount, numel(theta2)));
        updater.reset();
        lastRearrangementStep = step;
        rearrangementTimes(end + 1, 1) = t3; %#ok<AGROW>
        sbfControlPointCounts(end + 1, 1) = sbfInfo.controlPointCount; %#ok<AGROW>
        sbfFillDistances(end + 1, 1) = sbfInfo.fillDistance; %#ok<AGROW>
        qualityValues(step + 1) = pointQuality(x2);
    end
end

cex = exactConcentration(T, theta2, phi2);
result = emptyRunResult();
if useRearrangement
    result.mode = "rearranged";
else
    result.mode = "lagrangian";
end
result.relerr = norm(c2 - cex) / max(norm(cex), 1.0e-14);
result.numRearrangements = numel(rearrangementTimes);
result.initialQuality = qualityValues(1);
result.finalQuality = qualityValues(nsteps + 1);
result.maxQuality = max(qualityValues);
result.rearrangementTimes = {rearrangementTimes};
if isempty(sbfControlPointCounts)
    result.sbfControlPointCount = NaN;
    result.sbfControlPointFraction = NaN;
    result.sbfFillDistance = NaN;
else
    result.sbfControlPointCount = max(sbfControlPointCounts);
    result.sbfControlPointFraction = result.sbfControlPointCount / N;
    result.sbfFillDistance = max(sbfFillDistances);
end
result.qualityTimes = {qualityTimes};
result.qualityValues = {qualityValues};
result.finalTheta = {theta2};
result.finalPhi = {phi2};
result.finalX = {x2};
result.finalNumerical = {c2};
result.finalExact = {cex};
result.operatorStats = {operatorStats};
result.initialLaplaceRelError = operatorDiagnostics.laplaceRelError;
result.initialDivergenceRelError = operatorDiagnostics.divergenceRelError;
end

function diagnostics = computeOperatorDiagnostics(recordDiagnostics, xi, theta, ...
    useParallel, geom, c, t, thetaNodes, phiNodes, flow)
diagnostics = struct('laplaceRelError', NaN, 'divergenceRelError', NaN);
if ~recordDiagnostics
    return;
end

diagnosticUpdater = kp.manifold.TangentPlaneDiffMatUpdater(xi, ...
    'Theta', theta, ...
    'UseParallel', useParallel);
[L, Gx, Gy, Gz] = diagnosticUpdater.assemble(geom.X, geom.normals);
lapExact = parametricSurfaceLaplacian(t, thetaNodes, phiNodes);
diagnostics.laplaceRelError = norm(L * c - lapExact) / max(norm(lapExact), eps);

[thetaDot, phiDot] = labelVelocity(thetaNodes, phiNodes, flow);
velocity = parameterVelocity(thetaNodes, phiNodes, thetaDot, phiDot);
divExact = parametricSurfaceDivergence(thetaNodes, phiNodes, flow);
diagnostics.divergenceRelError = norm(surfaceDivergence(Gx, Gy, Gz, velocity) - divExact) / ...
    max(norm(divExact), eps);
end

function [c, stats] = solveBDFStep(order, updater, geom, velocity, t, theta, phi, ...
    history, nu, dt, flow, divergenceMode)
needsDiffusion = nu ~= 0;
needsRbfDivergence = lower(string(divergenceMode)) == "rbffd";
if needsDiffusion || needsRbfDivergence
    [L, Gx, Gy, Gz, stats] = updater.assemble(geom.X, geom.normals);
else
    L = [];
    Gx = [];
    Gy = [];
    Gz = [];
    stats = emptyOperatorStats();
end

switch lower(string(divergenceMode))
    case "rbffd"
        divv = surfaceDivergence(Gx, Gy, Gz, velocity);
    case "parametric"
        divv = parametricSurfaceDivergence(theta, phi, flow);
    otherwise
        error('kp:examples:BadDivergenceMode', ...
            'Unknown divergence mode "%s".', divergenceMode);
end
f = manufacturedForcing(t, theta, phi, nu, flow);
I = speye(size(geom.X, 1));

switch order
    case 1
        c0 = history{1};
        lhs = diffusionLhs(I, L, nu * dt);
        rhs = c0 + dt * f - dt * c0 .* divv;
    case 2
        c0 = history{1};
        c1 = history{2};
        beta = 2 / 3;
        extrap = 2 * c1 - c0;
        lhs = diffusionLhs(I, L, beta * nu * dt);
        rhs = (4 / 3) * c1 - (1 / 3) * c0 + beta * dt * f ...
            - beta * dt * extrap .* divv;
    case 3
        c0 = history{1};
        c1 = history{2};
        c2 = history{3};
        beta = 6 / 11;
        extrap = 3 * c2 - 3 * c1 + c0;
        lhs = diffusionLhs(I, L, beta * nu * dt);
        rhs = (18 / 11) * c2 - (9 / 11) * c1 + (2 / 11) * c0 ...
            + beta * dt * f - beta * dt * extrap .* divv;
    otherwise
        error('kp:examples:BadBDFOrder', 'Unsupported BDF order %d.', order);
end
if needsDiffusion
    c = lhs \ rhs;
else
    c = rhs;
end
end

function lhs = diffusionLhs(I, L, scale)
if scale == 0
    lhs = I;
else
    lhs = I - scale * L;
end
end

function [theta0New, phi0New, x0New, c0New, theta1New, phi1New, x1New, c1New, ...
    theta2New, phi2New, x2New, c2New] = rearrangeHistory( ...
    theta0, phi0, c0, theta1, phi1, c1, theta2, phi2, x2, c2, dt, flow, ...
    sbfControlPointCount)
geomCurrent = dupinGeometry(theta2, phi2);
geomCurrent.X = x2;
model = kp.manifold.fitGlobalParametricSBFGeometry(geomCurrent, ...
    'Degree', 7, ...
    'ControlPointCount', sbfControlPointCount);
[theta2New, phi2New] = dupinQuasiUniformMaterial(numel(theta2), model);

[theta1New, phi1New] = advanceLabels(theta2New, phi2New, -dt, flow);
[theta0New, phi0New] = advanceLabels(theta2New, phi2New, -2 * dt, flow);

c2New = periodicScalarInterpolate(theta2, phi2, c2, theta2New, phi2New, model.controlIds);
c1New = periodicScalarInterpolate(theta1, phi1, c1, theta1New, phi1New, model.controlIds);
c0New = periodicScalarInterpolate(theta0, phi0, c0, theta0New, phi0New, model.controlIds);

x2New = kp.manifold.evaluateGlobalParametricSBFGeometry(model, theta2New, phi2New);
x1New = kp.manifold.evaluateGlobalParametricSBFGeometry(model, theta1New, phi1New);
x0New = kp.manifold.evaluateGlobalParametricSBFGeometry(model, theta0New, phi0New);
end

function f = manufacturedForcing(t, theta, phi, nu, flow)
[thetaDot, phiDot] = labelVelocity(theta, phi, flow);
c = exactConcentration(t, theta, phi);
[ct, ctheta, cphi] = exactDerivatives(t, theta, phi);
materialDerivative = ct + thetaDot .* ctheta + phiDot .* cphi;
divv = parametricSurfaceDivergence(theta, phi, flow);
if nu == 0
    lapc = 0;
else
    lapc = parametricSurfaceLaplacian(t, theta, phi);
end
f = materialDerivative + c .* divv - nu .* lapc;
end

function c = exactConcentration(t, theta, phi)
base = concentrationBase(theta, phi);
c = exp(-t) .* base;
end

function [ct, ctheta, cphi, cthetaTheta, cthetaPhi, cphiPhi] = exactDerivatives(t, theta, phi)
base = concentrationBase(theta, phi);
baseTheta = -0.20 * sin(theta) - 0.10 * sin(theta - 2 * phi);
basePhi = 0.15 * cos(phi) + 0.20 * sin(theta - 2 * phi);
baseThetaTheta = -0.20 * cos(theta) - 0.10 * cos(theta - 2 * phi);
baseThetaPhi = 0.20 * cos(theta - 2 * phi);
basePhiPhi = -0.15 * sin(phi) - 0.40 * cos(theta - 2 * phi);
decay = exp(-t);
ct = -exp(-t) .* base;
ctheta = decay .* baseTheta;
cphi = decay .* basePhi;
cthetaTheta = decay .* baseThetaTheta;
cthetaPhi = decay .* baseThetaPhi;
cphiPhi = decay .* basePhiPhi;
end

function base = concentrationBase(theta, phi)
base = 1.0 + 0.20 * cos(theta) + 0.15 * sin(phi) + ...
    0.10 * cos(theta - 2 * phi);
end

function lap = parametricSurfaceLaplacian(t, theta, phi)
[E, F, G, sqrtg, ~, Ephi, Ftheta, Fphi, Gtheta, ~, detg, ...
    detTheta, detPhi] = surfaceMetric(theta, phi);
g11 = G ./ detg;
g12 = -F ./ detg;
g22 = E ./ detg;
sqrtTheta = detTheta ./ (2 * sqrtg);
sqrtPhi = detPhi ./ (2 * sqrtg);

a11 = sqrtg .* g11;
a12 = sqrtg .* g12;
a22 = sqrtg .* g22;
a11Theta = (Gtheta .* sqrtg - G .* sqrtTheta) ./ (sqrtg .^ 2);
a12Theta = -(Ftheta .* sqrtg - F .* sqrtTheta) ./ (sqrtg .^ 2);
a12Phi = -(Fphi .* sqrtg - F .* sqrtPhi) ./ (sqrtg .^ 2);
a22Phi = (Ephi .* sqrtg - E .* sqrtPhi) ./ (sqrtg .^ 2);

[~, ctheta, cphi, cthetaTheta, cthetaPhi, cphiPhi] = exactDerivatives(t, theta, phi);
lap = (a11Theta .* ctheta + a11 .* cthetaTheta + ...
    a12Theta .* cphi + a12 .* cthetaPhi + ...
    a12Phi .* ctheta + a12 .* cthetaPhi + ...
    a22Phi .* cphi + a22 .* cphiPhi) ./ sqrtg;
end

function divv = parametricSurfaceDivergence(theta, phi, flow)
[~, ~, ~, ~, ~, ~, ~, ~, ~, ~, detg, detTheta, detPhi] = ...
    surfaceMetric(theta, phi);
[thetaDot, phiDot] = labelVelocity(theta, phi, flow);
thetaDotTheta = 2 * flow.A * cos(2 * theta) .* cos(phi);
phiDotPhi = 3 * flow.C * cos(theta) .* cos(3 * phi);
logSqrtTheta = detTheta ./ (2 * detg);
logSqrtPhi = detPhi ./ (2 * detg);
divv = thetaDotTheta + phiDotPhi + ...
    logSqrtTheta .* thetaDot + logSqrtPhi .* phiDot;
end

function [thetaDot, phiDot] = labelVelocity(theta, phi, flow)
thetaDot = flow.A * sin(2 * theta) .* cos(phi);
phiDot = flow.B + flow.C * cos(theta) .* sin(3 * phi);
end

function velocity = parameterVelocity(theta, phi, thetaDot, phiDot)
[~, xtheta, xphi] = dupinDifferentialGeometry(theta, phi);
velocity = xtheta .* thetaDot(:) + xphi .* phiDot(:);
end

function [theta, phi] = advanceLabels(theta, phi, dt, flow)
maxStep = 2.5e-3;
nsub = max(1, ceil(abs(dt) / maxStep));
subdt = dt / nsub;
for k = 1:nsub
    [k1t, k1p] = labelVelocity(theta, phi, flow);
    [k2t, k2p] = labelVelocity(theta + 0.5 * subdt * k1t, ...
        phi + 0.5 * subdt * k1p, flow);
    [k3t, k3p] = labelVelocity(theta + 0.5 * subdt * k2t, ...
        phi + 0.5 * subdt * k2p, flow);
    [k4t, k4p] = labelVelocity(theta + subdt * k3t, ...
        phi + subdt * k3p, flow);
    theta = theta + (subdt / 6) * (k1t + 2 * k2t + 2 * k3t + k4t);
    phi = phi + (subdt / 6) * (k1p + 2 * k2p + 2 * k3p + k4p);
end
theta = mod(theta, 2 * pi);
phi = mod(phi, 2 * pi);
end

function geom = dupinGeometry(theta, phi)
geom = struct();
geom.X = dupinPosition(theta, phi);
geom.normals = dupinNormals(theta, phi);
geom.h = sqrt(dupinAreaApprox() / numel(theta));
geom.material = struct('theta', theta(:), 'phi', phi(:));
end

function X = dupinPosition(theta, phi)
[X, ~, ~] = dupinDifferentialGeometry(theta, phi);
end

function [X, xtheta, xphi, xthetaTheta, xthetaPhi, xphiPhi] = ...
    dupinDifferentialGeometry(theta, phi)
p = dupinParameters();
theta = theta(:);
phi = phi(:);
ct = cos(theta);
st = sin(theta);
cp = cos(phi);
sp = sin(phi);
den = p.a - p.c * ct .* cp;
n = [ ...
    p.d * (p.c - p.a * ct .* cp) + p.b^2 * ct, ...
    p.b * st .* (p.a - p.d * cp), ...
    p.b * sp .* (p.c * ct - p.d)];
nTheta = [ ...
    st .* (p.d * p.a * cp - p.b^2), ...
    p.b * ct .* (p.a - p.d * cp), ...
    -p.b * p.c * sp .* st];
nPhi = [ ...
    p.d * p.a * ct .* sp, ...
    p.b * p.d * st .* sp, ...
    p.b * cp .* (p.c * ct - p.d)];
nThetaTheta = [ ...
    ct .* (p.d * p.a * cp - p.b^2), ...
    -p.b * st .* (p.a - p.d * cp), ...
    -p.b * p.c * sp .* ct];
nThetaPhi = [ ...
    -p.d * p.a * st .* sp, ...
    p.b * p.d * ct .* sp, ...
    -p.b * p.c * cp .* st];
nPhiPhi = [ ...
    p.d * p.a * ct .* cp, ...
    p.b * p.d * st .* cp, ...
    -p.b * sp .* (p.c * ct - p.d)];

denTheta = p.c * st .* cp;
denPhi = p.c * ct .* sp;
denThetaTheta = p.c * ct .* cp;
denThetaPhi = -p.c * st .* sp;
denPhiPhi = p.c * ct .* cp;

[X, xtheta, xphi, xthetaTheta, xthetaPhi, xphiPhi] = quotientDerivatives( ...
    n, nTheta, nPhi, nThetaTheta, nThetaPhi, nPhiPhi, ...
    den, denTheta, denPhi, denThetaTheta, denThetaPhi, denPhiPhi);
end

function [x, xTheta, xPhi, xThetaTheta, xThetaPhi, xPhiPhi] = quotientDerivatives( ...
    n, nTheta, nPhi, nThetaTheta, nThetaPhi, nPhiPhi, ...
    q, qTheta, qPhi, qThetaTheta, qThetaPhi, qPhiPhi)
q2 = q .^ 2;
q3 = q .^ 3;
x = n ./ q;
xTheta = nTheta ./ q - n .* qTheta ./ q2;
xPhi = nPhi ./ q - n .* qPhi ./ q2;
xThetaTheta = nThetaTheta ./ q - ...
    (2 * nTheta .* qTheta + n .* qThetaTheta) ./ q2 + ...
    2 * n .* qTheta .^ 2 ./ q3;
xThetaPhi = nThetaPhi ./ q - ...
    (nTheta .* qPhi + nPhi .* qTheta + n .* qThetaPhi) ./ q2 + ...
    2 * n .* qTheta .* qPhi ./ q3;
xPhiPhi = nPhiPhi ./ q - ...
    (2 * nPhi .* qPhi + n .* qPhiPhi) ./ q2 + ...
    2 * n .* qPhi .^ 2 ./ q3;
end

function normals = dupinNormals(theta, phi)
[~, xtheta, xphi] = dupinDifferentialGeometry(theta, phi);
normals = kp.geometry.normalizeRows(cross(xtheta, xphi, 2));
X = dupinPosition(theta, phi);
center = mean(dupinPosition(linspace(0, 2 * pi, 64).', zeros(64, 1)), 1);
flip = sum((X - center) .* normals, 2) < 0;
normals(flip, :) = -normals(flip, :);
end

function [E, F, G, sqrtg] = firstFundamentalForm(theta, phi)
[E, F, G, sqrtg] = surfaceMetric(theta, phi);
end

function [E, F, G, sqrtg, Etheta, Ephi, Ftheta, Fphi, Gtheta, Gphi, detg, ...
    detTheta, detPhi] = surfaceMetric(theta, phi)
[~, xtheta, xphi, xthetaTheta, xthetaPhi, xphiPhi] = ...
    dupinDifferentialGeometry(theta, phi);
E = sum(xtheta .* xtheta, 2);
F = sum(xtheta .* xphi, 2);
G = sum(xphi .* xphi, 2);
Etheta = 2 * sum(xtheta .* xthetaTheta, 2);
Ephi = 2 * sum(xtheta .* xthetaPhi, 2);
Ftheta = sum(xthetaTheta .* xphi, 2) + sum(xtheta .* xthetaPhi, 2);
Fphi = sum(xthetaPhi .* xphi, 2) + sum(xtheta .* xphiPhi, 2);
Gtheta = 2 * sum(xphi .* xthetaPhi, 2);
Gphi = 2 * sum(xphi .* xphiPhi, 2);
detg = E .* G - F .^ 2;
detg = max(detg, realmin);
detTheta = Etheta .* G + E .* Gtheta - 2 * F .* Ftheta;
detPhi = Ephi .* G + E .* Gphi - 2 * F .* Fphi;
sqrtg = sqrt(detg);
end

function sqrtg = surfaceJacobian(theta, phi)
[~, ~, ~, sqrtg] = firstFundamentalForm(theta, phi);
end

function [theta, phi] = dupinQuasiUniformMaterial(N, model)
numCandidates = max(8 * N, 1600);
m = ceil(sqrt(numCandidates));
vals = 2 * pi * ((0:m-1).' + 0.5) / m;
[Theta, Phi] = ndgrid(vals, vals);
thetaCandidates = Theta(:);
phiCandidates = Phi(:);
if isempty(model)
    Xcand = dupinPosition(thetaCandidates, phiCandidates);
else
    Xcand = kp.manifold.evaluateGlobalParametricSBFGeometry(model, ...
        thetaCandidates, phiCandidates);
end
ids = farthestPointSubset(Xcand, N);
theta = thetaCandidates(ids);
phi = phiCandidates(ids);
end

function ids = farthestPointSubset(X, N)
numCandidates = size(X, 1);
ids = zeros(N, 1);
[~, ids(1)] = max(sum((X - mean(X, 1)) .^ 2, 2));
dist2 = inf(numCandidates, 1);
for k = 2:N
    diff = X - X(ids(k - 1), :);
    dist2 = min(dist2, sum(diff .^ 2, 2));
    [~, ids(k)] = max(dist2);
end
end

function q = pointQuality(X)
tree = KDTreeSearcher(X);
[~, d] = knnsearch(tree, X, 'K', 2);
nearest = d(:, 2);
q = max(nearest) / max(min(nearest), eps);
end

function values = periodicScalarInterpolate(theta, phi, data, thetaq, phiq, controlIds)
degree = 7;
theta = theta(:);
phi = phi(:);
data = data(:);
if nargin < 6 || isempty(controlIds)
    controlIds = (1:numel(theta)).';
else
    controlIds = controlIds(:);
end
thetaControl = theta(controlIds);
phiControl = phi(controlIds);
dataControl = data(controlIds);
r = periodicChordDistance(thetaControl, phiControl, thetaControl, phiControl);
kernel = kp.geometry.phsKernel(r, degree);
reg = 1.0e-12 * max(1.0, max(abs(kernel), [], 'all'));
coefficients = (kernel + reg * eye(numel(thetaControl))) \ dataControl;
rq = periodicChordDistance(thetaq(:), phiq(:), thetaControl, phiControl);
values = kp.geometry.phsKernel(rq, degree) * coefficients;
end

function count = normalizeControlPointCount(count, n)
if ~isfinite(count)
    count = n;
else
    count = min(n, max(1, round(count)));
end
end

function r = periodicChordDistance(thetaA, phiA, thetaB, phiB)
dtheta = thetaA(:) - thetaB(:).';
dphi = phiA(:) - phiB(:).';
r2 = 2 - 2 * cos(dtheta) + 2 - 2 * cos(dphi);
r = sqrt(max(r2, 0));
end

function divv = surfaceDivergence(Gx, Gy, Gz, velocity)
divv = Gx * velocity(:, 1) + Gy * velocity(:, 2) + Gz * velocity(:, 3);
end

function area = dupinAreaApprox()
persistent cachedArea
if isempty(cachedArea)
    m = 220;
    vals = 2 * pi * ((0:m-1).' + 0.5) / m;
    [Theta, Phi] = ndgrid(vals, vals);
    cachedArea = sum(surfaceJacobian(Theta(:), Phi(:))) * (2 * pi / m)^2;
end
area = cachedArea;
end

function p = dupinParameters()
p = struct();
p.a = 1.20;
p.c = 0.72;
p.b = sqrt(p.a^2 - p.c^2);
p.d = 0.34;
end

function flow = dupinFlowParameters(scale)
flow = struct('A', scale * 1.75, 'B', scale * 0.80, 'C', scale * 0.65);
end

function result = emptyRunResult()
result = struct( ...
    'mode', "", ...
    'relerr', NaN, ...
    'numRearrangements', 0, ...
    'initialQuality', NaN, ...
    'finalQuality', NaN, ...
    'maxQuality', NaN, ...
    'sbfControlPointCount', NaN, ...
    'sbfControlPointFraction', NaN, ...
    'sbfFillDistance', NaN, ...
    'rearrangementTimes', {cell(1, 1)}, ...
    'qualityTimes', {cell(1, 1)}, ...
    'qualityValues', {cell(1, 1)}, ...
    'finalTheta', {cell(1, 1)}, ...
    'finalPhi', {cell(1, 1)}, ...
    'finalX', {cell(1, 1)}, ...
    'finalNumerical', {cell(1, 1)}, ...
    'finalExact', {cell(1, 1)}, ...
    'operatorStats', {cell(1, 1)}, ...
    'initialLaplaceRelError', NaN, ...
    'initialDivergenceRelError', NaN);
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
fig = figure('Color', 'w', 'Position', [100, 100, 1360, 520]);
tiledlayout(fig, 1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
sgtitle(fig, sprintf('Forced advection with rearrangement on a sliding Dupin cyclide (N=%d, \\xi=%d)', ...
    study.N, study.xi));

nexttile;
hold on;
for k = 1:numel(study.results)
    plot(study.results(k).qualityTimes{1}, study.results(k).qualityValues{1}, ...
        'LineWidth', 1.6, 'DisplayName', study.results(k).mode);
    rt = study.results(k).rearrangementTimes{1};
    if ~isempty(rt)
        yl = ylim;
        for j = 1:numel(rt)
            plot([rt(j), rt(j)], yl, ':', 'Color', [0.3, 0.3, 0.3], ...
                'HandleVisibility', 'off');
        end
    end
end
hold off;
grid on;
xlabel('time');
ylabel('nearest-neighbor quality');
legend('Location', 'northwest');
title('Marker quality');

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
