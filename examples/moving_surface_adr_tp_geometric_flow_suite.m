function allResults = moving_surface_adr_tp_geometric_flow_suite(caseNames, xiVals, Nvals, dtScale, varargin)
%MOVING_SURFACE_ADR_TP_GEOMETRIC_FLOW_SUITE Analytic moving-surface ADR tests.
%   These are prescribed geometric flows used to validate the Lagrangian
%   moving-surface ADR discretization. The solver does not evolve the
%   geometry law; it solves ADR on analytically prescribed moving nodes.

if nargin < 1 || isempty(caseNames)
    caseNames = ["mcf_sphere", "imcf_sphere", "rotating_breathing_sphere", ...
        "anisotropic_ellipsoid", "breathing_torus"];
end
if nargin < 2
    xiVals = [];
end
if nargin < 3
    Nvals = [];
end
if nargin < 4 || isempty(dtScale)
    dtScale = 0.05;
end

parser = inputParser();
parser.addParameter('DiffMatUpdateMethod', "direct");
parser.addParameter('DefectTolerance', 1.0e-4);
parser.addParameter('MaxDefectIterations', 4);
parser.addParameter('NeighborUpdateMode', "periodic");
parser.addParameter('NeighborSearchInterval', 5);
parser.addParameter('SpectrumCheck', true);
parser.addParameter('HyperviscosityPower', NaN, ...
    @(x) isnumeric(x) && isscalar(x) && (isnan(x) || (isfinite(x) && x >= 1 && x == round(x))));
parser.addParameter('HyperviscosityUpdateMode', "everyStep");
parser.addParameter('HyperviscosityDriftTolerance', 0.05);
parser.addParameter('MaxHyperviscositySkippedSteps', 5);
parser.addParameter('NormalMode', "exact");
parser.addParameter('NormalNeighborCount', 32);
parser.addParameter('GlobalSBFNormalDegree', 7);
parser.addParameter('GlobalSBFControlPointCount', Inf);
parser.addParameter('GlobalSBFNormalOrder', NaN);
parser.addParameter('GlobalSBFBalanceSafety', 0.1);
parser.addParameter('GlobalSBFControlPointScale', 1 / 3);
parser.addParameter('GlobalSBFMinControlPointCount', 48);
parser.addParameter('RecordMassDiagnostics', true);
parser.addParameter('MassCorrectionMode', "balance", @(x) isstring(x) || ischar(x));
parser.addParameter('MassCorrectionMaxRelativeCorrection', Inf, ...
    @(x) isnumeric(x) && isscalar(x) && x > 0);
parser.addParameter('ReturnSolution', false);
parser.addParameter('ExactStartup', true, @(x) islogical(x) && isscalar(x));
parser.addParameter('FinalTime', []);
parser.addParameter('GlobalSolveMethod', "gmresIluDefect");
parser.addParameter('GlobalGMRESTolerance', 1.0e-12);
parser.addParameter('GlobalGMRESRestart', 40);
parser.addParameter('GlobalGMRESMaxIterations', 30);
parser.addParameter('GlobalDefectSweeps', 4);
parser.addParameter('GlobalILUDropTolerance', 1.0e-4);
parser.addParameter('GlobalILURefreshInterval', Inf);
parser.addParameter('GlobalILURefreshIterationThreshold', 30);
parser.addParameter('UseRearrangement', true, @(x) islogical(x) && isscalar(x));
parser.addParameter('QualityThreshold', 1.75, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('PredictiveLookaheadSteps', 3, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('MinStepsBetweenRearrangements', 8, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('RearrangementTransferMode', "localTp", @(x) isstring(x) || ischar(x));
parser.addParameter('WriteOutputs', true);
parser.addParameter('WriteFigures', true, @(x) islogical(x) && isscalar(x));
parser.parse(varargin{:});

caseNames = string(caseNames);
allResults = repmat(struct('caseName', "", 'results', []), 1, numel(caseNames));

for k = 1:numel(caseNames)
    problem = localProblem(caseNames(k));
    if ~isempty(parser.Results.FinalTime)
        problem.finalTime = parser.Results.FinalTime;
    end
    fprintf('\n%s\n', repmat('=', 1, strlength(problem.label)));
    fprintf('%s\n', problem.label);
    fprintf('%s\n\n', repmat('=', 1, strlength(problem.label)));

    if isempty(Nvals)
        caseNvals = problem.defaultNVals;
    else
        caseNvals = Nvals;
    end
    if isempty(xiVals)
        caseXiVals = problem.defaultXiVals;
    else
        caseXiVals = xiVals;
    end

    if parser.Results.WriteOutputs
        if parser.Results.WriteFigures
            imgPath = fullfile(pwd, 'docs', 'figures', problem.imageName);
        else
            imgPath = '';
        end
        matPath = fullfile(pwd, problem.resultsName);
    else
        imgPath = '';
        matPath = '';
    end
    results = kp.manifold.runLagrangianMovingADRConvergence( ...
        problem, caseXiVals, caseNvals, dtScale, ...
        'Mu', 0.1, ...
        'FinalTime', problem.finalTime, ...
        'HyperviscosityPower', parser.Results.HyperviscosityPower, ...
        'DiffMatUpdateMethod', parser.Results.DiffMatUpdateMethod, ...
        'DefectTolerance', parser.Results.DefectTolerance, ...
        'MaxDefectIterations', parser.Results.MaxDefectIterations, ...
        'NeighborUpdateMode', parser.Results.NeighborUpdateMode, ...
        'NeighborSearchInterval', parser.Results.NeighborSearchInterval, ...
        'SpectrumCheck', parser.Results.SpectrumCheck, ...
        'HyperviscosityUpdateMode', parser.Results.HyperviscosityUpdateMode, ...
        'HyperviscosityDriftTolerance', parser.Results.HyperviscosityDriftTolerance, ...
        'MaxHyperviscositySkippedSteps', parser.Results.MaxHyperviscositySkippedSteps, ...
        'NormalMode', parser.Results.NormalMode, ...
        'NormalNeighborCount', parser.Results.NormalNeighborCount, ...
        'GlobalSBFNormalDegree', parser.Results.GlobalSBFNormalDegree, ...
        'GlobalSBFControlPointCount', parser.Results.GlobalSBFControlPointCount, ...
        'GlobalSBFNormalOrder', parser.Results.GlobalSBFNormalOrder, ...
        'GlobalSBFBalanceSafety', parser.Results.GlobalSBFBalanceSafety, ...
        'GlobalSBFControlPointScale', parser.Results.GlobalSBFControlPointScale, ...
        'GlobalSBFMinControlPointCount', parser.Results.GlobalSBFMinControlPointCount, ...
        'RecordMassDiagnostics', parser.Results.RecordMassDiagnostics, ...
        'MassCorrectionMode', parser.Results.MassCorrectionMode, ...
        'MassCorrectionMaxRelativeCorrection', ...
        parser.Results.MassCorrectionMaxRelativeCorrection, ...
        'ReturnSolution', parser.Results.ReturnSolution, ...
        'ExactStartup', parser.Results.ExactStartup, ...
        'GlobalSolveMethod', parser.Results.GlobalSolveMethod, ...
        'GlobalGMRESTolerance', parser.Results.GlobalGMRESTolerance, ...
        'GlobalGMRESRestart', parser.Results.GlobalGMRESRestart, ...
        'GlobalGMRESMaxIterations', parser.Results.GlobalGMRESMaxIterations, ...
        'GlobalDefectSweeps', parser.Results.GlobalDefectSweeps, ...
        'GlobalILUDropTolerance', parser.Results.GlobalILUDropTolerance, ...
        'GlobalILURefreshInterval', parser.Results.GlobalILURefreshInterval, ...
        'GlobalILURefreshIterationThreshold', parser.Results.GlobalILURefreshIterationThreshold, ...
        'UseRearrangement', parser.Results.UseRearrangement, ...
        'QualityThreshold', parser.Results.QualityThreshold, ...
        'PredictiveLookaheadSteps', parser.Results.PredictiveLookaheadSteps, ...
        'MinStepsBetweenRearrangements', parser.Results.MinStepsBetweenRearrangements, ...
        'RearrangementTransferMode', parser.Results.RearrangementTransferMode, ...
        'ImagePath', imgPath, ...
        'ResultsPath', matPath);

    allResults(k).caseName = caseNames(k);
    allResults(k).results = results;
end
end

function problem = localProblem(caseName)
problem = moving_surface_adr_tp_problem(caseName);
end

function problem = sphereProblem(label, titleText, imageName, resultsName, radius, radiusPrime, rotation, finalTime)
problem.label = label;
problem.title = titleText;
problem.imageName = imageName;
problem.resultsName = resultsName;
problem.finalTime = finalTime;
problem.defaultXiVals = [2, 4, 6];
problem.defaultNVals = paperNVals();
problem.geometry = @(N, t) sphereGeometry(N, t, radius, rotation);
problem.exact = @(t, material) exp(-t) .* sphereBase(material.U);
problem.forcing = @(t, material, mu) sphereForcing(t, material.U, mu, radius, radiusPrime);
end

function geom = sphereGeometry(N, t, radius, rotation)
U = kp.geometry.fibonacciSphere(N);
U = U ./ vecnorm(U, 2, 2);
Q = rotation(t);
Y = U * Q.';
R = radius(t);
geom.X = R .* Y;
geom.normals = Y;
geom.h = R * sqrt(4 * pi / N);
geom.material = struct('U', U);
end

function f = sphereForcing(t, U, mu, radius, radiusPrime)
R = radius(t);
Rp = radiusPrime(t);
b = sphereBase(U);
lapb = sphereBaseLaplacian(U);
c = exp(-t) .* b;
materialDerivative = -c;
divv = 2 * Rp / R;
lapc = exp(-t) .* lapb ./ (R^2);
f = materialDerivative + c .* divv - mu .* lapc;
end

function b = sphereBase(U)
ux = U(:, 1);
uy = U(:, 2);
b = 0.8 + 0.15 * ux + 0.05 * (ux.^2 - uy.^2);
end

function lapb = sphereBaseLaplacian(U)
ux = U(:, 1);
uy = U(:, 2);
lapb = -0.30 * ux - 0.30 * (ux.^2 - uy.^2);
end

function Q = rotationZ(theta)
c = cos(theta);
s = sin(theta);
Q = [c, -s, 0; s, c, 0; 0, 0, 1];
end

function problem = anisotropicEllipsoidProblem()
problem.label = "Anisotropic ellipsoid ADR TP convergence";
problem.title = "Anisotropic ellipsoid ADR tangent-plane convergence";
problem.imageName = "moving_surface_adr_tp_convergence_anisotropic_ellipsoid.png";
problem.resultsName = "moving_surface_adr_tp_convergence_anisotropic_ellipsoid_results.mat";
problem.finalTime = 0.12;
problem.defaultXiVals = [2, 4, 6];
problem.defaultNVals = paperNVals();
problem.geometry = @ellipsoidGeometry;
problem.exact = @ellipsoidExact;
problem.forcing = @ellipsoidForcing;
end

function geom = ellipsoidGeometry(N, t)
U = kp.geometry.fibonacciSphere(N);
U = U ./ vecnorm(U, 2, 2);
axesLengths = ellipsoidAxes(t);
X = U .* axesLengths;
g = X ./ (axesLengths.^2);
normals = g ./ vecnorm(g, 2, 2);
geom.X = X;
geom.normals = normals;
geom.h = sqrt(ellipsoidAreaApprox(axesLengths) / N);
geom.meanCurvature = ellipsoidMeanCurvature(X, axesLengths);
geom.material = struct('U', U);
end

function c = ellipsoidExact(t, material)
U = material.U;
axesLengths = ellipsoidAxes(t);
X = U .* axesLengths;
alpha = [0.15, -0.08, 0.05];
b = 0.8 + X * alpha.';
c = exp(-t) .* b;
end

function f = ellipsoidForcing(t, material, mu)
U = material.U;
axesLengths = ellipsoidAxes(t);
axesPrime = ellipsoidAxesPrime(t);
X = U .* axesLengths;
V = U .* axesPrime;
geom = ellipsoidGeometry(size(U, 1), t);
alpha = [0.15, -0.08, 0.05];
b = 0.8 + X * alpha.';
c = exp(-t) .* b;
materialDerivative = exp(-t) .* (-b + V * alpha.');
divv = ellipsoidSurfaceDivergence(U, axesLengths, axesPrime);
lapX = -2 * geom.meanCurvature .* geom.normals;
lapc = exp(-t) .* (lapX * alpha.');
f = materialDerivative + c .* divv - mu .* lapc;
end

function axesLengths = ellipsoidAxes(t)
axesLengths = [ ...
    1.30 + 0.08 * sin(t), ...
    0.90 + 0.06 * cos(1.5 * t), ...
    0.70 + 0.05 * sin(2.0 * t)];
end

function axesPrime = ellipsoidAxesPrime(t)
axesPrime = [ ...
    0.08 * cos(t), ...
    -0.09 * sin(1.5 * t), ...
    0.10 * cos(2.0 * t)];
end

function divv = ellipsoidSurfaceDivergence(U, axesLengths, axesPrime)
logDetRate = sum(axesPrime ./ axesLengths);
s2 = sum((U.^2) ./ (axesLengths.^2), 2);
logStretchRate = -sum((U.^2) .* (axesPrime ./ (axesLengths.^3)), 2) ./ s2;
divv = logDetRate + logStretchRate;
end

function H = ellipsoidMeanCurvature(X, axesLengths)
A = 1 ./ (axesLengths.^2);
g = X .* A;
s = vecnorm(g, 2, 2);
trA = sum(A);
gAg = sum((g.^2) .* A, 2);
divn = trA ./ s - gAg ./ (s.^3);
H = 0.5 * divn;
end

function area = ellipsoidAreaApprox(axesLengths)
p = 1.6075;
a = axesLengths(1);
b = axesLengths(2);
c = axesLengths(3);
area = 4 * pi * ((a^p * b^p + a^p * c^p + b^p * c^p) / 3)^(1 / p);
end

function problem = breathingTorusProblem()
problem.label = "Breathing torus ADR TP convergence";
problem.title = "Breathing torus ADR tangent-plane convergence";
problem.imageName = "moving_surface_adr_tp_convergence_breathing_torus.png";
problem.resultsName = "moving_surface_adr_tp_convergence_breathing_torus_results.mat";
problem.finalTime = 0.12;
problem.defaultXiVals = [2, 4, 6];
problem.defaultNVals = paperNVals();
problem.geometry = @torusGeometry;
problem.exact = @(t, material) exp(-t) .* torusBase(material.theta, material.phi);
problem.forcing = @torusForcing;
end

function Nvals = paperNVals()
% Six levels up to the largest 2D scale used in the 2021 JCP paper.
Nvals = [576, 1024, 1600, 2500, 3600, 4900];
end

function geom = torusGeometry(N, t)
[theta, phi] = torusMaterialGrid(N);
[majorRadius, minorRadius] = torusRadii(t);
ct = cos(theta);
st = sin(theta);
cp = cos(phi);
sp = sin(phi);
q = majorRadius + minorRadius * ct;
geom.X = [q .* cp, q .* sp, minorRadius .* st];
geom.normals = [ct .* cp, ct .* sp, st];
geom.h = sqrt(4 * pi^2 * majorRadius * minorRadius / N);
geom.material = struct('theta', theta, 'phi', phi);
end

function [theta, phi] = torusMaterialGrid(N)
m = round(sqrt(N));
if m * m == N
    thetaVals = areaUniformTorusTheta(m);
    phiVals = 2 * pi * ((0:m-1).' + 0.5) / m;
else
    nTheta = m;
    nPhi = ceil(N / nTheta);
    thetaVals = areaUniformTorusTheta(nTheta);
    phiVals = 2 * pi * ((0:nPhi-1).' + 0.5) / nPhi;
end

[Theta, Phi] = ndgrid(thetaVals, phiVals);
theta = Theta(:);
phi = Phi(:);
theta = theta(1:N);
phi = phi(1:N);

% Break tensor-product stencil symmetries while retaining analytic material
% coordinates. The perturbation is a small fixed fraction of the grid spacing.
k = (1:N).';
dtheta = 2 * pi / numel(thetaVals);
dphi = 2 * pi / numel(phiVals);
theta = mod(theta + 0.13 * dtheta * sin(12.9898 * k), 2 * pi);
phi = mod(phi + 0.13 * dphi * sin(78.2330 * k), 2 * pi);
end

function theta = areaUniformTorusTheta(nTheta)
% Equal-area theta nodes for the initial torus radii.
[majorRadius, minorRadius] = torusRadii(0.0);
rho = minorRadius / majorRadius;
target = 2 * pi * ((0:nTheta-1).' + 0.5) / nTheta;
theta = target;
for it = 1:8
    theta = theta - (theta + rho * sin(theta) - target) ./ ...
        (1 + rho * cos(theta));
end
theta = mod(theta, 2 * pi);
end

function f = torusForcing(t, material, mu)
theta = material.theta;
phi = material.phi;
[majorRadius, minorRadius] = torusRadii(t);
[majorPrime, minorPrime] = torusRadiiPrime(t);
b = torusBase(theta, phi);
c = exp(-t) .* b;
materialDerivative = -c;
divv = minorPrime / minorRadius + ...
    (majorPrime + minorPrime * cos(theta)) ./ ...
    (majorRadius + minorRadius * cos(theta));
lapc = exp(-t) .* torusBaseLaplacian(theta, phi, majorRadius, minorRadius);
f = materialDerivative + c .* divv - mu .* lapc;
end

function [majorRadius, minorRadius] = torusRadii(t)
majorRadius = 1.35 + 0.08 * sin(t);
minorRadius = 0.38 + 0.04 * cos(1.3 * t);
end

function [majorPrime, minorPrime] = torusRadiiPrime(t)
majorPrime = 0.08 * cos(t);
minorPrime = -0.052 * sin(1.3 * t);
end

function b = torusBase(theta, phi)
b = 0.8 + 0.12 * cos(theta) + 0.08 * sin(phi) + ...
    0.05 * cos(theta - 2 * phi);
end

function lapb = torusBaseLaplacian(theta, phi, majorRadius, minorRadius)
q = majorRadius + minorRadius * cos(theta);
bTheta = -0.12 * sin(theta) - 0.05 * sin(theta - 2 * phi);
bThetaTheta = -0.12 * cos(theta) - 0.05 * cos(theta - 2 * phi);
bPhiPhi = -0.08 * sin(phi) - 0.20 * cos(theta - 2 * phi);
lapb = bThetaTheta ./ (minorRadius^2) ...
    - sin(theta) .* bTheta ./ (minorRadius .* q) ...
    + bPhiPhi ./ (q.^2);
end
