function problem = moving_surface_adr_tp_problem(caseName)
%MOVING_SURFACE_ADR_TP_PROBLEM Manufactured moving-surface ADR test problem.
%   This helper exposes the same geometry/forcing definitions used by the
%   moving-surface tangent-plane convergence suite so diagnostics can sample
%   the exact geometries without running a PDE solve.

switch lower(string(caseName))
    case "mcf_sphere"
        problem = sphereProblem( ...
            "Mean-curvature-flow sphere ADR TP convergence", ...
            "Mean-curvature-flow sphere ADR tangent-plane convergence", ...
            "moving_surface_adr_tp_convergence_mcf_sphere.png", ...
            "moving_surface_adr_tp_convergence_mcf_sphere_results.mat", ...
            @(t) sqrt(1.4^2 - 4 * t), ...
            @(t) -2 ./ sqrt(1.4^2 - 4 * t), ...
            @(t) eye(3), 0.12);

    case "imcf_sphere"
        problem = sphereProblem( ...
            "Inverse-mean-curvature-flow sphere ADR TP convergence", ...
            "Inverse-MCF sphere ADR tangent-plane convergence", ...
            "moving_surface_adr_tp_convergence_imcf_sphere.png", ...
            "moving_surface_adr_tp_convergence_imcf_sphere_results.mat", ...
            @(t) exp(0.5 * t), ...
            @(t) 0.5 * exp(0.5 * t), ...
            @(t) eye(3), 0.12);

    case "rotating_breathing_sphere"
        omega = 2.0;
        problem = sphereProblem( ...
            "Rotating breathing sphere ADR TP convergence", ...
            "Rotating breathing sphere ADR tangent-plane convergence", ...
            "moving_surface_adr_tp_convergence_rotating_breathing_sphere.png", ...
            "moving_surface_adr_tp_convergence_rotating_breathing_sphere_results.mat", ...
            @(t) 1 + 0.15 * sin(1.5 * t), ...
            @(t) 0.225 * cos(1.5 * t), ...
            @(t) rotationZ(omega * t), 0.12);

    case "anisotropic_ellipsoid"
        problem = anisotropicEllipsoidProblem();

    case "breathing_torus"
        problem = breathingTorusProblem();

    otherwise
        error('kp:examples:UnknownMovingSurfaceCase', ...
            'Unknown moving-surface case "%s".', caseName);
end
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
problem.geometryFromMaterial = @(material, t) sphereGeometryFromMaterial(material, t, radius, rotation);
problem.sampleMaterial = @(N, t) sphereMaterial(N);
problem.backtraceMaterial = @(material, tNow, tPast) material;
problem.exact = @(t, material) exp(-t) .* sphereBase(material.U);
problem.forcing = @(t, material, mu) sphereForcing(t, material.U, mu, radius, radiusPrime);
problem.exactMass = @(t, material, geom, ~) ...
    sum(geom.weights(:) .* exp(-t) .* sphereBase(material.U));
end

function geom = sphereGeometry(N, t, radius, rotation)
geom = sphereGeometryFromMaterial(sphereMaterial(N), t, radius, rotation);
end

function material = sphereMaterial(N)
U = kp.geometry.fibonacciSphere(N);
U = U ./ vecnorm(U, 2, 2);
material = struct('U', U);
end

function geom = sphereGeometryFromMaterial(material, t, radius, rotation)
U = material.U;
N = size(U, 1);
Q = rotation(t);
Y = U * Q.';
R = radius(t);
geom.X = R .* Y;
geom.normals = Y;
geom.weights = repmat(4 * pi * R^2 / N, N, 1);
geom.area = sum(geom.weights);
geom.h = sqrt(geom.area / N);
geom.material = struct('U', U, 'X', geom.X);
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
problem.geometryFromMaterial = @ellipsoidGeometryFromMaterial;
problem.sampleMaterial = @(N, t) sphereMaterial(N);
problem.backtraceMaterial = @(material, tNow, tPast) material;
problem.exact = @ellipsoidExact;
problem.forcing = @ellipsoidForcing;
problem.exactMass = @(t, material, geom, ~) ...
    sum(geom.weights(:) .* ellipsoidExact(t, material));
end

function geom = ellipsoidGeometry(N, t)
geom = ellipsoidGeometryFromMaterial(sphereMaterial(N), t);
end

function geom = ellipsoidGeometryFromMaterial(material, t)
U = material.U;
N = size(U, 1);
axesLengths = ellipsoidAxes(t);
X = U .* axesLengths;
g = X ./ (axesLengths.^2);
normals = g ./ vecnorm(g, 2, 2);
geom.X = X;
geom.normals = normals;
geom.weights = ellipsoidSurfaceWeights(U, axesLengths);
geom.area = sum(geom.weights);
geom.h = sqrt(geom.area / N);
geom.meanCurvature = ellipsoidMeanCurvature(X, axesLengths);
geom.material = struct('U', U, 'X', X);
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
geom = ellipsoidGeometryFromMaterial(material, t);
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

function weights = ellipsoidSurfaceWeights(U, axesLengths)
baseWeight = 4 * pi / size(U, 1);
jacobian = prod(axesLengths) .* vecnorm(U ./ axesLengths, 2, 2);
weights = baseWeight .* jacobian;
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
problem.geometryFromMaterial = @torusGeometryFromMaterial;
problem.sampleMaterial = @(N, t) torusMaterial(N);
problem.backtraceMaterial = @(material, tNow, tPast) material;
problem.exact = @(t, material) exp(-t) .* torusBase(material.theta, material.phi);
problem.forcing = @torusForcing;
problem.exactMass = @(t, material, geom, ~) ...
    sum(geom.weights(:) .* exp(-t) .* torusBase(material.theta, material.phi));
end

function Nvals = paperNVals()
Nvals = [576, 1024, 1600, 2500, 3600, 4900];
end

function geom = torusGeometry(N, t)
geom = torusGeometryFromMaterial(torusMaterial(N), t);
end

function material = torusMaterial(N)
[theta, phi] = torusMaterialGrid(N);
material = struct('theta', theta, 'phi', phi);
end

function geom = torusGeometryFromMaterial(material, t)
theta = material.theta;
phi = material.phi;
N = numel(theta);
[majorRadius, minorRadius] = torusRadii(t);
ct = cos(theta);
st = sin(theta);
cp = cos(phi);
sp = sin(phi);
q = majorRadius + minorRadius * ct;
geom.X = [q .* cp, q .* sp, minorRadius .* st];
geom.normals = [ct .* cp, ct .* sp, st];
geom.weights = torusSurfaceWeights(theta, majorRadius, minorRadius);
geom.area = sum(geom.weights);
geom.h = sqrt(geom.area / N);
geom.material = struct('theta', theta, 'phi', phi, 'X', geom.X);
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

k = (1:N).';
dtheta = 2 * pi / numel(thetaVals);
dphi = 2 * pi / numel(phiVals);
theta = mod(theta + 0.13 * dtheta * sin(12.9898 * k), 2 * pi);
phi = mod(phi + 0.13 * dphi * sin(78.2330 * k), 2 * pi);
end

function theta = areaUniformTorusTheta(nTheta)
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

function weights = torusSurfaceWeights(theta, majorRadius, minorRadius)
[majorRadius0, minorRadius0] = torusRadii(0.0);
area0 = 4 * pi^2 * majorRadius0 * minorRadius0;
jacobian0 = minorRadius0 .* (majorRadius0 + minorRadius0 .* cos(theta));
jacobian = minorRadius .* (majorRadius + minorRadius .* cos(theta));
weights = (area0 / numel(theta)) .* jacobian ./ jacobian0;
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
