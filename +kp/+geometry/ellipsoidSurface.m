function geom = ellipsoidSurface(n, axesLengths)
%ELLIPSOIDSURFACE Quasi-uniform sampled ellipsoid geometry.
%   geom = kp.geometry.ellipsoidSurface(n, axesLengths) maps Fibonacci
%   sphere nodes to the ellipsoid with semi-axes axesLengths = [a, b, c].
%   The returned struct contains nodes, outward normals, mean curvature,
%   an analytic-area approximation, and nominal node spacing h.

arguments
    n (1,1) double {mustBeInteger, mustBeNonnegative}
    axesLengths (1,3) double {mustBePositive}
end

U = kp.geometry.fibonacciSphere(n);
if n > 0
    U = U ./ vecnorm(U, 2, 2);
end

X = U .* axesLengths;
nr = ellipsoidNormals(X, axesLengths);
H = ellipsoidMeanCurvature(X, axesLengths);
area = ellipsoidAreaApprox(axesLengths);

geom = struct( ...
    'X', X, ...
    'normals', nr, ...
    'meanCurvature', H, ...
    'area', area, ...
    'h', sqrt(area / max(n, 1)), ...
    'axesLengths', axesLengths);
end

function nr = ellipsoidNormals(X, axesLengths)
if isempty(X)
    nr = zeros(0, 3);
    return;
end

g = X ./ (axesLengths.^2);
nr = g ./ vecnorm(g, 2, 2);
end

function H = ellipsoidMeanCurvature(X, axesLengths)
if isempty(X)
    H = zeros(0, 1);
    return;
end

A = 1 ./ (axesLengths.^2);
g = X .* A;
s = vecnorm(g, 2, 2);
trA = sum(A);
gAg = sum((g.^2) .* A, 2);
divn = trA ./ s - gAg ./ (s.^3);
H = 0.5 * divn;
end

function area = ellipsoidAreaApprox(axesLengths)
% Knud Thomsen's formula is accurate enough for nominal h selection.
p = 1.6075;
a = axesLengths(1);
b = axesLengths(2);
c = axesLengths(3);
area = 4 * pi * ((a^p * b^p + a^p * c^p + b^p * c^p) / 3)^(1 / p);
end
