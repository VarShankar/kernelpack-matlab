function [normals, info, state] = estimateNormalsGlobalSBF(geom, state, varargin)
%ESTIMATENORMALSGLOBALSBF Normals from a cached global parametric SBF fit.
%   The material parameter sites are fixed for a Lagrangian moving-surface
%   run, so the SBF interpolation matrix and parameter-derivative kernels can
%   be cached once.  At each time level only the coordinate coefficients are
%   updated by applying the cached LU factorization to the current XYZ data.
%
%   By default every surface node is used as an SBF center.  Passing
%   ControlPointCount=M uses a deterministic farthest-point subset of M
%   material sites, forms the M-by-M interpolation matrix there, and evaluates
%   the parameter derivatives at all N surface nodes.

parser = inputParser();
parser.addParameter('Degree', 7, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('ControlPointCount', Inf, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('ControlPointIds', [], @(x) isempty(x) || isnumeric(x));
parser.parse(varargin{:});
degree = parser.Results.Degree;
controlPointCount = parser.Results.ControlPointCount;
controlPointIds = parser.Results.ControlPointIds;

if nargin < 2 || isempty(state) || ...
        ~isReusableState(state, geom, degree, controlPointCount, controlPointIds)
    state = buildState(geom, degree, controlPointCount, controlPointIds);
end

coeffs = state.factor \ geom.X(state.controlIds, :);
tangent1 = state.derivativeKernel1 * coeffs;
tangent2 = state.derivativeKernel2 * coeffs;
normals = kp.geometry.normalizeRows(cross(tangent1, tangent2, 2));

reference = [];
if isfield(geom, 'normals') && ~isempty(geom.normals)
    reference = geom.normals;
elseif isfield(state, 'previousNormals') && ~isempty(state.previousNormals)
    reference = state.previousNormals;
end

[normals, info] = orientAndMeasure(normals, reference);
info.controlPointCount = numel(state.controlIds);
info.controlPointFraction = numel(state.controlIds) / size(geom.X, 1);
info.controlFillDistance = state.controlInfo.fillDistance;
state.previousNormals = normals;
end

function tf = isReusableState(state, geom, degree, controlPointCount, controlPointIds)
n = size(geom.X, 1);
[requestedCount, hasExplicitIds] = normalizeControlRequest(controlPointCount, controlPointIds, n);
tf = isfield(state, 'degree') && state.degree == degree ...
    && isfield(state, 'numPoints') && state.numPoints == n ...
    && isfield(state, 'controlPointCount') && state.controlPointCount == requestedCount ...
    && isfield(state, 'kind') && strcmp(state.kind, materialKind(geom));
if tf && hasExplicitIds
    tf = isequal(state.controlIds(:), controlPointIds(:));
end
end

function state = buildState(geom, degree, controlPointCount, controlPointIds)
kind = materialKind(geom);
n = size(geom.X, 1);
[requestedCount, hasExplicitIds] = normalizeControlRequest(controlPointCount, controlPointIds, n);
if hasExplicitIds
    [controlIds, controlInfo] = kp.manifold.selectSBFControlPoints(geom, ...
        'ControlPointIds', controlPointIds);
else
    [controlIds, controlInfo] = kp.manifold.selectSBFControlPoints(geom, ...
        'Count', requestedCount);
end
switch kind
    case 'sphere'
        [targets, uv, dpar1, dpar2] = sphereParameters(geom.material.U);
        centers = targets(controlIds, :);
        [r, ~] = kp.geometry.sphereChordDistance(centers, centers);
        [rtarget, ~] = kp.geometry.sphereChordDistance(targets, centers);
        targetFactor = kp.geometry.RBFLevelSet.radialDerivativeFactor(rtarget, degree);
        diff = kp.geometry.RBFLevelSet.differenceTensor(targets, centers);
        dot1 = diff(:, :, 1) .* dpar1(:, 1) ...
            + diff(:, :, 2) .* dpar1(:, 2) ...
            + diff(:, :, 3) .* dpar1(:, 3);
        dot2 = diff(:, :, 1) .* dpar2(:, 1) ...
            + diff(:, :, 2) .* dpar2(:, 2) ...
            + diff(:, :, 3) .* dpar2(:, 3);
        state.parameters = uv;
        state.derivativeKernel1 = targetFactor .* dot1;
        state.derivativeKernel2 = targetFactor .* dot2;

    case 'torus'
        theta = geom.material.theta(:);
        phi = geom.material.phi(:);
        thetaControl = theta(controlIds);
        phiControl = phi(controlIds);
        dtheta = thetaControl - thetaControl.';
        dphi = phiControl - phiControl.';
        r2 = 2 - 2 * cos(dtheta) + 2 - 2 * cos(dphi);
        r = sqrt(max(r2, 0));
        dthetaTarget = theta - thetaControl.';
        dphiTarget = phi - phiControl.';
        r2Target = 2 - 2 * cos(dthetaTarget) + 2 - 2 * cos(dphiTarget);
        rtarget = sqrt(max(r2Target, 0));
        targetFactor = kp.geometry.RBFLevelSet.radialDerivativeFactor(rtarget, degree);
        state.parameters = [theta, phi];
        state.derivativeKernel1 = targetFactor .* sin(dthetaTarget);
        state.derivativeKernel2 = targetFactor .* sin(dphiTarget);

    otherwise
        error('kp:manifold:BadSBFMaterial', ...
            'Global SBF normals need material.U or material.theta/material.phi.');
end

kernel = kp.geometry.phsKernel(r, degree);
reg = 1.0e-12 * max(1.0, max(abs(kernel), [], 'all'));
interp = kernel + reg * eye(size(kernel, 1));

state.kind = kind;
state.degree = degree;
state.numPoints = n;
state.controlPointCount = requestedCount;
state.controlIds = controlIds(:);
state.controlInfo = controlInfo;
state.factor = decomposition(interp, 'lu');
state.previousNormals = [];
end

function [count, hasExplicitIds] = normalizeControlRequest(controlPointCount, controlPointIds, n)
hasExplicitIds = ~isempty(controlPointIds);
if hasExplicitIds
    count = numel(unique(controlPointIds(:)));
else
    if ~isfinite(controlPointCount)
        count = n;
    else
        count = min(n, max(1, round(controlPointCount)));
    end
end
end

function kind = materialKind(geom)
if ~isfield(geom, 'material') || isempty(geom.material)
    kind = '';
elseif isfield(geom.material, 'U') && size(geom.material.U, 2) == 3
    kind = 'sphere';
elseif isfield(geom.material, 'theta') && isfield(geom.material, 'phi')
    kind = 'torus';
else
    kind = '';
end
end

function [centers, uv, du, dv] = sphereParameters(U)
centers = U ./ max(vecnorm(U, 2, 2), eps);
uvw = kp.geometry.cart2sphRows(centers);
uv = uvw(:, 1:2);
az = uv(:, 1);
el = uv(:, 2);
du = [-cos(el) .* sin(az), cos(el) .* cos(az), zeros(size(az))];
dv = [-sin(el) .* cos(az), -sin(el) .* sin(az), cos(el)];
end

function [normals, info] = orientAndMeasure(normals, reference)
info = struct('rmsAngleDegrees', NaN, 'maxAngleDegrees', NaN);
if isempty(reference)
    return;
end

reference = kp.geometry.normalizeRows(reference);
dotprod = sum(normals .* reference, 2);
flip = dotprod < 0;
normals(flip, :) = -normals(flip, :);
dotprod = abs(sum(normals .* reference, 2));
dotprod = min(max(dotprod, -1), 1);
angles = acosd(dotprod);
info.rmsAngleDegrees = sqrt(mean(angles .^ 2, 'omitnan'));
info.maxAngleDegrees = max(angles, [], 'omitnan');
end
