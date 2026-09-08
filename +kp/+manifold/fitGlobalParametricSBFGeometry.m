function model = fitGlobalParametricSBFGeometry(geom, varargin)
%FITGLOBALPARAMETRICSBFGEOMETRY Fit a periodic global SBF surface model.
%   MODEL = FITGLOBALPARAMETRICSBFGEOMETRY(GEOM) fits coordinate functions
%   X(theta,phi) using the periodic material coordinates stored in
%   GEOM.material.theta and GEOM.material.phi.  The returned model can be
%   evaluated with kp.manifold.evaluateGlobalParametricSBFGeometry.

parser = inputParser();
parser.addParameter('Degree', 7, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('ControlPointCount', Inf, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('ControlPointIds', [], @(x) isempty(x) || isnumeric(x));
parser.parse(varargin{:});
opts = parser.Results;

if ~isfield(geom, 'material') || ~isfield(geom.material, 'theta') || ...
        ~isfield(geom.material, 'phi')
    error('kp:manifold:BadSBFMaterial', ...
        'Global parametric SBF geometry needs material.theta and material.phi.');
end

if isempty(opts.ControlPointIds)
    [controlIds, controlInfo] = kp.manifold.selectSBFControlPoints(geom, ...
        'Count', normalizeControlPointCount(opts.ControlPointCount, size(geom.X, 1)));
else
    [controlIds, controlInfo] = kp.manifold.selectSBFControlPoints(geom, ...
        'ControlPointIds', opts.ControlPointIds);
end

theta = geom.material.theta(:);
phi = geom.material.phi(:);
thetaControl = theta(controlIds);
phiControl = phi(controlIds);
r = periodicChordDistance(thetaControl, phiControl, thetaControl, phiControl);
kernel = kp.geometry.phsKernel(r, opts.Degree);
reg = 1.0e-12 * max(1.0, max(abs(kernel), [], 'all'));
interp = kernel + reg * eye(size(kernel, 1));

model = struct();
model.kind = "periodic2";
model.degree = opts.Degree;
model.controlIds = controlIds(:);
model.controlInfo = controlInfo;
model.thetaControl = thetaControl(:);
model.phiControl = phiControl(:);
model.coefficients = interp \ geom.X(controlIds, :);
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
