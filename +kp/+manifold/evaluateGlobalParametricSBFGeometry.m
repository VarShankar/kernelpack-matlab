function [X, normals, tangentTheta, tangentPhi] = evaluateGlobalParametricSBFGeometry(model, theta, phi)
%EVALUATEGLOBALPARAMETRICSBFGEOMETRY Evaluate a periodic SBF surface model.
%   [X,N,TTH,TPH] = EVALUATEGLOBALPARAMETRICSBFGEOMETRY(MODEL,THETA,PHI)
%   returns surface positions, normals, and parameter tangents at the query
%   sites. MODEL is produced by kp.manifold.fitGlobalParametricSBFGeometry.

if ~isfield(model, 'kind') || string(model.kind) ~= "periodic2"
    error('kp:manifold:BadSBFModel', ...
        'Expected a periodic two-parameter SBF geometry model.');
end

theta = theta(:);
phi = phi(:);
[r, dtheta, dphi] = periodicChordDistance(theta, phi, ...
    model.thetaControl, model.phiControl);
kernel = kp.geometry.phsKernel(r, model.degree);
X = kernel * model.coefficients;

factor = kp.geometry.RBFLevelSet.radialDerivativeFactor(r, model.degree);
tangentTheta = (factor .* sin(dtheta)) * model.coefficients;
tangentPhi = (factor .* sin(dphi)) * model.coefficients;
normals = kp.geometry.normalizeRows(cross(tangentTheta, tangentPhi, 2));
end

function [r, dtheta, dphi] = periodicChordDistance(thetaA, phiA, thetaB, phiB)
dtheta = thetaA(:) - thetaB(:).';
dphi = phiA(:) - phiB(:).';
r2 = 2 - 2 * cos(dtheta) + 2 - 2 * cos(dphi);
r = sqrt(max(r2, 0));
end
