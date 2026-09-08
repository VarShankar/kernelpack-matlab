function [A, B, geom] = tangentPlaneLocalSystem(x, normal, rbf, drbfor, d2rbf, polyIndices, recurrence)
%TANGENTPLANELOCALSYSTEM Build one tangent-plane augmented RBF-FD system.

geom = kp.manifold.detail.tangentPlaneLocalGeometry(x, normal);

P = kp.poly.mpoly_eval(geom.xc, polyIndices, recurrence);
Lp = (kp.poly.mpoly_eval(geom.xc(1, :), polyIndices, recurrence, [2, 0]) + ...
      kp.poly.mpoly_eval(geom.xc(1, :), polyIndices, recurrence, [0, 2])) ./ (geom.w2^2);
Gpx = kp.poly.mpoly_eval(geom.xc(1, :), polyIndices, recurrence, [1, 0]) ./ geom.w2;
Gpy = kp.poly.mpoly_eval(geom.xc(1, :), polyIndices, recurrence, [0, 1]) ./ geom.w2;

% Scale the PHS block with the same local width used by the polynomial
% basis. Derivative right-hand sides are then mapped back to physical
% coordinates, preserving the exact scale invariance of RBF-FD weights.
Ar = rbf(1, geom.rdc);
dphidror = drbfor(1, geom.rdc(:, 1));
Lr = (d2rbf(1, geom.rdc(:, 1)) + dphidror) ./ (geom.w2^2);
dxPhi = -geom.xc .* dphidror ./ geom.w2;

A = [[Ar, P]; [P.', zeros(size(P, 2))]];
B = [[Lr, dxPhi]; [Lp.', Gpx.', Gpy.']];
end
