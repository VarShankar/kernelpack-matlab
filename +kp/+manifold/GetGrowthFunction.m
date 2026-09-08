function [eta, inveta] = GetGrowthFunction(L, X, h, k)
%GETGROWTHFUNCTION Hyperviscosity growth calibration for moving-surface ADR.

khat = 2 / h;
real_eig = (-1)^k * (3^k) * (khat^(2 * k));
if size(X, 2) < 3
    error('kp:manifold:GetGrowthFunction3D', 'GetGrowthFunction expects a 3D point cloud.');
end
f = exp(1i * khat * X(:, 1)) .* exp(1i * khat * X(:, 2)) .* exp(1i * khat * X(:, 3));
eta = kp.manifold.applyMatrixPower(L, f, k) ./ (real_eig * f);
eta(isnan(eta)) = 0;
eta = abs(real(eta));
inveta = 1 ./ eta;
end
