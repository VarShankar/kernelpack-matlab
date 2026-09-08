function [gamma, info] = hyperviscosityCoefficient( ...
    L, Gx, Gy, Gz, X, nr, k, h, targetOrder, minimumPower)
%HYPERVISCOSITYCOEFFICIENT Paper/source hyperviscosity scaling for surfaces.
%   Pass k=NaN and a targetOrder to select the smallest power satisfying
%   2*k-max(q) >= targetOrder. minimumPower prevents a moving-surface run
%   from decreasing a previously accepted power.

if nargin < 8 || isempty(h)
    h = sqrt(1 / size(X, 1));
end
if nargin < 9
    targetOrder = [];
end
if nargin < 10 || isempty(minimumPower)
    minimumPower = 1;
end

opts.tol = 1.0e-2;
taux = largestRealEigenvalue(Gx, opts);
tauy = largestRealEigenvalue(Gy, opts);
tauz = largestRealEigenvalue(Gz, opts);

[q1, q2, q3] = kp.manifold.EstimateGrowthSurfaceGrad(Gx, Gy, Gz, taux, tauy, tauz, X, nr, h);
q = [q1, q2, q3];
automaticPower = isempty(k) || (isscalar(k) && isnan(k));
if automaticPower
    if isempty(targetOrder) || ~isscalar(targetOrder) || ...
            ~isfinite(targetOrder) || targetOrder <= 0
        error('kp:manifold:MissingHyperviscosityTargetOrder', ...
            'Automatic hyperviscosity power selection requires a positive target order.');
    end
    k = max(round(minimumPower), ceil((targetOrder + max(q)) / 2));
else
    validateattributes(k, {'numeric'}, {'scalar', 'integer', 'positive', 'finite'});
end

hyptermComponents = [ ...
    taux * 2^(q1 - 2 * k) * h^(2 * k - q1), ...
    tauy * 2^(q2 - 2 * k) * h^(2 * k - q2), ...
    tauz * 2^(q3 - 2 * k) * h^(2 * k - q3)];
hypterm = sum(hyptermComponents);
[eta, ~] = kp.manifold.GetGrowthFunction(L, X, h, k);
gammaComponents = 3^(-k) * (-1)^(1 - k) * ...
    hyptermComponents ./ mean(eta);
gamma = sum(gammaComponents);
componentWeights = abs(hyptermComponents) ./ max(sum(abs(hyptermComponents)), eps);
geometryExponent = sum(componentWeights .* (2 * k - q - 1));

info = struct( ...
    'tau', [taux, tauy, tauz], ...
    'q', q, ...
    'etaMean', mean(eta), ...
    'hyptermComponents', hyptermComponents, ...
    'hypterm', hypterm, ...
    'gammaComponents', gammaComponents, ...
    'geometryExponent', geometryExponent, ...
    'h', h, ...
    'k', k, ...
    'targetOrder', targetOrder, ...
    'automaticPower', automaticPower, ...
    'consistencyMargin', 2 * k - max(q) - valueOrNaN(targetOrder));
end

function value = valueOrNaN(value)
if isempty(value)
    value = NaN;
end
end

function val = largestRealEigenvalue(A, opts)
val = tryLargestReal(A, opts);
while isnan(val)
    opts.tol = 2 * opts.tol;
    val = tryLargestReal(A, opts);
end
end

function val = tryLargestReal(A, opts)
warnState = warning();
cleanup = onCleanup(@() warning(warnState));
warning('off', 'all');
val = real(eigs(A, 1, 'LR', opts));
end
