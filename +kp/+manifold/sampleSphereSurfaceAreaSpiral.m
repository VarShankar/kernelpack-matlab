function [U, info] = sampleSphereSurfaceAreaSpiral(N, areaDensityFunction, varargin)
%SAMPLESPHERESURFACEAREASPIRAL Area-corrected low-discrepancy sphere samples.
%   U = SAMPLESPHERESURFACEAREASPIRAL(N, JFUN) returns N sites on the
%   reference sphere whose density is proportional to a sphere-like
%   surface-area density JFUN(U).  For a geometric or SBF map X(U), JFUN
%   should satisfy dA_X = JFUN(U(theta,z))*dtheta*dz.
%
%   Unlike embedded FPS, this keeps the smooth spiral ordering/structure that
%   high-order RBF-FD stencils tend to like, while correcting the density for
%   nonuniform surface stretch.

parser = inputParser();
parser.addParameter('NumEta', 2001, @(x) isnumeric(x) && isscalar(x) && x >= 16);
parser.addParameter('NumTheta', 721, @(x) isnumeric(x) && isscalar(x) && x >= 16);
parser.parse(varargin{:});
opts = parser.Results;

N = max(1, round(N));
numEta = round(opts.NumEta);
numTheta = round(opts.NumTheta);
etaGrid = linspace(-1 + 1.0e-10, 1 - 1.0e-10, numEta).';
thetaGrid = linspace(0, 2 * pi, numTheta);

marginal = zeros(numEta, 1);
for j = 1:numEta
    eta = etaGrid(j);
    r = sqrt(max(1 - eta^2, 0));
    Utheta = [r * cos(thetaGrid(:)), r * sin(thetaGrid(:)), ...
        eta * ones(numTheta, 1)];
    density = max(areaDensityFunction(Utheta), 0);
    marginal(j) = trapz(thetaGrid, density);
end

cdfEta = cumtrapz(etaGrid, marginal);
cdfEta = normalizeCDF(cdfEta);
levelsEta = ((0:N - 1).' + 0.5) / N;
eta = interp1(cdfEta, etaGrid, levelsEta, 'linear', 'extrap');

golden = (sqrt(5) - 1) / 2;
levelsTheta = mod((0:N - 1).' * golden + 0.5 / N, 1);
theta = zeros(N, 1);
for i = 1:N
    r = sqrt(max(1 - eta(i)^2, 0));
    Utheta = [r * cos(thetaGrid(:)), r * sin(thetaGrid(:)), ...
        eta(i) * ones(numTheta, 1)];
    density = max(areaDensityFunction(Utheta), 0);
    cdfTheta = normalizeCDF(cumtrapz(thetaGrid, density));
    theta(i) = interp1(cdfTheta, thetaGrid, levelsTheta(i), 'linear', 'extrap');
end

r = sqrt(max(1 - eta .^ 2, 0));
U = [r .* cos(theta), r .* sin(theta), eta];
U = U ./ max(vecnorm(U, 2, 2), eps);

info = struct();
info.numSamples = N;
info.numEta = numEta;
info.numTheta = numTheta;
info.areaMarginalMin = min(marginal);
info.areaMarginalMax = max(marginal);
end

function cdf = normalizeCDF(cdf)
cdf = cdf(:);
span = cdf(end) - cdf(1);
if abs(span) <= eps
    cdf = linspace(0, 1, numel(cdf)).';
else
    cdf = (cdf - cdf(1)) ./ span;
end
cdf = max(cdf, 0);
cdf = min(cdf, 1);
for k = 2:numel(cdf)
    if cdf(k) <= cdf(k - 1)
        cdf(k) = min(1, cdf(k - 1) + eps);
    end
end
cdf(end) = 1;
end
