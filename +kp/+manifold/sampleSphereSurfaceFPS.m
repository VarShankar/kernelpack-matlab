function [U, info] = sampleSphereSurfaceFPS(N, surfaceMap, varargin)
%SAMPLESPHERESURFACEFPS Quasi-uniform samples on a sphere-like surface.
%   U = SAMPLESPHERESURFACEFPS(N, SURFACEMAP) builds a dense deterministic
%   candidate set on the reference sphere, maps it to the embedded surface,
%   and returns N material sites selected by farthest-point sampling in the
%   embedded coordinates.  SURFACEMAP is a function handle X = SURFACEMAP(U).
%
%   An optional AreaDensityFunction handle J = JFUN(U) can be supplied to
%   preselect candidates with density proportional to the surface area
%   element before embedded FPS is applied.  This avoids uniform-reference
%   bias when the geometric/SBF map has nonuniform area stretch.

parser = inputParser();
parser.addParameter('CandidateFactor', 8, @(x) isnumeric(x) && isscalar(x) && x >= 1);
parser.addParameter('RawCandidateFactor', 3, @(x) isnumeric(x) && isscalar(x) && x >= 1);
parser.addParameter('AreaDensityFunction', [], @(x) isempty(x) || isa(x, 'function_handle'));
parser.addParameter('UseAreaWeightedCandidates', true, @(x) islogical(x) || isnumeric(x));
parser.parse(varargin{:});
opts = parser.Results;

N = max(1, round(N));
candidateCount = max(N, ceil(opts.CandidateFactor * N));
rawCount = max(candidateCount, ceil(opts.RawCandidateFactor * candidateCount));

Uraw = kp.geometry.fibonacciSphere(rawCount);
Uraw = Uraw ./ max(vecnorm(Uraw, 2, 2), eps);

if opts.UseAreaWeightedCandidates && ~isempty(opts.AreaDensityFunction)
    density = opts.AreaDensityFunction(Uraw);
    idsCandidate = systematicWeightedSubset(density(:), candidateCount);
else
    idsCandidate = (1:rawCount).';
end

Ucand = Uraw(idsCandidate, :);
Xcand = surfaceMap(Ucand);
ids = embeddedFarthestPointSubset(Xcand, min(N, size(Xcand, 1)));
U = Ucand(ids, :);
U = U ./ max(vecnorm(U, 2, 2), eps);

info = struct();
info.numSamples = size(U, 1);
info.numCandidates = size(Ucand, 1);
info.numRawCandidates = rawCount;
info.candidateFactor = opts.CandidateFactor;
info.rawCandidateFactor = opts.RawCandidateFactor;
info.usedAreaWeightedCandidates = opts.UseAreaWeightedCandidates && ...
    ~isempty(opts.AreaDensityFunction);
info.embeddedQuality = nearestNeighborQuality(surfaceMap(U));
end

function ids = systematicWeightedSubset(weights, count)
weights = max(weights(:), 0);
n = numel(weights);
if n == 0
    ids = zeros(0, 1);
    return;
end
if ~any(weights > 0)
    ids = round(linspace(1, n, min(count, n))).';
    ids = unique(ids, 'stable');
    return;
end

weights = weights ./ sum(weights);
cdf = cumsum(weights);
levels = ((0:count - 1).' + 0.5) / count;
ids = zeros(count, 1);
cursor = 1;
for k = 1:count
    while cursor < n && cdf(cursor) < levels(k)
        cursor = cursor + 1;
    end
    ids(k) = cursor;
end
ids = unique(ids, 'stable');

if numel(ids) < min(count, n)
    missing = setdiff((1:n).', ids, 'stable');
    [~, order] = sort(weights(missing), 'descend');
    fillCount = min(numel(order), min(count, n) - numel(ids));
    ids = [ids; missing(order(1:fillCount))]; %#ok<AGROW>
end
end

function ids = embeddedFarthestPointSubset(X, count)
n = size(X, 1);
count = min(n, max(1, round(count)));
ids = zeros(count, 1);
[~, ids(1)] = max(sum((X - mean(X, 1)) .^ 2, 2));
dist2 = sum((X - X(ids(1), :)) .^ 2, 2);
dist2(ids(1)) = 0;

for k = 2:count
    [~, ids(k)] = max(dist2);
    d2 = sum((X - X(ids(k), :)) .^ 2, 2);
    dist2 = min(dist2, d2);
    dist2(ids(1:k)) = 0;
end
ids = sort(ids);
end

function quality = nearestNeighborQuality(X)
if size(X, 1) < 2
    quality = NaN;
    return;
end
tree = KDTreeSearcher(X);
[~, dist] = knnsearch(tree, X, 'K', 2);
nearest = dist(:, 2);
quality = max(nearest) / max(min(nearest), eps);
end
