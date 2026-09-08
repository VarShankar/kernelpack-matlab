function [ids, info] = selectSBFControlPoints(geom, varargin)
%SELECTSBFCONTROLPOINTS Deterministic material-space SBF control subset.
%   IDS = SELECTSBFCONTROLPOINTS(GEOM,'Count',M) returns M farthest-point
%   material sites.  IDS = SELECTSBFCONTROLPOINTS(GEOM,'FillDistanceTarget',H)
%   keeps adding sites until the material fill distance is at most H.

parser = inputParser();
parser.addParameter('Count', Inf, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('ControlPointIds', [], @(x) isempty(x) || isnumeric(x));
parser.addParameter('FillDistanceTarget', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
parser.addParameter('MaxCount', Inf, @(x) isnumeric(x) && isscalar(x));
parser.parse(varargin{:});
opts = parser.Results;

[embedding, kind] = materialEmbedding(geom);
n = size(embedding, 1);

if ~isempty(opts.ControlPointIds)
    ids = unique(opts.ControlPointIds(:), 'stable');
    validateIds(ids, n);
    minDist2 = minSquaredDistances(embedding, ids);
else
    maxCount = min(n, max(1, round(opts.MaxCount)));
    if isempty(opts.FillDistanceTarget)
        count = requestedCount(opts.Count, n);
        ids = farthestPointSubset(embedding, min(count, maxCount));
        minDist2 = minSquaredDistances(embedding, ids);
    else
        target2 = max(opts.FillDistanceTarget, 0) .^ 2;
        [ids, minDist2] = farthestPointSubsetToFill(embedding, target2, maxCount);
    end
end

info = struct();
info.kind = kind;
info.controlPointCount = numel(ids);
info.controlPointFraction = numel(ids) / n;
info.fillDistance = sqrt(max(minDist2));
info.fillDistanceSquared = max(minDist2);
end

function [embedding, kind] = materialEmbedding(geom)
if ~isfield(geom, 'material') || isempty(geom.material)
    error('kp:manifold:BadSBFMaterial', ...
        'SBF control selection needs material.U or material.theta/material.phi.');
elseif isfield(geom.material, 'U') && size(geom.material.U, 2) == 3
    U = geom.material.U;
    embedding = U ./ max(vecnorm(U, 2, 2), eps);
    kind = 'sphere';
elseif isfield(geom.material, 'theta') && isfield(geom.material, 'phi')
    theta = geom.material.theta(:);
    phi = geom.material.phi(:);
    embedding = [cos(theta), sin(theta), cos(phi), sin(phi)];
    kind = 'torus';
else
    error('kp:manifold:BadSBFMaterial', ...
        'SBF control selection needs material.U or material.theta/material.phi.');
end
end

function count = requestedCount(count, n)
if ~isfinite(count)
    count = n;
else
    count = min(n, max(1, round(count)));
end
end

function validateIds(ids, n)
if any(ids < 1) || any(ids > n)
    error('kp:manifold:BadSBFControlIds', ...
        'ControlPointIds must be valid 1-based row indices.');
end
end

function ids = farthestPointSubset(Y, count)
[ids, ~] = farthestPointSubsetToFill(Y, -Inf, count);
ids = sort(ids);
end

function [ids, minDist2] = farthestPointSubsetToFill(Y, target2, maxCount)
n = size(Y, 1);
maxCount = min(n, max(1, round(maxCount)));
ids = zeros(maxCount, 1);
[~, first] = max(sum(Y .^ 2, 2));
ids(1) = first;
minDist2 = sum((Y - Y(first, :)) .^ 2, 2);
minDist2(first) = 0;
count = 1;

while max(minDist2) > target2 && count < maxCount
    [~, next] = max(minDist2);
    count = count + 1;
    ids(count) = next;
    d2 = sum((Y - Y(next, :)) .^ 2, 2);
    minDist2 = min(minDist2, d2);
end

ids = sort(ids(1:count));
minDist2 = minSquaredDistances(Y, ids);
end

function minDist2 = minSquaredDistances(Y, ids)
minDist2 = inf(size(Y, 1), 1);
for k = 1:numel(ids)
    d2 = sum((Y - Y(ids(k), :)) .^ 2, 2);
    minDist2 = min(minDist2, d2);
end
end
