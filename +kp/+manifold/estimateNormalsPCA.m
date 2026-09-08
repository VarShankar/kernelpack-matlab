function [normals, info] = estimateNormalsPCA(X, varargin)
%ESTIMATENORMALSPCA Estimate oriented point-cloud normals by local PCA.
%   NORMALS = ESTIMATENORMALSPCA(X) estimates one unit normal per point from
%   the smallest principal component of each local neighborhood.  If reference
%   normals are supplied, signs are chosen to agree with the reference.  If no
%   reference is supplied, signs are oriented away from the point-cloud
%   centroid.

parser = inputParser();
parser.addParameter('NumNeighbors', 32, @(x) isnumeric(x) && isscalar(x) && x >= 4);
parser.addParameter('ReferenceNormals', [], @(x) isempty(x) || isnumeric(x));
parser.addParameter('OrientationPoint', [], @(x) isempty(x) || isnumeric(x));
parser.parse(varargin{:});

n = size(X, 1);
dim = size(X, 2);
if dim ~= 3
    error('kp:manifold:BadNormalDimension', ...
        'PCA normal estimation currently expects 3D embedded points.');
end

k = min(n, max(4, round(parser.Results.NumNeighbors)));
ids = localNeighborIds(X, k);
normals = zeros(n, 3);

for i = 1:n
    Xi = X(ids(i, :), :);
    Xi = Xi - mean(Xi, 1);
    C = Xi.' * Xi;
    [V, D] = eig((C + C.') * 0.5);
    [~, j] = min(diag(D));
    ni = V(:, j).';
    normals(i, :) = ni ./ max(norm(ni), eps);
end

ref = parser.Results.ReferenceNormals;
if ~isempty(ref)
    dots = sum(normals .* ref, 2);
    normals(dots < 0, :) = -normals(dots < 0, :);
elseif ~isempty(parser.Results.OrientationPoint)
    center = parser.Results.OrientationPoint(:).';
    dots = sum(normals .* (X - center), 2);
    normals(dots < 0, :) = -normals(dots < 0, :);
else
    center = mean(X, 1);
    dots = sum(normals .* (X - center), 2);
    normals(dots < 0, :) = -normals(dots < 0, :);
end

info = struct();
info.numNeighbors = k;
info.rmsVectorError = NaN;
info.maxAngleDegrees = NaN;
if ~isempty(ref)
    ref = ref ./ max(vecnorm(ref, 2, 2), eps);
    diff = normals - ref;
    info.rmsVectorError = sqrt(mean(sum(diff.^2, 2)));
    cosang = min(max(sum(normals .* ref, 2), -1), 1);
    info.maxAngleDegrees = max(acosd(cosang));
    info.rmsAngleDegrees = sqrt(mean(acosd(cosang).^2));
end
end

function ids = localNeighborIds(X, k)
if exist('KDTreeSearcher', 'class') == 8 && exist('knnsearch', 'file') == 2
    tree = KDTreeSearcher(X);
    ids = knnsearch(tree, X, 'K', k);
    return;
end

n = size(X, 1);
ids = zeros(n, k);
for i = 1:n
    d2 = sum((X - X(i, :)).^2, 2);
    [~, order] = sort(d2, 'ascend');
    ids(i, :) = order(1:k);
end
end
