function out = puPatchGeometry(domain, xi, patch_spacing_factor, patch_radius_factor)
%PUPATCHGEOMETRY Shared PU patch geometry used by the MATLAB PU solvers.

if nargin < 3 || isempty(patch_spacing_factor)
    patch_spacing_factor = 0.0;
end
if nargin < 4 || isempty(patch_radius_factor)
    patch_radius_factor = 0.0;
end

domain.buildStructs();
Xi = domain.getInteriorNodes();
Xb = domain.getBdryNodes();
Xf = domain.getAllNodes();
h = domain.getSepRad();

% Match the C++ PU collocation helper: tune radius/spacing from the
% physical interior/boundary cloud, then build the actual patch covering on
% the full ghost-node cloud.
tuning = autoTunePatchGeometry(Xi, Xb, xi, h, patch_radius_factor, patch_spacing_factor);
spacing = tuning.spacing;
radius = tuning.radius;
min_nodes = tuning.min_required_nodes;
centers = choosePatchCenters(Xf, spacing);
node_ids = buildPatchNodeIds(domain, centers, radius, min_nodes);

out = struct();
out.centers = centers;
out.radius = radius;
out.spacing = spacing;
out.min_patch_nodes = min_nodes;
out.node_ids = node_ids;
out.center_tree = buildCenterTree(centers);
out.stencil_props = buildStencilProperties(size(Xf, 2), xi);
out.stencils = buildPatchStencils(Xf, node_ids, out.stencil_props);
end

function tuning = autoTunePatchGeometry(Xi, Xb, xi, h, patch_radius_factor, patch_spacing_factor)
dim = size(Xi, 2);
if dim == 0
    dim = size(Xb, 2);
end
npoly = nchoosek(xi + dim, dim);
min_required = max(2 * npoly + 1, ceil(2.5 * npoly));

if patch_radius_factor > 0
    radius = patch_radius_factor * h;
else
    radius = h * max(2.5, sqrt(double(min_required) / 3.0));
end
if patch_spacing_factor > 0
    spacing = patch_spacing_factor * h;
else
    spacing = 0.9 * radius;
end

X = [Xi; Xb];
tuning = struct('radius', radius, 'spacing', spacing, 'min_required_nodes', min_required);
best = struct();
haveBest = false;
for iter = 1:20
    if patch_spacing_factor > 0
        spacingCandidates = patch_spacing_factor * h;
    else
        spacingCandidates = [1.35, 1.20, 1.05, 0.9] * radius;
    end

    for spacingTry = spacingCandidates
        centerIds = choosePatchCenterIds(X, spacingTry);
        minPatchNodes = inf;
        maxPatchNodes = 0;
        for cid = centerIds(:).'
            d = sqrt(sum((X - X(cid, :)).^2, 2));
            count = nnz(d <= radius);
            minPatchNodes = min(minPatchNodes, count);
            maxPatchNodes = max(maxPatchNodes, count);
        end
        if isempty(centerIds)
            minPatchNodes = 0;
        end

        candidate = tuning;
        candidate.radius = radius;
        candidate.spacing = spacingTry;
        candidate.min_patch_nodes = minPatchNodes;
        candidate.max_patch_nodes = maxPatchNodes;
        candidate.num_centers = numel(centerIds);
        [candidate.avg_nnz_per_row, candidate.max_nnz_per_row] = estimateCollocationSparsity(X, centerIds, radius);
        if isempty(X)
            candidate.avg_nnz_fraction = 1.0;
            candidate.max_nnz_fraction = 1.0;
        else
            candidate.avg_nnz_fraction = candidate.avg_nnz_per_row / size(X, 1);
            candidate.max_nnz_fraction = candidate.max_nnz_per_row / size(X, 1);
        end

        if candidate.min_patch_nodes >= min_required && ...
                (~haveBest || candidate.avg_nnz_fraction < best.avg_nnz_fraction || ...
                (abs(candidate.avg_nnz_fraction - best.avg_nnz_fraction) < 1e-14 && candidate.radius < best.radius))
            best = candidate;
            haveBest = true;
        end

        if ~isempty(centerIds) && candidate.min_patch_nodes >= min_required && candidate.avg_nnz_fraction <= 0.10
            tuning = candidate;
            return;
        end
    end

    if patch_radius_factor > 0 && patch_spacing_factor > 0
        error('kp:solvers:BadPatchGeometry', 'Manual PU patch geometry does not provide enough nodes per patch.');
    end

    radius = 1.12 * radius;
end

if haveBest && best.min_patch_nodes >= min_required
    tuning = best;
    return;
end

error('kp:solvers:PatchTuningFailed', 'Failed to auto-tune PU patch geometry.');
end

function centers = choosePatchCenters(X, spacing)
if isempty(X)
    centers = zeros(0, size(X, 2));
    return;
end
centers = X(choosePatchCenterIds(X, spacing), :);
end

function ids = choosePatchCenterIds(X, spacing)
if isempty(X)
    ids = zeros(0, 1);
    return;
end
remaining = true(size(X, 1), 1);
ids = zeros(0, 1);
for i = 1:size(X, 1)
    if ~remaining(i)
        continue;
    end
    ids = [ids; i]; %#ok<AGROW>
    d = sqrt(sum((X - X(i, :)).^2, 2));
    remaining(d <= spacing) = false;
end
end

function [avgNnz, maxNnz] = estimateCollocationSparsity(X, centerIds, radius)
if isempty(X) || isempty(centerIds)
    avgNnz = 0.0;
    maxNnz = 0.0;
    return;
end
counts = zeros(numel(centerIds), 1);
for k = 1:numel(centerIds)
    d = sqrt(sum((X - X(centerIds(k), :)).^2, 2));
    mask = d <= radius;
    membership = sqrt(sum((X(centerIds, :) - X(centerIds(k), :)).^2, 2)) <= radius;
    if any(membership)
        centerSubset = centerIds(membership);
        cover = false(size(X, 1), 1);
        for j = 1:numel(centerSubset)
            dj = sqrt(sum((X - X(centerSubset(j), :)).^2, 2));
            cover = cover | (dj <= radius);
        end
        counts(k) = nnz(mask & cover);
    else
        counts(k) = nnz(mask);
    end
end
avgNnz = mean(counts);
maxNnz = max(counts);
end

function patch_node_ids = buildPatchNodeIds(domain, centers, radius, min_patch_nodes)
patch_node_ids = cell(size(centers, 1), 1);
num_all = domain.getNumTotalNodes();
for p = 1:size(centers, 1)
    ids = domain.queryBall("all", centers(p, :), radius);
    patch_node_ids{p} = ids{1}(:);
    if numel(patch_node_ids{p}) < min_patch_nodes
        [idsKnn, ~] = domain.queryKnn("all", centers(p, :), min(min_patch_nodes, num_all));
        patch_node_ids{p} = idsKnn(1, :).';
    end
end
end

function sp = buildStencilProperties(dim, xi)
sp = kp.rbffd.StencilProperties();
sp.dim = dim;
sp.ell = max(xi + 1, 2);
sp.npoly = size(kp.poly.total_degree_indices(dim, sp.ell), 1);
sp.spline_degree = max(5, sp.ell);
if mod(sp.spline_degree, 2) == 0
    sp.spline_degree = sp.spline_degree - 1;
end
end

function stencils = buildPatchStencils(Xf, node_ids, sp)
stencils = cell(size(node_ids));
for p = 1:numel(node_ids)
    stencil = kp.rbffd.RBFStencil();
    stencil.InitializeGeometry(Xf(node_ids{p}, :), sp);
    stencils{p} = stencil;
end
end

function tree = buildCenterTree(centers)
tree = struct('Points', centers, 'Searcher', [], 'HasSearcher', false);
if isempty(centers)
    return;
end
if exist('KDTreeSearcher', 'class') == 8 && exist('rangesearch', 'file') == 2
    tree.Searcher = KDTreeSearcher(centers);
    tree.HasSearcher = true;
end
end
