function out = puPatchGeometry(domain, xi, patch_spacing_factor, patch_radius_factor)
%PUPATCHGEOMETRY Shared PU patch geometry used by the MATLAB PU solvers.

if nargin < 3 || isempty(patch_spacing_factor)
    patch_spacing_factor = 0.0;
end
if nargin < 4 || isempty(patch_radius_factor)
    patch_radius_factor = 0.0;
end

domain.buildStructs();
Xf = domain.getAllNodes();
h = domain.getSepRad();

spacing = choosePatchSpacing(h, patch_spacing_factor);
radius = choosePatchRadius(h, patch_radius_factor);
min_nodes = chooseMinimumPatchNodes(size(Xf, 2), xi);
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

function spacing = choosePatchSpacing(h, patch_spacing_factor)
if patch_spacing_factor > 0
    spacing = patch_spacing_factor * h;
else
    spacing = 2.0 * h;
end
end

function radius = choosePatchRadius(h, patch_radius_factor)
if patch_radius_factor > 0
    radius = patch_radius_factor * h;
else
    radius = 3.0 * h;
end
end

function centers = choosePatchCenters(X, spacing)
if isempty(X)
    centers = zeros(0, size(X, 2));
    return;
end

remaining = true(size(X, 1), 1);
centers = zeros(0, size(X, 2));
for i = 1:size(X, 1)
    if ~remaining(i)
        continue;
    end
    xi = X(i, :);
    centers = [centers; xi]; %#ok<AGROW>
    d = sqrt(sum((X - xi).^2, 2));
    remaining(d <= spacing) = false;
end
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

function min_nodes = chooseMinimumPatchNodes(dim, xi)
ell = max(xi + 1, 2);
npoly = size(kp.poly.total_degree_indices(dim, ell), 1);
min_nodes = 2 * npoly + 1;
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
