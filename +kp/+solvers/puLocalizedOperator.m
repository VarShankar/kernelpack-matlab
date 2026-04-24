function A = puLocalizedOperator(domain, patchData, xi, Xq, opName, varargin)
%PULOCALIZEDOPERATOR Assemble a PU-collocation operator on arbitrary targets.

if isempty(Xq)
    A = sparse(0, domain.getNumTotalNodes());
    return;
end

Xnodes = domain.getAllNodes();
centers = patchData.centers;
node_ids = patchData.node_ids;
cachedStencils = patchData.stencils;
radius = patchData.radius;
dim = size(Xnodes, 2);
num_targets = size(Xq, 1);
num_all = size(Xnodes, 1);
patch_ids_per_query = queryPatchIds(patchData, Xq, radius);

parser = inputParser();
parser.addParameter('Normals', zeros(num_targets, dim));
parser.addParameter('NeuCoeff', zeros(num_targets, 1));
parser.addParameter('DirCoeff', zeros(num_targets, 1));
parser.parse(varargin{:});
normals = parser.Results.Normals;
neuCoeff = parser.Results.NeuCoeff(:);
dirCoeff = parser.Results.DirCoeff(:);

sp = patchData.stencil_props;
switch lower(string(opName))
    case {"interp", "interpolation"}
        theta = 0;
    case {"lap", "laplacian"}
        theta = 2;
    case {"bc", "boundary"}
        theta = 1;
    otherwise
        error('kp:solvers:BadPUOperator', 'Unknown PU operator "%s".', string(opName));
end
if theta ~= 1
    sp.ell = max(xi + theta - 1, 2);
    sp.npoly = size(kp.poly.total_degree_indices(dim, sp.ell), 1);
    sp.spline_degree = max(5, sp.ell);
    if mod(sp.spline_degree, 2) == 0
        sp.spline_degree = sp.spline_degree - 1;
    end
end

rows = cell(num_targets, 1);
cols = cell(num_targets, 1);
vals = cell(num_targets, 1);
for q = 1:num_targets
    xq = Xq(q, :);
    patch_ids = patch_ids_per_query{q};
    center_dist = sqrt(sum((centers - xq).^2, 2));

    alpha = puPatchWeight(center_dist(patch_ids) ./ radius);
    alpha_sum = sum(alpha);
    if alpha_sum <= 1.0e-14
        alpha = ones(size(alpha));
        alpha_sum = sum(alpha);
    end
    alpha = alpha / alpha_sum;

    row = zeros(1, num_all);
    for k = 1:numel(patch_ids)
        p = patch_ids(k);
        ids = node_ids{p};
        if theta == 2
            stencil = cachedStencils{p};
        else
            stencil = kp.rbffd.RBFStencil();
            stencil.InitializeGeometry(Xnodes(ids, :), sp);
        end
        wloc = localOperatorWeights(stencil, ids, xq, sp, opName, normals(q, :), neuCoeff(q), dirCoeff(q));
        row(ids) = row(ids) + alpha(k) * wloc;
    end

    nz = find(abs(row) > 0);
    rows{q} = q * ones(numel(nz), 1);
    cols{q} = nz(:);
    vals{q} = row(nz(:)).';
end

A = sparse(vertcat(rows{:}), vertcat(cols{:}), vertcat(vals{:}), num_targets, num_all);
end

function patch_ids_per_query = queryPatchIds(patchData, Xq, radius)
tree = patchData.center_tree;
centers = patchData.centers;
if tree.HasSearcher
    patch_ids_per_query = rangesearch(tree.Searcher, Xq, radius);
else
    D = kp.geometry.distanceMatrix(Xq, centers);
    patch_ids_per_query = cell(size(Xq, 1), 1);
    for q = 1:size(Xq, 1)
        patch_ids_per_query{q} = find(D(q, :) < radius);
    end
end

for q = 1:numel(patch_ids_per_query)
    if isempty(patch_ids_per_query{q})
        if isempty(centers)
            continue;
        end
        d = sqrt(sum((centers - Xq(q, :)).^2, 2));
        [~, nearest_patch] = min(d);
        patch_ids_per_query{q} = nearest_patch;
    else
        patch_ids_per_query{q} = patch_ids_per_query{q}(:).';
    end
end
end

function w = localOperatorWeights(stencil, ids, xq, sp, opName, nr, neuCoeff, dirCoeff)
Xloc = stencil.getStencilNodes();
xc = (xq - stencil.getCentroid()) / stencil.getWidth();
r = kp.geometry.distanceMatrix(xq, Xloc);
op = struct('nosolve', false, 'selectdim', 0);

switch lower(string(opName))
    case {"interp", "interpolation"}
        Bpoly = stencil.basis.evaluate(xc, zeros(1, size(Xloc, 2)), true).';
        B = [phsRbf(r, sp.spline_degree).'; Bpoly];
    case {"lap", "laplacian"}
        B = stencil.LapOp(sp, op, r, xq, Xloc, xc, stencil.getScaledStencilNodes());
    case {"bc", "boundary"}
        B = stencil.BCOp(sp, op, neuCoeff, dirCoeff, r, xq, Xloc, xc, stencil.getScaledStencilNodes(), nr);
    otherwise
        error('kp:solvers:BadLocalOperator', 'Unknown local operator "%s".', string(opName));
end

Wfull = stableSolve(stencil.solve_lhs, B);
w = Wfull(1:numel(ids), :).';
end

function Phi = phsRbf(r, degree)
Phi = r .^ degree;
if mod(degree, 2) == 0
    Phi = Phi .* log(r + 2e-16);
end
end

function X = stableSolve(A, B)
warnNear = warning('query', 'MATLAB:nearlySingularMatrix');
warnSing = warning('query', 'MATLAB:singularMatrix');
cleanupObj = onCleanup(@() restoreWarnings(warnNear, warnSing));
warning('off', 'MATLAB:nearlySingularMatrix');
warning('off', 'MATLAB:singularMatrix');
X = A \ B;
if any(~isfinite(X), 'all')
    X = pinv(A) * B;
end
X(~isfinite(X)) = 0;
end

function restoreWarnings(warnNear, warnSing)
warning(warnNear.state, 'MATLAB:nearlySingularMatrix');
warning(warnSing.state, 'MATLAB:singularMatrix');
end

function w = puPatchWeight(r)
w = zeros(size(r));
mask = r < 1;
t = 1 - r(mask);
rm = r(mask);
w(mask) = t.^8 .* (32 * rm.^3 + 25 * rm.^2 + 8 * rm + 1);
end
