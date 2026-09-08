function A = puLocalizedOperator(domain, patchData, xi, Xq, opName, varargin)
%PULOCALIZEDOPERATOR Assemble a PU-collocation operator on arbitrary targets.

if isempty(Xq)
    A = sparse(0, domain.getNumTotalNodes());
    return;
end

Xnodes = domain.getAllNodes();
node_ids = patchData.node_ids;
cachedStencils = patchData.stencils;
dim = size(Xnodes, 2);
num_targets = size(Xq, 1);
num_all = size(Xnodes, 1);
[~, alphaMat] = kp.solvers.puQueryPatchWeights(patchData, Xq);

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
patchStencilCache = cell(size(node_ids));
for p = 1:numel(node_ids)
    qIdx = find(alphaMat(:, p));
    if isempty(qIdx)
        continue;
    end
    ids = node_ids{p};
    if theta == 2
        stencil = cachedStencils{p};
    else
        stencil = patchStencilCache{p};
        if isempty(stencil)
            stencil = kp.rbffd.RBFStencil();
            stencil.InitializeGeometry(Xnodes(ids, :), sp);
            patchStencilCache{p} = stencil;
        end
    end

    wloc = localOperatorWeightsBatch(stencil, ids, Xq(qIdx, :), sp, opName, normals(qIdx, :), neuCoeff(qIdx), dirCoeff(qIdx));
    weighted = wloc .* full(alphaMat(qIdx, p));
    [ii, jj, vv] = find(weighted);
    for nzIdx = 1:numel(vv)
        q = qIdx(ii(nzIdx));
        rows{q}(end + 1, 1) = q;
        cols{q}(end + 1, 1) = ids(jj(nzIdx));
        vals{q}(end + 1, 1) = vv(nzIdx);
    end
end

A = sparse(vertcat(rows{:}), vertcat(cols{:}), vertcat(vals{:}), num_targets, num_all);
end
function w = localOperatorWeightsBatch(stencil, ids, Xq, sp, opName, nr, neuCoeff, dirCoeff)
Xloc = stencil.getStencilNodes();
xc = (Xq - stencil.getCentroid()) / stencil.getWidth();
r = kp.geometry.distanceMatrix(Xq, Xloc);
op = struct('nosolve', false, 'selectdim', 0);

switch lower(string(opName))
    case {"interp", "interpolation"}
        Bpoly = stencil.basis.evaluate(xc, zeros(1, size(Xloc, 2)), true).';
        B = [phsRbf(r, sp.spline_degree).'; Bpoly];
    case {"lap", "laplacian"}
        B = stencil.LapOp(sp, op, r, Xq, Xloc, xc, stencil.getScaledStencilNodes());
    case {"bc", "boundary"}
        w = zeros(size(Xq, 1), numel(ids));
        for q = 1:size(Xq, 1)
            B = stencil.BCOp(sp, op, neuCoeff(q), dirCoeff(q), r(q, :), Xq(q, :), Xloc, xc(q, :), stencil.getScaledStencilNodes(), nr(q, :));
            Wfull = stableSolve(stencil.solve_lhs, B);
            w(q, :) = Wfull(1:numel(ids), :).';
        end
        return;
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
