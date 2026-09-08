function rbffd_checks()
%RBFFD_CHECKS Lightweight checks for stencil and assembler classes.

[Xg, Yg] = ndgrid(linspace(-1, 1, 5), linspace(-1, 1, 5));
X = [Xg(:), Yg(:)];
interiorMask = abs(X(:, 1)) < 0.999 & abs(X(:, 2)) < 0.999;
activeRows = find(interiorMask);

domain = kp.domain.DomainDescriptor();
domain.setNodes(X, zeros(0, 2), zeros(0, 2));
domain.setSepRad(0.5);
domain.buildStructs();
assert(domain.getIntBdryTree().HasSearcher || ~isempty(domain.getIntBdryTree().Points), ...
    'DomainDescriptor should build a usable neighbor-search structure for interior-boundary nodes.');
[nnIdx, nnDist] = domain.queryKnn("interior_boundary", X(1, :), 4);
assert(size(nnIdx, 2) == 4 && nnDist(1, 1) < 1e-14, ...
    'DomainDescriptor knn queries should return the requested number of neighbors with self distance first.');

sp = kp.rbffd.StencilProperties('n', 9, 'dim', 2, 'ell', 2, 'spline_degree', 3, ...
    'treeMode', 'interior_boundary', 'pointSet', 'interior_boundary');
op = kp.rbffd.OpProperties('recordStencils', true);

specFromAccuracy = kp.rbffd.StencilProperties.fromAccuracy( ...
    'Operator', 'lap', ...
    'ConvergenceOrder', 3, ...
    'Dimension', 2, ...
    'Approximation', 'rbf', ...
    'treeMode', 'interior_boundary', ...
    'pointSet', 'interior_boundary');
assert(specFromAccuracy.ell == 4, 'fromAccuracy should choose ell = p + q - 1.');
assert(specFromAccuracy.npoly == size(kp.poly.total_degree_indices(2, 4), 1), ...
    'fromAccuracy should choose a total-degree polynomial space consistent with ell.');
assert(specFromAccuracy.n >= specFromAccuracy.npoly + 1, ...
    'fromAccuracy should choose a stencil size larger than the polynomial basis.');

f = X(:, 1).^2 + X(:, 2).^2;
g = 2 * X(:, 1) - 3 * X(:, 2) + 1;
gx_exact = 2 * ones(size(X, 1), 1);
gy_exact = -3 * ones(size(X, 1), 1);

fd_rbf = kp.rbffd.FDDiffOp(@() kp.rbffd.RBFStencil());
fd_rbf.AssembleOp(domain, 'lap', sp, op, 'ActiveRows', activeRows);
Lrbf = fd_rbf.getOp();
lap_rbf = Lrbf * f;
assert(all(abs(lap_rbf(activeRows) - 4) < 5e-1), 'RBF FDDiffOp should reproduce the Laplacian of a quadratic on interior grid points.');

fd_wls = kp.rbffd.FDDiffOp(@() kp.rbffd.WeightedLeastSquaresStencil());
fd_wls.AssembleOp(domain, 'lap', sp, op, 'ActiveRows', activeRows);
Lwls = fd_wls.getOp();
lap_wls = Lwls * f;
assert(all(abs(lap_wls(activeRows) - 4) < 1e-10), 'Weighted least-squares FDDiffOp should reproduce quadratic Laplacians exactly.');

gradxOp = kp.rbffd.OpProperties('recordStencils', true, 'selectdim', 0);
gradyOp = kp.rbffd.OpProperties('recordStencils', true, 'selectdim', 1);

fd_gradx_rbf = kp.rbffd.FDDiffOp(@() kp.rbffd.RBFStencil());
fd_gradx_rbf.AssembleOp(domain, 'grad', sp, gradxOp, 'ActiveRows', activeRows);
Gx_rbf = fd_gradx_rbf.getOp();
gx_rbf = Gx_rbf * g;
assert(all(abs(gx_rbf(activeRows) - gx_exact(activeRows)) < 5e-1), ...
    'RBF FDDiffOp should reproduce x-derivatives of linear functions on interior grid points.');

fd_grady_rbf = kp.rbffd.FDDiffOp(@() kp.rbffd.RBFStencil());
fd_grady_rbf.AssembleOp(domain, 'grad', sp, gradyOp, 'ActiveRows', activeRows);
Gy_rbf = fd_grady_rbf.getOp();
gy_rbf = Gy_rbf * g;
assert(all(abs(gy_rbf(activeRows) - gy_exact(activeRows)) < 5e-1), ...
    'RBF FDDiffOp should reproduce y-derivatives of linear functions on interior grid points.');

fd_gradx_wls = kp.rbffd.FDDiffOp(@() kp.rbffd.WeightedLeastSquaresStencil());
fd_gradx_wls.AssembleOp(domain, 'grad', sp, gradxOp, 'ActiveRows', activeRows);
Gx_wls = fd_gradx_wls.getOp();
gx_wls = Gx_wls * g;
assert(all(abs(gx_wls(activeRows) - gx_exact(activeRows)) < 1e-10), ...
    'Weighted least-squares FDDiffOp should reproduce x-derivatives of linear functions exactly.');

fd_grady_wls = kp.rbffd.FDDiffOp(@() kp.rbffd.WeightedLeastSquaresStencil());
fd_grady_wls.AssembleOp(domain, 'grad', sp, gradyOp, 'ActiveRows', activeRows);
Gy_wls = fd_grady_wls.getOp();
gy_wls = Gy_wls * g;
assert(all(abs(gy_wls(activeRows) - gy_exact(activeRows)) < 1e-10), ...
    'Weighted least-squares FDDiffOp should reproduce y-derivatives of linear functions exactly.');

fdo = kp.rbffd.FDODiffOp(@() kp.rbffd.RBFStencil());
fdo.AssembleOp(domain, 'lap', sp, kp.rbffd.OpProperties('recordStencils', true, 'OverlapLoad', 0.5), 'ActiveRows', activeRows);
Lfdo = fdo.getOp();
coveredRows = unique(find(any(Lfdo ~= 0, 2)));
assert(isequal(coveredRows(:), activeRows(:)), 'FDODiffOp should cover each requested active row exactly once.');
lap_fdo = Lfdo * f;
assert(all(abs(lap_fdo(activeRows) - 4) < 5e-1), ...
    'FDODiffOp should reproduce the Laplacian of a quadratic on covered interior rows.');

localStencil = X(abs(X(:, 1)) <= 0.5 & abs(X(:, 2)) <= 0.5, :);
spLocal = kp.rbffd.StencilProperties('n', 9, 'dim', 2, 'ell', 1, 'spline_degree', 3, ...
    'treeMode', 'interior_boundary', 'pointSet', 'interior_boundary');
stencil = kp.rbffd.WeightedLeastSquaresStencil();
W = stencil.ComputeWeights(localStencil, spLocal, op, 'interp', 1);
assert(size(W, 1) == 9 && size(W, 2) == 1, 'WeightedLeastSquaresStencil should return one weight vector per rhs point.');

leg = kp.rbffd.RBFStencil();
Wi = leg.ComputeWeights(localStencil, spLocal, op, 'interp', 1);
assert(size(Wi, 1) == 9 && size(Wi, 2) == 1, 'RBFStencil should return one interpolation weight vector per rhs point.');

centerIndex = ceil(size(localStencil, 1) / 2);
centerPoint = localStencil(centerIndex, :);
spInterp = kp.rbffd.StencilProperties('n', size(localStencil, 1), 'dim', 2, 'ell', 2, 'spline_degree', 3, ...
    'treeMode', 'interior_boundary', 'pointSet', 'interior_boundary');
interpStencilWls = kp.rbffd.WeightedLeastSquaresStencil();
WinterpWls = interpStencilWls.ComputeWeights(localStencil, spInterp, op, 'interp', centerIndex);
fLocal = localStencil(:, 1).^2 + localStencil(:, 2).^2;
assert(abs(WinterpWls(:, 1).' * fLocal - sum(centerPoint.^2)) < 1e-10, ...
    'Weighted least-squares interpolation weights should reproduce quadratic data at the query point.');

interpStencilRbf = kp.rbffd.RBFStencil();
WinterpRbf = interpStencilRbf.ComputeWeights(localStencil, spInterp, op, 'interp', centerIndex);
assert(abs(WinterpRbf(:, 1).' * fLocal - sum(centerPoint.^2)) < 5e-1, ...
    'RBF interpolation weights should reproduce quadratic data at the query point within a modest tolerance.');

t = linspace(0, 2*pi, 60).';
t(end) = [];
curve = [cos(t), 0.7 * sin(t)];
surface = kp.geometry.EmbeddedSurface();
surface.setDataSites(curve);
surface.buildClosedGeometricModelPS(2, 0.05, size(curve, 1));
surface.buildLevelSetFromGeometricModel([]);
generator = kp.nodes.DomainNodeGenerator();
geomDomain = generator.buildDomainDescriptorFromGeometry(surface, 0.08, ...
    'Seed', 29, 'StripCount', 5);
nb = geomDomain.getNumBdryNodes();
nf = geomDomain.getNumTotalNodes();
bcSp = kp.rbffd.StencilProperties.fromAccuracy( ...
    'Operator', 'bc', ...
    'ConvergenceOrder', 3, ...
    'Dimension', 2, ...
    'Approximation', 'rbf', ...
    'treeMode', 'all', ...
    'pointSet', 'boundary');
bcAssembler = kp.rbffd.FDDiffOp(@() kp.rbffd.RBFStencil());
bcAssembler.AssembleOp(geomDomain, 'bc', bcSp, op, ...
    'NeuCoeff', zeros(nb, 1), ...
    'DirCoeff', ones(nb, 1));
BC = bcAssembler.getOp();
assert(isequal(size(BC), [nb, nf]), ...
    'Boundary operators should have one row per boundary point and one column per total node.');

disp('rbffd checks passed');
end
