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

sp = kp.rbffd.StencilProperties('n', 9, 'dim', 2, 'ell', 2, 'spline_degree', 3, ...
    'treeMode', 'interior_boundary', 'pointSet', 'interior_boundary');
op = kp.rbffd.OpProperties('recordStencils', true);

f = X(:, 1).^2 + X(:, 2).^2;

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

fdo = kp.rbffd.FDODiffOp(@() kp.rbffd.RBFStencil());
fdo.AssembleOp(domain, 'lap', sp, kp.rbffd.OpProperties('recordStencils', true, 'OverlapLoad', 0.5), 'ActiveRows', activeRows);
Lfdo = fdo.getOp();
coveredRows = unique(find(any(Lfdo ~= 0, 2)));
assert(isequal(coveredRows(:), activeRows(:)), 'FDODiffOp should cover each requested active row exactly once.');

localStencil = X(abs(X(:, 1)) <= 0.5 & abs(X(:, 2)) <= 0.5, :);
spLocal = kp.rbffd.StencilProperties('n', 9, 'dim', 2, 'ell', 1, 'spline_degree', 3, ...
    'treeMode', 'interior_boundary', 'pointSet', 'interior_boundary');
stencil = kp.rbffd.WeightedLeastSquaresStencil();
W = stencil.ComputeWeights(localStencil, spLocal, op, 'interp', 1);
assert(size(W, 1) == 9 && size(W, 2) == 1, 'WeightedLeastSquaresStencil should return one weight vector per rhs point.');

leg = kp.rbffd.RBFStencil();
Wi = leg.ComputeWeights(localStencil, spLocal, op, 'interp', 1);
assert(size(Wi, 1) == 9 && size(Wi, 2) == 1, 'RBFStencil should return one interpolation weight vector per rhs point.');

disp('rbffd checks passed');
end
