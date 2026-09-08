function manifold_tangent_plane_checks()
%MANIFOLD_TANGENT_PLANE_CHECKS Targeted checks for tangent-plane manifold RBF-FD.

N = 642;
X = kp.geometry.fibonacciSphere(N);
X = X ./ vecnorm(X, 2, 2);
nr = X;
tree = KDTreeSearcher(X);
op = kp.manifold.rbffdop(2, 4, 2, 0);

[L, Gx, Gy, Gz] = kp.manifold.FormSurfaceDiffOpsTP( ...
    X, op.rbf, op.drbfor, op.d2rbf, nr, tree, op.stencilSize, op.ell);

x = X(:, 1);
y = X(:, 2);
z = X(:, 3);

lap_x_rel = relerr(L * x, -2 * x);
lap_y_rel = relerr(L * y, -2 * y);
lap_z_rel = relerr(L * z, -2 * z);

grad_x_true = [ones(N, 1) - x.^2, -x .* y, -x .* z];
grad_y_true = [-x .* y, ones(N, 1) - y.^2, -y .* z];
grad_z_true = [-x .* z, -y .* z, ones(N, 1) - z.^2];

grad_x_num = [Gx * x, Gy * x, Gz * x];
grad_y_num = [Gx * y, Gy * y, Gz * y];
grad_z_num = [Gx * z, Gy * z, Gz * z];

grad_x_rel = relerr(grad_x_num, grad_x_true);
grad_y_rel = relerr(grad_y_num, grad_y_true);
grad_z_rel = relerr(grad_z_num, grad_z_true);

assert(isfinite(lap_x_rel) && isfinite(lap_y_rel) && isfinite(lap_z_rel), ...
    'Surface Laplacian checks returned non-finite values.');
assert(lap_x_rel < 5.0e-2 && lap_y_rel < 5.0e-2 && lap_z_rel < 5.0e-2, ...
    'Surface Laplacian tangent-plane errors are too large.');
assert(grad_x_rel < 8.0e-2 && grad_y_rel < 8.0e-2 && grad_z_rel < 8.0e-2, ...
    'Surface gradient tangent-plane errors are too large.');

% Local coordinate normalization must make the assembled operators covariant
% under a uniform change of physical length units.
op6 = kp.manifold.rbffdop(2, 6, 2, 0);
neighborIds = knnsearch(tree, X, 'K', op6.stencilSize);
[L1, Gx1, Gy1, Gz1] = kp.manifold.FormSurfaceDiffOpsTP( ...
    X, op6.rbf, op6.drbfor, op6.d2rbf, nr, tree, ...
    op6.stencilSize, op6.ell, neighborIds);
scale = 0.1;
scaledX = scale * X;
scaledTree = KDTreeSearcher(scaledX);
[Ls, Gxs, Gys, Gzs] = kp.manifold.FormSurfaceDiffOpsTP( ...
    scaledX, op6.rbf, op6.drbfor, op6.d2rbf, nr, scaledTree, ...
    op6.stencilSize, op6.ell, neighborIds);
lapScaleError = relerr(scale^2 * Ls, L1);
gradScaleError = max([relerr(scale * Gxs, Gx1), ...
    relerr(scale * Gys, Gy1), relerr(scale * Gzs, Gz1)]);
assert(lapScaleError < 1.0e-9 && gradScaleError < 1.0e-9, ...
    'Tangent-plane RBF-FD operators are not scale covariant.');

disp('manifold tangent plane checks passed');
fprintf('  lap rel errors: [%.3e, %.3e, %.3e]\n', lap_x_rel, lap_y_rel, lap_z_rel);
fprintf('  grad rel errors: [%.3e, %.3e, %.3e]\n', grad_x_rel, grad_y_rel, grad_z_rel);
fprintf('  scale covariance errors (L / G): %.3e / %.3e\n', ...
    lapScaleError, gradScaleError);
end

function r = relerr(u, v)
r = norm(u(:) - v(:)) / max(norm(v(:)), 1.0e-14);
end
