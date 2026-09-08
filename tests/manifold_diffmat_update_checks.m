function manifold_diffmat_update_checks()
%MANIFOLD_DIFFMAT_UPDATE_CHECKS Defect-corrected TP matrix update checks.

N = 128;
xi = 4;
U = kp.geometry.fibonacciSphere(N);
U = U ./ vecnorm(U, 2, 2);
X0 = U;
X1 = 1.00001 .* U;
nr = U;
op = kp.manifold.rbffdop(2, xi, 2, 0);

updater = kp.manifold.TangentPlaneDiffMatUpdater(xi, ...
    'UseParallel', false, ...
    'DefectTolerance', 1.0e-7, ...
    'MaxDefectIterations', 4);
[L0u, Gx0u, Gy0u, Gz0u, stats0] = updater.assemble(X0, nr);

tree0 = KDTreeSearcher(X0);
[L0d, Gx0d, Gy0d, Gz0d] = kp.manifold.FormSurfaceDiffOpsTP( ...
    X0, op.rbf, op.drbfor, op.d2rbf, nr, tree0, op.stencilSize, op.ell);
initialRel = combinedRelativeDifference( ...
    {L0u, Gx0u, Gy0u, Gz0u}, {L0d, Gx0d, Gy0d, Gz0d});
assert(initialRel < 1.0e-10, 'Initial updater assembly must match direct assembly.');
assert(stats0.direct == N, 'Initial updater assembly should directly factor all stencils.');

[L1u, Gx1u, Gy1u, Gz1u, stats1] = updater.assemble(X1, nr);
tree1 = KDTreeSearcher(X1);
[L1d, Gx1d, Gy1d, Gz1d] = kp.manifold.FormSurfaceDiffOpsTP( ...
    X1, op.rbf, op.drbfor, op.d2rbf, nr, tree1, op.stencilSize, op.ell);
updateRel = combinedRelativeDifference( ...
    {L1u, Gx1u, Gy1u, Gz1u}, {L1d, Gx1d, Gy1d, Gz1d});
assert(updateRel < 2.0e-6, 'Defect-corrected update must match direct assembly.');
assert(stats1.defectCorrected > 0, 'Moved surface should use defect correction.');
assert(stats1.neighborReuses == 1 && stats1.neighborSearches == 0, ...
    'Small bounded motion should reuse cached neighbor indices.');

fprintf('manifold diff-mat update checks passed\n');
fprintf('  initial relative matrix difference: %.3e\n', initialRel);
fprintf('  update relative matrix difference: %.3e\n', updateRel);
fprintf('  update direct/defect/fallback: %d / %d / %d\n', ...
    stats1.direct, stats1.defectCorrected, stats1.defectFailedRefactored);
fprintf('  neighbor searches/reuses/rejections: %d / %d / %d\n', ...
    stats1.neighborSearches, stats1.neighborReuses, stats1.neighborReuseRejected);
fprintf('  update max residual: %.3e\n', stats1.maxRelativeResidual);
end

function rel = combinedRelativeDifference(A, B)
num = 0;
den = 0;
for k = 1:numel(A)
    num = num + norm(A{k} - B{k}, 'fro')^2;
    den = den + norm(B{k}, 'fro')^2;
end
rel = sqrt(num / max(den, eps));
end
