function manifold_mean_curvature_flow_checks()
%MANIFOLD_MEAN_CURVATURE_FLOW_CHECKS Focused checks for MCF geometry update.

N = 256;
R0 = 1.4;
dt = 1.0e-3;
U = kp.geometry.fibonacciSphere(N);
U = U ./ vecnorm(U, 2, 2);
X0 = R0 .* U;

solver = kp.manifold.MeanCurvatureFlowSolver( ...
    'Mode', "explicit", ...
    'ProjectNormal', true, ...
    'NormalProvider', @(X, ~, ~) kp.geometry.normalizeRows(X));
solver.init(X0, U, 4, dt);

V = solver.meanCurvatureVelocity();
Vex = -(2 / R0^2) .* X0;
relVelErr = norm(V - Vex) / norm(Vex);
assert(relVelErr < 5.0e-3, 'MCF coordinate Laplacian velocity is inaccurate.');

solver.step();
Rex = sqrt(R0^2 - 4 * dt);
r = vecnorm(solver.nodes(), 2, 2);
relRadiusErr = abs(mean(r) - Rex) / Rex;
assert(relRadiusErr < 2.0e-3, 'One MCF step gives the wrong sphere radius.');

fprintf('manifold mean-curvature-flow checks passed\n');
fprintf('  velocity relative error: %.3e\n', relVelErr);
fprintf('  one-step radius relative error: %.3e\n', relRadiusErr);
end
