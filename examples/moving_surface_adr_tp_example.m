function moving_surface_adr_tp_example()
%MOVING_SURFACE_ADR_TP_EXAMPLE Tangent-plane RBF-FD ADR example on a moving sphere.

N = 900;
xi = 4;
mu = 0.1;
dt = 1.0e-2;
T = 0.12;
nsteps = max(1, round(T / dt));
dt = T / nsteps;
hyppow = 3;

op = kp.manifold.rbffdop(2, xi, 2, 0);
[x0, ~, ~] = deformSphereSmoothLobedWarp(N, 0.0, 5, 5);
c0 = initialCondition(x0);

fprintf('Moving surface ADR tangent-plane example\n');
fprintf('  nodes: %d\n', N);
fprintf('  xi: %d\n', xi);
fprintf('  dt: %.4f\n', dt);
fprintf('  steps: %d\n', nsteps);

[x1, ~, nr1] = deformSphereSmoothLobedWarp(N, dt, 5, 5);
v1 = (x1 - x0) / dt;
[L, Gx, Gy, Gz, H, hyp_gamma] = surfaceOpsAndHyperviscosity(x1, nr1, op, hyppow);
I = speye(N, N);
lhs = I - mu * dt * L;
rhs = c0 + dt * reaction(0.0, x0, c0) ...
    - dt * c0 .* (Gx * v1(:, 1) + Gy * v1(:, 2) + Gz * v1(:, 3)) ...
    + dt * hyp_gamma * c0 .* (H * (v1(:, 1) + v1(:, 2) + v1(:, 3)));
c1 = lhs \ rhs;

[x2, ~, nr2] = deformSphereSmoothLobedWarp(N, 2 * dt, 5, 5);
v2 = (3 * x2 - 4 * x1 + x0) / (2 * dt);
[L, Gx, Gy, Gz, H, hyp_gamma] = surfaceOpsAndHyperviscosity(x2, nr2, op, hyppow);
lhs = I - (2 / 3) * mu * dt * L;
c2a = 2 * c1 - c0;
rhs = (4 / 3) * c1 - (1 / 3) * c0 + (2 / 3) * dt * (2 * reaction(dt, x1, c1) - reaction(0.0, x0, c0)) ...
    - (2 / 3) * dt * c2a .* (Gx * v2(:, 1) + Gy * v2(:, 2) + Gz * v2(:, 3)) ...
    + (2 / 3) * dt * hyp_gamma * c2a .* (H * (v2(:, 1) + v2(:, 2) + v2(:, 3)));
c2 = lhs \ rhs;
c3 = c2;

for step = 3:nsteps
    [x3, ~, nr3] = deformSphereSmoothLobedWarp(N, step * dt, 5, 5);
    v3 = (11 * x3 - 18 * x2 + 9 * x1 - 2 * x0) / (6 * dt);
    [L, Gx, Gy, Gz, H, hyp_gamma] = surfaceOpsAndHyperviscosity(x3, nr3, op, hyppow);
    lhs = I - (6 / 11) * mu * dt * L;

    c3a = 3 * c2 - 3 * c1 + c0;
    source = -(6 / 11) * dt * c3a .* (Gx * v3(:, 1) + Gy * v3(:, 2) + Gz * v3(:, 3)) ...
        + (6 / 11) * dt * hyp_gamma * c3a .* (H * (v3(:, 1) + v3(:, 2) + v3(:, 3)));
    rhs = (18 / 11) * c2 - (9 / 11) * c1 + (2 / 11) * c0 ...
        + (6 / 11) * dt * (3 * reaction(step * dt, x2, c2) ...
        - 3 * reaction((step - 1) * dt, x1, c1) + reaction((step - 2) * dt, x0, c0)) ...
        + source;
    c3 = lhs \ rhs;

    c0 = c1; c1 = c2; c2 = c3;
    x0 = x1; x1 = x2; x2 = x3;
end

finalField = c3;
tri = kp.geometry.triangulateClosedSurface(x2);

figure('Color', 'w', 'Position', [100, 100, 760, 620]);
trisurf(tri, x2(:, 1), x2(:, 2), x2(:, 3), finalField, 'EdgeColor', 'none');
camlight headlight;
lighting gouraud;
material dull;
shading interp;
daspect([1, 1, 1]);
view([-51.4487, 66]);
axis off;
title(sprintf('Moving surface ADR at T = %.2f', T));
colorbar;

fprintf('  final field min/max: %.6e / %.6e\n', min(finalField), max(finalField));
end

function [L, Gx, Gy, Gz, H, hyp_gamma] = surfaceOpsAndHyperviscosity(X, nr, op, hyppow)
tree = KDTreeSearcher(X);
[L, Gx, Gy, Gz] = kp.manifold.FormSurfaceDiffOpsTP( ...
    X, op.rbf, op.drbfor, op.d2rbf, nr, tree, op.stencilSize, op.ell);
assertNegativeLargestReal(L);

H = L^hyppow;
hx = sqrt(1 / size(X, 1));
hyp_gamma = kp.manifold.hyperviscosityCoefficient(H, Gx, Gy, Gz, X, nr, hyppow, hx);
end

function assertNegativeLargestReal(L)
tol = 1.0e-2;
[~, sig] = eigs(L, 1, 'largestreal', 'Tolerance', tol);
if sig > tol
    error('kp:manifold:BadSpectrum', ...
        ['+ve real eigenvalue detected, smooth normals with more nearest ' ...
         'neighbors or check points.']);
end
end

function c = initialCondition(X)
c = 0.45 + 0.20 * sin(pi * X(:, 1)) .* cos(pi * X(:, 2)) .* sin(pi * X(:, 3));
c = max(c, 0);
end

function out = reaction(~, ~, u)
out = u .* (1 - u) .* (u - 0.3);
end

function [X, A_est, N] = deformSphereSmoothLobedWarp(N_pts, t, l, m)
idx = (0:N_pts-1).';
golden = pi * (3 - sqrt(5));
theta = acos(1 - 2 * (idx + 0.5) / N_pts);
phi = mod(idx * golden, 2 * pi);

u = [sin(theta) .* cos(phi), sin(theta) .* sin(phi), cos(theta)];

eps_v = 0.4;
omega_v = 1;
S = sin(theta).^l;
Cphi = cos(m * phi);
Sphi = sin(m * phi);
amp_v = eps_v * sin(omega_v * t);
r = 1 + amp_v .* (S .* Cphi);
dr_dtheta = amp_v .* (l * sin(theta).^(l - 1) .* cos(theta) .* Cphi);
dr_dphi = amp_v .* (S .* (-m * Sphi));

e_theta = [cos(theta) .* cos(phi), cos(theta) .* sin(phi), -sin(theta)];
e_phi = [-sin(theta) .* sin(phi), sin(theta) .* cos(phi), zeros(N_pts, 1)];

eps_t = 0.4;
omega_t = 1;
amp_t = eps_t * sin(omega_t * t);
k_t = 3;

f = sin(k_t * phi) .* sin(theta);
g = cos(k_t * phi) .* sin(theta) .* cos(theta);
f_theta = sin(k_t * phi) .* cos(theta);
f_phi = k_t * cos(k_t * phi) .* sin(theta);
g_theta = cos(k_t * phi) .* (cos(theta).^2 - sin(theta).^2);
g_phi = -k_t * sin(k_t * phi) .* sin(theta) .* cos(theta);

e_theta_theta = [-sin(theta) .* cos(phi), -sin(theta) .* sin(phi), -cos(theta)];
e_theta_phi = [-cos(theta) .* sin(phi), cos(theta) .* cos(phi), zeros(N_pts, 1)];
e_phi_theta = [-cos(theta) .* sin(phi), cos(theta) .* cos(phi), zeros(N_pts, 1)];
e_phi_phi = [-sin(theta) .* cos(phi), -sin(theta) .* sin(phi), zeros(N_pts, 1)];

w = f .* e_theta + g .* e_phi;
w_theta = f_theta .* e_theta + f .* e_theta_theta + g_theta .* e_phi + g .* e_phi_theta;
w_phi = f_phi .* e_theta + f .* e_theta_phi + g_phi .* e_phi + g .* e_phi_phi;

X_rad = [r .* u(:, 1), r .* u(:, 2), r .* u(:, 3)];
X = X_rad + amp_t * w;

X_theta = dr_dtheta .* u + r .* e_theta + amp_t .* w_theta;
X_phi = dr_dphi .* u + r .* e_phi + amp_t .* w_phi;
N_raw = cross(X_theta, X_phi, 2);
N = N_raw ./ vecnorm(N_raw, 2, 2);
R = vecnorm(X, 2, 2);
A_est = (4 * pi / N_pts) * sum(R.^2);
end
