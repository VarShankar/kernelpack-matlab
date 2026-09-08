function results = moving_surface_adr_tp_convergence_study(xiVals, Nvals, dtScale)
%MOVING_SURFACE_ADR_TP_CONVERGENCE_STUDY Manufactured convergence on a moving surface.

if nargin < 1 || isempty(xiVals)
    xiVals = [2, 4, 6];
end
if nargin < 2 || isempty(Nvals)
    Nvals = [256, 576, 1024, 1600];
end
if nargin < 3 || isempty(dtScale)
    dtScale = 0.5;
end

mu = 0.1;
T = 0.12;
hyppow = 3;
l = 5;
m = 5;

Nvals = Nvals(:).';
xiVals = xiVals(:).';
numLevels = numel(Nvals);
numXi = numel(xiVals);
results = repmat(struct( ...
    'xi', 0, ...
    'N', zeros(1, numLevels), ...
    'h', zeros(1, numLevels), ...
    'dt', zeros(1, numLevels), ...
    'nsteps', zeros(1, numLevels), ...
    'relerr', zeros(1, numLevels), ...
    'rate', nan(1, numLevels)), 1, numXi);

for ix = 1:numXi
    xi = xiVals(ix);
    results(ix).xi = xi;
    for level = 1:numLevels
        N = Nvals(level);
        h = sqrt(4 * pi / N);
        dtTarget = dtScale * h^(xi / 3);
        nsteps = max(3, ceil(T / dtTarget));
        dt = T / nsteps;

        fprintf('Moving surface ADR TP convergence\n');
        fprintf('  xi: %d\n', xi);
        fprintf('  level %d / %d\n', level, numLevels);
        fprintf('  N: %d\n', N);
        fprintf('  h: %.6e\n', h);
        fprintf('  dt: %.6e\n', dt);
        fprintf('  steps: %d\n', nsteps);

        err = runOneLevel(N, xi, mu, T, dt, nsteps, hyppow, l, m);
        results(ix).N(level) = N;
        results(ix).h(level) = h;
        results(ix).dt(level) = dt;
        results(ix).nsteps(level) = nsteps;
        results(ix).relerr(level) = err;
        fprintf('  relative L2 error: %.6e\n', err);
    end

    for k = 2:numLevels
        results(ix).rate(k) = log(results(ix).relerr(k - 1) / results(ix).relerr(k)) / ...
            log(results(ix).h(k - 1) / results(ix).h(k));
    end
end

fig = figure('Color', 'w', 'Position', [100, 100, 760, 560]);
hold on;
colors = lines(numXi);
for ix = 1:numXi
    plot(results(ix).h, results(ix).relerr, '-o', 'LineWidth', 1.6, ...
        'MarkerSize', 7, 'Color', colors(ix, :), ...
        'DisplayName', sprintf('\\xi = %d', results(ix).xi));
end
hold off;
set(gca, 'XScale', 'log', 'YScale', 'log');
grid on;
xlabel('h');
ylabel('Relative L2 error');
title('Moving surface ADR tangent-plane convergence');
legend('Location', 'best');
set(gca, 'XDir', 'reverse');

imgPath = fullfile(pwd, 'docs', 'figures', 'moving_surface_adr_tp_convergence.png');
exportgraphics(fig, imgPath, 'Resolution', 180);

save(fullfile(pwd, 'moving_surface_adr_tp_convergence_results.mat'), 'results');
end

function relerr = runOneLevel(N, xi, mu, T, dt, nsteps, hyppow, l, m)
op = kp.manifold.rbffdop(2, xi, 2, 0);

[x0, ~, theta, phi] = movingSurfaceGeometry(N, 0.0, l, m);
c0 = exactField(0.0, theta, phi);

[x1, nr1] = movingSurfaceGeometry(N, dt, l, m);
v1 = (x1 - x0) / dt;
[L, Gx, Gy, Gz, H, hyp_gamma] = surfaceOpsAndHyperviscosity(x1, nr1, op, hyppow);
I = speye(N, N);
lhs = I - mu * dt * L;
rhs = c0 + dt * exactReaction(dt, theta, phi, mu, l, m) ...
    - dt * c0 .* (Gx * v1(:, 1) + Gy * v1(:, 2) + Gz * v1(:, 3)) ...
    + dt * hyp_gamma * c0 .* (H * (v1(:, 1) + v1(:, 2) + v1(:, 3)));
c1 = lhs \ rhs;

[x2, nr2] = movingSurfaceGeometry(N, 2 * dt, l, m);
v2 = (3 * x2 - 4 * x1 + x0) / (2 * dt);
[L, Gx, Gy, Gz, H, hyp_gamma] = surfaceOpsAndHyperviscosity(x2, nr2, op, hyppow);
lhs = I - (2 / 3) * mu * dt * L;
c2a = 2 * c1 - c0;
rhs = (4 / 3) * c1 - (1 / 3) * c0 + (2 / 3) * dt * ...
    exactReaction(2 * dt, theta, phi, mu, l, m) ...
    - (2 / 3) * dt * c2a .* (Gx * v2(:, 1) + Gy * v2(:, 2) + Gz * v2(:, 3)) ...
    + (2 / 3) * dt * hyp_gamma * c2a .* (H * (v2(:, 1) + v2(:, 2) + v2(:, 3)));
c2 = lhs \ rhs;
c3 = c2;

for step = 3:nsteps
    tnow = step * dt;
    [x3, nr3] = movingSurfaceGeometry(N, tnow, l, m);
    v3 = (11 * x3 - 18 * x2 + 9 * x1 - 2 * x0) / (6 * dt);
    [L, Gx, Gy, Gz, H, hyp_gamma] = surfaceOpsAndHyperviscosity(x3, nr3, op, hyppow);
    lhs = I - (6 / 11) * mu * dt * L;

    c3a = 3 * c2 - 3 * c1 + c0;
    rhs = (18 / 11) * c2 - (9 / 11) * c1 + (2 / 11) * c0 ...
        + (6 / 11) * dt * exactReaction(step * dt, theta, phi, mu, l, m) ...
        - (6 / 11) * dt * c3a .* (Gx * v3(:, 1) + Gy * v3(:, 2) + Gz * v3(:, 3)) ...
        + (6 / 11) * dt * hyp_gamma * c3a .* (H * (v3(:, 1) + v3(:, 2) + v3(:, 3)));
    c3 = lhs \ rhs;

    c0 = c1; c1 = c2; c2 = c3;
    x0 = x1; x1 = x2; x2 = x3;
end

cex = exactField(T, theta, phi);
relerr = norm(c3 - cex) / max(norm(cex), 1.0e-14);
end

function c = exactField(t, theta, phi)
base = 0.8 + 0.15 * sin(theta) .* cos(phi) + 0.05 * sin(2 * theta) .* sin(2 * phi);
c = exp(-t) .* base;
end

function r = exactReaction(t, theta, phi, mu, l, m)
c = exactField(t, theta, phi);
ct = -c;
divu = surfaceDivergence(theta, phi, t, l, m);
lapc = exp(-t) .* laplaceBeltramiBase(theta, phi, t, l, m);
r = ct + c .* divu - mu .* lapc;
end

function divu = surfaceDivergence(theta, phi, t, l, m)
tau = 1.0e-6;
[~, ~, sqrtgP2] = geometryMetric(theta, phi, t + 2 * tau, l, m);
[~, ~, sqrtgP1] = geometryMetric(theta, phi, t + tau, l, m);
[~, ~, sqrtgM1] = geometryMetric(theta, phi, t - tau, l, m);
[~, ~, sqrtgM2] = geometryMetric(theta, phi, t - 2 * tau, l, m);
[~, ~, sqrtg0] = geometryMetric(theta, phi, t, l, m);
divu = ((-sqrtgP2 + 8 * sqrtgP1 - 8 * sqrtgM1 + sqrtgM2) ./ (12 * tau)) ./ sqrtg0;
end

function lapf = laplaceBeltramiBase(theta, phi, t, l, m)
epsp = 1.0e-6;

[~, ~, ~, sqrtg] = geometryMetric(theta, phi, t, l, m);

AthetaP2 = fluxTheta(theta + 2 * epsp, phi, t, l, m);
AthetaP1 = fluxTheta(theta + epsp, phi, t, l, m);
AthetaM1 = fluxTheta(theta - epsp, phi, t, l, m);
AthetaM2 = fluxTheta(theta - 2 * epsp, phi, t, l, m);
AphiP2 = fluxPhi(theta, phi + 2 * epsp, t, l, m);
AphiP1 = fluxPhi(theta, phi + epsp, t, l, m);
AphiM1 = fluxPhi(theta, phi - epsp, t, l, m);
AphiM2 = fluxPhi(theta, phi - 2 * epsp, t, l, m);

lapf = ((-AthetaP2 + 8 * AthetaP1 - 8 * AthetaM1 + AthetaM2) ./ (12 * epsp) ...
    + (-AphiP2 + 8 * AphiP1 - 8 * AphiM1 + AphiM2) ./ (12 * epsp)) ./ sqrtg;
end

function val = fluxTheta(theta, phi, t, l, m)
[E, F, G, sqrtg] = geometryMetric(theta, phi, t, l, m);
detg = E .* G - F.^2;
[fth, fph] = baseDerivatives(theta, phi);
val = sqrtg .* (G .* fth - F .* fph) ./ detg;
end

function val = fluxPhi(theta, phi, t, l, m)
[E, F, G, sqrtg] = geometryMetric(theta, phi, t, l, m);
detg = E .* G - F.^2;
[fth, fph] = baseDerivatives(theta, phi);
val = sqrtg .* (-F .* fth + E .* fph) ./ detg;
end

function [fth, fph] = baseDerivatives(theta, phi)
fth = 0.15 * cos(theta) .* cos(phi) + 0.10 * cos(2 * theta) .* sin(2 * phi);
fph = -0.15 * sin(theta) .* sin(phi) + 0.10 * sin(2 * theta) .* cos(2 * phi);
end

function [X, nr, theta, phi] = movingSurfaceGeometry(N_pts, t, l, m)
idx = (0:N_pts-1).';
golden = pi * (3 - sqrt(5));
theta = acos(1 - 2 * (idx + 0.5) / N_pts);
phi = mod(idx * golden, 2 * pi);
[X, nr] = surfaceFromMaterial(theta, phi, t, l, m);
end

function [E, F, G, sqrtg] = geometryMetric(theta, phi, t, l, m)
[~, ~, Xtheta, Xphi] = surfaceFromMaterial(theta, phi, t, l, m);
E = sum(Xtheta .* Xtheta, 2);
F = sum(Xtheta .* Xphi, 2);
G = sum(Xphi .* Xphi, 2);
sqrtg = sqrt(max(E .* G - F.^2, 0));
end

function [X, nr, Xtheta, Xphi] = surfaceFromMaterial(theta, phi, t, l, m)
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
e_phi = [-sin(theta) .* sin(phi), sin(theta) .* cos(phi), zeros(size(theta))];

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
e_theta_phi = [-cos(theta) .* sin(phi), cos(theta) .* cos(phi), zeros(size(theta))];
e_phi_theta = [-cos(theta) .* sin(phi), cos(theta) .* cos(phi), zeros(size(theta))];
e_phi_phi = [-sin(theta) .* cos(phi), -sin(theta) .* sin(phi), zeros(size(theta))];

w = f .* e_theta + g .* e_phi;
w_theta = f_theta .* e_theta + f .* e_theta_theta + g_theta .* e_phi + g .* e_phi_theta;
w_phi = f_phi .* e_theta + f .* e_theta_phi + g_phi .* e_phi + g .* e_phi_phi;

X_rad = [r .* u(:, 1), r .* u(:, 2), r .* u(:, 3)];
X = X_rad + amp_t * w;

Xtheta = dr_dtheta .* u + r .* e_theta + amp_t .* w_theta;
Xphi = dr_dphi .* u + r .* e_phi + amp_t .* w_phi;
N_raw = cross(Xtheta, Xphi, 2);
nr = N_raw ./ vecnorm(N_raw, 2, 2);
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
[~, sig] = tryLargestAlgebraic(L, tol);
if sig > tol
    error('kp:manifold:BadSpectrum', ...
        ['+ve real eigenvalue detected, smooth normals with more nearest ' ...
         'neighbors or check points.']);
end
end

function [V, D] = tryLargestAlgebraic(A, tol)
warnState = warning();
cleanup = onCleanup(@() warning(warnState));
warning('off', 'all');
[V, D] = eigs(A, 1, 'largestreal', 'Tolerance', tol);
end
