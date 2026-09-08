function results = moving_surface_adr_tp_operator_diagnostics(xiVals, Nvals)
%MOVING_SURFACE_ADR_TP_OPERATOR_DIAGNOSTICS Check tangent-plane operator rates.
%   This diagnostic isolates the RBF-FD tangent-plane differentiation
%   matrices from time stepping, rearrangement, transfer, and forcing.

if nargin < 1 || isempty(xiVals)
    xiVals = [2, 4, 6];
end
if nargin < 2 || isempty(Nvals)
    Nvals = [256, 576, 1024, 1600];
end

xiVals = xiVals(:).';
Nvals = Nvals(:).';
results = table();

for xi = xiVals
    op = kp.manifold.rbffdop(2, xi, 2, 0);
    for N = Nvals
        [X, nr] = sphereNodes(N);
        tree = KDTreeSearcher(X);
        [L, Gx, Gy, Gz] = kp.manifold.FormSurfaceDiffOpsTP( ...
            X, op.rbf, op.drbfor, op.d2rbf, nr, tree, op.stencilSize, op.ell);

        u = sphereTestField(X);
        gradExact = sphereTestGradient(X);
        lapExact = sphereTestLaplacian(X);
        gradNum = [Gx * u, Gy * u, Gz * u];
        lapNum = L * u;

        row = table();
        row.xi = xi;
        row.N = N;
        row.sqrtN = sqrt(N);
        row.h = sqrt(4 * pi / N);
        row.gradRelerr = norm(gradNum - gradExact, 'fro') / norm(gradExact, 'fro');
        row.lapRelerr = norm(lapNum - lapExact) / norm(lapExact);
        results = [results; row]; %#ok<AGROW>

        fprintf('xi=%d N=%d grad=%.6e lap=%.6e\n', ...
            xi, N, row.gradRelerr, row.lapRelerr);
    end
end

results.gradRate = nan(height(results), 1);
results.lapRate = nan(height(results), 1);
for xi = xiVals
    ids = find(results.xi == xi);
    for k = 2:numel(ids)
        i0 = ids(k - 1);
        i1 = ids(k);
        results.gradRate(i1) = log(results.gradRelerr(i0) / results.gradRelerr(i1)) / ...
            log(results.h(i0) / results.h(i1));
        results.lapRate(i1) = log(results.lapRelerr(i0) / results.lapRelerr(i1)) / ...
            log(results.h(i0) / results.h(i1));
    end
end

writetable(results, fullfile(pwd, 'moving_surface_adr_tp_operator_diagnostics.csv'));
end

function [X, nr] = sphereNodes(N)
i = (0:N-1).';
golden = pi * (3 - sqrt(5));
z = 1 - 2 * (i + 0.5) / N;
r = sqrt(max(1 - z .^ 2, 0));
phi = mod(i * golden, 2 * pi);
X = [r .* cos(phi), r .* sin(phi), z];
nr = X;
end

function u = sphereTestField(X)
x = X(:, 1);
y = X(:, 2);
z = X(:, 3);
u = x + 0.5 * y + 0.25 * (3 * z .^ 2 - 1) + 0.1 * x .* y;
end

function grad = sphereTestGradient(X)
x = X(:, 1);
y = X(:, 2);
z = X(:, 3);
gradAmbient = [1 + 0.1 * y, 0.5 + 0.1 * x, 1.5 * z];
dotNormal = sum(gradAmbient .* X, 2);
grad = gradAmbient - dotNormal .* X;
end

function lap = sphereTestLaplacian(X)
x = X(:, 1);
y = X(:, 2);
z = X(:, 3);
% x and y are degree-one spherical harmonics, 3z^2-1 is degree two, and xy
% is degree two, so Delta_Gamma Y_l = -l(l+1)Y_l on the unit sphere.
lap = -2 * (x + 0.5 * y) - 6 * (0.25 * (3 * z .^ 2 - 1) + 0.1 * x .* y);
end
