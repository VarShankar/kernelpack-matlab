function report = moving_surface_adr_tp_local_kernel_diagnostic(N, xi, dt, varargin)
%MOVING_SURFACE_ADR_TP_LOCAL_KERNEL_DIAGNOSTIC Time local direct vs defect solves.
%   This isolates the local weight-solve kernels. It prebuilds the same
%   moved-surface stencil matrices and right-hand sides for both methods,
%   then times only direct factor-and-solve against defect correction using
%   cached old inverse actions.

if nargin < 1 || isempty(N)
    N = 1024;
end
if nargin < 2 || isempty(xi)
    xi = 4;
end
if nargin < 3 || isempty(dt)
    dt = 4.13793103448276e-3;
end

parser = inputParser();
parser.addParameter('DefectTolerance', 1.0e-4);
parser.addParameter('MaxDefectIterations', 2);
parser.addParameter('Repeats', 3);
parser.parse(varargin{:});

op = kp.manifold.rbffdop(2, xi, 2, 0);
n = op.stencilSize;
m = n + op.polyM;
recurrence = @(K) kp.poly.jacobi_recurrence(K, 0, 0);
polyIndices = kp.poly.total_degree_indices(2, op.ell);

U = kp.geometry.fibonacciSphere(N);
U = U ./ vecnorm(U, 2, 2);
radius0 = 1.4;
radius1 = sqrt(radius0^2 - 4 * dt);
X0 = radius0 .* U;
X1 = radius1 .* U;
normals = U;

tree = KDTreeSearcher(X1);
IDX = knnsearch(tree, X1, 'k', n);

Anew = zeros(m, m, N);
Bnew = zeros(m, 3, N);
oldInv = zeros(m, m, N);

fprintf('Prebuilding %d local systems, xi=%d, stencil=%d, matrix=%d\n', N, xi, n, m);
for i = 1:N
    idx = IDX(i, :);
    normal = normals(idx(1), :).';
    [A0, ~] = kp.manifold.detail.tangentPlaneLocalSystem( ...
        X0(idx, :), normal, op.rbf, op.drbfor, op.d2rbf, polyIndices, recurrence);
    [A1, B1] = kp.manifold.detail.tangentPlaneLocalSystem( ...
        X1(idx, :), normal, op.rbf, op.drbfor, op.d2rbf, polyIndices, recurrence);
    oldInv(:, :, i) = A0 \ eye(m);
    Anew(:, :, i) = A1;
    Bnew(:, :, i) = B1;
end

directTimes = zeros(parser.Results.Repeats, 1);
defectTimes = zeros(parser.Results.Repeats, 1);
maxRelDiff = 0;
maxResidual = 0;
defectIters = zeros(N, 1);
defectOk = false(N, 1);

for r = 1:parser.Results.Repeats
    tic;
    Wdirect = zeros(m, 3, N);
    for i = 1:N
        Wdirect(:, :, i) = Anew(:, :, i) \ Bnew(:, :, i);
    end
    directTimes(r) = toc;

    tic;
    Wdefect = zeros(m, 3, N);
    for i = 1:N
        [Wdefect(:, :, i), defectOk(i), defectIters(i), rr] = defectCorrect( ...
            Anew(:, :, i), Bnew(:, :, i), oldInv(:, :, i), ...
            parser.Results.DefectTolerance, parser.Results.MaxDefectIterations);
        maxResidual = max(maxResidual, rr);
    end
    defectTimes(r) = toc;

    maxRelDiff = max(maxRelDiff, norm(Wdefect(:) - Wdirect(:)) / max(norm(Wdirect(:)), eps));
end

report = struct();
report.N = N;
report.xi = xi;
report.dt = dt;
report.stencilSize = n;
report.matrixSize = m;
report.defectTolerance = parser.Results.DefectTolerance;
report.maxDefectIterations = parser.Results.MaxDefectIterations;
report.directTimes = directTimes;
report.defectTimes = defectTimes;
report.bestDirectTime = min(directTimes);
report.bestDefectTime = min(defectTimes);
report.speedup = report.bestDirectTime / report.bestDefectTime;
report.defectAccepted = nnz(defectOk);
report.defectRejected = N - report.defectAccepted;
report.meanDefectIterations = mean(defectIters(defectOk), 'omitnan');
report.maxResidual = maxResidual;
report.maxRelativeDifference = maxRelDiff;

fprintf('Local kernel diagnostic\n');
fprintf('  best direct factor+solve: %.6g s\n', report.bestDirectTime);
fprintf('  best defect correction:  %.6g s\n', report.bestDefectTime);
fprintf('  speedup direct/defect:   %.3f\n', report.speedup);
fprintf('  accepted/rejected:       %d / %d\n', report.defectAccepted, report.defectRejected);
fprintf('  mean defect iterations:  %.3f\n', report.meanDefectIterations);
fprintf('  max residual:            %.3e\n', report.maxResidual);
fprintf('  max relative W diff:     %.3e\n', report.maxRelativeDifference);
end

function [W, ok, iter, rr] = defectCorrect(A, B, invA, defectTol, maxIters)
W = invA * B;
rr = norm(B - A * W, 'fro') / max(norm(B, 'fro'), eps);
ok = rr <= defectTol;
iter = 0;

while ~ok && iter < maxIters
    residual = B - A * W;
    W = W + invA * residual;
    iter = iter + 1;
    rr = norm(B - A * W, 'fro') / max(norm(B, 'fro'), eps);
    ok = rr <= defectTol;
end
end
