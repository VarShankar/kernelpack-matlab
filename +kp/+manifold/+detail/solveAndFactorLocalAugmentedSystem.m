function [X, factor, relres] = solveAndFactorLocalAugmentedSystem(A, B, tol)
%SOLVEANDFACTORLOCALAUGMENTEDSYSTEM Direct local solve with LU->QR fallback.
%   LU is the intended fast path. Column-pivoted QR is used only if the LU
%   triangular factor is numerically singular or if the LU solve leaves an
%   unacceptable residual.

if nargin < 3 || isempty(tol)
    tol = 1.0e-10;
end

factor = kp.manifold.detail.factorLocalAugmentedSystem(A);
if ~isfinite(factor.reciprocalCondition) || factor.reciprocalCondition < 1.0e-12
    factor = qrFactor(A);
    X = kp.manifold.detail.solveLocalFactor(factor, B);
    relres = localRelativeResidual(A, X, B);
    return;
end
X = kp.manifold.detail.solveLocalFactor(factor, B);
relres = localRelativeResidual(A, X, B);
if relres <= tol
    return;
end

factor = qrFactor(A);
X = kp.manifold.detail.solveLocalFactor(factor, B);
relres = localRelativeResidual(A, X, B);
end

function relres = localRelativeResidual(A, X, B)
relres = norm(B - A * X, 'fro') / max(norm(B, 'fro'), eps);
end

function factor = qrFactor(A)
[Q, R, E] = qr(A, 0);
diagR = abs(diag(R));
if isempty(diagR)
    rankA = 0;
else
    tol = max(size(A)) * eps(max(diagR));
    rankA = sum(diagR > tol);
end

factor = struct( ...
    'kind', "qr", ...
    'L', [], ...
    'U', [], ...
    'P', [], ...
    'Q', Q, ...
    'R', R, ...
    'E', E, ...
    'rank', rankA, ...
    'reciprocalCondition', rcond(R));
end
