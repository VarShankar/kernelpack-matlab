function X = solveLocalAugmentedSystem(A, B)
%SOLVELOCALAUGMENTEDSYSTEM Solve local RBF/poly systems robustly.
%   The usual local systems are small enough that a condition estimate is
%   cheap. Well-conditioned systems use MATLAB's decomposition path. Poorly
%   conditioned systems fall back to column-pivoted QR, which is more robust
%   for nearly dependent polynomial/RBF columns on difficult surface stencils.

[X, ~] = kp.manifold.detail.solveAndFactorLocalAugmentedSystem(A, B);
end
