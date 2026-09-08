function p = mpoly_eval(x, alpha, recurrenceHandle, d)
%MPOLY_EVAL Thin wrapper for shared tensor-product polynomial evaluation.
    if nargin < 4
        d = zeros(1, size(x, 2));
    end
    p = kp.poly.JacobiPolynomials.tensorEvaluate(x, alpha, recurrenceHandle, d);
end
