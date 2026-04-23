function p = poly_eval(a, b, x, N, d)
%POLY_EVAL Thin wrapper for shared orthonormal polynomial evaluation.
    if nargin < 5
        d = 0;
    end
    p = kp.poly.JacobiPolynomials.evaluate(a, b, x, N, d);
end
