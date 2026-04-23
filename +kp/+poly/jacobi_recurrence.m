function [a, b] = jacobi_recurrence(N, alph, bet)
%JACOBI_RECURRENCE Thin wrapper for shared Jacobi recurrence coefficients.
    [a, b] = kp.poly.JacobiPolynomials.recurrence(N, alph, bet);
end
