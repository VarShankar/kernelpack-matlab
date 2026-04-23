function [a, b] = legendre_recurrence(N)
%LEGENDRE_RECURRENCE Orthonormal Legendre recurrence coefficients.
    [a, b] = kp.poly.jacobi_recurrence(N, 0, 0);
end
