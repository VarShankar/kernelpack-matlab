function [a, b] = chebyshev_recurrence(N)
%CHEBYSHEV_RECURRENCE Orthonormal Chebyshev recurrence for the arcsine measure.

a = zeros(N, 1);
b = ones(N, 1);

if N > 1
    b(2) = 1 / 2;
    b(3:end) = 1 / 4;
end
end
