function y = applyMatrixPower(A, x, k)
%APPLYMATRIXPOWER Apply A^k to x without explicitly forming A^k.

y = x;
for j = 1:k
    y = A * y;
end
end
