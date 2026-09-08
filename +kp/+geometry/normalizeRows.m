function Y = normalizeRows(X)
%NORMALIZEROWS Normalize each row of an array.

nrm = sqrt(sum(X.^2, 2));
nrm = max(nrm, eps);
Y = X ./ nrm;

