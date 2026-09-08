function D = distanceMatrix(X, Y)
%DISTANCEMATRIX Euclidean distance matrix.

xx = sum(X.^2, 2);
yy = sum(Y.^2, 2).';
D2 = bsxfun(@plus, xx, yy) - 2 * (X * Y.');
D2 = max(D2, 0);
D = sqrt(D2);

