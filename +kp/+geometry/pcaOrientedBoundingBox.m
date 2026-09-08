function box = pcaOrientedBoundingBox(X)
%PCAORIENTEDBOUNDINGBOX PCA-based oriented bounding box.

dim = size(X, 2);
if isempty(X)
    box = struct('p', zeros(0, dim), 'V', eye(dim), 'D', zeros(dim, 1));
    return;
end

if size(X, 1) == 1
    V = eye(dim);
    D = zeros(dim, 1);
    Xo = X;
else
    [V, S] = eig(cov(X));
    D = diag(S);
    Xo = (V.' * X.').';
end

xmin = min(Xo, [], 1);
xmax = max(Xo, [], 1);

if dim == 2
    p = [xmin(1) xmin(2); ...
         xmax(1) xmin(2); ...
         xmax(1) xmax(2); ...
         xmin(1) xmax(2)];
else
    p = [xmin(1) xmin(2) xmin(3); ...
         xmax(1) xmin(2) xmin(3); ...
         xmax(1) xmax(2) xmin(3); ...
         xmin(1) xmax(2) xmin(3); ...
         xmin(1) xmin(2) xmax(3); ...
         xmax(1) xmin(2) xmax(3); ...
         xmax(1) xmax(2) xmax(3); ...
         xmin(1) xmax(2) xmax(3)];
end

box = struct('p', (V * p.').', 'V', V, 'D', D);

