function [UV, origin, basis] = projectToBestFitPlane(X)
%PROJECTTOBESTFITPLANE Project 3D points to a best-fit plane.

origin = mean(X, 1);
Xc = X - origin;

if size(X, 1) < 3
    basis = eye(3, 2);
else
    [~, ~, V] = svd(Xc, 'econ');
    basis = V(:, 1:2);
end

UV = Xc * basis;

