function inside = pointsInConvexHull2D(P, X)
%POINTSINCONVEXHULL2D Inside test against the convex hull of a 2D cloud.

if isempty(P)
    inside = false(0, 1);
    return;
end

if size(X, 1) < 3
    inside = true(size(P, 1), 1);
    return;
end

k = convhull(X(:, 1), X(:, 2));
inside = inpolygon(P(:, 1), P(:, 2), X(k, 1), X(k, 2));

