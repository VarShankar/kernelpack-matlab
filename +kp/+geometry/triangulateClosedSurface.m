function triangles = triangulateClosedSurface(points, shrinkFactor)
%TRIANGULATECLOSEDSURFACE Reconstruct a closed surface for visualization.
%   TRIANGLES = TRIANGULATECLOSEDSURFACE(POINTS) uses MATLAB's boundary
%   reconstruction with a moderate shrink factor. This helper is intended
%   for plotting sampled closed surfaces, not for numerical PDE assembly.

arguments
    points (:, 3) double {mustBeFinite}
    shrinkFactor (1, 1) double {mustBeGreaterThanOrEqual(shrinkFactor, 0), ...
        mustBeLessThanOrEqual(shrinkFactor, 1)} = 0.8
end

if size(points, 1) < 4
    error('kp:geometry:triangulateClosedSurface:TooFewPoints', ...
        'At least four points are required to reconstruct a closed surface.');
end

triangles = boundary(points(:, 1), points(:, 2), points(:, 3), shrinkFactor);
end
