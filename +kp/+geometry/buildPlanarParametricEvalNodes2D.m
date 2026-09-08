function Xe = buildPlanarParametricEvalNodes2D(interpNodes, N)
%BUILDPLANARPARAMETRICEVALNODES2D Poisson evaluation nodes inside a projected hull.

[X, ~, ~] = kp.geometry.projectToBestFitPlane(interpNodes);
x0 = min(X(:, 1));
xm = max(X(:, 1));
y0 = min(X(:, 2));
ym = max(X(:, 2));
w = max(abs(xm - x0), 1e-12);
h = max(abs(ym - y0), 1e-12);
radius = sqrt((w * h) / max(2 * N, 1));

Xe = zeros(0, 2);
for attempt = 0:4
    samples = kp.nodes.generatePoissonNodesInBox(radius, [x0, y0], [xm, ym], ...
        'UseParallel', false, 'Seed', attempt);
    if ~isempty(samples)
        Xe = samples(kp.geometry.pointsInConvexHull2D(samples, X), :);
    else
        Xe = zeros(0, 2);
    end

    if size(Xe, 1) >= max(1, floor(N / 2)) || attempt == 4
        break;
    end
    radius = 0.85 * radius;
end
