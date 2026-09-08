function geometry_examples()
%GEOMETRY_EXAMPLES Minimal geometry examples with visualization.

fprintf('Running 2D smooth closed example...\n');
t = linspace(0, 2*pi, 40).';
t(end) = [];
x = [cos(t), 0.7 * sin(t)];

surf2 = kp.geometry.EmbeddedSurface();
surf2.setDataSites(x);
surf2.buildClosedGeometricModelPS(2, 0.05, size(x, 1));
surf2.buildLevelSetFromGeometricModel([]);

fprintf('2D sample sites: %d\n', surf2.getN());
plotSmoothClosed2D(surf2, x);

fprintf('Running piecewise 2D example...\n');
seg1 = [linspace(0, 1, 20).', zeros(20, 1)];
seg2 = [ones(20, 1), linspace(0, 1, 20).'];
seg3 = [linspace(1, 0, 20).', ones(20, 1)];
seg4 = [zeros(20, 1), linspace(1, 0, 20).'];

piece2 = kp.geometry.PiecewiseSmoothEmbeddedSurface();
piece2.generatePiecewiseSmoothSurfaceBySegment( ...
    {seg1, seg2, seg3, seg4}, [false false false false], 0.05, 1, 2, 2);
piece2.buildLevelSet();

fprintf('Piecewise 2D boundary nodes: %d\n', size(piece2.getBdryNodes(), 1));
plotPiecewise2D(piece2, {seg1, seg2, seg3, seg4});

fprintf('Running 3D smooth closed example...\n');
X = kp.geometry.fibonacciSphere(120);
uv = kp.geometry.cart2sphRows(X);
r = 1 + 0.15 * cos(3 * uv(:, 1)) .* cos(2 * uv(:, 2));
pts = X .* r;

surf3 = kp.geometry.EmbeddedSurface();
surf3.setDataSites(pts);
surf3.buildClosedGeometricModelPS(3, 0.2, size(pts, 1));
surf3.buildLevelSetFromGeometricModel([]);

fprintf('3D sample sites: %d\n', surf3.getN());
plotSmoothClosed3D(surf3, pts, 'Smooth Closed Surface (3D)');

fprintf('Running piecewise 3D seam example...\n');
faces = makeCubeFaces(13);

piece3 = kp.geometry.PiecewiseSmoothEmbeddedSurface();
piece3.generatePiecewiseSmoothSurfaceBySegment( ...
    faces, false(1, numel(faces)), 0.18, 1, 2, 2, true, 10);
piece3.buildLevelSet();

fprintf('Piecewise 3D boundary nodes: %d\n', size(piece3.getBdryNodes(), 1));
fprintf('Piecewise 3D corner flags: %d\n', nnz(piece3.getCornerFlags()));
plotPiecewise3D(piece3, faces, 'Piecewise Closed Surface (3D)');
end

function plotSmoothClosed2D(surface, dataSites)
xb = surface.getSampleSites();
nrmls = surface.getNrmls();

figure('Name', 'Smooth Closed Curve 2D', 'Color', 'w');
tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
plot(dataSites(:, 1), dataSites(:, 2), 'ko', 'MarkerFaceColor', [0.2 0.2 0.2]);
axis equal;
grid on;
title('Data sites');
xlabel('x');
ylabel('y');

nexttile;
plot(xb(:, 1), xb(:, 2), 'b.', 'MarkerSize', 10);
axis equal;
grid on;
title('Boundary samples');
xlabel('x');
ylabel('y');

nexttile;
plot(xb(:, 1), xb(:, 2), 'b.', 'MarkerSize', 10);
hold on;
drawNormalSegments2D(xb, nrmls, 1, 0.035, [0.85 0.2 0.2]);
axis equal;
grid on;
title('Boundary-sample normals');
xlabel('x');
ylabel('y');
end

function plotPiecewise2D(surface, segments)
xb = surface.getBdryNodes();
nrmls = surface.getBdryNrmls();
cornerMask = logical(surface.getCornerFlags());

figure('Name', 'Piecewise Boundary 2D', 'Color', 'w');
tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
hold on;
for k = 1:numel(segments)
    plot(segments{k}(:, 1), segments{k}(:, 2), 'k.-');
end
axis equal;
grid on;
title('Input segments');
xlabel('x');
ylabel('y');

nexttile;
plot(xb(:, 1), xb(:, 2), 'b.', 'MarkerSize', 10);
if any(cornerMask)
    hold on;
    plot(xb(cornerMask, 1), xb(cornerMask, 2), 'mo', 'MarkerSize', 7, 'LineWidth', 1.2);
end
axis equal;
grid on;
title('Assembled boundary');
xlabel('x');
ylabel('y');

nexttile;
plot(xb(:, 1), xb(:, 2), 'b.', 'MarkerSize', 10);
hold on;
drawNormalSegments2D(xb, nrmls, 1, 0.035, [0.85 0.2 0.2]);
if any(cornerMask)
    plot(xb(cornerMask, 1), xb(cornerMask, 2), 'mo', 'MarkerSize', 7, 'LineWidth', 1.2);
end
axis equal;
grid on;
title('Boundary normals and corners');
xlabel('x');
ylabel('y');
end

function plotSmoothClosed3D(surface, dataSites, figName)
xb = surface.getSampleSites();

figure('Name', figName, 'Color', 'w');
tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
plot3(dataSites(:, 1), dataSites(:, 2), dataSites(:, 3), 'k.', 'MarkerSize', 10);
axis equal;
grid on;
view(3);
title('Data cloud');
xlabel('x');
ylabel('y');
zlabel('z');

nexttile;
plot3(xb(:, 1), xb(:, 2), xb(:, 3), 'b.', 'MarkerSize', 10);
axis equal;
grid on;
view(3);
title('Boundary cloud');
xlabel('x');
ylabel('y');
zlabel('z');

nexttile;
tri = kp.geometry.triangulateClosedSurface(xb);
trisurf(tri, xb(:, 1), xb(:, 2), xb(:, 3), ...
    'FaceColor', [0.2 0.6 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
hold on;
plot3(xb(:, 1), xb(:, 2), xb(:, 3), 'k.', 'MarkerSize', 9);
axis equal;
grid on;
view(3);
camlight headlight;
lighting gouraud;
title('Triangulated surface + boundary cloud');
xlabel('x');
ylabel('y');
zlabel('z');
end

function plotPiecewise3D(surface, faces, figName)
xb = surface.getBdryNodes();
cornerMask = logical(surface.getCornerFlags());

raw = vertcat(faces{:});

figure('Name', figName, 'Color', 'w');
tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
plot3(raw(:, 1), raw(:, 2), raw(:, 3), 'k.', 'MarkerSize', 9);
axis equal;
grid on;
view(3);
title('Patch data cloud');
xlabel('x');
ylabel('y');
zlabel('z');

nexttile;
plot3(xb(:, 1), xb(:, 2), xb(:, 3), 'b.', 'MarkerSize', 9);
hold on;
if any(cornerMask)
    plot3(xb(cornerMask, 1), xb(cornerMask, 2), xb(cornerMask, 3), ...
        'mo', 'MarkerSize', 6, 'LineWidth', 1.1);
end
axis equal;
grid on;
view(3);
title('Boundary cloud');
xlabel('x');
ylabel('y');
zlabel('z');

nexttile;
tri = kp.geometry.triangulateClosedSurface(xb);
trisurf(tri, xb(:, 1), xb(:, 2), xb(:, 3), ...
    'FaceColor', [0.85 0.45 0.25], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
hold on;
plot3(xb(:, 1), xb(:, 2), xb(:, 3), 'k.', 'MarkerSize', 8);
axis equal;
grid on;
view(3);
camlight headlight;
lighting gouraud;
title('Triangulated surface + boundary cloud');
xlabel('x');
ylabel('y');
zlabel('z');
end

function faces = makeCubeFaces(n)
u = linspace(-1, 1, n);
[a, b] = ndgrid(u, u);

faces = {
    [ones(numel(a), 1), a(:), b(:)]
    [-ones(numel(a), 1), a(:), b(:)]
    [a(:), ones(numel(a), 1), b(:)]
    [a(:), -ones(numel(a), 1), b(:)]
    [a(:), b(:), ones(numel(a), 1)]
    [a(:), b(:), -ones(numel(a), 1)]
    };
end

function drawNormalSegments2D(points, normals, stride, scale, color)
idx = 1:stride:size(points, 1);
tips = points(idx, :) + scale * normals(idx, :);

for k = 1:numel(idx)
    line([points(idx(k), 1), tips(k, 1)], [points(idx(k), 2), tips(k, 2)], ...
        'Color', color, 'LineWidth', 1);
end
end
