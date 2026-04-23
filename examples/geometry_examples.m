function geometry_examples()
%GEOMETRY_EXAMPLES Minimal examples for current geometry objects.

fprintf('Running 2D smooth closed example...\n');
t = linspace(0, 2*pi, 40).';
t(end) = [];
x = [cos(t), 0.7 * sin(t)];

surf2 = kp.geometry.EmbeddedSurface();
surf2.setDataSites(x);
surf2.buildClosedGeometricModelPS(2, 0.05, size(x, 1), 120);
surf2.buildLevelSetFromGeometricModel([]);

fprintf('2D sample sites: %d\n', surf2.getN());

fprintf('Running 3D smooth closed example...\n');
X = kp.geometry.fibonacciSphere(120);
uv = kp.geometry.cart2sphRows(X);
r = 1 + 0.15 * cos(3 * uv(:, 1)) .* cos(2 * uv(:, 2));
pts = X .* r;

surf3 = kp.geometry.EmbeddedSurface();
surf3.setDataSites(pts);
surf3.buildClosedGeometricModelPS(3, 0.2, size(pts, 1), 160);
surf3.buildLevelSetFromGeometricModel([]);

fprintf('3D sample sites: %d\n', surf3.getN());

fprintf('Running piecewise 3D seam example...\n');
p1 = [rand(100, 2) * 2 - 1, zeros(100, 1)];
p2 = [rand(100, 2) * 2 - 1, ones(100, 1)];
p3 = [ones(100, 1), rand(100, 1) * 2 - 1, rand(100, 1)];
p4 = [ones(100, 1), rand(100, 1) * 2 - 1, rand(100, 1)];

piece3 = kp.geometry.PiecewiseSmoothEmbeddedSurface();
piece3.generatePiecewiseSmoothSurfaceBySegment({p1, p2, p3, p4}, ...
    [false false false false], 0.25, 1, 2, 2, true, 8);
piece3.buildLevelSet();

fprintf('Piecewise 3D boundary nodes: %d\n', size(piece3.getBdryNodes(), 1));
fprintf('Piecewise 3D corner flags: %d\n', nnz(piece3.getCornerFlags()));
