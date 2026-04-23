function geometry_checks()
%GEOMETRY_CHECKS Lightweight regression checks for current geometry layer.

t = linspace(0, 2*pi, 60).';
t(end) = [];
x = [cos(t), 0.7 * sin(t)];
s2 = kp.geometry.EmbeddedSurface();
s2.setDataSites(x);
s2.buildClosedGeometricModelPS(2, 0.05, size(x, 1), 120);
s2.buildLevelSetFromGeometricModel([]);
xb2 = s2.getSampleSites();
nr2 = s2.getNrmls();
ctr2 = mean(xb2, 1);
assert(mean(sum((xb2 - ctr2) .* nr2, 2)) > 0, '2D closed normals should point outward.');

X = kp.geometry.fibonacciSphere(160);
uv = kp.geometry.cart2sphRows(X);
r = 1 + 0.15 * cos(3 * uv(:, 1)) .* cos(2 * uv(:, 2));
pts = X .* r;
s3 = kp.geometry.EmbeddedSurface();
s3.setDataSites(pts);
s3.buildClosedGeometricModelPS(3, 0.18, size(pts, 1), 180);
s3.buildLevelSetFromGeometricModel([]);
assert(size(s3.getSampleSites(), 2) == 3, '3D sample sites should be 3D.');
assert(s3.getN() > 20, '3D WSE thinning should keep a reasonable sample set.');

p1 = [rand(120, 2) * 2 - 1, zeros(120, 1)];
p2 = [rand(120, 2) * 2 - 1, ones(120, 1)];
p3 = [ones(120, 1), rand(120, 1) * 2 - 1, rand(120, 1)];
p4 = [ones(120, 1), rand(120, 1) * 2 - 1, rand(120, 1)];

piece3 = kp.geometry.PiecewiseSmoothEmbeddedSurface();
piece3.generatePiecewiseSmoothSurfaceBySegment({p1, p2, p3, p4}, ...
    [false false false false], 0.25, 1, 2, 2, true, 8);
piece3.buildLevelSet();

assert(size(piece3.getBdryNodes(), 2) == 3, 'Piecewise 3D nodes should be 3D.');
assert(nnz(piece3.getCornerFlags()) > 0, 'Repeated seam patches should produce seam or corner flags.');
assert(size(piece3.getUniformBdryNodes(), 1) >= size(piece3.getBdryNodes(), 1), ...
    'Uniform boundary cloud should be at least as large as the assembled cloud.');

disp('geometry checks passed');
