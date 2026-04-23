function nodes_checks()
%NODES_CHECKS Lightweight checks for box Poisson disk sampling.

[x2a, info2a] = kp.nodes.generatePoissonNodesInBox(0.075, [0 0], [1 1], ...
    'Seed', 19, 'StripCount', 5);
[x2b, info2b] = kp.nodes.generatePoissonNodesInBox(0.075, [0 0], [1 1], ...
    'Seed', 19, 'StripCount', 5);
[x2c, ~] = kp.nodes.generatePoissonNodesInBox(0.075, [0 0], [1 1], ...
    'Seed', 23, 'StripCount', 5);

assert(isequal(x2a, x2b), 'Seeded 2D box sampling should be deterministic.');
assert(~isequal(x2a, x2c), 'Different seeds should usually produce different 2D clouds.');
assert(all(x2a >= 0, 'all') && all(x2a <= 1, 'all'), '2D nodes must stay in the box.');
assert(minPairDistance(x2a) >= 0.075 * (1 - 1e-10), '2D nodes must respect the exclusion radius.');
assert(info2a.deterministic && info2b.deterministic, 'Seeded runs should report deterministic mode.');

[x3, info3] = kp.nodes.generatePoissonNodesInBox(0.16, [0 0 0], [1 1 1], ...
    'Seed', 31, 'StripCount', 5);
assert(size(x3, 2) == 3, '3D nodes should have three coordinates.');
assert(minPairDistance(x3) >= 0.16 * (1 - 1e-10), '3D nodes must respect the exclusion radius.');
assert(info3.strip_count == 5, 'Explicit strip count should be preserved.');

[x4, ~] = kp.nodes.generatePoissonNodesInBox(0.35, zeros(1, 4), ones(1, 4), ...
    'Seed', 5, 'StripCount', 3, 'UseParallel', false);
assert(size(x4, 2) == 4, 'Sampler should work in dimensions beyond 3.');
assert(minPairDistance(x4) >= 0.35 * (1 - 1e-10), '4D nodes must respect the exclusion radius.');

gen = kp.nodes.DomainNodeGenerator();
gen.generatePoissonNodes(0.1, [0 0], [1 1], 'Seed', 13, 'StripCount', 4);
assert(~isempty(gen.getRawPoissonInteriorNodes()), 'DomainNodeGenerator should store the raw Poisson cloud.');

t = linspace(0, 2*pi, 60).';
t(end) = [];
curve = [cos(t), 0.7 * sin(t)];
surface = kp.geometry.EmbeddedSurface();
surface.setDataSites(curve);
surface.buildClosedGeometricModelPS(2, 0.05, size(curve, 1), 120);
surface.buildLevelSetFromGeometricModel([]);

gen2 = kp.nodes.DomainNodeGenerator();
gen2.generateInteriorNodesFromGeometry(surface, 0.08, 'Seed', 29, 'StripCount', 5);
raw2 = gen2.getRawPoissonInteriorNodes();
int2 = gen2.getInteriorNodes();
assert(size(int2, 1) < size(raw2, 1), 'Geometry clipping should remove points outside the surface bounding box fill.');
phi2 = surface.getLevelSet().Evaluate(int2);
assert(all(phi2 <= 1e-10), 'Interior nodes should lie on the negative side of the level set.');

seg1 = [linspace(0,1,25).', zeros(25,1)];
seg2 = [ones(25,1), linspace(0,1,25).'];
seg3 = [linspace(1,0,25).', ones(25,1)];
seg4 = [zeros(25,1), linspace(1,0,25).'];
piece = kp.geometry.PiecewiseSmoothEmbeddedSurface();
piece.generatePiecewiseSmoothSurfaceBySegment({seg1, seg2, seg3, seg4}, ...
    [false false false false], 0.05, 1, 2, 2);
piece.buildLevelSet();
[clippedPiece, ~, phiPiece] = kp.nodes.clipPointsByGeometry(raw2, piece, 'Keep', 'outside');
assert(~isempty(clippedPiece), 'Outside clipping against a piecewise geometry should keep some points.');
assert(all(phiPiece(phiPiece >= 0) >= 0), 'Outside clipping should use the positive level-set side.');

disp('nodes checks passed');
end

function dmin = minPairDistance(X)
    if size(X, 1) < 2
        dmin = inf;
        return;
    end
    D = kp.geometry.distanceMatrix(X, X);
    D(1:size(D, 1)+1:end) = inf;
    dmin = min(D, [], 'all');
end
