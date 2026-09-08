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
assert(info3.strip_count == 1, 'Deterministic sampling should collapse to one canonical strip, matching C++.');

[x4, ~] = kp.nodes.generatePoissonNodesInBox(0.35, zeros(1, 4), ones(1, 4), ...
    'Seed', 5, 'StripCount', 3, 'UseParallel', false);
assert(size(x4, 2) == 4, 'Sampler should work in dimensions beyond 3.');
assert(minPairDistance(x4) >= 0.35 * (1 - 1e-10), '4D nodes must respect the exclusion radius.');

radFunc = @(p, hmin) hmin * (1 + 0.5 * (p(1) > 0.5));
[xVar, infoVar] = kp.nodes.generatePoissonNodesInBox(radFunc, [0 0], [1 1], ...
    'MinRadius', 0.06, 'Seed', 41, 'StripCount', 5, 'UseParallel', false);
assert(strcmp(infoVar.mode, 'variable_radius'), 'Variable-density mode should be reported correctly.');
assert(size(xVar, 2) == 2, 'Variable-density nodes should keep the box dimension.');
assert(minPairDistance(xVar) >= 0.06 * (1 - 1e-10), 'Variable-density nodes must respect the minimum radius.');

gen = kp.nodes.DomainNodeGenerator();
gen.generatePoissonNodes(0.1, [0 0], [1 1], 'Seed', 13, 'StripCount', 4);
assert(~isempty(gen.getRawPoissonInteriorNodes()), 'DomainNodeGenerator should store the raw Poisson cloud.');

t = linspace(0, 2*pi, 60).';
t(end) = [];
curve = [cos(t), 0.7 * sin(t)];
surface = kp.geometry.EmbeddedSurface();
surface.setDataSites(curve);
surface.buildClosedGeometricModelPS(2, 0.05, size(curve, 1));
surface.buildLevelSetFromGeometricModel([]);

gen2 = kp.nodes.DomainNodeGenerator();
gen2.generateInteriorNodesFromGeometry(surface, 0.08, 'Seed', 29, 'StripCount', 5);
raw2 = gen2.getRawPoissonInteriorNodes();
int2 = gen2.getInteriorNodes();
assert(size(int2, 1) < size(raw2, 1), 'Geometry clipping should remove points outside the surface bounding box fill.');
phi2 = surface.getLevelSet().Evaluate(int2);
assert(all(phi2 >= 0.08 - 1e-10), 'Interior nodes should stay at least h away from the boundary.');

gen3 = kp.nodes.DomainNodeGenerator();
gen3.generateInteriorNodesFromGeometry(surface, 0.08, ...
    'Seed', 29, 'StripCount', 5, ...
    'DoOuterRefinement', true, ...
    'OuterFractionOfh', 0.5, ...
    'OuterRefinementZoneSizeAsMultipleOfh', 2.0);
int3 = gen3.getInteriorNodes();
phi3 = surface.getLevelSet().Evaluate(int3);
assert(any(phi3 < 0.16 - 1e-10), 'Outer refinement should populate the near-boundary band.');
assert(all(phi3 >= 0.04 - 1e-10), 'Clipping should use the smallest active radius in the refined node set.');
assert(size(int3, 1) >= size(int2, 1), 'Outer refinement should not reduce the interior node count.');

seg1 = [linspace(0,1,25).', zeros(25,1)];
seg2 = [ones(25,1), linspace(0,1,25).'];
seg3 = [linspace(1,0,25).', ones(25,1)];
seg4 = [zeros(25,1), linspace(1,0,25).'];
piece = kp.geometry.PiecewiseSmoothEmbeddedSurface();
piece.generatePiecewiseSmoothSurfaceBySegment({seg1, seg2, seg3, seg4}, ...
    [false false false false], 0.05, 1, 2, 2);
piece.buildLevelSet();
[clippedPiece, keepPiece, phiPiece] = kp.nodes.clipPointsByGeometry(raw2, piece, 'Keep', 'outside');
assert(~isempty(clippedPiece), 'Outside clipping against a piecewise geometry should keep some points.');
assert(all(phiPiece(keepPiece) <= 1e-10), 'Outside clipping should use the negative level-set side.');

descriptor = gen3.buildDomainDescriptorFromGeometry(surface, 0.08, ...
    'Seed', 29, 'StripCount', 5, ...
    'DoOuterRefinement', true, ...
    'OuterFractionOfh', 0.5, ...
    'OuterRefinementZoneSizeAsMultipleOfh', 2.0);
assert(~isempty(descriptor.getInteriorNodes()), 'DomainDescriptor should store interior nodes.');
assert(size(descriptor.getBdryNodes(), 1) == size(descriptor.getGhostNodes(), 1), ...
    'Each boundary node should receive one ghost node.');
assert(size(descriptor.getNrmls(), 1) == size(descriptor.getBdryNodes(), 1), ...
    'Boundary normals should match the boundary node count.');
assert(size(descriptor.getAllNodes(), 1) == ...
    size(descriptor.getInteriorNodes(), 1) + size(descriptor.getBdryNodes(), 1) + size(descriptor.getGhostNodes(), 1), ...
    'All-node cloud should concatenate interior, boundary, and ghost nodes.');

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
