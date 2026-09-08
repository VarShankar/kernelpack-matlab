function nodes_examples()
%NODES_EXAMPLES Minimal examples for box Poisson disk sampling.

fprintf('Running 2D seeded box Poisson example...\n');
[x2, info2] = kp.nodes.generatePoissonNodesInBox(0.06, [0 0], [1 1], ...
    'Seed', 7, 'StripCount', 5);
fprintf('2D nodes: %d (seed %d, strips %d)\n', size(x2, 1), info2.seed, info2.strip_count);

fprintf('Running 4D seeded box Poisson example...\n');
[x4, info4] = kp.nodes.generatePoissonNodesInBox(0.22, zeros(1, 4), ones(1, 4), ...
    'Seed', 11, 'StripCount', 5, 'UseParallel', false);
fprintf('4D nodes: %d (seed %d, strips %d)\n', size(x4, 1), info4.seed, info4.strip_count);

fprintf('Running geometry-clipped 2D interior-node example...\n');
t = linspace(0, 2*pi, 50).';
t(end) = [];
curve = [cos(t), 0.7 * sin(t)];
surface = kp.geometry.EmbeddedSurface();
surface.setDataSites(curve);
surface.buildClosedGeometricModelPS(2, 0.05, size(curve, 1));
surface.buildLevelSetFromGeometricModel([]);

generator = kp.nodes.DomainNodeGenerator();
generator.generateInteriorNodesFromGeometry(surface, 0.08, 'Seed', 17, 'StripCount', 5);
fprintf('Clipped interior nodes: %d from raw box nodes: %d\n', ...
    size(generator.getInteriorNodes(), 1), size(generator.getRawPoissonInteriorNodes(), 1));
fprintf('Boundary clearance used in interior generation: %.3f\n', 0.08);

fprintf('Running explicit parallel clipping example...\n');
[xin, ~, ~] = kp.nodes.clipPointsByGeometry(generator.getRawPoissonInteriorNodes(), ...
    surface, 'UseParallel', true, 'ChunkSize', 2000, 'MinParallelPoints', 1000, ...
    'BoundaryClearance', 0.08);
fprintf('Parallel-clipped interior nodes: %d\n', size(xin, 1));

fprintf('Running outer-refined interior-node example...\n');
refined = kp.nodes.DomainNodeGenerator();
refined.generateInteriorNodesFromGeometry(surface, 0.08, ...
    'Seed', 17, 'StripCount', 5, ...
    'DoOuterRefinement', true, ...
    'OuterFractionOfh', 0.5, ...
    'OuterRefinementZoneSizeAsMultipleOfh', 2.0);
fprintf('Outer-refined interior nodes: %d\n', size(refined.getInteriorNodes(), 1));

fprintf('Building DomainDescriptor from geometry...\n');
descriptor = refined.buildDomainDescriptorFromGeometry(surface, 0.08, ...
    'Seed', 17, 'StripCount', 5, ...
    'DoOuterRefinement', true, ...
    'OuterFractionOfh', 0.5, ...
    'OuterRefinementZoneSizeAsMultipleOfh', 2.0);
fprintf('Descriptor counts: Xi=%d, Xb=%d, Xg=%d\n', ...
    size(descriptor.getInteriorNodes(), 1), size(descriptor.getBdryNodes(), 1), size(descriptor.getGhostNodes(), 1));
