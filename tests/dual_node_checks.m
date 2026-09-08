function dual_node_checks()
%DUAL_NODE_CHECKS Targeted checks for dual velocity/pressure domain setup.

t = linspace(0, 2*pi, 60).';
t(end) = [];
curve = [cos(t), 0.7 * sin(t)];
surface = kp.geometry.EmbeddedSurface();
surface.setDataSites(curve);
surface.buildClosedGeometricModelPS(2, 0.05, size(curve, 1));
surface.buildLevelSetFromGeometricModel([]);

dualGen = kp.nodes.DualNodeDomainGenerator();
pressureInfo = dualGen.generateSmoothDomainNodesAutoPressure(surface, 0.08, 0.35, 20, ...
    'Seed', 29, 'StripCount', 5, ...
    'DoOuterRefinement', true, ...
    'OuterFractionOfh', 0.5, ...
    'OuterRefinementZoneSizeAsMultipleOfh', 2.0);
dual = dualGen.createDualNodeDomainDescriptor();

assert(dual.hasVelocityDomain(), 'Dual descriptor should contain a velocity domain.');
assert(dual.hasPressureDomain(), 'Dual descriptor should contain a pressure domain.');
assert(dual.getVelocityDomain().getNumIntBdryNodes() == pressureInfo.num_velocity_nodes, ...
    'Velocity node count should match the auto-pressure report.');
assert(dual.getPressureDomain().getNumIntBdryNodes() == pressureInfo.num_pressure_nodes, ...
    'Pressure node count should match the auto-pressure report.');
assert(dual.getPressureDomain().getNumGhostNodes() == 0, ...
    'Pressure ghosts should be stripped by default in the dual descriptor.');
assert(dual.getPressureDomain().getNumIntBdryNodes() <= floor(0.35 * dual.getVelocityDomain().getNumIntBdryNodes()), ...
    'Pressure cloud should satisfy the requested density cap.');

disp('dual node checks passed');
end
