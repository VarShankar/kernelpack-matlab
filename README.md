# kernelpack-matlab

`kernelpack-matlab` is a MATLAB companion project to
[KernelPack](https://github.com/VarShankar/kernelpack).

The aim is to bring the core geometry, node-generation, polynomial, and
meshfree discretization ingredients of KernelPack into a MATLAB codebase that
is easier to inspect, prototype with, and extend.

## Current focus

The repository now has a first pass of the main layers that sit underneath a
KernelPack-style RBF-FD workflow:

- `EmbeddedSurface`
- `PiecewiseSmoothEmbeddedSurface`
- `RBFLevelSet`
- `DomainNodeGenerator`
- `DomainDescriptor`
- `JacobiPolynomials`
- `PolynomialBasis`
- `RBFStencil`
- `WeightedLeastSquaresStencil`
- `FDDiffOp`
- `FDODiffOp`
- `PoissonSolver`
- `DiffusionSolver`
- total-degree and hyperbolic-cross index helpers
- Chebyshev recurrence and evaluation helpers

These classes live in [`+kp/+geometry`](+kp/+geometry) and
[`+kp/+nodes`](+kp/+nodes) together with [`+kp/+domain`](+kp/+domain),
[`+kp/+poly`](+kp/+poly), [`+kp/+rbffd`](+kp/+rbffd), and
[`+kp/+solvers`](+kp/+solvers).

## What is here now

The current codebase includes:

- `EmbeddedSurface` stores the data sites, sampled boundary points, normals,
  tangents, bounding boxes, thickened copies, and level-set representation for
  a single boundary object.
- `PiecewiseSmoothEmbeddedSurface` stores the segment list together with the
  assembled boundary cloud, normals, segment map, corner flags, bounding boxes,
  and level-set representation for a piecewise-smooth boundary.
- `RBFLevelSet` provides an implicit boundary representation with evaluation,
  gradient evaluation, inside-outside tests, and Newton projection routines.
- `JacobiPolynomials` provides shared Jacobi recurrence, univariate orthonormal
  polynomial evaluation, derivatives, and tensor-product evaluation.
- `PolynomialBasis` provides a higher-level basis object. It defaults to
  Legendre polynomials and stores a center and scalar scale factor for mapping
  stencil points into a unit disk or sphere before evaluation.
- `kp.poly` also includes multi-index builders and related helpers such as
  `total_degree_indices`, `hyperbolic_cross_indices`, `chebyshev_recurrence`,
  `chebyshev_eval`, `ratio_eval`, and `hyphelper`.
- `kp.rbffd` now includes KernelPack-shaped local stencil and assembler
  classes for RBF-FD and polynomial weighted least-squares work:
  `RBFStencil`, `WeightedLeastSquaresStencil`, `FDDiffOp`, `FDODiffOp`,
  `StencilProperties`, and `OpProperties`.
- `PoissonSolver` mirrors KernelPack's fixed-domain Poisson solver with cached
  Laplacian assembly, separate boundary assembly, and callback-driven solves.
- `DiffusionSolver` mirrors KernelPack's fixed-domain diffusion stepper with
  cached Laplacian assembly, state history, and BDF1/BDF2/BDF3 time steps.
- `DomainDescriptor` stores the pared-down domain state: interior nodes,
  boundary nodes, ghost nodes, boundary normals, KD-tree-backed search
  structures for the available node sets, and the separation radius used by
  the current node set.
- `DomainNodeGenerator` provides seeded fixed-radius Poisson disk sampling on
  axis-aligned boxes, geometry-aware clipping from raw box clouds to interior
  node sets, outer refinement near boundaries, and assembly into a
  `DomainDescriptor`.

At the moment, the implemented geometry and node-generation paths include:

- smooth closed curves
- open curve segments
- piecewise-smooth planar boundaries assembled from segments
- smooth closed 3D surfaces built from spherical-coordinate SBF interpolation
- open 3D surface patches built from best-fit-plane chart coordinates
- piecewise 3D surfaces assembled from patch segments
- fixed-radius Poisson disk sampling on axis-aligned boxes in any dimension
- level-set clipping of box Poisson clouds against `EmbeddedSurface` and
  `PiecewiseSmoothEmbeddedSurface`
- chunked `parfor` evaluation for large geometry-clipping jobs
- outer refinement bands near the boundary using a smaller local spacing

The current `rbffd` layer includes:

- local PHS-plus-polynomial `RBFStencil` construction
- local Legendre weighted least-squares reconstruction in
  `WeightedLeastSquaresStencil`
- one-row-per-center assembly through `FDDiffOp`
- overlapped row acceptance and assembly through `FDODiffOp`
- tree-backed nearest-neighbor selection through `DomainDescriptor`
- operator metadata via `StencilProperties` and `OpProperties`

## Basic use

The repository also includes:

- [`examples/geometry_examples.m`](examples/geometry_examples.m) for a quick run
  through the current 2D and 3D geometry paths, including plots of the
  source clouds, boundary samples, and 3D triangulated surfaces
- [`tests/geometry_checks.m`](tests/geometry_checks.m) for lightweight geometry
  checks on normals, assembled boundary clouds, and piecewise seam handling
- [`examples/nodes_examples.m`](examples/nodes_examples.m) for seeded box
  Poisson node generation examples, including parallel clipping
- [`tests/nodes_checks.m`](tests/nodes_checks.m) for basic node-generation
  determinism and spacing checks
- [`tests/poly_checks.m`](tests/poly_checks.m) for shared polynomial checks
- [`tests/rbffd_checks.m`](tests/rbffd_checks.m) for the current stencil and
  assembler checks
- [`examples/poisson_solver_example.m`](examples/poisson_solver_example.m) for
  a full Poisson solve on a geometry-built domain
- [`tests/poisson_solver_checks.m`](tests/poisson_solver_checks.m) for
  lightweight Poisson solver checks
- [`examples/diffusion_solver_example.m`](examples/diffusion_solver_example.m)
  for a fixed-domain diffusion time-stepping run on a geometry-built domain
- [`tests/diffusion_solver_checks.m`](tests/diffusion_solver_checks.m) for
  lightweight diffusion solver checks

### Seeded box Poisson nodes

```matlab
% Generate a seeded fixed-radius Poisson cloud directly in a box.
[x, info] = kp.nodes.generatePoissonNodesInBox(0.08, [0 0 0], [1 1 1], ...
    'Seed', 17, 'StripCount', 5);

% Use the higher-level generator wrapper for the same task.
generator = kp.nodes.DomainNodeGenerator();
generator.generatePoissonNodes(0.08, [0 0 0], [1 1 1], 'Seed', 17, 'StripCount', 5);

% Retrieve the raw interior cloud from the generator state.
raw_nodes = generator.getRawPoissonInteriorNodes();
```

### Geometry-clipped interior nodes

```matlab
% Define a smooth closed planar curve from scattered boundary data.
t = linspace(0, 2*pi, 50).';
t(end) = [];
curve = [cos(t), 0.7*sin(t)];

% Build the stacked parametric-plus-level-set representation.
surface = kp.geometry.EmbeddedSurface();
surface.setDataSites(curve);
surface.buildClosedGeometricModelPS(2, 0.05, size(curve,1));
surface.buildLevelSetFromGeometricModel([]);

% Generate a box cloud and clip it against the geometry.
generator = kp.nodes.DomainNodeGenerator();
generator.generateInteriorNodesFromGeometry(surface, 0.08, ...
    'Seed', 17, 'StripCount', 5);

% Pull out the clipped interior nodes, the raw box cloud, or a full descriptor.
interior_nodes = generator.getInteriorNodes();
raw_box_nodes = generator.getRawPoissonInteriorNodes();
descriptor = generator.buildDomainDescriptorFromGeometry(surface, 0.08, ...
    'Seed', 17, 'StripCount', 5);
```

The interior generation path keeps a clearance band of width `h` from the
boundary, so nodes within one local spacing of the level set are removed.

Near-boundary outer refinement is also available. In that mode the coarse
interior cloud stays outside the refinement band, a finer cloud is generated
inside the band, and the level-set clearance is based on the smallest radius
being used.

```matlab
% Enable a finer spacing in a band near the boundary.
generator = kp.nodes.DomainNodeGenerator();
generator.generateInteriorNodesFromGeometry(surface, 0.08, ...
    'Seed', 17, 'StripCount', 5, ...
    'DoOuterRefinement', true, ...
    'OuterFractionOfh', 0.5, ...
    'OuterRefinementZoneSizeAsMultipleOfh', 2.0);
```

For larger clouds, clipping can also evaluate the level set in parallel:

```matlab
% Evaluate the level set in chunks and farm large jobs out to parfor.
[interior_nodes, keep_mask, phi] = kp.nodes.clipPointsByGeometry( ...
    raw_box_nodes, surface, 'UseParallel', true, ...
    'ChunkSize', 5000, 'MinParallelPoints', 20000, ...
    'BoundaryClearance', 0.08);
```

### Smooth closed curve

```matlab
% Start from data sites on a smooth closed planar boundary.
t = linspace(0, 2*pi, 40).';
t(end) = [];
x = [cos(t), 0.7*sin(t)];

% Build the geometric model and its matching level set.
surface = kp.geometry.EmbeddedSurface();
surface.setDataSites(x);
surface.buildClosedGeometricModelPS(2, 0.05, size(x,1));
surface.buildLevelSetFromGeometricModel([]);

% Extract boundary samples, normals, and the implicit values at those samples.
xb = surface.getSampleSites();
nrmls = surface.getNrmls();
phi = surface.getLevelSet().Evaluate(xb);

% Plot the input sites, the sampled boundary, and the sampled normals.
figure('Color', 'w');
tiledlayout(1, 3);

nexttile;
plot(x(:,1), x(:,2), 'ko');
axis equal;
grid on;
title('Data sites');

nexttile;
plot(xb(:,1), xb(:,2), 'b.');
axis equal;
grid on;
title('Boundary samples');

nexttile;
plot(xb(:,1), xb(:,2), 'b.');
hold on;
idx = 1:size(xb,1);
tips = xb(idx,:) + 0.035*nrmls(idx,:);
for k = 1:numel(idx)
    line([xb(idx(k),1), tips(k,1)], [xb(idx(k),2), tips(k,2)], ...
        'Color', [0.85 0.2 0.2], 'LineWidth', 1);
end
axis equal;
grid on;
title('Boundary-sample normals');
```

### Open curve segment

```matlab
% Sample an open planar segment and build a nonperiodic parametric model.
u = linspace(0, 1, 30).';
x = [u, u.^2];

segment = kp.geometry.EmbeddedSurface();
segment.setDataSites(x);
segment.buildGeometricModelPS(2, 0.05, size(x,1));

xb = segment.getSampleSites();
nrmls = segment.getNrmls();
```

### Piecewise-smooth planar boundary

```matlab
% Build a closed piecewise-smooth boundary from four segments.
seg1 = [linspace(0,1,20).', zeros(20,1)];
seg2 = [ones(20,1), linspace(0,1,20).'];
seg3 = [linspace(1,0,20).', ones(20,1)];
seg4 = [zeros(20,1), linspace(1,0,20).'];

surface = kp.geometry.PiecewiseSmoothEmbeddedSurface();
surface.generatePiecewiseSmoothSurfaceBySegment( ...
    {seg1, seg2, seg3, seg4}, [false false false false], 0.05, 1, 2, 2);
surface.buildLevelSet();

xb = surface.getBdryNodes();
nrmls = surface.getBdryNrmls();
corners = surface.getCornerFlags();
cornerMask = logical(corners);

% Plot the original segments, the assembled boundary cloud, and the corner-aware normals.
figure('Color', 'w');
tiledlayout(1, 3);

nexttile;
plot(seg1(:,1), seg1(:,2), 'k.-');
hold on;
plot(seg2(:,1), seg2(:,2), 'k.-');
plot(seg3(:,1), seg3(:,2), 'k.-');
plot(seg4(:,1), seg4(:,2), 'k.-');
axis equal;
grid on;
title('Input segments');

nexttile;
plot(xb(:,1), xb(:,2), 'b.');
hold on;
plot(xb(cornerMask,1), xb(cornerMask,2), 'mo');
axis equal;
grid on;
title('Assembled boundary');

nexttile;
plot(xb(:,1), xb(:,2), 'b.');
hold on;
idx = 1:size(xb,1);
tips = xb(idx,:) + 0.035*nrmls(idx,:);
for k = 1:numel(idx)
    line([xb(idx(k),1), tips(k,1)], [xb(idx(k),2), tips(k,2)], ...
        'Color', [0.85 0.2 0.2], 'LineWidth', 1);
end
plot(xb(cornerMask,1), xb(cornerMask,2), 'mo');
axis equal;
grid on;
title('Boundary normals and corners');
```

### Smooth closed surface in 3D

```matlab
% Start from a smooth closed 3D point cloud.
X = kp.geometry.fibonacciSphere(120);
uv = kp.geometry.cart2sphRows(X);
r = 1 + 0.15*cos(3*uv(:,1)).*cos(2*uv(:,2));
pts = X .* r;

% Build the surface model and evaluate a boundary cloud on it.
surface = kp.geometry.EmbeddedSurface();
surface.setDataSites(pts);
surface.buildClosedGeometricModelPS(3, 0.2, size(pts,1));
surface.buildLevelSetFromGeometricModel([]);

xb = surface.getSampleSites();
nrmls = surface.getNrmls();

% Show the input cloud, the sampled boundary cloud, and a reconstructed triangulation.
figure('Color', 'w');
tiledlayout(1, 3);

nexttile;
plot3(pts(:,1), pts(:,2), pts(:,3), 'k.');
axis equal;
grid on;
view(3);
title('Data cloud');

nexttile;
plot3(xb(:,1), xb(:,2), xb(:,3), 'b.');
axis equal;
grid on;
view(3);
title('Boundary cloud');

nexttile;
tri = kp.geometry.MyRobustCrust(xb);
trisurf(tri, xb(:,1), xb(:,2), xb(:,3), ...
    'FaceColor', [0.2 0.6 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
hold on;
plot3(xb(:,1), xb(:,2), xb(:,3), 'k.');
axis equal;
grid on;
view(3);
camlight headlight;
lighting gouraud;
title('Triangulated surface');
```

### RBF-FD assembly on a descriptor

```matlab
% Build a descriptor directly from a small Cartesian grid.
[Xg, Yg] = ndgrid(linspace(-1, 1, 5), linspace(-1, 1, 5));
X = [Xg(:), Yg(:)];

domain = kp.domain.DomainDescriptor();
domain.setNodes(X, zeros(0, 2), zeros(0, 2));
domain.setSepRad(0.5);
domain.buildStructs();

% Describe the stencil and assemble a Laplacian on the descriptor.
sp = kp.rbffd.StencilProperties( ...
    'n', 9, 'dim', 2, 'ell', 2, ...
    'spline_degree', 3, ...
    'treeMode', 'interior_boundary', ...
    'pointSet', 'interior_boundary');
op = kp.rbffd.OpProperties('recordStencils', true);

assembler = kp.rbffd.FDDiffOp(@() kp.rbffd.RBFStencil());
assembler.AssembleOp(domain, 'lap', sp, op);
L = assembler.getOp();
```

### Full RBF-FD workflow from geometry to operator

High-level stencil selection:

```matlab
% Build a smooth geometry and convert it into a full domain descriptor.
t = linspace(0, 2*pi, 80).';
t(end) = [];
curve = [cos(t), 0.7*sin(t)];

surface = kp.geometry.EmbeddedSurface();
surface.setDataSites(curve);
surface.buildClosedGeometricModelPS(2, 0.05, size(curve, 1));
surface.buildLevelSetFromGeometricModel([]);

generator = kp.nodes.DomainNodeGenerator();
domain = generator.buildDomainDescriptorFromGeometry(surface, 0.08, ...
    'Seed', 17, ...
    'StripCount', 5, ...
    'DoOuterRefinement', true, ...
    'OuterFractionOfh', 0.5, ...
    'OuterRefinementZoneSizeAsMultipleOfh', 2.0);

% Ask the code to choose stencil parameters from the target accuracy.
sp = kp.rbffd.StencilProperties.fromAccuracy( ...
    'Operator', 'lap', ...
    'ConvergenceOrder', 3, ...
    'Dimension', 2, ...
    'Approximation', 'rbf', ...
    'treeMode', 'all', ...
    'pointSet', 'interior_boundary');
op = kp.rbffd.OpProperties('recordStencils', true);

% Assemble the interior Laplacian.
lapAssembler = kp.rbffd.FDDiffOp(@() kp.rbffd.RBFStencil());
lapAssembler.AssembleOp(domain, 'lap', sp, op);
L = lapAssembler.getOp();

% Assemble a boundary operator with one row per boundary node.
bcSp = kp.rbffd.StencilProperties.fromAccuracy( ...
    'Operator', 'bc', ...
    'ConvergenceOrder', 3, ...
    'Dimension', 2, ...
    'Approximation', 'rbf', ...
    'treeMode', 'all', ...
    'pointSet', 'boundary');

bcAssembler = kp.rbffd.FDDiffOp(@() kp.rbffd.RBFStencil());
nb = domain.getNumBdryNodes();
bcAssembler.AssembleOp(domain, 'bc', bcSp, op, ...
    'NeuCoeff', zeros(nb, 1), ...
    'DirCoeff', ones(nb, 1));
BC = bcAssembler.getOp();
```

That path starts from a geometric model, generates interior, boundary, and
ghost nodes through `DomainNodeGenerator`, packs them into a
`DomainDescriptor`, and then assembles interior and boundary RBF-FD operators
from the descriptor. The interior operator `L` uses the descriptor's
interior-plus-boundary row space, while the boundary operator `BC` is stored
with one row per boundary point and one column per total node.

Low-level stencil selection:

```matlab
% Supply stencil size and polynomial degree explicitly when you want manual control.
sp = kp.rbffd.StencilProperties( ...
    'n', 25, ...
    'dim', 2, ...
    'ell', 4, ...
    'spline_degree', 5, ...
    'treeMode', 'all', ...
    'pointSet', 'interior_boundary');

bcSp = kp.rbffd.StencilProperties( ...
    'n', 25, ...
    'dim', 2, ...
    'ell', 4, ...
    'spline_degree', 5, ...
    'treeMode', 'all', ...
    'pointSet', 'boundary');
```

### Fixed-domain Poisson solve on a geometry-backed domain

```matlab
% Build a geometry and convert it into a DomainDescriptor.
t = linspace(0, 2*pi, 80).';
t(end) = [];
curve = [cos(t), 0.8*sin(t)];

surface = kp.geometry.EmbeddedSurface();
surface.setDataSites(curve);
surface.buildClosedGeometricModelPS(2, 0.06, size(curve, 1));
surface.buildLevelSetFromGeometricModel([]);

generator = kp.nodes.DomainNodeGenerator();
domain = generator.buildDomainDescriptorFromGeometry(surface, 0.1, ...
    'Seed', 17, ...
    'StripCount', 5, ...
    'DoOuterRefinement', true, ...
    'OuterFractionOfh', 0.5, ...
    'OuterRefinementZoneSizeAsMultipleOfh', 2.0);

% Set up the fixed-domain Poisson solver with separate Laplacian and BC pieces.
solver = kp.solvers.PoissonSolver( ...
    'LapAssembler', 'fd', ...
    'BCAssembler', 'fd', ...
    'LapStencil', 'wls', ...
    'BCStencil', 'wls');
solver.init(domain, 3);

% Solve -Delta u = f with callback pieces for the forcing, BC coefficients,
% and boundary values.
uExact = @(X) X(:,1).^2 + X(:,2).^2;
forcing = @(Xeq) -4 * ones(size(Xeq, 1), 1);
neuCoeff = @(Xb) zeros(size(Xb, 1), 1);
dirCoeff = @(Xb) ones(size(Xb, 1), 1);
bc = @(NeuCoeffs, DirCoeffs, nr, Xb) uExact(Xb);

result = solver.solve(forcing, neuCoeff, dirCoeff, bc);
u = result.u;
```

### End-to-end Poisson solve with pure Neumann data

```matlab
% Build a smooth closed domain from boundary data.
t = linspace(0, 2*pi, 120).';
t(end) = [];
curve = [cos(t), sin(t)];

surface = kp.geometry.EmbeddedSurface();
surface.setDataSites(curve);
surface.buildClosedGeometricModelPS(2, 0.06, size(curve,1));
surface.buildLevelSetFromGeometricModel([]);

% Generate interior, boundary, and ghost nodes from the geometry.
generator = kp.nodes.DomainNodeGenerator();
domain = generator.buildDomainDescriptorFromGeometry(surface, 0.08, ...
    'Seed', 17, 'StripCount', 5);

% Set up an RBF-FD Poisson solve on the assembled domain.
solver = kp.solvers.PoissonSolver( ...
    'LapAssembler', 'fd', ...
    'BCAssembler', 'fd', ...
    'LapStencil', 'rbf', ...
    'BCStencil', 'rbf');
solver.init(domain, 4);

% Use a pure-Neumann manufactured solution on the unit disk.
uExact = @(X) (X(:,1).^2 + X(:,2).^2).^2 - (X(:,1).^2 + X(:,2).^2) + 1/6;
forcing = @(X) 4 - 16*(X(:,1).^2 + X(:,2).^2);
neuCoeff = @(Xb) ones(size(Xb,1), 1);
dirCoeff = @(Xb) zeros(size(Xb,1), 1);
bc = @(NeuCoeffs, DirCoeffs, nr, Xb) ...
    sum(([4*Xb(:,1).*(Xb(:,1).^2 + Xb(:,2).^2) - 2*Xb(:,1), ...
          4*Xb(:,2).*(Xb(:,1).^2 + Xb(:,2).^2) - 2*Xb(:,2)]).*nr, 2);

% Solve the square ghost-node system with the usual nullspace augmentation.
result = solver.solve(forcing, neuCoeff, dirCoeff, bc);

% Pure Neumann solutions are only defined up to a constant, so align the means.
Xphys = domain.getIntBdryNodes();
u = result.u;
uTrue = uExact(Xphys);
u = u - mean(u - uTrue);
```

The same solver can be switched to a different Laplacian or boundary backend
just by changing the constructor options:

```matlab
solver = kp.solvers.PoissonSolver( ...
    'LapAssembler', 'fd', ...
    'BCAssembler', 'fd', ...
    'LapStencil', 'rbf', ...
    'BCStencil', 'rbf');
```

Like the C++ solver, `PoissonSolver.solve(...)` returns the physical target
state on the interior-plus-boundary cloud, not the full all-node vector.

### Fixed-domain diffusion stepping on a geometry-backed domain

```matlab
% Set up the fixed-domain diffusion stepper with separate Laplacian and BC pieces.
solver = kp.solvers.DiffusionSolver( ...
    'LapAssembler', 'fd', ...
    'BCAssembler', 'fd', ...
    'LapStencil', 'wls', ...
    'BCStencil', 'wls');

nu = 0.25;
dt = 0.02;
solver.init(domain, 3, dt, nu);

% Seed the state history with the physical state on the interior-plus-boundary cloud.
uExact = @(time, X) exp(-time) .* (X(:,1).^2 + X(:,2).^2);
forcing = @(nuValue, time, X) -exp(-time) .* (X(:,1).^2 + X(:,2).^2) ...
    - 4 * nuValue * exp(-time);
neuCoeff = @(Xb) zeros(size(Xb, 1), 1);
dirCoeff = @(Xb) ones(size(Xb, 1), 1);
bc = @(NeuCoeffs, DirCoeffs, nr, time, Xb) uExact(time, Xb);

solver.setInitialState(uExact(0, domain.getIntBdryNodes()));
u1 = solver.bdf1Step(dt, forcing, neuCoeff, dirCoeff, bc);
u2 = solver.bdf2Step(2 * dt, forcing, neuCoeff, dirCoeff, bc);
u3 = solver.bdf3Step(3 * dt, forcing, neuCoeff, dirCoeff, bc);
```

## Project direction

The goal is still to keep building outward from the KernelPack contracts
rather than inventing a separate MATLAB abstraction stack. The next major steps
are now:

1. strengthen the smooth closed 3D and piecewise 3D geometry paths further
2. make the node-generation layer richer around boundary and ghost-node logic
3. expand the `rbffd` operator library beyond the current interpolation,
   gradient, Laplacian, and boundary-condition paths
4. add higher-level solver workflows on top of `DomainDescriptor` and `rbffd`

## Notes

This repository is under active construction. The current contents should be
viewed as an early but growing MATLAB implementation of the main KernelPack
ingredients rather than a full port of every KernelPack path.
