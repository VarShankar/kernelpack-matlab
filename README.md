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
- total-degree and hyperbolic-cross index helpers
- Chebyshev recurrence and evaluation helpers

These classes live in [`+kp/+geometry`](+kp/+geometry) and
[`+kp/+nodes`](+kp/+nodes) together with [`+kp/+domain`](+kp/+domain),
[`+kp/+poly`](+kp/+poly), and [`+kp/+rbffd`](+kp/+rbffd).

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
- `DomainDescriptor` stores the pared-down domain state: interior nodes,
  boundary nodes, ghost nodes, boundary normals, tree placeholders, and the
  separation radius used by the current node set.
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
- operator metadata via `StencilProperties` and `OpProperties`

## Basic use

The repository also includes:

- [`examples/geometry_examples.m`](examples/geometry_examples.m) for a quick run
  through the current 2D and 3D geometry paths
- [`tests/geometry_checks.m`](tests/geometry_checks.m) for lightweight geometry
  checks on normals, assembled boundary clouds, and piecewise seam handling
- [`examples/nodes_examples.m`](examples/nodes_examples.m) for seeded box
  Poisson node generation examples, including parallel clipping
- [`tests/nodes_checks.m`](tests/nodes_checks.m) for basic node-generation
  determinism and spacing checks
- [`tests/poly_checks.m`](tests/poly_checks.m) for shared polynomial checks
- [`tests/rbffd_checks.m`](tests/rbffd_checks.m) for the current stencil and
  assembler checks

### Seeded box Poisson nodes

```matlab
[x, info] = kp.nodes.generatePoissonNodesInBox(0.08, [0 0 0], [1 1 1], ...
    'Seed', 17, 'StripCount', 5);

generator = kp.nodes.DomainNodeGenerator();
generator.generatePoissonNodes(0.08, [0 0 0], [1 1 1], 'Seed', 17, 'StripCount', 5);

raw_nodes = generator.getRawPoissonInteriorNodes();
```

### Geometry-clipped interior nodes

```matlab
t = linspace(0, 2*pi, 50).';
t(end) = [];
curve = [cos(t), 0.7*sin(t)];

surface = kp.geometry.EmbeddedSurface();
surface.setDataSites(curve);
surface.buildClosedGeometricModelPS(2, 0.05, size(curve,1), 120);
surface.buildLevelSetFromGeometricModel([]);

generator = kp.nodes.DomainNodeGenerator();
generator.generateInteriorNodesFromGeometry(surface, 0.08, ...
    'Seed', 17, 'StripCount', 5);

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
generator = kp.nodes.DomainNodeGenerator();
generator.generateInteriorNodesFromGeometry(surface, 0.08, ...
    'Seed', 17, 'StripCount', 5, ...
    'DoOuterRefinement', true, ...
    'OuterFractionOfh', 0.5, ...
    'OuterRefinementZoneSizeAsMultipleOfh', 2.0);
```

For larger clouds, clipping can also evaluate the level set in parallel:

```matlab
[interior_nodes, keep_mask, phi] = kp.nodes.clipPointsByGeometry( ...
    raw_box_nodes, surface, 'UseParallel', true, ...
    'ChunkSize', 5000, 'MinParallelPoints', 20000, ...
    'BoundaryClearance', 0.08);
```

### Smooth closed curve

```matlab
t = linspace(0, 2*pi, 40).';
t(end) = [];
x = [cos(t), 0.7*sin(t)];

surface = kp.geometry.EmbeddedSurface();
surface.setDataSites(x);
surface.buildClosedGeometricModelPS(2, 0.05, size(x,1), 120);
surface.buildLevelSetFromGeometricModel([]);

xb = surface.getSampleSites();
nrmls = surface.getNrmls();
phi = surface.getLevelSet().Evaluate(xb);
```

### Open curve segment

```matlab
u = linspace(0, 1, 30).';
x = [u, u.^2];

segment = kp.geometry.EmbeddedSurface();
segment.setDataSites(x);
segment.buildGeometricModelPS(2, 0.05, size(x,1), 80);

xb = segment.getSampleSites();
nrmls = segment.getNrmls();
```

### Piecewise-smooth planar boundary

```matlab
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
```

### Smooth closed surface in 3D

```matlab
X = kp.geometry.fibonacciSphere(120);
uv = kp.geometry.cart2sphRows(X);
r = 1 + 0.15*cos(3*uv(:,1)).*cos(2*uv(:,2));
pts = X .* r;

surface = kp.geometry.EmbeddedSurface();
surface.setDataSites(pts);
surface.buildClosedGeometricModelPS(3, 0.2, size(pts,1), 160);
surface.buildLevelSetFromGeometricModel([]);

xb = surface.getSampleSites();
nrmls = surface.getNrmls();
```

### RBF-FD assembly on a descriptor

```matlab
[Xg, Yg] = ndgrid(linspace(-1, 1, 5), linspace(-1, 1, 5));
X = [Xg(:), Yg(:)];

domain = kp.domain.DomainDescriptor();
domain.setNodes(X, zeros(0, 2), zeros(0, 2));
domain.setSepRad(0.5);
domain.buildStructs();

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
