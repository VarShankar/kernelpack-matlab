# kernelpack-matlab

`kernelpack-matlab` is a MATLAB companion project to
[KernelPack](https://github.com/VarShankar/kernelpack).

The aim is to bring the core geometry and meshfree building blocks of
KernelPack into a MATLAB codebase that is easier to inspect, prototype with,
and extend.

## Current focus

The repository is currently centered on the geometry layer. The first pieces in
place are modeled after the main geometry and node-generation objects used in
KernelPack:

- `EmbeddedSurface`
- `PiecewiseSmoothEmbeddedSurface`
- `RBFLevelSet`
- `JacobiPolynomials`
- total-degree and hyperbolic-cross index helpers
- Chebyshev recurrence and evaluation helpers
- `DomainDescriptor`
- `DomainNodeGenerator`

These classes live in [`+kp/+geometry`](+kp/+geometry) and
[`+kp/+nodes`](+kp/+nodes) together with [`+kp/+domain`](+kp/+domain).

## What is here now

The current code establishes a KernelPack-shaped starting point for geometry:

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
- `kp.poly` also includes multi-index builders and related helpers such as
  `total_degree_indices`, `hyperbolic_cross_indices`, `chebyshev_recurrence`,
  `chebyshev_eval`, `ratio_eval`, and `hyphelper`.
- `DomainDescriptor` stores the pared-down domain state: interior nodes,
  boundary nodes, ghost nodes, boundary normals, and simple tree placeholders.
- `DomainNodeGenerator` provides seeded fixed-radius Poisson disk sampling on
  axis-aligned boxes, geometry-aware clipping from raw box clouds to interior
  node sets, and assembly into a `DomainDescriptor`.

At the moment, the implemented geometric-model construction path is the 2D
and early 3D cases:

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

## Project direction

The goal is to keep building outward from the KernelPack geometry contracts
rather than inventing a separate MATLAB abstraction stack. The next major steps
are:

1. strengthen the current geometry objects
2. improve the smooth closed 3D surface path further
3. extend piecewise geometry handling further in 3D
4. build node-generation and discretization layers on top of that geometry

## Notes

This repository is under active construction. The current contents should be
viewed as an early geometry foundation rather than a full MATLAB port of
KernelPack.
