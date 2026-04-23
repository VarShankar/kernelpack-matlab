# kernelpack-matlab

`kernelpack-matlab` is a MATLAB companion project to
[KernelPack](https://github.com/VarShankar/kernelpack).

The aim is to bring the core geometry and meshfree building blocks of
KernelPack into a MATLAB codebase that is easier to inspect, prototype with,
and extend.

## Current focus

The repository is currently centered on the geometry layer. The first pieces in
place are modeled after the main geometry objects used in KernelPack:

- `EmbeddedSurface`
- `PiecewiseSmoothEmbeddedSurface`
- `RBFLevelSet`

These classes live in [`+kp/+geometry`](+kp/+geometry).

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

At the moment, the implemented geometric-model construction path is the 2D
case:

- smooth closed curves
- open curve segments
- piecewise-smooth planar boundaries assembled from segments

## Project direction

The goal is to keep building outward from the KernelPack geometry contracts
rather than inventing a separate MATLAB abstraction stack. The next major steps
are:

1. strengthen the current geometry objects
2. add the smooth closed 3D surface path
3. extend piecewise geometry handling further in 3D
4. build node-generation and discretization layers on top of that geometry

## Notes

This repository is under active construction. The current contents should be
viewed as an early geometry foundation rather than a full MATLAB port of
KernelPack.
