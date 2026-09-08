# 3D IBAMR RBC capstone

This external IBAMR application produces the material trajectory used by the
KernelPack moving-surface ADR capstone. It evolves a closed biconcave membrane
in a three-dimensional, streamwise/spanwise-periodic Poiseuille microchannel.
The membrane is never tethered: edge elasticity, hinge bending, and global
area/volume penalties are recomputed from the current midpoint geometry before
each IBAMR force assembly.

Configure against an IBAMR 0.19 installation:

```bash
cmake -S . -B build \
  -DIBAMR_ROOT=/path/to/IBAMR-0.19.0
cmake --build build -j 4
mpiexec -n 2 ./build/kernelpack_rbc_3d input3d
```

The `trajectory` directory contains fixed connectivity, reference material
coordinates, per-frame positions and velocities, and area/volume diagnostics.
The MATLAB capstone consumes these files offline; it does not modify the IBAMR
structural mesh.
