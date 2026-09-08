function K = phsKernel(r, degree)
%PHSKERNEL Polyharmonic spline kernel.

K = r .^ degree;

