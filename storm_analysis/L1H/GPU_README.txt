
This was an attempt to further improve the computation speed by
using a GPU. It did not work, in that it turned out to be slower
than the our CPU version. The problem I believe is due to the
approach rather than any GPU specific issues. Profiling suggested
that the approach we followed involved too much data transfer
and too little computation. Basically the calculations we
were giving the GPU were too small so a lot of time got eaten up
shuffling data onto and off of the GPU. We were using the GPU to
handle the matrix / vector multiplications involved in solving a
single instance. A better approach would perhaps be to have each
GPU handle solving a single instance (i.e. solving for the
solution of a 7x7 pixel grid). Obviously this is much more
complicated to implement.

Additionally our experience was that the floating point math
that most GPUs are capable of is not accurate enough for L1H, and
you really need to use double precision math. This can be done
with some GPUs, but this issue combined with the likely
complexity of implementing what we thought was a better approach
put us off from any further attempts.

The relevant files are:
compile_gpu.bat
homotopy_gpu.c
homotopy_gpu.cl
homotopy_gpu.h

This depends on OpenCL.

We've added this in the hopes that others might find it useful
and we'd be happy to include any GPU related contributions to
this project.

Hazen Babcock
March 2014.

