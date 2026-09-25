#include "fdtd_cuda_internal.h"

#include <cuda_runtime.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void fdtd_cuda_free_clones(fdtd_cuda_ctx *ctx)
{
   if (!ctx) return;
   if (ctx->h_clone_jobs) {
      free(ctx->h_clone_jobs);
      ctx->h_clone_jobs = NULL;
   }
   ctx->clone_n = 0;
   ctx->clone_ready = 0;
}

/*
 * One ghost plane. sign=-1 is PMC (negated adjacent plane); sign=+1 is periodic
 * (copy of the opposite interior plane). Free axes a,b are absolute indices.
 */
__global__ void k_clone(fdtd_real *__restrict__ F, int nx, int ny,
                        int wall, int ghost, int source, int sign,
                        int a0, int b0, int na, int nb,
                        int ox, int oy, int oz)
{
   int ia = blockIdx.x * blockDim.x + threadIdx.x;
   int ib = blockIdx.y * blockDim.y + threadIdx.y;
   if (ia >= na || ib >= nb) return;
   int a = a0 + ia;
   int b = b0 + ib;
   int gi, gj, gk, si, sj, sk;
   if (wall == 0) {
      gi = ghost; gj = a; gk = b;
      si = source; sj = a; sk = b;
   } else if (wall == 1) {
      gj = ghost; gi = a; gk = b;
      sj = source; si = a; sk = b;
   } else {
      gk = ghost; gi = a; gj = b;
      sk = source; si = a; sj = b;
   }
   int ig = gi - ox, jg = gj - oy, kg = gk - oz;
   int is = si - ox, js = sj - oy, ks = sk - oz;
   F[fdtd_idx3(ig, jg, kg, nx, ny)] = (fdtd_real)sign * F[fdtd_idx3(is, js, ks, nx, ny)];
}

int fdtd_cuda_set_clone_jobs(fdtd_cuda_ctx *ctx, const fdtd_clone_job *jobs, int n)
{
   if (!ctx || n < 0) return 0;
   fdtd_cuda_free_clones(ctx);
   if (n == 0) {
      ctx->clone_ready = 1;
      return 1;
   }
   if (!jobs) return 0;
   ctx->h_clone_jobs = (fdtd_clone_job *)malloc((size_t)n * sizeof(fdtd_clone_job));
   if (!ctx->h_clone_jobs) return 0;
   memcpy(ctx->h_clone_jobs, jobs, (size_t)n * sizeof(fdtd_clone_job));
   ctx->clone_n = n;
   ctx->clone_ready = 1;
   return 1;
}

int fdtd_cuda_clone_ready(const fdtd_cuda_ctx *ctx)
{
   return ctx && ctx->clone_ready;
}

int fdtd_cuda_advance_clones(fdtd_cuda_ctx *ctx)
{
   if (!ctx || !ctx->clone_ready) return 0;
   for (int j = 0; j < ctx->clone_n; ++j) {
      const fdtd_clone_job job = ctx->h_clone_jobs[j];
      fdtd_real *F = fdtd_cuda_field_ptr(ctx, job.field_comp);
      fdtd_dims3 d = fdtd_cuda_field_dim(ctx, job.field_comp);
      int na = job.a1 - job.a0 + 1;
      int nb = job.b1 - job.b0 + 1;
      if (!F || na < 1 || nb < 1) return 0;
      dim3 block(16, 16);
      dim3 grid((na + block.x - 1) / block.x, (nb + block.y - 1) / block.y);
      k_clone<<<grid, block>>>(F, d.nx, d.ny, job.wall_axis, job.ghost, job.source, job.sign,
                               job.a0, job.b0, na, nb, job.e_xi, job.e_yi, job.e_zi);
      if (cudaGetLastError() != cudaSuccess) return 0;
   }
   return cudaDeviceSynchronize() == cudaSuccess;
}
