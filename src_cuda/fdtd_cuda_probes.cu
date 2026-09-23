#include "fdtd_cuda_internal.h"

#include <cuda_runtime.h>
#include <stdio.h>
#include <stdlib.h>

static int check_pr(cudaError_t e, const char *what)
{
   if (e != cudaSuccess) {
      fprintf(stderr, "FDTD CUDA error in %s: %s\n", what, cudaGetErrorString(e));
      return 0;
   }
   return 1;
}

static void free_probe_bufs(fdtd_cuda_ctx *ctx)
{
#define FREE(p) do { if (p) cudaFree(p); p = NULL; } while (0)
   FREE(ctx->d_probe_comp);
   FREE(ctx->d_probe_i);
   FREE(ctx->d_probe_j);
   FREE(ctx->d_probe_k);
   FREE(ctx->d_probe_out);
#undef FREE
   ctx->n_probes = 0;
}

static int validate_probe(fdtd_cuda_ctx *ctx, int comp, int i, int j, int k)
{
   if (comp < 0 || comp > 5) return 0;
   fdtd_dims3 d = fdtd_cuda_field_dim(ctx, comp);
   if (d.nx <= 0 || d.ny <= 0 || d.nz <= 0) return 0;
   if (i < 0 || i >= d.nx || j < 0 || j >= d.ny || k < 0 || k >= d.nz) return 0;
   return 1;
}

static __global__ void k_gather_probes(
   const fdtd_real *__restrict__ Ex, const fdtd_real *__restrict__ Ey, const fdtd_real *__restrict__ Ez,
   const fdtd_real *__restrict__ Hx, const fdtd_real *__restrict__ Hy, const fdtd_real *__restrict__ Hz,
   int ex_nx, int ex_ny, int ey_nx, int ey_ny, int ez_nx, int ez_ny,
   int hx_nx, int hx_ny, int hy_nx, int hy_ny, int hz_nx, int hz_ny,
   const int *__restrict__ comp, const int *__restrict__ ii, const int *__restrict__ jj,
   const int *__restrict__ kk, fdtd_real *__restrict__ out, int n)
{
   int t = blockIdx.x * blockDim.x + threadIdx.x;
   if (t >= n) return;
   int c = comp[t];
   int i = ii[t], j = jj[t], k = kk[t];
   fdtd_real v = 0;
   switch (c) {
   case 0: v = Ex[fdtd_idx3(i, j, k, ex_nx, ex_ny)]; break;
   case 1: v = Ey[fdtd_idx3(i, j, k, ey_nx, ey_ny)]; break;
   case 2: v = Ez[fdtd_idx3(i, j, k, ez_nx, ez_ny)]; break;
   case 3: v = Hx[fdtd_idx3(i, j, k, hx_nx, hx_ny)]; break;
   case 4: v = Hy[fdtd_idx3(i, j, k, hy_nx, hy_ny)]; break;
   case 5: v = Hz[fdtd_idx3(i, j, k, hz_nx, hz_ny)]; break;
   default: break;
   }
   out[t] = v;
}

int fdtd_cuda_set_point_probes(fdtd_cuda_ctx *ctx, int n,
                               const int *comp, const int *i, const int *j, const int *k)
{
   if (!ctx) return 0;
   free_probe_bufs(ctx);
   if (n < 0) return 0;
   if (n == 0) return 1;
   if (!comp || !i || !j || !k) return 0;
   if (!ctx->d_Ex) return 0; /* fields must be allocated */

   for (int t = 0; t < n; ++t) {
      if (!validate_probe(ctx, comp[t], i[t], j[t], k[t])) return 0;
   }

   size_t ni = (size_t)n * sizeof(int);
   size_t nr = (size_t)n * sizeof(fdtd_real);
   if (!check_pr(cudaMalloc((void **)&ctx->d_probe_comp, ni), "probe_comp") ||
       !check_pr(cudaMalloc((void **)&ctx->d_probe_i, ni), "probe_i") ||
       !check_pr(cudaMalloc((void **)&ctx->d_probe_j, ni), "probe_j") ||
       !check_pr(cudaMalloc((void **)&ctx->d_probe_k, ni), "probe_k") ||
       !check_pr(cudaMalloc((void **)&ctx->d_probe_out, nr), "probe_out")) {
      free_probe_bufs(ctx);
      return 0;
   }
   if (!check_pr(cudaMemcpy(ctx->d_probe_comp, comp, ni, cudaMemcpyHostToDevice), "H2D comp") ||
       !check_pr(cudaMemcpy(ctx->d_probe_i, i, ni, cudaMemcpyHostToDevice), "H2D i") ||
       !check_pr(cudaMemcpy(ctx->d_probe_j, j, ni, cudaMemcpyHostToDevice), "H2D j") ||
       !check_pr(cudaMemcpy(ctx->d_probe_k, k, ni, cudaMemcpyHostToDevice), "H2D k")) {
      free_probe_bufs(ctx);
      return 0;
   }
   ctx->n_probes = n;
   return 1;
}

int fdtd_cuda_gather_point_probes(fdtd_cuda_ctx *ctx, fdtd_real *out)
{
   if (!ctx) return 0;
   if (ctx->n_probes == 0) return 1;
   if (!out || !ctx->d_probe_out) return 0;

   fdtd_dims3 ex = ctx->dex, ey = ctx->dey, ez = ctx->dez;
   fdtd_dims3 hx = ctx->dhx, hy = ctx->dhy, hz = ctx->dhz;
   int n = ctx->n_probes;
   int threads = 256;
   int blocks = (n + threads - 1) / threads;
   k_gather_probes<<<blocks, threads>>>(
      ctx->d_Ex, ctx->d_Ey, ctx->d_Ez, ctx->d_Hx, ctx->d_Hy, ctx->d_Hz,
      ex.nx, ex.ny, ey.nx, ey.ny, ez.nx, ez.ny,
      hx.nx, hx.ny, hy.nx, hy.ny, hz.nx, hz.ny,
      ctx->d_probe_comp, ctx->d_probe_i, ctx->d_probe_j, ctx->d_probe_k,
      ctx->d_probe_out, n);
   if (!check_pr(cudaGetLastError(), "gather kernel") ||
       !check_pr(cudaDeviceSynchronize(), "gather sync"))
      return 0;
   return check_pr(cudaMemcpy(out, ctx->d_probe_out, (size_t)n * sizeof(fdtd_real),
                              cudaMemcpyDeviceToHost),
                   "gather D2H");
}
