#include "fdtd_cuda_internal.h"

#include <cuda_runtime.h>

static __global__ void k_advance_ex(
   fdtd_real *__restrict__ Ex, const fdtd_real *__restrict__ Hy, const fdtd_real *__restrict__ Hz,
   const fdtd_media *__restrict__ Mi, const fdtd_real *__restrict__ g1, const fdtd_real *__restrict__ g2,
   const fdtd_real *__restrict__ Idyh, const fdtd_real *__restrict__ Idzh,
   int nx, int ny, int nz, int is, int ie, int js, int je, int ks, int ke)
{
   int i = is + blockIdx.x * blockDim.x + threadIdx.x;
   int j = js + blockIdx.y * blockDim.y + threadIdx.y;
   int k = ks + blockIdx.z * blockDim.z + threadIdx.z;
   if (i > ie || j > je || k > ke) return;
   size_t id = fdtd_idx3(i, j, k, nx, ny);
   fdtd_media medio = Mi[id];
   fdtd_real Idzhk = Idzh[k];
   fdtd_real Idyhj = Idyh[j];
   size_t id_jm = fdtd_idx3(i, j - 1, k, nx, ny);
   size_t id_km = fdtd_idx3(i, j, k - 1, nx, ny);
   /* Hz uses Hz dims — caller ensures Hx/Hy/Hz dims match stencil access;
      here Hz/Hy share Ex spatial indexing in the remapped 0-based view only
      if Alloc dims match. We receive Hy/Hz with their own dims via separate
      pitches — for Yee, Fortran remaps each to its own 0:N*-1, and stencil
      uses the same (i,j,k) into each remapped array. So we need Hy/Hz pitches. */
   (void)nz;
   Ex[id] = g1[medio] * Ex[id] + g2[medio] *
            ((Hz[id] - Hz[id_jm]) * Idyhj - (Hy[id] - Hy[id_km]) * Idzhk);
}

/* Hy and Hz may have different nx,ny than Ex. Pass hy/hz pitches. */
static __global__ void k_advance_ex2(
   fdtd_real *__restrict__ Ex, const fdtd_real *__restrict__ Hy, const fdtd_real *__restrict__ Hz,
   const fdtd_media *__restrict__ Mi, const fdtd_real *__restrict__ g1, const fdtd_real *__restrict__ g2,
   const fdtd_real *__restrict__ Idyh, const fdtd_real *__restrict__ Idzh,
   int nx, int ny,
   int hy_nx, int hy_ny, int hz_nx, int hz_ny,
   int is, int ie, int js, int je, int ks, int ke)
{
   int i = is + blockIdx.x * blockDim.x + threadIdx.x;
   int j = js + blockIdx.y * blockDim.y + threadIdx.y;
   int k = ks + blockIdx.z * blockDim.z + threadIdx.z;
   if (i > ie || j > je || k > ke) return;
   size_t id = fdtd_idx3(i, j, k, nx, ny);
   fdtd_media medio = Mi[id];
   fdtd_real curl =
      (Hz[fdtd_idx3(i, j, k, hz_nx, hz_ny)] - Hz[fdtd_idx3(i, j - 1, k, hz_nx, hz_ny)]) * Idyh[j] -
      (Hy[fdtd_idx3(i, j, k, hy_nx, hy_ny)] - Hy[fdtd_idx3(i, j, k - 1, hy_nx, hy_ny)]) * Idzh[k];
   Ex[id] = g1[medio] * Ex[id] + g2[medio] * curl;
}

static __global__ void k_advance_ey2(
   fdtd_real *__restrict__ Ey, const fdtd_real *__restrict__ Hz, const fdtd_real *__restrict__ Hx,
   const fdtd_media *__restrict__ Mi, const fdtd_real *__restrict__ g1, const fdtd_real *__restrict__ g2,
   const fdtd_real *__restrict__ Idzh, const fdtd_real *__restrict__ Idxh,
   int nx, int ny, int hz_nx, int hz_ny, int hx_nx, int hx_ny,
   int is, int ie, int js, int je, int ks, int ke)
{
   int i = is + blockIdx.x * blockDim.x + threadIdx.x;
   int j = js + blockIdx.y * blockDim.y + threadIdx.y;
   int k = ks + blockIdx.z * blockDim.z + threadIdx.z;
   if (i > ie || j > je || k > ke) return;
   size_t id = fdtd_idx3(i, j, k, nx, ny);
   fdtd_media medio = Mi[id];
   fdtd_real curl =
      (Hx[fdtd_idx3(i, j, k, hx_nx, hx_ny)] - Hx[fdtd_idx3(i, j, k - 1, hx_nx, hx_ny)]) * Idzh[k] -
      (Hz[fdtd_idx3(i, j, k, hz_nx, hz_ny)] - Hz[fdtd_idx3(i - 1, j, k, hz_nx, hz_ny)]) * Idxh[i];
   Ey[id] = g1[medio] * Ey[id] + g2[medio] * curl;
}

static __global__ void k_advance_ez2(
   fdtd_real *__restrict__ Ez, const fdtd_real *__restrict__ Hx, const fdtd_real *__restrict__ Hy,
   const fdtd_media *__restrict__ Mi, const fdtd_real *__restrict__ g1, const fdtd_real *__restrict__ g2,
   const fdtd_real *__restrict__ Idyh, const fdtd_real *__restrict__ Idxh,
   int nx, int ny, int hx_nx, int hx_ny, int hy_nx, int hy_ny,
   int is, int ie, int js, int je, int ks, int ke)
{
   int i = is + blockIdx.x * blockDim.x + threadIdx.x;
   int j = js + blockIdx.y * blockDim.y + threadIdx.y;
   int k = ks + blockIdx.z * blockDim.z + threadIdx.z;
   if (i > ie || j > je || k > ke) return;
   size_t id = fdtd_idx3(i, j, k, nx, ny);
   fdtd_media medio = Mi[id];
   fdtd_real curl =
      (Hy[fdtd_idx3(i, j, k, hy_nx, hy_ny)] - Hy[fdtd_idx3(i - 1, j, k, hy_nx, hy_ny)]) * Idxh[i] -
      (Hx[fdtd_idx3(i, j, k, hx_nx, hx_ny)] - Hx[fdtd_idx3(i, j - 1, k, hx_nx, hx_ny)]) * Idyh[j];
   Ez[id] = g1[medio] * Ez[id] + g2[medio] * curl;
}

static __global__ void k_advance_hx2(
   fdtd_real *__restrict__ Hx, const fdtd_real *__restrict__ Ey, const fdtd_real *__restrict__ Ez,
   const fdtd_media *__restrict__ Mi, const fdtd_real *__restrict__ gm1, const fdtd_real *__restrict__ gm2,
   const fdtd_real *__restrict__ Idze, const fdtd_real *__restrict__ Idye,
   int nx, int ny, int ey_nx, int ey_ny, int ez_nx, int ez_ny,
   int is, int ie, int js, int je, int ks, int ke)
{
   int i = is + blockIdx.x * blockDim.x + threadIdx.x;
   int j = js + blockIdx.y * blockDim.y + threadIdx.y;
   int k = ks + blockIdx.z * blockDim.z + threadIdx.z;
   if (i > ie || j > je || k > ke) return;
   size_t id = fdtd_idx3(i, j, k, nx, ny);
   fdtd_media medio = Mi[id];
   fdtd_real curl =
      (Ey[fdtd_idx3(i, j, k + 1, ey_nx, ey_ny)] - Ey[fdtd_idx3(i, j, k, ey_nx, ey_ny)]) * Idze[k] -
      (Ez[fdtd_idx3(i, j + 1, k, ez_nx, ez_ny)] - Ez[fdtd_idx3(i, j, k, ez_nx, ez_ny)]) * Idye[j];
   Hx[id] = gm1[medio] * Hx[id] + gm2[medio] * curl;
}

static __global__ void k_advance_hy2(
   fdtd_real *__restrict__ Hy, const fdtd_real *__restrict__ Ez, const fdtd_real *__restrict__ Ex,
   const fdtd_media *__restrict__ Mi, const fdtd_real *__restrict__ gm1, const fdtd_real *__restrict__ gm2,
   const fdtd_real *__restrict__ Idxe, const fdtd_real *__restrict__ Idze,
   int nx, int ny, int ez_nx, int ez_ny, int ex_nx, int ex_ny,
   int is, int ie, int js, int je, int ks, int ke)
{
   int i = is + blockIdx.x * blockDim.x + threadIdx.x;
   int j = js + blockIdx.y * blockDim.y + threadIdx.y;
   int k = ks + blockIdx.z * blockDim.z + threadIdx.z;
   if (i > ie || j > je || k > ke) return;
   size_t id = fdtd_idx3(i, j, k, nx, ny);
   fdtd_media medio = Mi[id];
   fdtd_real curl =
      (Ez[fdtd_idx3(i + 1, j, k, ez_nx, ez_ny)] - Ez[fdtd_idx3(i, j, k, ez_nx, ez_ny)]) * Idxe[i] -
      (Ex[fdtd_idx3(i, j, k + 1, ex_nx, ex_ny)] - Ex[fdtd_idx3(i, j, k, ex_nx, ex_ny)]) * Idze[k];
   Hy[id] = gm1[medio] * Hy[id] + gm2[medio] * curl;
}

static __global__ void k_advance_hz2(
   fdtd_real *__restrict__ Hz, const fdtd_real *__restrict__ Ex, const fdtd_real *__restrict__ Ey,
   const fdtd_media *__restrict__ Mi, const fdtd_real *__restrict__ gm1, const fdtd_real *__restrict__ gm2,
   const fdtd_real *__restrict__ Idye, const fdtd_real *__restrict__ Idxe,
   int nx, int ny, int ex_nx, int ex_ny, int ey_nx, int ey_ny,
   int is, int ie, int js, int je, int ks, int ke)
{
   int i = is + blockIdx.x * blockDim.x + threadIdx.x;
   int j = js + blockIdx.y * blockDim.y + threadIdx.y;
   int k = ks + blockIdx.z * blockDim.z + threadIdx.z;
   if (i > ie || j > je || k > ke) return;
   size_t id = fdtd_idx3(i, j, k, nx, ny);
   fdtd_media medio = Mi[id];
   fdtd_real curl =
      (Ex[fdtd_idx3(i, j + 1, k, ex_nx, ex_ny)] - Ex[fdtd_idx3(i, j, k, ex_nx, ex_ny)]) * Idye[j] -
      (Ey[fdtd_idx3(i + 1, j, k, ey_nx, ey_ny)] - Ey[fdtd_idx3(i, j, k, ey_nx, ey_ny)]) * Idxe[i];
   Hz[id] = gm1[medio] * Hz[id] + gm2[medio] * curl;
}

static dim3 grid_for(fdtd_ibox s, dim3 block)
{
   int nx = s.ie - s.is + 1;
   int ny = s.je - s.js + 1;
   int nz = s.ke - s.ks + 1;
   if (nx < 1) nx = 1;
   if (ny < 1) ny = 1;
   if (nz < 1) nz = 1;
   return dim3((nx + block.x - 1) / block.x,
               (ny + block.y - 1) / block.y,
               (nz + block.z - 1) / block.z);
}

static int launch_ok(void)
{
   return cudaPeekAtLastError() == cudaSuccess && cudaDeviceSynchronize() == cudaSuccess;
}

int fdtd_cuda_advance_ex(fdtd_cuda_ctx *ctx, fdtd_ibox s)
{
   if (!ctx) return 0;
   fdtd_dims3 ex = fdtd_cuda_dim_Ex(ctx), hy = fdtd_cuda_dim_Hy(ctx), hz = fdtd_cuda_dim_Hz(ctx);
   dim3 block(8, 8, 4);
   k_advance_ex2<<<grid_for(s, block), block>>>(
      fdtd_cuda_ptr_Ex(ctx), fdtd_cuda_ptr_Hy(ctx), fdtd_cuda_ptr_Hz(ctx),
      fdtd_cuda_ptr_mEx(ctx), fdtd_cuda_ptr_g1(ctx), fdtd_cuda_ptr_g2(ctx),
      fdtd_cuda_ptr_Idyh(ctx), fdtd_cuda_ptr_Idzh(ctx),
      ex.nx, ex.ny, hy.nx, hy.ny, hz.nx, hz.ny,
      s.is, s.ie, s.js, s.je, s.ks, s.ke);
   return launch_ok();
}

int fdtd_cuda_advance_ey(fdtd_cuda_ctx *ctx, fdtd_ibox s)
{
   if (!ctx) return 0;
   fdtd_dims3 ey = fdtd_cuda_dim_Ey(ctx), hz = fdtd_cuda_dim_Hz(ctx), hx = fdtd_cuda_dim_Hx(ctx);
   dim3 block(8, 8, 4);
   k_advance_ey2<<<grid_for(s, block), block>>>(
      fdtd_cuda_ptr_Ey(ctx), fdtd_cuda_ptr_Hz(ctx), fdtd_cuda_ptr_Hx(ctx),
      fdtd_cuda_ptr_mEy(ctx), fdtd_cuda_ptr_g1(ctx), fdtd_cuda_ptr_g2(ctx),
      fdtd_cuda_ptr_Idzh(ctx), fdtd_cuda_ptr_Idxh(ctx),
      ey.nx, ey.ny, hz.nx, hz.ny, hx.nx, hx.ny,
      s.is, s.ie, s.js, s.je, s.ks, s.ke);
   return launch_ok();
}

int fdtd_cuda_advance_ez(fdtd_cuda_ctx *ctx, fdtd_ibox s)
{
   if (!ctx) return 0;
   fdtd_dims3 ez = fdtd_cuda_dim_Ez(ctx), hx = fdtd_cuda_dim_Hx(ctx), hy = fdtd_cuda_dim_Hy(ctx);
   dim3 block(8, 8, 4);
   k_advance_ez2<<<grid_for(s, block), block>>>(
      fdtd_cuda_ptr_Ez(ctx), fdtd_cuda_ptr_Hx(ctx), fdtd_cuda_ptr_Hy(ctx),
      fdtd_cuda_ptr_mEz(ctx), fdtd_cuda_ptr_g1(ctx), fdtd_cuda_ptr_g2(ctx),
      fdtd_cuda_ptr_Idyh(ctx), fdtd_cuda_ptr_Idxh(ctx),
      ez.nx, ez.ny, hx.nx, hx.ny, hy.nx, hy.ny,
      s.is, s.ie, s.js, s.je, s.ks, s.ke);
   return launch_ok();
}

int fdtd_cuda_advance_hx(fdtd_cuda_ctx *ctx, fdtd_ibox s)
{
   if (!ctx) return 0;
   fdtd_dims3 hx = fdtd_cuda_dim_Hx(ctx), ey = fdtd_cuda_dim_Ey(ctx), ez = fdtd_cuda_dim_Ez(ctx);
   dim3 block(8, 8, 4);
   k_advance_hx2<<<grid_for(s, block), block>>>(
      fdtd_cuda_ptr_Hx(ctx), fdtd_cuda_ptr_Ey(ctx), fdtd_cuda_ptr_Ez(ctx),
      fdtd_cuda_ptr_mHx(ctx), fdtd_cuda_ptr_gm1(ctx), fdtd_cuda_ptr_gm2(ctx),
      fdtd_cuda_ptr_Idze(ctx), fdtd_cuda_ptr_Idye(ctx),
      hx.nx, hx.ny, ey.nx, ey.ny, ez.nx, ez.ny,
      s.is, s.ie, s.js, s.je, s.ks, s.ke);
   return launch_ok();
}

int fdtd_cuda_advance_hy(fdtd_cuda_ctx *ctx, fdtd_ibox s)
{
   if (!ctx) return 0;
   fdtd_dims3 hy = fdtd_cuda_dim_Hy(ctx), ez = fdtd_cuda_dim_Ez(ctx), ex = fdtd_cuda_dim_Ex(ctx);
   dim3 block(8, 8, 4);
   k_advance_hy2<<<grid_for(s, block), block>>>(
      fdtd_cuda_ptr_Hy(ctx), fdtd_cuda_ptr_Ez(ctx), fdtd_cuda_ptr_Ex(ctx),
      fdtd_cuda_ptr_mHy(ctx), fdtd_cuda_ptr_gm1(ctx), fdtd_cuda_ptr_gm2(ctx),
      fdtd_cuda_ptr_Idxe(ctx), fdtd_cuda_ptr_Idze(ctx),
      hy.nx, hy.ny, ez.nx, ez.ny, ex.nx, ex.ny,
      s.is, s.ie, s.js, s.je, s.ks, s.ke);
   return launch_ok();
}

int fdtd_cuda_advance_hz(fdtd_cuda_ctx *ctx, fdtd_ibox s)
{
   if (!ctx) return 0;
   fdtd_dims3 hz = fdtd_cuda_dim_Hz(ctx), ex = fdtd_cuda_dim_Ex(ctx), ey = fdtd_cuda_dim_Ey(ctx);
   dim3 block(8, 8, 4);
   k_advance_hz2<<<grid_for(s, block), block>>>(
      fdtd_cuda_ptr_Hz(ctx), fdtd_cuda_ptr_Ex(ctx), fdtd_cuda_ptr_Ey(ctx),
      fdtd_cuda_ptr_mHz(ctx), fdtd_cuda_ptr_gm1(ctx), fdtd_cuda_ptr_gm2(ctx),
      fdtd_cuda_ptr_Idye(ctx), fdtd_cuda_ptr_Idxe(ctx),
      hz.nx, hz.ny, ex.nx, ex.ny, ey.nx, ey.ny,
      s.is, s.ie, s.js, s.je, s.ks, s.ke);
   return launch_ok();
}
