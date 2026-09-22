#pragma once
#include "fdtd_cuda.h"

#ifdef __cplusplus
extern "C" {
#endif

fdtd_real *fdtd_cuda_ptr_Ex(fdtd_cuda_ctx *c);
fdtd_real *fdtd_cuda_ptr_Ey(fdtd_cuda_ctx *c);
fdtd_real *fdtd_cuda_ptr_Ez(fdtd_cuda_ctx *c);
fdtd_real *fdtd_cuda_ptr_Hx(fdtd_cuda_ctx *c);
fdtd_real *fdtd_cuda_ptr_Hy(fdtd_cuda_ctx *c);
fdtd_real *fdtd_cuda_ptr_Hz(fdtd_cuda_ctx *c);
fdtd_media *fdtd_cuda_ptr_mEx(fdtd_cuda_ctx *c);
fdtd_media *fdtd_cuda_ptr_mEy(fdtd_cuda_ctx *c);
fdtd_media *fdtd_cuda_ptr_mEz(fdtd_cuda_ctx *c);
fdtd_media *fdtd_cuda_ptr_mHx(fdtd_cuda_ctx *c);
fdtd_media *fdtd_cuda_ptr_mHy(fdtd_cuda_ctx *c);
fdtd_media *fdtd_cuda_ptr_mHz(fdtd_cuda_ctx *c);
fdtd_real *fdtd_cuda_ptr_g1(fdtd_cuda_ctx *c);
fdtd_real *fdtd_cuda_ptr_g2(fdtd_cuda_ctx *c);
fdtd_real *fdtd_cuda_ptr_gm1(fdtd_cuda_ctx *c);
fdtd_real *fdtd_cuda_ptr_gm2(fdtd_cuda_ctx *c);
fdtd_real *fdtd_cuda_ptr_Idxh(fdtd_cuda_ctx *c);
fdtd_real *fdtd_cuda_ptr_Idyh(fdtd_cuda_ctx *c);
fdtd_real *fdtd_cuda_ptr_Idzh(fdtd_cuda_ctx *c);
fdtd_real *fdtd_cuda_ptr_Idxe(fdtd_cuda_ctx *c);
fdtd_real *fdtd_cuda_ptr_Idye(fdtd_cuda_ctx *c);
fdtd_real *fdtd_cuda_ptr_Idze(fdtd_cuda_ctx *c);
fdtd_dims3 fdtd_cuda_dim_Ex(fdtd_cuda_ctx *c);
fdtd_dims3 fdtd_cuda_dim_Ey(fdtd_cuda_ctx *c);
fdtd_dims3 fdtd_cuda_dim_Ez(fdtd_cuda_ctx *c);
fdtd_dims3 fdtd_cuda_dim_Hx(fdtd_cuda_ctx *c);
fdtd_dims3 fdtd_cuda_dim_Hy(fdtd_cuda_ctx *c);
fdtd_dims3 fdtd_cuda_dim_Hz(fdtd_cuda_ctx *c);
fdtd_real *fdtd_cuda_ptr_psi(fdtd_cuda_ctx *c, int s);
fdtd_real *fdtd_cuda_ptr_cpml1d(fdtd_cuda_ctx *c, int w);
fdtd_real *fdtd_cuda_field_ptr(fdtd_cuda_ctx *c, int comp);
fdtd_dims3 fdtd_cuda_field_dim(fdtd_cuda_ctx *c, int comp);
fdtd_media *fdtd_cuda_media_ptr(fdtd_cuda_ctx *c, int comp);
fdtd_dims3 fdtd_cuda_media_dim(fdtd_cuda_ctx *c, int comp);

#ifdef __cplusplus
}
#endif

/* Fortran column-major 0-based linear index */
__host__ __device__ inline size_t fdtd_idx3(int i, int j, int k, int nx, int ny)
{
   return (size_t)i + (size_t)j * (size_t)nx + (size_t)k * (size_t)nx * (size_t)ny;
}
