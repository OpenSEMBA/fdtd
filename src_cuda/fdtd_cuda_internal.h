#pragma once
#include "fdtd_cuda.h"

#define FDTD_CUDA_PSI_SLOTS 64
#define FDTD_CUDA_CPML1D    24

struct fdtd_cuda_ctx {
   int ok;
   fdtd_dims3 dex, dey, dez, dhx, dhy, dhz;
   fdtd_dims3 mex, mey, mez, mhx, mhy, mhz;
   fdtd_real *d_Ex, *d_Ey, *d_Ez, *d_Hx, *d_Hy, *d_Hz;
   fdtd_media *d_mEx, *d_mEy, *d_mEz, *d_mHx, *d_mHy, *d_mHz;
   fdtd_real *d_g1, *d_g2, *d_gm1, *d_gm2;
   int num_media;
   fdtd_real *d_Idxh, *d_Idyh, *d_Idzh, *d_Idxe, *d_Idye, *d_Idze;
   int n_idxh, n_idyh, n_idzh, n_idxe, n_idye, n_idze;
   fdtd_real *d_psi[FDTD_CUDA_PSI_SLOTS];
   int psi_n[FDTD_CUDA_PSI_SLOTS];
   fdtd_real *d_cpml1d[FDTD_CUDA_CPML1D];
   int cpml1d_n[FDTD_CUDA_CPML1D];
   /* Sparse point-probe gather (device SoA + output scratch). */
   int n_probes;
   int *d_probe_comp;
   int *d_probe_i;
   int *d_probe_j;
   int *d_probe_k;
   fdtd_real *d_probe_out;

   /* Device Huygens planewave (face kernels + Incid tables). */
   int pw_ready;
   int pw_n_waves, pw_max_modes, pw_max_numus;
   int pw_n_faces_e, pw_n_faces_h;
   fdtd_real pw_cluz;
   fdtd_real *d_pw_evol;
   fdtd_real *d_pw_deltaevol;
   int *d_pw_numus;
   int *d_pw_num_modes;
   fdtd_real *d_pw_px, *d_pw_py, *d_pw_pz, *d_pw_d0;
   fdtd_real *d_pw_fpw;
   fdtd_pw_face *d_pw_faces_e, *d_pw_faces_h;
   fdtd_pw_face *h_pw_faces_e, *h_pw_faces_h; /* host copies for launch */
   int *d_pw_still;
   fdtd_real *d_pw_phys_x[6], *d_pw_phys_y[6], *d_pw_phys_z[6];
   int pw_phys_x0[6], pw_phys_y0[6], pw_phys_z0[6];
   int pw_phys_nx[6], pw_phys_ny[6], pw_phys_nz[6];

   /* Device first-order magnetic Mur ABC. */
   int mur_ready;
   int mur_n_jobs;
   int mur_cab_n[FDTD_CUDA_MUR_CAB];
   fdtd_real *d_mur_cab[FDTD_CUDA_MUR_CAB];
   fdtd_mur_job h_mur_jobs[FDTD_CUDA_MUR_JOBS];
   fdtd_real *d_mur_past[FDTD_CUDA_MUR_JOBS];
   int mur_past_n[FDTD_CUDA_MUR_JOBS];
};

void fdtd_cuda_free_planewave(fdtd_cuda_ctx *ctx);
void fdtd_cuda_free_mur(fdtd_cuda_ctx *ctx);

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
