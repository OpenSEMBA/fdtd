#include "fdtd_cuda_internal.h"

#include <cuda_runtime.h>

/*
 * CPML face update matching AdvanceelectricCPML / AdvanceMagneticCPML loops:
 *   Psi = P_b(abs_free) * Psi + (H_a - H_b) * P_c(abs_free)
 *   Field += h_sign * G2(medio) * Psi
 *
 * Absolute (i,j,k) in [xi..xe]x[yi..ye]x[zi..ze].
 * Field index: (i - e_xi, j - e_yi, k - e_zi) in field array.
 * Psi index:   (i - psi_xi, ...) in psi array.
 * Neighbour H_b is H_a shifted -1 along h_diff_axis.
 */
static __global__ void k_cpml(
   fdtd_real *__restrict__ Field,
   fdtd_real *__restrict__ Psi,
   const fdtd_real *__restrict__ Ha,
   const fdtd_real *__restrict__ Hb_plane, /* same buffer as Ha; offset via axis */
   const fdtd_media *__restrict__ Mi,
   const fdtd_real *__restrict__ G,
   const fdtd_real *__restrict__ P_b,
   const fdtd_real *__restrict__ P_c,
   int f_nx, int f_ny, int h_nx, int h_ny, int m_nx, int m_ny,
   int nx_psi, int ny_psi,
   int xi, int xe, int yi, int ye, int zi, int ze,
   int e_xi, int e_yi, int e_zi,
   int psi_xi, int psi_yi, int psi_zi,
   int h_diff_axis, int free_axis, int p_base,
   int h_sign, int use_fixed_medio, int medio_fixed)
{
   int i = xi + blockIdx.x * blockDim.x + threadIdx.x;
   int j = yi + blockIdx.y * blockDim.y + threadIdx.y;
   int k = zi + blockIdx.z * blockDim.z + threadIdx.z;
   if (i > xe || j > ye || k > ze) return;

   int i_m = i - e_xi;
   int j_m = j - e_yi;
   int k_m = k - e_zi;

   int i_h0 = i_m, j_h0 = j_m, k_h0 = k_m;
   int i_h1 = i_m, j_h1 = j_m, k_h1 = k_m;
   if (h_diff_axis == 0) i_h1 = i_m - 1;
   else if (h_diff_axis == 1) j_h1 = j_m - 1;
   else k_h1 = k_m - 1;

   size_t id_f = fdtd_idx3(i_m, j_m, k_m, f_nx, f_ny);
   size_t id_h0 = fdtd_idx3(i_h0, j_h0, k_h0, h_nx, h_ny);
   size_t id_h1 = fdtd_idx3(i_h1, j_h1, k_h1, h_nx, h_ny);
   size_t id_p = fdtd_idx3(i - psi_xi, j - psi_yi, k - psi_zi, nx_psi, ny_psi);
   size_t id_m = fdtd_idx3(i_m, j_m, k_m, m_nx, m_ny);

   int abs_free = (free_axis == 0) ? i : (free_axis == 1) ? j : k;
   int pidx = abs_free - p_base;

   fdtd_real dH = Ha[id_h0] - Hb_plane[id_h1];
   Psi[id_p] = P_b[pidx] * Psi[id_p] + dH * P_c[pidx];

   int medio = use_fixed_medio ? medio_fixed : (int)Mi[id_m];
   Field[id_f] = Field[id_f] + (fdtd_real)h_sign * G[medio] * Psi[id_p];
}

static int launch_cpml(fdtd_cuda_ctx *ctx, const fdtd_cpml_job *job)
{
   fdtd_real *Field = fdtd_cuda_field_ptr(ctx, job->field_comp);
   fdtd_real *Ha = fdtd_cuda_field_ptr(ctx, job->h_comp_a);
   /* h_comp_b selects which field buffer; neighbour offset is via axis */
   (void)job->h_comp_b;
   fdtd_dims3 fd = fdtd_cuda_field_dim(ctx, job->field_comp);
   fdtd_dims3 hd = fdtd_cuda_field_dim(ctx, job->h_comp_a);
   fdtd_dims3 md = fdtd_cuda_media_dim(ctx, job->media_comp);
   fdtd_media *Mi = fdtd_cuda_media_ptr(ctx, job->media_comp);
   fdtd_real *G = (job->field_comp < 3) ? fdtd_cuda_ptr_g2(ctx) : fdtd_cuda_ptr_gm2(ctx);
   fdtd_real *Psi = fdtd_cuda_ptr_psi(ctx, job->psi_slot);
   fdtd_real *Pb = fdtd_cuda_ptr_cpml1d(ctx, job->p_b_which);
   fdtd_real *Pc = fdtd_cuda_ptr_cpml1d(ctx, job->p_c_which);
   if (!Field || !Ha || !Psi || !Pb || !Pc || !Mi || !G) return 0;

   int nx = job->xe - job->xi + 1;
   int ny = job->ye - job->yi + 1;
   int nz = job->ze - job->zi + 1;
   if (nx < 1 || ny < 1 || nz < 1) return 1; /* empty region OK */

   dim3 block(8, 8, 4);
   dim3 grid((nx + block.x - 1) / block.x,
             (ny + block.y - 1) / block.y,
             (nz + block.z - 1) / block.z);

   k_cpml<<<grid, block>>>(
      Field, Psi, Ha, Ha, Mi, G, Pb, Pc,
      fd.nx, fd.ny, hd.nx, hd.ny, md.nx, md.ny,
      job->nx_psi, job->ny_psi,
      job->xi, job->xe, job->yi, job->ye, job->zi, job->ze,
      job->e_xi, job->e_yi, job->e_zi,
      job->psi_xi, job->psi_yi, job->psi_zi,
      job->h_diff_axis, job->free_axis, job->p_base,
      job->h_sign, job->use_fixed_medio, job->medio_fixed);

   return cudaPeekAtLastError() == cudaSuccess;
}

int fdtd_cuda_cpml_apply(fdtd_cuda_ctx *ctx, const fdtd_cpml_job *job)
{
   if (!ctx || !job || !launch_cpml(ctx, job)) return 0;
   return cudaDeviceSynchronize() == cudaSuccess;
}

int fdtd_cuda_cpml_apply_n(fdtd_cuda_ctx *ctx, const fdtd_cpml_job *jobs, int n)
{
   if (!ctx || n < 0) return 0;
   if (n == 0) return 1;
   if (!jobs) return 0;
   for (int i = 0; i < n; ++i) {
      if (!launch_cpml(ctx, &jobs[i])) return 0;
   }
   return cudaDeviceSynchronize() == cudaSuccess;
}
