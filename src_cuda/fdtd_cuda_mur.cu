#include "fdtd_cuda_internal.h"

#include <cuda_runtime.h>
#include <stdio.h>
#include <string.h>

static int check_mur(cudaError_t e, const char *what)
{
   if (e != cudaSuccess) {
      fprintf(stderr, "FDTD CUDA Mur error in %s: %s\n", what, cudaGetErrorString(e));
      return 0;
   }
   return 1;
}

static void free_mur(fdtd_cuda_ctx *ctx)
{
#define FREE(p) do { if (p) cudaFree(p); p = NULL; } while (0)
   for (int i = 0; i < FDTD_CUDA_MUR_CAB; ++i) {
      FREE(ctx->d_mur_cab[i]);
      ctx->mur_cab_n[i] = 0;
   }
   for (int i = 0; i < FDTD_CUDA_MUR_JOBS; ++i) {
      FREE(ctx->d_mur_past[i]);
      ctx->mur_past_n[i] = 0;
   }
#undef FREE
   ctx->mur_n_jobs = 0;
   ctx->mur_ready = 0;
   memset(ctx->h_mur_jobs, 0, sizeof(ctx->h_mur_jobs));
}

void fdtd_cuda_free_mur(fdtd_cuda_ctx *ctx)
{
   if (!ctx) return;
   free_mur(ctx);
}

/*
 * First-order Mur on one ghost plane (host AdvanceMagneticMUR, mur_second=.false.):
 *   H(ghost) = Past(interior) + CAB1(medio) * (H(interior) - Past(ghost))
 * medio is taken at the interior neighbour (same as host sggMiH*(i_m, j_m±1, k_m)).
 *
 * Free axes: wall_axis 0 → a=j, b=k; 1 → a=i, b=k; 2 → a=i, b=j.
 */
static __global__ void k_mur_apply(
   fdtd_real *__restrict__ H,
   const fdtd_real *__restrict__ Past,
   const fdtd_media *__restrict__ Mi,
   const fdtd_real *__restrict__ cab,
   int n_cab,
   int f_nx, int f_ny, int f_nz, int m_nx, int m_ny, int nx_p, int ny_p, int nz_p,
   int wall_axis, int neigh_sign, int plane_abs,
   int a0, int a1, int b0, int b1,
   int e_xi, int e_yi, int e_zi,
   int p_xi, int p_yi, int p_zi)
{
   int a = a0 + blockIdx.x * blockDim.x + threadIdx.x;
   int b = b0 + blockIdx.y * blockDim.y + threadIdx.y;
   if (a > a1 || b > b1) return;

   int i, j, k;
   if (wall_axis == 0) {
      i = plane_abs;
      j = a;
      k = b;
   } else if (wall_axis == 1) {
      i = a;
      j = plane_abs;
      k = b;
   } else {
      i = a;
      j = b;
      k = plane_abs;
   }

   int i_m = i - e_xi;
   int j_m = j - e_yi;
   int k_m = k - e_zi;
   int i_n = i_m, j_n = j_m, k_n = k_m;
   int i_pn = i, j_pn = j, k_pn = k;
   if (wall_axis == 0) {
      i_n += neigh_sign;
      i_pn += neigh_sign;
   } else if (wall_axis == 1) {
      j_n += neigh_sign;
      j_pn += neigh_sign;
   } else {
      k_n += neigh_sign;
      k_pn += neigh_sign;
   }

   int pi_g = i - p_xi, pj_g = j - p_yi, pk_g = k - p_zi;
   int pi_n = i_pn - p_xi, pj_n = j_pn - p_yi, pk_n = k_pn - p_zi;
   if (i_m < 0 || j_m < 0 || k_m < 0 || i_n < 0 || j_n < 0 || k_n < 0) return;
   if (i_m >= f_nx || j_m >= f_ny || k_m >= f_nz) return;
   if (i_n >= f_nx || j_n >= f_ny || k_n >= f_nz) return;
   if (pi_g < 0 || pj_g < 0 || pk_g < 0 || pi_n < 0 || pj_n < 0 || pk_n < 0) return;
   if (pi_g >= nx_p || pj_g >= ny_p || pk_g >= nz_p) return;
   if (pi_n >= nx_p || pj_n >= ny_p || pk_n >= nz_p) return;

   size_t id_g = fdtd_idx3(i_m, j_m, k_m, f_nx, f_ny);
   size_t id_n = fdtd_idx3(i_n, j_n, k_n, f_nx, f_ny);
   size_t id_pg = fdtd_idx3(pi_g, pj_g, pk_g, nx_p, ny_p);
   size_t id_pn = fdtd_idx3(pi_n, pj_n, pk_n, nx_p, ny_p);
   size_t id_m = fdtd_idx3(i_n, j_n, k_n, m_nx, m_ny);

   int medio = (int)Mi[id_m];
   if (medio < 0 || medio >= n_cab) medio = 0;
   H[id_g] = Past[id_pn] + cab[medio] * (H[id_n] - Past[id_pg]);
}

static __global__ void k_mur_store(
   fdtd_real *__restrict__ Past,
   const fdtd_real *__restrict__ H,
   int f_nx, int f_ny, int f_nz, int nx_p, int ny_p, int nz_p,
   int xi, int xe, int yi, int ye, int zi, int ze,
   int e_xi, int e_yi, int e_zi,
   int p_xi, int p_yi, int p_zi)
{
   int i = xi + blockIdx.x * blockDim.x + threadIdx.x;
   int j = yi + blockIdx.y * blockDim.y + threadIdx.y;
   int k = zi + blockIdx.z * blockDim.z + threadIdx.z;
   if (i > xe || j > ye || k > ze) return;
   int i_m = i - e_xi, j_m = j - e_yi, k_m = k - e_zi;
   int pi = i - p_xi, pj = j - p_yi, pk = k - p_zi;
   if (i_m < 0 || j_m < 0 || k_m < 0 || pi < 0 || pj < 0 || pk < 0) return;
   if (i_m >= f_nx || j_m >= f_ny || k_m >= f_nz) return;
   if (pi >= nx_p || pj >= ny_p || pk >= nz_p) return;
   size_t id_p = fdtd_idx3(pi, pj, pk, nx_p, ny_p);
   size_t id_h = fdtd_idx3(i_m, j_m, k_m, f_nx, f_ny);
   Past[id_p] = H[id_h];
}

int fdtd_cuda_upload_mur_cab(fdtd_cuda_ctx *ctx, int which, int n, const fdtd_real *host)
{
   if (!ctx || which < 0 || which >= FDTD_CUDA_MUR_CAB || n <= 0 || !host) return 0;
   if (ctx->d_mur_cab[which]) cudaFree(ctx->d_mur_cab[which]);
   ctx->mur_cab_n[which] = n;
   return check_mur(cudaMalloc((void **)&ctx->d_mur_cab[which], (size_t)n * sizeof(fdtd_real)),
                    "mur cab") &&
          check_mur(cudaMemcpy(ctx->d_mur_cab[which], host, (size_t)n * sizeof(fdtd_real),
                               cudaMemcpyHostToDevice),
                    "mur cab H2D");
}

int fdtd_cuda_set_mur_jobs(fdtd_cuda_ctx *ctx, const fdtd_mur_job *jobs, int n_jobs)
{
   if (!ctx) return 0;
   for (int i = 0; i < FDTD_CUDA_MUR_JOBS; ++i) {
      if (ctx->d_mur_past[i]) {
         cudaFree(ctx->d_mur_past[i]);
         ctx->d_mur_past[i] = NULL;
      }
      ctx->mur_past_n[i] = 0;
   }
   ctx->mur_n_jobs = 0;
   ctx->mur_ready = 0;
   if (n_jobs < 0 || n_jobs > FDTD_CUDA_MUR_JOBS) return 0;
   if (n_jobs == 0) return 1;
   if (!jobs) return 0;
   memcpy(ctx->h_mur_jobs, jobs, (size_t)n_jobs * sizeof(fdtd_mur_job));
   for (int s = 0; s < n_jobs; ++s) {
      const fdtd_mur_job *j = &ctx->h_mur_jobs[s];
      int n = j->nx_p * j->ny_p * j->nz_p;
      if (n <= 0) return 0;
      if (!check_mur(cudaMalloc((void **)&ctx->d_mur_past[s], (size_t)n * sizeof(fdtd_real)),
                     "mur past") ||
          !check_mur(cudaMemset(ctx->d_mur_past[s], 0, (size_t)n * sizeof(fdtd_real)),
                     "mur past 0")) {
         return 0;
      }
      ctx->mur_past_n[s] = n;
   }
   ctx->mur_n_jobs = n_jobs;
   ctx->mur_ready = 1;
   return 1;
}

int fdtd_cuda_upload_mur_past(fdtd_cuda_ctx *ctx, int slot, const fdtd_real *host, int n_elem)
{
   if (!ctx || slot < 0 || slot >= ctx->mur_n_jobs || !host) return 0;
   if (!ctx->d_mur_past[slot] || n_elem != ctx->mur_past_n[slot]) return 0;
   return check_mur(cudaMemcpy(ctx->d_mur_past[slot], host, (size_t)n_elem * sizeof(fdtd_real),
                               cudaMemcpyHostToDevice),
                    "mur past H2D");
}

int fdtd_cuda_download_mur_past(fdtd_cuda_ctx *ctx, int slot, fdtd_real *host, int n_elem)
{
   if (!ctx || slot < 0 || slot >= ctx->mur_n_jobs || !host) return 0;
   if (!ctx->d_mur_past[slot] || n_elem != ctx->mur_past_n[slot]) return 0;
   return check_mur(cudaMemcpy(host, ctx->d_mur_past[slot], (size_t)n_elem * sizeof(fdtd_real),
                               cudaMemcpyDeviceToHost),
                    "mur past D2H");
}

int fdtd_cuda_mur_ready(const fdtd_cuda_ctx *ctx)
{
   return ctx && ctx->mur_ready && ctx->mur_n_jobs > 0;
}

static int launch_apply(fdtd_cuda_ctx *ctx, int slot)
{
   const fdtd_mur_job *j = &ctx->h_mur_jobs[slot];
   fdtd_real *H = fdtd_cuda_field_ptr(ctx, j->field_comp);
   fdtd_media *Mi = fdtd_cuda_media_ptr(ctx, j->media_comp);
   fdtd_dims3 fd = fdtd_cuda_field_dim(ctx, j->field_comp);
   fdtd_dims3 md = fdtd_cuda_media_dim(ctx, j->media_comp);
   fdtd_real *cab = ctx->d_mur_cab[j->cab_which];
   int n_cab = ctx->mur_cab_n[j->cab_which];
   if (!H || !Mi || !cab || !ctx->d_mur_past[slot] || n_cab <= 0) return 0;

   int na = j->a1 - j->a0 + 1;
   int nb = j->b1 - j->b0 + 1;
   if (na < 1 || nb < 1) return 1;

   dim3 block(16, 16);
   dim3 grid((na + block.x - 1) / block.x, (nb + block.y - 1) / block.y);
   k_mur_apply<<<grid, block>>>(
      H, ctx->d_mur_past[slot], Mi, cab, n_cab,
      fd.nx, fd.ny, fd.nz, md.nx, md.ny, j->nx_p, j->ny_p, j->nz_p,
      j->wall_axis, j->neigh_sign, j->plane_abs,
      j->a0, j->a1, j->b0, j->b1,
      j->e_xi, j->e_yi, j->e_zi,
      j->p_xi, j->p_yi, j->p_zi);
   return check_mur(cudaGetLastError(), "k_mur_apply");
}

static int launch_store(fdtd_cuda_ctx *ctx, int slot)
{
   const fdtd_mur_job *j = &ctx->h_mur_jobs[slot];
   fdtd_real *H = fdtd_cuda_field_ptr(ctx, j->field_comp);
   fdtd_dims3 fd = fdtd_cuda_field_dim(ctx, j->field_comp);
   if (!H || !ctx->d_mur_past[slot]) return 0;

   int nx = j->store_xe - j->store_xi + 1;
   int ny = j->store_ye - j->store_yi + 1;
   int nz = j->store_ze - j->store_zi + 1;
   if (nx < 1 || ny < 1 || nz < 1) return 1;

   dim3 block(8, 8, 4);
   dim3 grid((nx + block.x - 1) / block.x,
             (ny + block.y - 1) / block.y,
             (nz + block.z - 1) / block.z);
   k_mur_store<<<grid, block>>>(
      ctx->d_mur_past[slot], H,
      fd.nx, fd.ny, fd.nz, j->nx_p, j->ny_p, j->nz_p,
      j->store_xi, j->store_xe, j->store_yi, j->store_ye, j->store_zi, j->store_ze,
      j->e_xi, j->e_yi, j->e_zi,
      j->p_xi, j->p_yi, j->p_zi);
   return check_mur(cudaGetLastError(), "k_mur_store");
}

int fdtd_cuda_advance_mur(fdtd_cuda_ctx *ctx)
{
   if (!ctx || !ctx->mur_ready) return 0;
   for (int s = 0; s < ctx->mur_n_jobs; ++s) {
      if (!launch_apply(ctx, s)) return 0;
   }
   for (int s = 0; s < ctx->mur_n_jobs; ++s) {
      if (!launch_store(ctx, s)) return 0;
   }
   return check_mur(cudaDeviceSynchronize(), "mur sync");
}
