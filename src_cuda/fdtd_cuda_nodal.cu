#include "fdtd_cuda_internal.h"

#include <cuda_runtime.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int check_nodal(cudaError_t e, const char *what)
{
   if (e != cudaSuccess) {
      fprintf(stderr, "FDTD CUDA nodal error in %s: %s\n", what, cudaGetErrorString(e));
      return 0;
   }
   return 1;
}

static void free_nodal(fdtd_cuda_ctx *ctx)
{
   if (ctx->h_nodal_jobs) {
      free(ctx->h_nodal_jobs);
      ctx->h_nodal_jobs = NULL;
   }
   if (ctx->d_nodal_samples) {
      cudaFree(ctx->d_nodal_samples);
      ctx->d_nodal_samples = NULL;
   }
   if (ctx->d_nodal_skip_e) {
      cudaFree(ctx->d_nodal_skip_e);
      ctx->d_nodal_skip_e = NULL;
   }
   if (ctx->d_nodal_skip_h) {
      cudaFree(ctx->d_nodal_skip_h);
      ctx->d_nodal_skip_h = NULL;
   }
   ctx->nodal_n_jobs = 0;
   ctx->nodal_n_skip = 0;
   ctx->nodal_ready = 0;
}

void fdtd_cuda_free_nodal(fdtd_cuda_ctx *ctx)
{
   if (!ctx) return;
   free_nodal(ctx);
}

/* Host evolucion(): samples are evol(0:numus). Out of range returns 0. */
__device__ static fdtd_real nodal_wave(const fdtd_nodal_job &job, const fdtd_real *samples, fdtd_real t)
{
   if (job.initial_only) return samples[job.evol_off];
   fdtd_real delta = job.deltaevol;
   int nprev = (int)(t / delta);
   if ((nprev + 1 > job.numus) || (nprev + 1 <= 0)) return (fdtd_real)0;
   fdtd_real e0 = samples[job.evol_off + nprev];
   fdtd_real e1 = samples[job.evol_off + nprev + 1];
   return (e1 - e0) / delta * (t - (fdtd_real)nprev * delta) + e0;
}

__device__ static int axis_index(int is_i, int is_j, int i_m, int j_m, int k_m)
{
   if (is_i) return i_m;
   if (is_j) return j_m;
   return k_m;
}

/*
 * One box. Hard assigns amp*wave. Soft subtracts coeff(medio)*ax*ay*amp*wave.
 * ax/ay are the host metric pair for this component (null on the hard path).
 * E skips PEC (skip_e), H skips PMC (skip_h).
 */
__global__ void k_nodal(
   fdtd_real *__restrict__ F,
   const fdtd_media *__restrict__ Mi,
   const fdtd_real *__restrict__ coeff,
   const int *__restrict__ skip,
   int n_skip,
   const fdtd_real *__restrict__ samples,
   const fdtd_real *__restrict__ ax,
   const fdtd_real *__restrict__ ay,
   int n_ax, int n_ay,
   int ax_i, int ax_j, int ay_i, int ay_j,
   fdtd_nodal_job job,
   int f_nx, int f_ny, int f_nz, int m_nx, int m_ny,
   fdtd_real time)
{
   int nxb = job.xe - job.xi + 1;
   int nyb = job.ye - job.yi + 1;
   int nzb = job.ze - job.zi + 1;
   if (nxb <= 0 || nyb <= 0 || nzb <= 0) return;
   int ncell = nxb * nyb * nzb;
   int tid = (int)(blockIdx.x * blockDim.x + threadIdx.x);
   if (tid >= ncell) return;
   int i = job.xi + (tid % nxb);
   int j = job.yi + ((tid / nxb) % nyb);
   int k = job.zi + (tid / (nxb * nyb));
   int i_m = i - job.e_xi;
   int j_m = j - job.e_yi;
   int k_m = k - job.e_zi;
   if (i_m < 0 || j_m < 0 || k_m < 0) return;
   if (i_m >= f_nx || j_m >= f_ny || k_m >= f_nz) return;

   size_t id = fdtd_idx3(i_m, j_m, k_m, f_nx, f_ny);
   int medio = (int)Mi[fdtd_idx3(i_m, j_m, k_m, m_nx, m_ny)];
   if (medio < 0 || medio >= n_skip) return;
   if (skip[medio]) return;

   fdtd_real wave = nodal_wave(job, samples, time);
   if (job.hard) {
      F[id] = job.amplitude * wave;
      return;
   }
   int ia = axis_index(ax_i, ax_j, i_m, j_m, k_m);
   int ib = axis_index(ay_i, ay_j, i_m, j_m, k_m);
   if (!ax || !ay || ia < 0 || ia >= n_ax || ib < 0 || ib >= n_ay) return;
   F[id] = F[id] - coeff[medio] * ax[ia] * ay[ib] * job.amplitude * wave;
}

static int launch_job(fdtd_cuda_ctx *ctx, const fdtd_nodal_job &job, fdtd_real time, int step)
{
   if (job.initial_only && step != 0) return 1;
   if (job.field_comp < 0 || job.field_comp > 5) return 0;

   fdtd_real *F = fdtd_cuda_field_ptr(ctx, job.field_comp);
   fdtd_media *Mi = fdtd_cuda_media_ptr(ctx, job.field_comp);
   fdtd_dims3 fd = fdtd_cuda_field_dim(ctx, job.field_comp);
   fdtd_dims3 md = fdtd_cuda_media_dim(ctx, job.field_comp);
   if (!F || !Mi) return 0;

   int electric = job.field_comp < 3;
   const fdtd_real *coeff = electric ? ctx->d_g2 : ctx->d_gm2;
   const int *skip = electric ? ctx->d_nodal_skip_e : ctx->d_nodal_skip_h;
   if (!job.hard && !coeff) return 0;
   if (!skip) return 0;

   const fdtd_real *ax = nullptr;
   const fdtd_real *ay = nullptr;
   int n_ax = 1, n_ay = 1;
   int ax_i = 0, ax_j = 0, ay_i = 0, ay_j = 0;
   if (!job.hard) {
      switch (job.field_comp) {
      case 0: /* Ex: Idyh(j) * Idzh(k) */
         ax = ctx->d_Idyh; n_ax = ctx->n_idyh; ax_j = 1;
         ay = ctx->d_Idzh; n_ay = ctx->n_idzh;
         break;
      case 1: /* Ey: Idxh(i) * Idzh(k) */
         ax = ctx->d_Idxh; n_ax = ctx->n_idxh; ax_i = 1;
         ay = ctx->d_Idzh; n_ay = ctx->n_idzh;
         break;
      case 2: /* Ez: Idxh(i) * Idyh(j) */
         ax = ctx->d_Idxh; n_ax = ctx->n_idxh; ax_i = 1;
         ay = ctx->d_Idyh; n_ay = ctx->n_idyh; ay_j = 1;
         break;
      case 3: /* Hx: Idye(j) * Idze(k) */
         ax = ctx->d_Idye; n_ax = ctx->n_idye; ax_j = 1;
         ay = ctx->d_Idze; n_ay = ctx->n_idze;
         break;
      case 4: /* Hy: Idxe(i) * Idze(k) */
         ax = ctx->d_Idxe; n_ax = ctx->n_idxe; ax_i = 1;
         ay = ctx->d_Idze; n_ay = ctx->n_idze;
         break;
      default: /* Hz: Idye(j) * Idxe(i) */
         ax = ctx->d_Idye; n_ax = ctx->n_idye; ax_j = 1;
         ay = ctx->d_Idxe; n_ay = ctx->n_idxe; ay_i = 1;
         break;
      }
      if (!ax || !ay) return 0;
   }

   int ncell = (job.xe - job.xi + 1) * (job.ye - job.yi + 1) * (job.ze - job.zi + 1);
   if (ncell <= 0) return 0;
   int threads = 128;
   int blocks = (ncell + threads - 1) / threads;
   k_nodal<<<blocks, threads>>>(
      F, Mi, coeff, skip, ctx->nodal_n_skip, ctx->d_nodal_samples,
      ax, ay, n_ax, n_ay, ax_i, ax_j, ay_i, ay_j,
      job, fd.nx, fd.ny, fd.nz, md.nx, md.ny, time);
   return check_nodal(cudaGetLastError(), "nodal launch");
}

int fdtd_cuda_upload_nodal(fdtd_cuda_ctx *ctx,
                           const fdtd_nodal_job *jobs, int n_jobs,
                           const fdtd_real *samples, int n_samples,
                           const int *skip_e, const int *skip_h, int n_skip)
{
   if (!ctx || n_skip <= 0 || !skip_e || !skip_h) return 0;
   free_nodal(ctx);
   if (n_jobs < 0) return 0;
   if (n_jobs == 0) return 1;
   if (!jobs || n_samples <= 0 || !samples) return 0;

   ctx->h_nodal_jobs = (fdtd_nodal_job *)malloc((size_t)n_jobs * sizeof(fdtd_nodal_job));
   if (!ctx->h_nodal_jobs) return 0;
   memcpy(ctx->h_nodal_jobs, jobs, (size_t)n_jobs * sizeof(fdtd_nodal_job));
   ctx->nodal_n_jobs = n_jobs;
   ctx->nodal_n_skip = n_skip;

   if (!check_nodal(cudaMalloc((void **)&ctx->d_nodal_samples, (size_t)n_samples * sizeof(fdtd_real)),
                    "nodal samples") ||
       !check_nodal(cudaMemcpy(ctx->d_nodal_samples, samples, (size_t)n_samples * sizeof(fdtd_real),
                               cudaMemcpyHostToDevice),
                    "nodal samples H2D") ||
       !check_nodal(cudaMalloc((void **)&ctx->d_nodal_skip_e, (size_t)n_skip * sizeof(int)), "skip e") ||
       !check_nodal(cudaMalloc((void **)&ctx->d_nodal_skip_h, (size_t)n_skip * sizeof(int)), "skip h") ||
       !check_nodal(cudaMemcpy(ctx->d_nodal_skip_e, skip_e, (size_t)n_skip * sizeof(int), cudaMemcpyHostToDevice),
                    "skip e H2D") ||
       !check_nodal(cudaMemcpy(ctx->d_nodal_skip_h, skip_h, (size_t)n_skip * sizeof(int), cudaMemcpyHostToDevice),
                    "skip h H2D")) {
      free_nodal(ctx);
      return 0;
   }
   ctx->nodal_ready = 1;
   return 1;
}

int fdtd_cuda_nodal_ready(const fdtd_cuda_ctx *ctx)
{
   return ctx && ctx->nodal_ready;
}

static int advance_phase(fdtd_cuda_ctx *ctx, fdtd_real time, int step, int electric)
{
   if (!ctx || !ctx->nodal_ready) return 0;
   int queued = 0;
   for (int s = 0; s < ctx->nodal_n_jobs; ++s) {
      const fdtd_nodal_job &job = ctx->h_nodal_jobs[s];
      int is_e = job.field_comp < 3;
      if (is_e != electric) continue;
      if (job.initial_only && step != 0) continue;
      if (!launch_job(ctx, job, time, step)) return 0;
      queued = 1;
   }
   if (!queued) return 1;
   return check_nodal(cudaDeviceSynchronize(), "nodal sync");
}

int fdtd_cuda_advance_nodal_e(fdtd_cuda_ctx *ctx, fdtd_real time, int step)
{
   return advance_phase(ctx, time, step, 1);
}

int fdtd_cuda_advance_nodal_h(fdtd_cuda_ctx *ctx, fdtd_real time, int step)
{
   return advance_phase(ctx, time, step, 0);
}
