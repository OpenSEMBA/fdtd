#include "fdtd_cuda_internal.h"

#include <cuda_runtime.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int check_wires(cudaError_t e, const char *what)
{
   if (e != cudaSuccess) {
      fprintf(stderr, "FDTD CUDA wires error in %s: %s\n", what, cudaGetErrorString(e));
      return 0;
   }
   return 1;
}

static void free_wires(fdtd_cuda_ctx *ctx)
{
   if (ctx->d_wire_segs) cudaFree(ctx->d_wire_segs);
   if (ctx->d_wire_nodes) cudaFree(ctx->d_wire_nodes);
   if (ctx->d_wire_samples) cudaFree(ctx->d_wire_samples);
   if (ctx->d_wire_current) cudaFree(ctx->d_wire_current);
   if (ctx->d_wire_current_past) cudaFree(ctx->d_wire_current_past);
   if (ctx->d_wire_charge) cudaFree(ctx->d_wire_charge);
   if (ctx->d_wire_charge_past) cudaFree(ctx->d_wire_charge_past);
   ctx->d_wire_segs = NULL;
   ctx->d_wire_nodes = NULL;
   ctx->d_wire_samples = NULL;
   ctx->d_wire_current = NULL;
   ctx->d_wire_current_past = NULL;
   ctx->d_wire_charge = NULL;
   ctx->d_wire_charge_past = NULL;
   ctx->wires_nseg = 0;
   ctx->wires_nnode = 0;
   ctx->wires_ready = 0;
}

void fdtd_cuda_free_wires(fdtd_cuda_ctx *ctx)
{
   if (!ctx) return;
   free_wires(ctx);
}

/* Host evolucion(): samples are evol(0:numus). Out of range returns 0. */
__device__ static fdtd_real wire_wave(const fdtd_real *samples, int off, int numus,
                                      fdtd_real delta, fdtd_real t)
{
   if (delta == (fdtd_real)0) return (fdtd_real)0;
   int nprev = (int)(t / delta);
   if ((nprev + 1 > numus) || (nprev + 1 <= 0)) return (fdtd_real)0;
   fdtd_real e0 = samples[off + nprev];
   fdtd_real e1 = samples[off + nprev + 1];
   return (e1 - e0) / delta * (t - (fdtd_real)nprev * delta) + e0;
}

__global__ void k_wire_charge(fdtd_wire_node *nodes, int nnode, const fdtd_real *current,
                              fdtd_real *charge, fdtd_real *charge_past,
                              const fdtd_real *samples, fdtd_real time_q)
{
   int n = blockIdx.x * blockDim.x + threadIdx.x;
   if (n >= nnode) return;
   const fdtd_wire_node nd = nodes[n];
   if (!nd.exists) return;
   charge_past[n] = charge[n];
   if (!nd.is_mur) {
      fdtd_real Iplus = (fdtd_real)0, Iminus = (fdtd_real)0;
      for (int k = 0; k < nd.n_plus && k < FDTD_WIRE_NEIGH; ++k) {
         int s = nd.plus_seg[k];
         if (s >= 0) Iplus += current[s];
      }
      for (int k = 0; k < nd.n_minus && k < FDTD_WIRE_NEIGH; ++k) {
         int s = nd.minus_seg[k];
         if (s >= 0) Iminus += current[s];
      }
      if (nd.n_minus == 1 && nd.n_plus == 0)
         Iplus = nd.is_periodic ? Iminus : -Iminus;
      if (nd.n_minus == 0 && nd.n_plus == 1)
         Iminus = nd.is_periodic ? Iplus : -Iplus;
      charge[n] = nd.cte_prop * charge_past[n] - nd.cte_plain * (Iplus - Iminus);
   }
   if (nd.has_i) {
      fdtd_real Iinc = wire_wave(samples, nd.evol_off, nd.numus, nd.deltaevol, time_q);
      charge[n] = charge[n] + nd.cte_plain * Iinc;
   }
}

__global__ void k_wire_mur(const fdtd_wire_node *nodes, int nnode,
                           fdtd_real *charge, const fdtd_real *charge_past)
{
   int n = blockIdx.x * blockDim.x + threadIdx.x;
   if (n >= nnode) return;
   const fdtd_wire_node nd = nodes[n];
   if (!nd.exists || !nd.is_mur || nd.node_inside < 0) return;
   int inn = nd.node_inside;
   charge[n] = charge_past[inn] + nd.cte_mur * (charge[inn] - charge_past[n]);
}

__global__ void k_wire_inject(fdtd_real *Ex, fdtd_real *Ey, fdtd_real *Ez,
                              int nx0, int ny0, int nz0, int ox0, int oy0, int oz0,
                              int nx1, int ny1, int nz1, int ox1, int oy1, int oz1,
                              int nx2, int ny2, int nz2, int ox2, int oy2, int oz2,
                              const fdtd_wire_seg *segs, int nseg, const fdtd_real *current)
{
   int s = blockIdx.x * blockDim.x + threadIdx.x;
   if (s >= nseg) return;
   const fdtd_wire_seg seg = segs[s];
   if (!seg.coupled) return;
   fdtd_real *F;
   int nx, ny, ox, oy, oz;
   if (seg.field_comp == 0) {
      F = Ex; nx = nx0; ny = ny0; ox = ox0; oy = oy0; oz = oz0;
   } else if (seg.field_comp == 1) {
      F = Ey; nx = nx1; ny = ny1; ox = ox1; oy = oy1; oz = oz1;
   } else {
      F = Ez; nx = nx2; ny = ny2; ox = ox2; oy = oy2; oz = oz2;
   }
   int ii = seg.i - ox, jj = seg.j - oy, kk = seg.k - oz;
   int nz = (seg.field_comp == 0) ? nz0 : (seg.field_comp == 1) ? nz1 : nz2;
   if (ii < 0 || jj < 0 || kk < 0 || ii >= nx || jj >= ny || kk >= nz) return;
   atomicAdd(&F[fdtd_idx3(ii, jj, kk, nx, ny)], -seg.cte5 * current[s]);
}

__global__ void k_wire_current(fdtd_real *Ex, fdtd_real *Ey, fdtd_real *Ez,
                               int nx0, int ny0, int nz0, int ox0, int oy0, int oz0,
                               int nx1, int ny1, int nz1, int ox1, int oy1, int oz1,
                               int nx2, int ny2, int nz2, int ox2, int oy2, int oz2,
                               const fdtd_wire_seg *segs, int nseg,
                               const fdtd_real *charge, fdtd_real *current, fdtd_real *current_past,
                               const fdtd_real *samples, fdtd_real time_i)
{
   int s = blockIdx.x * blockDim.x + threadIdx.x;
   if (s >= nseg) return;
   const fdtd_wire_seg seg = segs[s];
   current_past[s] = current[s];
   fdtd_real I = current[s];
   if (seg.is_pmc) {
      I = (fdtd_real)0;
   } else {
      fdtd_real qplus = (fdtd_real)0, qminus = (fdtd_real)0;
      if (seg.charge_plus >= 0) qplus = charge[seg.charge_plus];
      if (seg.charge_minus >= 0) qminus = charge[seg.charge_minus];
      fdtd_real q = seg.fraction_plus * qplus - seg.fraction_minus * qminus;
      I = seg.cte1 * I - seg.cte3 * q;
      if (seg.coupled) {
         fdtd_real *F;
         int nx, ny, ox, oy, oz, nz;
         if (seg.field_comp == 0) {
            F = Ex; nx = nx0; ny = ny0; nz = nz0; ox = ox0; oy = oy0; oz = oz0;
         } else if (seg.field_comp == 1) {
            F = Ey; nx = nx1; ny = ny1; nz = nz1; ox = ox1; oy = oy1; oz = oz1;
         } else {
            F = Ez; nx = nx2; ny = ny2; nz = nz2; ox = ox2; oy = oy2; oz = oz2;
         }
         int ii = seg.i - ox, jj = seg.j - oy, kk = seg.k - oz;
         if (ii >= 0 && jj >= 0 && kk >= 0 && ii < nx && jj < ny && kk < nz)
            I = I + seg.cte2 * F[fdtd_idx3(ii, jj, kk, nx, ny)];
      }
   }
   if (seg.has_v) {
      fdtd_real Vinc = wire_wave(samples, seg.evol_off, seg.numus, seg.deltaevol, time_i);
      I = I + seg.vscale * Vinc;
   }
   current[s] = I;
}

int fdtd_cuda_upload_wires(fdtd_cuda_ctx *ctx,
                           const fdtd_wire_seg *segs, int nseg,
                           const fdtd_wire_node *nodes, int nnode,
                           const fdtd_real *samples, int n_samples,
                           const fdtd_real *current, const fdtd_real *charge,
                           const fdtd_real *charge_past,
                           int ex_xi, int ex_yi, int ex_zi,
                           int ey_xi, int ey_yi, int ey_zi,
                           int ez_xi, int ez_yi, int ez_zi)
{
   if (!ctx || nseg < 0 || nnode < 0) return 0;
   free_wires(ctx);
   ctx->wire_ex[0] = ex_xi; ctx->wire_ex[1] = ex_yi; ctx->wire_ex[2] = ex_zi;
   ctx->wire_ey[0] = ey_xi; ctx->wire_ey[1] = ey_yi; ctx->wire_ey[2] = ey_zi;
   ctx->wire_ez[0] = ez_xi; ctx->wire_ez[1] = ez_yi; ctx->wire_ez[2] = ez_zi;
   ctx->wires_nseg = nseg;
   ctx->wires_nnode = nnode;
   if (nseg == 0 && nnode == 0) {
      ctx->wires_ready = 1;
      return 1;
   }
   if (nseg > 0) {
      if (!segs || !current) return 0;
      if (!check_wires(cudaMalloc((void **)&ctx->d_wire_segs, (size_t)nseg * sizeof(fdtd_wire_seg)), "segs") ||
          !check_wires(cudaMemcpy(ctx->d_wire_segs, segs, (size_t)nseg * sizeof(fdtd_wire_seg), cudaMemcpyHostToDevice), "segs h2d") ||
          !check_wires(cudaMalloc((void **)&ctx->d_wire_current, (size_t)nseg * sizeof(fdtd_real)), "I") ||
          !check_wires(cudaMalloc((void **)&ctx->d_wire_current_past, (size_t)nseg * sizeof(fdtd_real)), "Ipast") ||
          !check_wires(cudaMemcpy(ctx->d_wire_current, current, (size_t)nseg * sizeof(fdtd_real), cudaMemcpyHostToDevice), "I h2d") ||
          !check_wires(cudaMemcpy(ctx->d_wire_current_past, current, (size_t)nseg * sizeof(fdtd_real), cudaMemcpyHostToDevice), "Ipast h2d"))
         return 0;
   }
   if (nnode > 0) {
      if (!nodes || !charge || !charge_past) return 0;
      if (!check_wires(cudaMalloc((void **)&ctx->d_wire_nodes, (size_t)nnode * sizeof(fdtd_wire_node)), "nodes") ||
          !check_wires(cudaMemcpy(ctx->d_wire_nodes, nodes, (size_t)nnode * sizeof(fdtd_wire_node), cudaMemcpyHostToDevice), "nodes h2d") ||
          !check_wires(cudaMalloc((void **)&ctx->d_wire_charge, (size_t)nnode * sizeof(fdtd_real)), "Q") ||
          !check_wires(cudaMalloc((void **)&ctx->d_wire_charge_past, (size_t)nnode * sizeof(fdtd_real)), "Qpast") ||
          !check_wires(cudaMemcpy(ctx->d_wire_charge, charge, (size_t)nnode * sizeof(fdtd_real), cudaMemcpyHostToDevice), "Q h2d") ||
          !check_wires(cudaMemcpy(ctx->d_wire_charge_past, charge_past, (size_t)nnode * sizeof(fdtd_real), cudaMemcpyHostToDevice), "Qpast h2d"))
         return 0;
   }
   if (n_samples > 0) {
      if (!samples) return 0;
      if (!check_wires(cudaMalloc((void **)&ctx->d_wire_samples, (size_t)n_samples * sizeof(fdtd_real)), "samples") ||
          !check_wires(cudaMemcpy(ctx->d_wire_samples, samples, (size_t)n_samples * sizeof(fdtd_real), cudaMemcpyHostToDevice), "samples h2d"))
         return 0;
   }
   ctx->wires_ready = 1;
   return 1;
}

int fdtd_cuda_wires_ready(const fdtd_cuda_ctx *ctx)
{
   return ctx && ctx->wires_ready;
}

int fdtd_cuda_advance_wires_e(fdtd_cuda_ctx *ctx, fdtd_real time_i, fdtd_real time_q)
{
   if (!ctx || !ctx->wires_ready) return 0;
   const int BS = 128;
   if (ctx->wires_nnode > 0) {
      int g = (ctx->wires_nnode + BS - 1) / BS;
      k_wire_charge<<<g, BS>>>(ctx->d_wire_nodes, ctx->wires_nnode, ctx->d_wire_current,
                               ctx->d_wire_charge, ctx->d_wire_charge_past,
                               ctx->d_wire_samples, time_q);
      k_wire_mur<<<g, BS>>>(ctx->d_wire_nodes, ctx->wires_nnode,
                            ctx->d_wire_charge, ctx->d_wire_charge_past);
   }
   if (ctx->wires_nseg > 0) {
      int g = (ctx->wires_nseg + BS - 1) / BS;
      k_wire_inject<<<g, BS>>>(ctx->d_Ex, ctx->d_Ey, ctx->d_Ez,
                               ctx->dex.nx, ctx->dex.ny, ctx->dex.nz, ctx->wire_ex[0], ctx->wire_ex[1], ctx->wire_ex[2],
                               ctx->dey.nx, ctx->dey.ny, ctx->dey.nz, ctx->wire_ey[0], ctx->wire_ey[1], ctx->wire_ey[2],
                               ctx->dez.nx, ctx->dez.ny, ctx->dez.nz, ctx->wire_ez[0], ctx->wire_ez[1], ctx->wire_ez[2],
                               ctx->d_wire_segs, ctx->wires_nseg, ctx->d_wire_current);
      k_wire_current<<<g, BS>>>(ctx->d_Ex, ctx->d_Ey, ctx->d_Ez,
                                ctx->dex.nx, ctx->dex.ny, ctx->dex.nz, ctx->wire_ex[0], ctx->wire_ex[1], ctx->wire_ex[2],
                                ctx->dey.nx, ctx->dey.ny, ctx->dey.nz, ctx->wire_ey[0], ctx->wire_ey[1], ctx->wire_ey[2],
                                ctx->dez.nx, ctx->dez.ny, ctx->dez.nz, ctx->wire_ez[0], ctx->wire_ez[1], ctx->wire_ez[2],
                                ctx->d_wire_segs, ctx->wires_nseg, ctx->d_wire_charge,
                                ctx->d_wire_current, ctx->d_wire_current_past,
                                ctx->d_wire_samples, time_i);
   }
   return check_wires(cudaGetLastError(), "wires launch") &&
          check_wires(cudaDeviceSynchronize(), "wires sync");
}

int fdtd_cuda_download_wires(fdtd_cuda_ctx *ctx,
                             fdtd_real *current, fdtd_real *current_past,
                             fdtd_real *charge, fdtd_real *charge_past)
{
   if (!ctx || !ctx->wires_ready) return 0;
   if (ctx->wires_nseg > 0) {
      if (!current || !current_past) return 0;
      if (!check_wires(cudaMemcpy(current, ctx->d_wire_current, (size_t)ctx->wires_nseg * sizeof(fdtd_real), cudaMemcpyDeviceToHost), "I d2h") ||
          !check_wires(cudaMemcpy(current_past, ctx->d_wire_current_past, (size_t)ctx->wires_nseg * sizeof(fdtd_real), cudaMemcpyDeviceToHost), "Ipast d2h"))
         return 0;
   }
   if (ctx->wires_nnode > 0) {
      if (!charge || !charge_past) return 0;
      if (!check_wires(cudaMemcpy(charge, ctx->d_wire_charge, (size_t)ctx->wires_nnode * sizeof(fdtd_real), cudaMemcpyDeviceToHost), "Q d2h") ||
          !check_wires(cudaMemcpy(charge_past, ctx->d_wire_charge_past, (size_t)ctx->wires_nnode * sizeof(fdtd_real), cudaMemcpyDeviceToHost), "Qpast d2h"))
         return 0;
   }
   return 1;
}
