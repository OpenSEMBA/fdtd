#include "fdtd_cuda_internal.h"

#include <cuda_runtime.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int check_pw(cudaError_t e, const char *what)
{
   if (e != cudaSuccess) {
      fprintf(stderr, "FDTD CUDA PW error in %s: %s\n", what, cudaGetErrorString(e));
      return 0;
   }
   return 1;
}

static void free_pw_tables(fdtd_cuda_ctx *ctx)
{
#define FREE(p) do { if (p) cudaFree(p); p = NULL; } while (0)
   FREE(ctx->d_pw_evol);
   FREE(ctx->d_pw_deltaevol);
   FREE(ctx->d_pw_numus);
   FREE(ctx->d_pw_num_modes);
   FREE(ctx->d_pw_px);
   FREE(ctx->d_pw_py);
   FREE(ctx->d_pw_pz);
   FREE(ctx->d_pw_d0);
   FREE(ctx->d_pw_fpw);
   FREE(ctx->d_pw_faces_e);
   FREE(ctx->d_pw_faces_h);
   FREE(ctx->d_pw_still);
#undef FREE
   free(ctx->h_pw_faces_e);
   free(ctx->h_pw_faces_h);
   ctx->h_pw_faces_e = NULL;
   ctx->h_pw_faces_h = NULL;
   ctx->pw_n_waves = 0;
   ctx->pw_max_modes = 0;
   ctx->pw_max_numus = 0;
   ctx->pw_n_faces_e = 0;
   ctx->pw_n_faces_h = 0;
   ctx->pw_ready = 0;
}

static void free_pw_phys(fdtd_cuda_ctx *ctx)
{
#define FREE(p) do { if (p) cudaFree(p); p = NULL; } while (0)
   for (int f = 0; f < 6; ++f) {
      FREE(ctx->d_pw_phys_x[f]);
      FREE(ctx->d_pw_phys_y[f]);
      FREE(ctx->d_pw_phys_z[f]);
      ctx->pw_phys_nx[f] = ctx->pw_phys_ny[f] = ctx->pw_phys_nz[f] = 0;
   }
#undef FREE
}

static void free_pw(fdtd_cuda_ctx *ctx)
{
   free_pw_tables(ctx);
   free_pw_phys(ctx);
}

void fdtd_cuda_free_planewave(fdtd_cuda_ctx *ctx)
{
   if (!ctx) return;
   free_pw(ctx);
}

/* evol[wave * (max_numus+1) + sample], sample in [0, numus] */
__device__ inline fdtd_real pw_evolucion(int wave, fdtd_real t, fdtd_real d, fdtd_real cluz,
                                         const fdtd_real *evol, const fdtd_real *deltaevol,
                                         const int *numus, int max_numus, int *still)
{
   fdtd_real de = deltaevol[wave];
   int nu = numus[wave];
   long long nprev = (long long)((t - d / cluz) / de);
   if (nprev + 1 <= (long long)nu) {
      atomicOr(still, 1);
      if (nprev > 0) {
         const fdtd_real *row = evol + (size_t)wave * (size_t)(max_numus + 1);
         fdtd_real e0 = row[nprev];
         fdtd_real e1 = row[nprev + 1];
         return (e1 - e0) / de * ((t - d / cluz) - (fdtd_real)nprev * de) + e0;
      }
   }
   return (fdtd_real)0;
}

/* fpw[wave * 6 * max_modes + (nfield-1)*max_modes + mode], nfield 1..6 */
__device__ inline fdtd_real pw_incid(
   int wave, int nfield, fdtd_real time, int i, int j, int k, fdtd_real cluz,
   const fdtd_real *px, const fdtd_real *py, const fdtd_real *pz, const fdtd_real *d0,
   const fdtd_real *fpw, const int *num_modes, int max_modes,
   const fdtd_real *evol, const fdtd_real *deltaevol, const int *numus, int max_numus,
   const fdtd_real *phys_x, const fdtd_real *phys_y, const fdtd_real *phys_z,
   int x0, int y0, int z0, int *still)
{
   fdtd_real xf = phys_x[i - x0];
   fdtd_real yf = phys_y[j - y0];
   fdtd_real zf = phys_z[k - z0];
   fdtd_real ehi = 0;
   int nm = num_modes[wave];
   for (int mode = 0; mode < nm; ++mode) {
      size_t mi = (size_t)wave * (size_t)max_modes + (size_t)mode;
      fdtd_real dd = xf * px[mi] + yf * py[mi] + zf * pz[mi] - d0[mi];
      size_t fi = (size_t)wave * 6u * (size_t)max_modes +
                  (size_t)(nfield - 1) * (size_t)max_modes + (size_t)mode;
      ehi += fpw[fi] * pw_evolucion(wave, time, dd, cluz, evol, deltaevol, numus, max_numus, still);
   }
   return ehi;
}

__global__ void k_pw_face(
   fdtd_real *Field, int fnx, int fny,
   const fdtd_pw_face face, fdtd_real G, fdtd_real cluz,
   const fdtd_real *Id, /* 1D metric along id_axis, remapped */
   const fdtd_real *px, const fdtd_real *py, const fdtd_real *pz, const fdtd_real *d0,
   const fdtd_real *fpw, const int *num_modes, int max_modes,
   const fdtd_real *evol, const fdtd_real *deltaevol, const int *numus, int max_numus,
   const fdtd_real *phys_x, const fdtd_real *phys_y, const fdtd_real *phys_z,
   int x0, int y0, int z0, fdtd_real time, int *still)
{
   int a = blockIdx.x * blockDim.x + threadIdx.x;
   int b = blockIdx.y * blockDim.y + threadIdx.y;
   int na = face.a1 - face.a0 + 1;
   int nb = face.b1 - face.b0 + 1;
   if (a >= na || b >= nb) return;

   int i, j, k, id_idx;
   if (face.free_mode == 0) {
      i = face.fixed_abs;
      j = face.a0 + a;
      k = face.b0 + b;
      id_idx = i - face.field_xi;
   } else if (face.free_mode == 1) {
      j = face.fixed_abs;
      i = face.a0 + a;
      k = face.b0 + b;
      id_idx = j - face.field_yi;
   } else {
      k = face.fixed_abs;
      i = face.a0 + a;
      j = face.b0 + b;
      id_idx = k - face.field_zi;
   }

   int i_m = i - face.field_xi;
   int j_m = j - face.field_yi;
   int k_m = k - face.field_zi;
   fdtd_real Idv = Id[id_idx];
   fdtd_real inc = pw_incid(face.wave, face.incid_nfield, time,
                            i + face.incid_di, j + face.incid_dj, k + face.incid_dk, cluz,
                            px, py, pz, d0, fpw, num_modes, max_modes,
                            evol, deltaevol, numus, max_numus,
                            phys_x, phys_y, phys_z, x0, y0, z0, still);
   size_t idx = fdtd_idx3(i_m, j_m, k_m, fnx, fny);
   Field[idx] = Field[idx] + (fdtd_real)face.sign * G * inc * Idv;
}

static const fdtd_real *metric_ptr(fdtd_cuda_ctx *ctx, const fdtd_pw_face *f)
{
   if (f->use_e_metric) {
      if (f->id_axis == 0) return ctx->d_Idxe;
      if (f->id_axis == 1) return ctx->d_Idye;
      return ctx->d_Idze;
   }
   if (f->id_axis == 0) return ctx->d_Idxh;
   if (f->id_axis == 1) return ctx->d_Idyh;
   return ctx->d_Idzh;
}

static int launch_faces(fdtd_cuda_ctx *ctx, const fdtd_pw_face *hfaces, int n_faces,
                        fdtd_real time, fdtd_real G)
{
   if (n_faces <= 0) return 1;
   if (!hfaces) {
      fprintf(stderr, "FDTD CUDA PW: launch_faces null hfaces\n");
      return 0;
   }

   for (int f = 0; f < n_faces; ++f) {
      fdtd_pw_face face = hfaces[f];
      int na = face.a1 - face.a0 + 1;
      int nb = face.b1 - face.b0 + 1;
      if (na <= 0 || nb <= 0) continue;

      fdtd_real *Field = fdtd_cuda_field_ptr(ctx, face.field_comp);
      fdtd_dims3 fd = fdtd_cuda_field_dim(ctx, face.field_comp);
      if (!Field || fd.nx <= 0) {
         fprintf(stderr, "FDTD CUDA PW: bad field_comp=%d\n", face.field_comp);
         return 0;
      }
      int nf = face.incid_nfield - 1; /* 0..5 */
      if (nf < 0 || nf > 5) {
         fprintf(stderr, "FDTD CUDA PW: bad incid_nfield=%d\n", face.incid_nfield);
         return 0;
      }
      if (!ctx->d_pw_phys_x[nf] || !ctx->d_pw_phys_y[nf] || !ctx->d_pw_phys_z[nf]) {
         fprintf(stderr, "FDTD CUDA PW: missing phys for incid field %d (face %d)\n",
                 face.incid_nfield, f);
         return 0;
      }
      const fdtd_real *Id = metric_ptr(ctx, &face);
      if (!Id) {
         fprintf(stderr, "FDTD CUDA PW: null metric face %d\n", f);
         return 0;
      }
      dim3 block(16, 16);
      dim3 grid((na + 15) / 16, (nb + 15) / 16);
      k_pw_face<<<grid, block>>>(
         Field, fd.nx, fd.ny, face, G, ctx->pw_cluz, Id,
         ctx->d_pw_px, ctx->d_pw_py, ctx->d_pw_pz, ctx->d_pw_d0,
         ctx->d_pw_fpw, ctx->d_pw_num_modes, ctx->pw_max_modes,
         ctx->d_pw_evol, ctx->d_pw_deltaevol, ctx->d_pw_numus, ctx->pw_max_numus,
         ctx->d_pw_phys_x[nf], ctx->d_pw_phys_y[nf], ctx->d_pw_phys_z[nf],
         ctx->pw_phys_x0[nf], ctx->pw_phys_y0[nf], ctx->pw_phys_z0[nf],
         time, ctx->d_pw_still);
      if (!check_pw(cudaGetLastError(), "k_pw_face")) {
         fprintf(stderr, "FDTD CUDA PW: face %d field=%d free=%d fixed=%d a=%d..%d b=%d..%d\n",
                 f, face.field_comp, face.free_mode, face.fixed_abs, face.a0, face.a1, face.b0,
                 face.b1);
         return 0;
      }
   }
   return 1;
}

int fdtd_cuda_upload_planewave(
   fdtd_cuda_ctx *ctx, int n_waves, int max_modes, int max_numus, fdtd_real cluz,
   const int *num_modes, const int *numus, const fdtd_real *deltaevol, const fdtd_real *evol,
   const fdtd_real *px, const fdtd_real *py, const fdtd_real *pz, const fdtd_real *d0,
   const fdtd_real *fpw, const fdtd_pw_face *faces_e, int n_faces_e,
   const fdtd_pw_face *faces_h, int n_faces_h)
{
   if (!ctx || !ctx->ok) return 0;
   free_pw_tables(ctx);
   if (n_waves <= 0) {
      ctx->pw_ready = 1; /* noop success */
      return 1;
   }
   if (!num_modes || !numus || !deltaevol || !evol || !px || !py || !pz || !d0 || !fpw)
      return 0;

   ctx->pw_n_waves = n_waves;
   ctx->pw_max_modes = max_modes;
   ctx->pw_max_numus = max_numus;
   ctx->pw_cluz = cluz;
   ctx->pw_n_faces_e = n_faces_e;
   ctx->pw_n_faces_h = n_faces_h;

   size_t n_evol = (size_t)n_waves * (size_t)(max_numus + 1);
   size_t n_mode = (size_t)n_waves * (size_t)max_modes;
   size_t n_fpw = (size_t)n_waves * 6u * (size_t)max_modes;

   if (!check_pw(cudaMalloc((void **)&ctx->d_pw_evol, n_evol * sizeof(fdtd_real)), "evol") ||
       !check_pw(cudaMalloc((void **)&ctx->d_pw_deltaevol, (size_t)n_waves * sizeof(fdtd_real)),
                 "delta") ||
       !check_pw(cudaMalloc((void **)&ctx->d_pw_numus, (size_t)n_waves * sizeof(int)), "numus") ||
       !check_pw(cudaMalloc((void **)&ctx->d_pw_num_modes, (size_t)n_waves * sizeof(int)),
                 "nmodes") ||
       !check_pw(cudaMalloc((void **)&ctx->d_pw_px, n_mode * sizeof(fdtd_real)), "px") ||
       !check_pw(cudaMalloc((void **)&ctx->d_pw_py, n_mode * sizeof(fdtd_real)), "py") ||
       !check_pw(cudaMalloc((void **)&ctx->d_pw_pz, n_mode * sizeof(fdtd_real)), "pz") ||
       !check_pw(cudaMalloc((void **)&ctx->d_pw_d0, n_mode * sizeof(fdtd_real)), "d0") ||
       !check_pw(cudaMalloc((void **)&ctx->d_pw_fpw, n_fpw * sizeof(fdtd_real)), "fpw") ||
       !check_pw(cudaMalloc((void **)&ctx->d_pw_still, sizeof(int)), "still"))
      return 0;

   if (!check_pw(cudaMemcpy(ctx->d_pw_evol, evol, n_evol * sizeof(fdtd_real), cudaMemcpyHostToDevice),
                 "evol H2D") ||
       !check_pw(cudaMemcpy(ctx->d_pw_deltaevol, deltaevol, (size_t)n_waves * sizeof(fdtd_real),
                            cudaMemcpyHostToDevice),
                 "delta H2D") ||
       !check_pw(cudaMemcpy(ctx->d_pw_numus, numus, (size_t)n_waves * sizeof(int),
                            cudaMemcpyHostToDevice),
                 "numus H2D") ||
       !check_pw(cudaMemcpy(ctx->d_pw_num_modes, num_modes, (size_t)n_waves * sizeof(int),
                            cudaMemcpyHostToDevice),
                 "nmodes H2D") ||
       !check_pw(cudaMemcpy(ctx->d_pw_px, px, n_mode * sizeof(fdtd_real), cudaMemcpyHostToDevice),
                 "px H2D") ||
       !check_pw(cudaMemcpy(ctx->d_pw_py, py, n_mode * sizeof(fdtd_real), cudaMemcpyHostToDevice),
                 "py H2D") ||
       !check_pw(cudaMemcpy(ctx->d_pw_pz, pz, n_mode * sizeof(fdtd_real), cudaMemcpyHostToDevice),
                 "pz H2D") ||
       !check_pw(cudaMemcpy(ctx->d_pw_d0, d0, n_mode * sizeof(fdtd_real), cudaMemcpyHostToDevice),
                 "d0 H2D") ||
       !check_pw(cudaMemcpy(ctx->d_pw_fpw, fpw, n_fpw * sizeof(fdtd_real), cudaMemcpyHostToDevice),
                 "fpw H2D"))
      return 0;

   if (n_faces_e > 0) {
      if (!faces_e) return 0;
      ctx->h_pw_faces_e = (fdtd_pw_face *)malloc((size_t)n_faces_e * sizeof(fdtd_pw_face));
      if (!ctx->h_pw_faces_e) return 0;
      memcpy(ctx->h_pw_faces_e, faces_e, (size_t)n_faces_e * sizeof(fdtd_pw_face));
      if (!check_pw(cudaMalloc((void **)&ctx->d_pw_faces_e, (size_t)n_faces_e * sizeof(fdtd_pw_face)),
                    "faces_e") ||
          !check_pw(cudaMemcpy(ctx->d_pw_faces_e, faces_e,
                               (size_t)n_faces_e * sizeof(fdtd_pw_face), cudaMemcpyHostToDevice),
                    "faces_e H2D"))
         return 0;
   }
   if (n_faces_h > 0) {
      if (!faces_h) return 0;
      ctx->h_pw_faces_h = (fdtd_pw_face *)malloc((size_t)n_faces_h * sizeof(fdtd_pw_face));
      if (!ctx->h_pw_faces_h) return 0;
      memcpy(ctx->h_pw_faces_h, faces_h, (size_t)n_faces_h * sizeof(fdtd_pw_face));
      if (!check_pw(cudaMalloc((void **)&ctx->d_pw_faces_h, (size_t)n_faces_h * sizeof(fdtd_pw_face)),
                    "faces_h") ||
          !check_pw(cudaMemcpy(ctx->d_pw_faces_h, faces_h,
                               (size_t)n_faces_h * sizeof(fdtd_pw_face), cudaMemcpyHostToDevice),
                    "faces_h H2D"))
         return 0;
   }

   ctx->pw_ready = 1;
   return 1;
}

int fdtd_cuda_upload_planewave_phys(fdtd_cuda_ctx *ctx, int field /*0..5*/, int axis /*0=x,1=y,2=z*/,
                                    int base_abs, int n, const fdtd_real *host)
{
   if (!ctx || !ctx->ok || field < 0 || field > 5 || axis < 0 || axis > 2 || n <= 0 || !host)
      return 0;
   fdtd_real **dst = (axis == 0)   ? &ctx->d_pw_phys_x[field]
                     : (axis == 1) ? &ctx->d_pw_phys_y[field]
                                   : &ctx->d_pw_phys_z[field];
   int *nslot = (axis == 0)   ? &ctx->pw_phys_nx[field]
                : (axis == 1) ? &ctx->pw_phys_ny[field]
                              : &ctx->pw_phys_nz[field];
   int *b0 = (axis == 0)   ? &ctx->pw_phys_x0[field]
             : (axis == 1) ? &ctx->pw_phys_y0[field]
                           : &ctx->pw_phys_z0[field];
   if (*dst) {
      cudaFree(*dst);
      *dst = NULL;
   }
   if (!check_pw(cudaMalloc((void **)dst, (size_t)n * sizeof(fdtd_real)), "phys") ||
       !check_pw(cudaMemcpy(*dst, host, (size_t)n * sizeof(fdtd_real), cudaMemcpyHostToDevice),
                 "phys H2D"))
      return 0;
   *nslot = n;
   *b0 = base_abs;
   return 1;
}

int fdtd_cuda_planewave_ready(const fdtd_cuda_ctx *ctx)
{
   return ctx && ctx->ok && ctx->pw_ready;
}

static int advance_pw(fdtd_cuda_ctx *ctx, fdtd_real time, int is_h, int *still_out)
{
   if (!ctx || !ctx->ok) {
      fprintf(stderr, "FDTD CUDA PW: advance bad ctx\n");
      return 0;
   }
   if (!ctx->pw_ready) {
      fprintf(stderr, "FDTD CUDA PW: advance not ready\n");
      return 0;
   }
   if (still_out) *still_out = 0;
   int n_faces = is_h ? ctx->pw_n_faces_h : ctx->pw_n_faces_e;
   if (n_faces <= 0 || ctx->pw_n_waves <= 0) {
      if (still_out) *still_out = 0;
      return 1;
   }
   if (!ctx->d_g2 || !ctx->d_gm2) {
      fprintf(stderr, "FDTD CUDA PW: missing G coeffs\n");
      return 0;
   }
   int zero = 0;
   if (!check_pw(cudaMemcpy(ctx->d_pw_still, &zero, sizeof(int), cudaMemcpyHostToDevice), "still clr"))
      return 0;

   /* G2(1) / Gm2(1): Fortran coeffs indexed 0..NumMedia */
   fdtd_real Ghost = 0;
   const fdtd_real *dG = is_h ? ctx->d_gm2 : ctx->d_g2;
   if (!check_pw(cudaMemcpy(&Ghost, dG + 1, sizeof(fdtd_real), cudaMemcpyDeviceToHost), "G1"))
      return 0;

   const fdtd_pw_face *hfaces = is_h ? ctx->h_pw_faces_h : ctx->h_pw_faces_e;
   if (!launch_faces(ctx, hfaces, n_faces, time, Ghost)) return 0;
   if (!check_pw(cudaDeviceSynchronize(), "pw sync")) return 0;

   int still = 0;
   if (!check_pw(cudaMemcpy(&still, ctx->d_pw_still, sizeof(int), cudaMemcpyDeviceToHost), "still D2H"))
      return 0;
   if (still_out) *still_out = still;
   return 1;
}

int fdtd_cuda_advance_planewave_e(fdtd_cuda_ctx *ctx, fdtd_real time, int *still_out)
{
   return advance_pw(ctx, time, 0, still_out);
}

int fdtd_cuda_advance_planewave_h(fdtd_cuda_ctx *ctx, fdtd_real time, int *still_out)
{
   return advance_pw(ctx, time, 1, still_out);
}
