#include "fdtd_cuda_internal.h"

#include <cuda_runtime.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int check(cudaError_t e, const char *what)
{
   if (e != cudaSuccess) {
      fprintf(stderr, "FDTD CUDA error in %s: %s\n", what, cudaGetErrorString(e));
      return 0;
   }
   return 1;
}

static size_t nelem(fdtd_dims3 d) { return (size_t)d.nx * (size_t)d.ny * (size_t)d.nz; }

fdtd_cuda_ctx *fdtd_cuda_create(void)
{
   fdtd_cuda_ctx *ctx = (fdtd_cuda_ctx *)calloc(1, sizeof(*ctx));
   if (!ctx) return NULL;
   int ndev = 0;
   if (!check(cudaGetDeviceCount(&ndev), "cudaGetDeviceCount") || ndev < 1) {
      free(ctx);
      return NULL;
   }
   if (!check(cudaSetDevice(0), "cudaSetDevice")) {
      free(ctx);
      return NULL;
   }
   ctx->ok = 1;
   return ctx;
}

void fdtd_cuda_destroy(fdtd_cuda_ctx *ctx)
{
   if (!ctx) return;
#define FREE(p) do { if (p) cudaFree(p); p = NULL; } while (0)
   FREE(ctx->d_Ex); FREE(ctx->d_Ey); FREE(ctx->d_Ez);
   FREE(ctx->d_Hx); FREE(ctx->d_Hy); FREE(ctx->d_Hz);
   FREE(ctx->d_mEx); FREE(ctx->d_mEy); FREE(ctx->d_mEz);
   FREE(ctx->d_mHx); FREE(ctx->d_mHy); FREE(ctx->d_mHz);
   FREE(ctx->d_g1); FREE(ctx->d_g2); FREE(ctx->d_gm1); FREE(ctx->d_gm2);
   FREE(ctx->d_Idxh); FREE(ctx->d_Idyh); FREE(ctx->d_Idzh);
   FREE(ctx->d_Idxe); FREE(ctx->d_Idye); FREE(ctx->d_Idze);
   for (int i = 0; i < FDTD_CUDA_PSI_SLOTS; ++i) FREE(ctx->d_psi[i]);
   for (int i = 0; i < FDTD_CUDA_CPML1D; ++i) FREE(ctx->d_cpml1d[i]);
   FREE(ctx->d_probe_comp); FREE(ctx->d_probe_i); FREE(ctx->d_probe_j); FREE(ctx->d_probe_k);
   FREE(ctx->d_probe_out);
#undef FREE
   fdtd_cuda_free_planewave(ctx);
   fdtd_cuda_free_mur(ctx);
   free(ctx);
}

int fdtd_cuda_ok(const fdtd_cuda_ctx *ctx) { return ctx && ctx->ok; }

static int alloc_real3(fdtd_real **p, fdtd_dims3 d)
{
   return check(cudaMalloc((void **)p, nelem(d) * sizeof(fdtd_real)), "cudaMalloc field");
}
static int alloc_media3(fdtd_media **p, fdtd_dims3 d)
{
   return check(cudaMalloc((void **)p, nelem(d) * sizeof(fdtd_media)), "cudaMalloc media");
}

int fdtd_cuda_alloc_fields(fdtd_cuda_ctx *ctx,
                           fdtd_dims3 ex, fdtd_dims3 ey, fdtd_dims3 ez,
                           fdtd_dims3 hx, fdtd_dims3 hy, fdtd_dims3 hz)
{
   if (!ctx) return 0;
   ctx->dex = ex; ctx->dey = ey; ctx->dez = ez;
   ctx->dhx = hx; ctx->dhy = hy; ctx->dhz = hz;
   return alloc_real3(&ctx->d_Ex, ex) && alloc_real3(&ctx->d_Ey, ey) &&
          alloc_real3(&ctx->d_Ez, ez) && alloc_real3(&ctx->d_Hx, hx) &&
          alloc_real3(&ctx->d_Hy, hy) && alloc_real3(&ctx->d_Hz, hz) &&
          check(cudaMemset(ctx->d_Ex, 0, nelem(ex) * sizeof(fdtd_real)), "memset Ex") &&
          check(cudaMemset(ctx->d_Ey, 0, nelem(ey) * sizeof(fdtd_real)), "memset Ey") &&
          check(cudaMemset(ctx->d_Ez, 0, nelem(ez) * sizeof(fdtd_real)), "memset Ez") &&
          check(cudaMemset(ctx->d_Hx, 0, nelem(hx) * sizeof(fdtd_real)), "memset Hx") &&
          check(cudaMemset(ctx->d_Hy, 0, nelem(hy) * sizeof(fdtd_real)), "memset Hy") &&
          check(cudaMemset(ctx->d_Hz, 0, nelem(hz) * sizeof(fdtd_real)), "memset Hz");
}

int fdtd_cuda_alloc_media(fdtd_cuda_ctx *ctx,
                          fdtd_dims3 mex, fdtd_dims3 mey, fdtd_dims3 mez,
                          fdtd_dims3 mhx, fdtd_dims3 mhy, fdtd_dims3 mhz)
{
   if (!ctx) return 0;
   ctx->mex = mex; ctx->mey = mey; ctx->mez = mez;
   ctx->mhx = mhx; ctx->mhy = mhy; ctx->mhz = mhz;
   return alloc_media3(&ctx->d_mEx, mex) && alloc_media3(&ctx->d_mEy, mey) &&
          alloc_media3(&ctx->d_mEz, mez) && alloc_media3(&ctx->d_mHx, mhx) &&
          alloc_media3(&ctx->d_mHy, mhy) && alloc_media3(&ctx->d_mHz, mhz);
}

int fdtd_cuda_alloc_coeffs(fdtd_cuda_ctx *ctx, int num_media)
{
   if (!ctx || num_media < 0) return 0;
   ctx->num_media = num_media;
   size_t n = (size_t)num_media + 1;
   return check(cudaMalloc((void **)&ctx->d_g1, n * sizeof(fdtd_real)), "g1") &&
          check(cudaMalloc((void **)&ctx->d_g2, n * sizeof(fdtd_real)), "g2") &&
          check(cudaMalloc((void **)&ctx->d_gm1, n * sizeof(fdtd_real)), "gm1") &&
          check(cudaMalloc((void **)&ctx->d_gm2, n * sizeof(fdtd_real)), "gm2");
}

int fdtd_cuda_alloc_metrics(fdtd_cuda_ctx *ctx,
                            int n_idxh, int n_idyh, int n_idzh,
                            int n_idxe, int n_idye, int n_idze)
{
   if (!ctx) return 0;
   ctx->n_idxh = n_idxh; ctx->n_idyh = n_idyh; ctx->n_idzh = n_idzh;
   ctx->n_idxe = n_idxe; ctx->n_idye = n_idye; ctx->n_idze = n_idze;
   return check(cudaMalloc((void **)&ctx->d_Idxh, (size_t)n_idxh * sizeof(fdtd_real)), "Idxh") &&
          check(cudaMalloc((void **)&ctx->d_Idyh, (size_t)n_idyh * sizeof(fdtd_real)), "Idyh") &&
          check(cudaMalloc((void **)&ctx->d_Idzh, (size_t)n_idzh * sizeof(fdtd_real)), "Idzh") &&
          check(cudaMalloc((void **)&ctx->d_Idxe, (size_t)n_idxe * sizeof(fdtd_real)), "Idxe") &&
          check(cudaMalloc((void **)&ctx->d_Idye, (size_t)n_idye * sizeof(fdtd_real)), "Idye") &&
          check(cudaMalloc((void **)&ctx->d_Idze, (size_t)n_idze * sizeof(fdtd_real)), "Idze");
}

#define UP(dst, src, n, t) check(cudaMemcpy(dst, src, (size_t)(n) * sizeof(t), cudaMemcpyHostToDevice), "H2D")
#define DN(dst, src, n, t) check(cudaMemcpy(dst, src, (size_t)(n) * sizeof(t), cudaMemcpyDeviceToHost), "D2H")

int fdtd_cuda_upload_fields(fdtd_cuda_ctx *ctx,
                            const fdtd_real *Ex, const fdtd_real *Ey, const fdtd_real *Ez,
                            const fdtd_real *Hx, const fdtd_real *Hy, const fdtd_real *Hz)
{
   if (!ctx) return 0;
   return UP(ctx->d_Ex, Ex, nelem(ctx->dex), fdtd_real) &&
          UP(ctx->d_Ey, Ey, nelem(ctx->dey), fdtd_real) &&
          UP(ctx->d_Ez, Ez, nelem(ctx->dez), fdtd_real) &&
          UP(ctx->d_Hx, Hx, nelem(ctx->dhx), fdtd_real) &&
          UP(ctx->d_Hy, Hy, nelem(ctx->dhy), fdtd_real) &&
          UP(ctx->d_Hz, Hz, nelem(ctx->dhz), fdtd_real);
}

int fdtd_cuda_download_fields(fdtd_cuda_ctx *ctx,
                              fdtd_real *Ex, fdtd_real *Ey, fdtd_real *Ez,
                              fdtd_real *Hx, fdtd_real *Hy, fdtd_real *Hz)
{
   if (!ctx) return 0;
   return DN(Ex, ctx->d_Ex, nelem(ctx->dex), fdtd_real) &&
          DN(Ey, ctx->d_Ey, nelem(ctx->dey), fdtd_real) &&
          DN(Ez, ctx->d_Ez, nelem(ctx->dez), fdtd_real) &&
          DN(Hx, ctx->d_Hx, nelem(ctx->dhx), fdtd_real) &&
          DN(Hy, ctx->d_Hy, nelem(ctx->dhy), fdtd_real) &&
          DN(Hz, ctx->d_Hz, nelem(ctx->dhz), fdtd_real);
}

/* Fortran column-major: index = i + j*nx + k*nx*ny. One cudaMemcpy3D for the box. */
static int copy_field_box(fdtd_real *dst, const fdtd_real *src, fdtd_dims3 d,
                          fdtd_ibox b, int host_to_device)
{
   int is = b.is < 0 ? 0 : b.is;
   int ie = b.ie >= d.nx ? d.nx - 1 : b.ie;
   int js = b.js < 0 ? 0 : b.js;
   int je = b.je >= d.ny ? d.ny - 1 : b.je;
   int ks = b.ks < 0 ? 0 : b.ks;
   int ke = b.ke >= d.nz ? d.nz - 1 : b.ke;
   if (is > ie || js > je || ks > ke) return 1;

   /* Full array → single bulk copy */
   if (is == 0 && js == 0 && ks == 0 && ie == d.nx - 1 && je == d.ny - 1 && ke == d.nz - 1) {
      enum cudaMemcpyKind kind = host_to_device ? cudaMemcpyHostToDevice : cudaMemcpyDeviceToHost;
      return check(cudaMemcpy(dst, src, nelem(d) * sizeof(fdtd_real), kind), "box full memcpy");
   }

   size_t width = (size_t)(ie - is + 1) * sizeof(fdtd_real);
   size_t height = (size_t)(je - js + 1);
   size_t depth = (size_t)(ke - ks + 1);
   size_t pitch = (size_t)d.nx * sizeof(fdtd_real);
   size_t spitch = pitch; /* same layout host and device */
   size_t offset = ((size_t)is + (size_t)js * (size_t)d.nx +
                    (size_t)ks * (size_t)d.nx * (size_t)d.ny) * sizeof(fdtd_real);

   cudaMemcpy3DParms p;
   memset(&p, 0, sizeof(p));
   p.extent = make_cudaExtent(width, height, depth);
   p.kind = host_to_device ? cudaMemcpyHostToDevice : cudaMemcpyDeviceToHost;
   if (host_to_device) {
      p.srcPtr = make_cudaPitchedPtr((void *)((char *)src + offset), spitch, (size_t)d.nx, (size_t)d.ny);
      p.dstPtr = make_cudaPitchedPtr((char *)dst + offset, pitch, (size_t)d.nx, (size_t)d.ny);
   } else {
      p.srcPtr = make_cudaPitchedPtr((char *)src + offset, pitch, (size_t)d.nx, (size_t)d.ny);
      p.dstPtr = make_cudaPitchedPtr((void *)((char *)dst + offset), spitch, (size_t)d.nx, (size_t)d.ny);
   }
   return check(cudaMemcpy3D(&p), "box memcpy3D");
}

int fdtd_cuda_upload_fields_box(fdtd_cuda_ctx *ctx,
                                const fdtd_real *Ex, const fdtd_real *Ey, const fdtd_real *Ez,
                                const fdtd_real *Hx, const fdtd_real *Hy, const fdtd_real *Hz,
                                const fdtd_ibox *ex, const fdtd_ibox *ey, const fdtd_ibox *ez,
                                const fdtd_ibox *hx, const fdtd_ibox *hy, const fdtd_ibox *hz)
{
   if (!ctx || !ex || !ey || !ez || !hx || !hy || !hz) return 0;
   return copy_field_box(ctx->d_Ex, Ex, ctx->dex, *ex, 1) &&
          copy_field_box(ctx->d_Ey, Ey, ctx->dey, *ey, 1) &&
          copy_field_box(ctx->d_Ez, Ez, ctx->dez, *ez, 1) &&
          copy_field_box(ctx->d_Hx, Hx, ctx->dhx, *hx, 1) &&
          copy_field_box(ctx->d_Hy, Hy, ctx->dhy, *hy, 1) &&
          copy_field_box(ctx->d_Hz, Hz, ctx->dhz, *hz, 1);
}

int fdtd_cuda_download_fields_box(fdtd_cuda_ctx *ctx,
                                  fdtd_real *Ex, fdtd_real *Ey, fdtd_real *Ez,
                                  fdtd_real *Hx, fdtd_real *Hy, fdtd_real *Hz,
                                  const fdtd_ibox *ex, const fdtd_ibox *ey, const fdtd_ibox *ez,
                                  const fdtd_ibox *hx, const fdtd_ibox *hy, const fdtd_ibox *hz)
{
   if (!ctx || !ex || !ey || !ez || !hx || !hy || !hz) return 0;
   return copy_field_box(Ex, ctx->d_Ex, ctx->dex, *ex, 0) &&
          copy_field_box(Ey, ctx->d_Ey, ctx->dey, *ey, 0) &&
          copy_field_box(Ez, ctx->d_Ez, ctx->dez, *ez, 0) &&
          copy_field_box(Hx, ctx->d_Hx, ctx->dhx, *hx, 0) &&
          copy_field_box(Hy, ctx->d_Hy, ctx->dhy, *hy, 0) &&
          copy_field_box(Hz, ctx->d_Hz, ctx->dhz, *hz, 0);
}

int fdtd_cuda_upload_media(fdtd_cuda_ctx *ctx,
                           const fdtd_media *mex, const fdtd_media *mey, const fdtd_media *mez,
                           const fdtd_media *mhx, const fdtd_media *mhy, const fdtd_media *mhz)
{
   if (!ctx) return 0;
   return UP(ctx->d_mEx, mex, nelem(ctx->mex), fdtd_media) &&
          UP(ctx->d_mEy, mey, nelem(ctx->mey), fdtd_media) &&
          UP(ctx->d_mEz, mez, nelem(ctx->mez), fdtd_media) &&
          UP(ctx->d_mHx, mhx, nelem(ctx->mhx), fdtd_media) &&
          UP(ctx->d_mHy, mhy, nelem(ctx->mhy), fdtd_media) &&
          UP(ctx->d_mHz, mhz, nelem(ctx->mhz), fdtd_media);
}

int fdtd_cuda_upload_coeffs(fdtd_cuda_ctx *ctx,
                            const fdtd_real *g1, const fdtd_real *g2,
                            const fdtd_real *gm1, const fdtd_real *gm2, int num_media)
{
   if (!ctx) return 0;
   size_t n = (size_t)num_media + 1;
   return UP(ctx->d_g1, g1, n, fdtd_real) && UP(ctx->d_g2, g2, n, fdtd_real) &&
          UP(ctx->d_gm1, gm1, n, fdtd_real) && UP(ctx->d_gm2, gm2, n, fdtd_real);
}

int fdtd_cuda_upload_metrics(fdtd_cuda_ctx *ctx,
                             const fdtd_real *Idxh, const fdtd_real *Idyh, const fdtd_real *Idzh,
                             const fdtd_real *Idxe, const fdtd_real *Idye, const fdtd_real *Idze)
{
   if (!ctx) return 0;
   return UP(ctx->d_Idxh, Idxh, ctx->n_idxh, fdtd_real) &&
          UP(ctx->d_Idyh, Idyh, ctx->n_idyh, fdtd_real) &&
          UP(ctx->d_Idzh, Idzh, ctx->n_idzh, fdtd_real) &&
          UP(ctx->d_Idxe, Idxe, ctx->n_idxe, fdtd_real) &&
          UP(ctx->d_Idye, Idye, ctx->n_idye, fdtd_real) &&
          UP(ctx->d_Idze, Idze, ctx->n_idze, fdtd_real);
}

int fdtd_cuda_alloc_psi(fdtd_cuda_ctx *ctx, int slot, int n_elem)
{
   if (!ctx || slot < 0 || slot >= FDTD_CUDA_PSI_SLOTS || n_elem <= 0) return 0;
   if (ctx->d_psi[slot]) cudaFree(ctx->d_psi[slot]);
   ctx->psi_n[slot] = n_elem;
   return check(cudaMalloc((void **)&ctx->d_psi[slot], (size_t)n_elem * sizeof(fdtd_real)), "psi") &&
          check(cudaMemset(ctx->d_psi[slot], 0, (size_t)n_elem * sizeof(fdtd_real)), "psi0");
}

int fdtd_cuda_upload_psi(fdtd_cuda_ctx *ctx, int slot, const fdtd_real *psi, int n_elem)
{
   if (!ctx || slot < 0 || slot >= FDTD_CUDA_PSI_SLOTS || !ctx->d_psi[slot]) return 0;
   return UP(ctx->d_psi[slot], psi, n_elem, fdtd_real);
}

int fdtd_cuda_download_psi(fdtd_cuda_ctx *ctx, int slot, fdtd_real *psi, int n_elem)
{
   if (!ctx || slot < 0 || slot >= FDTD_CUDA_PSI_SLOTS || !ctx->d_psi[slot]) return 0;
   return DN(psi, ctx->d_psi[slot], n_elem, fdtd_real);
}

int fdtd_cuda_alloc_cpml_1d(fdtd_cuda_ctx *ctx, int which, int n, const fdtd_real *host)
{
   if (!ctx || which < 0 || which >= FDTD_CUDA_CPML1D || n <= 0) return 0;
   if (ctx->d_cpml1d[which]) cudaFree(ctx->d_cpml1d[which]);
   ctx->cpml1d_n[which] = n;
   return check(cudaMalloc((void **)&ctx->d_cpml1d[which], (size_t)n * sizeof(fdtd_real)), "cpml1d") &&
          UP(ctx->d_cpml1d[which], host, n, fdtd_real);
}

/* Exposed for yee/cpml translation units */
extern "C" {
fdtd_real *fdtd_cuda_ptr_Ex(fdtd_cuda_ctx *c) { return c->d_Ex; }
fdtd_real *fdtd_cuda_ptr_Ey(fdtd_cuda_ctx *c) { return c->d_Ey; }
fdtd_real *fdtd_cuda_ptr_Ez(fdtd_cuda_ctx *c) { return c->d_Ez; }
fdtd_real *fdtd_cuda_ptr_Hx(fdtd_cuda_ctx *c) { return c->d_Hx; }
fdtd_real *fdtd_cuda_ptr_Hy(fdtd_cuda_ctx *c) { return c->d_Hy; }
fdtd_real *fdtd_cuda_ptr_Hz(fdtd_cuda_ctx *c) { return c->d_Hz; }
fdtd_media *fdtd_cuda_ptr_mEx(fdtd_cuda_ctx *c) { return c->d_mEx; }
fdtd_media *fdtd_cuda_ptr_mEy(fdtd_cuda_ctx *c) { return c->d_mEy; }
fdtd_media *fdtd_cuda_ptr_mEz(fdtd_cuda_ctx *c) { return c->d_mEz; }
fdtd_media *fdtd_cuda_ptr_mHx(fdtd_cuda_ctx *c) { return c->d_mHx; }
fdtd_media *fdtd_cuda_ptr_mHy(fdtd_cuda_ctx *c) { return c->d_mHy; }
fdtd_media *fdtd_cuda_ptr_mHz(fdtd_cuda_ctx *c) { return c->d_mHz; }
fdtd_real *fdtd_cuda_ptr_g1(fdtd_cuda_ctx *c) { return c->d_g1; }
fdtd_real *fdtd_cuda_ptr_g2(fdtd_cuda_ctx *c) { return c->d_g2; }
fdtd_real *fdtd_cuda_ptr_gm1(fdtd_cuda_ctx *c) { return c->d_gm1; }
fdtd_real *fdtd_cuda_ptr_gm2(fdtd_cuda_ctx *c) { return c->d_gm2; }
fdtd_real *fdtd_cuda_ptr_Idxh(fdtd_cuda_ctx *c) { return c->d_Idxh; }
fdtd_real *fdtd_cuda_ptr_Idyh(fdtd_cuda_ctx *c) { return c->d_Idyh; }
fdtd_real *fdtd_cuda_ptr_Idzh(fdtd_cuda_ctx *c) { return c->d_Idzh; }
fdtd_real *fdtd_cuda_ptr_Idxe(fdtd_cuda_ctx *c) { return c->d_Idxe; }
fdtd_real *fdtd_cuda_ptr_Idye(fdtd_cuda_ctx *c) { return c->d_Idye; }
fdtd_real *fdtd_cuda_ptr_Idze(fdtd_cuda_ctx *c) { return c->d_Idze; }
fdtd_dims3 fdtd_cuda_dim_Ex(fdtd_cuda_ctx *c) { return c->dex; }
fdtd_dims3 fdtd_cuda_dim_Ey(fdtd_cuda_ctx *c) { return c->dey; }
fdtd_dims3 fdtd_cuda_dim_Ez(fdtd_cuda_ctx *c) { return c->dez; }
fdtd_dims3 fdtd_cuda_dim_Hx(fdtd_cuda_ctx *c) { return c->dhx; }
fdtd_dims3 fdtd_cuda_dim_Hy(fdtd_cuda_ctx *c) { return c->dhy; }
fdtd_dims3 fdtd_cuda_dim_Hz(fdtd_cuda_ctx *c) { return c->dhz; }
fdtd_real *fdtd_cuda_ptr_psi(fdtd_cuda_ctx *c, int s) { return c->d_psi[s]; }
fdtd_real *fdtd_cuda_ptr_cpml1d(fdtd_cuda_ctx *c, int w) { return c->d_cpml1d[w]; }

fdtd_real *fdtd_cuda_field_ptr(fdtd_cuda_ctx *c, int comp)
{
   switch (comp) {
   case 0: return c->d_Ex; case 1: return c->d_Ey; case 2: return c->d_Ez;
   case 3: return c->d_Hx; case 4: return c->d_Hy; case 5: return c->d_Hz;
   default: return NULL;
   }
}
fdtd_dims3 fdtd_cuda_field_dim(fdtd_cuda_ctx *c, int comp)
{
   switch (comp) {
   case 0: return c->dex; case 1: return c->dey; case 2: return c->dez;
   case 3: return c->dhx; case 4: return c->dhy; case 5: return c->dhz;
   default: { fdtd_dims3 z = {0,0,0}; return z; }
   }
}
fdtd_media *fdtd_cuda_media_ptr(fdtd_cuda_ctx *c, int comp)
{
   switch (comp) {
   case 0: return c->d_mEx; case 1: return c->d_mEy; case 2: return c->d_mEz;
   case 3: return c->d_mHx; case 4: return c->d_mHy; case 5: return c->d_mHz;
   default: return NULL;
   }
}
fdtd_dims3 fdtd_cuda_media_dim(fdtd_cuda_ctx *c, int comp)
{
   switch (comp) {
   case 0: return c->mex; case 1: return c->mey; case 2: return c->mez;
   case 3: return c->mhx; case 4: return c->mhy; case 5: return c->mhz;
   default: { fdtd_dims3 z = {0,0,0}; return z; }
   }
}
} /* extern "C" */
