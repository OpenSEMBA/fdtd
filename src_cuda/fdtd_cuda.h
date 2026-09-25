#pragma once
/*
 * Device-resident Yee + CPML API for semba-fdtd (P0/P1 CUDA slice).
 * Fields stay on the GPU across timesteps; host sync is explicit.
 */
#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>
#include <stdint.h>

/* Match CompileWithInt2 / CompileWithReal4 defaults. */
typedef float fdtd_real;
typedef int16_t fdtd_media;

typedef struct {
   int nx, ny, nz; /* full allocated dims (0-based size) */
} fdtd_dims3;

typedef struct {
   int is, ie, js, je, ks, ke; /* inclusive sweep in remapped indices */
} fdtd_ibox;

typedef struct fdtd_cuda_ctx fdtd_cuda_ctx;

/* Lifecycle */
fdtd_cuda_ctx *fdtd_cuda_create(void);
void fdtd_cuda_destroy(fdtd_cuda_ctx *ctx);
int fdtd_cuda_ok(const fdtd_cuda_ctx *ctx);

/* Allocate / free device buffers (sizes in elements). */
int fdtd_cuda_alloc_fields(fdtd_cuda_ctx *ctx,
                           fdtd_dims3 ex, fdtd_dims3 ey, fdtd_dims3 ez,
                           fdtd_dims3 hx, fdtd_dims3 hy, fdtd_dims3 hz);
int fdtd_cuda_alloc_media(fdtd_cuda_ctx *ctx,
                          fdtd_dims3 mex, fdtd_dims3 mey, fdtd_dims3 mez,
                          fdtd_dims3 mhx, fdtd_dims3 mhy, fdtd_dims3 mhz);
int fdtd_cuda_alloc_coeffs(fdtd_cuda_ctx *ctx, int num_media);
int fdtd_cuda_alloc_metrics(fdtd_cuda_ctx *ctx,
                            int n_idxh, int n_idyh, int n_idzh,
                            int n_idxe, int n_idye, int n_idze);

/* Host <-> device bulk transfers (column-major Fortran layout). */
int fdtd_cuda_upload_fields(fdtd_cuda_ctx *ctx,
                            const fdtd_real *Ex, const fdtd_real *Ey, const fdtd_real *Ez,
                            const fdtd_real *Hx, const fdtd_real *Hy, const fdtd_real *Hz);
int fdtd_cuda_download_fields(fdtd_cuda_ctx *ctx,
                              fdtd_real *Ex, fdtd_real *Ey, fdtd_real *Ez,
                              fdtd_real *Hx, fdtd_real *Hy, fdtd_real *Hz);
/* Boxed H2D/D2H: each ibox is inclusive remapped indices into that field array. */
int fdtd_cuda_upload_fields_box(fdtd_cuda_ctx *ctx,
                                const fdtd_real *Ex, const fdtd_real *Ey, const fdtd_real *Ez,
                                const fdtd_real *Hx, const fdtd_real *Hy, const fdtd_real *Hz,
                                const fdtd_ibox *ex, const fdtd_ibox *ey, const fdtd_ibox *ez,
                                const fdtd_ibox *hx, const fdtd_ibox *hy, const fdtd_ibox *hz);
int fdtd_cuda_download_fields_box(fdtd_cuda_ctx *ctx,
                                  fdtd_real *Ex, fdtd_real *Ey, fdtd_real *Ez,
                                  fdtd_real *Hx, fdtd_real *Hy, fdtd_real *Hz,
                                  const fdtd_ibox *ex, const fdtd_ibox *ey, const fdtd_ibox *ez,
                                  const fdtd_ibox *hx, const fdtd_ibox *hy, const fdtd_ibox *hz);
int fdtd_cuda_upload_media(fdtd_cuda_ctx *ctx,
                           const fdtd_media *mex, const fdtd_media *mey, const fdtd_media *mez,
                           const fdtd_media *mhx, const fdtd_media *mhy, const fdtd_media *mhz);
int fdtd_cuda_upload_coeffs(fdtd_cuda_ctx *ctx,
                            const fdtd_real *g1, const fdtd_real *g2,
                            const fdtd_real *gm1, const fdtd_real *gm2, int num_media);
int fdtd_cuda_upload_metrics(fdtd_cuda_ctx *ctx,
                             const fdtd_real *Idxh, const fdtd_real *Idyh, const fdtd_real *Idzh,
                             const fdtd_real *Idxe, const fdtd_real *Idye, const fdtd_real *Idze);

/* P0: bulk Yee updates (operate on device-resident buffers). */
int fdtd_cuda_advance_ex(fdtd_cuda_ctx *ctx, fdtd_ibox sweep);
int fdtd_cuda_advance_ey(fdtd_cuda_ctx *ctx, fdtd_ibox sweep);
int fdtd_cuda_advance_ez(fdtd_cuda_ctx *ctx, fdtd_ibox sweep);
int fdtd_cuda_advance_hx(fdtd_cuda_ctx *ctx, fdtd_ibox sweep);
int fdtd_cuda_advance_hy(fdtd_cuda_ctx *ctx, fdtd_ibox sweep);
int fdtd_cuda_advance_hz(fdtd_cuda_ctx *ctx, fdtd_ibox sweep);

/*
 * P1: one CPML face contribution on a field component.
 * Updates Psi in place and adds ±G2(medio)*Psi into E (or ±Gm2 into H).
 *
 * Layout:
 *   E/H arrays: Fortran column-major, sizes from alloc_fields.
 *   Psi: sized (nx_psi * ny_psi * nz_psi) with absolute-ish indexing remapped
 *        so element (i,j,k) with i in [xi,xe] maps to
 *        (i-xi) + (j-yi)*nx_psi + (k-zi)*nx_psi*ny_psi.
 *   P_b, P_c: 1D coeffs indexed by the free coordinate (x, y, or z) using
 *             the absolute index used in the Fortran loops (passed as base
 *             offset p_base so device index = abs_index - p_base).
 *
 * h_diff_axis: 0=x, 1=y, 2=z — which neighbour direction for (H_a - H_b).
 * h_sign: +1 or -1 applied to G2*Psi when updating the field.
 * use_fixed_medio: if non-zero, medio_fixed is used instead of media lookup.
 */
int fdtd_cuda_alloc_psi(fdtd_cuda_ctx *ctx, int slot, int n_elem);
int fdtd_cuda_upload_psi(fdtd_cuda_ctx *ctx, int slot, const fdtd_real *psi, int n_elem);
int fdtd_cuda_download_psi(fdtd_cuda_ctx *ctx, int slot, fdtd_real *psi, int n_elem);
int fdtd_cuda_alloc_cpml_1d(fdtd_cuda_ctx *ctx, int which, int n, const fdtd_real *host);
/* which: 0=P_be_x .. 5=P_be_z, 6=P_ce_x .. 11=P_ce_z,
          12=P_bm_x .. 17=P_bm_z, 18=P_cm_x .. 23=P_cm_z */

typedef struct {
   int field_comp;   /* 0=Ex..2=Ez, 3=Hx..5=Hz */
   int psi_slot;
   int media_comp;   /* 0=Ex..5=Hz media plane matching field */
   int h_comp_a;     /* 0=Ex..5=Hz for the + neighbour sample */
   int h_comp_b;     /* neighbour for minus side (same plane, offset) */
   int h_diff_axis;  /* 0=x,1=y,2=z */
   int p_b_which;    /* index into CPML 1D table (0..23) */
   int p_c_which;
   int p_base;       /* subtract from absolute free index for 1D coeff */
   int free_axis;    /* which absolute index feeds P_b/P_c (0=x,1=y,2=z) */
   int xi, xe, yi, ye, zi, ze; /* absolute region inclusive */
   int e_xi, e_yi, e_zi;       /* field array absolute origin (Alloc lower) */
   int psi_xi, psi_yi, psi_zi; /* psi array absolute origin */
   int nx_psi, ny_psi, nz_psi;
   int h_sign;       /* +1 or -1 */
   int use_fixed_medio;
   int medio_fixed;
} fdtd_cpml_job;

int fdtd_cuda_cpml_apply(fdtd_cuda_ctx *ctx, const fdtd_cpml_job *job);

/*
 * Sparse point-probe gather: register remapped (comp,i,j,k) once, then each
 * step download only those N floats (avoids full-field D2H for point probes).
 * comp: 0=Ex .. 5=Hz. n==0 is a no-op success.
 */
int fdtd_cuda_set_point_probes(fdtd_cuda_ctx *ctx, int n,
                               const int *comp, const int *i, const int *j, const int *k);
int fdtd_cuda_gather_point_probes(fdtd_cuda_ctx *ctx, fdtd_real *out);

/*
 * Device Huygens planewave: upload tables/faces once, advance E/H faces each step.
 * field_comp 0=Ex..5=Hz; incid_nfield 1..6 (Fortran iEx..iHz); free_mode 0=fixed i,
 * 1=fixed j, 2=fixed k; a/b loop the two free absolute indices; id_axis 0=x,1=y,2=z;
 * use_e_metric 0=Idxh.., 1=Idxe..
 */
typedef struct {
   int field_comp;
   int incid_nfield;
   int wave;
   int free_mode;
   int fixed_abs;
   int a0, a1, b0, b1;
   int incid_di, incid_dj, incid_dk;
   int field_xi, field_yi, field_zi;
   int id_axis;
   int use_e_metric;
   int sign;
} fdtd_pw_face;

int fdtd_cuda_upload_planewave_phys(fdtd_cuda_ctx *ctx, int field, int axis, int base_abs, int n,
                                    const fdtd_real *host);
int fdtd_cuda_upload_planewave(fdtd_cuda_ctx *ctx, int n_waves, int max_modes, int max_numus,
                               fdtd_real cluz, const int *num_modes, const int *numus,
                               const fdtd_real *deltaevol, const fdtd_real *evol,
                               const fdtd_real *px, const fdtd_real *py, const fdtd_real *pz,
                               const fdtd_real *d0, const fdtd_real *fpw,
                               const fdtd_pw_face *faces_e, int n_faces_e,
                               const fdtd_pw_face *faces_h, int n_faces_h);
int fdtd_cuda_planewave_ready(const fdtd_cuda_ctx *ctx);
int fdtd_cuda_advance_planewave_e(fdtd_cuda_ctx *ctx, fdtd_real time, int *still_out);
int fdtd_cuda_advance_planewave_h(fdtd_cuda_ctx *ctx, fdtd_real time, int *still_out);

/*
 * Device first-order Mur ABC on magnetic faces (host AdvanceMagneticMUR).
 * Each job is one (face, H-component): update ghost plane, then copy the
 * 2-cell Past slab from current H. cab_which indexes CAB1[medio] tables
 * (0=left .. 5=front). wall_axis 0=x,1=y,2=z; neigh_sign +1 (min face) or -1.
 */
#define FDTD_CUDA_MUR_JOBS 12
#define FDTD_CUDA_MUR_CAB  6

typedef struct {
   int field_comp;     /* 3=Hx .. 5=Hz */
   int media_comp;
   int cab_which;      /* 0=left,1=right,2=down,3=up,4=back,5=front */
   int wall_axis;      /* 0=x, 1=y, 2=z */
   int neigh_sign;     /* +1 inward from min face, -1 from max face */
   int plane_abs;      /* absolute ghost index along wall_axis */
   int a0, a1, b0, b1; /* inclusive free axes (see wall_axis) */
   int e_xi, e_yi, e_zi;
   int p_xi, p_yi, p_zi;
   int nx_p, ny_p, nz_p;
   int store_xi, store_xe, store_yi, store_ye, store_zi, store_ze;
} fdtd_mur_job;

int fdtd_cuda_upload_mur_cab(fdtd_cuda_ctx *ctx, int which, int n, const fdtd_real *host);
int fdtd_cuda_set_mur_jobs(fdtd_cuda_ctx *ctx, const fdtd_mur_job *jobs, int n_jobs);
int fdtd_cuda_upload_mur_past(fdtd_cuda_ctx *ctx, int slot, const fdtd_real *host, int n_elem);
int fdtd_cuda_download_mur_past(fdtd_cuda_ctx *ctx, int slot, fdtd_real *host, int n_elem);
int fdtd_cuda_mur_ready(const fdtd_cuda_ctx *ctx);
int fdtd_cuda_advance_mur(fdtd_cuda_ctx *ctx);

/*
 * Device nodal sources (host AdvanceNodalE / AdvanceNodalH).
 * One job is one hard or soft box on one field component.
 * Samples are concatenated, 0-based, length numus+1 per job (evol(0:numus)).
 * field_comp 0=Ex..5=Hz. hard=1 assigns amp*wave; hard=0 subtracts
 * G2 or Gm2 * metric pair * amp * wave. initial_only applies only at step 0
 * and uses samples[evol_off] (evol(0)). E skips PEC media, H skips PMC.
 */
typedef struct {
   int field_comp;
   int hard;
   int initial_only;
   int xi, xe, yi, ye, zi, ze;
   int e_xi, e_yi, e_zi;
   int evol_off;
   int numus;
   fdtd_real amplitude;
   fdtd_real deltaevol;
} fdtd_nodal_job;

int fdtd_cuda_upload_nodal(fdtd_cuda_ctx *ctx,
                           const fdtd_nodal_job *jobs, int n_jobs,
                           const fdtd_real *samples, int n_samples,
                           const int *skip_e, const int *skip_h, int n_skip);
int fdtd_cuda_nodal_ready(const fdtd_cuda_ctx *ctx);
int fdtd_cuda_advance_nodal_e(fdtd_cuda_ctx *ctx, fdtd_real time, int step);
int fdtd_cuda_advance_nodal_h(fdtd_cuda_ctx *ctx, fdtd_real time, int step);

#ifdef __cplusplus
}
#endif
