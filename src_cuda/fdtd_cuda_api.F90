! Fortran ISO_C bindings for device-resident Yee + CPML (CompileWithCUDA).
module fdtd_cuda_m
   use, intrinsic :: iso_c_binding
   use FDETYPES_m, only: RKIND
   implicit none
   private

   public :: fdtd_cuda_enabled, fdtd_cuda_ctx_c
   public :: fdtd_cuda_create_f, fdtd_cuda_destroy_f, fdtd_cuda_ok_f
   public :: fdtd_cuda_alloc_fields_f, fdtd_cuda_alloc_media_f
   public :: fdtd_cuda_alloc_coeffs_f, fdtd_cuda_alloc_metrics_f
   public :: fdtd_cuda_upload_fields_f, fdtd_cuda_download_fields_f
   public :: fdtd_cuda_upload_fields_box_f, fdtd_cuda_download_fields_box_f
   public :: fdtd_cuda_upload_media_f, fdtd_cuda_upload_coeffs_f
   public :: fdtd_cuda_upload_metrics_f
   public :: fdtd_cuda_advance_ex_f, fdtd_cuda_advance_ey_f, fdtd_cuda_advance_ez_f
   public :: fdtd_cuda_advance_hx_f, fdtd_cuda_advance_hy_f, fdtd_cuda_advance_hz_f
   public :: fdtd_cuda_alloc_psi_f, fdtd_cuda_upload_psi_f, fdtd_cuda_download_psi_f
   public :: fdtd_cuda_alloc_cpml_1d_f, fdtd_cuda_cpml_apply_f
   public :: fdtd_cuda_set_point_probes_f, fdtd_cuda_gather_point_probes_f
   public :: fdtd_cuda_upload_planewave_phys_f, fdtd_cuda_upload_planewave_f
   public :: fdtd_cuda_planewave_ready_f
   public :: fdtd_cuda_advance_planewave_e_f, fdtd_cuda_advance_planewave_h_f
   public :: fdtd_cuda_upload_mur_cab_f, fdtd_cuda_set_mur_jobs_f
   public :: fdtd_cuda_upload_mur_past_f, fdtd_cuda_mur_ready_f
   public :: fdtd_cuda_advance_mur_f
   public :: fdtd_cuda_upload_nodal_f, fdtd_cuda_nodal_ready_f
   public :: fdtd_cuda_advance_nodal_e_f, fdtd_cuda_advance_nodal_h_f
   public :: fdtd_dims3_c, fdtd_ibox_c, fdtd_cpml_job_c, fdtd_pw_face_c, fdtd_mur_job_c, fdtd_nodal_job_c

   logical, save :: fdtd_cuda_enabled = .false.
   type(c_ptr), save :: fdtd_cuda_ctx_c = c_null_ptr

   type, bind(C) :: fdtd_dims3_c
      integer(c_int) :: nx, ny, nz
   end type

   type, bind(C) :: fdtd_ibox_c
      integer(c_int) :: is, ie, js, je, ks, ke
   end type

   type, bind(C) :: fdtd_cpml_job_c
      integer(c_int) :: field_comp, psi_slot, media_comp
      integer(c_int) :: h_comp_a, h_comp_b, h_diff_axis
      integer(c_int) :: p_b_which, p_c_which, p_base, free_axis
      integer(c_int) :: xi, xe, yi, ye, zi, ze
      integer(c_int) :: e_xi, e_yi, e_zi
      integer(c_int) :: psi_xi, psi_yi, psi_zi
      integer(c_int) :: nx_psi, ny_psi, nz_psi
      integer(c_int) :: h_sign, use_fixed_medio, medio_fixed
   end type

   type, bind(C) :: fdtd_pw_face_c
      integer(c_int) :: field_comp
      integer(c_int) :: incid_nfield
      integer(c_int) :: wave
      integer(c_int) :: free_mode
      integer(c_int) :: fixed_abs
      integer(c_int) :: a0, a1, b0, b1
      integer(c_int) :: incid_di, incid_dj, incid_dk
      integer(c_int) :: field_xi, field_yi, field_zi
      integer(c_int) :: id_axis
      integer(c_int) :: use_e_metric
      integer(c_int) :: sign
   end type

   type, bind(C) :: fdtd_mur_job_c
      integer(c_int) :: field_comp
      integer(c_int) :: media_comp
      integer(c_int) :: cab_which
      integer(c_int) :: wall_axis
      integer(c_int) :: neigh_sign
      integer(c_int) :: plane_abs
      integer(c_int) :: a0, a1, b0, b1
      integer(c_int) :: e_xi, e_yi, e_zi
      integer(c_int) :: p_xi, p_yi, p_zi
      integer(c_int) :: nx_p, ny_p, nz_p
      integer(c_int) :: store_xi, store_xe, store_yi, store_ye, store_zi, store_ze
   end type

   type, bind(C) :: fdtd_nodal_job_c
      integer(c_int) :: field_comp
      integer(c_int) :: hard
      integer(c_int) :: initial_only
      integer(c_int) :: xi, xe, yi, ye, zi, ze
      integer(c_int) :: e_xi, e_yi, e_zi
      integer(c_int) :: evol_off
      integer(c_int) :: numus
      real(c_float) :: amplitude
      real(c_float) :: deltaevol
   end type

   interface
      function fdtd_cuda_create() bind(C, name="fdtd_cuda_create")
         import :: c_ptr
         type(c_ptr) :: fdtd_cuda_create
      end function
      subroutine fdtd_cuda_destroy(ctx) bind(C, name="fdtd_cuda_destroy")
         import :: c_ptr
         type(c_ptr), value :: ctx
      end subroutine
      function fdtd_cuda_ok(ctx) bind(C, name="fdtd_cuda_ok")
         import :: c_ptr, c_int
         type(c_ptr), value :: ctx
         integer(c_int) :: fdtd_cuda_ok
      end function
      function fdtd_cuda_alloc_fields(ctx, ex, ey, ez, hx, hy, hz) &
         bind(C, name="fdtd_cuda_alloc_fields")
         import :: c_ptr, c_int, fdtd_dims3_c
         type(c_ptr), value :: ctx
         type(fdtd_dims3_c), value :: ex, ey, ez, hx, hy, hz
         integer(c_int) :: fdtd_cuda_alloc_fields
      end function
      function fdtd_cuda_alloc_media(ctx, mex, mey, mez, mhx, mhy, mhz) &
         bind(C, name="fdtd_cuda_alloc_media")
         import :: c_ptr, c_int, fdtd_dims3_c
         type(c_ptr), value :: ctx
         type(fdtd_dims3_c), value :: mex, mey, mez, mhx, mhy, mhz
         integer(c_int) :: fdtd_cuda_alloc_media
      end function
      function fdtd_cuda_alloc_coeffs(ctx, num_media) bind(C, name="fdtd_cuda_alloc_coeffs")
         import :: c_ptr, c_int
         type(c_ptr), value :: ctx
         integer(c_int), value :: num_media
         integer(c_int) :: fdtd_cuda_alloc_coeffs
      end function
      function fdtd_cuda_alloc_metrics(ctx, n_idxh, n_idyh, n_idzh, n_idxe, n_idye, n_idze) &
         bind(C, name="fdtd_cuda_alloc_metrics")
         import :: c_ptr, c_int
         type(c_ptr), value :: ctx
         integer(c_int), value :: n_idxh, n_idyh, n_idzh, n_idxe, n_idye, n_idze
         integer(c_int) :: fdtd_cuda_alloc_metrics
      end function
      function fdtd_cuda_upload_fields(ctx, Ex, Ey, Ez, Hx, Hy, Hz) &
         bind(C, name="fdtd_cuda_upload_fields")
         import :: c_ptr, c_int, c_float
         type(c_ptr), value :: ctx
         real(c_float), intent(in) :: Ex(*), Ey(*), Ez(*), Hx(*), Hy(*), Hz(*)
         integer(c_int) :: fdtd_cuda_upload_fields
      end function
      function fdtd_cuda_download_fields(ctx, Ex, Ey, Ez, Hx, Hy, Hz) &
         bind(C, name="fdtd_cuda_download_fields")
         import :: c_ptr, c_int, c_float
         type(c_ptr), value :: ctx
         real(c_float), intent(out) :: Ex(*), Ey(*), Ez(*), Hx(*), Hy(*), Hz(*)
         integer(c_int) :: fdtd_cuda_download_fields
      end function
      function fdtd_cuda_upload_fields_box(ctx, Ex, Ey, Ez, Hx, Hy, Hz, &
                                           bex, bey, bez, bhx, bhy, bhz) &
         bind(C, name="fdtd_cuda_upload_fields_box")
         import :: c_ptr, c_int, c_float, fdtd_ibox_c
         type(c_ptr), value :: ctx
         real(c_float), intent(in) :: Ex(*), Ey(*), Ez(*), Hx(*), Hy(*), Hz(*)
         type(fdtd_ibox_c), intent(in) :: bex, bey, bez, bhx, bhy, bhz
         integer(c_int) :: fdtd_cuda_upload_fields_box
      end function
      function fdtd_cuda_download_fields_box(ctx, Ex, Ey, Ez, Hx, Hy, Hz, &
                                             bex, bey, bez, bhx, bhy, bhz) &
         bind(C, name="fdtd_cuda_download_fields_box")
         import :: c_ptr, c_int, c_float, fdtd_ibox_c
         type(c_ptr), value :: ctx
         real(c_float), intent(out) :: Ex(*), Ey(*), Ez(*), Hx(*), Hy(*), Hz(*)
         type(fdtd_ibox_c), intent(in) :: bex, bey, bez, bhx, bhy, bhz
         integer(c_int) :: fdtd_cuda_download_fields_box
      end function
      function fdtd_cuda_upload_media(ctx, mex, mey, mez, mhx, mhy, mhz) &
         bind(C, name="fdtd_cuda_upload_media")
         import :: c_ptr, c_int, c_short
         type(c_ptr), value :: ctx
         integer(c_short), intent(in) :: mex(*), mey(*), mez(*), mhx(*), mhy(*), mhz(*)
         integer(c_int) :: fdtd_cuda_upload_media
      end function
      function fdtd_cuda_upload_coeffs(ctx, g1, g2, gm1, gm2, num_media) &
         bind(C, name="fdtd_cuda_upload_coeffs")
         import :: c_ptr, c_int, c_float
         type(c_ptr), value :: ctx
         real(c_float), intent(in) :: g1(*), g2(*), gm1(*), gm2(*)
         integer(c_int), value :: num_media
         integer(c_int) :: fdtd_cuda_upload_coeffs
      end function
      function fdtd_cuda_upload_metrics(ctx, Idxh, Idyh, Idzh, Idxe, Idye, Idze) &
         bind(C, name="fdtd_cuda_upload_metrics")
         import :: c_ptr, c_int, c_float
         type(c_ptr), value :: ctx
         real(c_float), intent(in) :: Idxh(*), Idyh(*), Idzh(*), Idxe(*), Idye(*), Idze(*)
         integer(c_int) :: fdtd_cuda_upload_metrics
      end function
      function fdtd_cuda_advance_ex(ctx, sweep) bind(C, name="fdtd_cuda_advance_ex")
         import :: c_ptr, c_int, fdtd_ibox_c
         type(c_ptr), value :: ctx
         type(fdtd_ibox_c), value :: sweep
         integer(c_int) :: fdtd_cuda_advance_ex
      end function
      function fdtd_cuda_advance_ey(ctx, sweep) bind(C, name="fdtd_cuda_advance_ey")
         import :: c_ptr, c_int, fdtd_ibox_c
         type(c_ptr), value :: ctx
         type(fdtd_ibox_c), value :: sweep
         integer(c_int) :: fdtd_cuda_advance_ey
      end function
      function fdtd_cuda_advance_ez(ctx, sweep) bind(C, name="fdtd_cuda_advance_ez")
         import :: c_ptr, c_int, fdtd_ibox_c
         type(c_ptr), value :: ctx
         type(fdtd_ibox_c), value :: sweep
         integer(c_int) :: fdtd_cuda_advance_ez
      end function
      function fdtd_cuda_advance_hx(ctx, sweep) bind(C, name="fdtd_cuda_advance_hx")
         import :: c_ptr, c_int, fdtd_ibox_c
         type(c_ptr), value :: ctx
         type(fdtd_ibox_c), value :: sweep
         integer(c_int) :: fdtd_cuda_advance_hx
      end function
      function fdtd_cuda_advance_hy(ctx, sweep) bind(C, name="fdtd_cuda_advance_hy")
         import :: c_ptr, c_int, fdtd_ibox_c
         type(c_ptr), value :: ctx
         type(fdtd_ibox_c), value :: sweep
         integer(c_int) :: fdtd_cuda_advance_hy
      end function
      function fdtd_cuda_advance_hz(ctx, sweep) bind(C, name="fdtd_cuda_advance_hz")
         import :: c_ptr, c_int, fdtd_ibox_c
         type(c_ptr), value :: ctx
         type(fdtd_ibox_c), value :: sweep
         integer(c_int) :: fdtd_cuda_advance_hz
      end function
      function fdtd_cuda_alloc_psi(ctx, slot, n_elem) bind(C, name="fdtd_cuda_alloc_psi")
         import :: c_ptr, c_int
         type(c_ptr), value :: ctx
         integer(c_int), value :: slot, n_elem
         integer(c_int) :: fdtd_cuda_alloc_psi
      end function
      function fdtd_cuda_upload_psi(ctx, slot, psi, n_elem) bind(C, name="fdtd_cuda_upload_psi")
         import :: c_ptr, c_int, c_float
         type(c_ptr), value :: ctx
         integer(c_int), value :: slot, n_elem
         real(c_float), intent(in) :: psi(*)
         integer(c_int) :: fdtd_cuda_upload_psi
      end function
      function fdtd_cuda_download_psi(ctx, slot, psi, n_elem) bind(C, name="fdtd_cuda_download_psi")
         import :: c_ptr, c_int, c_float
         type(c_ptr), value :: ctx
         integer(c_int), value :: slot, n_elem
         real(c_float), intent(out) :: psi(*)
         integer(c_int) :: fdtd_cuda_download_psi
      end function
      function fdtd_cuda_alloc_cpml_1d(ctx, which, n, host) bind(C, name="fdtd_cuda_alloc_cpml_1d")
         import :: c_ptr, c_int, c_float
         type(c_ptr), value :: ctx
         integer(c_int), value :: which, n
         real(c_float), intent(in) :: host(*)
         integer(c_int) :: fdtd_cuda_alloc_cpml_1d
      end function
      function fdtd_cuda_cpml_apply(ctx, job) bind(C, name="fdtd_cuda_cpml_apply")
         import :: c_ptr, c_int, fdtd_cpml_job_c
         type(c_ptr), value :: ctx
         type(fdtd_cpml_job_c), intent(in) :: job
         integer(c_int) :: fdtd_cuda_cpml_apply
      end function
      function fdtd_cuda_set_point_probes(ctx, n, comp, i, j, k) &
         bind(C, name="fdtd_cuda_set_point_probes")
         import :: c_ptr, c_int
         type(c_ptr), value :: ctx
         integer(c_int), value :: n
         integer(c_int), intent(in) :: comp(*), i(*), j(*), k(*)
         integer(c_int) :: fdtd_cuda_set_point_probes
      end function
      function fdtd_cuda_gather_point_probes(ctx, out) bind(C, name="fdtd_cuda_gather_point_probes")
         import :: c_ptr, c_int, c_float
         type(c_ptr), value :: ctx
         real(c_float), intent(out) :: out(*)
         integer(c_int) :: fdtd_cuda_gather_point_probes
      end function
      function fdtd_cuda_upload_planewave_phys(ctx, field, axis, base_abs, n, host) &
         bind(C, name="fdtd_cuda_upload_planewave_phys")
         import :: c_ptr, c_int, c_float
         type(c_ptr), value :: ctx
         integer(c_int), value :: field, axis, base_abs, n
         real(c_float), intent(in) :: host(*)
         integer(c_int) :: fdtd_cuda_upload_planewave_phys
      end function
      function fdtd_cuda_upload_planewave(ctx, n_waves, max_modes, max_numus, cluz, &
         num_modes, numus, deltaevol, evol, px, py, pz, d0, fpw, &
         faces_e, n_faces_e, faces_h, n_faces_h) bind(C, name="fdtd_cuda_upload_planewave")
         import :: c_ptr, c_int, c_float, fdtd_pw_face_c
         type(c_ptr), value :: ctx
         integer(c_int), value :: n_waves, max_modes, max_numus, n_faces_e, n_faces_h
         real(c_float), value :: cluz
         integer(c_int), intent(in) :: num_modes(*), numus(*)
         real(c_float), intent(in) :: deltaevol(*), evol(*), px(*), py(*), pz(*), d0(*), fpw(*)
         type(fdtd_pw_face_c), intent(in) :: faces_e(*), faces_h(*)
         integer(c_int) :: fdtd_cuda_upload_planewave
      end function
      function fdtd_cuda_planewave_ready(ctx) bind(C, name="fdtd_cuda_planewave_ready")
         import :: c_ptr, c_int
         type(c_ptr), value :: ctx
         integer(c_int) :: fdtd_cuda_planewave_ready
      end function
      function fdtd_cuda_advance_planewave_e(ctx, time, still_out) &
         bind(C, name="fdtd_cuda_advance_planewave_e")
         import :: c_ptr, c_int, c_float
         type(c_ptr), value :: ctx
         real(c_float), value :: time
         integer(c_int), intent(out) :: still_out
         integer(c_int) :: fdtd_cuda_advance_planewave_e
      end function
      function fdtd_cuda_advance_planewave_h(ctx, time, still_out) &
         bind(C, name="fdtd_cuda_advance_planewave_h")
         import :: c_ptr, c_int, c_float
         type(c_ptr), value :: ctx
         real(c_float), value :: time
         integer(c_int), intent(out) :: still_out
         integer(c_int) :: fdtd_cuda_advance_planewave_h
      end function
      function fdtd_cuda_upload_mur_cab(ctx, which, n, host) bind(C, name="fdtd_cuda_upload_mur_cab")
         import :: c_ptr, c_int, c_float
         type(c_ptr), value :: ctx
         integer(c_int), value :: which, n
         real(c_float), intent(in) :: host(*)
         integer(c_int) :: fdtd_cuda_upload_mur_cab
      end function
      function fdtd_cuda_set_mur_jobs(ctx, jobs, n_jobs) bind(C, name="fdtd_cuda_set_mur_jobs")
         import :: c_ptr, c_int, fdtd_mur_job_c
         type(c_ptr), value :: ctx
         type(fdtd_mur_job_c), intent(in) :: jobs(*)
         integer(c_int), value :: n_jobs
         integer(c_int) :: fdtd_cuda_set_mur_jobs
      end function
      function fdtd_cuda_upload_mur_past(ctx, slot, host, n_elem) &
         bind(C, name="fdtd_cuda_upload_mur_past")
         import :: c_ptr, c_int, c_float
         type(c_ptr), value :: ctx
         integer(c_int), value :: slot, n_elem
         real(c_float), intent(in) :: host(*)
         integer(c_int) :: fdtd_cuda_upload_mur_past
      end function
      function fdtd_cuda_mur_ready(ctx) bind(C, name="fdtd_cuda_mur_ready")
         import :: c_ptr, c_int
         type(c_ptr), value :: ctx
         integer(c_int) :: fdtd_cuda_mur_ready
      end function
      function fdtd_cuda_advance_mur(ctx) bind(C, name="fdtd_cuda_advance_mur")
         import :: c_ptr, c_int
         type(c_ptr), value :: ctx
         integer(c_int) :: fdtd_cuda_advance_mur
      end function
      function fdtd_cuda_upload_nodal(ctx, jobs, n_jobs, samples, n_samples, skip_e, skip_h, n_skip) &
         bind(C, name="fdtd_cuda_upload_nodal")
         import :: c_ptr, c_int, c_float, fdtd_nodal_job_c
         type(c_ptr), value :: ctx
         type(fdtd_nodal_job_c), intent(in) :: jobs(*)
         integer(c_int), value :: n_jobs, n_samples, n_skip
         real(c_float), intent(in) :: samples(*)
         integer(c_int), intent(in) :: skip_e(*), skip_h(*)
         integer(c_int) :: fdtd_cuda_upload_nodal
      end function
      function fdtd_cuda_nodal_ready(ctx) bind(C, name="fdtd_cuda_nodal_ready")
         import :: c_ptr, c_int
         type(c_ptr), value :: ctx
         integer(c_int) :: fdtd_cuda_nodal_ready
      end function
      function fdtd_cuda_advance_nodal_e(ctx, time, step) bind(C, name="fdtd_cuda_advance_nodal_e")
         import :: c_ptr, c_int, c_float
         type(c_ptr), value :: ctx
         real(c_float), value :: time
         integer(c_int), value :: step
         integer(c_int) :: fdtd_cuda_advance_nodal_e
      end function
      function fdtd_cuda_advance_nodal_h(ctx, time, step) bind(C, name="fdtd_cuda_advance_nodal_h")
         import :: c_ptr, c_int, c_float
         type(c_ptr), value :: ctx
         real(c_float), value :: time
         integer(c_int), value :: step
         integer(c_int) :: fdtd_cuda_advance_nodal_h
      end function
   end interface

contains

   subroutine fdtd_cuda_create_f(ok)
      logical, intent(out) :: ok
      fdtd_cuda_ctx_c = fdtd_cuda_create()
      ok = (fdtd_cuda_ok(fdtd_cuda_ctx_c) /= 0)
      fdtd_cuda_enabled = ok
   end subroutine

   subroutine fdtd_cuda_destroy_f()
      if (c_associated(fdtd_cuda_ctx_c)) call fdtd_cuda_destroy(fdtd_cuda_ctx_c)
      fdtd_cuda_ctx_c = c_null_ptr
      fdtd_cuda_enabled = .false.
   end subroutine

   logical function fdtd_cuda_ok_f()
      fdtd_cuda_ok_f = fdtd_cuda_enabled .and. (fdtd_cuda_ok(fdtd_cuda_ctx_c) /= 0)
   end function

   integer function fdtd_cuda_alloc_fields_f(ex, ey, ez, hx, hy, hz)
      type(fdtd_dims3_c), intent(in) :: ex, ey, ez, hx, hy, hz
      fdtd_cuda_alloc_fields_f = fdtd_cuda_alloc_fields(fdtd_cuda_ctx_c, ex, ey, ez, hx, hy, hz)
   end function

   integer function fdtd_cuda_alloc_media_f(mex, mey, mez, mhx, mhy, mhz)
      type(fdtd_dims3_c), intent(in) :: mex, mey, mez, mhx, mhy, mhz
      fdtd_cuda_alloc_media_f = fdtd_cuda_alloc_media(fdtd_cuda_ctx_c, mex, mey, mez, mhx, mhy, mhz)
   end function

   integer function fdtd_cuda_alloc_coeffs_f(num_media)
      integer, intent(in) :: num_media
      fdtd_cuda_alloc_coeffs_f = fdtd_cuda_alloc_coeffs(fdtd_cuda_ctx_c, int(num_media, c_int))
   end function

   integer function fdtd_cuda_alloc_metrics_f(n_idxh, n_idyh, n_idzh, n_idxe, n_idye, n_idze)
      integer, intent(in) :: n_idxh, n_idyh, n_idzh, n_idxe, n_idye, n_idze
      fdtd_cuda_alloc_metrics_f = fdtd_cuda_alloc_metrics(fdtd_cuda_ctx_c, &
         int(n_idxh, c_int), int(n_idyh, c_int), int(n_idzh, c_int), &
         int(n_idxe, c_int), int(n_idye, c_int), int(n_idze, c_int))
   end function

   integer function fdtd_cuda_upload_fields_f(Ex, Ey, Ez, Hx, Hy, Hz)
      real(kind=RKIND), intent(in), target :: Ex(*), Ey(*), Ez(*), Hx(*), Hy(*), Hz(*)
      fdtd_cuda_upload_fields_f = fdtd_cuda_upload_fields(fdtd_cuda_ctx_c, Ex, Ey, Ez, Hx, Hy, Hz)
   end function

   integer function fdtd_cuda_download_fields_f(Ex, Ey, Ez, Hx, Hy, Hz)
      real(kind=RKIND), intent(out), target :: Ex(*), Ey(*), Ez(*), Hx(*), Hy(*), Hz(*)
      fdtd_cuda_download_fields_f = fdtd_cuda_download_fields(fdtd_cuda_ctx_c, Ex, Ey, Ez, Hx, Hy, Hz)
   end function

   integer function fdtd_cuda_upload_fields_box_f(Ex, Ey, Ez, Hx, Hy, Hz, &
                                                  bex, bey, bez, bhx, bhy, bhz)
      real(kind=RKIND), intent(in), target :: Ex(*), Ey(*), Ez(*), Hx(*), Hy(*), Hz(*)
      type(fdtd_ibox_c), intent(in) :: bex, bey, bez, bhx, bhy, bhz
      fdtd_cuda_upload_fields_box_f = fdtd_cuda_upload_fields_box(fdtd_cuda_ctx_c, &
         Ex, Ey, Ez, Hx, Hy, Hz, bex, bey, bez, bhx, bhy, bhz)
   end function

   integer function fdtd_cuda_download_fields_box_f(Ex, Ey, Ez, Hx, Hy, Hz, &
                                                     bex, bey, bez, bhx, bhy, bhz)
      real(kind=RKIND), intent(out), target :: Ex(*), Ey(*), Ez(*), Hx(*), Hy(*), Hz(*)
      type(fdtd_ibox_c), intent(in) :: bex, bey, bez, bhx, bhy, bhz
      fdtd_cuda_download_fields_box_f = fdtd_cuda_download_fields_box(fdtd_cuda_ctx_c, &
         Ex, Ey, Ez, Hx, Hy, Hz, bex, bey, bez, bhx, bhy, bhz)
   end function

   integer function fdtd_cuda_upload_media_f(mex, mey, mez, mhx, mhy, mhz)
      integer(kind=2), intent(in), target :: mex(*), mey(*), mez(*), mhx(*), mhy(*), mhz(*)
      fdtd_cuda_upload_media_f = fdtd_cuda_upload_media(fdtd_cuda_ctx_c, mex, mey, mez, mhx, mhy, mhz)
   end function

   integer function fdtd_cuda_upload_coeffs_f(g1, g2, gm1, gm2, num_media)
      real(kind=RKIND), intent(in), target :: g1(*), g2(*), gm1(*), gm2(*)
      integer, intent(in) :: num_media
      fdtd_cuda_upload_coeffs_f = fdtd_cuda_upload_coeffs(fdtd_cuda_ctx_c, g1, g2, gm1, gm2, int(num_media, c_int))
   end function

   integer function fdtd_cuda_upload_metrics_f(Idxh, Idyh, Idzh, Idxe, Idye, Idze)
      real(kind=RKIND), intent(in), target :: Idxh(*), Idyh(*), Idzh(*), Idxe(*), Idye(*), Idze(*)
      fdtd_cuda_upload_metrics_f = fdtd_cuda_upload_metrics(fdtd_cuda_ctx_c, Idxh, Idyh, Idzh, Idxe, Idye, Idze)
   end function

   integer function fdtd_cuda_advance_ex_f(sweep)
      type(fdtd_ibox_c), intent(in) :: sweep
      fdtd_cuda_advance_ex_f = fdtd_cuda_advance_ex(fdtd_cuda_ctx_c, sweep)
   end function
   integer function fdtd_cuda_advance_ey_f(sweep)
      type(fdtd_ibox_c), intent(in) :: sweep
      fdtd_cuda_advance_ey_f = fdtd_cuda_advance_ey(fdtd_cuda_ctx_c, sweep)
   end function
   integer function fdtd_cuda_advance_ez_f(sweep)
      type(fdtd_ibox_c), intent(in) :: sweep
      fdtd_cuda_advance_ez_f = fdtd_cuda_advance_ez(fdtd_cuda_ctx_c, sweep)
   end function
   integer function fdtd_cuda_advance_hx_f(sweep)
      type(fdtd_ibox_c), intent(in) :: sweep
      fdtd_cuda_advance_hx_f = fdtd_cuda_advance_hx(fdtd_cuda_ctx_c, sweep)
   end function
   integer function fdtd_cuda_advance_hy_f(sweep)
      type(fdtd_ibox_c), intent(in) :: sweep
      fdtd_cuda_advance_hy_f = fdtd_cuda_advance_hy(fdtd_cuda_ctx_c, sweep)
   end function
   integer function fdtd_cuda_advance_hz_f(sweep)
      type(fdtd_ibox_c), intent(in) :: sweep
      fdtd_cuda_advance_hz_f = fdtd_cuda_advance_hz(fdtd_cuda_ctx_c, sweep)
   end function

   integer function fdtd_cuda_alloc_psi_f(slot, n_elem)
      integer, intent(in) :: slot, n_elem
      fdtd_cuda_alloc_psi_f = fdtd_cuda_alloc_psi(fdtd_cuda_ctx_c, int(slot, c_int), int(n_elem, c_int))
   end function
   integer function fdtd_cuda_upload_psi_f(slot, psi, n_elem)
      integer, intent(in) :: slot, n_elem
      real(kind=RKIND), intent(in), target :: psi(*)
      fdtd_cuda_upload_psi_f = fdtd_cuda_upload_psi(fdtd_cuda_ctx_c, int(slot, c_int), psi, int(n_elem, c_int))
   end function
   integer function fdtd_cuda_download_psi_f(slot, psi, n_elem)
      integer, intent(in) :: slot, n_elem
      real(kind=RKIND), intent(out), target :: psi(*)
      fdtd_cuda_download_psi_f = fdtd_cuda_download_psi(fdtd_cuda_ctx_c, int(slot, c_int), psi, int(n_elem, c_int))
   end function
   integer function fdtd_cuda_alloc_cpml_1d_f(which, n, host)
      integer, intent(in) :: which, n
      real(kind=RKIND), intent(in), target :: host(*)
      fdtd_cuda_alloc_cpml_1d_f = fdtd_cuda_alloc_cpml_1d(fdtd_cuda_ctx_c, int(which, c_int), int(n, c_int), host)
   end function
   integer function fdtd_cuda_cpml_apply_f(job)
      type(fdtd_cpml_job_c), intent(in) :: job
      fdtd_cuda_cpml_apply_f = fdtd_cuda_cpml_apply(fdtd_cuda_ctx_c, job)
   end function

   integer function fdtd_cuda_set_point_probes_f(n, comp, i, j, k)
      integer, intent(in) :: n
      integer(c_int), intent(in), target :: comp(*), i(*), j(*), k(*)
      fdtd_cuda_set_point_probes_f = fdtd_cuda_set_point_probes(fdtd_cuda_ctx_c, int(n, c_int), &
         comp, i, j, k)
   end function

   integer function fdtd_cuda_gather_point_probes_f(out)
      real(kind=RKIND), intent(out), target :: out(*)
      fdtd_cuda_gather_point_probes_f = fdtd_cuda_gather_point_probes(fdtd_cuda_ctx_c, out)
   end function

   integer function fdtd_cuda_upload_planewave_phys_f(field, axis, base_abs, n, host)
      integer, intent(in) :: field, axis, base_abs, n
      real(kind=RKIND), intent(in), target :: host(*)
      fdtd_cuda_upload_planewave_phys_f = fdtd_cuda_upload_planewave_phys(fdtd_cuda_ctx_c, &
         int(field, c_int), int(axis, c_int), int(base_abs, c_int), int(n, c_int), host)
   end function

   integer function fdtd_cuda_upload_planewave_f(n_waves, max_modes, max_numus, cluz, &
      num_modes, numus, deltaevol, evol, px, py, pz, d0, fpw, faces_e, n_faces_e, faces_h, n_faces_h)
      integer, intent(in) :: n_waves, max_modes, max_numus, n_faces_e, n_faces_h
      real(kind=RKIND), intent(in) :: cluz
      integer(c_int), intent(in), target :: num_modes(*), numus(*)
      real(kind=RKIND), intent(in), target :: deltaevol(*), evol(*), px(*), py(*), pz(*), d0(*), fpw(*)
      type(fdtd_pw_face_c), intent(in), target :: faces_e(*), faces_h(*)
      fdtd_cuda_upload_planewave_f = fdtd_cuda_upload_planewave(fdtd_cuda_ctx_c, &
         int(n_waves, c_int), int(max_modes, c_int), int(max_numus, c_int), real(cluz, c_float), &
         num_modes, numus, deltaevol, evol, px, py, pz, d0, fpw, &
         faces_e, int(n_faces_e, c_int), faces_h, int(n_faces_h, c_int))
   end function

   logical function fdtd_cuda_planewave_ready_f()
      fdtd_cuda_planewave_ready_f = fdtd_cuda_enabled .and. (fdtd_cuda_planewave_ready(fdtd_cuda_ctx_c) /= 0)
   end function

   integer function fdtd_cuda_advance_planewave_e_f(time, still_out)
      real(kind=RKIND), intent(in) :: time
      logical, intent(out) :: still_out
      integer(c_int) :: still_c
      fdtd_cuda_advance_planewave_e_f = fdtd_cuda_advance_planewave_e(fdtd_cuda_ctx_c, &
         real(time, c_float), still_c)
      still_out = (still_c /= 0)
   end function

   integer function fdtd_cuda_advance_planewave_h_f(time, still_out)
      real(kind=RKIND), intent(in) :: time
      logical, intent(out) :: still_out
      integer(c_int) :: still_c
      fdtd_cuda_advance_planewave_h_f = fdtd_cuda_advance_planewave_h(fdtd_cuda_ctx_c, &
         real(time, c_float), still_c)
      still_out = (still_c /= 0)
   end function

   integer function fdtd_cuda_upload_mur_cab_f(which, n, host)
      integer, intent(in) :: which, n
      real(kind=RKIND), intent(in), target :: host(*)
      fdtd_cuda_upload_mur_cab_f = fdtd_cuda_upload_mur_cab(fdtd_cuda_ctx_c, &
         int(which, c_int), int(n, c_int), host)
   end function

   integer function fdtd_cuda_set_mur_jobs_f(jobs, n_jobs)
      type(fdtd_mur_job_c), intent(in), target :: jobs(*)
      integer, intent(in) :: n_jobs
      fdtd_cuda_set_mur_jobs_f = fdtd_cuda_set_mur_jobs(fdtd_cuda_ctx_c, jobs, int(n_jobs, c_int))
   end function

   integer function fdtd_cuda_upload_mur_past_f(slot, host, n_elem)
      integer, intent(in) :: slot, n_elem
      real(kind=RKIND), intent(in), target :: host(*)
      fdtd_cuda_upload_mur_past_f = fdtd_cuda_upload_mur_past(fdtd_cuda_ctx_c, &
         int(slot, c_int), host, int(n_elem, c_int))
   end function

   logical function fdtd_cuda_mur_ready_f()
      fdtd_cuda_mur_ready_f = fdtd_cuda_enabled .and. (fdtd_cuda_mur_ready(fdtd_cuda_ctx_c) /= 0)
   end function

   integer function fdtd_cuda_advance_mur_f()
      fdtd_cuda_advance_mur_f = fdtd_cuda_advance_mur(fdtd_cuda_ctx_c)
   end function

   integer function fdtd_cuda_upload_nodal_f(jobs, n_jobs, samples, n_samples, skip_e, skip_h, n_skip)
      type(fdtd_nodal_job_c), intent(in), target :: jobs(*)
      integer, intent(in) :: n_jobs, n_samples, n_skip
      real(kind=RKIND), intent(in), target :: samples(*)
      integer(c_int), intent(in), target :: skip_e(*), skip_h(*)
      fdtd_cuda_upload_nodal_f = fdtd_cuda_upload_nodal(fdtd_cuda_ctx_c, jobs, &
         int(n_jobs, c_int), samples, int(n_samples, c_int), skip_e, skip_h, int(n_skip, c_int))
   end function

   logical function fdtd_cuda_nodal_ready_f()
      fdtd_cuda_nodal_ready_f = fdtd_cuda_enabled .and. (fdtd_cuda_nodal_ready(fdtd_cuda_ctx_c) /= 0)
   end function

   integer function fdtd_cuda_advance_nodal_e_f(time, step)
      real(kind=RKIND), intent(in) :: time
      integer, intent(in) :: step
      fdtd_cuda_advance_nodal_e_f = fdtd_cuda_advance_nodal_e(fdtd_cuda_ctx_c, &
         real(time, c_float), int(step, c_int))
   end function

   integer function fdtd_cuda_advance_nodal_h_f(time, step)
      real(kind=RKIND), intent(in) :: time
      integer, intent(in) :: step
      fdtd_cuda_advance_nodal_h_f = fdtd_cuda_advance_nodal_h(fdtd_cuda_ctx_c, &
         real(time, c_float), int(step, c_int))
   end function

end module fdtd_cuda_m
