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
   public :: fdtd_dims3_c, fdtd_ibox_c, fdtd_cpml_job_c

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

end module fdtd_cuda_m
