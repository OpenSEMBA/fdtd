module cuda_field_updates_m
   use, intrinsic :: iso_c_binding, only: c_int, c_ptr, c_loc
   use FDETYPES_m, only: RKIND, INTEGERSIZEOFMEDIAMATRICES
   implicit none
   private

   integer(c_int), parameter, public :: CUDA_EX=1, CUDA_EY=2, CUDA_EZ=3
   integer(c_int), parameter, public :: CUDA_HX=4, CUDA_HY=5, CUDA_HZ=6

   type, bind(C) :: cuda_field_view_t
      type(c_ptr) :: data
      integer(c_int) :: nx, ny, nz
      integer(c_int) :: lo_x, lo_y, lo_z
   end type

   type, bind(C) :: cuda_vector_view_t
      type(c_ptr) :: data
      integer(c_int) :: length
      integer(c_int) :: lo
   end type

   type, bind(C) :: cuda_sweep_t
      integer(c_int) :: nx, ny, nz
   end type

   type, bind(C), public :: cuda_cpml_job_t
      type(cuda_field_view_t) :: output, media, source, psi
      type(cuda_vector_view_t) :: coefficient_b, coefficient_c, coefficient_g
      integer(c_int) :: region_xi, region_xe, region_yi, region_ye, region_zi, region_ze
      integer(c_int) :: axis, source_direction, update_sign
   end type

   interface set_cuda_view
      module procedure set_field_view_real, set_field_view_integer, set_vector_view
   end interface

   interface
      function cuda_advance_field_c(operation, output, media, field_a, field_b, &
            metric_a, metric_b, coefficient_a, coefficient_b, sweep) &
            bind(C, name="semba_cuda_advance_field") result(status)
         import :: c_int, cuda_field_view_t, cuda_vector_view_t, cuda_sweep_t
         integer(c_int), value :: operation
         type(cuda_field_view_t), intent(in) :: output, media, field_a, field_b
         type(cuda_vector_view_t), intent(in) :: metric_a, metric_b
         type(cuda_vector_view_t), intent(in) :: coefficient_a, coefficient_b
         type(cuda_sweep_t), intent(in) :: sweep
         integer(c_int) :: status
      end function

      function cuda_advance_cpml_jobs_c(jobs, count) &
            bind(C, name="semba_cuda_advance_cpml_jobs") result(status)
         import :: c_int, cuda_cpml_job_t
         type(cuda_cpml_job_t), intent(in) :: jobs(*)
         integer(c_int), value :: count
         integer(c_int) :: status
      end function

      function cuda_advance_pml_body_c(count, psi, coefficient_b, coefficient_c, &
            difference, coefficient_g, field, update_sign) &
            bind(C, name="semba_cuda_advance_pml_body") result(status)
         import :: c_int, cuda_vector_view_t
         integer(c_int), value :: count
         type(cuda_vector_view_t), intent(in) :: psi, coefficient_b, coefficient_c
         type(cuda_vector_view_t), intent(in) :: difference, coefficient_g, field
         integer(c_int), value :: update_sign
         integer(c_int) :: status
      end function
   end interface

   public :: cuda_advance_field, cuda_advance_cpml_jobs, cuda_advance_pml_body
   public :: cuda_field_view_t, cuda_vector_view_t, set_cuda_view

contains

   subroutine cuda_advance_field(operation, output, media, field_a, field_b, &
         metric_a, metric_b, coefficient_a, coefficient_b, sweep_x, sweep_y, sweep_z)
      integer(c_int), intent(in) :: operation
      real(kind=RKIND), target, contiguous, intent(inout) :: output(:,:,:)
      integer(kind=INTEGERSIZEOFMEDIAMATRICES), target, contiguous, intent(in) :: media(:,:,:)
      real(kind=RKIND), target, contiguous, intent(in) :: field_a(:,:,:), field_b(:,:,:)
      real(kind=RKIND), target, contiguous, intent(in) :: metric_a(:), metric_b(:)
      real(kind=RKIND), target, contiguous, intent(in) :: coefficient_a(:), coefficient_b(:)
      integer, intent(in) :: sweep_x, sweep_y, sweep_z
      type(cuda_field_view_t) :: output_view, media_view, field_a_view, field_b_view
      type(cuda_vector_view_t) :: metric_a_view, metric_b_view, coefficient_a_view, coefficient_b_view
      type(cuda_sweep_t) :: sweep
      integer(c_int) :: status

      call set_cuda_view(output_view, output)
      call set_cuda_view(media_view, media)
      call set_cuda_view(field_a_view, field_a)
      call set_cuda_view(field_b_view, field_b)
      call set_cuda_view(metric_a_view, metric_a)
      call set_cuda_view(metric_b_view, metric_b)
      call set_cuda_view(coefficient_a_view, coefficient_a)
      call set_cuda_view(coefficient_b_view, coefficient_b)
      sweep%nx = int(sweep_x, c_int)
      sweep%ny = int(sweep_y, c_int)
      sweep%nz = int(sweep_z, c_int)

      status = cuda_advance_field_c(operation, output_view, media_view, field_a_view, &
            field_b_view, metric_a_view, metric_b_view, coefficient_a_view, &
            coefficient_b_view, sweep)
      if (status /= 0_c_int) then
         error stop "Native CUDA Yee update failed; see the CUDA runtime diagnostic."
      end if
   end subroutine cuda_advance_field

   subroutine cuda_advance_cpml_jobs(jobs)
      type(cuda_cpml_job_t), intent(in) :: jobs(:)
      integer(c_int) :: status

      if (size(jobs) == 0) return
      status = cuda_advance_cpml_jobs_c(jobs, int(size(jobs), c_int))
      if (status /= 0_c_int) then
         error stop "Native CUDA CPML update failed; see the CUDA runtime diagnostic."
      end if
   end subroutine cuda_advance_cpml_jobs

   subroutine cuda_advance_pml_body(psi, coefficient_b, coefficient_c, difference, &
         coefficient_g, field, update_sign)
      real(kind=RKIND), target, contiguous, intent(inout) :: psi(:), field(:)
      real(kind=RKIND), target, contiguous, intent(in) :: coefficient_b(:), coefficient_c(:)
      real(kind=RKIND), target, contiguous, intent(in) :: difference(:), coefficient_g(:)
      integer(c_int), intent(in) :: update_sign
      type(cuda_vector_view_t) :: psi_view, coefficient_b_view, coefficient_c_view
      type(cuda_vector_view_t) :: difference_view, coefficient_g_view, field_view
      integer(c_int) :: status

      call set_cuda_view(psi_view, psi)
      call set_cuda_view(coefficient_b_view, coefficient_b)
      call set_cuda_view(coefficient_c_view, coefficient_c)
      call set_cuda_view(difference_view, difference)
      call set_cuda_view(coefficient_g_view, coefficient_g)
      call set_cuda_view(field_view, field)
      status = cuda_advance_pml_body_c(int(size(psi), c_int), psi_view, &
            coefficient_b_view, coefficient_c_view, difference_view, &
            coefficient_g_view, field_view, update_sign)
      if (status /= 0_c_int) then
         error stop "Native CUDA PML-body update failed; see the CUDA runtime diagnostic."
      end if
   end subroutine cuda_advance_pml_body

   subroutine set_field_view_real(view, values, lo_x, lo_y, lo_z)
      type(cuda_field_view_t), intent(out) :: view
      real(kind=RKIND), target, contiguous, intent(in) :: values(:,:,:)
      integer, intent(in), optional :: lo_x, lo_y, lo_z

      view%data = c_loc(values(1,1,1))
      view%nx = int(size(values,1), c_int)
      view%ny = int(size(values,2), c_int)
      view%nz = int(size(values,3), c_int)
      view%lo_x = 0_c_int
      view%lo_y = 0_c_int
      view%lo_z = 0_c_int
      if (present(lo_x)) view%lo_x = int(lo_x, c_int)
      if (present(lo_y)) view%lo_y = int(lo_y, c_int)
      if (present(lo_z)) view%lo_z = int(lo_z, c_int)
   end subroutine set_field_view_real

   subroutine set_field_view_integer(view, values, lo_x, lo_y, lo_z)
      type(cuda_field_view_t), intent(out) :: view
      integer(kind=INTEGERSIZEOFMEDIAMATRICES), target, contiguous, intent(in) :: values(:,:,:)
      integer, intent(in), optional :: lo_x, lo_y, lo_z

      view%data = c_loc(values(1,1,1))
      view%nx = int(size(values,1), c_int)
      view%ny = int(size(values,2), c_int)
      view%nz = int(size(values,3), c_int)
      view%lo_x = 0_c_int
      view%lo_y = 0_c_int
      view%lo_z = 0_c_int
      if (present(lo_x)) view%lo_x = int(lo_x, c_int)
      if (present(lo_y)) view%lo_y = int(lo_y, c_int)
      if (present(lo_z)) view%lo_z = int(lo_z, c_int)
   end subroutine set_field_view_integer

   subroutine set_vector_view(view, values, lo)
      type(cuda_vector_view_t), intent(out) :: view
      real(kind=RKIND), target, contiguous, intent(in) :: values(:)
      integer, intent(in), optional :: lo

      view%data = c_loc(values(1))
      view%length = int(size(values), c_int)
      view%lo = 0_c_int
      if (present(lo)) view%lo = int(lo, c_int)
   end subroutine set_vector_view

end module cuda_field_updates_m
