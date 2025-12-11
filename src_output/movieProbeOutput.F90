module mod_movieProbeOutput
   use FDETYPES
   use outputTypes
   use mod_outputUtils
   implicit none
   private

   !===========================
   !  Public interface summary
   !===========================
   public :: init_movie_probe_output
   public :: update_movie_probe_output
   public :: flush_movie_probe_output
   !===========================

   !===========================
   !  Private interface summary
   !===========================
   private :: get_measurements_count

   !===========================

contains

   subroutine init_movie_probe_output(this, lowerBound, upperBound, field, domain, geometryMedia, registeredMedia, sinpml_fullsize, outputTypeExtension, mpidir)
      type(movie_probe_output_t), intent(inout) :: this
      type(cell_coordinate_t), intent(in) :: lowerBound, upperBound
      integer(kind=SINGLE), intent(in) :: mpidir, field
      character(len=BUFSIZE), intent(in) :: outputTypeExtension

      type(MediaData_t), pointer, dimension(:) :: registeredMedia
      type(media_matrices_t), pointer, intent(in) :: geometryMedia
      type(limit_t), pointer, dimension(:), intent(in)  :: sinpml_fullsize

      type(domain_t), intent(in) :: domain

      this%lowerBound = lowerBound
      this%upperBound = upperBound
      this%fieldComponent = field !This can refer to field or currentDensity
      this%domain = domain
      this%path = get_output_path()

      numberOfRequiredMeasures = get_measurements_count()

   contains
      function get_output_path() result(outputPath)
         character(len=BUFSIZE)  :: probeBoundsExtension, prefixFieldExtension
         character(len=BUFSIZE) :: outputPath
         probeBoundsExtension = get_coordinates_extension(this%lowerBound, this%upperBound, mpidir)
         prefixFieldExtension = get_prefix_extension(field, mpidir)
         outputPath = &
            trim(adjustl(outputTypeExtension))//'_'//trim(adjustl(prefixFieldExtension))//'_'//trim(adjustl(probeBoundsExtension))
         return
      end function get_output_path

   end subroutine init_movie_probe_output

   subroutine update_movie_probe_output()

   end subroutine update_movie_probe_output

   subroutine flush_movie_probe_output()

   end subroutine flush_movie_probe_output

   function get_measurements_count(this)
      type(movie_probe_output_t), intent(in) :: this
      integer(kind=SINGLE) :: i, j, k, field
      integer(kind=SINGLE) :: n

      

      n = 0_SINGLE
      do i = this%lowerBound%x, this%upperBound%x
      do j = this%lowerBound%y, this%upperBound%y
      do k = this%lowerBound%z, this%upperBound%z
        if 
      end do
      end do
      end do
   end function

end module mod_movieProbeOutput
