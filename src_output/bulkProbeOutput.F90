module mod_bulkProbeOutput
   use FDETYPES
   use outputTypes
   use FDETYPES_TOOLS
   use mod_outputUtils
   implicit none

contains

   subroutine init_bulk_probe_output(this, lowerBound, upperBound, field, domain, outputTypeExtension, mpidir)
      type(bulk_current_probe_output_t), intent(out) :: this
      type(cell_coordinate_t), intent(in) :: lowerBound, upperBound
      integer(kind=SINGLE), intent(in) :: mpidir, field
      character(len=BUFSIZE), intent(in) :: outputTypeExtension
      type(domain_t), intent(in) :: domain

      integer(kind=SINGLE) :: i

      this%lowerBound = lowerBound
      this%upperBound = upperBound
      this%fieldComponent = field

      this%domain = domain
      this%path = get_output_path()

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

   end subroutine init_bulk_probe_output

   subroutine create_bulk_probe_output(this)
      type(bulk_current_probe_output_t), intent(inout) :: this
      character(len=BUFSIZE) :: file_time
      integer(kind=SINGLE) :: err
      err = 0

      file_time = trim(adjustl(this%path))//'_'// &
                  trim(adjustl(timeExtension))//'_'// &
                  trim(adjustl(datFileExtension))
      call create_or_clear_file(file_time, this%fileUnitTime, err)
   end subroutine create_bulk_probe_output

   subroutine update_bulk_probe_output(this, step, field)
      type(bulk_current_probe_output_t), intent(out) :: this
      real(kind=RKIND_tiempo), intent(in) :: step
      type(field_data_t), intent(in) :: field

      integer(kind=SINGLE) :: i1_m, i2_m, j1_m, j2_m, k1_m, k2_m
      integer(kind=SINGLE) :: i1, i2, j1, j2, k1, k2
      integer(kind=SINGLE) :: iii, jjj, kkk

      real(kind=RKIND), pointer, dimension(:, :, :) :: xF, yF, zF
      real(kind=RKIND), pointer, dimension(:) :: dx, dy, dz

      i1_m = this%lowerBound%x
      j1_m = this%lowerBound%y
      k1_m = this%lowerBound%z
      i2_m = this%upperBound%x
      j2_m = this%upperBound%y
      k2_m = this%upperBound%z

      i1 = i1_m
      j1 = i2_m
      k1 = j1_m
      i2 = j2_m
      j2 = k1_m
      k2 = k2_m

      xF => field%x
      yF => field%y
      zF => field%z
      dx => field%deltaX
      dy => field%deltaY
      dz => field%deltaZ

      this%serializedTimeSize = this%serializedTimeSize + 1
      this%timeStep(this%serializedTimeSize) = step
      this%valueForTime(this%serializedTimeSize) = 0.0_RKIND !Clear uninitialized value
      selectcase (this%fieldComponent)
      case (iBloqueJx)
         do JJJ = j1, j2
            this%valueForTime(this%serializedTimeSize) = &
               this%valueForTime(this%serializedTimeSize) + &
               (yF(i1_m, JJJ, k1_m - 1) - yF(i1_m, JJJ, k2_m))*dy(JJJ)
         end do
         do KKK = k1, k2
            this%valueForTime(this%serializedTimeSize) = &
               this%valueForTime(this%serializedTimeSize) + &
               (-zF(i1_m, j1_m - 1, KKK) + zF(i1_m, j2_m, KKK))*dz(KKK)
         end do

      case (iBloqueJy)
         do KKK = k1, k2
            this%valueForTime(this%serializedTimeSize) = &
               this%valueForTime(this%serializedTimeSize) + &
               (-zF(i2_m, j1_m, KKK) + zF(i1_m - 1, j1_m, KKK))*dz(KKK)
         end do
         do III = i1, i2
            this%valueForTime(this%serializedTimeSize) = &
               this%valueForTime(this%serializedTimeSize) + &
               (xF(III, j1_m, k2_m) - xF(III, j1_m, k1_m - 1))*dx(III)
         end do

      case (iBloqueJz)
         do III = i1, i2
            this%valueForTime(this%serializedTimeSize) = &
               this%valueForTime(this%serializedTimeSize) + &
               (xF(III, j1_m - 1, k1_m) - xF(III, j2_m, k1_m))*dx(III)
         end do
         do JJJ = j1, j2
            this%valueForTime(this%serializedTimeSize) = &
               this%valueForTime(this%serializedTimeSize) + &
               (-yF(i1_m - 1, JJJ, k1_m) + yF(i2_m, JJJ, k1_m))*dy(JJJ)
         end do

      case (iBloqueMx)
         do JJJ = j1, j2
            this%valueForTime(this%serializedTimeSize) = &
               this%valueForTime(this%serializedTimeSize) + &
               (-yF(i1_m, JJJ, k1_m) + yF(i1_m, JJJ, k2_m + 1))*dy(JJJ)
         end do
         do KKK = k1, k2
            this%valueForTime(this%serializedTimeSize) = &
               this%valueForTime(this%serializedTimeSize) + &
               (zF(i1_m, j1_m, KKK) - zF(i1_m, j2_m + 1, KKK))*dz(KKK)
         end do

      case (iBloqueMy)
         do KKK = k1, k2
            this%valueForTime(this%serializedTimeSize) = &
               this%valueForTime(this%serializedTimeSize) + &
               (zF(i2_m + 1, j1_m, KKK) - zF(i1_m, j1_m, KKK))*dz(KKK)
         end do
         do III = i1, i2
            this%valueForTime(this%serializedTimeSize) = &
               this%valueForTime(this%serializedTimeSize) + &
               (-xF(III, j1_m, k2_m + 1) + xF(III, j1_m, k1_m))*dx(III)
         end do

      case (iBloqueMz)
         do III = i1, i2
            this%valueForTime(this%serializedTimeSize) = &
               this%valueForTime(this%serializedTimeSize) + &
               (-xF(III, j1_m, k1_m) + xF(III, j2_m + 1, k1_m))*dx(III)
         end do
         do JJJ = j1, j2
            this%valueForTime(this%serializedTimeSize) = &
               this%valueForTime(this%serializedTimeSize) + &
               (yF(i1_m, JJJ, k1_m) - yF(i2_m + 1, JJJ, k1_m))*dy(JJJ)
         end do

      end select

   end subroutine update_bulk_probe_output

   subroutine flush_bulk_probe_output(this)
      type(bulk_current_probe_output_t), intent(inout) :: this
      character(len=BUFSIZE) :: filename
      integer :: i
      if (this%serializedTimeSize <= 0) then
         print *, "No data to write."
         return
      end if

      filename = trim(adjustl(this%path))//'_'//trim(adjustl(timeExtension))//'_'//trim(adjustl(datFileExtension))
      open (unit=this%fileUnitTime, file=filename, status="old", action="write", position="append")

      do i = 1, this%serializedTimeSize
         write (this%fileUnitTime, fmt) this%timeStep(i), this%valueForTime(i)
      end do

      close (this%fileUnitTime)
      call clear_time_data()
   contains
      subroutine clear_time_data()
         this%timeStep = 0.0_RKIND_tiempo
         this%valueForTime = 0.0_RKIND

         this%serializedTimeSize = 0
      end subroutine clear_time_data
   end subroutine flush_bulk_probe_output

end module mod_bulkProbeOutput
