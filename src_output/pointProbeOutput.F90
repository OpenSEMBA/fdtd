module mod_pointProbeOutput
   use FDETYPES
   use outputTypes
   use mod_domain
   use mod_outputUtils

   implicit none

contains
   subroutine init_point_probe_output(this, coordinates, field, domain, outputTypeExtension, mpidir, timeInterval)
      type(point_probe_output_t), intent(out) :: this
      type(cell_coordinate_t) :: coordinates
      integer(kind=SINGLE), intent(in) :: mpidir, field
      character(len=*), intent(in) :: outputTypeExtension
      type(domain_t), intent(in) :: domain

      real(kind=RKIND_tiempo), intent(in) :: timeInterval

      integer(kind=SINGLE) :: i

      this%coordinates = coordinates

      this%fieldComponent = field

      this%domain = domain
      this%path = get_output_path()

      if (any(this%domain%domainType == (/FREQUENCY_DOMAIN, BOTH_DOMAIN/))) then
         this%nFreq = this%domain%fnum
         allocate (this%frequencySlice(this%domain%fnum))
         allocate (this%valueForFreq(this%domain%fnum))
         do i = 1, this%nFreq
            call init_frequency_slice(this%frequencySlice, this%domain)
         end do
         this%valueForFreq = (0.0_RKIND, 0.0_RKIND)

         allocate (this%auxExp_E(this%nFreq))
         allocate (this%auxExp_H(this%nFreq))
         do i = 1, this%nFreq
            this%auxExp_E(i) = timeInterval*(1.0E0_RKIND, 0.0E0_RKIND)*Exp(mcpi2*this%frequencySlice(i))   !el dt deberia ser algun tipo de promedio
            this%auxExp_H(i) = this%auxExp_E(i)*Exp(mcpi2*this%frequencySlice(i)*timeInterval*0.5_RKIND)
         end do
      end if

   contains
      function get_output_path() result(outputPath)
         character(len=BUFSIZE)  :: probeBoundsExtension, prefixFieldExtension
         character(len=BUFSIZE) :: outputPath
         probeBoundsExtension = get_coordinates_extension(this%coordinates, mpidir)
         prefixFieldExtension = get_prefix_extension(field, mpidir)
         outputPath = &
            trim(adjustl(outputTypeExtension))//'_'//trim(adjustl(prefixFieldExtension))//'_'//trim(adjustl(probeBoundsExtension))
         return
      end function get_output_path

   end subroutine init_point_probe_output

   subroutine create_point_probe_output_files(this)
      implicit none
      type(point_probe_output_t), intent(inout) :: this
      character(len=BUFSIZE) :: file_time, file_freq
      integer(kind=SINGLE) :: err
      err = 0

      file_time = trim(adjustl(this%path))//'_'// &
                  trim(adjustl(timeExtension))//'_'// &
                  trim(adjustl(datFileExtension))

      file_freq = trim(adjustl(this%path))//'_'// &
                  trim(adjustl(timeExtension))//'_'// &
                  trim(adjustl(datFileExtension))

      call create_or_clear_file(file_time, this%fileUnitTime, err)
      call create_or_clear_file(file_freq, this%fileUnitFreq, err)

   end subroutine create_point_probe_output_files

   subroutine update_point_probe_output(this, step, field)
      type(point_probe_output_t), intent(inout) :: this
      real(kind=RKIND), pointer, dimension(:, :, :), intent(in) :: field
      real(kind=RKIND_tiempo), intent(in) :: step

      integer(kind=SINGLE) :: iter

      if (any(this%domain%domainType == (/TIME_DOMAIN, BOTH_DOMAIN/))) then
         this%serializedTimeSize = this%serializedTimeSize + 1
         this%timeStep(this%serializedTimeSize) = step
         this%valueForTime(this%serializedTimeSize) = field(this%coordinates%x, this%coordinates%y, this%coordinates%z)
      end if

      if (any(this%domain%domainType == (/FREQUENCY_DOMAIN, BOTH_DOMAIN/))) then
         select case (this%fieldComponent)
         case (iEx, iEy, iEz)
            do iter = 1, this%nFreq
               this%valueForFreq(iter) = &
                  this%valueForFreq(iter) + field(this%coordinates%x, this%coordinates%y, this%coordinates%z)*(this%auxExp_E(iter)**step)
            end do
         case (iHx, iHy, iHz)
            do iter = 1, this%nFreq
               this%valueForFreq(iter) = &
                  this%valueForFreq(iter) + field(this%coordinates%x, this%coordinates%y, this%coordinates%z)*(this%auxExp_H(iter)**step)
            end do
         end select

      end if
   end subroutine update_point_probe_output

   subroutine flush_point_probe_output(this)
      type(point_probe_output_t), intent(inout) :: this
      if (any(this%domain%domainType == (/TIME_DOMAIN, BOTH_DOMAIN/))) then
         call flush_time_domain(this)
         call clear_time_data()
      end if
      if (any(this%domain%domainType == (/FREQUENCY_DOMAIN, BOTH_DOMAIN/))) then
         call flush_frequency_domain(this)
      end if
   contains

      subroutine flush_time_domain(this)
         type(point_probe_output_t), intent(in) :: this
         integer :: i
         character(len=BUFSIZE) :: filename

         if (this%serializedTimeSize <= 0) then
            print *, "No data to write."
            return
         end if

         filename = trim(adjustl(this%path))//'_'//trim(adjustl(timeExtension))//'_'//trim(adjustl(datFileExtension))
         open (unit=this%fileUnitTime, file=filename, status="old", action="write", position="append")

         do i = 1, this%serializedTimeSize
            write (this%fileUnitTime, '(F12.6,1X,F12.6)') this%timeStep(i), this%valueForTime(i)
         end do

         close (this%fileUnitTime)
      end subroutine flush_time_domain

      subroutine flush_frequency_domain(this)
         type(point_probe_output_t), intent(in) :: this
         integer ::i
         character(len=BUFSIZE) :: filename

         if (.not. allocated(this%frequencySlice) .or. .not. allocated(this%valueForFreq)) then
            print *, "Error: arrays not allocated."
            return
         end if

         if (this%nFreq <= 0) then
            print *, "No data to write."
            return
         end if
         filename = trim(adjustl(this%path))//'_'//trim(adjustl(frequencyExtension))//'_'//trim(adjustl(datFileExtension))
         open (unit=this%fileUnitFreq, file=filename, status="replace", action="write")

         do i = 1, this%nFreq
            write (this%fileUnitFreq, '(F12.6,1X,F12.6)') this%frequencySlice(i), this%valueForFreq(i)
         end do

         close (this%fileUnitFreq)
      end subroutine flush_frequency_domain

      subroutine clear_time_data()
         this%timeStep = 0.0_RKIND_tiempo
         this%valueForTime = 0.0_RKIND

         this%serializedTimeSize = 0
      end subroutine clear_time_data

   end subroutine flush_point_probe_output
end module
