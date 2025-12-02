module mod_bulkProbe
   use FDETYPES
   use FDETYPES_TOOLS
   use mod_domain
   use mod_outputUtils
   implicit none

   type bulk_probe_output_t
      integer(kind=SINGLE) :: columnas = 2_SINGLE !reference and field
      type(domain_t) :: domain
      integer(kind=SINGLE) :: xCoord, yCoord, zCoord
      integer(kind=SINGLE) :: x2Coord, y2Coord, z2Coord
      character(len=BUFSIZE) :: path
      integer(kind=SINGLE) :: fieldComponent
      integer(kind=SINGLE) :: serializedTimeSize = 0_SINGLE
      real(kind=RKIND_tiempo), dimension(BuffObse) :: timeStep = 0.0_RKIND
      real(kind=RKIND), dimension(BuffObse) :: valueForTime = 0.0_RKIND

   end type bulk_probe_output_t

contains

   subroutine init_bulk_probe_output(this, iCoord, jCoord, kCoord, i2Coord, j2Coord, k2Coord, field, domain, outputTypeExtension, mpidir)
      type(bulk_probe_output_t), intent(out) :: this
      integer(kind=SINGLE), intent(in) :: iCoord, jCoord, kCoord
      integer(kind=SINGLE), intent(in) :: i2Coord, j2Coord, k2Coord
      integer(kind=SINGLE), intent(in) :: mpidir, field
      character(len=BUFSIZE), intent(in) :: outputTypeExtension
      type(domain_t), intent(in) :: domain

      integer(kind=SINGLE) :: i

      this%xCoord = iCoord
      this%yCoord = jCoord
      this%zCoord = kCoord

      this%x2Coord = i2Coord
      this%y2Coord = j2Coord
      this%z2Coord = k2Coord

      this%fieldComponent = field

      this%domain = domain
      this%path = get_output_path()

   contains

      function get_output_path() result(outputPath)
         character(len=BUFSIZE)  :: probeBoundsExtension, prefixFieldExtension
         character(len=BUFSIZE) :: outputPath
         probeBoundsExtension = get_probe_bounds_extension()
         prefixFieldExtension = get_prefix_extension(field, mpidir)
         outputPath = &
            trim(adjustl(outputTypeExtension))//'_'//trim(adjustl(prefixFieldExtension))//trim(adjustl(probeBoundsExtension))
         return
      end function get_output_path

      function get_probe_bounds_extension() result(ext)
         character(len=BUFSIZE) :: ext
         character(len=BUFSIZE)  ::  chari, charj, chark, chari2, charj2, chark2

         write (chari, '(i7)') iCoord
         write (charj, '(i7)') jCoord
         write (chark, '(i7)') kCoord

         write (chari2, '(i7)') i2Coord
         write (charj2, '(i7)') j2Coord
         write (chark2, '(i7)') k2Coord

#if CompileWithMPI
         if (mpidir == 3) then
            ext = trim(adjustl(chari))//'_'//trim(adjustl(charj))//'_'//trim(adjustl(chark))//'__'// &
                  trim(adjustl(chari2))//'_'//trim(adjustl(charj2))//'_'//trim(adjustl(chark2))
         elseif (mpidir == 2) then
            ext = trim(adjustl(charj))//'_'//trim(adjustl(chark))//'_'//trim(adjustl(chari))//'__'// &
                  trim(adjustl(charj2))//'_'//trim(adjustl(chark2))//'_'//trim(adjustl(chari2))
         elseif (mpidir == 1) then
            ext = trim(adjustl(chark))//'_'//trim(adjustl(chari))//'_'//trim(adjustl(charj))//'__'// &
                  trim(adjustl(chark2))//'_'//trim(adjustl(chari2))//'_'//trim(adjustl(charj2))
         else
            call stoponerror('Buggy error in mpidir. ')
         end if
#else
         ext = trim(adjustl(chari))//'_'//trim(adjustl(charj))//'_'//trim(adjustl(chark))//'__'// &
               trim(adjustl(chari2))//'_'//trim(adjustl(charj2))//'_'//trim(adjustl(chark2))
#endif

         return
      end function get_probe_bounds_extension

   end subroutine init_bulk_probe_output

   subroutine update_bulk_probe_output(this, step, field)
      type(bulk_probe_output_t), intent(out) :: this
      real(kind=RKIND_tiempo), intent(in) :: step
      type(field_data_t), intent(in) :: field

      integer(kind=SINGLE) :: i1_m, i2_m, j1_m, j2_m, k1_m, k2_m
      integer(kind=SINGLE) :: i1, i2, j1, j2, k1, k2

      real(kind=RKIND), pointer, dimension(:,:,:) :: xF, yF, zF
      real(kind=RKIND), pointer, dimension(:) :: dx, dy, dz

      i1_m = this%xCoord
      i2_m = this%x2Coord
      j1_m = this%yCoord
      j2_m = this%y2Coord
      k1_m = this%zCoord
      k2_m = this%z2Coord

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
      selectcase (field)
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
               (zF(i1_m, j1_m, KKK_m) - zF(i1_m, j2_m + 1, KKK_m))*dz(KKK_m)
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

end module mod_bulkProbe
