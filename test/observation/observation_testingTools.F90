module observation_testingTools
   use FDETYPES
   type :: dummyFields_t
      real(kind=RKIND),allocatable, dimension(:,:,:) :: Ex, Ey, Ez, Hx, Hy, Hz
      real(kind=RKIND),allocatable, dimension(:) :: dxe, dye, dze, dxh, dyh, dzh
      contains
      procedure, public :: createDummyFields => create_dummy_fields
   end type dummyFields_t

   public dummyFields
   public check_shape_real, check_shape_complex, check_size
   public approx_equal
contains
   subroutine create_dummy_fields(this, lower, upper, delta)
      class(dummyFields_t), intent(inout) :: this
      integer, intent(in) :: lower, upper
      real(kind=rkind), intent(in) :: delta
      allocate(&
         this%Ex(lower:upper, lower:upper, lower:upper),&
         this%Ey(lower:upper, lower:upper, lower:upper),&
         this%Ez(lower:upper, lower:upper, lower:upper),&
         this%Hx(lower:upper, lower:upper, lower:upper),&
         this%Hy(lower:upper, lower:upper, lower:upper),&
         this%Hz(lower:upper, lower:upper, lower:upper)&
      )

      this%Ex = 0.0_RKIND
      this%Ey = 0.0_RKIND
      this%Ez = 0.0_RKIND
      this%Hx = 0.0_RKIND
      this%Hy = 0.0_RKIND
      this%Hz = 0.0_RKIND

      allocate(&
         this%dxh(lower:upper), &
         this%dyh(lower:upper), &
         this%dzh(lower:upper), &
         this%dxe(lower:upper), &
         this%dye(lower:upper), &
         this%dze(lower:upper)&
      )
      this%dxh = delta
      this%dyh = delta
      this%dzh = delta
      this%dxe = delta
      this%dye = delta
      this%dze = delta
   end subroutine create_dummy_fields

   subroutine check_shape_real(arr, n_expected, test_err, name)
      use Observa
      use FDETYPES
      real(kind=RKIND), intent(in), dimension(:, :) :: arr
      integer, intent(in) :: n_expected
      integer, intent(inout) :: test_err
      character(len=*), intent(in), optional :: name

      integer :: rank_arr
      integer, allocatable :: shp(:)
      character(len=:), allocatable :: nm

      if (present(name)) then
          nm = trim(adjustl(name))
      else
          nm = "array"
      end if

      rank_arr = rank(arr)
      shp = shape(arr)

      ! Expect exactly two dimensions [1, n_expected]
      if (rank_arr /= 2 .or. any(shp /= [1, n_expected])) then
         test_err = test_err + 1
      end if
   end subroutine check_shape_real

   subroutine check_shape_complex(arr, n_expected, test_err, name)
      use Observa
      use FDETYPES
      complex(kind=CKIND), intent(in), dimension(:, :) :: arr
      integer, intent(in) :: n_expected
      integer, intent(inout) :: test_err
      character(len=*), intent(in), optional :: name

      integer :: rank_arr
      integer, allocatable :: shp(:)
      character(len=:), allocatable :: nm

      if (present(name)) then
          nm = trim(adjustl(name))
      else
          nm = "array"
      end if

      rank_arr = rank(arr)
      shp = shape(arr)

      ! Expect exactly two dimensions [1, n_expected]
      if (rank_arr /= 2 .or. any(shp /= [1, n_expected])) then
         test_err = test_err + 1
      end if
   end subroutine check_shape_complex

   subroutine check_size(arr, n_expected, test_err, name)
      use Observa
      use FDETYPES
      integer, intent(in), dimension(:) :: arr
      integer, intent(in) :: n_expected
      integer, intent(inout) :: test_err
      character(len=*), intent(in), optional :: name

      integer :: rank_arr
      integer :: siz
      character(len=:), allocatable :: nm

      if (present(name)) then
          nm = trim(adjustl(name))
      else
          nm = "array"
      end if

      rank_arr = rank(arr)
      siz = size(arr)

      if (rank_arr /= 1 .or. siz /= n_expected) then
         test_err = test_err + 1
      end if
   end subroutine check_size

   logical function approx_equal(a, b, tol) result(equal)
      use FDETYPES
      real(kind=RKIND), intent(in) :: a, b, tol
      equal = abs(a - b) <= tol
   end function approx_equal

   function create_time_array(array_size, interval) result(arr)
      use FDETYPES
      integer, intent(in) :: array_size
      integer(kind=4) :: i
      real(kind=RKIND_tiempo) :: interval

      real(kind=RKIND_tiempo), pointer, dimension(:) :: arr
      allocate (arr(array_size))

      DO i = 1, array_size
         arr(i) = (i - 1)*interval
      END DO
   end function create_time_array

   function create_limit_type() result(r)
      use FDETYPES
      type(limit_t) :: r
   end function

   function create_xyz_limit_array(XI,YI,ZI,XE,YE,ZE) result(arr)
      use FDETYPES
      type(XYZlimit_t), dimension(1:6) :: arr
      integer (kind=4), intent(in) :: XI,YI,ZI,XE,YE,ZE
      integer :: i
      do i = 1, 6
         arr(i)%XI = XI
         arr(i)%XE = XE
         arr(i)%YI = YI
         arr(i)%YE = YE
         arr(i)%ZI = ZI
         arr(i)%ZE = ZE
      end do
   end function create_xyz_limit_array

   
   function create_facesNF2FF(tr, fr, iz, de, ab, ar) result(faces)
      use FDETYPES
      type(nf2ff_t) :: faces
      logical :: tr, fr, iz, de, ab, ar

      faces%tr = tr
      faces%fr = fr
      faces%iz = iz
      faces%de = de
      faces%ab = ab
      faces%ar = ar
   end function create_facesNF2FF

   function create_control_flags(layoutnumber, size, mpidir, finaltimestep, &
                                       nEntradaRoot, wiresflavor, &
                                       resume, saveall, NF2FFDecim, simu_devia, singlefilewrite, &
                                       facesNF2FF) result(control)
      use FDETYPES
      type(sim_control_t) :: control
      integer(kind=4), intent(in) :: layoutnumber, size, mpidir, finaltimestep
      character(len=*), intent(in) :: nEntradaRoot, wiresflavor
      logical, intent(in) :: resume, saveall, NF2FFDecim, simu_devia, singlefilewrite
      type(nf2ff_t), intent(in) :: facesNF2FF

      control%layoutnumber  = layoutnumber
      control%size = size
      control%mpidir  = mpidir
      control%finaltimestep = finaltimestep
      control%nEntradaRoot  = nEntradaRoot
      control%wiresflavor = wiresflavor
      control%resume  = resume
      control%saveall = saveall
      control%NF2FFDecim = NF2FFDecim
      control%simu_devia = simu_devia
      control%singlefilewrite  = singlefilewrite
      control%facesNF2FF = facesNF2FF

   end function create_control_flags

   function create_base_sgg() result(sgg)
      use FDETYPES
      type(SGGFDTDINFO) :: sgg
      
      sgg%NumMedia = 3
      allocate(sgg%Med(0:sgg%NumMedia))
      sgg%Med = create_basic_media()
      sgg%NumberRequest = 1
      sgg%dt = 0.1_RKIND_tiempo
      sgg%tiempo => create_time_array(100, sgg%dt)
      sgg%Sweep = create_xyz_limit_array(0,0,0,6,6,6)
      sgg%SINPMLSweep = create_xyz_limit_array(1,1,1,5,5,5)
      sgg%NumPlaneWaves = 1
      sgg%alloc = create_xyz_limit_array(0,0,0,6,6,6)
 
   end function create_base_sgg

   function create_basic_media () result(media)
      use FDETYPES
      type(MediaData_t) :: media
   end function create_basic_media
end module
