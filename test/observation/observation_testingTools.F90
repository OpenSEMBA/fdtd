module observation_testingTools
   implicit none
   public check_shape_real, check_shape_complex, check_size
   public approx_equal
contains
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

      nm = merge(name, "array", present(name))

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

      nm = merge(name, "array", present(name))

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

      nm = merge(name, "array", present(name))

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

   function create_sgg_sweep() result(arr)
      use FDETYPES
      type(XYZlimit_t), dimension(1:6) :: arr
      integer :: i
      do i = 1, 6
         arr(i)%XI = 0
         arr(i)%XE = 0
         arr(i)%YI = 0
         arr(i)%YE = 0
         arr(i)%ZI = 0
         arr(i)%ZE = 0
      end do
   end function create_sgg_sweep

   function create_sgg_alloc() result(arr)
      use FDETYPES
      type(XYZlimit_t), dimension(1:6) :: arr
      integer :: i
      do i = 1, 6
         arr(i)%XI = 0
         arr(i)%XE = 0
         arr(i)%YI = 0
         arr(i)%YE = 5
         arr(i)%ZI = 5
         arr(i)%ZE = 5
      end do
   end function create_sgg_alloc
   
   function create_sgg_SINPMLSweep() result(arr)
      use FDETYPES
      type(XYZlimit_t), dimension(1:6) :: arr
      integer :: i
      do i = 1, 6
         arr(i)%XI = 0
         arr(i)%XE = 0
         arr(i)%YI = 0
         arr(i)%YE = 5
         arr(i)%ZI = 5
         arr(i)%ZE = 5
      end do
   end function create_sgg_SINPMLSweep

   function create_facesNF2FF(tr, fr, iz, de, ab, ar) result(faces)
      use FDETYPES
      type(nf2ff_t) :: faces
      logical :: tr, fr, iz, de, ab, ar

      faces%tr = .false.
      faces%fr = .false.
      faces%iz = .false.
      faces%de = .false.
      faces%ab = .false.
      faces%ar = .false.
   end function create_facesNF2FF

   subroutine initialize_control_flags(control, &
                                       layoutnumber, size, mpidir, finaltimestep, &
                                       nEntradaRoot, wiresflavor, &
                                       resume, saveall, NF2FFDecim, simu_devia, singlefilewrite, &
                                       facesNF2FF)
      use FDETYPES
      type(sim_control_t), intent(out) :: control
      integer(kind=4), intent(in) :: layoutnumber, size, mpidir, finaltimestep
      character(len=bufsize), intent(in) :: nEntradaRoot, wiresflavor
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

   end subroutine initialize_control_flags

   function create_base_sgg() result(sgg)
      use FDETYPES
      type(SGGFDTDINFO) :: sgg
      
      sgg%NumMedia = 3
      allocate(sgg%Med(1:sgg%NumMedia))
      sgg%Med = create_basic_media()
      sgg%NumberRequest = 1
      sgg%dt = 0.1_RKIND_tiempo
      sgg%tiempo => create_time_array(100, sgg%dt)
      sgg%Sweep = create_sgg_sweep()
      sgg%SINPMLSweep = create_sgg_SINPMLSweep()
      sgg%NumPlaneWaves = 1
      sgg%alloc = create_sgg_alloc()
 
   end function create_base_sgg

   function create_basic_media () result(media)
      use FDETYPES
      type(MediaData_t) :: media
   end function create_basic_media
end module
