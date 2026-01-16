module mod_testOutputUtils
   use FDETYPES
   use FDETYPES_TOOLS
   use outputTypes
   implicit none
   private

   !===========================
   !  Public interface summary
   !===========================
   public :: dummyFields_t
   public :: create_point_probe_observation
   public :: create_volumic_probe_observation
   public :: create_movie_observation
   public :: create_frequency_slice_observation
   public :: create_dummy_fields
   public :: fillGradient
   public :: setup_dummy_problem_info
   public :: clean_dummy_problem_info
   !===========================
   
   ! Storage for dummy targets
   type(media_matrices_t), target :: dummyGeometry
   type(limit_t), allocatable, target :: dummyProblemDim(:)
   type(MediaData_t), allocatable, target :: dummyMaterialList(:)
   type(bounds_t), target :: dummyBounds

   !===========================
   !  Private interface summary
   !===========================

   !===========================

   type :: dummyFields_t
      real(kind=RKIND), allocatable, dimension(:, :, :) :: Ex, Ey, Ez, Hx, Hy, Hz
      real(kind=RKIND), allocatable, dimension(:) :: dxe, dye, dze, dxh, dyh, dzh
   contains
      procedure, public :: createDummyFields => create_dummy_fields
   end type dummyFields_t
contains
   function create_point_probe_observation(x, y, z) result(obs)
      integer, intent(in) :: x, y, z
      type(Obses_t) :: obs

      type(observable_t), dimension(:), allocatable :: P
      type(observation_domain_t) :: domain

      allocate (P(1))
      P(1) = create_observable(x, y, z, x, y, z, iEx)
      call initialize_observation_time_domain(domain, 0.0_RKIND, 10.0_RKIND, 0.1_RKIND)

      call set_observation(obs, P, 'poinProbe', domain, 'DummyFileNormalize')
   end function

   function create_volumic_probe_observation(xi, yi, zi, xe, ye, ze) result(obs)
      integer, intent(in) :: xi, yi, zi, xe, ye, ze
      type(Obses_t) :: obs

      type(observable_t), dimension(:), allocatable :: P
      type(observation_domain_t) :: domain

      allocate (P(1))
      P(1) = create_observable(xi, yi, zi, xe, ye, ze, iCurX)

      call initialize_observation_time_domain(domain, 0.0_RKIND, 10.0_RKIND, 0.1_RKIND)
      call initialize_observation_frequency_domain(domain, 0.0_RKIND, 1000.0_RKIND, 50.0_RKIND)

      call set_observation(obs, P, 'volumicProbe', domain, 'DummyFileNormalize')
   end function create_volumic_probe_observation

   function create_movie_observation(xi, yi, zi, xe, ye, ze, request) result(observation)
      integer, intent(in) :: xi, yi, zi, xe, ye, ze, request
      type(Obses_t) :: observation

      type(observable_t), dimension(:), allocatable :: P
      type(observation_domain_t) :: domain

      allocate (P(1))
      P(1) = create_observable(xi, yi, zi, xe, ye, ze, request)
      call initialize_observation_time_domain(domain, 0.0_RKIND, 10.0_RKIND, 0.1_RKIND)

      call set_observation(observation, P, 'movieProbe', domain, 'DummyFileNormalize')
   end function create_movie_observation

   function create_frequency_slice_observation(xi, yi, zi, xe, ye, ze, request) result(observation)
      integer, intent(in) :: xi, yi, zi, xe, ye, ze, request
      type(Obses_t) :: observation

      type(observable_t), dimension(:), allocatable :: P
      type(observation_domain_t) :: domain

      allocate (P(1))
      P(1) = create_observable(xi, yi, zi, xe, ye, ze, request)
      call initialize_observation_frequency_domain(domain, 0.0_RKIND, 100.0_RKIND, 20.0_RKIND)

      call set_observation(observation, P, 'frequency_sliceProbe', domain, 'DummyFileNormalize')
   end function create_frequency_slice_observation

   subroutine create_dummy_fields(this, lower, upper, delta)
      class(dummyFields_t), intent(inout) :: this
      integer, intent(in) :: lower, upper
      real(kind=rkind), intent(in) :: delta
      allocate ( &
         this%Ex(lower:upper, lower:upper, lower:upper), &
         this%Ey(lower:upper, lower:upper, lower:upper), &
         this%Ez(lower:upper, lower:upper, lower:upper), &
         this%Hx(lower:upper, lower:upper, lower:upper), &
         this%Hy(lower:upper, lower:upper, lower:upper), &
         this%Hz(lower:upper, lower:upper, lower:upper) &
         )

      this%Ex = 0.0_RKIND
      this%Ey = 0.0_RKIND
      this%Ez = 0.0_RKIND
      this%Hx = 0.0_RKIND
      this%Hy = 0.0_RKIND
      this%Hz = 0.0_RKIND

      allocate ( &
         this%dxh(lower:upper), &
         this%dyh(lower:upper), &
         this%dzh(lower:upper), &
         this%dxe(lower:upper), &
         this%dye(lower:upper), &
         this%dze(lower:upper) &
         )
      this%dxh = delta
      this%dyh = delta
      this%dzh = delta
      this%dxe = delta
      this%dye = delta
      this%dze = delta
   end subroutine create_dummy_fields

   subroutine fillGradient(dummyFields, direction, minVal, maxVal)
      !--------------------------------------------
      ! Fills dummyFields%Hx, Hy, Hz with a linear gradient
      ! along the specified direction (1=x, 2=y, 3=z)
      !--------------------------------------------
      implicit none
      type(dummyFields_t), intent(inout) :: dummyFields
      integer, intent(in) :: direction       ! 1=x, 2=y, 3=z
      real(RKIND), intent(in) :: minVal, maxVal

      integer :: i, j, k
      integer :: nx, ny, nz
      real(RKIND) :: factor

      ! Get array sizes
      nx = size(dummyFields%Hx, 1)
      ny = size(dummyFields%Hx, 2)
      nz = size(dummyFields%Hx, 3)

      select case (direction)
      case (1)  ! x-direction
         do i = 1, nx
            factor = real(i - 1, RKIND)/real(nx - 1, RKIND)
            dummyFields%Hx(i, :, :) = minVal + factor*(maxVal - minVal)
            dummyFields%Hy(i, :, :) = minVal + factor*(maxVal - minVal)
            dummyFields%Hz(i, :, :) = minVal + factor*(maxVal - minVal)
         end do
      case (2)  ! y-direction
         do j = 1, ny
            factor = real(j - 1, RKIND)/real(ny - 1, RKIND)
            dummyFields%Hx(:, j, :) = minVal + factor*(maxVal - minVal)
            dummyFields%Hy(:, j, :) = minVal + factor*(maxVal - minVal)
            dummyFields%Hz(:, j, :) = minVal + factor*(maxVal - minVal)
         end do
      case (3)  ! z-direction
         do k = 1, nz
            factor = real(k - 1, RKIND)/real(nz - 1, RKIND)
            dummyFields%Hx(:, :, k) = minVal + factor*(maxVal - minVal)
            dummyFields%Hy(:, :, k) = minVal + factor*(maxVal - minVal)
            dummyFields%Hz(:, :, k) = minVal + factor*(maxVal - minVal)
         end do
      case default
         print *, "Error: direction must be 1, 2, or 3."
      end select

   end subroutine fillGradient

   !--------------------------------------------------------------------------------
   ! Setup/Teardown
   !--------------------------------------------------------------------------------
   subroutine setup_dummy_problem_info(problemInfo)
      type(problem_info_t), intent(out) :: problemInfo

      integer :: i
      
      ! Create a 11x11x11 grid (0..10)
      if (allocated(dummyProblemDim)) deallocate(dummyProblemDim)
      allocate(dummyProblemDim(6))
      do i = 1,6 
         dummyProblemDim(i) = create_limit_t(0, 10, 0, 10, 0, 10, 1, 1, 1)
      end do
      problemInfo%problemDimension => dummyProblemDim

      call create_geometry_media(dummyGeometry, 0, 10, 0, 10, 0, 10)
      problemInfo%geometryToMaterialData => dummyGeometry

      call init_simulation_material_list(dummyMaterialList)
      problemInfo%materialList => dummyMaterialList
      
      problemInfo%simulationBounds => dummyBounds

   end subroutine setup_dummy_problem_info

   subroutine clean_dummy_problem_info(problemInfo)
      type(problem_info_t), intent(inout) :: problemInfo
      
      if (allocated(dummyProblemDim)) deallocate(dummyProblemDim)
      if (allocated(dummyMaterialList)) deallocate(dummyMaterialList)
      
      nullify(problemInfo%problemDimension)
      nullify(problemInfo%geometryToMaterialData)
      nullify(problemInfo%materialList)
      nullify(problemInfo%simulationBounds)
   end subroutine clean_dummy_problem_info

end module mod_testOutputUtils
