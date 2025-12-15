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
   public :: create_dummy_fields
   !===========================

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

   function create_movie_observation(xi, yi, zi, xe, ye, ze) result(observation)
      integer, intent(in) :: xi, yi, zi, xe, ye, ze
      type(Obses_t) :: observation

      type(observable_t), dimension(:), allocatable :: P
      type(observation_domain_t) :: domain

      allocate (P(1))
      P(1) = create_observable(xi, yi, zi, xe, ye, ze, iCur)
      call initialize_observation_time_domain(domain, 0.0_RKIND, 10.0_RKIND, 0.1_RKIND)

      call set_observation(observation, P, 'movieProbe', domain, 'DummyFileNormalize')
   end function create_movie_observation

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

end module mod_testOutputUtils
