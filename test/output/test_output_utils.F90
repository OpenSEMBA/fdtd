module mod_testOutputUtils
   use FDETYPES
   use FDETYPES_TOOLS

   implicit none
   type :: dummyFields_t
      real(kind=RKIND),allocatable, dimension(:,:,:) :: Ex, Ey, Ez, Hx, Hy, Hz
      real(kind=RKIND),allocatable, dimension(:) :: dxe, dye, dze, dxh, dyh, dzh
      contains
      procedure, public :: createDummyFields => create_dummy_fields
   end type dummyFields_t
contains
   function create_point_probe_observable() result(obs)
      type(Obses_t) :: obs

      type(observable_t), dimension(:), allocatable :: P
      type(observation_domain_t) :: domain

      call initialize_time_domain(domain, 0.0_RKIND, 10.0_RKIND, 0.1_RKIND)
      allocate(P(1))
      P(1) = create_observable(4, 4, 4, 6, 6, 6, iEx)
      call set_observable(obs, P, 'poinProbe', domain, 'DummyFileNormalize')

   end function

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

end module mod_testOutputUtils
