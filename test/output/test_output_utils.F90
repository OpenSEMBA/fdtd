module mod_testOutputUtils
   use FDETYPES
   use FDETYPES_TOOLS

   implicit none
contains
   function create_point_probe_observable() result(obs)
      type(Obses_t) :: obs

      type(observable_t), dimension(:), allocatable :: P
      type(observation_domain_t) :: domain

      call initialize_time_domain(domain, 0.0_RKIND, 10.0_RKIND, 0.1_RKIND)
      allocate(P(1))
      P(1) = create_observable(0, 0, 0, 6, 6, 6, iEx)
      call set_observable(obs, P, 'poinProbe', domain, 'DummyFileNormalize')

   end function
end module mod_testOutputUtils
