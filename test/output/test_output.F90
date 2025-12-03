function test_initialize() bind(C) result(err)
   use FDETYPES
   use FDETYPES_TOOLS
   use output

   type(SGGFDTDINFO) :: dummysgg
   type(sim_control_t) :: dummyControl
   type(media_matrices_t), pointer:: dummymedia
   type(limit_t), pointer, dimension(:) :: dummysinpml_fullsize
   type(solver_output_t), dimension(:), allocatable :: outputs
   logical :: ThereAreWires = .true.

   integer(kind=SINGLE) :: test_err = 0

   !Set requested observables
   dummysgg = create_base_sgg(nummedia=5, dt=0.1_RKIND_tiempo, time_steps=100)
   dummysgg%NumberRequest = 3
   allocate (dummysgg%Observation(3))
   dummysgg%Observation(1) = define_point_observation()
   dummysgg%Observation(2) = define_wire_current_observation()
   dummysgg%Observation(3) = define_wire_charge_observation()

   !Set control flags
   dummyControl = create_control_flags(mpidir=3, nEntradaRoot='entradaRoot', wiresflavor='holland')

   call init_outputs(dummysgg, dummymedia, dummysinpml_fullsize, dummyControl, outputs, ThereAreWires)

   deallocate (dummysgg%Observation)
   deallocate (outputs)
   err = test_err
end function test_initialize

