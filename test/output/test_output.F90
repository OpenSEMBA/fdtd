function test_initialize() bind(C) result(err)
   use FDETYPES
   use FDETYPES_TOOLS
   use output
   use mod_testOutputUtils

   type(SGGFDTDINFO) :: dummysgg
   type(sim_control_t) :: dummyControl
   type(media_matrices_t), pointer:: dummymedia
   type(limit_t), pointer, dimension(:) :: dummysinpml_fullsize
   type(solver_output_t), dimension(:), allocatable :: outputs
   type(MediaData_t) :: defaultMaterial, pecMaterial
   logical :: ThereAreWires = .true.

   type(Obses_t) :: pointProbeObservable

   integer(kind=SINGLE) :: test_err = 0
   if (allocated(outputs)) deallocate(outputs)

   !Set requested observables
   dummysgg = create_base_sgg(dt=0.1_RKIND_tiempo, time_steps=100)
   
   pointProbeObservable = create_point_probe_observable()
   call add_observation_to_sgg(dummysgg, pointProbeObservable)

   !Set dummymedia

   !set dummysinpml_fullsize

   !Set control flags
   dummyControl = create_control_flags(mpidir=3, nEntradaRoot='entradaRoot', wiresflavor='holland')

   call init_outputs(dummysgg, dummymedia, dummysinpml_fullsize, dummyControl, outputs, ThereAreWires)

   call clean_solver_output_array(outputs)
   err = test_err
end function test_initialize


function test_init_point_probe() bind(c) result(err)
   use FDETYPES
   use FDETYPES_TOOLS
   use output
   use mod_testOutputUtils

   type(SGGFDTDINFO) :: dummysgg
   type(sim_control_t) :: dummyControl
   type(media_matrices_t), pointer:: dummymedia
   type(limit_t), pointer, dimension(:) :: dummysinpml_fullsize
   type(solver_output_t), dimension(:), allocatable :: outputs
   type(MediaData_t) :: defaultMaterial, pecMaterial
   logical :: ThereAreWires = .true.

   type(Obses_t) :: pointProbeObservable

   integer(kind=SINGLE) :: test_err = 0

   !Cleanup 
   if (allocated(outputs)) deallocate(outputs)

   !Set requested observables
   dummysgg = create_base_sgg(dt=0.1_RKIND_tiempo, time_steps=100)
   
   pointProbeObservable = create_point_probe_observable()
   call add_observation_to_sgg(dummysgg, pointProbeObservable)

   !Set dummymedia

   !set dummysinpml_fullsize

   !Set control flags
   dummyControl = create_control_flags(mpidir=3, nEntradaRoot='entradaRoot', wiresflavor='holland')

   call init_outputs(dummysgg, dummymedia, dummysinpml_fullsize, dummyControl, outputs, ThereAreWires)

   test_err = assert_integer_equal(outputs(1)%outputID, POINT_PROBE_ID, 'Unexpected probe id')
   test_err = assert_integer_equal(outputs(1)%pointProbe%columnas, 2, 'Unexpected number of columns')
   test_err = assert_string_equal(outputs(1)%pointProbe%path, 'test', 'Unexpected path')
   
   deallocate (dummysgg%Observation)
   deallocate (outputs)
   err = test_err
end function test_init_point_probe
