integer function test_init_point_probe() bind(c) result(err)
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
   logical :: ThereAreWires = .false.

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

   test_err = test_err + assert_integer_equal(outputs(1)%outputID, POINT_PROBE_ID, 'Unexpected probe id')
   test_err = test_err + assert_integer_equal(outputs(1)%pointProbe%columnas, 2, 'Unexpected number of columns')
   test_err = test_err + assert_string_equal(outputs(1)%pointProbe%path, 'entradaRoot_poinProbe_Ex_4_4_4', 'Unexpected path')
   
   deallocate (dummysgg%Observation)
   deallocate (outputs)
   err = test_err
end function test_init_point_probe

integer function test_update_point_probe() bind(c) result(err)
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
   logical :: ThereAreWires = .false.

   type(Obses_t) :: pointProbeObservable
   type(dummyFields_t), target :: dummyfields
   type(fields_reference_t) :: fields
   type(XYZlimit_t), dimension(6) :: alloc
   REAL(KIND=RKIND), DIMENSION(:,:,:), POINTER :: temp_ptr => NULL()

   real(kind=rkind) :: fieldValue
   real(kind=RKIND_tiempo) :: timestep


   integer(kind=SINGLE) :: test_err = 0

   dummysgg = create_base_sgg(dt=0.1_RKIND_tiempo, time_steps=100)
   pointProbeObservable = create_point_probe_observable()
   call add_observation_to_sgg(dummysgg, pointProbeObservable)
   dummyControl = create_control_flags(mpidir=3, nEntradaRoot='entradaRoot', wiresflavor='holland')

   call init_outputs(dummysgg, dummymedia, dummysinpml_fullsize, dummyControl, outputs, ThereAreWires)

   call create_dummy_fields(dummyfields, 1, 10, 0.01)
   
   dummyfields%Ex(4,4,4) = 5

   fields%E%x => dummyfields%Ex
   fields%E%y => dummyfields%Ey
   fields%E%z => dummyfields%Ez
   fields%E%deltax => dummyfields%dxe
   fields%E%deltaY => dummyfields%dye
   fields%E%deltaZ => dummyfields%dze
   fields%H%x => dummyfields%Hx
   fields%H%y => dummyfields%Hy
   fields%H%z => dummyfields%Hz
   fields%H%deltax => dummyfields%dxh
   fields%H%deltaY => dummyfields%dyh
   fields%H%deltaZ => dummyfields%dzh


   call update_outputs(outputs, dummyControl, 0.5_RKIND_tiempo, fields)

   test_err = test_err + assert_real_equal(outputs%pointProbe(1)%timeStep(1), 0.5_RKIND_tiempo, 0.00001_RKIND_tiempo, 'Unexpected timestep')
   test_err = test_err + assert_real_equal(outputs%pointProbe(1)%valueForTime(1), 5,  0.00001_RKIND_tiempo, 'Unexpected field')

   dummyfields%Ex(4,4,4) = -4

   call update_outputs(outputs, dummyControl, 0.8_RKIND_tiempo, fields)

   test_err = test_err + assert_real_equal(outputs%pointProbe%timeStep(2), 0.8_RKIND_tiempo, 0.00001_RKIND_tiempo, 'Unexpected timestep')
   test_err = test_err + assert_real_equal(outputs%pointProbe%valueForTime(2), -4,  0.00001_RKIND_tiempo, 'Unexpected field')


   err = test_err
end function test_update_point_probe
