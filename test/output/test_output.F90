integer function test_init_point_probe() bind(c) result(err)
   use FDETYPES
   use FDETYPES_TOOLS
   use output
   use mod_testOutputUtils

   type(SGGFDTDINFO) :: dummysgg
   type(sim_control_t) :: dummyControl
   type(media_matrices_t), pointer:: dummymedia => NULL()
   type(limit_t), pointer, dimension(:) :: dummysinpml_fullsize => NULL()
   type(solver_output_t), dimension(:), allocatable :: outputs
   type(MediaData_t) :: defaultMaterial, pecMaterial
   logical :: ThereAreWires = .false.

   type(Obses_t) :: pointProbeObservable

   integer(kind=SINGLE) :: test_err = 0

   !Cleanup
   if (allocated(outputs)) deallocate (outputs)

   !Set requested observables
   dummysgg = create_base_sgg(dt=0.1_RKIND_tiempo, time_steps=100)

   pointProbeObservable = create_point_probe_observable()
   call add_observation_to_sgg(dummysgg, pointProbeObservable)

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
   type(media_matrices_t), pointer:: dummymedia => NULL()
   type(limit_t), pointer, dimension(:) :: dummysinpml_fullsize => NULL()
   type(solver_output_t), dimension(:), allocatable :: outputs
   logical :: ThereAreWires = .false.

   type(Obses_t) :: pointProbeObservable
   type(dummyFields_t), target :: dummyfields
   type(fields_reference_t) :: fields

   integer(kind=SINGLE) :: test_err = 0

   dummysgg = create_base_sgg(dt=0.1_RKIND_tiempo, time_steps=100)
   pointProbeObservable = create_point_probe_observable()
   call add_observation_to_sgg(dummysgg, pointProbeObservable)
   dummyControl = create_control_flags(mpidir=3, nEntradaRoot='entradaRoot', wiresflavor='holland')

   call init_outputs(dummysgg, dummymedia, dummysinpml_fullsize, dummyControl, outputs, ThereAreWires)

   call create_dummy_fields(dummyfields, 1, 10, 0.01)

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

   dummyfields%Ex(4, 4, 4) = 5.0_RKIND
   call update_outputs(outputs, dummyControl, 0.5_RKIND_tiempo, fields)

   test_err = test_err + assert_real_time_equal(outputs(1)%pointProbe%timeStep(1), 0.5_RKIND_tiempo, 0.00001_RKIND_tiempo, 'Unexpected timestep 1')
   test_err = test_err + assert_real_equal(outputs(1)%pointProbe%valueForTime(1), 5.0_RKIND, 0.00001_RKIND, 'Unexpected field 1')

   dummyfields%Ex(4, 4, 4) = -4.0_RKIND
   call update_outputs(outputs, dummyControl, 0.8_RKIND_tiempo, fields)

   test_err = test_err + assert_real_time_equal(outputs(1)%pointProbe%timeStep(2), 0.8_RKIND_tiempo, 0.00001_RKIND_tiempo, 'Unexpected timestep 2')
   test_err = test_err + assert_real_equal(outputs(1)%pointProbe%valueForTime(2), -4.0_RKIND, 0.00001_RKIND, 'Unexpected field 2')

   if (associated(dummymedia)) deallocate (dummymedia)
   if (associated(dummysinpml_fullsize)) deallocate (dummysinpml_fullsize)

   err = test_err
end function test_update_point_probe

integer function test_flush_point_probe() bind(c) result(err)
   use output
   use mod_domain
   use mod_testOutputUtils
   type(point_probe_output_t) :: probe
   type(domain_t):: domain
   character(len=BUFSIZE) :: file_time, file_freq
   character(len=27) :: test_extension
   integer :: n, i
   err = 1 !If test_err is not updated at the end it will be shown
   test_err = 0
   test_extension = 'tmp_cases/flush_point_probe'
   domain = create_domain(0.0_RKIND_tiempo, 10.0_RKIND_tiempo, 0.1_RKIND_tiempo, 10.0_RKIND, 100.0_RKIND, 10, .false.)
   call init_point_probe_output(probe, 2, 2, 2, iEx, domain, test_extension, 3)
   call create_point_probe_output_files(probe)

   n = 10
   do i = 1, n
      probe%timeStep(i) = real(i)
      probe%valueForTime(i) = 10.0*i
      probe%frequencySlice(i) = 0.1*i
      probe%valueForFreq(i) = 0.2*i
   end do
   probe%serializedTimeSize = n
   probe%nFreq = n

   file_time = trim(adjustl(probe%path))//'_'// &
               trim(adjustl(timeExtension))//'_'// &
               trim(adjustl(datFileExtension))

   file_freq = trim(adjustl(probe%path))//'_'// &
               trim(adjustl(timeExtension))//'_'// &   ! <-- SAME naming in your code
               trim(adjustl(datFileExtension))

   call flush_point_probe_output(probe)

   test_err = test_err + assert_written_output_file(file_time)
   test_err = test_err + assert_written_output_file(file_freq)

 test_err = test_err + assert_integer_equal(probe%serializedTimeSize, 0, "ERROR: clear_time_data did not reset serializedTimeSize!")
 test_err = test_err + assert_integer_equal(probe%serializedTimeSize, 0, "ERROR: clear_time_data did not reset serializedTimeSize!")

   if (all(probe%timeStep == 0.0) .and. all(probe%valueForTime == 0.0)) then
      print *, "Time arrays cleared correctly."
   else
      print *, "ERROR: time arrays not cleared!"
      test_err = test_err + 1
   end if

   if (probe%nFreq == 0) then
      print *, "ERROR: Destroyed frequency reference!"
      test_err = test_err + 1
   end if

   err = test_err
end function test_flush_point_probe

integer function test_multiple_flush_point_probe() bind(c) result(err)
   use output
   use mod_domain
   use mod_testOutputUtils
   type(point_probe_output_t) :: probe
   type(domain_t):: domain
   character(len=BUFSIZE) :: file_time, file_freq
   real(kind=RKIND), allocatable :: expectedTime(:,:), expectedFreq(:,:)
   character(len=36) :: test_extension
   integer :: n, i
   err = 1 !If test_err is not updated at the end it will be shown
   test_err = 0
   test_extension = 'tmp_cases/multiple_flush_point_probe'

   domain = create_domain(0.0_RKIND_tiempo, 10.0_RKIND_tiempo, 0.1_RKIND_tiempo, 10.0_RKIND, 100.0_RKIND, 10, .false.)
   call init_point_probe_output(probe, 2, 2, 2, iEx, domain, test_extension, 3)
   call create_point_probe_output_files(probe)

   file_time = trim(adjustl(probe%path))//'_'// &
               trim(adjustl(timeExtension))//'_'// &
               trim(adjustl(datFileExtension))

   file_freq = trim(adjustl(probe%path))//'_'// &
               trim(adjustl(frequencyExtension))//'_'// &
               trim(adjustl(datFileExtension))

   n = 10

   allocate(expectedTime(2*n,2))
   allocate(expectedFreq(n,2))
   !Simulate updates in probe
   do i = 1, n
      probe%timeStep(i) = real(i)
      probe%valueForTime(i) = 10.0*i
      probe%frequencySlice(i) = 0.1*i
      probe%valueForFreq(i) = 0.2*i

      expectedTime(i,1) = real(i)
      expectedTime(i,2) = 10.0*i

      expectedFreq(i,1) = 0.1*i
      expectedFreq(i,2) = 0.2*i
   end do
   probe%serializedTimeSize = n
   probe%nFreq = n
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   call flush_point_probe_output(probe)

   !Simulate new updates in probe
   do i = 1, n
      probe%timeStep(i) = real(i+10)
      probe%valueForTime(i) = 10.0*(i+10)
      probe%valueForFreq(i) = -0.5*i

      expectedTime(i+10,1) = real(i+10)
      expectedTime(i+10,2) = 10.0*(i+10)

      expectedFreq(i,1) = 0.1*i  ! frequency file overwrites, so expectedFreq(i,1) remains 0.1*i ?
      expectedFreq(i,2) = -0.5*i
   end do
   probe%serializedTimeSize = n
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   call flush_point_probe_output(probe)

   open(unit=probe%fileUnitTime, file=file_time, status="old", action="read")
   test_err = test_err + assert_file_content(probe%fileUnitTime, expectedTime, 2*n, 2)
   close(probe%fileUnitTime)

   open(unit=probe%fileUnitFreq, file=file_freq, status="old", action="read")
   test_err = test_err + assert_file_content(probe%fileUnitFreq, expectedFreq, n, 2)
   close(probe%fileUnitFreq)

   err = test_err


end function test_multiple_flush_point_probe
