integer function test_init_point_probe() bind(c) result(err)
   use FDETYPES
   use FDETYPES_TOOLS
   use output
   use mod_testOutputUtils
   use mod_sggMethods
   use mod_assertionTools

   type(SGGFDTDINFO) :: dummysgg
   type(sim_control_t) :: dummyControl
   type(media_matrices_t), pointer:: dummymedia => NULL()
   type(MediaData_t), dimension(:), allocatable, target :: simulationMedia
   type(MediaData_t), dimension(:), pointer :: simulationMediaPtr
   type(limit_t), pointer, dimension(:) :: dummysinpml_fullsize => NULL()
   type(solver_output_t), dimension(:), allocatable :: outputs
   type(MediaData_t) :: defaultMaterial, pecMaterial
   logical :: ThereAreWires = .false.

   type(Obses_t) :: pointProbeObservable

   real(kind=RKIND_tiempo) :: dt = 0.1_RKIND_tiempo
   integer(kind=SINGLE) :: nTimeSteps = 100_SINGLE
   real(kind=RKIND_tiempo), pointer, dimension(:) :: timeArray

   integer(kind=SINGLE) :: test_err = 0

   !Cleanup
   if (allocated(outputs)) deallocate (outputs)

   !Set requested observables
   call sgg_init(dummysgg)
   call init_time_array(timeArray, nTimeSteps, dt)
   call sgg_set_tiempo(dummysgg, timeArray)
   call sgg_set_dt(dummysgg, dt)

   call init_simulation_material_list(simulationMedia)
   simulationMediaPtr => simulationMedia
   call sgg_set_Med(dummysgg, simulationMediaPtr)

   pointProbeObservable = create_point_probe_observation(4, 4, 4)
   call sgg_add_observation(dummysgg, pointProbeObservable)

   !Set control flags
   dummyControl = create_control_flags(mpidir=3, nEntradaRoot='entradaRoot', wiresflavor='holland')

   call init_outputs(dummysgg, dummymedia, dummysinpml_fullsize, dummyControl, outputs, ThereAreWires)

   test_err = test_err + assert_integer_equal(outputs(1)%outputID, POINT_PROBE_ID, 'Unexpected probe id')
   test_err = test_err + assert_integer_equal(outputs(1)%pointProbe%columnas, 2, 'Unexpected number of columns')
   test_err = test_err + assert_string_equal(outputs(1)%pointProbe%path, 'entradaRoot_poinProbe_Ex_4_4_4', 'Unexpected path')

   call close_outputs(outputs)

   deallocate (dummysgg%Observation)
   deallocate (outputs)
   err = test_err
end function test_init_point_probe

integer function test_update_point_probe() bind(c) result(err)
   use FDETYPES
   use FDETYPES_TOOLS
   use output
   use mod_testOutputUtils
   use mod_sggMethods
   use mod_assertionTools

   type(SGGFDTDINFO) :: dummysgg
   type(sim_control_t) :: dummyControl
   type(media_matrices_t), pointer:: dummymedia => NULL()
   type(MediaData_t), dimension(:), allocatable, target :: simulationMedia
   type(MediaData_t), dimension(:), pointer :: simulationMediaPtr
   type(limit_t), pointer, dimension(:) :: dummysinpml_fullsize => NULL()
   type(solver_output_t), dimension(:), allocatable :: outputs
   logical :: ThereAreWires = .false.

   type(Obses_t) :: pointProbeObservable
   type(dummyFields_t), target :: dummyfields
   type(fields_reference_t) :: fields

   real(kind=RKIND_tiempo) :: dt = 0.1_RKIND_tiempo
   integer(kind=SINGLE) :: nTimeSteps = 100_SINGLE
   real(kind=RKIND_tiempo), pointer, dimension(:) :: timeArray

   integer(kind=SINGLE) :: test_err = 0

   call sgg_init(dummysgg)
   call init_time_array(timeArray, nTimeSteps, dt)
   call sgg_set_tiempo(dummysgg, timeArray)
   call sgg_set_dt(dummysgg, dt)
   pointProbeObservable = create_point_probe_observation(4, 4, 4)
   call sgg_add_observation(dummysgg, pointProbeObservable)
   dummyControl = create_control_flags(mpidir=3, nEntradaRoot='entradaRoot', wiresflavor='holland')

   call init_simulation_material_list(simulationMedia)
   simulationMediaPtr => simulationMedia
   call sgg_set_Med(dummysgg, simulationMediaPtr)

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
   call update_outputs(outputs, dummyMedia, simulationMediaPtr, dummysinpml_fullsize, dummyControl, 0.5_RKIND_tiempo, fields)

   test_err = test_err + assert_real_time_equal(outputs(1)%pointProbe%timeStep(1), 0.5_RKIND_tiempo, 0.00001_RKIND_tiempo, 'Unexpected timestep 1')
   test_err = test_err + assert_real_equal(outputs(1)%pointProbe%valueForTime(1), 5.0_RKIND, 0.00001_RKIND, 'Unexpected field 1')

   dummyfields%Ex(4, 4, 4) = -4.0_RKIND
   call update_outputs(outputs, dummyMedia, simulationMediaPtr, dummysinpml_fullsize, dummyControl, 0.8_RKIND_tiempo, fields)

   test_err = test_err + assert_real_time_equal(outputs(1)%pointProbe%timeStep(2), 0.8_RKIND_tiempo, 0.00001_RKIND_tiempo, 'Unexpected timestep 2')
   test_err = test_err + assert_real_equal(outputs(1)%pointProbe%valueForTime(2), -4.0_RKIND, 0.00001_RKIND, 'Unexpected field 2')

   call close_outputs(outputs)

   if (associated(dummymedia)) deallocate (dummymedia)
   if (associated(dummysinpml_fullsize)) deallocate (dummysinpml_fullsize)

   err = test_err
end function test_update_point_probe

integer function test_flush_point_probe() bind(c) result(err)
   use output
   use mod_domain
   use mod_testOutputUtils
   use mod_assertionTools
   type(point_probe_output_t) :: probe
   type(domain_t):: domain
   type(cell_coordinate_t) :: coordinates
   character(len=BUFSIZE) :: file_time, file_freq
   character(len=27) :: test_extension
   integer :: n, i
   err = 1 !If test_err is not updated at the end it will be shown
   test_err = 0
   test_extension = 'tmp_cases/flush_point_probe'
   domain = domain_t(0.0_RKIND_tiempo, 10.0_RKIND_tiempo, 0.1_RKIND_tiempo, 10.0_RKIND, 100.0_RKIND, 10, .false.)

   coordinates%x = 2
   coordinates%y = 2
   coordinates%z = 2
   call init_point_probe_output(probe, coordinates, iEx, domain, test_extension, 3, 0.1_RKIND_tiempo)
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
   use mod_assertionTools
   type(point_probe_output_t) :: probe
   type(domain_t):: domain
   type(cell_coordinate_t) :: coordinates
   character(len=BUFSIZE) :: file_time, file_freq
   real(kind=RKIND), allocatable :: expectedTime(:, :), expectedFreq(:, :)
   character(len=36) :: test_extension
   integer :: n, i
   err = 1 !If test_err is not updated at the end it will be shown
   test_err = 0
   test_extension = 'tmp_cases/multiple_flush_point_probe'

   domain = domain_t(0.0_RKIND_tiempo, 10.0_RKIND_tiempo, 0.1_RKIND_tiempo, 10.0_RKIND, 100.0_RKIND, 10, .false.)
   coordinates%x = 2
   coordinates%y = 2
   coordinates%z = 2
   call init_point_probe_output(probe, coordinates, iEx, domain, test_extension, 3, 0.1_RKIND_tiempo)
   call create_point_probe_output_files(probe)

   file_time = trim(adjustl(probe%path))//'_'// &
               trim(adjustl(timeExtension))//'_'// &
               trim(adjustl(datFileExtension))

   file_freq = trim(adjustl(probe%path))//'_'// &
               trim(adjustl(frequencyExtension))//'_'// &
               trim(adjustl(datFileExtension))

   n = 10

   allocate (expectedTime(2*n, 2))
   allocate (expectedFreq(n, 2))
   !Simulate updates in probe
   do i = 1, n
      probe%timeStep(i) = real(i)
      probe%valueForTime(i) = 10.0*i
      probe%frequencySlice(i) = 0.1*i
      probe%valueForFreq(i) = 0.2*i

      expectedTime(i, 1) = real(i)
      expectedTime(i, 2) = 10.0*i

      expectedFreq(i, 1) = 0.1*i
      expectedFreq(i, 2) = 0.2*i
   end do
   probe%serializedTimeSize = n
   probe%nFreq = n
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   call flush_point_probe_output(probe)

   !Simulate new updates in probe
   do i = 1, n
      probe%timeStep(i) = real(i + 10)
      probe%valueForTime(i) = 10.0*(i + 10)
      probe%valueForFreq(i) = -0.5*i

      expectedTime(i + 10, 1) = real(i + 10)
      expectedTime(i + 10, 2) = 10.0*(i + 10)

      expectedFreq(i, 1) = 0.1*i  ! frequency file overwrites, so expectedFreq(i,1) remains 0.1*i ?
      expectedFreq(i, 2) = -0.5*i
   end do
   probe%serializedTimeSize = n
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   call flush_point_probe_output(probe)

   open (unit=probe%fileUnitTime, file=file_time, status="old", action="read")
   test_err = test_err + assert_file_content(probe%fileUnitTime, expectedTime, 2*n, 2)
   close (probe%fileUnitTime)

   open (unit=probe%fileUnitFreq, file=file_freq, status="old", action="read")
   test_err = test_err + assert_file_content(probe%fileUnitFreq, expectedFreq, n, 2)
   close (probe%fileUnitFreq)

   err = test_err

end function test_multiple_flush_point_probe

integer function test_volumic_probe_count_relevant_surfaces() bind(c) result(err)
   use output
   use mod_testOutputUtils
   use FDETYPES_TOOLS
   use mod_sggMethods
   use mod_assertionTools

   integer(kind=RKIND) :: iter
   type(media_matrices_t), target :: media
   type(media_matrices_t), pointer :: mediaPtr
   type(MediaData_t), dimension(:), allocatable, target :: simulationMaterials
   type(MediaData_t), dimension(:), pointer :: simulationMaterialsPtr
   type(limit_t), dimension(1:6), target :: sinpml_fullsize
   type(limit_t), dimension(:), pointer :: sinpml_fullsizePtr
   type(Obses_t) :: volumicProbeObservable
   type(SGGFDTDINFO) :: dummysgg
   type(sim_control_t) :: dummyControl
   type(solver_output_t), dimension(:), allocatable :: outputs

   type(MediaData_t) :: thinWireSimulationMaterial
   character(len=BUFSIZE) :: test_extension = trim(adjustl('tmp_cases/flush_point_probe'))
   integer(kind=SINGLE) :: mpidir = 3
   logical :: ThereAreWires = .false.

   real(kind=RKIND_tiempo) :: dt = 0.1_RKIND_tiempo
   integer(kind=SINGLE) :: nTimeSteps = 100_SINGLE
   real(kind=RKIND_tiempo), pointer, dimension(:) :: timeArray

   integer(kind=SINGLE) :: test_err = 0

   call sgg_init(dummysgg)
   call init_time_array(timeArray, nTimeSteps, dt)

   call sgg_set_tiempo(dummysgg, timeArray)
   call sgg_set_dt(dummysgg, dt)

   err = 1 !If test_err is not updated at the end it will be shown
   test_err = 0

   do iter = 1, 6
      sinpml_fullsize(iter) = create_limit_t(0, 8, 0, 8, 0, 8, 10, 10, 10)
   end do
   sinpml_fullsizePtr => sinpml_fullsize

   call init_simulation_material_list(simulationMaterials)
   thinWireSimulationMaterial = create_thinWire_simulation_material(size(simulationMaterials))
   call add_simulation_material(simulationMaterials, thinWireSimulationMaterial)

   call create_geometry_media(media, 0, 8, 0, 8, 0, 8)
   call assing_material_id_to_media_matrix_coordinate(media, iEy, 1, 1, 1, simulationMaterials(0)%Id)
   call assing_material_id_to_media_matrix_coordinate(media, iHz, 1, 1, 1, simulationMaterials(2)%Id)
   call assing_material_id_to_media_matrix_coordinate(media, iEx, 2, 2, 2, thinWireSimulationMaterial%Id)
   mediaPtr => media

   call sgg_set_NumMedia(dummysgg, size(simulationMaterials))
   simulationMaterialsPtr => simulationMaterials
   call sgg_set_Med(dummysgg, simulationMaterialsPtr)

   volumicProbeObservable = create_volumic_probe_observation(4, 4, 4, 6, 6, 6)
   call sgg_add_observation(dummysgg, volumicProbeObservable)

   dummyControl = create_control_flags(mpidir=mpidir, nEntradaRoot='entradaRoot', wiresflavor='holland')

   call init_outputs(dummysgg, mediaPtr, sinpml_fullsizePtr, dummyControl, outputs, ThereAreWires)

   test_err = test_err + assert_integer_equal(outputs(1)%outputID, VOLUMIC_CURRENT_PROBE_ID, 'Unexpected probe id')
   test_err = test_err + assert_integer_equal(outputs(1)%volumicCurrentProbe%columnas, 4, 'Unexpected number of columns')
   test_err = test_err + assert_string_equal(outputs(1)%volumicCurrentProbe%path, 'entradaRoot_volumicProbe_BCX_4_4_4__6_6_6', 'Unexpected path')

   call close_outputs(outputs)

   err = test_err
end function

integer function test_init_movie_probe() bind(c) result(err)
   use output
   use mod_testOutputUtils
   use FDETYPES_TOOLS
   use mod_sggMethods
   use mod_assertionTools

   ! Init inputs
   type(SGGFDTDINFO) :: dummysgg
   type(media_matrices_t), pointer :: mediaPtr
   type(MediaData_t), dimension(:), pointer :: simulationMaterialsPtr
   type(limit_t), dimension(:), pointer :: sinpml_fullsizePtr
   type(sim_control_t) :: dummyControl
   type(solver_output_t), dimension(:), allocatable :: outputs
   logical :: ThereAreWires = .false.

   !Auxiliar variables
   type(media_matrices_t), target :: media
   type(MediaData_t), dimension(:), allocatable, target :: simulationMaterials
   type(limit_t), dimension(1:6), target :: sinpml_fullsize
   type(Obses_t) :: movieObservable
   type(cell_coordinate_t) :: lowerBoundMovieProbe, upperBoundMovieProbe

   real(kind=RKIND_tiempo), pointer, dimension(:) :: timeArray
   real(kind=RKIND_tiempo) :: dt = 0.1_RKIND_tiempo
   integer(kind=SINGLE) :: nTimeSteps = 100_SINGLE

   integer(kind=SINGLE) :: expectedNumMeasurments
   integer(kind=SINGLE) :: mpidir = 3
   character(len=BUFSIZE) :: test_folder_path = trim(adjustl('tmp_cases/'))

   err = 1 !If test_err is not updated at the end it will be shown
   test_err = 0

   lowerBoundMovieProbe%x = 2
   lowerBoundMovieProbe%y = 2
   lowerBoundMovieProbe%z = 2

   upperBoundMovieProbe%x = 5
   upperBoundMovieProbe%y = 5
   upperBoundMovieProbe%z = 5

   ! Setup sgg
   call sgg_init(dummysgg)

   call init_time_array(timeArray, nTimeSteps, dt)
   call sgg_set_tiempo(dummysgg, timeArray)
   call sgg_set_dt(dummysgg, dt)

   call init_simulation_material_list(simulationMaterials)
   call sgg_set_NumMedia(dummysgg, size(simulationMaterials))
   simulationMaterialsPtr => simulationMaterials
   call sgg_set_Med(dummysgg, simulationMaterialsPtr)

   call sgg_set_Sweep(dummysgg, create_xyz_limit_array(0, 0, 0, 6, 6, 6))
   call sgg_set_SINPMLSweep(dummysgg, create_xyz_limit_array(1, 1, 1, 5, 5, 5))
   call sgg_set_NumPlaneWaves(dummysgg, 1)
   call sgg_set_Alloc(dummysgg, create_xyz_limit_array(0, 0, 0, 6, 6, 6))

   ! Define movie observation on sgg
   movieObservable = create_movie_observation(2, 2, 2, 5, 5, 5)
   call sgg_add_observation(dummysgg, movieObservable)

   call create_geometry_media(media, 0, 8, 0, 8, 0, 8)
   !----- Defining PEC surface -----!
   call assing_material_id_to_media_matrix_coordinate(media, iEy, 3, 3, 3, simulationMaterials(0)%Id)
   call assing_material_id_to_media_matrix_coordinate(media, iEy, 4, 3, 3, simulationMaterials(0)%Id)
   call assing_material_id_to_media_matrix_coordinate(media, iEy, 4, 4, 3, simulationMaterials(0)%Id)
   call assing_material_id_to_media_matrix_coordinate(media, iEy, 3, 4, 3, simulationMaterials(0)%Id)
   expectedNumMeasurments = 4_SINGLE
   !----- -------------------- -----!
   mediaPtr => media

   do iter = 1, 6
      sinpml_fullsize(iter) = create_limit_t(0, 8, 0, 8, 0, 8, 10, 10, 10)
   end do
   sinpml_fullsizePtr => sinpml_fullsize

   dummyControl = create_control_flags(nEntradaRoot=test_folder_path, mpidir=mpidir)

   call init_outputs(dummysgg, mediaPtr, sinpml_fullsizePtr, dummyControl, outputs, ThereAreWires)

   test_err = test_err + assert_integer_equal(outputs(1)%outputID, MOVIE_PROBE_ID, 'Unexpected probe id')
   test_err = test_err + assert_integer_equal(outputs(1)%movieProbe%columnas, 4, 'Unexpected number of columns')
   test_err = test_err + assert_integer_equal(outputs(1)%movieProbe%nMeasuredElements, expectedNumMeasurments, 'Unexpected number of measurments')
   test_err = test_err + assert_integer_equal(size(outputs(1)%movieProbe%xValueForTime), expectedNumMeasurments * BuffObse, 'Unexpected allocation size')
   if (size(outputs(1)%movieProbe%timeStep) /= BuffObse) then
      test_err = 1
   end if

   call close_outputs(outputs)

   err = test_err
end function

integer function test_update_movie_probe() bind(c) result(err)
   use output
   use mod_testOutputUtils
   use FDETYPES_TOOLS
   use mod_sggMethods
   use mod_assertionTools

   ! Init inputs
   type(SGGFDTDINFO) :: dummysgg
   type(media_matrices_t), pointer :: mediaPtr
   type(MediaData_t), dimension(:), pointer :: simulationMaterialsPtr
   type(limit_t), dimension(:), pointer :: sinpml_fullsizePtr
   type(sim_control_t) :: dummyControl
   type(solver_output_t), dimension(:), allocatable :: outputs
   logical :: ThereAreWires = .false.

   !Auxiliar variables
   type(media_matrices_t), target :: media
   type(MediaData_t), dimension(:), allocatable, target :: simulationMaterials
   type(limit_t), dimension(1:6), target :: sinpml_fullsize
   type(Obses_t) :: movieObservable
   type(cell_coordinate_t) :: lowerBoundMovieProbe, upperBoundMovieProbe

   real(kind=RKIND_tiempo), pointer, dimension(:) :: timeArray
   real(kind=RKIND_tiempo) :: dt = 0.1_RKIND_tiempo
   integer(kind=SINGLE) :: nTimeSteps = 100_SINGLE

   integer(kind=SINGLE) :: expectedNumMeasurments
   integer(kind=SINGLE) :: mpidir = 3
   character(len=BUFSIZE) :: test_folder_path = trim(adjustl('tmp_cases/'))

   !DummyField required variables
   type(dummyFields_t), target :: dummyfields
   type(fields_reference_t) :: fields

   err = 1 !If test_err is not updated at the end it will be shown
   test_err = 0

   lowerBoundMovieProbe%x = 2
   lowerBoundMovieProbe%y = 2
   lowerBoundMovieProbe%z = 2

   upperBoundMovieProbe%x = 5
   upperBoundMovieProbe%y = 5
   upperBoundMovieProbe%z = 5

   ! Setup sgg
   call sgg_init(dummysgg)
   call init_time_array(timeArray, nTimeSteps, dt)
   call sgg_set_tiempo(dummysgg, timeArray)
   call sgg_set_dt(dummysgg, dt)
   call init_simulation_material_list(simulationMaterials)
   call sgg_set_NumMedia(dummysgg, size(simulationMaterials))
   simulationMaterialsPtr => simulationMaterials
   call sgg_set_Med(dummysgg, simulationMaterialsPtr)
   call sgg_set_Sweep(dummysgg, create_xyz_limit_array(0, 0, 0, 6, 6, 6))
   call sgg_set_SINPMLSweep(dummysgg, create_xyz_limit_array(1, 1, 1, 5, 5, 5))
   call sgg_set_NumPlaneWaves(dummysgg, 1)
   call sgg_set_Alloc(dummysgg, create_xyz_limit_array(0, 0, 0, 6, 6, 6))
   ! Define movie observation on sgg
   movieObservable = create_movie_observation(2, 2, 2, 5, 5, 5)
   call sgg_add_observation(dummysgg, movieObservable)
   call create_geometry_media(media, 0, 8, 0, 8, 0, 8)
   !----- Defining PEC surface -----!
   call assing_material_id_to_media_matrix_coordinate(media, iEy, 3, 3, 3, simulationMaterials(0)%Id)
   call assing_material_id_to_media_matrix_coordinate(media, iEy, 4, 3, 3, simulationMaterials(0)%Id)
   call assing_material_id_to_media_matrix_coordinate(media, iEy, 4, 4, 3, simulationMaterials(0)%Id)
   call assing_material_id_to_media_matrix_coordinate(media, iEy, 3, 4, 3, simulationMaterials(0)%Id)
   expectedNumMeasurments = 4_SINGLE
   !----- -------------------- -----!
   mediaPtr => media
   do iter = 1, 6
      sinpml_fullsize(iter) = create_limit_t(0, 8, 0, 8, 0, 8, 10, 10, 10)
   end do
   sinpml_fullsizePtr => sinpml_fullsize

   dummyControl = create_control_flags(nEntradaRoot=test_folder_path, mpidir=mpidir)

   call init_outputs(dummysgg, mediaPtr, sinpml_fullsizePtr, dummyControl, outputs, ThereAreWires)

   ! Set dummy field status

   call create_dummy_fields(dummyfields, 1, 5, 0.1_RKIND)
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

   dummyfields%Hx(3, 3, 3) = 2.0_RKIND
   dummyfields%Hy(3, 3, 3) = 5.0_RKIND
   dummyfields%Hz(3, 3, 3) = 4.0_RKIND

   call update_outputs(outputs, mediaPtr, simulationMaterialsPtr, sinpml_fullsizePtr, dummyControl, 0.5_RKIND_tiempo, fields)

   test_err = test_err + assert_integer_equal(outputs(1)%outputID, MOVIE_PROBE_ID, 'Unexpected probe id')
   test_err = test_err + assert_integer_equal(outputs(1)%movieProbe%columnas, 4, 'Unexpected number of columns')
   test_err = test_err + assert_integer_equal(outputs(1)%movieProbe%nMeasuredElements, expectedNumMeasurments, 'Unexpected number of measurments')
   test_err = test_err + assert_integer_equal(size(outputs(1)%movieProbe%xValueForTime), expectedNumMeasurments * BuffObse, 'Unexpected allocation size')
   test_err = test_err + assert_real_equal(outputs(1)%movieProbe%yValueForTime(1, 1), 0.2_RKIND, 0.00001_RKIND, 'Value error')
   test_err = test_err + assert_real_equal(outputs(1)%movieProbe%yValueForTime(1, 2), 0.0_RKIND, 0.00001_RKIND, 'Value error')
   test_err = test_err + assert_real_equal(outputs(1)%movieProbe%yValueForTime(1, 3), 0.2_RKIND, 0.00001_RKIND, 'Value error')
   test_err = test_err + assert_real_equal(outputs(1)%movieProbe%yValueForTime(1, 4), 0.0_RKIND, 0.00001_RKIND, 'Value error')
   if (size(outputs(1)%movieProbe%timeStep) /= BuffObse) then
      test_err = 1
   end if

   call close_outputs(outputs)

   err = test_err
end function

integer function test_flush_movie_probe() bind(c) result(err)
   use output
   use mod_testOutputUtils
   use FDETYPES_TOOLS
   use mod_sggMethods
   use mod_assertionTools

   ! Init inputs
   type(SGGFDTDINFO) :: dummysgg
   type(media_matrices_t), pointer :: mediaPtr
   type(MediaData_t), dimension(:), pointer :: simulationMaterialsPtr
   type(limit_t), dimension(:), pointer :: sinpml_fullsizePtr
   type(sim_control_t) :: dummyControl
   type(solver_output_t), dimension(:), allocatable :: outputs
   logical :: ThereAreWires = .false.

   !Auxiliar variables
   type(media_matrices_t), target :: media
   type(MediaData_t), dimension(:), allocatable, target :: simulationMaterials
   type(limit_t), dimension(1:6), target :: sinpml_fullsize
   type(Obses_t) :: movieObservable
   type(cell_coordinate_t) :: lowerBoundMovieProbe, upperBoundMovieProbe

   real(kind=RKIND_tiempo), pointer, dimension(:) :: timeArray
   real(kind=RKIND_tiempo) :: dt = 0.1_RKIND_tiempo
   integer(kind=SINGLE) :: nTimeSteps = 100_SINGLE

   integer(kind=SINGLE) :: expectedNumMeasurments
   integer(kind=SINGLE) :: mpidir = 3
   character(len=BUFSIZE) :: test_folder_path = trim(adjustl('tmp_cases/'))
   character(len=BUFSIZE) :: expectedPath

   err = 1 !If test_err is not updated at the end it will be shown
   test_err = 0

   lowerBoundMovieProbe%x = 2
   lowerBoundMovieProbe%y = 2
   lowerBoundMovieProbe%z = 2

   upperBoundMovieProbe%x = 5
   upperBoundMovieProbe%y = 5
   upperBoundMovieProbe%z = 5

   ! Setup sgg
   call sgg_init(dummysgg)
   call init_time_array(timeArray, nTimeSteps, dt)
   call sgg_set_tiempo(dummysgg, timeArray)
   call sgg_set_dt(dummysgg, dt)
   call init_simulation_material_list(simulationMaterials)
   call sgg_set_NumMedia(dummysgg, size(simulationMaterials))
   simulationMaterialsPtr => simulationMaterials
   call sgg_set_Med(dummysgg, simulationMaterialsPtr)
   call sgg_set_Sweep(dummysgg, create_xyz_limit_array(0, 0, 0, 6, 6, 6))
   call sgg_set_SINPMLSweep(dummysgg, create_xyz_limit_array(1, 1, 1, 5, 5, 5))
   call sgg_set_NumPlaneWaves(dummysgg, 1)
   call sgg_set_Alloc(dummysgg, create_xyz_limit_array(0, 0, 0, 6, 6, 6))
   ! Define movie observation on sgg
   movieObservable = create_movie_observation(2, 2, 2, 5, 5, 5)
   call sgg_add_observation(dummysgg, movieObservable)
   call create_geometry_media(media, 0, 8, 0, 8, 0, 8)
   !----- Defining PEC surface -----!
   call assing_material_id_to_media_matrix_coordinate(media, iEy, 3, 3, 3, simulationMaterials(0)%Id)
   call assing_material_id_to_media_matrix_coordinate(media, iEy, 4, 3, 3, simulationMaterials(0)%Id)
   call assing_material_id_to_media_matrix_coordinate(media, iEy, 4, 4, 3, simulationMaterials(0)%Id)
   call assing_material_id_to_media_matrix_coordinate(media, iEy, 3, 4, 3, simulationMaterials(0)%Id)
   expectedNumMeasurments = 4_SINGLE
   !----- -------------------- -----!
   mediaPtr => media
   do iter = 1, 6
      sinpml_fullsize(iter) = create_limit_t(0, 8, 0, 8, 0, 8, 10, 10, 10)
   end do
   sinpml_fullsizePtr => sinpml_fullsize

   dummyControl = create_control_flags(nEntradaRoot=test_folder_path, mpidir=mpidir)

   call init_outputs(dummysgg, mediaPtr, sinpml_fullsizePtr, dummyControl, outputs, ThereAreWires)

   !Dummy first update
   outputs(1)%movieProbe%serializedTimeSize = 1
   outputs(1)%movieProbe%timeStep(1) = 0.5_RKIND_tiempo

   outputs(1)%movieProbe%xValueForTime(1, 1) = 0.0_RKIND
   outputs(1)%movieProbe%xValueForTime(1, 2) = 0.0_RKIND
   outputs(1)%movieProbe%xValueForTime(1, 3) = 0.0_RKIND
   outputs(1)%movieProbe%xValueForTime(1, 4) = 0.0_RKIND

   outputs(1)%movieProbe%yValueForTime(1, 1) = 0.1_RKIND
   outputs(1)%movieProbe%yValueForTime(1, 2) = 0.2_RKIND
   outputs(1)%movieProbe%yValueForTime(1, 3) = 0.3_RKIND
   outputs(1)%movieProbe%yValueForTime(1, 4) = 0.4_RKIND

   outputs(1)%movieProbe%zValueForTime(1, 1) = 0.0_RKIND
   outputs(1)%movieProbe%zValueForTime(1, 2) = 0.0_RKIND
   outputs(1)%movieProbe%zValueForTime(1, 3) = 0.0_RKIND
   outputs(1)%movieProbe%zValueForTime(1, 4) = 0.0_RKIND

   !Dummy second update
   outputs(1)%movieProbe%serializedTimeSize = 2
   outputs(1)%movieProbe%timeStep(2) = 1.0_RKIND_tiempo

   outputs(1)%movieProbe%xValueForTime(2, 1) = 0.0_RKIND
   outputs(1)%movieProbe%xValueForTime(2, 2) = 0.0_RKIND
   outputs(1)%movieProbe%xValueForTime(2, 3) = 0.0_RKIND
   outputs(1)%movieProbe%xValueForTime(2, 4) = 0.0_RKIND

   outputs(1)%movieProbe%yValueForTime(2, 1) = 0.11_RKIND
   outputs(1)%movieProbe%yValueForTime(2, 2) = 0.22_RKIND
   outputs(1)%movieProbe%yValueForTime(2, 3) = 0.33_RKIND
   outputs(1)%movieProbe%yValueForTime(2, 4) = 0.44_RKIND

   outputs(1)%movieProbe%zValueForTime(2, 1) = 0.0_RKIND
   outputs(1)%movieProbe%zValueForTime(2, 2) = 0.0_RKIND
   outputs(1)%movieProbe%zValueForTime(2, 3) = 0.0_RKIND
   outputs(1)%movieProbe%zValueForTime(2, 4) = 0.0_RKIND

   call flush_outputs(outputs)

   expectedPath = trim(adjustl(outputs(1)%movieProbe%path))//'_ts'//'0001'//'.vtu'
   test_err = test_err + assert_file_exists(expectedPath)

   expectedPath = trim(adjustl(outputs(1)%movieProbe%path))//'_ts'//'0002'//'.vtu'
   test_err = test_err + assert_file_exists(expectedPath)

   call close_outputs(outputs)

   expectedPath = trim(adjustl(outputs(1)%movieProbe%path))//'.pvd'
   test_err = test_err + assert_file_exists(expectedPath)

   err = test_err
end function

integer function test_init_frequency_slice_probe() bind(c) result(err)
   use output
   use mod_testOutputUtils
   use FDETYPES_TOOLS
   use mod_sggMethods
   use mod_assertionTools

   ! Init inputs
   type(SGGFDTDINFO) :: dummysgg
   type(media_matrices_t), pointer :: mediaPtr
   type(MediaData_t), dimension(:), pointer :: simulationMaterialsPtr
   type(limit_t), dimension(:), pointer :: sinpml_fullsizePtr
   type(sim_control_t) :: dummyControl
   type(solver_output_t), dimension(:), allocatable :: outputs
   logical :: ThereAreWires = .false.

   !Auxiliar variables
   type(media_matrices_t), target :: media
   type(MediaData_t), dimension(:), allocatable, target :: simulationMaterials
   type(limit_t), dimension(1:6), target :: sinpml_fullsize
   type(Obses_t) :: frequencySliceObservation

   real(kind=RKIND_tiempo), pointer, dimension(:) :: timeArray
   real(kind=RKIND_tiempo) :: dt = 0.1_RKIND_tiempo
   integer(kind=SINGLE) :: nTimeSteps = 100_SINGLE

   integer(kind=SINGLE) :: expectedNumMeasurments
   integer(kind=SINGLE) :: expectedTotalFrequnecies
   integer(kind=SINGLE) :: mpidir = 3
   character(len=BUFSIZE) :: test_folder_path = trim(adjustl('tmp_cases/'))

   err = 1 !If test_err is not updated at the end it will be shown
   test_err = 0

   ! Setup sgg
   call sgg_init(dummysgg)

   call init_time_array(timeArray, nTimeSteps, dt)
   call sgg_set_tiempo(dummysgg, timeArray)
   call sgg_set_dt(dummysgg, dt)

   call init_simulation_material_list(simulationMaterials)
   call sgg_set_NumMedia(dummysgg, size(simulationMaterials))
   simulationMaterialsPtr => simulationMaterials
   call sgg_set_Med(dummysgg, simulationMaterialsPtr)

   call sgg_set_Sweep(dummysgg, create_xyz_limit_array(0, 0, 0, 6, 6, 6))
   call sgg_set_SINPMLSweep(dummysgg, create_xyz_limit_array(1, 1, 1, 5, 5, 5))
   call sgg_set_NumPlaneWaves(dummysgg, 1)
   call sgg_set_Alloc(dummysgg, create_xyz_limit_array(0, 0, 0, 6, 6, 6))

   ! Define movie observation on sgg
   frequencySliceObservation = create_frequency_slice_observation(2, 2, 2, 5, 5, 5)
   call sgg_add_observation(dummysgg, frequencySliceObservation)
   expectedTotalFrequnecies = 6_SINGLE

   call create_geometry_media(media, 0, 8, 0, 8, 0, 8)
   !----- Defining PEC surface -----!
   call assing_material_id_to_media_matrix_coordinate(media, iEy, 3, 3, 3, simulationMaterials(0)%Id)
   call assing_material_id_to_media_matrix_coordinate(media, iEy, 4, 3, 3, simulationMaterials(0)%Id)
   call assing_material_id_to_media_matrix_coordinate(media, iEy, 4, 4, 3, simulationMaterials(0)%Id)
   call assing_material_id_to_media_matrix_coordinate(media, iEy, 3, 4, 3, simulationMaterials(0)%Id)
   expectedNumMeasurments = 4_SINGLE
   !----- -------------------- -----!
   mediaPtr => media

   do iter = 1, 6
      sinpml_fullsize(iter) = create_limit_t(0, 8, 0, 8, 0, 8, 10, 10, 10)
   end do
   sinpml_fullsizePtr => sinpml_fullsize

   dummyControl = create_control_flags(nEntradaRoot=test_folder_path, mpidir=mpidir)

   call init_outputs(dummysgg, mediaPtr, sinpml_fullsizePtr, dummyControl, outputs, ThereAreWires)

   test_err = test_err + assert_integer_equal(outputs(1)%outputID, FREQUENCY_SLICE_PROBE_ID, 'Unexpected probe id')
   test_err = test_err + assert_integer_equal(outputs(1)%frequencySliceProbe%columnas, 4, 'Unexpected number of columns')
   test_err = test_err + assert_integer_equal(outputs(1)%frequencySliceProbe%nMeasuredElements, expectedNumMeasurments, 'Unexpected number of measurments')
   test_err = test_err + assert_integer_equal(size(outputs(1)%frequencySliceProbe%xValueForFreq), expectedNumMeasurments * expectedTotalFrequnecies, 'Unexpected allocation size')
   if (size(outputs(1)%frequencySliceProbe%frequencySlice) /= expectedTotalFrequnecies) then
      test_err = 1
   end if

   call close_outputs(outputs)

   err = test_err
end function

integer function test_update_frequency_slice_probe() bind(c) result(err)
   use output
   use mod_testOutputUtils
   use FDETYPES_TOOLS
   use mod_sggMethods
   use mod_assertionTools

   ! Init inputs
   type(SGGFDTDINFO) :: dummysgg
   type(media_matrices_t), pointer :: mediaPtr
   type(MediaData_t), dimension(:), pointer :: simulationMaterialsPtr
   type(limit_t), dimension(:), pointer :: sinpml_fullsizePtr
   type(sim_control_t) :: dummyControl
   type(solver_output_t), dimension(:), allocatable :: outputs
   logical :: ThereAreWires = .false.

   !Auxiliar variables
   type(media_matrices_t), target :: media
   type(MediaData_t), dimension(:), allocatable, target :: simulationMaterials
   type(limit_t), dimension(1:6), target :: sinpml_fullsize
   type(Obses_t) :: frequencySliceObservation

   real(kind=RKIND_tiempo), pointer, dimension(:) :: timeArray
   real(kind=RKIND_tiempo) :: dt = 0.1_RKIND_tiempo
   integer(kind=SINGLE) :: nTimeSteps = 100_SINGLE

   integer(kind=SINGLE) :: expectedNumMeasurments
   integer(kind=SINGLE) :: mpidir = 3
   character(len=BUFSIZE) :: test_folder_path = trim(adjustl('tmp_cases/'))

   !DummyField required variables
   type(dummyFields_t), target :: dummyfields
   type(fields_reference_t) :: fields

   err = 1 !If test_err is not updated at the end it will be shown
   test_err = 0

   ! Setup sgg
   call sgg_init(dummysgg)
   call init_time_array(timeArray, nTimeSteps, dt)
   call sgg_set_tiempo(dummysgg, timeArray)
   call sgg_set_dt(dummysgg, dt)
   call init_simulation_material_list(simulationMaterials)
   call sgg_set_NumMedia(dummysgg, size(simulationMaterials))
   simulationMaterialsPtr => simulationMaterials
   call sgg_set_Med(dummysgg, simulationMaterialsPtr)
   call sgg_set_Sweep(dummysgg, create_xyz_limit_array(0, 0, 0, 6, 6, 6))
   call sgg_set_SINPMLSweep(dummysgg, create_xyz_limit_array(1, 1, 1, 5, 5, 5))
   call sgg_set_NumPlaneWaves(dummysgg, 1)
   call sgg_set_Alloc(dummysgg, create_xyz_limit_array(0, 0, 0, 6, 6, 6))
   ! Define movie observation on sgg
   frequencySliceObservation = create_frequency_slice_observation(2, 2, 2, 5, 5, 5)
   call sgg_add_observation(dummysgg, frequencySliceObservation)
   call create_geometry_media(media, 0, 8, 0, 8, 0, 8)
   !----- Defining PEC surface -----!
   call assing_material_id_to_media_matrix_coordinate(media, iEy, 3, 3, 3, simulationMaterials(0)%Id)
   call assing_material_id_to_media_matrix_coordinate(media, iEy, 4, 3, 3, simulationMaterials(0)%Id)
   call assing_material_id_to_media_matrix_coordinate(media, iEy, 4, 4, 3, simulationMaterials(0)%Id)
   call assing_material_id_to_media_matrix_coordinate(media, iEy, 3, 4, 3, simulationMaterials(0)%Id)
   expectedNumMeasurments = 4_SINGLE
   !----- -------------------- -----!
   mediaPtr => media
   do iter = 1, 6
      sinpml_fullsize(iter) = create_limit_t(0, 8, 0, 8, 0, 8, 10, 10, 10)
   end do
   sinpml_fullsizePtr => sinpml_fullsize

   dummyControl = create_control_flags(nEntradaRoot=test_folder_path, mpidir=mpidir)

   call init_outputs(dummysgg, mediaPtr, sinpml_fullsizePtr, dummyControl, outputs, ThereAreWires)

   ! Set dummy field status

   call create_dummy_fields(dummyfields, 1, 5, 0.1_RKIND)
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

   dummyfields%Hx(3, 3, 3) = 2.0_RKIND
   dummyfields%Hy(3, 3, 3) = 5.0_RKIND
   dummyfields%Hz(3, 3, 3) = 4.0_RKIND

   call update_outputs(outputs, mediaPtr, simulationMaterialsPtr, sinpml_fullsizePtr, dummyControl, 0.5_RKIND_tiempo, fields)

   test_err = test_err + assert_integer_equal(outputs(1)%outputID, MOVIE_PROBE_ID, 'Unexpected probe id')
   test_err = test_err + assert_integer_equal(outputs(1)%frequencySliceProbe%columnas, 4, 'Unexpected number of columns')
   test_err = test_err + assert_integer_equal(outputs(1)%frequencySliceProbe%nMeasuredElements, expectedNumMeasurments, 'Unexpected number of measurments')
   test_err = test_err + assert_integer_equal(size(outputs(1)%frequencySliceProbe%frequencySlice), expectedNumMeasurments * BuffObse, 'Unexpected allocation size')
   test_err = test_err + assert_complex_equal(outputs(1)%frequencySliceProbe%yValueForFreq(3, 1), (0.2_CKIND ,0.2_CKIND), 0.00001_RKIND, 'Value error')
   test_err = test_err + assert_complex_equal(outputs(1)%frequencySliceProbe%yValueForFreq(3, 2), (0.2_CKIND ,0.2_CKIND), 0.00001_RKIND, 'Value error')
   test_err = test_err + assert_complex_equal(outputs(1)%frequencySliceProbe%yValueForFreq(3, 3), (0.2_CKIND ,0.2_CKIND), 0.00001_RKIND, 'Value error')
   test_err = test_err + assert_complex_equal(outputs(1)%frequencySliceProbe%yValueForFreq(3, 4), (0.2_CKIND ,0.2_CKIND), 0.00001_RKIND, 'Value error')
   test_err = test_err + assert_complex_equal(outputs(1)%frequencySliceProbe%yValueForFreq(3, 5), (0.2_CKIND ,0.2_CKIND), 0.00001_RKIND, 'Value error')

   call close_outputs(outputs)

   err = test_err
end function

integer function test_flush_frequency_slice_probe() bind(c) result(err)
   use output
   use mod_testOutputUtils
   use FDETYPES_TOOLS
   use mod_sggMethods
   use mod_assertionTools

   ! Init inputs
   type(SGGFDTDINFO) :: dummysgg
   type(media_matrices_t), pointer :: mediaPtr
   type(MediaData_t), dimension(:), pointer :: simulationMaterialsPtr
   type(limit_t), dimension(:), pointer :: sinpml_fullsizePtr
   type(sim_control_t) :: dummyControl
   type(solver_output_t), dimension(:), allocatable :: outputs
   logical :: ThereAreWires = .false.

   !Auxiliar variables
   type(media_matrices_t), target :: media
   type(MediaData_t), dimension(:), allocatable, target :: simulationMaterials
   type(limit_t), dimension(1:6), target :: sinpml_fullsize
   type(Obses_t) :: movieObservable
   type(cell_coordinate_t) :: lowerBoundMovieProbe, upperBoundMovieProbe

   real(kind=RKIND_tiempo), pointer, dimension(:) :: timeArray
   real(kind=RKIND_tiempo) :: dt = 0.1_RKIND_tiempo
   integer(kind=SINGLE) :: nTimeSteps = 100_SINGLE

   integer(kind=SINGLE) :: expectedNumMeasurments
   integer(kind=SINGLE) :: mpidir = 3
   character(len=BUFSIZE) :: test_folder_path = trim(adjustl('tmp_cases/'))
   character(len=BUFSIZE) :: expectedPath

   err = 1 !If test_err is not updated at the end it will be shown
   test_err = 0

   lowerBoundMovieProbe%x = 2
   lowerBoundMovieProbe%y = 2
   lowerBoundMovieProbe%z = 2

   upperBoundMovieProbe%x = 5
   upperBoundMovieProbe%y = 5
   upperBoundMovieProbe%z = 5

   ! Setup sgg
   call sgg_init(dummysgg)
   call init_time_array(timeArray, nTimeSteps, dt)
   call sgg_set_tiempo(dummysgg, timeArray)
   call sgg_set_dt(dummysgg, dt)
   call init_simulation_material_list(simulationMaterials)
   call sgg_set_NumMedia(dummysgg, size(simulationMaterials))
   simulationMaterialsPtr => simulationMaterials
   call sgg_set_Med(dummysgg, simulationMaterialsPtr)
   call sgg_set_Sweep(dummysgg, create_xyz_limit_array(0, 0, 0, 6, 6, 6))
   call sgg_set_SINPMLSweep(dummysgg, create_xyz_limit_array(1, 1, 1, 5, 5, 5))
   call sgg_set_NumPlaneWaves(dummysgg, 1)
   call sgg_set_Alloc(dummysgg, create_xyz_limit_array(0, 0, 0, 6, 6, 6))
   ! Define movie observation on sgg
   movieObservable = create_movie_observation(2, 2, 2, 5, 5, 5)
   call sgg_add_observation(dummysgg, movieObservable)
   call create_geometry_media(media, 0, 8, 0, 8, 0, 8)
   !----- Defining PEC surface -----!
   call assing_material_id_to_media_matrix_coordinate(media, iEy, 3, 3, 3, simulationMaterials(0)%Id)
   call assing_material_id_to_media_matrix_coordinate(media, iEy, 4, 3, 3, simulationMaterials(0)%Id)
   call assing_material_id_to_media_matrix_coordinate(media, iEy, 4, 4, 3, simulationMaterials(0)%Id)
   call assing_material_id_to_media_matrix_coordinate(media, iEy, 3, 4, 3, simulationMaterials(0)%Id)
   expectedNumMeasurments = 4_SINGLE
   !----- -------------------- -----!
   mediaPtr => media
   do iter = 1, 6
      sinpml_fullsize(iter) = create_limit_t(0, 8, 0, 8, 0, 8, 10, 10, 10)
   end do
   sinpml_fullsizePtr => sinpml_fullsize

   dummyControl = create_control_flags(nEntradaRoot=test_folder_path, mpidir=mpidir)

   call init_outputs(dummysgg, mediaPtr, sinpml_fullsizePtr, dummyControl, outputs, ThereAreWires)

   !Dummy first update
   outputs(1)%movieProbe%serializedTimeSize = 1
   outputs(1)%movieProbe%timeStep(1) = 0.5_RKIND_tiempo

   outputs(1)%movieProbe%xValueForTime(1, 1) = 0.0_RKIND
   outputs(1)%movieProbe%xValueForTime(1, 2) = 0.0_RKIND
   outputs(1)%movieProbe%xValueForTime(1, 3) = 0.0_RKIND
   outputs(1)%movieProbe%xValueForTime(1, 4) = 0.0_RKIND

   outputs(1)%movieProbe%yValueForTime(1, 1) = 0.1_RKIND
   outputs(1)%movieProbe%yValueForTime(1, 2) = 0.2_RKIND
   outputs(1)%movieProbe%yValueForTime(1, 3) = 0.3_RKIND
   outputs(1)%movieProbe%yValueForTime(1, 4) = 0.4_RKIND

   outputs(1)%movieProbe%zValueForTime(1, 1) = 0.0_RKIND
   outputs(1)%movieProbe%zValueForTime(1, 2) = 0.0_RKIND
   outputs(1)%movieProbe%zValueForTime(1, 3) = 0.0_RKIND
   outputs(1)%movieProbe%zValueForTime(1, 4) = 0.0_RKIND

   !Dummy second update
   outputs(1)%movieProbe%serializedTimeSize = 2
   outputs(1)%movieProbe%timeStep(2) = 1.0_RKIND_tiempo

   outputs(1)%movieProbe%xValueForTime(2, 1) = 0.0_RKIND
   outputs(1)%movieProbe%xValueForTime(2, 2) = 0.0_RKIND
   outputs(1)%movieProbe%xValueForTime(2, 3) = 0.0_RKIND
   outputs(1)%movieProbe%xValueForTime(2, 4) = 0.0_RKIND

   outputs(1)%movieProbe%yValueForTime(2, 1) = 0.11_RKIND
   outputs(1)%movieProbe%yValueForTime(2, 2) = 0.22_RKIND
   outputs(1)%movieProbe%yValueForTime(2, 3) = 0.33_RKIND
   outputs(1)%movieProbe%yValueForTime(2, 4) = 0.44_RKIND

   outputs(1)%movieProbe%zValueForTime(2, 1) = 0.0_RKIND
   outputs(1)%movieProbe%zValueForTime(2, 2) = 0.0_RKIND
   outputs(1)%movieProbe%zValueForTime(2, 3) = 0.0_RKIND
   outputs(1)%movieProbe%zValueForTime(2, 4) = 0.0_RKIND

   call flush_outputs(outputs)

   expectedPath = trim(adjustl(outputs(1)%movieProbe%path))//'_ts'//'0001'//'.vtu'
   test_err = test_err + assert_file_exists(expectedPath)

   expectedPath = trim(adjustl(outputs(1)%movieProbe%path))//'_ts'//'0002'//'.vtu'
   test_err = test_err + assert_file_exists(expectedPath)

   call close_outputs(outputs)

   expectedPath = trim(adjustl(outputs(1)%movieProbe%path))//'.pvd'
   test_err = test_err + assert_file_exists(expectedPath)

   err = test_err
end function
