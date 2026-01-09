integer function test_init_point_probe() bind(c) result(err)
   use FDETYPES
   use FDETYPES_TOOLS
   use output
   use mod_testOutputUtils
   use mod_sggMethods
   use mod_assertionTools

   type(SGGFDTDINFO)              :: sgg
   type(sim_control_t)            :: control
   type(bounds_t)                 :: bounds
   type(media_matrices_t)         :: media
   type(limit_t), allocatable     :: sinpml(:)
   type(Obses_t)                  :: probe
   type(solver_output_t), pointer :: outputs(:)
   type(MediaData_t), allocatable, target :: materials(:)
   type(MediaData_t), pointer     :: materialsPtr(:)

   real(kind=RKIND_tiempo), pointer :: timeArray(:)
   real(kind=RKIND_tiempo)          :: dt = 0.1_RKIND_tiempo
   integer(kind=SINGLE)             :: nSteps = 100_SINGLE
   logical                          :: outputRequested, hasWires = .false.
   integer(kind=SINGLE)             :: test_err = 0

   call sgg_init(sgg)
   call init_time_array(timeArray, nSteps, dt)
   call sgg_set_tiempo(sgg, timeArray)
   call sgg_set_dt(sgg, dt)

   call init_simulation_material_list(materials)
   materialsPtr => materials
   call sgg_set_Med(sgg, materialsPtr)

   probe = create_point_probe_observation(4, 4, 4)
   call sgg_add_observation(sgg, probe)

   control = create_control_flags(mpidir=3, nEntradaRoot='entradaRoot', wiresflavor='holland')

   call init_outputs(sgg, media, sinpml, bounds, control, outputRequested, hasWires)

   outputs => GetOutputs()
   test_err = test_err + assert_true(outputRequested, 'Valid probes not found')

   test_err = test_err + assert_integer_equal(outputs(1)%outputID, POINT_PROBE_ID, 'Unexpected probe id')
   test_err = test_err + assert_string_equal(outputs(1)%pointProbe%path, &
                                             'entradaRoot_poinProbe_Ex_4_4_4', 'Unexpected path')

   call close_outputs()
   deallocate (sgg%Observation, outputs)

   err = test_err
end function

integer function test_update_point_probe() bind(c) result(err)
   use FDETYPES
   use FDETYPES_TOOLS
   use output
   use outputTypes
   use mod_testOutputUtils
   use mod_sggMethods
   use mod_assertionTools

   type(SGGFDTDINFO)              :: sgg
   type(sim_control_t)            :: control
   type(bounds_t)                 :: bounds
   type(media_matrices_t)         :: media
   type(limit_t), allocatable     :: sinpml(:)
   type(Obses_t)                  :: probe
   type(solver_output_t), pointer :: outputs(:)
   type(MediaData_t), allocatable, target :: materials(:)
   type(MediaData_t), pointer     :: materialsPtr(:)

   type(dummyFields_t), target    :: dummyFields
   type(fields_reference_t)       :: fields

   real(kind=RKIND_tiempo), pointer :: timeArray(:)
   real(kind=RKIND_tiempo)          :: dt = 0.1_RKIND_tiempo
   integer(kind=SINGLE)             :: nSteps = 100_SINGLE
   logical                          :: outputRequested, hasWires = .false.
   integer(kind=SINGLE)             :: test_err = 0

   call sgg_init(sgg)
   call init_time_array(timeArray, nSteps, dt)
   call sgg_set_tiempo(sgg, timeArray)
   call sgg_set_dt(sgg, dt)

   probe = create_point_probe_observation(4, 4, 4)
   call sgg_add_observation(sgg, probe)

   call init_simulation_material_list(materials)
   materialsPtr => materials
   call sgg_set_Med(sgg, materialsPtr)

   control = create_control_flags(mpidir=3, nEntradaRoot='entradaRoot', wiresflavor='holland')

   call init_outputs(sgg, media, sinpml, bounds, control, outputRequested, hasWires)

   call create_dummy_fields(dummyFields, 1, 10, 0.01_RKIND)

   fields%E%x => dummyFields%Ex
   fields%E%y => dummyFields%Ey
   fields%E%z => dummyFields%Ez
   fields%E%deltax => dummyFields%dxe
   fields%E%deltaY => dummyFields%dye
   fields%E%deltaZ => dummyFields%dze
   fields%H%x => dummyFields%Hx
   fields%H%y => dummyFields%Hy
   fields%H%z => dummyFields%Hz
   fields%H%deltax => dummyFields%dxh
   fields%H%deltaY => dummyFields%dyh
   fields%H%deltaZ => dummyFields%dzh

   dummyFields%Ex(4, 4, 4) = 5.0_RKIND
   call update_outputs(control, sgg%tiempo, 1_SINGLE, fields)

   outputs => GetOutputs()

   test_err = test_err + assert_real_time_equal(outputs(1)%pointProbe%timeStep(1), 0.0_RKIND_tiempo, 1e-5_RKIND_tiempo, 'Unexpected timestep 1')
   test_err = test_err + assert_real_equal(outputs(1)%pointProbe%valueForTime(1), 5.0_RKIND, 1e-5_RKIND, 'Unexpected field 1')

   dummyFields%Ex(4, 4, 4) = -4.0_RKIND
   call update_outputs(control, sgg%tiempo, 2_SINGLE, fields)

   test_err = test_err + assert_real_time_equal(outputs(1)%pointProbe%timeStep(2), 0.1_RKIND_tiempo, 1e-5_RKIND_tiempo, 'Unexpected timestep 2')
   test_err = test_err + assert_real_equal(outputs(1)%pointProbe%valueForTime(2), -4.0_RKIND, 1e-5_RKIND, 'Unexpected field 2')

   call close_outputs()

   err = test_err
end function

integer function test_flush_point_probe() bind(c) result(err)
   use output
   use mod_pointProbeOutput
   use mod_domain
   use mod_testOutputUtils
   use mod_assertionTools

   type(point_probe_output_t) :: probe
   type(domain_t)             :: domain
   type(cell_coordinate_t)    :: coordinates

   character(len=BUFSIZE) :: file_time, file_freq
   character(len=27)      :: test_extension

   integer :: n, i
   integer :: test_err = 0

   err = 1
   test_extension = 'tmp_cases/flush_point_probe'

   domain = domain_t( &
            0.0_RKIND_tiempo, 10.0_RKIND_tiempo, 0.1_RKIND_tiempo, &
            10.0_RKIND, 100.0_RKIND, 10, .false.)

   coordinates%x = 2
   coordinates%y = 2
   coordinates%z = 2

   call init_point_probe_output(probe, coordinates, iEx, domain, &
                                test_extension, 3, 0.1_RKIND_tiempo)
   call create_point_probe_output_files(probe)

   n = 10
   do i = 1, n
      probe%timeStep(i) = real(i)
      probe%valueForTime(i) = 10.0*i
      probe%frequencySlice(i) = 0.1*i
      probe%valueForFreq(i) = 0.2*i
   end do

   probe%nTime = n
   probe%nFreq = n

   file_time = trim(adjustl(probe%path))//'_'// &
               trim(adjustl(timeExtension))//'_'// &
               trim(adjustl(datFileExtension))

   file_freq = trim(adjustl(probe%path))//'_'// &
               trim(adjustl(timeExtension))//'_'// &  ! intentional: mirrors implementation
               trim(adjustl(datFileExtension))

   call flush_point_probe_output(probe)

   test_err = test_err + assert_written_output_file(file_time)
   test_err = test_err + assert_written_output_file(file_freq)

   test_err = test_err + assert_integer_equal( &
              probe%nTime, 0, &
              'ERROR: clear_time_data did not reset serializedTimeSize!')

   if (.not. all(probe%timeStep == 0.0) .or. &
       .not. all(probe%valueForTime == 0.0)) then
      print *, 'ERROR: time arrays not cleared!'
      test_err = test_err + 1
   end if

   if (probe%nFreq == 0) then
      print *, 'ERROR: Destroyed frequency reference!'
      test_err = test_err + 1
   end if

   err = test_err
end function test_flush_point_probe

integer function test_multiple_flush_point_probe() bind(c) result(err)
   use output
   use mod_pointProbeOutput
   use mod_domain
   use mod_testOutputUtils
   use mod_assertionTools

   type(point_probe_output_t) :: probe
   type(domain_t)             :: domain
   type(cell_coordinate_t)    :: coordinates

   character(len=BUFSIZE) :: file_time, file_freq
   character(len=36)      :: test_extension

   real(kind=RKIND), allocatable :: expectedTime(:, :)
   real(kind=RKIND), allocatable :: expectedFreq(:, :)

   integer :: n, i
   integer :: test_err = 0

   err = 1
   test_extension = 'tmp_cases/multiple_flush_point_probe'

   domain = domain_t( &
            0.0_RKIND_tiempo, 10.0_RKIND_tiempo, 0.1_RKIND_tiempo, &
            10.0_RKIND, 100.0_RKIND, 10, .false.)

   coordinates%x = 2
   coordinates%y = 2
   coordinates%z = 2

   call init_point_probe_output(probe, coordinates, iEx, domain, &
                                test_extension, 3, 0.1_RKIND_tiempo)
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

   probe%nTime = n
   probe%nFreq = n

   call flush_point_probe_output(probe)

   do i = 1, n
      probe%timeStep(i) = real(i + 10)
      probe%valueForTime(i) = 10.0*(i + 10)
      probe%valueForFreq(i) = -0.5*i

      expectedTime(i + n, 1) = real(i + 10)
      expectedTime(i + n, 2) = 10.0*(i + 10)

      expectedFreq(i, 1) = 0.1*i
      expectedFreq(i, 2) = -0.5*i
   end do

   probe%nTime = n

   call flush_point_probe_output(probe)

   open (unit=probe%fileUnitTime, file=file_time, status='old', action='read')
   test_err = test_err + assert_file_content(probe%fileUnitTime, expectedTime, 2*n, 2)
   close (probe%fileUnitTime)

   open (unit=probe%fileUnitFreq, file=file_freq, status='old', action='read')
   test_err = test_err + assert_file_content(probe%fileUnitFreq, expectedFreq, n, 2)
   close (probe%fileUnitFreq)

   err = test_err
end function test_multiple_flush_point_probe

integer function test_volumic_probe_count_relevant_surfaces() bind(c) result(err)
!   use output
!   use mod_testOutputUtils
!   use FDETYPES_TOOLS
!   use mod_sggMethods
!   use mod_assertionTools
!
!   type(SGGFDTDINFO)              :: dummysgg
!   type(sim_control_t)            :: dummyControl
!   type(bounds_t)                 :: dummyBound
!   type(solver_output_t), pointer :: outputs(:)
!
!   type(media_matrices_t), target :: media
!   type(media_matrices_t), pointer :: mediaPtr
!
!   type(MediaData_t), allocatable, target :: simulationMaterials(:)
!   type(MediaData_t), pointer     :: simulationMaterialsPtr(:)
!   type(MediaData_t)              :: thinWireSimulationMaterial
!
!   type(limit_t), target          :: sinpml_fullsize(6)
!   type(limit_t), pointer         :: sinpml_fullsizePtr(:)
!
!   type(Obses_t)                  :: volumicProbeObservable
!
!   real(kind=RKIND_tiempo), pointer :: timeArray(:)
!   real(kind=RKIND_tiempo)          :: dt = 0.1_RKIND_tiempo
!   integer(kind=SINGLE)             :: nTimeSteps = 100_SINGLE
!
!   integer(kind=RKIND)              :: iter
!   integer(kind=SINGLE)             :: mpidir = 3
!   logical                          :: ThereAreWires = .false.
!   logical                          :: outputRequested
!   integer(kind=SINGLE)             :: test_err = 0
!
!   err = 1
!
!   call sgg_init(dummysgg)
!   call init_time_array(timeArray, nTimeSteps, dt)
!   call sgg_set_tiempo(dummysgg, timeArray)
!   call sgg_set_dt(dummysgg, dt)
!
!   do iter = 1, 6
!      sinpml_fullsize(iter) = create_limit_t(0, 8, 0, 8, 0, 8, 10, 10, 10)
!   end do
!   sinpml_fullsizePtr => sinpml_fullsize
!
!   call init_simulation_material_list(simulationMaterials)
!
!   thinWireSimulationMaterial = create_thinWire_simulation_material(size(simulationMaterials))
!   call add_simulation_material(simulationMaterials, thinWireSimulationMaterial)
!
!   simulationMaterialsPtr => simulationMaterials
!   call sgg_set_NumMedia(dummysgg, size(simulationMaterials))
!   call sgg_set_Med(dummysgg, simulationMaterialsPtr)
!
!   call create_geometry_media(media, 0, 8, 0, 8, 0, 8)
!   call assing_material_id_to_media_matrix_coordinate(media, iEy, 1, 1, 1, simulationMaterials(0)%Id)
!   call assing_material_id_to_media_matrix_coordinate(media, iHz, 1, 1, 1, simulationMaterials(2)%Id)
!   call assing_material_id_to_media_matrix_coordinate(media, iEx, 2, 2, 2, thinWireSimulationMaterial%Id)
!   mediaPtr => media
!
!   volumicProbeObservable = create_volumic_probe_observation(4, 4, 4, 6, 6, 6)
!   call sgg_add_observation(dummysgg, volumicProbeObservable)
!
!   dummyControl = create_control_flags(mpidir=mpidir, nEntradaRoot='entradaRoot', wiresflavor='holland')
!
!   call init_outputs(dummysgg, media, sinpml_fullsize, dummyBound, dummyControl, &
!                     outputRequested, ThereAreWires)
!
!   outputs => GetOutputs()
!
!   test_err = test_err + assert_integer_equal(outputs(1)%outputID, &
!      VOLUMIC_CURRENT_PROBE_ID, 'Unexpected probe id')
!
!   test_err = test_err + assert_integer_equal(outputs(1)%volumicCurrentProbe%columnas, &
!      4, 'Unexpected number of columns')
!
!   test_err = test_err + assert_string_equal(outputs(1)%volumicCurrentProbe%path, &
!      'entradaRoot_volumicProbe_BCX_4_4_4__6_6_6', 'Unexpected path')
!
!   call close_outputs()
!
!   err = test_err
end function

integer function test_init_movie_probe() bind(c) result(err)
   use output
   use outputTypes
   use mod_testOutputUtils
   use FDETYPES_TOOLS
   use mod_sggMethods
   use mod_assertionTools

   type(SGGFDTDINFO)              :: dummysgg
   type(sim_control_t)            :: dummyControl
   type(bounds_t)                 :: dummyBound
   type(solver_output_t), pointer :: outputs(:)

   type(media_matrices_t), target :: media
   type(media_matrices_t), pointer :: mediaPtr

   type(MediaData_t), allocatable, target :: simulationMaterials(:)
   type(MediaData_t), pointer     :: simulationMaterialsPtr(:)

   type(limit_t)                  :: sinpml(6)

   type(Obses_t)                  :: movieObservable
   type(cell_coordinate_t)        :: lowerBoundMovieProbe
   type(cell_coordinate_t)        :: upperBoundMovieProbe

   real(kind=RKIND_tiempo), pointer :: timeArray(:)
   real(kind=RKIND_tiempo)          :: dt = 0.1_RKIND_tiempo
   integer(kind=SINGLE)             :: nTimeSteps = 100_SINGLE

   integer(kind=SINGLE)             :: expectedNumMeasurments
   integer(kind=SINGLE)             :: mpidir = 3
   logical                          :: ThereAreWires = .false.
   logical                          :: outputRequested
   integer(kind=SINGLE)             :: iter
   integer(kind=SINGLE)             :: test_err = 0

   character(len=BUFSIZE) :: test_folder_path = 'tmp_cases/'

   err = 1

   lowerBoundMovieProbe = cell_coordinate_t(2, 2, 2)
   upperBoundMovieProbe = cell_coordinate_t(5, 5, 5)

   call sgg_init(dummysgg)
   call init_time_array(timeArray, nTimeSteps, dt)
   call sgg_set_tiempo(dummysgg, timeArray)
   call sgg_set_dt(dummysgg, dt)

   call init_simulation_material_list(simulationMaterials)
   simulationMaterialsPtr => simulationMaterials
   call sgg_set_NumMedia(dummysgg, size(simulationMaterials))
   call sgg_set_Med(dummysgg, simulationMaterialsPtr)

   call sgg_set_Sweep(dummysgg, create_xyz_limit_array(0, 0, 0, 6, 6, 6))
   call sgg_set_SINPMLSweep(dummysgg, create_xyz_limit_array(1, 1, 1, 5, 5, 5))
   call sgg_set_NumPlaneWaves(dummysgg, 1)
   call sgg_set_Alloc(dummysgg, create_xyz_limit_array(0, 0, 0, 6, 6, 6))

   movieObservable = create_movie_observation(2, 2, 2, 5, 5, 5, iCur)
   call sgg_add_observation(dummysgg, movieObservable)

   call create_geometry_media(media, 0, 8, 0, 8, 0, 8)

   call assing_material_id_to_media_matrix_coordinate(media, iEy, 3, 3, 3, simulationMaterials(0)%Id)
   call assing_material_id_to_media_matrix_coordinate(media, iEy, 4, 3, 3, simulationMaterials(0)%Id)
   call assing_material_id_to_media_matrix_coordinate(media, iEy, 4, 4, 3, simulationMaterials(0)%Id)
   call assing_material_id_to_media_matrix_coordinate(media, iEy, 3, 4, 3, simulationMaterials(0)%Id)

   expectedNumMeasurments = 4_SINGLE
   mediaPtr => media

   do iter = 1, 6
      sinpml(iter) = create_limit_t(0, 8, 0, 8, 0, 8, 10, 10, 10)
   end do

   dummyControl = create_control_flags(nEntradaRoot=test_folder_path, mpidir=mpidir)

   call init_outputs(dummysgg, media, sinpml, dummyBound, dummyControl, &
                     outputRequested, ThereAreWires)

   outputs => GetOutputs()

   test_err = test_err + assert_integer_equal(outputs(1)%outputID, &
                                              MOVIE_PROBE_ID, 'Unexpected probe id')

   test_err = test_err + assert_integer_equal(outputs(1)%movieProbe%nPoints, &
                                              expectedNumMeasurments, 'Unexpected number of measurements')

   test_err = test_err + assert_integer_equal( &
              size(outputs(1)%movieProbe%xValueForTime), &
              expectedNumMeasurments*BuffObse, 'Unexpected allocation size')

   test_err = test_err + assert_integer_equal( &
              size(outputs(1)%movieProbe%timeStep), BuffObse, 'Unexpected timestep buffer size')

   call close_outputs()

   err = test_err
end function

integer function test_update_movie_probe() bind(c) result(err)
   use output
   use outputTypes
   use mod_testOutputUtils
   use FDETYPES_TOOLS
   use mod_sggMethods
   use mod_assertionTools

   type(SGGFDTDINFO)              :: dummysgg
   type(sim_control_t)            :: dummyControl
   type(bounds_t)                 :: dummyBound
   type(solver_output_t), pointer :: outputs(:)

   type(media_matrices_t), target :: media
   type(media_matrices_t), pointer :: mediaPtr

   type(MediaData_t), allocatable, target :: simulationMaterials(:)
   type(MediaData_t), pointer     :: simulationMaterialsPtr(:)

   type(limit_t), target          :: sinpml_fullsize(6)
   type(limit_t), pointer         :: sinpml_fullsizePtr(:)

   type(Obses_t)                  :: movieObservable

   real(kind=RKIND_tiempo), pointer :: timeArray(:)
   real(kind=RKIND_tiempo)          :: dt = 0.1_RKIND_tiempo
   integer(kind=SINGLE)             :: nTimeSteps = 100_SINGLE

   type(dummyFields_t), target     :: dummyFields
   type(fields_reference_t)        :: fields

   integer(kind=SINGLE)             :: expectedNumMeasurments
   integer(kind=SINGLE)             :: mpidir = 3
   integer(kind=SINGLE)             :: iter
   integer(kind=SINGLE)             :: test_err = 0
   logical                          :: ThereAreWires = .false.
   logical                          :: outputRequested

   character(len=BUFSIZE) :: test_folder_path = 'tmp_cases/'

   err = 1

   call sgg_init(dummysgg)
   call init_time_array(timeArray, nTimeSteps, dt)
   call sgg_set_tiempo(dummysgg, timeArray)
   call sgg_set_dt(dummysgg, dt)

   call init_simulation_material_list(simulationMaterials)
   simulationMaterialsPtr => simulationMaterials
   call sgg_set_NumMedia(dummysgg, size(simulationMaterials))
   call sgg_set_Med(dummysgg, simulationMaterialsPtr)

   call sgg_set_Sweep(dummysgg, create_xyz_limit_array(0, 0, 0, 6, 6, 6))
   call sgg_set_SINPMLSweep(dummysgg, create_xyz_limit_array(1, 1, 1, 5, 5, 5))
   call sgg_set_NumPlaneWaves(dummysgg, 1)
   call sgg_set_Alloc(dummysgg, create_xyz_limit_array(0, 0, 0, 6, 6, 6))

   movieObservable = create_movie_observation(2, 2, 2, 5, 5, 5, iCur)
   call sgg_add_observation(dummysgg, movieObservable)

   call create_geometry_media(media, 0, 8, 0, 8, 0, 8)

   call assing_material_id_to_media_matrix_coordinate(media, iEy, 3, 3, 3, simulationMaterials(0)%Id)
   call assing_material_id_to_media_matrix_coordinate(media, iEy, 4, 3, 3, simulationMaterials(0)%Id)
   call assing_material_id_to_media_matrix_coordinate(media, iEy, 4, 4, 3, simulationMaterials(0)%Id)
   call assing_material_id_to_media_matrix_coordinate(media, iEy, 3, 4, 3, simulationMaterials(0)%Id)

   expectedNumMeasurments = 4_SINGLE
   mediaPtr => media

   do iter = 1, 6
      sinpml_fullsize(iter) = create_limit_t(0, 8, 0, 8, 0, 8, 10, 10, 10)
   end do
   sinpml_fullsizePtr => sinpml_fullsize

   dummyControl = create_control_flags(nEntradaRoot=test_folder_path, mpidir=mpidir)

   call init_outputs(dummysgg, media, sinpml_fullsize, dummyBound, dummyControl, &
                     outputRequested, ThereAreWires)

   outputs => GetOutputs()

   call create_dummy_fields(dummyFields, 1, 5, 0.1_RKIND)

   fields%E%x => dummyFields%Ex
   fields%E%y => dummyFields%Ey
   fields%E%z => dummyFields%Ez
   fields%E%deltax => dummyFields%dxe
   fields%E%deltaY => dummyFields%dye
   fields%E%deltaZ => dummyFields%dze
   fields%H%x => dummyFields%Hx
   fields%H%y => dummyFields%Hy
   fields%H%z => dummyFields%Hz
   fields%H%deltax => dummyFields%dxh
   fields%H%deltaY => dummyFields%dyh
   fields%H%deltaZ => dummyFields%dzh

   dummyFields%Hx(3, 3, 3) = 2.0_RKIND
   dummyFields%Hy(3, 3, 3) = 5.0_RKIND
   dummyFields%Hz(3, 3, 3) = 4.0_RKIND

   call update_outputs(dummyControl, dummysgg%tiempo, 1_SINGLE, fields)

   test_err = test_err + assert_real_equal(outputs(1)%movieProbe%yValueForTime(1, 1), &
                                           0.2_RKIND, 1e-5_RKIND, 'Value error')

   test_err = test_err + assert_real_equal(outputs(1)%movieProbe%yValueForTime(1, 2), &
                                           0.0_RKIND, 1e-5_RKIND, 'Value error')

   test_err = test_err + assert_real_equal(outputs(1)%movieProbe%yValueForTime(1, 3), &
                                           0.2_RKIND, 1e-5_RKIND, 'Value error')

   test_err = test_err + assert_real_equal(outputs(1)%movieProbe%yValueForTime(1, 4), &
                                           0.0_RKIND, 1e-5_RKIND, 'Value error')

   test_err = test_err + assert_integer_equal( &
              size(outputs(1)%movieProbe%timeStep), BuffObse, 'Unexpected timestep buffer size')

   call close_outputs()

   err = test_err
end function

integer function test_flush_movie_probe() bind(c) result(err)
   use output
   use outputTypes
   use mod_testOutputUtils
   use FDETYPES_TOOLS
   use mod_sggMethods
   use mod_assertionTools

   type(SGGFDTDINFO)              :: dummysgg
   type(sim_control_t)            :: dummyControl
   type(bounds_t)                 :: dummyBound
   type(solver_output_t), pointer :: outputs(:)

   type(media_matrices_t), target :: media
   type(media_matrices_t), pointer :: mediaPtr

   type(MediaData_t), allocatable, target :: simulationMaterials(:)
   type(MediaData_t), pointer     :: simulationMaterialsPtr(:)

   type(limit_t), target          :: sinpml_fullsize(6)
   type(limit_t), pointer         :: sinpml_fullsizePtr(:)

   type(Obses_t)                  :: movieCurrentObservable
   type(Obses_t)                  :: movieElectricXObservable
   type(Obses_t)                  :: movieMagneticYObservable
   type(fields_reference_t)       :: fields

   real(kind=RKIND_tiempo), pointer :: timeArray(:)
   real(kind=RKIND_tiempo)          :: dt = 0.1_RKIND_tiempo
   integer(kind=SINGLE)             :: nTimeSteps = 100_SINGLE

   integer(kind=SINGLE)             :: expectedNumMeasurments
   integer(kind=SINGLE)             :: mpidir = 3
   integer(kind=SINGLE)             :: iter
   integer(kind=SINGLE)             :: test_err = 0
   logical                          :: ThereAreWires = .false.
   logical                          :: outputRequested

   character(len=BUFSIZE) :: test_folder_path = 'tmp_cases/'
   character(len=BUFSIZE) :: expectedPath

   err = 1

   call sgg_init(dummysgg)
   call init_time_array(timeArray, nTimeSteps, dt)
   call sgg_set_tiempo(dummysgg, timeArray)
   call sgg_set_dt(dummysgg, dt)

   call init_simulation_material_list(simulationMaterials)
   simulationMaterialsPtr => simulationMaterials
   call sgg_set_NumMedia(dummysgg, size(simulationMaterials))
   call sgg_set_Med(dummysgg, simulationMaterialsPtr)

   call sgg_set_Sweep(dummysgg, create_xyz_limit_array(0, 0, 0, 6, 6, 6))
   call sgg_set_SINPMLSweep(dummysgg, create_xyz_limit_array(1, 1, 1, 5, 5, 5))
   call sgg_set_NumPlaneWaves(dummysgg, 1)
   call sgg_set_Alloc(dummysgg, create_xyz_limit_array(0, 0, 0, 6, 6, 6))

   movieCurrentObservable = create_movie_observation(2, 2, 2, 5, 5, 5, iCur)
   call sgg_add_observation(dummysgg, movieCurrentObservable)

   movieElectricXObservable = create_movie_observation(2, 2, 2, 5, 5, 5, iExC)
   call sgg_add_observation(dummysgg, movieElectricXObservable)

   movieMagneticYObservable = create_movie_observation(2, 2, 2, 5, 5, 5, iHyC)
   call sgg_add_observation(dummysgg, movieMagneticYObservable)

   call create_geometry_media(media, 0, 8, 0, 8, 0, 8)

   call assing_material_id_to_media_matrix_coordinate(media, iEx, 3, 3, 3, simulationMaterials(0)%Id)
   call assing_material_id_to_media_matrix_coordinate(media, iEy, 3, 3, 3, simulationMaterials(0)%Id)
   call assing_material_id_to_media_matrix_coordinate(media, iHy, 3, 3, 3, simulationMaterials(0)%Id)

   call assing_material_id_to_media_matrix_coordinate(media, iEx, 3, 4, 3, simulationMaterials(0)%Id)
   call assing_material_id_to_media_matrix_coordinate(media, iEy, 3, 4, 3, simulationMaterials(0)%Id)
   call assing_material_id_to_media_matrix_coordinate(media, iHy, 3, 4, 3, simulationMaterials(0)%Id)

   call assing_material_id_to_media_matrix_coordinate(media, iEx, 4, 4, 3, simulationMaterials(0)%Id)
   call assing_material_id_to_media_matrix_coordinate(media, iEy, 4, 4, 3, simulationMaterials(0)%Id)
   call assing_material_id_to_media_matrix_coordinate(media, iHy, 4, 4, 3, simulationMaterials(0)%Id)

   call assing_material_id_to_media_matrix_coordinate(media, iEx, 4, 3, 3, simulationMaterials(0)%Id)
   call assing_material_id_to_media_matrix_coordinate(media, iEy, 4, 3, 3, simulationMaterials(0)%Id)
   call assing_material_id_to_media_matrix_coordinate(media, iHy, 4, 3, 3, simulationMaterials(0)%Id)

   expectedNumMeasurments = 4_SINGLE
   mediaPtr => media

   do iter = 1, 6
      sinpml_fullsize(iter) = create_limit_t(0, 8, 0, 8, 0, 8, 10, 10, 10)
   end do
   sinpml_fullsizePtr => sinpml_fullsize

   dummyControl = create_control_flags(nEntradaRoot=test_folder_path, mpidir=mpidir)

   call init_outputs(dummysgg, media, sinpml_fullsize, dummyBound, dummyControl, &
                     outputRequested, ThereAreWires)

   outputs => GetOutputs()

   !--- Dummy first update ---
   !movieCurrentObservable
   outputs(1)%movieProbe%nTime = 1
   outputs(1)%movieProbe%timeStep(1) = 0.5_RKIND_tiempo
   outputs(1)%movieProbe%xValueForTime(1, :) = [0.1_RKIND, 0.2_RKIND, 0.3_RKIND, 0.4_RKIND]
   outputs(1)%movieProbe%yValueForTime(1, :) = [0.3_RKIND, 0.4_RKIND, 0.5_RKIND, 0.6_RKIND]
   outputs(1)%movieProbe%zValueForTime(1, :) = [0.7_RKIND, 0.8_RKIND, 0.9_RKIND, 1.0_RKIND]

   !movieElectricXObservable
   outputs(2)%movieProbe%nTime = 1
   outputs(2)%movieProbe%timeStep(1) = 0.5_RKIND_tiempo
   outputs(2)%movieProbe%xValueForTime(1, :) = [0.1_RKIND, 0.2_RKIND, 0.3_RKIND, 0.4_RKIND]

   !movieMagneticYObservable
   outputs(3)%movieProbe%nTime = 1
   outputs(3)%movieProbe%timeStep(1) = 0.5_RKIND_tiempo
   outputs(3)%movieProbe%yValueForTime(1, :) = [0.1_RKIND, 0.2_RKIND, 0.3_RKIND, 0.4_RKIND]

   !--- Dummy second update ---
   !movieCurrentObservable
   outputs(1)%movieProbe%nTime = 2
   outputs(1)%movieProbe%timeStep(2) = 0.75_RKIND_tiempo
   outputs(1)%movieProbe%xValueForTime(2, :) = [1.1_RKIND, 1.2_RKIND, 1.3_RKIND, 1.4_RKIND]
   outputs(1)%movieProbe%yValueForTime(2, :) = [1.3_RKIND, 1.4_RKIND, 1.5_RKIND, 1.6_RKIND]
   outputs(1)%movieProbe%zValueForTime(2, :) = [1.7_RKIND, 1.8_RKIND, 1.9_RKIND, 2.0_RKIND]

   !movieElectricXObservable
   outputs(2)%movieProbe%nTime = 2
   outputs(2)%movieProbe%timeStep(2) = 0.75_RKIND_tiempo
   outputs(2)%movieProbe%xValueForTime(2, :) = [1.1_RKIND, 1.2_RKIND, 1.3_RKIND, 1.4_RKIND]

   !movieMagneticYObservable
   outputs(3)%movieProbe%nTime = 2
   outputs(3)%movieProbe%timeStep(2) = 0.75_RKIND_tiempo
   outputs(3)%movieProbe%yValueForTime(2, :) = [1.1_RKIND, 1.2_RKIND, 1.3_RKIND, 1.4_RKIND]

   call flush_outputs(dummysgg%tiempo, 1_SINGLE, dummyControl, fields, dummyBound, .false.)

   ! --- Assert file existance
   do outputIdx = 1, 3
      expectedPath = trim(adjustl(outputs(outputIdx)%movieProbe%path))//'_ts0001.vtu'
      test_err = test_err + assert_file_exists(expectedPath)

      expectedPath = trim(adjustl(outputs(outputIdx)%movieProbe%path))//'_ts0002.vtu'
      test_err = test_err + assert_file_exists(expectedPath)
   end do

   call close_outputs()

   expectedPath = trim(adjustl(outputs(1)%movieProbe%path))//'.pvd'
   test_err = test_err + assert_file_exists(expectedPath)

   err = test_err
end function

integer function test_init_frequency_slice_probe() bind(c) result(err)
   use output
   use outputTypes
   use mod_testOutputUtils
   use FDETYPES_TOOLS
   use mod_sggMethods
   use mod_assertionTools

   type(SGGFDTDINFO)              :: dummysgg
   type(sim_control_t)            :: dummyControl
   type(bounds_t)                 :: dummyBound
   type(solver_output_t), pointer :: outputs(:)

   type(media_matrices_t), target :: media
   type(media_matrices_t), pointer :: mediaPtr

   type(MediaData_t), allocatable, target :: simulationMaterials(:)
   type(MediaData_t), pointer     :: simulationMaterialsPtr(:)

   type(limit_t), target          :: sinpml_fullsize(6)
   type(limit_t), pointer         :: sinpml_fullsizePtr(:)

   type(Obses_t)                  :: frequencySliceObservation

   real(kind=RKIND_tiempo), pointer :: timeArray(:)
   real(kind=RKIND_tiempo)          :: dt = 0.1_RKIND_tiempo
   integer(kind=SINGLE)             :: nTimeSteps = 100_SINGLE

   integer(kind=SINGLE)             :: expectedNumMeasurments
   integer(kind=SINGLE)             :: expectedTotalFrequnecies
   integer(kind=SINGLE)             :: mpidir = 3
   integer(kind=SINGLE)             :: iter
   integer(kind=SINGLE)             :: test_err = 0
   logical                          :: ThereAreWires = .false.
   logical                          :: outputRequested

   character(len=BUFSIZE) :: test_folder_path = 'tmp_cases/'

   err = 1

   call sgg_init(dummysgg)

   call init_time_array(timeArray, nTimeSteps, dt)
   call sgg_set_tiempo(dummysgg, timeArray)
   call sgg_set_dt(dummysgg, dt)

   call init_simulation_material_list(simulationMaterials)
   simulationMaterialsPtr => simulationMaterials
   call sgg_set_NumMedia(dummysgg, size(simulationMaterials))
   call sgg_set_Med(dummysgg, simulationMaterialsPtr)

   call sgg_set_Sweep(dummysgg, create_xyz_limit_array(0, 0, 0, 6, 6, 6))
   call sgg_set_SINPMLSweep(dummysgg, create_xyz_limit_array(1, 1, 1, 5, 5, 5))
   call sgg_set_NumPlaneWaves(dummysgg, 1)
   call sgg_set_Alloc(dummysgg, create_xyz_limit_array(0, 0, 0, 6, 6, 6))

   frequencySliceObservation = create_frequency_slice_observation(2, 2, 2, 5, 5, 5, iCur)
   call sgg_add_observation(dummysgg, frequencySliceObservation)

   expectedTotalFrequnecies = 6_SINGLE

   call create_geometry_media(media, 0, 8, 0, 8, 0, 8)

   call assing_material_id_to_media_matrix_coordinate(media, iEy, 3, 3, 3, simulationMaterials(0)%Id)
   call assing_material_id_to_media_matrix_coordinate(media, iEy, 4, 3, 3, simulationMaterials(0)%Id)
   call assing_material_id_to_media_matrix_coordinate(media, iEy, 4, 4, 3, simulationMaterials(0)%Id)
   call assing_material_id_to_media_matrix_coordinate(media, iEy, 3, 4, 3, simulationMaterials(0)%Id)

   expectedNumMeasurments = 4_SINGLE
   mediaPtr => media

   do iter = 1, 6
      sinpml_fullsize(iter) = create_limit_t(0, 8, 0, 8, 0, 8, 10, 10, 10)
   end do
   sinpml_fullsizePtr => sinpml_fullsize

   dummyControl = create_control_flags(nEntradaRoot=test_folder_path, mpidir=mpidir)

   call init_outputs(dummysgg, media, sinpml_fullsize, dummyBound, dummyControl, &
                     outputRequested, ThereAreWires)

   outputs => GetOutputs()

   test_err = test_err + assert_integer_equal(outputs(1)%outputID, &
                                              FREQUENCY_SLICE_PROBE_ID, 'Unexpected probe id')

   test_err = test_err + assert_integer_equal(outputs(1)%frequencySliceProbe%nFreq, &
                                              6, 'Unexpected number of frequencies')

   test_err = test_err + assert_integer_equal(outputs(1)%frequencySliceProbe%nPoints, &
                                              expectedNumMeasurments, 'Unexpected number of measurements')

   test_err = test_err + assert_integer_equal( &
              size(outputs(1)%frequencySliceProbe%xValueForFreq), &
              expectedNumMeasurments*expectedTotalFrequnecies, 'Unexpected allocation size')

   test_err = test_err + assert_integer_equal( &
              size(outputs(1)%frequencySliceProbe%frequencySlice), &
              expectedTotalFrequnecies, 'Unexpected frequency count')

   call close_outputs()

   err = test_err
end function

integer function test_update_frequency_slice_probe() bind(c) result(err)
   use output
   use outputTypes
   use mod_testOutputUtils
   use FDETYPES_TOOLS
   use mod_sggMethods
   use mod_assertionTools

   type(SGGFDTDINFO)              :: dummysgg
   type(sim_control_t)            :: dummyControl
   type(bounds_t)                 :: dummyBound
   type(solver_output_t), pointer :: outputs(:)

   type(media_matrices_t), target :: media
   type(media_matrices_t), pointer :: mediaPtr

   type(MediaData_t), allocatable, target :: simulationMaterials(:)
   type(MediaData_t), pointer     :: simulationMaterialsPtr(:)

   type(limit_t), target          :: sinpml_fullsize(6)
   type(limit_t), pointer         :: sinpml_fullsizePtr(:)

   type(Obses_t)                  :: frequencySliceObservation

   real(kind=RKIND_tiempo), pointer :: timeArray(:)
   real(kind=RKIND_tiempo)          :: dt = 0.1_RKIND_tiempo
   integer(kind=SINGLE)             :: nTimeSteps = 100_SINGLE

   type(dummyFields_t), target     :: dummyFields
   type(fields_reference_t)        :: fields

   integer(kind=SINGLE)             :: expectedNumMeasurments
   integer(kind=SINGLE)             :: expectedNumberFrequencies
   integer(kind=SINGLE)             :: mpidir = 3
   integer(kind=SINGLE)             :: iter
   integer(kind=SINGLE)             :: test_err = 0
   logical                          :: ThereAreWires = .false.
   logical                          :: outputRequested

   character(len=BUFSIZE) :: test_folder_path = 'tmp_cases/'

   err = 1

   call sgg_init(dummysgg)
   call init_time_array(timeArray, nTimeSteps, dt)
   call sgg_set_tiempo(dummysgg, timeArray)
   call sgg_set_dt(dummysgg, dt)

   call init_simulation_material_list(simulationMaterials)
   simulationMaterialsPtr => simulationMaterials
   call sgg_set_NumMedia(dummysgg, size(simulationMaterials))
   call sgg_set_Med(dummysgg, simulationMaterialsPtr)

   call sgg_set_Sweep(dummysgg, create_xyz_limit_array(0, 0, 0, 6, 6, 6))
   call sgg_set_SINPMLSweep(dummysgg, create_xyz_limit_array(1, 1, 1, 5, 5, 5))
   call sgg_set_NumPlaneWaves(dummysgg, 1)
   call sgg_set_Alloc(dummysgg, create_xyz_limit_array(0, 0, 0, 6, 6, 6))

   frequencySliceObservation = create_frequency_slice_observation(2, 2, 2, 5, 5, 5, iCur)
   call sgg_add_observation(dummysgg, frequencySliceObservation)

   call create_geometry_media(media, 0, 8, 0, 8, 0, 8)

   call assing_material_id_to_media_matrix_coordinate(media, iEy, 3, 3, 3, simulationMaterials(0)%Id)
   call assing_material_id_to_media_matrix_coordinate(media, iEy, 4, 3, 3, simulationMaterials(0)%Id)
   call assing_material_id_to_media_matrix_coordinate(media, iEy, 4, 4, 3, simulationMaterials(0)%Id)
   call assing_material_id_to_media_matrix_coordinate(media, iEy, 3, 4, 3, simulationMaterials(0)%Id)
   expectedNumberFrequencies = 6_SINGLE
   expectedNumMeasurments = 4_SINGLE
   mediaPtr => media

   do iter = 1, 6
      sinpml_fullsize(iter) = create_limit_t(0, 8, 0, 8, 0, 8, 10, 10, 10)
   end do
   sinpml_fullsizePtr => sinpml_fullsize

   dummyControl = create_control_flags(nEntradaRoot=test_folder_path, mpidir=mpidir)

   call init_outputs(dummysgg, media, sinpml_fullsize, dummyBound, dummyControl, &
                     outputRequested, ThereAreWires)

   outputs => GetOutputs()

   call create_dummy_fields(dummyFields, 1, 5, 0.1_RKIND)

   fields%E%x => dummyFields%Ex
   fields%E%y => dummyFields%Ey
   fields%E%z => dummyFields%Ez
   fields%E%deltax => dummyFields%dxe
   fields%E%deltaY => dummyFields%dye
   fields%E%deltaZ => dummyFields%dze
   fields%H%x => dummyFields%Hx
   fields%H%y => dummyFields%Hy
   fields%H%z => dummyFields%Hz
   fields%H%deltax => dummyFields%dxh
   fields%H%deltaY => dummyFields%dyh
   fields%H%deltaZ => dummyFields%dzh

   call fillGradient(dummyFields, 1, 0.0_RKIND, 10.0_RKIND)

   call update_outputs(dummyControl, dummysgg%tiempo, 2_SINGLE, fields)

   test_err = test_err + assert_integer_equal(outputs(1)%outputID, &
                                              FREQUENCY_SLICE_PROBE_ID, 'Unexpected probe id')

   test_err = test_err + assert_integer_equal(outputs(1)%frequencySliceProbe%nPoints, &
                                              expectedNumMeasurments, 'Unexpected number of measurements')

   test_err = test_err + assert_integer_equal( &
              size(outputs(1)%frequencySliceProbe%frequencySlice), &
              expectedNumberFrequencies, 'Unexpected allocation size')

   !This test generates X Gradient for H. It is expected to detect none Current accros X axis and Opposite values for Y and Z

   test_err = test_err + assert_array_value(outputs(1)%frequencySliceProbe%xValueForFreq, (0.0_CKIND , 0.0_CKIND), errormessage='Detected Current on X Axis for Hx gradient')
   test_err = test_err + assert_arrays_equal(outputs(1)%frequencySliceProbe%yValueForFreq, &
                                             -1.0_RKIND * outputs(1)%frequencySliceProbe%zValueForFreq, errormessage='Unequal values for Y and -Z')
   
   call close_outputs()

   err = test_err
end function

integer function test_flush_frequency_slice_probe() bind(c) result(err)
   use output
   use outputTypes
   use mod_testOutputUtils
   use FDETYPES_TOOLS
   use mod_sggMethods
   use mod_assertionTools

   type(SGGFDTDINFO)              :: dummysgg
   type(sim_control_t)            :: dummyControl
   type(bounds_t)                 :: dummyBound
   type(solver_output_t), pointer :: outputs(:)

   type(media_matrices_t), target :: media
   type(media_matrices_t), pointer :: mediaPtr

   type(MediaData_t), allocatable, target :: simulationMaterials(:)
   type(MediaData_t), pointer     :: simulationMaterialsPtr(:)

   type(limit_t), target          :: sinpml_fullsize(6)
   type(limit_t), pointer         :: sinpml_fullsizePtr(:)

   type(Obses_t)                  :: frequencySliceCurrentObservable
   type(Obses_t)                  :: frequencySliceElectricXObservable
   type(Obses_t)                  :: frequencySliceMagneticHObservable
   type(fields_reference_t)       :: fields
   type(dummyFields_t), target    :: dummyFields

   real(kind=RKIND_tiempo), pointer :: timeArray(:)
   real(kind=RKIND_tiempo)          :: dt = 0.1_RKIND_tiempo
   integer(kind=SINGLE)             :: nTimeSteps = 100_SINGLE
   integer(kind=SINGLE)             :: expectedNumMeasurments
   integer(kind=SINGLE)             :: expectedNumFrequencies
   integer(kind=SINGLE)             :: mpidir = 3
   integer(kind=SINGLE)             :: iter
   integer(kind=SINGLE)             :: test_err = 0
   logical                          :: ThereAreWires = .false.
   logical                          :: outputRequested
   character(len=BUFSIZE)           :: test_folder_path = 'tmp_cases/'
   character(len=BUFSIZE)           :: expectedPath
   character(len=3) :: freqIdName

   err = 1

   !--- Setup SGG ---
   call sgg_init(dummysgg)
   call init_time_array(timeArray, nTimeSteps, dt)
   call sgg_set_tiempo(dummysgg, timeArray)
   call sgg_set_dt(dummysgg, dt)

   call init_simulation_material_list(simulationMaterials)
   simulationMaterialsPtr => simulationMaterials
   call sgg_set_NumMedia(dummysgg, size(simulationMaterials))
   call sgg_set_Med(dummysgg, simulationMaterialsPtr)

   call sgg_set_Sweep(dummysgg, create_xyz_limit_array(0, 0, 0, 6, 6, 6))
   call sgg_set_SINPMLSweep(dummysgg, create_xyz_limit_array(1, 1, 1, 5, 5, 5))
   call sgg_set_NumPlaneWaves(dummysgg, 1)
   call sgg_set_Alloc(dummysgg, create_xyz_limit_array(0, 0, 0, 6, 6, 6))

   frequencySliceCurrentObservable = create_frequency_slice_observation(2, 2, 2, 5, 5, 5, iCur)
   call sgg_add_observation(dummysgg, frequencySliceCurrentObservable)

   frequencySliceElectricXObservable = create_frequency_slice_observation(2, 2, 2, 5, 5, 5, iExC)
   call sgg_add_observation(dummysgg, frequencySliceElectricXObservable)

   frequencySliceMagneticHObservable = create_frequency_slice_observation(2, 2, 2, 5, 5, 5, iHyC)
   call sgg_add_observation(dummysgg, frequencySliceMagneticHObservable)

   call create_geometry_media(media, 0, 8, 0, 8, 0, 8)
   call assing_material_id_to_media_matrix_coordinate(media, iEy, 3, 3, 3, simulationMaterials(0)%Id)
   call assing_material_id_to_media_matrix_coordinate(media, iEy, 4, 3, 3, simulationMaterials(0)%Id)
   call assing_material_id_to_media_matrix_coordinate(media, iEy, 4, 4, 3, simulationMaterials(0)%Id)
   call assing_material_id_to_media_matrix_coordinate(media, iEy, 3, 4, 3, simulationMaterials(0)%Id)

   expectedNumFrequencies = 6_SINGLE
   expectedNumMeasurments = 4_SINGLE

   mediaPtr => media

   do iter = 1, 6
      sinpml_fullsize(iter) = create_limit_t(0, 8, 0, 8, 0, 8, 10, 10, 10)
   end do
   sinpml_fullsizePtr => sinpml_fullsize

   dummyControl = create_control_flags(nEntradaRoot=test_folder_path, mpidir=mpidir)

   call init_outputs(dummysgg, media, sinpml_fullsize, dummyBound, dummyControl, &
                     outputRequested, ThereAreWires)
   outputs => GetOutputs()

   !--- Dummy update ---
   !frequencySliceObservable
   do freq = 1, expectedNumFrequencies
      outputs(1)%frequencySliceProbe%xvalueForFreq(freq, :) = [(0.1_RKIND, 0.1_RKIND), (0.2_RKIND, 0.2_RKIND), (0.3_RKIND, 0.3_RKIND), (0.4_RKIND, 0.4_RKIND)]
      outputs(1)%frequencySliceProbe%yvalueForFreq(freq, :) = [(0.5_RKIND, 0.5_RKIND), (0.6_RKIND, 0.6_RKIND), (0.7_RKIND, 0.7_RKIND), (0.8_RKIND, 0.8_RKIND)]
      outputs(1)%frequencySliceProbe%zvalueForFreq(freq, :) = [(0.9_RKIND, 0.9_RKIND), (1.0_RKIND, 1.0_RKIND), (1.1_RKIND, 1.1_RKIND), (1.2_RKIND, 1.2_RKIND)]
   end do
   !frequencySliceXObservable 
   do freq = 1, expectedNumFrequencies
      outputs(2)%frequencySliceProbe%xvalueForFreq(freq, :) = [(0.1_RKIND, 0.1_RKIND), (0.2_RKIND, 0.2_RKIND), (0.3_RKIND, 0.3_RKIND), (0.4_RKIND, 0.4_RKIND)]
   end do
   !frequencySliceYObservable 
   do freq = 1, expectedNumFrequencies
      outputs(3)%frequencySliceProbe%yvalueForFreq(freq, :) = [(0.1_RKIND, 0.1_RKIND), (0.2_RKIND, 0.2_RKIND), (0.3_RKIND, 0.3_RKIND), (0.4_RKIND, 0.4_RKIND)]
   end do

   call flush_outputs(dummysgg%tiempo, 1_SINGLE, dummyControl, fields, dummyBound, .false.)

   !--- Assert generated files ---
   do iter = 1, expectedNumFrequencies
      write(freqIdName, '(i3)') iter
      expectedPath = trim(adjustl(outputs(1)%frequencySliceProbe%path))//'_fq'//'000'//trim(adjustl(freqIdName))//'.vtu'
      test_err = test_err + assert_file_exists(expectedPath)
   end do

   call close_outputs()

   expectedPath = trim(adjustl(outputs(1)%frequencySliceProbe%path))//'.pvd'
   test_err = test_err + assert_file_exists(expectedPath)

   err = test_err
end function

