integer function test_init_point_probe() bind(c) result(err)
! This test initializes a single point probe at coordinates (4,4,4).
! It verifies that the probe is correctly registered in the simulation outputs,
! that the output ID matches POINT_PROBE_ID, and that the output paths for
! both the probe and its time data file are correctly set and exist.
   use FDETYPES_m
   use FDETYPES_TOOLS
   use output_m
   use outputTypes_m
   use testOutputUtils_m
   use sggMethods_m
   use assertionTools_m
   use directoryUtils_m
   implicit none

   ! Parameters
   character(len=14), parameter :: test_folder = 'testing_folder'
   character(len=18), parameter :: test_name = 'initPointProbeTest'

   ! Local variables
   character(len=1) :: sep
   character(len=BUFSIZE) :: nEntrada
   character(len=BUFSIZE) :: expectedProbePath
   character(len=BUFSIZE) :: expectedDataPath

   type(SGGFDTDINFO_t)              :: sgg
   type(sim_control_t)            :: control
   type(bounds_t)                 :: bounds
   type(media_matrices_t)         :: media
   type(limit_t), allocatable     :: sinpml(:)
   type(Obses_t)                  :: probe
   type(solver_output_t), pointer :: outputs(:)
   type(MediaData_t), allocatable, target :: materials(:)
   type(MediaData_t), pointer     :: materialsPtr(:)
   type(taglist_t)                :: tagNumbers

   real(kind=RKIND_tiempo), pointer :: timeArray(:)
   real(kind=RKIND_tiempo)          :: dt = 0.1_RKIND_tiempo
   integer(kind=SINGLE)             :: nSteps = 100_SINGLE

   logical :: outputRequested
   logical :: hasWires = .false.
   integer(kind=SINGLE) :: test_err = 0
   integer :: ios

   ! Setup
   sep = get_path_separator()
   nEntrada = test_folder//sep//test_name

   call sgg_init(sgg)
   call init_time_array(timeArray, nSteps, dt)
   call sgg_set_tiempo(sgg, timeArray)
   call sgg_set_dt(sgg, dt)

   call init_simulation_material_list(materials)
   materialsPtr => materials
   call sgg_set_Med(sgg, materialsPtr)

   probe = create_point_probe_observation(4, 4, 4)
   call sgg_add_observation(sgg, probe)

   control = create_control_flags(mpidir=3, nEntradaRoot=trim(nEntrada), wiresflavor='holland')

   ! Action
   call init_outputs(sgg, media, sinpml, tagNumbers, bounds, control, outputRequested, hasWires)
   outputs => GetOutputs()

   ! Assertions
   test_err = test_err + assert_true(outputRequested, 'Valid probes not found')
   test_err = test_err + assert_integer_equal(outputs(1)%outputID, POINT_PROBE_ID, 'Unexpected probe id')

   expectedProbePath = trim(nEntrada)//wordSeparation//'pointProbe_Ex_4_4_4'
   expectedDataPath = trim(expectedProbePath)//wordSeparation//timeExtension//datFileExtension

   test_err = test_err + assert_string_equal(outputs(1)%pointProbe%path, expectedProbePath, 'Unexpected path')
   test_err = test_err + assert_string_equal(outputs(1)%pointProbe%filePathTime, expectedDataPath, 'Unexpected path')
   test_err = test_err + assert_true(file_exists(expectedDataPath), 'Time data file do not exist')

   ! Cleanup
   call remove_folder(test_folder, ios)
   deallocate (sgg%Observation, outputs)

   err = test_err
end function

integer function test_update_point_probe() bind(c) result(err)
! This test updates the values recorded by a single point probe at (4,4,4)
! over two timesteps. It verifies that the probe correctly stores the time
! steps and associated field values, ensuring proper temporal tracking of
! measured quantities.
   use FDETYPES_m
   use FDETYPES_TOOLS
   use output_m
   use outputTypes_m
   use testOutputUtils_m
   use sggMethods_m
   use assertionTools_m
   use directoryUtils_m
   implicit none

   ! Parameters
   character(len=14), parameter :: test_folder = 'testing_folder'
   character(len=20), parameter :: test_name = 'updatePointProbeTest'

   ! Local variables
   character(len=1) :: sep
   character(len=BUFSIZE) :: nEntrada

   type(SGGFDTDINFO_t)              :: sgg
   type(sim_control_t)            :: control
   type(bounds_t)                 :: bounds
   type(media_matrices_t)         :: media
   type(limit_t), allocatable     :: sinpml(:)
   type(Obses_t)                  :: probe
   type(solver_output_t), pointer :: outputs(:)
   type(MediaData_t), allocatable, target :: materials(:)
   type(MediaData_t), pointer     :: materialsPtr(:)
   type(taglist_t)                :: tagNumbers

   type(dummyFields_t), target    :: dummyFields
   type(fields_reference_t)       :: fields

   real(kind=RKIND_tiempo), pointer :: timeArray(:)
   real(kind=RKIND_tiempo)          :: dt = 0.1_RKIND_tiempo
   integer(kind=SINGLE)             :: nSteps = 100_SINGLE

   logical :: outputRequested
   logical :: hasWires = .false.
   integer(kind=SINGLE) :: test_err = 0
   integer :: ios

   ! Setup
   sep = get_path_separator()
   nEntrada = test_folder//sep//test_name

   call sgg_init(sgg)
   call init_time_array(timeArray, nSteps, dt)
   call sgg_set_tiempo(sgg, timeArray)
   call sgg_set_dt(sgg, dt)

   probe = create_point_probe_observation(4, 4, 4)
   call sgg_add_observation(sgg, probe)

   call init_simulation_material_list(materials)
   materialsPtr => materials
   call sgg_set_Med(sgg, materialsPtr)

   control = create_control_flags(mpidir=3, nEntradaRoot=nEntrada, wiresflavor='holland')
   call init_outputs(sgg, media, sinpml, tagNumbers, bounds, control, outputRequested, hasWires)

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

   ! Action
   dummyFields%Ex(4, 4, 4) = 5.0_RKIND
   call update_outputs(control, sgg%tiempo, 1_SINGLE, fields)
   outputs => GetOutputs()

   ! Assertions
   test_err = test_err + assert_real_equal(outputs(1)%pointProbe%timeStep(1), 0.0_RKIND_tiempo, 1e-5_RKIND_tiempo, 'Unexpected timestep 1')
   test_err = test_err + assert_real_equal(outputs(1)%pointProbe%valueForTime(1), 5.0_RKIND, 1e-5_RKIND, 'Unexpected field 1')

   dummyFields%Ex(4, 4, 4) = -4.0_RKIND
   call update_outputs(control, sgg%tiempo, 2_SINGLE, fields)

   test_err = test_err + assert_real_equal(outputs(1)%pointProbe%timeStep(2), 0.1_RKIND_tiempo, 1e-5_RKIND_tiempo, 'Unexpected timestep 2')
   test_err = test_err + assert_real_equal(outputs(1)%pointProbe%valueForTime(2), -4.0_RKIND, 1e-5_RKIND, 'Unexpected field 2')

   !Cleanup
   call remove_folder(test_folder, ios)

   err = test_err
end function

integer function test_flush_point_probe() bind(c) result(err)
! This test validates the flush operation for a point probe. It ensures
! that time and frequency data are properly written to files, and that
! internal arrays are cleared/reset after flushing, preserving data integrity.
   use output_m
   use outputTypes_m
   use pointProbeOutput_m
   use domain_m
   use testOutputUtils_m
   use assertionTools_m
   use directoryUtils_m
   implicit none

   ! Parameters
   character(len=14), parameter :: test_folder = 'testing_folder'
   character(len=19), parameter :: test_name = 'flushPointProbeTest'

   ! Local variables
   character(len=1) :: sep
   character(len=BUFSIZE) :: nEntrada

   type(point_probe_output_t) :: probe
   type(domain_t)             :: domain
   type(cell_coordinate_t)    :: coordinates

   integer :: n, i
   integer :: test_err = 0
   integer :: ios

   ! Setup
   sep = get_path_separator()
   nEntrada = test_folder//sep//test_name

   domain = domain_t( &
            0.0_RKIND_tiempo, 10.0_RKIND_tiempo, 0.1_RKIND_tiempo, &
            10.0_RKIND, 100.0_RKIND, 10, .false.)

   coordinates%x = 2
   coordinates%y = 2
   coordinates%z = 2

   call init_point_probe_output(probe, coordinates, iEx, domain, nEntrada, 3, 0.1_RKIND_tiempo)

   ! Action
   n = 10
   do i = 1, n
      probe%timeStep(i) = real(i)
      probe%valueForTime(i) = 10.0*i
      probe%frequencySlice(i) = 0.1*i
      probe%valueForFreq(i) = 0.2*i
   end do

   probe%nTime = n
   probe%nFreq = n

   call flush_point_probe_output(probe)

   ! Assertions
   test_err = test_err + assert_written_output_file(probe%filePathTime)
   test_err = test_err + assert_written_output_file(probe%filePathFreq)

   test_err = test_err + assert_integer_equal(probe%nTime, 0, 'ERROR: clear_time_data did not reset serializedTimeSize!')

   if (.not. all(probe%timeStep == 0.0) .or. .not. all(probe%valueForTime == 0.0)) then
      print *, 'ERROR: time arrays not cleared!'
      test_err = test_err + 1
   end if

   if (probe%nFreq == 0) then
      print *, 'ERROR: Destroyed frequency reference!'
      test_err = test_err + 1
   end if

   !Cleanup
   call remove_folder(test_folder, ios)

   err = test_err
end function

integer function test_multiple_flush_point_probe() bind(c) result(err)
! This test verifies that multiple consecutive flushes of a point probe
! correctly append or overwrite data files without losing previous data.
! It ensures consistency of both time and frequency outputs across multiple flushes.
   use output_m
   use outputTypes_m
   use pointProbeOutput_m
   use domain_m
   use testOutputUtils_m
   use assertionTools_m
   use directoryUtils_m
   implicit none

   ! Parameters
   character(len=14), parameter :: test_folder = 'testing_folder'
   character(len=27), parameter :: test_name = 'flushMultiplePointProbeTest'

   ! Local variables
   character(len=1) :: sep
   character(len=BUFSIZE) :: nEntrada

   type(point_probe_output_t) :: probe
   type(domain_t)             :: domain
   type(cell_coordinate_t)    :: coordinates

   real(kind=RKIND), allocatable :: expectedTime(:, :)
   real(kind=RKIND), allocatable :: expectedFreq(:, :)

   integer :: n, i, unit
   integer :: test_err = 0
   integer :: ios

   ! Setup
   sep = get_path_separator()
   nEntrada = test_folder//sep//test_name

   domain = domain_t( &
            0.0_RKIND_tiempo, 10.0_RKIND_tiempo, 0.1_RKIND_tiempo, &
            10.0_RKIND, 100.0_RKIND, 10, .false.)

   coordinates%x = 2
   coordinates%y = 2
   coordinates%z = 2

   call init_point_probe_output(probe, coordinates, iEx, domain, nEntrada, 3, 0.1_RKIND_tiempo)

   n = 10
   allocate (expectedTime(2*n, 2))
   allocate (expectedFreq(n, 2))

   ! Action - first flush
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

   ! Action - second flush
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

   ! Assertions
   open (unit=unit, file=probe%filePathTime, status='old', action='read')
   test_err = test_err + assert_file_content(unit, expectedTime, 2*n, 2, 1e-06_RKIND)
   close (unit)

   open (unit=unit, file=probe%filePathFreq, status='old', action='read')
   test_err = test_err + assert_file_content(unit, expectedFreq, n, 2, 1e-06_RKIND)
   close (unit)

   !Cleanup
   call remove_folder(test_folder, ios)

   err = test_err
end function

integer function test_init_movie_probe() bind(c) result(err)
! This test initializes a movie probe over a 3D region from (2,2,2) to (5,5,5).
! It checks that the probe is correctly allocated, that the number of measurement
! points and buffer sizes are correct, and that the output folder and PVD file
! for the movie are properly created.
   use output_m
   use outputTypes_m
   use testOutputUtils_m
   use FDETYPES_TOOLS
   use sggMethods_m
   use assertionTools_m
   use directoryUtils_m
   implicit none

   type(SGGFDTDINFO_t)              :: dummysgg
   type(sim_control_t)            :: dummyControl
   type(bounds_t)                 :: dummyBound
   type(XYZlimit_t)               :: dummySweep(6)
   type(XYZlimit_t)               :: dummySinpmlSweep(6)
   type(XYZlimit_t)               :: allocationRange(6)
   type(solver_output_t), pointer :: outputs(:)

   type(media_matrices_t), target :: media
   type(media_matrices_t), pointer :: mediaPtr

   type(MediaData_t), allocatable, target :: simulationMaterials(:)
   type(MediaData_t), pointer     :: simulationMaterialsPtr(:)

   type(taglist_t)                :: tagNumbers

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

   character(len=14), parameter :: test_folder = 'testing_folder'
   character(len=9), parameter :: test_name = 'initMovie'

   character(len=BUFSIZE) :: nEntrada
   character(len=1) :: sep
   character(len=BUFSIZE) :: expectedProbePath
   character(len=BUFSIZE) :: pdvFileName
   integer :: ios

   sep = get_path_separator()
   nEntrada = test_folder//sep//test_name

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
   dummySweep = create_xyz_limit_array(0, 0, 0, 6, 6, 6)
   call sgg_set_Sweep(dummysgg, dummySweep)
   dummySinpmlSweep = create_xyz_limit_array(1, 1, 1, 5, 5, 5)
   call sgg_set_SINPMLSweep(dummysgg, dummySinpmlSweep)
   call sgg_set_NumPlaneWaves(dummysgg, 1)
   allocationRange = create_xyz_limit_array(0, 0, 0, 6, 6, 6)
   call sgg_set_Alloc(dummysgg, allocationRange)

   movieObservable = create_movie_observation(2, 2, 2, 5, 5, 5, iCur)
   call sgg_add_observation(dummysgg, movieObservable)

   call create_geometry_media(media, 0, 8, 0, 8, 0, 8)

   call assign_material_id_to_media_matrix_coordinate(media, iEy, 3, 3, 3, simulationMaterials(0)%Id)
   call assign_material_id_to_media_matrix_coordinate(media, iEy, 4, 3, 3, simulationMaterials(0)%Id)
   call assign_material_id_to_media_matrix_coordinate(media, iEy, 4, 4, 3, simulationMaterials(0)%Id)
   call assign_material_id_to_media_matrix_coordinate(media, iEy, 3, 4, 3, simulationMaterials(0)%Id)

   expectedNumMeasurments = 4_SINGLE
   mediaPtr => media

   do iter = 1, 6
      sinpml(iter) = create_limit_t(0, 8, 0, 8, 0, 8, 10, 10, 10)
   end do

   dummyControl = create_control_flags(nEntradaRoot=nEntrada, mpidir=mpidir)

   call init_outputs(dummysgg, media, sinpml, tagNumbers, dummyBound, dummyControl, &
                     outputRequested, ThereAreWires)

   outputs => GetOutputs()

   test_err = test_err + assert_integer_equal(outputs(1)%outputID, MOVIE_PROBE_ID, 'Unexpected probe id')
   test_err = test_err + assert_integer_equal(outputs(1)%movieProbe%nPoints, expectedNumMeasurments, 'Unexpected number of measurements')
   test_err = test_err + assert_integer_equal(size(outputs(1)%movieProbe%xValueForTime), expectedNumMeasurments*BuffObse, 'Unexpected allocation size')
   test_err = test_err + assert_integer_equal(size(outputs(1)%movieProbe%timeStep), BuffObse, 'Unexpected timestep buffer size')

   expectedProbePath = trim(nEntrada)//wordSeparation//'movieProbe_BC_2_2_2__5_5_5'
   pdvFileName = trim(get_last_component(expectedProbePath))//pvdExtension

   test_err = test_err + assert_string_equal(outputs(1)%movieProbe%path, expectedProbePath, 'Unexpected path')
   test_err = test_err + assert_true(folder_exists(expectedProbePath), 'Movie folder do not exist')

   !Cleanup
   call remove_folder(test_folder, ios)

   err = test_err
end function

integer function test_update_movie_probe() bind(c) result(err)
! This test updates a movie probe with field values at a single timestep.
! It verifies that the probe correctly stores the measured values in the x, y,
! and z components for all points, and that the timestep buffer is properly populated.
   use output_m
   use outputTypes_m
   use testOutputUtils_m
   use FDETYPES_TOOLS
   use sggMethods_m
   use assertionTools_m
   use directoryUtils_m
   implicit none

   type(SGGFDTDINFO_t)              :: dummysgg
   type(sim_control_t)            :: dummyControl
   type(bounds_t)                 :: dummyBound
   type(solver_output_t), pointer :: outputs(:)

   type(media_matrices_t), target :: media
   type(media_matrices_t), pointer :: mediaPtr

   type(MediaData_t), allocatable, target :: simulationMaterials(:)
   type(MediaData_t), pointer     :: simulationMaterialsPtr(:)

   type(limit_t), target          :: sinpml_fullsize(6)
   type(limit_t), pointer         :: sinpml_fullsizePtr(:)

   type(taglist_t)                :: tagNumbers

   type(XYZlimit_t)               :: dummySweep(6)
   type(XYZlimit_t)               :: dummySinpmlSweep(6)
   type(XYZlimit_t)               :: allocationRange(6)

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

   character(len=14), parameter :: test_folder = 'testing_folder'
   character(len=11), parameter :: test_name = 'updateMovie'

   character(len=BUFSIZE) :: nEntrada
   integer :: ios

   nEntrada = join_path(test_folder, test_name)

   err = 1

   call sgg_init(dummysgg)
   call init_time_array(timeArray, nTimeSteps, dt)
   call sgg_set_tiempo(dummysgg, timeArray)
   call sgg_set_dt(dummysgg, dt)

   call init_simulation_material_list(simulationMaterials)
   simulationMaterialsPtr => simulationMaterials
   call sgg_set_NumMedia(dummysgg, size(simulationMaterials))
   call sgg_set_Med(dummysgg, simulationMaterialsPtr)

   dummySweep = create_xyz_limit_array(0, 0, 0, 6, 6, 6)
   call sgg_set_Sweep(dummysgg, dummySweep)
   dummySinpmlSweep = create_xyz_limit_array(1, 1, 1, 5, 5, 5)
   call sgg_set_SINPMLSweep(dummysgg, dummySinpmlSweep)
   call sgg_set_NumPlaneWaves(dummysgg, 1)
   allocationRange = create_xyz_limit_array(0, 0, 0, 6, 6, 6)
   call sgg_set_Alloc(dummysgg, allocationRange)

   movieObservable = create_movie_observation(2, 2, 2, 5, 5, 5, iCur)
   call sgg_add_observation(dummysgg, movieObservable)

   call create_geometry_media(media, 0, 8, 0, 8, 0, 8)

   call assign_material_id_to_media_matrix_coordinate(media, iEy, 3, 3, 3, simulationMaterials(0)%Id)
   call assign_material_id_to_media_matrix_coordinate(media, iEy, 4, 3, 3, simulationMaterials(0)%Id)
   call assign_material_id_to_media_matrix_coordinate(media, iEy, 4, 4, 3, simulationMaterials(0)%Id)
   call assign_material_id_to_media_matrix_coordinate(media, iEy, 3, 4, 3, simulationMaterials(0)%Id)

   expectedNumMeasurments = 4_SINGLE
   mediaPtr => media

   do iter = 1, 6
      sinpml_fullsize(iter) = create_limit_t(0, 8, 0, 8, 0, 8, 10, 10, 10)
   end do
   sinpml_fullsizePtr => sinpml_fullsize

   dummyControl = create_control_flags(nEntradaRoot=nEntrada, mpidir=mpidir)

   call init_outputs(dummysgg, media, sinpml_fullsize, tagNumbers, dummyBound, dummyControl, &
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

   !Cleanup
   call remove_folder(test_folder, ios)

   err = test_err
end function

integer function test_flush_movie_probe() bind(c) result(err)
! This test validates flushing movie probes to disk. It ensures that
! VTU files for each timestep and the PVD file are created, confirming that
! the temporal sequence of the movie probe is correctly serialized and saved.
   use output_m
   use outputTypes_m
   use testOutputUtils_m
   use FDETYPES_TOOLS
   use sggMethods_m
   use assertionTools_m
   use directoryUtils_m
   implicit none

   type(SGGFDTDINFO_t)              :: dummysgg
   type(sim_control_t)            :: dummyControl
   type(bounds_t)                 :: dummyBound
   type(solver_output_t), pointer :: outputs(:)

   type(media_matrices_t), target :: media
   type(media_matrices_t), pointer :: mediaPtr

   type(MediaData_t), allocatable, target :: simulationMaterials(:)
   type(MediaData_t), pointer     :: simulationMaterialsPtr(:)

   type(limit_t), target          :: sinpml_fullsize(6)
   type(limit_t), pointer         :: sinpml_fullsizePtr(:)

   type(taglist_t)                :: tagNumbers

   type(XYZlimit_t)               :: dummySweep(6)
   type(XYZlimit_t)               :: dummySinpmlSweep(6)
   type(XYZlimit_t)               :: allocationRange(6)

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

   character(len=14), parameter :: test_folder = 'testing_folder'
   character(len=10), parameter :: test_name = 'flushMovie'

   character(len=BUFSIZE) :: nEntrada
   character(len=BUFSIZE) :: expectedPath
   integer :: outputIdx
   integer :: ios

   nEntrada = join_path(test_folder, test_name)

   err = 1

   call sgg_init(dummysgg)
   call init_time_array(timeArray, nTimeSteps, dt)
   call sgg_set_tiempo(dummysgg, timeArray)
   call sgg_set_dt(dummysgg, dt)

   call init_simulation_material_list(simulationMaterials)
   simulationMaterialsPtr => simulationMaterials
   call sgg_set_NumMedia(dummysgg, size(simulationMaterials))
   call sgg_set_Med(dummysgg, simulationMaterialsPtr)

   dummySweep = create_xyz_limit_array(0, 0, 0, 6, 6, 6)
   call sgg_set_Sweep(dummysgg, dummySweep)
   dummySinpmlSweep = create_xyz_limit_array(1, 1, 1, 5, 5, 5)
   call sgg_set_SINPMLSweep(dummysgg, dummySinpmlSweep)
   call sgg_set_NumPlaneWaves(dummysgg, 1)
   allocationRange = create_xyz_limit_array(0, 0, 0, 6, 6, 6)
   call sgg_set_Alloc(dummysgg, allocationRange)

   movieCurrentObservable = create_movie_observation(2, 2, 2, 5, 5, 5, iCur)
   call sgg_add_observation(dummysgg, movieCurrentObservable)

   movieElectricXObservable = create_movie_observation(2, 2, 2, 5, 5, 5, iExC)
   call sgg_add_observation(dummysgg, movieElectricXObservable)

   movieMagneticYObservable = create_movie_observation(2, 2, 2, 5, 5, 5, iHyC)
   call sgg_add_observation(dummysgg, movieMagneticYObservable)

   call create_geometry_media(media, 0, 8, 0, 8, 0, 8)

   call assign_material_id_to_media_matrix_coordinate(media, iEx, 3, 3, 3, simulationMaterials(0)%Id)
   call assign_material_id_to_media_matrix_coordinate(media, iEy, 3, 3, 3, simulationMaterials(0)%Id)
   call assign_material_id_to_media_matrix_coordinate(media, iHy, 3, 3, 3, simulationMaterials(0)%Id)

   call assign_material_id_to_media_matrix_coordinate(media, iEx, 3, 4, 3, simulationMaterials(0)%Id)
   call assign_material_id_to_media_matrix_coordinate(media, iEy, 3, 4, 3, simulationMaterials(0)%Id)
   call assign_material_id_to_media_matrix_coordinate(media, iHy, 3, 4, 3, simulationMaterials(0)%Id)

   call assign_material_id_to_media_matrix_coordinate(media, iEx, 4, 4, 3, simulationMaterials(0)%Id)
   call assign_material_id_to_media_matrix_coordinate(media, iEy, 4, 4, 3, simulationMaterials(0)%Id)
   call assign_material_id_to_media_matrix_coordinate(media, iHy, 4, 4, 3, simulationMaterials(0)%Id)

   call assign_material_id_to_media_matrix_coordinate(media, iEx, 4, 3, 3, simulationMaterials(0)%Id)
   call assign_material_id_to_media_matrix_coordinate(media, iEy, 4, 3, 3, simulationMaterials(0)%Id)
   call assign_material_id_to_media_matrix_coordinate(media, iHy, 4, 3, 3, simulationMaterials(0)%Id)

   expectedNumMeasurments = 4_SINGLE
   mediaPtr => media

   do iter = 1, 6
      sinpml_fullsize(iter) = create_limit_t(0, 8, 0, 8, 0, 8, 10, 10, 10)
   end do
   sinpml_fullsizePtr => sinpml_fullsize

   dummyControl = create_control_flags(nEntradaRoot=nEntrada, mpidir=mpidir)

   call init_outputs(dummysgg, media, sinpml_fullsize, tagNumbers, dummyBound, dummyControl, &
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

   call close_outputs()

   call remove_folder(test_folder, ios)

   err = test_err
end function

integer function test_init_frequency_slice_probe() bind(c) result(err)
! This test initializes a frequency slice probe over a 3D region (2,2,2) to (5,5,5).
! It verifies that the probe is correctly set up, that the expected number of measurement
! points and frequencies are allocated, and that the output folder and PVD file exist.
   use output_m
   use outputTypes_m
   use testOutputUtils_m
   use FDETYPES_TOOLS
   use sggMethods_m
   use assertionTools_m
   use directoryUtils_m
   implicit none

   type(SGGFDTDINFO_t)              :: dummysgg
   type(sim_control_t)            :: dummyControl
   type(bounds_t)                 :: dummyBound
   type(solver_output_t), pointer :: outputs(:)

   type(media_matrices_t), target :: media
   type(media_matrices_t), pointer :: mediaPtr

   type(MediaData_t), allocatable, target :: simulationMaterials(:)
   type(MediaData_t), pointer     :: simulationMaterialsPtr(:)

   type(taglist_t)                :: tagNumbers

   type(limit_t), target          :: sinpml_fullsize(6)
   type(limit_t), pointer         :: sinpml_fullsizePtr(:)

   type(XYZlimit_t)               :: dummySweep(6)
   type(XYZlimit_t)               :: dummySinpmlSweep(6)
   type(XYZlimit_t)               :: allocationRange(6)

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


   character(len=14), parameter :: test_folder = 'testing_folder'
   character(len=13), parameter :: test_name = 'initFrequency'

   character(len=BUFSIZE) :: nEntrada
   character(len=BUFSIZE) :: expectedProbePath
   character(len=BUFSIZE) :: pdvFileName
   integer :: ios

   nEntrada = join_path(test_folder, test_name)
   err = 1

   call sgg_init(dummysgg)

   call init_time_array(timeArray, nTimeSteps, dt)
   call sgg_set_tiempo(dummysgg, timeArray)
   call sgg_set_dt(dummysgg, dt)

   call init_simulation_material_list(simulationMaterials)
   simulationMaterialsPtr => simulationMaterials
   call sgg_set_NumMedia(dummysgg, size(simulationMaterials))
   call sgg_set_Med(dummysgg, simulationMaterialsPtr)

   dummySweep = create_xyz_limit_array(0, 0, 0, 6, 6, 6)
   call sgg_set_Sweep(dummysgg, dummySweep)
   dummySinpmlSweep = create_xyz_limit_array(1, 1, 1, 5, 5, 5)
   call sgg_set_SINPMLSweep(dummysgg, dummySinpmlSweep)
   call sgg_set_NumPlaneWaves(dummysgg, 1)
   allocationRange = create_xyz_limit_array(0, 0, 0, 6, 6, 6)
   call sgg_set_Alloc(dummysgg, allocationRange)

   frequencySliceObservation = create_frequency_slice_observation(2, 2, 2, 5, 5, 5, iCur)
   call sgg_add_observation(dummysgg, frequencySliceObservation)

   expectedTotalFrequnecies = 6_SINGLE

   call create_geometry_media(media, 0, 8, 0, 8, 0, 8)

   call assign_material_id_to_media_matrix_coordinate(media, iEy, 3, 3, 3, simulationMaterials(0)%Id)
   call assign_material_id_to_media_matrix_coordinate(media, iEy, 4, 3, 3, simulationMaterials(0)%Id)
   call assign_material_id_to_media_matrix_coordinate(media, iEy, 4, 4, 3, simulationMaterials(0)%Id)
   call assign_material_id_to_media_matrix_coordinate(media, iEy, 3, 4, 3, simulationMaterials(0)%Id)

   expectedNumMeasurments = 4_SINGLE
   mediaPtr => media

   do iter = 1, 6
      sinpml_fullsize(iter) = create_limit_t(0, 8, 0, 8, 0, 8, 10, 10, 10)
   end do
   sinpml_fullsizePtr => sinpml_fullsize

   dummyControl = create_control_flags(nEntradaRoot=nEntrada, mpidir=mpidir)

   call init_outputs(dummysgg, media, sinpml_fullsize, tagNumbers, dummyBound, dummyControl, &
                     outputRequested, ThereAreWires)

   outputs => GetOutputs()

   test_err = test_err + assert_integer_equal(outputs(1)%outputID, FREQUENCY_SLICE_PROBE_ID, 'Unexpected probe id')

   test_err = test_err + assert_integer_equal(outputs(1)%frequencySliceProbe%nFreq, 6, 'Unexpected number of frequencies')

   test_err = test_err + assert_integer_equal(outputs(1)%frequencySliceProbe%nPoints, expectedNumMeasurments, 'Unexpected number of measurements')

   test_err = test_err + assert_integer_equal( &
              size(outputs(1)%frequencySliceProbe%xValueForFreq), &
              expectedNumMeasurments*expectedTotalFrequnecies, 'Unexpected allocation size')

   test_err = test_err + assert_integer_equal( &
              size(outputs(1)%frequencySliceProbe%frequencySlice), &
              expectedTotalFrequnecies, 'Unexpected frequency count')

   expectedProbePath = trim(nEntrada)//wordSeparation//'frequencySliceProbe_BC_2_2_2__5_5_5'
   pdvFileName = trim(get_last_component(expectedProbePath))//pvdExtension

   test_err = test_err + assert_string_equal(outputs(1)%frequencySliceProbe%path, expectedProbePath, 'Unexpected path')
   test_err = test_err + assert_true(folder_exists(expectedProbePath), 'Frequency Slice folder do not exist')

   call remove_folder(test_folder, ios)

   err = test_err
end function

integer function test_update_frequency_slice_probe() bind(c) result(err)
! This test updates a frequency slice probe with field gradients.
! It checks that no current is detected along the X-axis for H-field gradients
! and verifies the correct relation between Y and Z values (Y = -Z), ensuring
! correct handling of field distributions across the frequency slice.
   use output_m
   use outputTypes_m
   use testOutputUtils_m
   use FDETYPES_TOOLS
   use sggMethods_m
   use assertionTools_m
   use directoryUtils_m
   implicit none

   type(SGGFDTDINFO_t)              :: dummysgg
   type(sim_control_t)            :: dummyControl
   type(bounds_t)                 :: dummyBound
   type(solver_output_t), pointer :: outputs(:)

   type(media_matrices_t), target :: media
   type(media_matrices_t), pointer :: mediaPtr

   type(MediaData_t), allocatable, target :: simulationMaterials(:)
   type(MediaData_t), pointer     :: simulationMaterialsPtr(:)

   type(limit_t), target          :: sinpml_fullsize(6)
   type(limit_t), pointer         :: sinpml_fullsizePtr(:)

   type(taglist_t)                :: tagNumbers

   type(XYZlimit_t)               :: dummySweep(6)
   type(XYZlimit_t)               :: dummySinpmlSweep(6)
   type(XYZlimit_t)               :: allocationRange(6)

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


   character(len=14), parameter :: test_folder = 'testing_folder'
   character(len=13), parameter :: test_name = 'initFrequency'

   character(len=BUFSIZE) :: nEntrada
   integer :: ios

   nEntrada = join_path(test_folder, test_name)
   
   err = 1

   call sgg_init(dummysgg)
   call init_time_array(timeArray, nTimeSteps, dt)
   call sgg_set_tiempo(dummysgg, timeArray)
   call sgg_set_dt(dummysgg, dt)

   call init_simulation_material_list(simulationMaterials)
   simulationMaterialsPtr => simulationMaterials
   call sgg_set_NumMedia(dummysgg, size(simulationMaterials))
   call sgg_set_Med(dummysgg, simulationMaterialsPtr)

   dummySweep = create_xyz_limit_array(0, 0, 0, 6, 6, 6)
   call sgg_set_Sweep(dummysgg, dummySweep)
   dummySinpmlSweep = create_xyz_limit_array(1, 1, 1, 5, 5, 5)
   call sgg_set_SINPMLSweep(dummysgg, dummySinpmlSweep)
   call sgg_set_NumPlaneWaves(dummysgg, 1)
   allocationRange = create_xyz_limit_array(0, 0, 0, 6, 6, 6)
   call sgg_set_Alloc(dummysgg, allocationRange)

   frequencySliceObservation = create_frequency_slice_observation(2, 2, 2, 5, 5, 5, iCur)
   call sgg_add_observation(dummysgg, frequencySliceObservation)

   call create_geometry_media(media, 0, 8, 0, 8, 0, 8)

   call assign_material_id_to_media_matrix_coordinate(media, iEy, 3, 3, 3, simulationMaterials(0)%Id)
   call assign_material_id_to_media_matrix_coordinate(media, iEy, 4, 3, 3, simulationMaterials(0)%Id)
   call assign_material_id_to_media_matrix_coordinate(media, iEy, 4, 4, 3, simulationMaterials(0)%Id)
   call assign_material_id_to_media_matrix_coordinate(media, iEy, 3, 4, 3, simulationMaterials(0)%Id)
   expectedNumberFrequencies = 6_SINGLE
   expectedNumMeasurments = 4_SINGLE
   mediaPtr => media

   do iter = 1, 6
      sinpml_fullsize(iter) = create_limit_t(0, 8, 0, 8, 0, 8, 10, 10, 10)
   end do
   sinpml_fullsizePtr => sinpml_fullsize

   dummyControl = create_control_flags(nEntradaRoot=nEntrada, mpidir=mpidir)

   call init_outputs(dummysgg, media, sinpml_fullsize, tagNumbers, dummyBound, dummyControl, &
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
                                -1.0_RKIND*outputs(1)%frequencySliceProbe%zValueForFreq, errormessage='Unequal values for Y and -Z')


   call remove_folder(test_folder, ios)

   err = test_err
end function
