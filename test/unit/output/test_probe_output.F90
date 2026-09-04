integer function test_init_point_probe() bind(c) result(err)
   ! Verifies point probes publish one flat text file without metadata sidecars.
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
   character(len=*), parameter :: test_folder = 'testing folder'
   character(len=*), parameter :: test_name = 'nested output/initPointProbeTest'

   ! Local variables
   character(len=BUFSIZE) :: testPath, nEntrada
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
   testPath = join_path(get_temp_folder(), test_folder)
   nEntrada = join_path(testPath, test_name)

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
   test_err = test_err + assert_string_equal(outputs(1)%pointProbe%artifacts(1)%relative_path, &
                                              expectedDataPath, 'Declared payload path changed')
   test_err = test_err + assert_true(file_exists(expectedDataPath), 'Time data file do not exist')
   test_err = test_err + assert_true(.not. file_exists(trim(expectedProbePath)//'.json'), &
                                      'Point probe created a descriptor sidecar')
   test_err = test_err + assert_true(.not. folder_exists(expectedProbePath), &
                                     'Point probe created a probe-specific folder')

   ! Cleanup
   call remove_folder(testPath, ios)
   deallocate (sgg%Observation, outputs)

   err = test_err
end function

integer function test_init_point_probe_with_incident() bind(c) result(err)
   ! Verifies incident point probes preserve their text header without a binary sidecar.
   use FDETYPES_m
   use FDETYPES_TOOLS
   use outputTypes_m, only: point_probe_output_t, domain_t, cell_coordinate_t, TIME_DOMAIN, OUTPUT_ARTIFACT_UNDEFINED
   use pointProbeOutput_m, only: init_point_probe_output
   use assertionTools_m, only: assert_true, assert_string_equal
   use directoryUtils_m, only: delete_file, file_exists, join_path
   use testOutputUtils_m, only: get_temp_folder
   implicit none

   type(point_probe_output_t) :: probe
   type(domain_t) :: domain
   type(cell_coordinate_t) :: coordinates
   character(len=4096) :: header, path
   integer :: ios, unit

   err = 0
   path = join_path(get_temp_folder(), 'incidentPointProbe')
   domain%domainType = TIME_DOMAIN
   coordinates = cell_coordinate_t(1, 1, 1)
   call init_point_probe_output(probe, coordinates, iEx, domain, path, 3, 0.1_RKIND_tiempo, .true.)

   err = err + assert_true(probe%hasIncident .and. allocated(probe%incidentForTime), &
                           'Incident point probe did not allocate incident samples')
   err = err + assert_true(size(probe%artifacts) == 1 .and. &
                           probe%artifacts(1)%kind /= OUTPUT_ARTIFACT_UNDEFINED, &
                           'Incident point probe declared a non-text sidecar')
   open (newunit=unit, file=probe%filePathTime, status='old', action='read', iostat=ios)
   read (unit, '(A)', iostat=ios) header
   close (unit)
   err = err + assert_true(ios == 0 .and. trim(header) == 't field incident', 'Incident text header is incorrect')
   err = err + assert_true(.not. file_exists(trim(path)//'_Ex_1_1_1_tm.bin'), &
                           'Incident point probe created a binary sidecar')
   call delete_file(probe%filePathTime, ios)
end function test_init_point_probe_with_incident

integer function test_scalar_probe_has_no_manifest() bind(c) result(err)
   ! Verifies a scalar-only run does not publish descriptors or a root manifest.
   use FDETYPES_m
   use FDETYPES_TOOLS
   use output_m
   use testOutputUtils_m
   use sggMethods_m
   use assertionTools_m
   use directoryUtils_m
   implicit none

   type(SGGFDTDINFO_t) :: sgg
   type(sim_control_t) :: control
   type(bounds_t) :: bounds
   type(media_matrices_t) :: media
   type(limit_t), allocatable :: sinpml(:)
   type(Obses_t) :: probe
   type(MediaData_t), allocatable, target :: materials(:)
   type(MediaData_t), pointer :: materials_ptr(:)
   type(taglist_t) :: tag_numbers
   real(kind=RKIND_tiempo), pointer :: time_array(:)
   character(len=BUFSIZE) :: path, probe_path
   logical :: observations_exist, wires_exist

   err = 0
   wires_exist = .false.
   path = join_path(get_temp_folder(), 'rootManifest')
   probe_path = trim(path)//'_pointProbe_Ex_4_4_4'
   call sgg_init(sgg)
   call init_time_array(time_array, 2_SINGLE, 0.1_RKIND_tiempo)
   call sgg_set_tiempo(sgg, time_array)
   call sgg_set_dt(sgg, 0.1_RKIND_tiempo)
   call init_simulation_material_list(materials)
   materials_ptr => materials
   call sgg_set_Med(sgg, materials_ptr)
   probe = create_point_probe_observation(4, 4, 4)
   call sgg_add_observation(sgg, probe)
   control = create_control_flags(nEntradaRoot=path, mpidir=3, size=1)

   call init_outputs(sgg, media, sinpml, tag_numbers, bounds, control, observations_exist, wires_exist)

   call close_outputs()
   err = err + assert_true(.not. file_exists(trim(path)//'_output_manifest.json'), &
                           'Scalar-only run published a root manifest')
   err = err + assert_true(file_exists(trim(probe_path)//'_tm.dat'), &
                           'Scalar-only run did not retain its flat data file')
   err = err + assert_true(.not. file_exists(trim(probe_path)//'.json'), &
                           'Scalar-only run published a probe descriptor')
   call delete_outputs(0)
   err = err + assert_true(.not. file_exists(trim(probe_path)//'_tm.dat'), &
                           'Direct output cleanup retained scalar data')
end function test_scalar_probe_has_no_manifest

integer function test_line_probe_integral() bind(c) result(err)
   ! Verifies the legacy signed E.dl line-integral convention in isolation.
   use FDETYPES_m, only: RKIND, direction_t, iEx, iEy, iEz
   use outputTypes_m, only: field_data_t
   use lineProbeOutput_m, only: calculate_line_integral
   use assertionTools_m, only: assert_real_equal
   implicit none

   type(direction_t) :: segments(3), reversed_segments(3)
   type(field_data_t) :: electric_field
   real(kind=RKIND), target :: ex(3, 3, 3), ey(3, 3, 3), ez(3, 3, 3)
   real(kind=RKIND), target :: dx(3), dy(3), dz(3)
   real(kind=RKIND) :: value

   err = 0
   ex = 0.0_RKIND
   ey = 0.0_RKIND
   ez = 0.0_RKIND
   dx = 1.0_RKIND
   dy = 1.0_RKIND
   dz = 1.0_RKIND

   ex(1, 1, 1) = 2.0_RKIND
   ey(2, 1, 1) = 3.0_RKIND
   ez(2, 2, 1) = 4.0_RKIND
   dx(1) = 0.5_RKIND
   dy(1) = 2.0_RKIND
   dz(1) = 1.5_RKIND

   segments(1) = direction_t(1, 1, 1, iEx)
   segments(2) = direction_t(2, 1, 1, -iEy)
   segments(3) = direction_t(2, 2, 1, iEz)
   reversed_segments = segments
   reversed_segments%orientation = -reversed_segments%orientation

   electric_field%x => ex
   electric_field%y => ey
   electric_field%z => ez
   electric_field%deltaX => dx
   electric_field%deltaY => dy
   electric_field%deltaZ => dz

   value = calculate_line_integral(segments, electric_field)
   err = err + assert_real_equal(value, 1.0_RKIND - 6.0_RKIND + 6.0_RKIND, 1.0e-6_RKIND, &
                                 'Mixed-direction line integral is incorrect')
   value = calculate_line_integral(reversed_segments, electric_field)
   err = err + assert_real_equal(value, -1.0_RKIND, 1.0e-6_RKIND, &
                                 'Reversed line orientation did not reverse the integral sign')
end function test_line_probe_integral

integer function test_line_probe_empty_path() bind(c) result(err)
   ! Verifies an empty line remains a valid zero-sample probe.
   use FDETYPES_m, only: RKIND, RKIND_tiempo, direction_t
   use outputTypes_m, only: line_probe_output_t, field_data_t, domain_t, TIME_DOMAIN
   use lineProbeOutput_m, only: init_line_probe_output, update_line_probe_output
   use assertionTools_m, only: assert_integer_equal
   use directoryUtils_m, only: delete_file
   use testOutputUtils_m, only: get_temp_folder
   use directoryUtils_m, only: join_path
   implicit none

   type(line_probe_output_t) :: probe
   type(domain_t) :: domain
   type(direction_t), allocatable :: segments(:)
   type(field_data_t) :: electric_field
   real(kind=RKIND), target :: field(1, 1, 1), spacing(1)
   character(len=4096) :: path
   integer :: ios

   err = 0
   path = join_path(get_temp_folder(), 'line-probe-empty')
   allocate (segments(0))
   domain%domainType = TIME_DOMAIN
   call init_line_probe_output(probe, segments, domain, path)
   field = 0.0_RKIND
   spacing = 1.0_RKIND
   electric_field%x => field
   electric_field%y => field
   electric_field%z => field
   electric_field%deltaX => spacing
   electric_field%deltaY => spacing
   electric_field%deltaZ => spacing
   call update_line_probe_output(probe, 0.0_RKIND_tiempo, electric_field)
   err = err + assert_integer_equal(probe%nTime, 0, 'Empty line probe recorded a fabricated sample')
   call delete_file(trim(path)//'_tm.dat', ios)
end function test_line_probe_empty_path

integer function test_line_probe_dat_output() bind(c) result(err)
   ! Verifies line-probe flush retains every sample in one text file.
   use FDETYPES_m, only: RKIND, RKIND_tiempo, direction_t, iEx
   use outputTypes_m, only: line_probe_output_t, field_data_t, domain_t, TIME_DOMAIN
   use lineProbeOutput_m, only: init_line_probe_output, update_line_probe_output, flush_line_probe_output
   use assertionTools_m, only: assert_integer_equal, assert_true
   use directoryUtils_m, only: delete_file, file_exists
   use testOutputUtils_m, only: get_temp_folder
   use directoryUtils_m, only: join_path
   implicit none

   type(line_probe_output_t) :: probe
   type(domain_t) :: domain
   type(direction_t) :: segments(1)
   type(field_data_t) :: electric_field
   real(kind=RKIND), target :: ex(1, 1, 1), ey(1, 1, 1), ez(1, 1, 1), spacing(1)
   character(len=4096) :: path
   integer :: ios, text_unit, text_records
   character(len=128) :: line

   err = 0
   path = join_path(get_temp_folder(), 'line-probe-artifacts')
   domain%domainType = TIME_DOMAIN
   segments(1) = direction_t(1, 1, 1, iEx)
   call init_line_probe_output(probe, segments, domain, path)
   ex = 2.0_RKIND
   ey = 0.0_RKIND
   ez = 0.0_RKIND
   spacing = 0.5_RKIND
   electric_field%x => ex
   electric_field%y => ey
   electric_field%z => ez
   electric_field%deltaX => spacing
   electric_field%deltaY => spacing
   electric_field%deltaZ => spacing
   call update_line_probe_output(probe, 0.0_RKIND_tiempo, electric_field)
   call update_line_probe_output(probe, 0.1_RKIND_tiempo, electric_field)
   call flush_line_probe_output(probe)

   text_records = 0
   open (newunit=text_unit, file=trim(path)//'_tm.dat', status='old', action='read', iostat=ios)
   do while (ios == 0)
      read (text_unit, '(A)', iostat=ios) line
      if (ios == 0) text_records = text_records + 1
   end do
   close (text_unit)
   err = err + assert_integer_equal(text_records, 3, 'Line text artifact lost its header or samples')
   err = err + assert_true(.not. file_exists(trim(path)//'_tm.bin'), &
                           'Line probe created a binary sidecar')
   call delete_file(trim(path)//'_tm.dat', ios)
end function test_line_probe_dat_output

integer function test_line_probe_shared_interface_owner() bind(c) result(err)
   use FDETYPES_m, only: direction_t, xyzlimit_t, iEx
   use lineProbeOutput_m, only: line_segment_is_local
   use assertionTools_m, only: assert_true
   implicit none

   type(direction_t) :: segment
   type(xyzlimit_t) :: lower_rank_sweeps(6), upper_rank_sweeps(6)

   err = 0
   lower_rank_sweeps = xyzlimit_t(0, 4, 0, 4, 0, 2)
   upper_rank_sweeps = xyzlimit_t(0, 4, 0, 4, 2, 4)
   segment = direction_t(1, 1, 2, iEx)
   err = err + assert_true(.not. line_segment_is_local(segment, lower_rank_sweeps, 0, 2), &
                           'Lower rank retained a shared Ex interface segment')
   err = err + assert_true(line_segment_is_local(segment, upper_rank_sweeps, 1, 2), &
                           'Higher rank did not own a shared Ex interface segment')
end function test_line_probe_shared_interface_owner

integer function test_nested_output_path() bind(c) result(err)
   ! Verifies nested output directories are created and removed correctly.
   use directoryUtils_m, only: create_file_with_path, file_exists, remove_folder
   use assertionTools_m, only: assert_integer_equal, assert_true
   use testOutputUtils_m, only: get_temp_folder
   use directoryUtils_m, only: join_path
   implicit none

   character(len=4096) :: path
   integer :: ios

   err = 0
   path = join_path(get_temp_folder(), 'testing folder/nested output/result.dat')
   call create_file_with_path(path, ios)
   err = err + assert_integer_equal(ios, 0, 'Nested path creation failed')
   err = err + assert_true(file_exists(path), 'Nested output file does not exist')
   call remove_folder(join_path(get_temp_folder(), 'testing folder'), ios)
   err = err + assert_integer_equal(ios, 0, 'Nested path cleanup failed')
end function

integer function test_output_artifact_contract() bind(c) result(err)
   ! Verifies binary artifact encoding and declared lifecycle values.
   use outputTypes_m, only: output_artifact_t, output_lifecycle_t, OUTPUT_ARTIFACT_BINARY, &
                            OUTPUT_LIFECYCLE_DECLARED, BINARY_ENDIAN_LITTLE, BINARY_COMPLEX_REAL_IMAG
   use assertionTools_m, only: assert_integer_equal
   implicit none

   type(output_artifact_t) :: artifact
   type(output_lifecycle_t) :: lifecycle

   err = 0
   artifact%kind = OUTPUT_ARTIFACT_BINARY
   artifact%byte_order = BINARY_ENDIAN_LITTLE
   artifact%complex_representation = BINARY_COMPLEX_REAL_IMAG
   lifecycle%state = OUTPUT_LIFECYCLE_DECLARED

   err = err + assert_integer_equal(artifact%kind, OUTPUT_ARTIFACT_BINARY, 'Binary artifact kind')
   err = err + assert_integer_equal(artifact%byte_order, BINARY_ENDIAN_LITTLE, 'Binary byte order')
   err = err + assert_integer_equal(artifact%complex_representation, BINARY_COMPLEX_REAL_IMAG, &
                                    'Complex representation')
   err = err + assert_integer_equal(lifecycle%state, OUTPUT_LIFECYCLE_DECLARED, 'Declared lifecycle state')
end function

integer function test_portable_binary_output() bind(c) result(err)
   ! Verifies little-endian real64 output and rejects unsupported byte order.
   use, intrinsic :: iso_fortran_env, only: int8, real64
   use outputBinary_m, only: append_binary_real64, validate_binary_layout, BINARY_WRITER_SUCCESS, &
                              BINARY_WRITER_INVALID_LAYOUT
   use outputTypes_m, only: output_artifact_t, OUTPUT_ARTIFACT_BINARY, BINARY_ENDIAN_LITTLE, &
                             BINARY_ENDIAN_BIG, BINARY_NUMERIC_REAL64, BINARY_COMPLEX_UNSPECIFIED
   use assertionTools_m, only: assert_integer_equal, assert_true
   use directoryUtils_m, only: remove_folder
   use testOutputUtils_m, only: get_temp_folder
   use directoryUtils_m, only: join_path
   implicit none

   type(output_artifact_t) :: artifact
   integer(int8) :: bytes(16)
   integer :: expected_bytes(16), ios, status, unit
   character(len=4096) :: folder, path

   err = 0
   folder = join_path(get_temp_folder(), 'testing binary')
   path = join_path(folder, 'payload.bin')
   artifact%kind = OUTPUT_ARTIFACT_BINARY
   artifact%byte_order = BINARY_ENDIAN_LITTLE
   artifact%numeric_representation = BINARY_NUMERIC_REAL64
   artifact%complex_representation = BINARY_COMPLEX_UNSPECIFIED
   artifact%record_bytes = 16
   artifact%component_order = 'time,value'
   expected_bytes = [0, 0, 0, 0, 0, 0, -16, 63, 0, 0, 0, 0, 0, 0, 4, -64]

   call validate_binary_layout(artifact, status)
   err = err + assert_integer_equal(status, BINARY_WRITER_SUCCESS, 'Valid real64 layout was rejected')
   call append_binary_real64(path, artifact, [1.0_real64, -2.5_real64], status)
   err = err + assert_integer_equal(status, BINARY_WRITER_SUCCESS, 'Portable binary write failed')

   open (newunit=unit, file=path, access='stream', form='unformatted', status='old', &
         action='read', iostat=ios)
   err = err + assert_integer_equal(ios, 0, 'Cannot read portable binary payload')
   if (ios == 0) then
      read (unit, iostat=ios) bytes
      close (unit)
      err = err + assert_integer_equal(ios, 0, 'Cannot read portable binary bytes')
      err = err + assert_true(all(int(bytes) == expected_bytes), &
                              'Payload is not little-endian IEEE real32')
   end if

   artifact%byte_order = BINARY_ENDIAN_BIG
   call validate_binary_layout(artifact, status)
   err = err + assert_integer_equal(status, BINARY_WRITER_INVALID_LAYOUT, 'Unsupported byte order was accepted')
   call remove_folder(folder, ios)
end function

integer function test_declared_output_artifacts() bind(c) result(err)
   ! Verifies output artifact kinds and relative paths are declared.
   use outputTypes_m, only: probe_metadata_t, OUTPUT_ARTIFACT_TEXT, OUTPUT_ARTIFACT_GEOMETRY, &
                            declare_output_artifacts
   use assertionTools_m, only: assert_integer_equal, assert_string_equal, assert_true
   implicit none

   type(probe_metadata_t) :: metadata

   err = 0
   call declare_output_artifacts(metadata, [character(len=16) :: 'probe_tm.dat', 'geometry.vtu'], &
                                 [OUTPUT_ARTIFACT_TEXT, OUTPUT_ARTIFACT_GEOMETRY])

   err = err + assert_true(allocated(metadata%artifacts), 'Artifacts were not declared')
   err = err + assert_integer_equal(size(metadata%artifacts), 2, 'Unexpected artifact count')
   err = err + assert_integer_equal(metadata%artifacts(1)%kind, OUTPUT_ARTIFACT_TEXT, 'Text kind was not retained')
   err = err + assert_string_equal(metadata%artifacts(1)%relative_path, 'probe_tm.dat', 'Text path was not retained')
   err = err + assert_integer_equal(metadata%artifacts(2)%kind, OUTPUT_ARTIFACT_GEOMETRY, 'Geometry kind was not retained')
   err = err + assert_string_equal(metadata%artifacts(2)%relative_path, 'geometry.vtu', 'Geometry path was not retained')
end function

integer function test_output_lifecycle_contract() bind(c) result(err)
   ! Verifies lifecycle terminal states and metadata completeness.
   use outputTypes_m, only: probe_metadata_t, output_artifact_t, output_lifecycle_is_terminal, &
                            probe_metadata_is_complete, OUTPUT_ARTIFACT_BINARY, OUTPUT_LIFECYCLE_DECLARED, &
                            OUTPUT_LIFECYCLE_ACTIVE, OUTPUT_LIFECYCLE_COMPLETE, OUTPUT_LIFECYCLE_FAILED
   use assertionTools_m, only: assert_true
   implicit none

   type(probe_metadata_t) :: metadata

   err = 0
   metadata%probe_id = 'lifecycle-001'
   metadata%quantity = 'Ex'
   allocate (output_artifact_t :: metadata%artifacts(1))
   metadata%artifacts(1)%kind = OUTPUT_ARTIFACT_BINARY
   metadata%artifacts(1)%relative_path = 'lifecycle.bin'

   metadata%lifecycle%state = OUTPUT_LIFECYCLE_DECLARED
   err = err + assert_true(.not. output_lifecycle_is_terminal(metadata%lifecycle), 'Declared lifecycle is terminal')
   err = err + assert_true(.not. probe_metadata_is_complete(metadata), 'Declared metadata is complete')

   metadata%lifecycle%state = OUTPUT_LIFECYCLE_ACTIVE
   err = err + assert_true(.not. output_lifecycle_is_terminal(metadata%lifecycle), 'Active lifecycle is terminal')

   metadata%lifecycle%state = OUTPUT_LIFECYCLE_COMPLETE
   err = err + assert_true(output_lifecycle_is_terminal(metadata%lifecycle), 'Complete lifecycle is not terminal')
   err = err + assert_true(probe_metadata_is_complete(metadata), 'Complete zero-sample metadata is incomplete')

   metadata%lifecycle%state = OUTPUT_LIFECYCLE_FAILED
   err = err + assert_true(output_lifecycle_is_terminal(metadata%lifecycle), 'Failed lifecycle is not terminal')
   err = err + assert_true(.not. probe_metadata_is_complete(metadata), 'Failed metadata is complete')
end function

integer function test_output_serial_distributed_equivalence() bind(c) result(err)
   ! Verifies serial and distributed artifacts have equivalent coverage.
   use FDETYPES_m, only: iEx, limit_t
   use outputTypes_m, only: cell_coordinate_t, output_artifact_t, OUTPUT_ARTIFACT_BINARY
   use outputDecomposition_m, only: output_partition_t, build_output_partition, OUTPUT_PARTITION_SUCCESS
   use assertionTools_m, only: assert_integer_equal, assert_string_equal, assert_true
   implicit none

   type(output_artifact_t) :: serial_artifact, distributed_artifact
   type(cell_coordinate_t) :: lower_bound, upper_bound
   type(limit_t) :: global_bounds, local_sweep
   type(output_partition_t) :: partition
   integer :: coverage(0:5), rank, z, status

   err = 0
   serial_artifact%kind = OUTPUT_ARTIFACT_BINARY
   serial_artifact%relative_path = 'probe.bin'
   distributed_artifact = serial_artifact
   err = err + assert_integer_equal(distributed_artifact%kind, serial_artifact%kind, 'Artifact kind differs')
   err = err + assert_string_equal(distributed_artifact%relative_path, serial_artifact%relative_path, &
                                   'Artifact path differs')

   lower_bound = cell_coordinate_t(0, 0, 0)
   upper_bound = cell_coordinate_t(0, 0, 5)
   global_bounds = limit_t(0, 0, 0, 0, 0, 5, 1, 1, 6)
   coverage = 0
   do rank = 0, 1
      if (rank == 0) then
         local_sweep = limit_t(0, 0, 0, 0, 0, 3, 1, 1, 4)
      else
         local_sweep = limit_t(0, 0, 0, 0, 3, 5, 1, 1, 3)
      end if
      call build_output_partition(lower_bound, upper_bound, global_bounds, local_sweep, iEx, rank, 2, partition, status)
      err = err + assert_integer_equal(status, OUTPUT_PARTITION_SUCCESS, 'Partition construction failed')
      do z = partition%local_lower%z, partition%local_upper%z
         coverage(z) = coverage(z) + 1
      end do
   end do
   err = err + assert_true(all(coverage == 1), 'Distributed partition does not match serial coverage')
end function

integer function test_volumetric_output_partition_attachment() bind(c) result(err)
   ! Verifies volumetric partitions attach to outputs and select serial fallback.
   use FDETYPES_m
   use FDETYPES_TOOLS, only: create_limit_t, create_control_flags, init_time_array, &
                             init_simulation_material_list, create_geometry_media, create_xyz_limit_array, create_tag_list
   use output_m, only: init_outputs, GetOutputs, GetOutputPartition, solver_output_t
   use outputDecomposition_m, only: output_partition_t, OUTPUT_PARTITION_SUCCESS
   use outputCollective_m, only: OUTPUT_PUBLICATION_ROOT_AGGREGATION
   use testOutputUtils_m, only: create_movie_observation, get_temp_folder
   use sggMethods_m, only: sgg_init, sgg_set_tiempo, sgg_set_dt, sgg_set_Med, sgg_set_NumMedia, &
                           sgg_set_Sweep, sgg_set_SINPMLSweep, sgg_set_NumPlaneWaves, sgg_set_Alloc, &
                           sgg_set_LineX, sgg_set_LineY, sgg_set_LineZ, sgg_add_observation
   use assertionTools_m, only: assert_integer_equal, assert_true
   use directoryUtils_m, only: delete_file, file_exists, join_path, remove_folder
   use testOutputUtils_m, only: get_temp_folder
   implicit none

   type(SGGFDTDINFO_t) :: sgg
   type(sim_control_t) :: control
   type(media_matrices_t) :: media
   type(limit_t) :: sinpml(6)
   type(bounds_t) :: bounds
   type(taglist_t) :: material_tags
   type(Obses_t) :: observation
   type(MediaData_t), allocatable, target :: materials(:)
   type(MediaData_t), pointer :: materials_ptr(:)
   type(XYZlimit_t) :: sweep(6)
   type(output_partition_t) :: partition
   type(solver_output_t), pointer :: outputs(:)
   real(kind=RKIND_tiempo), pointer :: time_array(:)
   real(kind=RKIND), pointer :: x_steps(:), y_steps(:), z_steps(:)
   logical :: observations_exist, wires_exist
   integer :: status, i, ios
   character(len=BUFSIZE) :: path

   err = 0
   path = join_path(get_temp_folder(), 'partitionAttachment')
   wires_exist = .false.
   call sgg_init(sgg)
   call init_time_array(time_array, 2_SINGLE, 0.1_RKIND_tiempo)
   call sgg_set_tiempo(sgg, time_array)
   call sgg_set_dt(sgg, 0.1_RKIND_tiempo)
   call init_simulation_material_list(materials)
   materials_ptr => materials
   call sgg_set_NumMedia(sgg, size(materials))
   call sgg_set_Med(sgg, materials_ptr)
   sweep = create_xyz_limit_array(0, 0, 0, 6, 6, 6)
   call sgg_set_Sweep(sgg, sweep)
   call sgg_set_SINPMLSweep(sgg, sweep)
   call sgg_set_NumPlaneWaves(sgg, 1)
   call sgg_set_Alloc(sgg, sweep)
   material_tags = create_tag_list(sweep)
   allocate (x_steps(0:8), y_steps(0:8), z_steps(0:8), source=1.0_RKIND)
   call sgg_set_LineX(sgg, x_steps)
   call sgg_set_LineY(sgg, y_steps)
   call sgg_set_LineZ(sgg, z_steps)
   do i = 1, size(sinpml)
      sinpml(i) = create_limit_t(0, 8, 0, 8, 0, 8, 9, 9, 9)
   end do
   call create_geometry_media(media, 0, 8, 0, 8, 0, 8)
   observation = create_movie_observation(2, 2, 2, 5, 5, 5, iCur)
   call sgg_add_observation(sgg, observation)
   control = create_control_flags(nEntradaRoot=path, mpidir=3, size=1)

   call init_outputs(sgg, media, sinpml, material_tags, bounds, control, observations_exist, wires_exist)
   call GetOutputPartition(1, partition, status)
   outputs => GetOutputs()

   err = err + assert_integer_equal(status, OUTPUT_PARTITION_SUCCESS, 'Volumetric partition was not retained')
   err = err + assert_true(partition%has_data, 'Serial volumetric partition has no data')
   err = err + assert_integer_equal(partition%global_lower%x, 2, 'Unexpected global lower x')
   err = err + assert_integer_equal(partition%global_upper%z, 5, 'Unexpected global upper z')
   err = err + assert_integer_equal(partition%local_lower%z, 2, 'Unexpected local lower z')
   err = err + assert_integer_equal(partition%local_upper%z, 5, 'Unexpected local upper z')
    err = err + assert_integer_equal(outputs(1)%movieProbe%publication%mode, OUTPUT_PUBLICATION_ROOT_AGGREGATION, &
                                     'Serial movie output did not select root aggregation fallback')
    err = err + assert_true(outputs(1)%movieProbe%publication%local_participates, &
                            'Serial movie output was excluded from publication')
   err = err + assert_true(.not. file_exists(trim(outputs(1)%movieProbe%filesPath)//'.json'), &
                           'Movie probe published a JSON descriptor during initialisation')
   call remove_folder(trim(path)//'_movieProbe_BC_2_2_2__5_5_5', ios)
end function

integer function test_update_point_probe() bind(c) result(err)
   ! Verifies time-frequency point probes honour their time window without decimating frequency updates.
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
   character(len=BUFSIZE) :: testPath, nEntrada

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
   integer :: i, ios

   ! Setup
   sep = get_path_separator()
   testPath = join_path(get_temp_folder(), test_folder)
   nEntrada = testPath//sep//test_name

   call sgg_init(sgg)
   call init_time_array(timeArray, nSteps, dt)
   call sgg_set_tiempo(sgg, timeArray)
   call sgg_set_dt(sgg, dt)

   probe = create_point_probe_observation(4, 4, 4)
   probe%InitialTime = 0.25_RKIND
   probe%FinalTime = 0.65_RKIND
   probe%TimeStep = 0.2_RKIND
   probe%FreqDomain = .true.
   probe%InitialFreq = 0.0_RKIND
   probe%FinalFreq = 0.0_RKIND
   probe%FreqStep = 0.0_RKIND
   call sgg_add_observation(sgg, probe)

   call init_simulation_material_list(materials)
   materialsPtr => materials
   call sgg_set_Med(sgg, materialsPtr)

   control = create_control_flags(mpidir=3, finaltimestep=nSteps - 2, nEntradaRoot=nEntrada, &
                                  wiresflavor='holland')
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
   do i = 0, 8
      dummyFields%Ex(4, 4, 4) = real(i, RKIND)
      call update_outputs(control, sgg%tiempo(i + 1), int(i, SINGLE), fields)
   end do
   outputs => GetOutputs()

   ! Assertions
   test_err = test_err + assert_integer_equal(outputs(1)%pointProbe%nTime, 2, 'Unexpected time sample count')
   test_err = test_err + assert_real_equal(outputs(1)%pointProbe%timeStep(1), 0.4_RKIND_tiempo, &
                                           1e-5_RKIND_tiempo, 'Unexpected timestep 1')
   test_err = test_err + assert_real_equal(outputs(1)%pointProbe%valueForTime(1), 4.0_RKIND, &
                                           1e-5_RKIND, 'Unexpected field 1')
   test_err = test_err + assert_real_equal(outputs(1)%pointProbe%timeStep(2), 0.6_RKIND_tiempo, &
                                           1e-5_RKIND_tiempo, 'Unexpected timestep 2')
   test_err = test_err + assert_real_equal(outputs(1)%pointProbe%valueForTime(2), 6.0_RKIND, &
                                           1e-5_RKIND, 'Unexpected field 2')
   test_err = test_err + assert_real_equal(real(outputs(1)%pointProbe%valueForFreq(1), RKIND), 3.6_RKIND, &
                                           1e-5_RKIND, 'Frequency accumulation was decimated')

   !Cleanup
   call remove_folder(testPath, ios)

   err = test_err
end function

integer function test_update_time_probe_ranges() bind(c) result(err)
   ! Verifies each scalar time output honours an explicit time window and sampling period.
   use FDETYPES_m
   use FDETYPES_TOOLS
   use output_m
   use outputTypes_m
   use testOutputUtils_m
   use sggMethods_m
   use assertionTools_m
   use directoryUtils_m
   use HollandWires_m, only: GetHwires
   use wiresHolland_constants_m, only: ThinWires_t
   implicit none

   character(len=22), parameter :: test_name = 'updateTimeProbeRanges'
   integer, parameter :: current_output = 1, charge_output = 2, bulk_output = 3, line_output = 4

   type(SGGFDTDINFO_t) :: sgg
   type(sim_control_t) :: control
   type(bounds_t) :: bounds
   type(media_matrices_t) :: media
   type(limit_t), allocatable :: sinpml(:)
   type(taglist_t) :: tagNumbers
   type(MediaData_t), allocatable, target :: materials(:)
   type(MediaData_t), pointer :: materialsPtr(:)
   type(observation_domain_t) :: domain
   type(observable_t) :: requestedOutput(1)
   type(Obses_t) :: observations(4)
   type(direction_t) :: lineSegments(1)
   type(XYZlimit_t) :: sweep(6)
   type(solver_output_t), pointer :: outputs(:)
   type(ThinWires_t), pointer :: wires
   type(dummyFields_t), target :: dummyFields
   type(fields_reference_t) :: fields
   real(kind=RKIND_tiempo), pointer :: timeArray(:)
   real(kind=RKIND_tiempo), parameter :: dt = 0.1_RKIND_tiempo
   real(kind=RKIND), target :: wireField
   character(len=BUFSIZE) :: testPath, nEntrada
   integer(kind=SINGLE), parameter :: nSteps = 100_SINGLE
   integer(kind=SINGLE) :: test_err
   integer :: i, ios
   logical :: outputRequested, hasWires

   test_err = 0
   testPath = join_path(get_temp_folder(), 'testing_folder')
   nEntrada = join_path(testPath, test_name)

   call sgg_init(sgg)
   call init_time_array(timeArray, nSteps, dt)
   call sgg_set_tiempo(sgg, timeArray)
   call sgg_set_dt(sgg, dt)
   sweep = create_xyz_limit_array(1, 1, 1, 5, 5, 5)
   call sgg_set_Sweep(sgg, sweep)

   call init_simulation_material_list(materials)
   materialsPtr => materials
   call sgg_set_Med(sgg, materialsPtr)

   call initialize_observation_time_domain(domain, 0.25_RKIND, 0.65_RKIND, 0.2_RKIND)

   requestedOutput(1) = create_observable(3, 3, 3, 3, 3, 3, iJx)
   requestedOutput(1)%Node = 1
   call set_observation(observations(current_output), requestedOutput, 'wireCurrentRange', domain, '')

   requestedOutput(1) = create_observable(3, 3, 3, 3, 3, 3, iQx)
   requestedOutput(1)%Node = 1
   call set_observation(observations(charge_output), requestedOutput, 'wireChargeRange', domain, '')

   requestedOutput(1) = create_observable(2, 2, 2, 3, 3, 3, iBloqueMx)
   call set_observation(observations(bulk_output), requestedOutput, 'bulkRange', domain, '')

   lineSegments(1) = direction_t(3, 3, 3, iEx)
   requestedOutput(1) = create_observable(3, 3, 3, 3, 3, 3, lineIntegral, lineSegments)
   call set_observation(observations(line_output), requestedOutput, 'lineRange', domain, '')

   do i = 1, size(observations)
      call sgg_add_observation(sgg, observations(i))
   end do

   wires => GetHwires()
   if (associated(wires%CurrentSegment)) deallocate (wires%CurrentSegment)
   if (associated(wires%ChargeNode)) deallocate (wires%ChargeNode)
   wires%NumDifferentWires = 0
   wires%NumCurrentSegments = 1
   wires%NumChargeNodes = 2
   allocate (wires%CurrentSegment(1), wires%ChargeNode(2))
   wires%CurrentSegment(1)%OrigIndex = 1
   wires%CurrentSegment(1)%i = 3
   wires%CurrentSegment(1)%j = 3
   wires%CurrentSegment(1)%k = 3
   wires%CurrentSegment(1)%tipofield = iEx
   wires%CurrentSegment(1)%indexmed = 0
   wires%CurrentSegment(1)%orientadoalreves = .false.
   wires%CurrentSegment(1)%delta = 1.0_RKIND_wires
   wires%CurrentSegment(1)%Lind = 1.0_RKIND_wires
   wires%CurrentSegment(1)%ChargePlus => wires%ChargeNode(1)
   wires%CurrentSegment(1)%ChargeMinus => wires%ChargeNode(2)
   wireField = 0.0_RKIND
   wires%CurrentSegment(1)%Efield_wire2main => wireField
   wires%ChargeNode%ChargePresent = 0.0_RKIND_wires
   wires%ChargeNode%ChargePast = 0.0_RKIND_wires

   control = create_control_flags(mpidir=3, finaltimestep=nSteps - 2, nEntradaRoot=nEntrada, &
                                  wiresflavor='holland')
   hasWires = .true.
   call init_outputs(sgg, media, sinpml, tagNumbers, bounds, control, outputRequested, hasWires)
   outputs => GetOutputs()

   call create_dummy_fields(dummyFields, 1, 5, 0.01_RKIND)
   fields%E%x => dummyFields%Ex
   fields%E%y => dummyFields%Ey
   fields%E%z => dummyFields%Ez
   fields%E%deltaX => dummyFields%dxe
   fields%E%deltaY => dummyFields%dye
   fields%E%deltaZ => dummyFields%dze
   fields%H%x => dummyFields%Hx
   fields%H%y => dummyFields%Hy
   fields%H%z => dummyFields%Hz
   fields%H%deltaX => dummyFields%dxh
   fields%H%deltaY => dummyFields%dyh
   fields%H%deltaZ => dummyFields%dzh

   do i = 0, 8
      wires%CurrentSegment(1)%CurrentPast = real(i, RKIND_wires)
      wires%ChargeNode(2)%ChargePresent = 10.0_RKIND_wires + real(i, RKIND_wires)
      dummyFields%Ex(3, 3, 3) = real(i, RKIND)
      call update_outputs(control, sgg%tiempo(i + 1), int(i, SINGLE), fields)
   end do

   test_err = test_err + assert_integer_equal(size(outputs), 4, 'Unexpected time output count')
   test_err = test_err + assert_integer_equal(outputs(current_output)%outputID, WIRE_CURRENT_PROBE_ID, &
                                               'Unexpected wire-current output id')
   test_err = test_err + assert_integer_equal(outputs(charge_output)%outputID, WIRE_CHARGE_PROBE_ID, &
                                               'Unexpected wire-charge output id')
   test_err = test_err + assert_integer_equal(outputs(bulk_output)%outputID, BULK_PROBE_ID, &
                                               'Unexpected bulk output id')
   test_err = test_err + assert_integer_equal(outputs(line_output)%outputID, LINE_PROBE_ID, &
                                               'Unexpected line output id')

   test_err = test_err + assert_integer_equal(outputs(current_output)%wireCurrentProbe%nTime, 2, &
                                               'Wire-current probe ignored its time range')
   test_err = test_err + assert_integer_equal(outputs(charge_output)%wireChargeProbe%nTime, 2, &
                                               'Wire-charge probe ignored its time range')
   test_err = test_err + assert_integer_equal(outputs(bulk_output)%bulkCurrentProbe%nTime, 2, &
                                               'Bulk probe ignored its time range')
   test_err = test_err + assert_integer_equal(outputs(line_output)%lineProbe%nTime, 2, &
                                               'Line probe ignored its time range')

   test_err = test_err + assert_real_equal(outputs(current_output)%wireCurrentProbe%timeStep(1), &
                                           0.4_RKIND_tiempo, 1e-5_RKIND_tiempo, &
                                           'Wire-current probe stored the wrong first time')
   test_err = test_err + assert_real_equal(outputs(current_output)%wireCurrentProbe%timeStep(2), &
                                           0.6_RKIND_tiempo, 1e-5_RKIND_tiempo, &
                                           'Wire-current probe stored the wrong second time')
   test_err = test_err + assert_real_equal(outputs(charge_output)%wireChargeProbe%timeStep(1), &
                                           0.4_RKIND_tiempo, 1e-5_RKIND_tiempo, &
                                           'Wire-charge probe stored the wrong first time')
   test_err = test_err + assert_real_equal(outputs(charge_output)%wireChargeProbe%timeStep(2), &
                                           0.6_RKIND_tiempo, 1e-5_RKIND_tiempo, &
                                           'Wire-charge probe stored the wrong second time')
   test_err = test_err + assert_real_equal(outputs(bulk_output)%bulkCurrentProbe%timeStep(1), &
                                           0.4_RKIND_tiempo, 1e-5_RKIND_tiempo, &
                                           'Bulk probe stored the wrong first time')
   test_err = test_err + assert_real_equal(outputs(bulk_output)%bulkCurrentProbe%timeStep(2), &
                                           0.6_RKIND_tiempo, 1e-5_RKIND_tiempo, &
                                           'Bulk probe stored the wrong second time')
   test_err = test_err + assert_real_equal(outputs(line_output)%lineProbe%timeStep(1), &
                                           0.4_RKIND_tiempo, 1e-5_RKIND_tiempo, &
                                           'Line probe stored the wrong first time')
   test_err = test_err + assert_real_equal(outputs(line_output)%lineProbe%timeStep(2), &
                                           0.6_RKIND_tiempo, 1e-5_RKIND_tiempo, &
                                           'Line probe stored the wrong second time')

   test_err = test_err + assert_real_equal(outputs(current_output)%wireCurrentProbe%currentValues(1)%current, &
                                           4.0_RKIND, 1e-5_RKIND, 'Wire-current probe stored an out-of-range value')
   test_err = test_err + assert_real_equal(outputs(current_output)%wireCurrentProbe%currentValues(2)%current, &
                                           6.0_RKIND, 1e-5_RKIND, 'Wire-current probe applied the wrong cadence')
   test_err = test_err + assert_real_equal(outputs(charge_output)%wireChargeProbe%chargeValue(1), &
                                           14.0_RKIND, 1e-5_RKIND, 'Wire-charge probe stored an out-of-range value')
   test_err = test_err + assert_real_equal(outputs(charge_output)%wireChargeProbe%chargeValue(2), &
                                           16.0_RKIND, 1e-5_RKIND, 'Wire-charge probe applied the wrong cadence')
   test_err = test_err + assert_real_equal(outputs(line_output)%lineProbe%valueForTime(1), &
                                           0.04_RKIND, 1e-5_RKIND, 'Line probe stored an out-of-range value')
   test_err = test_err + assert_real_equal(outputs(line_output)%lineProbe%valueForTime(2), &
                                           0.06_RKIND, 1e-5_RKIND, 'Line probe applied the wrong cadence')

   call delete_outputs(0)
   call remove_folder(testPath, ios)
   nullify (wires%CurrentSegment(1)%Efield_wire2main)
   deallocate (wires%CurrentSegment, wires%ChargeNode)
   wires%NumCurrentSegments = 0
   wires%NumChargeNodes = 0

   err = test_err
end function test_update_time_probe_ranges

integer function test_flush_point_probe() bind(c) result(err)
   ! Verifies point-probe time and frequency data are written and reset on flush.
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
   character(len=BUFSIZE) :: testPath, nEntrada

   type(point_probe_output_t) :: probe
   type(domain_t)             :: domain
   type(cell_coordinate_t)    :: coordinates

   integer :: file_size, n, i
   integer :: test_err = 0
   integer :: ios

   ! Setup
   sep = get_path_separator()
   testPath = join_path(get_temp_folder(), test_folder)
   nEntrada = testPath//sep//test_name

   domain = domain_t( &
            0.0_RKIND_tiempo, 10.0_RKIND_tiempo, 0.1_RKIND_tiempo, &
            10.0_RKIND, 100.0_RKIND, 10, .false.)

   coordinates%x = 2
   coordinates%y = 2
   coordinates%z = 2

   call init_point_probe_output(probe, coordinates, iEx, domain, nEntrada, 3, 0.1_RKIND_tiempo)

   test_err = test_err + assert_integer_equal(size(probe%artifacts), 2, &
                                               'Point probe did not declare both text domains')

   ! A declared scalar artifact remains discoverable when no sample is recorded.
   call flush_point_probe_output(probe)
   inquire (file=probe%filePathTime, size=file_size)
   test_err = test_err + assert_true(file_size > 0, 'Zero-sample point artifact is missing its header')
   test_err = test_err + assert_true(.not. file_exists(trim(probe%path)//'_tm.bin'), &
                                     'Zero-sample point probe created a time binary sidecar')
   test_err = test_err + assert_true(.not. file_exists(trim(probe%path)//'_fq.bin'), &
                                     'Zero-sample point probe created a frequency binary sidecar')

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
   test_err = test_err + assert_true(.not. file_exists(trim(probe%path)//'_tm.bin'), &
                                     'Point probe created a time binary sidecar')
   test_err = test_err + assert_true(.not. file_exists(trim(probe%path)//'_fq.bin'), &
                                     'Point probe created a frequency binary sidecar')

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
   call remove_folder(testPath, ios)

   err = test_err
end function

integer function test_flush_wire_probe_dat() bind(c) result(err)
   use FDETYPES_m, only: RKIND, RKIND_tiempo
   use outputTypes_m, only: wire_current_probe_output_t, wire_charge_probe_output_t, OUTPUT_ARTIFACT_TEXT
   use wireProbeOutput_m, only: flush_wire_current_probe_output, flush_wire_charge_probe_output
   use assertionTools_m, only: assert_true
   use directoryUtils_m, only: create_file_with_path, delete_file, file_exists
   use testOutputUtils_m, only: get_temp_folder
   use directoryUtils_m, only: join_path
   implicit none

   type(wire_current_probe_output_t) :: current_probe
   type(wire_charge_probe_output_t) :: charge_probe
   integer :: current_size, charge_size, ios
   character(len=4096) :: folder

   err = 0
   folder = join_path(get_temp_folder(), 'testing wire')
   current_probe%filePathTime = join_path(folder, 'current_tm.dat')
   current_probe%artifacts(1)%kind = OUTPUT_ARTIFACT_TEXT
   current_probe%artifacts(1)%relative_path = current_probe%filePathTime
   allocate (current_probe%timeStep(1))
   current_probe%nTime = 1
   current_probe%timeStep(1) = 1.0_RKIND_tiempo
   current_probe%currentValues(1)%current = 2.0_RKIND
   call create_file_with_path(current_probe%filePathTime, ios)
   call flush_wire_current_probe_output(current_probe)
   inquire (file=current_probe%filePathTime, size=current_size, iostat=ios)
   err = err + assert_true(ios == 0 .and. current_size > 0, 'Wire-current text record was not published')
   err = err + assert_true(.not. file_exists(join_path(folder, 'current_tm.bin')), &
                           'Wire-current probe created a binary sidecar')

   charge_probe%filePathTime = join_path(folder, 'charge_tm.dat')
   charge_probe%artifacts(1)%kind = OUTPUT_ARTIFACT_TEXT
   charge_probe%artifacts(1)%relative_path = charge_probe%filePathTime
   allocate (charge_probe%timeStep(1), charge_probe%chargeValue(1))
   charge_probe%nTime = 1
   charge_probe%timeStep(1) = 1.0_RKIND_tiempo
   charge_probe%chargeValue(1) = 2.0_RKIND
   call create_file_with_path(charge_probe%filePathTime, ios)
   call flush_wire_charge_probe_output(charge_probe)
   inquire (file=charge_probe%filePathTime, size=charge_size, iostat=ios)
   err = err + assert_true(ios == 0 .and. charge_size > 0, 'Wire-charge text record was not published')
   err = err + assert_true(.not. file_exists(join_path(folder, 'charge_tm.bin')), &
                           'Wire-charge probe created a binary sidecar')

   call delete_file(current_probe%filePathTime, ios)
   call delete_file(charge_probe%filePathTime, ios)
end function test_flush_wire_probe_dat

integer function test_flush_bulk_probe_dat() bind(c) result(err)
   use FDETYPES_m, only: RKIND, RKIND_tiempo
   use outputTypes_m, only: bulk_current_probe_output_t, OUTPUT_ARTIFACT_TEXT
   use bulkProbeOutput_m, only: flush_bulk_probe_output
   use assertionTools_m, only: assert_true
   use directoryUtils_m, only: create_file_with_path, delete_file, file_exists
   use testOutputUtils_m, only: get_temp_folder
   use directoryUtils_m, only: join_path
   implicit none

   type(bulk_current_probe_output_t) :: probe
   integer :: text_size, ios
   character(len=4096) :: folder

   err = 0
   folder = join_path(get_temp_folder(), 'testing bulk')
   probe%filePathTime = join_path(folder, 'probe_tm.dat')
   probe%artifacts(1)%kind = OUTPUT_ARTIFACT_TEXT
   probe%artifacts(1)%relative_path = probe%filePathTime
   allocate (probe%timeStep(1), probe%valueForTime(1))
   probe%nTime = 1
   probe%timeStep(1) = 1.0_RKIND_tiempo
   probe%valueForTime(1) = 2.0_RKIND
   call create_file_with_path(probe%filePathTime, ios)
   call flush_bulk_probe_output(probe)
   inquire (file=probe%filePathTime, size=text_size, iostat=ios)
   err = err + assert_true(ios == 0 .and. text_size > 0, 'Bulk text record was not published')
   err = err + assert_true(.not. file_exists(join_path(folder, 'probe_tm.bin')), &
                           'Bulk probe created a binary sidecar')
   call delete_file(probe%filePathTime, ios)
end function test_flush_bulk_probe_dat

integer function test_multiple_flush_point_probe() bind(c) result(err)
   ! Verifies consecutive point-probe flushes preserve time and frequency data.
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
   character(len=BUFSIZE) :: testPath, nEntrada

   type(point_probe_output_t) :: probe
   type(domain_t)             :: domain
   type(cell_coordinate_t)    :: coordinates

   real(kind=RKIND), allocatable :: expectedTime(:, :)
   real(kind=RKIND), allocatable :: expectedFreq(:, :)

   integer :: n, i, unit
   integer :: test_err = 0
   integer :: ios
   character(len=BUFSIZE) :: header

   ! Setup
   sep = get_path_separator()
   testPath = join_path(get_temp_folder(), test_folder)
   nEntrada = testPath//sep//test_name

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
   unit = 1
   open (unit=unit, file=probe%filePathTime, status='old', action='read')
   read (unit, '(A)', iostat=ios) header
   test_err = test_err + assert_true(ios == 0 .and. trim(header) == 't field', 'Point time header is incorrect')
   test_err = test_err + assert_file_content(unit, expectedTime, 2*n, 2, 1e-06_RKIND)
   close (unit)
   unit = 2
   open (unit=unit, file=probe%filePathFreq, status='old', action='read')
   read (unit, '(A)', iostat=ios) header
   test_err = test_err + assert_true(ios == 0 .and. trim(header) == 'frequency real imaginary', &
                                     'Point frequency header is incorrect')
   test_err = test_err + assert_file_content(unit, expectedFreq, n, 2, 1e-06_RKIND)
   close (unit)

   !Cleanup
   call remove_folder(testPath, ios)

   err = test_err
end function

integer function test_init_movie_probe() bind(c) result(err)
   ! Verifies movie-probe allocation, measurement sizes, and output directory.
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

   real(kind=RKIND), dimension(:), pointer   :: x_steps, y_steps, z_steps

   integer(kind=SINGLE)             :: expectedNumMeasurments
   integer(kind=SINGLE)             :: mpidir = 3
   logical                          :: ThereAreWires = .false.
   logical                          :: outputRequested
   integer(kind=SINGLE)             :: iter
   integer(kind=SINGLE)             :: test_err = 0

   character(len=14), parameter :: test_folder = 'testing_folder'
   character(len=9), parameter :: test_name = 'initMovie'

   character(len=BUFSIZE) :: testPath, nEntrada
   character(len=BUFSIZE) :: expectedProbePath
   integer :: ios

   testPath = join_path(get_temp_folder(), test_folder)
   nEntrada = join_path(testPath, test_name)

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
   tagNumbers = create_tag_list(allocationRange)

   allocate (x_steps(6), source=1.0_RKIND)
   allocate (y_steps(6), source=1.0_RKIND)
   allocate (z_steps(6), source=1.0_RKIND)
   call sgg_set_LineX(dummysgg, x_steps)
   call sgg_set_LineY(dummysgg, y_steps)
   call sgg_set_LineZ(dummysgg, z_steps)

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

   if (size(outputs) == 0 .or. outputs(1)%outputID /= MOVIE_PROBE_ID) then
      err = assert_true(.false., 'Movie probe initialisation did not create a movie output')
      call remove_folder(testPath, ios)
      return
   end if

   test_err = test_err + assert_integer_equal(outputs(1)%outputID, MOVIE_PROBE_ID, 'Unexpected probe id')
   test_err = test_err + assert_integer_equal(outputs(1)%movieProbe%nPoints, expectedNumMeasurments, 'Unexpected number of measurements')
   test_err = test_err + assert_integer_equal(size(outputs(1)%movieProbe%xValueForTime), expectedNumMeasurments*OUTPUT_TIME_BUFFER_SIZE, 'Unexpected allocation size')
   test_err = test_err + assert_integer_equal(size(outputs(1)%movieProbe%timeStep), OUTPUT_TIME_BUFFER_SIZE, 'Unexpected timestep buffer size')

   expectedProbePath = trim(nEntrada)//wordSeparation//'movieProbe_BC_2_2_2__5_5_5'
   test_err = test_err + assert_string_equal(outputs(1)%movieProbe%path, expectedProbePath, 'Unexpected path')
   test_err = test_err + assert_true(folder_exists(expectedProbePath), 'Movie folder do not exist')
   test_err = test_err + assert_true(file_exists(trim(outputs(1)%movieProbe%filesPath)//'.bin'), &
                                     'Movie binary payload was not created')
   test_err = test_err + assert_true(.not. file_exists(trim(outputs(1)%movieProbe%filesPath)//'.json'), &
                                     'Movie JSON descriptor was created')
   test_err = test_err + assert_true(file_exists(trim(outputs(1)%movieProbe%filesPath)//'_geometry.xdmf'), &
                                     'Movie geometry XDMF was not created')
   test_err = test_err + assert_true(file_exists(trim(outputs(1)%movieProbe%filesPath)//'_geometry.h5'), &
                                     'Movie geometry HDF5 was not created')
   !Cleanup
   call remove_folder(testPath, ios)

   err = test_err
end function

integer function test_update_movie_probe() bind(c) result(err)
   ! Verifies movie probes honour an explicit time window and sampling period.
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

   real(kind=RKIND), dimension(:), pointer   :: x_steps, y_steps, z_steps

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

   character(len=BUFSIZE) :: testPath, nEntrada
   integer :: ios

   testPath = join_path(get_temp_folder(), test_folder)
   nEntrada = join_path(testPath, test_name)

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
   tagNumbers = create_tag_list(allocationRange)

   allocate (x_steps(6), source=1.0_RKIND)
   allocate (y_steps(6), source=1.0_RKIND)
   allocate (z_steps(6), source=1.0_RKIND)
   call sgg_set_LineX(dummysgg, x_steps)
   call sgg_set_LineY(dummysgg, y_steps)
   call sgg_set_LineZ(dummysgg, z_steps)

   movieObservable = create_movie_observation(2, 2, 2, 5, 5, 5, iCur)
   movieObservable%InitialTime = 0.15_RKIND
   movieObservable%FinalTime = 0.35_RKIND
   movieObservable%TimeStep = 0.2_RKIND
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

   dummyControl = create_control_flags(nEntradaRoot=nEntrada, mpidir=mpidir, &
                                       finaltimestep=nTimeSteps - 2)

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

   do iter = 1, 4
      call update_outputs(dummyControl, dummysgg%tiempo(iter + 1), iter, fields)
   end do

   test_err = test_err + assert_real_equal(outputs(1)%movieProbe%yValueForTime(1, 1), &
                                           -0.2_RKIND, 1e-5_RKIND, 'Value error')

   test_err = test_err + assert_real_equal(outputs(1)%movieProbe%yValueForTime(1, 2), &
                                           0.4_RKIND, 1e-5_RKIND, 'Value error')

   test_err = test_err + assert_real_equal(outputs(1)%movieProbe%yValueForTime(1, 3), &
                                           0.0_RKIND, 1e-5_RKIND, 'Value error')

   test_err = test_err + assert_real_equal(outputs(1)%movieProbe%yValueForTime(1, 4), &
                                           0.0_RKIND, 1e-5_RKIND, 'Value error')

   test_err = test_err + assert_integer_equal( &
              size(outputs(1)%movieProbe%timeStep), OUTPUT_TIME_BUFFER_SIZE, 'Unexpected timestep buffer size')
   test_err = test_err + assert_integer_equal(outputs(1)%movieProbe%nTime, 1, 'Movie update did not buffer a timestep')
   test_err = test_err + assert_true(outputs(1)%movieProbe%timeStep(1) == dummysgg%tiempo(3), &
                                      'Movie update stored an incorrect timestep')

   !Cleanup
   call remove_folder(testPath, ios)

   err = test_err
end function

integer function test_flush_movie_probe() bind(c) result(err)
   ! Verifies movie flushes binary samples and clears its in-memory buffers.
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

   real(kind=RKIND), dimension(:), pointer   :: x_steps, y_steps, z_steps

   integer(kind=SINGLE)             :: expectedNumMeasurments
   integer(kind=SINGLE)             :: mpidir = 3
   integer(kind=SINGLE)             :: iter
   integer(kind=SINGLE)             :: test_err = 0
   logical                          :: ThereAreWires = .false.
   logical                          :: outputRequested

   character(len=14), parameter :: test_folder = 'testing_folder'
   character(len=10), parameter :: test_name = 'flushMovie'
   character(len=BUFSIZE) :: testPath
   character(len=BUFSIZE) :: nEntrada
   character(len=BUFSIZE) :: expectedPath
   integer :: binaryBytes, ios

   testPath = join_path(get_temp_folder(), test_folder)
   nEntrada = join_path(testPath, test_name)

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
   tagNumbers = create_tag_list(allocationRange)

   allocate (x_steps(6), source=1.0_RKIND)
   allocate (y_steps(6), source=1.0_RKIND)
   allocate (z_steps(6), source=1.0_RKIND)
   call sgg_set_LineX(dummysgg, x_steps)
   call sgg_set_LineY(dummysgg, y_steps)
   call sgg_set_LineZ(dummysgg, z_steps)

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
   outputs(2)%movieProbe%xValueForTime(1, :4) = [0.1_RKIND, 0.2_RKIND, 0.3_RKIND, 0.4_RKIND]

   !movieMagneticYObservable
   outputs(3)%movieProbe%nTime = 1
   outputs(3)%movieProbe%timeStep(1) = 0.5_RKIND_tiempo
   outputs(3)%movieProbe%yValueForTime(1, :4) = [0.1_RKIND, 0.2_RKIND, 0.3_RKIND, 0.4_RKIND]

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
   outputs(2)%movieProbe%xValueForTime(2, :4) = [1.1_RKIND, 1.2_RKIND, 1.3_RKIND, 1.4_RKIND]

   !movieMagneticYObservable
   outputs(3)%movieProbe%nTime = 2
   outputs(3)%movieProbe%timeStep(2) = 0.75_RKIND_tiempo
   outputs(3)%movieProbe%yValueForTime(2, :4) = [1.1_RKIND, 1.2_RKIND, 1.3_RKIND, 1.4_RKIND]

   call flush_outputs(dummysgg%tiempo, 1_SINGLE, dummyControl, fields, dummyBound, .false.)

   expectedPath = trim(outputs(1)%movieProbe%filesPath)
   test_err = test_err + assert_true(file_exists(trim(expectedPath)//'.bin'), 'Movie binary payload does not exist')
   inquire (file=trim(expectedPath)//'.bin', size=binaryBytes)
   test_err = test_err + assert_integer_equal(binaryBytes, 8*56, 'Movie binary record layout changed')
   test_err = test_err + assert_integer_equal(outputs(1)%movieProbe%nTime, 0, 'Movie flush did not clear current samples')
   test_err = test_err + assert_integer_equal(outputs(1)%movieProbe%nTimesFlushed, 2, 'Movie flush did not count samples')
   test_err = test_err + assert_integer_equal(outputs(2)%movieProbe%nTime, 0, 'Electric movie flush did not clear samples')
   test_err = test_err + assert_integer_equal(outputs(3)%movieProbe%nTime, 0, 'Magnetic movie flush did not clear samples')

   call close_outputs()
   call cleanup_test_artifacts(testPath, ios)

   err = test_err
end function

integer function test_close_movie_probe() bind(c) result(err)
   ! Verifies movie close finalises visualisation artifacts without a JSON sidecar.
   use output_m
   use outputTypes_m
   use testOutputUtils_m
   use FDETYPES_TOOLS
   use sggMethods_m
   use assertionTools_m
   use directoryUtils_m
   implicit none

   type(SGGFDTDINFO_t) :: dummysgg
   type(sim_control_t) :: dummyControl
   type(bounds_t) :: dummyBound
   type(solver_output_t), pointer :: outputs(:)
   type(media_matrices_t) :: media
   type(MediaData_t), allocatable, target :: simulationMaterials(:)
   type(MediaData_t), pointer :: simulationMaterialsPtr(:)
   type(limit_t) :: sinpml(6)
   type(taglist_t) :: tagNumbers
   type(XYZlimit_t) :: sweep(6)
   type(Obses_t) :: movieObservable
   real(kind=RKIND_tiempo), pointer :: timeArray(:)
   real(kind=RKIND), pointer :: x_steps(:), y_steps(:), z_steps(:)
   real(kind=RKIND_tiempo) :: dt = 0.1_RKIND_tiempo
   integer(kind=SINGLE) :: iter, mpidir = 3_SINGLE, test_err = 0_SINGLE
   logical :: outputRequested, thereAreWires = .false.
   character(len=14), parameter :: test_folder = 'testing_folder'
   character(len=10), parameter :: test_name = 'closeMovie'
   character(len=BUFSIZE) :: testPath, inputPath, expectedPath
   integer :: ios

   err = 1
   testPath = join_path(get_temp_folder(), test_folder)
   inputPath = join_path(testPath, test_name)
   call sgg_init(dummysgg)
   call init_time_array(timeArray, 100_SINGLE, dt)
   call sgg_set_tiempo(dummysgg, timeArray)
   call sgg_set_dt(dummysgg, dt)
   call init_simulation_material_list(simulationMaterials)
   simulationMaterialsPtr => simulationMaterials
   call sgg_set_NumMedia(dummysgg, size(simulationMaterials))
   call sgg_set_Med(dummysgg, simulationMaterialsPtr)
   sweep = create_xyz_limit_array(0, 0, 0, 6, 6, 6)
   call sgg_set_Sweep(dummysgg, sweep)
   call sgg_set_SINPMLSweep(dummysgg, sweep)
   call sgg_set_NumPlaneWaves(dummysgg, 1)
   call sgg_set_Alloc(dummysgg, sweep)
   tagNumbers = create_tag_list(sweep)
   allocate (x_steps(6), source=1.0_RKIND)
   allocate (y_steps(6), source=1.0_RKIND)
   allocate (z_steps(6), source=1.0_RKIND)
   call sgg_set_LineX(dummysgg, x_steps)
   call sgg_set_LineY(dummysgg, y_steps)
   call sgg_set_LineZ(dummysgg, z_steps)
   movieObservable = create_movie_observation(2, 2, 2, 5, 5, 5, iCur)
   call sgg_add_observation(dummysgg, movieObservable)
   call create_geometry_media(media, 0, 8, 0, 8, 0, 8)
   do iter = 1, 6
      sinpml(iter) = create_limit_t(0, 8, 0, 8, 0, 8, 10, 10, 10)
   end do
   dummyControl = create_control_flags(nEntradaRoot=inputPath, mpidir=mpidir)
   call init_outputs(dummysgg, media, sinpml, tagNumbers, dummyBound, dummyControl, outputRequested, thereAreWires)
   outputs => GetOutputs()

   call close_outputs()
   expectedPath = trim(outputs(1)%movieProbe%filesPath)
   test_err = test_err + assert_true(file_exists(trim(expectedPath)//'.xdmf'), 'Movie XDMF metadata was not finalised')
   test_err = test_err + assert_true(file_exists(trim(expectedPath)//'.h5'), 'Movie HDF5 payload was not finalised')
   test_err = test_err + assert_true(.not. file_exists(trim(expectedPath)//'.json'), &
                                     'Movie JSON descriptor was published')

   call cleanup_test_artifacts(testPath, ios)
   err = test_err
end function test_close_movie_probe

integer function test_init_frequency_slice_probe() bind(c) result(err)
   ! Verifies frequency-slice allocation, dimensions, paths, and visualisation output.
   use output_m
   use outputTypes_m
   use frequencySliceProbeOutput_m, only: flush_frequency_slice_probe_output, close_frequency_slice_probe_output
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
   real(kind=RKIND), pointer        :: x_steps(:), y_steps(:), z_steps(:)
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

   character(len=BUFSIZE) :: testPath, nEntrada
   character(len=BUFSIZE) :: expectedProbePath
   integer :: ios

   testPath = join_path(get_temp_folder(), test_folder)
   nEntrada = join_path(testPath, test_name)
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
   tagNumbers = create_tag_list(allocationRange)

   allocate (x_steps(6), source=1.0_RKIND)
   allocate (y_steps(6), source=1.0_RKIND)
   allocate (z_steps(6), source=1.0_RKIND)
   call sgg_set_LineX(dummysgg, x_steps)
   call sgg_set_LineY(dummysgg, y_steps)
   call sgg_set_LineZ(dummysgg, z_steps)

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
   test_err = test_err + assert_string_equal(outputs(1)%frequencySliceProbe%path, expectedProbePath, 'Unexpected path')
   test_err = test_err + assert_true(folder_exists(expectedProbePath), 'Frequency Slice folder do not exist')

   outputs(1)%frequencySliceProbe%xValueForFreq = (0.0_CKIND, 0.0_CKIND)
   outputs(1)%frequencySliceProbe%yValueForFreq = (0.0_CKIND, 0.0_CKIND)
   outputs(1)%frequencySliceProbe%zValueForFreq = (0.0_CKIND, 0.0_CKIND)
   outputs(1)%frequencySliceProbe%xValueForFreq(1, 1) = (1.0_CKIND, 2.0_CKIND)
   outputs(1)%frequencySliceProbe%yValueForFreq(1, 1) = (3.0_CKIND, 4.0_CKIND)
   outputs(1)%frequencySliceProbe%zValueForFreq(1, 1) = (5.0_CKIND, 6.0_CKIND)
   call flush_frequency_slice_probe_output(outputs(1)%frequencySliceProbe)
   call close_frequency_slice_probe_output(outputs(1)%frequencySliceProbe)

   expectedProbePath = trim(outputs(1)%frequencySliceProbe%filesPath)
   test_err = test_err + assert_true(.not. file_exists(trim(expectedProbePath)//'.bin'), &
                                     'Frequency binary payload was published')
   test_err = test_err + assert_true(file_exists(trim(expectedProbePath)//'.xdmf'), 'Frequency XDMF metadata does not exist')
   test_err = test_err + assert_true(file_exists(trim(expectedProbePath)//'.h5'), 'Frequency HDF5 payload does not exist')
   test_err = test_err + assert_true(file_exists(trim(expectedProbePath)//'_geometry.xdmf'), &
                                     'Frequency geometry XDMF metadata does not exist')
   test_err = test_err + assert_true(file_exists(trim(expectedProbePath)//'_geometry.h5'), &
                                     'Frequency geometry HDF5 payload does not exist')
   test_err = test_err + assert_true(.not. file_exists(trim(expectedProbePath)//'.json'), &
                                     'Frequency JSON descriptor was published')
   call remove_folder(testPath, ios)

   err = test_err
end function

integer function test_update_frequency_slice_probe() bind(c) result(err)
   ! Verifies frequency-slice gradients accumulate on every simulation step.
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
   real(kind=RKIND), pointer        :: x_steps(:), y_steps(:), z_steps(:)
   real(kind=RKIND_tiempo)          :: dt = 0.1_RKIND_tiempo
   integer(kind=SINGLE)             :: nTimeSteps = 100_SINGLE

   type(dummyFields_t), target     :: dummyFields
   type(fields_reference_t)        :: fields
   complex(kind=CKIND), allocatable :: firstFrequencyUpdate(:)

   integer(kind=SINGLE)             :: expectedNumMeasurments
   integer(kind=SINGLE)             :: expectedNumberFrequencies
   integer(kind=SINGLE)             :: mpidir = 3
   integer(kind=SINGLE)             :: iter
   integer(kind=SINGLE)             :: test_err = 0
   logical                          :: ThereAreWires = .false.
   logical                          :: outputRequested

   character(len=14), parameter :: test_folder = 'testing_folder'
   character(len=13), parameter :: test_name = 'initFrequency'

   character(len=BUFSIZE) :: testPath, nEntrada
   integer :: ios

   testPath = join_path(get_temp_folder(), test_folder)
   nEntrada = join_path(testPath, test_name)

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
   tagNumbers = create_tag_list(allocationRange)

   allocate (x_steps(6), source=1.0_RKIND)
   allocate (y_steps(6), source=1.0_RKIND)
   allocate (z_steps(6), source=1.0_RKIND)
   call sgg_set_LineX(dummysgg, x_steps)
   call sgg_set_LineY(dummysgg, y_steps)
   call sgg_set_LineZ(dummysgg, z_steps)

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

   call update_outputs(dummyControl, dummysgg%tiempo(3), 2_SINGLE, fields)
   firstFrequencyUpdate = outputs(1)%frequencySliceProbe%yValueForFreq(1, :)
   do iter = 3, 5
      call update_outputs(dummyControl, dummysgg%tiempo(iter + 1), iter, fields)
   end do

   test_err = test_err + assert_integer_equal(outputs(1)%outputID, &
                                              FREQUENCY_SLICE_PROBE_ID, 'Unexpected probe id')

   test_err = test_err + assert_integer_equal(outputs(1)%frequencySliceProbe%nPoints, &
                                              expectedNumMeasurments, 'Unexpected number of measurements')

   test_err = test_err + assert_integer_equal( &
               size(outputs(1)%frequencySliceProbe%frequencySlice), &
               expectedNumberFrequencies, 'Unexpected allocation size')
   test_err = test_err + assert_real_equal(outputs(1)%frequencySliceProbe%frequencySlice(1), &
                                           0.0_RKIND, 1e-5_RKIND, 'Unexpected initial frequency')
   test_err = test_err + assert_real_equal(outputs(1)%frequencySliceProbe%frequencySlice(6), &
                                           100.0_RKIND, 1e-5_RKIND, 'Unexpected final frequency')

   !This test generates X Gradient for H. It is expected to detect none Current accros X axis and Opposite values for Y and Z

   test_err = test_err + assert_array_value(outputs(1)%frequencySliceProbe%xValueForFreq, (0.0_CKIND , 0.0_CKIND), errormessage='Detected Current on X Axis for Hx gradient')
   test_err = test_err + assert_arrays_equal(outputs(1)%frequencySliceProbe%yValueForFreq, &
                                 -1.0_RKIND*outputs(1)%frequencySliceProbe%zValueForFreq, errormessage='Unequal values for Y and -Z')
   test_err = test_err + assert_true(any(abs(firstFrequencyUpdate) > 1e-6_RKIND), &
                                         'First frequency update produced no measurable value')
   test_err = test_err + assert_arrays_equal(outputs(1)%frequencySliceProbe%yValueForFreq(1, :), &
                                              4.0_RKIND*firstFrequencyUpdate, &
                                              errormessage='Frequency slice was not updated on every step')

   call remove_folder(testPath, ios)

   err = test_err
end function
