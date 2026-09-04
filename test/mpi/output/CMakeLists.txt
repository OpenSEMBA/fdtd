message(STATUS "Creating build system for test/mpi/output")

add_dependencies(fdtd-output vtkAPI)
add_executable(output_mpi_integration_test "test_output_mpi.F90")
target_link_libraries(output_mpi_integration_test PRIVATE fdtd-output vtkAPI fdtd-utils MPI::MPI_Fortran)
set_target_properties(output_mpi_integration_test PROPERTIES LINKER_LANGUAGE Fortran)
add_executable(output_transport_mpi_test "test_output_transport.F90")
target_link_libraries(output_transport_mpi_test PRIVATE fdtd-output MPI::MPI_Fortran)
set_target_properties(output_transport_mpi_test PROPERTIES LINKER_LANGUAGE Fortran)

add_test(
    NAME output_mpi_collective_and_root_aggregation
    COMMAND ${MPIEXEC_EXECUTABLE}
        ${MPIEXEC_NUMPROC_FLAG} ${SEMBA_FDTD_OUTPUT_MPI_TEST_RANKS}
        ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:output_mpi_integration_test>
        ${MPIEXEC_POSTFLAGS}
)
add_test(
    NAME output_transport
    COMMAND ${MPIEXEC_EXECUTABLE}
        ${MPIEXEC_NUMPROC_FLAG} ${SEMBA_FDTD_OUTPUT_MPI_TEST_RANKS}
        ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:output_transport_mpi_test>
        ${MPIEXEC_POSTFLAGS}
)

set_tests_properties(
    output_mpi_collective_and_root_aggregation
    output_transport
    PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1"
)
