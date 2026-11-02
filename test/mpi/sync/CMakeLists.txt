message(STATUS "Creating build system for test/mpi/sync")

add_executable(sync_timeout_mpi_test "test_sync_timeout.F90")
target_link_libraries(sync_timeout_mpi_test PRIVATE semba-outputs semba-reports MPI::MPI_Fortran)
set_target_properties(sync_timeout_mpi_test PROPERTIES LINKER_LANGUAGE Fortran)

add_test(
    NAME sync_timeout_normal_exchange
    COMMAND ${MPIEXEC_EXECUTABLE}
        ${MPIEXEC_NUMPROC_FLAG} 2
        ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:sync_timeout_mpi_test> ok
        ${MPIEXEC_POSTFLAGS}
)
add_test(
    NAME sync_timeout_aborts_on_stalled_neighbour
    COMMAND ${MPIEXEC_EXECUTABLE}
        ${MPIEXEC_NUMPROC_FLAG} 2
        ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:sync_timeout_mpi_test> stall
        ${MPIEXEC_POSTFLAGS}
)
set_tests_properties(
    sync_timeout_normal_exchange
    PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1" TIMEOUT 120
)
set_tests_properties(
    sync_timeout_aborts_on_stalled_neighbour
    PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1;SEMBA_FDTD_SYNC_TIMEOUT_S=5" TIMEOUT 120
    PASS_REGULAR_EXPRESSION "No answer from a neighbour rank"
)
