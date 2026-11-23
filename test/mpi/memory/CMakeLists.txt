message(STATUS "Creating build system for test/mpi/memory")

add_executable(memory_budget_mpi_test "test_memory_budget_mpi.F90")
target_link_libraries(memory_budget_mpi_test PRIVATE semba-reports fdtd-utils MPI::MPI_Fortran)
set_target_properties(memory_budget_mpi_test PROPERTIES LINKER_LANGUAGE Fortran)

add_test(
    NAME memory_budget_mpi
    COMMAND ${MPIEXEC_EXECUTABLE}
        ${MPIEXEC_NUMPROC_FLAG} 2
        ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:memory_budget_mpi_test>
        ${MPIEXEC_POSTFLAGS}
)
set_tests_properties(memory_budget_mpi PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1" TIMEOUT 60)
