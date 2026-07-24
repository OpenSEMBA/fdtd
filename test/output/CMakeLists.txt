message(STATUS "Creating build system for test/output")

add_library(
    output_test_fortran
    "unit/test_probe_output.F90"
    "unit/test_output_utils.F90"
    "unit/test_volumic_probe_utils.F90"
    "unit/test_output_decomposition.F90"
    "unit/test_output_collective.F90"
    "unit/test_output_publication_contract.F90"
    "unit/test_output_lifecycle.F90"
    "unit/test_output_transport.F90"
    "unit/test_output_metadata.F90"
    "unit/test_output_binary.F90"
)

add_library(vtkAPI_test_fortran
    "unit/test_vtk_api.F90"
)

add_library(output_tests "output_tests.cpp")
add_library(vtkAPI_tests "vtkAPI_tests.cpp")

target_link_libraries(output_test_fortran
    fdtd-output
    test_utils_fortran
)
target_link_libraries(output_tests
    output_test_fortran
    GTest::gtest
)
target_link_libraries(vtkAPI_test_fortran
    vtkAPI
    test_utils_fortran
)
target_link_libraries(vtkAPI_tests
    vtkAPI_test_fortran
    GTest::gtest
)
add_test(
    NAME output_implementation_wiring
    COMMAND ${CMAKE_COMMAND}
        -DPROJECT_SOURCE_DIR=${CMAKE_SOURCE_DIR}
        -P ${CMAKE_CURRENT_SOURCE_DIR}/unit/test_output_implementation_wiring.cmake
)

if(SEMBA_FDTD_ENABLE_MPI)
    add_dependencies(fdtd-output vtkAPI)
    add_executable(output_mpi_integration_test "test_output_mpi.F90")
    target_link_libraries(output_mpi_integration_test PRIVATE fdtd-output vtkAPI fdtd-utils MPI::MPI_Fortran)
    add_executable(output_transport_mpi_test "mpi/test_output_transport.F90")
    target_link_libraries(output_transport_mpi_test PRIVATE fdtd-output MPI::MPI_Fortran)

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
endif()
