message(STATUS "Creating build system for test/output")

add_library(
    output_test_fortran
    "test_output.F90"
    "test_output_utils.F90"
    "test_volumic_utils.F90"
    "test_output_decomposition.F90"
    "test_output_collective.F90"
)

add_library(vtkAPI_test_fortran
    "test_vtkAPI.F90"
)

add_library(xdmfAPI_test_fortran
    "test_xdmfAPI.F90"
)

add_library(output_tests "output_tests.cpp")
add_library(vtkAPI_tests "vtkAPI_tests.cpp")
add_library(xdmfAPI_tests "xdmfAPI_tests.cpp")

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

target_link_libraries(xdmfAPI_test_fortran
    xdmfAPI
    test_utils_fortran
    ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES}
)

target_link_libraries(xdmfAPI_tests
    xdmfAPI_test_fortran
    GTest::gtest
)

target_include_directories(xdmfAPI_test_fortran PRIVATE ${HDF5_INCLUDE_DIRS})

add_test(
    NAME output_implementation_wiring
    COMMAND ${CMAKE_COMMAND}
        -DPROJECT_SOURCE_DIR=${CMAKE_SOURCE_DIR}
        -P ${CMAKE_CURRENT_SOURCE_DIR}/test_output_implementation_wiring.cmake
)

if(SEMBA_FDTD_ENABLE_MPI)
    add_dependencies(fdtd-output vtkAPI xdmfAPI)
    add_executable(output_mpi_integration_test "test_output_mpi.F90")
    target_link_libraries(output_mpi_integration_test PRIVATE fdtd-output MPI::MPI_Fortran)

    add_test(
        NAME output_mpi_collective_and_root_aggregation
        COMMAND ${MPIEXEC_EXECUTABLE}
            ${MPIEXEC_NUMPROC_FLAG} ${SEMBA_FDTD_OUTPUT_MPI_TEST_RANKS}
            ${MPIEXEC_PREFLAGS}
            $<TARGET_FILE:output_mpi_integration_test>
            ${MPIEXEC_POSTFLAGS}
    )
endif()
