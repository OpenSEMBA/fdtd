message(STATUS "Creating build system for test/unit/output")

add_library(
    output_test_fortran
    # Shared fixtures and output utilities.
    "test_output_utils.F90"
    "test_volumic_probe_utils.F90"

    # Distributed publication contracts.
    "test_output_decomposition.F90"
    "test_output_collective.F90"
    "test_output_publication_contract.F90"
    "test_output_lifecycle.F90"
    "test_output_transport.F90"

    # Artifact adapters.
    "test_output_metadata.F90"
    "test_output_binary.F90"

    # Probe lifecycle behaviour.
    "test_probe_output.F90"
)

add_library(vtkAPI_test_fortran
    "test_vtk_api.F90"
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
        -P ${CMAKE_CURRENT_SOURCE_DIR}/test_output_implementation_wiring.cmake
)
