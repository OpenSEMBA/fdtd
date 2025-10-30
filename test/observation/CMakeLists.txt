message(STATUS "Creating build system for test/vtk")

add_library(
    observation_test_fortran
    "test_observation.F90"
)

target_link_libraries(observation_test_fortran
    semba-outputs
)

add_library(observation_tests "observation_tests.cpp")

target_link_libraries(observation_tests
    observation_test_fortran
    GTest::gtest
)