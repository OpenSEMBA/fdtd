message(STATUS "Creating build system for test/observation")

add_library(
    observation_test_fortran
    "observation_testingTools.F90"
    "test_observation.F90"
    "test_preprocess.F90"
    "test_observation_init.F90"
)

target_link_libraries(observation_test_fortran
    semba-outputs
    test_utils_fortran
)

add_library(observation_tests "observation_tests.cpp")

target_link_libraries(observation_tests
    observation_test_fortran
    GTest::gtest
)