message(STATUS "Creating build system for test/lorentz")

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

add_library(
    lorentz_test_fortran
    "test_lorentz.F90"
)

target_link_libraries(lorentz_test_fortran
    semba-main
)

add_library(lorentz_tests "lorentz_tests.cpp")

target_link_libraries(lorentz_tests
    lorentz_test_fortran
    GTest::gtest
)
