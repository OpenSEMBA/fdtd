message(STATUS "Creating build system for test/healer")

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

add_library(healer_test_fortran
    "test_healer.F90"
)

target_link_libraries(healer_test_fortran
    semba-main
    semba-reports
    fhash
)

add_library(healer_tests "healer_tests.cpp")

target_link_libraries(healer_tests
    healer_test_fortran
    GTest::gtest
)
