message(STATUS "Creating build system for test/maloneySheet")

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

add_library(maloneySheet_test_fortran
    "test_maloneySheet.F90"
)

target_link_libraries(maloneySheet_test_fortran
    semba-components
)

add_library(maloneySheet_tests "maloneySheet_tests.cpp")

target_link_libraries(maloneySheet_tests
    maloneySheet_test_fortran
    GTest::gtest
)
