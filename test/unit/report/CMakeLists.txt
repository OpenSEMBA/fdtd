message(STATUS "Creating build system for test/unit/report")

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

add_library(report_test_fortran
    "test_report_log.F90"
)

target_link_libraries(report_test_fortran
    semba-reports
)

add_library(report_tests "report_tests.cpp")

target_link_libraries(report_tests
    report_test_fortran
    GTest::gtest
)
