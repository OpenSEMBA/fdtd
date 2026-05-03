message(STATUS "Creating build system for test/semba_fdtd")

add_library(semba_fdtd_test_fortran
    "test_nfde2sgg.F90"
)

target_link_libraries(semba_fdtd_test_fortran
    semba-main
    smbjson
    fhash
    ${JSONFORTRAN_LIB}
)

add_library(semba_fdtd_tests "semba_fdtd_tests.cpp")

target_link_libraries(semba_fdtd_tests
    semba_fdtd_test_fortran
    GTest::gtest
)
