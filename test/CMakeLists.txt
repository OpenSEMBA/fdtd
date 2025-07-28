message(STATUS "Creating build system for fdtd_tests")

enable_testing ()
project(fdtd_tests)
enable_language(C CXX Fortran)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

if (SEMBA_FDTD_ENABLE_MTLN)
    add_subdirectory(mtln)
    set(MTLN_TESTS_LIBRARY mtln_tests)
    if (NOT SEMBA_FDTD_ENABLE_MPI)
        add_subdirectory(system)
        set(SYSTEM_TESTS_LIBRARY system_tests)
    endif()
endif()
if (SEMBA_FDTD_ENABLE_SMBJSON)
    add_subdirectory(smbjson)
    set(SMBJSON_TESTS_LIBRARY smbjson_tests)
    add_subdirectory(rotate)
    set(ROTATE_TESTS_LIBRARY rotate_tests)
endif()

add_executable(fdtd_tests 
    "fdtd_tests.cpp"
)

target_link_libraries(fdtd_tests
    ${MTLN_TESTS_LIBRARY}
    ${SMBJSON_TESTS_LIBRARY}
    ${ROTATE_TESTS_LIBRARY}
    ${SYSTEM_TESTS_LIBRARY}
    ${HDF_TESTS_LIBRARY}
    GTest::gtest_main
)