message(STATUS "Creating build system for fdtd_tests")

enable_testing ()

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

if (SEMBA_FDTD_ENABLE_MTLN)
    add_subdirectory(mtln)
    add_subdirectory(system)
    set(MTLN_TESTS_LIBRARY mtln_tests)
    set(SYSTEM_TESTS_LIBRARY system_tests)
endif()
if (SEMBA_FDTD_ENABLE_SMBJSON)
    add_subdirectory(smbjson)
    set(SMBJSON_TESTS_LIBRARY smbjson_tests)
endif()

add_executable(fdtd_tests 
    "fdtd_tests.cpp"
)

target_link_libraries(fdtd_tests
    ${MTLN_TESTS_LIBRARY}
    ${SMBJSON_TESTS_LIBRARY}
    ${SYSTEM_TESTS_LIBRARY}
    GTest::gtest_main
)


if (${CMAKE_SYSTEM_NAME} STREQUAL "Windows" AND SEMBA_FDTD_ENABLE_MTLN)
    add_custom_command(TARGET fdtd_tests POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${NGSPICE_DIR}ngspice.dll"      
            $<TARGET_FILE_DIR:fdtd_tests>
    )
endif()