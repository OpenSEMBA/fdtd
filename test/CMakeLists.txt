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
    add_subdirectory(utils)
    add_subdirectory(smbjson)
    set(SMBJSON_TESTS_LIBRARY smbjson_tests)
    add_subdirectory(unit/rotate)
    set(ROTATE_TESTS_LIBRARY rotate_tests)
    add_subdirectory(unit/output)
    set(OUTPUT_TESTS_LIBRARY output_tests)
    set(VTK_API_TESTS_LIBRARY vtkAPI_tests)
endif()

if(SEMBA_FDTD_ENABLE_MPI)
    add_subdirectory(mpi/output)
endif()

add_subdirectory(conformal)
set(CONFORMAL_TESTS_LIBRARY conformal_tests)

add_subdirectory(unit/preprocess)
set(PREPROCESS_TESTS_LIBRARY preprocess_tests)

add_executable(fdtd_tests 
    "fdtd_tests.cpp"
)

target_link_libraries(fdtd_tests
    ${MTLN_TESTS_LIBRARY}
    ${SMBJSON_TESTS_LIBRARY}
    ${CONFORMAL_TESTS_LIBRARY}
    ${PREPROCESS_TESTS_LIBRARY}
    ${ROTATE_TESTS_LIBRARY}
    ${OUTPUT_TESTS_LIBRARY}
    ${VTK_API_TESTS_LIBRARY}
    ${SYSTEM_TESTS_LIBRARY}
    GTest::gtest_main
)

add_test(NAME fdtd_unit COMMAND $<TARGET_FILE:fdtd_tests>)
set_tests_properties(fdtd_unit PROPERTIES WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

if(SEMBA_FDTD_ENABLE_MPI)
    find_package(MPI REQUIRED COMPONENTS CXX)
    target_link_libraries(fdtd_tests MPI::MPI_Fortran MPI::MPI_CXX)
endif()
