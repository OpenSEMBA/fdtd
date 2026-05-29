message(STATUS "Creating build system for test/cpp")

add_subdirectory(support)

add_executable(cpp_tests
    cpp_tests.cpp
)

target_include_directories(cpp_tests PRIVATE
    ${CPP_SOURCE_ROOT}/src_cpp/json_parser
    ${CPP_SOURCE_ROOT}/src_cpp/conformal
    ${CPP_SOURCE_ROOT}/src_cpp/main
    ${CPP_SOURCE_ROOT}/src_cpp/mtln
    ${CPP_SOURCE_ROOT}/test/cpp/support
    ${CPP_SOURCE_ROOT}/test/cpp
)

target_link_libraries(cpp_tests
    gtest_main
    nlohmann_json::nlohmann_json
    smbjson
    conformal
    semba-main
    semba-test-support
)

if(SEMBA_FDTD_ENABLE_MTLN)
    target_compile_definitions(cpp_tests PRIVATE CompileWithMTLN)
    target_link_libraries(cpp_tests mtlnsolver ngspice_interface ngspice ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
    target_include_directories(cpp_tests PRIVATE
        ${CPP_SOURCE_ROOT}/external/ngspice/src/include
    )
endif()

if(SEMBA_FDTD_ENABLE_MPI)
    target_compile_definitions(cpp_tests PRIVATE CompileWithMPI)
else()
    target_compile_definitions(cpp_tests PRIVATE CompileWithoutMPI)
endif()

if(SEMBA_FDTD_ENABLE_HDF)
    target_compile_definitions(cpp_tests PRIVATE SEMBA_CPP_ENABLE_HDF5)
endif()

add_test(NAME cpp_tests COMMAND cpp_tests)
set_tests_properties(cpp_tests PROPERTIES WORKING_DIRECTORY ${CPP_SOURCE_ROOT})

add_test(NAME cpp_nomtln_system COMMAND cpp_tests --gtest_filter=SystemInitSolver.*)
set_tests_properties(cpp_nomtln_system PROPERTIES WORKING_DIRECTORY ${CPP_SOURCE_ROOT})

add_test(NAME cpp_planewave_mur COMMAND cpp_tests --gtest_filter=Planewave*:BordersMur*)
set_tests_properties(cpp_planewave_mur PROPERTIES WORKING_DIRECTORY ${CPP_SOURCE_ROOT})

add_test(NAME cpp_mpi_one_axis COMMAND cpp_tests --gtest_filter=MpiOneAxis.*)
set_tests_properties(cpp_mpi_one_axis PROPERTIES WORKING_DIRECTORY ${CPP_SOURCE_ROOT})

add_test(NAME cpp_pml_boundary COMMAND cpp_tests --gtest_filter=PmlBoundary.*)
set_tests_properties(cpp_pml_boundary PROPERTIES WORKING_DIRECTORY ${CPP_SOURCE_ROOT})

if(SEMBA_FDTD_ENABLE_MPI AND SEMBA_FDTD_ENABLE_MTLN)
    add_test(
        NAME cpp_mtln_mpi_bundle
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 $<TARGET_FILE:cpp_tests> --gtest_filter=MtlnMpiBundle.*
    )
    set_tests_properties(cpp_mtln_mpi_bundle PROPERTIES WORKING_DIRECTORY ${CPP_SOURCE_ROOT})
endif()

if(SEMBA_FDTD_ENABLE_MPI)
    add_test(
        NAME cpp_mpi_one_axis_exchange
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 $<TARGET_FILE:cpp_tests> --gtest_filter=MpiOneAxis.Exchanges*
    )
    set_tests_properties(cpp_mpi_one_axis_exchange PROPERTIES WORKING_DIRECTORY ${CPP_SOURCE_ROOT})

    add_test(
        NAME cpp_mpi_one_axis_exchange_4ranks
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 $<TARGET_FILE:cpp_tests> --gtest_filter=MpiOneAxis.Exchanges*
    )
    set_tests_properties(cpp_mpi_one_axis_exchange_4ranks PROPERTIES WORKING_DIRECTORY ${CPP_SOURCE_ROOT})
endif()
