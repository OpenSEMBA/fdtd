message(STATUS "Creating build system for test/cpp")

add_executable(cpp_tests
    cpp_tests.cpp
)

target_sources(cpp_tests PRIVATE
    ${CPP_SOURCE_ROOT}/src_cpp/src_main_pub/maloney_nostoch.cpp
)

target_include_directories(cpp_tests PRIVATE
    ${CPP_SOURCE_ROOT}/src_cpp/src_json_parser
    ${CPP_SOURCE_ROOT}/src_cpp/src_conformal
    ${CPP_SOURCE_ROOT}/src_cpp/src_main_pub
    ${CPP_SOURCE_ROOT}/src_cpp/src_mtln
    ${CPP_SOURCE_ROOT}/src_cpp/src_test_support
    ${CPP_SOURCE_ROOT}/test/cpp
)

target_link_libraries(cpp_tests
    gtest_main
    nlohmann_json::nlohmann_json
    smbjson
    conformal
    semba-main
    ${SEMBA_TEST_SUPPORT}
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

add_test(NAME cpp_tests COMMAND cpp_tests)
set_tests_properties(cpp_tests PROPERTIES WORKING_DIRECTORY ${CPP_SOURCE_ROOT})

add_test(NAME cpp_nomtln_system COMMAND cpp_tests --gtest_filter=SystemInitSolver.*)
set_tests_properties(cpp_nomtln_system PROPERTIES WORKING_DIRECTORY ${CPP_SOURCE_ROOT})

add_test(NAME cpp_planewave_mur COMMAND cpp_tests --gtest_filter=Planewave*:BordersMur*)
set_tests_properties(cpp_planewave_mur PROPERTIES WORKING_DIRECTORY ${CPP_SOURCE_ROOT})
