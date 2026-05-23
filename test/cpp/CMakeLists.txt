message(STATUS "Creating build system for test/cpp")

add_executable(cpp_tests
    cpp_tests.cpp
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
    ${SEMBA_TEST_SUPPORT}
)

if(SEMBA_FDTD_ENABLE_MTLN)
    target_compile_definitions(cpp_tests PRIVATE CompileWithMTLN)
    target_link_libraries(cpp_tests mtlnsolver ngspice ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
endif()

if(SEMBA_FDTD_ENABLE_MPI)
    target_compile_definitions(cpp_tests PRIVATE CompileWithMPI)
else()
    target_compile_definitions(cpp_tests PRIVATE CompileWithoutMPI)
endif()

add_test(NAME cpp_tests COMMAND cpp_tests)
set_tests_properties(cpp_tests PROPERTIES WORKING_DIRECTORY ${CPP_SOURCE_ROOT})
