message(STATUS "Creating build system for test/cpp")

add_executable(cpp_tests
    cpp_tests.cpp
)

target_include_directories(cpp_tests PRIVATE
    ${CPP_SOURCE_ROOT}/src_cpp/src_json_parser
    ${CPP_SOURCE_ROOT}/src_cpp/src_conformal
    ${CPP_SOURCE_ROOT}/src_cpp/src_main_pub
    ${CPP_SOURCE_ROOT}/test
)

target_link_libraries(cpp_tests
    gtest_main
    nlohmann_json::nlohmann_json
    smbjson
)

# Add test to CTest
add_test(NAME cpp_tests COMMAND cpp_tests)
