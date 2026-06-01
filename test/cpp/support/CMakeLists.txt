add_library(semba-test-support STATIC
    observation_preprocess.cpp
    observation_movie_test.cpp
)

target_include_directories(semba-test-support PUBLIC
    ${CPP_SOURCE_ROOT}/src_cpp/main
    ${CPP_SOURCE_ROOT}/src_cpp/conformal
    ${CPP_SOURCE_ROOT}/src_cpp/mtln
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(semba-test-support PUBLIC semba-types)
