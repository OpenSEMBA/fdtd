add_library(semba-test-support STATIC
    observation_preprocess.cpp
    observation_movie_test.cpp
)

target_include_directories(semba-test-support PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/../src_main_pub
    ${CMAKE_CURRENT_SOURCE_DIR}/../src_conformal
    ${CMAKE_CURRENT_SOURCE_DIR}/../src_mtln
)

target_link_libraries(semba-test-support PUBLIC semba-types)
