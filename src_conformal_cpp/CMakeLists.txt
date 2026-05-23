message(STATUS "Creating build system for conformal (C++)")

add_library(conformal
    "${CMAKE_CURRENT_SOURCE_DIR}/../src_cpp/src_conformal/cell_map.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/../src_cpp/src_conformal/conformal_build.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/../src_cpp/src_conformal/geometry.cpp"
)

target_include_directories(conformal PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/../src_cpp/src_main_pub
    ${CMAKE_CURRENT_SOURCE_DIR}/../src_cpp/src_conformal
)
target_link_libraries(conformal PUBLIC semba-types)
