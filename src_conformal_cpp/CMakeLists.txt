message(STATUS "Creating build system for conformal (C++)")

add_library(conformal
    "${CMAKE_CURRENT_SOURCE_DIR}/../src_cpp/src_conformal/cell_map.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/../src_cpp/src_conformal/conformal.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/../src_cpp/src_conformal/geometry.cpp"
)

target_link_libraries(conformal PRIVATE semba-types)
