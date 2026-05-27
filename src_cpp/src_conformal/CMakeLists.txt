message(STATUS "Creating build system for conformal (C++)")

add_library(conformal
    cell_map.cpp
    conformal_build.cpp
    geometry.cpp
)

target_include_directories(conformal PUBLIC
    ${CPP_SOURCE_ROOT}/src_cpp/src_main_pub
    ${CMAKE_CURRENT_SOURCE_DIR}
)
target_link_libraries(conformal PUBLIC semba-types)
target_compile_features(conformal PUBLIC cxx_std_17)
