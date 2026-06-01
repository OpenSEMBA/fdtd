cmake_minimum_required(VERSION 3.15)

add_library(smbjson
    cells.cpp
    smbjson_labels.cpp
    smbjson.cpp
)

if(SEMBA_FDTD_ENABLE_MTLN)
    target_sources(smbjson PRIVATE smbjson_mtln.cpp)
    target_compile_definitions(smbjson PUBLIC CompileWithMTLN)
endif()

target_include_directories(smbjson PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CPP_SOURCE_ROOT}/src_cpp/main
    ${CPP_SOURCE_ROOT}/src_cpp/conformal
    ${CPP_SOURCE_ROOT}/src_cpp/mtln
    ${CPP_SOURCE_ROOT}/src_cpp/wires
)

target_compile_features(smbjson PUBLIC cxx_std_17)

find_package(nlohmann_json QUIET)
if(nlohmann_json_FOUND)
    target_link_libraries(smbjson PUBLIC nlohmann_json::nlohmann_json)
endif()
