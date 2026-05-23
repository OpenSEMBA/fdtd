message(STATUS "Creating build system for smbjson (C++)")

set(SMBJSON_CPP_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/../src_cpp/src_json_parser/smbjson_labels.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/../src_cpp/src_json_parser/cells.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/../src_cpp/src_json_parser/smbjson.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/../src_cpp/src_json_parser/idchildtable.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/../src_cpp/src_json_parser/mesh.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/../src_cpp/src_json_parser/parser_tools.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/../src_cpp/src_json_parser/nfdetypes_extension.cpp"
)
if(SEMBA_FDTD_ENABLE_MTLN)
    list(APPEND SMBJSON_CPP_SOURCES
        "${CMAKE_CURRENT_SOURCE_DIR}/../src_cpp/src_json_parser/smbjson_mtln.cpp")
endif()
add_library(smbjson ${SMBJSON_CPP_SOURCES})

target_link_libraries(smbjson PUBLIC
    nlohmann_json::nlohmann_json
    semba-types
    semba-reports
    conformal
)
