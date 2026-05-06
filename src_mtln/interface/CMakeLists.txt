message(STATUS "Creating build system for ngspice interface")

project(ngspice-interface)
enable_language(C CXX)

add_library(ngspice_interface
    "ngspice_interface.c"
)
target_link_libraries(ngspice_interface PRIVATE ${NGSPICE_LIB})
target_include_directories(ngspice_interface PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/../../external/ngspice/src/include"
)
if (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
	target_link_libraries(ngspice_interface PRIVATE	stdc++)
endif ()


