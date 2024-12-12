message(STATUS "Creating build system for ngspice interface")

project(ngspice-interface)
enable_language(C CXX)

add_library(ngspice_interface
    "ngspice_interface.c"
)

target_link_libraries(ngspice_interface PRIVATE ${NGSPICE_LIB})

