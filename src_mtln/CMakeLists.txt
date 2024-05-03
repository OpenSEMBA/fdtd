cmake_minimum_required (VERSION 3.7)

project (mtlnsolver)
enable_language (Fortran)
enable_language (C)

message(STATUS "Creating build system for mtlnsolver")

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

add_library(ngspice_interface
    "ngspice_interface.c"
)

target_link_libraries(ngspice_interface PRIVATE	${NGSPICE_LIB})

if (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
	target_link_libraries(ngspice_interface PRIVATE	stdc++)
endif ()

add_library(mtlnsolver
	"mtl.F90"
	"mtl_bundle.F90"
	"mtln_solver.F90"
	"network.F90"
	"network_manager.F90"
	"utils.F90"
	"probes.F90"
	"dispersive.F90"
	"rational_approximation.F90"
	"circuit.F90"
	"ngspice_interface.F90"
	"preprocess.F90"
)

target_link_libraries(mtlnsolver PRIVATE
	${LAPACK_LIB} 
	${BLAS_LIB} 
	fhash
	ngspice_interface
	semba-types
)