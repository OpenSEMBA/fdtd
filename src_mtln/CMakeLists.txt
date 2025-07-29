message(STATUS "Creating build system for mtlnsolver")

project(mtln)
enable_language(C CXX Fortran)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

add_subdirectory(interface)

add_library(mtlnsolver
	"mtl.F90"
	"mtl_bundle.F90"
	"mtln_solver.F90"
	"multipolar_expansion.F90"
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
	${MPI_Fortran_LIBRARIES}
)