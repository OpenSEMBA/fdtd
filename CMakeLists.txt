cmake_minimum_required (VERSION 3.18)

enable_language (CXX Fortran)
project(semba-fdtd Fortran)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

message(STATUS ${CMAKE_Fortran_COMPILER_ID})

option(CompileWithHDF "Use HDF" OFF)
option(CompileWithMPI "Use MPI" OFF)

include("${CMAKE_CURRENT_SOURCE_DIR}/set_precompiled_libraries.cmake")


if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
	set(CMAKE_Fortran_FLAGS "-ffree-form -ffree-line-length-none -fdec -fcheck=bounds -fopenmp -fallow-argument-mismatch")
	add_definitions(-DCompileWithIncludeMpifh)
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "^Intel")
	set(CMAKE_Fortran_FLAGS "-fpp -check all,nouninit -debug full -traceback -diag-disable=10441")
	set(CMAKE_CXX_FLAGS "-diag-disable=10448")
endif()



add_library(semba-types 
	"src_main_pub/nfde_types.F90"
	"src_main_pub/fdetypes.F90"
	"src_mtln/mtln_types.F90"
)
add_subdirectory(external/fhash/)
set(FHASH_INCLUDES  ${CMAKE_BINARY_DIR}/external/fhash/include/fortran_modules/${CMAKE_Fortran_COMPILER_ID}/${CMAKE_Fortran_COMPILER_VERSION})

add_subdirectory(external/googletest/)

add_subdirectory(src_json_parser)
add_subdirectory(src_mtln)

add_subdirectory(test)

add_library(semba-main
	"src_main_pub/errorreport.F90"
	"src_main_pub/anisotropic.F90"
	"src_main_pub/borderscpml.F90"
	"src_main_pub/bordersmur.F90"
	"src_main_pub/bordersother.F90"
	"src_main_pub/calc_constants.F90"
	"src_main_pub/dmma_thin_slot.F90"
	"src_main_pub/electricdispersive.F90"
	"src_main_pub/EpsMuTimeScale.F90"
	"src_main_pub/farfield.F90"
	"src_main_pub/getargs.F90"
	"src_main_pub/healer.F90"
	"src_main_pub/lumped_types.F90"
	"src_main_pub/lumped.F90"
	"src_main_pub/magneticdispersive.F90"
	"src_main_pub/maloney_nostoch.F90"
	"src_main_pub/mpicomm.F90"
	"src_main_pub/nfde_rotate.F90"
	"src_main_pub/nodalsources.F90"
	"src_main_pub/observation.F90"
	"src_main_pub/planewaves.F90"
	"src_main_pub/pml_bodies.F90"
	"src_main_pub/postprocess.F90"
	"src_main_pub/preprocess_geom.F90"
	"src_main_pub/resuming.F90"
	"src_main_pub/interpreta_switches.F90"
	"src_main_pub/storegeom.F90"
	"src_main_pub/timestepping.F90"
	"src_main_pub/version.F90"
	"src_wires_pub/wires_types.F90"
	"src_wires_pub/wires.F90"
	"src_wires_pub/wires_mtln.F90"
)

target_link_libraries(semba-main semba-types smbjson mtlnsolver ngspice_interface)
target_link_libraries(semba-main ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})
target_link_libraries(semba-main ${MPI_Fortran_LIBRARIES})

add_executable(semba-fdtd 
	"src_main_pub/semba_fdtd.F90"
)
target_link_libraries(semba-fdtd semba-main)

include_directories(${CMAKE_BINARY_DIR}/mod)
include_directories(${HDF5_INCLUDE_DIRS})
include_directories(${FHASH_INCLUDES})

add_definitions(		
-DCompileWithInt2 
-DCompileWithReal4 
-DCompileWithOpenMP 
-DCompileWithAnisotropic 
-DCompileWithEDispersives
-DCompileWithNF2FF
-DCompileWithNodalSources
-DCompileWithDMMA
-DCompileWithSGBC
-DCompileWithWires
-DCompileWithMTLN
)


