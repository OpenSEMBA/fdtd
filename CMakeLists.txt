cmake_minimum_required (VERSION 3.15)

project(semba-fdtd Fortran)
enable_language (CXX C Fortran)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

message(STATUS ${CMAKE_Fortran_COMPILER_ID})

option(SEMBA_FDTD_ENABLE_MPI "Use MPI" OFF)
option(SEMBA_FDTD_ENABLE_HDF "Use HDF" ON)
option(SEMBA_FDTD_ENABLE_MTLN "Use MTLN" ON)
option(SEMBA_FDTD_ENABLE_SMBJSON "Use smbjson" ON)
option(SEMBA_FDTD_ENABLE_TEST "Compile tests" ON)

option(SEMBA_FDTD_MAIN_LIB   "Compiles main library" ON)
option(SEMBA_FDTD_EXECUTABLE "Compiles executable" ON)
if (SEMBA_FDTD_EXECUTABLE)
	set(SEMBA_FDTD_MAIN_LIB ON)
else()
	set(SEMBA_FDTD_MAIN_LIB OFF)
endif()

include("${CMAKE_CURRENT_SOURCE_DIR}/set_precompiled_libraries.cmake")

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
	set(CMAKE_Fortran_FLAGS "-ffree-form -ffree-line-length-none -fdec -fcheck=bounds -fopenmp -fallow-argument-mismatch")
	set(CMAKE_CXX_FLAGS "-fopenmp")
	add_definitions(-DCompileWithIncludeMpifh)
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "^Intel")
	set(CMAKE_Fortran_FLAGS "-fpp")
	if (CMAKE_SYSTEM_NAME MATCHES "Linux")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qopenmp")
	endif()
	if(CMAKE_BUILD_TYPE MATCHES "Debug")
		set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -check all,nouninit -debug full -traceback")
	endif()
endif()

if (SEMBA_FDTD_ENABLE_MTLN OR SEMBA_FDTD_ENABLE_SMBJSON)
set (FHASH_SOURCES_FOLDER "external/fhash/src/")
add_library(fhash 
	${FHASH_SOURCES_FOLDER}/fhash.f90 
	${FHASH_SOURCES_FOLDER}/fhash_sll.f90 
	${FHASH_SOURCES_FOLDER}/fhash_data_container.f90 
	${FHASH_SOURCES_FOLDER}/fhash_tbl.f90 
	${FHASH_SOURCES_FOLDER}/fhash_tbl_iter.f90 
	${FHASH_SOURCES_FOLDER}/fhash_fnv.f90 
	${FHASH_SOURCES_FOLDER}/fhash_key/base.f90 
	${FHASH_SOURCES_FOLDER}/fhash_key/int32.f90 
	${FHASH_SOURCES_FOLDER}/fhash_key/int64.f90 
	${FHASH_SOURCES_FOLDER}/fhash_key/char.f90 
	${FHASH_SOURCES_FOLDER}/fhash_key/int32_1d.f90 
	${FHASH_SOURCES_FOLDER}/fhash_key/int64_1d.f90
)
endif()

add_library(semba-types 
	"src_main_pub/nfde_types.F90"
	"src_main_pub/fdetypes.F90"
	"src_mtln/mtln_types.F90"
	"src_wires_pub/wires_types.F90"
	"src_main_pub/lumped_types.F90"
)

if(SEMBA_FDTD_ENABLE_SMBJSON)
	add_definitions(-DCompileWithSMBJSON)
	add_subdirectory(src_json_parser)
	set(SMBJSON_LIBRARIES smbjson)
endif()
if (SEMBA_FDTD_ENABLE_MTLN)
	add_definitions(-DCompileWithMTLN)
	add_subdirectory(src_mtln)
	set(MTLN_LIBRARIES mtlnsolver ngspice_interface)
endif()

if (SEMBA_FDTD_ENABLE_TEST)
	add_subdirectory(external/googletest/)
	add_subdirectory(test)
endif()

if (SEMBA_FDTD_MAIN_LIB)
	add_library(semba-main
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
		"src_main_pub/lumped.F90"
		"src_main_pub/magneticdispersive.F90"
		"src_main_pub/maloney_nostoch.F90"
		"src_main_pub/mpicomm.F90"
		"src_main_pub/nfde_rotate.F90"
		"src_main_pub/nodalsources.F90"
		"src_main_pub/planewaves.F90"
		"src_main_pub/pml_bodies.F90"
		"src_main_pub/preprocess_geom.F90"
		"src_main_pub/storegeom.F90"
		"src_main_pub/snapxdmf.F90"
		"src_main_pub/version.F90"
		"src_main_pub/postprocess.F90"
		"src_main_pub/errorreport.F90"
		"src_main_pub/interpreta_switches.F90"
		"src_main_pub/resuming.F90"	
		"src_main_pub/timestepping.F90"
		"src_main_pub/observation.F90"
		"src_wires_pub/wires.F90"
		"src_wires_pub/wires_mtln.F90"
	)
	target_link_libraries(semba-main semba-types ${SMBJSON_LIBRARIES} ${MTLN_LIBRARIES})
	target_link_libraries(semba-main ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})
endif()

if (SEMBA_FDTD_EXECUTABLE)
	add_executable(semba-fdtd 
		"src_main_pub/semba_fdtd.F90"
	)
	target_link_libraries(semba-fdtd semba-main)
	target_link_libraries(semba-fdtd ${MPI_Fortran_LIBRARIES})
endif()


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
)


