cmake_minimum_required (VERSION 3.15)

project(semba-fdtd Fortran)
enable_language (Fortran)


set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

message(STATUS ${CMAKE_Fortran_COMPILER_ID})

option(SEMBA_FDTD_ENABLE_MPI "Use MPI" OFF)
option(SEMBA_FDTD_ENABLE_HDF "Use HDF" ON)
option(SEMBA_FDTD_ENABLE_MTLN "Use MTLN" ON)
option(SEMBA_FDTD_ENABLE_SMBJSON "Use smbjson" ON)
option(SEMBA_FDTD_ENABLE_DOUBLE_PRECISION "Use double precision (CompileWithReal8)" OFF)
option(SEMBA_FDTD_ENABLE_TEST "Compile tests" ON)

option(SEMBA_FDTD_EXECUTABLE "Compiles executable" ON)
option(SEMBA_FDTD_MAIN_LIB "Compiles main library" ON)
option(SEMBA_FDTD_COMPONENTS_LIB "Compiles components library" ON)
option(SEMBA_FDTD_OUTPUTS_LIB "Compiles outputs library" ON)

if(SEMBA_FDTD_ENABLE_SMBJSON)
	add_definitions(-DCompileWithSMBJSON)
endif()
if (SEMBA_FDTD_ENABLE_MTLN)
	add_definitions(-DCompileWithMTLN)
endif()
if (SEMBA_FDTD_ENABLE_DOUBLE_PRECISION)
	add_definitions(-DCompileWithReal8)
else()
	add_definitions(-DCompileWithReal4)
endif()
add_definitions(		
	-DCompileWithInt2 
	-DCompileWithOpenMP
)
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
	add_definitions(-DCompileWithDebug)
endif()

include("${CMAKE_CURRENT_SOURCE_DIR}/set_precompiled_libraries.cmake")

if (CMAKE_SYSTEM_NAME MATCHES "Linux")
	message(STATUS "Using flags for Linux system")
	if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
		add_definitions(-DCompileWithIncludeMpifh)
		set(CMAKE_CXX_FLAGS "-fopenmp")
		set(CMAKE_Fortran_FLAGS "-fopenmp -ffree-form -ffree-line-length-none -fdec -fallow-argument-mismatch")
		set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -ffast-math -fPIC")
		set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0")
	elseif(CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM")
		set(CMAKE_CXX_FLAGS "-qopenmp -qmkl")
		# set(CMAKE_CXX_FLAGS_RELEASE "-fast") # This causes problems because it activates -ipo.
		set(CMAKE_CXX_FLAGS_RELEASE "-O3 -xHost")
		set(CMAKE_Fortran_FLAGS "-fpp -qopenmp -qmkl")
		# set(CMAKE_Fortran_FLAGS_RELEASE "-fast") # This causes problems because it activates -ipo.
		set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -no-prec-div -fp-model fast=2w -xHost")
		set(CMAKE_Fortran_FLAGS_DEBUG "-check all,nouninit -debug full -traceback")
	elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
		set(CMAKE_CXX_FLAGS "-qopenmp -mkl")
		set(CMAKE_CXX_FLAGS_RELEASE "-O3 -xHost")
		set(CMAKE_Fortran_FLAGS "-fpp -qopenmp -mkl")
		set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -xHost")
		set(CMAKE_Fortran_FLAGS_DEBUG "-check all,nouninit -debug full -traceback")
	elseif(CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC")
		# TODO: Tune flags for NVHPC
	else()
		message(FATAL_ERROR "Unrecognized compiler:" ${CMAKE_Fortran_COMPILER_ID})
	endif()
elseif(CMAKE_SYSTEM_NAME MATCHES "Windows")
	if(CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM")
		set(CMAKE_CXX_FLAGS "/Qiopenmp /qmkl")
		set(CMAKE_Fortran_FLAGS "/fpp /qopenmp /qmkl")
		set(CMAKE_Fortran_FLAGS_RELEASE "/fast")
		set(CMAKE_Fortran_FLAGS_DEBUG "/check all,nouninit /debug full -traceback")
	elseif((CMAKE_Fortran_COMPILER_ID MATCHES "Intel"))
		set(CMAKE_CXX_FLAGS "/Qiopenmp /mkl")
		set(CMAKE_Fortran_FLAGS "/fpp /qopenmp /mkl")
		set(CMAKE_Fortran_FLAGS_RELEASE "/fast")
		set(CMAKE_Fortran_FLAGS_DEBUG "/check all,nouninit /debug full -traceback")
	else()
		message(FATAL_ERROR "Unrecognized compiler:" ${CMAKE_Fortran_COMPILER_ID})
	endif()
else()
	message(FATAL_ERROR "Unrecognized system name")
endif()

add_subdirectory(external)

add_library(semba-types 
	"src_main_pub/nfde_types.F90"
	"src_main_pub/fdetypes.F90"
	"src_mtln/mtln_types.F90"
	"src_wires_pub/wires_types.F90"
	"src_main_pub/lumped_types.F90"
)

if(SEMBA_FDTD_ENABLE_SMBJSON)
	add_subdirectory(src_json_parser)
	if(PROJECT_IS_TOP_LEVEL)
		set(SMBJSON_LIBRARIES smbjson)
	else()
		set(SMBJSON_LIBRARIES smbjson PARENT_SCOPE)
	endif()
endif()
if (SEMBA_FDTD_ENABLE_MTLN)
	set(NGSPICE_LIB ngspice)
	add_definitions(-DCompileWithMTLN)
	add_subdirectory(src_mtln)
	if (PROJECT_IS_TOP_LEVEL)
		set(MTLN_LIBRARIES mtlnsolver ngspice_interface)
	else()
		set(MTLN_LIBRARIES mtlnsolver ngspice_interface PARENT_SCOPE)
	endif()
endif()

if (SEMBA_FDTD_ENABLE_TEST)
	add_subdirectory(external/googletest/)
	add_subdirectory(test)
endif()

if(SEMBA_FDTD_COMPONENTS_LIB)
	add_library(semba-components
		"src_main_pub/errorreport.F90"
		"src_main_pub/anisotropic.F90"
		"src_main_pub/borderscpml.F90"
		"src_main_pub/bordersmur.F90"
		"src_main_pub/bordersother.F90"
		"src_main_pub/electricdispersive.F90"
		"src_main_pub/magneticdispersive.F90"
		"src_main_pub/nodalsources.F90"
		"src_main_pub/planewaves.F90"
		"src_main_pub/pml_bodies.F90"
		"src_main_pub/maloney_nostoch.F90"
		"src_main_pub/lumped.F90"
		"src_main_pub/dmma_thin_slot.F90"
		"src_main_pub/farfield.F90"
		"src_wires_pub/wires.F90"
		"src_wires_pub/wires_mtln.F90"
	)
	target_link_libraries(semba-components semba-types ${MTLN_LIBRARIES})
endif()

if(SEMBA_FDTD_OUTPUTS_LIB)
	add_library(semba-outputs
		"src_main_pub/mpicomm.F90"
		"src_main_pub/observation.F90"
		"src_main_pub/vtk.F90"
		"src_main_pub/snapxdmf.F90"
		"src_main_pub/xdmf.F90"
		"src_main_pub/xdmf_h5.F90"
	)
	target_link_libraries(semba-outputs semba-components ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})
endif()

if(SEMBA_FDTD_MAIN_LIB)
	add_library(semba-main
		"src_main_pub/calc_constants.F90"
		"src_main_pub/nfde_rotate.F90"
		"src_main_pub/EpsMuTimeScale.F90"
		"src_main_pub/getargs.F90"
		"src_main_pub/healer.F90"
		"src_main_pub/preprocess_geom.F90"
		"src_main_pub/storegeom.F90"
		"src_main_pub/version.F90"
		"src_main_pub/postprocess.F90"
		"src_main_pub/interpreta_switches.F90"
		"src_main_pub/resuming.F90"	
		"src_main_pub/timestepping.F90"
	)
	target_link_libraries(semba-main semba-outputs ${SMBJSON_LIBRARIES} ${MTLN_LIBRARIES})
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