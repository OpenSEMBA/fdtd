cmake_minimum_required (VERSION 3.15)

# Allow choosing between Fortran and C++ builds
option(SEMBA_FDTD_BUILD_CXX "Build C++ version instead of Fortran" OFF)

if(SEMBA_FDTD_BUILD_CXX)
    message(STATUS "Building C++ version of semba-fdtd")
    include("${CMAKE_CURRENT_SOURCE_DIR}/cpp_build.cmake")
    return()
endif()

project(semba-fdtd Fortran)
enable_language (Fortran)


set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

message(STATUS "Compiler Id is: ${CMAKE_Fortran_COMPILER_ID}")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

option(SEMBA_FDTD_ENABLE_MPI "Use MPI" OFF)
option(SEMBA_FDTD_ENABLE_HDF "Use HDF" ON)
option(SEMBA_FDTD_ENABLE_MTLN "Use MTLN" ON)
option(SEMBA_FDTD_ENABLE_SMBJSON "Use smbjson" ON)
option(SEMBA_FDTD_ENABLE_DOUBLE_PRECISION "Use double precision (CompileWithReal8)" OFF)

option(SEMBA_FDTD_ENABLE_TEST "Compile tests" ON)

option(SEMBA_FDTD_ENABLE_INTEL_XHOST_OPTIMIZATION "When compiling in Release, enables the -xHost optimization flag (not supported in github actions)" OFF) 
option(SEMBA_FDTD_ENABLE_INTEL_IPO "When compiling in Release, enables interprocedural optimization" OFF)

option(SEMBA_FDTD_EXECUTABLE "Compiles executable" ON)
option(SEMBA_FDTD_MAIN_LIB "Compiles main library" ON)
option(SEMBA_FDTD_COMPONENTS_LIB "Compiles components library" ON)
option(SEMBA_FDTD_OUTPUTS_LIB "Compiles outputs library" ON)
# Compilation defines.
if(CMAKE_BUILD_TYPE MATCHES "Release" OR CMAKE_BUILD_TYPE MATCHES "release" )
	add_definitions(-DCompileWithRelease)
else()
	add_definitions(-DCompileWithDebug)
endif()
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

include("${CMAKE_CURRENT_SOURCE_DIR}/set_precompiled_libraries.cmake")

if (CMAKE_SYSTEM_NAME MATCHES "Linux")
	message(STATUS "Using Linux flags")

	if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
		message(STATUS "Using GNU flags")

		set(CMAKE_CXX_FLAGS "-fopenmp")
		set(CMAKE_Fortran_FLAGS "-fopenmp -ffree-form -ffree-line-length-none -fdec -fallow-argument-mismatch")
		
		set(CMAKE_C_FLAGS_RELEASE "-Ofast")
		set(CMAKE_CXX_FLAGS_RELEASE "-Ofast")
		set(CMAKE_Fortran_FLAGS_RELEASE "-Ofast")

		set(CMAKE_C_FLAGS_DEBUG "-g -O0")
		set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
		set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -fno-inline -fcheck=all -fbacktrace")

	elseif(CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM")
		message(STATUS "Using IntelLLVM (ifx) flags")
		
		set(CMAKE_CXX_FLAGS "-qopenmp")
		set(CMAKE_Fortran_FLAGS "-qopenmp -fpp -static-intel")
		
		if (CMAKE_BUILD_TYPE STREQUAL "Release")
			set(CMAKE_C_FLAGS_RELEASE "-O3")
			set(CMAKE_CXX_FLAGS_RELEASE "-O3")
			set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -fp-model fast=2w")

			if(SEMBA_FDTD_ENABLE_INTEL_IPO)
				include(CheckIPOSupported)
				check_ipo_supported(RESULT iporesult)
				if(iporesult)
				  message (STATUS "IPO is supported and being enabled")
				else()
				  message(FATAL_ERROR "IPO is not supported")
				endif()
				set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ON)
				set(CMAKE_C_FLAGS_RELEASE "-ipo")
				set(CMAKE_CXX_FLAGS_RELEASE "-ipo")
				set(CMAKE_Fortran_FLAGS_RELEASE "-ipo")
			else()
				set(CMAKE_INTERPROCEDURAL_OPTIMIZATION OFF)
			endif()
		
			if(SEMBA_FDTD_ENABLE_INTEL_XHOST_OPTIMIZATION) 
				set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -xHost")
				set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -xHost")
				set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -xHost")
			endif()
		endif()
		
		set(CMAKE_Fortran_FLAGS_DEBUG "-check all,nouninit -debug full -traceback")

	
	elseif(CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC")
		message(STATUS "Using nvhpc flags")
		
		# TODO: Tune flags for NVHPC
	else()
		message(FATAL_ERROR "Unrecognized compiler:" ${CMAKE_Fortran_COMPILER_ID})
	endif()
elseif(CMAKE_SYSTEM_NAME MATCHES "Windows")
	message(STATUS "Using Windows flags")
	if(CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM")
		message(STATUS "Using IntelLLVM (ifx) flags")

		# Debug flags - enable full debugging support for breakpoints
		set(CMAKE_Fortran_FLAGS_DEBUG "/debug:full /Od /traceback /check:all,nouninit /Zi /debug-parameters:all /warn:interfaces")
		set(CMAKE_CXX_FLAGS_DEBUG "/Zi /Od /RTC1")
		set(CMAKE_C_FLAGS_DEBUG "/Zi /Od /RTC1")
		
		# Release flags
		set(CMAKE_Fortran_FLAGS_RELEASE "/O3")
		set(CMAKE_CXX_FLAGS_RELEASE "/O2")
		set(CMAKE_C_FLAGS_RELEASE "/O2")

	else()
		message(FATAL_ERROR "Unrecognized compiler id:" ${CMAKE_Fortran_COMPILER_ID})
	endif()
else()
	message(FATAL_ERROR "Unrecognized system name")
endif()

# =========================================

add_subdirectory(external)

add_library(semba-types 
	"src/main/nfde_types.F90"
	"src/main/fdetypes.F90"
	"src/mtln/mtln_types.F90"
	"src/conformal/conformal_types.F90"
	"src/wires/wires_types.F90"
	"src/main/lumped_types.F90"
)
target_link_libraries(semba-types ${MPI_Fortran_LIBRARIES})

add_library(semba-reports
	"src/main/errorreport.F90"
	"src/main/snapxdmf.F90"
)
target_link_libraries(semba-reports semba-types ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})

if(SEMBA_FDTD_ENABLE_SMBJSON)
	add_subdirectory(src/json_parser)
	if(PROJECT_IS_TOP_LEVEL)
		set(SMBJSON_LIBRARIES smbjson)
	else()
		set(SMBJSON_LIBRARIES smbjson PARENT_SCOPE)
	endif()
endif()
if (SEMBA_FDTD_ENABLE_MTLN)
	set(NGSPICE_LIB ngspice)
	add_definitions(-DCompileWithMTLN)
	add_subdirectory(src/mtln)
	if (PROJECT_IS_TOP_LEVEL)
		set(MTLN_LIBRARIES mtlnsolver ngspice_interface)
	else()
		set(MTLN_LIBRARIES mtlnsolver ngspice_interface PARENT_SCOPE)
	endif()
endif()

add_subdirectory(src/conformal)
set(CONFORMAL_LIBRARIES conformal)

if (SEMBA_FDTD_ENABLE_TEST)
	add_subdirectory(external/googletest/)
	add_subdirectory(test/fortran)
endif()

if(SEMBA_FDTD_COMPONENTS_LIB)
	add_library(semba-components
		"src/main/anisotropic.F90"
		"src/main/borderscpml.F90"
		"src/main/bordersmur.F90"
		"src/main/bordersother.F90"
		"src/main/electricdispersive.F90"
		"src/main/magneticdispersive.F90"
		"src/main/nodalsources.F90"
		"src/main/planewaves.F90"
		"src/main/pml_bodies.F90"
		"src/main/maloney_nostoch.F90"
		"src/main/lumped.F90"
		"src/main/dmma_thin_slot.F90"
		"src/main/farfield.F90"
		"src/wires/wires.F90"
		"src/wires/wires_mtln.F90"
	)
	target_link_libraries(semba-components semba-types semba-reports ${MTLN_LIBRARIES})
endif()

if(SEMBA_FDTD_OUTPUTS_LIB)
	add_library(semba-outputs
		"src/main/mpicomm.F90"
		"src/main/observation.F90"
		"src/main/vtk.F90"
		"src/main/xdmf.F90"
		"src/main/xdmf_h5.F90"
	)
	target_link_libraries(semba-outputs 
		semba-components 
		${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES} 
		${MPI_Fortran_LIBRARIES})
endif()

if(SEMBA_FDTD_MAIN_LIB)
	set(PROGRAM_NAME "semba-fdtd")
	include("${CMAKE_CURRENT_SOURCE_DIR}/get_commit_info.cmake")
	configure_file(
	${CMAKE_SOURCE_DIR}/src/main/version.F90.in
	${CMAKE_SOURCE_DIR}/src/main/version.F90
	@ONLY
	)

	add_library(semba-main
		"src/main/semba_fdtd.F90"
		"src/main/calc_constants.F90"
		"src/main/nfde_rotate.F90"
		"src/main/EpsMuTimeScale.F90"
		"src/main/getargs.F90"
		"src/main/healer.F90"
		"src/main/preprocess_geom.F90"
		"src/main/storegeom.F90"
		"src/main/version.F90"
		"src/main/postprocess.F90"
		"src/main/interpreta_switches.F90"
		"src/main/resuming.F90"	
		"src/main/timestepping.F90"
	)
	target_link_libraries(semba-main
		semba-outputs 
		${SMBJSON_LIBRARIES} 
		${MTLN_LIBRARIES})
endif()

if (SEMBA_FDTD_EXECUTABLE)
	add_executable(semba-fdtd 
		"src/main/launcher.F90"
	)
	target_link_libraries(semba-fdtd semba-main semba-reports)
	target_link_libraries(semba-fdtd ${MPI_Fortran_LIBRARIES})
endif()



include_directories(${CMAKE_BINARY_DIR}/mod)
include_directories(${HDF5_INCLUDE_DIRS})
include_directories(${FHASH_INCLUDES})
