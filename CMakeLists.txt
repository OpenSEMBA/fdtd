cmake_minimum_required (VERSION 3.15)

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
option(SEMBA_FDTD_ENABLE_INTEL_IPO "When compiling in Release, enables the interprocedural optimization" ON)

option(SEMBA_FDTD_EXECUTABLE "Compiles executable" ON)
option(SEMBA_FDTD_MAIN_LIB "Compiles main library" ON)
option(SEMBA_FDTD_COMPONENTS_LIB "Compiles components library" ON)
option(SEMBA_FDTD_OUTPUTS_LIB "Compiles outputs library" ON)

# Compilation defines.
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
	message(STATUS "Using Linux flags")

	if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
		message(STATUS "Using GNU flags")

		set(CMAKE_CXX_FLAGS "-fopenmp")
		set(CMAKE_Fortran_FLAGS "-fopenmp -ffree-form -ffree-line-length-none -fdec -fallow-argument-mismatch")
		
		set(CMAKE_C_FLAGS_RELEASE "-Ofast")
		set(CMAKE_CXX_FLAGS_RELEASE "-Ofast")
		set(CMAKE_Fortran_FLAGS_RELEASE "-Ofast")
		
		set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0")

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

	elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
		message(STATUS "Using Intel (ifort) flags")

		set(CMAKE_C_FLAGS "-diag-disable=10441")
		set(CMAKE_CXX_FLAGS "-qopenmp -diag-disable=10441")
		set(CMAKE_Fortran_FLAGS "-fpp -qopenmp -diag-disable=8889 -static-intel")

		set(CMAKE_C_FLAGS_RELEASE "-O3")
		set(CMAKE_CXX_FLAGS_RELEASE "-O3")
		set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
		set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ON)
		if(SEMBA_FDTD_ENABLE_INTEL_XHOST_OPTIMIZATION) 
			set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -xHost")
			set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xHost")
			set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -xHost")
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

		# This prevents a compilation error happening due to visual studio headers being too strict.
		# It can possibly can be removed in the future.
		add_compile_definitions("$<$<COMPILE_LANGUAGE:C,CXX>:-D_ALLOW_COMPILER_AND_STL_VERSION_MISMATCH>")

		set(CMAKE_CXX_FLAGS "/Qopenmp")
		set(CMAKE_Fortran_FLAGS "/Qopenmp /fpp")

		# Not tested

		# set(CMAKE_Fortran_FLAGS_RELEASE "/fast")
		
		set(CMAKE_CXX_FLAGS_FLAGS_DEBUG "/Qopenmp-simd- /Wno-debug-option-simd")
		# set(CMAKE_Fortran_FLAGS_DEBUG "/check all,nouninit /debug full -traceback")

	elseif((CMAKE_Fortran_COMPILER_ID MATCHES "Intel"))
		message(STATUS "Using Intel (ifort) flags")
		set(CMAKE_CXX_FLAGS "/Qopenmp")

		# Not tested
		# set(CMAKE_Fortran_FLAGS "/fpp /qopenmp") 

		# set(CMAKE_Fortran_FLAGS_RELEASE "/fast")

		# set(CMAKE_Fortran_FLAGS_DEBUG "/check all,nouninit /debug full -traceback")
	else()
		message(FATAL_ERROR "Unrecognized compiler id:" ${CMAKE_Fortran_COMPILER_ID})
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
target_link_libraries(semba-types ${MPI_Fortran_LIBRARIES})

add_library(semba-reports
	"src_main_pub/errorreport.F90"
	"src_main_pub/snapxdmf.F90"
)
target_link_libraries(semba-reports semba-types ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})

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
	target_link_libraries(semba-components semba-types semba-reports ${MTLN_LIBRARIES})
endif()

if(SEMBA_FDTD_OUTPUTS_LIB)
	add_library(semba-outputs
		"src_main_pub/mpicomm.F90"
		"src_main_pub/observation.F90"
		"src_main_pub/vtk.F90"
		"src_main_pub/xdmf.F90"
		"src_main_pub/xdmf_h5.F90"
	)
	target_link_libraries(semba-outputs 
		semba-components 
		${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES} 
		${MPI_Fortran_LIBRARIES})
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
	target_link_libraries(semba-main
		semba-outputs 
		${SMBJSON_LIBRARIES} 
		${MTLN_LIBRARIES})
endif()

if (SEMBA_FDTD_EXECUTABLE)
	add_executable(semba-fdtd 
		"src_main_pub/semba_fdtd.F90"
	)
	target_link_libraries(semba-fdtd semba-main semba-reports)
	target_link_libraries(semba-fdtd ${MPI_Fortran_LIBRARIES})
endif()


include_directories(${CMAKE_BINARY_DIR}/mod)
include_directories(${HDF5_INCLUDE_DIRS})
include_directories(${FHASH_INCLUDES})
