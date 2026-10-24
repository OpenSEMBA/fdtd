cmake_minimum_required (VERSION 3.18)

project(semba-fdtd LANGUAGES Fortran C)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

message(STATUS "Compiler Id is: ${CMAKE_Fortran_COMPILER_ID}")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

option(SEMBA_FDTD_ENABLE_MPI "Use MPI" OFF)
option(SEMBA_FDTD_ENABLE_OPENMP "Use OpenMP" ON)
option(SEMBA_FDTD_ENABLE_CUDA "Use native CUDA kernels for the Yee field updates" OFF)
set(SEMBA_FDTD_CUDA_ARCHITECTURES "89" CACHE STRING "CUDA architectures to compile for")
option(SEMBA_FDTD_ENABLE_MTLN "Use MTLN" ON)
option(SEMBA_FDTD_ENABLE_SMBJSON "Use smbjson" ON)
option(SEMBA_FDTD_ENABLE_DOUBLE_PRECISION "Use double precision (CompileWithReal8)" OFF)
option(SEMBA_FDTD_ENABLE_GNU_SAFE_MATH "Use GNU Release -O3 without unsafe floating-point optimizations or contraction" OFF)

option(SEMBA_FDTD_ENABLE_TEST "Compile tests" ON)
option(SEMBA_FDTD_OUTPUTS_LIB "Build standalone outputs library" ON)
option(SEMBA_FDTD_MAIN_LIB "Build standalone main library" ON)
option(SEMBA_FDTD_EXECUTABLE "Build standalone executable" ON)

if(SEMBA_FDTD_ENABLE_MPI AND SEMBA_FDTD_ENABLE_TEST)
	set(SEMBA_FDTD_OUTPUT_MPI_TEST_RANKS 3 CACHE STRING
		"Number of ranks used by distributed output integration tests")
endif()

if(SEMBA_FDTD_ENABLE_TEST OR XDMF_HDF5_BUILD_TESTING)
	enable_testing()
endif()

option(SEMBA_FDTD_ENABLE_INTEL_XHOST_OPTIMIZATION "When compiling in Release, enables the -xHost optimization flag (not supported in github actions)" OFF)
option(SEMBA_FDTD_ENABLE_INTEL_IPO "When compiling in Release, enables the interprocedural optimization" OFF)
set(SEMBA_FDTD_EXTRA_FORTRAN_FLAGS "" CACHE STRING "Extra Fortran flags appended to the Windows Intel Release flags")

if(SEMBA_FDTD_ENABLE_CUDA)
	if(NOT SEMBA_FDTD_MAIN_LIB)
		message(FATAL_ERROR "SEMBA_FDTD_ENABLE_CUDA requires SEMBA_FDTD_MAIN_LIB.")
	endif()
	enable_language(CXX)
	enable_language(CUDA)
	find_package(CUDAToolkit REQUIRED)
endif()

if(SEMBA_FDTD_ENABLE_MPI)
	set(HDF5_PREFER_PARALLEL TRUE)
endif()

set(XDMF_HDF5_ENABLE_MPI "${SEMBA_FDTD_ENABLE_MPI}")

# Compilation defines.
if(CMAKE_BUILD_TYPE MATCHES "Release" OR CMAKE_BUILD_TYPE MATCHES "release" )
	add_definitions(-DCompileWithRelease)
else()
	add_definitions(-DCompileWithDebug)
endif()
if (CMAKE_SYSTEM_NAME STREQUAL "Windows")
	add_compile_definitions(__WIN32__)
elseif (CMAKE_SYSTEM_NAME STREQUAL "Darwin")
	add_compile_definitions(__APPLE__)
elseif (CMAKE_SYSTEM_NAME STREQUAL "Linux")
	add_compile_definitions(__linux__)
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
)
if(SEMBA_FDTD_ENABLE_OPENMP)
	add_definitions(-DCompileWithOpenMP)
endif()
if(SEMBA_FDTD_ENABLE_CUDA)
	add_definitions(-DCompileWithCUDA)
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
	add_definitions(-DGNUCompiler)
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM")
	add_definitions(-DIFXCompiler)
endif()

include("${CMAKE_CURRENT_SOURCE_DIR}/set_precompiled_libraries.cmake")
include("${CMAKE_CURRENT_SOURCE_DIR}/cmake/solver_configuration_summary.cmake")

if (CMAKE_SYSTEM_NAME MATCHES "Linux")
	message(STATUS "Using Linux flags")

	if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
		message(STATUS "Using GNU flags")

		set(CMAKE_C_STANDARD 17)
		set(CMAKE_C_STANDARD_REQUIRED ON)
		set(CMAKE_C_EXTENSIONS ON)
		if(SEMBA_FDTD_ENABLE_OPENMP)
			set(CMAKE_CXX_FLAGS "-fopenmp")
			set(CMAKE_Fortran_FLAGS "-fopenmp -ffree-form -ffree-line-length-none -fdec -fallow-argument-mismatch")
		else()
			set(CMAKE_Fortran_FLAGS "-ffree-form -ffree-line-length-none -fdec -fallow-argument-mismatch")
		endif()
		
		set(CMAKE_C_FLAGS_RELEASE "-Ofast")
		set(CMAKE_CXX_FLAGS_RELEASE "-Ofast")
		set(CMAKE_Fortran_FLAGS_RELEASE "-Ofast")
		if(SEMBA_FDTD_ENABLE_GNU_SAFE_MATH)
			set(CMAKE_C_FLAGS_RELEASE "-O3 -fno-fast-math -ffp-contract=off")
			set(CMAKE_CXX_FLAGS_RELEASE "-O3 -fno-fast-math -ffp-contract=off")
			set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -fno-fast-math -ffp-contract=off")
		endif()

		set(CMAKE_C_FLAGS_DEBUG "-g -O0")
		set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
		set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -fno-inline -fcheck=all -fbacktrace")

		if(SEMBA_FDTD_ENABLE_COVERAGE)
			message(STATUS "Enabling coverage instrumentation")
			set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} --coverage")
			set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
			set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} --coverage")
			set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage")
			set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} --coverage")
		endif()

	elseif(CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM")
		message(STATUS "Using IntelLLVM (ifx) flags")

		if(SEMBA_FDTD_ENABLE_OPENMP)
			set(CMAKE_CXX_FLAGS "-qopenmp")
			set(CMAKE_Fortran_FLAGS "-qopenmp -fpp")
		else()
			set(CMAKE_Fortran_FLAGS "-fpp")
		endif()

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
		# MSVC links executables with a 1 MiB stack by default, which is too small
		# for HDF output initialisation with the Intel Fortran runtime.
		add_link_options(/STACK:67108864)
		# HDF output creates size-dependent Fortran temporaries. Keep them off the
		# limited Windows thread stack even when the final executable uses a larger stack.
		add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:/heap-arrays>)

		# This prevents a compilation error happening due to visual studio headers being too strict.
		# It can possibly can be removed in the future.
		#add_compile_definitions("$<$<COMPILE_LANGUAGE:C,CXX>:-D_ALLOW_COMPILER_AND_STL_VERSION_MISMATCH>")
		#add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:/nofpp>)
		#set(CMAKE_Fortran_FLAGS "/Qopenmp")
		#set(CMAKE_CXX_FLAGS "/Qopenmp")
		if(SEMBA_FDTD_ENABLE_OPENMP)
			set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} /Qopenmp")
		endif()

		# Debug flags - enable full debugging support for breakpoints
		set(CMAKE_Fortran_FLAGS_DEBUG "/debug:full /Od /traceback /check:all,nouninit /Zi /debug-parameters:all /warn:interfaces")
		set(CMAKE_CXX_FLAGS_DEBUG "/Zi /Od /RTC1")
		set(CMAKE_C_FLAGS_DEBUG "/Zi /Od /RTC1")
		
		# Release flags
		set(CMAKE_Fortran_FLAGS_RELEASE "/O3")
		set(CMAKE_CXX_FLAGS_RELEASE "/O2")
		set(CMAKE_C_FLAGS_RELEASE "/O2")
		if(SEMBA_FDTD_ENABLE_INTEL_XHOST_OPTIMIZATION)
			set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} /QxHost /fp:fast")
		endif()
		if(SEMBA_FDTD_EXTRA_FORTRAN_FLAGS)
			set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${SEMBA_FDTD_EXTRA_FORTRAN_FLAGS}")
		endif()

	else()
		message(FATAL_ERROR "Unrecognized compiler id:" ${CMAKE_Fortran_COMPILER_ID})
	endif()
else()
	message(FATAL_ERROR "Unrecognized system name")
endif()

# =========================================

add_subdirectory(external)

if(SEMBA_FDTD_ENABLE_CUDA)
	add_library(semba-cuda src_main_pub/cuda_field_updates.cu)
	set_target_properties(semba-cuda PROPERTIES
		CUDA_ARCHITECTURES "${SEMBA_FDTD_CUDA_ARCHITECTURES}"
		CUDA_STANDARD 17
		CUDA_STANDARD_REQUIRED ON
	)
	target_link_libraries(semba-cuda PRIVATE CUDA::cudart)
	if(SEMBA_FDTD_ENABLE_DOUBLE_PRECISION)
		target_compile_definitions(semba-cuda PRIVATE SEMBA_CUDA_DOUBLE_PRECISION)
	endif()
endif()

add_library(semba-types 
	"src_main_pub/nfde_types.F90"
	"src_main_pub/fdetypes.F90"
	"src_mtln/mtln_types.F90"
	"src_conformal/conformal_types.F90"
	"src_wires_pub/wires_types.F90"
	"src_main_pub/lumped_types.F90"
)
target_link_libraries(semba-types ${MPI_Fortran_LIBRARIES})

add_library(semba-reports
	"src_main_pub/errorreport.F90"
	"src_main_pub/snapxdmf.F90"
)
target_link_libraries(semba-reports semba-types XDMF::HDF5)

add_subdirectory(src_utils)

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

add_subdirectory(src_output)
set(OUTPUT_LIBRARIES fdtd-output)
set(VTK_API_LIBRARIES vtkAPI)

add_subdirectory(src_conformal)
set(CONFORMAL_LIBRARIES conformal)

if (SEMBA_FDTD_ENABLE_TEST)
	add_subdirectory(external/googletest/)
	add_subdirectory(test)
endif()



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
if(SEMBA_FDTD_ENABLE_CUDA)
	target_sources(semba-components PRIVATE "src_main_pub/cuda_field_updates.F90")
	target_link_libraries(semba-components semba-cuda CUDA::cudart)
endif()
target_link_libraries(semba-components semba-types semba-reports ${MTLN_LIBRARIES})
if (SEMBA_FDTD_ENABLE_MTLN)
	add_dependencies(semba-components mtlnsolver)
endif()

if(SEMBA_FDTD_OUTPUTS_LIB)
	add_library(semba-outputs
		"src_main_pub/mpicomm.F90"
	)
	target_link_libraries(semba-outputs
		semba-components
		${MPI_Fortran_LIBRARIES})
endif()

if(SEMBA_FDTD_MAIN_LIB)
	# Prepares variables containing information on commit-hash and build flags.
	# This information is only available at the configuration stage.
	set(PROGRAM_NAME "semba-fdtd")
	include("${CMAKE_CURRENT_SOURCE_DIR}/get_commit_info.cmake")
	configure_file(
		${CMAKE_CURRENT_SOURCE_DIR}/src_main_pub/version.F90.in
		${CMAKE_CURRENT_SOURCE_DIR}/src_main_pub/version.F90
		@ONLY
	)

add_library(semba-main
	"src_main_pub/semba_fdtd.F90"
	"src_main_pub/calc_constants.F90"
	"src_main_pub/nfde_rotate.F90"
	"src_main_pub/EpsMuTimeScale.F90"
	"src_main_pub/getargs.F90"
	"src_main_pub/healer.F90"
	"src_main_pub/preprocess_geom.F90"
	"src_main_pub/storegeom.F90"
	"src_main_pub/version.F90"
	"src_main_pub/interpret_switches.F90"
	"src_main_pub/resuming.F90"
	"src_main_pub/timestepping.F90"
)
 target_link_libraries(semba-main
		semba-outputs
		fdtd-utils
		${OUTPUT_LIBRARIES}
		${VTK_API_LIBRARIES}
		${SMBJSON_LIBRARIES}
		${MTLN_LIBRARIES})
endif()

if(SEMBA_FDTD_EXECUTABLE)
	add_executable(semba-fdtd
		"src_main_pub/launcher.F90"
	)
	target_link_libraries(semba-fdtd semba-main semba-reports)
	target_link_libraries(semba-fdtd ${MPI_Fortran_LIBRARIES})
	if(SEMBA_FDTD_EXECUTABLE_NAME)
		set_target_properties(semba-fdtd PROPERTIES OUTPUT_NAME ${SEMBA_FDTD_EXECUTABLE_NAME})
	endif()
endif()



include_directories(${CMAKE_BINARY_DIR}/mod)
include_directories(${HDF5_INCLUDE_DIRS})
include_directories(${FHASH_INCLUDES})

semba_print_configuration_summary()
