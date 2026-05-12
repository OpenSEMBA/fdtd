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
option(SEMBA_FDTD_ENABLE_INTEL_IPO "When compiling in Release, enables the interprocedural optimization" OFF)
option(SEMBA_FDTD_ENABLE_CUDA_FORTRAN "Enable CUDA Fortran GPU acceleration (NVHPC compiler)" OFF)

if(SEMBA_FDTD_ENABLE_CUDA_FORTRAN AND NOT CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC")
	message(FATAL_ERROR "SEMBA_FDTD_ENABLE_CUDA_FORTRAN requires NVHPC Fortran compiler")
endif()

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

if (SEMBA_FDTD_ENABLE_CUDA_FORTRAN)
	add_definitions(-DCompileWithACC)
	add_definitions(-DSEMBA_FDTD_ENABLE_CUDA_FORTRAN)
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
 		
 		# Detect NVHPC SDK location - try standard paths, then NVCOMPILERS env var
 		if(DEFINED ENV{NVCOMPILERS})
 			set(NVHPC_SDK_ROOT "$ENV{NVCOMPILERS}")
 		elseif(EXISTS "/opt/nvidia/hpc_sdk/Linux_x86_64")
 			file(GLOB NVHPC_VERSION_DIRS "/opt/nvidia/hpc_sdk/Linux_x86_64/[0-9]*.[0-9]*")
 			list(SORT NVHPC_VERSION_DIRS)
 			list(GET NVHPC_VERSION_DIRS -1 NVHPC_VERSION_DIR)
 			get_filename_component(NVHPC_SDK_ROOT "${NVHPC_VERSION_DIR}" DIRECTORY)
 		else()
 			set(NVHPC_SDK_ROOT "/opt/nvidia/hpc_sdk/Linux_x86_64/24.7")
 		endif()
 		message(STATUS "NVHPC SDK root: ${NVHPC_SDK_ROOT}")
 		
# Base flags: Fortran free-form preprocessor
		set(CMAKE_Fortran_FLAGS "-Minfo=accel -Mpreprocess -Mbyteswapio")
		# Only set C/CXX accel flags if using NVHPC C/C++ compiler
		if(CMAKE_C_COMPILER_ID MATCHES "NVHPC")
			set(CMAKE_C_FLAGS "-Minfo=accel")
			set(CMAKE_CXX_FLAGS "-Minfo=accel")
		endif()
		
		# Release: aggressive optimization with architecture-specific tuning
		set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -Minfo=accel -Mprefetch -pg")
		if(CMAKE_C_COMPILER_ID MATCHES "NVHPC")
			set(CMAKE_C_FLAGS_RELEASE "-O3 -Minfo=accel")
			set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Minfo=accel")
		endif()
		
		# Debug: no optimization + bounds checking
		set(CMAKE_Fortran_FLAGS_DEBUG "-g -C -Minfo=accel")
		set(CMAKE_C_FLAGS_DEBUG "-g")
		set(CMAKE_CXX_FLAGS_DEBUG "-g")
		
		# CUDA Fortran is enabled per-source (gpu_kernels_cuf.F90) to avoid impacting all Fortran units.
 		
 		# OpenACC requires CUDA runtime libraries at link time
 		set(NVHPC_COMPILER_LIB "${NVHPC_SDK_ROOT}/compilers/lib")
 		# Try CUDA 11.8 first (older SDKs), then 12.x, then 13.x
 		set(CUDA_LIB_DIR "")
 		set(NVVM_LIB_DIR "")
 		foreach(CUDA_VER "11.8" "12.0" "12.8" "13.1")
 			if(EXISTS "${NVHPC_SDK_ROOT}/cuda/${CUDA_VER}/lib64")
 				set(CUDA_LIB_DIR "${NVHPC_SDK_ROOT}/cuda/${CUDA_VER}/lib64")
 				set(NVVM_LIB_DIR "${NVHPC_SDK_ROOT}/cuda/${CUDA_VER}/nvvm/lib64")
 				break()
 			endif()
 		endforeach()
 		
 		if(EXISTS "${CUDA_LIB_DIR}")
 			link_directories("${CUDA_LIB_DIR}" "${NVVM_LIB_DIR}" "${NVHPC_COMPILER_LIB}")
 			list(APPEND CMAKE_EXE_LINKER_FLAGS "-L${CUDA_LIB_DIR} -L${NVVM_LIB_DIR} -L${NVHPC_COMPILER_LIB} -laccstub_static -lnvvm")
 			list(APPEND CMAKE_SHARED_LINKER_FLAGS "-L${CUDA_LIB_DIR} -L${NVVM_LIB_DIR} -L${NVHPC_COMPILER_LIB} -laccstub_static -lnvvm")
 			message(STATUS "NVHPC CUDA lib: ${CUDA_LIB_DIR}")
 			message(STATUS "NVHPC NVVM lib: ${NVVM_LIB_DIR}")
 		else()
 			message(WARNING "NVHPC CUDA toolkit not found at ${NVHPC_SDK_ROOT}/cuda/")
 		endif()
		
		# LAPACK from NVHPC SDK (already configured in set_precompiled_libraries.cmake)
	else()
		message(FATAL_ERROR "Unrecognized compiler:" ${CMAKE_Fortran_COMPILER_ID})
	endif()
elseif(CMAKE_SYSTEM_NAME MATCHES "Windows")
	message(STATUS "Using Windows flags")
	if(CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM")
		message(STATUS "Using IntelLLVM (ifx) flags")

		# This prevents a compilation error happening due to visual studio headers being too strict.
		# It can possibly can be removed in the future.
		#add_compile_definitions("$<$<COMPILE_LANGUAGE:C,CXX>:-D_ALLOW_COMPILER_AND_STL_VERSION_MISMATCH>")
		#add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:/nofpp>)
		#set(CMAKE_Fortran_FLAGS "/Qopenmp")
		#set(CMAKE_CXX_FLAGS "/Qopenmp")

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

add_subdirectory(src_conformal)
set(CONFORMAL_LIBRARIES conformal)

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
	# Prepares variables containing information on commit-hash and build flags.
	# This information is only available at the configuration stage.
	set(PROGRAM_NAME "semba-fdtd")
	include("${CMAKE_CURRENT_SOURCE_DIR}/get_commit_info.cmake")
	configure_file(
	${CMAKE_SOURCE_DIR}/src_main_pub/version.F90.in
	${CMAKE_SOURCE_DIR}/src_main_pub/version.F90
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
		"src_main_pub/postprocess.F90"
		"src_main_pub/interpreta_switches.F90"
		"src_main_pub/resuming.F90"	
		"src_main_pub/timestepping.F90"
	)
	if (SEMBA_FDTD_ENABLE_CUDA_FORTRAN)
		target_sources(semba-main PRIVATE
			"src_main_pub/gpu_core_m.F90"
			"src_main_pub/gpu_yee_m.F90"
			"src_main_pub/gpu_cpml_m.F90"
			"src_main_pub/gpu_sgbc_core_m.F90"
			"src_main_pub/gpu_sgbc_e_m.F90"
			"src_main_pub/gpu_sgbc_h_m.F90"
		)
		set_source_files_properties("src_main_pub/gpu_core_m.F90" "src_main_pub/gpu_yee_m.F90" "src_main_pub/gpu_cpml_m.F90" "src_main_pub/gpu_sgbc_core_m.F90" "src_main_pub/gpu_sgbc_e_m.F90" "src_main_pub/gpu_sgbc_h_m.F90" PROPERTIES COMPILE_OPTIONS "-cuda;-gpu=pinned")
		target_link_options(semba-main PRIVATE -cuda)
	endif()
	target_link_libraries(semba-main
		semba-outputs 
		${SMBJSON_LIBRARIES} 
		${MTLN_LIBRARIES})
endif()

if (SEMBA_FDTD_EXECUTABLE)
	add_executable(semba-fdtd 
		"src_main_pub/launcher.F90"
	)
	target_link_libraries(semba-fdtd semba-main semba-reports)
	target_link_libraries(semba-fdtd ${MPI_Fortran_LIBRARIES})
	
	# Add CUDA Fortran flags to link options
	if (SEMBA_FDTD_ENABLE_CUDA_FORTRAN)
		target_link_options(semba-fdtd PRIVATE -cuda)
	endif()
endif()



include_directories(${CMAKE_BINARY_DIR}/mod)
include_directories(${HDF5_INCLUDE_DIRS})
include_directories(${FHASH_INCLUDES})
