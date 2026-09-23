# CUDA support for device-resident Yee + CPML (P0/P1).
# Requires CUDA toolkit; first product slice assumes MPI OFF and MTLN OFF.

enable_language(CUDA)

find_package(CUDAToolkit REQUIRED)

# Match RTX 5080 (sm_120) when available; fall back to native.
if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES OR CMAKE_CUDA_ARCHITECTURES STREQUAL "")
	set(CMAKE_CUDA_ARCHITECTURES "90;100;120" CACHE STRING "CUDA archs" FORCE)
endif()

add_library(fdtd-cuda
	fdtd_cuda_context.cu
	fdtd_cuda_yee.cu
	fdtd_cuda_cpml.cu
	fdtd_cuda_probes.cu
	fdtd_cuda_planewave.cu
	fdtd_cuda_api.F90
)
target_include_directories(fdtd-cuda PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(fdtd-cuda PUBLIC CUDA::cudart semba-types stdc++)

set_target_properties(fdtd-cuda PROPERTIES
	CUDA_SEPARABLE_COMPILATION OFF
	CUDA_RESOLVE_DEVICE_SYMBOLS ON
	POSITION_INDEPENDENT_CODE ON
)
target_compile_definitions(fdtd-cuda PUBLIC CompileWithCUDA)
