# CUDA support for device-resident Yee + CPML (P0/P1).
# Toolkit, host compiler, and architecture are selected by cmake/semba_cuda.cmake.

add_library(fdtd-cuda
	fdtd_cuda_context.cu
	fdtd_cuda_yee.cu
	fdtd_cuda_cpml.cu
	fdtd_cuda_probes.cu
	fdtd_cuda_planewave.cu
	fdtd_cuda_mur.cu
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
