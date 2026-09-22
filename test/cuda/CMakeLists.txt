message(STATUS "Creating build system for test/cuda (SEMBA_FDTD_ENABLE_CUDA)")

# This subdirectory is only added when SEMBA_FDTD_ENABLE_CUDA is ON
# (see test/CMakeLists.txt). Do not add it to CPU or MPI CPU builds.
add_library(cuda_tests "cuda_tests.cpp")

target_compile_definitions(cuda_tests PUBLIC CompileWithCUDA)

target_include_directories(cuda_tests PRIVATE
    ${CMAKE_SOURCE_DIR}/src_cuda
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(cuda_tests
    fdtd-cuda
    GTest::gtest
    stdc++
)
