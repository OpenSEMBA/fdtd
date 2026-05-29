# cpp_build.cmake - Standalone C++ build for semba-fdtd

if(NOT DEFINED CPP_SOURCE_ROOT)
    set(CPP_SOURCE_ROOT ${CMAKE_CURRENT_SOURCE_DIR})
endif()

cmake_minimum_required(VERSION 3.15)
if(NOT DEFINED PROJECT_NAME)
    project(semba-fdtd-cpp CXX)
    enable_language(CXX)
endif()

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

message(STATUS "Compiler Id is: ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

option(SEMBA_FDTD_ENABLE_MPI "Use MPI" OFF)
option(SEMBA_FDTD_ENABLE_HDF "Use HDF" ON)
option(SEMBA_FDTD_ENABLE_MTLN "Use MTLN" ON)
option(SEMBA_FDTD_ENABLE_SMBJSON "Use smbjson" ON)
option(SEMBA_FDTD_ENABLE_DOUBLE_PRECISION "Use double precision (CompileWithReal8)" OFF)
option(SEMBA_FDTD_ENABLE_TEST "Compile tests" ON)
option(SEMBA_FDTD_ENABLE_INTEL_XHOST_OPTIMIZATION "When compiling in Release, enables the -xHost optimization flag" OFF)
option(SEMBA_FDTD_ENABLE_INTEL_IPO "When compiling in Release, enables interprocedural optimization" OFF)
option(SEMBA_FDTD_EXECUTABLE "Compiles executable" ON)
option(SEMBA_FDTD_MAIN_LIB "Compiles main library" ON)
option(SEMBA_FDTD_COMPONENTS_LIB "Compiles components library" ON)
option(SEMBA_FDTD_OUTPUTS_LIB "Compiles outputs library" ON)

if(CMAKE_BUILD_TYPE MATCHES "Release" OR CMAKE_BUILD_TYPE MATCHES "release")
    add_definitions(-DCompileWithRelease)
else()
    add_definitions(-DCompileWithDebug)
endif()

if(SEMBA_FDTD_ENABLE_SMBJSON)
    add_definitions(-DCompileWithSMBJSON)
endif()

if(SEMBA_FDTD_ENABLE_MTLN)
    add_definitions(-DCompileWithMTLN)
endif()

if(SEMBA_FDTD_ENABLE_DOUBLE_PRECISION)
    add_definitions(-DCompileWithReal8)
else()
    add_definitions(-DCompileWithReal4)
endif()

add_definitions(-DCompileWithInt2 -DCompileWithOpenMP)

if(SEMBA_FDTD_ENABLE_MPI)
    find_package(MPI REQUIRED COMPONENTS CXX)
    message(STATUS "MPI found: ${MPI_CXX_LIBRARIES}")
    add_definitions(-DCompileWithMPI)
endif()

if(SEMBA_FDTD_ENABLE_HDF)
    find_package(HDF5 REQUIRED COMPONENTS CXX HL)
    message(STATUS "HDF5 found: ${HDF5_LIBRARIES}")
endif()

find_package(nlohmann_json REQUIRED)
message(STATUS "nlohmann/json found: ${nlohmann_json_VERSION}")

if(SEMBA_FDTD_ENABLE_MTLN)
    find_package(LAPACK REQUIRED)
    find_package(BLAS REQUIRED)
    message(STATUS "LAPACK/BLAS found")
endif()

if(CMAKE_SYSTEM_NAME MATCHES "Linux")
    message(STATUS "Using Linux flags")

    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        message(STATUS "Using GNU flags")
        set(CMAKE_CXX_STANDARD 17)
        set(CMAKE_CXX_FLAGS "-std=c++17 -fopenmp")
        set(CMAKE_CXX_FLAGS_RELEASE "-Ofast")
        set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -fno-inline")
        set(CMAKE_C_FLAGS "-fopenmp")
        set(CMAKE_C_FLAGS_RELEASE "-Ofast")
        set(CMAKE_C_FLAGS_DEBUG "-g -O0")

    elseif(CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
        message(STATUS "Using IntelLLVM (ifx) flags")
        set(CMAKE_CXX_STANDARD 17)
        set(CMAKE_CXX_FLAGS "-std=c++17 -qopenmp")

        if(CMAKE_BUILD_TYPE STREQUAL "Release")
            set(CMAKE_CXX_FLAGS_RELEASE "-O3")

            if(SEMBA_FDTD_ENABLE_INTEL_IPO)
                include(CheckIPOSupported)
                check_ipo_supported(RESULT iporesult)
                if(iporesult)
                    message(STATUS "IPO is supported and being enabled")
                    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ON)
                    set(CMAKE_CXX_FLAGS_RELEASE "-ipo")
                else()
                    message(FATAL_ERROR "IPO is not supported")
                endif()
            endif()

            if(SEMBA_FDTD_ENABLE_INTEL_XHOST_OPTIMIZATION)
                set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -xHost")
            endif()
        endif()

        set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -check all,nouninit -debug full -traceback")
        set(CMAKE_C_FLAGS "-qopenmp")
        set(CMAKE_C_FLAGS_RELEASE "-O3")
        set(CMAKE_C_FLAGS_DEBUG "-g -O0")

    else()
        message(FATAL_ERROR "Unrecognized compiler: ${CMAKE_CXX_COMPILER_ID}")
    endif()

elseif(CMAKE_SYSTEM_NAME MATCHES "Windows")
    message(STATUS "Using Windows flags")
    if(CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
        message(STATUS "Using IntelLLVM (ifx) flags")
        set(CMAKE_CXX_STANDARD 17)
        set(CMAKE_CXX_FLAGS "/Qopenmp /std:c++17")
        set(CMAKE_CXX_FLAGS_DEBUG "/Zi /Od /RTC1")
        set(CMAKE_CXX_FLAGS_RELEASE "/O2")
        set(CMAKE_C_FLAGS "/Qopenmp")
        set(CMAKE_C_FLAGS_DEBUG "/Zi /Od /RTC1")
        set(CMAKE_C_FLAGS_RELEASE "/O2")
    else()
        message(FATAL_ERROR "Unrecognized compiler id: ${CMAKE_CXX_COMPILER_ID}")
    endif()
else()
    message(FATAL_ERROR "Unrecognized system name")
endif()

include_directories(${CPP_SOURCE_ROOT}/src_cpp)
include_directories(${HDF5_INCLUDE_DIRS})

add_library(semba-types INTERFACE)
target_include_directories(semba-types INTERFACE
    ${CPP_SOURCE_ROOT}/src_cpp/main
    ${CPP_SOURCE_ROOT}/src_cpp/mtln
    ${CPP_SOURCE_ROOT}/src_cpp/conformal
    ${CPP_SOURCE_ROOT}/src_cpp/wires
    ${CPP_SOURCE_ROOT}/src_cpp/json_parser
)

add_library(semba-reports
    "${CPP_SOURCE_ROOT}/src_cpp/main/errorreport_core.cpp"
    "${CPP_SOURCE_ROOT}/src_cpp/main/snapxdmf.cpp"
)
target_link_libraries(semba-reports PUBLIC semba-types)
if(SEMBA_FDTD_ENABLE_HDF)
    target_link_libraries(semba-reports PUBLIC ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})
endif()

if(SEMBA_FDTD_ENABLE_SMBJSON)
    add_subdirectory("${CPP_SOURCE_ROOT}/src_cpp/json_parser" "${CMAKE_CURRENT_BINARY_DIR}/smbjson")
    if(PROJECT_IS_TOP_LEVEL)
        set(SMBJSON_LIBRARIES smbjson)
    else()
        set(SMBJSON_LIBRARIES smbjson PARENT_SCOPE)
    endif()
endif()

add_subdirectory("${CPP_SOURCE_ROOT}/external/ngspice" "${CMAKE_CURRENT_BINARY_DIR}/ngspice")

if(SEMBA_FDTD_ENABLE_MTLN)
    set(NGSPICE_LIB ngspice)
    add_definitions(-DCompileWithMTLN)
    add_subdirectory("${CPP_SOURCE_ROOT}/src/mtln/interface" "${CMAKE_CURRENT_BINARY_DIR}/ngspice_interface")
    add_subdirectory("${CPP_SOURCE_ROOT}/src_cpp/mtln" "${CMAKE_CURRENT_BINARY_DIR}/mtln")
    if(PROJECT_IS_TOP_LEVEL)
        set(MTLN_LIBRARIES mtlnsolver)
    else()
        set(MTLN_LIBRARIES mtlnsolver PARENT_SCOPE)
    endif()
endif()

add_subdirectory("${CPP_SOURCE_ROOT}/src_cpp/conformal" "${CMAKE_CURRENT_BINARY_DIR}/conformal")
set(CONFORMAL_LIBRARIES conformal)

if((SEMBA_FDTD_EXECUTABLE OR SEMBA_FDTD_ENABLE_TEST) AND NOT SEMBA_FDTD_MAIN_LIB)
    message(FATAL_ERROR "C++ executable and tests require SEMBA_FDTD_MAIN_LIB=ON")
endif()

if(SEMBA_FDTD_COMPONENTS_LIB)
    add_library(semba-components
        "${CPP_SOURCE_ROOT}/src_cpp/main/anisotropic.cpp"
        "${CPP_SOURCE_ROOT}/src_cpp/main/borderscpml.cpp"
        "${CPP_SOURCE_ROOT}/src_cpp/main/bordersmur.cpp"
        "${CPP_SOURCE_ROOT}/src_cpp/main/bordersother.cpp"
        "${CPP_SOURCE_ROOT}/src_cpp/main/electricdispersive.cpp"
        "${CPP_SOURCE_ROOT}/src_cpp/main/magneticdispersive.cpp"
        "${CPP_SOURCE_ROOT}/src_cpp/main/nodalsources.cpp"
        "${CPP_SOURCE_ROOT}/src_cpp/main/planewaves.cpp"
        "${CPP_SOURCE_ROOT}/src_cpp/main/pml_bodies.cpp"
        "${CPP_SOURCE_ROOT}/src_cpp/main/maloney_nostoch.cpp"
        "${CPP_SOURCE_ROOT}/src_cpp/main/lumped.cpp"
        "${CPP_SOURCE_ROOT}/src_cpp/main/dmma_thin_slot.cpp"
        "${CPP_SOURCE_ROOT}/src_cpp/main/farfield.cpp"
        "${CPP_SOURCE_ROOT}/src_cpp/wires/wires.cpp"
        "${CPP_SOURCE_ROOT}/src_cpp/wires/wires_mtln.cpp"
    )
    target_link_libraries(semba-components PUBLIC semba-types semba-reports ${MTLN_LIBRARIES})
endif()

if(SEMBA_FDTD_OUTPUTS_LIB)
    add_library(semba-outputs
        "${CPP_SOURCE_ROOT}/src_cpp/main/observation_stub.cpp"
        "${CPP_SOURCE_ROOT}/src_cpp/main/vtk_stub.cpp"
        "${CPP_SOURCE_ROOT}/src_cpp/main/xdmf.cpp"
        "${CPP_SOURCE_ROOT}/src_cpp/main/xdmf_h5.cpp"
        "${CPP_SOURCE_ROOT}/src_cpp/main/snapxdmf.cpp"
    )
    if(SEMBA_FDTD_ENABLE_HDF)
        target_compile_definitions(semba-outputs PRIVATE SEMBA_CPP_ENABLE_HDF5)
    endif()
    target_link_libraries(semba-outputs PUBLIC semba-components)
    if(SEMBA_FDTD_ENABLE_HDF)
        target_link_libraries(semba-outputs PUBLIC ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})
    endif()
    if(SEMBA_FDTD_ENABLE_MPI)
        target_link_libraries(semba-outputs PUBLIC ${MPI_CXX_LIBRARIES})
    endif()
endif()

if(SEMBA_FDTD_MAIN_LIB)
    add_library(semba-main
        "${CPP_SOURCE_ROOT}/src_cpp/main/semba_fdtd.cpp"
        "${CPP_SOURCE_ROOT}/src_cpp/main/maloney_nostoch.cpp"
        "${CPP_SOURCE_ROOT}/src_cpp/main/lumped.cpp"
        "${CPP_SOURCE_ROOT}/src_cpp/main/mapvtk_writer.cpp"
    )
    target_include_directories(semba-main PUBLIC
        "${CPP_SOURCE_ROOT}/src_cpp/main"
        "${CPP_SOURCE_ROOT}/external/json/single_include/nlohmann"
    )
    target_link_libraries(semba-main PUBLIC
        semba-reports
        nlohmann_json::nlohmann_json
        ${SMBJSON_LIBRARIES}
    )
    if(SEMBA_FDTD_ENABLE_MPI)
        target_compile_definitions(semba-main PUBLIC CompileWithMPI)
        target_include_directories(semba-main PUBLIC ${MPI_CXX_INCLUDE_DIRS})
        if(TARGET MPI::MPI_CXX)
            target_link_libraries(semba-main PUBLIC MPI::MPI_CXX)
        else()
            target_link_libraries(semba-main PUBLIC ${MPI_CXX_LIBRARIES})
        endif()
    endif()
    if(SEMBA_FDTD_ENABLE_HDF)
        target_sources(semba-main PRIVATE
            "${CPP_SOURCE_ROOT}/src_cpp/main/xdmf_h5.cpp")
        target_compile_definitions(semba-main PRIVATE SEMBA_CPP_ENABLE_HDF5)
        target_link_libraries(semba-main PUBLIC ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})
    endif()
    if(SEMBA_FDTD_ENABLE_MTLN AND SEMBA_FDTD_ENABLE_SMBJSON)
        target_sources(semba-main PRIVATE
            "${CPP_SOURCE_ROOT}/src_cpp/wires/wires_mtln.cpp")
        target_include_directories(semba-main PUBLIC
            "${CPP_SOURCE_ROOT}/src_cpp/mtln"
            "${CPP_SOURCE_ROOT}/src_cpp/wires"
            "${CPP_SOURCE_ROOT}/src_cpp/json_parser"
            "${CPP_SOURCE_ROOT}/external/ngspice/src/include")
        target_link_libraries(semba-main PUBLIC
            smbjson mtlnsolver ngspice_interface ngspice
            ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
    endif()
endif()

if(SEMBA_FDTD_EXECUTABLE)
    add_executable(semba-fdtd-cpp
        "${CPP_SOURCE_ROOT}/src_cpp/main/launcher.cpp"
    )
    target_link_libraries(semba-fdtd-cpp PRIVATE semba-main semba-reports)
    if(SEMBA_FDTD_ENABLE_MPI)
        target_compile_definitions(semba-fdtd-cpp PRIVATE CompileWithMPI)
        if(TARGET MPI::MPI_CXX)
            target_link_libraries(semba-fdtd-cpp PRIVATE MPI::MPI_CXX)
        else()
            target_link_libraries(semba-fdtd-cpp PRIVATE ${MPI_CXX_LIBRARIES})
        endif()
    endif()
endif()

if(SEMBA_FDTD_ENABLE_TEST)
    enable_testing()
    find_package(GTest REQUIRED)
    add_subdirectory("${CPP_SOURCE_ROOT}/external/googletest" "${CMAKE_CURRENT_BINARY_DIR}/googletest")
    add_subdirectory("${CPP_SOURCE_ROOT}/test/cpp" "${CMAKE_CURRENT_BINARY_DIR}/cpp_tests")
endif()
