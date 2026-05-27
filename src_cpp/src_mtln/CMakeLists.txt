message(STATUS "Creating build system for mtlnsolver (C++)")

add_library(mtlnsolver
    mtl.cpp
    mtl_bundle.cpp
    mtln_solver.cpp
    multipolar_expansion.cpp
    network.cpp
    network_manager.cpp
    utils.cpp
    probes.cpp
    generators.cpp
    dispersive.cpp
    rational_approximation.cpp
    circuit.cpp
    preprocess.cpp
)

target_link_libraries(mtlnsolver PRIVATE
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES}
    ngspice_interface
    semba-types
    ${MPI_CXX_LIBRARIES}
)

target_include_directories(mtlnsolver PRIVATE
    ${CPP_SOURCE_ROOT}/external/ngspice/src/include
)
target_compile_features(mtlnsolver PUBLIC cxx_std_17)
