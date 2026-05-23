message(STATUS "Creating build system for mtlnsolver (C++)")

add_library(mtlnsolver
    "${CMAKE_CURRENT_SOURCE_DIR}/../src_cpp/src_mtln/mtl.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/../src_cpp/src_mtln/mtl_bundle.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/../src_cpp/src_mtln/mtln_solver.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/../src_cpp/src_mtln/multipolar_expansion.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/../src_cpp/src_mtln/network.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/../src_cpp/src_mtln/network_manager.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/../src_cpp/src_mtln/utils.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/../src_cpp/src_mtln/probes.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/../src_cpp/src_mtln/generators.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/../src_cpp/src_mtln/dispersive.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/../src_cpp/src_mtln/rational_approximation.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/../src_cpp/src_mtln/circuit.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/../src_cpp/src_mtln/preprocess.cpp"
)

target_link_libraries(mtlnsolver PRIVATE
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES}
    ngspice
    semba-types
    ${MPI_CXX_LIBRARIES}
)
