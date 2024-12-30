
message(STATUS "Creating build system for test/hdf")

add_library(
    hdf_test_fortran
    "test_hdf.F90"
)

target_link_libraries(hdf_test_fortran
    ${HDF5_LIBRARIES} 
    ${HDF5_HL_LIBRARIES}
)

add_library(hdf_tests "hdf_tests.cpp")

target_link_libraries(hdf_tests
    hdf_test_fortran
    GTest::gtest
)

include_directories(${HDF5_INCLUDE_DIRS})