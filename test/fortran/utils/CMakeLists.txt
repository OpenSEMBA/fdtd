message(STATUS "Creating build system for test/observation")

add_library(
    test_utils_fortran
    "fdetypes_tools.F90"
)

target_link_libraries(test_utils_fortran
    semba-types
)
