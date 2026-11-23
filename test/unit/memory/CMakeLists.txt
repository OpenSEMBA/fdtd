message(STATUS "Creating build system for test/unit/memory")

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

add_library(memory_test_fortran
    "test_memory_budget.F90"
)

target_link_libraries(memory_test_fortran
    semba-reports
    semba-components
    conformal
)

add_library(memory_tests "memory_tests.cpp")

target_link_libraries(memory_tests
    memory_test_fortran
    GTest::gtest
)
