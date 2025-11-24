message(STATUS "Creating build system for conformal")

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

add_library(conformal 
	"cell_map.F90"
    "conformal.F90"
    "geometry.F90"
)

target_link_libraries(conformal PRIVATE semba-types fhash smbjson)