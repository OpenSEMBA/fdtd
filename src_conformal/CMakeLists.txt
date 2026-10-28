message(STATUS "Creating build system for conformal")

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

add_library(conformal
	"cell_map.F90"
    "conformal_media.F90"
    "conformal_fields.F90"
    "conformal_advance.F90"
    "conformal_sgbc.F90"
    "geometry.F90"
)

target_link_libraries(conformal PRIVATE semba-types semba-reports fhash)
