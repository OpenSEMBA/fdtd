message(STATUS "Creating build system for smbjson")

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

set (JSONFORTRAN_SOURCES_FOLDER "${PROJECT_SOURCE_DIR}/external/json-fortran/src/")
add_library (jsonfortran  
	${JSONFORTRAN_SOURCES_FOLDER}/json_kinds.F90
    ${JSONFORTRAN_SOURCES_FOLDER}/json_parameters.F90
    ${JSONFORTRAN_SOURCES_FOLDER}/json_string_utilities.F90
    ${JSONFORTRAN_SOURCES_FOLDER}/json_value_module.F90
    ${JSONFORTRAN_SOURCES_FOLDER}/json_file_module.F90
    ${JSONFORTRAN_SOURCES_FOLDER}/json_module.F90
)

add_library(smbjson 
	"smbjson_labels.F90"
	"cells.F90"
	"smbjson.F90"
	"idchildtable.F90"
	"mesh.F90"
	"parser_tools.F90"
	"nfdetypes_extension.F90" 
	"geometry.F90"
	"cell_map.F90"
)

target_link_libraries(smbjson PRIVATE jsonfortran semba-types fhash)