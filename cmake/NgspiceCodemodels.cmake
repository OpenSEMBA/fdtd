# Builds the ngspice XSPICE code model libraries (*.cm) from the ngspice
# submodule sources, with the same compiler and the same ngspice config.h as
# the ngspice library linked into semba-fdtd.
#
# Code models are loaded at runtime by ngspice (`codemodel <file>.cm` in the
# spinit file) and share C structures with the ngspice core (SPICEdev,
# Mif_Private_t, coreInfo_t, ...). Prebuilt .cm files from another ngspice
# release or built with a different config.h are not ABI compatible and make
# XSPICE devices such as `filesource` fail to parse.
#
# Result: <build>/bin/codemodels/{spinit, *.cm}, with the *.cm also copied to
# the executable folder (<build>/bin).
#
# Usage: include(NgspiceCodemodels) after add_subdirectory(external/ngspice).

set(NGSPICE_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/../external/ngspice")
get_filename_component(NGSPICE_SOURCE_DIR "${NGSPICE_SOURCE_DIR}" ABSOLUTE)
set(NGSPICE_ICM_DIR "${NGSPICE_SOURCE_DIR}/src/xspice/icm")
set(NGSPICE_CMPP_DIR "${NGSPICE_SOURCE_DIR}/src/xspice/cmpp")

set(NGSPICE_CODEMODELS_OUTPUT_DIR "${CMAKE_BINARY_DIR}/bin/codemodels"
	CACHE PATH "Folder where the ngspice code models (*.cm) and spinit are written")
set(NGSPICE_CODEMODELS spice2poly analog digital xtradev xtraevt table)

# --- cmpp: the XSPICE code model preprocessor (host tool).
# The OpenSEMBA ngspice fork removed cmpp/main.c, so a copy of the upstream
# file is kept in this folder.
add_executable(ngspice-cmpp
	"${CMAKE_CURRENT_LIST_DIR}/ngspice_codemodels/cmpp_main.c"
	"${NGSPICE_CMPP_DIR}/file_buffer.c"
	"${NGSPICE_CMPP_DIR}/ifs_lex.c"
	"${NGSPICE_CMPP_DIR}/ifs_yacc.c"
	"${NGSPICE_CMPP_DIR}/mod_lex.c"
	"${NGSPICE_CMPP_DIR}/mod_yacc.c"
	"${NGSPICE_CMPP_DIR}/pp_ifs.c"
	"${NGSPICE_CMPP_DIR}/pp_lst.c"
	"${NGSPICE_CMPP_DIR}/pp_mod.c"
	"${NGSPICE_CMPP_DIR}/read_ifs.c"
	"${NGSPICE_CMPP_DIR}/util.c"
	"${NGSPICE_CMPP_DIR}/writ_ifs.c"
)
target_include_directories(ngspice-cmpp PRIVATE "${NGSPICE_CMPP_DIR}")
if (WIN32)
	target_compile_definitions(ngspice-cmpp PRIVATE YY_NO_UNISTD_H)
endif()

function(_ngspice_read_lst lst_file out_var)
	set(entries "")
	if (EXISTS "${lst_file}")
		file(STRINGS "${lst_file}" lines)
		foreach(line IN LISTS lines)
			string(STRIP "${line}" line)
			if (NOT line STREQUAL "" AND NOT line MATCHES "^#")
				list(APPEND entries "${line}")
			endif()
		endforeach()
	endif()
	set(${out_var} "${entries}" PARENT_SCOPE)
endfunction()

set(_ngspice_codemodel_targets "")
foreach(cm IN LISTS NGSPICE_CODEMODELS)
	set(cm_src "${NGSPICE_ICM_DIR}/${cm}")
	set(cm_gen "${CMAKE_CURRENT_BINARY_DIR}/ngspice_codemodels/${cm}")
	_ngspice_read_lst("${cm_src}/modpath.lst" cm_models)
	_ngspice_read_lst("${cm_src}/udnpath.lst" cm_udns)
	set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
		"${cm_src}/modpath.lst" "${cm_src}/udnpath.lst")

	# Library description headers (cminfo.h, udninfo.h, ...).
	set(cm_descr
		"${cm_gen}/cmextrn.h" "${cm_gen}/cminfo.h" "${cm_gen}/cminfo2.h"
		"${cm_gen}/udnextrn.h" "${cm_gen}/udninfo.h" "${cm_gen}/udninfo2.h"
		"${cm_gen}/objects.inc")
	set(cm_descr_deps "${cm_src}/modpath.lst" "${cm_src}/udnpath.lst")
	foreach(model IN LISTS cm_models)
		list(APPEND cm_descr_deps "${cm_src}/${model}/ifspec.ifs")
	endforeach()
	foreach(udn IN LISTS cm_udns)
		list(APPEND cm_descr_deps "${cm_src}/${udn}/udnfunc.c")
	endforeach()
	add_custom_command(
		OUTPUT ${cm_descr}
		COMMAND "${CMAKE_COMMAND}" -E make_directory "${cm_gen}"
		COMMAND "${CMAKE_COMMAND}" -E env "CMPP_IDIR=${cm_src}" "CMPP_ODIR=${cm_gen}"
			$<TARGET_FILE:ngspice-cmpp> -lst
		DEPENDS ngspice-cmpp ${cm_descr_deps}
		COMMENT "cmpp: generating ${cm} code model library description"
		VERBATIM)

	set(cm_sources
		"${NGSPICE_ICM_DIR}/dlmain.c"
		"${NGSPICE_SOURCE_DIR}/src/misc/dstring.c"
		${cm_descr})

	foreach(model IN LISTS cm_models)
		set(model_src "${cm_src}/${model}")
		set(model_gen "${cm_gen}/${model}")
		add_custom_command(
			OUTPUT "${model_gen}/ifspec.c" "${model_gen}/cfunc.c"
			COMMAND "${CMAKE_COMMAND}" -E make_directory "${model_gen}"
			COMMAND "${CMAKE_COMMAND}" -E env "CMPP_IDIR=${model_src}" "CMPP_ODIR=${model_gen}"
				$<TARGET_FILE:ngspice-cmpp> -ifs
			COMMAND "${CMAKE_COMMAND}" -E env "CMPP_IDIR=${model_src}" "CMPP_ODIR=${model_gen}"
				$<TARGET_FILE:ngspice-cmpp> -mod
			DEPENDS ngspice-cmpp "${model_src}/ifspec.ifs" "${model_src}/cfunc.mod"
			COMMENT "cmpp: generating ${cm}/${model}"
			VERBATIM)
		list(APPEND cm_sources "${model_gen}/ifspec.c" "${model_gen}/cfunc.c")
		# cfunc.c includes headers and sources that live next to cfunc.mod.
		set_source_files_properties("${model_gen}/ifspec.c" "${model_gen}/cfunc.c"
			PROPERTIES INCLUDE_DIRECTORIES "${model_src}")
	endforeach()

	foreach(udn IN LISTS cm_udns)
		list(APPEND cm_sources "${cm_src}/${udn}/udnfunc.c")
	endforeach()

	set(target "ngspice-codemodel-${cm}")
	add_library(${target} MODULE ${cm_sources})
	target_include_directories(${target} PRIVATE
		"${cm_gen}"
		"${cm_src}"
		"${NGSPICE_SOURCE_DIR}/src/include")
	if (WIN32)
		target_compile_definitions(${target} PRIVATE)
	else()
		target_link_libraries(${target} PRIVATE m)
	endif()
	set_target_properties(${target} PROPERTIES
		PREFIX ""
		SUFFIX ".cm"
		OUTPUT_NAME "${cm}"
		POSITION_INDEPENDENT_CODE ON
		LIBRARY_OUTPUT_DIRECTORY "${NGSPICE_CODEMODELS_OUTPUT_DIR}"
		RUNTIME_OUTPUT_DIRECTORY "${NGSPICE_CODEMODELS_OUTPUT_DIR}")
	# Multi-config generators append the configuration name otherwise.
	foreach(config IN ITEMS DEBUG RELEASE RELWITHDEBINFO MINSIZEREL)
		set_target_properties(${target} PROPERTIES
			LIBRARY_OUTPUT_DIRECTORY_${config} "${NGSPICE_CODEMODELS_OUTPUT_DIR}"
			RUNTIME_OUTPUT_DIRECTORY_${config} "${NGSPICE_CODEMODELS_OUTPUT_DIR}")
	endforeach()
	list(APPEND _ngspice_codemodel_targets ${target})
endforeach()

# spinit that loads the code models from its own folder.
configure_file("${CMAKE_CURRENT_LIST_DIR}/ngspice_codemodels/spinit"
	"${NGSPICE_CODEMODELS_OUTPUT_DIR}/spinit" COPYONLY)

# semba-fdtd loads the code models found next to its executable when no
# spinit loaded them, and the release archives ship them from there. Copy
# them on every build so stale, ABI incompatible ones are always replaced.
set(_ngspice_codemodel_files "")
foreach(target IN LISTS _ngspice_codemodel_targets)
	list(APPEND _ngspice_codemodel_files "$<TARGET_FILE:${target}>")
endforeach()
add_custom_target(ngspice-codemodels ALL
	COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}"
	COMMAND ${CMAKE_COMMAND} -E copy_if_different
		${_ngspice_codemodel_files} "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}"
	DEPENDS ${_ngspice_codemodel_targets}
	COMMENT "Copying ngspice code models next to semba-fdtd")
