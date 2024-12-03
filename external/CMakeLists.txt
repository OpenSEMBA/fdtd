
if (SEMBA_FDTD_ENABLE_MTLN OR SEMBA_FDTD_ENABLE_SMBJSON)
	set (FHASH_SOURCES_FOLDER "fhash/src/")
	add_library(fhash 
		${FHASH_SOURCES_FOLDER}/fhash.f90 
		${FHASH_SOURCES_FOLDER}/fhash_sll.f90 
		${FHASH_SOURCES_FOLDER}/fhash_data_container.f90 
		${FHASH_SOURCES_FOLDER}/fhash_tbl.f90 
		${FHASH_SOURCES_FOLDER}/fhash_tbl_iter.f90 
		${FHASH_SOURCES_FOLDER}/fhash_fnv.f90 
		${FHASH_SOURCES_FOLDER}/fhash_key/base.f90 
		${FHASH_SOURCES_FOLDER}/fhash_key/int32.f90 
		${FHASH_SOURCES_FOLDER}/fhash_key/int64.f90 
		${FHASH_SOURCES_FOLDER}/fhash_key/char.f90 
		${FHASH_SOURCES_FOLDER}/fhash_key/int32_1d.f90 
		${FHASH_SOURCES_FOLDER}/fhash_key/int64_1d.f90
	)
endif()
