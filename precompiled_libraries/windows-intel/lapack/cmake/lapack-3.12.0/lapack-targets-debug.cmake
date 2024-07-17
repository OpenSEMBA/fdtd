#----------------------------------------------------------------
# Generated CMake target import file for configuration "Debug".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "blas" for configuration "Debug"
set_property(TARGET blas APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(blas PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "Fortran"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib/libblas.lib"
  )

list(APPEND _cmake_import_check_targets blas )
list(APPEND _cmake_import_check_files_for_blas "${_IMPORT_PREFIX}/lib/libblas.lib" )

# Import target "lapack" for configuration "Debug"
set_property(TARGET lapack APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(lapack PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "Fortran"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib/liblapack.lib"
  )

list(APPEND _cmake_import_check_targets lapack )
list(APPEND _cmake_import_check_files_for_lapack "${_IMPORT_PREFIX}/lib/liblapack.lib" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
