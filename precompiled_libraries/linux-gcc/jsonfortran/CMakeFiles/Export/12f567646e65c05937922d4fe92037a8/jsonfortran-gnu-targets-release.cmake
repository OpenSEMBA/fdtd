#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "jsonfortran-gnu::jsonfortran" for configuration "Release"
set_property(TARGET jsonfortran-gnu::jsonfortran APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(jsonfortran-gnu::jsonfortran PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "Fortran"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/jsonfortran-gnu-8.3.0/lib/libjsonfortran.a"
  )

list(APPEND _cmake_import_check_targets jsonfortran-gnu::jsonfortran )
list(APPEND _cmake_import_check_files_for_jsonfortran-gnu::jsonfortran "${_IMPORT_PREFIX}/jsonfortran-gnu-8.3.0/lib/libjsonfortran.a" )

# Import target "jsonfortran-gnu::jsonfortran-static" for configuration "Release"
set_property(TARGET jsonfortran-gnu::jsonfortran-static APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(jsonfortran-gnu::jsonfortran-static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "Fortran"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/jsonfortran-gnu-8.3.0/lib/libjsonfortran.a"
  )

list(APPEND _cmake_import_check_targets jsonfortran-gnu::jsonfortran-static )
list(APPEND _cmake_import_check_files_for_jsonfortran-gnu::jsonfortran-static "${_IMPORT_PREFIX}/jsonfortran-gnu-8.3.0/lib/libjsonfortran.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
