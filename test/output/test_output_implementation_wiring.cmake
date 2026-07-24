file(READ "${PROJECT_SOURCE_DIR}/CMakeLists.txt" top_level_cmake)
file(READ "${PROJECT_SOURCE_DIR}/src_main_pub/timestepping.F90" timestepping)

foreach(legacy_source
    "src_main_pub/observation.F90"
    "src_main_pub/vtk.F90"
    "src_main_pub/xdmf.F90")
  string(FIND "${top_level_cmake}" "\"${legacy_source}\"" legacy_source_index)
  if(NOT legacy_source_index EQUAL -1)
    message(FATAL_ERROR "Legacy output source remains wired: ${legacy_source}")
  endif()
endforeach()

string(FIND "${top_level_cmake}" "-DCompileWithNewOutputModule" output_module_index)
if(output_module_index EQUAL -1)
  message(FATAL_ERROR "The sole configured output implementation is not enabled")
endif()
