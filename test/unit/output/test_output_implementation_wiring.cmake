foreach(legacy_source
    "src_main_pub/observation.F90"
    "src_main_pub/postprocess.F90"
    "src_main_pub/vtk.F90"
    "src_main_pub/xdmf.F90")
  if(EXISTS "${PROJECT_SOURCE_DIR}/${legacy_source}")
    message(FATAL_ERROR "Legacy output source remains: ${legacy_source}")
  endif()
endforeach()

file(READ "${PROJECT_SOURCE_DIR}/CMakeLists.txt" top_level_cmake)
file(READ "${PROJECT_SOURCE_DIR}/src_main_pub/timestepping.F90" timestepping)
string(FIND "${top_level_cmake}${timestepping}" "CompileWithNewOutputModule" output_toggle_index)
if(NOT output_toggle_index EQUAL -1)
  message(FATAL_ERROR "Retired output implementation toggle remains")
endif()
