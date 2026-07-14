find_package(Python3 QUIET COMPONENTS Interpreter)

add_executable(xdmf_hdf5_generate_cases generate_cases.F90)
target_link_libraries(xdmf_hdf5_generate_cases PRIVATE
  OpenSEMBA::xdmf-hdf5
)

set(_xdmf_hdf5_output_dir "${CMAKE_CURRENT_BINARY_DIR}/generated")
file(MAKE_DIRECTORY "${_xdmf_hdf5_output_dir}")

add_test(
  NAME xdmf-hdf5-generate
  COMMAND xdmf_hdf5_generate_cases "${_xdmf_hdf5_output_dir}"
)
if(Python3_Interpreter_FOUND)
  execute_process(
    COMMAND ${Python3_EXECUTABLE} -c "import h5py"
    RESULT_VARIABLE _xdmf_hdf5_h5py_result
    ERROR_QUIET
  )
  if(_xdmf_hdf5_h5py_result EQUAL 0)
    add_test(
      NAME xdmf-hdf5-verify
      COMMAND ${Python3_EXECUTABLE}
        "${CMAKE_CURRENT_SOURCE_DIR}/verify_pairs.py"
        "${_xdmf_hdf5_output_dir}"
    )
    set_tests_properties(xdmf-hdf5-verify PROPERTIES
      DEPENDS xdmf-hdf5-generate
    )
  else()
    message(WARNING
      "h5py was not found; skipping XDMF/HDF5 black-box verification")
  endif()
else()
  message(WARNING
    "Python was not found; skipping XDMF/HDF5 black-box verification")
endif()

find_program(XDMF_HDF5_PVPYTHON_EXECUTABLE pvpython)
if(XDMF_HDF5_PVPYTHON_EXECUTABLE)
  add_test(
    NAME xdmf-hdf5-paraview
    COMMAND ${XDMF_HDF5_PVPYTHON_EXECUTABLE}
      "${CMAKE_CURRENT_SOURCE_DIR}/verify_paraview.py"
      "${_xdmf_hdf5_output_dir}"
  )
  set_tests_properties(xdmf-hdf5-paraview PROPERTIES
    DEPENDS xdmf-hdf5-generate
  )
endif()
