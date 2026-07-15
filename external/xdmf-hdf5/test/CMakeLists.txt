find_package(Python3 QUIET COMPONENTS Interpreter)

add_executable(xdmf_hdf5_generate_cases generate_cases.F90)
target_link_libraries(xdmf_hdf5_generate_cases PRIVATE
  OpenSEMBA::xdmf-hdf5
)

set(_xdmf_hdf5_output_dir "${CMAKE_CURRENT_BINARY_DIR}/generated")
set(_xdmf_hdf5_example_dir
  "${CMAKE_CURRENT_SOURCE_DIR}/../examples/generated"
)
set(_xdmf_hdf5_example_names
  uniform
  volume
  curvilinear
  unstructured
  mixed
  time-series
)
file(MAKE_DIRECTORY "${_xdmf_hdf5_output_dir}")

add_custom_target(xdmf-hdf5-refresh-examples
  COMMAND ${CMAKE_COMMAND} -E make_directory "${_xdmf_hdf5_example_dir}"
  COMMAND xdmf_hdf5_generate_cases
    "${_xdmf_hdf5_example_dir}"
    --examples
  DEPENDS xdmf_hdf5_generate_cases
  COMMENT "Refreshing committed XDMF/HDF5 examples"
  VERBATIM
)

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
    add_test(
      NAME xdmf-hdf5-examples-verify
      COMMAND ${Python3_EXECUTABLE}
        "${CMAKE_CURRENT_SOURCE_DIR}/verify_pairs.py"
        "${_xdmf_hdf5_example_dir}"
        ${_xdmf_hdf5_example_names}
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
  add_test(
    NAME xdmf-hdf5-examples-paraview
    COMMAND ${XDMF_HDF5_PVPYTHON_EXECUTABLE}
      "${CMAKE_CURRENT_SOURCE_DIR}/verify_paraview.py"
      "${_xdmf_hdf5_example_dir}"
  )
endif()
