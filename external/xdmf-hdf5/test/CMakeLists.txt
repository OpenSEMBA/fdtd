# Prefer an activated virtual environment, or the parent project's .venv when
# CMake is launched by an IDE. CMake resolves bin/python or Scripts/python.exe.
set(_xdmf_hdf5_venv)
if(DEFINED ENV{VIRTUAL_ENV})
  set(_xdmf_hdf5_venv "$ENV{VIRTUAL_ENV}")
elseif(EXISTS "${CMAKE_SOURCE_DIR}/.venv")
  set(_xdmf_hdf5_venv "${CMAKE_SOURCE_DIR}/.venv")
endif()

if(_xdmf_hdf5_venv)
  get_filename_component(_xdmf_hdf5_venv "${_xdmf_hdf5_venv}" REALPATH)
  if(DEFINED Python3_EXECUTABLE)
    get_filename_component(_xdmf_hdf5_cached_python "${Python3_EXECUTABLE}" REALPATH)
  elseif(DEFINED _Python3_EXECUTABLE)
    get_filename_component(_xdmf_hdf5_cached_python "${_Python3_EXECUTABLE}" REALPATH)
  endif()
  if(DEFINED _xdmf_hdf5_cached_python)
    file(RELATIVE_PATH _xdmf_hdf5_python_relative
      "${_xdmf_hdf5_venv}"
      "${_xdmf_hdf5_cached_python}"
    )
    if(_xdmf_hdf5_python_relative MATCHES "^\\.\\.")
      unset(Python3_EXECUTABLE CACHE)
      unset(Python3_EXECUTABLE)
      unset(_Python3_EXECUTABLE CACHE)
      unset(_Python3_EXECUTABLE)
      unset(_Python3_INTERPRETER_PROPERTIES CACHE)
      unset(_Python3_INTERPRETER_SIGNATURE CACHE)
    endif()
  endif()
  set(Python3_FIND_VIRTUALENV ONLY)
  set(Python3_ROOT_DIR "${_xdmf_hdf5_venv}")
endif()

find_package(Python3 QUIET COMPONENTS Interpreter)

enable_language(C)

add_executable(xdmf_hdf5_generate_cases
  generate_cases.F90
  create_directory.c
)
target_link_libraries(xdmf_hdf5_generate_cases PRIVATE
  XDMF::HDF5
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
    --replace
  DEPENDS xdmf_hdf5_generate_cases
  COMMENT "Refreshing committed XDMF/HDF5 examples"
  VERBATIM
)

add_test(
  NAME xdmf-hdf5-generate
  COMMAND xdmf_hdf5_generate_cases "${_xdmf_hdf5_output_dir}" --replace
)
add_test(
  NAME xdmf-hdf5-generate-help
  COMMAND xdmf_hdf5_generate_cases --help
)
add_test(
  NAME xdmf-hdf5-generate-no-arguments
  COMMAND xdmf_hdf5_generate_cases
)
if(Python3_Interpreter_FOUND)
  execute_process(
    COMMAND ${Python3_EXECUTABLE} -c "import h5py"
    RESULT_VARIABLE _xdmf_hdf5_h5py_result
    ERROR_QUIET
  )
  if(_xdmf_hdf5_h5py_result EQUAL 0)
    message(STATUS "Found h5py with ${Python3_EXECUTABLE}; enabling XDMF/HDF5 black-box verification")
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

if(XDMF_HDF5_PARALLEL_AVAILABLE)
  add_executable(xdmf_hdf5_generate_parallel_case generate_parallel_case.F90)
  target_link_libraries(xdmf_hdf5_generate_parallel_case PRIVATE
    XDMF::HDF5
    MPI::MPI_Fortran
  )

  set(_xdmf_hdf5_mpi_output_dir
    "${CMAKE_CURRENT_BINARY_DIR}/generated-mpi"
  )
  file(MAKE_DIRECTORY "${_xdmf_hdf5_mpi_output_dir}")
  add_test(
    NAME xdmf-hdf5-mpi-generate
    COMMAND
      ${MPIEXEC_EXECUTABLE}
      ${MPIEXEC_NUMPROC_FLAG} 3
      ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:xdmf_hdf5_generate_parallel_case>
      ${MPIEXEC_POSTFLAGS}
      "${_xdmf_hdf5_mpi_output_dir}"
  )
  set_tests_properties(xdmf-hdf5-mpi-generate PROPERTIES
    PROCESSORS 3
    TIMEOUT 30
    LABELS "xdmf-hdf5;mpi;parallel-hdf5"
  )

  if(Python3_Interpreter_FOUND AND _xdmf_hdf5_h5py_result EQUAL 0)
    add_test(
      NAME xdmf-hdf5-mpi-verify
      COMMAND ${Python3_EXECUTABLE}
        "${CMAKE_CURRENT_SOURCE_DIR}/verify_pairs.py"
        "${_xdmf_hdf5_mpi_output_dir}"
        parallel-hyperslab
    )
    set_tests_properties(xdmf-hdf5-mpi-verify PROPERTIES
      DEPENDS xdmf-hdf5-mpi-generate
      LABELS "xdmf-hdf5;mpi;parallel-hdf5"
    )
  endif()
endif()
