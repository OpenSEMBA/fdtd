# Config file for the INSTALLED package
# Allow other CMake projects to find this package if it is installed
# Requires the use of the standard CMake module CMakePackageConfigHelpers

set ( jsonfortran_VERSION 8.3.0 )

# Check that the correct compiler is in use. Mod files and object files/archives
# are NOT compatible across different Fortran compilers when modules are present
set ( jsonfortran-gnu_Fortran_COMPILER_ID GNU )
set ( jsonfortran-gnu_COMPATIBLE_COMPILER TRUE )
if ( NOT ("GNU" MATCHES "${CMAKE_Fortran_COMPILER_ID}") )
  message ( SEND_ERROR "Incompatible Fortran compilers detected! jsonfortran-gnu was compiled with the GNU Fortran compiler, but the current project is trying to use the ${CMAKE_Fortran_COMPILER_ID} Fortran compiler! In general, Fortran modules and libraries can only link against other projects built using the same compiler." )
  set ( jsonfortran-gnu_COMPATIBLE_COMPILER FALSE )
endif ( NOT ("GNU" MATCHES "${CMAKE_Fortran_COMPILER_ID}") )


####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was jsonfortran-config.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../" ABSOLUTE)

macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

macro(check_required_components _NAME)
  foreach(comp ${${_NAME}_FIND_COMPONENTS})
    if(NOT ${_NAME}_${comp}_FOUND)
      if(${_NAME}_FIND_REQUIRED_${comp})
        set(${_NAME}_FOUND FALSE)
      endif()
    endif()
  endforeach()
endmacro()

####################################################################################

# Provide the targets
set_and_check ( jsonfortran-gnu_CONFIG_INSTALL_DIR "${PACKAGE_PREFIX_DIR}/jsonfortran-gnu-8.3.0/cmake" )
include ( "${jsonfortran-gnu_CONFIG_INSTALL_DIR}/jsonfortran-gnu-targets.cmake" )

# Make the module files available via include
set_and_check ( jsonfortran_INCLUDE_DIRS "${PACKAGE_PREFIX_DIR}/jsonfortran-gnu-8.3.0/lib" )
