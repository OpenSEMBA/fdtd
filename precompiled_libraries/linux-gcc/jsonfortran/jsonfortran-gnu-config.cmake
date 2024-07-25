# config file for the build tree
# Allow other CMake projects to find this one without installing jsonfortran
# (e.g. only configuring and building it)
# No need to use CMakePackageConfigHelpers since we know all the paths with
# certainty in the build tree.

set ( jsonfortran_VERSION 8.3.0 )

# Check that the correct compiler is in use. Mod files and object files/archives
# are NOT compatible across different Fortran compilers when modules are present
set ( jsonfortran-gnu_Fortran_COMPILER_ID GNU )
set ( jsonfortran-gnu_COMPATIBLE_COMPILER TRUE )
if ( NOT ("GNU" MATCHES "${CMAKE_Fortran_COMPILER_ID}") )
  message ( SEND_ERROR "Incompatible Fortran compilers detected! jsonfortran-gnu was compiled with the GNU Fortran compiler, but the current project is trying to use the ${CMAKE_Fortran_COMPILER_ID} Fortran compiler! In general, Fortran modules and libraries can only link against other projects built using the same compiler." )
  set ( jsonfortran-gnu_COMPATIBLE_COMPILER FALSE )
endif ( NOT ("GNU" MATCHES "${CMAKE_Fortran_COMPILER_ID}") )

# Make targets available to be built
include ( "/mnt/c/users/alberto/Work/ugr-fdtd/fdtd/precompiled_libraries/linux-gcc-rls/jsonfortran/jsonfortran-gnu-targets.cmake" )

# Tell the compiler where to find the mod files
set ( jsonfortran_INCLUDE_DIRS "/mnt/c/users/alberto/Work/ugr-fdtd/fdtd/precompiled_libraries/linux-gcc-rls/jsonfortran/include" )
