# Install script for directory: /mnt/c/users/alberto/Work/ugr-fdtd/fdtd/external/json-fortran

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/jsonfortran-gnu-8.3.0/lib" TYPE STATIC_LIBRARY FILES "/mnt/c/users/alberto/Work/ugr-fdtd/fdtd/precompiled_libraries/linux-gcc-rls/jsonfortran/libjsonfortran.a")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/jsonfortran-gnu-8.3.0/lib" TYPE STATIC_LIBRARY FILES "/mnt/c/users/alberto/Work/ugr-fdtd/fdtd/precompiled_libraries/linux-gcc-rls/jsonfortran/lib/libjsonfortran.a")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(GLOB_RECURSE MODULE_FILES "/mnt/c/users/alberto/Work/ugr-fdtd/fdtd/precompiled_libraries/linux-gcc-rls/jsonfortran/include/*.mod")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(GLOB_RECURSE SUBMOD_FILES "/mnt/c/users/alberto/Work/ugr-fdtd/fdtd/precompiled_libraries/linux-gcc-rls/jsonfortran/include/*.smod")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL ${MODULE_FILES} DESTINATION "${CMAKE_INSTALL_PREFIX}/jsonfortran-gnu-8.3.0/lib")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL ${SUBMOD_FILES} DESTINATION "${CMAKE_INSTALL_PREFIX}/jsonfortran-gnu-8.3.0/lib")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/jsonfortran-gnu-8.3.0/cmake/jsonfortran-gnu-targets.cmake")
    file(DIFFERENT _cmake_export_file_changed FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/jsonfortran-gnu-8.3.0/cmake/jsonfortran-gnu-targets.cmake"
         "/mnt/c/users/alberto/Work/ugr-fdtd/fdtd/precompiled_libraries/linux-gcc-rls/jsonfortran/CMakeFiles/Export/12f567646e65c05937922d4fe92037a8/jsonfortran-gnu-targets.cmake")
    if(_cmake_export_file_changed)
      file(GLOB _cmake_old_config_files "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/jsonfortran-gnu-8.3.0/cmake/jsonfortran-gnu-targets-*.cmake")
      if(_cmake_old_config_files)
        string(REPLACE ";" ", " _cmake_old_config_files_text "${_cmake_old_config_files}")
        message(STATUS "Old export file \"$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/jsonfortran-gnu-8.3.0/cmake/jsonfortran-gnu-targets.cmake\" will be replaced.  Removing files [${_cmake_old_config_files_text}].")
        unset(_cmake_old_config_files_text)
        file(REMOVE ${_cmake_old_config_files})
      endif()
      unset(_cmake_old_config_files)
    endif()
    unset(_cmake_export_file_changed)
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/jsonfortran-gnu-8.3.0/cmake" TYPE FILE FILES "/mnt/c/users/alberto/Work/ugr-fdtd/fdtd/precompiled_libraries/linux-gcc-rls/jsonfortran/CMakeFiles/Export/12f567646e65c05937922d4fe92037a8/jsonfortran-gnu-targets.cmake")
  if(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/jsonfortran-gnu-8.3.0/cmake" TYPE FILE FILES "/mnt/c/users/alberto/Work/ugr-fdtd/fdtd/precompiled_libraries/linux-gcc-rls/jsonfortran/CMakeFiles/Export/12f567646e65c05937922d4fe92037a8/jsonfortran-gnu-targets-release.cmake")
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/jsonfortran-gnu-8.3.0/cmake" TYPE FILE FILES
    "/mnt/c/users/alberto/Work/ugr-fdtd/fdtd/precompiled_libraries/linux-gcc-rls/jsonfortran/pkg/jsonfortran-gnu-config.cmake"
    "/mnt/c/users/alberto/Work/ugr-fdtd/fdtd/precompiled_libraries/linux-gcc-rls/jsonfortran/jsonfortran-gnu-config-version.cmake"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/jsonfortran-gnu-8.3.0/lib/pkgconfig" TYPE FILE FILES "/mnt/c/users/alberto/Work/ugr-fdtd/fdtd/precompiled_libraries/linux-gcc-rls/jsonfortran/json-fortran.pc")
endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/mnt/c/users/alberto/Work/ugr-fdtd/fdtd/precompiled_libraries/linux-gcc-rls/jsonfortran/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
