add_library(fdtd-utils
  "utils.F90"
  "valueReplacer.F90"
  "allocationUtils.F90"
  "directoryUtils.F90"
  "logUtils.F90"
  "systemMemory.F90"
  "systemMemory.c"
  "filesystemAdapter.c"
)
target_link_libraries(fdtd-utils
 semba-types
)
if(WIN32)
  target_link_libraries(fdtd-utils psapi)
endif()
