add_library(fdtd-utils
  "utils.F90"
  "valueReplacer.F90"
  "allocationUtils.F90"
)
target_link_libraries(fdtd-utils
 semba-types
)