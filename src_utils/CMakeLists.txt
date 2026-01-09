add_library(fdtd-utils
  "valueReplacer.f90"
)
target_link_libraries(fdtd-utils
 semba-types
)