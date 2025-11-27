add_library(fdtd-output 
  "output.F90"
  "domain.F90"
  "outputUtils.F90"
)
target_link_libraries(fdtd-output semba-types )