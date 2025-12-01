add_library(fdtd-output 
  "output.F90"
  "domain.F90"
  "outputUtils.F90"
  "point_probe_output.F90"
  "wire_current_probe_output.F90"
  "wire_charge_probe_output.F90"
)
target_link_libraries(fdtd-output semba-types )