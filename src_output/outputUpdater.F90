module mod_outputUpdater
  implicit none
  use FDETYPES
contains
  subroutine save_next_scalar(scalar, idx, val)
    real, intent(inout) :: scalar(:)
    integer, intent(in) :: idx
    real, intent(in) :: val
    scalar(idx) = val
  end subroutine save_next_scalar

  subroutine save_next_vector(xVector, yVector, zVector, idx, xVal, yVal, zVal)
    real, intent(inout) :: xVector(:), yVector(:), zVector(:)
    integer, intent(in) :: idx
    real, intent(in) :: xVal, yVal, zVal
    xVector(idx) = xVal
    yVector(idx) = yVal
    zVector(idx) = zVal
  end subroutine save_next_vector

  subroutine add_value(scalar, idx, val)
    complex, intent(inout) :: scalar(:)
    integer, intent(in) :: idx
    complex, intent(in) :: val
    scalar(idx) = val + scalar(idx)
  end subroutine update_scalar_value_freq

  subroutine update_vector_value_freq(xVector, yVector, zVector, idx, xVal, yVal, zVal)
    real, intent(inout) :: xVector(:), yVector(:), zVector(:)
    integer, intent(in) :: idx
    real, intent(in) :: xVal, yVal, zVal
    xVector(idx) = xVal + xVector(idx)
    yVector(idx) = yVal + yVector(idx)
    zVector(idx) = zVal + zVector(idx)
  end subroutine update_vector_value_freq

  subroutine save_scalar_timestep_for_valid_points(scalar, lowerCoord, upperCoord, idx)

end module mod_outputUpdater