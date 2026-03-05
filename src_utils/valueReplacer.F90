module mod_valueReplacer
  use FDETYPES, only: RKIND, CKIND, SINGLE, RKIND_tiempo
  implicit none
  private

  public :: replace_value

  interface replace_value
     ! Scalars
     module procedure replace_scalar_int
     module procedure replace_scalar_real
     module procedure replace_scalar_complex
     
     ! 1D arrays
     module procedure replace_1d_int
     module procedure replace_1d_real
     module procedure replace_1d_complex
     
     ! 2D arrays
     module procedure replace_2d_int
     module procedure replace_2d_real
     module procedure replace_2d_complex
     
     ! 3D arrays
     module procedure replace_3d_int
     module procedure replace_3d_real
     module procedure replace_3d_complex
  end interface

contains
  !=====================
  ! Scalar replacements
  !=====================
  subroutine replace_scalar_int(x, val)
    integer(SINGLE), intent(inout) :: x
    integer(SINGLE), intent(in) :: val
    x = val
  end subroutine

  subroutine replace_scalar_real(x, val)
    real(RKIND), intent(inout) :: x
    real(RKIND), intent(in) :: val
    x = val
  end subroutine

  subroutine replace_scalar_real_t(x, val)
    real(RKIND_tiempo), intent(inout) :: x
    real(RKIND_tiempo), intent(in) :: val
    x = val
  end subroutine

  subroutine replace_scalar_complex(x, val)
    complex(CKIND), intent(inout) :: x
    complex(CKIND), intent(in) :: val
    x = val
  end subroutine

  !=====================
  ! 1D array replacements
  !=====================
  subroutine replace_1d_int(x, idx1, val)
    integer(SINGLE), intent(inout) :: x(:)
    integer(SINGLE), intent(in) :: idx1
    integer(SINGLE), intent(in) :: val
    x(idx1) = val
  end subroutine

  subroutine replace_1d_real(x, idx1, val)
    real(RKIND), intent(inout) :: x(:)
    integer(SINGLE), intent(in) :: idx1
    real(RKIND), intent(in) :: val
    x(idx1) = val
  end subroutine

  subroutine replace_1d_complex(x, idx1, val)
    complex(CKIND), intent(inout) :: x(:)
    integer(SINGLE), intent(in) :: idx1
    complex(CKIND), intent(in) :: val
    x(idx1) = val
  end subroutine

  !=====================
  ! 2D array replacements
  !=====================
  subroutine replace_2d_int(x, idx1, idx2, val)
    integer(SINGLE), intent(inout) :: x(:,:)
    integer(SINGLE), intent(in) :: idx1, idx2
    integer(SINGLE), intent(in) :: val
    x(idx1, idx2) = val
  end subroutine

  subroutine replace_2d_real(x, idx1, idx2, val)
    real(RKIND), intent(inout) :: x(:,:)
    integer(SINGLE), intent(in) :: idx1, idx2
    real(RKIND), intent(in) :: val
    x(idx1, idx2) = val
  end subroutine

  subroutine replace_2d_complex(x, idx1, idx2, val)
    complex(CKIND), intent(inout) :: x(:,:)
    integer(SINGLE), intent(in) :: idx1, idx2
    complex(CKIND), intent(in) :: val
    x(idx1, idx2) = val
  end subroutine

  !=====================
  ! 3D array replacements
  !=====================
  subroutine replace_3d_int(x, idx1, idx2, idx3, val)
    integer(SINGLE), intent(inout) :: x(:,:,:)
    integer(SINGLE), intent(in) :: idx1, idx2, idx3
    integer(SINGLE), intent(in) :: val
    x(idx1, idx2, idx3) = val
  end subroutine

  subroutine replace_3d_real(x, idx1, idx2, idx3, val)
    real(RKIND), intent(inout) :: x(:,:,:)
    integer(SINGLE), intent(in) :: idx1, idx2, idx3
    real(RKIND), intent(in) :: val
    x(idx1, idx2, idx3) = val
  end subroutine

  subroutine replace_3d_complex(x, idx1, idx2, idx3, val)
    complex(CKIND), intent(inout) :: x(:,:,:)
    integer(SINGLE), intent(in) :: idx1, idx2, idx3
    complex(CKIND), intent(in) :: val
    x(idx1, idx2, idx3) = val
  end subroutine

end module mod_valueReplacer
