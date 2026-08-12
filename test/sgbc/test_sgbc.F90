! Unit tests for src_main_pub/maloney_nostoch.F90
! Tests the Thomas-algorithm tridiagonal solver solve_tridiag_iguales.

! Solve the 1D Poisson tridiagonal system of size 3:
!   [2 -1  0] [x1]   [1]
!   [-1 2 -1] [x2] = [0]
!   [0 -1  2] [x3]   [1]
! Exact solution: x1=x2=x3=1
integer function test_solve_tridiag_3x3_poisson() bind(C, name="test_solve_tridiag_3x3_poisson") result(status)
    use SGBC_nostoch_m, only: solve_tridiag_iguales
    use FDETYPES_m, only: RKIND
    implicit none
    integer, parameter :: n = 3
    real(kind=RKIND) :: d(n), x(n)
    real(kind=RKIND) :: aa, bb, cc       ! interior row values (sub, main, super)
    real(kind=RKIND) :: a1, b1, c1      ! first row
    real(kind=RKIND) :: an, bn, cn      ! last row
    real(kind=RKIND), parameter :: tol = 1.0e-5_RKIND
    integer :: i

    status = 0
    ! Tridiagonal: -1, 2, -1 interior; 0,2,-1 first; -1,2,0 last
    a1 = 0.0_RKIND;  b1 = 2.0_RKIND;  c1 = -1.0_RKIND
    aa = -1.0_RKIND; bb = 2.0_RKIND;  cc = -1.0_RKIND
    an = -1.0_RKIND; bn = 2.0_RKIND;  cn = 0.0_RKIND
    d  = [1.0_RKIND, 0.0_RKIND, 1.0_RKIND]

    call solve_tridiag_iguales(aa, bb, cc, a1, b1, c1, an, bn, cn, d, x, n)

    do i = 1, n
        if (abs(x(i) - 1.0_RKIND) > tol) then
            print *, "test_solve_tridiag_3x3_poisson FAILED: x(", i, ")=", x(i), " expected 1.0"
            status = 1
        end if
    end do
end function test_solve_tridiag_3x3_poisson

! Solve the 1D Poisson tridiagonal system of size 5:
!   Same stencil, d=[1,0,0,0,1]; exact solution is all ones.
integer function test_solve_tridiag_5x5_poisson() bind(C, name="test_solve_tridiag_5x5_poisson") result(status)
    use SGBC_nostoch_m, only: solve_tridiag_iguales
    use FDETYPES_m, only: RKIND
    implicit none
    integer, parameter :: n = 5
    real(kind=RKIND) :: d(n), x(n)
    real(kind=RKIND) :: aa, bb, cc
    real(kind=RKIND) :: a1, b1, c1
    real(kind=RKIND) :: an, bn, cn
    real(kind=RKIND), parameter :: tol = 1.0e-5_RKIND
    integer :: i

    status = 0
    a1 = 0.0_RKIND;  b1 = 2.0_RKIND;  c1 = -1.0_RKIND
    aa = -1.0_RKIND; bb = 2.0_RKIND;  cc = -1.0_RKIND
    an = -1.0_RKIND; bn = 2.0_RKIND;  cn = 0.0_RKIND
    d  = [1.0_RKIND, 0.0_RKIND, 0.0_RKIND, 0.0_RKIND, 1.0_RKIND]

    call solve_tridiag_iguales(aa, bb, cc, a1, b1, c1, an, bn, cn, d, x, n)

    do i = 1, n
        if (abs(x(i) - 1.0_RKIND) > tol) then
            print *, "test_solve_tridiag_5x5_poisson FAILED: x(", i, ")=", x(i), " expected 1.0"
            status = 1
        end if
    end do
end function test_solve_tridiag_5x5_poisson

! Solve a diagonal system of size 4:
!   All off-diagonals zero, main diagonal = 3.0; d=[6,9,-3,12] => x=[2,3,-1,4]
integer function test_solve_tridiag_diagonal_system() bind(C, name="test_solve_tridiag_diagonal_system") result(status)
    use SGBC_nostoch_m, only: solve_tridiag_iguales
    use FDETYPES_m, only: RKIND
    implicit none
    integer, parameter :: n = 4
    real(kind=RKIND) :: d(n), x(n)
    real(kind=RKIND) :: aa, bb, cc
    real(kind=RKIND) :: a1, b1, c1
    real(kind=RKIND) :: an, bn, cn
    real(kind=RKIND), parameter :: tol = 1.0e-5_RKIND
    real(kind=RKIND), dimension(n), parameter :: expected = [2.0_RKIND, 3.0_RKIND, -1.0_RKIND, 4.0_RKIND]
    integer :: i

    status = 0
    a1 = 0.0_RKIND; b1 = 3.0_RKIND; c1 = 0.0_RKIND
    aa = 0.0_RKIND; bb = 3.0_RKIND; cc = 0.0_RKIND
    an = 0.0_RKIND; bn = 3.0_RKIND; cn = 0.0_RKIND
    d  = [6.0_RKIND, 9.0_RKIND, -3.0_RKIND, 12.0_RKIND]

    call solve_tridiag_iguales(aa, bb, cc, a1, b1, c1, an, bn, cn, d, x, n)

    do i = 1, n
        if (abs(x(i) - expected(i)) > tol) then
            print *, "test_solve_tridiag_diagonal_system FAILED: x(", i, ")=", x(i), " expected", expected(i)
            status = 1
        end if
    end do
end function test_solve_tridiag_diagonal_system

