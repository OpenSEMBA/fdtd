!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  SEMBA_FDTD SOLVER MODULE
!  creation date Date :  June,  6, 2025
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module solver_mod

    !use ...
    implicit none


    type, public :: solver_t
        ! real (kind=rkind), pointer, dimension (:,:,:)  ::  Ex,Ey,Ez,Hx,Hy,Hz
        integer(kind=4) :: numberOfSteps
        private

    contains
        procedure :: run => solver_run
    end type

    interface solver_t
        module procedure solver_ctor 
    end interface

contains

    function solver_ctor() result(res)
        type(solver_t) :: res
    end function

    subroutine solver_run(this)
        class(solver_t) :: this
        integer :: i
        do i = 1, this%numberOfSteps
            ! call this%stepper%step()
            ! call this%step()
        end do

    end subroutine

end module