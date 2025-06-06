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
        

    contains
        procedure :: init => solver_init
        procedure :: run => solver_run
        procedure :: end => solver_end
    end type

    interface solver_t
        module procedure solver_ctor 
    end interface

contains

    function solver_ctor() result(res)
        type(solver_t) :: res
    end function

    subroutine solver_init(this)
        class(solver_t) :: this
        ! findBounds
        ! readFields
        ! crea_timevector
        ! calc_g1g2gm1gm2
        ! init reporting
        ! init other borders
        ! init PML bodies
        ! init lumped
        ! init wires
        ! init mur abc slanted
        ! init anisotropic
        ! init sgbcs
        ! init multiports
        ! init edispersives
        ! init mdispersives
        ! init planewave
        ! init nodalsources
        ! init hopf
        ! fill mtag
        ! init observation
        ! init timing
    end subroutine

    subroutine solver_run(this)
        class(solver_t) :: this
        integer :: i
        do i = 1, this%numberOfSteps
            ! call step from timeStepping class?
            ! call this%stepper%step()
            ! call this%step()
        end do

    end subroutine

    subroutine solver_end(this)
        class(solver_t) :: this
        ! end confornal simulation
        ! timing
        ! flush observation files
        ! close observation files
        ! postprocess
        ! create vtk
        ! create xdmf
    end subroutine

end module solver_mod