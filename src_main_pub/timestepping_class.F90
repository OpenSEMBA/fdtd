!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  SEMBA_FDTD TIMESTEPPING MODULE
!  creation date Date :  June,  6, 2025
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module timestepping_mod

!    use fdetypes
!    use report
!    use PostProcessing
!    use Ilumina
!    use Observa
!    use BORDERS_other
!    use Borders_CPML
!    use Borders_MUR
!    use Resuming
!    use nodalsources
!    use Lumped
!    use PMLbodies
!    use xdmf
!    use vtk
! #ifdef CompileWithMPI
!    use MPIcomm
! #endif
! #ifdef CompileWithStochastic
!    use MPI_stochastic
! #endif
! #ifdef CompileWithNIBC
!    use Multiports
! #endif

! #ifdef CompileWithStochastic
!    use sgbc_stoch
! #else
!    use sgbc_NOstoch
! #endif  
!    use EDispersives
!    use MDispersives
!    use Anisotropic
!    use HollandWires     

! #ifdef CompileWithMTLN  
!    use Wire_bundles_mtln_mod             
! #endif       

! #ifdef CompileWithBerengerWires
!    use WiresBerenger
! #ifdef CompileWithMPI
!    use WiresBerenger_MPI
! #endif
! #endif

! #ifdef CompileWithSlantedWires
!    use WiresSlanted
!    use estructura_slanted_m
! #endif


! #ifdef CompileWithConformal
!    USE conformal_time_stepping_m
!    USE CONFORMAL_MAPPED
! #endif
!    USE EpsMuTimeScale_m
!    USE CALC_CONSTANTS
! #ifdef CompileWithPrescale
!    USE P_rescale
! #endif              
! #ifdef CompileWithMTLN
!    ! use mtln_solver_mod, mtln_solver_t => mtln_t
!    use mtln_types_mod, only: mtln_t
!    use Wire_bundles_mtln_mod
! #endif
! !!
! #ifdef CompileWithProfiling
!    use nvtx
! #endif
   implicit none


    type, public :: stepper_t
        real (kind=rkind), pointer, dimension (:,:,:)  ::  Ex,Ey,Ez,Hx,Hy,Hz
        
        private

    contains
        procedure :: preStep => stepper_preStep
        procedure :: step => atepper_step
        procedure :: endStepping => stepper_endStepping
        procedure :: advanceE
        procedure :: advanceH
        procedure :: advanceEx
        procedure :: advanceEy
        procedure :: advanceEz
        procedure :: advanceHx
        procedure :: advanceHy
        procedure :: advanceHz
    end type

    interface stepper_t
        module procedure stepper_ctor 
    end interface
!    private

!    public launch_simulation
! #ifdef CompileWithMTLN
!    public launch_mtln_simulation
! #endif

contains

    function stepper_ctor() result(res)
        type(stepper_t) :: res


    end function

    ! 
    subroutine stepper_preStep(this)
        class(stepper_t) :: this


    end subroutine

    ! this is the time stepping process
    subroutine stepper_step(this)
        class(stepper_t) :: this

        ! flush planewave off
        ! advance anisotropic E
        call this%advanceE()
        ! advance conformal E
        ! freespace advance E
        ! advance wires
        ! advance PML E
        ! advance multiports E
        ! advance sgbc E
        ! advance lumped E
        ! advance dispersives E
        ! advance planewave E
        ! advance nodalE
        ! advance anisotropic H
        call this%advanceH()
        ! freespace advance H
        ! advance magnetic PML
        ! advance sgbc H
        ! advance dispersives H
        ! advance multiports H
        ! advance planewave H
        ! advance nodal H
        ! advance wires H
        ! advance conformal H
        ! advance magnetic MUR
    end subroutine

    ! the following is not part of the timestepping process
    subroutine endStepping(this)
        class(stepper_t) :: this

        ! end confornal simulation
        ! timing
        ! flush observation files
        ! close observation files
        ! postprocess
        ! create vtk
        ! create xdmf
    end subroutine

    subroutine advanceE(this)
        class(stepper_t) :: this
        call this%advanceEx()
        call this%advanceEy()
        call this%advanceEZ()
    end subroutine

    subroutine advanceH(this)
        class(stepper_t) :: this
        call this%advanceHx()
        call this%advanceHy()
        call this%advanceHZ()
    end subroutine

end module stepper_mod