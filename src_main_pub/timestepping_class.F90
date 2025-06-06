!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  SEMBA_FDTD TIMESTEPPING MODULE
!  creation date Date :  June,  6, 2025
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module timestepping_mod

    use fdetypes
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
        real (kind=rkind), pointer, dimension (:) ::  g1,g2,gM1,gM2
        real (kind=rkind), pointer, dimension (:) :: Idxh, Idyh, Idzh
        type (bounds_t) :: bounds

    contains
        procedure :: step => stepper_step
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
        ! bounds are defined in solver%findbounds
        ! no need to pass sgg

        ! type(sggfdtdinfo) :: sgg
        ! allocate(this%Idxh(sgg%alloc(iEx)%XI:sgg%alloc(iEx)%XE), &
        !          this%Idyh(sgg%alloc(iEy)%YI:sgg%alloc(iEy)%YE), &
        !          this%Idzh(sgg%alloc(iEz)%ZI:sgg%alloc(iEz)%ZE))
        allocate(this%Idxh(this%bounds%dxh%XI:this%bounds%dxh%XE), &
                 this%Idyh(this%bounds%dyh%YI:this%bounds%dyh%YE), &
                 this%Idzh(this%bounds%dzh%ZI:this%bounds%dzh%ZE))


    end function


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

    ! index shifting is needed for MPI
    ! in each slice arrays start on 1
    ! global indexing is lost
    subroutine advanceEx(this, Ex, Hy, Hz)
        class(stepper_t) :: this
        integer(kind = INTEGERSIZEOFMEDIAMATRICES) :: medio
        integer(kind=4) :: i, j, k
        real (kind=rkind), dimension(0:this%bounds%Ex%NX-1 , 0:this%bounds%Ex%NY-1 , 0:this%bounds%Ex%NZ-1), intent(inout) :: Ex
        real (kind=rkind), dimension(0:this%bounds%Hy%NX-1 , 0:this%bounds%Hy%NY-1 , 0:this%bounds%Hy%NZ-1), intent(inout) :: Hy
        real (kind=rkind), dimension(0:this%bounds%Hz%NX-1 , 0:this%bounds%Hz%NY-1 , 0:this%bounds%Hz%NZ-1), intent(inout) :: Hz
        real (kind=rkind), dimension(:), pointer  ::  Idyh
        real (kind=rkind), dimension(:), pointer  ::  Idzh
        real (kind=INTEGERSIZEOFMEDIAMATRICES), dimension(:,:,:), pointer  ::  sggmiEx
        
        allocate(Idyh(0:this%bounds%dyh%NY-1))
        allocate(Idzh(0:this%bounds%dzh%NZ-1))
        allocate(sggmiEx(0:this%bounds%sggMiEx%NX-1,0:this%bounds%sggMiEx%NY-1,0:this%bounds%sggMiEx%NZ-1 ))
        Idhy => this%Idhy
        Idhz => this%Idhz
        sggmiEx => this%sggmiEx
        do k=1,this%bounds%sweepEx%NZ
            do j=1,this%bounds%sweepEx%NY
                do i=1,this%bounds%sweepEx%NX
                    medio = sggMiEx(i,j,k)
                    Ex(i,j,k) = this%g1(medio)*Ex(i,j,k) + & 
                                this%g2(medio)*((Hz(i,j,k)-Hz(i,j-1,k))*Idyh(j) - & 
                                                (Hy(i,j,k)-Hy(i,j,k-1))*Idzh(k))
                end do
            end do
        end do

    end subroutine


end module timestepping_mod