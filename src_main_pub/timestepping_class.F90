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
        procedure :: advanceFreeSpaceE
        procedure :: advanceWires
        procedure :: advancePlaneWaveE
        
        procedure :: FreeSpace_Advance_Ex
        procedure :: FreeSpace_Advance_Ey
        procedure :: FreeSpace_Advance_Ez
        procedure :: FreeSpace_Advance_Hx
        procedure :: FreeSpace_Advance_Hy
        procedure :: FreeSpace_Advance_Hz
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
        if (Thereare%Anisotropic) call advanceAnisotropicE(sgg%alloc,ex,ey,ez,hx,hy,hz,Idxe,Idye,Idze,Idxh,Idyh,Idzh)
        call this%advanceE()
        if (planewavecorr) call this%advanceFreeSpaceE()
#ifdef CompileWithConformal
        if(input_conformal_flag) call conformal_advance_E()
#endif
        call this%advanceWires()
        call this%advancePMLE()
#ifdef CompileWithNIBC
        IF (Thereare%Multiports.and.(mibc)) call AdvanceMultiportE(sgg%alloc,Ex, Ey, Ez)
#endif
        IF (Thereare%sgbcs.and.(sgbc))  then
            call AdvancesgbcE(real(sgg%dt,RKIND),sgbcDispersive,simu_devia,stochastic)
        endif
        if (ThereAre%Lumpeds) call AdvanceLumpedE(sgg,n,simu_devia,stochastic)
        if (Thereare%Edispersives) call AdvanceEDispersiveE(sgg)
        call this%advancePlaneWaveE()
        If (Thereare%NodalE) call AdvanceNodalE(sgg,sggMiEx,sggMiEy,sggMiEz,sgg%NumMedia,n, b,G2,Idxh,Idyh,Idzh,Ex,Ey,Ez,simu_devia)

#ifdef CompileWithMPI
        if (size>1) then
            call MPI_Barrier(SUBCOMM_MPI,ierr)
            call   FlushMPI_E_Cray
        endif
#endif


        IF (Thereare%Anisotropic) call AdvanceAnisotropicH(sgg%alloc,ex,ey,ez,hx,hy,hz,Idxe,Idye,Idze,Idxh,Idyh,Idzh)
        call this%advanceH()
        if (planewavecorr) call this%advanceFreeSpaceH()
        call this%advancePMLH()
        call this%advancePMC()
        if (Thereare%sgbcs.and.(sgbc)) call AdvancesgbcH()
        if (Thereare%Mdispersives) call AdvanceMDispersiveH(sgg)
#ifdef CompileWithNIBC
         !Multiports H-field advancing
         IF (Thereare%Multiports    .and.(mibc))  &
         call AdvanceMultiportH    (sgg%alloc,Hx,Hy,Hz,Ex,Ey,Ez,Idxe,Idye,Idze,sggMiHx,sggMiHy,sggMiHz,gm2,sgg%nummedia,conformalskin)
#endif
        call this%AdvancePlaneWaveH()
        If (Thereare%NodalH) call AdvanceNodalH(sgg,sggMiHx,sggMiHy,sggMiHz,sgg%NumMedia,n, b ,GM2,Idxe,Idye,Idze,Hx,Hy,Hz,simu_devia)

        ! advance wires H: no hace nada salvo thick wires, y no implementado
        call this%advancePMC()

#ifdef CompileWithConformal                      
        if(input_conformal_flag) call conformal_advance_H()
#endif
        ! advance magnetic MUR
    end subroutine


    subroutine advanceE(this, Ex, Ey, Ez, Hx, Hy, Hz)
        class(stepper_t) :: this
        call this%advanceEx(Ex, Hy, Hz)
        call this%advanceEy(Ey, Hz, Hx)
        call this%advanceEZ(Ez, Hx, Hy)
    end subroutine

    subroutine advanceH(this, Ex, Ey, Ez, Hx, Hy, Hz)
        class(stepper_t) :: this
        call this%advanceHx(Hx, Ey, Ez)
        call this%advanceHy(Hy, Ez, Ex)
        call this%advanceHZ(Hz, Ex, Ey)
    end subroutine

    subroutine advanceFreeSpaceE(this)
        class(stepper_t) :: this
        call FreeSpace_Advance_Ex()
        call FreeSpace_Advance_Ey()
        call FreeSpace_Advance_Ez()
    end subroutine

    subroutine advanceFreeSpaceH(this)
        class(stepper_t) :: this
        call FreeSpace_Advance_Hx()
        call FreeSpace_Advance_Hy()
        call FreeSpace_Advance_Hz()
    end subroutine

    subroutine advanceWires(this)
        class(stepper_t) :: this

        if (((trim(adjustl(wiresflavor))=='holland') .or. &
              (trim(adjustl(wiresflavor))=='transition')) .and. .not. use_mtln_wires) then
            if (Thereare%Wires) then
                if (wirecrank) then
                    call AdvanceWiresEcrank(sgg,n, layoutnumber,wiresflavor,simu_devia,stochastic)
                else
#ifdef CompileWithMTLN
                    if (mtln_parsed%has_multiwires) then
                        write(buff, *) 'ERROR: Multiwires in simulation but -mtlnwires flag has not been selected'
                        call WarnErrReport(buff)
                    end if
#endif
                    call AdvanceWiresE(sgg,n, layoutnumber,wiresflavor,simu_devia,stochastic,experimentalVideal,wirethickness,eps0,mu0)                 
                endif
            endif
        endif
#ifdef CompileWithBerengerWires
        if (trim(adjustl(wiresflavor))=='berenger' .and. Thereare%Wires) then
            call AdvanceWiresE_Berenger(sgg,n)
        endif
#endif
#ifdef CompileWithSlantedWires
        if((trim(adjustl(wiresflavor))=='slanted').or.(trim(adjustl(wiresflavor))=='semistructured')) then
            call AdvanceWiresE_Slanted(sgg,n) 
        endif
#endif
        if (use_mtln_wires) then
#ifdef CompileWithMTLN
            call AdvanceWiresE_mtln(sgg,Idxh,Idyh,Idzh,eps0,mu0)
#else
            write(buff,'(a)') 'WIR_ERROR: Executable was not compiled with MTLN modules.'
#endif   
        endif
    end subroutine

    subroutine advancePMLE(this)
        class(stepper_t) :: this
        If (Thereare%PMLbodies) call AdvancePMLbodyE()
        If (Thereare%PMLBorders) then
            call AdvanceelectricCPML(sgg%NumMedia, b       ,sggMiEx,sggMiEy,sggMiEz,G2,Ex,Ey,Ez,Hx,Hy,Hz)
        endif
        if (planewavecorr) then
            If (Thereare%PMLBorders) then
                call AdvanceelectricCPML_freespace (sgg%NumMedia, b       ,sggMiEx,sggMiEy,sggMiEz,G2,Exvac,Eyvac,Ezvac,Hxvac,Hyvac,Hzvac)
            endif
        endif
    end subroutine

    subroutine advancePMLH(this)
        class(stepper_t) :: this
        if (Thereare%PMLbodies) call AdvancePMLbodyH
        If (Thereare%PMLBorders) then
            call AdvanceMagneticCPML(sgg%NumMedia, b, sggMiHx, sggMiHy, sggMiHz, gm2, Hx, Hy, Hz, Ex, Ey, Ez)
        endif
        if (planewavecorr) then
            If (Thereare%PMLBorders) then
                if (sgg%therearePMLMagneticMedia) then
                    call AdvanceMagneticCPML_freespace(sgg%NumMedia, b, sggMiHx, sggMiHy, sggMiHz, gm2, Hxvac, Hyvac, Hzvac, Exvac, Eyvac, Ezvac)
                else
                    continue
                endif
            endif
        endif
    end subroutine

    subroutine advancePMC(this)
        class(stepper_t) :: this
        if (Thereare%PMCBorders) call MinusCloneMagneticPMC(sgg%alloc,sgg%Border,Hx,Hy,Hz,sgg%sweep,layoutnumber,size)
        If (Thereare%PeriodicBorders) call CloneMagneticPeriodic(sgg%alloc,sgg%Border,Hx,Hy,Hz,sgg%sweep,layoutnumber,size)
    end subroutine

    subroutine advancePlaneWaveE(this)
        class(stepper_t) :: this
        !check logic is the same
        if (Thereare%PlaneWaveBoxes.and.still_planewave_time .and. .not.simu_devia) then
            call AdvancePlaneWaveE(sgg,n,b,G2,Idxh,Idyh,Idzh,Ex,Ey,Ez,still_planewave_time)
            if (planewavecorr) call AdvancePlaneWaveE(sgg,n, b,G2,Idxh,Idyh,Idzh,Exvac,Eyvac,Ezvac,still_planewave_time)
        endif
    end subroutine

    subroutine advancePlaneWaveH(this)
        class(stepper_t) :: this
        If (Thereare%PlaneWaveBoxes.and.still_planewave_time .and. .not.simu_devia)  then
            call AdvancePlaneWaveH(sgg,n, b, GM2, Idxe,Idye, Idze, Hx, Hy, Hz,still_planewave_time)
            if (planewavecorr) call AdvancePlaneWaveH(sgg,n, b , GM2, Idxe,Idye, Idze, Hxvac, Hyvac, Hzvac,still_planewave_time)
        endif
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

#ifdef CompileWithProfiling
        call nvtxStartRange("Antes del bucle EX")
#endif        
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
#ifdef CompileWithProfiling
        call nvtxEndRange
#endif

    end subroutine

    subroutine advanceEy(this, Ey, Hz, Hx)
        class(stepper_t) :: this
        integer(kind = INTEGERSIZEOFMEDIAMATRICES) :: medio
        integer(kind=4) :: i, j, k
        real (kind=rkind), dimension(0:this%bounds%Ey%NX-1 , 0:this%bounds%Ey%NY-1 , 0:this%bounds%Ey%NZ-1), intent(inout) :: Ey
        real (kind=rkind), dimension(0:this%bounds%Hz%NX-1 , 0:this%bounds%Hz%NY-1 , 0:this%bounds%Hz%NZ-1), intent(inout) :: Hz
        real (kind=rkind), dimension(0:this%bounds%Hx%NX-1 , 0:this%bounds%Hx%NY-1 , 0:this%bounds%Hx%NZ-1), intent(inout) :: Hx
        real (kind=rkind), dimension(:), pointer  ::  Idzh
        real (kind=rkind), dimension(:), pointer  ::  Idxh
        real (kind=INTEGERSIZEOFMEDIAMATRICES), dimension(:,:,:), pointer  ::  sggmiEy

#ifdef CompileWithProfiling
        call nvtxStartRange("Antes del bucle EY")
#endif        
        allocate(Idzh(0:this%bounds%dzh%NZ-1))
        allocate(Idyh(0:this%bounds%dxh%NX-1))
        allocate(sggmiEy(0:this%bounds%sggMiEy%NX-1,0:this%bounds%sggMiEy%NY-1,0:this%bounds%sggMiEy%NZ-1 ))
        Idhz => this%Idhz
        Idhx => this%Idhx
        sggmiEy => this%sggmiEy
        do k=1,this%bounds%sweepEy%NZ
            do j=1,this%bounds%sweepEy%NY
                do i=1,this%bounds%sweepEy%NX
                    medio = sggMiEy(i,j,k)
                    Ey(i,j,k) = this%g1(medio)*Ey(i,j,k) + & 
                                this%g2(medio)*((Hx(i,j,k)-Hx(i,j,k-1))*Idzh(k) - & 
                                                (Hz(i,j,k)-Hz(i-1,j,k))*Idxh(i))
                end do
            end do
        end do
#ifdef CompileWithProfiling
        call nvtxEndRange
#endif
    end subroutine

    subroutine advanceEz(this, Ez, Hx, Hy)
        class(stepper_t) :: this
        integer(kind = INTEGERSIZEOFMEDIAMATRICES) :: medio
        integer(kind=4) :: i, j, k
        real (kind=rkind), dimension(0:this%bounds%Ez%NX-1 , 0:this%bounds%Ez%NY-1 , 0:this%bounds%Ez%NZ-1), intent(inout) :: Ez
        real (kind=rkind), dimension(0:this%bounds%Hx%NX-1 , 0:this%bounds%Hx%NY-1 , 0:this%bounds%Hx%NZ-1), intent(inout) :: Hx
        real (kind=rkind), dimension(0:this%bounds%Hy%NX-1 , 0:this%bounds%Hy%NY-1 , 0:this%bounds%Hy%NZ-1), intent(inout) :: Hy
        real (kind=rkind), dimension(:), pointer  ::  Idxh
        real (kind=rkind), dimension(:), pointer  ::  Idzy
        real (kind=INTEGERSIZEOFMEDIAMATRICES), dimension(:,:,:), pointer  ::  sggmiEz
#ifdef CompileWithProfiling
        call nvtxStartRange("Antes del bucle EZ")
#endif        
        allocate(Idyh(0:this%bounds%dxh%NX-1))
        allocate(Idzh(0:this%bounds%dyh%NY-1))
        allocate(sggmiEz(0:this%bounds%sggMiEz%NX-1,0:this%bounds%sggMiEz%NY-1,0:this%bounds%sggMiEz%NZ-1 ))
        Idhx => this%Idhx
        Idhy => this%Idhy
        sggmiEz => this%sggmiEz
        do k=1,this%bounds%sweepEz%NZ
            do j=1,this%bounds%sweepEz%NY
                do i=1,this%bounds%sweepEz%NX
                    medio = sggMiEz(i,j,k)
                    Ez(i,j,k) = this%g1(medio)*Ez(i,j,k) + & 
                                this%g2(medio)*((Hy(i,j,k)-Hy(i-1,j,k))*Idzx(i) - & 
                                                (Hx(i,j,k)-Hx(i,j-1,k))*Idyh(j))
                end do
            end do
        end do
#ifdef CompileWithProfiling
        call nvtxEndRange
#endif
    end subroutine


end module timestepping_mod