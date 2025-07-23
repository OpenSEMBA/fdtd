!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module Lumped
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Elimino el tratamiento de los campos magneticos de Lumped para programar un multiLumped 
!!!solo teniendo en cuenta los parametros efectivos y sin actualizar los magneticos.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module Lumped

   use Report
   use fdetypes   
   use lumped_vars
#ifdef CompileWithStochastic
   use lumped_devia
#endif

   implicit none
   
   type (LumpedElem_t), save, target   ::  LumpElem
   
   REAL (KIND=RKIND), save           ::  eps0,mu0,zvac,cluz
   
   private
!!!
!!!
   public LumpedElem_t,Nodes_t 
   public AdvanceLumpedE,InitLumped,DestroyLumped,StoreFieldsLumpeds,calc_lumpedconstants,Getlumped

contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Subroutine to initialize the parameters
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine InitLumped(sgg,media,Ex,Ey,Ez,Hx,Hy,Hz,&
                         IDxe,IDye,IDze,IDxh,IDyh,IDzh, control, &
                         ThereAreLumped,eps00,mu00)
      REAL (KIND=RKIND)           ::  eps00,mu00
      type(media_matrices_t), intent(in) :: media

      type (SGGFDTDINFO), intent(IN)     ::  sgg
      REAL (KIND=RKIND)   , intent(in) , target     :: &
      Ex(sgg%alloc(iEx)%XI : sgg%alloc(iEx)%XE,sgg%alloc(iEx)%YI : sgg%alloc(iEx)%YE,sgg%alloc(iEx)%ZI : sgg%alloc(iEx)%ZE),&
      Ey(sgg%alloc(iEy)%XI : sgg%alloc(iEy)%XE,sgg%alloc(iEy)%YI : sgg%alloc(iEy)%YE,sgg%alloc(iEy)%ZI : sgg%alloc(iEy)%ZE),&
      Ez(sgg%alloc(iEz)%XI : sgg%alloc(iEz)%XE,sgg%alloc(iEz)%YI : sgg%alloc(iEz)%YE,sgg%alloc(iEz)%ZI : sgg%alloc(iEz)%ZE),&
      Hx(sgg%alloc(iHx)%XI : sgg%alloc(iHx)%XE,sgg%alloc(iHx)%YI : sgg%alloc(iHx)%YE,sgg%alloc(iHx)%ZI : sgg%alloc(iHx)%ZE),&
      Hy(sgg%alloc(iHy)%XI : sgg%alloc(iHy)%XE,sgg%alloc(iHy)%YI : sgg%alloc(iHy)%YE,sgg%alloc(iHy)%ZI : sgg%alloc(iHy)%ZE),&
      Hz(sgg%alloc(iHz)%XI : sgg%alloc(iHz)%XE,sgg%alloc(iHz)%YI : sgg%alloc(iHz)%YE,sgg%alloc(iHz)%ZI : sgg%alloc(iHz)%ZE)
      REAL (KIND=RKIND) , dimension (:)   , intent(in)   :: Idxh(sgg%ALLOC(iEx)%XI : sgg%ALLOC(iEx)%XE), &
                                                         &  Idyh(sgg%ALLOC(iEy)%YI : sgg%ALLOC(iEy)%YE), &
                                                         &  Idzh(sgg%ALLOC(iEz)%ZI : sgg%ALLOC(iEz)%ZE), &
                                                            Idxe(sgg%alloc(iHx)%XI : sgg%alloc(iHx)%XE), &
                                                            Idye(sgg%alloc(iHy)%YI : sgg%alloc(iHy)%YE), &
                                                            Idze(sgg%alloc(iHz)%ZI : sgg%alloc(iHz)%ZE)

      logical, INTENT(OUT)  ::  ThereAreLumped
      type(sim_control_t), intent(in) :: control
      integer (kind=4)  ::  jmed,j1,conta,k1,i1
      character(len=BUFSIZE) :: buff
      character (LEN=BUFSIZE)  ::  whoami
      type (Nodes_t), pointer :: lumped_
      logical :: unstable
!
      eps0=eps00; mu0=mu00; !chapuz para convertir la variables de paso en globales
!
!!!
      write(whoami,'(a,i5,a,i5,a)') '(',control%layoutnumber+1,'/',control%size,') '
      unstable=.false.
!
!
      ThereAreLumped=.FALSE.

      !precontaje

      conta=0
      Do k1=sgg%SINPMLSweep(iEx)%ZI,sgg%SINPMLSweep(iEx)%ZE
         Do j1=sgg%SINPMLSweep(iEx)%YI,sgg%SINPMLSweep(iEx)%YE
            Do i1=sgg%SINPMLSweep(iEx)%XI,sgg%SINPMLSweep(iEx)%XE
               jmed=media%sggMiEx(i1,j1,k1)
               if (SGG%Med(jmed)%Is%Lumped)  conta=conta+1
            end do
         end do
      end do
      Do k1=sgg%SINPMLSweep(iEy)%ZI,sgg%SINPMLSweep(iEy)%ZE
         Do j1=sgg%SINPMLSweep(iEy)%YI,sgg%SINPMLSweep(iEy)%YE
            Do i1=sgg%SINPMLSweep(iEy)%XI,sgg%SINPMLSweep(iEy)%XE
               jmed=media%sggMiEy(i1,j1,k1)
               if (SGG%Med(jmed)%Is%Lumped)  conta=conta+1
            end do
         end do
      end do

      Do k1=sgg%SINPMLSweep(iEz)%ZI,sgg%SINPMLSweep(iEz)%ZE
         Do j1=sgg%SINPMLSweep(iEz)%YI,sgg%SINPMLSweep(iEz)%YE
            Do i1=sgg%SINPMLSweep(iEz)%XI,sgg%SINPMLSweep(iEz)%XE
               jmed=media%sggMiEz(i1,j1,k1)
               if (SGG%Med(jmed)%Is%Lumped) conta=conta+1
            end do
         end do
      end do
      !!!!!!!!!!!!!!!!!!!!!!
      ThereAreLumped=(conta /=0)
      if (.not.thereareLumped) then
         return
      endif
      LumpElem%NumNodes=conta
      allocate (LumpElem%Nodes(1 : LumpElem%NumNodes))
      !!!!!!!!
      conta=0
      Do k1=sgg%SINPMLSweep(iEx)%ZI,sgg%SINPMLSweep(iEx)%ZE
         Do j1=sgg%SINPMLSweep(iEx)%YI,sgg%SINPMLSweep(iEx)%YE
            Do i1=sgg%SINPMLSweep(iEx)%XI,sgg%SINPMLSweep(iEx)%XE
               jmed=media%sggMiEx(i1,j1,k1)
               if (SGG%Med(jmed)%Is%Lumped)  then
                  conta=conta+1
                  lumped_ => LumpElem%Nodes(conta)
                  lumped_%alignedDeltaE      =    1.0_RKIND/IDXe(i1      )
                  lumped_%transversalDeltaHa =    1.0_RKIND/IDYh(   j1   )
                  lumped_%transversalDeltaHb    = 1.0_RKIND/IDzh(      k1)
                  lumped_%Orient=           SGG%Med(jmed)%Lumped(1)%Orient
                  lumped_%jmed       = jmed
                  lumped_%Efield    => Ex(i1,j1  ,k1  )
                  lumped_%Ha_Plus   => Hz(i1,j1  ,k1)
                  lumped_%Ha_Minu   => Hz(i1,j1-1,k1)
                  lumped_%Hb_Plus   => Hy(i1,j1  ,k1)
                  lumped_%Hb_Minu   => Hy(i1,j1  ,k1-1)
               endif
            end do
         end do
      end do
      !!!!!!!!!!
      Do k1=sgg%SINPMLSweep(iEy)%ZI,sgg%SINPMLSweep(iEy)%ZE
         Do j1=sgg%SINPMLSweep(iEy)%YI,sgg%SINPMLSweep(iEy)%YE
            Do i1=sgg%SINPMLSweep(iEy)%XI,sgg%SINPMLSweep(iEy)%XE
               jmed=media%sggMiEy(i1,j1,k1)
               if (SGG%Med(jmed)%Is%Lumped)  then
                  conta=conta+1
                  lumped_ => LumpElem%Nodes(conta)
                  lumped_%alignedDeltaE      = 1.0_RKIND/IDye(  j1   )
                  lumped_%transversalDeltaHa = 1.0_RKIND/IDzh(     k1)
                  lumped_%transversalDeltaHb = 1.0_RKIND/IDxh(i1     )
                  lumped_%Orient=           SGG%Med(jmed)%Lumped(1)%Orient
                  lumped_%jmed       = jmed
                  lumped_%Efield    => Ey(i1  ,j1  ,k1  )
                  lumped_%Ha_Plus   => Hx(i1  ,j1  ,k1  )
                  lumped_%Ha_Minu   => Hx(i1  ,j1  ,k1-1)
                  lumped_%Hb_Plus   => Hz(i1  ,j1  ,k1  )
                  lumped_%Hb_Minu   => Hz(i1-1,j1  ,k1  )
               endif
            end do
         end do
      end do

      Do k1=sgg%SINPMLSweep(iEz)%ZI,sgg%SINPMLSweep(iEz)%ZE
         Do j1=sgg%SINPMLSweep(iEz)%YI,sgg%SINPMLSweep(iEz)%YE
            Do i1=sgg%SINPMLSweep(iEz)%XI,sgg%SINPMLSweep(iEz)%XE
               jmed=media%sggMiEz(i1,j1,k1)
               if (SGG%Med(jmed)%Is%Lumped) then
                  conta=conta+1
                  lumped_ => LumpElem%Nodes(conta)
                  lumped_%alignedDeltaE      = 1.0_RKIND/IDzE(        k1)
                  lumped_%transversalDeltaHa = 1.0_RKIND/IDxh(i1        )
                  lumped_%transversalDeltaHb = 1.0_RKIND/IDyh(    j1    )
                  lumped_%Orient=           SGG%Med(jmed)%Lumped(1)%Orient
                  lumped_%jmed  = jmed
                  lumped_%Efield  => Ez(i1  ,j1  ,k1  )
                  lumped_%Ha_Plus => Hy(i1  ,j1  ,k1  )
                  lumped_%Ha_Minu => Hy(i1-1,j1  ,k1)
                  lumped_%Hb_Plus => Hx(i1  ,j1  ,k1  )
                  lumped_%Hb_Minu => Hx(i1  ,j1-1,k1  )
               endif
            end do
         end do
      end do

      call calc_lumpedconstants(sgg,eps0,mu0)   

      !!!!!!!!!resuming
      if (.not.control%resume) then  
         do conta=1,LumpElem%numnodes
            lumped_ => LumpElem%Nodes(conta)
            lumped_%EfieldPrevPrev=0.0_RKIND
            lumped_%EfieldPrev    =0.0_RKIND
            lumped_%Jcur          =0.0_RKIND
         end do
#ifdef CompileWithStochastic
         if (control%stochastic) then
             do conta=1,LumpElem%numnodes
                lumped_ => LumpElem%Nodes(conta)
            lumped_%EfieldPrevPrev_for_devia=0.0_RKIND
            lumped_%EfieldPrev_for_devia    =0.0_RKIND
            lumped_%Jcur_for_devia          =0.0_RKIND 
             end do
         endif
#endif
      else  
        do conta=1,LumpElem%numnodes
            lumped_ => LumpElem%Nodes(conta)
            read(14) lumped_%EfieldPrevPrev,lumped_%EfieldPrev,lumped_%Jcur 
        end do  
#ifdef CompileWithStochastic
         if (control%stochastic) then
             do conta=1,LumpElem%numnodes
                lumped_ => LumpElem%Nodes(conta)
                read(14) lumped_%EfieldPrevPrev_for_devia,lumped_%EfieldPrev_for_devia,lumped_%Jcur_for_devia
             end do
         endif
#endif
      endif
      return
   end subroutine InitLumped


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Subroutine to advance the E field in the Lumped: Usual Yee
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   subroutine AdvanceLumpedE(sgg,timestep,simu_devia,stochastic)
      logical :: simu_devia,stochastic 
      type (SGGFDTDINFO), intent(IN) ::  sgg
      integer (kind=4)  ::  conta,timestep
      real (KIND=RKIND) :: Enplus1, fieldC,A
      type (Nodes_t), pointer :: lumped_
      do conta=1,LumpElem%numnodes
         lumped_ => LumpElem%Nodes(conta)
         if (sgg%med(lumped_%jmed)%lumped(1)%inductor) then
             lumped_%EfieldPrevPrev = lumped_%EfieldPrev
             lumped_%EfieldPrev     = lumped_%Efield
             !!! The evolution of the current must be modified by a coefficient different from that shown in Mittra pag65.
             lumped_%Jcur = lumped_%currentCoeff * lumped_%Jcur + lumped_%sigmaEffResistInduct * (lumped_%EfieldPrev + lumped_%EfieldPrevPrev)
         else
             lumped_%Jcur=0.0_RKIND
         endif
         if (sgg%med(lumped_%jmed)%lumped(1)%diodo) then
             fieldC= - lumped_%G1 * lumped_%Efield -  (lumped_%G2a *(lumped_%Ha_Plus   - lumped_%Ha_Minu    ) - lumped_%G2b *(lumped_%Hb_Plus     - lumped_%Hb_Minu  ) ) - lumped_%diodepreA
             A= lumped_%diodepreA * exp(lumped_%diodeB * lumped_%Efield) 
             Enplus1 = newton_raphson(A,lumped_%diodeB,fieldC)
             lumped_%Efield = Enplus1
         else !debe entrar aqui si es un resistor, inductor o capacitor
            if (sgg%med(lumped_%jmed)%lumped(1)%resistor) then
                if ((timestep*sgg%dt >= sgg%Med(lumped_%jmed)%Lumped(1)%Rtime_on).and.(timestep*sgg%dt <= sgg%Med(lumped_%jmed)%Lumped(1)%Rtime_off)) then
                   lumped_%Efield = lumped_%G1 * lumped_%Efield +  (lumped_%G2a *(lumped_%Ha_Plus   - lumped_%Ha_Minu    ) - lumped_%G2b *(lumped_%Hb_Plus     - lumped_%Hb_Minu  ) ) - &
                                    lumped_%GJ * lumped_%Jcur
                else
                   lumped_%Efield = lumped_%G1_usual * lumped_%Efield +(lumped_%G2a_usual *(lumped_%Ha_Plus - lumped_%Ha_Minu) - &
                                                                        lumped_%G2b_usual *(lumped_%Hb_Plus - lumped_%Hb_Minu))
                endif
            else !inductor o capacitor
                lumped_%Efield = lumped_%G1 * lumped_%Efield +  (lumped_%G2a *(lumped_%Ha_Plus   - lumped_%Ha_Minu    ) - lumped_%G2b *(lumped_%Hb_Plus     - lumped_%Hb_Minu  ) ) - &
                                 lumped_%GJ * lumped_%Jcur
            endif
         endif                     
#ifdef CompileWithStochastic
         call inject_devialumped(sgg,timestep,simu_devia,stochastic,lumped_)
!inyecta las devia como sources
#endif
!!!
      end do
      return
   end subroutine AdvanceLumpedE


   subroutine calc_lumpedconstants(sgg,eps00,mu00)      
      REAL (KIND=RKIND)           ::  eps00,mu00
      type (SGGFDTDINFO), intent(IN)     ::  sgg
      integer (kind=4)  ::  conta
      type (Nodes_t), pointer :: lumped_
!!!variables locales
      integer (kind=4) :: jmed
      integer (kind=4) :: orient
      real (kind=RKIND) :: epsilon,sigma,g1,g2,Resist,Induct,Capaci,sigmaeff,epsiloneff,DiodB,DiodIsat
      real (kind=RKIND) :: g1_usual,g2_usual
      real (kind=RKIND) :: epsilonEffCapac    ,sigmaEffResistInduct,sigmaEffResist   ,sigmaEffResistCapac ,sigmaEffResistDiode, &
                           alignedDeltaE,transversalDeltaHa,transversalDeltaHb,  currentCoeff
      
      character(len=BUFSIZE) :: buff
!
      eps0=eps00; mu0=mu00; !chapuz para convertir la variables de paso en globales
      zvac=sqrt(mu0/eps0)
      cluz=1.0_RKIND/sqrt(eps0*mu0)
!
!!!

        do conta=1,LumpElem%numnodes
            lumped_ => LumpElem%Nodes(conta)
            jmed=lumped_%jmed

            orient=               sgg%Med(jmed)%Lumped(1)%orient
            Resist=               sgg%Med(jmed)%Lumped(1)%R
            Induct=               sgg%Med(jmed)%Lumped(1)%L
            Capaci=               sgg%Med(jmed)%Lumped(1)%C
            DiodB=                sgg%Med(jmed)%Lumped(1)%DiodB
            DiodIsat=             sgg%Med(jmed)%Lumped(1)%DiodIsat
            epsilon            =  sgg%Med(jmed)%epr * eps0
            sigma              =  sgg%Med(jmed)%sigma
            alignedDeltaE      =lumped_%alignedDeltaE      
            transversalDeltaHa =lumped_%transversalDeltaHa 
            transversalDeltaHb =lumped_%transversalDeltaHb     
!
            epsilonEffCapac      = alignedDeltaE * Capaci / (            transversalDeltaHa * transversalDeltaHb)
            sigmaEffResistInduct = alignedDeltaE * sgg%dt / (2.0_RKIND * transversalDeltaHa * transversalDeltaHb * (Induct + Resist * sgg%dt /2.0_RKIND))
            sigmaEffResist       = alignedDeltaE          / (   Resist * transversalDeltaHa * transversalDeltaHb)
            sigmaEffResistCapac  = sigmaEffResist
            sigmaEffResistDiode  = sigmaEffResist
            currentCoeff         = (Induct - Resist * sgg%dt /2.0_RKIND) / (Induct + Resist * sgg%dt /2.0_RKIND)

!!!! Mittra pag65   Parallel Finite-Difference Time-Domain Method
            if (sgg%med(jmed)%lumped(1)%resistor) then
                sigmaeff= sigma     + sigmaEffResist                     
                epsiloneff= epsilon 
            elseif (sgg%med(jmed)%lumped(1)%inductor) then
                sigmaeff= sigma     + sigmaEffResistInduct                     
                epsiloneff= epsilon 
            elseif (sgg%med(jmed)%lumped(1)%capacitor) then
                sigmaeff= sigma     + sigmaEffResistCapac                    
                epsiloneff= epsilon + epsilonEffCapac
            elseif (sgg%med(jmed)%lumped(1)%diodo) then
                sigmaeff= sigma     + sigmaEffResistDiode                    
                epsiloneff= epsilon 
            endif 
            if (.not.sgg%Med(jmed)%sigmareasignado) then
                sgg%Med(jmed)%sigma = sigmaeff !devuelve al principal el efectivo para hacer bien las conexiones lossy con thin wires 120123
                sgg%Med(jmed)%sigmareasignado=.true.
            else
                print *,'error buggy: reasignando sigma en un lumped'
                stop
            endif
!ORIGINALES COMENTADOS 151222 por otros equivalentes para simplificar implem stoch            
!           G1=(1.0_RKIND  - SigmaEff * sgg%dt / (2.0_RKIND * epsilonEff) ) / &
!              (1.0_RKIND  + SigmaEff * sgg%dt / (2.0_RKIND * epsilonEff) ) 
!           G2=  sgg%dt / epsilonEff                        / &
!              (1.0_RKIND  + SigmaEFF * sgg%dt / (2.0_RKIND * epsilonEff) ) 
!            
            G1=(epsilonEff/ sgg%dt  - SigmaEff / 2.0_RKIND ) / &
               (epsilonEff/ sgg%dt  + SigmaEff / 2.0_RKIND ) 
            G2=  1.0_RKIND  / &
               (epsilonEff/sgg%dt   + SigmaEff/2.0_RKIND  ) 
            
            
            !!!lo he comentado a 050122 por consistencia con stochastic
            !!if (g1 < 0.0_RKIND) then !exponential time stepping
            !!    g1=exp(- SigmaEff * sgg%dt / (epsilonEff ))
            !!    g2=(1.0_RKIND-g1)/ SigmaEff
            !!endif
            lumped_%g1=g1 
            lumped_%G2a= G2 / transversalDeltaHa 
            lumped_%G2b= G2 / transversalDeltaHb
            !!! The current derivation on page 65 of Mittra is incorrect and leads to saturation. 
            !!! The correct coefficient is the one shown here.
            lumped_%GJ = G2 * (1 + currentCoeff) / 2 !!!!Mittra pag65
            lumped_%sigmaEffResistInduct = sigmaEffResistInduct
            lumped_%currentCoeff = currentCoeff
  
            !!!!usual para resistencia que se encienden/apagan 200319
            G1_usual=(1.0_RKIND  - Sigma * sgg%dt / (2.0_RKIND * epsilon) ) / &
                (1.0_RKIND  + Sigma * sgg%dt / (2.0_RKIND * epsilon) ) 
            G2_usual=  sgg%dt / epsilon                        / &
                (1.0_RKIND  + Sigma * sgg%dt / (2.0_RKIND * epsilon) ) 

            if (g1_usual < 0.0_RKIND) then !exponential time stepping
                g1_usual=exp(- Sigma * sgg%dt / (epsilon ))
                g2_usual=(1.0_RKIND-g1_usual)/ Sigma
            endif
            lumped_%g1_usual=g1_usual
            lumped_%G2a_usual= G2_usual / transversalDeltaHa 
            lumped_%G2b_usual= G2_usual / transversalDeltaHb
                  
           !!!only for diodes 
            if (orient>0.0) then
                lumped_%diodeB    = lumped_%diodeB * alignedDeltaE / 2.0_RKIND
                lumped_%diodepreA = DiodIsat * G2 / ( transversalDeltaHa * transversalDeltaHb)  !DiodIsat es la corriente de saturacions
            elseif (orient<0.0) then
                lumped_%diodeB    = -lumped_%diodeB * alignedDeltaE / 2.0_RKIND
                lumped_%diodepreA = -DiodIsat * G2 / ( transversalDeltaHa * transversalDeltaHb) 
            else
                WRITE (buff, *)    'Buggy ERROR: In lumped orientations. '
                CALL WarnErrReport (buff,.TRUE.)
            endif   
        
#ifdef CompileWithStochastic
            call calc_lumped_deviaconsts(sgg,lumped_,jmed,Resist,Induct,Capaci,sigmaeff,epsiloneff,eps0,mu0)

#endif   
        end do
    return
!!!!!
   end subroutine calc_lumpedconstants

   subroutine StoreFieldsLumpeds(stochastic)
         logical :: stochastic
         integer (kind=4) :: conta
         type (Nodes_t), pointer :: lumped_
         do conta=1,LumpElem%numnodes
            lumped_ => LumpElem%Nodes(conta)
            write(14,err=634) lumped_%EfieldPrevPrev,lumped_%EfieldPrev,lumped_%Jcur !olvide almacenar jcur 071118
         end do
#ifdef CompileWithStochastic
         if (stochastic) then
             do conta=1,LumpElem%numnodes
                lumped_ => LumpElem%Nodes(conta)
                write(14,err=634) lumped_%EfieldPrevPrev_for_devia,lumped_%EfieldPrev_for_devia,lumped_%Jcur_for_devia !olvide almacenar jcur 071118
             end do
         endif
#endif
      goto 635
634   call print11(0,SEPARADOR//separador//separador)
      call print11(0,'LUMPED: ERROR WRITING RESTARTING FIELDS. IGNORING AND CONTINUING')
      call print11(0,SEPARADOR//separador//separador)          
635   return
   end subroutine StoreFieldsLumpeds


   subroutine DestroyLumped(sgg)

      type (SGGFDTDINFO), intent(INOUT)         ::  sgg
      integer (kind=4)  ::  i

      !free up memory
      do i=1,sgg%NumMedia
         if ((sgg%Med(i)%Is%Lumped).and.(.not.sgg%Med(i)%Is%PML)) then
            if (associated(sgg%Med(i)%Lumped)) deallocate (sgg%Med(i)%Lumped)
         endif
      end do


      if(allocated(LumpElem%nodes)) deallocate (LumpElem%nodes)

   end subroutine

   function newton_raphson(A,B,C) result (x)
   real (KIND=RKIND), intent (in) :: A,B,C
   real (KIND=RKIND) :: x
   real (KIND=RKIND) :: x0, xx0, fxx0, dfxx0
   real (KIND=RKIND) :: tol ! Tolerancia error relativo
   INTEGER, parameter :: nmax=1024  !limite de iteraciones/iteraciones realizadas
   INTEGER :: clave,i,n      ! Clave de exito

    clave = 1
    xx0 = x0
    busca: DO i=1,nmax
        fxx0=A*exp(B*xx0)+x+C
        dfxx0=A*B*exp(B*xx0)+1.0_RKIND
        x = xx0 - fxx0/dfxx0
        IF (ABS(x-xx0) < tol*ABS(x) ) THEN
            clave = 0
            n = i
            EXIT busca
        ENDIF
        xx0 = x
    END DO busca
    if (clave /= 0) CALL print11 (0, 'Error convergencia en Newton-Raphson')
    return
   end function newton_raphson

   function Getlumped() result(r)
      type(LumpedElem_t), pointer  ::  r
      r=>LumpElem
      return
   end function
      
   
   
   
end module Lumped
