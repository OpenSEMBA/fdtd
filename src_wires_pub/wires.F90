!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MIT License
! 
! Copyright (c) 2023 University of Granada
! 
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
! 
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
! 
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module thin wires from Wires paper
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module HollandWires
   !
#ifdef CompileWithWires

   use report
   use fdetypes
   use wiresHolland_constants

#ifdef CompileWithStochastic   
   use wiresHolland_devia
#endif


#ifdef CompileWithThickWires  
    use Thick_m
#endif    
   !
   implicit none
   
   
   real (KIND=RKIND_wires), PARAMETER             :: HEUR_RADIUSOVERDELTA=10.0
   !local variables

   logical                                         , save :: thereAreVsources,thereAreIsources,thereAreMurConditions
   type(Thinwires_t)     , target                  ,save  ::  HWires
   REAL (KIND=RKIND_wires)     , pointer, dimension ( : ),save  ::  InvEps  ,InvMu, OldInvEps  ,OldInvMu
   
!!!variables globales del modulo
   REAL (KIND=RKIND_wires), save           ::  eps0,mu0
!!!
   private

   public InitWires,AdvanceWiresE,AdvanceWiresH,AdvanceWiresEcrank,StoreFieldsWires,DestroyWires, GetHwires,ReportWireJunctions,calc_wirehollandconstants



contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Subroutine to initialize the parameters
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine InitWires(sgg,sggMiNo,sggMiEx,sggMiEy,sggMiEz,sggMiHx,sggMiHy,sggMiHz,layoutnumber,size,ThereAreWires,resume,makeholes,connectendings,isolategroupgroups,stableradholland,fieldtotl, &
   Ex,Ey,Ez,Hx,Hy,Hz,Idxe,Idye,Idze,Idxh,Idyh,Idzh, &
   inductance_model,wirethickness,groundwires,strictOLD,TAPARRABOS,g2,wiresflavor,SINPML_fullsize,fullsize,wirecrank,dtcritico, &
   eps00,mu00,simu_devia,stochastic,verbose,factorradius,factordelta)     
      integer, intent(IN) :: wirethickness
      logical :: simu_devia,stochastic,verbose
      REAL (KIND=RKIND), intent(in)           ::  eps00,mu00
      REAL (KIND=RKIND)           ::  eps000,mu000  !son dummies
      integer (kind=4), intent(in) :: layoutnumber,size
      type (limit_t), dimension(1:6), intent(in)  :: SINPML_fullsize,fullsize
      character(len=*), INTENT(in) :: inductance_model
      character(len=*), INTENT(in) :: wiresflavor
      type (SGGFDTDINFO), intent(INOUT) , target    ::  sgg
      REAL (KIND=RKIND) , pointer, dimension (:), intent(in)      :: G2
      REAL (KIND=RKIND), intent(out) :: dtcritico
      REAL (KIND=RKIND) , dimension (:)   , intent(in)      ::  &
           Idxe(sgg%ALLOC(iHx)%XI : sgg%ALLOC(iHx)%XE), &
           Idye(sgg%ALLOC(iHy)%YI : sgg%ALLOC(iHy)%YE), &
           Idze(sgg%ALLOC(iHz)%ZI : sgg%ALLOC(iHz)%ZE), &
           Idxh(sgg%ALLOC(iEx)%XI : sgg%ALLOC(iEx)%XE), &
           Idyh(sgg%ALLOC(iEy)%YI : sgg%ALLOC(iEy)%YE), &
           Idzh(sgg%ALLOC(iEz)%ZI : sgg%ALLOC(iEz)%ZE)
                                                           
        REAL (KIND=RKIND)   , intent(in), target      :: &
            Ex(sgg%Alloc(iEx)%XI : sgg%Alloc(iEx)%XE,sgg%Alloc(iEx)%YI : sgg%Alloc(iEx)%YE,sgg%Alloc(iEx)%ZI : sgg%Alloc(iEx)%ZE),&
            Ey(sgg%Alloc(iEy)%XI : sgg%Alloc(iEy)%XE,sgg%Alloc(iEy)%YI : sgg%Alloc(iEy)%YE,sgg%Alloc(iEy)%ZI : sgg%Alloc(iEy)%ZE),&
            Ez(sgg%Alloc(iEz)%XI : sgg%Alloc(iEz)%XE,sgg%Alloc(iEz)%YI : sgg%Alloc(iEz)%YE,sgg%Alloc(iEz)%ZI : sgg%Alloc(iEz)%ZE),&
            Hx(sgg%alloc(iHx)%XI : sgg%alloc(iHx)%XE,sgg%alloc(iHx)%YI : sgg%alloc(iHx)%YE,sgg%alloc(iHx)%ZI : sgg%alloc(iHx)%ZE),&
            Hy(sgg%alloc(iHy)%XI : sgg%alloc(iHy)%XE,sgg%alloc(iHy)%YI : sgg%alloc(iHy)%YE,sgg%alloc(iHy)%ZI : sgg%alloc(iHy)%ZE),&
            Hz(sgg%alloc(iHz)%XI : sgg%alloc(iHz)%XE,sgg%alloc(iHz)%YI : sgg%alloc(iHz)%YE,sgg%alloc(iHz)%ZI : sgg%alloc(iHz)%ZE)


      integer (KIND=INTEGERSIZEOFMEDIAMATRICES), intent(in) :: &
              sggMiNo(sgg%alloc(iHx)%XI : sgg%alloc(iHx)%XE,sgg%alloc(iHy)%YI : sgg%alloc(iHy)%YE,sgg%alloc(iHz)%ZI : sgg%alloc(iHz)%ZE)
      integer (KIND=INTEGERSIZEOFMEDIAMATRICES), intent(inout)   ::  &
              sggMiEx(sgg%Alloc(iEx)%XI : sgg%Alloc(iEx)%XE,sgg%Alloc(iEx)%YI : sgg%Alloc(iEx)%YE,sgg%Alloc(iEx)%ZI : sgg%Alloc(iEx)%ZE), &
              sggMiEy(sgg%Alloc(iEy)%XI : sgg%Alloc(iEy)%XE,sgg%Alloc(iEy)%YI : sgg%Alloc(iEy)%YE,sgg%Alloc(iEy)%ZI : sgg%Alloc(iEy)%ZE), &
              sggMiEz(sgg%Alloc(iEz)%XI : sgg%Alloc(iEz)%XE,sgg%Alloc(iEz)%YI : sgg%Alloc(iEz)%YE,sgg%Alloc(iEz)%ZI : sgg%Alloc(iEz)%ZE), &
              sggMiHx(sgg%Alloc(iHx)%XI : sgg%Alloc(iHx)%XE,sgg%Alloc(iHx)%YI : sgg%Alloc(iHx)%YE,sgg%Alloc(iHx)%ZI : sgg%Alloc(iHx)%ZE), &
              sggMiHy(sgg%Alloc(iHy)%XI : sgg%Alloc(iHy)%XE,sgg%Alloc(iHy)%YI : sgg%Alloc(iHy)%YE,sgg%Alloc(iHy)%ZI : sgg%Alloc(iHy)%ZE), &
              sggMiHz(sgg%Alloc(iHz)%XI : sgg%Alloc(iHz)%XE,sgg%Alloc(iHz)%YI : sgg%Alloc(iHz)%YE,sgg%Alloc(iHz)%ZI : sgg%Alloc(iHz)%ZE)
      integer (KIND=INTEGERSIZEOFMEDIAMATRICES) :: sggMiE 
      
      
      logical , dimension (:), pointer :: LindProb

      LOGICAL, INTENT(OUT)  ::  ThereAreWires
      LOGICAL, INTENT(in)  ::  resume,makeholes,connectendings,isolategroupgroups,groundwires,stableradholland,strictOLD,TAPARRABOS,fieldtotl
      logical :: proceed,proceed1,proceed2,NodeExists,Is_LeftEnd,Is_RightEnd,IsEnd_norLeft_norRight,repetido,conectado,conectado1,conectado2,asignado
      logical ::  IsPEC , islossy ,IsLossyPlus,IsLossyMinu,IsPecPlus,IsPECminu
      logical, intent (in) :: wirecrank
      real (kind=RKIND_wires) :: rlossy,newr0, factorradius,factordelta
      integer (KIND=INTEGERSIZEOFMEDIAMATRICES) :: med(0:11)=-1

      type (adyacc)  ::  adj
      integer (kind=4)  ::  conta,i1,j1,k1,i2,j2,k2,iwi,iwj,iwjjj,jmed,nn,nnn,i1libre,j1libre,k1libre, &
      whatfield,whatfield2,origIndex,OrigIndex2,LeftEnd_index,RightEnd_index,nm, &
      i,j,k,indexnode,kmenos1,kmasoffk,kmas1,tipofield,i22,j22,k22,i11,j11,k11,primernorabo,Jprimernorabo=-1
      REAL (KIND=RKIND_wires)   ::  r0, desp, deltadummy1 ,deltadummy2, deltadummy, oldr0,a,b, &
      despT1,despT2,DenominatorFractionMinusDummy,  &
      DenominatorFractionPlusDummy,givenautoin,resist,givenautoin_devia,resist_devia, &
      mindt,maxA,dt0,sigt,totalLind,capaci,autoin,deltax,sigtPlus,sigtMinu
      type (CurrentSegments), pointer  ::  dummy , org  ,fin ,orgmenos1 ,orgmas1 ,finmenos1,finmas1 ,segmento
      type (ChargeNodes), pointer  ::  nodo

      !dama
      integer   (KIND=4)                          ::  NumMultilines, NumParallel
      integer   (KIND=4)                          ::  iw1, is1, in1
      integer   (KIND=4)                          ::  iw2, is2, in2
      integer   (KIND=4)                          ::  contsgm, contmtln, contprll
      integer   (KIND=4)                          ::  contnds, contchgm, contchgp
      integer   (KIND=4)                          ::  imed, N
      real      (kind=RKIND_wires)							::  dl, dx1, dx2,Ceq
      real      (kind=RKIND_wires)							::  dist, phi
      real      (KIND=RKIND_wires), pointer, dimension(:,:) ::  Den
      !dama fin
      logical :: esPML
      integer (kind=4)  ::  int1, int2,ierr,mediox,medioy,medioz,ZI,ZE,offset,offi,offj,offk,NUMESEG,dummy1,dummy2,dummy3,multirabos,dummyfin,medio1,medio2,medio3,medio1m,medio2m,medio3m
      character(len=BUFSIZE) :: buff

      REAL (KIND=RKIND_wires)           ::  df1,df3,df2,Ddf1,Ddf3,Ddf2,vf1,vf3,vf2,runit
      character (LEN=BUFSIZE)  ::  whoami
      character (len=3), dimension(1:3) :: DIR
      
      !!!
!
      eps0=eps00; mu0=mu00; !chapuz para convertir la variables de paso en globales
!!!guarda el dt original (para permit scaling)
      Hwires%olddt=sgg%dt
!

      HWires%WireTipoMedio   => null()
      HWires%CurrentSegment  => null()
      HWires%ChargeNode      => null()
            !!!!!!!!!!!
      dtcritico=sgg%dt !es un parametro de vuelta

      ZI=sgg%Sweep(iHz)%ZI
      ZE=sgg%Sweep(iHz)%ZE


      dir(iEx)=' X '
      dir(iEy)=' Y '
      dir(iEz)=' Z '
      thereAreVsources=.false.
      thereAreIsources=.false.
      thereAreMurConditions=.false.

      write(whoami,'(a,i5,a,i5,a)') '(',layoutnumber+1,'/',size,') '

#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif


      adj%YESsegment(1:2) = -1 !default

      allocate (InvEps(0 : sgg%NumMedia),InvMu(0 : sgg%NumMedia))
      InvEps(0 : sgg%NumMedia)=1.0_RKIND_wires/(Eps0*sgg%Med(0 : sgg%NumMedia)%Epr) 
      InvMu(0 : sgg%NumMedia)=1.0_RKIND_wires/(Mu0*sgg%Med(0 : sgg%NumMedia)%Mur) 
      allocate (OldInvEps(0 : sgg%NumMedia),OldInvMu(0 : sgg%NumMedia))
      OldInvEps(0 : sgg%NumMedia)=InvEps(0 : sgg%NumMedia)
      OldInvMu(0 : sgg%NumMedia)=InvMu(0 : sgg%NumMedia)

#ifdef CompileWithMPI
      !Initialize MPI counters (later written to disk)
      Hwires%NumNeededCurrentUpMPI=0
      Hwires%NumNeededCurrentDownMPI=0
#endif
!!inicializaciones triviales

         !initialization
         !el trivial
         !
         HWires%NullNode%indexnode           =-1
         HWires%NullNode%ChargePresent       =0.0_RKIND_wires
         HWires%NullNode%ChargePast          =0.0_RKIND_wires
         HWires%NullNode%IsAttachedtoVoltage =.false.
         HWires%NullNode%IsMur               =.false.
         HWires%NullNode%IsBackDownLeftMur =.false.
         HWires%NullNode%IsFrontUpRightMur =.false.
         HWires%NullNode%IsPeriodic               =.false.
         HWires%NullNode%IsPEC               =.false.
         HWires%NullNode%already_YEEadvanced_byconformal_changedtoPECfield1 => null()
         HWires%NullNode%already_YEEadvanced_byconformal_changedtoPECfield2 => null()
         HWires%NullNode%already_YEEadvanced_byconformal_changedtoPECfield3 => null()
         HWires%NullNode%already_YEEadvanced_byconformal_changedtoPECfield4 => null()
         HWires%NullNode%already_YEEadvanced_byconformal_changedtoPECfield5 => null()
         HWires%NullNode%already_YEEadvanced_byconformal_changedtoPECfield6 => null()
         HWires%NullNode%IsLossy             =.false.
         HWires%NullNode%HasIsource          =.false.
         HWires%NullNode%IsHeterogeneousJunction          =.false.
         !just for informative !not implemented in MPI unsure behaviour under mpi 2011 \E7
         HWires%NullNode%Exists              =.false.
         HWires%NullNode%proc   =.false.
         HWires%NullNode%Is_LeftEnd   =.false.
         HWires%NullNode%Is_RightEnd   =.false.
         HWires%NullNode%IsInSingleRLCsegment   =.false.
         HWires%NullNode%NumCurrentMinus =0
         HWires%NullNode%NumCurrentPlus  =0
         HWires%NullNode%i                   =-1
         HWires%NullNode%j                   =-1
         HWires%NullNode%k                   =-1
         HWires%NullNode%cteMur              =0.0_RKIND_wires
         HWires%NullNode%cteProp            =0.0_RKIND_wires
         HWires%NullNode%oRIGctePlain            =0.0_RKIND_wires
         HWires%NullNode%ctePlain            =0.0_RKIND_wires
         !
         HWires%NullSegment%R      = 0.0_RKIND_wires
         HWires%NullSegment%Resist      = 0.0_RKIND_wires
         HWires%NullSegment%Resist_devia      = 0.0_RKIND_wires
         HWires%NullSegment%C      = 0.0_RKIND_wires
         HWires%NullSegment%L      = 0.0_RKIND_wires
         HWires%NullSegment%Lintrinsic    = 0.0_RKIND_wires
         HWires%NullSegment%NumParallel    =1
         HWires%NullSegment%origindex       =i1
         HWires%NullSegment%indexsegment    =i1
         HWires%NullSegment%currentpast         =0.0_RKIND_wires
         HWires%NullSegment%current         =0.0_RKIND_wires
         HWires%NullSegment%qplus_qminus         =0.0_RKIND_wires
         HWires%NullSegment%current_for_devia =0.0_RKIND_wires
         HWires%NullSegment%qplus_qminus_for_devia =0.0_RKIND_wires 
         HWires%NullSegment%Efield_main2wire_for_devia         =0.0_RKIND_wires
         HWires%NullSegment%inv_Lind_acum            =0.0_RKIND_wires
         HWires%NullSegment%Lind_acum            =0.0_RKIND_wires
         HWires%NullSegment%Lind            =0.0_RKIND_wires
         HWires%NullSegment%Lind_devia            =0.0_RKIND_wires
         HWires%NullSegment%HEUR_safety =0.0_RKIND_wires
         !!! HWires%NullSegment%logRoverR0      =0.0_RKIND_wires
         HWires%NullSegment%delta           =0.0_RKIND_wires
         HWires%NullSegment%deltaTransv1    =0.0_RKIND_wires
         HWires%NullSegment%deltaTransv2    =0.0_RKIND_wires
         HWires%NullSegment%cte1            =0.0_RKIND_wires
         HWires%NullSegment%cte2            =0.0_RKIND_wires
         HWires%NullSegment%cte3            =0.0_RKIND_wires
         HWires%NullSegment%cte1_for_devia            =0.0_RKIND_wires
         HWires%NullSegment%cte2_for_devia            =0.0_RKIND_wires
         HWires%NullSegment%cte3_for_devia            =0.0_RKIND_wires
         HWires%NullSegment%cte5            =0.0_RKIND_wires
         HWires%NullSegment%ilibre               =-1
         HWires%NullSegment%jlibre               =-1
         HWires%NullSegment%klibre               =-1
         HWires%NullSegment%i               =-1
         HWires%NullSegment%j               =-1
         HWires%NullSegment%k               =-1
         HWires%NullSegment%tipofield       =-1
         HWires%NullSegment%IsPMC           =.false.
         HWires%NullSegment%orientadoalreves =.false.
         HWires%NullSegment%HasVsource      =.false.
         HWires%NullSegment%IsShielded      =.false.
         HWires%NullSegment%HasAbsorbing_RightEnd =.false.
         HWires%NullSegment%HasAbsorbing_LeftEnd =.false.
         HWires%NullSegment%HasParallel_RightEnd =.false.
         HWires%NullSegment%HasParallel_LeftEnd =.false.
         HWires%NullSegment%HasSeries_RightEnd   =.false.
         HWires%NullSegment%HasSeries_LeftEnd   =.false.
         HWires%NullSegment%IsEnd_norLeft_norRight   =.false.
         HWires%NullSegment%Is_LeftEnd   =.false.
         HWires%NullSegment%Is_RightEnd   =.false.
         HWires%NullSegment%chargePlus => HWires%NullNode
         HWires%NullSegment%chargeMinus => HWires%NullNode
         !!!!!!!!!!
         !!!


      ThereAreWires=.FALSE.

      !detect thin wires :  same radius implies same medium independently of its orientation
      conta=0
      do jmed=1,sgg%NumMedia
         if (sgg%Med(jmed)%Is%ThinWire) conta=conta+1
      end do

      HWires%NumDifferentWires=conta
      allocate (HWires%WireTipoMedio(1 : HWires%NumDifferentWires))
      conta=0
      do jmed=1,sgg%NumMedia
         if (sgg%Med(jmed)%Is%ThinWire) then
            ThereAreWires=.TRUE.
            conta=conta+1
            HWires%WireTipoMedio(conta)=jmed
         endif
      end do

!!!chequeo previo 210323 para que no haya conectores dispersivos
         
      do iwi=1,HWires%NumDifferentWires
            if ((allocated(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%disp)).or. &
                (allocated(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%disp_RightEnd )).or. &
                (allocated(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%disp_LeftEnd  ))) then
                write (buff,'(a)')  'Dispersive wire or connectors unsupported in Holland wires'
                call WarnErrReport(buff,.true.)
            endif
      end do

      
      do iwi=1,HWires%NumDifferentWires
         do iwj=1,sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%numsegmentos
            sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%multirabo=.false.
         end do
      end do


      if (.not.therearewires) Return

      write (buff,'(a)')  '----------------------------------------------------------------'
      call WarnErrReport(buff)

      ! it directly reads the segments specified in the .nfde file

         !detects endings and set ending=.true.
         !esto implica que podra haber LeftEnd y RightEnd declarados y ademas ending 


         do iwi=1,HWires%NumDifferentWires
            do iwj=1,sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%numsegmentos
               if (.not.strictOLD) then
                  !
                  sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%ilibre = -1
                  sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%jlibre = -1
                  sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%klibre = -1
                  !(en formatos antiguas esta info solo viene bien si hay cargas (md me las pone al final y al principio) pero en general puede estar mal y hay que resetearla
                  if (connectendings) then
                     !
                     if (sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%numsegmentos == 1) then
                        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%Is_LeftEnd = .true.
                        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%Is_RightEnd = .true.
                        !ojooo esto esta bien? parecen los TL y TR intercambiados sgg 251019 pero no lo toco
                        if (sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%HasAbsorbing_RightEnd) sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%HasAbsorbing_LeftEnd=.true.
                        if (sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%HasAbsorbing_LeftEnd) sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%HasAbsorbing_RightEnd=.true.
                        if (sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%HasParallel_RightEnd) sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%HasParallel_LeftEnd=.true.
                        if (sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%HasParallel_LeftEnd) sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%HasParallel_RightEnd=.true.
                        if (sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%HasSeries_RightEnd  ) sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%HasSeries_LeftEnd  =.true.
                        if (sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%HasSeries_LeftEnd  ) sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%HasSeries_RightEnd  =.true.
                     endif
                     !
                     sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%Is_LeftEnd = sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%Is_LeftEnd .and. &
                     (sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%HasParallel_LeftEnd .or. sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%HasSeries_LeftEnd.or. &
                                                                                  sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%HasAbsorbing_LeftEnd)
                     sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%Is_RightEnd = sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%Is_RightEnd .and. &
                     (sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%HasParallel_RightEnd .or. sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%HasSeries_RightEnd .or. &
                                                                                     sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%HasAbsorbing_RightEnd)
                  endif
                  sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%IsEnd_norLeft_norRight= .not. &
                  (sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%Is_LeftEnd.or.sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%Is_RightEnd)
                  !fin reajuste
                  i1=       sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%i
                  j1=       sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%j
                  k1=       sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%k
                  whatfield=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%ori
                  origindex= sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%origIndex
                  !
                  conectado1 = .false.
                  conectado2 = .false.
                  buskk: do iwjjj=1,sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%numsegmentos
                     i2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%i
                     j2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%j
                     k2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%k
                     whatfield2=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%ori
                     origindex2= sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%origIndex
                     !
                     select case (whatfield) !cyclic
                      case (iEx)
                        select case (whatfield2)
                         case (iEx)
                           conectado1 = conectado1 .or. ((i1   == i2+1).and.(j1 == j2  ).and.(k1 == k2  ))
                           conectado2 = conectado2 .or. ((i1 +1== i2  ).and.(j1 == j2  ).and.(k1 == k2  ))
                         case (iEy)
                           conectado1 = conectado1 .or. ((i1   == i2  ).and.(j1 == j2  ).and.(k1 == k2  ))
                           conectado1 = conectado1 .or. ((i1   == i2  ).and.(j1 == j2+1).and.(k1 == k2  ))
                           !
                           conectado2 = conectado2 .or. ((i1 +1== i2  ).and.(j1 == j2  ).and.(k1 == k2  ))
                           conectado2 = conectado2 .or. ((i1 +1== i2  ).and.(j1 == j2+1).and.(k1 == k2  ))
                         case (iEz)
                           conectado1 = conectado1 .or. ((i1   == i2  ).and.(j1 == j2  ).and.(k1 == k2  ))
                           conectado1 = conectado1 .or. ((i1   == i2  ).and.(j1 == j2  ).and.(k1 == k2+1))
                           !
                           conectado2 = conectado2 .or. ((i1 +1== i2  ).and.(j1 == j2  ).and.(k1 == k2  ))
                           conectado2 = conectado2 .or. ((i1 +1== i2  ).and.(j1 == j2  ).and.(k1 == k2+1))
                        end select
                      case (iEy)
                        select case (whatfield2)
                         case (iEy)
                           conectado1 = conectado1 .or. ((j1   == j2+1).and.(k1 == k2  ).and.(i1 == i2  ))
                           conectado2 = conectado2 .or. ((j1 +1== j2  ).and.(k1 == k2  ).and.(i1 == i2  ))
                         case (iEz)
                           conectado1 = conectado1 .or. ((j1   == j2  ).and.(k1 == k2  ).and.(i1 == i2  ))
                           conectado1 = conectado1 .or. ((j1   == j2  ).and.(k1 == k2+1).and.(i1 == i2  ))
                           !
                           conectado2 = conectado2 .or. ((j1 +1== j2  ).and.(k1 == k2  ).and.(i1 == i2  ))
                           conectado2 = conectado2 .or. ((j1 +1== j2  ).and.(k1 == k2+1).and.(i1 == i2  ))
                         case (iEx)
                           conectado1 = conectado1 .or. ((j1   == j2  ).and.(k1 == k2  ).and.(i1 == i2  ))
                           conectado1 = conectado1 .or. ((j1   == j2  ).and.(k1 == k2  ).and.(i1 == i2+1))
                           !
                           conectado2 = conectado2 .or. ((j1 +1== j2  ).and.(k1 == k2  ).and.(i1 == i2  ))
                           conectado2 = conectado2 .or. ((j1 +1== j2  ).and.(k1 == k2  ).and.(i1 == i2+1))
                        end select
                      case (iEz)
                        select case (whatfield2)
                         case (iEz)
                           conectado1 = conectado1 .or. ((k1   == k2+1).and.(i1 == i2  ).and.(j1 == j2  ))
                           conectado2 = conectado2 .or. ((k1 +1== k2  ).and.(i1 == i2  ).and.(j1 == j2  ))
                         case (iEx)
                           conectado1 = conectado1 .or. ((k1   == k2  ).and.(i1 == i2  ).and.(j1 == j2  ))
                           conectado1 = conectado1 .or. ((k1   == k2  ).and.(i1 == i2+1).and.(j1 == j2  ))
                           !
                           conectado2 = conectado2 .or. ((k1 +1== k2  ).and.(i1 == i2  ).and.(j1 == j2  ))
                           conectado2 = conectado2 .or. ((k1 +1== k2  ).and.(i1 == i2+1).and.(j1 == j2  ))
                         case (iEy)
                           conectado1 = conectado1 .or. ((k1   == k2  ).and.(i1 == i2  ).and.(j1 == j2  ))
                           conectado1 = conectado1 .or. ((k1   == k2  ).and.(i1 == i2  ).and.(j1 == j2+1))
                           !
                           conectado2 = conectado2 .or. ((k1 +1== k2  ).and.(i1 == i2  ).and.(j1 == j2  ))
                           conectado2 = conectado2 .or. ((k1 +1== k2  ).and.(i1 == i2  ).and.(j1 == j2+1))
                        end select
                     end select

                     conectado=conectado1.and.conectado2
                     if (conectado) then
                        exit buskk
                     endif
                  end do buskk

                  if (connectendings) then
                     if ((.not.sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%Is_LeftEnd).and. &
                     (.not.sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%Is_RightEnd)) then
                        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%IsEnd_norLeft_norRight = sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%IsEnd_norLeft_norRight .or. &
                        ((.not.conectado).and.conectado2).and. &
                        (.not.(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%HasParallel_LeftEnd .or. sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%HasSeries_LeftEnd.or. &
                                                                                           sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%HasAbsorbing_LeftEnd))
                     endif
                     if ((.not.sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%Is_LeftEnd).and. &
                     (.not.sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%Is_RightEnd)) then
                        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%IsEnd_norLeft_norRight = sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%IsEnd_norLeft_norRight.or. &
                        ((.not.conectado).and.conectado1).and. &
                        (.not.(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%HasParallel_RightEnd .or. sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%HasSeries_RightEnd .or. &
                                                                                           sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%HasAbsorbing_RightEnd))
                     endif
                  endif
                  sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%IsEnd_norLeft_norRight = &
                  (sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%IsEnd_norLeft_norRight) .and. &
                  (.not.conectado).and.(conectado1.or.conectado2) .and. &
                  (.not.(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%Is_LeftEnd.or. &
                  sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%Is_RightEnd)) !si hay mas de uno este se pone a .true.
                  !detecta cual es el extremo libre
                  if ((.not.conectado).and.conectado1) then
                     select case (whatfield)
                      case (iEx)
                        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%ilibre = i1 +1
                        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%jlibre = j1
                        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%klibre = k1
                      case (iEy)
                        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%ilibre = i1
                        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%jlibre = j1 +1
                        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%klibre = k1
                      case (iEz)
                        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%ilibre = i1
                        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%jlibre = j1
                        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%klibre = k1 +1
                     end select
                  elseif ((.not.conectado).and.conectado2) then
                     sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%ilibre = i1
                     sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%jlibre = j1
                     sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%klibre = k1
                  endif

                  !caso especial
                  if (sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%numsegmentos == 1) then
                     sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%IsEnd_norLeft_norRight   = .false.
                     sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%Is_LeftEnd              = .true.
                     sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%Is_RightEnd              = .true.
                     select case (whatfield)
                      case (iEx)
                        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%ilibre = i1 +1
                        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%jlibre = j1
                        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%klibre = k1
                      case (iEy)
                        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%ilibre = i1
                        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%jlibre = j1 +1
                        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%klibre = k1
                      case (iEz)
                        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%ilibre = i1
                        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%jlibre = j1
                        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%klibre = k1 +1
                     end select
                  endif
                  !check for intermediate RLC error

                  if ( (sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%ilibre == -1).or. &
                  (sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%jlibre == -1).or. &
                  (sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%klibre == -1)) then
                     if (.not.conectado) then
                        write (buff,'(a,i7,3I7,a)')  'wir0_BUGGYERROR: Non-Intermediate multi-segment WIRE. ',origIndex,i1,j1,k1,' '//dir(whatfield)
                        if ((k1 >= ZI).and.(k1 <= ZE)) call WarnErrReport(buff,.true.)
                     endif
                     if ( ((sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%Is_RightEnd).and. &
                     (sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%HasParallel_RightEnd .or. sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%HasSeries_RightEnd)).or. &
                     ((sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%Is_LeftEnd).and. &
                     (sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%HasParallel_LeftEnd .or. sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%HasSeries_LeftEnd .or. &
                                                                                     sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%HasAbsorbing_LeftEnd)) ) then
                        if (conectado) then
                           sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%Is_RightEnd = .false.
                           sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%Is_LeftEnd = .false.
                           write (buff,'(a,i7,3I7,a)')  'wir0_WARNING: Intermediate segment with RLC. Neglecting RLC ',origIndex,i1,j1,k1,' '//dir(whatfield)
                           if ((k1 >= ZI).and.(k1 <= ZE)) call WarnErrReport(buff)
                        else
                           write (buff,'(a,i7,3I7,a)')  'wir0_BUGGYERROR: Non-Intermediate multi-segment WIRE with RLC. ',origIndex,i1,j1,k1,' '//dir(whatfield)
                           if ((k1 >= ZI).and.(k1 <= ZE)) call WarnErrReport(buff,.true.)
                        endif
                     endif
                  endif
                  !
               else !del strictOLD

                  sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%ilibre = -1
                  sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%jlibre = -1
                  sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%klibre = -1
                  sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%IsEnd_norLeft_norRight=.false. !irrelevante en strictOLD
                  i1=       sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%i
                  j1=       sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%j
                  k1=       sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%k
                  whatfield=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%ori
                  origindex= sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%origIndex
                  !
                  if (sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%Is_LeftEnd) then
                     dummy1=1
                     dummyfin=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%numsegmentos
                  elseif (sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%Is_RightEnd) then
                     dummy1=-1
                     dummyfin=1
                  endif
                  if ( ((sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%Is_LeftEnd).or.(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%Is_RightEnd)) ) THEN
                     IF (sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%numsegmentos/=1) then
                        dummy2=-1
                        buscakk2: do iwjjj=iwj+dummy1,dummyfin,dummy1 !atras o adelante
                           i2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%i
                           j2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%j
                           k2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%k
                           whatfield2=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%ori
                           if ((i1==i2).and.(j1==j2).and.(k1==k2).and.(whatfield==whatfield2)) then
                              dummy2=-dummy2 !detecta numero de rabitos para o impar aunque yo luego en las uniones solo trato 2 rabitos como mucho
                           else
                              continue
                              exit buscakk2
                           endif
                        end do buscakk2
                        if (dummy2==-1) then
                           dummy3=0
                        else
                           dummy3=1
                        endif

                        !
                        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%ilibre = i1
                        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%jlibre = j1
                        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%klibre = k1
                        if (whatfield2==whatfield) then
                           if (abs(i1-i2)+abs(j1-j2)+abs(k1-k2)>1) then
                              write (buff,'(a,i7,3I7,a)')  'wir0_ERROR: strictOLD LeftEnd/RightEnd segment disconnected.', origIndex,i1,j1,k1,' '//dir(whatfield)
                              if ((k1 >= ZI).and.(k1 <= ZE)) call WarnErrReport(buff,.true.)
                           endif
                           if (i1>i2) then
                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%ilibre = i1+1-dummy3
                           elseif (i1<i2) then
                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%ilibre = i1+dummy3
                           endif
                           if (j1>j2) then
                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%jlibre = j1+1-dummy3
                           elseif (j1<j2) then
                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%jlibre = j1+dummy3
                           endif
                           if (k1>k2) then
                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%klibre = k1+1-dummy3
                           elseif (k1<k2) then
                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%klibre = k1+dummy3
                           endif
                        else
                           if (abs(i1-i2)+abs(j1-j2)+abs(k1-k2)>2) then
                              write (buff,'(a,i7,3I7,a)')  'wir0_ERROR: strictOLD LeftEnd/RightEnd segment disconnected.', origIndex,i1,j1,k1,' '//dir(whatfield)
                              if ((k1 >= ZI).and.(k1 <= ZE)) call WarnErrReport(buff,.true.)
                           endif
                           select case (whatfield)
                            case (iEx)
                              if (i2 == i1) sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%ilibre = i1+1
                            case (iEy)
                              if (j2 == j1) sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%jlibre = j1+1
                            case (iEz)
                              if (k2 == k1) sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%klibre = k1+1
                           end select
                        endif
                     ELSE !DEL NUMERO SEGMENTOS '2014 NO PORTADO A !CHECK
                        select case (whatfield)
                         case (iEx)
                           sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%ilibre = i1+1
                           sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%jlibre = j1
                           sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%klibre = k1
                         case (iEy)
                           sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%ilibre = i1
                           sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%jlibre = j1+1
                           sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%klibre = k1
                         case (iEz)
                           sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%ilibre = i1
                           sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%jlibre = j1
                           sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%klibre = k1+1
                        end select
                     endif !DEL NUMERO SEGMENTOS '2014 NO PORTADO A !CHECK
                  ENDIF !Del Left_End Right_End
                  !!!!!!!!!!!!!!!!!!!
               endif !del strictOLD

            end do
         end do


         !preprocesa para eliminar multirabos luego se utiliza en repetido

         if (strictOLD) then
            do iwi=1,HWires%NumDifferentWires
               do iwj=1,sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%numsegmentos
                  if (.NOT.sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%multirabo) then
                     multirabos=1
                     i1=       sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%i
                     j1=       sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%j
                     k1=       sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%k
                     whatfield=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%ori
                     ORIGINDEX=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%origindex
                     !precontaje
                     buscarabos: do iwjjj=iwj+1,sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%numsegmentos !el Right_End aunque no se tape si debe detectarse a efectos par/impar
                        i2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%i
                        j2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%j
                        k2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%k
                        whatfield2=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%ori
                        ORIGINDEX2=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%origindex
                        if ((i1==i2).and.(j1==j2).and.(k1==k2).and.(whatfield==whatfield2)) then
                           multirabos=multirabos+1
                        else
                           primernorabo=origindex2
                           Jprimernorabo=iwjjj
                           exit buscarabos
                        endif
                     end do buscarabos
                     !machaca rabos
                     if (multirabos/=1) then
                        !
                        taparabos: do iwjjj=iwj+(2-mod(multirabos,2)),sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%numsegmentos-1 !el Right_End no debe taparse ya se tapa el de dentro en el otro bucle
                           i2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%i
                           j2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%j
                           k2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%k
                           whatfield2=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%ori
                           ORIGINDEX2=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%origindex
                           if ((i1==i2).and.(j1==j2).and.(k1==k2).and.(whatfield==whatfield2)) then
                              write (buff,'(a,i7,3I7,a,a,i7)')  'wir0_WARNING: strictOLD Redundannt zig-zag rabito, will be eliminated to Mod(2)', origIndex2,i2,j2,k2,dir(whatfield2), &
                              'by segment ',primernorabo
                              if ((k1 >= ZI).and.(k1 <= ZE)) call WarnErrReport(buff)
                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%multirabo=.true.
                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%multiraboDE=primernorabo
                           else
                              exit taparabos
                           endif
                        end do taparabos
                     endif !del multirabos no nulo
                  ENDIF
               end do
            end do
            !ahora la vuelta
            do iwi=1,HWires%NumDifferentWires
               do iwj=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%numsegmentos,1,-1
                  if (.NOT.sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%multirabo) then
                     multirabos=1
                     i1=       sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%i
                     j1=       sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%j
                     k1=       sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%k
                     whatfield=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%ori
                     ORIGINDEX=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%origindex
                     buscarabos2: do iwjjj=iwj-1,1,-1 !el Left_End aunque no se tape si debe detectarse a efectos par/impar
                        i2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%i
                        j2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%j
                        k2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%k
                        whatfield2=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%ori
                        ORIGINDEX2=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%origindex
                        if ((i1==i2).and.(j1==j2).and.(k1==k2).and.(whatfield==whatfield2)) then
                           multirabos=multirabos+1
                        else
                           primernorabo=origindex2
                           Jprimernorabo=iwjjj
                           exit buscarabos2
                        endif
                     end do buscarabos2

                     !machaca rabos
                     if (multirabos/=1) then
                        !
                        taparabos2: do iwjjj=iwj-(2-mod(multirabos,2)),2,-1 !el Left_End no debe taparse
                           i2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%i
                           j2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%j
                           k2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%k
                           whatfield2=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%ori
                           ORIGINDEX2=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%origindex
                           if ((i1==i2).and.(j1==j2).and.(k1==k2).and.(whatfield==whatfield2)) then
                              !faltaria eliminar sondas en multirabos 7/2/14
                              write (buff,'(a,i7,3I7,a,a,i7)')  'wir0_WARNING: strictOLD Redundannt zig-zag rabito, will be eliminated to Mod(2)', origIndex2,i2,j2,k2,dir(whatfield2), &
                              'by segment ',primernorabo
                              if ((k1 >= ZI).and.(k1 <= ZE)) call WarnErrReport(buff)
                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%multirabo=.true.
                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%multiraboDE=primernorabo
                           else
                              exit taparabos2
                           endif
                        end do taparabos2
                     endif !del multirabos no nulo
                  ENDIF
               end do
            end do


            !!!!!!!!!!!!!!!!
            !segunda pasada para procesar TAPARRABOS
            !!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!

            if (TAPARRABOS) then
               do iwi=1,HWires%NumDifferentWires
                  do iwj=1,sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%numsegmentos
                     if (.NOT.sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%multirabo) then
                        multirabos=1
                        i1=       sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%i
                        j1=       sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%j
                        k1=       sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%k
                        whatfield=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%ori
                        ORIGINDEX=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%origindex
                        !precontaje
                        buscarabos6: do iwjjj=iwj+1,sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%numsegmentos !el Right_End aunque no se tape si debe detectarse a efectos par/impar
                           i2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%i
                           j2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%j
                           k2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%k
                           whatfield2=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%ori
                           ORIGINDEX2=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%origindex
                           if ((i1==i2).and.(j1==j2).and.(k1==k2).and.(whatfield==whatfield2)) then
                              multirabos=multirabos+1
                           else
                              primernorabo=origindex2
                              Jprimernorabo=IWjjj
                              exit buscarabos6
                           endif
                        end do buscarabos6
                        !machaca rabos
                        if ((mod(multirabos,2)/=1).and.(Jprimernorabo /= sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%numsegmentos)) then !no al RightEnd
                           if (sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj          )%Is_LeftEnd) then
                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj          )%Is_LeftEnd=.false.
                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%Is_LeftEnd=.true. !pasa caracter Left_End al primernorabo
                              !!!tocado esto por el problema de gra_powerline_simple.nfde 190916 el rabito se quedaba con el libre mal computado. Los ilibre,jlbre,klibre para Left_End son el primer punto directamente
                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%ilibre=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%i  !!!!sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%ilibre
                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%jlibre=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%j  !!!!sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%jlibre
                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%klibre=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%k  !!!!sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%klibre
!!!y retocado 200117 por el problema del ntc1 de la demo de getafe que tambien calculaba mal el libre (no tiene que ser por guevos el primer punto, habra que ver como se conecta con el siguiente!!!)
                                 if ((sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%i)==(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo+1)%i).and. &
                                     (sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%j)==(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo+1)%j).and. &
                                     (sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%k)==(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo+1)%k)) then
                                     select case (sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%ori)
                                     case (1)
                                        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%ilibre=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%ilibre+1
                                     case (2)
                                        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%jlibre=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%jlibre+1
                                     case (3)
                                        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%klibre=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%klibre+1
                                     end select
                                endif
                           endif
                           !
                           if ((.NOT.sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%Is_LeftEnd).AND.(.NOT.sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%Is_RightEnd)) then
                              i2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%i
                              j2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%j
                              k2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%k
                              whatfield2=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%ori
                              ORIGINDEX2=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%origindex
                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%multirabo=.true.
                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%multiraboDE=primernorabo !el ultimo cazado es el primer NOrabo
                              write (buff,'(a,i7,3I7,a,a,i7)')  'wir0_WARNING: strictOLD and taparrabos redundannt zig-zag rabito, will be ALSO eliminated ', origIndex2,i2,j2,k2,dir(whatfield2), &
                              'by segment ',primernorabo
                              if ((k1 >= ZI).and.(k1 <= ZE)) call WarnErrReport(buff)
                           ENDIF
                           !
                           if  ((.NOT.sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj+1)%Is_LeftEnd).AND.(.NOT.sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj+1)%Is_RightEnd)) then
                              i2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj+1)%i
                              j2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj+1)%j
                              k2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj+1)%k
                              whatfield2=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj+1)%ori
                              ORIGINDEX2=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj+1)%origindex
                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj+1)%multirabo=.true.
                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj+1)%multiraboDE=primernorabo !el ultimo cazado es el primer NOrabo
                              write (buff,'(a,i7,3I7,a,a,i7)')  'wir0_WARNING: strictOLD and taparrabos redundannt zig-zag rabito, will be ALSO eliminated ', origIndex2,i2,j2,k2,dir(whatfield2), &
                              'by segment ',primernorabo
                              if ((k1 >= ZI).and.(k1 <= ZE)) call WarnErrReport(buff)
                           endif
                        endif !del multirabos no nulo
                     ENDIF
                  end do
               end do
               !ahora la vuelta
               do iwi=1,HWires%NumDifferentWires
                  do iwj=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%numsegmentos,1,-1
                     if (.NOT.sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%multirabo) then
                        multirabos=1
                        i1=       sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%i
                        j1=       sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%j
                        k1=       sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%k
                        whatfield=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%ori
                        ORIGINDEX=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%origindex
                        buscarabos7: do iwjjj=iwj-1,1,-1 !el Left_End aunque no se tape si debe detectarse a efectos par/impar
                           i2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%i
                           j2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%j
                           k2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%k
                           whatfield2=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%ori
                           ORIGINDEX2=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%origindex
                           if ((i1==i2).and.(j1==j2).and.(k1==k2).and.(whatfield==whatfield2)) then
                              multirabos=multirabos+1
                           else
                              primernorabo=origindex2
                              Jprimernorabo=iwjjj
                              exit buscarabos7
                           endif
                        end do buscarabos7

                        !machaca rabos
                        if ((mod(multirabos,2)/=1).and.(Jprimernorabo /= 1)) then   !no el Left_End
                           if (sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj          )%Is_RightEnd) then
                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj          )%Is_RightEnd=.false.
                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%Is_RightEnd=.true. !pasa caracter RightEnd al primernorabo
                              !!!tocado esto por el problema de gra_powerline_simple.nfde 190916 el rabito se quedaba con el libre mal computado. Los ilibre,jlbre,klibre para Left_End son el primer punto directamente
                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%ilibre=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%i  !!!!sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%ilibre
                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%jlibre=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%j  !!!!sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%jlibre
                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%klibre=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%k  !!!!sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%klibre
!!!y retocado 200117 por el problema del ntc1 de la demo de getafe que tambien calculaba mal el libre (no tiene que ser por guevos el primer punto, habra que ver como se conecta con el siguiente!!!)
                                 if ((sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%i)==(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo-1)%i).and. &
                                     (sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%j)==(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo-1)%j).and. &
                                     (sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%k)==(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo-1)%k)) then
                                     select case (sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%ori)
                                     case (1)
                                        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%ilibre=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%ilibre+1
                                     case (2)
                                        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%jlibre=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%jlibre+1
                                     case (3)
                                        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%klibre=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%klibre+1
                                     end select
                                endif
                           
                           endif
                           !
                           if ((.NOT.sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%Is_LeftEnd).AND.(.NOT.sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%Is_RightEnd)) then
                              i2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%i
                              j2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%j
                              k2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%k
                              whatfield2=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%ori
                              ORIGINDEX2=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%origindex
                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%multirabo=.true.
                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%multiraboDE=primernorabo !el ultimo cazado es el primer NOrabo
                              write (buff,'(a,i7,3I7,a,a,i7)')  'wir0_WARNING: strictOLD and taparrabos redundannt zig-zag rabito, will be ALSO eliminated ', origIndex2,i2,j2,k2,dir(whatfield2), &
                              'by segment ',primernorabo
                              if ((k1 >= ZI).and.(k1 <= ZE)) call WarnErrReport(buff)
                           endif
                           !
                           if  ((.NOT.sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj-1)%Is_LeftEnd).AND.(.NOT.sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj-1)%Is_RightEnd)) then
                              i2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj-1)%i
                              j2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj-1)%j
                              k2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj-1)%k
                              whatfield2=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj-1)%ori
                              ORIGINDEX2=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj-1)%origindex
                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj-1)%multirabo=.true.
                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj-1)%multiraboDE=primernorabo !el ultimo cazado es el primer NOrabo
                              write (buff,'(a,i7,3I7,a,a,i7)')  'wir0_WARNING: strictOLD and taparrabos redundannt zig-zag rabito, will be ALSO eliminated ', origIndex2,i2,j2,k2,dir(whatfield2), &
                              'by segment ',primernorabo
                              if ((k1 >= ZI).and.(k1 <= ZE)) call WarnErrReport(buff)
                           endif
                        endif !del multirabos no nulo
                     ENDIF
                  end do
               end do
            endif !del TAPARRABOS
            !
         endif !del if strictOLD

         !chequeo de errores buggy multirabo
         do iwi=1,HWires%NumDifferentWires
            do iwj=1,sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%numsegmentos
               if ((sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%multirabo).and. &
               ((sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%Is_LeftEnd).or.(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%Is_RightEnd))) then
                  i1=       sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%i
                  j1=       sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%j
                  k1=       sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%k
                  whatfield=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%ori
                  ORIGINDEX=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%origindex
                  write (buff,'(a,i7,3I7,a)')  'wir0_BuggyERROR: strictOLD LeftEnd/RightEnd cannot be multirabo ', origIndex,i1,j1,k1,' '//dir(whatfield)
                  if ((k1 >= ZI).and.(k1 <= ZE)) call WarnErrReport(buff,.true.)
               endif
            end do
         end do





         !!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !report them
         do iwi=1,HWires%NumDifferentWires
            do iwj=1,      sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%numsegmentos
               i1=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%i
               j1=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%j
               k1=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%k
               i1libre=   sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%ilibre
               j1libre=   sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%jlibre
               k1libre=   sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%klibre
               whatfield= sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%ori
               origindex= sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%origIndex
               LeftEnd_index=  sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%LeftEnd
               RightEnd_index=  sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%RightEnd
               !
               if (sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%Is_LeftEnd.and.sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%Is_RightEnd) then
                  write (buff,'(a,4I7,a,3I7,a,2i7)')  'wir0_INFO: Ending segment (LeftEnd_and_Right)',origIndex,i1,j1,k1,'-', &
                  i1libre,j1libre,k1libre,' '//dir(whatfield),LeftEnd_index,RightEnd_index
                  if ((k1 >= ZI).and.(k1 <= ZE).and.verbose) call WarnErrReport(buff)
               elseif (sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%Is_LeftEnd) then
                  write (buff,'(a,4I7,a,3I7,a,i7)')  'wir0_INFO: Ending segment (LeftEnd   )',origIndex,i1,j1,k1,'-', &
                  i1libre,j1libre,k1libre,' '//dir(whatfield),LeftEnd_index
                  if ((k1 >= ZI).and.(k1 <= ZE).and.verbose) call WarnErrReport(buff)
               elseif (sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%Is_RightEnd) then
                  write (buff,'(a,4I7,a,3I7,a,i7)')  'wir0_INFO: Ending segment (RightEnd   )',origIndex,i1,j1,k1,'-', &
                  i1libre,j1libre,k1libre,' '//dir(whatfield),RightEnd_index
                  if ((k1 >= ZI).and.(k1 <= ZE).and.verbose) call WarnErrReport(buff)
               endif
               if (sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%IsEnd_norLeft_norRight) then
                  if (connectendings) then
                     write (buff,'(a,4I7,a,3I7,a)')  'wir0_INFO: Ending segment (other )',origIndex,i1,j1,k1,'-', &
                     i1libre,j1libre,k1libre,' '//dir(whatfield)
                     if ((k1 >= ZI).and.(k1 <= ZE).and.verbose) call WarnErrReport(buff)
                  else
                     sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%IsEnd_norLeft_norRight =.false.
                     write (buff,'(a,4I7,a,3I7,a)')  'wir0_WARNING: Resetting Ending segment (other ) to NON-ENDING',origIndex,i1,j1,k1,'-', &
                     i1libre,j1libre,k1libre,' '//dir(whatfield)
                     if ((k1 >  ZI).and.(k1 <= ZE).and.(whatfield /= iEz)) call WarnErrReport(buff)
                     if ((k1 >= ZI).and.(k1 <= ZE).and.(whatfield == iEz)) call WarnErrReport(buff)
                  endif
               endif
            end do
         end do

         !Segment pre-counting incluyendo deteccion de repetidos

         do iwi=1,HWires%NumDifferentWires
            do iwj=1,sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%numsegmentos
               sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%repetido=.false.
            end do
         enddo
         !
         conta=0
         do iwi=1,HWires%NumDifferentWires
            do iwj=1,sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%numsegmentos
               if (((.not.sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%repetido).or.strictOLD).and.(.not.sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%multirabo)) then
                  i1=       sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%i
                  j1=       sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%j
                  k1=       sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%k
                  whatfield=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%ori
                  !
                  do iwjjj=iwj+1,sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%numsegmentos
                     i2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%i
                     j2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%j
                     k2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%k
                     whatfield2=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%ori
                     repetido = (i1 == i2).and.(j1 == j2).and.(k1 == k2).and.(whatfield == whatfield2)
                     if (repetido) then
                        if     (.not.(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%Is_LeftEnd .or. &
                        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%Is_RightEnd )) THEN

                           sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%repetido=repetido.or. &
                           sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%repetido
                        elseif (.not.(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%Is_LeftEnd .or. &
                        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%Is_RightEnd )) THEN

                           sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%repetido=repetido.or. &
                           sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%repetido
                        else
                           !aviso pero tomo una decision. md 260213 a veces lo duplica en principio y final!!!!!!

                           if ( ((abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Parallel_R_RightEnd) < 1.0e-12_RKIND_wires).and. &
                           (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Series_R_RightEnd  ) < 1.0e-12_RKIND_wires).and. &
                           (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Parallel_C_RightEnd) < 1.0e-12_RKIND_wires).and. &
                           (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Series_C_RightEnd  ) > 1.0e7_RKIND_wires).and. &
                           (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Parallel_L_RightEnd) < 1.0e-12_RKIND_wires).and. &
                           (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Series_L_RightEnd  ) < 1.0e-12_RKIND_wires)).and. &
                           !
                           ((abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Parallel_R_LeftEnd) >= 1.0e-12_RKIND_wires).or. &
                           (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Series_R_LeftEnd  ) >= 1.0e-12_RKIND_wires).or. &
                           (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Parallel_C_LeftEnd) >= 1.0e-12_RKIND_wires).or. &
                           (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Series_C_LeftEnd  ) <= 1.0e7_RKIND_wires).or. &
                           (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Parallel_L_LeftEnd) >= 1.0e-12_RKIND_wires).or. &
                           (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Series_L_LeftEnd  ) >= 1.0e-12_RKIND_wires)) ) then
                              if (sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%Is_RightEnd) then
                                 sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%repetido=.true.
                                 if (strictOLD) then
                                    write (buff,'(a,2i7,3i7)')  'wir0_INFO: Duplicate terminal LeftEnd and RightEnd Parallel segment. Keeping both', &
                                    sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%origindex, &
                                    sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%origindex,i1,j1,k1
                                 else
                                    write (buff,'(a,2i7,3i7)')  'wir0_INFO: Duplicate terminal LeftEnd and RightEnd Parallel segment. Removing the second one (no RightEnd RLC)', &
                                    sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%origindex, &
                                    sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%origindex,i1,j1,k1
                                 endif
                                 if ((k1 >  ZI).and.(k1 <= ZE).and.(whatfield /= iEz).and.verbose) call WarnErrReport(buff)
                                 if ((k1 >= ZI).and.(k1 <= ZE).and.(whatfield == iEz).and.verbose) call WarnErrReport(buff)
                              elseif (sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%Is_RightEnd) then
                                 sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%repetido=.true.
                                 if (strictOLD) then
                                    write (buff,'(a,2i7,3i7)')  'wir0_INFO: Duplicate terminal LeftEnd and RightEnd Parallel segment. Keeping both', &
                                    sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%origindex, &
                                    sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%origindex,i1,j1,k1
                                 else
                                    write (buff,'(a,2i7,3i7)')  'wir0_INFO: Duplicate terminal LeftEnd and RightEnd Parallel segment. Removing the first one (no RightEnd RLC)', &
                                    sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%origindex, &
                                    sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%origindex,i1,j1,k1
                                 endif
                                 if ((k1 >  ZI).and.(k1 <= ZE).and.(whatfield /= iEz).and.verbose) call WarnErrReport(buff)
                                 if ((k1 >= ZI).and.(k1 <= ZE).and.(whatfield == iEz).and.verbose) call WarnErrReport(buff)
                              endif
                           endif

                           if ( ((abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Parallel_R_RightEnd) >= 1.0e-12_RKIND_wires).or. &
                           (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Series_R_RightEnd  ) >= 1.0e-12_RKIND_wires).or. &
                           (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Parallel_C_RightEnd) >= 1.0e-12_RKIND_wires).or. &
                           (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Series_C_RightEnd  ) <= 1.0e7_RKIND_wires).or. &
                           (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Parallel_L_RightEnd) >= 1.0e-12_RKIND_wires).or. &
                           (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Series_L_RightEnd  ) >= 1.0e-12_RKIND_wires)).and. &
                           !
                           ((abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Parallel_R_LeftEnd) < 1.0e-12_RKIND_wires).and. &
                           (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Series_R_LeftEnd  ) < 1.0e-12_RKIND_wires).and. &
                           (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Parallel_C_LeftEnd) < 1.0e-12_RKIND_wires).and. &
                           (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Series_C_LeftEnd  ) > 1.0e7_RKIND_wires).and. &
                           (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Parallel_L_LeftEnd) < 1.0e-12_RKIND_wires).and. &
                           (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Series_L_LeftEnd  ) < 1.0e-12_RKIND_wires)) ) then
                              if (sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%Is_LeftEnd) then
                                 sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%repetido=.true.
                                 if (strictOLD) then
                                    write (buff,'(a,2i7,3i7)')  'wir0_INFO: Duplicate terminal LeftEnd and RightEnd Parallel segment. Keeping both', &
                                    sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%origindex, &
                                    sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%origindex,i1,j1,k1
                                 else
                                    write (buff,'(a,2i7,3i7)')  'wir0_INFO: Duplicate terminal LeftEnd and RightEnd Parallel segment. Removing the second one (NO LeftEnd RLC, Only RightEnd RLC)', &
                                    sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%origindex, &
                                    sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%origindex,i1,j1,k1
                                 endif
                                 if ((k1 >  ZI).and.(k1 <= ZE).and.(whatfield /= iEz).and.verbose) call WarnErrReport(buff)
                                 if ((k1 >= ZI).and.(k1 <= ZE).and.(whatfield == iEz).and.verbose) call WarnErrReport(buff)
                              elseif (sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%Is_LeftEnd) then
                                 sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%repetido=.true.
                                 if (strictOLD) then
                                    write (buff,'(a,2i7,3i7)')  'wir0_INFO: Duplicate terminal LeftEnd and RightEnd Parallel segment. Keeping both', &
                                    sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%origindex, &
                                    sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%origindex,i1,j1,k1
                                 else
                                    write (buff,'(a,2i7,3i7)')  'wir0_INFO: Duplicate terminal LeftEnd and RightEnd Parallel segment. Removing the first one (NO LeftEnd RLC, Only RightEnd RLC)', &
                                    sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%origindex, &
                                    sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%origindex,i1,j1,k1
                                 endif
                                 if ((k1 >  ZI).and.(k1 <= ZE).and.(whatfield /= iEz).and.verbose) call WarnErrReport(buff)
                                 if ((k1 >= ZI).and.(k1 <= ZE).and.(whatfield == iEz).and.verbose) call WarnErrReport(buff)
                              endif
                           endif
                           !
                           if ( ((abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Parallel_R_LeftEnd) >= 1.0e-12_RKIND_wires).or. &
                           (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Series_R_LeftEnd  ) >= 1.0e-12_RKIND_wires).or. &
                           (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Parallel_C_LeftEnd) >= 1.0e-12_RKIND_wires).or. &
                           (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Series_C_LeftEnd  ) <= 1.0e7_RKIND_wires).or. &
                           (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Parallel_L_LeftEnd) >= 1.0e-12_RKIND_wires).or. &
                           (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Series_L_LeftEnd  ) >= 1.0e-12_RKIND_wires)).AND. &
                           !
                           ((abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Parallel_R_RightEnd) >= 1.0e-12_RKIND_wires).or. &
                           (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Series_R_RightEnd  ) >= 1.0e-12_RKIND_wires).or. &
                           (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Parallel_C_RightEnd) >= 1.0e-12_RKIND_wires).or. &
                           (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Series_C_RightEnd  ) <= 1.0e7_RKIND_wires).or. &
                           (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Parallel_L_RightEnd) >= 1.0e-12_RKIND_wires).or. &
                           (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Series_L_RightEnd  ) >= 1.0e-12_RKIND_wires)) ) then

                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%repetido=repetido.or. &
                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%repetido

                              if (  (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Parallel_R_RightEnd -    &
                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Parallel_R_LeftEnd) < 1.0e-12_RKIND_wires).and. &
                              (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Series_R_RightEnd -    &
                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Series_R_LeftEnd  ) < 1.0e-12_RKIND_wires).and. &
                              (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Parallel_C_RightEnd -    &
                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Parallel_C_LeftEnd) < 1.0e-12_RKIND_wires).and. &
                              (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Series_C_RightEnd -    &
                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Series_C_LeftEnd  ) < 1.0e-12_RKIND_wires).and. &
                              (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Parallel_L_RightEnd -    &
                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Parallel_L_LeftEnd) < 1.0e-12_RKIND_wires).and. &
                              (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Series_L_RightEnd -    &
                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Series_L_LeftEnd  ) < 1.0e-12_RKIND_wires)) then
                                 if (strictOLD) then
                                    write (buff,'(a,2i7,3i7)')  'wir0_INFO: Duplicate terminal LeftEnd and RightEnd Parallel segment with the same RLC. Keeping both', &
                                    sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%origindex, &
                                    sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%origindex,i1,j1,k1
                                 else
                                    write (buff,'(a,2i7,3i7)')  'wir0_INFO: Duplicate terminal LeftEnd and RightEnd Parallel segment with the same RLC. Will remove the second one', &
                                    sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%origindex, &
                                    sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%origindex,i1,j1,k1
                                 endif
                                 if ((k1 >  ZI).and.(k1 <= ZE).and.(whatfield /= iEz).and.verbose) call WarnErrReport(buff)
                                 if ((k1 >= ZI).and.(k1 <= ZE).and.(whatfield == iEz).and.verbose) call WarnErrReport(buff)
                              else
                                 write (buff,'(a,2i7,3i7)')  'wir0_ERROR: Duplicate terminal LeftEnd and RightEnd Parallel segment with non-null different RLC', &
                                 sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%origindex, &
                                 sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%origindex,i1,j1,k1
                                 if ((k1 >  ZI).and.(k1 <= ZE).and.(whatfield /= iEz)) call WarnErrReport(buff,.true.)
                                 if ((k1 >= ZI).and.(k1 <= ZE).and.(whatfield == iEz)) call WarnErrReport(buff,.true.)
                              endif
                           ENDIF
                        endif
                     endif
                  end do
                  !
                  !clipping: in case of direct .nfde reading the PREPROCESSor has not clipped this data
                  if ((i1 >= sgg%Alloc(whatfield)%XI).and. &
                  (i1 <= sgg%Alloc(whatfield)%XE).and. &
                  (j1 >= sgg%Alloc(whatfield)%YI).and. &
                  (j1 <= sgg%Alloc(whatfield)%YE).and. &
                  (k1 >= sgg%Alloc(whatfield)%ZI).and. &
                  (k1 <= sgg%Alloc(whatfield)%ZE)) then
                     conta=conta+1
                  endif
               endif
            end do
         end do
         if (conta==0) therearewires=.false.
         if (.not.therearewires) Return
         !Report duplicated segments
         do iwi=1,HWires%NumDifferentWires
            do iwj=1,sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%numsegmentos
               i1=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%i
               j1=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%j
               k1=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%k
               whatfield =sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%ori
               origindex= sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%origIndex
               if ((sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%repetido).and.(.not.sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%multirabo)) then
                  if (strictOLD) then
                     write (buff,'(a,4i7,a)')  'wir0_WARNING: Keeping duplicate (Parallel) intra-WIRE segment', &
                     origindex,i1,j1,k1,' '//dir(whatfield)
                  else
                     write (buff,'(a,4i7,a)')  'wir0_WARNING: Removing duplicate (Parallel) intra-WIRE segment and voiding ASSOCIATED probes ', &
                     origindex,i1,j1,k1,' '//dir(whatfield)
                  endif
                  if ((k1 >= ZI).and.(k1 <= ZE).and.verbose) call WarnErrReport(buff)
               endif
            end do
         end do

         HWires%NumCurrentSegments=conta
         !inicializa ctes segmentos
         allocate (HWires%CurrentSegment(1 : HWires%NumCurrentSegments))
         do i1=1,HWires%NumCurrentSegments
            nullify (HWires%CurrentSegment(i1)%ChargePlus,HWires%CurrentSegment(i1)%ChargeMinus,  &
            HWires%CurrentSegment(i1)%TipoWire, &
            HWires%CurrentSegment(i1)%Efield_main2wire,HWires%CurrentSegment(i1)%Efield_wire2main)
            !
            HWires%CurrentSegment(i1)%R      = 0.0_RKIND_wires
            HWires%CurrentSegment(i1)%Resist      = 0.0_RKIND_wires
            HWires%CurrentSegment(i1)%Resist_devia      = 0.0_RKIND_wires
            HWires%CurrentSegment(i1)%C      = 0.0_RKIND_wires
            HWires%CurrentSegment(i1)%L      = 0.0_RKIND_wires
            HWires%CurrentSegment(i1)%Lintrinsic    = 0.0_RKIND_wires
            HWires%CurrentSegment(i1)%proc   =.false.
            HWires%CurrentSegment(i1)%NumParallel    =1
            HWires%CurrentSegment(i1)%origindex       =i1
            HWires%CurrentSegment(i1)%indexsegment    =i1
            HWires%CurrentSegment(i1)%currentpast         =0.0_RKIND_wires
            HWires%CurrentSegment(i1)%current         =0.0_RKIND_wires
            HWires%CurrentSegment(i1)%inv_Lind_acum            =0.0_RKIND_wires
            HWires%CurrentSegment(i1)%Lind_acum            =0.0_RKIND_wires
            HWires%CurrentSegment(i1)%HEUR_safety =0.0_RKIND_wires
            HWires%CurrentSegment(i1)%Lind            =0.0_RKIND_wires
            HWires%CurrentSegment(i1)%Lind_devia            =0.0_RKIND_wires
            !!! HWires%CurrentSegment(i1)%logRoverR0      =0.0_RKIND_wires
            HWires%CurrentSegment(i1)%delta           =0.0_RKIND_wires
            HWires%CurrentSegment(i1)%deltaTransv1    =0.0_RKIND_wires
            HWires%CurrentSegment(i1)%deltaTransv2    =0.0_RKIND_wires
            HWires%CurrentSegment(i1)%cte1            =0.0_RKIND_wires
            HWires%CurrentSegment(i1)%cte2            =0.0_RKIND_wires
            HWires%CurrentSegment(i1)%cte3            =0.0_RKIND_wires
            HWires%CurrentSegment(i1)%cte1_for_devia            =0.0_RKIND_wires
            HWires%CurrentSegment(i1)%cte2_for_devia            =0.0_RKIND_wires
            HWires%CurrentSegment(i1)%cte3_for_devia            =0.0_RKIND_wires
            HWires%CurrentSegment(i1)%cte5            =0.0_RKIND_wires
            HWires%CurrentSegment(i1)%ilibre               =-1
            HWires%CurrentSegment(i1)%jlibre               =-1
            HWires%CurrentSegment(i1)%klibre               =-1
            HWires%CurrentSegment(i1)%i               =-1
            HWires%CurrentSegment(i1)%j               =-1
            HWires%CurrentSegment(i1)%k               =-1
            HWires%CurrentSegment(i1)%tipofield       =-1
            HWires%CurrentSegment(i1)%IsPMC           =.false.
            HWires%CurrentSegment(i1)%orientadoalreves =.false.
            HWires%CurrentSegment(i1)%HasVsource      =.false.
            HWires%CurrentSegment(i1)%IsShielded      =.false.
            HWires%CurrentSegment(i1)%HasAbsorbing_RightEnd =.false.
            HWires%CurrentSegment(i1)%HasAbsorbing_LeftEnd =.false.
            HWires%CurrentSegment(i1)%HasParallel_RightEnd =.false.
            HWires%CurrentSegment(i1)%HasParallel_LeftEnd =.false.
            HWires%CurrentSegment(i1)%HasSeries_RightEnd   =.false.
            HWires%CurrentSegment(i1)%HasSeries_LeftEnd   =.false.
            HWires%CurrentSegment(i1)%IsEnd_norLeft_norRight   =.false.
            HWires%CurrentSegment(i1)%Is_LeftEnd   =.false.
            HWires%CurrentSegment(i1)%Is_RightEnd   =.false.
         end do

         !assign segment info
         conta=0
         do iwi=1,HWires%NumDifferentWires
            do iwj=1,sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%numsegmentos
               i1libre=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%ilibre
               j1libre=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%jlibre
               k1libre=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%klibre
               i1=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%i
               j1=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%j
               k1=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%k
               whatfield= sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%ori
               origindex= sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%origIndex
               IsEnd_norLeft_norRight=     sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%IsEnd_norLeft_norRight
               Is_LeftEnd=     sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%Is_LeftEnd
               Is_RightEnd=     sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%Is_RightEnd
               if (((.not.(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%repetido)).or.strictOLD).and.(.not.sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%multirabo)) then
                  !clipping: in case of direct .nfde reading the PREPROCESSor has not clipped this data
                  if ((i1 >= sgg%Alloc(whatfield)%XI).and. &
                  (i1 <= sgg%Alloc(whatfield)%XE).and. &
                  (j1 >= sgg%Alloc(whatfield)%YI).and. &
                  (j1 <= sgg%Alloc(whatfield)%YE).and. &
                  (k1 >= sgg%Alloc(whatfield)%ZI).and. &
                  (k1 <= sgg%Alloc(whatfield)%ZE)) then
                     conta=conta+1
                     HWires%CurrentSegment(conta)%IsEnd_norLeft_norRight=IsEnd_norLeft_norRight
                     HWires%CurrentSegment(conta)%Is_LeftEnd =Is_LeftEnd
                     HWires%CurrentSegment(conta)%Is_RightEnd =Is_RightEnd
                     HWires%CurrentSegment(conta)%origindex=origindex
                     HWires%CurrentSegment(conta)%tipofield=whatfield
                     HWires%CurrentSegment(conta)%ilibre=i1libre
                     HWires%CurrentSegment(conta)%jlibre=j1libre
                     HWires%CurrentSegment(conta)%klibre=k1libre
                     HWires%CurrentSegment(conta)%i=i1
                     HWires%CurrentSegment(conta)%j=j1
                     HWires%CurrentSegment(conta)%k=k1
                     HWires%CurrentSegment(conta)%ie=i1
                     HWires%CurrentSegment(conta)%je=j1
                     HWires%CurrentSegment(conta)%ke=k1
                     HWires%CurrentSegment(conta)%indexmed=HWires%WireTipoMedio(iwi)
                     HWires%CurrentSegment(conta)%TipoWire=>sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)
                     !!only for the observation sign to match (not used in this routine)
                     if (iwj < sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%numsegmentos) then
                        if (.not.strictOLD) then
                           HWires%CurrentSegment(conta)%orientadoalreves = &
                           (i1 > sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj+1)%i).or. &
                           (j1 > sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj+1)%j).or. &
                           (k1 > sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj+1)%k)
                        else
                           HWires%CurrentSegment(conta)%orientadoalreves =.false. !later corrected
                        endif
                     elseif     (iwj > 1 ) then
                        !only for the observation sign to match (not used in this routine)
                        if (.not.strictOLD) then
                           HWires%CurrentSegment(conta)%orientadoalreves = &
                           (i1 < sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj-1)%i).or. &
                           (j1 < sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj-1)%j).or. &
                           (k1 < sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj-1)%k)
                        else
                           HWires%CurrentSegment(conta)%orientadoalreves =.false. !later corrected
                        endif
                     endif
                     !
                     select case (HWires%CurrentSegment(conta)%tipofield)
                     case (iEx)
                        ! default
                        HWires%CurrentSegment(conta)%Efield_wire2main => Ex(i1,j1,k1) 
                        HWires%CurrentSegment(conta)%Efield_main2wire => Ex(i1,j1,k1) 
!!al final se reapuntan los del thickness
                        HWires%CurrentSegment(conta)%delta=1.0_RKIND_wires / Idxe(i1)    !ojo esto de los delta habra que corregirlo  para uniones
                        HWires%CurrentSegment(conta)%deltaTransv1=1.0_RKIND_wires / Idyh(j1)
                        if (k1 <= sgg%ALLOC(iEz)%ZE) then !esta corriente en el limite de los alloc nunca se precisa
                           !updateo por exceso para que haya una celda de solapamiento y la union topologica de
                           !hilos no se pierda por culpa de la puta particion MPI
                           HWires%CurrentSegment(conta)%deltaTransv2=1.0_RKIND_wires / Idzh(k1)
                        else !no se comete error alguno
                           !es solo para que el indice no reviente 2012
                           HWires%CurrentSegment(conta)%deltaTransv2=1.0_RKIND_wires / Idzh(k1-1)
                        endif
                        !dama
                        HWires%CurrentSegment(conta)%ie     = i1+1
                        HWires%CurrentSegment(conta)%x      = j1+0.25_RKIND_wires
                        HWires%CurrentSegment(conta)%y      = k1+0.25_RKIND_wires 
                        sggmiE = sggmiEx(i1,j1,k1); call deembed_peclossyconformal_segments(sggmiE); sggmiEx(i1,j1,k1)= sggmiE !por si se ha modificado !ojo agresivo 180220 !ojo cambiado aqui 170323 de sitio pq no se habian puesto los deltatrans y salia division por cero
                        !fin dama
                      case (iEy)
                        ! default
                        HWires%CurrentSegment(conta)%Efield_wire2main => Ey(i1,j1,k1) 
                        HWires%CurrentSegment(conta)%Efield_main2wire => Ey(i1,j1,k1)      
                        !
                        HWires%CurrentSegment(conta)%delta=1.0_RKIND_wires / Idye(j1)
                        if (k1 <= sgg%ALLOC(iEz)%ZE) then !esta corriente en el limite de los alloc nunca se precisa
                           !updateo por exceso para que haya una celda de solapamiento y la union topologica de
                           !hilos no se pierda por culpa de la puta particion MPI
                           HWires%CurrentSegment(conta)%deltaTransv1=1.0_RKIND_wires / Idzh(k1)
                        else
                           HWires%CurrentSegment(conta)%deltaTransv1=1.0_RKIND_wires / Idzh(k1-1) !no se comete error alguno
                           !es solo para que el indice no reviente 2012
                        endif
                        HWires%CurrentSegment(conta)%deltaTransv2=1.0_RKIND_wires / Idxh(i1)
                        !dama
                        HWires%CurrentSegment(conta)%je     = j1+1
                        HWires%CurrentSegment(conta)%x      = k1+0.25_RKIND_wires
                        HWires%CurrentSegment(conta)%y      = i1+0.25_RKIND_wires
                        sggmiE = sggmiEy(i1,j1,k1); call deembed_peclossyconformal_segments(sggmiE); sggmiEy(i1,j1,k1)= sggmiE !por si se ha modificado !ojo agresivo 180220 !ojo cambiado aqui 170323 de sitio pq no se habian puesto los deltatrans y salia division por cero
                        !fin dama
                      case (iEz)
                        ! default
                        HWires%CurrentSegment(conta)%Efield_wire2main => Ez(i1,j1,k1) 
                        HWires%CurrentSegment(conta)%Efield_main2wire => Ez(i1,j1,k1)  
                        !
                        HWires%CurrentSegment(conta)%delta=1.0_RKIND_wires / Idze(k1)
                        HWires%CurrentSegment(conta)%deltaTransv1=1.0_RKIND_wires / Idxh(i1)
                        HWires%CurrentSegment(conta)%deltaTransv2=1.0_RKIND_wires / Idyh(j1)
                        !dama
                        HWires%CurrentSegment(conta)%ke     = k1+1
                        HWires%CurrentSegment(conta)%x      = i1+0.25_RKIND_wires
                        HWires%CurrentSegment(conta)%y      = j1+0.25_RKIND_wires
                        sggmiE = sggmiEz(i1,j1,k1); call deembed_peclossyconformal_segments(sggmiE); sggmiEz(i1,j1,k1)= sggmiE !por si se ha modificado !ojo agresivo 180220 !ojo cambiado aqui 170323 de sitio pq no se habian puesto los deltatrans y salia division por cero
                        !fin dama
                     end select
                  endif
               endif !del repetido
            end do
         end do
 
!!fin niapa 171216
      !
      HWires%NumCurrentSegments=conta
      !

      !hacer agujeros

      do i1=1,HWires%NumCurrentSegments
         segmento=>HWires%CurrentSegment(i1)
         i=segmento%i
         j=segmento%j
         k=segmento%k
         whatfield= segmento%tipofield
         IsEnd_norLeft_norRight=segmento%IsEnd_norLeft_norRight
         Is_LeftEnd=segmento%Is_LeftEnd
         Is_RightEnd=segmento%Is_RightEnd
         if ((i > SINPML_fullsize(whatfield)%XI).and. &
         (i < SINPML_fullsize(whatfield)%XE).and. &
         (j > SINPML_fullsize(whatfield)%YI).and. &
         (j < SINPML_fullsize(whatfield)%YE).and. &
         (k > SINPML_fullsize(whatfield)%ZI).and. &
         (k < SINPML_fullsize(whatfield)%ZE)) then
            if (makeholes.and.(.not.IsEnd_norLeft_norRight).and.(.not.Is_LeftEnd).and.(.not.Is_RightEnd)) then
                if (size==0) then 
                    call stoponerror(layoutnumber,size,'Makeholes not available for MPI. Stoppping. ')
                endif
               select case (whatfield)
                case (iEx)
                  sggmiHx(i  ,j  ,k      ) = 1
                  sggmiHx(i  ,j-1,k      ) = 1
                  sggmiHx(i  ,j  ,k -1   ) = 1
                  sggmiHx(i  ,j-1,k -1   ) = 1
                  sggmiHx(i+1,j  ,k      ) = 1
                  sggmiHx(i+1,j-1,k      ) = 1
                  sggmiHx(i+1,j  ,k -1   ) = 1
                  sggmiHx(i+1,j-1,k -1   ) = 1
                  if (.not.sgg%med(sggmiEx(i  ,j  ,k      ))%Is%ThinWire) then
                     sggmiEx(i  ,j  ,k      ) = 1
                     write (buff,'(a,3i7,a,i7,a)')  'wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at ', &
                     i,j,k,' for WIRE-segment ',segmento%origIndex,' '//dir(whatfield)
                     call WarnErrReport(buff)
                  endif
                  if (.not.sgg%med(sggmiEy(i  ,j  ,k      ))%Is%ThinWire)  then
                     sggmiEy(i  ,j  ,k      ) = 1
                     write (buff,'(a,3i7,a,i7,a)')  'wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at ', &
                     i,j,k,' for WIRE-segment ',segmento%origIndex,' '//dir(whatfield)
                     call WarnErrReport(buff)
                  endif
                  if (.not.sgg%med(sggmiEy(i  ,j-1,k      ))%Is%ThinWire)  then
                     sggmiEy(i  ,j-1,k      ) = 1
                     write (buff,'(a,3i7,a,i7,a)')  'wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at ', &
                     i,j,k,' for WIRE-segment ',segmento%origIndex,' '//dir(whatfield)
                     call WarnErrReport(buff)
                  endif
                  if       ((k <=  sgg%alloc(iEz)%ZE).and.(k >= sgg%alloc(iEz)%ZI)) then
                     if (.not.sgg%med(sggmiEz(i  ,j  ,k      ))%Is%ThinWire)  then
                        sggmiEz(i  ,j  ,k      ) = 1
                        write (buff,'(a,3i7,a,i7,a)')  'wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at ', &
                        i,j,k,' for WIRE-segment ',segmento%origIndex,' '//dir(whatfield)
                        call WarnErrReport(buff)
                     endif
                  endif
                  if   ((k-1 <=  sgg%alloc(iEz)%ZE).and.(k -1 >=sgg%alloc(iEz)%ZI)) then
                     if (.not.sgg%med(sggmiEz(i  ,j  ,k-1    ))%Is%ThinWire)  then
                        sggmiEz(i  ,j  ,k - 1  ) = 1
                        write (buff,'(a,3i7,a,i7,a)')  'wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at ', &
                        i,j,k,' for WIRE-segment ',segmento%origIndex,' '//dir(whatfield)
                        call WarnErrReport(buff)
                     endif
                  endif
                  if (.not.sgg%med(sggmiEy(i+1,j  ,k      ))%Is%ThinWire)  then
                     sggmiEy(i+1,j  ,k      ) = 1
                     write (buff,'(a,3i7,a,i7,a)')  'wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at ', &
                     i,j,k,' for WIRE-segment ',segmento%origIndex,' '//dir(whatfield)
                     call WarnErrReport(buff)
                  endif
                  if (.not.sgg%med(sggmiEy(i+1,j-1,k      ))%Is%ThinWire)  then
                     sggmiEy(i+1,j-1,k      ) = 1
                     write (buff,'(a,3i7,a,i7,a)')  'wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at ', &
                     i,j,k,' for WIRE-segment ',segmento%origIndex,' '//dir(whatfield)
                     call WarnErrReport(buff)
                  endif
                  if         ((k <=  sgg%alloc(iEz)%ZE).and.(k>=sgg%alloc(iEz)%ZI)) then
                     if (.not.sgg%med(sggmiEz(i+1,j  ,k      ))%Is%ThinWire)  then
                        sggmiEz(i+1,j  ,k      ) = 1
                        write (buff,'(a,3i7,a,i7,a)')  'wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at ', &
                        i,j,k,' for WIRE-segment ',segmento%origIndex,' '//dir(whatfield)
                        call WarnErrReport(buff)
                     endif
                  endif
                  if ((k -1 <=  sgg%alloc(iEz)%ZE).and.(k - 1 >=sgg%alloc(iEz)%ZI)) then
                     if (.not.sgg%med(sggmiEz(i+1,j  ,k - 1  ))%Is%ThinWire)  then
                        sggmiEz(i+1,j  ,k - 1  ) = 1
                        write (buff,'(a,3i7,a,i7,a)')  'wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at ', &
                        i,j,k,' for WIRE-segment ',segmento%origIndex,' '//dir(whatfield)
                        call WarnErrReport(buff)
                     endif
                  endif
                case (iEy)
                  sggmiHy(i   ,j   ,k      ) = 1
                  sggmiHy(i   ,j   ,k-1    ) = 1
                  sggmiHy(i -1,j   ,k      ) = 1
                  sggmiHy(i -1,j   ,k-1    ) = 1
                  sggmiHy(i   ,j+1 ,k      ) = 1
                  sggmiHy(i   ,j+1 ,k-1    ) = 1
                  sggmiHy(i -1,j+1 ,k      ) = 1
                  sggmiHy(i -1,j+1 ,k-1    ) = 1
                  if (.not.sgg%med(sggmiEy(i  ,j  ,k      ))%Is%ThinWire)  then
                     sggmiEy(i  ,j  ,k      ) = 1
                     write (buff,'(a,3i7,a,i7,a)')  'wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at ', &
                     i,j,k,' for WIRE-segment ',segmento%origIndex,' '//dir(whatfield)
                     call WarnErrReport(buff)
                  endif
                  if     ((k <=  sgg%alloc(iEz)%ZE).and.(k >= sgg%alloc(iEz)%ZI)) then
                     if (.not.sgg%med(sggmiEz(i  ,j  ,k      ))%Is%ThinWire)  then
                        sggmiEz(i  ,j  ,k      ) = 1
                        write (buff,'(a,3i7,a,i7,a)')  'wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at ', &
                        i,j,k,' for WIRE-segment ',segmento%origIndex,' '//dir(whatfield)
                        call WarnErrReport(buff)
                     endif
                  endif
                  if ((k-1 <=  sgg%alloc(iEz)%ZE).and.(k -1 >=sgg%alloc(iEz)%ZI)) then
                     if (.not.sgg%med(sggmiEz(i  ,j  ,k - 1  ))%Is%ThinWire)  then
                        sggmiEz(i  ,j  ,k - 1  ) = 1
                        write (buff,'(a,3i7,a,i7,a)')  'wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at ', &
                        i,j,k,' for WIRE-segment ',segmento%origIndex,' '//dir(whatfield)
                        call WarnErrReport(buff)
                     endif
                  endif
                  if (.not.sgg%med(sggmiEx(i  ,j  ,k      ))%Is%ThinWire)  then
                     sggmiEx(i  ,j  ,k      ) = 1
                     write (buff,'(a,3i7,a,i7,a)')  'wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at ', &
                     i,j,k,' for WIRE-segment ',segmento%origIndex,' '//dir(whatfield)
                     call WarnErrReport(buff)
                  endif
                  if (.not.sgg%med(sggmiEx(i-1,j  ,k      ))%Is%ThinWire)  then
                     sggmiEx(i-1,j  ,k      ) = 1
                     write (buff,'(a,3i7,a,i7,a)')  'wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at ', &
                     i,j,k,' for WIRE-segment ',segmento%origIndex,' '//dir(whatfield)
                     call WarnErrReport(buff)
                  endif
                  if ((k     <=  sgg%alloc(iEz)%ZE).and.(k >= sgg%alloc(iEz)%ZI)) then
                     if (.not.sgg%med(sggmiEz(i  ,j+1,k      ))%Is%ThinWire)  then
                        sggmiEz(i  ,j+1,k      ) = 1
                        write (buff,'(a,3i7,a,i7,a)')  'wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at ', &
                        i,j,k,' for WIRE-segment ',segmento%origIndex,' '//dir(whatfield)
                        call WarnErrReport(buff)
                     endif
                  endif
                  if ((k-1 <=  sgg%alloc(iEz)%ZE).and.(k -1 >=sgg%alloc(iEz)%ZI)) then
                     if (.not.sgg%med(sggmiEz(i  ,j+1,k - 1  ))%Is%ThinWire)  then
                        sggmiEz(i  ,j+1,k - 1  ) = 1
                        write (buff,'(a,3i7,a,i7,a)')  'wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at ', &
                        i,j,k,' for WIRE-segment ',segmento%origIndex,' '//dir(whatfield)
                        call WarnErrReport(buff)
                     endif
                  endif
                  if (.not.sgg%med(sggmiEx(i  ,j+1,k      ))%Is%ThinWire)  then
                     sggmiEx(i  ,j+1,k      ) = 1
                     write (buff,'(a,3i7,a,i7,a)')  'wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at ', &
                     i,j,k,' for WIRE-segment ',segmento%origIndex,' '//dir(whatfield)
                     call WarnErrReport(buff)
                  endif
                  if (.not.sgg%med(sggmiEx(i-1,j+1,k      ))%Is%ThinWire)  then
                     sggmiEx(i-1,j+1,k      ) = 1
                     write (buff,'(a,3i7,a,i7,a)')  'wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at ', &
                     i,j,k,' for WIRE-segment ',segmento%origIndex,' '//dir(whatfield)
                     call WarnErrReport(buff)
                  endif
                case (iEz)
                  sggmiHz(i    ,j   ,k    ) = 1
                  sggmiHz(i -1 ,j   ,k    ) = 1
                  sggmiHz(i    ,j-1 ,k    ) = 1
                  sggmiHz(i -1 ,j-1 ,k    ) = 1
                  sggmiHz(i    ,j   ,k+1 ) = 1
                  sggmiHz(i -1 ,j   ,k+1 ) = 1
                  sggmiHz(i    ,j-1 ,k+1 ) = 1
                  sggmiHz(i -1 ,j-1 ,k+1 ) = 1
                  if (.not.sgg%med(sggmiEz(i  ,j  ,k      ))%Is%ThinWire)  then
                     sggmiEz(i  ,j  ,k      ) = 1
                     write (buff,'(a,3i7,a,i7,a)')  'wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at ', &
                     i,j,k,' for WIRE-segment ',segmento%origIndex,' '//dir(whatfield)
                     call WarnErrReport(buff)
                  endif
                  if ((k   <=  sgg%alloc(iEx)%ZE).and.(k  >= sgg%alloc(iEx)%ZI)) then
                     if (.not.sgg%med(sggmiEx(i  ,j  ,k      ))%Is%ThinWire)  then
                        sggmiEx(i  ,j  ,k      ) = 1
                        write (buff,'(a,3i7,a,i7,a)')  'wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at ', &
                        i,j,k,' for WIRE-segment ',segmento%origIndex,' '//dir(whatfield)
                        call WarnErrReport(buff)
                     endif
                     if (.not.sgg%med(sggmiEx(i-1,j  ,k      ))%Is%ThinWire)  then
                        sggmiEx(i-1,j  ,k      ) = 1
                        write (buff,'(a,3i7,a,i7,a)')  'wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at ', &
                        i,j,k,' for WIRE-segment ',segmento%origIndex,' '//dir(whatfield)
                        call WarnErrReport(buff)
                     endif
                  endif
                  if ((k   <=  sgg%alloc(iEy)%ZE).and.(k  >= sgg%alloc(iEy)%ZI)) then
                     if (.not.sgg%med(sggmiEy(i  ,j  ,k      ))%Is%ThinWire)  then
                        sggmiEy(i  ,j  ,k      ) = 1
                        write (buff,'(a,3i7,a,i7,a)')  'wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at ', &
                        i,j,k,' for WIRE-segment ',segmento%origIndex,' '//dir(whatfield)
                        call WarnErrReport(buff)
                     endif
                     if (.not.sgg%med(sggmiEy(i  ,j-1,k      ))%Is%ThinWire)  then
                        sggmiEy(i  ,j-1,k      ) = 1
                        write (buff,'(a,3i7,a,i7,a)')  'wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at ', &
                        i,j,k,' for WIRE-segment ',segmento%origIndex,' '//dir(whatfield)
                        call WarnErrReport(buff)
                     endif
                  endif
                  if ((k + 1  <=  sgg%alloc(iEx)%ZE).and.(k + 1 >= sgg%alloc(iEx)%ZI)) then
                     if (.not.sgg%med(sggmiEx(i  ,j  ,k + 1  ))%Is%ThinWire)  then
                        sggmiEx(i  ,j  ,k + 1  ) = 1
                        write (buff,'(a,3i7,a,i7,a)')  'wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at ', &
                        i,j,k,' for WIRE-segment ',segmento%origIndex,' '//dir(whatfield)
                        call WarnErrReport(buff)
                     endif
                     if (.not.sgg%med(sggmiEx(i-1,j  ,k + 1  ))%Is%ThinWire)  then
                        sggmiEx(i-1,j  ,k + 1  ) = 1
                        write (buff,'(a,3i7,a,i7,a)')  'wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at ', &
                        i,j,k,' for WIRE-segment ',segmento%origIndex,' '//dir(whatfield)
                        call WarnErrReport(buff)
                     endif
                  endif
                  if ((k + 1  <=  sgg%alloc(iEy)%ZE).and.(k + 1 >= sgg%alloc(iEy)%ZI)) then
                     if (.not.sgg%med(sggmiEy(i  ,j  ,k + 1  ))%Is%ThinWire)  then
                        sggmiEy(i  ,j  ,k + 1  ) = 1
                        write (buff,'(a,3i7,a,i7,a)')  'wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at ', &
                        i,j,k,' for WIRE-segment ',segmento%origIndex,' '//dir(whatfield)
                        call WarnErrReport(buff)
                     endif
                     if (.not.sgg%med(sggmiEy(i  ,j-1,k + 1  ))%Is%ThinWire)  then
                        sggmiEy(i  ,j-1,k + 1  ) = 1
                        write (buff,'(a,3i7,a,i7,a)')  'wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at ', &
                        i,j,k,' for WIRE-segment ',segmento%origIndex,' '//dir(whatfield)
                        call WarnErrReport(buff)
                     endif
                  endif
               end select
            endif
         endif
      end do

!!!calculo y gestion autoinduccion
      do i1=1,HWires%NumCurrentSegments
         jmed=HWires%CurrentSegment(i1)%indexmed
         desp=HWires%CurrentSegment(i1)%delta
         despT1=HWires%CurrentSegment(i1)%deltaTransv1
         despT2=HWires%CurrentSegment(i1)%deltaTransv2
         ! Esto no debe ser preciso. la hipotesis es que vuelve promediando a una celda
         if (wirethickness>1) then        !solo tiene impacto en la %cte2 de acople y en la fuente de carga para que sea de voltaje
           despT1=wirethickness*despT1
           despT2=wirethickness*despT2
         endif
         r0=HWires%CurrentSegment(i1)%TipoWire%Radius   ! CON SU RADIO
         
          if (r0 < 1e-30) then 
              write (buff,*)  'wir0_ERROR: ire radius cannot be null'
              if ((k1 >= ZI).and.(k1 <= ZE))  call WarnErrReport(buff,.true.)
          endif   
         if ((r0 < 1e-9*desp)) then
            write(BUFF,'(a,e12.2e3)') 'wir0_WARNING: WIRE radius too small ',r0
            call WarnErrReport(buff)
         endif
         if ((r0 > 0.5_RKIND_wires *despT1).or.(r0 > 0.5_RKIND_wires *despT2)) then
            write(BUFF,'(a,e12.2e3)') 'wir0_WARNING: WIRE radius greater '// &
            'than half a space-step. Reduced to this limit',min(0.5_RKIND_wires *despT1,0.5_RKIND_wires *despT2)
            if ((HWires%CurrentSegment(i1)%k >= ZI).and.(HWires%CurrentSegment(i1)%k <= ZE)) call WarnErrReport(buff)
         endif
         select case (trim(adjustl(inductance_model)))
          case ('berenger')
            !---------------------------------------------------------------------------------
            !!Slanted=Berenger
            !The second best one for the edelvik PEC box (Boutayeb's PPT is the best one for this case)
            !menos ruidoso que el de Boutayeb en uniones de hilos paralelos
            HWires%CurrentSegment(i1)%Lind = &
            (1.0_RKIND_wires / (4.0_RKIND_wires * pi*InvMu(jmed)))*(log((despT1**2.0_RKIND_wires +despT2**2.0_RKIND_wires )/(4*r0**2.0_RKIND_wires ))  + &
            despT1/despT2*atan(despT2/despT1)     + &
            despT2/despT1*atan(despT1/despT2)     + &
            pi*r0**2.0_RKIND_wires /(despT2*despT1)-3.0_RKIND_wires)
            !!---------------------------------------------------------------------------------
            !Slanted corrected for wires of radius > 0.3_RKIND_wires  Delta
            ! Untested
            !just divides by a correction factor equal to that used by Boutayeb in his correction
            !proposed by Grando
            if ((r0 >0.3_RKIND_wires  *despT1).or.(r0 >0.3_RKIND_wires *despT2)) then
               HWires%CurrentSegment(i1)%Lind = HWires%CurrentSegment(i1)%Lind &
               /(1.0_RKIND_wires-pi*r0**2.0_RKIND_wires /(despT1*despT2))
            endif
          case ('ledfelt')
            !---------------------------------------------------------------------------------
            !Ledfelt thesis square symbols
            !The best for the antenna dipole, but it is not the best for edelvik's cavity
            !Dec'11 me ha dado muy bien al revisar el TC2_OF-SWAB-A2-B
            !warning: there is a mistake in formula 8.15 parenthesis
            HWires%CurrentSegment(i1)%Lind = &
            (1.0_RKIND_wires / (4.0_RKIND_wires * pi*InvMu(jmed)))*(log((despT1**2.0_RKIND_wires +despT2**2.0_RKIND_wires )/(  r0**2.0_RKIND_wires ))+ &
            despT1/despT2*atan(despT2/despT1)+ &
            despT2/despT1*atan(despT1/despT2)+ pi*r0**2.0_RKIND_wires /(16.0_RKIND_wires *despT2*despT1)-3.0_RKIND_wires)
            !Slanted makes this correction for radius>0.3_RKIND_wires  Delta.  a>0.3_RKIND_wires  Delta
            !I use it also in Ledleft !2012
            !never tested
            if ((r0 > 0.3_RKIND_wires  *despT1).or.(r0 > 0.3_RKIND_wires *despT2)) then
               HWires%CurrentSegment(i1)%Lind =  HWires%CurrentSegment(i1)%Lind &
               /(1.0_RKIND_wires-pi*r0**2.0_RKIND_wires /(despT1*despT2))
            endif
          case ('boutayeb')
            !!---------------------------------------------------------------------------------
            !Boutayeb's PPT. For radius>0.3_RKIND_wires  Delta coincides with that of Berenger
            !adds a correction for radiuos <0.3_RKIND_wires  Delta
            !and is equal to Guiffault for radius > 0.5_RKIND_wires  Delta
            !The BEST ONE for edelvik_box and for HIRF NTC2.0_RKIND_wires
            HWires%CurrentSegment(i1)%Lind = &
            (1.0_RKIND_wires / (4.0_RKIND_wires * pi*InvMu(jmed)))*(log((despT1**2.0_RKIND_wires +despT2**2.0_RKIND_wires )/(4.0_RKIND_wires *r0**2.0_RKIND_wires ))+     &
            despT1/despT2*atan(despT2/despT1)+        &
            despT2/despT1*atan(despT1/despT2)+ pi*r0**2.0_RKIND_wires /(despT2*despT1)-3.0_RKIND_wires)

            if ((r0 < 0.3_RKIND_wires  *despT1).or.(r0 < 0.3_RKIND_wires *despT2)) then
               HWires%CurrentSegment(i1)%Lind = HWires%CurrentSegment(i1)%Lind - 0.57_RKIND_wires/(4.0_RKIND_wires * pi*InvMu(jmed))
            endif
            !---------------------------------------------------------------------------------
            !Boutayeb's PPT for radius>0.3_RKIND_wires  Delta.  a>0.3_RKIND_wires  Delta.
            !(Slanted corrects for 0.3_RKIND_wires  Delta while boutayeb does it for 0.5_RKIND_wires  Delta, I take Slanted's)
            !Untested
            !just divides by a correction factor (warning becomes negative for r0/delta >0.56)
            if ((r0 > 0.3_RKIND_wires  *despT1).or.(r0 > 0.3_RKIND_wires *despT2)) then
               HWires%CurrentSegment(i1)%Lind =  HWires%CurrentSegment(i1)%Lind &
               /(1.0_RKIND_wires-pi*r0**2.0_RKIND_wires /(despT1*despT2))
            endif
            !!---------------------------------------------------------------------------------
          case default
            BUFF='wir0_ERROR: Incorrect inductance model'
            call WarnErrReport(buff,.TRUE.)
         end select
         if (HWires%CurrentSegment(i1)%Lind < 0.0_RKIND_wires) then
            BUFF='wir0_ERROR: Wrong self-inductance. '
            call WarnErrReport(buff,.TRUE.)
         endif
         !!!!for junctions after holland
         !!!    HWires%CurrentSegment(i1)%logRoverR0=log(((desp*despT1*despT2)**(1.0_RKIND_wires/3.0_RKIND_wires))/r0)  !HWires%CurrentSegment(i1)%Lind ! !holland pp90 (17)(18)
         !!!!22/07/15 el cambio de radio es equivalente a un cambio de delta.... uso algo mejor
         !!!  HWires%CurrentSegment(i1)%logRoverR0=HWires%CurrentSegment(i1)%Lind * InvMu(jmed) * InvEPS(jmed)
         !!!!

      end do
      !
      allocate (LindProb(1:HWires%NumCurrentSegments ))
      do i1=1,HWires%NumCurrentSegments
         LindProb(i1)=.true.
         HWires%CurrentSegment(i1)%inv_Lind_acum    = 1.0_RKIND_wires /(  HWires%CurrentSegment(i1)%Lind )
         HWires%CurrentSegment(i1)%HEUR_safety =   (sgg%dt**2.0_RKIND/ (eps0* HWires%CurrentSegment(i1)%deltaTransv1*HWires%CurrentSegment(i1)%deltaTransv2))
      end do
      do i1=1,HWires%NumCurrentSegments
         if (LindProb(i1)) then
            org=>HWires%CurrentSegment(i1)
            org%NumParallel=1
            org%Lind_acum     = org%Lind
            do j1=i1+1,HWires%NumCurrentSegments
               fin=>HWires%CurrentSegment(j1)
               if ( (org%i == fin%i).and.(org%j == fin%j).and.(org%k == fin%k).and.(org%tipofield == fin%tipofield)) then
                  org%NumParallel=org%NumParallel + 1
                  if (stableradholland) org%Lind_acum     = org%Lind_acum     + fin%Lind
               endif
            end do
            do j1=i1+1,HWires%NumCurrentSegments
               fin=>HWires%CurrentSegment(j1)
               if ( (org%i == fin%i).and.(org%j == fin%j).and.(org%k == fin%k).and.(org%tipofield == fin%tipofield)) then
                  fin%NumParallel=org%NumParallel
                  fin%Lind_acum     = org%Lind_acum
                  LindProb(j1)=.false.
               endif
            end do
         endif
      end do
      deallocate (LindProb)


      !!!!solo para pruebas
      !!!open (629,file='param.txt',form='formatted')
      !!!    read (629,*) deltadummy
      !!!close (629)
      !!!    print *,'---------deltadummy=',deltadummy
      !!!do i1=1,HWires%NumCurrentSegments
      !!!    HWires%CurrentSegment(i1)%inv_Lind_acum=deltadummy*mu0
      !!!end do
      !!!!!!!!
      !!!!!!!!!!!!!!!!!\E7\E7\E7\E7\E7\E7\E7\E7\E7\E7
      !!!!!!!! hago 120715 mi criterior
      dtcritico=sgg%dt
      do i1=1,HWires%NumCurrentSegments
         dummy => HWires%CurrentSegment(i1)
         jmed=dummy%indexmed
         desp=  dummy%delta
         despT1=dummy%deltaTransv1
         despT2=dummy%deltaTransv2
         r0=    dummy%TipoWire%Radius
         if (r0 < 1e-30) then 
              write (buff,*)  'wir0_ERROR: wire radius cannot be null'
              if ((k1 >= ZI).and.(k1 <= ZE))  call WarnErrReport(buff,.true.)
         endif   
         !!!!!!!!!!!!correccion bruta 13/07/15
         deltadummy=dummy%inv_Lind_acum * dummy%HEUR_safety/0.9_RKIND_wires  !EL 0.9 ES POR MI TRANQUILIDAD
         if (deltadummy > 1.0_RKIND_wires) then
            if (stableradholland) then
               if (trim(adjustl(inductance_model))=='boutayeb') then
                  b= -4.0 * pi*InvMu(jmed) * (dummy%Lind*deltadummy) + log(despT1**2.0_RKIND+despT2**2.0_RKIND)+ despT1/despT2*atan(despT2/despT1)+ despT2/despT1*atan(despT1/despT2)-3.0_RKIND_wires
                  a= pi/(despT2*despT1)
                  if ((r0 < 0.3_RKIND_wires  *despT1).or.(r0 < 0.3_RKIND_wires *despT2)) then
                     B= B-0.57
                  endif
                  newr0=Sqrt(-Lambert(-A*Exp(b)/4.0_RKIND_wires) /A )
                  !!!doublechecking
                  b= -4.0 * pi*InvMu(jmed) * (dummy%Lind)            + log(despT1**2.0_RKIND+despT2**2.0_RKIND)+ despT1/despT2*atan(despT2/despT1)+ despT2/despT1*atan(despT1/despT2)-3.0_RKIND_wires
                  a= pi/(despT2*despT1)
                  if ((r0 < 0.3_RKIND_wires  *despT1).or.(r0 < 0.3_RKIND_wires *despT2)) then
                     B=B-0.57
                  endif
                  OLDR0=Sqrt(-Lambert(-A*Exp(b)/4.0_RKIND_wires) /A )
                  !!!!!!!!!!
                  write (buff,'(a,e10.2e3,a,2e10.2e3,a,i9,a,3i9,a,e10.2e3,a,e10.2e3)')  'wir0_WARNING: AUTOMATIC CORRECTION OF L/mu0=', dummy%Lind/mu0, ' for r0=',r0,oldr0,&
                  ' ',dummy%NumParallel, ' wires at ',dummy%i,dummy%j,dummy%k,' to L/mu0=',dummy%Lind*deltadummy/mu0, ' for newr0=',newr0
               else
                  write (buff,'(a,e10.2e3,a,i9,a,3i9,a,e10.2e3)')  'wir0_WARNING: AUTOMATIC CORRECTION OF L/mu0=', dummy%Lind/mu0, &
                  ' ',dummy%NumParallel, ' wires at ',dummy%i,dummy%j,dummy%k,' to L/mu0=',dummy%Lind*deltadummy/mu0
               endif
               if ((dummy%k > ZI).and.(dummy%k <= ZE)) call WarnErrReport(buff)
               dummy%Lind = dummy%Lind* deltadummy !bajo repartiendo proporcialmente
            else
               write (buff,'(a,e10.2e3,a,i9,a,3i9,a,e10.2e3)')  'wir0_SEVEREWARNING: L/mu0=', dummy%Lind/mu0, &
               ' in ',dummy%NumParallel, ' wires at ',dummy%i,dummy%j,dummy%k,' smaller (posibly unstable) than L/mu0=',dummy%Lind*deltadummy/mu0
               if ((dummy%k > ZI).and.(dummy%k <= ZE)) call WarnErrReport(buff)
               dtcritico=min(sgg%dt/sqrt(deltadummy),dtcritico)
            endif
         endif
      end do
      !!if (dtcritico<sgg%dt) then
      !!            write(buff,'(a,e9.2e2,a,e9.2e2)') &
      !!            &    'wir0_ERROR: UNSTABLE sgg%dt, decrease wire radius, number of Parallel WIREs, or make sgg%dt < ',dtcritico
      !!            call WarnErrReport(buff,.true.)
      !!endif

      !!!!!!!!!!!fin !mi criterio 13/07/15


      !Grounding R_RightEnd and R_LeftEnd resistances info
      do i1=1,HWires%NumCurrentSegments
         segmento=>HWires%CurrentSegment(i1)
         !
         if (segmento%TipoWire%HasAbsorbing_LeftEnd) then
            if (segmento%Is_LeftEnd) then
               segmento%HasAbsorbing_LeftEnd=.true.
               write (buff,'(a,5i7)')  'wir1_INFO: Absorbing conditions in terminal LeftEnd segment ', &
               segmento%origIndex,segmento%i,segmento%j,segmento%k,segmento%tipofield
               if ((segmento%k >= ZI).and.(segmento%k <= ZE).and.verbose) call WarnErrReport(buff)
            endif
         endif
         if (segmento%TipoWire%HasAbsorbing_RightEnd) then
            if (segmento%Is_RightEnd) then
               segmento%HasAbsorbing_RightEnd=.true.
               write (buff,'(a,5i7)')  'wir1_WARNING: Absorbing conditions  in terminal RightEnd segment ', &
               segmento%origIndex,segmento%i,segmento%j,segmento%k,segmento%tipofield
               if ((segmento%k >= ZI).and.(segmento%k <= ZE)) call WarnErrReport(buff)
            endif
         endif
         !
         if (segmento%TipoWire%HasParallel_LeftEnd) then
            if (segmento%Is_LeftEnd) then
               segmento%HasParallel_LeftEnd=.true.
               write (buff,'(a,5i7)')  'wir1_WARNING: Parallel RLC in terminal LeftEnd segment ', &
               segmento%origIndex,segmento%i,segmento%j,segmento%k,segmento%tipofield
               if ((segmento%k >= ZI).and.(segmento%k <= ZE)) call WarnErrReport(buff)
            endif
         endif
         if (segmento%TipoWire%HasParallel_RightEnd) then
            if (segmento%Is_RightEnd) then
               segmento%HasParallel_RightEnd=.true.
               write (buff,'(a,5i7)')  'wir1_WARNING: Parallel RLC in terminal RightEnd segment ', &
               segmento%origIndex,segmento%i,segmento%j,segmento%k,segmento%tipofield
               if ((segmento%k >= ZI).and.(segmento%k <= ZE)) call WarnErrReport(buff)
            endif
         endif

         if (segmento%TipoWire%HasSeries_LeftEnd) then
            if (segmento%Is_LeftEnd) then
               segmento%HasSeries_LeftEnd=.true.
               write (buff,'(a,5i7)')  'wir1_WARNING: Series RLC in terminal LeftEnd segment ', &
               segmento%origIndex,segmento%i,segmento%j,segmento%k,segmento%tipofield
               if ((segmento%k >= ZI).and.(segmento%k <= ZE)) call WarnErrReport(buff)
            endif
         endif
         if (segmento%TipoWire%HasSeries_RightEnd) then
            if (segmento%Is_RightEnd) then
               segmento%HasSeries_RightEnd=.true.
               write (buff,'(a,5i7)')  'wir1_WARNING: Series RLC in terminal RightEnd segment ', &
               segmento%origIndex,segmento%i,segmento%j,segmento%k,segmento%tipofield
               if ((segmento%k >= ZI).and.(segmento%k <= ZE)) call WarnErrReport(buff)
            endif
         endif


      end do
      !
      !Create the final update constants for the advance of the currents
      !It takes into account the extra inductance and resistance per unit length specified in ORIGINAL
      !It also takes into account the Series/Parallel Grounding Inductance at and the end segments TR and TL !untested
      !Junctions do no affect to these constants (later taken into account by means of the fractionplus and
      !fractionminus constants)
      do i1=1,HWires%NumCurrentSegments
         !constantes de actualizacion
         dummy=> HWires%CurrentSegment(i1)
         resist=0.0_RKIND_wires

         !!!for lossy groundings
         i = dummy%i
         j = dummy%j
         k = dummy%k
         whatfield=dummy%tipofield    
         !
         rlossy=0.0_RKIND_wires
         sigt=0.0_RKIND_wires
         sigtPlus=0.0_RKIND_wires
         sigtMinu=0.0_RKIND_wires
         IsLossy=.false. ; IsPEC = .false. ; 
         IsLossyPlus=.false. ; IsPECPlus = .false. ;  
         IsLossyMinu=.false. ; IsPECMinu = .false. ;       
!!!! 
         select case (whatfield)
             case (iEx)
                esPML= sgg%med(sggmiEx(i,j,k))%is%PML
             case (iEy)
                esPML= sgg%med(sggmiEy(i,j,k))%is%PML
             case (iEz)
                esPML= sgg%med(sggmiEz(i,j,k))%is%PML
         end select
         if (esPML) then 
            continue
            !!\C7  dummy%IsShielded=.true.
         else
             if ((k <= sgg%alloc(iEZ)%ZE).and.(k >= sgg%alloc(iEZ)%ZI)) then
                kmenos1= k-1
                kmas1  = k+1
                !
                !esta informacion solo se utiliza si realmente luego hay un nodo terminal y se suma la resistencia. En cualquier otro caso no sirver para nada
                !de todos modos hay un bug en la deteccion. sgg 110815
                ! se usal la informacion nodal (lo que sigue algun d\EDa se mover\E1 a la rutina de generacion nodal que se creo en preprocess y se podra dejar solo lo que sigue 110815
                !!!bug 270224 gg Cuando se sobrepasa el MPI  kmenos1=k o kmas1=k da un buggy error. Pero no se toma decision alguna. se puede ignorar con -ignoreerrors
                select case (whatfield)
                 case (iEx)
                   med(0)  = sggMiEx(i + 1  , j     , k       )
                   med(1)  = sggMiEy(i + 1  , j     , k       )
                   med(2)  = sggMiEy(i + 1  , j - 1 , k       )
                   med(3)  = sggMiEz(i + 1  , j     , k       )
                   if (kmenos1 <  sgg%alloc(iEz)%ZI) kmenos1=k; med(4)  = sggMiEz(i + 1  , j     , kmenos1 )
                   med(5)  = sggMiNo(i + 1  , j     , k       )
                   !
                   med(6)  = sggMiEx(i  - 1 , j     , k       )
                   med(7)  = sggMiEy(i      , j     , k       )
                   med(8)  = sggMiEy(i      , j - 1 , k       )
                   med(9)  = sggMiEz(i      , j     , k       )
                   if (kmenos1 <  sgg%alloc(iEz)%ZI) kmenos1=k; med(10) = sggMiEz(i      , j     , kmenos1 )
                   med(11) = sggMiNo(i      , j     , k       )
                 case (iEy)
                   med(0)  = sggMiEy(i      , j + 1 , k       )
                   med(1)  = sggMiEz(i      , j + 1 , k       )
                   if (kmenos1 <  sgg%alloc(iEz)%ZI) kmenos1=k; med(2)  = sggMiEz(i      , j + 1 , kmenos1 )
                   med(3)  = sggMiEx(i      , j + 1 , k       )
                   med(4)  = sggMiEx(i - 1  , j + 1 , k       )
                   med(5)  = sggMiNo(i      , j + 1 , k       )
                   !
                   med(6)  = sggMiEy(i      , j  - 1 , k      )
                   med(7)  = sggMiEz(i      , j     , k       )
                   if (kmenos1 <  sgg%alloc(iEz)%ZI) kmenos1=k; med(8)  = sggMiEz(i      , j     , kmenos1 )
                   med(9)  = sggMiEx(i      , j     , k       )
                   med(10) = sggMiEx(i - 1  , j     , k       )
                   med(11) = sggMiNo(i      , j     , k       )
                 case (iEz)
!!!ojooooo 270224 esta logica esta mal porque machaco las variables kmas1 y kmenos1.... corregir.... no tiene impacto pq el bucle es informativo y no se toman decisiones. en todo caso solo afectaria a MPI!
                     !!!pero esta no es la razon 27024 por la que sucede el error gg   (    8/   40) wir1_BUGGYERROR:  Lossy, pec,  2.           329         180         160   15000.0000000000      F F
                   if (kmas1   >  sgg%alloc(iEz)%ZE) kmas1=k  ; med(0)  = sggMiEz(i     , j     , kmas1    )
                   if (kmas1   >  sgg%alloc(iEx)%ZE) kmas1=k  ; med(1)  = sggMiEx(i     , j     , kmas1    )
                   if (kmas1   >  sgg%alloc(iEx)%ZE) kmas1=k  ; med(2)  = sggMiEx(i - 1 , j     , kmas1    )
                   if (kmas1   >  sgg%alloc(iEy)%ZE) kmas1=k  ; med(3)  = sggMiEy(i     , j     , kmas1    )
                   if (kmas1   >  sgg%alloc(iEy)%ZE) kmas1=k  ; med(4)  = sggMiEy(i     , j - 1 , kmas1    )
                   if (kmas1   >  sgg%alloc(iHz)%ZE) kmas1=k  ; med(5)  = sggMiNo(i     , j     , kmas1    )
                   !
                   if (kmenos1 <  sgg%alloc(iEz)%ZI) kmenos1=k; med(6)  = sggMiEz(i     , j     , kmenos1  )
                   med(7)  = sggMiEx(i     , j     , k        )
                   med(8)  = sggMiEx(i - 1 , j     , k        )
                   med(9)  = sggMiEy(i     , j     , k        )
                   med(10) = sggMiEy(i     , j - 1 , k        )
                   med(11) = sggMiNo(i     , j     , k        )
                end select
                !hay que tratar cada extremo por separado
                do nm=0,5
                   IsLossyPlus             = IsLossyPlus             .or. sgg%Med(med(nm))%Is%Lossy
                   IsPECPlus                   = IsPECPlus           .or. sgg%med(med(nm))%is%PEC .or. (med(nm)==0)
                   if (.not.IsPECPlus) then
                      sigtPlus                    = max(sigtPlus,                    sgg%Med(med(nm))%sigma)
                   endif
                end do
                do nm=6,11
                   IsLossyMinu             = IsLossyMinu             .or. sgg%Med(med(nm))%Is%Lossy
                   IsPECMinu               = IsPECMinu               .or. sgg%med(med(nm))%is%PEC .or. (med(nm)==0)
                   if (.not.IsPECMInu) then
                      sigtMinu                    = max(sigtMinu,                    sgg%Med(med(nm))%sigma)
                   endif
                end do
!!!telaranias quitadas 060215
                if (IsPECPlus) then !domina el pec en un nodo con multiples medios
                   IsLossyPlus=.false.
                endif
                if (IsPECMinu) then !domina el pec en un nodo con multiples medios
                   IsLossyMinu=.false.
                endif
                sigt                  = max(sigtplus,sigtMinu) !con que uno de los dos extremos sea Lossy hay que corregir la resistencia de contacto
                IsLossy               = IsLossyPlus               .or. IsLossyMinu !con que uno de los dos extremos sea Lossy hay que corregir la resistencia de contacto
                IsPEC                 = IsPECPlus                 .and.  IsPECMinu  !!!importante para ser PEC tienen que serlo los dos (es un caso trivial de un segemento unido a pec por los dos sitios. pero se da en el siva)
                !!!checking de coherencia
                if (isPEC.and.IsLossy) then !es decir: los dos extremos PEC y alguno lossy es que hay un error
                      write (buff,*)  'wir1_BUGGYERROR:  Lossy, pec, 1.  ', i,j,k,sigt,ispec,isLossy
                      call WarnErrReport(buff,.true.)
                endif
                !!!
                if ((.not.ispec).and.(.not.isLossy)) then  !algun extremo no pec y ambos extremos no lossy debe haber error si la conductividad no es nula
                   if (abs(sigt) > 1.0e-19_RKIND_wires) then
                       !270224 creo que el bug gg 270224 es simplemente que conviven pec y lossy en plus o minus siendo el otro minus/plus vacio. Se voidea
                       !270224 el islosyplus/minus a false por ser pec. el ispec es false porque uno es vacio y el islossy tambien es vacio
                       !270224 pero se ha almacenado alguna sigt no nula. El arreglo es simplemente convertir esto en un warning de casuistica de conexionado en vez de un buggyerror
                      !!comentado a 010324 por los motivos anteriores    write (buff,*)  'wir1_BUGGYERROR:  Lossy, pec,  2.  ', i,j,k,sigt,ispec,isLossy
                      !!comentado a 010324 por los motivos anteriores    call WarnErrReport(buff,.true.)
                   endif
                   write (buff,*)  'wir1_WARNING:  Lossy, pec,  2.: A wire is connected at some ending both to a lossy and to a PEC edge. Assuming it PEC  ', i,j,k,sigt,ispec,isLossy,IsPECPlus,IsPECMinu
                   call WarnErrReport(buff)
                endif
                if (isLossy) then !alguno extremo lossy con conductividad desconocida (multiports de ss)
                   if (abs(sigt) < 1.0e-19_RKIND_wires) then
                      sigt=1e4 ! asignale una resistencia de contacto por defecto si es nula !tipico de los composites de ss !habra algun dia que afinar esto
                      write (buff,*)  'wir1_WARNING:  A Lossy segment with unknown conductivity. Assuming a STANDARD value of 1e4 S/m ', i,j,k,sigt,ispec,isLossy,IsLossyPlus,IsLossyMinu
                      call WarnErrReport(buff)
                   endif
                endif

                !!!!!!!!!!!!!!!!!!!!hasta aqui casuistica. Aniade ahora la resistencia si procede con la formula de tercero de fisicas
                if (isLossy) then
                   !rlossy=rlossy + 1.0_RKIND_wires/(2.0_RKIND_wires * pi*dummy%TipoWire%radius*sigt)/dummy%delta   !p.u.l.     !pruebas distintas 0223
                    !rlossy=rlossy + 1.0_RKIND_wires/(2.0_RKIND_wires * pi*dummy%DELTA          *sigt)/dummy%delta   !p.u.l.     !pruebas distintas 0223
                    rlossy=rlossy ! + 1.0_RKIND_wires/(2.0_RKIND_wires * pi * (factordelta*dummy%DELTA +factorradius*dummy%TipoWire%radius) *sigt)/dummy%delta   !p.u.l. !comentado a 290323 porque no me fio de esto hasta que no se valide bien agb
                endif
             endif
         endif !del esPML
         !
         resist=dummy%TipoWire%R
         givenautoin = HWires%CurrentSegment(i1)%TipoWire%L
         dummy%givenautoin = givenautoin
         dummy%resist =      resist
!
                    resist_devia =                     dummy%TipoWire%R_devia
               givenautoin_devia = HWires%CurrentSegment(i1)%TipoWire%L_devia
         dummy%givenautoin_devia =                          givenautoin_devia
              dummy%resist_devia =                               resist_devia

         !!bug'inest OLD 020413 !\E7 ahora solo se tratan resistencias y se aniaden a los segementos finales

         if ((dummy%Is_LeftEnd).and.(isLossy.or.isLossy)) then
            !no tengo en cuenta el caso particularisimo de un solo segmento conectado a lossy por los dos extremos !habria que sumarle la resistencia dos veces pero la casuistica se enfollona !\E7
            if ((.not.dummy%HasParallel_LeftEnd).and.(.not.dummy%HasSeries_LeftEnd).and.(.not.dummy%HasAbsorbing_LeftEnd)) then
               resist=resist+rlossy
               write (buff,'(a,4i7,a)')  'wir1_INFO: Adding Lossy material resistence to LeftEnd segment in contact with lossy without a terminal RLC ', &
               dummy%origIndex,dummy%i,dummy%j,dummy%k,' '//dir(dummy%tipofield)
            else
               resist=resist+rlossy
               write (buff,'(a,4i7,a)')  'wir1_INFO: Adding Lossy material resistence to LeftEnd segment grounded through RLC ', &
               dummy%origIndex,dummy%i,dummy%j,dummy%k,' '//dir(dummy%tipofield)
            endif
            if ((dummy%k >  ZI).and.(dummy%k <= ZE).and.(dummy%tipofield /= iEz).and.verbose) call WarnErrReport(buff)
            if ((dummy%k >= ZI).and.(dummy%k <= ZE).and.(dummy%tipofield == iEz).and.verbose) call WarnErrReport(buff)
         endif
         if ((dummy%Is_RightEnd).and.(isLossy.or.isLossy)) then
            !no tengo en cuenta el caso particularisimo de un solo segmento conectado a lossy por los dos extremos !habria que sumarle la resistencia dos veces pero la casuistica se enfollona !\E7
            if ((.not.dummy%HasParallel_RightEnd).and.(.not.dummy%HasSeries_RightEnd).and.(.not.dummy%HasAbsorbing_RightEnd)) then
               resist=resist+rlossy
               write (buff,'(a,4i7,a)')  'wir1_INFO: Adding Lossy material resistence to RightEnd segment in contact with lossy without a terminal RLC  ', &
               dummy%origIndex,dummy%i,dummy%j,dummy%k,' '//dir(dummy%tipofield)
            else
               resist=resist+rlossy
               write (buff,'(a,4i7,a)')  'wir1_INFO: Adding Lossy material resistence to RightEnd segment grounded through RLC ', &
               dummy%origIndex,dummy%i,dummy%j,dummy%k,' '//dir(dummy%tipofield)
            endif
            if ((dummy%k >  ZI).and.(dummy%k <= ZE).and.(dummy%tipofield /= iEz).and.verbose) call WarnErrReport(buff)
            if ((dummy%k >= ZI).and.(dummy%k <= ZE).and.(dummy%tipofield == iEz).and.verbose) call WarnErrReport(buff)
         endif
         !!! lOS QUE NO TENGAN RESITENCIAS Y ESTEN EN ABIERTO NO LOS CONECTO A LOSSY. sI SE QUIEREN HACER CONEXIONES A LOSSY
         !!! HAY QUE ESPECIFICAR UNA RESISTENCIA LUMPED
         !!! lUEGO SI SE CONECTARAN A pec DIRECTAMENTE SI LA TOPOLOGIA LO MANDA
         !!!
         if ((dummy%IsEnd_norLeft_norRight).and.(isLossy.or.isLossy)) then
            !no tengo en cuenta el caso particularisimo de un solo segmento conectado a lossy por los dos extremos !habria que sumarle la resistencia dos veces pero la casuistica se enfollona !\E7
            if ((.not.dummy%HasParallel_LeftEnd).and.(.not.dummy%HasSeries_LeftEnd).AND.(.not.dummy%HasParallel_RightEnd).and.(.not.dummy%HasSeries_RightEnd).and. &
                                                                                                              (.not.dummy%HasAbsorbing_RightEnd)) then
               resist=resist+rlossy
               write (buff,'(a,4i7,a)')  'wir1_INFO: Adding Lossy material resistence to Ending segment (other) segment in contact with lossy without a terminal RLC ', &
               dummy%origIndex,dummy%i,dummy%j,dummy%k,' '//dir(dummy%tipofield)
               if ((dummy%k >  ZI).and.(dummy%k <= ZE).and.(dummy%tipofield /= iEz).and.verbose) call WarnErrReport(buff)
               if ((dummy%k >= ZI).and.(dummy%k <= ZE).and.(dummy%tipofield == iEz).and.verbose) call WarnErrReport(buff)
            ELSE
               resist=resist+rlossy
               write (buff,'(a,4i7,a)')  'wir1_BUGGYERROR:  Lossy material resistence to Ending (other) segment grounded through RLC () ', &
               dummy%origIndex,dummy%i,dummy%j,dummy%k,' '//dir(dummy%tipofield)
               if ((dummy%k >  ZI).and.(dummy%k <= ZE).and.(dummy%tipofield /= iEz)) call WarnErrReport(buff,.true.)
               if ((dummy%k >= ZI).and.(dummy%k <= ZE).and.(dummy%tipofield == iEz)) call WarnErrReport(buff,.true.)
            ENDIF
         endif

         if (dummy%HasParallel_RightEnd) then
            givenautoin=givenautoin + dummy%TipoWire%Parallel_L_RightEnd/dummy%delta !se le suma la autoinduccion !2011 \E7 untested
            dummy%givenautoin=givenautoin
!stoch
                  givenautoin_devia=givenautoin_devia + dummy%TipoWire%Parallel_L_RightEnd_devia/dummy%delta !se le suma la autoinduccion !2011 \E7 untested
            dummy%givenautoin_devia=givenautoin_devia

            !!bug'inest OLD 020413 !\E7 ahora solo se tratan resistencias y se aniaden a los segmentos finales

            resist=resist +                dummy%TipoWire%Parallel_R_RightEnd/dummy%delta !se le suma la autoinduccion !2011 \E7 untested
            resist_devia=resist_devia +    dummy%TipoWire%Parallel_R_RightEnd_devia/dummy%delta !se le suma la autoinduccion !2011 \E7 untested
            if (dummy%TipoWire%Parallel_R_RightEnd /= 0.0_RKIND_wires) then
               write (buff,'(a,4i7,a)')  'wir1_INFO: Adding Parallel RightEnd Resistance in segment ', &
               dummy%origIndex,dummy%i,dummy%j,dummy%k,' '//dir(dummy%tipofield)
               if ((dummy%k >  ZI).and.(dummy%k <= ZE).and.(dummy%tipofield /= iEz).and.verbose) call WarnErrReport(buff)
               if ((dummy%k >= ZI).and.(dummy%k <= ZE).and.(dummy%tipofield == iEz).and.verbose) call WarnErrReport(buff)
            else
               write (buff,'(a,4i7,a)')  'wir1_INFO: Adding Parallel RightEnd null-Resistance in segment ', &
               dummy%origIndex,dummy%i,dummy%j,dummy%k,' '//dir(dummy%tipofield)
               if ((dummy%k >  ZI).and.(dummy%k <= ZE).and.(dummy%tipofield /= iEz).and.verbose) call WarnErrReport(buff)
               if ((dummy%k >= ZI).and.(dummy%k <= ZE).and.(dummy%tipofield == iEz).and.verbose) call WarnErrReport(buff)
            endif

            !(ojo que es per unit length la intrinsea)
            if (dummy%TipoWire%Parallel_L_RightEnd /= 0.0_RKIND_wires) then
               write (buff,'(a,4i7,a)')  'wir1_INFO: Adding Parallel RightEnd Inductance in segment ', &
               dummy%origIndex,dummy%i,dummy%j,dummy%k,' '//dir(dummy%tipofield)
               if ((dummy%k >  ZI).and.(dummy%k <= ZE).and.(dummy%tipofield /= iEz).and.verbose) call WarnErrReport(buff)
               if ((dummy%k >= ZI).and.(dummy%k <= ZE).and.(dummy%tipofield == iEz).and.verbose) call WarnErrReport(buff)
            endif
            !aniado tambien al ultimo segmento la resistencia y peto si hay capacitancias
            if (dummy%TipoWire%Parallel_C_RightEnd >= 1.0e-12_RKIND_wires) then
               write (buff,'(a,4i7,a)')  'wir1_ERROR: (Currently unsupported)  Capacitances in Parallel RightEnd at segment ', &
               dummy%origIndex,dummy%i,dummy%j,dummy%k,' '//dir(dummy%tipofield)
               if ((dummy%k >  ZI).and.(dummy%k <= ZE).and.(dummy%tipofield /= iEz)) call WarnErrReport(buff,.true.)
               if ((dummy%k >= ZI).and.(dummy%k <= ZE).and.(dummy%tipofield == iEz)) call WarnErrReport(buff,.true.)
            else
               dummy%TipoWire%Parallel_C_RightEnd=0.0_RKIND_wires
            endif
         endif
         if (dummy%HasParallel_LeftEnd) then
            givenautoin=givenautoin + dummy%TipoWire%Parallel_L_LeftEnd/dummy%delta !se le suma la autoinduccion !2011 \E7 untested
            dummy%givenautoin=givenautoin
!
                  givenautoin_devia=givenautoin_devia + dummy%TipoWire%Parallel_L_LeftEnd_devia/dummy%delta !se le suma la autoinduccion !2011 \E7 untested
            dummy%givenautoin_devia=givenautoin_devia
            !!bug'inest OLD 020413 !\E7 ahora solo se tratan resistencias y se aniaden a los segementos finales

            resist=resist +                  dummy%TipoWire%Parallel_R_LeftEnd      /dummy%delta !se le suma la autoinduccion !2011 \E7 untested
            resist_devia=resist_devia +      dummy%TipoWire%Parallel_R_LeftEnd_devia/dummy%delta !se le suma la autoinduccion !2011 \E7 untested
            if (dummy%TipoWire%Parallel_R_LeftEnd /= 0.0_RKIND_wires) then
               write (buff,'(a,4i7,a)')  'wir1_INFO: Adding Parallel LeftEnd Resistance in segment ', &
               dummy%origIndex,dummy%i,dummy%j,dummy%k,' '//dir(dummy%tipofield)
               if ((dummy%k >  ZI).and.(dummy%k <= ZE).and.(dummy%tipofield /= iEz).and.verbose) call WarnErrReport(buff)
               if ((dummy%k >= ZI).and.(dummy%k <= ZE).and.(dummy%tipofield == iEz).and.verbose) call WarnErrReport(buff)
            else
               write (buff,'(a,4i7,a)')  'wir1_INFO: Adding Parallel LeftEnd null-Resistance in segment ', &
               dummy%origIndex,dummy%i,dummy%j,dummy%k,' '//dir(dummy%tipofield)
               if ((dummy%k >  ZI).and.(dummy%k <= ZE).and.(dummy%tipofield /= iEz).and.verbose) call WarnErrReport(buff)
               if ((dummy%k >= ZI).and.(dummy%k <= ZE).and.(dummy%tipofield == iEz).and.verbose) call WarnErrReport(buff)
            endif

            if (dummy%TipoWire%Parallel_L_LeftEnd /= 0.0_RKIND_wires) then
               write (buff,'(a,4i7,a)')  'wir1_INFO: Adding Parallel LeftEnd Inductance in segment ', &
               dummy%origIndex,dummy%i,dummy%j,dummy%k,' '//dir(dummy%tipofield)
               if ((dummy%k >  ZI).and.(dummy%k <= ZE).and.(dummy%tipofield /= iEz).and.verbose) call WarnErrReport(buff)
               if ((dummy%k >= ZI).and.(dummy%k <= ZE).and.(dummy%tipofield == iEz).and.verbose) call WarnErrReport(buff)
            endif
            if (dummy%TipoWire%Parallel_C_LeftEnd >= 1.0e-12_RKIND_wires) then
               write (buff,'(a,4i7,a)')  'wir1_ERROR: (Currently unsupported)  Capacitances in Parallel LeftEnd at segment ', &
               dummy%origIndex,dummy%i,dummy%j,dummy%k,' '//dir(dummy%tipofield)
               if ((dummy%k >  ZI).and.(dummy%k <= ZE).and.(dummy%tipofield /= iEz)) call WarnErrReport(buff,.true.)
               if ((dummy%k >= ZI).and.(dummy%k <= ZE).and.(dummy%tipofield == iEz)) call WarnErrReport(buff,.true.)
            else
               dummy%TipoWire%Parallel_C_LeftEnd=0.0_RKIND_wires
            endif
         endif
         !
         if (dummy%HasSeries_RightEnd) then
            givenautoin=givenautoin + dummy%TipoWire%Series_L_RightEnd/dummy%delta !se le suma la autoinduccion !2011 \E7 untested
            resist=resist +           dummy%TipoWire%Series_R_RightEnd/dummy%delta 
            dummy%givenautoin=givenautoin
            dummy%resist=resist
!
                  givenautoin_devia=givenautoin_devia + dummy%TipoWire%Series_L_RightEnd_devia/dummy%delta !se le suma la autoinduccion !2011 \E7 untested
                       resist_devia=     resist_devia + dummy%TipoWire%Series_R_RightEnd_devia/dummy%delta 
            dummy%givenautoin_devia=givenautoin_devia
                 dummy%resist_devia=     resist_devia
            !(ojo que es per unit length la intrinsea)
            if (dummy%TipoWire%Series_L_RightEnd /= 0.0_RKIND_wires) then
               write (buff,'(a,4i7,a)')  'wir1_INFO: Adding Series RightEnd Inductance in segment ', &
               dummy%origIndex,dummy%i,dummy%j,dummy%k,' '//dir(dummy%tipofield)
               if ((dummy%k >  ZI).and.(dummy%k <= ZE).and.(dummy%tipofield /= iEz).and.verbose) call WarnErrReport(buff)
               if ((dummy%k >= ZI).and.(dummy%k <= ZE).and.(dummy%tipofield == iEz).and.verbose) call WarnErrReport(buff)
            endif
            if (dummy%TipoWire%Series_R_RightEnd /= 0.0_RKIND_wires) then
               write (buff,'(a,4i7,a)')  'wir1_INFO: Adding Series RightEnd Resistance in segment ', &
               dummy%origIndex,dummy%i,dummy%j,dummy%k,' '//dir(dummy%tipofield)
               if ((dummy%k >  ZI).and.(dummy%k <= ZE).and.(dummy%tipofield /= iEz).and.verbose) call WarnErrReport(buff)
               if ((dummy%k >= ZI).and.(dummy%k <= ZE).and.(dummy%tipofield == iEz).and.verbose) call WarnErrReport(buff)
            endif
            if (dummy%TipoWire%Series_C_RightEnd<= 1.0e7_RKIND_wires) then
               write (buff,'(a,4i7,a)')  'wir1_ERROR: (Currently unsupported)  Capacitances smaller than 1.0e7_RKIND_wires (inf) in Series RightEnd at segment ', &
               dummy%origIndex,dummy%i,dummy%j,dummy%k,' '//dir(dummy%tipofield)
               if ((dummy%k >  ZI).and.(dummy%k <= ZE).and.(dummy%tipofield /= iEz)) call WarnErrReport(buff,.true.)
               if ((dummy%k >= ZI).and.(dummy%k <= ZE).and.(dummy%tipofield == iEz)) call WarnErrReport(buff,.true.)
            else
               dummy%TipoWire%Series_C_RightEnd=2.0e7_RKIND_wires
            endif

         endif
         if (dummy%HasSeries_LeftEnd) then
            givenautoin=givenautoin + dummy%TipoWire%Series_L_LeftEnd/dummy%delta !se le suma la autoinduccion !2011 \E7 untested
            resist=resist +           dummy%TipoWire%Series_R_LeftEnd/dummy%delta !se le suma la autoinduccion !2011 \E7 untested
            dummy%givenautoin=givenautoin
            dummy%resist=resist
!
                   givenautoin_devia = givenautoin_devia + dummy%TipoWire%Series_L_LeftEnd_devia/dummy%delta !se le suma la autoinduccion !2011 \E7 untested
                        resist_devia =      resist_devia + dummy%TipoWire%Series_R_LeftEnd_devia/dummy%delta !se le suma la autoinduccion !2011 \E7 untested
             dummy%givenautoin_devia = givenautoin_devia
                  dummy%resist_devia =      resist_devia
            if (dummy%TipoWire%Series_L_LeftEnd /= 0.0_RKIND_wires) then
               write (buff,'(a,4i7,a)')  'wir1_INFO: Adding Series LeftEnd Inductance in segment ', &
               dummy%origIndex,dummy%i,dummy%j,dummy%k,' '//dir(dummy%tipofield)
               if ((dummy%k >  ZI).and.(dummy%k <= ZE).and.(dummy%tipofield /= iEz).and.verbose) call WarnErrReport(buff)
               if ((dummy%k >= ZI).and.(dummy%k <= ZE).and.(dummy%tipofield == iEz).and.verbose) call WarnErrReport(buff)
            endif
            if (dummy%TipoWire%Series_R_LeftEnd /= 0.0_RKIND_wires) then
               write (buff,'(a,4i7,a)')  'wir1_INFO: Adding Series LeftEnd Resistance in segment ', &
               dummy%origIndex,dummy%i,dummy%j,dummy%k,' '//dir(dummy%tipofield)
               if ((dummy%k >  ZI).and.(dummy%k <= ZE).and.(dummy%tipofield /= iEz).and.verbose) call WarnErrReport(buff)
               if ((dummy%k >= ZI).and.(dummy%k <= ZE).and.(dummy%tipofield == iEz).and.verbose) call WarnErrReport(buff)
            endif
            if (dummy%TipoWire%Series_C_LeftEnd <= 1.0e7_RKIND_wires) then
               write (buff,'(a,4i7,a)')  'wir1_ERROR: (Currently unsupported)  Capacitances  smaller than 1.0e7_RKIND_wires (inf) in Series LeftEnd atn segment ', &
               dummy%origIndex,dummy%i,dummy%j,dummy%k,' '//dir(dummy%tipofield)
               if ((dummy%k >  ZI).and.(dummy%k <= ZE).and.(dummy%tipofield /= iEz)) call WarnErrReport(buff,.true.)
               if ((dummy%k >= ZI).and.(dummy%k <= ZE).and.(dummy%tipofield == iEz)) call WarnErrReport(buff,.true.)
            else
               dummy%TipoWire%Series_C_LeftEnd=2.0e7_RKIND_wires
            endif

         endif 
    
         dummy%R = resist
         dummy%L = givenautoin
         dummy%givenautoin=givenautoin
         dummy%resist=resist
!
         dummy%R_devia = resist_devia
         dummy%L_devia = givenautoin_devia
         dummy%givenautoin_devia =givenautoin_devia
              dummy%resist_devia =     resist_devia
!!!ojooooo  acumulo en %lind toda la autoinduccion para que los calculos de capacidad la tengan en cuenta completa    
         if (.not.fieldtotl) then
             dummy%Lind=dummy%Lind+givenautoin
!
             dummy%Lind_devia=dummy%Lind_devia + givenautoin_devia
         else
             if (givenautoin <= tiny(1.0e-1) )  call StopOnError(0,0,'Fieldtotl not compatible with null given self inductance. ')
                 
             dummy%Lind=givenautoin
!
             dummy%Lind_devia= givenautoin_devia
         endif !para el fieldtotl no tento en cuenta mas que la autoin y no devuelvo nada !100517

         call wiresconstantes(fieldtotl,dummy,G2,sgg)


      end do

      !!!!
      !!!copiado 23/04/2014
      !
      if (trim(adjustl(wiresflavor))=='transition') then
         !contar multilines
         NumMultilines = 0
         do is1 = 1,HWires%NumCurrentSegments
            if(.not.(HWires%CurrentSegment(is1)%proc)) then
               NumMultilines = NumMultilines+1
               do is2 = is1+1,HWires%NumCurrentSegments
                  if ((HWires%CurrentSegment(is1)%i == HWires%CurrentSegment(is2)%i).and. &
                  (HWires%CurrentSegment(is1)%j == HWires%CurrentSegment(is2)%j).and. &
                  (HWires%CurrentSegment(is1)%k == HWires%CurrentSegment(is2)%k).and. &
                  (HWires%CurrentSegment(is1)%tipofield == HWires%CurrentSegment(is2)%tipofield).and. &
                  .not.(HWires%CurrentSegment(is1)%proc)) then

                     HWires%CurrentSegment(is2)%proc = .true.
                  end if
               end do
            end if
         end do

         do is1 = 1,HWires%NumCurrentSegments
            HWires%CurrentSegment(is1)%proc = .false.
         end do

         HWires%NumMultilines = NumMultilines
         write (buff,*)  'wir1_INFO: Numero de multilineas ', NumMultilines
         if (verbose) call WarnErrReport(buff)
         allocate(HWires%Multilines(1:NumMultilines))
         !contar paralelos
         contmtln = 0
         do is1 = 1,HWires%NumCurrentSegments
            if(.not.(HWires%CurrentSegment(is1)%proc)) then
               contmtln = contmtln + 1
               if(contmtln > NumMultilines) then
                  write (buff,*)  'wir0_BUGGYERROR: Demasiados multihilos'
                  if ((k1 >= ZI).and.(k1 <= ZE)) call WarnErrReport(buff,.true.)
                  ThereAreWires = .false.
               end if
               NumParallel = 1
               do is2 = is1+1,HWires%NumCurrentSegments
                  if ((HWires%CurrentSegment(is1)%i == HWires%CurrentSegment(is2)%i).and. &
                  (HWires%CurrentSegment(is1)%j == HWires%CurrentSegment(is2)%j).and. &
                  (HWires%CurrentSegment(is1)%k == HWires%CurrentSegment(is2)%k).and. &
                  (HWires%CurrentSegment(is1)%tipofield == HWires%CurrentSegment(is2)%tipofield)) then

                     NumParallel = NumParallel + 1
                     HWires%CurrentSegment(is2)%proc = .true.
                  end if
               end do
               HWires%Multilines(contmtln)%NumParallel = NumParallel
               allocate(HWires%Multilines(contmtln)%Segments(1:NumParallel))
            end if
         end do

         do is1 = 1,HWires%NumCurrentSegments
            HWires%CurrentSegment(is1)%proc = .false.
         end do

         !asignar multilines
         contmtln = 0
         do is1 = 1,HWires%NumCurrentSegments
            if(.not.(HWires%CurrentSegment(is1)%proc)) then
               contmtln = contmtln + 1
               if(contmtln > NumMultilines) then
                  write (buff,*)  'wir0_BUGGYERROR: Demasiados multihilos'
                  if ((k1 >= ZI).and.(k1 <= ZE)) call WarnErrReport(buff,.true.)
                  ThereAreWires = .false.
               end if
               contprll = 0
               do is2 = is1,HWires%NumCurrentSegments
                  if ((HWires%CurrentSegment(is1)%i == HWires%CurrentSegment(is2)%i).and. &
                  (HWires%CurrentSegment(is1)%j == HWires%CurrentSegment(is2)%j).and. &
                  (HWires%CurrentSegment(is1)%k == HWires%CurrentSegment(is2)%k).and. &
                  (HWires%CurrentSegment(is1)%tipofield == HWires%CurrentSegment(is2)%tipofield)) then

                     contprll = contprll + 1
                     if(contprll > HWires%Multilines(contmtln)%NumParallel) then
                        write (buff,*)  'wir0_BUGGYERROR: Demasiados hilos paralelos'
                        if ((k1 >= ZI).and.(k1 <= ZE)) call WarnErrReport(buff,.true.)
                        ThereAreWires = .false.
                     end if
                     HWires%CurrentSegment(is2)%proc = .true.
                     HWires%Multilines(contmtln)%Segments(contprll)%ptr => HWires%CurrentSegment(is2)
                  end if
               end do
            end if
         end do

         !asignamos posiciones en celda
         do iw1 = 1,NumMultilines
            N = HWires%Multilines(iw1)%NumParallel
            do i1 = 2,N
               HWires%Multilines(iw1)%Segments(i1)%ptr%x = HWires%Multilines(iw1)%Segments(i1)%ptr%x+cos(real(i1-2,KIND=RKIND_wires)/real(N-1,KIND=RKIND_wires)*2.0_RKIND_wires * pi)*0.25_RKIND_wires
               HWires%Multilines(iw1)%Segments(i1)%ptr%y = HWires%Multilines(iw1)%Segments(i1)%ptr%y+sin(real(i1-2,KIND=RKIND_wires)/real(N-1,KIND=RKIND_wires)*2.0_RKIND_wires * pi)*0.25_RKIND_wires
            end do
         end do

         !asignar autoinducciones
         !asignar constantes de evolucion

         do iw1 = 1,NumMultilines
            NumParallel = HWires%Multilines(iw1)%NumParallel
            allocate(HWires%Multilines(iw1)%C(1:NumParallel,1:NumParallel))
            allocate(HWires%Multilines(iw1)%R(1:NumParallel,1:NumParallel))
            allocate(HWires%Multilines(iw1)%L(1:NumParallel,1:NumParallel))
            do is1 = 1,NumParallel
               dl      = HWires%Multilines(iw1)%Segments(is1)%ptr%delta
               dx1 	= HWires%Multilines(iw1)%Segments(is1)%ptr%deltaTransv1
               dx2		= HWires%Multilines(iw1)%Segments(is1)%ptr%deltaTransv2
               r0 		= HWires%Multilines(iw1)%Segments(is1)%ptr%TipoWire%radius
               
               if (r0 < 1e-30) then 
                    write (buff,*)  'wir0_ERROR: Wire radius cannot be null'
                    if ((k1 >= ZI).and.(k1 <= ZE))  call WarnErrReport(buff,.true.)
               endif   
               imed 	= HWires%Multilines(iw1)%Segments(is1)%ptr%indexmed
               HWires%Multilines(iw1)%Segments(is1)%ptr%Lintrinsic  = &
               (1.0_RKIND_wires / (4.0_RKIND_wires * pi*InvMu(imed)))*(log((dx1**2.0_RKIND_wires +dx2**2.0_RKIND_wires)/(4.0_RKIND_wires * r0**2.0_RKIND_wires)) + &
               dx1/dx2*atan(dx2/dx1)     + &
               dx2/dx1*atan(dx1/dx2)     + &
               pi*r0**2.0_RKIND_wires / (dx2*dx1)-3.0_RKIND_wires)
               if ((r0 < 0.3_RKIND_wires  *dx1).or.(r0 < 0.3_RKIND_wires *dx2)) then
                  HWires%Multilines(iw1)%Segments(is1)%ptr%Lintrinsic = HWires%Multilines(iw1)%Segments(is1)%ptr%Lintrinsic - 0.57_RKIND_wires /(4.0_RKIND_wires * pi*InvMu(imed))

               endif
               if ((r0 > 0.3_RKIND_wires  *dx1).or.(r0 > 0.3_RKIND_wires *dx2)) then
                  HWires%Multilines(iw1)%Segments(is1)%ptr%Lintrinsic = HWires%Multilines(iw1)%Segments(is1)%ptr%Lintrinsic /(1.0_RKIND_wires-pi*r0**2.0_RKIND_wires /(dx2*dx1))
               endif

               HWires%Multilines(iw1)%Segments(is1)%ptr%L      = HWires%Multilines(iw1)%Segments(is1)%ptr%L + &
               HWires%Multilines(iw1)%Segments(is1)%ptr%Lintrinsic
               HWires%Multilines(iw1)%Segments(is1)%ptr%C      = HWires%Multilines(iw1)%Segments(is1)%ptr%C + &
               1.0_RKIND_wires / (InvMu(imed)*InvEps(imed)* &
               HWires%Multilines(iw1)%Segments(is1)%ptr%Lintrinsic)
               do is2 = 1,NumParallel
                  if(is1 == is2) then
                     HWires%Multilines(iw1)%C(is1,is2) = HWires%Multilines(iw1)%Segments(is1)%ptr%C
                     HWires%Multilines(iw1)%R(is1,is2) = HWires%Multilines(iw1)%Segments(is1)%ptr%R
                     HWires%Multilines(iw1)%L(is1,is2) = HWires%Multilines(iw1)%Segments(is1)%ptr%L
                  else
                     HWires%Multilines(iw1)%C(is1,is2) = 0
                     HWires%Multilines(iw1)%R(is1,is2) = 0

                     dist = sqrt(((HWires%Multilines(iw1)%Segments(is1)%ptr%x-HWires%Multilines(iw1)%Segments(is2)%ptr%x)*HWires%Multilines(iw1)%Segments(is2)%ptr%deltaTransv1)**2.0_RKIND_wires  + &
                     ((HWires%Multilines(iw1)%Segments(is1)%ptr%y-HWires%Multilines(iw1)%Segments(is2)%ptr%y)*HWires%Multilines(iw1)%Segments(is2)%ptr%deltaTransv2)**2.0_RKIND_wires )
                     phi = atan((HWires%Multilines(iw1)%Segments(is1)%ptr%y-HWires%Multilines(iw1)%Segments(is2)%ptr%y) / &
                     (HWires%Multilines(iw1)%Segments(is1)%ptr%x-HWires%Multilines(iw1)%Segments(is2)%ptr%x))
                     HWires%Multilines(iw1)%L(is1,is2) = (1.0_RKIND_wires / (2.0_RKIND_wires * pi*InvMu(imed)))*F(HWires%Multilines(iw1)%Segments(is1)%ptr%deltaTransv1/2.0_RKIND_wires, HWires%Multilines(iw1)%Segments(is1)%ptr%deltaTransv2/2.0_RKIND_wires, &
                     HWires%Multilines(iw1)%Segments(is1)%ptr%TipoWire%radius, HWires%Multilines(iw1)%Segments(is2)%ptr%TipoWire%radius, &
                     dist, phi)
                  end if
               end do
            end do
            write(buff,*) 'wir1_INFO: Multihilo ', iw1, ' con numero de hilos paralelos ', NumParallel
            if (verbose) call WarnErrReport(buff)
            allocate(HWires%Multilines(iw1)%b1I(1:NumParallel,1:NumParallel))
            allocate(HWires%Multilines(iw1)%b2I(1:NumParallel,1:NumParallel))
            allocate(HWires%Multilines(iw1)%b3I(1:NumParallel,1:NumParallel))
            allocate(Den                       (1:NumParallel,1:NumParallel))
            Den = HWires%Multilines(iw1)%L(1:NumParallel,1:NumParallel)+HWires%Multilines(iw1)%R(1:NumParallel,1:NumParallel)*sgg%dt/2.0_RKIND_wires
            call MatInv(NumParallel,Den)

            HWires%Multilines(iw1)%b1I = MatMul(Den, HWires%Multilines(iw1)%L(1:NumParallel,1:NumParallel)  - &
            HWires%Multilines(iw1)%R(1:NumParallel,1:NumParallel)*sgg%dt/2.0_RKIND_wires)
            !yo lo hago con cargas y no hay que multiplicar por C
            dl  =  HWires%Multilines(iw1)%Segments(1)%ptr%delta !tomo el de 1 !dama no lo tenia
            imed 	= HWires%Multilines(iw1)%Segments(1)%ptr%indexmed !lo aniado yo pq no estaba definido sgg 141118  !tomo el de 1 !dama no lo tenia
            HWires%Multilines(iw1)%b2I = -sgg%dt/dl*InvMu(imed)*InvEps(imed)*MatMul(Den,HWires%Multilines(iw1)%L)
            !   HWires%Multilines(iw1)%b2I = -sgg%dt/dl*InvMu(imed)*InvEps(imed)*MatMul(MatMul(Den,HWires%Multilines(iw1)%L),HWires%Multilines(iw1)%C)
            HWires%Multilines(iw1)%b3I = sgg%dt*Den

            deallocate(Den)

            do is1 = 1,NumParallel
               HWires%Multilines(iw1)%Segments(is1)%ptr%bI = HWires%Multilines(iw1)%b2I(is1,is1)
            end do
         end do
      ENDIF !DEL FLAVOR transition
      !!!!!!!!!!fin dama
      !!!


      !voltage sources
      do i1=1,HWires%NumCurrentSegments
         segmento=>HWires%CurrentSegment(i1)
         if (segmento%TipoWire%VsourceExists) then
            DO k1=1,segmento%Tipowire%numvoltagesources
               if ((segmento%i == segmento%Tipowire%Vsource(k1)%I).and. &
               (segmento%j == segmento%Tipowire%Vsource(k1)%J).and. &
               (segmento%k == segmento%Tipowire%Vsource(k1)%K)) then
                  HWires%CurrentSegment(i1)%Vsource=>segmento%TipoWire%Vsource(k1)
                  segmento%HasVsource=.true.
                  thereareVsources=.true.
                  if ((segmento%HasSeries_RightEnd).or.(segmento%HasSeries_LeftEnd)) then
                     write (buff,'(a,5i7)')  'wir1_INFO: Voltage source with Series RL (C neglected if present) impedance in segment ', &
                     segmento%origIndex,segmento%i,segmento%j,segmento%k,segmento%tipofield
                     if ((segmento%k >= ZI).and.(segmento%k <= ZE).and.verbose) call WarnErrReport(buff)
                  elseif ((segmento%HasParallel_RightEnd).or.(segmento%HasParallel_LeftEnd)) then
                     write (buff,'(a,5i7)')  'wir1_WARNING: Voltage source with Parallel RLC  in segment ', &
                     segmento%origIndex,segmento%i,segmento%j,segmento%k,segmento%tipofield
                     if ((segmento%k >= ZI).and.(segmento%k <= ZE)) call WarnErrReport(buff)
                     !
                  elseif ((segmento%HasAbsorbing_RightEnd).or.(segmento%HasAbsorbing_LeftEnd)) then
                     write (buff,'(a,5i7)')  'wir1_WARNING: Voltage source with Absorbing  in segment ', &
                     segmento%origIndex,segmento%i,segmento%j,segmento%k,segmento%tipofield
                     if ((segmento%k >= ZI).and.(segmento%k <= ZE)) call WarnErrReport(buff)
                     !
                  else
                     write (buff,'(a,5i7)')  'wir1_INFO: Voltage source with null internal resistence in segment ', &
                     segmento%origIndex,segmento%i,segmento%j,segmento%k,segmento%tipofield
                     if ((segmento%k >= ZI).and.(segmento%k <= ZE).and.verbose) call WarnErrReport(buff)
                  endif
                  if (segmento%Vsource%Fichero%DeltaSamples > sgg%dt) then
                     write (buff,'(a,a,a,e15.4e3)')  'wir1_WARNING: ',trim(adjustl(segmento%Vsource%Fichero%Name)), &
                     ' undersampled by a factor ', &
                     segmento%Vsource%Fichero%DeltaSamples/sgg%dt
                     if ((segmento%k >= ZI).and.(segmento%k <= ZE)) call WarnErrReport(buff)
                  endif
               endif
            end do
         endif
      end do

      !END OF THE PROCESSING OF CURRENT SEGMENT INFORMATION

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      !START OF CHARGE NODAL PROCESSING
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
      
      HWires%NumChargeNodes=2*HWires%NumCurrentSegments
      !There can be up to twice the ammount of segments (this allocation could be reduced later -I don't do it since memory is enough-)
      allocate (HWires%ChargeNode(HWires%NumChargeNodes))
      do i1=1,HWires%NumChargeNodes
         nullify ( HWires%ChargeNode(i1)%CurrentPlus_1, HWires%ChargeNode(i1)%CurrentPlus_2, HWires%ChargeNode(i1)%CurrentPlus_3, &
         HWires%ChargeNode(i1)%CurrentPlus_4, HWires%ChargeNode(i1)%CurrentPlus_5, HWires%ChargeNode(i1)%CurrentPlus_6, &
         HWires%ChargeNode(i1)%CurrentPlus_7, HWires%ChargeNode(i1)%CurrentPlus_8, HWires%ChargeNode(i1)%CurrentPlus_9, &
         HWires%ChargeNode(i1)%CurrentMinus_1,HWires%ChargeNode(i1)%CurrentMinus_2,HWires%ChargeNode(i1)%CurrentMinus_3, &
         HWires%ChargeNode(i1)%CurrentMinus_4,HWires%ChargeNode(i1)%CurrentMinus_5,HWires%ChargeNode(i1)%CurrentMinus_6, &
         HWires%ChargeNode(i1)%CurrentMinus_7,HWires%ChargeNode(i1)%CurrentMinus_8,HWires%ChargeNode(i1)%CurrentMinus_9, &
         HWires%ChargeNode(i1)%NodeInside)
         HWires%ChargeNode(i1)%indexnode           =-1
         HWires%ChargeNode(i1)%ChargePresent       =0.0_RKIND_wires
         HWires%ChargeNode(i1)%ChargePast          =0.0_RKIND_wires
         HWires%ChargeNode(i1)%IsAttachedtoVoltage =.false.
         HWires%ChargeNode(i1)%IsMur               =.false.
         HWires%ChargeNode(i1)%IsBackDownLeftMur =.false.
         HWires%ChargeNode(i1)%IsFrontUpRightMur =.false.
         HWires%ChargeNode(i1)%IsPeriodic               =.false.
         HWires%ChargeNode(i1)%IsPEC               =.false.
         HWires%ChargeNode(i1)%already_YEEadvanced_byconformal_changedtoPECfield1 => null()
         HWires%ChargeNode(i1)%already_YEEadvanced_byconformal_changedtoPECfield2 => null()
         HWires%ChargeNode(i1)%already_YEEadvanced_byconformal_changedtoPECfield3 => null()
         HWires%ChargeNode(i1)%already_YEEadvanced_byconformal_changedtoPECfield4 => null()
         HWires%ChargeNode(i1)%already_YEEadvanced_byconformal_changedtoPECfield5 => null()
         HWires%ChargeNode(i1)%already_YEEadvanced_byconformal_changedtoPECfield6 => null()
         HWires%ChargeNode(i1)%IsLossy             =.false.
         HWires%ChargeNode(i1)%HasIsource          =.false.
         HWires%ChargeNode(i1)%IsHeterogeneousJunction          =.false.
         !just for informative !not implemented in MPI unsure behaviour under mpi 2011 \E7
         HWires%ChargeNode(i1)%Exists              =.false.
         HWires%ChargeNode(i1)%Is_LeftEnd   =.false.
         HWires%ChargeNode(i1)%Is_RightEnd   =.false.
         HWires%ChargeNode(i1)%IsInSingleRLCsegment   =.false.
         HWires%ChargeNode(i1)%NumCurrentMinus =0
         HWires%ChargeNode(i1)%NumCurrentPlus  =0
         HWires%ChargeNode(i1)%i                   =-1
         HWires%ChargeNode(i1)%j                   =-1
         HWires%ChargeNode(i1)%k                   =-1
         HWires%ChargeNode(i1)%cteMur              =0.0_RKIND_wires
         HWires%ChargeNode(i1)%cteProp            =0.0_RKIND_wires
         HWires%ChargeNode(i1)%oRIGctePlain            =0.0_RKIND_wires
         HWires%ChargeNode(i1)%ctePlain            =0.0_RKIND_wires
         do j1=1,2*MaxNumCurrentMinusPlus
            HWires%ChargeNode(i1)%YESsegment(j1)      =-j1 !default !servira luego para contar nodos repetidos
         end do
      end do


      write (buff,'(a)')  '----------------------------------------------------------------'
      call WarnErrReport(buff)
      !Detect nodes which are Between two segments (posible duplicities are removed)
      conta=0
      NUMESEG=HWires%NumCurrentSegments
      !runtina de adjacencias
      do i1=1,HWires%NumCurrentSegments
         org=>HWires%CurrentSegment(i1)
         orgmenos1 => null()
         orgmas1 => null()
         if (i1>1) then
            IF (HWires%CurrentSegment(i1-1)%indexmed==HWires%CurrentSegment(i1)%indexmed) orgmenos1=>HWires%CurrentSegment(i1-1)
         endif
         if (i1<NUMESEG) then
            IF (HWires%CurrentSegment(i1+1)%indexmed==HWires%CurrentSegment(i1)%indexmed) orgmas1=>HWires%CurrentSegment(i1+1)
         endif
         do j1=1,HWires%NumCurrentSegments
            fin=>HWires%CurrentSegment(j1)
            if (J1>1) then
               FINmenos1=>HWires%CurrentSegment(J1-1)
            else
               FINmenos1 => null()
            endif
            if (J1<NUMESEG) then
               FINmas1=>HWires%CurrentSegment(J1+1)
            else
               FINmas1 => null()
            endif
            if (i1 /= j1) then
               !calls the adjacent routine to see wether the org and fin segments are connected (
               Adj=TestAdjacency(org,i1,fin,j1,connectendings,isolategroupgroups,strictOLD,ZI,ZE,NUMESEG,orgmenos1,orgmas1,finmenos1,finmas1,verbose)
               !
               If (Adj%Is) then
                  NodeExists=.false.
                  bprin: do nn=1,conta
                     if ((HWires%ChargeNode(nn)%i==adj%i).and. &
                     (HWires%ChargeNode(nn)%j==adj%j).and. &
                     (HWires%ChargeNode(nn)%k==adj%k)) then
                        if ((      adj%IsHeterogeneousJunction .and.      HWires%ChargeNode(nn)%IsHeterogeneousJunction).or. &
                        ((.not.adj%IsHeterogeneousJunction).and.(.not.HWires%ChargeNode(nn)%IsHeterogeneousJunction))) then
                           if ((adj%YESsegment(1) /= -1)) then
                              proceed1=.true.
                              proceed2=.true.
                              b1: do j2=1,2*MaxNumCurrentMinusPlus
                                 if ((HWires%ChargeNode(nn)%YESsegment(j2) == adj%YESsegment(1)))  then
                                    proceed1=.false.
                                    exit b1
                                 endif
                              end do b1
                              b11: do j2=1,2*MaxNumCurrentMinusPlus
                                 if ((HWires%ChargeNode(nn)%YESsegment(j2) == adj%YESsegment(2)))  then
                                    proceed2=.false.
                                    exit b11
                                 endif
                              end do b11
                              if ((.not.proceed1).or.(.not.proceed2).or.(.not.(proceed1.and.proceed2))) then
                                 !This node was already present. There is connectivity and the node actually exists
                                 NodeExists=.true.
                                 if (proceed1) then
                                    b2: do j2=1,2*MaxNumCurrentMinusPlus
                                       if (HWires%ChargeNode(nn)%YESsegment(j2) < 0 )   then
                                          HWires%ChargeNode(nn)%YESsegment(j2)  = adj%YESsegment(1)
                                          exit b2
                                       endif
                                    end do b2
                                 endif
                                 IF (PROCEED2) THEN
                                    b4: do j2=1,2*MaxNumCurrentMinusPlus
                                       if (HWires%ChargeNode(nn)%YESsegment(j2)  < 0 )   then
                                          HWires%ChargeNode(nn)%YESsegment(j2)  = adj%YESsegment(2)
                                          exit b4
                                       endif
                                    end do b4
                                 endif
                                 exit bprin
                              else
                                 NodeExists=.false.
                              endif
                           endif
                        endif
                     endif
                  end do bprin
                  if (.not.(NodeExists))   then
                     conta=conta+1
                     if ((adj%YESsegment(1) /= -1)) then
                        b6: do j2=1,2*MaxNumCurrentMinusPlus
                           if (HWires%ChargeNode(conta)%YESsegment(j2) < 0)   then
                              HWires%ChargeNode(conta)%YESsegment(j2)  = adj%YESsegment(1)
                              exit b6
                           endif
                        end do b6
                        b8: do j2=1,2*MaxNumCurrentMinusPlus
                           if (HWires%ChargeNode(conta)%YESsegment(j2) < 0)   then
                              HWires%ChargeNode(conta)%YESsegment(j2)  = adj%YESsegment(2)
                              exit b8
                           endif
                        end do b8
                     endif
                     !
                     HWires%ChargeNode(conta)%IsHeterogeneousJunction=adj%IsHeterogeneousJunction
                     HWires%ChargeNode(conta)%Exists=.true.
                     HWires%ChargeNode(conta)%i=adj%i
                     HWires%ChargeNode(conta)%j=adj%j
                     HWires%ChargeNode(conta)%k=adj%k
                  endif
               endif

               if (Adj%Is .and. ADj%BothEndingsConnected) then !busca el otro extremo y crealo tambien !SOLO NECESARIO EN WIRES DE UN SOLO SEGMENTO
                  NodeExists=.false.
                  if     ((adj%i /= org%i)) then
                     adj%i  = org%i
                  elseif ((adj%i /= org%ilibre)) then
                     adj%i  = org%ilibre
                  else
                     adj%i  = org%i
                  endif
                  if     ((adj%J /= org%J)) then
                     adj%J  = org%J
                  elseif ((adj%J /= org%Jlibre)) then
                     adj%J  = org%Jlibre
                  else
                     adj%J  = org%J
                  endif
                  if     ((adj%K /= org%K)) then
                     adj%K  = org%K
                  elseif ((adj%K /= org%Klibre)) then
                     adj%K  = org%Klibre
                  else
                     adj%K  = org%K
                  endif
                  bprin2: do nn=1,conta
                     if ((HWires%ChargeNode(nn)%i==adj%i).and. &
                     (HWires%ChargeNode(nn)%j==adj%j).and. &
                     (HWires%ChargeNode(nn)%k==adj%k)) then
                        if ((      adj%IsHeterogeneousJunction .and.      HWires%ChargeNode(nn)%IsHeterogeneousJunction).or. &
                        ((.not.adj%IsHeterogeneousJunction).and.(.not.HWires%ChargeNode(nn)%IsHeterogeneousJunction))) then
                           if ((adj%YESsegment(1) /= -1)) then
                              proceed1=.true.
                              proceed2=.true.
                              b12: do j2=1,2*MaxNumCurrentMinusPlus
                                 if ((HWires%ChargeNode(nn)%YESsegment(j2) == adj%YESsegment(1)))  then
                                    proceed1=.false.
                                    exit b12
                                 endif
                              end do b12
                              b112: do j2=1,2*MaxNumCurrentMinusPlus
                                 if ((HWires%ChargeNode(nn)%YESsegment(j2) == adj%YESsegment(2)))  then
                                    proceed2=.false.
                                    exit b112
                                 endif
                              end do b112
                              if ((.not.proceed1).or.(.not.proceed2).or.(.not.(proceed1.and.proceed2))) then
                                 !This node was already present. There is connectivity and the node actually exists
                                 NodeExists=.true.
                                 if (proceed1) then
                                    b22: do j2=1,2*MaxNumCurrentMinusPlus
                                       if (HWires%ChargeNode(nn)%YESsegment(j2) < 0 )   then
                                          HWires%ChargeNode(nn)%YESsegment(j2)  = adj%YESsegment(1)
                                          exit b22
                                       endif
                                    end do b22
                                 endif
                                 IF (PROCEED2) THEN
                                    b42: do j2=1,2*MaxNumCurrentMinusPlus
                                       if (HWires%ChargeNode(nn)%YESsegment(j2)  < 0 )   then
                                          HWires%ChargeNode(nn)%YESsegment(j2)  = adj%YESsegment(2)
                                          exit b42
                                       endif
                                    end do b42
                                 endif
                                 exit bprin2
                              else
                                 NodeExists=.false.
                              endif
                           endif
                        endif
                     endif
                  end do bprin2
                  if (.not.(NodeExists))   then
                     conta=conta+1
                     if ((adj%YESsegment(1) /= -1)) then
                        b62: do j2=1,2*MaxNumCurrentMinusPlus
                           if (HWires%ChargeNode(conta)%YESsegment(j2) < 0)   then
                              HWires%ChargeNode(conta)%YESsegment(j2)  = adj%YESsegment(1)
                              exit b62
                           endif
                        end do b62
                        b82: do j2=1,2*MaxNumCurrentMinusPlus
                           if (HWires%ChargeNode(conta)%YESsegment(j2) < 0)   then
                              HWires%ChargeNode(conta)%YESsegment(j2)  = adj%YESsegment(2)
                              exit b82
                           endif
                        end do b82
                     endif
                     !
                     HWires%ChargeNode(conta)%IsHeterogeneousJunction=adj%IsHeterogeneousJunction
                     HWires%ChargeNode(conta)%Exists=.true.
                     HWires%ChargeNode(conta)%i=adj%i
                     HWires%ChargeNode(conta)%j=adj%j
                     HWires%ChargeNode(conta)%k=adj%k
                  endif
               endif
            endif !del i1/=ji1
            !
         end do
      end do
      HWires%NumChargeNodes=conta
      !reactualize the number of nodes


      write (buff,'(a)')  '----------------------------------------------------------------'
      call WarnErrReport(buff)

      !Now I remove duplicate nodes, using info of the shared segments
      !More than one of such node may exist if there are commond adjacencies
      !They are not sequentially detected
      do nn=1,HWires%NumChargeNodes
         do nnn=nn+1,HWires%NumChargeNodes
            if ((HWires%ChargeNode(nn)%i==HWires%ChargeNode(nnn)%i).and. &
            (HWires%ChargeNode(nn)%j==HWires%ChargeNode(nnn)%j).and. &
            (HWires%ChargeNode(nn)%k==HWires%ChargeNode(nnn)%k)) then
               conta=0
               do j1=1,2*MaxNumCurrentMinusPlus
                  do k1=1,2*MaxNumCurrentMinusPlus
                     if (HWires%ChargeNode(nn)%YESsegment(k1) == HWires%ChargeNode(nnn)%YESsegment(j1)) then
                        conta=conta+1
                     endif
                  end do
               end do
               if (conta == 2*MaxNumCurrentMinusPlus) then !They coincide
                  HWires%ChargeNode(nn)%Exists=.false. !voids it
               endif
            endif
         end do
      end do

      !compress the authentic nodes
      conta=HWires%NumChargeNodes
      do nn=1,HWires%NumChargeNodes
         if (.not.HWires%ChargeNode(nn)%Exists) then
            conta=conta-1
            do nnn=nn,HWires%NumChargeNodes-1
               HWires%ChargeNode(nnn) = HWires%ChargeNode(nnn+1)
            end do
            continue
         endif
      end do
      !voids the rest 
      do i1=conta+1,HWires%NumChargeNodes

         nullify (HWires%ChargeNode(i1)%CurrentPlus_1, HWires%ChargeNode(i1)%CurrentPlus_2, HWires%ChargeNode(i1)%CurrentPlus_3, &
         HWires%ChargeNode(i1)%CurrentPlus_4, HWires%ChargeNode(i1)%CurrentPlus_5, HWires%ChargeNode(i1)%CurrentPlus_6, &
         HWires%ChargeNode(i1)%CurrentPlus_7, HWires%ChargeNode(i1)%CurrentPlus_8, HWires%ChargeNode(i1)%CurrentPlus_9, &
         HWires%ChargeNode(i1)%CurrentMinus_1,HWires%ChargeNode(i1)%CurrentMinus_2,HWires%ChargeNode(i1)%CurrentMinus_3, &
         HWires%ChargeNode(i1)%CurrentMinus_4,HWires%ChargeNode(i1)%CurrentMinus_5,HWires%ChargeNode(i1)%CurrentMinus_6, &
         HWires%ChargeNode(i1)%CurrentMinus_7,HWires%ChargeNode(i1)%CurrentMinus_8,HWires%ChargeNode(i1)%CurrentMinus_9, &
         HWires%ChargeNode(i1)%NodeInside)
         HWires%ChargeNode(i1)%ChargePresent       =0.0_RKIND_wires
         HWires%ChargeNode(i1)%ChargePast          =0.0_RKIND_wires
         HWires%ChargeNode(i1)%IsAttachedtoVoltage =.false.
         HWires%ChargeNode(i1)%IsMur               =.false.
         HWires%ChargeNode(i1)%IsBackDownLeftMur =.false.
         HWires%ChargeNode(i1)%IsFrontUpRightMur =.false.
         HWires%ChargeNode(i1)%Isperiodic               =.false.
         HWires%ChargeNode(i1)%IsPEC               =.false.
         HWires%ChargeNode(i1)%already_YEEadvanced_byconformal_changedtoPECfield1 => null()
         HWires%ChargeNode(i1)%already_YEEadvanced_byconformal_changedtoPECfield2 => null()
         HWires%ChargeNode(i1)%already_YEEadvanced_byconformal_changedtoPECfield3 => null()
         HWires%ChargeNode(i1)%already_YEEadvanced_byconformal_changedtoPECfield4 => null()
         HWires%ChargeNode(i1)%already_YEEadvanced_byconformal_changedtoPECfield5 => null()
         HWires%ChargeNode(i1)%already_YEEadvanced_byconformal_changedtoPECfield6 => null()
         HWires%ChargeNode(i1)%IsLossy             =.false.
         HWires%ChargeNode(i1)%HasIsource          =.false.
         HWires%ChargeNode(i1)%IsHeterogeneousJunction          =.false.
         HWires%ChargeNode(i1)%Exists              =.false.
         HWires%ChargeNode(i1)%Is_LeftEnd   =.false.
         HWires%ChargeNode(i1)%Is_RightEnd   =.false.
         HWires%ChargeNode(i1)%NumCurrentMinus =0
         HWires%ChargeNode(i1)%NumCurrentPlus  =0
         HWires%ChargeNode(i1)%i                   =-1
         HWires%ChargeNode(i1)%j                   =-1
         HWires%ChargeNode(i1)%k                   =-1
         HWires%ChargeNode(i1)%cteMur              =0.0_RKIND_wires
         HWires%ChargeNode(i1)%cteProp            =0.0_RKIND_wires
         HWires%ChargeNode(i1)%oRIGctePlain            =0.0_RKIND_wires
         HWires%ChargeNode(i1)%ctePlain            =0.0_RKIND_wires
         do j1=1,2*MaxNumCurrentMinusPlus
            HWires%ChargeNode(i1)%YESsegment(j1)      =-j1 !default !used for duplicate nodes
         end do
      end do
      !
      HWires%NumChargeNodes=conta !reactualize the number of authentic nodes


      !I point segments nodal info to the nodes connected to each segment
      do i1=1,HWires%NumChargeNodes
         nodo=>HWires%ChargeNode(i1)
         do j1=1,HWires%NumCurrentSegments
            segmento=>HWires%CurrentSegment(j1)
            proceed=.false.
            total: do k1=1,2*MaxNumCurrentMinusPlus
               if (j1 == nodo%YESsegment(k1))  then
                  proceed=.true.
                  exit total
               endif
            end do total
            if (proceed) then
               if  ((segmento%i  ==nodo%i).and.(segmento%j==nodo%j).and.(segmento%k==nodo%k)) then
                  segmento%ChargeMinus=>HWires%ChargeNode(i1)
               endif
               if  ((segmento%i+1==nodo%i).and.(segmento%j==nodo%j).and.(segmento%k==nodo%k).and.(segmento%tipofield==iEx)) then
                  segmento%ChargePlus=>HWires%ChargeNode(i1)
               endif
               if  ((segmento%i==nodo%i).and.(segmento%j+1==nodo%j).and.(segmento%k==nodo%k).and.(segmento%tipofield==iEy)) then
                  segmento%ChargePlus=>HWires%ChargeNode(i1)
               endif
               if  ((segmento%i==nodo%i).and.(segmento%j==nodo%j).and.(segmento%k+1==nodo%k).and.(segmento%tipofield==iEz)) then
                  segmento%ChargePlus=>HWires%ChargeNode(i1)
               endif
            else
               continue
            endif
         end do
      end do

      !I point nodes to their segments, detecting possible junctions (I set the IsJunction Flag)
      !Beware that I set the IsJunction flag even in the case of TWO segments which are both in the plus
      !or in the minus direction (later I set another flag for the authentic junctions (up to 2*MaxNumCurrentMinusPlus) segments
      do j1=1,HWires%NumCurrentSegments
         segmento=>HWires%CurrentSegment(j1)
         if (associated(segmento%ChargePlus)) then
            segmento%ChargePlus%NumCurrentMinus=segmento%ChargePlus%NumCurrentMinus+1
            if (segmento%ChargePlus%NumCurrentMinus >9) then
               BUFF='wir1_ERROR: More than 9 Minus WIREs joined'
               call WarnErrReport(buff,.true.)
            endif
            select case (segmento%ChargePlus%NumCurrentMinus)
             case (1)
               segmento%ChargePlus%CurrentMinus_1=>HWires%CurrentSegment(j1)
             case (2)
               segmento%ChargePlus%CurrentMinus_2=>HWires%CurrentSegment(j1)
             case (3)
               segmento%ChargePlus%CurrentMinus_3=>HWires%CurrentSegment(j1)
             case (4)
               segmento%ChargePlus%CurrentMinus_4=>HWires%CurrentSegment(j1)
             case (5)
               segmento%ChargePlus%CurrentMinus_5=>HWires%CurrentSegment(j1)
             case (6)
               segmento%ChargePlus%CurrentMinus_6=>HWires%CurrentSegment(j1)
             case (7)
               segmento%ChargePlus%CurrentMinus_7=>HWires%CurrentSegment(j1)
             case (8)
               segmento%ChargePlus%CurrentMinus_8=>HWires%CurrentSegment(j1)
             case (9)
               segmento%ChargePlus%CurrentMinus_9=>HWires%CurrentSegment(j1)
            end select
         endif
         if (associated(segmento%ChargeMinus)) then
            segmento%ChargeMinus%NumCurrentPlus=segmento%ChargeMinus%NumCurrentPlus+1
            if (segmento%ChargeMinus%NumCurrentPlus >9) then
               BUFF='wir1_ERROR: More than 9 Plus WIREs joined'
               call WarnErrReport(buff,.true.)
            endif
            select case (segmento%ChargeMinus%NumCurrentPlus)
             case (1)
               segmento%ChargeMinus%CurrentPlus_1=>HWires%CurrentSegment(j1)
             case (2)
               segmento%ChargeMinus%CurrentPlus_2=>HWires%CurrentSegment(j1)
             case (3)
               segmento%ChargeMinus%CurrentPlus_3=>HWires%CurrentSegment(j1)
             case (4)
               segmento%ChargeMinus%CurrentPlus_4=>HWires%CurrentSegment(j1)
             case (5)
               segmento%ChargeMinus%CurrentPlus_5=>HWires%CurrentSegment(j1)
             case (6)
               segmento%ChargeMinus%CurrentPlus_6=>HWires%CurrentSegment(j1)
             case (7)
               segmento%ChargeMinus%CurrentPlus_7=>HWires%CurrentSegment(j1)
             case (8)
               segmento%ChargeMinus%CurrentPlus_8=>HWires%CurrentSegment(j1)
             case (9)
               segmento%ChargeMinus%CurrentPlus_9=>HWires%CurrentSegment(j1)
            end select
         endif
      end do

      !Now I detect hanging nodes (at wire terminations) and increment the node info
      do i1=1,HWires%NumCurrentSegments
         if (.not.associated(HWires%CurrentSegment(i1)%ChargePlus)) then
            !segments must end in a charge node
            conta=conta+1
            HWires%CurrentSegment(i1)%ChargePlus=>HWires%ChargeNode(conta)
            HWires%ChargeNode(conta)%CurrentMinus_1=>HWires%CurrentSegment(i1)
            HWires%ChargeNode(conta)%NumCurrentMinus=1
            HWires%ChargeNode(conta)%i=HWires%ChargeNode(conta)%CurrentMinus_1%i
            HWires%ChargeNode(conta)%j=HWires%ChargeNode(conta)%CurrentMinus_1%j
            HWires%ChargeNode(conta)%k=HWires%ChargeNode(conta)%CurrentMinus_1%k
            HWires%ChargeNode(conta)%Exists=.true.
            If (HWires%ChargeNode(conta)%CurrentMinus_1%tipofield==iEx) then
               HWires%ChargeNode(conta)%i=HWires%ChargeNode(conta)%CurrentMinus_1%i+1
            endif
            If (HWires%ChargeNode(conta)%CurrentMinus_1%tipofield==iEy) then
               HWires%ChargeNode(conta)%j=HWires%ChargeNode(conta)%CurrentMinus_1%j+1
            endif
            If (HWires%ChargeNode(conta)%CurrentMinus_1%tipofield==iEz) then
               HWires%ChargeNode(conta)%k=HWires%ChargeNode(conta)%CurrentMinus_1%k+1
            endif
         endif
         if (.not.associated(HWires%CurrentSegment(i1)%ChargeMinus)) then
            !segments must end in a charge node
            conta=conta+1
            HWires%CurrentSegment(i1)%ChargeMinus=>HWires%ChargeNode(conta)
            HWires%ChargeNode(conta)%CurrentPlus_1=>HWires%CurrentSegment(i1)
            HWires%ChargeNode(conta)%NumCurrentPlus=1
            HWires%ChargeNode(conta)%i=HWires%ChargeNode(conta)%CurrentPlus_1%i
            HWires%ChargeNode(conta)%j=HWires%ChargeNode(conta)%CurrentPlus_1%j
            HWires%ChargeNode(conta)%k=HWires%ChargeNode(conta)%CurrentPlus_1%k
            HWires%ChargeNode(conta)%Exists=.true.
            If (HWires%ChargeNode(conta)%CurrentPlus_1%tipofield==iEx) then
               HWires%ChargeNode(conta)%i=HWires%ChargeNode(conta)%CurrentPlus_1%i
            endif
            If (HWires%ChargeNode(conta)%CurrentPlus_1%tipofield==iEy) then
               HWires%ChargeNode(conta)%j=HWires%ChargeNode(conta)%CurrentPlus_1%j
            endif
            If (HWires%ChargeNode(conta)%CurrentPlus_1%tipofield==iEz) then
               HWires%ChargeNode(conta)%k=HWires%ChargeNode(conta)%CurrentPlus_1%k
            endif
         endif
      end do
      HWires%NumChargeNodes=conta !reactualize POR ULTIMA VEZ the number of charge nodes
      !!!!!!!!!!!


      do i1=1,HWires%NumCurrentSegments
         segmento=>HWires%CurrentSegment(i1)
         asignado=.false.
         if (segmento%Is_LeftEnd.or.segmento%IsEnd_norLeft_norRight) then
            if ((segmento%chargeplus%i == segmento%ilibre).and. &
            (segmento%chargeplus%j == segmento%jlibre).and. &
            (segmento%chargeplus%k == segmento%klibre)) then
               segmento%chargeplus%Is_LeftEnd =.true.
               asignado=.true.
            endif
            if ((segmento%chargeMinus%i == segmento%ilibre).and. &
            (segmento%chargeMinus%j == segmento%jlibre).and. &
            (segmento%chargeMinus%k == segmento%klibre)) then
               segmento%chargeMinus%Is_LeftEnd =.true.
               asignado=.true.
            endif
         endif
         if (segmento%Is_RightEnd) then
            if ((segmento%chargeplus%i == segmento%ilibre).and. &
            (segmento%chargeplus%j == segmento%jlibre).and. &
            (segmento%chargeplus%k == segmento%klibre)) then
               segmento%chargeplus%Is_RightEnd =.true.
               asignado=.true.
            endif
            if ((segmento%chargeMinus%i == segmento%ilibre).and. &
            (segmento%chargeMinus%j == segmento%jlibre).and. &
            (segmento%chargeMinus%k == segmento%klibre)) then
               segmento%chargeMinus%Is_RightEnd =.true.
               asignado=.true.
            endif
         endif
         !caso particular hilos de 1  segmento
         if (segmento%Is_LeftEnd.and.segmento%Is_RightEnd) then
            segmento%chargeplus%Is_LeftEnd =.true.
            segmento%chargeMinus%Is_RightEnd =.true.
            asignado=.true.
            !bug OLD 060313 hilos de un solo segmento con RLC. Correcion en caso de formatos antiguos solo
            if (connectendings) then
               segmento%chargeplus%Is_LeftEnd  =.true.
               segmento%chargeplus%Is_RightEnd  =.true.
               segmento%chargeMinus%Is_RightEnd =.true.
               segmento%chargeMinus%Is_LeftEnd =.true.
            endif
            segmento%chargeplus%IsInSingleRLCsegment  = (segmento%HasParallel_LeftEnd.or.segmento%HasSeries_LeftEnd).or. &
            (segmento%HasParallel_RightEnd.or.segmento%HasSeries_RightEnd)
            segmento%chargeminus%IsInSingleRLCsegment = (segmento%HasParallel_LeftEnd.or.segmento%HasSeries_LeftEnd).or. &
            (segmento%HasParallel_RightEnd.or.segmento%HasSeries_RightEnd)
         endif
         if (.not.asignado) then
            if (segmento%HasParallel_LeftEnd .or. segmento%HasSeries_LeftEnd ) then
               if (segmento%chargeplus%isPEC .or. segmento%chargeplus%isLossy.or. segmento%chargeplus%isLossy.or. segmento%chargeplus%isLossy) then
                  if ((.not. segmento%chargeplus%Is_LeftEnd).and.(.not. segmento%chargeplus%Is_RightEnd)) then
                     segmento%chargeplus%Is_LeftEnd =.true.
                     write (buff,'(a,4i7)')  'wir1_INFO: Forcing non-terminal node to LeftEnd to attach RLC ', &
                     segmento%chargeplus%i,segmento%chargeplus%j,segmento%chargeplus%k
                     if ((segmento%chargeplus%k >  ZI).and.(segmento%chargeplus%k <= ZE).and.verbose) call WarnErrReport(buff)
                  endif
               elseif (segmento%chargeminus%isPEC .or. segmento%chargeminus%isLossy .or. segmento%chargeminus%isLossy .or. segmento%chargeminus%isLossy) then
                  if ((.not. segmento%chargeminus%Is_LeftEnd).and.(.not. segmento%chargeminus%Is_RightEnd)) then
                     segmento%chargeminus%Is_LeftEnd =.true.
                     write (buff,'(a,4i7)')  'wir1_INFO: Forcing non-terminal node to LeftEnd to attach RLC ', &
                     segmento%chargeminus%i,segmento%chargeminus%j,segmento%chargeminus%k
                     if ((segmento%chargeminus%k >  ZI).and.(segmento%chargeminus%k <= ZE).and.verbose) call WarnErrReport(buff)
                  endif
               else
                  if (segmento%TipoWire%Series_C_LeftEnd < 1.0e7_RKIND_wires) then
                     write (buff,'(a,4i7,a)')  'wir1_ERROR: Series LeftEnd Capacitance in INTERMEDIATE segment smaller than 1e7 (inf)', &
                     segmento%origIndex,segmento%i,segmento%j,segmento%k,dir(segmento%tipofield)
                     if ((segmento%k >  ZI).and.(segmento%k <= ZE).and.(segmento%tipofield /= iEz)) call WarnErrReport(buff,.TRUE.)
                     if ((segmento%k >= ZI).and.(segmento%k <= ZE).and.(segmento%tipofield == iEz)) call WarnErrReport(buff,.TRUE.)
                  endif
                  if (segmento%TipoWire%Parallel_C_LeftEnd /= 0.0_RKIND_wires) then
                     write (buff,'(a,4i7,a)')  'wir1_ERROR: Parallel LeftEnd Capacitance in INTERMEDIATE segment ', &
                     segmento%origIndex,segmento%i,segmento%j,segmento%k,dir(segmento%tipofield)
                     if ((segmento%k >  ZI).and.(segmento%k <= ZE).and.(segmento%tipofield /= iEz)) call WarnErrReport(buff,.TRUE.)
                     if ((segmento%k >= ZI).and.(segmento%k <= ZE).and.(segmento%tipofield == iEz)) call WarnErrReport(buff,.TRUE.)
                  endif
                  if (segmento%TipoWire%Parallel_R_LeftEnd /= 0.0_RKIND_wires) then
                     write (buff,'(a,4i7,a)')  'wir1_ERROR:  Parallel LeftEnd Resistance in INTERMEDIATE segment ', &
                     segmento%origIndex,segmento%i,segmento%j,segmento%k,dir(segmento%tipofield)
                     if ((segmento%k >  ZI).and.(segmento%k <= ZE).and.(segmento%tipofield /= iEz)) call WarnErrReport(buff,.TRUE.)
                     if ((segmento%k >= ZI).and.(segmento%k <= ZE).and.(segmento%tipofield == iEz)) call WarnErrReport(buff,.TRUE.)
                  endif
               endif
            endif
            if (segmento%HasParallel_RightEnd .or. segmento%HasSeries_RightEnd ) then
               if (segmento%chargeplus%isPEC .or. segmento%chargeplus%isLossy .or. segmento%chargeplus%isLossy .or. segmento%chargeplus%isLossy) then
                  if ((.not. segmento%chargeplus%Is_LeftEnd).and.(.not. segmento%chargeplus%Is_RightEnd)) then
                     segmento%chargeplus%Is_RightEnd =.true.
                     write (buff,'(a,4i7)')  'wir1_INFO: Forcing non-terminal node to RightEnd to attach RLC ', &
                     segmento%chargeplus%i,segmento%chargeplus%j,segmento%chargeplus%k
                     if ((segmento%chargeplus%k >  ZI).and.(segmento%chargeplus%k <= ZE).and.verbose) call WarnErrReport(buff)
                  endif
               elseif (segmento%chargeminus%isPEC .or. segmento%chargeminus%isLossy .or. segmento%chargeminus%isLossy .or. segmento%chargeminus%isLossy) then
                  if ((.not. segmento%chargeminus%Is_LeftEnd).and.(.not. segmento%chargeminus%Is_RightEnd)) then
                     segmento%chargeminus%Is_RightEnd =.true.
                     write (buff,'(a,4i7)')  'wir1_INFO: Forcing non-terminal node to RightEnd to attach RLC ', &
                     segmento%chargeminus%i,segmento%chargeminus%j,segmento%chargeminus%k
                     if ((segmento%chargeminus%k >  ZI).and.(segmento%chargeminus%k <= ZE).and.verbose) call WarnErrReport(buff)
                  endif
               else
                  if (segmento%TipoWire%Series_C_RightEnd < 1.0e7_RKIND_wires) then
                     write (buff,'(a,4i7,a)')  'wir1_ERROR: Series RightEnd Capacitance in INTERMEDIATE segment  smaller than 1e7 (inf)', &
                     segmento%origIndex,segmento%i,segmento%j,segmento%k,dir(segmento%tipofield)
                     if ((segmento%k >  ZI).and.(segmento%k <= ZE).and.(segmento%tipofield /= iEz)) call WarnErrReport(buff,.TRUE.)
                     if ((segmento%k >= ZI).and.(segmento%k <= ZE).and.(segmento%tipofield == iEz)) call WarnErrReport(buff,.TRUE.)
                  endif
                  if (segmento%TipoWire%Parallel_C_RightEnd /= 0.0_RKIND_wires) then
                     write (buff,'(a,4i7,a)')  'wir1_ERROR: Parallel RightEnd Capacitance in INTERMEDIATE segment ', &
                     segmento%origIndex,segmento%i,segmento%j,segmento%k,dir(segmento%tipofield)
                     if ((segmento%k >  ZI).and.(segmento%k <= ZE).and.(segmento%tipofield /= iEz)) call WarnErrReport(buff,.TRUE.)
                     if ((segmento%k >= ZI).and.(segmento%k <= ZE).and.(segmento%tipofield == iEz)) call WarnErrReport(buff,.TRUE.)
                  endif
                  if (segmento%TipoWire%Parallel_R_RightEnd /= 0.0_RKIND_wires) then
                     write (buff,'(a,4i7,a)')  'wir1_ERROR: Series RightEnd Resistance in INTERMEDIATE segment ', &
                     segmento%origIndex,segmento%i,segmento%j,segmento%k,dir(segmento%tipofield)
                     if ((segmento%k >  ZI).and.(segmento%k <= ZE).and.(segmento%tipofield /= iEz)) call WarnErrReport(buff,.TRUE.)
                     if ((segmento%k >= ZI).and.(segmento%k <= ZE).and.(segmento%tipofield == iEz)) call WarnErrReport(buff,.TRUE.)
                  endif
               endif
            endif
         endif
      end do

      call detect_peclossyconformal_nodes

      !END of geometrical PREPROCESSING of junctions and adjacencies

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Adjust constants and sort of nodes

      !find PMC segments
      do i1=1,HWires%NumCurrentSegments
         segmento=>HWires%CurrentSegment(i1)
         select case (segmento%tipofield)
          case (iEx)
            if ((segmento%i+1==sgg%Alloc(segmento%tipofield)%XE).and.(sgg%Border%IsFrontPMC)) then
               segmento%IsPMC =.true.
            endif
            if ((segmento%i==sgg%Alloc(segmento%tipofield)%XI).and.(sgg%Border%IsBackPMC)) then
               segmento%IsPMC =.true.
            endif
          case (iEy)
            if ((segmento%j+1==sgg%Alloc(segmento%tipofield)%YE).and.(sgg%Border%IsRightPMC)) then
               segmento%IsPMC =.true.
            endif
            if ((segmento%j==sgg%Alloc(segmento%tipofield)%YI).and.(sgg%Border%IsLeftPMC)) then
               segmento%IsPMC =.true.
            endif
          case (iEz)
            if ((segmento%k+1==sgg%Alloc(segmento%tipofield)%ZE).and.(sgg%Border%IsUpPMC)) then
               segmento%IsPMC =.true.
            endif
            if ((segmento%k==sgg%Alloc(segmento%tipofield)%ZI).and.(sgg%Border%IsDownPMC)) then
               segmento%IsPMC =.true.
            endif
         end select
         if (segmento%isPMC) then
            write (buff,'(a,i7,a)') 'wir1_WARNING: PMC endings in WIREs are UNTESTED'
            if ((segmento%k >= ZI).and.(segmento%k <= ZE)) call WarnErrReport(buff)
         endif
      end do

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !find nodes at the boundaries
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !ojo que lo del shielding siguiente no es valido si hay branching. habria que shieldear canda segmento del junction !\E7\E7
      do i1=1,HWires%NumChargeNodes
         nodo=>HWires%Chargenode(i1)
         if (nodo%i==SINPML_FULLSIZE(iHx)%XI) then
            if (sgg%Border%IsBackPML.or.sgg%Border%IsBackMur) then
               thereAreMurConditions=.true.
               nodo%IsMur =.true.
               nodo%IsBackDownLeftMur =.true.
               nodo%IsFrontUpRightMur =.false.
               dummy => nodo%CurrentPlus_1
               !adjust constants
               HWires%Chargenode(i1)%NodeInside =>  HWires%Chargenode(i1)%CurrentPlus_1%ChargePlus
               !!dummy%IsShielded                                                  =.true.
               !!dummy%ChargePlus%CurrentPlus_1%IsShielded                         =.true.
               !!dummy%ChargePlus%CurrentPlus_1%ChargePlus%CurrentPlus_1%IsShielded=.true.
               nodo%cteMur=(sqrt(InvMu (dummy%indexmed)*InvEps(dummy%indexmed))  *sgg%dt- &
               dummy%delta )/ &
               (sqrt(InvMu (dummy%indexmed)*InvEps(dummy%indexmed))  *sgg%dt+ &
               dummy%delta )

               !
               write (buff,'(a,3i7,a)')  'wir1_WARNING: Node ',nodo%I,nodo%J,nodo%k,'with BACK Mur conditions'
               if ((nodo%k >=  ZI).and.(nodo%k <= ZE)) call WarnErrReport(buff)
            elseif (sgg%Border%IsBackPeriodic) then
               nodo%IsPeriodic =.true.
               write (buff,'(a,3i7,a)')  'wir1_WARNING: Node ',nodo%I,nodo%J,nodo%k,'with BACK Periodic conditions'
               if ((nodo%k >=  ZI).and.(nodo%k <= ZE)) call WarnErrReport(buff)
            elseif (sgg%Border%IsBackPEC) then
               nodo%IsPec =.true.
               write (buff,'(a,3i7,a)')  'wir1_WARNING: Node ',nodo%I,nodo%J,nodo%k,'with BACK PEC boundary conditions'
               if ((nodo%k >=  ZI).and.(nodo%k <= ZE)) call WarnErrReport(buff)
            endif
         endif
         if (nodo%i==SINPML_FULLSIZe(iHx)%XE) then
            if (sgg%Border%IsFrontPML.or.sgg%Border%IsFrontMur) then
               thereAreMurConditions=.true.
               nodo%IsMur =.true.
               nodo%IsBackDownLeftMur =.false.
               nodo%IsFrontUpRightMur =.true.
               dummy => nodo%CurrentMInus_1
               !adjust constantas shielding 3 levels of segments (no junctions permitted) !must be tested!!!!
               HWires%Chargenode(i1)%NodeInside =>  HWires%Chargenode(i1)%CurrentMinus_1%ChargeMinus
               !!dummy%IsShielded                                                      =.true.  !One current back
               !!dummy%ChargeMinus%CurrentMinus_1%IsShielded                           =.true.  !Two currents back
               !!dummy%ChargeMinus%CurrentMinus_1%ChargeMinus%CurrentMinus_1%IsShielded=.true.  !Three current back
               nodo%cteMur=(sqrt(InvMu (dummy%indexmed)*InvEps(dummy%indexmed))*sgg%dt- &
               dummy%delta)/ &
               (sqrt(InvMu (dummy%indexmed)*InvEps(dummy%indexmed))*sgg%dt+ &
               dummy%delta)
               !
               write (buff,'(a,3i7,a)')  'wir1_WARNING: Node ',nodo%I,nodo%J,nodo%k,'with FRONT Mur conditions'
               if ((nodo%k >=  ZI).and.(nodo%k <= ZE)) call WarnErrReport(buff)
            elseif (sgg%Border%IsFrontPeriodic) then
               nodo%IsPeriodic =.true.
               write (buff,'(a,3i7,a)')  'wir1_WARNING: Node ',nodo%I,nodo%J,nodo%k,'with Front Periodic conditions'
               if ((nodo%k >=  ZI).and.(nodo%k <= ZE)) call WarnErrReport(buff)
            elseif (sgg%Border%IsFrontPEC) then
               nodo%IsPec =.true.
               write (buff,'(a,3i7,a)')  'wir1_WARNING: Node ',nodo%I,nodo%J,nodo%k,'with FRONT PEC boundary conditions'
               if ((nodo%k >=  ZI).and.(nodo%k <= ZE)) call WarnErrReport(buff)
            endif
         endif
         !
         if (nodo%j==SINPML_FULLSIZe(iHy)%YI) then
            if (sgg%Border%IsLeftPML.or.sgg%Border%IsLeftMur) then
               thereAreMurConditions=.true.
               nodo%IsMur =.true.
               nodo%IsBackDownLeftMur =.true.
               nodo%IsFrontUpRightMur =.false.
               dummy => nodo%CurrentPlus_1
               !adjust constants
               HWires%Chargenode(i1)%NodeInside =>  HWires%Chargenode(i1)%CurrentPlus_1%ChargePlus
               !!dummy%IsShielded                                                  =.true.
               !!dummy%ChargePlus%CurrentPlus_1%IsShielded                         =.true.
               !!dummy%ChargePlus%CurrentPlus_1%ChargePlus%CurrentPlus_1%IsShielded=.true.
               nodo%cteMur=(sqrt(InvMu (dummy%indexmed)*InvEps(dummy%indexmed))  *sgg%dt- &
               dummy%delta )/ &
               (sqrt(InvMu (dummy%indexmed)*InvEps(dummy%indexmed))  *sgg%dt+ &
               dummy%delta )

               !
               write (buff,'(a,3i7,a)')  'wir1_WARNING: Node ',nodo%I,nodo%J,nodo%k,'with LEFT Mur conditions'
               if ((nodo%k >=  ZI).and.(nodo%k <= ZE)) call WarnErrReport(buff)
            elseif (sgg%Border%IsLeftPeriodic) then
               nodo%IsPeriodic =.true.
               write (buff,'(a,3i7,a)')  'wir1_WARNING: Node ',nodo%I,nodo%J,nodo%k,'with Left Periodic conditions'
               if ((nodo%k >=  ZI).and.(nodo%k <= ZE)) call WarnErrReport(buff)
            elseif (sgg%Border%IsLeftPEC) then
               nodo%IsPec =.true.
               write (buff,'(a,3i7,a)')  'wir1_WARNING: Node ',nodo%I,nodo%J,nodo%k,'with LEFT PEC boundary conditions'
               if ((nodo%k >=  ZI).and.(nodo%k <= ZE)) call WarnErrReport(buff)
            endif
         endif
         if (nodo%j==SINPML_FULLSIZe(iHy)%YE) then
            if (sgg%Border%IsRightPML.or.sgg%Border%IsRightMur) then
               thereAreMurConditions=.true.
               nodo%IsMur =.true.
               nodo%IsBackDownLeftMur =.false.
               nodo%IsFrontUpRightMur =.true.
               dummy => nodo%CurrentMInus_1
               !adjust constantas shielding 3 levels of segments (no junctions permitted) !must be tested!!!!
               HWires%Chargenode(i1)%NodeInside =>  HWires%Chargenode(i1)%CurrentMinus_1%ChargeMinus
               !!dummy%IsShielded                                                      =.true.  !One current back
               !!dummy%ChargeMinus%CurrentMinus_1%IsShielded                           =.true.  !Two currents back
               !!dummy%ChargeMinus%CurrentMinus_1%ChargeMinus%CurrentMinus_1%IsShielded=.true.  !Three current back
               nodo%cteMur=(sqrt(InvMu (dummy%indexmed)*InvEps(dummy%indexmed))*sgg%dt- &
               dummy%delta)/ &
               (sqrt(InvMu (dummy%indexmed)*InvEps(dummy%indexmed))*sgg%dt+ &
               dummy%delta)
               !
               write (buff,'(a,3i7,a)')  'wir1_WARNING: Node ',nodo%I,nodo%J,nodo%k,'with RIGHT Mur conditions'
               if ((nodo%k >=  ZI).and.(nodo%k <= ZE)) call WarnErrReport(buff)
            elseif (sgg%Border%IsRightPeriodic) then
               nodo%IsPeriodic =.true.
               write (buff,'(a,3i7,a)')  'wir1_WARNING: Node ',nodo%I,nodo%J,nodo%k,'with Right Periodic conditions'
               if ((nodo%k >=  ZI).and.(nodo%k <= ZE)) call WarnErrReport(buff)
            elseif (sgg%Border%IsRightPEC) then
               nodo%IsPec =.true.
               write (buff,'(a,3i7,a)')  'wir1_WARNING: Node ',nodo%I,nodo%J,nodo%k,'with RIGHT PEC boundary conditions'
               if ((nodo%k >=  ZI).and.(nodo%k <= ZE)) call WarnErrReport(buff)
            endif
         endif
         !
         if (nodo%k==SINPML_FULLSIZe(iHz)%ZI) then
            if (sgg%Border%IsDownPML.or.sgg%Border%IsDownMur) then
               thereAreMurConditions=.true.
               nodo%IsMur =.true.
               nodo%IsBackDownLeftMur =.true.
               nodo%IsFrontUpRightMur =.false.
               dummy => nodo%CurrentPlus_1
               !adjust constants
               HWires%Chargenode(i1)%NodeInside =>  HWires%Chargenode(i1)%CurrentPlus_1%ChargePlus
               !!dummy%IsShielded                                                  =.true.
               !!dummy%ChargePlus%CurrentPlus_1%IsShielded                         =.true.
               !!dummy%ChargePlus%CurrentPlus_1%ChargePlus%CurrentPlus_1%IsShielded=.true.
               nodo%cteMur=(sqrt(InvMu (dummy%indexmed)*InvEps(dummy%indexmed))  *sgg%dt- &
               dummy%delta )/ &
               (sqrt(InvMu (dummy%indexmed)*InvEps(dummy%indexmed))  *sgg%dt+ &
               dummy%delta )

               !
               write (buff,'(a,3i7,a)')  'wir1_WARNING: Node ',nodo%I,nodo%J,nodo%k,'with DOWN Mur conditions'
               if ((nodo%k >=  ZI).and.(nodo%k <= ZE)) call WarnErrReport(buff)
            elseif (sgg%Border%IsDownPeriodic) then
               nodo%IsPeriodic =.true.
               write (buff,'(a,3i7,a)')  'wir1_WARNING: Node ',nodo%I,nodo%J,nodo%k,'with Down Periodic conditions'
               if ((nodo%k >=  ZI).and.(nodo%k <= ZE)) call WarnErrReport(buff)
            elseif (sgg%Border%IsDownPEC) then
               nodo%IsPec =.true.
               write (buff,'(a,3i7,a)')  'wir1_WARNING: Node ',nodo%I,nodo%J,nodo%k,'with DOWN PEC boundary conditions'
               if ((nodo%k >=  ZI).and.(nodo%k <= ZE)) call WarnErrReport(buff)
            endif
         endif
         if (nodo%k==SINPML_FULLSIZe(iHz)%ZE) then
            if (sgg%Border%IsUpPML.or.sgg%Border%IsUpMur) then
               thereAreMurConditions=.true.
               nodo%IsMur =.true.
               nodo%IsBackDownLeftMur =.false.
               nodo%IsFrontUpRightMur =.true.
               dummy => nodo%CurrentMInus_1
               !adjust constantas shielding 3 levels of segments (no junctions permitted) !must be tested!!!!
               HWires%Chargenode(i1)%NodeInside =>  HWires%Chargenode(i1)%CurrentMinus_1%ChargeMinus
               !!dummy%IsShielded                                                      =.true.  !One current back
               !!dummy%ChargeMinus%CurrentMinus_1%IsShielded                           =.true.  !Two currents back
               !!dummy%ChargeMinus%CurrentMinus_1%ChargeMinus%CurrentMinus_1%IsShielded=.true.  !Three current back
               nodo%cteMur=(sqrt(InvMu (dummy%indexmed)*InvEps(dummy%indexmed))*sgg%dt- &
               dummy%delta)/ &
               (sqrt(InvMu (dummy%indexmed)*InvEps(dummy%indexmed))*sgg%dt+ &
               dummy%delta)
               !
               write (buff,'(a,3i7,a)')  'wir1_WARNING: Node ',nodo%I,nodo%J,nodo%k,'with UP Mur conditions'
               if ((nodo%k >=  ZI).and.(nodo%k <= ZE)) call WarnErrReport(buff)
            elseif (sgg%Border%IsUpPeriodic) then
               nodo%IsPeriodic =.true.
               write (buff,'(a,3i7,a)')  'wir1_WARNING: Node ',nodo%I,nodo%J,nodo%k,'with Up Periodic conditions'
               if ((nodo%k >=  ZI).and.(nodo%k <= ZE)) call WarnErrReport(buff)
            elseif (sgg%Border%IsUpPEC) then
               nodo%IsPec =.true.
               write (buff,'(a,3i7,a)')  'wir1_WARNING: Node ',nodo%I,nodo%J,nodo%k,'with UP PEC boundary conditions'
               if ((nodo%k >=  ZI).and.(nodo%k <= ZE)) call WarnErrReport(buff)
            endif
         endif
         if ((nodo%isMur).and.(nodo%NumCurrentPlus+nodo%NumCurrentMinus >1)) then
            write (buff,'(a,3i7,a)')  'wir1_WARNING: Node ',nodo%I,nodo%J,nodo%k,'is both a non-open and a Mur'
            if ((nodo%k >=  ZI).and.(nodo%k<=ZE)) call WarnErrReport(buff)
         endif
      end do


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Constantes de actualizacion de nodos
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do i1=1,HWires%NumChargeNodes
         nodo => HWires%ChargeNode(i1)
         !asign constants taking junctions into account
         deltadummy=0.0_RKIND_wires
         indexnode = 0
         do j1=1,nodo%NumCurrentPlus + nodo%NumCurrentMinus
            if (j1 <= nodo%NumCurrentPlus) then
               select case (j1)
                case (1)
                  segmento=>nodo%CurrentPlus_1
                case (2)
                  segmento=>nodo%CurrentPlus_2
                case (3)
                  segmento=>nodo%CurrentPlus_3
                case (4)
                  segmento=>nodo%CurrentPlus_4
                case (5)
                  segmento=>nodo%CurrentPlus_5
                case (6)
                  segmento=>nodo%CurrentPlus_6
                case (7)
                  segmento=>nodo%CurrentPlus_7
                case (8)
                  segmento=>nodo%CurrentPlus_8
                case (9)
                  segmento=>nodo%CurrentPlus_9
               end select
            else
               select case (j1-nodo%NumCurrentPlus)
                case (1)
                  segmento=>nodo%CurrentMinus_1
                case (2)
                  segmento=>nodo%CurrentMinus_2
                case (3)
                  segmento=>nodo%CurrentMinus_3
                case (4)
                  segmento=>nodo%CurrentMinus_4
                case (5)
                  segmento=>nodo%CurrentMinus_5
                case (6)
                  segmento=>nodo%CurrentMinus_6
                case (7)
                  segmento=>nodo%CurrentMinus_7
                case (8)
                  segmento=>nodo%CurrentMinus_8
                case (9)
                  segmento=>nodo%CurrentMinus_9
               end select
            endif
            deltadummy=deltadummy + segmento%delta/2.0_RKIND_wires   !holland pag. 92
            indexnode = indexnode + segmento%OrigIndex
            !
         end do
         !
         nodo%indexnode = indexnode !just for normalization of numbers and debugging mpi/nonmpi
         nodo%CteProp = 1.0_RKIND_wires
         nodo%CtePlain = sgg%dt/deltadummy
         if (nodo%NumCurrentPlus + nodo%NumCurrentMinus == 1) then
            nodo%CtePlain = sgg%dt/(2.0_RKIND_wires * deltadummy) !correct the ending in case of open terminal
         endif
         NODO%oRIGctePlain =nodo%CtePlain !lo salvo para que al ungroundear retome esta constante
      end do

      !!!!!!!!!!
      !corrige cteplain para PEC lossy nodes
      do i1=1,HWires%NumChargeNodes
        nodo => HWires%ChargeNode(i1)
        if (nodo%IsLossy.or.nodo%isPEC) then
            nodo%CtePlain = 0.0_RKIND_wires
            nodo%cteprop  = 0.0_RKIND_wires
        endif
      end do
   
      !!!!FINAL reporting of unGROUNDING OF PEC/LOSSY NODES

      do i1=1,HWires%NumChargeNodes
         nodo => Hwires%ChargeNode(i1)
         if ((nodo%k >= sgg%Sweep(iHz)%ZI).and. &
         (nodo%k <= sgg%Sweep(iHz)%ZE)) then !check that the node is inside the layout (just for MPI)
            if (nodo%IsPec.or.nodo%IsLossy) then
               if ((nodo%Is_LeftEnd).or.(nodo%Is_RightEnd).or.(nodo%NumCurrentPlus + nodo%NumCurrentMinus < 2)) then !es un nodo terminal
                    continue
               else !no es un nodo terminal LeftEnd o RightEnd o abierto pero Si es un nodo material
                     if (.not.groundwires) then
                         if (nodo%IsPec) write (buff,'(a,i7,3i7,a,2i3,a,a)')  'wir1_BUGGYERROR: NON-terminal node lying on PEC ()', &
                        nodo%indexnode, nodo%i,nodo%j,nodo%k, &
                        ' (',nodo%numcurrentminus,nodo%numcurrentplus,')'
                         if (nodo%IsLossy) write (buff,'(a,i7,3i7,a,2i3,a,a)')  'wir1_BUGGYERROR: NON-terminal node lying on Lossy ()', &
                        nodo%indexnode, nodo%i,nodo%j,nodo%k, &
                        ' (',nodo%numcurrentminus,nodo%numcurrentplus,')'
                        if ((nodo%k >=  ZI).and.(nodo%k <= ZE)) call WarnErrReport(buff,.true.)
                        !!if (nodo%isPEC)   nodo%isPEC = .false.
                        !!if (nodo%isLossy) nodo%isLossy = .false.
                        !!!restaura las constantes
                        !!nodo%CteProp = 1.0_RKIND_wires
                        !!nodo%Cteplain = nodo%OrigCteplain
                     else
                        !lo deja como esta
                        if (nodo%IsPec) write (buff,'(a,i7,3i7,a,2i3,a,a)')  'wir1_INFO: Leaving grounded a  NON-Terminal node lying on PEC (LAUNCHED WITH -groundwires) ', &
                        nodo%indexnode, nodo%i,nodo%j,nodo%k, &
                        ' (',nodo%numcurrentminus,nodo%numcurrentplus,')'
                         if (nodo%IsLossy) write (buff,'(a,i7,3i7,a,2i3,a,a)')  'wir1_INFO: Leaving grounded a  NON-Terminal node lying on Lossy (LAUNCHED WITH -groundwires) ', &
                        nodo%indexnode, nodo%i,nodo%j,nodo%k, &
                        ' (',nodo%numcurrentminus,nodo%numcurrentplus,')'
                        if ((nodo%k >=  ZI).and.(nodo%k <= ZE).and.verbose) call WarnErrReport(buff)
                     endif
               endif
            endif
         endif
      end do

   
      
      !gestion condiciones mur en nodos frontera

      do i1=1,HWires%NumChargeNodes
         nodo => HWires%ChargeNode(i1)

         do j1=1,nodo%NumCurrentPlus + nodo%NumCurrentMinus
            if (j1 <= nodo%NumCurrentPlus) then
               select case (j1)
                case (1)
                  segmento=>nodo%CurrentPlus_1
                case (2)
                  segmento=>nodo%CurrentPlus_2
                case (3)
                  segmento=>nodo%CurrentPlus_3
                case (4)
                  segmento=>nodo%CurrentPlus_4
                case (5)
                  segmento=>nodo%CurrentPlus_5
                case (6)
                  segmento=>nodo%CurrentPlus_6
                case (7)
                  segmento=>nodo%CurrentPlus_7
                case (8)
                  segmento=>nodo%CurrentPlus_8
                case (9)
                  segmento=>nodo%CurrentPlus_9
               end select
            else
               select case (j1-nodo%NumCurrentPlus)
                case (1)
                  segmento=>nodo%CurrentMinus_1
                case (2)
                  segmento=>nodo%CurrentMinus_2
                case (3)
                  segmento=>nodo%CurrentMinus_3
                case (4)
                  segmento=>nodo%CurrentMinus_4
                case (5)
                  segmento=>nodo%CurrentMinus_5
                case (6)
                  segmento=>nodo%CurrentMinus_6
                case (7)
                  segmento=>nodo%CurrentMinus_7
                case (8)
                  segmento=>nodo%CurrentMinus_8
                case (9)
                  segmento=>nodo%CurrentMinus_9
               end select
            endif

            !HOLD 251019
            if     ((segmento%HasAbsorbing_RightEnd).and.(nodo%Is_RightEnd))  then
               thereAreMurConditions=.true.
               nodo%IsMur = .true.
               if (associated( nodo%CurrentPlus_1)) then 
                   nodo%NodeInside =>  nodo%CurrentPlus_1%ChargePlus
                              dummy => nodo%CurrentPlus_1
               elseif (associated(nodo%CurrentMinus_1)) then  
                              dummy => nodo%CurrentMinus_1
                   nodo%NodeInside =>  nodo%CurrentMinus_1%ChargeMinus
               else
                  write (buff,*)  'wir0_BUGGYERROR: Mur1 on wires wrong. '
                  call WarnErrReport(buff,.true.)
               endif
               nodo%cteMur=(sqrt(InvMu (dummy%indexmed)*InvEps(dummy%indexmed))*sgg%dt - dummy%delta)/ &
                     (sqrt(InvMu (dummy%indexmed)*InvEps(dummy%indexmed))*sgg%dt + dummy%delta)

            elseif ((segmento%HasAbsorbing_LeftEnd).and.(nodo%Is_LeftEnd))  then
               thereAreMurConditions=.true.
               nodo%IsMur = .true.
               if (associated( nodo%CurrentPlus_1)) then 
                   nodo%NodeInside =>  nodo%CurrentPlus_1%ChargePlus
                              dummy => nodo%CurrentPlus_1
               elseif (associated(nodo%CurrentMinus_1)) then  
                              dummy => nodo%CurrentMinus_1
                   nodo%NodeInside =>  nodo%CurrentMinus_1%ChargeMinus
               else
                  write (buff,*)  'wir0_BUGGYERROR: Mur1 on wires wrong. '
                  call WarnErrReport(buff,.true.)
               endif
               nodo%cteMur=(sqrt(InvMu (dummy%indexmed)*InvEps(dummy%indexmed))*sgg%dt - dummy%delta)/ &
                     (sqrt(InvMu (dummy%indexmed)*InvEps(dummy%indexmed))*sgg%dt + dummy%delta)
            endif
            !!!
         end do !del barrido de nodos !movido aqui 07/02/14 bug esn.nfde solo en reporte
      end do




      !I create the FranctionPlus and FractionMinus actualization constants needed to update the currents in junctions
      do i1=1,HWires%NumCurrentSegments
         segmento=>HWires%CurrentSegment(i1)
         nodo=>segmento%ChargeMinus
         segmento%FractionMinus=-1.0e30_RKIND_wires !valor absurdo tiene que entrar siempre en algun case
         !Junctions and plain nodes require this correction in case of two distinct radius segments meet
         DenominatorFractionMinusDummy=0.0_RKIND_wires
         deltadummy1=0.0_RKIND_wires
         if (nodo%NumCurrentPlus<=0)  call StopOnError(0,0,'Bug fractionminuss. ')
         do j1=1,nodo%NumCurrentPlus
            select case (j1)
             case (1)
               dummy=>nodo%CurrentPlus_1
             case (2)
               dummy=>nodo%CurrentPlus_2
             case (3)
               dummy=>nodo%CurrentPlus_3
             case (4)
               dummy=>nodo%CurrentPlus_4
             case (5)
               dummy=>nodo%CurrentPlus_5
             case (6)
               dummy=>nodo%CurrentPlus_6
             case (7)
               dummy=>nodo%CurrentPlus_7
             case (8)
               dummy=>nodo%CurrentPlus_8
             case (9)
               dummy=>nodo%CurrentPlus_9
            end select
            deltadummy1=deltadummy1 + dummy%delta
            DenominatorFractionMinusDummy=DenominatorFractionMinusDummy+dummy%delta/(dummy%Lind * InvMu(dummy%indexmed) * InvEPS(dummy%indexmed))
         end do
         do j1=1,nodo%NumCurrentMinus
            select case (j1)
             case (1)
               dummy=>nodo%CurrentMinus_1
             case (2)
               dummy=>nodo%CurrentMinus_2
             case (3)
               dummy=>nodo%CurrentMinus_3
             case (4)
               dummy=>nodo%CurrentMinus_4
             case (5)
               dummy=>nodo%CurrentMinus_5
             case (6)
               dummy=>nodo%CurrentMinus_6
             case (7)
               dummy=>nodo%CurrentMinus_7
             case (8)
               dummy=>nodo%CurrentMinus_8
             case (9)
               dummy=>nodo%CurrentMinus_9
            end select
            deltadummy1=deltadummy1 + dummy%delta
            DenominatorFractionMinusDummy=DenominatorFractionMinusDummy+dummy%delta/(dummy%Lind * InvMu(dummy%indexmed) * InvEPS(dummy%indexmed))
         end do
         segmento%FractionMinus=(deltadummy1/(segmento%Lind * InvMu(segmento%indexmed) * InvEPS(segmento%indexmed)))/(DenominatorFractionMinusDummy)  !Hollond paper on wires'81 (page 91) <--- esta mal. lo de berenger en -wiresflavor new es lo correcto.. yo lo hago ahora 22/07/15
         !
         nodo=>segmento%ChargePlus
         segmento%FractionPlus=-1.0e30_RKIND_wires !valor absurdo tiene que entrar siempre en algun case
         DenominatorFractionPlusDummy=0.0_RKIND_wires
         deltadummy2=0.0_RKIND_wires
         if (nodo%NumCurrentMinus<=0)  call StopOnError(0,0,'Bug fractionplus. ')
         do j1=1,nodo%NumCurrentMinus
            select case (j1)
             case (1)
               dummy=>nodo%CurrentMinus_1
             case (2)
               dummy=>nodo%CurrentMinus_2
             case (3)
               dummy=>nodo%CurrentMinus_3
             case (4)
               dummy=>nodo%CurrentMinus_4
             case (5)
               dummy=>nodo%CurrentMinus_5
             case (6)
               dummy=>nodo%CurrentMinus_6
             case (7)
               dummy=>nodo%CurrentMinus_7
             case (8)
               dummy=>nodo%CurrentMinus_8
             case (9)
               dummy=>nodo%CurrentMinus_9
            end select
            deltadummy2=deltadummy2 + dummy%delta
            DenominatorFractionPlusDummy=DenominatorFractionPlusDummy+dummy%delta/(dummy%Lind * InvMu(dummy%indexmed) * InvEPS(dummy%indexmed))
         end do
         do j1=1,nodo%NumCurrentPlus
            select case (j1)
             case (1)
               dummy=>nodo%CurrentPlus_1
             case (2)
               dummy=>nodo%CurrentPlus_2
             case (3)
               dummy=>nodo%CurrentPlus_3
             case (4)
               dummy=>nodo%CurrentPlus_4
             case (5)
               dummy=>nodo%CurrentPlus_5
             case (6)
               dummy=>nodo%CurrentPlus_6
             case (7)
               dummy=>nodo%CurrentPlus_7
             case (8)
               dummy=>nodo%CurrentPlus_8
             case (9)
               dummy=>nodo%CurrentPlus_9
            end select
            deltadummy2=deltadummy2 + dummy%delta
            DenominatorFractionPlusDummy=DenominatorFractionPlusDummy+dummy%delta/(dummy%Lind * InvMu(dummy%indexmed) * InvEPS(dummy%indexmed))
         end do
         segmento%FractionPlus=(deltadummy2/(segmento%Lind * InvMu(segmento%indexmed) * InvEPS(segmento%indexmed)))/(DenominatorFractionPlusDummy)
         continue
      end do


      !End of the adjusting of constants and sort of nodes


      !detect inverse oriented segments bug OLD segmentos al reves 11/03/15

      do conta=1,HWires%NumCurrentSegments
         !!only for the observation sign to match (not used in this routine)
         if (associated(HWires%CurrentSegment(conta)%chargeplus%currentplus_1)) then
            if       ((HWires%CurrentSegment(conta)%chargeplus%currentplus_1%origindex < HWires%CurrentSegment(conta)%origindex).and.&
            (HWires%CurrentSegment(conta)%chargeplus%currentplus_1%indexmed == HWires%CurrentSegment(conta)%indexmed)) then !tienen que estar en el mismo hilo
               HWires%CurrentSegment(conta)%orientadoalreves =.true. !relies on ORIGINAL orientation
            endif
         endif
         if (associated(HWires%CurrentSegment(conta)%chargeplus%currentplus_2)) then
            if       ((HWires%CurrentSegment(conta)%chargeplus%currentplus_2%origindex < HWires%CurrentSegment(conta)%origindex).and.&
            (HWires%CurrentSegment(conta)%chargeplus%currentplus_2%indexmed == HWires%CurrentSegment(conta)%indexmed)) then !tienen que estar en el mismo hilo
               HWires%CurrentSegment(conta)%orientadoalreves =.true. !relies on ORIGINAL orientation
            endif
         endif
         if (associated(HWires%CurrentSegment(conta)%chargeplus%currentminus_1)) then !puede que sea yo mismo, por eso miro tambien el currentminus2
            if       ((HWires%CurrentSegment(conta)%chargeplus%currentminus_1%origindex < HWires%CurrentSegment(conta)%origindex).and.&
            (HWires%CurrentSegment(conta)%chargeplus%currentminus_1%indexmed == HWires%CurrentSegment(conta)%indexmed)) then !tienen que estar en el mismo hilo
               HWires%CurrentSegment(conta)%orientadoalreves =.true. !relies on ORIGINAL orientation
            endif
         endif
         if (associated(HWires%CurrentSegment(conta)%chargeplus%currentminus_2)) then
            if       ((HWires%CurrentSegment(conta)%chargeplus%currentminus_2%origindex < HWires%CurrentSegment(conta)%origindex).and.&
            (HWires%CurrentSegment(conta)%chargeplus%currentminus_2%indexmed == HWires%CurrentSegment(conta)%indexmed)) then !tienen que estar en el mismo hilo
               HWires%CurrentSegment(conta)%orientadoalreves =.true. !relies on ORIGINAL orientation
            endif
         endif
         !!!!!!!!!
         if (associated(HWires%CurrentSegment(conta)%chargeminus%currentminus_1)) then
            if       ((HWires%CurrentSegment(conta)%chargeminus%currentminus_1%origindex > HWires%CurrentSegment(conta)%origindex).and. &
            (HWires%CurrentSegment(conta)%chargeminus%currentminus_1%indexmed == HWires%CurrentSegment(conta)%indexmed)) then
               HWires%CurrentSegment(conta)%orientadoalreves =.true. !relies on ORIGINAL orientation
            endif
         endif
         if (associated(HWires%CurrentSegment(conta)%chargeminus%currentminus_2)) then
            if       ((HWires%CurrentSegment(conta)%chargeminus%currentminus_2%origindex > HWires%CurrentSegment(conta)%origindex).and. &
            (HWires%CurrentSegment(conta)%chargeminus%currentminus_2%indexmed == HWires%CurrentSegment(conta)%indexmed)) then
               HWires%CurrentSegment(conta)%orientadoalreves =.true. !relies on ORIGINAL orientation
            endif
         endif
         if (associated(HWires%CurrentSegment(conta)%chargeminus%currentplus_1)) then  !puede que sea yo mismo, por eso miro tambien el currentplus2
            if       ((HWires%CurrentSegment(conta)%chargeminus%currentplus_1%origindex > HWires%CurrentSegment(conta)%origindex).and. &
            (HWires%CurrentSegment(conta)%chargeminus%currentplus_1%indexmed == HWires%CurrentSegment(conta)%indexmed)) then
               HWires%CurrentSegment(conta)%orientadoalreves =.true. !relies on ORIGINAL orientation
            endif
         endif
         if (associated(HWires%CurrentSegment(conta)%chargeminus%currentplus_2)) then
            if       ((HWires%CurrentSegment(conta)%chargeminus%currentplus_2%origindex > HWires%CurrentSegment(conta)%origindex).and. &
            (HWires%CurrentSegment(conta)%chargeminus%currentplus_2%indexmed == HWires%CurrentSegment(conta)%indexmed)) then
               HWires%CurrentSegment(conta)%orientadoalreves =.true. !relies on ORIGINAL orientation
            endif
         endif
      end do


      !Current sources (only permitted in Plain nodes -no hanging-)
      !Read the file with the time evolution of the source
      !I find ambiguity in the geometrical position (-untested with OLD-) which is given in .nfde
      !at at segment index (i,j,k). So  I assume that the current source
      !is located at the node wich is the minus direction of the segment. That is I detect if the the CurrentPlus_1
      !is there some current source

      do i1=1,HWires%NumChargeNodes
         If (associated(HWires%ChargeNode(i1)%CurrentPlus_1)) then
            dummy=>HWires%ChargeNode(i1)%CurrentPlus_1
            DO k1=1,dummy%Tipowire%numcurrentsources
               if (dummy%TipoWire%IsourceExists) then
                  if ((dummy%i == dummy%Tipowire%Isource(k1)%I).and. &
                  (dummy%j == dummy%Tipowire%Isource(k1)%J).and. &
                  (dummy%k == dummy%Tipowire%Isource(k1)%K)) then
                     if (HWires%ChargeNode(i1)%NumCurrentPlus+HWires%ChargeNode(i1)%NumCurrentMinus>2) then
                        write (buff,'(a,i7,a)') 'wir1_ERROR: Current source at ',HWIREs%ChargeNode(i1)%indexnode,' in junctions forbidden'
                        call WarnErrReport(buff,.true.)
                     endif
                     HWires%ChargeNode(i1)%Isource=>dummy%TipoWire%Isource(k1)
                     HWires%ChargeNode(i1)%HasIsource=.true.
                     thereareIsources=.true.
                     write (buff,'(a,i7,3i7)')  'wir1_INFO: Current source at node ', &
                     nodo%indexnode, nodo%i,nodo%j,nodo%k
                     if ((HWires%ChargeNode(i1)%k > ZI).and.(HWires%ChargeNode(i1)%k <= ZE).and.verbose) call WarnErrReport(buff)
                     !
                     if (HWires%ChargeNode(i1)%Isource%Fichero%DeltaSamples > sgg%dt) then
                        write (buff,'(a,e15.4e3)') 'wir1_WARNING: '//trim(adjustl(HWires%ChargeNode(i1)%Isource%Fichero%Name))// &
                        ' undersampled by a factor ',HWires%ChargeNode(i1)%Isource%Fichero%DeltaSamples/sgg%dt
                        if ((HWires%ChargeNode(i1)%k > ZI).and.(HWires%ChargeNode(i1)%k <= ZE)) call WarnErrReport(buff)
                     endif
                  endif
               endif
            end do
         endif
      end do


      !End of geometrical PREPROCESSING
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call resume_casuistics
    
      !!!!!!!!find crank-nicolson coefficients
       
      if (wirecrank) then
          call init_wirecrank
      endif


!!!!Se arregla lo del resuming permit scaling 221118 pero realmente no se por que. deberia ser lo mismo. Es redundante, pero....
       eps000=eps0;mu000=mu0; !niapa para evitar lo del rkind 020719
       call calc_wirehollandconstants(sgg,G2,fieldtotl,wiresflavor,mu000,eps000,simu_devia)
       
#ifdef CompileWithThickWires
!!!init extra thick stuff !hay que hacerlo aqui al final con toda la info
!!!ojo si permittivity scaling y si mpi.... habra que rehacer
      do conta=1,HWires%NumCurrentSegments
            if (wirethickness>1) &
                call init_thick(sgg,eps000,mu000,Ex,Ey,Ez,Hx,Hy,Hz,Idxe,Idye,Idze,Idxh,Idyh,Idzh,HWires%CurrentSegment(conta),wirethickness)
      end do
#endif      
      !!
      return

   contains
 
   
 !!!!!!!!!
   
   subroutine deembed_peclossyconformal_segments(sggmiE)
    
       integer (KIND=INTEGERSIZEOFMEDIAMATRICES) :: sggMiE
      
   !primero los conformal 130220 %Is%split_and_useless
       if ((sgg%Med(sggmiE)%Is%split_and_useless).and. &
                .not.(IsEnd_norLeft_norRight.or.Is_LeftEnd.or.Is_RightEnd)) then   !NO NO NO ES UN TERMINAL
                   call deembed_segment 
                   write (buff,'(a,6i9)') 'wir0_WARNING: YES de-embedding a NON-TERMINAL conformal split_and_useless WIRE segment: ', sggmiE, &
                            HWires%CurrentSegment(conta)%origIndex,HWires%CurrentSegment(conta)%i, &
                            HWires%CurrentSegment(conta)%j,HWires%CurrentSegment(conta)%k,HWires%CurrentSegment(conta)%tipofield
                   if ((k1 >= ZI).and.(k1 <= ZE)) call WarnErrReport(buff)
        elseif ((sgg%Med(sggmiE)%Is%split_and_useless).and. &
                (IsEnd_norLeft_norRight.or.Is_LeftEnd.or.Is_RightEnd)) then   !SI SI SI ES UN TERMINAL
                   call deembed_segment 
                   write (buff,'(a,6i9)') 'wir0_ERROR: YES-TERMINAL WIRE SEGMENT IN A CONFORMAL split_and_useless SURFACE (): ', sggmiE, &
                           HWires%CurrentSegment(conta)%origIndex,HWires%CurrentSegment(conta)%i, &
                          HWires%CurrentSegment(conta)%j,HWires%CurrentSegment(conta)%k,HWires%CurrentSegment(conta)%tipofield
                   if ((k1 >= ZI).and.(k1 <= ZE)) call WarnErrReport(buff,.true.)
        elseif ((sgg%Med(sggmiE)%Is%already_YEEadvanced_byconformal).and. &  !!!!!!!!!!!!already_YEEadvanced_byconformal
                .not.(IsEnd_norLeft_norRight.or.Is_LeftEnd.or.Is_RightEnd)) then   !NO NO NO ES UN TERMINAL                    
                if (.not.fieldtotl) then
                         HWires%CurrentSegment(conta)%cte5 = sgg%dt /eps0 /(HWires%CurrentSegment(conta)%deltaTransv1*HWires%CurrentSegment(conta)%deltaTransv2)
                else
                         HWires%CurrentSegment(conta)%cte5 = 0.0_RKIND_wires
                endif
                   write (buff,'(a,6i9)') 'wir0_WARNING: NO de-embedding a NON-TERMINAL conformal already_YEEadvanced_byconformal  WIRE segment: ', sggmiE, &
                                     HWires%CurrentSegment(conta)%origIndex,HWires%CurrentSegment(conta)%i, &
                                     HWires%CurrentSegment(conta)%j,HWires%CurrentSegment(conta)%k,HWires%CurrentSegment(conta)%tipofield
                   if ((k1 >= ZI).and.(k1 <= ZE)) call WarnErrReport(buff)
        elseif ((sgg%Med(sggmiE)%Is%already_YEEadvanced_byconformal).and. &
                (IsEnd_norLeft_norRight.or.Is_LeftEnd.or.Is_RightEnd)) then  !SI SI SI ES UN TERMINAL                    
                if (.not.fieldtotl) then
                         HWires%CurrentSegment(conta)%cte5 = sgg%dt /eps0 /(HWires%CurrentSegment(conta)%deltaTransv1*HWires%CurrentSegment(conta)%deltaTransv2)
                else
                         HWires%CurrentSegment(conta)%cte5=0.0_RKIND_wires
                endif
                   write (buff,'(a,6i9)') 'wir0_WARNING: NO de-embedding a YES-TERMINAL conformal already_YEEadvanced_byconformal  WIRE segment: ', sggmiE, &
                                  HWires%CurrentSegment(conta)%origIndex,HWires%CurrentSegment(conta)%i, &
                                  HWires%CurrentSegment(conta)%j,HWires%CurrentSegment(conta)%k,HWires%CurrentSegment(conta)%tipofield
                   if ((k1 >= ZI).and.(k1 <= ZE)) call WarnErrReport(buff)
        elseif ((sggmiE == 0).or.(sgg%med(sggmiE)%is%pec).or. &
            (abs(sgg%Med(sggmiE)%sigma) >= 1.0e-15_RKIND_wires).or.(abs(sgg%Med(sggmiE)%sigmam) >= 1.0e-15_RKIND_wires).or. &
            sgg%Med(sggmiE)%Is%Lossy ) then
                  call deembed_segment
            !
                  IF (.not.(IsEnd_norLeft_norRight.or.Is_LeftEnd.or.Is_RightEnd)) then !NO NO NO ES UN TERMINAL
                   if ((sggmiE == 0).or.(sgg%med(sggmiE)%is%pec)) then
                       write (buff,'(a,6i9)')        'wir0_WARNING: YES De-embedding a NON-TERMINAL struct segment from PEC ', sggmiE, &
                                  HWires%CurrentSegment(conta)%origIndex,HWires%CurrentSegment(conta)%i, &
                                  HWires%CurrentSegment(conta)%j,HWires%CurrentSegment(conta)%k,HWires%CurrentSegment(conta)%tipofield
                   else
                       write (buff,'(a,6i9)')        'wir0_WARNING: YES De-embedding a NON-TERMINAL struct segment from Lossy ', sggmiE, &
                                  HWires%CurrentSegment(conta)%origIndex,HWires%CurrentSegment(conta)%i, &
                                  HWires%CurrentSegment(conta)%j,HWires%CurrentSegment(conta)%k,HWires%CurrentSegment(conta)%tipofield
                   endif
                   
                  else  !SI SI SI ES UN TERMINAL
                   if ((sggmiE == 0).or.(sgg%med(sggmiE)%is%pec)) then
                       call deembed_segment !ojoo esto no estaba a 290323 y lo he aniadido porque parece que tiene sentido deembed si lo dice
                       write (buff,'(a,6i9)')  'wir0_SEVEREWARNING: YES de-embedding a YES-TERMINAL struct SEGMENT from PEC  ', sggmiE, &
                                  HWires%CurrentSegment(conta)%origIndex,HWires%CurrentSegment(conta)%i, &
                                  HWires%CurrentSegment(conta)%j,HWires%CurrentSegment(conta)%k,HWires%CurrentSegment(conta)%tipofield
                    else
                       write (buff,'(a,6i9)')  'wir0_SEVEREWARNING: YES de-embedding a YES-TERMINAL struct SEGMENT from Lossy  ', sggmiE, &
                                  HWires%CurrentSegment(conta)%origIndex,HWires%CurrentSegment(conta)%i, &
                                  HWires%CurrentSegment(conta)%j,HWires%CurrentSegment(conta)%k,HWires%CurrentSegment(conta)%tipofield
                     endif
                  endif  
                  if ((k1 >= ZI).and.(k1 <= ZE)) call WarnErrReport(buff)
!!!!  luego normales
          else
            if ((sggmiE /= 1).and.(.not.sgg%Med(sggmiE)%Is%ThinWire)) then
                write (buff,'(a,6i9)')'wir0_WARNING: NO de-embedding a terminal/non-terminal segment in Lossless medium ', sggmiE, &
                      HWires%CurrentSegment(conta)%origIndex,HWires%CurrentSegment(conta)%i, &
                      HWires%CurrentSegment(conta)%j,HWires%CurrentSegment(conta)%k,HWires%CurrentSegment(conta)%tipofield
                      call WarnErrReport(buff)
            endif
          endif
   

          
   end subroutine deembed_peclossyconformal_segments
   !!!
   
   
   subroutine deembed_segment
   
                !AJUSTAR LA CONSTANTE !no se precisa realmente porque luego se iran su campos a cero (ver nota antes) 180220
                    
                if (.not.fieldtotl) then
                         HWires%CurrentSegment(conta)%cte5=sgg%dt /eps0 /(HWires%CurrentSegment(conta)%deltaTransv1*HWires%CurrentSegment(conta)%deltaTransv2)
                else
                         HWires%CurrentSegment(conta)%cte5=0.0_RKIND_wires
                endif
                    
                HWires%CurrentSegment(conta)%Efield_wire2main => HWires%null_field  ! YES de-embedding
                HWires%CurrentSegment(conta)%Efield_main2wire => HWires%null_field  ! YES de-embedding
            
end subroutine deembed_segment
   
   subroutine detect_peclossyconformal_nodes
      logical embed
      
      !detect PEC and lossy nodes
      do i1=1,HWires%NumChargeNodes
         nodo =>  HWires%ChargeNode(i1)
         i = nodo%i
         j = nodo%j
         k = nodo%k

         if ((i >  SINPML_fullsize(iEx)%XI).and. &
         (i <= SINPML_fullsize(iEx)%XE).and. &
         (j >  SINPML_fullsize(iEy)%YI).and. &
         (j <= SINPML_fullsize(iEy)%YE).and. &
         (k >  SINPML_fullsize(iEz)%ZI).and. &
         (k <= SINPML_fullsize(iEz)%ZE)) then



            if ((k <= sgg%alloc(iEZ)%ZE).and.(k >= sgg%alloc(iEZ)%ZI)) then
               kmenos1= k-1
               if (k-1 <  sgg%alloc(iEz)%ZI) kmenos1=k
               !
               Nodo%IsLossy   = &
               sgg%Med(sggMino(i  ,j,k         ))%Is%Lossy .or. &
               sgg%Med(sggMiEx(i  ,j,k         ))%Is%Lossy .or. &
               sgg%Med(sggMiEy(i,j  ,k         ))%Is%Lossy .or. &
               sgg%Med(sggMiEz(i,j,k           ))%Is%Lossy .or. &
               sgg%Med(sggMiEx(i-1,j,k         ))%Is%Lossy .or. &
               sgg%Med(sggMiEy(i,j-1,k         ))%Is%Lossy .or. &
               sgg%Med(sggMiEz(i,j,kmenos1     ))%Is%Lossy

               
               sigt =MAX(max(max(sgg%Med(sggMiEx(i  ,j  ,k      ))%sigma , &
                                   sgg%Med(sggMiEy(i  ,j  ,k      ))%sigma), &
                                   MAX(sgg%Med(sggMiEz(i  ,j  ,k      ))%sigma,  &
                                   sgg%Med(sggMiEx(i-1,j  ,k      ))%sigma)),&
                                   MAX(    sgg%Med(sggMiEy(i  ,j-1,k      ))%sigma,  &
                                   sgg%Med(sggMiEz(i  ,j  ,kmenos1))%sigma) )
               nodo%islossy = nodo%islossy .or. (abs(sigt) > 1.0e-19_RKIND_wires)

            
               nodo%ispec = &
               ( sggMiNo(i  ,j,k) ==0).or.(sgg%med( sggMiNo(i  ,j,k) )%is%pec) .or. &
               ( sggMiEx(i  ,j,k) ==0).or.(sgg%med( sggMiEx(i  ,j,k) )%is%pec) .or. &
               ( sggMiEy(i,j  ,k) ==0).or.(sgg%med( sggMiEy(i,j  ,k) )%is%pec) .or. &
               ( sggMiEz(i,j,k  ) ==0).or.(sgg%med( sggMiEz(i,j,k  ) )%is%pec) .or. &
               ( sggMiEx(i-1,j,k) ==0).or.(sgg%med( sggMiEx(i-1,j,k) )%is%pec) .or. &
               ( sggMiEy(i,j-1,k) ==0).or.(sgg%med( sggMiEy(i,j-1,k) )%is%pec) .or. &
               ( sggMiEz(i,j,kmenos1) ==0).or.(sgg%med( sggMiEz(i,j,kmenos1) )%is%pec)
               
               if ((.not.(nodo%Is_RightEnd.or.nodo%Is_LeftEnd)).and.(nodo%ispec.or.nodo%islossy)) then
                   if (nodo%ispec) write (buff,*)  'wir1_WARNING: Un-grounding NON-TERMINAL node from PEC ', i,j,k
                   if (nodo%islossy) write (buff,*)  'wir1_WARNING: Un-grounding NON-TERMINAL node from Lossy ', i,j,k
                   if ((k >  ZI).and.(k <= ZE)) call WarnErrReport(buff)
                   nodo%ispec=.false.
                   nodo%islossy=.false.
               endif
               
               
!pedazo de niapa para poner los nodos conformal a voltage nulo y sus segmentos conformal tambien si son already_YEEadvanced_byconformal (old notouch o no_touch) y que Dios reparta suerte 140220
               if ((nodo%Is_RightEnd.or.nodo%Is_LeftEnd)) then !!!los busy_nodes se han puesto a pec cuando en realidad estan unidos a already_YEEadvanced_byconformal .and.(.not.(nodo%ispec.or.nodo%islossy))) then   !solo si es un extremo y no estaba ya puesto a pec            
                      medio1  =sggMiEx(i  ,j  ,k  )    
                      medio1m =sggMiEx(i-1,j  ,k  )
                      medio2  =sggMiEy(i  ,j  ,k  )
                      medio2m =sggMiEy(i  ,j-1,k  )
                      medio3  =sggMiEz(i  ,j  ,k  )    
                      medio3m =sggMiEz(i  ,j  ,kmenos1)
                      
                     !
                      if (sgg%med(medio1 )%is%split_and_useless .or. sgg%med(medio2 )%is%split_and_useless .or. sgg%med(medio3 )%is%split_and_useless.or. &
                          sgg%med(medio1m)%is%split_and_useless .or. sgg%med(medio2m)%is%split_and_useless .or. sgg%med(medio3m)%is%split_and_useless) then
                           write (buff,*)  'wir1_BUGGYERROR: Conformal node CONNECTED TO at least one split_and_useless conformal edge that cannot be safely set to PEC () ', i,j,k
                           if ((k >  ZI).and.(k <= ZE)) call WarnErrReport(buff,.true.)
                      endif
                     !
                      if (sgg%med(medio1 )%is%already_YEEadvanced_byconformal .or. sgg%med(medio2 )%is%already_YEEadvanced_byconformal .or. sgg%med(medio3 )%is%already_YEEadvanced_byconformal.or. &
                          sgg%med(medio1m)%is%already_YEEadvanced_byconformal .or. sgg%med(medio2m)%is%already_YEEadvanced_byconformal .or. sgg%med(medio3m)%is%already_YEEadvanced_byconformal) then
			               nodo%ispec=.true. !luego se pondra nodo%cteplain = 0 para todos los pec
                           write (buff,*)  'wir1_WARNING: Conformal already_YEEadvanced_byconformal node changed to PEC grounded node at ', i,j,k
                           if ((k >  ZI).and.(k <= ZE)) call WarnErrReport(buff)
                           !ademas anotar y poner a cero el efield correspondiente already_YEEadvanced_byconformal overrideando el conformal_advance_E() que ya se ha hecho antes a partir de esta version   
                           if (sgg%med(medio1)%is%already_YEEadvanced_byconformal) then
                               call check_embed(embed,iEx,i,j,k)
                               if (.not.embed) then
                                  ! sggmiEx(i,j,k)=0;!ojoo quitar luego solo para visualiz
                                   nodo%already_YEEadvanced_byconformal_changedtoPECfield1 => Ex(i,j,k)
                                   write (buff,*)  'wir1_WARNING: Conformal already_YEEadvanced_byconformal segment atached to PEC grounded node, also changed to PEC line at Ex ', i,j,k
                                   if ((k >  ZI).and.(k <= ZE)) call WarnErrReport(buff)
                               else
                                   write (buff,*)  'wir1_WARNING: Conformal already_YEEadvanced_byconformal segment atached to PEC grounded node, cannot be also changed to PEC line for embedding a wire at Ex ', i,j,k
                                   if ((k >  ZI).and.(k <= ZE)) call WarnErrReport(buff) 
                                endif
                           endif  
                           if (sgg%med(medio1m)%is%already_YEEadvanced_byconformal) then
                               call check_embed(embed,iEx,i-1,j,k)
                               if (.not.embed) then
                               !    sggmiEx(i-1,j,k)=0;!ojoo quitar luego solo para visualiz
                                   nodo%already_YEEadvanced_byconformal_changedtoPECfield2 => Ex(i-1,j,k)
                                   write (buff,*)  'wir1_WARNING: Conformal already_YEEadvanced_byconformal segment atached to PEC grounded node, also changed to PEC line at mEx ', i-1,j,k
                                   if ((k >  ZI).and.(k <= ZE)) call WarnErrReport(buff)
                               else
                                   write (buff,*)  'wir1_WARNING: Conformal already_YEEadvanced_byconformal segment atached to PEC grounded node, cannot be also changed to PEC line for embedding a wire at mEx ', i-1,j,k
                                   if ((k >  ZI).and.(k <= ZE)) call WarnErrReport(buff) 
                                endif
                           endif 
                           if (sgg%med(medio2)%is%already_YEEadvanced_byconformal) then
                               call check_embed(embed,iEy,i,j,k)
                               if (.not.embed) then
                            !       sggmiEy(i,j,k)=0;!ojoo quitar luego solo para visualiz
                                   nodo%already_YEEadvanced_byconformal_changedtoPECfield3 => Ey(i,j,k)
                                   write (buff,*)  'wir1_WARNING: Conformal already_YEEadvanced_byconformal segment atached to PEC grounded node, also changed to PEC line at Ey ', i,j,k
                                   if ((k >  ZI).and.(k <= ZE)) call WarnErrReport(buff)
                               else
                                   write (buff,*)  'wir1_WARNING: Conformal already_YEEadvanced_byconformal segment atached to PEC grounded node, cannot be also changed to PEC line for embedding a wire at Ey ', i,j,k
                                   if ((k >  ZI).and.(k <= ZE)) call WarnErrReport(buff) 
                                endif
                           endif 
                           if (sgg%med(medio2m)%is%already_YEEadvanced_byconformal) then
                               call check_embed(embed,iEy,i,j-1,k)
                               if (.not.embed) then
                          !         sggmiEy(i,j-1,k)=0;!ojoo quitar luego solo para visualiz
                                   nodo%already_YEEadvanced_byconformal_changedtoPECfield4 => Ey(i,j-1,k)
                                   write (buff,*)  'wir1_WARNING: Conformal already_YEEadvanced_byconformal segment atached to PEC grounded node, also changed to PEC line at mEy ', i,j-1,k
                                   if ((k >  ZI).and.(k <= ZE)) call WarnErrReport(buff)
                               else
                                   write (buff,*)  'wir1_WARNING: Conformal already_YEEadvanced_byconformal segment atached to PEC grounded node, cannot be also changed to PEC line for embedding a wire at mEy ', i,j-1,k
                                   if ((k >  ZI).and.(k <= ZE)) call WarnErrReport(buff) 
                                endif
                           endif 
                           if (sgg%med(medio3)%is%already_YEEadvanced_byconformal) then
                               call check_embed(embed,iEz,i,j,k)
                               if (.not.embed) then
                          !         sggmiEz(i,j,k)=0;!ojoo quitar luego solo para visualiz
                                   nodo%already_YEEadvanced_byconformal_changedtoPECfield5 => Ez(i,j,k)
                                   write (buff,*)  'wir1_WARNING: Conformal already_YEEadvanced_byconformal segment atached to PEC grounded node, also changed to PEC line at Ez ', i,j,k
                                   if ((k >  ZI).and.(k <= ZE)) call WarnErrReport(buff)
                               else
                                   write (buff,*)  'wir1_WARNING: Conformal already_YEEadvanced_byconformal segment atached to PEC grounded node, cannot be also changed to PEC line for embedding a wire at Ez ', i,j,k
                                   if ((k >  ZI).and.(k <= ZE)) call WarnErrReport(buff) 
                                endif
                           endif
                           if (sgg%med(medio3m)%is%already_YEEadvanced_byconformal) then
                               call check_embed(embed,iEz,i,j,kmenos1)
                               if (.not.embed) then
                          !         sggmiEz(i,j,kmenos1)=0;!ojoo quitar luego solo para visualiz
                                   nodo%already_YEEadvanced_byconformal_changedtoPECfield6 => Ez(i,j,kmenos1)
                                   write (buff,*)  'wir1_WARNING: Conformal already_YEEadvanced_byconformal segment atached to PEC grounded node, also changed to PEC line at mEz ', i,j,kmenos1
                                   if ((k >  ZI).and.(k <= ZE)) call WarnErrReport(buff)
                               else
                                   write (buff,*)  'wir1_WARNING: Conformal already_YEEadvanced_byconformal segment atached to PEC grounded node, cannot be also changed to PEC line for embedding a wire at mEz ', i,j,kmenos1
                                   if ((k >  ZI).and.(k <= ZE)) call WarnErrReport(buff) 
                                endif
                           endif
                           
                      ELSEIF (sgg%med(medio1 )%is%split_and_useless .and. sgg%med(medio2 )%is%split_and_useless .and. sgg%med(medio3 )%is%split_and_useless .and. &
                              sgg%med(medio1m)%is%split_and_useless .and. sgg%med(medio2m)%is%split_and_useless .and. sgg%med(medio3m)%is%split_and_useless) then
                           write (buff,*)  'wir1_ERROR: Conformal split_and_useless node NOT changed (IMPOSIBLE) to PEC grounded node at ', i,j,k
                           if ((k >  ZI).and.(k <= ZE)) call WarnErrReport(buff,.true.)
                           !ojo no tendria arreglo porque aunque se updatee a 0, este campo luego esta partido en dos en conformal y es split_and_useless        
                      endif
                      !esta a PEC/lossy o al punietero aire
                              !!!lo he sacado del if. solo para reporte y que coincida con lo que da estructurado  280323
                      if (nodo%ispec)  then
                          write (buff,*)  'wir1_INFO: (SHOULD BE REDUNDANT) Terminal Node grounded to PEC ', i,j,k
                      elseif (nodo%islossy) then
                          write (buff,*)  'wir1_INFO: (SHOULD BE REDUNDANT) Terminal Node grounded to Lossy ', i,j,k
                      else
                          write (buff,*)  'wir1_INFO: (SHOULD BE REDUNDANT) Terminal Node embedded in air ', i,j,k
                      endif
                      if ((k >  ZI).and.(k <= ZE)) call WarnErrReport(buff)
                 endif
!!!!fin pedazo de niapa           
               

!!!!! ojo cambio AGRESIVO: sgg 08092016 a peticion OLD. si un nodo es heterogeneousjunction automaticamente no es ni pec ni lossy !esto es delicado !validado con gra_simple_conectado179_162_173.nfde
               if (nodo%isPEC.and.(nodo%IsHeterogeneousJunction)) then
                    nodo%isPEC   = .false.
			        write (buff,*)  'wir1_ERROR: (Deprecated 170220-010324: no longer an error, though I stop. If sure relaunch with -ignoreerrors) ENL/ENR PEC grounded node detached from ground-material for being a wire junction', i,j,k,nodo%ispec,nodo%islossy
                    if ((k >  ZI).and.(k <= ZE)) call WarnErrReport(buff,.true.)
               endif               
               if (nodo%isLossy.and.(nodo%IsHeterogeneousJunction)) then
                    nodo%isLossy   = .false.
                    write (buff,*) 'wir1_ERROR: (Deprecated 170220-010324: no longer an error, though I stop. If sure relaunch with -ignoreerrors)  ENL/ENR Lossy grounded node detached from ground-material for being a wire junction', i,j,k,nodo%ispec,nodo%islossy
                    if ((k >  ZI).and.(k <= ZE)) call WarnErrReport(buff,.true.)
               endif
!!!!!!!!! fin cambio sgg 080912016               
               
               
               if (nodo%isPec) then
                    nodo%isLossy=.false. !tocado 030615 para evitar conflictos en nodos frontera entre Lossy y pec
               endif
               !
            endif
         endif !del fullsize
      end do

   end subroutine detect_peclossyconformal_nodes
   
   subroutine check_embed(embed,tipofieldo,io,jo,ko)
   logical embed
   integer :: ib,jb,kb,tipofieldb
   integer :: io,jo,ko,tipofieldo
   type (CurrentSegments), pointer  ::  dummy
  !!!!!!!!!!!!!!!!! embed=.true.; return !!!!ojoooo sgg tocado a mano para ver bug conformal 220323 
   
    embed=.false.
    if (associated(nodo%CurrentMinus_1))  then
        dummy => nodo%CurrentMinus_1
        call auxem(embed,tipofieldo,io,jo,ko,dummy)  
    endif
    if (associated(nodo%CurrentMinus_2))  then
        dummy => nodo%CurrentMinus_2
        call auxem(embed,tipofieldo,io,jo,ko,dummy)  
    endif
    if (associated(nodo%CurrentMinus_3))  then
        dummy => nodo%CurrentMinus_3
        call auxem(embed,tipofieldo,io,jo,ko,dummy)  
    endif
    if (associated(nodo%CurrentMinus_4))  then
        dummy => nodo%CurrentMinus_4
        call auxem(embed,tipofieldo,io,jo,ko,dummy)  
    endif
    if (associated(nodo%CurrentMinus_5))  then
        dummy => nodo%CurrentMinus_5
        call auxem(embed,tipofieldo,io,jo,ko,dummy)  
    endif
    if (associated(nodo%CurrentMinus_6))  then
        dummy => nodo%CurrentMinus_6
        call auxem(embed,tipofieldo,io,jo,ko,dummy)  
    endif
    if (associated(nodo%CurrentMinus_7))  then
        dummy => nodo%CurrentMinus_7
        call auxem(embed,tipofieldo,io,jo,ko,dummy)  
    endif
    if (associated(nodo%CurrentMinus_8))  then
        dummy => nodo%CurrentMinus_8
        call auxem(embed,tipofieldo,io,jo,ko,dummy)  
    endif
    if (associated(nodo%CurrentMinus_9))  then
        dummy => nodo%CurrentMinus_9
        call auxem(embed,tipofieldo,io,jo,ko,dummy)  
    endif
    !
    if (associated(nodo%CurrentPlus_1))  then
        dummy => nodo%CurrentPlus_1
        call auxem(embed,tipofieldo,io,jo,ko,dummy)  
    endif
    if (associated(nodo%CurrentPlus_2))  then
        dummy => nodo%CurrentPlus_2
        call auxem(embed,tipofieldo,io,jo,ko,dummy)  
    endif
    if (associated(nodo%CurrentPlus_3))  then
        dummy => nodo%CurrentPlus_3
        call auxem(embed,tipofieldo,io,jo,ko,dummy)  
    endif
    if (associated(nodo%CurrentPlus_4))  then
        dummy => nodo%CurrentPlus_4
        call auxem(embed,tipofieldo,io,jo,ko,dummy)  
    endif
    if (associated(nodo%CurrentPlus_5))  then
        dummy => nodo%CurrentPlus_5
        call auxem(embed,tipofieldo,io,jo,ko,dummy)  
    endif
    if (associated(nodo%CurrentPlus_6))  then
        dummy => nodo%CurrentPlus_6
        call auxem(embed,tipofieldo,io,jo,ko,dummy)  
    endif
    if (associated(nodo%CurrentPlus_7))  then
        dummy => nodo%CurrentPlus_7
        call auxem(embed,tipofieldo,io,jo,ko,dummy)  
    endif
    if (associated(nodo%CurrentPlus_8))  then
        dummy => nodo%CurrentPlus_8
        call auxem(embed,tipofieldo,io,jo,ko,dummy)  
    endif
    if (associated(nodo%CurrentPlus_9))  then
        dummy => nodo%CurrentPlus_9
        call auxem(embed,tipofieldo,io,jo,ko,dummy)  
    endif
    
    return
   end subroutine check_embed
   
   subroutine auxem(embed,tipofieldo,io,jo,ko,dummy)
   logical embed
   integer :: ib,jb,kb,tipofieldb
   integer :: io,jo,ko,tipofieldo
   type (CurrentSegments), pointer  ::  dummy
        ib = dummy%i
        jb = dummy%j
        kb = dummy%k
        tipofieldb = dummy%tipofield
        embed = embed .or. ((ib == io).and.(jb == jo).and.(kb == ko).and.(tipofieldo == tipofieldb))
    return
    end subroutine
   
subroutine resume_casuistics

      !In case of resuming a problem, re-starting currents and charges must be read instead of initialized
      !resuming
      if (.not.resume) then
         do i1=1,HWires%NumChargeNodes
            HWires%ChargeNode(i1)%ChargePresent      =0.0_RKIND_wires
            HWires%ChargeNode(i1)%ChargePast         =0.0_RKIND_wires
         end do
         do i1=1,HWires%NumCurrentSegments
            HWires%CurrentSegment(i1)%Current        =0.0_RKIND_wires
            HWires%CurrentSegment(i1)%qplus_qminus         =0.0_RKIND_wires
            HWires%CurrentSegment(i1)%current_for_devia =0.0_RKIND_wires
            HWires%CurrentSegment(i1)%qplus_qminus_for_devia =0.0_RKIND_wires 
            HWires%CurrentSegment(i1)%Efield_main2wire_for_devia =0.0_RKIND_wires 
         end do
#ifdef CompileWithMPI
         !Initialize MPI counters (later written to disk)
         Hwires%NumNeededCurrentUpMPI=0
         Hwires%NumNeededCurrentDownMPI=0
#endif
      else
         do i1=1,HWires%NumChargeNodes
            READ (14) HWires%ChargeNode(i1)%ChargePresent
            READ (14) HWires%ChargeNode(i1)%ChargePast
         end do
         do i1=1,HWires%NumCurrentSegments
            READ (14) HWires%CurrentSegment(i1)%Current
            READ (14) HWires%CurrentSegment(i1)%qplus_qminus
            READ (14) HWires%CurrentSegment(i1)%current_for_devia
            READ (14) HWires%CurrentSegment(i1)%qplus_qminus_for_devia
            READ (14) HWires%CurrentSegment(i1)%Efield_main2wire_for_devia
         end do
#ifdef CompileWithMPI
         READ (14) Hwires%NumNeededCurrentUpMPI,Hwires%NumNeededCurrentDownMPI
         allocate (Hwires%MPIUpNeededCurrentSegment(1 : Hwires%NumNeededCurrentUpMPI))
         allocate (Hwires%MPIDownNeededCurrentSegment(1 : Hwires%NumNeededCurrentDownMPI))
         do i1=1,Hwires%NumNeededCurrentUpMPI
            READ (14) HWires%MPIUpNeededCurrentSegment(i1)%Current
         end do
         do i1=1,Hwires%NumNeededCurrentDownMPI
            READ (14) HWires%MPIDownNeededCurrentSegment(i1)%Current
         end do
         !for new wires MPI march'12 sgg bug en multiwires mpi

#endif
      endif

    end subroutine resume_casuistics
   
    subroutine  init_wirecrank
    
         if (layoutnumber /=0) then
            buff='Unsupported crank in MPI'
            call WarnErrReport(buff,.true.)
         endif
         do n=1,HWires%NumCurrentSegments
            Segmento => HWires%CurrentSegment(n)
            jmed   = Segmento%indexmed
            Resist = Segmento%TipoWire%R
            Autoin =  Segmento%Lind !!!+ Segmento%TipoWire%L !he comentado esto a 110517 porque con lo del fieldtotl ya he metido en %lind tambien la tipowire%L 
            Capaci = 1.0_RKIND_wires / (InvMu(jmed)*InvEps(jmed)* Autoin)
            deltax = Segmento%delta
            !!!!
            segmento%rightCHminus = 1.0_RKIND_wires/(deltax*Capaci)
            segmento%rightCHplus   = - segmento%rightCHminus
            segmento%rightCUminus = (sgg%dt)/(4.0_RKIND_wires*deltax**2*Capaci)
            segmento%rightCUplus  =  segmento%rightCUminus
            !!!!!!
            segmento%upperdiag = -sgg%dt/(4.0_RKIND_wires*deltax**2*Capaci) !comun a todos
            segmento%lowerdiag = segmento%upperdiag                !comun a todos
            !!!
            if ((Segmento%ChargeMinus%NumCurrentMinus == 1).and. &
            (Segmento%ChargeMinus%NumCurrentPlus  == 1).and. &
            (Segmento%ChargePlus%NumCurrentMinus  == 1).and. &
            (Segmento%ChargePlus%NumCurrentPlus   == 1) ) then !es un intermedio
               segmento%diag         =  (sgg%dt/(deltax**2_RKIND_wires*Capaci) + (2_RKIND_wires*Autoin)/sgg%dt + resist)/2.0_RKIND_wires
               segmento%rightCU      = (-sgg%dt/(2.0_RKIND_wires*deltax**2_RKIND_wires*Capaci) + Autoin/sgg%dt - resist/2.0_RKIND_wires)
            else
               segmento%diag    =   ((3_RKIND_wires*sgg%dt)/(2.0_RKIND_wires*deltax**2*Capaci) + (2_RKIND_wires*Autoin)/sgg%dt + resist)/2.0_RKIND_wires  !comun a todos
               segmento%rightCU =  ((-3_RKIND_wires*sgg%dt)/(4.0_RKIND_wires*deltax**2*Capaci) + Autoin/sgg%dt - resist/2.0_RKIND_wires)
               if ((Segmento%ChargeMinus%NumCurrentMinus+Segmento%ChargeMinus%NumCurrentPlus == 1).and. &
               (Segmento%ChargePlus%NumCurrentMinus                                      == 1).and. &
               (Segmento%ChargePlus%NumCurrentPlus                                       == 1)) then
                  segmento%lowerdiag = -1e30 !lo voideo pq no debe tener efecto
               elseif ((Segmento%ChargeMinus%NumCurrentMinus                                     == 1).and. &
               (Segmento%ChargeMinus%NumCurrentPlus                                      == 1).and. &
               (Segmento%ChargePlus%NumCurrentMinus+Segmento%ChargePlus%NumCurrentPlus   == 1)) then
                  segmento%upperdiag = -1e30_RKIND_wires !lo voideo pq no debe tener efecto
               else
                  buff='Unsupported crank'
                  call WarnErrReport(buff,.true.)
               endif
            endif

         end do
         end subroutine init_wirecrank
   
   end subroutine InitWires
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !end of initialization subroutine
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!141118 Permittivity scaling

    subroutine calc_wirehollandconstants(sgg,G2,fieldtotl,wiresflavor,mu00,eps00,simu_devia)

      logical :: simu_devia
      type (SGGFDTDINFO), intent(IN)      ::  sgg
      REAL (KIND=RKIND) , pointer, dimension (:),intent(in)      :: G2
      character(len=*), INTENT(in) :: wiresflavor
      logical, intent(in) :: fieldtotl
      REAL (KIND=RKIND),intent(in)            ::  eps00,mu00
      REAL (KIND=RKIND_wires)           ::  dl
      real      (KIND=RKIND_wires), pointer, dimension(:,:) ::  Den
      integer (kind=4)  :: n,jmed,layoutnumber,iw1,is1,is2,i1,NumParallel,imed

      REAL (KIND=RKIND_wires)           ::  df1,df3,df2,Ddf1,Ddf3,Ddf2,vf1,vf3,vf2,runit

      type (CurrentSegments), pointer  ::  dummy
      type (ChargeNodes), pointer  ::  Nodo
      type (TMultiline), pointer                      ::  Multiline
   

      eps0=eps00; mu0=mu00; !chapuz para convertir la variables de paso en globales

      OldInvEps(0 : sgg%NumMedia)=InvEps(0 : sgg%NumMedia)
      OldInvMu(0 : sgg%NumMedia)=InvMu(0 : sgg%NumMedia)
      InvEps(0 : sgg%NumMedia)=1.0_RKIND_wires/(Eps0*sgg%Med(0 : sgg%NumMedia)%Epr) 
      InvMu(0 : sgg%NumMedia)=1.0_RKIND_wires/(Mu0*sgg%Med(0 : sgg%NumMedia)%Mur) 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Constantes nodales
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do i1=1,HWires%NumChargeNodes
         nodo => HWires%ChargeNode(i1)
         nodo%CtePlain = nodo%CtePlain/Hwires%olddt*sgg%dt 
      end do

!ctes Mur
      do i1=1,HWires%NumChargeNodes
         nodo=>HWires%Chargenode(i1)
         if (nodo%IsBackDownLeftMur) then
            dummy => nodo%CurrentPlus_1
            nodo%cteMur=(sqrt(InvMu (dummy%indexmed)*InvEps(dummy%indexmed))*sgg%dt - dummy%delta)/ &
                     (sqrt(InvMu (dummy%indexmed)*InvEps(dummy%indexmed))*sgg%dt + dummy%delta)
         elseif (nodo%IsFrontUpRightMur) then
            dummy => nodo%CurrentMInus_1
            nodo%cteMur=(sqrt(InvMu (dummy%indexmed)*InvEps(dummy%indexmed))*sgg%dt - dummy%delta)/ &
                     (sqrt(InvMu (dummy%indexmed)*InvEps(dummy%indexmed))*sgg%dt + dummy%delta)
         endif
      end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Constantes de segmentos
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do i1=1,HWires%NumCurrentSegments
         !constantes de actualizacion
         dummy=> HWires%CurrentSegment(i1)
         jmed=dummy%indexmed
         dummy%Lind=dummy%Lind * OldInvMu(jmed)/InvMu(jmed) !escalo la autoinduccion por invmu (aunque el mu no se escala nunca)  
         dummy%Lind_devia=dummy%Lind_devia * OldInvMu(jmed)/InvMu(jmed) !escalo la autoinduccion por invmu (aunque el mu no se escala nunca)       
         call wiresconstantes(fieldtotl,dummy,G2,sgg)
      end do
!!!!!!!!!!!!!!!!!!!!
!Junctions segmento%FractionMinus & segmento%FractionPlus not affected by permittivity scaling

!Solapes con conformal (niapa 171216) no modificados en %cte5

!!!!!!!!!!!!!!! Berenger transition !untested
      !!!copiado 23/04/2014
    if (trim(adjustl(wiresflavor))=='transition') then
         do iw1 = 1,HWires%NumMultilines
            NumParallel = HWires%Multilines(iw1)%NumParallel
            do is1 = 1,NumParallel
               imed 	= HWires%Multilines(iw1)%Segments(is1)%ptr%indexmed
               HWires%Multilines(iw1)%Segments(is1)%ptr%Lintrinsic = HWires%Multilines(iw1)%Segments(is1)%ptr%Lintrinsic  * OldInvMu (imed)/InvMu (imed) 
               HWires%Multilines(iw1)%Segments(is1)%ptr%L      = HWires%Multilines(iw1)%Segments(is1)%ptr%L     * OldInvMu (imed)/InvMu (imed) 
               HWires%Multilines(iw1)%Segments(is1)%ptr%C      = HWires%Multilines(iw1)%Segments(is1)%ptr%C     * OldInvEps(imed)/InvEps(imed) 
               do is2 = 1,NumParallel
                     HWires%Multilines(iw1)%C(is1,is2) = HWires%Multilines(iw1)%Segments(is1)%ptr%C
                     HWires%Multilines(iw1)%L(is1,is2) = HWires%Multilines(iw1)%Segments(is1)%ptr%L
               end do
            end do
            allocate(Den                       (1:NumParallel,1:NumParallel))
            Den = HWires%Multilines(iw1)%L(1:NumParallel,1:NumParallel) + HWires%Multilines(iw1)%R(1:NumParallel,1:NumParallel)*sgg%dt/2.0_RKIND_wires
            call MatInv(NumParallel,Den)

            HWires%Multilines(iw1)%b1I = MatMul(Den, HWires%Multilines(iw1)%L(1:NumParallel,1:NumParallel)  - &
                                                     HWires%Multilines(iw1)%R(1:NumParallel,1:NumParallel)*sgg%dt/2.0_RKIND_wires)
            dl  =  HWires%Multilines(iw1)%Segments(1)%ptr%delta !tomo el de 1 !dama no lo tenia
            imed 	= HWires%Multilines(iw1)%Segments(1)%ptr%indexmed !lo aniado yo pq no estaba definido sgg 141118  !tomo el de 1 !dama no lo tenia
            HWires%Multilines(iw1)%b2I = -sgg%dt/dl*InvMu(imed)*InvEps(imed)*MatMul(Den,HWires%Multilines(iw1)%L)
            HWires%Multilines(iw1)%b3I = sgg%dt*Den

            deallocate(Den)

            do is1 = 1,NumParallel
               HWires%Multilines(iw1)%Segments(is1)%ptr%bI = HWires%Multilines(iw1)%b2I(is1,is1)
            end do
         end do

      ENDIF 

!!!guarda el dt anterior

         Hwires%olddt=sgg%dt


      return


    end subroutine calc_wirehollandconstants


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !Returns info on the adjacency of two segments (first and second)
   !all the wires (Parallel) are collapsed and everything is ohmicly connected if adjacent
   !If -connectendings wires are not collapsed, but the endings are connected even if not  specified by the .nfde
   !Parallel coincident wires with coincident nodes are
   !different wires with different nodes (later will be handled with the Berenger BerengerMTLN routine).
   !segments are collapsed into ONLY 1 -any- and reported accordingly in the warnings file. Furthermore coincident
   !nodes are connected ohmicly -are collapsed also- in this latter case.
   !

   function TestAdjacency(first,numfirst,second,numsecond,connectendings,isolategroupgroups,strictOLD,ZI,ZE,NUMESEG,firstmenos1,FIRSTMAS1,secondmenos1,secondmas1,verbose) RESULT(adj)
      type(CurrentSegments), pointer, intent(in)        ::  first ,second,firstmenos1 ,firstmas1 ,secondmenos1 ,secondmas1
      type(CurrentSegments), pointer     ::  firstprevio,secondprevio
      type (adyacc)  ::  adj
      logical :: verbose
      logical :: conexionados,RRConnected,RLConnected,LLConnected,LRConnected,EndingsConnected,connectendings,isolategroupgroups,RequestedConnection,strictOLD,success,entro1,entro2,entro3
      integer (kind=4) :: numfirst,numsecond,ZI,ZE,void,offx,offy,offz,NUMESEG
      character (len=3), dimension(1:3) :: DIR
      character(len=BUFSIZE) :: buff
      RequestedConnection=.false. ;

      dir(iEx)=' X '
      dir(iEy)=' Y '
      dir(iEz)=' Z '
      adj%YESsegment(1:2) = -1
      adj%Is=.false.
      adj%Parallel=.false.
      adj%IsHeterogeneousJunction=.false.
      adj%BothEndingsConnected = .false.
      if (numfirst == numsecond) return !trivial case: equal segment
      !
      if (first%tipofield==second%tipofield) then
         if    ((first%j==second%j).and.(first%k==second%k).and.(first%tipofield==iEx)) then
            adj%j=first%j
            adj%k=first%k
            if (first%i==second%i-1) then
               adj%Is=.true.
               adj%i=second%i
            endif
            if (first%i==second%i+1) then
               adj%Is=.true.
               adj%i=first%i
            endif
            if (first%i==second%i) then
               adj%Parallel=.true.
               adj%i=first%i
            endif
         elseif ((first%i==second%i).and.(first%k==second%k).and.(first%tipofield==iEy)) then
            adj%i=first%i
            adj%k=first%k
            if (first%j==second%j-1) then
               adj%Is=.true.
               adj%j=second%j
            endif
            if (first%j==second%j+1) then
               adj%Is=.true.
               adj%j=first%j
            endif
            if (first%j==second%j) then
               adj%Parallel=.true.
               adj%j=first%j
            endif
         elseif ((first%i==second%i).and.(first%j==second%j).and.(first%tipofield==iEz)) then
            adj%i=first%i
            adj%j=first%j
            if (first%k==second%k-1) then
               adj%Is=.true.
               adj%k=second%k
            endif
            if (first%k==second%k+1) then
               adj%Is=.true.
               adj%k=first%k
            endif
            if (first%k==second%k) then
               adj%Parallel=.true.
               adj%k=first%k
            endif
         endif
      elseif (((first%tipofield==iEx).and.(second%tipofield==iEy)).and. &
      (((first%i==second%i  ).and.(first%j==second%j  ).and.(first%k==second%k  )).or. &
      ((first%i==second%i  ).and.(first%j==second%j+1).and.(first%k==second%k  )))) then
         adj%Is=.true.
         adj%i=first%i
         adj%j=first%j
         adj%k=first%k
      elseif (((first%tipofield==iEx).and.(second%tipofield==iEy)).and. &
      (((first%i==second%i-1).and.(first%j==second%j  ).and.(first%k==second%k  )).or. &
      ((first%i==second%i-1).and.(first%j==second%j+1).and.(first%k==second%k  )))) then
         adj%Is=.true.
         adj%i=first%i+1
         adj%j=first%j
         adj%k=first%k
         !
      elseif (((first%tipofield==iEy).and.(second%tipofield==iEz)).and. &
      (((first%j==second%j  ).and.(first%k==second%k  ).and.(first%i==second%i  )).or. &
      ((first%j==second%j  ).and.(first%k==second%k+1).and.(first%i==second%i  )))) then
         adj%Is=.true.
         adj%j=first%j
         adj%k=first%k
         adj%i=first%i
      elseif (((first%tipofield==iEy).and.(second%tipofield==iEz)).and. &
      (((first%j==second%j-1).and.(first%k==second%k  ).and.(first%i==second%i  )).or. &
      ((first%j==second%j-1).and.(first%k==second%k+1).and.(first%i==second%i  )))) then
         adj%Is=.true.
         adj%j=first%j+1
         adj%k=first%k
         adj%i=first%i
         !
      elseif (((first%tipofield==iEz).and.(second%tipofield==iEx)).and. &
      (((first%k==second%k  ).and.(first%i==second%i  ).and.(first%j==second%j  )).or. &
      ((first%k==second%k  ).and.(first%i==second%i+1).and.(first%j==second%j  )))) then
         adj%Is=.true.
         adj%k=first%k
         adj%i=first%i
         adj%j=first%j
      elseif (((first%tipofield==iEz).and.(second%tipofield==iEx)).and. &
      (((first%k==second%k-1).and.(first%i==second%i  ).and.(first%j==second%j  )).or. &
      ((first%k==second%k-1).and.(first%i==second%i+1).and.(first%j==second%j  )))) then
         adj%Is=.true.
         adj%k=first%k+1
         adj%i=first%i
         adj%j=first%j
         !
         !
      elseif (((first%tipofield==iEy).and.(second%tipofield==iEx)).and. &
      (((first%i==second%i  ).and.(first%j==second%j  ).and.(first%k==second%k  )).or. &
      ((first%i==second%i+1).and.(first%j==second%j  ).and.(first%k==second%k  )))) then
         adj%Is=.true.
         adj%i=first%i
         adj%j=first%j
         adj%k=first%k
      elseif (((first%tipofield==iEy).and.(second%tipofield==iEx)).and. &
      (((first%i==second%i  ).and.(first%j==second%j-1).and.(first%k==second%k  )).or. &
      ((first%i==second%i+1).and.(first%j==second%j-1).and.(first%k==second%k  )))) then
         adj%Is=.true.
         adj%i=first%i
         adj%j=first%j+1
         adj%k=first%k
         !
      elseif (((first%tipofield==iEx).and.(second%tipofield==iEz)).and. &
      (((first%k==second%k  ).and.(first%i==second%i  ).and.(first%j==second%j  )).or. &
      ((first%k==second%k+1).and.(first%i==second%i  ).and.(first%j==second%j  )))) then
         adj%Is=.true.
         adj%k=first%k
         adj%i=first%i
         adj%j=first%j
      elseif (((first%tipofield==iEx).and.(second%tipofield==iEz)).and. &
      (((first%k==second%k  ).and.(first%i==second%i-1).and.(first%j==second%j  )).or. &
      ((first%k==second%k+1).and.(first%i==second%i-1).and.(first%j==second%j  )))) then
         adj%Is=.true.
         adj%k=first%k
         adj%i=first%i+1
         adj%j=first%j
         !
      elseif (((first%tipofield==iEz).and.(second%tipofield==iEy)).and. &
      (((first%j==second%j  ).and.(first%k==second%k  ).and.(first%i==second%i  )).or. &
      ((first%j==second%j+1).and.(first%k==second%k  ).and.(first%i==second%i  )))) then
         adj%Is=.true.
         adj%j=first%j
         adj%k=first%k
         adj%i=first%i
      elseif (((first%tipofield==iEz).and.(second%tipofield==iEy)).and. &
      (((first%j==second%j  ).and.(first%k==second%k-1).and.(first%i==second%i  )).or. &
      ((first%j==second%j+1).and.(first%k==second%k-1).and.(first%i==second%i  )))) then
         adj%Is=.true.
         adj%j=first%j
         adj%k=first%k+1
         adj%i=first%i
      endif
      !added July'12 to connect requested endings
      if (adj%Parallel) adj%Is = .true.
      !RightConnected adyacency  foreseen in the .nfde file      !added to solve  aNTh-THW_PW bug


      !aniadido 6/2/14 para tratar bien ESN.nfde (solo sirve para dar warnings a OLD con sus rabillos)
      if ((first%Is_LeftEnd).and.(second%Is_LeftEnd)) then
         RequestedConnection = (.not.connectendings).and.(first%indexmed /= second%indexmed).and. &
         ((first%tipowire%LeftEnd == second%tipowire%LeftEnd))
      endif
      if ((first%Is_LeftEnd).and.(second%Is_RightEnd)) then
         RequestedConnection = (.not.connectendings).and.(first%indexmed /= second%indexmed).and. &
         ((first%tipowire%LeftEnd == second%tipowire%RightEnd))
      endif
      if ((first%Is_RightEnd).and.(second%Is_RightEnd))  then
         RequestedConnection = (.not.connectendings).and.(first%indexmed /= second%indexmed).and. &
         ((first%tipowire%RightEnd == second%tipowire%RightEnd))
      endif
      if ((first%Is_RightEnd).and.(second%Is_LeftEnd))  then
         RequestedConnection = (.not.connectendings).and.(first%indexmed /= second%indexmed).and. &
         ((first%tipowire%RightEnd == second%tipowire%LeftEnd))
      endif

      !esto lo he quitado el 6/2/14 pq en ESN.nfde ambos extremos estaban adjancentes pero solo se pedia conexion en uno y por lo tanto hay que discernir con la condicional
      !de mas arriba
      !!!!     RequestedConnection = (.not.connectendings).and.(first%indexmed /= second%indexmed).and. &
      !!!!                                    ((first%tipowire%LeftEnd == second%tipowire%LeftEnd).or. &
      !!!!                                     (first%tipowire%LeftEnd == second%tipowire%RightEnd).or. &
      !!!!                                     (first%tipowire%RightEnd == second%tipowire%LeftEnd).or. &
      !!!!                                     (first%tipowire%RightEnd == second%tipowire%RightEnd))  !solo para reporte bug OLD 12/09/13 Model_unidos en .nfde

      !quito cualquier referencia a LextremoI,j,k porque yo calculo bien los extremos libres tanto si es strictOLD como si no

      LLConnected =  ( first%indexmed /= second%indexmed ).and. &
      ((second%Is_LeftEnd .and. first%Is_LeftEnd).or. &
      (first%IsEnd_norLeft_norRight .and. second%IsEnd_norLeft_norRight).or. &
      (first%IsEnd_norLeft_norRight .and. second%Is_LeftEnd)           .or. &
      (first%Is_LeftEnd .and. second%IsEnd_norLeft_norRight)                 )

      LRConnected = ( first%indexmed /= second%indexmed ).and. &
      ((second%Is_RightEnd .and. first%Is_LeftEnd).or. &
      (first%IsEnd_norLeft_norRight .and. second%IsEnd_norLeft_norRight))

      RLConnected =  ( first%indexmed /= second%indexmed ).and. &
      ((second%Is_LeftEnd .and. first%Is_RightEnd).or. &
      (first%IsEnd_norLeft_norRight .and. second%IsEnd_norLeft_norRight))

      RRConnected = ( first%indexmed /= second%indexmed ).and. &
      ((second%Is_RightEnd .and. first%Is_RightEnd).or. &
      (first%IsEnd_norLeft_norRight .and. second%IsEnd_norLeft_norRight).or. &
      (first%IsEnd_norLeft_norRight .and. second%Is_RightEnd)           .or. &
      (first%Is_RightEnd .and. second%IsEnd_norLeft_norRight)                 )

      if ((first%Is_LeftEnd .and. first%Is_RightEnd).and. (.not.((second%Is_LeftEnd .and. second%Is_RightEnd)))) then
         conexionados = &
         (((second%ilibre  ==  first%ilibre).and. &
         (second%jlibre  ==  first%jlibre).and. &
         (second%klibre  ==  first%klibre)).or. &
         ((second%ilibre  ==  first%i).and. &
         (second%jlibre  ==  first%j).and. &
         (second%klibre  ==  first%k)))

         LLConnected = LLConnected .and. conexionados
         LRConnected = LRConnected .and. conexionados
         RLConnected = RLConnected .and. conexionados
         RRConnected = RRConnected .and. conexionados
      elseif ((second%Is_LeftEnd .and. second%Is_RightEnd).and. (.not.((first%Is_LeftEnd .and. first%Is_RightEnd)))) then
         conexionados = &
         (((second%ilibre  ==  first%ilibre).and. &
         (second%jlibre  ==  first%jlibre).and. &
         (second%klibre  ==  first%klibre)).or. &
         ((second%i       ==  first%ilibre).and. &
         (second%j       ==  first%jlibre).and. &
         (second%k       ==  first%klibre)))

         LLConnected = LLConnected .and. conexionados
         LRConnected = LRConnected .and. conexionados
         RLConnected = RLConnected .and. conexionados
         RRConnected = RRConnected .and. conexionados
      elseif ((second%Is_LeftEnd .and. second%Is_RightEnd).and. (first%Is_LeftEnd .and. first%Is_RightEnd)) then
         conexionados = &
         (((second%ilibre  ==  first%ilibre).and. &
         (second%jlibre  ==  first%jlibre).and. &
         (second%klibre  ==  first%klibre)).or. &
         ((second%i       ==  first%ilibre).and. &
         (second%j       ==  first%jlibre).and. &
         (second%k       ==  first%klibre)).or. &
         ((second%ilibre  ==  first%i).and. &
         (second%jlibre  ==  first%j).and. &
         (second%klibre  ==  first%k)))

         LLConnected = LLConnected .and. conexionados
         LRConnected = LRConnected .and. conexionados
         RLConnected = RLConnected .and. conexionados
         RRConnected = RRConnected .and. conexionados
      else
         conexionados = &
         (second%ilibre  ==  first%ilibre).and. &
         (second%jlibre  ==  first%jlibre).and. &
         (second%klibre  ==  first%klibre)

         LLConnected = LLConnected .and. conexionados
         LRConnected = LRConnected .and. conexionados
         RLConnected = RLConnected .and. conexionados
         RRConnected = RRConnected .and. conexionados
      endif
      if (.not.connectendings) LLConnected=LLConnected.and.(first%TIPOWIRE%LeftEnd==second%TIPOWIRE%LeftEnd)
      if (.not.connectendings) LRConnected=LRConnected.and.(first%TIPOWIRE%LeftEnd==second%TIPOWIRE%RightEnd)
      if (.not.connectendings) RLConnected=RLConnected.and.(first%TIPOWIRE%RightEnd==second%TIPOWIRE%LeftEnd)
      if (.not.connectendings) RRConnected=RRConnected.and.(first%TIPOWIRE%RightEnd==second%TIPOWIRE%RightEnd)


      EndingsConnected =(LLConnected.or.LRConnected.or.RLConnected.or.RRConnected)
      !casos especiales con un hilo de un segmento

      if (adj%is .and. connectendings) then
         void=0
         if (second%Is_LeftEnd) void=void+1
         if (second%Is_RightEnd) void=void+1
         if (second%IsEnd_norLeft_norRight) void=void+1
         if ( (first%Is_LeftEnd .and. first%Is_RightEnd).and.(void==1)) then
            if ( &
            ((first%i      == second%ilibre).and.(first%j      == second%jlibre).and.(first%k      == second%klibre)).or. &
            ((first%ilibre == second%ilibre).and.(first%jlibre == second%jlibre).and.(first%klibre == second%klibre)) ) then
               EndingsConnected = .true.
               adj%i = second%ilibre
               adj%j=  second%jlibre
               adj%k=  second%klibre
            endif
         endif
         void=0
         if (first%Is_LeftEnd) void=void+1
         if (first%Is_RightEnd) void=void+1
         if (first%IsEnd_norLeft_norRight) void=void+1
         if ( (void ==1).and.(second%Is_LeftEnd .and. second%Is_RightEnd)) then
            if ( &
            ((first%ilibre == second%i     ).and.(first%jlibre == second%j     ).and.(first%klibre == second%k     )).or. &
            ((first%ilibre == second%ilibre).and.(first%jlibre == second%jlibre).and.(first%klibre == second%klibre)) ) then
               EndingsConnected = .true.
               adj%i = first%ilibre
               adj%j=  first%jlibre
               adj%k=  first%klibre
            endif
         endif
      endif
      !los cuatro conectados (puede pasar en amelet o en nfde)
      if ((first%Is_LeftEnd .and. first%Is_RightEnd).and.(second%Is_LeftEnd .and. second%Is_RightEnd)) then
         if (connectendings) then
            if ( &
            ((first%i      == second%ilibre).and.(first%j      == second%jlibre).and.(first%k      == second%klibre)).or. &
            ((first%ilibre == second%ilibre).and.(first%jlibre == second%jlibre).and.(first%klibre == second%klibre)) ) then
               EndingsConnected = .true.
               adj%i = second%ilibre
               adj%j=  second%jlibre
               adj%k=  second%klibre
            elseif ( &
            ((first%ilibre == second%i     ).and.(first%jlibre == second%j     ).and.(first%klibre == second%k     )).or. &
            ((first%ilibre == second%ilibre).and.(first%jlibre == second%jlibre).and.(first%klibre == second%klibre)) ) then
               EndingsConnected = .true.
               adj%i = first%ilibre
               adj%j=  first%jlibre
               adj%k=  first%klibre
            elseif ( &
            ((first%i == second%i     ).and.(first%j == second%j     ).and.(first%k == second%k     )) ) then
               EndingsConnected = .true.
               adj%i = first%i
               adj%j=  first%j
               adj%k=  first%k
            elseif ( &
            ((first%ilibre == second%ilibre).and.(first%jlibre == second%jlibre).and.(first%klibre == second%klibre)) ) then
               EndingsConnected = .true.
               adj%i = first%ilibre
               adj%j=  first%jlibre
               adj%k=  first%klibre
            endif
            !
            adj%BothEndingsConnected =adj%Parallel
         else
            EndingsConnected =(LLConnected.or.LRConnected.or.RLConnected.or.RRConnected) !redundante ya hecho antes
            adj%BothEndingsConnected = ((RRConnected.and.LLConnected).or.(RLConnected.and.LRConnected))
         endif
      endif


      !fin casos especiales
      if (adj%is) then
            if  (first%indexmed /= second%indexmed) then
               if  (.not.EndingsConnected) then
                  ADJ%IS=.FALSE.
                  adj%IsHeterogeneousJunction=.false.
                  if (adj%Parallel) then
                     if ((second%Is_RightEnd .or. second%Is_LeftEnd .or. second%IsEnd_norLeft_norRight).and. &
                     (first%Is_RightEnd .or. first%Is_LeftEnd .or. first%IsEnd_norLeft_norRight)) THEN
                        if (RequestedConnection) then
                           write (buff,'(a)')  'wir2_ERROR: Requested connection on non-connected Parallel Adjacent ENDING segments from multiWIREs:  '
                           if ((first%k >= ZI).and.(first%k <= ZE)) call WarnErrReport(buff,.true.)
                           write (buff,'(i7,3i7,a,i7,3i7,a)') first%origindex,first%i,first%j,first%k,' '//dir(first%tipofield),&
                           second%origindex,second%i,second%j,second%k,' '//dir(second%tipofield)
                           if ((first%k >= ZI).and.(first%k <= ZE)) call WarnErrReport(buff,.true.)
                        else
                           write (buff,'(a)')  'wir2_WARNING: DISCONNECTING Parallel Adjacent ENDING segments from multiWIREs:  '
                           if ((first%k >= ZI).and.(first%k <= ZE)) call WarnErrReport(buff)
                           write (buff,'(i7,3i7,a,i7,3i7,a)') first%origindex,first%i,first%j,first%k,' '//dir(first%tipofield),&
                           second%origindex,second%i,second%j,second%k,' '//dir(second%tipofield)
                           if ((first%k >= ZI).and.(first%k <= ZE)) call WarnErrReport(buff)
                        endif
                     ELSE
                        write (buff,'(a)')  'wir2_INFO: DISCONNECTING Parallel Adjacent intermediate segments from multiWIREs:  '
                        if ((first%k >= ZI).and.(first%k <= ZE).and.verbose) call WarnErrReport(buff)
                        write (buff,'(i7,3i7,a,i7,3i7,a)') first%origindex,first%i,first%j,first%k,' '//dir(first%tipofield),&
                        second%origindex,second%i,second%j,second%k,' '//dir(second%tipofield)
                        if ((first%k >= ZI).and.(first%k <= ZE).and.verbose) call WarnErrReport(buff)
                     ENDIF
                  else
                     if ((second%Is_RightEnd .or. second%Is_LeftEnd .or. second%IsEnd_norLeft_norRight).and. &
                     (first%Is_RightEnd .or. first%Is_LeftEnd .or. first%IsEnd_norLeft_norRight)) THEN
                        if (RequestedConnection) then
                           write (buff,'(a)')   'wir2_ERROR: Requested connection on non-connected Non-Parallel Adjacent ENDING segments from multiWIREs:  '
                           if ((first%k >= ZI).and.(first%k <= ZE)) call WarnErrReport(buff)
                           write (buff,'(i7,3i7,a,i7,3i7,a)') first%origindex,first%i,first%j,first%k,' '//dir(first%tipofield),&
                           second%origindex,second%i,second%j,second%k,' '//dir(second%tipofield)
                           if ((first%k >= ZI).and.(first%k <= ZE)) call WarnErrReport(buff,.true.)
                        else
                           write (buff,'(a)')  'wir2_WARNING: DISCONNECTING NON-Parallel Adjacent ENDING segments from multiWIREs:  '
                           if ((first%k >= ZI).and.(first%k <= ZE)) call WarnErrReport(buff)
                           write (buff,'(i7,3i7,a,i7,3i7,a)') first%origindex,first%i,first%j,first%k,' '//dir(first%tipofield),&
                           second%origindex,second%i,second%j,second%k,' '//dir(second%tipofield)
                           if ((first%k >= ZI).and.(first%k <= ZE)) call WarnErrReport(buff)
                        endif
                     else
                        write (buff,'(a)')  'wir2_INFO: DISCONNECTING Non-Parallel Adjacent intermediate segments from multiWIREs:'
                        if ((first%k >= ZI).and.(first%k <= ZE).and.verbose) call WarnErrReport(buff)
                        write (buff,'(i7,3i7,a,i7,3i7,a)') first%origindex,first%i,first%j,first%k,' '//dir(first%tipofield),&
                        second%origindex,second%i,second%j,second%k,' '//dir(second%tipofield)
                        if ((first%k >= ZI).and.(first%k <= ZE).and.verbose) call WarnErrReport(buff)
                     ENDIF
                  endif
               elseif (EndingsConnected) then
                  if (adj%Parallel) then
                     write (buff,'(a)')  'wir2_INFO: CONNECTING Parallel Adjacent ENDING segments from multiWIREs:  '
                  else
                     write (buff,'(a)')  'wir2_INFO: CONNECTING Non-Parallel Adjacent ENDING segments from multiWIREs:'
                  endif
                  if ((first%k >= ZI).and.(first%k <= ZE).and.verbose) call WarnErrReport(buff)
                  write (buff,'(i7,3i7,a,i7,3i7,a)') first%origindex,first%i,first%j,first%k,' '//dir(first%tipofield),&
                  second%origindex,second%i,second%j,second%k,' '//dir(second%tipofield)
                  if ((first%k >= ZI).and.(first%k <= ZE).and.verbose) call WarnErrReport(buff)
                  write (buff,'(a,3i7)')  '           AT :',adj%i,adj%j,adj%k
                  if ((first%k >= ZI).and.(first%k <= ZE).and.verbose) call WarnErrReport(buff)

                  adj%IsHeterogeneousJunction=.true.
                  ADJ%IS=.TRUE.
                  adj%YESsegment(1) = numfirst
                  adj%YESsegment(2) = numsecond
                  if (LLConnected.or.LRConnected.or.RLConnected.or.RRConnected) then
!!!!esto daba problemas con edel_wuking2.nfde y en general con el ultimo cable si solo tiene un solo segmento no unido, poniendo nodo%heterogeneousjunction=.true. al nodo libre en la vuelta (numfirst>numsecond)
!!!!comentado, pues, el 141218 y aniadido codigo 
                         !!adj%i = first%ilibre
                         !!adj%j=  first%jlibre
                         !!adj%k=  first%klibre
                    entro1=.false.; entro2=.false.; entro3=.false.
!!!lo de antes descomentado no daba problemas OLD. En correo 150219 entro el primer buggy de abajo. Por tanto aniado mas
!!!casuistica. Ahora puede que falle el edel_wuking2.nfde por este cambio, pero prefiero no descomentar lo anterior. Arreglar cuando pase !!
                    if ( &
                    ((first%i      == second%ilibre).and.(first%j      == second%jlibre).and.(first%k      == second%klibre))) then
                        adj%i = second%ilibre
                        adj%j=  second%jlibre
                        adj%k=  second%klibre
                        entro1=.true.
                    endif
                    if ( &
                    ((first%ilibre == second%i     ).and.(first%jlibre == second%j     ).and.(first%klibre == second%k     ))) then
                     adj%i = first%ilibre
                     adj%j=  first%jlibre
                     adj%k=  first%klibre
                        entro2=.true.
                    endif
!!!esta es la casuistica que aniado 150219                   
                    if ( &
                    ((first%ilibre     == second%ilibre).and.(first%jlibre      == second%jlibre).and.(first%klibre      == second%klibre))) then
                        adj%i = second%ilibre
                        adj%j=  second%jlibre
                        adj%k=  second%klibre
                        entro3=.true.
                    endif
!!!fin casuistica que anidado 150219
                    if ((.not.entro1).and.(.not.entro2).and.(.not.entro3)) then
                            write (buff,*)  'wir2_BUGGYERROR: bug in adjacencies entro. ',first%ilibre,first%jlibre,first%klibre,second%ilibre,second%jlibre,second%klibre
                            call WarnErrReport(buff,.true.)
                    endif
                    if (entro1.and.entro2) then
                        if (.not.entro3) then
                            write (buff,'(a)')  'wir2_BUGGYERROR: bug in adjacencies entro. ',first%ilibre,first%jlibre,first%klibre,second%ilibre,second%jlibre,second%klibre
                            call WarnErrReport(buff,.true.)
                        endif
                    endif
!!!!fina aniadido 141218
                  ENDIF
               endif
            ELSEIF (First%indexmed == second%indexmed) THEN
               if (adj%Parallel) then !paralelos del mismo hilo
                  if (.not.strictOLD) then
                     adj%is=.false.
                     !NUNCA DEBERIA ENTRAR AQUI porque los paralelos se han colapaso en la version no estricta
                     write (buff,'(a)')  'wir2_BUGGYERROR: DISCONNECTING Parallel segments from the same WIRE:'
                     call WarnErrReport(buff,.true.)
                     write (buff,'(i7,3i7,a,i7,3i7,a)') first%origindex,first%i,first%j,first%k,' '//dir(first%tipofield),&
                     second%origindex,second%i,second%j,second%k,' '//dir(second%tipofield)
                     call WarnErrReport(buff,.true.)
                  else
                     !!!!!!!!!!aqui es donde viene el meollo porque NO he quitado segmentos repetidos
                     !voy a suponer que forman una cadena ordenada

                     if  (abs(numfirst - numsecond)>1) then
                        adj%is=.false.
                        adj%IsHeterogeneousJunction=.false.
                        write (buff,'(a)')  'wir2_INFO: DISCONNECTING NON-CORRELATIVE Parallel segments from the same WIRE:'
                        if (verbose) call WarnErrReport(buff)
                        write (buff,'(i7,3i7,a,i7,3i7,a)') first%origindex,first%i,first%j,first%k,' '//dir(first%tipofield),&
                        second%origindex,second%i,second%j,second%k,' '//dir(second%tipofield)
                        if (verbose) call WarnErrReport(buff)
                     else
                        sucCess=.false.
                        firstprevio => null ()
                        secondprevio => null ()
                        if (associated(firstmenos1)) then
                           if (.not.((firstmenos1%i == first%i).and.(firstmenos1%j == first%j).and.(firstmenos1%k == first%k).and.(firstmenos1%tipofield == first%tipofield))) then
                              firstprevio=>firstmenos1
                           else
                              if (associated(firstmas1)) then
                                 if (.not.((firstmas1%i == first%i).and.(firstmas1%j == first%j).and.(firstmas1%k == first%k).and.(firstmas1%tipofield == first%tipofield))) then
                                    firstprevio=>firstmas1
                                 endif
                              endif
                           endif
                        else
                           if (associated(firstmas1)) then
                              if (.not.((firstmas1%i == first%i).and.(firstmas1%j == first%j).and.(firstmas1%k == first%k).and.(firstmas1%tipofield == first%tipofield))) then
                                 firstprevio=>firstmas1
                              endif
                           endif
                        endif
                        if (associated(secondmenos1)) then
                           if (.not.((secondmenos1%i == second%i).and.(secondmenos1%j == second%j).and.(secondmenos1%k == second%k).and.(secondmenos1%tipofield == second%tipofield))) then
                              secondprevio=>secondmenos1
                           else
                              if (associated(secondmas1)) then
                                 if (.not.((secondmas1%i == second%i).and.(secondmas1%j == second%j).and.(secondmas1%k == second%k).and.(secondmas1%tipofield == second%tipofield))) then
                                    secondprevio=>secondmas1
                                 endif
                              endif
                           endif
                        else
                           if (associated(secondmas1)) then
                              if (.not.((secondmas1%i == second%i).and.(secondmas1%j == second%j).and.(secondmas1%k == second%k).and.(secondmas1%tipofield == second%tipofield))) then
                                 secondprevio=>secondmas1
                              endif
                           endif
                        endif
                        !
                        offx=0; offy=0; offz=0
                        if ((first%tipofield==second%tipofield).and.(first%tipofield==iEx)) offx=1
                        if ((first%tipofield==second%tipofield).and.(first%tipofield==iEy)) offy=1
                        if ((first%tipofield==second%tipofield).and.(first%tipofield==iEz)) offz=1
                        if (associated(firstprevio)) then
                           select case(first%tipofield)
                            case  (iEx)
                              if (firstprevio%i == first%i) then
                                 success=.true.
                                 adj%i=first%i+1
                                 adj%j=first%j
                                 adj%k=first%k
                              elseif ((firstprevio%i - first%i)==1) then
                                 success=.true.
                                 adj%i=first%i
                                 adj%j=first%j
                                 adj%k=first%k
                              elseif ((firstprevio%i - first%i)==-1) then
                                 success=.true.
                                 adj%i=first%i+offx
                                 adj%j=first%j
                                 adj%k=first%k
                              endif
                            case  (iEy)
                              if (firstprevio%j == first%j) then
                                 success=.true.
                                 adj%i=first%i
                                 adj%j=first%j+1
                                 adj%k=first%k
                              elseif ((firstprevio%j - first%j)==1) then
                                 success=.true.
                                 adj%i=first%i
                                 adj%j=first%j
                                 adj%k=first%k
                              elseif ((firstprevio%j - first%j)==-1) then
                                 success=.true.
                                 adj%i=first%i
                                 adj%j=first%j+offy
                                 adj%k=first%k
                              endif
                            case  (iEz)
                              if (firstprevio%k == first%k) then
                                 success=.true.
                                 adj%i=first%i
                                 adj%j=first%j
                                 adj%k=first%k+1
                              elseif ((firstprevio%k - first%k)==1) then
                                 success=.true.
                                 adj%i=first%i
                                 adj%j=first%j
                                 adj%k=first%k
                              elseif ((firstprevio%k - first%k)==-1) then
                                 success=.true.
                                 adj%i=first%i
                                 adj%j=first%j
                                 adj%k=first%k+offz
                              endif
                           end select
                        elseif (associated(secondprevio)) then
                           select case(second%tipofield)
                            case  (iEx)
                              if (secondprevio%i == second%i) then
                                 success=.true.
                                 adj%i=second%i+1
                                 adj%j=second%j
                                 adj%k=second%k
                              elseif ((secondprevio%i - second%i)==1) then
                                 success=.true.
                                 adj%i=second%i
                                 adj%j=second%j
                                 adj%k=second%k
                              elseif ((secondprevio%i - second%i)==-1) then
                                 success=.true.
                                 adj%i=second%i+offx
                                 adj%j=second%j
                                 adj%k=second%k
                              endif
                            case  (iEy)
                              if (secondprevio%j == second%j) then
                                 success=.true.
                                 adj%i=second%i
                                 adj%j=second%j+1
                                 adj%k=second%k
                              elseif ((secondprevio%j - second%j)==1) then
                                 success=.true.
                                 adj%i=second%i
                                 adj%j=second%j
                                 adj%k=second%k
                              elseif ((secondprevio%j - second%j)==-1) then
                                 success=.true.
                                 adj%i=second%i
                                 adj%j=second%j+offy
                                 adj%k=second%k
                              endif
                            case  (iEz)
                              if (secondprevio%k == second%k) then
                                 success=.true.
                                 adj%i=second%i
                                 adj%j=second%j
                                 adj%k=second%k+1
                              elseif ((secondprevio%k - second%k)==1) then
                                 success=.true.
                                 adj%i=second%i
                                 adj%j=second%j
                                 adj%k=second%k
                              elseif ((secondprevio%k - second%k)==-1) then
                                 success=.true.
                                 adj%i=second%i
                                 adj%j=second%j
                                 adj%k=second%k+offz
                              endif
                           end select
                        else
                           if ((first%Is_LeftEnd .or. first%Is_RightEnd).and.(second%Is_LeftEnd .or. second%Is_RightEnd)) then
                              success=.true.
                              adj%i=second%i+offx
                              adj%j=second%j+offy
                              adj%k=second%k+offz
                           endif
                        endif
                        if (success) then
                           adj%IsHeterogeneousJunction=.false.
                           ADJ%IS=.TRUE.
                           adj%YESsegment(1) = numfirst
                           adj%YESsegment(2) = numsecond
                           write (buff,'(a)')  'wir2_INFO: CONNECTING CORRELATIVE Parallel segments from the same WIRE (rabitos):'
                           if (verbose) call WarnErrReport(buff)
                           write (buff,'(i7,3i7,a,i7,3i7,a)') first%origindex,first%i,first%j,first%k,' '//dir(first%tipofield),&
                           second%origindex,second%i,second%j,second%k,' '//dir(second%tipofield)
                           if (verbose) call WarnErrReport(buff)
                           write (buff,'(a,3i7)')  '           AT :',adj%i,adj%j,adj%k
                           if (verbose) call WarnErrReport(buff)
                        else
                           adj%IsHeterogeneousJunction=.false.
                           ADJ%IS=.false.
                           write (buff,'(a)')  'wir2_BUGGYERROR:  Cannot determine point of contact of Parallel intra-WIRE segment connection (mas de dos rabitos doblados?). '
                           call WarnErrReport(buff,.true.)
                           write (buff,'(i7,3i7,a,i7,3i7,a)') first%origindex,first%i,first%j,first%k,' '//dir(first%tipofield),&
                           second%origindex,second%i,second%j,second%k,' '//dir(second%tipofield)
                           call WarnErrReport(buff,.true.)
                        endif
                        !
                     endif
                  endif !DEL NO STRICTO
               else !no son paralelos pero si estan en el mismo hilko
                  if (.not.strictOLD) then
                     adj%IsHeterogeneousJunction=.false.
                     ADJ%IS=.TRUE.
                     adj%YESsegment(1) = numfirst
                     adj%YESsegment(2) = numsecond
                     if (abs(numfirst - numsecond)>1) then
                        write (buff,'(a)')  'wir2_INFO: CONNECTING NON-CORRELATIVE Non-Parallel segments from the same WIRE:'
                        if (verbose) call WarnErrReport(buff)
                     ELSE
                        write (buff,'(a)')  'wir2_INFO: CONNECTING CORRELATIVE Non-Parallel segments from the same WIRE:'
                        if (verbose) call WarnErrReport(buff)
                     ENDIF
                     write (buff,'(i7,3i7,a,i7,3i7,a)') first%origindex,first%i,first%j,first%k,' '//dir(first%tipofield),&
                     second%origindex,second%i,second%j,second%k,' '//dir(second%tipofield)
                     call WarnErrReport(buff)
                     write (buff,'(a,3i7)')  '           AT :',adj%i,adj%j,adj%k
                     call WarnErrReport(buff)
                  else
                     if (abs(numfirst - numsecond)>1) then !solo si estan leidos contiguamente se toman como adyacentes
                        adj%IsHeterogeneousJunction=.false.
                        ADJ%IS=.false.
                        write (buff,'(a)')  'wir2_INFO: DISCONNECTING NON-CORRELATIVE Non-Parallel segments from the same WIRE:'
                        if (verbose) call WarnErrReport(buff)  !demasiado verbose
                        write (buff,'(i7,3i7,a,i7,3i7,a)') first%origindex,first%i,first%j,first%k,' '//dir(first%tipofield),&
                        second%origindex,second%i,second%j,second%k,' '//dir(second%tipofield)
                        if (verbose) call WarnErrReport(buff)
                     else
                        adj%IsHeterogeneousJunction=.false.
                        ADJ%IS=.TRUE.
                        adj%YESsegment(1) = numfirst
                        adj%YESsegment(2) = numsecond
                        write (buff,'(a)')  'wir2_INFO: CONNECTING CORRELATIVE Non-Parallel segments from the same WIRE:'
                        if (verbose) call WarnErrReport(buff) !demasiado verbose
                        write (buff,'(i7,3i7,a,i7,3i7,a)') first%origindex,first%i,first%j,first%k,' '//dir(first%tipofield),&
                        second%origindex,second%i,second%j,second%k,' '//dir(second%tipofield)
                        if (verbose) call WarnErrReport(buff)
                        write (buff,'(a,3i7)')  '           AT :',adj%i,adj%j,adj%k
                        if (verbose) call WarnErrReport(buff)
                     endif
                  endif
               endif
            endif
         endif
      !
      !!!!!!!!!isolategroupsgroups
      if (isolategroupgroups) then
         if ((adj%is).AND.(first%tipowire%LeftEnd /=  second%tipowire%LeftEnd)) then
            write (buff,'(a)')  'wir2_WARNING: DISCONNECTING PREVIOUSLY Connected segments from multiWIREs being in DIFFERENT GROUPGROUPS:  '
            if ((first%k >= ZI).and.(first%k <= ZE)) call WarnErrReport(buff)
            write (buff,'(i7,3i7,a,i7,3i7,a)') first%tipowire%LeftEnd ,first%origindex,first%i,first%j,first%k,' '//dir(first%tipofield),&
            second%tipowire%LeftEnd,second%origindex,second%i,second%j,second%k,' '//dir(second%tipofield)
            if ((first%k >= ZI).and.(first%k <= ZE)) call WarnErrReport(buff)
            adj%is=.false.
            adj%BothEndingsConnected=.false.
         endif
      endif
      !!!!!!!!!!
      return
   end function TestAdjacency



   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !Advancing charge and current routine
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   subroutine AdvanceWiresE(sgg,timeinstant, layoutnumber,wiresflavor,simu_devia,stochastic,experimentalVideal,wirethickness,eps0,mu0)
                    
      REAL (KIND=RKIND), intent(IN)           ::  eps0,mu0
      integer, intent(IN) :: wirethickness
      logical :: simu_devia,stochastic,experimentalVideal
      type (SGGFDTDINFO), intent(IN)      ::  sgg

      integer (kind=4)  :: n,jmed,layoutnumber,iw1,is1,is2

      integer (kind=4), intent(IN)  ::  timeinstant
      REAL (KIND=RKIND_wires)   ::  Iplus,IMinus,Qplus,QMinus,timei
      REAL (KIND=RKIND_wires)   ::  Vincid,Iincid
      type (CurrentSegments), pointer  ::  Segmento, Segmento2
      type (ChargeNodes), pointer  ::  Nodo
      type (TMultiline), pointer                      ::  Multiline
      character(len=*), INTENT(in) :: wiresflavor
      timei = sgg%tiempo(timeinstant) 
      !!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !FIRST ADVANCE THE CHARGE from n+1.0_RKIND_wires / 2 to n+3/2 using the current known at n+1
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      !it is important that charge is first advanced for time-stepping coherency
      !Crucial also for MPI correctness!!!,
      !Beware than for MPI some currents outside the layout are also calculated with wrong information from charges, but
      !the MPI exchange call will overwrite them with the correct ones coming from the adjacent layers
      !WARNING: MPI does not handle correctly PERIODIC MIRRORING !\E7 20sept11 in currents !to dooooooooo

!!!!!!!!!!!!!!!051223 no pongo a cero los already_YEEadvanced_byconformal_changedtoPECfield !No son cero, sino lo que se updatee mas lo que se inyecte
!!!#ifdef CompileWithOpenMP
!!!!$OMP PARALLEL DO DEFAULT(SHARED)  private(Nodo)
!!!#endif
!!!      do n=1,HWires%NumChargeNodes
!!!         Nodo => HWires%ChargeNode(n)
!!!!!!!140220 pon a PEC los viejos notouch=already_YEEadvanced_byconformal_changedtoPECfield que tenga conectado un nodo
!!!         !debe ser lo primero que se hace para overridear el call conformal_advance_E
!!!         if (associated(nodo%already_YEEadvanced_byconformal_changedtoPECfield1)) then
!!!             nodo%already_YEEadvanced_byconformal_changedtoPECfield1=0.0_RKIND
!!!         endif
!!!         if (associated(nodo%already_YEEadvanced_byconformal_changedtoPECfield2)) then
!!!             nodo%already_YEEadvanced_byconformal_changedtoPECfield2=0.0_RKIND
!!!         endif
!!!         if (associated(nodo%already_YEEadvanced_byconformal_changedtoPECfield3)) then
!!!             nodo%already_YEEadvanced_byconformal_changedtoPECfield3=0.0_RKIND
!!!         endif
!!!         if (associated(nodo%already_YEEadvanced_byconformal_changedtoPECfield4)) then
!!!             nodo%already_YEEadvanced_byconformal_changedtoPECfield4=0.0_RKIND
!!!         endif
!!!         if (associated(nodo%already_YEEadvanced_byconformal_changedtoPECfield5)) then
!!!             nodo%already_YEEadvanced_byconformal_changedtoPECfield5=0.0_RKIND
!!!         endif
!!!         if (associated(nodo%already_YEEadvanced_byconformal_changedtoPECfield6)) then
!!!             nodo%already_YEEadvanced_byconformal_changedtoPECfield6=0.0_RKIND
!!!         endif
!!!      end do
!!!#ifdef CompileWithOpenMP
!!!!$OMP END PARALLEL DO
!!!#endif
!
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED)  private(Iplus,IMinus,Nodo)
#endif
      do n=1,HWires%NumChargeNodes
         Nodo => HWires%ChargeNode(n)
!
         if (nodo%exists) then
            Nodo%ChargePast=Nodo%ChargePresent
            !I sum up the plus and the minus currents accordingly in advance
            Iplus = 0.0_RKIND_wires
            IMinus = 0.0_RKIND_wires
            If (Nodo%NumCurrentPlus>=1) Iplus = Nodo%CurrentPlus_1%Current
            If (Nodo%NumCurrentPlus>=2) Iplus = Iplus+Nodo%CurrentPlus_2%Current
            If (Nodo%NumCurrentPlus>=3) Iplus = Iplus+Nodo%CurrentPlus_3%Current
            If (Nodo%NumCurrentPlus>=4) Iplus = Iplus+Nodo%CurrentPlus_4%Current
            If (Nodo%NumCurrentPlus>=5) Iplus = Iplus+Nodo%CurrentPlus_5%Current
            If (Nodo%NumCurrentPlus>=6) Iplus = Iplus+Nodo%CurrentPlus_6%Current
            If (Nodo%NumCurrentPlus>=7) Iplus = Iplus+Nodo%CurrentPlus_7%Current
            If (Nodo%NumCurrentPlus>=8) Iplus = Iplus+Nodo%CurrentPlus_8%Current
            If (Nodo%NumCurrentPlus>=9) Iplus = Iplus+Nodo%CurrentPlus_9%Current
            !
            If (Nodo%NumCurrentMinus>=1) IMinus = Nodo%CurrentMinus_1%Current
            If (Nodo%NumCurrentMinus>=2) IMinus = IMinus+Nodo%CurrentMinus_2%Current
            If (Nodo%NumCurrentMinus>=3) IMinus = IMinus+Nodo%CurrentMinus_3%Current
            If (Nodo%NumCurrentMinus>=4) IMinus = IMinus+Nodo%CurrentMinus_4%Current
            If (Nodo%NumCurrentMinus>=5) IMinus = IMinus+Nodo%CurrentMinus_5%Current
            If (Nodo%NumCurrentMinus>=6) IMinus = IMinus+Nodo%CurrentMinus_6%Current
            If (Nodo%NumCurrentMinus>=7) IMinus = IMinus+Nodo%CurrentMinus_7%Current
            If (Nodo%NumCurrentMinus>=8) IMinus = IMinus+Nodo%CurrentMinus_8%Current
            If (Nodo%NumCurrentMinus>=9) IMinus = IMinus+Nodo%CurrentMinus_9%Current

            if ((Nodo%NumCurrentMinus == 1).and.(Nodo%NumCurrentPlus == 0)) then
                if (Nodo%IsPeriodic) then
                    Iplus = + Iminus
                else
            !The node is a true terminal one and the mirror of the current is employed for updating (Edelvik's treatment similar to PMC)
                    Iplus = -Iminus
                endif
            endif
            if ((Nodo%NumCurrentMinus == 0).and.(Nodo%NumCurrentPlus == 1)) then
                if (Nodo%IsPeriodic) then
                    IMinus = +Iplus
                else
                    IMinus = -Iplus
                endif
            endif
            !Algoritmo comun a toda la casuistica !feb 13
            if (.not.nodo%IsMur)   then
               Nodo%ChargePresent =  Nodo%CteProp*Nodo%ChargePast   - Nodo%CtePlain*(Iplus-IMinus)
            else
                continue
            endif
            !        if (nodo%IsAttachedtoVoltage) Nodo%ChargePresent =  0.0_RKIND_wires

         endif !del nodo%exists !voideo algunos cuando spliteo
      end do
!!!no toco lo que esta comentado a 230323 para fuentes duras/blandas
      !!!las convierto en duras mas abajo. correos jag rayos junio'15
      !!!Transparent current source feeding in the charge nodes if necessary
      !!!if (.not.simu_devia) then              
      !!!if (thereAreIsources) then
      !!    do n=1,HWires%NumChargeNodes
      !!        if (HWires%ChargeNode(n)%exists) then
      !!        If  (HWires%ChargeNode(n)%HasIsource) then
      !!            Nodo => HWires%ChargeNode(n)
      !!            Iincid=evolucion(timei-unmedio*sgg%dt,Nodo%Isource%Fichero%Samples, &
      !!                              Nodo%Isource%Fichero%DeltaSamples,Nodo%Isource%Fichero%NumSamples)
      !!            Nodo%ChargePresent                         =Nodo%ChargePresent                          +  Nodo%CtePlain                          * Iincid
      !!!!!!!!no funciona esta idea previa de fuente de corriente  quizas porque violo kirchhoff
      !!            Nodo%CurrentPlus_1%Chargeplus%ChargePresent=Nodo%CurrentPlus_1%Chargeplus%ChargePresent -  Nodo%CurrentPlus_1%Chargeplus%CtePlain * Iincid
      !!        endif
      !!        endif
      !!    end do
      !!!endif
      !!!endif
      
!aniado eso a 230323 para fuentes de corrientes duras blandas. Me inspiro en lo que esta comentado antes   
      if (.not.simu_devia) then              
          if (thereAreIsources) then
              do n=1,HWires%NumChargeNodes
                  if (HWires%ChargeNode(n)%exists) then
                  If  (HWires%ChargeNode(n)%HasIsource) then
                      Nodo => HWires%ChargeNode(n)
                      Iincid=evolucion(timei-unmedio*sgg%dt,Nodo%Isource%Fichero%Samples, &
                                        Nodo%Isource%Fichero%DeltaSamples,Nodo%Isource%Fichero%NumSamples)
                      if (Nodo%Isource%soft) then
                            Nodo%ChargePresent = Nodo%ChargePresent    +  Nodo%CtePlain     * Iincid
                      else
                          Vincid=iincid !es realmente una fuente de carga
                          if (associated(Nodo%CurrentPlus_1)) then
                            Nodo%ChargePresent = Vincid/ (Nodo%CurrentPlus_1%Lind * InvMu(Nodo%CurrentPlus_1%indexmed)*InvEps(Nodo%CurrentPlus_1%indexmed))
                          elseif (associated(Nodo%CurrentMinus_1)) then !alguno estara asociado, solo es para sacar el L, eps, mu
                            Nodo%ChargePresent = Vincid/ (Nodo%CurrentMinus_1%Lind * InvMu(Nodo%CurrentMinus_1%indexmed)*InvEps(Nodo%CurrentMinus_1%indexmed))
                          else
                              print *,'bug in current sources '
                              stop
                          endif
                      endif
                  endif
                  endif
              end do
          endif
      endif

      !Absorbing Mur boundary conditions if necessary in the charges
      if (thereAreMurConditions) then
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED)  private(Nodo)
#endif
         do n=1,HWires%NumChargeNodes
            Nodo => HWires%ChargeNode(n)
            if (nodo%exists.and.nodo%IsMur) then
               Nodo%ChargePresent=Nodo%NodeInside%ChargePast + Nodo%cteMur*(Nodo%NodeInside%ChargePresent - Nodo%ChargePast)
            endif
         end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
      endif
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !END OF CHARGE ADVANCING from n+1.0_RKIND_wires / 2 to n+3/2
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Inject the current at n+1 into the electring FDTD field previously calculated at n+1
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !actualizo bien el campo con los globales (hilos paralelos)
      !I sum up all the currents (in case of unshielded segment)
      !this is the correct approach for Parallel wires (later to be corrected for BerengerMTLN)
      !

      !voids the null_field vale for embedded and paralel segments (if requested at run time with -groundwires and without the -stableradholland)
      HWires%null_field = 0.0_RKIND_wires
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do n=1,HWires%NumCurrentSegments
         Segmento => HWires%CurrentSegment(n)
         if (.not.Segmento%IsShielded) then
!171216quitado            Segmento%Efield_wire2main_past = real(Segmento%Efield_wire2main,KIND=RKIND_wires)
#ifdef CompileWithThickWires
             if (wirethickness>1) then
                call Advance_Thick_Efield_wire2main(sgg,segmento,eps0,mu0,wirethickness)
             endif 
#endif                      
             if (wirethickness==1) then
                segmento%efield_wire2main=real(segmento%efield_wire2main,kind=rkind_wires) &
                    - segmento%cte5 * segmento%current
             endif
         endif
      end do

      !revoids the null_field value (unvoided above) for embedded and paralel segments (if requested at run time with -groundwires)
      HWires%null_field = 0.0_RKIND_wires
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !CURRENT ADVANCING from n+1 to n+2
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !

      if (trim(adjustl(wiresflavor))=='transition') then
         !recuerdo intensidades
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED)   private(Segmento)
#endif
         do is1 = 1,HWires%NumCurrentSegments
            Segmento => HWires%CurrentSegment(is1)
            Segmento%CurrentPast = Segmento%Current
         end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif

#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED)   private(Multiline,Segmento,Segmento2,Qplus,QMinus)
#endif
         do iw1 = 1,HWires%NumMultilines
            Multiline => HWires%Multilines(iw1)
            do is1 = 1,Multiline%NumParallel
               Segmento => Multiline%Segments(is1)%ptr
               Segmento%Current = 0.0_RKIND_wires
               do is2 = 1,Multiline%NumParallel
                  Segmento2 => Multiline%Segments(is2)%ptr
                  Qplus  = Segmento2%ChargePlus%ChargePresent
                  QMinus = Segmento2%ChargeMinus%ChargePresent
                  Segmento2 => Multiline%Segments(is2)%ptr
                  Segmento%Current = Segmento%Current                                         + &
                  Multiline%b1I(is1,is2)* Segmento2%CurrentPast            + &
                  Multiline%b2I(is1,is2)*(Segmento2%fractionPlus*Qplus-Segmento2%fractionMinus*QMinus)
                  if(.not.(Segmento%IsShielded.and.Segmento2%IsShielded)) then
                      !!!lo he descomentado a 300523 porque creo que estaba mal. ojo si algun dia se usa el flavor transition. pongo un stop para avisar
                     Segmento%Current = Segmento%Current + Multiline%b3I(is1,is2)*real(Segmento2%Efield_main2wire,KIND=RKIND_wires)
                     stop
                  end if
               end do
            end do
         end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
         !!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!
      elseif (trim(adjustl(wiresflavor))=='holland') then
         !!!!!!!!!!!!!!!!!!!
         !!!!!!MY FLAVOR!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED)   private(Segmento,Qplus,QMinus,jmed)
#endif
         do n=1,HWires%NumCurrentSegments
            Segmento => HWires%CurrentSegment(n)
            Segmento%CurrentPast = Segmento%Current !save this data only for the probes observation right time positioning
            !--->
            if (Segmento%IsPMC) then
               Segmento%Current=0.0_RKIND_wires
            else
               Qplus  = Segmento%ChargePlus%ChargePresent
               QMinus = Segmento%ChargeMinus%ChargePresent
               Segmento%qplus_qminus=Segmento%fractionPlus*Qplus-Segmento%fractionMinus*QMinus
               Segmento%Current=Segmento%cte1*Segmento%Current - Segmento%cte3*(Segmento%qplus_qminus)
               if (.not.Segmento%IsShielded) then
#ifdef CompileWithThickWires
                   if (wirethickness>1) then
                       call Advance_Thick_Efield_main2wire(sgg,segmento,eps0,mu0,wirethickness)
                   endif 
#endif                            
                   if (wirethickness==1) then
                     Segmento%Current = Segmento%Current + &
                         Segmento%cte2*real(Segmento%Efield_main2wire,KIND=RKIND_wires)
                   endif           
               endif
            endif
         end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
         !!!!!!!!!!!!
         !!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!
      else
         call WarnErrReport('wir0_BUGGYERROR: this flavor not permitted in sgg routine',.true.)
      endif !wiresflavor
      !igual en ambos
      if (.not.simu_devia) then             
          if (thereAreVsources) then
             !Feed the transparent voltage sources as appropriate charges in the current segments (Like in Edelvik's model)
             do n=1,HWires%NumCurrentSegments
                If  (HWires%CurrentSegment(n)%HasVsource) then
                   Segmento => HWires%CurrentSegment(n)
                   Vincid=evolucion(timei,Segmento%Vsource%Fichero%Samples, &
                   Segmento%Vsource%Fichero%DeltaSamples,Segmento%Vsource%Fichero%NumSamples)
                   !!!if (experimentalVideal) then
                   !!!    if ((.not.Segmento%ChargePlus%isPEC).and.Segmento%ChargeMinus%isPEC) then
                   !!!         Segmento%ChargePlus%ChargePresent = Vincid  / (Segmento%Lind * InvMu(Segmento%indexmed)*InvEps(Segmento%indexmed))
                   !!!    elseif ((.not.Segmento%ChargeMinus%isPEC).and.Segmento%ChargePlus%isPEC) then
                   !!!         Segmento%ChargeMinus%ChargePresent =  Vincid  / (Segmento%Lind * InvMu(Segmento%indexmed)*InvEps(Segmento%indexmed))
                   !!!    else
                   !!!        print *,'error en experimentalVideal 200621'
                   !!!    endif
                   !!!else
                       !lo de siempre. aniado lo anterior para ver lo de las fuentes duras !C=Q/V-> C=c^-2/L-> Q=V*C=V/(L c^2)
                       if (Segmento%Vsource%soft) then !fuentes blandas usuales 230323 
                           !fuentes de voltaje. para que sean de carga hay que comentar el Lind
                           Segmento%Current = Segmento%Current + &
                                 Segmento%cte3 * Vincid  / (Segmento%Lind * InvMu(Segmento%indexmed)*InvEps(Segmento%indexmed))
                           !I use the capacitance to find the incident charge
                           !assuming that the evolution file contains a voltage, not a charge
                       else !nuevas fuentes duras 230323
                           Iincid=Vincid !realmente se trata de forzar la corriente 230323
                           Segmento%Current =  Iincid
                       endif
                   !!!endif
                endif
             end do
          endif
      endif
      
!!!230323 comento lo que sigue porque ya manejo fuentes duras y blandas antes, pero esto funciono con jag rayos junio'15
      !!!!!!080615  uso una fuente dura de corriente correos jag simulacion rayos Junio'15
      !!!if (.not.simu_devia) then             
      !!!    if (thereAreIsources) then
      !!!       do n=1,HWires%NumChargeNodes
      !!!          if (HWires%ChargeNode(n)%exists) then
      !!!             If  (HWires%ChargeNode(n)%HasIsource) then
      !!!                Nodo => HWires%ChargeNode(n)
      !!!                Iincid=evolucion(timei,Nodo%Isource%Fichero%Samples, &
      !!!                Nodo%Isource%Fichero%DeltaSamples,Nodo%Isource%Fichero%NumSamples)
      !!!                Nodo%CurrentPlus_1%Current = Iincid
      !!!             endif
      !!!          endif
      !!!       end do
      !!!    endif
      !!!endif
      !!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !END OF CURRENT ADVANCING from n+1 to n+2
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!
!stochastic
#ifdef CompileWithStochastic
      if (stochastic.and.simu_devia) call inject_deviasources(layoutnumber,Hwires) !solo son los segmentos los afectados
#endif
 
!!!!!!!!!!!!!!!051223 no pongo a cero los already_YEEadvanced_byconformal_changedtoPECfield !No son cero, sino lo que se updatee mas lo que se inyecte
!!!machaca el campo que haya metido el wires

!!!#ifdef CompileWithOpenMP
!!!!$OMP PARALLEL DO DEFAULT(SHARED)  private(Nodo)
!!!#endif
!!!      do n=1,HWires%NumChargeNodes
!!!         Nodo => HWires%ChargeNode(n)
!!!!!!!140220 pon a PEC los viejos notouch=already_YEEadvanced_byconformal_changedtoPECfield que tenga conectado un nodo
!!!         !debe ser lo primero que se hace para overridear el call conformal_advance_E 
!!!         if (associated(nodo%already_YEEadvanced_byconformal_changedtoPECfield1)) then
!!!             nodo%already_YEEadvanced_byconformal_changedtoPECfield1=0.0_RKIND
!!!         endif
!!!         if (associated(nodo%already_YEEadvanced_byconformal_changedtoPECfield2)) then
!!!             nodo%already_YEEadvanced_byconformal_changedtoPECfield2=0.0_RKIND
!!!         endif
!!!         if (associated(nodo%already_YEEadvanced_byconformal_changedtoPECfield3)) then
!!!             nodo%already_YEEadvanced_byconformal_changedtoPECfield3=0.0_RKIND
!!!         endif
!!!         if (associated(nodo%already_YEEadvanced_byconformal_changedtoPECfield4)) then
!!!             nodo%already_YEEadvanced_byconformal_changedtoPECfield4=0.0_RKIND
!!!         endif
!!!         if (associated(nodo%already_YEEadvanced_byconformal_changedtoPECfield5)) then
!!!             nodo%already_YEEadvanced_byconformal_changedtoPECfield5=0.0_RKIND
!!!         endif
!!!         if (associated(nodo%already_YEEadvanced_byconformal_changedtoPECfield6)) then
!!!             nodo%already_YEEadvanced_byconformal_changedtoPECfield6=0.0_RKIND
!!!         endif
!!!      end do
!!!#ifdef CompileWithOpenMP
!!!!$OMP END PARALLEL DO
!!!#endif
!
      
   return

   end subroutine AdvanceWiresE

!!!
!!!!!!!!!!!!!!!!
!!!!!!


   subroutine AdvanceWiresH(sgg,timeinstant, layoutnumber,wiresflavor,simu_devia,stochastic,experimentalVideal,wirethickness,eps0,mu0)
                    
      REAL (KIND=RKIND), intent(IN)           ::  eps0,mu0
      integer, intent(IN) :: wirethickness
      logical :: simu_devia,stochastic,experimentalVideal
      type (SGGFDTDINFO), intent(IN)      ::  sgg

      integer (kind=4)  :: n,jmed,layoutnumber,iw1,is1,is2

      integer (kind=4), intent(IN)  ::  timeinstant
      REAL (KIND=RKIND_wires)   ::  Iplus,IMinus,Qplus,QMinus,timei
      REAL (KIND=RKIND_wires)   ::  Vincid,Iincid
      type (CurrentSegments), pointer  ::  Segmento, Segmento2
      type (ChargeNodes), pointer  ::  Nodo
      type (TMultiline), pointer                      ::  Multiline
      character(len=*), INTENT(in) :: wiresflavor
      timei = sgg%tiempo(timeinstant) 
   
!

      do n=1,HWires%NumCurrentSegments
         Segmento => HWires%CurrentSegment(n)
#ifdef CompileWithThickWires
         if (wirethickness>1) then
            call Advance_Thick_Hfield_wire2main(sgg,Segmento,eps0,mu0)
         endif 
#endif                      
        if (wirethickness==1) then
            continue
        endif
      end do

      
   return

   end subroutine AdvanceWiresH

!!!


!!!!!!!
!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !Advancing charge and current routine
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine AdvanceWiresEcrank(sgg,timeinstant,layoutnumber,wiresflavor,simu_devia,stochastic)
      type (SGGFDTDINFO), intent(IN)      ::  sgg
      logical :: simu_devia,stochastic

      !!!

      integer (kind=4)  :: n,jmed,layoutnumber,iw1,is1,is2

      integer (kind=4), intent(IN)  ::  timeinstant
      REAL (KIND=RKIND_wires)   ::  Iplus,IMinus,IplusPast,IMinusPast,source,timei
      REAL (KIND=RKIND_wires)   ::  Vincid,Iincid
      type (CurrentSegments), pointer  ::  Segmento , Segmento2
      type (ChargeNodes), pointer  ::  Nodo
      type (TMultiline), pointer                      ::  Multiline
      character(len=*), INTENT(in) :: wiresflavor
      REAL (KIND=RKIND_wires) , dimension(1:HWires%NumCurrentSegments)  ::  a,b,c,d,x
      
      timei = sgg%tiempo(timeinstant) 
      !!!
      iplus=-1.0; iminus=-1.0;
      
      !!!! deprecado en pscale  110219 

      

#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED)  private(Nodo)
#endif
      do n=1,HWires%NumChargeNodes
         Nodo => HWires%ChargeNode(n)
!!!!140220 pon a PEC los viejos notouch=already_YEEadvanced_byconformal_changedtoPECfield que tenga conectado un nodo
         !debe ser lo primero que se hace para overridear el call conformal_advance_E 
         if (associated(nodo%already_YEEadvanced_byconformal_changedtoPECfield1)) then
             nodo%already_YEEadvanced_byconformal_changedtoPECfield1=0.0_RKIND
         endif
         if (associated(nodo%already_YEEadvanced_byconformal_changedtoPECfield2)) then
             nodo%already_YEEadvanced_byconformal_changedtoPECfield2=0.0_RKIND
         endif
         if (associated(nodo%already_YEEadvanced_byconformal_changedtoPECfield3)) then
             nodo%already_YEEadvanced_byconformal_changedtoPECfield3=0.0_RKIND
         endif
         if (associated(nodo%already_YEEadvanced_byconformal_changedtoPECfield4)) then
             nodo%already_YEEadvanced_byconformal_changedtoPECfield4=0.0_RKIND
         endif
         if (associated(nodo%already_YEEadvanced_byconformal_changedtoPECfield5)) then
             nodo%already_YEEadvanced_byconformal_changedtoPECfield5=0.0_RKIND
         endif
         if (associated(nodo%already_YEEadvanced_byconformal_changedtoPECfield6)) then
             nodo%already_YEEadvanced_byconformal_changedtoPECfield6=0.0_RKIND
         endif
      end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Inject the current at n+1 into the electring FDTD field previously calculated at n+1
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do n=1,HWires%NumCurrentSegments
         Segmento => HWires%CurrentSegment(n)
!171216quitado         Segmento%Efield_wire2main_past = Segmento%Efield_wire2main
         Segmento%Efield_wire2main=real(Segmento%Efield_wire2main,KIND=RKIND_wires) - Segmento%cte5 * Segmento%Current
      end do


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !CURRENT ADVANCING from n+1 to n+2
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !

      do n=1,HWires%NumCurrentSegments
         Segmento => HWires%CurrentSegment(n)
         b(n)=segmento%diag
         c(n)=segmento%upperdiag
         a(n)=segmento%lowerdiag
         d(n)=        segmento%rightCU       * Segmento%Current
         d(n)= d(n) + segmento%rightCHminus  * Segmento%ChargeMinus%ChargePresent
         d(n)= d(n) + segmento%rightCHplus   * Segmento%ChargePlus%ChargePresent
         d(n)= d(n) + real(Segmento%Efield_wire2main,KIND=RKIND_wires)
         if (Segmento%ChargeMinus%NumCurrentMinus==1) d(n)= d(n) + segmento%rightCUminus  * Segmento%ChargeMinus%CurrentMinus_1%Current
         if (Segmento%ChargePlus%NumCurrentPlus==1  ) d(n)= d(n) + segmento%rightCUplus   * Segmento%ChargePlus%CurrentPlus_1%Current
         if (.not.simu_devia) then             
             If  (Segmento%HasVsource) then
                source=evolucion(timei,Segmento%Vsource%Fichero%Samples, &
                Segmento%Vsource%Fichero%DeltaSamples,Segmento%Vsource%Fichero%NumSamples)
                d(n)=d(n) + source
             endif
         endif
         
      end do
      n=HWires%NumCurrentSegments
      call solve_tridiag_wires(a,b,c,d,x,n)
      !  a - sub-diagonal (means it is the diagonal below the main diagonal)
      !  b - the main diagonal
      !  c - sup-diagonal (means it is the diagonal above the main diagonal)
      !  d - right part
      !  x - the answer
      !  n - number of equations

#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED)   private(Segmento,jmed)
#endif
      do n=1,HWires%NumCurrentSegments
         Segmento => HWires%CurrentSegment(n)
         Segmento%CurrentPast = Segmento%Current !save this data only for the probes observation right time positioning
         Segmento%Current=x(n)
      end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif



!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Next ADVANCE THE CHARGE from n-1.0_RKIND_wires / 2 to n+1/2 using the current known at n+1/2
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED)  private(Iplus,IMinus,Nodo)
#endif
      do n=1,HWires%NumChargeNodes
         Nodo => HWires%ChargeNode(n)
         if (nodo%exists) then
            Nodo%ChargePast=Nodo%ChargePresent
            If (Nodo%NumCurrentPlus==1) Iplus = Nodo%CurrentPlus_1%Current
            If (Nodo%NumCurrentPlus==1) IplusPast = Nodo%CurrentPlus_1%CurrentPast
            !
            If (Nodo%NumCurrentMinus==1) IMinus = Nodo%CurrentMinus_1%Current
            If (Nodo%NumCurrentMinus==1) IMinusPast = Nodo%CurrentMinus_1%CurrentPast
            !
            if ((Nodo%NumCurrentMinus == 1).and.(Nodo%NumCurrentPlus == 0)) Iplus = -Iminus
            if ((Nodo%NumCurrentMinus == 0).and.(Nodo%NumCurrentPlus == 1)) IMinus = -Iplus
            if ((Nodo%NumCurrentMinus == 1).and.(Nodo%NumCurrentPlus == 0)) IplusPast = -IminusPast
            if ((Nodo%NumCurrentMinus == 0).and.(Nodo%NumCurrentPlus == 1)) IMinusPast = -IplusPast

            Nodo%ChargePresent =  Nodo%CteProp*Nodo%ChargePast   - Nodo%CtePlain*((Iplus+IplusPast)/2.0_RKIND_wires-(IMinus+IminusPast)/2.0_RKIND_wires)

         endif
      end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !END OF CHARGE ADVANCING from n+1.0_RKIND_wires / 2 to n+3/2
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED)  private(Nodo)
#endif
      do n=1,HWires%NumChargeNodes
         Nodo => HWires%ChargeNode(n)
!!!!140220 pon a PEC los viejos notouch=already_YEEadvanced_byconformal_changedtoPECfield que tenga conectado un nodo
         !debe ser lo primero que se hace para overridear el call conformal_advance_E 
         if (associated(nodo%already_YEEadvanced_byconformal_changedtoPECfield1)) then
             nodo%already_YEEadvanced_byconformal_changedtoPECfield1=0.0_RKIND
         endif
         if (associated(nodo%already_YEEadvanced_byconformal_changedtoPECfield2)) then
             nodo%already_YEEadvanced_byconformal_changedtoPECfield2=0.0_RKIND
         endif
         if (associated(nodo%already_YEEadvanced_byconformal_changedtoPECfield3)) then
             nodo%already_YEEadvanced_byconformal_changedtoPECfield3=0.0_RKIND
         endif
         if (associated(nodo%already_YEEadvanced_byconformal_changedtoPECfield4)) then
             nodo%already_YEEadvanced_byconformal_changedtoPECfield4=0.0_RKIND
         endif
         if (associated(nodo%already_YEEadvanced_byconformal_changedtoPECfield5)) then
             nodo%already_YEEadvanced_byconformal_changedtoPECfield5=0.0_RKIND
         endif
         if (associated(nodo%already_YEEadvanced_byconformal_changedtoPECfield6)) then
             nodo%already_YEEadvanced_byconformal_changedtoPECfield6=0.0_RKIND
         endif
      end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !END OF CURRENT ADVANCING from n+1 to n+2
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!
      return
   end subroutine AdvanceWiresEcrank


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Function to interpolate the evolution files at the desired time
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   REAL (KIND=RKIND_wires)  function evolucion(t,evol,deltaevol,numus)
      integer (kind=4)  ::  numus
      integer (kind=8)  ::  nprev
      REAL (KIND=RKIND_wires) , dimension(0 : numus)  ::  evol
      REAL (KIND=RKIND_wires)   ::  deltaevol, t
      !
      nprev=int((t)/deltaevol)
      if ((nprev+1 > numus).OR.(NPREV+1 <= 0)) then !SI NPREV<0 ES PORQUE SE HA DESBORADO EL ENTERO !BUG MIGEL 130614
         evolucion=0.0_RKIND_wires !if running out of samples, assume they are null
      else
         evolucion=(evol(nprev+1)-evol(nprev))/deltaevol*((t)-nprev*deltaevol)+evol(nprev) !linear interpolation
      endif

      return
   end function evolucion


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !Routine to store the fields for resumed problems
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine StoreFieldsWires
      integer (kind=4)  ::  i1

      !store data for resuming
      do i1=1,HWires%NumChargeNodes
         write(14,err=634) HWires%ChargeNode(i1)%ChargePresent
         write(14,err=634) HWires%ChargeNode(i1)%ChargePast
      end do
      do i1=1,HWires%NumCurrentSegments
         write(14,err=634) HWires%CurrentSegment(i1)%Current
         write(14,err=634) HWires%CurrentSegment(i1)%qplus_qminus
         write(14,err=634) HWires%CurrentSegment(i1)%current_for_devia
         write(14,err=634) HWires%CurrentSegment(i1)%qplus_qminus_for_devia
         write(14,err=634) HWires%CurrentSegment(i1)%Efield_main2wire_for_devia
      end do
      !
#ifdef CompileWithMPI
      write (14,err=634) Hwires%NumNeededCurrentUpMPI,Hwires%NumNeededCurrentDownMPI
      do i1=1,Hwires%NumNeededCurrentUpMPI
         write (14,err=634) HWires%MPIUpNeededCurrentSegment(i1)%Current
      end do
      do i1=1,Hwires%NumNeededCurrentDownMPI
         write (14,err=634) HWires%MPIDownNeededCurrentSegment(i1)%Current
      end do
#endif

      goto 635
634   call print11(0,SEPARADOR//separador//separador)
      call print11(0,'WIRES: ERROR WRITING RESTARTING FIELDS. IGNORING AND CONTINUING')
      call print11(0,SEPARADOR//separador//separador)          
635   return
   end subroutine StoreFieldsWires

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !Routine to free-up memory upon termination
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!

   Subroutine DestroyWires(sgg)
      type (SGGFDTDINFO), intent(INOUT)      ::  sgg
      integer (kind=4)  ::  i

      !free up memory !ojo no se como hacerlo
      do i=1,sgg%NumMedia
         if (sgg%Med(i)%Is%ThinWire) then
            if (associated(sgg%Med(i)%wire(1)%Vsource)) deallocate (sgg%Med(i)%wire(1)%Vsource)
            if (associated(sgg%Med(i)%wire(1)%Isource)) deallocate (sgg%Med(i)%wire(1)%Isource)
            if (associated(sgg%Med(i)%wire)) deallocate (sgg%Med(i)%wire)
         endif
      end do

      if (associated(HWires%WireTipoMedio )) deallocate (HWires%WireTipoMedio )
      if (associated(HWires%CurrentSegment)) deallocate (HWires%CurrentSegment)
      if (associated(HWires%ChargeNode    )) deallocate (HWires%ChargeNode    )
#ifdef CompileWithMPI
      if (Hwires%NumNeededCurrentUpMPI>0)   deallocate (HWires%MPIUpNeededCurrentSegment)
      if (Hwires%NumNeededCurrentDownMPI>0) deallocate (HWires%MPIDownNeededCurrentSegment)
#endif
   end subroutine

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !Function to publish the private wire data (used in observation and in MPI)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!

   function GetHwires() result(r)
      type(Thinwires_t), pointer  ::  r
      r=>Hwires
      return
   end function






   !!!!!!!!!!!!!!!!!!!!Just for reporting
   !!! It also set the IsHeterogeneousJunction flag correctly (not needed for anything else than reporting)
   subroutine ReportWireJunctions(layoutnumber,size,therearewires,ZI,ZE,groundwires,strictOLD,verbose)

      logical :: paralelos,groundwires,therearewires,Terminal,IsHeterogeneousJunction,paraErr,strictOLD,verbose
      character(len=BUFSIZE) :: buff
      integer (kind=4)  ::  i1,j1,layoutnumber,zi,ze,ierr,size,indio
      integer (kind=4)  ::  mini=1000000000,minj=1000000000,mink=1000000000,maxi=-1000000000,maxj=-1000000000,maxk=-1000000000
      type (CurrentSegments), pointer  ::  org,fin
      character (len=3), dimension(1:3) :: DIR
      character (LEN=BUFSIZE) :: ig
      type (ChargeNodes), pointer :: nodo
      type :: nodosopentoair_t 
        integer (kind=4) :: i,j,k,indexnode
      end type
      
      type (nodosopentoair_t), allocatable, dimension(:) :: nodosopentoair


      dir(iEx)=' X '
      dir(iEy)=' Y '
      dir(iEz)=' Z '




#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif
      if (therearewires) then
         write (buff,'(a)')  '----------------------------------------------------------------'
         call WarnErrReport(buff)
      endif

      do i1=1,HWires%NumCurrentSegments
         if (HWires%CurrentSegment(i1)%i < mini) mini=HWires%CurrentSegment(i1)%i
         if (HWires%CurrentSegment(i1)%j < minj) minj=HWires%CurrentSegment(i1)%j
         if (HWires%CurrentSegment(i1)%k < mink) mink=HWires%CurrentSegment(i1)%k
         if (HWires%CurrentSegment(i1)%i > maxi) maxi=HWires%CurrentSegment(i1)%i
         if (HWires%CurrentSegment(i1)%j > maxj) maxj=HWires%CurrentSegment(i1)%j
         if (HWires%CurrentSegment(i1)%k > maxk) maxk=HWires%CurrentSegment(i1)%k
         !
         if ((HWires%CurrentSegment(i1)%ChargeMinus%NumCurrentPlus  > MaxNumCurrentMinusPlus) .or. &
         (HWires%CurrentSegment(i1)%ChargeMinus%NumCurrentMinus > MaxNumCurrentMinusPlus)) then
            write (buff,'(a,i2,a,i7,3i7)')  'wir3_BUGGYERROR: More than ', MaxNumCurrentMinusPlus,' plus/minus junctions at ', &
            HWires%CurrentSegment(i1)%origindex,HWires%CurrentSegment(i1)%i, &
            HWires%CurrentSegment(i1)%j,HWires%CurrentSegment(i1)%k
            call WarnErrReport(buff,.true.)
            write (buff,'(a)')  'Contact with  to increase this limit.'
            call WarnErrReport(buff,.true.)
         endif
      end do
      if ((maxi >=mini).and.(maxj>=minj).and.(maxk>=mink)) write (buff,'(a,6i12)')  'wir3_INFO: BBOX for Holland WIREs ', mini,minj,mink,maxi,maxj,maxk
      if (HWires%NumCurrentSegments /= 0) call WarnErrReport(buff)


      do i1=1,HWires%NumCurrentSegments
         org=>HWires%CurrentSegment(i1)
         do j1=i1+1,HWires%NumCurrentSegments
            fin=>HWires%CurrentSegment(j1)
            paralelos = (org%i == fin%i).and.(org%j == fin%j).and.(org%k == fin%k).and.(org%tipofield == fin%tipofield)
            if (paralelos) then
               if (org%indexmed /= fin%indexmed) then
                  write (buff,'(a,2i7,a,3i7)')  'wir3_INFO: Parallel segments from different wires (multiWIRE) ', &
                  org%origindex,fin%origindex,' at ',org%i,org%j,org%k

                  if ((org%k >= ZI).and.(org%k <= ZE).and.verbose) call WarnErrReport(buff)
               else
                  if (.not.strictOLD) then
                     write (buff,'(a,2i7,a,3i7)')  'wir3_BUGGYERROR: Parallel segments from the same WIRE. ', &
                     org%origindex,fin%origindex,' at ',org%i,org%j,org%k
                     call WarnErrReport(buff,.true.)
                  else
                     write (buff,'(a,2i7,a,3i7)')  'wir3_INFO: Parallel segments from the same wire (multiWIRE) ', &
                     org%origindex,fin%origindex,' at ',org%i,org%j,org%k

                     if ((org%k >= ZI).and.(org%k <= ZE).and.verbose) call WarnErrReport(buff)
                  endif
               endif
            endif
         end do
      end do
      !

#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif

      do i1=1,HWires%NumChargeNodes
         nodo => HWires%Chargenode(i1)
         if (nodo%exists) then
            if ((nodo%NumCurrentMinus == 1).and.(nodo%NumCurrentPlus == 1)) then
               !special cases of 1plus and 1minus
               IsHeterogeneousJunction = (nodo%CurrentMinus_1%indexmed /= nodo%CurrentPlus_1%indexmed )
            elseif (nodo%NumCurrentMinus + nodo%NumCurrentPlus == 1) then
               IsHeterogeneousJunction =.false.
            else
               IsHeterogeneousJunction =.true.
               if (nodo%NumCurrentMInus >= 2) &
               IsHeterogeneousJunction = IsHeterogeneousJunction .and.  &
               (nodo%CurrentMinus_1%indexmed /= nodo%CurrentMinus_2%indexmed )
               if (nodo%NumCurrentMInus >= 3) &
               IsHeterogeneousJunction = IsHeterogeneousJunction .and.  &
               (nodo%CurrentMinus_1%indexmed /= nodo%CurrentMinus_3%indexmed ).and. &
               (nodo%CurrentMinus_2%indexmed /= nodo%CurrentMinus_3%indexmed )
               if (nodo%NumCurrentMInus >= 4) &
               IsHeterogeneousJunction = IsHeterogeneousJunction .and.  &
               (nodo%CurrentMinus_1%indexmed /= nodo%CurrentMinus_4%indexmed ).and. &
               (nodo%CurrentMinus_2%indexmed /= nodo%CurrentMinus_4%indexmed ).and. &
               (nodo%CurrentMinus_3%indexmed /= nodo%CurrentMinus_4%indexmed )
               if (nodo%NumCurrentMInus >= 5) &
               IsHeterogeneousJunction = IsHeterogeneousJunction .and.  &
               (nodo%CurrentMinus_1%indexmed /= nodo%CurrentMinus_5%indexmed ).and. &
               (nodo%CurrentMinus_2%indexmed /= nodo%CurrentMinus_5%indexmed ).and. &
               (nodo%CurrentMinus_3%indexmed /= nodo%CurrentMinus_5%indexmed ).and. &
               (nodo%CurrentMinus_4%indexmed /= nodo%CurrentMinus_5%indexmed )
               if (nodo%NumCurrentMInus >= 6) &
               IsHeterogeneousJunction = IsHeterogeneousJunction .and.  &
               (nodo%CurrentMinus_1%indexmed /= nodo%CurrentMinus_6%indexmed ).and. &
               (nodo%CurrentMinus_2%indexmed /= nodo%CurrentMinus_6%indexmed ).and. &
               (nodo%CurrentMinus_3%indexmed /= nodo%CurrentMinus_6%indexmed ).and. &
               (nodo%CurrentMinus_4%indexmed /= nodo%CurrentMinus_6%indexmed ).and. &
               (nodo%CurrentMinus_5%indexmed /= nodo%CurrentMinus_6%indexmed )
               if (nodo%NumCurrentMInus >= 7) &
               IsHeterogeneousJunction = IsHeterogeneousJunction .and.  &
               (nodo%CurrentMinus_1%indexmed /= nodo%CurrentMinus_7%indexmed ).and. &
               (nodo%CurrentMinus_2%indexmed /= nodo%CurrentMinus_7%indexmed ).and. &
               (nodo%CurrentMinus_3%indexmed /= nodo%CurrentMinus_7%indexmed ).and. &
               (nodo%CurrentMinus_4%indexmed /= nodo%CurrentMinus_7%indexmed ).and. &
               (nodo%CurrentMinus_5%indexmed /= nodo%CurrentMinus_7%indexmed ).and. &
               (nodo%CurrentMinus_6%indexmed /= nodo%CurrentMinus_7%indexmed )
               if (nodo%NumCurrentMInus >= 8) &
               IsHeterogeneousJunction = IsHeterogeneousJunction .and.  &
               (nodo%CurrentMinus_1%indexmed /= nodo%CurrentMinus_8%indexmed ).and. &
               (nodo%CurrentMinus_2%indexmed /= nodo%CurrentMinus_8%indexmed ).and. &
               (nodo%CurrentMinus_3%indexmed /= nodo%CurrentMinus_8%indexmed ).and. &
               (nodo%CurrentMinus_4%indexmed /= nodo%CurrentMinus_8%indexmed ).and. &
               (nodo%CurrentMinus_5%indexmed /= nodo%CurrentMinus_8%indexmed ).and. &
               (nodo%CurrentMinus_6%indexmed /= nodo%CurrentMinus_8%indexmed ).and. &
               (nodo%CurrentMinus_7%indexmed /= nodo%CurrentMinus_8%indexmed )
               if (nodo%NumCurrentMInus >= 9) &
               IsHeterogeneousJunction = IsHeterogeneousJunction .and.  &
               (nodo%CurrentMinus_1%indexmed /= nodo%CurrentMinus_9%indexmed ).and. &
               (nodo%CurrentMinus_2%indexmed /= nodo%CurrentMinus_9%indexmed ).and. &
               (nodo%CurrentMinus_3%indexmed /= nodo%CurrentMinus_9%indexmed ).and. &
               (nodo%CurrentMinus_4%indexmed /= nodo%CurrentMinus_9%indexmed ).and. &
               (nodo%CurrentMinus_5%indexmed /= nodo%CurrentMinus_9%indexmed ).and. &
               (nodo%CurrentMinus_6%indexmed /= nodo%CurrentMinus_9%indexmed ).and. &
               (nodo%CurrentMinus_7%indexmed /= nodo%CurrentMinus_9%indexmed ).and. &
               (nodo%CurrentMinus_8%indexmed /= nodo%CurrentMinus_9%indexmed )
               !
               if (nodo%NumCurrentPlus >= 2) &
               IsHeterogeneousJunction =  IsHeterogeneousJunction .and.  &
               (nodo%CurrentPlus_1%indexmed /= nodo%CurrentPlus_2%indexmed )
               if (nodo%NumCurrentPlus >= 3) &
               IsHeterogeneousJunction = IsHeterogeneousJunction .and.  &
               (nodo%CurrentPlus_1%indexmed /= nodo%CurrentPlus_3%indexmed ).and. &
               (nodo%CurrentPlus_2%indexmed /= nodo%CurrentPlus_3%indexmed )
               if (nodo%NumCurrentPlus >= 4) &
               IsHeterogeneousJunction = IsHeterogeneousJunction .and.  &
               (nodo%CurrentPlus_1%indexmed /= nodo%CurrentPlus_4%indexmed ).and. &
               (nodo%CurrentPlus_2%indexmed /= nodo%CurrentPlus_4%indexmed ).and. &
               (nodo%CurrentPlus_3%indexmed /= nodo%CurrentPlus_4%indexmed )
               if (nodo%NumCurrentPlus >= 5) &
               IsHeterogeneousJunction = IsHeterogeneousJunction .and.  &
               (nodo%CurrentPlus_1%indexmed /= nodo%CurrentPlus_5%indexmed ).and. &
               (nodo%CurrentPlus_2%indexmed /= nodo%CurrentPlus_5%indexmed ).and. &
               (nodo%CurrentPlus_3%indexmed /= nodo%CurrentPlus_5%indexmed ).and. &
               (nodo%CurrentPlus_4%indexmed /= nodo%CurrentPlus_5%indexmed )
               if (nodo%NumCurrentPlus >= 6) &
               IsHeterogeneousJunction = IsHeterogeneousJunction .and.  &
               (nodo%CurrentPlus_1%indexmed /= nodo%CurrentPlus_6%indexmed ).and. &
               (nodo%CurrentPlus_2%indexmed /= nodo%CurrentPlus_6%indexmed ).and. &
               (nodo%CurrentPlus_3%indexmed /= nodo%CurrentPlus_6%indexmed ).and. &
               (nodo%CurrentPlus_4%indexmed /= nodo%CurrentPlus_6%indexmed ).and. &
               (nodo%CurrentPlus_5%indexmed /= nodo%CurrentPlus_6%indexmed )
               if (nodo%NumCurrentPlus >= 7) &
               IsHeterogeneousJunction = IsHeterogeneousJunction .and.  &
               (nodo%CurrentPlus_1%indexmed /= nodo%CurrentPlus_7%indexmed ).and. &
               (nodo%CurrentPlus_2%indexmed /= nodo%CurrentPlus_7%indexmed ).and. &
               (nodo%CurrentPlus_3%indexmed /= nodo%CurrentPlus_7%indexmed ).and. &
               (nodo%CurrentPlus_4%indexmed /= nodo%CurrentPlus_7%indexmed ).and. &
               (nodo%CurrentPlus_5%indexmed /= nodo%CurrentPlus_7%indexmed ).and. &
               (nodo%CurrentPlus_6%indexmed /= nodo%CurrentPlus_7%indexmed )
               if (nodo%NumCurrentPlus >= 8) &
               IsHeterogeneousJunction = IsHeterogeneousJunction .and.  &
               (nodo%CurrentPlus_1%indexmed /= nodo%CurrentPlus_8%indexmed ).and. &
               (nodo%CurrentPlus_2%indexmed /= nodo%CurrentPlus_8%indexmed ).and. &
               (nodo%CurrentPlus_3%indexmed /= nodo%CurrentPlus_8%indexmed ).and. &
               (nodo%CurrentPlus_4%indexmed /= nodo%CurrentPlus_8%indexmed ).and. &
               (nodo%CurrentPlus_5%indexmed /= nodo%CurrentPlus_8%indexmed ).and. &
               (nodo%CurrentPlus_6%indexmed /= nodo%CurrentPlus_8%indexmed ).and. &
               (nodo%CurrentPlus_7%indexmed /= nodo%CurrentPlus_8%indexmed )
               if (nodo%NumCurrentPlus >= 9) &
               IsHeterogeneousJunction = IsHeterogeneousJunction .and.  &
               (nodo%CurrentPlus_1%indexmed /= nodo%CurrentPlus_9%indexmed ).and. &
               (nodo%CurrentPlus_2%indexmed /= nodo%CurrentPlus_9%indexmed ).and. &
               (nodo%CurrentPlus_3%indexmed /= nodo%CurrentPlus_9%indexmed ).and. &
               (nodo%CurrentPlus_4%indexmed /= nodo%CurrentPlus_9%indexmed ).and. &
               (nodo%CurrentPlus_5%indexmed /= nodo%CurrentPlus_9%indexmed ).and. &
               (nodo%CurrentPlus_6%indexmed /= nodo%CurrentPlus_9%indexmed ).and. &
               (nodo%CurrentPlus_7%indexmed /= nodo%CurrentPlus_9%indexmed ).and. &
               (nodo%CurrentPlus_8%indexmed /= nodo%CurrentPlus_9%indexmed )
            endif

            if ((IsHeterogeneousJunction       .and. nodo%IsHeterogeneousJunction).or. &
            (.not.IsHeterogeneousJunction).and.(.not.nodo%IsHeterogeneousJunction)) then
               continue
            else
               write (buff,'(a,i7,3i7,a,2i3,a)')  'wir3_BUGGYERROR: Heterogeneous Junctions mismatch ', &
               nodo%indexnode, nodo%i,nodo%j,nodo%k, &
               ' (',nodo%numcurrentminus,nodo%numcurrentplus,'). '
               if ((nodo%k >=  ZI).and.(nodo%k<=ZE)) call WarnErrReport(buff,.true.)
            endif

         endif !del if exist
      end do


#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif


      !!!!!!!!!!!!!!!!
      do i1=1,HWires%NumChargeNodes
         nodo => Hwires%ChargeNode(i1)
         if (nodo%exists) then
            if (nodo%isLossy.or.nodo%isPEC) then
               if (nodo%isLossy) ig=' to Lossy'
               if (nodo%isLossy)   ig=' to Lossy'
               if (nodo%ispec)   ig=' to PEC'
                  if ((nodo%Is_LeftEnd).or.(nodo%Is_RightEnd)) then
                     write (buff,'(a,i7,3i7,a,2i3,a,3e12.2e3)')  'wir3_INFO: Terminal (LeftEnd/RightEnd) node directly GROUNDED  ',nodo%indexnode, &
                     nodo%i,nodo%j,nodo%k, &
                     ' (',nodo%numcurrentminus,nodo%numcurrentplus,')'//trim(adjustl(ig)),Nodo%CteProp,Nodo%CtePlain
                     if ((nodo%k >=  ZI).and.(nodo%k <= ZE).and.verbose) call WarnErrReport(buff)
                  elseif ((nodo%NumCurrentPlus + nodo%NumCurrentMinus < 2)) then
                     write (buff,'(a,i7,3i7,a,2i3,a,3e12.2e3)')  'wir3_WARNING: Terminal (other) node  direcly GROUNDED  ',nodo%indexnode, &
                     nodo%i,nodo%j,nodo%k, &
                     ' (',nodo%numcurrentminus,nodo%numcurrentplus,')'//trim(adjustl(ig)),Nodo%CteProp,Nodo%CtePlain
                     if ((nodo%k >=  ZI).and.(nodo%k <= ZE)) call WarnErrReport(buff)
                  else
                     if (groundwires) then
                        write (buff,'(a,i7,3i7,a,2i3,a,3e12.2e3)')  'wir3_INFO: NON-terminal node directly GROUNDED (-groundWIREs)', &
                        nodo%indexnode, nodo%i,nodo%j,nodo%k, &
                        ' (',nodo%numcurrentminus,nodo%numcurrentplus,')'//trim(adjustl(ig)),Nodo%CteProp,Nodo%CtePlain
                     if ((nodo%k >=  ZI).and.(nodo%k <= ZE).and.verbose) call WarnErrReport(buff)
                     else
                        write (buff,'(a,i7,3i7,a,2i3,a,3e12.2e3)')  'wir3_WARNING: NON-terminal node directly GROUNDED ', &
                        nodo%indexnode, nodo%i,nodo%j,nodo%k, &
                        ' (',nodo%numcurrentminus,nodo%numcurrentplus,')'//trim(adjustl(ig)),Nodo%CteProp,Nodo%CtePlain
                     if ((nodo%k >=  ZI).and.(nodo%k <= ZE)) call WarnErrReport(buff)
                     endif
                  endif
            endif
         endif
      end do


      do i1=1,HWires%NumChargeNodes
         nodo => Hwires%ChargeNode(i1)
         if (nodo%exists) then
            if (nodo%IsHeterogeneousJunction) then
               !
               write (buff,'(a,i7,3i7,a,2i3,a)')  'wir3_INFO: MultiWIRE Junction made at node ',nodo%indexnode, &
               nodo%i,nodo%j,nodo%k, &
               ' (',nodo%numcurrentminus,nodo%numcurrentplus,')'
               if ((nodo%k >=  ZI).and.(nodo%k<=ZE).and.verbose) call WarnErrReport(buff)
               Terminal = .false. !7/2/14 esti estaba mal
               if (nodo%NumCurrentMinus >=1 ) then
                  Terminal = Terminal .or. &
                  (nodo%CurrentMinus_1%Is_LeftEnd.or.nodo%CurrentMinus_1%Is_RightEnd.or.nodo%CurrentMinus_1%IsEnd_norLeft_norRight)
                  write (buff,'(a,i7,3i7,a,i7)')  '         Segment -1: ', &
                  nodo%CurrentMinus_1%origindex, &
                  nodo%CurrentMinus_1%i, &
                  nodo%CurrentMinus_1%j, &
                  nodo%CurrentMinus_1%k, &
                  dir(nodo%CurrentMinus_1%tipofield), &
                  nodo%CurrentMinus_1%indexmed
                  if ((nodo%k >=  ZI).and.(nodo%k<=ZE).and.verbose) call WarnErrReport(buff)
               endif
               if (nodo%NumCurrentMinus >= 2) then
                  Terminal = Terminal .or. &
                  (nodo%CurrentMinus_2%Is_LeftEnd.or.nodo%CurrentMinus_2%Is_RightEnd.or.nodo%CurrentMinus_2%IsEnd_norLeft_norRight)
                  write (buff,'(a,i7,3i7,a,i7)')  '         Segment -2: ',  &
                  nodo%CurrentMinus_2%origindex, &
                  nodo%CurrentMinus_2%i, &
                  nodo%CurrentMinus_2%j, &
                  nodo%CurrentMinus_2%k, &
                  dir(nodo%CurrentMinus_2%tipofield), &
                  nodo%CurrentMinus_2%indexmed
                  if ((nodo%k >=  ZI).and.(nodo%k<=ZE).and.verbose) call WarnErrReport(buff)
               endif
               if (nodo%NumCurrentMinus >= 3) then
                  Terminal = Terminal .or. &
                  (nodo%CurrentMinus_3%Is_LeftEnd.or.nodo%CurrentMinus_3%Is_RightEnd.or.nodo%CurrentMinus_3%IsEnd_norLeft_norRight)
                  write (buff,'(a,i7,3i7,a,i7)')  '         Segment -3: ',  &
                  nodo%CurrentMinus_3%origindex, &
                  nodo%CurrentMinus_3%i, &
                  nodo%CurrentMinus_3%j, &
                  nodo%CurrentMinus_3%k, &
                  dir(nodo%CurrentMinus_3%tipofield), &
                  nodo%CurrentMinus_3%indexmed
                  if ((nodo%k >=  ZI).and.(nodo%k<=ZE).and.verbose) call WarnErrReport(buff)
               endif
               if (nodo%NumCurrentMinus >= 4) then
                  Terminal = Terminal .or. &
                  (nodo%CurrentMinus_4%Is_LeftEnd.or.nodo%CurrentMinus_4%Is_RightEnd.or.nodo%CurrentMinus_4%IsEnd_norLeft_norRight)
                  write (buff,'(a,i7,3i7,a,i7)')  '         Segment -4: ',  &
                  nodo%CurrentMinus_4%origindex, &
                  nodo%CurrentMinus_4%i, &
                  nodo%CurrentMinus_4%j, &
                  nodo%CurrentMinus_4%k, &
                  dir(nodo%CurrentMinus_4%tipofield), &
                  nodo%CurrentMinus_4%indexmed
                  if ((nodo%k >=  ZI).and.(nodo%k<=ZE).and.verbose) call WarnErrReport(buff)
               endif
               if (nodo%NumCurrentMinus >= 5) then
                  Terminal = Terminal .or. &
                  (nodo%CurrentMinus_5%Is_LeftEnd.or.nodo%CurrentMinus_5%Is_RightEnd.or.nodo%CurrentMinus_5%IsEnd_norLeft_norRight)
                  write (buff,'(a,i7,3i7,a,i7)')  '         Segment -5: ',  &
                  nodo%CurrentMinus_5%origindex, &
                  nodo%CurrentMinus_5%i, &
                  nodo%CurrentMinus_5%j, &
                  nodo%CurrentMinus_5%k, &
                  dir(nodo%CurrentMinus_5%tipofield), &
                  nodo%CurrentMinus_5%indexmed
                  if ((nodo%k >=  ZI).and.(nodo%k<=ZE).and.verbose) call WarnErrReport(buff)
               endif
               if (nodo%NumCurrentMinus >= 6) then
                  Terminal = Terminal .or. &
                  (nodo%CurrentMinus_6%Is_LeftEnd.or.nodo%CurrentMinus_6%Is_RightEnd.or.nodo%CurrentMinus_6%IsEnd_norLeft_norRight)
                  write (buff,'(a,i7,3i7,a,i7)')  '         Segment -6: ',  &
                  nodo%CurrentMinus_6%origindex, &
                  nodo%CurrentMinus_6%i, &
                  nodo%CurrentMinus_6%j, &
                  nodo%CurrentMinus_6%k, &
                  dir(nodo%CurrentMinus_6%tipofield), &
                  nodo%CurrentMinus_6%indexmed
                  if ((nodo%k >=  ZI).and.(nodo%k<=ZE).and.verbose) call WarnErrReport(buff)
               endif
               if (nodo%NumCurrentMinus >= 7) then
                  Terminal = Terminal .or. &
                  (nodo%CurrentMinus_7%Is_LeftEnd.or.nodo%CurrentMinus_7%Is_RightEnd.or.nodo%CurrentMinus_7%IsEnd_norLeft_norRight)
                  write (buff,'(a,i7,3i7,a,i7)')  '         Segment -7: ',  &
                  nodo%CurrentMinus_7%origindex, &
                  nodo%CurrentMinus_7%i, &
                  nodo%CurrentMinus_7%j, &
                  nodo%CurrentMinus_7%k, &
                  dir(nodo%CurrentMinus_7%tipofield), &
                  nodo%CurrentMinus_7%indexmed
                  if ((nodo%k >=  ZI).and.(nodo%k<=ZE).and.verbose) call WarnErrReport(buff)
               endif
               if (nodo%NumCurrentMinus >= 8) then
                  Terminal = Terminal .or. &
                  (nodo%CurrentMinus_8%Is_LeftEnd.or.nodo%CurrentMinus_8%Is_RightEnd.or.nodo%CurrentMinus_8%IsEnd_norLeft_norRight)
                  write (buff,'(a,i7,3i7,a,i7)')  '         Segment -8: ',  &
                  nodo%CurrentMinus_8%origindex, &
                  nodo%CurrentMinus_8%i, &
                  nodo%CurrentMinus_8%j, &
                  nodo%CurrentMinus_8%k, &
                  dir(nodo%CurrentMinus_8%tipofield), &
                  nodo%CurrentMinus_8%indexmed
                  if ((nodo%k >=  ZI).and.(nodo%k<=ZE).and.verbose) call WarnErrReport(buff)
               endif
               if (nodo%NumCurrentMinus >= 9) then
                  Terminal = Terminal .or. &
                  (nodo%CurrentMinus_9%Is_LeftEnd.or.nodo%CurrentMinus_9%Is_RightEnd.or.nodo%CurrentMinus_9%IsEnd_norLeft_norRight)
                  write (buff,'(a,i7,3i7,a,i7)')  '         Segment -9: ',  &
                  nodo%CurrentMinus_9%origindex, &
                  nodo%CurrentMinus_9%i, &
                  nodo%CurrentMinus_9%j, &
                  nodo%CurrentMinus_9%k, &
                  dir(nodo%CurrentMinus_9%tipofield), &
                  nodo%CurrentMinus_9%indexmed
                  if ((nodo%k >=  ZI).and.(nodo%k<=ZE)) call WarnErrReport(buff)
               endif
               !
               if (nodo%NumCurrentplus >= 1) then
                  Terminal = Terminal .or. &
                  (nodo%CurrentPlus_1%Is_LeftEnd.or.nodo%CurrentPlus_1%Is_RightEnd.or.nodo%CurrentPlus_1%IsEnd_norLeft_norRight)
                  write (buff,'(a,i7,3i7,a,i7)')  '         Segment +1: ',  &
                  nodo%Currentplus_1%origindex, &
                  nodo%Currentplus_1%i, &
                  nodo%Currentplus_1%j, &
                  nodo%Currentplus_1%k, &
                  dir(nodo%Currentplus_1%tipofield), &
                  nodo%Currentplus_1%indexmed
                  if ((nodo%k >=  ZI).and.(nodo%k<=ZE).and.verbose) call WarnErrReport(buff)
               endif
               if (nodo%NumCurrentplus >= 2) then
                  Terminal = Terminal .or. &
                  (nodo%CurrentPlus_2%Is_LeftEnd.or.nodo%CurrentPlus_2%Is_RightEnd.or.nodo%CurrentPlus_2%IsEnd_norLeft_norRight)
                  write (buff,'(a,i7,3i7,a,i7)')  '         Segment +2: ',  &
                  nodo%Currentplus_2%origindex, &
                  nodo%Currentplus_2%i, &
                  nodo%Currentplus_2%j, &
                  nodo%Currentplus_2%k, &
                  dir(nodo%Currentplus_2%tipofield), &
                  nodo%Currentplus_2%indexmed
                  if ((nodo%k >=  ZI).and.(nodo%k<=ZE).and.verbose) call WarnErrReport(buff)
               endif
               if (nodo%NumCurrentplus >= 3) then
                  Terminal = Terminal .or. &
                  (nodo%CurrentPlus_3%Is_LeftEnd.or.nodo%CurrentPlus_3%Is_RightEnd.or.nodo%CurrentPlus_3%IsEnd_norLeft_norRight)
                  write (buff,'(a,i7,3i7,a,i7)')  '         Segment +3: ',  &
                  nodo%Currentplus_3%origindex, &
                  nodo%Currentplus_3%i, &
                  nodo%Currentplus_3%j, &
                  nodo%Currentplus_3%k, &
                  dir(nodo%Currentplus_3%tipofield), &
                  nodo%Currentplus_3%indexmed
                  if ((nodo%k >=  ZI).and.(nodo%k<=ZE)) call WarnErrReport(buff)
               endif
               if (nodo%NumCurrentplus >= 4) then
                  Terminal = Terminal .or. &
                  (nodo%CurrentPlus_4%Is_LeftEnd.or.nodo%CurrentPlus_4%Is_RightEnd.or.nodo%CurrentPlus_4%IsEnd_norLeft_norRight)
                  write (buff,'(a,i7,3i7,a,i7)')  '         Segment +4: ',  &
                  nodo%Currentplus_4%origindex, &
                  nodo%Currentplus_4%i, &
                  nodo%Currentplus_4%j, &
                  nodo%Currentplus_4%k, &
                  dir(nodo%Currentplus_4%tipofield), &
                  nodo%Currentplus_4%indexmed
                  if ((nodo%k >=  ZI).and.(nodo%k<=ZE).and.verbose) call WarnErrReport(buff)
               endif
               if (nodo%NumCurrentplus >= 5) then
                  Terminal = Terminal .or. &
                  (nodo%CurrentPlus_5%Is_LeftEnd.or.nodo%CurrentPlus_5%Is_RightEnd.or.nodo%CurrentPlus_5%IsEnd_norLeft_norRight)
                  write (buff,'(a,i7,3i7,a,i7)')  '         Segment +5: ',  &
                  nodo%Currentplus_5%origindex, &
                  nodo%Currentplus_5%i, &
                  nodo%Currentplus_5%j, &
                  nodo%Currentplus_5%k, &
                  dir(nodo%Currentplus_5%tipofield), &
                  nodo%Currentplus_5%indexmed
                  if ((nodo%k >=  ZI).and.(nodo%k<=ZE).and.verbose) call WarnErrReport(buff)
               endif
               if (nodo%NumCurrentplus >= 6) then
                  Terminal = Terminal .or. &
                  (nodo%CurrentPlus_6%Is_LeftEnd.or.nodo%CurrentPlus_6%Is_RightEnd.or.nodo%CurrentPlus_6%IsEnd_norLeft_norRight)
                  write (buff,'(a,i7,3i7,a,i7)')  '         Segment +6: ',  &
                  nodo%Currentplus_6%origindex, &
                  nodo%Currentplus_6%i, &
                  nodo%Currentplus_6%j, &
                  nodo%Currentplus_6%k, &
                  dir(nodo%Currentplus_6%tipofield), &
                  nodo%Currentplus_6%indexmed
                  if ((nodo%k >=  ZI).and.(nodo%k<=ZE).and.verbose) call WarnErrReport(buff)
               endif
               if (nodo%NumCurrentplus >= 7) then
                  Terminal = Terminal .or. &
                  (nodo%CurrentPlus_7%Is_LeftEnd.or.nodo%CurrentPlus_7%Is_RightEnd.or.nodo%CurrentPlus_7%IsEnd_norLeft_norRight)
                  write (buff,'(a,i7,3i7,a,i7)')  '         Segment +7: ',  &
                  nodo%Currentplus_7%origindex, &
                  nodo%Currentplus_7%i, &
                  nodo%Currentplus_7%j, &
                  nodo%Currentplus_7%k, &
                  dir(nodo%Currentplus_7%tipofield), &
                  nodo%Currentplus_7%indexmed
                  if ((nodo%k >=  ZI).and.(nodo%k<=ZE).and.verbose) call WarnErrReport(buff)
               endif
               if (nodo%NumCurrentplus >= 8) then
                  Terminal = Terminal .or. &
                  (nodo%CurrentPlus_8%Is_LeftEnd.or.nodo%CurrentPlus_8%Is_RightEnd.or.nodo%CurrentPlus_8%IsEnd_norLeft_norRight)
                  write (buff,'(a,i7,3i7,a,i7)')  '         Segment +8: ',  &
                  nodo%Currentplus_8%origindex, &
                  nodo%Currentplus_8%i, &
                  nodo%Currentplus_8%j, &
                  nodo%Currentplus_8%k, &
                  dir(nodo%Currentplus_8%tipofield), &
                  nodo%Currentplus_8%indexmed
                  if ((nodo%k >=  ZI).and.(nodo%k<=ZE).and.verbose) call WarnErrReport(buff)
               endif
               if (nodo%NumCurrentplus >= 9) then
                  Terminal = Terminal .or. &
                  (nodo%CurrentPlus_9%Is_LeftEnd.or.nodo%CurrentPlus_9%Is_RightEnd.or.nodo%CurrentPlus_9%IsEnd_norLeft_norRight)
                  write (buff,'(a,i7,3i7,a,i7)')  '         Segment +9: ',  &
                  nodo%Currentplus_9%origindex, &
                  nodo%Currentplus_9%i, &
                  nodo%Currentplus_9%j, &
                  nodo%Currentplus_9%k, &
                  dir(nodo%Currentplus_9%tipofield), &
                  nodo%Currentplus_9%indexmed
                  if ((nodo%k >=  ZI).and.(nodo%k<=ZE).and.verbose) call WarnErrReport(buff)
               endif
               if (.not.Terminal) then
                  write (buff,'(a)')  'wir3_BUGGYERROR: Some of the segments are not terminal'
                  call WarnErrReport(buff,.true.)
               endif

            endif
         endif
      end do

      !!!!!!!!!!!!!!!!

      do i1=1,HWires%NumChargeNodes
         nodo => Hwires%ChargeNode(i1)
         if (nodo%exists) then
            if ((.not.nodo%IsHeterogeneousJunction).and. &
            (nodo%NumCurrentMinus + nodo%NumCurrentPlus >= 3)) then !homogeneo y verdadera union
               if (.not.strictOLD) then
                  write (buff,'(a)')  'wir3_INFO: Intra-WIRE Junction made between: '
                  paraerr=.false.
                  if ((nodo%k >=  ZI).and.(nodo%k<=ZE).and.verbose) call WarnErrReport(buff,paraerr)
               else
                  write (buff,'(a)')  'wir3_ERROR: Intra-WIRE Junction of more than 2 segments is forbidden: '
                  paraerr=.true.
                  if ((nodo%k >=  ZI).and.(nodo%k<=ZE)) call WarnErrReport(buff,paraerr)
               endif
               write (buff,'(a,3i7)')  '         At node: ', &
               nodo%i, &
               nodo%j, &
               nodo%k
               if ((nodo%k >=  ZI).and.(nodo%k<=ZE)) call WarnErrReport(buff,paraerr)
               if (nodo%NumCurrentMinus >=1 ) then
                  write (buff,'(a,i7,3i7,a,i7)')  '         Segment -1: ', &
                  nodo%CurrentMinus_1%origindex, &
                  nodo%CurrentMinus_1%i, &
                  nodo%CurrentMinus_1%j, &
                  nodo%CurrentMinus_1%k, &
                  dir(nodo%CurrentMinus_1%tipofield), &
                  nodo%CurrentMinus_1%indexmed
                  if ((nodo%k >=  ZI).and.(nodo%k<=ZE)) call WarnErrReport(buff,paraerr)
               endif
               if (nodo%NumCurrentMinus >= 2) then
                  write (buff,'(a,i7,3i7,a,i7)')  '         Segment -2: ',  &
                  nodo%CurrentMinus_2%origindex, &
                  nodo%CurrentMinus_2%i, &
                  nodo%CurrentMinus_2%j, &
                  nodo%CurrentMinus_2%k, &
                  dir(nodo%CurrentMinus_2%tipofield), &
                  nodo%CurrentMinus_2%indexmed
                  if ((nodo%k >=  ZI).and.(nodo%k<=ZE)) call WarnErrReport(buff,paraerr)
               endif
               if (nodo%NumCurrentMinus >= 3) then
                  write (buff,'(a,i7,3i7,a,i7)')  '         Segment -3: ',  &
                  nodo%CurrentMinus_3%origindex, &
                  nodo%CurrentMinus_3%i, &
                  nodo%CurrentMinus_3%j, &
                  nodo%CurrentMinus_3%k, &
                  dir(nodo%CurrentMinus_3%tipofield), &
                  nodo%CurrentMinus_3%indexmed
                  if ((nodo%k >=  ZI).and.(nodo%k<=ZE)) call WarnErrReport(buff,paraerr)
               endif
               if (nodo%NumCurrentMinus >= 4) then
                  write (buff,'(a,i7,3i7,a,i7)')  '         Segment -4: ',  &
                  nodo%CurrentMinus_4%origindex, &
                  nodo%CurrentMinus_4%i, &
                  nodo%CurrentMinus_4%j, &
                  nodo%CurrentMinus_4%k, &
                  dir(nodo%CurrentMinus_4%tipofield), &
                  nodo%CurrentMinus_4%indexmed
                  if ((nodo%k >=  ZI).and.(nodo%k<=ZE)) call WarnErrReport(buff,paraerr)
               endif
               if (nodo%NumCurrentMinus >= 5) then
                  write (buff,'(a,i7,3i7,a,i7)')  '         Segment -5: ',  &
                  nodo%CurrentMinus_5%origindex, &
                  nodo%CurrentMinus_5%i, &
                  nodo%CurrentMinus_5%j, &
                  nodo%CurrentMinus_5%k, &
                  dir(nodo%CurrentMinus_5%tipofield), &
                  nodo%CurrentMinus_5%indexmed
                  if ((nodo%k >=  ZI).and.(nodo%k<=ZE)) call WarnErrReport(buff,paraerr)
               endif
               if (nodo%NumCurrentMinus >= 6) then
                  write (buff,'(a,i7,3i7,a,i7)')  '         Segment -6: ',  &
                  nodo%CurrentMinus_6%origindex, &
                  nodo%CurrentMinus_6%i, &
                  nodo%CurrentMinus_6%j, &
                  nodo%CurrentMinus_6%k, &
                  dir(nodo%CurrentMinus_6%tipofield), &
                  nodo%CurrentMinus_6%indexmed
                  if ((nodo%k >=  ZI).and.(nodo%k<=ZE)) call WarnErrReport(buff,paraerr)
               endif
               if (nodo%NumCurrentMinus >= 7) then
                  write (buff,'(a,i7,3i7,a,i7)')  '         Segment -7: ',  &
                  nodo%CurrentMinus_7%origindex, &
                  nodo%CurrentMinus_7%i, &
                  nodo%CurrentMinus_7%j, &
                  nodo%CurrentMinus_7%k, &
                  dir(nodo%CurrentMinus_7%tipofield), &
                  nodo%CurrentMinus_7%indexmed
                  if ((nodo%k >=  ZI).and.(nodo%k<=ZE)) call WarnErrReport(buff,paraerr)
               endif
               if (nodo%NumCurrentMinus >= 8) then
                  write (buff,'(a,i7,3i7,a,i7)')  '         Segment -8: ',  &
                  nodo%CurrentMinus_8%origindex, &
                  nodo%CurrentMinus_8%i, &
                  nodo%CurrentMinus_8%j, &
                  nodo%CurrentMinus_8%k, &
                  dir(nodo%CurrentMinus_8%tipofield), &
                  nodo%CurrentMinus_8%indexmed
                  if ((nodo%k >=  ZI).and.(nodo%k<=ZE)) call WarnErrReport(buff,paraerr)
               endif
               if (nodo%NumCurrentMinus >= 9) then
                  write (buff,'(a,i7,3i7,a,i7)')  '         Segment -9: ',  &
                  nodo%CurrentMinus_9%origindex, &
                  nodo%CurrentMinus_9%i, &
                  nodo%CurrentMinus_9%j, &
                  nodo%CurrentMinus_9%k, &
                  dir(nodo%CurrentMinus_9%tipofield), &
                  nodo%CurrentMinus_9%indexmed
                  if ((nodo%k >=  ZI).and.(nodo%k<=ZE)) call WarnErrReport(buff,paraerr)
               endif
               !
               if (nodo%NumCurrentplus >= 1) then
                  write (buff,'(a,i7,3i7,a,i7)')  '         Segment +1: ',  &
                  nodo%Currentplus_1%origindex, &
                  nodo%Currentplus_1%i, &
                  nodo%Currentplus_1%j, &
                  nodo%Currentplus_1%k, &
                  dir(nodo%Currentplus_1%tipofield), &
                  nodo%Currentplus_1%indexmed
                  if ((nodo%k >=  ZI).and.(nodo%k<=ZE)) call WarnErrReport(buff,paraerr)
               endif
               if (nodo%NumCurrentplus >= 2) then
                  write (buff,'(a,i7,3i7,a,i7)')  '         Segment +2: ',  &
                  nodo%Currentplus_2%origindex, &
                  nodo%Currentplus_2%i, &
                  nodo%Currentplus_2%j, &
                  nodo%Currentplus_2%k, &
                  dir(nodo%Currentplus_2%tipofield), &
                  nodo%Currentplus_2%indexmed
                  if ((nodo%k >=  ZI).and.(nodo%k<=ZE)) call WarnErrReport(buff,paraerr)
               endif
               if (nodo%NumCurrentplus >= 3) then
                  write (buff,'(a,i7,3i7,a,i7)')  '         Segment +3: ',  &
                  nodo%Currentplus_3%origindex, &
                  nodo%Currentplus_3%i, &
                  nodo%Currentplus_3%j, &
                  nodo%Currentplus_3%k, &
                  dir(nodo%Currentplus_3%tipofield), &
                  nodo%Currentplus_3%indexmed
                  if ((nodo%k >=  ZI).and.(nodo%k<=ZE)) call WarnErrReport(buff,paraerr)
               endif
               if (nodo%NumCurrentplus >= 4) then
                  write (buff,'(a,i7,3i7,a,i7)')  '         Segment +4: ',  &
                  nodo%Currentplus_4%origindex, &
                  nodo%Currentplus_4%i, &
                  nodo%Currentplus_4%j, &
                  nodo%Currentplus_4%k, &
                  dir(nodo%Currentplus_4%tipofield), &
                  nodo%Currentplus_4%indexmed
                  if ((nodo%k >=  ZI).and.(nodo%k<=ZE)) call WarnErrReport(buff,paraerr)
               endif
               if (nodo%NumCurrentplus >= 5) then
                  write (buff,'(a,i7,3i7,a,i7)')  '         Segment +5: ',  &
                  nodo%Currentplus_5%origindex, &
                  nodo%Currentplus_5%i, &
                  nodo%Currentplus_5%j, &
                  nodo%Currentplus_5%k, &
                  dir(nodo%Currentplus_5%tipofield), &
                  nodo%Currentplus_5%indexmed
                  if ((nodo%k >=  ZI).and.(nodo%k<=ZE)) call WarnErrReport(buff,paraerr)
               endif
               if (nodo%NumCurrentplus >= 6) then
                  write (buff,'(a,i7,3i7,a,i7)')  '         Segment +6: ',  &
                  nodo%Currentplus_6%origindex, &
                  nodo%Currentplus_6%i, &
                  nodo%Currentplus_6%j, &
                  nodo%Currentplus_6%k, &
                  dir(nodo%Currentplus_6%tipofield), &
                  nodo%Currentplus_6%indexmed
                  if ((nodo%k >=  ZI).and.(nodo%k<=ZE)) call WarnErrReport(buff,paraerr)
               endif
               if (nodo%NumCurrentplus >= 7) then
                  write (buff,'(a,i7,3i7,a,i7)')  '         Segment +7: ',  &
                  nodo%Currentplus_7%origindex, &
                  nodo%Currentplus_7%i, &
                  nodo%Currentplus_7%j, &
                  nodo%Currentplus_7%k, &
                  dir(nodo%Currentplus_7%tipofield), &
                  nodo%Currentplus_7%indexmed
                  if ((nodo%k >=  ZI).and.(nodo%k<=ZE)) call WarnErrReport(buff,paraerr)
               endif
               if (nodo%NumCurrentplus >= 8) then
                  write (buff,'(a,i7,3i7,a,i7)')  '         Segment +8: ',  &
                  nodo%Currentplus_8%origindex, &
                  nodo%Currentplus_8%i, &
                  nodo%Currentplus_8%j, &
                  nodo%Currentplus_8%k, &
                  dir(nodo%Currentplus_8%tipofield), &
                  nodo%Currentplus_8%indexmed
                  if ((nodo%k >=  ZI).and.(nodo%k<=ZE)) call WarnErrReport(buff,paraerr)
               endif
               if (nodo%NumCurrentplus >= 9) then
                  write (buff,'(a,i7,3i7,a,i7)')  '         Segment +9: ',  &
                  nodo%Currentplus_9%origindex, &
                  nodo%Currentplus_9%i, &
                  nodo%Currentplus_9%j, &
                  nodo%Currentplus_9%k, &
                  dir(nodo%Currentplus_9%tipofield), &
                  nodo%Currentplus_9%indexmed
                  if ((nodo%k >=  ZI).and.(nodo%k<=ZE)) call WarnErrReport(buff,paraerr)
               endif
            endif
         endif
      end do



#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif


      !!!!!!!!!!!!!!!!

      allocate (nodosopentoair(1:HWires%NumChargeNodes))
      indio=0
      do i1=1,HWires%NumChargeNodes
         nodo => Hwires%ChargeNode(i1)
         if (nodo%exists) then
            if ((nodo%numcurrentPlus+nodo%numcurrentMinus < 2)) then
               if ((.not.(nodo%IsPec.or.nodo%IsLossy))) then
                  indio=indio+1
                  nodosopentoair(indio)%indexnode=nodo%indexnode
                  nodosopentoair(indio)%i=nodo%i
                  nodosopentoair(indio)%j=nodo%j
                  nodosopentoair(indio)%k=nodo%k
                  write (buff,'(a,i7,3i7)')  'wir3_WARNING: Node open to air ',nodo%indexnode, &
                  nodo%i,nodo%j,nodo%k
                  if ((nodo%k >=  ZI).and.(nodo%k <= ZE)) call WarnErrReport(buff)
               else
                     write (buff,'(a,i7,3i7)')  'wir3_INFO: NON-JUNCTION Node  GROUNDED ',nodo%indexnode, &
                     nodo%i,nodo%j,nodo%k
                     if ((nodo%k >=  ZI).and.(nodo%k <= ZE).and.verbose) call WarnErrReport(buff)
               endif
            endif
         endif
      end do
      
      do i1=1,indio
        do j1=i1+1,indio
            if ( (nodosopentoair(i1)%i==nodosopentoair(j1)%i).and.(nodosopentoair(i1)%j==nodosopentoair(j1)%j).and.(nodosopentoair(i1)%k==nodosopentoair(j1)%k) ) then
                  write (buff,'(a,8i7)')  'wir3_ERROR: TWO nodes at the same location open to air. Should be connected? ',indio,i1,j1, &
                  nodosopentoair(i1)%indexnode, nodosopentoair(j1)%indexnode, nodosopentoair(i1)%i,nodosopentoair(i1)%j,nodosopentoair(i1)%k
                  if ((nodosopentoair(i1)%k >=  ZI).and.(nodosopentoair(i1)%k <= ZE)) call WarnErrReport(buff,.true.)
            endif
        end do
      end do
      deallocate (nodosopentoair)



#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif
      !more testing

      do i1=1,HWires%NumCurrentSegments
         if (.not.(HWires%CurrentSegment(i1)%chargeplus%exists.and.HWires%CurrentSegment(i1)%chargeminus%exists)) then
            write (buff,'(a,i7,3i7,a,2i3,a)') 'wir3_BUGGYERROR: Bug in WIRE node assignment. ', &
            nodo%indexnode,nodo%i,nodo%j,nodo%k, &
            ' (',nodo%numcurrentminus,nodo%numcurrentplus,')'
            call WarnErrReport(buff,.true.)
         endif
      end do


      !!!!writes the lines in a DXF file if requested with -map
      !!!    do i1=1,HWires%NumCurrentSegments
      !!!        write(dxfbuff,'(a)') 'LINE'
      !!!        CALL DXFWRITE(DXFBUFF)
      !!!        write(dxfbuff,'(a)') '8'
      !!!        CALL DXFWRITE(DXFBUFF)
      !!!        write(dxfbuff,'(i7)') HWires%CurrentSegment(i1)%indexmed+20
      !!!        CALL DXFWRITE(DXFBUFF)
      !!!        write(dxfbuff,'(a)') '62'
      !!!        CALL DXFWRITE(DXFBUFF)
      !!!        write(dxfbuff,'(i7)') HWires%CurrentSegment(i1)%indexmed+20
      !!!        CALL DXFWRITE(DXFBUFF)
      !!!    select case(HWires%CurrentSegment(i1)%tipofield)
      !!!    case(iEx)
      !!!        write(dxfbuff,'(a)') '10'
      !!!        CALL DXFWRITE(DXFBUFF)
      !!!        write(dxfbuff,'(i7)') HWires%CurrentSegment(i1)%i
      !!!        CALL DXFWRITE(DXFBUFF)
      !!!        write(dxfbuff,'(a)') '20'
      !!!        CALL DXFWRITE(DXFBUFF)
      !!!        write(dxfbuff,'(i7)') HWires%CurrentSegment(i1)%J
      !!!        CALL DXFWRITE(DXFBUFF)
      !!!        write(dxfbuff,'(a)') '30'
      !!!        CALL DXFWRITE(DXFBUFF)
      !!!        write(dxfbuff,'(i7)') HWires%CurrentSegment(i1)%K
      !!!        CALL DXFWRITE(DXFBUFF)
      !!!        write(dxfbuff,'(a)') '11'
      !!!        CALL DXFWRITE(DXFBUFF)
      !!!        write(dxfbuff,'(i7)') HWires%CurrentSegment(i1)%i+1
      !!!        CALL DXFWRITE(DXFBUFF)
      !!!        write(dxfbuff,'(a)') '21'
      !!!        CALL DXFWRITE(DXFBUFF)
      !!!        write(dxfbuff,'(i7)') HWires%CurrentSegment(i1)%J
      !!!        CALL DXFWRITE(DXFBUFF)
      !!!        write(dxfbuff,'(a)') '31'
      !!!        CALL DXFWRITE(DXFBUFF)
      !!!        write(dxfbuff,'(i7)') HWires%CurrentSegment(i1)%K
      !!!        CALL DXFWRITE(DXFBUFF)
      !!!        write(dxfbuff,'(a)') '0'
      !!!        CALL DXFWRITE(DXFBUFF)
      !!!    case(iEy)
      !!!        write(dxfbuff,'(a)') '10'
      !!!        CALL DXFWRITE(DXFBUFF)
      !!!        write(dxfbuff,'(i7)') HWires%CurrentSegment(i1)%i
      !!!        CALL DXFWRITE(DXFBUFF)
      !!!        write(dxfbuff,'(a)') '20'
      !!!        CALL DXFWRITE(DXFBUFF)
      !!!        write(dxfbuff,'(i7)') HWires%CurrentSegment(i1)%J
      !!!        CALL DXFWRITE(DXFBUFF)
      !!!        write(dxfbuff,'(a)') '30'
      !!!        CALL DXFWRITE(DXFBUFF)
      !!!        write(dxfbuff,'(i7)') HWires%CurrentSegment(i1)%K
      !!!        CALL DXFWRITE(DXFBUFF)
      !!!        write(dxfbuff,'(a)') '11'
      !!!        CALL DXFWRITE(DXFBUFF)
      !!!        write(dxfbuff,'(i7)') HWires%CurrentSegment(i1)%i
      !!!        CALL DXFWRITE(DXFBUFF)
      !!!        write(dxfbuff,'(a)') '21'
      !!!        CALL DXFWRITE(DXFBUFF)
      !!!        write(dxfbuff,'(i7)') HWires%CurrentSegment(i1)%J+1
      !!!        CALL DXFWRITE(DXFBUFF)
      !!!        write(dxfbuff,'(a)') '31'
      !!!        CALL DXFWRITE(DXFBUFF)
      !!!        write(dxfbuff,'(i7)') HWires%CurrentSegment(i1)%K
      !!!        CALL DXFWRITE(DXFBUFF)
      !!!        write(dxfbuff,'(a)') '0'
      !!!        CALL DXFWRITE(DXFBUFF)
      !!!    case(iEz)
      !!!        write(dxfbuff,'(a)') '10'
      !!!        CALL DXFWRITE(DXFBUFF)
      !!!        write(dxfbuff,'(i7)') HWires%CurrentSegment(i1)%i
      !!!        CALL DXFWRITE(DXFBUFF)
      !!!        write(dxfbuff,'(a)') '20'
      !!!        CALL DXFWRITE(DXFBUFF)
      !!!        write(dxfbuff,'(i7)') HWires%CurrentSegment(i1)%J
      !!!        CALL DXFWRITE(DXFBUFF)
      !!!        write(dxfbuff,'(a)') '30'
      !!!        CALL DXFWRITE(DXFBUFF)
      !!!        write(dxfbuff,'(i7)') HWires%CurrentSegment(i1)%K
      !!!        CALL DXFWRITE(DXFBUFF)
      !!!        write(dxfbuff,'(a)') '11'
      !!!        CALL DXFWRITE(DXFBUFF)
      !!!        write(dxfbuff,'(i7)') HWires%CurrentSegment(i1)%i
      !!!        CALL DXFWRITE(DXFBUFF)
      !!!        write(dxfbuff,'(a)') '21'
      !!!        CALL DXFWRITE(DXFBUFF)
      !!!        write(dxfbuff,'(i7)') HWires%CurrentSegment(i1)%J
      !!!        CALL DXFWRITE(DXFBUFF)
      !!!        write(dxfbuff,'(a)') '31'
      !!!        CALL DXFWRITE(DXFBUFF)
      !!!        write(dxfbuff,'(i7)') HWires%CurrentSegment(i1)%K+1
      !!!        CALL DXFWRITE(DXFBUFF)
      !!!        write(dxfbuff,'(a)') '0'
      !!!        CALL DXFWRITE(DXFBUFF)
      !!!    end select
      !!!    end do
      !!!!END DXFOUT by LAYER 0 and COLOR 0
      !!!

#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif



      return
   end subroutine ReportWireJunctions

   


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Functions to get the values of the
   !!! non diagonal elementos of the autoinduction
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   real (KIND=RKIND_wires) function F(A, B, a1, a2, d, phi)
      integer (KIND=4)     :: i, j
      real    (KIND=RKIND_wires) :: A, B, a1, a2, d, phi, frac

      frac = 1.0_RKIND_wires / (8.0_RKIND_wires * A*B);
      F = pi*a2*a2*(1.0_RKIND_wires-2.0_RKIND_wires * log(a2/d))-3.0_RKIND_wires / (2.0_RKIND_wires * frac);
      do i = 1,2
         do j = 1,2
            F = F+Gkl(i,j,A,B,d,phi)+Hkl(i,j,A,B,d,phi)
         end do
      end do

      F = F*frac;
      return
   end function

   real (KIND=RKIND_wires) function Gkl(k, l, A, B, d, phi)
      integer (KIND=4)     :: k, l
      real    (KIND=RKIND_wires) :: A, B, d, phi, Ak, Bl
      Ak = Ai(k, A, d, phi)
      Bl = Bi(l, B, d, phi)
      Gkl = Ak*Bl*log((Ak*Ak+Bl*Bl)/(d*d))
      return
   end function

   real (KIND=RKIND_wires) function Hkl(k, l, A, B, d, phi)
      integer (KIND=4)     :: k, l
      real    (KIND=RKIND_wires) :: A, B, d, phi, Ak, Bl
      Ak = Ai(k, A, d, phi)
      Bl = Bi(l, B, d, phi)
      Hkl = Ak*Ak*atan(Bl/Ak)+Bl*Bl*atan(Ak/Bl)
      return
   end function

   real (KIND=RKIND_wires) function Ai(i, A, d, phi)
      integer (KIND=4)     :: i
      real    (KIND=RKIND_wires) :: A, d, phi
      ai=-1.0
      select case(i)
       case (1)
         Ai=A-d*cos(phi);
       case (2)
         Ai=A+d*cos(phi);
      end select

      return
   end function

   real (KIND=RKIND_wires) function Bi(i, B, d, phi)
      integer (KIND=4)     :: i
      real    (KIND=RKIND_wires) :: B, d, phi
      bi=-1.0
      select case(i)
       case (1)
         Bi = B-d*sin(phi);
       case (2)
         Bi = B+d*sin(phi);
      end select

      return
   end function

   !Rutina de inversionde matrices
   subroutine MatInv(N, M)
      !soubroutine parameters
      integer (kind=4), intent(in)                ::  N
      real    (kind=RKIND_wires), intent(inout)             ::  M(1:N, 1:N)

      !local variables
      real      (kind=RKIND_wires), pointer, dimension(:,:) ::  B, eye
      integer   (kind=4), pointer, dimension(:)   ::  P
      real      (kind=RKIND_wires), pointer, dimension(:)   ::  y
      integer   (kind=4)                          ::  i, j, k
      integer   (kind=4)                          ::  pivot, tmpi
      real      (kind=RKIND_wires)                          ::  tmpr, val
      character (len=BUFSIZE)                         ::  buff
      pivot=-1
      !function body
      allocate(P(1:N),y(1:N))
      allocate(B(1:N,1:N))
      allocate(eye(1:N,1:N))
      B(1:N,1:N) = M(1:N,1:N)

      do i = 1,N
         do j = 1,N
            eye(i,j) = 0.0_RKIND_wires
         end do
         eye(i,i) = 1.0_RKIND_wires
      end do

      do i = 1,N
         P(i) = i
      end do

      do k=1,N
         val = 0.0_RKIND_wires
         do i = k,N
            if(abs(B(i,k)) > val) then
               val = abs(B(i,k))
               pivot = i
            end if
         end do

         if(val == 0.0_RKIND_wires) then
            write(buff,*) 'WIR1_ERROR: Inversion de matriz fallida'
            call WarnErrReport(buff,.true.)
            deallocate(P,B,y,eye)
            return
         end if

         tmpi     = P(k)
         P(k)     = P(pivot)
         P(pivot) = tmpi

         do i = 1,N
            tmpr        = B(k,i)
            B(k,i)      = B(pivot,i)
            B(pivot,i)  = tmpr
         end do

         do i = k+1,N
            B(i,k) = B(i,k)/B(k,k)
            do j = k+1,N
               B(i,j) = B(i,j)-B(i,k)*B(k,j)
            end do
         end do
      end do

      do k = 1,N
         do i = 1,N
            y(i) = eye(P(i),k)
            do j = 1,i-1
               y(i) = y(i) - B(i,j)*y(j)
            end do
         end do
         do i=N,1,-1
            M(i,k) = y(i)
            do j=i+1,N
               M(i,k) = M(i,k) - B(i,j)*M(j,k)
            end do
            M(i,k) = M(i,k)/B(i,i)
         end do
      end do

      deallocate(P,B,y,eye)
   end subroutine

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !End of auxiliar functions
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 

   function Lambert(z)  result(x)
      implicit none
      real (kind=RKIND_wires) :: x,expx,newx,z
      real (kind=RKIND_wires) :: eps
      integer(kind=4) :: maxn=2**12,n
      x=z
      eps=abs(x/100000.0)
      if (z<-1.0_RKIND_wires/3.0_RKIND_wires) then
         X=-1
         return !no va a converger nunca
      endif
      n=1
      loop:do while ((abs(x*Exp(x)-z)>eps).and.(n<maxn))
         expx = exp(x)
         !!!newton        newx = x - (x * expx - z) / (expx * x + expx)
         !!!otro mejor
         newx = x - (x * expx - z)/(expx * (x + 1.0_RKIND_wires) - (x + 2.0_RKIND_wires)*(x* expx - z)/(2.0_RKIND_wires * x + 2.0_RKIND_wires))
         x=newx
         n=n+1
      end do loop
      expx = exp(x)
      if (n>=maxn) then
         X=-1
      endif
      return
   end function Lambert

  !!!!!!tridiagonal solver

      subroutine solve_tridiag_wires(a,b,c,d,x,n)
      implicit none
      !  a - sub-diagonal (means it is the diagonal below the main diagonal)
      !  b - the main diagonal
      !  c - sup-diagonal (means it is the diagonal above the main diagonal)
      !  d - right part
      !  x - the answer
      !  n - number of equations

      integer,intent(in) :: n
      real (kind=RKIND_wires) ,dimension(n),intent(in) :: a,b,c,d
      real (kind=RKIND_wires) ,dimension(n),intent(out) :: x
      real (kind=RKIND_wires) ,dimension(n) :: cp,dp
      real (kind=RKIND_wires)  :: m
      integer i

      !  initialize c-prime and d-prime
      cp(1) = c(1)/b(1)
      dp(1) = d(1)/b(1)
      ! solve for vectors c-prime and d-prime
      do i = 2,n
         m = b(i)-cp(i-1)*a(i)
         cp(i) = c(i)/m
         dp(i) = (d(i)-dp(i-1)*a(i))/m
      enddo
      ! initialize x
      x(n) = dp(n)
      ! solve for x from the vectors c-prime and d-prime
      do i = n-1, 1, -1
         x(i) = dp(i)-cp(i)*x(i+1)
      end do
      return
      end subroutine solve_tridiag_wires

      subroutine wiresconstantes(fieldtotl,dummy,G2,sgg)
      
          type (SGGFDTDINFO), intent(IN) , target    ::  sgg
          REAL (KIND=RKIND) , pointer, dimension (:), intent(in)      :: G2
          
          logical, intent(in) :: fieldtotl
          type (CurrentSegments), pointer  ::  dummy
!!!ojooooo 110517 acumulo en %lind toda la autoinduccion para que los calculos de capacidad la tengan en cuenta completa    
         if (.not.fieldtotl) then
             dummy%cte5 = G2(dummy%indexmed)/(dummy%deltaTransv1*dummy%deltaTransv2)
         else
             dummy%cte5 = 0.0_RKIND_wires !esta es la cte de acoplo wire2main 
         endif !para el fieldtotl no tento en cuenta mas que la autoin y no devuelvo nada !100517
!!!!
!
!Lind ya contiene el givenautoin y a los autoin LeftEnd y RightEnd, que entiendo que de escalarse el mu tambien afectaria a todos ellos !ojo
         dummy%cte1= ((dummy%Lind )/sgg%dt-dummy%Resist/2.0_RKIND_wires) &
                    /((dummy%Lind )/sgg%dt+dummy%Resist/2.0_RKIND_wires)
         dummy%cte3= InvMu(dummy%indexmed) * InvEps(dummy%indexmed) / (dummy%delta) * (dummy%Lind) &
                    /((dummy%Lind )/sgg%dt+dummy%Resist/2.0_RKIND_wires)
         dummy%cte2= 1.0_RKIND_wires  &
                    /((dummy%Lind )/sgg%dt+dummy%Resist/2.0_RKIND_wires)
!!!stochastic

#ifdef CompileWithStochastic
         call calc_wirehollandconstants_for_devia(sgg,dummy,InvEps,InvMu)
#endif         
!         continue


    end subroutine wiresconstantes

#endif

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module HollandWires
