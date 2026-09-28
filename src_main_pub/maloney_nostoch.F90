
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module SGBCs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!17/08/15 update!!!!!!!!!!!
!!!Removed the treatment of the SGBC magnetic fields to program a multiSGBC 
!!!only taking into account the effective parameters and without updating the magnetics.
!!!I keep the old version in the file SGBC_pre170815_noupdateababienH.F90
!!!
!!! 211115 THE AVERAGING OF filo_placaS THAT IS DONE IN COMPOSITES IS NOT DONE. 
!!!!       SGBCS ARE ASSIGNED ON A FIRST-COME-FIRST-SERVE BASIS. 
!!!!       THE filo_placaS TREATMENT DETECTS EDGES AND DOES SOMETHING SIMILAR TO THE SHARED ONE
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module SGBC_nostoch_m

use Report_m

use FDETYPES_m
implicit none
private
!structures needed by the SGBC
type :: val_t
    complex(kind=CKIND), allocatable, dimension(:) :: val
end type


type :: MDfield_t
   real(kind=RKIND), pointer                 :: FieldPresent !points to the background field
   real(kind=RKIND)                          :: FieldPrevious
   complex(kind=CKIND), pointer, dimension(:) :: Current
end type


type  :: SGBCSurface_t
   real(kind=RKIND), allocatable, dimension(:) :: E,H,E_past
   real(kind=RKIND), pointer :: Efield,Ha_Plus,Ha_Minu,Hb_Plus,Hb_Minu
   real(kind=RKIND), allocatable, dimension(:) :: delta_entreEinterno
   real(kind=RKIND), dimension(0:1) :: g1,g2a,g2b
!!!SGBC dispersive 12/05/16
   type(MDfield_t), allocatable, dimension(:) :: EDis
   integer(kind=4) :: numpolres
   type(val_t) :: Beta,Kappa,G3
!!!!!!!      
   logical :: correct_ha, correct_hb, es_unfilo_placa
   
   integer(kind=4) :: depth,jmed
   integer(kind=4), allocatable, dimension(:) ::layerIndex !!!0121
   real(kind=RKIND) , allocatable, dimension(:) :: G2_interno,GM2_interno,G1_interno,GM1_interno   
   real(kind=RKIND) :: GM2_externo   !gm1_externo is not needed because outside there is no magnetic conductivity and it is trivially 1. Storing gm2_externo makes sense because, even without conductivity, it is not unity
   real(kind=RKIND) :: Hyee__left, Hyee_right      
!!!!! Crank-Nicolson 311015
   real(kind=RKIND) , allocatable, dimension(:) :: a,b,c,rb,rh,rhm1
   real(kind=RKIND)                                :: a1,b1,c1,rb1,rh1,an,bn,cn,rbn,rhn 
   real(kind=RKIND) , allocatable, dimension(:) :: D !independent term CRANK-NICOLSON
   logical :: SGBCCrank

   real(kind=RKIND) :: transversalDeltaE,transversalDeltaH,alignedlDeltaH
   integer(kind=4), dimension(0:1) :: med
   complex(kind=ckind), allocatable, dimension(:) :: a11, c11
end type SGBCSurface_t


type :: MalDisp_t
    integer(kind=4) :: numpolres
    complex(kind=ckind), allocatable, dimension(:) :: a11, c11
end type

type  :: Malon_t
    logical :: SGBCdispersive
   integer(kind=4) :: NumNodes
   type(SGBCSurface_t), allocatable, dimension(:) :: nodes
   type(MalDisp_t), allocatable, dimension(:) :: dispersiveMedia
end type Malon_t



!!!module global variables  
type(Malon_t), save, target   :: malon
!
real(kind=RKIND), save           :: eps0,mu0,zvac,cluz
logical, save  :: SGBCcrank,SGBCDispersive
real(kind=RKIND), save  :: SGBCFreq,SGBCresol
integer(kind=4), save:: SGBCdepth
!!!
public Malon_t,SGBCSurface_t !the type is public
public AdvanceSGBCE,AdvanceSGBCH,InitSGBCs,DestroySGBCs,StoreFieldsSGBCs,calc_SGBCconstants,GetSGBCs
public solve_tridiag_iguales

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine to initialize the parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine InitSGBCs(sgg,media,Ex,Ey,Ez,Hx,Hy,Hz,IDxe,IDye,IDze,IDxh,IDyh,IDzh, &
                     layoutnumber,num_procs,g,ThereAreSGBCs,resume, &
                     temp_SGBCcrank,temp_SGBCFreq,temp_SGBCresol,temp_SGBCDepth,temp_SGBCDispersive, &
                     eps00,mu00,simu_devia,stochastic)
   type(media_matrices_t), intent(in) :: media
   logical :: simu_devia,stochastic 
   real(kind=RKIND) :: eps00,mu00
   logical :: temp_SGBCcrank,temp_SGBCDispersive
   type(constants_t), intent(inout) :: g
   type(SGGFDTDINFO_t), intent(inout) :: sgg !careful because epr, mur, sigma, sigmam are overwritten for dispersive materials
   real(kind=RKIND)   , intent(in) , target     :: &
   Ex(sgg%alloc(iEx)%XI : sgg%alloc(iEx)%XE,sgg%alloc(iEx)%YI : sgg%alloc(iEx)%YE,sgg%alloc(iEx)%ZI : sgg%alloc(iEx)%ZE),&
   Ey(sgg%alloc(iEy)%XI : sgg%alloc(iEy)%XE,sgg%alloc(iEy)%YI : sgg%alloc(iEy)%YE,sgg%alloc(iEy)%ZI : sgg%alloc(iEy)%ZE),&
   Ez(sgg%alloc(IEZ)%XI : sgg%alloc(IEZ)%XE,sgg%alloc(IEZ)%YI : sgg%alloc(IEZ)%YE,sgg%alloc(IEZ)%ZI : sgg%alloc(IEZ)%ZE),&
   Hx(sgg%alloc(IHX)%XI : sgg%alloc(IHX)%XE,sgg%alloc(IHX)%YI : sgg%alloc(IHX)%YE,sgg%alloc(IHX)%ZI : sgg%alloc(IHX)%ZE),&
   Hy(sgg%alloc(IHY)%XI : sgg%alloc(IHY)%XE,sgg%alloc(IHY)%YI : sgg%alloc(IHY)%YE,sgg%alloc(IHY)%ZI : sgg%alloc(IHY)%ZE),&
   Hz(sgg%alloc(IHZ)%XI : sgg%alloc(IHZ)%XE,sgg%alloc(IHZ)%YI : sgg%alloc(IHZ)%YE,sgg%alloc(IHZ)%ZI : sgg%alloc(IHZ)%ZE)
   real(kind=RKIND) , dimension(:)   , intent(in) :: Idxh(sgg%ALLOC(iEx)%XI : sgg%ALLOC(iEx)%XE), &
                                                      &  Idyh(sgg%ALLOC(iEy)%YI : sgg%ALLOC(iEy)%YE), &
                                                      &  Idzh(sgg%ALLOC(IEZ)%ZI : sgg%ALLOC(IEZ)%ZE), &
                                                         Idxe(sgg%alloc(IHX)%XI : sgg%alloc(IHX)%XE), &
                                                         Idye(sgg%alloc(IHY)%YI : sgg%alloc(IHY)%YE), &
                                                         Idze(sgg%alloc(IHZ)%ZI : sgg%alloc(IHZ)%ZE)

   real(kind=RKIND) :: temp_SGBCFreq,temp_SGBCresol, rra,rrb,rrc,rrd
   real(kind=RKIND) :: signo,g1eff_0,g1eff_1,g2eff_0,g2eff_1,Sigmam,epsilonValue,Mu,Sigma   
   real(kind=RKIND) :: factor
   real(kind=RKIND) , allocatable, dimension(:,:) :: derivcte

   integer(kind=4), intent(in) :: layoutnumber,num_procs,temp_SGBCdepth
   logical, intent(in) :: resume
   logical, intent(out) :: ThereAreSGBCs
   integer(kind=4) :: jmed,j1,conta,k1,i1,SGBCdir,i,filo_placas,idummy,numpolres,ient,incert,maxnumcapas,ii
   character(len=BUFSIZE) :: buff
   character(len=BUFSIZE) :: whoami
   type(SGBCSurface_t), pointer :: compo,compo_temp
   logical :: unstable, errnofile,es_unfilo_placa
   complex(kind=ckind) :: value1, value2
   character(len=BUFSIZE)                            :: filePoles

   eps0=eps00; mu0=mu00; !hack to turn the step variables into globals
   SGBCcrank        = temp_SGBCcrank     
   SGBCDispersive   = temp_SGBCDispersive
   SGBCFreq         = temp_SGBCFreq      
   SGBCresol        = temp_SGBCresol     
   SGBCdepth        = temp_SGBCdepth     
!
!!!
   write(whoami,'(a,i5,a,i5,a)') '(',layoutnumber+1,'/',num_procs,') '
   unstable=.false.

!

!   
   malon%SGBCDispersive=SGBCDispersive
   ThereAreSGBCs=.FALSE.
   do jmed=1,sgg%NumMedia
      if (SGG%Med(jmed)%Is%SGBC) then
         ThereAreSGBCs=.true.
      end if
   end do

!pre-count the media
   conta=0
   do k1=sgg%SINPMLSweep(iEx)%ZI,sgg%SINPMLSweep(iEx)%ZE
      do j1=sgg%SINPMLSweep(iEx)%YI,sgg%SINPMLSweep(iEx)%YE
         do i1=sgg%SINPMLSweep(iEx)%XI,sgg%SINPMLSweep(iEx)%XE
            jmed=media%sggMiEx(i1,j1,k1)
            if (SGG%Med(jmed)%Is%SGBC)  conta=conta+1
         end do
      end do
   end do
   do k1=sgg%SINPMLSweep(iEy)%ZI,sgg%SINPMLSweep(iEy)%ZE
      do j1=sgg%SINPMLSweep(iEy)%YI,sgg%SINPMLSweep(iEy)%YE
         do i1=sgg%SINPMLSweep(iEy)%XI,sgg%SINPMLSweep(iEy)%XE
            jmed=media%sggMiEy(i1,j1,k1)
            if (SGG%Med(jmed)%Is%SGBC)  conta=conta+1
         end do
      end do
   end do

   do k1=sgg%SINPMLSweep(IEZ)%ZI,sgg%SINPMLSweep(IEZ)%ZE
      do j1=sgg%SINPMLSweep(IEZ)%YI,sgg%SINPMLSweep(IEZ)%YE
         do i1=sgg%SINPMLSweep(IEZ)%XI,sgg%SINPMLSweep(IEZ)%XE
            jmed=media%sggMiEz(i1,j1,k1)
            if (SGG%Med(jmed)%Is%SGBC) conta=conta+1
         end do
      end do
   end do
   !!!!!!!!!!!!!!!!!!!!!!
   ThereAreSGBCs=ThereAreSGBCs.and.(conta /=0)
   if (.not.thereareSGBCs) then
      return
   end if
   malon%NumNodes=conta
   allocate (malon%Nodes(1 : malon%NumNodes))
   !!!!DISPERSIVE
   allocate (malon%dispersiveMedia(1:sgg%NumMedia))
   malon%dispersiveMedia(:)%numpolres=0
   !
   !!!!!!!! dispersive SGBC sgg 12/05/15   
!070717
!!!first I check whether the poles file exists. If it exists, it means SGBCDispersive stopped and reported
!!!!note that for now the sgbc thing is global (070717). If the switch is on all are dispersive sgbc. but it is half-prepared to be medium by medium. change someday....
!!!!0121 I remove this because it looks at the poles file generated by ugrmat_multilayer and enabled dispersive by default just like that
   !! do jmed=1,sgg%NumMedia
   !!      if ((.not.SGG%Med(jmed)%Is%SGBCDispersive).and.(SGG%Med(jmed)%Is%SGBC).and.(.not.(SGG%Med(jmed)%Is%PML))) then    
   !!          ficheropolos=SGG%Med(jmed)%multiport(1)%multiportFileZ11 !although I call it Z it has the syntax of an Edispersive ISOTROPIC WITH THE NEW STANDARD (SEE LINES 6749 OF NFDEPARSER). 
   !!          ! I ONLY READ THE FIRST POLES. THE REST OF THE DATA I DISCARD (SECOND-ORDER POLE INFORMATION, MAGNETIC POLES, ANISOTROPIES...)
   !!          !new file style without the _z11
   !!          i1=index(ficheropolos,'_z11.txt')
   !!          ficheropolos=trim(adjustl(ficheropolos(1:i1-1)))
   !!  !
   !!          errnofile=.false.
   !!          inquire(FILE=trim(adjustl(ficheropolos)), EXIST=errnofile)
   !!          if (errnofile) then
   !!               write(buff, *)    'ERROR: -sgbcdispersive not used and poles files exist. Correct .nfde or issue -sgbcdispersive'
   !!               call WarnErrReport (buff,.false.)
   !!               stop
   !!               SGG%Med(jmed)%Is%SGBCDispersive=.true.
   !!               SGBCdispersive=.true. !in case of ignoreerrors it can just continue with sgbcdispersive set to true
   !!          end if
   !!      end if
   !!end do 
   if (SGBCDispersive) then
       do jmed=1,sgg%NumMedia
          if ((SGG%Med(jmed)%Is%SGBCDispersive).and.(.not.(SGG%Med(jmed)%Is%PML))) then    
!!!only one dispersive layer
              if (sgg%Med(jmed)%multiport(1)%numLayers>1) then
                 buff='No more than 1 layer of dispersive SGBC currently supported'
                 call StopOnError(layoutnumber,num_procs,buff)
              end if
!!!!!!!!
              filePoles=SGG%Med(jmed)%multiport(1)%multiportFileZ11 !although I call it Z it has the syntax of an Edispersive ISOTROPIC WITH THE NEW STANDARD (SEE LINES 6749 OF NFDEPARSER). 
              ! I ONLY READ THE FIRST POLES. THE REST OF THE DATA I DISCARD (SECOND-ORDER POLE INFORMATION, MAGNETIC POLES, ANISOTROPIES...)
              !new file style without the _z11
              i1=index(filePoles,'_z11.txt')
              filePoles=trim(adjustl(filePoles(1:i1-1)))
      !
              errnofile=.false.
              inquire(FILE=trim(adjustl(filePoles)), EXIST=errnofile)
              if (.not.errnofile) then
                 buff='FILE '//trim(adjustl(filePoles))//' DOES NOT EXIST'
                 call StopOnError(layoutnumber,num_procs,buff)
              end if
              open (7345,file=trim(adjustl(filePoles)),form='formatted')
              read (7345,*) rra,rrb,rrc,rrd 
              rrb= rrb/eps0 ;  rrc = rrc/mu0 !permit scaling does not affect them I think 071118 because they are relative to the program input which MUST COME WITH the genuine eps0 and mu0
              SGG%Med(jmed)%multiport(1)%sigma(1)=rra; 
              SGG%Med(jmed)%multiport(1)%epr(1)=rrb; 
              SGG%Med(jmed)%multiport(1)%mur(1)=rrc; 
              SGG%Med(jmed)%multiport(1)%sigmam(1)=rrd;
              SGG%Med(jmed)%sigma=rra; 
              SGG%Med(jmed)%epr=rrb; 
              SGG%Med(jmed)%mur=rrc; 
              SGG%Med(jmed)%sigmam=rrd;
              !
              read (7345,*) numpolres, IDUMMY, IDUMMY, IDUMMY
              malon%dispersiveMedia(jmed)%numpolres = numpolres
              allocate (malon%dispersiveMedia(jmed)%a11(1:numpolres)) 
              allocate (malon%dispersiveMedia(jmed)%c11(1:numpolres)) 
              do i = 1, numpolres
                read(7345,*) value1, value2
                malon%dispersiveMedia(jmed)%c11 (i) = (value1) 
                malon%dispersiveMedia(jmed)%a11 (i) = - (value2) !the EM pole has its sign flipped !see also preprocess
              end do          
              close (7345)
!!!moved 071118 to the calculation of constants for permit scaling
!!!!!!!!070717  I recalculate and overwrite G1,G2,GM1, and GM2 with those appearing in the dispersive file although I think it is not used at all
!!!!!                  Sigmam  =      SGG%Med(jmed)%multiport(1)%sigmam(1)
!!!!!                  Epsilon = Eps0*SGG%Med(jmed)%multiport(1)%epr(1)
!!!!!                  Mu      = Mu0* SGG%Med(jmed)%multiport(1)%mur(1)
!!!!!                  Sigma   =      SGG%Med(jmed)%multiport(1)%sigma(1)
!!!!!                  G1(jmed)=(1 -  Sigma * sgg%dt / (2.0_RKIND * Epsilon ) ) / (1.0_RKIND + Sigma * sgg%dt / (2.0_RKIND * Epsilon ))
!!!!!                  G2(jmed)=sgg%dt /Epsilon                        / (1.0_RKIND + Sigma * sgg%dt / (2.0_RKIND * Epsilon ))
!!!!!                  if (g1(jmed) < 0.0_RKIND) then !exponential time stepping
!!!!!                     g1(jmed)=exp(- Sigma * sgg%dt / (Epsilon ))
!!!!!                     g2(jmed)=(1.0_RKIND-g1(jmed))/ Sigma
!!!!!                  end if
!!!!!                  GM1(jmed)=(1- SigmaM*sgg%dt/(2.0_RKIND *  Mu )) /(1.0_RKIND + SigmaM*sgg%dt/(2.0_RKIND *  Mu ))
!!!!!                  GM2(jmed)=sgg%dt/ Mu                   /(1.0_RKIND + SigmaM*sgg%dt/(2.0_RKIND *  Mu ))
!!!!!                  if (gm1(jmed) < 0.0_RKIND) then !exponential time stepping
!!!!!                     gm1(jmed)=exp(- Sigmam * sgg%dt / (Mu ))
!!!!!                     gm2(jmed)=(1.0_RKIND-gm1(jmed))/ Sigmam
!!!!!                  end if
!!!!!!!!!end 070717
         end if
      end do
   end if
   !!!!!!!! end of SGBCdispersive
   conta=0
!assigns the H curl signs and the transversal deltas. 0->left plate edge, 1->right plate edge and the variable es_unfilo_placa for the edges of the SGBC sheet      
do k1=sgg%SINPMLSweep(iEx)%ZI,sgg%SINPMLSweep(iEx)%ZE
    do j1=sgg%SINPMLSweep(iEx)%YI,sgg%SINPMLSweep(iEx)%YE
         do i1=sgg%SINPMLSweep(iEx)%XI,sgg%SINPMLSweep(iEx)%XE
            jmed=media%sggMiEx(i1,j1,k1)
            if (SGG%Med(jmed)%Is%SGBC)  then
!!!
               filo_placas=0
               es_unfilo_placa=.false.
               if (SGG%Med(media%sggMiEx(i1,j1,k1+1))%Is%SGBC) filo_placas=filo_placas+1
               if (SGG%Med(media%sggMiEx(i1,j1,k1-1))%Is%SGBC) filo_placas=filo_placas+1
               if (SGG%Med(media%sggMiEx(i1,j1+1,k1))%Is%SGBC) filo_placas=filo_placas+1
               if (SGG%Med(media%sggMiEx(i1,j1-1,k1))%Is%SGBC) filo_placas=filo_placas+1
               if (filo_placas < 2) then
                     es_unfilo_placa=.true.
               end if
!!!!!
               conta=conta+1
               compo => malon%Nodes(conta)
               compo%es_unfilo_placa = es_unfilo_placa
               SGBCdir=abs(SGG%Med(jmed)%Multiport(1)%Multiportdir)
               select case (SGBCdir)
                case (iEy)
                  compo%transversalDeltaE =    1.0_RKIND/IDYe(j1)
                  compo%transversalDeltaH =    1.0_RKIND/IDYh(j1)
                  compo%alignedlDeltaH     =    1.0_RKIND/Idzh(k1)
                  compo%med(1) =                     media%sggMiHz(i1,j1  ,k1)
                  compo%med(0) =                     media%sggMiHz(i1,j1-1,k1)
                  compo%Correct_Ha=.true. !they are cyclic a,b -> x,y,z
                  compo%Correct_Hb=.false.
                case (IEZ)
                  compo%transversalDeltaE    = 1.0_RKIND/IDze(k1)
                  compo%transversalDeltaH    = 1.0_RKIND/IDzh(k1)
                  compo%alignedlDeltaH        = 1.0_RKIND/IDYh(j1)
                  compo%med(1) =                     media%sggMiHy(i1,j1,k1)
                  compo%med(0) =                     media%sggMiHy(i1,j1,k1-1)
                  compo%Correct_Hb=.true.
                  compo%Correct_Ha=.false.
                case DEFAULT
                  write(buff, *)    'Buggy ERROR: In SGBCs. '
                  call WarnErrReport (buff,.TRUE.)
               end select
               if ((compo%med(1) == 0).or.(compo%med(0) == 0).or.(sgg%med(compo%med(1))%Is%PEC).or.(sgg%med(compo%med(0))%Is%PEC)) then
                  write(buff, *)    'Buggy ERROR: In SGBCs. '
                  call WarnErrReport (buff,.TRUE.)
               end if
!!!!!!!!!
               compo%jmed       = jmed
               compo%Efield    => Ex(i1,j1  ,k1)
               compo%Ha_Plus   => Hz(i1,j1  ,k1)
               compo%Ha_Minu   => Hz(i1,j1-1,k1)
               compo%Hb_Plus   => Hy(i1,j1  ,k1)
               compo%Hb_Minu   => Hy(i1,j1  ,k1-1)
               call  depth(compo,sgg,jmed,SGBCFreq,SGBCresol,SGBCdepth)
               if (compo%depth==0) then
                    compo%SGBCCrank=.false.
               else if (compo%depth<0) then
                  write(buff, *)    'Buggy ERROR: In SGBCs compo%depth<0. ',compo%depth
                  call WarnErrReport (buff,.TRUE.)
               else
                    compo%SGBCCrank=SGBCCrank
               end if
               allocate (compo%E         (-compo%depth:compo%depth))
               if (compo%depth>0) allocate (compo%H         (-compo%depth:compo%depth-1)) 
               allocate(compo%E_past(-compo%depth:compo%depth)) !not needed in yee but communicated in mpi_stochastic. moved outside the following if 170519
               if (compo%SGBCCrank)  then
                    allocate(compo%d     (-compo%depth:compo%depth)) 
               end if

               compo%numpolres=malon%dispersiveMedia(compo%jmed)%numpolres !I duplicate this info
               if (SGBCDispersive) then 
                 allocate (compo%a11(1:compo%numpolres)) 
                 allocate (compo%c11(1:compo%numpolres)) 
                 compo%a11 = malon%dispersiveMedia(compo%jmed)%a11
                 compo%c11 = malon%dispersiveMedia(compo%jmed)%c11
                 allocate (compo%beta%val(1:compo%numpolres))
                 allocate (compo%kappa%val(1:compo%numpolres))
                 allocate (compo%G3%val(1:compo%numpolres))
             !!  call calc_g1g2(sgg,GM2,compo)   ! permit scal 071118 
                 allocate (compo%EDis   (-compo%depth:compo%depth))
                 do ient=-compo%depth , compo%depth
                     allocate (compo%EDis(ient)%Current(1 : malon%dispersiveMedia(jmed)%numpolres))
                     compo%EDis(ient)%FieldPresent => compo%E(ient)
                 end do
               end if
             end if
         end do
      end do
   end do
   !!!!!!!!!!
   do k1=sgg%SINPMLSweep(iEy)%ZI,sgg%SINPMLSweep(iEy)%ZE
      do j1=sgg%SINPMLSweep(iEy)%YI,sgg%SINPMLSweep(iEy)%YE
         do i1=sgg%SINPMLSweep(iEy)%XI,sgg%SINPMLSweep(iEy)%XE
            jmed=media%sggMiEy(i1,j1,k1)
            if (SGG%Med(jmed)%Is%SGBC)  then
!!!
               filo_placas=0
               es_unfilo_placa=.false.
               if (SGG%Med(media%sggMiEy(i1+1,j1,k1))%Is%SGBC) filo_placas=filo_placas+1
               if (SGG%Med(media%sggMiEy(i1-1,j1,k1))%Is%SGBC) filo_placas=filo_placas+1
               if (SGG%Med(media%sggMiEy(i1,j1,k1+1))%Is%SGBC) filo_placas=filo_placas+1
               if (SGG%Med(media%sggMiEy(i1,j1,k1-1))%Is%SGBC) filo_placas=filo_placas+1
               if (filo_placas < 2) then
                     es_unfilo_placa=.true.
               end if
!!!!!
               conta=conta+1
               compo => malon%Nodes(conta)
               compo%es_unfilo_placa = es_unfilo_placa
               SGBCdir=abs(SGG%Med(jmed)%Multiport(1)%Multiportdir)
               select case (SGBCdir)
                case (IEZ)
                  compo%transversalDeltaE = 1.0_RKIND/IDze(k1)
                  compo%transversalDeltaH = 1.0_RKIND/IDzh(k1)
                  compo%alignedlDeltaH     = 1.0_RKIND/Idxh(i1)
                  compo%med(1) =                  media%sggMiHx(i1,j1,k1)
                  compo%med(0) =                  media%sggMiHx(i1,j1,k1-1)
                  compo%Correct_Ha=.true.
                  compo%Correct_Hb=.false.
                case (iEx)
                  compo%transversalDeltaE = 1.0_RKIND/IDxe(i1)
                  compo%transversalDeltaH = 1.0_RKIND/IDxh(i1)
                  compo%alignedlDeltaH     = 1.0_RKIND/IDzh(k1)
                  compo%med(1) =                  media%sggMiHz(i1  ,j1,k1)
                  compo%med(0) =                  media%sggMiHz(i1-1,j1,k1)
                  compo%Correct_Ha=.false.
                  compo%Correct_Hb=.true.
                case DEFAULT
                  write(buff, *)    'Buggy ERROR: In SGBCs. '
                  call WarnErrReport (buff,.TRUE.)
               end select
               if ((compo%med(1) == 0).or.(compo%med(0) == 0).or.(sgg%med(compo%med(1))%Is%PEC).or.(sgg%med(compo%med(0))%Is%PEC)) then
                  write(buff, *)    'Buggy ERROR: In SGBCs. '
                  call WarnErrReport (buff,.TRUE.)
               end if
!!!!   
               compo%jmed       = jmed
               compo%Efield    => Ey(i1  ,j1  ,k1)
               compo%Ha_Plus   => Hx(i1  ,j1  ,k1)
               compo%Ha_Minu   => Hx(i1  ,j1  ,k1-1)
               compo%Hb_Plus   => Hz(i1  ,j1  ,k1)
               compo%Hb_Minu   => Hz(i1-1,j1  ,k1)
               call  depth(compo,sgg,jmed,SGBCFreq,SGBCresol,SGBCdepth)
               if (compo%depth==0) then
                    compo%SGBCCrank=.false.
               else if (compo%depth<0) then
                  write(buff, *)    'Buggy ERROR: In SGBCs compo%depth<0. ',compo%depth
                  call WarnErrReport (buff,.TRUE.)
               else
                    compo%SGBCCrank=SGBCcrank
               end if
               allocate (compo%E         (-compo%depth:compo%depth))
               if (compo%depth>0) allocate (compo%H         (-compo%depth:compo%depth-1)) 
               allocate(compo%E_past(-compo%depth:compo%depth))
               if (compo%SGBCcrank)  then
                    allocate(compo%d     (-compo%depth:compo%depth)) 
               end if

               compo%numpolres=malon%dispersiveMedia(compo%jmed)%numpolres !I duplicate this info
               if (SGBCDispersive) then 
                 allocate (compo%a11(1:compo%numpolres)) 
                 allocate (compo%c11(1:compo%numpolres)) 
                 compo%a11 = malon%dispersiveMedia(compo%jmed)%a11
                 compo%c11 = malon%dispersiveMedia(compo%jmed)%c11
                 allocate (compo%beta%val(1:compo%numpolres))
                 allocate (compo%kappa%val(1:compo%numpolres))
                 allocate (compo%G3%val(1:compo%numpolres))
                 !! call calc_g1g2(sgg,GM2,compo)     ! permit scal 071118
                 allocate (compo%EDis   (-compo%depth:compo%depth))
                 do ient=-compo%depth , compo%depth
                     allocate (compo%EDis(ient)%Current(1 : malon%dispersiveMedia(jmed)%numpolres))
                     compo%EDis(ient)%FieldPresent => compo%E(ient)
                 end do
               end if
            end if
         end do
      end do
   end do

   do k1=sgg%SINPMLSweep(IEZ)%ZI,sgg%SINPMLSweep(IEZ)%ZE
      do j1=sgg%SINPMLSweep(IEZ)%YI,sgg%SINPMLSweep(IEZ)%YE
         do i1=sgg%SINPMLSweep(IEZ)%XI,sgg%SINPMLSweep(IEZ)%XE
            jmed=media%sggMiEz(i1,j1,k1)
            if (SGG%Med(jmed)%Is%SGBC) then
!!!
               filo_placas=0
               es_unfilo_placa=.false.
               if (SGG%Med(media%sggMiEz(i1,j1+1,k1))%Is%SGBC) filo_placas=filo_placas+1
               if (SGG%Med(media%sggMiEz(i1,j1-1,k1))%Is%SGBC) filo_placas=filo_placas+1
               if (SGG%Med(media%sggMiEz(i1+1,j1,k1))%Is%SGBC) filo_placas=filo_placas+1
               if (SGG%Med(media%sggMiEz(i1-1,j1,k1))%Is%SGBC) filo_placas=filo_placas+1
               if (filo_placas < 2) then
                     es_unfilo_placa=.true.
               end if
!!!!!
               conta=conta+1
               compo => malon%Nodes(conta)
               compo%es_unfilo_placa = es_unfilo_placa
               SGBCdir=abs(SGG%Med(jmed)%Multiport(1)%Multiportdir)
               select case (SGBCdir)
                case (iEx)
                  compo%transversalDeltaE = 1.0_RKIND/IDxE(i1)
                  compo%transversalDeltaH = 1.0_RKIND/IDxh(i1)
                  compo%alignedlDeltaH     = 1.0_RKIND/IDyh(j1)
                  compo%med(1) =                  media%sggMiHy(i1  ,j1,k1)
                  compo%med(0) =                  media%sggMiHy(i1-1,j1,k1)
                  compo%Correct_Ha=.true.
                  compo%Correct_Hb=.false.
                case (iEy)
                  compo%transversalDeltaE = 1.0_RKIND/IDyE(j1)
                  compo%transversalDeltaH = 1.0_RKIND/IDyh(j1)
                  compo%alignedlDeltaH     = 1.0_RKIND/IDxh(i1)
                  compo%med(1) =                    media%sggMiHx(i1,j1  ,k1)
                  compo%med(0) =                    media%sggMiHx(i1,j1-1,k1)
                  compo%Correct_Ha=.false.
                  compo%Correct_Hb=.true.
                case DEFAULT
                  write(buff, *)    'Buggy ERROR: In SGBCs. '
                  call WarnErrReport (buff,.TRUE.)
               end select
               if ((compo%med(1) == 0).or.(compo%med(0) == 0).or.(sgg%med(compo%med(1))%Is%PEC).or.(sgg%med(compo%med(0))%Is%PEC)) then
                  write(buff, *)    'Buggy ERROR: In SGBCs. '
                  call WarnErrReport (buff,.TRUE.)
               end if
!!!!!!!!!
               compo%jmed  = jmed
               compo%Efield  => Ez(i1  ,j1  ,k1)
               compo%Ha_Plus => Hy(i1  ,j1  ,k1)
               compo%Ha_Minu => Hy(i1-1,j1  ,k1)
               compo%Hb_Plus => Hx(i1  ,j1  ,k1)
               compo%Hb_Minu => Hx(i1  ,j1-1,k1)
               call  depth(compo,sgg,jmed,SGBCFreq,SGBCresol,SGBCdepth)
               if (compo%depth==0) then
                    compo%SGBCCrank=.false.
               else if (compo%depth<0) then
                  write(buff, *)    'Buggy ERROR: In SGBCs compo%depth<0. ',compo%depth
                  call WarnErrReport (buff,.TRUE.)
               else
                    compo%SGBCCrank=SGBCcrank
               end if
               allocate (compo%E         (-compo%depth:compo%depth))
               if (compo%depth>0) allocate (compo%H         (-compo%depth:compo%depth-1)) 
               allocate(compo%E_past(-compo%depth:compo%depth))
               if (compo%SGBCcrank)  then
                    allocate(compo%d     (-compo%depth:compo%depth)) 
               end if

               compo%numpolres=malon%dispersiveMedia(compo%jmed)%numpolres
               if (SGBCDispersive) then  !I duplicate this info
                 allocate (compo%a11(1:compo%numpolres)) 
                 allocate (compo%c11(1:compo%numpolres)) 
                 compo%a11 = malon%dispersiveMedia(compo%jmed)%a11
                 compo%c11 = malon%dispersiveMedia(compo%jmed)%c11
                 allocate (compo%beta%val(1:compo%numpolres))
                 allocate (compo%kappa%val(1:compo%numpolres))
                 allocate (compo%G3%val(1:compo%numpolres))
                 !! call calc_g1g2(sgg,GM2,compo)     ! permit scal 071118
                 allocate (compo%EDis   (-compo%depth:compo%depth))
                 do ient=-compo%depth , compo%depth
                     allocate (compo%EDis(ient)%Current(1 : malon%dispersiveMedia(jmed)%numpolres))
                     compo%EDis(ient)%FieldPresent => compo%E(ient)
                 end do
               end if
            end if
         end do
      end do
   end do

   call calc_SGBCconstants(sgg,g,eps0,mu0,stochastic)

!!!depth reporting
    i=-100
    do conta=1,malon%numnodes
      compo => malon%Nodes(conta)
      if (compo%depth>i) then
         i=compo%depth
         jmed=compo%jmed
      end if
    end do

    write(buff, *)  ' Maximum SGBC depth= ',2*i,' for medium jmed= ',jmed
    call WarnErrReport (buff)



!!! I don't use it. it is not a rigorous stability criterion 311015
!!!      call test_stab(G2,GM2)

   !!!!!!!!!resuming
   if (.not.resume) then  
      do conta=1,malon%numnodes
         compo => malon%Nodes(conta)
         compo%E     =0.0_RKIND
         if (compo%SGBCcrank)  then 
             compo%E_past=0.0_RKIND
         end if
         compo%Hyee__left=0.0_RKIND
         compo%Hyee_right=0.0_RKIND
         compo%H     =0.0_RKIND
         !
         if (malon%SGBCDispersive) then 
             do i=-compo%depth,compo%depth
                compo%EDis(i)%fieldPresent=0.0_rkind
                compo%EDis(i)%fieldPrevious=0.0_rkind
                compo%EDis(i)%current=0.0_rkind
             end do
         end if
      end do
   else  
      do conta=1,malon%numnodes
         compo => malon%Nodes(conta)
         read (14) (compo%E     (i),i=-compo%depth,compo%depth)
         read (14) (compo%E_past(i),i=-compo%depth,compo%depth)
         read (14) compo%Hyee__left
         read (14) compo%Hyee_right
         read (14) (compo%H     (i),i=-compo%depth,compo%depth-1)
         !
         if (malon%SGBCDispersive) then 
             read(14) (compo%EDis(i)%fieldPrevious, i=-compo%depth,compo%depth)
             do k1=1,compo%NumPolRes
                read(14) (compo%EDis(i)%current(k1), i=-compo%depth,compo%depth)
             end do
         end if
      end do
   end if
   return

end subroutine InitSGBCs

subroutine calc_SGBCconstants(sgg,g,eps00,mu00,stochastic)
   real(kind=RKIND), intent(in) :: eps00,mu00
   type(SGGFDTDINFO_t), intent(in) :: sgg 
   real(kind=RKIND), pointer, dimension(:) :: gm1,g1,gm2,g2
   type(constants_t) :: g
 integer :: jmed,conta,i
   real(kind=RKIND) :: sigmam,sigma,mu,epsilonValue,signo,g1eff_0,g2eff_0,g1eff_1,g2eff_1
   type(SGBCSurface_t), pointer :: compo
   character(len=BUFSIZE) :: buFF
   logical :: stochastic
!
   g1 => g%g1; g2 => g%g2; gm1 => g%gm1; gm2 => g%gm2; 
   eps0=eps00; mu0=mu00; !hack to turn the step variables into globals
   zvac=sqrt(mu0/eps0)
   cluz=1.0_RKIND/sqrt(mu0*eps0)
!!!I allocate all the constant matrices
 do conta=1,malon%numnodes
     compo => malon%Nodes(conta)  
     !!!I call depth again so it recalculates deltaentreEinterno correctly !110523 needed for stochastic
     jmed=compo%jmed
     call depth(compo,sgg,jmed,SGBCFreq,SGBCresol,SGBCdepth)
               
     if (.not.allocated(compo%GM1_interno)) then !!!there are spare extremes of g and gm but I leave it like this to avoid tinkering more 0121
         allocate(compo%GM1_interno (-compo%depth  :compo%depth-1) ,&   !the extremes are adjusted correctly        
                  compo%GM2_interno (-compo%depth  :compo%depth-1) ,&
                  compo%G1_interno  (-compo%depth+1:compo%depth-1) ,&
                  compo%G2_interno  (-compo%depth+1:compo%depth-1) ,&            
                  compo%a           (-compo%depth  :compo%depth) ,& !one more is added at the beginning and end but I don't touch it 0121
                  compo%b           (-compo%depth  :compo%depth) ,&
                  compo%c           (-compo%depth  :compo%depth) ,&
                  compo%rb          (-compo%depth  :compo%depth) ,&
                  compo%rh          (-compo%depth  :compo%depth) ,&
                  compo%rhm1        (-compo%depth  :compo%depth))
     end if

 end do
 
 !absurd defaults
 compo%GM1_interno = -1e23 
 compo%G1_interno  = +2e24
 compo%GM2_interno = -1e26 
 compo%G2_interno  = +2e26
 compo%a  =+1.3e24        
 compo%b  =+2.3e24        
 compo%c  =-1.3e24        
 compo%rb =+1.3e24        
 compo%rh =+4.3e24     
 compo%rhm1 =+2.7e25 
 
!!!!CALCULATION OF THE YEE-FDTD COEFFICIENTS (and seed of the CN-FDTD)

   if (SGBCDispersive) then
       do jmed=1,sgg%NumMedia
          if ((SGG%Med(jmed)%Is%SGBCDispersive).and.(.not.(SGG%Med(jmed)%Is%PML))) then   
!!!071118 for permit scaling
!!!070717  I recalculate and overwrite G1,G2,GM1, and GM2 with those appearing in the dispersive file although I think it is not used at all
               Sigmam  =      SGG%Med(jmed)%multiport(1)%sigmam(1)
               epsilonValue = Eps0*SGG%Med(jmed)%multiport(1)%epr(1)
               Mu      = Mu0* SGG%Med(jmed)%multiport(1)%mur(1)
               Sigma   =      SGG%Med(jmed)%multiport(1)%sigma(1)
               G1(jmed)=(1 -  Sigma * sgg%dt / (2.0_RKIND * epsilonValue ) ) / (1.0_RKIND + Sigma * sgg%dt / (2.0_RKIND * epsilonValue ))
               G2(jmed)=sgg%dt /epsilonValue                        / (1.0_RKIND + Sigma * sgg%dt / (2.0_RKIND * epsilonValue ))
               if (g1(jmed) < 0.0_RKIND) then !exponential time stepping
                  g1(jmed)=exp(- Sigma * sgg%dt / (epsilonValue ))
                  g2(jmed)=(1.0_RKIND-g1(jmed))/ Sigma
               end if
               GM1(jmed)=(1- SigmaM*sgg%dt/(2.0_RKIND *  Mu )) /(1.0_RKIND + SigmaM*sgg%dt/(2.0_RKIND *  Mu ))
               GM2(jmed)=sgg%dt/ Mu                   /(1.0_RKIND + SigmaM*sgg%dt/(2.0_RKIND *  Mu ))
               if (gm1(jmed) < 0.0_RKIND) then !exponential time stepping
                  gm1(jmed)=exp(- Sigmam * sgg%dt / (Mu ))
                  gm2(jmed)=(1.0_RKIND-gm1(jmed))/ Sigmam
               end if
!!!!end 070717
          end if
      end do
   end if    

#ifdef CompileWithOpenMP
!$OMP  PARALLEL do  DEFAULT(none) private (compo,buff) shared(malon,sgg,eps00,mu00,GM2,SGBCDispersive)
#endif
 do conta=1,malon%numnodes
     compo => malon%Nodes(conta)
     call calc_g1g2gm1gm2_compo(sgg,compo,eps00,mu00,SGBCDispersive)
!update constant for the H external to the multilayer!!! 0121
     compo%Gm2_externo=Gm2(compo%jmed) / compo%transversalDeltaE !note there was compo%transversalDeltaE before sgg 130516 ! but I think it is  compo%transversalDeltaH! 0121 No. it is deltaE because it is used to update the external H
 end do

#ifdef CompileWithOpenMP
!$OMP  END PARALLEL DO
#endif
!
!

#ifdef CompileWithOpenMP
!$OMP  PARALLEL do  DEFAULT(none) private (compo,signo,g1eff_0,g1eff_1,g2eff_0,g2eff_1) shared(malon,gm1,gm2)
#endif
    do conta=1,malon%numnodes
      compo => malon%Nodes(conta) 
      if (compo%depth>0) then !one of those in the curl
          if (compo%Correct_Ha) then
             signo=+1.0_RKIND
             g1eff_0=   compo%G1(0)  
             g1eff_1=   compo%G1(1)   
             g2eff_0=   signo *compo%G2a(0)  
             g2eff_1=   signo *compo%G2a(1)
          else if (compo%Correct_Hb) then
             signo=-1.0_RKIND
             g1eff_0=   compo%G1(0)  
             g1eff_1=   compo%G1(1)   
             g2eff_0=   signo *compo%G2b(0)  
             g2eff_1=   signo *compo%G2b(1)
          end if
!   
!
          !
          do i=-compo%depth , compo%depth-1
              compo%GM2_interno (i) = signo * compo%Gm2_interno (i) 
          end do
          do i=-compo%depth+1 , compo%depth-1
              compo%G2_interno  (i) = signo * compo%G2_interno (i) 
          end do
!!!Crank               
          do i=-compo%depth+1 , compo%depth-1 !the first and last are not used, but I leave them to not alter the algorithm pre 0121
!                 compo%a  (i)        =                      - compo%G2_interno (i  ) * compo%GM2_interno (i  ) /4.0_RKIND
               compo%a  (i)        =                      - compo%G2_interno (i) * compo%GM2_interno (i-1) /4.0_RKIND   !!!!the minus 1
               
!                 compo%b  (i)        = 1.0_RKIND            + compo%G2_interno (i  ) * compo%GM2_interno (i  ) /2.0_RKIND   
               compo%b  (i)        = 1.0_RKIND            + compo%G2_interno (i) * compo%GM2_interno (i-1) /4.0_RKIND  + compo%G2_interno (i) * compo%GM2_interno (i) /4.0_RKIND   
               
!                 compo%c  (i)        =                      - compo%G2_interno (i  ) * compo%GM2_interno (i  ) /4.0_RKIND     
               compo%c  (i)        =                      - compo%G2_interno (i) * compo%GM2_interno (i) /4.0_RKIND   !same because the +1 is indexed in i
               
!                 compo%rb (i)        = compo%G1_interno (i) - compo%G2_interno (i  ) * compo%GM2_interno (i  ) /2.0_RKIND
               compo%rb (i)        = compo%G1_interno (i) - compo%G2_interno (i) * compo%GM2_interno (i-1) /4.0_RKIND  - compo%G2_interno (i) * compo%GM2_interno (i) /4.0_RKIND   !!the minus 1
               
!                 compo%rh  (i)        =(compo%G2_interno (i  ) * compo%GM1_interno(i  ) + compo%G2_interno  (i  ))/2.0_RKIND
               compo%rh  (i)        =(compo%G2_interno (i) * compo%GM1_interno(i) + compo%G2_interno  (i))/2.0_RKIND  !it must be broken down
               compo%rhm1(i)        =(compo%G2_interno (i) * compo%GM1_interno(i-1) + compo%G2_interno  (i))/2.0_RKIND
!!!        
          end do
          i=-compo%depth
          compo%a1          =  0.0
          compo%c1          =                 - g2eff_0 * compo%GM2_interno (i) /4.0_RKIND !the one to its right, internal, together with the yee
          compo%b1           = 1.0_RKIND      + g2eff_0 * compo%GM2_interno (i) /4.0_RKIND
          compo%rb1          = g1eff_0 - g2eff_0        * compo%GM2_interno (i) /4.0_RKIND
          compo%rh1          =  (g2eff_0                * compo%GM1_interno (i)+ g2eff_0)/2.0_RKIND
          i=compo%depth
          compo%cn          =  0.0
          compo%an          =                 - g2eff_1 * compo%GM2_interno (i-1) /4.0_RKIND !the one to its left, internal, together with the yee
          compo%bn           = 1.0_RKIND      + g2eff_1 * compo%GM2_interno (i-1) /4.0_RKIND
          compo%rbn          = g1eff_1 - g2eff_1        * compo%GM2_interno (i-1) /4.0_RKIND
          compo%rhn          =  (g2eff_1                * compo%GM1_interno (i-1) + g2eff_1)/2.0_RKIND
!!!!end of new formulation
!!pre-jav 310116
!!!                 compo%a1=0.0; compo%c1=0.0; compo%b1=1.0; compo%an=0.0; compo%cn=0.0; compo%bn=1.0; compo%rb1=0.0; compo%rbn=0.0; compo%rh1=0.0; compo%rhn=0.0    
!!!!post-jav new formulation jav 310116       
      end if
    end do 
#ifdef CompileWithOpenMP
!$OMP  END PARALLEL DO
#endif
    return
end subroutine calc_SGBCconstants


                        
 subroutine YeeAdvanceSGBCDispersive (tempnode,numpolres,G3,kappa,beta,dt)
 
   real(kind=RKIND), intent(in) :: dt
   type(MDfield_t), pointer  :: tempnode
   integer(kind=4) :: k1,numpolres
   type(val_t) :: Beta,Kappa,G3
   
         do k1=1,NumPolRes
            tempnode%fieldPresent=tempnode%FieldPresent-real(G3%val(k1)*tempnode%current(k1))
         end do
         do k1=1,NumPolRes
            tempnode%current(k1)= Kappa%val(k1)  *tempnode%current(k1) + &
                                  Beta%val(k1)*(tempnode%fieldPresent-tempnode%fieldPrevious) /dt
         end do
         tempnode%fieldPrevious=tempnode%fieldPresent
         !stores previous field (careful, it is not a pointer but a value assignment)
         !before the background algorithm starts calculating it again

 end subroutine YeeAdvanceSGBCDispersive
 
                         
 subroutine primero_CNAdvanceSGBCDispersive (tempnode,tempD,numpolres,G3,kappa,beta,dt)
 
   real(kind=RKIND), intent(in) :: dt
   real(kind=RKIND),  pointer  :: tempD
   type(MDfield_t), pointer  :: tempnode
   integer(kind=4) :: k1,numpolres
   type(val_t) :: Beta,Kappa,G3
         do k1=1,NumPolRes
            tempD=tempD-real(G3%val(k1)*tempnode%current(k1))
         end do
 end subroutine primero_CNAdvanceSGBCDispersive
 
 subroutine segundo_CNAdvanceSGBCDispersive (campocalculado,tempnode,numpolres,G3,kappa,beta,dt)
 
   real(kind=RKIND) :: campocalculado 
   real(kind=RKIND), intent(in) :: dt
   type(MDfield_t), pointer  :: tempnode
   integer(kind=4) :: k1,numpolres
   type(val_t) :: Beta,Kappa,G3

         do k1=1,NumPolRes
            tempnode%fieldPresent=campocalculado
            tempnode%current(k1)= Kappa%val(k1)  *tempnode%current(k1) + &
                                  Beta%val(k1)*(tempnode%fieldPresent-tempnode%fieldPrevious) /dt
         end do
         tempnode%fieldPrevious=tempnode%fieldPresent
         !stores previous field (careful, it is not a pointer but a value assignment)
         !before the background algorithm starts calculating it again

 end subroutine segundo_CNAdvanceSGBCDispersive
                        

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine to advance the E field in the SGBC: Usual Yee
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine AdvanceSGBCE(dt,SGBCDispersive,simu_devia,stochastic)
   logical :: simu_devia,stochastic 
   real(kind=RKIND), intent(in) :: dt
   logical :: SGBCDispersive
   integer(kind=4) :: conta
!       call AdvanceSGBCE_single_node(dt,SGBCDispersive)
#ifdef CompileWithOpenMP
!$OMP  PARALLEL do DEFAULT(SHARED) private (conta) schedule(guided) 
#endif
do conta=1,malon%numnodes
   !call AdvanceSGBCE_single_node(malon%Nodes(conta), dt,SGBCDispersive)
   call AdvanceSGBCE_single_node(conta, dt,SGBCDispersive)
end do
#ifdef CompileWithOpenMP
!$OMP  END PARALLEL DO
#endif

contains
   !subroutine AdvanceSGBCE_single_node(compo,dt,SGBCDispersive) !this one would also be valid
      !type(SGBCSurface_t),target   :: compo
   subroutine AdvanceSGBCE_single_node(conta,dt,SGBCDispersive)
      !non argument arguments
      integer(kind=4),  intent(in) :: conta
      real(kind=RKIND), intent(in) :: dt
      logical,           intent(in) :: SGBCDispersive
      !local
      integer(kind=4) :: i
      type(SGBCSurface_t),pointer   :: compo
      type(MDfield_t) , pointer :: EDIS
      real(kind=RKIND), pointer :: dDIS
      

      
      !do conta=1,malon%numnodes
      compo => malon%Nodes(conta)

!!!the extremes of the internal E fields
         if (compo%depth>0) then
!this only does something in the yee case.
            if (.not.compo%SGBCcrank)  then !yee
               if (compo%Correct_Ha) then
!!!!THE EXTREMES
!!!!note THE filo_placaS ARE NOT CORRECTED WITH FDTD IN DISPERSIVE CASES. SO IT IS TO BE EXPECTED THAT IT IS UNSTABLE WITH DISPERSIVE SGBCYEEE
                  compo%E(compo%depth) = compo%G1(1) *compo%E(compo%depth) +  &
                                       (compo%G2a(1) *(compo%Ha_Plus   - compo%Hyee_Right) - compo%G2b(1) *(compo%Hb_Plus   - compo%Hb_Minu  ) )
                  !
                  compo%E(-compo%depth) = compo%G1(0) *compo%E(-compo%depth) +  &
                                          (compo%G2a(0) *(compo%Hyee__Left - compo%Ha_Minu  ) - compo%G2b(0) *(compo%Hb_Plus   - compo%Hb_Minu  ) )

               else if (compo%Correct_Hb) then
                  compo%E(compo%depth) = compo%G1(1) *compo%E(compo%depth) +  &
                                          (compo%G2a(1) *(compo%Ha_Plus   - compo%Ha_Minu  ) - compo%G2b(1) *(compo%Hb_Plus   - compo%Hyee_Right ) )
                  compo%E(-compo%depth) = compo%G1(0) *compo%E(-compo%depth) +  &
                                          (compo%G2a(0) *(compo%Ha_Plus   - compo%Ha_Minu  ) - compo%G2b(0) *(compo%Hyee__Left - compo%Hb_Minu ) )
                     
               end if
            end if !OF MALONECRANK
         else !if it has depth=0
            compo%E(compo%depth) = compo%G1(0) *compo%E(compo%depth) +  (compo%G2a(0) *(compo%Ha_Plus   - compo%Ha_Minu    ) - compo%G2b(0) *(compo%Hb_Plus     - compo%Hb_Minu  ) )
         end if !OF COMODEPTH
!!!the FDTD1D internal E fields
         if (compo%SGBCcrank)  then !below 2 it makes no sense
            do i=-compo%depth , compo%depth
               compo%E_past (i) = compo%E(i)
            end do
!this is not necessary. it is already passed via MPI
!!!pre-jav 310116 (only these lines
            !compo%d( compo%depth  )      = compo%E( compo%depth  )
            !compo%d(-compo%depth  )      = compo%E(-compo%depth  )
!!post-jav 310116 boundaries 
!!!!note these coefficients would not work well with magnetic materials because I do not take into account the gm1 (in jav's crank-nicolson) 
            if (compo%Correct_Ha) then !!!compo%G2a(0) is equal to compo%rh1 !!!compo%G2a(1) is equal to compo%rhn
               i=compo%depth
               compo%d(i) =    - compo%an * (compo%E_past(i-1)) + compo%rbn * compo%E_past(i) + &
                                          (compo%G2a(1) *(compo%Ha_Plus   - compo%Hyee_Right) - compo%G2b(1) *(compo%Hb_Plus   - compo%Hb_Minu  ) )
               !
               !
               i=-compo%depth
               compo%d(i) =    - compo%c1 * (compo%E_past(i+1)) + compo%rb1 * compo%E_past(i) + &
                                          (compo%G2a(0) *(compo%Hyee__Left - compo%Ha_Minu  ) - compo%G2b(0) *(compo%Hb_Plus   - compo%Hb_Minu  ) )
               !
               !
            else if (compo%Correct_Hb) then !!!-compo%G2b(0) is equal to compo%rh1 !!!-compo%G2b(1) is equal to compo%rhn
               i=compo%depth
               compo%d(i) =   - compo%an * (compo%E_past(i-1)) + compo%rbn * compo%E_past(i) + &
                                          (compo%G2a(1) *(compo%Ha_Plus   - compo%Ha_Minu  ) - compo%G2b(1) *(compo%Hb_Plus   - compo%Hyee_Right ) )
               i=-compo%depth
               compo%d(i) =    - compo%c1 * (compo%E_past(i+1)) + compo%rb1 * compo%E_past(i) + &
                    (compo%G2a(0) *(compo%Ha_Plus   - compo%Ha_Minu  ) - compo%G2b(0) *(compo%Hyee__Left - compo%Hb_Minu ) )
            end if
!end of new boundaries jav 310116 (must match by commenting the previous and leaving the pre-jav
            do i=-compo%depth+1 , compo%depth-1
!!!                  compo%d(i) = - compo%a(i) * (compo%E_past(i-1)+compo%E_past(i+1)) &!!!+ compo%rb(i) * compo%E_past(i ) + compo%rh(i) * (compo%H(i)-compo%H(i-1)) !! compo%a is equal to compo%c
!!!!   mmmmmmm it is not true that a = c for multilayer because it is asymmetric 0121: I fix the above roughly
                compo%d(i)=-compo%a(i)*compo%E_past(i-1) &
                           -compo%c(i)*compo%E_past(i+1) &
                           +compo%rb(i)*compo%E_past(i)  &
                           +compo%rh  (i)*compo%H(i)   & !!roughly!!!0121
                           -compo%rhm1(i)*compo%H(i-1)    !!roughly!!!0121
            end do
            
            !
            if (SGBCDispersive) then 
               do i=-compo%depth + 1 , compo%depth -1 !without the damn filo_placas
               EDIS=>compo%EDis(i)
               dDIS=>compo%d(i)
               call primero_CNAdvanceSGBCDispersive (EDIS,dDIS,compo%numpolres,compo%G3,Compo%kappa,compo%beta,dt)
               end do
            end if
            call solve_tridiag_distintos(compo%a ,compo%b ,compo%c , &
                                       compo%a1,compo%b1,compo%c1, &
                                       compo%an,compo%bn,compo%cn, &
                                       compo%d,compo%E,2*compo%depth+1)
            if (SGBCDispersive) then 
               do i=-compo%depth + 1 , compo%depth -1 !without the damn filo_placas
               EDIS=>compo%EDis(i)
               call segundo_CNAdvanceSGBCDispersive (compo%E(i),EDIS,compo%numpolres,compo%G3,Compo%kappa,compo%beta,dt)
               end do
            end if
         else !YEE
            do i=-compo%depth+1 , compo%depth-1
               compo%E(i) = compo%G1_interno (i)  *compo%E(i) + compo%G2_interno (i)  *( compo%H(i) - compo%H(i-1) )               
               ! THE DAMN filo_placaS WERE NOT CORRECTED BEFORE. HERE ONLY THE INTERIOR IS CORRECTED
               if (SGBCDispersive) then 
                  EDIS=>compo%EDis(i)
                  call YeeAdvanceSGBCDispersive (EDIS,compo%numpolres,compo%G3,Compo%kappa,compo%beta,dt)
               end if
            end do
                                    
         end if
!!!the 1D internal H fields 
         if (compo%SGBCcrank)  then
            do i=-compo%depth , compo%depth-1
               compo%H(i) = compo%GM1_interno(i) *compo%H(i) + compo%GM2_interno(i)   /2.0_RKIND *( compo%E(i+1)      - compo%E(i)      + &
                                                                                               compo%E_past(i+1) - compo%E_past(i) )
            end do                 
!for crank-nicolson a half-step advance is necessary for the H field since E and H are synchronous. Therefore two H fields coexist at (-compo%depth) and (compo%depth-1): one at n and another at n+1/2
!!pre-jav only 310116 (only these lines
!!!!this correction is analytically correct but a source of possible instabilities since it is a Yee. I comment it out 02115 and leave a backwards approx which is more than enough for metals where the speed is very low compared to vacuum
!            compo%Hyee__Left = compo%GM1_interno *compo%Hyee__Left + compo%GM2_interno *( compo%E(-compo%depth+1) - compo%E(-compo%depth  ) )
!            compo%Hyee_Right = compo%GM1_interno *compo%Hyee_Right + compo%GM2_interno *( compo%E(   compo%depth) - compo%E( compo%depth-1) )
!!!!post-jav 310116 boundaries only what follows is correct. It is no longer a TB approximation but the result of the new CN_jav
            compo%Hyee__Left = compo%H(-compo%depth) 
            compo%Hyee_Right = compo%H(compo%depth-1)              
            
         else !yee
            if (compo%depth/=0) then
               do i=-compo%depth , compo%depth-1
                  compo%H(i) = compo%GM1_interno(i) *compo%H(i) +  compo%GM2_interno(i) *( compo%E(i+1) - compo%E(i) )  !E and H are always used reciprocally with the same sign
               end do
               compo%Hyee__Left = compo%H(-compo%depth)
               compo%Hyee_Right = compo%H(compo%depth-1)
            end if
         end if
!!!!I copy the average into its Efield for probe request purposes and calculation at the borders (then when advancing H the main one on both sides will use this Efield, but advanceSGBCH will correct it with the correct one
         compo%Efield =(compo%E(-compo%depth)+compo%E(compo%depth))/2. 
      !end do
   end subroutine AdvanceSGBCE_single_node
end subroutine AdvanceSGBCE

subroutine AdvanceSGBCH
   integer(kind=4) :: conta
   type(SGBCSurface_t), pointer :: compo
   character(len=BUFSIZE) :: buFF
   !NOTE: This cannot be optimized 
   !      because two or more compo%H{a,b} can point
   !      to the same field an access conflict can occur
   do conta=1,malon%numnodes
      compo => malon%Nodes(conta)
!!!!note: it is a correction to what the main one does using the correct electric field
      if (compo%Correct_Ha) then
         compo%Ha_Plus = compo%Ha_Plus +  compo%gm2_externo* (compo%Efield - compo%E(compo%depth)) !I insist: it is a correction: the main one has added/removed Efield and must remove/add E of the corresponding extreme
         compo%Ha_Minu = compo%Ha_Minu -  compo%gm2_externo* (compo%Efield - compo%E(-compo%depth))
      else if (compo%Correct_Hb) then                                                       
         compo%Hb_Plus = compo%Hb_Plus -  compo%gm2_externo* (compo%Efield - compo%E(compo%depth))
         compo%Hb_Minu = compo%Hb_Minu +  compo%gm2_externo* (compo%Efield - compo%E(-compo%depth))
      else     
         write(buff, *)    'Buggy ERROR: In SGBCs. '
         call StopOnError (0,0,buff)
      end if
   end do
   return
end subroutine AdvanceSGBCH

subroutine calc_g1g2gm1gm2_compo(sgg,compo,eps00,mu00,SGBCDispersive)
   real(kind=RKIND), intent(in) :: eps00,mu00
   complex(kind=CKIND), pointer, dimension(:) :: Beta,Kappa,G3
!!!!      
   type(SGGFDTDINFO_t), intent(in) :: sgg
   type(SGBCSurface_t), pointer, intent(inout) :: compo
   character(len=BUFSIZE) :: buff
!!!local variables
 real(kind=RKIND) :: width,sigmatemp,eprtemp,sigmamtemp,murtemp,epsilonValue,sigma,mu,sigmam,g1,g2,gm1,gm2,delta_entreEinterno_temp,epr_adyacentei,sig_adyacentei
 real(kind=RKIND), dimension(0:1) :: epr_adyacente,sig_adyacente
 integer(kind=4) :: i,ib,ib_ady
 logical :: SGBCDispersive
 eps0=eps00; mu0=mu00; !hack to turn the step variables into globals
   


   if (compo%depth==0) then
         !!!!!!! compo%delta_entreeinterno=0.0 !!never used !note check case 0121 
         !!!!!!!averagefactor  = width / compo%transversaldeltah /factor !!! sgg well-averaged filo_placas 201115
         !!!!!!!epsilon = (1.0_rkind - averagefactor ) * (epr_adyacente(0)+epr_adyacente(1))/2.0_rkind  * eps0 + &
         !!!!!!!                       averagefactor   * eprtemp           * eps0
         !!!!!!!sigma =   (1.0_rkind - averagefactor ) * (sig_adyacente(0)+sig_adyacente(1))/2.0_rkind  + &
         !!!!!!!                       averagefactor   * sigmatemp
         compo%delta_entreEinterno=0.0 !!never used
         do i=0,1
 !0121 I am going to take vacuum BECAUSE AT THE CORNERS BETWEEN sgbc IT DETECTS THE ADJACENT MEDIUM WRONG. IT WAS DONE BEFORE 0121 LIKE THIS TOO
 !050421 I return it to vacuum because I cannot quite see the case of the corners between SGBC
           !  epr_adyacente(i) = Sgg%Med(compo%med(i))%epr   
           !  sig_adyacente(i) = Sgg%Med(compo%med(i))%sigma 
           !  if (((epr_adyacente(i)-1.0_RKIND>1E-3).OR.(sig_adyacente(i)>1E-3)).AND.(.NOT.(SGG%Med(compo%med(i))%Is%SGBC)) )  then 
           !           write(buff, *)    '(WARNING) Collision of composite with non free-space medium. Assuming free-space, instead ', compo%med(i), epr_adyacente(i), sig_adyacente(i)
                      epr_adyacente(i) = 1.0_RKIND
                      sig_adyacente(i) = 0.0_RKIND
           !           call WarnErrReport (buff,.false.)
           !  end if
         end do
         width=sgg%med(compo%jmed)%Multiport(1)%width(1)
         sigmatemp=sgg%Med(compo%jmed)%multiport(1)%sigma(1)
         eprtemp= sgg%Med(compo%jmed)%multiport(1)%epr(1)   
       !I do without the filo_placa 0121 because at the PEC boundaries it detects them incorrectly !anyway I have never liked this 0121
         !I cannot do without the filo_placas as of 040523 SinSTOCH_antiguou_th0.0001
         if (compo%es_unfilo_placa) then
             epsilonValue = ((epr_adyacente(0)+epr_adyacente(1))/2.0_rkind  * eps0 *(compo%transversaldeltah - width/2.0_rkind)   + &
                         eprtemp                                       * eps0 * width /2.0_rkind) / &
                        (compo%transversaldeltah)
             sigma =   ((sig_adyacente(0)+sig_adyacente(1))/2.0_rkind         *(compo%transversaldeltah - width/2.0_rkind)   + &
                         sigmatemp                                            *width  /2.0_rkind) / &
                        (compo%transversaldeltah)
         else
         epsilonValue = ((epr_adyacente(0)+epr_adyacente(1))/2.0_rkind  * eps0 *(compo%transversaldeltah-width)   + &
                     eprtemp                                       * eps0 *width             ) / &
                        (compo%transversaldeltah)
         sigma =   ((sig_adyacente(0)+sig_adyacente(1))/2.0_rkind         *(compo%transversaldeltah-width)   + &
                     sigmatemp                                            *width             ) / &
                        (compo%transversaldeltah)
         end if
!!!!!first adjusts the g1 and g2 of the layer thickness edges !note in sgbcdispersive the constants kappa, beta, g3 are not used in the damn filo_placas. only in the interior
         call g1g2(sgg%dt,epsilonValue,sigma,g1,g2)
         compo%g1   (0)=g1 
         if (compo%correct_ha) then
             compo%g2a(0)= g2 / compo%transversaldeltah 
             compo%g2b(0)= g2 / compo%alignedldeltah
 !!!!the indices (1) are not used for the particular case depth=0
             compo%g2a(1)= compo%g2a(0)
             compo%g2b(1)= compo%g2b(0)
         else if (compo%correct_hb) then
             compo%g2a(0)= g2 / compo%alignedldeltah
             compo%g2b(0)= g2 / compo%transversaldeltah
             compo%g2a(1)= compo%g2a(0)
             compo%g2b(1)= compo%g2b(0)
         end if
!!!!!!!!!!!!now the interior  
        epsilonValue =  eprtemp * eps0
        sigma =    sigmatemp 
        if (sgbcdispersive) then
             beta => compo%beta %val
             kappa => compo%kappa %val
             g3 => compo%g3 %val
             call g1g2_dispersive(sgg%dt,epsilonValue,sigma,g1,g2,beta,kappa,g3,compo%numpolres,compo%a11,compo%c11)
         else
             call g1g2(sgg%dt,epsilonValue,sigma,g1,g2)
         end if
         compo%g1_interno(0)=g1
         compo%g2_interno(0)=g2
         sigmamtemp= 0.
         murtemp=    1.
         mu=murtemp * mu0
         Sigmam = sigmamtemp
         call gm1gm2(sgg%dt,mu,sigmam,gm1,gm2)
         compo%gm1_interno(0)=gm1
         compo%gm2_interno(0)=gm2
 else !OF MALONYEDEPTH
!FIRST the constants at the boundaries in the thickness dimension
     do i=0,1
 !0121 I am going to take vacuum BECAUSE AT THE CORNERS BETWEEN sgbc IT DETECTS THE ADJACENT MEDIUM WRONG. IT WAS DONE BEFORE 0121 LIKE THIS TOO
 !050421 I return it to vacuum because I cannot quite see the case of the corners between SGBC
         ! epr_adyacente(i) = Sgg%Med(compo%med(i))%epr   
         ! sig_adyacente(i) = Sgg%Med(compo%med(i))%sigma 
         ! if (((epr_adyacente(i)-1.0_RKIND>1E-3).OR.(sig_adyacente(i)>1E-3)).AND.(.NOT.(SGG%Med(compo%med(i))%Is%SGBC)) )  then 
         !         write(buff, *)    '(WARNING) Collision of composite with non free-space medium. Assuming free-space, instead ', compo%med(i), epr_adyacente(i), sig_adyacente(i)
                  epr_adyacente(i) = 1.0_RKIND
                  sig_adyacente(i) = 0.0_RKIND
         !         call WarnErrReport (buff,.false.)
         ! end if
         if (i==0) then 
             ib=1 !first layer
             delta_entreEinterno_temp=compo%delta_entreEinterno(-compo%depth)
         else
             ib=sgg%Med(compo%jmed)%multiport(1)%numLayers !last layer
             delta_entreEinterno_temp=compo%delta_entreEinterno(compo%depth-1)     
         end if
         width=sgg%med(compo%jmed)%Multiport(1)%width(ib)
         sigmatemp= sgg%Med(compo%jmed)%multiport(1)%sigma(ib)
         eprtemp=   sgg%Med(compo%jmed)%multiport(1)%epr(ib)  
       !I do without the filo_placa 0121 because at the PEC boundaries it detects them incorrectly !anyway I have never liked this 0121
                     !I cannot do without the filo_placas as of 040523 SinSTOCH_antiguou_th0.0001
         if (compo%es_unfilo_placa) then
             epsilonValue = (epr_adyacente(i)* eps0 *(compo%transversalDeltaH + delta_entreEinterno_temp /2.0_RKIND)   + &
                                eprtemp         * eps0 *               (delta_entreEinterno_temp /2.0_RKIND)) / &
                               (compo%transversalDeltaH +               delta_entreEinterno_temp)
             Sigma =   (sig_adyacente(i)       *(compo%transversalDeltaH + delta_entreEinterno_temp /2.0_RKIND)   + &
                                sigmatemp               *(delta_entreEinterno_temp  /2.0_RKIND)) / &
                               (compo%transversalDeltaH + delta_entreEinterno_temp)
         else
             epsilonValue = (epr_adyacente(i)* eps0 *compo%transversalDeltaH   + &
                        eprtemp         * eps0 *delta_entreEinterno_temp    ) / &
                       (compo%transversalDeltaH + delta_entreEinterno_temp)                
             Sigma =   (sig_adyacente(i)       *compo%transversalDeltaH   + &
                        sigmatemp              *delta_entreEinterno_temp    ) / &
                       (compo%transversalDeltaH + delta_entreEinterno_temp)
             
         end if

       
         
         
 !first adjusts the g1 AND G2      !I do not need gm1 or gm2 in the filo_placas
         call g1g2(sgg%dt,epsilonValue,sigma,g1,g2)
         compo%g1(i)=g1 
         if (compo%Correct_Ha) then
             compo%G2a(i)= G2 / (0.5_RKIND * compo%transversalDeltaH + 0.5_RKIND*delta_entreEinterno_temp)
             compo%G2b(i)= G2 / compo%alignedlDeltaH
         else if (compo%Correct_Hb) then
             compo%G2a(i)= G2 / compo%alignedlDeltaH
             compo%G2b(i)= G2 / (0.5_RKIND * compo%transversalDeltaH + 0.5_RKIND*delta_entreEinterno_temp)
         end if
     end do !OF THE SWEEP i=0,1 OF the two sheet boundaries in the thickness dimension
     
 !now the interior !The first and last G are not used 0121
     compo%G2_interno =2e31 !absurd default values to detect errors
     compo%G1_interno =-2e21 
     barridoporcapas: do i=-compo%depth+1,compo%depth-1   !0121
         ib=compo%layerIndex(i)         
         ib_ady=compo%layerIndex(i-1)         
         if ((ib<1).or.(ib>sgg%Med(compo%jmed)%multiport(1)%numLayers)) then
             write(buff, *)   'Buggy error in ib fuera de rango en compo numcapas. Contact '
             call StopOnError (0,0,buff)
             stop
         end if
         eprtemp=    sgg%Med(compo%jmed)%multiport(1)%epr(ib)  
         sigmatemp=  sgg%Med(compo%jmed)%multiport(1)%sigma(ib)
         epr_adyacentei = sgg%Med(compo%jmed)%multiport(1)%epr(ib_ady)
         sig_adyacentei = sgg%Med(compo%jmed)%multiport(1)%sigma(ib_ady)
         !!! Interpolatory average done well 0121
         eprtemp = (epr_adyacentei     * compo%delta_entreEinterno(i-1)   + &
                    eprtemp            * compo%delta_entreEinterno(i) ) / &
                                        (compo%delta_entreEinterno(i-1) + compo%delta_entreEinterno(i))
         sigmatemp =   (sig_adyacentei  *compo%delta_entreEinterno(i-1)   + &
                        sigmatemp       *compo%delta_entreEinterno(i) ) / &
                                        (compo%delta_entreEinterno(i-1) + compo%delta_entreEinterno(i))
         !!!
         epsilonValue =  eprtemp * eps0
         Sigma =    sigmatemp 
         call g1g2(sgg%dt,epsilonValue,sigma,g1,g2)
         compo%g1_interno(i)=g1
         compo%g2_interno(i)=g2 /((compo%delta_entreEinterno (i)+compo%delta_entreEinterno (i-1))/2.0_RKIND) !half-sum diff  not centered between layers
     end do barridoporcapas
     
     compo%GM2_interno=-1e30 !absurd default values to detect errors
     compo%GM1_interno=3e22 
 !now the interior !the last GM which is not used 0121
     barridoporcapasH: do i=-compo%depth,compo%depth-1   !0121
         ib=compo%layerIndex(i)         
         if ((ib<1).or.(ib>sgg%Med(compo%jmed)%multiport(1)%numLayers)) then
             write(buff, *)   'Buggy error in ib fuera de rango en compo numcapas. '
             call StopOnError (0,0,buff)
             stop
         end if
         sigmamtemp= sgg%Med(compo%jmed)%multiport(1)%sigmam(ib)
         murtemp=    sgg%Med(compo%jmed)%multiport(1)%mur(ib)
         mu=murtemp * mu0
         Sigmam = sigmamtemp
         call gm1gm2(sgg%dt,mu,sigmam,gm1,gm2)
         compo%gm1_interno(i)=gm1
         compo%gm2_interno(i)=gm2  /compo%delta_entreEinterno (i) !there is no need to half-sum because it is internal
     end do barridoporcapasH
     
     
 end if !of compodepth

   return
end subroutine calc_g1g2gm1gm2_compo

!!!!!!!
subroutine g1g2(dt,epsilonValue,sigma,G1,G2)
   real(kind=RKIND_TIME), intent(in) :: dt
   real(kind=RKIND), intent(in) :: epsilonValue,sigma
   real(kind=RKIND), intent(out) :: g1,g2

   G1=(1.0_RKIND  - Sigma * dt / (2.0_RKIND * epsilonValue ) ) / &
      (1.0_RKIND  + Sigma * dt / (2.0_RKIND * epsilonValue ) )
   G2=dt / epsilonValue                       / &
      (1.0_RKIND  + Sigma * dt / (2.0_RKIND * epsilonValue))

   if (g1 < 0.0_RKIND) then !exponential time stepping
      g1=exp(- Sigma * dt / (epsilonValue ))
      g2=(1.0_RKIND-g1)/ Sigma
   else
      continue
   end if   
   return
end subroutine g1g2

!!!!!!!
subroutine gm1gm2(dt,mu,sigmam,Gm1,Gm2)
   real(kind=RKIND_TIME), intent(in) :: dt
   real(kind=RKIND), intent(in) :: mu,sigmam
   real(kind=RKIND), intent(out) :: gm1,gm2

   Gm1=(1.0_RKIND  - Sigmam * dt / (2.0_RKIND * mu) ) / &
      (1.0_RKIND  + Sigmam * dt / (2.0_RKIND * mu ) )
   Gm2=dt / mu                       / &
      (1.0_RKIND  + Sigmam * dt / (2.0_RKIND * mu))

   if (gm1 < 0.0_RKIND) then !exponential time stepping
      gm1=exp(- Sigmam * dt / (mu ))
      gm2=(1.0_RKIND-gm1)/ Sigmam
   else
      continue
   end if
   return
end subroutine gm1gm2

!!!!!!! dispersive media sgg 12/05/16 
subroutine g1g2_Dispersive(dt,epsilonValue,sigma,G1,G2,Beta,Kappa,G3,numpolres,a11,c11)
   real(kind=RKIND), intent(in) :: epsilonValue,sigma
   real(kind=RKIND_TIME), intent(in) :: dt
   real(kind=RKIND), intent(out) :: g1,g2
   complex(kind=ckind), intent(in), allocatable, dimension(:) :: a11, c11
   integer(kind=4) :: numpolres, i1
   real(kind=RKIND) :: tempo
!!!SGBC dispersive 12/05/16
   complex(kind=CKIND), pointer, dimension(:) :: Beta,Kappa,G3

     do i1=1,numpolres
         Kappa(i1) =(1.0_RKIND + a11(i1)*dt/2.0_RKIND)/&
                    (1.0_RKIND - a11(i1)*dt/2.0_RKIND)
         Beta(i1)=  (C11(i1)*dt) /&
                    (1.0_RKIND - a11(i1)*dt/2.0_RKIND)
     end do
     tempo=0.0_RKIND
     do i1=1,NumPolRes
         tempo=tempo+real(Beta(i1))
     end do
     G1=                        (2.0_RKIND * epsilonValue + tempo - sigma*dt) / & !note Do not be tempted to change this sign. It matches han dutton 130516 and edispersives 
                                (2.0_RKIND * epsilonValue + tempo + sigma*dt)
     G2=         2.0_RKIND *dt/ (2.0_RKIND * epsilonValue + tempo + sigma*dt)
!!!! here exponential time stepping does not fit
     do i1=1,NumPolRes
         G3(i1)=G2/2.0_RKIND * (1.0_RKIND+Kappa(i1))
     end do
   return
end subroutine g1g2_Dispersive

subroutine StoreFieldsSGBCs(stochastic)
      integer(kind=4) :: conta,i,k1
      logical :: SGBCDispersive,stochastic
      type(SGBCSurface_t), pointer :: compo
      do conta=1,malon%numnodes
         compo => malon%Nodes(conta)
         write(14,err=634) (compo%E     (i),i=-compo%depth,compo%depth)

         if (compo%SGBCcrank)  then 
             write(14,err=634) (compo%E_past(i),i=-compo%depth,compo%depth)
         end if
         write(14,err=634) compo%Hyee__left
         write(14,err=634) compo%Hyee_right
         write(14,err=634) (compo%H     (i),i=-compo%depth,compo%depth-1)
         !
         if (malon%SGBCDispersive) then 
             write(14,err=634) (compo%EDis(i)%fieldPrevious, i=-compo%depth,compo%depth)
             do k1=1,compo%NumPolRes
                write(14,err=634) (compo%EDis(i)%current(k1), i=-compo%depth,compo%depth)
             end do
         end if
         
      end do

   goto 635
634   call print11(0,SEPARADOR//separador//separador)
   call print11(0,'SGBC: ERROR WRITING RESTARTING FIELDS. IGNORING AND CONTINUING')
   call print11(0,SEPARADOR//separador//separador)          
635   return
end subroutine StoreFieldsSGBCs


subroutine DestroySGBCs(sgg)

   type(SGGFDTDINFO_t), intent(inout) :: sgg
   integer(kind=4) :: i,conta

   !free up memory
   do i=1,sgg%NumMedia
      if (allocated(malon%dispersiveMedia)) then
          if (allocated(malon%dispersiveMedia(i)%a11)) deallocate(malon%dispersiveMedia(i)%a11)
          if (allocated(malon%dispersiveMedia(i)%c11)) deallocate(malon%dispersiveMedia(i)%c11)
      end if
      if ((sgg%Med(i)%Is%SGBC).and.(.not.sgg%Med(i)%Is%PML))  deallocate(sgg%Med(i)%Multiport)      
   end do
   if (allocated(malon%dispersiveMedia)) deallocate(malon%dispersiveMedia)
   !
   do conta=1,malon%numnodes
      if (allocated(malon%Nodes(conta)%d))  deallocate(malon%Nodes(conta)%d) !CRANK-NICOLSON AUXILIARY
      if (allocated(malon%Nodes(conta)%beta%val))  deallocate(malon%Nodes(conta)%beta%val) !CRANK-NICOLSON dispersive AUXILIARY
      if (allocated(malon%Nodes(conta)%kappa%val))  deallocate(malon%Nodes(conta)%kappa%val) !CRANK-NICOLSON dispersive AUXILIARY
      if (allocated(malon%Nodes(conta)%G3%val))  deallocate(malon%Nodes(conta)%G3%val) !CRANK-NICOLSON dispersive AUXILIARY
      if (allocated(malon%Nodes(conta)%Edis))  deallocate(malon%Nodes(conta)%Edis) !CRANK-NICOLSON dispersive AUXILIARY
     deallocate(malon%Nodes(conta)%GM1_interno ,&           
                malon%Nodes(conta)%GM2_interno ,&
                malon%Nodes(conta)%G1_interno  ,&
                malon%Nodes(conta)%G2_interno  ,&            
                malon%Nodes(conta)%a           ,&
                malon%Nodes(conta)%b           ,&
                malon%Nodes(conta)%c           ,&
                malon%Nodes(conta)%rb          ,&
                malon%Nodes(conta)%rh          ,&
                malon%Nodes(conta)%rhm1       )
      
   end do

   if(allocated(malon%nodes)) deallocate(malon%nodes)
end subroutine


subroutine test_stab(G2,GM2)
   real(kind=RKIND)     , pointer, dimension(:) :: g2, gm2
   integer(kind=4) :: conta
   logical :: unstable
   type(SGBCSurface_t), pointer :: compo
   real(kind=RKIND) :: heur
   character(len=BUFSIZE) :: buff

   heur=1.0_RKIND/sqrt(3.0_RKIND)
   unstable = .false.
   do conta=1,malon%numnodes
      compo => malon%Nodes(conta)
!!!the extremes of the internal E fields
      unstable= unstable.or. &
             (G2(compo%jmed) * Gm2(compo%jmed)  > heur) .or. &
             (compo%G2a(1)   * Gm2(compo%jmed)  > heur) .or. &
             (compo%G2b(1)   * Gm2(compo%jmed)  > heur) .or. &
             (compo%G2a(0)   * Gm2(compo%jmed)  > heur) .or. &
             (compo%G2a(0)   * Gm2(compo%jmed)  > heur)
   end do

   if (unstable) then
        write(buff, *)    'ERROR: SGBCs may become unstable. Reduce cfl'
        call WarnErrReport (buff,.true.)
   end if

   return

end subroutine test_stab

subroutine depth(compo,sgg,jmed,SGBCFreq,SGBCresol,SGBCdepth) 
 type(SGGFDTDINFO_t), intent(in) :: sgg
 real(kind=rkind) :: SGBCFreq,SGBCresol,sigma, epr,epsilonValue,skin_depth,width,widthtotal
 integer(kind=4) :: jmed,i,SGBCdepth,numLayers,precuenta,celdafinal,celdainicial,layerWidth
 integer(kind=4) , pointer, dimension(:) :: layerIndex
 logical :: ultimacapamas1
 character(len=BUFSIZE) :: buff
    
 type(SGBCSurface_t), pointer :: compo

!!!0121 multilayers
 numLayers = sgg%Med(jmed)%multiport(1)%numLayers
 compo%depth=0
 do precuenta=0,1
     if (precuenta==1) then
         if (mod(compo%depth,2)/=0) then
             compo%depth=compo%depth+1 !rounds the total number of layers to an even number
!fills the remainder with the last layer
             ultimacapamas1=.true.
         else 
             ultimacapamas1=.false.
         end if
         compo%depth=int(compo%depth/2.0_RKIND) !divides by 2 because it starts at -compo%depth and reaches +compo%depth 
         if (compo%depth>0) then
             if (.not.allocated(compo%layerIndex))                allocate (compo%layerIndex(-compo%depth:compo%depth-1))
             if (.not.allocated(compo%delta_entreEinterno)) allocate (compo%delta_entreEinterno(-compo%depth:compo%depth-1))
         else
             if (.not.allocated(compo%layerIndex))                allocate (compo%layerIndex(0:0))
             if (.not.allocated(compo%delta_entreEinterno)) allocate (compo%delta_entreEinterno(0:0))
         end if
         
         celdafinal=-compo%depth-1
     end if
     widthtotal=0.; width=0.; sigma=0.; epr=0.; 
     do i=1,numLayers
         width=      sgg%Med(jmed)%multiport(1)%width(i) 
         sigma=      sgg%Med(jmed)%multiport(1)%sigma(i) 
         epr=        sgg%Med(jmed)%multiport(1)%epr(i)  
         epsilonValue=epr * eps0
         widthtotal=widthtotal +     sgg%Med(compo%jmed)%multiport(1)%width(i) 
         skin_depth=1.0_RKIND / (Sqrt(2.0_RKIND)*SGBCFreq*Pi*(Mu0**2*(4*epsilonValue**2.0_RKIND + Sigma**2/(SGBCFreq**2*Pi**2.0_RKIND )))**0.25_RKIND * &
                                 Sin(atan2(2*Pi*epsilonValue*Mu0, -(Mu0*Sigma)/SGBCFreq)/2.0_RKIND))
         if (SGBCdepth==0) then !numlayers must necessarily be 1
             if (numLayers > 1) then
                write(buff, *)   'SGBCDepth=0 and numcapas>1 not compatible. Please, relaunch'
                call StopOnError (0,0,buff)
             else
                 layerWidth=1 !numlayers is necessarily 1 if it continues
             end if
         else if (SGBCdepth>0) then
             layerWidth=SGBCdepth
         else !if it is negative it is calculated with the resol
             layerWidth=1+int(SGBCresol*width/skin_depth)
         end if
         if (layerWidth<2) layerWidth=2 !it is reasonable to never leave it at 1
         !end layers
         if (precuenta==0) then 
             if (SGBCDepth==0) then 
                 compo%depth=0
             else
                 compo%depth=compo%depth+layerWidth
             end if
         else if (precuenta==1) then
             if (SGBCDepth==0) then      !!!bug fixed on 040523
                     celdainicial=0
                     celdafinal=0
                     layerWidth=1
                     compo%layerIndex(celdainicial:celdafinal) = i
                     compo%delta_entreEinterno(celdainicial:celdafinal)=width/layerWidth
                     continue
             else  
             celdainicial=celdafinal+1
             celdafinal=celdainicial+layerWidth-1
             if ((i==numLayers).and.ultimacapamas1) then
!fills the remainder with the last layer if it is not an exact division
                     layerWidth=layerWidth+1
                     celdafinal=celdafinal+1
             end if
             compo%layerIndex(celdainicial:celdafinal) = i
             compo%delta_entreEinterno(celdainicial:celdafinal)=width/layerWidth
             continue
         end if
         end if
     end do
     if (precuenta==1) then
         if ((celdafinal/=compo%depth-1).and.(compo%depth/=0)) then
                write(buff, *)   'Buggy error redondeo final ultima capa. '
                call StopOnError (0,0,buff)
         end if
     end if
 end do
!!!end 02121
 return
end subroutine depth

function GetSGBCs() result(r)
   type(Malon_t), pointer  :: r
   r=>malon
   return
end function


   
!!!!!!tridiagonal solver

subroutine solve_tridiag_distintos(aa,bb,cc,a1,b1,c1,an,bn,cn,d,x,n)
   implicit none
   !  a - sub-diagonal (means it is the diagonal below the main diagonal)
   !  b - the main diagonal
   !  c - sup-diagonal (means it is the diagonal above the main diagonal)
   !  d - right part
   !  x - the answer
   !  n - number of equations

   integer,intent(in) :: n
 real(kind=RKIND) ,intent(in),dimension(n) :: aa,bb,cc
 real(kind=RKIND) ,intent(in) :: a1,b1,c1,an,bn,cn
 real(kind=RKIND) ,dimension(1:n) :: a,b,c
   real(kind=RKIND) ,dimension(n),intent(in) :: d
   real(kind=RKIND) ,dimension(n),intent(out) :: x
   real(kind=RKIND) ,dimension(n) :: cp,dp
   real(kind=RKIND) :: m
   integer i

   a(1)=a1
   b(1)=b1
   c(1)=c1
   a(n)=an
   b(n)=bn
   c(n)=cn
 a(2:n-1)=aa(2:n-1)
 b(2:n-1)=bb(2:n-1)
 c(2:n-1)=cc(2:n-1)
    !  initialize c-prime and d-prime
   cp(1) = c(1)/b(1)
   dp(1) = d(1)/b(1)
   ! solve for vectors c-prime and d-prime
   do i = 2,n
      m = b(i)-cp(i-1)*a(i)
      cp(i) = c(i)/m
      dp(i) = (d(i)-dp(i-1)*a(i))/m
   end do
   ! initialize x
   x(n) = dp(n)
   ! solve for x from the vectors c-prime and d-prime
   do i = n-1, 1, -1
      x(i) = dp(i)-cp(i)*x(i+1)
   end do
   return
end subroutine solve_tridiag_distintos

!!!!!!tridiagonal solver

   subroutine solve_tridiag_iguales(aa,bb,cc,a1,b1,c1,an,bn,cn,d,x,n)
   implicit none
   !  a - sub-diagonal (means it is the diagonal below the main diagonal)
   !  b - the main diagonal
   !  c - sup-diagonal (means it is the diagonal above the main diagonal)
   !  d - right part
   !  x - the answer
   !  n - number of equations

   integer,intent(in) :: n
   real(kind=RKIND) ,intent(in) :: aa,bb,cc,a1,b1,c1,an,bn,cn
   real(kind=RKIND) ,dimension(n) :: a,b,c
   real(kind=RKIND) ,dimension(n),intent(in) :: d
   real(kind=RKIND) ,dimension(n),intent(out) :: x
   real(kind=RKIND) ,dimension(n) :: cp,dp
   real(kind=RKIND) :: m
   integer i

   a(1)=a1
   b(1)=b1
   c(1)=c1
   a(n)=an
   b(n)=bn
   c(n)=cn
   a(2:n-1)=aa
   b(2:n-1)=bb
   c(2:n-1)=cc
    !  initialize c-prime and d-prime
   cp(1) = c(1)/b(1)
   dp(1) = d(1)/b(1)
   ! solve for vectors c-prime and d-prime
   do i = 2,n
      m = b(i)-cp(i-1)*a(i)
      cp(i) = c(i)/m
      dp(i) = (d(i)-dp(i-1)*a(i))/m
   end do
   ! initialize x
   x(n) = dp(n)
   ! solve for x from the vectors c-prime and d-prime
   do i = n-1, 1, -1
      x(i) = dp(i)-cp(i)*x(i+1)
   end do
   return
   end subroutine solve_tridiag_iguales           


end module SGBC_nostoch_m

