
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module SGBCs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!17/08/15 update!!!!!!!!!!!
!!!Elimino el tratamiento de los campos magneticos de SGBC para programar un multiSGBC 
!!!solo teniendo en cuenta los parametros efectivos y sin actualizar los magneticos.
!!!Mantengo en el fichero SGBC_pre170815_noupdateababienH.F90 la version antigua
!!!
!!! 211115 NO ESTA HECHO EL PROMEDIADO DE filo_placaS QUE SE HACE EN COMPOSITES. 
!!!!       LOS SGBCS SE ASIGNAN EN UN FIRST-COME-FIRST-SERVE BASIS. 
!!!!       EL TRATAMIENTO DE filo_placaS DETECTA ARISTAS Y HACE ALGO SIMILAR A LO DEL SHARED
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module SGBC_nostoch

use Report

use fdetypes
implicit none
private
!structures needed by the SGBC
type :: val_t
    Complex (Kind=CKIND), allocatable, dimension ( : )      ::  val
end type


type :: MDfield_t
   REAL (KIND=RKIND), pointer                 ::  FieldPresent !apunta al campo del background
   REAL (KIND=RKIND)                          ::  FieldPrevious
   Complex (Kind=CKIND), pointer, dimension ( : )   ::  Current
end type


TYPE  ::  SGBCSurface_t
   REAL (KIND=RKIND), allocatable, dimension (:) :: E,H,E_past
   REAL (KIND=RKIND), pointer ::  Efield,Ha_Plus,Ha_Minu,Hb_Plus,Hb_Minu
   REAL (KIND=RKIND), allocatable, dimension (:)          ::  delta_entreEinterno
   REAL (KIND=RKIND), dimension (0:1)  ::  g1,g2a,g2b
!!!SGBC dispersivos 12/05/16
   type (MDfield_t), allocatable, dimension (:) :: EDis
   integer (kind=4) :: numpolres
   type (val_t) :: Beta,Kappa,G3
!!!!!!!      
   logical :: correct_ha, correct_hb, es_unfilo_placa
   
   integer (kind=4) :: depth,jmed
   INTEGER (KIND=4), allocatable, dimension (:) ::capa !!!0121
   REAL (KIND=RKIND) , allocatable, dimension (:)   :: G2_interno,GM2_interno,G1_interno,GM1_interno   
   REAL (KIND=RKIND)          :: GM2_externo   !no se precisa gm1_externo porque fuera no hay conductividad magnetica y es trivialmente 1. El gm2_externo tiene sentido almecnarlo porque aun no habiendo conductividad, no es la unidad
   REAL (KIND=RKIND)          :: Hyee__left, Hyee_right      
!!!!! Crank-Nicolson 311015
   REAL (KIND=RKIND) , allocatable, dimension (:)   :: a,b,c,rb,rh,rhm1
   REAL (KIND=RKIND)                                :: a1,b1,c1,rb1,rh1,an,bn,cn,rbn,rhn 
   REAL (KIND=RKIND) , allocatable, dimension (:) :: D !termino independiente CRANK-NICOLSON
   logical :: SGBCCrank

   REAL (KIND=RKIND)  :: transversalDeltaE,transversalDeltaH,alignedlDeltaH
   integer (kind=4), dimension (0:1)  :: med
   COMPLEX (kind=ckind), allocatable, dimension (:) :: a11, c11
END TYPE SGBCSurface_t


type :: MalDisp_t
    integer (kind=4) :: numpolres
    COMPLEX (kind=ckind), allocatable, dimension (:) :: a11, c11
end type

TYPE  ::  Malon_t
    logical :: SGBCdispersive
   integer (kind=4)   ::   NumNodes
   type (SGBCSurface_t), allocatable, dimension (:) :: nodes
   type (MalDisp_t), allocatable, dimension(:) :: mediosDis
end type Malon_t



!!!variables globales del modulo  
type (Malon_t), save, target   ::  malon
!
REAL (KIND=RKIND), save           ::  eps0,mu0,zvac,cluz
logical, save  :: SGBCcrank,SGBCDispersive
REAL (KIND=RKIND), save  :: SGBCFreq,SGBCresol
integer (kind=4), save:: SGBCdepth
!!!
public Malon_t,SGBCSurface_t !el tipo es publico
public AdvanceSGBCE,AdvanceSGBCH,InitSGBCs,DestroySGBCs,StoreFieldsSGBCs,calc_SGBCconstants,GetSGBCs

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine to initialize the parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine InitSGBCs(sgg,media,Ex,Ey,Ez,Hx,Hy,Hz,IDxe,IDye,IDze,IDxh,IDyh,IDzh, &
                     layoutnumber,size,g,ThereAreSGBCs,resume,temp_SGBCcrank,temp_SGBCFreq,temp_SGBCresol,temp_SGBCDepth,temp_SGBCDispersive, &
                     eps00,mu00,simu_devia,stochastic)
   type(media_matrices_t), intent(in) :: media
   logical :: simu_devia,stochastic 
   REAL (KIND=RKIND)           ::  eps00,mu00
   logical :: temp_SGBCcrank,temp_SGBCDispersive
   type(constants_t), intent(inout) :: g
   type (SGGFDTDINFO), intent(INOUT)     ::  sgg !ojo pq se machacan los epr, mur, sigma, sigmam en caso de materiales dispersivos
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

   REAL (KIND=RKIND)  :: temp_SGBCFreq,temp_SGBCresol, rra,rrb,rrc,rrd
   REAL (KIND=RKIND)  :: signo,g1eff_0,g1eff_1,g2eff_0,g2eff_1,Sigmam,Epsilon,Mu,Sigma   
   REAL (KIND=RKIND)  :: factor
   REAL (KIND=RKIND) , allocatable, dimension(:,:) :: derivcte

   integer (kind=4), intent(in) :: layoutnumber,size,temp_SGBCdepth
   logical, INTENT(IN)  :: resume
   logical, INTENT(OUT)  ::  ThereAreSGBCs
   integer (kind=4)  ::  jmed,j1,conta,k1,i1,SGBCdir,i,filo_placas,idummy,numpolres,ient,incert,maxnumcapas,ii
   character(len=BUFSIZE) :: buff
   character (LEN=BUFSIZE)  ::  whoami
   type (SGBCSurface_t), pointer :: compo,compo_temp
   logical :: unstable, errnofile,es_unfilo_placa
   COMPLEX (kind=ckind) :: value1, value2
   character (LEN=BUFSIZE)                            ::   ficheropolos

   eps0=eps00; mu0=mu00; !chapuz para convertir la variables de paso en globales
   SGBCcrank        = temp_SGBCcrank     
   SGBCDispersive   = temp_SGBCDispersive
   SGBCFreq         = temp_SGBCFreq      
   SGBCresol        = temp_SGBCresol     
   SGBCdepth        = temp_SGBCdepth     
!
!!!
   write(whoami,'(a,i5,a,i5,a)') '(',layoutnumber+1,'/',size,') '
   unstable=.false.

!

!   
   malon%SGBCDispersive=SGBCDispersive
   ThereAreSGBCs=.FALSE.
   do jmed=1,sgg%NumMedia
      if (SGG%Med(jmed)%Is%SGBC) then
         ThereAreSGBCs=.true.
      endif
   end do

!pre-cuenta los medios
   conta=0
   Do k1=sgg%SINPMLSweep(iEx)%ZI,sgg%SINPMLSweep(iEx)%ZE
      Do j1=sgg%SINPMLSweep(iEx)%YI,sgg%SINPMLSweep(iEx)%YE
         Do i1=sgg%SINPMLSweep(iEx)%XI,sgg%SINPMLSweep(iEx)%XE
            jmed=media%sggMiEx(i1,j1,k1)
            if (SGG%Med(jmed)%Is%SGBC)  conta=conta+1
         end do
      end do
   end do
   Do k1=sgg%SINPMLSweep(iEy)%ZI,sgg%SINPMLSweep(iEy)%ZE
      Do j1=sgg%SINPMLSweep(iEy)%YI,sgg%SINPMLSweep(iEy)%YE
         Do i1=sgg%SINPMLSweep(iEy)%XI,sgg%SINPMLSweep(iEy)%XE
            jmed=media%sggMiEy(i1,j1,k1)
            if (SGG%Med(jmed)%Is%SGBC)  conta=conta+1
         end do
      end do
   end do

   Do k1=sgg%SINPMLSweep(iEz)%ZI,sgg%SINPMLSweep(iEz)%ZE
      Do j1=sgg%SINPMLSweep(iEz)%YI,sgg%SINPMLSweep(iEz)%YE
         Do i1=sgg%SINPMLSweep(iEz)%XI,sgg%SINPMLSweep(iEz)%XE
            jmed=media%sggMiEz(i1,j1,k1)
            if (SGG%Med(jmed)%Is%SGBC) conta=conta+1
         end do
      end do
   end do
   !!!!!!!!!!!!!!!!!!!!!!
   ThereAreSGBCs=ThereAreSGBCs.and.(conta /=0)
   if (.not.thereareSGBCs) then
      return
   endif
   malon%NumNodes=conta
   allocate (malon%Nodes(1 : malon%NumNodes))
   !!!!DISPERSIVOS
   allocate (malon%mediosDis(1:sgg%NumMedia))
   malon%mediosDis(:)%numpolres=0
   !
   !!!!!!!! dispersivos SGBC sgg 12/05/15   
!070717
!!!primero chequeo si existe el fichero de polos. si existe es que SGBCDispersive paro e informo
!!!!ojo lo del sgbc de momento es global (070717). Si el switch esta todos son sgbc dispersive. pero esta semipreparado para que sea medio a medio. cambiar algun dia....
!!!!0121 quito esto porque mira el fichero polos que genera ugrmat_multilayer y activaba dispersive por defecto sin mas
   !! do jmed=1,sgg%NumMedia
   !!      if ((.not.SGG%Med(jmed)%Is%SGBCDispersive).and.(SGG%Med(jmed)%Is%SGBC).and.(.not.(SGG%Med(jmed)%Is%PML))) then    
   !!          ficheropolos=SGG%Med(jmed)%multiport(1)%multiportFileZ11 !aunque le llamo Z tiene la sintaxis de un Edispersive ISOTROPO CON EL NUEVO STANDARD (VER LINEAS 6749 DE NFDEPARSER). 
   !!          ! SOLO LEO LOS PRIMEROS POLOS. eL RESTO DE DATOS LOS TIRO (INFORMACION DE POLOS DE SEGUNDO ORDEN, POLOS MAGNETICOS, ANISOTROPIAS...)
   !!          !nuevo estio del ficheros sin el _z11
   !!          i1=index(ficheropolos,'_z11.txt')
   !!          ficheropolos=trim(adjustl(ficheropolos(1:i1-1)))
   !!  !
   !!          errnofile=.false.
   !!          INQUIRE (FILE=trim(adjustl(ficheropolos)), EXIST=errnofile)
   !!          if (errnofile) then
   !!               WRITE (buff, *)    'ERROR: -sgbcdispersive not used and poles files exist. Correct .nfde or issue -sgbcdispersive'
   !!               CALL WarnErrReport (buff,.false.)
   !!               stop
   !!               SGG%Med(jmed)%Is%SGBCDispersive=.true.
   !!               SGBCdispersive=.true. !en caso de ignoreerrors puede seguir simplemente con sgbcdispersive a true
   !!          endif
   !!      endif
   !!end do 
   if (SGBCDispersive) then
       do jmed=1,sgg%NumMedia
          if ((SGG%Med(jmed)%Is%SGBCDispersive).and.(.not.(SGG%Med(jmed)%Is%PML))) then    
!!!solo una capa de dispersivo
              if (sgg%Med(jmed)%multiport(1)%numcapas>1) then
                 buff='No more than 1 layer of dispersive SGBC currently supported'
                 call StopOnError(layoutnumber,size,buff)
              endif
!!!!!!!!
              ficheropolos=SGG%Med(jmed)%multiport(1)%multiportFileZ11 !aunque le llamo Z tiene la sintaxis de un Edispersive ISOTROPO CON EL NUEVO STANDARD (VER LINEAS 6749 DE NFDEPARSER). 
              ! SOLO LEO LOS PRIMEROS POLOS. eL RESTO DE DATOS LOS TIRO (INFORMACION DE POLOS DE SEGUNDO ORDEN, POLOS MAGNETICOS, ANISOTROPIAS...)
              !nuevo estio del ficheros sin el _z11
              i1=index(ficheropolos,'_z11.txt')
              ficheropolos=trim(adjustl(ficheropolos(1:i1-1)))
      !
              errnofile=.false.
              INQUIRE (FILE=trim(adjustl(ficheropolos)), EXIST=errnofile)
              if (.not.errnofile) then
                 buff='FILE '//trim(adjustl(ficheropolos))//' DOES NOT EXIST'
                 call StopOnError(layoutnumber,size,buff)
              endif
              open (7345,file=trim(adjustl(ficheropolos)),form='formatted')
              READ (7345,*) rra,rrb,rrc,rrd 
              rrb= rrb/eps0 ;  rrc = rrc/mu0 !no les afecta el permit scaling creo 071118 pq son relativos a la entrada del programa que DEBE ENTRAR CON los eps0 y mu0 autenticos
              SGG%Med(jmed)%multiport(1)%sigma(1)=rra; 
              SGG%Med(jmed)%multiport(1)%epr(1)=rrb; 
              SGG%Med(jmed)%multiport(1)%mur(1)=rrc; 
              SGG%Med(jmed)%multiport(1)%sigmam(1)=rrd;
              SGG%Med(jmed)%sigma=rra; 
              SGG%Med(jmed)%epr=rrb; 
              SGG%Med(jmed)%mur=rrc; 
              SGG%Med(jmed)%sigmam=rrd;
              !
              READ (7345,*) numpolres, IDUMMY, IDUMMY, IDUMMY
              malon%mediosDis(jmed)%numpolres = numpolres
              allocate (malon%mediosDis(jmed)%a11(1:numpolres)) 
              allocate (malon%mediosDis(jmed)%c11(1:numpolres)) 
              DO i = 1, numpolres
                read(7345,*) value1, value2
                malon%mediosDis(jmed)%c11 (i) = (value1) 
                malon%mediosDis(jmed)%a11 (i) = - (value2) !el polo de EM esta cambiado de signo !ver tambien preprocess
              END DO          
              close (7345)
!!!movido 071118 al calculo de constantes para permit scaling
!!!!!!!!070717  recalculo y machaco los G1,G2,GM1, y GM2 con los que aparecen en el fichero de dispersivos aunque creo que no se usa para nada
!!!!!                  Sigmam  =      SGG%Med(jmed)%multiport(1)%sigmam(1)
!!!!!                  Epsilon = Eps0*SGG%Med(jmed)%multiport(1)%epr(1)
!!!!!                  Mu      = Mu0* SGG%Med(jmed)%multiport(1)%mur(1)
!!!!!                  Sigma   =      SGG%Med(jmed)%multiport(1)%sigma(1)
!!!!!                  G1(jmed)=(1 -  Sigma * sgg%dt / (2.0_RKIND * Epsilon ) ) / (1.0_RKIND + Sigma * sgg%dt / (2.0_RKIND * Epsilon ))
!!!!!                  G2(jmed)=sgg%dt /Epsilon                        / (1.0_RKIND + Sigma * sgg%dt / (2.0_RKIND * Epsilon ))
!!!!!                  if (g1(jmed) < 0.0_RKIND) then !exponential time stepping
!!!!!                     g1(jmed)=exp(- Sigma * sgg%dt / (Epsilon ))
!!!!!                     g2(jmed)=(1.0_RKIND-g1(jmed))/ Sigma
!!!!!                  endif
!!!!!                  GM1(jmed)=(1- SigmaM*sgg%dt/(2.0_RKIND *  Mu )) /(1.0_RKIND + SigmaM*sgg%dt/(2.0_RKIND *  Mu ))
!!!!!                  GM2(jmed)=sgg%dt/ Mu                   /(1.0_RKIND + SigmaM*sgg%dt/(2.0_RKIND *  Mu ))
!!!!!                  if (gm1(jmed) < 0.0_RKIND) then !exponential time stepping
!!!!!                     gm1(jmed)=exp(- Sigmam * sgg%dt / (Mu ))
!!!!!                     gm2(jmed)=(1.0_RKIND-gm1(jmed))/ Sigmam
!!!!!                  endif
!!!!!!!!!fin 070717
         endif
      end do
   endif
   !!!!!!!! del SGBCdispersive
   conta=0
!asigna los signos del rotacional de H y los transversaldelta.0->filo_placaizquierdo, 1->filo_placaderecho y las variable es_unfilo_placa para los filos de la lamina SGBC      
Do k1=sgg%SINPMLSweep(iEx)%ZI,sgg%SINPMLSweep(iEx)%ZE
    Do j1=sgg%SINPMLSweep(iEx)%YI,sgg%SINPMLSweep(iEx)%YE
         Do i1=sgg%SINPMLSweep(iEx)%XI,sgg%SINPMLSweep(iEx)%XE
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
               endif
!!!!!
               conta=conta+1
               compo => malon%Nodes(conta)
               compo%es_unfilo_placa = es_unfilo_placa
               SGBCdir=abs(SGG%Med(jmed)%Multiport(1)%Multiportdir)
               select case (SGBCdir)
                case (iEy)
                  compo%transversalDeltaE =    1.0_RKIND/IDYe(   j1     )
                  compo%transversalDeltaH =    1.0_RKIND/IDYh(   j1     )
                  compo%alignedlDeltaH     =    1.0_RKIND/Idzh(        k1)
                  compo%med(1) =                     media%sggMiHz(i1,j1  ,k1)
                  compo%med(0) =                     media%sggMiHz(i1,j1-1,k1)
                  compo%Correct_Ha=.true. !son ciclicos a,b -> x,y,z
                  compo%Correct_Hb=.false.
                case (iEz)
                  compo%transversalDeltaE    = 1.0_RKIND/IDze(      k1  )
                  compo%transversalDeltaH    = 1.0_RKIND/IDzh(      k1  )
                  compo%alignedlDeltaH        = 1.0_RKIND/IDYh(   j1     )
                  compo%med(1) =                     media%sggMiHy(i1,j1,k1  )
                  compo%med(0) =                     media%sggMiHy(i1,j1,k1-1)
                  compo%Correct_Hb=.true.
                  compo%Correct_Ha=.false.
                case DEFAULT
                  WRITE (buff, *)    'Buggy ERROR: In SGBCs. '
                  CALL WarnErrReport (buff,.TRUE.)
               end select
               if ((compo%med(1) == 0).or.(compo%med(0) == 0).or.(sgg%med(compo%med(1))%Is%PEC).or.(sgg%med(compo%med(0))%Is%PEC)) then
                  WRITE (buff, *)    'Buggy ERROR: In SGBCs. '
                  CALL WarnErrReport (buff,.TRUE.)
               endif
!!!!!!!!!
               compo%jmed       = jmed
               compo%Efield    => Ex(i1,j1  ,k1  )
               compo%Ha_Plus   => Hz(i1,j1  ,k1)
               compo%Ha_Minu   => Hz(i1,j1-1,k1)
               compo%Hb_Plus   => Hy(i1,j1  ,k1)
               compo%Hb_Minu   => Hy(i1,j1  ,k1-1)
               call  depth(compo,sgg,jmed,SGBCFreq,SGBCresol,SGBCdepth)
               if (compo%depth==0) then
                    compo%SGBCCrank=.false.
               elseif (compo%depth<0) then
                  WRITE (buff, *)    'Buggy ERROR: In SGBCs compo%depth<0. ',compo%depth
                  CALL WarnErrReport (buff,.TRUE.)
               else
                    compo%SGBCCrank=SGBCCrank
               endif
               allocate (compo%E         (-compo%depth:compo%depth  ))
               if (compo%depth>0) allocate (compo%H         (-compo%depth:compo%depth-1)) 
               allocate(compo%E_past(-compo%depth:compo%depth  )) !no se precisa en yee pero se comunica en mpi_stochastic. movido fuera del if sigugiente 170519
               if (compo%SGBCCrank)  then
                    allocate(compo%d     (-compo%depth:compo%depth  )) 
               endif

               compo%numpolres=malon%MediosDis(compo%jmed)%numpolres !duplico esta info
               if (SGBCDispersive) then 
                 allocate (compo%a11(1:compo%numpolres)) 
                 allocate (compo%c11(1:compo%numpolres)) 
                 compo%a11 = malon%MediosDis(compo%jmed)%a11
                 compo%c11 = malon%MediosDis(compo%jmed)%c11
                 allocate (compo%beta%val(1:compo%numpolres))
                 allocate (compo%kappa%val(1:compo%numpolres))
                 allocate (compo%G3%val(1:compo%numpolres))
             !!  call calc_g1g2(sgg,GM2,compo)   ! permit scal 071118 
                 allocate (compo%EDis   (-compo%depth:compo%depth  ))
                 do ient=-compo%depth , compo%depth
                     allocate (compo%EDis(ient)%Current(1 : malon%MediosDis(jmed)%numpolres))
                     compo%EDis(ient)%FieldPresent => compo%E(ient)
                 end do
               endif
             endif
         end do
      end do
   end do
   !!!!!!!!!!
   Do k1=sgg%SINPMLSweep(iEy)%ZI,sgg%SINPMLSweep(iEy)%ZE
      Do j1=sgg%SINPMLSweep(iEy)%YI,sgg%SINPMLSweep(iEy)%YE
         Do i1=sgg%SINPMLSweep(iEy)%XI,sgg%SINPMLSweep(iEy)%XE
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
               endif
!!!!!
               conta=conta+1
               compo => malon%Nodes(conta)
               compo%es_unfilo_placa = es_unfilo_placa
               SGBCdir=abs(SGG%Med(jmed)%Multiport(1)%Multiportdir)
               select case (SGBCdir)
                case (iEz)
                  compo%transversalDeltaE = 1.0_RKIND/IDze(      k1  )
                  compo%transversalDeltaH = 1.0_RKIND/IDzh(      k1  )
                  compo%alignedlDeltaH     = 1.0_RKIND/Idxh(i1        )
                  compo%med(1) =                  media%sggMiHx(i1,j1,k1  )
                  compo%med(0) =                  media%sggMiHx(i1,j1,k1-1)
                  compo%Correct_Ha=.true.
                  compo%Correct_Hb=.false.
                case (iEx)
                  compo%transversalDeltaE = 1.0_RKIND/IDxe(i1        )
                  compo%transversalDeltaH = 1.0_RKIND/IDxh(i1        )
                  compo%alignedlDeltaH     = 1.0_RKIND/IDzh(        k1)
                  compo%med(1) =                  media%sggMiHz(i1  ,j1,k1)
                  compo%med(0) =                  media%sggMiHz(i1-1,j1,k1)
                  compo%Correct_Ha=.false.
                  compo%Correct_Hb=.true.
                case DEFAULT
                  WRITE (buff, *)    'Buggy ERROR: In SGBCs. '
                  CALL WarnErrReport (buff,.TRUE.)
               end select
               if ((compo%med(1) == 0).or.(compo%med(0) == 0).or.(sgg%med(compo%med(1))%Is%PEC).or.(sgg%med(compo%med(0))%Is%PEC)) then
                  WRITE (buff, *)    'Buggy ERROR: In SGBCs. '
                  CALL WarnErrReport (buff,.TRUE.)
               endif
!!!!   
               compo%jmed       = jmed
               compo%Efield    => Ey(i1  ,j1  ,k1  )
               compo%Ha_Plus   => Hx(i1  ,j1  ,k1  )
               compo%Ha_Minu   => Hx(i1  ,j1  ,k1-1)
               compo%Hb_Plus   => Hz(i1  ,j1  ,k1  )
               compo%Hb_Minu   => Hz(i1-1,j1  ,k1  )
               call  depth(compo,sgg,jmed,SGBCFreq,SGBCresol,SGBCdepth)
               if (compo%depth==0) then
                    compo%SGBCCrank=.false.
               elseif (compo%depth<0) then
                  WRITE (buff, *)    'Buggy ERROR: In SGBCs compo%depth<0. ',compo%depth
                  CALL WarnErrReport (buff,.TRUE.)
               else
                    compo%SGBCCrank=SGBCcrank
               endif
               allocate (compo%E         (-compo%depth:compo%depth  ))
               if (compo%depth>0) allocate (compo%H         (-compo%depth:compo%depth-1)) 
               allocate(compo%E_past(-compo%depth:compo%depth  ))
               if (compo%SGBCcrank)  then
                    allocate(compo%d     (-compo%depth:compo%depth  )) 
               endif

               compo%numpolres=malon%MediosDis(compo%jmed)%numpolres !duplico esta info
               if (SGBCDispersive) then 
                 allocate (compo%a11(1:compo%numpolres)) 
                 allocate (compo%c11(1:compo%numpolres)) 
                 compo%a11 = malon%MediosDis(compo%jmed)%a11
                 compo%c11 = malon%MediosDis(compo%jmed)%c11
                 allocate (compo%beta%val(1:compo%numpolres))
                 allocate (compo%kappa%val(1:compo%numpolres))
                 allocate (compo%G3%val(1:compo%numpolres))
                 !! call calc_g1g2(sgg,GM2,compo)     ! permit scal 071118
                 allocate (compo%EDis   (-compo%depth:compo%depth  ))
                 do ient=-compo%depth , compo%depth
                     allocate (compo%EDis(ient)%Current(1 : malon%MediosDis(jmed)%numpolres))
                     compo%EDis(ient)%FieldPresent => compo%E(ient)
                 end do
               endif
            endif
         end do
      end do
   end do

   Do k1=sgg%SINPMLSweep(iEz)%ZI,sgg%SINPMLSweep(iEz)%ZE
      Do j1=sgg%SINPMLSweep(iEz)%YI,sgg%SINPMLSweep(iEz)%YE
         Do i1=sgg%SINPMLSweep(iEz)%XI,sgg%SINPMLSweep(iEz)%XE
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
               endif
!!!!!
               conta=conta+1
               compo => malon%Nodes(conta)
               compo%es_unfilo_placa = es_unfilo_placa
               SGBCdir=abs(SGG%Med(jmed)%Multiport(1)%Multiportdir)
               select case (SGBCdir)
                case (iEx)
                  compo%transversalDeltaE = 1.0_RKIND/IDxE(i1        )
                  compo%transversalDeltaH = 1.0_RKIND/IDxh(i1        )
                  compo%alignedlDeltaH     = 1.0_RKIND/IDyh(     j1   )
                  compo%med(1) =                  media%sggMiHy(i1  ,j1,k1)
                  compo%med(0) =                  media%sggMiHy(i1-1,j1,k1)
                  compo%Correct_Ha=.true.
                  compo%Correct_Hb=.false.
                case (iEy)
                  compo%transversalDeltaE = 1.0_RKIND/IDyE(    j1     )
                  compo%transversalDeltaH = 1.0_RKIND/IDyh(    j1     )
                  compo%alignedlDeltaH     = 1.0_RKIND/IDxh(i1        )
                  compo%med(1) =                    media%sggMiHx(i1,j1  ,k1)
                  compo%med(0) =                    media%sggMiHx(i1,j1-1,k1)
                  compo%Correct_Ha=.false.
                  compo%Correct_Hb=.true.
                case DEFAULT
                  WRITE (buff, *)    'Buggy ERROR: In SGBCs. '
                  CALL WarnErrReport (buff,.TRUE.)
               end select
               if ((compo%med(1) == 0).or.(compo%med(0) == 0).or.(sgg%med(compo%med(1))%Is%PEC).or.(sgg%med(compo%med(0))%Is%PEC)) then
                  WRITE (buff, *)    'Buggy ERROR: In SGBCs. '
                  CALL WarnErrReport (buff,.TRUE.)
               endif
!!!!!!!!!
               compo%jmed  = jmed
               compo%Efield  => Ez(i1  ,j1  ,k1  )
               compo%Ha_Plus => Hy(i1  ,j1  ,k1  )
               compo%Ha_Minu => Hy(i1-1,j1  ,k1)
               compo%Hb_Plus => Hx(i1  ,j1  ,k1  )
               compo%Hb_Minu => Hx(i1  ,j1-1,k1  )
               call  depth(compo,sgg,jmed,SGBCFreq,SGBCresol,SGBCdepth)
               if (compo%depth==0) then
                    compo%SGBCCrank=.false.
               elseif (compo%depth<0) then
                  WRITE (buff, *)    'Buggy ERROR: In SGBCs compo%depth<0. ',compo%depth
                  CALL WarnErrReport (buff,.TRUE.)
               else
                    compo%SGBCCrank=SGBCcrank
               endif
               allocate (compo%E         (-compo%depth:compo%depth  ))
               if (compo%depth>0) allocate (compo%H         (-compo%depth:compo%depth-1)) 
               allocate(compo%E_past(-compo%depth:compo%depth  ))
               if (compo%SGBCcrank)  then
                    allocate(compo%d     (-compo%depth:compo%depth  )) 
               endif

               compo%numpolres=malon%MediosDis(compo%jmed)%numpolres
               if (SGBCDispersive) then  !duplico esta info
                 allocate (compo%a11(1:compo%numpolres)) 
                 allocate (compo%c11(1:compo%numpolres)) 
                 compo%a11 = malon%MediosDis(compo%jmed)%a11
                 compo%c11 = malon%MediosDis(compo%jmed)%c11
                 allocate (compo%beta%val(1:compo%numpolres))
                 allocate (compo%kappa%val(1:compo%numpolres))
                 allocate (compo%G3%val(1:compo%numpolres))
                 !! call calc_g1g2(sgg,GM2,compo)     ! permit scal 071118
                 allocate (compo%EDis   (-compo%depth:compo%depth  ))
                 do ient=-compo%depth , compo%depth
                     allocate (compo%EDis(ient)%Current(1 : malon%MediosDis(jmed)%numpolres))
                     compo%EDis(ient)%FieldPresent => compo%E(ient)
                 end do
               endif
            endif
         end do
      end do
   end do

   call calc_SGBCconstants(sgg,g,eps0,mu0,stochastic)

!!!reporting de depth
    i=-100
    do conta=1,malon%numnodes
      compo => malon%Nodes(conta)
      if (compo%depth>i) then
         i=compo%depth
         jmed=compo%jmed
      endif
    end do

    WRITE (buff, *)  ' Maximum SGBC depth= ',2*i,' for medium jmed= ',jmed
    CALL WarnErrReport (buff)



!!! no lo uso. no es un criterio riguroso de estabilidad 311015
!!!      call test_stab(G2,GM2)

   !!!!!!!!!resuming
   if (.not.resume) then  
      do conta=1,malon%numnodes
         compo => malon%Nodes(conta)
         compo%E     =0.0_RKIND
         if (compo%SGBCcrank)  then 
             compo%E_past=0.0_RKIND
         endif
         compo%Hyee__left=0.0_RKIND
         compo%Hyee_right=0.0_RKIND
         compo%H     =0.0_RKIND
         !
         if (malon%SGBCDispersive) then 
             Do i=-compo%depth,compo%depth
                compo%EDis(i)%fieldPresent=0.0_rkind
                compo%EDis(i)%fieldPrevious=0.0_rkind
                compo%EDis(i)%current=0.0_rkind
             enddo
         endif
      end do
   else  
      do conta=1,malon%numnodes
         compo => malon%Nodes(conta)
         READ (14) (compo%E     (i),i=-compo%depth,compo%depth)
         READ (14) (compo%E_past(i),i=-compo%depth,compo%depth)
         READ (14) compo%Hyee__left
         READ (14) compo%Hyee_right
         READ (14) (compo%H     (i),i=-compo%depth,compo%depth-1)
         !
         if (malon%SGBCDispersive) then 
             read(14) (compo%EDis(i)%fieldPrevious, i=-compo%depth,compo%depth)
             Do k1=1,compo%NumPolRes
                read(14) (compo%EDis(i)%current(k1), i=-compo%depth,compo%depth)
             enddo
         endif
      end do
   endif
   return

end subroutine InitSGBCs

subroutine calc_SGBCconstants(sgg,g,eps00,mu00,stochastic)
   REAL (KIND=RKIND), intent(IN)           ::  eps00,mu00
   type (SGGFDTDINFO), intent(IN)     ::  sgg 
   REAL (KIND=RKIND), pointer, dimension ( : )   ::  gm1,g1,gm2,g2
   type(constants_t) :: g
 integer :: jmed,conta,i
   REAL (KIND=RKIND) :: sigmam,sigma,mu,epsilon,signo,g1eff_0,g2eff_0,g1eff_1,g2eff_1
   type (SGBCSurface_t), pointer :: compo
   character(len=BUFSIZE) :: buFF
   logical :: stochastic
!
   g1 => g%g1; g2 => g%g2; gm1 => g%gm1; gm2 => g%gm2; 
   eps0=eps00; mu0=mu00; !chapuz para convertir la variables de paso en globales
   zvac=sqrt(mu0/eps0)
   cluz=1.0_RKIND/sqrt(mu0*eps0)
!!!allocateo todas las matrices de constantes
 do conta=1,malon%numnodes
     compo => malon%Nodes(conta)  
     !!!rellamo a depth para que recalcule bien el deltaentreEinterno !110523 necesario para stochastic
     jmed=compo%jmed
     call depth(compo,sgg,jmed,SGBCFreq,SGBCresol,SGBCdepth)
               
     if (.not.allocated(compo%GM1_interno)) then !!!sobran extremos de g y gm pero lo dejo asi para no toquetear mas 0121
         allocate(compo%GM1_interno (-compo%depth  :compo%depth-1) ,&   !se ajustan bien los extremos        
                  compo%GM2_interno (-compo%depth  :compo%depth-1) ,&
                  compo%G1_interno  (-compo%depth+1:compo%depth-1) ,&
                  compo%G2_interno  (-compo%depth+1:compo%depth-1) ,&            
                  compo%a           (-compo%depth  :compo%depth) ,& !se echan 1 mas al ppio y al final pero no lo toco 0121
                  compo%b           (-compo%depth  :compo%depth) ,&
                  compo%c           (-compo%depth  :compo%depth) ,&
                  compo%rb          (-compo%depth  :compo%depth) ,&
                  compo%rh          (-compo%depth  :compo%depth) ,&
                  compo%rhm1        (-compo%depth  :compo%depth))
     endif

 end do
 
 !default absurdos
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
 
!!!!CALCULO DE LOS COEFICIENTES YEE-FDTD (y semilla de los CN-FDTD)

   if (SGBCDispersive) then
       do jmed=1,sgg%NumMedia
          if ((SGG%Med(jmed)%Is%SGBCDispersive).and.(.not.(SGG%Med(jmed)%Is%PML))) then   
!!!071118 para permit scaling
!!!070717  recalculo y machaco los G1,G2,GM1, y GM2 con los que aparecen en el fichero de dispersivos aunque creo que no se usa para nada
               Sigmam  =      SGG%Med(jmed)%multiport(1)%sigmam(1)
               Epsilon = Eps0*SGG%Med(jmed)%multiport(1)%epr(1)
               Mu      = Mu0* SGG%Med(jmed)%multiport(1)%mur(1)
               Sigma   =      SGG%Med(jmed)%multiport(1)%sigma(1)
               G1(jmed)=(1 -  Sigma * sgg%dt / (2.0_RKIND * Epsilon ) ) / (1.0_RKIND + Sigma * sgg%dt / (2.0_RKIND * Epsilon ))
               G2(jmed)=sgg%dt /Epsilon                        / (1.0_RKIND + Sigma * sgg%dt / (2.0_RKIND * Epsilon ))
               if (g1(jmed) < 0.0_RKIND) then !exponential time stepping
                  g1(jmed)=exp(- Sigma * sgg%dt / (Epsilon ))
                  g2(jmed)=(1.0_RKIND-g1(jmed))/ Sigma
               endif
               GM1(jmed)=(1- SigmaM*sgg%dt/(2.0_RKIND *  Mu )) /(1.0_RKIND + SigmaM*sgg%dt/(2.0_RKIND *  Mu ))
               GM2(jmed)=sgg%dt/ Mu                   /(1.0_RKIND + SigmaM*sgg%dt/(2.0_RKIND *  Mu ))
               if (gm1(jmed) < 0.0_RKIND) then !exponential time stepping
                  gm1(jmed)=exp(- Sigmam * sgg%dt / (Mu ))
                  gm2(jmed)=(1.0_RKIND-gm1(jmed))/ Sigmam
               endif
!!!!fin 070717
          endif
      end do
   endif    

#ifdef CompileWithOpenMP
!$OMP  PARALLEL DO  DEFAULT(none) private (compo,buff) shared(malon,sgg,eps00,mu00,GM2,SGBCDispersive)
#endif
 do conta=1,malon%numnodes
     compo => malon%Nodes(conta)
     call calc_g1g2gm1gm2_compo(sgg,compo,eps00,mu00,SGBCDispersive)
!cte de actualizacion de los H externos a la multicapa!!! 0121
     compo%Gm2_externo=Gm2(compo%jmed) / compo%transversalDeltaE !ojo habia compo%transversalDeltaE antes sgg 130516 ! pero yo creo que es  compo%transversalDeltaH! 0121 No. es deltaE pq se usa para actualizar H externo
 end do

#ifdef CompileWithOpenMP
!$OMP  END PARALLEL DO
#endif
!
!

#ifdef CompileWithOpenMP
!$OMP  PARALLEL DO  DEFAULT(none) private (compo,signo,g1eff_0,g1eff_1,g2eff_0,g2eff_1) shared(malon,gm1,gm2)
#endif
    do conta=1,malon%numnodes
      compo => malon%Nodes(conta) 
      if (compo%depth>0) then !uno de los del rotacional
          if (compo%Correct_Ha) then
             signo=+1.0_RKIND
             g1eff_0=   compo%G1(0)  
             g1eff_1=   compo%G1(1)   
             g2eff_0=   signo *compo%G2a(0)  
             g2eff_1=   signo *compo%G2a(1)
          elseif (compo%Correct_Hb) then
             signo=-1.0_RKIND
             g1eff_0=   compo%G1(0)  
             g1eff_1=   compo%G1(1)   
             g2eff_0=   signo *compo%G2b(0)  
             g2eff_1=   signo *compo%G2b(1)
          endif
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
          do i=-compo%depth+1 , compo%depth-1 !los primeros y ultimos no se usan, pero los dejo para no alterar el algoritmo pre 0121
!                 compo%a  (i)        =                      - compo%G2_interno (i  ) * compo%GM2_interno (i  ) /4.0_RKIND
               compo%a  (i)        =                      - compo%G2_interno (i  ) * compo%GM2_interno (i-1) /4.0_RKIND   !!!!el menos 1
               
!                 compo%b  (i)        = 1.0_RKIND            + compo%G2_interno (i  ) * compo%GM2_interno (i  ) /2.0_RKIND   
               compo%b  (i)        = 1.0_RKIND            + compo%G2_interno (i  ) * compo%GM2_interno (i-1) /4.0_RKIND  + compo%G2_interno (i) * compo%GM2_interno (i) /4.0_RKIND   
               
!                 compo%c  (i)        =                      - compo%G2_interno (i  ) * compo%GM2_interno (i  ) /4.0_RKIND     
               compo%c  (i)        =                      - compo%G2_interno (i  ) * compo%GM2_interno (i  ) /4.0_RKIND   !igual porque el +1 esta indexado en i
               
!                 compo%rb (i)        = compo%G1_interno (i) - compo%G2_interno (i  ) * compo%GM2_interno (i  ) /2.0_RKIND
               compo%rb (i)        = compo%G1_interno (i) - compo%G2_interno (i  ) * compo%GM2_interno (i-1) /4.0_RKIND  - compo%G2_interno (i) * compo%GM2_interno (i) /4.0_RKIND   !!el menos 1
               
!                 compo%rh  (i)        =(compo%G2_interno (i  ) * compo%GM1_interno(i  ) + compo%G2_interno  (i  ))/2.0_RKIND
               compo%rh  (i)        =(compo%G2_interno (i  ) * compo%GM1_interno(i  ) + compo%G2_interno  (i  ))/2.0_RKIND  !hay que desglosarlo
               compo%rhm1(i)        =(compo%G2_interno (i  ) * compo%GM1_interno(i-1) + compo%G2_interno  (i  ))/2.0_RKIND
!!!        
          end do
          i=-compo%depth
          compo%a1          =  0.0
          compo%c1          =                 - g2eff_0 * compo%GM2_interno (i) /4.0_RKIND !el que esta a su derecha interno junto con el yee
          compo%b1           = 1.0_RKIND      + g2eff_0 * compo%GM2_interno (i) /4.0_RKIND
          compo%rb1          = g1eff_0 - g2eff_0        * compo%GM2_interno (i) /4.0_RKIND
          compo%rh1          =  (g2eff_0                * compo%GM1_interno (i)+ g2eff_0)/2.0_RKIND
          i=compo%depth
          compo%cn          =  0.0
          compo%an          =                 - g2eff_1 * compo%GM2_interno (i-1) /4.0_RKIND !el que esta a su izquierda interno junto con el yee
          compo%bn           = 1.0_RKIND      + g2eff_1 * compo%GM2_interno (i-1) /4.0_RKIND
          compo%rbn          = g1eff_1 - g2eff_1        * compo%GM2_interno (i-1) /4.0_RKIND
          compo%rhn          =  (g2eff_1                * compo%GM1_interno (i-1) + g2eff_1)/2.0_RKIND
!!!!fin nueva formulacion
!!pre-jav 310116
!!!                 compo%a1=0.0; compo%c1=0.0; compo%b1=1.0; compo%an=0.0; compo%cn=0.0; compo%bn=1.0; compo%rb1=0.0; compo%rbn=0.0; compo%rh1=0.0; compo%rhn=0.0    
!!!!post-jav nueva formulacion jav 310116       
      endif
    end do 
#ifdef CompileWithOpenMP
!$OMP  END PARALLEL DO
#endif
    return
end subroutine calc_SGBCconstants


                        
 subroutine YeeAdvanceSGBCDispersive (tempnode,numpolres,G3,kappa,beta,dt)
 
   REAL (KIND=RKIND), intent(IN)     :: dt
   type (MDfield_t), pointer  ::  tempnode
   integer (kind=4) :: k1,numpolres
   type (val_t) :: Beta,Kappa,G3
   
         Do k1=1,NumPolRes
            tempnode%fieldPresent=tempnode%FieldPresent-REAL (G3%val(k1)*tempnode%current(k1) )
         enddo
         Do k1=1,NumPolRes
            tempnode%current(k1)= Kappa%val(k1)  *tempnode%current(k1) + &
                                  Beta%val(k1)*(tempnode%fieldPresent-tempnode%fieldPrevious) /dt
         enddo
         tempnode%fieldPrevious=tempnode%fieldPresent
         !stores previous field (cuidado no es un apuntamiento sino una igualdad de valores)
         !antes de que re-empieze a calcularlo el algoritmo del background

 end subroutine YeeAdvanceSGBCDispersive
 
                         
 subroutine primero_CNAdvanceSGBCDispersive (tempnode,tempD,numpolres,G3,kappa,beta,dt)
 
   REAL (KIND=RKIND), intent(IN)     :: dt
   REAL (KIND=RKIND),  pointer  ::  tempD
   type (MDfield_t), pointer  ::  tempnode
   integer (kind=4) :: k1,numpolres
   type (val_t) :: Beta,Kappa,G3
         Do k1=1,NumPolRes
            tempD=tempD-REAL (G3%val(k1)*tempnode%current(k1) )
         enddo
 end subroutine primero_CNAdvanceSGBCDispersive
 
 subroutine segundo_CNAdvanceSGBCDispersive (campocalculado,tempnode,numpolres,G3,kappa,beta,dt)
 
   REAL (KIND=RKIND) :: campocalculado 
   REAL (KIND=RKIND), intent(IN)     :: dt
   type (MDfield_t), pointer  ::  tempnode
   integer (kind=4) :: k1,numpolres
   type (val_t) :: Beta,Kappa,G3

         Do k1=1,NumPolRes
            tempnode%fieldPresent=campocalculado
            tempnode%current(k1)= Kappa%val(k1)  *tempnode%current(k1) + &
                                  Beta%val(k1)*(tempnode%fieldPresent-tempnode%fieldPrevious) /dt
         enddo
         tempnode%fieldPrevious=tempnode%fieldPresent
         !stores previous field (cuidado no es un apuntamiento sino una igualdad de valores)
         !antes de que re-empieze a calcularlo el algoritmo del background

 end subroutine segundo_CNAdvanceSGBCDispersive
                        

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine to advance the E field in the SGBC: Usual Yee
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine AdvanceSGBCE(dt,SGBCDispersive,simu_devia,stochastic)
   logical :: simu_devia,stochastic 
   REAL (KIND=RKIND), intent(IN)     :: dt
   logical :: SGBCDispersive
   integer (kind=4)  ::  conta
!       call AdvanceSGBCE_single_node(dt,SGBCDispersive)
#ifdef CompileWithOpenMP
!$OMP  PARALLEL DO DEFAULT(SHARED) private (conta) schedule(guided) 
#endif
do conta=1,malon%numnodes
   !call AdvanceSGBCE_single_node(malon%Nodes(conta), dt,SGBCDispersive)
   call AdvanceSGBCE_single_node(conta, dt,SGBCDispersive)
end do
#ifdef CompileWithOpenMP
!$OMP  END PARALLEL DO
#endif

contains
   !subroutine AdvanceSGBCE_single_node(compo,dt,SGBCDispersive) !esta tb seria valida
      !type (SGBCSurface_t),target   :: compo
   subroutine AdvanceSGBCE_single_node(conta,dt,SGBCDispersive)
      !non argument arguments
      integer (kind=4),  intent(IN) ::  conta
      REAL (KIND=RKIND), intent(IN) :: dt
      logical,           intent(IN) :: SGBCDispersive
      !local
      integer (kind=4)  ::  i
      type (SGBCSurface_t),pointer   :: compo
      type (MDfield_t) , pointer :: EDIS
      REAL (KIND=RKIND), pointer :: dDIS
      

      
      !do conta=1,malon%numnodes
      compo => malon%Nodes(conta)

!!!los extremos de los E internos
         if (compo%depth>0) then
!esto solo hace algo en el caso de yee.
            if (.not.compo%SGBCcrank)  then !yee
               if (compo%Correct_Ha) then
!!!!lOS EXTREMOS
!!!!ojo NO SE CORRIGEN LOS filo_placaS CON FDTD EN DISPERSIVOS. ASI QUE CABE ESPERAR QUE SEA INESTABLE CON DISPERSIVOS SGBCYEEE
                  compo%E( compo%depth  ) = compo%G1(1) *compo%E( compo%depth ) +  &
                                       (compo%G2a(1) *(compo%Ha_Plus   - compo%Hyee_Right) - compo%G2b(1) *(compo%Hb_Plus   - compo%Hb_Minu  ) )
                  !
                  compo%E(-compo%depth  ) = compo%G1(0) *compo%E(-compo%depth ) +  &
                                          (compo%G2a(0) *(compo%Hyee__Left - compo%Ha_Minu  ) - compo%G2b(0) *(compo%Hb_Plus   - compo%Hb_Minu  ) )

               elseif (compo%Correct_Hb) then
                  compo%E( compo%depth  ) = compo%G1(1) *compo%E( compo%depth ) +  &
                                          (compo%G2a(1) *(compo%Ha_Plus   - compo%Ha_Minu  ) - compo%G2b(1) *(compo%Hb_Plus   - compo%Hyee_Right ) )
                  compo%E(-compo%depth  ) = compo%G1(0) *compo%E(-compo%depth ) +  &
                                          (compo%G2a(0) *(compo%Ha_Plus   - compo%Ha_Minu  ) - compo%G2b(0) *(compo%Hyee__Left - compo%Hb_Minu ) )
                     
               endif
            endif !DEL MALONECRANK
         else !si tiene depth=0
            compo%E( compo%depth  ) = compo%G1(0) *compo%E( compo%depth ) +  (compo%G2a(0) *(compo%Ha_Plus   - compo%Ha_Minu    ) - compo%G2b(0) *(compo%Hb_Plus     - compo%Hb_Minu  ) )
         endif !DEL COMODEPTH
!!!los E internos FDTD1D
         if (compo%SGBCcrank)  then !por debajo de 2 no tiene ningun sentido
            do i=-compo%depth , compo%depth
               compo%E_past (i) = compo%E(i)
            end do
!esto no es necesario. ya se pasa por MPI
!!!pre-jav 310116 (solo estas lineas
            !compo%d( compo%depth  )      = compo%E( compo%depth  )
            !compo%d(-compo%depth  )      = compo%E(-compo%depth  )
!!post-jav 310116 boundaries 
!!!!ojo estos coeficientes no irian bien con materiales magneticos porque no tengo en cuenta los gm1 (en crank-nicolson de jav) 
            if (compo%Correct_Ha) then !!!compo%G2a(0) es igual a compo%rh1 !!!compo%G2a(1) es igual a compo%rhn
               i=compo%depth
               compo%d( i ) =    - compo%an * (compo%E_past(i-1)) + compo%rbn * compo%E_past(i ) + &
                                          (compo%G2a(1) *(compo%Ha_Plus   - compo%Hyee_Right) - compo%G2b(1) *(compo%Hb_Plus   - compo%Hb_Minu  ) )
               !
               !
               i=-compo%depth
               compo%d( i ) =    - compo%c1 * (compo%E_past(i+1)) + compo%rb1 * compo%E_past(i ) + &
                                          (compo%G2a(0) *(compo%Hyee__Left - compo%Ha_Minu  ) - compo%G2b(0) *(compo%Hb_Plus   - compo%Hb_Minu  ) )
               !
               !
            elseif (compo%Correct_Hb) then !!!-compo%G2b(0) es igual a compo%rh1 !!!-compo%G2b(1) es igual a compo%rhn
               i=compo%depth
               compo%d( i  ) =   - compo%an * (compo%E_past(i-1)) + compo%rbn * compo%E_past(i ) + &
                                          (compo%G2a(1) *(compo%Ha_Plus   - compo%Ha_Minu  ) - compo%G2b(1) *(compo%Hb_Plus   - compo%Hyee_Right ) )
               i=-compo%depth
               compo%d(i  ) =    - compo%c1 * (compo%E_past(i+1)) + compo%rb1 * compo%E_past(i ) + &
                    (compo%G2a(0) *(compo%Ha_Plus   - compo%Ha_Minu  ) - compo%G2b(0) *(compo%Hyee__Left - compo%Hb_Minu ) )
            endif
!fin nuevas boundaries jav 310116 (debe coincidir comentando lo anterior y dejando el pre-jav
            do i=-compo%depth+1 , compo%depth-1
!!!                  compo%d(i) = - compo%a(i) * (compo%E_past(i-1)+compo%E_past(i+1)) &!!!+ compo%rb(i) * compo%E_past(i ) + compo%rh(i) * (compo%H(i)-compo%H(i-1)) !! compo%a es igual a compo%c
!!!!   mmmmmmm esto no es cierto que a = c para multicapa porque es asimetrico 0121: corrijo lo anteiror a lo bruto
                compo%d(i)=-compo%a(i)*compo%E_past(i-1) &
                           -compo%c(i)*compo%E_past(i+1) &
                           +compo%rb(i)*compo%E_past(i)  &
                           +compo%rh  (i)*compo%H(i  )   & !!alobruto!!!0121
                           -compo%rhm1(i)*compo%H(i-1)    !!alobruto!!!0121
            end do
            
            !
            if (SGBCDispersive) then 
               do i=-compo%depth + 1 , compo%depth -1 !sin los putos filo_placas
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
               do i=-compo%depth + 1 , compo%depth -1 !sin los putos filo_placas
               EDIS=>compo%EDis(i)
               call segundo_CNAdvanceSGBCDispersive (compo%E(i),EDIS,compo%numpolres,compo%G3,Compo%kappa,compo%beta,dt)
               end do
            end if
         else !YEE
            do i=-compo%depth+1 , compo%depth-1
               compo%E(i) = compo%G1_interno (i)  *compo%E(i) + compo%G2_interno (i)  *( compo%H(i) - compo%H(i-1) )               
               ! NO SE CORRIGIERON LOS PUTOS filo_placaS ANTES. AQUI SOLO SE CORRIGE EL INTERIOR
               if (SGBCDispersive) then 
                  EDIS=>compo%EDis(i)
                  call YeeAdvanceSGBCDispersive (EDIS,compo%numpolres,compo%G3,Compo%kappa,compo%beta,dt)
               endif
            end do
                                    
         endif
!!!los H internos 1D 
         if (compo%SGBCcrank)  then
            do i=-compo%depth , compo%depth-1
               compo%H(i) = compo%GM1_interno(i) *compo%H(i) + compo%GM2_interno(i)   /2.0_RKIND *( compo%E(i+1)      - compo%E(i)      + &
                                                                                               compo%E_past(i+1) - compo%E_past(i) )
            end do                 
!for crank-nicolson a half-step advance is necessary for the H field since E and H are synchronous. Coexisten por tanto dos campos H en (-compo%depth) y (compo%depth-1): uno en n y otro en n+1/2
!!solo pre-jav 310116 (solo estas lineas
!!!!esta correccion es analiticamente correcta pero fuente de posibles intestabilidades pues es un Yee. Lo comento 02115 y dejo una aprox backwards que es mas que suficiente para metales donde la velocidad es muy baja comparada con el vacio
!            compo%Hyee__Left = compo%GM1_interno *compo%Hyee__Left + compo%GM2_interno *( compo%E(-compo%depth+1) - compo%E(-compo%depth  ) )
!            compo%Hyee_Right = compo%GM1_interno *compo%Hyee_Right + compo%GM2_interno *( compo%E(   compo%depth) - compo%E( compo%depth-1) )
!!!!post-jav 310116 boundaries solo lo que sigue es correcto. Ya no es una aproximacion TB sino el resultado del nuevo CN_jav
            compo%Hyee__Left = compo%H(-compo%depth  ) 
            compo%Hyee_Right = compo%H( compo%depth-1)              
            
         else !yee
            IF (compo%depth/=0) THEN
               do i=-compo%depth , compo%depth-1
                  compo%H(i) = compo%GM1_interno(i) *compo%H(i) +  compo%GM2_interno(i) *( compo%E(i+1) - compo%E(i) )  !E y H se usan reciprocamente siempre con el mismo signo
               end do
               compo%Hyee__Left = compo%H( -compo%depth)
               compo%Hyee_Right = compo%H(compo%depth-1)
            ENDIF
         endif
!!!!copio en su Efield el promedio a efectos de peticion de sondas y calculo en los bordoes (luego al avanzar el H el principal de ambos lados, utilizara este Efield, pero el advanceSGBCH lo recorregira por el correcto
         compo%Efield =(compo%E(-compo%depth)+compo%E(compo%depth))/2. 
      !end do
   end subroutine AdvanceSGBCE_single_node
end subroutine AdvanceSGBCE

subroutine AdvanceSGBCH
   integer (kind=4)  ::  conta
   type (SGBCSurface_t), pointer :: compo
   character(len=BUFSIZE) :: buFF
   !NOTE: Esto no se puede optimizar 
   !      porque los dos o mas compo%H{a,b} pueden apuntar
   !      al mismo campo se puede produce un conflicto de acceso
   do conta=1,malon%numnodes
      compo => malon%Nodes(conta)
!!!!ojo: es una correccion a lo que hace el principal utilizando el campo electrico correcto
      if (compo%Correct_Ha) then
         compo%Ha_Plus = compo%Ha_Plus +  compo%gm2_externo* (compo%Efield - compo%E( compo%depth)) !insisto: es una correccion: el principal ha aniadido/quitado Efield y debe quitar/aniadir E del extremo correspondiente
         compo%Ha_Minu = compo%Ha_Minu -  compo%gm2_externo* (compo%Efield - compo%E(-compo%depth))
      elseif (compo%Correct_Hb) then                                                       
         compo%Hb_Plus = compo%Hb_Plus -  compo%gm2_externo* (compo%Efield - compo%E( compo%depth))
         compo%Hb_Minu = compo%Hb_Minu +  compo%gm2_externo* (compo%Efield - compo%E(-compo%depth))
      else     
         WRITE (buff, *)    'Buggy ERROR: In SGBCs. '
         CALL StopOnError (0,0,buff)
      endif
   end do
   return
end subroutine AdvanceSGBCH

subroutine calc_g1g2gm1gm2_compo(sgg,compo,eps00,mu00,SGBCDispersive)
   REAL (KIND=RKIND), intent(IN)           ::  eps00,mu00
   Complex (Kind=CKIND), pointer, dimension ( : )      ::  Beta,Kappa,G3
!!!!      
   type (SGGFDTDINFO), intent(IN) ::  sgg
   type (SGBCSurface_t), pointer, intent(INOUT) :: compo
   character(len=BUFSIZE) :: buff
!!!variables locales
 real (kind=RKIND) :: width,sigmatemp,eprtemp,sigmamtemp,murtemp,epsilon,sigma,mu,sigmam,g1,g2,gm1,gm2,delta_entreEinterno_temp,epr_adyacentei,sig_adyacentei
 real (kind=RKIND), dimension(0:1) :: epr_adyacente,sig_adyacente
 integer (kind=4) :: i,ib,ib_ady
 logical :: SGBCDispersive
 eps0=eps00; mu0=mu00; !chapuz para convertir la variables de paso en globales
   


   IF (compo%depth==0) then
         !!!!!!! compo%delta_entreeinterno=0.0 !!no se usa nunca !ojoooo revisar caso 0121 
         !!!!!!!averagefactor  = width / compo%transversaldeltah /factor !!! sgg promediados buenos filo_placas 201115
         !!!!!!!epsilon = (1.0_rkind - averagefactor ) * (epr_adyacente(0)+epr_adyacente(1))/2.0_rkind  * eps0 + &
         !!!!!!!                       averagefactor   * eprtemp           * eps0
         !!!!!!!sigma =   (1.0_rkind - averagefactor ) * (sig_adyacente(0)+sig_adyacente(1))/2.0_rkind  + &
         !!!!!!!                       averagefactor   * sigmatemp
         compo%delta_entreEinterno=0.0 !!no se usa nunca
         do i=0,1
 !0121 Voy a tomar vacio PORQUE EN LAS ESQUINAS ENTRE sgbc DETECTA MAL EL MEDIO ADYACENTE. sE HACIA ANTES DE 0121 ASI TAMBIEN
 !050421 lo devuelvo a vacio porque no acabo de ver el caso de las esquinas entre SGBC
           !  epr_adyacente(i) = Sgg%Med(compo%med(i))%epr   
           !  sig_adyacente(i) = Sgg%Med(compo%med(i))%sigma 
           !  if (((epr_adyacente(i)-1.0_RKIND>1E-3).OR.(sig_adyacente(i)>1E-3)).AND.(.NOT.(SGG%Med(compo%med(i))%Is%SGBC)) )  then 
           !           WRITE (buff, *)    '(WARNING) Collision of composite with non free-space medium. Assuming free-space, instead ', compo%med(i), epr_adyacente(i), sig_adyacente(i)
                      epr_adyacente(i) = 1.0_RKIND
                      sig_adyacente(i) = 0.0_RKIND
           !           CALL WarnErrReport (buff,.false.)
           !  endif
         end do
         width=sgg%med(compo%jmed)%Multiport(1)%width(1)
         sigmatemp=sgg%Med(compo%jmed)%multiport(1)%sigma(1)
         eprtemp= sgg%Med(compo%jmed)%multiport(1)%epr(1)   
       !prescindo de los filo_placa 0121 poque en las pec de boundaries las detecta incorrectamente !de todos modos esto nunca me ha gustado 0121
         !no puedo prescindir de los filo_placas a 040523 SinSTOCH_antiguou_th0.0001
         if (compo%es_unfilo_placa) then
             epsilon = ((epr_adyacente(0)+epr_adyacente(1))/2.0_rkind  * eps0 *(compo%transversaldeltah - width/2.0_rkind)   + &
                         eprtemp                                       * eps0 * width /2.0_rkind) / &
                        (compo%transversaldeltah)
             sigma =   ((sig_adyacente(0)+sig_adyacente(1))/2.0_rkind         *(compo%transversaldeltah - width/2.0_rkind)   + &
                         sigmatemp                                            *width  /2.0_rkind) / &
                        (compo%transversaldeltah)
         else
         epsilon = ((epr_adyacente(0)+epr_adyacente(1))/2.0_rkind  * eps0 *(compo%transversaldeltah-width)   + &
                     eprtemp                                       * eps0 *width             ) / &
                        (compo%transversaldeltah)
         sigma =   ((sig_adyacente(0)+sig_adyacente(1))/2.0_rkind         *(compo%transversaldeltah-width)   + &
                     sigmatemp                                            *width             ) / &
                        (compo%transversaldeltah)
         endif
!!!!!ajusta primero los g1 y g2 de los bordes del espesor de la capa !ojo en sgbcdispersive no se utilizan las constantes kappa, beta, g3 en lo putos filo_placas. solo en el interior
         call g1g2(sgg%dt,epsilon,sigma,g1,g2)
         compo%g1   (0)=g1 
         if (compo%correct_ha) then
             compo%g2a(0)= g2 / compo%transversaldeltah 
             compo%g2b(0)= g2 / compo%alignedldeltah
 !!!!los indices (1) no se utilizan para el caso particular depth=0
             compo%g2a(1)= compo%g2a(0)
             compo%g2b(1)= compo%g2b(0)
         elseif (compo%correct_hb) then
             compo%g2a(0)= g2 / compo%alignedldeltah
             compo%g2b(0)= g2 / compo%transversaldeltah
             compo%g2a(1)= compo%g2a(0)
             compo%g2b(1)= compo%g2b(0)
         endif
!!!!!!!!!!!!ahora el interior  
        epsilon =  eprtemp * eps0
        sigma =    sigmatemp 
        if (sgbcdispersive) then
             beta => compo%beta %val
             kappa => compo%kappa %val
             g3 => compo%g3 %val
             call g1g2_dispersive(sgg%dt,epsilon,sigma,g1,g2,beta,kappa,g3,compo%numpolres,compo%a11,compo%c11)
         else
             call g1g2(sgg%dt,epsilon,sigma,g1,g2)
         endif
         compo%g1_interno(0)=g1
         compo%g2_interno(0)=g2
         sigmamtemp= 0.
         murtemp=    1.
         mu=murtemp * mu0
         Sigmam = sigmamtemp
         call gm1gm2(sgg%dt,mu,sigmam,gm1,gm2)
         compo%gm1_interno(0)=gm1
         compo%gm2_interno(0)=gm2
 else !DEL MALONYEDEPTH
!PRIMERO las constantes en las fronteras en la dimension del espesor
     do i=0,1
 !0121 Voy a tomar vacio PORQUE EN LAS ESQUINAS ENTRE sgbc DETECTA MAL EL MEDIO ADYACENTE. sE HACIA ANTES DE 0121 ASI TAMBIEN
 !050421 lo devuelvo a vacio porque no acabo de ver el caso de las esquinas entre SGBC
         ! epr_adyacente(i) = Sgg%Med(compo%med(i))%epr   
         ! sig_adyacente(i) = Sgg%Med(compo%med(i))%sigma 
         ! if (((epr_adyacente(i)-1.0_RKIND>1E-3).OR.(sig_adyacente(i)>1E-3)).AND.(.NOT.(SGG%Med(compo%med(i))%Is%SGBC)) )  then 
         !         WRITE (buff, *)    '(WARNING) Collision of composite with non free-space medium. Assuming free-space, instead ', compo%med(i), epr_adyacente(i), sig_adyacente(i)
                  epr_adyacente(i) = 1.0_RKIND
                  sig_adyacente(i) = 0.0_RKIND
         !         CALL WarnErrReport (buff,.false.)
         ! endif
         if (i==0) then 
             ib=1 !primera capa
             delta_entreEinterno_temp=compo%delta_entreEinterno(-compo%depth)
         else
             ib=sgg%Med(compo%jmed)%multiport(1)%numcapas !ultima capa
             delta_entreEinterno_temp=compo%delta_entreEinterno(compo%depth-1)     
         endif
         width=sgg%med(compo%jmed)%Multiport(1)%width(ib)
         sigmatemp= sgg%Med(compo%jmed)%multiport(1)%sigma(ib)
         eprtemp=   sgg%Med(compo%jmed)%multiport(1)%epr(ib)  
       !prescindo de los filo_placa 0121 poque en las pec de boundaries las detecta incorrectamente !de todos modos esto nunca me ha gustado 0121
                     !no puedo prescindir de los filo_placas a 040523 SinSTOCH_antiguou_th0.0001
         if (compo%es_unfilo_placa) then
             epsilon = (epr_adyacente(i)* eps0 *(compo%transversalDeltaH + delta_entreEinterno_temp /2.0_RKIND)   + &
                                eprtemp         * eps0 *               (delta_entreEinterno_temp /2.0_RKIND)) / &
                               (compo%transversalDeltaH +               delta_entreEinterno_temp)
             Sigma =   (sig_adyacente(i)       *(compo%transversalDeltaH + delta_entreEinterno_temp /2.0_RKIND)   + &
                                sigmatemp               *(delta_entreEinterno_temp  /2.0_RKIND)) / &
                               (compo%transversalDeltaH + delta_entreEinterno_temp)
         else
             epsilon = (epr_adyacente(i)* eps0 *compo%transversalDeltaH   + &
                        eprtemp         * eps0 *delta_entreEinterno_temp    ) / &
                       (compo%transversalDeltaH + delta_entreEinterno_temp)                
             Sigma =   (sig_adyacente(i)       *compo%transversalDeltaH   + &
                        sigmatemp              *delta_entreEinterno_temp    ) / &
                       (compo%transversalDeltaH + delta_entreEinterno_temp)
             
         endif

       
         
         
 !ajusta primero los g1 Y G2      !no preciso gm1 ni gm2 en los filo_placas
         call g1g2(sgg%dt,epsilon,sigma,g1,g2)
         compo%g1(i)=g1 
         if (compo%Correct_Ha) then
             compo%G2a(i)= G2 / (0.5_RKIND * compo%transversalDeltaH + 0.5_RKIND*delta_entreEinterno_temp)
             compo%G2b(i)= G2 / compo%alignedlDeltaH
         elseif (compo%Correct_Hb) then
             compo%G2a(i)= G2 / compo%alignedlDeltaH
             compo%G2b(i)= G2 / (0.5_RKIND * compo%transversalDeltaH + 0.5_RKIND*delta_entreEinterno_temp)
         endif
     end do !DEL BARRIDO i=0,1 De las dos fronteras de la lamina en la dimension del espesor
     
 !ahora el interior !El primero y ultimo G no se usa 0121
     compo%G2_interno =2e31 !valores default absurdos para detectar errores
     compo%G1_interno =-2e21 
     barridoporcapas: do i=-compo%depth+1,compo%depth-1   !0121
         ib=compo%capa(i)         
         ib_ady=compo%capa(i-1)         
         if ((ib<1).or.(ib>sgg%Med(compo%jmed)%multiport(1)%numcapas)) then
             WRITE (buff, *)   'Buggy error in ib fuera de rango en compo numcapas. Contact '
             CALL StopOnError (0,0,buff)
             stop
         endif
         eprtemp=    sgg%Med(compo%jmed)%multiport(1)%epr(ib)  
         sigmatemp=  sgg%Med(compo%jmed)%multiport(1)%sigma(ib)
         epr_adyacentei = sgg%Med(compo%jmed)%multiport(1)%epr(ib_ady)
         sig_adyacentei = sgg%Med(compo%jmed)%multiport(1)%sigma(ib_ady)
         !!! Promedio interpolatoriamente bien 0121
         eprtemp = (epr_adyacentei     * compo%delta_entreEinterno(i-1)   + &
                    eprtemp            * compo%delta_entreEinterno(i) ) / &
                                        (compo%delta_entreEinterno(i-1) + compo%delta_entreEinterno(i))
         sigmatemp =   (sig_adyacentei  *compo%delta_entreEinterno(i-1)   + &
                        sigmatemp       *compo%delta_entreEinterno(i) ) / &
                                        (compo%delta_entreEinterno(i-1) + compo%delta_entreEinterno(i))
         !!!
         epsilon =  eprtemp * eps0
         Sigma =    sigmatemp 
         call g1g2(sgg%dt,epsilon,sigma,g1,g2)
         compo%g1_interno(i)=g1
         compo%g2_interno(i)=g2 /((compo%delta_entreEinterno (i)+compo%delta_entreEinterno (i-1))/2.0_RKIND) !semisuma diff  no centrada entre capas
     end do barridoporcapas
     
     compo%GM2_interno=-1e30 !valores default absurdos para detectar errores
     compo%GM1_interno=3e22 
 !ahora el interior !el ultimo GM que no se usa 0121
     barridoporcapasH: do i=-compo%depth,compo%depth-1   !0121
         ib=compo%capa(i)         
         if ((ib<1).or.(ib>sgg%Med(compo%jmed)%multiport(1)%numcapas)) then
             WRITE (buff, *)   'Buggy error in ib fuera de rango en compo numcapas. '
             CALL StopOnError (0,0,buff)
             stop
         endif
         sigmamtemp= sgg%Med(compo%jmed)%multiport(1)%sigmam(ib)
         murtemp=    sgg%Med(compo%jmed)%multiport(1)%mur(ib)
         mu=murtemp * mu0
         Sigmam = sigmamtemp
         call gm1gm2(sgg%dt,mu,sigmam,gm1,gm2)
         compo%gm1_interno(i)=gm1
         compo%gm2_interno(i)=gm2  /compo%delta_entreEinterno (i) !hay que no hay que semisumar porque es interno
     end do barridoporcapasH
     
     
 endif !del compodepth

   return
end subroutine calc_g1g2gm1gm2_compo

!!!!!!!
subroutine g1g2(dt,epsilon,sigma,G1,G2)
   real (kind=RKIND_tiempo), intent(in) :: dt
   real (kind=RKIND), intent(in) :: epsilon,sigma
   real (kind=RKIND), intent(out) :: g1,g2

   G1=(1.0_RKIND  - Sigma * dt / (2.0_RKIND * epsilon ) ) / &
      (1.0_RKIND  + Sigma * dt / (2.0_RKIND * epsilon ) )
   G2=dt / epsilon                       / &
      (1.0_RKIND  + Sigma * dt / (2.0_RKIND * epsilon))

   if (g1 < 0.0_RKIND) then !exponential time stepping
      g1=exp(- Sigma * dt / (epsilon ))
      g2=(1.0_RKIND-g1)/ Sigma
   else
      continue
   endif   
   return
end subroutine g1g2

!!!!!!!
subroutine gm1gm2(dt,mu,sigmam,Gm1,Gm2)
   real (kind=RKIND_tiempo), intent(in) :: dt
   real (kind=RKIND), intent(in) :: mu,sigmam
   real (kind=RKIND), intent(out) :: gm1,gm2

   Gm1=(1.0_RKIND  - Sigmam * dt / (2.0_RKIND * mu) ) / &
      (1.0_RKIND  + Sigmam * dt / (2.0_RKIND * mu ) )
   Gm2=dt / mu                       / &
      (1.0_RKIND  + Sigmam * dt / (2.0_RKIND * mu))

   if (gm1 < 0.0_RKIND) then !exponential time stepping
      gm1=exp(- Sigmam * dt / (mu ))
      gm2=(1.0_RKIND-gm1)/ Sigmam
   else
      continue
   endif
   return
end subroutine gm1gm2

!!!!!!! medios dispersivos sgg 12/05/16 
subroutine g1g2_Dispersive(dt,epsilon,sigma,G1,G2,Beta,Kappa,G3,numpolres,a11,c11)
   real (kind=RKIND), intent(in) :: epsilon,sigma
   real (kind=RKIND_tiempo), intent(in) :: dt
   real (kind=RKIND), intent(out) :: g1,g2
   COMPLEX (kind=ckind), intent(in), allocatable, dimension(:) :: a11, c11
   integer (kind=4) :: numpolres, i1
   real (kind=RKIND) :: tempo
!!!SGBC dispersivos 12/05/16
   Complex (Kind=CKIND), pointer, dimension ( : )      ::  Beta,Kappa,G3

     do i1=1,numpolres
         Kappa(i1) =(1.0_RKIND + a11(i1)*dt/2.0_RKIND)/&
                    (1.0_RKIND - a11(i1)*dt/2.0_RKIND)
         Beta(i1)=  (C11(i1)*dt) /&
                    (1.0_RKIND - a11(i1)*dt/2.0_RKIND)
     end do
     tempo=0.0_RKIND
     Do i1=1,NumPolRes
         tempo=tempo+REAL (Beta(i1))
     end do
     G1=                        (2.0_RKIND * epsilon + tempo - sigma*dt) / & !ojo No estes tentado de cambiar este signo. Cuadra con han dutton 130516 y con edispersives 
                                (2.0_RKIND * epsilon + tempo + sigma*dt)
     G2=         2.0_RKIND *dt/ (2.0_RKIND * epsilon + tempo + sigma*dt)
!!!! aqui no cabe exponential time stepping
     Do i1=1,NumPolRes
         G3(i1)=G2/2.0_RKIND * (1.0_RKIND+Kappa(i1))
     end do
   return
end subroutine g1g2_Dispersive

subroutine StoreFieldsSGBCs(stochastic)
      integer (kind=4) :: conta,i,k1
      logical :: SGBCDispersive,stochastic
      type (SGBCSurface_t), pointer :: compo
      do conta=1,malon%numnodes
         compo => malon%Nodes(conta)
         write(14,err=634) (compo%E     (i),i=-compo%depth,compo%depth)

         if (compo%SGBCcrank)  then 
             write(14,err=634) (compo%E_past(i),i=-compo%depth,compo%depth)
         endif
         write(14,err=634) compo%Hyee__left
         write(14,err=634) compo%Hyee_right
         write(14,err=634) (compo%H     (i),i=-compo%depth,compo%depth-1)
         !
         if (malon%SGBCDispersive) then 
             write(14,err=634) (compo%EDis(i)%fieldPrevious, i=-compo%depth,compo%depth)
             Do k1=1,compo%NumPolRes
                write(14,err=634) (compo%EDis(i)%current(k1), i=-compo%depth,compo%depth)
             enddo
         endif
         
      end do

   goto 635
634   call print11(0,SEPARADOR//separador//separador)
   call print11(0,'SGBC: ERROR WRITING RESTARTING FIELDS. IGNORING AND CONTINUING')
   call print11(0,SEPARADOR//separador//separador)          
635   return
end subroutine StoreFieldsSGBCs


subroutine DestroySGBCs(sgg)

   type (SGGFDTDINFO), intent(INOUT)         ::  sgg
   integer (kind=4)  ::  i,conta

   !free up memory
   do i=1,sgg%NumMedia
      if (allocated(malon%mediosDis)) then
          if (allocated(malon%mediosDis(i)%a11)) deallocate(malon%mediosDis(i)%a11)
          if (allocated(malon%mediosDis(i)%c11)) deallocate(malon%mediosDis(i)%c11)
      endif
      if ((sgg%Med(i)%Is%SGBC).and.(.not.sgg%Med(i)%Is%PML))  deallocate (sgg%Med(i)%Multiport)      
   end do
   if (allocated(malon%mediosDis)) deallocate(malon%mediosDis)
   !
   do conta=1,malon%numnodes
      if (allocated(malon%Nodes(conta)%d))  deallocate (malon%Nodes(conta)%d) !AUXILIAR DE CRANK-NICOLSON
      if (allocated(malon%Nodes(conta)%beta%val))  deallocate (malon%Nodes(conta)%beta%val) !AUXILIAR DE CRANK-NICOLSON dispersivo
      if (allocated(malon%Nodes(conta)%kappa%val))  deallocate (malon%Nodes(conta)%kappa%val) !AUXILIAR DE CRANK-NICOLSON dispersivo
      if (allocated(malon%Nodes(conta)%G3%val))  deallocate (malon%Nodes(conta)%G3%val) !AUXILIAR DE CRANK-NICOLSON dispersivo
      if (allocated(malon%Nodes(conta)%Edis))  deallocate (malon%Nodes(conta)%Edis) !AUXILIAR DE CRANK-NICOLSON dispersivo
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

   if(allocated(malon%nodes)) deallocate (malon%nodes)
end subroutine


subroutine test_stab(G2,GM2)
   REAL (KIND=RKIND)     , pointer, dimension ( : )   ::   g2, gm2
   integer (kind=4)  ::  conta
   logical :: unstable
   type (SGBCSurface_t), pointer :: compo
   REAL (KIND=RKIND) :: heur
   character(len=BUFSIZE) :: buff

   heur=1.0_RKIND/sqrt(3.0_RKIND)
   unstable = .false.
   do conta=1,malon%numnodes
      compo => malon%Nodes(conta)
!!!los extremos de los E internos
      unstable= unstable.or. &
             (G2(compo%jmed) * Gm2(compo%jmed)  > heur) .or. &
             (compo%G2a(1)   * Gm2(compo%jmed)  > heur) .or. &
             (compo%G2b(1)   * Gm2(compo%jmed)  > heur) .or. &
             (compo%G2a(0)   * Gm2(compo%jmed)  > heur) .or. &
             (compo%G2a(0)   * Gm2(compo%jmed)  > heur)
   end do

   if (unstable) then
        WRITE (buff, *)    'ERROR: SGBCs may become unstable. Reduce cfl'
        CALL WarnErrReport (buff,.true.)
   endif

   return

end subroutine test_stab

subroutine depth(compo,sgg,jmed,SGBCFreq,SGBCresol,SGBCdepth) 
 type (SGGFDTDINFO), intent(IN)     ::  sgg
 real (kind=rkind) :: SGBCFreq,SGBCresol,sigma, epr,epsilon,skin_depth,width,widthtotal
 integer (kind=4) :: jmed,i,SGBCdepth,numcapas,precuenta,celdafinal,celdainicial,anchocapa
 integer (kind=4) , pointer, dimension(:) :: capa
 logical :: ultimacapamas1
 character(len=BUFSIZE) :: buff
    
 type (SGBCSurface_t), pointer :: compo

!!!0121 multicapas
 numcapas = sgg%Med(jmed)%multiport(1)%numcapas
 compo%depth=0
 do precuenta=0,1
     if (precuenta==1) then
         if (mod(compo%depth,2)/=0) then
             compo%depth=compo%depth+1 !redondea el numero de capas total a un numero par
!rellena el sobrante con la ultima capa
             ultimacapamas1=.true.
         else 
             ultimacapamas1=.false.
         endif
         compo%depth=int(compo%depth/2.0_RKIND) !divide por 2 porque arranca en -compo%depth y llega a +compo%depth 
         if (compo%depth>0) then
             if (.not.allocated(compo%capa))                allocate (compo%capa(-compo%depth:compo%depth-1))
             if (.not.allocated(compo%delta_entreEinterno)) allocate (compo%delta_entreEinterno(-compo%depth:compo%depth-1))
         else
             if (.not.allocated(compo%capa))                allocate (compo%capa(0:0))
             if (.not.allocated(compo%delta_entreEinterno)) allocate (compo%delta_entreEinterno(0:0))
         endif
         
         celdafinal=-compo%depth-1
     endif
     widthtotal=0.; width=0.; sigma=0.; epr=0.; 
     do i=1,numcapas
         width=      sgg%Med(jmed)%multiport(1)%width(i) 
         sigma=      sgg%Med(jmed)%multiport(1)%sigma(i) 
         epr=        sgg%Med(jmed)%multiport(1)%epr(i)  
         epsilon=epr * eps0
         widthtotal=widthtotal +     sgg%Med(compo%jmed)%multiport(1)%width(i) 
         skin_depth=1.0_RKIND / (Sqrt(2.0_RKIND)*SGBCFreq*Pi*(Mu0**2*(4*Epsilon**2.0_RKIND + Sigma**2/(SGBCFreq**2*Pi**2.0_RKIND )))**0.25_RKIND * &
                                 Sin(atan2(2*Pi*Epsilon*Mu0, -(Mu0*Sigma)/SGBCFreq)/2.0_RKIND))
         if (SGBCdepth==0) then !numcapas debe ser necesariamente 1
             if (numcapas > 1) then
                WRITE (buff, *)   'SGBCDepth=0 and numcapas>1 not compatible. Please, relaunch'
                CALL StopOnError (0,0,buff)
             else
                 anchocapa=1 !numcapas es necesariamente 1 si continua
             endif
         elseif (SGBCdepth>0) then
             anchocapa=SGBCdepth
         else !si es negativo se calcula con la resol
             anchocapa=1+int(SGBCresol*width/skin_depth)
         endif
         if (anchocapa<2) anchocapa=2 !es razonable no dejarlo nunca en 1
         !fin niapas
         if (precuenta==0) then 
             if (SGBCDepth==0) then 
                 compo%depth=0
             else
                 compo%depth=compo%depth+anchocapa
             endif
         elseif (precuenta==1) then
             if (SGBCDepth==0) then      !!!bug corregido a 040523
                     celdainicial=0
                     celdafinal=0
                     anchocapa=1
                     compo%capa(celdainicial:celdafinal) = i
                     compo%delta_entreEinterno(celdainicial:celdafinal)=width/anchocapa
                     continue
             else  
             celdainicial=celdafinal+1
             celdafinal=celdainicial+anchocapa-1
             if ((i==numcapas).and.ultimacapamas1) then
!rellena el sobrante con la ultima capa si no es una division cabal
                     anchocapa=anchocapa+1
                     celdafinal=celdafinal+1
             endif
             compo%capa(celdainicial:celdafinal) = i
             compo%delta_entreEinterno(celdainicial:celdafinal)=width/anchocapa
             continue
         endif
         endif
     end do
     if (precuenta==1) then
         if ((celdafinal/=compo%depth-1).and.(compo%depth/=0)) then
                WRITE (buff, *)   'Buggy error redondeo final ultima capa. '
                CALL StopOnError (0,0,buff)
         endif
     endif
 end do
!!!end 02121
 return
end subroutine depth

function GetSGBCs() result(r)
   type(Malon_t), pointer  ::  r
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
 real (kind=RKIND) ,intent(in),dimension(n) :: aa,bb,cc
 real (kind=RKIND) ,intent(in)   :: a1,b1,c1,an,bn,cn
 real (kind=RKIND) ,dimension(1:n) :: a,b,c
   real (kind=RKIND) ,dimension(n),intent(in) :: d
   real (kind=RKIND) ,dimension(n),intent(out) :: x
   real (kind=RKIND) ,dimension(n) :: cp,dp
   real (kind=RKIND)  :: m
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
   enddo
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
   real (kind=RKIND) ,intent(in) :: aa,bb,cc,a1,b1,c1,an,bn,cn
   real (kind=RKIND) ,dimension(n) :: a,b,c
   real (kind=RKIND) ,dimension(n),intent(in) :: d
   real (kind=RKIND) ,dimension(n),intent(out) :: x
   real (kind=RKIND) ,dimension(n) :: cp,dp
   real (kind=RKIND)  :: m
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
   enddo
   ! initialize x
   x(n) = dp(n)
   ! solve for x from the vectors c-prime and d-prime
   do i = n-1, 1, -1
      x(i) = dp(i)-cp(i)*x(i+1)
   end do
   return
   end subroutine solve_tridiag_iguales           


end module SGBC_nostoch

