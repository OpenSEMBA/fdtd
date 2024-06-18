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
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Observation module to store the observed data
!  Creation date Date :  April, 8, 2010
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module Observa
   use fdetypes

#ifdef CompileWithMPI
   use MPIcomm
#endif

#ifdef CompileWithWires
   use wiresHolland_constants
   use HollandWires
   use Wire_bundles_mtln_mod
   use mtln_solver_mod , mtln_solver_t => mtln_t 
#endif
#ifdef CompileWithBerengerWires
   use WiresBerenger
#endif   
#ifdef CompileWithSlantedWires
   use WiresSlanted
   use WiresSlanted_Types
   use WiresSlanted_Constants
#endif
   use report

#ifdef CompileWithNF2FF
   use farfield_m
#endif
#ifdef CompileWithNodalSources
   use nodalsources
#endif
!
   IMPLICIT NONE
   private

   type Serialized_t
      REAL (KIND=RKIND), pointer, dimension(:,:)   ::  valor,valor_x,valor_y,valor_z ! (step, valor)
      INTEGER (kind=4), POINTER, DIMENSION(:) :: eI,eJ,eK,currentType,sggmtag
      complex( kind = CKIND), dimension( :,:), allocatable  :: valorComplex_x,valorComplex_y,valorComplex_z
   end type Serialized_t
   type item_t
#ifdef CompileWithWires
      type (CurrentSegments), pointer  ::  segmento !segmento de hilo que se observa si lo hubiere
#endif
#ifdef CompileWithBerengerWires
      type (TSegment)       , pointer  ::  segmento_Berenger !segmento de hilo que se observa si lo hubiere
#endif
#ifdef CompileWithSlantedWires
      class(Segment)        , pointer  ::  segmento_Slanted !segmento de hilo que se observa si lo hubiere
#endif
      character (LEN=BUFSIZE)  ::  path
      integer (kind=4) :: unit,unitmaster !to store the unit of the file y en caso de singlefileginario el unitmaster que escribe
      integer (kind=4) :: columnas !number of columns in the output file
      REAL (KIND=RKIND), pointer, dimension( : )   ::  valor,valor2,valor3,valor4,valor5 !stored values at each time step !not read but calculate !210521 also store -edl+vdrop
      REAL (KIND=RKIND)    ::  valorsigno !just to store the sign of the current for the wires
      REAL (KIND=RKIND), pointer, dimension( :, :, :, : )   ::  valor3D !stored values at each time step !not read but calculate
      type (Serialized_t)  ::  Serialized !para almecenar valores serializados en volumenes en vez de bulk
      !freqdomain probles
      complex( kind = CKIND), dimension( :, :, :, :,: ), allocatable  :: valor3DComplex !freqdomain probes
#ifdef CompileWithMPI
      integer (kind=4)      ::  MPISubcomm,MPIRoot,MPIGroupIndex
      integer (kind=4)      :: ZIorig,ZEorig
#endif
      integer (kind=4)      :: Xtrancos,Ytrancos,Ztrancos
      integer (kind=4)      :: XItrancos,YItrancos,ZItrancos
      integer (kind=4)      :: XEtrancos,YEtrancos,ZEtrancos
   end type

   type output_t
      type (item_t), dimension( : ), pointer  ::  item !path con el output y sus valores
      integer (kind=4)      ::  Trancos
      logical  ::  SaveAll
      integer (kind=4) :: TimesWritten !to control the volumic probes
      !freqdomain probles
      integer (KIND=4) :: NumFreqs
      real (kind=Rkind) , dimension( :), allocatable  :: Freq
      real (kind=Rkind) :: InitialFreq,FinalFreq,FreqStep
      complex( kind = CKIND), dimension( :), allocatable  :: auxExp_E,auxExp_H,dftEntrada   !para sondas freqdomain
   end type output_t



#ifdef CompileWithWires
   type(Thinwires_t), pointer  ::  Hwireslocal
#endif
#ifdef CompileWithBerengerWires
   type(TWires)     , pointer  ::  Hwireslocal_Berenger
#endif
#ifdef CompileWithSlantedWires
   type(WiresData)  , pointer  ::  Hwireslocal_Slanted
#endif

#ifdef CompileWithMPI
   REAL (KIND=RKIND), pointer, dimension( : )   ::  valores ,newvalores  !auxiliary for Bloque currents sync
#endif
!!!variables globales del modulo
   REAL (KIND=RKIND), save           ::  eps0,mu0
!!!
   !!!!!!!!!variables local

   REAL (KIND=RKIND)     , pointer, dimension ( : ), save  ::  InvEps ,InvMu
   type (output_t)       , pointer, dimension ( : ), save  ::  output

   public InitObservation,FlushObservationFiles,UpdateObservation,DestroyObservation,CloseObservationFiles,unpacksinglefiles, &
   GetOutput
   public output_t,item_t,Serialized_t,dtft
   public FlushMTLNObservationFiles
contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Initializes observation stuff
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine InitObservation(sgg,sggMiEx,sggMiEy,sggMiEz,sggMiHx,sggMiHy,sggMiHz,sggMtag,&
                              ThereAreObservation,ThereAreWires,ThereAreFarFields,resume,initialtimestep, finaltimestep,lastexecutedtime, &
                              nEntradaRoot,layoutnumber,size, saveall, b,singlefilewrite,wiresflavor, &
                              SINPML_fullsize,facesNF2FF,NF2FFDecim,eps00,mu00,simu_devia,mpidir,niapapostprocess)
      type (SGGFDTDINFO), intent(IN)         ::  sgg
      INTEGER (KIND=IKINDMTAG), intent(in) :: sggMtag  (sgg%Alloc(iHx)%XI:sgg%Alloc(iHx)%XE, sgg%Alloc(iHy)%YI:sgg%Alloc(iHy)%YE, sgg%Alloc(iHz)%ZI:sgg%Alloc(iHz)%ZE)
      logical :: simu_devia,niapapostprocess
      REAL (KIND=RKIND)           ::  eps00,mu00
      !---------------------------> inputs <----------------------------------------------------------
      integer (kind=4), intent(in) :: layoutnumber,size,mpidir
      type (nf2ff_t) :: facesNF2FF
      type (limit_t), dimension(1:6), intent(in)  ::  SINPML_fullsize
      type( bounds_t), intent( IN)  ::  b !needed by far field
      integer(kind = INTEGERSIZEOFMEDIAMATRICES), dimension ( 0 : b%sggMiEx%NX-1 , 0 : b%sggMiEx%NY-1 , 0 : b%sggMiEx%NZ-1 )  , intent( IN)     ::  sggMiEx
      integer(kind = INTEGERSIZEOFMEDIAMATRICES), dimension ( 0 : b%sggMiEy%NX-1 , 0 : b%sggMiEy%NY-1 , 0 : b%sggMiEy%NZ-1 )  , intent( IN)     ::  sggMiEy
      integer(kind = INTEGERSIZEOFMEDIAMATRICES), dimension ( 0 : b%sggMiEz%NX-1 , 0 : b%sggMiEz%NY-1 , 0 : b%sggMiEz%NZ-1 )  , intent( IN)     ::  sggMiEz
      integer(kind = INTEGERSIZEOFMEDIAMATRICES), dimension ( 0 : b%sggMiHx%NX-1 , 0 : b%sggMiHx%NY-1 , 0 : b%sggMiHx%NZ-1 )  , intent( IN)     ::  sggMiHx
      integer(kind = INTEGERSIZEOFMEDIAMATRICES), dimension ( 0 : b%sggMiHy%NX-1 , 0 : b%sggMiHy%NY-1 , 0 : b%sggMiHy%NZ-1 )  , intent( IN)     ::  sggMiHy
      integer(kind = INTEGERSIZEOFMEDIAMATRICES), dimension ( 0 : b%sggMiHz%NX-1 , 0 : b%sggMiHz%NY-1 , 0 : b%sggMiHz%NZ-1 )  , intent( IN)     ::  sggMiHz
      !
      !!!
      character(len=*), INTENT(in) :: wiresflavor
      logical  ::  saveall,singlefilewrite,NF2FFDecim, INIT,GEOM,ASIGNA,electric,magnetic
      character (LEN=BUFSIZE)  ::  p1,p2
      real (kind=RKIND_tiempo) :: lastexecutedtime
      
      character (len=*), intent(in)  ::  nEntradaRoot

      integer (kind=4)  ::  i,field,ii,i1,j1,k1,n,i2,j2,k2,initialtimestep, finaltimestep,NO,NO2,iwi,iwj,compo,ntime,ntimeforvolumic,iff1,i0t

      logical, intent(inout)   ::  ThereAreObservation,ThereAreFarFields
      logical, intent(in)      ::  ThereAreWires,resume
      character (LEN=BUFSIZE)  ::  chari,charj,chark,chari2,charj2,chark2,charNO
      character (LEN=BUFSIZE)  ::  ext,extpoint,adum,prefix_field
      logical  ::  incident,errnofile,first
      REAL (KIND=RKIND)    ::  rdum,field1,field2
      REAL (KIND=RKIND_tiempo)    ::  at,dtevol,tiempo1,tiempo2
      integer (kind=4)  ::  unit,ndum,unitmaster,conta,III,JJJ,KKK,pozi,i1t,j1t,k1t
      character (LEN=BUFSIZE)  ::  whoami,whoamishort
      logical :: ok,existe,wrotemaster,found
      integer (kind=8)  :: memo,ntini,ntfin
      character(LEN=BUFSIZE) :: buff,path,buff2
#ifdef CompileWithMPI
      integer(kind=MPI_OFFSET_KIND) disp
      integer (kind=4)  ::  ierr
#endif
      logical :: Esborde
      integer (kind=4)  ::  imed,imed1,imed2,imed3,imed4
      integer (kind=4)  ::  thefile !for file management
!for dft
      REAL (KIND=RKIND), allocatable, dimension(:) :: signal,fqPos
      REAL (KIND=RKIND_tiempo), allocatable, dimension(:) :: samplingtime
      complex (kind=CKIND), allocatable, dimension(:) :: fqValues
      integer (kind=4) :: timesteps,klk,fqlength
      integer :: my_iostat
!
      
   if (size.gt.maxcores) then
       print *,'Maximum cores ',maxcores,' reached.  to recompile'
       stop
   endif
      eps0=eps00; mu0=mu00; !chapuz para convertir la variables de paso en globales
!
!!
      output=> null()
#ifdef CompileWithMPI
      valores=>null()
      newvalores => null()  !auxiliary for Bloque currents sync
#endif

      unitmaster=-1000 !!!no se bien. Lo pongo absurdo ����
      unit=1000 !initial
      if (unit >= 2.0_RKIND**31.0_RKIND-1.0_RKIND) then
         call stoponerror(layoutnumber,size,'Excesive number of probes')
      endif
      !
      write(whoamishort,'(i5)') layoutnumber+1
      write(whoami,'(a,i5,a,i5,a)') '(',layoutnumber+1,'/',size,') '

      call crea_gnuplot

      allocate (InvEps(0 : sgg%NumMedia),InvMu(0 : sgg%NumMedia))

      incident=.false.
      InvEps(0 : sgg%NumMedia)=1.0_RKIND/(Eps0*sgg%Med(0 : sgg%NumMedia)%Epr)
      InvMu (0 : sgg%NumMedia)=1.0_RKIND/(Mu0* sgg%Med(0 : sgg%NumMedia)%Mur) 


      allocate (output(1 : sgg%NumberRequest))
      output(1 : sgg%NumberRequest)%Trancos = -1
      output(1 : sgg%NumberRequest)%SaveAll = .false.
      output(1 : sgg%NumberRequest)%TimesWritten = -1

      !preprocesa par evitar overflows en integers
      !!!!!
      !corrige los time step en funcion del final (bug OLD 1 solo time step pero con trancos mayor que 0) !31/03/2014
      do ii=1,sgg%NumberRequest
         sgg%Observation(ii)%done = .false.
         sgg%Observation(ii)%begun = .false.
         sgg%Observation(ii)%flushed = .false.

         sgg%observation(ii)%TimeStep = max(sgg%observation(ii)%TimeStep,sgg%dt)
         if (10.0_RKIND*(sgg%observation(ii)%FinalTime-sgg%observation(ii)%InitialTime)/min(sgg%dt,sgg%observation(ii)%TimeStep ) >= huge (1_4)) then
            sgg%observation(ii)%FinalTime = sgg%observation(ii)%InitialTime+min(sgg%dt,sgg%observation(ii)%TimeStep )*huge(1_4)/10.0_RKIND
         endif

         if (sgg%observation(ii)%InitialTime < sgg%observation(ii)%TimeStep ) then
            sgg%observation(ii)%InitialTime = 0.0_RKIND !!! sgg%observation(ii)%TimeStep  !para que saque tambien el instante inicial
         endif
         if (( sgg%observation(ii)%TimeStep > sgg%observation(ii)%FinalTime   - sgg%observation(ii)%InitialTime)) then !!bug 260321 sondas Frequency Transfer .or.(sgg%observation(ii)%TimeStep ==0))  then
            ! sgg%observation(ii)%TimeStep = sgg%observation(ii)%FinalTime   - sgg%observation(ii)%InitialTime !!bug 260321 sondas Frequency Transfer !ya esta puesto antes bien
            if (sgg%observation(ii)%P(1)%what == mapvtk) then
                sgg%observation(ii)%FinalTime=0.
                sgg%observation(ii)%InitialTime=0.
            else                
              sgg%observation(ii)%FinalTime= sgg%observation(ii)%InitialTime + sgg%observation(ii)%TimeStep !al menos guarda uno !!bug 260321 sondas Frequency Transfer
            endif
            
         endif
         sgg%observation(ii)%InitialTime = int(sgg%observation(ii)%InitialTime/sgg%dt)*sgg%dt
         sgg%observation(ii)%FinalTime =   int(sgg%observation(ii)%FinalTime  /sgg%dt)*sgg%dt
         !!!!!sgg%observation(ii)%InitialTime = max(sgg%observation(ii)%InitialTime,sgg%dt) !no empezar en 0*sgg%dt
         sgg%observation(ii)%FreqStep = min(sgg%observation(ii)%FreqStep,2.0_RKIND / sgg%dt)
         if ((sgg%observation(ii)%FreqStep > sgg%observation(ii)%FinalFreq   - sgg%observation(ii)%InitialFreq).or.(sgg%observation(ii)%FreqStep ==0)) then
            sgg%observation(ii)%FreqStep = sgg%observation(ii)%FinalFreq   - sgg%observation(ii)%InitialFreq
            sgg%observation(ii)%FinalFreq= sgg%observation(ii)%InitialFreq + sgg%observation(ii)%FreqStep
         endif

         !if (sgg%Observation(ii)%nP/=0) then
            if (.not.SGG%Observation(ii)%Volumic) then  !pasate por los guevos los tiempos iniciales finales y step de las sondas excepto las volumicas porque pueden ser grandes para evitar lo que me paso con el SIVA!!! 051115   
                    sgg%observation(ii)%Saveall=sgg%observation(ii)%Saveall.or.saveall 
                    output(ii)%SaveAll=sgg%observation(ii)%Saveall     
    !!!!!!salvo la volumicas, el resto lo guardan todo           !!!repensado y cambiado para que lo respete 070916
               !!!     output(ii)%SaveAll=.true.
               !!!     sgg%observation(ii)%Saveall=.true.
            else              
                    output(ii)%SaveAll=.false.
                    sgg%observation(ii)%Saveall=.false.
            endif
#ifdef miguelConformalStandAlone
            output(ii)%SaveAll=.false.
#endif
            !!!        if (SGG%Observation(ii)%TimeDomain) then !CREO QUE NO SE PRECISA PQ LAS FREQ. DOMAIN PROBES TAMBIEN DEBEN PASAR POR AQUI
            if (sgg%observation(ii)%Saveall)  then
                output(ii)%Trancos       = 1
                sgg%observation(ii)%InitialTime=0.0_RKIND
                sgg%observation(ii)%FinalTime=sgg%tiempo(FINALTIMESTEP+2)
            else
                output(ii)%Trancos       = max(1,int(sgg%Observation(ii)%TimeStep/sgg%dt))
                sgg%observation(ii)%InitialTime=max(0.0_RKIND,sgg%observation(ii)%InitialTime)
                sgg%observation(ii)%FinalTime=min(sgg%tiempo(FINALTIMESTEP+2),sgg%observation(ii)%FinalTime) !CLIPEA
                if (sgg%observation(ii)%FinalTime < sgg%observation(ii)%InitialTime) sgg%observation(ii)%FinalTime = sgg%observation(ii)%InitialTime
            endif
         !endif
         if (sgg%Observation(ii)%nP/=0) then
            if (sgg%observation(ii)%P(1)%what == mapvtk) then    
                    output(ii)%SaveAll=.false.
                    sgg%observation(ii)%Saveall=.false.
            endif              
         endif
!!!!
      END DO



      do ii=1,sgg%NumberRequest
         allocate (output(ii)%item(1 : sgg%Observation(ii)%nP))
#ifdef CompileWithMPI
         do i=1,sgg%Observation(ii)%nP
            output(ii)%item(i)%ZIorig=sgg%Observation(ii)%P(i)%ZI
            output(ii)%item(i)%ZEorig=sgg%Observation(ii)%P(i)%ZE
         end do
#endif
         output(ii)%TimesWritten=0 !for volumic probes
      END DO

      !eliminate unnecesary observation points
      do ii=1,sgg%NumberRequest
         do i=1,sgg%Observation(ii)%nP
         !trancos
         output(ii)%item(i)%Xtrancos=sgg%Observation(ii)%P(i)%Xtrancos
         output(ii)%item(i)%Ytrancos=sgg%Observation(ii)%P(i)%Ytrancos
         output(ii)%item(i)%Ztrancos=sgg%Observation(ii)%P(i)%Ztrancos      
         output(ii)%item(i)%XItrancos = int(sgg%Observation(ii)%P(i)%XI/output(ii)%item(i)%Xtrancos)
         output(ii)%item(i)%XEtrancos = int(sgg%Observation(ii)%P(i)%XE/output(ii)%item(i)%Xtrancos)
         output(ii)%item(i)%YItrancos = int(sgg%observation(ii)%P(i)%YI/output(ii)%item(i)%Ytrancos)
         output(ii)%item(i)%YEtrancos = int(sgg%observation(ii)%P(i)%YE/output(ii)%item(i)%Ytrancos)
         output(ii)%item(i)%ZItrancos = int(sgg%observation(ii)%P(i)%ZI/output(ii)%item(i)%Ztrancos)
         output(ii)%item(i)%ZEtrancos = int(sgg%observation(ii)%P(i)%ZE/output(ii)%item(i)%Ztrancos)
         !
         if (mod(sgg%Observation(ii)%P(i)%XI,output(ii)%item(i)%Xtrancos) /= 0) output(ii)%item(i)%XItrancos=output(ii)%item(i)%XItrancos+1
         if (mod(sgg%Observation(ii)%P(i)%YI,output(ii)%item(i)%Ytrancos) /= 0) output(ii)%item(i)%YItrancos=output(ii)%item(i)%YItrancos+1
         if (mod(sgg%Observation(ii)%P(i)%ZI,output(ii)%item(i)%Ztrancos) /= 0) output(ii)%item(i)%ZItrancos=output(ii)%item(i)%ZItrancos+1
                        
         !fin trancos
#ifdef CompileWithMPI
            output(ii)%item(i)%MPISubComm=-1 !just to void it
#endif
            field=sgg%observation(ii)%P(i)%What
            select case (field)
             case (iBloqueJx,iBloqueJy,iBloqueMx,iBloqueMy)     !the end point info is only relevant for Bloque
               if ((sgg%Observation(ii)%P(i)%ZI >  sgg%Sweep(fieldo(field,'Z'))%ZE).or. &
               (sgg%Observation(ii)%P(i)%ZE <  sgg%Sweep(fieldo(field,'Z'))%ZI)) then
                  sgg%observation(ii)%P(i)%What=nothing
#ifdef CompileWithMPI
                  output(ii)%item(i)%MPISubComm=-1 !just to void it
               else
                  output(ii)%item(i)%MPISubComm=1
               endif
               output(ii)%item(i)%MPIRoot=0
               If ((sgg%observation(ii)%P(i)%ZI >= sgg%Sweep(fieldo(field,'Z'))%ZI).and. &
               (sgg%observation(ii)%P(i)%ZI <= sgg%Sweep(fieldo(field,'Z'))%ZE)) then
                  output(ii)%item(i)%MPIRoot=layoutnumber
               endif
               !all of them must call the init routine even if they do not sync
               call MPIinitSubcomm(layoutnumber,size, &
               output(ii)%item(i)%MPISubComm,output(ii)%item(i)%MPIRoot,output(ii)%item(i)%MPIGroupIndex)

#else
               endif
#endif
             case (iEx,iVx,iEy,iVy,iHz,iBloqueMz,iJx,iJy) !in case of MPI the flushing is only cared by one of the sharing layouts
               !in case of MPI the flushing is only cared by one of the sharing layouts
               !este es el unico caso en el que un punto es susceptible de ser escrito por dos layouts. Por eso se lo echo
               ! solo a uno de ellos: al de abajo (a menos que que sea el layout de mas arriba, en cuyo caso tiene que tratarlo el) !bug del itc2 con el pathx hasta el borde
               if (((sgg%Observation(ii)%P(i)%ZI >= sgg%Sweep(fieldo(field,'Z'))%ZE).and.(layoutnumber /= size-1)).or. &
               (sgg%observation(ii)%P(i)%ZI <  sgg%Sweep(fieldo(field,'Z'))%ZI)                             ) then
                  sgg%observation(ii)%P(i)%What=Nothing !do not observe anything
               endif
             case (iEz,iVz,iJz,iBloqueJz,iHx,iHy)
               if ((sgg%Observation(ii)%P(i)%ZI > sgg%Sweep(fieldo(field,'Z'))%ZE).or. &
               (sgg%Observation(ii)%P(i)%ZI < sgg%Sweep(fieldo(field,'Z'))%ZI)) then
                  sgg%observation(ii)%P(i)%What=Nothing !do not observe anything
               endif
               !the end point info is only relevant for Volumic and Bloque
             case (iExC,iEyC,iHzC,iMhC,iEzC,iHxC,iHyC,iMeC)
               if ((sgg%Observation(ii)%P(i)%ZI >  sgg%Sweep(fieldo(field,'Z'))%ZE).or. &
               (sgg%Observation(ii)%P(i)%ZE <  sgg%Sweep(fieldo(field,'Z'))%ZI)) then
                  sgg%Observation(ii)%P(i)%What=nothing
#ifdef CompileWithMPI
                  output(ii)%item(i)%MPISubComm=-1 !just to void it
               else
                  output(ii)%item(i)%MPISubComm=1
               endif
               output(ii)%item(i)%MPIRoot=0
               If ((sgg%observation(ii)%P(i)%ZI >= sgg%Sweep(fieldo(field,'Z'))%ZI).and. &
               (sgg%observation(ii)%P(i)%ZI <= sgg%Sweep(fieldo(field,'Z'))%ZE)) then
                  output(ii)%item(i)%MPIRoot=layoutnumber
               endif
               !all of them must call the init routine even if they do not sync
               call MPIinitSubcomm(layoutnumber,size, &
               output(ii)%item(i)%MPISubComm,output(ii)%item(i)%MPIRoot,output(ii)%item(i)%MPIGroupIndex)

#else
               endif
#endif
               !the end point info is only relevant for Volumic and Bloque
             case (iCur,iCurX,iCurY,iCurZ,mapvtk)
               if ((sgg%Observation(ii)%P(i)%ZI >=  sgg%Sweep(iHz)%ZE).or. &
               (sgg%Observation(ii)%P(i)%ZE <  sgg%Sweep(iHZ)%ZI)) then
                  sgg%Observation(ii)%P(i)%What=nothing
#ifdef CompileWithMPI
                  output(ii)%item(i)%MPISubComm=-1 !just to void it
               else
                  output(ii)%item(i)%MPISubComm=1
               endif
               !clipeo los finales porque luego tengo que interpolar y el MPI me puede molestar 06/07/15
               if ((field==icur).or.(field==icurX).or.(field==icurY).or.(field==mapvtk)) then
                  sgg%Observation(ii)%P(i)%ZE=Min(sgg%Observation(ii)%P(i)%ZE,sgg%Sweep(iHx)%ZE)
               endif
               !!!!!!!!!!!!
               output(ii)%item(i)%MPIRoot=0
               If ((sgg%observation(ii)%P(i)%ZI >= sgg%Sweep(fieldo(field,'Z'))%ZI).and. &
               (sgg%observation(ii)%P(i)%ZI <= sgg%Sweep(fieldo(field,'Z'))%ZE)) then
                  output(ii)%item(i)%MPIRoot=layoutnumber
               endif
               !all of them must call the init routine even if they do not sync
               call MPIinitSubcomm(layoutnumber,size, &
               output(ii)%item(i)%MPISubComm,output(ii)%item(i)%MPIRoot,output(ii)%item(i)%MPIGroupIndex)

#else
               endif
#endif
             case (FarField)
#ifdef CompileWithNF2FF
               if ( (sgg%observation(ii)%P(1)%ZI >  sgg%SINPMLSweep(IHz)%ZE).or. &  !MPI NO DUPLICAR CALCULOS
               (sgg%observation(ii)%P(1)%ZE < sgg%SINPMLSweep(iHz)%ZI)) then

                  sgg%Observation(ii)%P(i)%What=nothing
#ifdef CompileWithMPI
                  output(ii)%item(i)%MPISubComm=-1 !just to void it
               else
                  output(ii)%item(i)%MPISubComm=1
               endif
               output(ii)%item(i)%MPIRoot=0
               If ((sgg%observation(ii)%P(i)%ZI >= sgg%SINPMLSweep(iHz)%ZI).and. &
               (sgg%observation(ii)%P(i)%ZI < sgg%SINPMLSweep(iHz)%ZE)) then
                  output(ii)%item(i)%MPIRoot=layoutnumber
               endif
               !all of them must call the init routine even if they do not sync
               call MPIinitSubcomm(layoutnumber,size, &
               output(ii)%item(i)%MPISubComm,output(ii)%item(i)%MPIRoot,output(ii)%item(i)%MPIGroupIndex)
#else
               endif
#endif
               !
#else
               call stoponerror(layoutnumber,size,'Current version does not support Near to Far field. Recompile')
#endif
            end select
         end do
      end do

      !
      ThereAreObservation=.false.
      ThereAreFarFields=.false.
      do ii=1,sgg%NumberRequest
         do i=1,sgg%Observation(ii)%nP
            field=sgg%observation(ii)%P(i)%what
            if (field/=nothing) ThereAreObservation=.true.
         end do
      end do
      !
      memo=0
      !
      IF (ThereAreObservation) then
#ifdef CompileWithMPI
         allocate (valores(0 : BuffObse),newvalores(0 : BuffObse))
         valores = 0.0_RKIND
         newvalores = 0.0_RKIND
#endif
         if (sgg%NumPlaneWaves>=1) incident=.true. !150419 creo que no se ha dado nunca >1 porque los RC tocan el numero de modos pero solo hay una planewave
         !Write also the incident fields in case there are plane waves (useful in SE calculations)
         open (119,file=trim(adjustl(nEntradaRoot))//'_Outputrequests_'//trim(adjustl(whoamishort))//'.txt' )
         write (119,'(a)') '!END'
         close (119,status='delete')
         my_iostat=0
9138     if(my_iostat /= 0) write(*,fmt='(a)',advance='no') '.' !!if(my_iostat /= 0) print '(i5,a1,i4,2x,a)',9138,'.',layoutnumber,trim(adjustl(nEntradaRoot))//'_Outputrequests_'//trim(adjustl(whoamishort))//'.txt'
         open (19,file=trim(adjustl(nEntradaRoot))//'_Outputrequests_'//trim(adjustl(whoamishort))//'.txt',err=9138,iostat=my_iostat,status='new',action='write')

#ifdef CompileWithWires
         if ((trim(adjustl(wiresflavor))=='holland') .or. &
             (trim(adjustl(wiresflavor))=='transition')) then
            if (Therearewires) Hwireslocal => GetHwires()
         endif
#endif
#ifdef CompileWithBerengerWires
         if (trim(adjustl(wiresflavor))=='berenger') then
            if (Therearewires) Hwireslocal_Berenger => GetHwires_Berenger()
         endif
#endif
#ifdef CompileWithSlantedWires
         if ((trim(adjustl(wiresflavor))=='slanted').or.(trim(adjustl(wiresflavor))=='semistructured')) then
            if (Therearewires) Hwireslocal_Slanted => GetHwires_Slanted()
         endif
#endif


         !!!!!!!!!Comun a todas las sondas freqdomain
         do ii=1,sgg%NumberRequest
            if (SGG%Observation(ii)%FreqDomain) then
               !
               !Read the time evolution file common to all iten i
               !
               output(ii)%InitialFreq=sgg%observation(ii)%InitialFreq
               output(ii)%FinalFreq=sgg%observation(ii)%FinalFreq
               output(ii)%FreqStep=sgg%observation(ii)%FreqStep
               !
               if (SGG%Observation(ii)%FreqStep/=0) then
                  output(ii)%NumFreqs=int(abs(SGG%Observation(ii)%InitialFreq-SGG%Observation(ii)%FinalFreq)/SGG%Observation(ii)%FreqStep)+1
               else
                  output(ii)%NumFreqs=1 !default
               endif
               if ((output(ii)%NumFreqs < 0)) then
                  Buff='Freq. range for Freq. probes invalid'
                  call stoponerror(layoutnumber,size,Buff)
               endif
               if ((output(ii)%NumFreqs > 100000)) then
                  Buff='Too many Freqs requested (>100000)'
                  call stoponerror(layoutnumber,size,Buff)
               endif

               ALLOCATE ( output(ii)%Freq (1: output(ii)%NumFreqs), &
               output(ii)%auxExp_E(1: output(ii)%NumFreqs), &
               output(ii)%auxExp_H(1: output(ii)%NumFreqs))
!
               pozi=index(sgg%observation(ii)%outputrequest,'_log_')
               if (pozi == 0) then

                  do iff1=1,output(ii)%NumFreqs
                     output(ii)%Freq(iff1)=output(ii)%InitialFreq + (iff1-1) *output(ii)%FreqStep

                  END DO
               else !logaritmico
                  output(ii)%InitialFreq=log10(output(ii)%InitialFreq)
                  output(ii)%FinalFreq=log10(output(ii)%FinalFreq)
                  output(ii)%FreqStep=abs((output(ii)%InitialFreq-output(ii)%FinalFreq)/(output(ii)%NumFreqs))

                  do iff1=1,output(ii)%NumFreqs
                     output(ii)%Freq(iff1)=10.0_RKIND **(output(ii)%InitialFreq + (iff1-1) *output(ii)%FreqStep)

                  END DO
               endif
               !
               errnofile = .FALSE.

               IF (SGG%Observation(ii)%Transfer) then
                  ALLOCATE (output(ii)%dftEntrada(1:output(ii)%NumFreqs))
                  output(ii)%dftEntrada=0.0_RKIND
!
                  INQUIRE (file=trim(adjustl(sgg%observation(ii)%FileNormalize)), EXIST=errnofile)
                  IF ( .NOT. errnofile) THEN
                     Buff=trim(adjustl(sgg%observation(ii)%FileNormalize))//' NORMALIZATION FILE DOES NOT EXIST'
                     CALL STOPONERROR (layoutnumber,size,Buff)
                  END IF
!precuenta
                  timesteps=0
                  OPEN (15, file=trim(adjustl(sgg%observation(ii)%FileNormalize)))
                  do
                    READ (15,*, end=998) tiempo1, field1
                    timesteps=timesteps+1
                    continue
                  end do
998               CLOSE (15)
                  allocate (samplingTime(1:timesteps))
                  allocate (signal(1:timesteps))
                  signal=0.0_RKIND
                  samplingTime=0.0_RKIND_tiempo
!
                  !read the normalization file and find its DFT
                  OPEN (15, file=trim(adjustl(sgg%observation(ii)%FileNormalize)))
                  DO klk=1,timesteps
                     READ (15,*) samplingTime(klk),signal(klk)
                  END DO
                  CLOSE (15)
!
                  !niapa quitar 200120 ojooo
                  if (niapapostprocess) then
                      print *,'Correcting in observation ',timesteps,trim(adjustl(sgg%observation(ii)%FileNormalize))
                      DO klk=1,timesteps
                          samplingTime(klk)=real(klk*sgg%dt,RKIND_tiempo)
                      end do
                  endif
                  !fin niapa
                  
                  fqLength=output(ii)%NumFreqs
                  allocate (fqValues(1:fqLength),fqPos(1:fqLength))
                  fqValues(1:fqLength)=0.0_RKIND
                  fqPos(1:fqLength)=output(ii)%Freq(1:fqLength)
                  call dtft(fqValues, fqPos, fqLength, samplingTime, signal, timesteps)
                  output(ii)%dftEntrada=fqValues
                  deallocate (samplingTime,signal,fqValues,fqPos)
               ENDIF !DEL TRANSFER
            endif !del freqdomain
         end do
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!



         do ii=1,sgg%NumberRequest
            wrotemaster=.false.
            loop_ob :  do i=1,sgg%Observation(ii)%nP
               I1=sgg%observation(ii)%P(i)%XI
               J1=sgg%observation(ii)%P(i)%YI
               K1=sgg%observation(ii)%P(i)%ZI
               I2=sgg%observation(ii)%P(i)%XE
               J2=sgg%observation(ii)%P(i)%YE
               K2=sgg%observation(ii)%P(i)%ZE
               NO=sgg%observation(ii)%P(i)%NODE
               write(chari,'(i7)') i1
               write(charj,'(i7)') j1
               write(chark,'(i7)') k1
               field=sgg%observation(ii)%P(i)%what
               select case(field)
                  !
               case (iEx,iEy,iEz,iVx,iVy,iVz,iJx,iJy,iJz,iHx,iHy,iHz)
                  !
                  if (((field == iEx).or.(field == iEy).or.(field == iEz).or. &
                  (field == iHx).or.(field == iHy).or.(field == iHz)).and. &
                  (sgg%NumPlaneWaves>=1)) then
                     output(ii)%item(i)%columnas=2
                  elseif ((field == iJx).or.(field == iJy).or.(field == iJz)) then
                     output(ii)%item(i)%columnas=5 ! corriente -e*dl vplus vminus vplus-vminus
                  elseif ((field == iVx).or.(field == iVy).or.(field == iVz)) then
                     output(ii)%item(i)%columnas=1
                  else
                     output(ii)%item(i)%columnas=1
                  endif                
                  !mpdir 190319       !desrotacion para que los nombres sean correctos 
                  if      (mpidir==3) then
                      extpoint=trim(adjustl(chari)) //'_'//trim(adjustl(charj)) //'_'//trim(adjustl(chark))
                      select case (field)
                      case (iEx)
                        prefix_field=prefix(iEx)
                      case (iEy)
                        prefix_field=prefix(iEy)
                      case (iEz)
                        prefix_field=prefix(iEz)
                      case (iJx)
                        prefix_field=prefix(iJx)
                      case (iJy)
                        prefix_field=prefix(iJy)
                      case (iJz)
                        prefix_field=prefix(iJz)
                      case (iVx)
                        prefix_field=prefix(iVx)
                      case (iVy)
                        prefix_field=prefix(iVy)
                      case (iVz)
                        prefix_field=prefix(iVz)
                      case (iHx)
                        prefix_field=prefix(iHx)
                      case (iHy)
                        prefix_field=prefix(iHy)
                      case (iHz)
                        prefix_field=prefix(iHz)
                      case default
                        prefix_field=prefix(field)
                      end select
                  elseif  (mpidir==2) then
                      extpoint=trim(adjustl(charj)) //'_'//trim(adjustl(chark)) //'_'//trim(adjustl(chari))
                      select case (field)
                      case (iEx)
                        prefix_field=prefix(iEz)
                      case (iEy)
                        prefix_field=prefix(iEx)
                      case (iEz)
                        prefix_field=prefix(iEy)
                      case (iJx)
                        prefix_field=prefix(iJz)
                      case (iJy)
                        prefix_field=prefix(iJx)
                      case (iJz)
                        prefix_field=prefix(iJy)
                      case (iVx)
                        prefix_field=prefix(iVz)
                      case (iVy)
                        prefix_field=prefix(iVx)
                      case (iVz)
                        prefix_field=prefix(iVy)
                      case (iHx)
                        prefix_field=prefix(iHz)
                      case (iHy)
                        prefix_field=prefix(iHx)
                      case (iHz)
                        prefix_field=prefix(iHy)
                      case default
                        prefix_field=prefix(field)
                      end select
                  elseif  (mpidir==1) then
                      extpoint=trim(adjustl(chark)) //'_'//trim(adjustl(chari)) //'_'//trim(adjustl(charj))
                      select case (field)
                      case (iEx)
                        prefix_field=prefix(iEy)
                      case (iEy)
                        prefix_field=prefix(iEz)
                      case (iEz)
                        prefix_field=prefix(iEx)
                      case (iJx)
                        prefix_field=prefix(iJy)
                      case (iJy)
                        prefix_field=prefix(iJz)
                      case (iJz)
                        prefix_field=prefix(iJx)
                      case (iVx)
                        prefix_field=prefix(iVy)
                      case (iVy)
                        prefix_field=prefix(iVz)
                      case (iVz)
                        prefix_field=prefix(iVx)
                      case (iHx)
                        prefix_field=prefix(iHy)
                      case (iHy)
                        prefix_field=prefix(iHz)
                      case (iHz)
                        prefix_field=prefix(iHx)
                      case default
                        prefix_field=prefix(field)
                      end select
                  else
                      call stoponerror(layoutnumber,size,'Buggy error in mpidir. ')
                  endif
                  !     
                  if ((field == iJx).or.(field == iJy).or.(field == iJz)) then
                     write(charNO,'(i7)') NO
                     !append the number of the segment
                     extpoint=trim(adjustl(extpoint))//'_s'//trim(adjustl(charNO))
                  endif
                  ext=trim(adjustl(nEntradaRoot))//'_'//trim(adjustl(sgg%observation(ii)%outputrequest))
                  !do not use layername since no two observations from different layers will overlap
                  output(ii)%item(i)%path= &
                  trim(adjustl(ext))//'_'//trim(adjustl(prefix_field))//trim(adjustl(extpoint))//'.dat'
                  !
                  unit=unit+1
                  if (unit >= 2.0_RKIND  **31.0_RKIND-1.0_RKIND) then
                     call stoponerror(layoutnumber,size,'Excesive number of probes')
                  endif
                  output(ii)%item(i)%unit=unit
                  !

                  !!!busca nombres de ficheros por duplicado y resuelve la duplicidad
                  call checkduplicatenames
                  !!!!!!
                  
                  my_iostat=0
9235              if(my_iostat /= 0) write(*,fmt='(a)',advance='no'), '.' !!if(my_iostat /= 0) print '(i5,a1,i4,2x,a)',9235,layoutnumber,trim(adjustl(nEntradaRoot))//'_Outputrequests_'//trim(adjustl(whoamishort))//'.txt'
                  write(19,'(a)',err=9235,iostat=my_iostat) trim(adjustl(output(ii)%item(i)%path))
                  !
                  memo=memo+rkind*BuffObse
                  if (memo > MaxMemoryProbes) then
                     call stoponerror(layoutnumber,size,'Recompile: excesive memory for probes.'// &
                     &                                   'Increase MaxMemoryProbes')
                  endif
                  allocate (output(ii)%item(i)%valor(0 : BuffObse))
                  output(ii)%item(i)%valor(0 : BuffObse)=0.0_RKIND
#ifdef CompileWithWires
                  if ((trim(adjustl(wiresflavor))=='holland') .or. &
                      (trim(adjustl(wiresflavor))=='transition')) then
                     found=.false.
                     if ((Therearewires).and.((field==iJx).or.(field==iJy).or.(field==iJz))) then

                        memo=memo+3*4*BuffObse
                        if (memo > MaxMemoryProbes) then
                           call stoponerror(layoutnumber,size,'Recompile: excesive memory for probes.'// &
                           &                                   'Increase MaxMemoryProbes')
                        endif
                        allocate (&
                        output(ii)%item(i)%valor2(0 : BuffObse), &
                        output(ii)%item(i)%valor3(0 : BuffObse), &
                        output(ii)%item(i)%valor4(0 : BuffObse), &
                        output(ii)%item(i)%valor5(0 : BuffObse))
                        output(ii)%item(i)%valor2(0 : BuffObse)=0.0_RKIND
                        output(ii)%item(i)%valor3(0 : BuffObse)=0.0_RKIND
                        output(ii)%item(i)%valor4(0 : BuffObse)=0.0_RKIND
                        output(ii)%item(i)%valor5(0 : BuffObse)=0.0_RKIND
                        output(ii)%item(i)%valorsigno=+1
                        !en caso de hilos se necesitan
                        !parsea los hilos
                        found=.false.
                        output(ii)%item(i)%segmento => HWireslocal%NullSegment !por si no encuentra el segmento por estar repetido !bug serio que impide discernir para segmentos
                        do n=1,HWireslocal%NumCurrentSegments
                           if ((HWireslocal%CurrentSegment(n)%origindex==no).and. & !el nodo tambien debe coincidir !bug serio que impide discernir para segmentos paralelos 12/09/13 cazado gracias a Model_unidos.nfde
                           (HWireslocal%CurrentSegment(n)%i==i1).and. &
                           (HWireslocal%CurrentSegment(n)%j==j1).and. &
                           (HWireslocal%CurrentSegment(n)%k==k1).and. &
                           (HWireslocal%CurrentSegment(n)%tipofield*10==field))  then
                              !I have chosen IJx=10 IEx, etc. !do not change
                              output(ii)%item(i)%segmento => HWireslocal%CurrentSegment(n)
                              if (output(ii)%item(i)%segmento%orientadoalreves) output(ii)%item(i)%valorsigno=-1
                              found=.true.
                           endif
                        end do
                        if (.not.found) then
                           !busca por si fuera un multirabito
                           buscarabono: do iwi=1,Hwireslocal%NumDifferentWires
                              do iwj=1,sgg%Med(Hwireslocal%WireTipoMedio(iwi))%wire(1)%numsegmentos
                                 if ((no==sgg%Med(Hwireslocal%WireTipoMedio(iwi))%wire(1)%segm(iwj)%origindex).and.sgg%Med(Hwireslocal%WireTipoMedio(iwi))%wire(1)%segm(iwj)%multirabo) then
                                    no2=sgg%Med(Hwireslocal%WireTipoMedio(iwi))%wire(1)%segm(iwj)%multiraboDE
                                    do n=1,HWireslocal%NumCurrentSegments
                                       if (HWireslocal%CurrentSegment(n)%origindex==no2) then !el nodo tambien debe coincidir aunque no necesariamente coordenadas ni campo porque se ha cortado el rabo original
                                          !I have chosen IJx=10 IEx, etc. !do not change
                                          output(ii)%item(i)%segmento => HWireslocal%CurrentSegment(n)
                                          if (output(ii)%item(i)%segmento%orientadoalreves) output(ii)%item(i)%valorsigno=-1
                                          found=.true.
                                       endif
                                    end do
                                    found=.true.
                                    exit buscarabono
                                 endif
                              enddo
                           enddo buscarabono
                        endif
                     endif
                     if ((.not.found).and.((((field==iJx).or.(field==iJy).or.(field==iJz))))) then
                        sgg%Observation(ii)%P(i)%What=nothing
                          !ojoo 010423 para debugeo lbb1
                          write(buff,'(a,4i7,a)') 'ERROR: WIRE probe ',no,i1,j1,k1,' DOES NOT EXIST'
                          CALL WarnErrReport (buff,.true.)
                     endif
                  endif
#endif
#ifdef CompileWithBerengerWires
                  if (trim(adjustl(wiresflavor))=='berenger') then 
                     found=.false.
                     if ((Therearewires).and.((field==iJx).or.(field==iJy).or.(field==iJz))) then

                        memo=memo+3*4*BuffObse
                        if (memo > MaxMemoryProbes) then
                           call stoponerror(layoutnumber,size,'Recompile: excesive memory for probes.'// &
                           &                                   'Increase MaxMemoryProbes')
                        endif
                        allocate (&
                        output(ii)%item(i)%valor2(0 : BuffObse), &
                        output(ii)%item(i)%valor3(0 : BuffObse), &
                        output(ii)%item(i)%valor4(0 : BuffObse), &
                        output(ii)%item(i)%valor5(0 : BuffObse))
                        output(ii)%item(i)%valor2(0 : BuffObse)=0.0_RKIND
                        output(ii)%item(i)%valor3(0 : BuffObse)=0.0_RKIND
                        output(ii)%item(i)%valor4(0 : BuffObse)=0.0_RKIND
                        output(ii)%item(i)%valor5(0 : BuffObse)=0.0_RKIND
                        output(ii)%item(i)%valorsigno=+1
                        !en caso de hilos se necesitan
                        !parsea los hilos
                        found=.false.
                        do n=1,Hwireslocal_Berenger%NumSegments
                           if (Hwireslocal_Berenger%Segments(n)%IndexSegment==no) then !solo miro el nodo. dama ya corrige esto la observation de Berenger 
                           !bug dectectado por OLD 311019 con probe_issue.nfde
!                           if ((Hwireslocal_Berenger%Segments(n)%IndexSegment==no).and. & !el nodo tambien debe coincidir !bug serio que impide discernir para segmentos paralelos 12/09/13 cazado gracias a Model_unidos.nfde
!                           (Hwireslocal_Berenger%Segments(n)%ii==i1).and. &
!                           (Hwireslocal_Berenger%Segments(n)%ji==j1).and. &
!                           (Hwireslocal_Berenger%Segments(n)%ki==k1).and. &
!                           (Hwireslocal_Berenger%Segments(n)%orient*10==field))  then
                              !I have chosen IJx=10 IEx, etc. !do not change
                              output(ii)%item(i)%segmento_Berenger => Hwireslocal_Berenger%Segments(n)
                              if (output(ii)%item(i)%segmento_Berenger%orientadoalreves) output(ii)%item(i)%valorsigno=-1
                              found=.true.
                           endif
                        end do
                     endif
                     if ((.not.found).and.((((field==iJx).or.(field==iJy).or.(field==iJz))))) then
                        sgg%Observation(ii)%P(i)%What=nothing
                        write(buff,'(a,4i7,a)') 'ERROR: WIRE probe ',no,i1,j1,k1,' DOES NOT EXIST'
                        CALL WarnErrReport (buff,.TRUE.)
                     endif
                  endif
#endif
#ifdef CompileWithSlantedWires
                  if ((trim(adjustl(wiresflavor))=='slanted').or.(trim(adjustl(wiresflavor))=='semistructured')) then 
                     found=.false.
                     if ((Therearewires).and.((field==iJx).or.(field==iJy).or.(field==iJz))) then

                        memo=memo+3*4*BuffObse
                        if (memo > MaxMemoryProbes) then
                           call stoponerror(layoutnumber,size,'Recompile: excesive memory for probes.'// &
                           &                                   'Increase MaxMemoryProbes')
                        endif
                        allocate (&
                        output(ii)%item(i)%valor2(0 : BuffObse), &
                        output(ii)%item(i)%valor3(0 : BuffObse), &
                        output(ii)%item(i)%valor4(0 : BuffObse), &
                        output(ii)%item(i)%valor5(0 : BuffObse))
                        output(ii)%item(i)%valor2(0 : BuffObse)=0.0_RKIND
                        output(ii)%item(i)%valor3(0 : BuffObse)=0.0_RKIND
                        output(ii)%item(i)%valor4(0 : BuffObse)=0.0_RKIND
                        output(ii)%item(i)%valor5(0 : BuffObse)=0.0_RKIND
                        output(ii)%item(i)%valorsigno=+1
                        !en caso de hilos se necesitan
                        !parsea los hilos
                        found=.false.
                        do n=1,Hwireslocal_Slanted%NumSegments
                           if (Hwireslocal_Slanted%Segments(n)%ptr%Index==no)  then
                              !I have chosen IJx=10 IEx, etc. !do not change
                              output(ii)%item(i)%segmento_Slanted => Hwireslocal_Slanted%Segments(n)%ptr
                              found=.true.
                           endif
                        end do
                     endif 
                     !010423  creo que si no lo encuentra es porque el indice es el exterior bug lbb1 epg 0323
                     if (.not.found) then   
                         do n=1,Hwireslocal_Slanted%NumSegments
                           if (Hwireslocal_Slanted%Segments(n)%ptr%elotroindice ==no)  then
                              output(ii)%item(i)%segmento_Slanted => Hwireslocal_Slanted%Segments(n)%ptr
                              found=.true.
                           endif
                         end do
                     endif                       
                     if ((.not.found).and.((((field==iJx).or.(field==iJy).or.(field==iJz))))) then
                        sgg%Observation(ii)%P(i)%What=nothing
                        
                          !ojoo 010423 para debugeo lbb1
                          write(buff,'(a,4i7,a)') 'ERROR: WIRE probe ',no,i1,j1,k1,' DOES NOT EXIST'
                          CALL WarnErrReport (buff,.true.)
                     endif
                  endif   
#endif            
                
                  !!!!!!!!!!!!!!!
                  !erase pre-existing data unless this is a resuming simulation
                  if (.not.resume) then
                     !
                     if (singlefilewrite) then
                        if (.not.wrotemaster) then
                           wrotemaster=.true.
                           unitmaster=output(ii)%item(i)%unit
                           output(ii)%item(i)%unitmaster=unitmaster
                           open (unitmaster,recl=1000,file= trim(adjustl(output(ii)%item(i)%path))//'_'//trim(adjustl(whoamishort))//'_master.bin')
                           write(unitmaster,*) '!END'
                           close(unitmaster,status='delete')
                           
                           my_iostat=0 
9238                       if(my_iostat /= 0) write(*,fmt='(a)',advance='no'), '.' !!if(my_iostat /= 0) print '(i5,a1,i4,2x,a)',9238,'.',layoutnumber,trim(adjustl(output(ii)%item(i)%path))//'_'//trim(adjustl(whoamishort))//'_master.bin'  
                           open (unitmaster,file= trim(adjustl(output(ii)%item(i)%path))//'_'//trim(adjustl(whoamishort))//'_master.bin',form='unformatted',err=9238,iostat=my_iostat,status='new',action='write')
                        else
                           output(ii)%item(i)%unitmaster=unitmaster
                        endif
                     else
                        !
                        open (output(ii)%item(i)%unit,recl=1000,file= trim(adjustl(output(ii)%item(i)%path)))
                        write(output(ii)%item(i)%unit,*) '!END'
                        close (output(ii)%item(i)%unit,status='delete')
                        my_iostat=0
8766                    if(my_iostat /= 0) write(*,fmt='(a)',advance='no'), '.' !!if(my_iostat /= 0) print '(i5,a1,i4,2x,a)',8766,'.',layoutnumber,trim(adjustl(output(ii)%item(i)%path)) 
                        open (output(ii)%item(i)%unit,recl=1000,file= trim(adjustl(output(ii)%item(i)%path)),err=8766,iostat=my_iostat,status='new',action='write')
                        write(output(ii)%item(i)%unit,'(a)') trim(adjustl( ' t'//'              '// &
                        trim(adjustl(output(ii)%item(i)%path))//'       '//trim(adjustl(suffix(field,incident))) ))
                     endif
                     !
                     !                        write(output(ii)%item(i)%unit,'(5a)') ' t','              ', &
                     !!                            trim(adjustl(prefix(field)))//trim(adjustl(extpoint)),'   ',trim(adjustl(suffix(field,incident)))
                     !                            trim(adjustl(output(ii)%item(i)%path)),'       ',trim(adjustl(suffix(field,incident)))
                     !
                  else !wipe out duplicate data after non synchronous data and field resuming
                     !
                     if (singlefilewrite) then
                        if (.not.wrotemaster) then
                           wrotemaster=.true.
                           unitmaster=output(ii)%item(i)%unit
                           output(ii)%item(i)%unitmaster=unitmaster
                           open (unitmaster,recl=1000,file= trim(adjustl(output(ii)%item(i)%path))//'_master'//'_'//trim(adjustl(whoamishort))//'_master.bin')
                           write(unitmaster,*) '!END'
                           close(unitmaster,status='delete')
                           my_iostat=0
9239                       if(my_iostat /= 0) write(*,fmt='(a)',advance='no'), '.' !!if(my_iostat /= 0) print '(i5,a1,i4,2x,a)',9239,'.',layoutnumber,trim(adjustl(output(ii)%item(i)%path))//'_'//trim(adjustl(whoamishort))//'_master.bin' 
                           open (unitmaster,file= trim(adjustl(output(ii)%item(i)%path))//'_'//trim(adjustl(whoamishort))//'_master.bin',form='unformatted',err=9239,iostat=my_iostat,status='new',action='write')
                        else
                           output(ii)%item(i)%unitmaster=unitmaster
                        endif
                     else
                        inquire(file= trim(adjustl(output(ii)%item(i)%path)),exist=existe)
                        if (.not.existe) then
                           call stoponerror(layoutnumber,size,'Data files for resuming non existent (Ex, etc.) '//trim(adjustl(output(ii)%item(i)%path)))
                        endif
                        !
                        open (output(ii)%item(i)%unit,recl=1000,access='sequential',file= trim(adjustl(output(ii)%item(i)%path)))
                        read(output(ii)%item(i)%unit,'(a)') adum !first line contains characters
                        cutting :  do
                           read(output(ii)%item(i)%unit,*,end=678) at 
                           if (rdum > lastexecutedtime) then
                               print '(i4,a,a,2e19.9e3)',quienmpi,'Cutting 1 ',trim(adjustl(output(ii)%item(i)%path)),at,lastexecutedtime
                               backspace(output(ii)%item(i)%unit)
                               endfile(output(ii)%item(i)%unit)
                               exit cutting
                           endif
                        end do cutting
678                     continue
                        close (output(ii)%item(i)%unit)
                        open (output(ii)%item(i)%unit,recl=1000,file= trim(adjustl(output(ii)%item(i)%path)),position='append')
                     endif

                  endif
                case (iBloqueJx,iBloqueJy,iBloqueJz,iBloqueMx,iBloqueMy,iBloqueMz)
                  output(ii)%item(i)%columnas=1
                  write(chari2,'(i7)') i2
                  write(charj2,'(i7)') j2
                  write(chark2,'(i7)') k2
                  !mpidir 190319     !desrotacion para que los nombres sean correctos 
                  if      (mpidir==3) then
                      extpoint=trim(adjustl(chari)) //'_'//trim(adjustl(charj)) //'_'//trim(adjustl(chark))//'__'// &
                               trim(adjustl(chari2))//'_'//trim(adjustl(charj2))//'_'//trim(adjustl(chark2))
                      select case (field)
                      case (iBloqueJx)
                        prefix_field=prefix(iBloqueJx)
                      case (iBloqueJy)
                        prefix_field=prefix(iBloqueJy)
                      case (iBloqueJz)
                        prefix_field=prefix(iBloqueJz)
                      case (iBloqueMx)
                        prefix_field=prefix(iBloqueMx)
                      case (iBloqueMy)
                        prefix_field=prefix(iBloqueMy)
                      case (iBloqueMz)
                        prefix_field=prefix(iBloqueMz)
                      end select
                  elseif  (mpidir==2) then
                      extpoint=trim(adjustl(charj)) //'_'//trim(adjustl(chark)) //'_'//trim(adjustl(chari))//'__'// &
                               trim(adjustl(charj2))//'_'//trim(adjustl(chark2))//'_'//trim(adjustl(chari2))
                      select case (field)
                      case (iBloqueJx)
                        prefix_field=prefix(iBloqueJz)
                      case (iBloqueJy)
                        prefix_field=prefix(iBloqueJx)
                      case (iBloqueJz)
                        prefix_field=prefix(iBloqueJy)
                      case (iBloqueMx)
                        prefix_field=prefix(iBloqueMz)
                      case (iBloqueMy)
                        prefix_field=prefix(iBloqueMx)
                      case (iBloqueMz)
                        prefix_field=prefix(iBloqueMy)
                      end select
                  elseif  (mpidir==1) then
                      extpoint=trim(adjustl(chark)) //'_'//trim(adjustl(chari)) //'_'//trim(adjustl(charj))//'__'// &
                               trim(adjustl(chark2))//'_'//trim(adjustl(chari2))//'_'//trim(adjustl(charj2))
                      select case (field)
                      case (iBloqueJx)
                        prefix_field=prefix(iBloqueJy)
                      case (iBloqueJy)
                        prefix_field=prefix(iBloqueJz)
                      case (iBloqueJz)
                        prefix_field=prefix(iBloqueJx)
                      case (iBloqueMx)
                        prefix_field=prefix(iBloqueMy)
                      case (iBloqueMy)
                        prefix_field=prefix(iBloqueMz)
                      case (iBloqueMz)
                        prefix_field=prefix(iBloqueMx)
                      end select
                  else
                      call stoponerror(layoutnumber,size,'Buggy error in mpidir. ')
                  endif
                  !
                  ext=trim(adjustl(nEntradaRoot))//'_'//trim(adjustl(sgg%observation(ii)%outputrequest))
                  !
                  output(ii)%item(i)%path= &
                  trim(adjustl(ext))//'_'//trim(adjustl(prefix_field))//trim(adjustl(extpoint))//'.dat'

                  !
                  unit=unit+1
                  if (unit >= 2.0_RKIND  **31.0_RKIND-1.0_RKIND) then
                     call stoponerror(layoutnumber,size,'Excesive number of probes')
                  endif
                  output(ii)%item(i)%unit=unit
                  !
                  !!!busca nombres de ficheros por duplicado y resuelve la duplicidad
                  call checkduplicatenames
                  !!!!!!

                  memo=memo+rkind*BuffObse
                  if (memo > MaxMemoryProbes) then
                     call stoponerror(layoutnumber,size,'ERROR: Recompile: excesive memory for the probes.'// &
                     &                                   'Recompile increasing MaxMemoryProbes')
                  endif
                  allocate (output(ii)%item(i)%valor(0 : BuffObse))
                  output(ii)%item(i)%valor(0 : BuffObse)=0.0_RKIND
                  !readjust correctly the calculation region (for currents crossing the MPI region)
                  select case(field)
                   case (iBloqueJx,iBloqueJy)
                     sgg%observation(ii)%P(i)%ZI=max(sgg%Sweep(fieldo(field,'Z'))%ZI+1,sgg%observation(ii)%P(i)%ZI)
                     sgg%observation(ii)%P(i)%ZE=min(sgg%Sweep(fieldo(field,'Z'))%ZE  ,sgg%observation(ii)%P(i)%ZE)
                   case (iBloqueMx,iBloqueMy,iBloqueJz,iBloqueMz)
                     sgg%observation(ii)%P(i)%ZI=max(sgg%Sweep(fieldo(field,'Z'))%ZI  ,sgg%observation(ii)%P(i)%ZI)
                     sgg%observation(ii)%P(i)%ZE=min(sgg%Sweep(fieldo(field,'Z'))%ZE  ,sgg%observation(ii)%P(i)%ZE)
                  end select
                  !
#ifdef CompileWithMPI
                  if ((layoutnumber == output(ii)%item(i)%MPIRoot).or. &
                  (field==iBloqueJz).or.(field==iBloqueMz)) then !only the master
#endif
                  my_iostat=0
9837              if(my_iostat /= 0) write(*,fmt='(a)',advance='no'), '.' !!if(my_iostat /= 0) print '(i5,a1,i4,2x,a)',9837,layoutnumber,trim(adjustl(nEntradaRoot))//'_Outputrequests_'//trim(adjustl(whoamishort))//'.txt'
                  write(19,'(a)',err=9837,iostat=my_iostat) trim(adjustl(output(ii)%item(i)%path))
                  !erase pre-existing data unless this is a resuming simulation
                  if (.not.resume) then             
                     open (output(ii)%item(i)%unit,recl=1000,file= trim(adjustl(output(ii)%item(i)%path)))
                     write(output(ii)%item(i)%unit,*) '!END'
                     close (output(ii)%item(i)%unit,status='DELETE')
                     my_iostat=0
9838                 if(my_iostat /= 0) write(*,fmt='(a)',advance='no'), '.' !!if(my_iostat /= 0) print '(i5,a1,i4,2x,a)',9838,'.',layoutnumber,trim(adjustl(output(ii)%item(i)%path)) 
                     open (output(ii)%item(i)%unit,recl=1000,file= trim(adjustl(output(ii)%item(i)%path)),err=9838,iostat=my_iostat,status='new',action='write')
                     write(output(ii)%item(i)%unit,'(a)') trim(adjustl(' t'//'              '// &
                     trim(adjustl(output(ii)%item(i)%path))//'       '//&
                     trim(adjustl(suffix(field,incident))) ))
                                        !ojo a esto no le he corregido lo del mpidir de a 190319 porque estaba comentado
                     !
                     !                            write(output(ii)%item(i)%unit,'(5a)') ' t','              ', &
                     !                               !trim(adjustl(prefix(field)))//trim(adjustl(extpoint)),'   ',&
                     !                                trim(adjustl(output(ii)%item(i)%path)),'       ',&
                     !                                trim(adjustl(suffix(field,incident)))
                     !
                     !wipe out duplicate data after non synchronous data and field resuming
                  else
                     inquire(file= trim(adjustl(output(ii)%item(i)%path)),exist=existe)
                     if (.not.existe) then
                        call stoponerror(layoutnumber,size,'Data files for resuming non existent (Bloque, etc.) '//trim(adjustl(output(ii)%item(i)%path)))
                     endif
                     open (output(ii)%item(i)%unit,recl=1000,access='sequential',  &
                     file= trim(adjustl(output(ii)%item(i)%path)))
                     read(output(ii)%item(i)%unit,'(a)') adum !first line contains characters
                     cutting2 :  do
                        read(output(ii)%item(i)%unit,*,end=679) at 
                        if (at > lastexecutedtime) then
                               print '(i4,a,a,2e19.9e3)',quienmpi,'Cutting 2 ',trim(adjustl(output(ii)%item(i)%path)),at,lastexecutedtime
                               backspace(output(ii)%item(i)%unit)
                               endfile(output(ii)%item(i)%unit)
                           exit cutting2
                        endif
                     end do cutting2
679                  continue
                     close (output(ii)%item(i)%unit)
                     open (output(ii)%item(i)%unit,recl=1000,file= &
                     trim(adjustl(output(ii)%item(i)%path)),position='append')
                  endif
#ifdef CompileWithMPI
               endif
#endif
                  !voluminc Bloque current probes along edges (wires, surfaces
                case (iCur,iCurX,iCurY,iCurZ,mapvtk)
                  if (sgg%Observation(ii)%Volumic) then !they are necssaryly
                     if (sgg%Observation(ii)%nP /= 1) then
                        call stoponerror(layoutnumber,size,'ERROR! More than a volumic probe per group')
                     endif
                     !readjut correctly the calculation region
                     !de momento sere conservador 20/2/14 por lo que truene el MPI luego quitare el -1 si acaso !!!!����� !a priori puedo necesitar el HZ(alloc+1) para calcular las Bloque currents pero de momento me estoy quieto
                     sgg%observation(ii)%P(i)%ZI=max(sgg%Sweep(iEx)%ZI  ,sgg%observation(ii)%P(i)%ZI) !ojo estaba sweep(iEz) para ser conservador...puede dar problemas!! 03/07/15
                     sgg%observation(ii)%P(i)%ZE=min(sgg%Sweep(iEx)%ZE  ,sgg%observation(ii)%P(i)%ZE) !ojo estaba sweep(iEz) para ser conservador...puede dar problemas!! 03/07/15
                     !solo acepto que P(1:1) !!!
                     write(chari ,'(i7)') sgg%observation(ii)%P(i)%XI
                     write(charj ,'(i7)') sgg%observation(ii)%P(i)%YI
                     write(chark ,'(i7)') sgg%observation(ii)%P(i)%ZI
                     write(chari2,'(i7)') sgg%observation(ii)%P(i)%XE
                     write(charj2,'(i7)') sgg%observation(ii)%P(i)%YE
                     write(chark2,'(i7)') sgg%observation(ii)%P(i)%ZE
                      !mpidir 190319     !desrotacion para que los nombres sean correctos 
                     if      (mpidir==3) then
                        extpoint=trim(adjustl(chari)) //'_'//trim(adjustl(charj)) //'_'//trim(adjustl(chark))//'__'// &
                                 trim(adjustl(chari2))//'_'//trim(adjustl(charj2))//'_'//trim(adjustl(chark2))
                        select case (field)
                        case (iCurX)
                          prefix_field=prefix(iCurX)
                        case (iCurY)
                          prefix_field=prefix(iCurY)
                        case (iCurZ)
                          prefix_field=prefix(iCurZ)
                        case default
                          prefix_field=prefix(field)
                        end select
                     elseif  (mpidir==2) then
                        extpoint=trim(adjustl(charj)) //'_'//trim(adjustl(chark)) //'_'//trim(adjustl(chari))//'__'// &
                                 trim(adjustl(charj2))//'_'//trim(adjustl(chark2))//'_'//trim(adjustl(chari2))
                        select case (field)
                        case (iCurX)
                          prefix_field=prefix(iCurZ)
                        case (iCurY)
                          prefix_field=prefix(iCurX)
                        case (iCurZ)
                          prefix_field=prefix(iCurY)
                        case default
                          prefix_field=prefix(field)
                        end select
                     elseif  (mpidir==1) then
                          extpoint=trim(adjustl(chark)) //'_'//trim(adjustl(chari)) //'_'//trim(adjustl(charj))//'__'// &
                                   trim(adjustl(chark2))//'_'//trim(adjustl(chari2))//'_'//trim(adjustl(charj2))
                          select case (field)
                          case (iCurX)
                            prefix_field=prefix(iCurY)
                          case (iCurY)
                            prefix_field=prefix(iCurZ)
                          case (iCurZ)
                            prefix_field=prefix(iCurX)
                          case default
                            prefix_field=prefix(field)
                          end select
                     else
                          call stoponerror(layoutnumber,size,'Buggy error in mpidir. ')
                     endif
                     !
                     ext=trim(adjustl(nEntradaRoot))//'_'//trim(adjustl(sgg%observation(ii)%outputrequest))
                     !cada mpi su nombre
                     output(ii)%item(i)%path=trim(adjustl(ext))//'_'//trim(adjustl(prefix_field))//trim(adjustl(extpoint))//'.bin'

                     !
                     unit=unit+1
                     if (unit >= 2.0_RKIND  **31.0_RKIND-1.0_RKIND) then
                        call stoponerror(layoutnumber,size,'Excesive number of probes')
                     endif
                     output(ii)%item(i)%unit=unit
                     output(ii)%item(i)%columnas=0 !numero de items a escribir

                     !!!busca nombres de ficheros por duplicado y resuelve la duplicidad
                     call checkduplicatenames
                     !!!!!!
                     !precontaje
                     conta=0
                     do kkk=sgg%Observation(ii)%P(i)%ZI, sgg%Observation(ii)%P(i)%ZE
                        do jjj=sgg%observation(ii)%P(i)%YI, sgg%observation(ii)%P(i)%YE
                           do iii=sgg%observation(ii)%P(i)%XI, sgg%observation(ii)%P(i)%XE
                              if (field/=mapvtk) then
                                 if ((sgg%med(sggMiEx(III - b%Ex%XI, JJJ - b%Ex%YI, KKK - b%Ex%ZI))%Is%ThinWire).and.(iii <= SINPML_fullsize(iEx)%XE).and.(jjj <= SINPML_fullsize(iEx)%YE).and.(kkk <= SINPML_fullsize(iEx)%ZE)) then
                                    conta=conta+1
                                 endif
                                 if ((sgg%med(sggMiEy(III - b%Ey%XI, JJJ - b%Ey%YI, KKK - b%Ey%ZI))%Is%ThinWire).and.(iii <= SINPML_fullsize(iEy)%XE).and.(jjj <= SINPML_fullsize(iEy)%YE).and.(kkk <= SINPML_fullsize(iEy)%ZE)) then
                                    conta=conta+1
                                 endif
                                 if ((sgg%med(sggMiEz(III - b%Ez%XI, JJJ - b%Ez%YI, KKK - b%Ez%ZI))%Is%ThinWire).and.(iii <= SINPML_fullsize(iEz)%XE).and.(jjj <= SINPML_fullsize(iEz)%YE).and.(kkk <= SINPML_fullsize(iEz)%ZE)) then
                                    conta=conta+1
                                 endif
                                 if (((sgg%med(sggMiEx(III - b%Ex%XI, JJJ - b%Ex%YI, KKK - b%Ex%ZI))%Is%Line.AND.sgg%med(sggMiEx(III - b%Ex%XI, JJJ - b%Ex%YI, KKK - b%Ex%ZI))%Is%PEC)).and.(iii <= SINPML_fullsize(iEx)%XE).and.(jjj <= SINPML_fullsize(iEx)%YE).and.(kkk <= SINPML_fullsize(iEx)%ZE)) then
                                    if ((.not.sgg%med(sggMiHy(III - b%Hy%XI  , JJJ - b%Hy%YI  , KKK - b%Hy%ZI  ))%Is%PEC).and. &
                                    (.not.sgg%med(sggMiHy(III - b%Hy%XI  , JJJ - b%Hy%YI  , KKK - b%Hy%ZI-1))%Is%PEC).and. &
                                    (.not.sgg%med(sggMiHz(III - b%Hz%XI  , JJJ - b%Hz%YI  , KKK - b%Hz%ZI  ))%Is%PEC).and. &
                                    (.not.sgg%med(sggMiHz(III - b%Hz%XI  , JJJ - b%Hz%YI-1, KKK - b%Hz%ZI  ))%Is%PEC)) then
                                       conta=conta+1
                                    endif
                                 endif
                                 if (((sgg%med(sggMiEy(III - b%Ey%XI, JJJ - b%Ey%YI, KKK - b%Ey%ZI))%Is%Line.AND.sgg%med(sggMiEy(III - b%Ey%XI, JJJ - b%Ey%YI, KKK - b%Ey%ZI))%Is%PEC)).and.(iii <= SINPML_fullsize(iEy)%XE).and.(jjj <= SINPML_fullsize(iEy)%YE).and.(kkk <= SINPML_fullsize(iEy)%ZE)) then
                                    if ((.not.sgg%med(sggMiHz(III - b%Hz%XI  , JJJ - b%Hz%YI  , KKK - b%Hz%ZI  ))%Is%PEC).and. &
                                    (.not.sgg%med(sggMiHz(III - b%Hz%XI-1, JJJ - b%Hz%YI  , KKK - b%Hz%ZI  ))%Is%PEC).and. &
                                    (.not.sgg%med(sggMiHx(III - b%Hx%XI  , JJJ - b%Hx%YI  , KKK - b%Hx%ZI  ))%Is%PEC).and. &
                                    (.not.sgg%med(sggMiHx(III - b%Hx%XI  , JJJ - b%Hx%YI  , KKK - b%Hx%ZI-1))%Is%PEC)) then
                                       conta=conta+1
                                    endif
                                 endif
                                 if (((sgg%med(sggMiEz(III - b%Ez%XI, JJJ - b%Ez%YI, KKK - b%Ez%ZI))%Is%Line.AND.sgg%med(sggMiEz(III - b%Ez%XI, JJJ - b%Ez%YI, KKK - b%Ez%ZI))%Is%PEC)).and.(iii <= SINPML_fullsize(iEz)%XE).and.(jjj <= SINPML_fullsize(iEz)%YE).and.(kkk <= SINPML_fullsize(iEz)%ZE)) then
                                    if ((.not.sgg%med(sggMiHx(III - b%Hx%XI  , JJJ - b%Hx%YI  , KKK - b%Hx%ZI  ))%Is%PEC).and. &
                                    (.not.sgg%med(sggMiHx(III - b%Hx%XI  , JJJ - b%Hx%YI-1, KKK - b%Hx%ZI  ))%Is%PEC).and. &
                                    (.not.sgg%med(sggMiHy(III - b%Hy%XI  , JJJ - b%Hy%YI  , KKK - b%Hy%ZI  ))%Is%PEC).and. &
                                    (.not.sgg%med(sggMiHy(III - b%Hy%XI-1, JJJ - b%Hy%YI  , KKK - b%Hy%ZI  ))%Is%PEC)) then
                                       conta=conta+1
                                    endif
                                 endif
                              else !si es mapvtk
                                 !si es mapvtk y si no es vacio
                                 imed =sggMiEx(III -b%Ex%XI, JJJ- b%Ex%YI  , KKK- b%Ex%ZI  )
                                 imed1=sggMiHy(III -b%Hy%XI, JJJ- b%Hy%YI  , KKK- b%Hy%ZI  )
                                 imed2=sggMiHy(III -b%Hy%XI, JJJ- b%Hy%YI  , KKK- b%Hy%ZI-1)
                                 imed3=sggMiHz(III -b%Hz%XI, JJJ- b%Hz%YI  , KKK- b%Hz%ZI  )
                                 imed4=sggMiHz(III -b%Hz%XI, JJJ- b%Hz%YI-1, KKK- b%Hz%ZI  )
                                 call contabordes(sgg,imed,imed1,imed2,imed3,imed4,EsBorde,SINPML_fullsize,iEx,iii,jjj,kkk)
                                 if (EsBorde) then
                                    conta=conta+1
                                 endif
                                 !
                                 imed =sggMiEy(III- b%Ey%XI  , JJJ- b%Ey%YI  , KKK- b%Ey%ZI  )
                                 imed1=sggMiHz(III -b%Hz%XI  , JJJ- b%Hz%YI  , KKK- b%Hz%ZI  )
                                 imed2=sggMiHz(III -b%Hz%XI-1, JJJ- b%Hz%YI  , KKK- b%Hz%ZI  )
                                 imed3=sggMiHx(III -b%Hx%XI  , JJJ- b%Hx%YI  , KKK- b%Hx%ZI  )
                                 imed4=sggMiHx(III -b%Hx%XI  , JJJ- b%Hx%YI  , KKK- b%Hx%ZI-1)
                                 call contabordes(sgg,imed,imed1,imed2,imed3,imed4,EsBorde,SINPML_fullsize,iEy,iii,jjj,kkk)
                                 if (EsBorde) then
                                    conta=conta+1
                                 endif
                                 !
                                 imed =sggMiEz(III- b%Ez%XI  , JJJ- b%Ez%YI  , KKK- b%Ez%ZI  )
                                 imed1=sggMiHx(III -b%Hx%XI  , JJJ- b%Hx%YI  , KKK- b%Hx%ZI  )
                                 imed2=sggMiHx(III -b%Hx%XI  , JJJ- b%Hx%YI-1, KKK- b%Hx%ZI  )
                                 imed3=sggMiHy(III -b%Hy%XI  , JJJ- b%Hy%YI  , KKK- b%Hy%ZI  )
                                 imed4=sggMiHy(III -b%Hy%XI-1, JJJ- b%Hy%YI  , KKK- b%Hy%ZI  )
                                 call contabordes(sgg,imed,imed1,imed2,imed3,imed4,EsBorde,SINPML_fullsize,iEz,iii,jjj,kkk)
                                 if (EsBorde) then
                                    conta=conta+1
                                 endif
                              endif
                           end do
                        end do
                     end do

                     !!!
                     if (field==mapvtk) then
                        INIT=.TRUE.; geom=.false. ; asigna=.false.; magnetic=.false. ; electric=.true.

#ifdef CompileWithNodalSources
                        call nodalvtk(sgg,sggMiEx,sggMiEy,sggMiEz,sggMiHx,sggMiHy,sggMiHz,sggMtag, &
                        b,init,geom,asigna,electric,magnetic,conta,i,ii,output,Ntimeforvolumic)
#endif                              

#ifdef CompileWithWires
                        call wirebundlesvtk(sgg,init,geom,asigna,conta,i,ii,output,Ntimeforvolumic,wiresflavor,sggMtag)
#endif
                     endif
                     !!!
                     do kkk=sgg%Observation(ii)%P(i)%ZI, sgg%Observation(ii)%P(i)%ZE
                        do jjj=sgg%observation(ii)%P(i)%YI, sgg%observation(ii)%P(i)%YE
                           do iii=sgg%observation(ii)%P(i)%XI, sgg%observation(ii)%P(i)%XE
                              if (field/=mapvtk) then
                                 if (((sgg%med(sggMiHx(III -b%Hx%XI, JJJ- b%Hx%YI, KKK- b%Hx%ZI))%Is%PEC).or.  &
                                 (sgg%med(sggMiHx(III -b%Hx%XI, JJJ- b%Hx%YI, KKK- b%Hx%ZI))%Is%Surface).or.  &
                                 (field==icurX)).and.(iii <= SINPML_fullsize(iHx)%XE).and.(jjj <= SINPML_fullsize(iHx)%YE).and.(kkk <= SINPML_fullsize(iHx)%ZE)) then
                                    conta=conta+1
                                 endif
                                 if (((sgg%med(sggMiHy(III- b%Hy%XI, JJJ- b%Hy%YI, KKK- b%Hy%ZI))%Is%PEC).or.  &
                                 (sgg%med(sggMiHy(III- b%Hy%XI, JJJ- b%Hy%YI, KKK- b%Hy%ZI))%Is%Surface).or. &
                                 (field==icurY)).and.(iii <= SINPML_fullsize(iHy)%XE).and.(jjj <= SINPML_fullsize(iHy)%YE).and.(kkk <= SINPML_fullsize(iHy)%ZE)) then
                                    conta=conta+1
                                 endif
                                 if (((sgg%med(sggMiHz(III- b%Hz%XI, JJJ- b%Hz%YI, KKK- b%Hz%ZI))%Is%PEC).or.  &
                                 (sgg%med(sggMiHz(III- b%Hz%XI, JJJ- b%Hz%YI, KKK- b%Hz%ZI))%Is%Surface).or. &
                                 (field==icurZ)).and.(iii <= SINPML_fullsize(iHz)%XE).and.(jjj <= SINPML_fullsize(iHz)%YE).and.(kkk <= SINPML_fullsize(iHz)%ZE)) then
                                    conta=conta+1
                                 endif
                              else
                                 !si es mapvtk y si no es vacio
                                 if ((sggMiHx(III -b%Hx%XI, JJJ- b%Hx%YI, KKK- b%Hx%ZI)/=1).and. &
                                 (.not.sgg%med(sggMiHx(III -b%Hx%XI, JJJ- b%Hx%YI, KKK- b%Hx%ZI))%is%PML).and.(iii <= SINPML_fullsize(iHx)%XE).and.(jjj <= SINPML_fullsize(iHx)%YE).and.(kkk <= SINPML_fullsize(iHx)%ZE)) then
                                    conta=conta+1
                                 endif
                                 if ((sggMiHy(III- b%Hy%XI, JJJ- b%Hy%YI, KKK- b%Hy%ZI)/=1).and. &
                                 (.not.sgg%med(sggMiHy(III -b%Hy%XI, JJJ- b%Hy%YI, KKK- b%Hy%ZI))%is%PML).and.(iii <= SINPML_fullsize(iHy)%XE).and.(jjj <= SINPML_fullsize(iHy)%YE).and.(kkk <= SINPML_fullsize(iHy)%ZE)) then
                                    conta=conta+1
                                 endif
                                 if ((sggMiHz(III- b%Hz%XI, JJJ- b%Hz%YI, KKK- b%Hz%ZI)/=1).and. &
                                 (.not.sgg%med(sggMiHz(III -b%Hz%XI, JJJ- b%Hz%YI, KKK- b%Hz%ZI))%is%PML).and.(iii <= SINPML_fullsize(iHz)%XE).and.(jjj <= SINPML_fullsize(iHz)%YE).and.(kkk <= SINPML_fullsize(iHz)%ZE)) then
                                    conta=conta+1
                                 endif
                                 ! los tags de vacio negativos 141020 para mapvtk
                                 if (sggMtag(iii,jjj,kkk)<0) then
                                     if ( (btest(iabs(sggMtag(iii,jjj,kkk)),3)).and. & 
                                     (.not.sgg%med(sggMiHx(III -b%Hx%XI, JJJ- b%Hx%YI, KKK- b%Hx%ZI))%is%PML).and.(iii <= SINPML_fullsize(iHx)%XE).and.(jjj <= SINPML_fullsize(iHx)%YE).and.(kkk <= SINPML_fullsize(iHx)%ZE)) then
                                        conta=conta+1
                                     endif
                                     if ( (btest(iabs(sggMtag(iii,jjj,kkk)),4)).and. &
                                     (.not.sgg%med(sggMiHy(III -b%Hy%XI, JJJ- b%Hy%YI, KKK- b%Hy%ZI))%is%PML).and.(iii <= SINPML_fullsize(iHy)%XE).and.(jjj <= SINPML_fullsize(iHy)%YE).and.(kkk <= SINPML_fullsize(iHy)%ZE)) then
                                        conta=conta+1
                                     endif
                                     if ( (btest(iabs(sggMtag(iii,jjj,kkk)),5)).and. &
                                     (.not.sgg%med(sggMiHz(III -b%Hz%XI, JJJ- b%Hz%YI, KKK- b%Hz%ZI))%is%PML).and.(iii <= SINPML_fullsize(iHz)%XE).and.(jjj <= SINPML_fullsize(iHz)%YE).and.(kkk <= SINPML_fullsize(iHz)%ZE)) then
                                        conta=conta+1
                                     endif  
                                 endif
                              endif
                           end do
                        end do
                     end do
                     !!!
                     if (field==mapvtk) then
                        INIT=.false.; geom=.false. ; asigna=.false.; magnetic=.true. ; electric=.false.
#ifdef CompileWithNodalSources
                        call nodalvtk(sgg,sggMiEx,sggMiEy,sggMiEz,sggMiHx,sggMiHy,sggMiHz,sggMtag, &
                                      b,init,geom,asigna,electric,magnetic,conta,i,ii,output,Ntimeforvolumic)
#endif                              
                     endif
                     !!!
                     output(ii)%item(i)%columnas=conta
                     !print *,' ---- ALLOC DE output(ii)%item(i)%columnas=conta ', II,I,CONTA
                     !allocateo

                     !
                     IF (SGG%Observation(ii)%TimeDomain) THEN
                        !ojo por si algun dia esto molestara a Cray
                        !replico los ifs de transferencia y escritura
                        ntini=0
                        ntfin=0
                        first=.true.
                        do Ntime=initialtimestep,finaltimestep
                           at=sgg%tiempo(Ntime)
                           Ntimeforvolumic=Ntime !!!-nint(0.4999999+sgg%OBSERVATION(ii)%InitialTime/sgg%dt)
                           if (mod(Ntimeforvolumic,output(ii)%Trancos)==0) then
                              Ntimeforvolumic=Ntimeforvolumic/output(ii)%Trancos
                              if ( ((at >= sgg%OBSERVATION(ii)%InitialTime).and.(at <= sgg%OBSERVATION(ii)%FinalTime+sgg%dt/2.0_RKIND))) then
                                 if (first) then
                                    ntini=ntimeforvolumic
                                    first=.false.
                                 endif
                                 ntfin=ntimeforvolumic
                              endif
                           endif
                        end do
                        !!!                            ntinI=0
                        !!!                            ntfin=min(int((finaltimestep*sgg%dt-sgg%OBSERVATION(ii)%InitialTime)/sgg%dt/output(ii)%Trancos), &
                        !!!                                      int((sgg%Observation(ii)%FinalTIME-sgg%OBSERVATION(ii)%InitialTime)/sgg%dt/output(ii)%Trancos))+1


                        memo=memo+RKIND*(ntfin-ntini)*output(ii)%item(i)%columnas + 16 * output(ii)%item(i)%columnas ! 4 integers de 4 bytes

                        if (memo > MaxMemoryProbes) then
                           call stoponerror(layoutnumber,size,'ERROR: Recompile: excesive memory for the 3D probes.'// &
                           &                                   'Recompile increasing MaxMemoryProbes')
                        endif
                                                                                                                     
                        allocate (output(ii)%item(i)%Serialized%valor(ntinI : ntfin , 1:output(ii)%item(i)%columnas))      
                        !almaceno tambien los vectores
                        allocate (output(ii)%item(i)%Serialized%valor_x(ntinI : ntfin , 1:output(ii)%item(i)%columnas))
                        allocate (output(ii)%item(i)%Serialized%valor_y(ntinI : ntfin , 1:output(ii)%item(i)%columnas))
                        allocate (output(ii)%item(i)%Serialized%valor_z(ntinI : ntfin , 1:output(ii)%item(i)%columnas))   
                        output(ii)%item(i)%Serialized%Valor=0.0_RKIND
                        output(ii)%item(i)%Serialized%Valor_x=0.0_RKIND
                        output(ii)%item(i)%Serialized%Valor_y=0.0_RKIND
                        output(ii)%item(i)%Serialized%Valor_z=0.0_RKIND
                     ELSEIF (SGG%Observation(ii)%FreqDomain) THEN
                        memo=memo+RKIND*output(ii)%NumFreqs*output(ii)%item(i)%columnas + 16 * output(ii)%item(i)%columnas ! 4 integers de 4 bytes
                        if (memo > MaxMemoryProbes) then
                           call stoponerror(layoutnumber,size,'ERROR: Recompile: excesive memory for the probes.'// &
                           &                                   'Recompile increasing MaxMemoryProbes')
                        endif
                        !almaceno tambien los vectores          
                        allocate (output(ii)%item(i)%Serialized%valorComplex_x(1 : output(ii)%NumFreqs , 1:output(ii)%item(i)%columnas)) !dos posibles componentes
                        allocate (output(ii)%item(i)%Serialized%valorComplex_y(1 : output(ii)%NumFreqs , 1:output(ii)%item(i)%columnas)) !dos posibles componentes
                        allocate (output(ii)%item(i)%Serialized%valorComplex_z(1 : output(ii)%NumFreqs , 1:output(ii)%item(i)%columnas)) !dos posibles componentes
                               
                        output(ii)%item(i)%Serialized%ValorComplex_x=(0.0_RKIND,0.0_RKIND)
                        output(ii)%item(i)%Serialized%ValorComplex_y=(0.0_RKIND,0.0_RKIND)
                        output(ii)%item(i)%Serialized%ValorComplex_z=(0.0_RKIND,0.0_RKIND)
                     endif

                     allocate (output(ii)%item(i)%Serialized%eI(1:output(ii)%item(i)%columnas)) !he aprovechado la variable columnas!!jeje...comentario hecho el 120120 por bug vtk OLD
                     allocate (output(ii)%item(i)%Serialized%eJ(1:output(ii)%item(i)%columnas))
                     allocate (output(ii)%item(i)%Serialized%eK(1:output(ii)%item(i)%columnas))
                     allocate (output(ii)%item(i)%Serialized%currentType(1:output(ii)%item(i)%columnas))
                     allocate (output(ii)%item(i)%Serialized%sggMtag(1:output(ii)%item(i)%columnas))

                     !relleno info geometrica edge
                     conta=0
                     do kkk=sgg%Observation(ii)%P(i)%ZI, sgg%Observation(ii)%P(i)%ZE
                        do jjj=sgg%observation(ii)%P(i)%YI, sgg%observation(ii)%P(i)%YE
                           do iii=sgg%observation(ii)%P(i)%XI, sgg%observation(ii)%P(i)%XE
                              if (field/=mapvtk) then
                                 if ((sgg%med(sggMiEx(III - b%Ex%XI, JJJ - b%Ex%YI, KKK - b%Ex%ZI))%Is%ThinWire).and.(iii <= SINPML_fullsize(iEx)%XE).and.(jjj <= SINPML_fullsize(iEx)%YE).and.(kkk <= SINPML_fullsize(iEx)%ZE)) then
                                    conta=conta+1
                                    output(ii)%item(i)%Serialized%eI(conta)=iii
                                    output(ii)%item(i)%Serialized%eJ(conta)=jjj
                                    output(ii)%item(i)%Serialized%eK(conta)=kkk
                                    output(ii)%item(i)%Serialized%currentType(conta)=iJx
                                    output(ii)%item(i)%Serialized%sggMtag(conta)=iabs(sggMtag(iii,jjj,kkk))
                                 endif
                                 if ((sgg%med(sggMiEy(III - b%Ey%XI, JJJ - b%Ey%YI, KKK - b%Ey%ZI))%Is%ThinWire).and.(iii <= SINPML_fullsize(iEy)%XE).and.(jjj <= SINPML_fullsize(iEy)%YE).and.(kkk <= SINPML_fullsize(iEy)%ZE)) then
                                    conta=conta+1
                                    output(ii)%item(i)%Serialized%eI(conta)=iii
                                    output(ii)%item(i)%Serialized%eJ(conta)=jjj
                                    output(ii)%item(i)%Serialized%eK(conta)=kkk
                                    output(ii)%item(i)%Serialized%currentType(conta)=iJy
                                    output(ii)%item(i)%Serialized%sggMtag(conta)=iabs(sggMtag(iii,jjj,kkk))
                                 endif
                                 if ((sgg%med(sggMiEz(III - b%Ez%XI, JJJ - b%Ez%YI, KKK - b%Ez%ZI))%Is%ThinWire).and.(iii <= SINPML_fullsize(iEz)%XE).and.(jjj <= SINPML_fullsize(iEz)%YE).and.(kkk <= SINPML_fullsize(iEz)%ZE)) then
                                    conta=conta+1
                                    output(ii)%item(i)%Serialized%eI(conta)=iii
                                    output(ii)%item(i)%Serialized%eJ(conta)=jjj
                                    output(ii)%item(i)%Serialized%eK(conta)=kkk
                                    output(ii)%item(i)%Serialized%currentType(conta)=iJz
                                    output(ii)%item(i)%Serialized%sggMtag(conta)=iabs(sggMtag(iii,jjj,kkk))
                                 endif
                                 if (((sgg%med(sggMiEx(III - b%Ex%XI, JJJ - b%Ex%YI, KKK - b%Ex%ZI))%Is%Line.AND.sgg%med(sggMiEx(III - b%Ex%XI, JJJ - b%Ex%YI, KKK - b%Ex%ZI))%Is%PEC)).and.(iii <= SINPML_fullsize(iEx)%XE).and.(jjj <= SINPML_fullsize(iEx)%YE).and.(kkk <= SINPML_fullsize(iEx)%ZE)) then
                                    if ((.not.sgg%med(sggMiHy(III - b%Hy%XI  , JJJ - b%Hy%YI  , KKK - b%Hy%ZI  ))%Is%PEC).and. &
                                    (.not.sgg%med(sggMiHy(III - b%Hy%XI  , JJJ - b%Hy%YI  , KKK - b%Hy%ZI-1))%Is%PEC).and. &
                                    (.not.sgg%med(sggMiHz(III - b%Hz%XI  , JJJ - b%Hz%YI  , KKK - b%Hz%ZI  ))%Is%PEC).and. &
                                    (.not.sgg%med(sggMiHz(III - b%Hz%XI  , JJJ - b%Hz%YI-1, KKK - b%Hz%ZI  ))%Is%PEC)) then
                                       conta=conta+1
                                       output(ii)%item(i)%Serialized%eI(conta)=iii
                                       output(ii)%item(i)%Serialized%eJ(conta)=jjj
                                       output(ii)%item(i)%Serialized%eK(conta)=kkk
                                       output(ii)%item(i)%Serialized%currentType(conta)=iJx
                                    output(ii)%item(i)%Serialized%sggMtag(conta)=iabs(sggMtag(iii,jjj,kkk))
                                    endif
                                 endif
                                 if (((sgg%med(sggMiEy(III - b%Ey%XI, JJJ - b%Ey%YI, KKK - b%Ey%ZI))%Is%Line.AND.sgg%med(sggMiEy(III - b%Ey%XI, JJJ - b%Ey%YI, KKK - b%Ey%ZI))%Is%PEC)).and.(iii <= SINPML_fullsize(iEy)%XE).and.(jjj <= SINPML_fullsize(iEy)%YE).and.(kkk <= SINPML_fullsize(iEy)%ZE)) then
                                    if ((.not.sgg%med(sggMiHz(III - b%Hz%XI  , JJJ - b%Hz%YI  , KKK - b%Hz%ZI  ))%Is%PEC).and. &
                                    (.not.sgg%med(sggMiHz(III - b%Hz%XI-1, JJJ - b%Hz%YI  , KKK - b%Hz%ZI  ))%Is%PEC).and. &
                                    (.not.sgg%med(sggMiHx(III - b%Hx%XI  , JJJ - b%Hx%YI  , KKK - b%Hx%ZI  ))%Is%PEC).and. &
                                    (.not.sgg%med(sggMiHx(III - b%Hx%XI  , JJJ - b%Hx%YI  , KKK - b%Hx%ZI-1))%Is%PEC)) then
                                       conta=conta+1
                                       output(ii)%item(i)%Serialized%eI(conta)=iii
                                       output(ii)%item(i)%Serialized%eJ(conta)=jjj
                                       output(ii)%item(i)%Serialized%eK(conta)=kkk
                                       output(ii)%item(i)%Serialized%currentType(conta)=iJy
                                    output(ii)%item(i)%Serialized%sggMtag(conta)=iabs(sggMtag(iii,jjj,kkk))
                                    endif
                                 endif
                                 if (((sgg%med(sggMiEz(III - b%Ez%XI, JJJ - b%Ez%YI, KKK - b%Ez%ZI))%Is%Line.AND.sgg%med(sggMiEz(III - b%Ez%XI, JJJ - b%Ez%YI, KKK - b%Ez%ZI))%Is%PEC)).and.(iii <= SINPML_fullsize(iEz)%XE).and.(jjj <= SINPML_fullsize(iEz)%YE).and.(kkk <= SINPML_fullsize(iEz)%ZE)) then
                                    if ((.not.sgg%med(sggMiHx(III - b%Hx%XI  , JJJ - b%Hx%YI  , KKK - b%Hx%ZI  ))%Is%PEC).and. &
                                    (.not.sgg%med(sggMiHx(III - b%Hx%XI  , JJJ - b%Hx%YI-1, KKK - b%Hx%ZI  ))%Is%PEC).and. &
                                    (.not.sgg%med(sggMiHy(III - b%Hy%XI  , JJJ - b%Hy%YI  , KKK - b%Hy%ZI  ))%Is%PEC).and. &
                                    (.not.sgg%med(sggMiHy(III - b%Hy%XI-1, JJJ - b%Hy%YI  , KKK - b%Hy%ZI  ))%Is%PEC)) then
                                       conta=conta+1
                                       output(ii)%item(i)%Serialized%eI(conta)=iii
                                       output(ii)%item(i)%Serialized%eJ(conta)=jjj
                                       output(ii)%item(i)%Serialized%eK(conta)=kkk
                                       output(ii)%item(i)%Serialized%currentType(conta)=iJz
                                    output(ii)%item(i)%Serialized%sggMtag(conta)=iabs(sggMtag(iii,jjj,kkk))
                                    endif
                                 endif
                              else !si es mapvtk
                                 !si es mapvtk y si no es vacio

                                 imed =sggMiEx(III -b%Ex%XI, JJJ- b%Ex%YI  , KKK- b%Ex%ZI  )
                                 imed1=sggMiHy(III -b%Hy%XI, JJJ- b%Hy%YI  , KKK- b%Hy%ZI  )
                                 imed2=sggMiHy(III -b%Hy%XI, JJJ- b%Hy%YI  , KKK- b%Hy%ZI-1)
                                 imed3=sggMiHz(III -b%Hz%XI, JJJ- b%Hz%YI  , KKK- b%Hz%ZI  )
                                 imed4=sggMiHz(III -b%Hz%XI, JJJ- b%Hz%YI-1, KKK- b%Hz%ZI  )
                                 call contabordes(sgg,imed,imed1,imed2,imed3,imed4,EsBorde,SINPML_fullsize,iEx,iii,jjj,kkk)
                                 if (EsBorde) then
                                    conta=conta+1
                                    output(ii)%item(i)%Serialized%eI(conta)=iii
                                    output(ii)%item(i)%Serialized%eJ(conta)=jjj
                                    output(ii)%item(i)%Serialized%eK(conta)=kkk
                                    output(ii)%item(i)%Serialized%currentType(conta)=iJx !las lineas las asimilo a corrientes para que salgan en edges
                                    output(ii)%item(i)%Serialized%sggMtag(conta)=sggMtag(iii,jjj,kkk) !sin valor absoluto pq es mapvtk
                                 endif
                                 imed =sggMiEy(III- b%Ey%XI  , JJJ- b%Ey%YI  , KKK- b%Ey%ZI  )
                                 imed1=sggMiHz(III -b%Hz%XI  , JJJ- b%Hz%YI  , KKK- b%Hz%ZI  )
                                 imed2=sggMiHz(III -b%Hz%XI-1, JJJ- b%Hz%YI  , KKK- b%Hz%ZI  )
                                 imed3=sggMiHx(III -b%Hx%XI  , JJJ- b%Hx%YI  , KKK- b%Hx%ZI  )
                                 imed4=sggMiHx(III -b%Hx%XI  , JJJ- b%Hx%YI  , KKK- b%Hx%ZI-1)
                                 call contabordes(sgg,imed,imed1,imed2,imed3,imed4,EsBorde,SINPML_fullsize,iEy,iii,jjj,kkk)
                                 if (EsBorde) then
                                    conta=conta+1
                                    output(ii)%item(i)%Serialized%eI(conta)=iii
                                    output(ii)%item(i)%Serialized%eJ(conta)=jjj
                                    output(ii)%item(i)%Serialized%eK(conta)=kkk
                                    output(ii)%item(i)%Serialized%currentType(conta)=iJy !las lineas las asimilo a corrientes para que salgan en edges
                                    output(ii)%item(i)%Serialized%sggMtag(conta)=sggMtag(iii,jjj,kkk)
                                 endif
                                 imed =sggMiEz(III- b%Ez%XI  , JJJ- b%Ez%YI  , KKK- b%Ez%ZI  )
                                 imed1=sggMiHx(III -b%Hx%XI  , JJJ- b%Hx%YI  , KKK- b%Hx%ZI  )
                                 imed2=sggMiHx(III -b%Hx%XI  , JJJ- b%Hx%YI-1, KKK- b%Hx%ZI  )
                                 imed3=sggMiHy(III -b%Hy%XI  , JJJ- b%Hy%YI  , KKK- b%Hy%ZI  )
                                 imed4=sggMiHy(III -b%Hy%XI-1, JJJ- b%Hy%YI  , KKK- b%Hy%ZI  )
                                 call contabordes(sgg,imed,imed1,imed2,imed3,imed4,EsBorde,SINPML_fullsize,iEz,iii,jjj,kkk)
                                 if (EsBorde) then
                                    conta=conta+1
                                    output(ii)%item(i)%Serialized%eI(conta)=iii
                                    output(ii)%item(i)%Serialized%eJ(conta)=jjj
                                    output(ii)%item(i)%Serialized%eK(conta)=kkk
                                    output(ii)%item(i)%Serialized%currentType(conta)=iJz !las lineas las asimilo a corrientes para que salgan en edges
                                    output(ii)%item(i)%Serialized%sggMtag(conta)=sggMtag(iii,jjj,kkk)
                                 endif
                              endif
                              !
                           end do
                        end do
                     end do

                     !!!
                     if (field==mapvtk) then
                        INIT=.false.; geom=.true. ; asigna=.false.; magnetic=.false. ; electric=.true.
#ifdef CompileWithNodalSources
                        call nodalvtk(sgg,sggMiEx,sggMiEy,sggMiEz,sggMiHx,sggMiHy,sggMiHz,sggMtag,&
                                      b,init,geom,asigna,electric,magnetic,conta,i,ii,output,Ntimeforvolumic)
#endif                              
#ifdef CompileWithWires
                        call wirebundlesvtk(sgg,init,geom,asigna,conta,i,ii,output,Ntimeforvolumic,wiresflavor,sggMtag)
#endif
                     endif
                     !!!
                     do kkk=sgg%Observation(ii)%P(i)%ZI, sgg%Observation(ii)%P(i)%ZE
                        do jjj=sgg%observation(ii)%P(i)%YI, sgg%observation(ii)%P(i)%YE
                           do iii=sgg%observation(ii)%P(i)%XI, sgg%observation(ii)%P(i)%XE
                              if (field/=mapvtk) then
                                 if (((sgg%med(sggMiHx(III -b%Hx%XI, JJJ- b%Hx%YI, KKK- b%Hx%ZI))%Is%PEC).or. &
                                 (sgg%med(sggMiHx(III -b%Hx%XI, JJJ- b%Hx%YI, KKK- b%Hx%ZI))%Is%Surface).or. &
                                 (field==icurX)).and.(iii <= SINPML_fullsize(iHx)%XE).and.(jjj <= SINPML_fullsize(iHx)%YE).and.(kkk <= SINPML_fullsize(iHx)%ZE)) then
                                    conta=conta+1
                                    output(ii)%item(i)%Serialized%eI(conta)=iii
                                    output(ii)%item(i)%Serialized%eJ(conta)=jjj
                                    output(ii)%item(i)%Serialized%eK(conta)=kkk
                                    output(ii)%item(i)%Serialized%currentType(conta)=iBloqueJx
                                    output(ii)%item(i)%Serialized%sggMtag(conta)=iabs(sggMtag(iii,jjj,kkk))
                                 endif
                                 if (((sgg%med(sggMiHy(III- b%Hy%XI, JJJ- b%Hy%YI, KKK- b%Hy%ZI))%Is%PEC).or. &
                                 (sgg%med(sggMiHy(III- b%Hy%XI, JJJ- b%Hy%YI, KKK- b%Hy%ZI))%Is%Surface).or. &
                                 (field==icurY)).and.(iii <= SINPML_fullsize(iHy)%XE).and.(jjj <= SINPML_fullsize(iHy)%YE).and.(kkk <= SINPML_fullsize(iHy)%ZE)) then
                                    conta=conta+1
                                    output(ii)%item(i)%Serialized%eI(conta)=iii
                                    output(ii)%item(i)%Serialized%eJ(conta)=jjj
                                    output(ii)%item(i)%Serialized%eK(conta)=kkk
                                    output(ii)%item(i)%Serialized%currentType(conta)=iBloqueJy
                                    output(ii)%item(i)%Serialized%sggMtag(conta)=iabs(sggMtag(iii,jjj,kkk))
                                 endif
                                 if (((sgg%med(sggMiHz(III- b%Hz%XI, JJJ- b%Hz%YI, KKK- b%Hz%ZI))%Is%PEC).or. &
                                 (sgg%med(sggMiHz(III- b%Hz%XI, JJJ- b%Hz%YI, KKK- b%Hz%ZI))%Is%Surface).or. &
                                 (field==icurZ)).and.(iii <= SINPML_fullsize(iHz)%XE).and.(jjj <= SINPML_fullsize(iHz)%YE).and.(kkk <= SINPML_fullsize(iHz)%ZE)) then
                                    conta=conta+1
                                    output(ii)%item(i)%Serialized%eI(conta)=iii
                                    output(ii)%item(i)%Serialized%eJ(conta)=jjj
                                    output(ii)%item(i)%Serialized%eK(conta)=kkk
                                    output(ii)%item(i)%Serialized%currentType(conta)=iBloqueJz
                                    output(ii)%item(i)%Serialized%sggMtag(conta)=iabs(sggMtag(iii,jjj,kkk))
                                 endif
                              else !mapvtk y si no es vacio, asimilo la salida a corrientes iBloqueJ? para que vtk.f90 los escriba en quads
                                 if ((sggMiHx(III -b%Hx%XI, JJJ- b%Hx%YI, KKK- b%Hx%ZI)/=1).and. &
                                 (.not.sgg%med(sggMiHx(III -b%Hx%XI, JJJ- b%Hx%YI, KKK- b%Hx%ZI))%is%PML).and.(iii <= SINPML_fullsize(iHx)%XE).and.(jjj <= SINPML_fullsize(iHx)%YE).and.(kkk <= SINPML_fullsize(iHx)%ZE)) then
                                    conta=conta+1
                                    output(ii)%item(i)%Serialized%eI(conta)=iii
                                    output(ii)%item(i)%Serialized%eJ(conta)=jjj
                                    output(ii)%item(i)%Serialized%eK(conta)=kkk
                                    output(ii)%item(i)%Serialized%currentType(conta)=iBloqueJx
                                    output(ii)%item(i)%Serialized%sggMtag(conta)=sggMtag(iii,jjj,kkk)  !sin valor absoluto pq es mapvtk
                                 endif
                                 if ((sggMiHy(III- b%Hy%XI, JJJ- b%Hy%YI, KKK- b%Hy%ZI)/=1).and. &
                                 (.not.sgg%med(sggMiHy(III -b%Hy%XI, JJJ- b%Hy%YI, KKK- b%Hy%ZI))%is%PML).and.(iii <= SINPML_fullsize(iHy)%XE).and.(jjj <= SINPML_fullsize(iHy)%YE).and.(kkk <= SINPML_fullsize(iHy)%ZE)) then
                                    conta=conta+1
                                    output(ii)%item(i)%Serialized%eI(conta)=iii
                                    output(ii)%item(i)%Serialized%eJ(conta)=jjj
                                    output(ii)%item(i)%Serialized%eK(conta)=kkk
                                    output(ii)%item(i)%Serialized%currentType(conta)=iBloqueJy
                                    output(ii)%item(i)%Serialized%sggMtag(conta)=sggMtag(iii,jjj,kkk)
                                 endif
                                 if ((sggMiHz(III- b%Hz%XI, JJJ- b%Hz%YI, KKK- b%Hz%ZI)/=1).and. &
                                 (.not.sgg%med(sggMiHz(III -b%Hz%XI, JJJ- b%Hz%YI, KKK- b%Hz%ZI))%is%PML).and.(iii <= SINPML_fullsize(iHz)%XE).and.(jjj <= SINPML_fullsize(iHz)%YE).and.(kkk <= SINPML_fullsize(iHz)%ZE)) then
                                    conta=conta+1
                                    output(ii)%item(i)%Serialized%eI(conta)=iii
                                    output(ii)%item(i)%Serialized%eJ(conta)=jjj
                                    output(ii)%item(i)%Serialized%eK(conta)=kkk
                                    output(ii)%item(i)%Serialized%currentType(conta)=iBloqueJz
                                    output(ii)%item(i)%Serialized%sggMtag(conta)=sggMtag(iii,jjj,kkk)
                                 endif
!                                 ! los tags 141020 para mapvtk
                                 if (sggMtag(iii,jjj,kkk)<0) then
                                     if ( (btest(iabs(sggMtag(iii,jjj,kkk)),3)).and. & 
                                     (.not.sgg%med(sggMiHx(III -b%Hx%XI, JJJ- b%Hx%YI, KKK- b%Hx%ZI))%is%PML).and.(iii <= SINPML_fullsize(iHx)%XE).and.(jjj <= SINPML_fullsize(iHx)%YE).and.(kkk <= SINPML_fullsize(iHx)%ZE)) then
                                        conta=conta+1
                                        output(ii)%item(i)%Serialized%eI(conta)=iii
                                        output(ii)%item(i)%Serialized%eJ(conta)=jjj
                                        output(ii)%item(i)%Serialized%eK(conta)=kkk
                                        output(ii)%item(i)%Serialized%currentType(conta)=iBloqueJx
                                        output(ii)%item(i)%Serialized%sggMtag(conta)=sggMtag(iii,jjj,kkk) !sin valor absoluto pq es mapvtk
                                     endif
                                     if ( (btest(iabs(sggMtag(iii,jjj,kkk)),4)).and. & 
                                     (.not.sgg%med(sggMiHy(III -b%Hy%XI, JJJ- b%Hy%YI, KKK- b%Hy%ZI))%is%PML).and.(iii <= SINPML_fullsize(iHy)%XE).and.(jjj <= SINPML_fullsize(iHy)%YE).and.(kkk <= SINPML_fullsize(iHy)%ZE)) then
                                        conta=conta+1
                                        output(ii)%item(i)%Serialized%eI(conta)=iii
                                        output(ii)%item(i)%Serialized%eJ(conta)=jjj
                                        output(ii)%item(i)%Serialized%eK(conta)=kkk
                                        output(ii)%item(i)%Serialized%currentType(conta)=iBloqueJy
                                        output(ii)%item(i)%Serialized%sggMtag(conta)=sggMtag(iii,jjj,kkk)
                                     endif
                                     if ( (btest(iabs(sggMtag(iii,jjj,kkk)),5)).and. &
                                     (.not.sgg%med(sggMiHz(III -b%Hz%XI, JJJ- b%Hz%YI, KKK- b%Hz%ZI))%is%PML).and.(iii <= SINPML_fullsize(iHz)%XE).and.(jjj <= SINPML_fullsize(iHz)%YE).and.(kkk <= SINPML_fullsize(iHz)%ZE)) then
                                        conta=conta+1
                                        output(ii)%item(i)%Serialized%eI(conta)=iii
                                        output(ii)%item(i)%Serialized%eJ(conta)=jjj
                                        output(ii)%item(i)%Serialized%eK(conta)=kkk
                                        output(ii)%item(i)%Serialized%currentType(conta)=iBloqueJz
                                        output(ii)%item(i)%Serialized%sggMtag(conta)=sggMtag(iii,jjj,kkk)
                                     endif
                                endif
                              endif
                              !
                           end do
                        end do
                     end do

                     !!!
                     if (field==mapvtk) then
                        INIT=.false.; geom=.true. ; asigna=.false.; magnetic=.true. ; electric=.false.
#ifdef CompileWithNodalSources
                        call nodalvtk(sgg,sggMiEx,sggMiEy,sggMiEz,sggMiHx,sggMiHy,sggMiHz,sggMtag,&
                                      b,init,geom,asigna,electric,magnetic,conta,i,ii,output,Ntimeforvolumic)
#endif                              
                     endif
                     !!!
                     my_iostat=0
9137                 if(my_iostat /= 0) write(*,fmt='(a)',advance='no'), '.' !!if(my_iostat /= 0) print '(i5,a1,i4,2x,a)',9137,layoutnumber,trim(adjustl(nEntradaRoot))//'_Outputrequests_'//trim(adjustl(whoamishort))//'.txt'
                     write(19,'(a)',err=9137,iostat=my_iostat) trim(adjustl(output(ii)%item(i)%path))

                     !erase pre-existing data unless this is a resuming simulation

                     if (.not.resume) then
                        if (SGG%Observation(ii)%TimeDomain) then
                           open ( output(ii)%item(i)%unit,file= trim(adjustl(output(ii)%item(i)%path)),form='unformatted')
                           write(output(ii)%item(i)%unit) '!END'
                           close (output(ii)%item(i)%unit,status='DELETE')
                           my_iostat=0
9240                       if(my_iostat /= 0) write(*,fmt='(a)',advance='no'), '.' !!if(my_iostat /= 0) print '(i5,a1,i4,2x,a,2x,i5)',9240,'.',layoutnumber,trim(adjustl(output(ii)%item(i)%path)),output(ii)%item(i)%unit         
                           open ( output(ii)%item(i)%unit,file= trim(adjustl(output(ii)%item(i)%path)),form='unformatted',err=9240,iostat=my_iostat,status='new',action='write') !����
                               
                           write(output(ii)%item(i)%unit) output(ii)%item(i)%columnas
                           do conta=1,output(ii)%item(i)%columnas
                              write(output(ii)%item(i)%unit) &
                              output(ii)%item(i)%Serialized%eI(conta), &
                              output(ii)%item(i)%Serialized%eJ(conta), &
                              output(ii)%item(i)%Serialized%eK(conta), &
                              output(ii)%item(i)%Serialized%currentType(conta), &
                              output(ii)%item(i)%Serialized%sggMtag(conta) !added to resuming file 121020  
                           end do
                        elseif (SGG%Observation(ii)%FreqDomain) then
                           open ( output(ii)%item(i)%unit,file= trim(adjustl(output(ii)%item(i)%path)),form='unformatted')
                           write(output(ii)%item(i)%unit) '!END'
                           close (output(ii)%item(i)%unit,status='DELETE')
                        endif !no need to keep it open
                        !wipe out duplicate data after non synchronous data and field resuming !later
                     else !SE RESUMEA
                        inquire(file= trim(adjustl(output(ii)%item(i)%path)),exist=existe)
                        if (.not.existe) then
                           call stoponerror(layoutnumber,size,'Data files for resuming non existent (Volume) '//trim(adjustl(output(ii)%item(i)%path)))
                        endif
                        open (output(ii)%item(i)%unit,access='sequential',file= trim(adjustl(output(ii)%item(i)%path)),  &
                        form='unformatted')
     
                        if ((SGG%Observation(ii)%TimeDomain).and.(sgg%observation(ii)%P(1)%what /= mapvtk)) then
                           !
                           read(output(ii)%item(i)%unit) ndum
                           if (output(ii)%item(i)%columnas/=ndum) call stoponerror(layoutnumber,size,'BUGGYError reading resuming files () ')
                           do conta=1,output(ii)%item(i)%columnas
                              read(output(ii)%item(i)%unit) ndum, ndum, ndum, ndum, ndum
                           end do
                           cutting3 :  do
                              read(output(ii)%item(i)%unit,end=699) at 
                              if (output(ii)%item(i)%columnas/=0) read(output(ii)%item(i)%unit,end=699) (rdum,conta=1,output(ii)%item(i)%columnas)
                              output(ii)%TimesWritten=output(ii)%TimesWritten+1
                              if (at > lastexecutedtime) then
                                 print '(i4,a,a,2e19.9e3)',quienmpi,'Cutting 3 ',trim(adjustl(output(ii)%item(i)%path)),at,lastexecutedtime
                                 backspace(output(ii)%item(i)%unit)
                                 endfile(output(ii)%item(i)%unit)
                                 exit cutting3
                              endif
                           end do cutting3
699                        continue
                           !!            backspace(output(ii)%item(i)%unit) !machaco el timeswritten proque solo puede haber uno al final
                           close (output(ii)%item(i)%unit)
                           open (output(ii)%item(i)%unit,file= trim(adjustl(output(ii)%item(i)%path)),position='append',  &
                           form='unformatted')
                        elseif ((SGG%Observation(ii)%FreqDomain).and.(sgg%observation(ii)%P(1)%what /= mapvtk)) then
                           !
                           read(output(ii)%item(i)%unit) ndum
                           do conta=1,output(ii)%item(i)%columnas
                              read(output(ii)%item(i)%unit) ndum, ndum, ndum, ndum, ndum
                           end do
                           read(output(ii)%item(i)%unit) at
                           if (int(at/sgg%dt) /= (initialtimestep-1)) then
                              write (buff,*) nint(at/sgg%dt), initialtimestep-1,' Data files for resuming 3D freq domain probes might be corrupt. Continuing....'
                              call print11(layoutnumber,buff)
                           endif
                           do N=1,output(ii)%NumFreqs
                              read(output(ii)%item(i)%unit,end=6919) rdum
                              if (output(ii)%item(i)%columnas/=0) then
                                     do conta=1,output(ii)%item(i)%columnas
                                         read(output(ii)%item(i)%unit,end=6919) &   
                                           output(ii)%item(i)%Serialized%ValorComplex_x(N,conta),&
                                           output(ii)%item(i)%Serialized%ValorComplex_y(N,conta),&
                                           output(ii)%item(i)%Serialized%ValorComplex_z(N,conta) 
                                     end do
                              endif
                              if (SGG%Observation(ii)%transfer)then   
                                  output(ii)%item(i)%Serialized%ValorComplex_x = output(ii)%item(i)%Serialized%ValorComplex_x * output(ii)%dftEntrada(N) !desnormaliza
                                  output(ii)%item(i)%Serialized%ValorComplex_y = output(ii)%item(i)%Serialized%ValorComplex_y * output(ii)%dftEntrada(N) !desnormaliza
                                  output(ii)%item(i)%Serialized%ValorComplex_z = output(ii)%item(i)%Serialized%ValorComplex_z * output(ii)%dftEntrada(N) !desnormaliza
                              endif
                           end do
6919                       continue
                           close (output(ii)%item(i)%unit,status='delete')
                        endif
                     endif
                  endif
                  !!!!!!!!!!!!!!!!!!fin vtk
                  !Volumic probes
                case (iMEC,iMHC,iExC,iEyC,iEzC,iHxC,iHyC,iHzC)
                  if (sgg%Observation(ii)%Volumic) then !they are necssaryly
                     if (sgg%Observation(ii)%nP /= 1) then
                        call stoponerror(layoutnumber,size,'ERROR! More than a volumic probe per group')
                     endif
                     !readjust correctly the calculation region
                     select case(field)
                      case (iExC,iEyC,iHzC,iMhC)
                        sgg%observation(ii)%P(i)%ZI=max(sgg%Sweep(fieldo(field,'Z'))%ZI  ,sgg%observation(ii)%P(i)%ZI)
                        sgg%observation(ii)%P(i)%ZE=min(sgg%Sweep(fieldo(field,'Z'))%ZE-1  ,sgg%observation(ii)%P(i)%ZE)
                      case (iEzC,iHxC,iHyC,iMeC)
                        sgg%observation(ii)%P(i)%ZI=max(sgg%Sweep(fieldo(field,'Z'))%ZI  ,sgg%observation(ii)%P(i)%ZI)
                        sgg%observation(ii)%P(i)%ZE=min(sgg%Sweep(fieldo(field,'Z'))%ZE  ,sgg%observation(ii)%P(i)%ZE)
                     end select
                     !solo acepto que P(1:1) !!!
                     write(chari ,'(i7)') sgg%observation(ii)%P(i)%XI
                     write(charj ,'(i7)') sgg%observation(ii)%P(i)%YI
                     write(chark ,'(i7)') sgg%observation(ii)%P(i)%ZI
                     write(chari2,'(i7)') sgg%observation(ii)%P(i)%XE
                     write(charj2,'(i7)') sgg%observation(ii)%P(i)%YE
                     write(chark2,'(i7)') sgg%observation(ii)%P(i)%ZE
                     
                  !mpidir 190319      !desrotacion para que los nombres sean correctos 
                  if      (mpidir==3) then
                      extpoint=trim(adjustl(chari)) //'_'//trim(adjustl(charj)) //'_'//trim(adjustl(chark))//'__'// &
                               trim(adjustl(chari2))//'_'//trim(adjustl(charj2))//'_'//trim(adjustl(chark2))
                      select case (field)
                      case (iExC)
                        prefix_field=prefix(iExC)
                      case (iEyC)
                        prefix_field=prefix(iEyC)
                      case (iEzC)
                        prefix_field=prefix(iEzC)
                      case (iHxC)
                        prefix_field=prefix(iHxC)
                      case (iHyC)
                        prefix_field=prefix(iHyC)
                      case (iHzC)
                        prefix_field=prefix(iHzC)
                      case default
                        prefix_field=prefix(field)
                      end select
                  elseif  (mpidir==2) then
                      extpoint=trim(adjustl(charj)) //'_'//trim(adjustl(chark)) //'_'//trim(adjustl(chari))//'__'// &
                               trim(adjustl(charj2))//'_'//trim(adjustl(chark2))//'_'//trim(adjustl(chari2))
                      select case (field)
                      case (iExC)
                        prefix_field=prefix(iEzC)
                      case (iEyC)
                        prefix_field=prefix(iExC)
                      case (iEzC)
                        prefix_field=prefix(iEyC)
                      case (iHxC)
                        prefix_field=prefix(iHzC)
                      case (iHyC)
                        prefix_field=prefix(iHxC)
                      case (iHzC)
                        prefix_field=prefix(iHyC)
                      case default
                        prefix_field=prefix(field)
                      end select
                  elseif  (mpidir==1) then
                      extpoint=trim(adjustl(chark)) //'_'//trim(adjustl(chari)) //'_'//trim(adjustl(charj))//'__'// &
                               trim(adjustl(chark2))//'_'//trim(adjustl(chari2))//'_'//trim(adjustl(charj2))
                      select case (field)
                      case (iExC)
                        prefix_field=prefix(iEyC)
                      case (iEyC)
                        prefix_field=prefix(iEzC)
                      case (iEzC)
                        prefix_field=prefix(iExC)
                      case (iHxC)
                        prefix_field=prefix(iHyC)
                      case (iHyC)
                        prefix_field=prefix(iHzC)
                      case (iHzC)
                        prefix_field=prefix(iHxC)
                      case default
                        prefix_field=prefix(field)
                      end select
                  else
                      call stoponerror(layoutnumber,size,'Buggy error in mpidir. ')
                  endif
                  !
                     !
                     ext=trim(adjustl(nEntradaRoot))//'_'//trim(adjustl(sgg%observation(ii)%outputrequest))
                     !cada mpi su nombre
                     output(ii)%item(i)%path=trim(adjustl(ext))//'_'//trim(adjustl(prefix_field))//trim(adjustl(extpoint))//'.bin'

                     !
                     unit=unit+1
                     if (unit >= 2.0_RKIND  **31.0_RKIND-1.0_RKIND) then
                        call stoponerror(layoutnumber,size,'Excesive number of probes')
                     endif
                     output(ii)%item(i)%unit=unit

                     !!!busca nombres de ficheros por duplicado y resuelve la duplicidad
                     call checkduplicatenames
                     !!!!!!
                     !
                     output(ii)%item(i)%columnas=(sgg%Observation(ii)%P(i)%XE-sgg%Observation(ii)%P(i)%XI+1)* &
                     (sgg%observation(ii)%P(i)%YE-sgg%observation(ii)%P(i)%YI+1)* &
                     (sgg%observation(ii)%P(i)%ZE-sgg%observation(ii)%P(i)%ZI+1)
                     !
                     !ojo por si algun dia esto molestara a Cray
                     IF (SGG%Observation(ii)%TimeDomain) THEN

                        !replico los ifs de transferencia y escritura
                        ntini=0
                        ntfin=0
                        first=.true.
                        do Ntime=initialtimestep,finaltimestep
                           at=sgg%tiempo(Ntime)
                           Ntimeforvolumic=Ntime !!!-nint(0.4999999+sgg%OBSERVATION(ii)%InitialTime/sgg%dt)
                           if (mod(Ntimeforvolumic,output(ii)%Trancos)==0) then
                              Ntimeforvolumic=Ntimeforvolumic/output(ii)%Trancos
                              if ( ((at >= sgg%OBSERVATION(ii)%InitialTime).and.(at <= sgg%OBSERVATION(ii)%FinalTime+sgg%dt/2))) then
                                 if (first) then
                                    ntini=ntimeforvolumic
                                    first=.false.
                                 endif
                                 ntfin=ntimeforvolumic
                              endif
                           endif
                        end do

                        !!!                             ntinI=0
                        !!!                             ntfin=min(int((finaltimestep*sgg%dt-sgg%OBSERVATION(ii)%InitialTime)/sgg%dt/output(ii)%Trancos), &
                        !!!                                       int((sgg%Observation(ii)%FinalTIME-sgg%OBSERVATION(ii)%InitialTime)/sgg%dt/output(ii)%Trancos))+1
                        
                        memo=memo+RKIND* &
                            (ntfin-ntini+1) * &
                        (output(ii)%item(i)%XEtrancos-output(ii)%item(i)%XItrancos+1)* &
                        (output(ii)%item(i)%YEtrancos-output(ii)%item(i)%YItrancos+1)* &
                        (output(ii)%item(i)%ZEtrancos-output(ii)%item(i)%ZItrancos+1)

                        if (memo > MaxMemoryProbes) then
                           call stoponerror(layoutnumber,size,'ERROR: Recompile: excesive memory for the 3D probes.'// &
                           &                                   'Recompile increasing MaxMemoryProbes')
                        endif

                       
                        allocate (output(ii)%item(i)%valor3D(ntinI : ntfin, &
                            output(ii)%item(i)%XItrancos : output(ii)%item(i)%XEtrancos, &
                            output(ii)%item(i)%YItrancos : output(ii)%item(i)%YEtrancos, &
                            output(ii)%item(i)%ZItrancos : output(ii)%item(i)%ZEtrancos) )
                        output(ii)%item(i)%valor3D=0.0_RKIND

                     ELSEIF (SGG%Observation(ii)%FreqDomain) THEN
                        memo=memo+RKIND*output(ii)%NumFreqs*output(ii)%item(i)%columnas + 16 * output(ii)%item(i)%columnas ! 4 integers de 4 bytes
                        if (memo > MaxMemoryProbes) then
                           call stoponerror(layoutnumber,size,'ERROR: Recompile: excesive memory for the probes.'// &
                           &                                   'Recompile increasing MaxMemoryProbes')
                        endif

                        allocate (output(ii)%item(i)%valor3DComplex(1 : output(ii)%NumFreqs,1:3, & !tres posibles componentes
                            output(ii)%item(i)%XItrancos : output(ii)%item(i)%XEtrancos, &
                            output(ii)%item(i)%YItrancos : output(ii)%item(i)%YEtrancos, &
                            output(ii)%item(i)%ZItrancos : output(ii)%item(i)%ZEtrancos) )
                        output(ii)%item(i)%valor3DComplex=(0.0_RKIND,0.0_RKIND)
                     endif
                     !            
                     my_iostat=0
9234                 if(my_iostat /= 0) write(*,fmt='(a)',advance='no'), '.' !!if(my_iostat /= 0) write(*,fmt='(a)',advance='no'), '.' !!if(my_iostat /= 0) print '(i5,a1,i4,2x,a)',9234,layoutnumber,trim(adjustl(nEntradaRoot))//'_Outputrequests_'//trim(adjustl(whoamishort))//'.txt'
                     write(19,'(a)',err=9234,iostat=my_iostat) trim(adjustl(output(ii)%item(i)%path))
                     !erase pre-existing data unless this is a resuming simulation

                     if (.not.resume) then
                        if (SGG%Observation(ii)%TimeDomain) then
                           open ( output(ii)%item(i)%unit,file= trim(adjustl(output(ii)%item(i)%path)),form='unformatted')
                           write(output(ii)%item(i)%unit) '!END'
                           close (output(ii)%item(i)%unit,status='DELETE')
                           my_iostat=0
9271                       if(my_iostat /= 0) write(*,fmt='(a)',advance='no'), '.' !!if(my_iostat /= 0) print '(i5,a1,i4,2x,a)',9271,'.',layoutnumber,trim(adjustl(output(ii)%item(i)%path)) 
                           open ( output(ii)%item(i)%unit,file= trim(adjustl(output(ii)%item(i)%path)),form='unformatted',err=9271,iostat=my_iostat,status='new',action='write')
                           write(output(ii)%item(i)%unit) &
                                output(ii)%item(i)%XItrancos,output(ii)%item(i)%XEtrancos, &
                                output(ii)%item(i)%YItrancos,output(ii)%item(i)%YEtrancos, &
                                output(ii)%item(i)%ZItrancos,output(ii)%item(i)%ZEtrancos
                          !!!&      sgg%observation(ii)%P(i)%xI,sgg%observation(ii)%P(i)%xE, &
                          !!!&      sgg%observation(ii)%P(i)%YI,sgg%observation(ii)%P(i)%YE, &
                          !!!&      sgg%observation(ii)%P(i)%zI,sgg%observation(ii)%P(i)%ZE
                           !wipe out duplicate data after non synchronous data and field resuming !later

                        elseif (SGG%Observation(ii)%FreqDomain) then
                           open ( output(ii)%item(i)%unit,file= trim(adjustl(output(ii)%item(i)%path)),form='unformatted')
                           write(output(ii)%item(i)%unit) '!END'
                           close (output(ii)%item(i)%unit,status='DELETE')
                        endif !no need to keep it open
                     else
                        if (SGG%Observation(ii)%TimeDomain) then
                           inquire(file= trim(adjustl(output(ii)%item(i)%path)),exist=existe)
                           if (.not.existe) then
                              call stoponerror(layoutnumber,size,'Data files for resuming non existent (volume xdmf...) '//trim(adjustl(output(ii)%item(i)%path)))
                           endif
                           open (output(ii)%item(i)%unit,access='sequential',file= trim(adjustl(output(ii)%item(i)%path)),  &
                           form='unformatted')
                           read(output(ii)%item(i)%unit) ndum,ndum,ndum,ndum,ndum,ndum
                           cutting4 :  do
                              read(output(ii)%item(i)%unit,end=6999) at 
                              do k1=output(ii)%item(i)%ZItrancos,output(ii)%item(i)%ZEtrancos !sgg%Observation(ii)%P(i)%ZI,sgg%Observation(ii)%P(i)%ZE
                                 do j1=output(ii)%item(i)%YItrancos,output(ii)%item(i)%YEtrancos !sgg%Observation(ii)%P(i)%YI,sgg%Observation(ii)%P(i)%YE
                                    read(output(ii)%item(i)%unit,end=6999) (rdum, &
                                    &            i1=output(ii)%item(i)%XItrancos,output(ii)%item(i)%XEtrancos) ! sgg%Observation(ii)%P(i)%XI,sgg%Observation(ii)%P(i)%XE)
                                 end do
                              end do
                              output(ii)%TimesWritten=output(ii)%TimesWritten+1
                              if (at > lastexecutedtime) then
                                 print '(i4,a,a,2e19.9e3)',quienmpi,'Cutting 4 ',trim(adjustl(output(ii)%item(i)%path)),at,lastexecutedtime
                                 backspace(output(ii)%item(i)%unit)
                                 endfile(output(ii)%item(i)%unit)
                                 exit cutting4
                              endif
                           end do cutting4
6999                       continue
                           !!!    backspace(output(ii)%item(i)%unit) !machaco el timeswritten proque solo puede haber uno al final
                           close (output(ii)%item(i)%unit)
                           open (output(ii)%item(i)%unit,file= trim(adjustl(output(ii)%item(i)%path)),position='append',  &
                           form='unformatted')
                        elseif (SGG%Observation(ii)%FreqDomain) then
                           open ( output(ii)%item(i)%unit,file= trim(adjustl(output(ii)%item(i)%path)),form='unformatted')
                           read(output(ii)%item(i)%unit) ndum,ndum,ndum,ndum,ndum,ndum
                           read(output(ii)%item(i)%unit) at
                           if (nint(at/sgg%dt) /= (initialtimestep-1)) then  !estaba mal? ponia initialstep) sin el -1 !261119. lo he cambiado
                              write (buff,*) nint(at/sgg%dt), initialtimestep-1,' Data files for resuming 3D freq domain probes might be corrupt. Continuing....'
                              call print11(layoutnumber,buff)
                           endif

                           DO N=1,output(ii)%NumFreqs
                              read(output(ii)%item(i)%unit) rdum
                              do compo=1,3
                              do k1t=output(ii)%item(i)%ZItrancos,output(ii)%item(i)%ZEtrancos
                                 do j1t=output(ii)%item(i)%YItrancos,output(ii)%item(i)%YEtrancos
                                                    read(output(ii)%item(i)%unit) (output(ii)%item(i)%valor3DComplex(N,compo,i1t,j1t,k1t), &
                                    &                                i1t=output(ii)%item(i)%XItrancos,output(ii)%item(i)%XEtrancos)
                                        end do
                                 end do
                              end do
                              if (sgg%Observation(ii)%Transfer) output(ii)%item(i)%valor3DComplex = output(ii)%item(i)%valor3DComplex * output(ii)%dftEntrada(n) !desnormaliza
                           end do
                           close (output(ii)%item(i)%unit,status='delete')
                        endif
                     endif
                  endif
#ifdef CompileWithNF2FF
                case (farfield)
                  ThereAreFarFields=.true.
                  !
                  write(chari ,'(i7)') sgg%observation(ii)%P(1)%XI
                  write(charj ,'(i7)') sgg%observation(ii)%P(1)%YI
                  write(chark ,'(i7)') sgg%observation(ii)%P(1)%ZI
                  write(chari2,'(i7)') sgg%observation(ii)%P(1)%XE
                  write(charj2,'(i7)') sgg%observation(ii)%P(1)%YE
                  write(chark2,'(i7)') sgg%observation(ii)%P(1)%ZE
                  !mpidir 190319      !desrotacion para que los nombres sean correctos 
                  if      (mpidir==3) then
                      extpoint=trim(adjustl(chari)) //'_'//trim(adjustl(charj)) //'_'//trim(adjustl(chark))//'__'// &
                               trim(adjustl(chari2))//'_'//trim(adjustl(charj2))//'_'//trim(adjustl(chark2))
                      prefix_field=prefix(field)
                  elseif  (mpidir==2) then
                      extpoint=trim(adjustl(charj)) //'_'//trim(adjustl(chark)) //'_'//trim(adjustl(chari))//'__'// &
                               trim(adjustl(charj2))//'_'//trim(adjustl(chark2))//'_'//trim(adjustl(chari2))
                      prefix_field=prefix(field)
                  elseif  (mpidir==1) then
                      extpoint=trim(adjustl(chark)) //'_'//trim(adjustl(chari)) //'_'//trim(adjustl(charj))//'__'// &
                               trim(adjustl(chark2))//'_'//trim(adjustl(chari2))//'_'//trim(adjustl(charj2))
                      prefix_field=prefix(field)
                  else
                      call stoponerror(layoutnumber,size,'Buggy error in mpidir. ')
                  endif
                  !
                  !
                  ext=trim(adjustl(nEntradaRoot))//'_'//trim(adjustl(sgg%observation(ii)%outputrequest))
                  output(ii)%item(i)%path=trim(adjustl(ext))//'_'//trim(adjustl(prefix_field))//&
                  trim(adjustl(extpoint))//'.dat'

                  output(ii)%item(i)%columnas=1
                  !
                  unit=unit+1
                  if (unit >= 2.0_RKIND  **31.0_RKIND-1.0_RKIND) then
                     call stoponerror(layoutnumber,size,'Excesive number of probes')
                  endif
                  output(ii)%item(i)%unit=unit
                  !!!busca nombres de ficheros por duplicado y resuelve la duplicidad
                  call checkduplicatenames
                  !!!!!!
                  !inicializacion especifica del farfield
                  call InitFarField(sgg,sggMiEx,sggMiEy,sggMiEz,sggMiHx,sggMiHy,sggMiHz,layoutnumber,size, &
                  b,resume, &
                  output(ii)%item(i)%unit        ,   &
                  output(ii)%item(i)%path          , &
                  sgg%observation(ii)%P(1)%XI,       &
                  sgg%observation(ii)%P(1)%XE,       &
                  sgg%observation(ii)%P(1)%YI,       &
                  sgg%observation(ii)%P(1)%YE,       &
                  sgg%observation(ii)%P(1)%ZI,       &
                  sgg%observation(ii)%P(1)%ZE,       &
                  sgg%observation(ii)%InitialFreq  , &
                  sgg%observation(ii)%FinalFreq    , &
                  sgg%observation(ii)%FreqStep     , &
                  sgg%observation(ii)%phiStart     , &
                  sgg%observation(ii)%phiStop      , &
                  sgg%observation(ii)%phiStep      , &
                  sgg%observation(ii)%thetaStart   , &
                  sgg%observation(ii)%thetaStop    , &
                  sgg%observation(ii)%thetaStep    , &
                  sgg%observation(ii)%FileNormalize, SINPML_fullsize,facesNF2FF,NF2FFDecim &
#ifdef CompileWithMPI
                  ,output(ii)%item(i)%MPISubComm,output(ii)%item(i)%MPIRoot &
#endif
                  ,eps0,mu0)
#endif
                  !no es necesario hacer wipe out pq en DF se van machacando
               end select
            end do loop_ob
            !!!!        endif !del time domain !NO ES PRECISO 25/02/14
         end do !del ii=1,numberrequest

         write(19,'(a)') '!END '
         close (19)

      endif


      return

   contains

      subroutine checkduplicatenames
         integer(Kind=4) :: n_ii, n_i,off
         p1=output(ii)%item(i)%path
         do n_ii=1,ii
            off=sgg%Observation(n_ii)%nP
            if (n_ii==ii) off=i-1
            do n_i=1,off
               if ( sgg%Observation(n_ii)%P(n_i)%What /= nothing) then
                  p2=output(n_ii)%item(n_i)%path
                  if (trim(adjustl(p1))==trim(adjustl(p2))) then
                     write(charNO,'(i7)') output(ii)%item(i)%unit
                     output(ii)%item(i)%path=trim(adjustl(output(ii)%item(i)%path))//'_duplicate_'//trim(adjustl(charNO))//'.dat'
                  endif
               endif
            end do
         end do
         return
      end subroutine

      !*************************

      subroutine crea_gnuplot

         buff2=trim(adjustl(nEntradaRoot))//'_gnuplot.pl'
         thefile=openfile_mpi (layoutnumber,buff2)

         !!!!!
         conta=0
         do ii=1,sgg%NumberRequest
            do i=1,sgg%Observation(ii)%nP
               I1=sgg%observation(ii)%P(i)%XI
               J1=sgg%observation(ii)%P(i)%YI
               K1=sgg%observation(ii)%P(i)%ZI
               I2=sgg%observation(ii)%P(i)%XE
               J2=sgg%observation(ii)%P(i)%YE
               K2=sgg%observation(ii)%P(i)%ZE
               NO=sgg%observation(ii)%P(i)%NODE
               write(chari,'(i7)') i1
               write(charj,'(i7)') j1
               write(chark,'(i7)') k1
               field=sgg%observation(ii)%P(i)%What
               select case(field)
                case (iEx,iEy,iEz,iVx,iVy,iVz,iJx,iJy,iJz,iHx,iHy,iHz)
                  conta=conta+1                  
                  !mpidir 190319 !desrotacion para que los nombres sean correctos 
                  if      (mpidir==3)then
                      extpoint=trim(adjustl(chari)) //'_'//trim(adjustl(charj)) //'_'//trim(adjustl(chark))
                      select case (field)
                      case (iEx)
                        prefix_field=prefix(iEx)
                      case (iEy)
                        prefix_field=prefix(iEy)
                      case (iEz)
                        prefix_field=prefix(iEz)
                      case (iJx)
                        prefix_field=prefix(iJx)
                      case (iJy)
                        prefix_field=prefix(iJy)
                      case (iJz)
                        prefix_field=prefix(iJz)
                      case (iVx)
                        prefix_field=prefix(iVx)
                      case (iVy)
                        prefix_field=prefix(iVy)
                      case (iVz)
                        prefix_field=prefix(iVz)
                      case (iHx)
                        prefix_field=prefix(iHx)
                      case (iHy)
                        prefix_field=prefix(iHy)
                      case (iHz)
                        prefix_field=prefix(iHz)
                      case default
                        prefix_field=prefix(field)
                      end select
                  elseif  (mpidir==2) then
                      extpoint=trim(adjustl(charj)) //'_'//trim(adjustl(chark)) //'_'//trim(adjustl(chari))
                      select case (field)
                      case (iEx)
                        prefix_field=prefix(iEz)
                      case (iEy)
                        prefix_field=prefix(iEx)
                      case (iEz)
                        prefix_field=prefix(iEy)
                      case (iJx)
                        prefix_field=prefix(iJz)
                      case (iJy)
                        prefix_field=prefix(iJx)
                      case (iJz)
                        prefix_field=prefix(iJy)
                      case (iVx)
                        prefix_field=prefix(iVz)
                      case (iVy)
                        prefix_field=prefix(iVx)
                      case (iVz)
                        prefix_field=prefix(iVy)
                      case (iHx)
                        prefix_field=prefix(iHz)
                      case (iHy)
                        prefix_field=prefix(iHx)
                      case (iHz)
                        prefix_field=prefix(iHy)
                      case default
                        prefix_field=prefix(field)
                      end select
                  elseif  (mpidir==1) then
                      extpoint=trim(adjustl(chark)) //'_'//trim(adjustl(chari)) //'_'//trim(adjustl(charj))
                      select case (field)
                      case (iEx)
                        prefix_field=prefix(iEy)
                      case (iEy)
                        prefix_field=prefix(iEz)
                      case (iEz)
                        prefix_field=prefix(iEx)
                      case (iJx)
                        prefix_field=prefix(iJy)
                      case (iJy)
                        prefix_field=prefix(iJz)
                      case (iJz)
                        prefix_field=prefix(iJx)
                      case (iVx)
                        prefix_field=prefix(iVy)
                      case (iVy)
                        prefix_field=prefix(iVz)
                      case (iVz)
                        prefix_field=prefix(iVx)
                      case (iHx)
                        prefix_field=prefix(iHy)
                      case (iHy)
                        prefix_field=prefix(iHz)
                      case (iHz)
                        prefix_field=prefix(iHx)
                      case default
                        prefix_field=prefix(field)
                      end select
                  else
                      call stoponerror(layoutnumber,size,'Buggy error in mpidir. ')
                  endif
                  !                
                  if ((field == iJx).or.(field == iJy).or.(field == iJz)) then
                     write(charNO,'(i7)') NO
                     !append the number of the segment
                     extpoint=trim(adjustl(extpoint))//'_s'//trim(adjustl(charNO))
                  endif
                  ext=trim(adjustl(nEntradaRoot))//'_'//trim(adjustl(sgg%observation(ii)%outputrequest))
                  !do not use layername since no two observations from different layers will overlap
                  path=CHAR(39)//trim(adjustl(ext))//'_'//trim(adjustl(prefix_field))//trim(adjustl(extpoint))//'.dat'//CHAR(39)
                  write(buff2,*) 'set term x11 persist',conta+10000*layoutnumber
                  call writefile_mpi(layoutnumber,thefile, buff2)
                  write(buff2,*) 'plot ',trim(adjustl(path)),' using 1:2 every 1::2 with lines '
                  call writefile_mpi(layoutnumber,thefile, buff2)
                case (iBloqueJx,iBloqueJy,iBloqueJz,iBloqueMx,iBloqueMy,iBloqueMz)
                  conta=conta+1
                  write(chari2,'(i7)') i2
                  write(charj2,'(i7)') j2
                  write(chark2,'(i7)') k2
                  !mpidir 190319    !desrotacion para que los nombres sean correctos 
                  if      (mpidir==3) then
                      extpoint=trim(adjustl(chari)) //'_'//trim(adjustl(charj)) //'_'//trim(adjustl(chark))//'__'// &
                               trim(adjustl(chari2))//'_'//trim(adjustl(charj2))//'_'//trim(adjustl(chark2))
                      select case (field)
                      case (iBloqueJX)
                        prefix_field=prefix(iBloqueJX)
                      case (iBloqueJY)
                        prefix_field=prefix(iBloqueJY)
                      case (iBloqueJZ)
                        prefix_field=prefix(iBloqueJZ)
                      case (iBloqueMX)
                        prefix_field=prefix(iBloqueMX)
                      case (iBloqueMY)
                        prefix_field=prefix(iBloqueMY)
                      case (iBloqueMZ)
                        prefix_field=prefix(iBloqueMZ)
                      case default
                        prefix_field=prefix(field)
                      end select
                  elseif  (mpidir==2) then
                      extpoint=trim(adjustl(charj)) //'_'//trim(adjustl(chark)) //'_'//trim(adjustl(chari))//'__'// &
                               trim(adjustl(charj2))//'_'//trim(adjustl(chark2))//'_'//trim(adjustl(chari2))
                      select case (field)
                      case (iBloqueJX)
                        prefix_field=prefix(iBloqueJZ)
                      case (iBloqueJY)
                        prefix_field=prefix(iBloqueJX)
                      case (iBloqueJZ)
                        prefix_field=prefix(iBloqueJY)
                      case (iBloqueMX)
                        prefix_field=prefix(iBloqueMZ)
                      case (iBloqueMY)
                        prefix_field=prefix(iBloqueMX)
                      case (iBloqueMZ)
                        prefix_field=prefix(iBloqueMY)
                      case default
                        prefix_field=prefix(field)
                      end select
                  elseif  (mpidir==1) then
                      extpoint=trim(adjustl(chark)) //'_'//trim(adjustl(chari)) //'_'//trim(adjustl(charj))//'__'// &
                               trim(adjustl(chark2))//'_'//trim(adjustl(chari2))//'_'//trim(adjustl(charj2))
                      select case (field)
                      case (iBloqueJX)
                        prefix_field=prefix(iBloqueJY)
                      case (iBloqueJY)
                        prefix_field=prefix(iBloqueJZ)
                      case (iBloqueJZ)
                        prefix_field=prefix(iBloqueJX)
                      case (iBloqueMX)
                        prefix_field=prefix(iBloqueMY)
                      case (iBloqueMY)
                        prefix_field=prefix(iBloqueMZ)
                      case (iBloqueMZ)
                        prefix_field=prefix(iBloqueMX)
                      case default
                        prefix_field=prefix(field)
                      end select
                  else
                      call stoponerror(layoutnumber,size,'Buggy error in mpidir. ')
                  endif
                     !
                  !
                  ext=trim(adjustl(nEntradaRoot))//'_'//trim(adjustl(sgg%observation(ii)%outputrequest))
                  !
                  path=CHAR(39)//trim(adjustl(ext))//'_'//trim(adjustl(prefix_field))//trim(adjustl(extpoint))//'.dat'//CHAR(39)
                  write(buff2,*) 'set term x11 persist',conta+10000*layoutnumber
                  call writefile_mpi(layoutnumber,thefile, buff2)
                  write(buff2,*) 'plot ',trim(adjustl(path)),' using 1:2 every 1::2 with lines '
                  call writefile_mpi(layoutnumber,thefile, buff2)
               end select
            end do
         end do !del ii=1,numberrequest

         buff2=trim(adjustl(nEntradaRoot))//'_gnuplot.pl'
         call closefile_mpi(layoutnumber,size,buff2,thefile)

         return
      end subroutine crea_gnuplot

      
      
   end subroutine InitObservation


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Closes observation stuff
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine CloseObservationFiles(sgg,layoutnumber,size,singlefilewrite,initialtimestep,lastexecutedtime,resume)
      type (SGGFDTDINFO), intent(IN)         ::  sgg
      integer (kind=4)  ::  i,ii,layoutnumber,field,initialtimestep,unidad,size,idum
      logical :: singlefilewrite,resume,incident,existe,wrotemaster
      REAL (KIND=RKIND)    ::  rdum1,rdum2,rdum3,rdum4,rdum5,rdum6,rdum
      REAL (KIND=RKIND_tiempo)    ::  lastexecutedtime
      character (LEN=BUFSIZE) :: chdum
      character (LEN=BUFSIZE)  ::  whoamishort
      integer :: my_iostat
      real (kind=RKIND_tiempo) :: at
      !!!

      write(whoamishort,'(i5)') layoutnumber+1
      !
      if (sgg%NumPlaneWaves>=1) incident=.true.
      do ii=1,sgg%NumberRequest
         if (SGG%Observation(ii)%TimeDomain) then
            if (.not. SGG%Observation(ii)%Volumic) then
               wrotemaster=.false.
               do i=1,sgg%Observation(ii)%nP
                  if ((sgg%observation(ii)%P(i)%What /=nothing).and. &
                  (SGG%Observation(ii)%P(1)%what /= farfield)) then !el farfield se cierra y abre de su propio modo
                     field=sgg%observation(ii)%P(i)%what
                     if (singlefilewrite.and.((field==iEx).or.(field==iEy).or.(field==iEz).or. &
                     (field==iVx).or.(field==iVy).or.(field==iVz).or. &
                     (field==iJx).or.(field==iJy).or.(field==iJz).or. &
                     (field==iHx).or.(field==iHy).or.(field==iHz))) then
                        if (.not.wrotemaster) then
                           wrotemaster=.true.
                           close (output(ii)%item(i)%unitmaster)
                           open  (output(ii)%item(i)%unitmaster-1,file= trim(adjustl(output(ii)%item(i)%path))//'_'//trim(adjustl(whoamishort))//'_master.bin',form='unformatted')
                        else
                           rewind (output(ii)%item(i)%unitmaster-1)
                        endif
                        if (.not.resume) then
                           unidad=output(ii)%item(i)%unit
                           open (unidad,recl=1000,file= trim(adjustl(output(ii)%item(i)%path)))
                           write(unidad,*) '!END'
                           close(unidad,status='delete')
                           my_iostat=0
9242                       if(my_iostat /= 0) write(*,fmt='(a)',advance='no'), '.' !!if(my_iostat /= 0) print '(i5,a1,i4,2x,a)',9242,'.',layoutnumber,trim(adjustl(output(ii)%item(i)%path)) 
                           open (unidad,recl=1000,file= trim(adjustl(output(ii)%item(i)%path)),err=9242,iostat=my_iostat,status='new',action='write')
                           write(unidad,'(a)') trim(adjustl( ' t'//'              '// &
                           trim(adjustl(output(ii)%item(i)%path))//'       '//trim(adjustl(suffix(field,incident))) ))
                        else
                           inquire(file= trim(adjustl(output(ii)%item(i)%path)),exist=existe)
                           if (.not.existe) then
                              call stoponerror(layoutnumber,size,'Data files for resuming non existent (generic closing) '//trim(adjustl(output(ii)%item(i)%path)))
                           endif
                           !
                           open (output(ii)%item(i)%unit,recl=1000,access='sequential',file= trim(adjustl(output(ii)%item(i)%path)))
                           read(output(ii)%item(i)%unit,'(a)') chdum !first line contains characters
                           cutting :  do
                              read(output(ii)%item(i)%unit,*,end=678) at 
                              if (at > lastexecutedtime) then
                                 print '(i4,a,a,2e19.9e3)',quienmpi,'Cutting 5 ',trim(adjustl(output(ii)%item(i)%path)),at,lastexecutedtime
                                 backspace(output(ii)%item(i)%unit)
                                 endfile(output(ii)%item(i)%unit)
                                 exit cutting
                              endif
                           end do cutting
678                        continue
                           close (output(ii)%item(i)%unit)
                           !
                           open (output(ii)%item(i)%unit,recl=1000,file= trim(adjustl(output(ii)%item(i)%path)),position='append')
                        endif
                        !
                        do
                           select case(field)
                            case (iHx,iHy,iHz,iEx,iEy,iEz)
                              if (incident) then
                                 read  (output(ii)%item(i)%unitmaster-1,end=777) idum,rdum1,rdum2,rdum3
                                 if (idum == output(ii)%item(i)%unit) write (output(ii)%item(i)%unit,fmt) rdum1,rdum2,rdum3
                              else
                                 read  (output(ii)%item(i)%unitmaster-1,end=777) idum,rdum1,rdum2
                                 if (idum == output(ii)%item(i)%unit) write (output(ii)%item(i)%unit,fmt) rdum1,rdum2
                              endif
                            case(iJx,iJy,iJz)
                              read  (output(ii)%item(i)%unitmaster-1,end=777) idum,rdum1,rdum2,rdum3,rdum4,rdum5,rdum6
                              if (idum == output(ii)%item(i)%unit) write (output(ii)%item(i)%unit,fmt) rdum1,rdum2,rdum3,rdum4,rdum5,rdum6
                            case default
                              read  (output(ii)%item(i)%unitmaster-1,end=777) idum,rdum1,rdum2 !caso de los votalges ivx, etc
                              if (idum == output(ii)%item(i)%unit) write (output(ii)%item(i)%unit,fmt) rdum1,rdum2
                           end select
                        end do
777                     close (output(ii)%item(i)%unit)
                     else
                        close (output(ii)%item(i)%unit)
                     endif
                  endif
               end do
            else !sondas volumicas
               do i=1,sgg%Observation(ii)%nP
                  if ((sgg%observation(ii)%P(i)%What /=nothing).and. &
                  (SGG%Observation(ii)%P(1)%what /= farfield)) then
                     !write(output(ii)%item(i)%unit) output(ii)%TimesWritten !el farfield se cierra y abre de su propio modo
                     endfile(output(ii)%item(i)%unit)
                     close (output(ii)%item(i)%unit)
                  endif
               end do
               !!!!
            endif
         elseif (SGG%Observation(ii)%FreqDomain) then
            continue !nothing to do since the freq domain is updated by opening, writing and closing each time !esot hace que no se pueda hacer cutting (si hay fallos y restarteos)
         endif
         wrotemaster=.false.
         do i=1,sgg%Observation(ii)%nP
            if (SGG%Observation(ii)%TimeDomain) then
               if ((sgg%observation(ii)%P(i)%What /=nothing).and. &
               (SGG%Observation(ii)%P(1)%what /= farfield)) then !el farfield se cierra y abre de su propio modo
                  field=sgg%observation(ii)%P(i)%what
                  if (singlefilewrite.and.((field==iEx).or.(field==iEy).or.(field==iEz).or. &
                  (field==iVx).or.(field==iVy).or.(field==iVz).or. &
                  (field==iJx).or.(field==iJy).or.(field==iJz).or. &
                  (field==iHx).or.(field==iHy).or.(field==iHz))) then
                     if (.not.wrotemaster) then
                        wrotemaster=.true.
                        close (output(ii)%item(i)%unitmaster-1,status='delete') !el binario no se precisa para nada
                     endif
                  endif
               endif
            endif !no hay singlewrite para sondas freqdomain
         end do
      end do


      return
   end subroutine CloseObservationFiles


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Upacks .bin files observation stuff
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine UnpackSingleFiles(sgg,layoutnumber,size,singlefilewrite,initialtimestep,resume)
      type (SGGFDTDINFO), intent(IN)         ::  sgg
      integer (kind=4)  ::  i,ii,layoutnumber,field,initialtimestep,unidad,size,idum
      logical :: singlefilewrite,resume,incident,existe,wrotemaster
      REAL (KIND=RKIND)    ::  rdum1,rdum2,rdum3,rdum4,rdum5,rdum6,rdum
      character (LEN=BUFSIZE) :: chdum
      character (LEN=BUFSIZE)  ::  whoamishort
      integer :: my_iostat
      !!!
      write(whoamishort,'(i5)') layoutnumber+1
      !
      if (sgg%NumPlaneWaves>=1) incident=.true.
      do ii=1,sgg%NumberRequest
         if (SGG%Observation(ii)%TimeDomain) then
            if (.not. SGG%Observation(ii)%Volumic) then
               wrotemaster=.false.
               do i=1,sgg%Observation(ii)%nP
                  if ((sgg%observation(ii)%P(i)%What /=nothing).and. &
                  (SGG%Observation(ii)%P(1)%what /= farfield)) then !el farfield se cierra y abre de su propio modo
                     field=sgg%observation(ii)%P(i)%what
                     if (singlefilewrite.and.((field==iEx).or.(field==iEy).or.(field==iEz).or. &
                     (field==iVx).or.(field==iVy).or.(field==iVz).or. &
                     (field==iJx).or.(field==iJy).or.(field==iJz).or. &
                     (field==iHx).or.(field==iHy).or.(field==iHz))) then
                        rewind (output(ii)%item(i)%unitmaster)
                        unidad=35
                        open (unidad,recl=1000,file= trim(adjustl(output(ii)%item(i)%path)))
                        write(unidad,*) '!END'
                        close(unidad,status='delete')
                        my_iostat=0
9243                    if(my_iostat /= 0) write(*,fmt='(a)',advance='no'), '.' !!if(my_iostat /= 0) print '(i5,a1,i4,2x,a)',9243,'.',layoutnumber,trim(adjustl(output(ii)%item(i)%path)) 
                        open (unidad,recl=1000,file= trim(adjustl(output(ii)%item(i)%path)),err=9243,iostat=my_iostat,status='new',action='write')
                        write(unidad,'(a)') trim(adjustl(' t'//'              '// &
                        trim(adjustl(output(ii)%item(i)%path))//'       '//trim(adjustl(suffix(field,incident))) ))
                        !
                        do
                           select case(field)
                            case (iHx,iHy,iHz,iEx,iEy,iEz)
                              if (incident) then
                                 read  (output(ii)%item(i)%unitmaster,end=7778) idum,rdum1,rdum2,rdum3
                                 if (idum == output(ii)%item(i)%unit) write (unidad,fmt) rdum1,rdum2,rdum3
                              else
                                 read  (output(ii)%item(i)%unitmaster,end=7778) idum,rdum1,rdum2
                                 if (idum == output(ii)%item(i)%unit) write (unidad,fmt) rdum1,rdum2
                              endif
                            case(iJx,iJy,iJz)
                              read  (output(ii)%item(i)%unitmaster,end=7778) idum,rdum1,rdum2,rdum3,rdum4,rdum5,rdum6
                              if (idum == output(ii)%item(i)%unit) write (unidad,fmt) rdum1,rdum2,rdum3,rdum4,rdum5,rdum6
                            case default
                              read  (output(ii)%item(i)%unitmaster,end=7778) idum,rdum1,rdum2 !caso de los votalges ivx, etc
                              if (idum == output(ii)%item(i)%unit) write (unidad,fmt) rdum1,rdum2
                           end select
                        end do
7778                    close (unidad)
                     endif
                  endif
               end do
            endif
         endif
      end do


      return
   end subroutine UnpackSingleFiles



   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Updates the observed values. A nodal average is used for each field
   !!! The Wire modules uses its own updating procedure
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine UpdateObservation(sgg,sggMiEx,sggMiEy,sggMiEz,sggMiHx,sggMiHy,sggMiHz,sggMtag, &
      nTime,nInit, b, Ex, Ey, Ez, Hx, Hy, Hz, dxe, dye, dze, dxh, dyh, dzh,wiresflavor,SINPML_fullsize,wirecrank, &
       Exvac, Eyvac, Ezvac, Hxvac, Hyvac, Hzvac,Excor, Eycor, Ezcor, Hxcor, Hycor, Hzcor,planewavecorr,noconformalmapvtk)
      logical :: noconformalmapvtk
      type (SGGFDTDINFO), intent(IN)         ::  sgg
      INTEGER (KIND=IKINDMTAG), intent(in) :: sggMtag  (sgg%Alloc(iHx)%XI:sgg%Alloc(iHx)%XE, sgg%Alloc(iHy)%YI:sgg%Alloc(iHy)%YE, sgg%Alloc(iHz)%ZI:sgg%Alloc(iHz)%ZE)
      !---------------------------> inputs <----------------------------------------------------------
      logical :: planewavecorr
      type (limit_t), dimension(1:6), intent(in)  ::  SINPML_fullsize
      integer, intent( IN)  ::  nTime,nInit
      type( bounds_t), intent( IN)  ::  b
      !--->
      real (kind = RKIND), dimension( 0 :  b%Ex%NX-1, 0 :  b%Ex%NY-1, 0 :  b%Ex%NZ-1), intent( IN)  ::  Ex
      real (kind = RKIND), dimension( 0 :  b%Ey%NX-1, 0 :  b%Ey%NY-1, 0 :  b%Ey%NZ-1), intent( IN)  ::  Ey
      real (kind = RKIND), dimension( 0 :  b%Ez%NX-1, 0 :  b%Ez%NY-1, 0 :  b%Ez%NZ-1), intent( IN)  ::  Ez
      !--->
      real (kind = RKIND), dimension( 0 :  b%Hx%NX-1, 0 :  b%Hx%NY-1, 0 :  b%Hx%NZ-1), intent( IN)  ::  Hx
      real (kind = RKIND), dimension( 0 :  b%Hy%NX-1, 0 :  b%Hy%NY-1, 0 :  b%Hy%NZ-1), intent( IN)  ::  Hy
      real (kind = RKIND), dimension( 0 :  b%Hz%NX-1, 0 :  b%Hz%NY-1, 0 :  b%Hz%NZ-1), intent( IN)  ::  Hz
      !--->
      !--->
      real (kind = RKIND), dimension( 0 :  b%Ex%NX-1, 0 :  b%Ex%NY-1, 0 :  b%Ex%NZ-1), intent( IN)  ::  Exvac
      real (kind = RKIND), dimension( 0 :  b%Ey%NX-1, 0 :  b%Ey%NY-1, 0 :  b%Ey%NZ-1), intent( IN)  ::  Eyvac
      real (kind = RKIND), dimension( 0 :  b%Ez%NX-1, 0 :  b%Ez%NY-1, 0 :  b%Ez%NZ-1), intent( IN)  ::  Ezvac
      !--->                                                                                               
      real (kind = RKIND), dimension( 0 :  b%Hx%NX-1, 0 :  b%Hx%NY-1, 0 :  b%Hx%NZ-1), intent( IN)  ::  Hxvac
      real (kind = RKIND), dimension( 0 :  b%Hy%NX-1, 0 :  b%Hy%NY-1, 0 :  b%Hy%NZ-1), intent( IN)  ::  Hyvac
      real (kind = RKIND), dimension( 0 :  b%Hz%NX-1, 0 :  b%Hz%NY-1, 0 :  b%Hz%NZ-1), intent( IN)  ::  Hzvac
      !--->
      !--->
      real (kind = RKIND), dimension( 0 :  b%Ex%NX-1, 0 :  b%Ex%NY-1, 0 :  b%Ex%NZ-1), intent( INout)  ::  Excor
      real (kind = RKIND), dimension( 0 :  b%Ey%NX-1, 0 :  b%Ey%NY-1, 0 :  b%Ey%NZ-1), intent( INout)  ::  Eycor
      real (kind = RKIND), dimension( 0 :  b%Ez%NX-1, 0 :  b%Ez%NY-1, 0 :  b%Ez%NZ-1), intent( INout)  ::  Ezcor
      !--->                                                                                               
      real (kind = RKIND), dimension( 0 :  b%Hx%NX-1, 0 :  b%Hx%NY-1, 0 :  b%Hx%NZ-1), intent( INout)  ::  Hxcor
      real (kind = RKIND), dimension( 0 :  b%Hy%NX-1, 0 :  b%Hy%NY-1, 0 :  b%Hy%NZ-1), intent( INout)  ::  Hycor
      real (kind = RKIND), dimension( 0 :  b%Hz%NX-1, 0 :  b%Hz%NY-1, 0 :  b%Hz%NZ-1), intent( INout)  ::  Hzcor
      !--->
      
      real (kind = RKIND), dimension( 0 :  b%dxe%NX-1), intent( IN)  ::  dxe
      real (kind = RKIND), dimension( 0 :  b%dye%NY-1), intent( IN)  ::  dye
      real (kind = RKIND), dimension( 0 :  b%dze%NZ-1), intent( IN)  ::  dze
      !--->
      real (kind = RKIND), dimension( 0 :  b%dxh%NX-1), intent( IN)  ::  dxh
      real (kind = RKIND), dimension( 0 :  b%dyh%NY-1), intent( IN)  ::  dyh
      real (kind = RKIND), dimension( 0 :  b%dzh%NZ-1), intent( IN)  ::  dzh
      !!!
      !----->
      integer(kind = INTEGERSIZEOFMEDIAMATRICES), dimension ( 0 : b%sggMiEx%NX-1 , 0 : b%sggMiEx%NY-1 , 0 : b%sggMiEx%NZ-1 )  , intent( IN)     ::  sggMiEx
      integer(kind = INTEGERSIZEOFMEDIAMATRICES), dimension ( 0 : b%sggMiEy%NX-1 , 0 : b%sggMiEy%NY-1 , 0 : b%sggMiEy%NZ-1 )  , intent( IN)     ::  sggMiEy
      integer(kind = INTEGERSIZEOFMEDIAMATRICES), dimension ( 0 : b%sggMiEz%NX-1 , 0 : b%sggMiEz%NY-1 , 0 : b%sggMiEz%NZ-1 )  , intent( IN)     ::  sggMiEz
      integer(kind = INTEGERSIZEOFMEDIAMATRICES), dimension ( 0 : b%sggMiHx%NX-1 , 0 : b%sggMiHx%NY-1 , 0 : b%sggMiHx%NZ-1 )  , intent( IN)     ::  sggMiHx
      integer(kind = INTEGERSIZEOFMEDIAMATRICES), dimension ( 0 : b%sggMiHy%NX-1 , 0 : b%sggMiHy%NY-1 , 0 : b%sggMiHy%NZ-1 )  , intent( IN)     ::  sggMiHy
      integer(kind = INTEGERSIZEOFMEDIAMATRICES), dimension ( 0 : b%sggMiHz%NX-1 , 0 : b%sggMiHz%NY-1 , 0 : b%sggMiHz%NZ-1 )  , intent( IN)     ::  sggMiHz

      !---------------------------> variables locales <-----------------------------------------------
      integer( kind = 4)  ::  i, ii, i1, i2, j1, j2, k1, k2, i1_m, i2_m, j1_m, j2_m, k1_m, k2_m, field,jjx,jjy,jjz,if1,i1t,j1t,k1t,iff1
      integer( kind = 4)  ::  iii, kkk, jjj, jjj_m, iii_m, kkk_m,NtimeforVolumic,imed,imed1,imed2,imed3,imed4
      logical :: esborde,wirecrank
      REAL (KIND=RKIND_tiempo)    ::  at
      real (kind = RKIND) :: jx,jy,jz
      integer(kind=4) :: conta !para realmente dar tangenciales de campos en los medios superficiales
      character(len=*), INTENT(in) :: wiresflavor
#ifdef CompileWithWires
      type( CurrentSegments), pointer  ::  segmDumm !segmento de hilo que se observa si lo hubiere
#endif
      !
#ifdef CompileWithBerengerWires      
      type(TSegment)        , pointer  ::  segmDumm_Berenger !segmento de hilo que se observa si lo hubiere
#endif
      !
#ifdef CompileWithSlantedWires      
      class(Segment)        , pointer  ::  segmDumm_Slanted !segmento de hilo que se observa si lo hubiere
#endif

      logical ::  INIT,GEOM,ASIGNA,electric,magnetic
      !!!
      at=-1; jx=-1; jy=-1;jz=-1;  !para que gfortran no me diga que no las inicializo

      !---------------------------> empieza UpdateObservation <---------------------------------------
      do ii = 1, sgg%NumberRequest
         loop_obser :  do i = 1, sgg%Observation( ii)%nP
            field = SGG%Observation( ii)%P( i)%what
            if( field /= nothing) then
               I1 = SGG%Observation( ii)%P( i)%XI
               J1 = SGG%Observation( ii)%P( i)%YI
               K1 = SGG%Observation( ii)%P( i)%ZI
               I2 = SGG%Observation( ii)%P( i)%XE !ojo estos no se usan salvo en Bloques y Volumics
               J2 = SGG%Observation( ii)%P( i)%YE
               K2 = SGG%Observation( ii)%P( i)%ZE
               !--->
               if (SGG%Observation(ii)%TimeDomain) then
                  selectcase( field)
                   case( iEx)
                     output( ii)%item( i)%valor(nTime-nInit) = 0.0_RKIND !wipe value
                     i1_m = I1 - b%Ex%XI
                     j1_m = J1 - b%Ex%YI
                     k1_m = K1 - b%Ex%ZI
                     output( ii)%item( i)%valor(nTime-nInit) = Ex( i1_m, j1_m, k1_m)
                   case( iEy)
                     output( ii)%item( i)%valor(nTime-nInit) = 0.0_RKIND !wipe value
                     i1_m = I1 - b%Ey%XI
                     j1_m = J1 - b%Ey%YI
                     k1_m = K1 - b%Ey%ZI
                     output( ii)%item( i)%valor(nTime-nInit) = Ey( i1_m, j1_m, k1_m)
                   case( iEz)
                     output( ii)%item( i)%valor(nTime-nInit) = 0.0_RKIND !wipe value
                     i1_m = I1 - b%Ez%XI
                     j1_m = J1 - b%Ez%YI
                     k1_m = K1 - b%Ez%ZI
                     output( ii)%item( i)%valor(nTime-nInit) = Ez( i1_m, j1_m, k1_m)
                   case( iHx)
                     output( ii)%item( i)%valor(nTime-nInit) = 0.0_RKIND !wipe value
                     i1_m = I1 - b%Hx%XI
                     j1_m = J1 - b%Hx%YI
                     k1_m = K1 - b%Hx%ZI
                     output( ii)%item( i)%valor(nTime-nInit) = Hx( i1_m, j1_m, k1_m)
                   case( iHy)
                     output( ii)%item( i)%valor(nTime-nInit) = 0.0_RKIND !wipe value
                     i1_m = I1 - b%Hy%XI
                     j1_m = J1 - b%Hy%YI
                     k1_m = K1 - b%Hy%ZI
                     output( ii)%item( i)%valor(nTime-nInit) = Hy( i1_m, j1_m, k1_m)
                   case( iHz)
                     output( ii)%item( i)%valor(nTime-nInit) = 0.0_RKIND !wipe value
                     i1_m = I1 - b%Hz%XI
                     j1_m = J1 - b%Hz%YI
                     k1_m = K1 - b%Hz%ZI
                     output( ii)%item( i)%valor(nTime-nInit) = Hz( i1_m, j1_m, k1_m)
                   case( iBloqueJx)
                     output( ii)%item( i)%valor(nTime-nInit) = 0.0_RKIND !wipe value
                     i1_m = I1 - b%Hy%XI
                     j1_m = J1 - b%Hy%YI
                     k1_m = K1 - b%Hy%ZI
                     k2_m = K2 - b%Hy%ZI
                     do JJJ = j1, j2
                        JJJ_m = JJJ - b%Hy%YI
                        !--->
                        output( ii)%item( i)%valor(nTime-nInit) =   &
                        output( ii)%item( i)%valor(nTime-nInit) +   &
                        (Hy( i1_m, JJJ_m, k1_m-1) - Hy(i1_m, JJJ_m, k2_m)) * dyh( JJJ_m )
                     enddo
                     !--->
                     i1_m = I1 - b%Hz%XI
                     j1_m = J1 - b%Hz%YI
                     j2_m = J2 - b%Hz%YI
                     do KKK = k1, k2
                        KKK_m = KKK - b%Hz%ZI
                        !--->
                        output( ii)%item( i)%valor(nTime-nInit) =  &
                        output( ii)%item( i)%valor(nTime-nInit) +   &
                        (-Hz( i1_m, j1_m-1, KKK_m) + Hz( i1_m, j2_m, KKK_m)) * dzh( KKK_m )
                     enddo
                   case( iBloqueJy)
                     output( ii)%item( i)%valor(nTime-nInit) = 0.0_RKIND !wipe value
                     i1_m = I1 - b%Hz%XI
                     j1_m = J1 - b%Hz%YI
                     i2_m = I2 - b%Hz%XI
                     do KKK = k1,k2
                        KKK_m = KKK - b%Hz%ZI
                        !--->
                        output( ii)%item( i)%valor(nTime-nInit) =  &
                        output( ii)%item( i)%valor(nTime-nInit) +   &
                        (-Hz( i2_m, j1_m, KKK_m) + Hz( i1_m-1, j1_m, KKK_m)) * dzh( KKK_m )
                     enddo
                     !--->
                     j1_m = J1 - b%Hx%YI
                     k1_m = K1 - b%Hx%ZI
                     k2_m = K2 - b%Hx%ZI
                     do III = i1, i2
                        III_m = III - b%Hx%XI
                        !--->
                        output( ii)%item( i)%valor(nTime-nInit) =  &
                        output( ii)%item( i)%valor(nTime-nInit) +   &
                        (Hx( III_m, j1_m, k2_m) - Hx( III_m, j1_m, k1_m-1)) * dxh( III_m )
                     enddo
                   case( iBloqueJz)
                     output( ii)%item( i)%valor(nTime-nInit) = 0.0_RKIND !wipe value
                     j1_m = J1 - b%Hx%YI
                     k1_m = K1 - b%Hx%ZI
                     j2_m = J2 - b%Hx%YI
                     do III = i1, i2
                        III_m = III - b%Hx%XI
                        !--->
                        output( ii)%item( i)%valor(nTime-nInit) =  &
                        output( ii)%item( i)%valor(nTime-nInit) +   &
                        (Hx( III_m, j1_m-1, k1_m) - Hx( III_m, j2_m, k1_m)) * dxh( III_m )
                     enddo
                     !--->
                     i1_m = I1 - b%Hy%XI
                     k1_m = K1 - b%Hy%ZI
                     i2_m = I2 - b%Hy%XI
                     do JJJ = j1, j2
                        JJJ_m = JJJ - b%Hy%YI
                        !--->
                        output( ii)%item( i)%valor(nTime-nInit) =  &
                        output( ii)%item( i)%valor(nTime-nInit) +   &
                        (-Hy( i1_m-1, JJJ_m, k1_m) + Hy( i2_m, JJJ_m, k1_m)) * dyh( JJJ_m )
                     enddo
                   case( iBloqueMx)
                     output( ii)%item( i)%valor(nTime-nInit) = 0.0_RKIND !wipe value
                     i1_m = I1 - b%Ey%XI
                     k1_m = K1 - b%Ey%ZI
                     k2_m = K2 - b%Ey%ZI
                     do JJJ = j1, j2
                        JJJ_m = JJJ - b%Ey%YI
                        !--->
                        output( ii)%item( i)%valor(nTime-nInit) =  &
                        output( ii)%item( i)%valor(nTime-nInit) +   &
                        (-Ey( i1_m, JJJ_m, k1_m) + Ey( i1_m, JJJ_m, k2_m+1)) * dye( JJJ_m )
                     enddo
                     !--->
                     i1_m = I1 - b%Ez%XI
                     j1_m = J1 - b%Ez%YI
                     j2_m = J2 - b%Ez%YI
                     do KKK = k1,k2
                        KKK_m = KKK - b%Ez%ZI
                        !--->
                        output( ii)%item( i)%valor(nTime-nInit) =  &
                        output( ii)%item( i)%valor(nTime-nInit) +   &
                        (Ez( i1_m, j1_m, KKK_m) - Ez( i1_m, j2_m+1, KKK_m)) * dze(KKK_m )
                     enddo
                   case( iBloqueMy)
                     output( ii)%item( i)%valor(nTime-nInit) = 0.0_RKIND !wipe value
                     i1_m = I1 - b%Ez%XI
                     j1_m = J1 - b%Ez%YI
                     i2_m = I2 - b%Ez%XI
                     do KKK = k1, k2
                        KKK_m = KKK - b%Ez%ZI
                        !--->
                        output( ii)%item( i)%valor(nTime-nInit) =  &
                        output( ii)%item( i)%valor(nTime-nInit) +   &
                        (Ez( i2_m+1, j1_m, KKK_m) - Ez( i1_m, j1_m, KKK_m)) * dze( KKK_m )
                     enddo
                     !--->
                     j1_m = J1 - b%Ex%YI
                     k1_m = K1 - b%Ex%ZI
                     k2_m = K2 - b%Ex%ZI
                     do III = i1, i2
                        III_m = III - b%Ex%XI
                        !--->
                        output( ii)%item( i)%valor(nTime-nInit) =  &
                        output( ii)%item( i)%valor(nTime-nInit) +   &
                        (-Ex( III_m, j1_m, k2_m+1) + Ex( III_m, j1_m, k1_m)) * dxe( III_m )
                     enddo
                   case( iBloqueMz)
                     output( ii)%item( i)%valor(nTime-nInit) = 0.0_RKIND !wipe value
                     j1_m = J1 - b%Ex%YI
                     k1_m = K1 - b%Ex%ZI
                     j2_m = J2 - b%Ex%YI
                     do III = i1, i2
                        III_m = III - b%Ex%XI
                        !--->
                        output( ii)%item( i)%valor(nTime-nInit) =  &
                        output( ii)%item( i)%valor(nTime-nInit) +   &
                        (-Ex( III_m, j1_m, k1_m) + Ex( III_m, j2_m+1, k1_m)) * dxe( III_m )
                     enddo
                     !--->
                     i1_m = I1 - b%Ey%XI
                     k1_m = K1 - b%Ey%ZI
                     i2_m = I2 - b%Ey%XI
                     do JJJ = j1, j2
                        JJJ_m = JJJ - b%Ey%YI
                        !--->
                        output( ii)%item( i)%valor(nTime-nInit) =  &
                        output( ii)%item( i)%valor(nTime-nInit) +   &
                        (Ey( i1_m, JJJ_m, k1_m) - Ey( i2_m+1, JJJ_m, k1_m)) * dye( JJJ_m )
                     enddo
                   case( iJx, iJy, iJz)
#ifdef CompileWithWires
                     if ((trim(adjustl(wiresflavor))=='holland') .or. &
                         (trim(adjustl(wiresflavor))=='transition')) then
                        output( ii)%item( i)%valor(nTime-nInit) = 0.0_RKIND !wipe value
                        SegmDumm => output( ii)%item( i)%Segmento
                        if (wirecrank) then !no hay que promediar nada porque estan co-locados en tiempo
                           output( ii)%item(i)%valor(nTime-nInit)= output(ii)%item(i)%valorsigno* &
                                 SegmDumm%Currentpast
                           output( ii)%item(i)%valor2(nTime-nInit)= -SegmDumm%Efield_wire2main * SegmDumm%delta                           
                           output( ii)%item(i)%valor3(nTime-nInit)= output(ii)%item(i)%valorsigno* &
                                 (((SegmDumm%ChargePlus%ChargePresent)))* SegmDumm%Lind * ( InvMu(SegmDumm%indexmed) * InvEps( SegmDumm%indexmed))
                           output( ii)%item(i)%valor4(nTime-nInit)= output(ii)%item(i)%valorsigno*&
                                 (((SegmDumm%ChargeMinus%ChargePresent)))*SegmDumm%Lind * (InvMu( SegmDumm%indexmed) * InvEps( SegmDumm%indexmed))
                           output( ii)%item(i)%valor5(nTime-nInit)=output( ii)%item(i)%valor3(nTime-nInit)-output( ii)%item(i)%valor4(nTime-nInit)
    
                        else
!!saco el potencial calculado con E*delta !051115 !!y aniado el vdrop antinugo porque la Z se hacien bien con este 030719
                           output( ii)%item(i)%valor(nTime-nInit)= output(ii)%item(i)%valorsigno* &
                                 SegmDumm%currentpast
                           output( ii)%item(i)%valor2(nTime-nInit)= -SegmDumm%Efield_wire2main * SegmDumm%delta
                           output( ii)%item(i)%valor3(nTime-nInit)= output(ii)%item(i)%valorsigno*&
                                 (((SegmDumm%ChargePlus%ChargePresent + SegmDumm%ChargePlus%ChargePast))/2.0_RKIND)  * &
                                                                              SegmDumm%Lind * (InvMu(SegmDumm%indexmed) * InvEps( SegmDumm%indexmed)) 
                           output( ii)%item(i)%valor4(nTime-nInit)= output(ii)%item(i)%valorsigno*&
                                 (((SegmDumm%ChargeMinus%ChargePresent + SegmDumm%ChargeMinus%ChargePast))/2.0_RKIND)* &
                                                                              SegmDumm%Lind * (InvMu(SegmDumm%indexmed) * InvEps( SegmDumm%indexmed))
                           output( ii)%item(i)%valor5(nTime-nInit)=output( ii)%item(i)%valor3(nTime-nInit)-output( ii)%item(i)%valor4(nTime-nInit)

                        endif
                        !!!!!!!!!!!!!!!!!!
                     endif   
#endif                    
#ifdef CompileWithBerengerWires
                     if (trim(adjustl(wiresflavor))=='berenger') then 
                        SegmDumm_Berenger => output( ii)%item( i)%Segmento_Berenger
                        output( ii)%item(i)%valor(nTime-nInit)= output(ii)%item(i)%valorsigno* &
                                 SegmDumm_Berenger%Currentpast
                        output( ii)%item(i)%valor2(nTime-nInit)= -SegmDumm_Berenger%field * SegmDumm_Berenger%dl                          
                        output( ii)%item(i)%valor3(nTime-nInit)= output(ii)%item(i)%valorsigno*&
                                  (((SegmDumm_Berenger%ChargePlus  + SegmDumm_Berenger%ChargePlusPast))/2.0_RKIND)  * &
                                   SegmDumm_Berenger%L *    (InvMu(SegmDumm_Berenger%imed) * InvEps( SegmDumm_Berenger%imed))                                                                  
                        output( ii)%item(i)%valor4(nTime-nInit)= output(ii)%item(i)%valorsigno*&
                                 (((SegmDumm_Berenger%ChargeMinus + SegmDumm_Berenger%ChargeMinusPast))/2.0_RKIND) * &
                                   SegmDumm_Berenger%L *    (InvMu(SegmDumm_Berenger%imed) * InvEps( SegmDumm_Berenger%imed))
                        output( ii)%item(i)%valor5(nTime-nInit)=output( ii)%item(i)%valor3(nTime-nInit)-output( ii)%item(i)%valor4(nTime-nInit)
                     endif
#endif
#ifdef CompileWithSlantedWires
                     if((trim(adjustl(wiresflavor))=='slanted').or.(trim(adjustl(wiresflavor))=='semistructured')) then !del wiresflavor
                        SegmDumm_Slanted => output( ii)%item( i)%Segmento_Slanted
                        output( ii)%item(i)%valor(nTime-nInit)= & !ojo: slanted ya los orienta bien y no hay que multiplicar por valorsigno
                                 SegmDumm_Slanted%Currentpast
                        output( ii)%item(i)%valor2(nTime-nInit)= -SegmDumm_Slanted%field * SegmDumm_Slanted%dl                        
                        output( ii)%item(i)%valor3(nTime-nInit)= &
                                 (((SegmDumm_Slanted%Voltage(iPlus)%ptr%Voltage + SegmDumm_Slanted%Voltage(iPlus)%ptr%VoltagePast))/2.0_RKIND)                  
                        output( ii)%item(i)%valor4(nTime-nInit)= &
                                 (((SegmDumm_Slanted%Voltage(iMinus)%ptr%Voltage + SegmDumm_Slanted%Voltage(iMinus)%ptr%VoltagePast))/2.0_RKIND)
                        output( ii)%item(i)%valor5(nTime-nInit)=output( ii)%item(i)%valor3(nTime-nInit)-output( ii)%item(i)%valor4(nTime-nInit)
                     endif
#endif
                     !Volumic probes
                   case( iExC)
                     at=sgg%tiempo(ntime)
                        if (at >  sgg%OBSERVATION(ii)%FinalTime+sgg%dt/2.0_RKIND) sgg%OBSERVATION(ii)%Done=.true.
                        if (at >= sgg%OBSERVATION(ii)%InitialTime)     sgg%OBSERVATION(ii)%Begun=.true.
                     Ntimeforvolumic=Ntime !!!-nint(0.4999999+sgg%OBSERVATION(ii)%InitialTime/sgg%dt)
                     if (mod(Ntimeforvolumic,output(ii)%Trancos)==0) then
                        Ntimeforvolumic=Ntimeforvolumic/output(ii)%Trancos
                        if ( ((at >= sgg%OBSERVATION(ii)%InitialTime).and.(at <= sgg%OBSERVATION(ii)%FinalTime+sgg%dt/2.0_RKIND))) then
                           do KKK = k1, k2
                           if (mod(KKK,output(ii)%item(i)%Ztrancos)==0) then
                           k1t=int(kkk/output(ii)%item(i)%Ztrancos)
                             KKK_m = KKK - b%Ex%ZI
                             do JJJ = j1, j2
                             if (mod(jjj,output(ii)%item(i)%Ytrancos)==0) then
                             j1t=int(jjj/output(ii)%item(i)%Ytrancos)
                                 JJJ_m = JJJ - b%Ex%YI
                                 do III = i1, i2
                                 if (mod(iii,output(ii)%item(i)%Xtrancos)==0) then
                                 i1t=int(iii/output(ii)%item(i)%Xtrancos)
                                    III_m = III - b%Ex%XI
                                    output( ii)%item( i)%valor3D(Ntimeforvolumic,i1t,j1t,k1t) =  Ex( III_m, JJJ_m, KKK_m)
                                 endif
                                 enddo
                              endif
                              enddo
                           endif
                           enddo
                        endif
                     endif
                   case( iEyC)
                     at=sgg%tiempo(ntime)
                        if (at > sgg%OBSERVATION(ii)%FinalTime+sgg%dt/2.0_RKIND) sgg%OBSERVATION(ii)%Done=.true.
                        if (at >= sgg%OBSERVATION(ii)%InitialTime)     sgg%OBSERVATION(ii)%Begun=.true.
                     Ntimeforvolumic=Ntime !!!-nint(0.4999999+sgg%OBSERVATION(ii)%InitialTime/sgg%dt)
                     if (mod(Ntimeforvolumic,output(ii)%Trancos)==0) then
                        Ntimeforvolumic=Ntimeforvolumic/output(ii)%Trancos
                        if ( ((at >= sgg%OBSERVATION(ii)%InitialTime).and.(at <= sgg%OBSERVATION(ii)%FinalTime+sgg%dt/2.0_RKIND))) then
                           do KKK = k1, k2
                           if (mod(KKK,output(ii)%item(i)%Ztrancos)==0) then
                           k1t=int(kkk/output(ii)%item(i)%Ztrancos)
                              KKK_m = KKK - b%Ey%ZI
                              do JJJ = j1, j2
                              if (mod(jjj,output(ii)%item(i)%Ytrancos)==0) then
                              j1t=int(jjj/output(ii)%item(i)%Ytrancos)
                                 JJJ_m = JJJ - b%Ey%YI
                                 do III = i1, i2
                                 if (mod(iii,output(ii)%item(i)%Xtrancos)==0) then
                                 i1t=int(iii/output(ii)%item(i)%Xtrancos)
                                    III_m = III - b%Ey%XI
                                    output( ii)%item( i)%valor3D(Ntimeforvolumic,i1t,j1t,k1t) =  Ey( III_m, JJJ_m, KKK_m)
                                 endif
                                 enddo
                              endif
                              enddo
                           endif
                           enddo
                        endif
                     endif
                   case( iEzC)
                     at=sgg%tiempo(ntime)
                        if (at > sgg%OBSERVATION(ii)%FinalTime+sgg%dt/2.0_RKIND) sgg%OBSERVATION(ii)%Done=.true.
                        if (at >= sgg%OBSERVATION(ii)%InitialTime)     sgg%OBSERVATION(ii)%Begun=.true.
                     Ntimeforvolumic=Ntime !!!-nint(0.4999999+sgg%OBSERVATION(ii)%InitialTime/sgg%dt)
                     if (mod(Ntimeforvolumic,output(ii)%Trancos)==0) then
                        Ntimeforvolumic=Ntimeforvolumic/output(ii)%Trancos
                        if ( ((at >= sgg%OBSERVATION(ii)%InitialTime).and.(at <= sgg%OBSERVATION(ii)%FinalTime+sgg%dt/2.0_RKIND))) then
                           do KKK = k1, k2
                           if (mod(KKK,output(ii)%item(i)%Ztrancos)==0) then
                           k1t=int(kkk/output(ii)%item(i)%Ztrancos)
                              KKK_m = KKK - b%Ez%ZI
                              do JJJ = j1, j2
                              if (mod(jjj,output(ii)%item(i)%Ytrancos)==0) then
                              j1t=int(jjj/output(ii)%item(i)%Ytrancos)
                                 JJJ_m = JJJ - b%Ez%YI
                                 do III = i1, i2
                                 if (mod(iii,output(ii)%item(i)%Xtrancos)==0) then
                                 i1t=int(iii/output(ii)%item(i)%Xtrancos)
                                    III_m = III - b%Ez%XI
                                    output( ii)%item( i)%valor3D(Ntimeforvolumic,i1t,j1t,k1t) =  Ez( III_m, JJJ_m, KKK_m)
                                 endif
                                 enddo
                              endif
                              enddo
                           endif
                           enddo
                        endif
                     endif
!por aqui voy con los i1t, j1t, k1t
                   case( iHxC)
                     at=sgg%tiempo(ntime)
                        if (at > sgg%OBSERVATION(ii)%FinalTime+sgg%dt/2.0_RKIND) sgg%OBSERVATION(ii)%Done=.true.
                        if (at >= sgg%OBSERVATION(ii)%InitialTime)     sgg%OBSERVATION(ii)%Begun=.true.
                     Ntimeforvolumic=Ntime !!!-nint(0.4999999+sgg%OBSERVATION(ii)%InitialTime/sgg%dt)
                     if (mod(Ntimeforvolumic,output(ii)%Trancos)==0) then
                        Ntimeforvolumic=Ntimeforvolumic/output(ii)%Trancos
                        if ( ((at >= sgg%OBSERVATION(ii)%InitialTime).and.(at <= sgg%OBSERVATION(ii)%FinalTime+sgg%dt/2.0_RKIND))) then
                           do KKK = k1, k2
                           if (mod(KKK,output(ii)%item(i)%Ztrancos)==0) then
                           k1t=int(kkk/output(ii)%item(i)%Ztrancos)
                              KKK_m = KKK - b%Hx%ZI
                              do JJJ = j1, j2
                              if (mod(jjj,output(ii)%item(i)%Ytrancos)==0) then
                              j1t=int(jjj/output(ii)%item(i)%Ytrancos)
                                 JJJ_m = JJJ - b%Hx%YI
                                 do III = i1, i2
                                 if (mod(iii,output(ii)%item(i)%Xtrancos)==0) then
                                 i1t=int(iii/output(ii)%item(i)%Xtrancos)
                                    III_m = III - b%Hx%XI
                                    output( ii)%item( i)%valor3D(Ntimeforvolumic,i1t,j1t,k1t) =  Hx( III_m, JJJ_m, KKK_m)
                                 endif
                                 enddo
                              endif
                              enddo
                           endif
                           enddo
                        endif
                     endif
                   case( iHyC)
                     AT=sgg%tiempo(ntime)
                        if (at > sgg%OBSERVATION(ii)%FinalTime+sgg%dt/2.0_RKIND) sgg%OBSERVATION(ii)%Done=.true.
                        if (at >= sgg%OBSERVATION(ii)%InitialTime)     sgg%OBSERVATION(ii)%Begun=.true.
                     Ntimeforvolumic=Ntime !!!-nint(0.4999999+sgg%OBSERVATION(ii)%InitialTime/sgg%dt)
                     if (mod(Ntimeforvolumic,output(ii)%Trancos)==0) then
                        Ntimeforvolumic=Ntimeforvolumic/output(ii)%Trancos
                        if ( ((at >= sgg%OBSERVATION(ii)%InitialTime).and.(at <= sgg%OBSERVATION(ii)%FinalTime+sgg%dt/2.0_RKIND))) then
                           do KKK = k1, k2
                           if (mod(KKK,output(ii)%item(i)%Ztrancos)==0) then
                           k1t=int(kkk/output(ii)%item(i)%Ztrancos)
                              KKK_m = KKK - b%Hy%ZI
                              do JJJ = j1, j2
                              if (mod(jjj,output(ii)%item(i)%Ytrancos)==0) then
                              j1t=int(jjj/output(ii)%item(i)%Ytrancos)
                                 JJJ_m = JJJ - b%Hy%YI
                                 do III = i1, i2
                                 if (mod(iii,output(ii)%item(i)%Xtrancos)==0) then
                                 i1t=int(iii/output(ii)%item(i)%Xtrancos)
                                    III_m = III - b%Hy%XI
                                    output( ii)%item( i)%valor3D(Ntimeforvolumic,i1t,j1t,k1t) =  Hy( III_m, JJJ_m, KKK_m)
                                 endif
                                 enddo
                              endif
                              enddo
                           endif
                           enddo
                        endif
                     endif
                   case( iHzC)
                     at=sgg%tiempo(ntime)
                        if (at > sgg%OBSERVATION(ii)%FinalTime+sgg%dt/2.0_RKIND) sgg%OBSERVATION(ii)%Done=.true.
                        if (at >= sgg%OBSERVATION(ii)%InitialTime)     sgg%OBSERVATION(ii)%Begun=.true.
                     Ntimeforvolumic=Ntime !!!-nint(0.4999999+sgg%OBSERVATION(ii)%InitialTime/sgg%dt)
                     if (mod(Ntimeforvolumic,output(ii)%Trancos)==0) then
                        Ntimeforvolumic=Ntimeforvolumic/output(ii)%Trancos
                        if ( ((at >= sgg%OBSERVATION(ii)%InitialTime).and.(at <= sgg%OBSERVATION(ii)%FinalTime+sgg%dt/2.0_RKIND))) then
                           do KKK = k1, k2
                           if (mod(KKK,output(ii)%item(i)%Ztrancos)==0) then
                           k1t=int(kkk/output(ii)%item(i)%Ztrancos)
                              KKK_m = KKK - b%Hz%ZI
                              do JJJ = j1, j2
                              if (mod(jjj,output(ii)%item(i)%Ytrancos)==0) then
                              j1t=int(jjj/output(ii)%item(i)%Ytrancos)
                                 JJJ_m = JJJ - b%Hz%YI
                                 do III = i1, i2
                                 if (mod(iii,output(ii)%item(i)%Xtrancos)==0) then
                                 i1t=int(iii/output(ii)%item(i)%Xtrancos)
                                    III_m = III - b%Hz%XI
                                    output( ii)%item( i)%valor3D(Ntimeforvolumic,i1t,j1t,k1t) =  Hz( III_m, JJJ_m, KKK_m)
                                 endif
                                 enddo
                              endif
                              enddo
                           endif
                           enddo
                        endif
                     endif
                     !
                   case( iMEC)
                     at=sgg%tiempo(ntime)
                     if (at > sgg%OBSERVATION(ii)%FinalTime+sgg%dt/2.0_RKIND) sgg%OBSERVATION(ii)%Done=.true.
                     if (at >= sgg%OBSERVATION(ii)%InitialTime)     sgg%OBSERVATION(ii)%Begun=.true.
                     Ntimeforvolumic=Ntime !!!-nint(0.4999999+sgg%OBSERVATION(ii)%InitialTime/sgg%dt)
                     if (mod(Ntimeforvolumic,output(ii)%Trancos)==0) then
                        Ntimeforvolumic=Ntimeforvolumic/output(ii)%Trancos
                        if ( ((at >= sgg%OBSERVATION(ii)%InitialTime).and.(at <= sgg%OBSERVATION(ii)%FinalTime+sgg%dt/2.0_RKIND))) then
                           do KKK = k1, k2
                           if (mod(KKK,output(ii)%item(i)%Ztrancos)==0) then
                           k1t=int(kkk/output(ii)%item(i)%Ztrancos)
                              KKK_m = KKK - b%Ez%ZI
                              do JJJ = j1, j2
                              if (mod(jjj,output(ii)%item(i)%Ytrancos)==0) then
                              j1t=int(jjj/output(ii)%item(i)%Ytrancos)
                                 JJJ_m = JJJ - b%Ey%YI
                                 do III = i1, i2
                                 if (mod(iii,output(ii)%item(i)%Xtrancos)==0) then
                                 i1t=int(iii/output(ii)%item(i)%Xtrancos)
                                    III_m = III - b%Ex%XI
                                    output( ii)%item( i)%valor3D(Ntimeforvolumic,i1t,j1t,k1t) =  &
                                    &    sqrt(Ex( III_m, JJJ_m, KKK_m)**2.0_RKIND + Ey( III_m, JJJ_m, KKK_m)**2.0_RKIND + &
                                    Ez( III_m, JJJ_m, KKK_m)**2.0_RKIND )
                                 endif
                                 enddo
                              endif
                              enddo
                           endif
                           enddo
                        endif
                     endif
                   case( iMHC)
                     at=sgg%tiempo(ntime)
                     if (at > sgg%OBSERVATION(ii)%FinalTime+sgg%dt/2.0_RKIND) sgg%OBSERVATION(ii)%Done=.true.
                     if (at >= sgg%OBSERVATION(ii)%InitialTime)     sgg%OBSERVATION(ii)%Begun=.true.
                     Ntimeforvolumic=Ntime !!!-nint(0.4999999+sgg%OBSERVATION(ii)%InitialTime/sgg%dt)
                     if (mod(Ntimeforvolumic,output(ii)%Trancos)==0) then
                        Ntimeforvolumic=Ntimeforvolumic/output(ii)%Trancos
                        if ( ((at >= sgg%OBSERVATION(ii)%InitialTime).and.(at <= sgg%OBSERVATION(ii)%FinalTime+sgg%dt/2.0_RKIND))) then
                           do KKK = k1, k2
                           if (mod(KKK,output(ii)%item(i)%Ztrancos)==0) then
                           k1t=int(kkk/output(ii)%item(i)%Ztrancos)
                              KKK_m = KKK - b%Hz%ZI
                              do JJJ = j1, j2
                              if (mod(jjj,output(ii)%item(i)%Ytrancos)==0) then
                              j1t=int(jjj/output(ii)%item(i)%Ytrancos)
                                 JJJ_m = JJJ - b%Hy%YI
                                 do III = i1, i2
                                 if (mod(iii,output(ii)%item(i)%Xtrancos)==0) then
                                 i1t=int(iii/output(ii)%item(i)%Xtrancos)
                                    III_m = III - b%Hx%XI
                                    output( ii)%item( i)%valor3D(Ntimeforvolumic,i1t,j1t,k1t) =  &
                                    &    sqrt(Hx( III_m, JJJ_m, KKK_m)**2.0_RKIND + Hy( III_m, JJJ_m, KKK_m)**2.0_RKIND + &
                                    Hz( III_m, JJJ_m, KKK_m)**2.0_RKIND )
                                 endif
                                 enddo
                              endif
                              enddo
                           endif
                           enddo
                        endif
                     endif
                     !sondas corriente 20/2/14
                   case( iCur,iCurX,iCurY,iCurZ,mapvtk)
                     at=sgg%tiempo(ntime)
                     if (at > sgg%OBSERVATION(ii)%FinalTime+sgg%dt/2.0_RKIND) sgg%OBSERVATION(ii)%Done=.true.
                     if (at >= sgg%OBSERVATION(ii)%InitialTime)     sgg%OBSERVATION(ii)%Begun=.true.
                     Ntimeforvolumic=Ntime !!!-nint(0.4999999+sgg%OBSERVATION(ii)%InitialTime/sgg%dt)
                     if (mod(Ntimeforvolumic,output(ii)%Trancos)==0) then
                        Ntimeforvolumic=Ntimeforvolumic/output(ii)%Trancos
                        if ( ((at >= sgg%OBSERVATION(ii)%InitialTime).and.(at <= sgg%OBSERVATION(ii)%FinalTime+sgg%dt/2.0_RKIND))) then
                           conta=0 !en el mismo orden que se allocatearon
                           do KKK = k1, k2
                              do JJJ = j1, j2
                                 do III = i1, i2
                                    !saca  current a lo largo del edge con las sondas icur
                                    if (field/=mapvtk) then
                                       if ((sgg%med(sggMiEx(III - b%Ex%XI, JJJ - b%Ex%YI, KKK - b%Ex%ZI))%Is%ThinWire).and.(iii <= SINPML_fullsize(iEx)%XE).and.(jjj <= SINPML_fullsize(iEx)%YE).and.(kkk <= SINPML_fullsize(iEx)%ZE)) then
                                          conta=conta+1
                                          Jx= dyh(JJJ - b%Hy%YI) * (- Hy( III - b%Hy%XI, JJJ - b%Hy%YI, KKK - b%Ex%ZI)+Hy( III - b%Hy%XI  , JJJ - b%Hy%YI  , KKK - b%Hy%ZI-1))  + &  
                                              dzh(KKK - b%Hz%ZI) * (  Hz( III - b%Hz%XI, JJJ - b%Hz%YI, KKK - b%Hz%ZI)-Hz( III - b%Hz%XI  , JJJ - b%Hz%YI-1, KKK - b%Hz%ZI  ))  
                                           output( ii)%item( i)%Serialized%valor_x(Ntimeforvolumic,conta) = Jx
                                           output( ii)%item( i)%Serialized%valor_y(Ntimeforvolumic,conta) = 0.0_RKIND
                                           output( ii)%item( i)%Serialized%valor_z(Ntimeforvolumic,conta) = 0.0_RKIND
                                       endif                                                                
                                       if ((sgg%med(sggMiEy(III - b%Ey%XI, JJJ - b%Ey%YI, KKK - b%Ey%ZI))%Is%ThinWire).and.(iii <= SINPML_fullsize(iEy)%XE).and.(jjj <= SINPML_fullsize(iEy)%YE).and.(kkk <= SINPML_fullsize(iEy)%ZE)) then
                                          conta=conta+1
                                          Jy=dxh(III - b%Hx%XI) * (  Hx( III - b%Hx%XI, JJJ - b%Hx%YI, KKK - b%Hx%ZI)-Hx( III - b%Hx%XI  , JJJ - b%Hx%YI  , KKK - b%Hx%ZI-1)) + &
                                             dzh(KKK - b%Hz%ZI) * ( -Hz( III - b%Hz%XI, JJJ - b%Hz%YI, KKK - b%Hz%ZI)+Hz( III - b%Hz%XI-1, JJJ - b%Hz%YI  , KKK - b%Hz%ZI  ))    
                                           output( ii)%item( i)%Serialized%valor_x(Ntimeforvolumic,conta) = 0.0_RKIND  
                                           output( ii)%item( i)%Serialized%valor_y(Ntimeforvolumic,conta) = Jy 
                                           output( ii)%item( i)%Serialized%valor_z(Ntimeforvolumic,conta) = 0.0_RKIND
                                       endif
                                       if ((sgg%med(sggMiEz(III - b%Ez%XI, JJJ - b%Ez%YI, KKK - b%Ez%ZI))%Is%ThinWire).and.(iii <= SINPML_fullsize(iEz)%XE).and.(jjj <= SINPML_fullsize(iEz)%YE).and.(kkk <= SINPML_fullsize(iEz)%ZE)) then
                                          conta=conta+1
                                          Jz=dyh(JJJ - b%Hy%YI) * (  Hy( III - b%Hy%XI, JJJ - b%Hy%YI, KKK - b%Hy%ZI) - Hy( III - b%Hy%XI-1, JJJ - b%Hy%YI  , KKK - b%Hy%ZI  ))  + &
                                             dxh(III - b%Hx%XI) * ( -Hx( III - b%Hx%XI, JJJ - b%Hx%YI, KKK - b%Hx%ZI) + Hx( III - b%Hx%XI  , JJJ - b%Hx%YI-1, KKK - b%Hx%ZI  ))     
                                          output( ii)%item( i)%Serialized%valor_x(Ntimeforvolumic,conta) = 0.0_RKIND 
                                          output( ii)%item( i)%Serialized%valor_y(Ntimeforvolumic,conta) = 0.0_RKIND
                                          output( ii)%item( i)%Serialized%valor_z(Ntimeforvolumic,conta) = Jz
                                       endif
                                       if (((sgg%med(sggMiEx(III - b%Ex%XI, JJJ - b%Ex%YI, KKK - b%Ex%ZI))%Is%Line.AND.sgg%med(sggMiEx(III - b%Ex%XI, JJJ - b%Ex%YI, KKK - b%Ex%ZI))%Is%PEC)).and.&
                                           (iii <= SINPML_fullsize(iEx)%XE).and.(jjj <= SINPML_fullsize(iEx)%YE).and.(kkk <= SINPML_fullsize(iEx)%ZE)) then
                                          if ((.not.sgg%med(sggMiHy(III - b%Hy%XI  , JJJ - b%Hy%YI  , KKK - b%Hy%ZI  ))%Is%PEC).and. &
                                          (.not.sgg%med(sggMiHy(III - b%Hy%XI  , JJJ - b%Hy%YI  , KKK - b%Hy%ZI-1))%Is%PEC).and. &
                                          (.not.sgg%med(sggMiHz(III - b%Hz%XI  , JJJ - b%Hz%YI  , KKK - b%Hz%ZI  ))%Is%PEC).and. &
                                          (.not.sgg%med(sggMiHz(III - b%Hz%XI  , JJJ - b%Hz%YI-1, KKK - b%Hz%ZI  ))%Is%PEC)) then
                                             conta=conta+1
                                           Jx= dyh(JJJ - b%Hy%YI) * (- Hy( III - b%Hy%XI, JJJ - b%Hy%YI, KKK - b%Ex%ZI)+Hy( III - b%Hy%XI  , JJJ - b%Hy%YI  , KKK - b%Hy%ZI-1)) + &
                                               dzh(KKK - b%Hz%ZI) * (  Hz( III - b%Hz%XI, JJJ - b%Hz%YI, KKK - b%Hz%ZI)-Hz( III - b%Hz%XI  , JJJ - b%Hz%YI-1, KKK - b%Hz%ZI  ))  
                                           output( ii)%item( i)%Serialized%valor_x(Ntimeforvolumic,conta) = Jx
                                           output( ii)%item( i)%Serialized%valor_y(Ntimeforvolumic,conta) = 0.0_RKIND
                                           output( ii)%item( i)%Serialized%valor_z(Ntimeforvolumic,conta) = 0.0_RKIND
                                          endif
                                       endif
                                       if (((sgg%med(sggMiEy(III - b%Ey%XI, JJJ - b%Ey%YI, KKK - b%Ey%ZI))%Is%Line.AND.sgg%med(sggMiEy(III - b%Ey%XI, JJJ - b%Ey%YI, KKK - b%Ey%ZI))%Is%PEC)).and. &
                                           (iii <= SINPML_fullsize(iEy)%XE).and.(jjj <= SINPML_fullsize(iEy)%YE).and.(kkk <= SINPML_fullsize(iEy)%ZE)) then
                                          if ((.not.sgg%med(sggMiHz(III - b%Hz%XI  , JJJ - b%Hz%YI  , KKK - b%Hz%ZI  ))%Is%PEC).and. &
                                          (.not.sgg%med(sggMiHz(III - b%Hz%XI-1, JJJ - b%Hz%YI  , KKK - b%Hz%ZI  ))%Is%PEC).and. &
                                          (.not.sgg%med(sggMiHx(III - b%Hx%XI  , JJJ - b%Hx%YI  , KKK - b%Hx%ZI  ))%Is%PEC).and. &
                                          (.not.sgg%med(sggMiHx(III - b%Hx%XI  , JJJ - b%Hx%YI  , KKK - b%Hx%ZI-1))%Is%PEC)) then
                                             conta=conta+1
                                             Jy=dxh(III - b%Hx%XI) * (  Hx( III - b%Hx%XI, JJJ - b%Hx%YI, KKK - b%Hx%ZI)-Hx( III - b%Hx%XI  , JJJ - b%Hx%YI  , KKK - b%Hx%ZI-1))  + &
                                                dzh(KKK - b%Hz%ZI) * ( -Hz( III - b%Hz%XI, JJJ - b%Hz%YI, KKK - b%Hz%ZI)+Hz( III - b%Hz%XI-1, JJJ - b%Hz%YI  , KKK - b%Hz%ZI  ))   
                                           output( ii)%item( i)%Serialized%valor_x(Ntimeforvolumic,conta) = 0.0_RKIND  
                                           output( ii)%item( i)%Serialized%valor_y(Ntimeforvolumic,conta) = Jy 
                                           output( ii)%item( i)%Serialized%valor_z(Ntimeforvolumic,conta) = 0.0_RKIND
                                          endif
                                       endif
                                       if (((sgg%med(sggMiEz(III - b%Ez%XI, JJJ - b%Ez%YI, KKK - b%Ez%ZI))%Is%Line.AND.sgg%med(sggMiEz(III - b%Ez%XI, JJJ - b%Ez%YI, KKK - b%Ez%ZI))%Is%PEC)).and. &
                                           (iii <= SINPML_fullsize(iEz)%XE).and.(jjj <= SINPML_fullsize(iEz)%YE).and.(kkk <= SINPML_fullsize(iEz)%ZE)) then
                                          if ((.not.sgg%med(sggMiHx(III - b%Hx%XI  , JJJ - b%Hx%YI  , KKK - b%Hx%ZI  ))%Is%PEC).and. &
                                          (.not.sgg%med(sggMiHx(III - b%Hx%XI  , JJJ - b%Hx%YI-1, KKK - b%Hx%ZI  ))%Is%PEC).and. &
                                          (.not.sgg%med(sggMiHy(III - b%Hy%XI  , JJJ - b%Hy%YI  , KKK - b%Hy%ZI  ))%Is%PEC).and. &
                                          (.not.sgg%med(sggMiHy(III - b%Hy%XI-1, JJJ - b%Hy%YI  , KKK - b%Hy%ZI  ))%Is%PEC)) then
                                             conta=conta+1
                                             Jz=dyh(JJJ - b%Hy%YI) * (  Hy( III - b%Hy%XI, JJJ - b%Hy%YI, KKK - b%Hy%ZI) - Hy( III - b%Hy%XI-1, JJJ - b%Hy%YI  , KKK - b%Hy%ZI  ))  + &
                                                dxh(III - b%Hx%XI) * ( -Hx( III - b%Hx%XI, JJJ - b%Hx%YI, KKK - b%Hx%ZI) + Hx( III - b%Hx%XI  , JJJ - b%Hx%YI-1, KKK - b%Hx%ZI  ))    
                                             output( ii)%item( i)%Serialized%valor(Ntimeforvolumic,conta)  = Jz
                                             output( ii)%item( i)%Serialized%valor_x(Ntimeforvolumic,conta) = 0.0_RKIND
                                             output( ii)%item( i)%Serialized%valor_y(Ntimeforvolumic,conta) = 0.0_RKIND
                                             output( ii)%item( i)%Serialized%valor_z(Ntimeforvolumic,conta) = Jz
                                          endif
                                       endif
                                    else !si es mapvtk
                                       imed =sggMiEx(III -b%Ex%XI, JJJ- b%Ex%YI  , KKK- b%Ex%ZI  )
                                       imed1=sggMiHy(III -b%Hy%XI, JJJ- b%Hy%YI  , KKK- b%Hy%ZI  )
                                       imed2=sggMiHy(III -b%Hy%XI, JJJ- b%Hy%YI  , KKK- b%Hy%ZI-1)
                                       imed3=sggMiHz(III -b%Hz%XI, JJJ- b%Hz%YI  , KKK- b%Hz%ZI  )
                                       imed4=sggMiHz(III -b%Hz%XI, JJJ- b%Hz%YI-1, KKK- b%Hz%ZI  )
                                       call contabordes(sgg,imed,imed1,imed2,imed3,imed4,EsBorde,SINPML_fullsize,iEx,iii,jjj,kkk)
                                       if (EsBorde) then
                                          conta=conta+1
                                          jJx=sggMiEx(III -b%Ex%XI, JJJ- b%Ex%YI, KKK- b%Ex%ZI)
                                          !!!discretizo los colores para saber mejor que son (27/06/15)
                                          if ((jJx==0).or.(sgg%Med(jJx)%is%Pec)) then
                                             jx=0.5_RKIND
                                          elseif (sgg%Med(jJx)%is%thinwire) then
                                             if (((sggMiEy(III -b%Ey%XI  , JJJ- b%Ey%YI  , KKK- b%Ey%ZI  )/=1).AND.(.not.sgg%med(sggMiEy(III -b%Ey%XI  , JJJ- b%Ey%YI  , KKK- b%Ey%ZI  ))%is%thinwire)).or. &
                                             ((sggMiEy(III -b%Ey%XI  , JJJ- b%Ey%YI-1, KKK- b%Ey%ZI  )/=1).AND.(.not.sgg%med(sggMiEy(III -b%Ey%XI  , JJJ- b%Ey%YI-1, KKK- b%Ey%ZI  ))%is%thinwire)).or. &
                                             ((sggMiEz(III -b%Ez%XI  , JJJ- b%Ez%YI  , KKK- b%Ez%ZI  )/=1).AND.(.not.sgg%med(sggMiEz(III -b%Ez%XI  , JJJ- b%Ez%YI  , KKK- b%Ez%ZI  ))%is%thinwire)).or. &
                                             ((sggMiEz(III -b%Ez%XI  , JJJ- b%Ez%YI  , KKK- b%Ez%ZI-1)/=1).AND.(.not.sgg%med(sggMiEz(III -b%Ez%XI  , JJJ- b%Ez%YI  , KKK- b%Ez%ZI-1))%is%thinwire)).or. &
                                             ((sggMiEy(III -b%Ey%XI+1, JJJ- b%Ey%YI  , KKK- b%Ey%ZI  )/=1).AND.(.not.sgg%med(sggMiEy(III -b%Ey%XI+1, JJJ- b%Ey%YI  , KKK- b%Ey%ZI  ))%is%thinwire)).or. &
                                             ((sggMiEy(III -b%Ey%XI+1, JJJ- b%Ey%YI-1, KKK- b%Ey%ZI  )/=1).AND.(.not.sgg%med(sggMiEy(III -b%Ey%XI+1, JJJ- b%Ey%YI-1, KKK- b%Ey%ZI  ))%is%thinwire)).or. &
                                             ((sggMiEz(III -b%Ez%XI+1, JJJ- b%Ez%YI  , KKK- b%Ez%ZI  )/=1).AND.(.not.sgg%med(sggMiEz(III -b%Ez%XI+1, JJJ- b%Ez%YI  , KKK- b%Ez%ZI  ))%is%thinwire)).or. &
                                             ((sggMiEz(III -b%Ez%XI+1, JJJ- b%Ez%YI  , KKK- b%Ez%ZI-1)/=1).AND.(.not.sgg%med(sggMiEz(III -b%Ez%XI+1, JJJ- b%Ez%YI  , KKK- b%Ez%ZI-1))%is%thinwire))) then
                                                jx=8
                                             else  !no hay una colision
                                                jx=7
                                             endif
                                          elseif ((sgg%Med(jJx)%is%SGBC).or.(sgg%Med(jJx)%is%multiport).or.(sgg%Med(jJx)%is%anismultiport)) then
                                             jx=3.5
                                          elseif ((sgg%Med(jJx)%is%edispersive).or.(sgg%Med(jJx)%is%EDispersiveANIS).or.(sgg%Med(jJx)%is%mDispersive).or.(sgg%Med(jJx)%is%mDispersiveANIS)) then
                                             jx=1.5
                                          elseif ((sgg%Med(jJx)%is%Dielectric).or.(sgg%Med(jJx)%is%Anisotropic)) then
                                             jx=2.5
                                          elseif (sgg%Med(jJx)%is%thinslot) then
                                             jx=4.5
                                          elseif ((sgg%Med(jJx)%is%already_YEEadvanced_byconformal).and.(.not.noconformalmapvtk)) then
                                             jx=5.5
                                          elseif ((sgg%Med(jJx)%is%split_and_useless).and.(.not.noconformalmapvtk)) then 
                                             jx=6.5
                                          else
                                             jx=-0.5_RKIND
                                          endif
                                          output( ii)%item( i)%Serialized%valor(Ntimeforvolumic,conta) = jx    
                                          !!!fin discretizo los colores para saber mejor que son
                                       endif
                                       imed =sggMiEy(III- b%Ey%XI  , JJJ- b%Ey%YI  , KKK- b%Ey%ZI  )
                                       imed1=sggMiHz(III -b%Hz%XI  , JJJ- b%Hz%YI  , KKK- b%Hz%ZI  )
                                       imed2=sggMiHz(III -b%Hz%XI-1, JJJ- b%Hz%YI  , KKK- b%Hz%ZI  )
                                       imed3=sggMiHx(III -b%Hx%XI  , JJJ- b%Hx%YI  , KKK- b%Hx%ZI  )
                                       imed4=sggMiHx(III -b%Hx%XI  , JJJ- b%Hx%YI  , KKK- b%Hx%ZI-1)
                                       call contabordes(sgg,imed,imed1,imed2,imed3,imed4,EsBorde,SINPML_fullsize,iEy,iii,jjj,kkk)
                                       if (EsBorde) then
                                          conta=conta+1
                                          jJy=sggMiEy(III -b%Ey%XI, JJJ- b%Ey%YI, KKK- b%Ey%ZI)
                                          !!!discretizo los colores para saber mejor que son (27/06/15)
                                          if ((jJy==0).or.(sgg%Med(jJy)%is%Pec)) then
                                             jy=0.5_RKIND
                                          elseif (sgg%Med(jJy)%is%thinwire) then
                                             if (((sggMiEz(III -b%Ez%XI  , JJJ- b%Ez%YI    , KKK- b%Ez%ZI  )/=1).AND.(.not.sgg%med(sggMiEz(III -b%Ez%XI  , JJJ- b%Ez%YI    , KKK- b%Ez%ZI  ))%is%thinwire)).or. &
                                             ((sggMiEz(III -b%Ez%XI  , JJJ- b%Ez%YI    , KKK- b%Ez%ZI-1)/=1).AND.(.not.sgg%med(sggMiEz(III -b%Ez%XI  , JJJ- b%Ez%YI    , KKK- b%Ez%ZI-1))%is%thinwire)).or. &
                                             ((sggMiEx(III -b%Ex%XI  , JJJ- b%Ex%YI    , KKK- b%Ex%ZI  )/=1).AND.(.not.sgg%med(sggMiEx(III -b%Ex%XI  , JJJ- b%Ex%YI    , KKK- b%Ex%ZI  ))%is%thinwire)).or. &
                                             ((sggMiEx(III -b%Ex%XI-1, JJJ- b%Ex%YI    , KKK- b%Ex%ZI  )/=1).AND.(.not.sgg%med(sggMiEx(III -b%Ex%XI-1, JJJ- b%Ex%YI    , KKK- b%Ex%ZI  ))%is%thinwire)).or. &
                                             ((sggMiEz(III -b%Ez%XI  , JJJ- b%Ez%YI+1  , KKK- b%Ez%ZI  )/=1).AND.(.not.sgg%med(sggMiEz(III -b%Ez%XI  , JJJ- b%Ez%YI+1  , KKK- b%Ez%ZI  ))%is%thinwire)).or. &
                                             ((sggMiEz(III -b%Ez%XI  , JJJ- b%Ez%YI+1  , KKK- b%Ez%ZI-1)/=1).AND.(.not.sgg%med(sggMiEz(III -b%Ez%XI  , JJJ- b%Ez%YI+1  , KKK- b%Ez%ZI-1))%is%thinwire)).or. &
                                             ((sggMiEx(III -b%Ex%XI  , JJJ- b%Ex%YI+1  , KKK- b%Ex%ZI  )/=1).AND.(.not.sgg%med(sggMiEx(III -b%Ex%XI  , JJJ- b%Ex%YI+1  , KKK- b%Ex%ZI  ))%is%thinwire)).or. &
                                             ((sggMiEx(III -b%Ex%XI-1, JJJ- b%Ex%YI+1  , KKK- b%Ex%ZI  )/=1).AND.(.not.sgg%med(sggMiEx(III -b%Ex%XI-1, JJJ- b%Ex%YI+1  , KKK- b%Ex%ZI  ))%is%thinwire))) then
                                                jy=8
                                             else  !no hay una colision
                                                jy=7
                                             endif
                                          elseif ((sgg%Med(jJy)%is%SGBC).or.(sgg%Med(jJy)%is%multiport).or.(sgg%Med(jJy)%is%anismultiport)) then
                                             jy=3.5
                                          elseif ((sgg%Med(jJy)%is%edispersive).or.(sgg%Med(jJy)%is%EDispersiveANIS).or.(sgg%Med(jJy)%is%mDispersive).or.(sgg%Med(jJy)%is%mDispersiveANIS)) then
                                             jy=1.5
                                          elseif ((sgg%Med(jJy)%is%Dielectric).or.(sgg%Med(jJy)%is%Anisotropic)) then
                                             jy=2.5
                                          elseif (sgg%Med(jJy)%is%thinslot) then
                                             jy=4.5
                                          elseif ((sgg%Med(jJy)%is%already_YEEadvanced_byconformal).and.(.not.noconformalmapvtk)) then 
                                             jy=5.5
                                          elseif ((sgg%Med(jJy)%is%split_and_useless).and.(.not.noconformalmapvtk)) then
                                             jy=6.5
                                          else
                                             jy=-0.5_RKIND
                                          endif
                                          output( ii)%item( i)%Serialized%valor(Ntimeforvolumic,conta) = jy    
                                          !!!fin discretizo los colores para saber mejor que son
                                       endif
                                       imed =sggMiEz(III- b%Ez%XI  , JJJ- b%Ez%YI  , KKK- b%Ez%ZI  )
                                       imed1=sggMiHx(III -b%Hx%XI  , JJJ- b%Hx%YI  , KKK- b%Hx%ZI  )
                                       imed2=sggMiHx(III -b%Hx%XI  , JJJ- b%Hx%YI-1, KKK- b%Hx%ZI  )
                                       imed3=sggMiHy(III -b%Hy%XI  , JJJ- b%Hy%YI  , KKK- b%Hy%ZI  )
                                       imed4=sggMiHy(III -b%Hy%XI-1, JJJ- b%Hy%YI  , KKK- b%Hy%ZI  )
                                       call contabordes(sgg,imed,imed1,imed2,imed3,imed4,EsBorde,SINPML_fullsize,iEz,iii,jjj,kkk)
                                       if (EsBorde) then
                                          conta=conta+1
                                          jJz=sggMiEz(III -b%Ez%XI, JJJ- b%Ez%YI, KKK- b%Ez%ZI)
                                          !!!discretizo los colores para saber mejor que son (27/06/15)
                                          if ((jJz==0).or.(sgg%Med(jJz)%is%Pec)) then
                                             jz=0.5_RKIND
                                          elseif (sgg%Med(jJz)%is%thinwire) then
                                             if (((sggMiEy(III -b%Ey%XI  , JJJ- b%Ey%YI  , KKK- b%Ey%ZI    )/=1).AND.(.not.sgg%med(sggMiEy(III -b%Ey%XI  , JJJ- b%Ey%YI  , KKK- b%Ey%ZI    ))%is%thinwire)).or. &
                                             ((sggMiEy(III -b%Ey%XI  , JJJ- b%Ey%YI-1, KKK- b%Ey%ZI    )/=1).AND.(.not.sgg%med(sggMiEy(III -b%Ey%XI  , JJJ- b%Ey%YI-1, KKK- b%Ey%ZI    ))%is%thinwire)).or. &
                                             ((sggMiEx(III -b%Ex%XI  , JJJ- b%Ex%YI  , KKK- b%Ex%ZI    )/=1).AND.(.not.sgg%med(sggMiEx(III -b%Ex%XI  , JJJ- b%Ex%YI  , KKK- b%Ex%ZI    ))%is%thinwire)).or. &
                                             ((sggMiEx(III -b%Ex%XI-1, JJJ- b%Ex%YI  , KKK- b%Ex%ZI    )/=1).AND.(.not.sgg%med(sggMiEx(III -b%Ex%XI-1, JJJ- b%Ex%YI  , KKK- b%Ex%ZI    ))%is%thinwire)).or. &
                                             ((sggMiEy(III -b%Ey%XI  , JJJ- b%Ey%YI  , KKK- b%Ey%ZI+1  )/=1).AND.(.not.sgg%med(sggMiEy(III -b%Ey%XI  , JJJ- b%Ey%YI  , KKK- b%Ey%ZI+1  ))%is%thinwire)).or. &
                                             ((sggMiEy(III -b%Ey%XI  , JJJ- b%Ey%YI-1, KKK- b%Ey%ZI+1  )/=1).AND.(.not.sgg%med(sggMiEy(III -b%Ey%XI  , JJJ- b%Ey%YI-1, KKK- b%Ey%ZI+1  ))%is%thinwire)).or. &
                                             ((sggMiEx(III -b%Ex%XI  , JJJ- b%Ex%YI  , KKK- b%Ex%ZI+1  )/=1).AND.(.not.sgg%med(sggMiEx(III -b%Ex%XI  , JJJ- b%Ex%YI  , KKK- b%Ex%ZI+1  ))%is%thinwire)).or. &
                                             ((sggMiEx(III -b%Ex%XI-1, JJJ- b%Ex%YI  , KKK- b%Ex%ZI+1  )/=1).AND.(.not.sgg%med(sggMiEx(III -b%Ex%XI-1, JJJ- b%Ex%YI  , KKK- b%Ex%ZI+1  ))%is%thinwire))) then
                                                jz=8
                                             else  !no hay una colision
                                                jz=7
                                             endif
                                          elseif ((sgg%Med(jJz)%is%SGBC).or.(sgg%Med(jJz)%is%multiport).or.(sgg%Med(jJz)%is%anismultiport)) then
                                             jz=3.5
                                          elseif ((sgg%Med(jJz)%is%edispersive).or.(sgg%Med(jJz)%is%EDispersiveANIS).or.(sgg%Med(jJz)%is%mDispersive).or.(sgg%Med(jJz)%is%mDispersiveANIS)) then
                                             jz=1.5
                                          elseif ((sgg%Med(jJz)%is%Dielectric).or.(sgg%Med(jJz)%is%Anisotropic)) then
                                             jz=2.5
                                          elseif (sgg%Med(jJz)%is%thinslot) then
                                             jz=4.5
                                          elseif ((sgg%Med(jJz)%is%already_YEEadvanced_byconformal).and.(.not.noconformalmapvtk)) then
                                              if (jjj==10.and.iii==16) then
                                                  continue !solo para debugeo 0624 ojoo deshacer
                                              endif
                                             jz=5.5
                                          elseif ((sgg%Med(jJz)%is%split_and_useless).and.(.not.noconformalmapvtk)) then 
                                             jz=6.5
                                          else
                                             jz=-0.5_RKIND
                                          endif
                                          !!!fin discretizo los colores para saber mejor que son
                                          output( ii)%item( i)%Serialized%valor(Ntimeforvolumic,conta) = jz  
                                       endif

                                    endif
                                 end do
                              end do
                           end do

                           !!!
                           if (field==mapvtk) then
                              INIT=.false.; geom=.false. ; asigna=.true.; magnetic=.false. ; electric=.true.
#ifdef CompileWithNodalSources
                              call nodalvtk(sgg,sggMiEx,sggMiEy,sggMiEz,sggMiHx,sggMiHy,sggMiHz,sggMtag,b,&
                                            init,geom,asigna,electric,magnetic,conta,i,ii,output,Ntimeforvolumic)
#endif                              
#ifdef CompileWithWires
                              call wirebundlesvtk(sgg,init,geom,asigna,conta,i,ii,output,Ntimeforvolumic,wiresflavor,sggMtag)
#endif
                           endif
                           !!!
                           do KKK = k1, k2
                              do JJJ = j1, j2
                                 do III = i1, i2      
                                   !saca current en surfaces 0124
                                    if (field/=mapvtk) then
                                       if (((sgg%med(sggMiHx(III -b%Hx%XI, JJJ- b%Hx%YI, KKK- b%Hx%ZI))%Is%PEC).or. &
                                       (sgg%med(sggMiHx(III -b%Hx%XI, JJJ- b%Hx%YI, KKK- b%Hx%ZI))%Is%Surface).or.  &
                                       (field==icurX)).and.(iii <= SINPML_fullsize(iHx)%XE).and.(jjj <= SINPML_fullsize(iHx)%YE).and.(kkk <= SINPML_fullsize(iHx)%ZE)) then
                                          conta=conta+1
                                          Jy=(dzh(KKK - b%Hz%ZI) * Hz( III - b%Hz%XI-1, JJJ - b%Hz%YI  , KKK - b%Hz%ZI  ) + dzh(KKK - b%Hz%ZI+1) *Hz( III - b%Hz%XI-1, JJJ - b%Hz%YI  , KKK - b%Hz%ZI+1) )/2.0_RKIND -  &
                                             (dzh(KKK - b%Hz%ZI) * Hz( III - b%Hz%XI  , JJJ - b%Hz%YI  , KKK - b%Hz%ZI  ) + dzh(KKK - b%Hz%ZI+1) *Hz( III - b%Hz%XI  , JJJ - b%Hz%YI  , KKK - b%Hz%ZI+1) )/2.0_RKIND +  &
                                              dxh(III - b%Hx%XI)*( Hx( III - b%Hx%XI  , JJJ - b%Hx%YI  , KKK - b%Hx%ZI+1) -                       Hx( III - b%Hx%XI  , JJJ - b%Hx%YI  , KKK - b%Hx%ZI-1) )/2.0_RKIND
                                          !el Hx al promediarlo con el suyo (i,j,k) a ambos lados pierde su componente y solo quedan las adyancentes
                                          Jz=(dyh(JJJ - b%Hy%YI) * Hy( III - b%Hy%XI  , JJJ - b%Hy%YI  , KKK - b%Hy%ZI  ) + dyh(JJJ - b%Hy%YI+1) *Hy( III - b%Hy%XI  , JJJ - b%Hy%YI+1, KKK - b%Hy%ZI  ) )/2.0_RKIND -  &
                                             (dyh(JJJ - b%Hy%YI) * Hy( III - b%Hy%XI-1, JJJ - b%Hy%YI  , KKK - b%Hy%ZI  ) + dyh(JJJ - b%Hy%YI+1) *Hy( III - b%Hy%XI-1, JJJ - b%Hy%YI+1, KKK - b%Hy%ZI  ) )/2.0_RKIND +  &
                                              dxh(III - b%Hx%XI)*( Hx( III - b%Hx%XI  , JJJ - b%Hx%YI-1, KKK - b%Hx%ZI  ) -                       Hx( III - b%Hx%XI  , JJJ - b%Hx%YI+1, KKK - b%Hx%ZI  ) )/2.0_RKIND 
                                          output( ii)%item( i)%Serialized%valor_x(Ntimeforvolumic,conta) = 0.0_RKIND
                                          output( ii)%item( i)%Serialized%valor_y(Ntimeforvolumic,conta) = Jy
                                          output( ii)%item( i)%Serialized%valor_z(Ntimeforvolumic,conta) = Jz
                                       endif
                                       if (((sgg%med(sggMiHy(III- b%Hy%XI, JJJ- b%Hy%YI, KKK- b%Hy%ZI))%Is%PEC).or. &
                                       (sgg%med(sggMiHy(III- b%Hy%XI, JJJ- b%Hy%YI, KKK- b%Hy%ZI))%Is%Surface).or. &
                                       (field==icurY)).and.(iii <= SINPML_fullsize(iHy)%XE).and.(jjj <= SINPML_fullsize(iHy)%YE).and.(kkk <= SINPML_fullsize(iHy)%ZE)) then
                                          conta=conta+1
                                          Jz=(dxh(III - b%Hx%XI) * Hx( III - b%Hx%XI  , JJJ - b%Hx%YI-1, KKK - b%Hx%ZI  ) + dxh(III - b%Hx%XI+1) *Hx( III - b%Hx%XI+1, JJJ - b%Hx%YI-1, KKK - b%Hx%ZI  ) )/2.0_RKIND -  &
                                             (dxh(III - b%Hx%XI) * Hx( III - b%Hx%XI  , JJJ - b%Hx%YI  , KKK - b%Hx%ZI  ) + dxh(III - b%Hx%XI+1) *Hx( III - b%Hx%XI+1, JJJ - b%Hx%YI  , KKK - b%Hx%ZI  ) )/2.0_RKIND +  &
                                              dyh(JJJ - b%Hy%YI)*( Hy( III - b%Hy%XI+1, JJJ - b%Hy%YI  , KKK - b%Hy%ZI  ) -                       Hy( III - b%Hy%XI-1, JJJ - b%Hy%YI  , KKK - b%Hy%ZI  ) )/2.0_RKIND
                                          !
                                          Jx=(dzh(KKK - b%Hz%ZI) * Hz( III - b%Hz%XI  , JJJ - b%Hz%YI  , KKK - b%Hz%ZI  ) + dzh(KKK - b%Hz%ZI+1) *Hz( III - b%Hz%XI  , JJJ - b%Hz%YI  , KKK - b%Hz%ZI+1) )/2.0_RKIND -  &
                                             (dzh(KKK - b%Hz%ZI) * Hz( III - b%Hz%XI  , JJJ - b%Hz%YI-1, KKK - b%Hz%ZI  ) + dzh(KKK - b%Hz%ZI+1) *Hz( III - b%Hz%XI  , JJJ - b%Hz%YI-1, KKK - b%Hz%ZI+1) )/2.0_RKIND +  &
                                              dyh(JJJ - b%Hy%YI)*( Hy( III - b%Hy%XI  , JJJ - b%Hy%YI  , KKK - b%Hy%ZI-1) -                       Hy( III - b%Hy%XI  , JJJ - b%Hy%YI  , KKK - b%Hy%ZI+1) )/2.0_RKIND
                                          !
                                          output( ii)%item( i)%Serialized%valor_x(Ntimeforvolumic,conta) = Jx
                                          output( ii)%item( i)%Serialized%valor_y(Ntimeforvolumic,conta) = 0.0_RKIND
                                          output( ii)%item( i)%Serialized%valor_z(Ntimeforvolumic,conta) = Jz
                                       endif
                                       if (((sgg%med(sggMiHz(III- b%Hz%XI, JJJ- b%Hz%YI, KKK- b%Hz%ZI))%Is%PEC).or. &
                                       (sgg%med(sggMiHz(III- b%Hz%XI, JJJ- b%Hz%YI, KKK- b%Hz%ZI))%Is%Surface).or. &
                                       (field==icurZ)).and.(iii <= SINPML_fullsize(iHz)%XE).and.(jjj <= SINPML_fullsize(iHz)%YE).and.(kkk <= SINPML_fullsize(iHz)%ZE)) then
                                          conta=conta+1
                                          Jx=(dyh(JJJ - b%Hy%YI) * Hy( III - b%Hy%XI  , JJJ - b%Hy%YI  , KKK - b%Hy%ZI-1) + dyh(JJJ - b%Hy%YI+1) *Hy( III - b%Hy%XI  , JJJ - b%Hy%YI+1, KKK - b%Hy%ZI-1) )/2.0_RKIND -  &
                                             (dyh(JJJ - b%Hy%YI) * Hy( III - b%Hy%XI  , JJJ - b%Hy%YI  , KKK - b%Hy%ZI  ) + dyh(JJJ - b%Hy%YI+1) *Hy( III - b%Hy%XI  , JJJ - b%Hy%YI+1, KKK - b%Hy%ZI  ) )/2.0_RKIND +  &
                                              dzh(KKK - b%Hz%ZI)*( Hz( III - b%Hz%XI  , JJJ - b%Hz%YI+1, KKK - b%Hz%ZI  ) -                       Hz( III - b%Hz%XI  , JJJ - b%Hz%YI-1, KKK - b%Hz%ZI  ) )/2.0_RKIND
                                          !
                                          Jy=(dxh(III - b%Hx%XI) * Hx( III - b%Hx%XI  , JJJ - b%Hx%YI  , KKK - b%Hx%ZI  ) + dxh(III - b%Hx%XI+1) *Hx( III - b%Hx%XI+1, JJJ - b%Hx%YI  , KKK - b%Hx%ZI  ) )/2.0_RKIND -  &
                                             (dxh(III - b%Hx%XI) * Hx( III - b%Hx%XI  , JJJ - b%Hx%YI  , KKK - b%Hx%ZI-1) + dxh(III - b%Hx%XI+1) *Hx( III - b%Hx%XI+1, JJJ - b%Hx%YI  , KKK - b%Hx%ZI-1) )/2.0_RKIND +  &
                                              dzh(KKK - b%Hz%ZI)*( Hz( III - b%Hz%XI-1, JJJ - b%Hz%YI  , KKK - b%Hz%ZI  ) -                       Hz( III - b%Hz%XI+1, JJJ - b%Hz%YI  , KKK - b%Hz%ZI  ) )/2.0_RKIND
                                          !
                                           output( ii)%item( i)%Serialized%valor_x(Ntimeforvolumic,conta) = Jx
                                           output( ii)%item( i)%Serialized%valor_y(Ntimeforvolumic,conta) = Jy
                                           output( ii)%item( i)%Serialized%valor_z(Ntimeforvolumic,conta) = 0.0_RKIND
                                       endif
                                    else                                       !si es mapvtk y si no es vacio, asimilo la salida a corrientes iBloqueJ? para que vtk.f90 los escriba en quads
                                       if ((sggMiHx(III -b%Hx%XI, JJJ- b%Hx%YI, KKK- b%Hx%ZI)/=1).and. &
                                       (.not.sgg%med(sggMiHx(III -b%Hx%XI, JJJ- b%Hx%YI, KKK- b%Hx%ZI))%is%PML).and.(iii <= SINPML_fullsize(iHx)%XE).and.(jjj <= SINPML_fullsize(iHx)%YE).and.(kkk <= SINPML_fullsize(iHx)%ZE)) then
                                          conta=conta+1
                                          jJx=sggMiHx(III -b%Hx%XI, JJJ- b%Hx%YI, KKK- b%Hx%ZI)
                                          !!!discretizo los colores para saber mejor que son (27/06/15)
                                          if ((jJx==0).or.(sgg%Med(jJx)%is%Pec)) then
                                             jx=0
                                          elseif (sgg%Med(jJx)%is%thinwire) then
                                             CALL StopOnError (0,1,'ERROR: A magnetic field cannot be a thin-wire')
                                          elseif ((sgg%Med(jJx)%is%SGBC).or.(sgg%Med(jJx)%is%multiport).or.(sgg%Med(jJx)%is%anismultiport)) then
                                             jx=300+jJx
                                          elseif ((sgg%Med(jJx)%is%edispersive).or.(sgg%Med(jJx)%is%EDispersiveANIS).or.(sgg%Med(jJx)%is%mDispersive).or.(sgg%Med(jJx)%is%mDispersiveANIS)) then
                                             jx=100+jJx
                                          elseif ((sgg%Med(jJx)%is%Dielectric).or.(sgg%Med(jJx)%is%Anisotropic)) then
                                             jx=200+jJx
                                          elseif (sgg%Med(jJx)%is%thinslot) then
                                             jx=400+jJx
                                          elseif ((sgg%Med(jJx)%is%already_YEEadvanced_byconformal).and.(.not.noconformalmapvtk)) then
                                             jx=5
                                          elseif ((sgg%Med(jJx)%is%split_and_useless).and.(.not.noconformalmapvtk)) then
                                             jx=6
                                          else
                                             jx=-1
                                          endif
                                          !!!fin discretizo los colores para saber mejor que son
                                          output( ii)%item( i)%Serialized%valor(Ntimeforvolumic,conta) = Jx    
                                       endif
                                       if ((sggMiHy(III- b%Hy%XI, JJJ- b%Hy%YI, KKK- b%Hy%ZI)/=1).and. &
                                       (.not.sgg%med(sggMiHy(III -b%Hy%XI, JJJ- b%Hy%YI, KKK- b%Hy%ZI))%is%PML).and.(iii <= SINPML_fullsize(iHy)%XE).and.(jjj <= SINPML_fullsize(iHy)%YE).and.(kkk <= SINPML_fullsize(iHy)%ZE)) then
                                          conta=conta+1
                                          jJy=sggMiHy(III- b%Hy%XI, JJJ- b%Hy%YI, KKK- b%Hy%ZI)
                                          !!!discretizo los colores para saber mejor que son (27/06/15)
                                          if ((jJy==0).or.(sgg%Med(jJy)%is%Pec)) then
                                             jy=0
                                          elseif (sgg%Med(jJy)%is%thinwire) then
                                             CALL StopOnError (0,1,'ERROR: A magnetic field cannot be a thin-wire')
                                          elseif ((sgg%Med(jJy)%is%SGBC).or.(sgg%Med(jJy)%is%multiport).or.(sgg%Med(jJy)%is%anismultiport)) then
                                             jy=300+jJy
                                          elseif ((sgg%Med(jJy)%is%edispersive).or.(sgg%Med(jJy)%is%EDispersiveANIS).or.(sgg%Med(jJy)%is%mDispersive).or.(sgg%Med(jJy)%is%mDispersiveANIS)) then
                                             jy=100+jJy
                                          elseif ((sgg%Med(jJy)%is%Dielectric).or.(sgg%Med(jJy)%is%Anisotropic)) then
                                             jy=200+jJy
                                          elseif (sgg%Med(jJy)%is%thinslot) then
                                             jy=400+jJy
                                          elseif ((sgg%Med(jJy)%is%already_YEEadvanced_byconformal).and.(.not.noconformalmapvtk)) then
                                             jy=5
                                          elseif ((sgg%Med(jJy)%is%split_and_useless) .and.(.not.noconformalmapvtk)) then
                                             jy=6
                                          else
                                             jy=-1
                                          endif
                                          !!!fin discretizo los colores para saber mejor que son
                                          output( ii)%item( i)%Serialized%valor(Ntimeforvolumic,conta) = Jy     
                                       endif
                                       if ((sggMiHz(III- b%Hz%XI, JJJ- b%Hz%YI, KKK- b%Hz%ZI)/=1).and. &
                                       (.not.sgg%med(sggMiHz(III -b%Hz%XI, JJJ- b%Hz%YI, KKK- b%Hz%ZI))%is%PML).and.(iii <= SINPML_fullsize(iHz)%XE).and.(jjj <= SINPML_fullsize(iHz)%YE).and.(kkk <= SINPML_fullsize(iHz)%ZE)) then
                                          conta=conta+1
                                          jJz=sggMiHz(III- b%Hz%XI, JJJ- b%Hz%YI, KKK- b%Hz%ZI)
                                          !!!discretizo los colores para saber mejor que son (27/06/15)
                                          if ((jJz==0).or.(sgg%Med(jJz)%is%Pec)) then
                                             jz=0
                                          elseif (sgg%Med(jJz)%is%thinwire) then
                                             CALL StopOnError (0,1,'ERROR: A magnetic field cannot be a thin-wire')
                                          elseif ((sgg%Med(jJz)%is%SGBC).or.(sgg%Med(jJz)%is%multiport).or.(sgg%Med(jJz)%is%anismultiport)) then
                                             jz=300+jJz
                                          elseif ((sgg%Med(jJz)%is%edispersive).or.(sgg%Med(jJz)%is%EDispersiveANIS).or.(sgg%Med(jJz)%is%mDispersive).or.(sgg%Med(jJz)%is%mDispersiveANIS)) then
                                             jz=100+jJz
                                          elseif ((sgg%Med(jJz)%is%Dielectric).or.(sgg%Med(jJz)%is%Anisotropic)) then
                                             jz=200+jJz
                                          elseif (sgg%Med(jJz)%is%thinslot) then
                                             jz=400+jJz
                                          elseif ((sgg%Med(jJz)%is%already_YEEadvanced_byconformal).and.(.not.noconformalmapvtk)) then
                                             jz=5
                                          elseif ((sgg%Med(jJz)%is%split_and_useless).and.(.not.noconformalmapvtk)) then
                                             jz=6
                                          else
                                             jz=-1
                                          endif
                                          !!!fin discretizo los colores para saber mejor que son
                                          output( ii)%item( i)%Serialized%valor(Ntimeforvolumic,conta) = Jz      
                                       endif
                                     ! los tags 141020 para mapvtk se quedan con el medio -100: es una forma de voidearlos para visualizacion
                                       if (sggMtag(iii,jjj,kkk)<0) then
                                           if ( (btest(iabs(sggMtag(iii,jjj,kkk)),3)).and. & 
                                           (.not.sgg%med(sggMiHx(III -b%Hx%XI, JJJ- b%Hx%YI, KKK- b%Hx%ZI))%is%PML).and.(iii <= SINPML_fullsize(iHx)%XE).and.(jjj <= SINPML_fullsize(iHx)%YE).and.(kkk <= SINPML_fullsize(iHx)%ZE)) then
                                              conta=conta+1
                                              jx=-100 !ojo pq el vacio que es 1 pasa a ser -100 si es un candidato a slot
                                              output( ii)%item( i)%Serialized%valor(Ntimeforvolumic,conta) = Jx  
                                           endif
                                           if ( (btest(iabs(sggMtag(iii,jjj,kkk)),4)).and. & 
                                           (.not.sgg%med(sggMiHy(III -b%Hy%XI, JJJ- b%Hy%YI, KKK- b%Hy%ZI))%is%PML).and.(iii <= SINPML_fullsize(iHy)%XE).and.(jjj <= SINPML_fullsize(iHy)%YE).and.(kkk <= SINPML_fullsize(iHy)%ZE)) then
                                              conta=conta+1
                                              jy=-100
                                              output( ii)%item( i)%Serialized%valor(Ntimeforvolumic,conta) = Jy
                                           endif
                                           if ( (btest(iabs(sggMtag(iii,jjj,kkk)),5)).and. &
                                           (.not.sgg%med(sggMiHz(III -b%Hz%XI, JJJ- b%Hz%YI, KKK- b%Hz%ZI))%is%PML).and.(iii <= SINPML_fullsize(iHz)%XE).and.(jjj <= SINPML_fullsize(iHz)%YE).and.(kkk <= SINPML_fullsize(iHz)%ZE)) then
                                              conta=conta+1
                                              jz=-100
                                              output( ii)%item( i)%Serialized%valor(Ntimeforvolumic,conta) = Jz
                                           endif
                                       endif !de los tags negativos
                                    endif !del if mapvtk
                                    !
                                 end do
                              end do
                           end do

                           !!!
                           if (field==mapvtk) then
                              INIT=.false.; geom=.false. ; asigna=.true.; magnetic=.true. ; electric=.false.
#ifdef CompileWithNodalSources
                              call nodalvtk(sgg,sggMiEx,sggMiEy,sggMiEz,sggMiHx,sggMiHy,sggMiHz,sggMtag, &
                                            b,init,geom,asigna,electric,magnetic,conta,i,ii,output,Ntimeforvolumic)
#endif                              
                           endif
                           !!!
                           !!!!!!!!!!!!esto dara problemas en los angulos y aristas donde porque ahi sacara la Bloque current en Hx!!!! 19/2/14
                        endif
                     endif
                     !!!!!!!!fin sondas corriente
#ifdef CompileWithNF2FF
                   case( FarField)
                     if (planewavecorr) then
                         Excor=Ex-Exvac; Eycor=Ey-Eyvac; Ezcor=Ez-Ezvac;
                         Hxcor=Hx-Hxvac; Hycor=Hy-Hyvac; Hzcor=Hz-Hzvac;
                         call UpdateFarField(ntime, b, Excor, Eycor, Ezcor,Hxcor,Hycor,Hzcor)
                     else
                         call UpdateFarField(ntime, b, Ex, Ey, Ez,Hx,Hy,Hz)
                     endif
                     
#endif
                  endselect
                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FREQMAIN
                  !!!!!!!!!!!!!!!!!!!!
                  !!!!!!!!!!!!!!!!!!!!
               elseif (SGG%Observation(ii)%FreqDomain) then
                  at=sgg%tiempo(ntime)
!!!!! permit scaling

                  do iff1=1,output(ii)%NumFreqs
                    output(ii)%auxExp_E(iff1)= sgg%dt * (1.0E0_RKIND, 0.0E0_RKIND) * Exp(mcpi2* output(ii)%Freq(iff1) * at )   !el dt deberia ser algun tipo de promedio pero no me complico permit scaling 211118
                    output(ii)%auxExp_H(iff1)= output(ii)%auxExp_E(iff1)    * Exp(mcpi2* output(ii)%Freq(iff1) * sgg%dt * 0.5_RKIND  )   

                  end do
!!!
                  SElect case (field)
!!!las freqdomain NUNCA ESTAN DONE
                   case( iMEC,iExC,iEyC,iEzC) !como los tengo que guardar todas las componentes cartesianas solo variara el output final en el .xdmf
 !                    if (at > sgg%OBSERVATION(ii)%FinalTime+sgg%dt/2.0_RKIND) sgg%OBSERVATION(ii)%Done=.true.
                     if (at >= sgg%OBSERVATION(ii)%InitialTime)     sgg%OBSERVATION(ii)%Begun=.true.
                     do KKK = k1, k2
                     if (mod(KKK,output(ii)%item(i)%Ztrancos)==0) then
                     k1t=int(kkk/output(ii)%item(i)%Ztrancos)
                        KKK_m = KKK - b%Ez%ZI
                        do JJJ = j1, j2
                        if (mod(jjj,output(ii)%item(i)%Ytrancos)==0) then
                        j1t=int(jjj/output(ii)%item(i)%Ytrancos)
                           JJJ_m = JJJ - b%Ey%YI
                           do III = i1, i2
                           if (mod(iii,output(ii)%item(i)%Xtrancos)==0) then
                           i1t=int(iii/output(ii)%item(i)%Xtrancos)
                              III_m = III - b%Ex%XI
                              do if1=1,output( ii)%NumFreqs
                                 !!!������
                                 output( ii)%item( i)%valor3DComplex(if1,1,i1t,j1t,k1t) = output( ii)%item( i)%valor3DComplex(if1,1,i1t,j1t,k1t)+  &
                                                       output( ii)%auxExp_E(if1) * Ex( III_m, JJJ_m, KKK_m)
                                 output( ii)%item( i)%valor3DComplex(if1,2,i1t,j1t,k1t) = output( ii)%item( i)%valor3DComplex(if1,2,i1t,j1t,k1t)+  &
                                                       output( ii)%auxExp_E(if1) * Ey( III_m, JJJ_m, KKK_m)
                                 output( ii)%item( i)%valor3DComplex(if1,3,i1t,j1t,k1t) = output( ii)%item( i)%valor3DComplex(if1,3,i1t,j1t,k1t)+  &
                                                       output( ii)%auxExp_E(if1) * Ez( III_m, JJJ_m, KKK_m)
                              enddo
                           endif
                           enddo
                        endif
                        enddo
                    endif
                    enddo
                   case( iMHC,iHxC,iHyC,iHzC)
!                     if (at > sgg%OBSERVATION(ii)%FinalTime+sgg%dt/2.0_RKIND) sgg%OBSERVATION(ii)%Done=.true.
                     if (at >= sgg%OBSERVATION(ii)%InitialTime)     sgg%OBSERVATION(ii)%Begun=.true.
                     do KKK = k1, k2
                     if (mod(KKK,output(ii)%item(i)%Ztrancos)==0) then
                     k1t=int(kkk/output(ii)%item(i)%Ztrancos)
                        KKK_m = KKK - b%Hz%ZI
                        do JJJ = j1, j2
                        if (mod(jjj,output(ii)%item(i)%Ytrancos)==0) then
                        j1t=int(jjj/output(ii)%item(i)%Ytrancos)
                           JJJ_m = JJJ - b%Hy%YI
                           do III = i1, i2
                           if (mod(iii,output(ii)%item(i)%Xtrancos)==0) then
                           i1t=int(iii/output(ii)%item(i)%Xtrancos)
                              III_m = III - b%Hx%XI
                              do if1=1,output( ii)%NumFreqs
                                 !!!������
                                 output( ii)%item( i)%valor3DComplex(if1,1,i1t,j1t,k1t) = output( ii)%item( i)%valor3DComplex(if1,1,i1t,j1t,k1t)+  &
                                                    output( ii)%auxExp_H(if1) * Hx( III_m, JJJ_m, KKK_m)
                                 output( ii)%item( i)%valor3DComplex(if1,2,i1t,j1t,k1t) = output( ii)%item( i)%valor3DComplex(if1,2,i1t,j1t,k1t)+  &
                                                    output( ii)%auxExp_H(if1) * Hy( III_m, JJJ_m, KKK_m)
                                 output( ii)%item( i)%valor3DComplex(if1,3,i1t,j1t,k1t) = output( ii)%item( i)%valor3DComplex(if1,3,i1t,j1t,k1t)+  &
                                                    output( ii)%auxExp_H(if1) * Hz( III_m, JJJ_m, KKK_m)
                              enddo
                           endif
                           enddo
                        endif
                        enddo
                    endif
                    enddo
                     !sondas corriente 20/2/14
                   case( iCur,iCurX,iCurY,iCurZ)
!                     if (at > sgg%OBSERVATION(ii)%FinalTime+sgg%dt/2.0_RKIND) sgg%OBSERVATION(ii)%Done=.true.
                     if (at >= sgg%OBSERVATION(ii)%InitialTime)     sgg%OBSERVATION(ii)%Begun=.true.
                     conta=0 !en el mismo orden que se allocatearon
                     do KKK = k1, k2
                        do JJJ = j1, j2
                           do III = i1, i2
                              !saca bul current a lo largo del edgje con las sondas icur
                              if ((sgg%med(sggMiEx(III - b%Ex%XI, JJJ - b%Ex%YI, KKK - b%Ex%ZI))%Is%ThinWire).and.(iii <= SINPML_fullsize(iEx)%XE).and.(jjj <= SINPML_fullsize(iEx)%YE).and.(kkk <= SINPML_fullsize(iEx)%ZE)) then
                                 conta=conta+1
                                 Jx= dyh(JJJ - b%Hy%YI) * (- Hy( III - b%Hy%XI, JJJ - b%Hy%YI, KKK - b%Ex%ZI)+Hy( III - b%Hy%XI  , JJJ - b%Hy%YI  , KKK - b%Hy%ZI-1))  + &
                                     dzh(KKK - b%Hz%ZI) * (  Hz( III - b%Hz%XI, JJJ - b%Hz%YI, KKK - b%Hz%ZI)-Hz( III - b%Hz%XI  , JJJ - b%Hz%YI-1, KKK - b%Hz%ZI  ))
                                 do if1=1,output( ii)%NumFreqs  
                                    output( ii)%item( i)%Serialized%valorComplex_x(if1,conta) = output( ii)%item( i)%Serialized%valorComplex_x(if1,conta) +  &
                                                       output( ii)%auxExp_H(if1) * Jx      
                                 end do
                              endif
                              if ((sgg%med(sggMiEy(III - b%Ey%XI, JJJ - b%Ey%YI, KKK - b%Ey%ZI))%Is%ThinWire).and.(iii <= SINPML_fullsize(iEy)%XE).and.(jjj <= SINPML_fullsize(iEy)%YE).and.(kkk <= SINPML_fullsize(iEy)%ZE)) then
                                 conta=conta+1
                                 Jy=dxh(III - b%Hx%XI) * (  Hx( III - b%Hx%XI, JJJ - b%Hx%YI, KKK - b%Hx%ZI)-Hx( III - b%Hx%XI  , JJJ - b%Hx%YI  , KKK - b%Hx%ZI-1))  + &
                                    dzh(KKK - b%Hz%ZI) * ( -Hz( III - b%Hz%XI, JJJ - b%Hz%YI, KKK - b%Hz%ZI)+Hz( III - b%Hz%XI-1, JJJ - b%Hz%YI  , KKK - b%Hz%ZI  ))
                                 do if1=1,output( ii)%NumFreqs                                   
                                    output( ii)%item( i)%Serialized%valorComplex_y(if1,conta) = output( ii)%item( i)%Serialized%valorComplex_y(if1,conta) +  &
                                                       output( ii)%auxExp_H(if1) * Jy  
                                 end do
                              endif
                              if ((sgg%med(sggMiEz(III - b%Ez%XI, JJJ - b%Ez%YI, KKK - b%Ez%ZI))%Is%ThinWire).and.(iii <= SINPML_fullsize(iEz)%XE).and.(jjj <= SINPML_fullsize(iEz)%YE).and.(kkk <= SINPML_fullsize(iEz)%ZE)) then
                                 conta=conta+1
                                 Jz=dyh(JJJ - b%Hy%YI) * (  Hy( III - b%Hy%XI, JJJ - b%Hy%YI, KKK - b%Hy%ZI) - Hy( III - b%Hy%XI-1, JJJ - b%Hy%YI  , KKK - b%Hy%ZI  )) + &
                                    dxh(III - b%Hx%XI) * ( -Hx( III - b%Hx%XI, JJJ - b%Hx%YI, KKK - b%Hx%ZI) + Hx( III - b%Hx%XI  , JJJ - b%Hx%YI-1, KKK - b%Hx%ZI  ))
                                 do if1=1,output( ii)%NumFreqs
                                    output( ii)%item( i)%Serialized%valorComplex_z(if1,conta) = output( ii)%item( i)%Serialized%valorComplex_z(if1,conta) +  &
                                                       output( ii)%auxExp_H(if1) * Jz
                                 end do
                              endif
                              if (((sgg%med(sggMiEx(III - b%Ex%XI, JJJ - b%Ex%YI, KKK - b%Ex%ZI))%Is%Line.and.sgg%med(sggMiEx(III - b%Ex%XI, JJJ - b%Ex%YI, KKK - b%Ex%ZI))%Is%PEC)).and.(iii <= SINPML_fullsize(iEx)%XE).and.(jjj <= SINPML_fullsize(iEx)%YE).and.(kkk <= SINPML_fullsize(iEx)%ZE)) then
                                 if ((.not.sgg%med(sggMiHy(III - b%Hy%XI  , JJJ - b%Hy%YI  , KKK - b%Hy%ZI  ))%Is%PEC).and. &
                                 (.not.sgg%med(sggMiHy(III - b%Hy%XI  , JJJ - b%Hy%YI  , KKK - b%Hy%ZI-1))%Is%PEC).and. &
                                 (.not.sgg%med(sggMiHz(III - b%Hz%XI  , JJJ - b%Hz%YI  , KKK - b%Hz%ZI  ))%Is%PEC).and. &
                                 (.not.sgg%med(sggMiHz(III - b%Hz%XI  , JJJ - b%Hz%YI-1, KKK - b%Hz%ZI  ))%Is%PEC)) then
                                    conta=conta+1
                                    Jx= dyh(JJJ - b%Hy%YI) * (- Hy( III - b%Hy%XI, JJJ - b%Hy%YI, KKK - b%Ex%ZI)+Hy( III - b%Hy%XI  , JJJ - b%Hy%YI  , KKK - b%Hy%ZI-1))    + &
                                        dzh(KKK - b%Hz%ZI) * (  Hz( III - b%Hz%XI, JJJ - b%Hz%YI, KKK - b%Hz%ZI)-Hz( III - b%Hz%XI  , JJJ - b%Hz%YI-1, KKK - b%Hz%ZI  ))
                                    do if1=1,output( ii)%NumFreqs
                                       output( ii)%item( i)%Serialized%valorComplex_x(if1,conta) = output( ii)%item( i)%Serialized%valorComplex_x(if1,conta) +  &
                                                          output( ii)%auxExp_H(if1) * Jx    
                                    end do
                                 endif
                              endif
                              if (((sgg%med(sggMiEy(III - b%Ey%XI, JJJ - b%Ey%YI, KKK - b%Ey%ZI))%Is%Line.and.sgg%med(sggMiEy(III - b%Ey%XI, JJJ - b%Ey%YI, KKK - b%Ey%ZI))%Is%PEC)).and.(iii <= SINPML_fullsize(iEy)%XE).and.(jjj <= SINPML_fullsize(iEy)%YE).and.(kkk <= SINPML_fullsize(iEy)%ZE)) then
                                 if ((.not.sgg%med(sggMiHz(III - b%Hz%XI  , JJJ - b%Hz%YI  , KKK - b%Hz%ZI  ))%Is%PEC).and. &
                                 (.not.sgg%med(sggMiHz(III - b%Hz%XI-1, JJJ - b%Hz%YI  , KKK - b%Hz%ZI  ))%Is%PEC).and. &
                                 (.not.sgg%med(sggMiHx(III - b%Hx%XI  , JJJ - b%Hx%YI  , KKK - b%Hx%ZI  ))%Is%PEC).and. &
                                 (.not.sgg%med(sggMiHx(III - b%Hx%XI  , JJJ - b%Hx%YI  , KKK - b%Hx%ZI-1))%Is%PEC)) then
                                    conta=conta+1
                                    Jy=dxh(III - b%Hx%XI) * (  Hx( III - b%Hx%XI, JJJ - b%Hx%YI, KKK - b%Hx%ZI)-Hx( III - b%Hx%XI  , JJJ - b%Hx%YI  , KKK - b%Hx%ZI-1)) + &
                                       dzh(KKK - b%Hz%ZI) * ( -Hz( III - b%Hz%XI, JJJ - b%Hz%YI, KKK - b%Hz%ZI)+Hz( III - b%Hz%XI-1, JJJ - b%Hz%YI  , KKK - b%Hz%ZI  ))
                                    do if1=1,output( ii)%NumFreqs
                                       output( ii)%item( i)%Serialized%valorComplex_y(if1,conta) = output( ii)%item( i)%Serialized%valorComplex_y(if1,conta) +  &
                                                          output( ii)%auxExp_H(if1) * Jy        
                                       
                                    end do
                                 endif
                              endif
                              if (((sgg%med(sggMiEz(III - b%Ez%XI, JJJ - b%Ez%YI, KKK - b%Ez%ZI))%Is%Line.and.sgg%med(sggMiEz(III - b%Ez%XI, JJJ - b%Ez%YI, KKK - b%Ez%ZI))%Is%PEC)).and.(iii <= SINPML_fullsize(iEz)%XE).and.(jjj <= SINPML_fullsize(iEz)%YE).and.(kkk <= SINPML_fullsize(iEz)%ZE)) then
                                 if ((.not.sgg%med(sggMiHx(III - b%Hx%XI  , JJJ - b%Hx%YI  , KKK - b%Hx%ZI  ))%Is%PEC).and. &
                                 (.not.sgg%med(sggMiHx(III - b%Hx%XI  , JJJ - b%Hx%YI-1, KKK - b%Hx%ZI  ))%Is%PEC).and. &
                                 (.not.sgg%med(sggMiHy(III - b%Hy%XI  , JJJ - b%Hy%YI  , KKK - b%Hy%ZI  ))%Is%PEC).and. &
                                 (.not.sgg%med(sggMiHy(III - b%Hy%XI-1, JJJ - b%Hy%YI  , KKK - b%Hy%ZI  ))%Is%PEC)) then
                                    conta=conta+1
                                    Jz=dyh(JJJ - b%Hy%YI) * (  Hy( III - b%Hy%XI, JJJ - b%Hy%YI, KKK - b%Hy%ZI) - Hy( III - b%Hy%XI-1, JJJ - b%Hy%YI  , KKK - b%Hy%ZI  ))  + &
                                       dxh(III - b%Hx%XI) * ( -Hx( III - b%Hx%XI, JJJ - b%Hx%YI, KKK - b%Hx%ZI) + Hx( III - b%Hx%XI  , JJJ - b%Hx%YI-1, KKK - b%Hx%ZI  ))
                                    do if1=1,output( ii)%NumFreqs     
                                       output( ii)%item( i)%Serialized%valorComplex_z(if1,conta) = output( ii)%item( i)%Serialized%valorComplex_z(if1,conta) +  &
                                                          output( ii)%auxExp_H(if1) * Jz
                                    end do
                                 endif
                              endif
                           end do
                        end do
                     end do
                     do KKK = k1, k2
                        do JJJ = j1, j2
                           do III = i1, i2
                              !saca current en surfaces 0124
                              if (((sgg%med(sggMiHx(III -b%Hx%XI, JJJ- b%Hx%YI, KKK- b%Hx%ZI))%Is%PEC).or. &
                              (sgg%med(sggMiHx(III -b%Hx%XI, JJJ- b%Hx%YI, KKK- b%Hx%ZI))%Is%Surface).or.  &
                              (field==icurX)).and.(iii <= SINPML_fullsize(iHx)%XE).and.(jjj <= SINPML_fullsize(iHx)%YE).and.(kkk <= SINPML_fullsize(iHx)%ZE)) then
                                 conta=conta+1
                                 Jy=(dzh(KKK - b%Hz%ZI) * Hz( III - b%Hz%XI-1, JJJ - b%Hz%YI  , KKK - b%Hz%ZI  ) + dzh(KKK - b%Hz%ZI+1) *Hz( III - b%Hz%XI-1, JJJ - b%Hz%YI  , KKK - b%Hz%ZI+1) )/2.0_RKIND -  &
                                    (dzh(KKK - b%Hz%ZI) * Hz( III - b%Hz%XI  , JJJ - b%Hz%YI  , KKK - b%Hz%ZI  ) + dzh(KKK - b%Hz%ZI+1) *Hz( III - b%Hz%XI  , JJJ - b%Hz%YI  , KKK - b%Hz%ZI+1) )/2.0_RKIND +  &
                                     dxh(III - b%Hx%XI)*( Hx( III - b%Hx%XI  , JJJ - b%Hx%YI  , KKK - b%Hx%ZI+1) -                       Hx( III - b%Hx%XI  , JJJ - b%Hx%YI  , KKK - b%Hx%ZI-1) )/2.0_RKIND
                                 !el Hx al promediarlo con el suyo (i,j,k) a ambos lados pierde su componente y solo quedan las adyancentes
                                 Jz=(dyh(JJJ - b%Hy%YI) * Hy( III - b%Hy%XI  , JJJ - b%Hy%YI  , KKK - b%Hy%ZI  ) + dyh(JJJ - b%Hy%YI+1) *Hy( III - b%Hy%XI  , JJJ - b%Hy%YI+1, KKK - b%Hy%ZI  ) )/2.0_RKIND -  &
                                    (dyh(JJJ - b%Hy%YI) * Hy( III - b%Hy%XI-1, JJJ - b%Hy%YI  , KKK - b%Hy%ZI  ) + dyh(JJJ - b%Hy%YI+1) *Hy( III - b%Hy%XI-1, JJJ - b%Hy%YI+1, KKK - b%Hy%ZI  ) )/2.0_RKIND +  &
                                     dxh(III - b%Hx%XI)*( Hx( III - b%Hx%XI  , JJJ - b%Hx%YI-1, KKK - b%Hx%ZI  ) -                       Hx( III - b%Hx%XI  , JJJ - b%Hx%YI+1, KKK - b%Hx%ZI  ) )/2.0_RKIND
                                 do if1=1,output( ii)%NumFreqs
                                    output( ii)%item( i)%Serialized%valorComplex_y(if1,conta) = output( ii)%item( i)%Serialized%valorComplex_y(if1,conta) +  &
                                                       output( ii)%auxExp_H(if1) * Jy
                                    output( ii)%item( i)%Serialized%valorComplex_z(if1,conta) = output( ii)%item( i)%Serialized%valorComplex_z(if1,conta) +  &
                                                       output( ii)%auxExp_H(if1) * Jz
                                 end do
                              endif
                              if (((sgg%med(sggMiHy(III- b%Hy%XI, JJJ- b%Hy%YI, KKK- b%Hy%ZI))%Is%PEC).or. &
                              (sgg%med(sggMiHy(III- b%Hy%XI, JJJ- b%Hy%YI, KKK- b%Hy%ZI))%Is%Surface).or. &
                              (field==icurY)).and.(iii <= SINPML_fullsize(iHy)%XE).and.(jjj <= SINPML_fullsize(iHy)%YE).and.(kkk <= SINPML_fullsize(iHy)%ZE)) then
                                 conta=conta+1
                                 Jz=(dxh(III - b%Hx%XI) * Hx( III - b%Hx%XI  , JJJ - b%Hx%YI-1, KKK - b%Hx%ZI  ) + dxh(III - b%Hx%XI+1) *Hx( III - b%Hx%XI+1, JJJ - b%Hx%YI-1, KKK - b%Hx%ZI  ) )/2.0_RKIND -  &
                                    (dxh(III - b%Hx%XI) * Hx( III - b%Hx%XI  , JJJ - b%Hx%YI  , KKK - b%Hx%ZI  ) + dxh(III - b%Hx%XI+1) *Hx( III - b%Hx%XI+1, JJJ - b%Hx%YI  , KKK - b%Hx%ZI  ) )/2.0_RKIND +  &
                                     dyh(JJJ - b%Hy%YI)*( Hy( III - b%Hy%XI+1, JJJ - b%Hy%YI  , KKK - b%Hy%ZI  ) -                       Hy( III - b%Hy%XI-1, JJJ - b%Hy%YI  , KKK - b%Hy%ZI  ) )/2.0_RKIND
                                 !
                                 Jx=(dzh(KKK - b%Hz%ZI) * Hz( III - b%Hz%XI  , JJJ - b%Hz%YI  , KKK - b%Hz%ZI  ) + dzh(KKK - b%Hz%ZI+1) *Hz( III - b%Hz%XI  , JJJ - b%Hz%YI  , KKK - b%Hz%ZI+1) )/2.0_RKIND -  &
                                    (dzh(KKK - b%Hz%ZI) * Hz( III - b%Hz%XI  , JJJ - b%Hz%YI-1, KKK - b%Hz%ZI  ) + dzh(KKK - b%Hz%ZI+1) *Hz( III - b%Hz%XI  , JJJ - b%Hz%YI-1, KKK - b%Hz%ZI+1) )/2.0_RKIND +  &
                                     dyh(JJJ - b%Hy%YI)*( Hy( III - b%Hy%XI  , JJJ - b%Hy%YI  , KKK - b%Hy%ZI-1) -                       Hy( III - b%Hy%XI  , JJJ - b%Hy%YI  , KKK - b%Hy%ZI+1) )/2.0_RKIND
                                 !
                                 do if1=1,output( ii)%NumFreqs
                                    output( ii)%item( i)%Serialized%valorComplex_z(if1,conta) = output( ii)%item( i)%Serialized%valorComplex_z(if1,conta) +  &
                                                       output( ii)%auxExp_H(if1) * Jz
                                    output( ii)%item( i)%Serialized%valorComplex_x(if1,conta) = output( ii)%item( i)%Serialized%valorComplex_x(if1,conta) +  &
                                                       output( ii)%auxExp_H(if1) * Jx
                                 end do
                              endif
                              if (((sgg%med(sggMiHz(III- b%Hz%XI, JJJ- b%Hz%YI, KKK- b%Hz%ZI))%Is%PEC).or. &
                              (sgg%med(sggMiHz(III- b%Hz%XI, JJJ- b%Hz%YI, KKK- b%Hz%ZI))%Is%Surface).or. &
                              (field==icurZ)).and.(iii <= SINPML_fullsize(iHz)%XE).and.(jjj <= SINPML_fullsize(iHz)%YE).and.(kkk <= SINPML_fullsize(iHz)%ZE)) then
                                 conta=conta+1
                                 Jx=(dyh(JJJ - b%Hy%YI) * Hy( III - b%Hy%XI  , JJJ - b%Hy%YI  , KKK - b%Hy%ZI-1) + dyh(JJJ - b%Hy%YI+1) *Hy( III - b%Hy%XI  , JJJ - b%Hy%YI+1, KKK - b%Hy%ZI-1) )/2.0_RKIND -  &
                                    (dyh(JJJ - b%Hy%YI) * Hy( III - b%Hy%XI  , JJJ - b%Hy%YI  , KKK - b%Hy%ZI  ) + dyh(JJJ - b%Hy%YI+1) *Hy( III - b%Hy%XI  , JJJ - b%Hy%YI+1, KKK - b%Hy%ZI  ) )/2.0_RKIND +  &
                                     dzh(KKK - b%Hz%ZI)*( Hz( III - b%Hz%XI  , JJJ - b%Hz%YI+1, KKK - b%Hz%ZI  ) -                       Hz( III - b%Hz%XI  , JJJ - b%Hz%YI-1, KKK - b%Hz%ZI  ) )/2.0_RKIND
                                 !
                                 Jy=(dxh(III - b%Hx%XI) * Hx( III - b%Hx%XI  , JJJ - b%Hx%YI  , KKK - b%Hx%ZI  ) + dxh(III - b%Hx%XI+1) *Hx( III - b%Hx%XI+1, JJJ - b%Hx%YI  , KKK - b%Hx%ZI  ) )/2.0_RKIND -  &
                                    (dxh(III - b%Hx%XI) * Hx( III - b%Hx%XI  , JJJ - b%Hx%YI  , KKK - b%Hx%ZI-1) + dxh(III - b%Hx%XI+1) *Hx( III - b%Hx%XI+1, JJJ - b%Hx%YI  , KKK - b%Hx%ZI-1) )/2.0_RKIND +  &
                                     dzh(KKK - b%Hz%ZI)*( Hz( III - b%Hz%XI-1, JJJ - b%Hz%YI  , KKK - b%Hz%ZI  ) -                       Hz( III - b%Hz%XI+1, JJJ - b%Hz%YI  , KKK - b%Hz%ZI  ) )/2.0_RKIND
                                 do if1=1,output( ii)%NumFreqs  
                                    output( ii)%item( i)%Serialized%valorComplex_x(if1,conta) = output( ii)%item( i)%Serialized%valorComplex_x(if1,conta) +  &
                                                       output(ii)%auxExp_H(if1) * Jx
                                    output( ii)%item( i)%Serialized%valorComplex_y(if1,conta) = output( ii)%item( i)%Serialized%valorComplex_y(if1,conta) +  &
                                                       output(ii)%auxExp_H(if1) * Jy
                                 end do
                              endif
                              !nunca va a haver mapvtk en frequencia!!!
                              !
                           end do
                        end do
                     end do
                     !!!!!!!!!!!!esto dara problemas en los angulos y aristas donde porque ahi sacara la Bloque current en Hx!!!! 19/2/14

                     !!!!!!!!fin sondas corriente
                  end select

               endif !time domain
            endif !del nothing
         enddo   loop_obser
      enddo
      !---------------------------> acaba UpdateObservation <-----------------------------------------
      return
   endsubroutine UpdateObservation


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Flushes the observed magnitudes to disk
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine FlushObservationFiles(sgg,nInit,FinalInstant,layoutnumber,size, dxe,dye,dze,dxh,dyh,dzh,b,singlefilewrite,facesNF2FF,flushff)
      USE ILUMINA !is needed to also calculate the incident field in the observed points
      !
      type (nf2ff_t) :: facesNF2FF
      !!!
      type( bounds_t), intent( IN)      ::  b
      real (kind = RKIND), dimension( 0 :  b%dxe%NX-1), intent( IN)  ::  dxe
      real (kind = RKIND), dimension( 0 :  b%dye%NY-1), intent( IN)  ::  dye
      real (kind = RKIND), dimension( 0 :  b%dze%NZ-1), intent( IN)  ::  dze
      !--->
      real (kind = RKIND), dimension( 0 :  b%dxh%NX-1), intent( IN)  ::  dxh
      real (kind = RKIND), dimension( 0 :  b%dyh%NY-1), intent( IN)  ::  dyh
      real (kind = RKIND), dimension( 0 :  b%dzh%NZ-1), intent( IN)  ::  dzh
      !
      type (SGGFDTDINFO), intent(IN)         ::  sgg
      integer (kind=4), intent(in) :: layoutnumber,size
      integer (kind=4)  ::  nInit,FinalInstant,unidad,compo,conta
      integer (kind=4)  ::  i,field,N,ii,i1,j1,k1,Ntimeforvolumic,dummy_jjj,i1t,j1t,k1t,i0t
      logical  ::  incident,singlefilewrite,flushff,ISyaopen
      REAL (KIND=RKIND_tiempo)  ::  at

      logical :: ok
      logical :: called_fromobservation,dummy_logical
      integer :: my_iostat

      character (LEN=BUFSIZE)  ::  whoami
      !!!
      write(whoami,'(a,i5,a,i5,a)') '(',layoutnumber+1,'/',size,') '
      called_fromobservation=.true.
      !!!ojo dummy_logical lo dejo que incid( lo fije a still_planewave_time sin tocarlo aqui para no afectar al principal 210419
      dummy_jjj=1  !no es preciso fijarlo, porque a incid( se le pasa  called_fromobservation
      !Write also the incident fields in case there are plane waves (useful in SE calculations)
      incident=.false.
      if (sgg%NumPlaneWaves>=1) incident=.true.
      do ii=1,sgg%NumberRequest
         loop_obb :  do i=1,sgg%Observation(ii)%nP
            field=sgg%observation(ii)%P(i)%what
            if (SGG%Observation(ii)%TimeDomain) then
               select case(field)
                case (iEx,iEy,iEz,iJx,iJy,iJz,iHx,iHy,iHz)
                  I1=sgg%OBSERVATION(ii)%P(i)%XI
                  J1=sgg%OBSERVATION(ii)%P(i)%YI
                  K1=sgg%OBSERVATION(ii)%P(i)%ZI
                  !
                  !                nInit=max(nint(0.4999999+sgg%OBSERVATION(ii)%InitialTime/sgg%dt),nInit)
                  DO N=nInit,FinalInstant !!! bug octubre'14 mur1.nfde. Quitado --->>>,output(ii)%Trancos  !save only for the requested data
                     if (mod(n,output(ii)%Trancos)==0) then  !save only for the requested data
                        select case(field)
                         case (iHx,iHy,iHz)
                           at=sgg%tiempo(N)+sgg%dt/2.0_RKIND_tiempo
                         case(iEx,iEy,iEz,iJx,iJy,iJz)
                           at=sgg%tiempo(N)
                        end select
                        if ( ((at >= sgg%OBSERVATION(ii)%InitialTime).and.(at <= sgg%OBSERVATION(ii)%FinalTime+sgg%dt/2.0_RKIND)).or. &
                        output(ii)%saveall) then

                           select case(field)
                            case (iEx,iEy,iEz)
                              !
                              if (singlefilewrite) then
                                 unidad=output(ii)%item(i)%unitmaster
                                 if (incident) then
                                    WRITE (unidad) output(ii)%item(i)%unit,at, output(ii)%item(i)%valor(n-nInit), &
                                    Incid(sgg,dummy_jjj,field,real(at+0.0_RKIND * sgg%dt,RKIND),i1,j1,k1,dummy_logical,called_fromobservation) !quitado el 3 de ORIGINAL sync para pscale bien sincronizado
                                    !by hand para clavarlo
                                 else
                                    WRITE (unidad) output(ii)%item(i)%unit,at, output(ii)%item(i)%valor(n-nInit)
                                 endif
                              else
                                 unidad=output(ii)%item(i)%unit
                                 if (incident) then
#ifdef CompileWithReal16
                                    WRITE (unidad,*) at, output(ii)%item(i)%valor(n-nInit), &
                                    Incid(sgg,dummy_jjj,field,at*1.0_RKIND+0.0_RKIND * sgg%dt,i1,j1,k1,dummy_logical,called_fromobservation) !quitado el 3 de ORIGINAL sync para pscale bien sincronizado
#else
#ifdef CompileWithmMiguelStandaloneObservation
                                    WRITE (unidad,*) at, output(ii)%item(i)%valor(n-nInit), &
                                    Incid(sgg,dummy_jjj,field,real(at+0.0_RKIND * sgg%dt,RKIND),i1,j1,k1,dummy_logical,called_fromobservation) !quitado el 3 de ORIGINAL sync para pscale bien sincronizado
#else
                                    WRITE (unidad,fmt) at, output(ii)%item(i)%valor(n-nInit), &
                                    Incid(sgg,dummy_jjj,field,real(at+0.0_RKIND * sgg%dt,RKIND),i1,j1,k1,dummy_logical,called_fromobservation) !quitado el 3 de ORIGINAL sync para pscale bien sincronizado
                                    !by hand para clavarlo
#endif
#endif
                                 else
                                    WRITE (unidad,fmt) at, output(ii)%item(i)%valor(n-nInit)
                                 endif
                              endif
                              !
                            case (iHx,iHy,iHz)
                              !
                              if (singlefilewrite) then
                                 unidad=output(ii)%item(i)%unitmaster
                                 if (incident) then
                                    WRITE (unidad) output(ii)%item(i)%unit,at-sgg%dt/2.0_RKIND_tiempo, output(ii)%item(i)%valor(n-nInit), &  !SOLO A EFECTOS DE SALIDA EN FICHERO CHAPUZ SGG MAIL OLD 070916 �
                                    Incid(sgg,dummy_jjj,field,real(at+0.0_RKIND * sgg%dt,RKIND),i1,j1,k1,dummy_logical,called_fromobservation) !quitado el 3 de ORIGINAL sync para pscale bien sincronizado
                                    !by hand para clavarlo
                                 else
                                    WRITE (unidad) output(ii)%item(i)%unit,at-sgg%dt/2.0_RKIND_tiempo, output(ii)%item(i)%valor(n-nInit) !SOLO A EFECTOS DE SALIDA EN FICHERO CHAPUZ SGG MAIL OLD 070916 �
                                 endif
                              else
                                 unidad=output(ii)%item(i)%unit
                                 if (incident) then
#ifdef CompileWithReal16
                                    WRITE (unidad,*) at-sgg%dt/2.0_RKIND_tiempo, output(ii)%item(i)%valor(n-nInit), &  !SOLO A EFECTOS DE SALIDA EN FICHERO CHAPUZ SGG MAIL OLD 070916 �
                                    Incid(sgg,dummy_jjj,field,real(at+0.0_RKIND * sgg%dt,RKIND),i1,j1,k1,dummy_logical,called_fromobservation) !quitado el 3 de ORIGINAL sync para pscale bien sincronizado
#else
#ifdef CompileWithmMiguelStandaloneObservation
                                    WRITE (unidad,*) at-sgg%dt/2.0_RKIND_tiempo, output(ii)%item(i)%valor(n-nInit), &  !SOLO A EFECTOS DE SALIDA EN FICHERO CHAPUZ SGG MAIL OLD 070916 �
                                    Incid(sgg,dummy_jjj,field,real(at+0.0_RKIND * sgg%dt,RKIND),i1,j1,k1,dummy_logical,called_fromobservation) !quitado el 3 de ORIGINAL sync para pscale bien sincronizado
#else
                                    WRITE (unidad,fmt) at-sgg%dt/2.0_RKIND_tiempo, output(ii)%item(i)%valor(n-nInit), &  !SOLO A EFECTOS DE SALIDA EN FICHERO CHAPUZ SGG MAIL OLD 070916 �
                                    Incid(sgg,dummy_jjj,field,real(at+0.0_RKIND * sgg%dt,RKIND),i1,j1,k1,dummy_logical,called_fromobservation) !quitado el 3 de ORIGINAL sync para pscale bien sincronizado
                                    !by hand para clavarlo
#endif
#endif
                                 else
                                    WRITE (unidad,fmt) at-sgg%dt/2.0_RKIND_tiempo, output(ii)%item(i)%valor(n-nInit) !SOLO A EFECTOS DE SALIDA EN FICHERO CHAPUZ SGG MAIL OLD 070916 �
                                 endif
                              endif
                              !

#ifdef CompileWithWires
                            case(iJx,iJy,iJz)
                              if (singlefilewrite) then
                                 unidad=output(ii)%item(i)%unitmaster
                                 WRITE (unidad) output(ii)%item(i)%unit,at, &
                                                output(ii)%item(i)%valor(n-nInit), &
                                                output(ii)%item(i)%valor2(n-nInit), & !saco el valor2 -e*dl
                                                output(ii)%item(i)%valor3(n-nInit), & ! VpluS
                                                output(ii)%item(i)%valor4(n-nInit), & ! Vminus
                                                output(ii)%item(i)%valor5(n-nInit) ! vplus-vminus
                              else
                                 unidad=output(ii)%item(i)%unit
                                 WRITE (unidad,fmt) at, output(ii)%item(i)%valor(n-nInit), &
                                                        output(ii)%item(i)%valor2(n-nInit) , & !saco el valor2 -e*dl
                                                        output(ii)%item(i)%valor3(n-nInit) , & ! VPLUS 
                                                        output(ii)%item(i)%valor4(n-nInit) , & ! Vminus 
                                                        output(ii)%item(i)%valor5(n-nInit) ! vplus-vminus
                              endif
#endif
                           end select
                        endif
                     endif
                  END DO
                  if (singlefilewrite.and.((field==iEx).or.(field==iEy).or.(field==iEz).or. &
                  (field==iVx).or.(field==iVy).or.(field==iVz).or. &
                  (field==iJx).or.(field==iJy).or.(field==iJz).or. &
                  (field==iHx).or.(field==iHy).or.(field==iHz))) then
                     unidad=output(ii)%item(i)%unitmaster
                  else
                     unidad=output(ii)%item(i)%unit
                  endif
                  call flush(unidad)
                  !
                case (iBloqueMx,iBloqueMz,iBloqueMy,iBloqueJx,iBloqueJz,iBloqueJy)
#ifdef CompileWithMPI
                  if ((field==iBloqueMx).or.(field==iBloqueMy).or.(field==iBloqueJx).or.(field==iBloqueJy)) then
                     if (output(ii)%item(i)%MPISubComm /= -1) then !solo si alguien tiene que hacerlo
                        valores(0 : BuffObse)=output(ii)%item(i)%valor(0 : BuffObse)
                        call MPIupdateBloques(layoutnumber,valores,newvalores, &
                        output(ii)%item(i)%MPISubComm)
                        output(ii)%item(i)%valor(0 : BuffObse)=newvalores(0 : BuffObse)
                     endif
                  endif

                  if ((layoutnumber == output(ii)%item(i)%MPIRoot).or. &
                  (field==iBloqueJz).or.(field==iBloqueMz)) then !only the master
#endif
                  !                    nInit=max(nint(0.4999999+sgg%OBSERVATION(ii)%InitialTime/sgg%dt),nInit)
                  DO N=nInit,FinalInstant
                     if (mod(n,output(ii)%Trancos)==0) then  !save only for the requested data
                        select case(field)
                         case (iBloqueMx,iBloqueMz,iBloqueMy)
                           at=sgg%tiempo(N)+sgg%dt/2.0_RKIND_tiempo
                         case (iBloqueJx,iBloqueJz,iBloqueJy)
                           at=sgg%tiempo(N)
                        end select
                        if ( ((at >= sgg%OBSERVATION(ii)%InitialTime).and.(at <= sgg%OBSERVATION(ii)%FinalTime+sgg%dt/2.0_RKIND)).or. &
                        output(ii)%saveall) then
                        select case(field)
                         case (iBloqueMx,iBloqueMz,iBloqueMy)
#ifdef CompileWithReal16
                           WRITE (output(ii)%item(i)%unit,*) at-sgg%dt/2.0_RKIND_tiempo, output(ii)%item(i)%valor(n-nInit) !SOLO A EFECTOS DE SALIDA EN FICHERO CHAPUZ SGG MAIL OLD 070916 �
                         case (iBloqueJx,iBloqueJz,iBloqueJy)
                           WRITE (output(ii)%item(i)%unit,*) at, output(ii)%item(i)%valor(n-nInit)
#else
#ifdef CompileWithmMiguelStandaloneObservation
                           WRITE (output(ii)%item(i)%unit,*) at-sgg%dt/2.0_RKIND_tiempo, output(ii)%item(i)%valor(n-nInit) !SOLO A EFECTOS DE SALIDA EN FICHERO CHAPUZ SGG MAIL OLD 070916 �
                         case (iBloqueJx,iBloqueJz,iBloqueJy)
                           WRITE (output(ii)%item(i)%unit,*) at, output(ii)%item(i)%valor(n-nInit)
#else
                           WRITE (output(ii)%item(i)%unit,fmt) at-sgg%dt/2.0_RKIND_tiempo, output(ii)%item(i)%valor(n-nInit) !SOLO A EFECTOS DE SALIDA EN FICHERO CHAPUZ SGG MAIL OLD 070916 �
                         case (iBloqueJx,iBloqueJz,iBloqueJy)
                           WRITE (output(ii)%item(i)%unit,fmt) at, output(ii)%item(i)%valor(n-nInit)
#endif
#endif
                        end select
                        endif
                     endif
                  END DO
                  call flush(output(ii)%item(i)%unit)
                  !--->

                  !--->
#ifdef CompileWithMPI
               endif
#endif
#ifdef CompileWithNF2FF
                CASe (FarField) !no emplear tiempo calculando rcs por el camino solo al final
                  at=sgg%tiempo(FinalInstant)
                  if (flushFF) call FlushFarfield(layoutnumber,size,  b, dxe, dye, dze, dxh, dyh, dzh,facesNF2FF,at)
#endif
                case (iMHC,iHxC,iHyC,iHzC,iMEC,iExC,iEyC,iEzC,icur,iCurX,iCurY,iCurZ,mapvtk)
                  DO N=nInit,FinalInstant
                     at=sgg%tiempo(N)
                     Ntimeforvolumic=N !!!-nint(0.4999999+sgg%OBSERVATION(ii)%InitialTime/sgg%dt)
                     if (mod(Ntimeforvolumic,output(ii)%Trancos)==0) then
                        Ntimeforvolumic=Ntimeforvolumic/output(ii)%Trancos
                        if ( ((at >= sgg%OBSERVATION(ii)%InitialTime).and.(at <= sgg%OBSERVATION(ii)%FinalTime+sgg%dt/2.0_RKIND))) then
                           !assumo que todos son electricos o magneticos en una probe Volumic para calcular el tiempo !logico
                           output(ii)%TimesWritten=output(ii)%TimesWritten+1
                           select case(field)
                            case (iMHC,iHxC,iHyC,iHzC,iMEC,iExC,iEyC,iEzC)
                              write(output(ii)%item(i)%unit) at
                              do k1t=output(ii)%item(i)%ZItrancos,output(ii)%item(i)%ZEtrancos
                                 do j1t=output(ii)%item(i)%YItrancos,output(ii)%item(i)%YEtrancos
                                    write(output(ii)%item(i)%unit) (output(ii)%item(i)%valor3D(Ntimeforvolumic,i1t,j1t,k1t), &
                                    &                                i1t=output(ii)%item(i)%XItrancos,output(ii)%item(i)%XEtrancos)
                                 end do
                              end do
                            case(icur,iCurX,iCurY,iCurZ,mapvtk)
                              write(output(ii)%item(i)%unit) at
                              if (output(ii)%item(i)%columnas /=0) then
                                  do i1=1,output(ii)%item(i)%columnas
                                    write(output(ii)%item(i)%unit) output(ii)%item(i)%Serialized%valor(Ntimeforvolumic,i1), &
                                                                   output(ii)%item(i)%Serialized%valor_x(Ntimeforvolumic,i1), &
                                                                   output(ii)%item(i)%Serialized%valor_y(Ntimeforvolumic,i1), &
                                                                   output(ii)%item(i)%Serialized%valor_z(Ntimeforvolumic,i1)
                                  end do  
                              endif
                           end select
                        endif
                     endif
                  end do
                  !!!!!!!!!!!!!!!!!!!
                  !
                  call flush(output(ii)%item(i)%unit)
               end select

            elseif (SGG%Observation(ii)%FreqDomain) then !only volumic probes are handled in freq domain in this way
               select case (field)
                case (iMHC,iHxC,iHyC,iHzC,iMEC,iExC,iEyC,iEzC)
                  at=sgg%tiempo(FinalInstant)
                  !!!!assumo que todos son electricos o magneticos en una probe Volumic para calcular el tiempo !logico
                  output(ii)%TimesWritten=output(ii)%NumFreqs  !util para leer el numero exacto de freq points
                  INQUIRE(file= trim(adjustl(output(ii)%item(i)%path)),OPENED = ISyaopen)
                  if (isyaopen) close (output(ii)%item(i)%unit,status='delete')   
                  !            
                  my_iostat=0
9244              if(my_iostat /= 0) write(*,fmt='(a)',advance='no'), '.' !!if(my_iostat /= 0) print '(i5,a1,i4,2x,a)',9244,'.',layoutnumber,trim(adjustl(output(ii)%item(i)%path)) 
                  open (output(ii)%item(i)%unit,file= trim(adjustl(output(ii)%item(i)%path)),form='unformatted',err=9244,iostat=my_iostat)
                  write(output(ii)%item(i)%unit) &
                            output(ii)%item(i)%XItrancos , output(ii)%item(i)%XEtrancos, &
                            output(ii)%item(i)%YItrancos , output(ii)%item(i)%YEtrancos, &
                            output(ii)%item(i)%ZItrancos , output(ii)%item(i)%ZEtrancos
                !!! &      sgg%observation(ii)%P(i)%xI,sgg%observation(ii)%P(i)%xE, &
                !!! &      sgg%observation(ii)%P(i)%YI,sgg%observation(ii)%P(i)%YE, &
                !!! &      sgg%observation(ii)%P(i)%zI,sgg%observation(ii)%P(i)%ZE
                  write(output(ii)%item(i)%unit) at !deteccion errores dft si se resumea a partir de instantes posteriores al ultimo escrito
                  DO N=1,output(ii)%NumFreqs
                     write(output(ii)%item(i)%unit) output(ii)%Freq(n)
                     do compo=1,3
                              do k1t=output(ii)%item(i)%ZItrancos,output(ii)%item(i)%ZEtrancos
                                 do j1t=output(ii)%item(i)%YItrancos,output(ii)%item(i)%YEtrancos
                                  IF (SGG%Observation(ii)%Transfer) then
                                     write(output(ii)%item(i)%unit) (output(ii)%item(i)%valor3DComplex(N,compo,i1t,j1t,k1t)/output(ii)%dftEntrada(n), &
                                     &                                             i1t=output(ii)%item(i)%XItrancos,output(ii)%item(i)%XEtrancos)
                                  else !solo la transformada sin normalizar
                                     write(output(ii)%item(i)%unit) (output(ii)%item(i)%valor3DComplex(N,compo,i1t,j1t,k1t), &
                                     &                                             i1t=output(ii)%item(i)%XItrancos,output(ii)%item(i)%XEtrancos)
                                  endif
                                end do
                        end do
                     end do
                  end do
                  close (output(ii)%item(i)%unit)
                case(icur,iCurX,iCurY,iCurZ)  !!!quitadp de aqui el mapvtk porque nunca puede estar en frecuencia!!!! 050216���
                  at=sgg%tiempo(FinalInstant)
                  output(ii)%TimesWritten=output(ii)%NumFreqs  !util para leer el numero exacto de freq points
                  INQUIRE(file= trim(adjustl(output(ii)%item(i)%path)),OPENED = ISyaopen)
                  if (isyaopen) close (output(ii)%item(i)%unit,status='delete')   
                  !               
                  my_iostat=0
9245              if(my_iostat /= 0) write(*,fmt='(a)',advance='no'), '.' !!if(my_iostat /= 0) print '(i5,a1,i4,2x,a)',9245,'.',layoutnumber,trim(adjustl(output(ii)%item(i)%path)) 
                  open (output(ii)%item(i)%unit,file= trim(adjustl(output(ii)%item(i)%path)),form='unformatted',err=9245,iostat=my_iostat)
                  write(output(ii)%item(i)%unit) output(ii)%item(i)%columnas
                  do conta=1,output(ii)%item(i)%columnas
                     write(output(ii)%item(i)%unit) &
                     output(ii)%item(i)%Serialized%eI(conta), &
                     output(ii)%item(i)%Serialized%eJ(conta), &
                     output(ii)%item(i)%Serialized%eK(conta), &
                     output(ii)%item(i)%Serialized%currentType(conta), &
                     output(ii)%item(i)%Serialized%sggMtag(conta) !added to resuming file 121020
                  end do
                  write(output(ii)%item(i)%unit) at !deteccion errores dft si se resumea a partir de instantes posteriores al ultimo escrito
                  DO N=1,output(ii)%NumFreqs
                     write(output(ii)%item(i)%unit) output(ii)%Freq(n)
                        IF (SGG%Observation(ii)%Transfer) then
                           if (output(ii)%item(i)%columnas /=0) then    
                               do i1=1,output(ii)%item(i)%columnas
                                  write(output(ii)%item(i)%unit) &     
                                          output(ii)%item(i)%Serialized%valorComplex_x(N,i1)/output(ii)%dftEntrada(n), &
                                          output(ii)%item(i)%Serialized%valorComplex_y(N,i1)/output(ii)%dftEntrada(n), &
                                          output(ii)%item(i)%Serialized%valorComplex_z(N,i1)/output(ii)%dftEntrada(n) 
                               end do
                            endif
                        ELSE
                           if (output(ii)%item(i)%columnas /=0) then
                               
                                do i1=1,output(ii)%item(i)%columnas
                                    write(output(ii)%item(i)%unit) & 
                                       output(ii)%item(i)%Serialized%valorComplex_x(N,i1), &
                                       output(ii)%item(i)%Serialized%valorComplex_y(N,i1), &
                                       output(ii)%item(i)%Serialized%valorComplex_z(N,i1) 
                                end do
                            endif

                        ENDIF
                  end do
                  close (output(ii)%item(i)%unit)
               end select
               !
            endif !del time domain
         end do  loop_obb
      end do

      nInit=FinalInstant+1
      !voids valor
      do ii=1,sgg%NumberRequest
         do i=1,sgg%Observation(ii)%nP
            if (SGG%Observation(ii)%TimeDomain) then
               field=sgg%observation(ii)%P(i)%what
               select case(field)
                  !estas sondas no se anulan tras escribir
                  !case (iCur,iCurX,iCurY,iCurZ,mapvtk)
                  !    output(ii)%item(i)%Serialized%valor = 0.0_RKIND
                  !case (iMHC,iHxC,iHyC,iHzC,iMEC,iExC,iEyC,iEzC)
                  !    output(ii)%item(i)%valor3D = 0.0_RKIND
                case (iJx,iJy,iJz)
                  output(ii)%item(i)%valor(0:BuffObse)=0.0_RKIND
                  output(ii)%item(i)%valor2(0:BuffObse)=0.0_RKIND   
                  output(ii)%item(i)%valor3(0:BuffObse)=0.0_RKIND
                  output(ii)%item(i)%valor4(0:BuffObse)=0.0_RKIND
                  output(ii)%item(i)%valor5(0:BuffObse)=0.0_RKIND
                case (iBloqueMx,iBloqueMz,iBloqueMy,iBloqueJx,iBloqueJz,iBloqueJy,iEx,iEy,iEz,iHx,iHy,iHz)
                  output(ii)%item(i)%valor(0:BuffObse)=0.0_RKIND
               end select
            endif
         end do
      end do


      return
   end subroutine FlushObservationFiles

   subroutine FlushMTLNObservationFiles(nEntradaRoot)
      character (len=*), intent(in)  ::  nEntradaRoot
      type(mtln_solver_t), pointer :: mtln_solver
      integer :: i,j,k,n
      integer :: unit 
      character (len=bufsize)  ::  path, temp
      character (len=:), allocatable :: buffer

      mtln_solver => GetSolverPtr()
      unit = 2000
      do i = 1, size(mtln_solver%bundles)
         do j = 1, size(mtln_solver%bundles(i)%probes)
            path = trim(nEntradaRoot)//"_"//trim(mtln_solver%bundles(i)%probes(j)%name)//".dat"
            open(unit = unit , file = trim(path))

            buffer = "time"

            do k = 1, size(mtln_solver%bundles(i)%probes(j)%val,2)
               write(temp,*) k
               buffer = buffer//" "//"conductor_"//trim(adjustl(temp))
            end do
            write(unit,*) trim(buffer)
            do k = 1, size(mtln_solver%bundles(i)%probes(j)%t)
               buffer = ""
               write(temp,*) mtln_solver%bundles(i)%probes(j)%t(k)
               buffer = buffer//trim(temp)
               do n = 1, size(mtln_solver%bundles(i)%probes(j)%val,2)
                  write(temp,*) mtln_solver%bundles(i)%probes(j)%val(k,n)
                  buffer = buffer//" "//trim(temp)
               end do
               write(unit,*) trim(buffer)
            end do
            close(unit)
            unit = unit + 1
         end do
      end do
      
   end subroutine





   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Free up memory
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine DestroyObservation(sgg)
      type (SGGFDTDINFO), intent(INOUT)         ::  sgg
      integer (kind=4)  ::  ii,i,field

#ifdef CompileWithMPI
      integer (kind=4) :: ierr
#endif

#ifdef CompileWithMPI
      if (associated(valores)) deallocate (valores,newvalores)
#endif
      do ii=1,sgg%NumberRequest
         if (SGG%Observation(ii)%Transfer) deallocate (output(ii)%dftEntrada)
         if (SGG%Observation(ii)%FreqDomain) deallocate (output(ii)%auxExp_E,output(ii)%auxExp_H,output(ii)%Freq)
         DO i=1,sgg%Observation(ii)%nP
            field=sgg%observation(ii)%P(i)%what
            select case(field)
             case (iJx,iJy,iJz)
#ifdef CompileWithWires
               deallocate (output(ii)%item(i)%valor)
               deallocate (output(ii)%item(i)%valor2,output(ii)%item(i)%valor3,output(ii)%item(i)%valor4,output(ii)%item(i)%valor5)  !en caso de hilos se necesitan
#endif
             case (iBloqueJx,iBloqueJy,iBloqueMx,iBloqueMy)
               deallocate (output(ii)%item(i)%valor)
#ifdef CompileWithMPI
               if (output(ii)%item(i)%MPISubComm /= -1) then
                  call MPI_Group_free(output(ii)%item(i)%MPIgroupindex,ierr)
               endif
#endif
             case (iMHC,iHxC,iHyC,iHzC,iMEC,iExC,iEyC,iEzC)
#ifdef CompileWithMPI
               if (output(ii)%item(i)%MPISubComm /= -1) then
                  call MPI_Group_free(output(ii)%item(i)%MPIgroupindex,ierr)
               endif
#endif
               if (SGG%Observation(ii)%TimeDomain) deallocate (output(ii)%item(i)%valor3D)
               if (SGG%Observation(ii)%FreqDomain) deallocate (output(ii)%item(i)%valor3DComplex)
             case (iCur,iCurX,iCurY,iCurZ,mapvtk) !!!���
#ifdef CompileWithMPI
               if (output(ii)%item(i)%MPISubComm /= -1) then
                  call MPI_Group_free(output(ii)%item(i)%MPIgroupindex,ierr)
               endif
#endif
               if (SGG%Observation(ii)%TimeDomain) deallocate (output(ii)%item(i)%Serialized%valor,   &
                                                               output(ii)%item(i)%Serialized%valor_x, &
                                                               output(ii)%item(i)%Serialized%valor_y, &
                                                               output(ii)%item(i)%Serialized%valor_z)
               if (SGG%Observation(ii)%FreqDomain) then         
                   deallocate (output(ii)%item(i)%Serialized%valorComplex_x)
                   deallocate (output(ii)%item(i)%Serialized%valorComplex_y)
                   deallocate (output(ii)%item(i)%Serialized%valorComplex_z)
               endif
               deallocate (output(ii)%item(i)%Serialized%eI)
               deallocate (output(ii)%item(i)%Serialized%eJ)
               deallocate (output(ii)%item(i)%Serialized%eK)
               deallocate (output(ii)%item(i)%Serialized%currentType)
               deallocate (output(ii)%item(i)%Serialized%sggMtag)

             case (iBloqueMz,iBloqueJz,iEx,iEy,iEz,iHx,iHy,iHz)
               deallocate (output(ii)%item(i)%valor)
#ifdef CompileWithNF2FF
             case (farfield)
               call DestroyFarField
#ifdef CompileWithMPI
               if (output(ii)%item(i)%MPISubComm /= -1) then
                  call MPI_Group_free(output(ii)%item(i)%MPIgroupindex,ierr)
               endif
#endif
#endif
            end select
         end do
         if (associated(sgg%Observation(ii)%P)) deallocate (sgg%Observation(ii)%P)
         if (associated(output(ii)%item)) deallocate (output(ii)%item)

      end do

      if (associated(sgg%Observation)) deallocate (sgg%Observation)
      if (associated(output)) deallocate (output)

   end subroutine


   !!!!!!!!!!!!!!!!!!!!!

   function prefix(campo) result(ext)
      integer (kind=4)  ::  campo
      character (len=BUFSIZE)  ::  ext

      select case (campo)
       case (iEx)
         ext='Ex_'
       case (iEy)
         ext='Ey_'
       case (iEz)
         ext='Ez_'
       case (iVx)
         ext='Vx_'
       case (iVy)
         ext='Vy_'
       case (iVz)
         ext='Vz_'
       case (iHx)
         ext='Hx_'
       case (iHy)
         ext='Hy_'
       case (iHz)
         ext='Hz_'
       case (iBloqueJx)
         ext='Jx_'
       case (iBloqueJy)
         ext='Jy_'
       case (iBloqueJz)
         ext='Jz_'
       case (iBloqueMx)
         ext='Mx_'
       case (iBloqueMy)
         ext='My_'
       case (iBloqueMz)
         ext='Mz_'
       case (iJx)
         ext='Wx_'
       case (iJy)
         ext='Wy_'
       case (iJz)
         ext='Wz_'
       case (iExC)
         ext='ExC_'
       case (iEyC)
         ext='EyC_'
       case (iEzC)
         ext='EzC_'
       case (iHxC)
         ext='HxC_'
       case (iHyC)
         ext='HyC_'
       case (iHzC)
         ext='HzC_'
       case (iMEC)
         ext='ME_'
       case (iMHC)
         ext='MH_'
       case (iCur)
         ext='BC_'
       case (mapvtk)
         ext='MAP_'
       case (iCurX)
         ext='BCX_'
       case (iCurY)
         ext='BCY_'
       case (iCurZ)
         ext='BCZ_'
#ifdef CompileWithNF2FF
       case (farfield)
         ext='FF_'
#endif
      end select

      return

   end function prefix

   function suffix(campo,incid) result(ext)
      integer (kind=4)  ::  campo
      character (LEN=BUFSIZE)  ::  ext
      logical  ::  incid

      ext=' '

      select case (campo)
       case (iEx,iEy,iEz,iHx,iHy,iHz)
         if (incid) ext='incid'
       case (iJx,iJy,iJz)
         ext=' -E*dl Vplus Vminus Vplus-Vminus' 
      end select


      return

   end function suffix


   function fieldo(field,dir) result(fieldo2)
      integer  ::  fieldo2,field
      character (len=1) :: dir
      fieldo2=-1
      select case(field)
       case(iEx,iEy,iEz,iHx,iHy,iHz)
         fieldo2=field
       case (iJx,iVx,iBloqueJx,iExC)
         fieldo2=iEx
       case (iJy,iVy,iBloqueJy,iEyC)
         fieldo2=iEy
       case (iJz,iVz,iBloqueJz,iEzC)
         fieldo2=iEz
       case (iBloqueMx,iHxC)
         fieldo2=iHx
       case (iBloqueMy,iHyC)
         fieldo2=iHy
       case (iBloqueMz,iHzC)
         fieldo2=iHz
       case(iMEC)
         select case (dir)
          CASE ('X','x')
            fieldo2=iEx
          CASE ('Y','y')
            fieldo2=iEY
          CASE ('Z','z')
            fieldo2=iEz
         END SELECT
       case(iMHC)
         select case (dir)
          CASE ('X','x')
            fieldo2=ihx
          CASE ('Y','y')
            fieldo2=iHY
          CASE ('Z','z')
            fieldo2=iHz
         END SELECT
       case(iCur,iCurX,icurY,icurZ,mapvtk)  !los pongo en efield para evitar problemas con el MPI
         select case (dir)
          CASE ('X','x')
            fieldo2=iEx
          CASE ('Y','y')
            fieldo2=iEY
          CASE ('Z','z')
            fieldo2=iEz
         END SELECT
      end select
   end function

   !!!cuenta los bordes adyacentes
   subroutine contabordes(sgg,imed,imed1,imed2,imed3,imed4,EsBorde,SINPML_fullsize,campo,iii,jjj,kkk)
      type (SGGFDTDINFO), intent(IN)       ::  sgg
      type (limit_t), dimension(1:6), intent(in)  ::  SINPML_fullsize
      integer(Kind=4) imed,imed1,imed2,imed3,imed4,contaborde,campo,iii,jjj,kkk
      logical :: esborde
      !!!!
      esborde=.false.
      contaborde=0
      !esta primera opcion solo considera bordes los externos
      IF (imed/=1) THEN
         !    if     (sgg%med(imed )%is%SGBC) then
         !        if (sgg%med(imed1)%is%SGBC) THEN
         !            if (trim(adjustl(sgg%Med(imed )%Multiport(1)%MultiportFileZ11))==trim(adjustl(sgg%Med(imed1)%Multiport(1)%MultiportFileZ11)) ) contaborde=contaborde+1
         !        endif
         !        if (sgg%med(imed2)%is%SGBC) THEN
         !            if (trim(adjustl(sgg%Med(imed )%Multiport(1)%MultiportFileZ11))==trim(adjustl(sgg%Med(imed2)%Multiport(1)%MultiportFileZ11)) ) contaborde=contaborde+1
         !        endif
         !        if (sgg%med(imed3)%is%SGBC) THEN
         !            if (trim(adjustl(sgg%Med(imed )%Multiport(1)%MultiportFileZ11))==trim(adjustl(sgg%Med(imed3)%Multiport(1)%MultiportFileZ11)) ) contaborde=contaborde+1
         !        endif
         !        if (sgg%med(imed4)%is%SGBC) THEN
         !            if (trim(adjustl(sgg%Med(imed )%Multiport(1)%MultiportFileZ11))==trim(adjustl(sgg%Med(imed4)%Multiport(1)%MultiportFileZ11)) ) contaborde=contaborde+1
         !        endif
         !   elseif  (sgg%med(imed )%is%Multiport) then
         !        if (sgg%med(imed1)%is%Multiport) THEN
         !            if (trim(adjustl(sgg%Med(imed )%Multiport(1)%MultiportFileZ11))==trim(adjustl(sgg%Med(imed1)%Multiport(1)%MultiportFileZ11)) ) contaborde=contaborde+1
         !        endif
         !        if (sgg%med(imed2)%is%Multiport) THEN
         !            if (trim(adjustl(sgg%Med(imed )%Multiport(1)%MultiportFileZ11))==trim(adjustl(sgg%Med(imed2)%Multiport(1)%MultiportFileZ11)) ) contaborde=contaborde+1
         !        endif
         !        if (sgg%med(imed3)%is%Multiport) THEN
         !            if (trim(adjustl(sgg%Med(imed )%Multiport(1)%MultiportFileZ11))==trim(adjustl(sgg%Med(imed3)%Multiport(1)%MultiportFileZ11)) ) contaborde=contaborde+1
         !        endif
         !        if (sgg%med(imed4)%is%Multiport) THEN
         !            if (trim(adjustl(sgg%Med(imed )%Multiport(1)%MultiportFileZ11))==trim(adjustl(sgg%Med(imed4)%Multiport(1)%MultiportFileZ11)) ) contaborde=contaborde+1
         !        endif
         !    elseif (sgg%med(imed )%is%AnisMultiport) then
         !        if (sgg%med(imed1)%is%AnisMultiport) THEN
         !            if (trim(adjustl(sgg%Med(imed )%AnisMultiport(1)%MultiportFileZ11))==trim(adjustl(sgg%Med(imed1)%AnisMultiport(1)%MultiportFileZ11)) ) contaborde=contaborde+1
         !        endif
         !        if (sgg%med(imed2)%is%AnisMultiport) THEN
         !            if (trim(adjustl(sgg%Med(imed )%AnisMultiport(1)%MultiportFileZ11))==trim(adjustl(sgg%Med(imed2)%AnisMultiport(1)%MultiportFileZ11)) ) contaborde=contaborde+1
         !        endif
         !        if (sgg%med(imed3)%is%AnisMultiport) THEN
         !            if (trim(adjustl(sgg%Med(imed )%AnisMultiport(1)%MultiportFileZ11))==trim(adjustl(sgg%Med(imed3)%AnisMultiport(1)%MultiportFileZ11)) ) contaborde=contaborde+1
         !        endif
         !        if (sgg%med(imed4)%is%AnisMultiport) THEN
         !            if (trim(adjustl(sgg%Med(imed )%AnisMultiport(1)%MultiportFileZ11))==trim(adjustl(sgg%Med(imed4)%AnisMultiport(1)%MultiportFileZ11)) ) contaborde=contaborde+1
         !        endif
         !    else
         !        if (imed==imed1) contaborde=contaborde+1
         !        if (imed==imed2) contaborde=contaborde+1
         !        if (imed==imed3) contaborde=contaborde+1
         !        if (imed==imed4) contaborde=contaborde+1
         !    endif
         !    if (contaborde <=1) esborde=.true.
         !!!!alternativa
         if     (sgg%med(imed )%is%SGBC) then
            if (sgg%med(imed1)%is%SGBC) THEN
               if (trim(adjustl(sgg%Med(imed )%Multiport(1)%MultiportFileZ11))/=trim(adjustl(sgg%Med(imed1)%Multiport(1)%MultiportFileZ11)) ) contaborde=contaborde+1
            elseif (imed1/=1)  THEN
               contaborde=contaborde+1
            endif
            if (sgg%med(imed2)%is%SGBC) THEN
               if (trim(adjustl(sgg%Med(imed )%Multiport(1)%MultiportFileZ11))/=trim(adjustl(sgg%Med(imed2)%Multiport(1)%MultiportFileZ11)) ) contaborde=contaborde+1
            elseif (imed2/=1)  THEN
               contaborde=contaborde+1
            endif
            if (sgg%med(imed3)%is%SGBC) THEN
               if (trim(adjustl(sgg%Med(imed )%Multiport(1)%MultiportFileZ11))/=trim(adjustl(sgg%Med(imed3)%Multiport(1)%MultiportFileZ11)) ) contaborde=contaborde+1
            elseif (imed3/=1)  THEN
               contaborde=contaborde+1
            endif
            if (sgg%med(imed4)%is%SGBC) THEN
               if (trim(adjustl(sgg%Med(imed )%Multiport(1)%MultiportFileZ11))/=trim(adjustl(sgg%Med(imed4)%Multiport(1)%MultiportFileZ11)) ) contaborde=contaborde+1
            elseif (imed4/=1)  THEN
               contaborde=contaborde+1
            endif
         elseif  (sgg%med(imed )%is%Multiport) then
            if (sgg%med(imed1)%is%Multiport) THEN
               if (trim(adjustl(sgg%Med(imed )%Multiport(1)%MultiportFileZ11))/=trim(adjustl(sgg%Med(imed1)%Multiport(1)%MultiportFileZ11)) ) contaborde=contaborde+1
            elseif (imed1/=1)  THEN
               contaborde=contaborde+1
            endif
            if (sgg%med(imed2)%is%Multiport) THEN
               if (trim(adjustl(sgg%Med(imed )%Multiport(1)%MultiportFileZ11))/=trim(adjustl(sgg%Med(imed2)%Multiport(1)%MultiportFileZ11)) ) contaborde=contaborde+1
            elseif (imed2/=1)  THEN
               contaborde=contaborde+1
            endif
            if (sgg%med(imed3)%is%Multiport) THEN
               if (trim(adjustl(sgg%Med(imed )%Multiport(1)%MultiportFileZ11))/=trim(adjustl(sgg%Med(imed3)%Multiport(1)%MultiportFileZ11)) ) contaborde=contaborde+1
            elseif (imed3/=1)  THEN
               contaborde=contaborde+1
            endif
            if (sgg%med(imed4)%is%Multiport) THEN
               if (trim(adjustl(sgg%Med(imed )%Multiport(1)%MultiportFileZ11))/=trim(adjustl(sgg%Med(imed4)%Multiport(1)%MultiportFileZ11)) ) contaborde=contaborde+1
            elseif (imed4/=1)  THEN
               contaborde=contaborde+1
            endif
         elseif (sgg%med(imed )%is%AnisMultiport) then
            if (sgg%med(imed1)%is%AnisMultiport) THEN
               if (trim(adjustl(sgg%Med(imed )%AnisMultiport(1)%MultiportFileZ11))/=trim(adjustl(sgg%Med(imed1)%AnisMultiport(1)%MultiportFileZ11)) ) contaborde=contaborde+1
            elseif (imed1/=1)  THEN
               contaborde=contaborde+1
            endif
            if (sgg%med(imed2)%is%AnisMultiport) THEN
               if (trim(adjustl(sgg%Med(imed )%AnisMultiport(1)%MultiportFileZ11))/=trim(adjustl(sgg%Med(imed2)%AnisMultiport(1)%MultiportFileZ11)) ) contaborde=contaborde+1
            elseif (imed2/=1)  THEN
               contaborde=contaborde+1
            endif
            if (sgg%med(imed3)%is%AnisMultiport) THEN
               if (trim(adjustl(sgg%Med(imed )%AnisMultiport(1)%MultiportFileZ11))/=trim(adjustl(sgg%Med(imed3)%AnisMultiport(1)%MultiportFileZ11)) ) contaborde=contaborde+1
            elseif (imed3/=1)  THEN
               contaborde=contaborde+1
            endif
            if (sgg%med(imed4)%is%AnisMultiport) THEN
               if (trim(adjustl(sgg%Med(imed )%AnisMultiport(1)%MultiportFileZ11))/=trim(adjustl(sgg%Med(imed4)%AnisMultiport(1)%MultiportFileZ11)) ) contaborde=contaborde+1
            elseif (imed4/=1)  THEN
               contaborde=contaborde+1
            endif
         else
            if ((imed/=imed1).and.(imed1/=1)) contaborde=contaborde+1
            if ((imed/=imed2).and.(imed2/=1)) contaborde=contaborde+1
            if ((imed/=imed3).and.(imed3/=1)) contaborde=contaborde+1
            if ((imed/=imed4).and.(imed4/=1)) contaborde=contaborde+1
         endif
         if ((imed1==1).and.(imed2==1).and.(imed3==1).and.(imed4/=1)) esborde=.true. !un borde con vacion
         if ((imed2==1).and.(imed3==1).and.(imed4==1).and.(imed1/=1)) esborde=.true. !un borde con vacion
         if ((imed3==1).and.(imed4==1).and.(imed1==1).and.(imed2/=1)) esborde=.true. !un borde con vacion
         if ((imed4==1).and.(imed1==1).and.(imed2==1).and.(imed3/=1)) esborde=.true. !un borde con vacion
         if (contaborde > 0) esborde=.true.
         if ((imed1==1).and.(imed2==1).and.(imed3==1).and.(imed4==1)) esborde=.true. !un segmento aislado
         !!!!!! Fin de la alternativa para que considere bordes tambi�n ejes internos
         !!!!!!!!!!!!!!
         !!!!!!!!!!!!!
         if ((iii > SINPML_fullsize(campo)%XE).or.(jjj > SINPML_fullsize(campo)%YE).or.(kkk > SINPML_fullsize(campo)%ZE)) then
            esborde=.false.
         endif
         if ((iii < SINPML_fullsize(campo)%XI).or.(jjj < SINPML_fullsize(campo)%Yi).or.(kkk < SINPML_fullsize(campo)%Zi)) then
            esborde=.false.
         endif
      else
         esborde=.false.
      ENDIF !DEL IMED1
      return
    end subroutine contabordes

    


#ifdef CompileWithNodalSources
   subroutine nodalvtk(sgg,sggMiEx,sggMiEy,sggMiEz,sggMiHx,sggMiHy,sggMiHz,sggMtag, &
                         b,init,geom,asigna,electric,magnetic,conta,i,ii,output,Ntimeforvolumic)
      type (SGGFDTDINFO), intent(IN)         ::  sgg
      INTEGER (KIND=IKINDMTAG), intent(in) :: sggMtag  (sgg%Alloc(iHx)%XI:sgg%Alloc(iHx)%XE, sgg%Alloc(iHy)%YI:sgg%Alloc(iHy)%YE, sgg%Alloc(iHz)%ZI:sgg%Alloc(iHz)%ZE)

      type (output_t), pointer, dimension( : )  ::  output
      integer(kind=4), intent(IN) :: i,ii,Ntimeforvolumic
      type( bounds_t), intent( IN)  ::  b
      logical geom,asigNa,init,electric,magnetic
      integer(kind=4) conta,sweep,ni,nj,nk,i_m,j_m,k_m,IMED

      integer(kind = INTEGERSIZEOFMEDIAMATRICES), dimension ( 0 : b%sggMiEx%NX-1 , 0 : b%sggMiEx%NY-1 , 0 : b%sggMiEx%NZ-1 )  , intent( IN)     ::  sggMiEx
      integer(kind = INTEGERSIZEOFMEDIAMATRICES), dimension ( 0 : b%sggMiEy%NX-1 , 0 : b%sggMiEy%NY-1 , 0 : b%sggMiEy%NZ-1 )  , intent( IN)     ::  sggMiEy
      integer(kind = INTEGERSIZEOFMEDIAMATRICES), dimension ( 0 : b%sggMiEz%NX-1 , 0 : b%sggMiEz%NY-1 , 0 : b%sggMiEz%NZ-1 )  , intent( IN)     ::  sggMiEz
      integer(kind = INTEGERSIZEOFMEDIAMATRICES), dimension ( 0 : b%sggMiHx%NX-1 , 0 : b%sggMiHx%NY-1 , 0 : b%sggMiHx%NZ-1 )  , intent( IN)     ::  sggMiHx
      integer(kind = INTEGERSIZEOFMEDIAMATRICES), dimension ( 0 : b%sggMiHy%NX-1 , 0 : b%sggMiHy%NY-1 , 0 : b%sggMiHy%NZ-1 )  , intent( IN)     ::  sggMiHy
      integer(kind = INTEGERSIZEOFMEDIAMATRICES), dimension ( 0 : b%sggMiHz%NX-1 , 0 : b%sggMiHz%NY-1 , 0 : b%sggMiHz%NZ-1 )  , intent( IN)     ::  sggMiHz


      !to fetch info of nodal sources for the vtkmap
      type (nodsou), pointer, save :: rNodal_Ex,rNodal_Ey,rNodal_Ez
      type (nodsou), pointer, save :: rNodal_Hx,rNodal_Hy,rNodal_Hz
      !!!!!!!!!!!!


      if (init) call getnodal(rNodal_Ex,rNodal_Ey,rNodal_Ez,rNodal_Hx,rNodal_Hy,rNodal_Hz)
      IF (ELECTRIC) THEN
         do sweep=1,rNodal_Ex%numHard
            do nk=rNodal_Ex%nodHard(sweep)%punto%zi,rNodal_Ex%nodHard(sweep)%punto%ze
               k_m = nk - b%Ex%ZI
               do nj=rNodal_Ex%nodHard(sweep)%punto%yi,rNodal_Ex%nodHard(sweep)%punto%ye
                  j_m = nj - b%Ex%YI
                  do ni=rNodal_Ex%nodHard(sweep)%punto%xi,rNodal_Ex%nodHard(sweep)%punto%xe
                     i_m = ni - b%Ex%XI
                     imed = sggMiEx(i_m,j_m,k_m)
                     if (.not.sgg%Med(imed)%Is%PEC) then
                        conta=conta+1
                        if (geom) then
                           ! print *,'-------antes de GEOM ',II,I, conta
                           output(ii)%item(i)%Serialized%eI(conta)=ni
                           output(ii)%item(i)%Serialized%eJ(conta)=nj
                           output(ii)%item(i)%Serialized%eK(conta)=nk
                           output(ii)%item(i)%Serialized%currentType(conta)=iJx
                           output(ii)%item(i)%Serialized%sggMtag(conta)=iabs(sggMtag(ni,nj,nk))
                           ! print *,'-------tras  GEOM',output(ii)%item(i)%Serialized%eI(conta),output(ii)%item(i)%Serialized%currentType(conta)
                        endif
                        if (asigna) then
                           ! print *,'-------antes de asigna ', Ntimeforvolumic,conta
                           output( ii)%item( i)%Serialized%valor(Ntimeforvolumic,conta) = 8.5
                           ! print *,'-------tras  de asigna ',output( ii)%item( i)%Serialized%valor(Ntimeforvolumic,conta)
                        endif
                     endif
                  end do
               End do
            end do
         end do
         !
         do sweep=1,rNodal_Ex%numSoft
            do nk=rNodal_Ex%nodSoft(sweep)%punto%zi,rNodal_Ex%nodSoft(sweep)%punto%ze
               k_m = nk - b%Ex%ZI
               do nj=rNodal_Ex%nodSoft(sweep)%punto%yi,rNodal_Ex%nodSoft(sweep)%punto%ye
                  j_m = nj - b%Ex%YI
                  do ni=rNodal_Ex%nodSoft(sweep)%punto%xi,rNodal_Ex%nodSoft(sweep)%punto%xe
                     i_m = ni - b%Ex%XI
                     imed = sggMiEx(i_m,j_m,k_m)
                     if (.not.sgg%Med(imed)%Is%PEC) then
                        conta=conta+1
                        if (geom) then
                           ! print *,'-------antes de GEOM ',II,I, conta
                           output(ii)%item(i)%Serialized%eI(conta)=ni
                           output(ii)%item(i)%Serialized%eJ(conta)=nj
                           output(ii)%item(i)%Serialized%eK(conta)=nk
                           output(ii)%item(i)%Serialized%currentType(conta)=iJx
                           output(ii)%item(i)%Serialized%sggMtag(conta)=iabs(sggMtag(ni,nj,nk))
                           ! print *,'-------tras  GEOM',output(ii)%item(i)%Serialized%eI(conta),output(ii)%item(i)%Serialized%currentType(conta)
                        endif
                        if (asigna) then
                           ! print *,'-------antes de asigna ', Ntimeforvolumic,conta
                           output( ii)%item( i)%Serialized%valor(Ntimeforvolumic,conta) = 8.5
                           ! print *,'-------tras  de asigna ',output( ii)%item( i)%Serialized%valor(Ntimeforvolumic,conta)
                        endif
                     endif
                  end do
               End do
            end do
         end do
         !
         !
         do sweep=1,rNodal_Ey%numHard
            do nk=rNodal_Ey%nodHard(sweep)%punto%zi,rNodal_Ey%nodHard(sweep)%punto%ze
               k_m = nk - b%Ey%ZI
               do nj=rNodal_Ey%nodHard(sweep)%punto%yi,rNodal_Ey%nodHard(sweep)%punto%ye
                  j_m = nj - b%Ey%YI
                  do ni=rNodal_Ey%nodHard(sweep)%punto%xi,rNodal_Ey%nodHard(sweep)%punto%xe
                     i_m = ni - b%Ey%XI
                     imed = sggMiEy(i_m,j_m,k_m)
                     if (.not.sgg%Med(imed)%Is%PEC) then
                        conta=conta+1
                        if (geom) then
                           ! print *,'-------antes de GEOM ',II,I, conta
                           output(ii)%item(i)%Serialized%eI(conta)=ni
                           output(ii)%item(i)%Serialized%eJ(conta)=nj
                           output(ii)%item(i)%Serialized%eK(conta)=nk
                           output(ii)%item(i)%Serialized%currentType(conta)=iJy
                           output(ii)%item(i)%Serialized%sggMtag(conta)=iabs(sggMtag(ni,nj,nk))
                           ! print *,'-------tras  GEOM',output(ii)%item(i)%Serialized%eI(conta),output(ii)%item(i)%Serialized%currentType(conta)
                        endif
                        if (asigna) then
                           ! print *,'-------antes de asigna ', Ntimeforvolumic,conta
                           output( ii)%item( i)%Serialized%valor(Ntimeforvolumic,conta) = 8.5
                           ! print *,'-------tras  de asigna ',output( ii)%item( i)%Serialized%valor(Ntimeforvolumic,conta)
                        endif
                     endif
                  end do
               End do
            end do
         end do
         !
         do sweep=1,rNodal_Ey%numSoft
            do nk=rNodal_Ey%nodSoft(sweep)%punto%zi,rNodal_Ey%nodSoft(sweep)%punto%ze
               k_m = nk - b%Ey%ZI
               do nj=rNodal_Ey%nodSoft(sweep)%punto%yi,rNodal_Ey%nodSoft(sweep)%punto%ye
                  j_m = nj - b%Ey%YI
                  do ni=rNodal_Ey%nodSoft(sweep)%punto%xi,rNodal_Ey%nodSoft(sweep)%punto%xe
                     i_m = ni - b%Ey%XI
                     imed = sggMiEy(i_m,j_m,k_m)
                     if (.not.sgg%Med(imed)%Is%PEC) then
                        conta=conta+1
                        if (geom) then
                           ! print *,'-------antes de GEOM ',II,I, conta
                           output(ii)%item(i)%Serialized%eI(conta)=ni
                           output(ii)%item(i)%Serialized%eJ(conta)=nj
                           output(ii)%item(i)%Serialized%eK(conta)=nk
                           output(ii)%item(i)%Serialized%currentType(conta)=iJy
                           output(ii)%item(i)%Serialized%sggMtag(conta)=iabs(sggMtag(ni,nj,nk))
                           ! print *,'-------tras  GEOM',output(ii)%item(i)%Serialized%eI(conta),output(ii)%item(i)%Serialized%currentType(conta)
                        endif
                        if (asigna) then
                           ! print *,'-------antes de asigna ', Ntimeforvolumic,conta
                           output( ii)%item( i)%Serialized%valor(Ntimeforvolumic,conta) = 8.5
                           ! print *,'-------tras  de asigna ',output( ii)%item( i)%Serialized%valor(Ntimeforvolumic,conta)
                        endif
                     endif
                  end do
               End do
            end do
         end do

         do sweep=1,rNodal_Ez%numHard
            do nk=rNodal_Ez%nodHard(sweep)%punto%zi,rNodal_Ez%nodHard(sweep)%punto%ze
               k_m = nk - b%Ez%ZI
               do nj=rNodal_Ez%nodHard(sweep)%punto%yi,rNodal_Ez%nodHard(sweep)%punto%ye
                  j_m = nj - b%Ez%YI
                  do ni=rNodal_Ez%nodHard(sweep)%punto%xi,rNodal_Ez%nodHard(sweep)%punto%xe
                     i_m = ni - b%Ez%XI
                     imed = sggMiEz(i_m,j_m,k_m)
                     if (.not.sgg%Med(imed)%Is%PEC) then
                        conta=conta+1
                        if (geom) then
                           ! print *,'-------antes de GEOM ',II,I, conta
                           output(ii)%item(i)%Serialized%eI(conta)=ni
                           output(ii)%item(i)%Serialized%eJ(conta)=nj
                           output(ii)%item(i)%Serialized%eK(conta)=nk
                           output(ii)%item(i)%Serialized%currentType(conta)=iJz
                           output(ii)%item(i)%Serialized%sggMtag(conta)=iabs(sggMtag(ni,nj,nk))
                           ! print *,'-------tras  GEOM',output(ii)%item(i)%Serialized%eI(conta),output(ii)%item(i)%Serialized%currentType(conta)
                        endif
                        if (asigna) then
                           ! print *,'-------antes de asigna ', Ntimeforvolumic,conta
                           output( ii)%item( i)%Serialized%valor(Ntimeforvolumic,conta) = 8.5
                           ! print *,'-------tras  de asigna ',output( ii)%item( i)%Serialized%valor(Ntimeforvolumic,conta)
                        endif
                     endif
                  end do
               End do
            end do
         end do
         !
         do sweep=1,rNodal_Ez%numSoft
            do nk=rNodal_Ez%nodSoft(sweep)%punto%zi,rNodal_Ez%nodSoft(sweep)%punto%ze
               k_m = nk - b%Ez%ZI
               do nj=rNodal_Ez%nodSoft(sweep)%punto%yi,rNodal_Ez%nodSoft(sweep)%punto%ye
                  j_m = nj - b%Ez%YI
                  do ni=rNodal_Ez%nodSoft(sweep)%punto%xi,rNodal_Ez%nodSoft(sweep)%punto%xe
                     i_m = ni - b%Ez%XI
                     imed = sggMiEz(i_m,j_m,k_m)
                     if (.not.sgg%Med(imed)%Is%PEC) then
                        conta=conta+1
                        if (geom) then
                           ! print *,'-------antes de GEOM ',II,I, conta
                           output(ii)%item(i)%Serialized%eI(conta)=ni
                           output(ii)%item(i)%Serialized%eJ(conta)=nj
                           output(ii)%item(i)%Serialized%eK(conta)=nk
                           output(ii)%item(i)%Serialized%currentType(conta)=iJz
                           output(ii)%item(i)%Serialized%sggMtag(conta)=iabs(sggMtag(ni,nj,nk))
                           ! print *,'-------tras  GEOM',output(ii)%item(i)%Serialized%eI(conta),output(ii)%item(i)%Serialized%currentType(conta)
                        endif
                        if (asigna) then
                           ! print *,'-------antes de asigna ', Ntimeforvolumic,conta
                           output( ii)%item( i)%Serialized%valor(Ntimeforvolumic,conta) = 8.5
                           ! print *,'-------tras  de asigna ',output( ii)%item( i)%Serialized%valor(Ntimeforvolumic,conta)
                        endif
                     endif
                  end do
               End do
            end do
         end do
      ENDIF !DEL ELECTRIC

      IF (MAGNETIC) THEN

         do sweep=1,rNodal_Hx%numHard
            do nk=rNodal_Hx%nodHard(sweep)%punto%zi,rNodal_Hx%nodHard(sweep)%punto%ze
               k_m = nk - b%Hx%ZI
               do nj=rNodal_Hx%nodHard(sweep)%punto%yi,rNodal_Hx%nodHard(sweep)%punto%ye
                  j_m = nj - b%Hx%YI
                  do ni=rNodal_Hx%nodHard(sweep)%punto%xi,rNodal_Hx%nodHard(sweep)%punto%xe
                     i_m = ni - b%Hx%XI
                     imed = sggMiHx(i_m,j_m,k_m)
                     if (.not.sgg%Med(imed)%Is%PMC) then
                        conta=conta+1
                        if (geom) then
                           output(ii)%item(i)%Serialized%eI(conta)=ni
                           output(ii)%item(i)%Serialized%eJ(conta)=nj
                           output(ii)%item(i)%Serialized%eK(conta)=nk
                           output(ii)%item(i)%Serialized%currentType(conta)=iBloqueJx
                           output(ii)%item(i)%Serialized%sggMtag(conta)=iabs(sggMtag(ni,nj,nk))
                        endif
                        if (asigna) output( ii)%item( i)%Serialized%valor(Ntimeforvolumic,conta) = 9.0
                     endif
                  end do
               End do
            end do
         end do
         !
         do sweep=1,rNodal_Hx%numSoft
            do nk=rNodal_Hx%nodSoft(sweep)%punto%zi,rNodal_Hx%nodSoft(sweep)%punto%ze
               k_m = nk - b%Hx%ZI
               do nj=rNodal_Hx%nodSoft(sweep)%punto%yi,rNodal_Hx%nodSoft(sweep)%punto%ye
                  j_m = nj - b%Hx%YI
                  do ni=rNodal_Hx%nodSoft(sweep)%punto%xi,rNodal_Hx%nodSoft(sweep)%punto%xe
                     i_m = ni - b%Hx%XI
                     imed = sggMiHx(i_m,j_m,k_m)
                     if (.not.sgg%Med(imed)%Is%PMC) then
                        conta=conta+1
                        if (geom) then
                           output(ii)%item(i)%Serialized%eI(conta)=ni
                           output(ii)%item(i)%Serialized%eJ(conta)=nj
                           output(ii)%item(i)%Serialized%eK(conta)=nk
                           output(ii)%item(i)%Serialized%currentType(conta)=iBloqueJx
                           output(ii)%item(i)%Serialized%sggMtag(conta)=iabs(sggMtag(ni,nj,nk))
                        endif
                        if (asigna) output( ii)%item( i)%Serialized%valor(Ntimeforvolumic,conta) = 9.0
                     endif
                  end do
               End do
            end do
         end do
         !
         !
         do sweep=1,rNodal_Hy%numHard
            do nk=rNodal_Hy%nodHard(sweep)%punto%zi,rNodal_Hy%nodHard(sweep)%punto%ze
               k_m = nk - b%Hy%ZI
               do nj=rNodal_Hy%nodHard(sweep)%punto%yi,rNodal_Hy%nodHard(sweep)%punto%ye
                  j_m = nj - b%Hy%YI
                  do ni=rNodal_Hy%nodHard(sweep)%punto%xi,rNodal_Hy%nodHard(sweep)%punto%xe
                     i_m = ni - b%Hy%XI
                     imed = sggMiHx(i_m,j_m,k_m)
                     if (.not.sgg%Med(imed)%Is%PMC) then
                        conta=conta+1
                        if (geom) then
                           output(ii)%item(i)%Serialized%eI(conta)=ni
                           output(ii)%item(i)%Serialized%eJ(conta)=nj
                           output(ii)%item(i)%Serialized%eK(conta)=nk
                           output(ii)%item(i)%Serialized%currentType(conta)=iBloqueJy
                           output(ii)%item(i)%Serialized%sggMtag(conta)=iabs(sggMtag(ni,nj,nk))
                        endif
                        if (asigna) output( ii)%item( i)%Serialized%valor(Ntimeforvolumic,conta) = 9.0
                     endif
                  end do
               End do
            end do
         end do
         !
         do sweep=1,rNodal_Hy%numSoft
            do nk=rNodal_Hy%nodSoft(sweep)%punto%zi,rNodal_Hy%nodSoft(sweep)%punto%ze
               k_m = nk - b%Hy%ZI
               do nj=rNodal_Hy%nodSoft(sweep)%punto%yi,rNodal_Hy%nodSoft(sweep)%punto%ye
                  j_m = nj - b%Hy%YI
                  do ni=rNodal_Hy%nodSoft(sweep)%punto%xi,rNodal_Hy%nodSoft(sweep)%punto%xe
                     i_m = ni - b%Hy%XI
                     imed = sggMiHy(i_m,j_m,k_m)
                     if (.not.sgg%Med(imed)%Is%PMC) then
                        conta=conta+1
                        if (geom) then
                           output(ii)%item(i)%Serialized%eI(conta)=ni
                           output(ii)%item(i)%Serialized%eJ(conta)=nj
                           output(ii)%item(i)%Serialized%eK(conta)=nk
                           output(ii)%item(i)%Serialized%currentType(conta)=iBloqueJy
                           output(ii)%item(i)%Serialized%sggMtag(conta)=iabs(sggMtag(ni,nj,nk))
                        endif
                        if (asigna) output( ii)%item( i)%Serialized%valor(Ntimeforvolumic,conta) = 9.0
                     endif
                  end do
               End do
            end do
         end do

         do sweep=1,rNodal_Hz%numHard
            do nk=rNodal_Hz%nodHard(sweep)%punto%zi,rNodal_Hz%nodHard(sweep)%punto%ze
               k_m = nk - b%Hz%ZI
               do nj=rNodal_Hz%nodHard(sweep)%punto%yi,rNodal_Hz%nodHard(sweep)%punto%ye
                  j_m = nj - b%Hz%YI
                  do ni=rNodal_Hz%nodHard(sweep)%punto%xi,rNodal_Hz%nodHard(sweep)%punto%xe
                     i_m = ni - b%Hz%XI
                     imed = sggMiHx(i_m,j_m,k_m)
                     if (.not.sgg%Med(imed)%Is%PMC) then
                        conta=conta+1
                        if (geom) then
                           output(ii)%item(i)%Serialized%eI(conta)=ni
                           output(ii)%item(i)%Serialized%eJ(conta)=nj
                           output(ii)%item(i)%Serialized%eK(conta)=nk
                           output(ii)%item(i)%Serialized%currentType(conta)=iBloqueJz
                           output(ii)%item(i)%Serialized%sggMtag(conta)=iabs(sggMtag(ni,nj,nk))
                        endif
                        if (asigna) output( ii)%item( i)%Serialized%valor(Ntimeforvolumic,conta) = 9.0
                     endif
                  end do
               End do
            end do
         end do
         !
         do sweep=1,rNodal_Hz%numSoft
            do nk=rNodal_Hz%nodSoft(sweep)%punto%zi,rNodal_Hz%nodSoft(sweep)%punto%ze
               k_m = nk - b%Hz%ZI
               do nj=rNodal_Hz%nodSoft(sweep)%punto%yi,rNodal_Hz%nodSoft(sweep)%punto%ye
                  j_m = nj - b%Hz%YI
                  do ni=rNodal_Hz%nodSoft(sweep)%punto%xi,rNodal_Hz%nodSoft(sweep)%punto%xe
                     i_m = ni - b%Hz%XI
                     imed = sggMiHz(i_m,j_m,k_m)
                     if (.not.sgg%Med(imed)%Is%PMC) then
                        conta=conta+1
                        if (geom) then
                           output(ii)%item(i)%Serialized%eI(conta)=ni
                           output(ii)%item(i)%Serialized%eJ(conta)=nj
                           output(ii)%item(i)%Serialized%eK(conta)=nk
                           output(ii)%item(i)%Serialized%currentType(conta)=iBloqueJz
                           output(ii)%item(i)%Serialized%sggMtag(conta)=iabs(sggMtag(ni,nj,nk))
                        endif
                        if (asigna) output( ii)%item( i)%Serialized%valor(Ntimeforvolumic,conta) = 9.0
                     endif
                  end do
               End do
            end do
         end do
      ENDIF


      !!!!!!

      !print *,'----tras nodalvtk init,geom,asigna,electric,magnetic,conta,i,ii ',init,geom,asigna,electric,magnetic,conta,i,ii
      return
   end subroutine

#endif
!del CompileWithNodalSources


#ifdef CompileWithWires
   subroutine wirebundlesvtk(sgg,init,geom,asigna,conta,i,ii,output,Ntimeforvolumic,wiresflavor,sggMtag)
   
      type (SGGFDTDINFO), intent(IN)   :: sgg
      INTEGER (KIND=IKINDMTAG), intent(in) :: sggMtag  (sgg%Alloc(iHx)%XI:sgg%Alloc(iHx)%XE, sgg%Alloc(iHy)%YI:sgg%Alloc(iHy)%YE, sgg%Alloc(iHz)%ZI:sgg%Alloc(iHz)%ZE)
      type (output_t), pointer, dimension( : )  ::  output
      integer(kind=4), intent(IN) :: i,ii,Ntimeforvolumic
      logical, intent( IN) :: geom,asigNa,init
      integer(kind=4) conta,ni,nj,nk,n
      character(len=*), INTENT(in) :: wiresflavor
      integer(kind=4), SAVE :: MINIMED


      type(Thinwires_t), pointer, save  ::  Hwireslocal
      !
#ifdef CompileWithBerengerWires
      type(TWires)     , pointer, save  ::  Hwireslocal_Berenger
#endif
#ifdef CompileWithSlantedWires
      !
      type(WiresData)  , pointer, save  ::  Hwireslocal_Slanted
#endif


      !print *,'----antes wires init,geom,asigna,conta,i,ii',init,geom,asigna,conta,i,ii
      if (init) then
#ifdef CompileWithWires
         if ((trim(adjustl(wiresflavor))=='holland') .or. &
             (trim(adjustl(wiresflavor))=='transition')) then
            Hwireslocal => GetHwires()
         endif
#endif
#ifdef CompileWithBerengerWires
         if (trim(adjustl(wiresflavor))=='berenger') then
            Hwireslocal_Berenger => GetHwires_Berenger()
         endif
#endif
#ifdef CompileWithSlantedWires
         if((trim(adjustl(wiresflavor))=='slanted').or.(trim(adjustl(wiresflavor))=='semistructured')) then
            Hwireslocal_Slanted => GetHwires_Slanted()
         endif
#endif
      endif
#ifdef CompileWithBerengerWires
      if (trim(adjustl(wiresflavor))=='berenger') then

         !parsea los hilos
         if (geom) then
            MINIMED=2**12
            do n=1,Hwireslocal_Berenger%NumSegments
               conta=conta+1
               minimed=MIN(MINIMED, Hwireslocal_Berenger%Segments(n)%imeD)
               ni=Hwireslocal_Berenger%Segments(n)%ii
               nj=Hwireslocal_Berenger%Segments(n)%ji
               nk=Hwireslocal_Berenger%Segments(n)%ki
               select case (Hwireslocal_Berenger%Segments(n)%orient)
                case (iEx)
                  output(ii)%item(i)%Serialized%currentType(conta)=iJx
                case (iEy)
                  output(ii)%item(i)%Serialized%currentType(conta)=iJy
                case (iEz)
                  output(ii)%item(i)%Serialized%currentType(conta)=iJz
               end select
               output(ii)%item(i)%Serialized%eI(conta)=ni
               output(ii)%item(i)%Serialized%eJ(conta)=nj
               output(ii)%item(i)%Serialized%eK(conta)=nk
               output(ii)%item(i)%Serialized%sggMtag(conta)=iabs(sggMtag(ni,nj,nk))
            end do

         elseif (asigna) then
            do n=1,Hwireslocal_Berenger%NumSegments
               conta=conta+1
               if ((Hwireslocal_Berenger%Segments(n)%Is_LeftEnd).or. &
               (Hwireslocal_Berenger%Segments(n)%Is_RightEnd)) then
                  output( ii)%item( i)%Serialized%valor(Ntimeforvolumic,conta) = 10
               else
                  output( ii)%item( i)%Serialized%valor(Ntimeforvolumic,conta) = 20 + Hwireslocal_Berenger%Segments(n)%imed-MINIMED
               endif
            end do
         else
            do n=1,Hwireslocal_Berenger%NumSegments
               conta=conta+1
            end do
         endif
      endif
#endif
#ifdef CompileWithWires
      if ((trim(adjustl(wiresflavor))=='holland') .or. &
              (trim(adjustl(wiresflavor))=='transition')) then
         if (geom) then
            MINIMED=2**30
            do n=1,Hwireslocal%NumCurrentSegments
               conta=conta+1
               minimed=MIN(MINIMED, Hwireslocal%CurrentSegment(n)%indexmed)
               ni=Hwireslocal%CurrentSegment(n)%i
               nj=Hwireslocal%CurrentSegment(n)%j
               nk=Hwireslocal%CurrentSegment(n)%k
               select case (Hwireslocal%CurrentSegment(n)%tipofield)
                case (iEx)
                  output(ii)%item(i)%Serialized%currentType(conta)=iJx
                case (iEy)
                  output(ii)%item(i)%Serialized%currentType(conta)=iJy
                case (iEz)
                  output(ii)%item(i)%Serialized%currentType(conta)=iJz
               end select
               output(ii)%item(i)%Serialized%eI(conta)=ni
               output(ii)%item(i)%Serialized%eJ(conta)=nj
               output(ii)%item(i)%Serialized%eK(conta)=nk
               output(ii)%item(i)%Serialized%sggMtag(conta)=iabs(sggMtag(ni,nj,nk))
            end do

         elseif (asigna) then
            do n=1,Hwireslocal%NumCurrentSegments
               conta=conta+1
               if ((Hwireslocal%CurrentSegment(n)%Is_LeftEnd).or. &
               (Hwireslocal%CurrentSegment(n)%Is_RightEnd)) then
                  output( ii)%item( i)%Serialized%valor(Ntimeforvolumic,conta) = 10
               elseif (Hwireslocal%CurrentSegment(n)%IsEnd_norLeft_norRight) then
                  output( ii)%item( i)%Serialized%valor(Ntimeforvolumic,conta) = 11
               else
                  !  output( ii)%item( i)%Serialized%valor(Ntimeforvolumic,conta) = 20 + Hwireslocal%CurrentSegment(n)%indexmed-MINIMED
                  output( ii)%item( i)%Serialized%valor(Ntimeforvolumic,conta) = 20 + Hwireslocal%CurrentSegment(n)%NumParallel
               endif
            end do
         else
            do n=1,Hwireslocal%NumCurrentSegments
               conta=conta+1
            end do
         endif
      endif
#endif
#ifdef CompileWithSlantedWires
      if ((trim(adjustl(wiresflavor))=='slanted').or.(trim(adjustl(wiresflavor))=='semistructured')) then
         !parsea los hilos
         if (geom) then
            MINIMED=2**12
            do n=1,Hwireslocal_Slanted%NumSegmentsStr
               conta=conta+1
               minimed=MIN(MINIMED, Hwireslocal_Slanted%SegmentsStr(n)%imeD)
               ni=Hwireslocal_Slanted%SegmentsStr(n)%cell(iX)
               nj=Hwireslocal_Slanted%SegmentsStr(n)%cell(iY)
               nk=Hwireslocal_Slanted%SegmentsStr(n)%cell(iZ)
               select case (Hwireslocal_Slanted%SegmentsStr(n)%orient)
                case (iEx)
                  output(ii)%item(i)%Serialized%currentType(conta)=iJx
                case (iEy)
                  output(ii)%item(i)%Serialized%currentType(conta)=iJy
                case (iEz)
                  output(ii)%item(i)%Serialized%currentType(conta)=iJz
               end select
               output(ii)%item(i)%Serialized%eI(conta)=ni
               output(ii)%item(i)%Serialized%eJ(conta)=nj
               output(ii)%item(i)%Serialized%eK(conta)=nk
               output(ii)%item(i)%Serialized%sggMtag(conta)=iabs(sggMtag(ni,nj,nk))
            end do

         elseif (asigna) then
            do n=1,Hwireslocal_Slanted%NumSegmentsStr
               conta=conta+1
               if ((Hwireslocal_Slanted%SegmentsStr(n)%IsExt(iBeg)).or. &
                   (Hwireslocal_Slanted%SegmentsStr(n)%IsExt(iEnd))) then
                  output( ii)%item( i)%Serialized%valor(Ntimeforvolumic,conta) = 10
               else
                  output( ii)%item( i)%Serialized%valor(Ntimeforvolumic,conta) = 20 + Hwireslocal_Slanted%SegmentsStr(n)%imed-MINIMED
               endif
            end do
         else
            do n=1,Hwireslocal_Slanted%NumSegmentsStr
               conta=conta+1
            end do
         endif
    
      endif
#endif     
      
      !print *,'----tras wires init,geom,asigna,conta,i,ii',init,geom,asigna,conta,i,ii

      return
   end subroutine

#endif


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !Function to publish the private output data (used in postprocess)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!

   function GetOutput() result(r)
      type(output_t), pointer, dimension( : )  ::  r

      r => output
      return
   end function


   ! ========================================================================
   ! PURPOSE:
   ! Performs Discrete Time Fourier Transform in a signal given in sig
   ! sampled at times st.
   ! The frequency values are stored in fqVal for frequencies given in fq.
   ! This subroutine is efficient when fqSize << sigSize.
   ! ========================================================================
   subroutine dtft(fqVal, fq, fqSize, st, sig, sigSize)
      implicit none
      ! ----------- Variables --------------------------------------------------
      ! In & out vars.
      ! Size of the input frequencies and signal vectors.
      integer, intent(in) :: fqSize, sigSize
      ! Fourier transform value.
      complex (kind=CKIND), intent(out), dimension(fqSize) :: fqVal
      ! Vector of frequencies to compute the values.
      real (kind=RKIND), intent(in), dimension(fqSize) ::  fq
      ! Input signal.
      real (kind=RKIND), intent(in), dimension(sigSize) :: sig
      ! Input signal sampling time.
      real (kind=RKIND_tiempo), intent(in), dimension(sigSize) :: st
      ! Aux. vars.
      integer :: i, j;
      ! ----------- Body -------------------------------------------------------
      ! Checks that frequencies are below the sf.
      fqval=0.0_RKIND
      
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j)
#endif
      do i=1, fqSize
         do j = 1, sigSize-1 !algun delta promedio habra que tomar permit scaling 211118
            fqVal(i) = fqVal(i) + abs(st(j+1)-st(j))/2.0_RKIND * sig(j) * exp(mcPI2 * fq(i) * st(j)); !nosisosampleada 200120 !ojo valor absoluto delta 
!!            fqVal(i) = fqVal(i) + abs(st(2)-st(1))/2.0_RKIND * sig(j) * exp(mcPI2 * fq(i) * (j-1)*(st(2)-st(1))); !isosampleada 200120 !ojo valor absoluto delta    
!!            if ((i==1).and.(j==sigsize-1)) write (634,*) j,abs(st(2)-st(1)),abs(st(j+1)-st(j))
         end do
      end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif

#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j)
#endif
      do i=1, fqSize
         j=sigSize  !algun delta promedio habra que tomar permit scaling 211118
         !no debia mucho tener impacto por ser el ultimo timpicamente pequenio, pero....         
         fqVal(i) = fqVal(i) + abs(st(j-1)-st(j))/2.0_RKIND * sig(j) * exp(mcPI2 * fq(i) * st(j)); !nosisosampleada 200120 !ojo valor absoluto delta  !ojo valor absoluto delta pq el ultimo cambiaba elsigno 14111' 
!!         fqVal(i) = fqVal(i) + abs(st(2)-st(1))/2.0_RKIND * sig(j) * exp(mcPI2 * fq(i) * (j-1)*(st(2)-st(1))); !isosampleada 200120 !ojo valor absoluto delta  
     
      end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif


      fqVal=fqVal*2.0_RKIND !BUG HIRAI ENERGIA DOBLE PARSEVAL  mail 24/07/19

   end subroutine


end module Observa
