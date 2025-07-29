

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This module defines a class to deal with the transfering information
! from the parser to a class where it can be stored into memory and
! used in any further development of the FDTD simulation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!#define DetectAdj              : call the calculus of
!!!!!!!!!                                 since the routine is not deeply tested)
MODULE Preprocess_m
#undef DetectAdj
   !
   USE Report
   USE NFDETypes
   !healer sgg10
   USE CreateMatrices
   !typos que leo desde mi FDE
   USE FDEtypes
   USE DMMA
#ifdef CompileWithConformal
   USE CONFORMAL_INI_CLASS
   USE CONFORMAL_TOOLS
   USE CONFORMAL_MAPPED
   USE CONFORMAL_TYPES
   USE Conformal_TimeSteps_m
#endif
   USE conformal_mod, F_X => FACE_X, F_Y => FACE_Y, F_Z => FACE_Z, E_X => EDGE_X, E_Y => EDGE_Y, E_Z => EDGE_Z
   IMPLICIT NONE
!!!variables globales del modulo
   REAL (KIND=RKIND), save           ::  cluz,zvac
   REAL (KIND=RKIND), save           ::  eps0,mu0
!!!
   PRIVATE
   !
   PUBLIC read_geomData, read_limits_nogeom,AssigLossyOrPECtoNodes
   !
CONTAINS
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE read_geomData (sgg,media,tag_numbers, fichin, layoutnumber, size, SINPML_fullsize, fullsize, this, &
      groundwires,attfactor,mibc,SGBC,SGBCDispersive,MEDIOEXTRA,maxSourceValue,skindepthpre,createmapvtk,input_conformal_flag,CLIPREGION,boundwireradius,maxwireradius,updateshared,run_with_dmma, &
      eps00,mu00,simu_devia,hay_slanted_wires,verbose,ignoresamplingerrors,tagtype,wiresflavor)
      type(media_matrices_t), intent(inout) :: media
      logical :: simu_devia,verbose,hay_slanted_wires
      REAL (KIND=RKIND)           ::  eps00,mu00

      TYPE (MedioExtra_t), INTENT (INout) :: MEDIOEXTRA
      !
      CHARACTER (LEN=BUFSIZE), intent(in) :: wiresflavor
      logical, intent (in) :: updateshared,run_with_dmma,ignoresamplingerrors
      LOGICAL, INTENT (INout) :: mibc,SGBC,CLIPREGION,boundwireradius,SGBCDispersive,skindepthpre
      LOGICAL, INTENT (INout) :: createmapvtk
      TYPE (limit_t), DIMENSION (1:6) :: SINPML_fullsize, fullsize
      type (SGGFDTDINFO), intent(INOUT)    :: sgg
      character(len=BUFSIZE) :: extraswitches

      type(taglist_t) :: tag_numbers
      TYPE (Parseador), INTENT (INOUT) :: this
      INTEGER (KIND=4) :: tama, tama2, tama3, tama4, tama5, tama6, i, j, k, tipotemp, tamaSonda,  &
      &      tamaoldSONDA, tamaBloquePrb, tamaScrPrb,pozi,tama2bis,numeroasignaciones,ci
      CHARACTER (LEN=*), INTENT (IN) :: fichin
      !
      CHARACTER (LEN=BUFSIZE) :: probenumber
      REAL (KIND=RKIND) :: ex, ey, ez, px, py, pz, amplitud,attfactor,maxSourceValue,minSpaceStep,maxwireradius

      CHARACTER (LEN=BUFSIZE) :: tag

      TYPE (XYZlimit_t) :: punto, BoundingBox
      TYPE (XYZlimit_t_scaled) :: punto_s
      INTEGER (KIND=4) :: orientacion,orientacionL,orientacionR, direccion, contamedia,oldcontamedia, maxcontamedia, mincontamedia, inicontamedia, &
         i1, j1, field, k1, pecmedio, ii, medio1, medio2, sondas,CONTACURR,CONTAVOLT,I_,J_
      !
      LOGICAL ::  isathinwire, VALIDO, existia,medioespecial,input_conformal_flag,nodo_cazado
      LOGICAL :: errnofile,errnofile1,errnofile2,errnofile3,errnofile4
      REAL (KIND=RKIND) :: tiempo1, tiempo2, field1, field2,rdummy
      INTEGER (KIND=4) :: nsurfs, numus, layoutnumber, size,OrigIndex,numminus
      REAL (KIND=RKIND) :: delta,del,sig_max
      INTEGER (KIND=4), DIMENSION (:), ALLOCATABLE :: contapuntos
      INTEGER (KIND=4) :: conta1, conta2, MEDIO,imenos1,jmenos1,kmenos1,o,p,puntoxi,puntoyi,puntozi, &
         bboxwirXI,dummy_bboxwirXI,bboxwirYI,dummy_bboxwirYI,bboxwirzI,dummy_bboxwirzI, &
         bboxwirXE,dummy_bboxwirXE,bboxwirYE,dummy_bboxwirYE,bboxwirZE,dummy_bboxwirZE,IERR
      INTEGER (KIND=8) :: memo
      CHARACTER (LEN=BUFSIZE) :: MultiportFile,MultiportFile2
      CHARACTER (LEN=BUFSIZE) :: buff
      REAL (KIND=RKIND), allocatable, dimension (:) :: dummy_px,dummy_py,dummy_pz,dummy_ex,dummy_ey,dummy_ez,dummy_INCERT
      !
      CHARACTER (LEN=BUFSIZE) :: whoami
      character (LEN=BUFSIZE)  ::  whoamishort
      character (LEN=BUFSIZE)  ::  ext,extpoint
      character (LEN=BUFSIZE)  ::  chari,charj,chark,chari2,charj2,chark2
      !
      logical :: paraerrhilo,groundwires,islossy,DENTRO
      REAL (KIND=RKIND) :: width, dir (1:3), epr1, mur1
      LOGICAL :: oriX, oriY, oriZ, oriX2, oriY2, oriZ2, oriX3, oriY3, oriZ3, iguales
      LOGICAL :: oriX4, oriY4, oriZ4
      REAL (KIND=RKIND), DIMENSION (3, 3) :: EprSlot, MurSlot
      INTEGER (KIND=4) :: indicemedio
      INTEGER (KIND=4) :: i11, j11
      !
      type (tagtype_t) :: tagtype
      TYPE (FreqDepenMaterial), POINTER :: fdgeom
      !
      INTEGER (KIND=4) :: numertag
      INTEGER (KIND=4) :: Alloc_iEx_XI, Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, &
      & Alloc_iEy_XE, Alloc_iEy_YI, Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, &
      & Alloc_iEz_YE, Alloc_iEz_ZI, Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, &
      & Alloc_iHx_ZE, Alloc_iHy_XI, Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, &
      & Alloc_iHz_XE, Alloc_iHz_YI, Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE
      !
      !
      type(ConformalMedia_t) :: conformal_media
      eps0=eps00; mu0=mu00; !chapuz para convertir la variables de paso en globales
      cluz=1.0_RKIND/sqrt(eps0*mu0)
      zvac=sqrt(mu0/eps0)
!
      call cuentatags(this,tagtype,layoutnumber,fichin)
!
      delta=-1 !para que no se queje gfortran de variables sin inicializar

      sgg%thereAreMagneticMedia=.true.  !caso mas general
      sgg%thereArePMLMagneticMedia=.true. !caso mas general
      !antes de hacer nada preprocesa si es preciso y no allocatear nada
      !09/07/13 !los SGBCs con skindepth se deben preprocesar
      if (skindepthpre) then
         if (layoutnumber == 0) then
            call print11(layoutnumber,'Preprocessing SGBC materials to include skin-depth effects....')
            call prepro_skindepth(this,fichin)
            call print11(layoutnumber,'Finished preprocessing for skin-depth.')
         endif
#ifdef CompileWithMPI
         call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
         return
      endif
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      WRITE (whoami, '(a,i5,a,i5,a)') '(', layoutnumber + 1, '/', size, ') '
      write(whoamishort,'(i5)') layoutnumber+1
      !create space for the etangential shared info
      sgg%EShared%Conta = 0
      sgg%EShared%MaxConta = 10
      ALLOCATE (sgg%EShared%elem(1:sgg%EShared%MaxConta))


      if (associated(this%conformalRegs%volumes)) then 
         conformal_media = buildConformalMedia(this%conformalRegs)
         sgg%dt = sgg%dt*conformal_media%cfl
      end if


      ! Cuenta los medios
      !!!!!calcula tamanios
      !reserva espacio
      !el medio 0 se reserva para PEC
      !regiones PEC
      !el medio 1 se reserva para sustrato  y saltamos
      contamedia = 1
      IF ((this%pmcregs%nvols)+(this%pmcregs%nsurfs)+(this%pmcregs%nLINS) /= 0) THEN
         !los PMC empiezan en 2
         contamedia = 2
         !fin regions PMC
      END IF
      !materialList
      !NonMetalREgions   and frequencydependent media
      !Anisotropic
      !
      contamedia = contamedia + (this%DielRegs%nvols) + (this%DielRegs%nsurfs) + (this%DielRegs%nLINS) + this%FRQDEPMATS%nvols + &
      & this%FRQDEPMATS%nsurfs + this%FRQDEPMATS%nLINS + (this%ANIMATS%nvols+this%ANIMATS%nsurfs+this%ANIMATS%nLINS)
      !Multiports
      !worst case 6 orientations per surface plus the the lossy padding
      contamedia = contamedia + this%LossyThinSurfs%length * 7
      !wires
      !nueva formulacion que almacena also the lenghts
      contamedia = contamedia + this%twires%n_tw
      contamedia = contamedia + this%swires%n_sw
      !echo por demas, habria que precontar pero es complicado porque depende del procesamiento
      !thin Slots

      if (run_with_dmma) then
         DO j = 1, this%tSlots%n_tg
            contamedia = contamedia + this%tSlots%Tg(j)%N_tgc
         END DO
      endif

      !end thin Slots
      !PARA LA CAPA EXTRA 2013
      if (medioextra%exists) then
         CONTAMEDIA = CONTAMEDIA+1
         MEDIOEXTRA%index=CONTAMEDIA
      endif
      !para modulos que necesiten senialar con already_YEEadvanced_byconformal y split_and_useless (eg. conformal)
      !se crea siempre por defecto
      contamedia = contamedia+2 !para acomodar los no_use no_use_notouch
      !!!!!!!!!!!!!
      contamedia = contamedia +1 !para acomodar los nodal sources como caso especial de linea vacia

      ! contamedia = contamedia + this%conformalRegs%nEdges + this%conformalRegs%nFaces
      contamedia = contamedia + conformal_media%n_edges_media + conformal_media%n_faces_media

      sgg%NumMedia = contamedia
      sgg%AllocMed = contamedia
      !reserva espacio
      ALLOCATE (sgg%Med(0:sgg%NumMedia))
      !comienzo barrido resto :  medios y observaciones
      BoundingBox%XI = sgg%Alloc(iHx)%XI
      BoundingBox%XE = sgg%Alloc(iHx)%XE
      BoundingBox%YI = sgg%Alloc(iHy)%YI
      BoundingBox%YE = sgg%Alloc(iHy)%YE
      BoundingBox%ZI = sgg%Alloc(iHz)%ZI
      BoundingBox%ZE = sgg%Alloc(iHz)%ZE
      !
      Alloc_iEx_XI = sgg%Alloc(iEx)%XI
      Alloc_iEx_XE = sgg%Alloc(iEx)%XE
      Alloc_iEx_YI = sgg%Alloc(iEx)%YI
      Alloc_iEx_YE = sgg%Alloc(iEx)%YE
      Alloc_iEx_ZI = sgg%Alloc(iEx)%ZI
      Alloc_iEx_ZE = sgg%Alloc(iEx)%ZE
      Alloc_iEy_XI = sgg%Alloc(iEy)%XI
      Alloc_iEy_XE = sgg%Alloc(iEy)%XE
      Alloc_iEy_YI = sgg%Alloc(iEy)%YI
      Alloc_iEy_YE = sgg%Alloc(iEy)%YE
      Alloc_iEy_ZI = sgg%Alloc(iEy)%ZI
      Alloc_iEy_ZE = sgg%Alloc(iEy)%ZE
      Alloc_iEz_XI = sgg%Alloc(iEz)%XI
      Alloc_iEz_XE = sgg%Alloc(iEz)%XE
      Alloc_iEz_YI = sgg%Alloc(iEz)%YI
      Alloc_iEz_YE = sgg%Alloc(iEz)%YE
      Alloc_iEz_ZI = sgg%Alloc(iEz)%ZI
      Alloc_iEz_ZE = sgg%Alloc(iEz)%ZE
      Alloc_iHx_XI = sgg%Alloc(iHx)%XI
      Alloc_iHx_XE = sgg%Alloc(iHx)%XE
      Alloc_iHx_YI = sgg%Alloc(iHx)%YI
      Alloc_iHx_YE = sgg%Alloc(iHx)%YE
      Alloc_iHx_ZI = sgg%Alloc(iHx)%ZI
      Alloc_iHx_ZE = sgg%Alloc(iHx)%ZE
      Alloc_iHy_XI = sgg%Alloc(iHy)%XI
      Alloc_iHy_XE = sgg%Alloc(iHy)%XE
      Alloc_iHy_YI = sgg%Alloc(iHy)%YI
      Alloc_iHy_YE = sgg%Alloc(iHy)%YE
      Alloc_iHy_ZI = sgg%Alloc(iHy)%ZI
      Alloc_iHy_ZE = sgg%Alloc(iHy)%ZE
      Alloc_iHz_XI = sgg%Alloc(iHz)%XI
      Alloc_iHz_XE = sgg%Alloc(iHz)%XE
      Alloc_iHz_YI = sgg%Alloc(iHz)%YI
      Alloc_iHz_YE = sgg%Alloc(iHz)%YE
      Alloc_iHz_ZI = sgg%Alloc(iHz)%ZI
      Alloc_iHz_ZE = sgg%Alloc(iHz)%ZE
      !
      !
      field = 1
      !
      numertag = 0
      ALLOCATE (media%sggMtag(Alloc_iHx_XI:Alloc_iHx_XE, Alloc_iHy_YI:Alloc_iHy_YE, Alloc_iHz_ZI:Alloc_iHz_ZE))
      ALLOCATE (media%sggMiNo(Alloc_iHx_XI:Alloc_iHx_XE, Alloc_iHy_YI:Alloc_iHy_YE, Alloc_iHz_ZI:Alloc_iHz_ZE))

      ALLOCATE (tag_numbers%edge%x(Alloc_iEx_XI:Alloc_iEx_XE, Alloc_iEx_YI:Alloc_iEx_YE, Alloc_iEx_ZI:Alloc_iEx_ZE))
      ALLOCATE (tag_numbers%edge%y(Alloc_iEy_XI:Alloc_iEy_XE, Alloc_iEy_YI:Alloc_iEy_YE, Alloc_iEy_ZI:Alloc_iEy_ZE))
      ALLOCATE (tag_numbers%edge%z(Alloc_iEz_XI:Alloc_iEz_XE, Alloc_iEz_YI:Alloc_iEz_YE, Alloc_iEz_ZI:Alloc_iEz_ZE))
      ALLOCATE (tag_numbers%face%x(Alloc_iHx_XI:Alloc_iHx_XE, Alloc_iHx_YI:Alloc_iHx_YE, Alloc_iHx_ZI:Alloc_iHx_ZE))
      ALLOCATE (tag_numbers%face%y(Alloc_iHy_XI:Alloc_iHy_XE, Alloc_iHy_YI:Alloc_iHy_YE, Alloc_iHy_ZI:Alloc_iHy_ZE))
      ALLOCATE (tag_numbers%face%z(Alloc_iHz_XI:Alloc_iHz_XE, Alloc_iHz_YI:Alloc_iHz_YE, Alloc_iHz_ZI:Alloc_iHz_ZE))

      !!!nodos materiales: se precisan para el conformal !sgg310715
      ALLOCATE (media%sggMiEx(Alloc_iEx_XI:Alloc_iEx_XE, Alloc_iEx_YI:Alloc_iEx_YE, Alloc_iEx_ZI:Alloc_iEx_ZE))
      ALLOCATE (media%sggMiEy(Alloc_iEy_XI:Alloc_iEy_XE, Alloc_iEy_YI:Alloc_iEy_YE, Alloc_iEy_ZI:Alloc_iEy_ZE))
      ALLOCATE (media%sggMiEz(Alloc_iEz_XI:Alloc_iEz_XE, Alloc_iEz_YI:Alloc_iEz_YE, Alloc_iEz_ZI:Alloc_iEz_ZE))
      ALLOCATE (media%sggMiHx(Alloc_iHx_XI:Alloc_iHx_XE, Alloc_iHx_YI:Alloc_iHx_YE, Alloc_iHx_ZI:Alloc_iHx_ZE))
      ALLOCATE (media%sggMiHy(Alloc_iHy_XI:Alloc_iHy_XE, Alloc_iHy_YI:Alloc_iHy_YE, Alloc_iHy_ZI:Alloc_iHy_ZE))
      ALLOCATE (media%sggMiHz(Alloc_iHz_XI:Alloc_iHz_XE, Alloc_iHz_YI:Alloc_iHz_YE, Alloc_iHz_ZI:Alloc_iHz_ZE))
      !el tag esta voided porque luego el numero va con el del tag
      media%sggMtag (:, :, :) = 0 !LO VOIDEO A 0 EN VEZ DE A -1 PORQUE EL TAG 0 NO VA A EXISTIR NUNCA 141020
      tag_numbers%edge%x(:,:,:) = 0
      tag_numbers%edge%y(:,:,:) = 0
      tag_numbers%edge%z(:,:,:) = 0
      tag_numbers%face%x(:,:,:) = 0
      tag_numbers%face%y(:,:,:) = 0
      tag_numbers%face%z(:,:,:) = 0
      !todo sustrato por defecto
      media%sggMiNo (:, :, :) = 1
      media%sggMiEx (:, :, :) = 1
      media%sggMiEy (:, :, :) = 1
      media%sggMiEz (:, :, :) = 1
      media%sggMiHx (:, :, :) = 1
      media%sggMiHy (:, :, :) = 1
      media%sggMiHz (:, :, :) = 1

      !planeWaves
      !
      tama = (this%plnSrc%nc)
!!!      WRITE (buff,*) 'More than 1 Huygens box unsupported'
!!!      IF (tama > 1) CALL STOPONERROR(layoutnumber,size,buff)
      !LO PONGO A MANO ojo
      amplitud = 1.0_RKIND
      sgg%NumPlaneWaves = tama
      ALLOCATE (sgg%PlaneWave(1:sgg%NumPlaneWaves))
      DO i = 1, sgg%NumPlaneWaves
         punto%XI = Min (this%plnSrc%collection(i)%coor1(1), this%plnSrc%collection(i)%coor2(1))
         punto%XE = Max (this%plnSrc%collection(i)%coor1(1), this%plnSrc%collection(i)%coor2(1))
         punto%YI = Min (this%plnSrc%collection(i)%coor1(2), this%plnSrc%collection(i)%coor2(2))
         punto%YE = Max (this%plnSrc%collection(i)%coor1(2), this%plnSrc%collection(i)%coor2(2))
         punto%ZI = Min (this%plnSrc%collection(i)%coor1(3), this%plnSrc%collection(i)%coor2(3))
         punto%ZE = Max (this%plnSrc%collection(i)%coor1(3), this%plnSrc%collection(i)%coor2(3))
         !just for the sake of peace of my mind
         !readjust Huygens surface CLEARLY in/out in case of coincidente
         IF ((punto%XI == SINPML_fullsize(iHx)%XI)) THEN
            punto%XI = SINPML_fullsize(iHx)%XI - 5
         END IF
         IF ((punto%XE == SINPML_fullsize(iHx)%XE)) THEN
            punto%XE = SINPML_fullsize(iHx)%XE + 5
         END IF
         IF ((punto%YI == SINPML_fullsize(iHy)%YI)) THEN
            punto%YI = SINPML_fullsize(iHy)%YI - 5
         END IF
         IF ((punto%YE == SINPML_fullsize(iHy)%YE)) THEN
            punto%YE = SINPML_fullsize(iHy)%YE + 5
         END IF
         IF ((punto%ZI == SINPML_fullsize(iHz)%ZI)) THEN
            punto%ZI = SINPML_fullsize(iHz)%ZI - 5
         END IF
         IF ((punto%ZE == SINPML_fullsize(iHz)%ZE)) THEN
            punto%ZE = SINPML_fullsize(iHz)%ZE + 5
         END IF
         !
         sgg%PlaneWave(i)%isRC    = this%plnSrc%collection(i)%isRC
         sgg%PlaneWave(i)%numModes= this%plnSrc%collection(i)%numModes
         sgg%PlaneWave(i)%incertmax= this%plnSrc%collection(i)%incertmax
         if (sgg%PlaneWave(i)%isRC) then
            allocate (sgg%PlaneWave(i)%px(1:sgg%PlaneWave(i)%numModes))
            allocate (sgg%PlaneWave(i)%py(1:sgg%PlaneWave(i)%numModes))
            allocate (sgg%PlaneWave(i)%pz(1:sgg%PlaneWave(i)%numModes))
            allocate (sgg%PlaneWave(i)%ex(1:sgg%PlaneWave(i)%numModes))
            allocate (sgg%PlaneWave(i)%ey(1:sgg%PlaneWave(i)%numModes))
            allocate (sgg%PlaneWave(i)%ez(1:sgg%PlaneWave(i)%numModes))
            allocate (sgg%PlaneWave(i)%INCERT(1:sgg%PlaneWave(i)%numModes))
            allocate (dummy_px(1:sgg%PlaneWave(i)%numModes))
            allocate (dummy_py(1:sgg%PlaneWave(i)%numModes))
            allocate (dummy_pz(1:sgg%PlaneWave(i)%numModes))
            allocate (dummy_ex(1:sgg%PlaneWave(i)%numModes))
            allocate (dummy_ey(1:sgg%PlaneWave(i)%numModes))
            allocate (dummy_ez(1:sgg%PlaneWave(i)%numModes))
            allocate (dummy_INCERT(1:sgg%PlaneWave(i)%numModes))
            sgg%PlaneWave(i)%px=0.0
            sgg%PlaneWave(i)%py=0.0
            sgg%PlaneWave(i)%pz=0.0
            sgg%PlaneWave(i)%ex=0.0
            sgg%PlaneWave(i)%ey=0.0
            sgg%PlaneWave(i)%ez=0.0
            sgg%PlaneWave(i)%INCERT=0.0
            if (layoutnumber==0) call populatePlaneWaveRC(sgg%PlaneWave(i),simu_devia) !only the master populates
#ifdef CompileWithMPI
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            call MPI_AllReduce( sgg%PlaneWave(i)%px, dummy_px, sgg%PlaneWave(i)%numModes, REALSIZE, MPI_SUM, MPI_COMM_WORLD, ierr)
            call MPI_AllReduce( sgg%PlaneWave(i)%py, dummy_py, sgg%PlaneWave(i)%numModes, REALSIZE, MPI_SUM, MPI_COMM_WORLD, ierr)
            call MPI_AllReduce( sgg%PlaneWave(i)%pz, dummy_pz, sgg%PlaneWave(i)%numModes, REALSIZE, MPI_SUM, MPI_COMM_WORLD, ierr)
            call MPI_AllReduce( sgg%PlaneWave(i)%ex, dummy_ex, sgg%PlaneWave(i)%numModes, REALSIZE, MPI_SUM, MPI_COMM_WORLD, ierr)
            call MPI_AllReduce( sgg%PlaneWave(i)%ey, dummy_ey, sgg%PlaneWave(i)%numModes, REALSIZE, MPI_SUM, MPI_COMM_WORLD, ierr)
            call MPI_AllReduce( sgg%PlaneWave(i)%ez, dummy_ez, sgg%PlaneWave(i)%numModes, REALSIZE, MPI_SUM, MPI_COMM_WORLD, ierr)
            call MPI_AllReduce( sgg%PlaneWave(i)%INCERT, dummy_INCERT, sgg%PlaneWave(i)%numModes, REALSIZE, MPI_SUM, MPI_COMM_WORLD, ierr)
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            sgg%PlaneWave(i)%px=dummy_px
            sgg%PlaneWave(i)%py=dummy_py
            sgg%PlaneWave(i)%pz=dummy_pz
            sgg%PlaneWave(i)%ex=dummy_ex
            sgg%PlaneWave(i)%ey=dummy_ey
            sgg%PlaneWave(i)%ez=dummy_ez
            sgg%PlaneWave(i)%INCERT=dummy_INCERT
#endif
            deallocate (dummy_px)
            deallocate (dummy_py)
            deallocate (dummy_pz)
            deallocate (dummy_ex)
            deallocate (dummy_ey)
            deallocate (dummy_ez)
            deallocate (dummy_INCERT)
         else
            sgg%PlaneWave(i)%numModes=1
            allocate (sgg%PlaneWave(i)%px(1:sgg%PlaneWave(i)%numModes))
            allocate (sgg%PlaneWave(i)%py(1:sgg%PlaneWave(i)%numModes))
            allocate (sgg%PlaneWave(i)%pz(1:sgg%PlaneWave(i)%numModes))
            allocate (sgg%PlaneWave(i)%ex(1:sgg%PlaneWave(i)%numModes))
            allocate (sgg%PlaneWave(i)%ey(1:sgg%PlaneWave(i)%numModes))
            allocate (sgg%PlaneWave(i)%ez(1:sgg%PlaneWave(i)%numModes))
            allocate (sgg%PlaneWave(i)%INCERT(1:sgg%PlaneWave(i)%numModes))
            !
            ez = amplitud * Cos (this%plnSrc%collection(i)%alpha)
            ey = amplitud * Sin (this%plnSrc%collection(i)%alpha) * Sin (this%plnSrc%collection(i)%beta)
            ex = amplitud * Sin (this%plnSrc%collection(i)%alpha) * Cos (this%plnSrc%collection(i)%beta)
            pz = Cos (this%plnSrc%collection(i)%theta)
            py = Sin (this%plnSrc%collection(i)%theta) * Sin (this%plnSrc%collection(i)%phi)
            px = Sin (this%plnSrc%collection(i)%theta) * Cos (this%plnSrc%collection(i)%phi)
            !ojo con estos redondeos.
            !!!IF (Abs(ex/amplitud) < 1e-4) ex = 0.0_RKIND
            !!!IF (Abs(ey/amplitud) < 1e-4) ey = 0.0_RKIND
            !!!IF (Abs(ez/amplitud) < 1e-4) ez = 0.0_RKIND
            !!!IF (Abs(px) < 1e-4) px = 0.0_RKIND
            !!!IF (Abs(py) < 1e-4) py = 0.0_RKIND
            !!!IF (Abs(pz) < 1e-4) pz = 0.0_RKIND
            IF (Abs(px*ex+py*ey+pz*ez) >= 1e-4) THEN
               WRITE (buff,*) 'NO TEM PLANEWAVE',ex,ey,ez,px,py,pz,(px*ex+py*ey+pz*ez),this%plnSrc%collection(i)%alpha,  &
                  this%plnSrc%collection(i)%beta,this%plnSrc%collection(i)%theta,this%plnSrc%collection(i)%phi
               CALL STOPONERROR(layoutnumber,size,buff)
            END IF
            !
            sgg%PlaneWave(i)%px(1) = px
            sgg%PlaneWave(i)%py(1) = py
            sgg%PlaneWave(i)%pz(1) = pz
            sgg%PlaneWave(i)%ex(1) = ex
            sgg%PlaneWave(i)%ey(1) = ey
            sgg%PlaneWave(i)%ez(1) = ez
            sgg%PlaneWave(i)%INCERT(1)=0.0_RKIND
         endif
         sgg%PlaneWave(i)%fichero%name = trim (adjustl(this%plnSrc%collection(i)%nombre_fichero))
         sgg%PlaneWave(i)%esqx1 = Min (punto%XI, punto%XE)
         sgg%PlaneWave(i)%esqy1 = Min (punto%YI, punto%YE)
         sgg%PlaneWave(i)%esqz1 = Min (punto%ZI, punto%ZE)
         sgg%PlaneWave(i)%esqx2 = Max (punto%XI, punto%XE)
         sgg%PlaneWave(i)%esqy2 = Max (punto%YI, punto%YE)
         sgg%PlaneWave(i)%esqz2 = Max (punto%ZI, punto%ZE)
      END DO
      !Media parsing
      !Default
      !background
      sgg%Med%Priority = prior_BV
      sgg%Med%Epr = 1.0
      sgg%Med%Sigma = 0.0
      sgg%Med%Sigmareasignado = .false. !solo afecta a un chequeo de errores en lumped 120123
      sgg%Med%Mur = 1.0
      sgg%Med%SigmaM = 0.0
      sgg%Med%Is%Interfase = .FALSE.
      sgg%Med%Is%PMLbody = .false.
      sgg%Med%Is%Needed = .TRUE.
      sgg%Med%Is%Anisotropic = .FALSE.
      sgg%Med%Is%Dielectric = .FALSE.
      sgg%Med%Is%EDispersive = .FALSE.
      sgg%Med%Is%EDispersiveAnis = .FALSE.
      sgg%Med%Is%MDispersive = .FALSE.
      sgg%Med%Is%MDispersiveAnis = .FALSE.
      sgg%Med%Is%Lumped = .FALSE.
      sgg%Med%Is%SGBC = .FALSE.
      sgg%Med%Is%SGBCDispersive = .FALSE.
      sgg%Med%Is%Lossy = .FALSE.
      sgg%Med%Is%multiport = .FALSE.
      sgg%Med%Is%multiportpadding = .FALSE.
      sgg%Med%Is%AnisMultiport = .FALSE.
      sgg%Med%Is%ThinWire = .FALSE.
      sgg%Med%Is%SlantedWire = .FALSE.
      sgg%Med%Is%ThinSlot = .FALSE.
      sgg%Med%Is%PEC = .FALSE.
      sgg%Med%Is%ConformalPEC = .FALSE.
      sgg%Med%Is%PMC = .FALSE.
      sgg%Med%Is%PML = .FALSE.
      sgg%Med%Is%Volume = .FALSE.
      sgg%Med%Is%Surface = .FALSE.
      sgg%Med%Is%Line = .FALSE.
      sgg%Med%Is%already_YEEadvanced_byconformal = .FALSE.
      sgg%Med%Is%split_and_useless = .FALSE.
      !ojo tocar tambien en el readjust de healing si se crean nuevos flags
      !
      !medio PEC y PML es intrascendente si es surface o volume
      !son los de prioridad mas alta y siempre contienen a sus campos tangenciales electricos
      !Background    only differences from default are needed
      sgg%Med(1)%Priority = prior_BV
      sgg%Med(1)%Epr = this%mats%mats(1)%eps / Eps0
      sgg%Med(1)%Sigma = this%mats%mats(1)%Sigma
      sgg%Med(1)%Mur = this%mats%mats(1)%mu / Mu0
      sgg%Med(1)%SigmaM = this%mats%mats(1)%SigmaM
      sgg%Med(1)%Is%Dielectric = .false. !considero el vacio como NO dielectrico '251114
      sgg%Med(1)%Is%Volume = .false.  !considero el vacio como no volumic false '251114
      !
      sgg%Med(0)%Is%PEC = .TRUE.
      sgg%Med(0)%Is%Needed = .TRUE.
      sgg%Med(0)%Priority = prior_PEC
      sgg%Med(0)%Epr = this%mats%mats(1)%eps / Eps0
      sgg%Med(0)%Sigma = 1.0e29_RKIND
      sgg%Med(0)%Mur = this%mats%mats(1)%mu / Mu0
      sgg%Med(0)%SigmaM = 0.0_RKIND
      !CAPA EXTRA
      !Background    only differences from default are needed

      if (medioextra%exists) then
         !!!!estimate in terms of percentage of the maximum PML conductivity the conductivity of the extra medium
         !!!This info is available from read_limits_nogeom
         !the calculus is taken from borderscpml.F90
         sig_max=0.0_RKIND
         do o=1,3
            do p=1,2
               if ((o == 1).and.(p == 1)) del=sgg%dx(SINPML_fullsize(iHx)%XI)
               if ((o == 1).and.(p == 2)) del=sgg%dx(SINPML_fullsize(iHx)%XE-1)
               if ((o == 2).and.(p == 1)) del=sgg%dy(SINPML_fullsize(iHy)%YI)
               if ((o == 2).and.(p == 2)) del=sgg%dy(SINPML_fullsize(iHy)%YE-1)
               if ((o == 3).and.(p == 1)) del=sgg%dz(SINPML_fullsize(iHz)%ZI)
               if ((o == 3).and.(p == 2)) del=sgg%dz(SINPML_fullsize(iHz)%ZE-1)
               if (sgg%PML%NumLayers(o,p) /= 0) then
                  if ((sgg%PML%NumLayers(o,p) == 10).or.(sgg%PML%NumLayers(o,p) == 5)) then
                     sig_max = max( sig_max , 0.8*(sgg%PML%orden(o,p)+1)/(zvac*del))
                  else
                     if (sgg%PML%CoeffReflPML(o,p)==1.0_RKIND) then
                        !realmente en el borderscpml
                        !sig_max(sig_max,-((log( 0.99999d0                 )*(sgg%PML%orden(o,p)+1))/ &
                        !    (2.0_RKIND *sqrt(Mu0/eps0)*sgg%PML%NumLayers(o,p)*del)))
                        !trampa para que entonces tome la conductividad autentica que se especifique y poder anular las PML y solo dejar capa fisica !!?!?
                        sig_max = 1.0_RKIND
                     else
                        sig_max = max(sig_max,-((log( sgg%PML%CoeffReflPML(o,p) )*(sgg%PML%orden(o,p)+1))/ &
                           (2.0_RKIND *sqrt(Mu0/eps0)*sgg%PML%NumLayers(o,p)*del)))
                     endif
                  endif
               endif
            end do
         end do
         MEDIOEXTRA%sigma = MEDIOEXTRA%sigma * sig_max !la especificacion se da en terminos de tanto por uno en la linea de comandos
         !
         sgg%Med(MEDIOEXTRA%index)%Epr = this%mats%mats(1)%eps / Eps0 !luego se machaca este valor
         sgg%Med(MEDIOEXTRA%index)%Sigma = MEDIOEXTRA%sigma !luego se machaca este valor
         sgg%Med(MEDIOEXTRA%index)%Mur = this%mats%mats(1)%mu / Mu0 !luego se machaca este valor
         sgg%Med(MEDIOEXTRA%index)%SigmaM = 0.0_RKIND !solo lo creo para las tangenciales electricas
         sgg%Med(MEDIOEXTRA%index)%Priority = prior_PEC
         sgg%Med(MEDIOEXTRA%index)%Is%Dielectric = .TRUE.
         sgg%Med(MEDIOEXTRA%index)%Is%Volume = .TRUE.
         sgg%Med(MEDIOEXTRA%index)%Is%PML = .TRUE.
      endif
      !
      !barre los medios
      !Primero todos los pec
      !PECRegions
      !volumenes
      !el medio 0 se reserva para PEC
      !regiones PEC
      !
      IF ((this%pecregs%nvols)+(this%pecregs%nsurfs)+(this%pecregs%nLINS) /= 0) THEN
         pecmedio = 0
         tama = (this%pecregs%nvols)
         !BODYes
         DO i = 1, tama
            punto%XI = this%pecregs%vols(i)%XI
            punto%XE = this%pecregs%vols(i)%XE
            punto%YI = this%pecregs%vols(i)%YI
            punto%YE = this%pecregs%vols(i)%YE
            punto%ZI = this%pecregs%vols(i)%ZI
            punto%ZE = this%pecregs%vols(i)%ZE
            numertag = searchtag(tagtype,this%pecregs%vols(i)%tag )
            CALL CreateVolumeMM (layoutnumber, media%sggMtag, tag_numbers, numertag, media%sggMiEx, media%sggMiEy, media%sggMiEz, &
            & media%sggMiHx, media%sggMiHy, media%sggMiHz,  Alloc_iEx_XI, &
            & Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI, &
            & Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI, &
            & Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI, &
            & Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI, &
            & Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg%Med, sgg%NumMedia, sgg%EShared, BoundingBox, punto, pecmedio)
         END DO
         !SURFs
         tama = (this%pecregs%nsurfs)
         DO i = 1, tama
            punto%XI = this%pecregs%surfs(i)%XI
            punto%XE = this%pecregs%surfs(i)%XE
            punto%YI = this%pecregs%surfs(i)%YI
            punto%YE = this%pecregs%surfs(i)%YE
            punto%ZI = this%pecregs%surfs(i)%ZI
            punto%ZE = this%pecregs%surfs(i)%ZE
            orientacion = this%pecregs%surfs(i)%or
            numertag = searchtag(tagtype,this%pecregs%surfs(i)%tag)
            CALL CreateSurfaceMM (layoutnumber, media%sggMtag, tag_numbers, numertag, media%sggMiEx, media%sggMiEy, media%sggMiEz, &
            & media%sggMiHx, media%sggMiHy, media%sggMiHz, &
            & Alloc_iEx_XI, Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, &
            & Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI, Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, &
            & Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI, Alloc_iEz_ZE, &
            & Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, &
            & Alloc_iHy_XI, Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, &
            & Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI, Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, &
            & sgg%Med, sgg%NumMedia, sgg%EShared, BoundingBox, punto, orientacion, pecmedio)
         END DO
         !LINs
         tama = (this%pecregs%nLINS)
         DO i = 1, tama
            punto%XI = this%pecregs%lins(i)%XI
            punto%XE = this%pecregs%lins(i)%XE
            punto%YI = this%pecregs%lins(i)%YI
            punto%YE = this%pecregs%lins(i)%YE
            punto%ZI = this%pecregs%lins(i)%ZI
            punto%ZE = this%pecregs%lins(i)%ZE
            orientacion = this%pecregs%lins(i)%or
            isathinwire = .FALSE.
            numertag = searchtag(tagtype,this%pecregs%lins(i)%tag)
            CALL CreateLineMM (layoutnumber, media%sggMtag, tag_numbers, numertag, media%sggMiEx, media%sggMiEy, media%sggMiEz, &
            & media%sggMiHx, media%sggMiHy, media%sggMiHz, Alloc_iEx_XI, &
            & Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI, &
            & Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI, &
            & Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI, &
            & Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI, &
            & Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg%Med, sgg%NumMedia, sgg%EShared, BoundingBox, punto, orientacion, &
            & pecmedio, isathinwire, verbose,numeroasignaciones)
         END DO
         !regiones PEC
      END IF

      !el medio 1 se reserva para sustrato  y saltamos
      contamedia = 1

      !para el conformal !debe ser tipicamente contamedia =1+1=2 pq el 0 es pec y el 1 es vacio. Ojo cambiado de sitio el PMC porque podia hacer que fuesen 3 y 4. 130220!!! y puede haber error pq por ahi se comprueba el 2 y el 3
      contamedia = contamedia + 1
      sgg%Med(contamedia)%Is%already_YEEadvanced_byconformal = .TRUE.
      !debe ser contamedia =2+1=3
      contamedia = contamedia + 1
      sgg%Med(contamedia)%Is%split_and_useless = .TRUE.

!!!!cambiado aqui 130220

      !materialList
      !regiones PMC
      IF ((this%pmcregs%nvols)+(this%pmcregs%nsurfs)+(this%pmcregs%nLINS) /= 0) THEN
         !los PMC de existir tienen todos indice 2
         contamedia =contamedia+1      !!!!contamedia = 2 !!!ufff. cambiado a 130220 por posible bug con conformal si algun dia habia regiones PMC
         sgg%Med(contamedia)%Epr = sgg%Med(1)%Epr
         sgg%Med(contamedia)%Mur = sgg%Med(1)%Mur
         sgg%Med(contamedia)%Sigma = 0.0_RKIND
         sgg%Med(contamedia)%SigmaM = 1.0e29_RKIND
         sgg%Med(contamedia)%Priority = prior_PMC
         sgg%Med(contamedia)%Is%PMC = .TRUE.
         !BODYes
         tama = (this%pmcregs%nvols)
         DO i = 1, tama
            punto%XI = this%pmcregs%vols(i)%XI
            punto%XE = this%pmcregs%vols(i)%XE
            punto%YI = this%pmcregs%vols(i)%YI
            punto%YE = this%pmcregs%vols(i)%YE
            punto%ZI = this%pmcregs%vols(i)%ZI
            punto%ZE = this%pmcregs%vols(i)%ZE
            !
            !
            numertag = searchtag(tagtype,this%pmcregs%vols(i)%tag)
            CALL CreateVolumeMM (layoutnumber, media%sggMtag, tag_numbers, numertag, media%sggMiEx, media%sggMiEy, media%sggMiEz, &
            & media%sggMiHx, media%sggMiHy, media%sggMiHz, Alloc_iEx_XI, &
            & Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI, &
            & Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI, &
            & Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI, &
            & Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI, &
            & Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg%Med, sgg%NumMedia, sgg%EShared, BoundingBox, punto, contamedia)
         END DO
         !SURFs
         tama = (this%pmcregs%nsurfs)
         DO i = 1, tama
            punto%XI = this%pmcregs%surfs(i)%XI
            punto%XE = this%pmcregs%surfs(i)%XE
            punto%YI = this%pmcregs%surfs(i)%YI
            punto%YE = this%pmcregs%surfs(i)%YE
            punto%ZI = this%pmcregs%surfs(i)%ZI
            punto%ZE = this%pmcregs%surfs(i)%ZE
            orientacion = this%pmcregs%surfs(i)%or
            numertag = searchtag(tagtype,this%pmcregs%surfs(i)%tag)
            CALL CreateSurfaceMM (layoutnumber, media%sggMtag, tag_numbers, numertag, media%sggMiEx, media%sggMiEy, media%sggMiEz, &
            & media%sggMiHx, media%sggMiHy, media%sggMiHz, Alloc_iEx_XI, &
            & Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI, &
            & Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI, &
            & Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI, &
            & Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI, &
            & Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg%Med, sgg%NumMedia, sgg%EShared, BoundingBox, punto, orientacion, &
            & contamedia)
         END DO
         !LINs
         tama = (this%pmcregs%nLINS)
         DO i = 1, tama
            punto%XI = this%pmcregs%lins(i)%XI
            punto%XE = this%pmcregs%lins(i)%XE
            punto%YI = this%pmcregs%lins(i)%YI
            punto%YE = this%pmcregs%lins(i)%YE
            punto%ZI = this%pmcregs%lins(i)%ZI
            punto%ZE = this%pmcregs%lins(i)%ZE
            orientacion = this%pmcregs%lins(i)%or
            isathinwire = .FALSE.
            numertag = searchtag(tagtype,this%pmcregs%lins(i)%tag)
            CALL CreateLineMM (layoutnumber, media%sggMtag, tag_numbers, numertag, media%sggMiEx, media%sggMiEy, media%sggMiEz, &
            & media%sggMiHx, media%sggMiHy, media%sggMiHz, Alloc_iEx_XI, &
            & Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI, &
            & Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI, &
            & Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI, &
            & Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI, &
            & Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg%Med, sgg%NumMedia, sgg%EShared, BoundingBox, punto, orientacion, &
            & contamedia, isathinwire,verbose,numeroasignaciones)
            !
         END DO
         !fin regions PMC
      END IF
!!!!fin cambiado 130220

      !NonMetalREgions
      !BODYes
      tama = (this%DielRegs%nvols)
      DO i = 1, tama
         contamedia = contamedia + 1
         sgg%Med(contamedia)%Is%Dielectric = .TRUE.
         sgg%Med(contamedia)%Priority = prior_IB
         sgg%Med(contamedia)%Epr = this%DielRegs%vols(i)%eps / Eps0
         sgg%Med(contamedia)%Sigma = this%DielRegs%vols(i)%Sigma
         sgg%Med(contamedia)%Mur =    this%DielRegs%vols(i)%mu / Mu0
         sgg%Med(contamedia)%SigmaM = this%DielRegs%vols(i)%SigmaM
!!!!pmlbody
         if (this%DielRegs%vols(i)%PMLbody) then
            sgg%Med(contamedia)%Priority = prior_pmlbody !machaca con una prioridad superior a la de thin wires y backgroud !prueba HOLD 251019 coax
            sgg%Med(contamedia)%Is%PMLbody = .true.
            ALLOCATE (sgg%Med(contamedia)%PMLbody(1))
            sgg%Med(contamedia)%PMLbody(1)%orient    = this%DielRegs%vols(i)%orient
         endif
!!!!!
         tama2 = (this%DielRegs%vols(i)%n_c2P)
         DO j = 1, tama2
            IF ((J==1).and.(this%DielRegs%vols(i)%PMLbody)) sgg%Med(contamedia)%PMLbody(1)%orient = this%DielRegs%vols(i)%c2P(j)%OR !ES IGUAL PARA TODOS
            punto%XI = this%DielRegs%vols(i)%c2P(j)%XI
            punto%XE = this%DielRegs%vols(i)%c2P(j)%XE
            punto%YI = this%DielRegs%vols(i)%c2P(j)%YI
            punto%YE = this%DielRegs%vols(i)%c2P(j)%YE
            punto%ZI = this%DielRegs%vols(i)%c2P(j)%ZI
            punto%ZE = this%DielRegs%vols(i)%c2P(j)%ZE
            numertag = searchtag(tagtype,this%DielRegs%vols(i)%c2P(j)%tag)
            CALL CreateVolumeMM (layoutnumber, media%sggMtag, tag_numbers, numertag, media%sggMiEx, media%sggMiEy, media%sggMiEz, &
            & media%sggMiHx, media%sggMiHy, media%sggMiHz, Alloc_iEx_XI, &
            & Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI, &
            & Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI, &
            & Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI, &
            & Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI, &
            & Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg%Med, sgg%NumMedia, sgg%EShared, BoundingBox, punto, contamedia)
         END DO
         tama3 = (this%DielRegs%vols(i)%n_c1P)
         DO j = 1, tama3
            IF ((J==1).and.(this%DielRegs%vols(i)%PMLbody)) sgg%Med(contamedia)%PMLbody(1)%orient = this%DielRegs%vols(i)%c1P(j)%OR !ES IGUAL PARA TODOS
            punto%XI = this%DielRegs%vols(i)%c1P(j)%XI
            punto%XE = this%DielRegs%vols(i)%c1P(j)%XI
            punto%YI = this%DielRegs%vols(i)%c1P(j)%YI
            punto%YE = this%DielRegs%vols(i)%c1P(j)%YI
            punto%ZI = this%DielRegs%vols(i)%c1P(j)%ZI
            punto%ZE = this%DielRegs%vols(i)%c1P(j)%ZI
            !
            numertag = searchtag(tagtype,this%DielRegs%vols(i)%c1P(j)%tag)
            CALL CreateVolumeMM (layoutnumber, media%sggMtag, tag_numbers, numertag, media%sggMiEx, media%sggMiEy, media%sggMiEz, &
            & media%sggMiHx, media%sggMiHy, media%sggMiHz, Alloc_iEx_XI, &
            & Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI, &
            & Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI, &
            & Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI, &
            & Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI, &
            & Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg%Med, sgg%NumMedia, sgg%EShared, BoundingBox, punto, contamedia)
         END DO
      END DO
      !SURFs
      tama = (this%DielRegs%nsurfs)
      DO i = 1, tama
         contamedia = contamedia + 1
         sgg%Med(contamedia)%Is%Dielectric = .TRUE.
         sgg%Med(contamedia)%Priority = prior_IS
         sgg%Med(contamedia)%Epr = this%DielRegs%surfs(i)%eps / Eps0
         sgg%Med(contamedia)%Sigma = this%DielRegs%surfs(i)%Sigma
         sgg%Med(contamedia)%Mur = this%DielRegs%surfs(i)%mu / Mu0
         sgg%Med(contamedia)%SigmaM = this%DielRegs%surfs(i)%SigmaM
         tama2 = (this%DielRegs%surfs(i)%n_c2P)
         DO j = 1, tama2
            punto%XI = this%DielRegs%surfs(i)%c2P(j)%XI
            punto%XE = this%DielRegs%surfs(i)%c2P(j)%XE
            punto%YI = this%DielRegs%surfs(i)%c2P(j)%YI
            punto%YE = this%DielRegs%surfs(i)%c2P(j)%YE
            punto%ZI = this%DielRegs%surfs(i)%c2P(j)%ZI
            punto%ZE = this%DielRegs%surfs(i)%c2P(j)%ZE
            orientacion = this%DielRegs%surfs(i)%c2P(j)%or
            numertag = searchtag(tagtype,this%DielRegs%surfs(i)%c2P(j)%tag)
            CALL CreateSurfaceMM (layoutnumber, media%sggMtag, tag_numbers, numertag, media%sggMiEx, media%sggMiEy, media%sggMiEz, &
            & media%sggMiHx, media%sggMiHy, media%sggMiHz, Alloc_iEx_XI, &
            & Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI, &
            & Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI, &
            & Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI, &
            & Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI, &
            & Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg%Med, sgg%NumMedia, sgg%EShared, BoundingBox, punto, orientacion, &
            & contamedia)
         END DO
         tama3 = (this%DielRegs%surfs(i)%n_c1P)

         DO j = 1, tama3
            punto%XI = this%DielRegs%surfs(i)%c1P(j)%XI
            punto%XE = this%DielRegs%surfs(i)%c1P(j)%XI
            punto%YI = this%DielRegs%surfs(i)%c1P(j)%YI
            punto%YE = this%DielRegs%surfs(i)%c1P(j)%YI
            punto%ZI = this%DielRegs%surfs(i)%c1P(j)%ZI
            punto%ZE = this%DielRegs%surfs(i)%c1P(j)%ZI
            orientacion = this%DielRegs%surfs(i)%c1P(j)%or
            numertag = searchtag(tagtype,this%DielRegs%surfs(i)%c1P(j)%tag)
            CALL CreateSurfaceMM (layoutnumber, media%sggMtag, tag_numbers, numertag, media%sggMiEx, media%sggMiEy, media%sggMiEz, &
            & media%sggMiHx, media%sggMiHy, media%sggMiHz, Alloc_iEx_XI, &
            & Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI, &
            & Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI, &
            & Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI, &
            & Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI, &
            & Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg%Med, sgg%NumMedia, sgg%EShared, BoundingBox, punto, orientacion, &
            & contamedia)
         END DO
      END DO
      !LINs
      tama = (this%DielRegs%nLINS)
      DO i = 1, tama
         numeroasignaciones=0 !solo lo usa lumped para echarselo al primer y el resto ponerlo a PEC
         contamedia = contamedia + 1
         sgg%Med(contamedia)%Is%Dielectric = .TRUE.
         sgg%Med(contamedia)%Priority = prior_IL
         sgg%Med(contamedia)%Epr = this%DielRegs%lins(i)%eps / Eps0
         sgg%Med(contamedia)%Sigma = this%DielRegs%lins(i)%Sigma
         sgg%Med(contamedia)%Mur = this%DielRegs%lins(i)%mu / Mu0
         sgg%Med(contamedia)%SigmaM = this%DielRegs%lins(i)%SigmaM
!!!!lumped
         if (this%DielRegs%lins(i)%resistor) then
            sgg%Med(contamedia)%Is%Lumped = .true.
            sgg%Med(contamedia)%Is%lossy = .true. !importante que si es lumped esto se ponga a lossy para que thin-wires haga bien el bonding !bug agb 120123 test_GGGbugresis_wire_stoch_foragasconbug
            ALLOCATE (sgg%Med(contamedia)%lumped(1))
            sgg%Med(contamedia)%lumped(1)%resistor =.true.
            sgg%Med(contamedia)%lumped(1)%inductor =.false.
            sgg%Med(contamedia)%lumped(1)%capacitor=.false.
            sgg%Med(contamedia)%lumped(1)%diodo    =.false.
            sgg%Med(contamedia)%lumped(1)%R = this%DielRegs%lins(i)%R
            sgg%Med(contamedia)%lumped(1)%L = 0.0_RKIND
            sgg%Med(contamedia)%lumped(1)%C = 0.0_RKIND
            sgg%Med(contamedia)%lumped(1)%R_devia = this%DielRegs%lins(i)%R_devia
            sgg%Med(contamedia)%lumped(1)%L_devia = 0.0_RKIND
            sgg%Med(contamedia)%lumped(1)%C_devia = 0.0_RKIND
            sgg%Med(contamedia)%lumped(1)%Rtime_on = this%DielRegs%lins(i)%Rtime_on
            sgg%Med(contamedia)%lumped(1)%Rtime_off = this%DielRegs%lins(i)%Rtime_off
            sgg%Med(contamedia)%lumped(1)%DiodB = 0.0_RKIND
            sgg%Med(contamedia)%lumped(1)%DiodIsat = 0.0_RKIND
            sgg%Med(contamedia)%lumped(1)%orient = this%DielRegs%lins(i)%DiodOri
         elseif (this%DielRegs%lins(i)%inductor) then
            sgg%Med(contamedia)%Is%Lumped = .true.
            sgg%Med(contamedia)%Is%lossy = .true. !importante que si es lumped esto se ponga a lossy para que thin-wires haga bien el bonding !bug agb 120123 test_GGGbugresis_wire_stoch_foragasconbug
            ALLOCATE (sgg%Med(contamedia)%lumped(1))
            sgg%Med(contamedia)%lumped(1)%resistor =.false.
            sgg%Med(contamedia)%lumped(1)%inductor =.true.
            sgg%Med(contamedia)%lumped(1)%capacitor=.false.
            sgg%Med(contamedia)%lumped(1)%diodo    =.false.
            sgg%Med(contamedia)%lumped(1)%R = this%DielRegs%lins(i)%R
            sgg%Med(contamedia)%lumped(1)%L = this%DielRegs%lins(i)%L
            sgg%Med(contamedia)%lumped(1)%C = 0.0_RKIND
            sgg%Med(contamedia)%lumped(1)%R_devia = this%DielRegs%lins(i)%R_devia
            sgg%Med(contamedia)%lumped(1)%L_devia = this%DielRegs%lins(i)%L_devia
            sgg%Med(contamedia)%lumped(1)%C_devia = 0.0_RKIND
            sgg%Med(contamedia)%lumped(1)%Rtime_on = 0.0 !irrelevant
            sgg%Med(contamedia)%lumped(1)%Rtime_off = 0.0 !irrelevant
            sgg%Med(contamedia)%lumped(1)%DiodB = 0.0_RKIND
            sgg%Med(contamedia)%lumped(1)%DiodIsat = 0.0_RKIND
            sgg%Med(contamedia)%lumped(1)%orient = this%DielRegs%lins(i)%DiodOri
         elseif (this%DielRegs%lins(i)%capacitor) then
            sgg%Med(contamedia)%Is%Lumped = .true.
            sgg%Med(contamedia)%Is%lossy = .true. !importante que si es lumped esto se ponga a lossy para que thin-wires haga bien el bonding !bug agb 120123 test_GGGbugresis_wire_stoch_foragasconbug
            ALLOCATE (sgg%Med(contamedia)%lumped(1))
            sgg%Med(contamedia)%lumped(1)%resistor =.false.
            sgg%Med(contamedia)%lumped(1)%inductor =.false.
            sgg%Med(contamedia)%lumped(1)%capacitor=.true.
            sgg%Med(contamedia)%lumped(1)%diodo    =.false.
            sgg%Med(contamedia)%lumped(1)%R = this%DielRegs%lins(i)%R
            sgg%Med(contamedia)%lumped(1)%L = 0.0_RKIND
            sgg%Med(contamedia)%lumped(1)%C = this%DielRegs%lins(i)%C
            sgg%Med(contamedia)%lumped(1)%R_devia = this%DielRegs%lins(i)%R_devia
            sgg%Med(contamedia)%lumped(1)%L_devia = 0.0_RKIND
            sgg%Med(contamedia)%lumped(1)%C_devia = this%DielRegs%lins(i)%C_devia
            sgg%Med(contamedia)%lumped(1)%Rtime_on = 0.0 !irrelevant
            sgg%Med(contamedia)%lumped(1)%Rtime_off = 0.0 !irrelevant
            sgg%Med(contamedia)%lumped(1)%DiodB = 0.0_RKIND
            sgg%Med(contamedia)%lumped(1)%DiodIsat = 0.0_RKIND
            sgg%Med(contamedia)%lumped(1)%orient = this%DielRegs%lins(i)%DiodOri
         elseif (this%DielRegs%lins(i)%diodo) then
!!!27/08/15 diodos aun no soportados
            WRITE (buff, '(a)')    'Lumped Diodes currently unsupported. .'
            CALL STOPONERROR(layoutnumber,size,buff)
!!!
            sgg%Med(contamedia)%Is%Lumped = .true.
            sgg%Med(contamedia)%Is%lossy = .true. !importante que si es lumped esto se ponga a lossy para que thin-wires haga bien el bonding !bug agb 120123 test_GGGbugresis_wire_stoch_foragasconbug
            ALLOCATE (sgg%Med(contamedia)%lumped(1))
            sgg%Med(contamedia)%lumped(1)%resistor =.false.
            sgg%Med(contamedia)%lumped(1)%inductor =.false.
            sgg%Med(contamedia)%lumped(1)%capacitor=.false.
            sgg%Med(contamedia)%lumped(1)%diodo    =.true.
            sgg%Med(contamedia)%lumped(1)%R = this%DielRegs%lins(i)%R
            sgg%Med(contamedia)%lumped(1)%Rtime_on = 0.0 !irrelevant
            sgg%Med(contamedia)%lumped(1)%Rtime_off = 0.0 !irrelevant
            sgg%Med(contamedia)%lumped(1)%L = 0.0_RKIND
            sgg%Med(contamedia)%lumped(1)%C = 0.0_RKIND
            sgg%Med(contamedia)%lumped(1)%DiodB = this%DielRegs%lins(i)%DiodB
            sgg%Med(contamedia)%lumped(1)%DiodIsat = this%DielRegs%lins(i)%DiodIsat
            sgg%Med(contamedia)%lumped(1)%orient = this%DielRegs%lins(i)%DiodOri
         else
            sgg%Med(contamedia)%Is%Lumped = .false.
            if (.not. this%DielRegs%lins(i)%plain) then
               WRITE (buff, '(a)')    'Buggy error 1 in preprocess lumped. .'
               CALL STOPONERROR(layoutnumber,size,buff)
            endif
         endif
!!!fin lumped
         tama2 = (this%DielRegs%lins(i)%n_c2P)
         DO j = 1, tama2
            punto%XI = this%DielRegs%lins(i)%c2P(j)%XI
            punto%XE = this%DielRegs%lins(i)%c2P(j)%XE
            punto%YI = this%DielRegs%lins(i)%c2P(j)%YI
            punto%YE = this%DielRegs%lins(i)%c2P(j)%YE
            punto%ZI = this%DielRegs%lins(i)%c2P(j)%ZI
            punto%ZE = this%DielRegs%lins(i)%c2P(j)%ZE
            orientacion = this%DielRegs%lins(i)%c2P(j)%or
            isathinwire = .FALSE.
            numertag = searchtag(tagtype,this%DielRegs%lins(i)%c2P(j)%tag)
            CALL CreateLineMM (layoutnumber, media%sggMtag, tag_numbers, numertag, media%sggMiEx, media%sggMiEy, media%sggMiEz, &
            & media%sggMiHx, media%sggMiHy, media%sggMiHz, Alloc_iEx_XI, &
            & Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI, &
            & Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI, &
            & Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI, &
            & Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI, &
            & Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg%Med, sgg%NumMedia, sgg%EShared, BoundingBox, punto, orientacion, &
            & contamedia, isathinwire,verbose,numeroasignaciones)
         END DO
         tama3 = (this%DielRegs%lins(i)%n_c1P)
         DO j = 1, tama3
            punto%XI = this%DielRegs%lins(i)%c1P(j)%XI
            punto%XE = this%DielRegs%lins(i)%c1P(j)%XI
            punto%YI = this%DielRegs%lins(i)%c1P(j)%YI
            punto%YE = this%DielRegs%lins(i)%c1P(j)%YI
            punto%ZI = this%DielRegs%lins(i)%c1P(j)%ZI
            punto%ZE = this%DielRegs%lins(i)%c1P(j)%ZI
            orientacion = this%DielRegs%lins(i)%c1P(j)%or
            isathinwire = .FALSE.
            numertag = searchtag(tagtype,this%DielRegs%lins(i)%c1P(j)%tag)
            CALL CreateLineMM (layoutnumber, media%sggMtag, tag_numbers, numertag, media%sggMiEx, media%sggMiEy, media%sggMiEz, &
            & media%sggMiHx, media%sggMiHy, media%sggMiHz, Alloc_iEx_XI, &
            & Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI, &
            & Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI, &
            & Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI, &
            & Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI, &
            & Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg%Med, sgg%NumMedia, sgg%EShared, BoundingBox, punto, orientacion, &
            & contamedia, isathinwire,verbose,numeroasignaciones)
         END DO
      END DO


      !Anisotropic materials
      !materialList
      !BODYes
      tama = (this%ANIMATS%nvols)
      DO i = 1, tama
         contamedia = contamedia + 1
         ALLOCATE (sgg%Med(contamedia)%Anisotropic(1))
         sgg%Med(contamedia)%Is%Anisotropic = .TRUE.
         sgg%Med(contamedia)%Priority = prior_AB
         sgg%Med(contamedia)%Anisotropic(1)%Epr = this%ANIMATS%vols(i)%eps / Eps0
         sgg%Med(contamedia)%Anisotropic(1)%Sigma = this%ANIMATS%vols(i)%Sigma
         sgg%Med(contamedia)%Anisotropic(1)%Mur = this%ANIMATS%vols(i)%mu / Mu0
         sgg%Med(contamedia)%Anisotropic(1)%SigmaM = this%ANIMATS%vols(i)%SigmaM
         tama2 = (this%ANIMATS%vols(i)%n_c2P)
         DO j = 1, tama2
            punto%XI = this%ANIMATS%vols(i)%c2P(j)%XI
            punto%XE = this%ANIMATS%vols(i)%c2P(j)%XE
            punto%YI = this%ANIMATS%vols(i)%c2P(j)%YI
            punto%YE = this%ANIMATS%vols(i)%c2P(j)%YE
            punto%ZI = this%ANIMATS%vols(i)%c2P(j)%ZI
            punto%ZE = this%ANIMATS%vols(i)%c2P(j)%ZE
            numertag = searchtag(tagtype,this%ANIMATS%vols(i)%c2P(j)%tag)
            CALL CreateVolumeMM (layoutnumber, media%sggMtag, tag_numbers, numertag, media%sggMiEx, media%sggMiEy, media%sggMiEz, &
            & media%sggMiHx, media%sggMiHy, media%sggMiHz, Alloc_iEx_XI, &
            & Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI, &
            & Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI, &
            & Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI, &
            & Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI, &
            & Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg%Med, sgg%NumMedia, sgg%EShared, BoundingBox, punto, contamedia)
         END DO
         tama3 = (this%ANIMATS%vols(i)%n_c1P)
         DO j = 1, tama3
            punto%XI = this%ANIMATS%vols(i)%c1P(j)%XI
            punto%XE = this%ANIMATS%vols(i)%c1P(j)%XI
            punto%YI = this%ANIMATS%vols(i)%c1P(j)%YI
            punto%YE = this%ANIMATS%vols(i)%c1P(j)%YI
            punto%ZI = this%ANIMATS%vols(i)%c1P(j)%ZI
            punto%ZE = this%ANIMATS%vols(i)%c1P(j)%ZI
            !
            numertag = searchtag(tagtype,this%ANIMATS%vols(i)%c1P(j)%tag)
            CALL CreateVolumeMM (layoutnumber, media%sggMtag, tag_numbers, numertag, media%sggMiEx, media%sggMiEy, media%sggMiEz, &
            & media%sggMiHx, media%sggMiHy, media%sggMiHz, Alloc_iEx_XI, &
            & Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI, &
            & Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI, &
            & Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI, &
            & Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI, &
            & Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg%Med, sgg%NumMedia, sgg%EShared, BoundingBox, punto, contamedia)
         END DO
      END DO
      !SURFs
      tama = (this%ANIMATS%nsurfs)
      DO i = 1, tama
         contamedia = contamedia + 1
         ALLOCATE (sgg%Med(contamedia)%Anisotropic(1))
         sgg%Med(contamedia)%Is%Anisotropic = .TRUE.
         sgg%Med(contamedia)%Priority = prior_IS
         sgg%Med(contamedia)%Anisotropic(1)%Epr = this%ANIMATS%surfs(i)%eps / Eps0
         sgg%Med(contamedia)%Anisotropic(1)%Sigma = this%ANIMATS%surfs(i)%Sigma
         sgg%Med(contamedia)%Anisotropic(1)%Mur = this%ANIMATS%surfs(i)%mu / Mu0
         sgg%Med(contamedia)%Anisotropic(1)%SigmaM = this%ANIMATS%surfs(i)%SigmaM
         tama2 = (this%ANIMATS%surfs(i)%n_c2P)
         DO j = 1, tama2
            punto%XI = this%ANIMATS%surfs(i)%c2P(j)%XI
            punto%XE = this%ANIMATS%surfs(i)%c2P(j)%XE
            punto%YI = this%ANIMATS%surfs(i)%c2P(j)%YI
            punto%YE = this%ANIMATS%surfs(i)%c2P(j)%YE
            punto%ZI = this%ANIMATS%surfs(i)%c2P(j)%ZI
            punto%ZE = this%ANIMATS%surfs(i)%c2P(j)%ZE
            orientacion = this%ANIMATS%surfs(i)%c2P(j)%or
            numertag = searchtag(tagtype,this%ANIMATS%surfs(i)%c2P(j)%tag)
            CALL CreateSurfaceMM (layoutnumber, media%sggMtag, tag_numbers, numertag, media%sggMiEx, media%sggMiEy, media%sggMiEz, &
            & media%sggMiHx, media%sggMiHy, media%sggMiHz, Alloc_iEx_XI, &
            & Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI, &
            & Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI, &
            & Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI, &
            & Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI, &
            & Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg%Med, sgg%NumMedia, sgg%EShared, BoundingBox, punto, orientacion, &
            & contamedia)
         END DO
         tama3 = (this%ANIMATS%surfs(i)%n_c1P)
         DO j = 1, tama3
            punto%XI = this%ANIMATS%surfs(i)%c1P(j)%XI
            punto%XE = this%ANIMATS%surfs(i)%c1P(j)%XI
            punto%YI = this%ANIMATS%surfs(i)%c1P(j)%YI
            punto%YE = this%ANIMATS%surfs(i)%c1P(j)%YI
            punto%ZI = this%ANIMATS%surfs(i)%c1P(j)%ZI
            punto%ZE = this%ANIMATS%surfs(i)%c1P(j)%ZI
            orientacion = this%ANIMATS%surfs(i)%c1P(j)%or
            numertag = searchtag(tagtype,this%ANIMATS%surfs(i)%c1P(j)%tag)
            CALL CreateSurfaceMM (layoutnumber, media%sggMtag, tag_numbers, numertag, media%sggMiEx, media%sggMiEy, media%sggMiEz, &
            & media%sggMiHx, media%sggMiHy, media%sggMiHz, Alloc_iEx_XI, &
            & Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI, &
            & Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI, &
            & Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI, &
            & Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI, &
            & Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg%Med, sgg%NumMedia, sgg%EShared, BoundingBox, punto, orientacion, &
            & contamedia)
         END DO
      END DO
      !LINs
      tama = (this%ANIMATS%nLINS)
      DO i = 1, tama
         contamedia = contamedia + 1
         ALLOCATE (sgg%Med(contamedia)%Anisotropic(1))
         sgg%Med(contamedia)%Is%Anisotropic = .TRUE.
         sgg%Med(contamedia)%Priority = prior_IL
         sgg%Med(contamedia)%Anisotropic(1)%Epr = this%ANIMATS%lins(i)%eps / Eps0
         sgg%Med(contamedia)%Anisotropic(1)%Sigma = this%ANIMATS%lins(i)%Sigma
         sgg%Med(contamedia)%Anisotropic(1)%Mur = this%ANIMATS%lins(i)%mu / Mu0
         sgg%Med(contamedia)%Anisotropic(1)%SigmaM = this%ANIMATS%lins(i)%SigmaM
         tama2 = (this%ANIMATS%lins(i)%n_c2P)
         DO j = 1, tama2
            punto%XI = this%ANIMATS%lins(i)%c2P(j)%XI
            punto%XE = this%ANIMATS%lins(i)%c2P(j)%XE
            punto%YI = this%ANIMATS%lins(i)%c2P(j)%YI
            punto%YE = this%ANIMATS%lins(i)%c2P(j)%YE
            punto%ZI = this%ANIMATS%lins(i)%c2P(j)%ZI
            punto%ZE = this%ANIMATS%lins(i)%c2P(j)%ZE
            orientacion = this%ANIMATS%lins(i)%c2P(j)%or
            isathinwire = .FALSE.
            numertag = searchtag(tagtype,this%ANIMATS%lins(i)%c2P(j)%tag)
            CALL CreateLineMM (layoutnumber, media%sggMtag, tag_numbers, numertag, media%sggMiEx, media%sggMiEy, media%sggMiEz, &
            & media%sggMiHx, media%sggMiHy, media%sggMiHz, Alloc_iEx_XI, &
            & Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI, &
            & Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI, &
            & Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI, &
            & Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI, &
            & Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg%Med, sgg%NumMedia, sgg%EShared, BoundingBox, punto, orientacion, &
            & contamedia, isathinwire,verbose,numeroasignaciones)
         END DO
         tama3 = (this%ANIMATS%lins(i)%n_c1P)
         DO j = 1, tama3
            punto%XI = this%ANIMATS%lins(i)%c1P(j)%XI
            punto%XE = this%ANIMATS%lins(i)%c1P(j)%XI
            punto%YI = this%ANIMATS%lins(i)%c1P(j)%YI
            punto%YE = this%ANIMATS%lins(i)%c1P(j)%YI
            punto%ZI = this%ANIMATS%lins(i)%c1P(j)%ZI
            punto%ZE = this%ANIMATS%lins(i)%c1P(j)%ZI
            orientacion = this%ANIMATS%lins(i)%c1P(j)%or
            isathinwire = .FALSE.
            numertag = searchtag(tagtype,this%ANIMATS%lins(i)%c1P(j)%tag)
            CALL CreateLineMM (layoutnumber, media%sggMtag, tag_numbers, numertag, media%sggMiEx, media%sggMiEy, media%sggMiEz, &
            & media%sggMiHx, media%sggMiHy, media%sggMiHz, Alloc_iEx_XI, &
            & Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI, &
            & Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI, &
            & Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI, &
            & Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI, &
            & Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg%Med, sgg%NumMedia, sgg%EShared, BoundingBox, punto, orientacion, &
            & contamedia, isathinwire,verbose,numeroasignaciones)
         END DO
      END DO
      !frequency dependent materials
      !bodies
      !
      tama = this%FRQDEPMATS%nvols
      DO i = 1, tama
         contamedia = contamedia + 1
         NULLIFY(fdgeom)
         fdgeom=>this%FRQDEPMATS%vols(i)
         call asignadisper(fdgeom)
         !geometry
         !!!!!!!!!

         tama2 = this%FRQDEPMATS%vols(i)%n_C
         DO j = 1, tama2
            punto%XI = this%FRQDEPMATS%vols(i)%C(j)%XI
            punto%XE = this%FRQDEPMATS%vols(i)%C(j)%XE
            punto%YI = this%FRQDEPMATS%vols(i)%C(j)%YI
            punto%YE = this%FRQDEPMATS%vols(i)%C(j)%YE
            punto%ZI = this%FRQDEPMATS%vols(i)%C(j)%ZI
            punto%ZE = this%FRQDEPMATS%vols(i)%C(j)%ZE
            numertag = searchtag(tagtype,this%FRQDEPMATS%vols(i)%C(j)%tag)
            CALL CreateVolumeMM (layoutnumber, media%sggMtag, tag_numbers, numertag, media%sggMiEx, media%sggMiEy, media%sggMiEz, &
            & media%sggMiHx, media%sggMiHy, media%sggMiHz, Alloc_iEx_XI, &
            & Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI, &
            & Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI, &
            & Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI, &
            & Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI, &
            & Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg%Med, sgg%NumMedia, sgg%EShared, BoundingBox, punto, contamedia)
         END DO
      END DO


      !SURFs
      tama = this%FRQDEPMATS%nsurfs
      DO i = 1, tama
         contamedia = contamedia + 1
         NULLIFY(fdgeom)
         fdgeom=>this%FRQDEPMATS%surfs(i)
         call asignadisper(fdgeom)

         tama2 = this%FRQDEPMATS%surfs(i)%n_C
         DO j = 1, tama2
            punto%XI = this%FRQDEPMATS%surfs(i)%C(j)%XI
            punto%XE = this%FRQDEPMATS%surfs(i)%C(j)%XE
            punto%YI = this%FRQDEPMATS%surfs(i)%C(j)%YI
            punto%YE = this%FRQDEPMATS%surfs(i)%C(j)%YE
            punto%ZI = this%FRQDEPMATS%surfs(i)%C(j)%ZI
            punto%ZE = this%FRQDEPMATS%surfs(i)%C(j)%ZE
            orientacion = this%FRQDEPMATS%surfs(i)%C(j)%or
            numertag = searchtag(tagtype,this%FRQDEPMATS%surfs(i)%C(j)%tag)
            CALL CreateSurfaceMM (layoutnumber, media%sggMtag, tag_numbers, numertag, media%sggMiEx, media%sggMiEy, media%sggMiEz, &
            & media%sggMiHx, media%sggMiHy, media%sggMiHz, Alloc_iEx_XI, &
            & Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI, &
            & Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI, &
            & Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI, &
            & Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI, &
            & Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg%Med, sgg%NumMedia, sgg%EShared, BoundingBox, punto, orientacion, &
            & contamedia)
         END DO
      END DO
      !LINs
      tama = this%FRQDEPMATS%nLINS
      DO i = 1, tama
         contamedia = contamedia + 1
         NULLIFY(fdgeom)
         fdgeom=>this%FRQDEPMATS%lins(i)
         call asignadisper(fdgeom)

         tama2 = this%FRQDEPMATS%lins(i)%n_C
         DO j = 1, tama2
            punto%XI = this%FRQDEPMATS%lins(i)%C(j)%XI
            punto%XE = this%FRQDEPMATS%lins(i)%C(j)%XE
            punto%YI = this%FRQDEPMATS%lins(i)%C(j)%YI
            punto%YE = this%FRQDEPMATS%lins(i)%C(j)%YE
            punto%ZI = this%FRQDEPMATS%lins(i)%C(j)%ZI
            punto%ZE = this%FRQDEPMATS%lins(i)%C(j)%ZE
            orientacion = this%FRQDEPMATS%lins(i)%C(j)%or
            isathinwire = .FALSE.
            numertag = searchtag(tagtype,this%FRQDEPMATS%lins(i)%C(j)%tag)
            CALL CreateLineMM (layoutnumber, media%sggMtag, tag_numbers, numertag, media%sggMiEx, media%sggMiEy, media%sggMiEz, &
            & media%sggMiHx, media%sggMiHy, media%sggMiHz, Alloc_iEx_XI, &
            & Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI, &
            & Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI, &
            & Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI, &
            & Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI, &
            & Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg%Med, sgg%NumMedia, sgg%EShared, BoundingBox, punto, orientacion, &
            & contamedia, isathinwire,verbose,numeroasignaciones)
         END DO
      END DO
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !ISOTROPIC Multiports
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      inicontamedia = contamedia + 1
      maxcontamedia = contamedia
      tama = this%LossyThinSurfs%length
      DO j = 1, tama
         !carbon Multiports a guevo
         if (this%LossyThinSurfs%cs(j)%numcapas==0) then
            this%LossyThinSurfs%cs(j)%numcapas=1
            allocate(this%LossyThinSurfs%cs(j)%SigmaM(1))
            allocate(this%LossyThinSurfs%cs(j)%Sigma(1))
            allocate(this%LossyThinSurfs%cs(j)%EPS(1))
            allocate(this%LossyThinSurfs%cs(j)%MU(1))
            allocate(this%LossyThinSurfs%cs(j)%thk(1))

            !_for_devia 090519
            allocate(this%LossyThinSurfs%cs(j)%SigmaM_devia(1))
            allocate(this%LossyThinSurfs%cs(j)%Sigma_devia(1))
            allocate(this%LossyThinSurfs%cs(j)%EPS_devia(1))
            allocate(this%LossyThinSurfs%cs(j)%MU_devia(1))
            allocate(this%LossyThinSurfs%cs(j)%thk_devia(1))
            !!!

            this%LossyThinSurfs%cs(J)%SigmaM(1)=0.0_RKIND !TRUCO PARA QUE CUANDO NO TENGA CAPAS (LECTURA DESDE FICHERO DE POLOS /RESIDUOS) NO PETE
            this%LossyThinSurfs%cs(J)%Sigma(1)=0.0_RKIND
            this%LossyThinSurfs%cs(J)%EPS(1)=EPS0
            this%LossyThinSurfs%cs(J)%MU(1)=MU0
            this%LossyThinSurfs%cs(J)%thk(1)=-1.0

            !_for_devia 090519
            this%LossyThinSurfs%cs(J)%SigmaM_devia(1)=0.0_RKIND
            this%LossyThinSurfs%cs(J)%Sigma_devia(1)=0.0_RKIND
            this%LossyThinSurfs%cs(J)%EPS_devia(1)=0.0_RKIND
            this%LossyThinSurfs%cs(J)%MU_devia(1)=0.0_RKIND
            this%LossyThinSurfs%cs(J)%thk_devia(1)=0.0_RKIND
            !!!
            !!!comentado el 120219 pq no se lleva bien con Semba !no entiendo ahora el comentario de malonyedispersive!!!120219
            !!!     WRITE (buff, '(a)')    'pre1_Error:  SGBC materials must have at least one layyer even in dummy for malonyedispersive'
            !!!     CALL WarnErrReport (buff,.true.)
         ENDIF


         IF (abs(this%LossyThinSurfs%cs(j)%SigmaM(1)) <= 1.0e-2_RKIND ) THEN  !!!ojoooo a 210319 manda guevos que tengamos que estar con el flag de la conductidad magnetica para llamar a SGBC todavia en 2015!!!
            this%LossyThinSurfs%cs(j)%SigmaM = 0.0_RKIND
            if (.not.mibc) then
               !if (this%LossyThinSurfs%cs(j)%numcapas >1) then
               !   WRITE (buff, '(a)')    'pre1_Warning:  SGBC materials are just averaged for multilayered structures.'// &
               !   ' Use preferably -mibc instead.'
               !   CALL WarnErrReport (buff)
               !endif
               SGBC=.true. !si la conductividad es 0.0_RKIND (o casi) utiliza directamente SGBC
               mibc=.false.
            endif
         endif
         IF (this%LossyThinSurfs%cs(j)%SigmaM(1) >= 0.0_RKIND) THEN
            !SURFs (siempre son surfs)
            !
            tama2 = this%LossyThinSurfs%cs(j)%nc
            mincontamedia = maxcontamedia + 1
            MultiportFile = trim (adjustl(this%LossyThinSurfs%cs(j)%files)) // '_z11.txt'
            DO i = 1, tama2
               orientacion = this%LossyThinSurfs%cs(j)%C(i)%or
               punto%XI = this%LossyThinSurfs%cs(j)%C(i)%XI
               punto%XE = this%LossyThinSurfs%cs(j)%C(i)%XE
               punto%YI = this%LossyThinSurfs%cs(j)%C(i)%YI
               punto%YE = this%LossyThinSurfs%cs(j)%C(i)%YE
               punto%ZI = this%LossyThinSurfs%cs(j)%C(i)%ZI
               punto%ZE = this%LossyThinSurfs%cs(j)%C(i)%ZE
               existia = .FALSE.
               doexis: DO k = inicontamedia, maxcontamedia
                  IF (trim(adjustl(sgg%Med(k)%multiport(1)%multiportFileZ11)) == trim(adjustl(MultiportFile))) THEN
                     IF (sgg%Med(k)%multiport(1)%Multiportdir == orientacion) THEN
                        contamedia = k
                        existia = .TRUE.
                        EXIT doexis
                     END IF
                  END IF
               END DO doexis

               IF ( .NOT. existia) THEN
                  maxcontamedia = maxcontamedia + 1
                  contamedia = maxcontamedia
                  ALLOCATE (sgg%Med(contamedia)%multiport(1))
                  !
                  if ((this%LossyThinSurfs%cs(j)%numcapas >1).and.SGBCDispersive) then
                     WRITE (buff, *)    'ERROR in SGBCs Number of layers >1 still unsupported for SGBCDispersive. '
                     CALL StopOnError (0,0,buff)
                  endif
                  !
                  allocate(sgg%Med(contamedia)%Multiport(1)%epr(1:this%LossyThinSurfs%cs(j)%numcapas), &
                     sgg%Med(contamedia)%Multiport(1)%mur(1:this%LossyThinSurfs%cs(j)%numcapas), &
                     sgg%Med(contamedia)%Multiport(1)%sigma(1:this%LossyThinSurfs%cs(j)%numcapas), &
                     sgg%Med(contamedia)%Multiport(1)%sigmam(1:this%LossyThinSurfs%cs(j)%numcapas), &
                     sgg%Med(contamedia)%Multiport(1)%width(1:this%LossyThinSurfs%cs(j)%numcapas))
                  !_for_devia 090519
                  allocate(sgg%Med(contamedia)%Multiport(1)%epr_devia(1:this%LossyThinSurfs%cs(j)%numcapas), &
                     sgg%Med(contamedia)%Multiport(1)%mur_devia(1:this%LossyThinSurfs%cs(j)%numcapas), &
                     sgg%Med(contamedia)%Multiport(1)%sigma_devia(1:this%LossyThinSurfs%cs(j)%numcapas), &
                     sgg%Med(contamedia)%Multiport(1)%sigmaM_devia(1:this%LossyThinSurfs%cs(j)%numcapas), &
                     sgg%Med(contamedia)%Multiport(1)%width_devia(1:this%LossyThinSurfs%cs(j)%numcapas))
                  !!!
                  puntoXI = Max (punto%XI, Min(BoundingBox%XI, BoundingBox%XE))
                  puntoYI = Max (punto%YI, Min(BoundingBox%YI, BoundingBox%YE))
                  puntoZI = Max (punto%ZI, Min(BoundingBox%ZI, BoundingBox%ZE))

                  !!!!!!!!estaba antes  maaaal. bug 140815verano
                  if(.not.((puntoXI>=sgg%allocDxI).and.(puntoXI<=sgg%allocDxE))) then
                     puntoXI= sgg%allocDxI
                     WRITE (buff, *)    'ERROR: precompo 2: Readjusting composite init point. Only ignore if parts of the geometry fall out of the the domain deliberately (only if manual clipping)',puntoXI,puntoYI,puntoZI,sgg%allocDxI,sgg%allocDyI,sgg%allocDzI
                     CALL WarnErrReport (buff,.TRUE.)
                  endif
                  if(.not.((puntoYI>=sgg%allocDyI).and.(puntoYI<=sgg%allocDyE))) then
                     puntoYI= sgg%allocDyI
                     WRITE (buff, *)    'ERROR: precompo 2: Readjusting composite init point. Only ignore if parts of the geometry fall out of the the domain deliberately (only if manual clipping)',puntoXI,puntoYI,puntoZI,sgg%allocDxI,sgg%allocDyI,sgg%allocDzI
                     CALL WarnErrReport (buff,.TRUE.)
                  endif
                  if(.not.((puntoZI>=sgg%allocDzI).and.(puntoZI<=sgg%allocDzE))) then
                     puntoZI= sgg%allocDzI
                     WRITE (buff, *)    'ERROR: precompo 2: Readjusting composite init point. Only ignore if parts of the geometry fall out of the the domain deliberately (only if manual clipping)',puntoXI,puntoYI,puntoZI,sgg%allocDxI,sgg%allocDyI,sgg%allocDzI
                     CALL WarnErrReport (buff,.TRUE.)
                  endif
                  dentro = (puntoXI>=sgg%allocDxI).and.(puntoXI<=sgg%allocDxE).and. &
                     (puntoYI>=sgg%allocDyI).and.(puntoYI<=sgg%allocDyE).and. &
                     (puntoZI>=sgg%allocDzI).and.(puntoZI<=sgg%allocDzE)
                  delta=-1.0_RKIND
                  IF (DENTRO) THEN
                     select case (abs(this%LossyThinSurfs%cs(j)%C(i)%or))
                      case (iEx)
                        delta=(sgg%DX(puntoXI)+sgg%DX(puntoXI-1))/2.0_RKIND
                      case (iEy)
                        delta=(sgg%DY(puntoYI)+sgg%Dy(puntoYI-1))/2.0_RKIND
                      case (iEz)
                        delta=(sgg%DZ(puntoZI)+sgg%Dz(puntoZI-1))/2.0_RKIND
                      case default
                        WRITE (buff, '(a)')    'Buggy error 1 in preprocess composites. .'
                        CALL STOPONERROR(layoutnumber,size,buff)
                     end select
                  ELSE
                     WRITE (buff, '(a)')    'Buggy error 2 in preprocess composites. .'
                     CALL STOPONERROR(layoutnumber,size,buff)
                  ENDIF
                  sgg%Med(contamedia)%Multiport(1)%numcapas = this%LossyThinSurfs%cs(j)%numcapas
                  !el especificado
                  sgg%Med(contamedia)%Multiport(1)%Multiportdir = this%LossyThinSurfs%cs(j)%C(i)%or
                  do I_=1,sgg%Med(contamedia)%Multiport(1)%numcapas
                     if (sgg%Med(contamedia)%Multiport(1)%Multiportdir>0) then
                        j_=i_
                     else
                        j_=sgg%Med(contamedia)%Multiport(1)%numcapas-i_+1 !dale la vuelta (medios no simetricos) !0121
                     endif
                     sgg%Med(contamedia)%Multiport(1)%epr         (j_) =  this%LossyThinSurfs%cs(j)%eps               (i_)    / Eps0
                     sgg%Med(contamedia)%Multiport(1)%mur         (j_) =  this%LossyThinSurfs%cs(j)%mu                (i_)    / mu0
                     sgg%Med(contamedia)%Multiport(1)%sigma       (j_) =  this%LossyThinSurfs%cs(j)%Sigma             (i_)
                     sgg%Med(contamedia)%Multiport(1)%sigmam      (j_) =  abs(this%LossyThinSurfs%cs(j)%Sigmam        (i_)    )
                     sgg%Med(contamedia)%Multiport(1)%width       (j_) =  this%LossyThinSurfs%cs(j)%thk               (i_)

                     !_for_devia 090519
                     sgg%Med(contamedia)%Multiport(1)%epr_devia   (j_) =  this%LossyThinSurfs%cs(j)%eps_devia         (i_)    / Eps0
                     sgg%Med(contamedia)%Multiport(1)%mur_devia   (j_) =  this%LossyThinSurfs%cs(j)%MU_devia          (i_)    / mu0
                     sgg%Med(contamedia)%Multiport(1)%sigma_devia (j_) =  this%LossyThinSurfs%cs(j)%Sigma_devia       (i_)
                     sgg%Med(contamedia)%Multiport(1)%sigmaM_devia(j_) =  abs(this%LossyThinSurfs%cs(j)%SigmaM_devia  (i_)    )
                     sgg%Med(contamedia)%Multiport(1)%width_devia (j_) =  this%LossyThinSurfs%cs(j)%thk_devia         (i_)
                  end do

                  rdummy=maxval(abs(this%LossyThinSurfs%cs(j)%MU_devia))+maxval(abs(this%LossyThinSurfs%cs(j)%SigmaM_devia))
                  if (rdummy>1.0e-15_RKIND) then
                     WRITE (buff, '(a)')    'Non null deviations found in sigmam or mu in composites. Still unsupported.'
                     CALL STOPONERROR(layoutnumber,size,buff)
                  endif

                  !!!


!!old pre 17/07/15
!!!                           sgg%Med(contamedia)%Multiport(1)%transversalSpaceDelta=delta
                  sgg%Med(contamedia)%Priority = prior_CS
                  sgg%Med(contamedia)%Epr   = this%LossyThinSurfs%cs(j)%eps(1) / Eps0
                  sgg%Med(contamedia)%Sigma =this%LossyThinSurfs%cs(j)%Sigma(1)
                  sgg%Med(contamedia)%Mur = this%LossyThinSurfs%cs(j)%mu(1) / Mu0
                  sgg%Med(contamedia)%SigmaM = abs(this%LossyThinSurfs%cs(j)%SigmaM(1)) !may be negative
                  if (mibc) then
                     sgg%Med(contamedia)%Is%multiport = .TRUE.
                     sgg%Med(contamedia)%Is%Lossy = .true.
                  elseif (SGBC) then
                     sgg%Med(contamedia)%Is%SGBC = .TRUE.
                     sgg%Med(contamedia)%Is%Lossy = .true.
                     if (SGBCDispersive)   sgg%Med(contamedia)%Is%SGBCDispersive = .TRUE.
                  else
                     WRITE (buff, '(a)')    'Some -mibc -sgbc switch should be used for Composites.'
                     CALL STOPONERROR(layoutnumber,size,buff)
                  endif
                  sgg%Med(contamedia)%Is%Dielectric = .FALSE.

                  sgg%Med(contamedia)%multiport(1)%multiportFileZ11 =  trim &
                  & (adjustl(this%LossyThinSurfs%cs(j)%files)) // '_z11.txt'
                  sgg%Med(contamedia)%multiport(1)%multiportFileZ22 =  trim &
                  & (adjustl(this%LossyThinSurfs%cs(j)%files)) // '_z22.txt'
                  sgg%Med(contamedia)%multiport(1)%multiportFileZ12 =  trim &
                  & (adjustl(this%LossyThinSurfs%cs(j)%files)) // '_z12.txt'
                  sgg%Med(contamedia)%multiport(1)%multiportFileZ21 =  trim &
                  & (adjustl(this%LossyThinSurfs%cs(j)%files)) // '_z12.txt'
                  !
!!
                  if (mibc) then     !se trataran con MIBC !!!151161
                     sgg%Med(contamedia)%Is%SGBC = .false.
                     sgg%Med(contamedia)%Is%SGBCDispersive = .false.
                     sgg%Med(contamedia)%Is%Lossy = .true.

                     errnofile1=.false.
                     errnofile2=.false.
                     errnofile3=.false.
                     errnofile4=.false.
                     pozi=index(sgg%Med(contamedia)%Multiport(1)%multiportFileZ11,'_z11.txt')
                     multiportfile2=trim(adjustl(sgg%Med(contamedia)%Multiport(1)%multiportFileZ11(1:pozi-1)))
                     INQUIRE (file=trim(adjustl(multiportfile2)), EXIST=errnofile1)
                     INQUIRE (file=trim(adjustl(sgg%Med(contamedia)%Multiport(1)%multiportFileZ11)), EXIST=errnofile2)
                     INQUIRE (file=trim(adjustl(sgg%Med(contamedia)%Multiport(1)%multiportFileZ12)), EXIST=errnofile3)
                     INQUIRE (file=trim(adjustl(sgg%Med(contamedia)%Multiport(1)%multiportFileZ22)), EXIST=errnofile4)

                     if (.not.(errnofile1.or.(errnofile2.and.errnofile3.and.errnofile4))) then
                        buff='Neither New nor Old style mibc FILE '//trim(adjustl(multiportfile2))//' EXISTS.'
                        CALL WarnErrReport (buff,.TRUE.)
                     endif
                  endif

                  !!!!!!!!end 09/07/13
                  !
                  !
               END IF
               !
               !
               numertag = searchtag(tagtype,this%LossyThinSurfs%cs(j)%C(i)%tag)
               CALL CreateSurfaceMM (layoutnumber, media%sggMtag, tag_numbers, numertag, media%sggMiEx, media%sggMiEy, media%sggMiEz, &
               & media%sggMiHx, media%sggMiHy, media%sggMiHz, Alloc_iEx_XI, &
               & Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI, &
               & Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI, &
               & Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI, &
               & Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI, &
               & Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg%Med, sgg%NumMedia, sgg%EShared, BoundingBox, punto, orientacion, &
               & contamedia)
            END DO
         END IF
         !Multiports
      END DO
      contamedia = maxcontamedia
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !end ISOTROPIC multiports
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !ANISOTROPIC Multiports
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      inicontamedia = contamedia + 1
      maxcontamedia = contamedia
      tama = this%LossyThinSurfs%length
      DO j = 1, tama
         !carbon Multiports a guevo
         IF (this%LossyThinSurfs%cs(j)%SigmaM(1) < 0.0_RKIND) THEN
            !SURFs (siempre son surfs)
            tama2 = this%LossyThinSurfs%cs(j)%nc
            mincontamedia = maxcontamedia + 1
            MultiportFile =  trim (adjustl(this%LossyThinSurfs%cs(j)%files)) // '_z11.txt'
            DO i = 1, tama2
               orientacion = this%LossyThinSurfs%cs(j)%C(i)%or
               punto%XI = this%LossyThinSurfs%cs(j)%C(i)%XI
               punto%XE = this%LossyThinSurfs%cs(j)%C(i)%XE
               punto%YI = this%LossyThinSurfs%cs(j)%C(i)%YI
               punto%YE = this%LossyThinSurfs%cs(j)%C(i)%YE
               punto%ZI = this%LossyThinSurfs%cs(j)%C(i)%ZI
               punto%ZE = this%LossyThinSurfs%cs(j)%C(i)%ZE
               existia = .FALSE.
               doexis2: DO k = inicontamedia, maxcontamedia
                  IF (trim(adjustl(sgg%Med(k)%AnisMultiport(1)%multiportFileZ11)) == trim(adjustl(MultiportFile))) THEN
                     IF (sgg%Med(k)%AnisMultiport(1)%Multiportdir == orientacion) THEN
                        contamedia = k
                        existia = .TRUE.
                        EXIT doexis2
                     END IF
                  END IF
               END DO doexis2
               IF ( .NOT. existia) THEN
                  maxcontamedia = maxcontamedia + 1
                  contamedia = maxcontamedia
                  ALLOCATE (sgg%Med(contamedia)%AnisMultiport(1))
                  !
                  !

                  if (this%LossyThinSurfs%cs(j)%numcapas >1) then
                     WRITE (buff, '(a)')    'pre1_ERROR:  Anisotropic multiport materials unsupported for multilayered structures.'
                     CALL WarnErrReport (buff,.TRUE.)
                  endif
                  puntoXI = Max (punto%XI, Min(BoundingBox%XI, BoundingBox%XE)) !copiado de healer
                  puntoYI = Max (punto%YI, Min(BoundingBox%YI, BoundingBox%YE))
                  puntoZI = Max (punto%ZI, Min(BoundingBox%ZI, BoundingBox%ZE))


                  !!!!!!!!estaba antes  maaaal. bug 140815verano
                  if(.not.((puntoXI>=sgg%allocDxI).and.(puntoXI<=sgg%allocDxE))) puntoXI= sgg%allocDxI
                  if(.not.((puntoYI>=sgg%allocDyI).and.(puntoYI<=sgg%allocDyE))) puntoYI= sgg%allocDyI
                  if(.not.((puntoZI>=sgg%allocDzI).and.(puntoZI<=sgg%allocDzE))) puntoZI= sgg%allocDzI
                  WRITE (buff, '(a)')    'ERROR: precompo 2: Readjusting composite init point. Only ignore if parts of the geometry fall out of the the domain deliberately (only if manual clipping)'
                  CALL WarnErrReport (buff,.TRUE.)

                  dentro = (puntoXI>=sgg%allocDxI).and.(puntoXI<=sgg%allocDxE).and. &
                     (puntoYI>=sgg%allocDyI).and.(puntoYI<=sgg%allocDyE).and. &
                     (puntoZI>=sgg%allocDzI).and.(puntoZI<=sgg%allocDzE)
                  delta=-1.0_RKIND
                  IF (DENTRO) THEN
                     select case (abs(this%LossyThinSurfs%cs(j)%C(i)%or))
                      case (iEx)
                        delta=(sgg%DX(puntoXI)+sgg%DX(puntoXI-1))/2.0_RKIND
                      case (iEy)
                        delta=(sgg%DY(puntoYI)+sgg%Dy(puntoYI-1))/2.0_RKIND
                      case (iEz)
                        delta=(sgg%DZ(puntoZI)+sgg%Dz(puntoZI-1))/2.0_RKIND
                      case default
                        WRITE (buff, '(a)')    'Buggy error 1 in preprocess composites. .'
                        CALL STOPONERROR(layoutnumber,size,buff)
                     end select
                  ELSE
                     WRITE (buff, '(a)')    'Buggy error 2 in preprocess composites. .'
                     CALL STOPONERROR(layoutnumber,size,buff)
                  ENDIF
                  sgg%Med(contamedia)%AnisMultiport(1)%Multiportdir = this%LossyThinSurfs%cs(j)%C(i)%or
                  sgg%Med(contamedia)%AnisMultiport(1)%epr    = this%LossyThinSurfs%cs(j)%eps / Eps0
                  sgg%Med(contamedia)%AnisMultiport(1)%mur    =  this%LossyThinSurfs%cs(j)%mu / mu0
                  sgg%Med(contamedia)%AnisMultiport(1)%sigma  =  this%LossyThinSurfs%cs(j)%Sigma
                  sgg%Med(contamedia)%AnisMultiport(1)%sigmam =  abs(this%LossyThinSurfs%cs(j)%Sigmam)
                  sgg%Med(contamedia)%AnisMultiport(1)%width  =  this%LossyThinSurfs%cs(j)%thk
                  sgg%Med(contamedia)%Priority = prior_CS
                  sgg%Med(contamedia)%Epr    =this%LossyThinSurfs%cs(j)%eps(1) / Eps0
                  sgg%Med(contamedia)%Sigma =this%LossyThinSurfs%cs(j)%Sigma(1)
                  sgg%Med(contamedia)%Mur = this%LossyThinSurfs%cs(j)%mu(1) / Mu0
                  sgg%Med(contamedia)%SigmaM = abs(this%LossyThinSurfs%cs(j)%SigmaM(1)) !may be negative
                  !
                  if (mibc) then
                     sgg%Med(contamedia)%Is%Anismultiport = .TRUE.
                     sgg%Med(contamedia)%Is%Lossy = .TRUE.
                  else
                     WRITE (buff, '(a)')    'Some -mibc -sgbc switch should be used for Anisotropic Composites.'
                     CALL STOPONERROR(layoutnumber,size,buff)
                  endif

                  sgg%Med(contamedia)%Is%Dielectric = .FALSE.

                  sgg%Med(contamedia)%AnisMultiport(1)%multiportFileZ11 =  trim &
                  & (adjustl(this%LossyThinSurfs%cs(j)%files)) // '_z11.txt'
                  sgg%Med(contamedia)%AnisMultiport(1)%multiportFileZ22 =  trim &
                  & (adjustl(this%LossyThinSurfs%cs(j)%files)) // '_z22.txt'
                  sgg%Med(contamedia)%AnisMultiport(1)%multiportFileZ12 =  trim &
                  & (adjustl(this%LossyThinSurfs%cs(j)%files)) // '_z12.txt'
                  sgg%Med(contamedia)%AnisMultiport(1)%multiportFileZ21 =  trim &
                  & (adjustl(this%LossyThinSurfs%cs(j)%files)) // '_z12.txt'
                  !
               END IF
               !
               !
               !
               !
               numertag = searchtag(tagtype,this%LossyThinSurfs%cs(j)%C(i)%tag)
               CALL CreateSurfaceMM (layoutnumber, media%sggMtag, tag_numbers, numertag, media%sggMiEx, media%sggMiEy, media%sggMiEz, &
               & media%sggMiHx, media%sggMiHy, media%sggMiHz, Alloc_iEx_XI, &
               & Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI, &
               & Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI, &
               & Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI, &
               & Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI, &
               & Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg%Med, sgg%NumMedia, sgg%EShared, BoundingBox, punto, orientacion, &
               & contamedia)
            END DO
         END IF
         !Multiports
      END DO
      contamedia = maxcontamedia
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !end ANISOTROPIC multiports
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Multiports lossy padding
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (abs(attfactor-1.0_RKIND) > 1.0e-12_RKIND) then
         tama = this%LossyThinSurfs%length
         DO j = 1, tama
            tama2 = this%LossyThinSurfs%cs(j)%nc
            contamedia = contamedia + 1
            DO i = 1, tama2
               orientacion = this%LossyThinSurfs%cs(j)%C(i)%or
               !ES UN free-space multiportpadding CON LA PRIORIDAD DE UN MULTIPORT con conductividad magnetica que luego se desanulara
               sgg%Med(contamedia)%Priority = prior_CS
               sgg%Med(contamedia)%Is%multiport = .FALSE.
               sgg%Med(contamedia)%Is%ANISmultiport = .FALSE.
               sgg%Med(contamedia)%Is%MultiportPadding = .TRUE.
               sgg%Med(contamedia)%Is%Lossy = .TRUE.
               sgg%Med(contamedia)%Is%Dielectric = .False.
               sgg%Med(contamedia)%Epr = 1.0_RKIND !this%LossyThinSurfs%cs(j)%eps / Eps0
               sgg%Med(contamedia)%Sigma = 0.0_RKIND !abs(this%LossyThinSurfs%cs(j)%Sigma) !may be negative
               sgg%Med(contamedia)%Mur = 1.0_RKIND !this%LossyThinSurfs%cs(j)%mu / Mu0
               sgg%Med(contamedia)%Is%Dielectric = .TRUE.
               !provisionalmente (luego se retocara con el sigmam correcto)
               sgg%Med(contamedia)%SigmaM = 0.0_RKIND !abs(this%LossyThinSurfs%cs(j)%Sigma) !may be negative

               punto%XI = this%LossyThinSurfs%cs(j)%C(i)%XI
               punto%XE = this%LossyThinSurfs%cs(j)%C(i)%XE
               punto%YI = this%LossyThinSurfs%cs(j)%C(i)%YI
               punto%YE = this%LossyThinSurfs%cs(j)%C(i)%YE
               punto%ZI = this%LossyThinSurfs%cs(j)%C(i)%ZI
               punto%ZE = this%LossyThinSurfs%cs(j)%C(i)%ZE
               !!

               SELECT CASE (Abs(orientacion))
                CASE (iEx)
                  punto%XI = this%LossyThinSurfs%cs(j)%C(i)%XI
                  punto%XE = this%LossyThinSurfs%cs(j)%C(i)%XI
                  numertag = searchtag(tagtype,this%LossyThinSurfs%cs(j)%C(i)%tag)
                  CALL CreateMagneticSurface (layoutnumber, media%sggMtag, tag_numbers, numertag, media%sggMiEx, media%sggMiEy,   &
                     media%sggMiEz, &
                  & media%sggMiHx, media%sggMiHy, media%sggMiHz,   &
                     Alloc_iEx_XI, &
                  & Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE,   &
                     Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI, &
                  & Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE,   &
                     Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI, &
                  & Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE,   &
                     Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI, &
                  & Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE,   &
                     Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI, &
                  & Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg%Med, sgg%NumMedia,   &
                     sgg%EShared, BoundingBox, punto, orientacion, &
                  & contamedia)
                  punto%XI = this%LossyThinSurfs%cs(j)%C(i)%XI-1
                  punto%XE = this%LossyThinSurfs%cs(j)%C(i)%XI-1
                  numertag = searchtag(tagtype,this%LossyThinSurfs%cs(j)%C(i)%tag)
                  CALL CreateMagneticSurface (layoutnumber, media%sggMtag, tag_numbers, numertag, media%sggMiEx, media%sggMiEy,  &
                     media%sggMiEz, &
                  & media%sggMiHx, media%sggMiHy, media%sggMiHz,   &
                     Alloc_iEx_XI, &
                  & Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI,   &
                     Alloc_iEy_XE, Alloc_iEy_YI, &
                  & Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI,   &
                     Alloc_iEz_YE, Alloc_iEz_ZI, &
                  & Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI,   &
                     Alloc_iHx_ZE, Alloc_iHy_XI, &
                  & Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI,   &
                     Alloc_iHz_XE, Alloc_iHz_YI, &
                  & Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg%Med, sgg%NumMedia, sgg%EShared,   &
                     BoundingBox, punto, orientacion, &
                  & contamedia)
                CASE (iEy)
                  punto%YI = this%LossyThinSurfs%cs(j)%C(i)%YI
                  punto%YE = this%LossyThinSurfs%cs(j)%C(i)%YI
                  numertag = searchtag(tagtype,this%LossyThinSurfs%cs(j)%C(i)%tag)
                  CALL CreateMagneticSurface (layoutnumber, media%sggMtag, tag_numbers, numertag, media%sggMiEx, media%sggMiEy,   &
                     media%sggMiEz, &
                  & media%sggMiHx, media%sggMiHy, media%sggMiHz,   &
                     Alloc_iEx_XI, &
                  & Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI,   &
                     Alloc_iEy_XE, Alloc_iEy_YI, &
                  & Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI,   &
                     Alloc_iEz_YE, Alloc_iEz_ZI, &
                  & Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI,   &
                     Alloc_iHx_ZE, Alloc_iHy_XI, &
                  & Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI,   &
                     Alloc_iHz_XE, Alloc_iHz_YI, &
                  & Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg%Med, sgg%NumMedia, sgg%EShared,   &
                     BoundingBox, punto, orientacion, &
                  & contamedia)
                  punto%YI = this%LossyThinSurfs%cs(j)%C(i)%YI-1
                  punto%YE = this%LossyThinSurfs%cs(j)%C(i)%YI-1
                  numertag = searchtag(tagtype,this%LossyThinSurfs%cs(j)%C(i)%tag)
                  CALL CreateMagneticSurface (layoutnumber, media%sggMtag, tag_numbers, numertag, media%sggMiEx, media%sggMiEy,   &
                     media%sggMiEz, &
                  & media%sggMiHx, media%sggMiHy, media%sggMiHz,   &
                     Alloc_iEx_XI, &
                  & Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI,   &
                     Alloc_iEy_XE, Alloc_iEy_YI, &
                  & Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI,   &
                     Alloc_iEz_YE, Alloc_iEz_ZI, &
                  & Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI,   &
                     Alloc_iHx_ZE, Alloc_iHy_XI, &
                  & Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI,   &
                     Alloc_iHz_XE, Alloc_iHz_YI, &
                  & Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg%Med, sgg%NumMedia, sgg%EShared,   &
                     BoundingBox, punto, orientacion, &
                  & contamedia)
                CASE (iEz)
                  punto%ZI = this%LossyThinSurfs%cs(j)%C(i)%ZI
                  punto%ZE = this%LossyThinSurfs%cs(j)%C(i)%ZI
                  numertag = searchtag(tagtype,this%LossyThinSurfs%cs(j)%C(i)%tag)
                  CALL CreateMagneticSurface (layoutnumber, media%sggMtag, tag_numbers, numertag, media%sggMiEx, media%sggMiEy,   &
                     media%sggMiEz, &
                  & media%sggMiHx, media%sggMiHy, media%sggMiHz,   &
                     Alloc_iEx_XI, &
                  & Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI,   &
                     Alloc_iEy_XE, Alloc_iEy_YI, &
                  & Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI,   &
                     Alloc_iEz_YE, Alloc_iEz_ZI, &
                  & Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI,   &
                     Alloc_iHx_ZE, Alloc_iHy_XI, &
                  & Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI,   &
                     Alloc_iHz_XE, Alloc_iHz_YI, &
                  & Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg%Med, sgg%NumMedia, sgg%EShared,   &
                     BoundingBox, punto, orientacion, &
                  & contamedia)
                  punto%ZI = this%LossyThinSurfs%cs(j)%C(i)%ZI-1
                  punto%ZE = this%LossyThinSurfs%cs(j)%C(i)%ZI-1
                  numertag = searchtag(tagtype,this%LossyThinSurfs%cs(j)%C(i)%tag)
                  CALL CreateMagneticSurface (layoutnumber, media%sggMtag, tag_numbers, numertag, media%sggMiEx, media%sggMiEy,   &
                     media%sggMiEz, &
                  & media%sggMiHx, media%sggMiHy, media%sggMiHz,   &
                     Alloc_iEx_XI, &
                  & Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI,   &
                     Alloc_iEy_XE, Alloc_iEy_YI, &
                  & Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI,   &
                     Alloc_iEz_YE, Alloc_iEz_ZI, &
                  & Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI,   &
                     Alloc_iHx_ZE, Alloc_iHy_XI, &
                  & Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI,   &
                     Alloc_iHz_XE, Alloc_iHz_YI, &
                  & Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg%Med, sgg%NumMedia, sgg%EShared,   &
                     BoundingBox, punto, orientacion, &
                  & contamedia)
               end select
            end do
         end do
      endif
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !wires
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      oldcontamedia = contamedIA
      tama = this%twires%n_tw
      DO j = 1, tama
         contamedia = contamedia + 1
         ALLOCATE (sgg%Med(contamedia)%wire(1))
         sgg%Med(contamedia)%Priority = prior_TW

         !
         !background
         !
         sgg%Med(contamedia)%Epr = sgg%Med(1)%Epr
         sgg%Med(contamedia)%Sigma = sgg%Med(1)%Sigma
         sgg%Med(contamedia)%Mur = sgg%Med(1)%Mur
         sgg%Med(contamedia)%SigmaM = sgg%Med(1)%SigmaM
         sgg%Med(contamedia)%Is%ThinWire = .TRUE.
         sgg%Med(contamedia)%Is%Dielectric = .FALSE.
         sgg%Med(contamedia)%wire(1)%radius = this%twires%TW(j)%RAD
         sgg%Med(contamedia)%wire(1)%radius_devia = this%twires%TW(j)%RAD_devia
         if (boundwireradius) then
            if (sgg%Med(contamedia)%wire(1)%radius > maxwireradius) sgg%Med(contamedia)%wire(1)%radius=maxwireradius
         endif
         sgg%Med(contamedia)%wire(1)%R = this%twires%TW(j)%RES
         sgg%Med(contamedia)%wire(1)%l = this%twires%TW(j)%IND
         sgg%Med(contamedia)%wire(1)%C = this%twires%TW(j)%CAP
         sgg%Med(contamedia)%wire(1)%P_R = this%twires%TW(j)%P_RES
         sgg%Med(contamedia)%wire(1)%P_l = this%twires%TW(j)%P_IND

         sgg%Med(contamedia)%wire(1)%P_C = this%twires%TW(j)%P_CAP
         if (this%twires%TW(j)%disp) then
            allocate (sgg%Med(contamedia)%wire(1)%disp(1))
            call asignawiredisper(sgg%Med(contamedia)%wire(1)%disp(1), &
               this%twires%TW(j)%dispfile)
         end if
         sgg%Med(contamedia)%wire(1)%LeftEnd = this%twires%TW(j)%LeftEnd
         sgg%Med(contamedia)%wire(1)%RightEnd = this%twires%TW(j)%RightEnd
         sgg%Med(contamedia)%wire(1)%VsourceExists = .FALSE.
         sgg%Med(contamedia)%wire(1)%IsourceExists = .FALSE.
         !
         sgg%Med(contamedia)%wire(1)%HasAbsorbing_RightEnd = .FALSE.
         sgg%Med(contamedia)%wire(1)%HasAbsorbing_LeftEnd = .FALSE.
         sgg%Med(contamedia)%wire(1)%HasParallel_RightEnd = .FALSE.
         sgg%Med(contamedia)%wire(1)%HasParallel_LeftEnd = .FALSE.
         sgg%Med(contamedia)%wire(1)%HasSeries_RightEnd = .FALSE.
         sgg%Med(contamedia)%wire(1)%HasSeries_LeftEnd = .FALSE.
         sgg%Med(contamedia)%wire(1)%Parallel_R_RightEnd = 0.0_RKIND
         sgg%Med(contamedia)%wire(1)%Parallel_R_LeftEnd = 0.0_RKIND
         sgg%Med(contamedia)%wire(1)%Series_R_RightEnd = 0.0_RKIND
         sgg%Med(contamedia)%wire(1)%Series_R_LeftEnd = 0.0_RKIND
         sgg%Med(contamedia)%wire(1)%Parallel_L_RightEnd = 0.0_RKIND
         sgg%Med(contamedia)%wire(1)%Parallel_L_LeftEnd = 0.0_RKIND
         sgg%Med(contamedia)%wire(1)%Series_L_RightEnd = 0.0_RKIND
         sgg%Med(contamedia)%wire(1)%Series_L_LeftEnd = 0.0_RKIND
         sgg%Med(contamedia)%wire(1)%Parallel_C_RightEnd = 0.0_RKIND
         sgg%Med(contamedia)%wire(1)%Parallel_C_LeftEnd = 0.0_RKIND
         sgg%Med(contamedia)%wire(1)%Series_C_RightEnd = 2.0e7_RKIND !en corto 14/2/14
         sgg%Med(contamedia)%wire(1)%Series_C_LeftEnd = 2.0e7_RKIND !en corto 14/2/14
!stoch
         sgg%Med(contamedia)%wire(1)%R_devia = this%twires%TW(j)%RES_devia
         sgg%Med(contamedia)%wire(1)%l_devia = this%twires%TW(j)%IND_devia
         sgg%Med(contamedia)%wire(1)%C_devia = this%twires%TW(j)%CAP_devia

         sgg%Med(contamedia)%wire(1)%Parallel_R_RightEnd_devia = 0.0_RKIND
         sgg%Med(contamedia)%wire(1)%Parallel_R_LeftEnd_devia = 0.0_RKIND
         sgg%Med(contamedia)%wire(1)%Parallel_L_RightEnd_devia = 0.0_RKIND
         sgg%Med(contamedia)%wire(1)%Parallel_L_LeftEnd_devia = 0.0_RKIND
         sgg%Med(contamedia)%wire(1)%Parallel_C_RightEnd_devia = 0.0_RKIND
         sgg%Med(contamedia)%wire(1)%Parallel_C_LeftEnd_devia = 0.0_RKIND
         sgg%Med(contamedia)%wire(1)%Series_R_RightEnd_devia = 0.0_RKIND
         sgg%Med(contamedia)%wire(1)%Series_R_LeftEnd_devia = 0.0_RKIND
         sgg%Med(contamedia)%wire(1)%Series_L_RightEnd_devia = 0.0_RKIND
         sgg%Med(contamedia)%wire(1)%Series_L_LeftEnd_devia = 0.0_RKIND
         sgg%Med(contamedia)%wire(1)%Series_C_RightEnd_devia = 0.0_RKIND
         sgg%Med(contamedia)%wire(1)%Series_C_LeftEnd_devia = 0.0_RKIND
!fin stoch
         !
         IF     (this%twires%TW(j)%TL ==MATERIAL_absorbing) then
            sgg%Med(contamedia)%wire(1)%HasAbsorbing_LeftEnd = .TRUE.
         elseIF (this%twires%TW(j)%TL == Parallel_CONS) THEN
            sgg%Med(contamedia)%wire(1)%HasParallel_LeftEnd = .TRUE.
            sgg%Med(contamedia)%wire(1)%Parallel_R_LeftEnd = this%twires%TW(j)%R_LeftEnd
            sgg%Med(contamedia)%wire(1)%Parallel_L_LeftEnd = this%twires%TW(j)%L_LeftEnd
            sgg%Med(contamedia)%wire(1)%Parallel_C_LeftEnd = this%twires%TW(j)%C_LeftEnd
!
            sgg%Med(contamedia)%wire(1)%Parallel_R_LeftEnd_devia = this%twires%TW(j)%R_LeftEnd_devia
            sgg%Med(contamedia)%wire(1)%Parallel_L_LeftEnd_devia = this%twires%TW(j)%L_LeftEnd_devia
            sgg%Med(contamedia)%wire(1)%Parallel_C_LeftEnd_devia = this%twires%TW(j)%C_LeftEnd_devia

         ELSE IF (this%twires%TW(j)%TL == SERIES_CONS) THEN
            sgg%Med(contamedia)%wire(1)%HasSeries_LeftEnd = .TRUE.
            sgg%Med(contamedia)%wire(1)%Series_R_LeftEnd = this%twires%TW(j)%R_LeftEnd
            sgg%Med(contamedia)%wire(1)%Series_L_LeftEnd = this%twires%TW(j)%L_LeftEnd
            sgg%Med(contamedia)%wire(1)%Series_C_LeftEnd = this%twires%TW(j)%C_LeftEnd
!
            sgg%Med(contamedia)%wire(1)%Series_R_LeftEnd_devia = this%twires%TW(j)%R_LeftEnd_devia
            sgg%Med(contamedia)%wire(1)%Series_L_LeftEnd_devia = this%twires%TW(j)%L_LeftEnd_devia
            sgg%Med(contamedia)%wire(1)%Series_C_LeftEnd_devia = this%twires%TW(j)%C_LeftEnd_devia
         ELSE IF (this%twires%TW(j)%TL == DISPERSIVE_CONS) THEN
            allocate (sgg%Med(contamedia)%wire(1)%disp_LeftEnd(1))
            call asignawiredisper(sgg%Med(contamedia)%wire(1)%disp_LeftEnd(1), &
               this%twires%TW(j)%dispfile_LeftEnd)
         END IF
         !

         IF     (this%twires%TW(j)%TR ==MATERIAL_absorbing) then
            sgg%Med(contamedia)%wire(1)%HasAbsorbing_RightEnd = .TRUE.
         elseIF (this%twires%TW(j)%TR == Parallel_CONS) THEN
            sgg%Med(contamedia)%wire(1)%HasParallel_RightEnd = .TRUE.
            sgg%Med(contamedia)%wire(1)%Parallel_R_RightEnd = this%twires%TW(j)%R_RightEnd
            sgg%Med(contamedia)%wire(1)%Parallel_L_RightEnd = this%twires%TW(j)%L_RightEnd
            sgg%Med(contamedia)%wire(1)%Parallel_C_RightEnd = this%twires%TW(j)%C_RightEnd
!
            sgg%Med(contamedia)%wire(1)%Parallel_R_RightEnd_devia = this%twires%TW(j)%R_RightEnd_devia
            sgg%Med(contamedia)%wire(1)%Parallel_L_RightEnd_devia = this%twires%TW(j)%L_RightEnd_devia
            sgg%Med(contamedia)%wire(1)%Parallel_C_RightEnd_devia = this%twires%TW(j)%C_RightEnd_devia
         ELSE IF (this%twires%TW(j)%TR == SERIES_CONS) THEN
            sgg%Med(contamedia)%wire(1)%HasSeries_RightEnd = .TRUE.
            sgg%Med(contamedia)%wire(1)%Series_R_RightEnd = this%twires%TW(j)%R_RightEnd
            sgg%Med(contamedia)%wire(1)%Series_L_RightEnd = this%twires%TW(j)%L_RightEnd
            sgg%Med(contamedia)%wire(1)%Series_C_RightEnd = this%twires%TW(j)%C_RightEnd
!
            sgg%Med(contamedia)%wire(1)%Series_R_RightEnd_devia = this%twires%TW(j)%R_RightEnd_devia
            sgg%Med(contamedia)%wire(1)%Series_L_RightEnd_devia = this%twires%TW(j)%L_RightEnd_devia
            sgg%Med(contamedia)%wire(1)%Series_C_RightEnd_devia = this%twires%TW(j)%C_RightEnd_devia

         ELSE IF (this%twires%TW(j)%TR == DISPERSIVE_CONS) THEN
            allocate (sgg%Med(contamedia)%wire(1)%disp_RightEnd(1))
            call asignawiredisper(sgg%Med(contamedia)%wire(1)%disp_RightEnd(1), &
               this%twires%TW(j)%dispfile_RightEnd)
         END IF

!stoch
         rdummy=abs(sgg%Med(contamedia)%wire(1)%radius_devia) &
            +abs(sgg%Med(contamedia)%wire(1)%l_devia) &
            +abs(sgg%Med(contamedia)%wire(1)%C_devia) &
            +abs(sgg%Med(contamedia)%wire(1)%Parallel_L_RightEnd_devia) &
            +abs(sgg%Med(contamedia)%wire(1)%Parallel_L_LeftEnd_devia) &
            +abs(sgg%Med(contamedia)%wire(1)%Parallel_C_RightEnd_devia) &
            +abs(sgg%Med(contamedia)%wire(1)%Parallel_C_LeftEnd_devia) &
            +abs(sgg%Med(contamedia)%wire(1)%Series_L_RightEnd_devia)   &
            +abs(sgg%Med(contamedia)%wire(1)%Series_L_LeftEnd_devia)   &
            +abs(sgg%Med(contamedia)%wire(1)%Series_C_RightEnd_devia)   &
            +abs(sgg%Med(contamedia)%wire(1)%Series_C_LeftEnd_devia)
         if (rdummy>1.0e-15_RKIND) then
            WRITE (buff, '(a)')    'Non null deviations found in L, C or radius in wires stoch. Still unsupported.'
            CALL STOPONERROR(layoutnumber,size,buff)
         endif
!fin stoch
         !
         !esto se soportaba desde versiones antiguas (hilos de un solo segmento. Por error se descomento en la R2417 cuando se trabajo en lo del strictnfde tras vuelta de madrid
         !vuelvo a comentarlo porque si que tenemos la capacidad de hilos de un solo segmento
         !
         !        tama2 = this%twires%TW(j)%N_TWC
         !        if (tama2 == 1) then
         !           call stoponerror(layoutnumber,size,'A WIRE must have at least two segments')
         !        endif
         !
         !esto no es ya necesario porque lo calculo yo luego en el wires
         !!record the LeftEnd and RightEnd coordinates (first and last points)
         !!
         !        sgg%Med(contamedia)%wire(1)%LextremoI = this%twires%TW(j)%TWC(1)%i
         !        sgg%Med(contamedia)%wire(1)%LextremoJ = this%twires%TW(j)%TWC(1)%j
         !        sgg%Med(contamedia)%wire(1)%LextremoK = this%twires%TW(j)%TWC(1)%k
         !        orientacionL                          = this%twires%TW(j)%TWC(1)%D
         !        sgg%Med(contamedia)%wire(1)%RextremoI = this%twires%TW(j)%TWC(tama2)%i
         !        sgg%Med(contamedia)%wire(1)%RextremoJ = this%twires%TW(j)%TWC(tama2)%j
         !        sgg%Med(contamedia)%wire(1)%RextremoK = this%twires%TW(j)%TWC(tama2)%k
         !        orientacionR                          = this%twires%TW(j)%TWC(tama2)%D
         !        !
         !!correct each ending
         !        numminus=0
         !        DO i = 2, tama2-1 !bug OLD 12/09/13  Model_unidos.nfde segmentos finales duplicados internamente
         !            punto%XI = this%twires%TW(j)%TWC(i)%i
         !            punto%YI = this%twires%TW(j)%TWC(i)%j
         !            punto%ZI = this%twires%TW(j)%TWC(i)%k
         !            orientacion = this%twires%TW(j)%TWC(i)%D
         !            select case (orientacion)
         !            case (iEx)
         !                if (                                      (punto%YI   == sgg%Med(contamedia)%wire(1)%LextremoJ).and.  &
         !                                                          (punto%ZI   == sgg%Med(contamedia)%wire(1)%LextremoK)) then
         !                    if ((orientacion /= orientacionL).and.(punto%XI   == sgg%Med(contamedia)%wire(1)%LextremoI)) numminus=numminus +1 !bug OLD 12/09/13  Model_unidos.nfde segmentos finales duplicados internamente
         !                    if                                    (punto%XI+1 == sgg%Med(contamedia)%wire(1)%LextremoI) numminus =numminus  +1
         !                endif
         !            case (iEy)
         !                if (                                      (punto%XI   == sgg%Med(contamedia)%wire(1)%LextremoI).and.  &
         !                                                          (punto%ZI   == sgg%Med(contamedia)%wire(1)%LextremoK)) then
         !                    if ((orientacion /= orientacionL).and.(punto%YI   == sgg%Med(contamedia)%wire(1)%LextremoJ)) numminus=numminus +1
         !                    if                                    (punto%YI+1 == sgg%Med(contamedia)%wire(1)%LextremoJ) numminus =numminus  +1
         !                endif
         !            case (iEz)
         !                if (                                      (punto%YI   == sgg%Med(contamedia)%wire(1)%LextremoJ).and.  &
         !                                                          (punto%XI   == sgg%Med(contamedia)%wire(1)%LextremoI)) then
         !                    if ((orientacion /= orientacionL).and.(punto%ZI   == sgg%Med(contamedia)%wire(1)%LextremoK)) numminus=numminus +1
         !                    if                                    (punto%ZI+1 == sgg%Med(contamedia)%wire(1)%LextremoK) numminus =numminus  +1
         !                endif
         !            end select
         !
         !        END DO
         !        if (numminus >= 1) then !bug OLD 12/09/13  Model_unidos.nfde segmentos finales duplicados internamente
         !              select case (this%twires%TW(j)%TWC(1)%D)
         !              case (iEx) !si son iguales a 2 es cerrado
         !                  sgg%Med(contamedia)%wire(1)%LextremoI = sgg%Med(contamedia)%wire(1)%LextremoI + 1
         !              case (iEy)
         !                  sgg%Med(contamedia)%wire(1)%LextremoJ = sgg%Med(contamedia)%wire(1)%LextremoJ + 1
         !              case (iEz)
         !                  sgg%Med(contamedia)%wire(1)%LextremoK = sgg%Med(contamedia)%wire(1)%LextremoK + 1
         !              end select
         !        endif
         !        !
         !!correct each ending
         !        numminus=0
         !        DO i = 2, tama2-1 !bug OLD 12/09/13  Model_unidos.nfde segmentos finales duplicados internamente
         !            punto%XI = this%twires%TW(j)%TWC(i)%i
         !            punto%YI = this%twires%TW(j)%TWC(i)%j
         !            punto%ZI = this%twires%TW(j)%TWC(i)%k
         !            orientacion = this%twires%TW(j)%TWC(i)%D
         !            select case (orientacion)
         !            case (iEx)
         !                if (                                      (punto%YI   == sgg%Med(contamedia)%wire(1)%RextremoJ).and.  &
         !                                                          (punto%ZI   == sgg%Med(contamedia)%wire(1)%RextremoK)) then
         !                    if ((orientacion /= orientacionR).and.(punto%XI   == sgg%Med(contamedia)%wire(1)%RextremoI)) numminus=numminus +1
         !                    if                                    (punto%XI+1 == sgg%Med(contamedia)%wire(1)%RextremoI) numminus =numminus  +1
         !                endif
         !            case (iEy)
         !                if (                                      (punto%XI   == sgg%Med(contamedia)%wire(1)%RextremoI).and.  &
         !                                                          (punto%ZI   == sgg%Med(contamedia)%wire(1)%RextremoK)) then
         !                    if ((orientacion /= orientacionR).and.(punto%YI   == sgg%Med(contamedia)%wire(1)%RextremoJ)) numminus=numminus +1
         !                    if                                    (punto%YI+1 == sgg%Med(contamedia)%wire(1)%RextremoJ) numminus =numminus  +1
         !                endif
         !            case (iEz)
         !                if (                                      (punto%YI   == sgg%Med(contamedia)%wire(1)%RextremoJ).and.  &
         !                                                          (punto%XI   == sgg%Med(contamedia)%wire(1)%RextremoI)) then
         !                    if ((orientacion /= orientacionR).and.(punto%ZI   == sgg%Med(contamedia)%wire(1)%RextremoK)) numminus=numminus +1
         !                    if                                    (punto%ZI+1 == sgg%Med(contamedia)%wire(1)%RextremoK) numminus =numminus  +1
         !                endif
         !            end select
         !        END DO
         !        if ((numminus >= 1).or.(tama2 == 1)) then  !bug ca295 !bug OLD 12/09/13  Model_unidos.nfde segmentos finales duplicados internamente
         !              select case (this%twires%TW(j)%TWC(tama2)%D)
         !              case (iEx)
         !                  sgg%Med(contamedia)%wire(1)%RextremoI = sgg%Med(contamedia)%wire(1)%RextremoI  + 1
         !!si son iguales a 2 es cerrado
         !              case (iEy)
         !                  sgg%Med(contamedia)%wire(1)%RextremoJ = sgg%Med(contamedia)%wire(1)%RextremoJ  + 1
         !              case (iEz)
         !                  sgg%Med(contamedia)%wire(1)%RextremoK = sgg%Med(contamedia)%wire(1)%RextremoK  + 1
         !              end select
         !        endif
      end do !del tama


      !preanalisis de hilos embeddeds en materiales  antes de asignarlos
      tama = this%twires%n_tw
      paraerrhilo=.false.
      do j1=1, tama
         tama2 = this%twires%TW(j1)%N_TWC
         DO i1 = 1, tama2
            i=this%twires%TW(j1)%TWC(i1)%i
            j=this%twires%TW(j1)%TWC(i1)%j
            k=this%twires%TW(j1)%TWC(i1)%k
            orientacion = this%twires%TW(j1)%TWC(i1)%D
            OrigIndex=    this%twires%TW(j1)%TWC(i1)%nd
            IF ((i >= BoundingBox%XI) .AND. (i < BoundingBox%XE) .AND. &
            &    (j >= BoundingBox%YI) .AND. (j < BoundingBox%YE) .AND. &
            &    (k >= BoundingBox%ZI) .AND. (k < BoundingBox%ZE)) THEN
               if (i > BoundingBox%XI) then
                  imenos1=i-1
               else
                  imenos1=i
               endif
               if (j > BoundingBox%YI) then
                  jmenos1=j-1
               else
                  jmenos1=j
               endif
               if (k > BoundingBox%ZI) then
                  kmenos1=k-1
               else
                  kmenos1=k
               endif
               select case (orientacion)
                case (iEx)
                  if ((media%sggMiEx(i,j,k) ==0).or.(sgg%med(media%sggMiEx(i,j,k) )%is%pec)) then
                     paraerrhilo=.true.
                     WRITE (buff, '(a,i7,3i5,a)')    'pre1_WARNING:   x-WIRE at ',OrigIndex, i, j, k,' embedded within PEC'
                     if (verbose) CALL WarnErrReport (buff)
                  elseif (media%sggMiEx(i,j,k) /= 1) then
                     islossy = (sgg%Med(media%sggMiEx(i,j,k))%Sigma /= 0.0_RKIND)
                     if (islossy) then
                        paraerrhilo=.true.
                        WRITE (buff, '(a,i7,3i5,a,i5)') 'pre1_WARNING: x-WIRE at ',OrigIndex, i, j, k, &
                           ' embedded within LOSSY medium ', &
                           media%sggMiEx(i,j,k)
                     else
                        WRITE (buff, '(a,i7,3i5,a,i5)') 'pre1_WARNING: x-WIRE at ',OrigIndex, i, j, k,' embedded within medium ', &
                           media%sggMiEx(i,j,k)
                     endif
                     if (verbose) CALL WarnErrReport (buff)
                  endif
                  if ((((media%sggMiEy(i  ,j,k) ==0).or.(sgg%med(media%sggMiEy(i  ,j,k) )%is%pec)).or. &
                     ((media%sggMiEz(i  ,j,k) ==0).or.(sgg%med(media%sggMiEz(i  ,j,k) )%is%pec)).or. &
                     ((media%sggMiEy(i  ,jmenos1,k) ==0).or.(sgg%med(media%sggMiEy(i  ,jmenos1,k) )%is%pec)).or. &
                     ((media%sggMiEz(i  ,j,kmenos1) ==0).or.(sgg%med(media%sggMiEz(i  ,j,kmenos1) )%is%pec))).and. &
                  &     ((media%sggMiEx(i  ,j,k) /=0).and.(.not.(sgg%med(media%sggMiEx(i  ,j,k) )%is%pec)))) then
                     if ((i1 /= 1) .and. (i1 /= tama2)) then !solo en LeftEnd y RightEnd pueden tocar
                        paraerrhilo=.true.
                        WRITE (buff, '(a,i7,3i5,a)')    'pre1_WARNING:   intermediate node of x-WIRE at  ',OrigIndex, i, j, k, &
                           ' touching PEC'
                        if (verbose) CALL WarnErrReport (buff)
                     else
                        continue
                        !WRITE (buff, '(a,i7,3i5,a)')    'A node of terminal x-WIRE at ',OrigIndex, i, j, k,' touching PEC'
                        !if (verbose) CALL WarnErrReport (buff)
                     endif
                  elseif ((((media%sggMiEy(i+1,j,k) ==0).or.(sgg%med(media%sggMiEy(i+1,j,k) )%is%pec)).or.&
                     ((media%sggMiEz(i+1,j,k) ==0).or.(sgg%med(media%sggMiEz(i+1,j,k) )%is%pec)).or. &
                     ((media%sggMiEy(i+1,jmenos1,k) ==0).or.(sgg%med(media%sggMiEy(i+1,jmenos1,k) )%is%pec)).or. &
                     ((media%sggMiEz(i+1,j,kmenos1) ==0).or.(sgg%med(media%sggMiEz(i+1,j,kmenos1) )%is%pec))).and. &
                  &         ((media%sggMiEx(i  ,j,k) /=0).and.(.not.(sgg%med(media%sggMiEx(i  ,j,k) )%is%pec)))) then
                     if ((i1 /= 1) .and. (i1 /= tama2)) then !solo en LeftEnd y RightEnd pueden tocar
                        paraerrhilo=.true.
                        WRITE (buff, '(a,i7,3i5,a)')    'pre1_WARNING:   intermediate node of x-WIRE at  ',OrigIndex, i+1, j, k, &
                           ' touching PEC'
                        if (verbose) CALL WarnErrReport (buff)
                     else
                        continue
                        !WRITE (buff, '(a,i7,3i5,a)')    'A node of terminal x-WIRE at ',OrigIndex, i+1, j, k,' touching PEC'
                        !if (verbose) CALL WarnErrReport (buff)
                     endif
                  elseif (((media%sggMiEy(i  ,j,k) /= 1)).and. &
                     (media%sggMiEx(i  ,j,k) == 1)) then
                     WRITE (buff, '(a,i7,3i5,a,i5)') 'pre1_WARNING: x-WIRE at ',OrigIndex, i, j, k,' touching medium ', &
                     &                                  media%sggMiEy(i,j,k)
                     if (verbose) CALL WarnErrReport (buff)
                  elseif (((media%sggMiEz(i  ,j,k) /= 1)).and. &
                     (media%sggMiEx(i  ,j,k) == 1)) then
                     WRITE (buff, '(a,i7,3i5,a,i5)') 'pre1_WARNING: x-WIRE at ',OrigIndex, i, j, k,' touching medium ', &
                     &                                  media%sggMiEz(i,j,k)
                     if (verbose) CALL WarnErrReport (buff)
                  elseif (((media%sggMiEy(i+1,j,k) /= 1)).and. &
                     (media%sggMiEx(i  ,j,k) == 1)) then
                     WRITE (buff, '(a,i7,3i5,a,i5)') 'pre1_WARNING: x-WIRE at ',OrigIndex, i+1, j, k,' touching medium ', &
                     &                                  media%sggMiEy(i+1,j,k)
                     if (verbose) CALL WarnErrReport (buff)
                  elseif (((media%sggMiEz(i+1,j,k) /= 1)).and. &
                     (media%sggMiEx(i  ,j,k) == 1)) then
                     WRITE (buff, '(a,i7,3i5,a,i5)') 'pre1_WARNING: x-WIRE at ',OrigIndex, i+1, j, k,' touching medium ', &
                     &                                  media%sggMiEz(i+1,j,k)
                     if (verbose) CALL WarnErrReport (buff)
                  endif
                case (iEy)
                  if ((media%sggMiEy(i,j,k) ==0).or.(sgg%med(media%sggMiEy(i,j,k) )%is%pec)) then
                     paraerrhilo=.true.
                     WRITE (buff, '(a,i7,3i5,a)')    'pre1_WARNING:   y-WIRE at ',OrigIndex, i, j, k,' embedded within PEC'
                     if (verbose) CALL WarnErrReport (buff)
                  elseif (media%sggMiEy(i,j,k) /= 1) then
                     islossy = (sgg%Med(media%sggMiEy(i,j,k))%Sigma /= 0.0_RKIND)
                     if (islossy) then
                        paraerrhilo=.true.
                        WRITE (buff, '(a,i7,3i5,a,i5)') 'pre1_WARNING: Y-WIRE at ',OrigIndex, i, j, k,' embedded within LOSSY medium ', &
                           media%sggMiEY(i,j,k)
                     else
                        WRITE (buff, '(a,i7,3i5,a,i5)') 'pre1_WARNING: y-WIRE at ',OrigIndex, i, j, k,' embedded within medium ', &
                           media%sggMiEy(i,j,k)
                     ENDIF
                     if (verbose) CALL WarnErrReport (buff)
                  endif
                  if ((((media%sggMiEx(i,j  ,k) ==0).or.(sgg%med(media%sggMiEx(i,j  ,k) )%is%pec)).or. &
                     ((media%sggMiEz(i,j,k  ) ==0).or.(sgg%med(media%sggMiEz(i,j,k  ) )%is%pec)).or. &
                     ((media%sggMiEx(imenos1,j  ,k) ==0).or.(sgg%med(media%sggMiEx(imenos1,j  ,k) )%is%pec)).or. &
                     ((media%sggMiEz(i,j,kmenos1  ) ==0).or.(sgg%med(media%sggMiEz(i,j,kmenos1  ) )%is%pec))).and. &
                  &     ((media%sggMiEy(i,j  ,k) /=0).and.(.not.(sgg%med(media%sggMiEy(i,j  ,k) )%is%pec)))) then
                     if ((i1 /= 1) .and. (i1 /= tama2)) then !solo en LeftEnd y RightEnd pueden tocar
                        paraerrhilo=.true.
                        WRITE (buff, '(a,i7,3i5,a)')    'pre1_WARNING:   intermediate node of y-WIRE at ',OrigIndex, i, j, k, &
                           ' touching PEC'
                        if (verbose) CALL WarnErrReport (buff)
                     else
                        continue
                        !WRITE (buff, '(a,i7,3i5,a)')    'A node of terminal x-WIRE at ',OrigIndex, i, j, k,' touching PEC'
                        !if (verbose) CALL WarnErrReport (buff)
                     endif
                  elseif ((((media%sggMiEx(i,j+1,k) ==0).or.(sgg%med(media%sggMiEx(i,j+1,k) )%is%pec)).or. &
                     ((media%sggMiEz(i,j+1,k) ==0).or.(sgg%med(media%sggMiEz(i,j+1,k) )%is%pec)).or. &
                     ((media%sggMiEx(imenos1,j+1,k) ==0).or.(sgg%med(media%sggMiEx(imenos1,j+1,k) )%is%pec)).or. &
                     ((media%sggMiEz(i,j+1,kmenos1) ==0).or.(sgg%med(media%sggMiEz(i,j+1,kmenos1) )%is%pec))).and. &
                  &         ((media%sggMiEy(i,j  ,k) /=0).and.(.not.(sgg%med(media%sggMiEy(i,j  ,k) )%is%pec)))) then
                     if ((i1 /= 1) .and. (i1 /= tama2)) then !solo en LeftEnd y RightEnd pueden tocar
                        paraerrhilo=.true.
                        WRITE (buff, '(a,i7,3i5,a)')    'pre1_WARNING:   intermediate node of y-WIRE at ',OrigIndex, i, j+1, k, &
                           ' touching PEC'
                        if (verbose) CALL WarnErrReport (buff)
                     else
                        continue
                        !WRITE (buff, '(a,i7,3i5,a)')    'A node of terminal x-WIRE at ',OrigIndex, i, j+1, k,' touching PEC'
                        !if (verbose) CALL WarnErrReport (buff)
                     endif
                  elseif (((media%sggMiEx(i,j  ,k) /= 1)).and. &
                  &         (media%sggMiEy(i,j  ,k) == 1)) then
                     WRITE (buff, '(a,i7,3i5,a,i5)') 'pre1_WARNING: y-WIRE at ',OrigIndex, i, j, k,' touching medium ', &
                     &                                  media%sggMiEx(i,j,k)
                     if (verbose) CALL WarnErrReport (buff)
                  elseif (((media%sggMiEz(i,j  ,k) /= 1)).and. &
                  &         (media%sggMiEy(i,j  ,k) == 1)) then
                     WRITE (buff, '(a,i7,3i5,a,i5)') 'pre1_WARNING: y-WIRE at ',OrigIndex, i, j, k,' touching medium ', &
                     &                                  media%sggMiEz(i,j,k)
                     if (verbose) CALL WarnErrReport (buff)
                  elseif (((media%sggMiEx(i,j+1,k) /= 1)).and. &
                  &         (media%sggMiEy(i,j  ,k) == 1)) then
                     WRITE (buff, '(a,i7,3i5,a,i5)') 'pre1_WARNING: y-WIRE at ',OrigIndex, i, j+1, k,' touching medium ', &
                     &                                  media%sggMiEx(i,j+1,k)
                     if (verbose) CALL WarnErrReport (buff)
                  elseif (((media%sggMiEz(i,j+1,k) /= 1)).and. &
                  &         (media%sggMiEy(i,j  ,k) == 1)) then
                     WRITE (buff, '(a,i7,3i5,a,i5)') 'pre1_WARNING: y-WIRE at ',OrigIndex, i, j+1, k,' touching medium ', &
                     &                                  media%sggMiEz(i,j+1,k)
                     if (verbose) CALL WarnErrReport (buff)
                  endif
                case (iEz)
                  if ((media%sggMiEz(i,j,k) ==0).or.(sgg%med(media%sggMiEz(i,j,k) )%is%pec)) then
                     paraerrhilo=.true.
                     WRITE (buff, '(a,i7,3i5,a)')    'pre1_WARNING:   z-WIRE at ',OrigIndex, i, j, k,' embedded within PEC'
                     if (verbose) CALL WarnErrReport (buff)
                  elseif (media%sggMiEz(i,j,k) /= 1) then
                     islossy = (sgg%Med(media%sggMiEz(i,j,k))%Sigma /= 0.0_RKIND)
                     if (islossy) then
                        paraerrhilo=.true.
                        WRITE (buff, '(a,i7,3i5,a,i5)') 'pre1_WARNING: Y-WIRE at ',OrigIndex, i, j, k, &
                           ' embedded within LOSSY medium ', &
                           media%sggMiEz(i,j,k)
                     else
                        WRITE (buff, '(a,i7,3i5,a,i5)') 'pre1_WARNING: z-WIRE at ',OrigIndex, i, j, k,' embedded within medium ', &
                           media%sggMiEz(i,j,k)
                     ENDIF
                     if (verbose) CALL WarnErrReport (buff)
                  endif
                  if ((((media%sggMiEx(i,j,k  ) ==0).or.(sgg%med(media%sggMiEx(i,j,k  ) )%is%pec)).or. &
                     ((media%sggMiEy(i,j,k  ) ==0).or.(sgg%med(media%sggMiEy(i,j,k  ) )%is%pec)).or. &
                     ((media%sggMiEx(imenos1,j,k  ) ==0).or.(sgg%med(media%sggMiEx(imenos1,j,k  ) )%is%pec)).or. &
                     ((media%sggMiEy(i,jmenos1,k  ) ==0).or.(sgg%med(media%sggMiEy(i,jmenos1,k  ) )%is%pec))).and. &
                  &     ((media%sggMiEz(i,j,k  ) /=0).and.(.not.(sgg%med(media%sggMiEz(i,j,k  ) )%is%pec)))) then
                     if ((i1 /= 1) .and. (i1 /= tama2)) then !solo en LeftEnd y RightEnd pueden tocar
                        paraerrhilo=.true.
                        WRITE (buff, '(a,i7,3i5,a)')    'pre1_WARNING:   intermediate node of z-WIRE at ',OrigIndex, i, j, k, &
                           ' touching PEC'
                        if (verbose)  CALL WarnErrReport (buff)
                     else
                        continue
                        !WRITE (buff, '(a,i7,3i5,a)')    'A node of terminal x-WIRE at ',OrigIndex, i, j, k,' touching PEC'
                        !if (verbose) CALL WarnErrReport (buff)
                     endif
                  elseif ((((media%sggMiEx(i,j  ,k+1) ==0).or.(sgg%med(media%sggMiEx(i,j  ,k+1) )%is%pec)).or. &
                     ((media%sggMiEy(i,j  ,k+1) ==0).or.(sgg%med(media%sggMiEy(i,j  ,k+1) )%is%pec)).or.   &
                  &         ((media%sggMiEx(imenos1,j,k+1) ==0).or.(sgg%med(media%sggMiEx(imenos1,j,k+1) )%is%pec)).or. &
                     ((media%sggMiEy(i,jmenos1,k+1) ==0).or.(sgg%med(media%sggMiEy(i,jmenos1,k+1) )%is%pec))).and. &
                  &         (((media%sggMiEz(i,j,k  ) /=0).and.(.not.(sgg%med(media%sggMiEz(i,j,k  ) )%is%pec))).or.(.not.(sgg%med(media%sggMiEz(i,j,k  ) )%is%pec)))) then
                     if ((i1 /= 1) .and. (i1 /= tama2)) then !solo en LeftEnd y RightEnd pueden tocar
                        paraerrhilo=.true.
                        WRITE (buff, '(a,i7,3i5,a)')    'pre1_WARNING:   intermediate node of z-WIRE at ',OrigIndex, i, j, k+1, &
                           ' touching PEC'
                        if (verbose) CALL WarnErrReport (buff)
                     else
                        continue
                        !WRITE (buff, '(a,i7,3i5,a)')    'A node of terminal x-WIRE at ',OrigIndex, i, j, k+1,' touching PEC'
                        !if (verbose) CALL WarnErrReport (buff)
                     endif
                  elseif (((media%sggMiEx(i,j,k  ) /= 1)).and. &
                  &         (media%sggMiEz(i,j,k  ) == 1)) then
                     WRITE (buff, '(a,i7,3i5,a,i5)') 'pre1_WARNING: z-WIRE at ',OrigIndex, i, j, k,' touching medium ', &
                     &                                  media%sggMiEx(i,j,k)
                     if (verbose) CALL WarnErrReport (buff)
                  elseif (((media%sggMiEy(i,j,k  ) /= 1)).and. &
                  &         (media%sggMiEz(i,j,k  ) == 1)) then
                     WRITE (buff, '(a,i7,3i5,a,i5)') 'pre1_WARNING: z-WIRE at ',OrigIndex, i, j, k,' touching medium ', &
                     &                                  media%sggMiEy(i,j,k)
                     if (verbose) CALL WarnErrReport (buff)
                  elseif (((media%sggMiEy(i,j,k+1) /= 1)).and. &
                  &         (media%sggMiEz(i,j,k  ) == 1)) then
                     WRITE (buff, '(a,i7,3i5,a,i5)') 'pre1_WARNING: z-WIRE at ',OrigIndex, i, j, k,' touching medium ', &
                     &                                  media%sggMiEy(i,j,k+1)
                     if (verbose) CALL WarnErrReport (buff)
                  elseif (((media%sggMiEx(i,j,k+1) /= 1)).and. &
                  &         (media%sggMiEz(i,j,k  ) == 1)) then
                     WRITE (buff, '(a,i7,3i5,a,i5)') 'pre1_WARNING: z-WIRE at ',OrigIndex, i, j, k,' touching medium ', &
                     &                                  media%sggMiEx(i,j,k+1)
                     if (verbose) CALL WarnErrReport (buff)
                  endif
               end select
            endif
         end do
      end do
      !lo dejo que siga !luego el wires parara
      !if (paraerrhilo.and.(.not.groundwires)) then
      !    buff='Revise WIRE intersections!'
      !    call STOPONERROR(layoutnumber,size,buff)
      !endif
      !end preanalisis
      !perform the assignments
      bboxwirxi=  2**20
      bboxwirxe=-(2**20)
      bboxwirYi=  2**20
      bboxwirYe=-(2**20)
      bboxwirZi=  2**20
      bboxwirZe=-(2**20)
      contamedia = oldcontamedia
      tama = this%twires%n_tw
      DO j = 1, tama
         contamedia = contamedia + 1
         tama2 = this%twires%TW(j)%N_TWC
         TAMA2BIS=0
         DO i = 1, tama2
            punto%XI = this%twires%TW(j)%TWC(i)%i
            punto%XE = this%twires%TW(j)%TWC(i)%i
            punto%YI = this%twires%TW(j)%TWC(i)%j
            punto%YE = this%twires%TW(j)%TWC(i)%j
            punto%ZI = this%twires%TW(j)%TWC(i)%k
            punto%ZE = this%twires%TW(j)%TWC(i)%k
            !!!!!!!!!los clipeo agresivamente si lo lanzo con -CLIPREGION para que no me den problema 06/07/15 (solo sirve para debugeo y con el -wiresflavor holland (old))
            if (((punto%XI-2 >  SINPML_fullsize(iHx)%XI).and.(punto%XI+2 < SINPML_fullsize(iHx)%XE).and. &
               (punto%YI-2 >  SINPML_fullsize(iHy)%YI).and.(punto%YI+2 < SINPML_fullsize(iHy)%YE).and. &
               (punto%ZI-2 >  SINPML_fullsize(iHz)%ZI).and.(punto%zI+2 < SINPML_fullsize(iHz)%ZE).and. &
               (punto%XE-2 >  SINPML_fullsize(iHx)%XI).and.(punto%Xe+2 < SINPML_fullsize(iHx)%XE).and. &
               (punto%YE-2 >  SINPML_fullsize(iHy)%YI).and.(punto%Ye+2 < SINPML_fullsize(iHy)%YE).and. &
               (punto%ZE-2 >  SINPML_fullsize(iHz)%ZI).and.(punto%Ze+2 < SINPML_fullsize(iHz)%ZE)).or.(.not.CLIPREGION)) TAMA2BIS=TAMA2bis+1
         end do
         !
         ALLOCATE (sgg%Med(contamedia)%wire(1)%segm(1:TAMA2BIS))
         sgg%Med(contamedia)%wire(1)%numsegmentos = TAMA2BIS
         ALLOCATE (sgg%Med(contamedia)%wire(1)%VSource(1:TAMA2BIS))
         ALLOCATE (sgg%Med(contamedia)%wire(1)%ISource(1:TAMA2BIS))
         CONTAVOLT=0
         CONTACURR=0

         TAMA2BIS=0
         hilosbarre: DO i = 1, tama2
            punto%XI = this%twires%TW(j)%TWC(i)%i
            punto%XE = this%twires%TW(j)%TWC(i)%i
            punto%YI = this%twires%TW(j)%TWC(i)%j
            punto%YE = this%twires%TW(j)%TWC(i)%j
            punto%ZI = this%twires%TW(j)%TWC(i)%k
            punto%ZE = this%twires%TW(j)%TWC(i)%k
!!!sgg250418
!!!bug 2018
!!bug que aparece cuando en hilos de dos segmentos con ambos identicos. Lo que hago es clipearlo directamente.
            if ((i==2).and.(tama2==2)) then
               if  (((this%twires%TW(j)%TWC(i)%i).eq.(this%twires%TW(j)%TWC(i-1)%i)).and. &
                  ((this%twires%TW(j)%TWC(i)%j).eq.(this%twires%TW(j)%TWC(i-1)%j)).and. &
                  ((this%twires%TW(j)%TWC(i)%k).eq.(this%twires%TW(j)%TWC(i-1)%k)).and. &
                  ((this%twires%TW(j)%TWC(i)%D).eq.(this%twires%TW(j)%TWC(i-1)%D))) then
                  sgg%Med(contamedia)%wire(1)%numsegmentos = 1
                  WRITE (buff, '(a,1i7,a)')    'pre1_SEVEREWARNING: removing a repeteated segment from a 2-segment wire. Index= ',this%twires%TW(j)%TWC(i)%nd,'. Double check that no wire probes was attached to it'
                  CALL WarnErrReport (buff)
                  exit hilosbarre
               endif
            endif
!!!fin 250418

            !!!!!!!!!los clipeo agresivamente si lo lanzo con -CLIPREGION para que no me den problema 06/07/15 (solo sirve para debugeo y con el -wiresflavor holland (old))
            if (.not. &
               (((punto%XI-2 >  SINPML_fullsize(iHx)%XI).and.(punto%XI+2 < SINPML_fullsize(iHx)%XE).and. &
               (punto%YI-2 >  SINPML_fullsize(iHy)%YI).and.(punto%YI+2 < SINPML_fullsize(iHy)%YE).and. &
               (punto%ZI-2 >  SINPML_fullsize(iHz)%ZI).and.(punto%zI+2 < SINPML_fullsize(iHz)%ZE).and. &
               (punto%XE-2 >  SINPML_fullsize(iHx)%XI).and.(punto%Xe+2 < SINPML_fullsize(iHx)%XE).and. &
               (punto%YE-2 >  SINPML_fullsize(iHy)%YI).and.(punto%Ye+2 < SINPML_fullsize(iHy)%YE).and. &
               (punto%ZE-2 >  SINPML_fullsize(iHz)%ZI).and.(punto%Ze+2 < SINPML_fullsize(iHz)%ZE)).or.(.not.CLIPREGION)) ) CYCLE hilosbarre

            TAMA2BIS=TAMA2BIS+1
            orientacion = this%twires%TW(j)%TWC(i)%D
            origindex=    this%twires%TW(j)%TWC(i)%nd
            !
            !
            sgg%Med(contamedia)%wire(1)%SEGM(TAMA2BIS)%i = punto%XI
            sgg%Med(contamedia)%wire(1)%SEGM(TAMA2BIS)%j = punto%YI
            sgg%Med(contamedia)%wire(1)%SEGM(TAMA2BIS)%k = punto%ZI

            !!!2014 para informacion bbox hilos
            if (punto%XI < bboxwirxi) bboxwirXI=punto%XI
            if (punto%XE > bboxwirxE) bboxwirXE=punto%XE
            if (punto%YI < bboxwirYi) bboxwirYI=punto%YI
            if (punto%YE > bboxwirYE) bboxwirYE=punto%YE
            if (punto%ZI < bboxwirZi) bboxwirZI=punto%ZI
            if (punto%ZE > bboxwirZE) bboxwirZE=punto%ZE
            !!!!!!

            sgg%Med(contamedia)%wire(1)%SEGM(TAMA2BIS)%ori = orientacion
            sgg%Med(contamedia)%wire(1)%SEGM(TAMA2BIS)%origindex = origindex
            sgg%Med(contamedia)%wire(1)%SEGM(TAMA2BIS)%Is_LeftEnd = .false.
            sgg%Med(contamedia)%wire(1)%SEGM(TAMA2BIS)%Is_RightEnd = .false.
            sgg%Med(contamedia)%wire(1)%SEGM(TAMA2BIS)%repetido = .false. !luego el preprocesador del wires cambia esto
            if (i==1) sgg%Med(contamedia)%wire(1)%SEGM(TAMA2BIS)%Is_LeftEnd = .true.
            if (i==tama2) sgg%Med(contamedia)%wire(1)%SEGM(TAMA2BIS)%Is_RightEnd = .true.
            !
            isathinwire = .TRUE.
            numertag = searchtag(tagtype,this%twires%TW(j)%TWC(i)%tag)
            CALL CreateLineMM (layoutnumber, media%sggMtag, tag_numbers, numertag, media%sggMiEx, media%sggMiEy, media%sggMiEz, &
            & media%sggMiHx, media%sggMiHy, media%sggMiHz, Alloc_iEx_XI, &
            & Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI, &
            & Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI, &
            & Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI, &
            & Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI, &
            & Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg%Med, sgg%NumMedia, sgg%EShared, BoundingBox, punto, orientacion, &
            & contamedia, isathinwire,verbose,numeroasignaciones)
            IF ((trim(adjustl(this%twires%TW(j)%TWC(i)%SRCTYPE))   == 'VOLT') .OR.   &
            &     (trim(adjustl(this%twires%TW(j)%TWC(i)%SRCTYPE)) == 'CURR').OR.   &
            &     (trim(adjustl(this%twires%TW(j)%TWC(i)%SRCTYPE)) == 'SOFTVOLT').OR.   &
            &     (trim(adjustl(this%twires%TW(j)%TWC(i)%SRCTYPE)) == 'SOFTCURR').OR.   &
            &     (trim(adjustl(this%twires%TW(j)%TWC(i)%SRCTYPE)) == 'HARDVOLT').OR.   &
            &     (trim(adjustl(this%twires%TW(j)%TWC(i)%SRCTYPE)) == 'HARDCURR')) THEN
               !!!!!!!!!!!!!
               IF ((trim(adjustl(this%twires%TW(j)%TWC(i)%SRCTYPE)) == 'VOLT').OR. &
                  (trim(adjustl(this%twires%TW(j)%TWC(i)%SRCTYPE)) == 'SOFTVOLT')) THEN
                  CONTAVOLT=CONTAVOLT+1
                  sgg%Med(contamedia)%wire(1)%VsourceExists = .TRUE.
                  sgg%Med(contamedia)%wire(1)%VSource(CONTAVOLT)%SOFT=.TRUE.
                  sgg%Med(contamedia)%wire(1)%VSource(CONTAVOLT)%Multiplier = this%twires%TW(j)%TWC(i)%m
                  sgg%Med(contamedia)%wire(1)%VSource(CONTAVOLT)%Resistance = 0.0_RKIND
                  !not provided by .nfde but supported by the simulation though untested
                  !
                  sgg%Med(contamedia)%wire(1)%VSource(CONTAVOLT)%fichero%name = trim (adjustl(this%twires%TW(j)%TWC(i)%SRCFILE))
                  sgg%Med(contamedia)%wire(1)%VSource(CONTAVOLT)%i = this%twires%TW(j)%TWC(i)%i
                  sgg%Med(contamedia)%wire(1)%VSource(CONTAVOLT)%j = this%twires%TW(j)%TWC(i)%j
                  sgg%Med(contamedia)%wire(1)%VSource(CONTAVOLT)%k = this%twires%TW(j)%TWC(i)%k
               ELSE IF ((trim(adjustl(this%twires%TW(j)%TWC(i)%SRCTYPE)) == 'CURR').OR. &
                  (trim(adjustl(this%twires%TW(j)%TWC(i)%SRCTYPE)) == 'HARDCURR')) THEN
                  CONTAVOLT=CONTAVOLT+1
                  CALL WarnErrReport (buff)
                  !             CALL STOPONERROR(layoutnumber,size,buff)
                  sgg%Med(contamedia)%wire(1)%VsourceExists = .TRUE.
                  sgg%Med(contamedia)%wire(1)%VSource(CONTAVOLT)%SOFT=.FALSE. !230323 FUENTES DURAS Y FALSAS SGG
                  sgg%Med(contamedia)%wire(1)%VSource(CONTAVOLT)%Multiplier = this%twires%TW(j)%TWC(i)%m
                  sgg%Med(contamedia)%wire(1)%VSource(CONTAVOLT)%fichero%name = trim (adjustl(this%twires%TW(j)%TWC(i)%SRCFILE))
                  sgg%Med(contamedia)%wire(1)%VSource(CONTAVOLT)%Resistance = 0.0_RKIND
                  !not provided by .nfde but supported by the simulation though untested
                  !
                  sgg%Med(contamedia)%wire(1)%VSource(CONTAVOLT)%fichero%name = trim (adjustl(this%twires%TW(j)%TWC(i)%SRCFILE))
                  sgg%Med(contamedia)%wire(1)%VSource(CONTAVOLT)%i = this%twires%TW(j)%TWC(i)%i
                  sgg%Med(contamedia)%wire(1)%VSource(CONTAVOLT)%j = this%twires%TW(j)%TWC(i)%j
                  sgg%Med(contamedia)%wire(1)%VSource(CONTAVOLT)%k = this%twires%TW(j)%TWC(i)%k
               ELSEIF (trim(adjustl(this%twires%TW(j)%TWC(i)%SRCTYPE)) == 'HARDVOLT') THEN
                  CONTACURR=CONTACURR+1
                  sgg%Med(contamedia)%wire(1)%IsourceExists = .TRUE.
                  sgg%Med(contamedia)%wire(1)%ISource(CONTACURR)%SOFT=.FALSE. !230323 FUENTES DURAS Y FALSAS SGG
                  sgg%Med(contamedia)%wire(1)%ISource(CONTACURR)%Multiplier = this%twires%TW(j)%TWC(i)%m
                  sgg%Med(contamedia)%wire(1)%ISource(CONTACURR)%Resistance = 0.0_RKIND
                  !not provided by .nfde but supported by the simulation though untested
                  !
                  sgg%Med(contamedia)%wire(1)%ISource(CONTACURR)%fichero%name = trim (adjustl(this%twires%TW(j)%TWC(i)%SRCFILE))
                  sgg%Med(contamedia)%wire(1)%ISource(CONTACURR)%i = this%twires%TW(j)%TWC(i)%i
                  sgg%Med(contamedia)%wire(1)%ISource(CONTACURR)%j = this%twires%TW(j)%TWC(i)%j
                  sgg%Med(contamedia)%wire(1)%ISource(CONTACURR)%k = this%twires%TW(j)%TWC(i)%k
               ELSE IF ((trim(adjustl(this%twires%TW(j)%TWC(i)%SRCTYPE)) == 'SOFTCURR')) THEN
                  CONTACURR=CONTACURR+1
                  CALL WarnErrReport (buff)
                  !             CALL STOPONERROR(layoutnumber,size,buff)
                  sgg%Med(contamedia)%wire(1)%IsourceExists = .TRUE.
                  sgg%Med(contamedia)%wire(1)%ISource(CONTACURR)%SOFT=.true. !230323 FUENTES DURAS Y FALSAS SGG
                  sgg%Med(contamedia)%wire(1)%ISource(CONTACURR)%Multiplier = this%twires%TW(j)%TWC(i)%m
                  sgg%Med(contamedia)%wire(1)%ISource(CONTACURR)%fichero%name = trim (adjustl(this%twires%TW(j)%TWC(i)%SRCFILE))
                  sgg%Med(contamedia)%wire(1)%ISource(CONTACURR)%Resistance = 0.0_RKIND
                  !not provided by .nfde but supported by the simulation though untested
                  !
                  sgg%Med(contamedia)%wire(1)%ISource(CONTACURR)%fichero%name = trim (adjustl(this%twires%TW(j)%TWC(i)%SRCFILE))
                  sgg%Med(contamedia)%wire(1)%ISource(CONTACURR)%i = this%twires%TW(j)%TWC(i)%i
                  sgg%Med(contamedia)%wire(1)%ISource(CONTACURR)%j = this%twires%TW(j)%TWC(i)%j
                  sgg%Med(contamedia)%wire(1)%ISource(CONTACURR)%k = this%twires%TW(j)%TWC(i)%k
               END IF
            ELSEIF (trim(adjustl(this%twires%TW(j)%TWC(i)%SRCTYPE)) /= 'None') THEN
               write(buff,*) 'WRONG type of wire source '//trim(adjustl(this%twires%TW(j)%TWC(i)%SRCTYPE))
               call stoponerror (layoutnumber,size,buff)
            END IF
         END DO hilosbarre
         sgg%Med(contamedia)%wire(1)%NUMVOLTAGESOURCES=CONTAVOLT
         sgg%Med(contamedia)%wire(1)%NUMCURRENTSOURCES=CONTACURR
         !wires
      END DO

!!!ahora slanted

      if (this%swires%n_sw /=0) then
         hay_slanted_wires=.true.
      else
         hay_slanted_wires=.false.
      endif

      !SLANTED WIRES
      do j = 1,this%swires%n_sw
         contamedia = contamedia + 1
         ALLOCATE (sgg%Med(contamedia)%SlantedWire(1))

         sgg%Med(contamedia)%Epr    = sgg%Med(1)%Epr
         sgg%Med(contamedia)%Sigma  = sgg%Med(1)%Sigma
         sgg%Med(contamedia)%Mur    = sgg%Med(1)%Mur
         sgg%Med(contamedia)%SigmaM = sgg%Med(1)%SigmaM
         sgg%Med(contamedia)%Is%SlantedWire = .TRUE.
         sgg%Med(contamedia)%Is%Dielectric  = .FALSE.
         sgg%Med(contamedia)%SlantedWire(1)%radius = this%swires%SW(j)%RAD
         if (boundwireradius) then
            if (sgg%Med(contamedia)%SlantedWire(1)%radius > maxwireradius) sgg%Med(contamedia)%SlantedWire(1)%radius=maxwireradius
         endif
         sgg%Med(contamedia)%SlantedWire(1)%R = this%swires%SW(j)%res
         sgg%Med(contamedia)%SlantedWire(1)%L = this%swires%SW(j)%ind
         sgg%Med(contamedia)%SlantedWire(1)%C = this%swires%SW(j)%cap
         sgg%Med(contamedia)%SlantedWire(1)%P_R = this%swires%SW(j)%P_res
         sgg%Med(contamedia)%SlantedWire(1)%P_L = this%swires%SW(j)%P_ind
         sgg%Med(contamedia)%SlantedWire(1)%P_C = this%swires%SW(j)%P_cap

         if (this%swires%SW(j)%disp) then
            allocate (sgg%Med(contamedia)%SlantedWire(1)%disp(1))
            call asignawiredisper(sgg%Med(contamedia)%SlantedWire(1)%disp(1), &
               this%swires%SW(j)%dispfile)
         end if

         sgg%Med(contamedia)%SlantedWire(1)%LeftEnd = this%swires%SW(j)%LeftEnd
         sgg%Med(contamedia)%SlantedWire(1)%RightEnd = this%swires%SW(j)%RightEnd

         sgg%Med(contamedia)%SlantedWire(1)%HasParallel_LeftEnd = .FALSE.
         sgg%Med(contamedia)%SlantedWire(1)%Parallel_R_LeftEnd = 0.0_RKIND
         sgg%Med(contamedia)%SlantedWire(1)%Parallel_L_LeftEnd = 0.0_RKIND
         sgg%Med(contamedia)%SlantedWire(1)%Parallel_C_LeftEnd = 0.0_RKIND
         sgg%Med(contamedia)%SlantedWire(1)%HasParallel_RightEnd = .FALSE.
         sgg%Med(contamedia)%SlantedWire(1)%Parallel_R_RightEnd = 0.0_RKIND
         sgg%Med(contamedia)%SlantedWire(1)%Parallel_L_RightEnd = 0.0_RKIND
         sgg%Med(contamedia)%SlantedWire(1)%Parallel_C_RightEnd = 0.0_RKIND
         sgg%Med(contamedia)%SlantedWire(1)%HasSeries_LeftEnd = .FALSE.
         sgg%Med(contamedia)%SlantedWire(1)%Series_R_LeftEnd = 0.0_RKIND
         sgg%Med(contamedia)%SlantedWire(1)%Series_L_LeftEnd = 0.0_RKIND
         sgg%Med(contamedia)%SlantedWire(1)%Series_C_LeftEnd = 2.0e7_RKIND !en corto 14/2/14
         sgg%Med(contamedia)%SlantedWire(1)%HasSeries_RightEnd = .FALSE.
         sgg%Med(contamedia)%SlantedWire(1)%Series_R_RightEnd = 0.0_RKIND
         sgg%Med(contamedia)%SlantedWire(1)%Series_L_RightEnd = 0.0_RKIND
         sgg%Med(contamedia)%SlantedWire(1)%Series_C_RightEnd = 2.0e7_RKIND !en corto 14/2/14

         IF (this%swires%sw(j)%TL == Parallel_CONS) THEN
            sgg%Med(contamedia)%SlantedWire(1)%HasParallel_LeftEnd = .TRUE.
            sgg%Med(contamedia)%SlantedWire(1)%Parallel_R_LeftEnd = this%swires%sw(j)%R_LeftEnd
            sgg%Med(contamedia)%SlantedWire(1)%Parallel_L_LeftEnd = this%swires%sw(j)%L_LeftEnd
            sgg%Med(contamedia)%SlantedWire(1)%Parallel_C_LeftEnd = this%swires%sw(j)%C_LeftEnd
         ELSE IF (this%swires%sw(j)%TL == SERIES_CONS) THEN
            sgg%Med(contamedia)%SlantedWire(1)%HasSeries_LeftEnd = .TRUE.
            sgg%Med(contamedia)%SlantedWire(1)%Series_R_LeftEnd = this%swires%sw(j)%R_LeftEnd
            sgg%Med(contamedia)%SlantedWire(1)%Series_L_LeftEnd = this%swires%sw(j)%L_LeftEnd
            sgg%Med(contamedia)%SlantedWire(1)%Series_C_LeftEnd = this%swires%sw(j)%C_LeftEnd
         ELSE IF (this%swires%SW(j)%TL == DISPERSIVE_CONS) THEN
            allocate (sgg%Med(contamedia)%SlantedWire(1)%disp_LeftEnd(1))
            call asignawiredisper(sgg%Med(contamedia)%SlantedWire(1)%disp_LeftEnd(1), &
               this%swires%SW(j)%dispfile_LeftEnd)
         END IF
         IF (this%swires%sw(j)%TR == Parallel_CONS) THEN
            sgg%Med(contamedia)%SlantedWire(1)%HasParallel_RightEnd = .TRUE.
            sgg%Med(contamedia)%SlantedWire(1)%Parallel_R_RightEnd = this%swires%sw(j)%R_RightEnd
            sgg%Med(contamedia)%SlantedWire(1)%Parallel_L_RightEnd = this%swires%sw(j)%L_RightEnd
            sgg%Med(contamedia)%SlantedWire(1)%Parallel_C_RightEnd = this%swires%sw(j)%C_RightEnd
         ELSE IF (this%swires%sw(j)%TR == SERIES_CONS) THEN
            sgg%Med(contamedia)%SlantedWire(1)%HasSeries_RightEnd = .TRUE.
            sgg%Med(contamedia)%SlantedWire(1)%Series_R_RightEnd = this%swires%sw(j)%R_RightEnd
            sgg%Med(contamedia)%SlantedWire(1)%Series_L_RightEnd = this%swires%sw(j)%L_RightEnd
            sgg%Med(contamedia)%SlantedWire(1)%Series_C_RightEnd = this%swires%sw(j)%C_RightEnd
         ELSE IF (this%swires%SW(j)%TR == DISPERSIVE_CONS) THEN
            allocate (sgg%Med(contamedia)%SlantedWire(1)%disp_RightEnd(1))
            call asignawiredisper(sgg%Med(contamedia)%SlantedWire(1)%disp_RightEnd(1), &
               this%swires%SW(j)%dispfile_RightEnd)
         END IF

         sgg%Med(contamedia)%SlantedWire(1)%numNodes = this%swires%SW(j)%n_swc
         allocate(sgg%Med(contamedia)%SlantedWire(1)%nodes(1:this%swires%SW(j)%n_swc))
         do i = 1,this%swires%SW(j)%n_swc
            sgg%Med(contamedia)%SlantedWire(1)%nodes(i)%VsourceExists = .FALSE.
            sgg%Med(contamedia)%SlantedWire(1)%nodes(i)%IsourceExists = .FALSE.
            nullify(sgg%Med(contamedia)%SlantedWire(1)%nodes(i)%Vsource, &
               sgg%Med(contamedia)%SlantedWire(1)%nodes(i)%Isource)


            !!!2019 clipeo 2 celdas antes
            if (CLIPREGION) then
               if ((this%swires%SW(j)%swc(i)%x >= SINPML_Fullsize(iHx)%XE-2)) this%swires%SW(j)%swc(i)%x=SINPML_Fullsize(iHx)%XE-2
               if ((this%swires%SW(j)%swc(i)%x <= SINPML_Fullsize(iHx)%XI+2)) this%swires%SW(j)%swc(i)%x=SINPML_Fullsize(iHx)%XI+2
               if ((this%swires%SW(j)%swc(i)%y >= SINPML_Fullsize(iHy)%YE-2)) this%swires%SW(j)%swc(i)%y=SINPML_Fullsize(iHy)%YE-2
               if ((this%swires%SW(j)%swc(i)%y <= SINPML_Fullsize(iHy)%YI+2)) this%swires%SW(j)%swc(i)%y=SINPML_Fullsize(iHy)%YI+2
               if ((this%swires%SW(j)%swc(i)%z >= SINPML_Fullsize(iHz)%ZE-2)) this%swires%SW(j)%swc(i)%z=SINPML_Fullsize(iHz)%ZE-2
               if ((this%swires%SW(j)%swc(i)%z <= SINPML_Fullsize(iHz)%ZI+2)) this%swires%SW(j)%swc(i)%z=SINPML_Fullsize(iHz)%ZI+2
            endif

            !fin clipeo


            sgg%Med(contamedia)%SlantedWire(1)%nodes(i)%index = this%swires%SW(j)%swc(i)%nd
            sgg%Med(contamedia)%SlantedWire(1)%nodes(i)%x     = this%swires%SW(j)%swc(i)%x
            sgg%Med(contamedia)%SlantedWire(1)%nodes(i)%y     = this%swires%SW(j)%swc(i)%y
            sgg%Med(contamedia)%SlantedWire(1)%nodes(i)%z     = this%swires%SW(j)%swc(i)%z
            numertag = searchtag(tagtype,this%swires%SW(j)%swc(i)%tag)


            !!!2019 para informacion bbox hilos
            !!!!!!2020 retocado
            if ( int(this%swires%SW(j)%swc(i)%x)    < bboxwirxi) bboxwirXI=int(this%swires%SW(j)%swc(i)%x)
            if ( int(this%swires%SW(j)%swc(i)%x)+1  > bboxwirxE) bboxwirXE=int(this%swires%SW(j)%swc(i)%x)+1
            if ( int(this%swires%SW(j)%swc(i)%y)    < bboxwirYi) bboxwirYI=int(this%swires%SW(j)%swc(i)%y)
            if ( int(this%swires%SW(j)%swc(i)%y)+1  > bboxwirYE) bboxwirYE=int(this%swires%SW(j)%swc(i)%y)+1
            if ( int(this%swires%SW(j)%swc(i)%z)    < bboxwirZi) bboxwirZI=int(this%swires%SW(j)%swc(i)%z)
            if ( int(this%swires%SW(j)%swc(i)%z)+1  > bboxwirZE) bboxwirZE=int(this%swires%SW(j)%swc(i)%z)+1
            !!!!!!


            IF ((trim(adjustl(this%swires%SW(j)%swc(i)%SRCTYPE)) == 'VOLT').or. &
               (trim(adjustl(this%swires%SW(j)%swc(i)%SRCTYPE)) == 'SOFTVOLT')) THEN
               allocate (sgg%Med(contamedia)%SlantedWire(1)%nodes(i)%Vsource)
               sgg%Med(contamedia)%SlantedWire(1)%nodes(i)%VsourceExists = .TRUE.
               sgg%Med(contamedia)%SlantedWire(1)%nodes(i)%VSource%SOFT=.TRUE.  !fuentes duras sgg 230323. default blandas
               sgg%Med(contamedia)%SlantedWire(1)%nodes(i)%VSource%Multiplier = this%swires%SW(j)%swc(i)%m
               sgg%Med(contamedia)%SlantedWire(1)%nodes(i)%VSource%Resistance = 0.0_RKIND
               sgg%Med(contamedia)%SlantedWire(1)%nodes(i)%VSource%fichero%name = trim (adjustl(this%swires%SW(j)%swc(i)%SRCFILE))
            ELSE IF ((trim(adjustl(this%swires%SW(j)%swc(i)%SRCTYPE)) == 'CURR').OR. &
               (trim(adjustl(this%swires%SW(j)%swc(i)%SRCTYPE)) == 'HARDCURR')) THEN !EQUIVALENTES A IMPLEMENTAR UNA DE voltage
               allocate (sgg%Med(contamedia)%SlantedWire(1)%nodes(i)%Vsource)
               sgg%Med(contamedia)%SlantedWire(1)%nodes(i)%VsourceExists = .TRUE.
               sgg%Med(contamedia)%SlantedWire(1)%nodes(i)%Vsource%SOFT=.FALSE.
               sgg%Med(contamedia)%SlantedWire(1)%nodes(i)%Vsource%Multiplier = this%swires%SW(j)%swc(i)%m
               sgg%Med(contamedia)%SlantedWire(1)%nodes(i)%Vsource%Resistance = 0.0_RKIND
               sgg%Med(contamedia)%SlantedWire(1)%nodes(i)%Vsource%fichero%name = trim (adjustl(this%swires%SW(j)%swc(i)%SRCFILE))
            ELSEIF (trim(adjustl(this%swires%SW(j)%swc(i)%SRCTYPE)) == 'HARDVOLT') THEN !EQUIVALENTES A IMPLEMENTAR UNA DE Corriente
               allocate (sgg%Med(contamedia)%SlantedWire(1)%nodes(i)%Isource)
               sgg%Med(contamedia)%SlantedWire(1)%nodes(i)%IsourceExists = .TRUE.
               sgg%Med(contamedia)%SlantedWire(1)%nodes(i)%ISource%SOFT=.FALSE.
               sgg%Med(contamedia)%SlantedWire(1)%nodes(i)%ISource%Multiplier = this%swires%SW(j)%swc(i)%m
               sgg%Med(contamedia)%SlantedWire(1)%nodes(i)%ISource%Resistance = 0.0_RKIND
               sgg%Med(contamedia)%SlantedWire(1)%nodes(i)%ISource%fichero%name = trim (adjustl(this%swires%SW(j)%swc(i)%SRCFILE))
            ELSE IF ((trim(adjustl(this%swires%SW(j)%swc(i)%SRCTYPE)) == 'SOFTCURR')) THEN
               allocate (sgg%Med(contamedia)%SlantedWire(1)%nodes(i)%Isource)
               sgg%Med(contamedia)%SlantedWire(1)%nodes(i)%IsourceExists = .TRUE.
               sgg%Med(contamedia)%SlantedWire(1)%nodes(i)%Isource%SOFT=.TRUE.
               sgg%Med(contamedia)%SlantedWire(1)%nodes(i)%Isource%Multiplier = this%swires%SW(j)%swc(i)%m
               sgg%Med(contamedia)%SlantedWire(1)%nodes(i)%Isource%Resistance = 0.0_RKIND
               sgg%Med(contamedia)%SlantedWire(1)%nodes(i)%Isource%fichero%name = trim (adjustl(this%swires%SW(j)%swc(i)%SRCFILE))
            elseIF (trim(adjustl(this%swires%SW(j)%swc(i)%SRCTYPE)) /= 'None') THEN
               write(buff,*) 'WRONG type of wire source '//trim(adjustl(this%swires%SW(j)%swc(i)%SRCTYPE))
               call stoponerror (layoutnumber,size,buff)
            END IF
         end do
      end do
      !END SLANTED WIRES

      call addConformalMedia(sgg, media, conformal_media, contamedia)



      !reporta el bounding box

#ifdef CompileWithMPI
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_AllReduce( bboxwirXI, dummy_bboxwirXI, 1_4, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr)
      call MPI_AllReduce( bboxwirYI, dummy_bboxwirYI, 1_4, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr)
      call MPI_AllReduce( bboxwirZI, dummy_bboxwirzI, 1_4, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr)
      call MPI_AllReduce( bboxwirXE, dummy_bboxwirXE, 1_4, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
      call MPI_AllReduce( bboxwirYE, dummy_bboxwirYE, 1_4, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
      call MPI_AllReduce( bboxwirZE, dummy_bboxwirZE, 1_4, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      bboxwirXI=dummy_bboxwirXI
      bboxwirYI=dummy_bboxwirYI
      bboxwirzI=dummy_bboxwirzI
      bboxwirXE=dummy_bboxwirXE
      bboxwirYE=dummy_bboxwirYE
      bboxwirZE=dummy_bboxwirZE
#endif
      WRITE (buff, '(a, 6I12)')    'pre1_INFO:  Bounding Box for WIREs min_x,min_y,min_z, MAX_x,MAX_y,MAX_z',bboxwirXI,bboxwirYI,bboxwirZI,bboxwirXE,bboxwirYE,bboxwirZE
      if (((bboxwirXI<2**20).or.(bboxwirYI<2**20).or.(bboxwirZI<2**20).or.(bboxwirXE>-(2**20)).or.(bboxwirYE>-(2**20)).or.(bboxwirZE>-(2**20))).or.(VERBOSE)) then
         CALL WarnErrReport (buff)
      ENDIF
      !FIN WIRES

      if (run_with_dmma) then
         !always at the end since the orientation is found from the PEC one
         !thin Slots
         !the embedding material properties are also if also needed
         tama = this%tSlots%n_tg
         !
         DO j = 1, tama
            !
            !
            tama2 = this%tSlots%Tg(j)%N_tgc
            DO i = 1, tama2
               !del Slot
               direccion = this%tSlots%Tg(j)%TgC(i)%dir
               i1 = this%tSlots%Tg(j)%TgC(i)%i
               j1 = this%tSlots%Tg(j)%TgC(i)%j
               k1 = this%tSlots%Tg(j)%TgC(i)%k
               ORIX=.FALSE. ; ORIY=.FALSE. ; ORIZ=.FALSE. ;
               ORIX2=.FALSE. ; ORIY2=.FALSE. ; ORIZ2=.FALSE. ;
               ORIX3=.FALSE. ; ORIY3=.FALSE. ; ORIZ3=.FALSE. ;
               ORIX4=.FALSE. ; ORIY4=.FALSE. ; ORIZ4=.FALSE. ;
               !
               IF ((i1 >= BoundingBox%XI) .AND. (i1 < BoundingBox%XE) .AND. &
               &    (j1 >= BoundingBox%YI) .AND. (j1 < BoundingBox%YE) .AND. &
               &    (k1 >= BoundingBox%ZI) .AND. (k1 < BoundingBox%ZE)) THEN
                  !encuentra la orientacion del plano PEC que contiene al Slot
                  oriX = (direccion == iEy)  .AND.   &
                  &       (((media%sggMiHx(i1, j1, k1) ==0).or.(sgg%med(media%sggMiHx(i1, j1, k1) )%is%pec)) .OR. &
                  &       (sgg%Med(media%sggMiHx(i1, j1, k1))%Is%ThinSlot))          !&
                  !& .AND. (((sggmiHz(i1, j1, k1) /=0).and.(.not.(sgg%med(sggmiHz(i1, j1, k1) )%is%pec))) .AND.                   &
                  !&       ( .NOT. sgg%Med(sggmiHz(i1, j1,k1))%Is%ThinSlot))

                  oriX4 =(direccion == iEz) .AND.    &
                  &       (((               media%sggMiHx(i1, j1, k1) ==0).or.(sgg%med(               media%sggMiHx(i1, j1, k1) )%is%pec)) .OR.       &
                  &        (       sgg%Med(media%sggMiHx(i1, j1, k1))%Is%ThinSlot))   !&
                  !& .AND. (((               sggmiHy(i1, j1, k1) /=0).and.(.not.(sgg%med(               sggmiHy(i1, j1, k1) )%is%pec))) .AND.      &
                  !&       ( .NOT.  sgg%Med(sggmiHy(i1, j1, k1))%Is%ThinSlot))

                  oriY = (direccion == iEx) .AND.   &
                  &       (((               media%sggMiHy(i1, j1, k1) ==0).or.(sgg%med(               media%sggMiHy(i1, j1, k1) )%is%pec)) .OR.       &
                  &        (       sgg%Med(media%sggMiHy(i1, j1, k1))%Is%ThinSlot))   !&
                  !& .AND. (((               sggmiHz(i1, j1, k1) /=0).and.(.not.(sgg%med(               sggmiHz(i1, j1, k1) )%is%pec))) .AND.      &
                  !&        ( .NOT. sgg%Med(sggmiHz(i1, j1, k1))%Is%ThinSlot))

                  oriY4 =(direccion == iEz) .AND.   &
                  &       (((              media%sggMiHy(i1, j1, k1) ==0).or.(sgg%med(              media%sggMiHy(i1, j1, k1) )%is%pec)) .OR.        &
                  &        (      sgg%Med(media%sggMiHy(i1, j1, k1))%Is%ThinSlot))    !&
                  !&  .AND.(((               sggmiHx(i1, j1, k1) /=0).and.(.not.(sgg%med(              sggmiHx(i1, j1, k1) )%is%pec))) .AND.       &
                  !&       ( .NOT. sgg%Med(sggmiHx(i1, j1, k1))%Is%ThinSlot))

                  oriZ = (direccion == iEx) .AND.  &
                  &       (((               media%sggMiHz(i1, j1, k1) ==0).or.(sgg%med(               media%sggMiHz(i1, j1, k1) )%is%pec)) .OR.        &
                  &        (       sgg%Med(media%sggMiHz(i1, j1, k1))%Is%ThinSlot))    !&
                  !&       (((               sggmiHy(i1, j1, k1) /=0).and.(.not.(sgg%med(               sggmiHy(i1, j1, k1) )%is%pec))) .AND.       &
                  !&        (.NOT.  sgg%Med(sggmiHy(i1, j1, k1))%Is%ThinSlot))

                  oriZ4 = (direccion == iEy) .AND.   &
                  &        (((               media%sggMiHz(i1, j1, k1) ==0).or.(sgg%med(               media%sggMiHz(i1, j1, k1) )%is%pec)) .OR.       &
                  &         (       sgg%Med(media%sggMiHz(i1, j1, k1))%Is%ThinSlot))   !&
                  !& .AND.  (((              sggmiHx(i1, j1, k1) /=0).and.(.not.(sgg%med(               sggmiHx(i1, j1, k1) )%is%pec))) .AND.      &
                  !&         ( .NOT. sgg%Med(sggmiHx(i1, j1, k1))%Is%ThinSlot))

                  !encuentra la orientacion del plano PEC que contiene al Slot (considera los vecinos)
                  oriX2 = (direccion == iEy) .AND.   &
                  &        (((               media%sggMiHx(i1, j1, k1-1) ==0).or.(sgg%med(               media%sggMiHx(i1, j1, k1-1) )%is%pec)) .OR.     &
                  &         (       sgg%Med(media%sggMiHx(i1, j1, k1-1))%Is%ThinSlot)) !&
                  !& .AND.  (((              sggmiHz(i1, j1, k1-1) /=0).and.(.not.(sgg%med(               sggmiHz(i1, j1, k1-1) )%is%pec))) .AND.    &
                  !&         ( .NOT. sgg%Med(sggmiHz(i1, j1, k1-1))%Is%ThinSlot))

                  oriX3 = (direccion == iEz) .AND.   &
                  &        (((                media%sggMiHx(i1, j1-1, k1) ==0).or.(sgg%med(                media%sggMiHx(i1, j1-1, k1) )%is%pec)) .OR.     &
                  &         (        sgg%Med(media%sggMiHx(i1, j1-1, k1))%Is%ThinSlot)) !&
                  !& .AND.  (((              sggmiHy(i1, j1-1, k1) /=0).and.(.not.(sgg%med(                sggmiHy(i1, j1-1, k1) )%is%pec))) .AND.    &
                  !&         ( .NOT.  sgg%Med(sggmiHy(i1, j1-1, k1))%Is%ThinSlot))

                  oriY2 = (direccion == iEx) .AND.   &
                  &        (((               media%sggMiHy(i1, j1, k1-1) ==0).or.(sgg%med(               media%sggMiHy(i1, j1, k1-1) )%is%pec)) .OR.      &
                  &         (sgg%Med(       media%sggMiHy(i1, j1, k1-1))%Is%ThinSlot))  !&
                  !& .AND.  (((              sggmiHz(i1, j1, k1-1) /=0).and.(.not.(sgg%med(               sggmiHz(i1, j1, k1-1) )%is%pec))) .AND.     &
                  !&         ( .NOT. sgg%Med(sggmiHz(i1, j1, k1-1))%Is%ThinSlot))

                  oriY3 = (direccion == iEz) .AND.   &
                  &        (((               media%sggMiHy(i1-1, j1, k1) ==0).or.(sgg%med(               media%sggMiHy(i1-1, j1, k1) )%is%pec)) .OR.      &
                  &         (       sgg%Med(media%sggMiHy(i1-1, j1, k1))%Is%ThinSlot))  !&
                  !& .AND.  (((              sggmiHx(i1-1, j1, k1) /=0).and.(.not.(sgg%med(               sggmiHx(i1-1, j1, k1) )%is%pec))) .AND.     &
                  !&         ( .NOT. sgg%Med(sggmiHx(i1-1, j1, k1))%Is%ThinSlot))

                  oriZ2 = (direccion == iEx)  .AND.   &
                  &        (((               media%sggMiHz(i1, j1-1, k1) ==0).or.(sgg%med(               media%sggMiHz(i1, j1-1, k1) )%is%pec)) .OR.      &
                  &         (       sgg%Med(media%sggMiHz(i1, j1-1, k1))%Is%ThinSlot))  !&
                  !& .AND.  (((              sggmiHy(i1, j1-1, k1) /=0).and.(.not.(sgg%med(               sggmiHy(i1, j1-1, k1) )%is%pec))) .AND.     &
                  !&         ( .NOT. sgg%Med(sggmiHy(i1, j1-1, k1))%Is%ThinSlot))


                  oriZ3 = (direccion == iEy) .AND.   &
                  &        (((               media%sggMiHz(i1-1, j1, k1) ==0).or.(sgg%med(               media%sggMiHz(i1-1, j1, k1) )%is%pec)) .OR.      &
                  &         (       sgg%Med(media%sggMiHz(i1-1, j1, k1))%Is%ThinSlot))  !&
                  !& .AND.  (((              sggmiHx(i1-1, j1, k1) /=0).and.(.not.(sgg%med(               sggmiHx(i1-1, j1, k1) )%is%pec))) .AND.     &
                  !&         ( .NOT. sgg%Med(sggmiHx(i1-1, j1, k1))%Is%ThinSlot))

                  IF (oriX.or.oriX4) THEN
                     orientacion = iEx
                  ELSE IF (oriY.or.oriY4) THEN
                     orientacion = iEy
                  ELSE IF (oriZ.or.oriZ4) THEN
                     orientacion = iEz
                     !vecinos
                  else IF (oriX2) THEN
                     orientacion = iEx
                     k1 = k1-1
                  ELSE IF (oriY2) THEN
                     orientacion = iEy
                     k1 = k1-1
                  ELSE IF (oriZ2) THEN
                     orientacion = iEz
                     j1 = j1-1
                     !vecinos
                  else IF (oriX3) THEN
                     orientacion = iEx
                     j1 = j1-1
                  ELSE IF (oriY3) THEN
                     orientacion = iEy
                     i1 = i1-1
                  ELSE IF (oriZ3) THEN
                     orientacion = iEz
                     i1 = i1-1
                  ELSE
                     write(buff,*) 'Cannot determine ortientation of the PEC plane with the Slot',i1, j1, k1, direccion
                     call stoponerror (layoutnumber,size,buff)
                     !ojo con el nfde no se puede hacer Slots en escalera porque no se puede determinar la orientacion de los planos
                     !en los tramos comunes. Por tanto No he podido testear los shared electricos anisotropos. solo los magneticos
                  END IF

                  this%tSlots%Tg(j)%TgC(i)%or = orientacion
                  medio2=-1
                  medio1=-1
                  SELECT CASE (Abs(orientacion))
                   CASE (iEx)
                     medio1 = media%sggMiHx(i1,j1,k1) !!!sggmcen (i1, j1, k1) !tocaco 03/07/15 para lo eliminar lo de los media matrix !puede que me haya cargado los thin-slots en materialescon esto 03/07/15
                     IF (i1 > BoundingBox%XI) then
                        medio2 = media%sggMiHx(i1-1,j1,k1)  !!!sggmcen (i1-1, j1, k1) !tocaco 03/07/15 para lo eliminar lo de los media matrix
                     else
                        medio2=medio1
                     endif
                   CASE (iEy)
                     medio1 = media%sggMiHy(i1,j1,k1) !!!sggmcen (i1, j1, k1) !tocaco 03/07/15 para lo eliminar lo de los media matrix
                     IF (j1 > BoundingBox%YI) then
                        medio2 = media%sggMiHy(i1,j1-1,k1) !!!sggmcen (i1, j1-1, k1) !tocaco 03/07/15 para lo eliminar lo de los media matrix
                     else
                        medio2=medio1
                     endif
                   CASE (iEz)
                     medio1 = media%sggMiHz(i1,j1,k1) !!!sggmcen (i1, j1, k1) !tocaco 03/07/15 para lo eliminar lo de los media matrix
                     IF (k1 > BoundingBox%ZI) then
                        medio2 = media%sggMiHz(i1,j1,k1-1) !!! sggmcen (i1, j1, k1-1) !tocaco 03/07/15 para lo eliminar lo de los media matrix
                     else
                        medio2=medio1
                     endif
                  END SELECT

                  IF ( ((sgg%Med(medio1)%Is%Dielectric).or.(sgg%Med(medio1)%Is%Thinslot).or.(sgg%Med(medio1)%Is%Pec).or.(medio1 ==1 )).and. &
                     ((sgg%Med(medio2)%Is%Dielectric).or.(sgg%Med(medio2)%Is%Thinslot).or.(sgg%Med(medio2)%Is%Pec).or.(medio2 ==1 )) ) THEN
                     !average adjacent media
                     !
                     epr1 = 0.5_RKIND  * (sgg%Med(medio1)%Epr+sgg%Med(medio2)%Epr)
                     mur1 = 0.5_RKIND  * (sgg%Med(medio1)%Mur+sgg%Med(medio2)%Mur)
                  ELSE
                     write(buff,*) 'Media around the Slot are not plain media: ', medio1, medio2
                     CALL STOPONERROR(layoutnumber,size,buff)
                  END IF
                  width = this%tSlots%Tg(j)%width
                  IF (sgg%NumPlaneWaves == 1) THEN
                     !ojo no me gusta pq no es general
                     !assume the incident plane wave if there are planewaves
                     !
                     dir (1) = px
                     dir (2) = py
                     dir (3) = pz
                  ELSE
                     !assume normal incidence
                     !
                     SELECT CASE (Abs(orientacion))
                      CASE (iEx)
                        dir (1) = 1.0_RKIND
                        dir (2) = 0.0_RKIND
                        dir (3) = 0.0_RKIND
                      CASE (iEy)
                        dir (1) = 0.0_RKIND
                        dir (2) = 1.0_RKIND
                        dir (3) = 0.0_RKIND
                      CASE (iEz)
                        dir (1) = 0.0_RKIND
                        dir (2) = 0.0_RKIND
                        dir (3) = 1.0_RKIND
                     END SELECT
                  END IF
                  CALL dmma_thin_Slot (sgg%dx(i1), sgg%dy(j1), sgg%dz(k1), dir,   &
                  &        orientacion, direccion, width, epr1, mur1, EprSlot, MurSlot,eps0,mu0)
                  ! y tocar el precounting if so
                  !chequear que son distintos para incrementar contamedia
                  !
                  indicemedio = contamedia + 1
                  buscaiguales: DO ii = 1, contamedia
                     IF (sgg%Med(ii)%Is%ThinSlot) THEN
                        iguales = .TRUE.
                        DO j11 = 1, 3
                           DO i11 = 1, 3
                              iguales = iguales .AND. (sgg%Med(ii)%Anisotropic(1)%Epr(i11, j11) == EprSlot(i11, j11)) .AND. &
                              & (sgg%Med(ii)%Anisotropic(1)%Mur(i11, j11) == MurSlot(i11, j11))
                           END DO
                        END DO
                        IF (iguales) THEN
                           indicemedio = ii
                           EXIT buscaiguales
                        END IF
                     END IF
                  END DO buscaiguales
                  IF (indicemedio == contamedia+1) THEN
                     contamedia = indicemedio
                     ALLOCATE (sgg%Med(contamedia)%Anisotropic(1))
                     sgg%Med(contamedia)%Anisotropic(1)%Epr = EprSlot
                     sgg%Med(contamedia)%Anisotropic(1)%Mur = MurSlot
                     !lossless
                     sgg%Med(contamedia)%Anisotropic(1)%Sigma = 0.0_RKIND
                     sgg%Med(contamedia)%Anisotropic(1)%SigmaM = 0.0_RKIND
                     !
                     sgg%Med(contamedia)%Epr = epr1
                     sgg%Med(contamedia)%Mur = mur1
                     sgg%Med(contamedia)%Sigma = 0.0_RKIND
                     sgg%Med(contamedia)%SigmaM = 0.0_RKIND
                     sgg%Med(contamedia)%Is%Anisotropic = .TRUE.
                     !just for signaling
                     sgg%Med(contamedia)%Is%ThinSlot = .TRUE.
                     sgg%Med(contamedia)%Priority = prior_TG
                  END IF
                  !
                  !
                  !record coordinates
                  !
                  punto%XI = i1
                  punto%XE = i1
                  punto%YI = j1
                  punto%YE = j1
                  punto%ZI = k1
                  punto%ZE = k1
                  numertag = searchtag(tagtype,this%tSlots%Tg(j)%TgC(i)%tag)
                  CALL CreateSurfaceSlotMM (layoutnumber, media%sggMtag, tag_numbers, numertag, media%sggMiEx, media%sggMiEy, media%sggMiEz,&
                     media%sggMiHx, media%sggMiHy, media%sggMiHz, Alloc_iEx_XI,&
                     Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, &
                     Alloc_iEy_YI,&
                     Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, &
                     Alloc_iEz_ZI,&
                     Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, &
                     Alloc_iHy_XI,&
                     Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, &
                     Alloc_iHz_YI,&
                     Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg%Med, sgg%NumMedia, sgg%EShared, sgg%HShared, BoundingBox, &
                     punto, orientacion, direccion, indicemedio)
                  !del if esta dentro del bounding box
               END IF
            END DO
            !thin Slots
         END DO
         !
      endif !del run_with_dmma

      !debe ir al final para respetar el tipo de medio que haya SI SE TRATASE COMO A UN MEDIO
      !nodalsource
      !precounting
      tama = this%nodsrc%n_nodSrc
      !at most
      ALLOCATE (contapuntos(tama*(this%nodsrc%n_C2p_max+this%nodsrc%n_C1p_max)))
      contapuntos = 0
      conta1 = 0
      DO i = 1, tama
         conta2 = 0
         tama2 = this%nodsrc%NodalSource(i)%n_c1P
         tama3 = this%nodsrc%NodalSource(i)%n_c2P
         DO ii = 1, tama2
            punto_s%or = this%nodsrc%NodalSource(i)%c1P(ii)%or
            punto_s%XI = this%nodsrc%NodalSource(i)%c1P(ii)%XI
            punto_s%XE = this%nodsrc%NodalSource(i)%c1P(ii)%XE
            punto_s%YI = this%nodsrc%NodalSource(i)%c1P(ii)%YI
            punto_s%YE = this%nodsrc%NodalSource(i)%c1P(ii)%YE
            punto_s%ZI = this%nodsrc%NodalSource(i)%c1P(ii)%ZI
            punto_s%ZE = this%nodsrc%NodalSource(i)%c1P(ii)%ZE
            IF ((punto_s%XI <= punto_s%XE) .AND. (punto_s%YI <= punto_s%YE) .AND. (punto_s%ZI <= punto_s%ZE)) THEN
               conta2 = conta2 + 1
            END IF
         END DO
         !
         !
         DO ii = 1, tama3
            punto_s%or = this%nodsrc%NodalSource(i)%c2P(ii)%or
            punto_s%XI = this%nodsrc%NodalSource(i)%c2P(ii)%XI
            punto_s%XE = this%nodsrc%NodalSource(i)%c2P(ii)%XE
            punto_s%YI = this%nodsrc%NodalSource(i)%c2P(ii)%YI
            punto_s%YE = this%nodsrc%NodalSource(i)%c2P(ii)%YE
            punto_s%ZI = this%nodsrc%NodalSource(i)%c2P(ii)%ZI
            punto_s%ZE = this%nodsrc%NodalSource(i)%c2P(ii)%ZE
            IF ((punto_s%XI <= punto_s%XE) .AND. (punto_s%YI <= punto_s%YE) .AND. (punto_s%ZI <= punto_s%ZE)) THEN
               conta2 = conta2 + 1
            END IF
         END DO
         IF (conta2 /= 0) THEN
            conta1 = conta1 + 1
            contapuntos (conta1) = conta2
         END IF
      END DO
      sgg%NumNodalSources = conta1
      ALLOCATE (sgg%NodalSource(conta1))
      if (sgg%NumNodalSources /=0) contamedia =contamedia+1
      !
      conta1 = 0
      DO i = 1, tama
         IF (contapuntos(i) /= 0) THEN
            conta1 = conta1 + 1
            sgg%NodalSource(conta1)%numpuntos = contapuntos (conta1)
            ALLOCATE (sgg%NodalSource(conta1)%punto(contapuntos(conta1)))
            !initialization
            do ii=1,contapuntos(conta1)
               sgg%NodalSource(conta1)%punto(ii)%or = 0
               sgg%NodalSource(conta1)%punto(ii)%xc = 0.0_RKIND
               sgg%NodalSource(conta1)%punto(ii)%yc = 0.0_RKIND
               sgg%NodalSource(conta1)%punto(ii)%zc = 0.0_RKIND
               sgg%NodalSource(conta1)%punto(ii)%XI = -1
               sgg%NodalSource(conta1)%punto(ii)%XE = -1
               sgg%NodalSource(conta1)%punto(ii)%YI = -1
               sgg%NodalSource(conta1)%punto(ii)%YE = -1
               sgg%NodalSource(conta1)%punto(ii)%ZI = -1
               sgg%NodalSource(conta1)%punto(ii)%ZE = -1
            end do
         endif
      END DO
      !asignacion
      conta1 = 0
      DO i = 1, tama
         conta2 = 0
         IF (contapuntos(i) /= 0) then
            conta1 = conta1 + 1
            !
            sgg%NodalSource(conta1)%fichero%name = trim (adjustl(this%nodsrc%NodalSource(i)%nombre))
            sgg%NodalSource(conta1)%isElec = this%nodsrc%NodalSource(i)%isElec
            sgg%NodalSource(conta1)%IsHard = this%nodsrc%NodalSource(i)%isHard
            sgg%NodalSource(conta1)%IsInitialValue = this%nodsrc%NodalSource(i)%IsInitialValue 
         endif
         !
         tama2 = this%nodsrc%NodalSource(i)%n_c1P
         tama3 = this%nodsrc%NodalSource(i)%n_c2P
         DO ii = 1, tama2
            !!correct bounding box
            punto_s%or = this%nodsrc%NodalSource(i)%c1P(ii)%or
            punto_s%xc = this%nodsrc%NodalSource(i)%c1P(ii)%xc
            punto_s%yc = this%nodsrc%NodalSource(i)%c1P(ii)%yc
            punto_s%zc = this%nodsrc%NodalSource(i)%c1P(ii)%zc
            !
            punto_s%XI = Max (this%nodsrc%NodalSource(i)%c1P(ii)%XI, Min(BoundingBox%XI, BoundingBox%XE))
            punto_s%YI = Max (this%nodsrc%NodalSource(i)%c1P(ii)%YI, Min(BoundingBox%YI, BoundingBox%YE))
            punto_s%ZI = Max (this%nodsrc%NodalSource(i)%c1P(ii)%ZI, Min(BoundingBox%ZI, BoundingBox%ZE))
            !
            punto_s%XE = Min (this%nodsrc%NodalSource(i)%c1P(ii)%XE, Max(BoundingBox%XI, BoundingBox%XE))
            punto_s%YE = Min (this%nodsrc%NodalSource(i)%c1P(ii)%YE, Max(BoundingBox%YI, BoundingBox%YE))
            IF ((punto_s%zc /= 0).and.(this%nodsrc%NodalSource(i)%isElec))  then !only in case of Ez
               punto_s%ZE = Min (this%nodsrc%NodalSource(i)%c1P(ii)%ZE, Max(BoundingBox%ZI, BoundingBox%ZE-1))
            else
               punto_s%ZE = Min (this%nodsrc%NodalSource(i)%c1P(ii)%ZE, Max(BoundingBox%ZI, BoundingBox%ZE  ))
            endif
            !
            !
            DO k1 = punto_s%ZI, punto_s%ZE
               DO j1 = punto_s%YI, punto_s%YE
                  DO i1 = punto_s%XI, punto_s%XE
                     IF (punto_s%xc /= 0) THEN
                        !bug OLD 181214 sl_4_20mm_gli.nfde. Fuente nodal electrica embebida en pec y nodal magnetica en pmc se ignoraran sean hard or soft
                        MEDIO = media%sggMiEx (i1, j1, k1)
                        valido=.true.
                        !IF ( .NOT. this%nodsrc%NodalSource(i)%isHard) THEN
                        !  VALIDO = (sgg%Med(MEDIO)%Is%Dielectric) .OR. (sgg%Med(MEDIO)%Is%EDispersive) .OR. &
                        ! & (sgg%Med(MEDIO)%Is%MDispersive)
                        !ELSE
                        !  VALIDO = .TRUE.
                        !END IF
                        IF (this%nodsrc%NodalSource(i)%isElec) THEN
                           VALIDO = VALIDO .AND. ( .NOT. sgg%Med(MEDIO)%Is%PEC)
                        END IF
                        write (buff,*) 'WARNING: Ex Nodal source on PEC media will be ignored (', i1, j1, k1,')'
                        IF ( .NOT. VALIDO) CALL  WarnErrReport (buff)
                        !
                        ! COMENTADO 250816 PQ DA UN ERROR JUISTO CUANDO CAE LA FUENTE EN UN CORTE MPI. HABRIA QU TOCA LA CASUISTICA DE punto_s%ZE = Min (this%nodsrc%NodalSource(i)%c1P(ii)%ZE, Max(BoundingBox%ZI, BoundingBox%ZE  )) PERO NO LO HE QUERIDO HACER
                        !!!MEDIO = sggmiHx (i1, j1, k1)
                        !!!valido=.true.
                        !!!!IF ( .NOT. this%nodsrc%NodalSource(i)%isHard) THEN
                        !!!!  VALIDO = (sgg%Med(MEDIO)%Is%Dielectric) .OR. (sgg%Med(MEDIO)%Is%EDispersive) .OR. &
                        !!!! & (sgg%Med(MEDIO)%Is%MDispersive)
                        !!!!ELSE
                        !!!!  VALIDO = .TRUE.
                        !!!!END IF
                        !!!IF ( .NOT. this%nodsrc%NodalSource(i)%isElec) THEN
                        !!!   VALIDO = VALIDO .AND. ( .NOT. sgg%Med(MEDIO)%Is%PMC)
                        !!!END IF
                        !!!write (buff,*) 'WARNING: Hx Nodal source on PMC media will be ignored (', i1, j1, k1,')'
                        !!!IF ( .NOT. VALIDO) CALL  WarnErrReport (buff)
                     END IF
                     !
                     !
                     IF (punto_s%yc /= 0) THEN
                        MEDIO = media%sggMiEy (i1, j1, k1)
                        valido=.true.
                        !IF ( .NOT. this%nodsrc%NodalSource(i)%isHard) THEN
                        !  VALIDO = (sgg%Med(MEDIO)%Is%Dielectric) .OR. (sgg%Med(MEDIO)%Is%EDispersive) .OR. &
                        ! & (sgg%Med(MEDIO)%Is%MDispersive)
                        !ELSE
                        !  VALIDO = .TRUE.
                        !END IF
                        IF (this%nodsrc%NodalSource(i)%isElec) THEN
                           VALIDO = VALIDO .AND. ( .NOT. sgg%Med(MEDIO)%Is%PEC)
                        END IF
                        write (buff,*) 'WARNING: Ey Nodal source on PMC media will be ignored (', i1, j1, k1,')'
                        IF ( .NOT. VALIDO) CALL  WarnErrReport (buff)
                        !
                        !
                        !!!MEDIO = sggmiHy (i1, j1, k1)
                        !!!valido=.true.
                        !!!!IF ( .NOT. this%nodsrc%NodalSource(i)%isHard) THEN
                        !!!!  VALIDO = (sgg%Med(MEDIO)%Is%Dielectric) .OR. (sgg%Med(MEDIO)%Is%EDispersive) .OR. &
                        !!!! & (sgg%Med(MEDIO)%Is%MDispersive)
                        !!!!ELSE
                        !!!!  VALIDO = .TRUE.
                        !!!!END IF
                        !!!IF ( .NOT. this%nodsrc%NodalSource(i)%isElec) THEN
                        !!!   VALIDO = VALIDO .AND. ( .NOT. sgg%Med(MEDIO)%Is%PMC)
                        !!!END IF
                        !!!write (buff,*) 'WARNING: Hy Nodal source on PMC media will be ignored (', i1, j1, k1,')'
                        !!!IF ( .NOT. VALIDO) CALL  WarnErrReport (buff)
                     END IF
                     !
                     !
                     IF (punto_s%zc /= 0) THEN
                        MEDIO = media%sggMiEz (i1, j1, k1)
                        valido=.true.
                        !IF ( .NOT. this%nodsrc%NodalSource(i)%isHard) THEN
                        !  VALIDO = (sgg%Med(MEDIO)%Is%Dielectric) .OR. (sgg%Med(MEDIO)%Is%EDispersive) .OR. &
                        ! & (sgg%Med(MEDIO)%Is%MDispersive)
                        !ELSE
                        !  VALIDO = .TRUE.
                        !END IF
                        IF (this%nodsrc%NodalSource(i)%isElec) THEN
                           VALIDO = VALIDO .AND. ( .NOT. sgg%Med(MEDIO)%Is%PEC)
                        END IF
                        write (buff,*) 'WARNING: Ez Nodal source on PMC media will be ignored (', i1, j1, k1,')'
                        IF ( .NOT. VALIDO) CALL  WarnErrReport (buff)
                        !
                        !
                        !!!MEDIO = sggmiHz (i1, j1, k1)
                        !!!valido=.true.
                        !!!!IF ( .NOT. this%nodsrc%NodalSource(i)%isHard) THEN
                        !!!!  VALIDO = (sgg%Med(MEDIO)%Is%Dielectric) .OR. (sgg%Med(MEDIO)%Is%EDispersive) .OR. &
                        !!!! & (sgg%Med(MEDIO)%Is%MDispersive)
                        !!!!ELSE
                        !!!!  VALIDO = .TRUE.
                        !!!!END IF
                        !!!IF ( .NOT. this%nodsrc%NodalSource(i)%isElec) THEN
                        !!!   VALIDO = VALIDO .AND. ( .NOT. sgg%Med(MEDIO)%Is%PMC)
                        !!!END IF
                        !!!write (buff,*) 'WARNING: Hz Nodal source on PMC media will be ignored (', i1, j1, k1,')'
                        !!!IF ( .NOT. VALIDO) CALL  WarnErrReport (buff)
                     END IF
                  END DO
               END DO
            END DO
            !
            !
            IF ((punto_s%XI <= punto_s%XE) .AND. (punto_s%YI <= punto_s%YE) .AND. (punto_s%ZI <= punto_s%ZE)) THEN
               conta2 = conta2 + 1
               sgg%NodalSource(conta1)%punto(conta2)%or = punto_s%or
               sgg%NodalSource(conta1)%punto(conta2)%xc = punto_s%xc
               sgg%NodalSource(conta1)%punto(conta2)%yc = punto_s%yc
               sgg%NodalSource(conta1)%punto(conta2)%zc = punto_s%zc
               sgg%NodalSource(conta1)%punto(conta2)%XI = punto_s%XI
               sgg%NodalSource(conta1)%punto(conta2)%XE = punto_s%XE
               sgg%NodalSource(conta1)%punto(conta2)%YI = punto_s%YI
               sgg%NodalSource(conta1)%punto(conta2)%YE = punto_s%YE
               sgg%NodalSource(conta1)%punto(conta2)%ZI = punto_s%ZI
               sgg%NodalSource(conta1)%punto(conta2)%ZE = punto_s%ZE
               !PARA ACOMODAR LAS NODAL SOURCE COMO MEDIOS LINE Y PODER VISUALIZAR SONDAS 010824
               sgg%Med(contamedia)%Is%Dielectric = .TRUE.
               sgg%Med(contamedia)%Is%LINE = .TRUE.
               sgg%Med(contamedia)%Priority = prior_IL
               sgg%Med(contamedia)%Epr =  1.0
               sgg%Med(contamedia)%Sigma = 0.
               sgg%Med(contamedia)%Mur =  1.0
               sgg%Med(contamedia)%SigmaM = 0.
               punto%XI = punto_s%XI
               punto%XE = punto_s%XE
               punto%YI = punto_s%YI
               punto%YE = punto_s%YE
               punto%ZI = punto_s%ZI
               punto%ZE = punto_s%ZE
               orientacion = punto_s%or
               isathinwire = .FALSE.
               numertag = 37
               CALL CreateLineMM (layoutnumber, media%sggMtag, tag_numbers, numertag, media%sggMiEx, media%sggMiEy, media%sggMiEz, &
                  media%sggMiHx, media%sggMiHy, media%sggMiHz, Alloc_iEx_XI, &
                  Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI, &
                  Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI, &
                  Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI, &
                  Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI, &
                  Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg%Med, sgg%NumMedia, sgg%EShared, BoundingBox, punto, orientacion, &
                  contamedia, isathinwire,verbose,numeroasignaciones)
            END IF
            sgg%NodalSource(conta1)%numpuntos = conta2 !update with the correct value
         END DO
         !
         !
         DO ii = 1, tama3
            punto_s%or = this%nodsrc%NodalSource(i)%c2P(ii)%or
            punto_s%XI = this%nodsrc%NodalSource(i)%c2P(ii)%XI
            punto_s%XE = this%nodsrc%NodalSource(i)%c2P(ii)%XE
            punto_s%YI = this%nodsrc%NodalSource(i)%c2P(ii)%YI
            punto_s%YE = this%nodsrc%NodalSource(i)%c2P(ii)%YE
            punto_s%ZI = this%nodsrc%NodalSource(i)%c2P(ii)%ZI
            punto_s%ZE = this%nodsrc%NodalSource(i)%c2P(ii)%ZE
            punto_s%xc = this%nodsrc%NodalSource(i)%c2P(ii)%xc
            punto_s%yc = this%nodsrc%NodalSource(i)%c2P(ii)%yc
            punto_s%zc = this%nodsrc%NodalSource(i)%c2P(ii)%zc
            !!correct bounding box
            punto_s%XI = Max (this%nodsrc%NodalSource(i)%c2p(ii)%XI, Min(BoundingBox%XI, BoundingBox%XE))
            punto_s%YI = Max (this%nodsrc%NodalSource(i)%c2p(ii)%YI, Min(BoundingBox%YI, BoundingBox%YE))
            punto_s%ZI = Max (this%nodsrc%NodalSource(i)%c2p(ii)%ZI, Min(BoundingBox%ZI, BoundingBox%ZE))
            !
            punto_s%XE = Min (this%nodsrc%NodalSource(i)%c2p(ii)%XE, Max(BoundingBox%XI, BoundingBox%XE))
            punto_s%YE = Min (this%nodsrc%NodalSource(i)%c2p(ii)%YE, Max(BoundingBox%YI, BoundingBox%YE))
            IF ((punto_s%zc /= 0).and.(this%nodsrc%NodalSource(i)%isElec))  then !only in case of Ez
               punto_s%ZE = Min (this%nodsrc%NodalSource(i)%c2p(ii)%ZE, Max(BoundingBox%ZI, BoundingBox%ZE-1))
            else
               punto_s%ZE = Min (this%nodsrc%NodalSource(i)%c2p(ii)%ZE, Max(BoundingBox%ZI, BoundingBox%ZE  ))
            endif
            !
            punto_s%or = this%nodsrc%NodalSource(i)%c2p(ii)%or
            punto_s%xc = this%nodsrc%NodalSource(i)%c2p(ii)%xc
            punto_s%yc = this%nodsrc%NodalSource(i)%c2p(ii)%yc
            punto_s%zc = this%nodsrc%NodalSource(i)%c2p(ii)%zc
            !
            !
            DO k1 = punto_s%ZI, punto_s%ZE
               DO j1 = punto_s%YI, punto_s%YE
                  DO i1 = punto_s%XI, punto_s%XE
                     IF (punto_s%xc /= 0) THEN
                        MEDIO = media%sggMiEx (i1, j1, k1)
                        valido=.true.
                        !IF ( .NOT. this%nodsrc%NodalSource(i)%isHard) THEN
                        !  VALIDO = (sgg%Med(MEDIO)%Is%Dielectric) .OR. (sgg%Med(MEDIO)%Is%EDispersive) .OR. &
                        ! & (sgg%Med(MEDIO)%Is%MDispersive)
                        !ELSE
                        !  VALIDO = .TRUE.
                        !END IF
                        IF (this%nodsrc%NodalSource(i)%isElec) THEN
                           VALIDO = VALIDO .AND. ( .NOT. sgg%Med(MEDIO)%Is%PEC)
                        END IF
                        write (buff,*) 'WARNING: Ex Nodal source on PEC media will be ignored (', i1, j1, k1,')'
                        IF ( .NOT. VALIDO) CALL  WarnErrReport (buff)
                        !
                        !
                        !
                        !!!MEDIO = sggmiHx (i1, j1, k1)
                        !!!valido=.true.
                        !!!!IF ( .NOT. this%nodsrc%NodalSource(i)%isHard) THEN
                        !!!!  VALIDO = (sgg%Med(MEDIO)%Is%Dielectric) .OR. (sgg%Med(MEDIO)%Is%EDispersive) .OR. &
                        !!!! & (sgg%Med(MEDIO)%Is%MDispersive)
                        !!!!ELSE
                        !!!!  VALIDO = .TRUE.
                        !!!!END IF
                        !!!IF ( .NOT. this%nodsrc%NodalSource(i)%isElec) THEN
                        !!!   VALIDO = VALIDO .AND. ( .NOT. sgg%Med(MEDIO)%Is%PMC)
                        !!!END IF
                        !!!write (buff,*) 'WARNING: Hx Nodal source on PMC media will be ignored (', i1, j1, k1,')'
                        !!!IF ( .NOT. VALIDO) CALL  WarnErrReport (buff)
                     END IF
                     !
                     !
                     IF (punto_s%yc /= 0) THEN
                        MEDIO = media%sggMiEy (i1, j1, k1)
                        valido=.true.
                        !IF ( .NOT. this%nodsrc%NodalSource(i)%isHard) THEN
                        !  VALIDO = (sgg%Med(MEDIO)%Is%Dielectric) .OR. (sgg%Med(MEDIO)%Is%EDispersive) .OR. &
                        ! & (sgg%Med(MEDIO)%Is%MDispersive)
                        !ELSE
                        !  VALIDO = .TRUE.
                        !END IF
                        IF (this%nodsrc%NodalSource(i)%isElec) THEN
                           VALIDO = VALIDO .AND. ( .NOT. sgg%Med(MEDIO)%Is%PEC)
                        END IF
                        write (buff,*) 'WARNING: Ey Nodal source on PMC media will be ignored (', i1, j1, k1,')'
                        IF ( .NOT. VALIDO) CALL  WarnErrReport (buff)
                        !
                        !
                        !!!MEDIO = sggmiHy (i1, j1, k1)
                        !!!valido=.true.
                        !!!!IF ( .NOT. this%nodsrc%NodalSource(i)%isHard) THEN
                        !!!!  VALIDO = (sgg%Med(MEDIO)%Is%Dielectric) .OR. (sgg%Med(MEDIO)%Is%EDispersive) .OR. &
                        !!!! & (sgg%Med(MEDIO)%Is%MDispersive)
                        !!!!ELSE
                        !!!!  VALIDO = .TRUE.
                        !!!!END IF
                        !!!IF ( .NOT. this%nodsrc%NodalSource(i)%isElec) THEN
                        !!!   VALIDO = VALIDO .AND. ( .NOT. sgg%Med(MEDIO)%Is%PMC)
                        !!!END IF
                        !!!write (buff,*) 'WARNING: Hy Nodal source on PMC media will be ignored (', i1, j1, k1,')'
                        !!!IF ( .NOT. VALIDO) CALL  WarnErrReport (buff)
                     END IF
                     !
                     !
                     IF (punto_s%zc /= 0) THEN
                        MEDIO = media%sggMiEz (i1, j1, k1)
                        valido=.true.
                        !IF ( .NOT. this%nodsrc%NodalSource(i)%isHard) THEN
                        !  VALIDO = (sgg%Med(MEDIO)%Is%Dielectric) .OR. (sgg%Med(MEDIO)%Is%EDispersive) .OR. &
                        ! & (sgg%Med(MEDIO)%Is%MDispersive)
                        !ELSE
                        !  VALIDO = .TRUE.
                        !END IF
                        IF (this%nodsrc%NodalSource(i)%isElec) THEN
                           VALIDO = VALIDO .AND. ( .NOT. sgg%Med(MEDIO)%Is%PEC)
                        END IF
                        write (buff,*) 'WARNING: Ez Nodal source on PMC media will be ignored (', i1, j1, k1,')'
                        IF ( .NOT. VALIDO) CALL  WarnErrReport (buff)
                        !
                        !
                        !!!MEDIO = sggmiHz (i1, j1, k1)
                        !!!valido=.true.
                        !!!!IF ( .NOT. this%nodsrc%NodalSource(i)%isHard) THEN
                        !!!!  VALIDO = (sgg%Med(MEDIO)%Is%Dielectric) .OR. (sgg%Med(MEDIO)%Is%EDispersive) .OR. &
                        !!!! & (sgg%Med(MEDIO)%Is%MDispersive)
                        !!!!ELSE
                        !!!!  VALIDO = .TRUE.
                        !!!!END IF
                        !!!IF ( .NOT. this%nodsrc%NodalSource(i)%isElec) THEN
                        !!!   VALIDO = VALIDO .AND. ( .NOT. sgg%Med(MEDIO)%Is%PMC)
                        !!!END IF
                        !!!write (buff,*) 'WARNING: Hz Nodal source on PMC media will be ignored (', i1, j1, k1,')'
                        !!!IF ( .NOT. VALIDO) CALL  WarnErrReport (buff)
                     END IF
                  END DO
               END DO
            END DO
            !
            !
            IF ((punto_s%XI <= punto_s%XE) .AND. (punto_s%YI <= punto_s%YE) .AND. (punto_s%ZI <= punto_s%ZE)) THEN
               conta2 = conta2 + 1
               sgg%NodalSource(conta1)%punto(conta2)%or = punto_s%or
               sgg%NodalSource(conta1)%punto(conta2)%xc = punto_s%xc
               sgg%NodalSource(conta1)%punto(conta2)%yc = punto_s%yc
               sgg%NodalSource(conta1)%punto(conta2)%zc = punto_s%zc
               sgg%NodalSource(conta1)%punto(conta2)%XI = punto_s%XI
               sgg%NodalSource(conta1)%punto(conta2)%XE = punto_s%XE
               sgg%NodalSource(conta1)%punto(conta2)%YI = punto_s%YI
               sgg%NodalSource(conta1)%punto(conta2)%YE = punto_s%YE
               sgg%NodalSource(conta1)%punto(conta2)%ZI = punto_s%ZI
               sgg%NodalSource(conta1)%punto(conta2)%ZE = punto_s%ZE
               !PARA ACOMODAR LAS NODAL SOURCE COMO MEDIOS LINE Y PODER VISUALIZAR SONDAS 010824
               sgg%Med(contamedia)%Is%Dielectric = .TRUE.
               sgg%Med(contamedia)%Is%LINE = .TRUE.
               sgg%Med(contamedia)%Priority = prior_IL
               sgg%Med(contamedia)%Epr =  1.0
               sgg%Med(contamedia)%Sigma = 0.
               sgg%Med(contamedia)%Mur =  1.0
               sgg%Med(contamedia)%SigmaM = 0.
               punto%XI = punto_s%XI
               punto%XE = punto_s%XE
               punto%YI = punto_s%YI
               punto%YE = punto_s%YE
               punto%ZI = punto_s%ZI
               punto%ZE = punto_s%ZE
               orientacion = punto_s%or
               isathinwire = .FALSE.
               numertag = 37
               CALL CreateLineMM (layoutnumber, media%sggMtag, tag_numbers, numertag, media%sggMiEx, media%sggMiEy, media%sggMiEz, &
                  media%sggMiHx, media%sggMiHy, media%sggMiHz, Alloc_iEx_XI, &
                  Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI, &
                  Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI, &
                  Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI, &
                  Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI, &
                  Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg%Med, sgg%NumMedia, sgg%EShared, BoundingBox, punto, orientacion, &
                  contamedia, isathinwire,verbose,numeroasignaciones)
            END IF
            sgg%NodalSource(conta1)%numpuntos = conta2 !update with the correct value
         END DO
      END DO
      !
      IF (allocated(contapuntos)) deallocate (contapuntos)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!PROBES
      !!!!!!!!!!!!!!!!!!!!!PROBES
      !!!!!!!!!!!!!!!!!!!!!PROBES
      !!!!!!!!!!!!!!!!!!!!!PROBES
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!PRIMERO LAS CUENTO
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!MasSondas     fields
      !
      tamaScrPrb = this%VolPrb%length
      if (createmapvtk) tamaScrPrb=tamaScrPrb+1
      tamaScrPrb = (tamaScrPrb)*3 !!!210618 allocateo el triple por si las volumicas son de tipo tifr mezcladas
      tamaSonda = this%Sonda%length
      tamaoldSONDA = this%oldSONDA%n_probes
      tamaBloquePrb = this%BloquePRB%N_BP
      !probes totales
      sgg%NumberRequest = (tamaSonda) + (tamaoldSONDA) + (tamaBloquePrb) + (tamaScrPrb)
      ALLOCATE (sgg%observation(1:sgg%NumberRequest))
      !inicializacion
      sgg%observation(1:sgg%NumberRequest)%nP            =-1
      sgg%observation(1:sgg%NumberRequest)%InitialTime   =-1
      sgg%observation(1:sgg%NumberRequest)%FinalTime     =-1
      sgg%observation(1:sgg%NumberRequest)%TimeStep      =-1
      sgg%observation(1:sgg%NumberRequest)%InitialFreq   =-1.0_RKIND
      sgg%observation(1:sgg%NumberRequest)%FinalFreq     =-1.0_RKIND
      sgg%observation(1:sgg%NumberRequest)%FreqStep      =-1.0_RKIND
      sgg%observation(1:sgg%NumberRequest)%outputrequest =' '
      sgg%observation(1:sgg%NumberRequest)%FileNormalize =' '
      sgg%observation(1:sgg%NumberRequest)%FreqDomain    =.false.
      sgg%observation(1:sgg%NumberRequest)%TimeDomain    =.false.
      sgg%observation(1:sgg%NumberRequest)%Saveall       =.false.
      sgg%observation(1:sgg%NumberRequest)%TRANSFER      =.false.
      sgg%observation(1:sgg%NumberRequest)%Volumic       =.false.

      !
      !ahora las cuento por bloques
      !
      DO i = 1, tamaSonda
         ii = i
         sgg%observation(ii)%nP = 0
         tama2 = (this%Sonda%collection(i)%len_cor)
         DO j = 1, tama2
            tipotemp = this%Sonda%collection(i)%cordinates(j)%or
            punto%XI = this%Sonda%collection(i)%cordinates(j)%XI
            punto%YI = this%Sonda%collection(i)%cordinates(j)%YI
            punto%ZI = this%Sonda%collection(i)%cordinates(j)%ZI
            IF ((punto%XI >= BoundingBox%XI) .AND. (punto%XI <= BoundingBox%XE) .AND. (punto%YI >= BoundingBox%YI) .AND.   &
            &    (punto%YI <= BoundingBox%YE) .AND. (punto%ZI >= BoundingBox%ZI) .AND. (punto%ZI <= BoundingBox%ZE) .AND.   &
            &   ((tipotemp == NP_COR_EX) .OR. (tipotemp == NP_COR_EY) .OR. (tipotemp == NP_COR_EZ) .OR.   &
            &     (tipotemp == NP_COR_HX) .OR. (tipotemp ==  NP_COR_HY) .OR. (tipotemp == NP_COR_HZ))) THEN
               sgg%observation(ii)%nP = sgg%observation(ii)%nP + 1
            ELSE IF (tipotemp == NP_COR_WIRECURRENT .or. tipotemp == NP_COR_CHARGE) THEN
               nodo_cazado=.false.
               loop_busqueda1: DO j1 = 1, this%twires%n_tw
                  DO i1 = 1, this%twires%TW(j1)%N_TWC
                     !nodo cazado
                     IF (this%twires%TW(j1)%TWC(i1)%nd == this%Sonda%collection(i)%cordinates(j)%XI) THEN
                        punto%XI = this%twires%TW(j1)%TWC(i1)%i
                        punto%YI = this%twires%TW(j1)%TWC(i1)%j
                        punto%ZI = this%twires%TW(j1)%TWC(i1)%k
                        nodo_cazado=.true.
                        !
                        IF ((punto%XI >= BoundingBox%XI) .AND. (punto%XI <= BoundingBox%XE).AND.(punto%YI >= BoundingBox%YI) .AND. &
                        & (punto%YI <= BoundingBox%YE) .AND. (punto%ZI >= BoundingBox%ZI) .AND. (punto%ZI <= BoundingBox%ZE)) THEN
                           !
                           SELECT CASE (this%twires%TW(j1)%TWC(i1)%D)
                            CASE (iEx, iEy, iEz)
                              sgg%observation(ii)%nP = sgg%observation(ii)%nP + 1
                           END SELECT
                           EXIT loop_busqueda1
                        END IF
                     END IF
                  END DO
               END DO loop_busqueda1
               !si no lo ha cazado... probamos SLANTED
               IF (.not.nodo_cazado) then
                  loop_busqueda2: DO j1 = 1, this%swires%n_sw
                     DO i1 = 1, this%swires%SW(j1)%N_SWC
                        !nodo cazado
                        IF (this%swires%SW(j1)%SWC(i1)%nd == this%Sonda%collection(i)%cordinates(j)%XI) THEN
                           punto%XI = floor(this%swires%SW(j1)%SWC(i1)%x)
                           punto%YI = floor(this%swires%SW(j1)%SWC(i1)%y)
                           punto%ZI = floor(this%swires%SW(j1)%SWC(i1)%z)
                           nodo_cazado=.true.
                           !
                           IF ((punto%XI >= BoundingBox%XI) .AND. (punto%XI <= BoundingBox%XE).AND.(punto%YI >= BoundingBox%YI) .AND. &
                           & (punto%YI <= BoundingBox%YE) .AND. (punto%ZI >= BoundingBox%ZI) .AND. (punto%ZI <= BoundingBox%ZE)) THEN
                              sgg%observation(ii)%nP = sgg%observation(ii)%nP + 1
                              !
                              EXIT loop_busqueda2
                           END IF
                        END IF
                     END DO
                  END DO loop_busqueda2
               END IF
               !si no lo ha cazado
               IF (.not.nodo_cazado) then
                  write(buff,'(a,i9)') 'Current probe not found in WIRE segment ',this%Sonda%collection(i)%cordinates(j)%XI
                  call StopOnError(layoutnumber,size,buff)
               endif
            ELSE IF (this%Sonda%collection(i)%cordinates(j)%or == NP_COR_DDP) THEN
               if (run_with_dmma) then
                  nodo_cazado=.false.
                  do_loop_busquedatg1: DO j1 = 1, this%tSlots%n_tg
                     DO i1 = 1, this%tSlots%Tg(j1)%N_tgc
                        !nodo cazado
                        IF (this%tSlots%Tg(j1)%TgC(i1)%node == this%Sonda%collection(i)%cordinates(j)%XI) THEN
                           punto%XI = this%tSlots%Tg(j1)%TgC(i1)%i
                           punto%YI = this%tSlots%Tg(j1)%TgC(i1)%j
                           punto%ZI = this%tSlots%Tg(j1)%TgC(i1)%k
                           nodo_cazado=.true.
                           IF ((punto%XI >= BoundingBox%XI) .AND. (punto%XI <= BoundingBox%XE) .AND. (punto%YI >= BoundingBox%YI) .AND. &
                           & (punto%YI <= BoundingBox%YE) .AND. (punto%ZI >= BoundingBox%ZI) .AND. (punto%ZI <= BoundingBox%ZE)) THEN
                              sgg%observation(i)%nP = sgg%observation(i)%nP + 1
                              EXIT do_loop_busquedatg1
                           END IF
                        END IF
                     END DO
                  END DO do_loop_busquedatg1
                  !si no lo ha cazado
                  IF (.not.nodo_cazado) then
                     write(buff,'(a,i9)') 'Voltage probe not found ',this%Sonda%collection(i)%cordinates(j)%XI
                     call StopOnError(layoutnumber,size,buff)
                  endif
               else
                  write(buff,'(a,i9)') 'ERROR: Voltage probe in gaps only available under -dmma flag '
                  call StopOnError(layoutnumber,size,buff)
               endif !del run_with_dmma
            else if (abs(tipotemp) == NP_COR_LINE) then
               sgg%observation(ii)%nP = sgg%observation(ii)%nP + 1
            END IF
         END DO
      END DO

      !Sondas
      !
      DO i = 1, tamaoldSONDA
         !acumulador
         ii = i + tamaSonda
         !far fields (no es time domain pero una forma especial de ellos)
         tama2 = (this%oldSONDA%probes(i)%n_FarField)
         if (tama2 > 1) then
            buff='Only one Far Field probe allowed'
            CALL STOPONERROR(layoutnumber,size,buff)
         endif
         if (tama2 > 0) sgg%observation(ii)%nP = 1 !un punto para todo el farfield (es simbolico)
         !electric FIELDS
         tama2 = (this%oldSONDA%probes(i)%n_Electric)
         DO j = 1, tama2
            tama3 = (this%oldSONDA%probes(i)%Electric(j)%probe%n_cord)
            tama4 = (this%oldSONDA%probes(i)%Electric(j)%probe%n_cord)
            tama5 = (this%oldSONDA%probes(i)%Electric(j)%probe%n_cord)
            tama6 = (this%oldSONDA%probes(i)%Electric(j)%probe%n_cord)
            buff='TAMANIOS DE PROBES EH RAROS'
            IF ((tama3 /= tama4) .OR. (tama3 /= tama5) .OR. (tama3 /= tama6)) &
            &                CALL STOPONERROR(layoutnumber,size,buff)
            DO k = 1, tama3
               punto%XI = this%oldSONDA%probes(i)%Electric(j)%probe%i(k)
               punto%YI = this%oldSONDA%probes(i)%Electric(j)%probe%j(k)
               punto%ZI = this%oldSONDA%probes(i)%Electric(j)%probe%k(k)
               IF ((punto%XI >= BoundingBox%XI) .AND. (punto%XI <= BoundingBox%XE) .AND. (punto%YI >= BoundingBox%YI) .AND. &
               & (punto%YI <= BoundingBox%YE) .AND. (punto%ZI >= BoundingBox%ZI) .AND. (punto%ZI <= BoundingBox%ZE)) THEN
                  sgg%observation(ii)%nP = sgg%observation(ii)%nP + 3
               END IF
            END DO
         END DO
         !MAGNETIC FIELDS
         tama2 = (this%oldSONDA%probes(i)%n_Magnetic)
         DO j = 1, tama2
            tama3 = (this%oldSONDA%probes(i)%Magnetic(j)%probe%n_cord)
            tama4 = (this%oldSONDA%probes(i)%Magnetic(j)%probe%n_cord)
            tama5 = (this%oldSONDA%probes(i)%Magnetic(j)%probe%n_cord)
            tama6 = (this%oldSONDA%probes(i)%Magnetic(j)%probe%n_cord)
            buff='pre1_ERROR: TAMANIOS DE PROBES EH RAROS'
            IF ((tama3 /= tama4) .OR. (tama3 /= tama5) .OR. (tama3 /= tama6)) &
            &           CALL STOPONERROR(layoutnumber,size,buff)
            DO k = 1, tama3
               punto%XI = this%oldSONDA%probes(i)%Magnetic(j)%probe%i(k)
               punto%YI = this%oldSONDA%probes(i)%Magnetic(j)%probe%j(k)
               punto%ZI = this%oldSONDA%probes(i)%Magnetic(j)%probe%k(k)
               IF ((punto%XI >= BoundingBox%XI) .AND. (punto%XI <= BoundingBox%XE) .AND. (punto%YI >= BoundingBox%YI) .AND. &
               & (punto%YI <= BoundingBox%YE) .AND. (punto%ZI >= BoundingBox%ZI) .AND. (punto%ZI <= BoundingBox%ZE)) THEN
                  sgg%observation(ii)%nP = sgg%observation(ii)%nP + 3
               END IF
            END DO
         END DO
      END DO
      !Bloque Current Probes
      !
      DO i = 1, tamaBloquePrb
         ii = i + tamaSonda + tamaoldSONDA
         sgg%observation(ii)%nP = 0
         punto%XI = this%BloquePRB%BP(i)%i1
         punto%YI = this%BloquePRB%BP(i)%j1
         punto%ZI = this%BloquePRB%BP(i)%k1
         punto%XE = this%BloquePRB%BP(i)%I2
         punto%YE = this%BloquePRB%BP(i)%J2
         punto%ZE = this%BloquePRB%BP(i)%K2
         !!!
         IF (((punto%XI >= BoundingBox%XI) .OR. (punto%XI <= BoundingBox%XE)) .AND. ((punto%YI >= BoundingBox%YI).OR. (punto%YI <= BoundingBox%YE)) .AND. &
            ((punto%ZI >= BoundingBox%ZI) .OR. (punto%ZI <= BoundingBox%ZE)) .AND. ((punto%XE >= BoundingBox%XI).OR. (punto%XE <= BoundingBox%XE)) .AND. &
            ((punto%YE >= BoundingBox%YI) .OR. (punto%YE <= BoundingBox%YE)) .AND. ((punto%ZE >= BoundingBox%ZI).OR. (punto%ZE <= BoundingBox%ZE))) THEN
            SELECT CASE (this%BloquePRB%BP(i)%NML)
             CASE (iEx)
               DO k = this%BloquePRB%BP(i)%i1, this%BloquePRB%BP(i)%I2, this%BloquePRB%BP(i)%skip
                  sgg%observation(ii)%nP = sgg%observation(ii)%nP + 1
               END DO
             CASE (iEy)
               DO k = this%BloquePRB%BP(i)%j1, this%BloquePRB%BP(i)%J2, this%BloquePRB%BP(i)%skip
                  sgg%observation(ii)%nP = sgg%observation(ii)%nP + 1
               END DO
             CASE (iEz)
               DO k = this%BloquePRB%BP(i)%k1, this%BloquePRB%BP(i)%K2, this%BloquePRB%BP(i)%skip
                  sgg%observation(ii)%nP = sgg%observation(ii)%nP + 1
               END DO
            END SELECT
         END IF
         !DEL TAMA DEL Bloque CURRENT PROBES
      END DO

      !Volumic probes (similar to MasSondas PERO CON PUNTOS FINALES COMO LAS Bloque PROBES)
      !ahora las cuento por bloques
      DO i = 1, tamaScrPrb/3  !!!210618 En realidad hay un tercio
         ii = i + tamaSonda + tamaoldSONDA + tamaBloquePrb
         sgg%observation(ii)%nP = 0
         ! crea una sonda vtk vacion con el instante incial 0 a efectos de mapa
         if (createmapvtk.and.(i==tamaScrPrb/3)) then !!!210618 En realidad hay un tercio
            tama2=1
            DO j = 1, tama2
               tipotemp = mapvtk
               punto%XI = SINPML_fullsize(iHx)%XI !!! +1   !ojo si se cambia aqui tambien mas abajo !le quito 1 para que con condiciones PEC no las pinte !habria que manejar el dibujo de las condiciones aparte
               punto%YI = SINPML_fullsize(iHy)%YI !!! +1
               punto%ZI = SINPML_fullsize(iHz)%ZI !!! +1
               punto%XE = SINPML_fullsize(iHx)%XE !!! -1
               punto%YE = SINPML_fullsize(iHy)%YE !!! -1
               punto%ZE = SINPML_fullsize(iHz)%ZE !!! -1
!!!               print *,layoutnumber,punto%XI,punto%YI,punto%ZI,punto%XE,punto%YE,punto%ZE
               IF (((punto%XI >= BoundingBox%XI) .OR. (punto%XI <= BoundingBox%XE)) .AND. &
               &   ((punto%YI >= BoundingBox%YI) .OR. (punto%YI <= BoundingBox%YE)) .AND. &
               &   ((punto%ZI >= BoundingBox%ZI) .OR. (punto%ZI <= BoundingBox%ZE)) .AND. &
               &   ((punto%XE >= BoundingBox%XI) .OR. (punto%XE <= BoundingBox%XE)) .AND. &
               &   ((punto%YE >= BoundingBox%YI) .OR. (punto%YE <= BoundingBox%YE)) .AND. &
               &   ((punto%ZE >= BoundingBox%ZI) .OR. (punto%ZE <= BoundingBox%ZE))) THEN
!!!                   print *,'----Dentro->',layoutnumber,punto%XI,punto%YI,punto%ZI,punto%XE,punto%YE,punto%ZE
                  sgg%observation(ii)%nP = sgg%observation(ii)%nP + 1
               ENDIF
            end do
         else
            tama2 = (this%VolPrb%collection(i)%len_cor)
            DO j = 1, tama2
               tipotemp = this%VolPrb%collection(i)%cordinates(j)%or
               punto%XI = this%VolPrb%collection(i)%cordinates(j)%XI
               punto%YI = this%VolPrb%collection(i)%cordinates(j)%YI
               punto%ZI = this%VolPrb%collection(i)%cordinates(j)%ZI
               punto%XE = this%VolPrb%collection(i)%cordinates(j)%XE
               punto%YE = this%VolPrb%collection(i)%cordinates(j)%YE
               punto%ZE = this%VolPrb%collection(i)%cordinates(j)%ZE

               IF (((punto%XI >= BoundingBox%XI) .OR. (punto%XI <= BoundingBox%XE)) .AND. ((punto%YI >= BoundingBox%YI) .OR.  &
               &     (punto%YI <= BoundingBox%YE)) .AND. ((punto%ZI >= BoundingBox%ZI) .OR. (punto%ZI <= BoundingBox%ZE)) .AND.  &
               &    ((punto%XE >= BoundingBox%XI) .OR. (punto%XE <= BoundingBox%XE)) .AND. ((punto%YE >= BoundingBox%YI) .OR.   &
               &     (punto%YE <= BoundingBox%YE)) .AND. ((punto%ZE >= BoundingBox%ZI) .OR. (punto%ZE <= BoundingBox%ZE))) THEN

                  sgg%observation(ii)%nP = sgg%observation(ii)%nP + 1
               ENDIF
            END DO
         endif
      END DO
      ! si se lanza con -mapvtk se crea una slice probe para ver la estructura

      !Ahora creo los puntos de observacion
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      sondas = 0
      DO ii = 1, sgg%NumberRequest
         sondas = sondas + sgg%observation(ii)%nP
      END DO
      !
      IF (sondas > Maxprobes) THEN
         WRITE (buff,*) 'Too many probes= ', sondas, 'Either  or reduce the number of probes below ', &
         & Maxprobes
         CALL STOPONERROR(layoutnumber,size,buff)
      END IF
      IF (sondas > 1024) THEN
         WRITE (buff,*) 'Number of probes ', &
         &   ' is over 1024. UNIX limit is typically 1024'// &
            ' Check! It will crash if ulimit -n is not set > ', sondas
         call WarnErrReport(buff)
      END IF

      !luego chequeo las memoria de las sondas si se van de memoria en observation.f90
      !      IF (sondas*BuffObse*4 > MaxMemoryProbes) THEN
      !        WRITE (buff,*) 'Too much memory for the probes= ', sondas * BuffObse * 4, 'Probes= ', sondas,   &
      !         &     'Either reduce the number o&
      !       &f probes or recompile decreasing BuffObse ', BuffObse, 'or increasing ', MaxMemoryProbes
      !        CALL STOPONERROR(layoutnumber,size,buff)
      !      END IF
      !
      IF (sgg%NumberRequest /= 0) THEN
         !alocateo
         DO ii = 1, sgg%NumberRequest
            ALLOCATE (sgg%observation(ii)%P(1:sgg%observation(ii)%nP))
            sgg%observation(ii)%P(1:sgg%observation(ii)%nP)%what=nothing !bug peligroso 2012
            sgg%observation(ii)%TimeDomain = .false.
            sgg%observation(ii)%FreqDomain = .FALSE.
            sgg%observation(ii)%TRANSFER = .FALSE.
            sgg%observation(ii)%Volumic = .FALSE.
            sgg%observation(ii)%FileNormalize=' '
            !trancos
            sgg%observation(ii)%P(1:sgg%observation(ii)%nP)%Xtrancos = 1 !default
            sgg%observation(ii)%P(1:sgg%observation(ii)%nP)%Ytrancos = 1 !default
            sgg%observation(ii)%P(1:sgg%observation(ii)%nP)%Ztrancos = 1 !default
            !fin trancos
            !al final debe quedar igual
            !lo reseteo porque lo reutilizo de contador
            sgg%observation(ii)%nP=0 !reset it
            !
         END DO
         DO i = 1, tamaSonda
            ii = i
            tama2 = (this%Sonda%collection(i)%len_cor)
            !
            SELECT CASE (this%Sonda%collection(i)%type2)
             CASE (NP_T2_time)
               sgg%observation(ii)%TimeDomain = .TRUE.
               sgg%observation(ii)%FreqDomain = .FALSE.
               sgg%observation(ii)%TRANSFER = .FALSE.
             CASE (NP_T2_FREQ)
               !I will output everything in time and transform it later
               sgg%observation(ii)%TimeDomain = .TRUE.
               sgg%observation(ii)%FreqDomain = .TRUE.
               sgg%observation(ii)%TRANSFER = .FALSE.
               !                call STOPONERROR(layoutnumber,size,'ONLY TIME DOMAIN DATA IN NEW PROBE')
             CASE (NP_T2_TRANSFER)
               sgg%observation(ii)%TimeDomain = .TRUE.
               sgg%observation(ii)%FreqDomain = .TRUE.
               sgg%observation(ii)%TRANSFER = .TRUE.
               buff='Transfer function only in Frequency Domain'
               !!           CALL STOPONERROR(layoutnumber,size,buff)
             CASE (NP_T2_TIMEFREQ )
               sgg%observation(ii)%TimeDomain = .TRUE.
               sgg%observation(ii)%FreqDomain = .TRUE.
               sgg%observation(ii)%TRANSFER = .FALSE.
               !                call STOPONERROR(layoutnumber,size,'ONLY TIME DOMAIN DATA IN NEW PROBE')
             CASE (NP_T2_TIMETRANSF)
               sgg%observation(ii)%TimeDomain = .TRUE.
               sgg%observation(ii)%FreqDomain = .TRUE.
               sgg%observation(ii)%TRANSFER = .TRUE.
               !                call STOPONERROR(layoutnumber,size,'ONLY TIME DOMAIN DATA IN NEW PROBE')
             CASE (NP_T2_FREQTRANSF)
               sgg%observation(ii)%TimeDomain = .TRUE.
               sgg%observation(ii)%FreqDomain = .TRUE.
               sgg%observation(ii)%TRANSFER = .TRUE.
               !                call STOPONERROR(layoutnumber,size,'ONLY TIME DOMAIN DATA IN NEW PROBE')
             CASE (NP_T2_TIMEFRECTRANSF)
               sgg%observation(ii)%TimeDomain = .TRUE.
               sgg%observation(ii)%FreqDomain = .TRUE.
               sgg%observation(ii)%TRANSFER = .TRUE.
            END SELECT
            !repair info
            !
            IF (sgg%observation(ii)%FreqDomain) THEN
               !save everything for later transform
               sgg%observation(ii)%Saveall = .TRUE.
            END IF
            IF ((this%Sonda%collection(i)%type1 == NP_T1_AMBOS)) THEN
               sgg%observation(ii)%Saveall = .TRUE.
               sgg%observation(ii)%TimeDomain = .TRUE.
            ELSE
               continue
            END IF
            !
            WRITE (probenumber, '(i7)') ii
            !        sgg%observation(ii)%outputrequest=trim(adjustl(probenumber))//'_'// &
            !                                          trim(adjustl(this%Sonda%collection(i)%outputrequest))
            sgg%observation(ii)%outputrequest = trim (adjustl(this%Sonda%collection(i)%outputrequest))
            sgg%observation(ii)%InitialTime = this%Sonda%collection(i)%tstart
            sgg%observation(ii)%FinalTime = this%Sonda%collection(i)%tstop
            sgg%observation(ii)%TimeStep = this%Sonda%collection(i)%tstep
            IF ((sgg%observation(ii)%FinalTime<TINY(1.0_RKIND)).OR.(sgg%observation(ii)%TimeStep<TINY(1.0_RKIND))) sgg%observation(ii)%Saveall = .TRUE.
            sgg%observation(ii)%InitialFreq = this%Sonda%collection(i)%fstart
            sgg%observation(ii)%FinalFreq = this%Sonda%collection(i)%fstop
            sgg%observation(ii)%FreqStep = this%Sonda%collection(i)%fstep
            sgg%observation(ii)%FileNormalize = trim (adjustl(this%Sonda%collection(i)%filename))
            !!!
            IF ((sgg%observation(ii)%InitialFreq < 0.).or. &
               (sgg%observation(ii)%FinalFreq <= 1e-9).or. &
               (sgg%observation(ii)%FreqStep <= 1e-9)) then
               WRITE (buff,*) 'ERROR: Some incorrect frequency domain parameters (initial,final,step) ',sgg%observation(ii)%InitialFreq,sgg%observation(ii)%FinalFreq,sgg%observation(ii)%FreqStep
               if (sgg%observation(ii)%FreqDomain) CALL STOPONERROR(layoutnumber,size,buff)
            ENDIF
            !!!
            DO j = 1, tama2
               punto%XI = this%Sonda%collection(i)%cordinates(j)%XI
               punto%YI = this%Sonda%collection(i)%cordinates(j)%YI
               punto%ZI = this%Sonda%collection(i)%cordinates(j)%ZI
               IF (this%Sonda%collection(i)%cordinates(j)%or == NP_COR_EX) THEN
                  IF ((punto%XI >= BoundingBox%XI) .AND. (punto%XI <= BoundingBox%XE) .AND. (punto%YI >= BoundingBox%YI) .AND. &
                  & (punto%YI <= BoundingBox%YE) .AND. (punto%ZI >= BoundingBox%ZI) .AND. (punto%ZI <= BoundingBox%ZE)) THEN
                     sgg%observation(ii)%nP = sgg%observation(ii)%nP + 1
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%XI = punto%XI
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%YI = punto%YI
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%ZI = punto%ZI
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%What = iEx
                  END IF
               ELSE IF (this%Sonda%collection(i)%cordinates(j)%or == NP_COR_EY) THEN
                  IF ((punto%XI >= BoundingBox%XI) .AND. (punto%XI <= BoundingBox%XE) .AND. (punto%YI >= BoundingBox%YI) .AND. &
                  & (punto%YI <= BoundingBox%YE) .AND. (punto%ZI >= BoundingBox%ZI) .AND. (punto%ZI <= BoundingBox%ZE)) THEN
                     sgg%observation(ii)%nP = sgg%observation(ii)%nP + 1
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%XI = punto%XI
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%YI = punto%YI
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%ZI = punto%ZI
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%What = iEy
                  END IF
               ELSE IF (this%Sonda%collection(i)%cordinates(j)%or == NP_COR_EZ) THEN
                  IF ((punto%XI >= BoundingBox%XI) .AND. (punto%XI <= BoundingBox%XE) .AND. (punto%YI >= BoundingBox%YI) .AND. &
                  & (punto%YI <= BoundingBox%YE) .AND. (punto%ZI >= BoundingBox%ZI) .AND. (punto%ZI <= BoundingBox%ZE)) THEN
                     sgg%observation(ii)%nP = sgg%observation(ii)%nP + 1
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%XI = punto%XI
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%YI = punto%YI
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%ZI = punto%ZI
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%What = iEz
                  END IF
               ELSE IF (this%Sonda%collection(i)%cordinates(j)%or == NP_COR_HX) THEN
                  IF ((punto%XI >= BoundingBox%XI) .AND. (punto%XI <= BoundingBox%XE) .AND. (punto%YI >= BoundingBox%YI) .AND. &
                  & (punto%YI <= BoundingBox%YE) .AND. (punto%ZI >= BoundingBox%ZI) .AND. (punto%ZI <= BoundingBox%ZE)) THEN
                     sgg%observation(ii)%nP = sgg%observation(ii)%nP + 1
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%XI = punto%XI
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%YI = punto%YI
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%ZI = punto%ZI
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%What = iHx
                  END IF
               ELSE IF (this%Sonda%collection(i)%cordinates(j)%or == NP_COR_HY) THEN
                  IF ((punto%XI >= BoundingBox%XI) .AND. (punto%XI <= BoundingBox%XE) .AND. (punto%YI >= BoundingBox%YI) .AND. &
                  & (punto%YI <= BoundingBox%YE) .AND. (punto%ZI >= BoundingBox%ZI) .AND. (punto%ZI <= BoundingBox%ZE)) THEN
                     sgg%observation(ii)%nP = sgg%observation(ii)%nP + 1
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%XI = punto%XI
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%YI = punto%YI
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%ZI = punto%ZI
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%What = iHy
                  END IF
               ELSE IF (this%Sonda%collection(i)%cordinates(j)%or == NP_COR_HZ) THEN
                  IF ((punto%XI >= BoundingBox%XI) .AND. (punto%XI <= BoundingBox%XE) .AND. (punto%YI >= BoundingBox%YI) .AND. &
                  & (punto%YI <= BoundingBox%YE) .AND. (punto%ZI >= BoundingBox%ZI) .AND. (punto%ZI <= BoundingBox%ZE)) THEN
                     sgg%observation(ii)%nP = sgg%observation(ii)%nP + 1
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%XI = punto%XI
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%YI = punto%YI
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%ZI = punto%ZI
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%What = iHz
                  END IF
               ELSE IF (this%Sonda%collection(i)%cordinates(j)%or == NP_COR_WIRECURRENT .or. &
                        this%Sonda%collection(i)%cordinates(j)%or == NP_COR_CHARGE) THEN
                  nodo_cazado=.false.
                  do_loop_busqueda: DO j1 = 1, this%twires%n_tw
                     DO i1 = 1, this%twires%TW(j1)%N_TWC
                        !nodo cazado
                        IF (this%twires%TW(j1)%TWC(i1)%nd == this%Sonda%collection(i)%cordinates(j)%XI) THEN
                           nodo_cazado=.true.
                           punto%XI = this%twires%TW(j1)%TWC(i1)%i
                           punto%YI = this%twires%TW(j1)%TWC(i1)%j
                           punto%ZI = this%twires%TW(j1)%TWC(i1)%k
                           IF ((punto%XI >= BoundingBox%XI) .AND. (punto%XI <= BoundingBox%XE) .AND.   &
                           &    (punto%YI >= BoundingBox%YI) .AND. &
                           &    (punto%YI <= BoundingBox%YE) .AND. (punto%ZI >= BoundingBox%ZI) .AND. (punto%ZI <= BoundingBox%ZE)) THEN
                              SELECT CASE (this%twires%TW(j1)%TWC(i1)%D)
                               CASE (iEx)
                                 sgg%observation(i)%nP = sgg%observation(i)%nP + 1
                                 sgg%observation(i)%P(sgg%observation(i)%nP)%node = this%twires%TW(j1)%TWC(i1)%nd
                                 sgg%observation(i)%P(sgg%observation(i)%nP)%XI = punto%XI
                                 sgg%observation(i)%P(sgg%observation(i)%nP)%YI = punto%YI
                                 sgg%observation(i)%P(sgg%observation(i)%nP)%ZI = punto%ZI
                                 if (this%Sonda%collection(i)%cordinates(j)%or == NP_COR_WIRECURRENT) then 
                                    sgg%observation(i)%P(sgg%observation(i)%nP)%What = iJx
                                 else if (this%Sonda%collection(i)%cordinates(j)%or == NP_COR_CHARGE) then 
                                    sgg%observation(i)%P(sgg%observation(i)%nP)%What = iQx
                                 end if

                                 !se nota con un indice distinto
                                 !
                               CASE (iEy)
                                 sgg%observation(i)%nP = sgg%observation(i)%nP + 1
                                 sgg%observation(i)%P(sgg%observation(i)%nP)%node = this%twires%TW(j1)%TWC(i1)%nd
                                 sgg%observation(i)%P(sgg%observation(i)%nP)%XI = punto%XI
                                 sgg%observation(i)%P(sgg%observation(i)%nP)%YI = punto%YI
                                 sgg%observation(i)%P(sgg%observation(i)%nP)%ZI = punto%ZI
                                 if (this%Sonda%collection(i)%cordinates(j)%or == NP_COR_WIRECURRENT) then 
                                    sgg%observation(i)%P(sgg%observation(i)%nP)%What = iJy
                                 else if (this%Sonda%collection(i)%cordinates(j)%or == NP_COR_CHARGE) then 
                                    sgg%observation(i)%P(sgg%observation(i)%nP)%What = iQy
                                 end if
                               CASE (iEz)
                                 sgg%observation(i)%nP = sgg%observation(i)%nP + 1
                                 sgg%observation(i)%P(sgg%observation(i)%nP)%node = this%twires%TW(j1)%TWC(i1)%nd
                                 sgg%observation(i)%P(sgg%observation(i)%nP)%XI = punto%XI
                                 sgg%observation(i)%P(sgg%observation(i)%nP)%YI = punto%YI
                                 sgg%observation(i)%P(sgg%observation(i)%nP)%ZI = punto%ZI
                                 if (this%Sonda%collection(i)%cordinates(j)%or == NP_COR_WIRECURRENT) then 
                                    sgg%observation(i)%P(sgg%observation(i)%nP)%What = iJz
                                 else if (this%Sonda%collection(i)%cordinates(j)%or == NP_COR_CHARGE) then 
                                    sgg%observation(i)%P(sgg%observation(i)%nP)%What = iQz
                                 end if
                              END SELECT
                              EXIT do_loop_busqueda
                           END IF
                        END IF
                     END DO
                  END DO do_loop_busqueda
                  IF (.not.nodo_cazado) THEN
                     do_loop_busqueda3: DO j1 = 1, this%swires%n_sw
                        DO i1 = 1, this%swires%SW(j1)%N_SWC
                           IF (this%swires%SW(j1)%SWC(i1)%nd == this%Sonda%collection(i)%cordinates(j)%XI) THEN
                              nodo_cazado=.true.
                              punto%XI = floor(this%swires%SW(j1)%SWC(i1)%x)
                              punto%YI = floor(this%swires%SW(j1)%SWC(i1)%y)
                              punto%ZI = floor(this%swires%SW(j1)%SWC(i1)%z)
                              IF ((punto%XI >= BoundingBox%XI) .AND. (punto%XI <= BoundingBox%XE) .AND.  &
                              &    (punto%YI >= BoundingBox%YI) .AND. &
                              &    (punto%YI <= BoundingBox%YE) .AND. (punto%ZI >= BoundingBox%ZI) .AND. (punto%ZI <= BoundingBox%ZE)) THEN
                                 sgg%observation(i)%nP = sgg%observation(i)%nP + 1
                                 sgg%observation(i)%P(sgg%observation(i)%nP)%node = this%swires%SW(j1)%SWC(i1)%nd
                                 sgg%observation(i)%P(sgg%observation(i)%nP)%XI = floor(this%swires%SW(j1)%SWC(i1)%x)
                                 sgg%observation(i)%P(sgg%observation(i)%nP)%YI = floor(this%swires%SW(j1)%SWC(i1)%y)
                                 sgg%observation(i)%P(sgg%observation(i)%nP)%ZI = floor(this%swires%SW(j1)%SWC(i1)%z)
                                 if (this%Sonda%collection(i)%cordinates(j)%or == NP_COR_WIRECURRENT) then 
                                    sgg%observation(i)%P(sgg%observation(i)%nP)%What = iJx
                                 else if (this%Sonda%collection(i)%cordinates(j)%or == NP_COR_CHARGE) then 
                                    sgg%observation(i)%P(sgg%observation(i)%nP)%What = iQx
                                 end if
                                 EXIT do_loop_busqueda3
                              END IF
                           END IF
                        END DO
                     END DO do_loop_busqueda3
                  END IF
               ELSE IF (this%Sonda%collection(i)%cordinates(j)%or == NP_COR_DDP) THEN
                  do_loop_busquedatg: DO j1 = 1, this%tSlots%n_tg
                     DO i1 = 1, this%tSlots%Tg(j1)%N_tgc
                        !nodo cazado
                        IF (this%tSlots%Tg(j1)%TgC(i1)%node == this%Sonda%collection(i)%cordinates(j)%XI) THEN
                           punto%XI = this%tSlots%Tg(j1)%TgC(i1)%i
                           punto%YI = this%tSlots%Tg(j1)%TgC(i1)%j
                           punto%ZI = this%tSlots%Tg(j1)%TgC(i1)%k
                           IF ((punto%XI >= BoundingBox%XI) .AND. (punto%XI <= BoundingBox%XE) .AND.   &
                           &     (punto%YI >= BoundingBox%YI) .AND. &
                           & (punto%YI <= BoundingBox%YE) .AND. (punto%ZI >= BoundingBox%ZI) .AND. (punto%ZI <= BoundingBox%ZE)) THEN
                              sgg%observation(i)%nP = sgg%observation(i)%nP + 1
                              sgg%observation(i)%P(sgg%observation(i)%nP)%node = this%tSlots%Tg(j1)%TgC(i1)%node
                              sgg%observation(i)%P(sgg%observation(i)%nP)%XI = punto%XI
                              sgg%observation(i)%P(sgg%observation(i)%nP)%YI = punto%YI
                              sgg%observation(i)%P(sgg%observation(i)%nP)%ZI = punto%ZI
                              direccion = this%tSlots%Tg(j1)%TgC(i1)%dir
                              SELECT CASE (this%tSlots%Tg(j1)%TgC(i1)%or)
                               CASE (iEx)
                                 SELECT CASE (direccion)
                                  CASE (iEz)
                                    sgg%observation(i)%P(sgg%observation(i)%nP)%What = iVy
                                  CASE (iEy)
                                    sgg%observation(i)%P(sgg%observation(i)%nP)%What = iVz
                                 END SELECT
                               CASE (iEy)
                                 SELECT CASE (direccion)
                                  CASE (iEx)
                                    sgg%observation(i)%P(sgg%observation(i)%nP)%What = iVz
                                  CASE (iEz)
                                    sgg%observation(i)%P(sgg%observation(i)%nP)%What = iVx
                                 END SELECT
                               CASE (iEz)
                                 SELECT CASE (direccion)
                                  CASE (iEy)
                                    sgg%observation(i)%P(sgg%observation(i)%nP)%What = iVx
                                  CASE (iEx)
                                    sgg%observation(i)%P(sgg%observation(i)%nP)%What = iVy
                                 END SELECT
                              END SELECT
                              EXIT do_loop_busquedatg
                           END IF
                        END IF
                     END DO
                  END DO do_loop_busquedatg
               else if (abs(this%Sonda%collection(i)%cordinates(j)%or) == NP_COR_LINE) then 
                  block
                     integer (kind=4) :: line_size, obs_size, idx
                     !intrinsic size function coopted by size global variable...
                     line_size = ubound(this%Sonda%collection(i)%cordinates,1)-lbound(this%Sonda%collection(i)%cordinates,1)+1
                     sgg%observation(i)%nP = sgg%observation(i)%nP + 1
                     obs_size =  sgg%observation(i)%nP
                     allocate(sgg%observation(i)%P(obs_size)%line(line_size))
                     sgg%observation(i)%P(obs_size)%What = lineIntegral
                     do idx = 1, line_size
                        sgg%observation(i)%P(obs_size)%line(idx)%x = this%Sonda%collection(i)%cordinates(idx)%Xi
                        sgg%observation(i)%P(obs_size)%line(idx)%y = this%Sonda%collection(i)%cordinates(idx)%Yi
                        sgg%observation(i)%P(obs_size)%line(idx)%z = this%Sonda%collection(i)%cordinates(idx)%Zi
                        if (this%Sonda%collection(i)%cordinates(idx)%Xe /= -1) then 
                           sgg%observation(i)%P(obs_size)%line(idx)%orientation = sign(1, this%Sonda%collection(i)%cordinates(idx)%or)
                        else if (this%Sonda%collection(i)%cordinates(idx)%Ye /= -1) then 
                           sgg%observation(i)%P(obs_size)%line(idx)%orientation = sign(2, this%Sonda%collection(i)%cordinates(idx)%or)
                        else if (this%Sonda%collection(i)%cordinates(idx)%Ze /= -1) then 
                           sgg%observation(i)%P(obs_size)%line(idx)%orientation = sign(3, this%Sonda%collection(i)%cordinates(idx)%or)
                        end if
                     end do
                     sgg%observation(i)%P(obs_size)%XI = sgg%observation(i)%P(obs_size)%line(1)%x
                     sgg%observation(i)%P(obs_size)%YI = sgg%observation(i)%P(obs_size)%line(1)%y
                     sgg%observation(i)%P(obs_size)%ZI = sgg%observation(i)%P(obs_size)%line(1)%z
                  end block
               END IF
            END DO
         END DO
         !
         !Sondas propiamente dichas
         !
         DO i = 1, tamaoldSONDA
            ii = i + tamaSonda
            !only the MasSondas accept the freqdomain
            sgg%observation(ii)%TimeDomain = .TRUE. !NO  CONSIDERO EL FARFIELD FREQDOMAIN PQ LA TRATO BIEN COMO TIMEDOMAIN Y NO QUIERO JODERLA !26/02/14
            sgg%observation(ii)%FreqDomain = .FALSE.
            sgg%observation(ii)%TRANSFER = .FALSE.
            !farfields (no es time domain pero una forma especial de ellos)
            tama2 = (this%oldSONDA%probes(i)%n_FarField)
            WRITE (buff,*) 'More than 1 Far Field box unsupported'
            IF (tama2 > 1) CALL STOPONERROR(layoutnumber,size,buff)
            !
            DO j = 1, tama2
               tama3 = (this%oldSONDA%probes(i)%FarField(j)%probe%n_cord)
               buff='FAR FIELD PROBE REQUIRES TWO COORDINATES FOR THE BOX'
               IF (tama3 /= 2)  CALL STOPONERROR(layoutnumber,size,buff)
               !
               sgg%observation(ii)%nP = sgg%observation(ii)%nP + 1
               punto%XI = this%oldSONDA%probes(i)%FarField(j)%probe%i(1)
               punto%YI = this%oldSONDA%probes(i)%FarField(j)%probe%j(1)
               punto%ZI = this%oldSONDA%probes(i)%FarField(j)%probe%k(1)
               punto%XE = this%oldSONDA%probes(i)%FarField(j)%probe%i(2)
               punto%YE = this%oldSONDA%probes(i)%FarField(j)%probe%j(2)
               punto%ZE = this%oldSONDA%probes(i)%FarField(j)%probe%k(2)
               !
               sgg%observation(ii)%P(1)%XI = punto%XI
               sgg%observation(ii)%P(1)%YI = punto%YI
               sgg%observation(ii)%P(1)%ZI = punto%ZI
               sgg%observation(ii)%P(1)%XE = punto%XE
               sgg%observation(ii)%P(1)%YE = punto%YE
               sgg%observation(ii)%P(1)%ZE = punto%ZE
               sgg%observation(ii)%P(1)%what = farfield
               !
               !no se clipea porque se manejan como ondas planas
               !
               sgg%observation(ii)%InitialFreq   =this%oldSONDA%probes(i)%FarField(j)%probe%fstart
               sgg%observation(ii)%FinalFreq     =this%oldSONDA%probes(i)%FarField(j)%probe%fstop
               sgg%observation(ii)%FreqStep      =this%oldSONDA%probes(i)%FarField(j)%probe%fstep
               !
               sgg%observation(ii)%thetaStart    =this%oldSONDA%probes(i)%FarField(j)%probe%thetaStart
               sgg%observation(ii)%thetaStop     =this%oldSONDA%probes(i)%FarField(j)%probe%thetaStop
               sgg%observation(ii)%thetaStep     =this%oldSONDA%probes(i)%FarField(j)%probe%thetaStep
               !
               sgg%observation(ii)%phiStart   =this%oldSONDA%probes(i)%FarField(j)%probe%phiStart
               sgg%observation(ii)%phiStop    =this%oldSONDA%probes(i)%FarField(j)%probe%phiStop
               sgg%observation(ii)%phiStep    =this%oldSONDA%probes(i)%FarField(j)%probe%phiStep
               sgg%observation(ii)%FileNormalize    =this%oldSONDA%probes(i)%FarField(j)%probe%FileNormalize

               sgg%observation(ii)%outputrequest=trim(adjustl(this%oldSONDA%probes(i)%FarField(j)%probe%outputrequest(1:3)))

               IF ((sgg%observation(ii)%InitialFreq < 0.).or. &
                  (sgg%observation(ii)%FinalFreq <= 1e-9).or. &
                  (sgg%observation(ii)%FreqStep <= 1e-9)) then
                  WRITE (buff,*) 'ERROR: Some incorrect frequency domain parameters (initial,final,step) ',sgg%observation(ii)%InitialFreq,sgg%observation(ii)%FinalFreq,sgg%observation(ii)%FreqStep
                  if (sgg%observation(ii)%FreqDomain) CALL STOPONERROR(layoutnumber,size,buff)
               ENDIF
               !

            END DO
            !electric FIELDS
            tama2 = (this%oldSONDA%probes(i)%n_Electric)
            DO j = 1, tama2
               tama3 = (this%oldSONDA%probes(i)%Electric(j)%probe%n_cord)
               tama4 = (this%oldSONDA%probes(i)%Electric(j)%probe%n_cord)
               tama5 = (this%oldSONDA%probes(i)%Electric(j)%probe%n_cord)
               tama6 = (this%oldSONDA%probes(i)%Electric(j)%probe%n_cord)
               buff='TAMANIOS DE PROBES EH RAROS'
               IF ((tama3 /= tama4) .OR. (tama3 /= tama5) .OR. (tama3 /= tama6))  &
               &                CALL STOPONERROR(layoutnumber,size,buff)
               WRITE (probenumber, '(i7)') ii
               DO k = 1, tama3
                  punto%XI = this%oldSONDA%probes(i)%Electric(j)%probe%i(k)
                  punto%YI = this%oldSONDA%probes(i)%Electric(j)%probe%j(k)
                  punto%ZI = this%oldSONDA%probes(i)%Electric(j)%probe%k(k)
                  IF ((punto%XI >= BoundingBox%XI) .AND. (punto%XI <= BoundingBox%XE) .AND. (punto%YI >= BoundingBox%YI) .AND. &
                  & (punto%YI <= BoundingBox%YE) .AND. (punto%ZI >= BoundingBox%ZI) .AND. (punto%ZI <= BoundingBox%ZE)) THEN
                     sgg%observation(ii)%nP = sgg%observation(ii)%nP + 1
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%XI = punto%XI
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%YI = punto%YI
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%ZI = punto%ZI
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%What = iEx
                     !
                     sgg%observation(ii)%nP = sgg%observation(ii)%nP + 1
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%XI = punto%XI
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%YI = punto%YI
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%ZI = punto%ZI
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%What = iEy
                     !
                     sgg%observation(ii)%nP = sgg%observation(ii)%nP + 1
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%XI = punto%XI
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%YI = punto%YI
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%ZI = punto%ZI
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%What = iEz
                     !
                     !
                     sgg%observation(ii)%InitialTime = this%oldSONDA%probes(i)%Electric(j)%probe%tstart
                     sgg%observation(ii)%FinalTime = this%oldSONDA%probes(i)%Electric(j)%probe%tstop
                     sgg%observation(ii)%TimeStep = this%oldSONDA%probes(i)%Electric(j)%probe%tstep
                     IF ((sgg%observation(ii)%FinalTime<TINY(1.0_RKIND)).OR.(sgg%observation(ii)%TimeStep<TINY(1.0_RKIND))) sgg%observation(ii)%Saveall = .TRUE.
                     sgg%observation(ii)%outputrequest = trim (adjustl(this%oldSONDA%probes(i)%Electric(j)%probe%outputrequest))
                     !
                  END IF
               END DO
            END DO
            !MAGNETIC FIELDS
            tama2 = (this%oldSONDA%probes(i)%n_Magnetic)
            DO j = 1, tama2
               tama3 = (this%oldSONDA%probes(i)%Magnetic(j)%probe%n_cord)
               tama4 = (this%oldSONDA%probes(i)%Magnetic(j)%probe%n_cord)
               tama5 = (this%oldSONDA%probes(i)%Magnetic(j)%probe%n_cord)
               tama6 = (this%oldSONDA%probes(i)%Magnetic(j)%probe%n_cord)
               buff='TAMANIOS DE PROBES EH RAROS'
               IF ((tama3 /= tama4) .OR. (tama3 /= tama5) .OR. (tama3 /= tama6))  &
               &                CALL STOPONERROR(layoutnumber,size,buff)
               DO k = 1, tama3
                  punto%XI = this%oldSONDA%probes(i)%Magnetic(j)%probe%i(k)
                  punto%YI = this%oldSONDA%probes(i)%Magnetic(j)%probe%j(k)
                  punto%ZI = this%oldSONDA%probes(i)%Magnetic(j)%probe%k(k)
                  IF ((punto%XI >= BoundingBox%XI) .AND. (punto%XI <= BoundingBox%XE) .AND. (punto%YI >= BoundingBox%YI) .AND. &
                  & (punto%YI <= BoundingBox%YE) .AND. (punto%ZI >= BoundingBox%ZI) .AND. (punto%ZI <= BoundingBox%ZE)) THEN
                     sgg%observation(ii)%nP = sgg%observation(ii)%nP + 1
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%XI = punto%XI
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%YI = punto%YI
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%ZI = punto%ZI
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%What = iHx
                     !
                     sgg%observation(ii)%nP = sgg%observation(ii)%nP + 1
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%XI = punto%XI
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%YI = punto%YI
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%ZI = punto%ZI
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%What = iHy
                     !
                     sgg%observation(ii)%nP = sgg%observation(ii)%nP + 1
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%XI = punto%XI
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%YI = punto%YI
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%ZI = punto%ZI
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%What = iHz
                     !
                     !
                     sgg%observation(ii)%InitialTime = this%oldSONDA%probes(i)%Magnetic(j)%probe%tstart
                     sgg%observation(ii)%FinalTime = this%oldSONDA%probes(i)%Magnetic(j)%probe%tstop
                     sgg%observation(ii)%TimeStep = this%oldSONDA%probes(i)%Magnetic(j)%probe%tstep

                     IF ((sgg%observation(ii)%FinalTime<TINY(1.0_RKIND)).OR.(sgg%observation(ii)%TimeStep<TINY(1.0_RKIND))) sgg%observation(ii)%Saveall = .TRUE.
                     sgg%observation(ii)%outputrequest = trim (adjustl(this%oldSONDA%probes(i)%Magnetic(j)%probe%outputrequest))
                  END IF
               END DO
            END DO
            !ojo faltan por implementar
            !traditional probes
         END DO
         !Bloque current probes
         !
         DO i = 1, tamaBloquePrb
            ii = i + tamaSonda + tamaoldSONDA
            punto%XI = this%BloquePRB%BP(i)%i1
            punto%YI = this%BloquePRB%BP(i)%j1
            punto%ZI = this%BloquePRB%BP(i)%k1
            punto%XE = this%BloquePRB%BP(i)%I2
            punto%YE = this%BloquePRB%BP(i)%J2
            punto%ZE = this%BloquePRB%BP(i)%K2
            IF (((punto%XI >= BoundingBox%XI) .OR. (punto%XI <= BoundingBox%XE)) .AND. ((punto%YI >= BoundingBox%YI) .OR.   &
            &     (punto%YI <= BoundingBox%YE)) .AND. ((punto%ZI >= BoundingBox%ZI) .OR. (punto%ZI <= BoundingBox%ZE)) .AND.   &
            &     ((punto%XE >= BoundingBox%XI) .OR. (punto%XE <= BoundingBox%XE)) .AND. ((punto%YE >= BoundingBox%YI) .OR.   &
            &     (punto%YE <= BoundingBox%YE)) .AND. ((punto%ZE >= BoundingBox%ZI) .OR. (punto%ZE <= BoundingBox%ZE))) THEN
               !
               !!!ANIADIDO 15/07/15 PARA COMPATIBILIDD NEW PROBE EN FRECUENCIA
               SELECT CASE (this%BloquePRB%BP(i)%type2)
                CASE (NP_T2_time)
                  sgg%observation(ii)%TimeDomain = .TRUE.
                  sgg%observation(ii)%FreqDomain = .FALSE.
                  sgg%observation(ii)%TRANSFER = .FALSE.
                CASE (NP_T2_FREQ)
                  !I will output everything in time and transform it later
                  sgg%observation(ii)%TimeDomain = .TRUE.
                  sgg%observation(ii)%FreqDomain = .TRUE.
                  sgg%observation(ii)%TRANSFER = .FALSE.
                  !                call STOPONERROR(layoutnumber,size,'ONLY TIME DOMAIN DATA IN NEW PROBE')
                CASE (NP_T2_TRANSFER)
                  sgg%observation(ii)%TimeDomain = .TRUE.
                  sgg%observation(ii)%FreqDomain = .TRUE.
                  sgg%observation(ii)%TRANSFER = .TRUE.
                  buff='Transfer function only in Frequency Domain'
                  !!           CALL STOPONERROR(layoutnumber,size,buff)
                CASE (NP_T2_TIMEFREQ )
                  sgg%observation(ii)%TimeDomain = .TRUE.
                  sgg%observation(ii)%FreqDomain = .TRUE.
                  sgg%observation(ii)%TRANSFER = .FALSE.
                  !                call STOPONERROR(layoutnumber,size,'ONLY TIME DOMAIN DATA IN NEW PROBE')
                CASE (NP_T2_TIMETRANSF)
                  sgg%observation(ii)%TimeDomain = .TRUE.
                  sgg%observation(ii)%FreqDomain = .TRUE.
                  sgg%observation(ii)%TRANSFER = .TRUE.
                  !                call STOPONERROR(layoutnumber,size,'ONLY TIME DOMAIN DATA IN NEW PROBE')
                CASE (NP_T2_FREQTRANSF)
                  sgg%observation(ii)%TimeDomain = .TRUE.
                  sgg%observation(ii)%FreqDomain = .TRUE.
                  sgg%observation(ii)%TRANSFER = .TRUE.
                  !                call STOPONERROR(layoutnumber,size,'ONLY TIME DOMAIN DATA IN NEW PROBE')
                CASE (NP_T2_TIMEFRECTRANSF)
                  sgg%observation(ii)%TimeDomain = .TRUE.
                  sgg%observation(ii)%FreqDomain = .TRUE.
                  sgg%observation(ii)%TRANSFER = .TRUE.
               END SELECT
               !repair info
               !
               IF (sgg%observation(ii)%FreqDomain) THEN
                  !save everything for later transform
                  sgg%observation(ii)%Saveall = .TRUE.
               END IF
               !!!          END IF
               !
               WRITE (probenumber, '(i7)') ii
               !        sgg%observation(ii)%outputrequest=trim(adjustl(probenumber))//'_'// &
               !                                          trim(adjustl(this%BloquePRB%BP(i)%outputrequest))
               sgg%observation(ii)%outputrequest = trim (adjustl(this%BloquePRB%BP(i)%outputrequest))
               sgg%observation(ii)%InitialTime = this%BloquePRB%BP(i)%tstart
               sgg%observation(ii)%FinalTime = this%BloquePRB%BP(i)%tstop
               sgg%observation(ii)%TimeStep = this%BloquePRB%BP(i)%tstep
               IF ((sgg%observation(ii)%FinalTime<TINY(1.0_RKIND)).OR.(sgg%observation(ii)%TimeStep<TINY(1.0_RKIND))) sgg%observation(ii)%Saveall = .TRUE.
               sgg%observation(ii)%InitialFreq = this%BloquePRB%BP(i)%fstart
               sgg%observation(ii)%FinalFreq = this%BloquePRB%BP(i)%fstop
               sgg%observation(ii)%FreqStep = this%BloquePRB%BP(i)%fstep
               sgg%observation(ii)%FileNormalize = trim (adjustl(this%BloquePRB%BP(i)%FileNormalize))

               IF ((sgg%observation(ii)%InitialFreq < 0.).or. &
                  (sgg%observation(ii)%FinalFreq <= 1e-9).or. &
                  (sgg%observation(ii)%FreqStep <= 1e-9) ) then
                  WRITE (buff,*) 'ERROR: Some incorrect frequency domain parameters (initial,final,step) ',sgg%observation(ii)%InitialFreq,sgg%observation(ii)%FinalFreq,sgg%observation(ii)%FreqStep
                  if (sgg%observation(ii)%FreqDomain) CALL STOPONERROR(layoutnumber,size,buff)
               ENDIF
               !FIN COMPATIBILIDAD 15/07/15
               SELECT CASE (this%BloquePRB%BP(i)%NML)
                CASE (iEx)
                  DO k = this%BloquePRB%BP(i)%i1, this%BloquePRB%BP(i)%I2, this%BloquePRB%BP(i)%skip
                     sgg%observation(ii)%nP = sgg%observation(ii)%nP + 1
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%XI = k
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%YI = punto%YI
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%ZI = punto%ZI
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%XE = k
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%YE = punto%YE
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%ZE = punto%ZE
                     IF (this%BloquePRB%BP(i)%t) THEN
                        !electric TYPE
                        sgg%observation(ii)%P(sgg%observation(ii)%nP)%What = iBloqueJx
                     ELSE
                        !MAGNETIC TYPE
                        sgg%observation(ii)%P(sgg%observation(ii)%nP)%What = iBloqueMx
                     END IF
                  END DO
                CASE (iEy)
                  DO k = this%BloquePRB%BP(i)%j1, this%BloquePRB%BP(i)%J2, this%BloquePRB%BP(i)%skip
                     sgg%observation(ii)%nP = sgg%observation(ii)%nP + 1
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%XI = punto%XI
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%YI = k
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%ZI = punto%ZI
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%XE = punto%XE
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%YE = k
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%ZE = punto%ZE
                     IF (this%BloquePRB%BP(i)%t) THEN
                        !electric TYPE
                        sgg%observation(ii)%P(sgg%observation(ii)%nP)%What = iBloqueJy
                     ELSE
                        !MAGNETIC TYPE
                        sgg%observation(ii)%P(sgg%observation(ii)%nP)%What = iBloqueMy
                     END IF
                  END DO
                CASE (iEz)
                  DO k = this%BloquePRB%BP(i)%k1, this%BloquePRB%BP(i)%K2, this%BloquePRB%BP(i)%skip
                     sgg%observation(ii)%nP = sgg%observation(ii)%nP + 1
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%XI = punto%XI
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%YI = punto%YI
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%ZI = k
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%XE = punto%XE
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%YE = punto%YE
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%ZE = k
                     IF (this%BloquePRB%BP(i)%t) THEN
                        !electric TYPE
                        sgg%observation(ii)%P(sgg%observation(ii)%nP)%What = iBloqueJz
                     ELSE
                        !MAGNETIC TYPE
                        sgg%observation(ii)%P(sgg%observation(ii)%nP)%What = iBloqueMz
                     END IF
                  END DO
               END SELECT
            END IF
            !DEL TAMA DEL Bloque CURRENT PROBES
         END DO

         !Volumic probes (similar to MasSondas PERO CON PUNTOS FINALES COMO LAS Bloque PROBES)
         !ahora las cuento por bloques
         !
         memo=0
!!!210618
         DO i = 1, tamaScrPrb/3 !!!210618 En realidad hay un tercio
            ii = i + tamaSonda + tamaoldSONDA + tamaBloquePrb
            !crea sonda vtk al final de mapeo
            if (createmapvtk.and.(i==tamaScrPrb/3)) then !!!210618 En realidad hay un tercio
               sgg%observation(ii)%TimeDomain = .TRUE.
               sgg%observation(ii)%FreqDomain = .FALSE.
               sgg%observation(ii)%TRANSFER = .FALSE.
               sgg%observation(ii)%saveall = .FALSE.
               sgg%observation(ii)%nP = 0
               sgg%observation(ii)%Volumic = .true.
               sgg%observation(ii)%InitialTime =     sgg%dt
               sgg%observation(ii)%FinalTime =        sgg%dt+sgg%dt/300.0_RKIND
               sgg%observation(ii)%TimeStep =         sgg%dt !SACA SOLO UNO
               sgg%observation(ii)%outputrequest = ' '
               sgg%observation(ii)%InitialFreq = 0.0_RKIND
               sgg%observation(ii)%FinalFreq =  0.0_RKIND
               sgg%observation(ii)%FreqStep =  0.0_RKIND
               sgg%observation(ii)%FileNormalize = ''
               tama2 = 1
               IF (tama2 >1 ) THEN
                  WRITE (buff,*) 'Only 1 Volumic probe allown per section'
                  CALL STOPONERROR(layoutnumber,size,buff)
               END IF
               DO j = 1, tama2
                  !I clip these probes to allow out-of-the box snapshot probes !ojo si se cambia aqui tambien mas arriba
                  tipotemp = mapvtk
                  punto%XI = SINPML_fullsize(iHx)%XI   !!! +1
                  punto%YI = SINPML_fullsize(iHy)%YI   !!! +1
                  punto%ZI = SINPML_fullsize(iHz)%ZI   !!! +1
                  punto%XE = SINPML_fullsize(iHx)%XE   !!! -1
                  punto%YE = SINPML_fullsize(iHy)%YE   !!! -1
                  punto%ZE = SINPML_fullsize(iHz)%ZE   !!! -1

                  memo=memo+(punto%XE-punto%XI+1)*(punto%YE-punto%YI+1)*(punto%ZE-punto%ZI+1)

                  IF (( (punto%XI >= BoundingBox%XI) .OR. (punto%XI <= BoundingBox%XE)) .AND. &
                  &    ((punto%YI >= BoundingBox%YI) .OR. (punto%YI <= BoundingBox%YE)) .AND. &
                  &    ((punto%ZI >= BoundingBox%ZI) .OR. (punto%ZI <= BoundingBox%ZE)) .AND. &
                  &    ((punto%XE >= BoundingBox%XI) .OR. (punto%XE <= BoundingBox%XE)) .AND. &
                  &    ((punto%YE >= BoundingBox%YI) .OR. (punto%YE <= BoundingBox%YE)) .AND. &
                  &    ((punto%ZE >= BoundingBox%ZI) .OR. (punto%ZE <= BoundingBox%ZE))) THEN

                     !
                     WRITE (probenumber, '(i7)') ii
                     !            sgg%observation(ii)%outputrequest=trim(adjustl(probenumber))//'_'// &
                     !                                              trim(adjustl(this%BloquePrb%BP(i)%outputrequest))
                     !
                     sgg%observation(ii)%nP = sgg%observation(ii)%nP + 1
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%XI = punto%XI
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%YI = punto%YI
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%ZI = punto%ZI
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%XE = punto%XE
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%YE = punto%YE
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%ZE = punto%ZE
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%what = tipotemp
                  ENDIF
               END DO
!!!!210618 tambien se crrean extras dummy para los vtk
               sgg%observation(  tamaScrPrb/3+ii)%TimeDomain = .false.
               sgg%observation(  tamaScrPrb/3+ii)%FreqDomain = .false.
               sgg%observation(  tamaScrPrb/3+ii)%TRANSFER =   .false.
               if (associated(this%VolPrb%collection).and.(tamaScrPrb/=0).and.(i<=this%VolPrb%length)) then !280618 & 220319
                  sgg%observation(  tamaScrPrb/3+ii)%outputrequest=trim (adjustl( this%VolPrb%collection(i)%outputrequest))//'_df_'
               else
                  sgg%observation(  tamaScrPrb/3+ii)%outputrequest=' ' !es un dummy mapvtk sin nombre 280618
               endif
               sgg%observation(  tamaScrPrb/3+ii)%nP= sgg%observation(               ii)%np
               ALLOCATE (sgg%observation(  tamaScrPrb/3+ii)%P(1:sgg%observation(  tamaScrPrb/3+ii)%nP))
               sgg%observation(  tamaScrPrb/3+ii)%P(1:sgg%observation(  tamaScrPrb/3+ii)%nP)%what = nothing
!
               sgg%observation(2*tamaScrPrb/3+ii)%TimeDomain = .false.
               sgg%observation(2*tamaScrPrb/3+ii)%FreqDomain = .false.
               sgg%observation(2*tamaScrPrb/3+ii)%TRANSFER =   .false.
               if (associated(this%VolPrb%collection).and.(tamaScrPrb/=0).and.(i<=this%VolPrb%length)) then !280618 & 220319
                  sgg%observation(2*tamaScrPrb/3+ii)%outputrequest=trim (adjustl( this%VolPrb%collection(i)%outputrequest))//'_tr_'
               else
                  sgg%observation(2*tamaScrPrb/3+ii)%outputrequest=' ' !es un dummy mapvtk sin nombre 280618
               endif
               sgg%observation(2*tamaScrPrb/3+ii)%nP= sgg%observation(               ii)%np
               ALLOCATE (sgg%observation(2*tamaScrPrb/3+ii)%P(1:sgg%observation(2*tamaScrPrb/3+ii)%nP))
               sgg%observation(2*tamaScrPrb/3+ii)%P(1:sgg%observation(2*tamaScrPrb/3+ii)%nP)%what = nothing

!!!210618 triplica info sondas de frequencia
               sgg%observation(tamaScrPrb/3+ii)%Volumic                                            =  sgg%observation(ii)%Volumic
               sgg%observation(tamaScrPrb/3+ii)%InitialTime                                        =  sgg%observation(ii)%InitialTime
               sgg%observation(tamaScrPrb/3+ii)%FinalTime                                          =  sgg%observation(ii)%FinalTime
               sgg%observation(tamaScrPrb/3+ii)%TimeStep                                           =  sgg%observation(ii)%TimeStep
               sgg%observation(tamaScrPrb/3+ii)%InitialFreq                                        =  sgg%observation(ii)%InitialFreq
               sgg%observation(tamaScrPrb/3+ii)%FinalFreq                                          =  sgg%observation(ii)%FinalFreq
               sgg%observation(tamaScrPrb/3+ii)%FreqStep                                           =  sgg%observation(ii)%FreqStep
               sgg%observation(tamaScrPrb/3+ii)%FileNormalize                                      =  sgg%observation(ii)%FileNormalize
               sgg%observation(tamaScrPrb/3+ii)%P(1:sgg%observation(    tamaScrPrb/3+ii)%nP)%XI    =  sgg%observation(ii)%P(1:sgg%observation(ii)%nP)%XI
               sgg%observation(tamaScrPrb/3+ii)%P(1:sgg%observation(    tamaScrPrb/3+ii)%nP)%YI    =  sgg%observation(ii)%P(1:sgg%observation(ii)%nP)%YI
               sgg%observation(tamaScrPrb/3+ii)%P(1:sgg%observation(    tamaScrPrb/3+ii)%nP)%ZI    =  sgg%observation(ii)%P(1:sgg%observation(ii)%nP)%ZI
               sgg%observation(tamaScrPrb/3+ii)%P(1:sgg%observation(    tamaScrPrb/3+ii)%nP)%XE    =  sgg%observation(ii)%P(1:sgg%observation(ii)%nP)%XE
               sgg%observation(tamaScrPrb/3+ii)%P(1:sgg%observation(    tamaScrPrb/3+ii)%nP)%YE    =  sgg%observation(ii)%P(1:sgg%observation(ii)%nP)%YE
               sgg%observation(tamaScrPrb/3+ii)%P(1:sgg%observation(    tamaScrPrb/3+ii)%nP)%ZE    =  sgg%observation(ii)%P(1:sgg%observation(ii)%nP)%ZE
               !
               !trancos
               sgg%observation(tamaScrPrb/3+ii)%P(1:sgg%observation(    tamaScrPrb/3+ii)%nP)%Xtrancos = 1 !default
               sgg%observation(tamaScrPrb/3+ii)%P(1:sgg%observation(    tamaScrPrb/3+ii)%nP)%Ytrancos = 1 !default
               sgg%observation(tamaScrPrb/3+ii)%P(1:sgg%observation(    tamaScrPrb/3+ii)%nP)%Ztrancos = 1 !default
               !fin trancos

               !!!!!!!!!!!!!!!!!!
               sgg%observation(2*tamaScrPrb/3+ii)%Volumic                                          =  sgg%observation(ii)%Volumic
               sgg%observation(2*tamaScrPrb/3+ii)%InitialTime                                      =  sgg%observation(ii)%InitialTime
               sgg%observation(2*tamaScrPrb/3+ii)%FinalTime                                        =  sgg%observation(ii)%FinalTime
               sgg%observation(2*tamaScrPrb/3+ii)%TimeStep                                         =  sgg%observation(ii)%TimeStep
               sgg%observation(2*tamaScrPrb/3+ii)%InitialFreq                                      =  sgg%observation(ii)%InitialFreq
               sgg%observation(2*tamaScrPrb/3+ii)%FinalFreq                                        =  sgg%observation(ii)%FinalFreq
               sgg%observation(2*tamaScrPrb/3+ii)%FreqStep                                         =  sgg%observation(ii)%FreqStep
               sgg%observation(2*tamaScrPrb/3+ii)%FileNormalize                                    =  sgg%observation(ii)%FileNormalize
               sgg%observation(2*tamaScrPrb/3+ii)%P(1:sgg%observation(2*tamaScrPrb/3+ii)%nP)%XI    =  sgg%observation(ii)%P(1:sgg%observation(ii)%nP)%XI
               sgg%observation(2*tamaScrPrb/3+ii)%P(1:sgg%observation(2*tamaScrPrb/3+ii)%nP)%YI    =  sgg%observation(ii)%P(1:sgg%observation(ii)%nP)%YI
               sgg%observation(2*tamaScrPrb/3+ii)%P(1:sgg%observation(2*tamaScrPrb/3+ii)%nP)%ZI    =  sgg%observation(ii)%P(1:sgg%observation(ii)%nP)%ZI
               sgg%observation(2*tamaScrPrb/3+ii)%P(1:sgg%observation(2*tamaScrPrb/3+ii)%nP)%XE    =  sgg%observation(ii)%P(1:sgg%observation(ii)%nP)%XE
               sgg%observation(2*tamaScrPrb/3+ii)%P(1:sgg%observation(2*tamaScrPrb/3+ii)%nP)%YE    =  sgg%observation(ii)%P(1:sgg%observation(ii)%nP)%YE
               sgg%observation(2*tamaScrPrb/3+ii)%P(1:sgg%observation(2*tamaScrPrb/3+ii)%nP)%ZE    =  sgg%observation(ii)%P(1:sgg%observation(ii)%nP)%ZE
               !
               !trancos
               sgg%observation(2*tamaScrPrb/3+ii)%P(1:sgg%observation(2*tamaScrPrb/3+ii)%nP)%Xtrancos = 1 !default
               sgg%observation(2*tamaScrPrb/3+ii)%P(1:sgg%observation(2*tamaScrPrb/3+ii)%nP)%Ytrancos = 1 !default
               sgg%observation(2*tamaScrPrb/3+ii)%P(1:sgg%observation(2*tamaScrPrb/3+ii)%nP)%Ztrancos = 1 !default
               !fin trancos

            else !del mapvtk
!!!210618
               sgg%observation(ii)%Volumic = .true.
               sgg%observation(ii)%InitialTime =      this%VolPrb%collection(i)%tstart
               sgg%observation(ii)%FinalTime =        this%VolPrb%collection(i)%tstop
               sgg%observation(ii)%TimeStep =         this%VolPrb%collection(i)%tstep
               sgg%observation(ii)%outputrequest = trim (adjustl( this%VolPrb%collection(i)%outputrequest))
               sgg%observation(ii)%InitialFreq =      this%VolPrb%collection(i)%fstart
               sgg%observation(ii)%FinalFreq =        this%VolPrb%collection(i)%fstop
               sgg%observation(ii)%FreqStep =         this%VolPrb%collection(i)%fstep
               sgg%observation(ii)%FileNormalize = trim (adjustl(this%VolPrb%collection(i)%filename))
               sgg%observation(ii)%nP = 0
               tama2 = (this%VolPrb%collection(i)%len_cor)
               IF (tama2 >1 ) THEN
                  WRITE (buff,*) 'Only 1 Volumic probe allown per section'
                  CALL STOPONERROR(layoutnumber,size,buff)
               END IF
               DO j = 1, tama2
                  !I clip these probes to allow out-of-the box snapshot probes
                  tipotemp = this%VolPrb%collection(i)%cordinates(j)%or
                  punto%XI = max(this%VolPrb%collection(i)%cordinates(j)%XI,SINPML_fullsize(iEx)%XI)
                  punto%YI = max(this%VolPrb%collection(i)%cordinates(j)%YI,SINPML_fullsize(iEy)%YI)
                  punto%ZI = max(this%VolPrb%collection(i)%cordinates(j)%ZI,SINPML_fullsize(iEz)%ZI)
                  punto%XE = min(this%VolPrb%collection(i)%cordinates(j)%XE,SINPML_fullsize(iEx)%XE)
                  punto%YE = min(this%VolPrb%collection(i)%cordinates(j)%YE,SINPML_fullsize(iEy)%YE)
                  punto%ZE = min(this%VolPrb%collection(i)%cordinates(j)%ZE,SINPML_fullsize(iEz)%ZE)
                  memo=memo+(punto%XE-punto%XI+1)*(punto%YE-punto%YI+1)*(punto%ZE-punto%ZI+1)

                  IF (((punto%XI >= BoundingBox%XI) .OR. (punto%XI <= BoundingBox%XE)) .AND. &
                  &    ((punto%YI >= BoundingBox%YI) .OR. (punto%YI <= BoundingBox%YE)) .AND. &
                  &    ((punto%ZI >= BoundingBox%ZI) .OR. (punto%ZI <= BoundingBox%ZE)) .AND. &
                  &    ((punto%XE >= BoundingBox%XI) .OR. (punto%XE <= BoundingBox%XE)) .AND. &
                  &    ((punto%YE >= BoundingBox%YI) .OR. (punto%YE <= BoundingBox%YE)) .AND. &
                  &    ((punto%ZE >= BoundingBox%ZI) .OR. (punto%ZE <= BoundingBox%ZE))) THEN

                     !
                     WRITE (probenumber, '(i7)') ii
                     !            sgg%observation(ii)%outputrequest=trim(adjustl(probenumber))//'_'// &
                     !                                              trim(adjustl(this%BloquePrb%BP(i)%outputrequest))
                     !
                     sgg%observation(ii)%nP = sgg%observation(ii)%nP + 1
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%XI = punto%XI
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%YI = punto%YI
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%ZI = punto%ZI
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%XE = punto%XE
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%YE = punto%YE
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%ZE = punto%ZE
                     !trancos
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%Xtrancos = this%VolPrb%collection(i)%cordinates(sgg%observation(ii)%nP)%Xtrancos
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%Ytrancos = this%VolPrb%collection(i)%cordinates(sgg%observation(ii)%nP)%Ytrancos
                     sgg%observation(ii)%P(sgg%observation(ii)%nP)%Ztrancos = this%VolPrb%collection(i)%cordinates(sgg%observation(ii)%nP)%Ztrancos
                     !fin trancos

                  ENDIF
               END DO
!
               SELECT CASE (this%VolPrb%collection(i)%type2)
                CASE (NP_T2_time)

                  sgg%observation(ii)%TimeDomain = .TRUE.
                  sgg%observation(ii)%FreqDomain = .FALSE.
                  sgg%observation(ii)%TRANSFER = .FALSE.
                  sgg%observation(ii)%P(1:sgg%observation(ii)%nP)%what=tipotemp
!
                  sgg%observation(  tamaScrPrb/3+ii)%TimeDomain = .false.
                  sgg%observation(  tamaScrPrb/3+ii)%FreqDomain = .false.
                  sgg%observation(  tamaScrPrb/3+ii)%TRANSFER =   .false.
                  sgg%observation(  tamaScrPrb/3+ii)%outputrequest=trim (adjustl( this%VolPrb%collection(i)%outputrequest))//'_df_'
                  sgg%observation(  tamaScrPrb/3+ii)%nP= sgg%observation(               ii)%np
                  ALLOCATE (sgg%observation(  tamaScrPrb/3+ii)%P(1:sgg%observation(  tamaScrPrb/3+ii)%nP))
                  sgg%observation(  tamaScrPrb/3+ii)%P(1:sgg%observation(  tamaScrPrb/3+ii)%nP)%what = nothing
!
                  sgg%observation(2*tamaScrPrb/3+ii)%TimeDomain = .false.
                  sgg%observation(2*tamaScrPrb/3+ii)%FreqDomain = .false.
                  sgg%observation(2*tamaScrPrb/3+ii)%TRANSFER =   .false.
                  sgg%observation(2*tamaScrPrb/3+ii)%outputrequest=trim (adjustl( this%VolPrb%collection(i)%outputrequest))//'_tr_'
                  sgg%observation(2*tamaScrPrb/3+ii)%nP= sgg%observation(               ii)%np
                  ALLOCATE (sgg%observation(2*tamaScrPrb/3+ii)%P(1:sgg%observation(2*tamaScrPrb/3+ii)%nP))
                  sgg%observation(2*tamaScrPrb/3+ii)%P(1:sgg%observation(2*tamaScrPrb/3+ii)%nP)%what = nothing
!
                CASE (NP_T2_FREQ)
                  !I will TRANSFORM ON THE FLY
                  sgg%observation(ii)%TimeDomain = .false.
                  sgg%observation(ii)%FreqDomain = .false.
                  sgg%observation(ii)%TRANSFER = .false.
                  sgg%observation(ii)%P(1:sgg%observation(ii)%nP)%what=nothing  !el nothing debe predominar
!
                  sgg%observation(  tamaScrPrb/3+ii)%TimeDomain = .false.
                  sgg%observation(  tamaScrPrb/3+ii)%FreqDomain = .TRUE.
                  sgg%observation(  tamaScrPrb/3+ii)%TRANSFER =   .false.
                  sgg%observation(  tamaScrPrb/3+ii)%outputrequest=trim (adjustl( this%VolPrb%collection(i)%outputrequest))//'_df_'
                  sgg%observation(  tamaScrPrb/3+ii)%nP= sgg%observation(               ii)%np
                  ALLOCATE (sgg%observation(  tamaScrPrb/3+ii)%P(1:sgg%observation(  tamaScrPrb/3+ii)%nP))
                  sgg%observation(  tamaScrPrb/3+ii)%P(1:sgg%observation(  tamaScrPrb/3+ii)%nP)%what = tipotemp
!
                  sgg%observation(2*tamaScrPrb/3+ii)%TimeDomain = .false.
                  sgg%observation(2*tamaScrPrb/3+ii)%FreqDomain = .false.
                  sgg%observation(2*tamaScrPrb/3+ii)%TRANSFER =   .false.
                  sgg%observation(2*tamaScrPrb/3+ii)%outputrequest=trim (adjustl( this%VolPrb%collection(i)%outputrequest))//'_tr_'
                  sgg%observation(2*tamaScrPrb/3+ii)%nP= sgg%observation(               ii)%np
                  ALLOCATE (sgg%observation(2*tamaScrPrb/3+ii)%P(1:sgg%observation(2*tamaScrPrb/3+ii)%nP))
                  sgg%observation(2*tamaScrPrb/3+ii)%P(1:sgg%observation(2*tamaScrPrb/3+ii)%nP)%what = nothing
!
                CASE (NP_T2_TRANSFER)
                  !I will TRANSFORM ON THE FLY
                  sgg%observation(ii)%TimeDomain = .false.
                  sgg%observation(ii)%FreqDomain = .false.
                  sgg%observation(ii)%TRANSFER = .false.
                  sgg%observation(ii)%P(1:sgg%observation(ii)%nP)%what=nothing   !el nothing predomina sobre los true anteriores
!
                  sgg%observation(  tamaScrPrb/3+ii)%TimeDomain = .false.
                  sgg%observation(  tamaScrPrb/3+ii)%FreqDomain = .false.
                  sgg%observation(  tamaScrPrb/3+ii)%TRANSFER =   .false.
                  sgg%observation(  tamaScrPrb/3+ii)%outputrequest=trim (adjustl( this%VolPrb%collection(i)%outputrequest))//'_df_'
                  sgg%observation(  tamaScrPrb/3+ii)%nP= sgg%observation(               ii)%np
                  ALLOCATE (sgg%observation(  tamaScrPrb/3+ii)%P(1:sgg%observation(  tamaScrPrb/3+ii)%nP))
                  sgg%observation(  tamaScrPrb/3+ii)%P(1:sgg%observation(  tamaScrPrb/3+ii)%nP)%what = nothing
!
                  sgg%observation(2*tamaScrPrb/3+ii)%TimeDomain = .false.
                  sgg%observation(2*tamaScrPrb/3+ii)%FreqDomain = .true.
                  sgg%observation(2*tamaScrPrb/3+ii)%TRANSFER =   .true.
                  sgg%observation(2*tamaScrPrb/3+ii)%outputrequest=trim (adjustl( this%VolPrb%collection(i)%outputrequest))//'_tr_'
                  sgg%observation(2*tamaScrPrb/3+ii)%nP= sgg%observation(               ii)%np
                  ALLOCATE (sgg%observation(2*tamaScrPrb/3+ii)%P(1:sgg%observation(2*tamaScrPrb/3+ii)%nP))
                  sgg%observation(2*tamaScrPrb/3+ii)%P(1:sgg%observation(2*tamaScrPrb/3+ii)%nP)%what = tipotemp
                CASE (NP_T2_TIMEFREQ )
                  sgg%observation(ii)%TimeDomain = .true.
                  sgg%observation(ii)%FreqDomain = .false.
                  sgg%observation(ii)%TRANSFER = .false.
                  sgg%observation(ii)%P(1:sgg%observation(ii)%nP)%what=tipotemp   !el nothing predomina sobre los true anteriores
!
                  sgg%observation(  tamaScrPrb/3+ii)%TimeDomain = .false.
                  sgg%observation(  tamaScrPrb/3+ii)%FreqDomain = .TRUE.
                  sgg%observation(  tamaScrPrb/3+ii)%TRANSFER =   .false.
                  sgg%observation(  tamaScrPrb/3+ii)%outputrequest=trim (adjustl( this%VolPrb%collection(i)%outputrequest))//'_df_'
                  sgg%observation(  tamaScrPrb/3+ii)%nP= sgg%observation(               ii)%np
                  ALLOCATE (sgg%observation(  tamaScrPrb/3+ii)%P(1:sgg%observation(  tamaScrPrb/3+ii)%nP))
                  sgg%observation(  tamaScrPrb/3+ii)%P(1:sgg%observation(  tamaScrPrb/3+ii)%nP)%what = tipotemp
!
                  sgg%observation(2*tamaScrPrb/3+ii)%TimeDomain = .false.
                  sgg%observation(2*tamaScrPrb/3+ii)%FreqDomain = .false.
                  sgg%observation(2*tamaScrPrb/3+ii)%TRANSFER =   .false.
                  sgg%observation(2*tamaScrPrb/3+ii)%outputrequest=trim (adjustl( this%VolPrb%collection(i)%outputrequest))//'_tr_'
                  sgg%observation(2*tamaScrPrb/3+ii)%nP= sgg%observation(               ii)%np
                  ALLOCATE (sgg%observation(2*tamaScrPrb/3+ii)%P(1:sgg%observation(2*tamaScrPrb/3+ii)%nP))
                  sgg%observation(2*tamaScrPrb/3+ii)%P(1:sgg%observation(2*tamaScrPrb/3+ii)%nP)%what = nothing
                CASE (NP_T2_TIMETRANSF)
                  sgg%observation(ii)%TimeDomain = .true.
                  sgg%observation(ii)%FreqDomain = .false.
                  sgg%observation(ii)%TRANSFER = .false.
                  sgg%observation(ii)%P(1:sgg%observation(ii)%nP)%what=tipotemp   !el nothing predomina sobre los true anteriores
!
                  sgg%observation(  tamaScrPrb/3+ii)%TimeDomain = .false.
                  sgg%observation(  tamaScrPrb/3+ii)%FreqDomain = .false.
                  sgg%observation(  tamaScrPrb/3+ii)%TRANSFER =   .false.
                  sgg%observation(  tamaScrPrb/3+ii)%outputrequest=trim (adjustl( this%VolPrb%collection(i)%outputrequest))//'_df_'
                  sgg%observation(  tamaScrPrb/3+ii)%nP= sgg%observation(               ii)%np
                  ALLOCATE (sgg%observation(  tamaScrPrb/3+ii)%P(1:sgg%observation(  tamaScrPrb/3+ii)%nP))
                  sgg%observation(  tamaScrPrb/3+ii)%P(1:sgg%observation(  tamaScrPrb/3+ii)%nP)%what = nothing
!
                  sgg%observation(2*tamaScrPrb/3+ii)%TimeDomain = .false.
                  sgg%observation(2*tamaScrPrb/3+ii)%FreqDomain = .true.
                  sgg%observation(2*tamaScrPrb/3+ii)%TRANSFER =   .true.
                  sgg%observation(2*tamaScrPrb/3+ii)%outputrequest=trim (adjustl( this%VolPrb%collection(i)%outputrequest))//'_tr_'
                  sgg%observation(2*tamaScrPrb/3+ii)%nP= sgg%observation(               ii)%np
                  ALLOCATE (sgg%observation(2*tamaScrPrb/3+ii)%P(1:sgg%observation(2*tamaScrPrb/3+ii)%nP))
                  sgg%observation(2*tamaScrPrb/3+ii)%P(1:sgg%observation(2*tamaScrPrb/3+ii)%nP)%what = tipotemp
                CASE (NP_T2_FREQTRANSF)
                  sgg%observation(ii)%TimeDomain = .false.
                  sgg%observation(ii)%FreqDomain = .false.
                  sgg%observation(ii)%TRANSFER = .false.
                  sgg%observation(ii)%P(1:sgg%observation(ii)%nP)%what=nothing
!
                  sgg%observation(  tamaScrPrb/3+ii)%TimeDomain = .false.
                  sgg%observation(  tamaScrPrb/3+ii)%FreqDomain = .TRUE.
                  sgg%observation(  tamaScrPrb/3+ii)%TRANSFER =   .false.
                  sgg%observation(  tamaScrPrb/3+ii)%outputrequest=trim (adjustl( this%VolPrb%collection(i)%outputrequest))//'_df_'
                  sgg%observation(  tamaScrPrb/3+ii)%nP= sgg%observation(               ii)%np
                  ALLOCATE (sgg%observation(  tamaScrPrb/3+ii)%P(1:sgg%observation(  tamaScrPrb/3+ii)%nP))
                  sgg%observation(  tamaScrPrb/3+ii)%P(1:sgg%observation(  tamaScrPrb/3+ii)%nP)%what = tipotemp
!
                  sgg%observation(2*tamaScrPrb/3+ii)%TimeDomain = .false.
                  sgg%observation(2*tamaScrPrb/3+ii)%FreqDomain = .true.
                  sgg%observation(2*tamaScrPrb/3+ii)%TRANSFER =   .true.
                  sgg%observation(2*tamaScrPrb/3+ii)%outputrequest=trim (adjustl( this%VolPrb%collection(i)%outputrequest))//'_tr_'
                  sgg%observation(2*tamaScrPrb/3+ii)%nP= sgg%observation(               ii)%np
                  ALLOCATE (sgg%observation(2*tamaScrPrb/3+ii)%P(1:sgg%observation(2*tamaScrPrb/3+ii)%nP))
                  sgg%observation(2*tamaScrPrb/3+ii)%P(1:sgg%observation(2*tamaScrPrb/3+ii)%nP)%what = tipotemp
!
                CASE (NP_T2_TIMEFRECTRANSF)
                  sgg%observation(ii)%TimeDomain = .true.
                  sgg%observation(ii)%FreqDomain = .false.
                  sgg%observation(ii)%TRANSFER = .false.
                  sgg%observation(ii)%P(1:sgg%observation(ii)%nP)%what=tipotemp
!
                  sgg%observation(  tamaScrPrb/3+ii)%TimeDomain = .false.
                  sgg%observation(  tamaScrPrb/3+ii)%FreqDomain = .TRUE.
                  sgg%observation(  tamaScrPrb/3+ii)%TRANSFER =   .false.
                  sgg%observation(  tamaScrPrb/3+ii)%outputrequest=trim (adjustl( this%VolPrb%collection(i)%outputrequest))//'_df_'
                  sgg%observation(  tamaScrPrb/3+ii)%nP= sgg%observation(               ii)%np
                  ALLOCATE (sgg%observation(  tamaScrPrb/3+ii)%P(1:sgg%observation(  tamaScrPrb/3+ii)%nP))
                  sgg%observation(  tamaScrPrb/3+ii)%P(1:sgg%observation(  tamaScrPrb/3+ii)%nP)%what = tipotemp
!
                  sgg%observation(2*tamaScrPrb/3+ii)%TimeDomain = .false.
                  sgg%observation(2*tamaScrPrb/3+ii)%FreqDomain = .true.
                  sgg%observation(2*tamaScrPrb/3+ii)%TRANSFER =   .true.
                  sgg%observation(2*tamaScrPrb/3+ii)%outputrequest=trim (adjustl( this%VolPrb%collection(i)%outputrequest))//'_tr_'
                  sgg%observation(2*tamaScrPrb/3+ii)%nP= sgg%observation(               ii)%np
                  ALLOCATE (sgg%observation(2*tamaScrPrb/3+ii)%P(1:sgg%observation(2*tamaScrPrb/3+ii)%nP))
                  sgg%observation(2*tamaScrPrb/3+ii)%P(1:sgg%observation(2*tamaScrPrb/3+ii)%nP)%what = tipotemp
!
               END SELECT
!!!210618 triplica info sondas de frequencia
               sgg%observation(tamaScrPrb/3+ii)%Volumic                                            =  sgg%observation(ii)%Volumic
               sgg%observation(tamaScrPrb/3+ii)%InitialTime                                        =  sgg%observation(ii)%InitialTime
               sgg%observation(tamaScrPrb/3+ii)%FinalTime                                          =  sgg%observation(ii)%FinalTime
               sgg%observation(tamaScrPrb/3+ii)%TimeStep                                           =  sgg%observation(ii)%TimeStep
               sgg%observation(tamaScrPrb/3+ii)%InitialFreq                                        =  sgg%observation(ii)%InitialFreq
               sgg%observation(tamaScrPrb/3+ii)%FinalFreq                                          =  sgg%observation(ii)%FinalFreq
               sgg%observation(tamaScrPrb/3+ii)%FreqStep                                           =  sgg%observation(ii)%FreqStep
               sgg%observation(tamaScrPrb/3+ii)%FileNormalize                                      =  sgg%observation(ii)%FileNormalize
               sgg%observation(tamaScrPrb/3+ii)%P(1:sgg%observation(    tamaScrPrb/3+ii)%nP)%XI    =  sgg%observation(ii)%P(1:sgg%observation(ii)%nP)%XI
               sgg%observation(tamaScrPrb/3+ii)%P(1:sgg%observation(    tamaScrPrb/3+ii)%nP)%YI    =  sgg%observation(ii)%P(1:sgg%observation(ii)%nP)%YI
               sgg%observation(tamaScrPrb/3+ii)%P(1:sgg%observation(    tamaScrPrb/3+ii)%nP)%ZI    =  sgg%observation(ii)%P(1:sgg%observation(ii)%nP)%ZI
               sgg%observation(tamaScrPrb/3+ii)%P(1:sgg%observation(    tamaScrPrb/3+ii)%nP)%XE    =  sgg%observation(ii)%P(1:sgg%observation(ii)%nP)%XE
               sgg%observation(tamaScrPrb/3+ii)%P(1:sgg%observation(    tamaScrPrb/3+ii)%nP)%YE    =  sgg%observation(ii)%P(1:sgg%observation(ii)%nP)%YE
               sgg%observation(tamaScrPrb/3+ii)%P(1:sgg%observation(    tamaScrPrb/3+ii)%nP)%ZE    =  sgg%observation(ii)%P(1:sgg%observation(ii)%nP)%ZE
               !trancos
               sgg%observation(tamaScrPrb/3+ii)%P(1:sgg%observation(    tamaScrPrb/3+ii)%nP)%Xtrancos =sgg%observation(ii)%P(1:sgg%observation(ii)%nP)%Xtrancos
               sgg%observation(tamaScrPrb/3+ii)%P(1:sgg%observation(    tamaScrPrb/3+ii)%nP)%Ytrancos =sgg%observation(ii)%P(1:sgg%observation(ii)%nP)%Ytrancos
               sgg%observation(tamaScrPrb/3+ii)%P(1:sgg%observation(    tamaScrPrb/3+ii)%nP)%Ztrancos =sgg%observation(ii)%P(1:sgg%observation(ii)%nP)%Ztrancos
               !fin trancos
               !!!!!!!!!!!!!!!!!!
               sgg%observation(2*tamaScrPrb/3+ii)%Volumic                                          =  sgg%observation(ii)%Volumic
               sgg%observation(2*tamaScrPrb/3+ii)%InitialTime                                      =  sgg%observation(ii)%InitialTime
               sgg%observation(2*tamaScrPrb/3+ii)%FinalTime                                        =  sgg%observation(ii)%FinalTime
               sgg%observation(2*tamaScrPrb/3+ii)%TimeStep                                         =  sgg%observation(ii)%TimeStep
               sgg%observation(2*tamaScrPrb/3+ii)%InitialFreq                                      =  sgg%observation(ii)%InitialFreq
               sgg%observation(2*tamaScrPrb/3+ii)%FinalFreq                                        =  sgg%observation(ii)%FinalFreq
               sgg%observation(2*tamaScrPrb/3+ii)%FreqStep                                         =  sgg%observation(ii)%FreqStep
               sgg%observation(2*tamaScrPrb/3+ii)%FileNormalize                                    =  sgg%observation(ii)%FileNormalize
               sgg%observation(2*tamaScrPrb/3+ii)%P(1:sgg%observation(2*tamaScrPrb/3+ii)%nP)%XI    =  sgg%observation(ii)%P(1:sgg%observation(ii)%nP)%XI
               sgg%observation(2*tamaScrPrb/3+ii)%P(1:sgg%observation(2*tamaScrPrb/3+ii)%nP)%YI    =  sgg%observation(ii)%P(1:sgg%observation(ii)%nP)%YI
               sgg%observation(2*tamaScrPrb/3+ii)%P(1:sgg%observation(2*tamaScrPrb/3+ii)%nP)%ZI    =  sgg%observation(ii)%P(1:sgg%observation(ii)%nP)%ZI
               sgg%observation(2*tamaScrPrb/3+ii)%P(1:sgg%observation(2*tamaScrPrb/3+ii)%nP)%XE    =  sgg%observation(ii)%P(1:sgg%observation(ii)%nP)%XE
               sgg%observation(2*tamaScrPrb/3+ii)%P(1:sgg%observation(2*tamaScrPrb/3+ii)%nP)%YE    =  sgg%observation(ii)%P(1:sgg%observation(ii)%nP)%YE
               sgg%observation(2*tamaScrPrb/3+ii)%P(1:sgg%observation(2*tamaScrPrb/3+ii)%nP)%ZE    =  sgg%observation(ii)%P(1:sgg%observation(ii)%nP)%ZE
               !trancos
               sgg%observation(2*tamaScrPrb/3+ii)%P(1:sgg%observation(2*tamaScrPrb/3+ii)%nP)%Xtrancos =sgg%observation(ii)%P(1:sgg%observation(ii)%nP)%Xtrancos
               sgg%observation(2*tamaScrPrb/3+ii)%P(1:sgg%observation(2*tamaScrPrb/3+ii)%nP)%Ytrancos =sgg%observation(ii)%P(1:sgg%observation(ii)%nP)%Ytrancos
               sgg%observation(2*tamaScrPrb/3+ii)%P(1:sgg%observation(2*tamaScrPrb/3+ii)%nP)%Ztrancos =sgg%observation(ii)%P(1:sgg%observation(ii)%nP)%Ztrancos
               !fin trancos

!!!210618
!!!find 210618
            endif
            !DE LAS VolumicPROBLES
         END DO
!!!210618
         DO i = tamaScrPrb/3+1,tamaScrPrb
            ii = i + tamaSonda + tamaoldSONDA + tamaBloquePrb !Bug 040718
            if (sgg%observation(ii)%nP /=1) then
!para 040718
               WRITE (buff,*) '----> Volumic probe ii. np=',sgg%observation(ii)%nP
               CALL print11 (layoutnumber, buff)
               WRITE (buff,*) '----> Volumic probe ii. outputrequest=',trim(adjustl(sgg%observation(ii)%outputrequest))
               CALL print11 (layoutnumber, buff)
!fin para debugear
               WRITE (buff,*) 'Buggy error in Volumic probes. np/=1. , np=',sgg%observation(ii)%nP
               CALL STOPONERROR(layoutnumber,size,buff)
            endif
         end do
!!!find 210618

         !luego chequeo las sondas  si se van de memoria en observation.f90
         !        IF ((memo+sondas)*BuffObse*4 > MaxMemoryProbes) THEN
         !          WRITE (buff,*) 'Too much memory for the probes= ', (memo+sondas)*BuffObse*4, 'Probes= ', (memo+sondas), &
         !         & 'Either reduce the number of probes or recompile decreasing BuffObse ', BuffObse, 'or increasing ', MaxMemoryProbes
         !          CALL STOPONERROR(layoutnumber,size,buff)
         !        END IF

         !del if sgg%numberrequest
      END IF
      !las lineas goto 8 que sigue la comento a 27/10/14 porque "creo" que la informacion de shared es necesaria actualizarse
      !este bug aparece en bug_OLD221014_a400m_skindepth en Modelo.nfde
      !!!goto 8 !!!!?
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Update the number of the shared fields
      if (updateshared) then !!aqui se pierde mucho tiempo aniadido flag -noshared para evitarlo 040717
         WRITE (buff,*) 'INIT UPDATING SHARED INFO. This process may take time!'
         CALL print11 (layoutnumber, buff)
         WRITE (buff,*) 'Launch with -noshared to skip this process (just relevant for structured NIBC CFCs and Anisot.)'
         CALL print11 (layoutnumber, buff)
         DO i1 = 1, sgg%EShared%conta
            DO j1 = i1 + 1, sgg%EShared%conta
               IF ((sgg%Med(sgg%EShared%elem(i1)%SharedMed)%Priority == sgg%Med(sgg%EShared%elem(j1)%SharedMed)%Priority) .AND. &
               & (sgg%EShared%elem(i1)%field == sgg%EShared%elem(j1)%field) .AND. (sgg%EShared%elem(i1)%i == sgg%EShared%elem(j1)%i) &
               & .AND. (sgg%EShared%elem(i1)%j == sgg%EShared%elem(j1)%j) .AND. (sgg%EShared%elem(i1)%k == sgg%EShared%elem(j1)%k)) &
               & THEN
                  sgg%EShared%elem(i1)%times = sgg%EShared%elem(i1)%times + 1
                  sgg%EShared%elem(j1)%times = sgg%EShared%elem(j1)%times + 1
                  !!!!!!!!   write (17,*) i1,j1,sgg%EShared%elem(j1)%times,sgg%Med(sgg%EShared%elem(j1)%SharedMed)%Priority
                  if (sgg%EShared%elem(j1)%times > 4) then
                     write(buff,'(a,4i5)') 'WARNING: More than 4 media try to occupy ', sgg%EShared%elem(j1)%field,sgg%EShared%elem(j1)%i,   &
                        sgg%EShared%elem(j1)%j, sgg%EShared%elem(j1)%k
                     if ( ((.not.sgg%Med(sgg%EShared%elem(i1)%SharedMed)%is%thinwire).and.(.not.sgg%Med(sgg%EShared%elem(j1)%SharedMed)%is%thinwire)).and. &
                        ((.not.sgg%Med(sgg%EShared%elem(i1)%SharedMed)%is%SlantedWire).and.(.not.sgg%Med(sgg%EShared%elem(j1)%SharedMed)%is%SlantedWire)).or.verbose) CALL  WarnErrReport (buff)
                  endif
               END IF
            END DO
         END DO
         !Update the number of the shared fields
         DO i1 = 1, sgg%HShared%conta
            DO j1 = i1 + 1, sgg%HShared%conta
               IF ((sgg%Med(sgg%HShared%elem(i1)%SharedMed)%Priority == sgg%Med(sgg%HShared%elem(j1)%SharedMed)%Priority) .AND. &
               & (sgg%HShared%elem(i1)%field == sgg%HShared%elem(j1)%field) .AND. (sgg%HShared%elem(i1)%i == sgg%HShared%elem(j1)%i) &
               & .AND. (sgg%HShared%elem(i1)%j == sgg%HShared%elem(j1)%j) .AND. (sgg%HShared%elem(i1)%k == sgg%HShared%elem(j1)%k)) &
               & THEN
                  sgg%HShared%elem(i1)%times = sgg%HShared%elem(i1)%times + 1
                  sgg%HShared%elem(j1)%times = sgg%HShared%elem(j1)%times + 1
                  if (sgg%HShared%elem(j1)%times > 4) then
                     write(buff,'(a,4i5)') 'WARNING: More than 4 media try to occupy ', sgg%hShared%elem(j1)%field,sgg%HShared%elem(j1)%i,  &
                        sgg%HShared%elem(j1)%j, sgg%HShared%elem(j1)%k
                     if ( ((.not.sgg%Med(sgg%EShared%elem(i1)%SharedMed)%is%thinwire).and.(.not.sgg%Med(sgg%EShared%elem(j1)%SharedMed)%is%thinwire)).and. &
                        ((.not.sgg%Med(sgg%EShared%elem(i1)%SharedMed)%is%SlantedWire).and.(.not.sgg%Med(sgg%EShared%elem(j1)%SharedMed)%is%SlantedWire)).or.verbose) CALL  WarnErrReport (buff)
                  endif
               END IF
            END DO
         END DO
         !end shared
8        continue
         WRITE (buff,*) '[OK] END UPDATING SHARED INFO'
         CALL print11 (layoutnumber, buff)
      end if !del updateshared

      !PARA LA CAPA EXTRA 2013
      if (medioextra%exists) then
         CONTAMEDIA = CONTAMEDIA+1
         if  (MEDIOEXTRA%index /= contamedia) then !should be already done earlier
            CALL STOPONERROR(layoutnumber,size,'Bug in media count. ')
         endif
         MEDIOEXTRA%index=CONTAMEDIA
      endif
      !!!!!!!!!!!!!
      sgg%NumMedia = contamedia
      !el medio 0 no precisa compresion


      IF ((CLIPREGION)) THEN !ALLOW four cells OF AIR CELLS TO CLIP LARGE PROBLES WITH NO PROBLEMS WITH BOUNDARIES solo sin MPI
         do field=1,6
            DO K= sgg%Alloc(field)%ZI   , sgg%Alloc(field)%ZE
               DO J= sgg%Alloc(field)%YI   , sgg%Alloc(field)%YE
                  DO I= sgg%Alloc(field)%XI   , sgg%Alloc(field)%XE
                     if ( (i>=sinpml_FULLSIZE(field)%XI)  .AND.(i<=sinpml_FULLSIZE(field)%XI+4).OR. &
                        (i>=sinpml_FULLSIZE(field)%XE-4).AND.(i<=sinpml_FULLSIZE(field)%XE  ).OR. &
                        (J>=sinpml_FULLSIZE(field)%YI)  .AND.(j<=sinpml_FULLSIZE(field)%YI+4).OR. &
                        (J>=sinpml_FULLSIZE(field)%YE-4).AND.(j<=sinpml_FULLSIZE(field)%YE  ).OR. &
                        (K>=sinpml_FULLSIZE(field)%ZI)  .AND.(K<=sinpml_FULLSIZE(field)%ZI+4).OR. &
                        (K>=sinpml_FULLSIZE(field)%ZE-4).AND.(K<=sinpml_FULLSIZE(field)%ZE  )) THEN
                        select case (field)
                         case (iEx)
                           media%sggMIEX(I,J,K)=1
                         case (iEy)
                           media%sggMIEY(I,J,K)=1
                         case (iEz)
                           media%sggMIEZ(I,J,K)=1
                         case (iHx)
                           media%sggMiHX(I,J,K)=1
                         case (iHy)
                           media%sggMiHY(I,J,K)=1
                         case (iHz)
                           media%sggMiHZ(I,J,K)=1
                        end select
                     ENDIF
                  END DO
               END DO
            END DO
         end do !del field
      endif !del CLIPREGION

      !!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!fin clipeado
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CALL CreatePMLmatrix (layoutnumber, SIZE,sgg,media%sggMiEx,media%sggMiEy,media%sggMiEz,media%sggMiHx,media%sggMiHy,media%sggMiHz, SINPML_fullsize, fullsize, BoundingBox, sgg%Med, sgg%NumMedia, sgg%Border,MEDIOEXTRA)
      sgg%EndPMLMedia = sgg%NumMedia

      !
#ifdef CompileWithInt1
      IF (sgg%NumMedia > 127) THEN
         CLOSE (14)
         IF (sgg%NumMedia > 32767) THEN
            buff='Number of media>32767. Recompile with #define CompileWithInt4'
            CALL STOPONERROR(layoutnumber,size,buff)
         ELSE
            buff='Number of media>127. Recompile with #define CompileWithInt2'
            CALL STOPONERROR(layoutnumber,size,buff)
         END IF
      END IF
#endif
#ifdef CompileWithInt2
      IF (sgg%NumMedia > 32767) THEN
         CLOSE (14)
         buff='Number of media>32767. Recompile with #define CompileWithInt4'
         CALL STOPONERROR(layoutnumber,size,buff)
      END IF
#endif
#ifdef CompileWithInt4
      IF (sgg%NumMedia > 2.0e9) THEN
         CLOSE (14)
         buff='Number of media>2^31-1. Cannot continue. '
         CALL STOPONERROR(layoutnumber,size,buff)
      END IF
#endif
      !read the source files
      !
      call read_TIMEFRECTRANSFsourcefiles(simu_devia)
      !
      DO ii = 1, tamaSonda
         !Read the time normalization file
         !
         IF (sgg%observation(ii)%TRANSFER) THEN
            errnofile = .FALSE.
            INQUIRE (file=trim(adjustl(sgg%observation(ii)%FileNormalize)), EXIST=errnofile)
            IF ( .NOT. errnofile) THEN
               buff=trim(adjustl(sgg%observation(ii)%FileNormalize))//' DOES NOT EXIST'
               CALL STOPONERROR(layoutnumber,size,buff)
            END IF
         END IF
         !
      END DO
!!!!ajusta el flag lossy de los medios (261115) !aunque ya esta hecho por ahi arriba, lo rehago aqui
!!!!Ojo con que el usuario siempre ponga conductividad en !!compo para que esto no reviente
!!!mamma mia. comentado lo siguiente a 120123 para que las puestas a lossy thin wire conectado con resistor sean correctas. bug test_GGGbugresis_wire_stoch_foragasconbug
      !!!DO i = 1, sgg%NumMedia
      !!!      if ( (.not.(sgg%med(i)%is%PEC)).and.(sgg%med(i)%sigma >= 1e-4) ) then
      !!!         sgg%med(i)%is%lossy = .true.
      !!!      else
      !!!         sgg%med(i)%is%lossy = .false.
      !!!      endif
      !!!end do
      !!!fin 120123
!!!!!!
      !!!!!!!!!!do a final check if magnetic media are present
      !sgg jun'12 dejarlo siempre a .true. pq es lo mas seguro
      !sgg%thereAreMagneticMedia=.false.
      !sgg%thereArePMLMagneticMedia=.false.
      !medioespecial = .false.
      !do ii=1,sgg%NumMedia
      !    buff='PEC media can only have index 0'
      !    IF (sgg%Med(ii)%Is%PEC.or.sgg%Med(ii)%Is%Lumped) call STOPONERROR(layoutnumber,size,buff)
      !    medioespecial =medioespecial .or. &
      !                   sgg%Med(ii)%Is%EDispersive   .or. &
      !                   sgg%Med(ii)%Is%multiport     .or. &
      !                   sgg%Med(ii)%Is%AnisMultiport .or. &
      !                   sgg%Med(ii)%Is%Anisotropic   .or. &
      !                   sgg%Med(ii)%Is%ThinSlot      .or. &
      !                   sgg%Med(ii)%Is%MDispersive
      !end do
      !if (.not.medioespecial) then
      !    do ii=1,sgg%NumMedia
      !        if ((SGG%Med(ii)%mur >1.001).or.(SGG%Med(ii)%mur <0.999).or. &
      !            (abs(SGG%Med(ii)%sigmam) >1.0e-3_RKIND0)) then
      !            if (.not. sgg%Med(ii)%Is%PML ) then
      !                sgg%thereAreMagneticMedia=.true.
      !            else
      !                sgg%thereArePMLMagneticMedia=.true.
      !            endif
      !        endif
      !    end do
      !else
      sgg%thereAreMagneticMedia=.true.
      sgg%thereArePMLMagneticMedia=.true.
      !endif
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      RETURN

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine addConformalMedia(sgg, media, conformal_media, contamedia)
         type(sggfdtdinfo), intent(inout)    :: sgg
         type(media_matrices_t), intent(inout) :: media
         type(ConformalMedia_t), intent(in) :: conformal_media
         integer (kind=4), intent(inout) :: contamedia
         call addConformalEdgeMedia(sgg, media, conformal_media, contamedia)
         call addConformalFaceMedia(sgg, media, conformal_media, contamedia)
      end subroutine

      subroutine addConformalFaceMedia(sgg, media, conformal_media, contamedia)
         type (SGGFDTDINFO), intent(INOUT)    :: sgg
         type(media_matrices_t), intent(inout) :: media
         type(ConformalMedia_t), intent(in) :: conformal_media
         integer (kind=4), intent(inout) :: contamedia
         integer (kind=4) :: cell_i, cell_j, cell_k
         integer :: j, k
         do j = 1, conformal_media%n_faces_media
            contamedia = contamedia + 1
            sgg%Med(contamedia)%Is%ConformalPEC = .TRUE.
            sgg%Med(contamedia)%Is%Needed = .TRUE.
            sgg%Med(contamedia)%Priority = prior_PEC
            sgg%Med(contamedia)%Epr = this%mats%mats(1)%eps / Eps0
            sgg%Med(contamedia)%Sigma = 1.0e29_RKIND
            sgg%Med(contamedia)%Mur = conformal_media%face_media(j)%ratio * this%mats%mats(1)%mu / Mu0
            sgg%Med(contamedia)%SigmaM = 0.0_RKIND
            ! funcion al healer para tener en cuenta la prioridad
            do k = 1, conformal_media%face_media(j)%size
               cell_i = conformal_media%face_media(j)%faces(k)%cell(1)
               cell_j = conformal_media%face_media(j)%faces(k)%cell(2)
               cell_k = conformal_media%face_media(j)%faces(k)%cell(3)
               select case(conformal_media%face_media(j)%faces(k)%direction)
               case(F_X)
                  media%sggMiHx(cell_i, cell_j, cell_k) = contamedia
               case(F_Y)
                  media%sggMiHy(cell_i, cell_j, cell_k) = contamedia
               case(F_Z)
                  media%sggMiHz(cell_i, cell_j, cell_k) = contamedia
               end select
            end do
         end do
   
      end subroutine

      subroutine addConformalEdgeMedia(sgg, media, conformal_media, contamedia)
         type (SGGFDTDINFO), intent(INOUT)    :: sgg
         type(media_matrices_t), intent(inout) :: media
         type(ConformalMedia_t), intent(in) :: conformal_media
         integer (kind=4), intent(inout) :: contamedia
         integer (kind=4) :: edge_media
         integer (kind=4) :: cell_i, cell_j, cell_k
         integer :: j, k
         do j = 1, conformal_media%n_edges_media
            if (conformal_media%edge_media(j)%ratio /= 0) then 
               contamedia = contamedia + 1
               sgg%Med(contamedia)%Is%ConformalPEC = .TRUE.
               sgg%Med(contamedia)%Is%Needed = .TRUE.
               sgg%Med(contamedia)%Priority = prior_PEC
               sgg%Med(contamedia)%Epr = (this%mats%mats(1)%eps / conformal_media%edge_media(j)%ratio ) / Eps0
               sgg%Med(contamedia)%Sigma = 1.0e29_RKIND
               sgg%Med(contamedia)%Mur = this%mats%mats(1)%mu / Mu0
               sgg%Med(contamedia)%SigmaM = 0.0_RKIND
               edge_media = contamedia
            else
               edge_media = 0
            end if
            ! funcion al healer para tener en cuenta la prioridad
            do k = 1, conformal_media%edge_media(j)%size
               cell_i = conformal_media%edge_media(j)%edges(k)%cell(1)
               cell_j = conformal_media%edge_media(j)%edges(k)%cell(2)
               cell_k = conformal_media%edge_media(j)%edges(k)%cell(3)
               select case(conformal_media%edge_media(j)%edges(k)%direction)
               case(E_X)
                  media%sggMiEx(cell_i, cell_j, cell_k) = edge_media
               case(E_Y)
                  media%sggMiEy(cell_i, cell_j, cell_k) = edge_media
               case(E_Z)
                  media%sggMiEz(cell_i, cell_j, cell_k) = edge_media
               end select
            end do
         end do
      end subroutine

      subroutine read_TIMEFRECTRANSFsourcefiles(simu_devia)
         logical :: simu_devia
         real (kind=RKIND) :: deviafactor_Multiplier,unillo,tiempoant
         if (simu_devia) then
            deviafactor_Multiplier=1.0_RKIND !correccion de 160619. gestiono la iluminacion o no con la variable simu_devia, para que postprocese bien los ficheros devia_ !deprecar deviafactor_Multiplier algun dia
         ELSE
            deviafactor_Multiplier=1.0_RKIND
         ENDIF
!!!! SOURCES IN WIRES
         maxSourceValue=0.0_RKIND
         minSpaceStep=min(min(minval(sgg%dx),minval(sgg%dy)),minval(sgg%dz))
         DO i = 1, sgg%NumMedia
            IF (sgg%Med(i)%Is%ThinWire) THEN
               IF (sgg%Med(i)%wire(1)%VsourceExists) THEN
                  DO CONTAVOLT=1,sgg%Med(i)%wire(1)%NUMVOLTAGESOURCES
                     errnofile = .FALSE.
                     INQUIRE (file=trim(adjustl(sgg%Med(i)%wire(1)%VSource(CONTAVOLT)%fichero%NAME)), EXIST=errnofile)
                     IF ( .NOT. errnofile) THEN
                        buff=trim(adjustl(sgg%Med(i)%wire(1)%VSource(CONTAVOLT)%fichero%name))//' DOES NOT EXIST'
                        CALL STOPONERROR(layoutnumber,size,buff)
                     END IF
                     OPEN (15, file=trim(adjustl(sgg%Med(i)%wire(1)%VSource(CONTAVOLT)%fichero%NAME)),action='read')
                     READ (15,*) tiempo1, field1
                     READ (15,*) tiempo2, field2
                     sgg%Med(i)%wire(1)%VSource(CONTAVOLT)%fichero%deltaSamples = tiempo2 - tiempo1
                     nsurfs = 3
                     !problemas con multivac
                     ! while (.not.eof(15))
                     DO
                        READ (15,*, end=77) tiempo1, field1
                        if (field1/minspacestep > maxSourceValue) maxSourceValue=field1/minspacestep
                        nsurfs = nsurfs + 1
                     END DO
77                   CONTINUE
                     CLOSE (15)
                     numus = nsurfs - 2
                     sgg%Med(i)%wire(1)%VSource(CONTAVOLT)%fichero%NumSamples = numus
                     ALLOCATE (sgg%Med(i)%wire(1)%VSource(CONTAVOLT)%fichero%Samples(0:numus))
                     OPEN (15, file=trim(adjustl(sgg%Med(i)%wire(1)%VSource(CONTAVOLT)%fichero%NAME)),action='read')
                     DO k = 0, numus
                        tiempoant=tiempo1
                        READ (15,*) tiempo1, sgg%Med(i)%wire(1)%VSource(CONTAVOLT)%fichero%Samples(k)
                        sgg%Med(i)%wire(1)%VSource(CONTAVOLT)%fichero%Samples(k) = sgg%Med(i)%wire(1)%VSource(CONTAVOLT)%fichero%Samples(k) * &
                        & sgg%Med(i)%wire(1)%VSource(CONTAVOLT)%Multiplier * &
                           deviafactor_Multiplier
                        !!evitar sampleos no uniformes
                        if ((k>1).and.(k<numus-1)) then
                           unillo=(tiempo1-tiempoant)/sgg%Med(i)%wire(1)%VSource(CONTAVOLT)%fichero%deltaSamples
                           if ((unillo<0.9).or.(unillo>1.0_RKIND/0.9)) then
                              if (2.0_RKIND*(tiempo1-tiempoant)/(tiempo1+tiempoant)>1e-6_RKIND) then !a tiempos muy altos ignoro el redondeo
                                 buff=trim(adjustl(sgg%Med(i)%wire(1)%VSource(CONTAVOLT)%fichero%NAME))//' not uniformly sampled. Relaunch with -ignoresamplingerrors to ignore it.'
                                 if (.not.ignoresamplingerrors) then
                                    CLOSE(15)
                                    CALL STOPONERROR(layoutnumber,size,buff)
                                 endif
                              endif
                           endif
                        endif
                     END DO
                     CLOSE (15)
                     sgg%Med(i)%wire(1)%VSource(CONTAVOLT)%fichero%NumSamples = numus
                  END DO
               END IF
            END IF
            IF (sgg%Med(i)%Is%ThinWire) THEN
               IF (sgg%Med(i)%wire(1)%IsourceExists) THEN
                  DO CONTACURR=1,sgg%Med(i)%wire(1)%NUMCURRENTSOURCES
                     errnofile = .FALSE.
                     INQUIRE (file=trim(adjustl(sgg%Med(i)%wire(1)%ISource(CONTACURR)%fichero%NAME)), EXIST=errnofile)
                     IF ( .NOT. errnofile) THEN
                        buff=trim(adjustl(sgg%Med(i)%wire(1)%ISource(CONTACURR)%fichero%name))//' DOES NOT EXIST'
                        CALL STOPONERROR(layoutnumber,size,buff)
                     END IF
                     OPEN (15, file=trim(adjustl(sgg%Med(i)%wire(1)%ISource(CONTACURR)%fichero%NAME)),action='read')
                     READ (15,*) tiempo1, field1
                     READ (15,*) tiempo2, field2
                     sgg%Med(i)%wire(1)%ISource(CONTACURR)%fichero%deltaSamples = tiempo2 - tiempo1
                     nsurfs = 3
                     !problemas con multivac
                     ! while (.not.eof(15))
                     DO
                        READ (15,*, end=79) tiempo1, field1
                        if (field1/minspacestep**2.0_RKIND > maxSourceValue) maxSourceValue=field1/minspacestep**2.0_RKIND !aqui no tengo feeling, pero estas fuentes no se usan !!? repensar
                        nsurfs = nsurfs + 1
                     END DO
79                   CONTINUE
                     CLOSE (15)
                     numus = nsurfs - 2
                     sgg%Med(i)%wire(1)%ISource(CONTACURR)%fichero%NumSamples = numus
                     ALLOCATE (sgg%Med(i)%wire(1)%ISource(CONTACURR)%fichero%Samples(0:numus))
                     OPEN (15, file=trim(adjustl(sgg%Med(i)%wire(1)%ISource(CONTACURR)%fichero%NAME)),action='read')
                     DO k = 0, numus
                        tiempoant=tiempo1
                        READ (15,*) tiempo1, sgg%Med(i)%wire(1)%ISource(CONTACURR)%fichero%Samples(k)
                        sgg%Med(i)%wire(1)%ISource(CONTACURR)%fichero%Samples(k) = sgg%Med(i)%wire(1)%ISource(CONTACURR)%fichero%Samples(k) * &
                        & sgg%Med(i)%wire(1)%ISource(CONTACURR)%Multiplier * &
                           deviafactor_Multiplier
                        !!evitar sampleos no uniformes
                        if ((k>1).and.(k<numus-1)) then
                           unillo=(tiempo1-tiempoant)/sgg%Med(i)%wire(1)%ISource(CONTACURR)%fichero%deltaSamples
                           if ((unillo<0.9).or.(unillo>1.0_RKIND/0.9)) then
                              if (2.0_RKIND*(tiempo1-tiempoant)/(tiempo1+tiempoant)>1e-6_RKIND) then !a tiempos muy altos ignoro el redondeo
                                 buff=trim(adjustl(sgg%Med(i)%wire(1)%ISource(CONTACURR)%fichero%NAME))//' not uniformly sampled. Relaunch with -ignoresamplingerrors to ignore it.'
                                 if (.not.ignoresamplingerrors) then
                                    CLOSE(15)
                                    CALL STOPONERROR(layoutnumber,size,buff)
                                 endif
                              endif
                           endif
                        endif
                        !!
                     END DO
                     CLOSE (15)
                     sgg%Med(i)%wire(1)%ISource(CONTACURR)%fichero%NumSamples = numus
                  END DO
               END IF
            END IF
            !
            IF (sgg%Med(i)%Is%SlantedWire) THEN
               DO j = 1,sgg%Med(i)%SlantedWire(1)%numNodes
                  IF (sgg%Med(i)%SlantedWire(1)%nodes(j)%VsourceExists) then
                     errnofile = .FALSE.
                     INQUIRE (file=trim(adjustl(sgg%Med(i)%SlantedWire(1)%nodes(j)%Vsource%fichero%NAME)), EXIST=errnofile)
                     IF ( .NOT. errnofile) THEN
                        buff=trim(adjustl(sgg%Med(i)%SlantedWire(1)%nodes(j)%Vsource%fichero%name))//' DOES NOT EXIST'
                        CALL STOPONERROR(layoutnumber,size,buff)
                     END IF
                     OPEN (15, file=trim(adjustl(sgg%Med(i)%SlantedWire(1)%nodes(j)%Vsource%fichero%NAME)),action='read')
                     READ (15,*) tiempo1, field1
                     READ (15,*) tiempo2, field2
                     sgg%Med(i)%SlantedWire(1)%nodes(j)%Vsource%fichero%deltaSamples = tiempo2 - tiempo1
                     nsurfs = 3
                     !problemas con multivac
                     ! while (.not.eof(15))
                     DO
                        READ (15,*, end=179) tiempo1, field1
                        if (field1/minspacestep > maxSourceValue) maxSourceValue=field1/minspacestep
                        nsurfs = nsurfs + 1
                     END DO
179                  CONTINUE
                     CLOSE (15)
                     numus = nsurfs - 2
                     sgg%Med(i)%SlantedWire(1)%nodes(j)%Vsource%fichero%NumSamples = numus
                     ALLOCATE (sgg%Med(i)%SlantedWire(1)%nodes(j)%Vsource%fichero%Samples(0:numus))
                     OPEN (15, file=trim(adjustl(sgg%Med(i)%SlantedWire(1)%nodes(j)%Vsource%fichero%NAME)),action='read')
                     DO k = 0, numus
                        tiempoant=tiempo1
                        READ (15,*) tiempo1, sgg%Med(i)%SlantedWire(1)%nodes(j)%Vsource%fichero%Samples(k)
                        sgg%Med(i)%SlantedWire(1)%nodes(j)%Vsource%fichero%Samples(k) = sgg%Med(i)%SlantedWire(1)%nodes(j)%Vsource%fichero%Samples(k) * &
                           sgg%Med(i)%SlantedWire(1)%nodes(j)%Vsource%Multiplier * &
                           deviafactor_Multiplier
                        !!evitar sampleos no uniformes
                        if ((k>1).and.(k<numus-1)) then
                           unillo=(tiempo1-tiempoant)/sgg%Med(i)%SlantedWire(1)%nodes(j)%Vsource%fichero%deltaSamples
                           if ((unillo<0.9).or.(unillo>1.0_RKIND/0.9)) then
                              if (2.0_RKIND*(tiempo1-tiempoant)/(tiempo1+tiempoant)>1e-6_RKIND) then !a tiempos muy altos ignoro el redondeo
                                 buff=trim(adjustl(sgg%Med(i)%SlantedWire(1)%nodes(j)%Vsource%fichero%NAME))//' not uniformly sampled. Relaunch with -ignoresamplingerrors to ignore it.'
                                 if (.not.ignoresamplingerrors) then
                                    CLOSE(15)
                                    CALL STOPONERROR(layoutnumber,size,buff)
                                 endif
                              endif
                           endif
                        endif
                        !!
                     END DO
                     CLOSE (15)
                     sgg%Med(i)%SlantedWire(1)%nodes(j)%Vsource%fichero%NumSamples = numus
                  END IF
                  IF (sgg%Med(i)%SlantedWire(1)%nodes(j)%IsourceExists) then
                     errnofile = .FALSE.
                     INQUIRE (file=trim(adjustl(sgg%Med(i)%SlantedWire(1)%nodes(j)%Isource%fichero%NAME)), EXIST=errnofile)
                     IF ( .NOT. errnofile) THEN
                        buff=trim(adjustl(sgg%Med(i)%SlantedWire(1)%nodes(j)%Isource%fichero%name))//' DOES NOT EXIST'
                        CALL STOPONERROR(layoutnumber,size,buff)
                     END IF
                     OPEN (15, file=trim(adjustl(sgg%Med(i)%SlantedWire(1)%nodes(j)%Isource%fichero%NAME)),action='read')
                     READ (15,*) tiempo1, field1
                     READ (15,*) tiempo2, field2
                     sgg%Med(i)%SlantedWire(1)%nodes(j)%Isource%fichero%deltaSamples = tiempo2 - tiempo1
                     nsurfs = 3
                     !problemas con multivac
                     ! while (.not.eof(15))
                     DO
                        READ (15,*, end=279) tiempo1, field1
                        if (field1/minspacestep**2.0_RKIND > maxSourceValue) maxSourceValue=field1/minspacestep**2.0_RKIND !aqui no tengo feeling, pero estas fuentes no se usan !!? repensar
                        nsurfs = nsurfs + 1
                     END DO
279                  CONTINUE
                     CLOSE (15)
                     numus = nsurfs - 2
                     sgg%Med(i)%SlantedWire(1)%nodes(j)%Isource%fichero%NumSamples = numus
                     ALLOCATE (sgg%Med(i)%SlantedWire(1)%nodes(j)%Isource%fichero%Samples(0:numus))
                     OPEN (15, file=trim(adjustl(sgg%Med(i)%SlantedWire(1)%nodes(j)%Isource%fichero%NAME)),action='read')
                     DO k = 0, numus
                        tiempoant=tiempo1
                        READ (15,*) tiempo1, sgg%Med(i)%SlantedWire(1)%nodes(j)%Isource%fichero%Samples(k)
                        sgg%Med(i)%SlantedWire(1)%nodes(j)%Isource%fichero%Samples(k) = sgg%Med(i)%SlantedWire(1)%nodes(j)%Isource%fichero%Samples(k) * &
                           sgg%Med(i)%SlantedWire(1)%nodes(j)%Isource%Multiplier * &
                           deviafactor_Multiplier
                        !!evitar sampleos no uniformes
                        if ((k>1).and.(k<numus-1)) then
                           unillo=(tiempo1-tiempoant)/sgg%Med(i)%SlantedWire(1)%nodes(j)%Isource%fichero%deltaSamples
                           if ((unillo<0.9).or.(unillo>1.0_RKIND/0.9)) then
                              if (2.0_RKIND*(tiempo1-tiempoant)/(tiempo1+tiempoant)>1e-6_RKIND) then !a tiempos muy altos ignoro el redondeo
                                 buff=trim(adjustl(sgg%Med(i)%SlantedWire(1)%nodes(j)%Isource%fichero%NAME))//' not uniformly sampled. Relaunch with -ignoresamplingerrors to ignore it.'
                                 if (.not.ignoresamplingerrors) then
                                    CLOSE(15)
                                    CALL STOPONERROR(layoutnumber,size,buff)
                                 endif
                              endif
                           endif
                        endif
                        !!
                     END DO
                     CLOSE (15)
                     sgg%Med(i)%SlantedWire(1)%nodes(j)%Isource%fichero%NumSamples = numus
                  END IF
               END DO
            END IF
         END DO
         !nodal sources
         DO j = 1, sgg%NumNodalSources
            !Read the time evolution file
            !
            if (.not.sgg%NodalSource(j)%IsInitialValue) then
               errnofile = .FALSE.
               INQUIRE (file=trim(adjustl(sgg%NodalSource(j)%fichero%NAME)), EXIST=errnofile)
               IF ( .NOT. errnofile) THEN
                  buff=trim(adjustl(sgg%NodalSource(j)%fichero%name))//' DOES NOT EXIST'
                  CALL STOPONERROR(layoutnumber,size,buff)
               END IF
               OPEN (15, file=trim(adjustl(sgg%NodalSource(j)%fichero%NAME)),action='read')
               READ (15,*) tiempo1, field1
               READ (15,*) tiempo2, field2
               sgg%NodalSource(j)%fichero%deltaSamples = tiempo2 - tiempo1
               nsurfs = 3
               DO
                  READ (15,*, end=78) tiempo1, field1
                  if (sgg%NodalSource(j)%isElec) then
                     if (field1 > maxSourceValue) maxSourceValue=field1
                  else
                     if (field1*zvac > maxSourceValue) maxSourceValue=field1*zvac
                  endif
                  nsurfs = nsurfs + 1
               END DO
78             CONTINUE
               CLOSE (15)
               numus = nsurfs - 2
               ALLOCATE (sgg%NodalSource(j)%fichero%Samples(0:numus))
               OPEN (15, file=trim(adjustl(sgg%NodalSource(j)%fichero%NAME)),action='read')
               DO k = 0, numus
                  tiempoant=tiempo1
                  READ (15,*) tiempo1, sgg%NodalSource(j)%fichero%Samples(k)
                  sgg%NodalSource(j)%fichero%Samples(k) = sgg%NodalSource(j)%fichero%Samples(k) * deviafactor_Multiplier
                  !!evitar sampleos no uniformes
                  if ((k>1).and.(k<numus-1)) then
                     unillo=(tiempo1-tiempoant)/sgg%NodalSource(j)%fichero%deltaSamples
                     if ((unillo<0.9).or.(unillo>1.0_RKIND/0.9)) then
                        if (2.0_RKIND*(tiempo1-tiempoant)/(tiempo1+tiempoant)>1e-6_RKIND) then !a tiempos muy altos ignoro el redondeo
                           buff=trim(adjustl(sgg%NodalSource(j)%fichero%NAME))//' not uniformly sampled. Relaunch with -ignoresamplingerrors to ignore it.'
                           if (.not.ignoresamplingerrors) then
                              CLOSE(15)
                              CALL STOPONERROR(layoutnumber,size,buff)
                           endif
                        endif
                     endif
                  endif
                  !!
               END DO
               CLOSE (15)
               sgg%NodalSource(j)%fichero%NumSamples = numus
            else !es un initialvalue que no precisa fichero. incializar trivialmente
               numus = 0
               sgg%NodalSource(j)%fichero%deltaSamples = 1.0_RKIND !must be non-null
               ALLOCATE (sgg%NodalSource(j)%fichero%Samples(0:numus))
               DO k = 0, numus
                  sgg%NodalSource(j)%fichero%Samples(k) = 1.0_RKIND * deviafactor_Multiplier
               END DO
               sgg%NodalSource(j)%fichero%NumSamples = numus
            endif
         END DO
         !Plane wave sources
         DO j = 1, sgg%NumPlaneWaves
            !Read the time evolution file
            !
            errnofile = .FALSE.
            INQUIRE (file=trim(adjustl(sgg%PlaneWave(j)%fichero%NAME)), EXIST=errnofile)
            IF ( .NOT. errnofile) THEN
               buff=trim(adjustl(sgg%PlaneWave(j)%fichero%name))//' DOES NOT EXIST'
               CALL STOPONERROR(layoutnumber,size,buff)
            END IF
            OPEN (15, file=trim(adjustl(sgg%PlaneWave(j)%fichero%NAME)),action='read')
            READ (15,*) tiempo1, field1
            READ (15,*) tiempo2, field2
            sgg%PlaneWave(j)%fichero%deltaSamples = tiempo2 - tiempo1
            nsurfs = 3
            DO
               READ (15,*, end=98) tiempo1, field1
               if (field1 > maxSourceValue) maxSourceValue=field1
               nsurfs = nsurfs + 1
            END DO
98          CONTINUE
            CLOSE (15)
            numus = nsurfs - 2
            ALLOCATE (sgg%PlaneWave(j)%fichero%Samples(0:numus))
            OPEN (15, file=trim(adjustl(sgg%PlaneWave(j)%fichero%NAME)),action='read')
            DO k = 0, numus
               tiempoant=tiempo1
               READ (15,*) tiempo1, sgg%PlaneWave(j)%fichero%Samples(k)
               sgg%PlaneWave(j)%fichero%Samples(k) = sgg%PlaneWave(j)%fichero%Samples(k) * deviafactor_Multiplier
               !!evitar sampleos no uniformes
               if ((k>1).and.(k<numus-1)) then
                  unillo=(tiempo1-tiempoant)/sgg%PlaneWave(j)%fichero%deltaSamples
                  if ((unillo<0.9).or.(unillo>1.0_RKIND/0.9)) then
                     if (2.0_RKIND*(tiempo1-tiempoant)/(tiempo1+tiempoant)>1e-6_RKIND) then !a tiempos muy altos ignoro el redondeo
                        buff=trim(adjustl(sgg%PlaneWave(j)%fichero%NAME))//' not uniformly sampled. Relaunch with -ignoresamplingerrors to ignore it.'
                        if (.not.ignoresamplingerrors) then
                           CLOSE(15)
                           CALL STOPONERROR(layoutnumber,size,buff)
                        endif
                     endif
                  endif
               endif
               !!
            END DO
            CLOSE (15)
            sgg%PlaneWave(j)%fichero%NumSamples = numus
         END DO
         return
      end subroutine read_TIMEFRECTRANSFsourcefiles

!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!


      subroutine asignawiredisper(disp, file)
         type (WireDispersiveParams_t), intent(inout) :: disp
         CHARACTER (LEN=BUFSIZE), intent(in)        :: file

         integer (kind=4)  :: i, numPoles
         real (kind=RKIND) :: valreal, valimag

         open (1712, file=file, action='READ')
         read (1712, *) numPoles

         disp%NumPoles = numPoles
         ALLOCATE (disp%res(1:numPoles))
         ALLOCATE (disp%p  (1:numPoles))
         do i = 1, numPoles
            read (1712, *) valreal, valimag
            disp%res(i) = CMPLX(valreal, valimag, CKIND)
         end do
         do i = 1, numPoles
            read (1712, *) valreal, valimag
            disp%p(i) = CMPLX(valreal, valimag, CKIND)
         end do
         read (1712, *) valreal, valimag
         disp%d = CMPLX(valreal, valimag, CKIND)
         read (1712, *) valreal, valimag
         disp%e = CMPLX(valreal, valimag, CKIND)

         return
      end subroutine asignawiredisper

      subroutine asignadisper(fdgeom)
         TYPE (FreqDepenMaterial), POINTER :: fdgeom

         IF (fdgeom%l+fdgeom%LM /=0 ) THEN
            BUFF='ERROR: SECOND ORDER DISPERSIVE MEDIA UNSUPPORTED. TRANSLATE THEM TO FIRST ORDER ()'
            CALL WarnErrReport (buff,.TRUE.)
         endif
         !
         ALLOCATE (sgg%Med(contamedia)%EDispersive(1))
         ALLOCATE (sgg%Med(contamedia)%MDispersive(1))
         sgg%Med(contamedia)%Priority = prior_FDB
         sgg%Med(contamedia)%Epr = fdgeom%eps11 / Eps0
         sgg%Med(contamedia)%Sigma = fdgeom%Sigma11
         sgg%Med(contamedia)%Mur = fdgeom%mu11 / Mu0
         sgg%Med(contamedia)%SigmaM = fdgeom%SigmaM11
         !electric dispersion
         sgg%Med(contamedia)%Edispersive(1)%eps11=       fdgeom%eps11
         sgg%Med(contamedia)%Edispersive(1)%eps12=       fdgeom%eps12
         sgg%Med(contamedia)%Edispersive(1)%eps13=       fdgeom%eps13
         sgg%Med(contamedia)%Edispersive(1)%eps22=       fdgeom%eps22
         sgg%Med(contamedia)%Edispersive(1)%eps23=       fdgeom%eps23
         sgg%Med(contamedia)%Edispersive(1)%eps33=       fdgeom%eps33

         sgg%Med(contamedia)%Edispersive(1)%mu11=       fdgeom%MU11
         sgg%Med(contamedia)%Edispersive(1)%mu12=       fdgeom%MU12
         sgg%Med(contamedia)%Edispersive(1)%mu13=       fdgeom%MU13
         sgg%Med(contamedia)%Edispersive(1)%mu22=       fdgeom%MU22
         sgg%Med(contamedia)%Edispersive(1)%mu23=       fdgeom%MU23
         sgg%Med(contamedia)%Edispersive(1)%mu33=       fdgeom%MU33

         sgg%Med(contamedia)%EDispersive(1)%SIGMA11=  fdgeom%SIGMA11
         sgg%Med(contamedia)%EDispersive(1)%SIGMA12=  fdgeom%SIGMA12
         sgg%Med(contamedia)%EDispersive(1)%SIGMA13=  fdgeom%SIGMA13
         sgg%Med(contamedia)%EDispersive(1)%SIGMA22=  fdgeom%SIGMA22
         sgg%Med(contamedia)%EDispersive(1)%SIGMA23=  fdgeom%SIGMA23
         sgg%Med(contamedia)%EDispersive(1)%SIGMA33=  fdgeom%SIGMA33

         sgg%Med(contamedia)%EDispersive(1)%SIGMAM11=fdgeom%SIGMAM11
         sgg%Med(contamedia)%EDispersive(1)%SIGMAM12=fdgeom%SIGMAM12
         sgg%Med(contamedia)%EDispersive(1)%SIGMAM13=fdgeom%SIGMAM13
         sgg%Med(contamedia)%EDispersive(1)%SIGMAM22=fdgeom%SIGMAM22
         sgg%Med(contamedia)%EDispersive(1)%SIGMAM23=fdgeom%SIGMAM23
         sgg%Med(contamedia)%EDispersive(1)%SIGMAM33=fdgeom%SIGMAM33

         sgg%Med(contamedia)%Mdispersive(1)%eps11=        sgg%Med(contamedia)%Edispersive(1)%eps11
         sgg%Med(contamedia)%Mdispersive(1)%eps12=        sgg%Med(contamedia)%Edispersive(1)%eps12
         sgg%Med(contamedia)%Mdispersive(1)%eps13=        sgg%Med(contamedia)%Edispersive(1)%eps13
         sgg%Med(contamedia)%Mdispersive(1)%eps22=        sgg%Med(contamedia)%Edispersive(1)%eps22
         sgg%Med(contamedia)%Mdispersive(1)%eps23=        sgg%Med(contamedia)%Edispersive(1)%eps23
         sgg%Med(contamedia)%Mdispersive(1)%eps33=        sgg%Med(contamedia)%Edispersive(1)%eps33

         sgg%Med(contamedia)%Mdispersive(1)%mu11=        sgg%Med(contamedia)%Edispersive(1)%mu11
         sgg%Med(contamedia)%Mdispersive(1)%mu12=        sgg%Med(contamedia)%Edispersive(1)%mu12
         sgg%Med(contamedia)%Mdispersive(1)%mu13=        sgg%Med(contamedia)%Edispersive(1)%mu13
         sgg%Med(contamedia)%Mdispersive(1)%mu22=        sgg%Med(contamedia)%Edispersive(1)%mu22
         sgg%Med(contamedia)%Mdispersive(1)%mu23=        sgg%Med(contamedia)%Edispersive(1)%mu23
         sgg%Med(contamedia)%Mdispersive(1)%mu33=        sgg%Med(contamedia)%Edispersive(1)%mu33

         sgg%Med(contamedia)%MDISPERSIVE(1)%SIGMA11=    sgg%Med(contamedia)%EDispersive(1)%SIGMA11
         sgg%Med(contamedia)%MDISPERSIVE(1)%SIGMA12=    sgg%Med(contamedia)%EDispersive(1)%SIGMA12
         sgg%Med(contamedia)%MDISPERSIVE(1)%SIGMA13=    sgg%Med(contamedia)%EDispersive(1)%SIGMA13
         sgg%Med(contamedia)%MDISPERSIVE(1)%SIGMA22=    sgg%Med(contamedia)%EDispersive(1)%SIGMA22
         sgg%Med(contamedia)%MDISPERSIVE(1)%SIGMA23=    sgg%Med(contamedia)%EDispersive(1)%SIGMA23
         sgg%Med(contamedia)%MDISPERSIVE(1)%SIGMA33=    sgg%Med(contamedia)%EDispersive(1)%SIGMA33

         sgg%Med(contamedia)%MDISPERSIVE(1)%SIGMAM11=  sgg%Med(contamedia)%EDispersive(1)%SIGMAM11
         sgg%Med(contamedia)%MDISPERSIVE(1)%SIGMAM12=  sgg%Med(contamedia)%EDispersive(1)%SIGMAM12
         sgg%Med(contamedia)%MDISPERSIVE(1)%SIGMAM13=  sgg%Med(contamedia)%EDispersive(1)%SIGMAM13
         sgg%Med(contamedia)%MDISPERSIVE(1)%SIGMAM22=  sgg%Med(contamedia)%EDispersive(1)%SIGMAM22
         sgg%Med(contamedia)%MDISPERSIVE(1)%SIGMAM23=  sgg%Med(contamedia)%EDispersive(1)%SIGMAM23
         sgg%Med(contamedia)%MDISPERSIVE(1)%SIGMAM33=  sgg%Med(contamedia)%EDispersive(1)%SIGMAM33
         !
         !los de primer orden solo. Los de segundo no juegan
         sgg%Med(contamedia)%EDispersive(1)%NumPolRes11 = fdgeom%k11  !+ fdgeom%l
         sgg%Med(contamedia)%EDispersive(1)%NumPolRes12 = fdgeom%k12
         sgg%Med(contamedia)%EDispersive(1)%NumPolRes13 = fdgeom%k13
         sgg%Med(contamedia)%EDispersive(1)%NumPolRes22 = fdgeom%k22
         sgg%Med(contamedia)%EDispersive(1)%NumPolRes23 = fdgeom%k23
         sgg%Med(contamedia)%EDispersive(1)%NumPolRes33 = fdgeom%k33
         !magnetic
         sgg%Med(contamedia)%MDispersive(1)%NumPolRes11 = fdgeom%KM11 ! + fdgeom%LM
         sgg%Med(contamedia)%MDispersive(1)%NumPolRes12 = fdgeom%KM12
         sgg%Med(contamedia)%MDispersive(1)%NumPolRes13 = fdgeom%KM13
         sgg%Med(contamedia)%MDispersive(1)%NumPolRes22 = fdgeom%KM22
         sgg%Med(contamedia)%MDispersive(1)%NumPolRes23 = fdgeom%KM23
         sgg%Med(contamedia)%MDispersive(1)%NumPolRes33 = fdgeom%KM33
         !
         IF (sgg%Med(contamedia)%EDispersive(1)%NumPolRes11 /= 0) THEN
            sgg%Med(contamedia)%Is%EDispersive = .TRUE.
            sgg%Med(contamedia)%Is%EDispersiveANIS = .FALSE.
            sgg%Med(contamedia)%Is%Dielectric = .FALSE.
         END IF
         IF (sgg%Med(contamedia)%EDispersive(1)%NumPolRes12+sgg%Med(contamedia)%EDispersive(1)%NumPolRes13+ &
            sgg%Med(contamedia)%EDispersive(1)%NumPolRes22+sgg%Med(contamedia)%EDispersive(1)%NumPolRes23+ &
            sgg%Med(contamedia)%EDispersive(1)%NumPolRes33 /= 0) THEN
            sgg%Med(contamedia)%Is%EDispersive = .TRUE.
            sgg%Med(contamedia)%Is%EDispersiveAnis = .TRUE.
            print *, "Error: anisotropic dispersive unsupported"
            stop
            sgg%Med(contamedia)%Is%Dielectric = .FALSE.
         END IF
         IF (sgg%Med(contamedia)%MDispersive(1)%NumPolRes11 /= 0) THEN
            sgg%Med(contamedia)%Is%MDispersive = .TRUE.
            sgg%Med(contamedia)%Is%Dielectric = .FALSE.
         END IF
         IF (sgg%Med(contamedia)%MDISPERSIVE(1)%NumPolRes12+sgg%Med(contamedia)%MDISPERSIVE(1)%NumPolRes13+ &
            sgg%Med(contamedia)%MDISPERSIVE(1)%NumPolRes22+sgg%Med(contamedia)%MDISPERSIVE(1)%NumPolRes23+ &
            sgg%Med(contamedia)%MDISPERSIVE(1)%NumPolRes33 /= 0) THEN
            sgg%Med(contamedia)%Is%MDISPERSIVE = .TRUE.
            sgg%Med(contamedia)%Is%MDISPERSIVEAnis = .TRUE.
            print *, "Error: anisotropic dispersive unsupported"
            stop
            sgg%Med(contamedia)%Is%Dielectric = .FALSE.
         END IF
         !!!!
         ALLOCATE (sgg%Med(contamedia)%EDispersive(1)%C11(1:sgg%Med(contamedia)%EDispersive(1)%NumPolRes11), &
         &          sgg%Med(contamedia)%EDispersive(1)%a11(1:sgg%Med(contamedia)%EDispersive(1)%NumPolRes11))
         ALLOCATE (sgg%Med(contamedia)%EDispersive(1)%C12(1:sgg%Med(contamedia)%EDispersive(1)%NumPolRes12), &
         &          sgg%Med(contamedia)%EDispersive(1)%a12(1:sgg%Med(contamedia)%EDispersive(1)%NumPolRes12))
         ALLOCATE (sgg%Med(contamedia)%EDispersive(1)%C13(1:sgg%Med(contamedia)%EDispersive(1)%NumPolRes13), &
         &          sgg%Med(contamedia)%EDispersive(1)%a13(1:sgg%Med(contamedia)%EDispersive(1)%NumPolRes13))
         ALLOCATE (sgg%Med(contamedia)%EDispersive(1)%C22(1:sgg%Med(contamedia)%EDispersive(1)%NumPolRes22), &
         &          sgg%Med(contamedia)%EDispersive(1)%a22(1:sgg%Med(contamedia)%EDispersive(1)%NumPolRes22))
         ALLOCATE (sgg%Med(contamedia)%EDispersive(1)%C23(1:sgg%Med(contamedia)%EDispersive(1)%NumPolRes23), &
         &          sgg%Med(contamedia)%EDispersive(1)%a23(1:sgg%Med(contamedia)%EDispersive(1)%NumPolRes23))
         ALLOCATE (sgg%Med(contamedia)%EDispersive(1)%C33(1:sgg%Med(contamedia)%EDispersive(1)%NumPolRes33), &
         &          sgg%Med(contamedia)%EDispersive(1)%a33(1:sgg%Med(contamedia)%EDispersive(1)%NumPolRes33))
         !
         ALLOCATE (sgg%Med(contamedia)%MDispersive(1)%C11(1:sgg%Med(contamedia)%MDispersive(1)%NumPolRes11), &
         &          sgg%Med(contamedia)%MDispersive(1)%a11(1:sgg%Med(contamedia)%MDispersive(1)%NumPolRes11))
         ALLOCATE (sgg%Med(contamedia)%MDispersive(1)%C12(1:sgg%Med(contamedia)%MDispersive(1)%NumPolRes12), &
         &          sgg%Med(contamedia)%MDispersive(1)%a12(1:sgg%Med(contamedia)%MDispersive(1)%NumPolRes12))
         ALLOCATE (sgg%Med(contamedia)%MDispersive(1)%C13(1:sgg%Med(contamedia)%MDispersive(1)%NumPolRes13), &
         &          sgg%Med(contamedia)%MDispersive(1)%a13(1:sgg%Med(contamedia)%MDispersive(1)%NumPolRes13))
         ALLOCATE (sgg%Med(contamedia)%MDispersive(1)%C22(1:sgg%Med(contamedia)%MDispersive(1)%NumPolRes22), &
         &          sgg%Med(contamedia)%MDispersive(1)%a22(1:sgg%Med(contamedia)%MDispersive(1)%NumPolRes22))
         ALLOCATE (sgg%Med(contamedia)%MDispersive(1)%C23(1:sgg%Med(contamedia)%MDispersive(1)%NumPolRes23), &
         &          sgg%Med(contamedia)%MDispersive(1)%a23(1:sgg%Med(contamedia)%MDispersive(1)%NumPolRes23))
         ALLOCATE (sgg%Med(contamedia)%MDispersive(1)%C33(1:sgg%Med(contamedia)%MDispersive(1)%NumPolRes33), &
         &          sgg%Med(contamedia)%MDispersive(1)%a33(1:sgg%Med(contamedia)%MDispersive(1)%NumPolRes33))
         !
         DO k1 = 1,    (fdgeom%k11)
            sgg%Med(contamedia)%EDispersive(1)%C11(k1) =  (fdgeom%a11(k1))
            sgg%Med(contamedia)%EDispersive(1)%a11(k1) = (-fdgeom%b11(k1)) !el polo de ORIGINAL esta cambiado de signo
         END DO
         DO k1 = 1,    (fdgeom%k12)
            sgg%Med(contamedia)%EDispersive(1)%C12(k1) =  (fdgeom%a12(k1))
            sgg%Med(contamedia)%EDispersive(1)%a12(k1) = (-fdgeom%b12(k1))
         END DO
         DO k1 = 1,    (fdgeom%k13)
            sgg%Med(contamedia)%EDispersive(1)%C13(k1) =  (fdgeom%a13(k1))
            sgg%Med(contamedia)%EDispersive(1)%a13(k1) = (-fdgeom%b13(k1))
         END DO
         DO k1 = 1,    (fdgeom%k22)
            sgg%Med(contamedia)%EDispersive(1)%C22(k1) =  (fdgeom%a22(k1))
            sgg%Med(contamedia)%EDispersive(1)%a22(k1) = (-fdgeom%b22(k1))
         END DO
         DO k1 = 1,    (fdgeom%k23)
            sgg%Med(contamedia)%EDispersive(1)%C23(k1) =  (fdgeom%a23(k1))
            sgg%Med(contamedia)%EDispersive(1)%a23(k1) = (-fdgeom%b23(k1))
         END DO
         DO k1 = 1,    (fdgeom%k33)
            sgg%Med(contamedia)%EDispersive(1)%C33(k1) =  (fdgeom%a33(k1))
            sgg%Med(contamedia)%EDispersive(1)%a33(k1) = (-fdgeom%b33(k1))
         END DO
         !
         DO k1 = 1,   (fdgeom%KM11)
            sgg%Med(contamedia)%MDispersive(1)%C11(k1) = (fdgeom%aM11(k1))
            sgg%Med(contamedia)%MDispersive(1)%a11(k1) = (-fdgeom%bM11(k1))
         END DO
         DO k1 = 1,   (fdgeom%KM12)
            sgg%Med(contamedia)%MDispersive(1)%C12(k1) =  (fdgeom%aM12(k1))
            sgg%Med(contamedia)%MDispersive(1)%a12(k1) = (-fdgeom%bM12(k1))
         END DO
         DO k1 = 1,   (fdgeom%KM13)
            sgg%Med(contamedia)%MDispersive(1)%C13(k1) =  (fdgeom%aM13(k1))
            sgg%Med(contamedia)%MDispersive(1)%a13(k1) = (-fdgeom%bM13(k1))
         END DO
         DO k1 = 1,   (fdgeom%KM22)
            sgg%Med(contamedia)%MDispersive(1)%C22(k1) =  (fdgeom%aM22(k1))
            sgg%Med(contamedia)%MDispersive(1)%a22(k1) = (-fdgeom%bM22(k1))
         END DO
         DO k1 = 1,   (fdgeom%KM23)
            sgg%Med(contamedia)%MDispersive(1)%C23(k1) =  (fdgeom%aM23(k1))
            sgg%Med(contamedia)%MDispersive(1)%a23(k1) = (-fdgeom%bM23(k1))
         END DO
         DO k1 = 1,   (fdgeom%KM33)
            sgg%Med(contamedia)%MDispersive(1)%C33(k1) =  (fdgeom%aM33(k1))
            sgg%Med(contamedia)%MDispersive(1)%a33(k1) = (-fdgeom%bM33(k1))
         END DO

         return

      end subroutine asignadisper
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   END SUBROUTINE read_geomData

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE read_limits_nogeom (layoutnumber,size, sgg, fullsize, SINPML_fullsize, this,MurAfterPML,mur_exist)
      TYPE (limit_t), DIMENSION (1:6) :: fullsize, SINPML_fullsize
      type (SGGFDTDINFO), intent(INOUT)    :: sgg

      TYPE (Parseador), INTENT (IN) :: this
      INTEGER (KIND=4) :: tama, i, field,j,k
      character(len=BUFSIZE) :: buff
      logical MurAfterPML,mur_exist
      integer (kind=4), intent(in)           ::  layoutnumber,size
      REAL (KIND=RKIND), POINTER, DIMENSION (:) ::    DummyD
      REAL (KIND=RKIND) :: delta
      !
      REAL (KIND=RKIND), POINTER, DIMENSION (:) :: lineasX, lineasY, lineasZ
      !
      ALLOCATE (sgg%dx(this%despl%mx1:this%despl%mx2-1), sgg%dy(this%despl%my1:this%despl%my2-1), &
      & sgg%dz(this%despl%mz1:this%despl%mz2-1))
      !displacement
      !materialMatrix
      tama = (this%despl%nx)
      IF (tama == 1) THEN
         DO i = this%despl%mx1, this%despl%mx2 - 1
            sgg%dx (i) = this%despl%desx(1)
         END DO
      ELSE
         IF (tama /= this%despl%mx2-this%despl%mx1) THEN
            buff='Tamanio discretizacion distinto de la region'
            CALL STOPONERROR(layoutnumber,size,buff)
         END IF
         DO i = this%despl%mx1, this%despl%mx2 - 1
            sgg%dx (i) = this%despl%desx(i)
         END DO
      END IF
      !
      !
      tama = (this%despl%nY)
      IF (tama == 1) THEN
         DO i = this%despl%my1, this%despl%my2 - 1
            sgg%dy (i) = this%despl%desY(1)
         END DO
      ELSE
         IF (tama /= this%despl%my2-this%despl%my1) THEN
            buff='Tamanio discretizacion distinto de la region'
            CALL STOPONERROR(layoutnumber,size,buff)
         END IF
         DO i = this%despl%my1, this%despl%my2 - 1
            sgg%dy (i) = this%despl%desY(i)
         END DO
      END IF
      !
      !
      tama = (this%despl%nZ)
      IF (tama == 1) THEN
         DO i = this%despl%mz1, this%despl%mz2 - 1
            sgg%dz (i) = this%despl%desZ(1)
         END DO
      ELSE
         IF (tama /= this%despl%mz2-this%despl%mz1) THEN
            buff='Tamanio discretizacion distinto de la region'
            CALL STOPONERROR(layoutnumber,size,buff)
         END IF
         DO i = this%despl%mz1, this%despl%mz2 - 1
            sgg%dz (i) = this%despl%desZ(i)
         END DO
      END IF
      !DISCRETIZATION LINES
      !displacement
      !materialMatrix
      tama = (this%despl%nx)
      ALLOCATE (lineasX(this%despl%mx1:this%despl%mx2))
      lineasX (this%despl%mx1) = this%despl%originX + this%despl%mx1 * this%despl%desx(1)
      IF (tama == 1) THEN
         DO i = this%despl%mx1, this%despl%mx2 - 1
            lineasX (i+1) = this%despl%desx(1) * (i-this%despl%mx1+1) + lineasX (this%despl%mx1)
         END DO
      ELSE
         IF (tama /= this%despl%mx2-this%despl%mx1) THEN
            buff='Tamanio discretizacion distinto de la region'
            CALL STOPONERROR(layoutnumber,size,buff)
         END IF
         DO i = this%despl%mx1, this%despl%mx2 - 1
            lineasX (i+1) = this%despl%desx(i) + lineasX (i)
         END DO
      END IF
      !Y
      tama = (this%despl%nY)
      ALLOCATE (lineasY(this%despl%my1:this%despl%my2))
      lineasY (this%despl%my1) = this%despl%originy+this%despl%my1*this%despl%desY(1)
      IF (tama == 1) THEN
         DO i = this%despl%my1, this%despl%my2 - 1
            lineasY (i+1) = this%despl%desY(1) * (i-this%despl%my1+1) + lineasY (this%despl%my1)
         END DO
      ELSE
         IF (tama /= this%despl%my2-this%despl%my1) then
            buff='Tamanio discretizacion distinto de la region'
            CALL STOPONERROR(layoutnumber,size,buff)
         endif
         DO i = this%despl%my1, this%despl%my2 - 1
            lineasY (i+1) = this%despl%desY(i) + lineasY (i)
         END DO
      END IF
      !Z
      tama = (this%despl%nZ)
      ALLOCATE (lineasZ(this%despl%mz1:this%despl%mz2))
      lineasZ (this%despl%mz1) = this%despl%originZ+this%despl%mz1*this%despl%desZ(1)
      IF (tama == 1) THEN
         DO i = this%despl%mz1, this%despl%mz2 - 1
            lineasZ (i+1) = this%despl%desZ(1) * (i-this%despl%mz1+1) + lineasZ (this%despl%mz1)
         END DO
      ELSE
         IF (tama /= this%despl%mz2-this%despl%mz1) then
            buff='Tamanio discretizacion distinto de la region'
            CALL STOPONERROR(layoutnumber,size,buff)
         endif
         DO i = this%despl%mz1, this%despl%mz2 - 1
            lineasZ (i+1) = this%despl%desZ(i) + lineasZ (i)
         END DO
      END IF
      !
      ALLOCATE (sgg%LineX(this%despl%mx1:this%despl%mx2), sgg%LineY(this%despl%my1:this%despl%my2), &
      & sgg%LineZ(this%despl%mz1:this%despl%mz2))
      !
      sgg%LineX (this%despl%mx1:this%despl%mx2) = lineasX
      sgg%LineY (this%despl%my1:this%despl%my2) = lineasY
      sgg%LineZ (this%despl%mz1:this%despl%mz2) = lineasZ
      DEALLOCATE (lineasX, lineasY, lineasZ)
      !General Parameter
      sgg%InitialTimeStep = 0
      sgg%TimeSteps = this%general%nmax
      sgg%dt = this%general%dt
      if (sgg%dt < abs(tiny(1.0_RKIND)*2.0_RKIND)) sgg%dt=huge(1.0_RKIND)/2.0_RKIND
      !border
      !this%BORDER%PROPIEDADESPML(I)%ORDEN no lo considero porque en el interior de mi programa lo pongo (esta a 2 normalmente)
      sgg%Border%IsBackPEC = .FALSE.
      sgg%Border%IsFrontPEC = .FALSE.
      sgg%Border%IsLeftPEC = .FALSE.
      sgg%Border%IsRightPEC = .FALSE.
      sgg%Border%IsUpPEC = .FALSE.
      sgg%Border%IsDownPEC = .FALSE.
      sgg%Border%IsBackPMC = .FALSE.
      sgg%Border%IsFrontPMC = .FALSE.
      sgg%Border%IsLeftPMC = .FALSE.
      sgg%Border%IsRightPMC = .FALSE.
      sgg%Border%IsUpPMC = .FALSE.
      sgg%Border%IsDownPMC = .FALSE.
      sgg%Border%IsBackPML = .FALSE.
      sgg%Border%IsFrontPML = .FALSE.
      sgg%Border%IsLeftPML = .FALSE.
      sgg%Border%IsRightPML = .FALSE.
      sgg%Border%IsUpPML = .FALSE.
      sgg%Border%IsDownPML = .FALSE.
      sgg%Border%IsBackPeriodic = .FALSE.
      sgg%Border%IsFrontPeriodic = .FALSE.
      sgg%Border%IsLeftPeriodic = .FALSE.
      sgg%Border%IsRightPeriodic = .FALSE.
      sgg%Border%IsUpPeriodic = .FALSE.
      sgg%Border%IsDownPeriodic = .FALSE.
      sgg%Border%IsBackMUR = .FALSE.
      sgg%Border%IsFrontMUR = .FALSE.
      sgg%Border%IsLeftMUR = .FALSE.
      sgg%Border%IsRightMUR = .FALSE.
      sgg%Border%IsUpMUR = .FALSE.
      sgg%Border%IsDownMUR = .FALSE.
      sgg%PML%NumLayers = 0
      DO i = 1, 6
         IF (this%front%tipofrontera(i) == F_PML) THEN
            SELECT CASE (i)
               !xmin
             CASE (1)
               sgg%Border%IsBackPML = .TRUE.
               sgg%PML%NumLayers (icoord, comi) = this%front%PROPIEDADESPML(i)%NUMCAPAS
               sgg%PML%CoeffReflPML (icoord, comi) = this%front%PROPIEDADESPML(i)%REFL
               if (sgg%PML%CoeffReflPML (icoord, comi)>=1.0_RKIND) sgg%PML%CoeffReflPML(icoord, comi)=0.99999d0
               sgg%PML%orden (icoord, comi) = this%front%PROPIEDADESPML(i)%orden
               !xmax
             CASE (2)
               sgg%Border%IsFrontPML = .TRUE.
               sgg%PML%NumLayers (icoord, fine) = this%front%PROPIEDADESPML(i)%NUMCAPAS
               sgg%PML%CoeffReflPML (icoord, fine) = this%front%PROPIEDADESPML(i)%REFL
               if (sgg%PML%CoeffReflPML (icoord, fine)>=1.0_RKIND) sgg%PML%CoeffReflPML(icoord, fine)=0.99999d0
               sgg%PML%orden (icoord, fine) = this%front%PROPIEDADESPML(i)%orden
               !ymin
             CASE (3)
               sgg%Border%IsLeftPML = .TRUE.
               sgg%PML%NumLayers (jcoord, comi) = this%front%PROPIEDADESPML(i)%NUMCAPAS
               sgg%PML%CoeffReflPML (jcoord, comi) = this%front%PROPIEDADESPML(i)%REFL
               if (sgg%PML%CoeffReflPML (jcoord, comi)>=1.0_RKIND) sgg%PML%CoeffReflPML(jcoord, comi)=0.99999d0
               sgg%PML%orden (jcoord, comi) = this%front%PROPIEDADESPML(i)%orden
               !ymax
             CASE (4)
               sgg%Border%IsRightPML = .TRUE.
               sgg%PML%NumLayers (jcoord, fine) = this%front%PROPIEDADESPML(i)%NUMCAPAS
               sgg%PML%CoeffReflPML (jcoord, fine) = this%front%PROPIEDADESPML(i)%REFL
               if (sgg%PML%CoeffReflPML (jcoord, fine)>=1.0_RKIND) sgg%PML%CoeffReflPML(jcoord, fine)=0.99999d0
               sgg%PML%orden (jcoord, fine) = this%front%PROPIEDADESPML(i)%orden
               !zmin
             CASE (5)
               sgg%Border%IsDownPML = .TRUE.
               sgg%PML%NumLayers (kcoord, comi) = this%front%PROPIEDADESPML(i)%NUMCAPAS
               sgg%PML%CoeffReflPML (kcoord, comi) = this%front%PROPIEDADESPML(i)%REFL
               if (sgg%PML%CoeffReflPML (kcoord, comi)>=1.0_RKIND) sgg%PML%CoeffReflPML(kcoord, comi)=0.99999d0
               sgg%PML%orden (kcoord, comi) = this%front%PROPIEDADESPML(i)%orden
               !zmax
             CASE (6)
               sgg%Border%IsUpPML = .TRUE.
               sgg%PML%NumLayers (kcoord, fine) = this%front%PROPIEDADESPML(i)%NUMCAPAS
               sgg%PML%CoeffReflPML (kcoord, fine) = this%front%PROPIEDADESPML(i)%REFL
               if (sgg%PML%CoeffReflPML (kcoord, fine)>=1.0_RKIND) sgg%PML%CoeffReflPML(kcoord, fine)=0.99999d0
               sgg%PML%orden (kcoord, fine) = this%front%PROPIEDADESPML(i)%orden
            END SELECT
         elseIF (this%front%tipofrontera(i) == F_MUR) THEN
            mur_exist=.true.
            SELECT CASE (i)
               !xmin
             CASE (1)
               sgg%Border%IsBackMUR = .TRUE.
               !xmax
             CASE (2)
               sgg%Border%IsFrontMUR = .TRUE.
               !ymin
             CASE (3)
               sgg%Border%IsLeftMUR = .TRUE.
               !ymax
             CASE (4)
               sgg%Border%IsRightMUR = .TRUE.
               !zmin
             CASE (5)
               sgg%Border%IsDownMUR = .TRUE.
               !zmax
             CASE (6)
               sgg%Border%IsUpMUR = .TRUE.
            END SELECT
         ELSE IF (this%front%tipofrontera(i) == F_PEC) THEN
            SELECT CASE (i)
               !xmin
             CASE (1)
               sgg%Border%IsBackPEC = .TRUE.
               !xmax
             CASE (2)
               sgg%Border%IsFrontPEC = .TRUE.
               !ymin
             CASE (3)
               sgg%Border%IsLeftPEC = .TRUE.
               !ymax
             CASE (4)
               sgg%Border%IsRightPEC = .TRUE.
               !zmin
             CASE (5)
               sgg%Border%IsDownPEC = .TRUE.
               !zmax
             CASE (6)
               sgg%Border%IsUpPEC = .TRUE.
            END SELECT
         ELSE IF (this%front%tipofrontera(i) == F_PMC) THEN
            SELECT CASE (i)
               !xmin
             CASE (1)
               sgg%Border%IsBackPMC = .TRUE.
               !xmax
             CASE (2)
               sgg%Border%IsFrontPMC = .TRUE.
               !ymin
             CASE (3)
               sgg%Border%IsLeftPMC = .TRUE.
               !ymax
             CASE (4)
               sgg%Border%IsRightPMC = .TRUE.
               !zmin
             CASE (5)
               sgg%Border%IsDownPMC = .TRUE.
               !zmax
             CASE (6)
               sgg%Border%IsUpPMC = .TRUE.
            END SELECT
         ELSE IF (this%front%tipofrontera(i) == F_Per) THEN
            SELECT CASE (i)
               !xmin
             CASE (1)
               sgg%Border%IsBackPeriodic = .TRUE.
               !xmax
             CASE (2)
               sgg%Border%IsFrontPeriodic = .TRUE.
               !ymin
             CASE (3)
               sgg%Border%IsLeftPeriodic = .TRUE.
               !ymax
             CASE (4)
               sgg%Border%IsRightPeriodic = .TRUE.
               !zmin
             CASE (5)
               sgg%Border%IsDownPeriodic = .TRUE.
               !zmax
             CASE (6)
               sgg%Border%IsUpPeriodic = .TRUE.
            END SELECT
         END IF
      END DO
      !assign limits
      DO field = iEx, iHz
         SINPML_fullsize(field)%XI = this%despl%mx1
         SINPML_fullsize(field)%YI = this%despl%my1
         SINPML_fullsize(field)%ZI = this%despl%mz1
         SINPML_fullsize(field)%XE = this%despl%mx2
         SINPML_fullsize(field)%YE = this%despl%my2
         SINPML_fullsize(field)%ZE = this%despl%mz2
      END DO
      !adjust the endings
      SINPML_fullsize(iEx)%XE = SINPML_fullsize(iEx)%XE - 1
      SINPML_fullsize(iEy)%YE = SINPML_fullsize(iEy)%YE - 1
      SINPML_fullsize(iEz)%ZE = SINPML_fullsize(iEz)%ZE - 1
      !
      !
      SINPML_fullsize(iHx)%YE = SINPML_fullsize(iHx)%YE - 1
      SINPML_fullsize(iHx)%ZE = SINPML_fullsize(iHx)%ZE - 1
      SINPML_fullsize(iHy)%ZE = SINPML_fullsize(iHy)%ZE - 1
      SINPML_fullsize(iHy)%XE = SINPML_fullsize(iHy)%XE - 1
      SINPML_fullsize(iHz)%XE = SINPML_fullsize(iHz)%XE - 1
      SINPML_fullsize(iHz)%YE = SINPML_fullsize(iHz)%YE - 1
      !
      DO field = iEx, iHz
         fullsize(field)%XI = SINPML_fullsize(field)%XI - sgg%PML%NumLayers(icoord, comi)
         fullsize(field)%YI = SINPML_fullsize(field)%YI - sgg%PML%NumLayers(jcoord, comi)
         fullsize(field)%ZI = SINPML_fullsize(field)%ZI - sgg%PML%NumLayers(kcoord, comi)
         fullsize(field)%XE = SINPML_fullsize(field)%XE + sgg%PML%NumLayers(icoord, fine)
         fullsize(field)%YE = SINPML_fullsize(field)%YE + sgg%PML%NumLayers(jcoord, fine)
         fullsize(field)%ZE = SINPML_fullsize(field)%ZE + sgg%PML%NumLayers(kcoord, fine)
      END DO
      !
      !readjust mur boundaries if necessary

      sgg%Border%IsBackMUR  = (sgg%Border%IsBackMUR  ).or.(sgg%Border%IsBackPML   .and. MurAfterPML)
      sgg%Border%IsFrontMUR = (sgg%Border%IsFrontMUR ).or.(sgg%Border%IsFrontPML  .and. MurAfterPML)
      sgg%Border%IsLeftMUR  = (sgg%Border%IsLeftMUR  ).or.(sgg%Border%IsLeftPML   .and. MurAfterPML)
      sgg%Border%IsRightMUR = (sgg%Border%IsRightMUR ).or.(sgg%Border%IsRightPML  .and. MurAfterPML)
      sgg%Border%IsUpMUR    = (sgg%Border%IsUpMUR    ).or.(sgg%Border%IsUpPML     .and. MurAfterPML)
      sgg%Border%IsDownMUR  = (sgg%Border%IsDownMUR  ).or.(sgg%Border%IsDownPML   .and. MurAfterPML)


      !readjust space steps accordingly 140815 para que esten bien allocateados los dx
      ! Discretization Lines Matrix Resizing to accomodate PML regions
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ALLOCATE (DummyD(SINPML_fullsize(iHx)%XI:SINPML_fullsize(iHx)%XE-1))
      DummyD = sgg%dx
      DEALLOCATE (sgg%dx)
      !

      sgg%allocDxI=-sgg%PML%NumLayers(1, 1)+SINPML_fullsize(iHx)%XI-1-1
      sgg%allocDxE= SINPML_fullsize(iHx)%XE+sgg%PML%NumLayers(1, 2)+1+1
      ALLOCATE (sgg%dx(sgg%allocDxI:sgg%allocDxE))
      !
      sgg%dx (SINPML_fullsize(iHx)%XI:SINPML_fullsize(iHx)%XE-1) = DummyD (SINPML_fullsize(iHx)%XI:SINPML_fullsize(iHx)%XE-1)
      DEALLOCATE (DummyD)
      !
      ALLOCATE (DummyD(SINPML_fullsize(iHy)%YI:SINPML_fullsize(iHy)%YE-1))
      DummyD = sgg%dy
      DEALLOCATE (sgg%dy)
      !
      sgg%allocDyI=-sgg%PML%NumLayers(2, 1)+SINPML_fullsize(iHy)%YI-1-1
      sgg%allocDyE= SINPML_fullsize(iHy)%YE+sgg%PML%NumLayers(2, 2)+1+1
      ALLOCATE (sgg%dy(sgg%allocDyI:sgg%allocDyE))
      !
      sgg%dy (SINPML_fullsize(iHy)%YI:SINPML_fullsize(iHy)%YE-1) = DummyD (SINPML_fullsize(iHy)%YI:SINPML_fullsize(iHy)%YE-1)
      DEALLOCATE (DummyD)
      !
      ALLOCATE (DummyD(SINPML_fullsize(iHz)%ZI:SINPML_fullsize(iHz)%ZE-1))
      DummyD = sgg%dz
      DEALLOCATE (sgg%dz)
      !
      sgg%allocDzI=-sgg%PML%NumLayers(3, 1)+SINPML_fullsize(iHz)%ZI-1-1
      sgg%allocDzE= SINPML_fullsize(iHz)%ZE+sgg%PML%NumLayers(3, 2)+1+1
      ALLOCATE (sgg%dz(sgg%allocDzI:sgg%allocDzE))
      sgg%dz (SINPML_fullsize(iHz)%ZI:SINPML_fullsize(iHz)%ZE-1) = DummyD (SINPML_fullsize(iHz)%ZI:SINPML_fullsize(iHz)%ZE-1)
      DEALLOCATE (DummyD)
      !
      delta = sgg%dx (SINPML_fullsize(iHx)%XI)
      DO i = SINPML_fullsize(iHx)%XI - 1, SINPML_fullsize(iHx)%XI - 1 - sgg%PML%NumLayers(1, 1) - 1, - 1
         sgg%dx (i) = delta
      END DO
      delta = sgg%dx (SINPML_fullsize(iHx)%XE-1)
      DO i = SINPML_fullsize(iHx)%XE, SINPML_fullsize(iHx)%XE + sgg%PML%NumLayers(1, 2) + 1 + 1
         sgg%dx (i) = delta
      END DO
      !
      delta = sgg%dy (SINPML_fullsize(iHy)%YI)
      DO j = SINPML_fullsize(iHy)%YI - 1, SINPML_fullsize(iHy)%YI - 1 - sgg%PML%NumLayers(2, 1) - 1, - 1
         sgg%dy (j) = delta
      END DO
      !
      delta = sgg%dy (SINPML_fullsize(iHy)%YE-1)
      DO j = SINPML_fullsize(iHy)%YE, SINPML_fullsize(iHy)%YE + sgg%PML%NumLayers(2, 2) + 1 + 1
         sgg%dy (j) = delta
      END DO
      !
      delta = sgg%dz (SINPML_fullsize(iHz)%ZI)
      DO k = SINPML_fullsize(iHz)%ZI - 1, SINPML_fullsize(iHz)%ZI - 1 - sgg%PML%NumLayers(3, 1) - 1, - 1
         sgg%dz (k) = delta
      END DO
      !
      delta = sgg%dz (SINPML_fullsize(iHz)%ZE-1)
      DO k = SINPML_fullsize(iHz)%ZE, SINPML_fullsize(iHz)%ZE + sgg%PML%NumLayers(3, 2) + 1 + 1
         sgg%dz (k) = delta
      END DO
      !DISCRETIZATION LINES (TO BE DEPRECATED IN A NEAR FUTURE, ONLY NEEDED BY THE PLANEWAVE CORNER ROUTINE)
      !
      ALLOCATE (DummyD(SINPML_fullsize(iHx)%XI:SINPML_fullsize(iHx)%XE))
      DummyD = sgg%LineX
      DEALLOCATE (sgg%LineX)
      ALLOCATE (sgg%LineX(-sgg%PML%NumLayers(1, 1)+SINPML_fullsize(iHx)%XI-1:SINPML_fullsize(iHx)%XE+sgg%PML%NumLayers(1, 2)+1))
      sgg%LineX (SINPML_fullsize(iHx)%XI:SINPML_fullsize(iHx)%XE) = DummyD (SINPML_fullsize(iHx)%XI:SINPML_fullsize(iHx)%XE)
      DEALLOCATE (DummyD)
      !
      ALLOCATE (DummyD(SINPML_fullsize(iHy)%YI:SINPML_fullsize(iHy)%YE))
      DummyD = sgg%LineY
      DEALLOCATE (sgg%LineY)
      ALLOCATE (sgg%LineY(-sgg%PML%NumLayers(2, 1)+SINPML_fullsize(iHy)%YI-1:SINPML_fullsize(iHy)%YE+sgg%PML%NumLayers(2, 2)+1))
      sgg%LineY (SINPML_fullsize(iHy)%YI:SINPML_fullsize(iHy)%YE) = DummyD (SINPML_fullsize(iHy)%YI:SINPML_fullsize(iHy)%YE)
      DEALLOCATE (DummyD)
      !
      ALLOCATE (DummyD(SINPML_fullsize(iHz)%ZI:SINPML_fullsize(iHz)%ZE))
      DummyD = sgg%LineZ
      DEALLOCATE (sgg%LineZ)
      ALLOCATE (sgg%LineZ(-sgg%PML%NumLayers(3, 1)+SINPML_fullsize(iHz)%ZI-1:SINPML_fullsize(iHz)%ZE+sgg%PML%NumLayers(3, 2)+1))
      sgg%LineZ (SINPML_fullsize(iHz)%ZI:SINPML_fullsize(iHz)%ZE) = DummyD (SINPML_fullsize(iHz)%ZI:SINPML_fullsize(iHz)%ZE)
      DEALLOCATE (DummyD)
      !
      delta = sgg%LineX (SINPML_fullsize(iHx)%XI+1) - sgg%LineX(SINPML_fullsize(iHx)%XI)
      DO i = SINPML_fullsize(iHx)%XI - 1, SINPML_fullsize(iHx)%XI - 1 - sgg%PML%NumLayers(1, 1), - 1
         sgg%LineX (i) = sgg%LineX(i+1) - delta
      END DO
      delta = sgg%LineX (SINPML_fullsize(iHx)%XE) - sgg%LineX(SINPML_fullsize(iHx)%XE-1)
      DO i = SINPML_fullsize(iHx)%XE + 1, SINPML_fullsize(iHx)%XE + sgg%PML%NumLayers(1, 2) + 1
         sgg%LineX (i) = sgg%LineX(i-1) + delta
      END DO
      !
      delta = sgg%LineY (SINPML_fullsize(iHy)%YI+1) - sgg%LineY(SINPML_fullsize(iHy)%YI)
      DO j = SINPML_fullsize(iHy)%YI - 1, SINPML_fullsize(iHy)%YI - 1 - sgg%PML%NumLayers(2, 1), - 1
         sgg%LineY (j) = sgg%LineY(j+1) - delta
      END DO
      !
      delta = sgg%LineY (SINPML_fullsize(iHy)%YE) - sgg%LineY(SINPML_fullsize(iHy)%YE-1)
      DO j = SINPML_fullsize(iHy)%YE + 1, SINPML_fullsize(iHy)%YE + sgg%PML%NumLayers(2, 2) + 1
         sgg%LineY (j) = sgg%LineY(j-1) + delta
      END DO
      !
      delta = sgg%LineZ (SINPML_fullsize(iHz)%ZI+1) - sgg%LineZ(SINPML_fullsize(iHz)%ZI)
      DO k = SINPML_fullsize(iHz)%ZI - 1, SINPML_fullsize(iHz)%ZI - 1 - sgg%PML%NumLayers(3, 1), - 1
         sgg%LineZ (k) = sgg%LineZ(k+1) - delta
      END DO
      !
      delta = sgg%LineZ (SINPML_fullsize(iHz)%ZE) - sgg%LineZ(SINPML_fullsize(iHz)%ZE-1)
      DO k = SINPML_fullsize(iHz)%ZE + 1, SINPML_fullsize(iHz)%ZE + sgg%PML%NumLayers(3, 2) + 1
         sgg%LineZ (k) = sgg%LineZ(k-1) + delta
      END DO
      !2012
      !update actual number of media
      !!!!!!!

      !!      DO i = SINPML_fullsize(iHx)%XI - 1,  SINPML_fullsize(iHx)%XE + sgg%PML%NumLayers(1, 2) + 1
      !!        write (6678,*) 'lineX ',i, sgg%LineX (i)
      !!      END DO
      !!!
      !!      DO j = SINPML_fullsize(iHy)%YI - 1,  SINPML_fullsize(iHy)%YE + sgg%PML%NumLayers(2, 2) + 1
      !!        write (6678,*) 'lineY ',j,sgg%LineY (j)
      !!      END DO
      !!!
      !!      DO k = SINPML_fullsize(iHz)%ZI - 1, SINPML_fullsize(iHz)%ZE + sgg%PML%NumLayers(3, 2) + 1
      !!        write (6678,*) 'lineZ ',k,sgg%LineZ (k)
      !!      END DO


      RETURN
      !
   END SUBROUTINE
   !



   !!!!!!!!!!!!!!!!PREPROCESADOR PARA SKIN-DEPTH 09/07/13
   subroutine prepro_skindepth(this,fichin)
      integer pozi,tama,j,k
      character (LEN=BUFSIZE)       ::   multiportFile
      TYPE (Parseador), INTENT (IN) :: this
      CHARACTER (LEN=*), INTENT (IN) :: fichin
      character (LEN=BUFSIZE)  ::  restocadena
      integer :: my_iostat

      open (unit=7533,file='UGRskindepthmatlab.layers')
      close(7533,status='delete')
      my_iostat=0
9306  if(my_iostat /= 0) write(*,fmt='(a)',advance='no'), '.' !!if(my_iostat /= 0) print '(i5,a1,i4,2x,a)',9306,'.',quienmpi,'UGRskindepthmatlab.layers'
      open (unit=7533,file='UGRskindepthmatlab.layers',err=9306,iostat=my_iostat,status='new',action='write')

      tama = this%LossyThinSurfs%length
      DO j = 1, tama
         IF (abs(this%LossyThinSurfs%cs(j)%SigmaM(1)) <= 1.0e-2_RKIND ) THEN !SGBCs que hay que sustituir
            multiportFile =  trim(adjustl(this%LossyThinSurfs%cs(j)%files)) // '_z11.txt'
            !
            !09/07/13 !los SGBCs con skindepth se deben preprocesar

            !crea el fichero de entrada para usar con el compilado de Matlab
            pozi=index(multiportFile,'_z11.txt')
            write(7533,'(a)') trim(adjustl(multiportFile(1:pozi-1)))
            write(7533,*)     'layers    ',this%LossyThinSurfs%cs(j)%numcapas
            do k=1,this%LossyThinSurfs%cs(j)%numcapas
               write(7533,*) 'eps       ',k,this%LossyThinSurfs%cs(j)%eps(k)
               write(7533,*) 'mu        ',k,this%LossyThinSurfs%cs(j)%mu(k)
               write(7533,*) 'sigma     ',k,this%LossyThinSurfs%cs(j)%sigma(k)
               write(7533,*) 'thickness ',k,this%LossyThinSurfs%cs(j)%thk(k)
            end do
            write(7533,*) 'fmin      ',10**4
            write(7533,*) 'fmax      ',10**9
            write(7533,*) 'order     ',24
         ENDIF
      END DO
      close(7533)
      RETURN

   end subroutine !prepro_skindepth
   !!!!!!!!end 09/07/13

#ifdef CompileWithConformal
   subroutine AssigLossyOrPECtoNodes(sgg,media,conf_conflicts,input_conformal_flag)
#else
   subroutine AssigLossyOrPECtoNodes(sgg,media)
#endif
      type (SGGFDTDINFO), intent(INOUT) :: sgg
      type (media_matrices_t), intent(inout) :: media

      logical :: ispec, isSGBC, IsComposite, islossy, input_conformal_flag,NODALMENTEIGUALES,iguaSGM,iguaSIG,iguaMUR,iguaPEC,iguaLOS,iguaEPR,ISconformal
      REAL (KIND=RKIND)   :: sigt,epst,SIGMA,SIGMAM,EPR,MUR
      integer (KIND=4) i,j,k,n,kmenos1,jmenos1,imenos1,med(0:5),r,imed,i1
      character(len=BUFSIZE) :: buff

#ifdef CompileWithConformal
      type (conf_conflicts_t), pointer  :: conf_conflicts
      type(conf_node_t), dimension (:), pointer :: conf_busy_node
      integer (KIND=confIKIND) :: dims
      logical :: mediois1,mediois2,mediois3
      integer (KIND=INTEGERSIZEOFMEDIAMATRICES) :: medio1,medio2,medio3
#endif


      !!!lo dejo aqui para en un futuro hacerlo crecer y asignar la informacion de nodo con su tipo de material correctamente
      !!copiado de la casuistica de wires !310715  !no esta anulado alli, pero cuando se haga crecer habra que aumenta el sgg%med para acomodar a los nuevos materiales nodales....

      do k= sgg%Alloc(iEz)%ZI , sgg%Alloc(iEz)%ZE
         do j= sgg%Alloc(iEy)%YI , sgg%Alloc(iEy)%YE
            do i= sgg%Alloc(iEx)%XI , sgg%Alloc(iEx)%XE
               imenos1= i-1
               jmenos1= j-1
               kmenos1= k-1
               if (i-1 <  sgg%alloc(iEx)%XI) imenos1=i
               if (j-1 <  sgg%alloc(iEy)%YI) jmenos1=j
               if (k-1 <  sgg%alloc(iEz)%ZI) kmenos1=k

               med(0)  = media%sggMiEx(i       , j       , k       )
               med(1)  = media%sggMiEx(imenos1 , j       , k       )
               med(2)  = media%sggMiEy(i       , j       , k       )
               med(3)  = media%sggMiEy(i       , jmenos1 , k       )
               med(4)  = media%sggMiEz(i       , j       , k       )
               med(5)  = media%sggMiEz(i       , j       , kmenos1 )
               sigma                       = 0.0_RKIND
               sigmam                      = 0.0_RKIND
               epr                         = 0.0_RKIND
               mur                         = 0.0_RKIND
               ISPEC                       =.FALSE.
               ISLOSSY                     =.FALSE.
               do i1=0,5
                  imed=med(i1)
                  sigma                       = max(sigma,                      sgg%Med(imed)%sigma)
                  sigmam                      = max(sigmam,                     sgg%Med(imed)%sigmam)
                  epr                         =     epr                       + sgg%Med(imed)%epr / 6.0_RKIND
                  mur                         =     mur                       + sgg%Med(imed)%mur /6.0_RKIND
                  if ((sgg%med(imed)%is%PEC).or.(imed==0)) isPEC = .true.
               end do
               if ( (.not.isPEC).and.(sigma >= 1e-4) ) then
                  islossy = .true.
               else
                  islossy = .false.
               endif
               !  CREAR NUEVO MEDIO y asignarle sus propiedades de acuerdo a sus adyacencias
               IF (.not.(sgg%med(med(0))%is%PML.OR.sgg%med(med(1))%is%PML.OR.sgg%med(med(2))%is%PML.OR.sgg%med(med(3))%is%PML.OR.sgg%med(med(4))%is%PML.OR.sgg%med(med(5))%is%PML)) THEN
                  IF ((MED(0)/=MED(1)).OR.(MED(1)/=MED(2)).OR.(MED(2)/=MED(3)).OR.(MED(3)/=MED(4)).OR.(MED(4)/=MED(5)).OR.(MED(5)/=MED(0))) THEN
                     NODALMENTEIGUALES=.FALSE.
                     busqueda: DO I1=0,SGG%NUMMEDIA
!cambios 230817 bug milano borja en rutina iguales
                        iguaSGM=IGUALES(SGG%MED(I1)%SIGMAM,SIGMAM)
                        iguaSIG=IGUALES(SGG%MED(I1)%SIGMA,SIGMA) !sgg230817 al poner sigma 1e29 y pec ademas, la rutina de iguales fallaba
                        iguaEPR=IGUALES(SGG%MED(I1)%EPR,EPR)
                        iguaMUR=IGUALES(SGG%MED(I1)%MUR,MUR)
                        iguaPEC=(SGG%MED(I1)%iS%PEC.eqv.ISPEC)
                        iguaLOS=(SGG%MED(I1)%iS%LOSSY.eqv.ISLOSSY)
                        ISconformal=((SGG%MED(I1)%iS%already_YEEadvanced_byconformal).or.(SGG%MED(I1)%is%split_and_useless))
                        NODALMENTEIGUALES=NODALMENTEIGUALES.OR.(((iguaSGM.and.iguaSIG.and.iguaEPR.and.iguaMUR.and.iguaLOS).OR.iguaPEC).and.(.not.ISconformal))
                        if (nodalmenteiguales) exit busqueda
                     END DO busqueda
                     IF (.NOT.NODALMENTEIGUALES) THEN
                        IF (SGG%NUMMEDIA+1 > SGG%ALLOCmed) then
                           call READJUST(SGG%ALLOCmed,sgg%med,2*SGG%ALLOCmed) !LO HAgo REallocatando al doble. gENERO NUEVO PARAMETRO sgg%ALLOCmed. Pero esto es un guirigay.... 261115
                        endif
                        SGG%NUMMEDIA=SGG%NUMMEDIA+1
                        media%sggMiNo(i,j,k)=SGG%NUMMEDIA
                        r=media%sggMiNo(i,j,k)
                        sgg%med(r)%sigma =sigma
                        sgg%med(r)%sigmam=sigmam
                        sgg%med(r)%epr   =epr
                        sgg%med(r)%mur   =mur
                        sgg%med(r)%is%PEC = ISPEC !ojo con estos medios que el sistema de prioridades ya no les afecta porque esta rutina va despues del preprocess 03116
                        sgg%med(r)%is%LOSSY = ISLOSSY
                        sgg%med(r)%is%needed = .true.  !sgg 220817 por defecto lo he puesto en readjust a false
!write(113,*) '.NOT.NODALMENTEIGUALES--> ',i,j,k,' - ',med(0),med(1),med(2),med(3),med(4),med(5),' - ',SGG%NUMMEDIA
                     ELSE
                        media%sggMiNo(i,j,k)=i1  !PUEDE QUE NO SEAN IGUALES PERO NODALMENTE LO SON (SOLO A EFECTOS DE SIGMA,EPR,SIGMAM,MUR,ISLOSSY,ISPEC
                        !bug 060417 debo ponerlo al medio que ha encontrado igual (i1) y estaba a med(0)!!!!
!write(114,*) '.YES.NODALMENTEIGUALES--> ',i,j,k,' - ',med(0),med(1),med(2),med(3),med(4),med(5),' - ',SGG%NUMMEDIA
                     ENDIF
                  else
                     media%sggMiNo(i,j,k)=MED(0)  !todos iguales
                  endif
               endif !del no es pml
               !!!!aqui habra luego que ir creando y almacenando lo nuevos tipos de medio nodales en funcion de los sigt y epst para que wires use directamente esa info
            end do
         end do
      end do

!CORRIGE AHORA CON LA INFO CONFORMAL (HABIA UN BUG 211116 PQ EL NODALMENTE IGUALES LO MACHACABA SI ESTO SE HACIA ANTES)
#ifdef CompileWithConformal
      if(input_conformal_flag)then
         dims      =conf_conflicts%conf_busy_nodesGroup%dims
         conf_busy_node=> conf_conflicts%conf_busy_nodesGroup%conf_busy_node
         do n=1,dims
            i=conf_busy_node(n)%i
            j=conf_busy_node(n)%j
            k=conf_busy_node(n)%k
            if (conf_busy_node(n)%ispec) then
               media%sggMiNo(i,j,k)=0
            endif
         end do
         !!!conf_busy_node(n)%i = 0
         !!!conf_busy_node(n)%j = 0
         !!!conf_busy_node(n)%k = 0
         !!!conf_busy_node(n)%isfree  = .true.
         !!!conf_busy_node(n)%ispec   = .false.
         !!!conf_busy_node(n)%sigmaEquiv      = 0.0_RKIND
         !!!conf_busy_node(n)%epsilonRelEquiv = 1.0_RKIND
         !!!conf_busy_node(n)%muRelEquiv      = 1.0_RKIND

!!!      Movido a cada wires.F90. Aqui me parece inreportable y peligroso
!!!!pedazo de niapa para poner los sggmiNo conformal a voltage nulo y que Dios reparta suerte 130220
!!!         !barro el interior, por eso el shifing k,x - j,z - i,y
!!!         Do k=sgg%alloc(iHx)%ZI , sgg%alloc(iHx)%ZE
!!!            Do j=sgg%alloc(iHz)%YI , sgg%alloc(iHz)%YE
!!!               Do i=sgg%alloc(iHy)%XI , sgg%alloc(iHy)%XE
!!!                  medio1 =sggMiEx(i,j,k)
!!!                  medio2 =sggMiEy(i,j,k)
!!!                  medio3 =sggMiEz(i,j,k)
!!!                  mediois1= sgg%med(medio1)%is%already_YEEadvanced_byconformal .or. sgg%med(medio1)%is%split_and_useless
!!!                  mediois2= sgg%med(medio2)%is%already_YEEadvanced_byconformal .or. sgg%med(medio2)%is%split_and_useless
!!!                  mediois3= sgg%med(medio3)%is%already_YEEadvanced_byconformal .or. sgg%med(medio3)%is%split_and_useless
!!!                  if (mediois1.or.mediois2.or.mediois3)  then
!!!                      sggMiNo(i,j,k)=0 !ojo !esto no sirve para contactos conformal SGBC-wires, solo para conformal PEC-wires
!!!                  endif
!!!               End do
!!!            End do
!!!         End do
!!!!fin pedado de niapa


      endif
#endif

      return
   end subroutine AssigLossyOrPECtoNodes

   LOGICAL FUNCTION IGUALES(A,B) RESULT(IGUAL)
      REAL (KIND=RKIND) :: A,B,ERR
      igual=.false.
      IF (abs(A+B)>1e-20 ) THEN
         ERR=2.0_RKIND*ABS((A-B)/(A+B))
         if (err <1e-2_RKIND) igual=.true. !en tanto por ciento me apanio con un 1 por ciento
      ELSE
         ERR=ABS(A-B)
         if (err <1e-20_RKIND) igual=.true. !en valor absoluto para valores casi nulos le pido que el error sea casi nulo
      ENDIF
      RETURN
   END FUNCTION IGUALES


   subroutine populatePlaneWaveRC(Planewave,simu_devia)
      type (planeonde_t)           ::  PlaneWave
      integer :: kkk
      real (kind=8) :: theta, phi, alpha,beta,alpha1,alpha2,amplitud,FACTOR
      complex :: beta1,beta2
      logical :: primeravez,simu_devia
      CHARACTER (LEN=BUFSIZE) :: buff
!integer :: values(1:8), k
!integer, dimension(:), allocatable :: seed
!real(8) :: r
!call date_and_time(values=values)
!call random_seed(size=k)
!allocate(seed(1:k))
!seed(:) = values(8)
!call random_seed(put=seed)

      call random_seed()

      amplitud=1.0
      do kkk=1,PlaneWave%numModes
1        continue !punto de retorno si hay algun error de redondeo !vivan los gotos !!!!
         primeravez=.true.
         theta=0.; phi=0.;
         do while(((2.0_RKIND *pi*sin(theta) < phi)).or.primeravez) !moglie
            primeravez=.false.
            CALL RANDOM_NUMBER(theta)
            theta=pi*theta
            CALL RANDOM_NUMBER(phi)
            phi=2.0_RKIND *pi*phi
         end do !moglie
         phi=phi/sin(theta) !moglie
!ahora la polarizacion
!!!! si los hago asi hay apegotonamiento en los polos 281115 pero con los beta tampoco me sale. Seguir pensando y hacerlo con poincare algun dia 281115
!generado con ortogonalidad_teM_parafuentesRC.nb
!ojo que el atan de fortran y de mathematica estan invertidos!!!!
2        continue
         CALL RANDOM_NUMBER(beta)
         beta=2.0_RKIND *pi*beta
         alpha1=atan2(  Cos(theta)/Sqrt(Cos(theta)**2.0_RKIND+ Cos(beta - phi)**2.0*Sin(theta)**2.0),-((Cos(beta - phi)*Sin(theta))/Sqrt(Cos(theta)**2.0_RKIND+ Cos(beta - phi)**2.0*Sin(theta)**2.0)))
         alpha2=atan2(-(Cos(theta)/Sqrt(Cos(theta)**2.0_RKIND+ Cos(beta - phi)**2.0*Sin(theta)**2.0)), (Cos(beta - phi)*Sin(theta))/Sqrt(Cos(theta)**2.0_RKIND+ Cos(beta - phi)**2.0*Sin(theta)**2.0))
         if ((alpha1 <= pi).and.(alpha1 >= 0.0)) then
            alpha=alpha1
         elseif ((alpha2 <= pi).and.(alpha2 >= 0.0)) then
            alpha=alpha2
         else
            goto 2
            ! WRITE (buff,*) 'Error generando direcciones aleatorias para RC planewaves. '
            ! CALL STOPONERROR(0,0,buff)
         endif
         !!! beta1=atan2(-2.0*1.0/Tan(alpha)*1.0/Tan(theta)*Sin(phi) + Sqrt(2.0)*1.0/Tan(phi)*1.0/Sin(alpha)**2.0*1.0/Sin(theta)**2.0*Sqrt(-((Cos(2.0*alpha) + Cos(2.0*theta))*Sin(alpha)**2.0*Sin(phi)**2.0*Sin(theta)**2.0)), &
         !!!    -(1.0/Sin(alpha)*1.0/Sin(theta)*(2.0*Cos(alpha)*Cos(phi)*Cos(theta) + Sqrt(2.0)*1.0/Sin(alpha)*1.0/Sin(theta)*Sqrt(-((Cos(2.0*alpha) + Cos(2.0*theta))*Sin(alpha)**2.0*Sin(phi)**2.0*Sin(theta)**2.0)))))
         !!! beta2=atan2(-(1.0/Sin(alpha)*1.0/Sin(theta)*(2.0*Cos(alpha)*Cos(theta)*Sin(phi) + Sqrt(2.0)*1.0/Tan(phi)*1.0/Sin(alpha)*1.0/Sin(theta)*Sqrt(-((Cos(2.0*alpha) + Cos(2.0*theta))*Sin(alpha)**2.0*Sin(phi)**2.0*Sin(theta)**2.0)))), &
         !!!      -2.0*Cos(phi)*1.0/Tan(alpha)*1.0/Tan(theta) + Sqrt(2.0)*1.0/Sin(alpha)**2.0*1.0/Sin(theta)**2.0*Sqrt(-((Cos(2.0*alpha) + Cos(2.0*theta))*Sin(alpha)**2.0*Sin(phi)**2.0*Sin(theta)**2.0)))


!!!ahora la incertumbre en la posicion
         CALL RANDOM_NUMBER(factor)
         planewave%INCERT(kkk)=planewave%incertMax*factor
!!!
         !
         planewave%px(kkk) = Sin (theta) * Cos (phi)
         planewave%py(kkk) = Sin (theta) * Sin (phi)
         planewave%pz(kkk) = Cos (theta)
         planewave%ex(kkk) = amplitud * Sin (alpha) * Cos (beta)
         planewave%ey(kkk) = amplitud * Sin (alpha) * Sin (beta)
         planewave%ez(kkk) = amplitud * Cos (alpha)
         !ojo con estos redondeos.
         !!!IF (Abs(planewave%ex(KKK)/amplitud) < 1e-4) planewave%ex(KKK) = 0.0_RKIND
         !!!IF (Abs(planewave%ey(KKK)/amplitud) < 1e-4) planewave%ey(KKK) = 0.0_RKIND
         !!!IF (Abs(planewave%ez(KKK)/amplitud) < 1e-4) planewave%ez(KKK) = 0.0_RKIND
         !!!IF (Abs(planewave%px(KKK)) < 1e-4) planewave%px(KKK) = 0.0_RKIND
         !!!IF (Abs(planewave%py(KKK)) < 1e-4) planewave%py(KKK) = 0.0_RKIND
         !!!IF (Abs(planewave%pz(KKK)) < 1e-4) planewave%pz(KKK) = 0.0_RKIND
         IF (abs(planewave%px(KKK)**2.+planewave%py(KKK)**2.+planewave%pz(KKK)**2.-1.)>1e-4) THEN
!!!            WRITE (buff,*) 'NO TEM PLANEWAVE in RC routine'
            goto 1 !CALL STOPONERROR(0,0,buff)
         END IF
         IF (abs(planewave%ex(KKK)**2.+planewave%ey(KKK)**2.+planewave%ez(KKK)**2.-amplitud**2.)>1e-4) THEN
!!!            WRITE (buff,*) 'NO TEM PLANEWAVE in RC routine'
            goto 1 !CALL STOPONERROR(0,0,buff)
         END IF
         IF (Abs(planewave%px(KKK)*planewave%ex(KKK)+planewave%py(KKK)*planewave%ey(KKK)+planewave%pz(KKK)*planewave%ez(KKK)) >= 1e-4) THEN
!!!            WRITE (buff,*) 'NO TEM PLANEWAVE in RC routine'
            goto 1 !CALL STOPONERROR(0,0,buff)
         END IF
!!!         write (777,'(i5,12e19.9e3)') kkk,theta, phi, alpha, beta,factor
      end do

      open(888,file='rc_EP.dat', FORM='formatted')
      do
         read (888,'(i5,12e19.9e3)',end=888,err=888) kkk, planewave%px(kkk),planewave%py(kkk),planewave%pz(kkk), &
            planewave%ex(kkk),planewave%ey(kkk),planewave%ez(kkk), planewave%INCERT(kkk)
      end do
888   continue
      close(888)

      if (.not.simu_devia) then !solo lo escribe el principal
         open(888,file='rc_EP.dat', FORM='formatted')
         do kkk=1,PlaneWave%numModes
            write (888,'(i5,12e19.9e3)') kkk, planewave%px(kkk),planewave%py(kkk),planewave%pz(kkk), &
               planewave%ex(kkk),planewave%ey(kkk),planewave%ez(kkk),planewave%INCERT(kkk)
         end do
         close(888)
      endif

   end subroutine populatePlaneWaveRC

!!!!       subroutine init_random_seed()
!!!!           use iso_fortran_env, only: int64
!!!!           implicit none
!!!!           integer, allocatable :: seed(:)
!!!!           integer :: i, n, un, istat, dt(8), pid
!!!!           integer(int64) :: t
!!!!
!!!!           call random_seed(size = n)
!!!!           allocate(seed(n))
!!!!           ! First try if the OS provides a random number generator
!!!!           open(newunit=un, file="/dev/urandom", access="stream", &
!!!!                form="unformatted", action="read", status="old", iostat=istat)
!!!!           if (istat == 0) then
!!!!              read(un) seed
!!!!              close(un)
!!!!           else
!!!!              ! Fallback to XOR:ing the current time and pid. The PID is
!!!!              ! useful in case one launches multiple instances of the same
!!!!              ! program in Parallel.
!!!!              call system_clock(t)
!!!!              if (t == 0) then
!!!!                 call date_and_time(values=dt)
!!!!                 t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
!!!!                      + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
!!!!                      + dt(3) * 24_int64 * 60 * 60 * 1000 &
!!!!                      + dt(5) * 60 * 60 * 1000 &
!!!!                      + dt(6) * 60 * 1000 + dt(7) * 1000 &
!!!!                      + dt(8)
!!!!              end if
!!!!              pid = getpid()
!!!!              t = ieor(t, int(pid, kind(t)))
!!!!              do i = 1, n
!!!!                 seed(i) = lcg(t)
!!!!              end do
!!!!           end if
!!!!           call random_seed(put=seed)
!!!!         contains
!!!!           ! This simple PRNG might not be good enough for real work, but is
!!!!           ! sufficient for seeding a better PRNG.
!!!!           function lcg(s)
!!!!             integer :: lcg
!!!!             integer(int64) :: s
!!!!             if (s == 0) then
!!!!                s = 104729
!!!!             else
!!!!                s = mod(s, 4294967296_int64)
!!!!             end if
!!!!             s = mod(s * 279470273_int64, 4294967291_int64)
!!!!             lcg = int(mod(s, int(huge(0), int64)), kind(0))
!!!!           end function lcg
!!!!         end subroutine init_random_seed

   subroutine cuentatags(this,tagtype,layoutnumber,fichin)

      CHARACTER (LEN=*), INTENT (IN) :: fichin
      INTEGER (KIND=4), INTENT (IN) :: layoutnumber

      TYPE (Parseador), INTENT (INOUT) :: this
      LOGICAL :: foundDuplicate
      integer (Kind=4) :: numertag, i,j, k, m, tama,tama2,tama3,tama2p,tama3p,precounting,acum,thefileno
      CHARACTER(LEN=BUFSIZE) :: tagToCheck
      type (tagtype_t) :: tagtype

      !!!ojoo
!!!!return !ojooo

      do precounting=0,1
         numertag=0
         !
         tama = (this%pecregs%nvols)
         DO i = 1, tama
            numertag = numertag + 1;
            if ((i>1)) then
               if ((this%pecregs%vols(i)%tag == this%pecregs%vols(i-1)%tag)) then !do not increase
                  numertag=numertag-1
               endif
            endif
            if (precounting==1) tagtype%tag(numertag) = this%pecregs%vols(i)%tag
         end do
         tama = (this%pecregs%nsurfs)
         DO i = 1, tama
            numertag = numertag + 1;
            if ((i>1)) then
               if ((this%pecregs%surfs(i)%tag == this%pecregs%surfs(i-1)%tag)) then !do not increase
                  numertag=numertag-1
               endif
            endif
            if (precounting==1) tagtype%tag(numertag) = this%pecregs%surfs(i)%tag
         end do
         tama = (this%pecregs%nLINS)
         DO i = 1, tama
            numertag = numertag + 1;
            if ((i>1)) then
               if ((this%pecregs%lins(i)%tag == this%pecregs%lins(i-1)%tag)) then !do not increase
                  numertag=numertag-1
               endif
            endif
            if (precounting==1) tagtype%tag(numertag) = this%pecregs%lins(i)%tag
         enddo
         !
         tama = (this%pmcregs%nvols)
         DO i = 1, tama
            numertag = numertag + 1;
            if ((i>1)) then
               if ((this%pmcregs%vols(i)%tag == this%pmcregs%vols(i-1)%tag)) then !do not increase
                  numertag=numertag-1
               endif
            endif
            if (precounting==1) tagtype%tag(numertag) = this%pmcregs%vols(i)%tag
         end do
         tama = (this%pmcregs%nsurfs)
         DO i = 1, tama
            numertag = numertag + 1;
            if ((i>1)) then
               if ((this%pmcregs%surfs(i)%tag == this%pmcregs%surfs(i-1)%tag)) then !do not increase
                  numertag=numertag-1
               endif
            endif
            if (precounting==1) tagtype%tag(numertag) = this%pmcregs%surfs(i)%tag
         end do
         tama = (this%pmcregs%nLINS)
         DO i = 1, tama
            numertag = numertag + 1;
            if ((i>1)) then
               if ((this%pmcregs%lins(i)%tag == this%pmcregs%lins(i-1)%tag)) then !do not increase
                  numertag=numertag-1
               endif
            endif
            if (precounting==1) tagtype%tag(numertag) = this%pmcregs%lins(i)%tag
         enddo
         !
         !
         !

         tama = (this%DielRegs%nvols)
         DO i = 1, tama
            tama2 = (this%DielRegs%vols(i)%n_c1P)
            tama3 = (this%DielRegs%vols(i)%n_c2P)

            if ((tama2 == 0) .and. (tama3==0)) then
               print *,'Bug in Dielectric Tags. Missing coordinates'
               stop
            end if

            ! Check c1P coordinates
            call checkDielectricComponentTags(this%DielRegs%vols(i), this%DielRegs%vols(1:i-1), i-1, &
               'c1P', numertag, tagtype, precounting, &
               'Bug in Dielectric Volume Tags')

            ! Check c2P coordinates
            call checkDielectricComponentTags(this%DielRegs%vols(i), this%DielRegs%vols(1:i-1), i-1, &
               'c2P', numertag, tagtype, precounting, &
               'Bug in Dielectric Volume Tags')
         end do

         ! Similar for surfaces
         tama = (this%DielRegs%nSurfs)
         DO i = 1, tama
            tama2 = (this%DielRegs%surfs(i)%n_c1P)
            tama3 = (this%DielRegs%surfs(i)%n_c2P)

            if ((tama2 == 0) .and. (tama3==0)) then
               print *,'Bug in Dielectric Tags. Missing coordinates'
               stop
            end if

            ! Check c1P coordinates
            call checkDielectricComponentTags(this%DielRegs%surfs(i), this%DielRegs%surfs(1:i-1), i-1, &
               'c1P', numertag, tagtype, precounting, &
               'Bug in Dielectric Surface Tags')

            ! Check c2P coordinates
            call checkDielectricComponentTags(this%DielRegs%surfs(i), this%DielRegs%surfs(1:i-1), i-1, &
               'c2P', numertag, tagtype, precounting, &
               'Bug in Dielectric Surface Tags')
         end do

         ! Similar for surfaces
         tama = (this%DielRegs%nLins)
         DO i = 1, tama
            tama2 = (this%DielRegs%Lins(i)%n_c1P)
            tama3 = (this%DielRegs%Lins(i)%n_c2P)

            if ((tama2 == 0) .and. (tama3==0)) then
               print *,'Bug in Dielectric Tags. Missing coordinates'
               stop
            end if

            ! Check c1P coordinates
            call checkDielectricComponentTags(this%DielRegs%lins(i), this%DielRegs%lins(1:i-1), i-1, &
               'c2P', numertag, tagtype, precounting, &
               'Bug in Dielectric Surface Tags')

            ! Check c2P coordinates
            call checkDielectricComponentTags(this%DielRegs%lins(i), this%DielRegs%lins(1:i-1), i-1, &
               'c2P', numertag, tagtype, precounting, &
               'Bug in Dielectric Surface Tags')
         end do
!
         tama = (this%animats%nvols)
         DO i = 1, tama
            tama2 = (this%DielRegs%vols(i)%n_c1P)
            tama3 = (this%DielRegs%vols(i)%n_c2P)

            if ((tama2 == 0) .and. (tama3==0)) then
               print *,'Bug in Animat Tags. Missing coordinates'
               stop
            end if

            ! Check c1P coordinates
            call checkAnimatedComponentTags(this%aniMats%vols(i), this%aniMats%vols(1:i-1), i-1, &
               'c1P', numertag, tagtype, precounting, &
               'Bug in Animat Volume Tags')

            ! Check c2P coordinates
            call checkAnimatedComponentTags(this%aniMats%vols(i), this%aniMats%vols(1:i-1), i-1, &
               'c2P', numertag, tagtype, precounting, &
               'Bug in Animat Volume Tags')
         end do
         !
         tama = (this%animats%nSurfs)
         DO i = 1, tama
            tama2 = (this%DielRegs%Surfs(i)%n_c1P)
            tama3 = (this%DielRegs%Surfs(i)%n_c2P)

            if ((tama2 == 0) .and. (tama3==0)) then
               print *,'Bug in Animat Tags. Missing coordinates'
               stop
            end if

            ! Check c1P coordinates
            call checkAnimatedComponentTags(this%aniMats%Surfs(i), this%aniMats%Surfs(1:i-1), i-1, &
               'c1P', numertag, tagtype, precounting, &
               'Bug in Animat Surface Tags')

            ! Check c2P coordinates
            call checkAnimatedComponentTags(this%aniMats%Surfs(i), this%aniMats%Surfs(1:i-1), i-1, &
               'c2P', numertag, tagtype, precounting, &
               'Bug in Animat Surface Tags')
         end do
         !
         tama = (this%animats%nLins)
         DO i = 1, tama
            tama2 = (this%DielRegs%Lins(i)%n_c1P)
            tama3 = (this%DielRegs%Lins(i)%n_c2P)

            if ((tama2 == 0) .and. (tama3==0)) then
               print *,'Bug in Animat Tags. Missing coordinates'
               stop
            end if

            ! Check c1P coordinates
            call checkAnimatedComponentTags(this%aniMats%Lins(i), this%aniMats%Lins(1:i-1), i-1, &
               'c1P', numertag, tagtype, precounting, &
               'Bug in Animat Line Tags')

            ! Check c2P coordinates
            call checkAnimatedComponentTags(this%aniMats%Lins(i), this%aniMats%Lins(1:i-1), i-1, &
               'c2P', numertag, tagtype, precounting, &
               'Bug in Animat Line Tags')
         end do
!
!
!
         tama = (this%frqdepmats%nvols)
         DO i = 1, tama
            numertag = numertag + 1;
            tama2 = (this%frqdepmats%vols(i)%n_c)
            if (tama2/=0) then
               if ((i>1)) then
                  if ((this%frqdepmats%vols(i)%c(1)%tag == this%frqdepmats%vols(i-1)%c(1)%tag)) then !do not increase
                     numertag=numertag-1
                  endif
               endif
            endif
            if (precounting==1) then
               if (tama2/=0) then
                  tagtype%tag(numertag) =  this%frqdepmats%vols(i)%c(1)%tag
               else
                  print *,'bug in tags. '
                  stop
               endif
               DO j = 1, tama2
                  if (trim(adjustl(this%frqdepmats%vols(i)%c(j)%tag)) /= trim(adjustl(tagtype%tag(numertag)))) then
                     print *,'bug in tags. '
                     stop
                  endif
               end do
            endif
         end do
         !
         tama = (this%frqdepmats%nsurfs)
         DO i = 1, tama
            numertag = numertag + 1;
            tama2 = (this%frqdepmats%surfs(i)%n_c)
            if (tama2/=0) then
               if ((i>1)) then
                  if ((this%frqdepmats%surfs(i)%c(1)%tag == this%frqdepmats%surfs(i-1)%c(1)%tag)) then !do not increase
                     numertag=numertag-1
                  endif
               endif
            endif
            if (precounting==1) then
               if (tama2/=0) then
                  tagtype%tag(numertag) =  this%frqdepmats%surfs(i)%c(1)%tag
               else
                  print *,'bug in tags. '
                  stop
               endif
               DO j = 1, tama2
                  if (trim(adjustl(this%frqdepmats%surfs(i)%c(j)%tag)) /= trim(adjustl(tagtype%tag(numertag)))) then
                     print *,'bug in tags. '
                     stop
                  endif
               end do
            endif
         end do
         !
         tama = (this%frqdepmats%nlins)
         DO i = 1, tama
            numertag = numertag + 1;
            tama2 = (this%frqdepmats%lins(i)%n_c)
            if (tama2/=0) then
               if ((i>1)) then
                  if ((this%frqdepmats%lins(i)%c(1)%tag == this%frqdepmats%lins(i-1)%c(1)%tag)) then !do not increase
                     numertag=numertag-1
                  endif
               endif
            endif
            if (precounting==1) then
               if (tama2/=0) then
                  tagtype%tag(numertag) =  this%LossyThinSurfs%cs(i)%C(1)%tag
               else
                  print *,'bug in tags. '
                  stop
               endif
               DO j = 1, tama2
                  if (trim(adjustl(this%frqdepmats%lins(i)%c(j)%tag)) /= trim(adjustl(tagtype%tag(numertag)))) then
                     print *,'bug in tags. '
                     stop
                  endif
               end do
            endif
         end do
!!!
         tama = this%LossyThinSurfs%length
         DO i = 1, tama
            call checkLossyTags(this%LossyThinSurfs%cs(i), &
                              this%LossyThinSurfs%cs(1:i-1), &
                              i-1, numertag, tagtype, precounting)
         END DO

         tama = this%twires%n_tw
         do i=1, tama
            numertag = numertag + 1
            tama2 = this%twires%TW(i)%N_TWC
            if (tama2/=0) then
               if ((i>1)) then
                  if ((this%twires%TW(i)%TWC(1)%tag   == this%twires%TW(i-1)%TWC(1)%tag  )) then !do not increase
                     numertag=numertag-1
                  endif
               endif
            endif
            if (precounting==1) then
               if (tama2/=0) then
                  tagtype%tag(numertag) =  this%twires%TW(i)%TWC(1)%tag
               else
                  print *,'bug in tags. '
                  stop
               endif
               DO j = 1, tama2
                  if (trim(adjustl(this%twires%TW(i)%TWC(j)%tag)) /= trim(adjustl(tagtype%tag(numertag)))) then
                     print *,'bug in tags. '
                     stop
                  endif
               end do
            endif
         end do


         tama = this%swires%n_sw
         do i=1, tama
            numertag = numertag + 1
            tama2=this%swires%SW(i)%n_swc
            if (tama2/=0) then
               if ((i>1)) then
                  if ((this%swires%SW(i)%swc(1)%tag == this%swires%SW(i-1)%swc(1)%tag)) then !do not increase
                     numertag=numertag-1
                  endif
               endif
            endif
            if (precounting==1) then
               if (tama2/=0) then
                  tagtype%tag(numertag) =  this%swires%SW(i)%swc(1)%tag
               else
                  print *,'bug in tags. '
                  stop
               endif
               do j = 1,tama2
                  if (trim(adjustl(this%swires%SW(i)%swc(j)%tag)) /= trim(adjustl(tagtype%tag(numertag)))) then
                     print *,'bug in tags. '
                     stop
                  endif
               end do
            endif
         end do

         tama = this%tSlots%n_tg
         DO i = 1, tama
            numertag = numertag + 1
            tama2 = this%tSlots%Tg(i)%N_tgc
            if (tama2/=0) then
               if ((i>1)) then
                  if ((this%tSlots%Tg(i)%TgC(1)%tag == this%tSlots%Tg(i-1)%TgC(1)%tag)) then !do not increase
                     numertag=numertag-1
                  endif
               endif
            endif
            if (precounting==1) then
               if (tama2/=0) then
                  tagtype%tag(numertag) =  this%tSlots%Tg(i)%TgC(1)%tag
               else
                  print *,'bug in tags. '
                  stop
               endif
               DO j = 1, tama2
                  if (trim(adjustl(this%tSlots%Tg(i)%TgC(j)%tag)) /= trim(adjustl(tagtype%tag(numertag)))) then
                     print *,'bug in tags. '
                     stop
                  endif
               end do
            endif
         end do
!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!
         if (precounting==0) then
            tagtype%numertags = numertag
            allocate(tagtype%tag(1:numertag+1)) !uno mas para luego jugar
            tagtype%tag=''
         else !elimina repetidos
            do i=1,numertag
               do j=i+1,numertag
                  if ((trim(adjustl(tagtype%tag(i)))==trim(adjustl(tagtype%tag(j))))) then
                     tagtype%tag(j)=''
                  endif
               end do
            end do
            i=1
            acum=0
            do while ((i<=numertag).and.(acum<=numertag+1))
               if (trim(adjustl(tagtype%tag(i)))=='') then
                  tagtype%tag(i:numertag)=tagtype%tag(i+1:numertag+1)
                  acum=acum+1
               else
                  i=i+1
               endif
            end do
            numertag=i-1
            tagtype%numertags = numertag
         endif
      end do !del precounting


      return
   end subroutine cuentatags


   SUBROUTINE checkDielectricComponentTags (component, prev_components, n_prev, coord_type, numertag, tagtype, precounting, error_msg)
      type(Dielectric_t), intent(in) :: component        ! Current component
      type(Dielectric_t), intent(in) :: prev_components(:) ! Array of previous components
      integer, intent(in) :: n_prev                            ! Number of previous components
      character(len=*), intent(in) :: coord_type               ! 'c1P' or 'c2P'
      integer, intent(inout) :: numertag
      type(tagtype_t), intent(inout) :: tagtype
      integer, intent(in) :: precounting
      character(len=*), intent(in) :: error_msg

      logical :: foundDuplicate
      character(len=BUFSIZE) :: tagToCheck
      integer :: i, j, k, m, tama2, prev_size

      if (coord_type == 'c1P') then
         tama2 = component%n_c1P
      else
         tama2 = component%n_c2P
      endif

      check_tags: do j = 1, tama2
         numertag = numertag + 1
         foundDuplicate = .false.

         ! Get tag to check based on coord_type
         if (coord_type == 'c1P') then
            tagToCheck = trim(adjustl(component%c1P(j)%tag))
         else
            tagToCheck = trim(adjustl(component%c2P(j)%tag))
         endif

         ! Check current component up to j-1
         if (j > 1) then
            check_current: do k = 1, j-1
               if (coord_type == 'c1P') then
                  if (tagToCheck == trim(adjustl(component%c1P(k)%tag))) then
                     foundDuplicate = .true.
                     exit check_current
                  endif
               else
                  if (tagToCheck == trim(adjustl(component%c2P(k)%tag))) then
                     foundDuplicate = .true.
                     exit check_current
                  endif
               endif
            end do check_current
         endif

         ! If not found, check all previous components
         if (.not. foundDuplicate) then
            check_previous: do m = 1, n_prev
               ! Check c1P of previous component
               if (prev_components(m)%n_c1P > 0) then
                  do k = 1, prev_components(m)%n_c1P
                     if (tagToCheck == trim(adjustl(prev_components(m)%c1P(k)%tag))) then
                        print *, error_msg
                        print *, 'Duplicate tag found:', tagToCheck
                        stop
                     endif
                  end do
               endif

               ! Check c2P of previous component
               if (prev_components(m)%n_c2P > 0) then
                  do k = 1, prev_components(m)%n_c2P
                     if (tagToCheck == trim(adjustl(prev_components(m)%c2P(k)%tag))) then
                        print *, error_msg
                        print *, 'Duplicate tag found:', tagToCheck
                        stop
                     endif
                  end do
               endif
            end do check_previous
         endif

         if (foundDuplicate) then
            numertag = numertag - 1
         else if (precounting == 1) then
            tagtype%tag(numertag) = tagToCheck
         end if
      end do check_tags
   END SUBROUTINE


   SUBROUTINE checkAnimatedComponentTags (component, prev_components, n_prev, coord_type, numertag, tagtype, precounting, error_msg)
      type(ANISOTROPICbody_t), intent(in) :: component        ! Current component
      type(ANISOTROPICbody_t), intent(in) :: prev_components(:) ! Array of previous components
      integer, intent(in) :: n_prev                            ! Number of previous components
      character(len=*), intent(in) :: coord_type               ! 'c1P' or 'c2P'
      integer, intent(inout) :: numertag
      type(tagtype_t), intent(inout) :: tagtype
      integer, intent(in) :: precounting
      character(len=*), intent(in) :: error_msg

      logical :: foundDuplicate
      character(len=BUFSIZE) :: tagToCheck
      integer :: i, j, k, m, tama2, prev_size

      if (coord_type == 'c1P') then
         tama2 = component%n_c1P
      else
         tama2 = component%n_c2P
      endif

      check_tags: do j = 1, tama2
         numertag = numertag + 1
         foundDuplicate = .false.

         ! Get tag to check based on coord_type
         if (coord_type == 'c1P') then
            tagToCheck = trim(adjustl(component%c1P(j)%tag))
         else
            tagToCheck = trim(adjustl(component%c2P(j)%tag))
         endif

         ! Check current component up to j-1
         if (j > 1) then
            check_current: do k = 1, j-1
               if (coord_type == 'c1P') then
                  if (tagToCheck == trim(adjustl(component%c1P(k)%tag))) then
                     foundDuplicate = .true.
                     exit check_current
                  endif
               else
                  if (tagToCheck == trim(adjustl(component%c2P(k)%tag))) then
                     foundDuplicate = .true.
                     exit check_current
                  endif
               endif
            end do check_current
         endif

         ! If not found, check all previous components
         if (.not. foundDuplicate) then
            check_previous: do m = 1, n_prev
               ! Check c1P of previous component
               if (prev_components(m)%n_c1P > 0) then
                  do k = 1, prev_components(m)%n_c1P
                     if (tagToCheck == trim(adjustl(prev_components(m)%c1P(k)%tag))) then
                        print *, error_msg
                        print *, 'Duplicate tag found:', tagToCheck
                        stop
                     endif
                  end do
               endif

               ! Check c2P of previous component
               if (prev_components(m)%n_c2P > 0) then
                  do k = 1, prev_components(m)%n_c2P
                     if (tagToCheck == trim(adjustl(prev_components(m)%c2P(k)%tag))) then
                        print *, error_msg
                        print *, 'Duplicate tag found:', tagToCheck
                        stop
                     endif
                  end do
               endif
            end do check_previous
         endif

         if (foundDuplicate) then
            numertag = numertag - 1
         else if (precounting == 1) then
            tagtype%tag(numertag) = tagToCheck
         end if
      end do check_tags
   END SUBROUTINE

   SUBROUTINE checkLossyTags(component, prev_components, n_prev, numertag, tagtype, precounting)
      type(LossyThinSurface), intent(in) :: component        ! Current component
      type(LossyThinSurface), intent(in) :: prev_components(:) ! Array of previous components
      integer, intent(in) :: n_prev                         ! Number of previous components
      integer, intent(inout) :: numertag
      type(tagtype_t), intent(inout) :: tagtype
      integer, intent(in) :: precounting

      logical :: foundDuplicate
      character(len=BUFSIZE) :: tagToCheck
      integer :: i, j, k, m, tama2

      tama2 = component%nc
      if (tama2 == 0) then
         print *, 'Bug in LossyThinSurf Tags. Missing coordinates'
         stop
      endif

      check_tags: do j = 1, tama2
         numertag = numertag + 1
         foundDuplicate = .false.

         tagToCheck = trim(adjustl(component%C(j)%tag))

         ! Check current component up to j-1
         if (j > 1) then
            check_current: do k = 1, j-1
               if (tagToCheck == trim(adjustl(component%C(k)%tag))) then
                  foundDuplicate = .true.
                  exit check_current
               endif
            end do check_current
         endif

         ! If not found, check all previous components
         
         if ((.not. foundDuplicate) .and. (n_prev>0)) then
            check_previous: do m = 1, n_prev
               if (prev_components(m)%nc > 0) then
                  do k = 1, prev_components(m)%nc
                     if (tagToCheck == trim(adjustl(prev_components(m)%C(k)%tag))) then
                        foundDuplicate = .true.
                     endif
                  end do
               endif
            end do check_previous
         endif

         if (foundDuplicate) then
            numertag = numertag - 1
         else if (precounting == 1) then
            tagtype%tag(numertag) = tagToCheck
         end if
      end do check_tags
   END SUBROUTINE


   function searchtag(tagtype,tag) result(numertag)


      CHARACTER (LEN=BUFSIZE) :: tag
      integer (Kind=4) :: i,numertag
      type (tagtype_t) :: tagtype

      numertag=-1
      busca: do i=1,tagtype%numertags
         if (trim(adjustl(tagtype%tag(i)))==trim(adjustl(tag))) then
            numertag=i
            exit busca
         endif
      end do busca
      return
   end function searchtag


END MODULE Preprocess_m
