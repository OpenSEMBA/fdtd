MODULE VTK
      !
   USE fdetypes
   USE Observa
   use report
   !
   !
   !
   !
   IMPLICIT NONE
   !
   PRIVATE
   PUBLIC createVTK,createVTKOnTheFly
CONTAINS
   !Subrutine to parse the volumic probes to create VTK files on PEC and on wires
   !
   SUBROUTINE createVTK (layoutnumber, size, sgg,vtkindex,somethingdone,mpidir,tagtype,sggMtag,dontwritevtk)
   
   
      type (SGGFDTDINFO), intent(IN)   :: sgg
      INTEGER (KIND=IKINDMTAG), intent(in) :: sggMtag  (sgg%Alloc(iHx)%XI:sgg%Alloc(iHx)%XE, sgg%Alloc(iHy)%YI:sgg%Alloc(iHy)%YE, sgg%Alloc(iHz)%ZI:sgg%Alloc(iHz)%ZE)
      type (tagtype_t) :: tagtype
      integer (KIND=4) :: mpidir
      logical :: vtkindex,yacreado,dontwritevtk
      !------------------------>
      CHARACTER (LEN=BUFSIZE) :: filename ! File name
      CHARACTER (LEN=BUFSIZE) :: fichero,fichero_input,char_i_sub_time ! File name
      !
      !
      !

      INTEGER (KIND=4), INTENT (IN) :: layoutnumber, size
      INTEGER (KIND=4) ::  ierr,  posicionMPI,conta,ecurrentType,eei,eej,eek,esggMtag
      INTEGER (KIND=4) , allocatable , dimension(:) ::   sizeofvalores,NewsizeOfValores

      real (kind=RKIND) :: time,rdum
      !
      !
      TYPE (output_t), POINTER, DIMENSION (:) :: output
      INTEGER (KIND=4) :: iroot
      !
#ifdef CompileWithMPI
      type (Serialized_t)  ::  NewSerialized !para sondas Volumic
#endif
      type (Serialized_t)  ::  Serialized !para almecenar valores serializados en volumenes en vez de Bloque
      INTEGER (KIND=4) , dimension (:) , allocatable :: PosiMPI,NewPosiMPI
      INTEGER (KIND=4) :: indi,numberOfSerialized
      real (  KINd=RKIND), ALLOCATABLE, DIMENSION (:) :: att  
      real (  KINd=RKIND) :: att_rkind
      real (  KINd=RKIND_tiempo) :: att_rkind_tiempo
      !
      INTEGER (KIND=4) :: ii, i1, finalstep
      LOGICAL :: lexis,freqdomain,somethingdone
      character (LEN=BUFSIZE)     ::  dubuf
      INTEGER (KIND=4) :: minXabs, maxXabs, minYabs, maxYabs, minZabs, maxZabs
      CHARACTER (LEN=BUFSIZE) :: pathroot
      character (LEN=BUFSIZE)  ::  chari,charj,chark,chari2,charj2,chark2
      character (LEN=BUFSIZE)  ::  extpoint
      character(len=BUFSIZE) :: buff
      CHARACTER (LEN=BUFSIZE) :: charc
      CHARACTER (LEN=BUFSIZE) :: tag
      !
      CHARACTER (LEN=BUFSIZE) :: whoami, whoamishort
      INTEGER (kind=4)::  numNodes,numEdges,numQuads , iroot2,iroot1,i_sub_time, total_sub_times
      INTEGER (kind=4), parameter :: time_phases_param=35
      real (kind= RKIND), allocatable, dimension(:,:) :: Nodes
      integer (kind=4), allocatable, dimension(:,:) :: Elems
      integer (kind=4) :: coldummy
!      print *,'RKIND,CKIND,REALSIZE,COMPLEXSIZE,MPI_DOUBLE_PRECISION, MPI_DOUBLE_COMPLEX',RKIND,CKIND,REALSIZE,COMPLEXSIZE,MPI_DOUBLE_PRECISION, MPI_DOUBLE_COMPLEX
      yacreado=.false.
      numNodes=0; numEdges=0;numQuads=0;

      WRITE (whoamishort, '(i5)') layoutnumber + 1
      WRITE (whoami, '(a,i5,a,i5,a)') '(', layoutnumber + 1, '/', size, ') '
      !
      output => GetOutput ()!get the output private info from observation
      !
      somethingdone=.false.
      barridoprobes: DO ii = 1, sgg%NumberRequest

         IF (sgg%observation(ii)%Volumic) then
         if (sgg%observation(ii)%nP == 1) then
         if ((sgg%observation(ii)%P(1)%What == iCur).or.(sgg%observation(ii)%P(1)%What == iCurX).or. &
               (sgg%observation(ii)%P(1)%What == iCurY).or.(sgg%observation(ii)%P(1)%What == iCurZ).or. &
               (sgg%observation(ii)%P(1)%What == mapvtk)) THEN !solo corrientes volumicas
               if (sgg%Observation(ii)%done.and.(sgg%Observation(ii)%flushed)) then
                  cycle barridoprobes
               elseif (sgg%Observation(ii)%done) then
                  sgg%Observation(ii)%flushed=.true. !ultima que se flushea
                  continue
               elseif ((.not.(sgg%Observation(ii)%done)).and.(sgg%Observation(ii)%Begun)) then
                  continue
               elseif (.not.(sgg%Observation(ii)%begun))  then
                  cycle barridoprobes
               else !creo que tengo toda la casuistica, por si se me escapa algo continuo, y ya debajo se manejara
                  continue
               endif
         else 
               cycle barridoprobes
         endif
         endif 
         endif
         !sondas Volumic traducelas a VTK
         IF (sgg%observation(ii)%Volumic) then
            if (sgg%observation(ii)%nP == 1) then
               if ((sgg%observation(ii)%P(1)%What == iCur).or.(sgg%observation(ii)%P(1)%What == iCurX).or. &
               (sgg%observation(ii)%P(1)%What == iCurY).or.(sgg%observation(ii)%P(1)%What == iCurZ).or. &
               (sgg%observation(ii)%P(1)%What == mapvtk)) THEN !solo corrientes volumicas
                  INQUIRE (FILE=trim(adjustl(output(ii)%item(1)%path)), EXIST=lexis)
                  if ((lexis).and.(output(ii)%TimesWritten/=0)) then
                     !


                     minXabs = sgg%observation(ii)%P(1)%XI
                     maxXabs = sgg%observation(ii)%P(1)%XE
                     minYabs = sgg%observation(ii)%P(1)%YI
                     maxYabs = sgg%observation(ii)%P(1)%YE
#ifdef CompileWithMPI
                     minZabs = output(ii)%item(1)%ZIorig
                     maxZabs = output(ii)%item(1)%ZEorig
#else
                     minZabs = sgg%observation(ii)%P(1)%zI
                     maxZabs = sgg%observation(ii)%P(1)%zE
#endif
                     write(chari ,'(i7)') minXabs
                     write(charj ,'(i7)') minYabs
                     write(chark ,'(i7)') minZabs
                     write(chari2,'(i7)') maxXabs
                     write(charj2,'(i7)') maxYabs
                     write(chark2,'(i7)') maxZabs
                     !mpidir 190319    !desrotacion para que los nombres sean correctos 
                        if      (mpidir==3) then
                           extpoint=trim(adjustl(chari)) //'_'//trim(adjustl(charj)) //'_'//trim(adjustl(chark))//'__'// &
                                    trim(adjustl(chari2))//'_'//trim(adjustl(charj2))//'_'//trim(adjustl(chark2))
                        elseif  (mpidir==2) then
                           extpoint=trim(adjustl(charj)) //'_'//trim(adjustl(chark)) //'_'//trim(adjustl(chari))//'__'// &
                                    trim(adjustl(charj2))//'_'//trim(adjustl(chark2))//'_'//trim(adjustl(chari2))
                        elseif  (mpidir==1) then
                           extpoint=trim(adjustl(chark)) //'_'//trim(adjustl(chari)) //'_'//trim(adjustl(charj))//'__'// &
                                    trim(adjustl(chark2))//'_'//trim(adjustl(chari2))//'_'//trim(adjustl(charj2))
                        else
                           call stoponerror(layoutnumber,size,'Buggy error in mpidir. ')
                        endif
                     !fin mpidir
                     !
                     iroot=index(output(ii)%item(1)%path,'__',.true.)
                     pathroot=trim(adjustl(output(ii)%item(1)%path(1:iroot-1)))
                     iroot = index (pathroot, '_',.true.)
                     pathroot = trim (adjustl(pathroot(1:iroot-1)))
                     iroot = index (pathroot, '_',.true.)
                     pathroot = trim (adjustl(pathroot(1:iroot-1)))
                     iroot = index (pathroot, '_',.true.)
                     pathroot = trim(adjustl(pathroot(1:iroot-1)))//'_'//trim(adjustl(extpoint))
                     filename = trim (adjustl(pathroot))
                     !


#ifdef CompileWithMPI
                     if (size>1) then
                        if (output(ii)%item(1)%MPISubComm /= -1) then
                           !!! CALL print11 (layoutnumber, trim(adjustl(whoami))////' Init processing file '//trim(adjustl(filename)), .TRUE.) !enforces print
                           continue
                        endif
                     endif
#endif


                     !!!                    OPEN (output(ii)%item(1)%UNIT, FILE=trim(adjustl(output(ii)%item(1)%path)), &
                     !!!                            FORM='unformatted')
                     !!!                        DO
                     !!!                            READ (output(ii)%item(1)%UNIT, end=762) finalstep
                     !!!                        END DO
                     !!!762                     CONTINUE
                     !!!                    CLOSE (output(ii)%item(1)%UNIT)
                     finalstep=output(ii)%TimesWritten
                     allocate (att(1:finalstep))
                     !!!!!!!!!!!!!
                     numberOfSerialized=0
                     allocate (sizeOfValores(0:size-1))
                     sizeOfValores=0
                     sizeofvalores(layoutnumber) = output(ii)%item(1)%columnas
                     !SINCRONIZA EL TAMANIO DE CADA LAYER
#ifdef CompileWithMPI
                     if (size>1) then
                        allocate (NewsizeOfValores(0:size-1))
                        NewsizeOfValores=0
                        if (output(ii)%item(1)%MPISubComm /= -1) then
                           call MPI_Barrier(output(ii)%item(1)%MPISubComm,ierr)
                           CALL MPI_AllReduce (sizeofvalores, newSizeofvalores, size, MPI_INTEGER, MPI_SUM, &
                           &                     output(ii)%item(1)%MPISubComm, ierr)
                        endif
                        sizeofvalores = newSizeofvalores
                     endif
#endif
                     !
                     do i1=0,size-1
                        numberOfSerialized=numberOfSerialized + sizeofvalores(i1)
                     end do
                     !asumo solamente un time step por lectura
                     ALLOCATE (PosiMPI(1:numberOfSerialized))

                     if (SGG%Observation(ii)%TimeDomain) then            
                        ALLOCATE (Serialized%Valor(1,1:numberOfSerialized), &
                                    Serialized%Valor_x(1,1:numberOfSerialized), &
                                    Serialized%Valor_y(1,1:numberOfSerialized), &
                                    Serialized%Valor_z(1,1:numberOfSerialized) )
                        ALLOCATE (Serialized%ValorE(1,1:numberOfSerialized), &
                                    Serialized%Valor_Ex(1,1:numberOfSerialized), &
                                    Serialized%Valor_Ey(1,1:numberOfSerialized), &
                                    Serialized%Valor_Ez(1,1:numberOfSerialized) )
                        ALLOCATE (Serialized%ValorH(1,1:numberOfSerialized), &
                                    Serialized%Valor_Hx(1,1:numberOfSerialized), &
                                    Serialized%Valor_Hy(1,1:numberOfSerialized), &
                                    Serialized%Valor_Hz(1,1:numberOfSerialized) )  
                        Serialized%Valor = 0.   
                        Serialized%Valor_x= 0.
                        Serialized%Valor_y= 0.
                        Serialized%Valor_z= 0. 
                        Serialized%ValorE = 0.   
                        Serialized%Valor_Ex= 0.
                        Serialized%Valor_Ey= 0.
                        Serialized%Valor_Ez= 0. 
                        Serialized%ValorH = 0.   
                        Serialized%Valor_Hx= 0.
                        Serialized%Valor_Hy= 0.
                        Serialized%Valor_Hz= 0.
                        freqdomain=.false.
                     elseif (SGG%Observation(ii)%FreqDomain) then        
                        ALLOCATE (Serialized%Valor(1,1:numberOfSerialized), &
                                    Serialized%Valor_x(1,1:numberOfSerialized), &
                                    Serialized%Valor_y(1,1:numberOfSerialized), &
                                    Serialized%Valor_z(1,1:numberOfSerialized) ) !auxiliar para sincronizar MPI
                        ALLOCATE (Serialized%ValorE(1,1:numberOfSerialized), &
                                    Serialized%Valor_Ex(1,1:numberOfSerialized), &
                                    Serialized%Valor_Ey(1,1:numberOfSerialized), &
                                    Serialized%Valor_Ez(1,1:numberOfSerialized) ) !auxiliar para sincronizar MPI
                        ALLOCATE (Serialized%ValorH(1,1:numberOfSerialized), &
                                    Serialized%Valor_Hx(1,1:numberOfSerialized), &
                                    Serialized%Valor_Hy(1,1:numberOfSerialized), &
                                    Serialized%Valor_Hz(1,1:numberOfSerialized) ) !auxiliar para sincronizar MPI     
                        Serialized%Valor = 0.
                        Serialized%Valor_x = 0.
                        Serialized%Valor_y = 0.
                        Serialized%Valor_z = 0.   
                        Serialized%ValorE = 0.
                        Serialized%Valor_Ex = 0.
                        Serialized%Valor_Ey = 0.
                        Serialized%Valor_Ez = 0.   
                        Serialized%ValorH = 0.
                        Serialized%Valor_Hx = 0.
                        Serialized%Valor_Hy = 0.
                        Serialized%Valor_Hz = 0.      
                        ALLOCATE (Serialized%ValorComplex_x(1,1:numberOfSerialized), &
                                    Serialized%ValorComplex_y(1,1:numberOfSerialized), &
                                    Serialized%ValorComplex_z(1,1:numberOfSerialized) )
                        ALLOCATE (Serialized%ValorComplex_Ex(1,1:numberOfSerialized), &
                                    Serialized%ValorComplex_Ey(1,1:numberOfSerialized), &
                                    Serialized%ValorComplex_Ez(1,1:numberOfSerialized) )
                        ALLOCATE (Serialized%ValorComplex_Hx(1,1:numberOfSerialized), &
                                    Serialized%ValorComplex_Hy(1,1:numberOfSerialized), &
                                    Serialized%ValorComplex_Hz(1,1:numberOfSerialized) )  
                        Serialized%ValorComplex_x = 0.
                        Serialized%ValorComplex_y = 0.
                        Serialized%ValorComplex_z = 0.
                        Serialized%ValorComplex_Ex = 0.
                        Serialized%ValorComplex_Ey = 0.
                        Serialized%ValorComplex_Ez = 0.
                        Serialized%ValorComplex_Hx = 0.
                        Serialized%ValorComplex_Hy = 0.
                        Serialized%ValorComplex_Hz = 0.
                        freqdomain=.true.
                     endif
                     allocate (Serialized%eI(1:numberOfSerialized))
                     allocate (Serialized%eJ(1:numberOfSerialized))
                     allocate (Serialized%eK(1:numberOfSerialized))
                     allocate (Serialized%currentType(1:numberOfSerialized))
                     allocate (Serialized%sggMtag(1:numberOfSerialized))
                     PosiMPI=0
                     Serialized%eI = 0
                     Serialized%eJ = 0
                     Serialized%eK = 0
                     Serialized%currentType = 0
                     Serialized%sggMtag = 0


                     !!!BUSCA LA POSICION mpi E INICIALIZA LOS NUEVOS
                     posicionMPI=0
#ifdef CompileWithMPI
                     if (size>1) then
                        if (output(ii)%item(1)%MPISubComm /= -1) then
                           buscaMPI: do i1=0,layoutnumber-1
                              posicionMPI=posicionMPI+sizeofvalores(i1)
                           end do buscaMPI
                        endif
                        ALLOCATE (newPosiMPI(1:numberOfSerialized))
                        if (SGG%Observation(ii)%TimeDomain) then
                           ALLOCATE (NewSerialized%Valor(1,1:numberOfSerialized))
                           ALLOCATE (NewSerialized%Valor_x(1,1:numberOfSerialized))
                           ALLOCATE (NewSerialized%Valor_y(1,1:numberOfSerialized))
                           ALLOCATE (NewSerialized%Valor_z(1,1:numberOfSerialized))   
                           
                           ALLOCATE (NewSerialized%ValorE(1,1:numberOfSerialized))
                           ALLOCATE (NewSerialized%Valor_Ex(1,1:numberOfSerialized))
                           ALLOCATE (NewSerialized%Valor_Ey(1,1:numberOfSerialized))
                           ALLOCATE (NewSerialized%Valor_Ez(1,1:numberOfSerialized))   
                           ALLOCATE (NewSerialized%ValorH(1,1:numberOfSerialized))
                           ALLOCATE (NewSerialized%Valor_Hx(1,1:numberOfSerialized))
                           ALLOCATE (NewSerialized%Valor_Hy(1,1:numberOfSerialized))
                           ALLOCATE (NewSerialized%Valor_Hz(1,1:numberOfSerialized))   
                           
                           NewSerialized%Valor = 0.
                           NewSerialized%Valor_x = 0.
                           NewSerialized%Valor_y = 0.
                           NewSerialized%Valor_z = 0.
                           
                           
                           NewSerialized%ValorE = 0.
                           NewSerialized%Valor_Ex = 0.
                           NewSerialized%Valor_Ey = 0.
                           NewSerialized%Valor_Ez = 0.
                           
                           NewSerialized%ValorH = 0.
                           NewSerialized%Valor_Hx = 0.
                           NewSerialized%Valor_Hy = 0.
                           NewSerialized%Valor_Hz = 0.
                        elseif (SGG%Observation(ii)%FreqDomain) then
!!!lo hago con parte real e imaginaria
                           ALLOCATE (NewSerialized%Valor  (1,1:numberOfSerialized))
                           ALLOCATE (NewSerialized%Valor_x(1,1:numberOfSerialized))
                           ALLOCATE (NewSerialized%Valor_y(1,1:numberOfSerialized))
                           ALLOCATE (NewSerialized%Valor_z(1,1:numberOfSerialized))     
                           
                           ALLOCATE (NewSerialized%ValorE  (1,1:numberOfSerialized))
                           ALLOCATE (NewSerialized%Valor_Ex(1,1:numberOfSerialized))
                           ALLOCATE (NewSerialized%Valor_Ey(1,1:numberOfSerialized))
                           ALLOCATE (NewSerialized%Valor_Ez(1,1:numberOfSerialized))     
                           ALLOCATE (NewSerialized%ValorH  (1,1:numberOfSerialized))
                           ALLOCATE (NewSerialized%Valor_Hx(1,1:numberOfSerialized))
                           ALLOCATE (NewSerialized%Valor_Hy(1,1:numberOfSerialized))
                           ALLOCATE (NewSerialized%Valor_Hz(1,1:numberOfSerialized))     
                           NewSerialized%Valor = 0.
                           NewSerialized%Valor_x = 0.
                           NewSerialized%Valor_y = 0.
                           NewSerialized%Valor_z = 0.  
                              
                           NewSerialized%ValorE = 0.
                           NewSerialized%Valor_Ex = 0.
                           NewSerialized%Valor_Ey = 0.
                           NewSerialized%Valor_Ez = 0.    
                           NewSerialized%ValorH = 0.
                           NewSerialized%Valor_Hx = 0.
                           NewSerialized%Valor_Hy = 0.
                           NewSerialized%Valor_Hz = 0.  
                        endif
                        allocate (newSerialized%eI(1:numberOfSerialized))
                        allocate (newSerialized%eJ(1:numberOfSerialized))
                        allocate (newSerialized%eK(1:numberOfSerialized))
                        allocate (newSerialized%currentType(1:numberOfSerialized))
                        allocate (newSerialized%sggMtag(1:numberOfSerialized))
                        NewPosiMPI=0
                        newSerialized%eI = 0
                        newSerialized%eJ = 0
                        newSerialized%eK = 0
                        newSerialized%currentType = 0
                        newSerialized%sggMtag = 0
                        !
                     endif
#endif
                     !LEE INFO GEOMETRICA TENIENDO EN CUENTA POSICION MPI
                     OPEN (output(ii)%item(1)%UNIT, FILE=trim(adjustl(output(ii)%item(1)%path)), FORM='unformatted')
                     read(output(ii)%item(1)%unit) coldummy
                     if (coldummy/= output(ii)%item(1)%columnas) then
                           write (buff,'(a,2i9)') 'ERROR: Buggy error creating .vtk',coldummy, output(ii)%item(1)%columnas
                           CALL print11(0_4, buff)
                     endif
                     do conta=1,output(ii)%item(1)%columnas
                        read(output(ii)%item(1)%unit) eei,eej,eek,ecurrentType,esggMtag
                        PosiMPI(posicionMPI+conta)=posicionMPI+conta
                        Serialized%eI(posicionMPI+conta)=eeI
                        Serialized%eJ(posicionMPI+conta)=eeJ
                        Serialized%eK(posicionMPI+conta)=eeK
                        Serialized%currentType(posicionMPI+conta)=ecurrentType
                        Serialized%sggMtag(posicionMPI+conta)=esggMtag
                     end do
                     if (SGG%Observation(ii)%FreqDomain) read(output(ii)%item(1)%unit) rdum !instante en el que se ha escrito la info frequencial
                     !SINCRONIZA SUMPANDO LA INFO GEOMETRICA DE TODOS LOS LAYERS
#ifdef CompileWithMPI
                     if (size>1) then
                        newPosiMPI=-1
                        if (output(ii)%item(1)%MPISubComm /= -1) then
                           call MPI_Barrier(output(ii)%item(1)%MPISubComm,ierr)
                           CALL MPI_AllReduce (PosiMPI, newPosiMPI, numberOfSerialized, MPI_INTEGER, MPI_SUM, &
                           &                     output(ii)%item(1)%MPISubComm, ierr)
                        endif
                        PosiMPI = newPosiMPI
                        !
                        if (output(ii)%item(1)%MPISubComm /= -1) then
                           call MPI_Barrier(output(ii)%item(1)%MPISubComm,ierr)
                           CALL MPI_AllReduce (Serialized%eI, newSerialized%eI, numberOfSerialized, MPI_INTEGER, MPI_SUM, &
                           &                     output(ii)%item(1)%MPISubComm, ierr)
                        endif
                        Serialized%eI = newSerialized%eI
                        !
                        if (output(ii)%item(1)%MPISubComm /= -1) then
                           call MPI_Barrier(output(ii)%item(1)%MPISubComm,ierr)
                           CALL MPI_AllReduce (Serialized%eJ, newSerialized%eJ, numberOfSerialized, MPI_INTEGER, MPI_SUM, &
                           &                     output(ii)%item(1)%MPISubComm, ierr)
                        endif
                        Serialized%eJ = newSerialized%eJ
                        !
                        if (output(ii)%item(1)%MPISubComm /= -1) then
                           call MPI_Barrier(output(ii)%item(1)%MPISubComm,ierr)
                           CALL MPI_AllReduce (Serialized%eK, newSerialized%eK, numberOfSerialized, MPI_INTEGER, MPI_SUM, &
                           &                     output(ii)%item(1)%MPISubComm, ierr)
                        endif
                        Serialized%eK = newSerialized%eK
                        !
                        if (output(ii)%item(1)%MPISubComm /= -1) then
                           call MPI_Barrier(output(ii)%item(1)%MPISubComm,ierr)
                           CALL MPI_AllReduce (Serialized%currentType, newSerialized%currentType, numberOfSerialized, MPI_INTEGER, MPI_SUM, &
                           &                     output(ii)%item(1)%MPISubComm, ierr)
                        endif
                        Serialized%currentType = newSerialized%currentType
                        !
                        if (output(ii)%item(1)%MPISubComm /= -1) then
                           call MPI_Barrier(output(ii)%item(1)%MPISubComm,ierr)
                           CALL MPI_AllReduce (Serialized%sggMtag, newSerialized%sggMtag, numberOfSerialized, MPI_INTEGER, MPI_SUM, &
                           &                     output(ii)%item(1)%MPISubComm, ierr)
                        endif
                        Serialized%sggMtag = newSerialized%sggMtag
                     endif
#endif

                     !crea informacion unstruct y escribela en el fichero
#ifdef CompileWithMPI

                     IF (layoutnumber == output(ii)%item(1)%MPIRoot) THEN
#else
                     IF (layoutnumber == 0) THEN
#endif
                        call creaUnstructData(Serialized,  numberOfSerialized,sgg,Nodes,NumNodes,Elems,NumEdges,NumQuads,vtkindex)
                     endif

                     !LEE CADA TIME STEPO TENIENDO EN CUENTA POSICION MPI

                     bucleindi: DO indi = 1, finalstep
                        if (SGG%Observation(ii)%TimeDomain) then
                           Serialized%Valor = 0.
                           Serialized%Valor_x = 0.
                           Serialized%Valor_y = 0.
                           Serialized%Valor_z = 0.
                           
                           Serialized%ValorE = 0.
                           Serialized%Valor_Ex = 0.
                           Serialized%Valor_Ey = 0.
                           Serialized%Valor_Ez = 0.
                           Serialized%ValorH = 0.
                           Serialized%Valor_Hx = 0.
                           Serialized%Valor_Hy = 0.
                           Serialized%Valor_Hz = 0.
                           READ (output(ii)%item(1)%UNIT) att_rkind_tiempo  
                           att(indi)=att_rkind_tiempo
                           if (output(ii)%item(1)%columnas /=0) then
                                 do   conta=1,output(ii)%item(1)%columnas
                                    READ (output(ii)%item(1)%UNIT)           Serialized%valor  (1,posicionMPI+conta), &
                                                                              Serialized%valor_x(1,posicionMPI+conta), &
                                                                              Serialized%valor_y(1,posicionMPI+conta), &
                                                                              Serialized%valor_z(1,posicionMPI+conta) !lo meto en el unico step
                                    
                                    READ (output(ii)%item(1)%UNIT)           Serialized%valorE  (1,posicionMPI+conta), &
                                                                              Serialized%valor_Ex(1,posicionMPI+conta), &
                                                                              Serialized%valor_Ey(1,posicionMPI+conta), &
                                                                              Serialized%valor_Ez(1,posicionMPI+conta) !lo meto en el unico step
                                    
                                    READ (output(ii)%item(1)%UNIT)           Serialized%valorH  (1,posicionMPI+conta), &
                                                                              Serialized%valor_Hx(1,posicionMPI+conta), &
                                                                              Serialized%valor_Hy(1,posicionMPI+conta), &
                                                                              Serialized%valor_Hz(1,posicionMPI+conta) !lo meto en el unico step
                                 end do
                           endif
                        elseif (SGG%Observation(ii)%FreqDomain) then    
                           Serialized%ValorComplex_x = 0.
                           Serialized%ValorComplex_y = 0.
                           Serialized%ValorComplex_z = 0.
                           Serialized%ValorComplex_Ex = 0.
                           Serialized%ValorComplex_Ey = 0.
                           Serialized%ValorComplex_Ez = 0.
                           Serialized%ValorComplex_Hx = 0.
                           Serialized%ValorComplex_Hy = 0.
                           Serialized%ValorComplex_Hz = 0.
                           READ (output(ii)%item(1)%UNIT) att_rkind  
                           att(indi)=att_rkind
                           if (output(ii)%item(1)%columnas /=0) then
                                    do conta=1,output(ii)%item(1)%columnas   
                                       READ (output(ii)%item(1)%UNIT) &
                                          Serialized%valorComplex_x(1,posicionMPI+conta), &
                                          Serialized%valorComplex_y(1,posicionMPI+conta), &
                                          Serialized%valorComplex_z(1,posicionMPI+conta) !lo meto en el unico step
                                       READ (output(ii)%item(1)%UNIT) &
                                          Serialized%valorComplex_Ex(1,posicionMPI+conta), &
                                          Serialized%valorComplex_Ey(1,posicionMPI+conta), &
                                          Serialized%valorComplex_Ez(1,posicionMPI+conta) !lo meto en el unico step
                                       READ (output(ii)%item(1)%UNIT) &
                                          Serialized%valorComplex_Hx(1,posicionMPI+conta), &
                                          Serialized%valorComplex_Hy(1,posicionMPI+conta), &
                                          Serialized%valorComplex_Hz(1,posicionMPI+conta) !lo meto en el unico step
                                    end do
                           endif

                        endif
                        !SINCRONIZA TODOS LOS LAYERS Y SOLO EL ROOT LLAMA A LA RUTINA DE ESCRITURA
#ifdef CompileWithMPI
                        call MPI_Barrier(output(ii)%item(1)%MPISubComm,ierr)
                        if (size>1) then
                           if (output(ii)%item(1)%MPISubComm /= -1) then

                              if (SGG%Observation(ii)%TimeDomain) then
                                 CALL MPI_AllReduce (Serialized%Valor, newSerialized%Valor, numberOfSerialized, REALSIZE, MPI_SUM, &
                                 &                     output(ii)%item(1)%MPISubComm, ierr)
                                 call MPI_Barrier(output(ii)%item(1)%MPISubComm,ierr)
                                 Serialized%Valor = newSerialized%Valor
                                 !
                                 CALL MPI_AllReduce (Serialized%Valor_x, newSerialized%Valor_x, numberOfSerialized, REALSIZE, MPI_SUM, &
                                 &                     output(ii)%item(1)%MPISubComm, ierr)
                                 call MPI_Barrier(output(ii)%item(1)%MPISubComm,ierr)
                                 Serialized%Valor_x = newSerialized%Valor_x
                                 !
                                 CALL MPI_AllReduce (Serialized%Valor_y, newSerialized%Valor_y, numberOfSerialized, REALSIZE, MPI_SUM, &
                                 &                     output(ii)%item(1)%MPISubComm, ierr)
                                 call MPI_Barrier(output(ii)%item(1)%MPISubComm,ierr)
                                 Serialized%Valor_y = newSerialized%Valor_y
                                 !
                                 CALL MPI_AllReduce (Serialized%Valor_z, newSerialized%Valor_z, numberOfSerialized, REALSIZE, MPI_SUM, &
                                 &                     output(ii)%item(1)%MPISubComm, ierr)
                                 call MPI_Barrier(output(ii)%item(1)%MPISubComm,ierr)
                                 Serialized%Valor_z = newSerialized%Valor_z
                                 !electric
                                 
                                 CALL MPI_AllReduce (Serialized%ValorE, newSerialized%ValorE, numberOfSerialized, REALSIZE, MPI_SUM, &
                                 &                     output(ii)%item(1)%MPISubComm, ierr)
                                 call MPI_Barrier(output(ii)%item(1)%MPISubComm,ierr)
                                 Serialized%ValorE = newSerialized%ValorE
                                 !
                                 CALL MPI_AllReduce (Serialized%Valor_Ex, newSerialized%Valor_Ex, numberOfSerialized, REALSIZE, MPI_SUM, &
                                 &                     output(ii)%item(1)%MPISubComm, ierr)
                                 call MPI_Barrier(output(ii)%item(1)%MPISubComm,ierr)
                                 Serialized%Valor_Ex = newSerialized%Valor_Ex
                                 !
                                 CALL MPI_AllReduce (Serialized%Valor_Ey, newSerialized%Valor_Ey, numberOfSerialized, REALSIZE, MPI_SUM, &
                                 &                     output(ii)%item(1)%MPISubComm, ierr)
                                 call MPI_Barrier(output(ii)%item(1)%MPISubComm,ierr)
                                 Serialized%Valor_Ey = newSerialized%Valor_Ey
                                 !
                                 CALL MPI_AllReduce (Serialized%Valor_Ez, newSerialized%Valor_Ez, numberOfSerialized, REALSIZE, MPI_SUM, &
                                 &                     output(ii)%item(1)%MPISubComm, ierr)
                                 call MPI_Barrier(output(ii)%item(1)%MPISubComm,ierr)
                                 Serialized%Valor_Ez = newSerialized%Valor_Ez
                                 !magnetic
                                 
                                 CALL MPI_AllReduce (Serialized%ValorH, newSerialized%ValorH, numberOfSerialized, REALSIZE, MPI_SUM, &
                                 &                     output(ii)%item(1)%MPISubComm, ierr)
                                 call MPI_Barrier(output(ii)%item(1)%MPISubComm,ierr)
                                 Serialized%ValorH = newSerialized%ValorH
                                 !
                                 CALL MPI_AllReduce (Serialized%Valor_Hx, newSerialized%Valor_Hx, numberOfSerialized, REALSIZE, MPI_SUM, &
                                 &                     output(ii)%item(1)%MPISubComm, ierr)
                                 call MPI_Barrier(output(ii)%item(1)%MPISubComm,ierr)
                                 Serialized%Valor_Hx = newSerialized%Valor_Hx
                                 !
                                 CALL MPI_AllReduce (Serialized%Valor_Hy, newSerialized%Valor_Hy, numberOfSerialized, REALSIZE, MPI_SUM, &
                                 &                     output(ii)%item(1)%MPISubComm, ierr)
                                 call MPI_Barrier(output(ii)%item(1)%MPISubComm,ierr)
                                 Serialized%Valor_Hy = newSerialized%Valor_Hy
                                 !
                                 CALL MPI_AllReduce (Serialized%Valor_Hz, newSerialized%Valor_Hz, numberOfSerialized, REALSIZE, MPI_SUM, &
                                 &                     output(ii)%item(1)%MPISubComm, ierr)
                                 call MPI_Barrier(output(ii)%item(1)%MPISubComm,ierr)
                                 Serialized%Valor_Hz = newSerialized%Valor_Hz
                              elseif (SGG%Observation(ii)%FreqDomain) then
!parte real
                                    Serialized%Valor_x = 0.0_RKIND
                                    newSerialized%Valor_x = 0.0_RKIND
                                    do conta=1,output(ii)%item(1)%columnas
                                          Serialized%Valor_x(1,posicionMPI+conta)=real(Serialized%ValorComplex_x(1,posicionMPI+conta))
                                    end do
                                    CALL MPI_AllReduce (Serialized%Valor_x, newSerialized%Valor_x, numberOfSerialized, REALSIZE, MPI_SUM, &
                                    &                     output(ii)%item(1)%MPISubComm, ierr)
                                    call MPI_Barrier(output(ii)%item(1)%MPISubComm,ierr)
                                    do conta=1,numberOfSerialized
                                          Serialized%ValorComplex_x(1,conta)=cmplx(newSerialized%Valor_x(1,conta),0.0_RKIND)
                                    end do
!parte imaginaria
                                    Serialized%Valor_x = 0.0_RKIND
                                    newSerialized%Valor_x = 0.0_RKIND
                                    do conta=1,output(ii)%item(1)%columnas
                                          Serialized%Valor_x(1,posicionMPI+conta)=aimag(Serialized%ValorComplex_x(1,posicionMPI+conta))
                                    end do
                                    CALL MPI_AllReduce (Serialized%Valor_x, newSerialized%Valor_x, numberOfSerialized, REALSIZE, MPI_SUM, &
                                    &                     output(ii)%item(1)%MPISubComm, ierr)
                                    call MPI_Barrier(output(ii)%item(1)%MPISubComm,ierr)
                                    do conta=1,numberOfSerialized
                                          Serialized%ValorComplex_x(1,conta)=Serialized%ValorComplex_x(1,conta)+cmplx(0.0_RKIND,newSerialized%Valor_x(1,conta))
                                    end do   
!!!!   y
!parte real
                                    Serialized%Valor_y = 0.0_RKIND
                                    newSerialized%Valor_y = 0.0_RKIND
                                    do conta=1,output(ii)%item(1)%columnas
                                          Serialized%Valor_y(1,posicionMPI+conta)=real(Serialized%ValorComplex_y(1,posicionMPI+conta))
                                 end do                            
                                    CALL MPI_AllReduce (Serialized%Valor_y, newSerialized%Valor_y, numberOfSerialized, REALSIZE, MPI_SUM, &
                                    &                     output(ii)%item(1)%MPISubComm, ierr)
                                    call MPI_Barrier(output(ii)%item(1)%MPISubComm,ierr)
                                    do conta=1,numberOfSerialized
                                          Serialized%ValorComplex_y(1,conta)=cmplx(newSerialized%Valor_y(1,conta),0.0_RKIND)
                                    end do
!parte imaginaria
                                    Serialized%Valor_y = 0.0_RKIND
                                    newSerialized%Valor_y = 0.0_RKIND
                                    do conta=1,output(ii)%item(1)%columnas
                                          Serialized%Valor_y(1,posicionMPI+conta)=aimag(Serialized%ValorComplex_y(1,posicionMPI+conta))
                                    end do
                                    CALL MPI_AllReduce (Serialized%Valor_y, newSerialized%Valor_y, numberOfSerialized, REALSIZE, MPI_SUM, &
                                    &                     output(ii)%item(1)%MPISubComm, ierr)
                                    call MPI_Barrier(output(ii)%item(1)%MPISubComm,ierr)
                                    do conta=1,numberOfSerialized
                                          Serialized%ValorComplex_y(1,conta)=Serialized%ValorComplex_y(1,conta)+cmplx(0.0_RKIND,newSerialized%Valor_y(1,conta))
                                    end do   
!!!!   z
!parte real
                                    Serialized%valor_z = 0.0_RKIND
                                    newSerialized%valor_z = 0.0_RKIND
                                    do conta=1,output(ii)%item(1)%columnas
                                          Serialized%valor_z(1,posicionMPI+conta)=real(Serialized%ValorComplex_z(1,posicionMPI+conta))
                                    end do
                                    CALL MPI_AllReduce (Serialized%valor_z, newSerialized%valor_z, numberOfSerialized, REALSIZE, MPI_SUM, &
                                    &                     output(ii)%item(1)%MPISubComm, ierr)
                                    call MPI_Barrier(output(ii)%item(1)%MPISubComm,ierr)
                                    do conta=1,numberOfSerialized
                                          Serialized%ValorComplex_z(1,conta)=cmplx(newSerialized%valor_z(1,conta),0.0_RKIND)
                                    end do
!parte imaginaria
                                    Serialized%valor_z = 0.0_RKIND
                                    newSerialized%valor_z = 0.0_RKIND
                                    do conta=1,output(ii)%item(1)%columnas
                                          Serialized%valor_z(1,posicionMPI+conta)=aimag(Serialized%ValorComplex_z(1,posicionMPI+conta))
                                    end do
                                    CALL MPI_AllReduce (Serialized%valor_z, newSerialized%valor_z, numberOfSerialized, REALSIZE, MPI_SUM, &
                                    &                     output(ii)%item(1)%MPISubComm, ierr)
                                    call MPI_Barrier(output(ii)%item(1)%MPISubComm,ierr)
                                    do conta=1,numberOfSerialized
                                          Serialized%ValorComplex_z(1,conta)=Serialized%ValorComplex_z(1,conta)+cmplx(0.0_RKIND,newSerialized%valor_z(1,conta))
                                    end do   
!ELECTRIC
                                    
!parte real
                                    Serialized%Valor_Ex = 0.0_RKIND
                                    newSerialized%Valor_Ex = 0.0_RKIND
                                    do conta=1,output(ii)%item(1)%columnas
                                          Serialized%Valor_Ex(1,posicionMPI+conta)=real(Serialized%ValorComplex_Ex(1,posicionMPI+conta))
                                    end do
                                    CALL MPI_AllReduce (Serialized%Valor_Ex, newSerialized%Valor_Ex, numberOfSerialized, REALSIZE, MPI_SUM, &
                                    &                     output(ii)%item(1)%MPISubComm, ierr)
                                    call MPI_Barrier(output(ii)%item(1)%MPISubComm,ierr)
                                    do conta=1,numberOfSerialized
                                          Serialized%ValorComplex_Ex(1,conta)=cmplx(newSerialized%Valor_Ex(1,conta),0.0_RKIND)
                                    end do
!parte imaginaria
                                    Serialized%Valor_Ex = 0.0_RKIND
                                    newSerialized%Valor_Ex = 0.0_RKIND
                                    do conta=1,output(ii)%item(1)%columnas
                                          Serialized%Valor_Ex(1,posicionMPI+conta)=aimag(Serialized%ValorComplex_Ex(1,posicionMPI+conta))
                                    end do
                                    CALL MPI_AllReduce (Serialized%Valor_Ex, newSerialized%Valor_Ex, numberOfSerialized, REALSIZE, MPI_SUM, &
                                    &                     output(ii)%item(1)%MPISubComm, ierr)
                                    call MPI_Barrier(output(ii)%item(1)%MPISubComm,ierr)
                                    do conta=1,numberOfSerialized
                                          Serialized%ValorComplex_Ex(1,conta)=Serialized%ValorComplex_Ex(1,conta)+cmplx(0.0_RKIND,newSerialized%Valor_Ex(1,conta))
                                    end do   
!!!!   y
!parte real
                                    Serialized%Valor_Ey = 0.0_RKIND
                                    newSerialized%Valor_Ey = 0.0_RKIND
                                    do conta=1,output(ii)%item(1)%columnas
                                          Serialized%Valor_Ey(1,posicionMPI+conta)=real(Serialized%ValorComplex_Ey(1,posicionMPI+conta))
                                 end do                            
                                    CALL MPI_AllReduce (Serialized%Valor_Ey, newSerialized%Valor_Ey, numberOfSerialized, REALSIZE, MPI_SUM, &
                                    &                     output(ii)%item(1)%MPISubComm, ierr)
                                    call MPI_Barrier(output(ii)%item(1)%MPISubComm,ierr)
                                    do conta=1,numberOfSerialized
                                          Serialized%ValorComplex_Ey(1,conta)=cmplx(newSerialized%Valor_Ey(1,conta),0.0_RKIND)
                                    end do
!parte imaginaria
                                    Serialized%Valor_Ey = 0.0_RKIND
                                    newSerialized%Valor_Ey = 0.0_RKIND
                                    do conta=1,output(ii)%item(1)%columnas
                                          Serialized%Valor_Ey(1,posicionMPI+conta)=aimag(Serialized%ValorComplex_Ey(1,posicionMPI+conta))
                                    end do
                                    CALL MPI_AllReduce (Serialized%Valor_Ey, newSerialized%Valor_Ey, numberOfSerialized, REALSIZE, MPI_SUM, &
                                    &                     output(ii)%item(1)%MPISubComm, ierr)
                                    call MPI_Barrier(output(ii)%item(1)%MPISubComm,ierr)
                                    do conta=1,numberOfSerialized
                                          Serialized%ValorComplex_Ey(1,conta)=Serialized%ValorComplex_Ey(1,conta)+cmplx(0.0_RKIND,newSerialized%Valor_Ey(1,conta))
                                    end do   
!!!!   z
!parte real
                                    Serialized%valor_Ez = 0.0_RKIND
                                    newSerialized%valor_Ez = 0.0_RKIND
                                    do conta=1,output(ii)%item(1)%columnas
                                          Serialized%valor_Ez(1,posicionMPI+conta)=real(Serialized%ValorComplex_Ez(1,posicionMPI+conta))
                                    end do
                                    CALL MPI_AllReduce (Serialized%valor_Ez, newSerialized%valor_Ez, numberOfSerialized, REALSIZE, MPI_SUM, &
                                    &                     output(ii)%item(1)%MPISubComm, ierr)
                                    call MPI_Barrier(output(ii)%item(1)%MPISubComm,ierr)
                                    do conta=1,numberOfSerialized
                                          Serialized%ValorComplex_Ez(1,conta)=cmplx(newSerialized%valor_Ez(1,conta),0.0_RKIND)
                                    end do
!parte imaginaria
                                    Serialized%valor_Ez = 0.0_RKIND
                                    newSerialized%valor_Ez = 0.0_RKIND
                                    do conta=1,output(ii)%item(1)%columnas
                                          Serialized%valor_Ez(1,posicionMPI+conta)=aimag(Serialized%ValorComplex_Ez(1,posicionMPI+conta))
                                    end do
                                    CALL MPI_AllReduce (Serialized%valor_Ez, newSerialized%valor_Ez, numberOfSerialized, REALSIZE, MPI_SUM, &
                                    &                     output(ii)%item(1)%MPISubComm, ierr)
                                    call MPI_Barrier(output(ii)%item(1)%MPISubComm,ierr)
                                    do conta=1,numberOfSerialized
                                          Serialized%ValorComplex_Ez(1,conta)=Serialized%ValorComplex_Ez(1,conta)+cmplx(0.0_RKIND,newSerialized%valor_Ez(1,conta))
                                    end do   
!MAGNETIC
                                    
!parte real
                                    Serialized%Valor_Hx = 0.0_RKIND
                                    newSerialized%Valor_Hx = 0.0_RKIND
                                    do conta=1,output(ii)%item(1)%columnas
                                          Serialized%Valor_Hx(1,posicionMPI+conta)=real(Serialized%ValorComplex_Hx(1,posicionMPI+conta))
                                    end do
                                    CALL MPI_AllReduce (Serialized%Valor_Hx, newSerialized%Valor_Hx, numberOfSerialized, REALSIZE, MPI_SUM, &
                                    &                     output(ii)%item(1)%MPISubComm, ierr)
                                    call MPI_Barrier(output(ii)%item(1)%MPISubComm,ierr)
                                    do conta=1,numberOfSerialized
                                          Serialized%ValorComplex_Hx(1,conta)=cmplx(newSerialized%Valor_Hx(1,conta),0.0_RKIND)
                                    end do
!parte imaginaria
                                    Serialized%Valor_Hx = 0.0_RKIND
                                    newSerialized%Valor_Hx = 0.0_RKIND
                                    do conta=1,output(ii)%item(1)%columnas
                                          Serialized%Valor_Hx(1,posicionMPI+conta)=aimag(Serialized%ValorComplex_Hx(1,posicionMPI+conta))
                                    end do
                                    CALL MPI_AllReduce (Serialized%Valor_Hx, newSerialized%Valor_Hx, numberOfSerialized, REALSIZE, MPI_SUM, &
                                    &                     output(ii)%item(1)%MPISubComm, ierr)
                                    call MPI_Barrier(output(ii)%item(1)%MPISubComm,ierr)
                                    do conta=1,numberOfSerialized
                                          Serialized%ValorComplex_Hx(1,conta)=Serialized%ValorComplex_Hx(1,conta)+cmplx(0.0_RKIND,newSerialized%Valor_Hx(1,conta))
                                    end do   
!!!!   y
!parte real
                                    Serialized%Valor_Hy = 0.0_RKIND
                                    newSerialized%Valor_Hy = 0.0_RKIND
                                    do conta=1,output(ii)%item(1)%columnas
                                          Serialized%Valor_Hy(1,posicionMPI+conta)=real(Serialized%ValorComplex_Hy(1,posicionMPI+conta))
                                 end do                            
                                    CALL MPI_AllReduce (Serialized%Valor_Hy, newSerialized%Valor_Hy, numberOfSerialized, REALSIZE, MPI_SUM, &
                                    &                     output(ii)%item(1)%MPISubComm, ierr)
                                    call MPI_Barrier(output(ii)%item(1)%MPISubComm,ierr)
                                    do conta=1,numberOfSerialized
                                          Serialized%ValorComplex_Hy(1,conta)=cmplx(newSerialized%Valor_Hy(1,conta),0.0_RKIND)
                                    end do
!parte imaginaria
                                    Serialized%Valor_Hy = 0.0_RKIND
                                    newSerialized%Valor_Hy = 0.0_RKIND
                                    do conta=1,output(ii)%item(1)%columnas
                                          Serialized%Valor_Hy(1,posicionMPI+conta)=aimag(Serialized%ValorComplex_Hy(1,posicionMPI+conta))
                                    end do
                                    CALL MPI_AllReduce (Serialized%Valor_Hy, newSerialized%Valor_Hy, numberOfSerialized, REALSIZE, MPI_SUM, &
                                    &                     output(ii)%item(1)%MPISubComm, ierr)
                                    call MPI_Barrier(output(ii)%item(1)%MPISubComm,ierr)
                                    do conta=1,numberOfSerialized
                                          Serialized%ValorComplex_Hy(1,conta)=Serialized%ValorComplex_Hy(1,conta)+cmplx(0.0_RKIND,newSerialized%Valor_Hy(1,conta))
                                    end do   
!!!!   z
!parte real
                                    Serialized%valor_Hz = 0.0_RKIND
                                    newSerialized%valor_Hz = 0.0_RKIND
                                    do conta=1,output(ii)%item(1)%columnas
                                          Serialized%valor_Hz(1,posicionMPI+conta)=real(Serialized%ValorComplex_Hz(1,posicionMPI+conta))
                                    end do
                                    CALL MPI_AllReduce (Serialized%valor_Hz, newSerialized%valor_Hz, numberOfSerialized, REALSIZE, MPI_SUM, &
                                    &                     output(ii)%item(1)%MPISubComm, ierr)
                                    call MPI_Barrier(output(ii)%item(1)%MPISubComm,ierr)
                                    do conta=1,numberOfSerialized
                                          Serialized%ValorComplex_Hz(1,conta)=cmplx(newSerialized%valor_Hz(1,conta),0.0_RKIND)
                                    end do
!parte imaginaria
                                    Serialized%valor_Hz = 0.0_RKIND
                                    newSerialized%valor_Hz = 0.0_RKIND
                                    do conta=1,output(ii)%item(1)%columnas
                                          Serialized%valor_Hz(1,posicionMPI+conta)=aimag(Serialized%ValorComplex_Hz(1,posicionMPI+conta))
                                    end do
                                    CALL MPI_AllReduce (Serialized%valor_Hz, newSerialized%valor_Hz, numberOfSerialized, REALSIZE, MPI_SUM, &
                                    &                     output(ii)%item(1)%MPISubComm, ierr)
                                    call MPI_Barrier(output(ii)%item(1)%MPISubComm,ierr)
                                    do conta=1,numberOfSerialized
                                          Serialized%ValorComplex_Hz(1,conta)=Serialized%ValorComplex_Hz(1,conta)+cmplx(0.0_RKIND,newSerialized%valor_Hz(1,conta))
                                    end do   
                                 
                              endif

                           endif
                        endif
                        IF (layoutnumber == output(ii)%item(1)%MPIRoot) THEN
#else
                        IF (layoutnumber == 0) THEN
#endif
                           !
                           time=att(indi)
                           WRITE (charc, '(i10)') indi
                           fichero=trim(adjustl(filename))//'_'//trim (adjustl(charc))//'.vtk'


                           if ((.not.dontwritevtk).or.(sgg%observation(ii)%P(1)%What==mapvtk)) then !el mapvtk lo procesa siempre
                                 write(dubuf,'(a,i9,a,i9)')  ' ----> file '//trim(adjustl(fichero))//' ',indi,'/',finalstep
                                 call print11(layoutnumber,dubuf)

                                 iroot1 = index (trim(adjustl(fichero)), '.vtk',.true.)     
                                 iroot2 = index (trim(adjustl(fichero(1:iroot1))),'_',.true.)
                                 iroot2=iroot2-1
                                 if (indi==1) then 
                                    block
                                       logical :: dir_e 
                                       !inquire(DIRECTORY= trim(adjustl(fichero(1:iroot2))), exist=dir_e)
                                       dir_e=.false. !intento crearlo por defecto 0624: ya dara un warning el sistema 
                                       if ( dir_e ) then         
                                          continue
                                          ! write(*,*) "dir exists! "//trim(probeName)
                                       else
                                          ! workaround: it calls an extern program...  
                                          CALL SYSTEM('mkdir ' // trim(adjustl(fichero(1:iroot2))))  
                                       end if     
                                    end block
                                 endif 
                                 if (sgg%observation(ii)%P(1)%What==mapvtk) then
                                       fichero_input=fichero(1:iroot1-1)//'.vtk'  
                                       i_sub_time=-30 !cualquier cosa
                                       total_sub_times=-12 !cualquier cosa
                                       CALL write_VTKfile(sgg,fichero_input,iroot2, Serialized,  numberOfSerialized,Nodes,Numnodes,Elems,NumEdges,NumQuads,time, &
                                                            i_sub_time,total_sub_times,freqDomain,sgg%observation(ii)%P(1)%What,sggMtag,'vt')  
                                 else
                                    if (freqDomain) then
                                       total_sub_times=time_phases_param
                                       do i_sub_time=0,total_sub_times
                                             write (char_i_sub_time,'(i3)') i_sub_time
                                             fichero_input=fichero(1:iroot1-1)//'_n_'//trim(adjustl(char_i_sub_time))//'_current.vtk' 
                                             CALL write_VTKfile(sgg,fichero_input,iroot2, Serialized,  numberOfSerialized,Nodes,Numnodes,Elems,NumEdges,NumQuads,time, &
                                                               i_sub_time,total_sub_times,freqDomain,sgg%observation(ii)%P(1)%What,sggMtag,'cu')
                                             write(dubuf,'(a,i9,a,i9)')  trim(adjustl(whoamishort))//' -------> Dumped frequency phase file '//trim(adjustl(fichero_input))//', ',i_sub_time,'/',total_sub_times
                                             call print11(layoutnumber,dubuf,.true.)        
                                             !electric
                                       
                                             fichero_input=fichero(1:iroot1-1)//'_n_'//trim(adjustl(char_i_sub_time))//'_efield.vtk' 
                                             CALL write_VTKfile(sgg,fichero_input,iroot2, Serialized,  numberOfSerialized,Nodes,Numnodes,Elems,NumEdges,NumQuads,time, &
                                                               i_sub_time,total_sub_times,freqDomain,sgg%observation(ii)%P(1)%What,sggMtag,'ef')
                                             write(dubuf,'(a,i9,a,i9)')  trim(adjustl(whoamishort))//' -------> Dumped frequency phase file '//trim(adjustl(fichero_input))//', ',i_sub_time,'/',total_sub_times
                                             call print11(layoutnumber,dubuf,.true.) 
                                             !      
                                             !      magnetic
                                       
                                             fichero_input=fichero(1:iroot1-1)//'_n_'//trim(adjustl(char_i_sub_time))//'_hfield.vtk' 
                                             CALL write_VTKfile(sgg,fichero_input,iroot2, Serialized,  numberOfSerialized,Nodes,Numnodes,Elems,NumEdges,NumQuads,time, &
                                                               i_sub_time,total_sub_times,freqDomain,sgg%observation(ii)%P(1)%What,sggMtag,'hf')
                                             write(dubuf,'(a,i9,a,i9)')  trim(adjustl(whoamishort))//' -------> Dumped frequency phase file '//trim(adjustl(fichero_input))//', ',i_sub_time,'/',total_sub_times
                                             call print11(layoutnumber,dubuf,.true.) 
                                             !      
                                       end do
                                    else
                                       fichero_input=fichero(1:iroot1-1)//'_current.vtk'  
                                       i_sub_time=-30 !cualquier cosa
                                       total_sub_times=-12 !cualquier cosa
                                       CALL write_VTKfile(sgg,fichero_input,iroot2, Serialized,  numberOfSerialized,Nodes,Numnodes,Elems,NumEdges,NumQuads,time, &
                                                            i_sub_time,total_sub_times,freqDomain,sgg%observation(ii)%P(1)%What,sggMtag,'cu')  
                                       !electric
                                    
                                       fichero_input=fichero(1:iroot1-1)//'_efield.vtk'  
                                       i_sub_time=-30 !cualquier cosa
                                       total_sub_times=-12 !cualquier cosa
                                       CALL write_VTKfile(sgg,fichero_input,iroot2, Serialized,  numberOfSerialized,Nodes,Numnodes,Elems,NumEdges,NumQuads,time, &
                                                            i_sub_time,total_sub_times,freqDomain,sgg%observation(ii)%P(1)%What,sggMtag,'ef')
                                    
                                    
                                       !magnetic
                                    
                                       fichero_input=fichero(1:iroot1-1)//'_hfield.vtk'  
                                       i_sub_time=-30 !cualquier cosa
                                       total_sub_times=-12 !cualquier cosa
                                       CALL write_VTKfile(sgg,fichero_input,iroot2, Serialized,  numberOfSerialized,Nodes,Numnodes,Elems,NumEdges,NumQuads,time, &
                                                            i_sub_time,total_sub_times,freqDomain,sgg%observation(ii)%P(1)%What,sggMtag,'hf')
                                    
                                    
                                    endif        
                                    
                                    !!! CALL print11 (layoutnumber, trim(adjustl(whoami))////' Written into file '//trim(adjustl(fichero)), .TRUE.) !enforces print
                                 endif !DEL VTK
                           else
                                 write(dubuf,'(a,i9,a,i9)')  trim(adjustl(whoamishort))//' Requesting not to dump .vtk ----> file '//trim(adjustl(fichero))//' ',indi,'/',finalstep
                                 call print11(layoutnumber,dubuf,.true.)
                           endif
                                 
                        ENDIF

                        !
                     END DO bucleindi
                     CLOSE (output(ii)%item(1)%UNIT)
                     !
                     if (SGG%Observation(ii)%TimeDomain) then 
                        DEALLOCATE (Serialized%Valor)
                        DEALLOCATE (Serialized%Valor_x)
                        DEALLOCATE (Serialized%Valor_y)
                        DEALLOCATE (Serialized%Valor_z)
                        DEALLOCATE (Serialized%ValorE)
                        DEALLOCATE (Serialized%Valor_Ex)
                        DEALLOCATE (Serialized%Valor_Ey)
                        DEALLOCATE (Serialized%Valor_Ez)
                        DEALLOCATE (Serialized%ValorH)
                        DEALLOCATE (Serialized%Valor_Hx)
                        DEALLOCATE (Serialized%Valor_Hy)
                        DEALLOCATE (Serialized%Valor_Hz)
                     elseif (SGG%Observation(ii)%FreqDomain) then    
                        DEALLOCATE (Serialized%Valor)
                        DEALLOCATE (Serialized%Valor_x)
                        DEALLOCATE (Serialized%Valor_y)
                        DEALLOCATE (Serialized%Valor_z)    
                        DEALLOCATE (Serialized%ValorComplex_x)
                        DEALLOCATE (Serialized%ValorComplex_y)
                        DEALLOCATE (Serialized%ValorComplex_z)
                        
                        DEALLOCATE (Serialized%ValorE)
                        DEALLOCATE (Serialized%Valor_Ex)
                        DEALLOCATE (Serialized%Valor_Ey)
                        DEALLOCATE (Serialized%Valor_Ez)    
                        DEALLOCATE (Serialized%ValorComplex_Ex)
                        DEALLOCATE (Serialized%ValorComplex_Ey)
                        DEALLOCATE (Serialized%ValorComplex_Ez)
                        
                        DEALLOCATE (Serialized%ValorH)
                        DEALLOCATE (Serialized%Valor_Hx)
                        DEALLOCATE (Serialized%Valor_Hy)
                        DEALLOCATE (Serialized%Valor_Hz)    
                        DEALLOCATE (Serialized%ValorComplex_Hx)
                        DEALLOCATE (Serialized%ValorComplex_Hy)
                        DEALLOCATE (Serialized%ValorComplex_Hz)
                     endif
                     deallocate (Serialized%eI)
                     deallocate (Serialized%eJ)
                     deallocate (Serialized%eK)
                     deallocate (Serialized%currentType)
                     deallocate (Serialized%sggMtag)
#ifdef CompileWithMPI
                     if (size>1) then
                        if (SGG%Observation(ii)%TimeDomain) then  
                           DEALLOCATE (NewSerialized%Valor)
                           DEALLOCATE (NewSerialized%Valor_x)
                           DEALLOCATE (NewSerialized%Valor_y)
                           DEALLOCATE (NewSerialized%Valor_z)
                           DEALLOCATE (NewSerialized%ValorE)
                           DEALLOCATE (NewSerialized%Valor_Ex)
                           DEALLOCATE (NewSerialized%Valor_Ey)
                           DEALLOCATE (NewSerialized%Valor_Ez)
                           DEALLOCATE (NewSerialized%ValorH)
                           DEALLOCATE (NewSerialized%Valor_Hx)
                           DEALLOCATE (NewSerialized%Valor_Hy)
                           DEALLOCATE (NewSerialized%Valor_Hz)
                        elseif (SGG%Observation(ii)%FreqDomain) then  
                           DEALLOCATE (NewSerialized%Valor) !auxiliar
                           DEALLOCATE (NewSerialized%Valor_x) !auxiliar
                           DEALLOCATE (NewSerialized%Valor_y) !auxiliar
                           DEALLOCATE (NewSerialized%Valor_z) !auxiliar
   !!                          DEALLOCATE (NewSerialized%ValorComplex)
                           
                           DEALLOCATE (NewSerialized%ValorE) !auxiliar
                           DEALLOCATE (NewSerialized%Valor_Ex) !auxiliar
                           DEALLOCATE (NewSerialized%Valor_Ey) !auxiliar
                           DEALLOCATE (NewSerialized%Valor_Ez) !auxiliar
   !!                          DEALLOCATE (NewSerialized%ValorComplexE)
                           DEALLOCATE (NewSerialized%ValorH) !auxiliar
                           DEALLOCATE (NewSerialized%Valor_Hx) !auxiliar
                           DEALLOCATE (NewSerialized%Valor_Hy) !auxiliar
                           DEALLOCATE (NewSerialized%Valor_Hz) !auxiliar
   !!                          DEALLOCATE (NewSerialized%ValorComplexH)
                        endif
                        deallocate (newSerialized%eI)
                        deallocate (newSerialized%eJ)
                        deallocate (newSerialized%eK)
                        deallocate (newSerialized%currentType)
                        deallocate (newSerialized%sggMtag)
                     ENDIF
#endif

                     !deallocatea
#ifdef CompileWithMPI
                     IF (layoutnumber == output(ii)%item(1)%MPIRoot) THEN
#else
                     IF (layoutnumber == 0) THEN
#endif
                        if (numberOfSerialized/=0) deallocate (Nodes,Elems)
                     endif

#ifdef CompileWithMPI
                     if (size>1) then
                        DEALLOCATE (newSizeofvalores,newPosiMPI)
                     endif
#endif
                     DEALLOCATE (SIZEOFVALORES,PosiMPI)
                     DEALLOCATE (ATT)
#ifdef CompileWithMPI
                     if (size>1) then
                        if (output(ii)%item(1)%MPISubComm /= -1) then
                           call MPI_Barrier(output(ii)%item(1)%MPISubComm,ierr)
                        endif
                        !!! CALL print11 (layoutnumber, trim(adjustl(whoami))////' End processing file '//trim(adjustl(filename)), .TRUE.) !enforces print
                     endif
#endif
                  else !del lexis
                     buff='NOT PROCESSING: Ignoring: Inexistent or void file '//trim(adjustl(output(ii)%item(1)%path))
                     CALL print11(layoutnumber, buff,.true.)
                  endif !del lexis


               somethingdone=.true.

               ENDIF !DEL WHAT
            ENDIF
         ENDIF

      END DO  barridoprobes !barrido puntos de observacion



      RETURN
   END SUBROUTINE createVTK

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   SUBROUTINE createVTKOnTheFly (layoutnumber, size, sgg,vtkindex,somethingdone,mpidir,tagtype,sggMtag,dontwritevtk)
   
      type (SGGFDTDINFO), intent(IN)    :: sgg
      INTEGER (KIND=IKINDMTAG), intent(in) ::  sggMtag  (sgg%Alloc(iHx)%XI:sgg%Alloc(iHx)%XE, sgg%Alloc(iHy)%YI:sgg%Alloc(iHy)%YE, sgg%Alloc(iHz)%ZI:sgg%Alloc(iHz)%ZE)
   
      type (tagtype_t) :: tagtype
   
      integer (KIND=4) :: mpidir
      logical :: vtkindex,somethingdone

      INTEGER (KIND=4), INTENT (IN) :: layoutnumber, size
      TYPE (output_t), POINTER, DIMENSION (:) :: output
      INTEGER (KIND=4) :: ii
      logical :: lexis,dontwritevtk
      character(len=BUFSIZE) :: buff
      character (LEN=BUFSIZE) :: path


      !
      output => GetOutput ()!get the output private info from observation
      !

      DO ii = 1, sgg%NumberRequest
         !sondas Volumic traducelas a xdfm
         IF (sgg%observation(ii)%Volumic) then
            if (sgg%observation(ii)%nP == 1) then

               if ((sgg%observation(ii)%P(1)%What == iCur).or.(sgg%observation(ii)%P(1)%What == iCurX).or. &
               (sgg%observation(ii)%P(1)%What == iCurY).or.(sgg%observation(ii)%P(1)%What == iCurZ).or. &
               (sgg%observation(ii)%P(1)%What == mapvtk)) THEN !solo corrientes volumicas
                  !
                  INQUIRE (FILE=trim(adjustl(output(ii)%item(1)%path)), EXIST=lexis)
                  if (.not.lexis) then
                     buff='NOT PROCESSING: Inexistent file '//trim(adjustl(output(ii)%item(1)%path))
                     CALL print11(layoutnumber, buff,.true.)
                     return
                  ELSE
                     close (output(ii)%item(1)%unit)
                  ENDIF !DEL LEXIS
               ENDIF
            ENDIF
         ENDIF

      END DO  !barrido puntos de observacion
      call createVTK (layoutnumber, size, sgg,vtkindex,somethingdone,mpidir,tagtype,sggMtag,dontwritevtk)
      DO ii = 1, sgg%NumberRequest
         !sondas Volumic traducelas a xdfm
         IF (sgg%observation(ii)%Volumic) then
            if (sgg%observation(ii)%nP == 1) then
               if ((sgg%observation(ii)%P(1)%What == iCur).or.(sgg%observation(ii)%P(1)%What == iCurX).or.(sgg%observation(ii)%P(1)%What == iCurY).or.(sgg%observation(ii)%P(1)%What == iCurZ).or. &
               (sgg%observation(ii)%P(1)%What == mapvtk)) THEN !solo corrientes volumicas
                  !
                  INQUIRE (FILE=trim(adjustl(output(ii)%item(1)%path)), EXIST=lexis)
                  if (.not.lexis) then
                     buff='NOT PROCESSING: Inexistent file '//trim(adjustl(output(ii)%item(1)%path))
                     CALL print11(layoutnumber, buff,.true.)
                     return
                  ELSE
                     open (output(ii)%item(1)%unit,file=trim(adjustl(output(ii)%item(1)%path)),FORM='unformatted',position='append')
                  ENDIF !DEL LEXIS
               ENDIF
            ENDIF
         ENDIF

      END DO  !barrido puntos de observacion

      RETURN
   END SUBROUTINE createVTKOnTheFly


   !!!!!!!

   SUBROUTINE write_VTKfile(sgg,fichero,iroot2, Serialized,  numberOfSerialized,Nodes,Numnodes,Elems,NumEdges,NumQuads,time,  &
                              i_sub_time,total_sub_times,FreqDomain,what,sggMtag,que_saco)
   
      type (SGGFDTDINFO), intent(IN)   ::  sgg
      INTEGER (KIND=IKINDMTAG), intent(in) ::  sggMtag  (sgg%Alloc(iHx)%XI:sgg%Alloc(iHx)%XE, sgg%Alloc(iHy)%YI:sgg%Alloc(iHy)%YE, sgg%Alloc(iHz)%ZI:sgg%Alloc(iHz)%ZE)
      CHARACTER (LEN=BUFSIZE), intent(in) :: fichero

      type (Serialized_t), intent(in)::  Serialized            
      INTEGER (kind=4), intent(in):: numberOfSerialized,numNodes,numEdges,NumQuads,iroot2,i_sub_time,total_sub_times    
      real (kind=RKIND), intent(in) :: time
      real (kind=RKIND) :: phase_x,phase_y,phase_z,raa,rbb,rcc 
      real (kind=RKIND) :: phase_Ex,phase_Ey,phase_Ez
      real (kind=RKIND) :: phase_Hx,phase_Hy,phase_Hz
      LOGICAL, intent(in) :: FREQDOMAIN
      integer (kind=4), intent(in):: what
      INTEGER (kind=4) :: conta,myunit
      CHARACTER (LEN=BUFSIZE) :: buff,buff2 ! File name
      real (kind= RKIND), allocatable, dimension(:,:) :: Nodes
      integer (kind=4), allocatable, dimension(:,:) ::  Elems
      character*2, intent(in) :: que_saco
      
      !!!!!!!
      !if (what==mapvtk) THEN
      !call fillinparaviewstate
      !open(newunit=myunit,file=trim(adjustl(fichero))//'.pvsm',form='formatted')
      !    write (myunit,'(a)') 'Generador del .pvsm'
      !close(myunit)
      !endif

      !!!! 

      open(newunit=myunit,file=trim(adjustl(fichero(1:iroot2)))//'/'//trim(adjustl(fichero)),form='formatted')
      close(myunit,status='delete')
      open(newunit=myunit,file=trim(adjustl(fichero(1:iroot2)))//'/'//trim(adjustl(fichero)),form='formatted')
      write(myunit,'(a)') '# vtk DataFile Version 1.0'
      !a modo de ayuda saco en el fichero MAP el tipo de material en la segunda linea como manda el standard vtk
      if (what==mapvtk) THEN
         write(myunit,'(a)') 'PEC=0, already_YEEadvanced_byconformal=5, NOTOUCHNOUSE=6, WIRE=7, WIRE-COLISION=8, COMPO=3, DISPER=1, DIEL=2, SLOT=4, CONF=5/6, OTHER=-1 (ADD +0.5 for borders)'
      else                                      
         if (.not.Freqdomain) then
            write(myunit,'(a,e21.12e3)') 'Time= ',time
         else
            write(myunit,'(a,e21.12e3)') 'Frequency= ',time   
         endif
      endif
      write(myunit,'(a)') 'ASCII'
      write(myunit,'(a)') ' '
      write(myunit,'(a)') 'DATASET UNSTRUCTURED_GRID'
      write(myunit,'(a)') 'FIELD FieldData 1'
      write(myunit,'(a)') 'TIME 1 1 double'
      write(myunit,'(e21.12e3)')  time
      write (buff,'(a,i9,a)') 'POINTS ',numNodes+1,' float'
      write(myunit,'(a)') trim(adjustl(buff))
      do conta=0,NumNodes
         write (buff,'(3e21.12e3)')  Nodes(conta,1), Nodes(conta,2), Nodes(conta,3)
         write(myunit,'(a)') trim(adjustl(buff))
      end do
      write(myunit,'(a)') ' '
      write (buff,'(a,2i9)') 'CELLS ',(NumEdges+1)+(NumQuads+1),3*(NumEdges+1)+5*(NumQuads+1)
      write(myunit,'(a)') trim(adjustl(buff))
      do conta=1,numberOfSerialized
         if (Elems(conta,3)==-1) then !es un edge
            write (myunit,'(i2,2i9)')  2,Elems(conta,1), Elems(conta,2)
         else
            write (myunit,'(i2,4i9)')  4,Elems(conta,1), Elems(conta,2), Elems(conta,3), Elems(conta,4)
         endif
      end do
      write(myunit,'(a)') ' '
      write (buff,'(a,i9)') 'CELL_TYPES ',(NumEdges+1)+(NumQuads+1)
      write(myunit,'(a)') trim(adjustl(buff))
      do conta=1,numberOfSerialized
         if (Elems(conta,3)==-1) then !es un edge
            write (myunit,'(i2)')  3
         else
            write (myunit,'(i2)')  9
         endif
      end do
      write(myunit,'(a)') ' '
      write (buff,'(a,i9)') 'CELL_DATA ',numberOfSerialized
      write(myunit,'(a)') trim(adjustl(buff))
      write (buff2,'(e21.12e3)') time
      if ((what==mapvtk).AND.(que_saco=='vt')) THEN
            write (buff,'(a)') 'SCALARS mediatype float 1'
      else                   
         select case(que_saco)           
         case('cu')
         if (.not.Freqdomain) then
!probando2023 vectores        write (buff,'(a)') 'SCALARS current_t float 1'
               write (buff,'(a)') 'SCALARS current_t float 3'
         else
            write (buff,'(a)') 'SCALARS current_f float 3'
         endif     
         case('ef')    
         if (.not.Freqdomain) then
!probando2023 vectores        write (buff,'(a)') 'SCALARS current_t float 1'
               write (buff,'(a)') 'SCALARS efield_t float 3'
         else
            write (buff,'(a)') 'SCALARS efield_f float 3'
         endif
         case('hf')    
         if (.not.Freqdomain) then
!probando2023 vectores        write (buff,'(a)') 'SCALARS current_t float 1'
               write (buff,'(a)') 'SCALARS hfield_t float 3'
         else
            write (buff,'(a)') 'SCALARS hfield_f float 3'
         endif
         end select
      endif
      write(myunit,'(a)') trim(adjustl(buff))
      write(myunit,'(a)') 'LOOKUP_TABLE default'

      if (.not.Freqdomain) then
         do conta=1,numberOfSerialized
!  Vectorial 0124
               if (what==mapvtk) THEN      
                  write (myunit,'(1e21.12e3)')  Serialized%valor(1,conta)      !sin vectores
               else             
               select case(que_saco)          
               case('cu')          
                  raa=Serialized%valor_x(1,conta)
                  rbb=Serialized%valor_y(1,conta)
                  rcc=Serialized%valor_z(1,conta)   
               case('ef') 
                  raa=Serialized%valor_Ex(1,conta)
                  rbb=Serialized%valor_Ey(1,conta)
                  rcc=Serialized%valor_Ez(1,conta)   
               case('hf')  
                  raa=Serialized%valor_Hx(1,conta)
                  rbb=Serialized%valor_Hy(1,conta)
                  rcc=Serialized%valor_Hz(1,conta)   
               end select
               if (raa>1.e37)  raa=1.e37; if (raa<-1.e37) raa=-1.e37; if (abs(raa)<1e-37 ) raa=0.  
               if (rbb>1.e37)  rbb=1.e37; if (rbb<-1.e37) rbb=-1.e37; if (abs(rbb)<1e-37 ) rbb=0.
               if (rcc>1.e37)  rcc=1.e37; if (rcc<-1.e37) rcc=-1.e37; if (abs(rcc)<1e-37 ) rcc=0.
               write (myunit,'(3e21.12e3)')  raa,rbb,rcc
               endif
         end do
      else
         do conta=1,numberOfSerialized             
            select case(que_saco)          
            case('cu')                                      
               phase_x=atan2(AIMAG(Serialized%valorComplex_x(1,conta)),REAL(Serialized%valorComplex_x(1,conta)))
               phase_y=atan2(AIMAG(Serialized%valorComplex_y(1,conta)),REAL(Serialized%valorComplex_y(1,conta)))
               phase_z=atan2(AIMAG(Serialized%valorComplex_z(1,conta)),REAL(Serialized%valorComplex_z(1,conta)))
               raa=abs(Serialized%valorComplex_x(1,conta))*cos(real(i_sub_time)*2.*pi/real(total_sub_times)+phase_x) 
               rbb=abs(Serialized%valorComplex_y(1,conta))*cos(real(i_sub_time)*2.*pi/real(total_sub_times)+phase_y)
               rcc=abs(Serialized%valorComplex_z(1,conta))*cos(real(i_sub_time)*2.*pi/real(total_sub_times)+phase_z)          
               if (raa>1.e37)  raa=1.e37; if (raa<-1.e37) raa=-1.e37; if (abs(raa)<1e-37 ) raa=0.      !bug 1e-40 unsupported in paraview 
               if (rbb>1.e37)  rbb=1.e37; if (rbb<-1.e37) rbb=-1.e37; if (abs(rbb)<1e-37 ) rbb=0.
               if (rcc>1.e37)  rcc=1.e37; if (rcc<-1.e37) rcc=-1.e37; if (abs(rcc)<1e-37 ) rcc=0.
               write (myunit,'(3e21.12e3)')  raa,rbb,rcc      
               case('ef')                                  
               phase_Ex=atan2(AIMAG(Serialized%valorComplex_Ex(1,conta)),REAL(Serialized%valorComplex_Ex(1,conta)))
               phase_Ey=atan2(AIMAG(Serialized%valorComplex_Ey(1,conta)),REAL(Serialized%valorComplex_Ey(1,conta)))
               phase_Ez=atan2(AIMAG(Serialized%valorComplex_Ez(1,conta)),REAL(Serialized%valorComplex_Ez(1,conta)))
               raa=abs(Serialized%valorComplex_Ex(1,conta))*cos(real(i_sub_time)*2.*pi/real(total_sub_times)+phase_Ex) 
               rbb=abs(Serialized%valorComplex_Ey(1,conta))*cos(real(i_sub_time)*2.*pi/real(total_sub_times)+phase_Ey)
               rcc=abs(Serialized%valorComplex_Ez(1,conta))*cos(real(i_sub_time)*2.*pi/real(total_sub_times)+phase_Ez)          
               if (raa>1.e37)  raa=1.e37; if (raa<-1.e37) raa=-1.e37; if (abs(raa)<1e-37 ) raa=0.      !bug 1e-40 unsupported in paraview 
               if (rbb>1.e37)  rbb=1.e37; if (rbb<-1.e37) rbb=-1.e37; if (abs(rbb)<1e-37 ) rbb=0.
               if (rcc>1.e37)  rcc=1.e37; if (rcc<-1.e37) rcc=-1.e37; if (abs(rcc)<1e-37 ) rcc=0.
               write (myunit,'(3e21.12e3)')  raa,rbb,rcc      
               case('hf')                                  
               phase_Ex=atan2(AIMAG(Serialized%valorComplex_Ex(1,conta)),REAL(Serialized%valorComplex_Ex(1,conta)))
               phase_Ey=atan2(AIMAG(Serialized%valorComplex_Ey(1,conta)),REAL(Serialized%valorComplex_Ey(1,conta)))
               phase_Ez=atan2(AIMAG(Serialized%valorComplex_Ez(1,conta)),REAL(Serialized%valorComplex_Ez(1,conta)))
               raa=abs(Serialized%valorComplex_Ex(1,conta))*cos(real(i_sub_time)*2.*pi/real(total_sub_times)+phase_Ex) 
               rbb=abs(Serialized%valorComplex_Ey(1,conta))*cos(real(i_sub_time)*2.*pi/real(total_sub_times)+phase_Ey)
               rcc=abs(Serialized%valorComplex_Ez(1,conta))*cos(real(i_sub_time)*2.*pi/real(total_sub_times)+phase_Ez)          
               if (raa>1.e37)  raa=1.e37; if (raa<-1.e37) raa=-1.e37; if (abs(raa)<1e-37 ) raa=0.      !bug 1e-40 unsupported in paraview 
               if (rbb>1.e37)  rbb=1.e37; if (rbb<-1.e37) rbb=-1.e37; if (abs(rbb)<1e-37 ) rbb=0.
               if (rcc>1.e37)  rcc=1.e37; if (rcc<-1.e37) rcc=-1.e37; if (abs(rcc)<1e-37 ) rcc=0.
               write (myunit,'(3e21.12e3)')  raa,rbb,rcc      
               end select
                                             
         end do
      endif

      write(myunit,'(a)') ' '

      !!!info del tag 240220
      write (buff,'(a,i9)') 'SCALARS tagnumber float 1'
      write (myunit,'(a)') trim(adjustl(buff))
      write (buff,'(a)') 'LOOKUP_TABLE default'
      write (myunit,'(a)') trim(adjustl(buff))
      do conta=1,numberOfSerialized
            !!!escribo por exceso tags en sitios donde no hay realmente ese medio incluyendo en quads solo porque uno de sus lados tiene ese tag. arreglar algun dia aunque no es critiico porque los tags solo serviran para filtran luego visualizaciones !240220
               !if (Elems(conta,3)==-1) then !es un edge
               !   write (myunit,'(i4)')  sggMtag(Serialized%eI(conta),Serialized%eJ(conta),Serialized%eK(conta)) 
               !else !es un quad y hay que verlo mejor
               !    if ( ((sggMtag(Serialized%eI(conta),Serialized%eJ(conta),Serialized%eK(conta))==sggMtag(Serialized%eI(conta)+1,Serialized%eJ(conta)  ,Serialized%eK(conta)  )).and. &
               !          (sggMtag(Serialized%eI(conta),Serialized%eJ(conta),Serialized%eK(conta))==sggMtag(Serialized%eI(conta)  ,Serialized%eJ(conta)+1,Serialized%eK(conta)  )).and. &
               !          (sggMtag(Serialized%eI(conta),Serialized%eJ(conta),Serialized%eK(conta))==sggMtag(Serialized%eI(conta)+1,Serialized%eJ(conta)+1,Serialized%eK(conta)  )) ).or. &
               !        !
               !         ((sggMtag(Serialized%eI(conta),Serialized%eJ(conta),Serialized%eK(conta))==sggMtag(Serialized%eI(conta)+1,Serialized%eJ(conta)  ,Serialized%eK(conta)  )).and. &
               !          (sggMtag(Serialized%eI(conta),Serialized%eJ(conta),Serialized%eK(conta))==sggMtag(Serialized%eI(conta)  ,Serialized%eJ(conta)  ,Serialized%eK(conta)+1)).and. &
               !          (sggMtag(Serialized%eI(conta),Serialized%eJ(conta),Serialized%eK(conta))==sggMtag(Serialized%eI(conta)+1,Serialized%eJ(conta)+1,Serialized%eK(conta)+1)) ).or. &
               !        !
               !         ((sggMtag(Serialized%eI(conta),Serialized%eJ(conta),Serialized%eK(conta))==sggMtag(Serialized%eI(conta)  ,Serialized%eJ(conta)  ,Serialized%eK(conta)  )).and. &
               !          (sggMtag(Serialized%eI(conta),Serialized%eJ(conta),Serialized%eK(conta))==sggMtag(Serialized%eI(conta)  ,Serialized%eJ(conta)+1,Serialized%eK(conta)  )).and. &
               !          (sggMtag(Serialized%eI(conta),Serialized%eJ(conta),Serialized%eK(conta))==sggMtag(Serialized%eI(conta)  ,Serialized%eJ(conta)+1,Serialized%eK(conta)+1)) ) ) then
                     if (tamaniompi == 0) then !mantengo la dicotomia mpisize=0 nocero solo para degugeo comparando
                           print *,quienmpi,' writting original Mtag'
                           write (myunit,'(i7)')  sggMtag(Serialized%eI(conta),Serialized%eJ(conta),Serialized%eK(conta)) !!! esto estaba mal en MPI: bug OLD vtk 121090  !
                     else
                           write (myunit,'(i7)')  Serialized%sggMtag(conta)  
                     endif
               !      else
               !          write (myunit,'(i4)')  -1
               !      endif
                        
               !endif
         end do
         !!!fin info tag  

      write(myunit,'(a)') ' '
      
      CLOSE(myunit)




      return
   END SUBROUTINE write_VTKfile


   SUBROUTINE creaUnstructData(Serialized,  numberOfSerialized,sgg,Nodes,Numnodes,Elems,NumEdges,NumQuads,vtkindex)

      INTEGER (kind=4), intent(out):: numNodes,numQuads,numEdges
      real (kind= RKIND), allocatable, dimension(:,:), intent(out) :: Nodes
      integer (kind=4), allocatable, dimension(:,:), intent(out) ::  Elems
      
      logical, intent(IN) :: vtkindex
      type (SGGFDTDINFO), intent(IN)    :: sgg
      INTEGER (kind=4), intent(in):: numberOfSerialized
      type (Serialized_t), intent(in)  ::  Serialized
      
      CHARACTER (LEN=BUFSIZE) :: buff ! File name
      INTEGER (kind=4):: conta


      numNodes=-1
      numEdges=-1
      numQuads=-1
      !creo por demas
      if (numberOfSerialized/=0) then
         Allocate (Nodes(0:numberOfSerialized * 4,3) )
         allocate (Elems(1:numberOfSerialized, 4) )
      else
         return
      endif
      !





      do conta=1,numberOfSerialized
         select case (Serialized%currentType(conta))
            case (iJx)
            numNodes=numNodes+1
            if (vtkindex) then
               Nodes(numNodes,1)=(Serialized%eI(conta))*1.0_RKIND
               Nodes(numNodes,2)=(Serialized%eJ(conta))*1.0_RKIND
               Nodes(numNodes,3)=(Serialized%eK(conta))*1.0_RKIND
               numNodes=numNodes+1
               Nodes(numNodes,1)=(1 + Serialized%eI(conta))*1.0_RKIND
               Nodes(numNodes,2)=(    Serialized%eJ(conta))*1.0_RKIND
               Nodes(numNodes,3)=(    Serialized%eK(conta))*1.0_RKIND
            else
               Nodes(numNodes,1)=sgg%LineX(Serialized%eI(conta))
               Nodes(numNodes,2)=sgg%Liney(Serialized%eJ(conta))
               Nodes(numNodes,3)=sgg%Linez(Serialized%eK(conta))
               numNodes=numNodes+1
               Nodes(numNodes,1)=sgg%LineX(1 + Serialized%eI(conta))
               Nodes(numNodes,2)=sgg%Liney(    Serialized%eJ(conta))
               Nodes(numNodes,3)=sgg%Linez(    Serialized%eK(conta))
            endif
            !
            numEdges=numEdges + 1
            Elems(conta,1)=NumNodes - 1
            Elems(conta,2)=NumNodes
            Elems(conta,3)=-1 !marcar como edge para luego escribir bien
            case (iJy)
            numNodes=numNodes+1
            if (vtkindex) then
               Nodes(numNodes,1)=(Serialized%eI(conta))*1.0_RKIND
               Nodes(numNodes,2)=(Serialized%eJ(conta))*1.0_RKIND
               Nodes(numNodes,3)=(Serialized%eK(conta))*1.0_RKIND
               numNodes=numNodes+1
               Nodes(numNodes,1)=(    Serialized%eI(conta))*1.0_RKIND
               Nodes(numNodes,2)=(1 + Serialized%eJ(conta))*1.0_RKIND
               Nodes(numNodes,3)=(    Serialized%eK(conta))*1.0_RKIND
            else
               Nodes(numNodes,1)=sgg%LineX(Serialized%eI(conta))
               Nodes(numNodes,2)=sgg%Liney(Serialized%eJ(conta))
               Nodes(numNodes,3)=sgg%Linez(Serialized%eK(conta))
               numNodes=numNodes+1
               Nodes(numNodes,1)=sgg%LineX(    Serialized%eI(conta))
               Nodes(numNodes,2)=sgg%Liney(1 + Serialized%eJ(conta))
               Nodes(numNodes,3)=sgg%Linez(    Serialized%eK(conta))
            endif
            !

            numEdges=numEdges + 1
            Elems(conta,1)=NumNodes - 1
            Elems(conta,2)=NumNodes
            Elems(conta,3)=-1 !marcar como edge para luego escribir bien
            case (iJz)
            numNodes=numNodes+1
            if (vtkindex) then
               Nodes(numNodes,1)=(Serialized%eI(conta))*1.0_RKIND
               Nodes(numNodes,2)=(Serialized%eJ(conta))*1.0_RKIND
               Nodes(numNodes,3)=(Serialized%eK(conta))*1.0_RKIND
               numNodes=numNodes+1
               Nodes(numNodes,1)=(    Serialized%eI(conta))*1.0_RKIND
               Nodes(numNodes,2)=(    Serialized%eJ(conta))*1.0_RKIND
               Nodes(numNodes,3)=(1 + Serialized%eK(conta))*1.0_RKIND
            else
               Nodes(numNodes,1)=sgg%LineX(Serialized%eI(conta))
               Nodes(numNodes,2)=sgg%Liney(Serialized%eJ(conta))
               Nodes(numNodes,3)=sgg%Linez(Serialized%eK(conta))
               numNodes=numNodes+1
               Nodes(numNodes,1)=sgg%LineX(    Serialized%eI(conta))
               Nodes(numNodes,2)=sgg%Liney(    Serialized%eJ(conta))
               Nodes(numNodes,3)=sgg%Linez(1 + Serialized%eK(conta))
            endif
            !
            numEdges=numEdges + 1
            Elems(conta,1)=NumNodes - 1
            Elems(conta,2)=NumNodes
            Elems(conta,3)=-1 !marcar como edge para luego escribir bien
            !
            case (iBloqueJx)
            numNodes=numNodes+1
            if (vtkindex) then
               Nodes(numNodes,1)=(Serialized%eI(conta))*1.0_RKIND
               Nodes(numNodes,2)=(Serialized%eJ(conta))*1.0_RKIND
               Nodes(numNodes,3)=(Serialized%eK(conta))*1.0_RKIND
            else
               Nodes(numNodes,1)=sgg%LineX(Serialized%eI(conta))
               Nodes(numNodes,2)=sgg%Liney(Serialized%eJ(conta))
               Nodes(numNodes,3)=sgg%Linez(Serialized%eK(conta))
            endif
            numNodes=numNodes+1
            if (vtkindex) then
               Nodes(numNodes,1)=(    Serialized%eI(conta))*1.0_RKIND
               Nodes(numNodes,2)=(1 + Serialized%eJ(conta))*1.0_RKIND
               Nodes(numNodes,3)=(    Serialized%eK(conta))*1.0_RKIND
            else
               Nodes(numNodes,1)=sgg%LineX(    Serialized%eI(conta))
               Nodes(numNodes,2)=sgg%Liney(1 + Serialized%eJ(conta))
               Nodes(numNodes,3)=sgg%Linez(    Serialized%eK(conta))
            endif
            numNodes=numNodes+1
            if (vtkindex) then
               Nodes(numNodes,1)=(    Serialized%eI(conta))*1.0_RKIND
               Nodes(numNodes,2)=(1 + Serialized%eJ(conta))*1.0_RKIND
               Nodes(numNodes,3)=(1 + Serialized%eK(conta))*1.0_RKIND
            else
               Nodes(numNodes,1)=sgg%LineX(    Serialized%eI(conta))
               Nodes(numNodes,2)=sgg%Liney(1 + Serialized%eJ(conta))
               Nodes(numNodes,3)=sgg%Linez(1 + Serialized%eK(conta))
            endif
            numNodes=numNodes+1
            if (vtkindex) then
               Nodes(numNodes,1)=(    Serialized%eI(conta))*1.0_RKIND
               Nodes(numNodes,2)=(    Serialized%eJ(conta))*1.0_RKIND
               Nodes(numNodes,3)=(1 + Serialized%eK(conta))*1.0_RKIND
            else
               Nodes(numNodes,1)=sgg%LineX(    Serialized%eI(conta))
               Nodes(numNodes,2)=sgg%Liney(    Serialized%eJ(conta))
               Nodes(numNodes,3)=sgg%Linez(1 + Serialized%eK(conta))
            endif
            !

            numQuads=numQuads + 1
            Elems(conta,1)=NumNodes - 3
            Elems(conta,2)=NumNodes - 2
            Elems(conta,3)=NumNodes - 1
            Elems(conta,4)=NumNodes
            case (iBloqueJy)
            numNodes=numNodes+1
            if (vtkindex) then
               Nodes(numNodes,1)=(Serialized%eI(conta))*1.0_RKIND
               Nodes(numNodes,2)=(Serialized%eJ(conta))*1.0_RKIND
               Nodes(numNodes,3)=(Serialized%eK(conta))*1.0_RKIND
            else
               Nodes(numNodes,1)=sgg%LineX(Serialized%eI(conta))
               Nodes(numNodes,2)=sgg%Liney(Serialized%eJ(conta))
               Nodes(numNodes,3)=sgg%Linez(Serialized%eK(conta))
            endif
            numNodes=numNodes+1
            if (vtkindex) then
               Nodes(numNodes,1)=(1 + Serialized%eI(conta))*1.0_RKIND
               Nodes(numNodes,2)=(    Serialized%eJ(conta))*1.0_RKIND
               Nodes(numNodes,3)=(    Serialized%eK(conta))*1.0_RKIND
            else
               Nodes(numNodes,1)=sgg%LineX(1 + Serialized%eI(conta))
               Nodes(numNodes,2)=sgg%Liney(    Serialized%eJ(conta))
               Nodes(numNodes,3)=sgg%Linez(    Serialized%eK(conta))
            endif
            numNodes=numNodes+1
            if (vtkindex) then
               Nodes(numNodes,1)=(1 + Serialized%eI(conta))*1.0_RKIND
               Nodes(numNodes,2)=(    Serialized%eJ(conta))*1.0_RKIND
               Nodes(numNodes,3)=(1 + Serialized%eK(conta))*1.0_RKIND
            else
               Nodes(numNodes,1)=sgg%LineX(1 + Serialized%eI(conta))
               Nodes(numNodes,2)=sgg%Liney(    Serialized%eJ(conta))
               Nodes(numNodes,3)=sgg%Linez(1 + Serialized%eK(conta))
            endif
            numNodes=numNodes+1
            if (vtkindex) then
               Nodes(numNodes,1)=(    Serialized%eI(conta))*1.0_RKIND
               Nodes(numNodes,2)=(    Serialized%eJ(conta))*1.0_RKIND
               Nodes(numNodes,3)=(1 + Serialized%eK(conta))*1.0_RKIND
            else
               Nodes(numNodes,1)=sgg%LineX(    Serialized%eI(conta))
               Nodes(numNodes,2)=sgg%Liney(    Serialized%eJ(conta))
               Nodes(numNodes,3)=sgg%Linez(1 + Serialized%eK(conta))
            endif
            !
            numQuads=numQuads + 1
            Elems(conta,1)=NumNodes - 3
            Elems(conta,2)=NumNodes - 2
            Elems(conta,3)=NumNodes - 1
            Elems(conta,4)=NumNodes
            case (iBloqueJz)
            numNodes=numNodes+1
            if (vtkindex) then
               Nodes(numNodes,1)=(Serialized%eI(conta))*1.0_RKIND
               Nodes(numNodes,2)=(Serialized%eJ(conta))*1.0_RKIND
               Nodes(numNodes,3)=(Serialized%eK(conta))*1.0_RKIND
            else
               Nodes(numNodes,1)=sgg%LineX(Serialized%eI(conta))
               Nodes(numNodes,2)=sgg%Liney(Serialized%eJ(conta))
               Nodes(numNodes,3)=sgg%Linez(Serialized%eK(conta))
            endif
            numNodes=numNodes+1
            if (vtkindex) then
               Nodes(numNodes,1)=(1 + Serialized%eI(conta))*1.0_RKIND
               Nodes(numNodes,2)=(    Serialized%eJ(conta))*1.0_RKIND
               Nodes(numNodes,3)=(    Serialized%eK(conta))*1.0_RKIND
            else
               Nodes(numNodes,1)=sgg%LineX(1 + Serialized%eI(conta))
               Nodes(numNodes,2)=sgg%Liney(    Serialized%eJ(conta))
               Nodes(numNodes,3)=sgg%Linez(    Serialized%eK(conta))
            endif
            numNodes=numNodes+1
            if (vtkindex) then
               Nodes(numNodes,1)=(1 + Serialized%eI(conta))*1.0_RKIND
               Nodes(numNodes,2)=(1 + Serialized%eJ(conta))*1.0_RKIND
               Nodes(numNodes,3)=(    Serialized%eK(conta))*1.0_RKIND
            else
               Nodes(numNodes,1)=sgg%LineX(1 + Serialized%eI(conta))
               Nodes(numNodes,2)=sgg%Liney(1 + Serialized%eJ(conta))
               Nodes(numNodes,3)=sgg%Linez(    Serialized%eK(conta))
            endif
            numNodes=numNodes+1
            if (vtkindex) then
               Nodes(numNodes,1)=(    Serialized%eI(conta))*1.0_RKIND
               Nodes(numNodes,2)=(1 + Serialized%eJ(conta))*1.0_RKIND
               Nodes(numNodes,3)=(    Serialized%eK(conta))*1.0_RKIND
            else
               Nodes(numNodes,1)=sgg%LineX(    Serialized%eI(conta))
               Nodes(numNodes,2)=sgg%Liney(1 + Serialized%eJ(conta))
               Nodes(numNodes,3)=sgg%Linez(    Serialized%eK(conta))
            endif
            !

            numQuads=numQuads + 1
            Elems(conta,1)=NumNodes - 3
            Elems(conta,2)=NumNodes - 2
            Elems(conta,3)=NumNodes - 1
            Elems(conta,4)=NumNodes
         end select
      end do

      if ((NumEdges+1)+(NumQuads+1)/=numberofSerialized) then
         buff='ERROR: Buggy error sumas creating .vtk'
         CALL print11(0_4, buff)
      endif

      return
   END SUBROUTINE creaUnstructData

   !subroutine fillinparaviewstate
   !
   !return
   !end subroutine
END MODULE VTK
!
!
