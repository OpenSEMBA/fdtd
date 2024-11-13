MODULE xdmf
   !
   USE fdetypes
   USE Observa
   use report
   use xdmf_h5
   !
   !
   !
   !
   IMPLICIT NONE
   !
   PRIVATE
   PUBLIC createxdmf,createxdmfOnTheFly,createh5bintxt
   !!!        public create_interpreted_mesh
CONTAINS
   ! =================================================================================================
   !                       ====>>>> SALVADOR's CODE <<<<====
   ! =================================================================================================
   !
   !Subrutine to parse the volumic probes to create .xdmf and .h5 files
   !
   SUBROUTINE createxdmf (sgg,layoutnumber, size,vtkindex,createh5bin,somethingdone,mpidir)
      logical, save :: firsttimeenteringcreatexdmf=.true.
      integer (KIND=4) :: mpidir
      logical :: vtkindex,createh5bin
      !------------------------>
      CHARACTER (LEN=BUFSIZE) :: filename ! File name
      !
      type (SGGFDTDINFO), intent(IN)        :: sgg
      INTEGER (KIND=4), INTENT (IN) :: layoutnumber, size
      INTEGER (KIND=4) ::  ierr,  sizeofvalores,COMPO
      complex( kind = CKIND), dimension( :, :, :, :,: ), allocatable  :: valor3DComplex !freqdomain probes
      !
      TYPE (output_t), POINTER, DIMENSION (:) :: output
      INTEGER (KIND=4) :: iroot
      integer (KIND=4) :: myunit
      !
#ifdef CompileWithMPI
      REAL (KIND=RKIND), ALLOCATABLE, DIMENSION (:, :, :, :) :: newvalor3d !para sondas Volumic
#endif


      REAL (KIND=RKIND), ALLOCATABLE, DIMENSION (:, :, :, :) :: valor3d !para sondas Volumic
      real (  KINd=RKIND_TIEMPO), ALLOCATABLE, DIMENSION (:) :: att


      INTEGER (KIND=4) :: indi,fieldob
      !
      INTEGER (KIND=4) :: ii, i1, j1, k1, finalstep
      INTEGER (KIND=4) :: minx, maxx, miny, maxy, minz, maxz,pasadas,pasadastotales
      LOGICAL :: lexis,somethingdone
      character (LEN=BUFSIZE)     ::  dubuf
      INTEGER (KIND=4) :: minXabs, maxXabs, minYabs, maxYabs, minZabs, maxZabs 
      INTEGER (KIND=4) :: minXabs_primero,minYabs_primero,minZabs_primero,imdice
      CHARACTER (LEN=BUFSIZE) :: pathroot
      character (LEN=BUFSIZE)  ::  chari,charj,chark,chari2,charj2,chark2
      character (LEN=BUFSIZE)  ::  extpoint
      character(len=BUFSIZE) :: buff
      REAL (KIND=RKIND) :: linez_minZabs_primero,liney_minYabs_primero,linex_minXabs_primero, &
                             dz_minZabs,dy_minYabs,dx_minXabs                 
      !
      CHARACTER (LEN=BUFSIZE) :: whoami, whoamishort
      REAL (KIND=RKIND)  :: rdum
      integer :: my_iostat
      
      WRITE (whoamishort, '(i5)') layoutnumber + 1
      WRITE (whoami, '(a,i5,a,i5,a)') '(', layoutnumber + 1, '/', size, ') '
      !
      output => GetOutput ()!get the output private info from observation

      somethingdone=.false.
      barridoprobes: DO ii = 1, sgg%NumberRequest

         IF (sgg%observation(ii)%Volumic) then
         if (sgg%observation(ii)%nP == 1) then
         if ((sgg%observation(ii)%P(1)%What /= nothing).AND.(sgg%observation(ii)%P(1)%What /= iCur).AND.(sgg%observation(ii)%P(1)%What /= mapvtk).AND. &
             (sgg%observation(ii)%P(1)%What  /=  iCurX).AND.(sgg%observation(ii)%P(1)%What  /=  iCurY).AND.(sgg%observation(ii)%P(1)%What  /=  iCurZ)) THEN
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
         !
         !sondas Volumic traducelas a xdfm
         IF (sgg%observation(ii)%Volumic) then
            if (sgg%observation(ii)%nP == 1) then
               if ((sgg%observation(ii)%P(1)%What /= nothing).AND.(sgg%observation(ii)%P(1)%What /= iCur).AND.(sgg%observation(ii)%P(1)%What /= mapvtk).AND. &
               (sgg%observation(ii)%P(1)%What  /=  iCurX).AND.(sgg%observation(ii)%P(1)%What  /=  iCurY).AND.(sgg%observation(ii)%P(1)%What  /=  iCurZ)) THEN
                  INQUIRE (FILE=trim(adjustl(output(ii)%item(1)%path)), EXIST=lexis)
                  if ((lexis).and.(output(ii)%TimesWritten/=0)) then
                     fieldob=sgg%observation(ii)%P(1)%what
!inicializaciones varias
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
                    !mpidir 190319   !desrotacion para que los nombres sean correctos 
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
                      
                     !! CORREGIDO PARA TRANCOS   AHORA DESPUES DE HABER PUESTO BIEN EXTPOINT
!guarda el original   
                     minXabs_primero = minXabs
                     minYabs_primero = minYabs   
                     minZabs_primero = minZabs
                     im1: do imdice=minXabs,maxXabs
                         if (mod(imdice,output(ii)%item(1)%Xtrancos)==0) then
                            minXabs_primero=imdice
                            exit im1
                        endif
                     end do im1
                     im2: do imdice=minYabs,maxYabs
                         if (mod(imdice,output(ii)%item(1)%Ytrancos)==0) then
                            minYabs_primero=imdice
                            exit im2
                        endif
                     end do im2
                     im3: do imdice=minZabs,maxZabs
                         if (mod(imdice,output(ii)%item(1)%Ztrancos)==0) then
                            minZabs_primero=imdice
                            exit im3
                        endif
                     end do im3
!pufff hay mucha reduncancia minxabs = minx, etc. 021219 limpiar algun dia 
                     minXabs = int(sgg%Observation(ii)%P(1)%XI/output(ii)%item(1)%Xtrancos)
                     if (mod(sgg%Observation(ii)%P(1)%XI,output(ii)%item(1)%Xtrancos) /= 0) minXabs=minXabs+1
                     maxXabs = int(sgg%Observation(ii)%P(1)%XE/output(ii)%item(1)%Xtrancos)
                     minYabs = int(sgg%observation(ii)%P(1)%YI/output(ii)%item(1)%Ytrancos)
                     if (mod(sgg%Observation(ii)%P(1)%YI,output(ii)%item(1)%Ytrancos) /= 0) minYabs=minYabs+1
                     maxYabs = int(sgg%observation(ii)%P(1)%YE/output(ii)%item(1)%Ytrancos)    
                     
#ifdef CompileWithMPI
                     minZabs = int(output(ii)%item(1)%ZIorig/output(ii)%item(1)%Ztrancos)
                     if (mod(output(ii)%item(1)%ZIorig,output(ii)%item(1)%Ztrancos) /= 0) minZabs=minZabs+1
                     maxZabs = int(output(ii)%item(1)%ZEorig/output(ii)%item(1)%Ztrancos)
#else
                     minZabs = int(sgg%observation(ii)%P(1)%zI/output(ii)%item(1)%Ztrancos)
                     if (mod(sgg%Observation(ii)%P(1)%ZI,output(ii)%item(1)%Ztrancos) /= 0) minZabs=minZabs+1
                     maxZabs = int(sgg%observation(ii)%P(1)%zE/output(ii)%item(1)%Ztrancos)
#endif

                        
                     !fin trancos
                      
                     iroot=index(output(ii)%item(1)%path,'__',.true.)
                     pathroot=trim(adjustl(output(ii)%item(1)%path(1:iroot-1)))
                     iroot = index (pathroot, '_',.true.)
                     pathroot = trim (adjustl(pathroot(1:iroot-1)))
                     iroot = index (pathroot, '_',.true.)
                     pathroot = trim (adjustl(pathroot(1:iroot-1)))
                     iroot = index (pathroot, '_',.true.)
                     pathroot = trim(adjustl(pathroot(1:iroot-1)))//'_'//trim(adjustl(extpoint))


                     linez_minZabs_primero     = sgg%linez(minZabs_primero)
                     liney_minYabs_primero     = sgg%liney(minYabs_primero)
                     linex_minXabs_primero     = sgg%linex(minXabs_primero) 
                     dz_minZabs                = sgg%dz(minZabs)*output(ii)%item(1)%Ztrancos
                     dy_minYabs                = sgg%dy(minYabs)*output(ii)%item(1)%Ytrancos
                     dx_minXabs                = sgg%dx(minXabs)*output(ii)%item(1)%Xtrancos
                     
                     OPEN (output(ii)%item(1)%UNIT, FILE=trim(adjustl(output(ii)%item(1)%path)), FORM='unformatted')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                      
                     READ (output(ii)%item(1)%UNIT) minx, maxx , miny, maxy, minz, maxz    !ya deben venir bien escritos incluyendo correccion de TRANCOS
!!!allocate space                     
                     if (SGG%Observation(ii)%TimeDomain) then
                       finalstep=output(ii)%TimesWritten
                       allocate (att(1:finalstep))
                          att = 0.0_RKIND !aquiiiii
                       pasadastotales=1
                       ALLOCATE (valor3d(minXabs:maxXabs, minYabs:maxYabs, minZabs:maxZabs, 1))
                     elseif (SGG%Observation(ii)%FreqDomain) then
                       !este read solo se precisa para la frecuencial y es dummy 
                       read(output(ii)%item(1)%unit) rdum !instante en el que se ha escrito la info frequencial
                       finalstep= output(ii)%NumFreqs
                       allocate (att(1:finalstep))
                          att = 0.0_RKIND  !aquiiiii
                       pasadastotales=2
                       ALLOCATE (valor3d(minXabs:maxXabs, minYabs:maxYabs, minZabs:maxZabs, 1))
                       ALLOCATE (valor3dCOMPLEX(1,1:3,minXabs:maxXabs, minYabs:maxYabs, minZabs:maxZabs))
                     endif

#ifdef CompileWithMPI
                     ALLOCATE (newvalor3d(minXabs:maxXabs, minYabs:maxYabs, minZabs:maxZabs, 1))
#endif

#ifdef CompileWithMPI
                     IF (layoutnumber == output(ii)%item(1)%MPIRoot) THEN
#else                 
                     IF (layoutnumber == 0) THEN
#endif                
                         if (createh5bin) then
                            if (firsttimeenteringcreatexdmf) then         
                                open(newunit=myunit,file=trim(adjustl(sgg%nEntradaRoot))//'_'//trim(adjustl(whoamishort))//'_h5bin.txt',form='formatted') !lista de todos los .h5bin     
                                WRITE (myunit, '(a)') '!END'      
                                close(myunit,status='delete')
                                firsttimeenteringcreatexdmf=.false.
                            endif        
                            my_iostat=0
9138                        if(my_iostat /= 0) write(*,fmt='(a)',advance='no'), '.' !!if(my_iostat /= 0) print '(i5,a1,i4,2x,a)',9138,'.',layoutnumber,trim(adjustl(sgg%nEntradaRoot))//'_'//trim(adjustl(whoamishort))//'_h5bin.txt'
                            open(newunit=myunit,file=trim(adjustl(sgg%nEntradaRoot))//'_'//trim(adjustl(whoamishort))//'_h5bin.txt',form='formatted',position='append',err=9138,iostat=my_iostat,status='new',action='write') !lista de todos los .h5bin   
                ! !lista de todos los .h5bin
                            write (myunit,'(a)') trim(adjustl(pathroot))//'.h5bin'
                            close(myunit)
                            !
                            open(newunit=myunit,file=trim(adjustl(pathroot))//'.h5bin',form='unformatted')
                            write (myunit) finalstep,minXabs, maxXabs, minYabs, maxYabs, minZabs, maxZabs,fieldob,SGG%Observation(ii)%TimeDomain,pasadastotales
                         endif
                     endif !del layoutnumber
                        
#ifdef CompileWithMPI
                     call MPI_Barrier(output(ii)%item(1)%MPISubComm,ierr)
#endif                 
                     buclepasadas: do pasadas=1,pasadastotales  
!!!inicializa a cero en cada pasada
                       if (SGG%Observation(ii)%TimeDomain) then
                          valor3d = 0.0_RKIND
                        elseif (SGG%Observation(ii)%FreqDomain) then
                          valor3d = 0.0_RKIND
                          valor3dCOMPLEX = 0.0_RKIND
                        endif
                        
#ifdef CompileWithMPI
                        newvalor3d = 0.0_RKIND
#endif
!!!!abre fiecho escritura .h5 
                   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                                   
                        if (SGG%Observation(ii)%TimeDomain) then         
                            if (pasadas==1) then   
                                filename = trim (adjustl(pathroot))//'_time'
                            else
                                print *,'Buggy error in valor3d. '
                                stop
                            endif               
                            continue !ya se ha leido valor3d
                        else 
                            if (pasadas==1) then   
                                filename = trim (adjustl(pathroot))//'_mod'
                            elseif (pasadas==2) then  
                                filename = trim (adjustl(pathroot))//'_phase'
                            else
                                print *,'Buggy error in valor3d. '
                                stop
                            endif
                        endif
#ifdef CompileWithHDF
#ifdef CompileWithMPI
                        IF (layoutnumber == output(ii)%item(1)%MPIRoot) THEN
#else                 
                        IF (layoutnumber == 0) THEN
#endif                
                           
                           if (.not.(((fieldob == iMEC).or.(fieldob ==iMHC)).and.(pasadas ==2))) then ! no tiene sentido esccribir la fase del modulo
                              call openh5file(filename,finalstep,minXabs,maxXabs, minYabs,maxYabs, minZabs,maxZabs)
                           endif
                        endif
#endif               
                           
                        bucleindi: DO indi = 1, finalstep
                               if (pasadas == 1) then !solo es preciso leer los datos una vez
                                   READ (output(ii)%item(1)%UNIT) att(indi)
                                   write(dubuf,*)  ' ----> .xdmf file ',att(indi),'(',indi,'/',finalstep,')'
                                   call print11(layoutnumber,dubuf)
                              
                                   if (SGG%Observation(ii)%TimeDomain) then
                                     DO k1 = minz, maxz
                                        DO j1 = miny, maxy
                                           READ (output(ii)%item(1)%UNIT) (valor3d(i1, j1, k1, 1), i1=minx, maxx)
                                        END DO
                                     END DO
                                     !
                                   elseif (SGG%Observation(ii)%FreqDomain) then
                                     DO COMPO=1,3
                                        DO k1 = minz, maxz
                                           DO j1 = miny, maxy
                                              READ (output(ii)%item(1)%UNIT) (valor3dCOMPLEX(1,COMPO,i1, j1, k1), i1=minx, maxx)
                                           END DO
                                        END DO
                                     END DO
                                   endif
                               endif !del if (pasadas==1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                                   
                               if (SGG%Observation(ii)%TimeDomain) then                     
                                    continue !ya se ha leido valor3d
                               else !freqdomain construir valor3d
                                    select case (fieldob)
                                    case(iMEC,iMHC)
                                         !modulo 
                                         DO k1 = minz, maxz
                                            DO j1 = miny, maxy
                                               DO i1=minx, maxx
                                                  if (pasadas==1) then      !modulo
                                                     valor3d(i1, j1, k1, 1)=SQRT( ABS(valor3dCOMPLEX(1,1,i1, j1, k1))**2. + &
                                                                                  ABS(valor3dCOMPLEX(1,2,i1, j1, k1))**2. + &
                                                                                  ABS(valor3dCOMPLEX(1,3,i1, j1, k1))**2. )  !sgg  301119 faltaba este cuadrado creo
                                                  else !phase
                                                      valor3d=0.0_RKIND !LA fase no tiene sentido para el modulo del vector
                                                  endif
                                               END DO
                                            END DO
                                         END DO
                                    case(iExC,iHxC)
                                         DO k1 = minz, maxz
                                            DO j1 = miny, maxy
                                               DO i1=minx, maxx
                                                  if (pasadas==1) then    !modulo
                                                      valor3d(i1, j1, k1, 1)= ABS(valor3dCOMPLEX(1,1,i1, j1, k1)) 
                                                   else   !phase
                                                      valor3d(i1, j1, k1, 1)= ATAN2(AIMAG(valor3dCOMPLEX(1,1,i1, j1, k1)),REAL(valor3dCOMPLEX(1,1,i1, j1, k1)))
                                                   endif
                                               END DO
                                            END DO
                                         END DO
                                    case(iEyC,iHyC) 
                                         DO k1 = minz, maxz
                                            DO j1 = miny, maxy
                                               DO i1=minx, maxx
                                                  if (pasadas==1) then       !modulo 
                                                      valor3d(i1, j1, k1, 1)=ABS(valor3dCOMPLEX(1,2,i1, j1, k1))
                                                  else    !phase
                                                      valor3d(i1, j1, k1, 1)= ATAN2(AIMAG(valor3dCOMPLEX(1,2,i1, j1, k1)),REAL(valor3dCOMPLEX(1,2,i1, j1, k1)))
                                                  endif
                                               END DO
                                            END DO
                                         END DO
                                    case(iEzC,iHzC)
                                         DO k1 = minz, maxz
                                            DO j1 = miny, maxy
                                               DO i1=minx, maxx
                                                  if (pasadas==1) then      !modulo 
                                                     valor3d(i1, j1, k1, 1)=ABS(valor3dCOMPLEX(1,3,i1, j1, k1))
                                                  else                 !phase
                                                     valor3d(i1, j1, k1, 1)= ATAN2(AIMAG(valor3dCOMPLEX(1,3,i1, j1, k1)),REAL(valor3dCOMPLEX(1,3,i1, j1, k1)))
                                                  endif
                                               END DO
                                            END DO
                                         END DO
                                    case default
                                          print *,'Buggy error in valor3d. Not processing continuing. '
                                          continue
                                    end select                                                 
                               endif   !del time domain
!!!!!!!!!!!!!!!!sincroniza valor3d y aunalos en el root
#ifdef CompileWithMPI
                               if (size>1) then
                                  if (output(ii)%item(1)%MPISubComm /= -1) then
                                     sizeofvalores = (maxXabs-minXabs+1) * (maxYabs-minYabs+1) * (maxZabs-minZabs+1)
                                     call MPI_Barrier(output(ii)%item(1)%MPISubComm,ierr)
                                     CALL MPI_AllReduce (valor3d, newvalor3d, sizeofvalores, REALSIZE, MPI_SUM, &
                                     &                     output(ii)%item(1)%MPISubComm, ierr)
                                  endif
                                  valor3d = newvalor3d
                               endif
#endif
                    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!escribe los ficheros de salida
 
#ifdef CompileWithMPI
                              IF (layoutnumber == output(ii)%item(1)%MPIRoot) THEN
#else                   
                              IF (layoutnumber == 0) THEN
#endif                
#ifdef CompileWithHDF
                                   if (.not.(((fieldob == iMEC).or.(fieldob ==iMHC)).and.(pasadas ==2))) then ! no tiene sentido esccribir la fase del modulo
                                       call writeh5file(filename,valor3d,indi,att(indi),minXabs,maxXabs, minYabs,maxYabs, minZabs,maxZabs, &
                                                        linez_minZabs_primero,liney_minYabs_primero,linex_minXabs_primero, &
                                                        dz_minZabs,dy_minYabs,dx_minXabs,&
                                                        minZabs_primero,minYabs_primero,minXabs_primero,finalstep,vtkindex)
                                   endif
#endif                             
                                   if (createh5bin) then
                                       write (myunit) (minZabs_primero),(minYabs_primero), (minXabs_primero)      
                                       write (myunit) linez_minZabs_primero,liney_minYabs_primero,linex_minXabs_primero
                                       write (myunit) dz_minZabs,dy_minYabs,dx_minXabs
                                       WRITE (myunit) att(indi)
                                       DO k1 = minzabs, maxzabs
                                          DO j1 = minyabs, maxyabs
                                             WRITE (myunit) (valor3d(i1, j1, k1, 1), i1=minxabs, maxxabs)
                                          END DO
                                       END DO
                                   endif
                              endif                  
                          
                        END DO bucleindi
                              
#ifdef CompileWithHDF
#ifdef CompileWithMPI
                        IF (layoutnumber == output(ii)%item(1)%MPIRoot) THEN
#else                  
                        IF (layoutnumber == 0) THEN
#endif   
                           
                           if (.not.(((fieldob == iMEC).or.(fieldob ==iMHC)).and.(pasadas ==2))) then ! no tiene sentido esccribir la fase del modulo
                              call closeh5file(finalstep,att)
                              CALL print11 (layoutnumber, trim(adjustl(whoami))//' Written into '//trim(adjustl(filename))//'.h5', .TRUE.) !enforces print
                           endif
                        endif                        
#endif                 
                     end do buclepasadas
                   !
#ifdef CompileWithMPI
                     IF (layoutnumber == output(ii)%item(1)%MPIRoot) THEN
#else              
                     IF (layoutnumber == 0) THEN
#endif   
                        if (createh5bin) then
                             close(myunit)
                             CALL print11 (layoutnumber, trim(adjustl(whoami))//' Written into '//trim(adjustl(sgg%nEntradaRoot))//'.h5bin', .TRUE.)
                        endif
                     endif                     

                     DEALLOCATE (valor3d)
                     if (SGG%Observation(ii)%FreqDomain) then
                        DEALLOCATE (valor3dCOMPLEX)
                     ENDIF
#ifdef CompileWithMPI
                     DEALLOCATE (newvalor3d)
#endif
                     DEALLOCATE (ATT)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                     
                     CLOSE (output(ii)%item(1)%UNIT)
                  else !del lexis
                     buff='NOT PROCESSING: Ignoring: Inexistent or void file '//trim(adjustl(output(ii)%item(1)%path))
                     CALL print11(layoutnumber, buff)
                  ENDIF !DEL LEXIS
               somethingdone=.true.
               ENDIF
            ENDIF
         ENDIF
      END DO barridoprobes !barrido puntos de observacion
                     
      RETURN
   END SUBROUTINE createxdmf


   SUBROUTINE createh5bintxt(sgg,layoutnumber,size)
      type (SGGFDTDINFO), intent(IN)        :: sgg
      INTEGER (KIND=4), INTENT (IN) :: layoutnumber, size
      logical :: lexis,algoescrito
      INTEGER (KIND=4) :: ii,ierr
      integer (KIND=4) :: myunit,myunit2
      CHARACTER (LEN=BUFSIZE) ::  whoamishort
      CHARACTER (LEN=BUFSIZE) :: pathroot
      integer :: my_iostat
#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif              
      if (layoutnumber == 0) then !solo el root
          open(newunit=myunit,file=trim(adjustl(sgg%nEntradaRoot))//'_h5bin.txt',form='formatted') !lista de todos los .h5bin
          write (myunit,'(a)') '!END'
          close(myunit,status='delete')  
          my_iostat=0
9138      if(my_iostat /= 0) write(*,fmt='(a)',advance='no'), '.' !!if(my_iostat /= 0) print '(i5,a1,i4,2x,a)',9138,'.',layoutnumber,trim(adjustl(sgg%nEntradaRoot))//'_h5bin.txt'
          open(newunit=myunit,file=trim(adjustl(sgg%nEntradaRoot))//'_h5bin.txt',form='formatted',err=9138,iostat=my_iostat,status='new',action='write') !lista de todos los .h5bin   
          algoescrito=.false.
          do ii=0,size-1 !auna todos los _h5bin.txt 
             WRITE (whoamishort, '(i5)') ii + 1
             INQUIRE (FILE=trim(adjustl(trim(adjustl(sgg%nEntradaRoot))//'_'//trim(adjustl(whoamishort))//'_h5bin.txt')), EXIST=lexis)
             if (lexis) then
                 open(newunit=myunit2,file=trim(adjustl(sgg%nEntradaRoot))//'_'//trim(adjustl(whoamishort))//'_h5bin.txt',form='formatted')
                 do  
                     read  (myunit2, '(a)',end=9874) pathroot           
                     write (myunit,'(a)') trim(adjustl(pathroot))  
                     algoescrito=.true.
                 end do
9874            close (myunit2,status='delete')    
             endif      
          end do
          if (algoescrito) then          
            close(myunit)
          else
            close(myunit,status='delete')
          endif
      endif
#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif     
   end SUBROUTINE createh5bintxt

   SUBROUTINE createxdmfOnTheFly (sgg,layoutnumber,size,vtkindex,createh5bin,somethingdone,mpidir)
      integer (KIND=4) :: mpidir
      logical :: vtkindex,createh5bin
      !------------------------>

      type (SGGFDTDINFO), intent(IN)        :: sgg
      INTEGER (KIND=4), INTENT (IN) :: layoutnumber, size
      TYPE (output_t), POINTER, DIMENSION (:) :: output
      INTEGER (KIND=4) :: ii
      logical :: lexis,somethingdone
      character(len=BUFSIZE) :: buff
      !
      output => GetOutput ()!get the output private info from observation
      !

      DO ii = 1, sgg%NumberRequest
         !sondas Volumic traducelas a xdfm
         IF (sgg%observation(ii)%Volumic) then
            if (sgg%observation(ii)%nP == 1) then
               if ((sgg%observation(ii)%P(1)%What /= nothing).AND.(sgg%observation(ii)%P(1)%What /= iCur).AND.(sgg%observation(ii)%P(1)%What  /=  iCurX).AND.(sgg%observation(ii)%P(1)%What  /=  iCurY).AND.(sgg%observation(ii)%P(1)%What  /=  iCurZ)) THEN
                  !
                  INQUIRE (FILE=trim(adjustl(output(ii)%item(1)%path)), EXIST=lexis)
                  if (.not.lexis) then
                     buff='NOT PROCESSING: Inexistent file '//trim(adjustl(output(ii)%item(1)%path))
                     CALL print11(layoutnumber, buff)
                     return
                  ELSE
                     close (output(ii)%item(1)%unit)
                  ENDIF !DEL LEXIS
               ENDIF
            ENDIF
         ENDIF

      END DO !barrido puntos de observacion
      call createxdmf (sgg,layoutnumber, size,vtkindex,createh5bin,somethingdone,mpidir)
      DO ii = 1, sgg%NumberRequest
         !sondas Volumic traducelas a xdfm
         IF (sgg%observation(ii)%Volumic) then
            if (sgg%observation(ii)%nP == 1) then
               if ((sgg%observation(ii)%P(1)%What /= nothing).AND.(sgg%observation(ii)%P(1)%What /= iCur).AND.(sgg%observation(ii)%P(1)%What  /=  iCurX).AND.(sgg%observation(ii)%P(1)%What  /=  iCurY).AND.(sgg%observation(ii)%P(1)%What  /=  iCurZ)) THEN
                  !
                  INQUIRE (FILE=trim(adjustl(output(ii)%item(1)%path)), EXIST=lexis)
                  if (.not.lexis) then
                     buff='NOT PROCESSING: Inexistent file '//trim(adjustl(output(ii)%item(1)%path))
                     CALL print11(layoutnumber, buff)
                     return
                  ELSE
                     open (output(ii)%item(1)%unit,file=trim(adjustl(output(ii)%item(1)%path)),FORM='unformatted',position='append')
                  ENDIF !DEL LEXIS
               ENDIF
            ENDIF
         ENDIF

      END DO !barrido puntos de observacion

      RETURN
   END SUBROUTINE createxdmfOnTheFly
END MODULE xdmf
!
!
