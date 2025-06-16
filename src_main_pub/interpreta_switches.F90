module interpreta_switches_m  

   USE FDETYPES
   USE Getargs   
   use EpsMuTimeScale_m
   use Report
   use version
 
   IMPLICIT NONE
   PRIVATE
   !   
   type entrada_t

        logical ::                            &
            forcing                         , &
            singlefilewrite                 , &
            ignoresamplingerrors            , &
            ignoreerrors                    , &
            updateshared                    , &
            prioritizeISOTROPICBODYoverall  , &
            wirecrank                       , &
            CLIPREGION                      , &
            verbose                         , &
            resume                          , &
            forcesteps                      , &
            resume_fromold                  , &
            freshstart                      , &
            run                             , &
            createmap                       , &
            dontwritevtk                    , &
            vtkindex                        , &
            createmapvtk                    , &
            hopf                            , &
            run_with_dmma                   , &
            run_with_abrezanjas             , &
            input_conformal_flag            , &
            pausar                          , &
            relaunching                     , &
            forcestop                       , &
            l_aux                           , &
            flag_conf_sgg                   , &
            takeintcripte                   , &
            skindepthpre                    , &
            sgbc                            , &
            conformalskin                   , &
            ade                             , &
            mibc                            , &
            NOcompomur                      , &
            MurAfterPML                     , &
            sgbccrank                       , &
            sgbcDispersive                  , &
            saveall                         , &
            boundwireradius                 , &
            hay_slanted_wires               , &
            makeholes                       , &
            mur_first                       , &    
            mur_second                      , &
            mur_exist                       , &
            connectendings                  , &
            strictOLD                       , &
            mtlnberenger                    , &
            stableradholland                , &
            TAPARRABOS                      , &
            fieldtotl                       , &
            forceresampled                  , &
            isolategroupgroups              , &
            groundwires                     , &
            noSlantedcrecepelo              , &
            forcecfl                        , &
            niapapostprocess                , &
            permitscaling                   , &
            stochastic                      , &
            chosenyesornostochastic         , &
            prioritizeCOMPOoverPEC          , &
            prioritizeTHINWIRE              , &
            createh5bin                     , &
            deleteintermediates             , &
            existeNFDE                      , &
            file11isopen                    , &
            NF2FFDecim                      , &
            existeh5                        , &  
            fatalerror                      , & 
            fatalerrornfde2sgg              , & 
            existeconf                      , &
            existecmsh                      , &  
            thereare_stoch                  , &
            experimentalVideal              , &   
            simu_devia                      , &  
            noconformalmapvtk               , &
            createh5filefromsinglebin       , &
            creditosyaprinteados            , &
            use_mtln_wires
      
        integer (kind=4) ::                   &
            wirethickness                    ,&
            inductance_order                 ,&
            finaltimestep                    ,&
            ierr                             ,&
            layoutnumber                     ,&
            size                             ,&
            length                           ,&
            mpidir                           ,&
            flushminutesFields               ,&
            flushminutesData                 ,&
            flushsecondsFields               ,&
            flushsecondsData                 ,&
            forced                           ,&
            maxCPUtime                       ,&
            sgbcdepth                        ,&
            precision                        ,&
            statuse
                                           
        REAL (KIND=RKIND) ::                  &
            maxwireradius                    ,&
            mindistwires                     ,&
            attfactorc                       ,&
            attfactorw                       ,&
            cfltemp                          ,&
            cfl                              ,&
            sgbcfreq                         ,&
            sgbcresol                        ,&
            alphamaxpar                      ,&
            kappamaxpar                      ,&
            alphaOrden                     
                                           
        REAL (KIND=8) ::                      &
             time_begin                      ,&
             time_end                        
        real (kind=RKIND_wires) ::            &
             factorradius                    ,&
             factordelta                     
        
        type (nf2ff_T) ::                           facesNF2FF    
        TYPE (MedioExtra_t) ::                      MEDIOEXTRA    
       type (EpsMuTimeScale_input_parameters_t) :: EpsMuTimeScale_input_parameters
        type (tiempo_t)  ::                         time_out2 
        
!pgi        CHARACTER (LEN=BUFSIZE_LONG) ::   &          
            CHARACTER (LEN=BUFSIZE) ::   &  
             prefix              ,&
             prefixopci          ,&
             prefixopci1         ,&
             opcionespararesumeo ,&
             opcionesoriginales  ,&
             slicesoriginales    ,&
             chdummy  
       
   
        CHARACTER (LEN=BUFSIZE) :: chaininput,     &  
                                   chain,     &
                                   chain2,     &
                                   fichin,         &
                                   extension,         &
                                   nresumeable2,   &
                                   fileFDE,          &
                                   fileH5,           &
                                   inductance_model, &
                                   wiresflavor     , &
                                   nEntradaRoot    , &
                                   opcionestotales , &
                                   ficherohopf,      &
                                   conformal_file_input_name  , &
                                   geomfile         
                        
   end type entrada_t
       
   PUBLIC interpreta,insertalogtmp,print_help,print_basic_help,print_credits, &
       removeintraspaces,buscaswitchficheroinput,default_flags
   PUBLIC entrada_t
   !                                        
CONTAINS
    
   subroutine interpreta(l,statuse)
   
!!!!!!!!!!!!!           
   type (entrada_t), intent(INOUT) :: l
!!!!!!!!!      
   
   CHARACTER (LEN=BUFSIZE) :: chari,f,dubuf,buff
   logical :: existiarunningigual,mpidirset,resume3
   integer (kind=4) :: i,j,donde,n, newmpidir,statuse
   real (KIND=RKIND) :: pausetime

  
!!!
   l%input_conformal_flag=input_conformal_flag !ojooo 051223 es un flag global
   
   mpidirset=.false.
   existiarunningigual=.false.
   statuse=0
   !!!!!!!!!!!!!!!
   n = commandargumentcount (l%chaininput)
   IF (n == 0) THEN
      call print_basic_help(l) 
      call stoponerror(l%layoutnumber,l%size,'Error: NO arguments neither command line nor in launch file. Correct and remove pause...',.true.)
      statuse=-1
      !goto 668
   END IF
   l%opcionestotales=''
   do i=1,n
      CALL getcommandargument (l%chaininput,i,l%chain,l%length,statuse)
      IF (statuse /= 0) THEN
         CALL stoponerror (l%layoutnumber, l%size, 'Reading input',.true.)
          statuse=-1
          !goto 668
      END IF
      l%opcionestotales=trim(adjustl(l%opcionestotales))//' '//trim(adjustl(l%chain))
   end do
   CALL print11 (l%layoutnumber, 'Switches '//trim(adjustl(l%opcionestotales)))

 
   IF (n > 0) THEN
      i = 2  ! se empieza en 2 porque el primer argumento es siempre el nombre del ejecutable
      DO while (i <= n)
         CALL getcommandargument (l%chaininput, i, l%chain, l%length, statuse)
         IF (statuse /= 0) THEN
            CALL stoponerror (l%layoutnumber, l%size, 'Reading input',.true.)
          statuse=-1
          !goto 668
         END IF
         SELECT CASE (trim(adjustl(l%chain)))   
          CASE ('-i')
               i = i + 1
               CALL getcommandargument (l%chaininput, i, f, l%length,  statuse)
             continue !ya interpretado
          case ('-a')
               i = i + 1
               CALL getcommandargument (l%chaininput, i, f, l%length,  statuse)
             continue !ya interpretado
          CASE ('-mpidir')
            i = i + 1
            CALL getcommandargument (l%chaininput, i, f, l%length,  statuse)
            select case (trim (adjustl(f)))
             case ('x','X')
               newmpidir=1  !!!lo cambie por error !161018
             case ('y','Y')
               newmpidir=2   !!!lo cambie por error !161018
             case ('z','Z')
               newmpidir=3
             CASE DEFAULT
               GOTO 1762
            END SELECT 
            if (newmpidir.ne.l%mpidir) then
               GOTO 1762
            endif
            GO TO 2762
1762        CALL stoponerror (l%layoutnumber, l%size, 'Invalid or duplicate incoherent -l%mpidir option',.true.)
          statuse=-1
          goto 668
2762      CONTINUE
          if (.not.mpidirset) then
            l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain)) // ' ' // trim (adjustl(f))
            mpidirset=.true.
          endif
              
          case ('-pause')
            i = i + 1
            CALL getcommandargument (l%chaininput, i, f, l%length,  statuse)
            ! Converts the characters to integer
            READ (f,*, ERR=7312) pausetime
            GO TO 8312
7312        CALL stoponerror (l%layoutnumber, l%size, 'Invalid pause time',.true.)
          statuse=-1
          !goto 668
8312        IF (pausetime <= 0) THEN
               CALL stoponerror (l%layoutnumber, l%size, 'Invalid pause time',.true.)
          statuse=-1
          !goto 668
            END IF
            !
            l%pausar=.true.
#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,l%ierr)
#endif
            CALL get_secnds (l%time_out2)
            l%time_begin = l%time_out2%segundos
            WRITE (dubuf,*) 'Paused for (secs) ', pausetime
            CALL print11 (l%layoutnumber, dubuf)
            DO while (l%pausar)
#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,l%ierr)
#endif
               CALL get_secnds (l%time_out2)
               l%time_end = l%time_out2%segundos
               IF (l%time_end-l%time_begin > pausetime) THEN
                  l%pausar =.false.
               END IF
            END DO
#ifdef CompileWithMPI
            CALL MPI_Barrier (SUBCOMM_MPI, l%ierr)
            l%l_aux = l%pausar
            CALL MPI_AllReduce (l%l_aux, l%pausar, 1_4, MPI_LOGICAL, MPI_LOR, SUBCOMM_MPI, l%ierr)
#endif

            !!!        CASE ('-maxmessages')
            !!!            i = i + 1
            !!!          CALL getcommandargument (l%chaininput, i, f, l%length,  statuse)
            !!!          READ (f,*, ERR=1012) maxmessages
            !!!          GO TO 2012
            !!!1012      CALL stoponerror (l%layoutnumber, l%size, 'Invalid Number of maxmessages',.true.)
          !!!statuse=-1
          !!!!goto 668
            !!!2012      CONTINUE
            !!!          l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain)) // ' ' // trim (adjustl(f))
          CASE ('-NF2FFDecim')
            l%NF2FFDecim=.true.
            l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain)) // ' ' // trim (adjustl(f))
          CASE ('-noNF2FF')
            i = i + 1
            CALL getcommandargument (l%chaininput, i, f, l%length,  statuse)
            select case (trim (adjustl(f)))
             case ('back','BACK')
               l%facesNF2FF%TR=.FALSE.
             case ('front','FRONT')
               l%facesNF2FF%FR=.FALSE.
             case ('left','LEFT')
               l%facesNF2FF%IZ=.FALSE.
             case ('right','RIGHT')
               l%facesNF2FF%DE=.FALSE.
             case ('down','DOWN')
               l%facesNF2FF%AB=.FALSE.
             case ('up','UP')
               l%facesNF2FF%AR=.FALSE.
             CASE DEFAULT
               GOTO 1712
            END SELECT
            GO TO 2712
1712        CALL stoponerror (l%layoutnumber, l%size, 'Invalid -noNF2FF option',.true.)
          statuse=-1
          !goto 668
2712      CONTINUE
            !COMO LA RCS SE CALCULA SOLO AL FINAL NO OBLIGO A RESUMEAR CON IGUAL -NONFF2FF PARA PODER CALCULAR CON Y SIN ESTA OPCION resumeando
            !          l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain)) // ' ' // trim (adjustl(f))      
          CASE ('-force')      
            l%forcing = .TRUE. 
            i = i + 1
            CALL getcommandargument (l%chaininput, i, f, l%length,  statuse)
            READ (f,*, ERR=412) l%forced
            GO TO 312
412         CALL stoponerror (l%layoutnumber, l%size, 'Invalid cut',.true.)
          statuse=-1
          !goto 668
312         CONTINUE
            l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain)) // ' ' // trim (adjustl(f))
          CASE ('-singlefile')
            l%singlefilewrite = .TRUE.
            l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain))  
          CASE ('-ignoresamplingerrors')
            l%ignoresamplingerrors = .TRUE.      
          CASE ('-prioritizeTHINWIRE')
            l%prioritizeTHINWIRE=.true.
            l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain))
            l%ignoreerrors = .TRUE.
          CASE ('-prioritizeCOMPOoverPEC')
            l%prioritizeCOMPOoverPEC=.true.
            l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain))
            l%ignoreerrors = .TRUE.
          CASE ('-noshared')
            l%updateshared=.false.
            l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain))
          CASE ('-prioritizeISOTROPICBODYoverall')
            l%prioritizeISOTROPICBODYoverall=.true.
            l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain))
          CASE ('-wirecrank')
            l%wirecrank = .TRUE.
            l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain))
          CASE ('-clip')
            l%CLIPREGION = .TRUE.
            l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain))
!!!!#endif        
!     
  
          CASE ('-verbose')
            l%verbose = .TRUE.
          CASE ('-ignoreerrors')
            l%ignoreerrors = .TRUE.         
          CASE ('-r')
            l%resume = .TRUE.
            l%forcesteps=.true.
#ifdef CompileWithOldSaving
          CASE ('-old')
            l%resume_fromold = .TRUE.
#endif  
          CASE ('-cpumax')
            i = i + 1
            CALL getcommandargument (l%chaininput, i, f, l%length,  statuse)
            ! Converts the characters to integer
            READ (f,*, ERR=712) l%maxCPUtime
            GO TO 812
712         CALL stoponerror (l%layoutnumber, l%size, 'Invalid CPU maximum time',.true.)
          statuse=-1
          !goto 668
812         IF (l%maxCPUtime <= 0) THEN
               CALL stoponerror (l%layoutnumber, l%size, 'Invalid CPU maximum time',.true.)
          statuse=-1
          !goto 668
            END IF   

          CASE ('-s')
            l%freshstart = .TRUE.
          CASE ('-flush')
            i = i + 1
            CALL getcommandargument (l%chaininput, i, f, l%length,  statuse)
            ! Converts the characters to integer
            READ (f,*, ERR=300) l%flushminutesFields
            GO TO 400
300         CALL stoponerror (l%layoutnumber, l%size, 'Invalid flushing interval',.true.)
          statuse=-1
          !goto 668
400         IF (l%flushminutesFields <= 0) THEN
               CALL stoponerror (l%layoutnumber, l%size, 'Invalid flushing interval',.true.)
          statuse=-1
          !goto 668
            END IF
          CASE ('-flushdata')
            i = i + 1
            CALL getcommandargument (l%chaininput, i, f, l%length,  statuse)
            ! Converts the characters to integer
            READ (f,*, ERR=301) l%flushminutesData
            GO TO 401
301         CALL stoponerror (l%layoutnumber, l%size, 'Invalid flushing interval',.true.)
          statuse=-1
          !goto 668
401         IF (l%flushminutesData <= 0) THEN
               CALL stoponerror (l%layoutnumber, l%size, 'Invalid flushing interval',.true.)
          statuse=-1
          !goto 668
            END IF   
          CASE ('-run')
            l%run = .TRUE.                
          CASE ('-map')
            !dump the map files
            l%createmap = .TRUE.
          CASE ('-dontwritevtk')
            l%dontwritevtk=.true.
          CASE ('-vtkindex')
            l%vtkindex = .TRUE.  
          CASE ('-mapvtk')
            l%createmapvtk = .TRUE.
          CASE ('-hopf')
            l%hopf=.true.
            i = i + 1;
            l%ficherohopf = char(0);
            CALL getcommandargument (l%chaininput, i, f, l%length,  statuse)
            l%ficherohopf = trim(adjustl(f));
            INQUIRE (file=trim(adjustl(f)), EXIST=l%existeNFDE)
            IF ( .NOT. l%existeNFDE) THEN
               l%hopf = .FALSE.;
               buff = 'The l%hopf input file was not found '//trim(adjustl(l%ficherohopf));
               call WarnErrReport(Trim(buff),.true.)
            END IF
            l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain)) // ' ' // trim (adjustl(f))
!
          CASE ('-dmma')
              l%run_with_dmma = .TRUE.
              l%run_with_abrezanjas = .FALSE.
              l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain))
!!              i = i + 1;

#ifdef CompileWithConformal
          CASE ('-abrezanjas') !Provisional FEB-2018

            !NOTE: Comento lo de abajo para debugear 
            !   print *,'Abrezanjas not available (290521).  '
            !   stop

            l%run_with_dmma = .FALSE.
            l%run_with_abrezanjas = .true.
            if (.NOT.l%input_conformal_flag) then
                l%conformal_file_input_name = char(0)
                l%input_conformal_flag = .true. 
            end if
            l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain))

!!            i = i + 1;
          CASE ('-activateconf') !Provisional FEB-2018
            if (.NOT.l%input_conformal_flag) then
                l%conformal_file_input_name = char(0)
                l%input_conformal_flag = .true.
            end if
            i = i + 1;      
            l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain))      
#endif            
          CASE ('-conf')
            l%flag_conf_sgg=.true.
            i = i + 1;
#ifdef CompileWithConformal              
            l%input_conformal_flag = .true.;
            l%conformal_file_input_name = char(0);

            CALL getcommandargument (l%chaininput, i, f, l%length,  statuse)
            l%conformal_file_input_name = trim(adjustl(f));

            INQUIRE (file=trim(adjustl(f)), EXIST=l%existeNFDE)
            IF ( .NOT. l%existeNFDE) THEN
               l%input_conformal_flag = .FALSE.;
               buff = 'The conformal input file was not found '//trim(adjustl(l%fichin));
               call WarnErrReport(Trim(buff),.true.)
            END IF
            l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain))
#endif
#ifndef CompileWithConformal
            IF (l%input_conformal_flag)THEN
               buff='';buff='Conformal without conformal support. Recompile!';
               call WarnErrReport(Trim(buff),.true.)
            END IF
#endif

          CASE ('-takeintcripte')
            l%takeintcripte=.true.
            l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain))
#ifdef CompileWithNIBC
          CASE ('-skindepthpre')
            l%skindepthpre=.true.
!            l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain))
          CASE ('-mibc','-skindepth')
            l%mibc=.true.
            l%sgbc=.false.
            l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain))
          CASE ('-conformalskin')
            l%conformalskin=.true.
            l%mibc=.true.
            l%sgbc=.false.
            l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain))
          CASE ('-ade')
            l%ade=.true.
            l%mibc=.true.
            l%sgbc=.false.
            l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain))
          CASE ('-NOcompomur')
            l%NOcompomur=.true.
            l%mibc=.true.
            l%sgbc=.false.
            l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain))
#endif
          CASE ('-mur2')
            l%MurAfterPML=.true.
            !l%mur_second=.true.
            l%mur_first=.true.
            !arreglar cuando resuelva el bug en mur segundo orden
            l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain))
          CASE ('-mur1')
            l%MurAfterPML=.true.
            l%mur_first=.true.
            l%mur_second=.false.
            l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain))
          CASE ('-pmlalpha')
            i = i + 1
            CALL getcommandargument (l%chaininput, i, f, l%length,  statuse)
            ! Converts the characters to real
            READ (f,*, ERR=7621) l%alphamaxpar
            GO TO 8621
7621        CALL stoponerror (l%layoutnumber, l%size, 'Invalid CPML alpha factor',.true.)
          statuse=-1
          !goto 668
8621        IF (l%alphamaxpar < 0.0_RKIND) THEN
               CALL stoponerror (l%layoutnumber, l%size, 'Invalid CPML alpha factor',.true.)
          statuse=-1
          !goto 668
            END IF
            i = i + 1
            !          l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain))// ' ' // trim (adjustl(f))
            CALL getcommandargument (l%chaininput, i, f, l%length,  statuse)
            ! Converts the characters to real
            READ (f,*, ERR=7121) l%alphaOrden
            GO TO 8121
7121        CALL stoponerror (l%layoutnumber, l%size, 'Invalid CPML order factor',.true.)
          statuse=-1
          !goto 668
8121        IF (l%alphaOrden < 0.0_RKIND) THEN
               CALL stoponerror (l%layoutnumber, l%size, 'Invalid CPML alpha factor',.true.)
          statuse=-1
          !goto 668
            END IF
            !          l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain))// ' ' // trim (adjustl(f))
          CASE ('-pmlkappa')
            i = i + 1
            CALL getcommandargument (l%chaininput, i, f, l%length,  statuse)
            ! Converts the characters to real
            READ (f,*, ERR=7622) l%kappamaxpar
            GO TO 8622
7622        CALL stoponerror (l%layoutnumber, l%size, 'Invalid CPML kappa factor',.true.)
          statuse=-1
          !goto 668
8622        IF (l%kappamaxpar < 1.0_RKIND) THEN
               CALL stoponerror (l%layoutnumber, l%size, 'Invalid CPML kappa factor',.true.)
          statuse=-1
          !goto 668
            END IF
            !          l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain))// ' ' // trim (adjustl(f))
          CASE ('-pmlcorr')
            l%MEDIOEXTRA%exists=.true.
            i = i + 1
            CALL getcommandargument (l%chaininput, i, f, l%length,  statuse)
            ! Converts the characters to real
            READ (f,*, ERR=7672) l%MEDIOEXTRA%sigma
            GO TO 8672
7672        CALL stoponerror (l%layoutnumber, l%size, 'Invalid pmlcorr sigma factor',.true.)
          statuse=-1
          !goto 668
8672        IF (l%MEDIOEXTRA%sigma < 0.0_RKIND) THEN
               CALL stoponerror (l%layoutnumber, l%size, 'Invalid pmlcorr sigma factor',.true.)
          statuse=-1
          !goto 668
            END IF
            l%MEDIOEXTRA%sigmam=-1.0_RKIND!voids it. later overriden
            !          l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain))// ' ' // trim (adjustl(f))
            i = i + 1
            CALL getcommandargument (l%chaininput, i, f, l%length,  statuse)
            ! Converts the characters to real
            READ (f,*, ERR=7662) l%MEDIOEXTRA%size
            GO TO 8662
7662        CALL stoponerror (l%layoutnumber, l%size, 'Invalid pmlcorr depth factor',.true.)
          statuse=-1
          !goto 668
8662        IF (l%MEDIOEXTRA%size < 0) THEN
               CALL stoponerror (l%layoutnumber, l%size, 'Invalid pmlcorr depth factor',.true.); statuse=-1; !goto 668
            END IF
            !          l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain))// ' ' // trim (adjustl(f))
          CASE ('-attc')
            i = i + 1
            CALL getcommandargument (l%chaininput, i, f, l%length,  statuse)
            ! Converts the characters to real
            READ (f,*, ERR=766) l%attfactorc
            GO TO 866
766         CALL stoponerror (l%layoutnumber, l%size, 'Invalid dissipation factor',.true.); statuse=-1; !goto 668
866         IF ((l%attfactorc <= -1.0_RKIND ).or.(l%attfactorc > 1.0_RKIND)) THEN
               CALL stoponerror (l%layoutnumber, l%size, 'Invalid dissipation factor',.true.); statuse=-1; !goto 668
            END IF
            l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain))// ' ' // trim (adjustl(f))
          CASE ('-sgbcdepth')
            l%mibc=.false.
            l%sgbc=.true.
            i = i + 1
            CALL getcommandargument (l%chaininput, i, f, l%length,  statuse)
            ! Converts the characters to real
            READ (f,*, ERR=7466) l%sgbcdepth
            GO TO 8466
7466        CALL stoponerror (l%layoutnumber, l%size, 'Invalid sgbc depth ',.true.); statuse=-1; !goto 668
8466        IF (l%sgbcdepth < -1 ) THEN
               CALL stoponerror (l%layoutnumber, l%size, 'Invalid sgbc depth',.true.); statuse=-1; !goto 668
            END IF
            l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain))// ' ' // trim (adjustl(f))
          CASE ('-sgbcfreq')
            l%sgbc=.true.
            l%mibc=.false.
            i = i + 1
            CALL getcommandargument (l%chaininput, i, f, l%length,  statuse)
            ! Converts the characters to real
            READ (f,*, ERR=74616) l%sgbcfreq
            GO TO 84616
74616       CALL stoponerror (l%layoutnumber, l%size, 'Invalid sgbc freq ',.true.); statuse=-1; !goto 668
84616       IF (l%sgbcfreq < 0. ) THEN
            CALL stoponerror (l%layoutnumber, l%size, 'Invalid sgbc freq',.true.); statuse=-1; !goto 668
            END IF
            l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain))// ' ' // trim (adjustl(f))
          CASE ('-sgbcresol')
            l%mibc=.false.
            l%sgbc=.true.
            i = i + 1
            CALL getcommandargument (l%chaininput, i, f, l%length,  statuse)
            ! Converts the characters to real
            READ (f,*, ERR=74626) l%sgbcresol
            GO TO 84626
74626       CALL stoponerror (l%layoutnumber, l%size, 'Invalid sgbc decay ',.true.); statuse=-1; !goto 668
84626       IF (l%sgbcresol < 0.0 ) THEN
            CALL stoponerror (l%layoutnumber, l%size, 'Invalid sgbc decay',.true.); statuse=-1; !goto 668
            END IF
            l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain))// ' ' // trim (adjustl(f))
          CASE ('-sgbcyee')
            l%sgbc=.true.
            l%mibc=.false.
            l%sgbccrank=.false.
            l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain))
          CASE ('-sgbccrank') !es el default. Lo mantengo por compatibilidad con lanzamientos previos
            l%sgbccrank=.true.
            l%mibc=.false.
            l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain))
          CASE ('-nosgbc') !opcion generica que aglutina varios switches que estan den default (l%sgbcresol, l%sgbccrank, l%sgbcfreq)
            l%sgbc=.false.
            l%mibc=.true.
            l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain))
          CASE ('-sgbc') !opcion generica que aglutina varios switches que estan den default (l%sgbcresol, l%sgbccrank, l%sgbcfreq)
            l%sgbc=.true.
            l%mibc=.false.
            l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain))
          CASE ('-sgbcDispersive') !opcion generica que aglutina varios switches que estan den default (l%sgbcresol, l%sgbccrank, l%sgbcfreq)
            l%sgbc=.true.
            l%mibc=.false.
            l%sgbcDispersive=.true.
            l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain))
          CASE ('-saveall')
            l%saveall = .TRUE.
          CASE ('-attw')
            i = i + 1
            CALL getcommandargument (l%chaininput, i, f, l%length,  statuse)
            ! Converts the characters to real
            READ (f,*, ERR=732) l%attfactorw
            GO TO 832
732         CALL stoponerror (l%layoutnumber, l%size, 'Invalid dissipation factor',.true.); statuse=-1; !goto 668
832         IF ((l%attfactorw <= -1.0_RKIND ).or.(l%attfactorw > 1.0_RKIND)) THEN
               CALL stoponerror (l%layoutnumber, l%size, 'Invalid dissipation factor',.true.); statuse=-1; !goto 668
            END IF
            l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain))// ' ' // trim (adjustl(f))
          CASE ('-maxwireradius')
            l%boundwireradius=.true.
            i = i + 1
            CALL getcommandargument (l%chaininput, i, f, l%length,  statuse)
            ! Converts the characters to real
            READ (f,*, ERR=737) l%maxwireradius
            GO TO 837
737         CALL stoponerror (l%layoutnumber, l%size, 'Invalid dissipation factor',.true.); statuse=-1; !goto 668
837         IF ((l%maxwireradius <= 0.0_RKIND )) THEN
               CALL stoponerror (l%layoutnumber, l%size, 'Invalid maximumwireradius',.true.); statuse=-1; !goto 668
            END IF
            l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain))// ' ' // trim (adjustl(f))
          CASE ('-mindistwires')
            i = i + 1
            CALL getcommandargument (l%chaininput, i, f, l%length,  statuse)
            ! Converts the characters to real
            READ (f,*, ERR=1732) l%mindistwires
            GO TO 1832
1732        CALL stoponerror (l%layoutnumber, l%size, 'Invalid minimum distance between wires',.true.); statuse=-1; !goto 668
1832        IF (l%mindistwires <= 0.0_RKIND ) THEN
               CALL stoponerror (l%layoutnumber, l%size, 'Invalid minimum distance between wires',.true.); statuse=-1; !goto 668
            END IF
            l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain))// ' ' // trim (adjustl(f))
          CASE ('-makeholes')
            l%makeholes = .TRUE.
            l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain))
          CASE ('-connectendings')
            l%connectendings = .TRUE.
            l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain))
          CASE ('-nostrictOLD')
            l%strictOLD = .false.
            l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain))
          CASE ('-nomtlnberenger')
            l%mtlnberenger = .false.
            l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain))
          CASE ('-stableradholland')
            l%stableradholland = .true.
            l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain))
         !  CASE ('-mtln')
         !      buff='-mtln option deprecated and ignored. Check -nomtlnberenger or -l%stableradholland'
         !      call WarnErrReport(Trim(buff),.false.)
          CASE ('-intrawiresimplify')
            l%strictOLD = .false.
            l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain))
          CASE ('-notaparrabos')
            l%TAPARRABOS = .false.
            l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain))
          case ('-fieldtotl')
              l%fieldtotl=.true.
            l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain))
          !!case ('-experimentalVideal')
          !!    l%experimentalVideal=.true.
          !!    l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain))

          case ('-forceresampled') !a menos que se pida explicitamente, no se resamplea 120123
              l%forceresampled=.true.
              l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain))   
              
          case ('-mtlnwires')
            l%use_mtln_wires = .true.
          CASE ('-wirethickness')
            i = i + 1
            CALL getcommandargument (l%chaininput, i, f, l%length,  statuse)
            ! Converts the characters to real
            READ (f,*, ERR=7416) l%wirethickness
            GO TO 8416
7416        CALL stoponerror (l%layoutnumber, l%size, 'Invalid l%wirethickness ',.true.); statuse=-1; !goto 668
8416        IF (l%sgbcdepth < -1 ) THEN
               CALL stoponerror (l%layoutnumber, l%size, 'Invalid l%wirethickness',.true.); statuse=-1; !goto 668
            END IF
            l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain))// ' ' // trim (adjustl(f))              
          CASE ('-wiresflavor')
            i = i + 1
            CALL getcommandargument (l%chaininput, i, f, l%length,  statuse)
            l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain)) // ' ' // trim (adjustl(f))
            READ (f, '(a)', ERR=3621) l%wiresflavor
            if (trim(adjustl(l%wiresflavor(1:1)))=='g') l%wiresflavor='slanted' 
            select case (trim(adjustl(l%wiresflavor)))
            case ('holland','old')
                l%wiresflavor='holland'
            case ('berenger','new')
                l%wiresflavor='berenger'
            case ('slanted','experimental')
                l%wiresflavor='slanted'
            case ('transition')
                l%wiresflavor='transition'
            case ('semistructured')
                l%wiresflavor='semistructured'
                !
                i = i + 1
                CALL getcommandargument (l%chaininput, i, f, l%length,  statuse)
                l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(f))
            ! Converts the characters to real
                READ (f,*, ERR=2561) l%precision
                GO TO 2562
2561            CALL stoponerror (l%layoutnumber, l%size, 'Invalid l%precision for semistructured',.true.); statuse=-1; !goto 668
2562            IF (l%precision < 0 ) THEN
                    CALL stoponerror (l%layoutnumber, l%size, 'Invalid l%precision for semistructured',.true.); statuse=-1; !goto 668
                END IF
            !   
            end select   
            GO TO 4621
3621        CALL stoponerror (l%layoutnumber, l%size, 'Invalid wires flavor',.true.); statuse=-1; !goto 668
4621        IF ( ((trim(adjustl(l%wiresflavor)) /= 'holland')  .AND. &
                  (trim(adjustl(l%wiresflavor)) /= 'transition') .AND. &
                  (trim(adjustl(l%wiresflavor)) /= 'berenger') .AND.  &
                  (trim(adjustl(l%wiresflavor)) /= 'slanted').and. &
                  (trim(adjustl(l%wiresflavor))/='semistructured')) .or. &
             .not.((trim(adjustl(l%wiresflavor)) == 'holland')  .xor.  &
                  (trim(adjustl(l%wiresflavor)) == 'transition') .xor. &
                  (trim(adjustl(l%wiresflavor)) == 'berenger') .xor.  &
                  (trim(adjustl(l%wiresflavor)) == 'slanted').xor.  &
                  (trim(adjustl(l%wiresflavor)) =='semistructured')) )  THEN
               CALL stoponerror (l%layoutnumber, l%size, 'Invalid wires flavor->'//trim(adjustl(l%wiresflavor)),.true.); statuse=-1; !goto 668
            END IF    
#ifndef CompileWithThickWires
            select case (trim(adjustl(l%wiresflavor)))
            case ('holland','transition')
                if (l%wirethickness/=1) then
                    CALL stoponerror (l%layoutnumber, l%size, 'Holland wire flavor not available in this compilation',.true.); statuse=-1; !goto 668
                endif
            end select   
#endif    
#ifndef CompileWithThickWires
            select case (trim(adjustl(l%wiresflavor)))
            case ('holland')
                if (l%wirethickness/=1) then
                    CALL stoponerror (l%layoutnumber, l%size, 'Holland wire flavor thickness>1 requires recompiling',.true.); statuse=-1; !goto 668
                endif
            end select   
#endif
            select case (trim(adjustl(l%wiresflavor)))
            case ('berenger','slanted','experimental','transition')   
                if (l%wirethickness/=1) then
                    CALL stoponerror (l%layoutnumber, l%size, 'Thickness>1 unsupported for this wireflavor',.true.); statuse=-1; !goto 668
                endif    
            end select   
#ifndef CompileWithBerengerWires
            select case (trim(adjustl(l%wiresflavor)))
            case ('berenger')
                CALL stoponerror (l%layoutnumber, l%size, 'Berenger wire flavor not available in this compilation',.true.); statuse=-1; !goto 668
            end select   
#endif
#ifndef CompileWithSlantedWires
            select case (trim(adjustl(l%wiresflavor)))
            case ('slanted','experimental')
                CALL stoponerror (l%layoutnumber, l%size, 'Experimental wire flavor not available in this compilation',.true.); statuse=-1; !goto 668
            end select   
#endif
          CASE ('-isolategroupgroups')
            l%isolategroupgroups = .TRUE.
            l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain))
          CASE ('-groundwires')
            l%groundwires = .TRUE.
            l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain))
          CASE ('-noSlantedcrecepelo ') !opcion niapa excperimental 131219
            l%noSlantedcrecepelo  = .TRUE.
            l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain))
          CASE ('-inductance')
            i = i + 1
            CALL getcommandargument (l%chaininput, i, f, l%length,  statuse)
            READ (f, '(a)', ERR=361) l%inductance_model
            GO TO 461
361         CALL stoponerror (l%layoutnumber, l%size, 'Invalid inductance model',.true.); statuse=-1; !goto 668
461         IF ((l%inductance_model /= 'ledfelt') .AND. (l%inductance_model /= 'berenger') .AND. &
            &    (l%inductance_model /= 'boutayeb')) THEN
               CALL stoponerror (l%layoutnumber, l%size, 'Invalid inductance model',.true.); statuse=-1; !goto 668
            END IF
            l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain)) // ' ' // trim (adjustl(f))
          CASE ('-inductanceorder')
            i = i + 1
            CALL getcommandargument (l%chaininput, i, f, l%length,  statuse)
            READ (f,*, ERR=179) l%inductance_order
            GO TO 180
179         CALL stoponerror (l%layoutnumber, l%size, 'Invalid inductance order',.true.); statuse=-1; !goto 668
180         l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain)) // ' ' // trim (adjustl(f))
          CASE ('-prefix')
            i = i + 1
            CALL getcommandargument (l%chaininput, i, f, l%length,  statuse)
            l%prefix = '_' // trim (adjustl(f))
            l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain)) // ' ' // trim (adjustl(f))
          CASE ('-cfl')
            i = i + 1
            CALL getcommandargument (l%chaininput, i, f, l%length,  statuse)
            ! Converts the characters to real
            READ (f,*, ERR=3762) l%cfltemp
            GO TO 3862
3762        CALL stoponerror (l%layoutnumber, l%size, 'Invalid Courant Number',.true.); statuse=-1; !goto 668
3862        IF (l%cfltemp <= 0.0 ) THEN
               call print11(l%layoutnumber,'------> Ignoring negative or null l%cfl Courant Number')
!!!!!!!!!               CALL stoponerror (l%layoutnumber, l%size, 'Invalid negative or null l%cfl Courant Number',.true.); statuse=-1; !goto 668 !!!sgg 151216 para evitar el error l%cfl 0 del problem-type sigue como si no estuviera
               l%forcecfl=.false.
            else
               l%cfl=l%cfltemp
               l%forcecfl=.true.
               l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain))// ' ' // trim (adjustl(f))
            END IF           
          CASE ('-noconformalmapvtk')
            l%noconformalmapvtk=.true.
          CASE ('-niapapostprocess')
            l%niapapostprocess=.true.
#ifdef CompileWithPrescale
!!!!210918 permit scaling
          CASE ('-pscale')
            l%permitscaling=.true.
            l%saveall=.true. !lo salvo todo en permit scaling para evitar errores
            i = i + 1
            buff=""
            CALL getcommandargument (l%chaininput, i, f, l%length,  statuse)
            READ (f,*, ERR=33762) buff
            l%EpsMuTimeScale_input_parameters%electric=.False.
            l%EpsMuTimeScale_input_parameters%electric=.False.
            Select case (trim(adjustl(buff)))
            case("ee")
               l%EpsMuTimeScale_input_parameters%electric=.True.
            case("hh")
               l%EpsMuTimeScale_input_parameters%magnetic=.True.
            case("eh": "he")
               l%EpsMuTimeScale_input_parameters%electric=.True.
               l%EpsMuTimeScale_input_parameters%magnetic=.True.
            case default
               GO TO 33862
            end select
            i = i + 1
            CALL getcommandargument (l%chaininput, i, f, l%length,  statuse)
            l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain))// ' ' // trim (adjustl(f))
            ! Converts the characters to real
            READ (f,*, ERR=33762) l%EpsMuTimeScale_input_parameters%tini
            i = i + 1
            CALL getcommandargument (l%chaininput, i, f, l%length,  statuse)
            l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(f))
            READ (f,*, ERR=33762) l%EpsMuTimeScale_input_parameters%tend
            i = i + 1
            CALL getcommandargument (l%chaininput, i, f, l%length,  statuse)
            l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(f))
            READ (f,*, ERR=33762) l%EpsMuTimeScale_input_parameters%alpha_max
            GO TO 33862
33762       CALL stoponerror (l%layoutnumber, l%size, 'Invalid pscale parameters',.true.); statuse=-1; !goto 668
33862       continue
            IF (l%EpsMuTimeScale_input_parameters%checkError()/=0) THEN
                CALL stoponerror (l%layoutnumber, l%size, &
   &'Invalid -pscale parameters: some parameters have to be greater than 0.0: -pscale t0(>=0) tend slope(>0)'&
                  &,.true.); statuse=-1; !goto 668
            else
               l%EpsMuTimeScale_input_parameters%are_there = .true.
            endif
#endif       
          CASE ('-n')
            l%forcesteps = .TRUE.
            i = i + 1
            CALL getcommandargument (l%chaininput, i, f, l%length,  statuse)
            ! Converts the characters to integer
            READ (f,*, ERR=602) l%finaltimestep
            GO TO 702
602         CALL stoponerror (l%layoutnumber, l%size, 'Invalid time step',.true.); statuse=-1; !goto 668
702         IF (l%finaltimestep < -2) THEN
               CALL stoponerror (l%layoutnumber, l%size, 'Invalid time step',.true.); statuse=-1; !goto 668
            END IF      
!!!!!!     
          CASE ('-factorradius')
            i = i + 1
            CALL getcommandargument (l%chaininput, i, f, l%length,  statuse)
            ! Converts the characters to integer
            READ (f,*, ERR=6032) l%factorradius
            GO TO 7032
6032         CALL stoponerror (l%layoutnumber, l%size, 'Invalid l%factorradius',.true.); statuse=-1; !goto 668
7032         continue
          CASE ('-factordelta')
            i = i + 1
            CALL getcommandargument (l%chaininput, i, f, l%length,  statuse)
            ! Converts the characters to integer
            READ (f,*, ERR=6072) l%factordelta
            GO TO 7072
6072         CALL stoponerror (l%layoutnumber, l%size, 'Invalid l%factordelta',.true.); statuse=-1; !goto 668
7072         continue
!!!!!!!!!!!!!
          CASE ('-stoch')
            l%stochastic=.true.
            l%chosenyesornostochastic=.true.
            l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain))
#ifndef CompileWithMPI
            CALL stoponerror (l%layoutnumber, l%size, 'l%stochastic simulation unsupported without MPI compilation',.true.); statuse=-1; !goto 668
#endif
          CASE ('-nostoch')
            l%stochastic=.false.
            l%chosenyesornostochastic=.true.
            l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain))         
          case ('-forcecreateh5bin')         
            l%createh5bin=.true.     
          CASE ('') !100615 para evitar el crlf del .sh
            continue
          CASE DEFAULT
            CALL stoponerror (l%layoutnumber, l%size, 'Wrong switch '//trim(adjustl(l%chain)),.true.); statuse=-1; !goto 668
         END SELECT
         i = i + 1
      END DO

   END IF
   !some checkings
   !just to be sure that I do not have stupid errors
#ifdef CompileWithMPI
   IF ((sizeof(MPI_DOUBLE_PRECISION) /= 4) .OR. (sizeof(MPI_REAL) /= 4)) THEN
      CALL stoponerror (l%layoutnumber, l%size, 'SEVERE COMPILATION ERROR: MPI_REAL l%size is not 4. ')
   END IF
#endif

   IF (l%connectendings .AND. l%strictOLD) THEN
      CALL stoponerror (l%layoutnumber, l%size, 'l%strictOLD option not compatible with -l%connectendings',.true.); statuse=-1; !goto 668
   END IF
   IF (l%TAPARRABOS .AND. (.not.l%strictOLD)) THEN
      CALL stoponerror (l%layoutnumber, l%size, '-nostrictOLD option requires -notaparrabos ',.true.); statuse=-1; !goto 668
   END IF
   IF (l%isolategroupgroups .AND. l%strictOLD) THEN
      CALL stoponerror (l%layoutnumber, l%size, '-intrawiresimplify option not compatible with -l%isolategroupgroups',.true.); statuse=-1; !goto 668
   END IF

   IF ((l%sgbc .AND. l%mibc)) THEN
      CALL stoponerror (l%layoutnumber, l%size, 'Use only one of -sgbc -l%mibc',.true.); statuse=-1; !goto 668
   END IF
   IF (l%freshstart .AND. l%resume) THEN
      CALL stoponerror (l%layoutnumber, l%size, 'Fresh Start option -s not compatible with restarting -r',.true.); statuse=-1; !goto 668
   END IF
   IF (l%freshstart .AND. l%resume_fromold) THEN
      CALL stoponerror (l%layoutnumber, l%size, 'Fresh Start option -s not compatible with -old',.true.); statuse=-1; !goto 668
   END IF
   IF (( .NOT. l%resume).and.(.not.l%run) .AND. l%resume_fromold) THEN
      CALL stoponerror (l%layoutnumber, l%size, 'l%resume option -r must be used if issuing -old',.true.); statuse=-1; !goto 668
   END IF
   IF ((l%flushminutesFields /= 0) .AND. (l%deleteintermediates)) THEN
      CALL stoponerror (l%layoutnumber, l%size, '-delete is not compatible with -flush',.true.); statuse=-1; !goto 668
   END IF
   if (l%run_with_abrezanjas.and.l%run_with_dmma) then
      CALL stoponerror (l%layoutnumber, l%size, '-abrezanjas is not compatible with -dmma',.true.); statuse=-1; !goto 668
   END IF
   if (l%stochastic.and.(trim(adjustl(l%wiresflavor))/='holland')) then
      CALL stoponerror (l%layoutnumber, l%size, 'Old wires flavor is the only supported with l%stochastic',.true.); statuse=-1; !goto 668
   END IF
   if (l%stochastic.and.l%wirecrank) then
      CALL stoponerror (l%layoutnumber, l%size, 'Wires Crank Nicolson is unsupported with l%stochastic',.true.); statuse=-1; !goto 668
   END IF
   !!!si esta soportado 170719
   !! if (l%permitscaling.and.l%resume) then
   !!   CALL stoponerror (l%layoutnumber, l%size, 'Resuming with Permittivity scaling unsupported',.true.); statuse=-1; !goto 668
   !!END IF
   if (l%permitscaling.and.(l%kappamaxpar.gt.1.000001_rkind)) then
   !!!061118 no lo permito porque cpml toca los idxe, idye, idze en funcion del kappa y permittivity scaling conflicta
            CALL stoponerror (l%layoutnumber, l%size, 'Unsupported CPML kappa factor since 061118 because conflicts with Idxe...in permittivity scaling',.true.)
   endif
   if (l%stochastic) then
#ifndef CompileWithStochastic
        call StopOnError(l%layoutnumber,l%size,'l%stochastic without compilation support. Recompile')
#endif
#ifdef CompileWithStochastic
#ifndef CompileWithMPI
        call StopOnError(l%layoutnumber,l%size,'l%stochastic unsupported without MPI compilation. Recompile')
#endif
#endif 
        continue
   endif
   
!!!
       
   !
   l%prefixopci1=trim (adjustl(l%opcionespararesumeo))
   l%prefixopci=' '
   do i=1,len(trim (adjustl(l%prefixopci1)))
      l%prefixopci(i:i)=l%prefixopci1(i:i)
      j=ichar(l%prefixopci1(i:i))
      if (j <= 47) l%prefixopci(i:i)='_'
      if (j >= 123) l%prefixopci(i:i)='_'
      if ((j >= 58).and.(j <= 64)) l%prefixopci(i:i)='_'
      if ((j >= 91).and.(j <= 96)) l%prefixopci(i:i)='_'
      if (j == 46) l%prefixopci(i:i)='p'
   end do

   do i=1,len(trim (adjustl(l%prefixopci)))
      do while (l%prefixopci(i:i+1) =='__')
         l%prefixopci(i:)=l%prefixopci(i+1:)
      end do
   end do
   if (l%prefix(1:1)=='_') then
     !!!acortado 120219  l%nEntradaRoot = trim (adjustl(l%fichin)) // trim (adjustl(prefix))// trim (adjustl(l%prefixopci))
      l%nEntradaRoot = trim (adjustl(l%fichin)) //'_'// trim (adjustl(l%prefixopci))
   else
      l%nEntradaRoot = trim (adjustl(l%fichin))
   endif
!!!l%stochastic
#ifdef CompileWithStochastic
  if (l%stochastic) then  
    if (l%layoutnumber <= l%size/2-1) then !aun no se ha dividido el l%size
          l%nEntradaRoot=trim (adjustl(l%nEntradaRoot))
    else
          l%nEntradaRoot=trim (adjustl('devia_'//trim (adjustl(l%nEntradaRoot))))
    endif
  endif
#ifdef CompileWithMPI
  CALL MPI_Barrier (SUBCOMM_MPI, l%ierr)
#endif
#endif
!!!fin l%stochastic
!!!   sgg%nEntradaRoot=trim (adjustl(l%nEntradaRoot))
   !
   WRITE (chari, '(i5)') l%layoutnumber + 1
   l%nresumeable2 = trim (adjustl(l%nEntradaRoot)) // '_' // trim (adjustl(chari)) // '.fields'
   !

   l%geomfile = trim (adjustl(l%nEntradaRoot)) // '_' // trim (adjustl(chari))
   !warning file management
   if (statuse/=-1) then
       CALL CLOSEWARNINGFILE(l%layoutnumber,l%size,l%fatalerror,.false.,.false.) !!cierra el temporal !todavia no se ha dividido el l%size
       !!!borra el tmp si no hay l%fatalerror y reabre el de verdad
       if ((.not.l%fatalerror).and.(l%layoutnumber==0)) then 
           open(unit=1320,file=trim (adjustl(l%fichin))//'_tmpWarnings.txt_Warnings.txt')
           close(unit=1320,status='delete')
       endif
       CALL INITWARNINGFILE (l%layoutnumber, l%size, l%nEntradaRoot,l%verbose,l%ignoreerrors)
   endif
   !
   !
   IF (l%resume_fromold) THEN
      INQUIRE (file=trim(adjustl(l%nresumeable2))//'.old', EXIST=resume3)
   ELSE
      INQUIRE (file=trim(adjustl(l%nresumeable2)), EXIST=resume3)
   END IF
   IF (l%resume) THEN
      IF ( .NOT. resume3) THEN
         CALL stoponerror (l%layoutnumber, l%size, 'l%resume fields not present',.true.); statuse=-1; !goto 668
      END IF
      WRITE (dubuf,*) 'RESUMING simulation ', trim (adjustl(l%nEntradaRoot)), ' until n= ', l%finaltimestep
      CALL print11 (l%layoutnumber, dubuf)
   ELSE
      IF (resume3 .AND. ( .NOT. l%freshstart).and.(.not.l%run)) THEN
         CALL stoponerror (l%layoutnumber, l%size, 'Restarting file exists. Either specify -r to l%resume, -s to do a fresh START, or -run to run in whatever the case',.true.); statuse=-1; !goto 668
      ELSEIF (resume3.and.(l%run)) THEN
         l%resume=.true.
      ELSE
         OPEN (35, file=trim(adjustl(l%nresumeable2)))
         WRITE (35, '(a)') '!END'
         CLOSE (35, status='DELETE')
         OPEN (35, file=trim(adjustl(l%nresumeable2))//'.old')
         WRITE (35, '(a)') '!END'
         CLOSE (35, status='DELETE')
      END IF
   END IF
!
!
!
   if (((l%wiresflavor=='slanted').or.(l%wiresflavor=='semistructured')).AND.(l%mpidir/=3)) then
       continue !arreglado l%mpidir slanted 2019
       !         CALL stoponerror (l%layoutnumber, l%size, 'slanted wires unsupported with -l%mpidir {x,y}',.true.); statuse=-1; !goto 668
   endif
   if (l%input_conformal_flag.AND.(l%mpidir/=3)) then
        continue !arreglado l%mpidir conformal 2019
         !TODO: under test
         !26-sep-2018: lo comento 
         !CALL stoponerror (l%layoutnumber, l%size, 'CONFORMAL -conf  unsupported with -l%mpidir {x,y}',.true.); statuse=-1; !goto 668
   endif
   if (l%run_with_abrezanjas.AND.(l%mpidir/=3)) then
        continue !arreglado l%mpidir conformal 2019
         !under test
         !26-sep-2018: lo comento 
         !CALL stoponerror (l%layoutnumber, l%size, 'New abrezanjas thin gaps unsupported with -l%mpidir {x,y}',.true.); statuse=-1; !goto 668
   endif
   if (l%run_with_abrezanjas.AND.l%flag_conf_sgg) then
      !pass Mayo-2018
         !CALL stoponerror (l%layoutnumber, l%size, 'CONFORMAL -conf currently unsupported with new abrezanjas thin gaps (unsupported 2 simultaneous conformal meshes at this moment',.true.); statuse=-1; !goto 668
         !se hace en otro sitio
   endif


   !
   IF (((l%forcesteps) .AND. ( .NOT. l%freshstart)).and.(statuse/=-1)) THEN
      !in case of option -n withouth the l%freshstart option -s, it will l%resume or do a fresh start
      !depending on wether the resuming files are present or not
      IF (l%resume_fromold) THEN
         INQUIRE (file=trim(adjustl(l%nresumeable2))//'.old', EXIST=l%resume)
      ELSE
         INQUIRE (file=trim(adjustl(l%nresumeable2)), EXIST=l%resume)
      END IF
      IF (l%resume) THEN
         IF ((l%layoutnumber==0).or.((l%layoutnumber == l%size/2).and.l%stochastic)) THEN
            !the temporary
            CLOSE (11)
            l%file11isopen=.false.
            OPEN (11, file=trim(adjustl(l%nEntradaRoot))//'_Report.txt', FORM='formatted', POSITION='append')
            l%file11isopen=.true.
 !!!           if (l%layoutnumber==0) call insertalogtmp !ojo lo quito aqui porque borra el _log con la info de credits
            IF (l%resume_fromold) THEN
               CALL print11 (l%layoutnumber, 'Resuming from .fields.old files')
            ELSE
               CALL print11 (l%layoutnumber, 'Resuming from .fields files')
            END IF
         END IF
      ELSE
         !!!IF ((l%layoutnumber==0).or.((l%layoutnumber == l%size/2).and.l%stochastic)) THEN
         !!!   !the temporary
         !!!   CLOSE (11)
         !!!   l%file11isopen=.false.
         !!!   OPEN (11, file=trim(adjustl(l%nEntradaRoot))//'_Report.txt', FORM='formatted')
         !!!   l%file11isopen=.true.
         !!!   if (l%layoutnumber==0) call insertalogtmp
         !!!   CALL print11 (l%layoutnumber, 'Doing a new simulation from n=1')
         !!!END IF
         l%freshstart = .TRUE.
         l%resume_fromold = .FALSE.
      END IF
      IF ( ((l%layoutnumber==0).or.((l%layoutnumber == l%size/2).and.l%stochastic)).and.l%file11isopen) close (11)
      l%file11isopen=.false.
   END IF
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if ((l%run)) then 
#ifdef keeppause
          !!!solo para el cluster
                  INQUIRE (file='running', EXIST=hayinput)
#ifdef CompileWithMPI
                      CALL MPI_Barrier (SUBCOMM_MPI, l%ierr)
#endif
                  if (hayinput) then
                      OPEN (9, file='running', FORM='formatted',action='read')
                      READ (9, '(a)') chain4
                      chain4=trim(adjustl(chain4))
                      CLOSE (9)
                  !!!!!!!!
#ifdef CompileWithMPI
                      CALL MPI_Barrier (SUBCOMM_MPI, l%ierr)
#endif
                      call removeintraspaces(l%opcionespararesumeo)
                      call removeintraspaces(chain4)
                      if (trim(adjustl(l%opcionespararesumeo))==trim(adjustl(chain4))) then
                          existiarunningigual=.true.
                      endif
                  endif
#endif
#ifdef CompileWithMPI
                      CALL MPI_Barrier (SUBCOMM_MPI, l%ierr)
#endif
      if (l%layoutnumber==0) then
         OPEN (38, file='running')
         WRITE (38, '(a)') trim(adjustl(l%opcionespararesumeo))
         CLOSE (38)
      endif
   endif
#ifdef CompileWithMPI
   CALL MPI_Barrier (SUBCOMM_MPI, l%ierr)
#endif
   !open ReportingFiles
   !only the master

   IF (((l%layoutnumber==0).or.((l%layoutnumber == l%size/2).and.l%stochastic)).and.(statuse/=-1)) THEN
      print *,'Opening _Report.txt file'
      IF (l%resume) THEN
         !the temporary
         CLOSE (11)
         l%file11isopen=.false.
         OPEN (11, file=trim(adjustl(l%nEntradaRoot))//'_Report.txt', FORM='formatted')
         l%file11isopen=.true.
         donde = 0
         DO while (donde ==  0)
            !the first one is a dummy read
            READ (11, '(a)') l%chdummy
            donde = index (l%chdummy, 'mpirun -n')
         END DO
         CLOSE (11)
         l%file11isopen=.false.
         l%opcionesoriginales = l%chdummy
!
         call removeintraspaces(l%opcionespararesumeo)
         call removeintraspaces(l%opcionesoriginales)
         IF (trim(adjustl(l%opcionesoriginales)) /= trim(adjustl(l%opcionespararesumeo))) THEN
            CALL stoponerror (l%layoutnumber, l%size, 'Different resumed/original switches: '//trim(adjustl(l%opcionespararesumeo))//' <> '//&
            & trim(adjustl(l%opcionesoriginales)),.true.); statuse=-1; !goto 668
         END IF
         !
         !!!!!!!!!        CLOSE (11, status='delete')
         !!!!!!!!!        l%file11isopen=.false.
         OPEN (11, file=trim(adjustl(l%nEntradaRoot))//'_Report.txt', FORM='formatted')
         l%file11isopen=.true.
         donde = 0
         DO while (donde ==  0)
            !the first one is a dummy read
            READ (11, '(a)') l%chdummy
            donde = index (l%chdummy, '!SLICES')
         END DO
         l%slicesoriginales = trim (adjustl(l%chdummy))
         CLOSE (11)
         l%file11isopen=.false.
         OPEN (11, file=trim(adjustl(l%nEntradaRoot))//'_Report.txt', FORM='formatted', POSITION='append')
         l%file11isopen=.true.
         if (l%layoutnumber==0) call insertalogtmp(l)
      ELSE
         CLOSE (11)
         l%file11isopen=.false.
         OPEN (11, file=trim(adjustl(l%nEntradaRoot))//'_Report.txt', FORM='formatted')
         l%file11isopen=.true.
         if (l%layoutnumber==0) call insertalogtmp(l)
      END IF
      !
      CALL get_secnds (l%time_out2)
      call print_credits(l)
#ifdef CompileWithReal8
      WRITE (dubuf,*) 'Compiled with Double precision (real*8)'
      CALL print11 (l%layoutnumber, dubuf)
#endif      
#ifdef CompileWithReal4
      WRITE (dubuf,*) 'Compiled with Single precision (real*4)'
      CALL print11 (l%layoutnumber, dubuf)
#endif      
#ifdef CompileWithReal16
      WRITE (dubuf,*) 'Compiled with Quadruple precision (real*16)'
      CALL print11 (l%layoutnumber, dubuf)
#endif      
      WRITE (dubuf,*) SEPARADOR // SEPARADOR // SEPARADOR
      !!!CALL print11 (l%layoutnumber, dubuf,.true.)
      write (11,'(a)') trim(adjustl(dubuf)) !a capon para que el l%stochastic pueda resumear
      WRITE (dubuf,*) 'Launched on              ', l%time_out2%fecha(7:8), '/', l%time_out2%fecha(5:6), '/', &
      &                l%time_out2%fecha(1:4), ' ', l%time_out2%hora(1:2), ':', l%time_out2%hora(3:4)
      !!!CALL print11 (l%layoutnumber, dubuf,.true.)
      write (11,'(a)') trim(adjustl(dubuf)) !a capon para que el l%stochastic pueda resumear
      WRITE (dubuf,*) SEPARADOR // SEPARADOR // SEPARADOR
      !!!CALL print11 (l%layoutnumber, dubuf,.true.)
      write (11,'(a)') trim(adjustl(dubuf)) !a capon para que el l%stochastic pueda resumear
      WRITE (dubuf,'(a)') 'Launched with total options '
      !!!CALL print11 (l%layoutnumber, dubuf,.true.)
      write (11,'(a)') trim(adjustl(dubuf)) !a capon para que el l%stochastic pueda resumear
      WRITE (dubuf,*) trim (adjustl(l%opcionestotales))
      !!!CALL print11 (l%layoutnumber, dubuf,.true.)
      write (11,'(a)') trim(adjustl(dubuf)) !a capon para que el l%stochastic pueda resumear
      WRITE (dubuf,'(a)') 'If later resuming use compulsory options '
      !!!CALL print11 (l%layoutnumber, dubuf,.true.)
      write (11,'(a)') trim(adjustl(dubuf)) !a capon para que el l%stochastic pueda resumear
      WRITE (dubuf,*) trim (adjustl(l%opcionespararesumeo))
      !!!CALL print11 (l%layoutnumber, dubuf,.true.)
      write (11,'(a)') trim(adjustl(dubuf)) !a capon para que el l%stochastic pueda resumear
      WRITE (dubuf,*) SEPARADOR // SEPARADOR // SEPARADOR
      CALL print11 (l%layoutnumber, dubuf)
   END IF
   !
   !
   !in seconds
   l%flushsecondsFields = l%flushminutesFields * 60
   !in seconds
   l%flushsecondsData = l%flushminutesData * 60
   
   IF (( .NOT. l%existeNFDE) .AND. ( .NOT. l%existeh5)) THEN
      CALL stoponerror (l%layoutnumber, l%size, 'Some input file missing .h5/.nfde/.conf',.true.); statuse=-1; !goto 668
   END IF
   !
   !
#ifdef CompileWithMPI
   CALL MPI_Barrier (SUBCOMM_MPI, l%ierr)
#endif
   !
#ifdef CompileWithMPI
   CALL MPI_Barrier (SUBCOMM_MPI, l%ierr)
#endif
   !
    if (existiarunningigual) then !lo pongo aqui pq si no no se escribe en el report
        CALL stoponerror (l%layoutnumber, l%size, 'Running flag file with same options than requested exist. ',.true.); statuse=-1;
    endif
              
       
668 continue
                 
    input_conformal_flag=l%input_conformal_flag  !es un flag global!!!!ojooo 051223 !devolverlo correctamente
    return !el unico return que he dejado !240817
   
   end subroutine interpreta
   
   subroutine insertalogtmp(l) !para 100920              
           type (entrada_t), intent(INOUT) :: l  
           CHARACTER (LEN=BUFSIZE) ::  dubuf
           integer (kind=4) :: MYUNIT11
           CALL OffPrint !no reimprimas, esto ya estaba por pantalla
            OPEN (newunit=myunit11, file='SEMBA_FDTD_temp.log')
            do
                read(myunit11,'(1024a)',end=7211) dubuf
                dubuf='&'//dubuf !para respetar los espacios
                CALL print11 (l%layoutnumber, dubuf)
            end do
7211        CLOSE (myunit11, status='delete')   
            CALL OnPrint
            return
   end subroutine insertalogtmp
   
   subroutine print_basic_help(l)                 
      type (entrada_t), intent(INOUT) :: l
      call print_credits(l)
      CALL print11 (l%layoutnumber, '___________________________________________________________________________')
      CALL print11 (l%layoutnumber, 'Basic usage: ')
      CALL print11 (l%layoutnumber, '&   For help use          -h ')
      CALL print11 (l%layoutnumber, '&   For launching use                     ')
      CALL print11 (l%layoutnumber, '&                         -i inputfile (native)')
      CALL print11 (l%layoutnumber, '___________________________________________________________________________')
      return
   end subroutine print_basic_help
 
   subroutine print_credits(l)       
      type (entrada_t), intent(INOUT) :: l
      CHARACTER (LEN=BUFSIZE) ::  dubuf
      
      if (l%creditosyaprinteados) return
      l%creditosyaprinteados=.true.
      CALL print11 (l%layoutnumber, '=========================')
      CALL print11 (l%layoutnumber, 'SEMBA-FDTD SOLVER')
      CALL print11 (l%layoutnumber, '=========================')
      
      WRITE (dubuf,*) SEPARADOR // SEPARADOR // SEPARADOR
      CALL print11 (l%layoutnumber, dubuf)
      WRITE (dubuf,*) trim (adjustl(dataversion))
      CALL print11 (l%layoutnumber, dubuf)
      WRITE (dubuf,*) SEPARADOR // SEPARADOR // SEPARADOR
      CALL print11 (l%layoutnumber, dubuf)
      CALL print11 (l%layoutnumber, 'All rights reserved by the University of Granada (Spain)')
      CALL print11 (l%layoutnumber, '       Contact person: Salvador G. Garcia <salva@ugr.es>')
      CALL print11 (l%layoutnumber, ' ')
      !*******************************************************************************


      WRITE (dubuf,*) SEPARADOR // SEPARADOR // SEPARADOR
      CALL print11 (l%layoutnumber, dubuf)
#ifdef CompileWithMPI
      CALL print11 (l%layoutnumber, 'Compiled WITH MPI support')
#endif
#ifdef CompileWithHDF
      CALL print11 (l%layoutnumber, 'Compiled WITH .h5 HDF support')
#endif
#ifdef CompileWithConformal
      CALL print11 (l%layoutnumber, 'Compiled WITH Conformal support')
#endif
#ifdef CompileWithMTLN
      CALL print11 (l%layoutnumber, 'Compiled WITH MTLN support')
#endif
#ifdef CompileWithSMBJSON
      CALL print11 (l%layoutnumber, 'Compiled WITH SMBJSON support')
#endif
      WRITE (dubuf,*) SEPARADOR // SEPARADOR // SEPARADOR
      CALL print11 (l%layoutnumber, dubuf)
      CALL get_secnds (l%time_out2)
      WRITE (dubuf,*) 'Launched on              ', l%time_out2%fecha(7:8), '/', l%time_out2%fecha(5:6), '/', &
      &                l%time_out2%fecha(1:4), ' ', l%time_out2%hora(1:2), ':', l%time_out2%hora(3:4)
      CALL print11 (l%layoutnumber, dubuf)
      if (l%layoutnumber==0) print *, 'Highest integer ',huge(1_4)
      return
   end subroutine print_credits

   subroutine print_help(l)     
   type (entrada_t), intent(INOUT) :: l
      CHARACTER (LEN=BUFSIZE) :: buff
      CALL print11 (l%layoutnumber, '___________________________________________________________________________')
      CALL print11 (l%layoutnumber, 'Command line arguments: ')
      CALL print11 (l%layoutnumber, '___________________________________________________________________________')
      CALL print11 (l%layoutnumber, '-i geometryfile        : Simulates the Native format input file            ')
      CALL print11 (l%layoutnumber, '-r                     : Restarts a previous execution until a given step. ')
      CALL print11 (l%layoutnumber, '&                        Needs -n                                          ')
      CALL print11 (l%layoutnumber, '-run                   : Uses a semaphore running file and automatically   ')
      CALL print11 (l%layoutnumber, '&                        relaunches simulation if ended or aborted (cluter)')
#ifdef CompileWithOldSaving
      CALL print11 (l%layoutnumber, '-old                   : Jointly with -r restarts from .fields.old files   ')
      CALL print11 (l%layoutnumber, '&                        instead (for safety .fields.old fields are saved  ')
      CALL print11 (l%layoutnumber, '&                        too if -flush is issued)                          ')
#endif
      CALL print11 (l%layoutnumber, '-cfl number            : Courant number (suggested<=0.8)  overriding input ')
      CALL print11 (l%layoutnumber, '-n numberoftimesteps   : Run the simulation until a specified step         ')
      CALL print11 (l%layoutnumber, '&                        either restarting if the necessary files are      ')
      CALL print11 (l%layoutnumber, '&                        present, or starting a fresh new one otherwise    ')
      CALL print11 (l%layoutnumber, '&                        Special cases: n=-1 -> Run only .h5/.nfde preproc.')
      CALL print11 (l%layoutnumber, '&                        Special cases: n=-2 -> Run only .h5 preprocessing ')
      CALL print11 (l%layoutnumber, '-s                     : Forces a fresh new simulation, erasing the        ')
      CALL print11 (l%layoutnumber, '&                        restarting files if they are present              ')
      CALL print11 (l%layoutnumber, '&                        Jointly with -n, it enforces a fresh restart      ')
      CALL print11 (l%layoutnumber, '&                        (erases .fields files from previous simulations)  ')
      CALL print11 (l%layoutnumber, '___________________________________________________________________________')
      CALL print11 (l%layoutnumber, '-pause seconds         : Wait seconds to start simulation                  ')
      CALL print11 (l%layoutnumber, '-prefix string         : Adds a string to the output filenames             ')
      CALL print11 (l%layoutnumber, '-saveall               : Saves all the observation time steps              ')
      CALL print11 (l%layoutnumber, '&                        (default saves only the specified windows of time)')
      CALL print11 (l%layoutnumber, '-singlefile            : Compacts E, H, J probes in single files to        ')
      CALL print11 (l%layoutnumber, '&                        overcome a large number of file openings          ')
      !!#ifdef CompileWithMPI
      !!          CALL print11 (l%layoutnumber, '-maxmessages number    : Buffer of messages for MPI Warnings file. Just    ')
      !!          CALL print11 (l%layoutnumber, '&                        increase if requested at runtime                  ')
      !!#endif

      !*********************************************************************************************************************
      !***[conformal] **************************************************************************************
      !*********************************************************************************************************************
      !conformal -help printf line   ref:  ##Confhelp##
#ifdef CompileWithConformal
      CALL print11 (l%layoutnumber, '-conf                    : Adds the conformal file to the simulation.'        )
#endif
      !*********************************************************************************************************************
#ifdef CompileWithNIBC
      CALL print11 (l%layoutnumber, '-skindepthpre          : Pre-processor for sgbc metals including skin depth.')
      CALL print11 (l%layoutnumber, '-mibc                  : Uses pure l%mibc to deal with composites.  ')
      CALL print11 (l%layoutnumber, '-ade                   : Uses l%ade-l%mibc to deal with composites. ')
      CALL print11 (l%layoutnumber, '&                        Alternative to -l%mibc.'      )
      CALL print11 (l%layoutnumber, '-conformalskin         : Uses a conformal l%mibc to deal with skin-depth')
      CALL print11 (l%layoutnumber, '&                        Do not use this switch if the problem also involves ')
      CALL print11 (l%layoutnumber, '&                        traditional composites, since these do not hold the right ')
      CALL print11 (l%layoutnumber, '&                        thickness parameter. Only use it if the problem only ')
      CALL print11 (l%layoutnumber, '&                        contains metals for which both the conductivity and ')
      CALL print11 (l%layoutnumber, '&                        thickness are CORRECTLY specified in the .nfde file. ')
      CALL print11 (l%layoutnumber, '-NOcompomur            : Uses OLD (possibly unstable) upwinding scheme to deal with composites, ')
      CALL print11 (l%layoutnumber, '&                        instead of the NEW default, which uses a causal time-domain extrapolation ')
      CALL print11 (l%layoutnumber, '&                        of magnetic fields at the surface, by using the one-way ')
      CALL print11 (l%layoutnumber, '&                        advection equation (similar to 1D Mur ABCs) for its ')
      CALL print11 (l%layoutnumber, '&                        superior stability of the default new Mur formulation')
      CALL print11 (l%layoutnumber, '-attc   dissipation    : Positive factor (under 1) for stable composites,   ')
      CALL print11 (l%layoutnumber, '&                        permits to solve some instabilities in the simulation of l%mibc materials.')
      CALL print11 (l%layoutnumber, '&                        It just adds a 1 cell lossy magnetic coating to the l%mibc composite.')
      CALL print11 (l%layoutnumber, '&                        The dissipation factor is used to find the magnetic conductivity ')
      CALL print11 (l%layoutnumber, '&                        from the coefficient updating the current magnetic ')
      CALL print11 (l%layoutnumber, '&                        field from the previous one.  ')       
      write(buff,'(a,e10.2e3)')     '&                        Default= ',l%attfactorc
      CALL print11 (l%layoutnumber, buff)
#endif
      CALL print11 (l%layoutnumber, '-prioritizeCOMPOoverPEC: Uses Composites instead of PEC in conflicts.       ')
      CALL print11 (l%layoutnumber, '-prioritizeISOTROPICBODYoverall: Uses ISOTROPIC BODY FOR conflicts (JUST FOR SIVA).       ')
      CALL print11 (l%layoutnumber, '-sgbc               : Enables the defaults sgbc model for composites. Default sgbc:')
      CALL print11 (l%layoutnumber, '-nosgbc             : Disables the defaults sgbc model for composites. Default sgbc:')
      CALL print11 (l%layoutnumber, '&                        -sgbfreq 3e9 -sgbresol 1 -sgbcrank      ')
      CALL print11 (l%layoutnumber, '-sgbcfreq           : Maximum frequency to consider the skin-depth       ')
      CALL print11 (l%layoutnumber, '-sgbcresol          : Number of cells per skin-depth a the Maximum frequency')
      CALL print11 (l%layoutnumber, '-sgbcyee            : Uses pure Yee ETD sgbc instead of Crank-Nicolson')
      CALL print11 (l%layoutnumber, '-sgbccrank          : Uses sgbc Crank-Nicolson (default)        ')
      CALL print11 (l%layoutnumber, '-sgbcdepth number   : Overrides automatic calculation of number of cells ')
      CALL print11 (l%layoutnumber, '&                        within sgbc                              ')

      CALL print11 (l%layoutnumber, '-pmlalpha factor order : CPML Alpha factor (>=0, <1 sug.) & polyn. grading.')
      CALL print11 (l%layoutnumber, '&                        alpha=factor * maximum_PML_sigma , order=polynom. ')
      write(buff,'(a,2e10.2e3)')    '&                        Default= ',l%alphamaxpar,l%alphaOrden
      CALL print11 (l%layoutnumber, buff)
      write(buff,'(a,e10.2e3)')   '-pmlkappa number       : CPML Kappa (>=1). Default= ',l%kappamaxpar
      CALL print11 (l%layoutnumber, buff)
      CALL print11 (l%layoutnumber, '-pmlcorr factor depth  : Factor for CPML enhanced stability (default none).')
      CALL print11 (l%layoutnumber, '&                        sigma=factor * maximum_PML_sigma, depth= # layers ')
      CALL print11 (l%layoutnumber, '-mur1                  : Supplement PMLs with 1st order Mur ABCs           ')
      CALL print11 (l%layoutnumber, '-mur2                  : Supplement PMLs with 2nd order Mur ABCs           ')
      CALL print11 (l%layoutnumber, '-wiresflavor {holland.or.old} : model for the wires    ')
#ifdef CompileWithBerengerWires
      CALL print11 (l%layoutnumber, '-wiresflavor {berenger} : model for the wires    ')   
#endif
#ifdef CompileWithSlantedWires
      CALL print11 (l%layoutnumber, '-wiresflavor {new/Slanted.or.experimental.or.slanted/transition/semistructured l%precision} : model for the wires    ')   
#endif
      CALL print11 (l%layoutnumber, '&                        (default '//trim(adjustl(l%wiresflavor))//')   ')
      CALL print11 (l%layoutnumber, '-mtlnwires             : Use mtln solver to advance wires currents ')
      CALL print11 (l%layoutnumber, '-notaparrabos          : Do not remove extra double tails at the end of the wires ')
      CALL print11 (l%layoutnumber, '&                        only available for the native format.             ')
      CALL print11 (l%layoutnumber, '-intrawiresimplify     : Disable strict interpretation of .NFDE topology.  ')
      CALL print11 (l%layoutnumber, '&                        Collapse internal parallel wires and create       ')
      CALL print11 (l%layoutnumber, '&                        intra-wire junctions.                             ')
      CALL print11 (l%layoutnumber, '-nomtlnberenger        : Disables MTLN improvements for Berenger l%wiresflavor')
      CALL print11 (l%layoutnumber, '-stableradholland             : Automatic correction of radii for Holland l%wiresflavor')
      CALL print11 (l%layoutnumber, '&                        Use only in case of instabilities.  (experimental)')
      CALL print11 (l%layoutnumber, '-groundwires           : Ground wires touching/embedded/crossing PEC/Lossy.')
      CALL print11 (l%layoutnumber, '&                        Use with CAUTION. Revise *Warnings.txt file!      ')
      CALL print11 (l%layoutnumber, '-noSlantedcrecepelo : Ground open nodes. Experimental. Do not use.')
      CALL print11 (l%layoutnumber, '-connectendings        : Joins ohmicly endings nodes of adjacent segments  ')
      CALL print11 (l%layoutnumber, '&                        from multiwires (segments do no collapse).        ')
      CALL print11 (l%layoutnumber, '&                        regardless of whether they are actually connected ')
      CALL print11 (l%layoutnumber, '&                        through the LeftEnd/RightEnd numbering ')
      CALL print11 (l%layoutnumber, '&                        Automatic with -a                                 ')
      CALL print11 (l%layoutnumber, '&                        Use with CAUTION. Revise *Warnings.txt file!      ')
      CALL print11 (l%layoutnumber, '-isolategroupgroups    : Detach ohmicly endings nodes of adjacent segments ')
      CALL print11 (l%layoutnumber, '&                        from multiwires if they are in different          ')
      CALL print11 (l%layoutnumber, '-makeholes             : Create a void 2-cell area around wire segments    ')
      CALL print11 (l%layoutnumber, '&                        Use with CAUTION. Revise *Warnings.txt (experim.) ')
      CALL print11 (l%layoutnumber, '-mindistwires dist     : Specify the min distance between wires in a       ')
      CALL print11 (l%layoutnumber, '&                        multiwire in new and experimental wires flavors   ')
      write(buff,'(a,e10.2e3)')   '&                        Default= ',l%mindistwires
      CALL print11 (l%layoutnumber, buff)
      CALL print11 (l%layoutnumber, '-inductance {ledfelt/berenger/boutayeb} : model for the self-inductance    ')
      CALL print11 (l%layoutnumber, '&                        (default '//trim(adjustl(l%inductance_model))//')   ')
      CALL print11 (l%layoutnumber, '-inductanceorder order : order for the self-inductance calculation for     ')
      CALL print11 (l%layoutnumber, '&                        slanted wires in experimental l%wiresflavor         ')
      write(buff,'(a,i8)')   '&                        Default= ',l%inductance_order
      CALL print11 (l%layoutnumber, '-attw   dissipation    : Positive factor (under 1) for stability in wires, ')
      write(buff,'(a,e10.2e3)')   '&                        Default= ',l%attfactorw                              
      CALL print11 (l%layoutnumber, '-maxwireradius number  : Bounds globally the wire radius                   ')
      CALL print11 (l%layoutnumber, '-clip                  : Permits to clip a bigger problem truncating wires.')
      CALL print11 (l%layoutnumber, '-wirecrank             : Uses Crank-Nicolson for wires (development)       ')
      CALL print11 (l%layoutnumber, '-noNF2FF string        : Supress a NF2FF plane for calculation             ')
      CALL print11 (l%layoutnumber, '&                        String can be: up, down, left, right, back , front')
      CALL print11 (l%layoutnumber, '-NF2FFDecim            : Uses decimation in NF2FF calculation (faster).    ')
      CALL print11 (l%layoutnumber, '&                        WARNING: High-freq aliasing may occur             ')
      CALL print11 (l%layoutnumber, '-vtkindex              : Output index instead of real point in 3D slices.  ')
      CALL print11 (l%layoutnumber, '-ignoreerrors          : Run even if errors reported in *Warnings.txt file.')
      CALL print11 (l%layoutnumber, '___________________________________________________________________________')
      CALL print11 (l%layoutnumber, '-cpumax minutes        : CPU runtime (useful for limited CPU queuing       ')
      CALL print11 (l%layoutnumber, '-noshared              : Do not waste time with shared fields              ')
      CALL print11 (l%layoutnumber, '-flush minutes         : Minutes between data flush of restarting fields   ')
      CALL print11 (l%layoutnumber, '&                        (default 0=No flush)                              ')
      CALL print11 (l%layoutnumber, '-flushdata minutes     : Minutes between flushing observation data         ')
      CALL print11 (l%layoutnumber, '&                        (default is every 5 minutes)                      ')
      CALL print11 (l%layoutnumber, '-map                   : Creates map ASCII files of the geometry           ')
      CALL print11 (l%layoutnumber, '&                        with wires and PEC                ')
      CALL print11 (l%layoutnumber, '&                        (in conjunction with -n 0 only creates the maps)  ')
      CALL print11 (l%layoutnumber, '-mapvtk                : Creates .VTK map of the PEC/wires/Surface geometry')
#ifdef CompileWithConformal
      CALL print11 (l%layoutnumber, '-conf file             : conformal file  ')
      CALL print11 (l%layoutnumber, '-abrezanjas            : Thin-gaps treated in conformal manner  ')
#endif
      CALL print11 (l%layoutnumber, '-dmma                  : Thin-gaps treated in DMMA manner  ')
#ifdef CompileWithMPI
      CALL print11 (l%layoutnumber, '-mpidir {x,y,z}        : Rotate model to force MPI along z be the largest  ')
      CALL print11 (l%layoutnumber, '-force    cutplane     : Force a MPI layout to begin at cutplane (debug!)  ')
#endif
      CALL print11 (l%layoutnumber, '___________________________________________________________________________')
      CALL print11 (l%layoutnumber, 'Control through signaling files during the simulation: (after erased)      ')
      CALL print11 (l%layoutnumber, '&  stop         : (void) Forces a graceful end (it Cannot be resumed)      ')
      CALL print11 (l%layoutnumber, '&                 No restarting data is flushed, only observation data     ')
      CALL print11 (l%layoutnumber, '&  stopflushing : (void) Forces a graceful end (it can be resumed)         ')
      CALL print11 (l%layoutnumber, '&  flush        : (void) Forces a flush of resuming fields and observation ')
      CALL print11 (l%layoutnumber, '&                 data in 1 minute time approx.                            ')
      CALL print11 (l%layoutnumber, '&  flushdata    : (void) Forces a flush only of the observation data in    ')
      CALL print11 (l%layoutnumber, '&                 1 minute time approx.                                    ')
      CALL print11 (l%layoutnumber, '&                 Both restarting and observation data are flushed         ')
      CALL print11 (l%layoutnumber, '&  stop_only         : Forces a graceful end (cannot be resumed) only of a ')
      CALL print11 (l%layoutnumber, '&                      given problem name (without the .nfde extension)    ')
      CALL print11 (l%layoutnumber, '&                      No restarting data is flushed, only observation data')
      CALL print11 (l%layoutnumber, '&  stopflushing_only : Forces a graceful end (it can be resumed) only of a ')
      CALL print11 (l%layoutnumber, '&                      give problem name (without the .nfde extension)     ')
      CALL print11 (l%layoutnumber, '&                      Both restarting and observation data is flushed     ')
      CALL print11 (l%layoutnumber, '&  flush_only   : Forces flush of resuming fields and observation data only')
      CALL print11 (l%layoutnumber, '&                 of a given problem name (without the .nfde extension)    ')
      CALL print11 (l%layoutnumber, '&                 in 1 minute time approx.                                 ')
      CALL print11 (l%layoutnumber, '&  flushdata_only : Forces a flush only of the observation data only of a  ')
      CALL print11 (l%layoutnumber, '&                   given problem name (without the .nfde extension)       ')
      CALL print11 (l%layoutnumber, '&                   in 1 minute time approx.                               ')
      CALL print11 (l%layoutnumber, '&                   Both restarting and observation data are flushed       ')
      CALL print11 (l%layoutnumber, '&  pause        : (void) While this field exist no simulation is started   ')
      CALL print11 (l%layoutnumber, '&  unpack       : (void) Unpacks on-the-fly .bin probes files created      ')
      CALL print11 (l%layoutnumber, '&                 with the -singlefile packaging option                    ')
      CALL print11 (l%layoutnumber, '&  postprocess  : (void) Do frequency domain and transfer function         ')
      CALL print11 (l%layoutnumber, '&                 postprocess on-the-fly                                   ')
      CALL print11 (l%layoutnumber, '&  flushxdmf    : (void) Flush .xdmf animation probes on the fly           ')
      CALL print11 (l%layoutnumber, '&  flushvtk     : (void) Flush .vtk  animation probes on the fly           ')
      CALL print11 (l%layoutnumber, '&  snap         : Creates a .h5 and .xdmf snapshot per MPI layout if the   ')
      CALL print11 (l%layoutnumber, '&                 field value is over the first number found in this file  ')
      CALL print11 (l%layoutnumber, '&                 in space steps by the 2nd integer number                 ')
      CALL print11 (l%layoutnumber, '&                 in time steps by the 3rd integer number (1-minute lapse) ')
      CALL print11 (l%layoutnumber, '&  relaunch     : Relaunches the simulation upon termination with the      ')
      CALL print11 (l%layoutnumber, '&                 switches read from this file. Used jointly with a        ')
      CALL print11 (l%layoutnumber, '&                 stop file permits to launch simulations on-demand        ')
      CALL print11 (l%layoutnumber, '___________________________________________________________________________')
      !
      write (buff,'(a,i14,a)') 'Max CPU time is ',topCPUtime,' seconds (can be overriden by -cpumax)'
      CALL print11 (l%layoutnumber, buff)
#ifdef CompileWithOpenMP
      CALL print11 (l%layoutnumber, 'SUPPORTED:   MultiCPU parallel simulation (OpenMP)')
#else
      !CALL print11 (l%layoutnumber, 'UNSUPPORTED: MultiCPU parallel simulation (OpenMP)')
#endif
!
#ifdef CompileWithMPI
      CALL print11 (l%layoutnumber, 'SUPPORTED:   MultiCPU/Multinode parallel simulation (MPI)')
#else
      !CALL print11 (l%layoutnumber, 'UNSUPPORTED: MultiCPU/Multinode parallel simulation (MPI)')
#endif
#ifdef CompileWithConformal
      CALL print11 (l%layoutnumber, 'SUPPORTED:   Conformal algorithm')
#else
      !CALL print11 (l%layoutnumber, 'UNSUPPORTED: Conformal algorithm')
#endif
      CALL print11 (l%layoutnumber, 'SUPPORTED:   Near-to-Far field probes')
      CALL print11 (l%layoutnumber, 'SUPPORTED:   Lossy anistropic materials, both electric and magnetic')
      CALL print11 (l%layoutnumber, 'SUPPORTED:   Thin Slots ')
      CALL print11 (l%layoutnumber, 'SUPPORTED:   Electric and Magnetic Dispersive materials ')
      CALL print11 (l%layoutnumber, 'SUPPORTED:   Isotropic Multilayer Skin-depth Materials (sgbc)')
#ifdef CompileWithNIBC
      CALL print11 (l%layoutnumber, 'SUPPORTED:   Isotropic Multilayer Skin-depth Materials (l%mibc)')
#else
      !CALL print11 (l%layoutnumber, 'UNSUPPORTED: Isotropic Multilayer Skin-depth Materials (l%mibc)')
#endif
      CALL print11 (l%layoutnumber, 'SUPPORTED:   Loaded and grounded thin-wires with juntions')
      CALL print11 (l%layoutnumber, 'SUPPORTED:   Nodal hard/soft electric and magnetic sources')
#ifdef CompileWithHDF
      CALL print11 (l%layoutnumber, 'SUPPORTED:   .xdmf+.h5 probes ')
#else
      !CALL print11 (l%layoutnumber, 'UNSUPPORTED: .xdmf+.h5 probes ')
#endif
#ifdef CompileWithOldSaving
      CALL print11 (l%layoutnumber, 'SUPPORTED:   .fields.old files created (fail-safe)')
#else
      !CALL print11 (l%layoutnumber, 'UNSUPPORTED: .fields.old files created (fail-safe)')
#endif
#ifdef CompileWithStochastic
      CALL print11 (l%layoutnumber, 'SUPPORTED:   l%stochastic analysis')
#else
      !CALL print11 (l%layoutnumber, 'UNSUPPORTED: l%stochastic analysis')
#endif
#ifdef CompileWithPrescale
      CALL print11 (l%layoutnumber, 'SUPPORTED:   Permittivity scaling accelerations')
#else
      !CALL print11 (l%layoutnumber, 'UNSUPPORTED: Permittivity scaling accelerations')
#endif
      CALL print11 (l%layoutnumber, 'SUPPORTED:   Holland Wires')
#ifdef CompileWithBerengerWires
      CALL print11 (l%layoutnumber, 'SUPPORTED:   Multi-Wires')
#else
      !CALL print11 (l%layoutnumber, 'UNSUPPORTED: Multi-Wires')
#endif
#ifdef CompileWithSlantedWires
      CALL print11 (l%layoutnumber, 'SUPPORTED:   Slanted Wires')
#else
      !CALL print11 (l%layoutnumber, 'UNSUPPORTED: Slanted Wires')
#endif
!!!!!!!!!!!!!!!!!
#ifdef CompileWithReal4
      CALL print11 (l%layoutnumber, 'Single precission simulations (reals are 4-byte)')
#endif
#ifdef CompileWithReal8
      CALL print11 (l%layoutnumber, 'Double precission simulations (reals are 8-byte)')
#endif
#ifdef CompileWithInt4
      CALL print11 (l%layoutnumber, 'Media matrices are 4 bytes')
#endif
#ifdef CompileWithInt2
      CALL print11 (l%layoutnumber, 'Media matrices are 2 bytes')
#endif
#ifdef CompileWithInt1
      CALL print11 (l%layoutnumber, 'Media matrices are 1 byte')
#endif
#ifdef CompileWithMPI
      CALL MPI_FINALIZE (l%ierr)
#endif
      return
   end subroutine print_help
 
   subroutine removeintraspaces(a)
        CHARACTER (LEN=*), intent(inout):: a
        integer (Kind=4) :: i,longi
        logical correc
        correc=.true.
        do while(correc)
            correc=.false.
            a=trim(adjustl(a))
            longi=len(trim(adjustl(a)))
            buscae: do i=1,longi-1
              if ((a(i:i)==' ').and.(a(i+1:i+1)==' ')) then
                a=a(1:i-1)//a(i+1:longi) 
                correc=.true.
                exit buscae
              endif
            end do buscae
        end do
   return
   end subroutine removeintraspaces


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine buscaswitchficheroinput(l)     
 
!!!!!!!!!!!!!        
   type (entrada_t), intent(INOUT) :: l
!!!!!!!!!      
   
   CHARACTER (LEN=BUFSIZE) :: dato,buff,f
   integer (kind=4) :: i,n,statuse, NUM_NFDES,TEMP_NUMNFDES,p
   CHARACTER (LEN=5) :: NFDEEXTENSION, CONFEXTENSION, CMSHEXTENSION
    

!!!
   
   NFDEEXTENSION='.nfde'; CONFEXTENSION='.conf'; CMSHEXTENSION='.cmsh'
   statuse=0
   !!!!!!!!!!!!!!!
   n = commandargumentcount (l%chain2)
   IF (n == 0) THEN
      call print_basic_help(l)  
      call stoponerror(l%layoutnumber,l%size,'Error: NO arguments neither command line nor in launch file. Correct and remove pause...',.true.)
      statuse=-1
      goto 667
   END IF

   IF (n > 0) THEN
      num_nfdes=0
      i = 2
      DO while (i <= n)
         CALL getcommandargument (l%chain2, i, l%chain, l%length, statuse)
         IF (statuse /= 0) THEN
            CALL stoponerror (l%layoutnumber, l%size, 'Reading input',.true.)
            goto 667
         END IF
         !
         SELECT CASE (trim(adjustl(l%chain)))
          CASE ('-mpidir')
            i = i + 1
            CALL getcommandargument (l%chain2, i, f, l%length,  statuse)
            select case (trim (adjustl(f)))
             case ('x','X')
               l%mpidir=1  !!!lo cambie por error !161018
             case ('y','Y')
               l%mpidir=2   !!!lo cambie por error !161018
             case ('z','Z')
               l%mpidir=3
             CASE DEFAULT
               GOTO 1762
            END SELECT
            GO TO 2762
1762        CALL stoponerror (l%layoutnumber, l%size, 'Invalid -l%mpidir option',.true.)
            statuse=-1
            goto 667
2762        CONTINUE
          CASE ('-h')
            call print_credits(l)
            call print_help(l)
            call print_credits(l)
            STOP
          CASE ('-hh')
            call print_credits(l)
            call print_help(l)
            call print_credits(l)
            STOP
          CASE ('-i')
            num_nfdes=num_nfdes + 1
         end select
         i = i + 1
      END DO
      if (num_nfdes > 1) then
         temp_numnfdes=0
         i = 2 ! se empieza en 2 porque el primer argumento es siempre el nombre del ejecutable
         DO while (i <= n)
            CALL getcommandargument (l%chain2, i, l%chain, l%length, statuse)
             IF (statuse /= 0) THEN
                CALL stoponerror (l%layoutnumber, l%size, 'Reading input',.true.)
                goto 667
             END IF
            !
            SELECT CASE (trim(adjustl(l%chain)))
             CASE ('-i')
               temp_numnfdes=temp_numnfdes + 1
               i = i + 1
               CALL getcommandargument (l%chain2, i, f, l%length,  statuse)
               p = LEN_trim (adjustl(f))
               IF ((p-4) >= 1) THEN
                  IF (f((p-4) :(p-4)) == NFDEEXTENSION(1:1)) THEN
                     NFDEEXTENSION= f((p-4) :p)
                     l%extension=NFDEEXTENSION 
                     l%fichin = f (1:p-5)
                  ELSE
                     l%fichin = f (1:p)
                  END IF
               ELSE IF (p >= 1) THEN
                  l%fichin = f (1:p)
               ELSE
                  CALL stoponerror (l%layoutnumber, l%size, 'There is not a .nfde file for input',.true.)
                  statuse=-1
                  goto 667
               END IF
               INQUIRE (file=trim(adjustl(l%fichin))//NFDEEXTENSION, EXIST=l%existeNFDE)
               IF ( .NOT. l%existeNFDE) THEN
                  buff='The input file was not found '//trim(adjustl(l%fichin))//NFDEEXTENSION
                  CALL stoponerror (l%layoutnumber, l%size, buff,.true.)
                  statuse=-1
                  goto 667  
               END IF
!aniadido para chequear que no haya .conf sin haber invocado el -conf 15/12/16 sgg
               INQUIRE (file=trim(adjustl(l%fichin))//CONFEXTENSION, EXIST=l%existeconf)
               IF ((l%existeconf).AND.(.not.(l%input_conformal_flag))) THEN
                  buff='No -conf issued but existing file '//trim(adjustl(l%fichin))//confEXTENSION//' . Either remove file or relaunch with -conf'
                  CALL stoponerror (l%layoutnumber, l%size, buff,.true.)
                  statuse=-1
                  goto 667  
               END IF
               INQUIRE (file=trim(adjustl(l%fichin))//CMSHEXTENSION, EXIST=l%existecmsh)
               IF ((l%existecmsh).AND.(.not.(l%input_conformal_flag))) THEN
                  buff='No -conf issued but existing file '//trim(adjustl(l%fichin))//CMSHEXTENSION//' . Either remove file or relaunch with -conf'
                  CALL stoponerror (l%layoutnumber, l%size, buff,.true.)
                  statuse=-1
                  goto 667  
               END IF
!
               if (temp_numnfdes ==1) then !solo el primero
                  IF (l%layoutnumber == 0) open (194,file='multi_'//trim(adjustl(l%fichin))//NFDEEXTENSION,form='formatted')
               endif
               IF (l%layoutnumber == 0) then
                  open (196,file=trim(adjustl(l%fichin))//NFDEEXTENSION,form='formatted')
                  do
                     read (196,'(a)',end=197) dato
                     if (trim(adjustl(dato)) /= '!END') then
                        write(194,'(a)') trim(adjustl(dato))
                     else
                        dato='***** End merging file: '//trim(adjustl(l%fichin))//NFDEEXTENSION//' ********'
                        write(194,'(a)') trim(adjustl(dato))
                     endif
                  end do
197               close(196)
               endif
               if (temp_numnfdes == num_nfdes) then !solo el primero
                  IF (l%layoutnumber == 0) then
                     write(194,'(a)') '!END'
                     close (194)
                  endif
               endif
            end select
            i = i + 1
         END DO
      endif
   endif
#ifdef CompileWithMPI
   CALL MPI_Barrier (SUBCOMM_MPI, l%ierr)
#endif
   !
   !concatenado multiples ORIGINAL 26/06/14
   
   !fin concatenado

   temp_numnfdes=0
   IF (n > 0) THEN
   i = 2  ! se empieza en 2 porque el primer argumento es siempre el nombre del ejecutable
      DO while (i <= n)
         CALL getcommandargument (l%chain2, i, l%chain, l%length, statuse)
         IF (statuse /= 0) THEN
            CALL stoponerror (l%layoutnumber, l%size, 'Reading input',.true.)
            goto 667
         END IF
         !
         SELECT CASE (trim(adjustl(l%chain)))
            !
          CASE ('-i')
            temp_numnfdes=temp_numnfdes + 1
            i = i + 1
            if (temp_numnfdes == 1) then
               !
               CALL getcommandargument (l%chain2, i, f, l%length,  statuse)
               p = LEN_trim (adjustl(f))
               IF ((p-4) >= 1) THEN
                  IF (f((p-4) :(p-4)) == NFDEEXTENSION(1:1)) THEN
                     NFDEEXTENSION= f((p-4) :p)
                     l%extension=NFDEEXTENSION 
                     l%fichin = f (1:p-5)
                  ELSE
                     l%fichin = f (1:p)
                  END IF
               ELSE IF (p >= 1) THEN
                  l%fichin = f (1:p)
               ELSE
                  CALL stoponerror (l%layoutnumber, l%size, 'There is not a .nfde file for input',.true.)
                  statuse=-1
                  goto 667
               END IF
               INQUIRE (file=trim(adjustl(l%fichin))//NFDEEXTENSION, EXIST=l%existeNFDE)
               IF ( .NOT. l%existeNFDE) THEN
                  buff='The input file was not found '//trim(adjustl(l%fichin))//NFDEEXTENSION
                  CALL stoponerror (l%layoutnumber, l%size, buff,.true.)
                  statuse=-1
                  goto 667
               END IF
            elseif (temp_numnfdes == 2) then
               l%fichin='multi_'//trim(adjustl(l%fichin))
            else
               l%fichin=l%fichin
            endif
            !!!          l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain)) // ' ' // trim (adjustl(f))
         END SELECT
         i = i + 1
      END DO
   endif
      !
    ! If no input is present we stop
    IF (len(trim(adjustl(l%fichin))) <= 0) THEN
        CALL stoponerror (l%layoutnumber, l%size, 'ERROR! -> No input file was specified. Use -i ****.nfde',.true.); statuse=-1; goto 667
    END IF

   l%fileFDE = trim (adjustl(l%fichin)) // NFDEEXTENSION
   l%fileH5 = trim (adjustl(l%fichin)) // '.h5'
   CALL INITWARNINGFILE (l%layoutnumber, l%size, trim (adjustl(l%fichin))//'_tmpWarnings.txt',l%verbose,l%ignoreerrors)
   

!!!
667   return
  end subroutine buscaswitchficheroinput

!!!!!!!!!!!!!!!!!!   
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine default_flags(l)   
!!!!!!!!!!!!!        
   type (entrada_t), intent(INOUT) :: l    
      l%noconformalmapvtk=.false.
      l%forced=-1
      l%sgbcdepth=-1
      l%statuse=0
      l%time_begin=0
      l%ficherohopf=''
!
      l%hopf=.false.
      l%precision=0 !redondeo del semiestructurado
      l%stochastic=.false.
      l%chosenyesornostochastic=.false. !es un flag informativo que debe inicializarse a .false. a pesar de qu el sentido comun diga lo contrario
      l%simu_devia=.false. !puto bug semana santa 2019 cazado a 210419
   
#ifdef CompileWithHDF
      l%createh5bin=.false.
#else
      l%createh5bin=.true.
#endif     
      l%createh5filefromsinglebin=.false.
      l%permitscaling=.false.
      l%niapapostprocess=.false.
      l%prioritizeCOMPOoverPEC=.false.  !pec has default more priority than compo (para siva hay que cambiarlo)
      l%prioritizeTHINWIRE=.false. !solo para visualizacion y experimentacion 231024
      l%prioritizeISOTROPICBODYoverall=.FALSE. !PARA EL SIVA SE CAMBIA POR LINEA DE COMANDO
      l%mpidir=3 !DEFAULT DO NOT ROTATE GEOMETRY !JUST TO TAKE PROFIT OF MPI
      l%maxwireradius=-1.0_RKIND
      l%boundwireradius=.false.
      l%wirecrank=.FALSE.
      l%ignoreerrors=.false.
      l%ignoresamplingerrors=.false.
      l%vtkindex=.FALSE. !SOLO AFECTA A LOS VTK (SACA INDICES EN VEZ DE POSICION FISICA)
      l%CLIPREGION=.false.
      l%NF2FFDecim=.FALSE.
      l%facesNF2FF%tr=.true.
      l%facesNF2FF%fr=.true.
      l%facesNF2FF%iz=.true.
      l%facesNF2FF%de=.true.
      l%facesNF2FF%ab=.true.
      l%facesNF2FF%ar=.true.
      !defaults
      l%use_mtln_wires = .false.
      l%hay_slanted_wires=.false.
      l%forcing = .FALSE.
      l%resume_fromold = .FALSE.
      l%singlefilewrite = .FALSE.
      l%updateshared=.true. !040717 para que no pierda tiempo en el update shared de preprocess se crea flag -noshared que pone esta variable a false. 
                          !Esa info solo la usan anisotropic y composites. podria hacerse algo mas automatico. de momento manual.

      l%finaltimestep=0
      l%cfltemp=1.0 !dummy
      l%cfl=1.0 !default courant number !no tocarlo 310715 solo afecta si se usa -l%cfl
      l%forcecfl=.false.
      !PML default
      !cpml stretching maximum parameters !!l%alphamaxpar=StaticFrequency*2*pi*Eps0
      l%alphamaxpar=0.0_RKIND !0.24  !expresion 7.78 taflove 3 edic)
      l%alphaOrden=1.0_RKIND
      l%kappamaxpar=1.0_RKIND !15.0_RKIND !061118 mantener a 1 por conflictos cpml and permittivity scaling
      !and final layer electric sigma
      l%MEDIOEXTRA%exists=.false.
      l%MEDIOEXTRA%index=-7 !void
      l%MEDIOEXTRA%size=-1  !void
      l%MEDIOEXTRA%sigma=-1e20 !void
      !
      l%MurAfterPML=.false.
      l%mur_second=.false.
      l%mur_first=.false.
      l%mur_exist=.false.
      !!!!!!!!!!!!
      l%takeintcripte=.false. !a peticion de OLD, redondear los nodos de cripte a la baja
      l%attfactorc=1.0_RKIND !default dissipation factor for composites
      l%attfactorw=1.0_RKIND!default dissipation factor for wires

      l%mibc=.false.
      l%ade=.false. !auxiliary differential equation for composites
      l%conformalskin=.false.
      l%sgbc=.true. !default is false unless required
      l%sgbcDispersive=.false. !default is false unless required    
      l%skindepthpre=.false.
      l%sgbcdepth=-1  ! se calcula automaticamente a menos que se use el switch
      l%sgbcfreq=1e9 !default es cazar el skin depth hasta 1e9
      l%sgbcresol=1.0 !numero de celdas por skin depth (caida a exp(-1))
      l%sgbccrank=.true. !default es l%sgbccrank

      l%fatalerror=.false.
      l%fatalerrornfde2sgg=.false.
      !**************************************************************************************************
      !***[conformal] *******************************************************************
      !**************************************************************************************************
      !conformal existence flags   ref: ##Confflag##
      l%input_conformal_flag = .false.
      !**************************************************************************************************
      !**************************************************************************************************
      !**************************************************************************************************
      !added 2701418 para prevenir conflictos de dobles mallas conformal si se usa -l%run_with_abrezanjas
      l%flag_conf_sgg=.false.
      !

      l%dontwritevtk=.false.

      l%NOcompomur=.false. !DEFAULT mi formulacion

      !default no join the wires which are adjacent (ORIGINAL election)
      !do not connect endings unless specified in ORIGINAL
      l%makeholes = .FALSE.
      l%connectendings=.false.
      l%isolategroupgroups=.false.
      l%strictOLD=.true. !default is strict ORIGINAL overriden manually 
      l%TAPARRABOS=.true. !default since 101116 !cortar los multizigzag rabitos
      l%mtlnberenger=.true. !solo actua si se invoca con l%wiresflavor berenger esto va a ser siempre true a menos que tambien se invoque con -nomtlnberenger (solo para debugeo y que coincida con Holland) 020719
      l%stableradholland=.false. !solo actua si se invoca con l%wiresflavor holland 
      l%fieldtotl=.false.
      l%experimentalVideal=.false.
      l%forceresampled=.false.
      l%factorradius=1.0e+30 !para evitar division por cero 120123
      l%factordelta=1.0e+30 !para evitar division por cero 120123
      !default
      l%groundwires = .false.
      l%noSlantedcrecepelo =.false. !131219 experimental niapa ojoooo
      l%inductance_model = 'boutayeb'
      l%inductance_order = 8
      l%wiresflavor='holland'
      l%wirethickness=1
      l%mindistwires = 0.5_RKIND
      !
      l%MurAfterPML = .false.
      !
      l%createmap = .FALSE.
      l%createmapvtk = .FALSE.
      l%verbose = .FALSE.
      l%saveall = .FALSE.
      l%forcesteps = .FALSE.
      l%resume = .FALSE.
      l%freshstart = .FALSE.
      l%run = .FALSE.  !si hay .fields restartea y si no comienza
      l%deleteintermediates = .FALSE.
      !
      l%existeNFDE = .FALSE.
      l%existeh5  = .FALSE.
      !
      !default is NO flush fields
      l%flushminutesFields = 0
      !default is to flush data when the buffer is filled up
      !si se pone cada tantos minutos y se guardan las sondas en trancos!puede haber errores de redondeo porque el buffer se limpia tras cada flusheo
      l%flushminutesData = topCPUtime
      !
      !maximum runtime
      l%maxCPUtime = topCPUtime
      l%input_conformal_flag = .false.
      l%file11isopen=.false.
      l%relaunching=.false.
      l%forcestop=.false.
      l%input_conformal_flag = .false.
      l%run_with_dmma = .true.
#ifdef CompileWithConformal
      l%run_with_dmma = .false.
! todo esto para el abrezanjas. se precisa tambien el l%input_conformal_flag  
!!!!quitado sgg ojo 290521 esto no se ha arreglado aim... quito el abrezanjas !290521 bug     
  !!!    l%run_with_abrezanjas = .true. !OJO 0323 A VECES DA ERROR. PONER A FALSE SI SUCEDE
      l%run_with_abrezanjas = .false. !OJO 0323 A VECES DA ERROR. PONER A FALSE SI SUCEDE
      !!!!l%run_with_abrezanjas = .false.
#else
      l%run_with_abrezanjas = .false.
#endif

!fin thin gaps

      input_conformal_flag=l%input_conformal_flag    !ojooo 051223 es un flag globaaaallll
      return
   end subroutine default_flags




   
   end module interpreta_switches_m

