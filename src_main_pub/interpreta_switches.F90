module interpreta_switches_m

   use FDETYPES_m
   use Getargs_m
   use EpsMuTimeScale_m
   use Report_m
   use version_m

   implicit none
   private
   !
   type entrada_t

      logical :: &
         forcing, &
         singlefilewrite, &
         ignoresamplingerrors, &
         ignoreerrors, &
         updateshared, &
         prioritizeISOTROPICBODYoverall, &
         wirecrank, &
         CLIPREGION, &
         verbose, &
         resume, &
         forcesteps, &
         resume_fromold, &
         freshstart, &
         run, &
         createmap, &
         dontwritevtk, &
         vtkindex, &
         createmapvtk, &
         run_with_dmma, &
         run_with_abrezanjas, &
         input_conformal_flag, &
         pausar, &
         relaunching, &
         forcestop, &
         l_aux, &
         flag_conf_sgg, &
         takeintcripte, &
         skindepthpre, &
         sgbc, &
         conformalskin, &
         ade, &
         mibc, &
         NOcompomur, &
         MurAfterPML, &
         sgbccrank, &
         sgbcDispersive, &
         saveall, &
         boundwireradius, &
         hay_slanted_wires, &
         makeholes, &
         mur_first, &
         mur_second, &
         mur_exist, &
         connectendings, &
         strictOLD, &
         mtlnberenger, &
         stableradholland, &
         TAPARRABOS, &
         fieldtotl, &
         forceresampled, &
         isolategroupgroups, &
         groundwires, &
         noSlantedcrecepelo, &
         forcecfl, &
         niapapostprocess, &
         permitscaling, &
         stochastic, &
         chosenyesornostochastic, &
         prioritizeCOMPOoverPEC, &
         prioritizeTHINWIRE, &
         createh5bin, &
         deleteintermediates, &
         existeNFDE, &
         file11isopen, &
         NF2FFDecim, &
         existeh5, &
         fatalerror, &
         fatalerrornfde2sgg, &
         existeconf, &
         existecmsh, &
         thereare_stoch, &
         experimentalVideal, &
         simu_devia, &
         noconformalmapvtk, &
         createh5filefromsinglebin, &
         creditosyaprinteados, &
         read_command_line

      integer(kind=4) :: &
         wirethickness, &
         inductance_order, &
         finaltimestep, &
         ierr, &
         layoutnumber, &
         num_procs, &
         length, &
         mpidir, &
         flushminutesFields, &
         flushminutesData, &
         flushsecondsFields, &
         flushsecondsData, &
         forced, &
         maxCPUtime, &
         sgbcdepth, &
         precision, &
         statuse

      real(kind=RKIND) :: &
         maxwireradius, &
         mindistwires, &
         attfactorc, &
         attfactorw, &
         cfltemp, &
         cfl, &
         sgbcfreq, &
         sgbcresol, &
         alphamaxpar, &
         kappamaxpar, &
         alphaOrden

      real(kind=8) :: &
         time_begin, &
         time_end
      real(kind=RKIND_wires) :: &
         factorradius, &
         factordelta

      type(nf2ff_T) :: facesNF2FF
      type(MedioExtra_t) :: MEDIOEXTRA
      type(EpsMuTimeScale_input_parameters_t) :: EpsMuTimeScale_input_parameters
      type(tiempo_t) :: time_out2

!pgi        character(len=BUFSIZE_LONG) :: &
      character(len=BUFSIZE) :: &
         prefix, &
         prefixopci, &
         prefixopci1, &
         opcionespararesumeo, &
         opcionesoriginales, &
         slicesoriginales, &
         chdummy

      character(len=BUFSIZE) :: chaininput, &
                                chain, &
                                chain2, &
                                fichin, &
                                extension, &
                                nresumeable2, &
                                fileFDE, &
                                fileH5, &
                                inductance_model, &
                                wiresflavor, &
                                nEntradaRoot, &
                                opcionestotales, &
                                conformal_file_input_name, &
                                geomfile

   end type entrada_t

   public interpreta, insertalogtmp, print_help, print_basic_help, print_credits, &
      removeintraspaces, buscaswitchficheroinput, default_flags
   public entrada_t
   !
contains

   subroutine interpreta(l, statuse)

!!!!!!!!!!!!!
      type(entrada_t), intent(INOUT) :: l
      integer(kind=4), intent(out) :: statuse
!!!!!!!!!

      character(len=BUFSIZE) :: chari, f, dubuf, buff, binaryPath
      logical :: existiarunningigual, mpidirset, resume3
      integer(kind=4) :: i, j, donde, n, newmpidir
      real(kind=RKIND) :: pausetime
      integer(kind=4) :: iostatus = 0

!!!
      l%input_conformal_flag = input_conformal_flag !ojooo 051223 es un flag global

      mpidirset = .false.
      existiarunningigual = .false.
      statuse = 0
   !!!!!!!!!!!!!!!
      binaryPath = getBinaryPath()
      n = commandargumentcount(l%chaininput, binaryPath)
      if (n == 0) then
         call print_basic_help(l)
      call stoponerror(l%layoutnumber,l%num_procs,'Error: NO arguments neither command line nor in launch file. Correct and remove pause...',.true.)
         statuse = -1
      end if
      l%opcionestotales = ''
      do i = 1, n
         call getcommandargument(l%chaininput, i, l%chain, l%length, statuse, binaryPath)
         if (statuse /= 0) then
            call stoponerror(l%layoutnumber, l%num_procs, 'Reading input', .true.)
            statuse = -1
         end if
         l%opcionestotales = trim(adjustl(l%opcionestotales))//' '//trim(adjustl(l%chain))
      end do
      call print11(l%layoutnumber, 'Switches '//trim(adjustl(l%opcionestotales)))

      if (n > 0) then
         i = 2  ! se empieza en 2 porque el primer argumento es siempre el nombre del ejecutable
         do while (i <= n)
            call getcommandargument(l%chaininput, i, l%chain, l%length, statuse, binaryPath)
            if (statuse /= 0) then
               call stoponerror(l%layoutnumber, l%num_procs, 'Reading input', .true.)
               statuse = -1
            end if
            select case (trim(adjustl(l%chain)))
            case ('-i')
               i = i + 1
               call getcommandargument(l%chaininput, i, f, l%length, statuse, binaryPath)
               continue !ya interpretado
            case ('-a')
               i = i + 1
               call getcommandargument(l%chaininput, i, f, l%length, statuse, binaryPath)
               continue !ya interpretado
            case ('-mpidir')
               i = i + 1
               call getcommandargument(l%chaininput, i, f, l%length, statuse, binaryPath)
               select case (trim(adjustl(f)))
               case ('x', 'X')
                  l%mpidir = 1  
               case ('y', 'Y')
                  l%mpidir = 2   
               case ('z', 'Z')
                  l%mpidir = 3
               case default
                  call stoponerror(l%layoutnumber, l%num_procs, 'Invalid or duplicate incoherent -l%mpidir option', .true.)
                  statuse = -1
                  continue
               end select
               if (.not. mpidirset) then
                  l%opcionespararesumeo = trim(adjustl(l%opcionespararesumeo))//' '//trim(adjustl(l%chain))//' '//trim(adjustl(f))
                  mpidirset = .true.
               end if

            case ('-pause')
               i = i + 1
               call getcommandargument(l%chaininput, i, f, l%length, statuse, binaryPath)
               read (f, *, iostat=iostatus) pausetime
               if (iostatus /= 0) call stoponerror(l%layoutnumber, l%num_procs, 'Invalid pause time', .true.)
               if (pausetime <= 0) then
                  call stoponerror(l%layoutnumber, l%num_procs, 'Invalid pause time: zero or negative value', .true.)
                  statuse = -1
               end if
               !
               l%pausar = .true.
#ifdef CompileWithMPI
               call MPI_Barrier(SUBCOMM_MPI, l%ierr)
#endif
               call get_secnds(l%time_out2)
               l%time_begin = l%time_out2%segundos
               write(dubuf, *) 'Paused for (secs) ', pausetime
               call print11(l%layoutnumber, dubuf)
               do while (l%pausar)
#ifdef CompileWithMPI
                  call MPI_Barrier(SUBCOMM_MPI, l%ierr)
#endif
                  call get_secnds(l%time_out2)
                  l%time_end = l%time_out2%segundos
                  if (l%time_end - l%time_begin > pausetime) then
                     l%pausar = .false.
                  end if
               end do
#ifdef CompileWithMPI
               call MPI_Barrier(SUBCOMM_MPI, l%ierr)
               l%l_aux = l%pausar
               call MPI_AllReduce(l%l_aux, l%pausar, 1_4, MPI_LOGICAL, MPI_LOR, SUBCOMM_MPI, l%ierr)
#endif
            case ('-NF2FFDecim')
               l%NF2FFDecim = .true.
               l%opcionespararesumeo = trim(adjustl(l%opcionespararesumeo))//' '//trim(adjustl(l%chain))//' '//trim(adjustl(f))
            case ('-noNF2FF')
               i = i + 1
               call getcommandargument(l%chaininput, i, f, l%length, statuse, binaryPath)
               select case (trim(adjustl(f)))
               case ('back', 'BACK')
                  l%facesNF2FF%TR = .FALSE.
               case ('front', 'FRONT')
                  l%facesNF2FF%FR = .FALSE.
               case ('left', 'LEFT')
                  l%facesNF2FF%IZ = .FALSE.
               case ('right', 'RIGHT')
                  l%facesNF2FF%DE = .FALSE.
               case ('down', 'DOWN')
                  l%facesNF2FF%AB = .FALSE.
               case ('up', 'UP')
                  l%facesNF2FF%AR = .FALSE.
               case default
                  call stoponerror(l%layoutnumber, l%num_procs, 'Invalid -noNF2FF option', .true.)
                  statuse = -1
               end select
               continue
               !COMO LA RCS SE CALCULA SOLO AL FINAL NO OBLIGO A RESUMEAR CON IGUAL -NONFF2FF PARA PODER CALCULAR CON Y SIN ESTA OPCION resumeando
               !          l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain)) // ' ' // trim (adjustl(f))
            CASE ('-force')
               l%forcing = .TRUE.
               i = i + 1
               call getcommandargument(l%chaininput, i, f, l%length, statuse, binaryPath)
               READ (f, *, ERR=412) l%forced
               GO TO 312
412            call stoponerror(l%layoutnumber, l%num_procs, 'Invalid cut', .true.)
               statuse = -1
312            continue
               l%opcionespararesumeo = trim(adjustl(l%opcionespararesumeo))//' '//trim(adjustl(l%chain))//' '//trim(adjustl(f))
            CASE ('-singlefile')
               l%singlefilewrite = .TRUE.
               l%opcionespararesumeo = trim(adjustl(l%opcionespararesumeo))//' '//trim(adjustl(l%chain))
            CASE ('-ignoresamplingerrors')
               l%ignoresamplingerrors = .TRUE.
            CASE ('-prioritizeTHINWIRE')
               l%prioritizeTHINWIRE = .true.
               l%opcionespararesumeo = trim(adjustl(l%opcionespararesumeo))//' '//trim(adjustl(l%chain))
               l%ignoreerrors = .TRUE.
            CASE ('-prioritizeCOMPOoverPEC')
               l%prioritizeCOMPOoverPEC = .true.
               l%opcionespararesumeo = trim(adjustl(l%opcionespararesumeo))//' '//trim(adjustl(l%chain))
               l%ignoreerrors = .TRUE.
            CASE ('-noshared')
               l%updateshared = .false.
               l%opcionespararesumeo = trim(adjustl(l%opcionespararesumeo))//' '//trim(adjustl(l%chain))
            CASE ('-prioritizeISOTROPICBODYoverall')
               l%prioritizeISOTROPICBODYoverall = .true.
               l%opcionespararesumeo = trim(adjustl(l%opcionespararesumeo))//' '//trim(adjustl(l%chain))
            CASE ('-wirecrank')
               l%wirecrank = .TRUE.
               l%opcionespararesumeo = trim(adjustl(l%opcionespararesumeo))//' '//trim(adjustl(l%chain))
            CASE ('-clip')
               l%CLIPREGION = .TRUE.
               l%opcionespararesumeo = trim(adjustl(l%opcionespararesumeo))//' '//trim(adjustl(l%chain))
!!!!#endif
!

            CASE ('-verbose')
               l%verbose = .TRUE.
            CASE ('-ignoreerrors')
               l%ignoreerrors = .TRUE.
            CASE ('-r')
               l%resume = .TRUE.
               l%forcesteps = .true.
#ifdef CompileWithOldSaving
            CASE ('-old')
               l%resume_fromold = .TRUE.
#endif
            CASE ('-cpumax')
               i = i + 1
               call getcommandargument(l%chaininput, i, f, l%length, statuse, binaryPath)
               READ (f, *, iostat=iostatus) l%maxCPUtime
               if (iostatus /= 0) call stoponerror(l%layoutnumber, l%num_procs, 'Invalid CPU maximum time', .true.)
               if (l%maxCPUtime <= 0) then
                  call stoponerror(l%layoutnumber, l%num_procs, 'Invalid CPU maximum time', .true.)
                  statuse = -1
               end if

            CASE ('-s')
               l%freshstart = .TRUE.
            CASE ('-flush')
               i = i + 1
               call getcommandargument(l%chaininput, i, f, l%length, statuse, binaryPath)
               READ (f, *, iostat=iostatus) l%flushminutesFields
               if (iostatus /= 0) call stoponerror(l%layoutnumber, l%num_procs, 'Invalid flushing interval', .true.)
               if (l%flushminutesFields <= 0) then
                  call stoponerror(l%layoutnumber, l%num_procs, 'Invalid flushing interval', .true.)
                  statuse = -1
               end if
            CASE ('-flushdata')
               i = i + 1
               call getcommandargument(l%chaininput, i, f, l%length, statuse, binaryPath)
               READ (f, *, iostat=iostatus) l%flushminutesData
               if (iostatus /= 0) call stoponerror(l%layoutnumber, l%num_procs, 'Invalid flushing interval', .true.)
401            if (l%flushminutesData <= 0) then
                  call stoponerror(l%layoutnumber, l%num_procs, 'Invalid flushing interval', .true.)
                  statuse = -1
               end if
            CASE ('-run')
               l%run = .TRUE.
            CASE ('-map')
               l%createmap = .TRUE.
            CASE ('-dontwritevtk')
               l%dontwritevtk = .true.
            CASE ('-vtkindex')
               l%vtkindex = .TRUE.
            CASE ('-mapvtk')
               l%createmapvtk = .TRUE.
            CASE ('-dmma')
               l%run_with_dmma = .TRUE.
               l%run_with_abrezanjas = .FALSE.
               l%opcionespararesumeo = trim(adjustl(l%opcionespararesumeo))//' '//trim(adjustl(l%chain))

#ifdef CompileWithConformal
            CASE ('-abrezanjas') !Provisional FEB-2018

               !NOTE: Comento lo de abajo para debugear
               !   print *,'Abrezanjas not available (290521).  '
               !   stop

               l%run_with_dmma = .FALSE.
               l%run_with_abrezanjas = .true.
               if (.NOT. l%input_conformal_flag) then
                  l%conformal_file_input_name = char(0)
                  l%input_conformal_flag = .true.
               end if
               l%opcionespararesumeo = trim(adjustl(l%opcionespararesumeo))//' '//trim(adjustl(l%chain))

            CASE ('-activateconf') !Provisional FEB-2018
               if (.NOT. l%input_conformal_flag) then
                  l%conformal_file_input_name = char(0)
                  l%input_conformal_flag = .true.
               end if
               i = i + 1; 
               l%opcionespararesumeo = trim(adjustl(l%opcionespararesumeo))//' '//trim(adjustl(l%chain))
#endif
            CASE ('-conf')
               l%flag_conf_sgg = .true.
               i = i + 1; 
#ifdef CompileWithConformal
               l%input_conformal_flag = .true.; 
               l%conformal_file_input_name = char(0); 
               call getcommandargument(l%chaininput, i, f, l%length, statuse, binaryPath)
               l%conformal_file_input_name = trim(adjustl(f)); 
               inquire(file=trim(adjustl(f)), EXIST=l%existeNFDE)
               if (.NOT. l%existeNFDE) then
                  l%input_conformal_flag = .FALSE.; 
                  buff = 'The conformal input file was not found '//trim(adjustl(l%fichin)); 
                  call WarnErrReport(Trim(buff), .true.)
               end if
               l%opcionespararesumeo = trim(adjustl(l%opcionespararesumeo))//' '//trim(adjustl(l%chain))
#endif
#ifndef CompileWithConformal
               if (l%input_conformal_flag) then
                  buff = ''; buff = 'Conformal without conformal support. Recompile!'; 
                  call WarnErrReport(Trim(buff), .true.)
               end if
#endif

            CASE ('-takeintcripte')
               l%takeintcripte = .true.
               l%opcionespararesumeo = trim(adjustl(l%opcionespararesumeo))//' '//trim(adjustl(l%chain))
#ifdef CompileWithNIBC
            CASE ('-skindepthpre')
               l%skindepthpre = .true.
!            l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain))
            CASE ('-mibc', '-skindepth')
               l%mibc = .true.
               l%sgbc = .false.
               l%opcionespararesumeo = trim(adjustl(l%opcionespararesumeo))//' '//trim(adjustl(l%chain))
            CASE ('-conformalskin')
               l%conformalskin = .true.
               l%mibc = .true.
               l%sgbc = .false.
               l%opcionespararesumeo = trim(adjustl(l%opcionespararesumeo))//' '//trim(adjustl(l%chain))
            CASE ('-ade')
               l%ade = .true.
               l%mibc = .true.
               l%sgbc = .false.
               l%opcionespararesumeo = trim(adjustl(l%opcionespararesumeo))//' '//trim(adjustl(l%chain))
            CASE ('-NOcompomur')
               l%NOcompomur = .true.
               l%mibc = .true.
               l%sgbc = .false.
               l%opcionespararesumeo = trim(adjustl(l%opcionespararesumeo))//' '//trim(adjustl(l%chain))
#endif
            CASE ('-mur2')
               l%MurAfterPML = .true.
               !l%mur_second=.true.
               l%mur_first = .true.
               !arreglar cuando resuelva el bug en mur segundo orden
               l%opcionespararesumeo = trim(adjustl(l%opcionespararesumeo))//' '//trim(adjustl(l%chain))
            CASE ('-mur1')
               l%MurAfterPML = .true.
               l%mur_first = .true.
               l%mur_second = .false.
               l%opcionespararesumeo = trim(adjustl(l%opcionespararesumeo))//' '//trim(adjustl(l%chain))
            CASE ('-pmlalpha')
               i = i + 1
               call getcommandargument(l%chaininput, i, f, l%length, statuse, binaryPath)
               ! Converts the characters to real
               READ (f, *, ERR=7621) l%alphamaxpar
               GO TO 8621
7621           call stoponerror(l%layoutnumber, l%num_procs, 'Invalid CPML alpha factor', .true.)
               statuse = -1
               !goto 668
8621           if (l%alphamaxpar < 0.0_RKIND) then
                  call stoponerror(l%layoutnumber, l%num_procs, 'Invalid CPML alpha factor', .true.)
                  statuse = -1
                  !goto 668
               end if
               i = i + 1
               !          l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain))// ' ' // trim (adjustl(f))
               call getcommandargument(l%chaininput, i, f, l%length, statuse, binaryPath)
               ! Converts the characters to real
               READ (f, *, ERR=7121) l%alphaOrden
               GO TO 8121
7121           call stoponerror(l%layoutnumber, l%num_procs, 'Invalid CPML order factor', .true.)
               statuse = -1
               !goto 668
8121           if (l%alphaOrden < 0.0_RKIND) then
                  call stoponerror(l%layoutnumber, l%num_procs, 'Invalid CPML alpha factor', .true.)
                  statuse = -1
                  !goto 668
               end if
               !          l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain))// ' ' // trim (adjustl(f))
            CASE ('-pmlkappa')
               i = i + 1
               call getcommandargument(l%chaininput, i, f, l%length, statuse, binaryPath)
               ! Converts the characters to real
               READ (f, *, ERR=7622) l%kappamaxpar
               GO TO 8622
7622           call stoponerror(l%layoutnumber, l%num_procs, 'Invalid CPML kappa factor', .true.)
               statuse = -1
               !goto 668
8622           if (l%kappamaxpar < 1.0_RKIND) then
                  call stoponerror(l%layoutnumber, l%num_procs, 'Invalid CPML kappa factor', .true.)
                  statuse = -1
                  !goto 668
               end if
               !          l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain))// ' ' // trim (adjustl(f))
            CASE ('-pmlcorr')
               l%MEDIOEXTRA%exists = .true.
               i = i + 1
               call getcommandargument(l%chaininput, i, f, l%length, statuse, binaryPath)
               ! Converts the characters to real
               READ (f, *, ERR=7672) l%MEDIOEXTRA%sigma
               GO TO 8672
7672           call stoponerror(l%layoutnumber, l%num_procs, 'Invalid pmlcorr sigma factor', .true.)
               statuse = -1
               !goto 668
8672           if (l%MEDIOEXTRA%sigma < 0.0_RKIND) then
                  call stoponerror(l%layoutnumber, l%num_procs, 'Invalid pmlcorr sigma factor', .true.)
                  statuse = -1
                  !goto 668
               end if
               l%MEDIOEXTRA%sigmam = -1.0_RKIND!voids it. later overriden
               !          l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain))// ' ' // trim (adjustl(f))
               i = i + 1
               call getcommandargument(l%chaininput, i, f, l%length, statuse, binaryPath)
               ! Converts the characters to real
               READ (f, *, ERR=7662) l%MEDIOEXTRA%pml_size
               GO TO 8662
7662           call stoponerror(l%layoutnumber, l%num_procs, 'Invalid pmlcorr depth factor', .true.)
               statuse = -1
               !goto 668
8662           if (l%MEDIOEXTRA%pml_size < 0) then
                  call stoponerror(l%layoutnumber, l%num_procs, 'Invalid pmlcorr depth factor', .true.); statuse = -1; !goto 668
               end if
               !          l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain))// ' ' // trim (adjustl(f))
            CASE ('-attc')
               i = i + 1
               call getcommandargument(l%chaininput, i, f, l%length, statuse, binaryPath)
               ! Converts the characters to real
               READ (f, *, ERR=766) l%attfactorc
               GO TO 866
766            call stoponerror(l%layoutnumber, l%num_procs, 'Invalid dissipation factor', .true.); statuse = -1; !goto 668
866            if ((l%attfactorc <= -1.0_RKIND) .or. (l%attfactorc > 1.0_RKIND)) then
                  call stoponerror(l%layoutnumber, l%num_procs, 'Invalid dissipation factor', .true.); statuse = -1; !goto 668
               end if
               l%opcionespararesumeo = trim(adjustl(l%opcionespararesumeo))//' '//trim(adjustl(l%chain))//' '//trim(adjustl(f))
            CASE ('-sgbcdepth')
               l%mibc = .false.
               l%sgbc = .true.
               i = i + 1
               call getcommandargument(l%chaininput, i, f, l%length, statuse, binaryPath)
               ! Converts the characters to real
               READ (f, *, ERR=7466) l%sgbcdepth
               GO TO 8466
7466           call stoponerror(l%layoutnumber, l%num_procs, 'Invalid sgbc depth ', .true.); statuse = -1; !goto 668
8466           if (l%sgbcdepth < -1) then
                  call stoponerror(l%layoutnumber, l%num_procs, 'Invalid sgbc depth', .true.); statuse = -1; !goto 668
               end if
               l%opcionespararesumeo = trim(adjustl(l%opcionespararesumeo))//' '//trim(adjustl(l%chain))//' '//trim(adjustl(f))
            CASE ('-sgbcfreq')
               l%sgbc = .true.
               l%mibc = .false.
               i = i + 1
               call getcommandargument(l%chaininput, i, f, l%length, statuse, binaryPath)
               ! Converts the characters to real
               READ (f, *, ERR=74616) l%sgbcfreq
               GO TO 84616
74616          call stoponerror(l%layoutnumber, l%num_procs, 'Invalid sgbc freq ', .true.); statuse = -1; !goto 668
84616          if (l%sgbcfreq < 0.) then
                  call stoponerror(l%layoutnumber, l%num_procs, 'Invalid sgbc freq', .true.); statuse = -1; !goto 668
               end if
               l%opcionespararesumeo = trim(adjustl(l%opcionespararesumeo))//' '//trim(adjustl(l%chain))//' '//trim(adjustl(f))
            CASE ('-sgbcresol')
               l%mibc = .false.
               l%sgbc = .true.
               i = i + 1
               call getcommandargument(l%chaininput, i, f, l%length, statuse, binaryPath)
               ! Converts the characters to real
               READ (f, *, ERR=74626) l%sgbcresol
               GO TO 84626
74626          call stoponerror(l%layoutnumber, l%num_procs, 'Invalid sgbc decay ', .true.); statuse = -1; !goto 668
84626          if (l%sgbcresol < 0.0) then
                  call stoponerror(l%layoutnumber, l%num_procs, 'Invalid sgbc decay', .true.); statuse = -1; !goto 668
               end if
               l%opcionespararesumeo = trim(adjustl(l%opcionespararesumeo))//' '//trim(adjustl(l%chain))//' '//trim(adjustl(f))
            CASE ('-sgbcyee')
               l%sgbc = .true.
               l%mibc = .false.
               l%sgbccrank = .false.
               l%opcionespararesumeo = trim(adjustl(l%opcionespararesumeo))//' '//trim(adjustl(l%chain))
            CASE ('-sgbccrank') !es el default. Lo mantengo por compatibilidad con lanzamientos previos
               l%sgbccrank = .true.
               l%mibc = .false.
               l%opcionespararesumeo = trim(adjustl(l%opcionespararesumeo))//' '//trim(adjustl(l%chain))
            CASE ('-nosgbc') !opcion generica que aglutina varios switches que estan den default (l%sgbcresol, l%sgbccrank, l%sgbcfreq)
               l%sgbc = .false.
               l%mibc = .true.
               l%opcionespararesumeo = trim(adjustl(l%opcionespararesumeo))//' '//trim(adjustl(l%chain))
            CASE ('-sgbc') !opcion generica que aglutina varios switches que estan den default (l%sgbcresol, l%sgbccrank, l%sgbcfreq)
               l%sgbc = .true.
               l%mibc = .false.
               l%opcionespararesumeo = trim(adjustl(l%opcionespararesumeo))//' '//trim(adjustl(l%chain))
            CASE ('-sgbcDispersive') !opcion generica que aglutina varios switches que estan den default (l%sgbcresol, l%sgbccrank, l%sgbcfreq)
               l%sgbc = .true.
               l%mibc = .false.
               l%sgbcDispersive = .true.
               l%opcionespararesumeo = trim(adjustl(l%opcionespararesumeo))//' '//trim(adjustl(l%chain))
            CASE ('-saveall')
               l%saveall = .TRUE.
            CASE ('-attw')
               i = i + 1
               call getcommandargument(l%chaininput, i, f, l%length, statuse, binaryPath)
               ! Converts the characters to real
               READ (f, *, ERR=732) l%attfactorw
               GO TO 832
732            call stoponerror(l%layoutnumber, l%num_procs, 'Invalid dissipation factor', .true.); statuse = -1; !goto 668
832            if ((l%attfactorw <= -1.0_RKIND) .or. (l%attfactorw > 1.0_RKIND)) then
                  call stoponerror(l%layoutnumber, l%num_procs, 'Invalid dissipation factor', .true.); statuse = -1; !goto 668
               end if
               l%opcionespararesumeo = trim(adjustl(l%opcionespararesumeo))//' '//trim(adjustl(l%chain))//' '//trim(adjustl(f))
            CASE ('-maxwireradius')
               l%boundwireradius = .true.
               i = i + 1
               call getcommandargument(l%chaininput, i, f, l%length, statuse, binaryPath)
               ! Converts the characters to real
               READ (f, *, ERR=737) l%maxwireradius
               GO TO 837
737            call stoponerror(l%layoutnumber, l%num_procs, 'Invalid dissipation factor', .true.); statuse = -1; !goto 668
837            if ((l%maxwireradius <= 0.0_RKIND)) then
                  call stoponerror(l%layoutnumber, l%num_procs, 'Invalid maximumwireradius', .true.); statuse = -1; !goto 668
               end if
               l%opcionespararesumeo = trim(adjustl(l%opcionespararesumeo))//' '//trim(adjustl(l%chain))//' '//trim(adjustl(f))
            CASE ('-mindistwires')
               i = i + 1
               call getcommandargument(l%chaininput, i, f, l%length, statuse, binaryPath)
               ! Converts the characters to real
               READ (f, *, ERR=1732) l%mindistwires
               GO TO 1832
1732           call stoponerror(l%layoutnumber, l%num_procs, 'Invalid minimum distance between wires', .true.); statuse = -1; !goto 668
1832           if (l%mindistwires <= 0.0_RKIND) then
                  call stoponerror(l%layoutnumber, l%num_procs, 'Invalid minimum distance between wires', .true.); statuse = -1; !goto 668
               end if
               l%opcionespararesumeo = trim(adjustl(l%opcionespararesumeo))//' '//trim(adjustl(l%chain))//' '//trim(adjustl(f))
            CASE ('-makeholes')
               l%makeholes = .TRUE.
               l%opcionespararesumeo = trim(adjustl(l%opcionespararesumeo))//' '//trim(adjustl(l%chain))
            CASE ('-connectendings')
               l%connectendings = .TRUE.
               l%opcionespararesumeo = trim(adjustl(l%opcionespararesumeo))//' '//trim(adjustl(l%chain))
            CASE ('-nostrictOLD')
               l%strictOLD = .false.
               l%opcionespararesumeo = trim(adjustl(l%opcionespararesumeo))//' '//trim(adjustl(l%chain))
            CASE ('-nomtlnberenger')
               l%mtlnberenger = .false.
               l%opcionespararesumeo = trim(adjustl(l%opcionespararesumeo))//' '//trim(adjustl(l%chain))
            CASE ('-stableradholland')
               l%stableradholland = .true.
               l%opcionespararesumeo = trim(adjustl(l%opcionespararesumeo))//' '//trim(adjustl(l%chain))
               !  CASE ('-mtln')
               !      buff='-mtln option deprecated and ignored. Check -nomtlnberenger or -l%stableradholland'
               !      call WarnErrReport(Trim(buff),.false.)
            CASE ('-intrawiresimplify')
               l%strictOLD = .false.
               l%opcionespararesumeo = trim(adjustl(l%opcionespararesumeo))//' '//trim(adjustl(l%chain))
            CASE ('-notaparrabos')
               l%TAPARRABOS = .false.
               l%opcionespararesumeo = trim(adjustl(l%opcionespararesumeo))//' '//trim(adjustl(l%chain))
            case ('-fieldtotl')
               l%fieldtotl = .true.
               l%opcionespararesumeo = trim(adjustl(l%opcionespararesumeo))//' '//trim(adjustl(l%chain))
          !!case ('-experimentalVideal')
          !!    l%experimentalVideal=.true.
          !!    l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain))

            case ('-forceresampled') !a menos que se pida explicitamente, no se resamplea 120123
               l%forceresampled = .true.
               l%opcionespararesumeo = trim(adjustl(l%opcionespararesumeo))//' '//trim(adjustl(l%chain))

            CASE ('-wirethickness')
               i = i + 1
               call getcommandargument(l%chaininput, i, f, l%length, statuse, binaryPath)
               ! Converts the characters to real
               READ (f, *, ERR=7416) l%wirethickness
               GO TO 8416
7416           call stoponerror(l%layoutnumber, l%num_procs, 'Invalid l%wirethickness ', .true.); statuse = -1; !goto 668
8416           if (l%sgbcdepth < -1) then
                  call stoponerror(l%layoutnumber, l%num_procs, 'Invalid l%wirethickness', .true.); statuse = -1; !goto 668
               end if
               l%opcionespararesumeo = trim(adjustl(l%opcionespararesumeo))//' '//trim(adjustl(l%chain))//' '//trim(adjustl(f))
            CASE ('-wiresflavor')
               i = i + 1
               call getcommandargument(l%chaininput, i, f, l%length, statuse, binaryPath)
               l%opcionespararesumeo = trim(adjustl(l%opcionespararesumeo))//' '//trim(adjustl(l%chain))//' '//trim(adjustl(f))
               READ (f, '(a)', ERR=3621) l%wiresflavor
               if (trim(adjustl(l%wiresflavor(1:1))) == 'g') l%wiresflavor = 'slanted'
               select case (trim(adjustl(l%wiresflavor)))
               case ('holland', 'old')
                  l%wiresflavor = 'holland'
               case ('berenger', 'new')
                  l%wiresflavor = 'berenger'
               case ('slanted', 'experimental')
                  l%wiresflavor = 'slanted'
               case ('transition')
                  l%wiresflavor = 'transition'
               case ('semistructured')
                  l%wiresflavor = 'semistructured'
                  !
                  i = i + 1
                  call getcommandargument(l%chaininput, i, f, l%length, statuse, binaryPath)
                  l%opcionespararesumeo = trim(adjustl(l%opcionespararesumeo))//' '//trim(adjustl(f))
                  ! Converts the characters to real
                  READ (f, *, ERR=2561) l%precision
                  GO TO 2562
2561              call stoponerror(l%layoutnumber, l%num_procs, 'Invalid l%precision for semistructured', .true.); statuse = -1; !goto 668
2562              if (l%precision < 0) then
                     call stoponerror(l%layoutnumber, l%num_procs, 'Invalid l%precision for semistructured', .true.); statuse = -1; !goto 668
                  end if
                  !
               end select
               GO TO 4621
3621           call stoponerror(l%layoutnumber, l%num_procs, 'Invalid wires flavor', .true.); statuse = -1; !goto 668
4621           if (((trim(adjustl(l%wiresflavor)) /= 'holland') .AND. &
                    (trim(adjustl(l%wiresflavor)) /= 'transition') .AND. &
                    (trim(adjustl(l%wiresflavor)) /= 'berenger') .AND. &
                    (trim(adjustl(l%wiresflavor)) /= 'slanted') .and. &
                    (trim(adjustl(l%wiresflavor)) /= 'semistructured')) .or. &
                   .not. ((trim(adjustl(l%wiresflavor)) == 'holland') .xor. &
                          (trim(adjustl(l%wiresflavor)) == 'transition') .xor. &
                          (trim(adjustl(l%wiresflavor)) == 'berenger') .xor. &
                          (trim(adjustl(l%wiresflavor)) == 'slanted') .xor. &
                          (trim(adjustl(l%wiresflavor)) == 'semistructured'))) then
                  call stoponerror(l%layoutnumber, l%num_procs, 'Invalid wires flavor->'//trim(adjustl(l%wiresflavor)), .true.); statuse = -1; !goto 668
               end if
#ifndef CompileWithThickWires
               select case (trim(adjustl(l%wiresflavor)))
               case ('holland', 'transition')
                  if (l%wirethickness /= 1) then
                     call stoponerror(l%layoutnumber, l%num_procs, 'Holland wire flavor not available in this compilation', .true.); statuse = -1; !goto 668
                  end if
               end select
#endif
#ifndef CompileWithThickWires
               select case (trim(adjustl(l%wiresflavor)))
               case ('holland')
                  if (l%wirethickness /= 1) then
                     call stoponerror(l%layoutnumber, l%num_procs, 'Holland wire flavor thickness>1 requires recompiling', .true.); statuse = -1; !goto 668
                  end if
               end select
#endif
               select case (trim(adjustl(l%wiresflavor)))
               case ('berenger', 'slanted', 'experimental', 'transition')
                  if (l%wirethickness /= 1) then
                     call stoponerror(l%layoutnumber, l%num_procs, 'Thickness>1 unsupported for this wireflavor', .true.); statuse = -1; !goto 668
                  end if
               end select
#ifndef CompileWithBerengerWires
               select case (trim(adjustl(l%wiresflavor)))
               case ('berenger')
                  call stoponerror(l%layoutnumber, l%num_procs, 'Berenger wire flavor not available in this compilation', .true.); statuse = -1; !goto 668
               end select
#endif
#ifndef CompileWithSlantedWires
               select case (trim(adjustl(l%wiresflavor)))
               case ('slanted', 'experimental')
                  call stoponerror(l%layoutnumber, l%num_procs, 'Experimental wire flavor not available in this compilation', .true.); statuse = -1; !goto 668
               end select
#endif
            CASE ('-isolategroupgroups')
               l%isolategroupgroups = .TRUE.
               l%opcionespararesumeo = trim(adjustl(l%opcionespararesumeo))//' '//trim(adjustl(l%chain))
            CASE ('-groundwires')
               l%groundwires = .TRUE.
               l%opcionespararesumeo = trim(adjustl(l%opcionespararesumeo))//' '//trim(adjustl(l%chain))
            CASE ('-noSlantedcrecepelo ') !opcion niapa excperimental 131219
               l%noSlantedcrecepelo = .TRUE.
               l%opcionespararesumeo = trim(adjustl(l%opcionespararesumeo))//' '//trim(adjustl(l%chain))
            CASE ('-inductance')
               i = i + 1
               call getcommandargument(l%chaininput, i, f, l%length, statuse, binaryPath)
               READ (f, '(a)', ERR=361) l%inductance_model
               GO TO 461
361            call stoponerror(l%layoutnumber, l%num_procs, 'Invalid inductance model', .true.); statuse = -1; !goto 668
461            if ((l%inductance_model /= 'ledfelt') .AND. (l%inductance_model /= 'berenger') .AND. &
                       &    (l%inductance_model /= 'boutayeb')) then
                  call stoponerror(l%layoutnumber, l%num_procs, 'Invalid inductance model', .true.); statuse = -1; !goto 668
               end if
               l%opcionespararesumeo = trim(adjustl(l%opcionespararesumeo))//' '//trim(adjustl(l%chain))//' '//trim(adjustl(f))
            CASE ('-inductanceorder')
               i = i + 1
               call getcommandargument(l%chaininput, i, f, l%length, statuse, binaryPath)
               READ (f, *, ERR=179) l%inductance_order
               GO TO 180
179            call stoponerror(l%layoutnumber, l%num_procs, 'Invalid inductance order', .true.); statuse = -1; !goto 668
180            l%opcionespararesumeo = trim(adjustl(l%opcionespararesumeo))//' '//trim(adjustl(l%chain))//' '//trim(adjustl(f))
            CASE ('-prefix')
               i = i + 1
               call getcommandargument(l%chaininput, i, f, l%length, statuse, binaryPath)
               l%prefix = '_'//trim(adjustl(f))
               l%opcionespararesumeo = trim(adjustl(l%opcionespararesumeo))//' '//trim(adjustl(l%chain))//' '//trim(adjustl(f))
            CASE ('-cfl')
               i = i + 1
               call getcommandargument(l%chaininput, i, f, l%length, statuse, binaryPath)
               ! Converts the characters to real
               READ (f, *, ERR=3762) l%cfltemp
               GO TO 3862
3762           call stoponerror(l%layoutnumber, l%num_procs, 'Invalid Courant Number', .true.); statuse = -1; !goto 668
3862           if (l%cfltemp <= 0.0) then
                  call print11(l%layoutnumber, '------> Ignoring negative or null l%cfl Courant Number')
!!!!!!!!!               call stoponerror (l%layoutnumber, l%num_procs, 'Invalid negative or null l%cfl Courant Number',.true.); statuse=-1; !goto 668 !!!sgg 151216 para evitar el error l%cfl 0 del problem-type sigue como si no estuviera
                  l%forcecfl = .false.
               else
                  l%cfl = l%cfltemp
                  l%forcecfl = .true.
                  l%opcionespararesumeo = trim(adjustl(l%opcionespararesumeo))//' '//trim(adjustl(l%chain))//' '//trim(adjustl(f))
               end if
            CASE ('-noconformalmapvtk')
               l%noconformalmapvtk = .true.
            CASE ('-niapapostprocess')
               l%niapapostprocess = .true.
#ifdef CompileWithPrescale
!!!!210918 permit scaling
            CASE ('-pscale')
               l%permitscaling = .true.
               l%saveall = .true. !lo salvo todo en permit scaling para evitar errores
               i = i + 1
               buff = ""
               call getcommandargument(l%chaininput, i, f, l%length, statuse, binaryPath)
               READ (f, *, ERR=33762) buff
               l%EpsMuTimeScale_input_parameters%electric = .False.
               l%EpsMuTimeScale_input_parameters%electric = .False.
               Select case (trim(adjustl(buff)))
               case ("ee")
                  l%EpsMuTimeScale_input_parameters%electric = .True.
               case ("hh")
                  l%EpsMuTimeScale_input_parameters%magnetic = .True.
               case ("eh":"he")
                  l%EpsMuTimeScale_input_parameters%electric = .True.
                  l%EpsMuTimeScale_input_parameters%magnetic = .True.
               case default
                  GO TO 33862
               end select
               i = i + 1
               call getcommandargument(l%chaininput, i, f, l%length, statuse, binaryPath)
               l%opcionespararesumeo = trim(adjustl(l%opcionespararesumeo))//' '//trim(adjustl(l%chain))//' '//trim(adjustl(f))
               ! Converts the characters to real
               READ (f, *, ERR=33762) l%EpsMuTimeScale_input_parameters%tini
               i = i + 1
               call getcommandargument(l%chaininput, i, f, l%length, statuse, binaryPath)
               l%opcionespararesumeo = trim(adjustl(l%opcionespararesumeo))//' '//trim(adjustl(f))
               READ (f, *, ERR=33762) l%EpsMuTimeScale_input_parameters%tend
               i = i + 1
               call getcommandargument(l%chaininput, i, f, l%length, statuse, binaryPath)
               l%opcionespararesumeo = trim(adjustl(l%opcionespararesumeo))//' '//trim(adjustl(f))
               READ (f, *, ERR=33762) l%EpsMuTimeScale_input_parameters%alpha_max
               GO TO 33862
33762          call stoponerror(l%layoutnumber, l%num_procs, 'Invalid pscale parameters', .true.); statuse = -1; !goto 668
33862          continue
               if (l%EpsMuTimeScale_input_parameters%checkError() /= 0) then
                  call stoponerror(l%layoutnumber, l%num_procs, &
     &'Invalid -pscale parameters: some parameters have to be greater than 0.0: -pscale t0(>=0) tend slope(>0)'&
                    &, .true.); statuse = -1; !goto 668
               else
                  l%EpsMuTimeScale_input_parameters%are_there = .true.
               end if
#endif
            CASE ('-n')
               l%forcesteps = .TRUE.
               i = i + 1
               call getcommandargument(l%chaininput, i, f, l%length, statuse, binaryPath)
               ! Converts the characters to integer
               READ (f, *, ERR=602) l%finaltimestep
               GO TO 702
602            call stoponerror(l%layoutnumber, l%num_procs, 'Invalid time step', .true.); statuse = -1; !goto 668
702            if (l%finaltimestep < -2) then
                  call stoponerror(l%layoutnumber, l%num_procs, 'Invalid time step', .true.); statuse = -1; !goto 668
               end if
!!!!!!
            CASE ('-factorradius')
               i = i + 1
               call getcommandargument(l%chaininput, i, f, l%length, statuse, binaryPath)
               ! Converts the characters to integer
               READ (f, *, ERR=6032) l%factorradius
               GO TO 7032
6032           call stoponerror(l%layoutnumber, l%num_procs, 'Invalid l%factorradius', .true.); statuse = -1; !goto 668
7032           continue
            CASE ('-factordelta')
               i = i + 1
               call getcommandargument(l%chaininput, i, f, l%length, statuse, binaryPath)
               ! Converts the characters to integer
               READ (f, *, ERR=6072) l%factordelta
               GO TO 7072
6072           call stoponerror(l%layoutnumber, l%num_procs, 'Invalid l%factordelta', .true.); statuse = -1; !goto 668
7072           continue
!!!!!!!!!!!!!
            CASE ('-stoch')
               l%stochastic = .true.
               l%chosenyesornostochastic = .true.
               l%opcionespararesumeo = trim(adjustl(l%opcionespararesumeo))//' '//trim(adjustl(l%chain))
#ifndef CompileWithMPI
               call stoponerror(l%layoutnumber, l%num_procs, 'l%stochastic simulation unsupported without MPI compilation', .true.); statuse = -1; !goto 668
#endif
            CASE ('-nostoch')
               l%stochastic = .false.
               l%chosenyesornostochastic = .true.
               l%opcionespararesumeo = trim(adjustl(l%opcionespararesumeo))//' '//trim(adjustl(l%chain))
            case ('-forcecreateh5bin')
               l%createh5bin = .true.
            CASE ('') !100615 para evitar el crlf del .sh
               continue
            CASE DEFAULT
               call stoponerror(l%layoutnumber, l%num_procs, 'Wrong switch '//trim(adjustl(l%chain)), .true.); statuse = -1; !goto 668
            end select
            i = i + 1
         end do

      end if
      !some checkings
      !just to be sure that I do not have stupid errors
#ifdef CompileWithMPI
      if ((sizeof(MPI_DOUBLE_PRECISION) /= 4) .OR. (sizeof(MPI_REAL) /= 4)) then
         call stoponerror(l%layoutnumber, l%num_procs, 'SEVERE COMPILATION ERROR: MPI_REAL l%num_procs is not 4. ')
      end if
#endif

      if (l%connectendings .AND. l%strictOLD) then
         call stoponerror(l%layoutnumber, l%num_procs, 'l%strictOLD option not compatible with -l%connectendings', .true.); statuse = -1; !goto 668
      end if
      if (l%TAPARRABOS .AND. (.not. l%strictOLD)) then
         call stoponerror(l%layoutnumber, l%num_procs, '-nostrictOLD option requires -notaparrabos ', .true.); statuse = -1; !goto 668
      end if
      if (l%isolategroupgroups .AND. l%strictOLD) then
         call stoponerror(l%layoutnumber, l%num_procs, '-intrawiresimplify option not compatible with -l%isolategroupgroups', .true.); statuse = -1; !goto 668
      end if

      if ((l%sgbc .AND. l%mibc)) then
         call stoponerror(l%layoutnumber, l%num_procs, 'Use only one of -sgbc -l%mibc', .true.); statuse = -1; !goto 668
      end if
      if (l%freshstart .AND. l%resume) then
         call stoponerror(l%layoutnumber, l%num_procs, 'Fresh Start option -s not compatible with restarting -r', .true.); statuse = -1; !goto 668
      end if
      if (l%freshstart .AND. l%resume_fromold) then
         call stoponerror(l%layoutnumber, l%num_procs, 'Fresh Start option -s not compatible with -old', .true.); statuse = -1; !goto 668
      end if
      if ((.NOT. l%resume) .and. (.not. l%run) .AND. l%resume_fromold) then
         call stoponerror(l%layoutnumber, l%num_procs, 'l%resume option -r must be used if issuing -old', .true.); statuse = -1; !goto 668
      end if
      if ((l%flushminutesFields /= 0) .AND. (l%deleteintermediates)) then
         call stoponerror(l%layoutnumber, l%num_procs, '-delete is not compatible with -flush', .true.); statuse = -1; !goto 668
      end if
      if (l%run_with_abrezanjas .and. l%run_with_dmma) then
         call stoponerror(l%layoutnumber, l%num_procs, '-abrezanjas is not compatible with -dmma', .true.); statuse = -1; !goto 668
      end if
      if (l%stochastic .and. (trim(adjustl(l%wiresflavor)) /= 'holland')) then
         call stoponerror(l%layoutnumber, l%num_procs, 'Old wires flavor is the only supported with l%stochastic', .true.); statuse = -1; !goto 668
      end if
      if (l%stochastic .and. l%wirecrank) then
         call stoponerror(l%layoutnumber, l%num_procs, 'Wires Crank Nicolson is unsupported with l%stochastic', .true.); statuse = -1; !goto 668
      end if
   !!!si esta soportado 170719
   !! if (l%permitscaling.and.l%resume) then
   !!   call stoponerror (l%layoutnumber, l%num_procs, 'Resuming with Permittivity scaling unsupported',.true.); statuse=-1; !goto 668
   !!end if
      if (l%permitscaling .and. (l%kappamaxpar .gt. 1.000001_rkind)) then
   !!!061118 no lo permito porque cpml toca los idxe, idye, idze en funcion del kappa y permittivity scaling conflicta
            call stoponerror (l%layoutnumber, l%num_procs, 'Unsupported CPML kappa factor since 061118 because conflicts with Idxe...in permittivity scaling',.true.)
      end if
      if (l%stochastic) then
#ifndef CompileWithStochastic
         call StopOnError(l%layoutnumber, l%num_procs, 'l%stochastic without compilation support. Recompile')
#endif
#ifdef CompileWithStochastic
#ifndef CompileWithMPI
         call StopOnError(l%layoutnumber, l%num_procs, 'l%stochastic unsupported without MPI compilation. Recompile')
#endif
#endif
         continue
      end if

!!!

      !
      l%prefixopci1 = trim(adjustl(l%opcionespararesumeo))
      l%prefixopci = ' '
      do i = 1, len(trim(adjustl(l%prefixopci1)))
         l%prefixopci(i:i) = l%prefixopci1(i:i)
         j = ichar(l%prefixopci1(i:i))
         if (j <= 47) l%prefixopci(i:i) = '_'
         if (j >= 123) l%prefixopci(i:i) = '_'
         if ((j >= 58) .and. (j <= 64)) l%prefixopci(i:i) = '_'
         if ((j >= 91) .and. (j <= 96)) l%prefixopci(i:i) = '_'
         if (j == 46) l%prefixopci(i:i) = 'p'
      end do

      do i = 1, len(trim(adjustl(l%prefixopci)))
         do while (l%prefixopci(i:i + 1) == '__')
            l%prefixopci(i:) = l%prefixopci(i + 1:)
         end do
      end do
      if (l%prefix(1:1) == '_') then
     !!!acortado 120219  l%nEntradaRoot = trim (adjustl(l%fichin)) // trim (adjustl(prefix))// trim (adjustl(l%prefixopci))
         l%nEntradaRoot = trim(adjustl(l%fichin))//'_'//trim(adjustl(l%prefixopci))
      else
         l%nEntradaRoot = trim(adjustl(l%fichin))
      end if
!!!l%stochastic
#ifdef CompileWithStochastic
      if (l%stochastic) then
         if (l%layoutnumber <= l%num_procs/2 - 1) then !aun no se ha dividido el l%num_procs
            l%nEntradaRoot = trim(adjustl(l%nEntradaRoot))
         else
            l%nEntradaRoot = trim(adjustl('devia_'//trim(adjustl(l%nEntradaRoot))))
         end if
      end if
#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI, l%ierr)
#endif
#endif
!!!fin l%stochastic
!!!   sgg%nEntradaRoot=trim (adjustl(l%nEntradaRoot))
      !
      write(chari, '(i5)') l%layoutnumber + 1
      l%nresumeable2 = trim(adjustl(l%nEntradaRoot))//'_'//trim(adjustl(chari))//'.fields'
      !

      l%geomfile = trim(adjustl(l%nEntradaRoot))//'_'//trim(adjustl(chari))
      !warning file management
      if (statuse /= -1) then
         call CLOSEWARNINGFILE(l%layoutnumber, l%num_procs, l%fatalerror, .false., .false.) !!cierra el temporal !todavia no se ha dividido el l%num_procs
       !!!borra el tmp si no hay l%fatalerror y reabre el de verdad
         if ((.not. l%fatalerror) .and. (l%layoutnumber == 0)) then
            open (unit=1320, file=trim(adjustl(l%fichin))//'_tmpWarnings.txt_Warnings.txt')
            close (unit=1320, status='delete')
         end if
         call INITWARNINGFILE(l%layoutnumber, l%num_procs, l%nEntradaRoot, l%verbose, l%ignoreerrors)
      end if
      !
      !
      if (l%resume_fromold) then
         inquire(file=trim(adjustl(l%nresumeable2))//'.old', EXIST=resume3)
      ELSE
         inquire(file=trim(adjustl(l%nresumeable2)), EXIST=resume3)
      end if
      if (l%resume) then
         if (.NOT. resume3) then
            call stoponerror(l%layoutnumber, l%num_procs, 'l%resume fields not present', .true.); statuse = -1; !goto 668
         end if
         write(dubuf, *) 'RESUMING simulation ', trim(adjustl(l%nEntradaRoot)), ' until n= ', l%finaltimestep
         call print11(l%layoutnumber, dubuf)
      ELSE
         if (resume3 .AND. (.NOT. l%freshstart) .and. (.not. l%run)) then
         call stoponerror (l%layoutnumber, l%num_procs, 'Restarting file exists. Either specify -r to l%resume, -s to do a fresh START, or -run to run in whatever the case',.true.); statuse = -1; !goto 668
         ELSEIF (resume3 .and. (l%run)) then
            l%resume = .true.
         ELSE
            open(35, file=trim(adjustl(l%nresumeable2)))
            write(35, '(a)') '!END'
            CLOSE (35, status='DELETE')
            open(35, file=trim(adjustl(l%nresumeable2))//'.old')
            write(35, '(a)') '!END'
            CLOSE (35, status='DELETE')
         end if
      end if
!
!
!
      if (((l%wiresflavor == 'slanted') .or. (l%wiresflavor == 'semistructured')) .AND. (l%mpidir /= 3)) then
         continue !arreglado l%mpidir slanted 2019
         !         call stoponerror (l%layoutnumber, l%num_procs, 'slanted wires unsupported with -l%mpidir {x,y}',.true.); statuse=-1; !goto 668
      end if
      if (l%input_conformal_flag .AND. (l%mpidir /= 3)) then
         continue !arreglado l%mpidir conformal 2019
         !TODO: under test
         !26-sep-2018: lo comento
         !call stoponerror (l%layoutnumber, l%num_procs, 'CONFORMAL -conf  unsupported with -l%mpidir {x,y}',.true.); statuse=-1; !goto 668
      end if
      if (l%run_with_abrezanjas .AND. (l%mpidir /= 3)) then
         continue !arreglado l%mpidir conformal 2019
         !under test
         !26-sep-2018: lo comento
         !call stoponerror (l%layoutnumber, l%num_procs, 'New abrezanjas thin gaps unsupported with -l%mpidir {x,y}',.true.); statuse=-1; !goto 668
      end if
      if (l%run_with_abrezanjas .AND. l%flag_conf_sgg) then
         !pass Mayo-2018
         !call stoponerror (l%layoutnumber, l%num_procs, 'CONFORMAL -conf currently unsupported with new abrezanjas thin gaps (unsupported 2 simultaneous conformal meshes at this moment',.true.); statuse=-1; !goto 668
         !se hace en otro sitio
      end if

      !
      if (((l%forcesteps) .AND. (.NOT. l%freshstart)) .and. (statuse /= -1)) then
         !in case of option -n withouth the l%freshstart option -s, it will l%resume or do a fresh start
         !depending on wether the resuming files are present or not
         if (l%resume_fromold) then
            inquire(file=trim(adjustl(l%nresumeable2))//'.old', EXIST=l%resume)
         ELSE
            inquire(file=trim(adjustl(l%nresumeable2)), EXIST=l%resume)
         end if
         if (l%resume) then
            if ((l%layoutnumber == 0) .or. ((l%layoutnumber == l%num_procs/2) .and. l%stochastic)) then
               !the temporary
               CLOSE (11)
               l%file11isopen = .false.
               open(11, file=trim(adjustl(l%nEntradaRoot))//'_Report.txt', FORM='formatted', POSITION='append')
               l%file11isopen = .true.
 !!!           if (l%layoutnumber==0) call insertalogtmp !ojo lo quito aqui porque borra el _log con la info de credits
               if (l%resume_fromold) then
                  call print11(l%layoutnumber, 'Resuming from .fields.old files')
               ELSE
                  call print11(l%layoutnumber, 'Resuming from .fields files')
               end if
            end if
         ELSE
         !!!if ((l%layoutnumber==0).or.((l%layoutnumber == l%num_procs/2).and.l%stochastic)) then
         !!!   !the temporary
         !!!   CLOSE (11)
         !!!   l%file11isopen=.false.
         !!!   open(11, file=trim(adjustl(l%nEntradaRoot))//'_Report.txt', FORM='formatted')
         !!!   l%file11isopen=.true.
         !!!   if (l%layoutnumber==0) call insertalogtmp
         !!!   call print11 (l%layoutnumber, 'Doing a new simulation from n=1')
         !!!end if
            l%freshstart = .TRUE.
            l%resume_fromold = .FALSE.
         end if
         if (((l%layoutnumber == 0) .or. ((l%layoutnumber == l%num_procs/2) .and. l%stochastic)) .and. l%file11isopen) close (11)
         l%file11isopen = .false.
      end if
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if ((l%run)) then
#ifdef keeppause
          !!!solo para el cluster
         inquire(file='running', EXIST=hayinput)
#ifdef CompileWithMPI
         call MPI_Barrier(SUBCOMM_MPI, l%ierr)
#endif
         if (hayinput) then
            open(9, file='running', FORM='formatted', action='read')
            READ (9, '(a)') chain4
            chain4 = trim(adjustl(chain4))
            CLOSE (9)
                  !!!!!!!!
#ifdef CompileWithMPI
            call MPI_Barrier(SUBCOMM_MPI, l%ierr)
#endif
            call removeintraspaces(l%opcionespararesumeo)
            call removeintraspaces(chain4)
            if (trim(adjustl(l%opcionespararesumeo)) == trim(adjustl(chain4))) then
               existiarunningigual = .true.
            end if
         end if
#endif
#ifdef CompileWithMPI
         call MPI_Barrier(SUBCOMM_MPI, l%ierr)
#endif
         if (l%layoutnumber == 0) then
            open(38, file='running')
            write(38, '(a)') trim(adjustl(l%opcionespararesumeo))
            CLOSE (38)
         end if
      end if
#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI, l%ierr)
#endif
      !open ReportingFiles
      !only the master

      if (((l%layoutnumber == 0) .or. ((l%layoutnumber == l%num_procs/2) .and. l%stochastic)) .and. (statuse /= -1)) then
         print *, 'Opening _Report.txt file'
         if (l%resume) then
            !the temporary
            CLOSE (11)
            l%file11isopen = .false.
            open(11, file=trim(adjustl(l%nEntradaRoot))//'_Report.txt', FORM='formatted')
            l%file11isopen = .true.
            donde = 0
            do while (donde == 0)
               !the first one is a dummy read
               READ (11, '(a)') l%chdummy
               donde = index(l%chdummy, 'mpirun -n')
            end do
            CLOSE (11)
            l%file11isopen = .false.
            l%opcionesoriginales = l%chdummy
!
            call removeintraspaces(l%opcionespararesumeo)
            call removeintraspaces(l%opcionesoriginales)
            if (trim(adjustl(l%opcionesoriginales)) /= trim(adjustl(l%opcionespararesumeo))) then
   call stoponerror(l%layoutnumber, l%num_procs, 'Different resumed/original switches: '//trim(adjustl(l%opcionespararesumeo))//' <> '//&
                        & trim(adjustl(l%opcionesoriginales)), .true.); statuse = -1; !goto 668
            end if
            !
         !!!!!!!!!        CLOSE (11, status='delete')
         !!!!!!!!!        l%file11isopen=.false.
            open(11, file=trim(adjustl(l%nEntradaRoot))//'_Report.txt', FORM='formatted')
            l%file11isopen = .true.
            donde = 0
            do while (donde == 0)
               !the first one is a dummy read
               READ (11, '(a)') l%chdummy
               donde = index(l%chdummy, '!SLICES')
            end do
            l%slicesoriginales = trim(adjustl(l%chdummy))
            CLOSE (11)
            l%file11isopen = .false.
            open(11, file=trim(adjustl(l%nEntradaRoot))//'_Report.txt', FORM='formatted', POSITION='append')
            l%file11isopen = .true.
            if (l%layoutnumber == 0) call insertalogtmp(l)
         ELSE
            CLOSE (11)
            l%file11isopen = .false.
            open(11, file=trim(adjustl(l%nEntradaRoot))//'_Report.txt', FORM='formatted')
            l%file11isopen = .true.
            if (l%layoutnumber == 0) call insertalogtmp(l)
         end if
         !
         call get_secnds(l%time_out2)
         call print_credits(l)
#ifdef CompileWithReal8
         write(dubuf, *) 'Compiled with Double precision (real*8)'
         call print11(l%layoutnumber, dubuf)
#endif
#ifdef CompileWithReal4
         write(dubuf, *) 'Compiled with Single precision (real*4)'
         call print11(l%layoutnumber, dubuf)
#endif
#ifdef CompileWithReal16
         write(dubuf, *) 'Compiled with Quadruple precision (real*16)'
         call print11(l%layoutnumber, dubuf)
#endif
         write(dubuf, *) SEPARADOR//SEPARADOR//SEPARADOR
      !!!call print11 (l%layoutnumber, dubuf,.true.)
         write (11, '(a)') trim(adjustl(dubuf)) !a capon para que el l%stochastic pueda resumear
         write(dubuf, *) 'Launched on              ', l%time_out2%fecha(7:8), '/', l%time_out2%fecha(5:6), '/', &
         &                l%time_out2%fecha(1:4), ' ', l%time_out2%hora(1:2), ':', l%time_out2%hora(3:4)
      !!!call print11 (l%layoutnumber, dubuf,.true.)
         write (11, '(a)') trim(adjustl(dubuf)) !a capon para que el l%stochastic pueda resumear
         write(dubuf, *) SEPARADOR//SEPARADOR//SEPARADOR
      !!!call print11 (l%layoutnumber, dubuf,.true.)
         write (11, '(a)') trim(adjustl(dubuf)) !a capon para que el l%stochastic pueda resumear
         write(dubuf, '(a)') 'Launched with total options '
      !!!call print11 (l%layoutnumber, dubuf,.true.)
         write (11, '(a)') trim(adjustl(dubuf)) !a capon para que el l%stochastic pueda resumear
         write(dubuf, *) trim(adjustl(l%opcionestotales))
      !!!call print11 (l%layoutnumber, dubuf,.true.)
         write (11, '(a)') trim(adjustl(dubuf)) !a capon para que el l%stochastic pueda resumear
         write(dubuf, '(a)') 'If later resuming use compulsory options '
      !!!call print11 (l%layoutnumber, dubuf,.true.)
         write (11, '(a)') trim(adjustl(dubuf)) !a capon para que el l%stochastic pueda resumear
         write(dubuf, *) trim(adjustl(l%opcionespararesumeo))
      !!!call print11 (l%layoutnumber, dubuf,.true.)
         write (11, '(a)') trim(adjustl(dubuf)) !a capon para que el l%stochastic pueda resumear
         write(dubuf, *) SEPARADOR//SEPARADOR//SEPARADOR
         call print11(l%layoutnumber, dubuf)
      end if
      !
      !
      !in seconds
      l%flushsecondsFields = l%flushminutesFields*60
      !in seconds
      l%flushsecondsData = l%flushminutesData*60

      if ((.NOT. l%existeNFDE) .AND. (.NOT. l%existeh5)) then
         call stoponerror(l%layoutnumber, l%num_procs, 'Some input file missing .h5/.nfde/.conf', .true.); statuse = -1; !goto 668
      end if
      !
      !
#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI, l%ierr)
#endif
      if (existiarunningigual) then !lo pongo aqui pq si no no se escribe en el report
         call stoponerror(l%layoutnumber, l%num_procs, 'Running flag file with same options than requested exist. ', .true.); statuse = -1; 
      end if

668   continue

      input_conformal_flag = l%input_conformal_flag  !es un flag global!!!!ojooo 051223 !devolverlo correctamente
      return !el unico return que he dejado !240817

   end subroutine interpreta

   subroutine insertalogtmp(l) !para 100920
      type(entrada_t), intent(INOUT) :: l
      character(len=BUFSIZE) :: dubuf
      integer(kind=4) :: MYUNIT11
      call OffPrint !no reimprimas, esto ya estaba por pantalla
      open(newunit=myunit11, file='SEMBA_FDTD_temp.log')
      do
         read (myunit11, '(1024a)', end=7211) dubuf
         dubuf = '&'//dubuf !para respetar los espacios
         call print11(l%layoutnumber, dubuf)
      end do
7211  CLOSE (myunit11, status='delete')
      call OnPrint
      return
   end subroutine insertalogtmp

   subroutine print_basic_help(l)
      type(entrada_t), intent(INOUT) :: l
      call print_credits(l)
      call print11(l%layoutnumber, '___________________________________________________________________________')
      call print11(l%layoutnumber, 'Basic usage: ')
      call print11(l%layoutnumber, '&   For help use          -h ')
      call print11(l%layoutnumber, '&   For launching use                     ')
      call print11(l%layoutnumber, '&                         -i inputfile (native)')
      call print11(l%layoutnumber, '___________________________________________________________________________')
      return
   end subroutine print_basic_help

   subroutine print_credits(l)
      type(entrada_t), intent(INOUT) :: l
      character(len=BUFSIZE) :: dubuf

      if (l%creditosyaprinteados) return
      l%creditosyaprinteados = .true.
      call print11(l%layoutnumber, '=========================')
      call print11(l%layoutnumber, program_name)
      call print11(l%layoutnumber, '=========================')

      write(dubuf, *) SEPARADOR//SEPARADOR//SEPARADOR
      call print11(l%layoutnumber, dubuf)
      call print11(l%layoutnumber, 'Compilation date: '//compilation_date)
      call print11(l%layoutnumber, 'Compiler Id: '//compiler_id)
      call print11(l%layoutnumber, 'git commit: '//git_commit)
      call print11(l%layoutnumber, 'cmake build type: '//cmake_build_type)
      if (cmake_build_type == "Debug") then
         call print11(l%layoutnumber, 'cmake compilation flags: '//compilation_flags_debug)
      elseif (cmake_build_type == "Release") then
         call print11(l%layoutnumber, 'cmake compilation flags: '//compilation_flags_release)
      else
         call print11(l%layoutnumber, 'cmake compilation flags: '//compilation_flags)
      end if
      call print11(l%layoutnumber, dubuf)
      write(dubuf, *) SEPARADOR//SEPARADOR//SEPARADOR
      call print11(l%layoutnumber, dubuf)
      call print11(l%layoutnumber, 'All rights reserved by the University of Granada (Spain)')
      call print11(l%layoutnumber, '       Contact person: Luis D. Angulo <lmdiazangulo@ugr.es>')
      call print11(l%layoutnumber, ' ')
      !*******************************************************************************

      write(dubuf, *) SEPARADOR//SEPARADOR//SEPARADOR
      call print11(l%layoutnumber, dubuf)
#ifdef CompileWithMPI
      call print11(l%layoutnumber, 'Compiled WITH MPI support')
#endif
#ifdef CompileWithHDF
      call print11(l%layoutnumber, 'Compiled WITH .h5 HDF support')
#endif
#ifdef CompileWithConformal
      call print11(l%layoutnumber, 'Compiled WITH Conformal support')
#endif
#ifdef CompileWithMTLN
      call print11(l%layoutnumber, 'Compiled WITH MTLN support')
#endif
#ifdef CompileWithSMBJSON
      call print11(l%layoutnumber, 'Compiled WITH SMBJSON support')
#endif
      write(dubuf, *) SEPARADOR//SEPARADOR//SEPARADOR
      call print11(l%layoutnumber, dubuf)
      call get_secnds(l%time_out2)
      write(dubuf, *) 'Launched on              ', l%time_out2%fecha(7:8), '/', l%time_out2%fecha(5:6), '/', &
      &                l%time_out2%fecha(1:4), ' ', l%time_out2%hora(1:2), ':', l%time_out2%hora(3:4)
      call print11(l%layoutnumber, dubuf)
      if (l%layoutnumber == 0) print *, 'Highest integer ', huge(1_4)
      return
   end subroutine print_credits

   subroutine print_help(l)
      type(entrada_t), intent(INOUT) :: l
      character(len=BUFSIZE) :: buff
      call print11(l%layoutnumber, '___________________________________________________________________________')
      call print11(l%layoutnumber, 'Command line arguments: ')
      call print11(l%layoutnumber, '___________________________________________________________________________')
      call print11(l%layoutnumber, '-i geometryfile        : Simulates the Native format input file            ')
      call print11(l%layoutnumber, '-r                     : Restarts a previous execution until a given step. ')
      call print11(l%layoutnumber, '&                        Needs -n                                          ')
      call print11(l%layoutnumber, '-run                   : Uses a semaphore running file and automatically   ')
      call print11(l%layoutnumber, '&                        relaunches simulation if ended or aborted (cluter)')
#ifdef CompileWithOldSaving
      call print11(l%layoutnumber, '-old                   : Jointly with -r restarts from .fields.old files   ')
      call print11(l%layoutnumber, '&                        instead (for safety .fields.old fields are saved  ')
      call print11(l%layoutnumber, '&                        too if -flush is issued)                          ')
#endif
      call print11(l%layoutnumber, '-cfl number            : Courant number (suggested<=0.8)  overriding input ')
      call print11(l%layoutnumber, '-n numberoftimesteps   : Run the simulation until a specified step         ')
      call print11(l%layoutnumber, '&                        either restarting if the necessary files are      ')
      call print11(l%layoutnumber, '&                        present, or starting a fresh new one otherwise    ')
      call print11(l%layoutnumber, '&                        Special cases: n=-1 -> Run only .h5/.nfde preproc.')
      call print11(l%layoutnumber, '&                        Special cases: n=-2 -> Run only .h5 preprocessing ')
      call print11(l%layoutnumber, '-s                     : Forces a fresh new simulation, erasing the        ')
      call print11(l%layoutnumber, '&                        restarting files if they are present              ')
      call print11(l%layoutnumber, '&                        Jointly with -n, it enforces a fresh restart      ')
      call print11(l%layoutnumber, '&                        (erases .fields files from previous simulations)  ')
      call print11(l%layoutnumber, '___________________________________________________________________________')
      call print11(l%layoutnumber, '-pause seconds         : Wait seconds to start simulation                  ')
      call print11(l%layoutnumber, '-prefix string         : Adds a string to the output filenames             ')
      call print11(l%layoutnumber, '-saveall               : Saves all the observation time steps              ')
      call print11(l%layoutnumber, '&                        (default saves only the specified windows of time)')
      call print11(l%layoutnumber, '-singlefile            : Compacts E, H, J probes in single files to        ')
      call print11(l%layoutnumber, '&                        overcome a large number of file openings          ')
      !!#ifdef CompileWithMPI
      !!          call print11 (l%layoutnumber, '-maxmessages number    : Buffer of messages for MPI Warnings file. Just    ')
      !!          call print11 (l%layoutnumber, '&                        increase if requested at runtime                  ')
      !!#endif

      !*********************************************************************************************************************
      !***[conformal] **************************************************************************************
      !*********************************************************************************************************************
      !conformal -help printf line   ref:  ##Confhelp##
#ifdef CompileWithConformal
      call print11(l%layoutnumber, '-conf                    : Adds the conformal file to the simulation.')
#endif
      !*********************************************************************************************************************
#ifdef CompileWithNIBC
      call print11(l%layoutnumber, '-skindepthpre          : Pre-processor for sgbc metals including skin depth.')
      call print11(l%layoutnumber, '-mibc                  : Uses pure l%mibc to deal with composites.  ')
      call print11(l%layoutnumber, '-ade                   : Uses l%ade-l%mibc to deal with composites. ')
      call print11(l%layoutnumber, '&                        Alternative to -l%mibc.')
      call print11(l%layoutnumber, '-conformalskin         : Uses a conformal l%mibc to deal with skin-depth')
      call print11(l%layoutnumber, '&                        Do not use this switch if the problem also involves ')
      call print11(l%layoutnumber, '&                        traditional composites, since these do not hold the right ')
      call print11(l%layoutnumber, '&                        thickness parameter. Only use it if the problem only ')
      call print11(l%layoutnumber, '&                        contains metals for which both the conductivity and ')
      call print11(l%layoutnumber, '&                        thickness are CORRECTLY specified in the .nfde file. ')
    call print11(l%layoutnumber, '-NOcompomur            : Uses OLD (possibly unstable) upwinding scheme to deal with composites, ')
 call print11(l%layoutnumber, '&                        instead of the NEW default, which uses a causal time-domain extrapolation ')
      call print11(l%layoutnumber, '&                        of magnetic fields at the surface, by using the one-way ')
      call print11(l%layoutnumber, '&                        advection equation (similar to 1D Mur ABCs) for its ')
      call print11(l%layoutnumber, '&                        superior stability of the default new Mur formulation')
      call print11(l%layoutnumber, '-attc   dissipation    : Positive factor (under 1) for stable composites,   ')
 call print11(l%layoutnumber, '&                        permits to solve some instabilities in the simulation of l%mibc materials.')
      call print11(l%layoutnumber, '&                        It just adds a 1 cell lossy magnetic coating to the l%mibc composite.')
      call print11(l%layoutnumber, '&                        The dissipation factor is used to find the magnetic conductivity ')
      call print11(l%layoutnumber, '&                        from the coefficient updating the current magnetic ')
      call print11(l%layoutnumber, '&                        field from the previous one.  ')
      write (buff, '(a,e10.2e3)') '&                        Default= ', l%attfactorc
      call print11(l%layoutnumber, buff)
#endif
      call print11(l%layoutnumber, '-prioritizeCOMPOoverPEC: Uses Composites instead of PEC in conflicts.       ')
      call print11(l%layoutnumber, '-prioritizeISOTROPICBODYoverall: Uses ISOTROPIC BODY FOR conflicts (JUST FOR SIVA).       ')
      call print11(l%layoutnumber, '-sgbc               : Enables the defaults sgbc model for composites. Default sgbc:')
      call print11(l%layoutnumber, '-nosgbc             : Disables the defaults sgbc model for composites. Default sgbc:')
      call print11(l%layoutnumber, '&                        -sgbfreq 3e9 -sgbresol 1 -sgbcrank      ')
      call print11(l%layoutnumber, '-sgbcfreq           : Maximum frequency to consider the skin-depth       ')
      call print11(l%layoutnumber, '-sgbcresol          : Number of cells per skin-depth a the Maximum frequency')
      call print11(l%layoutnumber, '-sgbcyee            : Uses pure Yee ETD sgbc instead of Crank-Nicolson')
      call print11(l%layoutnumber, '-sgbccrank          : Uses sgbc Crank-Nicolson (default)        ')
      call print11(l%layoutnumber, '-sgbcdepth number   : Overrides automatic calculation of number of cells ')
      call print11(l%layoutnumber, '&                        within sgbc                              ')

      call print11(l%layoutnumber, '-pmlalpha factor order : CPML Alpha factor (>=0, <1 sug.) & polyn. grading.')
      call print11(l%layoutnumber, '&                        alpha=factor * maximum_PML_sigma , order=polynom. ')
      write (buff, '(a,2e10.2e3)') '&                        Default= ', l%alphamaxpar, l%alphaOrden
      call print11(l%layoutnumber, buff)
      write (buff, '(a,e10.2e3)') '-pmlkappa number       : CPML Kappa (>=1). Default= ', l%kappamaxpar
      call print11(l%layoutnumber, buff)
      call print11(l%layoutnumber, '-pmlcorr factor depth  : Factor for CPML enhanced stability (default none).')
      call print11(l%layoutnumber, '&                        sigma=factor * maximum_PML_sigma, depth= # layers ')
      call print11(l%layoutnumber, '-mur1                  : Supplement PMLs with 1st order Mur ABCs           ')
      call print11(l%layoutnumber, '-mur2                  : Supplement PMLs with 2nd order Mur ABCs           ')
      call print11(l%layoutnumber, '-wiresflavor {holland.or.old} : model for the wires    ')
#ifdef CompileWithBerengerWires
      call print11(l%layoutnumber, '-wiresflavor {berenger} : model for the wires    ')
#endif
#ifdef CompileWithSlantedWires
      call print11 (l%layoutnumber, '-wiresflavor {new/Slanted.or.experimental.or.slanted/transition/semistructured l%precision} : model for the wires    ')   
#endif
      call print11(l%layoutnumber, '&                        (default '//trim(adjustl(l%wiresflavor))//')   ')
      call print11(l%layoutnumber, '-notaparrabos          : Do not remove extra double tails at the end of the wires ')
      call print11(l%layoutnumber, '&                        only available for the native format.             ')
      call print11(l%layoutnumber, '-intrawiresimplify     : Disable strict interpretation of .NFDE topology.  ')
      call print11(l%layoutnumber, '&                        Collapse internal parallel wires and create       ')
      call print11(l%layoutnumber, '&                        intra-wire junctions.                             ')
      call print11(l%layoutnumber, '-nomtlnberenger        : Disables MTLN improvements for Berenger l%wiresflavor')
      call print11(l%layoutnumber, '-stableradholland             : Automatic correction of radii for Holland l%wiresflavor')
      call print11(l%layoutnumber, '&                        Use only in case of instabilities.  (experimental)')
      call print11(l%layoutnumber, '-groundwires           : Ground wires touching/embedded/crossing PEC/Lossy.')
      call print11(l%layoutnumber, '&                        Use with CAUTION. Revise *Warnings.txt file!      ')
      call print11(l%layoutnumber, '-noSlantedcrecepelo : Ground open nodes. Experimental. Do not use.')
      call print11(l%layoutnumber, '-connectendings        : Joins ohmicly endings nodes of adjacent segments  ')
      call print11(l%layoutnumber, '&                        from multiwires (segments do no collapse).        ')
      call print11(l%layoutnumber, '&                        regardless of whether they are actually connected ')
      call print11(l%layoutnumber, '&                        through the LeftEnd/RightEnd numbering ')
      call print11(l%layoutnumber, '&                        Automatic with -a                                 ')
      call print11(l%layoutnumber, '&                        Use with CAUTION. Revise *Warnings.txt file!      ')
      call print11(l%layoutnumber, '-isolategroupgroups    : Detach ohmicly endings nodes of adjacent segments ')
      call print11(l%layoutnumber, '&                        from multiwires if they are in different          ')
      call print11(l%layoutnumber, '-makeholes             : Create a void 2-cell area around wire segments    ')
      call print11(l%layoutnumber, '&                        Use with CAUTION. Revise *Warnings.txt (experim.) ')
      call print11(l%layoutnumber, '-mindistwires dist     : Specify the min distance between wires in a       ')
      call print11(l%layoutnumber, '&                        multiwire in new and experimental wires flavors   ')
      write (buff, '(a,e10.2e3)') '&                        Default= ', l%mindistwires
      call print11(l%layoutnumber, buff)
      call print11(l%layoutnumber, '-inductance {ledfelt/berenger/boutayeb} : model for the self-inductance    ')
      call print11(l%layoutnumber, '&                        (default '//trim(adjustl(l%inductance_model))//')   ')
      call print11(l%layoutnumber, '-inductanceorder order : order for the self-inductance calculation for     ')
      call print11(l%layoutnumber, '&                        slanted wires in experimental l%wiresflavor         ')
      write (buff, '(a,i8)') '&                        Default= ', l%inductance_order
      call print11(l%layoutnumber, '-attw   dissipation    : Positive factor (under 1) for stability in wires, ')
      write (buff, '(a,e10.2e3)') '&                        Default= ', l%attfactorw
      call print11(l%layoutnumber, '-maxwireradius number  : Bounds globally the wire radius                   ')
      call print11(l%layoutnumber, '-clip                  : Permits to clip a bigger problem truncating wires.')
      call print11(l%layoutnumber, '-wirecrank             : Uses Crank-Nicolson for wires (development)       ')
      call print11(l%layoutnumber, '-noNF2FF string        : Supress a NF2FF plane for calculation             ')
      call print11(l%layoutnumber, '&                        String can be: up, down, left, right, back , front')
      call print11(l%layoutnumber, '-NF2FFDecim            : Uses decimation in NF2FF calculation (faster).    ')
      call print11(l%layoutnumber, '&                        WARNING: High-freq aliasing may occur             ')
      call print11(l%layoutnumber, '-vtkindex              : Output index instead of real point in 3D slices.  ')
      call print11(l%layoutnumber, '-ignoreerrors          : Run even if errors reported in *Warnings.txt file.')
      call print11(l%layoutnumber, '___________________________________________________________________________')
      call print11(l%layoutnumber, '-cpumax minutes        : CPU runtime (useful for limited CPU queuing       ')
      call print11(l%layoutnumber, '-noshared              : Do not waste time with shared fields              ')
      call print11(l%layoutnumber, '-flush minutes         : Minutes between data flush of restarting fields   ')
      call print11(l%layoutnumber, '&                        (default 0=No flush)                              ')
      call print11(l%layoutnumber, '-flushdata minutes     : Minutes between flushing observation data         ')
      call print11(l%layoutnumber, '&                        (default is every 5 minutes)                      ')
      call print11(l%layoutnumber, '-map                   : Creates map ASCII files of the geometry           ')
      call print11(l%layoutnumber, '&                        with wires and PEC                ')
      call print11(l%layoutnumber, '&                        (in conjunction with -n 0 only creates the maps)  ')
      call print11(l%layoutnumber, '-mapvtk                : Creates .VTK map of the PEC/wires/Surface geometry')
#ifdef CompileWithConformal
      call print11(l%layoutnumber, '-conf file             : conformal file  ')
      call print11(l%layoutnumber, '-abrezanjas            : Thin-gaps treated in conformal manner  ')
#endif
      call print11(l%layoutnumber, '-dmma                  : Thin-gaps treated in DMMA manner  ')
#ifdef CompileWithMPI
      call print11(l%layoutnumber, '-mpidir {x,y,z}        : Rotate model to force MPI along z be the largest  ')
      call print11(l%layoutnumber, '-force    cutplane     : Force a MPI layout to begin at cutplane (debug!)  ')
#endif
      call print11(l%layoutnumber, '___________________________________________________________________________')
      call print11(l%layoutnumber, 'Control through signaling files during the simulation: (after erased)      ')
      call print11(l%layoutnumber, '&  stop         : (void) Forces a graceful end (it Cannot be resumed)      ')
      call print11(l%layoutnumber, '&                 No restarting data is flushed, only observation data     ')
      call print11(l%layoutnumber, '&  stopflushing : (void) Forces a graceful end (it can be resumed)         ')
      call print11(l%layoutnumber, '&  flush        : (void) Forces a flush of resuming fields and observation ')
      call print11(l%layoutnumber, '&                 data in 1 minute time approx.                            ')
      call print11(l%layoutnumber, '&  flushdata    : (void) Forces a flush only of the observation data in    ')
      call print11(l%layoutnumber, '&                 1 minute time approx.                                    ')
      call print11(l%layoutnumber, '&                 Both restarting and observation data are flushed         ')
      call print11(l%layoutnumber, '&  stop_only         : Forces a graceful end (cannot be resumed) only of a ')
      call print11(l%layoutnumber, '&                      given problem name (without the .nfde extension)    ')
      call print11(l%layoutnumber, '&                      No restarting data is flushed, only observation data')
      call print11(l%layoutnumber, '&  stopflushing_only : Forces a graceful end (it can be resumed) only of a ')
      call print11(l%layoutnumber, '&                      give problem name (without the .nfde extension)     ')
      call print11(l%layoutnumber, '&                      Both restarting and observation data is flushed     ')
      call print11(l%layoutnumber, '&  flush_only   : Forces flush of resuming fields and observation data only')
      call print11(l%layoutnumber, '&                 of a given problem name (without the .nfde extension)    ')
      call print11(l%layoutnumber, '&                 in 1 minute time approx.                                 ')
      call print11(l%layoutnumber, '&  flushdata_only : Forces a flush only of the observation data only of a  ')
      call print11(l%layoutnumber, '&                   given problem name (without the .nfde extension)       ')
      call print11(l%layoutnumber, '&                   in 1 minute time approx.                               ')
      call print11(l%layoutnumber, '&                   Both restarting and observation data are flushed       ')
      call print11(l%layoutnumber, '&  pause        : (void) While this field exist no simulation is started   ')
      call print11(l%layoutnumber, '&  unpack       : (void) Unpacks on-the-fly .bin probes files created      ')
      call print11(l%layoutnumber, '&                 with the -singlefile packaging option                    ')
      call print11(l%layoutnumber, '&  postprocess  : (void) Do frequency domain and transfer function         ')
      call print11(l%layoutnumber, '&                 postprocess on-the-fly                                   ')
      call print11(l%layoutnumber, '&  flushxdmf    : (void) Flush .xdmf animation probes on the fly           ')
      call print11(l%layoutnumber, '&  flushvtk     : (void) Flush .vtk  animation probes on the fly           ')
      call print11(l%layoutnumber, '&  snap         : Creates a .h5 and .xdmf snapshot per MPI layout if the   ')
      call print11(l%layoutnumber, '&                 field value is over the first number found in this file  ')
      call print11(l%layoutnumber, '&                 in space steps by the 2nd integer number                 ')
      call print11(l%layoutnumber, '&                 in time steps by the 3rd integer number (1-minute lapse) ')
      call print11(l%layoutnumber, '&  relaunch     : Relaunches the simulation upon termination with the      ')
      call print11(l%layoutnumber, '&                 switches read from this file. Used jointly with a        ')
      call print11(l%layoutnumber, '&                 stop file permits to launch simulations on-demand        ')
      call print11(l%layoutnumber, '___________________________________________________________________________')
      !
      write (buff, '(a,i14,a)') 'Max CPU time is ', topCPUtime, ' seconds (can be overriden by -cpumax)'
      call print11(l%layoutnumber, buff)
#ifdef CompileWithOpenMP
      call print11(l%layoutnumber, 'SUPPORTED:   MultiCPU parallel simulation (OpenMP)')
#else
      !call print11 (l%layoutnumber, 'UNSUPPORTED: MultiCPU parallel simulation (OpenMP)')
#endif
!
#ifdef CompileWithMPI
      call print11(l%layoutnumber, 'SUPPORTED:   MultiCPU/Multinode parallel simulation (MPI)')
#else
      !call print11 (l%layoutnumber, 'UNSUPPORTED: MultiCPU/Multinode parallel simulation (MPI)')
#endif
#ifdef CompileWithConformal
      call print11(l%layoutnumber, 'SUPPORTED:   Conformal algorithm')
#else
      !call print11 (l%layoutnumber, 'UNSUPPORTED: Conformal algorithm')
#endif
      call print11(l%layoutnumber, 'SUPPORTED:   Near-to-Far field probes')
      call print11(l%layoutnumber, 'SUPPORTED:   Lossy anistropic materials, both electric and magnetic')
      call print11(l%layoutnumber, 'SUPPORTED:   Thin Slots ')
      call print11(l%layoutnumber, 'SUPPORTED:   Electric and Magnetic Dispersive materials ')
      call print11(l%layoutnumber, 'SUPPORTED:   Isotropic Multilayer Skin-depth Materials (sgbc)')
#ifdef CompileWithNIBC
      call print11(l%layoutnumber, 'SUPPORTED:   Isotropic Multilayer Skin-depth Materials (l%mibc)')
#else
      !call print11 (l%layoutnumber, 'UNSUPPORTED: Isotropic Multilayer Skin-depth Materials (l%mibc)')
#endif
      call print11(l%layoutnumber, 'SUPPORTED:   Loaded and grounded thin-wires with juntions')
      call print11(l%layoutnumber, 'SUPPORTED:   Nodal hard/soft electric and magnetic sources')
#ifdef CompileWithHDF
      call print11(l%layoutnumber, 'SUPPORTED:   .xdmf+.h5 probes ')
#else
      !call print11 (l%layoutnumber, 'UNSUPPORTED: .xdmf+.h5 probes ')
#endif
#ifdef CompileWithOldSaving
      call print11(l%layoutnumber, 'SUPPORTED:   .fields.old files created (fail-safe)')
#else
      !call print11 (l%layoutnumber, 'UNSUPPORTED: .fields.old files created (fail-safe)')
#endif
#ifdef CompileWithStochastic
      call print11(l%layoutnumber, 'SUPPORTED:   l%stochastic analysis')
#else
      !call print11 (l%layoutnumber, 'UNSUPPORTED: l%stochastic analysis')
#endif
#ifdef CompileWithPrescale
      call print11(l%layoutnumber, 'SUPPORTED:   Permittivity scaling accelerations')
#else
      !call print11 (l%layoutnumber, 'UNSUPPORTED: Permittivity scaling accelerations')
#endif
      call print11(l%layoutnumber, 'SUPPORTED:   Holland Wires')
#ifdef CompileWithBerengerWires
      call print11(l%layoutnumber, 'SUPPORTED:   Multi-Wires')
#else
      !call print11 (l%layoutnumber, 'UNSUPPORTED: Multi-Wires')
#endif
#ifdef CompileWithSlantedWires
      call print11(l%layoutnumber, 'SUPPORTED:   Slanted Wires')
#else
      !call print11 (l%layoutnumber, 'UNSUPPORTED: Slanted Wires')
#endif
!!!!!!!!!!!!!!!!!
#ifdef CompileWithReal4
      call print11(l%layoutnumber, 'Single precission simulations (reals are 4-byte)')
#endif
#ifdef CompileWithReal8
      call print11(l%layoutnumber, 'Double precission simulations (reals are 8-byte)')
#endif
#ifdef CompileWithInt4
      call print11(l%layoutnumber, 'Media matrices are 4 bytes')
#endif
#ifdef CompileWithInt2
      call print11(l%layoutnumber, 'Media matrices are 2 bytes')
#endif
#ifdef CompileWithInt1
      call print11(l%layoutnumber, 'Media matrices are 1 byte')
#endif
#ifdef CompileWithMPI
      call MPI_FINALIZE(l%ierr)
#endif
      return
   end subroutine print_help

   subroutine removeintraspaces(a)
      character(len=*), intent(inout):: a
      integer(Kind=4) :: i, longi
      logical correc
      correc = .true.
      do while (correc)
         correc = .false.
         a = trim(adjustl(a))
         longi = len(trim(adjustl(a)))
         buscae: do i = 1, longi - 1
            if ((a(i:i) == ' ') .and. (a(i + 1:i + 1) == ' ')) then
               a = a(1:i - 1)//a(i + 1:longi)
               correc = .true.
               exit buscae
            end if
         end do buscae
      end do
      return
   end subroutine removeintraspaces

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine buscaswitchficheroinput(l)

!!!!!!!!!!!!!
      type(entrada_t), intent(INOUT) :: l
!!!!!!!!!

      character(len=BUFSIZE) :: dato, buff, f, binaryPath
      integer(kind=4) :: i, n, statuse, NUM_NFDES, TEMP_NUMNFDES, p
      character(len=5) :: NFDEEXTENSION, CONFEXTENSION, CMSHEXTENSION

!!!

      NFDEEXTENSION = '.nfde'; CONFEXTENSION = '.conf'; CMSHEXTENSION = '.cmsh'
      statuse = 0
   !!!!!!!!!!!!!!!
      binaryPath = getBinaryPath()
      n = commandargumentcount(l%chain2, binaryPath)
      if (n == 0) then
         call print_basic_help(l)
      call stoponerror(l%layoutnumber,l%num_procs,'Error: NO arguments neither command line nor in launch file. Correct and remove pause...',.true.)
         statuse = -1
         goto 667
      end if

      if (n > 0) then
         num_nfdes = 0
         i = 2
         do while (i <= n)
            call getcommandargument(l%chain2, i, l%chain, l%length, statuse, binaryPath)
            if (statuse /= 0) then
               call stoponerror(l%layoutnumber, l%num_procs, 'Reading input', .true.)
               goto 667
            end if
            !
            SELECT CASE (trim(adjustl(l%chain)))
            CASE ('-mpidir')
               i = i + 1
               call getcommandargument(l%chain2, i, f, l%length, statuse, binaryPath)
               select case (trim(adjustl(f)))
               case ('x', 'X')
                  l%mpidir = 1  
               case ('y', 'Y')
                  l%mpidir = 2   
               case ('z', 'Z')
                  l%mpidir = 3
               CASE DEFAULT
                  GOTO 1762
               end select
               GO TO 2762
1762           call stoponerror(l%layoutnumber, l%num_procs, 'Invalid -l%mpidir option', .true.)
               statuse = -1
               goto 667
2762           continue
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
               num_nfdes = num_nfdes + 1
            end select
            i = i + 1
         end do
         if (num_nfdes > 1) then
            temp_numnfdes = 0
            i = 2 ! se empieza en 2 porque el primer argumento es siempre el nombre del ejecutable
            do while (i <= n)
               call getcommandargument(l%chain2, i, l%chain, l%length, statuse, binaryPath)
               if (statuse /= 0) then
                  call stoponerror(l%layoutnumber, l%num_procs, 'Reading input', .true.)
                  goto 667
               end if
               !
               SELECT CASE (trim(adjustl(l%chain)))
               CASE ('-i')
                  temp_numnfdes = temp_numnfdes + 1
                  i = i + 1
                  call getcommandargument(l%chain2, i, f, l%length, statuse, binaryPath)
                  p = LEN_trim(adjustl(f))
                  if ((p - 4) >= 1) then
                     if (f((p - 4):(p - 4)) == NFDEEXTENSION(1:1)) then
                        NFDEEXTENSION = f((p - 4):p)
                        l%extension = NFDEEXTENSION
                        l%fichin = f(1:p - 5)
                     ELSE
                        l%fichin = f(1:p)
                     end if
                  ELSE if (p >= 1) then
                     l%fichin = f(1:p)
                  ELSE
                     call stoponerror(l%layoutnumber, l%num_procs, 'There is not a .nfde file for input', .true.)
                     statuse = -1
                     goto 667
                  end if
                  inquire(file=trim(adjustl(l%fichin))//NFDEEXTENSION, EXIST=l%existeNFDE)
                  if (.NOT. l%existeNFDE) then
                     buff = 'The input file was not found '//trim(adjustl(l%fichin))//NFDEEXTENSION
                     call stoponerror(l%layoutnumber, l%num_procs, buff, .true.)
                     statuse = -1
                     goto 667
                  end if
!aniadido para chequear que no haya .conf sin haber invocado el -conf 15/12/16 sgg
                  inquire(file=trim(adjustl(l%fichin))//CONFEXTENSION, EXIST=l%existeconf)
                  if ((l%existeconf) .AND. (.not. (l%input_conformal_flag))) then
 buff = 'No -conf issued but existing file '//trim(adjustl(l%fichin))//confEXTENSION//' . Either remove file or relaunch with -conf'
                     call stoponerror(l%layoutnumber, l%num_procs, buff, .true.)
                     statuse = -1
                     goto 667
                  end if
                  inquire(file=trim(adjustl(l%fichin))//CMSHEXTENSION, EXIST=l%existecmsh)
                  if ((l%existecmsh) .AND. (.not. (l%input_conformal_flag))) then
 buff = 'No -conf issued but existing file '//trim(adjustl(l%fichin))//CMSHEXTENSION//' . Either remove file or relaunch with -conf'
                     call stoponerror(l%layoutnumber, l%num_procs, buff, .true.)
                     statuse = -1
                     goto 667
                  end if
!
                  if (temp_numnfdes == 1) then !solo el primero
                     if (l%layoutnumber == 0) open (194, file='multi_'//trim(adjustl(l%fichin))//NFDEEXTENSION, form='formatted')
                  end if
                  if (l%layoutnumber == 0) then
                     open (196, file=trim(adjustl(l%fichin))//NFDEEXTENSION, form='formatted')
                     do
                        read (196, '(a)', end=197) dato
                        if (trim(adjustl(dato)) /= '!END') then
                           write (194, '(a)') trim(adjustl(dato))
                        else
                           dato = '***** End merging file: '//trim(adjustl(l%fichin))//NFDEEXTENSION//' ********'
                           write (194, '(a)') trim(adjustl(dato))
                        end if
                     end do
197                  close (196)
                  end if
                  if (temp_numnfdes == num_nfdes) then !solo el primero
                     if (l%layoutnumber == 0) then
                        write (194, '(a)') '!END'
                        close (194)
                     end if
                  end if
               end select
               i = i + 1
            end do
         end if
      end if
#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI, l%ierr)
#endif
      !
      !concatenado multiples ORIGINAL 26/06/14

      !fin concatenado

      temp_numnfdes = 0
      if (n > 0) then
         i = 2  ! se empieza en 2 porque el primer argumento es siempre el nombre del ejecutable
         do while (i <= n)
            call getcommandargument(l%chain2, i, l%chain, l%length, statuse, binaryPath)
            if (statuse /= 0) then
               call stoponerror(l%layoutnumber, l%num_procs, 'Reading input', .true.)
               goto 667
            end if
            !
            SELECT CASE (trim(adjustl(l%chain)))
               !
            CASE ('-i')
               temp_numnfdes = temp_numnfdes + 1
               i = i + 1
               if (temp_numnfdes == 1) then
                  !
                  call getcommandargument(l%chain2, i, f, l%length, statuse, binaryPath)
                  p = LEN_trim(adjustl(f))
                  if ((p - 4) >= 1) then
                     if (f((p - 4):(p - 4)) == NFDEEXTENSION(1:1)) then
                        NFDEEXTENSION = f((p - 4):p)
                        l%extension = NFDEEXTENSION
                        l%fichin = f(1:p - 5)
                     ELSE
                        l%fichin = f(1:p)
                     end if
                  ELSE if (p >= 1) then
                     l%fichin = f(1:p)
                  ELSE
                     call stoponerror(l%layoutnumber, l%num_procs, 'There is not a .nfde file for input', .true.)
                     statuse = -1
                     goto 667
                  end if
                  block
                     character(len=255) :: cwd
                     call getcwd(cwd)
                     write(*, *) TRIM(cwd)
                  end block
                  inquire(file=trim(adjustl(l%fichin))//NFDEEXTENSION, EXIST=l%existeNFDE)
                  if (.NOT. l%existeNFDE) then
                     buff = 'The input file was not found '//trim(adjustl(l%fichin))//NFDEEXTENSION
                     call stoponerror(l%layoutnumber, l%num_procs, buff, .true.)
                     statuse = -1
                     goto 667
                  end if
               elseif (temp_numnfdes == 2) then
                  l%fichin = 'multi_'//trim(adjustl(l%fichin))
               else
                  l%fichin = l%fichin
               end if
            !!!          l%opcionespararesumeo = trim (adjustl(l%opcionespararesumeo)) // ' ' // trim (adjustl(l%chain)) // ' ' // trim (adjustl(f))
            end select
            i = i + 1
         end do
      end if
      !
      ! If no input is present we stop
      if (len(trim(adjustl(l%fichin))) <= 0) then
         call stoponerror(l%layoutnumber, l%num_procs, 'ERROR! -> No input file was specified. Use -i ****.fdtd.json', .true.); statuse = -1; goto 667
      end if

      l%fileFDE = trim(adjustl(l%fichin))//NFDEEXTENSION
      l%fileH5 = trim(adjustl(l%fichin))//'.h5'
      call INITWARNINGFILE(l%layoutnumber, l%num_procs, trim(adjustl(l%fichin))//'_tmpWarnings.txt', l%verbose, l%ignoreerrors)

!!!
667   return
   end subroutine buscaswitchficheroinput

!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine default_flags(l)
!!!!!!!!!!!!!
      type(entrada_t), intent(INOUT) :: l
      l%noconformalmapvtk = .false.
      l%forced = -1
      l%sgbcdepth = -1
      l%statuse = 0
      l%time_begin = 0
!
      l%precision = 0 !redondeo del semiestructurado
      l%stochastic = .false.
      l%chosenyesornostochastic = .false. !es un flag informativo que debe inicializarse a .false. a pesar de qu el sentido comun diga lo contrario
      l%simu_devia = .false.

#ifdef CompileWithHDF
      l%createh5bin = .false.
#else
      l%createh5bin = .true.
#endif
      l%createh5filefromsinglebin = .false.
      l%permitscaling = .false.
      l%niapapostprocess = .false.
      l%prioritizeCOMPOoverPEC = .false.  !pec has default more priority than compo (para siva hay que cambiarlo)
      l%prioritizeTHINWIRE = .false. !solo para visualizacion y experimentacion 231024
      l%prioritizeISOTROPICBODYoverall = .FALSE. !PARA EL SIVA SE CAMBIA POR LINEA DE COMANDO
      l%mpidir = 3 !DEFAULT do NOT ROTATE GEOMETRY !JUST TO TAKE PROFIT OF MPI
      l%maxwireradius = -1.0_RKIND
      l%boundwireradius = .false.
      l%wirecrank = .FALSE.
      l%ignoreerrors = .false.
      l%ignoresamplingerrors = .false.
      l%vtkindex = .FALSE. !SOLO AFECTA A LOS VTK (SACA INDICES EN VEZ DE POSICION FISICA)
      l%CLIPREGION = .false.
      l%NF2FFDecim = .FALSE.
      l%facesNF2FF%tr = .true.
      l%facesNF2FF%fr = .true.
      l%facesNF2FF%iz = .true.
      l%facesNF2FF%de = .true.
      l%facesNF2FF%ab = .true.
      l%facesNF2FF%ar = .true.
      !defaults
      l%read_command_line = .true.
      l%hay_slanted_wires = .false.
      l%forcing = .FALSE.
      l%resume_fromold = .FALSE.
      l%singlefilewrite = .FALSE.
      l%updateshared = .true.

      l%finaltimestep = 0
      l%cfltemp = 1.0 !dummy
      l%cfl = 1.0 !default courant number !no tocarlo 310715 solo afecta si se usa -l%cfl
      l%forcecfl = .false.
      !PML default
      !cpml stretching maximum parameters !!l%alphamaxpar=StaticFrequency*2*pi*Eps0
      l%alphamaxpar = 0.0_RKIND !0.24  !expresion 7.78 taflove 3 edic)
      l%alphaOrden = 1.0_RKIND
      l%kappamaxpar = 1.0_RKIND !15.0_RKIND !061118 mantener a 1 por conflictos cpml and permittivity scaling
      !and final layer electric sigma
      l%MEDIOEXTRA%exists = .false.
      l%MEDIOEXTRA%index = -7 !void
      l%MEDIOEXTRA%pml_size = -1  !void
      l%MEDIOEXTRA%sigma = -1e20 !void
      !
      l%MurAfterPML = .false.
      l%mur_second = .false.
      l%mur_first = .false.
      l%mur_exist = .false.
      !!!!!!!!!!!!
      l%takeintcripte = .false. !a peticion de OLD, redondear los nodos de cripte a la baja
      l%attfactorc = 1.0_RKIND !default dissipation factor for composites
      l%attfactorw = 1.0_RKIND!default dissipation factor for wires

      l%mibc = .false.
      l%ade = .false. !auxiliary differential equation for composites
      l%conformalskin = .false.
      l%sgbc = .true. !default is false unless required
      l%sgbcDispersive = .false. !default is false unless required
      l%skindepthpre = .false.
      l%sgbcdepth = -1  ! se calcula automaticamente a menos que se use el switch
      l%sgbcfreq = 1e9 !default es cazar el skin depth hasta 1e9
      l%sgbcresol = 1.0 !numero de celdas por skin depth (caida a exp(-1))
      l%sgbccrank = .true. !default es l%sgbccrank

      l%fatalerror = .false.
      l%fatalerrornfde2sgg = .false.
      !**************************************************************************************************
      !***[conformal] *******************************************************************
      !**************************************************************************************************
      !conformal existence flags   ref: ##Confflag##
      l%input_conformal_flag = .false.
      l%flag_conf_sgg = .false.
      !

      l%dontwritevtk = .false.

      l%NOcompomur = .false. !DEFAULT mi formulacion

      !default no join the wires which are adjacent (ORIGINAL election)
      !do not connect endings unless specified in ORIGINAL
      l%makeholes = .FALSE.
      l%connectendings = .false.
      l%isolategroupgroups = .false.
      l%strictOLD = .true. !default is strict ORIGINAL overriden manually
      l%TAPARRABOS = .true. !default since 101116 !cortar los multizigzag rabitos
      l%mtlnberenger = .true. !solo actua si se invoca con l%wiresflavor berenger esto va a ser siempre true a menos que tambien se invoque con -nomtlnberenger (solo para debugeo y que coincida con Holland) 020719
      l%stableradholland = .false. !solo actua si se invoca con l%wiresflavor holland
      l%fieldtotl = .false.
      l%experimentalVideal = .false.
      l%thereare_stoch = .false.
      l%forceresampled = .false.
      l%factorradius = 1.0e+30 !para evitar division por cero 120123
      l%factordelta = 1.0e+30 !para evitar division por cero 120123
      !default
      l%groundwires = .false.
      l%noSlantedcrecepelo = .false. !131219 experimental niapa ojoooo
      l%inductance_model = 'boutayeb'
      l%inductance_order = 8
      l%wiresflavor = 'holland'
      l%wirethickness = 1
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
      l%existeh5 = .FALSE.
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
      l%file11isopen = .false.
      l%relaunching = .false.
      l%forcestop = .false.
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

      input_conformal_flag = l%input_conformal_flag    !ojooo 051223 es un flag globaaaallll
      return
   end subroutine default_flags

end module interpreta_switches_m

