integer function test_read_sphere() bind (C) result(err)
   use smbjson
   use smbjson_testingTools

   implicit none

   character(len=*),parameter :: filename = PATH_TO_TEST_DATA//INPUT_EXAMPLES//'sphere.fdtd.json'
   type(Parseador) :: pr, ex
   type(parser_t) :: parser
   logical :: areSame
   err = 0

   ex = expectedProblemDescription()
   parser = parser_t(filename)
   pr = parser%readProblemDescription()
   
   call expect_eq(err, ex, pr, ignoreRegions=.true.)
   
   if (err == 0) write(*,*) "Read and expected inputs are equal."   
   
contains
   function expectedProblemDescription() result (ex)
      type(Parseador) :: ex

      call initializeProblemDescription(ex)

      ! Expected general info.
      ex%general%dt = 3.85167e-11
      ex%general%nmax = 100

      ! Excected media matrix.
      ex%matriz%totalX = 80
      ex%matriz%totalY = 80
      ex%matriz%totalZ = 80

      ! Expected grid.
      ex%despl%nX = 80
      ex%despl%nY = 80
      ex%despl%nZ = 80

      allocate(ex%despl%desX(80))
      allocate(ex%despl%desY(80))
      allocate(ex%despl%desZ(80))
      ex%despl%desX = 0.025
      ex%despl%desY = 0.025
      ex%despl%desZ = 0.025
      ex%despl%mx1 = 0
      ex%despl%mx2 = 80
      ex%despl%my1 = 0
      ex%despl%my2 = 80
      ex%despl%mz1 = 0
      ex%despl%mz2 = 80

      ! Expected boundaries.
      ex%front%tipoFrontera(:) = F_PML
      ex%front%propiedadesPML(:)%numCapas = 10
      ex%front%propiedadesPML(:)%orden = 2
      ex%front%propiedadesPML(:)%refl = 0.001

      ! Expected material regions.
      ex%pecRegs%nSurfs = 1
      ex%pecRegs%nSurfs_max = 1
      allocate(ex%pecRegs%Surfs(1))
      ! -- specific surfs not included DO NOT use comparison --

      ! Expected sources.
      allocate(ex%plnSrc%collection(1))
      ex%plnSrc%collection(1)%nombre_fichero = "gauss.exc"
      ex%plnSrc%collection(1)%atributo = ""
      ex%plnSrc%collection(1)%coor1 = [0, 0, 0]
      ex%plnSrc%collection(1)%coor2 = [79, 79, 79]
      ex%plnSrc%collection(1)%theta = 1.5707963268
      ex%plnSrc%collection(1)%phi = 0.0
      ex%plnSrc%collection(1)%alpha = 1.5707963268
      ex%plnSrc%collection(1)%beta = 4.7123889804
      ex%plnSrc%collection(1)%isRC=.false.
      ex%plnSrc%collection(1)%nummodes=1
      ex%plnSrc%collection(1)%INCERTMAX=0.0
      ex%plnSrc%nc = 1
      ex%plnSrc%nC_max = 1

      ! Expected probes
      allocate(ex%oldSONDA)
      ex%oldSONDA%n_probes = 1
      ex%oldSONDA%n_probes_max = 1
      allocate(ex%oldSONDA%probes(1))
      ex%oldSonda%probes(1)%n_FarField = 1
      ex%oldSonda%probes(1)%n_FarField_max = 1
      allocate(ex%oldSonda%probes(1)%FarField(1))
      ex%oldSonda%probes(1)%FarField(1)%probe%grname = " "
      ex%oldSonda%probes(1)%FarField(1)%probe%outputrequest = "FarField_log_"
      ex%oldSonda%probes(1)%FarField(1)%probe%tstart = 0.0
      ex%oldSonda%probes(1)%FarField(1)%probe%tstop = 0.0
      ex%oldSonda%probes(1)%FarField(1)%probe%tstep = 0.0
      ex%oldSonda%probes(1)%FarField(1)%probe%fstart = 1e6
      ex%oldSonda%probes(1)%FarField(1)%probe%fstop = 1e9
      ex%oldSonda%probes(1)%FarField(1)%probe%fstep = (1e9-1e6)/5
      ex%oldSONDA%probes(1)%FarField(1)%probe%FileNormalize = "gauss.exc"
      allocate(ex%oldSonda%probes(1)%FarField(1)%probe%i(2))
      allocate(ex%oldSonda%probes(1)%FarField(1)%probe%j(2))
      allocate(ex%oldSonda%probes(1)%FarField(1)%probe%k(2))
      allocate(ex%oldSonda%probes(1)%FarField(1)%probe%node(0))
      ex%oldSonda%probes(1)%FarField(1)%probe%i = [2, 77]
      ex%oldSonda%probes(1)%FarField(1)%probe%j = [2, 77]
      ex%oldSonda%probes(1)%FarField(1)%probe%k = [2, 77]
      ex%oldSonda%probes(1)%FarField(1)%probe%n_cord = 2
      ex%oldSonda%probes(1)%FarField(1)%probe%n_cord_max = 2
      ex%oldSonda%probes(1)%FarField(1)%probe%thetastart = 0.0
      ex%oldSonda%probes(1)%FarField(1)%probe%thetastop  = 180.0
      ex%oldSonda%probes(1)%FarField(1)%probe%thetastep  = 90.0 
      ex%oldSonda%probes(1)%FarField(1)%probe%phistart   = 0.0
      ex%oldSonda%probes(1)%FarField(1)%probe%phistop    = 360.0
      ex%oldSonda%probes(1)%FarField(1)%probe%phistep    = 90.0   
      
      allocate(ex%VolPrb)
      ex%VolPrb%length = 1
      ex%VolPrb%length_max = 1
      ex%VolPrb%len_cor_max = 2
      allocate(ex%VolPrb%collection(1))
      allocate(ex%VolPrb%collection(1)%cordinates(1))
      ex%VolPrb%collection(1)%len_cor = 1
      ex%VolPrb%collection(1)%cordinates(1)%Xi = 2
      ex%VolPrb%collection(1)%cordinates(1)%Xe = 77
      ex%VolPrb%collection(1)%cordinates(1)%Yi = 2
      ex%VolPrb%collection(1)%cordinates(1)%Ye = 77
      ex%VolPrb%collection(1)%cordinates(1)%Zi = 2
      ex%VolPrb%collection(1)%cordinates(1)%Ze = 77
      ex%VolPrb%collection(1)%cordinates(1)%or = iExC
      ex%VolPrb%collection(1)%cordinates(1)%xtrancos = 1
      ex%VolPrb%collection(1)%cordinates(1)%ytrancos = 1
      ex%VolPrb%collection(1)%cordinates(1)%ztrancos = 1
      ex%VolPrb%collection(1)%cordinates(1)%tag = ""
      ex%VolPrb%collection(1)%tstart = 0.0
      ex%VolPrb%collection(1)%tstop = 0.0
      ex%VolPrb%collection(1)%tstep = 1e-9
      ex%VolPrb%collection(1)%fstart = 0.0
      ex%VolPrb%collection(1)%fstop = 0.0
      ex%VolPrb%collection(1)%fstep = 0.0
      ex%VolPrb%collection(1)%outputrequest = "electric_field_movie"
      ex%VolPrb%collection(1)%filename = " "
      ex%VolPrb%collection(1)%type2 = NP_T2_TIME
   end function
end function

