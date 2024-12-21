integer function test_read_planewave() bind (C) result(err)
   use smbjson
   use smbjson_testingTools

   implicit none

   character(len=*), parameter :: filename = PATH_TO_TEST_DATA//'cases/planewave.fdtd.json'
   type(Parseador) :: pr, ex
   type(parser_t) :: parser
   logical :: areSame
   err = 0

   ex = expectedProblemDescription()
   parser = parser_t(filename)
   pr = parser%readProblemDescription()
   call expect_eq(err, ex, pr)

contains
   function expectedProblemDescription() result (expected)
      type(Parseador) :: expected

      call initializeProblemDescription(expected)

      ! Expected general info.
      expected%general%dt = 10e-12
      expected%general%nmax = 2000

      ! Excected media matrix.
      expected%matriz%totalX = 10
      expected%matriz%totalY = 10
      expected%matriz%totalZ = 10

      ! Expected grid.
      expected%despl%nX = 10
      expected%despl%nY = 10
      expected%despl%nZ = 10

      allocate(expected%despl%desX(10))
      allocate(expected%despl%desY(10))
      allocate(expected%despl%desZ(10))
      expected%despl%desX = 0.1
      expected%despl%desY = 0.1
      expected%despl%desZ = 0.1
      expected%despl%mx1 = 0
      expected%despl%mx2 = 10
      expected%despl%my1 = 0
      expected%despl%my2 = 10
      expected%despl%mz1 = 0
      expected%despl%mz2 = 10

      ! Expected boundaries.
      expected%front%tipoFrontera(:) = F_MUR

      ! Expected sources.
      allocate(expected%plnSrc%collection(1))
      expected%plnSrc%collection(1)%nombre_fichero = "gauss.exc"
      expected%plnSrc%collection(1)%atributo = ""
      expected%plnSrc%collection(1)%coor1 = [1, 1, 1]
      expected%plnSrc%collection(1)%coor2 = [8, 8, 8]
      expected%plnSrc%collection(1)%theta = 0.0
      expected%plnSrc%collection(1)%phi = 0.0
      expected%plnSrc%collection(1)%alpha = 1.5708
      expected%plnSrc%collection(1)%beta = 0.0
      expected%plnSrc%collection(1)%isRC=.false.
      expected%plnSrc%collection(1)%nummodes=1
      expected%plnSrc%collection(1)%INCERTMAX=0.0
      expected%plnSrc%nc = 1
      expected%plnSrc%nC_max = 1


      ! Expected probes
      ! sonda
      expected%Sonda%len_cor_max = 0
      expected%Sonda%length = 2
      expected%Sonda%length_max = 2
      allocate(expected%Sonda%collection(2))
      expected%Sonda%collection(1)%outputrequest = "electric_field_point_probe"
      expected%Sonda%collection(1)%type1 = NP_T1_PLAIN
      expected%Sonda%collection(1)%type2 = NP_T2_TIME
      expected%Sonda%collection(1)%filename = ' '
      expected%Sonda%collection(1)%tstart = 0.0
      expected%Sonda%collection(1)%tstop = 0.0
      expected%Sonda%collection(1)%tstep = 0.0
      expected%Sonda%collection(1)%fstart = 0.0
      expected%Sonda%collection(1)%fstop = 0.0
      expected%Sonda%collection(1)%fstep = 0.0
      allocate(expected%Sonda%collection(1)%cordinates(3))
      expected%Sonda%collection(1)%len_cor = 3
      expected%Sonda%collection(1)%cordinates(1:3)%tag = "electric_field_point_probe"
      expected%Sonda%collection(1)%cordinates(1:3)%Xi = 4
      expected%Sonda%collection(1)%cordinates(1:3)%Yi = 4
      expected%Sonda%collection(1)%cordinates(1:3)%Zi = 4
      expected%Sonda%collection(1)%cordinates(1)%Or = NP_COR_EX
      expected%Sonda%collection(1)%cordinates(2)%Or = NP_COR_EY
      expected%Sonda%collection(1)%cordinates(3)%Or = NP_COR_EZ

      expected%Sonda%collection(2)%outputrequest = "magnetic_field_point_probe"
      expected%Sonda%collection(2)%type1 = NP_T1_PLAIN
      expected%Sonda%collection(2)%type2 = NP_T2_TIME
      expected%Sonda%collection(2)%filename = ' '
      expected%Sonda%collection(2)%tstart = 0.0
      expected%Sonda%collection(2)%tstop = 0.0
      expected%Sonda%collection(2)%tstep = 0.0
      expected%Sonda%collection(2)%fstart = 0.0
      expected%Sonda%collection(2)%fstop = 0.0
      expected%Sonda%collection(2)%fstep = 0.0
      allocate(expected%Sonda%collection(2)%cordinates(1))
      expected%Sonda%collection(2)%len_cor = 1
      expected%Sonda%collection(2)%cordinates(1)%tag = "magnetic_field_point_probe"
      expected%Sonda%collection(2)%cordinates(1)%Xi = 6
      expected%Sonda%collection(2)%cordinates(1)%Yi = 6
      expected%Sonda%collection(2)%cordinates(1)%Zi = 6
      expected%Sonda%collection(2)%cordinates(1)%Or = NP_COR_HX

   end function
end function

