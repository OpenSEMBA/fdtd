integer function test_read_thinSlot() bind (C) result(err)
   use smbjson
   use smbjson_testingTools

   implicit none

   character(len=*), parameter :: filename = PATH_TO_TEST_DATA//'cases/thinSlot.fdtd.json'
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
      integer :: i

      call initializeProblemDescription(expected)

      ! Expected general info.
      expected%general%dt = 10e-12
      expected%general%nmax = 2000

      ! Excected media matrix.
      expected%matriz%totalX = 4
      expected%matriz%totalY = 4
      expected%matriz%totalZ = 50

      ! Expected grid.
      expected%despl%nX = 4
      expected%despl%nY = 4
      expected%despl%nZ = 50

      allocate(expected%despl%desX(4))
      allocate(expected%despl%desY(4))
      allocate(expected%despl%desZ(50))
      expected%despl%desX = 0.1
      expected%despl%desY = 0.1
      expected%despl%desZ = 0.1
      expected%despl%mx1 = 0
      expected%despl%mx2 = 4
      expected%despl%my1 = 0
      expected%despl%my2 = 4
      expected%despl%mz1 = 0
      expected%despl%mz2 = 50

      ! Expected boundaries.
      expected%front%tipoFrontera(F_XL) = F_PEC
      expected%front%tipoFrontera(F_XU) = F_PEC
      expected%front%tipoFrontera(F_YL) = F_PMC
      expected%front%tipoFrontera(F_YU) = F_PMC
      expected%front%tipoFrontera(F_ZL) = F_PML
      expected%front%tipoFrontera(F_ZU) = F_PML

      ! Expected sources.
      allocate(expected%plnSrc%collection(1))
      expected%plnSrc%collection(1)%nombre_fichero = "gauss.exc"
      expected%plnSrc%collection(1)%atributo = ""
      expected%plnSrc%collection(1)%coor1 = [0, 0, 2]
      expected%plnSrc%collection(1)%coor2 = [3, 3, 47]
      expected%plnSrc%collection(1)%theta = 0.0
      expected%plnSrc%collection(1)%phi = 0.0
      expected%plnSrc%collection(1)%alpha = 1.5708
      expected%plnSrc%collection(1)%beta = 0.0
      expected%plnSrc%collection(1)%isRC=.false.
      expected%plnSrc%collection(1)%nummodes=1
      expected%plnSrc%collection(1)%INCERTMAX=0.0
      expected%plnSrc%nc = 1
      expected%plnSrc%nC_max = 1

      ! materials
      !! pec square
      expected%pecRegs%nVols = 0
      expected%pecRegs%nSurfs = 1
      expected%pecRegs%nLins = 0
      expected%pecRegs%nVols_max = 0
      expected%pecRegs%nSurfs_max = 1
      expected%pecRegs%nLins_max = 0
      allocate(expected%pecRegs%Vols(0))
      allocate(expected%pecRegs%Surfs(1))
      allocate(expected%pecRegs%Lins(0))
      expected%pecRegs%Surfs(1)%Or = +iEz
      expected%pecRegs%Surfs(1)%Xi = 0
      expected%pecRegs%Surfs(1)%Xe = 3
      expected%pecRegs%Surfs(1)%Yi = 0
      expected%pecRegs%Surfs(1)%Ye = 3
      expected%pecRegs%Surfs(1)%Zi = 25
      expected%pecRegs%Surfs(1)%Ze = 25
      expected%pecRegs%Surfs(1)%tag = 'copper@square'
      
      !! thin slot
      expected%tSlots%n_tg = 1
      expected%tSlots%n_tg_max = 1
      allocate(expected%tslots%tg(1))
      expected%tSlots%tg(1)%width = 3e-3
      expected%tSlots%tg(1)%n_tgc = 2
      expected%tSlots%tg(1)%n_tgc_max = 2
      allocate(expected%tSlots%tg(1)%tgc(2))
      expected%tSlots%tg(1)%tgc(1)%i = 1
      expected%tSlots%tg(1)%tgc(1)%j = 2
      expected%tSlots%tg(1)%tgc(1)%k = 25
      expected%tSlots%tg(1)%tgc(1)%node = 0
      expected%tSlots%tg(1)%tgc(1)%dir = iEx
      expected%tSlots%tg(1)%tgc(1)%Or = -1
      expected%tSlots%tg(1)%tgc(1)%tag = "3mm-gap@slot"
      expected%tSlots%tg(1)%tgc(2) = expected%tSlots%tg(1)%tgc(1)
      expected%tSlots%tg(1)%tgc(2)%i = 2

      ! Expected probes
      ! sonda
      expected%Sonda%len_cor_max = 0
      expected%Sonda%length = 2
      expected%Sonda%length_max = 2
      allocate(expected%Sonda%collection(2))
      ! common data
      do i = 1, 2
         expected%Sonda%collection(i)%type1 = NP_T1_PLAIN
         expected%Sonda%collection(i)%type2 = NP_T2_TIME
         expected%Sonda%collection(i)%filename = ' '
         expected%Sonda%collection(i)%tstart = 0.0
         expected%Sonda%collection(i)%tstop = 0.0
         expected%Sonda%collection(i)%tstep = 0.0
         expected%Sonda%collection(i)%fstart = 0.0
         expected%Sonda%collection(i)%fstop = 0.0
         expected%Sonda%collection(i)%fstep = 0.0
         allocate(expected%Sonda%collection(i)%cordinates(3))
         expected%Sonda%collection(i)%cordinates(1)%Or = NP_COR_EX
         expected%Sonda%collection(i)%cordinates(2)%Or = NP_COR_EY
         expected%Sonda%collection(i)%cordinates(3)%Or = NP_COR_EZ
         expected%Sonda%collection(i)%len_cor = 3
      end do
      ! point probe at front
      expected%Sonda%collection(1)%outputrequest = "front"
      expected%Sonda%collection(1)%cordinates(1:3)%tag = "front"
      expected%Sonda%collection(1)%cordinates(1:3)%Xi = 2
      expected%Sonda%collection(1)%cordinates(1:3)%Yi = 2
      expected%Sonda%collection(1)%cordinates(1:3)%Zi = 10
      ! point probe at back
      expected%Sonda%collection(2)%outputrequest = "back"
      expected%Sonda%collection(2)%cordinates(1:3)%tag = "back"
      expected%Sonda%collection(2)%cordinates(1:3)%Xi = 2
      expected%Sonda%collection(2)%cordinates(1:3)%Yi = 2
      expected%Sonda%collection(2)%cordinates(1:3)%Zi = 40
      
   end function
end function

