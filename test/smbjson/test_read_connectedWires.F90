integer function test_read_connectedwires() bind (C) result(err)
   use smbjson
   use smbjson_testingTools

   implicit none

   character(len=*),parameter :: filename = PATH_TO_TEST_DATA//'cases/connectedWires.fdtd.json'
   type(Parseador) :: problem, expected
   type(parser_t) :: parser
   logical :: areSame
   err = 0

   expected = expectedProblemDescription()
   parser = parser_t(filename)
   problem = parser%readProblemDescription()
   call expect_eq(err, expected, problem)

contains
   function expectedProblemDescription() result (expected)
      type(Parseador) :: expected

      integer :: i

      call initializeProblemDescription(expected)

      ! Expected general info.
      expected%general%dt = 1e-12
      expected%general%nmax = 1000

      ! Excected media matrix.
      expected%matriz%totalX = 60
      expected%matriz%totalY = 60
      expected%matriz%totalZ = 60

      ! Expected grid.
      expected%despl%nX = 60
      expected%despl%nY = 60
      expected%despl%nZ = 60

      allocate(expected%despl%desX(60))
      allocate(expected%despl%desY(60))
      allocate(expected%despl%desZ(60))
      expected%despl%desX = 0.01
      expected%despl%desY = 0.01
      expected%despl%desZ = 0.01
      expected%despl%mx1 = 0
      expected%despl%mx2 = 60
      expected%despl%my1 = 0
      expected%despl%my2 = 60
      expected%despl%mz1 = 0
      expected%despl%mz2 = 60

      ! Expected boundaries.
      expected%front%tipoFrontera(:) = F_PML
      expected%front%propiedadesPML(:)%numCapas = 6
      expected%front%propiedadesPML(:)%orden = 2.0
      expected%front%propiedadesPML(:)%refl = 0.001

      ! Expected PEC regions.
      expected%pecRegs%nLins = 0
      expected%pecRegs%nLins_max = 0
      expected%pecRegs%nSurfs = 1
      expected%pecRegs%nSurfs_max = 1
      expected%pecRegs%nVols = 0
      expected%pecRegs%nVols_max = 0
      allocate(expected%pecRegs%Surfs(1))
      expected%pecRegs%Surfs(1)%Xi = 25
      expected%pecRegs%Surfs(1)%Xe = 44
      expected%pecRegs%Surfs(1)%yi = 20
      expected%pecRegs%Surfs(1)%ye = 29
      expected%pecRegs%Surfs(1)%zi = 30
      expected%pecRegs%Surfs(1)%ze = 30
      expected%pecRegs%Surfs(1)%Xtrancos = 1
      expected%pecRegs%Surfs(1)%Ytrancos = 1
      expected%pecRegs%Surfs(1)%Ztrancos = 1
      expected%pecRegs%Surfs(1)%Or = 3
      expected%pecRegs%Surfs(1)%tag =  trim(adjustl(" "))

      ! Expected probes
      ! sonda
      expected%Sonda%length = 2
      expected%Sonda%length_max = 2
      allocate(expected%Sonda%collection(2))
      expected%Sonda%collection(1)%outputrequest = "wire_start"
      expected%Sonda%collection(1)%type1 = NP_T1_PLAIN
      expected%Sonda%collection(1)%type2 = NP_T2_TIME
      expected%Sonda%collection(1)%filename = ' '
      expected%Sonda%collection(1)%tstart = 0.0
      expected%Sonda%collection(1)%tstop = 0.0
      expected%Sonda%collection(1)%tstep = 0.0
      expected%Sonda%collection(1)%fstart = 0.0
      expected%Sonda%collection(1)%fstop = 0.0
      expected%Sonda%collection(1)%fstep = 0.0
      allocate(expected%Sonda%collection(1)%cordinates(1))
      expected%Sonda%collection(1)%len_cor = 1
      expected%Sonda%collection(1)%cordinates(1)%tag = "wire_start"
      expected%Sonda%collection(1)%cordinates(1)%Xi = 1 ! Coord id as tag.
      expected%Sonda%collection(1)%cordinates(1)%Yi = 0
      expected%Sonda%collection(1)%cordinates(1)%Zi = 0
      expected%Sonda%collection(1)%cordinates(1)%Or = NP_COR_WIRECURRENT
      
      expected%Sonda%collection(2)%outputrequest = "wire_end"
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
      expected%Sonda%collection(2)%cordinates(1)%tag = "wire_end"
      expected%Sonda%collection(2)%cordinates(1)%Xi = 4 ! Coord id as tag.
      expected%Sonda%collection(2)%cordinates(1)%Yi = 0
      expected%Sonda%collection(2)%cordinates(1)%Zi = 0
      expected%Sonda%collection(2)%cordinates(1)%Or = NP_COR_WIRECURRENT
      
      
      ! Expected thin wires
      allocate(expected%tWires%tw(2))
      ! thin wire 1
      expected%tWires%tw(1)%rad=0.1e-3
      expected%tWires%tw(1)%dispfile = trim(adjustl(" "))
      expected%tWires%tw(1)%dispfile_LeftEnd = trim(adjustl(" "))
      expected%tWires%tw(1)%dispfile_RightEnd = trim(adjustl(" "))
      expected%tWires%tw(1)%n_twc=10
      expected%tWires%tw(1)%n_twc_max=10
      allocate(expected%tWires%tw(1)%twc(10))
      expected%tWires%tw(1)%twc(1)%srcfile = 'ramp.exc'
      expected%tWires%tw(1)%twc(1)%srctype = 'VOLT'
      expected%tWires%tw(1)%twc(1)%m = 1.0
      expected%tWires%tw(1)%twc(2)%srcfile = 'None'
      expected%tWires%tw(1)%twc(2)%srctype = 'None'
      expected%tWires%tw(1)%twc(2)%m = 0.0
      expected%tWires%tw(1)%twc(1:2)%i = 27
      expected%tWires%tw(1)%twc(1:2)%j = 25
      expected%tWires%tw(1)%twc(1:2)%k = [(i, i = 30, 31)]
      expected%tWires%tw(1)%twc(1:2)%d = DIR_Z
      expected%tWires%tw(1)%twc(1:2)%nd = -1

      expected%tWires%tw(1)%twc(3:10)%srcfile = 'None'
      expected%tWires%tw(1)%twc(3:10)%srctype = 'None'
      expected%tWires%tw(1)%twc(3:10)%m = 0.0
      expected%tWires%tw(1)%twc(3:10)%i = [(i, i=27,34)]
      expected%tWires%tw(1)%twc(3:10)%j = 25
      expected%tWires%tw(1)%twc(3:10)%k = 32
      expected%tWires%tw(1)%twc(3:10)%d = DIR_X
      expected%tWires%tw(1)%twc(3:10)%nd = -1

      expected%tWires%tw(1)%twc(1)%nd  = 1
      expected%tWires%tw(1)%twc(3)%nd  = 2
      expected%tWires%tw(1)%twc(10)%nd  = 5
      
      expected%tWires%tw(1)%twc(1:10)%tag = trim(adjustl("2"))   ! The polyline id is used as tag.
      
      expected%tWires%tw(1)%tl = MATERIAL_CONS
      expected%tWires%tw(1)%C_LeftEnd = 1e22
      expected%tWires%tw(1)%tr = SERIES_CONS
      expected%tWires%tw(1)%R_RightEnd = 50
      expected%tWires%tw(1)%C_RightEnd = 1e22

      ! thin wire 2

      expected%tWires%tw(2)%rad=0.1e-3
      expected%tWires%tw(2)%dispfile = trim(adjustl(" "))
      expected%tWires%tw(2)%dispfile_LeftEnd = trim(adjustl(" "))
      expected%tWires%tw(2)%dispfile_RightEnd = trim(adjustl(" "))
      expected%tWires%tw(2)%n_twc=10
      expected%tWires%tw(2)%n_twc_max=10
      allocate(expected%tWires%tw(2)%twc(10))
      expected%tWires%tw(2)%twc(1:8)%srcfile = 'None'
      expected%tWires%tw(2)%twc(1:8)%srctype = 'None'
      expected%tWires%tw(2)%twc(1:8)%i = [(i, i=35,42)]
      expected%tWires%tw(2)%twc(1:8)%j = 25
      expected%tWires%tw(2)%twc(1:8)%k = 32
      expected%tWires%tw(2)%twc(1:8)%d = DIR_X
      expected%tWires%tw(2)%twc(1:8)%nd = -1
      

      expected%tWires%tw(2)%twc(9:10)%srcfile = 'None'
      expected%tWires%tw(2)%twc(9:10)%srctype = 'None'
      expected%tWires%tw(2)%twc(9:10)%m = 0.0
      expected%tWires%tw(2)%twc(9:10)%i = 43
      expected%tWires%tw(2)%twc(9:10)%j = 25
      expected%tWires%tw(2)%twc(9:10)%k = [(i, i=31,30,-1)]
      expected%tWires%tw(2)%twc(9:10)%d = DIR_Z
      expected%tWires%tw(2)%twc(9:10)%nd = -1


      expected%tWires%tw(2)%twc(1)%nd  = 5
      expected%tWires%tw(2)%twc(9)%nd  = 3
      expected%tWires%tw(2)%twc(10)%nd  = 4
      
      expected%tWires%tw(2)%twc(1:10)%tag = trim(adjustl("5"))   ! The polyline id is used as tag.
      
      expected%tWires%tw(2)%tl = MATERIAL_CONS
      expected%tWires%tw(2)%C_LeftEnd = 1e22
      expected%tWires%tw(2)%tr = MATERIAL_CONS
      expected%tWires%tw(2)%C_RightEnd = 1e22

      expected%tWires%n_tw = 2
      expected%tWires%n_tw_max = 2
   end function
end function

