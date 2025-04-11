integer function test_read_lumped_fixture() bind (C) result(err)
   use smbjson
   use smbjson_testingTools

   implicit none

   character(len=*),parameter :: filename = PATH_TO_TEST_DATA//INPUT_EXAMPLES//'lumped_fixture.fdtd.json'
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

      ! Expected general info
      expected%general%dt = 7.7033e-12
      expected%general%nmax = 389

      ! Expected media matrix
      expected%matriz%totalX = 20
      expected%matriz%totalY = 8
      expected%matriz%totalZ = 9

      ! Expected grid
      expected%despl%nX = 20
      expected%despl%nY = 8
      expected%despl%nZ = 9

      allocate(expected%despl%desX(20))
      allocate(expected%despl%desY(8))
      allocate(expected%despl%desZ(9))
      expected%despl%desX = 0.005
      expected%despl%desY = 0.005
      expected%despl%desZ = 0.005
      expected%despl%mx1 = 0
      expected%despl%my1 = 0
      expected%despl%mz1 = 0
      expected%despl%mx2 = 20
      expected%despl%my2 = 8
      expected%despl%mz2 = 9

      ! Expected boundaries
      expected%front%tipoFrontera(:) = F_MUR

      ! Expected material regions
      expected%pecRegs%nVols = 0
      expected%pecRegs%nSurfs = 6
      expected%pecRegs%nLins = 0
      expected%pecRegs%nVols_max = 0
      expected%pecRegs%nSurfs_max = 6
      expected%pecRegs%nLins_max = 0
      allocate(expected%pecRegs%Vols(0))
      allocate(expected%pecRegs%Surfs(6))
      allocate(expected%pecRegs%Lins(0))

      ! Left side - Surface 1
      expected%pecRegs%Surfs(1)%Or = +iEx
      expected%pecRegs%Surfs(1)%Xi = 2
      expected%pecRegs%Surfs(1)%Xe = 2
      expected%pecRegs%Surfs(1)%Yi = 2
      expected%pecRegs%Surfs(1)%Ye = 5
      expected%pecRegs%Surfs(1)%Zi = 2
      expected%pecRegs%Surfs(1)%Ze = 6
      expected%pecRegs%Surfs(1)%tag = 'pec@left_side'

      ! Left side - Surface 2
      expected%pecRegs%Surfs(2)%Or = +iEz
      expected%pecRegs%Surfs(2)%Xi = 2
      expected%pecRegs%Surfs(2)%Xe = 8
      expected%pecRegs%Surfs(2)%Yi = 2
      expected%pecRegs%Surfs(2)%Ye = 5
      expected%pecRegs%Surfs(2)%Zi = 2
      expected%pecRegs%Surfs(2)%Ze = 2
      expected%pecRegs%Surfs(2)%tag = 'pec@left_side'

      ! Left side - Surface 3
      expected%pecRegs%Surfs(3)%Or = +iEz
      expected%pecRegs%Surfs(3)%Xi = 2
      expected%pecRegs%Surfs(3)%Xe = 8
      expected%pecRegs%Surfs(3)%Yi = 2
      expected%pecRegs%Surfs(3)%Ye = 5
      expected%pecRegs%Surfs(3)%Zi = 7
      expected%pecRegs%Surfs(3)%Ze = 7
      expected%pecRegs%Surfs(3)%tag = 'pec@left_side'

      ! Right side - Surface 1
      expected%pecRegs%Surfs(4)%Or = +iEz
      expected%pecRegs%Surfs(4)%Xi = 11
      expected%pecRegs%Surfs(4)%Xe = 17
      expected%pecRegs%Surfs(4)%Yi = 2
      expected%pecRegs%Surfs(4)%Ye = 5
      expected%pecRegs%Surfs(4)%Zi = 2
      expected%pecRegs%Surfs(4)%Ze = 2
      expected%pecRegs%Surfs(4)%tag = 'pec@right_side'

      ! Right side - Surface 2
      expected%pecRegs%Surfs(5)%Or = +iEz
      expected%pecRegs%Surfs(5)%Xi = 11
      expected%pecRegs%Surfs(5)%Xe = 17
      expected%pecRegs%Surfs(5)%Yi = 2
      expected%pecRegs%Surfs(5)%Ye = 5
      expected%pecRegs%Surfs(5)%Zi = 7
      expected%pecRegs%Surfs(5)%Ze = 7
      expected%pecRegs%Surfs(5)%tag = 'pec@right_side'

      ! Right side - Surface 3
      expected%pecRegs%Surfs(6)%Or = +iEx
      expected%pecRegs%Surfs(6)%Xi = 18
      expected%pecRegs%Surfs(6)%Xe = 18
      expected%pecRegs%Surfs(6)%Yi = 2
      expected%pecRegs%Surfs(6)%Ye = 5
      expected%pecRegs%Surfs(6)%Zi = 2
      expected%pecRegs%Surfs(6)%Ze = 6
      expected%pecRegs%Surfs(6)%tag = 'pec@right_side'

      ! Expected dielectric regions (including lumped resistor)
      expected%dielRegs%nVols = 0
      expected%dielRegs%nSurfs = 0
      expected%dielRegs%nLins = 1
      expected%dielRegs%nVols_max = 0
      expected%dielRegs%nSurfs_max = 0
      expected%dielRegs%nLins_max = 1
      allocate(expected%dielRegs%Vols(0))
      allocate(expected%dielRegs%Surfs(0))
      allocate(expected%dielRegs%Lins(1))

      ! Lumped resistor line
      allocate(expected%dielRegs%Lins(1)%c1P(0))
      allocate(expected%dielRegs%Lins(1)%c2P(1))
      expected%dielRegs%Lins(1)%n_C1P = 0
      expected%dielRegs%Lins(1)%n_C2P = 1
      expected%dielRegs%Lins(1)%c2P%Or = iEx
      expected%dielRegs%Lins(1)%c2P%Xi = 9
      expected%dielRegs%Lins(1)%c2P%Xe = 10
      expected%dielRegs%Lins(1)%c2P%Yi = 4
      expected%dielRegs%Lins(1)%c2P%Ye = 4
      expected%dielRegs%Lins(1)%c2P%Zi = 7
      expected%dielRegs%Lins(1)%c2P%Ze = 7
      expected%dielRegs%Lins(1)%c2P%tag = '100ohm_resistor@lumped_line'

      expected%dielRegs%Lins(1)%sigma = 0.0
      expected%dielRegs%Lins(1)%eps = EPSILON_VACUUM
      expected%dielRegs%Lins(1)%mu = MU_VACUUM
      expected%dielRegs%Lins(1)%sigmam = 0.0

      expected%dielRegs%Lins(1)%R = 100.0
      expected%dielRegs%Lins(1)%Rtime_on = 0.0
      expected%dielRegs%Lins(1)%Rtime_off = 1.0

      expected%dielRegs%Lins(1)%resistor = .true.

      
      ! Expected sources
      expected%nodSrc%n_nodSrc = 1
      expected%nodSrc%n_nodSrc_max = 1
      expected%nodSrc%n_C1P_max = 0
      expected%nodSrc%n_C2P_max = 1
      allocate(expected%nodSrc%NodalSource(1))
      expected%nodSrc%NodalSource(1)%nombre = "predefinedExcitation.1.exc"
      expected%nodSrc%NodalSource(1)%isElec = .true.
      expected%nodSrc%NodalSource(1)%isHard = .false.
      expected%nodSrc%NodalSource(1)%isInitialValue = .false.
      allocate(expected%nodSrc%NodalSource(1)%c1P(0))
      allocate(expected%nodSrc%NodalSource(1)%c2P(1))
      expected%nodSrc%NodalSource(1)%n_C2P = 1
      expected%nodSrc%NodalSource(1)%c2P(1)%Or = iEx
      expected%nodSrc%NodalSource(1)%c2P(1)%Xi = 9
      expected%nodSrc%NodalSource(1)%c2P(1)%Xe = 10
      expected%nodSrc%NodalSource(1)%c2P(1)%Yi = 4
      expected%nodSrc%NodalSource(1)%c2P(1)%Ye = 4
      expected%nodSrc%NodalSource(1)%c2P(1)%Zi = 2
      expected%nodSrc%NodalSource(1)%c2P(1)%Ze = 2
      expected%nodSrc%NodalSource(1)%c2P(1)%tag = 'nodal_source'
      expected%nodSrc%NodalSource(1)%c2P(1)%xc = 1.0
      expected%nodSrc%NodalSource(1)%c2P(1)%yc = 0.0
      expected%nodSrc%NodalSource(1)%c2P(1)%zc = 0.0

      ! Expected probes
      
      ! Electric field point probe
      expected%Sonda%length = 1
      expected%Sonda%length_max = 1
      allocate(expected%Sonda%collection(1))

      expected%Sonda%collection(1)%outputrequest = "e_probe"
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
      expected%Sonda%collection(1)%cordinates(1:3)%Xi = 10
      expected%Sonda%collection(1)%cordinates(1:3)%Yi = 3
      expected%Sonda%collection(1)%cordinates(1:3)%Zi = 7
      expected%Sonda%collection(1)%cordinates(1)%Or = NP_COR_EX
      expected%Sonda%collection(1)%cordinates(2)%Or = NP_COR_EY
      expected%Sonda%collection(1)%cordinates(3)%Or = NP_COR_EZ
      expected%Sonda%collection(1)%cordinates(1:3)%tag = "e_probe"

      ! Bulk current probe
      expected%BloquePrb%n_bp = 1
      expected%BloquePrb%n_bp_max = 1
      allocate(expected%BloquePrb%bp(1))
      
      expected%BloquePrb%bp(1)%outputrequest = "Bulk probe"
      expected%BloquePrb%bp(1)%FileNormalize = ' '
      expected%BloquePrb%bp(1)%type2 = NP_T2_TIME
      expected%BloquePrb%bp(1)%tstart = 0.0
      expected%BloquePrb%bp(1)%tstop = 0.0
      expected%BloquePrb%bp(1)%tstep = 0.0
      expected%BloquePrb%bp(1)%fstart = 0.0
      expected%BloquePrb%bp(1)%fstop = 0.0
      expected%BloquePrb%bp(1)%fstep = 0.0
      expected%BloquePrb%bp(1)%i1 = 6
      expected%BloquePrb%bp(1)%i2 = 6
      expected%BloquePrb%bp(1)%j1 = 1
      expected%BloquePrb%bp(1)%j2 = 6
      expected%BloquePrb%bp(1)%k1 = 6
      expected%BloquePrb%bp(1)%k2 = 7
      expected%BloquePrb%bp(1)%skip = 1
      expected%BloquePrb%bp(1)%nml = iEx
      expected%BloquePrb%bp(1)%t = BcELECT
      expected%BloquePrb%bp(1)%tag = "Bulk probe"

   end function
end function 