integer function test_read_currentinjection() bind (C) result(err)
   use smbjson
   use smbjson_testingTools

   implicit none

   character(len=*),parameter :: filename = PATH_TO_TEST_DATA//INPUT_EXAMPLES//'currentInjection.fdtd.json'
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
      expected%matriz%totalX = 20
      expected%matriz%totalY = 20
      expected%matriz%totalZ = 20

      ! Expected grid.
      expected%despl%nX = 20
      expected%despl%nY = 20
      expected%despl%nZ = 20

      allocate(expected%despl%desX(20))
      allocate(expected%despl%desY(20))
      allocate(expected%despl%desZ(20))
      expected%despl%desX = 0.1
      expected%despl%desY = 0.1
      expected%despl%desZ = 0.1
      expected%despl%mx1 = 0
      expected%despl%my1 = 0
      expected%despl%mz1 = 0
      expected%despl%mx2 = 20
      expected%despl%my2 = 20
      expected%despl%mz2 = 20

      ! Expected boundaries.
      expected%front%tipoFrontera(:) = F_MUR

      ! Expected material regions.
      expected%pecRegs%nVols = 0
      expected%pecRegs%nSurfs = 1
      expected%pecRegs%nLins = 1
      expected%pecRegs%nVols_max = 0
      expected%pecRegs%nSurfs_max = 1
      expected%pecRegs%nLins_max = 1
      allocate(expected%pecRegs%Vols(0))
      allocate(expected%pecRegs%Surfs(1))
      allocate(expected%pecRegs%Lins(1))
      
      ! Body
      expected%pecRegs%Surfs(1)%Or = +iEz
      expected%pecRegs%Surfs(1)%Xi = 5
      expected%pecRegs%Surfs(1)%Xe = 14
      expected%pecRegs%Surfs(1)%Yi = 5
      expected%pecRegs%Surfs(1)%Ye = 14
      expected%pecRegs%Surfs(1)%Zi = 10
      expected%pecRegs%Surfs(1)%Ze = 10
      expected%pecRegs%Surfs(1)%tag = 'aluminum@body'

      ! Exit line
      expected%pecRegs%Lins(1)%Or = +iEy
      expected%pecRegs%Lins(1)%Xi = 10
      expected%pecRegs%Lins(1)%Xe = 10
      expected%pecRegs%Lins(1)%Yi = 15
      expected%pecRegs%Lins(1)%Ye = 19
      expected%pecRegs%Lins(1)%Zi = 10
      expected%pecRegs%Lins(1)%Ze = 10
      expected%pecRegs%Lins(1)%tag = 'aluminum@exit'

      ! Expected sources.
      expected%nodSrc%n_nodSrc = 1
      expected%nodSrc%n_nodSrc_max = 1
      expected%nodSrc%n_C1P_max = 0
      expected%nodSrc%n_C2P_max = 1
      allocate(expected%nodSrc%NodalSource(1))
      expected%nodSrc%NodalSource(1)%nombre = "gauss.exc"
      expected%nodSrc%NodalSource(1)%isElec = .true.
      expected%nodSrc%NodalSource(1)%isHard = .true.
      expected%nodSrc%NodalSource(1)%isInitialValue = .false.
      allocate(expected%nodSrc%NodalSource(1)%c1P(0))
      allocate(expected%nodSrc%NodalSource(1)%c2P(1))
      expected%nodSrc%NodalSource(1)%n_C2P = 1
      expected%nodSrc%NodalSource(1)%c2P(1)%Or = iEy
      expected%nodSrc%NodalSource(1)%c2P(1)%Xi = 10
      expected%nodSrc%NodalSource(1)%c2P(1)%Xe = 10
      expected%nodSrc%NodalSource(1)%c2P(1)%Yi =  0
      expected%nodSrc%NodalSource(1)%c2P(1)%Ye =  4
      expected%nodSrc%NodalSource(1)%c2P(1)%Zi = 10
      expected%nodSrc%NodalSource(1)%c2P(1)%Ze = 10
      expected%nodSrc%NodalSource(1)%c2P(1)%tag = ''
      expected%nodSrc%NodalSource(1)%c2P(1)%xc = 0.0
      expected%nodSrc%NodalSource(1)%c2P(1)%yc = 1.0
      expected%nodSrc%NodalSource(1)%c2P(1)%zc = 0.0

      ! Expected probes
      ! bloqueprobes
      expected%BloquePrb%n_bp = 2
      expected%BloquePrb%n_bp_max = 2
      allocate(expected%BloquePrb%bp(2))
      expected%BloquePrb%bp(1)%outputrequest = "bulk_current_at_entry"
      expected%BloquePrb%bp(1)%FileNormalize = ''
      expected%BloquePrb%bp(1)%type2 = NP_T2_TIME
      expected%BloquePrb%bp(1)%tstart = 0.0
      expected%BloquePrb%bp(1)%tstop = 0.0
      expected%BloquePrb%bp(1)%tstep = 0.0
      expected%BloquePrb%bp(1)%fstart = 0.0
      expected%BloquePrb%bp(1)%fstop = 0.0
      expected%BloquePrb%bp(1)%fstep = 0.0
      expected%BloquePrb%bp(1)%i1 = 9
      expected%BloquePrb%bp(1)%i2 = 10
      expected%BloquePrb%bp(1)%j1 =  2
      expected%BloquePrb%bp(1)%j2 =  2
      expected%BloquePrb%bp(1)%k1 = 9
      expected%BloquePrb%bp(1)%k2 = 10
      expected%BloquePrb%bp(1)%skip = 1
      expected%BloquePrb%bp(1)%nml = iEy
      expected%BloquePrb%bp(1)%t = BcELECT
      expected%BloquePrb%bp(1)%tag = "bulk_current_at_entry"

      expected%BloquePrb%bp(2)%outputrequest = "bulk_current_at_exit"
      expected%BloquePrb%bp(2)%FileNormalize = ' '
      expected%BloquePrb%bp(2)%type2 = NP_T2_TIME
      expected%BloquePrb%bp(2)%tstart = 0.0
      expected%BloquePrb%bp(2)%tstop = 0.0
      expected%BloquePrb%bp(2)%tstep = 0.0
      expected%BloquePrb%bp(2)%fstart = 0.0
      expected%BloquePrb%bp(2)%fstop = 0.0
      expected%BloquePrb%bp(2)%fstep = 0.0
      expected%BloquePrb%bp(2)%i1 = 9
      expected%BloquePrb%bp(2)%i2 = 10
      expected%BloquePrb%bp(2)%j1 = 17
      expected%BloquePrb%bp(2)%j2 = 17
      expected%BloquePrb%bp(2)%k1 = 9
      expected%BloquePrb%bp(2)%k2 = 10
      expected%BloquePrb%bp(2)%skip = 1
      expected%BloquePrb%bp(2)%nml = iEy
      expected%BloquePrb%bp(2)%t = BcELECT
      expected%BloquePrb%bp(2)%tag = "bulk_current_at_exit"
   end function
end function

