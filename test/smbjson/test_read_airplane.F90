integer function test_read_airplane() bind (C) result(err)
   use smbjson
   use smbjson_testingTools

   implicit none

   character(len=*),parameter :: filename = PATH_TO_TEST_DATA//INPUT_EXAMPLES//'airplane.fdtd.json'
   type(Parseador) :: problem, expected
   type(parser_t) :: parser
   logical :: areSame
   err = 0

   expected = expectedProblemDescription()
   parser = parser_t(filename)
   problem = parser%readProblemDescription()

   call expect_eq(err, expected, problem, ignoreRegions=.true.)

contains
   function expectedProblemDescription() result (expected)
      type(Parseador) :: expected

      integer :: i

      call initializeProblemDescription(expected)

      ! Expected general info.
      expected%general%dt = 1.2971876e-09
      expected%general%nmax = 23

      ! Excected media matrix.
      expected%matriz%totalX = 50
      expected%matriz%totalY = 50
      expected%matriz%totalZ = 50

      ! Expected grid.
      expected%despl%nX = 50
      expected%despl%nY = 50
      expected%despl%nZ = 50

      allocate(expected%despl%desX(50))
      allocate(expected%despl%desY(50))
      allocate(expected%despl%desZ(50))
      expected%despl%desX = 0.32419496007084553
      expected%despl%desY = 0.12839303226115248
      expected%despl%desZ = 0.3621442456908099
      expected%despl%mx1 = 0
      expected%despl%mx2 = 50
      expected%despl%my1 = 0
      expected%despl%my2 = 50
      expected%despl%mz1 = 0
      expected%despl%mz2 = 50

      ! Expected boundaries.
      expected%front%tipoFrontera(:) = F_PML
      expected%front%propiedadesPML(:)%numCapas = 10
      expected%front%propiedadesPML(:)%orden = 2.0
      expected%front%propiedadesPML(:)%refl = 0.001

      ! Expected sources.
      expected%nodSrc%n_nodSrc = 1
      expected%nodSrc%n_nodSrc_max = 1
      expected%nodSrc%n_C1P_max = 0
      expected%nodSrc%n_C2P_max = 1
      allocate(expected%nodSrc%NodalSource(1))
      expected%nodSrc%NodalSource(1)%nombre = "gauss.exc"
      expected%nodSrc%NodalSource(1)%isElec = .true.
      expected%nodSrc%NodalSource(1)%isMagnet = .false.
      expected%nodSrc%NodalSource(1)%isCurrent = .false.
      expected%nodSrc%NodalSource(1)%isField = .true.
      expected%nodSrc%NodalSource(1)%isInitialValue = .false.
      allocate(expected%nodSrc%NodalSource(1)%c1P(0))
      allocate(expected%nodSrc%NodalSource(1)%c2P(1))
      expected%nodSrc%NodalSource(1)%n_C2P = 1
      expected%nodSrc%NodalSource(1)%c2P(1)%Or = iEz
      expected%nodSrc%NodalSource(1)%c2P(1)%Xi = 5
      expected%nodSrc%NodalSource(1)%c2P(1)%Xe = 5
      expected%nodSrc%NodalSource(1)%c2P(1)%Yi = 30
      expected%nodSrc%NodalSource(1)%c2P(1)%Ye = 30
      expected%nodSrc%NodalSource(1)%c2P(1)%Zi = 39
      expected%nodSrc%NodalSource(1)%c2P(1)%Ze = 46
      expected%nodSrc%NodalSource(1)%c2P(1)%tag = ''
      expected%nodSrc%NodalSource(1)%c2P(1)%xc = 0.0
      expected%nodSrc%NodalSource(1)%c2P(1)%yc = 0.0
      expected%nodSrc%NodalSource(1)%c2P(1)%zc = 1.0
   end function
end function

