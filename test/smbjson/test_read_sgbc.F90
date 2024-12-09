integer function test_read_sgbc() bind (C) result(err)
   use smbjson
   use smbjson_testingTools

   implicit none

   character(len=*), parameter :: filename = PATH_TO_TEST_DATA//'cases/sgbc.fdtd.json'
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

      ! Expected materials
      !! PECs
      expected%pecRegs%nSurfs = 1
      expected%pecRegs%nLins = 0
      expected%pecRegs%nVols_max = 0
      expected%pecRegs%nSurfs_max = 1
      expected%pecRegs%nLins_max = 0
      allocate(expected%pecRegs%Vols(0))
      allocate(expected%pecRegs%Surfs(1))
      
      !!! 2x2 PEC square
      expected%pecRegs%Surfs(1)%Or = +iEz
      expected%pecRegs%Surfs(1)%Xi = 3
      expected%pecRegs%Surfs(1)%Xe = 4
      expected%pecRegs%Surfs(1)%Yi = 3
      expected%pecRegs%Surfs(1)%Ye = 4
      expected%pecRegs%Surfs(1)%Zi = 3
      expected%pecRegs%Surfs(1)%Ze = 3
      expected%pecRegs%Surfs(1)%tag = 'material1@layer1'

      !! Composites
      allocate(expected%lossyThinSurfs%cs(2))
      expected%lossyThinSurfs%length = 2
      expected%lossyThinSurfs%length_max = 2
      expected%lossyThinSurfs%nC_max = 2
      
      !!! 2-layer composite
      allocate(expected%lossyThinSurfs%cs(1)%c(1))
      expected%lossyThinSurfs%cs(1)%nc = 1
      expected%lossyThinSurfs%cs(1)%c(1)%tag = '2-layers-composite@layer2'
      expected%lossyThinSurfs%cs(1)%c(1)%Or = +iEx
      expected%lossyThinSurfs%cs(1)%c(1)%Xi = 3
      expected%lossyThinSurfs%cs(1)%c(1)%Xe = 4
      expected%lossyThinSurfs%cs(1)%c(1)%Yi = 3
      expected%lossyThinSurfs%cs(1)%c(1)%Ye = 3
      expected%lossyThinSurfs%cs(1)%c(1)%Zi = 3
      expected%lossyThinSurfs%cs(1)%c(1)%Ze = 4
      expected%lossyThinSurfs%cs(1)%numcapas = 2
      allocate(expected%lossyThinSurfs%cs(1)%thk(2))
      allocate(expected%lossyThinSurfs%cs(1)%sigma(2))
      allocate(expected%lossyThinSurfs%cs(1)%eps(2))
      allocate(expected%lossyThinSurfs%cs(1)%mu(2))
      allocate(expected%lossyThinSurfs%cs(1)%sigmam(2))
      expected%lossyThinSurfs%cs(1)%thk    = [              1e-3,               5e-3]
      expected%lossyThinSurfs%cs(1)%sigma  = [              2e-4,                0.0]
      expected%lossyThinSurfs%cs(1)%eps    = [1.3*EPSILON_VACUUM, 1.3*EPSILON_VACUUM]
      expected%lossyThinSurfs%cs(1)%mu     = [         MU_VACUUM,          MU_VACUUM]
      expected%lossyThinSurfs%cs(1)%sigmam = [               0.0,                0.0]
      
      !!! 3-layer composite
      allocate(expected%lossyThinSurfs%cs(2)%c(1))
      expected%lossyThinSurfs%cs(2)%nc = 1
      expected%lossyThinSurfs%cs(2)%c(1)%tag = '3-layers-composite@layer3'
      expected%lossyThinSurfs%cs(2)%c(1)%Or = +iEy
      expected%lossyThinSurfs%cs(2)%c(1)%Xi = 3
      expected%lossyThinSurfs%cs(2)%c(1)%Xe = 3
      expected%lossyThinSurfs%cs(2)%c(1)%Yi = 3
      expected%lossyThinSurfs%cs(2)%c(1)%Ye = 4
      expected%lossyThinSurfs%cs(2)%c(1)%Zi = 3
      expected%lossyThinSurfs%cs(2)%c(1)%Ze = 4
      expected%lossyThinSurfs%cs(2)%numcapas = 3
      allocate(expected%lossyThinSurfs%cs(2)%thk(3))
      allocate(expected%lossyThinSurfs%cs(2)%sigma(3))
      allocate(expected%lossyThinSurfs%cs(2)%eps(3))
      allocate(expected%lossyThinSurfs%cs(2)%mu(3))
      allocate(expected%lossyThinSurfs%cs(2)%sigmam(3))
      expected%lossyThinSurfs%cs(2)%thk    = [          1e-3,           5e-3,           1e-3]
      expected%lossyThinSurfs%cs(2)%sigma  = [          2e-4,            0.0,            0.0]
      expected%lossyThinSurfs%cs(2)%eps    = [EPSILON_VACUUM, EPSILON_VACUUM, EPSILON_VACUUM]
      expected%lossyThinSurfs%cs(2)%mu     = [     MU_VACUUM,  1.3*MU_VACUUM,      MU_VACUUM]
      expected%lossyThinSurfs%cs(2)%sigmam = [           0.0,            0.0,           1e-4]
   end function
end function

