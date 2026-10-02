integer function test_read_maloneysheet() bind (C) result(err)
   use smbjson_m
   use smbjson_testingTools

   implicit none

   character(len=*), parameter :: filename = PATH_TO_TEST_DATA//INPUT_EXAMPLES//'maloneySheet.fdtd.json'
   type(Parseador_t) :: pr, ex
   type(parser_t) :: parser
   logical :: areSame
   err = 0

   ex = expectedProblemDescription()
   parser = parser_t(filename)
   pr = parser%readProblemDescription()
   call expect_eq(err, ex, pr)

contains
   function expectedProblemDescription() result (expected)
      type(Parseador_t) :: expected

      call initializeProblemDescription(expected)

      ! Expected general info.
      expected%general%dt = 10e-12_RKIND
      expected%general%nmax = 2000

      ! Excected media matrix.
      expected%matriz%totalX = 11
      expected%matriz%totalY = 11
      expected%matriz%totalZ = 11

      ! Expected grid.
      expected%despl%nX = 1
      expected%despl%nY = 1
      expected%despl%nZ = 1

      allocate(expected%despl%desX(1:1))
      allocate(expected%despl%desY(1:1))
      allocate(expected%despl%desZ(1:1))
      expected%despl%desX = 0.1_RKIND
      expected%despl%desY = 0.1_RKIND
      expected%despl%desZ = 0.1_RKIND
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

      !! Maloney sheets
      allocate(expected%maloneySheets%cs(2))
      expected%maloneySheets%length = 2
      expected%maloneySheets%length_max = 2
      expected%maloneySheets%nC_max = 1

      !!! sheet-1
      allocate(expected%maloneySheets%cs(1)%c(1))
      expected%maloneySheets%cs(1)%nc = 1
      expected%maloneySheets%cs(1)%files = 'sheet-1'
      expected%maloneySheets%cs(1)%c(1)%tag = 'sheet-1@layer2'
      expected%maloneySheets%cs(1)%c(1)%Or = +iEy
      expected%maloneySheets%cs(1)%c(1)%Xi = 3
      expected%maloneySheets%cs(1)%c(1)%Xe = 4
      expected%maloneySheets%cs(1)%c(1)%Yi = 3
      expected%maloneySheets%cs(1)%c(1)%Ye = 3
      expected%maloneySheets%cs(1)%c(1)%Zi = 3
      expected%maloneySheets%cs(1)%c(1)%Ze = 4
      expected%maloneySheets%cs(1)%thk = 1e-3_RKIND
      expected%maloneySheets%cs(1)%sigma = 2e-4_RKIND
      expected%maloneySheets%cs(1)%eps = 1.3_RKIND*EPSILON_VACUUM
      expected%maloneySheets%cs(1)%mu = MU_VACUUM
      expected%maloneySheets%cs(1)%sigmam = 0.0_RKIND

      !!! sheet-2
      allocate(expected%maloneySheets%cs(2)%c(1))
      expected%maloneySheets%cs(2)%nc = 1
      expected%maloneySheets%cs(2)%files = 'sheet-2'
      expected%maloneySheets%cs(2)%c(1)%tag = 'sheet-2@layer3'
      expected%maloneySheets%cs(2)%c(1)%Or = +iEx
      expected%maloneySheets%cs(2)%c(1)%Xi = 3
      expected%maloneySheets%cs(2)%c(1)%Xe = 3
      expected%maloneySheets%cs(2)%c(1)%Yi = 3
      expected%maloneySheets%cs(2)%c(1)%Ye = 4
      expected%maloneySheets%cs(2)%c(1)%Zi = 3
      expected%maloneySheets%cs(2)%c(1)%Ze = 4
      expected%maloneySheets%cs(2)%thk = 5e-3_RKIND
      expected%maloneySheets%cs(2)%sigma = 0.0_RKIND
      expected%maloneySheets%cs(2)%eps = 1.5e-11_RKIND
      expected%maloneySheets%cs(2)%mu = MU_VACUUM
      expected%maloneySheets%cs(2)%sigmam = 0.0_RKIND
   end function
end function

