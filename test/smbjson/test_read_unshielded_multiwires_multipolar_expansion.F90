integer function test_read_unshielded_multiwires_multipolar_expansion() bind (C) result(err)
   use smbjson
   use smbjson_testingTools

   implicit none

   character(len=*),parameter :: filename = &
      PATH_TO_TEST_DATA//INPUT_EXAMPLES//'unshielded_multiwires_multipolar_expansion.fdtd.json'
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
      ex%general%dt = 3e-11
      ex%general%nmax = 1100

      ! Excected media matrix.
      ex%matriz%totalX = 30
      ex%matriz%totalY = 30
      ex%matriz%totalZ = 30

      ! Expected grid.
      ex%despl%nX = 30
      ex%despl%nY = 30
      ex%despl%nZ = 30

      allocate(ex%despl%desX(30))
      allocate(ex%despl%desY(30))
      allocate(ex%despl%desZ(30))
      ex%despl%desX = 0.2
      ex%despl%desY = 0.2
      ex%despl%desZ = 0.2
      ex%despl%mx1 = 0
      ex%despl%mx2 = 30
      ex%despl%my1 = 0
      ex%despl%my2 = 30
      ex%despl%mz1 = 0
      ex%despl%mz2 = 30

      ! Expected boundaries.
      ex%front%tipoFrontera(:) = F_PML
      ex%front%propiedadesPML(:)%numCapas = 8
      ex%front%propiedadesPML(:)%orden = 2
      ex%front%propiedadesPML(:)%refl = 0.001

      ! Expected material regions.
      ex%pecRegs%nSurfs = 1
      ex%pecRegs%nSurfs_max = 1
      allocate(ex%pecRegs%Surfs(1))
      ! -- specific surfs not included DO NOT use comparison --

      ! Expected sources.
      allocate(ex%plnSrc%collection(1))
      ex%plnSrc%collection(1)%nombre_fichero = "unshielded_50ns.exc"
      ex%plnSrc%collection(1)%atributo = ""
      ex%plnSrc%collection(1)%coor1 = [1, 1, 1]
      ex%plnSrc%collection(1)%coor2 = [28, 28, 28]
      ex%plnSrc%collection(1)%theta = 1.5708
      ex%plnSrc%collection(1)%phi = 0.0
      ex%plnSrc%collection(1)%alpha = 0.0
      ex%plnSrc%collection(1)%beta = 0.0
      ex%plnSrc%collection(1)%isRC=.false.
      ex%plnSrc%collection(1)%nummodes=1
      ex%plnSrc%collection(1)%INCERTMAX=0.0
      ex%plnSrc%nc = 1
      ex%plnSrc%nC_max = 1

      ! Expected Probe
      deallocate(expected%mtln%probes)
      allocate(expected%mtln%probes(1))
      expected%mtln%probes(1)%attached_to_cable => expected%mtln%cables(1)
      expected%mtln%probes(1)%index = 5
      expected%mtln%probes(1)%probe_type = PROBE_TYPE_CURRENT
      expected%mtln%probes(1)%probe_name = "test"
      expected%mtln%probes(1)%probe_position = [2,11,14]
   end function
end function

