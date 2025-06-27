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
      class(cable_t), pointer :: ptr
      integer :: i

      call initializeProblemDescription(ex)

      ! Expected general info.
      ex%general%dt = 3e-11
      ex%general%nmax = 1100

      ! Expected media matrix.
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

      ! ex sources.
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

      ! ex mtln type
      ex%mtln%has_multiwires = .true.
      ex%mtln%time_step = 3e-11
      ex%mtln%number_of_steps = 1100
      deallocate(ex%mtln%cables)
      allocate(ex%mtln%cables(1))
      
      allocate(unshielded_multiwire_t :: ex%mtln%cables(1)%ptr)
      call initializeCablePULParameters(ex%mtln%cables(1)%ptr, 2)
      ptr => ex%mtln%cables(1)%ptr
      select type(ptr)
      type is(unshielded_multiwire_t)
      ! cable 1 - wire
         ptr%name = "unshielded_pair"

         deallocate(ptr%multipolar_expansion)
         allocate(ptr%multipolar_expansion(1))
         ptr%multipolar_expansion(1)%inner_region%min = &
            [-0.0265000002, -0.0310000002] 
         ptr%multipolar_expansion(1)%inner_region%max =  &
            [ 0.03550000020000001, 0.0310000002] 
         allocate(ptr%multipolar_expansion(1)%electric(2))
         ! First conductor.
         allocate(ptr%multipolar_expansion(1)%electric(1)%conductor_potentials(2))
         ptr%multipolar_expansion(1)%electric(1)%conductor_potentials = [ 1.0, 0.5909272203987278 ]
         ptr%multipolar_expansion(1)%electric(1)%expansion_center = [ -0.004970886788455953, 6.610694092023349e-07 ]
         ptr%multipolar_expansion(1)%electric(1)%inner_region_average_potential = 0.5608636261599323
         allocate(ptr%multipolar_expansion(1)%electric(1)%ab(1))
         ptr%multipolar_expansion(1)%electric(1)%ab(1)%a = 0.9488836986256424
         ptr%multipolar_expansion(1)%electric(1)%ab(1)%b = 0.0
         ! Second conductor
         allocate(ptr%multipolar_expansion(1)%electric(2)%conductor_potentials(2))
         ptr%multipolar_expansion(1)%electric(2)%conductor_potentials = [  0.8497110567446987, 1.0 ]
         ptr%multipolar_expansion(1)%electric(2)%expansion_center = [ 0.009920513440028656, 6.949869591535922e-07 ]
         ptr%multipolar_expansion(1)%electric(2)%inner_region_average_potential = 0.8070848243572611
         allocate(ptr%multipolar_expansion(1)%electric(2)%ab(1))
         ptr%multipolar_expansion(1)%electric(2)%ab(1)%a = 1.3644011168458479
         ptr%multipolar_expansion(1)%electric(2)%ab(1)%b = 0.0
         
         allocate(ptr%multipolar_expansion(1)%magnetic(2))
         ptr%multipolar_expansion(1)%magnetic = &
            ptr%multipolar_expansion(1)%electric

         allocate(ptr%step_size(15))
         ptr%step_size(:) =  0.2

         allocate(ptr%segments(15))
         ptr%segments(:)%x = 2
         ptr%segments(:)%y = 11
         do i = 1, 15
            ptr%segments(i)%z = 6+i
            ptr%segments(i)%orientation = DIRECTION_Z_POS
         end do


         ptr%initial_connector => null()
         ptr%end_connector => null()
      end select
      ! Expected Probe
      deallocate(ex%mtln%probes)
      allocate(ex%mtln%probes(1))
      ex%mtln%probes(1)%attached_to_cable => ex%mtln%cables(1)%ptr
      ex%mtln%probes(1)%index = 8
      ex%mtln%probes(1)%probe_type = PROBE_TYPE_CURRENT
      ex%mtln%probes(1)%probe_name = "test"
      ex%mtln%probes(1)%probe_position = [2,11,14]


      ! networks
      deallocate(ex%mtln%networks)
      allocate(ex%mtln%networks(2))

      allocate(ex%mtln%networks(1)%connections(2))
      allocate(ex%mtln%networks(1)%connections(1)%nodes(1))
      allocate(ex%mtln%networks(1)%connections(2)%nodes(1))
      ex%mtln%networks(1)%connections(1)%nodes(1)%conductor_in_cable = 1
      ex%mtln%networks(1)%connections(1)%nodes(1)%side = TERMINAL_NODE_SIDE_INI
      ex%mtln%networks(1)%connections(1)%nodes(1)%belongs_to_cable =>  ex%mtln%cables(1)%ptr
      ex%mtln%networks(1)%connections(1)%nodes(1)%termination%termination_type = TERMINATION_OPEN
      ex%mtln%networks(1)%connections(2)%nodes(1)%conductor_in_cable = 2
      ex%mtln%networks(1)%connections(2)%nodes(1)%side = TERMINAL_NODE_SIDE_INI
      ex%mtln%networks(1)%connections(2)%nodes(1)%belongs_to_cable =>  ex%mtln%cables(1)%ptr
      ex%mtln%networks(1)%connections(2)%nodes(1)%termination%termination_type = TERMINATION_OPEN

      allocate(ex%mtln%networks(2)%connections(2))
      allocate(ex%mtln%networks(2)%connections(1)%nodes(1))
      allocate(ex%mtln%networks(2)%connections(2)%nodes(1))
      ex%mtln%networks(2)%connections(1)%nodes(1)%conductor_in_cable = 1
      ex%mtln%networks(2)%connections(1)%nodes(1)%side = TERMINAL_NODE_SIDE_END
      ex%mtln%networks(2)%connections(1)%nodes(1)%belongs_to_cable =>  ex%mtln%cables(1)%ptr
      ex%mtln%networks(2)%connections(1)%nodes(1)%termination%termination_type = TERMINATION_OPEN
      ex%mtln%networks(2)%connections(2)%nodes(1)%conductor_in_cable = 2
      ex%mtln%networks(2)%connections(2)%nodes(1)%side = TERMINAL_NODE_SIDE_END
      ex%mtln%networks(2)%connections(2)%nodes(1)%belongs_to_cable =>  ex%mtln%cables(1)%ptr
      ex%mtln%networks(2)%connections(2)%nodes(1)%termination%termination_type = TERMINATION_OPEN


   end function
end function

