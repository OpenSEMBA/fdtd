integer function test_read_shieldedpair() bind (C) result(err)
   use smbjson
   use smbjson_testingTools

   implicit none
   character(len=*),parameter :: filename = PATH_TO_TEST_DATA//INPUT_EXAMPLES//'shieldedPair.fdtd.json'
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
      class(cable_t), pointer :: ptr
      integer :: i

      call initializeProblemDescription(expected)

      ! Expected general info.
      expected%general%dt = 0.43e-10
      expected%general%nmax = 700

      ! Excected media matrix.
      expected%matriz%totalX = 150
      expected%matriz%totalY = 150
      expected%matriz%totalZ = 150

      ! Expected grid.
      expected%despl%nX = 150
      expected%despl%nY = 150
      expected%despl%nZ = 150

      allocate(expected%despl%desX(150))
      allocate(expected%despl%desY(150))
      allocate(expected%despl%desZ(150))
      expected%despl%desX = 0.180
      expected%despl%desY = 0.180
      expected%despl%desZ = 0.0504
      expected%despl%mx1 = 0
      expected%despl%mx2 = 150
      expected%despl%my1 = 0
      expected%despl%my2 = 150
      expected%despl%mz1 = 0
      expected%despl%mz2 = 150

      ! Expected boundaries.
      expected%front%tipoFrontera(:) = F_PML
      expected%front%propiedadesPML(:)%numCapas = 6
      expected%front%propiedadesPML(:)%orden = 2.0
      expected%front%propiedadesPML(:)%refl = 0.0001

      ! Expected PEC regions.
      expected%pecRegs%nLins = 0
      expected%pecRegs%nLins_max = 0
      expected%pecRegs%nSurfs = 1
      expected%pecRegs%nSurfs_max = 1
      expected%pecRegs%nVols = 0
      expected%pecRegs%nVols_max = 0
      allocate(expected%pecRegs%Surfs(1))
      expected%pecRegs%Surfs(1)%Xi = 20
      expected%pecRegs%Surfs(1)%Xe = 129
      expected%pecRegs%Surfs(1)%yi = 20
      expected%pecRegs%Surfs(1)%ye = 129
      expected%pecRegs%Surfs(1)%zi = 74
      expected%pecRegs%Surfs(1)%ze = 74
      expected%pecRegs%Surfs(1)%Xtrancos = 1
      expected%pecRegs%Surfs(1)%Ytrancos = 1
      expected%pecRegs%Surfs(1)%Ztrancos = 1
      expected%pecRegs%Surfs(1)%Or = 3
      expected%pecRegs%Surfs(1)%tag =  trim(adjustl("material5@layer5"))
   
      ! Expected sources.
      allocate(expected%plnSrc%collection(1))
      expected%plnSrc%collection(1)%nombre_fichero = "shielded_pair.exc"
      expected%plnSrc%collection(1)%atributo = ""
      expected%plnSrc%collection(1)%coor1 = [10, 10, 10]
      expected%plnSrc%collection(1)%coor2 = [139, 139, 139]
      expected%plnSrc%collection(1)%theta = 3.1416
      expected%plnSrc%collection(1)%phi = 0.0
      expected%plnSrc%collection(1)%alpha = 1.5708
      expected%plnSrc%collection(1)%beta = -1.5708
      expected%plnSrc%collection(1)%isRC=.false.
      expected%plnSrc%collection(1)%nummodes=1
      expected%plnSrc%collection(1)%INCERTMAX=0.0
      expected%plnSrc%nc = 1
      expected%plnSrc%nC_max = 1

      
      ! Expected mtln type
      expected%mtln%has_multiwires = .true.
      expected%mtln%time_step = 0.43e-10
      expected%mtln%number_of_steps = 700

      
      deallocate(expected%mtln%cables)
      allocate(expected%mtln%cables(2))
      ! cable 1 - wire
      allocate(unshielded_multiwire_t :: expected%mtln%cables(1)%ptr)
      call initializeCablePULParameters(expected%mtln%cables(1)%ptr)
      ptr => expected%mtln%cables(1)%ptr
      select type(ptr)
      type is(unshielded_multiwire_t)
         expected%mtln%cables(1)%ptr%name = "line_0"
         ptr%resistance_per_meter = reshape(source=[22.9e-3], shape=[1,1])

         allocate(ptr%step_size(5))
         ptr%step_size(1) =  0.0504
         ptr%step_size(2:4) =  [(0.180, i = 2, 4)]
         ptr%step_size(5) =  0.0504

         allocate(ptr%segments(5))
         ptr%segments(1)%x = 75
         ptr%segments(1)%y = 71
         ptr%segments(1)%z = 74
         ptr%segments(1)%orientation = DIRECTION_Z_POS
         do i = 2, 4
            ptr%segments(i)%x = 75
            ptr%segments(i)%y = 69+i
            ptr%segments(i)%z = 75
            ptr%segments(i)%orientation = DIRECTION_Y_POS
         end do
         ptr%segments(5)%x = 75
         ptr%segments(5)%y = 74
         ptr%segments(5)%z = 74
         ptr%segments(5)%orientation = DIRECTION_Z_NEG

         ptr%initial_connector => null()
         ptr%end_connector => null()
      end select
      ! cable 2 - shieldedMultiwire
      allocate(shielded_multiwire_t :: expected%mtln%cables(2)%ptr)
      call initializeCablePULParameters(expected%mtln%cables(2)%ptr)
      ptr => expected%mtln%cables(2)%ptr
      select type(ptr)
      type is(shielded_multiwire_t)
         ptr%name = "line_1"

         ptr%inductance_per_meter = & 
            reshape( source = [ 3.13182309e-07, 7.45674981e-08, 7.45674981e-08, 3.13182309e-07 ], shape = [ 2,2 ] )
         ptr%capacitance_per_meter = &
            reshape( source = [85.0e-12, -20.5e-12, -20.5e-12, 85.0e-12 ], shape = [ 2,2 ] )
         
         allocate(ptr%step_size(5))
         ptr%step_size(1) =  0.0504
         ptr%step_size(2:4) =  [(0.180, i = 2, 4)]
         ptr%step_size(5) =  0.0504

         allocate(ptr%segments(5))
         ptr%segments(1)%x = 75
         ptr%segments(1)%y = 71
         ptr%segments(1)%z = 74
         ptr%segments(1)%orientation = DIRECTION_Z_POS
         do i = 2, 4
            ptr%segments(i)%x = 75
            ptr%segments(i)%y = 69+i
            ptr%segments(i)%z = 75
            ptr%segments(i)%orientation = DIRECTION_Y_POS
         end do
         ptr%segments(5)%x = 75
         ptr%segments(5)%y = 74
         ptr%segments(5)%z = 74
         ptr%segments(5)%orientation = DIRECTION_Z_NEG

         ptr%transfer_impedance%direction = TRANSFER_IMPEDANCE_DIRECTION_BOTH
         ptr%transfer_impedance%resistive_term = 0.0
         ptr%transfer_impedance%inductive_term = 4.0e-9
         allocate(ptr%transfer_impedance%poles(0))
         allocate(ptr%transfer_impedance%residues(0))

         ptr%parent_cable => expected%mtln%cables(1)%ptr
         ptr%conductor_in_parent = 1
         ptr%initial_connector => null()
         ptr%end_connector => null()
      end select
      ! probes
      deallocate(expected%mtln%probes)
      allocate(expected%mtln%probes(4))
      expected%mtln%probes(1)%attached_to_cable => expected%mtln%cables(1)%ptr ! to which cable is the probe attached in mtln?
      expected%mtln%probes(1)%index = 1
      expected%mtln%probes(1)%probe_type = PROBE_TYPE_CURRENT
      expected%mtln%probes(1)%probe_name = "wire_end"
      expected%mtln%probes(1)%probe_position = [75,71,74]
      
      expected%mtln%probes(2)%attached_to_cable => expected%mtln%cables(1)%ptr
      expected%mtln%probes(2)%index = 1
      expected%mtln%probes(2)%probe_type = PROBE_TYPE_VOLTAGE
      expected%mtln%probes(2)%probe_name = "wire_end"
      expected%mtln%probes(2)%probe_position = [75,71,74]
      
      expected%mtln%probes(3)%attached_to_cable => expected%mtln%cables(1)%ptr ! to which cable is the probe attached in mtln?
      expected%mtln%probes(3)%index = 6
      expected%mtln%probes(3)%probe_type = PROBE_TYPE_CURRENT
      expected%mtln%probes(3)%probe_name = "wire_start"
      expected%mtln%probes(3)%probe_position = [75,74,74]
      
      expected%mtln%probes(4)%attached_to_cable => expected%mtln%cables(1)%ptr
      expected%mtln%probes(4)%index = 6
      expected%mtln%probes(4)%probe_type = PROBE_TYPE_VOLTAGE
      expected%mtln%probes(4)%probe_name = "wire_start"
      expected%mtln%probes(4)%probe_position = [75,74,74]

      ! networks
      deallocate(expected%mtln%networks)
      allocate(expected%mtln%networks(2))
      allocate(expected%mtln%networks(1)%connections(3))

      allocate(expected%mtln%networks(1)%connections(1)%nodes(1))
      expected%mtln%networks(1)%connections(1)%nodes(1)%conductor_in_cable = 1
      expected%mtln%networks(1)%connections(1)%nodes(1)%side = TERMINAL_NODE_SIDE_INI
      expected%mtln%networks(1)%connections(1)%nodes(1)%belongs_to_cable =>  expected%mtln%cables(1)%ptr
      expected%mtln%networks(1)%connections(1)%nodes(1)%termination%termination_type = TERMINATION_SERIES
      expected%mtln%networks(1)%connections(1)%nodes(1)%termination%resistance = 50.0

      allocate(expected%mtln%networks(1)%connections(2)%nodes(1))
      expected%mtln%networks(1)%connections(2)%nodes(1)%conductor_in_cable = 1
      expected%mtln%networks(1)%connections(2)%nodes(1)%side = TERMINAL_NODE_SIDE_INI
      expected%mtln%networks(1)%connections(2)%nodes(1)%belongs_to_cable => expected%mtln%cables(2)%ptr
      expected%mtln%networks(1)%connections(2)%nodes(1)%termination%termination_type = TERMINATION_SERIES
      expected%mtln%networks(1)%connections(2)%nodes(1)%termination%resistance = 50.0

      allocate(expected%mtln%networks(1)%connections(3)%nodes(1))
      expected%mtln%networks(1)%connections(3)%nodes(1)%conductor_in_cable = 2
      expected%mtln%networks(1)%connections(3)%nodes(1)%side = TERMINAL_NODE_SIDE_INI
      expected%mtln%networks(1)%connections(3)%nodes(1)%belongs_to_cable => expected%mtln%cables(2)%ptr
      expected%mtln%networks(1)%connections(3)%nodes(1)%termination%termination_type = TERMINATION_SERIES
      expected%mtln%networks(1)%connections(3)%nodes(1)%termination%resistance = 50.0

      allocate(expected%mtln%networks(2)%connections(3))

      allocate(expected%mtln%networks(2)%connections(1)%nodes(1))
      expected%mtln%networks(2)%connections(1)%nodes(1)%conductor_in_cable = 1
      expected%mtln%networks(2)%connections(1)%nodes(1)%side = TERMINAL_NODE_SIDE_END
      expected%mtln%networks(2)%connections(1)%nodes(1)%belongs_to_cable => expected%mtln%cables(1)%ptr
      expected%mtln%networks(2)%connections(1)%nodes(1)%termination%termination_type = TERMINATION_SERIES
      expected%mtln%networks(2)%connections(1)%nodes(1)%termination%resistance = 50.0

      allocate(expected%mtln%networks(2)%connections(2)%nodes(1))
      expected%mtln%networks(2)%connections(2)%nodes(1)%conductor_in_cable = 1
      expected%mtln%networks(2)%connections(2)%nodes(1)%side = TERMINAL_NODE_SIDE_END
      expected%mtln%networks(2)%connections(2)%nodes(1)%belongs_to_cable => expected%mtln%cables(2)%ptr
      expected%mtln%networks(2)%connections(2)%nodes(1)%termination%termination_type = TERMINATION_SERIES
      expected%mtln%networks(2)%connections(2)%nodes(1)%termination%resistance = 50.0

      allocate(expected%mtln%networks(2)%connections(3)%nodes(1))
      expected%mtln%networks(2)%connections(3)%nodes(1)%conductor_in_cable = 2
      expected%mtln%networks(2)%connections(3)%nodes(1)%side = TERMINAL_NODE_SIDE_END
      expected%mtln%networks(2)%connections(3)%nodes(1)%belongs_to_cable => expected%mtln%cables(2)%ptr
      expected%mtln%networks(2)%connections(3)%nodes(1)%termination%termination_type = TERMINATION_SERIES
      expected%mtln%networks(2)%connections(3)%nodes(1)%termination%resistance = 50.0

      !connectors
      allocate(expected%mtln%connectors(0))

   end function
end function

