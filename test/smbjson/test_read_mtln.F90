integer function test_read_mtln() bind (C) result(err)
   use smbjson
   use smbjson_testingTools

   implicit none

   character(len=*),parameter :: filename = PATH_TO_TEST_DATA//INPUT_EXAMPLES//'mtln.fdtd.json'
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
      integer :: i, j

      call initializeProblemDescription(expected)

      ! Expected general info.
      expected%general%dt = 1e-12
      expected%general%nmax = 1000

      ! Excected media matrix.
      expected%matriz%totalX = 100
      expected%matriz%totalY = 7
      expected%matriz%totalZ = 2

      ! Expected grid.
      expected%despl%nX = 100
      expected%despl%nY = 7
      expected%despl%nZ = 2

      allocate(expected%despl%desX(100))
      allocate(expected%despl%desY(7))
      allocate(expected%despl%desZ(2))
      expected%despl%desX = 0.1
      expected%despl%desY = 0.1
      expected%despl%desZ = 0.1
      expected%despl%mx1 = 0
      expected%despl%mx2 = 100
      expected%despl%my1 = 0
      expected%despl%my2 = 7
      expected%despl%mz1 = 0
      expected%despl%mz2 = 2

      ! Expected boundaries.
      expected%front%tipoFrontera(:) = F_MUR

      ! Expected sources.
      allocate(expected%nodSrc%NodalSource(1))
      expected%nodSrc%n_C1P_max = 0
      expected%nodSrc%n_C2P_max = 1
      expected%nodSrc%n_nodSrc = 1
      expected%nodSrc%n_nodSrc_max = 1
      expected%nodSrc%NodalSource(1)%nombre = trim(adjustl("gauss.exc"))
      expected%nodSrc%NodalSource(1)%isElec = .true.
      expected%nodSrc%NodalSource(1)%isHard = .false.
      expected%nodSrc%NodalSource(1)%isInitialValue = .false.
      expected%nodSrc%NodalSource(1)%n_C1P = 0
      allocate(expected%nodSrc%NodalSource(1)%c1P(0))
      expected%nodSrc%NodalSource(1)%n_C2P = 1
      allocate(expected%nodSrc%NodalSource(1)%c2P(1))
      expected%nodSrc%NodalSource(1)%c2P(1)%xi = 2
      expected%nodSrc%NodalSource(1)%c2P(1)%xe = 3
      expected%nodSrc%NodalSource(1)%c2P(1)%yi = 7
      expected%nodSrc%NodalSource(1)%c2P(1)%ye = 7
      expected%nodSrc%NodalSource(1)%c2P(1)%zi = 1
      expected%nodSrc%NodalSource(1)%c2P(1)%ze = 1
      expected%nodSrc%NodalSource(1)%c2P(1)%xc = 1
      expected%nodSrc%NodalSource(1)%c2P(1)%yc = 0
      expected%nodSrc%NodalSource(1)%c2P(1)%zc = 0
      expected%nodSrc%NodalSource(1)%c2P(1)%or = 1
      expected%nodSrc%NodalSource(1)%c2P(1)%tag = "DistributedSource"
      
            

      ! Expected mtln type
      !connectors
      ! id = 24
      expected%mtln%has_multiwires = .true.
      expected%mtln%time_step = 1e-12
      expected%mtln%number_of_steps = 1000

      allocate(expected%mtln%connectors(4))

      allocate(expected%mtln%connectors(1)%resistances(1))
      expected%mtln%connectors(1)%id = 24
      expected%mtln%connectors(1)%resistances = [100.0e-3]
      expected%mtln%connectors(1)%transfer_impedance_per_meter%direction = TRANSFER_IMPEDANCE_DIRECTION_INWARDS
      allocate(expected%mtln%connectors(1)%transfer_impedance_per_meter%poles(0))
      allocate(expected%mtln%connectors(1)%transfer_impedance_per_meter%residues(0))

      ! id = 25
      allocate(expected%mtln%connectors(2)%resistances(1))
      expected%mtln%connectors(2)%id = 25
      expected%mtln%connectors(2)%resistances = [19.0]
      expected%mtln%connectors(2)%transfer_impedance_per_meter%direction = TRANSFER_IMPEDANCE_DIRECTION_INWARDS
      allocate(expected%mtln%connectors(2)%transfer_impedance_per_meter%poles(0))
      allocate(expected%mtln%connectors(2)%transfer_impedance_per_meter%residues(0))

      ! id = 204
      allocate(expected%mtln%connectors(3)%resistances(1))
      expected%mtln%connectors(3)%id = 204
      expected%mtln%connectors(3)%resistances = [100.0e-3]
      expected%mtln%connectors(3)%transfer_impedance_per_meter%direction = TRANSFER_IMPEDANCE_DIRECTION_INWARDS
      expected%mtln%connectors(3)%transfer_impedance_per_meter%resistive_term = 3.33
      expected%mtln%connectors(3)%transfer_impedance_per_meter%inductive_term = 2.6e-9
      allocate(expected%mtln%connectors(3)%transfer_impedance_per_meter%poles(0))
      allocate(expected%mtln%connectors(3)%transfer_impedance_per_meter%residues(0))

      ! id = 205
      allocate(expected%mtln%connectors(4)%resistances(1))
      expected%mtln%connectors(4)%id = 205
      expected%mtln%connectors(4)%resistances = [19.0]
      expected%mtln%connectors(4)%transfer_impedance_per_meter%direction = TRANSFER_IMPEDANCE_DIRECTION_INWARDS
      expected%mtln%connectors(4)%transfer_impedance_per_meter%resistive_term = 609.3
      expected%mtln%connectors(4)%transfer_impedance_per_meter%inductive_term = 2.6e-9
      allocate(expected%mtln%connectors(4)%transfer_impedance_per_meter%poles(0))
      allocate(expected%mtln%connectors(4)%transfer_impedance_per_meter%residues(0))

      !cables
      deallocate(expected%mtln%cables)
      allocate(expected%mtln%cables(9))
      ! level 0
      ! cable 1 - wire
      allocate(unshielded_multiwire_t :: expected%mtln%cables(1)%ptr)
      ptr => expected%mtln%cables(1)%ptr
      call initializeCablePULParameters(ptr)
      select type(ptr)
      type is (unshielded_multiwire_t)
         ptr%name = "line_0_0"
         ptr%cell_inductance_per_meter = reshape( source = [5.481553487168089e-07], shape = [ 1,1 ] )
         ptr%cell_capacitance_per_meter = reshape( source = [2.0270004E-11], shape = [ 1,1 ] )
         ptr%resistance_per_meter =  reshape(source=[22.9e-3], shape=[1,1])
         
         deallocate(ptr%multipolar_expansion)
         allocate(ptr%multipolar_expansion(0))

         allocate(ptr%step_size(9))
         ptr%step_size = [(0.1, i = 1, 9)]
         allocate(ptr%segments(9))
         do i = 1, 9
            ptr%segments(i)%x = i
            ptr%segments(i)%y = 7
            ptr%segments(i)%z = 1
            ptr%segments(i)%orientation = DIRECTION_X_POS
         end do

         ptr%initial_connector => expected%mtln%connectors(1)
         ptr%end_connector => null()
      end select
      ! cable 2 - shieldedMultiwire
      allocate(shielded_multiwire_t :: expected%mtln%cables(2)%ptr)
      ptr => expected%mtln%cables(2)%ptr
      call initializeCablePULParameters(ptr)
      select type(ptr)
      type is (shielded_multiwire_t)
         ptr%name = "line_1_0"
         ptr%inductance_per_meter = reshape(source=[8.802075200000001e-08], shape=[1,1])
         ptr%capacitance_per_meter = reshape(source=[5.5840010E-10], shape=[1,1])
         ptr%resistance_per_meter = reshape(source=[3.9e-3], shape=[1,1])

         allocate(ptr%step_size(9))
         ptr%step_size =  [(0.1, i = 1, 9)]

         allocate(ptr%segments(9))
         do i = 1, 9
            ptr%segments(i)%x = i
            ptr%segments(i)%y = 7
            ptr%segments(i)%z = 1
            ptr%segments(i)%orientation = DIRECTION_X_POS
         end do

         ptr%transfer_impedance%direction = TRANSFER_IMPEDANCE_DIRECTION_INWARDS
         ptr%transfer_impedance%resistive_term = 0.0
         ptr%transfer_impedance%inductive_term = 8.9e-9
         allocate(ptr%transfer_impedance%poles(0))
         allocate(ptr%transfer_impedance%residues(0))

         ptr%parent_cable => expected%mtln%cables(1)%ptr
         ptr%conductor_in_parent = 1
         ptr%initial_connector => expected%mtln%connectors(3)
         ptr%end_connector => null()
      end select
      ! cable 3 - shieldedMultiwire
      allocate(shielded_multiwire_t :: expected%mtln%cables(3)%ptr)
      ptr => expected%mtln%cables(3)%ptr
      call initializeCablePULParameters(ptr,8)
      select type(ptr)
      type is (shielded_multiwire_t)
         ptr%name = "line_2_0"
         ptr%inductance_per_meter(1:2,1:2) = & 
            reshape(source=[2.4382084E-07, 4.7377505E-08, 4.7377508E-08, 2.4382081E-07], shape=[2,2], order =[2,1])
         ptr%inductance_per_meter(3:4,3:4) = & 
            reshape(source=[2.4382084E-07, 4.7377505E-08, 4.7377508E-08, 2.4382081E-07], shape=[2,2], order =[2,1])
         ptr%inductance_per_meter(5:6,5:6) = & 
            reshape(source=[2.4382084E-07, 4.7377505E-08, 4.7377508E-08, 2.4382081E-07], shape=[2,2], order =[2,1])
         ptr%inductance_per_meter(7:8,7:8) = & 
            reshape(source=[2.4382084E-07, 4.7377505E-08, 4.7377508E-08, 2.4382081E-07], shape=[2,2], order =[2,1])

         ptr%capacitance_per_meter(1:2,1:2) = &
            reshape(source=[105.5e-12, -20.5e-12, -20.5e-12, 105.5e-12], shape=[2,2], order =[2,1])
         ptr%capacitance_per_meter(3:4,3:4) = &
            reshape(source=[105.5e-12, -20.5e-12, -20.5e-12, 105.5e-12], shape=[2,2], order =[2,1])
         ptr%capacitance_per_meter(5:6,5:6) = &
            reshape(source=[105.5e-12, -20.5e-12, -20.5e-12, 105.5e-12], shape=[2,2], order =[2,1])
         ptr%capacitance_per_meter(7:8,7:8) = &
            reshape(source=[105.5e-12, -20.5e-12, -20.5e-12, 105.5e-12], shape=[2,2], order =[2,1])


         do i = 1, 8
            ptr%resistance_per_meter(i,i) =  62.0e-3
         end do

         allocate(ptr%step_size(9))
         ptr%step_size =  [(0.1, i = 1, 9)]      

         allocate(ptr%segments(9))
         do i = 1, 9
            ptr%segments(i)%x = i
            ptr%segments(i)%y = 7
            ptr%segments(i)%z = 1
            ptr%segments(i)%orientation = DIRECTION_X_POS

         end do


         ptr%transfer_impedance%direction = TRANSFER_IMPEDANCE_DIRECTION_INWARDS
         ptr%transfer_impedance%resistive_term = 0.0
         ptr%transfer_impedance%inductive_term = 4.2e-9
         allocate(ptr%transfer_impedance%poles(0))
         allocate(ptr%transfer_impedance%residues(0))

         ptr%parent_cable => expected%mtln%cables(2)%ptr
         ptr%conductor_in_parent = 1
         ptr%initial_connector => null()
         ptr%end_connector => null()
      end select

      ! cable 4 - wire
      allocate(unshielded_multiwire_t :: expected%mtln%cables(4)%ptr)
      ptr => expected%mtln%cables(4)%ptr
      call initializeCablePULParameters(ptr)
      select type(ptr)
      type is (unshielded_multiwire_t)
         ptr%name = "line_0_1"

         ptr%cell_inductance_per_meter = reshape( source = [6.482560773828984e-07], shape = [ 1,1 ] )
         ptr%cell_capacitance_per_meter = reshape( source = [1.7140003E-11], shape = [ 1,1 ] )
         ptr%resistance_per_meter =  reshape(source=[11.8e-3], shape=[1,1])
         
         deallocate(ptr%multipolar_expansion)
         allocate(ptr%multipolar_expansion(0))

         allocate(ptr%step_size(8))
         ptr%step_size = [(0.1, i = 1, 8)]

         allocate(ptr%segments(8))
         do i = 1, 8
            ptr%segments(i)%x = 9+i
            ptr%segments(i)%y = 7
            ptr%segments(i)%z = 1
            ptr%segments(i)%orientation = DIRECTION_X_POS
         end do
         ptr%initial_connector => expected%mtln%connectors(2)
         ptr%end_connector => null()
      end select
      ! cable 5 - shieldedMultiwire
      allocate(shielded_multiwire_t :: expected%mtln%cables(5)%ptr)
      ptr => expected%mtln%cables(5)%ptr
      call initializeCablePULParameters(ptr)
      select type(ptr)
      type is (shielded_multiwire_t)
         ptr%name = "line_1_1"
         ptr%inductance_per_meter = reshape(source=[1.37228e-07], shape=[1,1])
         ptr%capacitance_per_meter = reshape(source=[3.2310005E-10], shape=[1,1])
         ptr%resistance_per_meter = reshape(source=[12.2e-3], shape=[1,1])
         ptr%conductance_per_meter = reshape(source=[0.0], shape=[1,1])

         allocate(ptr%step_size(8))
         ptr%step_size =  [(0.1, i = 1, 8)]

         allocate(ptr%segments(8))
         do i = 1, 8
            ptr%segments(i)%x = 9+i
            ptr%segments(i)%y = 7
            ptr%segments(i)%z = 1
            ptr%segments(i)%orientation = DIRECTION_X_POS
         end do

         ptr%transfer_impedance%direction = TRANSFER_IMPEDANCE_DIRECTION_INWARDS
         ptr%transfer_impedance%resistive_term = 0.0
         ptr%transfer_impedance%inductive_term = 7.4e-9
         allocate(ptr%transfer_impedance%poles(0))
         allocate(ptr%transfer_impedance%residues(0))

         ptr%parent_cable => expected%mtln%cables(4)%ptr
         ptr%conductor_in_parent = 1
         ptr%initial_connector => expected%mtln%connectors(4)
         ptr%end_connector => null()
      end select
      ! cable 6 - shieldedMultiwire
      allocate(shielded_multiwire_t :: expected%mtln%cables(6)%ptr)
      ptr => expected%mtln%cables(6)%ptr
      call initializeCablePULParameters(ptr,2)
      select type(ptr)
      type is (shielded_multiwire_t)
         ptr%name = "line_2_4"
         ptr%inductance_per_meter(1:2,1:2) = & 
            reshape(source=[2.4382084E-07, 4.7377505E-08, 4.7377508E-08, 2.4382081E-07], shape=[2,2], order =[2,1])

         ptr%capacitance_per_meter(1:2,1:2) = &
            reshape(source=[105.5e-12, -20.5e-12, -20.5e-12, 105.5e-12], shape=[2,2], order =[2,1])

         do i = 1, 2
            ptr%resistance_per_meter(i,i) = 62.0e-3
         end do

         allocate(ptr%step_size(8))
         ptr%step_size =  [(0.1, i = 1, 8)]

         allocate(ptr%segments(8))
         do i = 1, 8
            ptr%segments(i)%x = 9+i
            ptr%segments(i)%y = 7
            ptr%segments(i)%z = 1
            ptr%segments(i)%orientation = DIRECTION_X_POS
         end do


         ptr%transfer_impedance%direction = TRANSFER_IMPEDANCE_DIRECTION_INWARDS
         ptr%transfer_impedance%resistive_term = 0.0
         ptr%transfer_impedance%inductive_term = 4.2e-9
         allocate(ptr%transfer_impedance%poles(0))
         allocate(ptr%transfer_impedance%residues(0))

         ptr%parent_cable => expected%mtln%cables(5)%ptr
         ptr%conductor_in_parent = 1
         ptr%initial_connector => null()
         ptr%end_connector => null()
      end select
      ! cable 7 - wire
      allocate(unshielded_multiwire_t :: expected%mtln%cables(7)%ptr)
      ptr => expected%mtln%cables(7)%ptr
      call initializeCablePULParameters(ptr)
      select type(ptr)
      type is (unshielded_multiwire_t)
         ptr%name = "line_0_2"
         ptr%cell_inductance_per_meter = reshape( source = [5.802145885361537e-07], shape = [ 1,1 ] )
         ptr%cell_capacitance_per_meter = reshape( source = [1.9150003E-11], shape = [ 1,1 ] )
         ptr%resistance_per_meter =  reshape(source=[17.3e-3], shape=[1,1])
         
         deallocate(ptr%multipolar_expansion)
         allocate(ptr%multipolar_expansion(0))

         allocate(ptr%step_size(7))
         ptr%step_size = [(0.1, i = 1, 7)]

         allocate(ptr%segments(7))
         do i = 1,7
            ptr%segments(i)%x = 10
            ptr%segments(i)%y = 7-i
            ptr%segments(i)%z = 1
            ptr%segments(i)%orientation = DIRECTION_Y_NEG
         end do

         ptr%initial_connector => null()
         ptr%end_connector => null()
      end select
      ! cable 8 - shieldedMultiwire
      allocate(shielded_multiwire_t :: expected%mtln%cables(8)%ptr)
      ptr => expected%mtln%cables(8)%ptr
      call initializeCablePULParameters(ptr)
      select type(ptr)
      type is (shielded_multiwire_t)
         ptr%name = "line_1_2"
         ptr%inductance_per_meter = reshape(source=[9.1890502e-08], shape=[1,1])
         ptr%capacitance_per_meter = reshape(source=[4.7190007E-10], shape=[1,1])
         ptr%resistance_per_meter = reshape(source=[6.5e-3], shape=[1,1])
         ptr%conductance_per_meter = reshape(source=[0.0], shape=[1,1])

         allocate(ptr%step_size(7))
         ptr%step_size =  [(0.1, i = 1, 7)]

         allocate(ptr%segments(7))
         do i = 1,7
            ptr%segments(i)%x = 10
            ptr%segments(i)%y = 7-i
            ptr%segments(i)%z = 1
            ptr%segments(i)%orientation = DIRECTION_Y_NEG
         end do

         ptr%transfer_impedance%direction = TRANSFER_IMPEDANCE_DIRECTION_INWARDS
         ptr%transfer_impedance%resistive_term = 0.0
         ptr%transfer_impedance%inductive_term = 3.0e-9
         allocate(ptr%transfer_impedance%poles(0))
         allocate(ptr%transfer_impedance%residues(0))

         ptr%parent_cable => expected%mtln%cables(7)%ptr
         ptr%conductor_in_parent = 1
         ptr%initial_connector => null()
         ptr%end_connector => null()
      end select
      ! cable 9 - shieldedMultiwire
      allocate(shielded_multiwire_t :: expected%mtln%cables(9)%ptr)
      ptr => expected%mtln%cables(9)%ptr
      call initializeCablePULParameters(ptr,6)
      select type(ptr)
      type is (shielded_multiwire_t)
         ptr%name = "line_2_5"
         ptr%inductance_per_meter(1:2,1:2) = & 
            reshape(source=[2.4382084E-07, 4.7377505E-08, 4.7377508E-08, 2.4382081E-07], shape=[2,2], order =[2,1])
         ptr%inductance_per_meter(3:4,3:4) = & 
            reshape(source=[2.4382084E-07, 4.7377505E-08, 4.7377508E-08, 2.4382081E-07], shape=[2,2], order =[2,1])
         ptr%inductance_per_meter(5:6,5:6) = & 
            reshape(source=[2.4382084E-07, 4.7377505E-08, 4.7377508E-08, 2.4382081E-07], shape=[2,2], order =[2,1])

         ptr%capacitance_per_meter(1:2,1:2) = &
            reshape(source=[105.5e-12, -20.5e-12, -20.5e-12, 105.5e-12], shape=[2,2], order =[2,1])
         ptr%capacitance_per_meter(3:4,3:4) = &
            reshape(source=[105.5e-12, -20.5e-12, -20.5e-12, 105.5e-12], shape=[2,2], order =[2,1])
         ptr%capacitance_per_meter(5:6,5:6) = &
            reshape(source=[105.5e-12, -20.5e-12, -20.5e-12, 105.5e-12], shape=[2,2], order =[2,1])

         do i = 1, 6
            ptr%resistance_per_meter(i,i) = 62.0e-3
         end do
         allocate(ptr%step_size(7))
         ptr%step_size =  [(0.1, i = 1, 7)]

         allocate(ptr%segments(7))
         do i = 1,7
            ptr%segments(i)%x = 10
            ptr%segments(i)%y = 7-i
            ptr%segments(i)%z = 1
            ptr%segments(i)%orientation = DIRECTION_Y_NEG
         end do

         ptr%transfer_impedance%direction = TRANSFER_IMPEDANCE_DIRECTION_INWARDS
         ptr%transfer_impedance%resistive_term = 0.0
         ptr%transfer_impedance%inductive_term = 4.2e-9
         allocate(ptr%transfer_impedance%poles(0))
         allocate(ptr%transfer_impedance%residues(0))

         ptr%parent_cable => expected%mtln%cables(8)%ptr
         ptr%conductor_in_parent = 1
         ptr%initial_connector => null()
         ptr%end_connector => null()
      end select

      ! probes
      deallocate(expected%mtln%probes)
      allocate(expected%mtln%probes(7))
      expected%mtln%probes(1)%attached_to_cable => expected%mtln%cables(1)%ptr ! to which cable is the probe attached in mtln?
      expected%mtln%probes(1)%index = 1
      expected%mtln%probes(1)%probe_type = PROBE_TYPE_VOLTAGE
      expected%mtln%probes(1)%probe_name = "b1_terminal_voltage"
      expected%mtln%probes(1)%probe_position = [1,7,1]

      expected%mtln%probes(2)%attached_to_cable => expected%mtln%cables(1)%ptr
      expected%mtln%probes(2)%index = 1
      expected%mtln%probes(2)%probe_type = PROBE_TYPE_CURRENT
      expected%mtln%probes(2)%probe_name = "b1_terminal_current"
      expected%mtln%probes(2)%probe_position = [1,7,1]

      expected%mtln%probes(3)%attached_to_cable => expected%mtln%cables(1)%ptr
      expected%mtln%probes(3)%index = 10
      expected%mtln%probes(3)%probe_type = PROBE_TYPE_CURRENT
      expected%mtln%probes(3)%probe_name = "junction_current"
      expected%mtln%probes(3)%probe_position = [10, 7, 1]

      expected%mtln%probes(4)%attached_to_cable => expected%mtln%cables(4)%ptr
      expected%mtln%probes(4)%index = 1
      expected%mtln%probes(4)%probe_type = PROBE_TYPE_CURRENT
      expected%mtln%probes(4)%probe_name = "junction_current"
      expected%mtln%probes(4)%probe_position = [10, 7, 1]

      expected%mtln%probes(5)%attached_to_cable => expected%mtln%cables(7)%ptr
      expected%mtln%probes(5)%index = 1
      expected%mtln%probes(5)%probe_type = PROBE_TYPE_CURRENT
      expected%mtln%probes(5)%probe_name = "junction_current"
      expected%mtln%probes(5)%probe_position = [10, 7, 1]

      expected%mtln%probes(6)%attached_to_cable => expected%mtln%cables(4)%ptr
      expected%mtln%probes(6)%index = 9
      expected%mtln%probes(6)%probe_type = PROBE_TYPE_CURRENT
      expected%mtln%probes(6)%probe_name = "b2_terminal_current"
      expected%mtln%probes(6)%probe_position = [ 18, 7, 1]

      expected%mtln%probes(7)%attached_to_cable => expected%mtln%cables(7)%ptr
      expected%mtln%probes(7)%index = 8
      expected%mtln%probes(7)%probe_type = PROBE_TYPE_CURRENT
      expected%mtln%probes(7)%probe_name = "b3_terminal_current"
      expected%mtln%probes(7)%probe_position = [10, 0, 1]



      ! networks
      deallocate(expected%mtln%networks)
      allocate(expected%mtln%networks(4))

      ! NETWORK 1
      allocate(expected%mtln%networks(1)%connections(10))

      allocate(expected%mtln%networks(1)%connections(1)%nodes(1))
      expected%mtln%networks(1)%connections(1)%nodes(1)%conductor_in_cable = 1
      expected%mtln%networks(1)%connections(1)%nodes(1)%side = TERMINAL_NODE_SIDE_INI
      expected%mtln%networks(1)%connections(1)%nodes(1)%belongs_to_cable =>  expected%mtln%cables(1)%ptr
      expected%mtln%networks(1)%connections(1)%nodes(1)%termination%termination_type = TERMINATION_SERIES
      expected%mtln%networks(1)%connections(1)%nodes(1)%termination%resistance = 0.7e-3

      allocate(expected%mtln%networks(1)%connections(2)%nodes(1))
      expected%mtln%networks(1)%connections(2)%nodes(1)%conductor_in_cable = 1
      expected%mtln%networks(1)%connections(2)%nodes(1)%side = TERMINAL_NODE_SIDE_INI
      expected%mtln%networks(1)%connections(2)%nodes(1)%belongs_to_cable => expected%mtln%cables(2)%ptr
      expected%mtln%networks(1)%connections(2)%nodes(1)%termination%termination_type = TERMINATION_SERIES
      expected%mtln%networks(1)%connections(2)%nodes(1)%termination%resistance = 1e-6

      do i = 3, 10
         allocate(expected%mtln%networks(1)%connections(i)%nodes(1))
         expected%mtln%networks(1)%connections(i)%nodes(1)%side = TERMINAL_NODE_SIDE_INI
         expected%mtln%networks(1)%connections(i)%nodes(1)%belongs_to_cable => expected%mtln%cables(3)%ptr
         expected%mtln%networks(1)%connections(i)%nodes(1)%conductor_in_cable = i-2
      end do

      expected%mtln%networks(1)%connections(4)%nodes(1)%termination%termination_type = TERMINATION_SERIES
      expected%mtln%networks(1)%connections(4)%nodes(1)%termination%resistance = 1e10
      expected%mtln%networks(1)%connections(6)%nodes(1)%termination%termination_type = TERMINATION_SERIES
      expected%mtln%networks(1)%connections(6)%nodes(1)%termination%resistance = 1e10
      expected%mtln%networks(1)%connections(8)%nodes(1)%termination%termination_type = TERMINATION_SERIES
      expected%mtln%networks(1)%connections(8)%nodes(1)%termination%resistance = 1e10
      expected%mtln%networks(1)%connections(10)%nodes(1)%termination%termination_type = TERMINATION_SERIES
      expected%mtln%networks(1)%connections(10)%nodes(1)%termination%resistance = 1e10

      expected%mtln%networks(1)%connections(3)%nodes(1)%termination%termination_type = TERMINATION_RsLCp
      expected%mtln%networks(1)%connections(3)%nodes(1)%termination%resistance = 50.0
      expected%mtln%networks(1)%connections(3)%nodes(1)%termination%inductance = 30e-12
      expected%mtln%networks(1)%connections(3)%nodes(1)%termination%capacitance = 60e-9

      expected%mtln%networks(1)%connections(5)%nodes(1)%termination%termination_type = TERMINATION_RsLCp
      expected%mtln%networks(1)%connections(5)%nodes(1)%termination%resistance = 50.0
      expected%mtln%networks(1)%connections(5)%nodes(1)%termination%inductance = 30e-12
      expected%mtln%networks(1)%connections(5)%nodes(1)%termination%capacitance = 60e-9

      expected%mtln%networks(1)%connections(7)%nodes(1)%termination%termination_type = TERMINATION_RsLCp
      expected%mtln%networks(1)%connections(7)%nodes(1)%termination%resistance = 50.0
      expected%mtln%networks(1)%connections(7)%nodes(1)%termination%inductance = 30e-12
      expected%mtln%networks(1)%connections(7)%nodes(1)%termination%capacitance = 60e-9

      expected%mtln%networks(1)%connections(9)%nodes(1)%termination%termination_type = TERMINATION_RsLCp
      expected%mtln%networks(1)%connections(9)%nodes(1)%termination%resistance = 50.0
      expected%mtln%networks(1)%connections(9)%nodes(1)%termination%inductance = 30e-12
      expected%mtln%networks(1)%connections(9)%nodes(1)%termination%capacitance = 60e-9

      ! NETWORK 2
      allocate(expected%mtln%networks(2)%connections(10))
      allocate(expected%mtln%networks(2)%connections(1)%nodes(3))
      expected%mtln%networks(2)%connections(1)%nodes(1)%conductor_in_cable = 1
      expected%mtln%networks(2)%connections(1)%nodes(1)%side = TERMINAL_NODE_SIDE_END
      expected%mtln%networks(2)%connections(1)%nodes(1)%belongs_to_cable =>  expected%mtln%cables(1)%ptr

      expected%mtln%networks(2)%connections(1)%nodes(2)%conductor_in_cable = 1
      expected%mtln%networks(2)%connections(1)%nodes(2)%side = TERMINAL_NODE_SIDE_INI
      expected%mtln%networks(2)%connections(1)%nodes(2)%belongs_to_cable =>  expected%mtln%cables(4)%ptr

      expected%mtln%networks(2)%connections(1)%nodes(3)%conductor_in_cable = 1
      expected%mtln%networks(2)%connections(1)%nodes(3)%side = TERMINAL_NODE_SIDE_INI
      expected%mtln%networks(2)%connections(1)%nodes(3)%belongs_to_cable =>  expected%mtln%cables(7)%ptr

      do i = 1, 3
         expected%mtln%networks(2)%connections(1)%nodes(i)%termination%termination_type = TERMINATION_SERIES
         expected%mtln%networks(2)%connections(1)%nodes(i)%termination%resistance = 1e-6
      end do

      allocate(expected%mtln%networks(2)%connections(2)%nodes(3))
      expected%mtln%networks(2)%connections(2)%nodes(1)%conductor_in_cable = 1
      expected%mtln%networks(2)%connections(2)%nodes(1)%side = TERMINAL_NODE_SIDE_END
      expected%mtln%networks(2)%connections(2)%nodes(1)%belongs_to_cable =>  expected%mtln%cables(2)%ptr

      expected%mtln%networks(2)%connections(2)%nodes(2)%conductor_in_cable = 1
      expected%mtln%networks(2)%connections(2)%nodes(2)%side = TERMINAL_NODE_SIDE_INI
      expected%mtln%networks(2)%connections(2)%nodes(2)%belongs_to_cable =>  expected%mtln%cables(5)%ptr

      expected%mtln%networks(2)%connections(2)%nodes(3)%conductor_in_cable = 1
      expected%mtln%networks(2)%connections(2)%nodes(3)%side = TERMINAL_NODE_SIDE_INI
      expected%mtln%networks(2)%connections(2)%nodes(3)%belongs_to_cable =>  expected%mtln%cables(8)%ptr

      do i = 1, 3
         expected%mtln%networks(2)%connections(2)%nodes(i)%termination%termination_type = TERMINATION_SERIES
         expected%mtln%networks(2)%connections(2)%nodes(i)%termination%resistance = 1e-6
      end do

      ! NETWORK 2 - CONNECTIONS 3-10
      do i = 3, 8
         allocate(expected%mtln%networks(2)%connections(i)%nodes(2))

         expected%mtln%networks(2)%connections(i)%nodes(1)%conductor_in_cable = i-2
         expected%mtln%networks(2)%connections(i)%nodes(1)%side = TERMINAL_NODE_SIDE_END
         expected%mtln%networks(2)%connections(i)%nodes(1)%belongs_to_cable =>  expected%mtln%cables(3)%ptr
   
         expected%mtln%networks(2)%connections(i)%nodes(2)%conductor_in_cable = i-2
         expected%mtln%networks(2)%connections(i)%nodes(2)%side = TERMINAL_NODE_SIDE_INI
         expected%mtln%networks(2)%connections(i)%nodes(2)%belongs_to_cable =>  expected%mtln%cables(9)%ptr

         do j = 1, 2
            expected%mtln%networks(2)%connections(i)%nodes(j)%termination%termination_type = TERMINATION_SERIES
            expected%mtln%networks(2)%connections(i)%nodes(j)%termination%resistance = 1e-6
         end do

      end do

      do i = 9, 10
         allocate(expected%mtln%networks(2)%connections(i)%nodes(2))

         expected%mtln%networks(2)%connections(i)%nodes(1)%conductor_in_cable = i-2
         expected%mtln%networks(2)%connections(i)%nodes(1)%side = TERMINAL_NODE_SIDE_END
         expected%mtln%networks(2)%connections(i)%nodes(1)%belongs_to_cable =>  expected%mtln%cables(3)%ptr
   
         expected%mtln%networks(2)%connections(i)%nodes(2)%conductor_in_cable = i-8
         expected%mtln%networks(2)%connections(i)%nodes(2)%side = TERMINAL_NODE_SIDE_INI
         expected%mtln%networks(2)%connections(i)%nodes(2)%belongs_to_cable =>  expected%mtln%cables(6)%ptr

         do j = 1, 2
            expected%mtln%networks(2)%connections(i)%nodes(j)%termination%termination_type = TERMINATION_SERIES
            expected%mtln%networks(2)%connections(i)%nodes(j)%termination%resistance = 1e-6
         end do

      end do


      ! NETWORK 3
      allocate(expected%mtln%networks(3)%connections(4))

      allocate(expected%mtln%networks(3)%connections(1)%nodes(1))
      expected%mtln%networks(3)%connections(1)%nodes(1)%conductor_in_cable = 1
      expected%mtln%networks(3)%connections(1)%nodes(1)%side = TERMINAL_NODE_SIDE_END
      expected%mtln%networks(3)%connections(1)%nodes(1)%belongs_to_cable =>  expected%mtln%cables(4)%ptr
      expected%mtln%networks(3)%connections(1)%nodes(1)%termination%termination_type = TERMINATION_SERIES
      expected%mtln%networks(3)%connections(1)%nodes(1)%termination%resistance = 1

      allocate(expected%mtln%networks(3)%connections(2)%nodes(1))
      expected%mtln%networks(3)%connections(2)%nodes(1)%conductor_in_cable = 1
      expected%mtln%networks(3)%connections(2)%nodes(1)%side = TERMINAL_NODE_SIDE_END
      expected%mtln%networks(3)%connections(2)%nodes(1)%belongs_to_cable => expected%mtln%cables(5)%ptr
      expected%mtln%networks(3)%connections(2)%nodes(1)%termination%termination_type = TERMINATION_SERIES
      expected%mtln%networks(3)%connections(2)%nodes(1)%termination%resistance = 1e-6

      do i = 3, 4
         allocate(expected%mtln%networks(3)%connections(i)%nodes(1))
         expected%mtln%networks(3)%connections(i)%nodes(1)%side = TERMINAL_NODE_SIDE_END
         expected%mtln%networks(3)%connections(i)%nodes(1)%belongs_to_cable => expected%mtln%cables(6)%ptr
         expected%mtln%networks(3)%connections(i)%nodes(1)%conductor_in_cable = i-2
      end do

      expected%mtln%networks(3)%connections(3)%nodes(1)%termination%termination_type = TERMINATION_SERIES
      expected%mtln%networks(3)%connections(3)%nodes(1)%termination%resistance = 50
      expected%mtln%networks(3)%connections(4)%nodes(1)%termination%termination_type = TERMINATION_SERIES
      expected%mtln%networks(3)%connections(4)%nodes(1)%termination%resistance = 50


      ! NETWORK 4
      allocate(expected%mtln%networks(4)%connections(8))

      allocate(expected%mtln%networks(4)%connections(1)%nodes(1))
      expected%mtln%networks(4)%connections(1)%nodes(1)%conductor_in_cable = 1
      expected%mtln%networks(4)%connections(1)%nodes(1)%side = TERMINAL_NODE_SIDE_END
      expected%mtln%networks(4)%connections(1)%nodes(1)%belongs_to_cable =>  expected%mtln%cables(7)%ptr
      expected%mtln%networks(4)%connections(1)%nodes(1)%termination%termination_type = TERMINATION_SERIES
      expected%mtln%networks(4)%connections(1)%nodes(1)%termination%resistance = 0.7e-3

      allocate(expected%mtln%networks(4)%connections(2)%nodes(1))
      expected%mtln%networks(4)%connections(2)%nodes(1)%conductor_in_cable = 1
      expected%mtln%networks(4)%connections(2)%nodes(1)%side = TERMINAL_NODE_SIDE_END
      expected%mtln%networks(4)%connections(2)%nodes(1)%belongs_to_cable => expected%mtln%cables(8)%ptr
      expected%mtln%networks(4)%connections(2)%nodes(1)%termination%termination_type = TERMINATION_SERIES
      expected%mtln%networks(4)%connections(2)%nodes(1)%termination%resistance = 1e-6

      do i = 3, 8
         allocate(expected%mtln%networks(4)%connections(i)%nodes(1))
         expected%mtln%networks(4)%connections(i)%nodes(1)%side = TERMINAL_NODE_SIDE_END
         expected%mtln%networks(4)%connections(i)%nodes(1)%belongs_to_cable => expected%mtln%cables(9)%ptr
         expected%mtln%networks(4)%connections(i)%nodes(1)%conductor_in_cable = i-2
      end do

      expected%mtln%networks(4)%connections(3)%nodes(1)%termination%termination_type = TERMINATION_SERIES
      expected%mtln%networks(4)%connections(3)%nodes(1)%termination%resistance = 50
      expected%mtln%networks(4)%connections(4)%nodes(1)%termination%termination_type = TERMINATION_SERIES
      expected%mtln%networks(4)%connections(4)%nodes(1)%termination%resistance = 50
      expected%mtln%networks(4)%connections(5)%nodes(1)%termination%termination_type = TERMINATION_SERIES
      expected%mtln%networks(4)%connections(5)%nodes(1)%termination%resistance = 50
      expected%mtln%networks(4)%connections(6)%nodes(1)%termination%termination_type = TERMINATION_SERIES
      expected%mtln%networks(4)%connections(6)%nodes(1)%termination%resistance = 50
      expected%mtln%networks(4)%connections(7)%nodes(1)%termination%termination_type = TERMINATION_SERIES
      expected%mtln%networks(4)%connections(7)%nodes(1)%termination%resistance = 50
      expected%mtln%networks(4)%connections(8)%nodes(1)%termination%termination_type = TERMINATION_SERIES
      expected%mtln%networks(4)%connections(8)%nodes(1)%termination%resistance = 50

   end function
end function

