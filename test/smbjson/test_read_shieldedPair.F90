integer function test_read_shieldedpair() bind (C) result(err)
   use smbjson
   use smbjson_testingTools

   implicit none

   character(len=*),parameter :: filename = PATH_TO_TEST_DATA//'cases/shieldedPair.fdtd.json'
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

      ! Expected sources.
      allocate(expected%plnSrc%collection(1))
      expected%plnSrc%collection(1)%nombre_fichero = "double_exp.exc"
      expected%plnSrc%collection(1)%atributo = ""
      expected%plnSrc%collection(1)%coor1 = [1, 1, 1]
      expected%plnSrc%collection(1)%coor2 = [2, 2, 18]
      expected%plnSrc%collection(1)%theta = 3.1416
      expected%plnSrc%collection(1)%phi = 1.5708
      expected%plnSrc%collection(1)%alpha = 1.5708
      expected%plnSrc%collection(1)%beta = 0.0
      expected%plnSrc%collection(1)%isRC=.false.
      expected%plnSrc%collection(1)%nummodes=1
      expected%plnSrc%collection(1)%INCERTMAX=0.0
      expected%plnSrc%nc = 1
      expected%plnSrc%nC_max = 1

      ! Expected probes
      ! oldSonda
      expected%oldSONDA%n_probes = 0
      expected%oldSONDA%n_probes_max = 0
      allocate(expected%oldSONDA%probes(0))

      ! sonda
      expected%Sonda%length = 2
      expected%Sonda%length_max = 2
      allocate(expected%Sonda%collection(2))
      
      expected%Sonda%collection(1)%outputrequest = "voltage_at_wire_end"
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
      expected%Sonda%collection(1)%cordinates(1)%tag = "voltage_at_wire_end"
      expected%Sonda%collection(1)%cordinates(1)%Xi = 2 ! Coord id as tag.
      expected%Sonda%collection(1)%cordinates(1)%Yi = 0
      expected%Sonda%collection(1)%cordinates(1)%Zi = 0
      expected%Sonda%collection(1)%cordinates(1)%Or = NP_COR_DDP
      
      expected%Sonda%collection(2)%outputrequest = "current_at_wire_end"
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
      expected%Sonda%collection(2)%cordinates(1)%tag = "current_at_wire_end"
      expected%Sonda%collection(2)%cordinates(1)%Xi = 2 ! Coord id as tag.
      expected%Sonda%collection(2)%cordinates(1)%Yi = 0
      expected%Sonda%collection(2)%cordinates(1)%Zi = 0
      expected%Sonda%collection(2)%cordinates(1)%Or = NP_COR_WIRECURRENT
            
      ! Expected thin wires
      allocate(expected%tWires%tw(1))
      expected%tWires%tw(1)%rad=0.0001
      expected%tWires%tw(1)%res=22.9e-3
      expected%tWires%tw(1)%dispfile = trim(adjustl(" "))
      expected%tWires%tw(1)%dispfile_LeftEnd = trim(adjustl(" "))
      expected%tWires%tw(1)%dispfile_RightEnd = trim(adjustl(" "))
      expected%tWires%tw(1)%n_twc=18
      expected%tWires%tw(1)%n_twc_max=18
      allocate(expected%tWires%tw(1)%twc(18))
      expected%tWires%tw(1)%twc(1:18)%srcfile = 'None'
      expected%tWires%tw(1)%twc(1:18)%srctype = 'None'
      expected%tWires%tw(1)%twc(1:18)%i = 2
      expected%tWires%tw(1)%twc(1:18)%j = 2
      expected%tWires%tw(1)%twc(1:18)%k = [(i, i=1, 18)]
      expected%tWires%tw(1)%twc(1:18)%d = DIR_Z
      expected%tWires%tw(1)%twc(1)%nd  = 1
      expected%tWires%tw(1)%twc(2:17)%nd = NO_TAG
      expected%tWires%tw(1)%twc(18)%nd  = 2
      
      expected%tWires%tw(1)%twc(1:18)%tag = trim(adjustl("1"))   ! The polyline id is used as tag.
      
      expected%tWires%tw(1)%tl = SERIES_CONS
      expected%tWires%tw(1)%R_LeftEnd = 50
      expected%tWires%tw(1)%C_LeftEnd = 1e22
      expected%tWires%tw(1)%tr = SERIES_CONS
      expected%tWires%tw(1)%R_RightEnd = 50
      expected%tWires%tw(1)%C_RightEnd = 1e22
            
      expected%tWires%n_tw = 1
      expected%tWires%n_tw_max = 1

      ! Expected mtln type
      allocate(expected%mtln%cables(2))
      ! cable 2 - wire
      expected%mtln%cables(2)%name = "line_0"
      allocate(expected%mtln%cables(2)%inductance_per_meter(1,1))
      allocate(expected%mtln%cables(2)%capacitance_per_meter(1,1))
      allocate(expected%mtln%cables(2)%resistance_per_meter(1,1))
      allocate(expected%mtln%cables(2)%conductance_per_meter(1,1))
      expected%mtln%cables(2)%inductance_per_meter = reshape(source=[5.362505362505362e-07], shape=[1,1])
      expected%mtln%cables(2)%capacitance_per_meter = reshape(source=[20.72e-12], shape=[1,1])
      expected%mtln%cables(2)%resistance_per_meter = reshape(source=[22.9e-3], shape=[1,1])
      expected%mtln%cables(2)%conductance_per_meter = reshape(source=[0.0], shape=[1,1])
      allocate(expected%mtln%cables(2)%step_size(18))
      expected%mtln%cables(2)%step_size =  [(0.03, i = 1, 18)]
      allocate(expected%mtln%cables(2)%external_field_segments(18))
      do i = 1, 18
         expected%mtln%cables(2)%external_field_segments(i)%position = (/2,2,i/)
         expected%mtln%cables(2)%external_field_segments(i)%direction = DIRECTION_Z_POS
         expected%mtln%cables(2)%external_field_segments(i)%field => null()

      end do

      allocate(expected%mtln%cables(2)%transfer_impedance%poles(0))
      allocate(expected%mtln%cables(2)%transfer_impedance%residues(0))

      expected%mtln%cables(2)%parent_cable => null()
      expected%mtln%cables(2)%conductor_in_parent = 0
      expected%mtln%cables(2)%initial_connector => null()
      expected%mtln%cables(2)%end_connector => null()

      ! cable 1 - multiwire
      expected%mtln%cables(1)%name = "line_1"
      allocate(expected%mtln%cables(1)%inductance_per_meter(2,2))
      allocate(expected%mtln%cables(1)%capacitance_per_meter(2,2))
      allocate(expected%mtln%cables(1)%resistance_per_meter(2,2))
      allocate(expected%mtln%cables(1)%conductance_per_meter(2,2))

      expected%mtln%cables(1)%inductance_per_meter = & 
         reshape( source = [ 3.13182309e-07, 7.45674981e-08, 7.45674981e-08, 3.13182309e-07 ], shape = [ 2,2 ] )
      expected%mtln%cables(1)%capacitance_per_meter = &
         reshape( source = [85.0e-12, -20.5e-12, -20.5e-12, 85.0e-12 ], shape = [ 2,2 ] )
      expected%mtln%cables(1)%resistance_per_meter =  reshape(source=[0.0, 0.0, 0.0, 0.0], shape=[2,2])
      expected%mtln%cables(1)%conductance_per_meter = reshape(source=[0.0, 0.0, 0.0, 0.0], shape=[2,2])
      
      allocate(expected%mtln%cables(1)%step_size(18))
      expected%mtln%cables(1)%step_size = [(0.03, i = 1, 18)]
      allocate(expected%mtln%cables(1)%external_field_segments(18))
      do i = 1, 18
         expected%mtln%cables(1)%external_field_segments(i)%position = (/2,2,i/)
         expected%mtln%cables(1)%external_field_segments(i)%direction = DIRECTION_Z_POS
         expected%mtln%cables(1)%external_field_segments(i)%field => null()
      end do
      expected%mtln%cables(1)%transfer_impedance%direction = TRANSFER_IMPEDANCE_DIRECTION_INWARDS
      expected%mtln%cables(1)%transfer_impedance%resistive_term = 0.0
      expected%mtln%cables(1)%transfer_impedance%inductive_term = 4.0e-9
      allocate(expected%mtln%cables(1)%transfer_impedance%poles(0))
      allocate(expected%mtln%cables(1)%transfer_impedance%residues(0))

      expected%mtln%cables(1)%parent_cable => expected%mtln%cables(2)
      expected%mtln%cables(1)%conductor_in_parent = 1
      expected%mtln%cables(1)%initial_connector => null()
      expected%mtln%cables(1)%end_connector => null()

      ! probes
      deallocate(expected%mtln%probes)
      allocate(expected%mtln%probes(2))
      expected%mtln%probes(1)%attached_to_cable => expected%mtln%cables(2) ! to which cable is the probe attached in mtln?
      expected%mtln%probes(1)%index = 19
      expected%mtln%probes(1)%probe_type = PROBE_TYPE_VOLTAGE

      expected%mtln%probes(2)%attached_to_cable => expected%mtln%cables(2)
      expected%mtln%probes(2)%index = 19
      expected%mtln%probes(2)%probe_type = PROBE_TYPE_CURRENT

      ! networks
      deallocate(expected%mtln%networks)
      allocate(expected%mtln%networks(2))
      allocate(expected%mtln%networks(1)%connections(3))

      allocate(expected%mtln%networks(1)%connections(3)%nodes(1))
      expected%mtln%networks(1)%connections(3)%nodes(1)%conductor_in_cable = 1
      expected%mtln%networks(1)%connections(3)%nodes(1)%side = TERMINAL_NODE_SIDE_INI
      expected%mtln%networks(1)%connections(3)%nodes(1)%belongs_to_cable =>  expected%mtln%cables(2)
      expected%mtln%networks(1)%connections(3)%nodes(1)%termination%termination_type = TERMINATION_SERIES
      expected%mtln%networks(1)%connections(3)%nodes(1)%termination%resistance = 50.0

      allocate(expected%mtln%networks(1)%connections(1)%nodes(1))
      expected%mtln%networks(1)%connections(1)%nodes(1)%conductor_in_cable = 1
      expected%mtln%networks(1)%connections(1)%nodes(1)%side = TERMINAL_NODE_SIDE_INI
      expected%mtln%networks(1)%connections(1)%nodes(1)%belongs_to_cable => expected%mtln%cables(1)
      expected%mtln%networks(1)%connections(1)%nodes(1)%termination%termination_type = TERMINATION_SERIES
      expected%mtln%networks(1)%connections(1)%nodes(1)%termination%resistance = 50.0

      allocate(expected%mtln%networks(1)%connections(2)%nodes(1))
      expected%mtln%networks(1)%connections(2)%nodes(1)%conductor_in_cable = 2
      expected%mtln%networks(1)%connections(2)%nodes(1)%side = TERMINAL_NODE_SIDE_INI
      expected%mtln%networks(1)%connections(2)%nodes(1)%belongs_to_cable => expected%mtln%cables(1)
      expected%mtln%networks(1)%connections(2)%nodes(1)%termination%termination_type = TERMINATION_SERIES
      expected%mtln%networks(1)%connections(2)%nodes(1)%termination%resistance = 50.0

      allocate(expected%mtln%networks(2)%connections(3))

      allocate(expected%mtln%networks(2)%connections(3)%nodes(1))
      expected%mtln%networks(2)%connections(3)%nodes(1)%conductor_in_cable = 1
      expected%mtln%networks(2)%connections(3)%nodes(1)%side = TERMINAL_NODE_SIDE_END
      expected%mtln%networks(2)%connections(3)%nodes(1)%belongs_to_cable => expected%mtln%cables(2)
      expected%mtln%networks(2)%connections(3)%nodes(1)%termination%termination_type = TERMINATION_SERIES
      expected%mtln%networks(2)%connections(3)%nodes(1)%termination%resistance = 50.0

      allocate(expected%mtln%networks(2)%connections(1)%nodes(1))
      expected%mtln%networks(2)%connections(1)%nodes(1)%conductor_in_cable = 1
      expected%mtln%networks(2)%connections(1)%nodes(1)%side = TERMINAL_NODE_SIDE_END
      expected%mtln%networks(2)%connections(1)%nodes(1)%belongs_to_cable => expected%mtln%cables(1)
      expected%mtln%networks(2)%connections(1)%nodes(1)%termination%termination_type = TERMINATION_SERIES
      expected%mtln%networks(2)%connections(1)%nodes(1)%termination%resistance = 50.0

      allocate(expected%mtln%networks(2)%connections(2)%nodes(1))
      expected%mtln%networks(2)%connections(2)%nodes(1)%conductor_in_cable = 2
      expected%mtln%networks(2)%connections(2)%nodes(1)%side = TERMINAL_NODE_SIDE_END
      expected%mtln%networks(2)%connections(2)%nodes(1)%belongs_to_cable => expected%mtln%cables(1)
      expected%mtln%networks(2)%connections(2)%nodes(1)%termination%termination_type = TERMINATION_SERIES
      expected%mtln%networks(2)%connections(2)%nodes(1)%termination%resistance = 50.0

      !connectors
      allocate(expected%mtln%connectors(0))

   end function
end function

