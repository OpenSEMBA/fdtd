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
      ! allocate(expected%nodSrc)
      allocate(expected%nodSrc%NodalSource(1))
      expected%nodSrc%n_C1P_max = 0
      expected%nodSrc%n_C2P_max = 1
      expected%nodSrc%n_nodSrc = 1
      expected%nodSrc%n_nodSrc_max = 1
      expected%nodSrc%NodalSource(1)%nombre = trim(adjustl("gauss.exc"))
      expected%nodSrc%NodalSource(1)%isElec = .true.
      expected%nodSrc%NodalSource(1)%isMagnet = .false.
      expected%nodSrc%NodalSource(1)%isCurrent = .false.
      expected%nodSrc%NodalSource(1)%isField = .true.
      expected%nodSrc%NodalSource(1)%isInitialValue = .false.
      expected%nodSrc%NodalSource(1)%n_C1P = 0
      allocate(expected%nodSrc%NodalSource(1)%c1P(0))
      expected%nodSrc%NodalSource(1)%n_C2P = 1
      allocate(expected%nodSrc%NodalSource(1)%c2P(1))
      expected%nodSrc%NodalSource(1)%c2P(1)%xi = 2
      expected%nodSrc%NodalSource(1)%c2P(1)%xe = 3
      expected%nodSrc%NodalSource(1)%c2P(1)%yi = 9
      expected%nodSrc%NodalSource(1)%c2P(1)%ye = 9
      expected%nodSrc%NodalSource(1)%c2P(1)%zi = 1
      expected%nodSrc%NodalSource(1)%c2P(1)%ze = 1
      expected%nodSrc%NodalSource(1)%c2P(1)%xc = 1
      expected%nodSrc%NodalSource(1)%c2P(1)%yc = 0
      expected%nodSrc%NodalSource(1)%c2P(1)%zc = 0
      expected%nodSrc%NodalSource(1)%c2P(1)%or = 1
      expected%nodSrc%NodalSource(1)%c2P(1)%tag = " "
      
      

      ! Expected probes
      ! oldSonda
      expected%oldSONDA%n_probes = 0
      expected%oldSONDA%n_probes_max = 0
      allocate(expected%oldSONDA%probes(0))

      ! sonda
      expected%Sonda%length = 4
      expected%Sonda%length_max = 4
      allocate(expected%Sonda%collection(4))
            
      expected%Sonda%collection(1)%outputrequest = "b1_terminal_current"
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
      expected%Sonda%collection(1)%cordinates(1)%tag = "b1_terminal_current"
      expected%Sonda%collection(1)%cordinates(1)%Xi = 11 ! Coord id as tag.
      expected%Sonda%collection(1)%cordinates(1)%Yi = 0
      expected%Sonda%collection(1)%cordinates(1)%Zi = 0
      expected%Sonda%collection(1)%cordinates(1)%Or = NP_COR_WIRECURRENT
            
      expected%Sonda%collection(2)%outputrequest = "junction_current"
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
      expected%Sonda%collection(2)%cordinates(1)%tag = "junction_current"
      expected%Sonda%collection(2)%cordinates(1)%Xi = 15 ! Coord id as tag.
      expected%Sonda%collection(2)%cordinates(1)%Yi = 0
      expected%Sonda%collection(2)%cordinates(1)%Zi = 0
      expected%Sonda%collection(2)%cordinates(1)%Or = NP_COR_WIRECURRENT
            
            
      expected%Sonda%collection(3)%outputrequest = "b2_terminal_current"
      expected%Sonda%collection(3)%type1 = NP_T1_PLAIN
      expected%Sonda%collection(3)%type2 = NP_T2_TIME
      expected%Sonda%collection(3)%filename = ' '
      expected%Sonda%collection(3)%tstart = 0.0
      expected%Sonda%collection(3)%tstop = 0.0
      expected%Sonda%collection(3)%tstep = 0.0
      expected%Sonda%collection(3)%fstart = 0.0
      expected%Sonda%collection(3)%fstop = 0.0
      expected%Sonda%collection(3)%fstep = 0.0
      allocate(expected%Sonda%collection(3)%cordinates(1))
      expected%Sonda%collection(3)%len_cor = 1
      expected%Sonda%collection(3)%cordinates(1)%tag = "b2_terminal_current"
      expected%Sonda%collection(3)%cordinates(1)%Xi = 23 ! Coord id as tag.
      expected%Sonda%collection(3)%cordinates(1)%Yi = 0
      expected%Sonda%collection(3)%cordinates(1)%Zi = 0
      expected%Sonda%collection(3)%cordinates(1)%Or = NP_COR_WIRECURRENT
            
           
      expected%Sonda%collection(4)%outputrequest = "b3_terminal_current"
      expected%Sonda%collection(4)%type1 = NP_T1_PLAIN
      expected%Sonda%collection(4)%type2 = NP_T2_TIME
      expected%Sonda%collection(4)%filename = ' '
      expected%Sonda%collection(4)%tstart = 0.0
      expected%Sonda%collection(4)%tstop = 0.0
      expected%Sonda%collection(4)%tstep = 0.0
      expected%Sonda%collection(4)%fstart = 0.0
      expected%Sonda%collection(4)%fstop = 0.0
      expected%Sonda%collection(4)%fstep = 0.0
      allocate(expected%Sonda%collection(4)%cordinates(1))
      expected%Sonda%collection(4)%len_cor = 1
      expected%Sonda%collection(4)%cordinates(1)%tag = "b3_terminal_current"
      expected%Sonda%collection(4)%cordinates(1)%Xi = 24 ! Coord id as tag.
      expected%Sonda%collection(4)%cordinates(1)%Yi = 0
      expected%Sonda%collection(4)%cordinates(1)%Zi = 0
      expected%Sonda%collection(4)%cordinates(1)%Or = NP_COR_WIRECURRENT
            
      ! Expected thin wires
      allocate(expected%tWires%tw(3))
      ! wire 1
      expected%tWires%tw(1)%rad=0.0001
      expected%tWires%tw(1)%res=22.9e-3
      expected%tWires%tw(1)%dispfile = trim(adjustl(" "))
      expected%tWires%tw(1)%dispfile_LeftEnd = trim(adjustl(" "))
      expected%tWires%tw(1)%dispfile_RightEnd = trim(adjustl(" "))
      expected%tWires%tw(1)%n_twc=9
      expected%tWires%tw(1)%n_twc_max=9
      allocate(expected%tWires%tw(1)%twc(9))
      expected%tWires%tw(1)%twc(1:9)%srcfile = 'None'
      expected%tWires%tw(1)%twc(1:9)%srctype = 'None'
      expected%tWires%tw(1)%twc(1:9)%i = [(i, i=1, 9)]
      expected%tWires%tw(1)%twc(1:9)%j = 9
      expected%tWires%tw(1)%twc(1:9)%k = 1
      expected%tWires%tw(1)%twc(1:9)%d = DIR_X
      expected%tWires%tw(1)%twc(1)%nd  = 1
      expected%tWires%tw(1)%twc(2:8)%nd = NO_TAG
      expected%tWires%tw(1)%twc(9)%nd  = 2
      
      expected%tWires%tw(1)%twc(1:9)%tag = trim(adjustl("1"))   ! The polyline id is used as tag.
      
      expected%tWires%tw(1)%tl = SERIES_CONS
      expected%tWires%tw(1)%R_LeftEnd = 0.7e-3
      expected%tWires%tw(1)%C_LeftEnd = 1e22
      expected%tWires%tw(1)%tr = SERIES_CONS
      expected%tWires%tw(1)%R_RightEnd = 1e-6
      expected%tWires%tw(1)%C_RightEnd = 1e22

      ! wire 2
      expected%tWires%tw(2)%rad=0.0001
      expected%tWires%tw(2)%res=11.8e-3
      expected%tWires%tw(2)%dispfile = trim(adjustl(" "))
      expected%tWires%tw(2)%dispfile_LeftEnd = trim(adjustl(" "))
      expected%tWires%tw(2)%dispfile_RightEnd = trim(adjustl(" "))
      expected%tWires%tw(2)%n_twc=8
      expected%tWires%tw(2)%n_twc_max=8
      allocate(expected%tWires%tw(2)%twc(8))
      expected%tWires%tw(2)%twc(1:8)%srcfile = 'None'
      expected%tWires%tw(2)%twc(1:8)%srctype = 'None'
      expected%tWires%tw(2)%twc(1:8)%i = [(i, i=10, 17)]
      expected%tWires%tw(2)%twc(1:8)%j = 9
      expected%tWires%tw(2)%twc(1:8)%k = 1
      expected%tWires%tw(2)%twc(1:8)%d = DIR_X
      expected%tWires%tw(2)%twc(1)%nd  = 2
      expected%tWires%tw(2)%twc(2:7)%nd = NO_TAG
      expected%tWires%tw(2)%twc(8)%nd  = 5
      
      expected%tWires%tw(2)%twc(1:8)%tag = trim(adjustl("2"))   ! The polyline id is used as tag.
      
      expected%tWires%tw(2)%tl = SERIES_CONS
      expected%tWires%tw(2)%R_LeftEnd = 1e-6
      expected%tWires%tw(2)%C_LeftEnd = 1e22
      expected%tWires%tw(2)%tr = SERIES_CONS
      expected%tWires%tw(2)%R_RightEnd = 1.0
      expected%tWires%tw(2)%C_RightEnd = 1e22

      ! wire 3
      expected%tWires%tw(3)%rad=0.0001
      expected%tWires%tw(3)%res= 17.3e-3
      expected%tWires%tw(3)%dispfile = trim(adjustl(" "))
      expected%tWires%tw(3)%dispfile_LeftEnd = trim(adjustl(" "))
      expected%tWires%tw(3)%dispfile_RightEnd = trim(adjustl(" "))
      expected%tWires%tw(3)%n_twc=7
      expected%tWires%tw(3)%n_twc_max=7
      allocate(expected%tWires%tw(3)%twc(7))
      expected%tWires%tw(3)%twc(1:7)%srcfile = 'None'
      expected%tWires%tw(3)%twc(1:7)%srctype = 'None'
      expected%tWires%tw(3)%twc(1:7)%i = 10
      expected%tWires%tw(3)%twc(1:7)%j = [(i, i=8, 2, -1)]
      expected%tWires%tw(3)%twc(1:7)%k = 1
      expected%tWires%tw(3)%twc(1:7)%d = DIR_Y
      expected%tWires%tw(3)%twc(1)%nd  = 2
      expected%tWires%tw(3)%twc(2:6)%nd = NO_TAG
      expected%tWires%tw(3)%twc(7)%nd  = 6
      
      expected%tWires%tw(3)%twc(1:7)%tag = trim(adjustl("3"))   ! The polyline id is used as tag.
      
      expected%tWires%tw(3)%tl = SERIES_CONS
      expected%tWires%tw(3)%R_LeftEnd = 1e-6
      expected%tWires%tw(3)%C_LeftEnd = 1e22
      expected%tWires%tw(3)%tr = SERIES_CONS
      expected%tWires%tw(3)%R_RightEnd = 0.7e-3
      expected%tWires%tw(3)%C_RightEnd = 1e22

      expected%tWires%n_tw = 3
      expected%tWires%n_tw_max = 3

      ! Expected mtln type
      !connectors
      ! id = 24
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
      allocate(expected%mtln%cables(9))
      ! level 0
      ! cable 1 - wire
      expected%mtln%cables(1)%name = "line_0_0"
      allocate(expected%mtln%cables(1)%inductance_per_meter(1,1))
      allocate(expected%mtln%cables(1)%capacitance_per_meter(1,1))
      allocate(expected%mtln%cables(1)%resistance_per_meter(1,1))
      allocate(expected%mtln%cables(1)%conductance_per_meter(1,1))

      expected%mtln%cables(1)%inductance_per_meter = reshape( source = [5.481553487168089e-07], shape = [ 1,1 ] )
      expected%mtln%cables(1)%capacitance_per_meter = reshape( source = [2.0270004E-11], shape = [ 1,1 ] )
      expected%mtln%cables(1)%resistance_per_meter =  reshape(source=[22.9e-3], shape=[1,1])
      expected%mtln%cables(1)%conductance_per_meter = reshape(source=[0.0], shape=[1,1])
      
      allocate(expected%mtln%cables(1)%step_size(9))
      expected%mtln%cables(1)%step_size = [(0.1, i = 1, 9)]
      allocate(expected%mtln%cables(1)%external_field_segments(9))
      do i = 1, 9
         expected%mtln%cables(1)%external_field_segments(i)%position = (/i,9,1/)
         expected%mtln%cables(1)%external_field_segments(i)%direction = DIRECTION_X_POS
         expected%mtln%cables(1)%external_field_segments(i)%field => null()
      end do

      allocate(expected%mtln%cables(1)%transfer_impedance%poles(0))
      allocate(expected%mtln%cables(1)%transfer_impedance%residues(0))

      expected%mtln%cables(1)%parent_cable => null()
      expected%mtln%cables(1)%conductor_in_parent = 0
      expected%mtln%cables(1)%initial_connector => expected%mtln%connectors(1)
      expected%mtln%cables(1)%end_connector => null()

      ! cable 2 - multiwire
      expected%mtln%cables(2)%name = "line_1_0"
      allocate(expected%mtln%cables(2)%inductance_per_meter(1,1))
      allocate(expected%mtln%cables(2)%capacitance_per_meter(1,1))
      allocate(expected%mtln%cables(2)%resistance_per_meter(1,1))
      allocate(expected%mtln%cables(2)%conductance_per_meter(1,1))
      expected%mtln%cables(2)%inductance_per_meter = reshape(source=[8.802075200000001e-08], shape=[1,1])
      expected%mtln%cables(2)%capacitance_per_meter = reshape(source=[5.5840010E-10], shape=[1,1])
      expected%mtln%cables(2)%resistance_per_meter = reshape(source=[3.9e-3], shape=[1,1])
      expected%mtln%cables(2)%conductance_per_meter = reshape(source=[0.0], shape=[1,1])
      allocate(expected%mtln%cables(2)%step_size(9))
      expected%mtln%cables(2)%step_size =  [(0.1, i = 1, 9)]

      allocate(expected%mtln%cables(2)%external_field_segments(9))
      do i = 1, 9
         expected%mtln%cables(2)%external_field_segments(i)%position = (/i,9,1/)
         expected%mtln%cables(2)%external_field_segments(i)%direction = DIRECTION_X_POS
         expected%mtln%cables(2)%external_field_segments(i)%field => null()

      end do

      expected%mtln%cables(2)%transfer_impedance%direction = TRANSFER_IMPEDANCE_DIRECTION_INWARDS
      expected%mtln%cables(2)%transfer_impedance%resistive_term = 0.0
      expected%mtln%cables(2)%transfer_impedance%inductive_term = 8.9e-9
      allocate(expected%mtln%cables(2)%transfer_impedance%poles(0))
      allocate(expected%mtln%cables(2)%transfer_impedance%residues(0))

      expected%mtln%cables(2)%parent_cable => expected%mtln%cables(1)
      expected%mtln%cables(2)%conductor_in_parent = 1
      expected%mtln%cables(2)%initial_connector => expected%mtln%connectors(3)
      expected%mtln%cables(2)%end_connector => null()

      ! cable 3 - multiwire
      expected%mtln%cables(3)%name = "line_2_0"
      allocate(expected%mtln%cables(3)%inductance_per_meter(8,8), source = 0.0)
      allocate(expected%mtln%cables(3)%capacitance_per_meter(8,8), source = 0.0)
      allocate(expected%mtln%cables(3)%resistance_per_meter(8,8), source = 0.0)
      allocate(expected%mtln%cables(3)%conductance_per_meter(8,8), source = 0.0)
      expected%mtln%cables(3)%inductance_per_meter(1:2,1:2) = & 
         reshape(source=[2.4382084E-07, 4.7377505E-08, 4.7377508E-08, 2.4382081E-07], shape=[2,2], order =[2,1])
      expected%mtln%cables(3)%inductance_per_meter(3:4,3:4) = & 
         reshape(source=[2.4382084E-07, 4.7377505E-08, 4.7377508E-08, 2.4382081E-07], shape=[2,2], order =[2,1])
      expected%mtln%cables(3)%inductance_per_meter(5:6,5:6) = & 
         reshape(source=[2.4382084E-07, 4.7377505E-08, 4.7377508E-08, 2.4382081E-07], shape=[2,2], order =[2,1])
      expected%mtln%cables(3)%inductance_per_meter(7:8,7:8) = & 
         reshape(source=[2.4382084E-07, 4.7377505E-08, 4.7377508E-08, 2.4382081E-07], shape=[2,2], order =[2,1])

      expected%mtln%cables(3)%capacitance_per_meter(1:2,1:2) = &
         reshape(source=[105.5e-12, -20.5e-12, -20.5e-12, 105.5e-12], shape=[2,2], order =[2,1])
      expected%mtln%cables(3)%capacitance_per_meter(3:4,3:4) = &
         reshape(source=[105.5e-12, -20.5e-12, -20.5e-12, 105.5e-12], shape=[2,2], order =[2,1])
      expected%mtln%cables(3)%capacitance_per_meter(5:6,5:6) = &
         reshape(source=[105.5e-12, -20.5e-12, -20.5e-12, 105.5e-12], shape=[2,2], order =[2,1])
      expected%mtln%cables(3)%capacitance_per_meter(7:8,7:8) = &
         reshape(source=[105.5e-12, -20.5e-12, -20.5e-12, 105.5e-12], shape=[2,2], order =[2,1])

      do i = 1, 8
         expected%mtln%cables(3)%resistance_per_meter(i,i) =  62.0e-3
      end do

      allocate(expected%mtln%cables(3)%step_size(9))
      expected%mtln%cables(3)%step_size =  [(0.1, i = 1, 9)]      

      allocate(expected%mtln%cables(3)%external_field_segments(9))
      do i = 1, 9
         expected%mtln%cables(3)%external_field_segments(i)%position = (/i,9,1/)
         expected%mtln%cables(3)%external_field_segments(i)%direction = DIRECTION_X_POS
         expected%mtln%cables(3)%external_field_segments(i)%field => null()

      end do


      expected%mtln%cables(3)%transfer_impedance%direction = TRANSFER_IMPEDANCE_DIRECTION_INWARDS
      expected%mtln%cables(3)%transfer_impedance%resistive_term = 0.0
      expected%mtln%cables(3)%transfer_impedance%inductive_term = 4.2e-9
      allocate(expected%mtln%cables(3)%transfer_impedance%poles(0))
      allocate(expected%mtln%cables(3)%transfer_impedance%residues(0))

      expected%mtln%cables(3)%parent_cable => expected%mtln%cables(2)
      expected%mtln%cables(3)%conductor_in_parent = 1
      expected%mtln%cables(3)%initial_connector => null()
      expected%mtln%cables(3)%end_connector => null()


      ! cable 4 - wire
      expected%mtln%cables(4)%name = "line_0_1"
      allocate(expected%mtln%cables(4)%inductance_per_meter(1,1))
      allocate(expected%mtln%cables(4)%capacitance_per_meter(1,1))
      allocate(expected%mtln%cables(4)%resistance_per_meter(1,1))
      allocate(expected%mtln%cables(4)%conductance_per_meter(1,1))

      expected%mtln%cables(4)%inductance_per_meter = reshape( source = [6.482560773828984e-07], shape = [ 1,1 ] )
      expected%mtln%cables(4)%capacitance_per_meter = reshape( source = [1.7140003E-11], shape = [ 1,1 ] )
      expected%mtln%cables(4)%resistance_per_meter =  reshape(source=[11.8e-3], shape=[1,1])
      expected%mtln%cables(4)%conductance_per_meter = reshape(source=[0.0], shape=[1,1])
      
      allocate(expected%mtln%cables(4)%step_size(8))
      expected%mtln%cables(4)%step_size = [(0.1, i = 1, 8)]

      allocate(expected%mtln%cables(4)%external_field_segments(8))
      do i = 1, 8
         expected%mtln%cables(4)%external_field_segments(i)%position = (/9+i,9,1/)
         expected%mtln%cables(4)%external_field_segments(i)%direction = DIRECTION_X_POS
         expected%mtln%cables(4)%external_field_segments(i)%field => null()
      end do


      allocate(expected%mtln%cables(4)%transfer_impedance%poles(0))
      allocate(expected%mtln%cables(4)%transfer_impedance%residues(0))

      expected%mtln%cables(4)%parent_cable => null()
      expected%mtln%cables(4)%conductor_in_parent = 0
      expected%mtln%cables(4)%initial_connector => expected%mtln%connectors(2)
      expected%mtln%cables(4)%end_connector => null()

      ! cable 5 - multiwire
      expected%mtln%cables(5)%name = "line_1_1"
      allocate(expected%mtln%cables(5)%inductance_per_meter(1,1))
      allocate(expected%mtln%cables(5)%capacitance_per_meter(1,1))
      allocate(expected%mtln%cables(5)%resistance_per_meter(1,1))
      allocate(expected%mtln%cables(5)%conductance_per_meter(1,1))
      expected%mtln%cables(5)%inductance_per_meter = reshape(source=[1.37228e-07], shape=[1,1])
      expected%mtln%cables(5)%capacitance_per_meter = reshape(source=[3.2310005E-10], shape=[1,1])
      expected%mtln%cables(5)%resistance_per_meter = reshape(source=[12.2e-3], shape=[1,1])
      expected%mtln%cables(5)%conductance_per_meter = reshape(source=[0.0], shape=[1,1])
      allocate(expected%mtln%cables(5)%step_size(8))
      expected%mtln%cables(5)%step_size =  [(0.1, i = 1, 8)]

      allocate(expected%mtln%cables(5)%external_field_segments(8))
      do i = 1, 8
         expected%mtln%cables(5)%external_field_segments(i)%position = (/9+i,9,1/)
         expected%mtln%cables(5)%external_field_segments(i)%direction = DIRECTION_X_POS
         expected%mtln%cables(5)%external_field_segments(i)%field => null()
      end do


      expected%mtln%cables(5)%transfer_impedance%direction = TRANSFER_IMPEDANCE_DIRECTION_INWARDS
      expected%mtln%cables(5)%transfer_impedance%resistive_term = 0.0
      expected%mtln%cables(5)%transfer_impedance%inductive_term = 7.4e-9
      allocate(expected%mtln%cables(5)%transfer_impedance%poles(0))
      allocate(expected%mtln%cables(5)%transfer_impedance%residues(0))

      expected%mtln%cables(5)%parent_cable => expected%mtln%cables(4)
      expected%mtln%cables(5)%conductor_in_parent = 1
      expected%mtln%cables(5)%initial_connector => expected%mtln%connectors(4)
      expected%mtln%cables(5)%end_connector => null()

      ! cable 6 - multiwire
      expected%mtln%cables(6)%name = "line_2_4"
      allocate(expected%mtln%cables(6)%inductance_per_meter(2,2), source = 0.0)
      allocate(expected%mtln%cables(6)%capacitance_per_meter(2,2), source = 0.0)
      allocate(expected%mtln%cables(6)%resistance_per_meter(2,2), source = 0.0)
      allocate(expected%mtln%cables(6)%conductance_per_meter(2,2), source = 0.0)
      expected%mtln%cables(6)%inductance_per_meter(1:2,1:2) = & 
         reshape(source=[2.4382084E-07, 4.7377505E-08, 4.7377508E-08, 2.4382081E-07], shape=[2,2], order =[2,1])

      expected%mtln%cables(6)%capacitance_per_meter(1:2,1:2) = &
         reshape(source=[105.5e-12, -20.5e-12, -20.5e-12, 105.5e-12], shape=[2,2], order =[2,1])

      do i = 1, 2
         expected%mtln%cables(6)%resistance_per_meter(i,i) = 62.0e-3
      end do

      allocate(expected%mtln%cables(6)%step_size(8))
      expected%mtln%cables(6)%step_size =  [(0.1, i = 1, 8)]

      allocate(expected%mtln%cables(6)%external_field_segments(8))
      do i = 1, 8
         expected%mtln%cables(6)%external_field_segments(i)%position = (/9+i,9,1/)
         expected%mtln%cables(6)%external_field_segments(i)%direction = DIRECTION_X_POS
         expected%mtln%cables(6)%external_field_segments(i)%field => null()
      end do


      expected%mtln%cables(6)%transfer_impedance%direction = TRANSFER_IMPEDANCE_DIRECTION_INWARDS
      expected%mtln%cables(6)%transfer_impedance%resistive_term = 0.0
      expected%mtln%cables(6)%transfer_impedance%inductive_term = 4.2e-9
      allocate(expected%mtln%cables(6)%transfer_impedance%poles(0))
      allocate(expected%mtln%cables(6)%transfer_impedance%residues(0))

      expected%mtln%cables(6)%parent_cable => expected%mtln%cables(5)
      expected%mtln%cables(6)%conductor_in_parent = 1
      expected%mtln%cables(6)%initial_connector => null()
      expected%mtln%cables(6)%end_connector => null()


      ! cable 7 - wire
      expected%mtln%cables(7)%name = "line_0_2"
      allocate(expected%mtln%cables(7)%inductance_per_meter(1,1))
      allocate(expected%mtln%cables(7)%capacitance_per_meter(1,1))
      allocate(expected%mtln%cables(7)%resistance_per_meter(1,1))
      allocate(expected%mtln%cables(7)%conductance_per_meter(1,1))

      expected%mtln%cables(7)%inductance_per_meter = reshape( source = [5.802145885361537e-07], shape = [ 1,1 ] )
      expected%mtln%cables(7)%capacitance_per_meter = reshape( source = [1.9150003E-11], shape = [ 1,1 ] )
      expected%mtln%cables(7)%resistance_per_meter =  reshape(source=[17.3e-3], shape=[1,1])
      expected%mtln%cables(7)%conductance_per_meter = reshape(source=[0.0], shape=[1,1])
      
      allocate(expected%mtln%cables(7)%step_size(7))
      expected%mtln%cables(7)%step_size = [(0.1, i = 1, 7)]

      allocate(expected%mtln%cables(7)%external_field_segments(7))
      do i = 1,7
         expected%mtln%cables(7)%external_field_segments(i)%position = (/10,9-i,1/)
         expected%mtln%cables(7)%external_field_segments(i)%direction = DIRECTION_Y_NEG
         expected%mtln%cables(7)%external_field_segments(i)%field => null()
      end do


      allocate(expected%mtln%cables(7)%transfer_impedance%poles(0))
      allocate(expected%mtln%cables(7)%transfer_impedance%residues(0))

      expected%mtln%cables(7)%parent_cable => null()
      expected%mtln%cables(7)%conductor_in_parent = 0
      expected%mtln%cables(7)%initial_connector => null()
      expected%mtln%cables(7)%end_connector => null()

      ! cable 8 - multiwire
      expected%mtln%cables(8)%name = "line_1_2"
      allocate(expected%mtln%cables(8)%inductance_per_meter(1,1))
      allocate(expected%mtln%cables(8)%capacitance_per_meter(1,1))
      allocate(expected%mtln%cables(8)%resistance_per_meter(1,1))
      allocate(expected%mtln%cables(8)%conductance_per_meter(1,1))
      expected%mtln%cables(8)%inductance_per_meter = reshape(source=[9.1890502e-08], shape=[1,1])
      expected%mtln%cables(8)%capacitance_per_meter = reshape(source=[4.7190007E-10], shape=[1,1])
      expected%mtln%cables(8)%resistance_per_meter = reshape(source=[6.5e-3], shape=[1,1])
      expected%mtln%cables(8)%conductance_per_meter = reshape(source=[0.0], shape=[1,1])
      allocate(expected%mtln%cables(8)%step_size(7))
      expected%mtln%cables(8)%step_size =  [(0.1, i = 1, 7)]

      allocate(expected%mtln%cables(8)%external_field_segments(7))
      do i = 1,7
         expected%mtln%cables(8)%external_field_segments(i)%position = (/10,9-i,1/)
         expected%mtln%cables(8)%external_field_segments(i)%direction = DIRECTION_Y_NEG
         expected%mtln%cables(8)%external_field_segments(i)%field => null()
      end do

      expected%mtln%cables(8)%transfer_impedance%direction = TRANSFER_IMPEDANCE_DIRECTION_INWARDS
      expected%mtln%cables(8)%transfer_impedance%resistive_term = 0.0
      expected%mtln%cables(8)%transfer_impedance%inductive_term = 3.0e-9
      allocate(expected%mtln%cables(8)%transfer_impedance%poles(0))
      allocate(expected%mtln%cables(8)%transfer_impedance%residues(0))

      expected%mtln%cables(8)%parent_cable => expected%mtln%cables(7)
      expected%mtln%cables(8)%conductor_in_parent = 1
      expected%mtln%cables(8)%initial_connector => null()
      expected%mtln%cables(8)%end_connector => null()

      ! cable 9 - multiwire
      expected%mtln%cables(9)%name = "line_2_5"
      allocate(expected%mtln%cables(9)%inductance_per_meter(6,6), source = 0.0)
      allocate(expected%mtln%cables(9)%capacitance_per_meter(6,6), source = 0.0)
      allocate(expected%mtln%cables(9)%resistance_per_meter(6,6), source = 0.0)
      allocate(expected%mtln%cables(9)%conductance_per_meter(6,6), source = 0.0)
      expected%mtln%cables(9)%inductance_per_meter(1:2,1:2) = & 
         reshape(source=[2.4382084E-07, 4.7377505E-08, 4.7377508E-08, 2.4382081E-07], shape=[2,2], order =[2,1])
      expected%mtln%cables(9)%inductance_per_meter(3:4,3:4) = & 
         reshape(source=[2.4382084E-07, 4.7377505E-08, 4.7377508E-08, 2.4382081E-07], shape=[2,2], order =[2,1])
      expected%mtln%cables(9)%inductance_per_meter(5:6,5:6) = & 
         reshape(source=[2.4382084E-07, 4.7377505E-08, 4.7377508E-08, 2.4382081E-07], shape=[2,2], order =[2,1])

      expected%mtln%cables(9)%capacitance_per_meter(1:2,1:2) = &
         reshape(source=[105.5e-12, -20.5e-12, -20.5e-12, 105.5e-12], shape=[2,2], order =[2,1])
      expected%mtln%cables(9)%capacitance_per_meter(3:4,3:4) = &
         reshape(source=[105.5e-12, -20.5e-12, -20.5e-12, 105.5e-12], shape=[2,2], order =[2,1])
      expected%mtln%cables(9)%capacitance_per_meter(5:6,5:6) = &
         reshape(source=[105.5e-12, -20.5e-12, -20.5e-12, 105.5e-12], shape=[2,2], order =[2,1])

      do i = 1, 6
         expected%mtln%cables(9)%resistance_per_meter(i,i) = 62.0e-3
      end do
      allocate(expected%mtln%cables(9)%step_size(7))
      expected%mtln%cables(9)%step_size =  [(0.1, i = 1, 7)]

      allocate(expected%mtln%cables(9)%external_field_segments(7))
      do i = 1,7
         expected%mtln%cables(9)%external_field_segments(i)%position = (/10,9-i,1/)
         expected%mtln%cables(9)%external_field_segments(i)%direction = DIRECTION_Y_NEG
         expected%mtln%cables(9)%external_field_segments(i)%field => null()
      end do

      expected%mtln%cables(9)%transfer_impedance%direction = TRANSFER_IMPEDANCE_DIRECTION_INWARDS
      expected%mtln%cables(9)%transfer_impedance%resistive_term = 0.0
      expected%mtln%cables(9)%transfer_impedance%inductive_term = 4.2e-9
      allocate(expected%mtln%cables(9)%transfer_impedance%poles(0))
      allocate(expected%mtln%cables(9)%transfer_impedance%residues(0))

      expected%mtln%cables(9)%parent_cable => expected%mtln%cables(8)
      expected%mtln%cables(9)%conductor_in_parent = 1
      expected%mtln%cables(9)%initial_connector => null()
      expected%mtln%cables(9)%end_connector => null()


      ! probes
      deallocate(expected%mtln%probes)
      allocate(expected%mtln%probes(7))
      expected%mtln%probes(1)%attached_to_cable => expected%mtln%cables(1) ! to which cable is the probe attached in mtln?
      expected%mtln%probes(1)%index = 1
      expected%mtln%probes(1)%probe_type = PROBE_TYPE_VOLTAGE
      expected%mtln%probes(1)%probe_name = "b1_terminal_voltage"
      expected%mtln%probes(1)%probe_position = [1,9,1]

      expected%mtln%probes(2)%attached_to_cable => expected%mtln%cables(1)
      expected%mtln%probes(2)%index = 1
      expected%mtln%probes(2)%probe_type = PROBE_TYPE_CURRENT
      expected%mtln%probes(2)%probe_name = "b1_terminal_current"
      expected%mtln%probes(2)%probe_position = [1,9,1]

      expected%mtln%probes(3)%attached_to_cable => expected%mtln%cables(1)
      expected%mtln%probes(3)%index = 10
      expected%mtln%probes(3)%probe_type = PROBE_TYPE_CURRENT
      expected%mtln%probes(3)%probe_name = "junction_current"
      expected%mtln%probes(3)%probe_position = [10, 9, 1]

      expected%mtln%probes(4)%attached_to_cable => expected%mtln%cables(4)
      expected%mtln%probes(4)%index = 1
      expected%mtln%probes(4)%probe_type = PROBE_TYPE_CURRENT
      expected%mtln%probes(4)%probe_name = "junction_current"
      expected%mtln%probes(4)%probe_position = [10, 9, 1]

      expected%mtln%probes(5)%attached_to_cable => expected%mtln%cables(7)
      expected%mtln%probes(5)%index = 1
      expected%mtln%probes(5)%probe_type = PROBE_TYPE_CURRENT
      expected%mtln%probes(5)%probe_name = "junction_current"
      expected%mtln%probes(5)%probe_position = [10, 9, 1]

      expected%mtln%probes(6)%attached_to_cable => expected%mtln%cables(4)
      expected%mtln%probes(6)%index = 9
      expected%mtln%probes(6)%probe_type = PROBE_TYPE_CURRENT
      expected%mtln%probes(6)%probe_name = "b2_terminal_current"
      expected%mtln%probes(6)%probe_position = [ 18, 9, 1]

      expected%mtln%probes(7)%attached_to_cable => expected%mtln%cables(7)
      expected%mtln%probes(7)%index = 8
      expected%mtln%probes(7)%probe_type = PROBE_TYPE_CURRENT
      expected%mtln%probes(7)%probe_name = "b3_terminal_current"
      expected%mtln%probes(7)%probe_position = [10, 2, 1]



      ! networks
      deallocate(expected%mtln%networks)
      allocate(expected%mtln%networks(4))

      ! NETWORK 1
      allocate(expected%mtln%networks(1)%connections(10))

      allocate(expected%mtln%networks(1)%connections(1)%nodes(1))
      expected%mtln%networks(1)%connections(1)%nodes(1)%conductor_in_cable = 1
      expected%mtln%networks(1)%connections(1)%nodes(1)%side = TERMINAL_NODE_SIDE_INI
      expected%mtln%networks(1)%connections(1)%nodes(1)%belongs_to_cable =>  expected%mtln%cables(1)
      expected%mtln%networks(1)%connections(1)%nodes(1)%termination%termination_type = TERMINATION_SERIES
      expected%mtln%networks(1)%connections(1)%nodes(1)%termination%resistance = 0.7e-3

      allocate(expected%mtln%networks(1)%connections(2)%nodes(1))
      expected%mtln%networks(1)%connections(2)%nodes(1)%conductor_in_cable = 1
      expected%mtln%networks(1)%connections(2)%nodes(1)%side = TERMINAL_NODE_SIDE_INI
      expected%mtln%networks(1)%connections(2)%nodes(1)%belongs_to_cable => expected%mtln%cables(2)
      expected%mtln%networks(1)%connections(2)%nodes(1)%termination%termination_type = TERMINATION_SERIES
      expected%mtln%networks(1)%connections(2)%nodes(1)%termination%resistance = 1e-6

      do i = 3, 10
         allocate(expected%mtln%networks(1)%connections(i)%nodes(1))
         expected%mtln%networks(1)%connections(i)%nodes(1)%side = TERMINAL_NODE_SIDE_INI
         expected%mtln%networks(1)%connections(i)%nodes(1)%belongs_to_cable => expected%mtln%cables(3)
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

      expected%mtln%networks(1)%connections(3)%nodes(1)%termination%termination_type = TERMINATION_LCpRs
      expected%mtln%networks(1)%connections(3)%nodes(1)%termination%resistance = 50.0
      expected%mtln%networks(1)%connections(3)%nodes(1)%termination%inductance = 30e-12
      expected%mtln%networks(1)%connections(3)%nodes(1)%termination%capacitance = 60e-9

      expected%mtln%networks(1)%connections(5)%nodes(1)%termination%termination_type = TERMINATION_LCpRs
      expected%mtln%networks(1)%connections(5)%nodes(1)%termination%resistance = 50.0
      expected%mtln%networks(1)%connections(5)%nodes(1)%termination%inductance = 30e-12
      expected%mtln%networks(1)%connections(5)%nodes(1)%termination%capacitance = 60e-9

      expected%mtln%networks(1)%connections(7)%nodes(1)%termination%termination_type = TERMINATION_LCpRs
      expected%mtln%networks(1)%connections(7)%nodes(1)%termination%resistance = 50.0
      expected%mtln%networks(1)%connections(7)%nodes(1)%termination%inductance = 30e-12
      expected%mtln%networks(1)%connections(7)%nodes(1)%termination%capacitance = 60e-9

      expected%mtln%networks(1)%connections(9)%nodes(1)%termination%termination_type = TERMINATION_LCpRs
      expected%mtln%networks(1)%connections(9)%nodes(1)%termination%resistance = 50.0
      expected%mtln%networks(1)%connections(9)%nodes(1)%termination%inductance = 30e-12
      expected%mtln%networks(1)%connections(9)%nodes(1)%termination%capacitance = 60e-9

      ! NETWORK 2
      allocate(expected%mtln%networks(2)%connections(10))
      allocate(expected%mtln%networks(2)%connections(1)%nodes(3))
      expected%mtln%networks(2)%connections(1)%nodes(1)%conductor_in_cable = 1
      expected%mtln%networks(2)%connections(1)%nodes(1)%side = TERMINAL_NODE_SIDE_END
      expected%mtln%networks(2)%connections(1)%nodes(1)%belongs_to_cable =>  expected%mtln%cables(1)

      expected%mtln%networks(2)%connections(1)%nodes(2)%conductor_in_cable = 1
      expected%mtln%networks(2)%connections(1)%nodes(2)%side = TERMINAL_NODE_SIDE_INI
      expected%mtln%networks(2)%connections(1)%nodes(2)%belongs_to_cable =>  expected%mtln%cables(4)

      expected%mtln%networks(2)%connections(1)%nodes(3)%conductor_in_cable = 1
      expected%mtln%networks(2)%connections(1)%nodes(3)%side = TERMINAL_NODE_SIDE_INI
      expected%mtln%networks(2)%connections(1)%nodes(3)%belongs_to_cable =>  expected%mtln%cables(7)

      do i = 1, 3
         expected%mtln%networks(2)%connections(1)%nodes(i)%termination%termination_type = TERMINATION_SERIES
         expected%mtln%networks(2)%connections(1)%nodes(i)%termination%resistance = 1e-6
      end do

      allocate(expected%mtln%networks(2)%connections(2)%nodes(3))
      expected%mtln%networks(2)%connections(2)%nodes(1)%conductor_in_cable = 1
      expected%mtln%networks(2)%connections(2)%nodes(1)%side = TERMINAL_NODE_SIDE_END
      expected%mtln%networks(2)%connections(2)%nodes(1)%belongs_to_cable =>  expected%mtln%cables(2)

      expected%mtln%networks(2)%connections(2)%nodes(2)%conductor_in_cable = 1
      expected%mtln%networks(2)%connections(2)%nodes(2)%side = TERMINAL_NODE_SIDE_INI
      expected%mtln%networks(2)%connections(2)%nodes(2)%belongs_to_cable =>  expected%mtln%cables(5)

      expected%mtln%networks(2)%connections(2)%nodes(3)%conductor_in_cable = 1
      expected%mtln%networks(2)%connections(2)%nodes(3)%side = TERMINAL_NODE_SIDE_INI
      expected%mtln%networks(2)%connections(2)%nodes(3)%belongs_to_cable =>  expected%mtln%cables(8)

      do i = 1, 3
         expected%mtln%networks(2)%connections(2)%nodes(i)%termination%termination_type = TERMINATION_SERIES
         expected%mtln%networks(2)%connections(2)%nodes(i)%termination%resistance = 1e-6
      end do

      ! NETWORK 2 - CONNECTIONS 3-10
      do i = 3, 8
         allocate(expected%mtln%networks(2)%connections(i)%nodes(2))

         expected%mtln%networks(2)%connections(i)%nodes(1)%conductor_in_cable = i-2
         expected%mtln%networks(2)%connections(i)%nodes(1)%side = TERMINAL_NODE_SIDE_END
         expected%mtln%networks(2)%connections(i)%nodes(1)%belongs_to_cable =>  expected%mtln%cables(3)
   
         expected%mtln%networks(2)%connections(i)%nodes(2)%conductor_in_cable = i-2
         expected%mtln%networks(2)%connections(i)%nodes(2)%side = TERMINAL_NODE_SIDE_INI
         expected%mtln%networks(2)%connections(i)%nodes(2)%belongs_to_cable =>  expected%mtln%cables(9)

         do j = 1, 2
            expected%mtln%networks(2)%connections(i)%nodes(j)%termination%termination_type = TERMINATION_SERIES
            expected%mtln%networks(2)%connections(i)%nodes(j)%termination%resistance = 1e-6
         end do

      end do

      do i = 9, 10
         allocate(expected%mtln%networks(2)%connections(i)%nodes(2))

         expected%mtln%networks(2)%connections(i)%nodes(1)%conductor_in_cable = i-2
         expected%mtln%networks(2)%connections(i)%nodes(1)%side = TERMINAL_NODE_SIDE_END
         expected%mtln%networks(2)%connections(i)%nodes(1)%belongs_to_cable =>  expected%mtln%cables(3)
   
         expected%mtln%networks(2)%connections(i)%nodes(2)%conductor_in_cable = i-8
         expected%mtln%networks(2)%connections(i)%nodes(2)%side = TERMINAL_NODE_SIDE_INI
         expected%mtln%networks(2)%connections(i)%nodes(2)%belongs_to_cable =>  expected%mtln%cables(6)

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
      expected%mtln%networks(3)%connections(1)%nodes(1)%belongs_to_cable =>  expected%mtln%cables(4)
      expected%mtln%networks(3)%connections(1)%nodes(1)%termination%termination_type = TERMINATION_SERIES
      expected%mtln%networks(3)%connections(1)%nodes(1)%termination%resistance = 1

      allocate(expected%mtln%networks(3)%connections(2)%nodes(1))
      expected%mtln%networks(3)%connections(2)%nodes(1)%conductor_in_cable = 1
      expected%mtln%networks(3)%connections(2)%nodes(1)%side = TERMINAL_NODE_SIDE_END
      expected%mtln%networks(3)%connections(2)%nodes(1)%belongs_to_cable => expected%mtln%cables(5)
      expected%mtln%networks(3)%connections(2)%nodes(1)%termination%termination_type = TERMINATION_SERIES
      expected%mtln%networks(3)%connections(2)%nodes(1)%termination%resistance = 1e-6

      do i = 3, 4
         allocate(expected%mtln%networks(3)%connections(i)%nodes(1))
         expected%mtln%networks(3)%connections(i)%nodes(1)%side = TERMINAL_NODE_SIDE_END
         expected%mtln%networks(3)%connections(i)%nodes(1)%belongs_to_cable => expected%mtln%cables(6)
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
      expected%mtln%networks(4)%connections(1)%nodes(1)%belongs_to_cable =>  expected%mtln%cables(7)
      expected%mtln%networks(4)%connections(1)%nodes(1)%termination%termination_type = TERMINATION_SERIES
      expected%mtln%networks(4)%connections(1)%nodes(1)%termination%resistance = 0.7e-3

      allocate(expected%mtln%networks(4)%connections(2)%nodes(1))
      expected%mtln%networks(4)%connections(2)%nodes(1)%conductor_in_cable = 1
      expected%mtln%networks(4)%connections(2)%nodes(1)%side = TERMINAL_NODE_SIDE_END
      expected%mtln%networks(4)%connections(2)%nodes(1)%belongs_to_cable => expected%mtln%cables(8)
      expected%mtln%networks(4)%connections(2)%nodes(1)%termination%termination_type = TERMINATION_SERIES
      expected%mtln%networks(4)%connections(2)%nodes(1)%termination%resistance = 1e-6

      do i = 3, 8
         allocate(expected%mtln%networks(4)%connections(i)%nodes(1))
         expected%mtln%networks(4)%connections(i)%nodes(1)%side = TERMINAL_NODE_SIDE_END
         expected%mtln%networks(4)%connections(i)%nodes(1)%belongs_to_cable => expected%mtln%cables(9)
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

