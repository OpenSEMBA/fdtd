integer function test_read_holland1981() bind (C) result(err)
   use smbjson
   use smbjson_testingTools

   implicit none

   character(len=*),parameter :: filename = PATH_TO_TEST_DATA//'cases/holland1981.fdtd.json'
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
      expected%general%dt = 30e-12
      expected%general%nmax = 1000

      ! Excected media matrix.
      expected%matriz%totalX = 20
      expected%matriz%totalY = 20
      expected%matriz%totalZ = 22

      ! Expected grid.
      expected%despl%nX = 20
      expected%despl%nY = 20
      expected%despl%nZ = 22

      allocate(expected%despl%desX(20))
      allocate(expected%despl%desY(20))
      allocate(expected%despl%desZ(22))
      expected%despl%desX = 0.1
      expected%despl%desY = 0.1
      expected%despl%desZ = 0.1
      expected%despl%mx1 = 0
      expected%despl%mx2 = 20
      expected%despl%my1 = 0
      expected%despl%my2 = 20
      expected%despl%mz1 = 0
      expected%despl%mz2 = 22

      ! Expected boundaries.
      expected%front%tipoFrontera(:) = F_PML
      expected%front%propiedadesPML(:)%numCapas = 6
      expected%front%propiedadesPML(:)%orden = 2.0
      expected%front%propiedadesPML(:)%refl = 0.001

      ! Expected sources.
      allocate(expected%plnSrc%collection(1))
      expected%plnSrc%collection(1)%nombre_fichero = "holland.exc"
      expected%plnSrc%collection(1)%atributo = ""
      expected%plnSrc%collection(1)%coor1 = [1, 1, 1]
      expected%plnSrc%collection(1)%coor2 = [18, 18, 20]
      expected%plnSrc%collection(1)%theta = 1.5708
      expected%plnSrc%collection(1)%phi = 0.0
      expected%plnSrc%collection(1)%alpha = 0.0
      expected%plnSrc%collection(1)%beta = 0.0
      expected%plnSrc%collection(1)%isRC=.false.
      expected%plnSrc%collection(1)%nummodes=1
      expected%plnSrc%collection(1)%INCERTMAX=0.0
      expected%plnSrc%nc = 1
      expected%plnSrc%nC_max = 1

      ! Expected probes
      ! sonda
      expected%Sonda%length = 1
      expected%Sonda%length_max = 1
      allocate(expected%Sonda%collection(1))
      expected%Sonda%collection(1)%outputrequest = "mid_point"
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
      expected%Sonda%collection(1)%cordinates(1)%tag = "mid_point"
      expected%Sonda%collection(1)%cordinates(1)%Xi = 2 ! Coord id as tag.
      expected%Sonda%collection(1)%cordinates(1)%Yi = 0
      expected%Sonda%collection(1)%cordinates(1)%Zi = 0
      expected%Sonda%collection(1)%cordinates(1)%Or = NP_COR_WIRECURRENT
      
      
      ! Expected thin wires
      allocate(expected%tWires%tw(1))
      expected%tWires%tw(1)%rad=0.02
      expected%tWires%tw(1)%dispfile = trim(adjustl(" "))
      expected%tWires%tw(1)%dispfile_LeftEnd = trim(adjustl(" "))
      expected%tWires%tw(1)%dispfile_RightEnd = trim(adjustl(" "))
      expected%tWires%tw(1)%n_twc=10
      expected%tWires%tw(1)%n_twc_max=10
      allocate(expected%tWires%tw(1)%twc(10))
      expected%tWires%tw(1)%twc(1:10)%srcfile = 'None'
      expected%tWires%tw(1)%twc(1:10)%srctype = 'None'
      expected%tWires%tw(1)%twc(1:10)%i = 11
      expected%tWires%tw(1)%twc(1:10)%j = 11
      expected%tWires%tw(1)%twc(1:10)%k = [(i, i=7, 16)]
      expected%tWires%tw(1)%twc(1:10)%d = DIR_Z
      expected%tWires%tw(1)%twc(1:10)%nd = -1
      expected%tWires%tw(1)%twc(1)%nd  = 1
      expected%tWires%tw(1)%twc(6)%nd  = 2
      expected%tWires%tw(1)%twc(10)%nd = 3
      
      expected%tWires%tw(1)%twc(1:10)%tag = trim(adjustl("2"))   ! The polyline id is used as tag.
      
      expected%tWires%tw(1)%tl = MATERIAL_CONS
      expected%tWires%tw(1)%tr = MATERIAL_CONS
      
      expected%tWires%n_tw = 1
      expected%tWires%n_tw_max = 1
  
      expected%mtln%time_step = 30e-12
      expected%mtln%number_of_steps = 1000

      allocate(expected%mtln%cables(1))

      expected%mtln%cables(1)%name = "single_wire"
      allocate(expected%mtln%cables(1)%inductance_per_meter(1,1))
      allocate(expected%mtln%cables(1)%capacitance_per_meter(1,1))
      allocate(expected%mtln%cables(1)%resistance_per_meter(1,1))
      allocate(expected%mtln%cables(1)%conductance_per_meter(1,1))
      expected%mtln%cables(1)%inductance_per_meter  = 0.0
      expected%mtln%cables(1)%capacitance_per_meter = 0.0
      expected%mtln%cables(1)%resistance_per_meter  = 0.0
      expected%mtln%cables(1)%conductance_per_meter = 0.0
      allocate(expected%mtln%cables(1)%step_size(10))
      expected%mtln%cables(1)%step_size =  [(0.1, i = 1, 10)]
      allocate(expected%mtln%cables(1)%external_field_segments(10))
      do i = 1, 10
         expected%mtln%cables(1)%external_field_segments(i)%position = (/11,11,i+6/)
         expected%mtln%cables(1)%external_field_segments(i)%direction = DIRECTION_Z_POS
         expected%mtln%cables(1)%external_field_segments(i)%field => null()

      end do

      allocate(expected%mtln%cables(1)%transfer_impedance%poles(0))
      allocate(expected%mtln%cables(1)%transfer_impedance%residues(0))

      expected%mtln%cables(1)%parent_cable => null()
      expected%mtln%cables(1)%conductor_in_parent = 0
      expected%mtln%cables(1)%initial_connector => null()
      expected%mtln%cables(1)%end_connector => null()

      ! probes
      deallocate(expected%mtln%probes)
      allocate(expected%mtln%probes(1))
      expected%mtln%probes(1)%attached_to_cable => expected%mtln%cables(1)
      expected%mtln%probes(1)%index = 6
      expected%mtln%probes(1)%probe_type = PROBE_TYPE_CURRENT
      expected%mtln%probes(1)%probe_name = "mid_point"
      expected%mtln%probes(1)%probe_position = [11,11,12]

      ! networks
      deallocate(expected%mtln%networks)
      allocate(expected%mtln%networks(2))

      allocate(expected%mtln%networks(1)%connections(1))
      allocate(expected%mtln%networks(1)%connections(1)%nodes(1))
      expected%mtln%networks(1)%connections(1)%nodes(1)%conductor_in_cable = 1
      expected%mtln%networks(1)%connections(1)%nodes(1)%side = TERMINAL_NODE_SIDE_INI
      expected%mtln%networks(1)%connections(1)%nodes(1)%belongs_to_cable =>  expected%mtln%cables(1)
      expected%mtln%networks(1)%connections(1)%nodes(1)%termination%termination_type = TERMINATION_OPEN

      allocate(expected%mtln%networks(2)%connections(1))
      allocate(expected%mtln%networks(2)%connections(1)%nodes(1))
      expected%mtln%networks(2)%connections(1)%nodes(1)%conductor_in_cable = 1
      expected%mtln%networks(2)%connections(1)%nodes(1)%side = TERMINAL_NODE_SIDE_END
      expected%mtln%networks(2)%connections(1)%nodes(1)%belongs_to_cable =>  expected%mtln%cables(1)
      expected%mtln%networks(2)%connections(1)%nodes(1)%termination%termination_type = TERMINATION_OPEN


   end function
end function

