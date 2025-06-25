integer function test_read_holland1981_mtl() bind (C) result(err)
   use smbjson
   use smbjson_testingTools

   implicit none

   character(len=*),parameter :: filename = PATH_TO_TEST_DATA//INPUT_EXAMPLES//'holland1981_mtl.fdtd.json'
   type(Parseador) :: problem, expected
   type(parser_t) :: parser
   logical :: areSame
   err = 0

   expected = expectedProblemDescription()
   parser = parser_t(filename)
   problem = parser%readProblemDescription()
   call expect_eq(err, expected, problem)

contains
   function expectedProblemDescription() result (ex)
      type(Parseador) :: ex
      integer :: i

      call initializeProblemDescription(ex)

      ! Expected general info.
      ex%general%dt = 30e-12
      ex%general%nmax = 1000

      ! Excected media matrix.
      ex%matriz%totalX = 20
      ex%matriz%totalY = 20
      ex%matriz%totalZ = 22

      ! Expected grid.
      ex%despl%nX = 20
      ex%despl%nY = 20
      ex%despl%nZ = 22

      allocate(ex%despl%desX(20))
      allocate(ex%despl%desY(20))
      allocate(ex%despl%desZ(22))
      ex%despl%desX = 0.1
      ex%despl%desY = 0.1
      ex%despl%desZ = 0.1
      ex%despl%mx1 = 0
      ex%despl%mx2 = 20
      ex%despl%my1 = 0
      ex%despl%my2 = 20
      ex%despl%mz1 = 0
      ex%despl%mz2 = 22

      ! Expected boundaries.
      ex%front%tipoFrontera(:) = F_PML
      ex%front%propiedadesPML(:)%numCapas = 6
      ex%front%propiedadesPML(:)%orden = 2.0
      ex%front%propiedadesPML(:)%refl = 0.001

      ! Expected sources.
      allocate(ex%plnSrc%collection(1))
      ex%plnSrc%collection(1)%nombre_fichero = "holland.exc"
      ex%plnSrc%collection(1)%atributo = ""
      ex%plnSrc%collection(1)%coor1 = [1, 1, 1]
      ex%plnSrc%collection(1)%coor2 = [18, 18, 20]
      ex%plnSrc%collection(1)%theta = 1.5708
      ex%plnSrc%collection(1)%phi = 0.0
      ex%plnSrc%collection(1)%alpha = 0.0
      ex%plnSrc%collection(1)%beta = 0.0
      ex%plnSrc%collection(1)%isRC=.false.
      ex%plnSrc%collection(1)%nummodes=1
      ex%plnSrc%collection(1)%INCERTMAX=0.0
      ex%plnSrc%nc = 1
      ex%plnSrc%nC_max = 1

      ! ! Expected probes
      ! ! sonda
      ! ex%Sonda%length = 1
      ! ex%Sonda%length_max = 1
      ! allocate(ex%Sonda%collection(1))
      ! ex%Sonda%collection(1)%outputrequest = "mid_point"
      ! ex%Sonda%collection(1)%type1 = NP_T1_PLAIN
      ! ex%Sonda%collection(1)%type2 = NP_T2_TIME
      ! ex%Sonda%collection(1)%filename = ' '
      ! ex%Sonda%collection(1)%tstart = 0.0
      ! ex%Sonda%collection(1)%tstop = 0.0
      ! ex%Sonda%collection(1)%tstep = 0.0
      ! ex%Sonda%collection(1)%fstart = 0.0
      ! ex%Sonda%collection(1)%fstop = 0.0
      ! ex%Sonda%collection(1)%fstep = 0.0
      ! allocate(ex%Sonda%collection(1)%cordinates(1))
      ! ex%Sonda%collection(1)%len_cor = 1
      ! ex%Sonda%collection(1)%cordinates(1)%tag = "mid_point"
      ! ex%Sonda%collection(1)%cordinates(1)%Xi = 2 ! Coord id as tag.
      ! ex%Sonda%collection(1)%cordinates(1)%Yi = 0
      ! ex%Sonda%collection(1)%cordinates(1)%Zi = 0
      ! ex%Sonda%collection(1)%cordinates(1)%Or = NP_COR_WIRECURRENT
      
      
      ! ! Expected thin wires
      ! allocate(ex%tWires%tw(1))
      ! ex%tWires%tw(1)%rad=0.02
      ! ex%tWires%tw(1)%dispfile = trim(adjustl(" "))
      ! ex%tWires%tw(1)%dispfile_LeftEnd = trim(adjustl(" "))
      ! ex%tWires%tw(1)%dispfile_RightEnd = trim(adjustl(" "))
      ! ex%tWires%tw(1)%n_twc=10
      ! ex%tWires%tw(1)%n_twc_max=10
      ! allocate(ex%tWires%tw(1)%twc(10))
      ! ex%tWires%tw(1)%twc(1:10)%srcfile = 'None'
      ! ex%tWires%tw(1)%twc(1:10)%srctype = 'None'
      ! ex%tWires%tw(1)%twc(1:10)%i = 11
      ! ex%tWires%tw(1)%twc(1:10)%j = 11
      ! ex%tWires%tw(1)%twc(1:10)%k = [(i, i=7, 16)]
      ! ex%tWires%tw(1)%twc(1:10)%d = DIR_Z
      ! ex%tWires%tw(1)%twc(1:10)%nd = -1
      ! ex%tWires%tw(1)%twc(1)%nd  = 1
      ! ex%tWires%tw(1)%twc(6)%nd  = 2
      ! ex%tWires%tw(1)%twc(10)%nd = 3
      
      ! ex%tWires%tw(1)%twc(1:10)%tag = trim(adjustl("2"))   ! The polyline id is used as tag.
      
      ! ex%tWires%tw(1)%tl = MATERIAL_CONS
      ! ex%tWires%tw(1)%tr = MATERIAL_CONS
      
      ! ex%tWires%n_tw = 1
      ! ex%tWires%n_tw_max = 1
  

      ex%mtln%time_step = 30e-12
      ex%mtln%number_of_steps = 1000

      
      deallocate(ex%mtln%cables)
      allocate(ex%mtln%cables(1))
      allocate(unshielded_multiwire_t :: ex%mtln%cables(1)%ptr)

      ex%mtln%cables(1)%ptr%name = "single_unshielded_multiwire"
      call initializeCablePULParameters(ex%mtln%cables(1)%ptr)
      allocate(ex%mtln%cables(1)%ptr%step_size(10))
      ex%mtln%cables(1)%ptr%step_size =  [(0.1, i = 1, 10)]
      allocate(ex%mtln%cables(1)%ptr%segments(10))
      do i = 1, 10
         ex%mtln%cables(1)%ptr%segments(i)%x = 11
         ex%mtln%cables(1)%ptr%segments(i)%y = 11
         ex%mtln%cables(1)%ptr%segments(i)%z = i+6
         ex%mtln%cables(1)%ptr%segments(i)%orientation = DIRECTION_Z_POS
      end do

      ex%mtln%cables(1)%ptr%initial_connector => null()
      ex%mtln%cables(1)%ptr%end_connector => null()

      ! probes
      deallocate(ex%mtln%probes)
      allocate(ex%mtln%probes(1))
      ex%mtln%probes(1)%attached_to_cable => ex%mtln%cables(1)%ptr
      ex%mtln%probes(1)%index = 6
      ex%mtln%probes(1)%probe_type = PROBE_TYPE_CURRENT
      ex%mtln%probes(1)%probe_name = "mid_point"
      ex%mtln%probes(1)%probe_position = [11,11,12]

      ! networks
      deallocate(ex%mtln%networks)
      allocate(ex%mtln%networks(2))

      allocate(ex%mtln%networks(1)%connections(1))
      allocate(ex%mtln%networks(1)%connections(1)%nodes(1))
      ex%mtln%networks(1)%connections(1)%nodes(1)%conductor_in_cable = 1
      ex%mtln%networks(1)%connections(1)%nodes(1)%side = TERMINAL_NODE_SIDE_INI
      ex%mtln%networks(1)%connections(1)%nodes(1)%belongs_to_cable =>  ex%mtln%cables(1)%ptr
      ex%mtln%networks(1)%connections(1)%nodes(1)%termination%termination_type = TERMINATION_OPEN

      allocate(ex%mtln%networks(2)%connections(1))
      allocate(ex%mtln%networks(2)%connections(1)%nodes(1))
      ex%mtln%networks(2)%connections(1)%nodes(1)%conductor_in_cable = 1
      ex%mtln%networks(2)%connections(1)%nodes(1)%side = TERMINAL_NODE_SIDE_END
      ex%mtln%networks(2)%connections(1)%nodes(1)%belongs_to_cable =>  ex%mtln%cables(1)%ptr
      ex%mtln%networks(2)%connections(1)%nodes(1)%termination%termination_type = TERMINATION_OPEN


   end function
end function

