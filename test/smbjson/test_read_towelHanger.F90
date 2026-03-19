integer function test_read_towelhanger() bind (C) result(err)
   use smbjson_m
   use smbjson_testingTools

   implicit none

   character(len=*),parameter :: filename = PATH_TO_TEST_DATA//INPUT_EXAMPLES//'towelHanger.fdtd.json'
   type(Parseador_t) :: problem, expected
   type(parser_t) :: parser
   logical :: areSame
   err = 0

   expected = expectedProblemDescription()
   parser = parser_t(filename)
   problem = parser%readProblemDescription()
   call expect_eq(err, expected, problem)

contains
   function expectedProblemDescription() result (expected)
      type(Parseador_t) :: expected

      integer :: i

      call initializeProblemDescription(expected)

      ! Expected general info.
      expected%general%dt = 1e-12
      expected%general%nmax = 2000

      ! Excected media matrix.
      expected%matriz%totalX = 60
      expected%matriz%totalY = 60
      expected%matriz%totalZ = 60

      ! Expected grid.
      expected%despl%nX = 60
      expected%despl%nY = 60
      expected%despl%nZ = 60

      allocate(expected%despl%desX(60))
      allocate(expected%despl%desY(60))
      allocate(expected%despl%desZ(60))
      expected%despl%desX = 0.01
      expected%despl%desY = 0.01
      expected%despl%desZ = 0.01
      expected%despl%mx1 = 0
      expected%despl%mx2 = 60
      expected%despl%my1 = 0
      expected%despl%my2 = 60
      expected%despl%mz1 = 0
      expected%despl%mz2 = 60

      ! Expected boundaries.
      expected%front%tipoFrontera(:) = F_PML
      expected%front%propiedadesPML(:)%numCapas = 6
      expected%front%propiedadesPML(:)%orden = 2.0
      expected%front%propiedadesPML(:)%refl = 0.001

      ! Expected PEC regions.
      expected%pecRegs%nLins = 0
      expected%pecRegs%nLins_max = 0
      expected%pecRegs%nSurfs = 1
      expected%pecRegs%nSurfs_max = 1
      expected%pecRegs%nVols = 0
      expected%pecRegs%nVols_max = 0
      allocate(expected%pecRegs%Surfs(1))
      expected%pecRegs%Surfs(1)%Xi = 25
      expected%pecRegs%Surfs(1)%Xe = 44
      expected%pecRegs%Surfs(1)%yi = 20
      expected%pecRegs%Surfs(1)%ye = 29
      expected%pecRegs%Surfs(1)%zi = 30
      expected%pecRegs%Surfs(1)%ze = 30
      expected%pecRegs%Surfs(1)%Xtrancos = 1
      expected%pecRegs%Surfs(1)%Ytrancos = 1
      expected%pecRegs%Surfs(1)%Ztrancos = 1
      expected%pecRegs%Surfs(1)%Or = 3
      expected%pecRegs%Surfs(1)%tag =  "copper@ground_plane"

      ! expected mtln bundles
      expected%mtln%time_step = 1e-12
      expected%mtln%number_of_steps = 2000
      deallocate(expected%mtln%cables)
      allocate(expected%mtln%cables(1))
      allocate(unshielded_multiwire_t :: expected%mtln%cables(1)%ptr)

      expected%mtln%cables(1)%ptr%name = "wire"
      call initializeCablePulParameters(expected%mtln%cables(1)%ptr)
      allocate(expected%mtln%cables(1)%ptr%step_size(20))
      expected%mtln%cables(1)%ptr%step_size =  [(0.01, i = 1, 20)]
      allocate(expected%mtln%cables(1)%ptr%segments(20))
      
      do i = 1,2
         expected%mtln%cables(1)%ptr%segments(i)%x = 27
         expected%mtln%cables(1)%ptr%segments(i)%y = 25
         expected%mtln%cables(1)%ptr%segments(i)%z = 29+i
         expected%mtln%cables(1)%ptr%segments(i)%orientation = DIRECTION_Z_POS
      end do

      do i = 3, 10
         expected%mtln%cables(1)%ptr%segments(i)%x = 24+i
         expected%mtln%cables(1)%ptr%segments(i)%y = 25
         expected%mtln%cables(1)%ptr%segments(i)%z = 32
         expected%mtln%cables(1)%ptr%segments(i)%orientation = DIRECTION_X_POS
      end do

      do i = 1,8
         expected%mtln%cables(1)%ptr%segments(10+i)%x = 34+i
         expected%mtln%cables(1)%ptr%segments(10+i)%y = 25
         expected%mtln%cables(1)%ptr%segments(10+i)%z = 32
         expected%mtln%cables(1)%ptr%segments(10+i)%orientation = DIRECTION_X_POS
      end do

      do i = 9,10
         expected%mtln%cables(1)%ptr%segments(10+i)%x = 43
         expected%mtln%cables(1)%ptr%segments(10+i)%y = 25
         expected%mtln%cables(1)%ptr%segments(10+i)%z = 40-i
         expected%mtln%cables(1)%ptr%segments(10+i)%orientation = DIRECTION_Z_NEG
      end do

      expected%mtln%cables(1)%ptr%initial_connector => null()
      expected%mtln%cables(1)%ptr%end_connector => null()



      ! probes
      deallocate(expected%mtln%probes)
      allocate(expected%mtln%probes(3))
      expected%mtln%probes(1)%attached_to_cable => expected%mtln%cables(1)%ptr
      expected%mtln%probes(1)%index = 1
      expected%mtln%probes(1)%probe_type = PROBE_TYPE_CURRENT
      expected%mtln%probes(1)%probe_name = "wire_start"
      expected%mtln%probes(1)%probe_position = [27,25,30]
      
      expected%mtln%probes(2)%attached_to_cable => expected%mtln%cables(1)%ptr
      expected%mtln%probes(2)%index = 21
      expected%mtln%probes(2)%probe_type = PROBE_TYPE_CURRENT
      expected%mtln%probes(2)%probe_name = "wire_end"
      expected%mtln%probes(2)%probe_position = [43, 25, 30]

      expected%mtln%probes(3)%attached_to_cable => expected%mtln%cables(1)%ptr
      expected%mtln%probes(3)%index = 11
      expected%mtln%probes(3)%probe_type = PROBE_TYPE_CURRENT
      expected%mtln%probes(3)%probe_name = "wire_mid"
      expected%mtln%probes(3)%probe_position = [35, 25, 32]

      ! networks
      deallocate(expected%mtln%networks)
      allocate(expected%mtln%networks(2))

      allocate(expected%mtln%networks(1)%connections(1))
      allocate(expected%mtln%networks(1)%connections(1)%nodes(1))
      expected%mtln%networks(1)%connections(1)%nodes(1)%conductor_in_cable = 1
      expected%mtln%networks(1)%connections(1)%nodes(1)%side = TERMINAL_NODE_SIDE_INI
      expected%mtln%networks(1)%connections(1)%nodes(1)%belongs_to_cable =>  expected%mtln%cables(1)%ptr
      expected%mtln%networks(1)%connections(1)%nodes(1)%termination%termination_type = TERMINATION_SERIES
      expected%mtln%networks(1)%connections(1)%nodes(1)%termination%resistance = 50.0
      expected%mtln%networks(1)%connections(1)%nodes(1)%termination%source%path_to_excitation = "towelHanger.exc"
      expected%mtln%networks(1)%connections(1)%nodes(1)%termination%source%source_type = SOURCE_TYPE_VOLTAGE

      allocate(expected%mtln%networks(2)%connections(1))
      allocate(expected%mtln%networks(2)%connections(1)%nodes(1))
      expected%mtln%networks(2)%connections(1)%nodes(1)%conductor_in_cable = 1
      expected%mtln%networks(2)%connections(1)%nodes(1)%side = TERMINAL_NODE_SIDE_END
      expected%mtln%networks(2)%connections(1)%nodes(1)%belongs_to_cable =>  expected%mtln%cables(1)%ptr
      expected%mtln%networks(2)%connections(1)%nodes(1)%termination%termination_type = TERMINATION_SHORT


   end function
end function

