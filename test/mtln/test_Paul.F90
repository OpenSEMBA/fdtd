! integer function test_termination_resistive() bind(C) result(error_cnt)
!     use mtln_solver_mod
!     use mtln_testingTools_mod
!     use preprocess_mod
!     implicit none

!     character(len=*), parameter :: square_excitation = PATH_TO_TEST_DATA//'excitations/termination_resistive_pulse.exc'

!     type(cable_t), target :: cable
!     type(terminal_node_t) :: node_left, node_right
!     type(terminal_network_t) :: nw_left, nw_right
!     type(parsed_mtln_t) :: parsed
!     type(terminal_connection_t) :: connection_left, connection_right
!     type(terminal_network_t) :: network_left, network_right

!     type(parsed_probe_t) :: probe_v, probe_i, probe_v_r

!     type(preprocess_t) :: pre
!     type(mtln_t) :: solver
!     real,dimension(1,1) :: lpul = 0.25e-6
!     real,dimension(1,1) :: gpul = 0.0
!     real,dimension(1,1) :: rpul = 0.0
!     real,dimension(1,1) :: cpul = 100.0e-12

!     error_cnt = 0

!     cable%name = "wire0"
!     cable%inductance_per_meter = lpul
!     cable%conductance_per_meter = gpul
!     cable%resistance_per_meter = rpul
!     cable%capacitance_per_meter = cpul
!     block
!         integer :: i
!         cable%step_size = [(4, i = 1, 50)]
!     end block

!     parsed%time_step = 1.0e-8
!     parsed%number_of_steps = 40e-6/parsed%time_step


!     node_left%belongs_to_cable => cable
!     node_left%conductor_in_cable = 1
!     node_left%side = TERMINAL_NODE_SIDE_INI

!     node_left%termination = termination_with_source_t(path_to_excitation=square_excitation, & 
!                                                 termination_type = TERMINATION_SERIES, &
!                                                 resistance = 150, &
!                                                 inductance = 0.0, &
!                                                 capacitance = 1e22)

!     connection_left%nodes = [node_left]

!     node_right%belongs_to_cable => cable
!     node_right%conductor_in_cable = 1
!     node_right%side = TERMINAL_NODE_SIDE_END
!     node_right%termination = termination_t(termination_type = TERMINATION_SERIES, &
!                                            resistance = 50, &
!                                            inductance = 0.0, &
!                                            capacitance = 1e22)

!     connection_right%nodes = [node_right]

!     network_left%connections = [connection_left]
!     network_right%connections = [connection_right]

!     probe_v%attached_to_cable => cable
!     probe_v%index = 1 
!     probe_v%probe_type = PROBE_TYPE_VOLTAGE

!     probe_i%attached_to_cable => cable
!     probe_i%index = 50
!     probe_i%probe_type = PROBE_TYPE_CURRENT

!     probe_v_r%attached_to_cable => cable
!     probe_v_r%index = 50
!     probe_v_r%probe_type = PROBE_TYPE_VOLTAGE


!     ! parsed%networks = [network_left]
!     parsed%networks = [network_left, network_right]
!     parsed%cables = [cable]
!     parsed%probes = [probe_v, probe_i, probe_v_r]    ! pre = preprocess(parsed)

!     solver = mtlnCtor(parsed)
!     call solver%runUntil(parsed%time_step * parsed%number_of_steps)
!     write(*,*) error_cnt
!     ! p = Parser(file)
!     ! p.run(finalTime = 18e-6)
!     block
!         integer :: i
!         open(unit = 1, file =  './probes/probes_R_termination.txt')
!         do i = 1, size(solver%bundles(1)%probes(1)%t)
!             write(1,*) solver%bundles(1)%probes(1)%t(i)," ", &
!                     solver%bundles(1)%probes(1)%val(i,1) ," ", &
!                     solver%bundles(1)%probes(2)%val(i,1) ," ", &
!                     solver%bundles(1)%probes(3)%val(i,1)
!         end do

!     end block

! end function

! integer function test_termination_resistive_inductive() bind(C) result(error_cnt)
!     use mtln_solver_mod
!     use mtln_testingTools_mod
!     use preprocess_mod
!     implicit none

!     character(len=*), parameter :: square_excitation = PATH_TO_TEST_DATA//'excitations/termination_resistive_pulse.exc'

!     type(cable_t), target :: cable
!     type(terminal_node_t) :: node_left, node_right
!     type(terminal_network_t) :: nw_left, nw_right
!     type(parsed_mtln_t) :: parsed
!     type(terminal_connection_t) :: connection_left, connection_right
!     type(terminal_network_t) :: network_left, network_right

!     type(parsed_probe_t) :: probe_v, probe_i, probe_v_R

!     type(preprocess_t) :: pre
!     type(mtln_t) :: solver
!     real,dimension(1,1) :: lpul = 0.25e-6
!     real,dimension(1,1) :: gpul = 0.0
!     real,dimension(1,1) :: rpul = 0.0
!     real,dimension(1,1) :: cpul = 100.0e-12

!     error_cnt = 0

!     cable%name = "wire0"
!     cable%inductance_per_meter = lpul
!     cable%conductance_per_meter = gpul
!     cable%resistance_per_meter = rpul
!     cable%capacitance_per_meter = cpul
!     block
!         integer :: i
!         cable%step_size = [(4, i = 1, 50)]
!     end block

!     parsed%time_step = 1.5e-8
!     parsed%number_of_steps = 40e-6/parsed%time_step


!     node_left%belongs_to_cable => cable
!     node_left%conductor_in_cable = 1
!     node_left%side = TERMINAL_NODE_SIDE_INI

!     node_left%termination = termination_with_source_t(path_to_excitation=square_excitation, & 
!                                                 termination_type = TERMINATION_SERIES, &
!                                                 resistance = 150, &
!                                                 inductance = 0.0, &
!                                                 capacitance = 1e22)

!     connection_left%nodes = [node_left]

!     node_right%belongs_to_cable => cable
!     node_right%conductor_in_cable = 1
!     node_right%side = TERMINAL_NODE_SIDE_END
!     node_right%termination = termination_t(termination_type = TERMINATION_SERIES, &
!                                             resistance = 50, &
!                                             inductance = 1e-5, &
!                                             capacitance = 1e22)

!     connection_right%nodes = [node_right]

!     network_left%connections = [connection_left]
!     network_right%connections = [connection_right]

!     probe_v%attached_to_cable => cable
!     probe_v%index = 1 
!     probe_v%probe_type = PROBE_TYPE_VOLTAGE

!     probe_i%attached_to_cable => cable
!     probe_i%index = 50
!     probe_i%probe_type = PROBE_TYPE_CURRENT

!     probe_v_R%attached_to_cable => cable
!     probe_v_R%index = 50
!     probe_v_R%probe_type = PROBE_TYPE_VOLTAGE


!     ! parsed%networks = [network_left]
!     parsed%networks = [network_left, network_right]
!     parsed%cables = [cable]
!     parsed%probes = [probe_v, probe_i, probe_v_R]    ! pre = preprocess(parsed)

!     solver = mtlnCtor(parsed)
!     call solver%runUntil(parsed%time_step * parsed%number_of_steps)
!     write(*,*) error_cnt
!     ! p = Parser(file)
!     ! p.run(finalTime = 18e-6)
!     block
!         integer :: i
!         open(unit = 1, file =  './testData/outputs/spice/RL_termination.txt')
!         do i = 1, size(solver%bundles(1)%probes(1)%t)
!             write(1,*) solver%bundles(1)%probes(1)%t(i)," ", &
!                     solver%bundles(1)%probes(1)%val(i,1) ," ", &
!                     solver%bundles(1)%probes(2)%val(i,1) ," ", &
!                     solver%bundles(1)%probes(3)%val(i,1)
!         end do

!     end block

! end function

! integer function test_termination_resistive_capacitive_parallel() bind(C) result(error_cnt)
!     use mtln_solver_mod
!     use mtln_testingTools_mod
!     use preprocess_mod
!     implicit none

!     character(len=*), parameter :: square_excitation = &
!         PATH_TO_TEST_DATA//'excitations/termination_resistive_pulse.exc'

!     type(cable_t), target :: cable
!     type(terminal_node_t) :: node_left, node_right
!     type(terminal_network_t) :: nw_left, nw_right
!     type(parsed_mtln_t) :: parsed
!     type(terminal_connection_t) :: connection_left, connection_right
!     type(terminal_network_t) :: network_left, network_right

!     type(parsed_probe_t) :: probe_v, probe_i, probe_v_r

!     type(preprocess_t) :: pre
!     type(mtln_t) :: solver
!     real,dimension(1,1) :: lpul = 0.25e-6
!     real,dimension(1,1) :: gpul = 0.0
!     real,dimension(1,1) :: rpul = 0.0
!     real,dimension(1,1) :: cpul = 100.0e-12

!     error_cnt = 0

!     cable%name = "wire0"
!     cable%inductance_per_meter = lpul
!     cable%conductance_per_meter = gpul
!     cable%resistance_per_meter = rpul
!     cable%capacitance_per_meter = cpul
!     block
!         integer :: i
!         cable%step_size = [(4, i = 1, 50)]
!     end block

!     parsed%time_step = 2.0e-8
!     parsed%number_of_steps = 40e-6/parsed%time_step


!     node_left%belongs_to_cable => cable
!     node_left%conductor_in_cable = 1
!     node_left%side = TERMINAL_NODE_SIDE_INI

!     node_left%termination = termination_with_source_t(path_to_excitation=square_excitation, & 
!                                                 termination_type = TERMINATION_SERIES, &
!                                                 resistance = 150, &
!                                                 inductance = 0.0, &
!                                                 capacitance = 1e22)

!     connection_left%nodes = [node_left]

!     node_right%belongs_to_cable => cable
!     node_right%conductor_in_cable = 1
!     node_right%side = TERMINAL_NODE_SIDE_END
!     node_right%termination = termination_t(termination_type = TERMINATION_RLsCp, &
!                                                   resistance = 0.1, &
!                                                   inductance = 0.0, &
!                                                   capacitance = 100e-12)

!     connection_right%nodes = [node_right]

!     network_left%connections = [connection_left]
!     network_right%connections = [connection_right]

!     probe_v%attached_to_cable => cable
!     probe_v%index = 1 
!     probe_v%probe_type = PROBE_TYPE_VOLTAGE

!     probe_i%attached_to_cable => cable
!     probe_i%index = 50
!     probe_i%probe_type = PROBE_TYPE_CURRENT

!     probe_v_r%attached_to_cable => cable
!     probe_v_r%index = 50
!     probe_v_r%probe_type = PROBE_TYPE_VOLTAGE


!     ! parsed%networks = [network_left]
!     parsed%networks = [network_left, network_right]
!     parsed%cables = [cable]
!     parsed%probes = [probe_v, probe_i, probe_v_r]    ! pre = preprocess(parsed)

!     solver = mtlnCtor(parsed)
!     call solver%runUntil(parsed%time_step * parsed%number_of_steps)
!     write(*,*) error_cnt
!     ! p = Parser(file)
!     ! p.run(finalTime = 18e-6)
!     block
!         integer :: i
!         open(unit = 1, file =  './testData/outputs/spice/RCp_termination.txt')
!         do i = 1, size(solver%bundles(1)%probes(1)%t)
!             write(1,*) solver%bundles(1)%probes(1)%t(i)," ", &
!                     solver%bundles(1)%probes(1)%val(i,1) ," ", &
!                     solver%bundles(1)%probes(2)%val(i,1) ," ", &
!                     solver%bundles(1)%probes(3)%val(i,1)
!         end do

!     end block

! end function

! integer function test_termination_rls_cp() bind(C) result(error_cnt)
!     use mtln_solver_mod
!     use mtln_testingTools_mod
!     use preprocess_mod
!     implicit none

!     character(len=*), parameter :: square_excitation = &
!         PATH_TO_TEST_DATA//'excitations/termination_resistive_pulse.exc'

!     type(cable_t), target :: cable
!     type(terminal_node_t) :: node_left, node_right
!     type(terminal_network_t) :: nw_left, nw_right
!     type(parsed_mtln_t) :: parsed
!     type(terminal_connection_t) :: connection_left, connection_right
!     type(terminal_network_t) :: network_left, network_right

!     type(parsed_probe_t) :: probe_v, probe_i, probe_v_r

!     type(preprocess_t) :: pre
!     type(mtln_t) :: solver
!     real,dimension(1,1) :: lpul = 0.25e-6
!     real,dimension(1,1) :: gpul = 0.0
!     real,dimension(1,1) :: rpul = 0.0
!     real,dimension(1,1) :: cpul = 100.0e-12

!     error_cnt = 0

!     cable%name = "wire0"
!     cable%inductance_per_meter = lpul
!     cable%conductance_per_meter = gpul
!     cable%resistance_per_meter = rpul
!     cable%capacitance_per_meter = cpul
!     block
!         integer :: i
!         cable%step_size = [(4, i = 1, 50)]
!     end block

!     parsed%time_step = 2.0e-8
!     parsed%number_of_steps = 40e-6/parsed%time_step


!     node_left%belongs_to_cable => cable
!     node_left%conductor_in_cable = 1
!     node_left%side = TERMINAL_NODE_SIDE_INI

!     node_left%termination = termination_with_source_t(path_to_excitation=square_excitation, & 
!                                                 termination_type = TERMINATION_SERIES, &
!                                                 resistance = 150, &
!                                                 inductance = 0.0, &
!                                                 capacitance = 1e22)

!     connection_left%nodes = [node_left]

!     node_right%belongs_to_cable => cable
!     node_right%conductor_in_cable = 1
!     node_right%side = TERMINAL_NODE_SIDE_END
!     node_right%termination = termination_t(termination_type = TERMINATION_RLsCP, &
!                                                   resistance = 10, &
!                                                   inductance = 10e-6, &
!                                                   capacitance = 100e-12)

!     connection_right%nodes = [node_right]

!     network_left%connections = [connection_left]
!     network_right%connections = [connection_right]

!     probe_v%attached_to_cable => cable
!     probe_v%index = 1 
!     probe_v%probe_type = PROBE_TYPE_VOLTAGE

!     probe_i%attached_to_cable => cable
!     probe_i%index = 50
!     probe_i%probe_type = PROBE_TYPE_CURRENT

!     probe_v_r%attached_to_cable => cable
!     probe_v_r%index = 50
!     probe_v_r%probe_type = PROBE_TYPE_VOLTAGE


!     ! parsed%networks = [network_left]
!     parsed%networks = [network_left, network_right]
!     parsed%cables = [cable]
!     parsed%probes = [probe_v, probe_i, probe_v_r]    ! pre = preprocess(parsed)

!     solver = mtlnCtor(parsed)
!     call solver%runUntil(parsed%time_step * parsed%number_of_steps)
!     write(*,*) error_cnt
!     ! p = Parser(file)
!     ! p.run(finalTime = 18e-6)
!     block
!         integer :: i
!         open(unit = 1, file =  './testData/outputs/spice/RLsCp_termination.txt')
!         do i = 1, size(solver%bundles(1)%probes(1)%t)
!             write(1,*) solver%bundles(1)%probes(1)%t(i)," ", &
!                     solver%bundles(1)%probes(1)%val(i,1) ," ", &
!                     solver%bundles(1)%probes(2)%val(i,1) ," ", &
!                     solver%bundles(1)%probes(3)%val(i,1)
!         end do

!     end block

! end function

! integer function test_termination_rls_cp_ns() bind(C) result(error_cnt)
!     use mtln_solver_mod
!     use mtln_testingTools_mod
!     use preprocess_mod
!     implicit none

!     character(len=*), parameter :: pulse_excitation = PATH_TO_TEST_DATA//'excitations/2_conductor_line_paul_9_6_gauss.exc'

!     type(cable_t), target :: cable
!     type(terminal_node_t) :: node_left, node_right
!     type(terminal_network_t) :: nw_left, nw_right
!     type(parsed_mtln_t) :: parsed
!     type(terminal_connection_t) :: connection_left, connection_right
!     type(terminal_network_t) :: network_left, network_right

!     type(parsed_probe_t) :: probe_v, probe_i, probe_v_r

!     type(preprocess_t) :: pre
!     type(mtln_t) :: solver
!     real,dimension(1,1) :: lpul = 0.25e-6
!     real,dimension(1,1) :: gpul = 0.0
!     real,dimension(1,1) :: rpul = 0.0
!     real,dimension(1,1) :: cpul = 100.0e-12

!     real :: final_time

!     error_cnt = 0

!     cable%name = "wire0"
!     cable%inductance_per_meter = lpul
!     cable%conductance_per_meter = gpul
!     cable%resistance_per_meter = rpul
!     cable%capacitance_per_meter = cpul
!     block
!         integer :: i
!         real :: length = 0.5
!         ! integer :: ndz = 1000
!         integer :: ndz = 795
!         cable%step_size = [(length/ndz, i = 1, ndz)]
!     end block

!     final_time = 40e-9
!     parsed%number_of_steps = 9410
!     parsed%time_step = final_time/parsed%number_of_steps


!     node_left%belongs_to_cable => cable
!     node_left%conductor_in_cable = 1
!     node_left%side = TERMINAL_NODE_SIDE_INI

!     node_left%termination = termination_with_source_t(path_to_excitation=pulse_excitation, & 
!                                                 termination_type = TERMINATION_SERIES, &
!                                                 resistance = 150, &
!                                                 inductance = 0.0, &
!                                                 capacitance = 1e22)

!     connection_left%nodes = [node_left]

!     node_right%belongs_to_cable => cable
!     node_right%conductor_in_cable = 1
!     node_right%side = TERMINAL_NODE_SIDE_END
!     node_right%termination = termination_t(termination_type = TERMINATION_RLsCP, &
!                                                   resistance = 10, &
!                                                   inductance = 10e-6, &
!                                                   capacitance = 100e-12)

!     connection_right%nodes = [node_right]

!     network_left%connections = [connection_left]
!     network_right%connections = [connection_right]

!     probe_v%attached_to_cable => cable
!     probe_v%index = 1 
!     probe_v%probe_type = PROBE_TYPE_VOLTAGE

!     probe_i%attached_to_cable => cable
!     probe_i%index = 50
!     probe_i%probe_type = PROBE_TYPE_CURRENT

!     probe_v_r%attached_to_cable => cable
!     probe_v_r%index = 50
!     probe_v_r%probe_type = PROBE_TYPE_VOLTAGE


!     ! parsed%networks = [network_left]
!     parsed%networks = [network_left, network_right]
!     parsed%cables = [cable]
!     parsed%probes = [probe_v, probe_i, probe_v_r]    ! pre = preprocess(parsed)

!     solver = mtlnCtor(parsed)
!     call solver%runUntil(parsed%time_step * parsed%number_of_steps)
!     write(*,*) error_cnt
!     ! p = Parser(file)
!     ! p.run(finalTime = 18e-6)
!     block
!         integer :: i
!         open(unit = 1, file =  './testData/outputs/spice/RLsCp_termination_ns.txt')
!         do i = 1, size(solver%bundles(1)%probes(1)%t)
!             write(1,*) solver%bundles(1)%probes(1)%t(i)," ", &
!                     solver%bundles(1)%probes(1)%val(i,1) ," ", &
!                     solver%bundles(1)%probes(2)%val(i,1) ," ", &
!                     solver%bundles(1)%probes(3)%val(i,1)
!         end do

!     end block

! end function

! integer function test_termination_rcp() bind(C) result(error_cnt)
!     use mtln_solver_mod
!     use mtln_testingTools_mod
!     use preprocess_mod
!     implicit none

!     character(len=*), parameter :: square_excitation = &
!         PATH_TO_TEST_DATA//'excitations/termination_resistive_pulse.exc'

!     type(cable_t), target :: cable
!     type(terminal_node_t) :: node_left, node_right
!     type(terminal_network_t) :: nw_left, nw_right
!     type(parsed_mtln_t) :: parsed
!     type(terminal_connection_t) :: connection_left, connection_right
!     type(terminal_network_t) :: network_left, network_right

!     type(parsed_probe_t) :: probe_v, probe_i, probe_v_r

!     type(preprocess_t) :: pre
!     type(mtln_t) :: solver
!     real,dimension(1,1) :: lpul = 0.25e-6
!     real,dimension(1,1) :: gpul = 0.0
!     real,dimension(1,1) :: rpul = 0.0
!     real,dimension(1,1) :: cpul = 100.0e-12

!     error_cnt = 0

!     cable%name = "wire0"
!     cable%inductance_per_meter = lpul
!     cable%conductance_per_meter = gpul
!     cable%resistance_per_meter = rpul
!     cable%capacitance_per_meter = cpul
!     block
!         integer :: i
!         cable%step_size = [(4, i = 1, 50)]
!     end block

!     parsed%time_step = 2.0e-8
!     parsed%number_of_steps = 40e-6/parsed%time_step


!     node_left%belongs_to_cable => cable
!     node_left%conductor_in_cable = 1
!     node_left%side = TERMINAL_NODE_SIDE_INI

!     node_left%termination = termination_with_source_t(path_to_excitation=square_excitation, & 
!                                                 termination_type = TERMINATION_SERIES, &
!                                                 resistance = 150, &
!                                                 inductance = 0.0, &
!                                                 capacitance = 1e22)

!     connection_left%nodes = [node_left]

!     node_right%belongs_to_cable => cable
!     node_right%conductor_in_cable = 1
!     node_right%side = TERMINAL_NODE_SIDE_END
!     node_right%termination = termination_t(termination_type = TERMINATION_RLsCp, &
!                                                   resistance = 10, &
!                                                   inductance = 0.0, &
!                                                   capacitance = 100e-12)

!     connection_right%nodes = [node_right]

!     network_left%connections = [connection_left]
!     network_right%connections = [connection_right]

!     probe_v%attached_to_cable => cable
!     probe_v%index = 1 
!     probe_v%probe_type = PROBE_TYPE_VOLTAGE

!     probe_i%attached_to_cable => cable
!     probe_i%index = 50
!     probe_i%probe_type = PROBE_TYPE_CURRENT

!     probe_v_r%attached_to_cable => cable
!     probe_v_r%index = 50
!     probe_v_r%probe_type = PROBE_TYPE_VOLTAGE


!     ! parsed%networks = [network_left]
!     parsed%networks = [network_left, network_right]
!     parsed%cables = [cable]
!     parsed%probes = [probe_v, probe_i, probe_v_r]    ! pre = preprocess(parsed)

!     solver = mtlnCtor(parsed)
!     call solver%runUntil(parsed%time_step * parsed%number_of_steps)
!     write(*,*) error_cnt
!     ! p = Parser(file)
!     ! p.run(finalTime = 18e-6)
!     block
!         integer :: i
!         open(unit = 1, file =  './testData/outputs/spice/RCp_termination.txt')
!         do i = 1, size(solver%bundles(1)%probes(1)%t)
!             write(1,*) solver%bundles(1)%probes(1)%t(i)," ", &
!                     solver%bundles(1)%probes(1)%val(i,1) ," ", &
!                     solver%bundles(1)%probes(2)%val(i,1) ," ", &
!                     solver%bundles(1)%probes(3)%val(i,1)
!         end do

!     end block

! end function

! integer function test_termination_resistive_capacitive() bind(C) result(error_cnt)
!     use mtln_solver_mod
!     use mtln_testingTools_mod
!     use preprocess_mod
!     implicit none

!     ! character(len=*), parameter :: square_excitation = PATH_TO_TEST_DATA//'excitations/termination_resistive_pulse.exc'
!     character(len=*), parameter :: pulse_excitation = PATH_TO_TEST_DATA//'excitations/5u_1u_gauss.exc'

!     type(cable_t), target :: cable
!     type(terminal_node_t) :: node_left, node_right
!     type(terminal_network_t) :: nw_left, nw_right
!     type(parsed_mtln_t) :: parsed
!     type(terminal_connection_t) :: connection_left, connection_right
!     type(terminal_network_t) :: network_left, network_right

!     type(parsed_probe_t) :: probe_v, probe_i, probe_v_r

!     type(preprocess_t) :: pre
!     type(mtln_t) :: solver
!     real,dimension(1,1) :: lpul = 0.25e-6
!     real,dimension(1,1) :: gpul = 0.0
!     real,dimension(1,1) :: rpul = 0.0
!     real,dimension(1,1) :: cpul = 100.0e-12

!     error_cnt = 0

!     cable%name = "wire0"
!     cable%inductance_per_meter = lpul
!     cable%conductance_per_meter = gpul
!     cable%resistance_per_meter = rpul
!     cable%capacitance_per_meter = cpul
!     block
!         integer :: i
!         cable%step_size = [(4, i = 1, 50)]
!     end block

!     parsed%time_step = 2.0e-8
!     parsed%number_of_steps = 40e-6/parsed%time_step


!     node_left%belongs_to_cable => cable
!     node_left%conductor_in_cable = 1
!     node_left%side = TERMINAL_NODE_SIDE_INI

!     node_left%termination = termination_with_source_t(path_to_excitation=pulse_excitation, & 
!                                                 termination_type = TERMINATION_SERIES, &
!                                                 resistance = 150, &
!                                                 inductance = 0.0, &
!                                                 capacitance = 1e22)

!     connection_left%nodes = [node_left]

!     node_right%belongs_to_cable => cable
!     node_right%conductor_in_cable = 1
!     node_right%side = TERMINAL_NODE_SIDE_END
!     node_right%termination = termination_t(termination_type = TERMINATION_SERIES, &
!                                                   resistance = 50, &
!                                                   inductance = 0.0, &
!                                                   capacitance = 100e-12)

!     connection_right%nodes = [node_right]

!     network_left%connections = [connection_left]
!     network_right%connections = [connection_right]

!     probe_v%attached_to_cable => cable
!     probe_v%index = 1 
!     probe_v%probe_type = PROBE_TYPE_VOLTAGE

!     probe_i%attached_to_cable => cable
!     probe_i%index = 50
!     probe_i%probe_type = PROBE_TYPE_CURRENT

!     probe_v_r%attached_to_cable => cable
!     probe_v_r%index = 50
!     probe_v_r%probe_type = PROBE_TYPE_VOLTAGE


!     ! parsed%networks = [network_left]
!     parsed%networks = [network_left, network_right]
!     parsed%cables = [cable]
!     parsed%probes = [probe_v, probe_i, probe_v_r]    ! pre = preprocess(parsed)

!     solver = mtlnCtor(parsed)
!     call solver%runUntil(parsed%time_step * parsed%number_of_steps)
!     write(*,*) error_cnt
!     ! p = Parser(file)
!     ! p.run(finalTime = 18e-6)
!     block
!         integer :: i
!         open(unit = 1, file =  './testData/outputs/spice/RC_termination_gauss.txt')
!         do i = 1, size(solver%bundles(1)%probes(1)%t)
!             write(1,*) solver%bundles(1)%probes(1)%t(i)," ", &
!                     solver%bundles(1)%probes(1)%val(i,1) ," ", &
!                     solver%bundles(1)%probes(2)%val(i,1) ," ", &
!                     solver%bundles(1)%probes(3)%val(i,1)
!         end do

!     end block

! end function

integer function test_coaxial_line_paul_8_6_square() bind(C) result(error_cnt)    
    use mtln_solver_mod
    use mtln_testingTools_mod
    use preprocess_mod
    implicit none

    ! character(len=*), parameter :: filename = PATH_TO_TEST_DATA//'excitations/coaxial_line_paul_8_6_0.5_square.smb.json'
    character(len=256), parameter :: square_excitation = trim(PATH_TO_TEST_DATA//'excitations/coaxial_line_paul_8_6_0.25_square.exc')
    
    type(cable_t), target :: cable
    type(terminal_node_t) :: node_left, node_right
    type(terminal_network_t) :: nw_left, nw_right
    type(parsed_mtln_t) :: parsed
    type(terminal_connection_t) :: connection_left, connection_right
    type(terminal_network_t) :: network_left, network_right

    type(parsed_probe_t) :: probe_v, probe_i

    type(preprocess_t) :: pre
    type(mtln_t) :: solver
    real,dimension(1,1) :: lpul = 0.25e-6
    real,dimension(1,1) :: gpul = 0.0
    real,dimension(1,1) :: rpul = 0.0
    real,dimension(1,1) :: cpul = 100.0e-12

    character(20) :: charR, charL, charC, lineC
    integer :: i
    type(external_field_segment_t), dimension(100) :: external_field_segments
    class(termination_t), allocatable :: d
   

    error_cnt = 0

    cable%name = "wire0"
    cable%inductance_per_meter = lpul
    cable%conductance_per_meter = gpul
    cable%resistance_per_meter = rpul
    cable%capacitance_per_meter = cpul
    do i = 1, 100
        external_field_segments(i)%position = (/i, 1, 1/)
        external_field_segments(i)%direction = 1
        external_field_segments(i)%Efield_main2wire => null()
        external_field_segments(i)%Efield_wire2main => null()
    end do
    cable%external_field_segments = external_field_segments
    cable%step_size = [(4.0, i = 1, 100)]

    parsed%time_step = 2e-8
    parsed%number_of_steps = 20e-6/parsed%time_step

    node_left%belongs_to_cable => cable
    node_left%conductor_in_cable = 1
    node_left%side = TERMINAL_NODE_SIDE_INI

    node_left%termination%path_to_excitation=square_excitation
    node_left%termination%termination_type = TERMINATION_SERIES
    node_left%termination%resistance = 150
    node_left%termination%inductance = 0.0 
    node_left%termination%capacitance = 1e22

                                             

    write(charR, '(F10.6)') node_left%termination%resistance

    ! write(*,*) node_left%termination%path_to_excitation

    connection_left%nodes = [node_left]

    node_right%belongs_to_cable => cable
    node_right%conductor_in_cable = 1
    node_right%side = TERMINAL_NODE_SIDE_END
    node_right%termination = termination_t(termination_type=TERMINATION_SHORT)

    allocate(connection_right%nodes(1))
    connection_right%nodes = [node_right]
    
    allocate(network_left%connections(1))
    network_left%connections = [connection_left]
    allocate(network_right%connections(1))
    network_right%connections = [connection_right]

    probe_v%attached_to_cable => cable
    probe_v%index = 1 
    probe_v%probe_type = PROBE_TYPE_VOLTAGE

    probe_i%attached_to_cable => cable
    probe_i%index = 100
    probe_i%probe_type = PROBE_TYPE_CURRENT

    allocate(parsed%networks(2))
    parsed%networks = [network_left, network_right]
    allocate(parsed%cables(1))
    parsed%cables = [cable]
    allocate(parsed%probes(2))
    parsed%probes = [probe_v, probe_i]

    solver = mtlnCtor(parsed)
    call solver%runUntil(parsed%time_step * parsed%number_of_steps)

    block
        integer :: i
        real, dimension(:), allocatable :: start_times, end_times, expected_voltages, aux_times
        integer :: j, start, end, idx

        ! open(unit = 1, file =  'testData/outputs/paul/paul_8.6_square.txt')
        ! do i = 1, size(solver%bundles(1)%probes(1)%t)
        !     write(1,*) solver%bundles(1)%probes(1)%t(i)," ", &
        !                solver%bundles(1)%probes(1)%val(i,1) ," ", &
        !                solver%bundles(1)%probes(2)%val(i,1)
        ! end do

        
        start_times = [0.1, 4.1, 6.1, 8.1, 10.1, 12.1, 14.1, 16.1]
        end_times = [3.9, 5.9, 7.9, 9.9, 11.9, 13.9, 15.9, 18.9]
        expected_voltages = [25.0, -12.5, -37.5, -18.75, 18.75, 9.375, -9.375, -4.6875]
        allocate(aux_times(size(solver%bundles(1)%probes(1)%t)), source = 0.0)
        do j = 1, size(start_times) 
            aux_times = 0.5*(start_times(j)*1e-6 + end_times(j)*1e-6)
            idx = minloc(abs(solver%bundles(1)%probes(1)%t - aux_times),1)
            if (.not. (abs((solver%bundles(1)%probes(1)%val(idx,1) - expected_voltages(j))/expected_voltages(j)) <= 5e-2)) then 
                error_cnt =  error_cnt + 1
            end if
        end do

        ! close(unit =1 , status='delete')
    end block



end function

integer function test_coaxial_line_paul_8_6_triangle() bind(C) result(error_cnt)    
    use mtln_solver_mod
    use mtln_testingTools_mod
    use preprocess_mod
    implicit none

    character(len=*), parameter :: square_excitation = PATH_TO_TEST_DATA//'excitations/coaxial_line_paul_8_6_0.05_triangle.exc'
    
    type(cable_t), target :: cable
    type(terminal_node_t) :: node_left, node_right
    type(terminal_network_t) :: nw_left, nw_right
    type(parsed_mtln_t) :: parsed
    type(terminal_connection_t) :: connection_left, connection_right
    type(terminal_network_t) :: network_left, network_right

    type(parsed_probe_t) :: probe_v, probe_i

    type(preprocess_t) :: pre
    type(mtln_t) :: solver
    real,dimension(1,1) :: lpul = 0.25e-6
    real,dimension(1,1) :: gpul = 0.0
    real,dimension(1,1) :: rpul = 0.0
    real,dimension(1,1) :: cpul = 100.0e-12

    character(20) :: charR, charL, charC, lineC
    integer :: i
    type(external_field_segment_t), dimension(100) :: external_field_segments
    error_cnt = 0

    cable%name = "wire0"
    cable%inductance_per_meter = lpul
    cable%conductance_per_meter = gpul
    cable%resistance_per_meter = rpul
    cable%capacitance_per_meter = cpul
    do i = 1, 100
        external_field_segments(i)%position = (/i, 1, 1/)
        external_field_segments(i)%direction = 1
        external_field_segments(i)%Efield_main2wire => null()
        external_field_segments(i)%Efield_wire2main => null()
    end do
    cable%external_field_segments = external_field_segments
    cable%step_size = [(4.0, i = 1, 100)]

    parsed%time_step = 2e-8
    parsed%number_of_steps = 20e-6/parsed%time_step


    node_left%belongs_to_cable => cable
    node_left%conductor_in_cable = 1
    node_left%side = TERMINAL_NODE_SIDE_INI

    node_left%termination%path_to_excitation=square_excitation
    node_left%termination%termination_type = TERMINATION_SERIES
    node_left%termination%resistance = 150
    node_left%termination%inductance = 0.0 
    node_left%termination%capacitance = 1e22


    connection_left%nodes = [node_left]

    node_right%belongs_to_cable => cable
    node_right%conductor_in_cable = 1
    node_right%side = TERMINAL_NODE_SIDE_END
    node_right%termination = termination_t(termination_type=TERMINATION_SHORT)

    connection_right%nodes = [node_right]
    
    network_left%connections = [connection_left]
    network_right%connections = [connection_right]

    probe_v%attached_to_cable => cable
    probe_v%index = 1 
    probe_v%probe_type = PROBE_TYPE_VOLTAGE

    probe_i%attached_to_cable => cable
    probe_i%index = 100
    probe_i%probe_type = PROBE_TYPE_CURRENT

    allocate(parsed%networks(2))
    parsed%networks = [network_left, network_right]
    allocate(parsed%cables(1))
    parsed%cables = [cable]
    allocate(parsed%probes(2))
    parsed%probes = [probe_v, probe_i]

    solver = mtlnCtor(parsed)
    call solver%runUntil(parsed%time_step * parsed%number_of_steps)
    block
        integer :: i
        real, dimension(:), allocatable :: times, expected_voltages, aux_times
        integer :: j, start, end, idx

        ! open(unit = 1, file =  'testData/outputs/paul/paul_8.6_triangle.txt')
        ! do i = 1, size(solver%bundles(1)%probes(1)%t)
        !     write(1,*) solver%bundles(1)%probes(1)%t(i)," ", &
        !                solver%bundles(1)%probes(1)%val(i,1) ," ", &
        !                solver%bundles(1)%probes(2)%val(i,1)
        ! end do
        
        times = [4.0, 5.9, 6.1, 8.0, 10.1, 12.0]
        expected_voltages = [16.67, 12.5, -12.5, -25.0, 6.25, 12.5]
        allocate(aux_times(size(solver%bundles(1)%probes(1)%t)), source = 0.0)
        do j = 1, size(times) 
            aux_times = times(j)*1e-6
            idx = minloc(abs(solver%bundles(1)%probes(1)%t - aux_times),1)
            ! write(*,*) 'expected: ', expected_voltages(j), ' simulated: ', solver%bundles(1)%probes(1)%val(idx,1)
            if (.not. (abs((solver%bundles(1)%probes(1)%val(idx,1) - expected_voltages(j))/expected_voltages(j)) <= 5e-2)) then 
                error_cnt =  error_cnt + 1
            end if
        end do

        ! close(unit =1 , status='delete')

    end block



end function

! integer function test_2_conductor_line_paul_9_6_1c() bind(C) result(error_cnt)    
!     use mtln_solver_mod
!     use mtln_testingTools_mod
!     use preprocess_mod
!     implicit none

!     ! character(len=*), parameter :: filename = PATH_TO_TEST_DATA//'excitations/coaxial_line_paul_8_6_square.smb.json'
!     ! character(len=*), parameter :: pulse_excitation = PATH_TO_TEST_DATA//'excitations/2_conductor_line_paul_9_6_pulse.exc'
!     character(len=*), parameter :: pulse_excitation = PATH_TO_TEST_DATA//'excitations/2_conductor_line_paul_9_6_gauss.exc'
    
!     type(cable_t), target :: cable
!     type(terminal_node_t) :: node_left_1, node_right_1
!     type(terminal_network_t) :: nw_left, nw_right
!     type(parsed_mtln_t) :: parsed
!     type(terminal_connection_t) :: connection_left_1, connection_right_1
!     type(terminal_network_t) :: network_left, network_right

!     type(parsed_probe_t) :: probe_v_left, probe_v_right, probe_v_mid, probe_i_left, probe_i_right, probe_i_mid
!     type(preprocess_t) :: pre
!     type(mtln_t) :: solver
!     real,dimension(1,1) :: lpul = 0.805756e-6
!     real,dimension(1,1) :: gpul = 0.0
!     real,dimension(1,1) :: rpul = 0.0
!     real,dimension(1,1) :: cpul = 117.791e-12


!     character(20) :: charR, charL, charC, lineC
!     real :: final_time
!     error_cnt = 0

!     cable%name = "wire0"
!     cable%inductance_per_meter = lpul
!     cable%conductance_per_meter = gpul
!     cable%resistance_per_meter = rpul
!     cable%capacitance_per_meter = cpul
!     block
!         integer :: i
!         real :: length = 0.5
!         integer :: ndz = 795
!         cable%step_size = [(length/ndz, i = 1, ndz)]
!     end block

!     final_time = 40e-9
!     parsed%number_of_steps = 9410
!     parsed%time_step = final_time/parsed%number_of_steps


                                                   

!     node_left_1%belongs_to_cable => cable
!     node_left_1%conductor_in_cable = 1
!     node_left_1%side = TERMINAL_NODE_SIDE_INI
!     node_left_1%termination = termination_with_source_t(path_to_excitation=pulse_excitation, & 
!                                                  termination_type = TERMINATION_SERIES, &
!                                                  resistance = 50, &
!                                                  inductance = 0.0, &
!                                                  capacitance = 1e22)


!     node_right_1%belongs_to_cable => cable
!     node_right_1%conductor_in_cable = 1
!     node_right_1%side = TERMINAL_NODE_SIDE_END
!     node_right_1%termination = termination_t(termination_type = TERMINATION_RLsCp, &
!                                              resistance = 50, &
!                                              inductance = 0.0, &
!                                              capacitance = 50e-12)
!     ! node_right_1%termination = termination_t(termination_type = TERMINATION_RLsCP, &
!     !                                          resistance = 100, &
!     !                                          inductance = 100e-6, &
!     !                                          capacitance = 10e-12)


    

!     connection_left_1%nodes = [node_left_1]
!     connection_right_1%nodes = [node_right_1]
    
!     network_left%connections = [connection_left_1]
!     network_right%connections = [connection_right_1]

!     probe_v_left%attached_to_cable => cable
!     probe_v_left%index = 1 
!     probe_v_left%probe_type = PROBE_TYPE_VOLTAGE

!     probe_v_right%attached_to_cable => cable
!     probe_v_right%index = 796
!     probe_v_right%probe_type = PROBE_TYPE_VOLTAGE

!     probe_i_left%attached_to_cable => cable
!     probe_i_left%index = 1 
!     probe_i_left%probe_type = PROBE_TYPE_CURRENT

!     probe_i_right%attached_to_cable => cable
!     probe_i_right%index = 795
!     probe_i_right%probe_type = PROBE_TYPE_CURRENT


!     parsed%networks = [network_left, network_right]
!     parsed%cables = [cable]
!     parsed%probes = [probe_v_left, probe_v_right, & 
!                      probe_i_left, probe_i_right]

!     solver = mtlnCtor(parsed)
!     call solver%runUntil(parsed%time_step * parsed%number_of_steps)
!     write(*,*) error_cnt

!     block
!         integer :: i
!         open(unit = 1, file =  'testData/outputs/paul/paul_9.6_gauss_1c_R_50_L_0_C_50p_par.txt')
!         ! open(unit = 1, file =  './probes/probes_9.6_gauss_1c_RLsCp_R_100_L_100u_C_10p.txt')
!         do i = 1, size(solver%bundles(1)%probes(1)%t)
!             write(1,*) solver%bundles(1)%probes(1)%t(i)," ", &
!                        solver%bundles(1)%probes(1)%val(i,1) ," ", &
!                        solver%bundles(1)%probes(2)%val(i,1) ," ", &
!                        solver%bundles(1)%probes(3)%val(i,1) ," ", &
!                        solver%bundles(1)%probes(4)%val(i,1)
!         end do

!     end block


! end function 

! integer function test_2_conductor_line_paul_9_6() bind(C) result(error_cnt)    
!     use mtln_solver_mod
!     use mtln_testingTools_mod
!     use preprocess_mod
!     implicit none

!     character(len=*), parameter :: pulse_excitation = PATH_TO_TEST_DATA//'excitations/2_conductor_line_paul_9_6_pulse.exc'
    
!     type(cable_t), target :: cable
!     type(terminal_node_t) :: node_left_1, node_right_1, node_left_2, node_right_2
!     type(terminal_network_t) :: nw_left, nw_right
!     type(parsed_mtln_t) :: parsed
!     type(terminal_connection_t) :: connection_left_1, connection_right_1, connection_left_2, connection_right_2
!     type(terminal_network_t) :: network_left, network_right

!     type(parsed_probe_t) :: probe_v_left, probe_v_right, probe_v_mid, probe_i_left, probe_i_right, probe_i_mid
!     type(preprocess_t) :: pre
!     type(mtln_t) :: solver
!     real,dimension(2,2) :: lpul = reshape( source = [0.805756e-6, 0.538771e-6, 0.538771e-6, 1.07754e-6], shape = [2,2])
!     real,dimension(2,2) :: gpul = reshape( source = [0.0, 0.0,0.0,0.0], shape = [2,2])
!     real,dimension(2,2) :: rpul = reshape( source = [0.0, 0.0,0.0,0.0], shape = [2,2])
!     real,dimension(2,2) :: cpul = reshape( source = [117.791e-12, -58.8956e-12, -58.8956e-12, 71.8544e-12], shape = [2,2])


!     character(20) :: charR, charL, charC, lineC
!     real :: final_time
!     error_cnt = 0

!     cable%name = "wire0"
!     cable%inductance_per_meter = lpul
!     cable%conductance_per_meter = gpul
!     cable%resistance_per_meter = rpul
!     cable%capacitance_per_meter = cpul
!     block
!         integer :: i
!         real :: length = 0.5
!         integer :: ndz = 795
!         cable%step_size = [(length/ndz, i = 1, ndz)]
!     end block

!     final_time = 40e-9
!     parsed%number_of_steps = 9410
!     parsed%time_step = final_time/parsed%number_of_steps


!     node_left_1%belongs_to_cable => cable
!     node_left_1%conductor_in_cable = 1
!     node_left_1%side = TERMINAL_NODE_SIDE_INI
!     node_left_1%termination = termination_t(termination_type = TERMINATION_SERIES, &
!                                             resistance = 5, &
!                                             inductance = 0.0, &
!                                             capacitance = 1e22)
                                                   

!     node_left_2%belongs_to_cable => cable
!     node_left_2%conductor_in_cable = 2
!     node_left_2%side = TERMINAL_NODE_SIDE_INI
!     node_left_2%termination = termination_with_source_t(path_to_excitation=pulse_excitation, & 
!                                                  termination_type = TERMINATION_SERIES, &
!                                                  resistance = 50, &
!                                                  inductance = 0.0, &
!                                                  capacitance = 1e22)


!     node_right_1%belongs_to_cable => cable
!     node_right_1%conductor_in_cable = 1
!     node_right_1%side = TERMINAL_NODE_SIDE_END
!     ! node_right_1%termination = termination_t(termination_type = TERMINATION_SERIES, &
!     !                                          resistance = 50, &
!     !                                          inductance = 0, &
!     !                                          capacitance = 1e22)
!     node_right_1%termination = termination_t(termination_type = TERMINATION_RLsCP, &
!                                              resistance = 10, &
!                                              inductance = 1e-6, &
!                                              capacitance = 100e-12)

!     node_right_2%belongs_to_cable => cable
!     node_right_2%conductor_in_cable = 2
!     node_right_2%side = TERMINAL_NODE_SIDE_END
!     node_right_2%termination = termination_t(termination_type = TERMINATION_SERIES, &
!                                              resistance = 50, &
!                                              inductance = 0.0, &
!                                              capacitance = 1e22)

    

!     connection_left_1%nodes = [node_left_1]
!     connection_left_2%nodes = [node_left_2]
!     connection_right_1%nodes = [node_right_1]
!     connection_right_2%nodes = [node_right_2]
    
!     network_left%connections = [connection_left_1,connection_left_2]
!     network_right%connections = [connection_right_1, connection_right_2]

!     probe_v_left%attached_to_cable => cable
!     probe_v_left%index = 1 
!     probe_v_left%probe_type = PROBE_TYPE_VOLTAGE

!     probe_v_mid%attached_to_cable => cable
!     probe_v_mid%index = 398
!     probe_v_mid%probe_type = PROBE_TYPE_VOLTAGE

!     probe_v_right%attached_to_cable => cable
!     probe_v_right%index = 796
!     probe_v_right%probe_type = PROBE_TYPE_VOLTAGE

!     probe_i_left%attached_to_cable => cable
!     probe_i_left%index = 1 
!     probe_i_left%probe_type = PROBE_TYPE_CURRENT

!     probe_i_mid%attached_to_cable => cable
!     probe_i_mid%index = 398
!     probe_i_mid%probe_type = PROBE_TYPE_CURRENT

!     probe_i_right%attached_to_cable => cable
!     probe_i_right%index = 795
!     probe_i_right%probe_type = PROBE_TYPE_CURRENT


!     ! parsed%networks = [network_left]
!     parsed%networks = [network_left, network_right]
!     parsed%cables = [cable]
!     parsed%probes = [probe_v_left, probe_v_mid, probe_v_right, & 
!                      probe_i_left, probe_i_mid, probe_i_right]    ! pre = preprocess(parsed)

!     solver = mtlnCtor(parsed)
!     call solver%runUntil(parsed%time_step * parsed%number_of_steps)
!     write(*,*) error_cnt
!     ! p = Parser(file)
!     ! p.run(finalTime = 18e-6)
!     block
!         integer :: i
!         open(unit = 1, file =  'testData/outputs/paul/paul_9.6_pulse.txt')
!         do i = 1, size(solver%bundles(1)%probes(1)%t)
!             write(1,*) solver%bundles(1)%probes(1)%t(i)," ", &
!                        solver%bundles(1)%probes(1)%val(i,1) ," ", &
!                        solver%bundles(1)%probes(1)%val(i,2) ," ", &
!                        solver%bundles(1)%probes(2)%val(i,1) ," ", &
!                        solver%bundles(1)%probes(2)%val(i,2) ," ", &
!                        solver%bundles(1)%probes(3)%val(i,1) ," ", &
!                        solver%bundles(1)%probes(3)%val(i,2) ," ", &
!                        solver%bundles(1)%probes(4)%val(i,1) ," ", &
!                        solver%bundles(1)%probes(4)%val(i,2) ," ", &
!                        solver%bundles(1)%probes(5)%val(i,1) ," ", &
!                        solver%bundles(1)%probes(5)%val(i,2) ," ", &
!                        solver%bundles(1)%probes(6)%val(i,1) ," ", &
!                        solver%bundles(1)%probes(6)%val(i,2)
!         end do

!     end block


! end function 

! integer function test_2_conductor_line_paul_9_11_20ns() bind(C) result(error_cnt)    
!     use mtln_solver_mod
!     use mtln_testingTools_mod
!     use preprocess_mod
!     implicit none

!     character(len=*), parameter :: pulse_excitation = PATH_TO_TEST_DATA//'excitations/2_conductor_line_paul_9_11_20ns.exc'
    
!     type(cable_t), target :: cable
!     type(terminal_node_t) :: node_left_1, node_right_1, node_left_2, node_right_2
!     type(terminal_network_t) :: nw_left, nw_right
!     type(parsed_mtln_t) :: parsed
!     type(terminal_connection_t) :: connection_left_1, connection_right_1, connection_left_2, connection_right_2
!     type(terminal_network_t) :: network_left, network_right

!     type(parsed_probe_t) :: probe_v_101, probe_v_201
!     type(preprocess_t) :: pre
!     type(mtln_t) :: solver
!     real,dimension(2,2) :: lpul = reshape( source = [0.7485e-6, 0.5077e-6, 0.5077e-6, 1.0154e-6], shape = [2,2])
!     real,dimension(2,2) :: gpul = reshape( source = [0.0, 0.0,0.0,0.0], shape = [2,2])
!     real,dimension(2,2) :: rpul = reshape( source = [0.0, 0.0,0.0,0.0], shape = [2,2])
!     real,dimension(2,2) :: cpul = reshape( source = [37.432e-12, -18.716e-12, -18.716e-12, 24.982e-12], shape = [2,2])

!     real :: final_time

!     character(20) :: charR, charL, charC, lineC

!     error_cnt = 0

!     cable%name = "wire0"
!     cable%inductance_per_meter = lpul
!     cable%conductance_per_meter = gpul
!     cable%resistance_per_meter = rpul
!     cable%capacitance_per_meter = cpul
!     block
!         integer :: i
!         cable%step_size = [(1.0, i = 1, 2)]
!     end block

!     final_time = 200e-9
!     parsed%number_of_steps = 54
!     parsed%time_step = final_time/parsed%number_of_steps


!     node_left_1%belongs_to_cable => cable
!     node_left_1%conductor_in_cable = 1
!     node_left_1%side = TERMINAL_NODE_SIDE_INI
!     node_left_1%termination = termination_t(termination_type = TERMINATION_SERIES, &
!                                             resistance = 50, &
!                                             inductance = 0.0, &
!                                             capacitance = 1e22)
                                                   

!     node_left_2%belongs_to_cable => cable
!     node_left_2%conductor_in_cable = 2
!     node_left_2%side = TERMINAL_NODE_SIDE_INI
!     node_left_2%termination = termination_with_source_t(path_to_excitation=pulse_excitation, & 
!                                                  termination_type = TERMINATION_SERIES, &
!                                                  resistance = 50, &
!                                                  inductance = 0.0, &
!                                                  capacitance = 1e22)


!     node_right_1%belongs_to_cable => cable
!     node_right_1%conductor_in_cable = 1
!     node_right_1%side = TERMINAL_NODE_SIDE_END
!     node_right_1%termination = termination_t(termination_type = TERMINATION_SERIES, &
!                                              resistance = 50, &
!                                              inductance = 0.0, &
!                                              capacitance = 1e22)

!     node_right_2%belongs_to_cable => cable
!     node_right_2%conductor_in_cable = 2
!     node_right_2%side = TERMINAL_NODE_SIDE_END
!     node_right_2%termination = termination_t(termination_type = TERMINATION_SERIES, &
!                                              resistance = 50, &
!                                              inductance = 0.0, &
!                                              capacitance = 1e22)

    

!     connection_left_1%nodes = [node_left_1]
!     connection_left_2%nodes = [node_left_2]
!     connection_right_1%nodes = [node_right_1]
!     connection_right_2%nodes = [node_right_2]
    
!     network_left%connections = [connection_left_1,connection_left_2]
!     network_right%connections = [connection_right_1, connection_right_2]

!     probe_v_101%attached_to_cable => cable
!     probe_v_101%index = 1 
!     probe_v_101%probe_type = PROBE_TYPE_VOLTAGE

!     parsed%networks = [network_left, network_right]
!     parsed%cables = [cable]
!     parsed%probes = [probe_v_101]

!        solver = mtlnCtor(parsed)
!     call solver%runUntil(parsed%time_step * parsed%number_of_steps)
!     write(*,*) error_cnt

!     block
!         integer :: i
!         open(unit = 1, file =  'testData/outputs/paul/paul_9.11_20ns.txt')
!         do i = 1, size(solver%bundles(1)%probes(1)%t)
!             write(1,*) solver%bundles(1)%probes(1)%t(i)," ", &
!                        solver%bundles(1)%probes(1)%val(i,1)," ", &
!                        solver%bundles(1)%probes(1)%val(i,2)
!         end do

!     end block


! end function 

! integer function test_2_conductor_line_paul_9_11_1ns() bind(C) result(error_cnt)    
!     use mtln_solver_mod
!     use mtln_testingTools_mod
!     use preprocess_mod
!     implicit none

!     character(len=*), parameter :: pulse_excitation = PATH_TO_TEST_DATA//'excitations/2_conductor_line_paul_9_11_1ns.exc'
    
!     type(cable_t), target :: cable
!     type(terminal_node_t) :: node_left_1, node_right_1, node_left_2, node_right_2
!     type(terminal_network_t) :: nw_left, nw_right
!     type(parsed_mtln_t) :: parsed
!     type(terminal_connection_t) :: connection_left_1, connection_right_1, connection_left_2, connection_right_2
!     type(terminal_network_t) :: network_left, network_right

!     type(parsed_probe_t) :: probe_v_101, probe_v_201
!     type(preprocess_t) :: pre
!     type(mtln_t) :: solver
!     real,dimension(2,2) :: lpul = reshape( source = [0.7485e-6, 0.5077e-6, 0.5077e-6, 1.0154e-6], shape = [2,2])
!     real,dimension(2,2) :: gpul = reshape( source = [0.0, 0.0,0.0,0.0], shape = [2,2])
!     real,dimension(2,2) :: rpul = reshape( source = [0.0, 0.0,0.0,0.0], shape = [2,2])
!     real,dimension(2,2) :: cpul = reshape( source = [37.432e-12, -18.716e-12, -18.716e-12, 24.982e-12], shape = [2,2])

!     real :: final_time

!     character(20) :: charR, charL, charC, lineC

!     error_cnt = 0

!     cable%name = "wire0"
!     cable%inductance_per_meter = lpul
!     cable%conductance_per_meter = gpul
!     cable%resistance_per_meter = rpul
!     cable%capacitance_per_meter = cpul
!     block
!         integer :: i
!         cable%step_size = [(0.02, i = 1, 100)]
!     end block

!     final_time = 200e-9
!     parsed%number_of_steps = 2600
!     parsed%time_step = final_time/parsed%number_of_steps


!     node_left_1%belongs_to_cable => cable
!     node_left_1%conductor_in_cable = 1
!     node_left_1%side = TERMINAL_NODE_SIDE_INI
!     node_left_1%termination = termination_t(termination_type = TERMINATION_SERIES, &
!                                             resistance = 50, &
!                                             inductance = 0.0, &
!                                             capacitance = 1e22)
                                                   

!     node_left_2%belongs_to_cable => cable
!     node_left_2%conductor_in_cable = 2
!     node_left_2%side = TERMINAL_NODE_SIDE_INI
!     node_left_2%termination = termination_with_source_t(path_to_excitation=pulse_excitation, & 
!                                                  termination_type = TERMINATION_SERIES, &
!                                                  resistance = 50, &
!                                                  inductance = 0.0, &
!                                                  capacitance = 1e22)


!     node_right_1%belongs_to_cable => cable
!     node_right_1%conductor_in_cable = 1
!     node_right_1%side = TERMINAL_NODE_SIDE_END
!     node_right_1%termination = termination_t(termination_type = TERMINATION_SERIES, &
!                                              resistance = 50, &
!                                              inductance = 0.0, &
!                                              capacitance = 1e22)

!     node_right_2%belongs_to_cable => cable
!     node_right_2%conductor_in_cable = 2
!     node_right_2%side = TERMINAL_NODE_SIDE_END
!     node_right_2%termination = termination_t(termination_type = TERMINATION_SERIES, &
!                                              resistance = 50, &
!                                              inductance = 0.0, &
!                                              capacitance = 1e22)

    

!     connection_left_1%nodes = [node_left_1]
!     connection_left_2%nodes = [node_left_2]
!     connection_right_1%nodes = [node_right_1]
!     connection_right_2%nodes = [node_right_2]
    
!     network_left%connections = [connection_left_1,connection_left_2]
!     network_right%connections = [connection_right_1, connection_right_2]

!     probe_v_101%attached_to_cable => cable
!     probe_v_101%index = 1 
!     probe_v_101%probe_type = PROBE_TYPE_VOLTAGE

!     parsed%networks = [network_left, network_right]
!     parsed%cables = [cable]
!     parsed%probes = [probe_v_101]

!     solver = mtlnCtor(parsed)
!     call solver%runUntil(parsed%time_step * parsed%number_of_steps)
!     write(*,*) error_cnt

!     block
!         integer :: i
!         open(unit = 1, file =  'testData/outputs/paul/paul_9.11_1ns.txt')
!         do i = 1, size(solver%bundles(1)%probes(1)%t)
!             write(1,*) solver%bundles(1)%probes(1)%t(i)," ", &
!                        solver%bundles(1)%probes(1)%val(i,1)," ", &
!                        solver%bundles(1)%probes(1)%val(i,2)
!         end do

!     end block


! end function 
