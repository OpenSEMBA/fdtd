integer function test_spice_read_message() bind(C) result(error_cnt)    

    use circuit_mod
    use mtln_testingTools_mod
    implicit none

    type(circuit_t) :: circuit
    real :: result(4)
    character(50), dimension(:), allocatable :: input
    integer :: i
    type(string_t), dimension(4) :: names
    names(1) = string_t("node1", 5)
    names(2) = string_t("node2", 5)
    names(3) = string_t("node3", 5)
    names(4) = string_t("v-sweep", 7)

    allocate(input(0))
    input = ["* Multiple dc sources"]
    input = [input, "vn1 node1 0 dc 24"]
    input = [input, "vn2 node3 0 dc 15"]
    input = [input, "rn1 node1 node2 10k"]
    input = [input, "rn2 node2 node3 8.1k"]
    input = [input, "rn3 node2 0 4.7k"]
    input = [input, ".dc vn1 24 24 1"]
    input = [input, ".save v(node3) v(node2) v(node1)"]
    input = [input, ".end"]
    input = [input, "NULL"]
    ! write(*,*) input
    error_cnt = 0

    call circuit%init(names=names)
    call circuit%readInput(input)
    call circuit%step()

    result = [24.000000000000000, 9.7469741675197206, 15.000000000000000, 24.000000000000000]
    if (size(circuit%nodes%values) /= 4) then 
        error_cnt = error_cnt + 1
    end if

    do i = 1, 4                      
        if (checkNear(circuit%nodes%values(i)%voltage, result(i), 0.01) .eqv. .false. ) then 
            error_cnt = error_cnt + 1
        end if
    end do

end function

integer function test_spice_dc() bind(C) result(error_cnt)    

    use circuit_mod
    use mtln_testingTools_mod
    implicit none

    type(circuit_t) :: circuit
    character(len=50) :: netlist
    real :: result(4)
    integer :: i
    type(string_t), dimension(4) :: names
    names(1) = string_t("node1", 5)
    names(2) = string_t("node2", 5)
    names(3) = string_t("node3", 5)
    ! names(3) = string_t("vn1#branch", 10)
    names(4) = string_t("v-sweep", 7)

    error_cnt = 0
    netlist = PATH_TO_TEST_DATA//'netlists/netlist_dc.cir'
    call circuit%init(names = names, netlist = netlist)
    call circuit%step()

    result = [24.000000000000000, 9.7469741675197206, 15.000000000000000, 24.000000000000000]
    if (size(circuit%nodes%values) /= 4) then 
        error_cnt = error_cnt + 1
    end if

    do i = 1, 4                      
        if (checkNear(circuit%nodes%values(i)%voltage, result(i), 0.01) .eqv. .false. ) then 
            error_cnt = error_cnt + 1
        end if
    end do

end function

integer function test_spice_tran() bind(C) result(error_cnt)    

    use circuit_mod
    use mtln_testingTools_mod
    implicit none

    type(circuit_t) :: circuit
    character(len=*, kind=c_char), parameter :: netlist= PATH_TO_TEST_DATA//c_char_'netlists/netlist_tran.cir'
    real :: finalTime
    real :: result(3)
    integer :: i
    type(string_t), dimension(4) :: names
    names(1) = string_t("in", 2)
    names(2) = string_t("int", 3)
    names(3) = string_t("out", 3)
    names(4) = string_t("time", 4)

    result = [5.0,0.092995181699999999,0.053166680000000001]

    circuit%time = 0.0
    circuit%dt = 50e-6
    finalTime = 200e-6

    error_cnt = 0
    call circuit%init(names=names,netlist=netlist)
    call circuit%setStopTimes(finalTime, circuit%dt)
    do while (circuit%time < finalTime)
        call circuit%step()
        circuit%time = circuit%time + circuit%dt
        if (checkNear(circuit%getTime(), circuit%time, 0.01) .eqv. .false. ) then 
            error_cnt = error_cnt + 1
        end if
    end do
    if (checkNear(circuit%getNodeVoltage("in"), result(1), 0.01) .eqv. .false. ) then 
        error_cnt = error_cnt + 1
    end if
    if (checkNear(circuit%getNodeVoltage("int"), result(2), 0.01) .eqv. .false. ) then 
        error_cnt = error_cnt + 1
    end if
    if (checkNear(circuit%getNodeVoltage("out"), result(3), 0.01) .eqv. .false. ) then 
        error_cnt = error_cnt + 1
    end if

    ! call circuit%quit()

end function

integer function test_spice_tran_2() bind(C) result(error_cnt)    

    use circuit_mod
    use mtln_testingTools_mod
    implicit none

    type(circuit_t) :: circuit
    character(len=*, kind=c_char), parameter :: netlist= PATH_TO_TEST_DATA//c_char_'netlists/netlist_tran_2.cir'
    real :: finalTime
    integer :: i
    real :: result(3)
    type(string_t), dimension(4) :: names
    names(1) = string_t("in", 2)
    names(2) = string_t("int", 3)
    names(3) = string_t("out", 3)
    names(4) = string_t("time", 4)
    
    result = [5.0,0.0039656539400000001,0.00069279532199999997]
    
    circuit%time = 0.0
    circuit%dt = 50e-6
    finalTime = 200e-6

    error_cnt = 0
    call circuit%init(names=names,netlist=netlist)
    call circuit%setStopTimes(finalTime, circuit%dt)
    do while (circuit%time < finalTime)
        call circuit%step()
        circuit%time = circuit%time + circuit%dt
        if (checkNear(circuit%getTime(), circuit%time, 0.01) .eqv. .false. ) then 
            error_cnt = error_cnt + 1
        end if
    end do
    if (checkNear(circuit%getNodeVoltage("in"), result(1), 0.01) .eqv. .false. ) then 
        error_cnt = error_cnt + 1
    end if
    if (checkNear(circuit%getNodeVoltage("int"), result(2), 0.01) .eqv. .false. ) then 
        error_cnt = error_cnt + 1
    end if
    if (checkNear(circuit%getNodeVoltage("out"), result(3), 0.01) .eqv. .false. ) then 
        error_cnt = error_cnt + 1
    end if

    ! call circuit%quit()
end function

integer function test_spice_current_source() bind(C) result(error_cnt)    

    use circuit_mod
    use mtln_testingTools_mod
    implicit none

    type(circuit_t) :: circuit
    character(len=50) :: netlist
    real :: finalTime, resistance
    integer :: i
    real :: current
    type(string_t), dimension(1) :: names
    names(1) = string_t("1_initial", 9)

    circuit%time = 0.0
    circuit%dt = 50e-6
    finalTime = 200e-6

    error_cnt = 0
    netlist = PATH_TO_TEST_DATA//'netlists/netlist_current_source.cir'

    resistance = 10
    current = 0.1
    call circuit%init(names = names, netlist = netlist)
    call circuit%setStopTimes(finalTime, circuit%dt)
    do while (circuit%time < finalTime)
        call circuit%updateNodeCurrent("1_initial", current)
        call circuit%step()
        circuit%time = circuit%time + circuit%dt
        if (checkNear(circuit%getNodeVoltage("1_initial"), current*resistance, 0.01) .eqv. .false. ) then 
            error_cnt = error_cnt + 1
        end if
        current = 2.0*current
    end do

end function

integer function test_spice_multiple() bind(C) result(error_cnt)

    use circuit_mod
    use mtln_testingTools_mod
    implicit none
    
    type(circuit_t) :: circuit
    character(len=50) :: netlist
    real :: dt = 50e-6
    real :: finalTime = 200e-6
    type(string_t), dimension(7) :: names
    names(1) = string_t("n1_in", 4)
    names(2) = string_t("n1_int", 5)
    names(3) = string_t("n1_out", 5)
    names(4) = string_t("time", 4)
    names(5) = string_t("n2_in", 5)
    names(6) = string_t("n2_int", 5)
    names(7) = string_t("n2_out", 5)
    circuit%time = 0.0
    circuit%dt = 50e-6
    finalTime = 200e-6

    error_cnt = 0

    netlist = PATH_TO_TEST_DATA//'netlists/netlist_multiple.cir'
    call circuit%init(names = names, netlist = netlist)
    call circuit%setStopTimes(finalTime, dt)
    do while (circuit%time < finalTime)
        call circuit%step()
        circuit%time = circuit%time + circuit%dt
        if (checkNear(circuit%getTime(), circuit%time, 0.01) .eqv. .false. ) then 
            error_cnt = error_cnt + 1
        end if
    end do
   


end function
