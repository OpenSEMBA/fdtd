integer function test_spice_read_message() bind(C) result(error_cnt)    

    use circuit_m
    use mtln_testingTools_m
    implicit none

    type(circuit_t) :: circuit
    real(kind=rkind) :: resultValue(4)
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
    call circuit%updateNodes()

    resultValue = [24.000000000000000, 9.7469741675197206, 15.000000000000000, 24.000000000000000]
    if (size(circuit%nodes%values) /= 4) then 
        error_cnt = error_cnt + 1
    end if

    do i = 1, 4                      
        if (checkNear(circuit%nodes%values(i)%voltage, resultValue(i), 0.01_rkind) .eqv. .false.) then 
            error_cnt = error_cnt + 1
        end if
    end do

end function

integer function test_spice_dc() bind(C) result(error_cnt)    

    use circuit_m
    use mtln_testingTools_m
    implicit none

    type(circuit_t) :: circuit
    character(len=50) :: netlist
    real(kind=rkind) :: resultValue(4)
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
    call circuit%updateNodes()

    resultValue = [24.000000000000000, 9.7469741675197206, 15.000000000000000, 24.000000000000000]
    if (size(circuit%nodes%values) /= 4) then 
        error_cnt = error_cnt + 1
    end if

    do i = 1, 4                      
        if (checkNear(circuit%nodes%values(i)%voltage, resultValue(i), 0.01_rkind) .eqv. .false.) then 
            error_cnt = error_cnt + 1
        end if
    end do

end function

integer function test_spice_tran() bind(C) result(error_cnt)    

    use circuit_m
    use mtln_testingTools_m
    implicit none

    type(circuit_t) :: circuit
    character(len=*, kind=c_char), parameter :: NETLIST= PATH_TO_TEST_DATA//c_char_'netlists/netlist_tran.cir'
    real(kind=RKIND_TIME) :: finalTime
    real(kind=rkind) :: resultValue(3)
    integer :: i
    type(string_t), dimension(4) :: names
    names(1) = string_t("in", 2)
    names(2) = string_t("int", 3)
    names(3) = string_t("out", 3)
    names(4) = string_t("time", 4)

    resultValue = [5.0,0.092995181699999999,0.053166680000000001]

    circuit%time = 0.0
    circuit%dt = 50e-6
    finalTime = 200e-6

    error_cnt = 0
    call circuit%init(names=names,netlist=netlist)
    call circuit%setStopTimes(finalTime, circuit%dt)
    do while (circuit%time < finalTime)
        call circuit%step()
        call circuit%updateNodes()
        if (checkNear_time(circuit%getTime(), circuit%time, 0.01_RKIND_TIME) .eqv. .false.) then 
            error_cnt = error_cnt + 1
        end if
    end do
    if (checkNear(circuit%getNodeVoltage("in"), resultValue(1), 0.01_rkind) .eqv. .false.) then 
        error_cnt = error_cnt + 1
    end if
    if (checkNear(circuit%getNodeVoltage("int"), resultValue(2), 0.01_rkind) .eqv. .false.) then 
        error_cnt = error_cnt + 1
    end if
    if (checkNear(circuit%getNodeVoltage("out"), resultValue(3), 0.01_rkind) .eqv. .false.) then 
        error_cnt = error_cnt + 1
    end if

    ! call circuit%quit()

end function

integer function test_spice_tran_2() bind(C) result(error_cnt)    

    use circuit_m
    use mtln_testingTools_m
    implicit none

    type(circuit_t) :: circuit
    character(len=*, kind=c_char), parameter :: NETLIST= PATH_TO_TEST_DATA//c_char_'netlists/netlist_tran_2.cir'
    real(kind=RKIND_TIME) :: finalTime
    integer :: i
    real(kind=rkind) :: resultValue(3)
    type(string_t), dimension(4) :: names
    names(1) = string_t("in", 2)
    names(2) = string_t("int", 3)
    names(3) = string_t("out", 3)
    names(4) = string_t("time", 4)
    
    resultValue = [5.0_rkind,0.0039656539400000001_rkind,0.00069279532199999997_rkind]
    
    circuit%time = 0.0
    circuit%dt = 50e-6
    finalTime = 200e-6

    error_cnt = 0
    call circuit%init(names=names,netlist=netlist)
    call circuit%setStopTimes(finalTime, circuit%dt)
    do while (circuit%time < finalTime)
        call circuit%step()
        call circuit%updateNodes()
        if (checkNear_time(circuit%getTime(), circuit%time, 0.01_RKIND_TIME) .eqv. .false.) then 
            error_cnt = error_cnt + 1
        end if
    end do
    if (checkNear(circuit%getNodeVoltage("in"), resultValue(1), 0.01_rkind) .eqv. .false.) then 
        error_cnt = error_cnt + 1
    end if
    if (checkNear(circuit%getNodeVoltage("int"), resultValue(2), 0.01_rkind) .eqv. .false.) then 
        error_cnt = error_cnt + 1
    end if
    if (checkNear(circuit%getNodeVoltage("out"), resultValue(3), 0.01_rkind) .eqv. .false.) then 
        error_cnt = error_cnt + 1
    end if

    ! call circuit%quit()
end function

integer function test_spice_current_source() bind(C) result(error_cnt)    

    use circuit_m
    use mtln_testingTools_m
    implicit none

    type(circuit_t) :: circuit
    character(len=50) :: netlist
    real(kind=RKIND_TIME) :: finalTime
    real(kind=rkind) ::resistance
    integer :: i
    real(kind=rkind) :: current, voltage
    character(50) :: sCurrent
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
        write(sCurrent, *) current
        call command("alter @I1_initial[dc]="//trim(sCurrent) // c_null_char)
        call circuit%step()
        call circuit%updateNodes()
        voltage = circuit%getNodeVoltage("1_initial")
        if (checkNear(voltage, current*resistance, 0.01_rkind) .eqv. .false.) then 
            error_cnt = error_cnt + 1
        end if
        current = 2.0*current
    end do

end function

integer function test_spice_multiple() bind(C) result(error_cnt)

    use circuit_m
    use mtln_testingTools_m
    implicit none
    
    type(circuit_t) :: circuit
    character(len=50) :: netlist
    real(kind=RKIND_TIME) :: dt= 50e-6
    real(kind=RKIND_TIME) :: finalTime = 200e-6
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
        call circuit%updateNodes()
        if (checkNear_time(circuit%getTime(), circuit%time, 0.01_RKIND_TIME) .eqv. .false.) then 
            error_cnt = error_cnt + 1
        end if
    end do
   


end function

integer function test_spice_stop_mod_times() bind(C) result(error_cnt)
    use circuit_m
    use mtln_testingTools_m
    implicit none

    type(circuit_t) :: circuit
    character(len=*, kind=c_char), parameter :: NETLIST= PATH_TO_TEST_DATA//c_char_'netlists/netlist_tran.cir'
    real(kind=RKIND_TIME) :: finalTime
    real(kind=rkind) :: resultValue(3)
    integer :: i
    type(string_t), dimension(4) :: names
    names(1) = string_t("in", 2)
    names(2) = string_t("int", 3)
    names(3) = string_t("out", 3)
    names(4) = string_t("time", 4)

    resultValue = [5.0,0.092995181699999999,0.053166680000000001]

    circuit%time = 0.0
    circuit%dt = 50e-6
    finalTime = 200e-6

    error_cnt = 0
    call circuit%init(names=names,netlist=netlist)
    call circuit%setModStopTimes(circuit%dt)
    do while (circuit%time < finalTime)
        call circuit%step()
        call circuit%updateNodes()
        if (checkNear_time(circuit%getTime(), circuit%time, 0.01_RKIND_TIME) .eqv. .false.) then 
            error_cnt = error_cnt + 1
        end if
    end do
    if (checkNear(circuit%getNodeVoltage("in"), resultValue(1), 0.01_rkind) .eqv. .false.) then 
        error_cnt = error_cnt + 1
    end if
    if (checkNear(circuit%getNodeVoltage("int"), resultValue(2), 0.01_rkind) .eqv. .false.) then 
        error_cnt = error_cnt + 1
    end if
    if (checkNear(circuit%getNodeVoltage("out"), resultValue(3), 0.01_rkind) .eqv. .false.) then 
        error_cnt = error_cnt + 1
    end if

    ! call circuit%quit()

end function


integer function test_load_codemodels() bind(C) result(error_cnt)
    ! Loads the XSPICE analog code models and solves a circuit with a 'gain'
    ! code model. Fails if the code models cannot be loaded or are not
    ! compatible with the ngspice library (the 'a1' line would not parse).
    use circuit_m
    use mtln_testingTools_m
    use ngspice_interface_m, only: command
    implicit none
    include 'ngspice_codemodels_dir.inc'
    type(circuit_t) :: circuit
    character(50), dimension(:), allocatable :: input
    type(string_t), dimension(2) :: names
    names(1) = string_t("in", 2)
    names(2) = string_t("out", 3)

    allocate(input(0))
    input = ["* XSPICE code model check"]
    input = [input, "vin in 0 dc 2"]
    input = [input, "a1 in out amp"]
    input = [input, ".model amp gain(gain=2.5)"]
    input = [input, "rl out 0 1k"]
    input = [input, ".dc vin 2 2 1"]
    input = [input, ".save v(in) v(out)"]
    input = [input, ".end"]
    input = [input, "NULL"]

    error_cnt = 0
    call circuit%init(names=names)
    ! ngspice passes quotes literally to dlopen, so the path must not contain spaces.
    call command(c_char_'codemodel '//NGSPICE_CODEMODELS_DIR//c_char_'analog.cm'//c_null_char)
    call circuit%readInput(input)
    call circuit%step()
    call circuit%updateNodes()

    if (.not. checkNear(circuit%getNodeVoltage("in"), 2.0_rkind, 0.01_rkind)) then
        error_cnt = error_cnt + 1
    end if
    if (.not. checkNear(circuit%getNodeVoltage("out"), 5.0_rkind, 0.01_rkind)) then
        error_cnt = error_cnt + 1
    end if
end function
