integer function test_networks_pointer_copy() bind(C) result(error_cnt)

    type node_t 
        real, pointer :: v
        real, pointer :: i
    end type

    type(node_t), dimension(2) :: nodes
    type(node_t), dimension(2) :: nodes_copy
    real, dimension(3), target :: v1, v2
    real, dimension(2), target :: i1, i2

    error_cnt = 0 

    nodes(1)%v => v1(3)
    nodes(1)%i => i1(2)
    
    nodes(2)%v => v2(1)
    nodes(2)%i => i2(1)

    v1 = [1,2,3]
    v2 = [4,5,6]

    i1 = [-1,-2]
    i2 = [-3,-4]

    ! write(*,*) 'v1: ', v1, ' i1: ', i1
    ! write(*,*) 'v2: ', v2, ' i2: ', i2

    nodes(1)%v = 2*nodes(1)%v
    nodes(2)%v = 2*nodes(2)%v
    i1(2) = 0.5*i1(2)
    i2(1) = 0.5*i2(1)

    ! write(*,*) 'nodes: '
    ! write(*,*) 'v: ', nodes(1)%v, ' ', nodes(2)%v
    ! write(*,*) 'i: ', nodes(1)%i, ' ', nodes(2)%i
    ! write(*,*) 'v1: ', v1, ' i1: ', i1
    ! write(*,*) 'v2: ', v2, ' i2: ', i2

    if ((v1(3) /= 6.0) .or. (v2(1) /= 8.0)) then 
        error_cnt = error_cnt + 1
    end if
    if ((nodes(1)%i /= -1.0) .or. (nodes(2)%i /=-1.5)) then 
        error_cnt = error_cnt + 1
    end if

    nodes_copy = nodes
    nodes_copy(1)%v = 2*nodes_copy(1)%v
    nodes_copy(2)%v = 2*nodes_copy(2)%v
    i1(2) = 0.5*i1(2)
    i2(1) = 0.5*i2(1)

    if ((v1(3) /= 12.0) .or. (v2(1) /= 16.0)) then 
        error_cnt = error_cnt + 1
    end if
    if ((nodes_copy(1)%i /= -0.5) .or. (nodes_copy(2)%i /=-0.75)) then 
        error_cnt = error_cnt + 1
    end if

end function

integer function test_networks_simple_manager() bind(C) result(error_cnt)

    use network_manager_mod
    use mtln_testingTools_mod, only: buildNetwork, checkNear
    implicit none

    type(network_manager_t) :: manager
    real :: final_time, dt
    character(100), dimension(:), allocatable :: nw_mng_desc, input
    character(20) :: sTime, sdt
    type(network_t) :: nw1, nw2
    type(network_t), dimension(:), allocatable :: nws(:)
    real, dimension(3), target :: tgt_v_1, tgt_v_2, tgt_i_1, tgt_i_2
    integer :: i
    final_time = 200e-6
    dt = 50e-6
    error_cnt = 0

    write(sTime, '(E10.2)') final_time
    write(sdt, '(E10.2)') dt
    
    nw1 = buildNetwork("n1","5k","2.5k","50u","200n", tgt_v_1, tgt_i_1)
    nw2 = buildNetwork("n2","10k","1k","1u","100n", tgt_v_2, tgt_i_2)


    allocate(nw_mng_desc(0))
    nw_mng_desc = [nw_mng_desc, trim("* manager with two networks")]
    nw_mng_desc = [nw_mng_desc, nw1%description]
    nw_mng_desc = [nw_mng_desc, nw2%description]
    nw_mng_desc = [nw_mng_desc, trim(".tran "//trim(sdt)//" "//sTime)]
    nw_mng_desc = [nw_mng_desc, ".save v(n1_int) v(n1_in) v(n1_out) v(n2_int) v(n2_in) v(n2_out) "]
    nw_mng_desc = [nw_mng_desc, ".end"]
    nw_mng_desc = [nw_mng_desc, "NULL"]
    
    allocate(nws(0))
    nws = [nw1, nw2]
    manager = network_managerCtor(nws, nw_mng_desc, final_time, dt)

    do while (manager%time < final_time)
        call manager%advanceVoltage()
        manager%time = manager%time + manager%dt
        if (checkNear(manager%circuit%getTime(), manager%circuit%time, 0.01) .eqv. .false. ) then 
            error_cnt = error_cnt + 1
        end if
    end do

    if ((checkNear(manager%networks(1)%nodes(1)%v, manager%circuit%nodes%values(1)%voltage, 0.01) .eqv. .false.) .or. &
        (checkNear(manager%networks(1)%nodes(2)%v, manager%circuit%nodes%values(2)%voltage, 0.01) .eqv. .false.) .or. &
        (checkNear(manager%networks(1)%nodes(3)%v, manager%circuit%nodes%values(3)%voltage, 0.01) .eqv. .false.) .or. &
        (checkNear(manager%networks(2)%nodes(1)%v, manager%circuit%nodes%values(4)%voltage, 0.01) .eqv. .false.) .or. &
        (checkNear(manager%networks(2)%nodes(2)%v, manager%circuit%nodes%values(5)%voltage, 0.01) .eqv. .false.) .or. &
        (checkNear(manager%networks(2)%nodes(3)%v, manager%circuit%nodes%values(6)%voltage, 0.01) .eqv. .false.)) then 
        error_cnt = error_cnt + 1
    end if

end function