module network_manager_m

    use network_m
    use circuit_m
    use mtln_types_m, only: node_source_t
    use FDETYPES_m, only: RKIND, RKIND_TIEMPO

    implicit none 

    type network_manager_t
        type(network_t), dimension(:), allocatable :: networks
        type(circuit_t) :: circuit
        type(nw_node_t), allocatable :: open_nodes(:)
        real(kind=rkind) :: time, dt
        logical :: has_active_node = .false.

    contains
        procedure :: advanceVoltage => network_advanceVoltage
        procedure :: updateCircuitCurrentsFromNetwork
        procedure :: updateNetworkVoltagesFromCircuit
    end type

    interface network_manager_t
        module procedure network_managerCtor
    end interface


contains

    subroutine appendToString_tArray(arr, str)
        ! This has been implemented because there seems to be a bug in gfortran: 
        ! https://fortran-lang.discourse.group/t/read-data-and-append-it-to-array-best-practice/1915
        ! and arr = [ arr, str ] can't be used.
        type(string_t), allocatable, intent(inout) :: arr(:)
        type(string_t), intent(in) :: str
        type(string_t), allocatable :: old_arr(:)
        
        old_arr = arr
        deallocate(arr)
        allocate(arr(size(old_arr)+1))
        arr(1:size(old_arr)) = old_arr 
        arr(size(old_arr)+1) = str
    end subroutine


    function copy_sources(networks) result(res)
        type(network_t), dimension(:), intent(in) :: networks
        type(node_source_t), dimension(:), allocatable :: res
        integer :: i,j,n
        type(string_t) :: temp
        n = 0
        do i = 1, size(networks)
            do j = 1, size(networks(i)%nodes)
                n = n + 1
            end do
        end do
        allocate(res(n))
        n = 1
        do i = 1, size(networks)
            do j = 1, size(networks(i)%nodes)
                res(n)%path_to_excitation = trim(networks(i)%nodes(j)%source%path_to_excitation)
                res(n)%source_type = networks(i)%nodes(j)%source%source_type
                res(n)%resistance = networks(i)%nodes(j)%source%resistance
                n = n + 1
            end do
        end do
    end function

    function copy_node_names(networks) result(res)
        type(network_t), dimension(:), intent(in) :: networks
        type(string_t), dimension(:), allocatable :: res
        integer :: i,j
        type(string_t) :: temp
        allocate(res(0))
        do i = 1, size(networks)
            do j = 1, size(networks(i)%nodes)
                temp = string_t(trim(networks(i)%nodes(j)%name), len(trim(networks(i)%nodes(j)%name)))
                call appendToString_tArray(res, temp)
            end do
        end do
        call appendToString_tArray(res, string_t("time",4))
    end function


    function network_managerCtor(networks, description, final_time, dt) result(res)
        type(network_t), dimension(:), intent(in) :: networks
        character(*), dimension(:), intent(in) :: description
        real(kind=RKIND_TIEMPO), intent(in) :: final_time, dt
        type(network_manager_t) :: res
        logical :: printInput = .true.
        res%dt = dt
        res%time = 0.0
        res%networks = networks

        res%open_nodes = collectOpenNodes(networks)

        call res%circuit%init(copy_node_names(networks), copy_sources(networks))
        res%circuit%dt = dt
#ifdef CompileWithRelease
        printInput = .false.
#endif        
        call res%circuit%readInput(description, printInput)
        call res%circuit%setModStopTimes(dt)

        contains
        
        function collectOpenNodes(nws) result(res)
            type(network_t), dimension(:), intent(in) :: nws
            integer :: i, j, n
            type(nw_node_t), allocatable :: res(:)
            n = 0
            do i = 1, size(nws)
                do j = 1, nws(i)%number_of_nodes
                    if (nws(i)%nodes(j)%open) n = n + 1
                end do
            end do
            allocate(res(n))
            if (n==0) return
            n = 0
            do i = 1, size(nws)
                do j = 1, nws(i)%number_of_nodes
                    if (nws(i)%nodes(j)%open) then 
                        n = n + 1
                        res(n) = nws(i)%nodes(j)
                    end if
                end do
            end do
            
        end function


    end function

    subroutine updateNetworkVoltages(this)
        class(network_manager_t) :: this
        integer :: i, j
        do i = 1, size(this%networks)
            do j = 1, this%networks(i)%number_of_nodes
                this%networks(i)%nodes(j)%v = this%circuit%getNodeVoltage(this%networks(i)%nodes(j)%name)
            end do
        end do

    end subroutine

    subroutine updateCircuitCurrentsFromNetwork(this)
        class(network_manager_t) :: this
        integer :: i, j
        do i = 1, size(this%networks)
            do j = 1, this%networks(i)%number_of_nodes
                call this%circuit%updateNodeCurrent(this%networks(i)%nodes(j)%name, this%networks(i)%nodes(j)%i)
            end do
        end do
    end subroutine

    subroutine network_advanceVoltage(this)
        class(network_manager_t) :: this
        call this%updateCircuitCurrentsFromNetwork()
        call this%circuit%step()
        ! this%circuit%time = this%circuit%time + this%circuit%dt
        call this%updateNetworkVoltagesFromCircuit()
    end subroutine

    subroutine updateNetworkVoltagesFromCircuit(this)
        class(network_manager_t) :: this
        integer :: i, j, idx
        type(vectorInfo_t), pointer :: info
        type(c_ptr) :: info_ptr
        real(kind=c_double), pointer :: values(:)
        type(string_t), allocatable :: names(:)

        do i = 1, size(this%networks)
            do j = 1, this%networks(i)%number_of_nodes
                info_ptr = get_vector_info(trim(this%networks(i)%nodes(j)%name)//c_null_char)
                if (.not. c_associated(info_ptr)) then
                    call WarnErrReport('Ngspice returned null vector info for '//trim(this%networks(i)%nodes(j)%name), .true.)
                    return
                end if

                call c_f_pointer(info_ptr, info)
                if (.not. c_associated(info%vRealData)) then
                    call WarnErrReport('Ngspice returned null vector data for '//trim(this%networks(i)%nodes(j)%name), .true.)
                    return
                end if
                if (info%vLength <= 0) then
                    call WarnErrReport('Ngspice returned empty vector for '//trim(this%networks(i)%nodes(j)%name), .true.)
                    return
                end if

                call c_f_pointer(info%vRealData, values,shape=[info%vLength])
                if (this%networks(i)%nodes(j)%name /= "time" .and. .not. this%networks(i)%nodes(j)%open) then 
                    this%networks(i)%nodes(j)%v = values(ubound(values,1))
                end if
            end do
        end do
    end subroutine
    
end module