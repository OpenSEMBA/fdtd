module network_manager_m

    use network_m
    use circuit_m
    use termination_handler_m
    use mtln_types_m, only: node_source_t
    use FDETYPES_m, only: RKIND, RKIND_TIEMPO

    implicit none 

    type network_manager_t
        type(network_t), dimension(:), allocatable :: networks
        type(circuit_t) :: circuit
        type(simple_termination_t), allocatable :: terminations(:)
        integer, allocatable :: ngspice_node_indices(:)
        integer :: num_simple
        integer :: num_ngspice
        real(kind=rkind) :: time, dt
    contains
        procedure :: advanceVoltage => network_advanceVoltage
        procedure :: updateCircuitCurrentsFromNetwork
        procedure :: updateNetworkVoltages
        procedure :: initTerminations

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
        call res%circuit%init(copy_node_names(networks), copy_sources(networks))
        res%circuit%dt = dt
#ifdef CompileWithRelease
        printInput = .false.
#endif        
        call res%circuit%readInput(description, printInput)
        call res%circuit%setModStopTimes(dt)
        call res%initTerminations()

    end function

    subroutine initTerminations(this)
        class(network_manager_t) :: this
        integer :: i, j, n_simple, n_ngspice
        integer :: ngspice_indices(10000)
        integer :: idx

        n_simple = 0
        n_ngspice = 0

        do i = 1, size(this%networks)
            do j = 1, this%networks(i)%number_of_nodes
                if (isSimpleTermination(this%networks(i)%nodes(j)%termination_type)) then
                    n_simple = n_simple + 1
                else
                    n_ngspice = n_ngspice + 1
                    ngspice_indices(n_ngspice) = (i - 1) * 1000 + j
                end if
            end do
        end do

        allocate(this%terminations(n_simple))
        allocate(this%ngspice_node_indices(n_ngspice))
        this%num_simple = n_simple
        this%num_ngspice = n_ngspice

        idx = 1
        do i = 1, size(this%networks)
            do j = 1, this%networks(i)%number_of_nodes
                if (isSimpleTermination(this%networks(i)%nodes(j)%termination_type)) then
                    call this%terminations(idx)%init( &
                         this%networks(i)%nodes(j)%termination_type, &
                         this%networks(i)%nodes(j)%R, &
                         this%networks(i)%nodes(j)%L, &
                         this%networks(i)%nodes(j)%C, &
                         this%dt)
                    idx = idx + 1
                end if
            end do
        end do

        this%ngspice_node_indices = ngspice_indices(1:n_ngspice)
    end subroutine initTerminations

    logical function isSimpleTermination(term_type)
        integer, intent(in) :: term_type
        isSimpleTermination = .false.
        select case (term_type)
        case (TERMINATION_SHORT, TERMINATION_OPEN, TERMINATION_SERIES, &
             TERMINATION_PARALLEL, TERMINATION_RsLCp, TERMINATION_RLsCp, &
             TERMINATION_LsRCp, TERMINATION_CsLRp, TERMINATION_RCsLp, &
             TERMINATION_LCsRp)
            isSimpleTermination = .true.
        end select
    end function isSimpleTermination

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
        integer :: i, j, idx

        ! Update simple terminations directly in Fortran
        idx = 1
        do i = 1, size(this%networks)
            do j = 1, this%networks(i)%number_of_nodes
                if (isSimpleTermination(this%networks(i)%nodes(j)%termination_type)) then
                    call this%terminations(idx)%step(this%networks(i)%nodes(j)%i, this%dt)
                    this%networks(i)%nodes(j)%v = this%terminations(idx)%v_node
                    idx = idx + 1
                end if
            end do
        end do

        ! Update complex terminations via ngspice
        if (this%num_ngspice > 0) then
            call this%updateCircuitCurrentsFromNetwork()
            call this%circuit%step()
            this%circuit%time = this%circuit%time + this%circuit%dt
            call this%updateNetworkVoltages()
        end if
    end subroutine

end module