module network_manager_mod

    use network_mod
    use circuit_mod
    implicit none 

    type network_manager_t
        type(network_t), dimension(:), allocatable :: networks
        type(circuit_t) :: circuit
        real :: time, dt
    contains
        procedure :: advanceVoltage => network_advanceVoltage
        procedure :: updateCircuitCurrentsFromNetwork
        procedure :: updateNetworkVoltages

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


    function copy_sources_names(networks) result(res)
        type(network_t), dimension(:), intent(in) :: networks
        type(string_t), dimension(:), allocatable :: res
        integer :: i,j
        type(string_t) :: temp
        allocate(res(0))
        do i = 1, size(networks)
            do j = 1, size(networks(i)%nodes)
                temp = string_t(trim(networks(i)%nodes(j)%source), len(trim(networks(i)%nodes(j)%source)))
                call appendToString_tArray(res, temp)
                ! res = [res, string_t(networks(i)%nodes(j)%source, len(networks(i)%nodes(j)%source))]
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
                ! res = [res, string_t(networks(i)%nodes(j)%name, len(networks(i)%nodes(j)%name))]
            end do
        end do
        call appendToString_tArray(res, string_t("time",4))
        ! res = [res, string_t("time",4)]
    end function


    function network_managerCtor(networks, description, final_time, dt) result(res)
        type(network_t), dimension(:), intent(in) :: networks
        character(*), dimension(:), intent(in) :: description
        real, intent(in) :: final_time, dt
        type(network_manager_t) :: res

        res%dt = dt
        res%time = 0.0
        res%networks = networks
        call res%circuit%init(copy_node_names(networks), copy_sources_names(networks))
        res%circuit%dt = dt
        call res%circuit%readInput(description)
        call res%circuit%setStopTimes(final_time, dt)

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
        this%circuit%time = this%circuit%time + this%circuit%dt
        call this%updateNetworkVoltages()
    end subroutine

end module