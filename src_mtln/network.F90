module network_mod

    use mtl_bundle_mod
    use mtln_types_mod, parsed_probe_t => probe_t, parsed_mtln_t => mtln_t
    use circuit_mod, only: string_t
    implicit none

    type nw_node_t
        character(len=:), allocatable :: name
        type(node_source_t) :: source
        integer :: source_type
        real :: line_c_per_meter, line_g_per_meter
        real :: step
        real :: v
        real :: i
        integer :: bundle_number, conductor_number, v_index, i_index
        integer :: side
    end type


    type, public :: network_t
        integer :: number_of_nodes = 0
        type(nw_node_t), dimension(:), allocatable :: nodes
        character(256), dimension(:), allocatable :: description
    contains

    end type network_t

    interface network_t
        module procedure networkCtor
    end interface


contains

    integer function countNodes(connections)
        type(terminal_connection_t), dimension(:), allocatable :: connections
        integer :: i
        countNodes = 0
        do i = 1, size(connections)
            countNodes = countNodes + size(connections(i)%nodes)
        end do
    end function

    function networkCtor(nodes, description) result(res)
        type(nw_node_t), dimension(:), intent(in) :: nodes
        character(*), dimension(:), intent(in) :: description
        type(string_t), dimension(:), allocatable :: names
        integer :: i
        type(network_t) :: res

        res%nodes = nodes
        res%description = description
        res%number_of_nodes = size(nodes)

        allocate(names(size(nodes)))
        do i = 1, size(nodes)
            names(i)%name = nodes(i)%name
            names(i)%length = len(nodes(i)%name)
        end do

    end function



end module