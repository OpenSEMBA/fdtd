module preprocess_mod

    use mtln_types_mod, parsed_probe_t => probe_t, parsed_mtln_t => mtln_t
    use mtl_bundle_mod
    use network_manager_mod
    use mtl_mod!, only: mtl_t, mtl_array_t, line_bundle_t,

    use fhash, only: fhash_tbl_t, key=>fhash_key, fhash_key_t

    
    implicit none


    type, public :: preprocess_t
        type(mtl_bundle_t), dimension(:), allocatable :: bundles
        type(network_manager_t) :: network_manager
        type(probe_t), dimension(:), allocatable :: probes
        type(fhash_tbl_t) :: conductors_before_cable
        type(fhash_tbl_t) :: cable_name_to_bundle_id
        real :: final_time, dt
    
    contains
        procedure :: buildMTLBundles
        procedure :: buildNetworkManager
        procedure :: buildNetwork
        procedure :: connectNodeToGround
        procedure :: connectNodes
        procedure :: addNodeWithId
        procedure :: addProbesWithId
    end type

    type, public :: cable_ptr_t
        type(cable_t), pointer :: p
    end type

    type, public :: cable_array_t
        type(cable_ptr_t), dimension(:), allocatable :: cables
    end type

    type, public :: cable_bundle_t
        type(cable_array_t), dimension(:), allocatable :: levels
    end type

    interface preprocess_t
        module procedure preprocess
    end interface

contains


    function preprocess(parsed) result(res)
        type(parsed_mtln_t), intent(in):: parsed
        type(preprocess_t) :: res
        type(fhash_tbl_t) :: cable_name_to_bundle_id
        type(line_bundle_t), dimension(:), allocatable :: line_bundles
        type(cable_bundle_t), dimension(:), allocatable :: cable_bundles


        res%final_time = parsed%time_step * parsed%number_of_steps
        res%dt = parsed%time_step

        cable_bundles = buildCableBundles(parsed%cables)
        line_bundles = buildLineBundles(cable_bundles)
        res%bundles = res%buildMTLBundles(line_bundles)
        res%cable_name_to_bundle_id = mapCablesToBundlesId(line_bundles, res%bundles)
        if (size(parsed%probes) /= 0) then
            res%probes = res%addProbesWithId(parsed%probes)
        end if
        res%network_manager = res%buildNetworkManager(parsed%networks)
        
    end function

    subroutine addConnector(line, connector, side)
        type(mtl_t), intent(inout) :: line
        type(connector_t) :: connector
        integer :: side

        integer :: i
        do i = 1, size(connector%resistances)
            line%rpul(side,i,i) = connector%resistances(i)
        end do

    end subroutine

    function conductorsInLevel(line) result(res)
        type(line_bundle_t), intent(in) :: line
        integer, dimension(:), allocatable :: res
        integer :: i,j

        allocate(res(size(line%levels)), source = 0)
        do i = 1, size(line%levels)
            do j = 1, size(line%levels(i)%lines)
                res(i) = res(i) + line%levels(i)%lines(j)%number_of_conductors
            end do
        end do
    end function

    function findConductorsBeforeCable(name, level) result(res)
        character(len=*), intent(in) :: name
        type(mtl_array_t), intent(in) :: level
        integer :: res 
        integer :: i
        res = 0
        do i = 1, size(level%lines)
            if (level%lines(i)%name /= name) then
                res = res + level%lines(i)%number_of_conductors
            else 
                return
            end if  
        end do
    end function

    function findOuterConductorNumber(line, level, conductors_in_level) result(res)
        type(mtl_t), intent(in) :: line
        type(mtl_array_t), intent(in) :: level
        integer, intent(in) :: conductors_in_level
        integer :: res
        res = findConductorsBeforeCable(line%parent_name, level) + &
              conductors_in_level + &
              line%conductor_in_parent
    end function

    function findInnerConductorRange(line, level, conductors_in_level) result(res)
        type(mtl_t), intent(in) :: line
        type(mtl_array_t), intent(in) :: level
        integer, intent(in) :: conductors_in_level
        integer, dimension(:), allocatable :: res
        integer :: k
        res = findConductorsBeforeCable(line%name, level) + & 
              conductors_in_level + &
              [(k, k = 1, line%number_of_conductors)]

    end function
    
    subroutine setBundleTransferImpedance(bundle, line)
        type(mtl_bundle_t), intent(inout) :: bundle
        type(line_bundle_t), intent(in) :: line
        integer :: i,j,k
        integer, dimension(:), allocatable :: range_in
        integer :: conductor_out
        type(transfer_impedance_per_meter_t) :: zt

        integer, dimension(:), allocatable :: conductors_in_level

        conductors_in_level = conductorsInLevel(line)
        bundle%conductors_in_level = conductors_in_level
        do i = 2, size(line%levels)
            do j = 1, size(line%levels(i)%lines)
                conductor_out = findOuterConductorNumber(line%levels(i)%lines(j), line%levels(i-1), sum(conductors_in_level(1:i-2)))
                range_in = findInnerConductorRange(line%levels(i)%lines(j), line%levels(i), sum(conductors_in_level(1:i-1)))
                call bundle%addTransferImpedance(conductor_out, range_in, line%levels(i)%lines(j)%transfer_impedance)
            end do
        end do  


    end subroutine

    subroutine mapConductorsBeforeCable(conductors_before_cable, line)
        type(fhash_tbl_t), intent(inout) :: conductors_before_cable
        type(line_bundle_t), intent(in) :: line
        integer, dimension(:), allocatable :: range_in
        integer, dimension(:), allocatable :: conductors_in_level
        integer :: i,j
        conductors_in_level = conductorsInLevel(line)
        call conductors_before_cable%set(key(line%levels(1)%lines(1)%name), 0)
        do i = 2, size(line%levels)
            do j = 1, size(line%levels(i)%lines)
                range_in = findInnerConductorRange(line%levels(i)%lines(j), line%levels(i), sum(conductors_in_level(1:i-1)))
                if (size(range_in) /= 0) then 
                    call conductors_before_cable%set(key(line%levels(i)%lines(j)%name), range_in(1) - 1)
                else
                    error stop 'range in cannot be empty'
                end if
            end do
        end do  

    end subroutine

    function buildMTLBundles(this, lines) result(res)
        class(preprocess_t) :: this
        type(line_bundle_t), dimension(:), intent(in) :: lines
        type(mtl_bundle_t), dimension(:), allocatable :: res
        type(fhash_tbl_t) :: conductors_before_cable
        integer :: i

        allocate(res(size(lines)))
        do i = 1, size(lines)
            res(i) = mtldCtor(lines(i)%levels, "bundle_"//lines(i)%levels(1)%lines(1)%name)
            call setBundleTransferImpedance(res(i), lines(i))
            call mapConductorsBeforeCable(conductors_before_cable, lines(i))
        end do  
        this%conductors_before_cable = conductors_before_cable
    end function    

    function buildLineFromCable(cable) result(res)
        type(cable_t), intent(in) :: cable
        type(mtl_t) :: res
        integer :: conductor_in_parent = 0
        character(len=:), allocatable :: parent_name
        if (associated(cable%parent_cable)) then 
            parent_name = cable%parent_cable%name
            conductor_in_parent = cable%conductor_in_parent
        end if  

        res = mtlHomogeneous(lpul = cable%inductance_per_meter, &
                             cpul = cable%capacitance_per_meter, &
                             rpul = cable%resistance_per_meter, &
                             gpul = cable%conductance_per_meter, &
                             step_size = cable%step_size, &
                             name = cable%name, &
                             parent_name = parent_name, &
                             conductor_in_parent = conductor_in_parent, & 
                             transfer_impedance = cable%transfer_impedance, &
                             external_field_segments = cable%external_field_segments)

        if (associated(cable%initial_connector)) call addConnector(res, cable%initial_connector, 0)
        if (associated(cable%end_connector))     call addConnector(res, cable%initial_connector, size(res%rpul,1))
                

    end function

    function buildLineBundles(cable_bundles) result(res)
        type(cable_bundle_t), dimension(:), allocatable :: cable_bundles
        type(line_bundle_t), dimension(:), allocatable :: res
        integer :: i, j, k
        integer :: nb, nl, nc
        nb = size(cable_bundles)

        allocate(res(nb))
        do i = 1, nb
            nl = size(cable_bundles(i)%levels)
            allocate(res(i)%levels(nl))
            do j = 1, nl
                nc = size(cable_bundles(i)%levels(j)%cables)
                allocate(res(i)%levels(j)%lines(nc))
                do k = 1, nc
                    res(i)%levels(j)%lines(k) = buildLineFromCable(cable_bundles(i)%levels(j)%cables(k)%p)
                end do
            end do
        end do

    end function

    function buildCableBundleFromParent(parent, cables) result(res)
        type(cable_ptr_t), intent(in) :: parent
        type(cable_t), dimension(:), intent(in), target :: cables
        type(cable_array_t) :: level
        type(cable_bundle_t) :: res

        allocate(res%levels(1))
        allocate(res%levels(1)%cables(1))

        allocate(level%cables(1))
        level%cables(1)%p => parent%p
        
        res%levels(1) = level

        do while (findNextLevel(level, cables) /= 0)
            res%levels = [res%levels, level]
        end do

        contains
            integer function findNextLevel(curr_level, c)
                type(cable_array_t), intent(inout) :: curr_level
                type(cable_t), dimension(:), intent(in), target :: c
                type(cable_t), target :: tgt
                type(cable_array_t) :: next_level
                integer :: i,j, next_level_size
                integer :: n
                ! allocate(next_level%cables(0))
                next_level_size = 0
                do i = 1, size(curr_level%cables) 
                    do j = 1, size(c)
                        if (associated(c(j)%parent_cable, curr_level%cables(i)%p)) then 
                            next_level_size = next_level_size + 1
                        end if
                    end do
                end do
                
                allocate(next_level%cables(next_level_size))
                n = 0
                do i = 1, size(curr_level%cables) 
                    ! tgt = curr_level%cables(i)%p
                    do j = 1, size(c)
                        if (associated(c(j)%parent_cable, curr_level%cables(i)%p)) then 
                        ! if (associated(c(j)%parent_cable, tgt)) then 
                            ! next_level%cables = [next_level%cables, c(j)]
                            n = n + 1
                            ! tgt = c(j)
                            next_level%cables(n)%p => c(j)
                        end if
                    end do
                end do
                curr_level = next_level
                findNextLevel = size(curr_level%cables)
            end function

    end function

    function findParentCables(cables) result(res)
        type(cable_t), dimension(:), intent(in), target :: cables
        type(cable_ptr_t), dimension(:), allocatable :: res
        integer :: i
        integer, dimension(:), allocatable :: parent_ids

        allocate(parent_ids(0))
        do i = 1, size(cables)
            if (associated(cables(i)%parent_cable) .eqv. .false.) then 
                parent_ids = [parent_ids, i]
            end if
        end do

        allocate(res(size(parent_ids)))
        do i = 1, size(parent_ids)
            res(i)%p => cables((parent_ids(i)))
        end do
    end function


    function buildCableBundles(cables) result(cable_bundles)
        type(cable_t), dimension(:), intent(in) :: cables
        type(cable_bundle_t), dimension(:), pointer :: cable_bundles
        type(cable_ptr_t), dimension(:), allocatable :: parents
        integer :: i

        parents = findParentCables(cables)
        allocate(cable_bundles(size(parents)))
        do i = 1, size(parents)
            cable_bundles(i) = buildCableBundleFromParent(parents(i), cables)
        end do

    end function

    function mapCablesToBundlesId(lines, bundles) result(res)
        type(line_bundle_t), dimension(:), allocatable :: lines
        type(mtl_bundle_t), dimension(:), allocatable :: bundles
        type(fhash_tbl_t) :: res
        integer :: i, j, k

        do i = 1, size(lines)
            do j = 1, size(lines(i)%levels)
                do k = 1, size(lines(i)%levels(j)%lines)
                    call res%set(key(lines(i)%levels(j)%lines(k)%name), value = i)
                end do
            end do
        end do

    end function

    function mapCablesToBundles(lines, bundles) result(res)
        type(line_bundle_t), dimension(:), allocatable :: lines
        type(mtl_bundle_t), dimension(:), allocatable :: bundles
        type(fhash_tbl_t) :: res
        integer :: i, j, k

        do i = 1, size(lines)
            do j = 1, size(lines(i)%levels)
                do k = 1, size(lines(i)%levels(j)%lines)
                    call res%set(key(lines(i)%levels(j)%lines(k)%name), value = bundles(i))
                end do
            end do
        end do

    end function

    function writeSeriesRLCnode(node, termination, end_node) result(res)
        type(nw_node_t), intent(in) :: node
        type(termination_t), intent(in) :: termination
        character(len=*), intent(in) :: end_node
        character(len=256), allocatable :: res(:)
        character(len=256) :: buff
        character(20) :: termination_r, termination_l, termination_c, line_c

        write(termination_c, *) termination%capacitance
        write(termination_r, *) termination%resistance
        write(termination_l, *) termination%inductance
        write(line_c, *) node%line_c_per_meter * node%step/2
        
        allocate(res(0))

        buff = trim("R" // node%name // " " // node%name // " "   // node%name //"_R " // termination_r)
        call appendToStringArray(res, buff)
        buff = trim("L" // node%name // " " // node%name // "_R " // node%name //"_L " // termination_l)
        call appendToStringArray(res, buff)
        if (termination%path_to_excitation /= "") then
            buff = trim("C" // node%name // " " // node%name // "_L " // node%name //"_V "// termination_c)
            call appendToStringArray(res, buff)
            buff = trim("V" // node%name // " " // node%name // "_V " // end_node //" dc 0" )
            call appendToStringArray(res, buff)
        else
            buff = trim("C" // node%name // " " // node%name // "_L " // end_node //" "// termination_c)
            call appendToStringArray(res, buff)
        end if
        buff = trim("I" // node%name // " " // node%name// " 0 " // " dc 0")
        call appendToStringArray(res, buff)
        buff = trim("CL" // node%name // " " // node%name // " 0 " // line_c)
        call appendToStringArray(res, buff)

    end function

    function writeSeriesRLnode(node, termination, end_node) result(res)
        type(nw_node_t), intent(in) :: node
        type(termination_t), intent(in) :: termination
        character(len=*), intent(in) :: end_node
        character(len=256), allocatable :: res(:)
        character(len=256) :: buff
        character(20) :: termination_r, termination_l, line_c

        write(termination_r, *) termination%resistance
        write(termination_l, *) termination%inductance
        write(line_c, *) node%line_c_per_meter * node%step/2

        allocate(res(0))

        res = [trim("R" // node%name // " " // node%name // "_R "   // node%name //" ")//" "//trim(termination_r)]
        if (termination%path_to_excitation /= "") then
            buff = trim("L" // node%name // " " // node%name // "_R " // node%name //"_L")//" "//trim(termination_l)
            call appendToStringArray(res, buff)
            buff = trim("V" // node%name // " " // node%name // "_L " // end_node //" dc 0" )
            call appendToStringArray(res, buff)
        else
            buff = trim("L" // node%name // " " // node%name // "_R " // end_node)//" "//trim(termination_l)
            call appendToStringArray(res, buff)
        end if
        buff = trim("I" // node%name // " " // node%name// " 0 " // " dc 0")
        call appendToStringArray(res, buff)
        buff = trim("CL" // node%name // " " // node%name // " 0 " // line_c)
        call appendToStringArray(res, buff)
        
    end function

    function writeRLsCpnode(node, termination, end_node) result(res)
        type(nw_node_t), intent(in) :: node
        type(termination_t), intent(in) :: termination
        character(len=*), intent(in) :: end_node
        character(len=256), allocatable :: res(:)
        character(len=256) :: buff
        character(20) :: termination_r, termination_l, termination_c, line_c

        write(termination_r, *) termination%resistance
        write(termination_c, *) termination%capacitance
        write(termination_l, *) termination%inductance
        write(line_c, *) node%line_c_per_meter * node%step/2

        allocate(res(0))

        buff = trim("R" // node%name // " " // node%name // " "   // node%name //"_R " // termination_r)
        call appendToStringArray(res, buff)
        if (termination%path_to_excitation /= "") then
            buff = trim("L" // node%name // " " // node%name // "_R " // node%name //"_V " // termination_l)
            call appendToStringArray(res, buff)
            buff = trim("C" // node%name // " " // node%name // " " // node%name //"_V " // termination_c)
            call appendToStringArray(res, buff)
            buff = trim("V" // node%name // " " // node%name // "_V " // end_node //" dc 0" )
            call appendToStringArray(res, buff)
        else 
            buff = trim("L" // node%name // " " // node%name // "_R " // end_node //" "// termination_l)
            call appendToStringArray(res, buff)
            buff = trim("C" // node%name // " " // node%name // " " // end_node //" "// termination_c)
            call appendToStringArray(res, buff)
        end if
        buff = trim("I" // node%name // " " // node%name// " 0 " // " dc 0")
        call appendToStringArray(res, buff)
        buff = trim("CL" // node%name // " " // node%name // " 0 " // line_c)
        call appendToStringArray(res, buff)


    end function

    function writeSeriesNode(node, termination, end_node) result(res)
        type(nw_node_t), intent(in) :: node
        type(termination_t), intent(in) :: termination
        character(len=*), intent(in) :: end_node
        character(len=256), allocatable :: res(:)

        if (termination%capacitance >= 1e22) then 
            res = writeSeriesRLnode(node, termination, end_node)
        else
            res = writeSeriesRLCnode(node, termination, end_node)
        end if

    end function

    subroutine appendToStringArray(arr, str)
        ! This has been implemented because there seems to be a bug in gfortran: 
        ! https://fortran-lang.discourse.group/t/read-data-and-append-it-to-array-best-practice/1915
        ! and arr = [ arr, str ] can't be used.
        character(len=256), allocatable, intent(inout) :: arr(:)
        character(len=256), intent(in) :: str
        character(len=256), allocatable :: old_arr(:)
        
        old_arr = arr
        deallocate(arr)
        allocate(arr(size(old_arr)+1))
        arr(1:size(old_arr)) = old_arr 
        arr(size(old_arr)+1) = str
    end subroutine

    function writeShortNode(node, termination, end_node) result(res)
        type(nw_node_t), intent(in) :: node
        type(termination_t), intent(in) :: termination
        character(len=*), intent(in) :: end_node
        character(len=256), allocatable :: res(:)
        character(len=256) :: buff
        character(20) :: short_R, line_c

        write(short_r, *) 1e-10
        write(line_c, *) node%line_c_per_meter*node%step/2

        allocate(res(0))
        if (termination%path_to_excitation /= "") then
            buff = trim("R" // node%name // " " // node%name // " " // node%name //"_R")//" "//trim(short_R)
            call appendToStringArray(res, buff)
            buff = trim("V" // node%name // " " // node%name // "_R " // end_node//" dc 0")
            call appendToStringArray(res, buff)
        else
            buff = trim("R" // node%name // " " // node%name // " " // end_node)//" "//trim(short_R)
            call appendToStringArray(res, buff)
        end if
        buff = trim("I" // node%name // " " // node%name// " 0 " // " dc 0")
        call appendToStringArray(res, buff)
        buff = trim("CL" // node%name // " " // node%name // " 0 " // line_c)
        call appendToStringArray(res, buff)
        
    end function

    function writeOpenNode(node, termination, end_node) result(res)
        type(nw_node_t), intent(in) :: node
        type(termination_t), intent(in) :: termination
        character(len=*), intent(in) :: end_node
        character(len=256), allocatable :: res(:)
        character(len=256) :: buff
        character(20) :: line_c

        write(line_c, *) node%line_c_per_meter*node%step/2

        allocate(res(0))
        buff = trim("I" // node%name // " " // node%name// " 0 " // " dc 0")
        call appendToStringArray(res, buff)
        buff = trim("CL" // node%name // " " // node%name // " 0 " // line_c)
        call appendToStringArray(res, buff)
        
    end function

    function writeLCpRsNode(node, termination, end_node) result(res)
        type(nw_node_t), intent(in) :: node
        type(termination_t), intent(in) :: termination
        character(len=*), intent(in) :: end_node
        character(len=256), allocatable :: res(:)
        character(len=256) :: buff
        character(len=:), allocatable :: node_name
        character(20) :: termination_r, termination_l, termination_c, line_c
        
        write(termination_r, *) termination%resistance
        write(termination_l, *) termination%inductance
        write(termination_c, *) termination%capacitance
        write(line_c, *) node%line_c_per_meter * node%step/2
       
        allocate(res(0))
        res = [trim("R" // node%name // " " // node%name // " "   // node%name //"_p " // termination_r)]
        if (termination%path_to_excitation /= "") then
            buff = trim("L" // node%name // " " // node%name // "_p " // node%name //"_V "// termination_l)
            call appendToStringArray(res, buff)
            buff = trim("C" // node%name // " " // node%name // "_p " // node%name //"_V "// termination_c)
            call appendToStringArray(res, buff)
            buff = trim("V" // node%name // " " // node%name // "_V " // end_node //" dc 0" )
            call appendToStringArray(res, buff)
        else
            buff =  trim("L" // node%name // " " // node%name // "_p " // end_node //" "// termination_l)
            call appendToStringArray(res, buff)
            buff =  trim("C" // node%name // " " // node%name // "_p " // end_node //" "// termination_c)
            call appendToStringArray(res, buff)
        end if
        buff =  trim("I" // node%name // " " // node%name// " 0 " // " dc 0")
        call appendToStringArray(res, buff)
        buff = trim("CL" // node%name // " " // node%name // " 0 " // line_c)
        call appendToStringArray(res, buff)

    end function

    function writeNodeDescription(node, termination, end_node) result(res)
        type(nw_node_t), intent(in) :: node
        type(termination_t), intent(in) :: termination
        character(len=256), allocatable :: res(:)
        character(len=*), intent(in) :: end_node

        if (termination%termination_type == TERMINATION_SERIES) then 
            res = writeSeriesNode(node, termination, end_node)
        else if (termination%termination_type == TERMINATION_LCpRs) then 
            res = writeLCpRsNode(node, termination, end_node)
        else if (termination%termination_type == TERMINATION_RLsCp) then 
            res = writeRLsCpNode(node, termination, end_node)
        else if (termination%termination_type == TERMINATION_SHORT) then 
            res = writeShortNode(node, termination , end_node)
        else if (termination%termination_type == TERMINATION_OPEN) then 
            res = writeOpenNode(node, termination , end_node)
        end if

    end function    

    function addNodeWithId(this, node) result(res)
        class(preprocess_t) :: this
        type(terminal_node_t) :: node
        integer :: stat
        type(mtl_bundle_t), target :: tbundle
        integer :: d
        type(nw_node_t) :: res
        character(len=4) :: sConductor
        integer :: conductor_number

        call this%conductors_before_cable%get(key(node%belongs_to_cable%name), conductor_number)
        conductor_number = conductor_number + node%conductor_in_cable
        
        call this%cable_name_to_bundle_id%get(key(node%belongs_to_cable%name), d, stat)
        if (stat /= 0) return
        tbundle = this%bundles(d)
        write(sConductor,'(I0)') node%conductor_in_cable
        res%name = trim(node%belongs_to_cable%name)//"_"//trim(sConductor)//"_"//nodeSideToString(node%side)
        write(*,*) res%name
        res%v = 0.0
        res%i = 0.0
        res%bundle_number = d
        res%conductor_number = conductor_number

        if (node%side == TERMINAL_NODE_SIDE_INI) then 
            res%v_index = lbound(tbundle%v,2)
            res%i_index = lbound(tbundle%i,2)
            res%line_c_per_meter = tbundle%cpul(lbound(tbundle%cpul,1), conductor_number, conductor_number)
            res%step = tbundle%du(lbound(tbundle%du,1), conductor_number, conductor_number)
            res%side = TERMINAL_NODE_SIDE_INI

        else if (node%side == TERMINAL_NODE_SIDE_END) then 
            res%v_index = ubound(tbundle%v,2)
            res%i_index = ubound(tbundle%i,2)
            res%line_c_per_meter = tbundle%cpul(ubound(tbundle%cpul,1), conductor_number, conductor_number)
            res%step = tbundle%du(ubound(tbundle%du,1), conductor_number, conductor_number)
            res%side = TERMINAL_NODE_SIDE_END
        end if
        
        res%source = node%termination%path_to_excitation

    contains
        function nodeSideToString(side) result(cSide)
            character (len=:), allocatable :: cSide
            integer, intent(in) :: side
            select case (side)
            case (TERMINAL_NODE_SIDE_INI)
                cSide = "initial"
            case (TERMINAL_NODE_SIDE_END)
                cSide = "end"
            end select
        end function

    end function

    subroutine connectNodeToGround(this, terminal_nodes, nodes, description)
        class(preprocess_t) :: this
        type(terminal_node_t), dimension(:), allocatable :: terminal_nodes
        type(nw_node_t),  dimension(:), allocatable, intent(inout) :: nodes
        character(256), dimension(:), allocatable, intent(inout) :: description
        character(256), dimension(:), allocatable :: node_description, old_description

        type(nw_node_t) :: new_node
        integer :: stat

        new_node = this%addNodeWithId(terminal_nodes(1))
        nodes = [nodes, new_node]

        node_description = writeNodeDescription(new_node, terminal_nodes(1)%termination, "0")
        old_description = description
        deallocate(description)
        allocate(description(size(old_description) + size(node_description)))
        description(1:size(old_description)) = old_description
        description((size(old_description)+1):size(description)) = node_description(:)
    end subroutine

    subroutine connectNodes(this, terminal_nodes, nodes, description)
        class(preprocess_t) :: this
        type(terminal_node_t), dimension(:), allocatable :: terminal_nodes
        type(nw_node_t),  dimension(:), intent(inout) :: nodes
        character(256), dimension(:), intent(inout) :: description
        type(nw_node_t) :: new_node
        integer :: i, stat
        character(len=:), allocatable :: interior_node

        interior_node = trim(terminal_nodes(1)%belongs_to_cable%name)//"_"//&
                        trim(terminal_nodes(2)%belongs_to_cable%name)//"_inter"
        do i = 1, 2
            new_node =this%addNodeWithId(terminal_nodes(i))
            nodes = [nodes, new_node]
            description = [description, writeNodeDescription(new_node, terminal_nodes(i)%termination, interior_node)]
        end do
    end subroutine


    function buildNetwork(this,terminal_network) result(res)
        class(preprocess_t) :: this
        type(terminal_network_t), intent(in) :: terminal_network
        type(nw_node_t), dimension(:), allocatable :: nodes
        character(256), dimension(:), allocatable :: description
        type(network_t) :: res
        integer :: i

        allocate(description(0))
        allocate(nodes(0))
        do i = 1, size(terminal_network%connections)
            if (size(terminal_network%connections(i)%nodes) == 1) then 
                call this%connectNodeToGround(terminal_network%connections(i)%nodes, nodes, description)
            else
                call this%connectNodes(terminal_network%connections(i)%nodes, nodes, description)
            end if
        end do

        res = networkCtor(nodes, description)
    end function

    subroutine endDescription(description)
        character(256), dimension(:), allocatable, intent(inout) :: description
        character(256) :: buff

        buff = ".end"
        call appendToStringArray(description, buff)
        
        buff = "NULL"
        call appendToStringArray(description, buff)
        
    end subroutine

    subroutine addNetworksDescription(description, networks)
        character(256), dimension(:), allocatable, intent(inout) :: description
        type(network_t), dimension(:), intent(in) :: networks
        integer :: i
        do i = 1, size(networks)
            description = [description, networks(i)%description]
        end do
    end subroutine

    subroutine addAnalysis(description, final_time, dt)
        character(256), dimension(:), allocatable, intent(inout) :: description
        character(256) :: buff
        real, intent(in) :: final_time, dt
        character(20) :: sTime, sdt
        
        write(sTime, '(E10.2)') final_time
        write(sdt, '(E10.2)') dt

        buff = trim(".tran "//sdt//" "//sTime//" 0 "//sdt)
        call appendToStringArray(description, buff)       

    end subroutine

    subroutine addSavedNodes(description, networks)
        character(256), dimension(:), allocatable, intent(inout) :: description
        character(256) :: buff
        type(network_t), dimension(:), intent(in) :: networks
        character(len=:), allocatable :: saved_nodes
        integer :: i,j
        saved_nodes = ".save "
        do j = 1, size(networks)
            do i = 1, size(networks(j)%nodes)
                saved_nodes = saved_nodes // "V1"//trim(networks(j)%nodes(i)%name)//"#branch "
                saved_nodes = saved_nodes // trim(networks(j)%nodes(i)%name) // " "
            end do
        end do
        
        buff = trim(saved_nodes)
        call appendToStringArray(description, buff)

    end subroutine


    function buildNetworkManager(this, terminal_networks) result(res)
        class(preprocess_t) :: this
        type(terminal_network_t), dimension(:), intent(in) :: terminal_networks
        type(network_t), dimension(:), allocatable :: networks
        type(network_manager_t) :: res
        character(256), dimension(:), allocatable :: description
        integer :: i

        allocate(networks(size(terminal_networks)))
        do i = 1, size(terminal_networks)
            networks(i) = this%buildNetwork(terminal_networks(i))
        end do
        
        allocate(description(0))
        description = ["* network description message"]
        call addNetworksDescription(description, networks)
        call addAnalysis(description, this%final_time, this%dt)
        call addSavedNodes(description, networks)
        call endDescription(description)        

        do i = 1, size(description)
            write(*,'(A)') trim(description(i))
        end do  
        res = network_managerCtor(networks, description, this%final_time, this%dt)

    end function


    function addProbesWithId(this, parsed_probes) result(res)
        class(preprocess_t) :: this
        type(parsed_probe_t), dimension(:), allocatable :: parsed_probes
        type(probe_t), dimension(:), allocatable :: res
        integer :: i, d
        integer :: stat
        type(mtl_bundle_t), target :: tbundle

        allocate(res(size(parsed_probes)))
        do i = 1, size(parsed_probes)
            call this%cable_name_to_bundle_id%get(key = key(parsed_probes(i)%attached_to_cable%name), &
                                               value = d, &
                                               stat=stat)

            if (stat /= 0) return
            res(i) =  this%bundles(d)%addProbe(index = parsed_probes(i)%index, probe_type = parsed_probes(i)%probe_type)
        end do
    end function

    
end module