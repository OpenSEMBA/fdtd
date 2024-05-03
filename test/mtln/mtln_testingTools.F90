module mtln_testingTools_mod
    use iso_c_binding
    use mtl_mod, only: mtl_t
    use network_mod
    use mtln_types_mod, only: terminal_node_t, termination_t
    implicit none
    
    character(len=*, kind=c_char), parameter :: PATH_TO_TEST_DATA = c_char_'./testData/'

contains
    
    

    type(network_t) function buildNetwork(name, r1, c1, r2, c2, target_v_nodes, target_i_nodes) result(res)
        character(len=*), intent(in) :: name, r1, c1, r2, c2
        character(100), dimension(:), allocatable :: description
        type(nw_node_t), dimension(3) :: nodes
        type(nw_node_t) :: node_int, node_out, node_in
        real, dimension(:),intent(in), target :: target_v_nodes, target_i_nodes

        ! node_int%values%v => target_v_nodes(1)
        ! node_int%values%i => target_i_nodes(1)
        ! node_in%values%v  => target_v_nodes(2)
        ! node_in%values%i  => target_i_nodes(2)
        ! node_out%values%v => target_v_nodes(3)
        ! node_out%values%i => target_i_nodes(3)

        ! node_int%name = name//"_int"
        ! node_int%values%v = 0.0
        ! node_int%values%i = 0.0
        ! node_int%line_c_per_meter = 0.0
    
        ! node_in%name = name//"_in"
        ! node_in%values%v = 0.0
        ! node_in%values%i = 0.0
        ! node_in%line_c_per_meter = 0.0
    
        ! node_out%name = name//"_out"
        ! node_out%values%v = 0.0
        ! node_out%values%i = 0.0
        ! node_out%line_c_per_meter = 0.0

        ! nodes = [node_in, node_int, node_out]
        ! allocate(description(0))
        ! description = [description, trim("R_"//name//"_1 "//name//"_int "//name//"_in "//r1)]
        ! description = [description, trim("V_"//name//"_1 "//name//"_in 0 dc 0 PULSE (0 5 1u 1u 1u 1 1)")]
        ! description = [description, trim("R_"//name//"_2 "//name//"_out "//name//"_int "//r2)]
        ! description = [description, trim("C_"//name//"_1 "//name//"_int 0 "//c1)]
        ! description = [description, trim("C_"//name//"_2 "//name//"_out 0 "//c2)]
        ! res = networkCtor(nodes, description)
    
    end function    

    type(mtl_t) function buildLineWithNConductors(n,name, parent_name, conductor_in_parent, dt) result(res)
    
        integer, intent(in) :: n
        character(len=*), intent(in) :: name
        real, intent(in), optional :: dt
        character(len=*), intent(in), optional :: parent_name
        integer, intent(in), optional :: conductor_in_parent
        real, allocatable, dimension(:,:) :: lpul, cpul, rpul, gpul
        real, dimension(5) :: step_size = [20.0, 20.0, 20.0, 20.0, 20.0]
        type(external_field_segment_t), dimension(5) :: external_field_segments
        integer :: i,j

        allocate(lpul(n,n), source = 0.0)
        allocate(cpul(n,n), source = 0.0)
        allocate(gpul(n,n), source = 0.0)
        allocate(rpul(n,n), source = 0.0)

        do i = 1, 5
            external_field_segments(i)%position =(/i,1,1/)            
            external_field_segments(i)%direction = DIRECTION_X_POS  
            external_field_segments(i)%Efield_main2wire => null()
            external_field_segments(i)%Efield_wire2main => null()
        end do

        do i = 1, n
            do j = 1, n
                rpul(i,j) = 0.0
                gpul(i,j) = 0.0
                if (i==j) then
                    lpul(i,j) = 4.4712610E-07
                    cpul(i,j) = 2.242e-10
                else 
                    lpul(i,j) = 1.4863653E-07
                    cpul(i,j) = -7.453e-11
                end if
            end do
        end do
        if (present(dt) .and. .not.present(parent_name)) then
            res = mtl_t(lpul, cpul, rpul, gpul, step_size, name, dt = dt, &
                        external_field_segments = external_field_segments)
        else if (.not.present(dt) .and. present(parent_name) ) then
            res = mtl_t(lpul, cpul, rpul, gpul, step_size, name, & 
                        parent_name= parent_name, &
                        conductor_in_parent=conductor_in_parent, &
                        external_field_segments = external_field_segments)
        else if (present(dt) .and. present(parent_name) ) then
            res = mtl_t(lpul, cpul, rpul, gpul, step_size, name, & 
                        parent_name= parent_name, &
                        conductor_in_parent=conductor_in_parent, &
                        dt = dt, external_field_segments = external_field_segments)
        else 
            res = mtl_t(lpul, cpul, rpul, gpul, step_size, name, &
                        external_field_segments = external_field_segments)
        end if
    end function    

    subroutine comparePULMatrices(error_cnt, m_line, m_input)
        integer, intent(inout) :: error_cnt
        real, intent(in), dimension(:,:,:) :: m_line
        real, intent(in), dimension(:,:) :: m_input
        integer :: i

        if (size(m_input, dim = 1) .ne. size(m_input, dim = 2)) then
            error_cnt = error_cnt + 1
            return
        end if   

        do i = 1, size(m_line, dim = 1)
            if (.not.ALL(m_line(i,:,:) == m_input(:,:))) then
                error_cnt = error_cnt + 1
            end if
        end do
        
    end subroutine 

    subroutine comparePULMatricesIH(error_cnt, m_line, m_input)
        integer, intent(inout) :: error_cnt
        real, intent(in), dimension(:,:,:) :: m_line
        real, intent(in), dimension(:,:,:) :: m_input
        integer :: i

        if (size(m_input, dim = 2) .ne. size(m_input, dim = 2)) then
            error_cnt = error_cnt + 1
            return
        end if   

        if (size(m_input, dim = 1) .ne. size(m_input, dim = 1)) then
            error_cnt = error_cnt + 1
            return
        end if   

        if (.not.ALL(m_line(:,:,:) == m_input(:,:,:))) then
            error_cnt = error_cnt + 1
        end if
        
    end subroutine 

    function checkNear_dp(target, number, rel_tol) result(is_near)
        double precision, intent(in) :: target, number
        real :: rel_tol
        logical :: is_near
        double precision :: abs_diff

        abs_diff = abs(target-number)
        if (abs_diff == 0.0) then
            is_near = .true.
        else 
            is_near = abs(target-number)/target < rel_tol
        endif

    end function 

    function checkNear(target, number, rel_tol) result(is_near)
        real, intent(in) :: target, number
        real :: rel_tol
        logical :: is_near
        real :: abs_diff

        abs_diff = abs(target-number)
        if (abs_diff == 0.0) then
            is_near = .true.
        else 
            is_near = abs(target-number)/target < rel_tol
        endif

    end function 

    function checkNear_real8(target, number, rel_tol) result(is_near)
        real(kind=8), intent(in) :: target, number
        real(kind=8) :: rel_tol
        logical :: is_near
        real(kind=8) :: abs_diff

        abs_diff = abs(target-number)
        if (abs_diff == 0.0) then
            is_near = .true.
        else 
            is_near = abs(target-number)/target < rel_tol
        endif

    end function 


 end module mtln_testingTools_mod
 