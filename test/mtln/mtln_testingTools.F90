module mtln_testingTools_mod
    use iso_c_binding
    use mtl_mod, only: mtl_t
    use network_mod
    use mtln_types_mod, only: terminal_node_t, termination_t
    implicit none
    
    character(len=*, kind=c_char), parameter :: PATH_TO_TEST_DATA = c_char_'./testData/'
    character(len=*), parameter :: MTL_TYPE_SHIELDED = "shielded"
    character(len=*), parameter :: MTL_TYPE_UNSHIELDED = "unshielded"

contains
    
    
    type(mtl_t) function buildLineWithNConductors(n,name, parent_name, conductor_in_parent, dt, type) result(res)
    
        integer, intent(in) :: n
        character(len=*), intent(in) :: name
        real, intent(in), optional :: dt
        character(len=*), intent(in), optional :: parent_name
        integer, intent(in), optional :: conductor_in_parent
        character(len=*), intent(in) :: type
        
        real, allocatable, dimension(:,:) :: lpul, cpul, rpul, gpul
        real, dimension(5) :: step_size = [20.0, 20.0, 20.0, 20.0, 20.0]
        type(segment_t), allocatable, dimension(:) :: segments
        integer :: i,j
        
        real :: time_step
        integer :: conductor
        character(len=:), allocatable :: parent

        type(transfer_impedance_per_meter_t):: Zt
        type(multipolar_expansion_t), dimension(:), allocatable:: mE

        Zt%inductive_term = 0.0
        Zt%resistive_term = 0.0
        allocate(Zt%poles(0), Zt%residues(0))
        allocate(mE(0))

        allocate(lpul(n,n), source = 0.0)
        allocate(cpul(n,n), source = 0.0)
        allocate(gpul(n,n), source = 0.0)
        allocate(rpul(n,n), source = 0.0)
        allocate(segments(5))
        do i = 1, 5
            segments(i)%x = 1
            segments(i)%y = i
            segments(i)%z = i
            segments(i)%orientation = DIRECTION_X_POS  
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
        if (.not. present(dt)) then 
            time_step = 1e-12
        else 
            time_step = dt
        end if
        
        if (type == MTL_TYPE_SHIELDED) then 
            if (.not. present(parent_name)) then 
                parent = "p" 
            else
                parent = parent_name
            end if
            if (.not. present(conductor_in_parent)) then 
                conductor = -1
            else 
                conductor = conductor_in_parent
            end if
            res = mtl_shielded(lpul, cpul, rpul, gpul, step_size, name, segments, time_step, parent, conductor, Zt)
        else if (type == MTL_TYPE_UNSHIELDED) then 
            res = mtl_unshielded(lpul, cpul, rpul, gpul, step_size, name, segments, time_step, mE)
        else
            write(*,*) 'Unrecognized line type'
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
 