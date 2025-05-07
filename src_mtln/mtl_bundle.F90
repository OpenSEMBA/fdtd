module mtl_bundle_mod

    use utils_mod
    use probes_mod
    use dispersive_mod
    use mtl_mod
#ifdef CompileWithMPI
    use FDETYPES, only: SUBCOMM_MPI, REALSIZE, MPI_STATUS_SIZE
#endif
!     use MPIcomm
     implicit none

    type, public :: mtl_bundle_t
        character (len=:), allocatable :: name
        real, allocatable, dimension(:,:,:) :: lpul, cpul, rpul, gpul
        integer  :: number_of_conductors = 0, number_of_divisions = 0
        real, dimension(:), allocatable :: step_size
        real, allocatable, dimension(:,:) :: v, i, e_L
        real, allocatable, dimension(:,:,:) :: du(:,:,:)
        real :: time = 0.0, dt = 1e10

        logical :: is_left_end = .true.
        logical :: is_right_end = .true.

        type(probe_t), allocatable, dimension(:) :: probes
        type(transfer_impedance_t) :: transfer_impedance
        integer, dimension(:), allocatable :: conductors_in_level
        
        real, dimension(:,:,:), allocatable :: v_term, i_term
        real, dimension(:,:,:), allocatable :: v_diff, i_diff

        type(external_field_segment_t), dimension(:), allocatable :: external_field_segments
        logical :: isPassthrough = .false.

    contains
        procedure :: mergePULMatrices
        procedure :: mergeDispersiveMatrices
        procedure :: initialAllocation
        procedure :: addProbe
        procedure :: updateLRTerms
        procedure :: updateCGTerms

        procedure :: updateSources => bundle_updateSources
        procedure :: advanceVoltage => bundle_advanceVoltage
        procedure :: advanceCurrent => bundle_advanceCurrent
        procedure :: addTransferImpedance => bundle_addTransferImpedance
        procedure :: setConnectorTransferImpedance => bundle_setConnectorTransferImpedance
        procedure :: setExternalLongitudinalField => bundle_setExternalLongitudinalField
#ifdef CompileWithMPI
        procedure :: Comm_MPI_MTL_V
#endif


    end type mtl_bundle_t

    interface mtl_bundle_t
        module procedure mtldCtor
    end interface

contains

    function mtldCtor(levels, name) result(res)
        type(mtl_bundle_t) :: res
        type(mtl_array_t), dimension(:), intent(in) :: levels
        character(len=*), intent(in), optional :: name
        
        res%name = ""
        if (present(name)) then
            res%name = name
        endif   
        allocate(res%probes(0))

        res%number_of_conductors = countNumberOfConductors(levels)
        res%dt = levels(1)%lines(1)%dt

        res%is_left_end = levels(1)%lines(1)%is_left_end
        res%is_right_end = levels(1)%lines(1)%is_right_end

        res%step_size = levels(1)%lines(1)%step_size
        res%number_of_divisions = size(res%step_size,1)
        res%external_field_segments = levels(1)%lines(1)%external_field_segments
        res%isPassthrough = levels(1)%lines(1)%isPassthrough
        call res%initialAllocation()
        call res%mergePULMatrices(levels)
        call res%mergeDispersiveMatrices(levels)


    end function

    subroutine initialAllocation(this)
        class(mtl_bundle_t) :: this
        allocate(this%lpul(this%number_of_divisions, this%number_of_conductors, this%number_of_conductors), source = 0.0)
        allocate(this%cpul(this%number_of_divisions + 1, this%number_of_conductors, this%number_of_conductors), source = 0.0)
        allocate(this%gpul(this%number_of_divisions + 1, this%number_of_conductors, this%number_of_conductors), source = 0.0)
        allocate(this%rpul(this%number_of_divisions, this%number_of_conductors, this%number_of_conductors), source = 0.0)
        allocate(this%du(this%number_of_divisions, this%number_of_conductors, this%number_of_conductors), source = 0.0)
        
        allocate(this%v(this%number_of_conductors, this%number_of_divisions + 1), source = 0.0)
        allocate(this%i(this%number_of_conductors, this%number_of_divisions), source = 0.0)
        allocate(this%e_L(this%number_of_conductors, this%number_of_divisions), source = 0.0)
        
        allocate(this%i_term(this%number_of_divisions,this%number_of_conductors,this%number_of_conductors), source = 0.0)
        allocate(this%v_diff(this%number_of_divisions,this%number_of_conductors,this%number_of_conductors), source = 0.0)

        allocate(this%v_term(this%number_of_divisions + 1,this%number_of_conductors,this%number_of_conductors), source = 0.0)
        allocate(this%i_diff(this%number_of_divisions + 1,this%number_of_conductors,this%number_of_conductors), source = 0.0)

    end subroutine

    function countNumberOfConductors(levels) result(res)
        type(mtl_array_t), dimension(:), intent(in) :: levels
        integer :: i,j, res
        res = 0
        do i = 1, size(levels)
            do j = 1, size(levels(i)%lines)
                res = res + levels(i)%lines(j)%number_of_conductors
            end do
        end do  
    end function

    subroutine mergePULMatrices(this, levels)
        class(mtl_bundle_t) :: this
        type(mtl_array_t), dimension(:), intent(in) :: levels
        integer :: i, j, n, n_sum
        n_sum = 0
        do i = 1, size(levels)
            do j = 1, size(levels(i)%lines)
                n = levels(i)%lines(j)%number_of_conductors
                this%lpul(:, n_sum + 1: n_sum+n , n_sum +1 : n_sum+n) = levels(i)%lines(j)%lpul(:,:,:)
                this%cpul(:, n_sum + 1: n_sum+n , n_sum +1 : n_sum+n) = levels(i)%lines(j)%cpul(:,:,:)
                this%rpul(:, n_sum + 1: n_sum+n , n_sum +1 : n_sum+n) = levels(i)%lines(j)%rpul(:,:,:)
                this%gpul(:, n_sum + 1: n_sum+n , n_sum +1 : n_sum+n) = levels(i)%lines(j)%gpul(:,:,:)
                this%du(:, n_sum + 1: n_sum+n , n_sum +1 : n_sum+n) = levels(i)%lines(j)%du(:,:,:)
                n_sum = n_sum+n
            end do
        end do
    end subroutine

    subroutine mergeDispersiveMatrices(this, levels)
        class(mtl_bundle_t) :: this
        type(mtl_array_t), dimension(:), intent(in) :: levels
        integer :: i, j, n, n_sum, number_of_poles
        n_sum = 0
        number_of_poles = 0
        do i = 1, size(levels)
            do j = 1, size(levels(i)%lines)
                number_of_poles = max(number_of_poles, levels(i)%lines(j)%lumped_elements%number_of_poles)
            end do
        end do
        this%transfer_impedance = &
            transfer_impedance_t(this%number_of_conductors, number_of_poles, this%number_of_divisions, this%dt)
        do i = 1, size(levels)
            do j = 1, size(levels(i)%lines)
                n = levels(i)%lines(j)%number_of_conductors

                this%transfer_impedance%q1(:,n_sum+1:n_sum+n,n_sum +1:n_sum+n,:) = &
                    levels(i)%lines(j)%lumped_elements%q1(:,:,:,:)
                
                this%transfer_impedance%q2(:,n_sum+1:n_sum+n,n_sum +1:n_sum+n,:) = &
                    levels(i)%lines(j)%lumped_elements%q2(:,:,:,:)
                
                this%transfer_impedance%q3(:,n_sum+1:n_sum+n,n_sum +1:n_sum+n,:) = &
                    levels(i)%lines(j)%lumped_elements%q3(:,:,:,:)
                
                this%transfer_impedance%q1_sum(:,n_sum+1:n_sum+n,n_sum +1:n_sum+n) = & 
                    levels(i)%lines(j)%lumped_elements%q1_sum(:,:,:)
                
                this%transfer_impedance%q2_sum(:,n_sum+1:n_sum+n,n_sum +1:n_sum+n) = & 
                    levels(i)%lines(j)%lumped_elements%q2_sum(:,:,:)
                
                this%transfer_impedance%q3_phi(:,n_sum+1:n_sum+n) = & 
                    levels(i)%lines(j)%lumped_elements%q3_phi(:,:)
                
                this%transfer_impedance%phi(:,n_sum+1:n_sum+n,:)  = & 
                    levels(i)%lines(j)%lumped_elements%phi(:,:,:)
                
                this%transfer_impedance%d(:,n_sum+1:n_sum+n,n_sum +1:n_sum+n) = & 
                    levels(i)%lines(j)%lumped_elements%d(:,:,:)
                
                this%transfer_impedance%e(:,n_sum+1:n_sum+n,n_sum +1:n_sum+n) = & 
                    levels(i)%lines(j)%lumped_elements%e(:,:,:)

                n_sum = n_sum + n
            end do
        end do

    end subroutine

    type(probe_t) function addProbe(this, index, probe_type, name, position) result(res)
        class(mtl_bundle_t) :: this
        integer, intent(in) :: index
        integer, intent(in) :: probe_type
        real, dimension(3), optional :: position
        character (len=:), allocatable, optional :: name
        type(probe_t), allocatable, dimension(:) :: aux_probes

        aux_probes = this%probes
        deallocate(this%probes)
        allocate(this%probes(size(aux_probes)+1))

        if (present(position) .and. present(name)) then
            res = probeCtor(index, probe_type, this%dt, name, position)
        else if (present(name) .and. .not. present(position)) then 
            res = probeCtor(index, probe_type, this%dt, name=name)
        else if (.not. present(name) .and. present(position)) then 
            res = probeCtor(index, probe_type, this%dt, position=position)
        else
            res = probeCtor(index, probe_type, this%dt)
        end if

        this%probes(1:size(this%probes)-1) = aux_probes
        this%probes(size(aux_probes)+1) = res
    end function

    subroutine bundle_setConnectorTransferImpedance(this, index, conductor_out, range_in, transfer_impedance)
        class(mtl_bundle_t) :: this
        integer, intent(in) :: index
        integer, intent(in) :: conductor_out
        integer, dimension(:), intent(in) :: range_in
        type(transfer_impedance_per_meter_t) :: transfer_impedance

        call this%transfer_impedance%setTransferImpedance(index, conductor_out, range_in, transfer_impedance)

    end subroutine

    subroutine bundle_addTransferImpedance(this, conductor_out, range_in, transfer_impedance)
        class(mtl_bundle_t) :: this
        integer, intent(in) :: conductor_out
        integer, dimension(:), intent(in) :: range_in
        type(transfer_impedance_per_meter_t) :: transfer_impedance

        call this%transfer_impedance%addTransferImpedance(conductor_out, range_in, transfer_impedance)

    end subroutine

    subroutine updateLRTerms(this)
        class(mtl_bundle_t) ::this
        real, dimension(this%number_of_divisions,this%number_of_conductors,this%number_of_conductors) :: F1, F2, IF1
        integer :: i

        F1 = reshape(source=[(matmul( &
            this%du(i,:,:), &
            this%lpul(i,:,:)/this%dt + &
                0.5*this%transfer_impedance%d(i,:,:) + &
                this%transfer_impedance%e(i,:,:)/this%dt + &
                0.5*this%rpul(i,:,:) + &
                this%transfer_impedance%q1_sum(i,:,:)), &
            i = 1,this%number_of_divisions)], & 
            shape=[this%number_of_divisions,this%number_of_conductors, this%number_of_conductors], &
            order=[2,3,1])

        F2 = reshape(source=[(matmul( &
            this%du(i,:,:), &
            this%lpul(i,:,:)/this%dt - &
            0.5*this%transfer_impedance%d(i,:,:) + &
            this%transfer_impedance%e(i,:,:)/this%dt - &
            0.5*this%rpul(i,:,:) - &
            this%transfer_impedance%q1_sum(i,:,:)), &
            i = 1,this%number_of_divisions)], & 
            shape=[this%number_of_divisions,this%number_of_conductors, this%number_of_conductors], &
            order=[2,3,1])

        IF1 = reshape(source=[(inv(F1(i,:,:)), i = 1, this%number_of_divisions)], &
                    shape=[this%number_of_divisions,this%number_of_conductors, this%number_of_conductors], &
                    order=[2,3,1])
        this%i_term = reshape(&
            source=[(matmul(IF1(i,:,:), F2(i,:,:)), i = 1, this%number_of_divisions)], &
            shape=[this%number_of_divisions,this%number_of_conductors, this%number_of_conductors], &
            order=[2,3,1])
        this%v_diff = IF1

    end subroutine



    subroutine updateCGTerms(this)
        class(mtl_bundle_t) ::this
        real, dimension(this%number_of_divisions + 1,this%number_of_conductors,this%number_of_conductors) :: F1, F2, IF1
        real, dimension(this%number_of_divisions + 1, this%number_of_conductors,this%number_of_conductors) :: extended_du
        integer :: i
        
        extended_du(1,:,:) = this%du(1,:,:)
        do i = 2, this%number_of_divisions
            extended_du(i,:,:)= 0.5*(this%du(i,:,:)+this%du(i-1,:,:))
        end do
        extended_du(this%number_of_divisions + 1,:,:) = this%du(this%number_of_divisions,:,:)

        F1 = reshape(&
            source=[(matmul(extended_du(i,:,:), &
            this%cpul(i,:,:)/this%dt) + 0.5*this%gpul(i,:,:), i = 1, this%number_of_divisions + 1)], &
            shape=[this%number_of_divisions + 1,this%number_of_conductors, this%number_of_conductors], &
            order=[2,3,1])
        F2 = reshape(&
            source=[(matmul(extended_du(i,:,:), &
            this%cpul(i,:,:)/this%dt) - 0.5*this%gpul(i,:,:), i = 1, this%number_of_divisions + 1)], &
            shape=[this%number_of_divisions + 1,this%number_of_conductors, this%number_of_conductors], &
            order=[2,3,1])

        IF1 = reshape(&
            source=[(inv(F1(i,:,:)), i = 1, this%number_of_divisions + 1)], &
            shape=[this%number_of_divisions + 1,this%number_of_conductors, this%number_of_conductors], &
            order=[2,3,1])

        this%v_term = reshape(&
            source=[(matmul(IF1(i,:,:), F2(i,:,:)), i = 1, this%number_of_divisions + 1)], &
            shape=[this%number_of_divisions + 1,this%number_of_conductors, this%number_of_conductors], &
            order=[2,3,1])
        this%i_diff = IF1

        do i = 2, this%number_of_divisions
            this%i_diff(i,1,1) = this%i_diff(i,1,1)/this%external_field_segments(i)%dielectric%effective_relative_permittivity
        end do
         
    end subroutine

    subroutine bundle_updateSources(this, time, dt)
        class(mtl_bundle_t) ::this
        real, intent(in) :: time, dt
        !TODO
    end subroutine

    subroutine bundle_advanceVoltage(this)
        class(mtl_bundle_t) ::this
        integer :: i


        do i = 2,this%number_of_divisions
            this%v(:, i) = matmul(this%v_term(i,:,:), this%v(:,i)) - &
                           matmul(this%i_diff(i,:,:), this%i(:,i) - this%i(:,i-1)  )
        end do
    end subroutine


    subroutine bundle_advanceCurrent(this)
        class(mtl_bundle_t) ::this
        real, dimension(:,:), allocatable :: i_prev, i_now
        integer :: i
        real :: eps_r

#ifdef CompileWithMPI
        call this%Comm_MPI_MTL_V()
#endif
        ! write(*,*)'advance'
        call this%transfer_impedance%updateQ3Phi()
        i_prev = this%i

        do i = 1, this%number_of_divisions 
            this%i(:,i) = matmul(this%i_term(i,:,:), this%i(:,i)) - &
                        matmul(this%v_diff(i,:,:), (this%v(:,i+1) - this%v(:,i)) - &
                                                    this%e_L(:,i) * this%step_size(i)) - &
                        matmul(this%v_diff(i,:,:), matmul(this%du(i,:,:), this%transfer_impedance%q3_phi(i,:)))
        enddo
        !TODO - revisar
        i_now = this%i
        call this%transfer_impedance%updatePhi(i_prev, i_now)
    end subroutine

    subroutine bundle_setExternalLongitudinalField(this)
        class(mtl_bundle_t) :: this
        integer :: i, j

        if (this%isPassthrough) then 
            do j = 2, 1 + this%conductors_in_level(2)
                do i = 1, size(this%e_L,2)
                    this%e_L(j,i) = this%external_field_segments(i)%field * &
                                    this%external_field_segments(i)%direction/abs(this%external_field_segments(i)%direction)
                end do
            end do
        else
            do i = 1, size(this%e_L,2)
                this%e_L(1,i) = this%external_field_segments(i)%field * &
                                this%external_field_segments(i)%direction/abs(this%external_field_segments(i)%direction)
            end do
        end if
    end subroutine

#ifdef CompileWithMPI
    subroutine Comm_MPI_MTL_V(this)
        class(mtl_bundle_t) :: this
        integer :: number_of_divisions, number_of_conductors, i
        integer :: ierr, rank, sizeof,status(MPI_STATUS_SIZE)

        ! call mpi_barrier(subcomm_mpi,ierr)
        call MPI_COMM_RANK(SUBCOMM_MPI, rank, ierr)
        call MPI_COMM_SIZE(SUBCOMM_MPI, sizeof, ierr)

        number_of_conductors = size(this%v,1)
        number_of_divisions = size(this%v,2)-1

        if (.not. this%is_left_end) then 
            do i = 1, number_of_conductors
                ! MPI_Send(buf, count, datatype, dest, tag, comm, ierror)
                ! MPI_Recv(buf, count, datatype, source, tag, comm, status, ierror)
                call MPI_send(this%v(i,3), 1, REALSIZE, rank-1, 1000*(rank-1)+i,  SUBCOMM_MPI, ierr)
                call MPI_recv(this%v(i,1), 1, REALSIZE, rank-1, 100*rank+i,  SUBCOMM_MPI, status, ierr)
            end do
        end if
        if (.not. this%is_right_end) then 
            do i = 1, number_of_conductors
                ! MPI_Send(buf, count, datatype, dest, tag, comm, ierror)
                ! MPI_Recv(buf, count, datatype, source, tag, comm, status, ierror)
                call MPI_send(this%v(i,number_of_divisions - 1), 1, REALSIZE, rank+1, 100*(rank+1)+i, SUBCOMM_MPI, ierr)
                call MPI_recv(this%v(i,number_of_divisions + 1), 1, REALSIZE, rank+1, 1000*rank+i, SUBCOMM_MPI, status, ierr)
            end do
        end if

        ! call mpi_barrier(subcomm_mpi,ierr)

    end subroutine
#endif
end module mtl_bundle_mod