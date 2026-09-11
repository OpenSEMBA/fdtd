module mtl_bundle_m

    use mtln_utils_m
    use probes_m, only: probe_t, probeCtor, PROBE_TYPE_CURRENT
    use generators_m
    use dispersive_m
    use mtl_m
#ifdef CompileWithMPI
    use FDETYPES_m, only: SUBCOMM_MPI, REALSIZE, INTEGERSIZE, MPI_STATUS_SIZE
    use mpi
#endif
    use mtln_types_m, only: SOURCE_TYPE_CURRENT, SOURCE_TYPE_VOLTAGE
    use FDETYPES_m, only: RKIND, RKIND_TIEMPO
    implicit none

    type, public :: mtl_bundle_t
        character(len=:), allocatable :: name
        real(kind=rkind), allocatable, dimension(:,:,:) :: lpul, cpul, rpul, gpul
        integer  :: number_of_conductors = 0, number_of_divisions = 0
        real(kind=RKIND), dimension(:), allocatable :: step_size
        real(kind=RKIND), pointer, dimension(:,:) :: v, i, i_prev
        real(kind=RKIND), allocatable, dimension(:,:) :: v_source, i_source, e_L
        real(kind=RKIND), allocatable, dimension(:,:,:) :: du(:,:,:)
        real(kind=RKIND_TIEMPO) :: time = 0.0, dt = 1e10

        ! type homogen
        type(generator_t), allocatable, dimension(:) :: generators

        type(probe_t), allocatable, dimension(:) :: probes
        type(transfer_impedance_t) :: transfer_impedance
        integer, dimension(:), allocatable :: conductors_in_level
        
        real(kind=rkind), dimension(:,:,:), allocatable :: v_term, i_term
        real(kind=rkind), dimension(:,:,:), allocatable :: v_diff, i_diff

        type(external_field_segment_t), dimension(:), allocatable :: external_field_segments
        logical :: bundle_in_layer = .true.

#ifdef CompileWithMPI
        integer(kind=4), allocatable, dimension(:,:) :: layer_indices
        type(comm_t) :: mpi_comm
        ! Owner rank of each physical current segment.
        integer(kind=4), allocatable :: segment_owner(:)
#endif

    contains
        procedure :: mergePULMatrices
        procedure :: mergeDispersiveMatrices
        procedure :: initialAllocation
        procedure :: addProbe
        procedure :: addGenerator
        procedure :: updateLRTerms
        procedure :: updateCGTerms
        procedure :: updateGenerators => bundle_updateGenerators
        procedure :: advanceVoltage => bundle_advanceVoltage
        procedure :: advanceCurrent => bundle_advanceCurrent
        procedure :: addTransferImpedance => bundle_addTransferImpedance
        procedure :: setConnectorTransferImpedance => bundle_setConnectorTransferImpedance
        procedure :: setExternalLongitudinalField => bundle_setExternalLongitudinalField
#ifdef CompileWithMPI
        procedure :: Comm_MPI_V
        procedure :: Comm_MPI_Fields
        procedure :: Comm_MPI_I
        procedure :: initMPIStateOwnership
#endif


    end type mtl_bundle_t

    interface mtl_bundle_t
        module procedure mtldCtor
    end interface

    type :: external_field_segment_t
        integer, dimension(3) ::position
        integer :: direction = 0
        real(kind=rkind) , pointer  :: field => null()
        ! True only on the rank that updates this FDTD field cell.
        logical :: field_is_local = .false.
        ! True on the geometric FDTD owner even when the cell is PEC/lossy.
        logical :: field_is_owned = .false.
    end type

contains

    function mtldCtor(levels, name) result(res)
        type(mtl_bundle_t) :: res
        type(transmission_line_level_t), dimension(:), intent(in) :: levels
        character(len=*), intent(in), optional :: name

        res%name = ""
        if (present(name)) then
            res%name = name
        end if
        allocate(res%probes(0))
        allocate(res%generators(0))

        res%number_of_conductors = countNumberOfConductors(levels)
        res%dt = levels(1)%lines(1)%dt

        res%step_size = levels(1)%lines(1)%step_size
        res%number_of_divisions = size(res%step_size,1)
        res%external_field_segments = buildExternalFieldSegments(levels)
        call res%initialAllocation()

        call res%mergePULMatrices(levels)
        call res%mergeDispersiveMatrices(levels)

#ifdef CompileWithMPI
        res%bundle_in_layer = levels(1)%lines(1)%bundle_in_layer
        res%layer_indices = levels(1)%lines(1)%layer_indices
        res%mpi_comm = levels(1)%lines(1)%mpi_comm
#endif

    end function

    subroutine initialAllocation(this)
        class(mtl_bundle_t) :: this
        allocate(this%lpul(this%number_of_divisions, this%number_of_conductors, this%number_of_conductors), source = 0.0_rkind)
        allocate(this%cpul(this%number_of_divisions + 1, this%number_of_conductors, this%number_of_conductors), source = 0.0_rkind)
        allocate(this%gpul(this%number_of_divisions + 1, this%number_of_conductors, this%number_of_conductors), source = 0.0_rkind)
        allocate(this%rpul(this%number_of_divisions, this%number_of_conductors, this%number_of_conductors), source = 0.0_rkind)
        allocate(this%du(this%number_of_divisions, this%number_of_conductors, this%number_of_conductors), source = 0.0_rkind)

        allocate(this%v(this%number_of_conductors, this%number_of_divisions + 1), source = 0.0_rkind)
        allocate(this%i(this%number_of_conductors, this%number_of_divisions), source = 0.0_rkind)
        allocate(this%i_prev(this%number_of_conductors, this%number_of_divisions), source = 0.0_rkind)
        allocate(this%e_L(this%number_of_conductors, this%number_of_divisions), source = 0.0_rkind)

        allocate(this%v_source(this%number_of_conductors, this%number_of_divisions + 1), source = 0.0_rkind)
        allocate(this%i_source(this%number_of_conductors, this%number_of_divisions), source = 0.0_rkind)

        allocate(this%i_term(this%number_of_divisions,this%number_of_conductors,this%number_of_conductors), source = 0.0_rkind)
        allocate(this%v_diff(this%number_of_divisions,this%number_of_conductors,this%number_of_conductors), source = 0.0_rkind)

        allocate(this%v_term(this%number_of_divisions + 1,this%number_of_conductors,this%number_of_conductors), source = 0.0_rkind)
        allocate(this%i_diff(this%number_of_divisions + 1,this%number_of_conductors,this%number_of_conductors), source = 0.0_rkind)

    end subroutine

    function countNumberOfConductors(levels) result(res)
        type(transmission_line_level_t), dimension(:), intent(in) :: levels
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
        type(transmission_line_level_t), dimension(:), intent(in) :: levels
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
        type(transmission_line_level_t), dimension(:), intent(in) :: levels
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

    function buildExternalFieldSegments(levels) result(res)
        type(transmission_line_level_t), dimension(:), intent(in) :: levels
        type(external_field_segment_t), dimension(:), allocatable :: res
        type(segment_t), dimension(:), allocatable :: segments
        integer :: i
        segments = levels(1)%lines(1)%segments
        allocate(res(size(segments)))
        do i = 1, size(segments)
            res(i)%position(1) = segments(i)%x
            res(i)%position(2) = segments(i)%y
            res(i)%position(3) = segments(i)%z
            res(i)%direction   = segments(i)%orientation
        end do
    end function

    subroutine addProbe(this, index, probe_type, name, position, layer_indices)
        class(mtl_bundle_t) :: this
        integer, intent(in) :: index
        integer, intent(in) :: probe_type
        real(kind=rkind), dimension(3) :: position
        character(len=:), allocatable :: name
        integer(kind=4), dimension(:,:), intent(in), optional :: layer_indices
        type(probe_t), allocatable, dimension(:) :: aux_probes
        type(probe_t) :: new_probe

        aux_probes = this%probes
        deallocate(this%probes)
        allocate(this%probes(size(aux_probes)+1))

#ifdef CompileWithMPI
        new_probe = probeCtor(index, probe_type, this%dt, name, position, layer_indices = layer_indices)
#else
        new_probe = probeCtor(index, probe_type, this%dt, name, position)
#endif
        this%probes(1:size(this%probes)-1) = aux_probes
        this%probes(size(aux_probes)+1) = new_probe
    end subroutine

    subroutine addGenerator(this, index, conductor, gen_type, resistance, path, layer_indices)
        class(mtl_bundle_t) :: this
        integer, intent(in) :: index, conductor, gen_type
        real(kind=rkind) :: resistance
        character(*), intent(in) :: path
        integer(kind=4), dimension(:,:), intent(in), optional :: layer_indices

        type(generator_t), allocatable, dimension(:) :: aux_generators
        type(generator_t) :: new_generator
        aux_generators = this%generators
#ifdef CompileWithMPI
        new_generator = generatorCtor(index, conductor, gen_type, resistance, path, layer_indices = layer_indices)
#else
        new_generator = generatorCtor(index, conductor, gen_type, resistance, path)
#endif
        ! Non-owning MPI ranks have no local PUL divisions to update.
        if (.not. new_generator%in_layer) return

        deallocate(this%generators)
        allocate(this%generators(size(aux_generators)+1))
        this%generators(1:size(this%generators)-1) = aux_generators
        this%generators(size(aux_generators)+1) = new_generator

        if (gen_type == SOURCE_TYPE_VOLTAGE) then
            this%rpul(new_generator%index, conductor, conductor) = &
                this%rpul(new_generator%index, conductor, conductor) + &
                resistance/this%du(new_generator%index, conductor, conductor)
        else if (new_generator%source_type == SOURCE_TYPE_CURRENT) then
            this%rpul(new_generator%index, conductor, conductor) = &
                this%rpul(new_generator%index, conductor, conductor) + &
                resistance/this%du(new_generator%index, conductor, conductor)

        end if

    end subroutine

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
        real(kind=rkind), dimension(this%number_of_divisions,this%number_of_conductors,this%number_of_conductors) :: F1, F2, IF1
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
        real(kind=rkind), dimension(this%number_of_divisions + 1,this%number_of_conductors,this%number_of_conductors) :: F1, F2, IF1
        real(kind=rkind), dimension(this%number_of_divisions + 1, this%number_of_conductors,this%number_of_conductors) :: extended_du
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


    end subroutine

    subroutine bundle_updateGenerators(this, time, dt)
        class(mtl_bundle_t) ::this
        real(kind=RKIND_TIEMPO), intent(in) :: time, dt
        real(kind=rkind) :: val
        integer :: i
        do i = 1, size(this%generators)
            if (this%generators(i)%in_layer) then 
                if (this%generators(i)%source_type == SOURCE_TYPE_VOLTAGE) then
                    val = 0.5*(this%generators(i)%interpolate(time+dt)+this%generators(i)%interpolate(time))
                    this%v_source(this%generators(i)%conductor, this%generators(i)%index) = & 
                        val/this%du(this%generators(i)%index, this%generators(i)%conductor, this%generators(i)%conductor)
                else if (this%generators(i)%source_type == SOURCE_TYPE_CURRENT) then
                    val = 0.5*(this%generators(i)%interpolate(time+dt)+this%generators(i)%interpolate(time))
                    this%v_source(this%generators(i)%conductor, this%generators(i)%index) = & 
                        val*this%generators(i)%resistance/this%du(this%generators(i)%index, this%generators(i)%conductor, this%generators(i)%conductor)
                end if
            end if
        end do

    end subroutine

#ifdef CompileWithMPI
    subroutine initMPIStateOwnership(this)
        class(mtl_bundle_t) :: this
        integer :: i, ierr, rank, n
        integer, allocatable :: candidate(:), owner(:)

        n = this%number_of_divisions
        allocate(candidate(n), owner(n), source = -1)
        call MPI_COMM_RANK(SUBCOMM_MPI, rank, ierr)
        do i = 1, n
            if (this%external_field_segments(i)%field_is_owned) candidate(i) = rank
        end do
        call MPI_Allreduce(candidate, owner, n, MPI_INTEGER, MPI_MAX, SUBCOMM_MPI, ierr)
        if (any(owner < 0)) error stop "MTLN MPI ownership missing for a current segment"
        if (allocated(this%segment_owner)) deallocate(this%segment_owner)
        call move_alloc(owner, this%segment_owner)
        do i = 1, n
            this%external_field_segments(i)%field_is_local = this%external_field_segments(i)%field_is_local .and. &
                this%segment_owner(i) == rank
        end do
        do i = 1, size(this%probes)
            if (this%probes(i)%type == PROBE_TYPE_CURRENT) then
                this%probes(i)%in_layer = this%segment_owner(this%probes(i)%index) == rank
            else
                this%probes(i)%in_layer = this%segment_owner(min(this%probes(i)%index, n)) == rank
            end if
        end do
        do i = 1, size(this%generators)
            this%generators(i)%in_layer = this%segment_owner(this%generators(i)%index) == rank
        end do
    end subroutine initMPIStateOwnership
#endif

    subroutine bundle_advanceVoltage(this)
        class(mtl_bundle_t) ::this
        integer :: i
#ifdef CompileWithMPI
        integer :: ierr, rank
        call MPI_COMM_RANK(SUBCOMM_MPI, rank, ierr)
#endif
        do i = 2,this%number_of_divisions
#ifdef CompileWithMPI
            if (allocated(this%segment_owner)) then
                if (this%segment_owner(i) /= rank) cycle
            end if
#endif
            this%v(:, i) = matmul(this%v_term(i,:,:), this%v(:,i)) - &
                           matmul(this%i_diff(i,:,:), (this%i(:,i) - this%i(:,i-1)) + matmul(this%du(i,:,:), this%i_source(:,i)))
        end do

    end subroutine


    subroutine bundle_advanceCurrent(this)
        class(mtl_bundle_t) ::this
        integer :: i
        real(kind=rkind) :: eps_r
#ifdef CompileWithMPI
        integer(kind=4) :: sizeof, ierr, rank
        call MPI_COMM_SIZE(SUBCOMM_MPI, sizeof, ierr)
        call MPI_COMM_RANK(SUBCOMM_MPI, rank, ierr)
        if (sizeof > 1) call this%Comm_MPI_V()
#endif
        call this%transfer_impedance%updateQ3Phi()
        this%i_prev = this%i

        do i = 1, this%number_of_divisions
#ifdef CompileWithMPI
            if (allocated(this%segment_owner)) then
                if (this%segment_owner(i) /= rank) cycle
            end if
#endif
            this%i(:,i) = matmul(this%i_term(i,:,:), this%i(:,i)) - &
                          matmul(this%v_diff(i,:,:), (this%v(:,i+1) - this%v(:,i)) - &
                                                      this%e_L(:,i) * this%step_size(i) - &
                                                      matmul(this%du(i,:,:),this%v_source(:,i))) - &
                          matmul(this%v_diff(i,:,:), matmul(this%du(i,:,:), this%transfer_impedance%q3_phi(i,:)))
        enddo
        call this%transfer_impedance%updatePhi(this%i_prev, this%i)
#ifdef CompileWithMPI
        if (sizeof > 1) call this%Comm_MPI_I()
#endif
    end subroutine

    subroutine bundle_setExternalLongitudinalField(this)
        class(mtl_bundle_t) :: this
        integer :: i, j
#ifdef CompileWithMPI
        integer :: ierr, rank
        call MPI_COMM_RANK(SUBCOMM_MPI, rank, ierr)
#endif
        this%e_L = 0.0_rkind
        do j = 1, this%conductors_in_level(1)
            do i = 1, size(this%e_L,2)
#ifdef CompileWithMPI
                if (allocated(this%segment_owner)) then
                    if (this%segment_owner(i) /= rank) cycle
                end if
#endif
                this%e_L(j,i) = this%external_field_segments(i)%field * &
                                  this%external_field_segments(i)%direction / abs(this%external_field_segments(i)%direction)
            end do
        end do
    end subroutine

#ifdef CompileWithMPI
    subroutine Comm_MPI_V(this)
        class(mtl_bundle_t) :: this
        integer :: i, ierr, rank, nreq, source_rank, destination_rank
        integer, allocatable :: requests(:), statuses(:,:)

        if (.not. allocated(this%segment_owner)) return
        call MPI_COMM_RANK(SUBCOMM_MPI, rank, ierr)
        allocate(requests(2 * max(0, this%number_of_divisions - 1)), &
                 statuses(MPI_STATUS_SIZE, 2 * max(0, this%number_of_divisions - 1)))
        nreq = 0
        do i = 2, this%number_of_divisions
            source_rank = this%segment_owner(i)
            destination_rank = this%segment_owner(i-1)
            if (source_rank == destination_rank) cycle
            if (rank == source_rank) then
                nreq = nreq + 1
                call MPI_Isend(this%v(:,i), size(this%v,1), REALSIZE, destination_rank, 501, SUBCOMM_MPI, requests(nreq), ierr)
            else if (rank == destination_rank) then
                nreq = nreq + 1
                call MPI_Irecv(this%v(:,i), size(this%v,1), REALSIZE, source_rank, 501, SUBCOMM_MPI, requests(nreq), ierr)
            end if
        end do
        if (nreq > 0) call MPI_Waitall(nreq, requests, statuses, ierr)
    end subroutine

    subroutine Comm_MPI_I(this)
        class(mtl_bundle_t) :: this
        integer :: i, ierr, rank, nreq, source_rank, destination_rank
        integer, allocatable :: requests(:), statuses(:,:)

        if (.not. allocated(this%segment_owner)) return
        call MPI_COMM_RANK(SUBCOMM_MPI, rank, ierr)
        allocate(requests(2 * max(0, this%number_of_divisions - 1)), &
                 statuses(MPI_STATUS_SIZE, 2 * max(0, this%number_of_divisions - 1)))
        nreq = 0
        do i = 2, this%number_of_divisions
            source_rank = this%segment_owner(i-1)
            destination_rank = this%segment_owner(i)
            if (source_rank == destination_rank) cycle
            if (rank == source_rank) then
                nreq = nreq + 1
                call MPI_Isend(this%i(:,i-1), size(this%i,1), REALSIZE, destination_rank, 502, SUBCOMM_MPI, requests(nreq), ierr)
            else if (rank == destination_rank) then
                nreq = nreq + 1
                call MPI_Irecv(this%i(:,i-1), size(this%i,1), REALSIZE, source_rank, 502, SUBCOMM_MPI, requests(nreq), ierr)
            end if
        end do
        if (nreq > 0) call MPI_Waitall(nreq, requests, statuses, ierr)
    end subroutine

        subroutine Comm_MPI_Fields(this)
        class(mtl_bundle_t) :: this
        integer :: i, ierr, rank, status(MPI_STATUS_SIZE)
        call MPI_COMM_RANK(SUBCOMM_MPI, rank, ierr)
        do i = 1, size(this%mpi_comm%comms)
            if (this%mpi_comm%comms(i)%comm_type == COMM_FIELD .or. this%mpi_comm%comms(i)%comm_type == COMM_BOTH) then
                if (this%mpi_comm%comms(i)%comm_task == COMM_SEND) then
                    call MPI_send(this%external_field_segments(this%mpi_comm%comms(i)%field_index)%field, 1, REALSIZE, &
                                    rank+this%mpi_comm%comms(i)%delta_rank, &
                                    100*(rank+this%mpi_comm%comms(i)%delta_rank+1), &
                                    SUBCOMM_MPI, ierr)
                end if
                if (this%mpi_comm%comms(i)%comm_task == COMM_RECV) then
                    call MPI_recv(this%external_field_segments(this%mpi_comm%comms(i)%field_index)%field,1, REALSIZE, &
                                    rank+this%mpi_comm%comms(i)%delta_rank, &
                                    100*(rank+1), &
                                    SUBCOMM_MPI, status, ierr)
                end if
            end if
        end do
    end subroutine

#endif
end module mtl_bundle_m
