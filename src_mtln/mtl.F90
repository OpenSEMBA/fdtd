module mtl_mod

    ! use NFDETypes
    use utils_mod
    use dispersive_mod
    use mtln_types_mod, only: external_field_segment_t, cable_t
#ifdef CompileWithMPI
        use FDETYPES, only: SUBCOMM_MPI, REALSIZE, INTEGERSIZE
#endif

    implicit none
#ifdef CompileWithMPI

    integer, parameter :: COMM_UNDEFINED = 0 
    integer, parameter :: COMM_NEXT = 1
    integer, parameter :: COMM_PREV = -1

    type, public :: side_communicator_t
        logical :: is_end = .true.
        integer :: index = -1
        integer :: comm_rank = -1
        integer :: direction = COMM_UNDEFINED
    end type

    type, public :: communicator_t
        type(side_communicator_t) :: left, right
    end type

    type, public :: comm_t
        type(communicator_t), dimension(:), allocatable :: comms
        integer :: rank
    end type
#endif

    type, public :: mtl_t
        character (len=:), allocatable :: name
        integer  :: number_of_conductors
        real, allocatable, dimension(:,:,:) :: lpul, cpul, rpul, gpul
        real, allocatable, dimension(:) :: step_size
        real, allocatable, dimension(:,:,:) :: du(:,:,:)
        type(lumped_t) :: lumped_elements
        real :: time = 0.0, dt = 0.0

        logical :: is_left_end = .true.
        logical :: is_right_end = .true.

        character (len=:), allocatable :: parent_name
        integer :: conductor_in_parent

        type(transfer_impedance_per_meter_t) :: transfer_impedance
        type(transfer_impedance_per_meter_t) :: initial_connector_transfer_impedance, end_connector_transfer_impedance


        type(external_field_segment_t), allocatable, dimension(:) :: external_field_segments
        logical :: isPassthrough = .false.

        integer (kind=4), allocatable, dimension(:,:) :: layer_indices
        ! integer (kind=4), dimension(1,2) :: layer_indices = [-1,-1]
        logical :: bundle_in_layer = .false.
#ifdef CompileWithMPI
        type(comm_t) :: mpi_comm
#endif

    contains
        procedure :: setTimeStep
        procedure :: initLCHomogeneous
        procedure :: initLCInhomogeneous
        procedure :: initRGHomogeneous
        procedure :: initRGInhomogeneous
        procedure :: initDirections
        procedure :: getMaxTimeStep
        procedure :: getPhaseVelocities
        !TODO
        ! procedure :: setResistanceInRegion
        ! procedure :: setResistanceAtPoint
        ! procedure :: setConductanceInRegion
        ! procedure :: setConductanceAtPoint
        ! procedure :: addDispersiveConnector
#ifdef CompileWithMPI
        procedure :: initCommunicators
#endif

    end type mtl_t

    interface mtl_t
        module procedure mtlHomogeneous
        module procedure mtlInhomogeneous
    end interface

    type, public :: mtl_array_t
        type(mtl_t), dimension(:), allocatable :: lines
    end type

    type, public :: line_bundle_t
        type(mtl_array_t), dimension(:), allocatable :: levels
    end type

contains


    subroutine initLCHomogeneous(this, lpul, cpul)
        class(mtl_t) :: this
        real, intent(in), dimension(:,:) :: lpul, cpul
        integer :: i
        allocate(this%lpul(size(this%step_size, 1), size(lpul, 1), size(lpul, 1)))
        allocate(this%cpul(size(this%step_size, 1) + 1, size(cpul,1), size(cpul, 1)))
        do i = 1, size(this%lpul, 1)
            this%lpul(i,:,:) = lpul(:,:)
        enddo
        do i = 1, size(this%cpul, 1)
            this%cpul(i,:,:) = cpul(:,:)
        enddo
    end subroutine

    subroutine initLCInhomogeneous(this, lpul, cpul)
        class(mtl_t) :: this
        real, intent(in), dimension(:, :, :) :: lpul, cpul
        ! integer :: i
        allocate(this%lpul(size(this%step_size, 1), size(lpul, 2), size(lpul, 2)))
        allocate(this%cpul(size(this%step_size, 1) + 1, size(cpul,2), size(cpul, 2)))
        this%lpul(:,:,:) = lpul(:,:,:)
        this%cpul(:,:,:) = cpul(:,:,:)
    end subroutine

    function mtlHomogeneous(lpul, cpul, rpul, gpul, &
                            step_size, name, &
                            dt, parent_name, conductor_in_parent, &
                            transfer_impedance, &
                            external_field_segments, &
                            isPassthrough, layer_indices, bundle_in_layer) result(res)
        type(mtl_t) :: res
        real, intent(in), dimension(:,:) :: lpul, cpul, rpul, gpul
        real, intent(in), dimension(:) :: step_size
        character(len=*), intent(in) :: name

        real, intent(in), optional :: dt
        real :: max_dt
        character(len=*), intent(in), optional :: parent_name
        integer, intent(in), optional :: conductor_in_parent
        type(transfer_impedance_per_meter_t), intent(in), optional :: transfer_impedance
        type(external_field_segment_t), intent(in), dimension(:), optional :: external_field_segments
        logical, optional :: isPassthrough
        integer (kind=4), allocatable, dimension(:,:), intent(in), optional :: layer_indices
        logical, optional :: bundle_in_layer

        res%name = name
        if (present(layer_indices)) then
#ifdef CompileWithMPI
            call res%initCommunicators(step_size, external_field_segments, layer_indices)
#else
            res%step_size =  step_size
            res%external_field_segments = external_field_segments
#endif
          res%layer_indices = layer_indices
        end if


        if (present(bundle_in_layer)) res%bundle_in_layer = bundle_in_layer
        call checkPULDimensionsHomogeneous(lpul, cpul, rpul, gpul)
        res%number_of_conductors = size(lpul, 1)

        call res%initDirections()
        call res%initLCHomogeneous(lpul, cpul)
        call res%initRGHomogeneous(rpul, gpul)
        

        if (present(dt)) then 
            if (lpul(1,1) /= 0.0) then 
                max_dt = res%getMaxTimeStep() 
                if (dt > max_dt) then
                    res%dt = max_dt
                    write(*,*) 'dt larger than maximum permitted. Changed to dt = ', max_dt 
                else 
                    res%dt = dt
                end if
            else 
                res%dt = dt
            end if
        else
            if (lpul(1,1) /= 0.0) then 
                res%dt = res%getMaxTimeStep() 
            end if
        end if

     
        res%lumped_elements = lumped_t(res%number_of_conductors, 0, size(res%step_size), res%dt)
        if (present(parent_name)) res%parent_name = parent_name
        if (present(conductor_in_parent)) res%conductor_in_parent = conductor_in_parent
        if (present(transfer_impedance)) res%transfer_impedance = transfer_impedance
        if (present(isPassthrough)) res%isPassthrough = isPassthrough

    end function

    function mtlInhomogeneous(lpul, cpul, rpul, gpul, &
                            step_size, name, &
                            dt, parent_name, conductor_in_parent, &
                            transfer_impedance, &
                            external_field_segments, &
                            isPassthrough) result(res)
        type(mtl_t) :: res
        real, intent(in), dimension(:,:,:) :: lpul, cpul, rpul, gpul
        real, intent(in), dimension(:) :: step_size
        character(len=*), intent(in) :: name
        real, intent(in), optional :: dt
        real :: max_dt
        character(len=*), intent(in), optional :: parent_name
        integer, intent(in), optional :: conductor_in_parent
        type(transfer_impedance_per_meter_t), intent(in), optional :: transfer_impedance
        type(external_field_segment_t), intent(in), dimension(:), optional :: external_field_segments
        logical, optional :: isPassthrough

        res%name = name
        res%step_size = step_size
        call checkPULDimensionsInhomogeneous(lpul, cpul, rpul, gpul)
        res%number_of_conductors = size(lpul, 2)

        call res%initDirections()
        call res%initLCInhomogeneous(lpul, cpul)
        call res%initRGInhomogeneous(rpul, gpul)

        if (present(dt)) then 
            if (lpul(1,1,1) /= 0.0) then 
                max_dt = res%getMaxTimeStep() 
                if (dt > max_dt) then
                    res%dt = max_dt
                    write(*,*) 'dt larger than maximum permitted. Changed to dt = ', max_dt 
                else 
                    res%dt = dt
                end if
            else 
                res%dt = dt
            end if
        else
            if (lpul(1,1,1) /= 0.0) then 
                res%dt = res%getMaxTimeStep() 
            end if
        end if

        res%lumped_elements = lumped_t(res%number_of_conductors, 0, size(step_size), res%dt)
        if (present(parent_name)) res%parent_name = parent_name
        if (present(conductor_in_parent)) res%conductor_in_parent = conductor_in_parent
        if (present(transfer_impedance))  res%transfer_impedance = transfer_impedance
        if (present(external_field_segments)) res%external_field_segments = external_field_segments
        if (present(isPassthrough)) res%isPassthrough = isPassthrough

    end function

    subroutine initDirections(this)
        class(mtl_t) :: this
        integer :: j

        this%du = reshape(source = [(this%step_size(j)*eye(this%number_of_conductors) , j = 1, size(this%step_size))], &
                              shape = [size(this%step_size), this%number_of_conductors, this%number_of_conductors], &
                              order=[2,3,1])
    end subroutine

    subroutine checkPULDimensionsHomogeneous(lpul, cpul, rpul, gpul)
        real, intent(in), dimension(:,:) :: lpul, cpul, rpul, gpul

        if ((size(lpul, 1) /= size(lpul, dim = 2)).or.&
            (size(cpul, 1) /= size(cpul, dim = 2)).or.&
            (size(rpul, 1) /= size(rpul, dim = 2)).or.&
            (size(gpul, 1) /= size(gpul, dim = 2))) then
            error stop 'PUL matrices are not square'
        endif

        if ((size(lpul, 1) /= size(cpul, 1)).or.&
            (size(lpul, 1) /= size(rpul, 1)).or.&
            (size(lpul, 1) /= size(gpul, 1))) then
            error stop 'PUL matrices do not have the same dimensions'
        endif

    end subroutine

    subroutine checkPULDimensionsInhomogeneous(lpul, cpul, rpul, gpul)
        real, intent(in), dimension(:, :,:) :: lpul, cpul, rpul, gpul

        if ((size(lpul, 2) /= size(lpul, dim = 3)).or.&
            (size(cpul, 2) /= size(cpul, dim = 3)).or.&
            (size(rpul, 2) /= size(rpul, dim = 3)).or.&
            (size(gpul, 2) /= size(gpul, dim = 3))) then
            error stop 'PUL matrices are not square'
        endif

        if ((size(lpul, 2) /= size(cpul, 2)).or.&
            (size(lpul, 2) /= size(rpul, 2)).or.&
            (size(lpul, 2) /= size(gpul, 2))) then
            error stop 'PUL matrices do not have the same dimensions'
        endif

    end subroutine

    function getPhaseVelocities(this) result(res)
        class(mtl_t) :: this
        real, dimension(size(this%step_size,1), this%number_of_conductors) :: res
        real, dimension(2*this%number_of_conductors) :: ev
        integer :: k

        do k = 1, size(this%step_size, 1)
            ev = getEigenValues(dble(matmul(this%lpul(k,:,:), this%cpul(k+1,:,:))))
            res(k,:) = 1.0/sqrt(ev(1:this%number_of_conductors))
        enddo

    end function

    function getMaxTimeStep(this) result(res)
        class(mtl_t) :: this
        real :: res
        res= minval(pack(this%du, this%du /= 0))/maxval(this%getPhaseVelocities())
    end function

    subroutine initRGHomogeneous(this, rpul, gpul)
        class(mtl_t) :: this
        real, intent(in), dimension(:,:) :: rpul, gpul
        integer :: i
        allocate(this%rpul(size(this%step_size, 1), size(rpul, 1), size(rpul, 1)))
        allocate(this%gpul(size(this%step_size, 1) + 1, size(gpul, 1), size(gpul, 1)))

        do i = 1, size(this%rpul, 1)
            this%rpul(i,:,:) = rpul(:,:)
        enddo
        do i = 1, size(this%gpul, 1)
            this%gpul(i,:,:) = gpul(:,:)
        enddo
    end subroutine

    subroutine initRGInhomogeneous(this, rpul, gpul)
        class(mtl_t) :: this
        real, intent(in), dimension(:, :,:) :: rpul, gpul
        allocate(this%rpul(size(this%step_size, 1), size(rpul, 2), size(rpul, 2)))
        allocate(this%gpul(size(this%step_size, 1) + 1, size(gpul, 2), size(gpul, 2)))
        this%rpul(:,:,:) = rpul(:,:,:)
        this%gpul(:,:,:) = gpul(:,:,:)
    end subroutine


    subroutine setTimeStep(this, numberOfSteps, finalTime)
        class(mtl_t) :: this
        integer, intent(in) :: numberOfSteps
        real, intent (in) ::finalTime
        this%dt = finalTime/numberOfSteps
    end subroutine


#ifdef CompileWithMPI
    subroutine initCommunicators(this, step_size, external_field_segments, layer_indices)
        class(mtl_t) :: this
        real, intent(in), dimension(:) :: step_size
        type(external_field_segment_t), intent(in), dimension(:) :: external_field_segments
        integer (kind=4), allocatable, dimension(:,:), intent(in) :: layer_indices
        integer :: j, n
        integer :: rank, ierr, inter, sizeof

            call MPI_COMM_RANK(SUBCOMM_MPI, rank, ierr)
            call MPI_COMM_SIZE(SUBCOMM_MPI, sizeof, ierr)

            write(*,*) 'rank: ',rank
            this%mpi_comm%rank = rank
            allocate(this%mpi_comm%comms(size(layer_indices,1)))
            write(*,*) layer_indices
            n =  0
            do j = 1, size(layer_indices, 1)
                n = n + layer_indices(j,2) - layer_indices(j,1) + 1
            end do
            n = n + size(layer_indices, 1) -1
            allocate(this%step_size(n))
            allocate(this%external_field_segments(n))

            n = 1
            inter = 0
            do j = 1, size(layer_indices, 1)
                this%step_size(n: n + layer_indices(j,2) - layer_indices(j,1)) = step_size(layer_indices(j, 1):layer_indices(j, 2))
                this%external_field_segments(n: n + layer_indices(j,2) - layer_indices(j,1)) = external_field_segments(layer_indices(j, 1):layer_indices(j, 2))
                
                ! if (layer_indices(j, 1) + inter /= 1) then 
                if (layer_indices(j, 1)  /= 1) then 
                    this%mpi_comm%comms(j)%left%is_end = .false.
                end if
                this%mpi_comm%comms(j)%left%index = n + 1
                ! if (layer_indices(j, 2) + inter /= size(step_size)) then 
                
                if (layer_indices(j, 2)  /= size(step_size)) then 
                    this%mpi_comm%comms(j)%right%is_end = .false.
                end if
                this%mpi_comm%comms(j)%right%index = n + layer_indices(j,2) - layer_indices(j,1)

                if (.not. this%mpi_comm%comms(j)%left%is_end) then 
                    if (external_field_segments(layer_indices(j, 1))%direction > 0) then 
                        this%mpi_comm%comms(j)%left%comm_rank = this%mpi_comm%rank - 1
                        this%mpi_comm%comms(j)%left%direction = COMM_PREV
                    else 
                        this%mpi_comm%comms(j)%left%comm_rank = this%mpi_comm%rank + 1
                        this%mpi_comm%comms(j)%left%direction = COMM_NEXT
                    end if
                else 
                    if (rank == sizeof -1) then 
                        this%mpi_comm%comms(j)%left%comm_rank = this%mpi_comm%rank - 1
                        this%mpi_comm%comms(j)%left%direction = COMM_PREV
                    else if (rank == 0) then
                        this%mpi_comm%comms(j)%left%comm_rank = this%mpi_comm%rank + 1
                        this%mpi_comm%comms(j)%left%direction = COMM_NEXT
                    else
                        if (external_field_segments(layer_indices(j, 2))%direction > 0) then 
                            this%mpi_comm%comms(j)%left%comm_rank = this%mpi_comm%rank - 1
                            this%mpi_comm%comms(j)%left%direction = COMM_PREV
                        else
                            this%mpi_comm%comms(j)%left%comm_rank = this%mpi_comm%rank + 1
                            this%mpi_comm%comms(j)%left%direction = COMM_NEXT
                        end if
                    end if
                end if

                
                if (.not. this%mpi_comm%comms(j)%right%is_end) then 
                    if (external_field_segments(layer_indices(j, 2))%direction > 0) then 
                        this%mpi_comm%comms(j)%right%comm_rank = this%mpi_comm%rank + 1
                        this%mpi_comm%comms(j)%right%direction = COMM_NEXT
                    else
                        this%mpi_comm%comms(j)%right%comm_rank = this%mpi_comm%rank - 1
                        this%mpi_comm%comms(j)%right%direction = COMM_PREV
                    end if
                else 
                    if (rank == sizeof -1) then 
                        this%mpi_comm%comms(j)%right%comm_rank = this%mpi_comm%rank - 1
                        this%mpi_comm%comms(j)%right%direction = COMM_PREV
                    else if (rank == 0) then
                        this%mpi_comm%comms(j)%right%comm_rank = this%mpi_comm%rank + 1
                        this%mpi_comm%comms(j)%right%direction = COMM_NEXT
                    else
                        if (external_field_segments(layer_indices(j, 2))%direction > 0) then 
                            this%mpi_comm%comms(j)%right%comm_rank = this%mpi_comm%rank + 1
                            this%mpi_comm%comms(j)%right%direction = COMM_NEXT
                        else
                            this%mpi_comm%comms(j)%right%comm_rank = this%mpi_comm%rank - 1
                            this%mpi_comm%comms(j)%right%direction = COMM_PREV
                        end if
                    end if
                end if

                if (j /= size(layer_indices,1)) then 
                    this%step_size(n + layer_indices(j,2) - layer_indices(j,1) + 1) = this%step_size(n + layer_indices(j,2) - layer_indices(j,1))
                    this%external_field_segments(n + layer_indices(j,2) - layer_indices(j,1) + 1) = this%external_field_segments(n + layer_indices(j,2) - layer_indices(j,1)) 
                    n = n + 1
                    inter = inter + 1
                end if

                n = n + layer_indices(j,2) - layer_indices(j,1) + 1
            end do
    end subroutine
#endif

end module
