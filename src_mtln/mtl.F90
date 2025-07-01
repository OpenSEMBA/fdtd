module mtl_mod

    ! use NFDETypes
    use utils_mod
    use dispersive_mod, dispersive_lumped_t => lumped_t
    use mtln_types_mod, only: segment_t, multipolar_expansion_t
    use multipolar_expansion_mod, only: getCellCapacitanceOnBox, getCellInductanceOnBox
#ifdef CompileWithMPI
    use fdetypes, only: SUBCOMM_MPI, REALSIZE, INTEGERSIZE
#endif
    implicit none
#ifdef CompileWithMPI

    integer, parameter :: COMM_SEND = 1
    integer, parameter :: COMM_RECV = -1
    integer, parameter :: COMM_NONE = 0

    type, public :: communicator_t
        integer :: field_index = -1, v_index = -1
        integer :: comm_task = COMM_NONE
        integer :: delta_rank = 0
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
        type(dispersive_lumped_t) :: lumped_elements
        real :: time = 0.0, dt = 0.0

        character (len=:), allocatable :: parent_name
        integer :: conductor_in_parent
        type(transfer_impedance_per_meter_t) :: transfer_impedance
        type(transfer_impedance_per_meter_t) :: initial_connector_transfer_impedance, end_connector_transfer_impedance
        type(segment_t), dimension(:), allocatable :: segments

#ifdef CompileWithMPI
        type(comm_t) :: mpi_comm
        integer (kind=4), allocatable, dimension(:,:) :: layer_indices
        logical :: bundle_in_layer = .true.
#endif

    contains
        procedure :: setTimeStep
        procedure :: checkTimeStep
        procedure :: allocatePULMatrices
        procedure :: computeLCParameters
        procedure :: initLC
        procedure :: initRG
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
        procedure :: initStepSizeAndFieldSegments
#endif

    end type mtl_t

    interface mtl_t
        module procedure mtl_shielded
        module procedure mtl_unshielded
    end interface

    type, public :: transmission_line_level_t
        type(mtl_t), dimension(:), allocatable :: lines
    end type

    type, public :: transmission_line_bundle_t
        type(transmission_line_level_t), dimension(:), allocatable :: levels
    end type

contains


    subroutine initLC(this, lpul, cpul)
        class(mtl_t) :: this
        real, intent(in), dimension(:,:) :: lpul, cpul
        integer :: i
        do i = 1, size(this%lpul, 1)
            this%lpul(i,:,:) = lpul(:,:)
        enddo
        do i = 1, size(this%cpul, 1)
            this%cpul(i,:,:) = cpul(:,:)
        enddo
    end subroutine

    function mtl_shielded(lpul, cpul, rpul, gpul, &
                            step_size, name, segments, &
                            dt, parent_name, conductor_in_parent, &
                            transfer_impedance, &
                            layer_indices, bundle_in_layer, alloc_z) result(res)
        real, intent(in), dimension(:,:) :: lpul, cpul, rpul, gpul
        real, intent(in), dimension(:) :: step_size
        character(len=*), intent(in) :: name
        type(segment_t), dimension(:), allocatable, intent(in) :: segments
        real, intent(in) :: dt
        character(len=*), intent(in) :: parent_name
        integer, intent(in) :: conductor_in_parent
        type(transfer_impedance_per_meter_t), intent(in) :: transfer_impedance

        integer (kind=4), allocatable, dimension(:,:), intent(in), optional :: layer_indices
        logical, optional :: bundle_in_layer
        integer(kind=4), dimension (2), intent(in), optional :: alloc_z
        
        type(mtl_t) :: res

        real :: max_dt

#ifdef CompileWithMPI
        integer (kind=4) :: sizeof, ierr
        if (present(layer_indices)) then 
            ! call res%initStepSizeAndFieldSegments(step_size, external_field_segments, layer_indices)
            ! => call res%initStepSizeAndFieldSegments(step_size, segments, layer_indices)
            call res%initCommunicators(alloc_z)
            res%layer_indices = layer_indices
            res%bundle_in_layer = bundle_in_layer
        else
            res%step_size =  step_size
            allocate(res%layer_indices(0,0))
            allocate(res%mpi_comm%comms(0))
            ! res%external_field_segments = external_field_segments
        end if
#else
        res%step_size =  step_size
        res%segments = segments
        ! res%external_field_segments = external_field_segments
#endif

        call checkPULDimensions(lpul, cpul, rpul, gpul)
        res%name = name
        res%number_of_conductors = size(lpul, 1)
        call res%initDirections()
        call res%allocatePULMatrices()
        call res%initLC(lpul, cpul)
        call res%initRG(rpul, gpul)
        call res%checkTimeStep(dt)
    
        res%parent_name = parent_name
        res%conductor_in_parent = conductor_in_parent
        res%transfer_impedance = transfer_impedance
        res%lumped_elements = dispersive_lumped_t(res%number_of_conductors, 0, size(res%step_size), res%dt)

    end function

    function mtl_unshielded(lpul, cpul, rpul, gpul, &
                            step_size, name, segments, &
                            dt, multipolar_expansion, &
                            layer_indices, bundle_in_layer, alloc_z) result(res)
        real, intent(in), dimension(:,:) :: lpul, cpul, rpul, gpul
        real, intent(in), dimension(:) :: step_size
        character(len=*), intent(in) :: name
        type(segment_t), dimension(:), allocatable, intent(in) :: segments
        real, intent(in) :: dt
        type(multipolar_expansion_t), dimension(:), allocatable :: multipolar_expansion

        integer (kind=4), allocatable, dimension(:,:), intent(in), optional :: layer_indices
        logical, optional :: bundle_in_layer
        integer(kind=4), dimension (2), intent(in), optional :: alloc_z
        
        type(mtl_t) :: res
        real :: max_dt
#ifdef CompileWithMPI
        integer (kind=4) :: sizeof, ierr
        if (present(layer_indices)) then 
            ! call res%initStepSizeAndFieldSegments(step_size, external_field_segments, layer_indices)
            ! => call res%initStepSizeAndFieldSegments(step_size, segments, layer_indices)
            call res%initCommunicators(alloc_z)
            res%layer_indices = layer_indices
            res%bundle_in_layer = bundle_in_layer
        else
            res%step_size =  step_size
            allocate(res%layer_indices(0,0))
            allocate(res%mpi_comm%comms(0))
            ! res%external_field_segments = external_field_segments
        end if
#else
        res%step_size =  step_size
        res%segments = segments
        ! res%external_field_segments = external_field_segments
#endif

        call checkPULDimensions(lpul, cpul, rpul, gpul)
        res%name = name
        res%number_of_conductors = size(lpul, 1)
        call res%initDirections()
        call res%allocatePULMatrices()
        if (size(multipolar_expansion) /= 0) then 
            call res%computeLCParameters(multipolar_expansion(1))
        else
            call res%initLC(lpul, cpul)
        end if
        call res%initRG(rpul, gpul)
        call res%checkTimeStep(dt)
        res%lumped_elements = dispersive_lumped_t(res%number_of_conductors, 0, size(res%step_size), res%dt)
    end function

    subroutine checkTimeStep(this, dt)
        class(mtl_t) :: this
        real, intent(in), optional :: dt
        real :: max_dt
        if (present(dt)) then 
            if (this%lpul(1,1,1) /= 0.0) then 
                max_dt = this%getMaxTimeStep() 
                if (dt > max_dt) then
                    this%dt = max_dt
                    write(*,*) 'dt larger than maximum permitted. Changed to dt = ', max_dt 
                else 
                    this%dt = dt
                end if
            else 
                this%dt = dt
            end if
        else
            if (this%lpul(1,1,1) /= 0.0) then 
                this%dt = this%getMaxTimeStep() 
            end if
        end if
    end subroutine

    subroutine allocatePULMatrices(this)
        class(mtl_t) :: this
        integer :: n
        n = this%number_of_conductors
        allocate(this%lpul(size(this%step_size, 1),     n, n))
        allocate(this%cpul(size(this%step_size, 1) + 1, n, n))
        allocate(this%rpul(size(this%step_size, 1),     n, n))
        allocate(this%gpul(size(this%step_size, 1) + 1, n, n))
    end subroutine

    subroutine computeLCParameters(this, multipolar_expansion)
        class(mtl_t) :: this
        type(multipolar_expansion_t), intent(in) :: multipolar_expansion
        integer :: i
        do i = 1, size(this%segments)
            this%lpul(i,:,:) = getCellInductanceOnBox(multipolar_expansion, this%segments(i)%dualBox)
            this%cpul(i,:,:) = getCellCapacitanceOnBox(multipolar_expansion, this%segments(i)%dualBox)
        end do
        this%cpul(size(this%segments)+1, :,:) = this%cpul(size(this%segments), :,:)
    end subroutine


    subroutine initDirections(this)
        class(mtl_t) :: this
        integer :: j

        this%du = reshape(source = [(this%step_size(j)*eye(this%number_of_conductors) , j = 1, size(this%step_size))], &
                              shape = [size(this%step_size), this%number_of_conductors, this%number_of_conductors], &
                              order=[2,3,1])
    end subroutine

    subroutine checkPULDimensions(lpul, cpul, rpul, gpul)
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

    subroutine initRG(this, rpul, gpul)
        class(mtl_t) :: this
        real, intent(in), dimension(:,:) :: rpul, gpul
        integer :: i
        do i = 1, size(this%rpul, 1)
            this%rpul(i,:,:) = rpul(:,:)
        enddo
        do i = 1, size(this%gpul, 1)
            this%gpul(i,:,:) = gpul(:,:)
        enddo
    end subroutine


    subroutine setTimeStep(this, numberOfSteps, finalTime)
        class(mtl_t) :: this
        integer, intent(in) :: numberOfSteps
        real, intent (in) ::finalTime
        this%dt = finalTime/numberOfSteps
    end subroutine


#ifdef CompileWithMPI

    subroutine initStepSizeAndFieldSegments(this, step_size, layer_indices)
    ! subroutine initStepSizeAndFieldSegments(this, step_size, external_field_segments, layer_indices)
        class(mtl_t) :: this
        real, intent(in), dimension(:) :: step_size
        ! type(external_field_segment_t), intent(in), dimension(:) :: external_field_segments
        integer (kind=4), allocatable, dimension(:,:), intent(in) :: layer_indices
        integer :: n, j
        n =  0
        do j = 1, size(layer_indices, 1)
            n = n + layer_indices(j,2) - layer_indices(j,1) + 1
        end do
        n = n + size(layer_indices, 1) -1
        allocate(this%step_size(n))
        ! allocate(this%external_field_segments(n))

        n = 1
        do j = 1, size(layer_indices, 1)
            this%step_size(n: n + layer_indices(j,2) - layer_indices(j,1)) = step_size(layer_indices(j, 1):layer_indices(j, 2))
            ! this%external_field_segments(n: n + layer_indices(j,2) - layer_indices(j,1)) = external_field_segments(layer_indices(j, 1):layer_indices(j, 2))

            if (j /= size(layer_indices,1)) then
                this%step_size(n + layer_indices(j,2) - layer_indices(j,1) + 1) = this%step_size(n + layer_indices(j,2) - layer_indices(j,1))
                ! this%external_field_segments(n + layer_indices(j,2) - layer_indices(j,1) + 1) = this%external_field_segments(n + layer_indices(j,2) - layer_indices(j,1))
                n = n + 1
            end if
            n = n + layer_indices(j,2) - layer_indices(j,1) + 1
        end do

    end subroutine


    subroutine initCommunicators(this, alloc_z)
        class(mtl_t) :: this
        integer (kind =4), dimension(2) :: alloc_z
        integer :: j, n
        integer :: rank, ierr
        integer (kind =4) :: z_init, z_end
        type(communicator_t), dimension(:), allocatable :: aux_comm

        call MPI_COMM_RANK(SUBCOMM_MPI, rank, ierr)
        this%mpi_comm%rank = rank
        allocate(this%mpi_comm%comms(0))
        allocate(aux_comm(0))
        z_init = alloc_z(1)
        z_end = alloc_z(2)

        ! do j = 1, size(this%external_field_segments)
            
        !     if (isSegmentZOriented(j) .and. &
        !        (isSegmentNextToLayerEnd(j,z_end) .or. isSegmentNextToLayerInit(j,z_init))) then
                
        !         n = size(this%mpi_comm%comms)
        !         deallocate(aux_comm)
        !         allocate(aux_comm(n+1))
        !         aux_comm(1:n) = this%mpi_comm%comms
        !         aux_comm(n+1)%field_index = j

        !         if (isSegmentNextToLayerEnd(j,z_end)) then 
        !             aux_comm(n+1)%delta_rank = 1

        !             if (isSegmentBeforeLayerEnd(j,z_end)) then 
        !                 aux_comm(n+1)%comm_task = COMM_SEND
        !                 if (isSegmentZPositive(j)) then 
        !                     aux_comm(n+1)%v_index = j
        !                 else 
        !                     aux_comm(n+1)%v_index = j+1
        !                 end if
        !             else if (isSegmentAfterLayerEnd(j,z_end)) then 
        !                 aux_comm(n+1)%comm_task = COMM_RECV
        !                 if (isSegmentZPositive(j)) then 
        !                     aux_comm(n+1)%v_index = j+1
        !                 else 
        !                     aux_comm(n+1)%v_index = j
        !                 end if
        !             end if

        !         else if  (isSegmentNextToLayerInit(j,z_init)) then 
        !             aux_comm(n+1)%delta_rank = -1

        !             if (isSegmentBeforeLayerInit(j,z_init)) then 
        !                 aux_comm(n+1)%comm_task = COMM_RECV
        !                 if (isSegmentZPositive(j)) then 
        !                     aux_comm(n+1)%v_index = j
        !                 else
        !                     aux_comm(n+1)%v_index = j+1
        !                 end if
        !             else if (isSegmentAfterLayerInit(j,z_init)) then 
        !                 aux_comm(n+1)%comm_task = COMM_SEND
        !                 if (isSegmentZPositive(j)) then 
        !                     aux_comm(n+1)%v_index = j+1
        !                 else
        !                     aux_comm(n+1)%v_index = j
        !                 end if
        !             end if

        !         end  if
        !         deallocate(this%mpi_comm%comms)
        !         allocate(this%mpi_comm%comms(n+1))
        !         this%mpi_comm%comms = aux_comm
        !     end if
        ! end do

    ! contains    

    ! logical function isSegmentZOriented(j)
    !     integer, intent(in) :: j
    !     isSegmentZOriented = (abs(this%external_field_segments(j)%direction) == 3)
    ! end function

    ! logical function isSegmentZPositive(j)
    !     integer, intent(in) :: j
    !     isSegmentZPositive = (this%external_field_segments(j)%direction > 0)
    ! end function

    ! logical function isSegmentBeforeLayerEnd(j, z_end)
    !     integer, intent(in) :: j, z_end
    !     integer :: z
    !     z = this%external_field_segments(j)%position(3)
    !     isSegmentBeforeLayerEnd = (z==z_end-1)
    ! end function

    ! logical function isSegmentAfterLayerEnd(j, z_end)
    !     integer, intent(in) :: j, z_end
    !     integer :: z
    !     z = this%external_field_segments(j)%position(3)
    !     isSegmentAfterLayerEnd = (z==z_end)
    ! end function
    
    ! logical function isSegmentBeforeLayerInit(j, z_init)
    !     integer, intent(in) :: j, z_init
    !     integer :: z
    !     z = this%external_field_segments(j)%position(3)
    !     isSegmentBeforeLayerInit = (z==z_init)
    ! end function
    
    ! logical function isSegmentAfterLayerInit(j, z_init)
    !     integer, intent(in) :: j, z_init
    !     integer :: z
    !     z = this%external_field_segments(j)%position(3)
    !     isSegmentAfterLayerInit = (z==z_init+1)
    ! end function

    ! logical function isSegmentNextToLayerEnd(j, z_end)
    !     integer, intent(in) :: j, z_end
    !     integer :: z
    !     z = this%external_field_segments(j)%position(3)
    !     isSegmentNextToLayerEnd = (abs(z-z_end)<= 1)
    ! end function

    ! logical function isSegmentNextToLayerInit(j, z_init)
    !     integer, intent(in) :: j, z_init
    !     integer :: z
    !     z = this%external_field_segments(j)%position(3)
    !     isSegmentNextToLayerInit = (abs(z-z_init-1) <= 1)
    ! end function



    end subroutine
#endif

end module
