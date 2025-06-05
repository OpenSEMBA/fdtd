module mtln_solver_mod 

    use mtl_bundle_mod
    use network_manager_mod
    use preprocess_mod
    use FDETYPES, only: XYZlimit_t
#ifdef CompileWithMPI
    use FDETYPES, only: SUBCOMM_MPI, REALSIZE, INTEGERSIZE, MPI_STATUS_SIZE
#endif
    implicit none


    type, public :: mtln_t
        real :: time, dt, final_time
        type(mtl_bundle_t), allocatable, dimension(:) :: bundles
        type(network_manager_t) :: network_manager
        type(probe_t), allocatable, dimension(:) :: probes
        integer :: number_of_bundles
        logical :: has_multiwires
        integer :: step_i
    contains

        procedure :: updateBundlesTimeStep
        procedure :: updatePULTerms
        procedure :: initNodes
        procedure :: getTimeRange
        procedure :: updateProbes
        procedure :: advanceNWVoltage
        procedure :: advanceBundlesVoltage
        procedure :: advanceBundlesCurrent
        procedure :: advanceTime
        procedure :: step => mtln_step
        procedure :: step_alone
        procedure :: setExternalLongitudinalField
        
        procedure :: runUntil
        procedure :: run => mtln_run

    end type mtln_t


    interface mtln_t
        module procedure mtlnCtor
    end interface
    

contains

    function mtlnCtor(parsed, alloc) result(res)
        type(parsed_mtln_t) :: parsed
        type (XYZlimit_t), dimension (1:6), intent(in), optional :: alloc
        type(mtln_t) :: res
        integer :: i
        type(preprocess_t) :: pre

#ifdef CompileWithMPI
        integer (kind=4) :: sizeof, ierr
#endif

#ifdef CompileWithMPI
        call mpi_barrier(subcomm_mpi, ierr)
        
#endif
        if (present(alloc)) then 
            pre = preprocess(parsed, alloc)
        else  
            pre = preprocess(parsed)
        endif

        if (size(pre%bundles) == 0) then
            res%number_of_bundles = 0
            res%network_manager%has_networks = .false.
            return
        end if

        res%has_multiwires = parsed%has_multiwires

        res%dt = pre%dt
        res%time  = 0.0
        res%final_time = pre%final_time

        res%bundles = pre%bundles
        res%number_of_bundles = size(res%bundles)
        
        res%network_manager = pre%network_manager
        res%probes = pre%probes
        call res%updateBundlesTimeStep(res%dt)
        call res%initNodes()

        res%step_i = 0

    end function

    subroutine initNodes(this)
        class(mtln_t) :: this
        integer :: i,j
        do i = 1, size(this%network_manager%networks)
            do j = 1, size(this%network_manager%networks(i)%nodes)
                this%network_manager%networks(i)%nodes(j)%v = 0.0
                this%network_manager%networks(i)%nodes(j)%i = 0.0
            end do
        end do
    end subroutine

    subroutine mtln_step(this)
        class(mtln_t) :: this
        call this%setExternalLongitudinalField()
        call this%advanceBundlesVoltage()
        call this%advanceNWVoltage()
        call this%advanceBundlesCurrent()
        call this%advanceTime()
        call this%updateProbes()

    end subroutine

    subroutine step_alone(this)
        class(mtln_t) :: this
        integer :: i 

        call this%advanceBundlesVoltage()
        call this%advanceNWVoltage()
        call this%advanceBundlesCurrent()

        call this%advanceTime()
        call this%updateProbes()

    end subroutine

    subroutine setExternalLongitudinalField(this)
        class(mtln_t) :: this
        integer :: i
        do i = 1, this%number_of_bundles
            if (this%bundles(i)%bundle_in_layer) call this%bundles(i)%setExternalLongitudinalField()
        end do

    end subroutine

    subroutine advanceBundlesVoltage(this)
        class(mtln_t) :: this
        integer :: i

        do i = 1, this%number_of_bundles
            if (this%bundles(i)%bundle_in_layer) then 
                call this%bundles(i)%updateSources(this%time, this%dt)
                call this%bundles(i)%advanceVoltage()
            end if
        end do

    end subroutine

    subroutine advanceNWVoltage(this)
        class(mtln_t) :: this
        integer :: i,j
        integer ::b, c, v_idx, i_idx
        integer :: n

        if (this%number_of_bundles /= 0) then 
            do i = 1, size(this%network_manager%networks)
                do j = 1, size(this%network_manager%networks(i)%nodes)
                    b = this%network_manager%networks(i)%nodes(j)%bundle_number
                    c = this%network_manager%networks(i)%nodes(j)%conductor_number
                    v_idx = this%network_manager%networks(i)%nodes(j)%v_index
                    i_idx = this%network_manager%networks(i)%nodes(j)%i_index
                    if (this%bundles(b)%bundle_in_layer) this%network_manager%networks(i)%nodes(j)%i = this%bundles(b)%i(c, i_idx)
                end do
            end do

            call this%network_manager%advanceVoltage()

            do i = 1, size(this%network_manager%networks)
                do j = 1, size(this%network_manager%networks(i)%nodes)
                    b = this%network_manager%networks(i)%nodes(j)%bundle_number
                    c = this%network_manager%networks(i)%nodes(j)%conductor_number
                    if (.not. this%network_manager%networks(i)%nodes(j)%open) then 
                        v_idx = this%network_manager%networks(i)%nodes(j)%v_index
                        i_idx = this%network_manager%networks(i)%nodes(j)%i_index
                        if (this%bundles(b)%bundle_in_layer) this%bundles(b)%v(c, v_idx) = this%network_manager%networks(i)%nodes(j)%v
                    else 
                        if (this%network_manager%networks(i)%nodes(j)%side == TERMINAL_NODE_SIDE_INI) then 
                            this%bundles(b)%v(c,1) = this%bundles(b)%v(c,1) - 2*dot_product(this%bundles(b)%i_diff(1,c,:), this%bundles(b)%i(:,1))
                        else if (this%network_manager%networks(i)%nodes(j)%side == TERMINAL_NODE_SIDE_END) then 
                            n = this%bundles(b)%number_of_divisions
                            this%bundles(b)%v(c,n+1) = this%bundles(b)%v(c,n+1) + 2*dot_product(this%bundles(b)%i_diff(n,c,:), this%bundles(b)%i(:,n))
                        end if
                    end if
                end do
            end do
        end if
    end subroutine

    subroutine advanceBundlesCurrent(this)
        class(mtln_t) :: this
        integer :: i
        do i = 1, this%number_of_bundles
            if (this%bundles(i)%bundle_in_layer) call this%bundles(i)%advanceCurrent()

        end do
    end subroutine

    subroutine advanceTime(this)
        class(mtln_t) :: this
        this%time = this%time + this%dt
    end subroutine

    subroutine updateProbes(this)
        class(mtln_t) :: this
        integer :: i, j
        do i = 1, this%number_of_bundles
            if (size(this%bundles(i)%probes) /= 0 .and. this%bundles(i)%bundle_in_layer) then 
                do j = 1, size(this%bundles(i)%probes)
                    if (this%bundles(i)%probes(j)%in_layer) call this%bundles(i)%probes(j)%update(this%time, this%bundles(i)%v, this%bundles(i)%i)
                end do 
            end if
        end do
    end subroutine

    function getTimeRange(this, time) result(res)
        class(mtln_t) :: this
        real, intent(in), optional :: time
        integer :: res
        if (present(time)) then 
            res =  floor(time / this%dt)
        else
            res =  floor(this%final_time / this%dt)
        end if
    end function

    subroutine updateBundlesTimeStep(this, dt)
        class(mtln_t) :: this
        real :: dt
        integer :: i
        do i = 1, this%number_of_bundles
            this%bundles(i)%dt = dt
        end do
    end subroutine

    subroutine updatePULTerms(this)
        class(mtln_t) :: this
        integer :: i, j 
        do i = 1, this%number_of_bundles
            if (this%bundles(i)%bundle_in_layer) then
                call this%bundles(i)%updateLRTerms()
                call this%bundles(i)%updateCGTerms()
                do j = 1, size(this%bundles(i)%probes)
                    call this%bundles(i)%probes(j)%resizeFrames(this%getTimeRange(this%final_time), & 
                                                                this%bundles(i)%number_of_conductors)
                end do
            end if
        end do
   
    end subroutine


    subroutine runUntil(this, final_time)
        class(mtln_t) :: this
        real, intent(in):: final_time
        real :: time
        integer :: i

        do i = 1, this%getTimeRange(final_time)
            call this%advanceBundlesVoltage()
            call this%advanceNWVoltage()
            call this%advanceBundlesCurrent()
            call this%advanceTime()
            call this%updateProbes()
        end do

    end subroutine

    subroutine mtln_run(this)
        class(mtln_t) :: this
        real :: time
        integer :: i

        call this%updatePULTerms()
        do i = 1, this%getTimeRange(this%final_time)
            call this%advanceBundlesVoltage()
            call this%advanceNWVoltage()
            call this%advanceBundlesCurrent()
            call this%advanceTime()
            call this%updateProbes()
        end do

    end subroutine

end module mtln_solver_mod