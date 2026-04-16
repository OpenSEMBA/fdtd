module mtln_solver_m 

    use mtl_bundle_m
    use network_manager_m
    use mtln_preprocess_m
    use FDETYPES_m, only: XYZlimit_t, RKIND_TIEMPO
#ifdef CompileWithMPI
    use FDETYPES_m, only: SUBCOMM_MPI, REALSIZE, INTEGERSIZE, MPI_STATUS_SIZE
#endif
    implicit none


    type, public :: mtln_t
        real(kind=RKIND_TIEMPO) :: time, dt, final_time
        type(mtl_bundle_t), allocatable, dimension(:) :: bundles
        type(network_manager_t) :: network_manager
        ! type(probe_t), allocatable, dimension(:) :: probes
        integer :: number_of_bundles
        integer :: number_of_steps
        real(kind=rkind) :: null_field

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
        procedure :: initObservation => mtln_initObservation
        procedure :: updateObservation => mtln_updateObservation
        procedure :: closeObservation => mtln_closeObservation
    end type mtln_t


    interface mtln_t
        module procedure mtlnCtor
    end interface
    

contains

    function mtlnCtor(parsed, alloc) result(res)
        type(parsed_mtln_t) :: parsed
        type(XYZlimit_t), dimension(1:6), intent(in), optional :: alloc
        type(mtln_t) :: res
        integer :: i
        type(preprocess_t) :: pre

#ifdef CompileWithMPI
        integer(kind=4) :: sizeof, ierr
#endif

#ifdef CompileWithMPI
        call mpi_barrier(subcomm_mpi, ierr)
#endif
        if (present(alloc)) then 
            pre = preprocess(parsed, alloc)
        else  
            pre = preprocess(parsed)
        end if
        if (size(pre%bundles) == 0) then
            res%number_of_bundles = 0
            return
        end if
        
        res%dt = pre%dt
        res%time  = 0.0
        res%final_time = pre%final_time
        
        res%bundles = pre%bundles
        res%number_of_bundles = size(res%bundles)
        
        res%network_manager = pre%network_manager
        ! res%probes = pre%probes
        call res%updateBundlesTimeStep(res%dt)
        call res%initNodes()

        res%number_of_steps = parsed%number_of_steps
        res%null_field = 0.0_rkind
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
        call this%updateProbes()
        call this%advanceTime()

    end subroutine

    subroutine step_alone(this)
        class(mtln_t) :: this
        integer :: i 

        call this%advanceBundlesVoltage()
        call this%advanceNWVoltage()
        call this%advanceBundlesCurrent()

        call this%updateProbes()
        call this%advanceTime()

    end subroutine

    subroutine setExternalLongitudinalField(this)
        class(mtln_t) :: this
        integer :: i
#ifdef CompileWithMPI
        integer :: ierr
      call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif
        do i = 1, this%number_of_bundles
            if (this%bundles(i)%bundle_in_layer) call this%bundles(i)%setExternalLongitudinalField()
        end do

    end subroutine

    subroutine advanceBundlesVoltage(this)
        class(mtln_t) :: this
        integer :: i

        do i = 1, this%number_of_bundles
            if (this%bundles(i)%bundle_in_layer) then
                call this%bundles(i)%updateGenerators(this%time, this%dt)
                call this%bundles(i)%advanceVoltage()
            end if
        end do

    end subroutine

    subroutine advanceNWVoltage(this)
        class(mtln_t) :: this
        integer :: i,j
        integer ::b, c, v_idx, i_idx
        integer :: n
! #ifdef CompileWithMPI
!         integer(kind=4) :: ierr
!         call mpi_barrier(subcomm_mpi, ierr)
! #endif
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
#ifdef CompileWithMPI
        integer :: ierr
        call mpi_barrier(subcomm_mpi, ierr)
#endif
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
                    if (this%bundles(i)%probes(j)%in_layer) then
                        call this%bundles(i)%probes(j)%update(this%time, this%bundles(i)%v, this%bundles(i)%i)
                    end if
                end do 
            end if
        end do
    end subroutine

    function getTimeRange(this, time) result(res)
        class(mtln_t) :: this
        real(kind=RKIND_TIEMPO), intent(in), optional :: time
        integer :: res
        if (present(time)) then 
            res =  floor(time / this%dt)
        else
            res =  floor(this%final_time / this%dt)
        end if
    end function

    subroutine updateBundlesTimeStep(this, dt)
        class(mtln_t) :: this
        real(kind=RKIND_TIEMPO) :: dt
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
        real(kind=RKIND_TIEMPO), intent(in):: final_time
        real(kind=RKIND_TIEMPO) :: time
        integer :: i

        do i = 0, this%getTimeRange(final_time)
            call this%advanceBundlesVoltage()
            call this%advanceNWVoltage()
            call this%advanceBundlesCurrent()
            call this%updateProbes()
            call this%advanceTime()
            call this%updateObservation(i)
        end do

    end subroutine

    subroutine mtln_run(this)
        class(mtln_t) :: this
        real(kind=RKIND_TIEMPO) :: time
        integer :: i

        do i = 0, this%getTimeRange(this%final_time)
            call this%advanceBundlesVoltage()
            call this%advanceNWVoltage()
            call this%advanceBundlesCurrent()
            call this%updateProbes()
            call this%advanceTime()
            call this%updateObservation(i)
        end do

    end subroutine

    subroutine mtln_initObservation(this, nEntradaRoot)
        class(mtln_t) :: this
        character(len=*), intent(in) :: nEntradaRoot
        integer :: i, j, k, unit
        character(len=bufsize) :: path
        character(len=bufsize) :: temp
        character(len=:), allocatable :: buffer

        if (.not. allocated(this%bundles)) return
        unit = 2000
        do i = 1, size(this%bundles)
            do j = 1, size(this%bundles(i)%probes)
#ifdef CompileWithMPI
                if (.not. this%bundles(i)%probes(j)%in_layer) cycle
#endif          
                this%bundles(i)%probes(j)%unit = unit
                path = trim(trim(nEntradaRoot)//"_"//trim(this%bundles(i)%probes(j)%name)//".dat")
                open (unit=unit, file=trim(path))
                write (*, *) 'name: ', trim(this%bundles(i)%probes(j)%name)
                buffer = "time"
                do k = 1, size(this%bundles(i)%probes(j)%val, 2)
                    write (temp, *) k
                    buffer = buffer//" "//"conductor_"//trim(adjustl(temp))
                end do
                write (unit, '(a)') trim(buffer)
                unit = unit + 1
            end do
        end do
    end subroutine

    subroutine mtln_updateObservation(this, step)
        class(mtln_t) :: this
        integer, intent(in) :: step
        integer :: i, j, k, n
        integer :: unit
        character(len=bufsize) :: temp
        character(len=bufsize) :: path
        character(len=:), allocatable :: buffer
#ifdef CompileWithMPI
        integer(kind=4) :: ierr
#endif
        if (.not. allocated(this%bundles)) return
        do i = 1, size(this%bundles)
            do j = 1, size(this%bundles(i)%probes)
#ifdef CompileWithMPI
                if (.not. this%bundles(i)%probes(j)%in_layer) cycle
#endif          
                buffer = ""
                write(temp, *) this%bundles(i)%probes(j)%t(step+1)
                buffer = buffer//trim(temp)
                do n = 1, size(this%bundles(i)%probes(j)%val, 2)
                    write (temp, *) this%bundles(i)%probes(j)%val(step+1, n)
                    buffer = buffer//" "//trim(temp)
                end do
                write (this%bundles(i)%probes(j)%unit, '(a)') trim(buffer)
            end do
        end do
    end subroutine

    subroutine mtln_closeObservation(this)
        class(mtln_t) :: this
        integer :: i, j, k
        if (.not. allocated(this%bundles)) return
        do i = 1, size(this%bundles)
            do j = 1, size(this%bundles(i)%probes)
#ifdef CompileWithMPI
                if (.not. this%bundles(i)%probes(j)%in_layer) cycle
#endif          
                close(this%bundles(i)%probes(k)%unit)
            end do
      end do
    end subroutine

    
end module mtln_solver_m