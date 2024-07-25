module mtln_solver_mod 

    use mtl_bundle_mod
    use network_manager_mod
    use preprocess_mod
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

    function mtlnCtor(parsed) result(res)
        type(parsed_mtln_t) :: parsed
        type(mtln_t) :: res
        integer :: i
        type(preprocess_t) :: pre

        pre = preprocess(parsed)
        if (size(pre%bundles) == 0) then
            res%number_of_bundles = 0
            return
        end if

        res%has_multiwires = parsed%has_multiwires

        res%dt = pre%dt
        res%time  = 0.0
        res%final_time = pre%final_time

        res%bundles = pre%bundles
        res%network_manager = pre%network_manager
        res%number_of_bundles = size(res%bundles)
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
            call this%bundles(i)%setExternalLongitudinalField()
        end do

    end subroutine

    subroutine advanceBundlesVoltage(this)
        class(mtln_t) :: this
        integer :: i
        do i = 1, this%number_of_bundles
            call this%bundles(i)%updateSources(this%time, this%dt)
            call this%bundles(i)%advanceVoltage()
        end do

    end subroutine

    subroutine advanceNWVoltage(this)
        class(mtln_t) :: this
        integer :: i,j
        integer ::b, c, v_idx, i_idx
            
        do i = 1, size(this%network_manager%networks)
            do j = 1, size(this%network_manager%networks(i)%nodes)
                b = this%network_manager%networks(i)%nodes(j)%bundle_number
                c = this%network_manager%networks(i)%nodes(j)%conductor_number
                v_idx = this%network_manager%networks(i)%nodes(j)%v_index
                i_idx = this%network_manager%networks(i)%nodes(j)%i_index

                this%network_manager%networks(i)%nodes(j)%i = this%bundles(b)%i(c, i_idx)
            end do
        end do

        call this%network_manager%advanceVoltage()

        do i = 1, size(this%network_manager%networks)
            do j = 1, size(this%network_manager%networks(i)%nodes)
                b = this%network_manager%networks(i)%nodes(j)%bundle_number
                c = this%network_manager%networks(i)%nodes(j)%conductor_number
                v_idx = this%network_manager%networks(i)%nodes(j)%v_index
                i_idx = this%network_manager%networks(i)%nodes(j)%i_index

                this%bundles(b)%v(c, v_idx) = this%network_manager%networks(i)%nodes(j)%v
            end do
        end do

    end subroutine

    subroutine advanceBundlesCurrent(this)
        class(mtln_t) :: this
        integer :: i
        do i = 1, this%number_of_bundles
            call this%bundles(i)%advanceCurrent()

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
            if (size(this%bundles(i)%probes) /= 0) then 
                do j = 1, size(this%bundles(i)%probes)
                    call this%bundles(i)%probes(j)%update(this%time, this%bundles(i)%v, this%bundles(i)%i)
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
            call this%bundles(i)%updateLRTerms()
            call this%bundles(i)%updateCGTerms()
            do j = 1, size(this%bundles(i)%probes)
                call this%bundles(i)%probes(j)%resizeFrames(this%getTimeRange(this%final_time), & 
                                                            this%bundles(i)%number_of_conductors)
            end do
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