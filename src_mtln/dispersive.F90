module dispersive_m
    use utils_m
    use rational_approximation_m
    use FDETYPES_m, only: RKIND, RKIND_TIEMPO
    implicit none
    
    type :: dispersive_t
        real(kind=RKIND_TIEMPO) :: dt
        integer :: number_of_divisions, number_of_conductors, number_of_poles
        real, allocatable :: u(:,:)
        real, allocatable, dimension(:,:,:) :: d, e
        complex, allocatable, dimension(:,:) :: q3_phi
        complex, allocatable, dimension(:,:,:) :: q1_sum, q2_sum
        complex, allocatable, dimension(:,:,:) :: phi
        complex, allocatable, dimension(:,:,:,:) :: q1, q2, q3
    contains
        procedure :: updateQ3Phi
        procedure :: updatePhi
        procedure, private :: increaseOrder
    end type dispersive_t
    
    interface dispersive_t
        module procedure dispersiveCtor
    end interface

    type, extends(dispersive_t) :: lumped_t
    contains
        procedure :: addDispersiveLumped
        procedure, private :: positionIsEmpty
        procedure, private :: addDispersiveLumpedInConductor
    end type lumped_t

    interface lumped_t
        module procedure lumpedCtor
    end interface

    type, extends(dispersive_t) :: transfer_impedance_t
    contains
        procedure :: addTransferImpedance
        procedure :: setTransferImpedance
        procedure, private :: addTransferImpedanceInConductors
        procedure, private :: setTransferImpedanceInConductors
    end type transfer_impedance_t

    interface transfer_impedance_t
        module procedure transferImpendaceCtor
    end interface


contains
    

    function dispersiveCtor(number_of_conductors, number_of_poles, number_of_divisions, dt) result(res)
        type(dispersive_t) :: res
        integer :: number_of_conductors, number_of_poles
        real(kind=RKIND_TIEMPO) :: dt
        integer, intent(in) :: number_of_divisions
        complex :: zero 
        res%dt = dt
        res%number_of_conductors = number_of_conductors
        res%number_of_poles = number_of_poles
        res%number_of_divisions = number_of_divisions
        zero%re = 0.0
        zero%im = 0.0
        allocate(res%phi(res%number_of_divisions, number_of_conductors,number_of_poles), source = zero)
        allocate(res%q1 (res%number_of_divisions, number_of_conductors,number_of_conductors, number_of_poles), source = zero)
        allocate(res%q2 (res%number_of_divisions, number_of_conductors,number_of_conductors, number_of_poles), source = zero)
        allocate(res%q3 (res%number_of_divisions, number_of_conductors,number_of_conductors, number_of_poles), source = zero)
        allocate(res%d  (res%number_of_divisions, number_of_conductors,number_of_conductors), source = 0.0)
        allocate(res%e  (res%number_of_divisions, number_of_conductors,number_of_conductors), source = 0.0)
        allocate(res%q1_sum (res%number_of_divisions, number_of_conductors,number_of_conductors), source = zero)
        allocate(res%q2_sum (res%number_of_divisions, number_of_conductors,number_of_conductors), source = zero)
        allocate(res%q3_phi (res%number_of_divisions, number_of_conductors), source = zero)

    end function

    subroutine updateQ3Phi(this)
        class(dispersive_t) :: this
        integer :: i_div,j,k

        this%q3_phi(:,:) = reshape(&
            source = [(dotmatrixmul(this%q3(i_div, :, :, :), this%phi(i_div, :, :)), i_div = 1, this%number_of_divisions)], &
            shape=[this%number_of_divisions, this%number_of_conductors], &
            order=[2,1])
    end subroutine

    subroutine updatePhi(this, i_prev, i_now)
        class(dispersive_t) :: this
        real, dimension(:,:), intent(in):: i_prev, i_now
        integer :: i_div,k

        do k = 1, this%number_of_poles
            do i_div = 1, this%number_of_divisions
                this%phi(i_div,:,k) = matmul(this%q1(i_div,:,:,k), i_now(:,i_div)) + &
                                      matmul(this%q2(i_div,:,:,k), i_prev(:,i_div)) + &
                                      matmul(this%q3(i_div,:,:,k), this%phi(i_div,:,k))
            end do
        end do

    end subroutine

    subroutine increaseOrder(this, number_of_poles)
        class(dispersive_t) :: this
        integer, intent(in) :: number_of_poles
        type(dispersive_t) :: new_dispersive

        new_dispersive = dispersive_t(this%number_of_conductors, number_of_poles, this%number_of_divisions, this%dt)
        new_dispersive%q1(:,:,:,1:this%number_of_poles) = this%q1
        new_dispersive%q2(:,:,:,1:this%number_of_poles) = this%q2
        new_dispersive%q3(:,:,:,1:this%number_of_poles) = this%q3
        new_dispersive%phi(:,:,1:this%number_of_poles) = this%phi
        call move_alloc(from=new_dispersive%q1, to=this%q1)
        call move_alloc(from=new_dispersive%q2, to=this%q2)
        call move_alloc(from=new_dispersive%q3, to=this%q3)
        call move_alloc(from=new_dispersive%phi, to=this%phi)
        this%number_of_poles = number_of_poles
    end subroutine

    
    function lumpedCtor(number_of_conductors, number_of_poles, number_of_divisions, dt) result(res)
        type(lumped_t) :: res
        integer :: number_of_conductors, number_of_poles
        real(kind=RKIND_TIEMPO) :: dt
        integer, intent(in) :: number_of_divisions
        res%dispersive_t = dispersiveCtor(number_of_conductors, number_of_poles, number_of_divisions, dt)
    end function 

    function positionIsEmpty(this, index, conductor) result(res)
        class(lumped_t) :: this
        integer, intent(in) :: index, conductor
        logical :: res
        res = .true.
        if ((this%d(index, conductor, conductor) /= 0.0).or.&
            (this%e(index, conductor, conductor) /= 0.0).or.&
            .not.(all(this%q1(index, conductor, conductor,:) == 0.0)) .or.&
            .not.(all(this%q2(index, conductor, conductor,:) == 0.0)) .or.&
            .not.(all(this%q3(index, conductor, conductor,:) == 0.0))) then
                res = .false.
        end if
    end function


    subroutine addDispersiveLumped(this, index, conductor, model)
        class(lumped_t) :: this
        integer, intent(in) :: index
        integer, intent(in) :: conductor
        type(transfer_impedance_per_meter_t), intent(in) :: model
        type(pol_res_t) :: connector

        if (.not.this%positionIsEmpty(index, conductor))then
            error stop 'Dispersive connector already in conductor at position'
        end if

        connector = pol_res_t(model, this%dt)
        if (connector%number_of_poles > this%number_of_poles) call this%increaseOrder(connector%number_of_poles)
        call this%addDispersiveLumpedInConductor(index, conductor, connector)
       
        this%q1_sum = sumQComponents(this%q1)
        this%q2_sum = sumQComponents(this%q2)
    end subroutine

    subroutine addDispersiveLumpedInConductor(this, index, conductor, connector)
        class(lumped_t) :: this
        integer, intent(in) :: index, conductor
        type(pol_res_t), intent(in) :: connector

        this%d(index, conductor, conductor) = this%d(index, conductor, conductor) + connector%r
        this%e(index, conductor, conductor) = this%e(index, conductor, conductor) + connector%l
        if (connector%number_of_poles /= 0) then
            this%q1(index, conductor, conductor,:) = this%q1(index, conductor, conductor,:) - connector%q1(:)
            this%q2(index, conductor, conductor,:) = this%q2(index, conductor, conductor,:) - connector%q2(:)
            this%q3(index, conductor, conductor,:) = this%q3(index, conductor, conductor,:) - connector%q3(:)
        end if
    end subroutine


    function transferImpendaceCtor(number_of_conductors, number_of_poles, number_of_divisions, dt) result(res)
        type(transfer_impedance_t) :: res
        integer :: number_of_conductors, number_of_poles
        real(kind=RKIND_TIEMPO) :: dt
        integer :: number_of_divisions
        res%dispersive_t = dispersiveCtor(number_of_conductors, number_of_poles, number_of_divisions, dt)
    end function 

    function isCouplingInwards(direction) result(res)
        integer, intent(in) :: direction
        logical :: res
        res = .false.
        if (direction == TRANSFER_IMPEDANCE_DIRECTION_INWARDS &
            .or. direction == TRANSFER_IMPEDANCE_DIRECTION_BOTH) res = .true.
    end function

    function isCouplingOutWards(direction) result(res)
        integer, intent(in) :: direction
        logical :: res
        res = .false.
        if (direction == TRANSFER_IMPEDANCE_DIRECTION_OUTWARDS &
            .or. direction == TRANSFER_IMPEDANCE_DIRECTION_BOTH) res = .true.
    end function
    
    subroutine addTransferImpedance(this, conductor_out, range_in, model)
        class(transfer_impedance_t) :: this
        integer, intent(in) :: conductor_out
        integer, dimension(:), intent(in) :: range_in
        type(transfer_impedance_per_meter_t) :: model
        integer :: i
        type(pol_res_t) :: connector

        connector = pol_res_t(model, this%dt)
        if (connector%number_of_poles > this%number_of_poles) call this%increaseOrder(connector%number_of_poles)

        do i = 1, size(range_in)
            if (isCouplingInwards(connector%direction)) then
                call this%addTransferImpedanceInConductors(range_in(i), conductor_out, connector)
            end if
            if (isCouplingOutwards(connector%direction)) then
                call this%addTransferImpedanceInConductors(conductor_out, range_in(i), connector)
            end if
        end do

        this%q1_sum = sumQComponents(this%q1)
        this%q2_sum = sumQComponents(this%q2)
    end subroutine

    subroutine setTransferImpedance(this, index, conductor_out, range_in, model)
        class(transfer_impedance_t) :: this
        integer, intent(in) :: index
        integer, intent(in) :: conductor_out
        integer, dimension(:), intent(in) :: range_in
        type(transfer_impedance_per_meter_t) :: model
        integer :: i
        type(pol_res_t) :: connector

        connector = pol_res_t(model, this%dt)
        if (connector%number_of_poles > this%number_of_poles) call this%increaseOrder(connector%number_of_poles)

        do i = 1, size(range_in)
            if (isCouplingInwards(connector%direction)) then
                call this%setTransferImpedanceInConductors(index,range_in(i), conductor_out, connector)
            end if  
            if (isCouplingOutwards(connector%direction)) then
                call this%setTransferImpedanceInConductors(index, conductor_out, range_in(i), connector)
            end if
        end do

        this%q1_sum = sumQComponents(this%q1)
        this%q2_sum = sumQComponents(this%q2)
    end subroutine

    subroutine setTransferImpedanceInConductors(this, index, conductor_1, conductor_2, connector)
        class(transfer_impedance_t) :: this
        integer, intent(in) :: index, conductor_1, conductor_2
        type(pol_res_t), intent(in) :: connector

        this%d(index, conductor_1, conductor_2) = - connector%r
        this%e(index, conductor_1, conductor_2) = - connector%l
        if (connector%number_of_poles /= 0) then
            this%q1(index, conductor_1, conductor_2,:) = connector%q1(:)
            this%q2(index, conductor_1, conductor_2,:) = connector%q2(:)
            this%q3(index, conductor_1, conductor_2,:) = connector%q3(:)
        end if
    end subroutine

    subroutine addTransferImpedanceInConductors(this, conductor_1, conductor_2, connector)
        class(transfer_impedance_t) :: this
        integer, intent(in) :: conductor_1, conductor_2
        type(pol_res_t), intent(in) :: connector
        integer :: i

        this%d(:, conductor_1, conductor_2) = this%d(:, conductor_1, conductor_2) - connector%r
        this%e(:, conductor_1, conductor_2) = this%e(:, conductor_1, conductor_2) - connector%l

        if (connector%number_of_poles /= 0) then
            do i = 1, this%number_of_divisions
                this%q1(i,conductor_1, conductor_2,:) = this%q1(i,conductor_1, conductor_2,:) + connector%q1(:)
                this%q2(i,conductor_1, conductor_2,:) = this%q2(i,conductor_1, conductor_2,:) + connector%q2(:)
                this%q3(i,conductor_1, conductor_2,:) = this%q3(i,conductor_1, conductor_2,:) + connector%q3(:)
            end do
        end if
    end subroutine



end module dispersive_m


