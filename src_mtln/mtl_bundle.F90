module mtl_bundle_mod

    use utils_mod
    use probes_mod
    use dispersive_mod
    use mtl_mod
    implicit none

    type, public :: mtl_bundle_t
        character (len=:), allocatable :: name
        real, allocatable, dimension(:,:,:) :: lpul, cpul, rpul, gpul
        integer  :: number_of_conductors = 0, number_of_divisions = 0
        real, dimension(:), allocatable :: step_size
        real, allocatable, dimension(:,:) :: v, i
        real, allocatable, dimension(:,:,:) :: du(:,:,:)
        real :: time = 0.0, dt = 1e10
        type(probe_t), allocatable, dimension(:) :: probes
        type(transfer_impedance_t) :: transfer_impedance
        integer, dimension(:), allocatable :: conductors_in_level
        
        real, dimension(:,:,:), allocatable :: v_term, i_term
        real, dimension(:,:,:), allocatable :: v_diff, i_diff

        type(external_field_segment_t), dimension(:), allocatable :: external_field_segments


    contains
        procedure :: mergePULMatrices
        procedure :: mergeDispersiveMatrices
        procedure :: initialAllocation
        procedure :: addProbe
        procedure :: updateLRTerms
        procedure :: updateCGTerms
        ! procedure :: get_phase_velocities
        procedure :: updateSources => bundle_updateSources
        procedure :: advanceVoltage => bundle_advanceVoltage
        procedure :: advanceCurrent => bundle_advanceCurrent
        procedure :: addTransferImpedance => bundle_addTransferImpedance
        ! procedure :: setConnectorTransferImpedance
        procedure :: setExternalVoltage => bundle_setExternalVoltage
        procedure :: updateExternalCurrent => bundle_updateExternalCurrent

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
        res%step_size = levels(1)%lines(1)%step_size
        res%number_of_divisions = size(res%step_size,1)
        res%external_field_segments = levels(1)%lines(1)%external_field_segments
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

    type(probe_t) function addProbe(this, index, probe_type) result(res)
        class(mtl_bundle_t) :: this
        integer, intent(in) :: index
        integer, intent(in) :: probe_type
        res = probeCtor(index, probe_type, this%dt)
        this%probes = [this%probes, res]
    end function

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

        do i = 1, this%number_of_divisions
            F1(i,:,:) = matmul(this%du(i,:,:), &
                               this%lpul(i,:,:)/this%dt + 0.5*this%transfer_impedance%d(i,:,:) + this%transfer_impedance%e(i,:,:)/this%dt + 0.5*this%rpul(i,:,:) + this%transfer_impedance%q1_sum(i,:,:))
            F2(i,:,:) = matmul(this%du(i,:,:), &
                               this%lpul(i,:,:)/this%dt - 0.5*this%transfer_impedance%d(i,:,:) + this%transfer_impedance%e(i,:,:)/this%dt - 0.5*this%rpul(i,:,:) - this%transfer_impedance%q1_sum(i,:,:))
            IF1(i,:,:) = inv(F1(i,:,:))
            ! this%i_term(i,:,:) = matmul(IF1(i,:,:),F2(i,:,:))
            ! this%v_diff(i,:,:) = IF1(i,:,:)
        enddo
        F1 = reshape(source=F1, shape=[this%number_of_divisions,this%number_of_conductors, this%number_of_conductors],order=[2,3,1])
        F2 = reshape(source=F2, shape=[this%number_of_divisions,this%number_of_conductors, this%number_of_conductors],order=[2,3,1])
        IF1 = reshape(source=IF1,shape=[this%number_of_divisions,this%number_of_conductors, this%number_of_conductors])
        do i = 1, this%number_of_divisions
            this%i_term(i,:,:) = matmul(IF1(i,:,:),F2(i,:,:))
            this%v_diff(i,:,:) = IF1(i,:,:)
        enddo
        this%i_term = reshape(source=this%i_term, shape=[this%number_of_divisions,this%number_of_conductors, this%number_of_conductors],order=[2,3,1])
        this%v_diff = reshape(source=this%v_diff, shape=[this%number_of_divisions,this%number_of_conductors, this%number_of_conductors],order=[2,3,1])
        write(*,*) size(F1,1), ' ',size(F1,2), ' ', size(F1,3)
        write(*,*) F1(1,1,1)
        write(*,*) size(F2,1), ' ',size(F2,2), ' ', size(F2,3)
        write(*,*) F2(1,1,1)
        write(*,*) size(IF1,1), ' ',size(IF1,2), ' ', size(IF1,3)
        write(*,*) IF1(1,1,1)
        write(*,*) size(this%i_term,1), ' ',size(this%i_term,2), ' ', size(this%i_term,3)
        write(*,*) this%i_term(1,1,1)
        write(*,*) size(this%v_diff,1), ' ',size(this%v_diff,2), ' ', size(this%v_diff,3)
        write(*,*) this%v_diff(1,1,1)

        ! ! IF1 = reshape(source=[(inv(F1(i,:,:)), i = 1, this%number_of_divisions)], &
        ! !         shape=[this%number_of_divisions,this%number_of_conductors, this%number_of_conductors], &
        ! !         order=[2,3,1])
        ! this%i_term = reshape(&
        !     source=[(matmul(IF1(i,:,:), F2(i,:,:)), i = 1, this%number_of_divisions)], &
        !     shape=[this%number_of_divisions,this%number_of_conductors, this%number_of_conductors], &
        !     order=[2,3,1])
        ! this%v_diff = IF1

        ! F1 = reshape(source=[(matmul( &
        !     this%du(i,:,:), &
        !     this%lpul(i,:,:)/this%dt + &
        !         0.5*this%transfer_impedance%d(i,:,:) + &
        !         this%transfer_impedance%e(i,:,:)/this%dt + &
        !         0.5*this%rpul(i,:,:) + &
        !         this%transfer_impedance%q1_sum(i,:,:)), &
        !     i = 1,this%number_of_divisions)], & 
        !     shape=[this%number_of_divisions,this%number_of_conductors, this%number_of_conductors], &
        !     order=[2,3,1])
        ! write(*,*) F1(1,1,1)
        ! F2 = reshape(source=[(matmul( &
        !     this%du(i,:,:), &
        !     this%lpul(i,:,:)/this%dt - &
        !     0.5*this%transfer_impedance%d(i,:,:) + &
        !     this%transfer_impedance%e(i,:,:)/this%dt - &
        !     0.5*this%rpul(i,:,:) - &
        !     this%transfer_impedance%q1_sum(i,:,:)), &
        !     i = 1,this%number_of_divisions)], & 
        !     shape=[this%number_of_divisions,this%number_of_conductors, this%number_of_conductors], &
        !     order=[2,3,1])
        ! write(*,*) F2(1,1,1)
        ! IF1 = reshape(source=[(inv(F1(i,:,:)), i = 1, 100)], &
        !             shape=[100,1, 1], &
        !             order=[2,3,1])
        ! ! IF1 = reshape(source=[(inv(F1(i,:,:)), i = 1, this%number_of_divisions)], &
        ! !             shape=[this%number_of_divisions,this%number_of_conductors, this%number_of_conductors], &
        ! !             order=[2,3,1])
        ! this%i_term = reshape(&
        !     source=[(matmul(IF1(i,:,:), F2(i,:,:)), i = 1, this%number_of_divisions)], &
        !     shape=[this%number_of_divisions,this%number_of_conductors, this%number_of_conductors], &
        !     order=[2,3,1])
        ! this%v_diff = IF1

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

    end subroutine

    subroutine bundle_updateSources(this, time, dt)
        class(mtl_bundle_t) ::this
        real, intent(in) :: time, dt
        !TODO
    end subroutine

    subroutine bundle_advanceVoltage(this)
        class(mtl_bundle_t) ::this
        integer :: i

        do i = 2, this%number_of_divisions
            this%v(:, i) = matmul(this%v_term(i,:,:), this%v(:,i)) - &
                           matmul(this%i_diff(i,:,:), this%i(:,i) - this%i(:,i-1)  )
        end do


        ! this%v(:,2:this%number_of_divisions) = &
        !          reshape(source=[(matmul(this%v_term(i,:,:), this%v(:,i)), i = 2, this%number_of_divisions - 1)], &
        !                  shape = [this%number_of_divisions -1,this%number_of_conductors], &
        !                  order = [2,1])
    end subroutine

    subroutine bundle_advanceCurrent(this)
        class(mtl_bundle_t) ::this
        real, dimension(:,:), allocatable :: i_prev, i_now
        integer :: i
        ! call this%transfer_impedance%updateQ3Phi()
        ! i_prev = this%i
        do i = 1, this%number_of_divisions 
            this%i(:,i) = matmul(this%i_term(i,:,:), this%i(:,i)) - &
                          matmul(this%v_diff(i,:,:), (this%v(:,i+1) - this%v(:,i))) !- &
                                !  matmul(0.5*this%du_length(i,:,:), this%el))
                        !   matmul(this%v_diff(i,:,:), matmul(this%du(i,:,:), this%transfer_impedance%q3_phi(i,:)))
        enddo
        !TODO - revisar
        ! i_now = this%i
        ! call this%transfer_impedance%updatePhi(i_prev, i_now)
    end subroutine

    subroutine bundle_setExternalVoltage(this)
        class(mtl_bundle_t) :: this
        integer :: i
        do i = 1, size(this%v,2) 
            this%v(1, i) = this%external_field_segments(i)%Efield_main2wire * this%step_size(i)
        end do
    end subroutine

    subroutine bundle_updateExternalCurrent(this, current)
        class(mtl_bundle_t) :: this
        real, dimension(:,:,:), intent(inout) :: current
        ! real, dimension(:), intent(inout) :: current
        integer, dimension(:), allocatable :: position
        integer :: i
        do i = 1, size(this%i,2)
            position = this%external_field_segments(i)%position
            current(position(1), position(2), position(3)) = this%i(1,i)
        end do
        ! current(:) =  this%i(1,:)
    end subroutine

end module mtl_bundle_mod