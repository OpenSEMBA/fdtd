module rational_approximation_mod
    use mtln_types_mod, only: &
        transfer_impedance_per_meter_t, &
        TRANSFER_IMPEDANCE_DIRECTION_BOTH, &
        TRANSFER_IMPEDANCE_DIRECTION_INWARDS, &
        TRANSFER_IMPEDANCE_DIRECTION_OUTWARDS
    implicit none

    type :: pol_res_t
        complex, allocatable, dimension(:) :: q1,q2,q3
        real :: r, l
        integer :: number_of_poles
        integer :: direction
    end type

    interface pol_res_t
        module procedure pol_resCtor
    end interface

contains

    function pol_resCtor(model, dt) result(res)
        type(transfer_impedance_per_meter_t), intent(in) :: model
        real, intent(in) :: dt
        type(pol_res_t) :: res
        
        res%r = model%resistive_term
        res%l = model%inductive_term
        res%number_of_poles = size(model%poles)
        
        allocate(res%q1(res%number_of_poles), res%q2(res%number_of_poles), res%q3(res%number_of_poles))
        if (res%number_of_poles /= 0) then
        block
            type(complex), dimension(:), allocatable :: alpha, beta
            alpha = model%residues / model%poles
            beta =  model%poles * dt
            res%q1 = - (alpha/beta)*(exp(beta)-beta-1)
            res%q2 = - (alpha/beta)*(1+exp(beta)*(beta-1))
            res%q3 = - exp(beta)
        end block
        end if
        res%direction = model%direction
    end function


end module