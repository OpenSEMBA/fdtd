module termination_handler_m

    use mtln_types_m, only: TERMINATION_SHORT, TERMINATION_OPEN, TERMINATION_SERIES, &
         TERMINATION_PARALLEL, TERMINATION_RsLCp, TERMINATION_RLsCp, &
         TERMINATION_LsRCp, TERMINATION_CsLRp, TERMINATION_RCsLp, TERMINATION_LCsRp, &
         TERMINATION_CIRCUIT, TERMINATION_NETWORK, TERMINATION_UNDEFINED
    use FDETYPES_m, only: RKIND
    implicit none

    ! Simple termination handler for RLC-type terminations that don't need ngspice
    ! These terminations can be computed directly in Fortran using lumped element models

    type, public :: simple_termination_t
        integer :: termination_type
        real(kind=rkind) :: R, L, C  ! Component values
        real(kind=rkind) :: v_node    ! Node voltage (output)
        real(kind=rkind) :: i_node    ! Node current from MTL (input)
        real(kind=rkind) :: i_prev    ! Previous current (for numerical differentiation)
        real(kind=rkind) :: v_prev    ! Previous voltage (for numerical integration)
        real(kind=rkind) :: dt        ! Time step
        real(kind=rkind) :: i_l       ! Inductor current state (for parallel L)
        real(kind=rkind) :: v_c       ! Capacitor voltage state (for series C)
    contains
        procedure :: init => termination_init
        procedure :: step => termination_step
    end type simple_termination_t

contains

    subroutine termination_init(this, term_type, R, L, C, dt)
        class(simple_termination_t) :: this
        integer, intent(in) :: term_type
        real(kind=rkind), intent(in) :: R, L, C
        real(kind=rkind), intent(in) :: dt

        this%termination_type = term_type
        this%R = R
        this%L = L
        this%C = C
        this%dt = dt
        this%v_node = 0.0_rkind
        this%i_node = 0.0_rkind
        this%i_prev = 0.0_rkind
        this%v_prev = 0.0_rkind
        this%i_l = 0.0_rkind
        this%v_c = 0.0_rkind
    end subroutine termination_init

    subroutine termination_step(this, i_in, dt)
        class(simple_termination_t) :: this
        real(kind=rkind), intent(in) :: i_in
        real(kind=rkind), intent(in), optional :: dt
        real(kind=rkind) :: di_dt, dv_dt
        real(kind=rkind) :: G, Y_eq, I_eq, V_new

        if (present(dt)) this%dt = dt

        ! Save previous values
        this%v_prev = this%v_node
        this%i_prev = this%i_node
        this%i_node = i_in

        ! Compute derivative of current (for inductor voltage)
        di_dt = (this%i_node - this%i_prev) / this%dt

        select case (this%termination_type)

        case (TERMINATION_SHORT)
            ! Short circuit: V = 0 (ideal) or V = I * R_parasitic
            this%v_node = this%i_node * this%R

        case (TERMINATION_OPEN)
            ! Open circuit with parasitic C and R
            ! dV/dt = I/C - V/(R*C) using trapezoidal
            if (this%C > 0.0_rkind .and. this%C < 1e20_rkind) then
                G = 1.0_rkind / this%R
                Y_eq = this%C / this%dt + 0.5_rkind * G
                I_eq = this%i_node + this%C * this%v_prev / this%dt - 0.5_rkind * G * this%v_prev
                this%v_node = I_eq / Y_eq
            else
                this%v_node = 0.0_rkind
            end if

        case (TERMINATION_SERIES)
            ! Series R-L or R-L-C: V = I*R + L*dI/dt + V_C
            this%v_node = this%i_node * this%R + this%L * di_dt
            if (this%C > 0.0_rkind .and. this%C < 1e20_rkind) then
                ! Add series capacitor voltage (trapezoidal integration)
                this%v_c = this%v_c + this%dt / (2.0_rkind * this%C) * &
                     (this%i_node + this%i_prev)
                this%v_node = this%v_node + this%v_c
            end if

        case (TERMINATION_PARALLEL)
            ! Parallel R-L-C using companion model (trapezoidal integration)
            ! I = V/R + I_L + C*dV/dt
            ! I_L,n = I_L,n-1 + dt/2 * (V_n + V_n-1)/L
            ! C*dV/dt = C/dt * (V_n - V_n-1)
            ! => I = V_n/R + I_L,n-1 + dt*V_n/(2L) + dt*V_n-1/(2L) + C*(V_n - V_n-1)/dt
            ! => V_n * (1/R + dt/(2L) + C/dt) = I - I_L,n-1 - dt*V_n-1/(2L) + C*V_n-1/dt
            call solve_parallel_rlc(this%R, this%L, this%C, this%dt, &
                 this%i_node, this%v_prev, this%i_l, this%v_node)

        case (TERMINATION_RsLCp)
            ! R series, (L||C) parallel
            call solve_parallel_lc(this%L, this%C, this%dt, &
                 this%i_node, this%v_prev, this%i_l, V_new)
            this%v_node = V_new + this%i_node * this%R

        case (TERMINATION_RLsCp)
            ! (R+L) series, C parallel
            call solve_parallel_c(this%C, this%dt, this%i_node, this%v_prev, V_new)
            this%v_node = V_new + this%i_node * this%R + this%L * di_dt

        case (TERMINATION_LsRCp)
            ! L series, (R||C) parallel
            call solve_parallel_rc(this%R, this%C, this%dt, this%i_node, this%v_prev, V_new)
            this%v_node = V_new + this%L * di_dt

        case (TERMINATION_CsLRp)
            ! C series, (L||R) parallel
            call solve_parallel_lr(this%L, this%R, this%dt, this%i_node, this%v_prev, this%i_l, V_new)
            ! Add series C voltage
            this%v_c = this%v_c + this%dt / (2.0_rkind * this%C) * &
                 (this%i_node + this%i_prev)
            this%v_node = V_new + this%v_c

        case (TERMINATION_RCsLp)
            ! (R+C) series, L parallel
            call solve_parallel_l(this%L, this%dt, this%i_node, this%v_prev, this%i_l, V_new)
            ! Add series R+C
            this%v_c = this%v_c + this%dt / (2.0_rkind * this%C) * &
                 (this%i_node + this%i_prev)
            this%v_node = V_new + this%i_node * this%R + this%v_c

        case (TERMINATION_LCsRp)
            ! (L+C) series, R parallel
            call solve_parallel_r(this%R, this%i_node, V_new)
            ! Add series L+C
            this%v_c = this%v_c + this%dt / (2.0_rkind * this%C) * &
                 (this%i_node + this%i_prev)
            this%v_node = V_new + this%L * di_dt + this%v_c

        case default
            ! Unknown type - should not happen
            this%v_node = 0.0_rkind

        end select

    end subroutine termination_step

    ! Solve parallel RLC using companion model (trapezoidal integration)
    subroutine solve_parallel_rlc(R, L, C, dt, I_in, V_prev, I_L_prev, V_out)
        real(kind=rkind), intent(in) :: R, L, C, dt, I_in, V_prev
        real(kind=rkind), intent(inout) :: I_L_prev
        real(kind=rkind), intent(out) :: V_out
        real(kind=rkind) :: G, Y_eq, I_eq, denom

        G = 0.0_rkind
        if (R > 0.0_rkind .and. R < 1e20_rkind) G = 1.0_rkind / R

        denom = G + C / dt
        if (L > 0.0_rkind .and. L < 1e20_rkind) denom = denom + dt / (2.0_rkind * L)

        I_eq = I_in + C * V_prev / dt - G * V_prev - I_L_prev
        if (L > 0.0_rkind .and. L < 1e20_rkind) then
            I_eq = I_eq - dt * V_prev / (2.0_rkind * L)
        end if

        V_out = I_eq / denom

        ! Update inductor current
        if (L > 0.0_rkind .and. L < 1e20_rkind) then
            I_L_prev = I_L_prev + dt / (2.0_rkind * L) * (V_out + V_prev)
        end if
    end subroutine solve_parallel_rlc

    ! Solve parallel LC
    subroutine solve_parallel_lc(L, C, dt, I_in, V_prev, I_L_prev, V_out)
        real(kind=rkind), intent(in) :: L, C, dt, I_in, V_prev
        real(kind=rkind), intent(inout) :: I_L_prev
        real(kind=rkind), intent(out) :: V_out
        real(kind=rkind) :: denom, I_eq

        denom = C / dt + dt / (2.0_rkind * L)
        I_eq = I_in + C * V_prev / dt - I_L_prev - dt * V_prev / (2.0_rkind * L)
        V_out = I_eq / denom
        I_L_prev = I_L_prev + dt / (2.0_rkind * L) * (V_out + V_prev)
    end subroutine solve_parallel_lc

    ! Solve parallel C only
    subroutine solve_parallel_c(C, dt, I_in, V_prev, V_out)
        real(kind=rkind), intent(in) :: C, dt, I_in, V_prev
        real(kind=rkind), intent(out) :: V_out
        V_out = (I_in + C * V_prev / dt) / (C / dt)
    end subroutine solve_parallel_c

    ! Solve parallel RC
    subroutine solve_parallel_rc(R, C, dt, I_in, V_prev, V_out)
        real(kind=rkind), intent(in) :: R, C, dt, I_in, V_prev
        real(kind=rkind), intent(out) :: V_out
        real(kind=rkind) :: G, Y_eq, I_eq

        G = 1.0_rkind / R
        Y_eq = G + C / dt
        I_eq = I_in + C * V_prev / dt - G * V_prev
        V_out = I_eq / Y_eq
    end subroutine solve_parallel_rc

    ! Solve parallel LR
    subroutine solve_parallel_lr(L, R, dt, I_in, V_prev, I_L_prev, V_out)
        real(kind=rkind), intent(in) :: L, R, dt, I_in, V_prev
        real(kind=rkind), intent(inout) :: I_L_prev
        real(kind=rkind), intent(out) :: V_out
        real(kind=rkind) :: G, Y_eq, I_eq

        G = 1.0_rkind / R
        Y_eq = G + dt / (2.0_rkind * L)
        I_eq = I_in - G * V_prev - I_L_prev - dt * V_prev / (2.0_rkind * L)
        V_out = I_eq / Y_eq
        I_L_prev = I_L_prev + dt / (2.0_rkind * L) * (V_out + V_prev)
    end subroutine solve_parallel_lr

    ! Solve parallel L only
    subroutine solve_parallel_l(L, dt, I_in, V_prev, I_L_prev, V_out)
        real(kind=rkind), intent(in) :: L, dt, I_in, V_prev
        real(kind=rkind), intent(inout) :: I_L_prev
        real(kind=rkind), intent(out) :: V_out
        real(kind=rkind) :: Y_eq, I_eq

        Y_eq = dt / (2.0_rkind * L)
        I_eq = I_in - I_L_prev - dt * V_prev / (2.0_rkind * L)
        V_out = I_eq / Y_eq
        I_L_prev = I_L_prev + dt / (2.0_rkind * L) * (V_out + V_prev)
    end subroutine solve_parallel_l

    ! Solve parallel R only
    subroutine solve_parallel_r(R, I_in, V_out)
        real(kind=rkind), intent(in) :: R, I_in
        real(kind=rkind), intent(out) :: V_out
        V_out = I_in * R
    end subroutine solve_parallel_r

end module termination_handler_m
