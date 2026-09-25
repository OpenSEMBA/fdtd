module conformal_sgbc_m

   use FDETYPES_m, only: RKIND, RKIND_tiempo, SGBCMaterialProfile_t, EPSILON_VACUUM, MU_VACUUM
   use Report_m, only: WarnErrReport
   implicit none
   private

   type, public :: conformal_sgbc_state_t
      logical :: initialized = .false.
      integer :: n = 0
      real(kind=RKIND), allocatable :: e(:), e_old(:), h(:), rhs(:)
      real(kind=RKIND), allocatable :: e_ga(:), e_gb(:), h_ga(:), h_gb(:)
      real(kind=RKIND), allocatable :: lower_diag(:), diagonal(:), upper_diag(:)
      real(kind=RKIND), allocatable :: rhs_coeff(:,:)
      real(kind=RKIND), allocatable :: upper_prime(:), rhs_prime(:)
   contains
      procedure :: init => conformal_sgbc_init
      procedure :: advance => conformal_sgbc_advance
      procedure :: lower_e => conformal_sgbc_lower_e
      procedure :: upper_e => conformal_sgbc_upper_e
      procedure :: destroy => conformal_sgbc_destroy
   end type conformal_sgbc_state_t

contains

   subroutine conformal_sgbc_init(this, profile, dt, lower_distance, upper_distance, maximum_frequency, resolution, depth, reverse_layers)
      class(conformal_sgbc_state_t), intent(inout) :: this
      type(SGBCMaterialProfile_t), intent(in) :: profile
      real(kind=RKIND_tiempo), intent(in) :: dt
      real(kind=RKIND), intent(in) :: lower_distance, upper_distance, maximum_frequency, resolution
      integer, intent(in) :: depth
      logical, intent(in) :: reverse_layers
      integer, allocatable :: subdivisions(:)
      real(kind=RKIND), allocatable :: e_ds(:), h_ds(:), e_xi(:), h_xi(:), e_sigma(:), h_sigma(:)
      real(kind=RKIND), allocatable :: cell_eps(:), cell_sigma(:)
      integer :: layer, source_layer, first, last, index
      real(kind=RKIND) :: fine_step

      call this%destroy()
      if (profile%num_layers <= 0) then
         call WarnErrReport('Cannot initialize a conformal SGBC with an empty material profile.', .true.)
         return
      end if
      if (depth == 0 .and. profile%num_layers > 1) then
         call WarnErrReport('Conformal SGBC depth 0 only supports a single material layer.', .true.)
         return
      end if

      allocate(subdivisions(profile%num_layers))
      do layer = 1, profile%num_layers
         subdivisions(layer) = layer_subdivisions(profile, layer, maximum_frequency, resolution, depth)
      end do
      this%n = sum(subdivisions)+1
      allocate(this%e(1:this%n), this%e_old(1:this%n), this%h(0:this%n), this%rhs(1:this%n))
      allocate(this%e_ga(1:this%n), this%e_gb(1:this%n), this%h_ga(0:this%n), this%h_gb(0:this%n))
      allocate(this%lower_diag(1:this%n), this%diagonal(1:this%n), this%upper_diag(1:this%n))
      allocate(this%rhs_coeff(3,1:this%n))
      allocate(this%upper_prime(1:this%n), this%rhs_prime(1:this%n))
      allocate(e_ds(1:this%n), h_ds(0:this%n), e_xi(1:this%n), h_xi(0:this%n), &
               e_sigma(1:this%n), h_sigma(0:this%n))
      allocate(cell_eps(0:this%n), cell_sigma(0:this%n))

      h_ds = 0.0_RKIND
      cell_eps = EPSILON_VACUUM
      cell_sigma = 0.0_RKIND
      h_xi = MU_VACUUM
      h_sigma = 0.0_RKIND
      last = 0
      do layer = 1, profile%num_layers
         if (reverse_layers) then
            source_layer = profile%num_layers-layer+1
         else
            source_layer = layer
         end if
         first = last+1
         last = first+subdivisions(source_layer)-1
         fine_step = profile%thickness(source_layer)/real(subdivisions(source_layer), RKIND)
         h_ds(first:last) = fine_step
         cell_eps(first:last) = profile%eps(source_layer)
         cell_sigma(first:last) = profile%sigma(source_layer)
         h_xi(first:last) = profile%mu(source_layer)
         h_sigma(first:last) = profile%sigmam(source_layer)
      end do
      h_ds(0) = lower_distance
      h_ds(this%n) = upper_distance
      h_xi(0) = MU_VACUUM
      h_xi(this%n) = MU_VACUUM
      h_sigma(0) = 0.0_RKIND
      h_sigma(this%n) = 0.0_RKIND
      do index = 1, this%n
         ! Each E node spans half of the magnetic cell on either side.
         ! The common factor one half cancels in this volume average.  This
         ! handles vacuum/material boundaries and dissimilar layer interfaces
         ! with the same expression.
         e_xi(index) = (cell_eps(index-1)*h_ds(index-1)+cell_eps(index)*h_ds(index))/ &
            (h_ds(index-1)+h_ds(index))
         e_sigma(index) = (cell_sigma(index-1)*h_ds(index-1)+cell_sigma(index)*h_ds(index))/ &
            (h_ds(index-1)+h_ds(index))
         e_ds(index) = 0.5_RKIND*(h_ds(index-1)+h_ds(index))
         this%e_ga(index) = advance_a(real(dt,RKIND), e_sigma(index), e_xi(index))
         this%e_gb(index) = advance_b(real(dt,RKIND), e_sigma(index), e_xi(index), e_ds(index))
      end do
      this%h_ga(0) = 1.0_RKIND
      this%h_ga(this%n) = 1.0_RKIND
      this%h_gb(0) = 0.0_RKIND
      this%h_gb(this%n) = 0.0_RKIND
      do index = 1, this%n-1
         this%h_ga(index) = advance_a(real(dt,RKIND), h_sigma(index), h_xi(index))
         this%h_gb(index) = advance_b(real(dt,RKIND), h_sigma(index), h_xi(index), h_ds(index))
      end do
      call initialize_cn_coefficients(this)
      this%e = 0.0_RKIND
      this%e_old = 0.0_RKIND
      this%h = 0.0_RKIND
      this%rhs = 0.0_RKIND
      this%initialized = .true.
   end subroutine conformal_sgbc_init

   integer function layer_subdivisions(profile, layer, maximum_frequency, resolution, depth) result(count)
      type(SGBCMaterialProfile_t), intent(in) :: profile
      integer, intent(in) :: layer, depth
      real(kind=RKIND), intent(in) :: maximum_frequency, resolution
      real(kind=RKIND) :: skin_depth, omega, argument

      if (depth == 0) then
         count = 1
      else if (depth > 0) then
         count = depth
      else if (maximum_frequency <= 0.0_RKIND) then
         count = 2
      else
         omega = 2.0_RKIND*acos(-1.0_RKIND)*maximum_frequency
         argument = sqrt((omega*profile%eps(layer))**2 + profile%sigma(layer)**2)- &
                    omega*profile%eps(layer)
         if (argument <= tiny(1.0_RKIND)) then
            skin_depth = huge(1.0_RKIND)
         else
            skin_depth = 1.0_RKIND/sqrt(0.5_RKIND*omega*profile%mu(layer)*argument)
         end if
         count = 1+int(resolution*profile%thickness(layer)/skin_depth)
         count = max(2, count)
      end if
   end function layer_subdivisions

   real(kind=RKIND) function advance_a(dt, sigma, xi) result(value)
      real(kind=RKIND), intent(in) :: dt, sigma, xi
      value = (2.0_RKIND*xi-dt*sigma)/(2.0_RKIND*xi+dt*sigma)
      if (value < 0.0_RKIND) value = exp(-sigma*dt/xi)
   end function advance_a

   real(kind=RKIND) function advance_b(dt, sigma, xi, step) result(value)
      real(kind=RKIND), intent(in) :: dt, sigma, xi, step
      real(kind=RKIND) :: trapezoidal_a
      trapezoidal_a = (2.0_RKIND*xi-dt*sigma)/(2.0_RKIND*xi+dt*sigma)
      if (trapezoidal_a < 0.0_RKIND) then
         value = (1.0_RKIND-exp(-sigma*dt/xi))/(sigma*step)
      else
         value = 2.0_RKIND*dt/(step*(2.0_RKIND*xi+dt*sigma))
      end if
   end function advance_b

   subroutine initialize_cn_coefficients(this)
      class(conformal_sgbc_state_t), intent(inout) :: this
      integer :: index
      this%lower_diag = 0.0_RKIND
      this%upper_diag = 0.0_RKIND
      do index = 2, this%n
         this%lower_diag(index) = -0.25_RKIND*this%e_gb(index)*this%h_gb(index-1)
      end do
      do index = 1, this%n-1
         this%upper_diag(index) = -0.25_RKIND*this%e_gb(index)*this%h_gb(index)
      end do
      do index = 1, this%n
         this%diagonal(index) = 1.0_RKIND+0.25_RKIND*this%e_gb(index)* &
            (this%h_gb(index-1)+this%h_gb(index))
         this%rhs_coeff(1,index) = 0.5_RKIND*this%e_gb(index)*(1.0_RKIND+this%h_ga(index-1))
         this%rhs_coeff(2,index) = 0.5_RKIND*this%e_gb(index)*(1.0_RKIND+this%h_ga(index))
         this%rhs_coeff(3,index) = this%e_ga(index)-0.25_RKIND*this%e_gb(index)* &
            (this%h_gb(index-1)+this%h_gb(index))
      end do
      this%diagonal(1) = 1.0_RKIND+0.25_RKIND*this%e_gb(1)*this%h_gb(1)
      this%rhs_coeff(1,1) = this%e_gb(1)
      this%rhs_coeff(3,1) = this%e_ga(1)-0.25_RKIND*this%e_gb(1)*this%h_gb(1)
      this%diagonal(this%n) = 1.0_RKIND+0.25_RKIND*this%e_gb(this%n)*this%h_gb(this%n-1)
      this%rhs_coeff(2,this%n) = this%e_gb(this%n)
      this%rhs_coeff(3,this%n) = this%e_ga(this%n)-0.25_RKIND*this%e_gb(this%n)*this%h_gb(this%n-1)
   end subroutine initialize_cn_coefficients

   subroutine conformal_sgbc_advance(this, lower_h, upper_h)
      class(conformal_sgbc_state_t), intent(inout) :: this
      real(kind=RKIND), intent(in) :: lower_h, upper_h
      integer :: index
      if (.not. this%initialized) return
      this%h(0) = lower_h
      this%h(this%n) = upper_h
      do index = 1, this%n
         this%rhs(index) = this%rhs_coeff(1,index)*this%h(index-1)- &
            this%rhs_coeff(2,index)*this%h(index)+this%rhs_coeff(3,index)*this%e(index)
         if (index > 1) this%rhs(index) = this%rhs(index)-this%lower_diag(index)*this%e(index-1)
         if (index < this%n) this%rhs(index) = this%rhs(index)-this%upper_diag(index)*this%e(index+1)
      end do
      this%e_old = this%e
      call solve_tridiagonal(this%lower_diag, this%diagonal, this%upper_diag, this%rhs, this%e, &
         this%upper_prime, this%rhs_prime)
      do index = 1, this%n-1
         this%h(index) = this%h_ga(index)*this%h(index)+0.5_RKIND*this%h_gb(index)* &
            (this%e(index)-this%e(index+1)+this%e_old(index)-this%e_old(index+1))
      end do
   end subroutine conformal_sgbc_advance

   subroutine solve_tridiagonal(lower, diagonal, upper, rhs, solution, upper_prime, rhs_prime)
      real(kind=RKIND), intent(in) :: lower(:), diagonal(:), upper(:), rhs(:)
      real(kind=RKIND), intent(out) :: solution(:)
      real(kind=RKIND), intent(inout) :: upper_prime(:), rhs_prime(:)
      real(kind=RKIND) :: pivot
      integer :: index, n
      n = size(diagonal)
      upper_prime(1) = upper(1)/diagonal(1)
      rhs_prime(1) = rhs(1)/diagonal(1)
      do index = 2, n
         pivot = diagonal(index)-lower(index)*upper_prime(index-1)
         upper_prime(index) = upper(index)/pivot
         rhs_prime(index) = (rhs(index)-lower(index)*rhs_prime(index-1))/pivot
      end do
      solution(n) = rhs_prime(n)
      do index = n-1, 1, -1
         solution(index) = rhs_prime(index)-upper_prime(index)*solution(index+1)
      end do
   end subroutine solve_tridiagonal

   real(kind=RKIND) function conformal_sgbc_lower_e(this) result(value)
      class(conformal_sgbc_state_t), intent(in) :: this
      value = 0.0_RKIND
      if (this%initialized) value = this%e(1)
   end function conformal_sgbc_lower_e

   real(kind=RKIND) function conformal_sgbc_upper_e(this) result(value)
      class(conformal_sgbc_state_t), intent(in) :: this
      value = 0.0_RKIND
      if (this%initialized) value = this%e(this%n)
   end function conformal_sgbc_upper_e

   subroutine conformal_sgbc_destroy(this)
      class(conformal_sgbc_state_t), intent(inout) :: this
      if (allocated(this%e)) deallocate(this%e)
      if (allocated(this%e_old)) deallocate(this%e_old)
      if (allocated(this%h)) deallocate(this%h)
      if (allocated(this%rhs)) deallocate(this%rhs)
      if (allocated(this%e_ga)) deallocate(this%e_ga)
      if (allocated(this%e_gb)) deallocate(this%e_gb)
      if (allocated(this%h_ga)) deallocate(this%h_ga)
      if (allocated(this%h_gb)) deallocate(this%h_gb)
      if (allocated(this%lower_diag)) deallocate(this%lower_diag)
      if (allocated(this%diagonal)) deallocate(this%diagonal)
      if (allocated(this%upper_diag)) deallocate(this%upper_diag)
      if (allocated(this%rhs_coeff)) deallocate(this%rhs_coeff)
      if (allocated(this%upper_prime)) deallocate(this%upper_prime)
      if (allocated(this%rhs_prime)) deallocate(this%rhs_prime)
      this%n = 0
      this%initialized = .false.
   end subroutine conformal_sgbc_destroy

end module conformal_sgbc_m
