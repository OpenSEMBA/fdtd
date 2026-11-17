module lorentz_testingTools_m
   use FDETYPES_m
   use EDispersives_m
   implicit none
   private
   public :: test_lorentz_ade_single_node

   ! Single-pole-pair medium used by the checks below. The values correspond to
   ! eps_r(w) = 2 + omegaP^2/(omega1^2 - w^2 + j*gamma*w) with omegaP = 3c,
   ! omega1 = 2c and gamma = 0.3c, already converted to the first-order form
   ! consumed by InitEDispersives.
   complex(kind=CKIND), parameter :: POLE1 = (4.49688687e7_CKIND, 5.97896205325447083e8_CKIND)
   complex(kind=CKIND), parameter :: RESIDUE1 = (0.0_CKIND, 5.98931083300393435e-3_CKIND)

   real(kind=RKIND), parameter :: TOL = 1.0e-5_RKIND

contains

   ! Unit test of the electric dispersive subroutines on a one-node mesh:
   ! InitEDispersives must produce the documented ADE coefficients and
   ! AdvanceEDispersiveE must follow the documented polarization recurrence.
   integer function test_lorentz_ade_single_node() bind(C) result(err)
      type(SGGFDTDINFO_t) :: sgg
      type(media_matrices_t) :: media
      real(kind=RKIND), target :: ex(1, 1, 1), ey(1, 1, 1), ez(1, 1, 1)
      real(kind=RKIND) :: g1(0:1), g2(0:1)
      real(kind=RKIND) :: dt, eps, sigma
      complex(kind=CKIND) :: kappa(2), beta(2), g3(2), current(2)
      real(kind=RKIND) :: sumBeta, g1Expected, g2Expected
      real(kind=RKIND) :: previous, expected
      logical :: thereAre
      integer :: step

      err = 0
      dt = 1.0e-11_RKIND
      eps = 2.0_RKIND*EPSILON_VACUUM
      sigma = 0.0_RKIND

      call buildSingleNodeMedium(sgg, media, eps, sigma, dt)
      ex = 0.0_RKIND
      ey = 0.0_RKIND
      ez = 0.0_RKIND
      g1 = 0.0_RKIND
      g2 = 0.0_RKIND

      call InitEDispersives(sgg, media, thereAre, .false., g1, g2, ex, ey, ez)
      if (.not. thereAre) call testFails(err, 'InitEDispersives should report electric dispersive nodes')

      ! Documented coefficients: Kappa = (1 + a*dt/2)/(1 - a*dt/2),
      ! Beta = C*dt/(1 - a*dt/2), G2 = 2*dt/(2*eps + sum(Re(Beta)) + sigma*dt)
      ! and G1 = (2*eps + sum(Re(Beta)) - sigma*dt)/(2*eps + sum(Re(Beta)) + sigma*dt).
      kappa(1) = (1.0_RKIND + POLE1*dt/2.0_RKIND)/(1.0_RKIND - POLE1*dt/2.0_RKIND)
      kappa(2) = conjg(kappa(1))
      beta(1) = RESIDUE1*dt/(1.0_RKIND - POLE1*dt/2.0_RKIND)
      beta(2) = conjg(beta(1))
      sumBeta = real(beta(1)) + real(beta(2))
      g2Expected = 2.0_RKIND*dt/(2.0_RKIND*eps + sumBeta + sigma*dt)
      g1Expected = (2.0_RKIND*eps + sumBeta - sigma*dt)/(2.0_RKIND*eps + sumBeta + sigma*dt)
      g3(1) = g2Expected/2.0_RKIND*(1.0_RKIND + kappa(1))
      g3(2) = g2Expected/2.0_RKIND*(1.0_RKIND + kappa(2))

      if (abs(g1(1) - g1Expected) > TOL*abs(g1Expected)) call testFails(err, 'InitEDispersives G1')
      if (abs(g2(1) - g2Expected) > TOL*abs(g2Expected)) call testFails(err, 'InitEDispersives G2')

      ! Advance the single node with a prescribed E and follow the recurrence:
      ! E <- E - Re(sum_k G3_k*J_k), J_k <- Kappa_k*J_k + Beta_k*(E - Eprevious)/dt.
      current = (0.0_RKIND, 0.0_RKIND)
      previous = 0.0_RKIND
      do step = 1, 3
         ex(1, 1, 1) = real(step, RKIND)
         call AdvanceEDispersiveE(sgg)
         expected = real(step, RKIND) - real(g3(1)*current(1) + g3(2)*current(2))
         if (abs(ex(1, 1, 1) - expected) > TOL*max(1.0_RKIND, abs(expected))) &
            call testFails(err, 'AdvanceEDispersiveE should subtract the polarization current')
         current(1) = kappa(1)*current(1) + beta(1)*(ex(1, 1, 1) - previous)/dt
         current(2) = kappa(2)*current(2) + beta(2)*(ex(1, 1, 1) - previous)/dt
         previous = ex(1, 1, 1)
      end do

      call DestroyEDispersives(sgg)
      deallocate(sgg%Med)
      deallocate(media%sggMiEx, media%sggMiEy, media%sggMiEz)
   end function

   subroutine buildSingleNodeMedium(sgg, media, eps, sigma, dt)
      type(SGGFDTDINFO_t), intent(inout) :: sgg
      type(media_matrices_t), intent(inout) :: media
      real(kind=RKIND), intent(in) :: eps, sigma, dt
      integer :: field

      allocate(sgg%Med(0:1))
      sgg%NumMedia = 1
      sgg%dt = dt
      sgg%Med(1)%Is%PML = .false.
      sgg%Med(1)%Is%EDispersive = .true.
      sgg%Med(1)%Is%EDispersiveAnis = .false.
      allocate(sgg%Med(1)%EDispersive(1))
      sgg%Med(1)%EDispersive(1)%numpolres11 = 2
      sgg%Med(1)%EDispersive(1)%eps11 = eps
      sgg%Med(1)%EDispersive(1)%Sigma11 = sigma
      allocate(sgg%Med(1)%EDispersive(1)%a11(2), sgg%Med(1)%EDispersive(1)%C11(2))
      sgg%Med(1)%EDispersive(1)%a11(1) = POLE1
      sgg%Med(1)%EDispersive(1)%a11(2) = conjg(POLE1)
      sgg%Med(1)%EDispersive(1)%C11(1) = RESIDUE1
      sgg%Med(1)%EDispersive(1)%C11(2) = conjg(RESIDUE1)

      do field = IEX, IHZ
         sgg%Alloc(field)%XI = 1
         sgg%Alloc(field)%XE = 1
         sgg%Alloc(field)%YI = 1
         sgg%Alloc(field)%YE = 1
         sgg%Alloc(field)%ZI = 1
         sgg%Alloc(field)%ZE = 1
         sgg%Sweep(field)%XI = 1
         sgg%Sweep(field)%XE = 1
         sgg%Sweep(field)%YI = 1
         sgg%Sweep(field)%YE = 1
         sgg%Sweep(field)%ZI = 1
         sgg%Sweep(field)%ZE = 1
      end do

      allocate(media%sggMiEx(1, 1, 1), media%sggMiEy(1, 1, 1), media%sggMiEz(1, 1, 1))
      media%sggMiEx = 0
      media%sggMiEy = 0
      media%sggMiEz = 0
      media%sggMiEx(1, 1, 1) = 1
   end subroutine

   subroutine testFails(err, message)
      integer, intent(inout) :: err
      character(len=*), intent(in) :: message
      err = err + 1
      write (*, '(a)') 'Lorentz ADE test: '//trim(message)
   end subroutine

end module
