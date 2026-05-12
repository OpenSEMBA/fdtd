!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  GPU NF2FF MODULE - CUDA Fortran (CUF) accelerated near-to-far-field
!  Infrastructure for GPU-accelerated far-field pattern computation.
!  NOTE: GPU kernel not yet implemented — falls through to CPU path.
!  Implementation pending: full face-specific indexing from farfield.F90.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module gpu_nf2ff_m

   use FDETYPES_m
   use Report_m
   use cudafor
   use gpu_core_m

   implicit none

contains

   !--------------------------------------------------------------------------------
   ! Device function: compute average of two complex values (geometric or arithmetic)
   !--------------------------------------------------------------------------------
   attributes(device) function AverageNF2FF(pasadas, z1, z2) result(z)
      complex(kind=rkind), intent(in) :: z1, z2
      integer, intent(in) :: pasadas
      complex(kind=rkind) :: z
      real(kind=rkind) :: phi1, phi2, nphi1, nphi2

      z = (0.0_rkind, 0.0_rkind)
      if (pasadas == 2) then
         phi1 = atan2(imag(z1), real(z1))
         phi2 = atan2(imag(z2), real(z2))
         nphi1 = phi1
         nphi2 = phi2
         if ((phi1 < -pi/2.0_rkind) .and. (phi2 > pi/2.0_rkind)) nphi1 = phi1 - 2.0_rkind * pi
         if ((phi2 < -pi/2.0_rkind) .and. (phi1 > pi/2.0_rkind)) nphi2 = phi2 - 2.0_rkind * pi
         z = sqrt(abs(z1*z2)) * exp((0.0_rkind, 1.0_rkind) * (nphi1 + nphi2) / 2.0_rkind)
      elseif (pasadas == 1) then
         z = (z1 + z2) / 2.0_rkind
      endif
   end function AverageNF2FF

   !--------------------------------------------------------------------------------
   ! Public interface: GPU far-field flush
   ! NOTE: GPU kernel not yet implemented. Falls through to CPU.
   ! The gpu_state_t%nf2ff_initialized flag is set during init,
   ! but the actual GPU computation is deferred until the full
   ! face-specific indexing from farfield.F90 is implemented.
   !--------------------------------------------------------------------------------
   subroutine gpu_flush_nf2ff(this, Etheta_h, Ephi_h, RCS_h)
      class(gpu_state_t), intent(inout) :: this
      real(kind=rkind), dimension(:), intent(out) :: Etheta_h, Ephi_h, RCS_h

      ! GPU kernel not yet implemented — CPU path handles this.
      ! The nf2ff_initialized flag is set during init, but the actual
      ! GPU flush is a no-op until full implementation.
      ! TODO: Implement gpu_flush_nf2ff_kernel when NVHPC argument
      ! limit is addressed (use derived type or split into multiple kernels).

      ! Zero output buffers (CPU fallback will overwrite)
      Etheta_h = 0.0_rkind
      Ephi_h = 0.0_rkind
      RCS_h = 0.0_rkind

   end subroutine gpu_flush_nf2ff

end module gpu_nf2ff_m
