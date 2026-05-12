!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  GPU NF2FF MODULE - CUDA Fortran (CUF) accelerated near-to-far-field
!  Computes far-field radiation patterns from Huygens box surface data.
!  Two phases: DFT accumulation (UpdateFarField) and far-field pattern (FlushFarfield).
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module gpu_nf2ff_m

   use FDETYPES_m
   use Report_m
   use cudafor
   use gpu_core_m

   implicit none

contains

   !--------------------------------------------------------------------------------
   ! GPU kernel: far-field pattern computation
   ! Each thread block handles one (freq, theta, phi) combination
   !--------------------------------------------------------------------------------
   attributes(global) subroutine gpu_flush_nf2ff_kernel(
        // DFT buffer arrays (read-only, 18 complex 3D arrays)
        ExIz_d, ExDe_d, ExAb_d, ExAr_d, EyFr_d, EyTr_d, EyAb_d, EyAr_d,
        EzIz_d, EzDe_d, EzFr_d, EzTr_d,
        HxIz_d, HxDe_d, HxAb_d, HxAr_d, HyFr_d, HyTr_d, HyAb_d, HyAr_d,
        HzIz_d, HzDe_d, HzFr_d, HzTr_d,
        HxIz2_d, HxDe2_d, HxAb2_d, HxAr2_d, HyFr2_d, HyTr2_d, HyAb2_d, HyAr2_d,
        HzIz2_d, HzDe2_d, HzFr2_d, HzTr2_d,

        // Geometry arrays (read-only)
        phys_x_Mx_d, phys_y_Mx_d, phys_z_Mx_d,
        phys_x_My_d, phys_y_My_d, phys_z_My_d,
        phys_x_Mz_d, phys_y_Mz_d, phys_z_Mz_d,
        phys_x_Jx_d, phys_y_Jx_d, phys_z_Jx_d,
        phys_x_Jy_d, phys_y_Jy_d, phys_z_Jy_d,
        phys_x_Jz_d, phys_y_Jz_d, phys_z_Jz_d,

        // Cell dimensions (read-only)
        dyh_d, dze_d, dye_d, dzh_d,

        // DFT phase factors (read-only)
        expIwdt_d, auxExp_E_d, auxExp_H_d,

        // Output buffers
        Etheta_d, Ephi_d, RCS_d,

        // Configuration
        numFreqs, Ntheta, Nphi,
        thetaStart, thetaStop, thetaStep, phiStart, phiStop, phiStep,
        freqStep, initialFreq, comun, cluz, z0,
        XDobleAncho, YDobleAncho, ZDobleAncho,
        sym_flags,
        // Face enable flags
        farfieldTr, farfieldFr, farfieldIz, farfieldDe, farfieldAb, farfieldAr,
        // Symmetry flags
        symTrPEC_Front, symTrPEC_Left, symTrPEC_Right, symTrPEC_Up, symTrPEC_Down,
        symTrPMC_Front, symTrPMC_Left, symTrPMC_Right, symTrPMC_Up, symTrPMC_Down,
        symFrPEC_Back, symFrPEC_Left, symFrPEC_Right, symFrPEC_Up, symFrPEC_Down,
        symFrPMC_Back, symFrPMC_Left, symFrPMC_Right, symFrPMC_Up, symFrPMC_Down,
        symIzPEC_Back, symIzPEC_Front, symIzPEC_Right, symIzPEC_Up, symIzPEC_Down,
        symIzPMC_Back, symIzPMC_Front, symIzPMC_Right, symIzPMC_Up, symIzPMC_Down,
        symDePEC_Back, symDePEC_Front, symDePEC_Left, symDePEC_Up, symDePEC_Down,
        symDePMC_Back, symDePMC_Front, symDePMC_Left, symDePMC_Up, symDePMC_Down,
        symArPEC_Back, symArPEC_Front, symArPEC_Left, symArPEC_Right, symArPEC_Down,
        symArPMC_Back, symArPMC_Front, symArPMC_Left, symArPMC_Right, symArPMC_Down,
        symAbPEC_Back, symAbPEC_Front, symAbPEC_Left, symAbPEC_Right, symAbPEC_Up,
        symAbPMC_Back, symAbPMC_Front, symAbPMC_Left, symAbPMC_Right, symAbPMC_Up
   )

      complex(kind=rkind), dimension(:,:,:), device, intent(in) :: ExIz, ExDe, ExAb, ExAr, EyFr, EyTr, EyAb, EyAr, EzIz, EzDe, EzFr, EzTr
      complex(kind=rkind), dimension(:,:,:), device, intent(in) :: HxIz, HxDe, HxAb, HxAr, HyFr, HyTr, HyAb, HyAr, HzIz, HzDe, HzFr, HzTr
      complex(kind=rkind), dimension(:,:,:), device, intent(in) :: HxIz2, HxDe2, HxAb2, HxAr2, HyFr2, HyTr2, HyAb2, HyAr2, HzIz2, HzDe2, HzFr2, HzTr2
      real(kind=rkind), dimension(:), device, intent(in) :: phys_x_Mx_d, phys_y_Mx_d, phys_z_Mx_d
      real(kind=rkind), dimension(:), device, intent(in) :: phys_x_My_d, phys_y_My_d, phys_z_My_d
      real(kind=rkind), dimension(:), device, intent(in) :: phys_x_Mz_d, phys_y_Mz_d, phys_z_Mz_d
      real(kind=rkind), dimension(:), device, intent(in) :: phys_x_Jx_d, phys_y_Jx_d, phys_z_Jx_d
      real(kind=rkind), dimension(:), device, intent(in) :: phys_x_Jy_d, phys_y_Jy_d, phys_z_Jy_d
      real(kind=rkind), dimension(:), device, intent(in) :: phys_x_Jz_d, phys_y_Jz_d, phys_z_Jz_d
      real(kind=rkind), dimension(:), device, intent(in) :: dyh_d, dze_d, dye_d, dzh_d
      complex(kind=rkind), dimension(:), device, intent(in) :: expIwdt_d, auxExp_E_d, auxExp_H_d
      real(kind=rkind), dimension(:), device, intent(out) :: Etheta_d, Ephi_d, RCS_d
      integer, value :: numFreqs, Ntheta, Nphi
      integer, value :: thetaStart, thetaStop, thetaStep, phiStart, phiStop, phiStep
      real(kind=rkind), value :: freqStep, initialFreq, comun, cluz, z0
      real(kind=rkind), value :: XDobleAncho, YDobleAncho, ZDobleAncho
      integer, value :: sym_flags
      logical, value :: farfieldTr, farfieldFr, farfieldIz, farfieldDe, farfieldAb, farfieldAr
      logical, value :: symTrPEC_Front, symTrPEC_Left, symTrPEC_Right, symTrPEC_Up, symTrPEC_Down
      logical, value :: symTrPMC_Front, symTrPMC_Left, symTrPMC_Right, symTrPMC_Up, symTrPMC_Down
      logical, value :: symFrPEC_Back, symFrPEC_Left, symFrPEC_Right, symFrPEC_Up, symFrPEC_Down
      logical, value :: symFrPMC_Back, symFrPMC_Left, symFrPMC_Right, symFrPMC_Up, symFrPMC_Down
      logical, value :: symIzPEC_Back, symIzPEC_Front, symIzPEC_Right, symIzPEC_Up, symIzPEC_Down
      logical, value :: symIzPMC_Back, symIzPMC_Front, symIzPMC_Right, symIzPMC_Up, symIzPMC_Down
      logical, value :: symDePEC_Back, symDePEC_Front, symDePEC_Left, symDePEC_Up, symDePEC_Down
      logical, value :: symDePMC_Back, symDePMC_Front, symDePMC_Left, symDePMC_Up, symDePMC_Down
      logical, value :: symArPEC_Back, symArPEC_Front, symArPEC_Left, symArPEC_Right, symArPEC_Down
      logical, value :: symArPMC_Back, symArPMC_Front, symArPMC_Left, symArPMC_Right, symArPMC_Down
      logical, value :: symAbPEC_Back, symAbPEC_Front, symAbPEC_Left, symAbPEC_Right, symAbPEC_Up
      logical, value :: symAbPMC_Back, symAbPMC_Front, symAbPMC_Left, symAbPMC_Right, symAbPMC_Up

      integer :: idx, freqIdx, thetaIdx, phiIdx, ii, pasadas, donde, cellIdx
      real(kind=rkind) :: theta, phi, freq, sintheta, costheta, sinphi, cosphi
      real(kind=rkind) :: sintheta_sinphi, sintheta_cosphi, costheta_cosphi, costheta_sinphi
      complex(kind=rkind) :: L_theta, L_phi, N_theta, N_phi
      complex(kind=rkind) :: L_theta_final, L_phi_final, N_theta_final, N_phi_final
      complex(kind=rkind) :: Mx, My, Mz, Jx, Jy, Jz
      complex(kind=rkind) :: comunMx, comunMy, comunMz, comunJx, comunJy, comunJz
      complex(kind=rkind) :: Etheta, Ephi
      real(kind=rkind) :: rcs
      real(kind=rkind) :: normal

      ! Decode thread index to (freq, theta, phi)
      idx = (blockidx%x - 1) * blockdim%x + threadidx%x
      if (idx > numFreqs * Ntheta * Nphi) return

      freqIdx = (idx - 1) / (Ntheta * Nphi) + 1
      phiIdx = mod((idx - 1) / Ntheta, Nphi) + 1
      thetaIdx = mod(idx - 1, Ntheta) + 1

      ! Compute frequency
      freq = initialFreq + (freqIdx - 1) * freqStep
      if (freqIdx > 1) freq = 10.0_rkind ** freq  ! logarithmic spacing

      ! Compute theta and phi
      theta = thetaStart + (thetaIdx - 1) * thetaStep
      phi = phiStart + (phiIdx - 1) * phiStep

      ! Precompute direction cosines (common for all cells at this direction)
      sintheta = Sin(theta)
      costheta = Cos(theta)
      sinphi = Sin(phi)
      cosphi = Cos(phi)
      sintheta_sinphi = sintheta * sinphi
      sintheta_cosphi = sintheta * cosphi
      costheta_cosphi = costheta * cosphi
      costheta_sinphi = costheta * sinphi

      ! Initialize accumulators
      L_theta_final = (0.0_rkind, 0.0_rkind)
      L_phi_final = (0.0_rkind, 0.0_rkind)
      N_theta_final = (0.0_rkind, 0.0_rkind)
      N_phi_final = (0.0_rkind, 0.0_rkind)

      ! Two passes: geometric (pasadas=2) then arithmetic (pasadas=1)
      do pasadas = 2, 1, -1

         ! Three face-pairs: Back/Front, Left/Right, Bottom/Top
         ! Pair 1: Back (Tr) / Front (Fr)
         if (farfieldTr) then
            normal = -1.0_rkind
            ! Process Back face cells
            ! Each thread handles one cell on the face
            ! Cell index mapping: cellIdx = (thetaIdx-1)*Nphi + phiIdx (unique per direction)
            ! We need to iterate over all cells on the face
            ! For simplicity, use a single thread per cell approach
            ! In practice, we'd launch more threads and have each thread handle one cell
            call process_face_cells(ExIz_d, ExDe_d, ExAb_d, ExAr_d, EyFr_d, EyTr_d, EyAb_d, EyAr_d, &
                                    EzIz_d, EzDe_d, EzFr_d, EzTr_d, &
                                    HxIz_d, HxDe_d, HxAb_d, HxAr_d, HyFr_d, HyTr_d, HyAb_d, HyAr_d, &
                                    HzIz_d, HzDe_d, HzFr_d, HzTr_d, &
                                    HxIz2_d, HxDe2_d, HxAb2_d, HxAr2_d, HyFr2_d, HyTr2_d, HyAb2_d, HyAr2_d, &
                                    HzIz2_d, HzDe2_d, HzFr2_d, HzTr2_d, &
                                    phys_x_Mx_d, phys_y_Mx_d, phys_z_Mx_d, &
                                    phys_x_My_d, phys_y_My_d, phys_z_My_d, &
                                    phys_x_Mz_d, phys_y_Mz_d, phys_z_Mz_d, &
                                    phys_x_Jx_d, phys_y_Jx_d, phys_z_Jx_d, &
                                    phys_x_Jy_d, phys_y_Jy_d, phys_z_Jy_d, &
                                    phys_x_Jz_d, phys_y_Jz_d, phys_z_Jz_d, &
                                    dyh_d, dze_d, dye_d, dzh_d, &
                                    expIwdt_d, auxExp_E_d, auxExp_H_d, &
                                    comun, sintheta_cosphi, sintheta_sinphi, costheta, &
                                    costheta_cosphi, costheta_sinphi, sintheta, sinphi, cosphi, &
                                    L_theta, L_phi, N_theta, N_phi, &
                                    pasadas, normal, freqIdx, &
                                    symTrPEC_Front, symTrPEC_Left, symTrPEC_Right, symTrPEC_Up, symTrPEC_Down, &
                                    symTrPMC_Front, symTrPMC_Left, symTrPMC_Right, symTrPMC_Up, symTrPMC_Down, &
                                    1, freq, cluz, z0)
            L_theta_final = L_theta_final + L_theta
            L_phi_final = L_phi_final + L_phi
            N_theta_final = N_theta_final + N_theta
            N_phi_final = N_phi_final + N_phi
         endif

         if (farfieldFr) then
            normal = +1.0_rkind
            call process_face_cells(ExIz_d, ExDe_d, ExAb_d, ExAr_d, EyFr_d, EyTr_d, EyAb_d, EyAr_d, &
                                    EzIz_d, EzDe_d, EzFr_d, EzTr_d, &
                                    HxIz_d, HxDe_d, HxAb_d, HxAr_d, HyFr_d, HyTr_d, HyAb_d, HyAr_d, &
                                    HzIz_d, HzDe_d, HzFr_d, HzTr_d, &
                                    HxIz2_d, HxDe2_d, HxAb2_d, HxAr2_d, HyFr2_d, HyTr2_d, HyAb2_d, HyAr2_d, &
                                    HzIz2_d, HzDe2_d, HzFr2_d, HzTr2_d, &
                                    phys_x_Mx_d, phys_y_Mx_d, phys_z_Mx_d, &
                                    phys_x_My_d, phys_y_My_d, phys_z_My_d, &
                                    phys_x_Mz_d, phys_y_Mz_d, phys_z_Mz_d, &
                                    phys_x_Jx_d, phys_y_Jx_d, phys_z_Jx_d, &
                                    phys_x_Jy_d, phys_y_Jy_d, phys_z_Jy_d, &
                                    phys_x_Jz_d, phys_y_Jz_d, phys_z_Jz_d, &
                                    dyh_d, dze_d, dye_d, dzh_d, &
                                    expIwdt_d, auxExp_E_d, auxExp_H_d, &
                                    comun, sintheta_cosphi, sintheta_sinphi, costheta, &
                                    costheta_cosphi, costheta_sinphi, sintheta, sinphi, cosphi, &
                                    L_theta, L_phi, N_theta, N_phi, &
                                    pasadas, normal, freqIdx, &
                                    symFrPEC_Back, symFrPEC_Left, symFrPEC_Right, symFrPEC_Up, symFrPEC_Down, &
                                    symFrPMC_Back, symFrPMC_Left, symFrPMC_Right, symFrPMC_Up, symFrPMC_Down, &
                                    2, freq, cluz, z0)
            L_theta_final = L_theta_final + L_theta
            L_phi_final = L_phi_final + L_phi
            N_theta_final = N_theta_final + N_theta
            N_phi_final = N_phi_final + N_phi
         endif

         ! Pair 2: Left (Iz) / Right (De)
         ! ... (similar pattern)

         ! Pair 3: Bottom (Ab) / Top (Ar)
         ! ... (similar pattern)

      end do ! pasadas

      ! Compute far-field patterns
      Etheta = -(j*freq/(2.0_rkind*cluz)) * (L_phi_final + z0 * N_theta_final)
      Ephi = (j*freq/(2.0_rkind*cluz)) * (L_theta_final - z0 * N_phi_final)

      ! Output indices
      ii = (freqIdx - 1) * Ntheta * Nphi + (thetaIdx - 1) * Nphi + phiIdx
      Etheta_d(ii) = abs(Etheta)
      Ephi_d(ii) = abs(Ephi)
      RCS_d(ii) = 4.0_rkind * pi * abs(Etheta)**2 / (freq/cluz)**2

   end subroutine gpu_flush_nf2ff_kernel

   !--------------------------------------------------------------------------------
   ! Helper: process cells on one Huygens face
   !--------------------------------------------------------------------------------
   attributes(global) subroutine process_face_cells(
        ExIz_d, ExDe_d, ExAb_d, ExAr_d, EyFr_d, EyTr_d, EyAb_d, EyAr_d,
        EzIz_d, EzDe_d, EzFr_d, EzTr_d,
        HxIz_d, HxDe_d, HxAb_d, HxAr_d, HyFr_d, HyTr_d, HyAb_d, HyAr_d,
        HzIz_d, HzDe_d, HzFr_d, HzTr_d,
        HxIz2_d, HxDe2_d, HxAb2_d, HxAr2_d, HyFr2_d, HyTr2_d, HyAb2_d, HyAr2_d,
        HzIz2_d, HzDe2_d, HzFr2_d, HzTr2_d,
        phys_x_Mx_d, phys_y_Mx_d, phys_z_Mx_d,
        phys_x_My_d, phys_y_My_d, phys_z_My_d,
        phys_x_Mz_d, phys_y_Mz_d, phys_z_Mz_d,
        phys_x_Jx_d, phys_y_Jx_d, phys_z_Jx_d,
        phys_x_Jy_d, phys_y_Jy_d, phys_z_Jy_d,
        phys_x_Jz_d, phys_y_Jz_d, phys_z_Jz_d,
        dyh_d, dze_d, dye_d, dzh_d,
        expIwdt_d, auxExp_E_d, auxExp_H_d,
        comun, sintheta_cosphi, sintheta_sinphi, costheta,
        costheta_cosphi, costheta_sinphi, sintheta, sinphi, cosphi,
        L_theta, L_phi, N_theta, N_phi,
        pasadas, normal, freqIdx,
        sym1, sym2, sym3, sym4, sym5,
        sym6, sym7, sym8, sym9, sym10,
        faceId, freq, cluz, z0
   )

      complex(kind=rkind), dimension(:,:,:), device, intent(in) :: ExIz, ExDe, ExAb, ExAr, EyFr, EyTr, EyAb, EyAr, EzIz, EzDe, EzFr, EzTr
      complex(kind=rkind), dimension(:,:,:), device, intent(in) :: HxIz, HxDe, HxAb, HxAr, HyFr, HyTr, HyAb, HyAr, HzIz, HzDe, HzFr, HzTr
      complex(kind=rkind), dimension(:,:,:), device, intent(in) :: HxIz2, HxDe2, HxAb2, HxAr2, HyFr2, HyTr2, HyAb2, HyAr2, HzIz2, HzDe2, HzFr2, HzTr2
      real(kind=rkind), dimension(:), device, intent(in) :: phys_x_Mx_d, phys_y_Mx_d, phys_z_Mx_d
      real(kind=rkind), dimension(:), device, intent(in) :: phys_x_My_d, phys_y_My_d, phys_z_My_d
      real(kind=rkind), dimension(:), device, intent(in) :: phys_x_Mz_d, phys_y_Mz_d, phys_z_Mz_d
      real(kind=rkind), dimension(:), device, intent(in) :: phys_x_Jx_d, phys_y_Jx_d, phys_z_Jx_d
      real(kind=rkind), dimension(:), device, intent(in) :: phys_x_Jy_d, phys_y_Jy_d, phys_z_Jy_d
      real(kind=rkind), dimension(:), device, intent(in) :: phys_x_Jz_d, phys_y_Jz_d, phys_z_Jz_d
      real(kind=rkind), dimension(:), device, intent(in) :: dyh_d, dze_d, dye_d, dzh_d
      complex(kind=rkind), dimension(:), device, intent(in) :: expIwdt_d, auxExp_E_d, auxExp_H_d
      real(kind=rkind), value :: comun, sintheta_cosphi, sintheta_sinphi, costheta
      real(kind=rkind), value :: costheta_cosphi, costheta_sinphi, sintheta, sinphi, cosphi
      real(kind=rkind), value :: freq, cluz, z0
      integer, value :: pasadas, normal, freqIdx, faceId
      logical, value :: sym1, sym2, sym3, sym4, sym5
      logical, value :: sym6, sym7, sym8, sym9, sym10
      complex(kind=rkind), intent(inout) :: L_theta, L_phi, N_theta, N_phi

      integer :: cellIdx
      complex(kind=rkind) :: Ez_cell, Ey_cell, Hz_cell, Hy_cell, Hz2_cell, Hy2_cell
      real(kind=rkind) :: dyh_val, dze_val, dye_val, dzh_val
      real(kind=rkind) :: x_My, y_My, z_My, x_Mz, y_Mz, z_Mz
      real(kind=rkind) :: x_Jy, y_Jy, z_Jy, x_Jz, y_Jz, z_Jz
      complex(kind=rkind) :: Mx, My, Mz, Jx, Jy, Jz
      complex(kind=rkind) :: comunMx, comunMy, comunMz, comunJx, comunJy, comunJz
      complex(kind=rkind) :: z1, z2, z_avg
      real(kind=rkind) :: phi1, phi2

      cellIdx = (blockidx%x - 1) * blockdim%x + threadidx%x

      ! Each thread processes one cell on the face
      ! Cell index mapping depends on face geometry
      ! For now, use a simple linear mapping (to be refined per face)
      if (cellIdx > 10000) return  ! Safety limit

      ! Read cell data from DFT buffers (simplified - actual indexing depends on face)
      ! This is a placeholder — actual implementation needs face-specific indexing
      if (faceId == 1) then
         ! Back face (Tr): uses EzTr, EyTr, HzTr, HyTr
         Ez_cell = EzTr(cellIdx, 0, freqIdx)
         Ey_cell = EyTr(cellIdx, 0, freqIdx)
         Hz_cell = HzTr(cellIdx, 0, freqIdx)
         Hy_cell = HyTr(cellIdx, 0, freqIdx)
         Hz2_cell = HzTr2(cellIdx, 0, freqIdx)
         Hy2_cell = HyTr2(cellIdx, 0, freqIdx)
      else
         ! Front face (Fr): uses EzFr, EyFr, HzFr, HyFr
         Ez_cell = EzFr(cellIdx, 0, freqIdx)
         Ey_cell = EyFr(cellIdx, 0, freqIdx)
         Hz_cell = HzFr(cellIdx, 0, freqIdx)
         Hy_cell = HyFr(cellIdx, 0, freqIdx)
         Hz2_cell = HzFr2(cellIdx, 0, freqIdx)
         Hy2_cell = HyFr2(cellIdx, 0, freqIdx)
      endif

      ! Read geometry and cell dimensions
      x_My = phys_x_My_d(cellIdx)
      y_My = phys_y_My_d(cellIdx)
      z_My = phys_z_My_d(cellIdx)
      x_Mz = phys_x_Mz_d(cellIdx)
      y_Mz = phys_y_Mz_d(cellIdx)
      z_Mz = phys_z_Mz_d(cellIdx)
      x_Jy = phys_x_Jy_d(cellIdx)
      y_Jy = phys_y_Jy_d(cellIdx)
      z_Jy = phys_z_Jy_d(cellIdx)
      x_Jz = phys_x_Jz_d(cellIdx)
      y_Jz = phys_y_Jz_d(cellIdx)
      z_Jz = phys_z_Jz_d(cellIdx)

      dyh_val = dyh_d(cellIdx)
      dze_val = dze_d(cellIdx)
      dye_val = dye_d(cellIdx)
      dzh_val = dzh_d(cellIdx)

      ! Compute equivalent currents (from farfield.F90 lines 2554-2557)
      My = -Ez_cell * dyh_val * dze_val * real(normal, RKIND)
      Mz = Ey_cell * dye_val * dzh_val * real(normal, RKIND)
      Jy = AverageNF2FF(pasadas, Hz_cell, Hz2_cell) * dye_val * dzh_val * real(normal, RKIND)
      Jz = -AverageNF2FF(pasadas, Hy_cell, Hy2_cell) * dyh_val * dze_val * real(normal, RKIND)
      Mx = (0.0_rkind, 0.0_rkind)
      Jx = (0.0_rkind, 0.0_rkind)

      ! Compute phase factors and update L, N (from update_LN, lines 3378-3389)
      comunMx = exp(comun * (x_My*sintheta_cosphi + y_My*sintheta_sinphi + z_My*costheta))
      comunMy = exp(comun * (x_Mz*sintheta_cosphi + y_Mz*sintheta_sinphi + z_Mz*costheta))
      comunMz = (0.0_rkind, 0.0_rkind)  ! Mz coordinate not used for Iz/De faces
      comunJx = (0.0_rkind, 0.0_rkind)
      comunJy = exp(comun * (x_Jy*sintheta_cosphi + y_Jy*sintheta_sinphi + z_Jy*costheta))
      comunJz = exp(comun * (x_Jz*sintheta_cosphi + y_Jz*sintheta_sinphi + z_Jz*costheta))

      L_theta = L_theta + (My * costheta_sinphi * comunMy - Mz * sintheta * comunMx)
      L_phi = L_phi + (-My * cosphi * comunMy)
      N_theta = N_theta + (Jy * costheta_sinphi * comunJy - Jz * sintheta * comunJz)
      N_phi = N_phi + (Jy * cosphi * comunJy)

      ! Apply symmetry clones (predicated execution)
      if (sym1 .or. sym2) then
         ! Clone with PEC symmetry
         call clone_and_update_LN(x_My, y_My, z_My, x_Mz, y_Mz, z_Mz, x_Jy, y_Jy, z_Jy, x_Jz, y_Jz, z_Jz, &
                                  My, Mz, Jy, Jz, comun, sintheta_cosphi, sintheta_sinphi, costheta, &
                                  sintheta, sinphi, cosphi, L_theta, L_phi, N_theta, N_phi, normal)
      endif
      ! ... (more symmetry clones)

   end subroutine process_face_cells

   !--------------------------------------------------------------------------------
   ! Helper: compute average (geometric or arithmetic) — matches farfield.F90:3498
   !--------------------------------------------------------------------------------
   attributes(device) function AverageNF2FF(pasadas, z1, z2) result(z)
      complex(kind=rkind), intent(in) :: z1, z2
      integer, intent(in) :: pasadas
      complex(kind=rkind) :: z
      real(kind=rkind) :: phi1, phi2, nphi1, nphi2

      z = (0.0_rkind, 0.0_rkind)
      if (pasadas == 2) then  ! geometric
         phi1 = atan2(imag(z1), real(z1))
         phi2 = atan2(imag(z2), real(z2))
         nphi1 = phi1
         nphi2 = phi2
         if ((phi1 < -pi/2.0_rkind) .and. (phi2 > pi/2.0_rkind)) nphi1 = phi1 - 2.0_rkind * pi
         if ((phi2 < -pi/2.0_rkind) .and. (phi1 > pi/2.0_rkind)) nphi2 = phi2 - 2.0_rkind * pi
         z = sqrt(abs(z1*z2)) * exp((0.0_rkind, 1.0_rkind) * (nphi1 + nphi2) / 2.0_rkind)
      elseif (pasadas == 1) then  ! arithmetic
         z = (z1 + z2) / 2.0_rkind
      endif
   end function AverageNF2FF

   !--------------------------------------------------------------------------------
   ! Public interface: GPU far-field flush
   !--------------------------------------------------------------------------------
   subroutine gpu_flush_nf2ff(this, Etheta_h, Ephi_h, RCS_h)
      class(gpu_state_t), intent(inout) :: this
      real(kind=rkind), dimension(:), intent(out) :: Etheta_h, Ephi_h, RCS_h

      integer(kind=4) :: totalAngles, blockSize, gridSize, cuda_status

      if (.not. this%nf2ff_initialized) return

      totalAngles = this%nf2ff_Ntheta * this%nf2ff_Nphi
      if (totalAngles == 0) return

      blockSize = 256
      gridSize = (totalAngles + blockSize - 1) / blockSize

      ! Launch kernel — each block handles one (freq, theta, phi)
      ! Note: This is a simplified launch configuration
      ! The full implementation would need face-specific cell iteration
      call gpu_flush_nf2ff_kernel<<<gridSize, blockSize>>>(
           this%nf2ff_ExIz_d, this%nf2ff_ExDe_d, this%nf2ff_ExAb_d, this%nf2ff_ExAr_d, &
           this%nf2ff_EyFr_d, this%nf2ff_EyTr_d, this%nf2ff_EyAb_d, this%nf2ff_EyAr_d, &
           this%nf2ff_EzIz_d, this%nf2ff_EzDe_d, this%nf2ff_EzFr_d, this%nf2ff_EzTr_d, &
           this%nf2ff_HxIz_d, this%nf2ff_HxDe_d, this%nf2ff_HxAb_d, this%nf2ff_HxAr_d, &
           this%nf2ff_HyFr_d, this%nf2ff_HyTr_d, this%nf2ff_HyAb_d, this%nf2ff_HyAr_d, &
           this%nf2ff_HzIz_d, this%nf2ff_HzDe_d, this%nf2ff_HzFr_d, this%nf2ff_HzTr_d, &
           this%nf2ff_HxIz2_d, this%nf2ff_HxDe2_d, this%nf2ff_HxAb2_d, this%nf2ff_HxAr2_d, &
           this%nf2ff_HyFr2_d, this%nf2ff_HyTr2_d, this%nf2ff_HyAb2_d, this%nf2ff_HyAr2_d, &
           this%nf2ff_HzIz2_d, this%nf2ff_HzDe2_d, this%nf2ff_HzFr2_d, this%nf2ff_HzTr2_d, &
           this%nf2ff_phys_x_Mx_d, this%nf2ff_phys_y_Mx_d, this%nf2ff_phys_z_Mx_d, &
           this%nf2ff_phys_x_My_d, this%nf2ff_phys_y_My_d, this%nf2ff_phys_z_My_d, &
           this%nf2ff_phys_x_Mz_d, this%nf2ff_phys_y_Mz_d, this%nf2ff_phys_z_Mz_d, &
           this%nf2ff_phys_x_Jx_d, this%nf2ff_phys_y_Jx_d, this%nf2ff_phys_z_Jx_d, &
           this%nf2ff_phys_x_Jy_d, this%nf2ff_phys_y_Jy_d, this%nf2ff_phys_z_Jy_d, &
           this%nf2ff_phys_x_Jz_d, this%nf2ff_phys_y_Jz_d, this%nf2ff_phys_z_Jz_d, &
           this%nf2ff_dyh_d, this%nf2ff_dze_d, this%nf2ff_dye_d, this%nf2ff_dzh_d, &
           this%nf2ff_expIwdt_d, this%nf2ff_auxExp_E_d, this%nf2ff_auxExp_H_d, &
           this%nf2ff_Etheta_d, this%nf2ff_Ephi_d, this%nf2ff_RCS_d, &
           this%nf2ff_num_freqs, this%nf2ff_Ntheta, this%nf2ff_Nphi, &
           0, 0, this%nf2ff_thetaStep, 0, 0, this%nf2ff_phiStep, &
           this%nf2ff_freqStep, this%nf2ff_initialFreq, 0.0_rkind, this%nf2ff_cluz, this%nf2ff_z0, &
           this%nf2ff_XDobleAncho, this%nf2ff_YDobleAncho, this%nf2ff_ZDobleAncho, &
           this%nf2ff_sym_flags, &
           .true., .true., .true., .true., .true., .true., &
           .false., .false., .false., .false., .false., &
           .false., .false., .false., .false., .false., &
           .false., .false., .false., .false., .false., &
           .false., .false., .false., .false., .false., &
           .false., .false., .false., .false., .false., &
           .false., .false., .false., .false., .false., &
           .false., .false., .false., .false., .false., &
           .false., .false., .false., .false., .false.
      )

      ! D2H transfer results
      cuda_status = cudaMemcpy(Etheta_h, this%nf2ff_Etheta_d, &
                               this%nf2ff_num_freqs * this%nf2ff_Ntheta * this%nf2ff_Nphi * 4, &
                               cudaMemcpyDeviceToHost)
      cuda_status = cudaMemcpy(Ephi_h, this%nf2ff_Ephi_d, &
                               this%nf2ff_num_freqs * this%nf2ff_Ntheta * this%nf2ff_Nphi * 4, &
                               cudaMemcpyDeviceToHost)
      cuda_status = cudaMemcpy(RCS_h, this%nf2ff_RCS_d, &
                               this%nf2ff_num_freqs * this%nf2ff_Ntheta * this%nf2ff_Nphi * 4, &
                               cudaMemcpyDeviceToHost)

   end subroutine gpu_flush_nf2ff

end module gpu_nf2ff_m
