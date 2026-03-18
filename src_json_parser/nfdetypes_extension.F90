module NFDETypes_extension_m      
#ifdef CompileWithSMBJSON
   use NFDETypes_m

   implicit none

   public

   interface operator(==)
      module procedure Parseador_eq

      module procedure NFDEGeneral_eq
      module procedure desplazamiento_eq

      module procedure MatrizMedios_eq
      module procedure Material_eq
      module procedure Materials_eq
      module procedure pecregions_eq
      module procedure dielectric_eq
      module procedure dielectricregions_eq
      module procedure freqdepenmaterial_eq
      module procedure freqdepenmaterials_eq
      module procedure anisotropicbody_eq
      module procedure anisotropicelements_eq
      module procedure LossyThinSurface_eq
      module procedure LossyThinSurfaces_eq

      module procedure frontera_eq
      module procedure fronteraPML_eq

      module procedure box_eq
      module procedure boxes_eq
      module procedure planewaves_eq
      module procedure planewave_eq
      module procedure curr_field_src_eq
      module procedure nodsource_eq

      module procedure coords_eq
      module procedure coords_scaled_eq
      module procedure abstractSonda_eq
      module procedure sonda_eq
      module procedure sondas_eq
      module procedure massonda_eq
      module procedure massondas_eq
      module procedure bloqueprobe_eq
      module procedure bloqueprobes_eq
      module procedure volprobe_eq
      module procedure volprobes_eq

      module procedure FarField_Sonda_eq
      module procedure Electric_Sonda_eq
      module procedure Magnetic_Sonda_eq
      module procedure NormalElectric_Sonda_eq
      module procedure NormalMagnetic_Sonda_eq
      module procedure SurfaceElectricCurrent_Sonda_eq
      module procedure SurfaceMagneticCurrent_Sonda_eq
      module procedure ThinWireComp_eq
      module procedure ThinWire_eq
      module procedure ThinWires_eq
      module procedure SlantedWireComp_eq
      module procedure SlantedWire_eq
      module procedure SlantedWires_eq
      module procedure ThinSlotComp_eq
      module procedure ThinSlot_eq
      module procedure ThinSlots_eq
   end interface

contains
   subroutine initializeProblemDescription(pD)
      type(Parseador_t), intent(inout) :: pD

      allocate(pD%general)
      allocate(pD%matriz)
      allocate(pD%despl)
      allocate(pD%front)

      allocate(pD%mats)
      pD%mats%n_Mats = 3
      pD%mats%n_Mats_max = 3
      allocate(pD%mats%mats(3))

      pD%mats%mats(1)%id = 1
      pD%mats%mats(1)%eps = EPSILON_VACUUM
      pD%mats%mats(1)%mu = MU_VACUUM
      pD%mats%mats(1)%sigma = 0.0
      pD%mats%mats(1)%sigmam = 0.0

      pD%mats%mats(2)%id = 2
      pD%mats%mats(2)%eps = EPSILON_VACUUM
      pD%mats%mats(2)%mu = MU_VACUUM
      pD%mats%mats(2)%sigma = SIGMA_PEC
      pD%mats%mats(2)%sigmam = 0.0

      pD%mats%mats(3)%id = 3
      pD%mats%mats(3)%eps = EPSILON_VACUUM
      pD%mats%mats(3)%mu = MU_VACUUM
      pD%mats%mats(3)%sigma = 0.0
      pD%mats%mats(3)%sigmam =SIGMA_PMC

      allocate(pD%pecRegs)
      allocate(pD%pecRegs%lins(0))
      allocate(pD%pecRegs%surfs(0))
      allocate(pD%pecRegs%vols(0))

      allocate(pD%pmcRegs)
      allocate(pD%pmcRegs%lins(0))
      allocate(pD%pmcRegs%surfs(0))
      allocate(pD%pmcRegs%vols(0))

      allocate(pD%DielRegs)
      allocate(pD%DielRegs%lins(0))
      allocate(pD%DielRegs%surfs(0))
      allocate(pD%DielRegs%vols(0))
      
      allocate(pD%LossyThinSurfs)
      allocate(pD%LossyThinSurfs%cs(0))

      allocate(pD%frqDepMats)
      allocate(pD%aniMats)
      !
      allocate(pD%plnSrc)
      allocate(pD%plnSrc%collection(0))
      allocate(pD%nodSrc)
      allocate(pD%nodSrc%NodalSource(0))
      allocate(pD%boxSrc)
      !
      allocate(pD%Sonda)
      pD%Sonda%len_cor_max = 0
      pD%Sonda%length = 0
      pD%Sonda%length_max = 0
      allocate(pD%Sonda%collection(0))
      !
      allocate(pD%oldSONDA)
      allocate(pD%oldSONDA%probes(0))
      pD%oldSONDA%n_probes = 0
      pD%oldSONDA%n_probes_max = 0

      allocate(pD%BloquePrb)
      allocate(pD%BloquePrb%bp(0))
      
      allocate(pD%VolPrb)
      allocate(pD%VolPrb%collection(0))
      !
      allocate(pD%tWires)
      allocate(pD%tWires%tw(0))
      pD%tWires%n_tw = 0
      pD%tWires%n_tw = 0

      allocate(pD%sWires)
      allocate(pD%tSlots)
      allocate(pD%tSlots%tg(0))

#ifdef CompileWithMTLN
      allocate(pD%mtln)
      allocate(pD%mtln%cables(0))
      allocate(pD%mtln%probes(0))
      allocate(pD%mtln%networks(0))
      allocate(pD%mtln%connectors(0))
#endif

      allocate(pD%conformalRegs)
      allocate(pD%conformalRegs%volumes(0))
      allocate(pD%conformalRegs%surfaces(0))


   end subroutine

   elemental logical function parseador_eq(a, b)
      type(Parseador_t), intent(in) :: a, b
      parseador_eq = &
         (a%switches == b%switches) .and. &
         (a%general        == b%general) .and. &
         (a%matriz         == b%matriz) .and. &
         (a%despl          == b%despl) .and. &
         (a%front          == b%front) .and. &
         (a%Mats           == b%Mats) .and. &
         (a%pecRegs        == b%pecRegs) .and. &
         (a%pmcRegs        == b%pmcRegs) .and. &
         (a%DielRegs       == b%DielRegs) .and. &
         (a%LossyThinSurfs == b%LossyThinSurfs) .and. &
         (a%frqDepMats     == b%frqDepMats) .and. &
         (a%aniMats        == b%aniMats) .and. &
         (a%boxSrc         == b%boxSrc) .and. &
         (a%plnSrc         == b%plnSrc) .and. &
         (a%nodSrc         == b%nodSrc) .and. &
         (a%oldSONDA       == b%oldSONDA) .and. &
         (a%Sonda          == b%Sonda) .and. &
         (a%BloquePrb      == b%BloquePrb) .and. &
         (a%VolPrb         == b%VolPrb) .and. &
         (a%tWires         == b%tWires) .and. &
         (a%sWires         == b%sWires) .and. &
         (a%tSlots         == b%tSlots)
   end function parseador_eq

   elemental logical function MatrizMedios_eq(a, b)
      type(MatrizMedios_t), intent(in) :: a, b
      MatrizMedios_eq = &
         a%totalX == b%totalX .and. &
         a%totalY == b%totalY .and. &
         a%totalZ == b%totalZ
   end function MatrizMedios_eq

   elemental logical function Material_eq(a, b)
      type(Material_t), intent(in) :: a, b
      Material_eq = &
         a%eps == b%eps .and. &
         a%mu == b%mu .and. &
         a%sigma == b%sigma .and. &
         a%sigmam == b%sigmam .and. &
         a%id == b%id
   end function Material_eq

   elemental logical function Materials_eq(a, b)
      type(Materials_t), intent(in) :: a, b

      Materials_eq = &
         (a%n_Mats == b%n_Mats) .and. &
         (a%n_Mats_max == b%n_Mats_max) .and. &
         all(a%Mats == b%Mats)
   end function Materials_eq

   elemental logical function pecregions_eq(a, b)
      type(PECRegions_t), intent(in) :: a, b
      logical :: lins_ok, surfs_ok, vols_ok

      if (a%nLins /= b%nLins .or. a%nLins_max /= b%nLins_max .or. &
          a%nSurfs /= b%nSurfs .or. a%nSurfs_max /= b%nSurfs_max .or. &
          a%nVols /= b%nVols .or. a%nVols_max /= b%nVols_max) then
         pecregions_eq = .false.
         return
      end if

      lins_ok = .false.
      if (associated(a%Lins) .eqv. associated(b%Lins)) then
         if (associated(a%Lins)) then
            lins_ok = all(a%Lins == b%Lins)
         else
            lins_ok = .true.
         end if
      else if (associated(a%Lins)) then
         if (size(a%Lins) == 0) lins_ok = .true.
      else
         if (size(b%Lins) == 0) lins_ok = .true.
      end if

      surfs_ok = .false.
      if (associated(a%Surfs) .eqv. associated(b%Surfs)) then
         if (associated(a%Surfs)) then
            surfs_ok = all(a%Surfs == b%Surfs)
         else
            surfs_ok = .true.
         end if
      else if (associated(a%Surfs)) then
         if (size(a%Surfs) == 0) surfs_ok = .true.
      else
         if (size(b%Surfs) == 0) surfs_ok = .true.
      end if

      vols_ok = .false.
      if (associated(a%Vols) .eqv. associated(b%Vols)) then
         if (associated(a%Vols)) then
            vols_ok = all(a%Vols == b%Vols)
         else
            vols_ok = .true.
         end if
      else if (associated(a%Vols)) then
         if (size(a%Vols) == 0) vols_ok = .true.
      else
         if (size(b%Vols) == 0) vols_ok = .true.
      end if

      pecregions_eq = lins_ok .and. surfs_ok .and. vols_ok

   end function pecregions_eq

   elemental logical function dielectric_eq(a, b)
      type(dielectric_t), intent(in) :: a, b
      logical :: allAssociated
      dielectric_eq = &
         (a%n_C1P == b%n_C1P) .and. &
         (a%n_C2P == b%n_C2P) .and. &
         (a%sigma == b%sigma) .and. &
         (a%eps == b%eps) .and. &
         (a%mu == b%mu) .and. &
         (a%sigmam == b%sigmam) .and. &
         (a%Rtime_on == b%Rtime_on) .and. &
         (a%Rtime_off == b%Rtime_off) .and. &
         (a%R == b%R) .and. &
         (a%L == b%L) .and. &
         (a%C == b%C) .and. &
         (a%R_devia == b%R_devia) .and. &
         (a%L_devia == b%L_devia) .and. &
         (a%C_devia == b%C_devia) .and. &
         (a%DiodB == b%DiodB) .and. &
         (a%DiodIsat == b%DiodIsat) .and. &
         (a%DiodOri == b%DiodOri) .and. &
         (a%orient == b%orient) .and. &
         (a%resistor  .eqv. b%resistor) .and. &
         (a%inductor  .eqv. b%inductor) .and. &
         (a%capacitor .eqv. b%capacitor) .and. &
         (a%diodo     .eqv. b%diodo) .and. &
         (a%plain     .eqv. b%plain) .and. &
         (a%PMLbody   .eqv. b%PMLbody)
      
      allAssociated = &
         associated(a%C1P) .and. associated(b%C1P) .and. &
         associated(a%C2P) .and. associated(b%C2P) 
      if (.not. allAssociated) then
         dielectric_eq = .false.
         return
      end if
      dielectric_eq = all(a%C1P == b%C1P) .and. all(a%C2P == b%C2P)
   end function dielectric_eq

   elemental logical function freqdepenmaterial_eq(a, b)
      type(FreqDepenMaterial_t), intent(in) :: a, b
      freqdepenmaterial_eq = &
         all(a%a11 == b%a11) .and. &
         all(a%b11 == b%b11) .and. &
         all(a%am11 == b%am11) .and. &
         all(a%bm11 == b%bm11) .and. &
         all(a%a12 == b%a12) .and. &
         all(a%b12 == b%b12) .and. &
         all(a%am12 == b%am12) .and. &
         all(a%bm12 == b%bm12) .and. &
         all(a%a13 == b%a13) .and. &
         all(a%b13 == b%b13) .and. &
         all(a%am13 == b%am13) .and. &
         all(a%bm13 == b%bm13) .and. &
         all(a%a22 == b%a22) .and. &
         all(a%b22 == b%b22) .and. &
         all(a%am22 == b%am22) .and. &
         all(a%bm22 == b%bm22) .and. &
         all(a%a23 == b%a23) .and. &
         all(a%b23 == b%b23) .and. &
         all(a%am23 == b%am23) .and. &
         all(a%bm23 == b%bm23) .and. &
         all(a%a33 == b%a33) .and. &
         all(a%b33 == b%b33) .and. &
         all(a%am33 == b%am33) .and. &
         all(a%bm33 == b%bm33) .and. &
         all(a%alpha == b%alpha) .and. &
         all(a%beta == b%beta) .and. &
         all(a%gamma == b%gamma) .and. &
         all(a%alpham == b%alpham) .and. &
         all(a%betam == b%betam) .and. &
         all(a%gammam == b%gammam) .and. &
         all(coords_eq(a%c(1:a%n_c), b%c(1:b%n_c))) .and. &
         (a%eps11 == b%eps11) .and. &
         (a%eps12 == b%eps12) .and. &
         (a%eps13 == b%eps13) .and. &
         (a%eps22 == b%eps22) .and. &
         (a%eps23 == b%eps23) .and. &
         (a%eps33 == b%eps33) .and. &
         (a%mu11 == b%mu11) .and. &
         (a%mu12 == b%mu12) .and. &
         (a%mu13 == b%mu13) .and. &
         (a%mu22 == b%mu22) .and. &
         (a%mu23 == b%mu23) .and. &
         (a%mu33 == b%mu33) .and. &
         (a%sigma11 == b%sigma11) .and. &
         (a%sigma12 == b%sigma12) .and. &
         (a%sigma13 == b%sigma13) .and. &
         (a%sigma22 == b%sigma22) .and. &
         (a%sigma23 == b%sigma23) .and. &
         (a%sigma33 == b%sigma33) .and. &
         (a%sigmam11 == b%sigmam11) .and. &
         (a%sigmam12 == b%sigmam12) .and. &
         (a%sigmam13 == b%sigmam13) .and. &
         (a%sigmam22 == b%sigmam22) .and. &
         (a%sigmam23 == b%sigmam23) .and. &
         (a%sigmam33 == b%sigmam33) .and. &
         (a%K11 == b%K11) .and. &
         (a%Km11 == b%Km11) .and. &
         (a%K12 == b%K12) .and. &
         (a%Km12 == b%Km12) .and. &
         (a%K13 == b%K13) .and. &
         (a%Km13 == b%Km13) .and. &
         (a%K22 == b%K22) .and. &
         (a%Km22 == b%Km22) .and. &
         (a%K23 == b%K23) .and. &
         (a%Km23 == b%Km23) .and. &
         (a%K33 == b%K33) .and. &
         (a%Km33 == b%Km33) .and. &
         (a%L == b%L) .and. &
         (a%Lm == b%Lm) .and. &
         (a%n_c == b%n_c) .and. &
         (a%files == b%files)
   end function freqdepenmaterial_eq

   elemental logical function freqdepenmaterials_eq(a, b)
      type(FreqDepenMaterials_t), intent(in) :: a, b
      freqdepenmaterials_eq = &
         (a%nVols == b%nVols) .and. &
         (a%nSurfs == b%nSurfs) .and. &
         (a%nLins == b%nLins) .and. &
         (a%nVols_max == b%nVols_max) .and. &
         (a%nSurfs_max == b%nSurfs_max) .and. &
         (a%nLins_max == b%nLins_max) .and. &
         (a%n_c_max == b%n_c_max) .and. &
         all(a%Vols  == b%Vols) .and. &
         all(a%Surfs == b%Surfs) .and. &
         all(a%Lins  == b%Lins)
   end function freqdepenmaterials_eq

   elemental logical function anisotropicbody_eq(a, b)
      type(ANISOTROPICbody_t), intent(in) :: a, b
      anisotropicbody_eq = &
         all(a%c1P == b%c1P) .and. &
         all(a%c2P == b%c2P) .and. &
         all(a%sigma == b%sigma) .and. &
         all(a%eps == b%eps) .and. &
         all(a%mu == b%mu) .and. &
         all(a%sigmam == b%sigmam) .and. &
         (a%n_C1P == b%n_C1P) .and. &
         (a%n_C2P == b%n_C2P)
   end function anisotropicbody_eq

   elemental logical function anisotropicelements_eq(a, b)
      type(ANISOTROPICelements_t), intent(in) :: a, b
      anisotropicelements_eq = &
         (a%nVols == b%nVols) .and. &
         (a%nSurfs == b%nSurfs) .and. &
         (a%nLins == b%nLins) .and. &
         (a%nVols_max == b%nVols_max) .and. &
         (a%nSurfs_max == b%nSurfs_max) .and. &
         (a%nLins_max == b%nLins_max) .and. &
         (a%n_C1P_max == b%n_C1P_max) .and. &
         (a%n_C2P_max == b%n_C2P_max) .and. &
         all(a%Vols  == b%Vols ) .and. &
         all(a%Surfs == b%Surfs) .and. &
         all(a%Lins  == b%Lins )
   end function anisotropicelements_eq

   elemental logical function dielectricregions_eq(a, b)
      type(DielectricRegions_t), intent(in) :: a, b
      logical :: allAssociated
      
      allAssociated = &
         associated(a%Lins) .and. associated(b%Lins) .and. &
         associated(a%Surfs) .and. associated(b%Surfs) .and. &
         associated(a%Vols) .and. associated(b%Vols)
         if (.not. allAssociated) then
            dielectricregions_eq = .false.
            return
         end if

      dielectricregions_eq = &
         a%nVols == b%nVols .and. &
         a%nSurfs == b%nSurfs .and. &
         a%nLins == b%nLins .and. &
         a%nVols_max == b%nVols_max .and. &
         a%nSurfs_max == b%nSurfs_max .and. &
         a%nLins_max == b%nLins_max .and. &
         a%n_C1P_max == b%n_C1P_max .and. &
         a%n_C2P_max == b%n_C2P_max

      dielectricregions_eq = dielectricregions_eq .and. all(a%Lins == b%Lins) 
      dielectricregions_eq = dielectricregions_eq .and. all(a%surfs == b%surfs) 
      dielectricregions_eq = dielectricregions_eq .and. all(a%vols == b%vols) 
   end function dielectricregions_eq

   elemental logical function LossyThinSurface_eq(a, b)
      type(LossyThinSurface_t), intent(in) :: a, b
      LossyThinSurface_eq = &
         all(a%c == b%c) .and. &
         all(a%sigma == b%sigma) .and. &
         all(a%eps == b%eps) .and. &
         all(a%mu == b%mu) .and. &
         all(a%sigmam == b%sigmam) .and. &
         all(a%thk == b%thk) .and. &
         all(a%sigma_devia == b%sigma_devia) .and. &
         all(a%eps_devia == b%eps_devia) .and. &
         all(a%mu_devia == b%mu_devia) .and. &
         all(a%sigmam_devia == b%sigmam_devia) .and. &
         all(a%thk_devia == b%thk_devia) .and. &
         (a%nc == b%nc) .and. &
         (a%files == b%files) .and. &
         (a%numcapas == b%numcapas)
   end function LossyThinSurface_eq

   elemental logical function LossyThinSurfaces_eq(a, b)
      type(LossyThinSurfaces_t), intent(in) :: a, b

      LossyThinSurfaces_eq = all(a%cs == b%cs) .and. &
         (a%length == b%length) .and. &
         (a%length_max == b%length_max) .and. &
         (a%nC_max == b%nC_max)
   end function LossyThinSurfaces_eq

   elemental logical function ThinWireComp_eq(a, b)
      type(ThinWireComp_t), intent(in) :: a, b

      ThinWireComp_eq = &
         (a%srctype == b%srctype) .and. &
         (a%srcfile == b%srcfile) .and. &
         (a%i == b%i) .and. &
         (a%j == b%j) .and. &
         (a%K == b%K) .and. &
         (a%nd == b%nd) .and. &
         (a%d == b%d) .and. &
         (a%m == b%m) .and. &
         (a%tag == b%tag)
   end function ThinWireComp_eq

   elemental logical function ThinWire_eq(a, b)
      type(ThinWire_t), intent(in) :: a, b

      ThinWire_eq = .false.
      if (a%rad /= b%rad) return
      if (a%disp .neqv. b%disp) return
      if (a%dispfile /= b%dispfile) return
      if (a%res /= b%res) return
      if (a%ind /= b%ind) return
      if (a%cap /= b%cap) return
      if (a%P_res /= b%P_res) return
      if (a%P_ind /= b%P_ind) return
      if (a%P_cap /= b%P_cap) return
      if (a%dispfile_LeftEnd /= b%dispfile_LeftEnd) return
      if (a%R_LeftEnd /= b%R_LeftEnd) return
      if (a%L_LeftEnd /= b%L_LeftEnd) return
      if (a%C_LeftEnd /= b%C_LeftEnd) return
      if (a%dispfile_RightEnd /= b%dispfile_RightEnd) return
      if (a%R_RightEnd /= b%R_RightEnd) return
      if (a%L_RightEnd /= b%L_RightEnd) return
      if (a%C_RightEnd /= b%C_RightEnd) return
      if (a%LeftEnd /= b%LeftEnd) return
      if (a%RightEnd /= b%RightEnd) return
      if (a%tl /= b%tl) return
      if (a%tr /= b%tr) return
      if (.not. all(a%twc == b%twc)) return
      if (a%n_twc /= b%n_twc) return
      if (a%n_twc_max /= b%n_twc_max) return
      ThinWire_eq = .true.
   end function ThinWire_eq

   elemental logical function ThinWires_eq(a, b)
      type(ThinWires_t), intent(in) :: a, b
      ThinWires_eq = .false.
      if (a%n_tw /= b%n_tw) return
      if (a%n_tw_max /= b%n_tw_max) return
      if (associated(a%tw) .eqv. associated(b%tw)) then
         if (associated(a%tw)) then
            if (.not. all(a%tw == b%tw)) return
         end if
      else if (associated(a%tw)) then
         if (size(a%tw) /= 0) return
      else
         if (size(b%tw) /= 0) return
      end if
      ThinWires_eq = .true.
   end function ThinWires_eq

   elemental logical function SlantedWireComp_eq(a, b)
      type(SlantedWireComp_t), intent(in) :: a, b

      SlantedWireComp_eq = &
         (a%srctype == b%srctype) .and. &
         (a%srcfile == b%srcfile) .and. &
         (a%x == b%x) .and. &
         (a%y == b%y) .and. &
         (a%z == b%z) .and. &
         (a%nd == b%nd) .and. &
         (a%m == b%m) .and. &
         (a%tag == b%tag)
   end function SlantedWireComp_eq

   elemental logical function SlantedWire_eq(a, b)
      type(SlantedWire_t), intent(in) :: a, b

      SlantedWire_eq = &
         (a%rad == b%rad) .and. &
         (a%disp .eqv. b%disp) .and. &
         (a%dispfile == b%dispfile) .and. &
         (a%res == b%res) .and. &
         (a%ind == b%ind) .and. &
         (a%cap == b%cap) .and. &
         (a%P_res == b%P_res) .and. &
         (a%P_ind == b%P_ind) .and. &
         (a%P_cap == b%P_cap) .and. &
         (a%dispfile_LeftEnd == b%dispfile_LeftEnd) .and. &
         (a%R_LeftEnd == b%R_LeftEnd) .and. &
         (a%L_LeftEnd == b%L_LeftEnd) .and. &
         (a%C_LeftEnd == b%C_LeftEnd) .and. &
         (a%dispfile_RightEnd == b%dispfile_RightEnd) .and. &
         (a%R_RightEnd == b%R_RightEnd) .and. &
         (a%L_RightEnd == b%L_RightEnd) .and. &
         (a%C_RightEnd == b%C_RightEnd) .and. &
         (a%LeftEnd == b%LeftEnd) .and. &
         (a%RightEnd == b%RightEnd) .and. &
         (a%tl == b%tl) .and. &
         (a%tr == b%tr) .and. &
         (a%n_swc == b%n_swc) .and. &
         (a%n_swc_max == b%n_swc_max)
   end function SlantedWire_eq

   elemental logical function SlantedWires_eq(a, b)
      type(SlantedWiresInfo_t), intent(in) :: a, b
      logical :: allAssociated
      allAssociated = &
         associated(a%sw) .and. associated(b%sw)
      if (.not. allAssociated) then
         SlantedWires_eq = .false.
         return
      end if
      SlantedWires_eq = &
         all(a%sw == b%sw) .and. &
         (a%n_sw == b%n_sw) .and. &
         (a%n_sw_max == b%n_sw_max)
   end function SlantedWires_eq

   elemental logical function ThinSlotComp_eq(a, b)
      type(ThinSlotComp_t), intent(in) :: a, b

      ThinSlotComp_eq = (a%i == b%i) .and. &
         (a%j == b%j) .and. &
         (a%K == b%K) .and. &
         (a%node == b%node) .and. &
         (a%dir == b%dir) .and. &
         (a%Or == b%Or) .and. &
         (a%tag == b%tag)
   end function ThinSlotComp_eq

   elemental logical function ThinSlot_eq(a, b)
      type(ThinSlot_t), intent(in) :: a, b
      ThinSlot_eq = .false.
      if (a%width /= b%width) return
      if (a%n_tgc /= b%n_tgc) return
      if (a%n_tgc_max /= b%n_tgc_max) return
      if (associated(a%tgc) .eqv. associated(b%tgc)) then
         if (associated(a%tgc)) then
            if (.not. all(a%tgc == b%tgc)) return
         end if
      else if (associated(a%tgc)) then
         if (size(a%tgc) /= 0) return
      else
         if (size(b%tgc) /= 0) return
      end if
      ThinSlot_eq = .true.
   end function ThinSlot_eq

   elemental logical function ThinSlots_eq(a, b)
      type(ThinSlots_t), intent(in) :: a, b
      ThinSlots_eq = .false.
      if (a%n_tg /= b%n_tg) return
      if (a%n_tg_max /= b%n_tg_max) return
      if (associated(a%tg) .eqv. associated(b%tg)) then
         if (associated(a%tg)) then
            if (.not. all(a%tg == b%tg)) return
         end if
      else if (associated(a%tg)) then
         if (size(a%tg) /= 0) return
      else
         if (size(b%tg) /= 0) return
      end if
      ThinSlots_eq = .true.
   end function ThinSlots_eq

   elemental logical function NFDEGeneral_eq(a, b) result (res)
      type(NFDEGeneral_t), intent(in) :: a, b
      res = .false.
      if (a%dt /= b%dt) return
      if (a%nmax /= b%nmax) return
      res = .true.
   end function

   elemental logical function box_eq(a, b)
      type(Box_t), intent(in) :: a, b
      box_eq = &
         (a%nombre_fichero == b%nombre_fichero) .and. &
         all(a%coor1 == b%coor1) .and. &
         all(a%coor2 == b%coor2)
   end function box_eq

   elemental logical function boxes_eq(a, b)
      type(Boxes_t), intent(in) :: a, b
      boxes_eq = &
         (a%nVols == b%nVols) .and. &
         (a%nVols_max == b%nVols_max) .and. &
         all(a%Vols == b%Vols)
   end function boxes_eq

   elemental logical function planewave_eq(a,b) result(res)
      type(PlaneWave_t), intent(in) :: a, b
      res = .false.
      if (a%nombre_fichero /= b%nombre_fichero) return
      if (a%atributo /= b%atributo) return
      if (any(a%coor1 /= b%coor1)) return
      if (any(a%coor2 /= b%coor2)) return
      if (a%theta /= b%theta) return
      if (a%phi /= b%phi) return
      if (a%alpha /= b%alpha) return
      if (a%beta /= b%beta) return
      if (a%isRC .neqv. b%isRC) return
      if (a%INCERTMAX /= b%INCERTMAX) return
      if (a%numModes /= b%numModes) return
      res = .true.
   end function

   elemental logical function planewaves_eq(a,b) result(res)
      type(PlaneWaves_t), intent(in) :: a, b
      res = .false.
      if (a%nc /= b%nc) return
      if (a%nC_max /= b%nC_max) return
      if (associated(a%collection) .eqv. associated(b%collection)) then
         if (associated(a%collection)) then
            if (any(.not. a%collection == b%collection)) return
         end if
      else if (associated(a%collection)) then
         if (size(a%collection) /= 0) return
      else
         if (size(b%collection) /= 0) return
      end if
      res = .true.
   end function

   elemental logical function curr_field_src_eq(a, b)
      type(Curr_Field_Src_t), intent(in) :: a, b
      curr_field_src_eq = &
         all(a%c1P == b%c1P) .and. &
         all(a%c2P == b%c2P) .and. &
         (a%n_C1P == b%n_C1P) .and. &
         (a%n_C2P == b%n_C2P) .and. &
         (a%nombre == b%nombre) .and. &
         (a%isElec .eqv. b%isElec) .and. &
         (a%isHard .eqv. b%isHard) .and. &
         (a%isInitialValue .eqv. b%isInitialValue)
   end function curr_field_src_eq

   elemental logical function nodsource_eq(a, b)
      type(NodSource_t), intent(in) :: a, b
      nodsource_eq = .false.
      if (a%n_nodSrc /= b%n_nodSrc) return
      if (a%n_nodSrc_max /= b%n_nodSrc_max) return
      if (a%n_C1P_max /= b%n_C1P_max) return
      if (a%n_C2P_max /= b%n_C2P_max) return
      if (associated(a%NodalSource) .eqv. associated(b%NodalSource)) then
         if (associated(a%NodalSource)) then
            if (any(.not. a%NodalSource == b%NodalSource)) return
         end if
      else if (associated(a%NodalSource)) then
         if (size(a%NodalSource) /= 0) return
      else
         if (size(b%NodalSource) /= 0) return
      end if
      nodsource_eq = .true.
   end function nodsource_eq

   elemental logical function fronteraPML_eq(a, b) result (res)
      type(FronteraPML_t), intent(in) :: a, b
      res = .false.
      if (a%orden    /= b%orden) return
      if (a%refl     /= b%refl) return
      if (a%numCapas /= b%numCapas) return
      res = .true.
   end function

   elemental logical function frontera_eq(a, b) result (res)
      type(Frontera_t), intent(in) :: a, b
      integer :: i
      res = .false.
      if (any(a%tipoFrontera /= b%tipoFrontera)) return
      if (any(.not. a%propiedadesPML == b%propiedadesPML)) return
      res = .true.
   end function

   elemental logical function desplazamiento_eq(a, b) result (res)
      type(Desplazamiento_t), intent(in) :: a, b

      res = .false.

      if (.not. associated(a%desX)) return
      if (.not. associated(a%desY)) return
      if (.not. associated(a%desZ)) return
      if (.not. associated(b%desX)) return
      if (.not. associated(b%desY)) return
      if (.not. associated(b%desZ)) return

      if (size(a%desX) /= size(b%desX)) return
      if (size(a%desY) /= size(b%desY)) return
      if (size(a%desZ) /= size(b%desZ)) return
      if (any(a%desX /= b%desX)) return
      if (any(a%desY /= b%desY)) return
      if (any(a%desZ /= b%desZ)) return

      if (a%nX /=  b%nX) return
      if (a%nY /=  b%nY) return
      if (a%nZ /=  b%nZ) return

      if (a%mx1 /=  b%mx1) return
      if (a%my1 /=  b%my1) return
      if (a%mz1 /=  b%mz1) return

      if (a%mx2 /=  b%mx2) return
      if (a%my2 /=  b%my2) return
      if (a%mz2 /=  b%mz2) return

      if (a%originx /=  b%originx) return
      if (a%originy /=  b%originy) return
      if (a%originz /=  b%originz) return

      res = .true.
   end function

   elemental logical function coords_eq(a, b) result(res)
      type(coords_t), intent(in) :: a, b
      res = .false.
      if (a%xi /= b%xi) return
      if (a%xe /= b%xe) return
      if (a%yi /= b%yi) return
      if (a%ye /= b%ye) return
      if (a%zi /= b%zi) return
      if (a%ze /= b%ze) return
      if (a%xtrancos /= b%xtrancos) return
      if (a%ytrancos /= b%ytrancos) return
      if (a%ztrancos /= b%ztrancos) return
      if (a%or /= b%or) return
      if (a%tag /= b%tag) return
      res = .true.
   end function

   elemental logical function coords_scaled_eq(a, b)
      type(coords_scaled_t), intent(in) :: a, b
      coords_scaled_eq = &
         (a%Xi == b%Xi) .and. &
         (a%Xe == b%Xe) .and. &
         (a%Yi == b%Yi) .and. &
         (a%Ye == b%Ye) .and. &
         (a%Zi == b%Zi) .and. &
         (a%Ze == b%Ze) .and. &
         (a%xc == b%xc) .and. &
         (a%yc == b%yc) .and. &
         (a%zc == b%zc) .and. &
         (a%Or == b%Or) .and. &
         (a%tag == b%tag)
   end function coords_scaled_eq

   elemental logical function FarField_Sonda_eq(a, b)
      type(FarField_Sonda_t), intent(in) :: a, b
      FarField_Sonda_eq = a%probe == b%probe
   end function FarField_Sonda_eq

   elemental logical function Electric_Sonda_eq(a, b)
      type(Electric_Sonda_t), intent(in) :: a, b
      Electric_Sonda_eq = a%probe == b%probe
   end function Electric_Sonda_eq

   elemental logical function Magnetic_Sonda_eq(a, b)
      type(Magnetic_Sonda_t), intent(in) :: a, b
      Magnetic_Sonda_eq = a%probe == b%probe
   end function Magnetic_Sonda_eq

   elemental logical function NormalElectric_Sonda_eq(a, b) result(res)
      type(NormalElectric_Sonda_t), intent(in) :: a, b
      res = .false.
      if (.not. (a%probe == b%probe)) return
      if (a%n_nml /= b%n_nml) return
      if (a%n_nml_max /= b%n_nml_max) return
      if (associated(a%nml) .eqv. associated(b%nml)) then
         if (associated(a%nml)) then
            if (any(a%nml /= b%nml)) return
         end if
      else if (associated(a%nml)) then
         if (size(a%nml) /= 0) return
      else
         if (size(b%nml) /= 0) return
      end if
      res = .true.
   end function NormalElectric_Sonda_eq

   elemental logical function NormalMagnetic_Sonda_eq(a, b) result(res)
      type(NormalMagnetic_Sonda_t), intent(in) :: a, b
      res = .false.
      if (.not. (a%probe == b%probe)) return
      if (a%n_nml /= b%n_nml) return
      if (a%n_nml_max /= b%n_nml_max) return
      if (associated(a%nml) .eqv. associated(b%nml)) then
         if (associated(a%nml)) then
            if (any(a%nml /= b%nml)) return
         end if
      else if (associated(a%nml)) then
         if (size(a%nml) /= 0) return
      else
         if (size(b%nml) /= 0) return
      end if
      res = .true.
   end function NormalMagnetic_Sonda_eq

   elemental logical function SurfaceElectricCurrent_Sonda_eq(a, b) result(res)
      type(SurfaceElectricCurrent_Sonda_t), intent(in) :: a, b
      res = .false.
      if (.not. (a%probe == b%probe)) return
      if (a%n_nml /= b%n_nml) return
      if (a%n_nml_max /= b%n_nml_max) return
      if (associated(a%nml) .eqv. associated(b%nml)) then
         if (associated(a%nml)) then
            if (any(a%nml /= b%nml)) return
         end if
      else if (associated(a%nml)) then
         if (size(a%nml) /= 0) return
      else
         if (size(b%nml) /= 0) return
      end if
      res = .true.
   end function SurfaceElectricCurrent_Sonda_eq

   elemental logical function SurfaceMagneticCurrent_Sonda_eq(a, b) result(res)
      type(SurfaceMagneticCurrent_Sonda_t), intent(in) :: a, b
      res = .false.
      if (.not. (a%probe == b%probe)) return
      if (a%n_nml /= b%n_nml) return
      if (a%n_nml_max /= b%n_nml_max) return
      if (associated(a%nml) .eqv. associated(b%nml)) then
         if (associated(a%nml)) then
            if (any(a%nml /= b%nml)) return
         end if
      else if (associated(a%nml)) then
         if (size(a%nml) /= 0) return
      else
         if (size(b%nml) /= 0) return
      end if
      res = .true.
   end function SurfaceMagneticCurrent_Sonda_eq

   elemental logical function abstractSonda_eq(a, b) result(res)
      type(abstractSonda_t), intent(in) :: a, b

      res = .false.

      if (a%n_FarField /= b%n_FarField) return
      if (a%n_Electric /= b%n_Electric) return
      if (a%n_Magnetic /= b%n_Magnetic) return
      if (a%n_NormalElectric /= b%n_NormalElectric) return
      if (a%n_NormalMagnetic /= b%n_NormalMagnetic) return
      if (a%n_SurfaceElectricCurrent /= b%n_SurfaceElectricCurrent) return
      if (a%n_SurfaceMagneticCurrent /= b%n_SurfaceMagneticCurrent) return
      if (a%n_FarField_max /= b%n_FarField_max) return
      if (a%n_Electric_max /= b%n_Electric_max) return
      if (a%n_Magnetic_max /= b%n_Magnetic_max) return
      if (a%n_NormalElectric_max /= b%n_NormalElectric_max) return
      if (a%n_NormalMagnetic_max /= b%n_NormalMagnetic_max) return
      if (a%n_SurfaceElectricCurrent_max /= b%n_SurfaceElectricCurrent_max) return
      if (a%n_SurfaceMagneticCurrent_max /= b%n_SurfaceMagneticCurrent_max) return

#define CHECK_PTR(field) \
      if (associated(a%field) .eqv. associated(b%field)) then ; \
         if (associated(a%field)) then ; \
            if (any(.not. a%field == b%field)) return ; \
         end if ; \
      else if (associated(a%field)) then ; \
         if (size(a%field) /= 0) return ; \
      else ; \
         if (size(b%field) /= 0) return ; \
      end if

      CHECK_PTR(FarField)
      CHECK_PTR(Electric)
      CHECK_PTR(Magnetic)
      CHECK_PTR(NormalElectric)
      CHECK_PTR(NormalMagnetic)
      CHECK_PTR(SurfaceElectricCurrent)
      CHECK_PTR(SurfaceMagneticCurrent)

#undef CHECK_PTR

      res = .true.
   end function abstractSonda_eq

   elemental logical function sondas_eq(a, b) result(res)
      type(Sondas_t), intent(in) :: a, b

      res = .false.
      if (a%n_probes /= b%n_probes) return
      if (a%n_probes_max /= b%n_probes_max) return
      if (associated(a%probes) .eqv. associated(b%probes)) then
         if (associated(a%probes)) then
            if (any(.not. a%probes == b%probes)) return
         end if
      else if (associated(a%probes)) then
         if (size(a%probes) /= 0) return
      else
         if (size(b%probes) /= 0) return
      end if
      res = .true.
   end function sondas_eq

   elemental logical function sonda_eq(a, b) result (res)
      type(Sonda_t), intent(in) :: a, b
      res = .false.

      if (a%grname /= b%grname) return

      if (associated(a%i) .eqv. associated(b%i)) then
         if (associated(a%i)) then
            if (any(a%i /= b%i)) return
         end if
      else if (associated(a%i)) then
         if (size(a%i) /= 0) return
      else
         if (size(b%i) /= 0) return
      end if

      if (associated(a%j) .eqv. associated(b%j)) then
         if (associated(a%j)) then
            if (any(a%j /= b%j)) return
         end if
      else if (associated(a%j)) then
         if (size(a%j) /= 0) return
      else
         if (size(b%j) /= 0) return
      end if

      if (associated(a%k) .eqv. associated(b%k)) then
         if (associated(a%k)) then
            if (any(a%k /= b%k)) return
         end if
      else if (associated(a%k)) then
         if (size(a%k) /= 0) return
      else
         if (size(b%k) /= 0) return
      end if

      if (associated(a%node) .eqv. associated(b%node)) then
         if (associated(a%node)) then
            if (any(a%node /= b%node)) return
         end if
      else if (associated(a%node)) then
         if (size(a%node) /= 0) return
      else
         if (size(b%node) /= 0) return
      end if

      if (a%n_cord /= b%n_cord) return
      if (a%n_cord_max /= b%n_cord_max) return
      if (a%tstart /= b%tstart) return
      if (a%tstop /= b%tstop) return
      if (a%tstep /= b%tstep) return
      if (a%outputrequest /= b%outputrequest) return
      if (a%fstart     /= b%fstart) return
      if (a%fstop      /= b%fstop) return
      if (a%fstep      /= b%fstep) return
      if (a%phistart   /= b%phistart) return
      if (a%phistop    /= b%phistop) return
      if (a%phistep    /= b%phistep) return
      if (a%thetastart /= b%thetastart) return
      if (a%thetastop  /= b%thetastop) return
      if (a%thetastep  /= b%thetastep) return
      if (a%FileNormalize /= b%FileNormalize) return
      res = .true.
   end function

   elemental logical function masSonda_eq(a, b) result (res)
      type(MasSonda_t), intent(in) :: a, b
      res = .false.

      if (a%filename /= b%filename) return
      if (a%type1 /= b%type1) return
      if (a%type2 /= b%type2) return
      if (a%outputrequest /= b%outputrequest) return
      if (a%len_cor /= b%len_cor) return
      if (a%tstart /= b%tstart) return
      if (a%tstop /= b%tstop) return
      if (a%tstep /= b%tstep) return
      if (a%fstart /= b%fstart) return
      if (a%fstop /= b%fstop) return
      if (a%fstep /= b%fstep) return

      if (associated(a%cordinates) .eqv. associated(b%cordinates)) then
         if (associated(a%cordinates)) then
            if (any(.not. a%cordinates == b%cordinates)) return
         end if
      else if (associated(a%cordinates)) then
         if (size(a%cordinates) /= 0) return
      else
         if (size(b%cordinates) /= 0) return
      end if
      res = .true.
   end function

   elemental logical function MasSondas_eq(a, b)
      type(MasSondas_t), intent(in) :: a, b

      MasSondas_eq = .false.
      if (a%length /= b%length) return
      if (a%length_max /= b%length_max) return
      if (a%len_cor_max /= b%len_cor_max) return
      if (associated(a%collection) .eqv. associated(b%collection)) then
         if (associated(a%collection)) then
            if (any(.not. a%collection == b%collection)) return
         end if
      else if (associated(a%collection)) then
         if (size(a%collection) /= 0) return
      else
         if (size(b%collection) /= 0) return
      end if
      MasSondas_eq = .true.
   end function MasSondas_eq

   elemental logical function bloqueprobe_eq(a, b)
      type(BloqueProbe_t), intent(in) :: a, b
      bloqueprobe_eq = .false.
      if (a%tstart /= b%tstart) return
      if (a%tstop /= b%tstop) return
      if (a%tstep /= b%tstep) return
      if (a%fstart /= b%fstart) return
      if (a%fstop /= b%fstop) return
      if (a%fstep /= b%fstep) return
      if (trim(a%FileNormalize) /= trim(b%FileNormalize)) return
      if (a%type2 /= b%type2) return
      if (a%i1 /= b%i1) return
      if (a%i2 /= b%i2) return
      if (a%j1 /= b%j1) return
      if (a%j2 /= b%j2) return
      if (a%k1 /= b%k1) return
      if (a%k2 /= b%k2) return
      if (a%skip /= b%skip) return
      if (a%nml /= b%nml) return
      if (a%t .neqv. b%t) return
      if (a%outputrequest /= b%outputrequest) return
      if (a%tag /= b%tag) return
      bloqueprobe_eq = .true.
   end function bloqueprobe_eq

   elemental logical function bloqueprobes_eq(a, b)
      type(BloqueProbes_t), intent(in) :: a, b
      bloqueprobes_eq = .false.
      if (a%n_bp /= b%n_bp) return
      if (a%n_bp_max /= b%n_bp_max) return
      if (associated(a%bp) .eqv. associated(b%bp)) then
         if (associated(a%bp)) then
            if (.not. all(a%bp == b%bp)) return
         end if
      else if (associated(a%bp)) then
         if (size(a%bp) /= 0) return
      else
         if (size(b%bp) /= 0) return
      end if
      bloqueprobes_eq = .true.
   end function bloqueprobes_eq

   elemental logical function volprobe_eq(a, b)
      type(VolProbe_t), intent(in) :: a, b
      volprobe_eq = .false.
      if (a%tstart /= b%tstart) return
      if (a%tstop /= b%tstop) return
      if (a%tstep /= b%tstep) return
      if (a%outputrequest /= b%outputrequest) return
      if (a%len_cor /= b%len_cor) return
      if (a%fstart /= b%fstart) return
      if (a%fstop /= b%fstop) return
      if (a%fstep /= b%fstep) return
      if (a%type2 /= b%type2) return
      if (a%filename /= b%filename) return
      if (associated(a%cordinates) .eqv. associated(b%cordinates)) then
         if (associated(a%cordinates)) then
            if (.not. all(a%cordinates == b%cordinates)) return
         end if
      else if (associated(a%cordinates)) then
         if (size(a%cordinates) /= 0) return
      else
         if (size(b%cordinates) /= 0) return
      end if
      volprobe_eq = .true.
   end function volprobe_eq

   elemental logical function volprobes_eq(a, b)
      type(VolProbes_t), intent(in) :: a, b
      volprobes_eq = .false.
      if (a%length /= b%length) return
      if (a%length_max /= b%length_max) return
      if (a%len_cor_max /= b%len_cor_max) return
      if (associated(a%collection) .eqv. associated(b%collection)) then
         if (associated(a%collection)) then
            if (.not. all(a%collection == b%collection)) return
         end if
      else if (associated(a%collection)) then
         if (size(a%collection) /= 0) return
      else
         if (size(b%collection) /= 0) return
      end if
      volprobes_eq = .true.
   end function volprobes_eq

#endif
end module
