! Device Holland thickness-1 wire step. Included only with CompileWithCUDA.

logical function cuda_wires_is_ready()
   cuda_wires_is_ready = cuda_wires_ready
end function

subroutine InitWires_cuda(control, b)
   type(sim_control_t), intent(in) :: control
   type(bounds_t), intent(in) :: b
   type(fdtd_wire_seg_c), allocatable :: segs(:)
   type(fdtd_wire_node_c), allocatable :: nodes(:)
   real(c_float), allocatable :: samples(:), current(:), charge(:), charge_past(:)
   integer :: nseg, nnode, ns, irc, s, nd

   cuda_wires_ready = .false.
   if (.not. fdtd_cuda_ok_f()) return
   if (HWires%NumCurrentSegments < 1) return
   if (control%wirethickness /= 1) return
   if (control%wirecrank .or. control%stochastic .or. control%simu_devia) return
   if (trim(adjustl(control%wiresflavor)) /= 'holland') return

   nseg = HWires%NumCurrentSegments
   nnode = HWires%NumChargeNodes
   allocate(segs(nseg), nodes(max(nnode, 1)))
   allocate(current(nseg), charge(max(nnode, 1)), charge_past(max(nnode, 1)))
   allocate(samples(1))
   samples(1) = 0.0_c_float
   current = 0.0_c_float
   charge = 0.0_c_float
   charge_past = 0.0_c_float
   ns = 0

   do nd = 1, nnode
      nodes(nd)%exists = 0
      nodes(nd)%is_mur = 0
      nodes(nd)%is_periodic = 0
      nodes(nd)%n_plus = 0
      nodes(nd)%n_minus = 0
      nodes(nd)%plus_seg = -1
      nodes(nd)%minus_seg = -1
      nodes(nd)%node_inside = -1
      nodes(nd)%has_i = 0
      nodes(nd)%evol_off = 0
      nodes(nd)%numus = 0
      nodes(nd)%deltaevol = 1.0_c_float
      nodes(nd)%cte_prop = 0.0_c_float
      nodes(nd)%cte_plain = 0.0_c_float
      nodes(nd)%cte_mur = 0.0_c_float
      if (.not. HWires%ChargeNode(nd)%exists) cycle
      nodes(nd)%exists = 1
      if (HWires%ChargeNode(nd)%IsMur) nodes(nd)%is_mur = 1
      if (HWires%ChargeNode(nd)%IsPeriodic) nodes(nd)%is_periodic = 1
      nodes(nd)%n_plus = HWires%ChargeNode(nd)%NumCurrentPlus
      nodes(nd)%n_minus = HWires%ChargeNode(nd)%NumCurrentMinus
      if (nodes(nd)%n_plus > 9 .or. nodes(nd)%n_minus > 9) return
      if (.not. push_plus(nd, 1, HWires%ChargeNode(nd)%CurrentPlus_1)) return
      if (.not. push_plus(nd, 2, HWires%ChargeNode(nd)%CurrentPlus_2)) return
      if (.not. push_plus(nd, 3, HWires%ChargeNode(nd)%CurrentPlus_3)) return
      if (.not. push_plus(nd, 4, HWires%ChargeNode(nd)%CurrentPlus_4)) return
      if (.not. push_plus(nd, 5, HWires%ChargeNode(nd)%CurrentPlus_5)) return
      if (.not. push_plus(nd, 6, HWires%ChargeNode(nd)%CurrentPlus_6)) return
      if (.not. push_plus(nd, 7, HWires%ChargeNode(nd)%CurrentPlus_7)) return
      if (.not. push_plus(nd, 8, HWires%ChargeNode(nd)%CurrentPlus_8)) return
      if (.not. push_plus(nd, 9, HWires%ChargeNode(nd)%CurrentPlus_9)) return
      if (.not. push_minus(nd, 1, HWires%ChargeNode(nd)%CurrentMinus_1)) return
      if (.not. push_minus(nd, 2, HWires%ChargeNode(nd)%CurrentMinus_2)) return
      if (.not. push_minus(nd, 3, HWires%ChargeNode(nd)%CurrentMinus_3)) return
      if (.not. push_minus(nd, 4, HWires%ChargeNode(nd)%CurrentMinus_4)) return
      if (.not. push_minus(nd, 5, HWires%ChargeNode(nd)%CurrentMinus_5)) return
      if (.not. push_minus(nd, 6, HWires%ChargeNode(nd)%CurrentMinus_6)) return
      if (.not. push_minus(nd, 7, HWires%ChargeNode(nd)%CurrentMinus_7)) return
      if (.not. push_minus(nd, 8, HWires%ChargeNode(nd)%CurrentMinus_8)) return
      if (.not. push_minus(nd, 9, HWires%ChargeNode(nd)%CurrentMinus_9)) return
      if (HWires%ChargeNode(nd)%IsMur) then
         nodes(nd)%node_inside = node_of(HWires%ChargeNode(nd)%NodeInside)
         if (nodes(nd)%node_inside < 0) return
         nodes(nd)%cte_mur = real(HWires%ChargeNode(nd)%cteMur, c_float)
      end if
      nodes(nd)%cte_prop = real(HWires%ChargeNode(nd)%cteprop, c_float)
      nodes(nd)%cte_plain = real(HWires%ChargeNode(nd)%ctePlain, c_float)
      charge(nd) = real(HWires%ChargeNode(nd)%ChargePresent, c_float)
      charge_past(nd) = real(HWires%ChargeNode(nd)%ChargePast, c_float)
      if (HWires%ChargeNode(nd)%HasIsource) then
         if (.not. associated(HWires%ChargeNode(nd)%Isource)) return
         nodes(nd)%has_i = 1
         call append_samples(HWires%ChargeNode(nd)%Isource%Fichero%Samples, &
            HWires%ChargeNode(nd)%Isource%Fichero%NumSamples, &
            HWires%ChargeNode(nd)%Isource%Fichero%DeltaSamples, &
            nodes(nd)%evol_off, nodes(nd)%numus, nodes(nd)%deltaevol)
      end if
   end do

   do s = 1, nseg
      segs(s)%coupled = 0
      segs(s)%field_comp = 0
      segs(s)%i = HWires%CurrentSegment(s)%i
      segs(s)%j = HWires%CurrentSegment(s)%j
      segs(s)%k = HWires%CurrentSegment(s)%k
      segs(s)%is_pmc = 0
      if (HWires%CurrentSegment(s)%IsPMC) segs(s)%is_pmc = 1
      segs(s)%charge_plus = node_of(HWires%CurrentSegment(s)%ChargePlus)
      segs(s)%charge_minus = node_of(HWires%CurrentSegment(s)%ChargeMinus)
      if (segs(s)%charge_plus < 0 .or. segs(s)%charge_minus < 0) return
      segs(s)%has_v = 0
      segs(s)%evol_off = 0
      segs(s)%numus = 0
      segs(s)%deltaevol = 1.0_c_float
      segs(s)%cte1 = real(HWires%CurrentSegment(s)%cte1, c_float)
      segs(s)%cte2 = real(HWires%CurrentSegment(s)%cte2, c_float)
      segs(s)%cte3 = real(HWires%CurrentSegment(s)%cte3, c_float)
      segs(s)%cte5 = real(HWires%CurrentSegment(s)%cte5, c_float)
      segs(s)%fraction_plus = real(HWires%CurrentSegment(s)%FractionPlus, c_float)
      segs(s)%fraction_minus = real(HWires%CurrentSegment(s)%FractionMinus, c_float)
      segs(s)%vscale = 0.0_c_float
      select case (HWires%CurrentSegment(s)%tipofield)
      case (iEx)
         segs(s)%field_comp = 0
      case (iEy)
         segs(s)%field_comp = 1
      case (iEz)
         segs(s)%field_comp = 2
      case default
         return
      end select
      if ((.not. HWires%CurrentSegment(s)%IsShielded) .and. &
          associated(HWires%CurrentSegment(s)%Efield_wire2main) .and. &
          .not. associated(HWires%CurrentSegment(s)%Efield_wire2main, HWires%null_field)) then
         segs(s)%coupled = 1
      end if
      current(s) = real(HWires%CurrentSegment(s)%Current, c_float)
      if (HWires%CurrentSegment(s)%HasVsource) then
         if (.not. associated(HWires%CurrentSegment(s)%Vsource)) return
         if (HWires%CurrentSegment(s)%Lind == 0.0_RKIND_wires) return
         segs(s)%has_v = 1
         segs(s)%vscale = real(HWires%CurrentSegment(s)%cte3 / (HWires%CurrentSegment(s)%Lind * &
            InvMu(HWires%CurrentSegment(s)%indexmed) * InvEps(HWires%CurrentSegment(s)%indexmed)), c_float)
         call append_samples(HWires%CurrentSegment(s)%Vsource%Fichero%Samples, &
            HWires%CurrentSegment(s)%Vsource%Fichero%NumSamples, &
            HWires%CurrentSegment(s)%Vsource%Fichero%DeltaSamples, &
            segs(s)%evol_off, segs(s)%numus, segs(s)%deltaevol)
      end if
   end do

   if (ns < 1) ns = 1
   irc = fdtd_cuda_upload_wires_f(segs, nseg, nodes, nnode, samples, ns, &
      current, charge, charge_past, &
      b%Ex%XI, b%Ex%YI, b%Ex%ZI, b%Ey%XI, b%Ey%YI, b%Ey%ZI, b%Ez%XI, b%Ez%YI, b%Ez%ZI)
   if (irc == 0) return
   cuda_wires_ready = .true.

contains
   integer function seg_of(ptr)
      type(CurrentSegments_t), pointer :: ptr
      integer :: n
      seg_of = -1
      if (.not. associated(ptr)) return
      do n = 1, HWires%NumCurrentSegments
         if (associated(ptr, HWires%CurrentSegment(n))) then
            seg_of = n - 1
            return
         end if
      end do
      seg_of = -2
   end function

   integer function node_of(ptr)
      type(ChargeNodes_t), pointer :: ptr
      integer :: n
      node_of = -1
      if (.not. associated(ptr)) return
      do n = 1, HWires%NumChargeNodes
         if (associated(ptr, HWires%ChargeNode(n))) then
            node_of = n - 1
            return
         end if
      end do
      node_of = -2
   end function

   logical function push_plus(nd, slot, ptr)
      integer, intent(in) :: nd, slot
      type(CurrentSegments_t), pointer :: ptr
      integer :: idx
      push_plus = .true.
      if (slot > HWires%ChargeNode(nd)%NumCurrentPlus) return
      idx = seg_of(ptr)
      if (idx < 0) then
         push_plus = .false.
         return
      end if
      nodes(nd)%plus_seg(slot) = idx
   end function

   logical function push_minus(nd, slot, ptr)
      integer, intent(in) :: nd, slot
      type(CurrentSegments_t), pointer :: ptr
      integer :: idx
      push_minus = .true.
      if (slot > HWires%ChargeNode(nd)%NumCurrentMinus) return
      idx = seg_of(ptr)
      if (idx < 0) then
         push_minus = .false.
         return
      end if
      nodes(nd)%minus_seg(slot) = idx
   end function

   subroutine append_samples(src, numus, delta, off, out_numus, out_delta)
      real(kind=RKIND_wires), pointer :: src(:)
      integer, intent(in) :: numus
      real(kind=RKIND_wires), intent(in) :: delta
      integer, intent(out) :: off, out_numus
      real(c_float), intent(out) :: out_delta
      real(c_float), allocatable :: grown(:)
      integer :: ncopy, i, lb
      off = ns
      out_numus = numus
      out_delta = real(delta, c_float)
      if (.not. associated(src)) return
      ncopy = numus + 1
      lb = lbound(src, 1)
      if (ns + ncopy > size(samples)) then
         allocate(grown(max(ns + ncopy, 2 * size(samples))))
         if (ns > 0) grown(1:ns) = samples(1:ns)
         call move_alloc(grown, samples)
      end if
      do i = 1, ncopy
         if (lb + i - 1 <= ubound(src, 1)) then
            samples(ns + i) = real(src(lb + i - 1), c_float)
         else
            samples(ns + i) = 0.0_c_float
         end if
      end do
      ns = ns + ncopy
   end subroutine
end subroutine InitWires_cuda

subroutine AdvanceWiresE_cuda(sgg, timeinstant)
   type(SGGFDTDINFO_t), intent(in) :: sgg
   integer, intent(in) :: timeinstant
   real(c_float), allocatable :: current(:), current_past(:), charge(:), charge_past(:)
   integer :: nseg, nnode, s, nd, irc
   real(kind=RKIND) :: time_i, time_q

   if (.not. cuda_wires_ready) return
   nseg = HWires%NumCurrentSegments
   nnode = HWires%NumChargeNodes
   allocate(current(max(nseg, 1)), current_past(max(nseg, 1)))
   allocate(charge(max(nnode, 1)), charge_past(max(nnode, 1)))
   time_i = sgg%tiempo(timeinstant)
   time_q = time_i - unmedio * sgg%dt
   irc = fdtd_cuda_advance_wires_e_f(time_i, time_q)
   if (irc == 0) call stoponerror(0, 0, 'CUDA Holland wire advance failed')
   irc = fdtd_cuda_download_wires_f(current, current_past, charge, charge_past)
   if (irc == 0) call stoponerror(0, 0, 'CUDA Holland wire download failed')
   do s = 1, nseg
      HWires%CurrentSegment(s)%CurrentPast = real(current_past(s), RKIND_wires)
      HWires%CurrentSegment(s)%Current = real(current(s), RKIND_wires)
   end do
   do nd = 1, nnode
      if (.not. HWires%ChargeNode(nd)%exists) cycle
      HWires%ChargeNode(nd)%ChargePast = real(charge_past(nd), RKIND_wires)
      HWires%ChargeNode(nd)%ChargePresent = real(charge(nd), RKIND_wires)
   end do
end subroutine AdvanceWiresE_cuda
