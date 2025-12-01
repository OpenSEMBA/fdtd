module FDETYPES_TOOLS
   use FDETYPES
contains
   function create_limit_t(XI, XE, YI, YE, ZI, ZE, NX, NY, NZ) result(r)
      type(limit_t) :: r
      integer(kind=4), intent(in) :: XI, XE, YI, YE, ZI, ZE, NX, NY, NZ
      r%XI = XI
      r%XE = XE
      r%YI = YI
      r%YE = YE
      r%ZI = ZI
      r%ZE = ZE
      r%NX = NX
      r%NY = NY
      r%NZ = NZ
   end function create_limit_t
   function create_tag_list(sggAlloc) result(r)
      type(XYZlimit_t), dimension(6), intent(in) :: sggAlloc
      type(taglist_t) :: r

      allocate (r%edge%x(sggAlloc(iEx)%XI:sggAlloc(iEx)%XE, sggAlloc(iEx)%YI:sggAlloc(iEx)%YE, sggAlloc(iEx)%ZI:sggAlloc(iEx)%ZE))
      allocate (r%edge%y(sggAlloc(iEy)%XI:sggAlloc(iEy)%XE, sggAlloc(iEy)%YI:sggAlloc(iEy)%YE, sggAlloc(iEy)%ZI:sggAlloc(iEy)%ZE))
      allocate (r%edge%z(sggAlloc(iEz)%XI:sggAlloc(iEz)%XE, sggAlloc(iEz)%YI:sggAlloc(iEz)%YE, sggAlloc(iEz)%ZI:sggAlloc(iEz)%ZE))
      allocate (r%face%x(sggAlloc(iHx)%XI:sggAlloc(iHx)%XE, sggAlloc(iHx)%YI:sggAlloc(iHx)%YE, sggAlloc(iHx)%ZI:sggAlloc(iHx)%ZE))
      allocate (r%face%y(sggAlloc(iHy)%XI:sggAlloc(iHy)%XE, sggAlloc(iHy)%YI:sggAlloc(iHy)%YE, sggAlloc(iHy)%ZI:sggAlloc(iHy)%ZE))
      allocate (r%face%z(sggAlloc(iHz)%XI:sggAlloc(iHz)%XE, sggAlloc(iHz)%YI:sggAlloc(iHz)%YE, sggAlloc(iHz)%ZI:sggAlloc(iHz)%ZE))

      r%edge%x(:, :, :) = 0
      r%edge%y(:, :, :) = 0
      r%edge%z(:, :, :) = 0
      r%face%x(:, :, :) = 0
      r%face%y(:, :, :) = 0
      r%face%z(:, :, :) = 0
   end function create_tag_list

   function create_media(sggAlloc) result(r)
      type(XYZlimit_t), dimension(6), intent(in) :: sggAlloc
      type(media_matrices_t) :: r

      allocate (r%sggMtag(sggAlloc(iHx)%XI:sggAlloc(iHx)%XE, sggAlloc(iHy)%YI:sggAlloc(iHy)%YE, sggAlloc(iHz)%ZI:sggAlloc(iHz)%ZE))
      allocate (r%sggMiNo(sggAlloc(iHx)%XI:sggAlloc(iHx)%XE, sggAlloc(iHy)%YI:sggAlloc(iHy)%YE, sggAlloc(iHz)%ZI:sggAlloc(iHz)%ZE))

      allocate (r%sggMiEx(sggAlloc(iEx)%XI:sggAlloc(iEx)%XE, sggAlloc(iEx)%YI:sggAlloc(iEx)%YE, sggAlloc(iEx)%ZI:sggAlloc(iEx)%ZE))
      allocate (r%sggMiEy(sggAlloc(iEy)%XI:sggAlloc(iEy)%XE, sggAlloc(iEy)%YI:sggAlloc(iEy)%YE, sggAlloc(iEy)%ZI:sggAlloc(iEy)%ZE))
      allocate (r%sggMiEz(sggAlloc(iEz)%XI:sggAlloc(iEz)%XE, sggAlloc(iEz)%YI:sggAlloc(iEz)%YE, sggAlloc(iEz)%ZI:sggAlloc(iEz)%ZE))

      allocate (r%sggMiHx(sggAlloc(iHx)%XI:sggAlloc(iHx)%XE, sggAlloc(iHx)%YI:sggAlloc(iHx)%YE, sggAlloc(iHx)%ZI:sggAlloc(iHx)%ZE))
      allocate (r%sggMiHy(sggAlloc(iHy)%XI:sggAlloc(iHy)%XE, sggAlloc(iHy)%YI:sggAlloc(iHy)%YE, sggAlloc(iHy)%ZI:sggAlloc(iHy)%ZE))
      allocate (r%sggMiHz(sggAlloc(iHz)%XI:sggAlloc(iHz)%XE, sggAlloc(iHz)%YI:sggAlloc(iHz)%YE, sggAlloc(iHz)%ZI:sggAlloc(iHz)%ZE))

      r%sggMtag(:, :, :) = 0
      r%sggMiNo(:, :, :) = 1
      r%sggMiEx(:, :, :) = 1
      r%sggMiEy(:, :, :) = 1
      r%sggMiEz(:, :, :) = 1
      r%sggMiHx(:, :, :) = 1
      r%sggMiHy(:, :, :) = 1
      r%sggMiHz(:, :, :) = 1
   end function create_media

   function create_control_flags(layoutnumber, size, mpidir, finaltimestep, &
                                 nEntradaRoot, wiresflavor, wirecrank, &
                                 resume, saveall, NF2FFDecim, simu_devia, singlefilewrite, &
                                 facesNF2FF) result(control)

      type(sim_control_t) :: control

      integer(kind=SINGLE), intent(in), optional :: layoutnumber, size, mpidir, finaltimestep
      character(len=*), intent(in), optional :: nEntradaRoot, wiresflavor
      logical, intent(in), optional :: wirecrank, resume, saveall, NF2FFDecim, simu_devia, singlefilewrite
      type(nf2ff_t), intent(in), optional :: facesNF2FF

      ! 1. Set explicit defaults for all components
      control%layoutnumber = 0
      control%size = 0
      control%mpidir = 0
      control%finaltimestep = 0
      control%nEntradaRoot = ""
      control%wiresflavor = ""
      control%wirecrank = .false.
      control%resume = .false.
      control%saveall = .false.
      control%NF2FFDecim = .false.
      control%simu_devia = .false.
      control%singlefilewrite = .false.
      ! Note: control%facesNF2FF retains its default initialized state

      ! 2. Overwrite defaults only if the optional argument is present
      if (present(layoutnumber)) control%layoutnumber = layoutnumber
      if (present(size)) control%size = size
      if (present(mpidir)) control%mpidir = mpidir
      if (present(finaltimestep)) control%finaltimestep = finaltimestep
      if (present(nEntradaRoot)) control%nEntradaRoot = nEntradaRoot
      if (present(wiresflavor)) control%wiresflavor = wiresflavor
      if (present(wiresflavor)) control%wirecrank = wirecrank
      if (present(resume)) control%resume = resume
      if (present(saveall)) control%saveall = saveall
      if (present(NF2FFDecim)) control%NF2FFDecim = NF2FFDecim
      if (present(simu_devia)) control%simu_devia = simu_devia
      if (present(singlefilewrite)) control%singlefilewrite = singlefilewrite
      if (present(facesNF2FF)) control%facesNF2FF = facesNF2FF

   end function create_control_flags

   function create_base_sgg(NumMedia, dt, time_steps) result(sgg)
      type(SGGFDTDINFO) :: sgg
      integer, optional, intent(in) :: NumMedia, time_steps
      real(kind=RKIND_tiempo), optional, intent(in) :: dt

      sgg%NumMedia = merge(NumMedia, 3, present(NumMedia))
      allocate (sgg%Med(1:sgg%NumMedia))
      sgg%Med = create_basic_media()
      sgg%NumberRequest = 1
      sgg%dt = merge(dt, 0.1_RKIND_tiempo, present(dt))

      ! Use the new optional-aware create_time_array
      sgg%tiempo = create_time_array(merge(time_steps, 100, present(time_steps)), sgg%dt)

      ! Hardcoded array limits now call the optional-aware function
      sgg%Sweep = create_xyz_limit_array(0, 0, 0, 6, 6, 6)
      sgg%SINPMLSweep = create_xyz_limit_array(1, 1, 1, 5, 5, 5)
      sgg%NumPlaneWaves = 1
      sgg%alloc = create_xyz_limit_array(0, 0, 0, 6, 6, 6)

   end function create_base_sgg

   function create_time_array(array_size, interval) result(arr)
      integer, intent(in), optional :: array_size
      real(kind=RKIND_tiempo), intent(in), optional :: interval
      integer(kind=4) :: i
      integer :: size_val
      real(kind=RKIND_tiempo) :: interval_val
      real(kind=RKIND_tiempo), allocatable, dimension(:) :: arr

      size_val = merge(array_size, 100, present(array_size))
      interval_val = merge(interval, 1.0_RKIND_tiempo, present(interval))

      

      allocate (arr(size_val))

      DO i = 1, size_val
         arr(i) = (i - 1)*interval_val
      END DO
   end function create_time_array

   function create_limit_type() result(r)
      type(limit_t) :: r
   end function create_limit_type

   function create_xyz_limit_array(XI, YI, ZI, XE, YE, ZE) result(arr)
      type(XYZlimit_t), dimension(1:6) :: arr
      integer(kind=4), intent(in), optional :: XI, YI, ZI, XE, YE, ZE
      integer :: i
      integer(kind=4) :: xi_val, yi_val, zi_val, xe_val, ye_val, ze_val

      ! Use merge for compact handling of optional inputs with defaults
      xi_val = merge(XI, 0, present(XI))
      yi_val = merge(YI, 0, present(YI))
      zi_val = merge(ZI, 0, present(ZI))
      xe_val = merge(XE, 6, present(XE))
      ye_val = merge(YE, 6, present(YE))
      ze_val = merge(ZE, 6, present(ZE))

      do i = 1, 6
         arr(i)%XI = xi_val
         arr(i)%XE = xe_val
         arr(i)%YI = yi_val
         arr(i)%YE = ye_val
         arr(i)%ZI = zi_val
         arr(i)%ZE = ze_val
      end do
   end function create_xyz_limit_array

   function create_facesNF2FF(tr, fr, iz, de, ab, ar) result(faces)
      type(nf2ff_t) :: faces
      logical, intent(in), optional :: tr, fr, iz, de, ab, ar

      faces%tr = .false.
      faces%fr = .false.
      faces%iz = .false.
      faces%de = .false.
      faces%ab = .false.
      faces%ar = .false.

      if (present(tr)) faces%tr = tr
      if (present(fr)) faces%fr = fr
      if (present(iz)) faces%iz = iz
      if (present(de)) faces%de = de
      if (present(ab)) faces%ab = ab
      if (present(ar)) faces%ar = ar
   end function create_facesNF2FF

   function create_basic_media() result(media)
      type(MediaData_t) :: media
   end function create_basic_media

   function define_point_observation() result(obs)
      type(Obses_t) :: obs

      obs%nP = 1
      allocate (obs%P(obs%nP))
      obs%P(1) = create_observable(1_SINGLE, 1_SINGLE, 1_SINGLE, 1_SINGLE, 1_SINGLE, 1_SINGLE, iEx)

      obs%InitialTime = 0.0_RKIND_tiempo
      obs%FinalTime = 1.0_RKIND_tiempo
      obs%TimeStep = 0.1_RKIND_tiempo

      obs%InitialFreq = 0.0_RKIND
      obs%FinalFreq = 0.0_RKIND
      obs%FreqStep = 0.0_RKIND

      obs%outputrequest = 'pointProbe'

      obs%FreqDomain = .false.
      obs%TimeDomain = .true.
      obs%Saveall = .false.
      obs%TransFer = .false.
      obs%Volumic = .false.
      obs%Done = .false.
      obs%Begun = .false.
      obs%Flushed = .false.

   end function define_point_observation

   function define_wire_current_observation() result(obs)
          type(Obses_t) :: obs

      obs%nP = 1
      allocate (obs%P(obs%nP))
      obs%P(1) = create_observable(3_SINGLE, 3_SINGLE, 3_SINGLE, 3_SINGLE, 3_SINGLE, 3_SINGLE, iJx)

      obs%InitialTime = 0.0_RKIND_tiempo
      obs%FinalTime = 1.0_RKIND_tiempo
      obs%TimeStep = 0.1_RKIND_tiempo

      obs%InitialFreq = 0.0_RKIND
      obs%FinalFreq = 0.0_RKIND
      obs%FreqStep = 0.0_RKIND

      obs%outputrequest = 'pointProbe'

      obs%FreqDomain = .false.
      obs%TimeDomain = .true.
      obs%Saveall = .false.
      obs%TransFer = .false.
      obs%Volumic = .false.
      obs%Done = .false.
      obs%Begun = .false.
      obs%Flushed = .false.
   end function define_wire_current_observation

   
   function define_wire_charge_observation() result(obs)
          type(Obses_t) :: obs

      obs%nP = 1
      allocate (obs%P(obs%nP))
      obs%P(1) = create_observable(3_SINGLE, 3_SINGLE, 3_SINGLE, 3_SINGLE, 3_SINGLE, 3_SINGLE, iQx)

      obs%InitialTime = 0.0_RKIND_tiempo
      obs%FinalTime = 1.0_RKIND_tiempo
      obs%TimeStep = 0.1_RKIND_tiempo

      obs%InitialFreq = 0.0_RKIND
      obs%FinalFreq = 0.0_RKIND
      obs%FreqStep = 0.0_RKIND

      obs%outputrequest = 'pointProbe'

      obs%FreqDomain = .false.
      obs%TimeDomain = .true.
      obs%Saveall = .false.
      obs%TransFer = .false.
      obs%Volumic = .false.
      obs%Done = .false.
      obs%Begun = .false.
      obs%Flushed = .false.
   end function define_wire_charge_observation

   function create_observable(XI,YI,ZI,XE,YE,ZE, what) result(observable)
                type(observable_t) :: observable
            integer (kind=4)  ::  XI,YI,ZI,XE,YE,ZE, what

            observable%XI = XI
            observable%YI = YI
            observable%ZI = ZI

            observable%XE = XE
            observable%YE = YE
            observable%ZE = ZE

            observable%Xtrancos = 1
            observable%Ytrancos = 1
            observable%Ztrancos = 1

            observable%What = what
    end function create_observable

end module FDETYPES_TOOLS
