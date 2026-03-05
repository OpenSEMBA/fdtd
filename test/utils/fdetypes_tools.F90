module FDETYPES_TOOLS
   use FDETYPES
   use mod_UTILS
   use NFDETypes
   implicit none
   private

   !===========================
   !  Public interface summary
   !===========================
   public :: observation_domain_t
   public :: initialize_observation_domain_logical_flags
   public :: initialize_observation_time_domain
   public :: initialize_observation_frequency_domain
   public :: initialize_observation_phi_domain
   public :: initialize_observation_theta_domain
   public :: create_observable
   public :: set_observation
   public :: init_time_array
   public :: create_limit_t
   public :: create_xyzlimit_t
   public :: create_xyz_limit_array
   public :: create_tag_list
   public :: create_geometry_media
   public :: create_geometry_media_from_sggAlloc
   public :: create_thinWire_simulation_material
   public :: init_simulation_material_list
   public :: create_facesNF2FF
   public :: create_control_flags
   public :: add_simulation_material
   public :: assing_material_id_to_media_matrix_coordinate
   !===========================

   !===========================
   !  Private interface summary
   !===========================

   !===========================

   real(kind=rkind) :: UTILEPS0 = 8.8541878176203898505365630317107502606083701665994498081024171524053950954599821142852891607182008932e-12
   real(kind=rkind) :: UTILMU0 = 1.2566370614359172953850573533118011536788677597500423283899778369231265625144835994512139301368468271e-6
   type :: observation_domain_t
      real(kind=RKIND) :: InitialTime = 0.0_RKIND
      real(kind=RKIND) :: FinalTime = 0.0_RKIND
      real(kind=RKIND) :: TimeStep = 0.0_RKIND

      real(kind=RKIND) :: InitialFreq = 0.0_RKIND
      real(kind=RKIND) :: FinalFreq = 0.0_RKIND
      real(kind=RKIND) :: FreqStep = 0.0_RKIND

      real(kind=RKIND) :: thetaStart = 0.0_RKIND
      real(kind=RKIND) :: thetaStop = 0.0_RKIND
      real(kind=RKIND) :: thetaStep = 0.0_RKIND

      real(kind=RKIND) :: phiStart = 0.0_RKIND
      real(kind=RKIND) :: phiStop = 0.0_RKIND
      real(kind=RKIND) :: phiStep = 0.0_RKIND

      logical :: FreqDomain = .FALSE.
      logical :: TimeDomain = .FALSE.
      logical :: Saveall = .FALSE.
      logical :: TransFer = .FALSE.
      logical :: Volumic = .FALSE.
   end type observation_domain_t

contains

   function create_xyzlimit_t(XI, XE, YI, YE, ZI, ZE) result(r)
      type(limit_t) :: r
      integer(kind=4), intent(in) :: XI, XE, YI, YE, ZI, ZE
      r%XI = XI
      r%XE = XE
      r%YI = YI
      r%YE = YE
      r%ZI = ZI
      r%ZE = ZE
   end function create_xyzlimit_t

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

   subroutine create_geometry_media(res, xi, xe, yi, ye, zi, ze)
      integer(kind=SINGLE), intent(in) :: xi, xe, yi, ye, zi, ze
      type(media_matrices_t), intent(inout) :: res

      ! Allocate each array with its own kind
      call alloc_and_init(res%sggMtag, xi, xe, yi, ye, zi, ze, 1_IKINDMTAG)
      call alloc_and_init(res%sggMiNo, xi, xe, yi, ye, zi, ze, 1_INTEGERSIZEOFMEDIAMATRICES)
      call alloc_and_init(res%sggMiEx, xi, xe, yi, ye, zi, ze, 1_INTEGERSIZEOFMEDIAMATRICES)
      call alloc_and_init(res%sggMiEy, xi, xe, yi, ye, zi, ze, 1_INTEGERSIZEOFMEDIAMATRICES)
      call alloc_and_init(res%sggMiEz, xi, xe, yi, ye, zi, ze, 1_INTEGERSIZEOFMEDIAMATRICES)
      call alloc_and_init(res%sggMiHx, xi, xe, yi, ye, zi, ze, 1_INTEGERSIZEOFMEDIAMATRICES)
      call alloc_and_init(res%sggMiHy, xi, xe, yi, ye, zi, ze, 1_INTEGERSIZEOFMEDIAMATRICES)
      call alloc_and_init(res%sggMiHz, xi, xe, yi, ye, zi, ze, 1_INTEGERSIZEOFMEDIAMATRICES)
   end subroutine create_geometry_media

   function create_geometry_media_from_sggAlloc(sggAlloc) result(r)
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

      r%sggMtag(:, :, :) = 1
      r%sggMiNo(:, :, :) = 1
      r%sggMiEx(:, :, :) = 1
      r%sggMiEy(:, :, :) = 1
      r%sggMiEz(:, :, :) = 1
      r%sggMiHx(:, :, :) = 1
      r%sggMiHy(:, :, :) = 1
      r%sggMiHz(:, :, :) = 1
   end function create_geometry_media_from_sggAlloc

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
      control%mpidir = 3
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
      if (present(wirecrank)) control%wirecrank = wirecrank
      if (present(resume)) control%resume = resume
      if (present(saveall)) control%saveall = saveall
      if (present(NF2FFDecim)) control%NF2FFDecim = NF2FFDecim
      if (present(simu_devia)) control%simu_devia = simu_devia
      if (present(singlefilewrite)) control%singlefilewrite = singlefilewrite
      if (present(facesNF2FF)) control%facesNF2FF = facesNF2FF

   end function create_control_flags

   subroutine init_simulation_material_list(simulationMaterials)
      implicit none
      type(MediaData_t), dimension(:), allocatable, intent(out) :: simulationMaterials
      if (allocated(simulationMaterials)) deallocate (simulationMaterials)
      allocate (simulationMaterials(0:2))
      simulationMaterials(0) = create_pec_simulation_material()
      simulationMaterials(1) = get_default_mediadata()
      simulationMaterials(2) = create_pmc_simulation_material()

   end subroutine init_simulation_material_list

   subroutine init_time_array(arr, array_size, interval)
      integer, intent(in), optional :: array_size
      real(kind=RKIND_tiempo), intent(in), optional :: interval
      integer(kind=4) :: i
      integer :: size_val
      real(kind=RKIND_tiempo) :: interval_val
      real(kind=RKIND_tiempo), pointer, dimension(:), intent(out) :: arr

      size_val = merge(array_size, 100, present(array_size))
      interval_val = merge(interval, 1.0_RKIND_tiempo, present(interval))

      allocate (arr(size_val))

      DO i = 1, size_val
         arr(i) = (i - 1)*interval_val
      END DO
   end subroutine init_time_array

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

   function create_observable(XI, YI, ZI, XE, YE, ZE, what, line_in) result(observable)
      type(observable_t) :: observable
      integer(kind=4), intent(in)  ::  XI, YI, ZI, XE, YE, ZE, what
      type(direction_t), dimension(:), optional, intent(in) :: line_in

      integer(kind=SINGLE) :: line_size

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

      if (present(line_in)) then
         line_size = size(line_in)

         if (line_size > 0) then
            allocate (observable%line(1:line_size))

            observable%line = line_in
         else

         end if
      end if
   end function create_observable

   subroutine add_simulation_material(simulationMaterials, newSimulationMaterial)
      type(MediaData_t), dimension(:), intent(inout), allocatable :: simulationMaterials
      type(MediaData_t), intent(in) :: newSimulationMaterial

      type(MediaData_t), dimension(:), target, allocatable :: tempSimulationMaterials
      integer(kind=SINGLE) :: oldSize, istat
      oldSize = size(simulationMaterials)
      allocate (tempSimulationMaterials(0:oldSize), stat=istat)
      if (istat /= 0) then
         stop "Allocation failed for temporary media array."
      end if

      if (oldSize > 0) then
         tempSimulationMaterials(0:oldSize - 1) = simulationMaterials
         deallocate (simulationMaterials)
      end if
      tempSimulationMaterials(oldSize) = newSimulationMaterial

      simulationMaterials = tempSimulationMaterials
   end subroutine add_simulation_material

   subroutine add_media_data_to_sgg(sgg, mediaData)
      implicit none

      type(SGGFDTDINFO), intent(inout) :: sgg
      type(MediaData_t), intent(in)    :: mediaData

      type(MediaData_t), dimension(:), target, allocatable :: temp_Med
      integer :: new_size, istat

      new_size = sgg%NumMedia + 1

      allocate (temp_Med(new_size), stat=istat)
      if (istat /= 0) then
         stop "Allocation failed for temporary media array."
      end if

      if (sgg%NumMedia > 0) then
         temp_Med(1:sgg%NumMedia) = sgg%Med

         deallocate (sgg%Med)
      end if

      temp_Med(new_size) = mediaData

      sgg%Med => temp_Med

      sgg%NumMedia = new_size

   end subroutine add_media_data_to_sgg

   subroutine assing_material_id_to_media_matrix_coordinate(media, fieldComponent, i, j, k, materialId)
      type(media_matrices_t), intent(inout) :: media
      integer(kind=SINGLE), intent(in) :: fieldComponent, i, j, k, materialId
      selectcase (fieldComponent)
      case (iEx); media%sggMiEx(i, j, k) = materialId
      case (iEy); media%sggMiEy(i, j, k) = materialId
      case (iEz); media%sggMiEz(i, j, k) = materialId
      case (iHx); media%sggMiHx(i, j, k) = materialId
      case (iHy); media%sggMiHy(i, j, k) = materialId
      case (iHz); media%sggMiHz(i, j, k) = materialId
      end select

   end subroutine assing_material_id_to_media_matrix_coordinate

   function get_default_mediadata() result(res)
      implicit none

      type(MediaData_t) :: res
      !Vacuum id
      res%Id = 0

      ! Reals
      res%Priority = 10
      res%Epr = 1.0_RKIND
      res%Sigma = 0.0_RKIND
      res%Mur = 1.0_RKIND
      res%SigmaM = 0.0_RKIND

      ! Logical
      res%sigmareasignado = .false.

      ! exists_t logicals
      res%Is%PML = .false.
      res%Is%PEC = .false.
      res%Is%PMC = .false.
      res%Is%ThinWire = .false.
      res%Is%SlantedWire = .false.
      res%Is%EDispersive = .false.
      res%Is%MDispersive = .false.
      res%Is%EDispersiveAnis = .false.
      res%Is%MDispersiveAnis = .false.
      res%Is%ThinSlot = .false.
      res%Is%PMLbody = .false.
      res%Is%SGBC = .false.
      res%Is%SGBCDispersive = .false.
      res%Is%Lumped = .false.
      res%Is%Lossy = .false.
      res%Is%AnisMultiport = .false.
      res%Is%Multiport = .false.
      res%Is%MultiportPadding = .false.
      res%Is%Dielectric = .false.
      res%Is%Anisotropic = .false.
      res%Is%Volume = .false.
      res%Is%Line = .false.
      res%Is%Surface = .false.
      res%Is%Needed = .true.
      res%Is%Interfase = .false.
      res%Is%already_YEEadvanced_byconformal = .false.
      res%Is%split_and_useless = .false.

      ! Pointers: They are automatically unassociated (nullified)
      ! when a function returns a type with pointer components,
      ! unless explicitly associated before return.
      ! For safety, we can explicitly nullify them, although Fortran often handles this.
      nullify (res%Wire)
      nullify (res%SlantedWire)
      nullify (res%PMLbody)
      nullify (res%Multiport)
      nullify (res%AnisMultiport)
      nullify (res%EDispersive)
      nullify (res%MDispersive)
      nullify (res%Anisotropic)
      nullify (res%Lumped)

   end function get_default_mediadata

   function create_pec_simulation_material() result(res)
      implicit none

      type(MediaData_t) :: res
      type(Material) :: mat

      mat = create_pec_material()
      res = get_default_mediadata()
      res%Id = mat%id
      res%Is%PEC = .TRUE.

      res%Priority = 150
      res%Epr = mat%eps/UTILEPS0
      res%Sigma = mat%sigma
      res%Mur = mat%mu/UTILMU0
      res%SigmaM = mat%sigmam

   end function create_pec_simulation_material

   function create_pmc_simulation_material() result(res)
      implicit none

      type(MediaData_t) :: res
      type(Material) :: mat

      mat = create_pmc_material()
      res = get_default_mediadata()

      res%Id = mat%id
      res%Is%PMC = .TRUE.

      res%Priority = 160
      res%Epr = mat%eps/UTILEPS0
      res%Sigma = mat%sigma
      res%Mur = mat%mu/UTILMU0
      res%SigmaM = mat%sigmam

   end function create_pmc_simulation_material

   function create_thinWire_simulation_material(materialId) result(res)
      implicit none
      integer(kind=SINGLE) :: materialId

      type(MediaData_t) :: res
      type(Material) :: mat

      type(Wires_t), target, dimension(1) :: wire

      res = get_default_mediadata()
      res%Id = materialId
      res%Is%ThinWire = .TRUE.

      allocate (res%Wire(1))
      wire(1) = get_default_wire()
      res%wire => wire

      res%Priority = 15

   end function create_thinWire_simulation_material

   function create_empty_material() result(mat)
      implicit none
      type(Material) :: mat
   end function create_empty_material

   function create_material(eps_in, mu_in, sigma_in, sigmam_in, id_in) result(mat)
      implicit none
      real(kind=RK), intent(in) :: eps_in, mu_in, sigma_in, sigmam_in
      integer(kind=4), intent(in) :: id_in
      type(Material) :: mat

      mat%eps = eps_in
      mat%mu = mu_in
      mat%sigma = sigma_in
      mat%sigmam = sigmam_in
      mat%id = id_in
   end function create_material

   function create_vacuum_material() result(mat)
      type(Material) :: mat
      mat = create_material(EPSILON_VACUUM, MU_VACUUM, 0.0_RKIND, 0.0_RKIND, 1)
   end function create_vacuum_material

   function create_pec_material() result(mat)
      type(Material) :: mat
      mat = create_material(EPSILON_VACUUM, MU_VACUUM, SIGMA_PEC, 0.0_RKIND, 0)
   end function create_pec_material

   function create_pmc_material() result(mat)
      type(Material) :: mat
      mat = create_material(EPSILON_VACUUM, MU_VACUUM, 0.0_RKIND, SIGMA_PMC, 2)
   end function create_pmc_material

   function create_empty_materials() result(mats)
      implicit none
      type(Materials) :: mats
   end function create_empty_materials

   subroutine add_material_to_materials(mats_collection, new_mat)
      implicit none
      type(Materials), intent(inout) :: mats_collection
      type(Material), intent(in) :: new_mat

      type(Material), dimension(:), target, allocatable :: temp_Mats
      integer :: old_size, new_size

      old_size = mats_collection%n_Mats
      new_size = old_size + 1

      allocate (temp_Mats(new_size))

      if (old_size > 0) then
         temp_Mats(1:old_size) = mats_collection%Mats

         deallocate (mats_collection%Mats)
      end if

      temp_Mats(new_size) = new_mat

      mats_collection%Mats => temp_Mats

      mats_collection%n_Mats = new_size
      mats_collection%n_Mats_max = new_size

   end subroutine add_material_to_materials

   function get_default_wire() result(wire)
      implicit none
      type(Wires_t) :: wire

      wire%Radius = 0.0_RKIND_wires
      wire%R = 0.0_RKIND_wires
      wire%L = 0.0_RKIND_wires
      wire%C = 0.0_RKIND_wires
      wire%P_R = 0.0_RKIND_wires
      wire%P_L = 0.0_RKIND_wires
      wire%P_C = 0.0_RKIND_wires
      wire%Radius_devia = 0.0_RKIND_wires
      wire%R_devia = 0.0_RKIND_wires
      wire%L_devia = 0.0_RKIND_wires
      wire%C_devia = 0.0_RKIND_wires

      wire%numsegmentos = 0
      wire%NUMVOLTAGESOURCES = 0
      wire%NUMCURRENTSOURCES = 0

      nullify (wire%segm)
      nullify (wire%Vsource)
      nullify (wire%Isource)

      wire%VsourceExists = .false.
      wire%IsourceExists = .false.
      wire%HasParallel_LeftEnd = .false.
      wire%HasParallel_RightEnd = .false.
      wire%HasSeries_LeftEnd = .false.
      wire%HasSeries_RightEnd = .false.
      wire%HasAbsorbing_LeftEnd = .false.
      wire%HasAbsorbing_RightEnd = .false.

      wire%Parallel_R_RightEnd = 0.0_RKIND_wires
      wire%Parallel_R_LeftEnd = 0.0_RKIND_wires
      wire%Series_R_RightEnd = 0.0_RKIND_wires
      wire%Series_R_LeftEnd = 0.0_RKIND_wires
      wire%Parallel_L_RightEnd = 0.0_RKIND_wires
      wire%Parallel_L_LeftEnd = 0.0_RKIND_wires
      wire%Series_L_RightEnd = 0.0_RKIND_wires
      wire%Series_L_LeftEnd = 0.0_RKIND_wires
      wire%Parallel_C_RightEnd = 0.0_RKIND_wires
      wire%Parallel_C_LeftEnd = 0.0_RKIND_wires
      wire%Series_C_RightEnd = 2.0e7_RKIND ! Valor por defecto de corto
      wire%Series_C_LeftEnd = 2.0e7_RKIND ! Valor por defecto de corto

      wire%Parallel_R_RightEnd_devia = 0.0_RKIND_wires
      wire%Parallel_R_LeftEnd_devia = 0.0_RKIND_wires
      wire%Series_R_RightEnd_devia = 0.0_RKIND_wires
      wire%Series_R_LeftEnd_devia = 0.0_RKIND_wires
      wire%Parallel_L_RightEnd_devia = 0.0_RKIND_wires
      wire%Parallel_L_LeftEnd_devia = 0.0_RKIND_wires
      wire%Series_L_RightEnd_devia = 0.0_RKIND_wires
      wire%Series_L_LeftEnd_devia = 0.0_RKIND_wires
      wire%Parallel_C_RightEnd_devia = 0.0_RKIND_wires
      wire%Parallel_C_LeftEnd_devia = 0.0_RKIND_wires
      wire%Series_C_RightEnd_devia = 0.0_RKIND_wires
      wire%Series_C_LeftEnd_devia = 0.0_RKIND_wires

      wire%LeftEnd = 0
      wire%RightEnd = 0
   end function get_default_wire

   subroutine set_observation(obs, P_in, outputrequest_in, domain_params, FileNormalize_in)
      implicit none

      type(observable_t), dimension(:), intent(in) :: P_in
      character(LEN=*), intent(in) :: outputrequest_in, FileNormalize_in
      type(observation_domain_t), intent(in) :: domain_params

      type(Obses_t), intent(out) :: obs
      integer(kind=4) :: n_count

      n_count = size(P_in)
      obs%nP = n_count

      allocate (obs%P(1:n_count))

      obs%P(1:n_count) = P_in(1:n_count)

      obs%outputrequest = outputrequest_in
      obs%FileNormalize = FileNormalize_in

      obs%InitialTime = domain_params%InitialTime
      obs%FinalTime = domain_params%FinalTime
      obs%TimeStep = domain_params%TimeStep

      obs%InitialFreq = domain_params%InitialFreq
      obs%FinalFreq = domain_params%FinalFreq
      obs%FreqStep = domain_params%FreqStep

      obs%thetaStart = domain_params%thetaStart
      obs%thetaStop = domain_params%thetaStop
      obs%thetaStep = domain_params%thetaStep

      obs%phiStart = domain_params%phiStart
      obs%phiStop = domain_params%phiStop
      obs%phiStep = domain_params%phiStep

      obs%FreqDomain = domain_params%FreqDomain
      obs%TimeDomain = domain_params%TimeDomain
      obs%Saveall = domain_params%Saveall
      obs%TransFer = domain_params%TransFer
      obs%Volumic = domain_params%Volumic

   end subroutine set_observation

   subroutine initialize_observation_time_domain(domain, InitialTime, FinalTime, TimeStep)
      implicit none

      type(observation_domain_t), intent(inout) :: domain
      real(kind=RKIND), intent(in) :: InitialTime, FinalTime, TimeStep

      domain%InitialTime = InitialTime
      domain%FinalTime = FinalTime
      domain%TimeStep = TimeStep

      domain%TimeDomain = .true.

   end subroutine initialize_observation_time_domain

   subroutine initialize_observation_frequency_domain(domain, InitialFreq, FinalFreq, FreqStep)
      implicit none

      type(observation_domain_t), intent(inout) :: domain
      real(kind=RKIND), intent(in) :: InitialFreq, FinalFreq, FreqStep

      domain%InitialFreq = InitialFreq
      domain%FinalFreq = FinalFreq
      domain%FreqStep = FreqStep

      domain%FreqDomain = .true.

   end subroutine initialize_observation_frequency_domain

   subroutine initialize_observation_theta_domain(domain, thetaStart, thetaStop, thetaStep)
      implicit none

      type(observation_domain_t), intent(inout) :: domain
      real(kind=RKIND), intent(in) :: thetaStart, thetaStop, thetaStep

      domain%thetaStart = thetaStart
      domain%thetaStop = thetaStop
      domain%thetaStep = thetaStep

   end subroutine initialize_observation_theta_domain

   subroutine initialize_observation_phi_domain(domain, phiStart, phiStop, phiStep)
      implicit none

      type(observation_domain_t), intent(inout) :: domain
      real(kind=RKIND), intent(in) :: phiStart, phiStop, phiStep

      domain%phiStart = phiStart
      domain%phiStop = phiStop
      domain%phiStep = phiStep

   end subroutine initialize_observation_phi_domain

   subroutine initialize_observation_domain_logical_flags(domain, Saveall_flag, TransFer_flag, Volumic_flag)
      implicit none

      type(observation_domain_t), intent(inout) :: domain
      logical, intent(in) :: Saveall_flag, TransFer_flag, Volumic_flag

      domain%Saveall = Saveall_flag
      domain%TransFer = TransFer_flag
      domain%Volumic = Volumic_flag

   end subroutine initialize_observation_domain_logical_flags

end module FDETYPES_TOOLS
