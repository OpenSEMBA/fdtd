module smbjson

#ifdef CompileWithSMBJSON
   use NFDETypes

   use NFDETypes_extension
   use smbjson_labels_mod
   use mesh_mod
   use parser_tools_mod
   use idchildtable_mod

   use json_module
   use json_kinds


   use, intrinsic :: iso_fortran_env , only: error_unit

   implicit none

   integer, private, parameter  :: MAX_LINE = BUFSIZE
   character (len=*), parameter :: TAG_MATERIAL = 'material'
   character (len=*), parameter :: TAG_LAYER = 'layer'

   type, public :: parser_t
      private
      character (len=:), allocatable :: filename
      type(json_file), pointer :: jsonfile => null()
      type(json_core), pointer :: core => null()
      type(json_value), pointer :: root => null()
      type(mesh_t) :: mesh
      type(IdChildTable_t) :: matTable, elementTable
      logical, public :: isInitialized = .false.

   contains
      procedure :: readProblemDescription
      procedure :: readMesh
   
#ifdef CompileWithMTLN
      procedure :: readMTLN
#endif

      ! private
      procedure, private :: readGeneral
      procedure, private :: readGrid
      procedure, private :: readMediaMatrix
      procedure, private :: readPECRegions
      procedure, private :: readPMCRegions
      procedure, private :: readDielectricRegions
      procedure, private :: readLossyThinSurfaces
      procedure, private :: readBoundary
      procedure, private :: readPlanewaves
      procedure, private :: readNodalSources
      procedure, private :: readProbes
      procedure, private :: readMoreProbes
      procedure, private :: readBlockProbes
      procedure, private :: readVolumicProbes
      procedure, private :: readThinWires
      procedure, private :: readThinSlots
      procedure, private :: readConformalRegions
      !
      !
      procedure, private :: getLogicalAt
      procedure, private :: getIntAt
      procedure, private :: getIntsAt
      procedure, private :: getRealAt
      procedure, private :: getRealsAt
      procedure, private :: getMatrixAt
      procedure, private :: getStrAt
      procedure, private :: existsAt
      procedure, private :: getDomain
      procedure, private :: buildPECPMCRegions
      procedure, private :: getMaterialAssociations
      procedure, private :: parseMaterialAssociation
      procedure, private :: matAssToCoords
      procedure, private :: buildTagName
      procedure, private :: jsonValueFilterByKeyValue
      procedure, private :: jsonValueFilterByKeyValues
      procedure, private :: getSingleVolumeInElementsIds
   end type
   interface parser_t
      module procedure parser_ctor
   end interface

   type, private :: thinwiretermination_t
      integer :: terminationType
      real :: r, l, c
   end type

   type, private :: generator_description_t
      character(:), allocatable :: srctype, srcfile
      real :: multiplier
   end type

   type, private :: materialAssociation_t
      ! Common fields
      character(:), allocatable :: name
      integer :: materialId
      integer, dimension(:), allocatable :: elementIds
      character(:), allocatable :: matAssType
      ! Cable specific fields.
      integer :: initialTerminalId = -1
      integer :: endTerminalId = -1
      integer :: initialConnectorId = -1
      integer :: endConnectorId = -1
      integer :: containedWithinElementId = -1
   end type

   type, private :: domain_t
      real :: tstart = 0.0, tstop = 0.0, tstep = 0.0
      real :: fstart = 0.0, fstop = 0.0
      integer :: fstep = 0
      character(len=:), allocatable :: filename
      integer :: type1 = NP_T1_PLAIN, type2 = NP_T2_TIME
      logical :: isLogarithmicFrequencySpacing = .false.
   end type
contains
   function parser_ctor(filename) result(res)
      type(parser_t) :: res
      character(len=*), intent(in) :: filename
      res%filename = filename
      
      allocate(res%jsonfile)
      call res%jsonfile%initialize()
      if (res%jsonfile%failed()) then
         call res%jsonfile%print_error_message(error_unit)
         return
      end if

      call res%jsonfile%load(filename = res%filename)
      if (res%jsonfile%failed()) then
         call res%jsonfile%print_error_message(error_unit)
         return
      end if

      allocate(res%core)
      call res%jsonfile%get_core(res%core)
      call res%jsonfile%get('.', res%root)

      res%isInitialized = .true.
   end function

   function readProblemDescription(this) result (res)
      class(parser_t) :: this
      type(Parseador) :: res
      integer :: stat 

      this%mesh = this%readMesh()
      this%matTable = IdChildTable_t(this%core, this%root, J_MATERIALS)
      this%elementTable = IdChildTable_t(this%core, this%root, J_MESH//'.'//J_ELEMENTS)
      
      call initializeProblemDescription(res)
      
      ! Basics
      res%general = this%readGeneral()
      res%matriz = this%readMediaMatrix()
      res%despl = this%readGrid()
      res%front = this%readBoundary()
      
      ! Materials
      res%pecRegs = this%readPECRegions()
      res%pmcRegs = this%readPMCRegions()
      res%dielRegs = this%readDielectricRegions()
      res%lossyThinSurfs = this%readLossyThinSurfaces()
      
      ! Sources
      res%plnSrc = this%readPlanewaves()
      res%nodSrc = this%readNodalSources()
      
      ! Probes
      res%oldSonda = this%readProbes()
      res%sonda = this%readMoreProbes()
      res%BloquePrb = this%readBlockProbes()
      res%VolPrb = this%readVolumicProbes()
      
      ! Thin elements
      res%tWires = this%readThinWires()
      res%tSlots = this%readThinSlots()

#ifdef CompileWithMTLN
      write(*,*) 'debug: mtl'
      res%mtln = this%readMTLN(res%despl)
#endif

      res%conformalRegs = this%readConformalRegions()

   end function

   function readMesh(this) result(res)
      class(parser_t) :: this
      type(Mesh_t) :: res
      call addCoordinates(res)
      call addElements(res)

   contains
      subroutine addCoordinates(mesh)
          type(mesh_t), intent(inout) :: mesh
          type(json_value), pointer :: jcs, jc
          integer :: id, i
          real, dimension(:), allocatable :: pos
          type(coordinate_t) :: c
          integer :: numberOfCoordinates 
          logical :: found
          
          call this%core%get(this%root, J_MESH//'.'//J_COORDINATES, jcs, found=found)
          if (found) then
             numberOfCoordinates = this%core%count(jcs)
             call res%allocateCoordinates(10*numberOfCoordinates)
             do i = 1, numberOfCoordinates
                call this%core%get_child(jcs, i, jc)
                id = this%getIntAt(jc, J_ID)
                pos = this%getRealsAt(jc, J_COORDINATE_POS)
                c%position = pos
                call mesh%addCoordinate(id, c)
             end do
          end if
      end subroutine
   
      subroutine addElements(mesh)
         type(mesh_t), intent(inout) :: mesh
         character (len=:), allocatable :: elementType
         type(json_value), pointer :: jes, je
         integer :: id, i
         type(node_t) :: node
         type(polyline_t) :: polyline
         integer, dimension(:), allocatable :: coordIds
         integer :: numberOfElements
         logical :: found
         
         call this%core%get(this%root, J_MESH//'.'//J_ELEMENTS, jes, found=found)
         numberOfElements = this%core%count(jes)
         call res%allocateElements(10*numberOfElements)
             
         if (found) then
            do i = 1, numberOfElements
               call this%core%get_child(jes, i, je)
               id = this%getIntAt(je, J_ID)
               elementType = this%getStrAt(je, J_TYPE)
               select case (elementType)
                case (J_ELEM_TYPE_NODE)
                  coordIds = this%getIntsAt(je, J_COORDINATE_IDS)
                  node%coordIds = coordIds
                  call mesh%addElement(id, node)
                case (J_ELEM_TYPE_POLYLINE)
                  coordIds = this%getIntsAt(je, J_COORDINATE_IDS)
                  polyline%coordIds = coordIds
                  call mesh%addElement(id, polyline)
                CASE (J_ELEM_TYPE_CELL)
                  block
                     type(cell_region_t) :: cR
                     type(cell_interval_t), dimension(:), allocatable :: intervals
                     cR%intervals = readCellIntervals(je, J_CELL_INTERVALS)
                     call mesh%addCellRegion(id, cR)
                     
                  end block
                case (J_ELEM_TYPE_CONF_VOLUME) 
                  block 
                     type(conformal_region_t) :: cV
                     cV%intervals = readCellIntervals(je, J_CELL_INTERVALS)
                     cV%triangles = readTriangles(je, J_CONF_VOLUME_TRIANGLES)
                     cV%offgrid_points = this%getIntAt(je,J_CONF_OFFGRID)
                     call mesh%addConformalRegion(id, cV)
                  end block
                case default
                  write (error_unit, *) 'Invalid element type'
               end select
            end do
         end if
      end subroutine

      function readCellIntervals(place, path) result(res)
         type(json_value), pointer, intent(in) :: place
         character (len=*), intent(in) :: path
         type(cell_interval_t), dimension(:), allocatable :: res

         type(json_value), pointer :: intervalsPlace, interval
         integer :: i, nIntervals
         real, dimension(:), allocatable :: cellIni, cellEnd
         logical :: containsInterval

         call this%core%get(place, path, intervalsPlace, found=containsInterval)
         if (.not. containsInterval) then
            allocate(res(0))
            return
         end if
         nIntervals = this%core%count(intervalsPlace)
         allocate(res(nIntervals))
         do i = 1, nIntervals
            call this%core%get_child(intervalsPlace, i, interval)
            cellIni = this%getRealsAt(interval, '(1)')
            cellEnd = this%getRealsAt(interval, '(2)')
            res(i)%ini%cell = cellIni(1:3)
            res(i)%end%cell = cellEnd(1:3)
         end do
      end function

      function readTriangles(place, path) result(res)
         type(json_value), pointer, intent(in) :: place
         character (len=*), intent(in) :: path
         type(triangle_t), dimension(:), allocatable :: res

         type(json_value), pointer :: triangles, triangle_ptr
         real, dimension(:), allocatable :: triangle
         integer :: i, j, nTriangles
         real, dimension(:), allocatable :: cellIni, cellEnd

         logical :: containsTriangles
         call this%core%get(place, path, triangles, found=containsTriangles)
         if (.not. containsTriangles) then
            allocate(res(0))
            return
         end if
         nTriangles = this%core%count(triangles)
         allocate(res(nTriangles))
         do i = 1, nTriangles
            call this%core%get_child(triangles, i, triangle_ptr)
            call this%core%get(triangle_ptr, triangle)
            do j = 1, 3
               res(i)%vertices(j)%id = triangle(j)
            end do
         end do

      end function

   end function

   function readGeneral(this) result (res)
      class(parser_t) :: this
      type(NFDEGeneral) :: res
      res%dt = this%getRealAt(this%root, J_GENERAL//'.'//J_GEN_TIME_STEP)
      res%nmax = this%getRealAt(this%root, J_GENERAL//'.'//J_GEN_NUMBER_OF_STEPS)
      res%mtlnProblem = this%getLogicalAt(this%root, J_GENERAL//'.'//J_GEN_MTLN_PROBLEM, default = .false.)
   end function

   function readMediaMatrix(this) result(res)
      class(parser_t) :: this
      type(MatrizMedios) :: res
      character (len=*), parameter :: P = J_MESH//'.'//J_GRID//'.'//J_GRID_NUMBER_OF_CELLS
      res%totalX = this%getIntAt(this%root, P//'(1)')
      res%totalY = this%getIntAt(this%root, P//'(2)')
      res%totalZ = this%getIntAt(this%root, P//'(3)')
   end function

   function readGrid(this) result (res)
      class(parser_t) :: this
      type(Desplazamiento) :: res
      real, dimension(:), allocatable :: vec
      character (len=*), parameter :: P = J_MESH//'.'//J_GRID

      res%nX = this%getIntAt(this%root, P//'.'//J_GRID_NUMBER_OF_CELLS//'(1)')
      res%nY = this%getIntAt(this%root, P//'.'//J_GRID_NUMBER_OF_CELLS//'(2)')
      res%nZ = this%getIntAt(this%root, P//'.'//J_GRID_NUMBER_OF_CELLS//'(3)')

      call assignDes(P//'.'//J_GRID_STEPS//'.x', res%desX, res%nX)
      call assignDes(P//'.'//J_GRID_STEPS//'.y', res%desY, res%nY)
      call assignDes(P//'.'//J_GRID_STEPS//'.z', res%desZ, res%nZ)

      res%mx1 = 0
      res%my1 = 0
      res%mz1 = 0
      res%mx2 = res%nX
      res%my2 = res%nY
      res%mz2 = res%nZ


   contains
      subroutine assignDes(path, dest, numberOfCells)
         character(kind=CK, len=*) :: path
         real (kind=rkind), dimension(:), pointer :: dest
         real, dimension(:), allocatable :: vec
         integer (kind=4), intent(in) :: numberOfCells
         logical :: found = .false.

         vec= this%getRealsAt(this%root, path, found)

         if (.not. found) then
            write(error_unit, *) 'Error reading grid: steps not found.'
         endif
         if (size(vec) /= 1 .and. size(vec) /= numberOfCells) then
            write(error_unit, *) &
               'Error reading grid: steps must be arrays of size 1 (for regular grids) or size equal to the number of cells.'
         end if

         allocate(dest(0 : numberOfCells-1))
         if (size(vec) == 1) then
            dest(:) = vec(1)
         else
            dest = vec
         end if
      end subroutine
   end function

   function readBoundary(this) result (res)
      class(parser_t) :: this
      type(Frontera) :: res
      character (len=:), allocatable :: bdrType
      type(json_value), pointer :: bdrs
      logical :: found
      character(len=*), parameter :: errorMsgInit = "ERROR reading boundary: "
      
      call this%core%get(this%root, J_BOUNDARY, bdrs, found)
      if (.not. found) then
         write(error_unit, * ) errorMsgInit, J_BOUNDARY, " object not found."
      end if
      
      block
         bdrType = this%getStrAt(bdrs, J_BND_ALL//'.'//J_TYPE, found)
         if (found) then
            res%tipoFrontera(:) = labelToBoundaryType(bdrType)
            if (all(res%tipoFrontera == F_PML)) then
               res%propiedadesPML(:) = readPMLProperties(J_BOUNDARY//"."//J_BND_ALL)
            end if
            return
         end if
      end block
         
      block
         character(len=*), dimension(6), parameter :: placeLabels = &
            [J_BND_XL, J_BND_XU, J_BND_YL, J_BND_YU, J_BND_ZL, J_BND_ZU]
         integer :: i, j
         do i = 1, 6
            bdrType = this%getStrAt(bdrs, placeLabels(i)//"."//J_TYPE, found)
            if (.not. found) then
               write(error_unit, *) errorMsgInit, placeLabels(i), " or ", J_BND_ALL, " not found."
            end if
            j = labelToBoundaryPlace(placeLabels(i))
            res%tipoFrontera(j) = labelToBoundaryType(bdrType)
            if (res%tipoFrontera(j) == F_PML) then
               res%propiedadesPML(j) = readPMLProperties(J_BOUNDARY//"."//placeLabels(i))
            end if
         end do
      end block

   contains
      function readPMLProperties(p) result(res)
         type(FronteraPML) :: res
         character(len=*), intent(in) :: p
         res%numCapas = this%getIntAt(this%root, p//'.'//J_BND_PML_LAYERS, default=8)
         res%orden = this%getRealAt(this%root, p//'.'//J_BND_PML_ORDER, default=2.0)
         res%refl = this%getRealAt(this%root, p//'.'//J_BND_PML_REFLECTION, default=0.001)
      end function

      function labelToBoundaryPlace(str) result (place)
         character(len=*), intent(in) :: str
         integer :: place
         select case (str)
            case (J_BND_XL)
               place = F_XL
            case (J_BND_XU)
               place = F_XU
            case (J_BND_YL)
               place = F_YL
            case (J_BND_YU)
               place = F_YU
            case (J_BND_ZL)
               place = F_ZL
            case (J_BND_ZU)
               place = F_ZU
         end select
      end function

      function labelToBoundaryType(str) result (bdrType)
         character(len=:), allocatable, intent(in) :: str
         integer :: bdrType
         select case (str)
          case (J_BND_TYPE_PEC)
            bdrType = F_PEC
          case (J_BND_TYPE_PMC)
            bdrType = F_PMC
          case (J_BND_TYPE_PERIODIC)
            bdrType = F_PER
          case (J_BND_TYPE_MUR)
            bdrType = F_MUR
          case (J_BND_TYPE_PML)
            bdrType = F_PML
         end select
      end function
   end function

   function readPECRegions(this) result (res)
      class(parser_t), intent(in) :: this
      type(PECRegions) :: res
      res = this%buildPECPMCRegions(J_MAT_TYPE_PEC)
   end function

   function readPMCRegions(this) result (res)
      class(parser_t), intent(in) :: this
      type(PECRegions) :: res
      res = this%buildPECPMCRegions(J_MAT_TYPE_PMC)
   end function

   function buildPECPMCRegions(this, matType) result(res)
      class(parser_t) :: this
      character (len=*), intent(in) :: matType 
      type(PECRegions) :: res
      type(materialAssociation_t), dimension(:), allocatable :: mAs
      type(coords), dimension(:), pointer :: cs
      integer :: i
      
      mAs = this%getMaterialAssociations([matType])
      
      block
         type(coords), dimension(:), pointer :: emptyCoords
         if (size(mAs) == 0) then 
            allocate(emptyCoords(0))
            call appendRegion(res%lins,  res%nLins,  res%nLins_max,  emptyCoords)
            call appendRegion(res%surfs, res%nSurfs, res%nSurfs_max, emptyCoords)
            call appendRegion(res%vols,  res%nVols,  res%nVols_max,  emptyCoords)
            return
         end if
      end block

      do i = 1, size(mAs)
         call this%matAssToCoords(cs, mAs(i), CELL_TYPE_LINEL)
         call appendRegion(res%lins,  res%nLins,  res%nLins_max, cs)
         call this%matAssToCoords(cs, mAs(i), CELL_TYPE_SURFEL)
         call appendRegion(res%surfs, res%nSurfs, res%nSurfs_max, cs)
         call this%matAssToCoords(cs, mAs(i), CELL_TYPE_VOXEL)
         call appendRegion(res%vols,  res%nVols,  res%nVols_max, cs)
         deallocate(cs)
      end do

   contains
      subroutine appendRegion(resCoords, resNCoords, resNCoordsMax, cs)
         type(coords), dimension(:), pointer :: resCoords
         integer, intent(out) :: resNCoords, resNCoordsMax
         type(coords), dimension(:), pointer, intent(in) :: cs
         type(coords) , dimension(:), allocatable :: auxCs
         integer :: i

         if (.not. associated(resCoords)) then
            allocate(resCoords(size(cs)))
            do i = 1, size(cs) ! Use do loop to prevent stack overflow in large cs
                resCoords(i) = cs(i)
            end do
            resNCoords = size(cs)
            resNCoordsMax = size(cs)
         else 
            allocate(auxCs(size(resCoords)))
            do i = 1, size(resCoords)
                auxCs(i) = resCoords(i)
            end do
            deallocate(resCoords)
            
            allocate(resCoords(size(auxCs) + size(cs)))
            resNCoords = size(resCoords)
            resNCoordsMax = size(resCoords)
            do i = 1, size(auxCs)
                resCoords(i) = auxCs(i)
            end do
            do i = 1, size(cs)
                resCoords(i + size(auxCs)) = cs(i)
            end do
         end if
      end subroutine         
   end function

   function readConformalRegions(this) result(res)
      class(parser_t) :: this
      type(ConformalPECRegions) :: res
      type(materialAssociation_t), dimension(:), allocatable :: mAs
      type(conformal_region_t) :: cR
      integer :: i, j, k
      type(cell_map_t) :: cell_map

      mAs = this%getMaterialAssociations([J_MAT_TYPE_PEC])
      if (size(mAs) == 0) then 
         do i = 1, size(mAs)
            do j = 1, size(mAs(i)%elementIds)
               cR = this%mesh%getConformalRegion(mAs(i)%elementIds(j))
               cell_map = buildCellTriMap(cR%triangles)
            end do
         end do
      end if

      ! buildEdgeRegions(cell_map)
      ! buildFaceRegions(cell_map)


   end function

   function buildCellTriMap(triangles) result(res)
      type(triangle_t), dimension(:), allocatable :: triangles
      type(cell_map_t) :: res
      integer :: i, stat
      integer (kind=4), dimension(3) :: cell
      do i = 1, size(triangles)
         call res%addTriangle(triangles(i))
      end do
   end function

   function buildEdgeRegions(cell_map) result(res)
      type(cell_map_t) :: cell_map
      type(triangle_t), dimension(:), allocatable :: triangles
      type(side_t), dimension(:), allocatable :: sides_X, sides_Y, sides_Z
      type(side_t), dimension(3) :: sides
      logical :: res
      integer :: i, j, k
      do i = 1, size(cell_map%keys) 
         triangles = getTrianglesOffFaces(cell_map%getTrianglesInCell(cell_map%keys(i)%cell))
         sides_X = getSidesOnFace(triangles, 1) !FACE_X
         sides_Y = getSidesOnFace(triangles, 2) !FACE_Y
         sides_Z = getSidesOnFace(triangles, 3) !FACE_Z
         ! do j = 1, size(triangles)
         !    sides = triangles(j)%getSides()
         !    do k = 1, 3
         !       side
         !    end do
         ! end do
      end do

   end function

   function getTrianglesOffFaces(triangles) result (res)
      type(triangle_t), dimension(:), allocatable, intent(in) :: triangles
      type(triangle_t), dimension(:), allocatable :: res
      integer :: i, j, n
      n = 0
      do i = 1, size(triangles)
         if (.not. (triangles(i)%isOnFace(1) .or. triangles(i)%isOnFace(1) .or. triangles(i)%isOnFace(1))) n = n + 1
      end do
      allocate(res(n))
      n = 0
      do i = 1, size(triangles)
         if (.not. (triangles(i)%isOnFace(1) .or. triangles(i)%isOnFace(1) .or. triangles(i)%isOnFace(1))) then 
            n = n + 1
            res(n) = triangles(j)
         end if
      end do

   end function   

   function getSidesOnFace(triangles, face) result (res)
      type(triangle_t), dimension(:), allocatable, intent(in) :: triangles
      type(side_t), dimension(3) :: sides
      integer, intent(in) :: face
      type(side_t), dimension(:), allocatable :: res
      integer :: i, j, n
      n = 0
      do i = 1, size(triangles)
         sides = triangles(i)%getSides()
         do j = 1, 3
            if (sides(j)%isOnFace(face) .and. all(sides(j)%getCell() .eq. triangles(i)%getCell())) n = n + 1
         end do
      end do
      allocate(res(n))
      n = 0
      do i = 1, size(triangles)
         sides = triangles(i)%getSides()
         do j = 1, 3
            if (sides(j)%isOnFace(face) .and. all(sides(j)%getCell() .eq. triangles(i)%getCell())) then 
               n = n + 1
               res(n) = sides(j)
            end if
         end do
      end do

   end function

   function readDielectricRegions(this) result (res)
      class(parser_t), intent(in) :: this
      type(DielectricRegions) :: res
      
      call fillDielectricsOfCellType(res%vols, CELL_TYPE_VOXEL)
      call fillDielectricsOfCellType(res%surfs, CELL_TYPE_SURFEL)
      call fillDielectricsOfCellType(res%lins, CELL_TYPE_LINEL)
      
      res%nVols = size(res%vols)
      res%nSurfs = size(res%Surfs)
      res%nLins = size(res%Lins)

      res%nVols_max = res%nVols
      res%nSurfs_max = res%nSurfs
      res%nLins_max = res%nLins
   contains
      subroutine fillDielectricsOfCellType(res, cellType)
         integer, intent(in) :: cellType
         type(dielectric_t), dimension(:), pointer :: res
         
         type(materialAssociation_t), dimension(:), allocatable :: mAs
         type(materialAssociation_t) :: mA
         type(cell_region_t) :: cR

         integer :: i, j
         integer :: nCs, nDielectrics
         
         mAs = this%getMaterialAssociations([J_MAT_TYPE_ISOTROPIC])
         if (size(mAs) == 0) then
            allocate(res(0))
            return
         end if

         ! Precounts
         nDielectrics = 0
         do i = 1, size(mAs)           
            if (containsCellRegionsWithType(mAs(i), cellType)) then
               nDielectrics = nDielectrics + 1
            end if 
         end do

         ! Fills
         allocate(res(nDielectrics))
         
         if (nDielectrics == 0) return

         j = 0
         do i = 1, size(mAs)       
            if (.not. containsCellRegionsWithType(mAs(i), cellType)) cycle
            j = j + 1
            res(j) = readDielectric(mAs(i), cellType)
         end do
      end subroutine

      function readDielectric(mA, cellType) result(res)
         type(materialAssociation_t), intent(in) :: mA
         integer, intent(in) :: cellType
         type(Dielectric_t) :: res
         type(cell_region_t) :: cR
         type (coords), dimension(:), allocatable :: coords
         type(json_value_ptr) :: matPtr
         integer :: e, j

         allocate(res%c1P(0))
         res%n_c1p = 0
         call this%matAssToCoords(res%c2p, mA, cellType)
         res%n_c2p = size(res%c2p)
         
         matPtr = this%matTable%getId(mA%materialId)
         ! Fills rest of dielectric data.
         res%sigma  = this%getRealAt(matPtr%p, J_MAT_ELECTRIC_CONDUCTIVITY, default=0.0)
         res%sigmam = this%getRealAt(matPtr%p, J_MAT_MAGNETIC_CONDUCTIVITY, default=0.0)
         res%eps    = this%getRealAt(matPtr%p, J_MAT_REL_PERMITTIVITY, default=1.0)*EPSILON_VACUUM
         res%mu     = this%getRealAt(matPtr%p, J_MAT_REL_PERMEABILITY, default=1.0)*MU_VACUUM

      end function

      logical function containsCellRegionsWithType(mA, cellType)
         integer, intent(in) :: cellType
         type(materialAssociation_t), intent(in) :: mA
         integer :: e
         type(cell_region_t) :: cR
         
         do e = 1, size(mA%elementIds)
            cR = this%mesh%getCellRegion(mA%elementIds(e))
            if (size(cellRegionToCoords(cR, cellType)) /= 0) then
               containsCellRegionsWithType = .true.
               return
            end if
         end do

         containsCellRegionsWithType = .false.
      end function
   end function

   subroutine matAssToCoords(this, res, mA, cellType)
      class(parser_t) :: this
      type(materialAssociation_t), intent(in) :: mA
      type (coords), dimension(:), pointer :: res
      integer, intent(in) :: cellType
      character (len=:), allocatable :: tagName
      type (coords), dimension(:), allocatable :: newCoords
      type (cell_region_t) :: cR
      integer :: nCs
      integer :: e, jIni, jEnd
      
      ! Precount
      nCs = 0
      do e = 1, size(mA%elementIds)
         cR = this%mesh%getCellRegion(mA%elementIds(e))
         newCoords = cellRegionToCoords(cR, cellType)
         nCs = nCs + size(newCoords)
      end do

      ! Fills coords
      jIni = 1
      allocate(res(nCs))
      do e = 1, size(mA%elementIds)
         cR = this%mesh%getCellRegion(mA%elementIds(e))
         tagName = this%buildTagName(mA%materialId, mA%elementIds(e))
         newCoords = cellRegionToCoords(cR, cellType, tag=tagName)
         if (size(newCoords) == 0) cycle
         jEnd = jIni + size(newCoords) - 1
         res(jIni:jEnd) = newCoords(:)
         jIni = jEnd + 1 
      end do
   end subroutine

   function readLossyThinSurfaces(this) result (res)
      class(parser_t), intent(in) :: this
      type(LossyThinSurfaces) :: res
      type(materialAssociation_t), dimension(:), allocatable :: mAs
      type(json_value_ptr) :: mat
      integer :: nLossySurfaces
      logical :: found
      integer :: i, j, k
      type(coords), dimension(:), pointer :: cs

      mAs = this%getMaterialAssociations([J_MAT_TYPE_MULTILAYERED_SURFACE])
      
      ! Precounts
      nLossySurfaces = 0
      do i = 1, size(mAs)
         call this%matAssToCoords(cs, mAs(i), CELL_TYPE_SURFEL)
         if (size(cs) > 0) nLossySurfaces = nLossySurfaces + 1
      end do

      ! Fills
      if (nLossySurfaces == 0) then
         res = emptyLossyThinSurfaces()
         return
      end if

      allocate(res%cs(nLossySurfaces))
      res%length = nLossySurfaces
      res%length_max = nLossySurfaces
      res%nC_max = nLossySurfaces
      k = 1
      do i = 1, size(mAs)
         call this%matAssToCoords(cs, mAs(i), CELL_TYPE_SURFEL)
         if (size(cs) == 0) cycle
         res%cs(k) = readLossyThinSurface(mAs(i))
         k = k + 1
      end do
      
   contains
      function readLossyThinSurface(mA) result(res)
         type(materialAssociation_t), intent(in) :: mA
         type(LossyThinSurface) :: res
         logical :: found
         character (len=*), parameter :: errorMsgInit = "ERROR reading lossy thin surface: "
         integer :: i
         type(json_value_ptr) :: mat
         type(json_value), pointer :: layer
         type(json_value), pointer :: layers
         
         call this%matAssToCoords(res%c, mA, CELL_TYPE_SURFEL)
         res%nc = size(res%c)

         mat = this%matTable%getId(mA%materialId)
         call this%core%get(mat%p, J_MAT_MULTILAYERED_SURF_LAYERS, layers)

         res%numcapas = this%core%count(layers)
         allocate(res%sigma( res%numcapas))
         allocate(res%eps(   res%numcapas))
         allocate(res%mu(    res%numcapas))
         allocate(res%sigmam(res%numcapas))
         allocate(res%thk(   res%numcapas))
         allocate(res%sigma_devia( res%numcapas))
         allocate(res%eps_devia(   res%numcapas))
         allocate(res%mu_devia(    res%numcapas))
         allocate(res%sigmam_devia(res%numcapas))
         allocate(res%thk_devia(   res%numcapas))
         do i = 1, res%numcapas
            call this%core%get_child(layers, i, layer)
            res%sigma(i)  = this%getRealAt(layer, J_MAT_ELECTRIC_CONDUCTIVITY, default=0.0)
            res%sigmam(i) = this%getRealAt(layer, J_MAT_MAGNETIC_CONDUCTIVITY, default=0.0)
            res%eps(i)    = this%getRealAt(layer, J_MAT_REL_PERMITTIVITY, default=1.0) * EPSILON_VACUUM
            res%mu(i)     = this%getRealAt(layer, J_MAT_REL_PERMEABILITY, default=1.0) * MU_VACUUM
            res%thk(i)    = this%getRealAt(layer, J_MAT_MULTILAYERED_SURF_THICKNESS, found)
            if (.not. found) then
               write(error_unit, *) errorMsgInit, J_MAT_MULTILAYERED_SURF_THICKNESS, " in layer not found."
            end if
            res%sigma_devia(i) = 0.0
            res%eps_devia(i) = 0.0
            res%mu_devia(i) = 0.0
            res%sigmam_devia(i) = 0.0
            res%thk_devia(i) = 0.0
         end do
      end function

      function emptyLossyThinSurfaces() result (res)
         type(LossyThinSurfaces) :: res
         allocate(res%cs(0))
         res%length = 0
         res%length_max = 0
         res%nC_max = 0
      end function
   end function

   function readPlanewaves(this) result (res)
      class(parser_t) :: this
      type(PlaneWaves) :: res
      type(json_value), pointer :: sources
      type(json_value_ptr), allocatable :: pws(:)
      integer :: i
      logical :: found

      call this%core%get(this%root, J_SOURCES, sources, found)

      if (.not. found) then
         allocate(res%collection(0))
         res%nc = size(res%collection)
         res%nc_max = size(res%collection)
         return
      end if

      pws = this%jsonValueFilterByKeyValue(sources, J_TYPE, J_SRC_TYPE_PW)

      allocate(res%collection(size(pws)))
      do i=1, size(pws)
         res%collection(i) = readPlanewave(pws(i)%p)
      end do
      res%nc = size(res%collection)
      res%nc_max = size(res%collection)

   contains
      function readPlanewave(pw) result (res)
         type(PlaneWave) :: res
         type(json_value), pointer :: pw

         character (len=:), allocatable :: label
         logical :: found

         res%nombre_fichero = trim(adjustl(this%getStrAt(pw,J_SRC_MAGNITUDE_FILE)))

         res%atributo = ""

         res%theta = this%getRealAt(pw, J_SRC_PW_DIRECTION//'.'//J_SRC_PW_THETA)
         res%phi = this%getRealAt(pw, J_SRC_PW_DIRECTION//'.'//J_SRC_PW_PHI)
         res%alpha = this%getRealAt(pw, J_SRC_PW_POLARIZATION//'.'//J_SRC_PW_THETA)
         res%beta = this%getRealAt(pw, J_SRC_PW_POLARIZATION//'.'//J_SRC_PW_PHI)

         block
            type(coords), dimension(:), allocatable :: nfdeCoords
            nfdeCoords = &
               cellIntervalsToCoords(this%getSingleVolumeInElementsIds(pw))
            res%coor1 = [nfdeCoords(1)%Xi, nfdeCoords(1)%Yi, nfdeCoords(1)%Zi]
            res%coor2 = [nfdeCoords(1)%Xe, nfdeCoords(1)%Ye, nfdeCoords(1)%Ze]
         end block

         res%isRC = .false.
         res%nummodes = 1
         res%incertmax = 0.0
      end function
   end function

   function readNodalSources(this) result (res)
      type(NodSource) :: res
      class(parser_t) :: this
      type(json_value), pointer :: sources
      type(json_value_ptr), dimension(:), allocatable :: nodSrcs
      logical :: found
      integer :: i

      call this%core%get(this%root, J_SOURCES, sources, found)
      if (.not. found) then
         allocate(res%NodalSource(0))
         return
      end if

      nodSrcs = this%jsonValueFilterByKeyValues(sources, J_TYPE, [J_SRC_TYPE_NS])
      if (size(nodSrcs) == 0) then
         allocate(res%NodalSource(0))
         return
      end if

      allocate(res%NodalSource(size(nodSrcs)))
      res%n_nodSrc = size(nodSrcs)
      res%n_nodSrc_max = size(nodSrcs)
      do i = 1, size(nodSrcs)
         res%NodalSource(i) = readField(nodSrcs(i)%p)
      end do
      do i = 1, size(nodSrcs)
         res%n_C2P_max = max(res%n_C2P_max, res%NodalSource(i)%n_C2P)
      end do

   contains
      function readField(jns) result(res)
         type(Curr_Field_Src) :: res
         type(json_value), pointer :: jns, entry
         integer, dimension(:), allocatable :: elementIds
         type(coords_scaled), dimension(:), allocatable :: coordsFromLinels

         select case (this%getStrAt(jns, J_FIELD))
          case (J_FIELD_ELECTRIC)
            res%isField = .true.
            res%isElec = .true.
            res%isMagnet = .false.
            res%isCurrent = .false.
          case (J_FIELD_MAGNETIC)
            res%isField = .true.
            res%isElec = .false.
            res%isMagnet = .true.
            res%isCurrent = .false.
          case default
            write(error_unit, *) 'Error reading current field source. Field label not recognized.'
         end select
         res%isInitialValue = .false.
         res%nombre = trim(adjustl(this%getStrAt(jns, J_SRC_MAGNITUDE_FILE)))

         allocate(res%c1P(0))
         res%n_C1P = 0

         elementIds = this%getIntsAt(jns, J_ELEMENTIDS)
         call cellRegionsToScaledCoords(res%C2P,  this%mesh%getCellRegions(elementIds) )
         res%n_C2P = size(res%C2p)

      end function
   end function

   function readProbes(this) result (res)
      class(parser_t) :: this
      type(Sondas) :: res
      type(json_value), pointer :: allProbes
      type(json_value_ptr), dimension(:), allocatable :: ps
      ! The only oldProbe present in the format is the far field.
      character (len=*), dimension(1), parameter :: validTypes = [J_PR_TYPE_FARFIELD]
      integer :: i
      logical :: found

      call this%core%get(this%root, J_PROBES, allProbes, found)
      if (.not. found) then
         allocate(res%probes(0))
         res%n_probes = size(res%probes)
         res%n_probes_max = size(res%probes)
         return
      end if

      ps = this%jsonValueFilterByKeyValues(allProbes, J_TYPE, validTypes)

      res%n_probes = size(ps)
      res%n_probes_max = size(ps)
      allocate(res%probes(size(ps)))
      do i=1, size(ps)
         res%probes(i) = readFarFieldProbe(ps(i)%p)
      end do

   contains
      function readFarFieldProbe(p) result (res)
         type(abstractSonda) :: res
         type(json_value), pointer :: p
         type(Sonda), pointer :: ff
         character (len=:), allocatable :: outputName
         logical :: transferFunctionFound
         type(domain_t) :: domain

         res%n_FarField = 1
         res%n_FarField_max = 1
         allocate(res%FarField(1))
         ff => res%FarField(1)%probe

         ff%grname = " "
         outputName = this%getStrAt(p, J_NAME)
         ff%outputrequest = trim(adjustl(outputName))

         ! Far fields only accept frequency domains.
         domain = this%getDomain(p, J_PR_DOMAIN)
         if (domain%type2 /= NP_T2_FREQ) then
            write(error_unit, *) "ERROR at far field probe: Only accepted domain is frequency."
         end if
         ff%tstart = 0.0
         ff%tstop = 0.0
         ff%tstep = 0.0
         ff%fstart = domain%fstart
         ff%fstop = domain%fstop
         ff%fstep = domain%fstep

         block
            logical :: sourcesFound
            type(json_value), pointer :: sources, src
            character (len=:), allocatable :: fn

            fn = this%getStrAt(p, J_PR_DOMAIN//J_PR_DOMAIN_MAGNITUDE_FILE, found=transferFunctionFound)
            if (.not. transferFunctionFound) then
               call this%core%get(this%root, J_SOURCES, sources, sourcesFound)
               if (sourcesFound) then
                  if (this%core%count(sources) == 1) then
                     call this%core%get_child(sources, 1, src)
                     fn = this%getStrAt(src, J_SRC_MAGNITUDE_FILE, found=transferFunctionFound)
                  end if
               end if
            end if

            if (transferFunctionFound) then
               ff%FileNormalize = trim(adjustl(fn))
            else
               ff%FileNormalize = " "
            end if

         end block

         if (domain%isLogarithmicFrequencySpacing) then
            call appendLogSufix(ff%outputrequest)
         end if

         block
            type(coords), dimension(:), allocatable :: nfdeCoords
            nfdeCoords = &
               cellIntervalsToCoords(this%getSingleVolumeInElementsIds(p))
            ff%n_cord = 2
            ff%n_cord_max = 2
            allocate(ff%i(2))
            allocate(ff%j(2))
            allocate(ff%k(2))
            allocate(ff%node(0))
            ff%i(1) = nfdeCoords(1)%Xi
            ff%i(2) = nfdeCoords(1)%Xe
            ff%j(1) = nfdeCoords(1)%Yi
            ff%j(2) = nfdeCoords(1)%Ye
            ff%k(1) = nfdeCoords(1)%Zi
            ff%k(2) = nfdeCoords(1)%Ze
         end block

         block
            call readDirection(&
               p, J_PR_FAR_FIELD_PHI, ff%phistart, ff%phistop, ff%phistep)
            call readDirection(&
               p, J_PR_FAR_FIELD_THETA, ff%thetastart, ff%thetastop, ff%thetastep)
         end block
      end function

      subroutine readDirection(p, label, initial, final, step)
         type(json_value), pointer :: p
         type(json_value), pointer :: dir
         character (len=*), intent(in) :: label
         logical :: found
         real (kind=rkind), intent(inout) :: initial, final, step

         call this%core%get(p, label, dir, found=found)
         if (.not. found) &
            write (error_unit, *) "Error reading far field probe. Direction label not found."
         initial = this%getRealAt(dir, J_PR_FAR_FIELD_DIR_INITIAL)
         final   = this%getRealAt(dir, J_PR_FAR_FIELD_DIR_FINAL)
         step    = this%getRealAt(dir, J_PR_FAR_FIELD_DIR_STEP)
      end subroutine
   end function

   function readMoreProbes(this) result (res)
      class(parser_t) :: this
      type(MasSondas) :: res
      type(json_value), pointer :: allProbes
      type(json_value_ptr), dimension(:), allocatable :: ps

      integer :: i
      character (len=*), dimension(2), parameter :: validTypes = &
         [J_PR_TYPE_POINT, J_PR_TYPE_WIRE]
      character (len=*), dimension(1), parameter :: validFields = &
         [J_FIELD_CURRENT]
      logical :: found
      character (len=:), allocatable :: fieldLbl
      integer :: filtered_size, n

      call this%core%get(this%root, J_PROBES, allProbes, found)
      if (.not. found) then
         allocate(res%collection(0))
         res%length = size(res%collection)
         res%length_max = size(res%collection)
         res%len_cor_max = 0
         return
      end if

      ps = this%jsonValueFilterByKeyValues(allProbes, J_TYPE, validTypes)
      
      filtered_size = 0
      do i=1, size(ps)
         fieldLbl = this%getStrAt(ps(i)%p, J_FIELD, default=J_FIELD_ELECTRIC)
         if (fieldLbl /= J_FIELD_VOLTAGE) then 
            filtered_size = filtered_size + 1
         end if
      end do

      n = 1
      allocate(res%collection(filtered_size))
      do i=1, size(ps)
         fieldLbl = this%getStrAt(ps(i)%p, J_FIELD, default=J_FIELD_ELECTRIC)
         if (fieldLbl /= J_FIELD_VOLTAGE) then 
            res%collection(n) = readPointProbe(ps(i)%p)
            n = n + 1
         end if
      end do

      res%length = size(res%collection)
      res%length_max = size(res%collection)
      res%len_cor_max = 0
   contains
      function readPointProbe(p) result (res)
         type(MasSonda) :: res
         type(json_value), pointer :: p, dirLabelPtr
         character(len=1), dimension(:), allocatable :: dirLabels
         integer :: i, j, k
         character (len=:), allocatable :: typeLabel, fieldLabel, outputName, dirLabel
         type(pixel_t) :: pixel

         integer, dimension(:), allocatable :: elemIds
         logical :: elementIdsFound, typeLabelFound, dirLabelsFound, fieldLabelFound, nameFound

         outputName = this%getStrAt(p, J_NAME, found=nameFound)
         if (.not. nameFound) then 
            write(error_unit, *) "ERROR: name entry not found for probe."
         end if
         res%outputrequest = trim(adjustl(outputName))

         call setDomain(res, this%getDomain(p, J_PR_DOMAIN))

         elemIds = this%getIntsAt(p, J_ELEMENTIDS, found=elementIdsFound)
         if (.not. elementIdsFound) then
            write(error_unit, *) "ERROR: element ids entry not found for probe."
         end if
         if (size(elemIds) /= 1) then
            write(error_unit, *) "ERROR: point probe must contain a single element id."
         end if

         pixel = getPixelFromElementId(this%mesh, elemIds(1))

         typeLabel = this%getStrAt(p, J_TYPE, found=typeLabelFound)
         if (.not. typeLabelFound) then
            write(error_unit, *) "ERROR: Point probe type label not found."
         end if
         select case (typeLabel)
          case (J_PR_TYPE_WIRE)
            allocate(res%cordinates(1))
            fieldLabel = this%getStrAt(p, J_FIELD, default=J_FIELD_VOLTAGE)
            res%cordinates(1)%tag = outputName
            res%cordinates(1)%Xi = pixel%tag
            res%cordinates(1)%Yi = 0
            res%cordinates(1)%Zi = 0
            res%cordinates(1)%Or = strToFieldType(fieldLabel)
          case (J_PR_TYPE_POINT)
            call this%core%get(p, J_PR_POINT_DIRECTIONS, dirLabelPtr, found=dirLabelsFound)
            if(dirLabelsFound) then
               dirLabels = buildDirLabels(dirLabelPtr)
            else 
               dirLabels = [J_DIR_X, J_DIR_Y, J_DIR_Z]
            end if
            fieldLabel = this%getStrAt(p, J_FIELD, default=J_FIELD_ELECTRIC, found=fieldLabelFound)          
            allocate(res%cordinates(size(dirLabels)))
            do j = 1, size(dirLabels)
               res%cordinates(j)%tag = outputName
               res%cordinates(j)%Xi = int (pixel%cell(1))
               res%cordinates(j)%Yi = int (pixel%cell(2))
               res%cordinates(j)%Zi = int (pixel%cell(3))
               res%cordinates(j)%Or = strToFieldType(fieldLabel, dirLabels(j))
            end do
         end select

         res%len_cor = size(res%cordinates)
      end function

      function buildDirLabels(dirLabelsPtr) result (res)
         type(json_value), pointer, intent(in) :: dirLabelsPtr
         character(len=1), dimension(:), allocatable :: res
         type(json_value), pointer :: child
         character(len=:), allocatable :: str
         integer :: i
         allocate(res(this%core%count(dirLabelsPtr)))
         do i = 1, this%core%count(dirLabelsPtr)
            call this%core%get_child(dirLabelsPtr, i, child)
            call this%core%get(child, str)
            res(i) = str
         end do
      end function

      subroutine setDomain(res, domain)
         type(MasSonda), intent(inout) :: res
         type(domain_t), intent(in) :: domain

         res%tstart = domain%tstart
         res%tstep  = domain%tstep
         res%tstop  = domain%tstop
         res%fstart = domain%fstart
         res%fstep  = domain%fstep
         res%fstop  = domain%fstop
         if (allocated(domain%filename)) then
            res%filename = domain%filename
         else
            res%filename = " "
         end if
         res%type1  = domain%type1
         res%type2  = domain%type2

         if (domain%isLogarithmicFrequencySpacing) then
            call appendLogSufix(res%outputrequest)
         end if
      end subroutine

      function strToFieldType(fieldLabel, dirLabel) result(res)
         integer (kind=4) :: res
         character (len=:), allocatable, intent(in) :: fieldLabel
         character (len=1), intent(in), optional :: dirLabel
         select case (fieldLabel)
          case (J_FIELD_ELECTRIC)
            if (.not. present(dirLabel)) then
               write(error_unit, *) "Dir label must be present"
            end if
            select case (dirLabel)
             case (J_DIR_X)
               res = NP_COR_EX
             case (J_DIR_Y)
               res = NP_COR_EY
             case (J_DIR_Z)
               res = NP_COR_EZ
             case default
               write(error_unit, *) "Invalid dir label"
            end select
          case (J_FIELD_MAGNETIC)
            if (.not. present(dirLabel)) then
               write(error_unit, *) "Dir label must be present"
            end if
            select case (dirLabel)
             case (J_DIR_X)
               res = NP_COR_HX
             case (J_DIR_Y)
               res = NP_COR_HY
             case (J_DIR_Z)
               res = NP_COR_HZ
             case default
               write(error_unit, *) "Invalid dir label"
            end select
          case (J_FIELD_CURRENT)
            res = NP_COR_WIRECURRENT
          case (J_FIELD_VOLTAGE)
            res = NP_COR_DDP
          case default
            write(error_unit,*) "Invalid field label for point/wire probe."
         end select

      end function
   end function

   function readBlockProbes(this) result (res)
      class(parser_t) :: this
      type(BloqueProbes) :: res
      type(json_value_ptr), dimension(:), allocatable :: bps
      type(json_value), pointer :: probes
      logical :: found
      integer :: i

      call this%core%get(this%root, J_PROBES, probes, found)
      if (.not. found) then
         allocate(res%bp(0))
         return
      end if

      bps = this%jsonValueFilterByKeyValues(probes, J_TYPE, [J_PR_TYPE_BULK_CURRENT])
      if (size(bps) == 0) then
         allocate(res%bp(0))
         return
      end if

      res%n_bp = size(bps)
      res%n_bp_max = size(bps)
      allocate(res%bp(size(bps)))
      do i = 1, size(bps)
         res%bp(i) = readBlockProbe(bps(i)%p)
      end do
   contains
      function readBlockProbe(bp) result(res)
         type(BloqueProbe) :: res
         type(json_value), pointer :: bp
         type(coords), dimension(:), allocatable :: cs
         type(cell_region_t), dimension(:), allocatable :: cRs

         cRs = this%mesh%getCellRegions(this%getIntsAt(bp, J_ELEMENTIDS))
         if (size(cRs) /= 1) write(error_unit, *) "Bulk current probe must be defined by a single cell region."

         if (size(cRs(1)%intervals) /= 1) write(error_unit, *) "Bulk current probe must be defined by a single cell interval."
         cs = cellIntervalsToCoords(cRs(1)%intervals)

         res%i1  = cs(1)%xi
         res%i2  = cs(1)%xe
         res%j1  = cs(1)%yi
         res%j2  = cs(1)%ye
         res%k1  = cs(1)%zi
         res%k2  = cs(1)%ze
         res%nml = cs(1)%Or

         res%outputrequest = trim(adjustl(this%getStrAt(bp, J_NAME)))
         call setDomain(res, this%getDomain(bp, J_PR_DOMAIN))

         res%skip = 1
         res%tag = trim(adjustl(this%getStrAt(bp, J_NAME, default=" ")))
         res%t = BcELECT

      end function

      subroutine setDomain(res, domain)
         type(BloqueProbe), intent(inout) :: res
         type(domain_t), intent(in) :: domain

         res%tstart = domain%tstart
         res%tstep = domain%tstep
         res%tstop = domain%tstop
         res%fstart = domain%fstart
         res%fstep = domain%fstep
         res%fstop = domain%fstop
         if (allocated(domain%filename)) then
            res%FileNormalize = domain%filename
         else
            res%fileNormalize = " "
         end if
         res%type2 = domain%type2

         if (domain%isLogarithmicFrequencySpacing) then
            call appendLogSufix(res%outputrequest)
         end if
      end subroutine
   end function

   function readVolumicProbes(this) result (res)
      class(parser_t) :: this
      type(VolProbes) :: res
      type(json_value_ptr), dimension(:), allocatable :: ps
      type(json_value), pointer :: probes
      logical :: found
      integer :: i

      call this%core%get(this%root, J_PROBES, probes, found)
      if (.not. found) then
         res = buildNoVolProbes()
         return
      end if

      ps = this%jsonValueFilterByKeyValues(probes, J_TYPE, [J_PR_TYPE_MOVIE])
      if (size(ps) == 0) then
         res = buildNoVolProbes()      
         return
      end if

      res%length = size(ps)
      res%length_max = size(ps)
      res%len_cor_max = 2*size(ps)
      allocate(res%collection(size(ps)))
      do i = 1, size(ps)
         res%collection(i) = readVolProbe(ps(i)%p)
      end do

   contains
      function buildNoVolProbes() result(res)
         type(VolProbes) :: res
         allocate(res%collection(0))
         res%length = 0
         res%length_max = 0
         res%len_cor_max = 0
      end function

      function readVolProbe(p) result(res)
         type(VolProbe) :: res
         type(json_value), pointer :: p, compsPtr, compPtr
         type(coords), dimension(:), allocatable :: cs
         type(cell_region_t), dimension(:), allocatable :: cRs
         character(len=:), allocatable :: fieldType, component
         integer :: i
         integer :: numberOfComponents
         logical :: componentsFound

         cRs = this%mesh%getCellRegions(this%getIntsAt(p, J_ELEMENTIDS))
         if (size(cRs) /= 1) then
            write(error_unit, *) "Movie probe must be defined over a single cell region."
         end if

         if (size(cRs(1)%intervals) /= 1) then
            write(error_unit, *) "Movie probe must be defined by a single cell interval."
         end if
         cs = cellIntervalsToCoords(cRs(1)%intervals)

         fieldType = this%getStrAt(p, J_FIELD, default=J_FIELD_ELECTRIC)
         call this%core%get(p, J_PR_MOVIE_COMPONENT, compsPtr, found=componentsFound)
         allocate(res%cordinates(1))
         if (componentsFound) then
            call this%core%get(compsPtr, component)
            res%cordinates(1) = cs(1)
            res%cordinates(1)%Or  = buildVolProbeType(fieldType, component)
         else 
            component = J_DIR_M
            res%cordinates(1)%Or  = buildVolProbeType(fieldType, component)
         endif
         res%len_cor = size(res%cordinates)
         
         res%outputrequest = trim(adjustl(this%getStrAt(p, J_NAME, default=" ")))
         call setDomain(res, this%getDomain(p, J_PR_DOMAIN))
      end function

      integer function buildVolProbeType(fieldType, component) result(res)
         character(len=:), allocatable, intent(in) :: fieldType, component
         select case (fieldType)
         case (J_FIELD_ELECTRIC)
            select case (component)
            case (J_DIR_X)
               res = iExC
            case (J_DIR_Y)
               res = iEyC
            case (J_DIR_Z)
               res = iEzC
            case (J_DIR_M)
               res = iMEC
            end select
         case (J_FIELD_MAGNETIC)
            select case (component)
            case (J_DIR_X)
               res = iHxC
            case (J_DIR_Y)
               res = iHyC
            case (J_DIR_Z)
               res = iHzC
            case (J_DIR_M)
               res = iMHC
            end select
         case (J_FIELD_CURRENT_DENSITY)
            select case (component)
            case (J_DIR_X)
               res = iCurX
            case (J_DIR_Y)
               res = iCurY
            case (J_DIR_Z)
               res = iCurZ
            case (J_DIR_M)
               res = iCur
            end select
         case default
            write(error_unit,*) "ERROR Determining vol probe type: invalid field type."
         end select
      end function

      subroutine setDomain(res, domain)
         type(VolProbe), intent(inout) :: res
         type(domain_t), intent(in) :: domain

         res%tstart = domain%tstart
         res%tstep = domain%tstep
         res%tstop = domain%tstop
         res%fstart = domain%fstart
         res%fstep = domain%fstep
         res%fstop = domain%fstop
         if (allocated(domain%filename)) then
            res%filename = domain%filename
         else
            res%filename = " "
         end if
         res%type2 = domain%type2

         if (domain%isLogarithmicFrequencySpacing) then
            call appendLogSufix(res%outputrequest)
         end if
      end subroutine
   end function

   subroutine appendLogSufix(fn) 
      character(len=BUFSIZE), intent(inout) :: fn
      character (len=*), parameter :: SMBJSON_LOG_SUFFIX = "_log_"
      fn = trim(fn) // SMBJSON_LOG_SUFFIX
   end subroutine

   function readThinSlots(this) result (res)
      class(parser_t) :: this
      type(ThinSlots) :: res
      
      type(materialAssociation_t), dimension(:), allocatable :: mAs
      integer :: i

      mAs = this%getMaterialAssociations([J_MAT_TYPE_SLOT])
      if (size(mAs) == 0) then
         allocate(res%tg(0))
         return
      end if

      res%n_tg = size(mAs)
      allocate(res%tg(res%n_tg))
      do i = 1, size(mAs)
         res%tg = readThinSlot(mAs(i))
      end do
   contains
      function readThinSlot(mA) result(res)
         type (materialAssociation_t), intent(in) :: mA
         type (thinSlot) :: res
         type (coords), dimension(:), pointer :: cs
         type(json_value_ptr) :: mat
         logical :: found
         
         mat = this%matTable%getId(mA%materialId)
         res%width = this%getRealAt(mat%p, J_MAT_THINSLOT_WIDTH, found)
         if (.not. found) then
            write(error_unit, *) "ERROR reading thin slot: ", &
               J_MAT_THINSLOT_WIDTH, "not found"
         end if

         call this%matAssToCoords(cs, mA, CELL_TYPE_LINEL)
         call coordsToThinSlotComp(res%tgc, cs)
         res%n_tgc = size(res%tgc)

      end function

      subroutine coordsToThinSlotComp(tc, cs)
         type(coords), dimension(:), pointer, intent(in) :: cs
         type(thinSlotComp), dimension(:), pointer :: tc
         integer :: i, j, k
         integer :: nTgc, nXYZ
         integer :: dir
         ! Precount
         nTgc = 0
         do i = 1, size(cs)
            nXYZ =  (cs(i)%xe - cs(i)%xi + 1) * &
                    (cs(i)%ye - cs(i)%yi + 1) * &
                    (cs(i)%ze - cs(i)%zi + 1) 
            nTgc = nTgc + nXYZ
         end do

         ! Fill
         j = 1
         allocate(tc(nTgc))
         do i = 1, size(cs)
            select case (abs(cs(i)%Or))
            case (iEx)
               do k = 1, (cs(i)%xe - cs(i)%xi + 1)
                  tc(j) = buildBaseThinSlotComponent(cs(i))
                  tc(j)%i = cs(i)%xi + k - 1
                  j = j + 1
               end do
            case (iEy)
               do k = 1, (cs(i)%xe - cs(i)%xi + 1)
                  tc(j) = buildBaseThinSlotComponent(cs(i))
                  tc(j)%j = cs(i)%yi + k - 1
                  j = j + 1
               end do
            case (iEz)
               do k = 1, (cs(i)%xe - cs(i)%xi + 1)
                  tc(j) = buildBaseThinSlotComponent(cs(i))
                  tc(j)%k = cs(i)%zi + k - 1
                  j = j + 1
               end do
            end select
         end do
      end subroutine

      function buildBaseThinSlotComponent(cs) result(res)
         type(coords), intent(in) :: cs
         type(thinSlotComp) :: res
         res%i = cs%xi
         res%j = cs%yi
         res%k = cs%zi
         res%dir = abs(cs%Or)
         res%tag = cs%tag
      end function
   end function

   function readThinWires(this) result (res)
      class(parser_t) :: this
      type(ThinWires) :: res
      type(materialAssociation_t), dimension(:), allocatable :: mAs
      integer :: i, j
      logical :: found

      mAs = this%getMaterialAssociations([J_MAT_TYPE_WIRE])

      ! Pre-allocates thin wires.
      block
         integer :: nTw
         nTw = 0
         if (size(mAs) /=0 ) then
            do i = 1, size(mAs)
               if (isThinWire(mAs(i))) nTw = nTw+1
            end do
         end if

         allocate(res%tw(nTw))
         res%n_tw = size(res%tw)
         res%n_tw_max = size(res%tw)
      end block

      j = 1
      if (size(mAs) /=0 ) then
         do i = 1, size(mAs)
            if (isThinWire(mAs(i))) then
               res%tw(j) = readThinWire(mAs(i))
               j = j+1
            end if
         end do
      end if

   contains
      function readThinWire(cable) result(res)
         type(ThinWire) :: res
         type(materialAssociation_t), intent(in) :: cable

         character (len=:), allocatable :: entry
         type(json_value), pointer :: je, je2
         integer :: i
         logical :: found
         real :: radius, resistance, inductance
         block
            type(json_value_ptr) :: m
            m = this%matTable%getId(cable%materialId)
            radius = this%getRealAt(m%p, J_MAT_WIRE_RADIUS, default = 0.0)
            resistance = this%getRealAt(m%p, J_MAT_WIRE_RESISTANCE, default = 0.0)
            inductance = this%getRealAt(m%p, J_MAT_WIRE_INDUCTANCE, default = 0.0)
            res%rad = radius 
            res%res = resistance
            res%ind = inductance
            res%dispfile = trim(adjustl(" "))
         end block

         block
            type(json_value_ptr) :: terminal
            type(thinwiretermination_t) :: term
            character (len=:), allocatable :: label
            terminal = this%matTable%getId(cable%initialTerminalId)
            term = readThinWireTermination(terminal%p)
            res%tl = term%terminationType
            res%R_LeftEnd = term%r
            res%L_LeftEnd = term%l
            res%C_LeftEnd = term%c
            res%dispfile_LeftEnd = trim(adjustl(" "))
         end block

         block
            type(json_value_ptr) :: terminal
            type(thinwiretermination_t) :: term
            terminal = this%matTable%getId(cable%endTerminalId)
            term = readThinWireTermination(terminal%p)
            res%tr = term%terminationType
            res%R_RightEnd = term%r
            res%L_RightEnd = term%l
            res%C_RightEnd = term%c
            res%dispfile_RightEnd = trim(adjustl(" "))
         end block

         block
            type(linel_t), dimension(:), allocatable :: linels
            type(polyline_t) :: polyline
            character (len=MAX_LINE) :: tagLabel
            type(generator_description_t), dimension(:), allocatable :: genDesc
            
            polyline = this%mesh%getPolyline(cable%elementIds(1))
            linels = this%mesh%polylineToLinels(polyline)

            write(tagLabel, '(i10)') cable%elementIds(1)

            genDesc = readGeneratorOnThinWire(linels, cable%elementIds)

            res%n_twc = size(linels)
            res%n_twc_max = size(linels)
            allocate(res%twc(size(linels)))
            do i = 1, size(linels)
               res%twc(i)%srcfile = genDesc(i)%srcfile
               res%twc(i)%srctype = genDesc(i)%srctype
               res%twc(i)%m = genDesc(i)%multiplier
               res%twc(i)%i = linels(i)%cell(1)
               res%twc(i)%j = linels(i)%cell(2)
               res%twc(i)%k = linels(i)%cell(3)
               res%twc(i)%d = abs(linels(i)%orientation)
               res%twc(i)%nd = linels(i)%tag
               res%twc(i)%tag = trim(adjustl(tagLabel))
            end do
         end block

      end function

      function readGeneratorOnThinWire(linels, plineElemIds) result(res)
         type(linel_t), dimension(:), intent(in) :: linels
         integer, dimension(:), intent(in) :: plineElemIds
         type(json_value), pointer :: sources
         type(json_value_ptr), dimension(:), allocatable :: genSrcs
         logical :: found
         type(generator_description_t), dimension(:), allocatable :: res
         integer :: i

         allocate(res(size(linels)))
         do i = 1, size(linels)
            res(i)%srcfile = 'None'
            res(i)%srctype = 'None'
            res(i)%multiplier = 0.0
         end do

         call this%core%get(this%root, J_SOURCES, sources, found)
         if (.not. found) then
            return
         end if

         genSrcs = this%jsonValueFilterByKeyValues(sources, J_TYPE, [J_SRC_TYPE_GEN])
         if (size(genSrcs) == 0) then
            return
         end if

         block
            integer, dimension(:), allocatable :: sourceElemIds
            integer :: position
            type(node_t) :: srcCoord
            type(polyline_t) :: polylineCoords
            do i = 1, size(genSrcs)
               sourceElemIds = this%getIntsAt(genSrcs(i)%p, J_ELEMENTIDS)
               srcCoord = this%mesh%getNode(sourceElemIds(1))
               polylineCoords = this%mesh%getPolyline(plineElemIds(1))
               if (.not. any(polylineCoords%coordIds == srcCoord%coordIds(1))) then
                  cycle ! generator is not in this polyline
               end if

               position = findSourcePositionInLinels(sourceElemIds, linels)

               if (.not. this%existsAt(genSrcs(i)%p, J_SRC_MAGNITUDE_FILE)) then
                  write(error_unit, *) 'magnitudeFile of source missing'
                  return
               end if

               select case(this%getStrAt(genSrcs(i)%p, J_FIELD))
                case (J_FIELD_VOLTAGE)
                  res(position)%srctype = "VOLT"
                  res(position)%srcfile = this%getStrAt(genSrcs(i)%p, J_SRC_MAGNITUDE_FILE)
                  res(position)%multiplier = 1.0
                case (J_FIELD_CURRENT)
                  res(position)%srctype = "CURR"
                  res(position)%srcfile = this%getStrAt(genSrcs(i)%p, J_SRC_MAGNITUDE_FILE)
                  res(position)%multiplier = 1.0
                case default
                  write(error_unit, *) 'Field block of source of type generator must be current or voltage'
               end select

            end do
         end block

      end function

      function findSourcePositionInLinels(srcElemIds, linels) result(res)
         integer, dimension(:), intent(in) :: srcElemIds
         type(linel_t), dimension(:), intent(in) :: linels
         type(pixel_t) :: pixel
         integer :: res
         integer :: i
         pixel = this%mesh%nodeToPixel(this%mesh%getNode(srcElemIds(1)))
         do i = 1, size(linels)
            if (linels(i)%tag == pixel%tag) then
               res = i
               return
            end if
         end do
         write (error_unit, * ) "ERROR: Source could not be found in linels."

      end function

      function readThinWireTermination(terminal) result(res)
         type(thinwiretermination_t) :: res
         type(json_value), pointer :: terminal, tms, tm
         character (len=:), allocatable :: label
         logical :: found

         call this%core%get(terminal, J_MAT_TERM_TERMINATIONS, tms, found)

         if (.not. found) then
            write(error_unit, *) "Error reading wire terminal. terminations not found."
         end if
         if (this%core%count(tms) /= 1) then
            write(error_unit, *) "Only terminals with a single termination are allowed for a wire."
         end if

         call this%core%get_child(tms, 1, tm)

         label = this%getStrAt(tm, J_TYPE, found)
         res%terminationType = strToTerminationType(label)
         if (.not. found) then
            write(error_unit, *) "Error reading wire terminal. termination must specify a type."
         end if

         select case(label)
          case(J_MAT_TERM_TYPE_OPEN)
            res%r = 0.0
            res%l = 0.0
            res%c = 0.0
          case default
            res%r = this%getRealAt(tm, J_MAT_TERM_RESISTANCE, default=0.0)
            res%l = this%getRealAt(tm, J_MAT_TERM_INDUCTANCE, default=0.0)
            res%c = this%getRealAt(tm, J_MAT_TERM_CAPACITANCE, default=1e22)
         end select

      end function

      function strToTerminationType(label) result(res)
         character (len=:), allocatable, intent(in) :: label
         integer :: res
         select case (label)
          case (J_MAT_TERM_TYPE_OPEN)
            res = MATERIAL_CONS
          case (J_MAT_TERM_TYPE_SERIES)
            res = SERIES_CONS
          case (J_MAT_TERM_TYPE_SHORT)
            res = MATERIAL_CONS
         end select
      end function

      logical function isThinWire(mA)
         type(materialAssociation_t) :: mA
         type(json_value_ptr) :: mat
         type(polyline_t) :: pl
         logical :: found
         isThinWire = .false.

         if (size(mA%elementIds) /= 1) then
            write(error_unit, *) "ERROR: Thin wires must be defined by a single element id."
         end if

         pl = this%mesh%getPolyline(mA%elementIds(1))
         if (.not. this%mesh%arePolylineSegmentsStructured(pl)) then
            write(error_unit, *) "ERROR: Thin wires must be defined by a structured polyline."
         end if
      
         isThinWire = .true.
      end function
   end function

   function getDomain(this, place, path) result(res)
      class(parser_t) :: this
      type(domain_t) :: res
      type(json_value), pointer :: place
      character(len=*), intent(in) :: path

      integer :: numberOfFrequencies
      type(json_value), pointer :: domain
      character (len=:), allocatable :: fn, domainType, freqSpacing
      logical :: found, transferFunctionFound
      real :: val

      call this%core%get(place, path, domain, found)
      if (.not. found) then
         res%filename = " "
         return
      end if

      fn = this%getStrAt(domain, J_PR_DOMAIN_MAGNITUDE_FILE, transferFunctionFound, default=" ")
      if (transferFunctionFound) then
         res%filename = trim(adjustl(fn))
      endif

      res%type1 = NP_T1_PLAIN

      domainType = this%getStrAt(domain, J_TYPE, default=J_PR_DOMAIN_TYPE_TIME)
      res%type2 = getNPDomainType(domainType, transferFunctionFound)

      res%tstart = this%getRealAt(domain, J_PR_DOMAIN_TIME_START, default=0.0)
      res%tstop = this%getRealAt(domain, J_PR_DOMAIN_TIME_STOP, default=0.0)
      res%tstep = this%getRealAt(domain, J_PR_DOMAIN_TIME_STEP, default=0.0)
      res%fstart = this%getRealAt(domain, J_PR_DOMAIN_FREQ_START, default=0.0)
      res%fstop = this%getRealAt(domain, J_PR_DOMAIN_FREQ_STOP, default=0.0)

      numberOfFrequencies = this%getIntAt(domain, J_PR_DOMAIN_FREQ_NUMBER, default=0)
      if (numberOfFrequencies == 0) then
         res%fstep = 0.0
      else
         res%fstep = res%fstart * numberOfFrequencies
      endif

      freqSpacing = &
         this%getStrAt(domain, J_PR_DOMAIN_FREQ_SPACING, default=J_PR_DOMAIN_FREQ_SPACING_LINEAR)
      select case (freqSpacing)
       case (J_PR_DOMAIN_FREQ_SPACING_LINEAR)
         res%isLogarithmicFrequencySpacing = .false.
       case (J_PR_DOMAIN_FREQ_SPACING_LOGARITHMIC)
         res%isLogarithmicFrequencySpacing = .true.
      end select

   contains
      function getNPDomainType(typeLabel, hasTransferFunction) result(res)
         integer (kind=4) :: res
         character (len=:), intent(in), allocatable :: typeLabel
         logical, intent(in) :: hasTransferFunction
         logical :: isTime, isFrequency
         select case(typeLabel)
          case (J_PR_DOMAIN_TYPE_TIME)
            isTime = .true.
            isFrequency = .false.
          case (J_PR_DOMAIN_TYPE_FREQ)
            isTime = .false.
            isFrequency = .true.
          case (J_PR_DOMAIN_TYPE_TIMEFREQ)
            isTime = .true.
            isFrequency = .true.
         end select

         if (           isTime .and. .not. isFrequency .and. .not. hasTransferFunction) then
            res = NP_T2_TIME
            return
         else if (.not. isTime .and.       isFrequency .and. .not. hasTransferFunction) then
            res = NP_T2_FREQ
            return
         else if (.not. isTime .and. .not. isFrequency .and.       hasTransferFunction) then
            res = NP_T2_TRANSFER
            return
         else if (      isTime .and.       isFrequency .and. .not. hasTransferFunction) then
            res = NP_T2_TIMEFREQ
            return
         else if (      isTime .and. .not. isFrequency .and.       hasTransferFunction) then
            res = NP_T2_TIMETRANSF
            return
         else if (.not. isTime .and.       isFrequency .and.       hasTransferFunction) then
            res = NP_T2_FREQTRANSF
            return
         else if (      isTime .and.       isFrequency .and.       hasTransferFunction) then
            res = NP_T2_TIMEFRECTRANSF
            return
         end if

         write(error_unit, *) "Error parsing domain."
      end function
   end function

   function parseMaterialAssociation(this, matAss) result(res)
      class(parser_t) :: this
      type(json_value), pointer, intent(in) :: matAss
      type(json_value_ptr) :: mat
      type(materialAssociation_t) :: res
      character (len=*), parameter :: errorMsgInit = "ERROR reading material association: "
      logical :: found
      logical :: isMultiwire, isWireOrMultiwire
      
      ! Fills material association.
      res%materialId = this%getIntAt(matAss, J_MATERIAL_ID, found)
      if (.not. found) call showLabelNotFoundError(J_MATERIAL_ID)
      
      res%elementIds = this%getIntsAt(matAss, J_ELEMENTIDS, found)
      if (.not. found) call showLabelNotFoundError(J_ELEMENTIDS)
      
      res%name = this%getStrAt(matAss, J_NAME, found)
      if (.not. found) then
         res%name = ""
      end if

      res%initialTerminalId  = this%getIntAt(matAss, J_MAT_ASS_CAB_INI_TERM_ID, default=-1)
      res%endTerminalId      = this%getIntAt(matAss, J_MAT_ASS_CAB_END_TERM_ID, default=-1)
      res%initialConnectorId = this%getIntAt(matAss, J_MAT_ASS_CAB_INI_CONN_ID, default=-1)
      res%endConnectorId     = this%getIntAt(matAss, J_MAT_ASS_CAB_END_CONN_ID, default=-1)
      res%containedWithinElementId = &
                  this%getIntAt(matAss, J_MAT_ASS_CAB_CONTAINED_WITHIN_ID, default=-1)

      ! Checks validity of associations.
      if (this%matTable%checkId(res%materialId) /= 0) then
         write(error_unit, *) errorMsgInit, "material with id ", res%materialId, " not found."
      endif
      
      if (size(res%elementIds) == 0) then
         write(error_unit, *) errorMsgInit, J_ELEMENTIDS, "must not be empty."
      end if
      block
         integer :: i
         do i = 1, size(res%elementIds)
            if (this%mesh%checkElementId(res%elementIds(i)) /= 0) then
               write(error_unit, *) errorMsgInit, "element with id ", res%elementIds(i), " not found."
            end if
         end do
      end block

      mat = this%matTable%getId(res%materialId)
      isMultiwire = this%getStrAt(mat%p, J_TYPE) == J_MAT_TYPE_MULTIWIRE
      isWireOrMultiwire = &
         this%getStrAt(mat%p, J_TYPE) == J_MAT_TYPE_WIRE .or. isMultiwire 
      
      if (isWireOrMultiwire) then
         if (res%initialTerminalId == -1 .or. res%endTerminalId == -1) then
            write(error_unit, *), errorMsgInit, "wire associations must include terminals."
         end if
         if (.not. isMaterialIdOfType(res%initialTerminalId, J_MAT_TYPE_TERMINAL)) then
            write(error_unit, *) errorMsgInit, "material with id ", res%materialId, " must be a terminal."
         end if
         if (.not. isMaterialIdOfType(res%endTerminalId, J_MAT_TYPE_TERMINAL)) then
            write(error_unit, *) errorMsgInit, "material with id ", res%materialId, " must be a terminal."
         end if
         if (res%initialConnectorId /= -1) then
            if (.not. isMaterialIdOfType(res%initialConnectorId, J_MAT_TYPE_CONNECTOR)) then
               write(error_unit, *) errorMsgInit, "material with id ", res%materialId, " must be a connector."
            end if
         end if 
         if (res%endConnectorId /= -1) then
            if (.not. isMaterialIdOfType(res%endConnectorId, J_MAT_TYPE_CONNECTOR)) then
               write(error_unit, *) errorMsgInit, "material with id ", res%materialId, " must be a connector."
            end if
         end if
      end if
      if (isMultiwire) then
         ! Not defininign a containedWithinElementId is an error if the simulation is a 3D-FDTD one. 
         ! For pure MTLN mode it is not an error.
         ! if (res%containedWithinElementId == -1) then
         !    write(error_unit, *) errorMsgInit, "multiwire associations must include: ", J_MAT_ASS_CAB_CONTAINED_WITHIN_ID
         ! end if
         if (.not. (this%getLogicalAt(this%root, J_GENERAL//'.'//J_GEN_MTLN_PROBLEM, default = .false.)) .and. &
            (this%mesh%checkElementId(res%containedWithinElementId) /= 0)) then
            write(error_unit, *) errorMsgInit, "element with id ", res%containedWithinElementId, " not found."
         end if
      end if
      
   contains 
      logical function isMaterialIdOfType(matId, matType)
         integer, intent(in) :: matId
         character (len=*), intent(in) :: matType
         type(json_value_ptr) :: mat
         logical :: materialFound
         if (this%matTable%checkId(matId) /= 0) then
            write(error_unit, *) "Material with id ", matId, " not found."
         end if
         mat = this%matTable%getId(matId)
         isMaterialIdOfType = this%getStrAt(mat%p, J_TYPE) == matType
      end function
      
      subroutine showLabelNotFoundError(label)
         character (len=*), intent(in) :: label
         
      end subroutine
   end function

   function getMaterialAssociations(this, materialTypes) result(res)
      class(parser_t) :: this
      character(len=*), intent(in) :: materialTypes(:)
      type(materialAssociation_t), dimension(:), allocatable :: res
      type(json_value), pointer :: allMatAss
      
      type(json_value), pointer :: mAPtr
      integer :: i, j, k
      integer :: nMaterials
      logical :: found

      call this%core%get(this%root, J_MATERIAL_ASSOCIATIONS, allMatAss, found)
      if (.not. found) then
         allocate(res(0))
         return
      end if

      nMaterials = 0
      do i = 1, this%core%count(allMatAss)
         call this%core%get_child(allMatAss, i, mAPtr)
         do j = 1, size(materialTypes)
            if (isAssociatedWithMaterial(mAPtr, trim(materialTypes(j)))) then
               nMaterials = nMaterials + 1
            end if
         end do
      end do

      allocate(res(nMaterials))
      j = 1
      do i = 1, this%core%count(allMatAss)
         call this%core%get_child(allMatAss, i, mAPtr)
         do k = 1, size(materialTypes)
            if (isAssociatedWithMaterial(mAPtr, trim(materialTypes(k)))) then
               res(j) = this%parseMaterialAssociation(mAPtr)
               j = j+1
            end if
         end do
      end do

   contains 
      logical function isAssociatedWithMaterial(mAPtr, materialType)
         type(json_value), pointer, intent(in) :: mAPtr
         character (len=*), intent(in) :: materialType
         
         type(materialAssociation_t) :: matAss
         type(json_value_ptr) :: mat

         matAss = this%parseMaterialAssociation(mAPtr)
         mat = this%matTable%getId(matAss%materialId)
         isAssociatedWithMaterial = this%getStrAt(mat%p, J_TYPE) == materialType
      end function
   end function

   function buildTagName(this, matId, elementId) result(res)
      class(parser_t) :: this
      integer, intent(in) :: matId, elementId
      character(len=BUFSIZE) :: res
      character(len=:), allocatable :: matName, layerName
      logical :: found
      
      block
         type(json_value_ptr) :: mat
         mat = this%matTable%getId(matId)
         matName = this%getStrAt(mat%p, J_NAME, found)
         if (.not. found) then
            deallocate(matName)
            allocate(character(len(TAG_MATERIAL) + 12) :: matName)
            write(matName, '(a,i0)') TAG_MATERIAL, matId
         end if
         matName = adaptName(matName)
      end block
      
      block
         type(json_value_ptr) :: elem
         elem = this%elementTable%getId(elementId)
         layerName = this%getStrAt(elem%p, J_NAME, found)
         if (.not. found) then
            deallocate(layerName)
            allocate(character(len(TAG_LAYER) + 12) :: layerName)
            write(layerName, '(a,i0)') TAG_LAYER, elementId 
         end if
         layerName = adaptName(layerName)
      end block
      
      call checkIsValidName(matName)
      call checkIsValidName(layerName)
      res = trim(matName // '@' // layerName)
   contains
      subroutine checkIsValidName(str)
         character (len=:), allocatable, intent(in) :: str
         character (len=*), parameter :: notAllowedChars = '@'
         integer :: i
         do i = 1, len((notAllowedChars))
            if (index(str, notAllowedChars(i:i)) /= 0) then
               write(error_unit, *) "ERROR in name: ", str, &
                  " contains invalid character ", notAllowedChars(i:i)
            end if 
         end do
      end subroutine

      function adaptName(str) result(res)
         character (len=:), allocatable, intent(in) :: str
         character (len=:), allocatable :: res
         integer :: i
         res = trim(adjustl(str))
         do i = 1, len(res)
            if (res(i:i) == ' ') then
               res(i:i) = '_'
            end if
         end do
      end function
   end function


#ifdef CompileWithMTLN
   function readMTLN(this, grid) result (mtln_res)
      class(parser_t) :: this
      type(Desplazamiento), intent(in) :: grid
      type(mtln_t) :: mtln_res
      type(fhash_tbl_t) :: elemIdToPosition, elemIdToCable, connIdToConnector
      type(materialAssociation_t), dimension(:), allocatable :: wires, multiwires, cables

      mtln_res%time_step = this%getRealAt(this%root, J_GENERAL//'.'//J_GEN_TIME_STEP)
      mtln_res%number_of_steps = this%getRealAt(this%root, J_GENERAL//'.'//J_GEN_NUMBER_OF_STEPS)

      wires = this%getMaterialAssociations([J_MAT_TYPE_WIRE])
      multiwires = this%getMaterialAssociations([J_MAT_TYPE_MULTIWIRE])
      cables = this%getMaterialAssociations(&
               [J_MAT_TYPE_WIRE//'     ',&
                J_MAT_TYPE_MULTIWIRE    ]) 
      ! 5 spaces are needed to make strings have same length. 
      ! Why? Because of FORTRAN! It only accepts fixed length strings for arrays.

      mtln_res%connectors => readConnectors()
      call addConnIdToConnectorMap(connIdToConnector, mtln_res%connectors)

      if (size(multiwires) /= 0) mtln_res%has_multiwires = .true.

      allocate (mtln_res%cables(size(wires) + size(multiwires)))
      block
         logical :: is_read
         integer :: i, j, ncc
         type(cable_t) :: read_cable
         ncc = 0
         do i = 1, size(cables)
            is_read = .true.
            read_cable = readMTLNCable(cables(i), is_read)
            if (.not. isCableNameUnique(read_cable, mtln_res%cables, ncc)) then
               error stop 'Cable name "'//read_cable%name//'" has already been used'
            end if
            ncc = ncc + 1
            mtln_res%cables(ncc) = read_cable
            call addElemIdToCableMap(elemIdToCable, cables(i)%elementIds, ncc)
            call addElemIdToPositionMap(elemIdToPosition, cables(i)%elementIds)
         end do
      end block

      block
         integer :: i, j, parentId, index
         type(json_value_ptr) :: mat
         j = 1
         do i = 1, size(cables)
            mat = this%matTable%getId(cables(i)%materialId)
            if (this%getStrAt(mat%p, J_TYPE) == J_MAT_TYPE_MULTIWIRE) then
               parentId = cables(i)%containedWithinElementId
               if (parentId == -1) then
                  mtln_res%cables(j)%parent_cable => null()
                  mtln_res%cables(j)%conductor_in_parent = 0
               else 
                  call elemIdToCable%get(key(parentId), value=index)
                  mtln_res%cables(j)%parent_cable => mtln_res%cables(index)
                  mtln_res%cables(j)%conductor_in_parent = getParentPositionInMultiwire(parentId)
               end if
            else if (this%getStrAt(mat%p, J_TYPE) == J_MAT_TYPE_WIRE) then 
               mtln_res%cables(j)%parent_cable => null()
               mtln_res%cables(j)%conductor_in_parent = 0
            else
               write(error_unit, *) 'ERROR: Material type not recognized'
            end if
            j = j + 1
         end do
      end block

      mtln_res%probes = readWireProbes()
      mtln_res%networks = buildNetworks()

   contains

      function isCableNameUnique(cable, cables, n) result(res)
         type(cable_t) :: cable
         type(cable_t), dimension(:), pointer :: cables
         integer :: n, i
         logical :: res
         res = .true.
         do i = 1, n
            if (cable%name == cables(i)%name) then
               res = .false.
               exit
            end if
         end do
      end function

      function readConnectors() result(res)
         type(connector_t), dimension(:), pointer :: res
         type(json_value), pointer :: mat, z

         type(json_value_ptr), dimension(:), allocatable :: connectors
         integer :: i, id
         call this%core%get(this%root, J_MATERIALS, mat)
         connectors = this%jsonValueFilterByKeyValue(mat, J_TYPE, J_MAT_TYPE_CONNECTOR)
         allocate(res(size(connectors)))
         if (size(connectors) /= 0) then
            do i = 1, size(connectors)
               res(i)%id = this%getIntAt(connectors(i)%p, J_ID)
               if (this%existsAt(connectors(i)%p, J_MAT_CONN_RESISTANCES)) then
                  res(i)%resistances = this%getRealsAt(connectors(i)%p, J_MAT_CONN_RESISTANCES)
               else
                  allocate(res(i)%resistances(0))
               end if

               if (this%existsAt(connectors(i)%p, J_MAT_CONN_TRANSFER_IMPEDANCE)) then
                  call this%core%get(connectors(i)%p, J_MAT_MULTIWIRE_TRANSFER_IMPEDANCE,  z)
                  res(i)%transfer_impedance_per_meter = readTransferImpedance(z)
               else
                  res(i)%transfer_impedance_per_meter = noTransferImpedance()
               end if

            end do
         end if

      end function

      function findMaxElemId(cables) result(res)
         type(json_value_ptr), dimension(:), intent(in) :: cables
         integer :: i, m
         integer, dimension(:), allocatable :: elemIds
         integer :: res
         res = 0

         if (size(cables) /= 0) then
            do i = 1, size(cables)
               elemIds = getCableElemIds(cables(i)%p)
               if (size(elemIds) == 0) return
               m = maxval(elemIds,dim = 1)
               if (m > res) then
                  res = m
               end if
            end do
         end if

      end function

      function buildNetworks() result(res)
         type(terminal_network_t), dimension(:), allocatable :: res
         type(aux_node_t), dimension(:), allocatable :: aux_nodes
         integer :: i,j
         integer, dimension(:), allocatable :: elemIds
         type(json_value), pointer :: terminations_ini, terminations_end
         type(coordinate_t), dimension(:), allocatable :: networks_coordinates
         type(subcircuit_t), dimension(:), allocatable :: subcircuits
         type(materialAssociation_t), dimension(:), allocatable :: cables
         subcircuits = readSubcircuits()
         
         allocate(aux_nodes(0))
         allocate(networks_coordinates(0))
         cables = [ this%getMaterialAssociations([J_MAT_TYPE_WIRE]), &
                     this%getMaterialAssociations([J_MAT_TYPE_MULTIWIRE]) ]
         do i = 1, size(cables)
            elemIds = cables(i)%elementIds
            terminations_ini => getTerminationsOnSide(cables(i)%initialTerminalId)
            terminations_end => getTerminationsOnSide(cables(i)%endTerminalId)
            do j = 1, size(elemIds)
               aux_nodes = [aux_nodes, buildNode(terminations_ini, TERMINAL_NODE_SIDE_INI, j, elemIds(j))]
               aux_nodes = [aux_nodes, buildNode(terminations_end, TERMINAL_NODE_SIDE_END, j, elemIds(j))]
               call updateListOfNetworksCoordinates(networks_coordinates, elemIds(j))
            end do

         end do
         allocate(res(size(networks_coordinates)))


         do i = 1, size(networks_coordinates)
            res(i) = buildNetwork(networks_coordinates(i), aux_nodes, subcircuits)
         end do

      end function

      function readSubcircuits() result(res_ckt)
         type(subcircuit_t), dimension(:), allocatable :: res_ckt
         type(json_value), pointer :: subCkt, ckt
         type(json_value_ptr) :: m
         integer :: i, j, id
         logical :: found
         type(coordinate_t) :: ports_coordinate
         type(node_t) :: node
         integer, dimension(:), allocatable :: elemIds

         if (this%existsAt(this%root,  J_SUBCIRCUITS)) then
            call this%core%get(this%root, J_SUBCIRCUITS, subCkt)
            allocate(res_ckt(this%core%count(subCkt)))
            do i = 1, this%core%count(subCkt)
               call this%core%get_child(subCkt, i, ckt)
               res_ckt(i)%subcircuit_name = trim(adjustl(this%getStrAt(ckt,J_SUBCKT_NAME)))
               
               elemIds = this%getIntsAt(ckt, J_ELEMENTIDS)
               node = this%mesh%getNode(elemIds(1))
               res_ckt(i)%nodeId = node%coordIds(1)

               id = this%getIntAt(ckt,J_MATERIAL_ID)
               m = this%matTable%getId(this%getIntAt(ckt, J_MATERIAL_ID, found))
               if (.not. found) &
                  write(error_unit, *) "Error reading material region: materialId label not found."
               res_ckt(i)%model_file = this%getStrAt(m%p, J_SUBCKT_FILE)
               res_ckt(i)%model_name = this%getStrAt(m%p, J_SUBCKT_NAME)
               res_ckt(i)%numberOfPorts = this%getIntAt(m%p, J_SUBCKT_PORTS)
            end do
            
         else
            allocate(res_ckt(0))
         end if


      end function

      function buildListOfCoordinates(elemIds) result(res)
         integer, dimension(:), intent(in) :: elemIds
         integer :: i
         type(coordinate_t), dimension(:), allocatable :: res
         allocate(res(0))
         do i = 1, size(elemIds)
            call updateListOfNetworksCoordinates(res, elemIds(i))
         end do
      end function

      function countSubcircuitsInNetwork(network_coordinate,subcircuits) result (res)
         type(coordinate_t) :: network_coordinate
         type(subcircuit_t), dimension(:), intent(in) :: subcircuits
         integer :: i, res
         res = 0
         do i = 1, size(subcircuits)
            if (network_coordinate == this%mesh%getCoordinate(subcircuits(i)%nodeId)) then 
               res = res + 1
            end if
         end do
      end function

      function buildNetwork(network_coordinate, aux_nodes, subcircuits) result(res)
         type(coordinate_t) :: network_coordinate
         type(aux_node_t), dimension(:), intent(in) :: aux_nodes
         type(subcircuit_t), dimension(:), intent(in) :: subcircuits

         type(aux_node_t), dimension(:), allocatable :: network_nodes
         integer, dimension(:), allocatable :: node_ids
         integer :: i, j
         type(terminal_network_t) :: res
         type(node_t) :: node
         
         network_nodes = filterNetworkNodes(network_coordinate, aux_nodes)
         node_ids = buildListOfNodeIds(network_nodes)

         do i = 1, size(node_ids)
            call res%add_connection(buildConnection(node_ids(i), network_nodes, subcircuits))
         end do


      end function

      function buildListOfNodeIds(network_nodes) result(res)
         type(aux_node_t), dimension(:), intent(in) :: network_nodes
         integer, dimension(:), allocatable :: res
         integer :: i
         allocate(res(0))
         do i = 1, size(network_nodes)
            if (findloc(res, network_nodes(i)%cId, 1) == 0) res = [res, network_nodes(i)%cId]
         end do
      end function

      function filterNetworkNodes(network_coordinate, aux_nodes) result(res)
         type(coordinate_t), intent(in) :: network_coordinate
         type(aux_node_t), dimension(:), intent(in) :: aux_nodes
         type(aux_node_t), dimension(:), allocatable :: res
         integer :: i
         allocate(res(0))
         do i = 1, size(aux_nodes)
            if (aux_nodes(i)%relPos == network_coordinate) then
               res = [res, aux_nodes(i)]
            end if
         end do
      end function

      function buildConnection(node_id, network_nodes, subcircuits) result (res)
         integer, intent(in) :: node_id
         type(aux_node_t), dimension(:), intent(in) :: network_nodes
         type(subcircuit_t), dimension(:), intent(in) :: subcircuits
         type(terminal_connection_t) :: res
         integer :: i
         do i = 1, size(network_nodes)
            if (network_nodes(i)%cId == node_id) then
               call res%add_node(network_nodes(i)%node)
            end if
         end do
         do i = 1, size(subcircuits)
            if (subcircuits(i)%nodeId == node_id) then
               res%subcircuit = subcircuits(i)
               res%has_subcircuit = .true.
            end if
         end do


      end function

      subroutine updateListOfConnectionIds(ids, id)
         integer, dimension(:), intent(inout) :: ids
         integer, intent(in) :: id
         if (findloc(ids, id, 1) == 0) ids = [ids, id]
      end subroutine

      subroutine updateListOfNetworksCoordinates(coordinates, conductor_index)
         type(coordinate_t), dimension(:), allocatable,  intent(inout) :: coordinates
         integer, intent(in) :: conductor_index
         type (polyline_t) ::polyline
         integer :: i
         logical :: found_ini, found_end
         type(coordinate_t) :: coord_ini, coord_end
         integer :: ub

         found_ini = .false.
         found_end = .false.
         polyline = this%mesh%getPolyline(conductor_index)
         coord_ini = this%mesh%getCoordinate(polyline%coordIds(1))

         ub = ubound(polyline%coordIds,1)
         coord_end = this%mesh%getCoordinate(polyline%coordIds(ub))

         if (size(coordinates) /= 0) then
            do i = 1, size(coordinates)
               if (coordinates(i) == coord_ini) then
                  found_ini = .true.
               end if
               if (coordinates(i) == coord_end) then
                  found_end = .true.
               end if
            end do
         end if

         if (.not. found_ini) then
            coordinates = [coordinates, coord_ini]
         end if

         if (.not. found_end) then
            coordinates = [coordinates, coord_end]
         end if


      end subroutine

      function getTerminationsOnSide(terminationId) result(res)
         integer, intent(in) :: terminationId
         type(json_value_ptr) :: terminal
         type(json_value), pointer :: res

         if (terminationId == -1) then
            write(error_unit, *) 'Error: missing terminal on cable side'
            res => null()
            return
         end if
         terminal = this%matTable%getId(terminationId)
         if (.not. this%existsAt(terminal%p, J_MAT_TERM_TERMINATIONS)) then
            write(error_unit, *) 'Error: missing terminations on terminal'
            res => null()
            return
         end if
         call this%core%get(terminal%p, J_MAT_TERM_TERMINATIONS, res)

      end function



      function buildNode(termination_list, label, index, id) result(res)
         type(json_value), pointer :: termination_list, termination
         integer, intent(in) :: label
         integer, intent(in) :: index, id
         type(polyline_t) :: polyline
         type(aux_node_t) :: res
         integer :: cable_index
         call this%core%get_child(termination_list, index, termination)
         
         res%node%termination%termination_type = readTerminationType(termination)
         res%node%termination%capacitance = readTerminationRLC(termination,J_MAT_TERM_CAPACITANCE, default = 1e22)
         res%node%termination%resistance = readTerminationRLC(termination, J_MAT_TERM_RESISTANCE, default = 0.0)
         res%node%termination%inductance = readTerminationRLC(termination, J_MAT_TERM_INDUCTANCE, default=0.0)
         res%node%termination%source = readGeneratorOnTermination(id,label)
         res%node%termination%model = readTerminationModel(termination)
         res%node%termination%subcircuitPort = readTerminationSubcircuitPort(termination, default = -1)
         
         res%node%side = label
         res%node%conductor_in_cable = index

         call elemIdToCable%get(key(id), value=cable_index)
         res%node%belongs_to_cable => mtln_res%cables(cable_index)

         polyline = this%mesh%getPolyline(id)

         if (label == TERMINAL_NODE_SIDE_INI) then
            res%cId = polyline%coordIds(1)
            res%relPos = this%mesh%getCoordinate(polyline%coordIds(1))
         else if (label == TERMINAL_NODE_SIDE_END) then
            res%cId = polyline%coordIds(ubound(polyline%coordIds,1))
            res%relPos = this%mesh%getCoordinate(polyline%coordIds(ubound(polyline%coordIds,1)))
         end if
      end function

      function readGeneratorOnTermination(id, label) result(res)
         integer, intent(in) :: id, label
         type(json_value), pointer :: sources
         type(json_value_ptr), dimension(:), allocatable :: genSrcs
         logical :: found
         type(node_source_t) :: res
         integer :: polylineId

         character (len=*), dimension(1), parameter :: validTypes = &
         [J_SRC_TYPE_GEN]

         call this%core%get(this%root, J_SOURCES, sources, found)
         if (.not. found) then
            res%path_to_excitation = trim("")
            res%source_type = SOURCE_TYPE_UNDEFINED
            return
         end if
         
         genSrcs = this%jsonValueFilterByKeyValues(sources, J_TYPE, validTypes)
         if (size(genSrcs) == 0) then
            res%path_to_excitation = trim("")
            res%source_type = SOURCE_TYPE_UNDEFINED
            return
         end if


         block
            type(node_t) :: srcCoord
            integer, dimension(:), allocatable :: sourceElemIds
            type(polyline_t) :: poly
            integer :: i
   
            poly = this%mesh%getPolyline(id)
            do i = 1, size(genSrcs)
               if (.not. this%existsAt(genSrcs(i)%p, J_SRC_MAGNITUDE_FILE)) then
                  write(error_unit, *) 'magnitudeFile of source missing'
                  res%path_to_excitation = trim("")
                  res%source_type = SOURCE_TYPE_UNDEFINED
                  return
               end if
               if (.not. this%existsAt(genSrcs(i)%p, J_FIELD)) then
                  write(error_unit, *) 'Type of generator is ambigous'
                  res%path_to_excitation = trim("")
                  res%source_type = SOURCE_TYPE_UNDEFINED
                  return
               end if
               if (this%getStrAt(genSrcs(i)%p, J_FIELD) /= J_FIELD_VOLTAGE .and. &
                   this%getStrAt(genSrcs(i)%p, J_FIELD) /= J_FIELD_CURRENT) then 
                  write(error_unit, *) 'Only voltage and current generators are supported'
                  res%path_to_excitation = trim("")
                  res%source_type = SOURCE_TYPE_UNDEFINED
                  return
               end if

               if (isSourceAttachedToLine(genSrcs(i)%p, poly, id, label)) then 
                  if (this%getStrAt(genSrcs(i)%p, J_FIELD) == J_FIELD_VOLTAGE) then 
                     res%source_type = SOURCE_TYPE_VOLTAGE 
                  else if (this%getStrAt(genSrcs(i)%p, J_FIELD) == J_FIELD_CURRENT) then 
                     res%source_type = SOURCE_TYPE_CURRENT
                  end if
                  res%path_to_excitation = trim(this%getStrAt(genSrcs(i)%p, J_SRC_MAGNITUDE_FILE))
                  return
               end if

            end do
            res%path_to_excitation = trim("")
            res%source_type = SOURCE_TYPE_UNDEFINED
         end block 
      end function

      function isSourceAttachedToLine(src, polyline, id, label) result(res)
         type(json_value), pointer, intent(in)  :: src
         type(polyline_t), intent(in) :: polyline
         integer, intent(in) :: id, label
         integer :: index
         integer, dimension(:), allocatable :: sourceElemIds
         type(node_t) :: srcCoord
         logical :: res

         sourceElemIds = this%getIntsAt(src, J_ELEMENTIDS)
         srcCoord = this%mesh%getNode(sourceElemIds(1))

         if (label == TERMINAL_NODE_SIDE_INI) then
            index = 1
         else if (label == TERMINAL_NODE_SIDE_END) then 
            index = ubound(polyline%coordIds,1)
         end if
         
         if (this%existsAt(src, J_SRC_ATTACHED_ID)) then 
            res = (srcCoord%coordIds(1) == polyline%coordIds(index)) .and. (this%getIntAt(src, J_SRC_ATTACHED_ID) == id)
         else
            res = (srcCoord%coordIds(1) == polyline%coordIds(index))
         end if
      
      end function

      function readTerminationType(termination) result(res)
         type(json_value), pointer :: termination
         integer :: res
         character(:), allocatable :: type
         type = this%getStrAt(termination, J_TYPE)
         if (type == J_MAT_TERM_TYPE_OPEN) then
            res = TERMINATION_OPEN
         else if (type == J_MAT_TERM_TYPE_SHORT) then
            res = TERMINATION_SHORT
         else if (type == J_MAT_TERM_TYPE_SERIES) then
            res = TERMINATION_SERIES
         else if (type == J_MAT_TERM_TYPE_LCpRs) then
            res = TERMINATION_LCpRs
         else if (type == J_MAT_TERM_TYPE_RLsCp) then
            res = TERMINATION_RLsCp
         else if (type == J_MAT_TERM_TYPE_CIRCUIT) then 
            res = TERMINATION_CIRCUIT
         else
            res = TERMINATION_UNDEFINED
         end if
      end function

      function readTerminationModel(termination) result(res)
         type(json_value), pointer :: termination
         type(terminal_circuit_t) :: res
         if (this%existsAt(termination, J_MAT_TERM_MODEL_FILE)) then
            res%file = this%getStrAt(termination, J_MAT_TERM_MODEL_FILE)
         else
            res%file = ""
         end if

         if (this%existsAt(termination, J_MAT_TERM_MODEL_NAME)) then
            res%model_name = this%getStrAt(termination, J_MAT_TERM_MODEL_NAME)
         else
            res%model_name = ""
         end if

      end function

      function readTerminationSubcircuitPort(termination, default) result(res)
         type(json_value), pointer :: termination
         integer, intent(in) :: default
         integer :: res
         if (this%existsAt(termination, J_MAT_TERM_MODEL_PORT)) then
            res = this%getIntAt(termination, J_MAT_TERM_MODEL_PORT)
         else
            res = default
         end if

      end function

      function readTerminationRLC(termination, label, default) result(res)
         type(json_value), pointer :: termination
         character(*), intent(in) :: label
         real, intent(in) :: default
         real :: res
         if (this%existsAt(termination, label)) then
            res = this%getRealAt(termination, label)
         else
            res = default
         end if

      end function

      function readWireProbes() result(res)
         type(probe_t), dimension(:), allocatable :: res
         type(json_value_ptr), dimension(:), allocatable :: wire_probes, polylines
         type(json_value), pointer :: probes, elements
         integer :: i,j, position, n_probes, k
         integer, dimension(:), allocatable :: elemIds, polylinecIds
         type(node_t) :: node
         type(cable_t), pointer :: cable_ptr
         integer :: index
         type(coordinate_t) :: node_coord
         if (this%existsAt(this%root, J_PROBES)) then
            call this%core%get(this%root, J_PROBES, probes)
         else
            allocate(res(0))
            return
         end if

         call this%core%get(this%root, J_MESH//'.'//J_ELEMENTS, elements)
         polylines = this%jsonValueFilterByKeyValue(elements, J_TYPE, J_ELEM_TYPE_POLYLINE)
         wire_probes = this%jsonValueFilterByKeyValue(probes, J_TYPE, J_MAT_TYPE_WIRE)

         n_probes = countProbes(wire_probes, polylines)

         allocate(res(n_probes))
         if (n_probes /= 0) then
            k = 1
            do i = 1, size(wire_probes)
               elemIds = this%getIntsAt(wire_probes(i)%p, J_ELEMENTIDS)
               node = this%mesh%getNode(elemIds(1))
               node_coord = this%mesh%getCoordinate(node%coordIds(1))
               do j = 1, size(polylines)
                  polylinecIds = this%getIntsAt(polylines(j)%p, J_COORDINATE_IDS)
                  position = findloc(polylinecIds, node%coordIds(1), dim=1)
                  if (position /= 0) then ! polyline found
                     res(k)%probe_type = readProbeType(wire_probes(i)%p)
                     res(k)%probe_name = readProbeName(wire_probes(i)%p)
                     res(k)%probe_position = node_coord%position
                     call elemIdToCable%get(key(this%getIntAt(polylines(j)%p, J_ID)), value=index)
                     cable_ptr => mtln_res%cables(index)
                     do while (associated(cable_ptr%parent_cable))
                        cable_ptr => cable_ptr%parent_cable
                     end do
                     res(k)%attached_to_cable => cable_ptr
                     res(k)%index = findProbeIndex(polylinecIds, position)
                     k = k + 1
                  end if
               end do
            end do
         end if


      end function

      function countProbes(probes, lines) result(res)
         type(json_value_ptr), dimension(:), allocatable :: probes, lines
         integer :: res
         integer, dimension(:), allocatable :: ids, lines_ids
         type(node_t) :: node
         integer :: i, j, position
         res = 0
         if (size(probes) /= 0) then
            do i = 1, size(probes)
               ids = this%getIntsAt(probes(i)%p, J_ELEMENTIDS)
               node = this%mesh%getNode(ids(1))
               do j = 1, size(lines)
                  lines_ids = this%getIntsAt(lines(j)%p, J_COORDINATE_IDS)
                  position = findloc(lines_ids, node%coordIds(1), dim=1)
                  if (position /= 0) then ! polyline found
                     res = res + 1
                  end if
               end do
            end do
         end if
      end function

      function readProbeType(probe) result(res)
         type(json_value), pointer :: probe
         character(:), allocatable :: probe_type
         integer :: res
         probe_type = this%getStrAt(probe, J_FIELD)
         if (probe_type == J_FIELD_VOLTAGE) then
            res = PROBE_TYPE_VOLTAGE
         else if (probe_type == J_FIELD_CURRENT) then
            res = PROBE_TYPE_CURRENT
         else
            write(error_unit,*) 'probe type '//probe_type//' not supported'
            res = PROBE_TYPE_UNDEFINED
         end if
      end function

      function readProbeName(probe) result(res)
         type(json_value), pointer :: probe
         character(:), allocatable :: res
         if (this%existsAt(probe, J_NAME)) then 
            res = this%getStrAt(probe, J_NAME)
         else 
            res = ""
         end if
      end function

      function findProbeIndex(polyline_cIds, node_position) result(res)
         integer, dimension(:), intent(in) :: polyline_cIds
         integer, intent(in) :: node_position
         integer :: k, res
         type(coordinate_t) :: c1, c2, delta
         res = 1
         do k=2, node_position
            c2 = this%mesh%getCoordinate(polyline_cIds(k))
            c1 = this%mesh%getCoordinate(polyline_cIds(k-1))
            delta = c2-c1
            res = res + abs(delta%position(findDirection(delta)))
         end do

      end function

      function getCableContainingElemId(id) result(res)
         integer, intent(in) :: id
         integer :: mStat
         class(*), pointer :: d
         type(cable_t), pointer :: res

         nullify(res)
         call elemIdToCable%check_key(key(id), mStat)
         if (mStat /= 0) then
            return
         end if

         call elemIdToCable%get_raw_ptr(key(id), d, mStat)
         if (mStat /= 0) then
            return
         end if
         select type(d)
          type is (cable_t)

            res => d
         end select
      end function

      function findConnectorWithId(conn_Id) result(res)
         integer, intent(in) :: conn_id
         integer :: conn_index
         type(connector_t), pointer :: res
         if (conn_id /= -1) then
            call connIdToConnector%get(key(conn_id), conn_index)
            res => mtln_res%connectors(conn_index)
         else
            res => null()
         end if
      end function

      function getConnectorWithIdFromMap(id) result(res)
         integer, intent(in) :: id
         integer :: mStat
         class(*), pointer :: d
         type(connector_t), pointer :: res

         nullify(res)
         call connIdToConnector%check_key(key(id), mStat)
         if (mStat /= 0) then
            res => null()
            return
         end if

         call connIdToConnector%get_raw_ptr(key(id), d, mStat)
         if (mStat /= 0) then
            res => null()
            return
         end if
         select type(d)
          type is (connector_t)
            res => d
         end select
      end function

      function getParentPositionInMultiwire(id) result(res)
         integer, intent(in) :: id
         integer :: mStat
         integer :: res

         call elemIdToPosition%check_key(key(id), mStat)
         if (mStat /= 0) then
            return
         end if
         call elemIdToPosition%get(key(id), value=res)
      end function

      subroutine addConnIdToConnectorMap(map, conn)
         type(fhash_tbl_t), intent(inout) :: map
         type(connector_t), dimension(:), intent(in) :: conn
         integer :: i
         if (size(conn) == 0) return
         do i = 1, size(conn)
            call map%set(key(conn(i)%id), i)
         end do
      end subroutine


      subroutine addElemIdToCableMap(map, elemIds, index)
         type(fhash_tbl_t), intent(inout) :: map
         integer, dimension(:), intent(in) :: elemIds
         integer :: index
         integer :: i
         do i = 1, size(elemIds)
            call map%set(key(elemIds(i)), index)
         end do
      end subroutine

      subroutine addElemIdToPositionMap(map, elemIds)
         type(fhash_tbl_t), intent(inout) :: map
         integer, dimension(:), allocatable, intent(in) :: elemIds
         integer :: i
         do i = 1, size(elemIds)
            call map%set(key(elemIds(i)), i)
         end do
      end subroutine

      function getCableElemIds(cable) result(res)
         type(json_value), pointer :: cable
         integer, dimension(:), allocatable :: res
         if (this%existsAt(cable, J_ELEMENTIDS)) then
            res = this%getIntsAt(cable, J_ELEMENTIDS)
         else
            allocate(res(0))
            write(error_unit,*) 'Error reading materialAssociation region: elementIds label not found'
         end if
      end function


      function readMTLNCable(j_cable, is_read) result(res)
         type(materialAssociation_t), intent(in) :: j_cable
         type(cable_t) :: res
         type(json_value_ptr) :: material
         logical, intent(inout) :: is_read
         integer :: nConductors
         logical :: found
         character(:), allocatable :: materialType

         res%name = j_cable%name

         res%step_size = buildStepSize(j_cable)
         res%external_field_segments = mapSegmentsToGridCoordinates(j_cable)

         material = this%matTable%getId(j_cable%materialId)
         materialType = this%getStrAt(material%p, J_TYPE)
         if (materialType == J_MAT_TYPE_WIRE) then
            call assignReferenceProperties(res, material)
            call assignExternalRadius(res, material)
            if (this%existsAt(material%p, J_MAT_WIRE_DIELECTRIC)) then
               call assignDielectricProperties(res, material)
            end if

            if (this%existsAt(material%p, J_MAT_WIRE_PASS)) then 
               res%isPassthrough = this%getLogicalAt(material%p, J_MAT_WIRE_PASS)
            end if

         else if (materialType == J_MAT_TYPE_MULTIWIRE) then
            call assignPULProperties(res, material, size(j_cable%elementIds))
         else
            write(error_unit, *) "Error reading cable: is neither wire nor multiwire"
         end if

         res%initial_connector => findConnectorWithId(j_cable%initialConnectorId)
         res%end_connector => findConnectorWithId(j_cable%endConnectorId)
         res%transfer_impedance = buildTransferImpedance(material)


      end function

      subroutine assignExternalRadius(res, mat)
         type(cable_t), intent(inout) :: res
         type(json_value_ptr) :: mat
         integer :: i

         if (this%existsAt(mat%p, J_MAT_WIRE_RADIUS)) then
            do i = 1, size(res%external_field_segments(:))
               res%external_field_segments(i)%radius = this%getRealAt(mat%p, J_MAT_WIRE_RADIUS)
            end do
         else
            write(error_unit, *) "Wire radius is missing"
         end if

      end subroutine

      subroutine assignDielectricProperties(res, mat)
         type(cable_t), intent(inout) :: res
         type(json_value_ptr) :: mat, diel
         type(json_value), pointer :: diel_ptr
         integer :: i

         call this%core%get(mat%p, J_MAT_WIRE_DIELECTRIC, diel_ptr)
         if (this%existsAt(diel_ptr, J_MAT_WIRE_DIELECTRIC_PERMITTIVITY) .and. & 
             this%existsAt(diel_ptr, J_MAT_WIRE_DIELECTRIC_RADIUS)) then

            do i = 1, size(res%external_field_segments(:))
               res%external_field_segments(i)%has_dielectric = .true.
               res%external_field_segments(i)%dielectric%relative_permittivity = this%getRealAt(diel_ptr, J_MAT_WIRE_DIELECTRIC_PERMITTIVITY)
               res%external_field_segments(i)%dielectric%radius = this%getRealAt(diel_ptr, J_MAT_WIRE_DIELECTRIC_RADIUS)
            end do
         else
            write(error_unit, *) "Dielectric permittivity and/of radius is missing"
         end if

      end subroutine

      function buildTransferImpedance(mat) result(res)
         type(json_value_ptr):: mat
         type(transfer_impedance_per_meter_t) :: res
         type(json_value), pointer :: z
         if (this%existsAt(mat%p, J_MAT_MULTIWIRE_TRANSFER_IMPEDANCE)) then
            call this%core%get(mat%p, J_MAT_MULTIWIRE_TRANSFER_IMPEDANCE,  z)
            res = readTransferImpedance(z)
         else
            res = noTransferImpedance()
         end if
      end function

      !needs correction
      function buildConnector(j_cable, side) result(res)
         type(json_value), pointer :: j_cable
         character(*), intent(in) :: side
         type(connector_t), pointer :: res
         type(connector_t), target :: res_conn
         type(json_value_ptr) :: conn
         type(json_value), pointer :: z
         type(json_value), pointer :: c_ptr

         logical :: found
         character(:), allocatable :: name, type
         integer :: id
         real, dimension(:), allocatable :: rs
         if (this%existsAt(j_cable, side)) then
            conn = this%matTable%getId(this%getIntAt(j_cable, side))

            if (this%existsAt(conn%p, J_MAT_CONN_RESISTANCES)) then
               res_conn%resistances = this%getRealsAt(conn%p, J_MAT_CONN_RESISTANCES)
            else
               allocate(res_conn%resistances(0))
               write(error_unit, *) "Error reading connector: no resistances label found"
            end if

            if (this%existsAt(conn%p, J_MAT_MULTIWIRE_TRANSFER_IMPEDANCE)) then
               call this%core%get(conn%p, J_MAT_MULTIWIRE_TRANSFER_IMPEDANCE,  z)
               res_conn%transfer_impedance_per_meter = readTransferImpedance(z)
            else
               res_conn%transfer_impedance_per_meter = noTransferImpedance()
               write(error_unit, *) "Error reading connector: no transferImpedancePerMeter label found"
            end if
            res => res_conn
         else
            res => null()
         end if
      end function

      subroutine assignReferenceProperties(res, mat)
         type(cable_t), intent(inout) :: res
         type(json_value_ptr) :: mat
         real, dimension(1,1) :: val
         allocate(res%capacitance_per_meter(1,1), source = 0.0)
         allocate(res%inductance_per_meter(1,1), source = 0.0)
         allocate(res%resistance_per_meter(1,1), source = 0.0)
         allocate(res%conductance_per_meter(1,1), source = 0.0)

         if (this%existsAt(mat%p, J_MAT_WIRE_REF_CAPACITANCE)) then
            res%capacitance_per_meter(1,1) = this%getRealAt(mat%p, J_MAT_WIRE_REF_CAPACITANCE)
         else
            write(error_unit, *) "Capacitance per meter will be assigned in module Wire_bundles_mtln"
            res%capacitance_per_meter(1,1) = 0.0
         end if
         
         if (this%existsAt(mat%p, J_MAT_WIRE_REF_INDUCTANCE)) then
            res%inductance_per_meter(1,1) = this%getRealAt(mat%p, J_MAT_WIRE_REF_INDUCTANCE)
         else
            write(error_unit, *) "Inductance per meter will be assigned in module Wire_bundles_mtln"
            res%inductance_per_meter(1,1) = 0.0
         end if

         if (this%existsAt(mat%p, J_MAT_WIRE_RESISTANCE)) then
            res%resistance_per_meter(1,1) = this%getRealAt(mat%p, J_MAT_WIRE_RESISTANCE)
         else
            res%resistance_per_meter(1,1) = 0.0
         end if

      end subroutine

      subroutine assignPULProperties(res, mat, n)
         type(cable_t), intent(inout) :: res
         type(json_value_ptr) :: mat
         integer, intent(in) :: n
         real, dimension(:,:), allocatable :: null_matrix
         logical :: found
         allocate(null_matrix(n,n), source = 0.0)
         if (this%existsAt(mat%p, J_MAT_MULTIWIRE_INDUCTANCE)) then
            res%inductance_per_meter = this%getMatrixAt(mat%p, J_MAT_MULTIWIRE_INDUCTANCE,found)
         else
            write(error_unit, *) "Error reading material region: inductancePerMeter label not found."
            res%inductance_per_meter = null_matrix
         end if

         if (this%existsAt(mat%p, J_MAT_MULTIWIRE_CAPACITANCE)) then
            res%capacitance_per_meter = this%getMatrixAt(mat%p, J_MAT_MULTIWIRE_CAPACITANCE,found)
         else
            write(error_unit, *) "Error reading material region: capacitancePerMeter label not found."
            res%capacitance_per_meter = null_matrix
         end if

         if (this%existsAt(mat%p, J_MAT_MULTIWIRE_RESISTANCE)) then
            res%resistance_per_meter = vectorToDiagonalMatrix(this%getRealsAt(mat%p, J_MAT_MULTIWIRE_RESISTANCE,found))
         else
            res%resistance_per_meter = null_matrix
         end if

         if (this%existsAt(mat%p, J_MAT_MULTIWIRE_CONDUCTANCE)) then
            res%conductance_per_meter = vectorToDiagonalMatrix(this%getRealsAt(mat%p, J_MAT_MULTIWIRE_CONDUCTANCE,found))
         else
            res%conductance_per_meter = null_matrix
         end if


      end subroutine

      function mapSegmentsToGridCoordinates(j_cable) result(res)
         type(materialAssociation_t), intent(in) :: j_cable
         type(external_field_segment_t), dimension(:), allocatable :: res
         integer, dimension(:), allocatable :: elemIds
         type(polyline_t) :: p_line

         elemIds = j_cable%elementIds
         if (size(elemIds) == 0) return

         p_line = this%mesh%getPolyline(elemIds(1))
         allocate(res(0))
         block
            type(coordinate_t) :: c1, c2
            integer :: i
            do i = 2, size(p_line%coordIds)
               c2 = this%mesh%getCoordinate(p_line%coordIds(i))
               c1 = this%mesh%getCoordinate(p_line%coordIds(i-1))

               if (findOrientation(c2-c1) > 0) then
                  res = [res, mapPositiveSegment(c1,c2)]
               else if (findOrientation(c2-c1) < 0) then
                  res = [res, mapNegativeSegment(c1,c2)]
               else
                  write(error_unit, *) 'Error: polyline first and last coordinate are identical'
               end if
            end do
         end block
      end function

      function mapNegativeSegment(c1, c2) result(res)
         type(coordinate_t), intent(in) :: c1, c2
         type(external_field_segment_t) :: curr_pos
         integer :: axis, i, n_segments
         type(external_field_segment_t), dimension(:), allocatable :: res

         axis = findDirection(c2-c1)
         n_segments = abs(ceiling(c2%position(axis)) - floor(c1%position(axis)))
         allocate(res(n_segments))
         curr_pos%position = [(c1%position(i), i = 1, 3)]
         curr_pos%field => null()
         curr_pos%radius = 0.0

         res = [(curr_pos, i = 1, n_segments)]
         res(:)%position(axis) = [(res(i)%position(axis) - i, i = 1, n_segments)]
         res(:)%direction = -axis
      end function

      function mapPositiveSegment(c1, c2) result(res)
         type(coordinate_t), intent(in) :: c1, c2
         type(external_field_segment_t) :: curr_pos
         integer :: axis, orientation, i, n_segments
         type(external_field_segment_t), dimension(:), allocatable :: res

         axis = findDirection(c2-c1)

         n_segments = abs(floor(c2%position(axis)) - ceiling(c1%position(axis)))
         allocate(res(n_segments))
         curr_pos%position = [(c1%position(i), i = 1, 3)]
         curr_pos%field => null()

         res = [(curr_pos, i = 1, n_segments)]
         res(:)%position(axis) = [(res(i)%position(axis) + (i-1), i = 1, n_segments)]
         res(:)%direction = axis
      end function

      function buildStepSize(j_cable) result(res)
         type(materialAssociation_t), intent(in) :: j_cable
         real, dimension(:), allocatable :: res
         integer, dimension(:), allocatable :: elemIds
         type(polyline_t) :: p_line
         type(Desplazamiento) :: desp

         desp = this%readGrid()

         elemIds = j_cable%elementIds
         if (size(elemIds) == 0) return

         p_line = this%mesh%getPolyline(elemIds(1))
         allocate(res(0))
         block
            type(coordinate_t) :: c1, c2
            integer :: axis, i, j
            integer :: index_1, index_2
            real :: f1, f2
            real, dimension(:), allocatable :: displacement
            do j = 2, size(p_line%coordIds)
               c2 = this%mesh%getCoordinate(p_line%coordIds(j))
               c1 = this%mesh%getCoordinate(p_line%coordIds(j-1))
               axis = findDirection(c2-c1)
               f1 = abs(ceiling(c1%position(axis))-c1%position(axis))
               f2 = abs(c2%position(axis)-floor(c2%position(axis)))
               displacement = assignDisplacement(desp, axis)
               if (f1 /= 0) then
                  res = [res, f1*displacement(floor(c1%position(axis)))]
               end if
               index_1 = ceiling(min(abs(c1%position(axis)), abs(c2%position(axis))))
               index_2 = floor(max(abs(c1%position(axis)), abs(c2%position(axis))))
               do i = 1, index_2 - index_1
                  res = [res, displacement(i)]
               enddo
               if (f2 /= 0) then
                  res = [res, f2*displacement(floor(c2%position(axis)))]
               end if
            end do
         end block
      end function

      function readTransferImpedance(z) result(res)
         type(json_value), pointer :: z
         type(transfer_impedance_per_meter_t) :: res
         character(len=:), allocatable :: direction
         if (this%existsAt(z, J_MAT_TRANSFER_IMPEDANCE_RESISTANCE)) then
            res%resistive_term = this%getRealAt(z,J_MAT_TRANSFER_IMPEDANCE_RESISTANCE)
         end if

         if (this%existsAt(z, J_MAT_TRANSFER_IMPEDANCE_INDUCTANCE)) then
            res%inductive_term = this%getRealAt(z,J_MAT_TRANSFER_IMPEDANCE_INDUCTANCE)
         end if

         if (this%existsAt(z, J_MAT_TRANSFER_IMPEDANCE_DIRECTION)) then
            direction = trim(adjustl(this%getStrAt(z,J_MAT_TRANSFER_IMPEDANCE_DIRECTION)))
         else
            write(error_unit,*) 'Error reading material: direction of transferImpedancePerMeter missing'
         end if

         if (direction == "inwards") then
            res%direction = TRANSFER_IMPEDANCE_DIRECTION_INWARDS
         else if (direction == "outwards") then
            res%direction = TRANSFER_IMPEDANCE_DIRECTION_OUTWARDS
         else if (direction == "both") then
            res%direction = TRANSFER_IMPEDANCE_DIRECTION_BOTH
         end if
         
         block
            type(json_value), pointer :: poles, residues, p, r
            integer :: i, n
            real :: read_value
            if (this%existsAt(z, J_MAT_TRANSFER_IMPEDANCE_POLES)) then 
               n = this%getIntAt(z, J_MAT_TRANSFER_IMPEDANCE_NUMBER_POLES)
               allocate(res%poles(n),res%residues(n))
               call this%core%get(z, J_MAT_TRANSFER_IMPEDANCE_POLES, poles)
               call this%core%get(z, J_MAT_TRANSFER_IMPEDANCE_RESIDUES, residues)
               do i = 1, n 
                  call this%core%get_child(poles, i, p)
                  call this%core%get(p, '(1)', read_value)
                  res%poles(i)%re = read_value
                  call this%core%get(p, '(2)', read_value)
                  res%poles(i)%im = read_value

                  call this%core%get_child(residues, i, r)
                  call this%core%get(r, '(1)', read_value)
                  res%residues(i)%re = read_value
                  call this%core%get(r, '(2)', read_value)
                  res%residues(i)%im = read_value
               end do
            else 
               allocate(res%poles(0),res%residues(0))
            end if
               
         end block
      end function

      function noTransferImpedance() result(res)
         type(transfer_impedance_per_meter_t) :: res
         character(len=:), allocatable :: direction
         ! res%resistive_term = 0.0
         ! res%inductive_term = 0.0
         ! res%direction = 0
         allocate(res%poles(0), res%residues(0))
      end function


      function assignDisplacement(desp, axis) result (res)
         type(Desplazamiento), intent(in) :: desp
         integer, intent(in) :: axis
         real, dimension(:), allocatable :: res

         if (axis == 1) then
            allocate(res(size(desp%desX)))
            res = desp%desX
         else if (axis == 2) then
            allocate(res(size(desp%desY)))
            res = desp%desY
         else if (axis == 3) then
            allocate(res(size(desp%desZ)))
            res = desp%desZ
         end if
      end function

      function findOrientation(coordDiference) result(res)
         type(coordinate_t), intent(in) :: coordDiference
         integer :: res
         integer :: i
         do i = 1, 3
            if (coordDiference%position(i) /= 0) then
               res = coordDiference%position(i)/abs(coordDiference%position(i))
            end if
         end do
      end function


      function findDirection(coordDiference) result(res)
         type(coordinate_t), intent(in) :: coordDiference
         integer :: res
         integer :: i
         do i = 1, 3
            if (coordDiference%position(i) /= 0) res = i
         end do
      end function


   end function
#endif

   subroutine handleFoundAndDefault(path, found, defaultPresent)
      character(len=*), intent(in) :: path
      logical, intent(in) :: found
      logical, intent(in) :: defaultPresent
      if (.not. found .and. .not. defaultPresent) then
         write(error_unit,*) 'ERROR expecting a value at: '//path
      end if
   end subroutine

   function getLogicalAt(this, place, path, found, default) result(res)
      logical :: res
      class(parser_t) :: this
      type(json_value), pointer :: place
      character(len=*) :: path
      logical, intent(out), optional :: found
      logical, optional :: default
      logical :: localFound
      
      call this%core%get(place, path, res, localFound, default)
      if (present(found)) then
         found = localFound
      else
         call handleFoundAndDefault(path, localFound, present(default))
      endif
   end function


   function getIntAt(this, place, path, found, default) result(res)
      integer :: res
      class(parser_t) :: this
      type(json_value), pointer :: place
      character(len=*) :: path
      logical, intent(out), optional :: found
      integer, optional :: default
      logical :: localFound
      
      call this%core%get(place, path, res, localFound, default)
      if (present(found)) then
         found = localFound
      else
         call handleFoundAndDefault(path, localFound, present(default))
      endif
   end function

   function getIntsAt(this, place, path, found) result(res)
      integer, dimension(:), allocatable :: res
      class(parser_t) :: this
      type(json_value), pointer :: place
      character(len=*) :: path
      logical, intent(out), optional :: found
      logical :: localFound
      
      call this%core%get(place, path, res, localFound)
      if (present(found)) then
         found = localFound
      else
         call handleFoundAndDefault(path, localFound, .false.)
      endif
   end function

   function getRealAt(this, place, path, found, default) result(res)
      real :: res
      class(parser_t) :: this
      type(json_value), pointer :: place
      character(len=*) :: path
      logical, intent(out), optional :: found
      real, optional :: default
      logical :: localFound

      call this%core%get(place, path, res, localFound, default)
      if (present(found)) then
         found = localFound
      else
         call handleFoundAndDefault(path, localFound, present(default))
      endif
   end function

   function getRealsAt(this, place, path, found) result(res)
      real, dimension(:), allocatable :: res
      class(parser_t) :: this
      type(json_value), pointer :: place
      character(len=*) :: path
      logical, intent(out), optional :: found
      logical :: localFound

      call this%core%get(place, path, res, localfound)
      if (present(found)) then
         found = localFound
      else
         call handleFoundAndDefault(path, localFound, .false.)
      endif
   end function

   function getMatrixAt(this, place, path, found) result(res)
      real, dimension(:,:), allocatable :: res
      class(parser_t) :: this
      type(json_value), pointer :: place, matrix, row
      character(len=*) :: path
      logical, intent(out), optional :: found
      integer :: i, vartype, nr
      real, dimension(:), allocatable :: res_row
      logical :: localFound

      call this%core%get(place, path,  matrix, localfound)
      if (present(found)) then
         found = localFound
      else
         call handleFoundAndDefault(path, localFound, .false.)
      endif
      call this%core%info(matrix, vartype, nr)
      allocate(res(nr,nr))

      do i = 1, nr
         call this%core%get_child(matrix, i, row)
         call this%core%get(row, res_row)
         res(i,:) = res_row
      end do
      ! need to check if not found
   end function


   function getStrAt(this, place, path, found, default) result(res)
      character (len=:), allocatable :: res
      class(parser_t) :: this
      type(json_value), pointer :: place
      character(len=*) :: path
      logical, intent(out), optional :: found
      character (len=*), optional :: default
      logical :: localFound
      
      call this%core%get(place, path, res, localFound, default)
      if (present(found)) then
         found = localFound
      else
         call handleFoundAndDefault(path, localFound, present(default))
      endif
   end function

   function existsAt(this, place, path) result(res)
      logical :: res
      class(parser_t) :: this
      type(json_value), pointer :: place
      character(len=*) :: path
      call this%core%info(place, path, found=res)
   end function

   function jsonValueFilterByKeyValues(this, srcs, key, values) result (res)
      class(parser_t) :: this
      type(json_value_ptr), dimension(:), allocatable :: res
      type(json_value), pointer :: srcs

      character (kind=JSON_CK, len=*) :: key
      character (kind=JSON_CK, len=*), dimension(:) :: values

      type(json_value_ptr), dimension (:), allocatable :: foundEntries
      integer :: i, lastEntry, nEntries

      allocate(res(0))
      do i = 1, size(values)
         foundEntries = this%jsonValueFilterByKeyValue(srcs, key, values(i))
         if (size(foundEntries) /= 0) then
            res = [res, foundEntries]
         end if
      end do
   end function

   function jsonValueFilterByKeyValue(this, place, key, value) result (res)
      class(parser_t) :: this
      type(json_value_ptr), allocatable :: res(:)
      character (kind=JSON_CK, len=*) :: key, value
      type(json_value), pointer :: place, src
      character (kind=JSON_CK, len=:), allocatable :: typeStr
      integer :: i, j, n
      logical :: found

      n = 0
      do i = 1, this%core%count(place)
         call this%core%get_child(place, i, src)
         typeStr = this%getStrAt(src, key, found)
         call this%core%get(src, key, typeStr, found)
         if(found .and. typeStr == trim(value)) then
            n = n + 1
         end if
      end do

      allocate(res(n))
      j = 1
      do i = 1, this%core%count(place)
         call this%core%get_child(place, i, src)
         typeStr = this%getStrAt(src, key, found)
         if(found .and. typeStr == value) then
            res(j)%p => src
            j = j + 1
         end if
      end do
   end function

   function getSingleVolumeInElementsIds(this, pw) result (res)
      class(parser_t) :: this
      type(json_value), pointer :: pw
      type(cell_region_t) :: cellRegion
      integer, dimension(:), allocatable :: elemIds
      type(cell_interval_t), dimension(:), allocatable :: res
      logical :: found

      elemIds = this%getIntsAt(pw, J_ELEMENTIDS, found=found)
      if (.not. found) then
         write(error_unit, *) "Error reading single volume elementIds label not found."
      end if
      if (size(elemIds) /= 1) then
         write(error_unit, *) "Entity must contain a single elementId."
      end if
      cellRegion = this%mesh%getCellRegion(elemIds(1), found)
      if (.not. found) then
         write(error_unit, *) "Entity elementId ", elemIds(1), " not found."
      end if
      res = cellRegion%getIntervalsOfType(CELL_TYPE_VOXEL)
      if (size(res) /= 1) &
         write(error_unit, *) "Entity must contain a single cell region defining a volume."
   end function

#endif   
end module