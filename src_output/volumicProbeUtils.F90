module mod_volumicProbeUtils
   use FDETYPES
   USE mod_UTILS
   use outputTypes
   use mod_outputUtils
   implicit none
   private

   ! Public interface
   public :: find_and_store_important_coords
   public :: isValidPointForCurrent
   public :: isValidPointForField
   public :: createUnstructuredDataForVTU

   abstract interface
      logical function logical_func(component, i, j, k, problemInfo)
         import :: problem_info_t, SINGLE
         type(problem_info_t), intent(in) :: problemInfo
         integer(kind=SINGLE), intent(in) :: component, i, j, k
      end function logical_func
   end interface

   interface registerNode
      module procedure &
         registerNodeByIndex, &
         registerNodeByCoordinate
   end interface

contains

   subroutine find_and_store_important_coords(lowerBound, upperBound, component, problemInfo, nPoints, coords)
      type(cell_coordinate_t), intent(in)             :: lowerBound, upperBound
      integer(kind=SINGLE), intent(in)                :: component
      type(problem_info_t), intent(in)                :: problemInfo
      integer(kind=SINGLE), intent(out)               :: nPoints
      integer(kind=SINGLE), allocatable, intent(inout) :: coords(:, :)

      call count_required_coords(lowerBound, upperBound, component, problemInfo, nPoints)
      call alloc_and_init(coords, 3, nPoints, 0_SINGLE)
      call store_required_coords(lowerBound, upperBound, component, problemInfo, coords)
   end subroutine find_and_store_important_coords

   subroutine count_required_coords(lowerBound, upperBound, requestComponent, problemInfo, count)
      type(cell_coordinate_t), intent(in) :: lowerBound, upperBound
      integer(kind=SINGLE), intent(in)    :: requestComponent
      type(problem_info_t), intent(in)    :: problemInfo
      integer(kind=SINGLE), intent(out)   :: count

      integer :: i, j, k
      procedure(logical_func), pointer :: checker => null()
      integer :: component

      call get_checker_and_component(requestComponent, checker, component)

      count = 0
      do k = lowerBound%z, upperBound%z
      do j = lowerBound%y, upperBound%y
      do i = lowerBound%x, upperBound%x
         if (checker(component, i, j, k, problemInfo)) count = count + 1
      end do
      end do
      end do
   end subroutine count_required_coords

   subroutine store_required_coords(lowerBound, upperBound, requestComponent, problemInfo, coords)
      type(cell_coordinate_t), intent(in) :: lowerBound, upperBound
      integer(kind=SINGLE), intent(in)    :: requestComponent
      type(problem_info_t), intent(in)    :: problemInfo
      integer(kind=SINGLE), intent(inout) :: coords(:, :)

      integer :: i, j, k, count
      procedure(logical_func), pointer :: checker => null()
      integer :: component

      call get_checker_and_component(requestComponent, checker, component)

      count = 0
      do k = lowerBound%z, upperBound%z
      do j = lowerBound%y, upperBound%y
      do i = lowerBound%x, upperBound%x
         if (checker(component, i, j, k, problemInfo)) then
            count = count + 1
            coords(1, count) = i
            coords(2, count) = j
            coords(3, count) = k
         end if
      end do
      end do
      end do
   end subroutine store_required_coords

   subroutine createUnstructuredDataForVTU(counter, coords, currentType, Nodes, Edges, Quads, numNodes, numEdges, numQuads, usevtkindex, realXGrid, realYGrid, realZGrid)
      integer, intent(in) :: counter
      integer(kind=SINGLE), intent(in) :: coords(:, :), currentType(:)
      logical, intent(in) :: usevtkindex
      real(KIND=RKIND), pointer, dimension(:), intent(in) :: realXGrid, realYGrid, realZGrid

      integer(kind=4), intent(out):: numNodes, numQuads, numEdges
      real(kind=RKIND), allocatable, dimension(:, :), intent(out) :: Nodes
      integer(kind=4), allocatable, dimension(:, :), intent(out) ::  Edges, Quads

      if (counter == 0) return

      call countElements(counter, currentType, numEdges, numQuads)

      allocate (Edges(2, numEdges))
      allocate (Quads(4, numQuads))
      allocate (Nodes(3, 2*numEdges + 4*numQuads))

      call registerElements(counter, coords, currentType, Nodes, Edges, Quads, usevtkindex, realXGrid, realYGrid, realZGrid)
      return
   end subroutine

   subroutine registerNodeByIndex(nodes, nodeIdx, x, y, z)
      real(kind=RKIND), dimension(:, :), intent(inout) :: nodes
      integer(kind=SINGLE), intent(in) :: nodeIdx, x, y, z
      !We need to avoid using idx 0
      nodes(1, nodeIdx + 1) = x*1.0_RKIND
      nodes(2, nodeIdx + 1) = y*1.0_RKIND
      nodes(3, nodeIdx + 1) = z*1.0_RKIND
   end subroutine

   subroutine registerNodeByCoordinate(nodes, nodeIdx, x, y, z)
      real(kind=RKIND), dimension(:, :), intent(inout) :: nodes
      integer(kind=SINGLE), intent(in) :: nodeIdx
      real(kind=RKIND), intent(in) :: x, y, z
      !We need to avoid using idx 0
      nodes(1, nodeIdx + 1) = x
      nodes(2, nodeIdx + 1) = y
      nodes(3, nodeIdx + 1) = z
   end subroutine

   subroutine registerEdge(edges, edgeIdx, startNodeIdx, endNodeIdx)
      integer(kind=SINGLE), dimension(:, :), intent(inout) :: edges
      integer(kind=SINGLE), intent(in) :: edgeIdx, startNodeIdx, endNodeIdx

      edges(1, edgeIdx + 1) = startNodeIdx
      edges(2, edgeIdx + 1) = endNodeIdx
   end subroutine

   subroutine registerQuad(quads, quadIdx, firstNodeIdx, secondNodeIdx, thirdNodeIdx, fourthNodeIdx)
      integer(kind=SINGLE), dimension(:, :), intent(inout) :: quads
      integer(kind=SINGLE), intent(in) :: quadIdx, firstNodeIdx, secondNodeIdx, thirdNodeIdx, fourthNodeIdx

      quads(1, quadIdx + 1) = firstNodeIdx
      quads(2, quadIdx + 1) = secondNodeIdx
      quads(3, quadIdx + 1) = thirdNodeIdx
      quads(4, quadIdx + 1) = fourthNodeIdx
   end subroutine

   subroutine countElements(counter, currentType, numEdges, numQuads)
      integer, intent(in) :: counter
      integer(kind=SINGLE), intent(in) :: currentType(:)
      integer(kind=4), intent(out) :: numEdges, numQuads
      integer :: i

      numEdges = 0
      numQuads = 0

      do i = 1, counter
         if ((currentType(i) == iJx) .or. (currentType(i) == iJy) .or. (currentType(i) == iJz)) numEdges = numEdges + 1
    if ((currentType(i) == iBloqueJx) .or. (currentType(i) == iBloqueJy) .or. (currentType(i) == iBloqueJz)) numQuads = numQuads + 1
      end do
   end subroutine

   subroutine registerElements(counter, coords, currentType, Nodes, Edges, Quads, usevtkindex, realXGrid, realYGrid, realZGrid)
      integer, intent(in) :: counter
      integer(kind=SINGLE), intent(in) :: coords(:, :), currentType(:)
      real(kind=RKIND), intent(inout) :: Nodes(:, :)
      integer(kind=4), intent(inout) :: Edges(:, :), Quads(:, :)
      logical :: usevtkindex
      real(KIND=RKIND), pointer, dimension(:), intent(in) :: realXGrid, realYGrid, realZGrid

      integer :: nodeIdx, quadIdx, edgeIdx
      integer :: xCoord, yCoord, zCoord
      integer :: i

      nodeIdx = -1
      quadIdx = -1
      edgeIdx = -1
      
      do i = 1, counter
         xCoord = coords(1, i)
         yCoord = coords(2, i)
         zCoord = coords(3, i)

         select case (currentType(i))
         case (iJx)
            nodeIdx = nodeIdx + 2
            if (usevtkindex) then
               call registerNode(Nodes, nodeIdx - 1, xCoord    , yCoord, zCoord    )
               call registerNode(Nodes, nodeIdx    , xCoord + 1, yCoord, zCoord    )
            else
               call registerNode(Nodes, nodeIdx - 1, realXGrid(xCoord)    , realYGrid(yCoord), realZGrid(zCoord)    )
               call registerNode(Nodes, nodeIdx    , realXGrid(xCoord + 1), realYGrid(yCoord), realZGrid(zCoord)    )
            endif
            edgeIdx = edgeIdx + 1
            call registerEdge(Edges, edgeIdx, nodeIdx - 1, nodeIdx)

         case (iJy)
            nodeIdx = nodeIdx + 2
            if (usevtkindex) then
               call registerNode(Nodes, nodeIdx - 1, xCoord    , yCoord    , zCoord    )
               call registerNode(Nodes, nodeIdx    , xCoord    , yCoord + 1, zCoord    )
            else
               call registerNode(Nodes, nodeIdx - 1, realXGrid(xCoord)    , realYGrid(yCoord)    , realZGrid(zCoord)    )
               call registerNode(Nodes, nodeIdx    , realXGrid(xCoord)    , realYGrid(yCoord + 1), realZGrid(zCoord)    )
            endif
            edgeIdx = edgeIdx + 1
            call registerEdge(Edges, edgeIdx, nodeIdx - 1, nodeIdx)

         case (iJz)
            nodeIdx = nodeIdx + 2
            if (usevtkindex) then
               call registerNode(Nodes, nodeIdx - 1, xCoord, yCoord    , zCoord    )
               call registerNode(Nodes, nodeIdx    , xCoord, yCoord    , zCoord + 1)
            else
               call registerNode(Nodes, nodeIdx - 1, realXGrid(xCoord)    , realYGrid(yCoord)    , realZGrid(zCoord)    )
               call registerNode(Nodes, nodeIdx    , realXGrid(xCoord)    , realYGrid(yCoord)    , realZGrid(zCoord + 1))
            endif
            edgeIdx = edgeIdx + 1
            call registerEdge(Edges, edgeIdx, nodeIdx - 1, nodeIdx)

         case (iBloqueJx)
            nodeIdx = nodeIdx + 4
            if (usevtkindex) then
               call registerNode(Nodes, nodeIdx - 3, xCoord, yCoord    , zCoord    )
               call registerNode(Nodes, nodeIdx - 2, xCoord, yCoord + 1, zCoord    )
               call registerNode(Nodes, nodeIdx - 1, xCoord, yCoord + 1, zCoord + 1)
               call registerNode(Nodes, nodeIdx    , xCoord, yCoord    , zCoord + 1)
            else
               call registerNode(Nodes, nodeIdx - 3, realXGrid(xCoord), realYGrid(yCoord)    , realZGrid(zCoord)    )
               call registerNode(Nodes, nodeIdx - 2, realXGrid(xCoord), realYGrid(yCoord + 1), realZGrid(zCoord)    )
               call registerNode(Nodes, nodeIdx - 1, realXGrid(xCoord), realYGrid(yCoord + 1), realZGrid(zCoord + 1))
               call registerNode(Nodes, nodeIdx    , realXGrid(xCoord), realYGrid(yCoord)    , realZGrid(zCoord + 1))
            endif
            quadIdx = quadIdx + 1
            call registerQuad(Quads, quadIdx, nodeIdx - 3, nodeIdx - 2, nodeIdx - 1, nodeIdx)

         case (iBloqueJy)
            nodeIdx = nodeIdx + 4
            if (usevtkindex) then
               call registerNode(Nodes, nodeIdx - 3, xCoord    , yCoord, zCoord    )
               call registerNode(Nodes, nodeIdx - 2, xCoord + 1, yCoord, zCoord    )
               call registerNode(Nodes, nodeIdx - 1, xCoord + 1, yCoord, zCoord + 1)
               call registerNode(Nodes, nodeIdx    , xCoord    , yCoord, zCoord + 1)
            else
               call registerNode(Nodes, nodeIdx - 3, realXGrid(xCoord)    , realYGrid(yCoord), realZGrid(zCoord)    )
               call registerNode(Nodes, nodeIdx - 2, realXGrid(xCoord + 1), realYGrid(yCoord), realZGrid(zCoord)    )
               call registerNode(Nodes, nodeIdx - 1, realXGrid(xCoord + 1), realYGrid(yCoord), realZGrid(zCoord + 1))
               call registerNode(Nodes, nodeIdx    , realXGrid(xCoord)    , realYGrid(yCoord), realZGrid(zCoord + 1))
            endif
            quadIdx = quadIdx + 1
            call registerQuad(Quads, quadIdx, nodeIdx - 3, nodeIdx - 2, nodeIdx - 1, nodeIdx)

         case (iBloqueJz)
            nodeIdx = nodeIdx + 4
            if (usevtkindex) then 
               call registerNode(Nodes, nodeIdx - 3, xCoord    , yCoord    , zCoord)
               call registerNode(Nodes, nodeIdx - 2, xCoord + 1, yCoord    , zCoord)
               call registerNode(Nodes, nodeIdx - 1, xCoord + 1, yCoord + 1, zCoord)
               call registerNode(Nodes, nodeIdx    , xCoord    , yCoord + 1, zCoord)
            else
               call registerNode(Nodes, nodeIdx - 3, realXGrid(xCoord)    , realYGrid(yCoord)    , realZGrid(zCoord))
               call registerNode(Nodes, nodeIdx - 2, realXGrid(xCoord + 1), realYGrid(yCoord)    , realZGrid(zCoord))
               call registerNode(Nodes, nodeIdx - 1, realXGrid(xCoord + 1), realYGrid(yCoord + 1), realZGrid(zCoord))
               call registerNode(Nodes, nodeIdx    , realXGrid(xCoord)    , realYGrid(yCoord + 1), realZGrid(zCoord))
            endif
            quadIdx = quadIdx + 1
            call registerQuad(Quads, quadIdx, nodeIdx - 3, nodeIdx - 2, nodeIdx - 1, nodeIdx)
         end select
      end do
   end subroutine

   subroutine get_checker_and_component(request, checker, component)
      integer(kind=SINGLE), intent(in)              :: request
      procedure(logical_func), pointer, intent(out) :: checker
      integer(kind=SINGLE), intent(out)             :: component

      select case (request)
      case (iCur); checker => volumicCurrentRequest; component = iCur
      case (iMEC); checker => volumicElectricRequest; component = iMEC
      case (iMHC); checker => volumicMagneticRequest; component = iMHC
      case (iCurx); checker => componentCurrentRequest; component = iEx
      case (iExC); checker => componentFieldRequest; component = iEx
      case (iHxC); checker => componentFieldRequest; component = iHx
      case (iCurY); checker => componentCurrentRequest; component = iEy
      case (iEyC); checker => componentFieldRequest; component = iEy
      case (iHyC); checker => componentFieldRequest; component = iHy
      case (iCurZ); checker => componentCurrentRequest; component = iEz
      case (iEzC); checker => componentFieldRequest; component = iEz
      case (iHzC); checker => componentFieldRequest; component = iHz
      end select
   end subroutine get_checker_and_component

   !--------------------------------------------------------------------------
   ! Logic Functions
   !--------------------------------------------------------------------------

   logical function isValidPointForCurrent(request, i, j, k, problemInfo)
      integer(kind=SINGLE), intent(in) :: request, i, j, k
      type(problem_info_t), intent(in) :: problemInfo
      select case (request)
      case (iCur)
         isValidPointForCurrent = volumicCurrentRequest(request, i, j, k, problemInfo)
      case (iEx, iEy, iEz)
         isValidPointForCurrent = componentCurrentRequest(request, i, j, k, problemInfo)
      case default
         isValidPointForCurrent = .false.
      end select
   end function isValidPointForCurrent

   logical function isValidPointForField(request, i, j, k, problemInfo)
      integer(kind=SINGLE), intent(in) :: request, i, j, k
      type(problem_info_t), intent(in) :: problemInfo
      select case (request)
      case (iMEC)
         isValidPointForField = volumicElectricRequest(request, i, j, k, problemInfo)
      case (iMHC)
         isValidPointForField = volumicMagneticRequest(request, i, j, k, problemInfo)
      case (iEx, iEy, iEz, iHx, iHy, iHz)
         isValidPointForField = componentFieldRequest(request, i, j, k, problemInfo)
      case default
         isValidPointForField = .false.
      end select
   end function isValidPointForField

   logical function volumicCurrentRequest(request, i, j, k, problemInfo)
      integer(kind=SINGLE), intent(in) :: request, i, j, k
      type(problem_info_t), intent(in) :: problemInfo
      volumicCurrentRequest = componentCurrentRequest(iEx, i, j, k, problemInfo) .or. &
                              componentCurrentRequest(iEy, i, j, k, problemInfo) .or. &
                              componentCurrentRequest(iEz, i, j, k, problemInfo)
   end function volumicCurrentRequest

   logical function volumicElectricRequest(request, i, j, k, problemInfo)
      integer(kind=SINGLE), intent(in) :: request, i, j, k
      type(problem_info_t), intent(in) :: problemInfo
      volumicElectricRequest = componentFieldRequest(iEx, i, j, k, problemInfo) .or. &
                               componentFieldRequest(iEy, i, j, k, problemInfo) .or. &
                               componentFieldRequest(iEz, i, j, k, problemInfo)
   end function volumicElectricRequest

   logical function volumicMagneticRequest(request, i, j, k, problemInfo)
      integer(kind=SINGLE), intent(in) :: request, i, j, k
      type(problem_info_t), intent(in) :: problemInfo
      volumicMagneticRequest = componentFieldRequest(iHx, i, j, k, problemInfo) .or. &
                               componentFieldRequest(iHy, i, j, k, problemInfo) .or. &
                               componentFieldRequest(iHz, i, j, k, problemInfo)
   end function volumicMagneticRequest

   logical function componentCurrentRequest(fieldDir, i, j, k, problemInfo)
      integer(kind=SINGLE), intent(in) :: fieldDir, i, j, k
      type(problem_info_t), intent(in) :: problemInfo
      componentCurrentRequest = isWithinBounds(fieldDir, i, j, k, problemInfo)
      if (componentCurrentRequest) then
         componentCurrentRequest = isPEC(fieldDir, i, j, k, problemInfo) .or. isThinWire(fieldDir, i, j, k, problemInfo)
      end if
   end function componentCurrentRequest

   logical function componentFieldRequest(fieldDir, i, j, k, problemInfo)
      integer(kind=SINGLE), intent(in) :: fieldDir, i, j, k
      type(problem_info_t), intent(in) :: problemInfo
      componentFieldRequest = isWithinBounds(fieldDir, i, j, k, problemInfo)
   end function componentFieldRequest

end module mod_volumicProbeUtils
