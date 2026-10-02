! Builder of slanted (arbitrary-orientation) MTLN wire divisions.
!
! A polyline that is not axis-aligned on the mesh grid is discretized into a
! sequence of divisions, each one fully contained in a single FDTD cell. Every
! division carries the list of E-field edges it couples to, following the
! arbitrary-orientation thin wire formalism: the wire tangential field is a
! bilinear interpolation of the three staggered E components, and the wire
! current is injected back into the same edges with the adjoint weights.
!
! Coordinates coming from the JSON input are cell-relative (cell index plus a
! fraction). Positions handled inside this module are physical.
module slanted_wire_segments_m

   use FDETYPES_m, only: rkind
   use NFDETypes_m, only: Desplazamiento_t
   use mtln_types_m, only: segment_t, segment_coupling_t, MAX_SEGMENT_COUPLINGS, &
                           DIRECTION_X_POS, DIRECTION_Y_POS, DIRECTION_Z_POS
   use Report_m, only: WarnErrReport

   implicit none
   private

   public :: buildSlantedSegments
   public :: findIndexInSlantedSegments

   integer, parameter :: X_AXIS = 1
   integer, parameter :: Y_AXIS = 2
   integer, parameter :: Z_AXIS = 3

   ! Tolerance along the segment parameter when locating face crossings.
   real(kind=rkind), parameter :: WALK_TOLERANCE = 1.0e-10_rkind
   ! Couplings whose weight is below this value are dropped.
   real(kind=rkind), parameter :: MIN_COUPLING_WEIGHT = 1.0e-12_rkind
   ! Tolerance used to decide whether a cell-relative coordinate is integer.
   real(kind=rkind), parameter :: INTEGER_TOLERANCE = 1.0e-9_rkind

   ! Uniform view of the graded grid: cell sizes and grid line positions.
   type :: cell_grid_t
      integer :: n(3) = 0
      real(kind=rkind), allocatable :: step(:, :)
      real(kind=rkind), allocatable :: line(:, :)
   end type

contains

   ! Builds the slanted divisions of a whole polyline. Points are given as
   ! cell-relative coordinates, column-wise: cell_positions(:, i).
   subroutine buildSlantedSegments(cell_positions, grid, segments)
      real(kind=rkind), dimension(:, :), intent(in) :: cell_positions
      type(Desplazamiento_t), intent(in) :: grid
      type(segment_t), dimension(:), allocatable, intent(out) :: segments

      type(cell_grid_t) :: g
      integer :: i, p, ndiv, idx, ndiv_piece
      integer :: cs(3), ce(3), maxtmp
      integer, allocatable :: cells(:,:)
      real(kind=rkind) :: p0(3), p1(3), tdir(3), length
      real(kind=rkind), allocatable :: points(:,:)
      real(kind=rkind), allocatable :: directions(:, :)
      logical :: ok

      if (size(cell_positions, 1) /= 3 .or. size(cell_positions, 2) < 2) then
         call WarnErrReport('buildSlantedSegments: polyline must have at least '// &
                            'two three-dimensional points.', .true.)
         allocate (segments(0))
         return
      end if

      g = buildCellGrid(grid)

      maxtmp = 1
      do i = 1, size(cell_positions, 2) - 1
         call validatePoint(cell_positions(:, i), g)
         call validatePoint(cell_positions(:, i + 1), g)
         p0 = positionOfCellCoordinate(g, cell_positions(:, i))
         p1 = positionOfCellCoordinate(g, cell_positions(:, i + 1))
         cs = startCellOf(g, cell_positions(:, i), p1 - p0)
         ce = endCellOf(g, cell_positions(:, i + 1), p1 - p0)
         maxtmp = max(maxtmp, 1 + sum(abs(ce - cs)))
      end do
      allocate (cells(3, maxtmp), points(3, maxtmp + 1))
      allocate (directions(3, size(cell_positions, 2) - 1))

      ! Count divisions.
      ndiv = 0
      do i = 1, size(cell_positions, 2) - 1
         p0 = positionOfCellCoordinate(g, cell_positions(:, i))
         p1 = positionOfCellCoordinate(g, cell_positions(:, i + 1))
         length = norm2(p1 - p0)
         if (length <= 0.0_rkind) then
            call WarnErrReport('buildSlantedSegments: zero length wire segment.', .true.)
            allocate (segments(0))
            return
         end if
         tdir = (p1 - p0)/length
         directions(:, i) = tdir
         cs = startCellOf(g, cell_positions(:, i), tdir)
         ce = endCellOf(g, cell_positions(:, i + 1), tdir)
         call walkSegment(g, p0, p1, cs, ce, cells, points, ndiv_piece, ok)
         if (.not. ok) then
            call WarnErrReport('buildSlantedSegments: could not walk the wire '// &
                               'through the mesh cells.', .true.)
            allocate (segments(0))
            return
         end if
         ndiv = ndiv + ndiv_piece
      end do

      allocate (segments(ndiv))

      ! Fill divisions.
      idx = 0
      do i = 1, size(cell_positions, 2) - 1
         p0 = positionOfCellCoordinate(g, cell_positions(:, i))
         p1 = positionOfCellCoordinate(g, cell_positions(:, i + 1))
         tdir = directions(:, i)
         cs = startCellOf(g, cell_positions(:, i), tdir)
         ce = endCellOf(g, cell_positions(:, i + 1), tdir)
         call walkSegment(g, p0, p1, cs, ce, cells, points, ndiv, ok)
         do p = 1, ndiv
            idx = idx + 1
            call fillSegment(g, tdir, cells(:, p), points(:, p), points(:, p + 1), segments(idx))
         end do
      end do

   end subroutine buildSlantedSegments

   ! Index of the closest division boundary to a cell-relative coordinate.
   ! Division i starts at point i, and the wire end is point n+1, so the
   ! result follows the same node numbering used by MTLN probes and generators.
   function findIndexInSlantedSegments(coord, grid, segments) result(res)
      real(kind=rkind), dimension(3), intent(in) :: coord
      type(Desplazamiento_t), intent(in) :: grid
      type(segment_t), dimension(:), intent(in) :: segments
      integer :: res

      type(cell_grid_t) :: g
      real(kind=rkind) :: p(3), dist, best
      integer :: i

      g = buildCellGrid(grid)
      p = positionOfCellCoordinate(g, coord)
      res = 1
      best = huge(best)
      do i = 1, size(segments)
         dist = norm2(p - segments(i)%position_begin)
         if (dist < best) then
            best = dist
            res = i
         end if
      end do
      if (size(segments) > 0) then
         dist = norm2(p - segments(size(segments))%position_end)
         if (dist < best) res = size(segments) + 1
      end if
   end function findIndexInSlantedSegments

   subroutine fillSegment(g, tdir, cell, p0, p1, seg)
      type(cell_grid_t), intent(in) :: g
      real(kind=rkind), dimension(3), intent(in) :: tdir, p0, p1
      integer, dimension(3), intent(in) :: cell
      type(segment_t), intent(inout) :: seg

      real(kind=rkind) :: chord, transverse(2)
      integer :: axis

      seg%is_slanted = .true.
      seg%x = cell(X_AXIS)
      seg%y = cell(Y_AXIS)
      seg%z = cell(Z_AXIS)
      seg%direction_cosines = tdir
      seg%position_begin = p0
      seg%position_end = p1
      axis = dominantAxis(tdir, g, cell)
      seg%orientation = axis*sign(1, nint(tdir(axis)))
      transverse = transverseSteps(g, cell, axis)
      seg%d1 = transverse(1)
      seg%d2 = transverse(2)
      chord = norm2(p1 - p0)
      call buildCouplings(g, cell, p0, p1, tdir, chord, seg%couplings, seg%n_couplings)
   end subroutine fillSegment

   ! Bilinear interpolation weights of the four E edges surrounding the wire
   ! inside one cell. The component with the largest projection over its cell
   ! step uses the average of the entry and exit positions (a trapezoidal
   ! approximation of the line integral), as in the reference implementation.
   subroutine buildCouplings(g, cell, p0, p1, tdir, chord, couplings, n_couplings)
      type(cell_grid_t), intent(in) :: g
      integer, dimension(3), intent(in) :: cell
      real(kind=rkind), dimension(3), intent(in) :: p0, p1, tdir
      real(kind=rkind), intent(in) :: chord
      type(segment_coupling_t), dimension(MAX_SEGMENT_COUPLINGS), intent(inout) :: couplings
      integer, intent(out) :: n_couplings

      real(kind=rkind) :: u0(3), u1(3), um(3), pond(3, 4), w, delta1, delta2
      integer :: d, d1, d2, offset, o1, o2, n, edge(3), axis

      axis = dominantAxis(tdir, g, cell)
      do d = 1, 3
         u0(d) = (p0(d) - g%line(d, cell(d)))/g%step(d, cell(d))
         u1(d) = (p1(d) - g%line(d, cell(d)))/g%step(d, cell(d))
      end do
      um = 0.5_rkind*(u0 + u1)
      do d = 1, 3
         d1 = mod(d, 3) + 1
         d2 = mod(d + 1, 3) + 1
         if (d == axis) then
            ! Average of the entry and exit interpolation weights, a
            ! trapezoidal approximation of the line integral along the wire.
            pond(d, 1) = 0.5_rkind*((1.0_rkind - u0(d1))*(1.0_rkind - u0(d2)) + &
                                    (1.0_rkind - u1(d1))*(1.0_rkind - u1(d2)))
            pond(d, 2) = 0.5_rkind*(u0(d1)*(1.0_rkind - u0(d2)) + &
                                    u1(d1)*(1.0_rkind - u1(d2)))
            pond(d, 3) = 0.5_rkind*((1.0_rkind - u0(d1))*u0(d2) + &
                                    (1.0_rkind - u1(d1))*u1(d2))
            pond(d, 4) = 0.5_rkind*(u0(d1)*u0(d2) + u1(d1)*u1(d2))
         else
            delta1 = um(d1)
            delta2 = um(d2)
            pond(d, 1) = (1.0_rkind - delta1)*(1.0_rkind - delta2)
            pond(d, 2) = delta1*(1.0_rkind - delta2)
            pond(d, 3) = (1.0_rkind - delta1)*delta2
            pond(d, 4) = delta1*delta2
         end if
      end do

      n = 0
      do d = 1, 3
         d1 = mod(d, 3) + 1
         d2 = mod(d + 1, 3) + 1
         do offset = 1, 4
            w = pond(d, offset)*tdir(d)
            if (abs(w) < MIN_COUPLING_WEIGHT) cycle
            o1 = mod(offset - 1, 2)
            o2 = (offset - 1)/2
            edge = cell
            edge(d1) = cell(d1) + o1
            edge(d2) = cell(d2) + o2
            if (any(edge < 0) .or. any(edge > g%n)) then
               call WarnErrReport('buildSlantedSegments: wire coupling falls '// &
                                  'outside the mesh.', .true.)
               cycle
            end if
            n = n + 1
            couplings(n)%i = edge(X_AXIS)
            couplings(n)%j = edge(Y_AXIS)
            couplings(n)%k = edge(Z_AXIS)
            couplings(n)%component = d
            couplings(n)%weight = w
            couplings(n)%chord = chord
         end do
      end do
      n_couplings = n
   end subroutine buildCouplings

   ! Walks a straight segment in physical space, splitting it at cell faces.
   ! cells(:, p) is the cell of division p, and points(:, p:p+1) its limits.
   subroutine walkSegment(g, p0, p1, cs, ce, cells, points, ndiv, ok)
      type(cell_grid_t), intent(in) :: g
      real(kind=rkind), dimension(3), intent(in) :: p0, p1
      integer, dimension(3), intent(in) :: cs, ce
      integer, dimension(:, :), intent(inout) :: cells
      real(kind=rkind), dimension(:, :), intent(inout) :: points
      integer, intent(out) :: ndiv
      logical, intent(out) :: ok

      integer :: c(3), d
      real(kind=rkind) :: dir(3), p(3), pnext(3), t, tface, face
      real(kind=rkind) :: td(3)
      logical :: tied(3)
      integer :: maxdiv

      maxdiv = size(cells, 2)
      dir = p1 - p0
      c = cs
      p = p0
      t = 0.0_rkind
      ndiv = 0
      ok = .true.
      do
         if (all(c == ce)) then
            call appendDivision(cells, points, ndiv, maxdiv, c, p, p1)
            exit
         end if

         td = -1.0_rkind
         do d = 1, 3
            if (dir(d) == 0.0_rkind) cycle
            if (dir(d) > 0.0_rkind) then
               face = g%line(d, c(d) + 1)
            else
               face = g%line(d, c(d))
            end if
            td(d) = (face - p(d))/dir(d)
            if (td(d) < 0.0_rkind .and. td(d) > -WALK_TOLERANCE) td(d) = 0.0_rkind
         end do
         tface = huge(tface)
         do d = 1, 3
            if (dir(d) /= 0.0_rkind) tface = min(tface, td(d))
         end do
         if (t + tface >= 1.0_rkind - WALK_TOLERANCE) then
            call appendDivision(cells, points, ndiv, maxdiv, c, p, p1)
            exit
         end if
         if (tface < 0.0_rkind) then
            ok = .false.
            return
         end if

         pnext = p + tface*dir
         tied = .false.
         do d = 1, 3
            if (dir(d) == 0.0_rkind) cycle
            if (abs(td(d) - tface) <= WALK_TOLERANCE*max(1.0_rkind, abs(tface))) then
               tied(d) = .true.
               if (dir(d) > 0.0_rkind) then
                  pnext(d) = g%line(d, c(d) + 1)
               else
                  pnext(d) = g%line(d, c(d))
               end if
            end if
         end do
         call appendDivision(cells, points, ndiv, maxdiv, c, p, pnext)
         if (ndiv < 0) then
            ok = .false.
            return
         end if
         do d = 1, 3
            if (tied(d)) c(d) = c(d) + nint(sign(1.0_rkind, dir(d)))
         end do
         p = pnext
         t = t + tface
      end do
      if (ndiv < 0) ok = .false.
   end subroutine walkSegment

   subroutine appendDivision(cells, points, ndiv, maxdiv, cell, p0, p1)
      integer, dimension(:, :), intent(inout) :: cells
      real(kind=rkind), dimension(:, :), intent(inout) :: points
      integer, intent(inout) :: ndiv
      integer, intent(in) :: maxdiv
      integer, dimension(3), intent(in) :: cell
      real(kind=rkind), dimension(3), intent(in) :: p0, p1

      if (ndiv >= maxdiv) then
         ndiv = -1
         return
      end if
      ndiv = ndiv + 1
      cells(:, ndiv) = cell
      points(:, ndiv) = p0
      points(:, ndiv + 1) = p1
   end subroutine appendDivision

   function buildCellGrid(grid) result(res)
      type(Desplazamiento_t), intent(in) :: grid
      type(cell_grid_t) :: res

      integer :: d, c, n, maxn

      ! The number of cells must not be read from nx/ny/nz: readGrid sets
      ! them to 1 for uniform grids, where a single step value is stored.
      ! The enlarged step arrays carry one value per cell, and mx2/my2/mz2
      ! keep the cell count of a raw (not enlarged) grid.
      res%n = [gridCellCount(grid%desX, grid%mx2), &
               gridCellCount(grid%desY, grid%my2), &
               gridCellCount(grid%desZ, grid%mz2)]
      if (any(res%n <= 0)) then
         call WarnErrReport('buildSlantedSegments: the mesh grid is not initialized.', .true.)
         res%n = max(res%n, 1)
      end if
      maxn = maxval(res%n)
      allocate (res%step(3, 0:maxn - 1), source=0.0_rkind)
      allocate (res%line(3, 0:maxn), source=0.0_rkind)
      do d = 1, 3
         n = res%n(d)
         do c = 0, n - 1
            res%step(d, c) = getGridStep(grid, d, c)
         end do
         do c = 1, n
            res%line(d, c) = res%line(d, c - 1) + res%step(d, c - 1)
         end do
      end do
   end function buildCellGrid

   function gridCellCount(des, numberOfCells) result(res)
      real(kind=RKIND), dimension(:), pointer, intent(in) :: des
      integer, intent(in) :: numberOfCells
      integer :: res

      res = max(size(des), numberOfCells)
   end function gridCellCount

   function getGridStep(grid, d, c) result(res)
      type(Desplazamiento_t), intent(in) :: grid
      integer, intent(in) :: d, c
      real(kind=rkind) :: res

      select case (d)
      case (X_AXIS)
         if (size(grid%desX) == 1) then
            res = grid%desX(1)
         else
            res = grid%desX(lbound(grid%desX, 1) + c)
         end if
      case (Y_AXIS)
         if (size(grid%desY) == 1) then
            res = grid%desY(1)
         else
            res = grid%desY(lbound(grid%desY, 1) + c)
         end if
      case (Z_AXIS)
         if (size(grid%desZ) == 1) then
            res = grid%desZ(1)
         else
            res = grid%desZ(lbound(grid%desZ, 1) + c)
         end if
      end select
   end function getGridStep

   function positionOfCellCoordinate(g, u) result(p)
      type(cell_grid_t), intent(in) :: g
      real(kind=rkind), dimension(3), intent(in) :: u
      real(kind=rkind), dimension(3) :: p

      integer :: d, c

      do d = 1, 3
         c = min(max(int(floor(u(d))), 0), g%n(d) - 1)
         p(d) = g%line(d, c) + (u(d) - real(c, rkind))*g%step(d, c)
      end do
   end function positionOfCellCoordinate

   subroutine validatePoint(u, g)
      real(kind=rkind), dimension(3), intent(in) :: u
      type(cell_grid_t), intent(in) :: g

      character(len=256) :: msg
      if (any(u < -INTEGER_TOLERANCE) .or. &
          any(u > real(g%n, rkind) + INTEGER_TOLERANCE)) then
         write(msg, '(a,3(1x,f10.3),a,3(1x,i5))') 'buildSlantedSegments: wire coordinate ', &
            u, ' outside the mesh ', g%n
         call WarnErrReport(msg, .true.)
      end if
   end subroutine validatePoint

   function startCellOf(g, u, dir) result(c)
      type(cell_grid_t), intent(in) :: g
      real(kind=rkind), dimension(3), intent(in) :: u, dir
      integer :: c(3)

      integer :: d

      do d = 1, 3
         c(d) = int(floor(u(d)))
         if (isInteger(u(d)) .and. dir(d) < 0.0_rkind) c(d) = c(d) - 1
         c(d) = min(max(c(d), 0), g%n(d) - 1)
      end do
   end function startCellOf

   function endCellOf(g, u, dir) result(c)
      type(cell_grid_t), intent(in) :: g
      real(kind=rkind), dimension(3), intent(in) :: u, dir
      integer :: c(3)

      integer :: d

      do d = 1, 3
         c(d) = int(floor(u(d)))
         if (isInteger(u(d)) .and. dir(d) > 0.0_rkind) c(d) = c(d) - 1
         c(d) = min(max(c(d), 0), g%n(d) - 1)
      end do
   end function endCellOf

   logical function isInteger(x)
      real(kind=rkind), intent(in) :: x
      isInteger = abs(x - real(nint(x), rkind)) <= INTEGER_TOLERANCE
   end function isInteger

   function dominantAxis(tdir, g, cell) result(res)
      real(kind=rkind), dimension(3), intent(in) :: tdir
      type(cell_grid_t), intent(in) :: g
      integer, dimension(3), intent(in) :: cell
      integer :: res

      real(kind=rkind) :: ratio(3)

      ratio(X_AXIS) = abs(tdir(X_AXIS))/g%step(X_AXIS, cell(X_AXIS))
      ratio(Y_AXIS) = abs(tdir(Y_AXIS))/g%step(Y_AXIS, cell(Y_AXIS))
      ratio(Z_AXIS) = abs(tdir(Z_AXIS))/g%step(Z_AXIS, cell(Z_AXIS))
      res = maxloc(ratio, dim=1)
   end function dominantAxis

   function transverseSteps(g, cell, axis) result(res)
      type(cell_grid_t), intent(in) :: g
      integer, dimension(3), intent(in) :: cell
      integer, intent(in) :: axis
      real(kind=rkind) :: res(2)

      select case (axis)
      case (X_AXIS)
         res = [g%step(Y_AXIS, cell(Y_AXIS)), g%step(Z_AXIS, cell(Z_AXIS))]
      case (Y_AXIS)
         res = [g%step(Z_AXIS, cell(Z_AXIS)), g%step(X_AXIS, cell(X_AXIS))]
      case (Z_AXIS)
         res = [g%step(X_AXIS, cell(X_AXIS)), g%step(Y_AXIS, cell(Y_AXIS))]
      end select
   end function transverseSteps

end module slanted_wire_segments_m
