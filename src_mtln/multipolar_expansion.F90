module multipolar_expansion_mod

   use mtln_types_mod
   use Report
   use FDETypes, only: pi, EPSILON_VACUUM, MU_VACUUM

   type, private :: integration_grid_t
      real, dimension(:), allocatable :: x, y
   end type

contains
   function multipolarExpansion2D(position, ab, expansionCenter) result(res)

      ! 2D multipolar expansion from:
      ! TSOGTGEREL GANTUMUR, MULTIPOLE EXPANSIONS IN THE PLANE.
      ! 2016-04-16 lecture notes.

      implicit none

      real, intent(in) :: position(2)
      real, intent(in) :: expansionCenter(2)
      type(multipolar_coefficient_t), dimension(:), intent(in) :: ab
      real :: res

      real :: rVec(2), r, phi
      integer :: n

      rVec = position - expansionCenter
      r = sqrt(rVec(1)**2 + rVec(2)**2)
      phi = atan2(rVec(2), rVec(1))

      res = 0.0
      do n = 1, size(ab)
         if (n == 1) then
            res = res - ab(n)%a * log(r)
            ! b0 should always be zero
         else
            res = res + (ab(n)%a * cos((n-1)*phi) + ab(n)%b * sin((n-1)*phi)) / r**(n-1)
         end if
      end do

      res = res / (2.0 * pi)
   end function

   function buildIntegrationGridForBox(integrationBox, innerRegionBox) result(res)
      ! Returns an integration_grid_t with sorted, unique, equispaced points for x and y
      type(box_2d_t), intent(in) :: integrationBox, innerRegionBox
      type(integration_grid_t) :: res

      integer, parameter :: GRID_INTEGRATION_SAMPLING_POINTS = 100
      real :: minval, maxval, step
      integer :: x, k
      real, allocatable :: controlPoints(:)
      real, allocatable :: allPoints(:)
      integer :: i, j

      ! Preconditions
      if (any(integrationBox%min >= innerRegionBox%min) .or. &
          any(integrationBox%max <= innerRegionBox%max)) then       
         call WarnErrReport( &
            "Error in mutipolar expansion innerRegion must be fully contained within the integration Box", .true.)
         return
      end if

      block
         real, dimension(2) :: innerRegionSize, integrationBoxSize         
         innerRegionSize = innerRegionBox%max - innerRegionBox%min
         integrationBoxSize = integrationBox%max - integrationBox%min
         if (any(integrationBoxSize < innerRegionSize * 1.25)) then
            call WarnErrReport( &
               "Error in multipolar expansion: integration box is too small for the inner region", .true.)
            return
         end if
      end block

      ! Builds grid
      allocate(allPoints( GRID_INTEGRATION_SAMPLING_POINTS * 3 + 1))
      do x = 1, 2 
         ! control points are ordered from min to max.
         controlPoints = [&
            integrationBox%min(x), &
            innerRegionBox%min(x), &
            innerRegionBox%max(x), &
            integrationBox%max(x)]
         
         do i = 2, size(controlPoints)
            minval = controlPoints(i-1)
            maxval = controlPoints(i)
            step = (maxval - minval) / GRID_INTEGRATION_SAMPLING_POINTS
            do k = 1, GRID_INTEGRATION_SAMPLING_POINTS
               allPoints((i-2)*GRID_INTEGRATION_SAMPLING_POINTS+k) = minval + (k-1)*step
            end do
         end do
         allPoints(size(allPoints)) = controlPoints(size(controlPoints))

         if (x == 1) then
            res%x = allPoints
         else
            res%y = allPoints
         end if

      end do
   end function

   real function boxArea(box) result(res)
      type(box_2d_t), intent(in) :: box
      res = (box%max(1) - box%min(1)) * (box%max(2) - box%min(2))
   end function

   logical function isWithinBox(box, point) result(res)
      type(box_2d_t), intent(in) :: box
      real, dimension(2), intent(in) :: point
      res = all(point >= box%min) .and. all(point <= box%max)
   end function

   function getAveragePotential(potential, innerBox, outerBox) result(avVj)
      type(field_reconstruction_t), intent(in) :: potential
      type(box_2d_t), intent(in) :: innerBox, outerBox
      real :: avVj

      type(integration_grid_t) :: integrationGrid
      real :: outerV, innerV, area
      real :: xMin, xMax, yMin, yMax
      real :: midPoint(2), center(2)
      integer :: m, n, nx, ny

      integrationGrid = buildIntegrationGridForBox(outerBox, innerBox)

      outerV = 0.0
      do m = 2, size(integrationGrid%x)
         do n = 2, size(integrationGrid%y)
            xMin = integrationGrid%x(m - 1)
            xMax = integrationGrid%x(m)
            yMin = integrationGrid%y(n - 1)
            yMax = integrationGrid%y(n)

            midPoint = [0.5 * (xMin + xMax), 0.5 * (yMin + yMax)]
            area = (xMax - xMin) * (yMax - yMin)

            if (isWithinBox(innerBox, midPoint)) then
               cycle ! skip if inside inner region
            end if

            center = potential%expansion_center
            outerV = outerV + area * multipolarExpansion2D(midPoint, potential%ab, center)
         end do
      end do

      innerV = potential%inner_region_average_potential * boxArea(innerBox)
      avVj = (innerV + outerV) / boxArea(outerBox)
   end function
  
   function getCellCapacitanceOnBox(multipolarExpansionParameters, cellBox) result (res)
      real, dimension(:,:), allocatable :: res

      type(multipolar_expansion_t), intent(in) :: multipolarExpansionParameters
      type(box_2d_t), intent(in) :: cellBox

      real :: Qj, avVj, ViWhenPrescribedVj

      integer :: i, j
      integer :: N

      N = size(multipolarExpansionParameters%electric)
      allocate(res(N,N))
      do i = 1, N
         do j = 1, N
            Qj = multipolarExpansionParameters%electric(j)%ab(1)%a
            avVj = getAveragePotential(&
               multipolarExpansionParameters%electric(j), &
               multipolarExpansionParameters%inner_region, &
               cellBox)
            ViWhenPrescribedVj = multipolarExpansionParameters%electric(j)%conductor_potentials(i)
            avVj = -avVj + ViWhenPrescribedVj

            res(i,j) = Qj / avVj * EPSILON_VACUUM
         end do
      end do

   end function
  
   function getCellInductanceOnBox(multipolarExpansionParameters, cellBox) result (res)
      real, dimension(:,:), allocatable :: res

      type(multipolar_expansion_t), intent(in) :: multipolarExpansionParameters
      type(box_2d_t), intent(in) :: cellBox

      real :: Ij, avAj, ViWhenPrescribedVj

      integer :: i, j
      integer :: N

      N = size(multipolarExpansionParameters%magnetic)
      allocate(res(N,N))
      do i = 1, N
         do j = 1, N
            Ij = multipolarExpansionParameters%magnetic(j)%ab(1)%a
            avAj = getAveragePotential(&
               multipolarExpansionParameters%magnetic(j), &
               multipolarExpansionParameters%inner_region, &
               cellBox)
            ViWhenPrescribedVj = multipolarExpansionParameters%magnetic(j)%conductor_potentials(i)
            avAj = -avAj + ViWhenPrescribedVj

            res(i,j) = avAj / Ij * MU_VACUUM
         end do
      end do

   end function

end module multipolar_expansion_mod
