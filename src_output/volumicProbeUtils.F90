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

   abstract interface
      logical function logical_func(component, i, j, k, problemInfo)
         import :: problem_info_t, SINGLE
         type(problem_info_t), intent(in) :: problemInfo
         integer(kind=SINGLE), intent(in) :: component, i, j, k
      end function logical_func
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
      do i = lowerBound%x, upperBound%x
      do j = lowerBound%y, upperBound%y
      do k = lowerBound%z, upperBound%z
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
      do i = lowerBound%x, upperBound%x
      do j = lowerBound%y, upperBound%y
      do k = lowerBound%z, upperBound%z
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
