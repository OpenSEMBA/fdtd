module outputDecomposition_m
   use iso_fortran_env, only: int64
   use FDETYPES_m, only: SINGLE, limit_t, iEx, iEy, iEz, iHx, iHy, iHz
   use outputTypes_m, only: cell_coordinate_t
   implicit none
   private

   integer, parameter, public :: OUTPUT_PARTITION_SUCCESS = 0
   integer, parameter, public :: OUTPUT_PARTITION_INVALID_ARGUMENT = 1
   integer, parameter, public :: OUTPUT_PARTITION_INVALID_COMPONENT = 2

   type, public :: output_partition_t
      ! Shapes and offsets use internal (x,y,z) order; file offsets are zero-based.
      type(cell_coordinate_t) :: global_lower = cell_coordinate_t(0_SINGLE, 0_SINGLE, 0_SINGLE)
      type(cell_coordinate_t) :: global_upper = cell_coordinate_t(0_SINGLE, 0_SINGLE, 0_SINGLE)
      type(cell_coordinate_t) :: local_lower = cell_coordinate_t(0_SINGLE, 0_SINGLE, 0_SINGLE)
      type(cell_coordinate_t) :: local_upper = cell_coordinate_t(0_SINGLE, 0_SINGLE, 0_SINGLE)
      integer(int64) :: global_shape(3) = 0_int64
      integer(int64) :: file_offset(3) = 0_int64
      integer(int64) :: local_shape(3) = 0_int64
      logical :: has_data = .false.
   end type output_partition_t

   public :: build_output_partition
   public :: point_is_in_partition

contains

   pure subroutine build_output_partition(request_lower, request_upper, global_bounds, local_sweep, &
                                          field_component, rank, rank_count, partition, status)
      type(cell_coordinate_t), intent(in) :: request_lower, request_upper
      type(limit_t), intent(in) :: global_bounds, local_sweep
      integer, intent(in) :: field_component, rank, rank_count
      type(output_partition_t), intent(out) :: partition
      integer, intent(out) :: status

      integer(kind=SINGLE) :: owned_upper_z

      partition = output_partition_t()
      status = OUTPUT_PARTITION_SUCCESS

      if (.not. is_supported_component(field_component)) then
         status = OUTPUT_PARTITION_INVALID_COMPONENT
         return
      end if

      if (rank_count < 1 .or. rank < 0 .or. rank >= rank_count .or. &
          .not. is_valid_box(request_lower, request_upper) .or. &
          .not. is_valid_limit(global_bounds) .or. .not. is_valid_limit(local_sweep)) then
         status = OUTPUT_PARTITION_INVALID_ARGUMENT
         return
      end if

      partition%global_lower = cell_coordinate_t( &
                               max(request_lower%x, global_bounds%XI), &
                               max(request_lower%y, global_bounds%YI), &
                               max(request_lower%z, global_bounds%ZI))
      partition%global_upper = cell_coordinate_t( &
                               min(request_upper%x, global_bounds%XE), &
                               min(request_upper%y, global_bounds%YE), &
                               min(request_upper%z, global_bounds%ZE))

      if (.not. is_valid_box(partition%global_lower, partition%global_upper)) return

      partition%global_shape = shape_of(partition%global_lower, partition%global_upper)

      owned_upper_z = local_sweep%ZE
      ! Ex, Ey, and Hz sweeps overlap at MPI interfaces; the higher rank owns the plane.
      if (rank < rank_count - 1 .and. has_shared_upper_plane(field_component)) then
         owned_upper_z = owned_upper_z - 1_SINGLE
      end if

      partition%local_lower = cell_coordinate_t( &
                              max(partition%global_lower%x, local_sweep%XI), &
                              max(partition%global_lower%y, local_sweep%YI), &
                              max(partition%global_lower%z, local_sweep%ZI))
      partition%local_upper = cell_coordinate_t( &
                              min(partition%global_upper%x, local_sweep%XE), &
                              min(partition%global_upper%y, local_sweep%YE), &
                              min(partition%global_upper%z, owned_upper_z))

      if (.not. is_valid_box(partition%local_lower, partition%local_upper)) return

      partition%local_shape = shape_of(partition%local_lower, partition%local_upper)
      partition%file_offset = [ &
                              int(partition%local_lower%x, int64) - int(partition%global_lower%x, int64), &
                              int(partition%local_lower%y, int64) - int(partition%global_lower%y, int64), &
                              int(partition%local_lower%z, int64) - int(partition%global_lower%z, int64)]
      partition%has_data = .true.
   end subroutine build_output_partition

   pure logical function point_is_in_partition(point, partition)
      type(cell_coordinate_t), intent(in) :: point
      type(output_partition_t), intent(in) :: partition

      point_is_in_partition = partition%has_data .and. &
                              point%x >= partition%local_lower%x .and. point%x <= partition%local_upper%x .and. &
                              point%y >= partition%local_lower%y .and. point%y <= partition%local_upper%y .and. &
                              point%z >= partition%local_lower%z .and. point%z <= partition%local_upper%z
   end function point_is_in_partition

   pure logical function is_supported_component(field_component)
      integer, intent(in) :: field_component

      is_supported_component = any(field_component == [iEx, iEy, iEz, iHx, iHy, iHz])
   end function is_supported_component

   pure logical function has_shared_upper_plane(field_component)
      integer, intent(in) :: field_component

      has_shared_upper_plane = any(field_component == [iEx, iEy, iHz])
   end function has_shared_upper_plane

   pure logical function is_valid_limit(bounds)
      type(limit_t), intent(in) :: bounds

      is_valid_limit = bounds%XI <= bounds%XE .and. &
                       bounds%YI <= bounds%YE .and. &
                       bounds%ZI <= bounds%ZE
   end function is_valid_limit

   pure logical function is_valid_box(lower, upper)
      type(cell_coordinate_t), intent(in) :: lower, upper

      is_valid_box = lower%x <= upper%x .and. &
                     lower%y <= upper%y .and. &
                     lower%z <= upper%z
   end function is_valid_box

   pure function shape_of(lower, upper) result(shape)
      type(cell_coordinate_t), intent(in) :: lower, upper
      integer(int64) :: shape(3)

      shape = [ &
              int(upper%x, int64) - int(lower%x, int64) + 1_int64, &
              int(upper%y, int64) - int(lower%y, int64) + 1_int64, &
              int(upper%z, int64) - int(lower%z, int64) + 1_int64]
   end function shape_of

end module outputDecomposition_m
