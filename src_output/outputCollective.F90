module outputCollective_m
   use outputDecomposition_m, only: output_partition_t
   implicit none
   private

   integer, parameter, public :: OUTPUT_COLLECTIVE_SUCCESS = 0
   integer, parameter, public :: OUTPUT_COLLECTIVE_INVALID_CONTEXT = 1
   integer, parameter, public :: OUTPUT_COLLECTIVE_INVALID_PARTICIPANTS = 2
   integer, parameter, public :: OUTPUT_COLLECTIVE_INVALID_OWNER = 3
   integer, parameter, public :: OUTPUT_COLLECTIVE_INVALID_PARTITION = 4

   integer, parameter, public :: OUTPUT_PUBLICATION_COLLECTIVE = 1
   integer, parameter, public :: OUTPUT_PUBLICATION_ROOT_AGGREGATION = 2

   type, public :: output_collective_t
      integer :: rank = 0
      integer :: rank_count = 1
      integer :: root_rank = 0
      logical :: collective_publication_available = .false.
   end type output_collective_t

   public :: init_output_collective
   public :: select_output_participants
   public :: validate_output_ownership
   public :: select_output_publication_mode
   public :: prepare_output_partition_publication

contains

   pure subroutine init_output_collective(collective, rank, rank_count, root_rank, &
                                          collective_publication_available, status)
      type(output_collective_t), intent(out) :: collective
      integer, intent(in) :: rank, rank_count, root_rank
      logical, intent(in) :: collective_publication_available
      integer, intent(out) :: status

      collective = output_collective_t()
      if (rank_count < 1 .or. rank < 0 .or. rank >= rank_count .or. &
          root_rank < 0 .or. root_rank >= rank_count) then
         status = OUTPUT_COLLECTIVE_INVALID_CONTEXT
         return
      end if

      collective%rank = rank
      collective%rank_count = rank_count
      collective%root_rank = root_rank
      collective%collective_publication_available = collective_publication_available
      status = OUTPUT_COLLECTIVE_SUCCESS
   end subroutine init_output_collective

   subroutine select_output_participants(collective, rank_has_data, participants, owner_rank, status)
      type(output_collective_t), intent(in) :: collective
      logical, intent(in) :: rank_has_data(:)
      integer, allocatable, intent(out) :: participants(:)
      integer, intent(out) :: owner_rank, status
      integer :: rank

      if (allocated(participants)) deallocate (participants)
      owner_rank = -1
      if (.not. valid_context(collective) .or. size(rank_has_data) /= collective%rank_count) then
         status = OUTPUT_COLLECTIVE_INVALID_CONTEXT
         return
      end if
      if (.not. any(rank_has_data)) then
         status = OUTPUT_COLLECTIVE_INVALID_PARTICIPANTS
         return
      end if

      participants = pack([(rank - 1, rank=1, collective%rank_count)], rank_has_data)
      owner_rank = participants(1)
      status = OUTPUT_COLLECTIVE_SUCCESS
   end subroutine select_output_participants

   pure subroutine validate_output_ownership(collective, participants, owner_rank, status)
      type(output_collective_t), intent(in) :: collective
      integer, intent(in) :: participants(:)
      integer, intent(in) :: owner_rank
      integer, intent(out) :: status
      integer :: i

      if (.not. valid_context(collective)) then
         status = OUTPUT_COLLECTIVE_INVALID_CONTEXT
         return
      end if
      if (size(participants) == 0) then
         status = OUTPUT_COLLECTIVE_INVALID_PARTICIPANTS
         return
      end if
      if (participants(1) < 0 .or. participants(size(participants)) >= collective%rank_count) then
         status = OUTPUT_COLLECTIVE_INVALID_PARTICIPANTS
         return
      end if
      do i = 1, size(participants) - 1
         if (participants(i) >= participants(i + 1)) then
            status = OUTPUT_COLLECTIVE_INVALID_PARTICIPANTS
            return
         end if
      end do
      if (owner_rank /= participants(1)) then
         status = OUTPUT_COLLECTIVE_INVALID_OWNER
         return
      end if

      status = OUTPUT_COLLECTIVE_SUCCESS
   end subroutine validate_output_ownership

   pure subroutine select_output_publication_mode(collective, publication_mode)
      type(output_collective_t), intent(in) :: collective
      integer, intent(out) :: publication_mode

      if (collective%rank_count > 1 .and. collective%collective_publication_available) then
         publication_mode = OUTPUT_PUBLICATION_COLLECTIVE
      else
         publication_mode = OUTPUT_PUBLICATION_ROOT_AGGREGATION
      end if
   end subroutine select_output_publication_mode

   pure subroutine prepare_output_partition_publication(collective, participants, owner_rank, partition, &
                                                        local_participates, publication_mode, status)
      type(output_collective_t), intent(in) :: collective
      integer, intent(in) :: participants(:)
      integer, intent(in) :: owner_rank
      type(output_partition_t), intent(in) :: partition
      logical, intent(out) :: local_participates
      integer, intent(out) :: publication_mode, status
      integer :: ownership_status

      call validate_output_ownership(collective, participants, owner_rank, ownership_status)
      if (ownership_status /= OUTPUT_COLLECTIVE_SUCCESS) then
         local_participates = .false.
         publication_mode = OUTPUT_PUBLICATION_ROOT_AGGREGATION
         status = ownership_status
         return
      end if

      local_participates = any(participants == collective%rank)
      if (local_participates .neqv. partition%has_data) then
         publication_mode = OUTPUT_PUBLICATION_ROOT_AGGREGATION
         status = OUTPUT_COLLECTIVE_INVALID_PARTITION
         return
      end if

      call select_output_publication_mode(collective, publication_mode)
      status = OUTPUT_COLLECTIVE_SUCCESS
   end subroutine prepare_output_partition_publication

   pure logical function valid_context(collective)
      type(output_collective_t), intent(in) :: collective

      valid_context = collective%rank_count >= 1 .and. collective%rank >= 0 .and. &
                      collective%rank < collective%rank_count .and. collective%root_rank >= 0 .and. &
                      collective%root_rank < collective%rank_count
   end function valid_context

end module outputCollective_m
