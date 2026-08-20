module outputMetadata_m
    use outputTypes_m, only: probe_metadata_t, output_lifecycle_is_terminal, &
                             output_artifact_identity_is_valid, output_fragment_descriptor_is_valid, &
                             OUTPUT_ARTIFACT_UNDEFINED, OUTPUT_ARTIFACT_ROLE_CANONICAL, &
                             OUTPUT_ARTIFACT_ROLE_FRAGMENT, OUTPUT_LIFECYCLE_DECLARED, &
                             OUTPUT_LIFECYCLE_FAILED, TIME_DOMAIN, FREQUENCY_DOMAIN, BOTH_DOMAIN, UNDEFINED_DOMAIN
   use directoryUtils_m, only: atomic_replace_file, create_file_with_path, delete_file
   implicit none
   private

   integer, parameter, public :: OUTPUT_METADATA_SUCCESS = 0
   integer, parameter, public :: OUTPUT_METADATA_INVALID = 1
   integer, parameter, public :: OUTPUT_METADATA_IO_ERROR = 2

   public :: publish_initial_probe_metadata, publish_final_probe_metadata, json_escape

contains

   subroutine publish_initial_probe_metadata(path, metadata, status)
      character(len=*), intent(in) :: path
      type(probe_metadata_t), intent(in) :: metadata
      integer, intent(out) :: status
      type(probe_metadata_t) :: descriptor

      descriptor = metadata
      descriptor%lifecycle%state = OUTPUT_LIFECYCLE_DECLARED
      descriptor%lifecycle%diagnostic = ''
      call publish_probe_metadata(path, descriptor, .false., status)
   end subroutine publish_initial_probe_metadata

   subroutine publish_final_probe_metadata(path, metadata, status)
      character(len=*), intent(in) :: path
      type(probe_metadata_t), intent(in) :: metadata
      integer, intent(out) :: status

      call publish_probe_metadata(path, metadata, .true., status)
   end subroutine publish_final_probe_metadata

   subroutine publish_probe_metadata(path, metadata, terminal, status)
      character(len=*), intent(in) :: path
      type(probe_metadata_t), intent(in) :: metadata
      logical, intent(in) :: terminal
      integer, intent(out) :: status
       integer :: artifact_index, close_ios, fragment_index, ios, unit, write_ios
       character(len=:), allocatable :: temporary_path

      if (.not. valid_metadata(metadata, terminal)) then
         status = OUTPUT_METADATA_INVALID
         return
      end if

       temporary_path = trim(path)//'.tmp'
       call create_file_with_path(temporary_path, ios)
       if (ios /= 0) then
          status = OUTPUT_METADATA_IO_ERROR
          return
       end if
       open(newunit=unit, file=temporary_path, status='replace', action='write', iostat=ios)
       if (ios /= 0) then
          status = OUTPUT_METADATA_IO_ERROR
          return
      end if

       write(unit, '(a)', iostat=ios) '{'
       if (ios == 0) write(unit, '(a,i0,a)', iostat=ios) '"schema_version":', metadata%schema_version, ','
       if (ios == 0) write(unit, '(a)', iostat=ios) '"probe_id":"'//json_escape(trim(metadata%probe_id))//'",'
       if (ios == 0) write(unit, '(a)', iostat=ios) '"parent_probe_id":"'// &
          json_escape(trim(metadata%parent_probe_id))//'",'
       if (ios == 0) write(unit, '(a,i0,a)', iostat=ios) '"contributor_rank":', metadata%contributor_rank, ','
       if (ios == 0) write(unit, '(a)', iostat=ios) '"quantity":"'//json_escape(trim(metadata%quantity))//'",'
      if (ios == 0) write(unit, '(a,i0,a,i0,a,i0,a)', iostat=ios) '"lower_bound":{"x":', metadata%lower_bound%x, &
         ',"y":', metadata%lower_bound%y, ',"z":', metadata%lower_bound%z, '},'
      if (ios == 0) write(unit, '(a,i0,a,i0,a,i0,a)', iostat=ios) '"upper_bound":{"x":', metadata%upper_bound%x, &
         ',"y":', metadata%upper_bound%y, ',"z":', metadata%upper_bound%z, '},'
      if (ios == 0) write(unit, '(a)', iostat=ios) '"domain":"'//domain_name(metadata%domain_type)//'",'
      if (ios == 0) write(unit, '(a)', iostat=ios) '"lifecycle":{"state":"'//lifecycle_name(metadata%lifecycle%state)// &
         '","diagnostic":"'//json_escape(trim(metadata%lifecycle%diagnostic))//'"},'
      if (ios == 0) write(unit, '(a)', advance='no', iostat=ios) '"artifacts":['
      do artifact_index = 1, size(metadata%artifacts)
         if (artifact_index > 1 .and. ios == 0) write(unit, '(a)', advance='no', iostat=ios) ','
         if (ios == 0) call write_artifact(unit, metadata, artifact_index, ios)
       end do
       if (ios == 0) write(unit, '(a)', iostat=ios) '],'
       if (ios == 0) write(unit, '(a)', advance='no', iostat=ios) '"fragment_descriptors":['
       if (allocated(metadata%fragment_descriptors)) then
          do fragment_index = 1, size(metadata%fragment_descriptors)
             if (fragment_index > 1 .and. ios == 0) write(unit, '(a)', advance='no', iostat=ios) ','
             if (ios == 0) call write_fragment_descriptor(unit, metadata, fragment_index, ios)
          end do
       end if
       if (ios == 0) write(unit, '(a)', iostat=ios) '],'
       if (ios == 0) write(unit, '(a)', iostat=ios) '"ownership":{"participant_ranks":'// &
         participant_ranks_json(metadata)//',"scalar_writer_rank":'//scalar_writer_rank_json(metadata)//'}'
      if (ios == 0) write(unit, '(a)', iostat=ios) '}'
      write_ios = ios
      close(unit, iostat=close_ios)
       if (write_ios /= 0 .or. close_ios /= 0) then
          call delete_file(temporary_path, ios)
          status = OUTPUT_METADATA_IO_ERROR
          return
       end if
       call atomic_replace_file(temporary_path, path, ios)
       if (ios == 0) then
          status = OUTPUT_METADATA_SUCCESS
       else
          call delete_file(temporary_path, close_ios)
          status = OUTPUT_METADATA_IO_ERROR
       end if
   end subroutine publish_probe_metadata

   subroutine write_artifact(unit, metadata, artifact_index, ios)
      integer, intent(in) :: unit, artifact_index
      type(probe_metadata_t), intent(in) :: metadata
      integer, intent(inout) :: ios
      character(len=:), allocatable :: artifact

       artifact = '{"kind":"'//artifact_kind_name(metadata%artifacts(artifact_index)%kind)//'","role":"'// &
          artifact_role_name(metadata%artifacts(artifact_index)%role)//'","relative_path":"'// &
          json_escape(trim(metadata%artifacts(artifact_index)%relative_path))//'","required":'// &
         logical_json(metadata%artifacts(artifact_index)%required)//',"byte_order":"'// &
         byte_order_name(metadata%artifacts(artifact_index)%byte_order)//'","numeric_representation":"'// &
          numeric_name(metadata%artifacts(artifact_index)%numeric_representation)//'","complex_representation":"'// &
          complex_name(metadata%artifacts(artifact_index)%complex_representation)//'","record_bytes":'// &
          integer_json(metadata%artifacts(artifact_index)%record_bytes)//',"component_order":"'// &
          json_escape(trim(metadata%artifacts(artifact_index)%component_order))//'"}'
       write(unit, '(a)', advance='no', iostat=ios) artifact
    end subroutine write_artifact

    subroutine write_fragment_descriptor(unit, metadata, descriptor_index, ios)
       integer, intent(in) :: unit, descriptor_index
       type(probe_metadata_t), intent(in) :: metadata
       integer, intent(inout) :: ios
       character(len=:), allocatable :: descriptor

       descriptor = '{"parent_probe_id":"'// &
          json_escape(trim(metadata%fragment_descriptors(descriptor_index)%identity%parent_probe_id))// &
          '","contributor_rank":'// &
          integer_json(int(metadata%fragment_descriptors(descriptor_index)%identity%contributor_rank, kind=8))// &
          ',"relative_path":"'//json_escape(trim(metadata%fragment_descriptors(descriptor_index)%relative_path))//'"}'
       write(unit, '(a)', advance='no', iostat=ios) descriptor
    end subroutine write_fragment_descriptor

   logical function valid_metadata(metadata, terminal)
      type(probe_metadata_t), intent(in) :: metadata
      logical, intent(in) :: terminal
       integer :: i, j

      valid_metadata = metadata%schema_version > 0 .and. len_trim(metadata%probe_id) > 0 .and. &
         len_trim(metadata%quantity) > 0 .and. allocated(metadata%artifacts)
      if (.not. valid_metadata) return
      valid_metadata = size(metadata%artifacts) > 0
      if (.not. valid_metadata) return
      if (terminal) valid_metadata = output_lifecycle_is_terminal(metadata%lifecycle)
      if (metadata%lifecycle%state == OUTPUT_LIFECYCLE_FAILED) then
         valid_metadata = valid_metadata .and. len_trim(metadata%lifecycle%diagnostic) > 0
      end if
      if (.not. valid_metadata) return
       do i = 1, size(metadata%artifacts)
          if (metadata%artifacts(i)%kind == OUTPUT_ARTIFACT_UNDEFINED .or. &
              .not. relative_path(metadata%artifacts(i)%relative_path) .or. &
              .not. output_artifact_identity_is_valid(metadata%artifacts(i))) then
             valid_metadata = .false.
             return
          end if
       end do

       if (len_trim(metadata%parent_probe_id) == 0) then
          if (metadata%contributor_rank /= -1 .or. any(metadata%artifacts%role /= OUTPUT_ARTIFACT_ROLE_CANONICAL)) then
             valid_metadata = .false.
             return
          end if
          if (allocated(metadata%fragment_descriptors)) then
             do i = 1, size(metadata%fragment_descriptors)
                if (.not. output_fragment_descriptor_is_valid(metadata%fragment_descriptors(i), metadata%probe_id) .or. &
                    .not. relative_path(metadata%fragment_descriptors(i)%relative_path)) then
                   valid_metadata = .false.
                   return
                end if
                do j = 1, i - 1
                   if (metadata%fragment_descriptors(j)%identity%contributor_rank == &
                       metadata%fragment_descriptors(i)%identity%contributor_rank) then
                      valid_metadata = .false.
                      return
                   end if
                end do
             end do
          end if
       else
          if (metadata%contributor_rank < 0 .or. any(metadata%artifacts%role /= OUTPUT_ARTIFACT_ROLE_FRAGMENT)) then
             valid_metadata = .false.
             return
          end if
          if (allocated(metadata%fragment_descriptors)) then
             if (size(metadata%fragment_descriptors) /= 0) then
                valid_metadata = .false.
                return
             end if
          end if
          do i = 1, size(metadata%artifacts)
             if (trim(metadata%artifacts(i)%fragment%parent_probe_id) /= trim(metadata%parent_probe_id) .or. &
                 metadata%artifacts(i)%fragment%contributor_rank /= metadata%contributor_rank) then
                valid_metadata = .false.
                return
             end if
          end do
       end if
    end function valid_metadata

   pure logical function relative_path(path)
      character(len=*), intent(in) :: path
      character(len=:), allocatable :: value

      value = trim(path)
      relative_path = len(value) > 0
      if (.not. relative_path) return
       relative_path = value(1:1) /= '/' .and. value(1:1) /= '\'
      if (len(value) > 1) relative_path = relative_path .and. value(2:2) /= ':'
   end function relative_path

   pure function json_escape(value) result(escaped)
      character(len=*), intent(in) :: value
      character(len=:), allocatable :: escaped
      integer :: i

      escaped = ''
      do i = 1, len_trim(value)
         select case (value(i:i))
         case ('"')
            escaped = escaped//'\"'
         case ('\')
            escaped = escaped//'\\'
         case (achar(9))
            escaped = escaped//'\t'
         case (achar(10))
            escaped = escaped//'\n'
         case (achar(13))
            escaped = escaped//'\r'
         case default
            if (iachar(value(i:i)) >= 32) escaped = escaped//value(i:i)
         end select
      end do
   end function json_escape

   pure function lifecycle_name(state) result(name)
      integer, intent(in) :: state
      character(len=:), allocatable :: name

      select case (state)
      case (0); name = 'declared'
      case (1); name = 'active'
      case (2); name = 'finalising'
      case (3); name = 'complete'
      case (4); name = 'failed'
      case default; name = 'unknown'
      end select
   end function lifecycle_name

   pure function domain_name(domain_type) result(name)
      integer, intent(in) :: domain_type
      character(len=:), allocatable :: name

      select case (domain_type)
      case (TIME_DOMAIN); name = 'time'
      case (FREQUENCY_DOMAIN); name = 'frequency'
      case (BOTH_DOMAIN); name = 'both'
      case (UNDEFINED_DOMAIN); name = 'undefined'
      case default; name = 'unknown'
      end select
   end function domain_name

    pure function artifact_kind_name(kind) result(name)
      integer, intent(in) :: kind
      character(len=:), allocatable :: name

      select case (kind)
      case (1); name = 'text'
      case (2); name = 'binary'
      case (3); name = 'metadata'
      case (4); name = 'visualisation_metadata'
      case (5); name = 'visualisation_data'
      case (6); name = 'geometry'
      case default; name = 'undefined'
      end select
    end function artifact_kind_name

    pure function artifact_role_name(role) result(name)
       integer, intent(in) :: role
       character(len=:), allocatable :: name

       select case (role)
       case (OUTPUT_ARTIFACT_ROLE_CANONICAL); name = 'canonical'
       case (OUTPUT_ARTIFACT_ROLE_FRAGMENT); name = 'fragment'
       case default; name = 'undefined'
       end select
    end function artifact_role_name

    pure function byte_order_name(value) result(name)
      integer, intent(in) :: value
      character(len=:), allocatable :: name

      select case (value)
      case (1); name = 'little_endian'
      case (2); name = 'big_endian'
      case default; name = 'unspecified'
      end select
   end function byte_order_name

   pure function numeric_name(value) result(name)
      integer, intent(in) :: value
      character(len=:), allocatable :: name

      select case (value)
      case (1); name = 'real32'
      case (2); name = 'real64'
      case (3); name = 'int32'
      case (4); name = 'int64'
      case default; name = 'unspecified'
      end select
   end function numeric_name

   pure function complex_name(value) result(name)
      integer, intent(in) :: value
      character(len=:), allocatable :: name

      select case (value)
      case (1); name = 'real_imag'
      case (2); name = 'magnitude_phase'
      case default; name = 'unspecified'
      end select
   end function complex_name

   pure function logical_json(value) result(text)
      logical, intent(in) :: value
      character(len=:), allocatable :: text

      if (value) then
         text = 'true'
      else
         text = 'false'
      end if
   end function logical_json

   function integer_json(value) result(text)
      integer(kind=8), intent(in) :: value
      character(len=:), allocatable :: text
      character(len=32) :: buffer

      write(buffer, '(i0)') value
      text = trim(buffer)
   end function integer_json

   function participant_ranks_json(metadata) result(text)
      type(probe_metadata_t), intent(in) :: metadata
      character(len=:), allocatable :: text
      character(len=32) :: buffer
      integer :: i

      text = '['
      if (allocated(metadata%ownership%participant_ranks)) then
         do i = 1, size(metadata%ownership%participant_ranks)
            if (i > 1) text = text//','
            write(buffer, '(i0)') metadata%ownership%participant_ranks(i)
            text = text//trim(buffer)
         end do
      end if
      text = text//']'
   end function participant_ranks_json

   function scalar_writer_rank_json(metadata) result(text)
      type(probe_metadata_t), intent(in) :: metadata
      character(len=:), allocatable :: text
      character(len=32) :: buffer

      write(buffer, '(i0)') metadata%ownership%scalar_writer_rank
      text = trim(buffer)
   end function scalar_writer_rank_json

end module outputMetadata_m
