module outputBinary_m
   use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
   use outputTypes_m, only: output_artifact_t, OUTPUT_ARTIFACT_BINARY, BINARY_ENDIAN_LITTLE, &
                            BINARY_NUMERIC_REAL32, BINARY_NUMERIC_REAL64, &
                            BINARY_COMPLEX_UNSPECIFIED, BINARY_COMPLEX_REAL_IMAG, &
                            BINARY_BYTES_REAL32, BINARY_BYTES_REAL64
    use directoryUtils_m, only: create_file_with_path, file_exists
   implicit none
   private

   integer, parameter, public :: BINARY_WRITER_SUCCESS = 0
   integer, parameter, public :: BINARY_WRITER_INVALID_LAYOUT = 1
   integer, parameter, public :: BINARY_WRITER_SIZE_MISMATCH = 2
   integer, parameter, public :: BINARY_WRITER_IO_ERROR = 3

    public :: validate_binary_layout, open_binary_append, write_binary_real32, write_binary_real64
    public :: write_binary_complex32, write_binary_complex64, write_binary_complex_record32

contains

    subroutine validate_binary_layout(artifact, status)
      type(output_artifact_t), intent(in) :: artifact
      integer, intent(out) :: status
      integer(kind=8) :: scalar_bytes

      status = BINARY_WRITER_INVALID_LAYOUT
      if (artifact%kind /= OUTPUT_ARTIFACT_BINARY) return
      if (artifact%byte_order /= BINARY_ENDIAN_LITTLE) return
      select case (artifact%numeric_representation)
      case (BINARY_NUMERIC_REAL32)
         scalar_bytes = BINARY_BYTES_REAL32
      case (BINARY_NUMERIC_REAL64)
         scalar_bytes = BINARY_BYTES_REAL64
      case default
         return
      end select
      select case (artifact%complex_representation)
      case (BINARY_COMPLEX_UNSPECIFIED)
      case (BINARY_COMPLEX_REAL_IMAG)
         scalar_bytes = 2 * scalar_bytes
      case default
         return
      end select
      if (artifact%record_bytes < scalar_bytes) return
      if (mod(artifact%record_bytes, scalar_bytes) /= 0) return
      status = BINARY_WRITER_SUCCESS
    end subroutine validate_binary_layout

    subroutine open_binary_append(path, artifact, unit, status)
       character(len=*), intent(in) :: path
       type(output_artifact_t), intent(in) :: artifact
       integer, intent(out) :: unit, status
       integer :: ios

       unit = -1
       call validate_binary_layout(artifact, status)
       if (status /= BINARY_WRITER_SUCCESS) return
        if (.not. file_exists(path)) then
           call create_file_with_path(path, ios)
           if (ios /= 0) then
              status = BINARY_WRITER_IO_ERROR
              return
           end if
        end if
       open(newunit=unit, file=trim(path), access='stream', form='unformatted', status='old', &
            position='append', action='write', iostat=ios)
       if (ios == 0) then
          status = BINARY_WRITER_SUCCESS
       else
          status = BINARY_WRITER_IO_ERROR
       end if
    end subroutine open_binary_append

   subroutine write_binary_real32(path, artifact, values, status)
      character(len=*), intent(in) :: path
      type(output_artifact_t), intent(in) :: artifact
      real(real32), intent(in) :: values(:)
      integer, intent(out) :: status

      call write_real32_values(path, artifact, values, BINARY_COMPLEX_UNSPECIFIED, status)
   end subroutine write_binary_real32

   subroutine write_binary_real64(path, artifact, values, status)
      character(len=*), intent(in) :: path
      type(output_artifact_t), intent(in) :: artifact
      real(real64), intent(in) :: values(:)
      integer, intent(out) :: status

      call write_real64_values(path, artifact, values, BINARY_COMPLEX_UNSPECIFIED, status)
   end subroutine write_binary_real64

    subroutine write_binary_complex32(path, artifact, values, status)
      character(len=*), intent(in) :: path
      type(output_artifact_t), intent(in) :: artifact
      complex(real32), intent(in) :: values(:)
      integer, intent(out) :: status
      real(real32), allocatable :: scalars(:)
      integer :: i

      allocate(scalars(2 * size(values)))
      do i = 1, size(values)
         scalars(2 * i - 1) = real(values(i), real32)
         scalars(2 * i) = aimag(values(i))
      end do
      call write_real32_values(path, artifact, scalars, BINARY_COMPLEX_REAL_IMAG, status)
    end subroutine write_binary_complex32

    subroutine write_binary_complex_record32(path, artifact, values, status)
       character(len=*), intent(in) :: path
       type(output_artifact_t), intent(in) :: artifact
       real(real32), intent(in) :: values(:)
       integer, intent(out) :: status

       call write_real32_values(path, artifact, values, BINARY_COMPLEX_REAL_IMAG, status)
    end subroutine write_binary_complex_record32

   subroutine write_binary_complex64(path, artifact, values, status)
      character(len=*), intent(in) :: path
      type(output_artifact_t), intent(in) :: artifact
      complex(real64), intent(in) :: values(:)
      integer, intent(out) :: status
      real(real64), allocatable :: scalars(:)
      integer :: i

      allocate(scalars(2 * size(values)))
      do i = 1, size(values)
         scalars(2 * i - 1) = real(values(i), real64)
         scalars(2 * i) = aimag(values(i))
      end do
      call write_real64_values(path, artifact, scalars, BINARY_COMPLEX_REAL_IMAG, status)
   end subroutine write_binary_complex64

   subroutine write_real32_values(path, artifact, values, complex_representation, status)
      character(len=*), intent(in) :: path
      type(output_artifact_t), intent(in) :: artifact
      real(real32), intent(in) :: values(:)
      integer, intent(in) :: complex_representation
      integer, intent(out) :: status
      integer :: i, ios, unit, write_ios

      call prepare_write(path, artifact, BINARY_NUMERIC_REAL32, complex_representation, &
                         int(size(values), kind=8) * BINARY_BYTES_REAL32, status)
      if (status /= BINARY_WRITER_SUCCESS) return
      open(newunit=unit, file=trim(path), access='stream', form='unformatted', status='replace', &
           action='write', iostat=ios)
      if (ios /= 0) then
         status = BINARY_WRITER_IO_ERROR
         return
      end if
      do i = 1, size(values)
         call write_int32_little_endian(unit, transfer(values(i), 0_int32), ios)
         if (ios /= 0) exit
      end do
      write_ios = ios
      close(unit, iostat=ios)
      if (write_ios == 0 .and. ios == 0) then
         status = BINARY_WRITER_SUCCESS
      else
         status = BINARY_WRITER_IO_ERROR
      end if
   end subroutine write_real32_values

   subroutine write_real64_values(path, artifact, values, complex_representation, status)
      character(len=*), intent(in) :: path
      type(output_artifact_t), intent(in) :: artifact
      real(real64), intent(in) :: values(:)
      integer, intent(in) :: complex_representation
      integer, intent(out) :: status
      integer :: i, ios, unit, write_ios

      call prepare_write(path, artifact, BINARY_NUMERIC_REAL64, complex_representation, &
                         int(size(values), kind=8) * BINARY_BYTES_REAL64, status)
      if (status /= BINARY_WRITER_SUCCESS) return
      open(newunit=unit, file=trim(path), access='stream', form='unformatted', status='replace', &
           action='write', iostat=ios)
      if (ios /= 0) then
         status = BINARY_WRITER_IO_ERROR
         return
      end if
      do i = 1, size(values)
         call write_int64_little_endian(unit, transfer(values(i), 0_int64), ios)
         if (ios /= 0) exit
      end do
      write_ios = ios
      close(unit, iostat=ios)
      if (write_ios == 0 .and. ios == 0) then
         status = BINARY_WRITER_SUCCESS
      else
         status = BINARY_WRITER_IO_ERROR
      end if
   end subroutine write_real64_values

   subroutine prepare_write(path, artifact, numeric_representation, complex_representation, byte_count, status)
      character(len=*), intent(in) :: path
      type(output_artifact_t), intent(in) :: artifact
      integer, intent(in) :: numeric_representation, complex_representation
      integer(kind=8), intent(in) :: byte_count
      integer, intent(out) :: status
      integer :: ios

      call validate_binary_layout(artifact, status)
      if (status /= BINARY_WRITER_SUCCESS) return
      if (artifact%numeric_representation /= numeric_representation .or. &
          artifact%complex_representation /= complex_representation) then
         status = BINARY_WRITER_INVALID_LAYOUT
         return
      end if
      if (mod(byte_count, artifact%record_bytes) /= 0) then
         status = BINARY_WRITER_SIZE_MISMATCH
         return
      end if
      call create_file_with_path(path, ios)
      if (ios /= 0) status = BINARY_WRITER_IO_ERROR
   end subroutine prepare_write

   subroutine write_int32_little_endian(unit, value, ios)
      integer, intent(in) :: unit
      integer(int32), intent(in) :: value
      integer, intent(inout) :: ios
      integer :: byte_index

      do byte_index = 0, BINARY_BYTES_REAL32 - 1
         write(unit, iostat=ios) achar(ibits(value, 8 * byte_index, 8))
         if (ios /= 0) return
      end do
   end subroutine write_int32_little_endian

   subroutine write_int64_little_endian(unit, value, ios)
      integer, intent(in) :: unit
      integer(int64), intent(in) :: value
      integer, intent(inout) :: ios
      integer :: byte_index

      do byte_index = 0, BINARY_BYTES_REAL64 - 1
         write(unit, iostat=ios) achar(ibits(value, 8 * byte_index, 8))
         if (ios /= 0) return
      end do
   end subroutine write_int64_little_endian

end module outputBinary_m
