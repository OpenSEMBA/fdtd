module mod_frequencySliceProbeOutput
   use FDETYPES_m
   use utils_m
   use Report
   use outputTypes
   use mod_outputUtils
   use mod_volumicProbeUtils
   use directoryUtils_m
   use HDF5
   use mod_xdmfAPI
   implicit none
   private

   !===========================
   !  Public interface summary
   !===========================
   public :: init_frequency_slice_probe_output
   public :: update_frequency_slice_probe_output
   public :: flush_frequency_slice_probe_output
   public :: close_frequency_slice_probe_output
   !===========================

   !===========================
   !  Private interface summary
   !===========================
   private :: save_field
   private :: save_field_module
   private :: save_field_component
   private :: save_current
   private :: save_current_module
   private :: save_current_component
   !===========================

   !===========================

contains

   subroutine init_frequency_slice_probe_output(this, lowerBound, upperBound, timeInterval, field, domain, outputTypeExtension, control, problemInfo)
      type(frequency_slice_probe_output_t), intent(out) :: this
      type(cell_coordinate_t), intent(in) :: lowerBound, upperBound
      real(kind=RKIND_tiempo), intent(in) :: timeInterval
      integer(kind=SINGLE), intent(in) :: field
      type(domain_t), intent(in) :: domain
      character(len=BUFSIZE), intent(in) :: outputTypeExtension
      type(sim_control_t), intent(in) :: control
      type(problem_info_t), intent(in) :: problemInfo

      integer :: i
      integer :: error
      character(len=BUFSIZE) :: filename

      this%mainCoords = lowerBound
      this%auxCoords = upperBound
      this%component = field !This can refer to electric, magnetic or currentDensity
      this%domain = domain
      this%nFreq = domain%fnum

      call alloc_and_init(this%frequencySlice, this%nFreq, 0.0_RKIND)
      do i = 1, this%nFreq
         call init_frequency_slice(this%frequencySlice, this%domain)
      end do

      call find_and_store_important_coords(this%mainCoords, this%auxCoords, this%component, problemInfo, this%nPoints, this%coords)

      call alloc_and_init(this%xValueForFreq, this%nFreq, this%nPoints, (0.0_CKIND, 0.0_CKIND))
      call alloc_and_init(this%yValueForFreq, this%nFreq, this%nPoints, (0.0_CKIND, 0.0_CKIND))
      call alloc_and_init(this%zValueForFreq, this%nFreq, this%nPoints, (0.0_CKIND, 0.0_CKIND))

      call alloc_and_init(this%auxExp_E, this%nFreq, (0.0_CKIND, 0.0_CKIND))
      call alloc_and_init(this%auxExp_H, this%nFreq, (0.0_CKIND, 0.0_CKIND))

      do i = 1, this%nFreq
         this%auxExp_E(i) = timeInterval*(1.0E0_RKIND, 0.0E0_RKIND)*Exp(mcpi2*this%frequencySlice(i))   ! the dt should be some kind of average
         this%auxExp_H(i) = this%auxExp_E(i)*Exp(mcpi2*this%frequencySlice(i)*timeInterval*0.5_RKIND)
      end do

      this%path = get_output_path_freq(this, outputTypeExtension, field, control)

      filename = get_last_component(this%path)
      this%filesPath = join_path(this%path, filename)

      call create_folder(this%path, error)
      call create_bin_file(this%filesPath, error)
   end subroutine init_frequency_slice_probe_output

   subroutine create_bin_file(filePath, error)
      character(len=*), intent(in) :: filePath
      integer, intent(out) :: error
      call create_file_with_path(add_extension(filePath, binaryExtension), error)
   end subroutine 

   subroutine write_bin_file(this)
      ! Check type definition for binary format
      type(frequency_slice_probe_output_t), intent(inout) :: this
      integer :: i, f, unit
      !We rewrite the binary as simulation continues
      open (unit=unit, file=add_extension(this%filesPath, binaryExtension), &
            status='old', form='unformatted', action='write', access='stream')
      do f = 1, this%nFreq
      do i = 1, this%nPoints
         write(unit) this%frequencySlice(f), this%coords(1,i), this%coords(2,i), this%coords(3,i), this%xValueForFreq(f,i), this%yValueForFreq(f,i), this%xValueForFreq(f,i)
      end do
      end do
      flush (unit)
      close (unit)
   end subroutine

   subroutine write_to_xdmf_h5(this)
      !If we call to this subrutine it will always replace old values
      type(frequency_slice_probe_output_t), intent(inout) :: this

      integer(HID_T) :: file_id
      integer :: f, error, xdmfunit
      real(dp), allocatable, dimension(:, :) :: coordsReal
      character(len=256) :: h5_filename, h5_filepath
      character(len=256) :: xdmf_filename
      character(len=10) :: dimension_string
      character(len=10) :: nCoordsString
      character(len=14) :: stepName
      h5_filepath = add_extension(this%filesPath, ".h5")
      h5_filename = get_last_component(h5_filepath)

      call H5open_f(error)
      call H5Fopen_f(trim(h5_filepath), H5F_ACC_RDWR_F, file_id, error)
      write(dimension_string, '(I0,1X,I0)') this%nPoints, this%nFreq
      write(nCoordsString, '(I0, I0)') this%nPoints, 1
      do f = 1, this%nFreq
         write(stepName, '((I5.5))') f
         call h5_append_rows_to_dataset(file_id, "xVal", reshape(real(abs(this%xValueForFreq(f, :)), dp), [1, this%nPoints]))
         call h5_append_rows_to_dataset(file_id, "yVal", reshape(real(abs(this%yValueForFreq(f, :)), dp), [1, this%nPoints]))
         call h5_append_rows_to_dataset(file_id, "zVal", reshape(real(abs(this%zValueForFreq(f, :)), dp), [1, this%nPoints]))
      end do
      call H5Fclose_f(file_id, error)
      call H5close_f(error)
   end subroutine

   function get_output_path_freq(this, outputTypeExtension, field, control) result(outputPath)
      type(frequency_slice_probe_output_t), intent(in) :: this
      character(len=*), intent(in) :: outputTypeExtension
      integer(kind=SINGLE), intent(in) :: field
      type(sim_control_t), intent(in) :: control
      character(len=BUFSIZE)  :: probeBoundsExtension, prefixFieldExtension
      character(len=BUFSIZE) :: outputPath
      probeBoundsExtension = get_coordinates_extension(this%mainCoords, this%auxCoords, control%mpidir)
      prefixFieldExtension = get_prefix_extension(field, control%mpidir)
      outputPath = &
         trim(adjustl(outputTypeExtension))//'_'//trim(adjustl(prefixFieldExtension))//'_'//trim(adjustl(probeBoundsExtension))
   end function get_output_path_freq

   subroutine update_frequency_slice_probe_output(this, step, fieldsReference, control, problemInfo)
      type(frequency_slice_probe_output_t), intent(inout) :: this
      real(kind=RKIND_tiempo), intent(in) :: step
      type(sim_control_t), intent(in) :: control
      type(problem_info_t), intent(in) :: problemInfo
      type(fields_reference_t), intent(in) :: fieldsReference

      integer(kind=4) :: request
      request = this%component

      if (any(VOLUMIC_M_MEASURE == request)) then
         select case (request)
         case (iCur); call save_current_module(this, fieldsReference, step, problemInfo)
         case (iMEC); call save_field_module(this, fieldsReference%E, step, request, problemInfo)
         case (iMHC); call save_field_module(this, fieldsReference%H, step, request, problemInfo)
         case default; call StopOnError(control%layoutnumber, control%size, "Volumic measure not supported")
         end select

      else if (any(VOLUMIC_X_MEASURE == request)) then
         select case (request)
         case (iCurX); call save_current_component(this, this%xValueForFreq, fieldsReference, problemInfo, iEx, this%auxExp_E, this%nFreq, step)
         case (iExC); call save_field_component(this, this%xValueForFreq, fieldsReference%E%x, step, problemInfo, iEx)
         case (iHxC); call save_field_component(this, this%xValueForFreq, fieldsReference%H%x, step, problemInfo, iHx)
         case default; call StopOnError(control%layoutnumber, control%size, "Volumic measure not supported")
         end select

      else if (any(VOLUMIC_Y_MEASURE == request)) then
         select case (request)
         case (iCurY); call save_current_component(this, this%yValueForFreq, fieldsReference, problemInfo, iEy, this%auxExp_E, this%nFreq, step)
         case (iEyC); call save_field_component(this, this%yValueForFreq, fieldsReference%E%y, step, problemInfo, iEy)
         case (iHyC); call save_field_component(this, this%yValueForFreq, fieldsReference%H%y, step, problemInfo, iHy)
         case default; call StopOnError(control%layoutnumber, control%size, "Volumic measure not supported")
         end select

      else if (any(VOLUMIC_Z_MEASURE == request)) then
         select case (request)
         case (iCurZ); call save_current_component(this, this%zValueForFreq, fieldsReference, problemInfo, iEz, this%auxExp_E, this%nFreq, step)
         case (iEzC); call save_field_component(this, this%zValueForFreq, fieldsReference%E%z, step, problemInfo, iEz)
         case (iHzC); call save_field_component(this, this%zValueForFreq, fieldsReference%H%z, step, problemInfo, iHz)
         case default; call StopOnError(control%layoutnumber, control%size, "Volumic measure not supported")
         end select
      end if
   end subroutine update_frequency_slice_probe_output

   subroutine save_current_module(this, fieldsReference, step, problemInfo)
      type(frequency_slice_probe_output_t), intent(inout) :: this
      type(fields_reference_t), intent(in) :: fieldsReference
      type(problem_info_t), intent(in) :: problemInfo
      real(kind=RKIND_tiempo), intent(in) :: step

      integer :: i, j, k, coordIdx

      coordIdx = 0
      do i = this%mainCoords%x, this%auxCoords%x
      do j = this%mainCoords%y, this%auxCoords%y
      do k = this%mainCoords%z, this%auxCoords%z
         if (isValidPointForCurrent(iCur, i, j, k, problemInfo)) then
            coordIdx = coordIdx + 1
            call save_current(this%xValueForFreq, iEx, coordIdx, i, j, k, fieldsReference, this%auxExp_E, this%nFreq, step)
            call save_current(this%yValueForFreq, iEy, coordIdx, i, j, k, fieldsReference, this%auxExp_E, this%nFreq, step)
            call save_current(this%zValueForFreq, iEz, coordIdx, i, j, k, fieldsReference, this%auxExp_E, this%nFreq, step)
         end if
      end do
      end do
      end do
   end subroutine

   subroutine save_current_component(this, currentData, fieldsReference, problemInfo, fieldDir, auxExp, nFreq, step)
      type(frequency_slice_probe_output_t), intent(inout) :: this
      complex(kind=CKIND), intent(inout) :: currentData(:, :)
      type(fields_reference_t), intent(in) :: fieldsReference
      type(problem_info_t), intent(in) :: problemInfo
      integer, intent(in) :: fieldDir, nFreq
      complex(kind=ckind), intent(in), dimension(:) :: auxExp
      real(kind=RKIND_tiempo), intent(in) :: step

      integer :: i, j, k, coordIdx

      coordIdx = 0
      do i = this%mainCoords%x, this%auxCoords%x
      do j = this%mainCoords%y, this%auxCoords%y
      do k = this%mainCoords%z, this%auxCoords%z
         if (isValidPointForCurrent(fieldDir, i, j, k, problemInfo)) then
            coordIdx = coordIdx + 1
            call save_current(currentData, fieldDir, coordIdx, i, j, k, fieldsReference, auxExp, nFreq, step)
         end if
      end do
      end do
      end do
   end subroutine

   subroutine save_current(valorComplex, direction, coordIdx, i, j, k, fieldsReference, auxExponential, nFreq, step)
      integer, intent(in) :: direction
      complex(kind=CKIND), intent(inout) :: valorComplex(:, :)
      complex(kind=CKIND), intent(in) :: auxExponential(:)
      integer, intent(in) :: i, j, k, coordIdx, nFreq
      type(fields_reference_t), intent(in) :: fieldsReference
      real(kind=RKIND_tiempo), intent(in) :: step

      integer :: iter
      complex(kind=CKIND) :: z_cplx = (0.0_RKIND, 0.0_RKIND)
      real(kind=rkind) :: jdir

      jdir = computej(direction, i, j, k, fieldsReference)

      do iter = 1, nFreq
         valorComplex(iter, coordIdx) = valorComplex(iter, coordIdx) + (auxExponential(iter)**step)*jdir
      end do
   end subroutine

   subroutine save_field_module(this, fieldInfo, simTime, request, problemInfo)
      type(frequency_slice_probe_output_t), intent(inout) :: this
      type(field_data_t), intent(in) :: fieldInfo
      real(kind=RKIND_tiempo), intent(in) :: simTime
      type(problem_info_t), intent(in) :: problemInfo
      integer, intent(in) :: request

      complex(kind=CKIND), dimension(this%nFreq) :: auxExponential
      integer :: i, j, k, coordIdx

      if (iMHC == request) auxExponential = this%auxExp_H**simTime
      if (iMEC == request) auxExponential = this%auxExp_E**simTime

      coordIdx = 0
      do i = this%mainCoords%x, this%auxCoords%x
      do j = this%mainCoords%y, this%auxCoords%y
      do k = this%mainCoords%z, this%auxCoords%z
         if (isValidPointForField(request, i, j, k, problemInfo)) then
            coordIdx = coordIdx + 1
            call save_field(this%xValueForFreq, auxExponential, fieldInfo%x(i, j, k), this%nFreq, coordIdx)
            call save_field(this%yValueForFreq, auxExponential, fieldInfo%y(i, j, k), this%nFreq, coordIdx)
            call save_field(this%zValueForFreq, auxExponential, fieldInfo%z(i, j, k), this%nFreq, coordIdx)
         end if
      end do
      end do
      end do

   end subroutine

   subroutine save_field_component(this, fieldData, fieldComponent, simTime, problemInfo, fieldDir)
      type(frequency_slice_probe_output_t), intent(inout) :: this
      complex(kind=CKIND), intent(inout) :: fieldData(:, :)
      real(kind=RKIND), intent(in) :: fieldComponent(:, :, :)
      real(kind=RKIND_tiempo), intent(in) :: simTime
      type(problem_info_t), intent(in) :: problemInfo
      integer, intent(in) :: fieldDir

      complex(kind=CKIND), dimension(this%nFreq) :: auxExponential
      integer :: i, j, k, coordIdx

      if (any(MAGNETIC_FIELD_DIRECTION == fieldDir)) auxExponential = this%auxExp_H**simTime
      if (any(ELECTRIC_FIELD_DIRECTION == fieldDir)) auxExponential = this%auxExp_E**simTime

      coordIdx = 0
      do i = this%mainCoords%x, this%auxCoords%x
      do j = this%mainCoords%y, this%auxCoords%y
      do k = this%mainCoords%z, this%auxCoords%z
         if (isValidPointForField(fieldDir, i, j, k, problemInfo)) then
            coordIdx = coordIdx + 1
            call save_field(fieldData, auxExponential, fieldComponent(i, j, k), this%nFreq, coordIdx)
         end if
      end do
      end do
      end do
   end subroutine

   subroutine save_field(valorComplex, auxExp, fieldValue, nFreq, coordIdx)
      complex(kind=CKIND), intent(inout) :: valorComplex(:, :)
      complex(kind=CKIND), intent(in) :: auxExp(:)
      real(KIND=RKIND), intent(in) :: fieldValue
      integer(KIND=SINGLE), intent(in) :: nFreq, coordIdx

      integer :: freq

      do freq = 1, nFreq
         valorComplex = valorComplex(freq, coordIdx) + auxExp(freq)*fieldValue
      end do
   end subroutine

   subroutine flush_frequency_slice_probe_output(this)
      type(frequency_slice_probe_output_t), intent(inout) :: this
      call write_bin_file(this)
   end subroutine flush_frequency_slice_probe_output

   subroutine close_frequency_slice_probe_output(this)
      type(frequency_slice_probe_output_t), intent(inout) :: this
      call write_to_xdmf_h5(this)
   end subroutine



end module mod_frequencySliceProbeOutput
