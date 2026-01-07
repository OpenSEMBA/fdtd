module mod_farFieldOutput
   use outputTypes
   use Report
   use mod_outputUtils
   use farfield_m
   implicit none
   private
   !===========================
   !  Public interface summary
   !===========================
   public :: init_farField_probe_output
   !public :: create_farField_probe_output
   public :: update_farField_probe_output
   public :: flush_farField_probe_output
   !===========================
contains

  subroutine init_farField_probe_output(this, sgg, lowerBound, upperBound, field, domain, sphericRange, control, outputTypeExtension, fileNormalize, eps0, mu0, geometricMedia, SINPML_fullsize, bounds)
      type(far_field_probe_output_t), intent(out) :: this
      type(domain_t), intent(in) :: domain
      type(SGGFDTDINFO), intent(in) :: sgg
      type(cell_coordinate_t), intent(in) :: lowerBound, upperBound
      integer(kind=SINGLE), intent(in) :: field
      type(spheric_domain_t), intent(in) :: sphericRange
      type(sim_control_t), intent(in) :: control
      type(media_matrices_t), intent(in) :: geometricMedia
      type(limit_t), dimension(:), intent(in)  ::  SINPML_fullsize
      character(len=*), intent(in) :: fileNormalize, outputTypeExtension
      real(kind=RKIND), intent(in) :: mu0, eps0
      type(bounds_t), intent(in) :: bounds

      if (domain%domainType /= TIME_DOMAIN) call StopOnError(0, 0, "Unexpected domain type for farField probe")

      this%domain = domain
      this%sphericRange = sphericRange
      this%component = field
      this%path = get_output_path()
      this%fileUnitFreq = 2025 !Dummy unit for now

      call InitFarField(sgg, &
         geometricMedia%sggMiEx,geometricMedia%sggMiEy,geometricMedia%sggMiEz,geometricMedia%sggMiHx,geometricMedia%sggMiHy,geometricMedia%sggMiHz, &
                        control%layoutnumber, control%size, bounds, control%resume, &
                        this%fileUnitFreq, this%path, &
                        lowerBound%x, upperBound%x, &
                        lowerBound%y, upperBound%y, &
                        lowerBound%z, upperBound%z, &
                        domain%fstart, domain%fstop, domain%fstep, &
                        sphericRange%phiStart, sphericRange%phiStop, sphericRange%phiStep, &
                        sphericRange%thetaStart, sphericRange%thetaStop, sphericRange%thetaStep, &
                        fileNormalize, SINPML_fullsize, &
                        control%facesNF2FF, control%NF2FFDecim, &
#ifdef CompileWithMPI
                        output(ii)%item(i)%MPISubComm, output(ii)%item(i)%MPIRoot, &
#endif
                        eps0, mu0)

   contains
      function get_output_path() result(outputPath)
         character(len=BUFSIZE)  :: probeBoundsExtension, prefixFieldExtension
         character(len=BUFSIZE) :: outputPath
         probeBoundsExtension = get_coordinates_extension(this%mainCoords, this%auxCoords, control%mpidir)
         prefixFieldExtension = get_prefix_extension(field, control%mpidir)
         outputPath = &
            trim(adjustl(outputTypeExtension))//'_'//trim(adjustl(prefixFieldExtension))//'_'//trim(adjustl(probeBoundsExtension))
         return
      end function get_output_path

   end subroutine init_farField_probe_output

   subroutine update_farField_probe_output(this, ntime, bounds, fieldsReference)
    type(far_field_probe_output_t), intent(inout) :: this
    type(fields_reference_t), intent(in) :: fieldsReference
    integer(kind=SINGLE), intent(in) :: ntime
    type(bounds_t), intent(in) :: bounds
    call UpdateFarField(ntime, bounds, &
      fieldsReference%E%x, fieldsReference%E%y, fieldsReference%E%z, &
      fieldsReference%H%x, fieldsReference%H%y, fieldsReference%H%z)
   end subroutine update_farField_probe_output

   subroutine flush_farField_probe_output(this, simlulationTimeArray, timeIndex, control, fieldsReference, bounds)
    type(far_field_probe_output_t), intent(out) :: this
    integer, intent(in) :: timeIndex
    real(KIND=RKIND_tiempo), pointer, dimension(:), intent(in) :: simlulationTimeArray
    type(sim_control_t), intent(in) :: control
    type(fields_reference_t), pointer, intent(in) :: fieldsReference
    type(bounds_t), intent(in) :: bounds

    real(kind=RKIND_tiempo) :: flushTime

    flushTime = simlulationTimeArray(timeIndex)
    call FlushFarfield(control%layoutnumber, control%size, bounds, & 
    fieldsReference%E%deltaX, fieldsReference%E%deltaY, fieldsReference%E%deltaZ, &
    fieldsReference%H%deltaX, fieldsReference%H%deltaY, fieldsReference%H%deltaZ, &
    control%facesNF2FF, flushTime)
    end subroutine flush_farfield_probe_output

end module mod_farFieldOutput
