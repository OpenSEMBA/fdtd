module output_m
   use FDETYPES_m
   use report_m
   use domain_m
   use outputUtils_m
   use directoryUtils_m, only: delete_file, file_exists, get_last_component, join_path, remove_folder
   use pointProbeOutput_m
   use wireProbeOutput_m
   use bulkProbeOutput_m
   use lineProbeOutput_m
   use movieProbeOutput_m
   use frequencySliceProbeOutput_m
   use farFieldOutput_m
   use mapVTKOutput_m
   use outputDecomposition_m, only: output_partition_t, build_output_partition, &
                                    OUTPUT_PARTITION_SUCCESS, OUTPUT_PARTITION_INVALID_ARGUMENT
   use outputCollective_m, only: output_collective_t, init_output_collective, select_point_owner, &
                                 select_output_publication_mode, OUTPUT_COLLECTIVE_SUCCESS, &
                                 OUTPUT_COLLECTIVE_UNOWNED_POINT
   use outputTransport_m, only: output_transport_t, gather_point_eligibility, OUTPUT_TRANSPORT_SUCCESS
   use outputMetadata_m, only: publish_initial_probe_metadata, publish_final_probe_metadata, json_escape, &
                               json_unescape, OUTPUT_METADATA_SUCCESS
   use outputTypes_m, only: probe_metadata_t, output_artifact_t, output_lifecycle_is_terminal, &
                            probe_metadata_is_complete, OUTPUT_ARTIFACT_UNDEFINED, &
                            probe_publication_plan_t, &
                            OUTPUT_LIFECYCLE_DECLARED, OUTPUT_LIFECYCLE_ACTIVE, &
                            OUTPUT_LIFECYCLE_FINALISING, OUTPUT_LIFECYCLE_COMPLETE, OUTPUT_LIFECYCLE_FAILED
#ifdef CompileWithMPI
    use mpi
#endif
#ifdef CompileWithMTLN
   use Wire_bundles_mtln_m, only: GetSolverPtr
   use mtln_solver_m, only: mtln_solver_t => mtln_t
#endif

   implicit none
   private

   !===========================
   !  Public interface summary
   !===========================
   public :: solver_output_t
   public :: GetOutputs
   public :: init_outputs
   public :: update_outputs
   public :: flush_outputs
   public :: close_outputs
   public :: GetOutputPartition
   public :: probe_output_t, run_output_manifest_t
   public :: init_run_output_manifest, declare_probe_output, begin_probe_output, finalise_probe_output, &
             fail_probe_output, finalise_run_outputs, select_probe_participants
   public :: prepare_point_publication_plan, publication_plan_allows_canonical_write
   public :: prepare_distributed_point_publication_plan, publish_planned_probe_metadata, &
             finalise_transport_run_outputs
   public :: publish_scalar_probe_metadata
   public :: delete_run_output_manifest
   public :: OUTPUT_COORDINATION_SUCCESS, OUTPUT_COORDINATION_INVALID_PROBE, &
             OUTPUT_COORDINATION_INVALID_STATE, OUTPUT_COORDINATION_INVALID_ARTIFACTS, &
             OUTPUT_COORDINATION_NOT_TERMINAL, OUTPUT_COORDINATION_NOT_ROOT, &
             OUTPUT_COORDINATION_INVALID_OWNERSHIP, OUTPUT_COORDINATION_UNOWNED_POINT, &
             OUTPUT_COORDINATION_TRANSPORT_ERROR

   public :: POINT_PROBE_ID, WIRE_CURRENT_PROBE_ID, WIRE_CHARGE_PROBE_ID, BULK_PROBE_ID, VOLUMIC_CURRENT_PROBE_ID, &
             MOVIE_PROBE_ID, FREQUENCY_SLICE_PROBE_ID, FAR_FIELD_PROBE_ID, LINE_PROBE_ID
   !===========================

   !===========================
   !  Private interface summary
   !===========================
   private :: get_required_output_count
   !===========================

   integer(kind=SINGLE), parameter :: UNDEFINED_PROBE = -1, &
                                      POINT_PROBE_ID = 0, &
                                      WIRE_CURRENT_PROBE_ID = 1, &
                                      WIRE_CHARGE_PROBE_ID = 2, &
                                      BULK_PROBE_ID = 3, &
                                      VOLUMIC_CURRENT_PROBE_ID = 4, &
                                      MOVIE_PROBE_ID = 5, &
                                      FREQUENCY_SLICE_PROBE_ID = 6, &
                                      FAR_FIELD_PROBE_ID = 7, &
                                      MAPVTK_ID = 8, LINE_PROBE_ID = 9

   integer, parameter :: OUTPUT_COORDINATION_SUCCESS = 0
   integer, parameter :: OUTPUT_COORDINATION_INVALID_PROBE = 1
   integer, parameter :: OUTPUT_COORDINATION_INVALID_STATE = 2
   integer, parameter :: OUTPUT_COORDINATION_INVALID_ARTIFACTS = 3
   integer, parameter :: OUTPUT_COORDINATION_NOT_TERMINAL = 4
   integer, parameter :: OUTPUT_COORDINATION_NOT_ROOT = 5
   integer, parameter :: OUTPUT_COORDINATION_INVALID_OWNERSHIP = 6
   integer, parameter :: OUTPUT_COORDINATION_UNOWNED_POINT = 7
   integer, parameter :: OUTPUT_COORDINATION_TRANSPORT_ERROR = 8

   type :: probe_output_t
      type(probe_metadata_t) :: metadata
      type(probe_publication_plan_t) :: publication_plan
   end type probe_output_t

   type :: run_output_manifest_t
      character(len=BUFSIZE) :: run_id = ''
      integer :: root_rank = 0
      logical :: published = .false.
      type(probe_output_t), allocatable :: probes(:)
   end type run_output_manifest_t

   REAL(KIND=RKIND), save           ::  eps0, mu0
   REAL(KIND=RKIND), pointer, dimension(:), save  ::  InvEps, InvMu
   type(solver_output_t), pointer, dimension(:), save  ::  outputs
   type(output_partition_t), allocatable, save :: outputPartitions(:)
   type(problem_info_t), save, target :: problemInfo

   interface init_solver_output
      module procedure &
         init_point_probe_output, &
         init_wire_current_probe_output, &
         init_wire_charge_probe_output, &
         init_bulk_probe_output, &
         init_line_probe_output, &
         init_movie_probe_output, &
         init_frequency_slice_probe_output, &
         init_farField_probe_output, &
         init_mapvtk_output
   end interface

   interface update_solver_output
      module procedure &
         update_point_probe_output, &
         update_wire_current_probe_output, &
         update_wire_charge_probe_output, &
         update_bulk_probe_output, &
         update_line_probe_output, &
         update_movie_probe_output, &
         update_frequency_slice_probe_output, &
         update_farField_probe_output
   end interface

   interface flush_solver_output
      module procedure &
         flush_point_probe_output, &
         flush_wire_current_probe_output, &
         flush_wire_charge_probe_output, &
         flush_bulk_probe_output, &
         flush_line_probe_output, &
         flush_movie_probe_output, &
         flush_frequency_slice_probe_output, &
         flush_farField_probe_output

   end interface

   interface close_solver_output
      module procedure close_movie_probe_output
      module procedure close_frequency_slice_probe_output
   end interface
contains

   function GetOutputs() result(r)
      type(solver_output_t), pointer, dimension(:)  ::  r
      r => outputs
      return
   end function

   function GetProblemInfo() result(r)
      type(problem_info_t), pointer ::  r
      r => problemInfo
      return
   end function

   subroutine prepare_point_publication_plan(plan, rank, rank_count, rank_is_eligible, status)
      type(probe_publication_plan_t), intent(out) :: plan
      integer, intent(in) :: rank, rank_count
      logical, intent(in) :: rank_is_eligible(:)
      integer, intent(out) :: status
      type(output_collective_t) :: collective

      plan = probe_publication_plan_t()
      call init_output_collective(collective, rank, rank_count, 0, .false., status)
      if (status /= OUTPUT_COLLECTIVE_SUCCESS) return

      call select_point_owner(collective, rank_is_eligible, plan%canonical_writer_rank, status)
      if (status /= OUTPUT_COLLECTIVE_SUCCESS) return

      plan%local_eligible = rank_is_eligible(rank + 1)
      plan%local_participates = plan%local_eligible
      plan%local_is_canonical_writer = plan%canonical_writer_rank == rank
   end subroutine prepare_point_publication_plan

   pure logical function publication_plan_allows_canonical_write(plan)
      type(probe_publication_plan_t), intent(in) :: plan

      publication_plan_allows_canonical_write = plan%canonical_writer_rank >= 0 .and. &
                                                plan%local_is_canonical_writer
   end function publication_plan_allows_canonical_write

   subroutine prepare_distributed_point_publication_plan(plan, transport, local_eligible, status)
      type(probe_publication_plan_t), intent(out) :: plan
      type(output_transport_t), intent(in) :: transport
      logical, intent(in) :: local_eligible
      logical, allocatable :: rank_is_eligible(:)
      integer, intent(out) :: status
      integer :: collective_status, transport_status

      plan = probe_publication_plan_t()
      call gather_point_eligibility(transport, local_eligible, rank_is_eligible, transport_status)
      if (transport_status /= OUTPUT_TRANSPORT_SUCCESS) then
         status = OUTPUT_COORDINATION_TRANSPORT_ERROR
         return
      end if

      call prepare_point_publication_plan(plan, transport%rank, transport%rank_count, rank_is_eligible, &
                                          collective_status)
      if (collective_status == OUTPUT_COLLECTIVE_SUCCESS) then
         status = OUTPUT_COORDINATION_SUCCESS
      else if (collective_status == OUTPUT_COLLECTIVE_UNOWNED_POINT) then
         status = OUTPUT_COORDINATION_UNOWNED_POINT
      else
         status = OUTPUT_COORDINATION_INVALID_OWNERSHIP
      end if
   end subroutine prepare_distributed_point_publication_plan

   subroutine publish_planned_probe_metadata(path, metadata, plan, status)
      character(len=*), intent(in) :: path
      type(probe_metadata_t), intent(in) :: metadata
      type(probe_publication_plan_t), intent(in) :: plan
      integer, intent(out) :: status
      integer :: metadata_status

      if (.not. publication_plan_allows_canonical_write(plan)) then
         status = OUTPUT_COORDINATION_SUCCESS
         return
      end if
      if (output_lifecycle_is_terminal(metadata%lifecycle)) then
         call publish_final_probe_metadata(path, metadata, metadata_status)
      else
         call publish_initial_probe_metadata(path, metadata, metadata_status)
      end if
      if (metadata_status == OUTPUT_METADATA_SUCCESS) then
         status = OUTPUT_COORDINATION_SUCCESS
      else
         status = OUTPUT_COORDINATION_INVALID_ARTIFACTS
      end if
   end subroutine publish_planned_probe_metadata

   subroutine publish_scalar_probe_metadata(path, probe_id, quantity, artifacts, status)
      character(len=*), intent(in) :: path, probe_id, quantity
      type(output_artifact_t), intent(in) :: artifacts(:)
      integer, intent(out) :: status
      type(probe_metadata_t) :: metadata
      type(output_artifact_t), allocatable :: normalised_artifacts(:)
      integer :: i, artifact_count, metadata_status

      status = OUTPUT_COORDINATION_INVALID_ARTIFACTS
      if (len_trim(probe_id) == 0 .or. len_trim(quantity) == 0) return
      artifact_count = count(artifacts%kind /= OUTPUT_ARTIFACT_UNDEFINED)
      if (artifact_count == 0) return

      metadata%probe_id = probe_id
      metadata%quantity = quantity
      allocate (metadata%artifacts(artifact_count))
      allocate (normalised_artifacts(artifact_count))
      artifact_count = 0
      do i = 1, size(artifacts)
         if (artifacts(i)%kind == OUTPUT_ARTIFACT_UNDEFINED) cycle
         artifact_count = artifact_count + 1
         normalised_artifacts(artifact_count) = artifacts(i)
         normalised_artifacts(artifact_count)%relative_path = &
            get_last_component(artifacts(i)%relative_path)
      end do
      metadata%artifacts = normalised_artifacts
      call publish_initial_probe_metadata(path, metadata, metadata_status)
      if (metadata_status == OUTPUT_METADATA_SUCCESS) then
         status = OUTPUT_COORDINATION_SUCCESS
      end if
   end subroutine publish_scalar_probe_metadata

   subroutine finalise_scalar_output_metadata(output_index)
      integer(kind=SINGLE), intent(in) :: output_index
      integer :: i, path_end, status
      logical :: complete
      character(len=:), allocatable :: artifact_path

      if (.not. allocated(outputs(output_index)%metadata%artifacts)) return
      complete = .true.
      do i = 1, size(outputs(output_index)%metadata%artifacts)
         path_end = scan(trim(outputs(output_index)%metadata_path), '/\', back=.true.)
         if (path_end > 0) then
            artifact_path = trim(outputs(output_index)%metadata_path(:path_end))// &
                            trim(outputs(output_index)%metadata%artifacts(i)%relative_path)
         else
            artifact_path = trim(outputs(output_index)%metadata%artifacts(i)%relative_path)
         end if
         if (outputs(output_index)%metadata%artifacts(i)%required .and. &
             .not. file_exists(artifact_path)) complete = .false.
      end do
      if (complete) then
         outputs(output_index)%metadata%lifecycle%state = OUTPUT_LIFECYCLE_COMPLETE
         outputs(output_index)%metadata%lifecycle%diagnostic = ''
      else
         outputs(output_index)%metadata%lifecycle%state = OUTPUT_LIFECYCLE_FAILED
         outputs(output_index)%metadata%lifecycle%diagnostic = 'Required scalar artifacts are incomplete'
      end if
      call publish_final_probe_metadata(outputs(output_index)%metadata_path, outputs(output_index)%metadata, status)
   end subroutine finalise_scalar_output_metadata

   subroutine GetOutputPartition(output_index, partition, status)
      integer, intent(in) :: output_index
      type(output_partition_t), intent(out) :: partition
      integer, intent(out) :: status

      partition = output_partition_t()
      if (.not. allocated(outputPartitions) .or. output_index < 1 .or. output_index > size(outputPartitions)) then
         status = OUTPUT_PARTITION_INVALID_ARGUMENT
         return
      end if

      partition = outputPartitions(output_index)
      status = OUTPUT_PARTITION_SUCCESS
   end subroutine GetOutputPartition

   subroutine init_run_output_manifest(manifest, run_id, root_rank)
      type(run_output_manifest_t), intent(out) :: manifest
      character(len=*), intent(in) :: run_id
      integer, intent(in) :: root_rank

      manifest%run_id = run_id
      manifest%root_rank = root_rank
   end subroutine init_run_output_manifest

   subroutine declare_probe_output(manifest, probe_id, quantity, artifacts, probe_index, status)
      type(run_output_manifest_t), intent(inout) :: manifest
      character(len=*), intent(in) :: probe_id
      character(len=*), intent(in) :: quantity
      type(output_artifact_t), intent(in) :: artifacts(:)
      integer, intent(out) :: probe_index, status
      type(probe_output_t), allocatable :: expanded_probes(:)
      integer :: i, probe_count

      probe_index = 0
      if (len_trim(probe_id) == 0 .or. len_trim(quantity) == 0 .or. &
          .not. artifacts_are_declared(artifacts)) then
         status = OUTPUT_COORDINATION_INVALID_ARTIFACTS
         return
      end if

      if (allocated(manifest%probes)) then
         do i = 1, size(manifest%probes)
            if (trim(manifest%probes(i)%metadata%probe_id) == trim(probe_id)) then
               status = OUTPUT_COORDINATION_INVALID_PROBE
               return
            end if
         end do
         probe_count = size(manifest%probes)
      else
         probe_count = 0
      end if

      allocate (expanded_probes(probe_count + 1))
      if (probe_count > 0) expanded_probes(:probe_count) = manifest%probes
      expanded_probes(probe_count + 1)%metadata%probe_id = probe_id
      expanded_probes(probe_count + 1)%metadata%quantity = quantity
      allocate (expanded_probes(probe_count + 1)%metadata%artifacts(size(artifacts)))
      expanded_probes(probe_count + 1)%metadata%artifacts = artifacts
      call move_alloc(expanded_probes, manifest%probes)
      probe_index = probe_count + 1
      status = OUTPUT_COORDINATION_SUCCESS
   end subroutine declare_probe_output

   subroutine begin_probe_output(manifest, probe_index, status)
      type(run_output_manifest_t), intent(inout) :: manifest
      integer, intent(in) :: probe_index
      integer, intent(out) :: status

      if (.not. valid_probe_index(manifest, probe_index)) then
         status = OUTPUT_COORDINATION_INVALID_PROBE
      else if (manifest%probes(probe_index)%metadata%lifecycle%state /= OUTPUT_LIFECYCLE_DECLARED) then
         status = OUTPUT_COORDINATION_INVALID_STATE
      else
         manifest%probes(probe_index)%metadata%lifecycle%state = OUTPUT_LIFECYCLE_ACTIVE
         status = OUTPUT_COORDINATION_SUCCESS
      end if
   end subroutine begin_probe_output

   subroutine select_probe_participants(manifest, probe_index, participant_ranks, scalar_writer_rank, status)
      type(run_output_manifest_t), intent(inout) :: manifest
      integer, intent(in) :: probe_index
      integer, intent(in) :: participant_ranks(:)
      integer, intent(in) :: scalar_writer_rank
      integer, intent(out) :: status

      if (.not. valid_probe_index(manifest, probe_index)) then
         status = OUTPUT_COORDINATION_INVALID_PROBE
         return
      end if
      if (manifest%probes(probe_index)%metadata%lifecycle%state /= OUTPUT_LIFECYCLE_DECLARED) then
         status = OUTPUT_COORDINATION_INVALID_STATE
         return
      end if
      if (.not. valid_probe_ownership(participant_ranks, scalar_writer_rank)) then
         status = OUTPUT_COORDINATION_INVALID_OWNERSHIP
         return
      end if

      associate (ownership => manifest%probes(probe_index)%metadata%ownership)
         ownership%participant_ranks = participant_ranks
         ownership%scalar_writer_rank = scalar_writer_rank
      end associate
      status = OUTPUT_COORDINATION_SUCCESS
   end subroutine select_probe_participants

   subroutine finalise_probe_output(manifest, probe_index, status)
      type(run_output_manifest_t), intent(inout) :: manifest
      integer, intent(in) :: probe_index
      integer, intent(out) :: status

      if (.not. valid_probe_index(manifest, probe_index)) then
         status = OUTPUT_COORDINATION_INVALID_PROBE
         return
      end if

      associate (metadata => manifest%probes(probe_index)%metadata)
         if (metadata%lifecycle%state /= OUTPUT_LIFECYCLE_DECLARED .and. &
             metadata%lifecycle%state /= OUTPUT_LIFECYCLE_ACTIVE) then
            status = OUTPUT_COORDINATION_INVALID_STATE
            return
         end if

         metadata%lifecycle%state = OUTPUT_LIFECYCLE_FINALISING
         metadata%lifecycle%state = OUTPUT_LIFECYCLE_COMPLETE
         if (probe_metadata_is_complete(metadata)) then
            status = OUTPUT_COORDINATION_SUCCESS
         else
            metadata%lifecycle%state = OUTPUT_LIFECYCLE_FAILED
            metadata%lifecycle%diagnostic = 'Required artifacts are incomplete'
            status = OUTPUT_COORDINATION_INVALID_ARTIFACTS
         end if
      end associate
   end subroutine finalise_probe_output

   subroutine fail_probe_output(manifest, probe_index, diagnostic, status)
      type(run_output_manifest_t), intent(inout) :: manifest
      integer, intent(in) :: probe_index
      character(len=*), intent(in) :: diagnostic
      integer, intent(out) :: status

      if (.not. valid_probe_index(manifest, probe_index)) then
         status = OUTPUT_COORDINATION_INVALID_PROBE
         return
      end if
      if (output_lifecycle_is_terminal(manifest%probes(probe_index)%metadata%lifecycle)) then
         status = OUTPUT_COORDINATION_INVALID_STATE
         return
      end if

      manifest%probes(probe_index)%metadata%lifecycle%state = OUTPUT_LIFECYCLE_FAILED
      manifest%probes(probe_index)%metadata%lifecycle%diagnostic = diagnostic
      status = OUTPUT_COORDINATION_SUCCESS
   end subroutine fail_probe_output

   subroutine finalise_run_outputs(manifest, writer_rank, status)
      type(run_output_manifest_t), intent(inout) :: manifest
      integer, intent(in) :: writer_rank
      integer, intent(out) :: status
      integer :: i

      if (writer_rank /= manifest%root_rank) then
         status = OUTPUT_COORDINATION_NOT_ROOT
         return
      end if
      if (allocated(manifest%probes)) then
         do i = 1, size(manifest%probes)
            if (.not. output_lifecycle_is_terminal(manifest%probes(i)%metadata%lifecycle)) then
               status = OUTPUT_COORDINATION_NOT_TERMINAL
               return
            end if
         end do
      end if

      manifest%published = .true.
      status = OUTPUT_COORDINATION_SUCCESS
   end subroutine finalise_run_outputs

   subroutine finalise_transport_run_outputs(manifest, transport, status)
      type(run_output_manifest_t), intent(inout) :: manifest
      type(output_transport_t), intent(in) :: transport
      integer, intent(out) :: status
      integer :: i

      if (transport%root_rank /= manifest%root_rank) then
         status = OUTPUT_COORDINATION_INVALID_OWNERSHIP
         return
      end if
      if (transport%rank == transport%root_rank) then
         call finalise_run_outputs(manifest, transport%rank, status)
         return
      end if
      if (allocated(manifest%probes)) then
         do i = 1, size(manifest%probes)
            if (.not. output_lifecycle_is_terminal(manifest%probes(i)%metadata%lifecycle)) then
               status = OUTPUT_COORDINATION_NOT_TERMINAL
               return
            end if
         end do
      end if
      status = OUTPUT_COORDINATION_SUCCESS
   end subroutine finalise_transport_run_outputs

   pure logical function valid_probe_index(manifest, probe_index)
      type(run_output_manifest_t), intent(in) :: manifest
      integer, intent(in) :: probe_index

      valid_probe_index = allocated(manifest%probes) .and. probe_index >= 1
      if (valid_probe_index) valid_probe_index = probe_index <= size(manifest%probes)
   end function valid_probe_index

   pure logical function artifacts_are_declared(artifacts)
      type(output_artifact_t), intent(in) :: artifacts(:)
      integer :: i

      artifacts_are_declared = size(artifacts) > 0
      if (.not. artifacts_are_declared) return
      do i = 1, size(artifacts)
         if (artifacts(i)%kind == OUTPUT_ARTIFACT_UNDEFINED .or. len_trim(artifacts(i)%relative_path) == 0) then
            artifacts_are_declared = .false.
            return
         end if
      end do
   end function artifacts_are_declared

   pure logical function valid_probe_ownership(participant_ranks, scalar_writer_rank)
      integer, intent(in) :: participant_ranks(:)
      integer, intent(in) :: scalar_writer_rank
      integer :: i

      valid_probe_ownership = size(participant_ranks) > 0 .and. any(participant_ranks == scalar_writer_rank)
      if (.not. valid_probe_ownership) return
      do i = 1, size(participant_ranks) - 1
         if (any(participant_ranks(i + 1:) == participant_ranks(i))) then
            valid_probe_ownership = .false.
            return
         end if
      end do
   end function valid_probe_ownership

   subroutine init_outputs(sgg, media, sinpml_fullsize, materialTags, bounds, control, observationsExists, &
                           wiresExists, eps0_input, mu0_input)

      type(SGGFDTDINFO_t), intent(in) ::  sgg
      type(media_matrices_t), target, intent(in) :: media
      type(limit_t), dimension(:), target, intent(in)  ::  SINPML_fullsize
      type(bounds_t), intent(in), target :: bounds
      type(taglist_t), intent(in), target :: materialTags
      type(sim_control_t), intent(in) :: control
      logical, intent(inout) :: wiresExists
      logical, intent(out) :: observationsExists
      real(kind=RKIND), intent(in), optional :: eps0_input, mu0_input

      type(domain_t) :: domain
      type(spheric_domain_t) :: sphericRange
      type(cell_coordinate_t) :: lowerBound, upperBound, localMapLower, localMapUpper
      integer(kind=SINGLE) :: i, ii, outputRequestType
      integer(kind=SINGLE) :: NODE
      integer(kind=SINGLE) :: outputCount
      integer(kind=SINGLE) :: requestedOutputs
      integer :: metadata_status
      logical :: mapHasData
      character(len=BUFSIZE) :: outputTypeExtension

#ifdef CompileWithMTLN
      logical :: thereAreMtlnObservations = .false.
#endif
      observationsExists = .false.
      eps0 = EPSILON_VACUUM
      mu0 = MU_VACUUM
      if (present(eps0_input)) eps0 = eps0_input
      if (present(mu0_input)) mu0 = mu0_input
      requestedOutputs = get_required_output_count(sgg)

      problemInfo%geometryToMaterialData => media
      problemInfo%materialList => sgg%Med
      problemInfo%simulationBounds => bounds
      problemInfo%problemDimension => SINPML_fullsize
      problemInfo%materialTag => materialTags
      problemInfo%xSteps => sgg%LineX
      problemInfo%ySteps => sgg%LineY
      problemInfo%zSteps => sgg%LineZ

      outputs => NULL()
      allocate (outputs(requestedOutputs))
      if (allocated(outputPartitions)) deallocate (outputPartitions)
      allocate (outputPartitions(requestedOutputs))

      allocate (InvEps(lbound(sgg%Med, 1):ubound(sgg%Med, 1)), &
                InvMu(lbound(sgg%Med, 1):ubound(sgg%Med, 1)))
      outputCount = 0

      InvEps = 1.0_RKIND/(eps0*sgg%Med%Epr)
      InvMu = 1.0_RKIND/(mu0*sgg%Med%Mur)

      !do ii = 1, sgg%NumberRequest
      !do i = 1, sgg%Observation(ii)%nP
      !   call eliminate_unnecesary_observation_points(sgg%Observation(ii)%P(i), output(ii)%item(i), &
      !     sgg%Sweep, sgg%SINPMLSweep, sgg%Observation(ii)%P(1)%ZI, sgg%Observation(ii)%P(1)%ZE, control%layoutnumber, control%num_procs)
      !end do
      !end do

#ifdef CompileWithMTLN
      block
         type(mtln_solver_t), pointer :: mtln_solver
         integer :: i, j
         mtln_solver => GetSolverPtr()
         if (.not. allocated(mtln_solver%bundles)) then
            if (control%layoutnumber == 0) then
               call print11(0, 'WARNING: mtln solver has not bundles allocated. No MTLN bundle requests were found.')
            end if
         else if (size(mtln_solver%bundles) == 0) then
            if (control%layoutnumber == 0) then
               call print11(0, 'WARNING: mtln solver has not bundles registered. No MTLN bundle requests were found.')
            end if
         else
            mtlnObservationExist: do i = 1, size(mtln_solver%bundles)
               if (.not. allocated(mtln_solver%bundles(i)%probes)) cycle
               do j = 1, size(mtln_solver%bundles(i)%probes)
                  if (mtln_solver%bundles(i)%probes(j)%in_layer) then
                     thereAreMtlnObservations = .true.
                     exit mtlnObservationExist
                  end if
               end do
            end do mtlnObservationExist
         end if
      end block
#endif

      do ii = 1, sgg%NumberRequest
         domain = preprocess_domain(sgg%Observation(ii), sgg%tiempo, sgg%dt, control%finaltimestep)
         if (domain%domainType == UNDEFINED_DOMAIN) cycle
         do i = 1, sgg%Observation(ii)%nP
            lowerBound%x = sgg%observation(ii)%P(i)%XI
            lowerBound%y = sgg%observation(ii)%P(i)%YI
            lowerBound%z = sgg%observation(ii)%P(i)%ZI

            upperBound%x = sgg%observation(ii)%P(i)%XE
            upperBound%y = sgg%observation(ii)%P(i)%YE
            upperBound%z = sgg%observation(ii)%P(i)%ZE
            NODE = sgg%observation(ii)%P(i)%NODE

            outputTypeExtension = trim(adjustl(control%nEntradaRoot))//'_'//trim(adjustl(sgg%observation(ii)%outputrequest))

            outputRequestType = sgg%observation(ii)%P(i)%what
            select case (outputRequestType)
            case (mapvtk)
                call get_local_map_bounds(lowerBound, upperBound, localMapLower, localMapUpper, mapHasData)
                if (mapHasData) then
                   outputCount = outputCount + 1
                   outputs(outputCount)%outputID = MAPVTK_ID

                   allocate (outputs(outputCount)%mapvtkOutput)
                   call init_solver_output(outputs(outputCount)%mapvtkOutput, localMapLower, localMapUpper, &
                                            outputRequestType, outputTypeExtension, control%mpidir, problemInfo)
                   call create_geometry_simulation_vtu(outputs(outputCount)%mapvtkOutput, control, sgg%LineX, sgg%LineY, &
                                                       sgg%LineZ, problemInfo)
                   call register_scalar_output_metadata(outputCount, &
                                                        join_path(outputs(outputCount)%mapvtkOutput%path, &
                                                                 get_last_component(outputs(outputCount)%mapvtkOutput%path)//'.json'), &
                                                        get_last_component(outputs(outputCount)%mapvtkOutput%path), &
                                                        get_prefix_extension(outputRequestType, control%mpidir), &
                                                        outputs(outputCount)%mapvtkOutput%artifacts, metadata_status)
                end if

            case (iEx, iEy, iEz, iHx, iHy, iHz)
               outputCount = outputCount + 1
               outputs(outputCount)%outputID = POINT_PROBE_ID

               allocate (outputs(outputCount)%pointProbe)
                call init_solver_output(outputs(outputCount)%pointProbe, lowerBound, outputRequestType, domain, outputTypeExtension, control%mpidir, sgg%dt, &
                                       sgg%NumPlaneWaves >= 1)
               call register_scalar_output_metadata(outputCount, trim(outputs(outputCount)%pointProbe%path)//'.json', &
                                                    get_last_component(outputs(outputCount)%pointProbe%path), &
                                                    get_prefix_extension(outputRequestType, control%mpidir), &
                                                    outputs(outputCount)%pointProbe%artifacts, metadata_status)
            case (iJx, iJy, iJz)
               if (wiresExists) then
                  outputCount = outputCount + 1
                  outputs(outputCount)%outputID = WIRE_CURRENT_PROBE_ID

                  allocate (outputs(outputCount)%wireCurrentProbe)
                  call init_solver_output(outputs(outputCount)%wireCurrentProbe, lowerBound, NODE, outputRequestType, domain, problemInfo%materialList, outputTypeExtension, control%mpidir, control%wiresflavor)
                  call register_scalar_output_metadata(outputCount, trim(outputs(outputCount)%wireCurrentProbe%path)//'.json', &
                                                       get_last_component(outputs(outputCount)%wireCurrentProbe%path), &
                                                       get_prefix_extension(outputRequestType, control%mpidir), &
                                                       outputs(outputCount)%wireCurrentProbe%artifacts, metadata_status)
               end if

            case (iQx, iQy, iQz)
               outputCount = outputCount + 1
               outputs(outputCount)%outputID = WIRE_CHARGE_PROBE_ID

               allocate (outputs(outputCount)%wireChargeProbe)
                call init_solver_output(outputs(outputCount)%wireChargeProbe, lowerBound, NODE, outputRequestType, domain, outputTypeExtension, control%mpidir, control%wiresflavor)
               call register_scalar_output_metadata(outputCount, trim(outputs(outputCount)%wireChargeProbe%path)//'.json', &
                                                    get_last_component(outputs(outputCount)%wireChargeProbe%path), &
                                                    get_prefix_extension(outputRequestType, control%mpidir), &
                                                    outputs(outputCount)%wireChargeProbe%artifacts, metadata_status)

            case (iBloqueJx, iBloqueJy, iBloqueJz, iBloqueMx, iBloqueMy, iBloqueMz)
               outputCount = outputCount + 1
               outputs(outputCount)%outputID = BULK_PROBE_ID

               allocate (outputs(outputCount)%bulkCurrentProbe)
                call init_solver_output(outputs(outputCount)%bulkCurrentProbe, lowerBound, upperBound, outputRequestType, domain, outputTypeExtension, control%mpidir)
               call register_scalar_output_metadata(outputCount, trim(outputs(outputCount)%bulkCurrentProbe%path)//'.json', &
                                                    get_last_component(outputs(outputCount)%bulkCurrentProbe%path), &
                                                    get_prefix_extension(outputRequestType, control%mpidir), &
                                                    outputs(outputCount)%bulkCurrentProbe%artifacts, metadata_status)
                !! call adjust_computation_range --- Required due to issues in mpi region edges

            case (lineIntegral)
               if (domain%domainType /= TIME_DOMAIN) then
                  call stoponerror(0, 0, 'Line probes only support the time domain')
               else
                  outputCount = outputCount + 1
                  outputs(outputCount)%outputID = LINE_PROBE_ID
                  allocate (outputs(outputCount)%lineProbe)
                  call init_line_probe_output(outputs(outputCount)%lineProbe, sgg%observation(ii)%P(i)%line, domain, &
                                              join_path(trim(outputTypeExtension)//'_LI', trim(outputTypeExtension)//'_LI'), &
                                              sgg%Sweep, control%layoutnumber, control%num_procs)
                  call register_scalar_output_metadata(outputCount, trim(outputs(outputCount)%lineProbe%path)//'.json', &
                                                       get_last_component(outputs(outputCount)%lineProbe%path), 'LI', &
                                                       outputs(outputCount)%lineProbe%artifacts, metadata_status)
                  outputs(outputCount)%metadata%domain_type = domain%domainType
                  if (allocated(outputs(outputCount)%metadata%ownership%participant_ranks)) then
                     deallocate (outputs(outputCount)%metadata%ownership%participant_ranks)
                  end if
                  allocate (outputs(outputCount)%metadata%ownership%participant_ranks(1))
                  outputs(outputCount)%metadata%ownership%participant_ranks(1) = control%layoutnumber
                  outputs(outputCount)%metadata%ownership%scalar_writer_rank = control%layoutnumber
                  call publish_initial_probe_metadata(outputs(outputCount)%metadata_path, outputs(outputCount)%metadata, &
                                                      metadata_status)
               end if

            case (iCur, iMEC, iMHC, iCurX, iCurY, iCurZ, iExC, iEyC, iEzC, iHxC, iHyC, iHzC)
               if (domain%domainType == TIME_DOMAIN) then

                  outputCount = outputCount + 1
                  outputs(outputCount)%outputID = MOVIE_PROBE_ID
                  call attach_output_partition(outputCount)
                  if (outputPartitions(outputCount)%has_data) then
                     allocate (outputs(outputCount)%movieProbe)
                     call init_solver_output(outputs(outputCount)%movieProbe, outputPartitions(outputCount)%local_lower, &
                                              outputPartitions(outputCount)%local_upper, outputRequestType, domain, control, &
                                              problemInfo, outputTypeExtension)
                     call configure_output_partition(outputCount)
                  else
                     outputs(outputCount)%outputID = UNDEFINED_PROBE
                     outputCount = outputCount - 1
                  end if
               else if (domain%domainType == FREQUENCY_DOMAIN) then

                  outputCount = outputCount + 1
                  outputs(outputCount)%outputID = FREQUENCY_SLICE_PROBE_ID
                  call attach_output_partition(outputCount)
                  if (outputPartitions(outputCount)%has_data) then
                     allocate (outputs(outputCount)%frequencySliceProbe)
                     call init_solver_output(outputs(outputCount)%frequencySliceProbe, outputPartitions(outputCount)%local_lower, &
                                              outputPartitions(outputCount)%local_upper, sgg%dt, outputRequestType, domain, &
                                              outputTypeExtension, control, problemInfo)
                     call configure_output_partition(outputCount)
                  else
                     outputs(outputCount)%outputID = UNDEFINED_PROBE
                     outputCount = outputCount - 1
                  end if
               end if
            case (farfield)
               sphericRange = preprocess_polar_range(sgg%Observation(ii))

               outputCount = outputCount + 1
               outputs(outputCount)%outputID = FAR_FIELD_PROBE_ID
               allocate (outputs(outputCount)%farFieldOutput)
#ifdef CompileWithMPI
                call configure_far_field_mpi(outputs(outputCount), lowerBound, upperBound)
#endif
                call init_solver_output(outputs(outputCount)%farFieldOutput, sgg, lowerBound, upperBound, outputRequestType, domain, sphericRange, outputTypeExtension, sgg%Observation(ii)%FileNormalize, control, problemInfo, &
#ifdef CompileWithMPI
                                         outputs(outputCount)%MPISubcomm, outputs(outputCount)%MPIRoot, &
#endif
                                         eps0, mu0)
               call register_scalar_output_metadata(outputCount, trim(outputs(outputCount)%farFieldOutput%path)//'.json', &
                                                    get_last_component(outputs(outputCount)%farFieldOutput%path), &
                                                    get_prefix_extension(outputRequestType, control%mpidir), &
                                                    outputs(outputCount)%farFieldOutput%artifacts, metadata_status)
            case default
               call stoponerror(0, 0, 'OutputRequestType type not implemented yet on new observations')
            end select
         end do
      end do
      if (outputCount /= requestedOutputs) then
         call remove_unused_outputs(outputs)
         outputCount = size(outputs)
      end if
      if (outputCount /= 0) observationsExists = .true.
#ifdef CompileWithMTLN
      observationsExists = observationsExists .or. thereAreMtlnObservations
#endif
      if (observationsExists) call write_run_output_manifest(control, outputCount)
      return
   contains
#ifdef CompileWithMPI
       subroutine configure_far_field_mpi(output, lower_bound, upper_bound)
          type(solver_output_t), intent(inout) :: output
          type(cell_coordinate_t), intent(in) :: lower_bound, upper_bound
          integer :: color, ierr, root_candidate

          output%MPISubcomm = -1
          if (lower_bound%z <= sgg%SINPMLSweep(iHz)%ZE .and. &
              upper_bound%z >= sgg%SINPMLSweep(iHz)%ZI) then
             output%MPISubcomm = 1
          end if
          root_candidate = -1
          if (lower_bound%z >= sgg%SINPMLSweep(iHz)%ZI .and. &
              lower_bound%z < sgg%SINPMLSweep(iHz)%ZE) then
             root_candidate = control%layoutnumber
          end if
          call MPI_AllReduce(root_candidate, output%MPIRoot, 1_4, MPI_INTEGER, MPI_MAX, SUBCOMM_MPI, ierr)

          color = MPI_UNDEFINED
          if (output%MPISubcomm == 1) color = 0
          call MPI_Comm_split(SUBCOMM_MPI, color, control%layoutnumber, output%MPISubcomm, ierr)
          if (output%MPISubcomm /= MPI_COMM_NULL) then
             call MPI_Comm_rank(output%MPISubcomm, output%MPIGroupIndex, ierr)
          end if
       end subroutine configure_far_field_mpi
#endif

       subroutine register_scalar_output_metadata(output_index, descriptor_path, probe_id, quantity, artifacts, status)
         integer(kind=SINGLE), intent(in) :: output_index
         character(len=*), intent(in) :: descriptor_path, probe_id, quantity
         type(output_artifact_t), intent(in) :: artifacts(:)
         integer, intent(out) :: status
         integer :: i, artifact_count

         call publish_scalar_probe_metadata(descriptor_path, probe_id, quantity, artifacts, status)
         if (status /= OUTPUT_COORDINATION_SUCCESS) return
         outputs(output_index)%metadata%probe_id = probe_id
         outputs(output_index)%metadata%quantity = quantity
         allocate (outputs(output_index)%metadata%artifacts(count(artifacts%kind /= OUTPUT_ARTIFACT_UNDEFINED)))
         artifact_count = 0
         do i = 1, size(artifacts)
            if (artifacts(i)%kind == OUTPUT_ARTIFACT_UNDEFINED) cycle
            artifact_count = artifact_count + 1
            outputs(output_index)%metadata%artifacts(artifact_count) = artifacts(i)
            outputs(output_index)%metadata%artifacts(artifact_count)%relative_path = &
               get_last_component(artifacts(i)%relative_path)
         end do
         outputs(output_index)%metadata_path = descriptor_path
      end subroutine register_scalar_output_metadata

       subroutine attach_output_partition(output_index)
          integer(kind=SINGLE), intent(in) :: output_index
          type(limit_t) :: local_sweep
          integer :: field_component, partition_status

         field_component = fieldo(outputRequestType, 'Z')
         local_sweep%XI = sgg%Sweep(field_component)%XI
         local_sweep%XE = sgg%Sweep(field_component)%XE
         local_sweep%YI = sgg%Sweep(field_component)%YI
         local_sweep%YE = sgg%Sweep(field_component)%YE
         local_sweep%ZI = sgg%Sweep(field_component)%ZI
         local_sweep%ZE = sgg%Sweep(field_component)%ZE
         local_sweep%NX = local_sweep%XE - local_sweep%XI + 1
         local_sweep%NY = local_sweep%YE - local_sweep%YI + 1
         local_sweep%NZ = local_sweep%ZE - local_sweep%ZI + 1
         call build_output_partition(lowerBound, upperBound, SINPML_fullsize(field_component), local_sweep, &
                                     field_component, control%layoutnumber, max(control%num_procs, 1), &
                                     outputPartitions(output_index), partition_status)
         if (partition_status /= OUTPUT_PARTITION_SUCCESS) return

       end subroutine attach_output_partition

       subroutine configure_output_partition(output_index)
          integer(kind=SINGLE), intent(in) :: output_index
          type(output_collective_t) :: collective
          integer :: collective_status, publication_mode

          ! Parallel visualisation output is available when MPI spans multiple ranks.
          call init_output_collective(collective, control%layoutnumber, max(control%num_procs, 1), 0, &
                                      control%num_procs > 1, &
                                      collective_status)
         if (collective_status /= OUTPUT_COLLECTIVE_SUCCESS) return
         call select_output_publication_mode(collective, publication_mode)
         select case (outputs(output_index)%outputID)
         case (MOVIE_PROBE_ID)
            call configure_movie_probe_publication(outputs(output_index)%movieProbe, publication_mode, &
                                                   outputPartitions(output_index)%has_data)
          case (FREQUENCY_SLICE_PROBE_ID)
             call configure_frequency_slice_probe_publication(outputs(output_index)%frequencySliceProbe, &
                                                              publication_mode, outputPartitions(output_index)%has_data)
          end select
       end subroutine configure_output_partition

       subroutine get_local_map_bounds(request_lower, request_upper, local_lower, local_upper, has_data)
          type(cell_coordinate_t), intent(in) :: request_lower, request_upper
          type(cell_coordinate_t), intent(out) :: local_lower, local_upper
          logical, intent(out) :: has_data
          integer :: field

          local_lower = request_lower
          local_upper = request_upper
          do field = iEx, iHz
             ! Restrict iteration to owned cells; isEdge uses Alloc halos for neighbours.
             local_lower%x = max(local_lower%x, sgg%Sweep(field)%XI)
             local_lower%y = max(local_lower%y, sgg%Sweep(field)%YI)
             local_lower%z = max(local_lower%z, sgg%Sweep(field)%ZI)
             local_upper%x = min(local_upper%x, sgg%Sweep(field)%XE)
             local_upper%y = min(local_upper%y, sgg%Sweep(field)%YE)
             local_upper%z = min(local_upper%z, sgg%Sweep(field)%ZE)
          end do
          has_data = local_lower%x <= local_upper%x .and. local_lower%y <= local_upper%y .and. &
                     local_lower%z <= local_upper%z
       end subroutine get_local_map_bounds

      function preprocess_domain(observation, timeArray, simulationTimeStep, finalStepIndex) result(newDomain)
         type(Obses_t), intent(in) :: observation
         real(kind=RKIND_tiempo), pointer, dimension(:), intent(in) :: timeArray
         real(kind=RKIND_tiempo), intent(in) :: simulationTimeStep
         integer(kind=4), intent(in) :: finalStepIndex
         type(domain_t) :: newDomain

         integer(kind=SINGLE) :: nFreq

         if (observation%TimeDomain .and. observation%FreqDomain) then
            nFreq = frequency_count(observation)
            newdomain = domain_t(real(observation%InitialTime, kind=RKIND_tiempo), &
                                 real(observation%FinalTime, kind=RKIND_tiempo), &
                                 real(observation%TimeStep, kind=RKIND_tiempo), &
                                 observation%InitialFreq, frequency_stop(observation, nFreq), nFreq, .false.)

         else if (observation%TimeDomain) then
            newdomain = domain_t(real(observation%InitialTime, kind=RKIND_tiempo), &
                                 real(observation%FinalTime, kind=RKIND_tiempo), &
                                 real(observation%TimeStep, kind=RKIND_tiempo))

            newdomain%tstep = max(newdomain%tstep, simulationTimeStep)

            if (10.0_RKIND*(newdomain%tstop - newdomain%tstart)/min(simulationTimeStep, newdomain%tstep) >= huge(1_4)) then
               newdomain%tstop = newdomain%tstart + min(simulationTimeStep, newdomain%tstep)*huge(1_4)/10.0_RKIND
            end if

            if (newDomain%tstart < newDomain%tstep) then
               newDomain%tstart = 0.0_RKIND_tiempo
            end if

            if (newDomain%tstep > (newdomain%tstop - newdomain%tstart)) then
               newDomain%tstop = newDomain%tstart + newDomain%tstep
            end if

         elseif (observation%FreqDomain) then
            nFreq = frequency_count(observation)
            newdomain = domain_t(observation%InitialFreq, frequency_stop(observation, nFreq), nFreq, &
                                 logarithmicspacing=.false.)

         else
            newDomain = domain_t()
         end if
         return
      end function preprocess_domain

      function preprocess_polar_range(observation) result(sphericDomain)
         type(spheric_domain_t) :: sphericDomain
         type(Obses_t), intent(in) :: observation

         sphericDomain%phiStart = observation%phiStart
         sphericDomain%phiStop = observation%phiStop
         sphericDomain%phiStep = observation%phiStep
         sphericDomain%thetaStart = observation%thetaStart
         sphericDomain%thetaStop = observation%thetaStop
         sphericDomain%thetaStep = observation%thetaStep
      end function preprocess_polar_range

   end subroutine init_outputs

   subroutine update_outputs(control, discreteTimeArray, timeIndx, fieldsReference, sgg)
      integer(kind=SINGLE), intent(in) :: timeIndx
      real(kind=RKIND_tiempo), dimension(:), intent(in) :: discreteTimeArray
      integer(kind=SINGLE) :: i, id
      type(sim_control_t), intent(in) :: control
      real(kind=RKIND), pointer, dimension(:, :, :) :: fieldComponent
      type(field_data_t) :: fieldReference
      type(fields_reference_t), intent(in) :: fieldsReference
      type(SGGFDTDINFO_t), intent(in), optional :: sgg
      real(kind=RKIND_tiempo) :: discreteTime

      ! Assumed-shape arguments are indexed from one, while the solver time
      ! vector begins at zero. Account for that remapping at the API boundary.
      discreteTime = discreteTimeArray(timeIndx + 1)

      do i = 1, size(outputs)
         select case (outputs(i)%outputID)
         case (POINT_PROBE_ID)
            fieldComponent => get_field_component(outputs(i)%pointProbe%component, fieldsReference) !Cada componente requiere de valores deiferentes pero estos valores no se como conseguirlos
            if (present(sgg)) then
               call update_solver_output(outputs(i)%pointProbe, discreteTime, fieldComponent, sgg)
            else
               call update_solver_output(outputs(i)%pointProbe, discreteTime, fieldComponent)
            end if
         case (WIRE_CURRENT_PROBE_ID)
            call update_solver_output(outputs(i)%wireCurrentProbe, discreteTime, control, InvEps, InvMu)
         case (WIRE_CHARGE_PROBE_ID)
            call update_solver_output(outputs(i)%wireChargeProbe, discreteTime)
         case (BULK_PROBE_ID)
            fieldReference = get_field_reference(outputs(i)%bulkCurrentProbe%component, fieldsReference)
            call update_solver_output(outputs(i)%bulkCurrentProbe, discreteTime, fieldReference)
         case (LINE_PROBE_ID)
            call update_solver_output(outputs(i)%lineProbe, discreteTime, fieldsReference%E)
         case (MOVIE_PROBE_ID)
            call update_solver_output(outputs(i)%movieProbe, discreteTime, fieldsReference, control, problemInfo)
         case (FREQUENCY_SLICE_PROBE_ID)
            call update_solver_output(outputs(i)%frequencySliceProbe, discreteTime, fieldsReference, control, problemInfo)
         case (FAR_FIELD_PROBE_ID)
            call update_solver_output(outputs(i)%farFieldOutput, timeIndx, problemInfo%simulationBounds, fieldsReference)
         case (MAPVTK_ID)
         case default
            call stoponerror(0, 0, 'Output update not implemented')
         end select
      end do

   end subroutine update_outputs

   subroutine flush_outputs(simulationTimeArray, simulationTimeIndex, control, fields, bounds, farFieldFlushRequested)
      implicit none
      type(fields_reference_t), target :: fields
      type(fields_reference_t), pointer :: fieldsPtr
      type(sim_control_t), intent(in) :: control
      type(bounds_t), intent(in) :: bounds
      logical, intent(in) :: farFieldFlushRequested
      real(KIND=RKIND_tiempo), pointer, dimension(:), intent(in) :: simulationTimeArray
      integer, intent(in) :: simulationTimeIndex
      integer :: outIdx

      fieldsPtr => fields

      do outIdx = 1, size(outputs)
         select case (outputs(outIdx)%outputID)
         case (POINT_PROBE_ID)
            call flush_solver_output(outputs(outIdx)%pointProbe)
         case (WIRE_CURRENT_PROBE_ID)
            call flush_solver_output(outputs(outIdx)%wireCurrentProbe)
         case (WIRE_CHARGE_PROBE_ID)
            call flush_solver_output(outputs(outIdx)%wireChargeProbe)
         case (BULK_PROBE_ID)
            call flush_solver_output(outputs(outIdx)%bulkCurrentProbe)
         case (LINE_PROBE_ID)
            call flush_solver_output(outputs(outIdx)%lineProbe)
         case (MOVIE_PROBE_ID)
            call flush_solver_output(outputs(outIdx)%movieProbe)
         case (FREQUENCY_SLICE_PROBE_ID)
            call flush_solver_output(outputs(outIdx)%frequencySliceProbe)
         case (FAR_FIELD_PROBE_ID)
            if (farFieldFlushRequested) call flush_solver_output(outputs(outIdx)%farFieldOutput, simulationTimeArray, simulationTimeIndex, control, fieldsPtr, bounds)
         case default
         end select
      end do
   end subroutine flush_outputs

   subroutine remove_unused_outputs(output_list)
      implicit none
      type(solver_output_t), pointer, intent(inout) :: output_list(:)

      type(solver_output_t), allocatable :: tmp(:)
      type(output_partition_t), allocatable :: compact_partitions(:)
      integer :: i, n, k

      n = count(output_list%outputID /= UNDEFINED_PROBE)

      allocate (tmp(n))
      allocate (compact_partitions(n))

      ! Copy valid elements
      k = 0
      do i = 1, size(output_list)
         if (output_list(i)%outputID /= UNDEFINED_PROBE) then
            k = k + 1
            tmp(k) = output_list(i)   ! deep copy of all allocatable components
            compact_partitions(k) = outputPartitions(i)
         end if
      end do

      ! Replace the saved pointer target safely
      if (associated(output_list)) deallocate (output_list)
      allocate (output_list(n))
      output_list = tmp
      call move_alloc(compact_partitions, outputPartitions)

   end subroutine remove_unused_outputs

   subroutine close_outputs()
      integer :: i
      do i = 1, size(outputs)
         select case (outputs(i)%outputID)
         case (POINT_PROBE_ID)
            call finalise_scalar_output_metadata(i)
         case (WIRE_CURRENT_PROBE_ID)
            call finalise_scalar_output_metadata(i)
         case (WIRE_CHARGE_PROBE_ID)
            call finalise_scalar_output_metadata(i)
         case (BULK_PROBE_ID)
            call finalise_scalar_output_metadata(i)
         case (LINE_PROBE_ID)
            call finalise_scalar_output_metadata(i)
         case (VOLUMIC_CURRENT_PROBE_ID)
            call finalise_scalar_output_metadata(i)
         case (FAR_FIELD_PROBE_ID)
            call finalise_scalar_output_metadata(i)
         case (MAPVTK_ID)
            call finalise_scalar_output_metadata(i)
         case (MOVIE_PROBE_ID)
            call close_solver_output(outputs(i)%movieProbe)
         case (FREQUENCY_SLICE_PROBE_ID)
            call close_solver_output(outputs(i)%frequencySliceProbe)
         end select
      end do
   end subroutine

   subroutine create_pvd(pvdPath)
      implicit none
      character(len=*), intent(in) :: pvdPath
      integer :: ios
      integer :: unit

      open (newunit=unit, file=trim(pvdPath), status="replace", action="write", iostat=ios)
      if (ios /= 0) stop "Error al crear archivo PVD"

      ! Escribimos encabezados XML
      write (unit, *) '<?xml version="1.0"?>'
      write (unit, *) '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">'
      write (unit, *) '  <Collection>'
      close (unit)
   end subroutine create_pvd

   subroutine close_pvd(pvdPath)
      implicit none
      character(len=*), intent(in) :: pvdPath
      integer :: unit
      integer :: ios
      open (newunit=unit, file=trim(pvdPath), status="old", action="write", iostat=ios)
      if (ios /= 0) stop "Error al abrir archivo PVD"
      write (unit, *) '  </Collection>'
      write (unit, *) '</VTKFile>'
      close (unit)
   end subroutine close_pvd

   function get_required_output_count(sgg) result(count)
      type(SGGFDTDINFO_t), intent(in) :: sgg
      integer(kind=SINGLE) ::i, count
      count = 0
      do i = 1, sgg%NumberRequest
         count = count + sgg%Observation(i)%nP
      end do
      return
   end function

   subroutine write_run_output_manifest(control, outputCount)
      type(sim_control_t), intent(in) :: control
      integer, intent(in) :: outputCount

      character(len=BUFSIZE) :: manifest_path
      integer :: i, iostat, unit
      logical :: first_artifact

      if (control%layoutnumber /= 0) return
      manifest_path = trim(control%nEntradaRoot)//'_output_manifest.json'
      open (newunit=unit, file=trim(manifest_path), status='replace', action='write', iostat=iostat)
      if (iostat /= 0) then
         call StopOnError(control%layoutnumber, control%num_procs, 'Error while creating output manifest')
         return
      end if

      first_artifact = .true.
      write (unit, '(a)') '{'
      write (unit, '(a)') '"artifacts":['
      do i = 1, outputCount
         select case (outputs(i)%outputID)
         case (POINT_PROBE_ID)
            call write_artifacts(outputs(i)%pointProbe%artifacts)
         case (WIRE_CURRENT_PROBE_ID)
            call write_artifacts(outputs(i)%wireCurrentProbe%artifacts)
         case (WIRE_CHARGE_PROBE_ID)
            call write_artifacts(outputs(i)%wireChargeProbe%artifacts)
         case (BULK_PROBE_ID)
            call write_artifacts(outputs(i)%bulkCurrentProbe%artifacts)
         case (LINE_PROBE_ID)
            call write_artifacts(outputs(i)%lineProbe%artifacts)
         case (MOVIE_PROBE_ID)
            call write_artifacts(outputs(i)%movieProbe%metadata%artifacts, outputs(i)%movieProbe%path)
         case (FREQUENCY_SLICE_PROBE_ID)
            call write_artifacts(outputs(i)%frequencySliceProbe%metadata%artifacts, outputs(i)%frequencySliceProbe%path)
         case (FAR_FIELD_PROBE_ID)
            call write_artifacts(outputs(i)%farFieldOutput%artifacts)
         case (MAPVTK_ID)
            call write_artifacts(outputs(i)%mapvtkOutput%artifacts)
         end select
      end do
      write (unit, '(a)') ']'
      write (unit, '(a)') '}'
      close (unit)
   contains
      subroutine write_artifacts(artifacts, base_path)
         type(output_artifact_t), intent(in) :: artifacts(:)
         character(len=*), intent(in), optional :: base_path
         character(len=:), allocatable :: artifact_path
         integer :: artifact_index

         do artifact_index = 1, size(artifacts)
            if (artifacts(artifact_index)%kind == OUTPUT_ARTIFACT_UNDEFINED) cycle
            artifact_path = trim(artifacts(artifact_index)%relative_path)
            if (present(base_path)) artifact_path = join_path(trim(base_path), artifact_path)
            if (.not. first_artifact) write (unit, '(a)') ','
            write (unit, '(a)') '{"path":"'//json_escape(trim(artifact_path))//'"}'
            first_artifact = .false.
         end do
      end subroutine write_artifacts
   end subroutine write_run_output_manifest

   subroutine delete_run_output_manifest(run_id, writer_rank)
      character(len=*), intent(in) :: run_id
      integer, intent(in) :: writer_rank
      character(len=BUFSIZE) :: artifact_path, artifact_folder, line, manifest_path
      integer :: ios, path_end, path_start, unit

      if (writer_rank /= 0) return
      manifest_path = trim(run_id)//'_output_manifest.json'
      if (.not. file_exists(manifest_path)) return
      open (newunit=unit, file=trim(manifest_path), status='old', action='read', iostat=ios)
      if (ios /= 0) return
      do
         read (unit, '(A)', iostat=ios) line
         if (ios /= 0) exit
         path_start = index(line, '"path":"')
         if (path_start == 0) cycle
         path_start = path_start + len('"path":"')
         path_end = index(line(path_start:), '"')
         if (path_end == 0) cycle
         artifact_path = json_unescape(line(path_start:path_start + path_end - 2))
         if (index(trim(artifact_path), trim(run_id)//'_') == 1) then
            call delete_file(trim(artifact_path), ios)
            path_end = scan(trim(artifact_path), '/\\', back=.true.)
            if (path_end > 0) then
               artifact_folder = artifact_path(:path_end - 1)
               call remove_folder(trim(artifact_folder), ios)
            end if
         end if
      end do
      close (unit)
      call delete_file(trim(manifest_path), ios)
   end subroutine delete_run_output_manifest

   function frequency_count(observation) result(count)
      type(Obses_t), intent(in) :: observation
      integer(kind=SINGLE) :: count

      if (observation%FreqStep <= 0.0_RKIND) then
         count = 1_SINGLE
      else
         count = int(floor((observation%FinalFreq - observation%InitialFreq)/observation%FreqStep), &
                     kind=SINGLE) + 1_SINGLE
      end if
   end function frequency_count

   function frequency_stop(observation, count) result(stop)
      type(Obses_t), intent(in) :: observation
      integer(kind=SINGLE), intent(in) :: count
      real(kind=RKIND) :: stop

      if (count == 1_SINGLE) then
         stop = observation%InitialFreq
      else
         stop = observation%InitialFreq + real(count - 1_SINGLE, RKIND)*observation%FreqStep
      end if
   end function frequency_stop

end module output_m
