module output_m
   use FDETYPES_m
   use report_m
   use domain_m
   use outputUtils_m
   use directoryUtils_m, only: delete_file, remove_folder
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
   use outputCollective_m, only: output_collective_t, init_output_collective, select_output_participants, &
                                 prepare_output_partition_publication, OUTPUT_COLLECTIVE_SUCCESS
   use outputTypes_m, only: solver_output_t, problem_info_t, output_artifact_t, &
                              spheric_domain_t, cell_coordinate_t, field_data_t, fields_reference_t, &
                              TIME_DOMAIN, FREQUENCY_DOMAIN, BOTH_DOMAIN, UNDEFINED_DOMAIN, volumetric_publication_t, &
                              OUTPUT_ARTIFACT_UNDEFINED
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
   public :: delete_outputs
   public :: GetOutputPartition
   public :: POINT_PROBE_ID, WIRE_CURRENT_PROBE_ID, WIRE_CHARGE_PROBE_ID, BULK_PROBE_ID, &
             MOVIE_PROBE_ID, FREQUENCY_SLICE_PROBE_ID, FAR_FIELD_PROBE_ID, LINE_PROBE_ID, MTLN_PROBE_ID
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
                                      MOVIE_PROBE_ID = 5, &
                                      FREQUENCY_SLICE_PROBE_ID = 6, &
                                      FAR_FIELD_PROBE_ID = 7, &
                                      MAPVTK_ID = 8, LINE_PROBE_ID = 9, MTLN_PROBE_ID = 10

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
      type(cell_coordinate_t) :: globalMapLower, globalMapUpper
      integer(kind=SINGLE) :: i, ii, outputRequestType
      integer(kind=SINGLE) :: NODE
      integer(kind=SINGLE) :: outputCount
      integer(kind=SINGLE) :: requestedOutputs
      logical :: mapHasData
      character(len=BUFSIZE) :: outputTypeExtension
#ifdef CompileWithMPI
      character(len=BUFSIZE) :: localMapPath
      character(len=BUFSIZE), allocatable :: mapPiecePaths(:)
      integer :: localMapMinimum(3), localMapMaximum(3), globalMapMinimum(3), globalMapMaximum(3)
      integer :: mpiError
#endif

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
                  if (mtln_solver%bundles(i)%probes(j)%output_writer_rank >= 0) then
                     thereAreMtlnObservations = .true.
                     exit mtlnObservationExist
                  end if
               end do
            end do mtlnObservationExist
         end if
      end block
#endif

      do ii = 1, sgg%NumberRequest
         domain = preprocess_domain(sgg%Observation(ii), sgg%tiempo, sgg%dt, control%finaltimestep, &
                                    control%saveall)
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
               globalMapLower = localMapLower
               globalMapUpper = localMapUpper
#ifdef CompileWithMPI
               if (control%num_procs > 1) then
                  localMapMinimum = huge(0)
                  localMapMaximum = -huge(0)
                  if (mapHasData) then
                     localMapMinimum = [localMapLower%x, localMapLower%y, localMapLower%z]
                     localMapMaximum = [localMapUpper%x, localMapUpper%y, localMapUpper%z]
                  end if
                  call MPI_Allreduce(localMapMinimum, globalMapMinimum, 3, MPI_INTEGER, MPI_MIN, &
                                     SUBCOMM_MPI, mpiError)
                  if (mpiError /= MPI_SUCCESS) then
                     call StopOnError(control%layoutnumber, control%num_procs, &
                                      'Unable to determine distributed geometry map lower bounds')
                  end if
                  call MPI_Allreduce(localMapMaximum, globalMapMaximum, 3, MPI_INTEGER, MPI_MAX, &
                                     SUBCOMM_MPI, mpiError)
                  if (mpiError /= MPI_SUCCESS) then
                     call StopOnError(control%layoutnumber, control%num_procs, &
                                      'Unable to determine distributed geometry map upper bounds')
                  end if
                  globalMapLower = cell_coordinate_t(globalMapMinimum(1), globalMapMinimum(2), globalMapMinimum(3))
                  globalMapUpper = cell_coordinate_t(globalMapMaximum(1), globalMapMaximum(2), globalMapMaximum(3))
               end if
#endif
               outputCount = outputCount + 1
               outputs(outputCount)%outputID = MAPVTK_ID
               allocate (outputs(outputCount)%mapvtkOutput)
               call init_solver_output(outputs(outputCount)%mapvtkOutput, localMapLower, localMapUpper, &
                                       globalMapLower, globalMapUpper, mapHasData, outputRequestType, &
                                       outputTypeExtension, control%mpidir, problemInfo)
               if (mapHasData) then
                  call create_geometry_simulation_vtu(outputs(outputCount)%mapvtkOutput, control, sgg%LineX, sgg%LineY, &
                                                      sgg%LineZ, problemInfo)
               end if
#ifdef CompileWithMPI
               if (control%num_procs > 1) then
                  localMapPath = ''
                  if (mapHasData) localMapPath = outputs(outputCount)%mapvtkOutput%artifacts(1)%relative_path
                  allocate (mapPiecePaths(control%num_procs))
                  call MPI_Allgather(localMapPath, BUFSIZE, MPI_CHARACTER, mapPiecePaths, BUFSIZE, &
                                     MPI_CHARACTER, SUBCOMM_MPI, mpiError)
                  if (mpiError /= MPI_SUCCESS) then
                     call StopOnError(control%layoutnumber, control%num_procs, &
                                      'Unable to collect distributed geometry map paths')
                  end if
                  call MPI_Barrier(SUBCOMM_MPI, mpiError)
                  if (mpiError /= MPI_SUCCESS) then
                     call StopOnError(control%layoutnumber, control%num_procs, &
                                      'Unable to synchronize distributed geometry map output')
                  end if
                  if (control%layoutnumber == 0) then
                     call create_parallel_geometry_vtu(outputs(outputCount)%mapvtkOutput, mapPiecePaths)
                  end if
                  deallocate (mapPiecePaths)
               end if
#endif

            case (iEx, iEy, iEz, iHx, iHy, iHz)
               outputCount = outputCount + 1
               outputs(outputCount)%outputID = POINT_PROBE_ID

               allocate (outputs(outputCount)%pointProbe)
                 call init_solver_output(outputs(outputCount)%pointProbe, lowerBound, outputRequestType, domain, outputTypeExtension, control%mpidir, sgg%dt, &
                                       sgg%NumPlaneWaves >= 1)
            case (iJx, iJy, iJz)
               if (wiresExists) then
                  outputCount = outputCount + 1
                  outputs(outputCount)%outputID = WIRE_CURRENT_PROBE_ID

                  allocate (outputs(outputCount)%wireCurrentProbe)
                   call init_solver_output(outputs(outputCount)%wireCurrentProbe, lowerBound, NODE, outputRequestType, domain, problemInfo%materialList, outputTypeExtension, control%mpidir, control%wiresflavor)
               end if

            case (iQx, iQy, iQz)
               outputCount = outputCount + 1
               outputs(outputCount)%outputID = WIRE_CHARGE_PROBE_ID

               allocate (outputs(outputCount)%wireChargeProbe)
                 call init_solver_output(outputs(outputCount)%wireChargeProbe, lowerBound, NODE, outputRequestType, domain, outputTypeExtension, control%mpidir, control%wiresflavor)

            case (iBloqueJx, iBloqueJy, iBloqueJz, iBloqueMx, iBloqueMy, iBloqueMz)
               outputCount = outputCount + 1
               outputs(outputCount)%outputID = BULK_PROBE_ID

               allocate (outputs(outputCount)%bulkCurrentProbe)
                 call init_solver_output(outputs(outputCount)%bulkCurrentProbe, lowerBound, upperBound, outputRequestType, domain, outputTypeExtension, control%mpidir)
                 !! call adjust_computation_range --- Required due to issues in mpi region edges

            case (lineIntegral)
               if (domain%domainType /= TIME_DOMAIN) then
                  call stoponerror(0, 0, 'Line probes only support the time domain')
               else
                  outputCount = outputCount + 1
                  outputs(outputCount)%outputID = LINE_PROBE_ID
                  allocate (outputs(outputCount)%lineProbe)
                  call init_line_probe_output(outputs(outputCount)%lineProbe, sgg%observation(ii)%P(i)%line, domain, &
                                              trim(outputTypeExtension)//'_LI', &
                                              sgg%Sweep, control%layoutnumber, control%num_procs)
               end if

            case (iCur, iMEC, iMHC, iCurX, iCurY, iCurZ, iExC, iEyC, iEzC, iHxC, iHyC, iHzC)
               if (domain%domainType == TIME_DOMAIN) then
                  block
                     type(volumetric_publication_t) :: publication
                     outputCount = outputCount + 1
                     outputs(outputCount)%outputID = MOVIE_PROBE_ID
                     call attach_output_partition(outputCount)
                     call configure_output_publication(outputCount, publication)
                     allocate (outputs(outputCount)%movieProbe)
                     call init_solver_output(outputs(outputCount)%movieProbe, publication, outputRequestType, domain, &
                                             control, problemInfo, outputTypeExtension)
                  end block
               else if (domain%domainType == FREQUENCY_DOMAIN) then
                  block
                     type(volumetric_publication_t) :: publication
                     outputCount = outputCount + 1
                     outputs(outputCount)%outputID = FREQUENCY_SLICE_PROBE_ID
                     call attach_output_partition(outputCount)
                     call configure_output_publication(outputCount, publication)
                     allocate (outputs(outputCount)%frequencySliceProbe)
                     call init_solver_output(outputs(outputCount)%frequencySliceProbe, publication, sgg%dt, &
                                             outputRequestType, domain, outputTypeExtension, control, problemInfo)
                  end block
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
      end subroutine configure_far_field_mpi
#endif

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

      subroutine configure_output_publication(output_index, publication)
         integer(kind=SINGLE), intent(in) :: output_index
         type(volumetric_publication_t), intent(out) :: publication
         type(output_collective_t) :: collective
         integer, allocatable :: participants(:)
         logical, allocatable :: rank_has_data(:)
         integer :: collective_status, owner_rank
         logical :: parallel_hdf5_available
#ifdef CompileWithMPI
         integer :: color, ierr
#endif

         publication = volumetric_publication_t()
         allocate (rank_has_data(max(control%num_procs, 1)))
#ifdef CompileWithMPI
         if (control%num_procs > 1) then
            call MPI_Allgather(outputPartitions(output_index)%has_data, 1, MPI_LOGICAL, rank_has_data, 1, &
                               MPI_LOGICAL, SUBCOMM_MPI, ierr)
            if (ierr /= MPI_SUCCESS) then
               call StopOnError(control%layoutnumber, control%num_procs, &
                                'Unable to identify distributed output participants')
            end if
         else
            rank_has_data(1) = outputPartitions(output_index)%has_data
         end if
#else
         rank_has_data(1) = outputPartitions(output_index)%has_data
#endif

         parallel_hdf5_available = .false.
#ifdef XDMF_HDF5_PARALLEL_AVAILABLE
         parallel_hdf5_available = .true.
#endif
         call init_output_collective(collective, control%layoutnumber, max(control%num_procs, 1), 0, &
                                     parallel_hdf5_available, collective_status)
         if (collective_status == OUTPUT_COLLECTIVE_SUCCESS) then
            call select_output_participants(collective, rank_has_data, participants, owner_rank, collective_status)
            collective%collective_publication_available = parallel_hdf5_available .and. size(participants) > 1
         end if
         if (collective_status == OUTPUT_COLLECTIVE_SUCCESS) then
            call prepare_output_partition_publication(collective, participants, owner_rank, &
                                                      outputPartitions(output_index), &
                                                      publication%local_participates, publication%mode, &
                                                      collective_status)
         end if
         if (collective_status /= OUTPUT_COLLECTIVE_SUCCESS) then
            call StopOnError(control%layoutnumber, control%num_procs, &
                             'Unable to configure distributed output publication')
         end if

         publication%participant_ranks = participants
         publication%owner_rank = owner_rank
         publication%local_is_owner = control%layoutnumber == owner_rank
         publication%file_offset = outputPartitions(output_index)%file_offset
         publication%local_shape = outputPartitions(output_index)%local_shape
         publication%global_lower = outputPartitions(output_index)%global_lower
         publication%global_upper = outputPartitions(output_index)%global_upper

#ifdef CompileWithMPI
         if (control%num_procs > 1) then
            color = MPI_UNDEFINED
            if (publication%local_participates) color = 0
            call MPI_Comm_split(SUBCOMM_MPI, color, control%layoutnumber, publication%communicator, ierr)
            if (ierr /= MPI_SUCCESS) then
               call StopOnError(control%layoutnumber, control%num_procs, &
                                'Unable to create distributed output communicator')
            end if
            if (publication%local_participates) then
               publication%owns_communicator = .true.
               call MPI_Comm_rank(publication%communicator, publication%communicator_rank, ierr)
               call MPI_Comm_size(publication%communicator, publication%communicator_size, ierr)
            end if
         else
            publication%communicator = MPI_COMM_WORLD
         end if
#endif
      end subroutine configure_output_publication

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

       function preprocess_domain(observation, timeArray, simulationTimeStep, finalStepIndex, globalSaveAll) &
          result(newDomain)
          type(Obses_t), intent(in) :: observation
          real(kind=RKIND_tiempo), pointer, dimension(:), intent(in) :: timeArray
          real(kind=RKIND_tiempo), intent(in) :: simulationTimeStep
          integer(kind=4), intent(in) :: finalStepIndex
          logical, intent(in) :: globalSaveAll
          type(domain_t) :: newDomain

          integer(kind=SINGLE) :: nFreq
          integer :: simulationEndIndex
          logical :: saveAllTimeSteps

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

          elseif (observation%FreqDomain) then
            nFreq = frequency_count(observation)
            newdomain = domain_t(observation%InitialFreq, frequency_stop(observation, nFreq), nFreq, &
                                 logarithmicspacing=.false.)

          else
             newDomain = domain_t()
          end if

          if (any(newDomain%domainType == [TIME_DOMAIN, BOTH_DOMAIN])) then
             simulationEndIndex = min(max(finalStepIndex + 2, lbound(timeArray, 1)), ubound(timeArray, 1))
             saveAllTimeSteps = globalSaveAll .or. observation%FinalTime < tiny(1.0_RKIND) .or. &
                                observation%TimeStep < tiny(1.0_RKIND)
             if (saveAllTimeSteps) then
                newDomain%tstart = 0.0_RKIND_tiempo
                newDomain%tstop = timeArray(simulationEndIndex)
                newDomain%tstep = simulationTimeStep
             else
                newDomain%tstart = max(0.0_RKIND_tiempo, newDomain%tstart)
                newDomain%tstop = min(timeArray(simulationEndIndex), newDomain%tstop)
                newDomain%tstop = max(newDomain%tstart, newDomain%tstop)
                newDomain%tstep = max(simulationTimeStep, newDomain%tstep)
             end if
             newDomain%tstride = max(1_SINGLE, int(newDomain%tstep/simulationTimeStep, kind=SINGLE))
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

   subroutine update_outputs(control, discreteTime, timeIndx, fieldsReference, sgg)
      integer(kind=SINGLE), intent(in) :: timeIndx
      real(kind=RKIND_tiempo), intent(in) :: discreteTime
      integer(kind=SINGLE) :: i, id
      type(sim_control_t), intent(in) :: control
      real(kind=RKIND), pointer, dimension(:, :, :) :: fieldComponent
      type(field_data_t) :: fieldReference
      type(fields_reference_t), intent(in) :: fieldsReference
      type(SGGFDTDINFO_t), intent(in), optional :: sgg
      logical :: timeSampleDue

      do i = 1, size(outputs)
         select case (outputs(i)%outputID)
         case (POINT_PROBE_ID)
            timeSampleDue = is_time_sample_due(outputs(i)%pointProbe%domain, timeIndx, discreteTime)
            fieldComponent => get_field_component(outputs(i)%pointProbe%component, fieldsReference) !Cada componente requiere de valores deiferentes pero estos valores no se como conseguirlos
            if (present(sgg)) then
               call update_solver_output(outputs(i)%pointProbe, discreteTime, fieldComponent, sgg, timeSampleDue)
            else
               call update_solver_output(outputs(i)%pointProbe, discreteTime, fieldComponent, &
                                         saveTimeSample=timeSampleDue)
            end if
         case (WIRE_CURRENT_PROBE_ID)
            if (.not. is_time_sample_due(outputs(i)%wireCurrentProbe%domain, timeIndx, discreteTime)) cycle
            call update_solver_output(outputs(i)%wireCurrentProbe, discreteTime, control, InvEps, InvMu)
         case (WIRE_CHARGE_PROBE_ID)
            if (.not. is_time_sample_due(outputs(i)%wireChargeProbe%domain, timeIndx, discreteTime)) cycle
            call update_solver_output(outputs(i)%wireChargeProbe, discreteTime)
         case (BULK_PROBE_ID)
            if (.not. is_time_sample_due(outputs(i)%bulkCurrentProbe%domain, timeIndx, discreteTime)) cycle
            fieldReference = get_field_reference(outputs(i)%bulkCurrentProbe%component, fieldsReference)
            call update_solver_output(outputs(i)%bulkCurrentProbe, discreteTime, fieldReference)
         case (LINE_PROBE_ID)
            if (.not. is_time_sample_due(outputs(i)%lineProbe%domain, timeIndx, discreteTime)) cycle
            call update_solver_output(outputs(i)%lineProbe, discreteTime, fieldsReference%E)
         case (MOVIE_PROBE_ID)
            if (.not. is_time_sample_due(outputs(i)%movieProbe%domain, timeIndx, discreteTime)) cycle
            call update_solver_output(outputs(i)%movieProbe, discreteTime, fieldsReference, control, problemInfo)
         case (FREQUENCY_SLICE_PROBE_ID)
            call update_solver_output(outputs(i)%frequencySliceProbe, discreteTime, fieldsReference, control, problemInfo)
         case (FAR_FIELD_PROBE_ID)
            call update_solver_output(outputs(i)%farFieldOutput, timeIndx, problemInfo%simulationBounds, fieldsReference)
         case (MAPVTK_ID)
         case (MTLN_PROBE_ID)
         case default
            call stoponerror(0, 0, 'Output update not implemented')
         end select
      end do

   end subroutine update_outputs

   logical function is_time_sample_due(domain, timeIndex, discreteTime) result(sampleDue)
      type(domain_t), intent(in) :: domain
      integer(kind=SINGLE), intent(in) :: timeIndex
      real(kind=RKIND_tiempo), intent(in) :: discreteTime
      real(kind=RKIND_tiempo) :: boundaryTolerance, timeScale

      sampleDue = .false.
      if (.not. any(domain%domainType == [TIME_DOMAIN, BOTH_DOMAIN])) return
      if (modulo(timeIndex, domain%tstride) /= 0) return

      timeScale = max(abs(discreteTime), abs(domain%tstart), abs(domain%tstop), domain%tstep)
      boundaryTolerance = 4.0_RKIND_tiempo*real(epsilon(1.0_RKIND), RKIND_tiempo)*timeScale
      sampleDue = discreteTime >= domain%tstart - boundaryTolerance .and. &
                  discreteTime <= domain%tstop + boundaryTolerance
   end function is_time_sample_due

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
         case (MOVIE_PROBE_ID)
            call close_solver_output(outputs(i)%movieProbe)
         case (FREQUENCY_SLICE_PROBE_ID)
            call close_solver_output(outputs(i)%frequencySliceProbe)
         end select
      end do
   end subroutine

   subroutine delete_outputs(writer_rank)
      integer, intent(in) :: writer_rank
      integer :: i, ios

      if (associated(outputs)) then
         do i = 1, size(outputs)
            if (outputs(i)%outputID /= MAPVTK_ID) cycle
            if (outputs(i)%mapvtkOutput%localParticipates) then
               call remove_folder(outputs(i)%mapvtkOutput%path, ios)
            end if
            if (writer_rank == 0) call delete_file(outputs(i)%mapvtkOutput%masterPath, ios)
         end do
      end if
      if (writer_rank /= 0) return
      if (associated(outputs)) then
         do i = 1, size(outputs)
            select case (outputs(i)%outputID)
            case (POINT_PROBE_ID)
               call delete_artifacts(outputs(i)%pointProbe%artifacts)
            case (WIRE_CURRENT_PROBE_ID)
               call delete_artifacts(outputs(i)%wireCurrentProbe%artifacts)
            case (WIRE_CHARGE_PROBE_ID)
               call delete_artifacts(outputs(i)%wireChargeProbe%artifacts)
            case (BULK_PROBE_ID)
               call delete_artifacts(outputs(i)%bulkCurrentProbe%artifacts)
            case (LINE_PROBE_ID)
               call delete_artifacts(outputs(i)%lineProbe%artifacts)
            case (FAR_FIELD_PROBE_ID)
               if (allocated(outputs(i)%farFieldOutput%artifacts)) then
                  call delete_artifacts(outputs(i)%farFieldOutput%artifacts)
               end if
            case (MOVIE_PROBE_ID)
               call remove_folder(outputs(i)%movieProbe%path, ios)
            case (FREQUENCY_SLICE_PROBE_ID)
               call remove_folder(outputs(i)%frequencySliceProbe%path, ios)
            end select
         end do
      end if
#ifdef CompileWithMTLN
      call delete_mtln_probe_outputs()
#endif
   end subroutine delete_outputs

   subroutine delete_artifacts(artifacts)
      type(output_artifact_t), intent(in) :: artifacts(:)
      integer :: i, ios

      do i = 1, size(artifacts)
         if (artifacts(i)%kind == OUTPUT_ARTIFACT_UNDEFINED .or. &
             len_trim(artifacts(i)%relative_path) == 0) cycle
         call delete_file(artifacts(i)%relative_path, ios)
      end do
   end subroutine delete_artifacts

#ifdef CompileWithMTLN
   subroutine delete_mtln_probe_outputs()
      type(mtln_solver_t), pointer :: mtln_solver
      integer :: i, ios, j, separator

      mtln_solver => GetSolverPtr()
      if (.not. allocated(mtln_solver%bundles)) return
      do i = 1, size(mtln_solver%bundles)
         if (.not. allocated(mtln_solver%bundles(i)%probes)) cycle
         do j = 1, size(mtln_solver%bundles(i)%probes)
            call delete_file(mtln_solver%bundles(i)%probes(j)%output_path, ios)
            separator = scan(trim(mtln_solver%bundles(i)%probes(j)%output_path), '/\', back=.true.)
            if (separator > 1) then
               call remove_folder(mtln_solver%bundles(i)%probes(j)%output_path(:separator - 1), ios)
            end if
         end do
      end do
   end subroutine delete_mtln_probe_outputs
#endif

   function get_required_output_count(sgg) result(count)
      type(SGGFDTDINFO_t), intent(in) :: sgg
      integer(kind=SINGLE) ::i, count
      count = 0
      do i = 1, sgg%NumberRequest
         count = count + sgg%Observation(i)%nP
      end do
      return
   end function

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
