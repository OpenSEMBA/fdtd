program test_output_mpi
   use mpi
    use FDETYPES_m, only: iEx, iEz, limit_t
   use outputCollective_m, only: output_collective_t, init_output_collective, &
                                 select_output_participants, &
                                 prepare_output_partition_publication, &
                                 OUTPUT_COLLECTIVE_SUCCESS, &
                                 OUTPUT_PUBLICATION_COLLECTIVE, &
                                 OUTPUT_PUBLICATION_ROOT_AGGREGATION
    use outputDecomposition_m, only: output_partition_t, build_output_partition, &
                                     OUTPUT_PARTITION_SUCCESS
    use outputTransport_m, only: output_transport_t, init_output_transport, OUTPUT_TRANSPORT_SUCCESS
    use outputTypes_m, only: cell_coordinate_t, output_artifact_t, OUTPUT_ARTIFACT_TEXT
    use output_m, only: run_output_manifest_t, init_run_output_manifest, declare_probe_output, begin_probe_output, &
                         finalise_probe_output, finalise_transport_run_outputs, OUTPUT_COORDINATION_SUCCESS
   implicit none

   integer, parameter :: root = 0
    integer :: ierr, rank, rank_count, status, owner, publication_mode
    integer :: local_count, total_count, i, probe_index, z, failures
   integer, allocatable :: participants(:), local_values(:), gathered_values(:), counts(:), displacements(:)
   integer, allocatable :: local_coverage(:), global_coverage(:)
   logical, allocatable :: rank_has_data(:)
     logical :: local_participates
    type(output_collective_t) :: collective
    type(output_transport_t) :: transport
    type(output_partition_t) :: partition
    type(cell_coordinate_t) :: request_lower, request_upper
    type(limit_t) :: global_bounds, local_sweep
    type(run_output_manifest_t) :: manifest
    type(output_artifact_t) :: artifacts(1)

   call MPI_Init(ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, rank_count, ierr)

   failures = 0
   if (rank_count < 2) then
      failures = 1
   else
      request_lower = cell_coordinate_t(0, 0, 0)
      request_upper = cell_coordinate_t(0, 0, 2 * rank_count)
      global_bounds = limit_t(0, 0, 0, 0, 0, 2 * rank_count, 1, 1, 2 * rank_count + 1)
      ! Adjacent rank sweeps share their upper/lower interface plane.
      local_sweep = limit_t(0, 0, 0, 0, 2 * rank, 2 * (rank + 1), 1, 1, 3)

      call build_output_partition(request_lower, request_upper, global_bounds, local_sweep, &
                                  iEx, rank, rank_count, partition, status)
      if (status /= OUTPUT_PARTITION_SUCCESS .or. .not. partition%has_data) failures = failures + 1

      allocate(rank_has_data(rank_count))
      rank_has_data = .true.
       call init_output_collective(collective, rank, rank_count, root, .true., status)
      if (status /= OUTPUT_COLLECTIVE_SUCCESS) failures = failures + 1
      call select_output_participants(collective, rank_has_data, participants, owner, status)
      if (status /= OUTPUT_COLLECTIVE_SUCCESS .or. owner /= root) failures = failures + 1
      call prepare_output_partition_publication(collective, participants, owner, partition, &
                                                local_participates, publication_mode, status)
       if (status /= OUTPUT_COLLECTIVE_SUCCESS .or. .not. local_participates .or. &
           publication_mode /= OUTPUT_PUBLICATION_COLLECTIVE) failures = failures + 1

       call init_output_transport(transport, root, status)
       if (status /= OUTPUT_TRANSPORT_SUCCESS) failures = failures + 1
       artifacts(1)%kind = OUTPUT_ARTIFACT_TEXT
       artifacts(1)%relative_path = 'point.dat'
       call init_run_output_manifest(manifest, 'mpi-run', root)
        call declare_probe_output(manifest, 'point-001', 'Ex', artifacts, probe_index, status)
       if (status /= OUTPUT_COORDINATION_SUCCESS) failures = failures + 1
       call begin_probe_output(manifest, probe_index, status)
       if (status /= OUTPUT_COORDINATION_SUCCESS) failures = failures + 1
       call finalise_probe_output(manifest, probe_index, status)
       if (status /= OUTPUT_COORDINATION_SUCCESS) failures = failures + 1
       call finalise_transport_run_outputs(manifest, transport, status)
       if (status /= OUTPUT_COORDINATION_SUCCESS) failures = failures + 1
       if (manifest%published .neqv. (rank == root)) failures = failures + 1

      allocate(local_coverage(0:2 * rank_count), global_coverage(0:2 * rank_count))
      local_coverage = 0
      do z = partition%local_lower%z, partition%local_upper%z
         local_coverage(z) = 1
      end do
      call MPI_Reduce(local_coverage, global_coverage, size(local_coverage), MPI_INTEGER, MPI_SUM, root, &
                      MPI_COMM_WORLD, ierr)
       if (rank == root) then
          do z = 0, 2 * rank_count
             if (global_coverage(z) /= 1) failures = failures + 1
          end do
       end if

       ! Electric-z ranges are already disjoint in MPIdivide and must retain
       ! every requested coordinate when used to initialise a movie fragment.
       local_sweep%ZI = 2 * rank
       local_sweep%ZE = 2 * (rank + 1) - 1
       if (rank == rank_count - 1) local_sweep%ZE = 2 * rank_count
       local_sweep%NZ = local_sweep%ZE - local_sweep%ZI + 1
       call build_output_partition(request_lower, request_upper, global_bounds, local_sweep, &
                                   iEz, rank, rank_count, partition, status)
       if (status /= OUTPUT_PARTITION_SUCCESS .or. .not. partition%has_data) failures = failures + 1
       local_coverage = 0
       do z = partition%local_lower%z, partition%local_upper%z
          local_coverage(z) = 1
       end do
       call MPI_Reduce(local_coverage, global_coverage, size(local_coverage), MPI_INTEGER, MPI_SUM, root, &
                       MPI_COMM_WORLD, ierr)
       if (rank == root) then
          do z = 0, 2 * rank_count
             if (global_coverage(z) /= 1) failures = failures + 1
          end do
       end if

      ! Root aggregation must reconstruct the same ordered values as a serial output.
      call init_output_collective(collective, rank, rank_count, root, .false., status)
      call prepare_output_partition_publication(collective, participants, owner, partition, &
                                                local_participates, publication_mode, status)
      if (status /= OUTPUT_COLLECTIVE_SUCCESS .or. publication_mode /= OUTPUT_PUBLICATION_ROOT_AGGREGATION) then
         failures = failures + 1
      end if

      local_count = int(partition%local_shape(3))
      allocate(local_values(local_count))
      do i = 1, local_count
         local_values(i) = partition%local_lower%z + i - 1
      end do
      allocate(counts(rank_count), displacements(rank_count), gathered_values(2 * rank_count + 1))
      call MPI_Gather(local_count, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
      if (rank == root) then
         displacements(1) = 0
         do i = 2, rank_count
            displacements(i) = displacements(i - 1) + counts(i - 1)
         end do
         total_count = sum(counts)
      end if
      call MPI_Gatherv(local_values, local_count, MPI_INTEGER, gathered_values, counts, displacements, &
                       MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
      if (rank == root) then
         if (total_count /= 2 * rank_count + 1) failures = failures + 1
         do i = 1, total_count
            if (gathered_values(i) /= i - 1) failures = failures + 1
         end do
      end if
   end if

   call MPI_Allreduce(MPI_IN_PLACE, failures, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
   call MPI_Finalize(ierr)
   if (failures /= 0) error stop failures
end program test_output_mpi
