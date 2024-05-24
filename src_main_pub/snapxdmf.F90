MODULE snapxdmf
   !

#ifdef CompileWithHDF
   USE fdetypes
   USE HDF5
   !
   !
   IMPLICIT NONE
   !
   PRIVATE
   PUBLIC WRITE_XDMFSNAP
CONTAINS

   !snaps

   SUBROUTINE write_xdmfsnap(ninstant,filename,minXabs,maxXabs,minYabs,maxYabs,minZabs,maxZabs,valor3D)

      !------------------------>



      INTEGER (KIND=4) :: minXabs, maxXabs, minYabs, maxYabs, minZabs, maxZabs
      real (kind=4), dimension(minXabs:maxXabs,minYabs:maxYabs,minZabs:maxZabs,1:1)  ::  valor3D
      CHARACTER (LEN=BUFSIZE) :: filename ! File name
      CHARACTER (LEN=BUFSIZE) :: dsetname ! Dataset name
      !
      INTEGER (HID_T) :: file_id  ! File identifier
      INTEGER (HID_T) :: dset_id ! Dataset identifier
      INTEGER (HID_T) :: dspace_id, slice2D_id ! Dataspace identifier
      !
      INTEGER :: error ! Error flag
      INTEGER :: rank ! Dataset rank
      INTEGER (HSIZE_T), ALLOCATABLE, DIMENSION (:) :: DATA_dims ! Dataset dimensions
      INTEGER (HSIZE_T), ALLOCATABLE, DIMENSION (:) :: offset
      INTEGER (HSIZE_T), ALLOCATABLE, DIMENSION (:) :: valor3d_dims ! slice dimensions
      !
      CHARACTER (LEN=BUFSIZE) :: charc
      INTEGER (KIND=4), INTENT (IN) :: ninstant
      !
      !
      !
      INTEGER (KIND=4) :: indi
      REAL (KIND=4), ALLOCATABLE, DIMENSION (:) :: att
      !
      INTEGER (KIND=4) :: finalstep
      !
      finalstep=1
      !
      allocate (att(1:finalstep))
      !
      rank = 4
      ALLOCATE(DATA_dims(1:RANK),valor3d_dims(1:RANK),offset(1:RANK))
      !
      DATA_dims (1) = maxXabs - minXabs + 1
      DATA_dims (2) = maxYabs - minYabs + 1
      DATA_dims (3) = maxZabs - minZabs + 1
      DATA_dims (4) = finalstep
      !
      valor3d_dims (1) = DATA_dims (1)
      valor3d_dims (2) = DATA_dims (2)
      valor3d_dims (3) = DATA_dims (3)
      valor3d_dims (4) = 1
      !
      dsetname = 'data'

      CALL h5open_f (error)
      CALL h5fcreate_f (trim(adjustl(filename))//'.h5', H5F_ACC_TRUNC_F, file_id, error)
      CALL h5screate_simple_f (rank, DATA_dims, dspace_id, error)
      CALL h5screate_simple_f (rank, valor3d_dims, slice2D_id, error)
      CALL h5dcreate_f (file_id, trim(adjustl(dsetname)), H5T_NATIVE_REAL , dspace_id, dset_id, error)
      OPEN (18, FILE=trim(adjustl(filename))//'.xdmf', FORM='formatted')
      WRITE (18,*) '<Xdmf>'
      WRITE (18,*) '<Domain>'
      WRITE (18,*) '<Grid Name="GridTime" GridType="Collection" CollectionType="Temporal">'

      !HDF5 transposes the dimensions !be careful
      indi=finalstep
      att(indi)=ninstant
      WRITE (charc, '(i10)') indi
      offset (1) = 0
      offset (2) = 0
      offset (3) = 0
      offset (4) = 0
      !
      CALL h5sselect_hyperslab_f (dspace_id, H5S_SELECT_SET_F, offset, valor3d_dims, error)
      CALL h5dwrite_f (dset_id, H5T_NATIVE_REAL, valor3d, valor3d_dims, error, slice2D_id,dspace_id)

      !
      !HDF5 transposes matrices
      WRITE (18, '(a)') '<Grid Name="IntGrid" GridType="Uniform"  CollectionType="Spatial">>'
      WRITE (18, '(a)') '<Time Value="' // trim (adjustl(charc)) // '" />'
      WRITE (18, '(a,3i5,a)') '<Topology TopologyType="3DCoRectMesh" Dimensions="', DATA_dims (3), &
      & DATA_dims (2), DATA_dims (1), '">'
      WRITE (18, '(a)') '</Topology>'
      WRITE (18, '(a)') '<Geometry Type="ORIGIN_DXDYDZ">'
      WRITE (18, '(a)') '<DataItem Format="XML" Dimensions="3">'
      WRITE (18, '(4f10.1)') 1.0_RKIND *minZabs, 1.0_RKIND *minYabs, 1.0_RKIND *minXabs
      WRITE (18, '(a)') '</DataItem>'
      WRITE (18, '(a)') '<DataItem Format="XML" Dimensions="3">'
      WRITE (18, '(a)') '1.0 1.0 1.0'
      WRITE (18, '(a)') '</DataItem>'
      WRITE (18, '(a)') '</Geometry>'
      WRITE (18, '(a)') '<Attribute Name="IntValues" Center="Node">'
      WRITE (18, '(a,4i5,a)') '<DataItem ItemType="HyperSlab" Dimensions="', 1, DATA_dims (3), &
      & DATA_dims (2), DATA_dims (1), '" Format="XML">'
      WRITE (18, '(a)') '<DataItem Dimensions="3 4" Format="XML">'
      WRITE (18, '(4i5)') offset (4), 0, 0, 0
      WRITE (18, '(4i5)') 1, 1, 1, 1
      WRITE (18, '(4i5)') 1, DATA_dims (3), DATA_dims (2), DATA_dims (1)
      WRITE (18, '(a)') '</DataItem>'
      WRITE (18, '(a,4i5,a)') '<DataItem Format="HDF" NumberType="Float" Precision="4" Dimensions="',&
      &  finalstep, DATA_dims (3), DATA_dims (2), DATA_dims (1), '">'
      WRITE (18, '(a)') trim (adjustl(filename)) // '.h5:/' // trim (adjustl(dsetname))
      WRITE (18, '(a)') '</DataItem>'
      WRITE (18, '(a)') '</DataItem>'
      WRITE (18, '(a)') '</Attribute>'
      WRITE (18, '(a)') '</Grid>'

      !timedata
      CALL h5dclose_f (dset_id, error)
      CALL h5sclose_f (slice2D_id, error)
      CALL h5sclose_f (dspace_id, error)

      DEALLOCATE(DATA_dims,valor3d_dims,offset)

      dsetname='Time'
      rank = 1
      ALLOCATE(DATA_dims(rank))
      data_dims(1) = finalstep

      CALL h5screate_simple_f (rank, DATA_dims, dspace_id, error)
      CALL h5dcreate_f (file_id, trim(adjustl(dsetname)), H5T_NATIVE_REAL, dspace_id, dset_id, error)
      CALL h5dwrite_f (dset_id, H5T_NATIVE_REAL, att, DATA_dims, error)
      !
      cALL h5dclose_f (dset_id, error)
      CALL h5sclose_f (dspace_id, error)
      !
      WRITE (18, '(a)') '</Grid>'
      WRITE (18, '(a)') '</Domain>'
      WRITE (18, '(a)') '</Xdmf>'


      CALL h5fclose_f (file_id, error)
      CALL h5close_f (error)

      CLOSE (18)

      DEALLOCATE (ATT)

      RETURN
   END SUBROUTINE write_xdmfsnap
#endif


END MODULE snapxdmf
!
!
