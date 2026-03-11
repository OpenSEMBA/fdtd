module xdmf_h5
#ifdef CompileWithHDF

    use fdetypes
    use HDF5
    implicit none
 
    integer(HID_T) :: file_id ! File identifier
    integer(HID_T) :: dset_id ! Dataset identifier
    integer(HID_T) :: dspace_id, slice2D_id ! Dataspace identifier
    integer(HSIZE_T), ALLOCATABLE, dimension(:) :: DATA_dims ! Dataset dimensions
    integer(HSIZE_T), ALLOCATABLE, dimension(:) :: offset
    integer(HSIZE_T), ALLOCATABLE, dimension(:) :: valor3d_dims ! slice dimensions
   
    !
    private
    public openh5file,writeh5file,closeh5file,createh5filefromsinglebin
    
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine openh5file(filename,finalstep,minXabs,maxXabs, minYabs,maxYabs, minZabs,maxZabs)
      
        integer :: error ! Error flag
        character(len=BUFSIZE) :: filename ! File name
        character(len=BUFSIZE) :: dsetname ! Dataset name
        !
        integer(kind=4) :: minXabs, maxXabs, minYabs, maxYabs, minZabs, maxZabs,finalstep

        integer :: rank ! Dataset rank
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
        
        dsetname = 'data'
        call h5open_f (error)
        call h5fcreate_f (trim(adjustl(filename))//'.h5', H5F_ACC_TRUNC_F, file_id, error)
        call h5screate_simple_f (rank, DATA_dims, dspace_id, error)
        call h5screate_simple_f (rank, valor3d_dims, slice2D_id, error)
#ifdef CompileWithReal8
        call h5dcreate_f (file_id, trim(adjustl(dsetname)), H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
#else
        call h5dcreate_f (file_id, trim(adjustl(dsetname)), H5T_NATIVE_REAL , dspace_id, dset_id, error)
#endif

!xdmf part
        open(18, FILE=trim(adjustl(filename))//'.xdmf', FORM='formatted')
        write(18,*) '<Xdmf>'
        write(18,*) '<Domain>'
        write(18,*) '<Grid Name="GridTime" GridType="Collection" CollectionType="Temporal">'


   end subroutine openh5file
   
   subroutine writeh5file(filename,valor3d,indi,attindi,minXabs,maxXabs, minYabs,maxYabs, minZabs,maxZabs, &
                          linez_minZabs_primero,liney_minYabs_primero,linex_minXabs_primero, &
                          dz_minZabs,dy_minYabs,dx_minXabs,&
                          minZabs_primero,minYabs_primero,minXabs_primero,finalstep,vtkindex)
        real(  KINd=RKIND_tiempo) :: attindi
        character(len=BUFSIZE) :: filename
        logical :: vtkindex
      !
        character(len=BUFSIZE) :: dsetname ! Dataset name
        integer(kind=4) :: indi
        real(kind=RKIND), dimension(:, :, :, :) :: valor3d
        integer :: error ! Error flag
        character(len=BUFSIZE) :: charc
        integer(kind=4) :: minXabs, maxXabs, minYabs, maxYabs, minZabs, maxZabs, &
                            minZabs_primero,minYabs_primero,minXabs_primero,finalstep
        real(kind=RKIND) :: linez_minZabs_primero,liney_minYabs_primero,linex_minXabs_primero, &
                             dz_minZabs,dy_minYabs,dx_minXabs                 
        
        
        offset (1) = 0
        offset (2) = 0
        offset (3) = 0
        offset (4) = indi - 1
        !
        call h5sselect_hyperslab_f (dspace_id, H5S_SELECT_SET_F, offset, valor3d_dims, error)
#ifdef CompileWithReal8
        call h5dwrite_f (dset_id, H5T_NATIVE_DOUBLE, valor3d, valor3d_dims, error, slice2D_id, &
        & dspace_id)
#elif CompileWithReal16
        call h5dwrite_f (dset_id, H5T_NATIVE_LDOUBLE, valor3d, valor3d_dims, error, slice2D_id, &
        & dspace_id)
#else
        call h5dwrite_f (dset_id, H5T_NATIVE_REAL, valor3d, valor3d_dims, error, slice2D_id, &
        & dspace_id)
#endif

        !el .xdmf como usualmente
        !HDF5 transposes matrices
        write(charc,'(e19.9e3)') attindi  !'(i9)') indi !
        dsetname = 'data'         
        DATA_dims(1) = maxXabs - minXabs + 1
        DATA_dims(2) = maxYabs - minYabs + 1
        DATA_dims(3) = maxZabs - minZabs + 1
        write(18, '(a)') '<Grid Name="IntGrid" GridType="Uniform"  CollectionType="Spatial">>'
        write(18, '(a)') '<Time Value="' // trim (adjustl(charc)) // '" />'
        write(18, '(a,3i5,a)') '<Topology TopologyType="3DCoRectMesh" Dimensions="', DATA_dims(3), &
        & DATA_dims(2), DATA_dims(1), '">'
        write(18, '(a)') '</Topology>'
        write(18, '(a)') '<Geometry Type="ORIGIN_DXDYDZ">'
        write(18, '(a)') '<DataItem Format="XML" Dimensions="3">'
        if (vtkindex) then
           write(18, *)  minZabs_primero,minYabs_primero,minXabs_primero !ojo solo funciona bien el escalado en mallados uniformes salva'oct14
        else
           write(18, *) linez_minZabs_primero,liney_minYabs_primero,linex_minXabs_primero !ojo solo funciona bien el escalado en mallados uniformes salva'oct14
        endif
        write(18, '(a)') '</DataItem>'
        write(18, '(a)') '<DataItem Format="XML" Dimensions="3">'
        write(18, *) dz_minZabs,dy_minYabs,dx_minXabs
        write(18, '(a)') '</DataItem>'
        write(18, '(a)') '</Geometry>'
        write(18, '(a)') '<Attribute Name="IntValues" Center="Node">'
        write(18, '(a,4i5,a)') '<DataItem ItemType="HyperSlab" Dimensions="', 1, DATA_dims(3), &
        & DATA_dims(2), DATA_dims(1), '" Format="XML">'
        write(18, '(a)') '<DataItem Dimensions="3 4" Format="XML">'
        write(18, '(4i5)') offset(4), 0, 0, 0
        write(18, '(4i5)') 1, 1, 1, 1
        write(18, '(4i5)') 1, DATA_dims(3), DATA_dims(2), DATA_dims(1)
        write(18, '(a)') '</DataItem>'
        write(18, '(a,4i5,a)') '<DataItem Format="HDF" NumberType="Float" Precision="4" Dimensions="',&
        &  finalstep, DATA_dims(3), DATA_dims(2), DATA_dims(1), '">'
        write(18, '(a)') trim (adjustl(filename)) // '.h5:/' // trim (adjustl(dsetname))
        write(18, '(a)') '</DataItem>'
        write(18, '(a)') '</DataItem>'
        write(18, '(a)') '</Attribute>'
        write(18, '(a)') '</Grid>'

                 
   end subroutine writeh5file
   
   
   subroutine closeh5file(finalstep,att)
      !
        integer :: rank ! Dataset rank
        real(  KINd=RKIND_tiempo), dimension(:) :: att
        integer :: error ! Error flag
        integer(kind=4) :: finalstep
        character(len=BUFSIZE) :: dsetname ! Dataset name

        
        DEALLOCATE(DATA_dims,valor3d_dims,offset)
        
        !timedata
        call h5dclose_f (dset_id, error)

        dsetname='Time'
        rank = 1
        ALLOCATE(DATA_dims(rank))
        data_dims(1) = finalstep

        call h5screate_simple_f (rank, DATA_dims, dspace_id, error)
#ifdef CompileWithReal8
        call h5dcreate_f (file_id, trim(adjustl(dsetname)), H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
        call h5dwrite_f (dset_id, H5T_NATIVE_DOUBLE, att, DATA_dims, error)
#elif CompileWithReal16
        call h5dcreate_f (file_id, trim(adjustl(dsetname)), H5T_NATIVE_LDOUBLE, dspace_id, dset_id, error)
        call h5dwrite_f (dset_id, H5T_NATIVE_REAL, att, DATA_dims, error)
#else
        call h5dcreate_f (file_id, trim(adjustl(dsetname)), H5T_NATIVE_REAL, dspace_id, dset_id, error)
        call h5dwrite_f (dset_id, H5T_NATIVE_REAL, att, DATA_dims, error)
#endif
        cALL h5dclose_f (dset_id, error)

        call h5sclose_f (slice2D_id, error)
        call h5sclose_f (dspace_id, error)
   
                        
                        
        call h5fclose_f (file_id, error)
        call h5close_f (error)
        !
                        !
        write(18, '(a)') '</Grid>'
        write(18, '(a)') '</Domain>'
        write(18, '(a)') '</Xdmf>'
        CLOSE (18)
        !
        DEALLOCATE(DATA_dims)
                        
                        
   end subroutine closeh5file

   subroutine createh5filefromsinglebin(filename,vtkindex)
        integer(kind=4) :: myunit,fieldob,pasadas,pasadastotales
        character(len=BUFSIZE) :: filename,fichin ! File name
        real(  KINd=RKIND_tiempo), ALLOCATABLE, dimension(:) :: att
        real(kind=RKIND), ALLOCATABLE, dimension(:, :, :, :) :: valor3d !para sondas Volumic
        logical :: vtkindex,SGGObservationiiTimeDomain
        integer(kind=4) :: minXabs, maxXabs, minYabs, maxYabs, minZabs, maxZabs, &
                            minZabs_primero,minYabs_primero,minXabs_primero,finalstep,indi,i1,j1,k1
        real(kind=RKIND) :: linez_minZabs_primero,liney_minYabs_primero,linex_minXabs_primero, &
                             dz_minZabs,dy_minYabs,dx_minXabs    
        character(len=BUFSIZE) :: dubuf

        filename=filename(1:index(filename,'.h5bin')-1); filename=trim(adjustl(filename))
   
        open(newunit=myunit,file=trim(adjustl(filename))//'.h5bin',form='unformatted')
        read (myunit)            finalstep,minXabs, maxXabs, minYabs, maxYabs, minZabs, maxZabs,fieldob, &
                                 SGGObservationiiTimeDomain,pasadastotales
        
       allocate(valor3d(minXabs:maxXabs, minYabs:maxYabs, minZabs:maxZabs, 1))
        allocate (att(1:finalstep))
        
        buclepasadas: do pasadas=1,pasadastotales                                     
            if (SGGObservationiiTimeDomain) then         
                if (pasadas==1) then   
                    fichin = trim (adjustl(filename))//'_time'
                else
                    print *,'Buggy error in valor3d. '
                    stop
                endif  
            else 
                if (pasadas==1) then   
                    fichin = trim (adjustl(filename))//'_mod'
                elseif (pasadas==2) then  
                    fichin = trim (adjustl(filename))//'_phase'
                else
                    print *,'Buggy error in valor3d. '
                    stop
                endif
            endif
            
            
            if (.not.(((fieldob == iMEC).or.(fieldob ==iMHC)).and.(pasadas ==2))) then ! no tiene sentido esccribir la fase
                call openh5file(fichin,finalstep,minXabs,maxXabs, minYabs,maxYabs, minZabs,maxZabs)
            endif
 
            valor3d = 0.0_RKIND
            att=0.0_RKIND
            do indi = 1, finalstep
                read(myunit) minZabs_primero,minYabs_primero,minXabs_primero      
                read(myunit) linez_minZabs_primero,liney_minYabs_primero,linex_minXabs_primero
                read(myunit) dz_minZabs,dy_minYabs,dx_minXabs                    
                read (myunit) att(indi)
                write(dubuf,*)  ' ----> .xdmf file ',att(indi),'(',indi,'/',finalstep,')'
                print *,trim(adjustl(dubuf))
                do k1 = minzabs, maxzabs
                    do j1 = minyabs, maxyabs
                        read (myunit) (valor3d(i1, j1, k1, 1), i1=minxabs, maxxabs)
                    end do
                end do
                
                if (.not.(((fieldob == iMEC).or.(fieldob ==iMHC)).and.(pasadas ==2))) then ! no tiene sentido esccribir la fase
                   call writeh5file(fichin,valor3d,indi,att(indi),minXabs,maxXabs, minYabs,maxYabs, minZabs,maxZabs, &
                                    linez_minZabs_primero,liney_minYabs_primero,linex_minXabs_primero, &
                                    dz_minZabs,dy_minYabs,dx_minXabs,&
                                    minZabs_primero,minYabs_primero,minXabs_primero,finalstep,vtkindex)
                endif
            end do
            
            if (.not.(((fieldob == iMEC).or.(fieldob ==iMHC)).and.(pasadas ==2))) then ! no tiene sentido esccribir la fase
                call closeh5file(finalstep,att)
            endif
        end do buclepasadas
    
        close(myunit)
        
        deallocate(valor3d)
        deallocate(ATT)
             
   
   end subroutine createh5filefromsinglebin


#endif               

end module xdmf_h5 