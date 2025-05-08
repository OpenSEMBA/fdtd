
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module thin wires from Wires paper
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module Wire_bundles_mtln_mod             
#ifdef CompileWithMTLN

   use report
   use fdetypes
   use mtln_solver_mod , mtln_solver_t => mtln_t 
   use mtln_types_mod, only: mtln_t
   use HollandWires
   use wiresHolland_constants
   use ilumina
   implicit none
   
   REAL (KIND=RKIND_wires)           ::  eps0,mu0
   private   

   public InitWires_mtln, AdvanceWiresE_mtln, GetSolverPtr, solveMTLNProblem, reportSimulationEnd
   type(mtln_solver_t), target :: mtln_solver
   integer, dimension(:,:), allocatable :: indexMap

contains


   subroutine InitWires_mtln(sgg,Ex,Ey,Ez, Idxh, Idyh, Idzh, eps00, mu00, mtln_parsed,thereAreMTLNbundles)
      type (SGGFDTDINFO), intent(IN), target    :: sgg 
      REAL (KIND=RKIND), intent(inout), target :: &
         Ex(sgg%Alloc(iEx)%XI : sgg%Alloc(iEx)%XE,  &
            sgg%Alloc(iEx)%YI : sgg%Alloc(iEx)%YE,  &
            sgg%Alloc(iEx)%ZI : sgg%Alloc(iEx)%ZE), &
         Ey(sgg%Alloc(iEy)%XI : sgg%Alloc(iEy)%XE,  &
            sgg%Alloc(iEy)%YI : sgg%Alloc(iEy)%YE,  &
            sgg%Alloc(iEy)%ZI : sgg%Alloc(iEy)%ZE), &
         Ez(sgg%Alloc(iEz)%XI : sgg%Alloc(iEz)%XE,  &
            sgg%Alloc(iEz)%YI : sgg%Alloc(iEz)%YE,  &
            sgg%Alloc(iEz)%ZI : sgg%Alloc(iEz)%ZE)
      real (kind=RKIND), dimension (:), intent(in) :: &
            Idxh(sgg%ALLOC(iEx)%XI : sgg%ALLOC(iEx)%XE),&
            Idyh(sgg%ALLOC(iEy)%YI : sgg%ALLOC(iEy)%YE),&
            Idzh(sgg%ALLOC(iEz)%ZI : sgg%ALLOC(iEz)%ZE)  
   
      real(KIND=RKIND) :: eps00,mu00
   
      type(mtln_t) :: mtln_parsed
      logical :: thereAreMTLNbundles
      type(Thinwires_t), pointer  ::  hwires
#ifdef CompileWithMPI
      integer(kind=4) :: ierr
#endif


      eps0 = eps00
      mu0 = mu00

#ifdef CompileWithMPI
      mtln_solver = mtlnCtor(mtln_parsed, sgg%alloc)
      call mpi_barrier(subcomm_mpi,ierr)
#else
      mtln_solver = mtlnCtor(mtln_parsed)
#endif

      if (mtln_solver%number_of_bundles>=1) then 
           thereAreMTLNbundles=.true.
        else    
           thereAreMTLNbundles=.false.
           return
      endif

      hwires => GetHwires()
      indexMap = mapFieldToCurrentSegments(hwires, mtln_solver%bundles)

! #ifdef CompileWithMPI
!       call mpi_barrier(subcomm_mpi,ierr)
! #endif

      call pointSegmentsToFields()

      call assignLCToExternalConductor()
! #ifdef CompileWithMPI
      ! call mpi_barrier(subcomm_mpi,ierr)
! #endif
      call updateNetworksLineCapacitors()
      call mtln_solver%updatePULTerms()

   contains

      subroutine pointSegmentsToFields()
         integer (kind=4) :: i, j, k, m, n, direction
         do m = 1, mtln_solver%number_of_bundles
            if (mtln_solver%bundles(m)%buildLine) then 
               do n = 1, ubound(mtln_solver%bundles(m)%external_field_segments,1)
                  call readGridIndices(i, j, k, mtln_solver%bundles(m)%external_field_segments(n))
                  direction = mtln_solver%bundles(m)%external_field_segments(n)%direction
                  select case (abs(mtln_solver%bundles(m)%external_field_segments(n)%direction))
                     case(1)                                
                     mtln_solver%bundles(m)%external_field_segments(n)%field => Ex(i, j, k) 
                     case(2)     
                     mtln_solver%bundles(m)%external_field_segments(n)%field => Ey(i, j, k) 
                     case(3)     
                     mtln_solver%bundles(m)%external_field_segments(n)%field => Ez(i, j, k) 
                     end select
               end do
            end if
         end do
      end subroutine

      subroutine assignLCToExternalConductor()
         integer(kind=4) :: m, n
         real (kind=rkind) :: l,c
! #ifdef CompileWithMPI
!          call mpi_barrier(subcomm_mpi,ierr)
! #endif

         do m = 1, mtln_solver%number_of_bundles
            if (mtln_solver%bundles(m)%buildLine) then 
               do n = 1, ubound(mtln_solver%bundles(m)%lpul,1)
                  l = hwires%CurrentSegment(indexMap(m,n))%Lind
                  c = mu0*eps0/l
                  if (mtln_solver%bundles(m)%lpul(n,1,1) == 0.0) then 
                  ! if (mtln_solver%bundles(m)%lpul(n,1,1) == 0.0 .and. mtln_solver%bundles(m)%isPassthrough .eqv. .false.) then 
                     mtln_solver%bundles(m)%lpul(n,1,1) = l 
                     mtln_solver%bundles(m)%cpul(n,1,1) = c 
                     if (mtln_solver%bundles(m)%external_field_segments(n)%has_dielectric) then
                        mtln_solver%bundles(m)%external_field_segments(n)%dielectric%effective_relative_permittivity = computeEffectivePermittivity(m,n,l,c)
                     end if
                  end if
               end do
               mtln_solver%bundles(m)%cpul(ubound(mtln_solver%bundles(m)%cpul,1),1,1) = &
                  mtln_solver%bundles(m)%cpul(ubound(mtln_solver%bundles(m)%cpul,1)-1,1,1)
         end if
         end do
      end subroutine

      function computeEffectivePermittivity(m,n,l0,c0) result(res)
         integer (kind=4) :: i, j, k, direction
         integer (kind=4) :: m,n
         real(kind=rkind) :: dS_inverse
         real(kind=rkind) :: l0, c0
         real(kind=rkind) :: shield_inductance, external_inductance
         real(kind=rkind) :: shield_capacitance, external_capacitance, effective_capacitance
         real(kind=rkind) :: r, r_diel
         real(kind=rkind) :: res

         call readGridIndices(i, j, k, mtln_solver%bundles(m)%external_field_segments(n))      
         direction = mtln_solver%bundles(m)%external_field_segments(n)%direction
         select case (abs(direction))  
         case(1)   
            dS_inverse = (idyh(j)*idzh(k))
         case(2)     
            dS_inverse = (idxh(i)*idzh(k))
         case(3)   
            dS_inverse = (idxh(i)*idyh(j))
         end select

         r = mtln_solver%bundles(m)%external_field_segments(n)%radius
         r_diel = mtln_solver%bundles(m)%external_field_segments(n)%dielectric%radius

         shield_inductance = (0.5*mu0/pi)*log(r_diel/r)*(1-pi*r**2*dS_inverse)
         external_inductance = l0 - shield_inductance

         shield_capacitance = mu0*eps0*mtln_solver%bundles(m)%external_field_segments(n)%dielectric%relative_permittivity/shield_inductance
         external_capacitance = mu0*eps0/external_inductance

         effective_capacitance = shield_capacitance*external_capacitance/(shield_capacitance+external_capacitance)
         res = effective_capacitance/c0

      end function

      subroutine updateNetworksLineCapacitors()
         integer(kind=4) :: m,init, end, sep
         do m = 1, mtln_solver%number_of_bundles
            if (mtln_solver%bundles(m)%buildLine) then 
               init = lbound(mtln_solver%bundles(m)%cpul,1)
               end = ubound(mtln_solver%bundles(m)%cpul,1)
               sep = index(mtln_solver%bundles(m)%name,"_")
               call mtln_solver%network_manager%circuit%modifyLineCapacitorValue(&
                  trim(mtln_solver%bundles(m)%name(sep+1:))//"_1_initial",&
                  mtln_solver%bundles(m)%cpul(init,1,1)*mtln_solver%bundles(m)%step_size(init)*0.5 )

               call mtln_solver%network_manager%circuit%modifyLineCapacitorValue(&
                  trim(mtln_solver%bundles(m)%name(sep+1:))//"_1_end",&
                  mtln_solver%bundles(m)%cpul(end,1,1)*mtln_solver%bundles(m)%step_size(end-1)*0.5)
            end if
         end do   
      end subroutine

      function mapFieldToCurrentSegments(wires, bundles) result (res)
         type(Thinwires_t), pointer  ::  wires
         type(mtl_bundle_t), allocatable, dimension(:) :: bundles
         integer, dimension(:,:), allocatable :: res
         integer :: m, n, nmax, iw
         integer :: i, j, k
         nmax = 0
         do m = 1, mtln_solver%number_of_bundles
            if (ubound(mtln_solver%bundles(m)%lpul,1) > nmax) then 
               nmax = ubound(mtln_solver%bundles(m)%lpul,1)
            end if 
         end do
         allocate(res(mtln_solver%number_of_bundles,nmax))
         res(:,:) = 0
         do m = 1, mtln_solver%number_of_bundles
            if (mtln_solver%bundles(m)%buildLine) then 
               do n = 1, ubound(mtln_solver%bundles(m)%lpul,1)
                  call readGridIndices(i, j, k, mtln_solver%bundles(m)%external_field_segments(n))                          
                  do iw = 1, wires%NumCurrentSegments
                     if ((i == wires%CurrentSegment(iw)%i) .and. &
                        (j == wires%CurrentSegment(iw)%j) .and. &
                        (k == wires%CurrentSegment(iw)%k)) then
                           res(m,n) = iw
                     end if
                  end do
               end do
            end if
         end do
      end function
   endsubroutine InitWires_mtln

   subroutine AdvanceWiresE_mtln(sgg,Idxh, Idyh, Idzh, eps00,mu00)  
      type (SGGFDTDINFO), intent(IN), target    :: sgg      
      real (kind=RKIND), dimension (:), intent(in) :: &
         Idxh(sgg%ALLOC(iEx)%XI : sgg%ALLOC(iEx)%XE),&
         Idyh(sgg%ALLOC(iEy)%YI : sgg%ALLOC(iEy)%YE),&
         Idzh(sgg%ALLOC(iEz)%ZI : sgg%ALLOC(iEz)%ZE)  
      real(KIND=RKIND) :: cte,eps00,mu00
      integer (kind=4) :: m, n
      REAL (KIND=RKIND),pointer:: punt
      type(Thinwires_t), pointer  ::  hwires
#ifdef CompileWithMPI      
      integer(kind=4) :: ierr, rank
#endif

      eps0 = eps00 
      mu0 = mu00
      
      hwires => GetHwires()
      do m = 1, mtln_solver%number_of_bundles
         if (mtln_solver%bundles(m)%buildLine) then 
            do n = 1, ubound(mtln_solver%bundles(m)%external_field_segments,1)
               punt => mtln_solver%bundles(m)%external_field_segments(n)%field
               punt = real(punt, kind=rkind_wires) - computeFieldFromCurrent(m,n)
               hwires%CurrentSegment(indexMap(m,n))%CurrentPast = getOrientedCurrent(m,n)
            end do
         end if
      end do

#ifdef CompileWithMPI      
!       call mpi_barrier(subcomm_mpi,ierr)
! if (mtln_solver%number_of_bundles /= 0) then 
      ! call MPI_COMM_RANK(SUBCOMM_MPI, rank, ierr)
      ! write(*,*) 'rank ', rank, ' mtln_solver%step()'
      ! call mtln_solver%step()
! else
!    call mtln_solver%advanceTime()
! end if
#endif
   call mtln_solver%step()


   contains

      function getOrientedCurrent(m, n) result(res)
         integer(kind=4), intent(in) :: m, n
         real(kind=rkind) :: res
         real(kind=rkind) :: curr
         integer (kind=4) :: direction, i
         direction = mtln_solver%bundles(m)%external_field_segments(n)%direction
         if (mtln_solver%bundles(m)%isPassthrough) then 
            curr = 0
            do i = 2, 1 + mtln_solver%bundles(m)%conductors_in_level(2)
               curr = curr + mtln_solver%bundles(m)%i(i, n)
            end do
            mtln_solver%bundles(m)%i(1, n) = curr
         end if
         
         res = mtln_solver%bundles(m)%i(1, n) * sign(1.0, real(direction))
      end function

      function computeFieldFromCurrent(m, n) result(res)
         integer(kind=4), intent(in) :: m, n
         real(kind=rkind) :: dS_inverse, factor
         real(kind=rkind) :: res
         integer (kind=4) :: i, j, k, direction
         direction = mtln_solver%bundles(m)%external_field_segments(n)%direction
         call readGridIndices(i, j, k, mtln_solver%bundles(m)%external_field_segments(n))      
         select case (abs(direction))  
         case(1)   
            dS_inverse = (idyh(j)*idzh(k))
         case(2)     
            dS_inverse = (idxh(i)*idzh(k))
         case(3)   
            dS_inverse = (idxh(i)*idyh(j))
         end select
         factor = (sgg%dt / (eps0)) * dS_inverse
         res = factor * getOrientedCurrent(m, n)
      end function

   end subroutine AdvanceWiresE_mtln


   subroutine readGridIndices(i, j, k, field_segment)
      type(external_field_segment_t), intent(in) :: field_segment
      integer, intent(inout) :: i, j, k
      i = field_segment%position(1)
      j = field_segment%position(2)
      k = field_segment%position(3)
   end subroutine

   function GetSolverPtr() result(res)
      type(mtln_solver_t), pointer :: res
      res => mtln_solver
   end function

   subroutine solveMTLNProblem(mtln_parsed)
      type(mtln_t) :: mtln_parsed
      mtln_solver = mtlnCtor(mtln_parsed)
      call mtln_solver%run()
   end subroutine

   subroutine reportSimulationEnd(layoutnumber)
      character (len=bufsize) :: dubuf
      integer (kind=4), intent(in) ::  layoutnumber
      write(dubuf, *) 'MTLN simulation finished. Init flusing probe data to output files'
      call print11(layoutnumber,dubuf)

   end subroutine


#endif

end module Wire_bundles_mtln_mod
