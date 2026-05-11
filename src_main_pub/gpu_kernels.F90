!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  GPU KERNELS MODULE - OpenACC accelerated YEE kernel subroutines
!  Standalone module for GPU execution with NVHPC compiler
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module gpu_kernels_m

   use FDETYPES_m
   use Report_m

   implicit none

   ! GPU state: holds device pointers and tracking
   type gpu_state_t
      real(kind=rkind), pointer, dimension(:,:,:), contiguous :: Ex, Ey, Ez, Hx, Hy, Hz
      real(kind=rkind), pointer, dimension(:), contiguous :: Idxe, Idye, Idze, Idxh, Idyh, Idzh, dxe, dye, dze, dxh, dyh, dzh
      integer(kind=integersizeofmediamatrices), pointer, dimension(:,:,:), contiguous :: sggMiEx, sggMiEy, sggMiEz, sggMiHx, sggMiHy, sggMiHz
      real(kind=rkind), pointer, dimension(:), contiguous :: g1, g2, gm1, gm2
      logical :: initialized = .false.
   end type

contains

   !--------------------------------------------------------------------------------
   ! Initialize GPU data regions - called once at simulation start
   !--------------------------------------------------------------------------------
   subroutine gpu_init(this, Ex, Ey, Ez, Hx, Hy, Hz, sggMiEx, sggMiEy, sggMiEz, &
                       sggMiHx, sggMiHy, sggMiHz, g1, g2, gm1, gm2, &
                       Idxe_in, Idye_in, Idze_in, Idxh_in, Idyh_in, Idzh_in, dxe_in, dye_in, dze_in, dxh_in, dyh_in, dzh_in)
      class(gpu_state_t), intent(inout) :: this
      real(kind=rkind), dimension(:,:,:), pointer, contiguous, intent(in) :: Ex, Ey, Ez, Hx, Hy, Hz
      integer(kind=integersizeofmediamatrices), dimension(:,:,:), pointer, contiguous, intent(in) :: sggMiEx, sggMiEy, sggMiEz, sggMiHx, sggMiHy, sggMiHz
      real(kind=rkind), dimension(:), pointer, intent(in) :: g1, g2, gm1, gm2
      real(kind=rkind), dimension(:), pointer, intent(in) :: Idxe_in, Idye_in, Idze_in, Idxh_in, Idyh_in, Idzh_in, dxe_in, dye_in, dze_in, dxh_in, dyh_in, dzh_in

#ifdef SEMBA_FDTD_ENABLE_ACC
   ! Data movement is handled by per-kernel implicit copies while debugging runtime stability.
#endif

      this%Ex => Ex
      this%Ey => Ey
      this%Ez => Ez
      this%Hx => Hx
      this%Hy => Hy
      this%Hz => Hz
      this%sggMiEx => sggMiEx
      this%sggMiEy => sggMiEy
      this%sggMiEz => sggMiEz
      this%sggMiHx => sggMiHx
      this%sggMiHy => sggMiHy
      this%sggMiHz => sggMiHz
      this%g1 => g1
      this%g2 => g2
      this%gm1 => gm1
      this%gm2 => gm2
      this%Idxe => Idxe_in
      this%Idye => Idye_in
      this%Idze => Idze_in
      this%Idxh => Idxh_in
      this%Idyh => Idyh_in
      this%Idzh => Idzh_in
      this%dxe => dxe_in
      this%dye => dye_in
      this%dze => dze_in
      this%dxh => dxh_in
      this%dyh => dyh_in
      this%dzh => dzh_in
      this%initialized = .true.

   end subroutine gpu_init

   !--------------------------------------------------------------------------------
   ! Destroy GPU data regions - called once at simulation end
   !--------------------------------------------------------------------------------
   subroutine gpu_destroy(this)
      class(gpu_state_t), intent(inout) :: this

#ifdef SEMBA_FDTD_ENABLE_ACC
   ! No persistent data region to tear down in debug mode.
#endif

      nullify(this%Ex)
      nullify(this%Ey)
      nullify(this%Ez)
      nullify(this%Hx)
      nullify(this%Hy)
      nullify(this%Hz)
      nullify(this%sggMiEx)
      nullify(this%sggMiEy)
      nullify(this%sggMiEz)
      nullify(this%sggMiHx)
      nullify(this%sggMiHy)
      nullify(this%sggMiHz)
      nullify(this%g1)
      nullify(this%g2)
      nullify(this%gm1)
      nullify(this%gm2)
      nullify(this%Idxe)
      nullify(this%Idye)
      nullify(this%Idze)
      nullify(this%Idxh)
      nullify(this%Idyh)
      nullify(this%Idzh)
      nullify(this%dxe)
      nullify(this%dye)
      nullify(this%dze)
      nullify(this%dxh)
      nullify(this%dyh)
      nullify(this%dzh)
      this%initialized = .false.

   end subroutine gpu_destroy

   !--------------------------------------------------------------------------------
   ! Advance Ex field - GPU accelerated YEE kernel (PENDING)
   !
   ! NOTE: This kernel is currently blocked at runtime due to NVHPC OpenACC compiler issues.
   ! Even a no-op kernel fails to launch, suggesting an incompatibility between NVHPC 26.3
   ! OpenACC runtime and the kernel's array bounds model. The kernel code is preserved below
   ! for future work with newer compiler versions or alternative GPU strategies (CUF, direct CUDA).
   !
   ! Call path: timestepping.F90::advanceEx() [currently disabled, falls through to CPU]
   !
   ! TODO: Revisit with NVHPC >= 24.7 or alternative GPU approaches
   !--------------------------------------------------------------------------------
   subroutine gpu_advanceEx(Ex, Hy, Hz, sggMiEx, Idyh, Idzh, g1, g2, nx, ny, nz)
      ! Preserved for future work; currently unused due to runtime blocker
      integer(kind=4), intent(in) :: nx, ny, nz
      real(kind=rkind), dimension(:,:,:), intent(inout) :: Ex
      real(kind=rkind), dimension(:,:,:), intent(in) :: Hy
      real(kind=rkind), dimension(:,:,:), intent(in) :: Hz
      integer(kind=integersizeofmediamatrices), dimension(:,:,:), intent(in) :: sggMiEx
      real(kind=rkind), dimension(:), intent(in) :: Idyh
      real(kind=rkind), dimension(:), intent(in) :: Idzh
      real(kind=rkind), dimension(:), intent(in) :: g1, g2

      real(kind=rkind) :: Idzhk, Idyhj
      integer(kind=4) :: i, j, k
      integer(kind=integersizeofmediamatrices) :: medio

#ifdef SEMBA_FDTD_ENABLE_ACC
   !$ACC parallel loop gang
#endif
      Do k=1,nz
         Do j=1,ny
            Do i=1,nx
               Idzhk=Idzh(k)
               Idyhj=Idyh(j)
               medio =sggMiEx(i,j,k)
               Ex(i,j,k)=g1(medio)*Ex(i,j,k)+g2(medio)* &
               ((Hz(i,j,k)-Hz(i,j-1,k))*Idyhj-(Hy(i,j,k)-Hy(i,j,k-1))*Idzhk)
            End do
         End do
      End do
#ifdef SEMBA_FDTD_ENABLE_ACC
   !$acc wait
#endif

   end subroutine gpu_advanceEx

   !--------------------------------------------------------------------------------
   ! Advance Ey field - GPU accelerated YEE kernel
   !--------------------------------------------------------------------------------
   subroutine gpu_advanceEy(this, bounds)
      class(gpu_state_t), intent(inout) :: this
      type(bounds_t), intent(in) :: bounds

      real(kind=rkind), dimension(:,:,:), pointer, contiguous :: Ey
      real(kind=rkind), dimension(:,:,:), pointer, contiguous :: Hz
      real(kind=rkind), dimension(:,:,:), pointer, contiguous :: Hx
      integer(kind=integersizeofmediamatrices), dimension(:,:,:), pointer, contiguous :: sggMiEy
      real(kind=rkind), dimension(:), pointer :: Idzh
      real(kind=rkind), dimension(:), pointer :: Idxh
      real(kind=rkind), dimension(:), pointer :: g1
      real(kind=rkind), dimension(:), pointer :: g2

      real(kind=rkind) :: Idzhk
      integer(kind=4) :: i, j, k
      integer(kind=integersizeofmediamatrices) :: medio

      Ey(0:bounds%Ey%NX-1,0:bounds%Ey%NY-1,0:bounds%Ey%NZ-1) => this%Ey
      Hz(0:bounds%Hz%NX-1,0:bounds%Hz%NY-1,0:bounds%Hz%NZ-1) => this%Hz
      Hx(0:bounds%Hx%NX-1,0:bounds%Hx%NY-1,0:bounds%Hx%NZ-1) => this%Hx
      sggMiEy(0:bounds%sggMiEy%NX-1,0:bounds%sggMiEy%NY-1,0:bounds%sggMiEy%NZ-1) => this%sggMiEy
      Idzh(0:bounds%dzh%NZ-1) => this%Idzh
      Idxh(0:bounds%dxh%NX-1) => this%Idxh
      g1 => this%g1
      g2 => this%g2

#ifdef SEMBA_FDTD_ENABLE_ACC
   !$ACC data copy(Ey) copyin(Hz,Hx,sggMiEy,Idzh,Idxh,g1,g2)
   !$ACC parallel loop gang vector collapse(3) independent private(i,j,k,medio,Idzhk)
#endif
      Do k=1,bounds%sweepEy%NZ
         Do j=1,bounds%sweepEy%NY
            Do i=1,bounds%sweepEy%NX
               Idzhk=Idzh(k)
               medio =sggMiEy(i,j,k)
               Ey(i,j,k)=g1(medio)*Ey(i,j,k)+g2(medio)*((Hx(i,j,k)-Hx(i,j,k-1))*Idzhk-(Hz(i,j,k)-Hz(i-1,j,k))*Idxh(i))
            End do
         End do
      End do
#ifdef SEMBA_FDTD_ENABLE_ACC
   !$acc wait
   !$acc end data
#endif

   end subroutine gpu_advanceEy

   !--------------------------------------------------------------------------------
   ! Advance Ez field - GPU accelerated YEE kernel
   !--------------------------------------------------------------------------------
   subroutine gpu_advanceEz(this, bounds)
      class(gpu_state_t), intent(inout) :: this
      type(bounds_t), intent(in) :: bounds

      real(kind=rkind), dimension(:,:,:), pointer, contiguous :: Ez
      real(kind=rkind), dimension(:,:,:), pointer, contiguous :: Hx
      real(kind=rkind), dimension(:,:,:), pointer, contiguous :: Hy
      integer(kind=integersizeofmediamatrices), dimension(:,:,:), pointer, contiguous :: sggMiEz
      real(kind=rkind), dimension(:), pointer :: Idyh
      real(kind=rkind), dimension(:), pointer :: Idxh
      real(kind=rkind), dimension(:), pointer :: g1
      real(kind=rkind), dimension(:), pointer :: g2

      real(kind = RKIND) :: Idyhj
      integer(kind = 4) :: i, j, k
      integer(kind = INTEGERSIZEOFMEDIAMATRICES) :: medio

      Ez(0:bounds%Ez%NX-1,0:bounds%Ez%NY-1,0:bounds%Ez%NZ-1) => this%Ez
      Hx(0:bounds%Hx%NX-1,0:bounds%Hx%NY-1,0:bounds%Hx%NZ-1) => this%Hx
      Hy(0:bounds%Hy%NX-1,0:bounds%Hy%NY-1,0:bounds%Hy%NZ-1) => this%Hy
      sggMiEz(0:bounds%sggMiEz%NX-1,0:bounds%sggMiEz%NY-1,0:bounds%sggMiEz%NZ-1) => this%sggMiEz
      Idyh(0:bounds%dyh%NY-1) => this%Idyh
      Idxh(0:bounds%dxh%NX-1) => this%Idxh
      g1 => this%g1
      g2 => this%g2

#ifdef SEMBA_FDTD_ENABLE_ACC
   !$ACC data copy(Ez) copyin(Hx,Hy,sggMiEz,Idyh,Idxh,g1,g2)
   !$ACC parallel loop gang vector collapse(3) independent private(i,j,k,medio,Idyhj)
#endif
      Do k=1,bounds%sweepEz%NZ
         Do j=1,bounds%sweepEz%NY
            Do i=1,bounds%sweepEz%NX
               Idyhj=Idyh(j)
               medio =sggMiEz(i,j,k)
               Ez(i,j,k)=g1(medio)*Ez(i,j,k)+g2(medio)*((Hy(i,j,k)-Hy(i-1,j,k))*Idxh(i)-(Hx(i,j,k)-Hx(i,j-1,k))*Idyhj)
            End do
         End do
      End do
#ifdef SEMBA_FDTD_ENABLE_ACC
   !$acc wait
   !$acc end data
#endif

   end subroutine gpu_advanceEz

   !--------------------------------------------------------------------------------
   ! Advance Hx field - GPU accelerated YEE kernel
   !--------------------------------------------------------------------------------
   subroutine gpu_advanceHx(this, bounds)
      class(gpu_state_t), intent(inout) :: this
      type(bounds_t), intent(in) :: bounds

      real(kind=rkind), dimension(:,:,:), pointer, contiguous  :: Hx
      real(kind=rkind), dimension(:,:,:), pointer, contiguous  :: Ey
      real(kind=rkind), dimension(:,:,:), pointer, contiguous  :: Ez
      integer(kind=integersizeofmediamatrices), dimension(:,:,:), pointer, contiguous :: sggMiHx
      real(kind=rkind), dimension(:), pointer :: Idye
      real(kind=rkind), dimension(:), pointer :: Idze
      real(kind=rkind), dimension(:), pointer :: gm1
      real(kind=rkind), dimension(:), pointer :: gm2

      real(kind=rkind) :: Idzek, Idyej
      integer(kind=4) :: i, j, k
      integer(kind=integersizeofmediamatrices) :: medio

      Hx(0:bounds%Hx%NX-1,0:bounds%Hx%NY-1,0:bounds%Hx%NZ-1) => this%Hx
      Ey(0:bounds%Ey%NX-1,0:bounds%Ey%NY-1,0:bounds%Ey%NZ-1) => this%Ey
      Ez(0:bounds%Ez%NX-1,0:bounds%Ez%NY-1,0:bounds%Ez%NZ-1) => this%Ez
      sggMiHx(0:bounds%sggMiHx%NX-1,0:bounds%sggMiHx%NY-1,0:bounds%sggMiHx%NZ-1) => this%sggMiHx
      Idye(0:bounds%dyE%NY-1) => this%Idye
      Idze(0:bounds%dzE%NZ-1) => this%Idze
      gm1 => this%gm1
      gm2 => this%gm2

#ifdef SEMBA_FDTD_ENABLE_ACC
   !$ACC data copy(Hx) copyin(Ey,Ez,sggMiHx,Idye,Idze,gm1,gm2)
   !$ACC parallel loop gang vector collapse(3) independent private(i,j,k,medio,Idzek,Idyej)
#endif
      Do k=1,bounds%sweepHx%NZ
         Do j=1,bounds%sweepHx%NY
            Do i=1,bounds%sweepHx%NX
               Idzek=Idze(k)
               Idyej=Idye(j)
               medio =sggMiHx(i,j,k)
               Hx(i,j,k)=gm1(medio)*Hx(i,j,k)+gm2(medio)*((Ey(i,j,k+1)-Ey(i,j,k))*Idzek-(Ez(i,j+1,k)-Ez(i,j,k))*Idyej)
            End do
         End do
      End do
#ifdef SEMBA_FDTD_ENABLE_ACC
   !$acc wait
   !$acc end data
#endif

   end subroutine gpu_advanceHx

   !--------------------------------------------------------------------------------
   ! Advance Hy field - GPU accelerated YEE kernel
   !--------------------------------------------------------------------------------
   subroutine gpu_advanceHy(this, bounds)
      class(gpu_state_t), intent(inout) :: this
      type(bounds_t), intent(in) :: bounds

      real(kind=rkind), dimension(:,:,:), pointer, contiguous :: Hy
      real(kind=rkind), dimension(:,:,:), pointer, contiguous :: Ez
      real(kind=rkind), dimension(:,:,:), pointer, contiguous :: Ex
      integer(kind=integersizeofmediamatrices), dimension(:,:,:), pointer, contiguous :: sggMiHy
      real(kind=rkind), dimension(:), pointer :: Idze
      real(kind=rkind), dimension(:), pointer :: Idxe
      real(kind=rkind), dimension(:), pointer :: gm1
      real(kind=rkind), dimension(:), pointer :: gm2

      real(kind=rkind) :: Idzek
      integer(kind=4) :: i, j, k
      integer(kind=integersizeofmediamatrices) :: medio

      Hy(0:bounds%Hy%NX-1,0:bounds%Hy%NY-1,0:bounds%Hy%NZ-1) => this%Hy
      Ez(0:bounds%Ez%NX-1,0:bounds%Ez%NY-1,0:bounds%Ez%NZ-1) => this%Ez
      Ex(0:bounds%Ex%NX-1,0:bounds%Ex%NY-1,0:bounds%Ex%NZ-1) => this%Ex
      sggMiHy(0:bounds%sggMiHy%NX-1,0:bounds%sggMiHy%NY-1,0:bounds%sggMiHy%NZ-1) => this%sggMiHy
      Idze(0:bounds%dzE%NZ-1) => this%Idze
      Idxe(0:bounds%dxE%NX-1) => this%Idxe
      gm1 => this%gm1
      gm2 => this%gm2

#ifdef SEMBA_FDTD_ENABLE_ACC
   !$ACC data copy(Hy) copyin(Ez,Ex,sggMiHy,Idze,Idxe,gm1,gm2)
   !$ACC parallel loop gang vector collapse(3) independent private(i,j,k,medio,Idzek)
#endif
      Do k=1,bounds%sweepHy%NZ
         Do j=1,bounds%sweepHy%NY
            Do i=1,bounds%sweepHy%NX
               Idzek=Idze(k)
               medio =sggMiHy(i,j,k)
               Hy(i,j,k)=gm1(medio)*Hy(i,j,k)+gm2(medio)*((Ez(i+1,j,k)-Ez(i,j,k))*Idxe(i)-(Ex(i,j,k+1)-Ex(i,j,k))*Idzek)
            End do
         End do
      End do
#ifdef SEMBA_FDTD_ENABLE_ACC
   !$acc wait
   !$acc end data
#endif

   end subroutine gpu_advanceHy

   !--------------------------------------------------------------------------------
   ! Advance Hz field - GPU accelerated YEE kernel
   !--------------------------------------------------------------------------------
   subroutine gpu_advanceHz(this, bounds)
      class(gpu_state_t), intent(inout) :: this
      type(bounds_t), intent(in) :: bounds

      real(kind=rkind), dimension(:,:,:), pointer, contiguous :: Hz
      real(kind=rkind), dimension(:,:,:), pointer, contiguous :: Ex
      real(kind=rkind), dimension(:,:,:), pointer, contiguous :: Ey
      integer(kind=integersizeofmediamatrices), dimension(:,:,:), pointer, contiguous :: sggMiHz
      real(kind=rkind), dimension(:), pointer :: Idye
      real(kind=rkind), dimension(:), pointer :: Idxe
      real(kind=rkind), dimension(:), pointer :: gm1
      real(kind=rkind), dimension(:), pointer :: gm2

      real(kind = RKIND) :: Idyej
      integer(kind = 4) :: i, j, k
      integer(kind = INTEGERSIZEOFMEDIAMATRICES) :: medio

      Hz(0:bounds%Hz%NX-1,0:bounds%Hz%NY-1,0:bounds%Hz%NZ-1) => this%Hz
      Ex(0:bounds%EX%NX-1,0:bounds%EX%NY-1,0:bounds%EX%NZ-1) => this%Ex
      Ey(0:bounds%Ey%NX-1,0:bounds%Ey%NY-1,0:bounds%Ey%NZ-1) => this%Ey
      sggMiHz(0:bounds%sggMiHz%NX-1,0:bounds%sggMiHz%NY-1,0:bounds%sggMiHz%NZ-1) => this%sggMiHz
      Idye(0:bounds%dyE%NY-1) => this%Idye
      Idxe(0:bounds%dxE%NX-1) => this%Idxe
      gm1 => this%gm1
      gm2 => this%gm2

#ifdef SEMBA_FDTD_ENABLE_ACC
   !$ACC data copy(Hz) copyin(Ex,Ey,sggMiHz,Idye,Idxe,gm1,gm2)
   !$ACC parallel loop gang vector collapse(3) independent private(i,j,k,medio,Idyej)
#endif
      Do k=1,bounds%sweepHz%NZ
         Do j=1,bounds%sweepHz%NY
            Do i=1,bounds%sweepHz%NX
               Idyej=Idye(j)
               medio =sggMiHz(i,j,k)
               Hz(i,j,k)=gm1(medio)*Hz(i,j,k)+gm2(medio)*((Ex(i,j+1,k)-Ex(i,j,k))*Idyej-(Ey(i+1,j,k)-Ey(i,j,k))*Idxe(i))
            End do
         End do
      End do
#ifdef SEMBA_FDTD_ENABLE_ACC
   !$acc wait
   !$acc end data
#endif

   end subroutine gpu_advanceHz

end module gpu_kernels_m
