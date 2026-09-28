integer function test_conformal_sgbc_zero_state() bind(C) result(err)
   use conformal_sgbc_m
   use FDETYPES_m
   implicit none
   type(SGBCMaterialProfile_t) :: profile
   type(conformal_sgbc_state_t) :: state

   err = 0
   call make_profile(profile, .false.)
   call state%init(profile, 1.0e-12_RKIND_tiempo, 0.005_RKIND, 0.005_RKIND, &
      1.0e9_RKIND, 1.0_RKIND, 2, .false.)
   call state%advance(0.0_RKIND, 0.0_RKIND)
   if (.not. state%initialized) err = err+1
   if (maxval(abs(state%e)) /= 0.0_RKIND) err = err+1
   if (maxval(abs(state%h)) /= 0.0_RKIND) err = err+1
contains
   subroutine make_profile(profile, asymmetric)
      type(SGBCMaterialProfile_t), intent(out) :: profile
      logical, intent(in) :: asymmetric
      profile%num_layers = 1
      allocate(profile%thickness(1), profile%eps(1), profile%mu(1), profile%sigma(1), profile%sigmam(1))
      profile%thickness = 0.001_RKIND
      profile%eps = EPSILON_VACUUM
      profile%mu = MU_VACUUM
      profile%sigma = 10.0_RKIND
      profile%sigmam = 0.0_RKIND
   end subroutine
end function test_conformal_sgbc_zero_state

integer function test_conformal_sgbc_layer_orientation() bind(C) result(err)
   use conformal_sgbc_m
   use FDETYPES_m
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   type(SGBCMaterialProfile_t) :: profile
   type(conformal_sgbc_state_t) :: forward, reverse
   integer :: iteration

   err = 0
   profile%num_layers = 2
   allocate(profile%thickness(2), profile%eps(2), profile%mu(2), profile%sigma(2), profile%sigmam(2))
   profile%thickness = [0.0005_RKIND, 0.0015_RKIND]
   profile%eps = [2.0_RKIND*EPSILON_VACUUM, EPSILON_VACUUM]
   profile%mu = MU_VACUUM
   profile%sigma = [0.0_RKIND, 100.0_RKIND]
   profile%sigmam = 0.0_RKIND
   call forward%init(profile, 1.0e-12_RKIND_tiempo, 0.005_RKIND, 0.005_RKIND, &
      1.0e9_RKIND, 1.0_RKIND, 2, .false.)
   call reverse%init(profile, 1.0e-12_RKIND_tiempo, 0.005_RKIND, 0.005_RKIND, &
      1.0e9_RKIND, 1.0_RKIND, 2, .true.)
   do iteration = 1, 8
      call forward%advance(1.0_RKIND, 0.0_RKIND)
      call reverse%advance(1.0_RKIND, 0.0_RKIND)
   end do
   if (.not. all(ieee_is_finite(forward%e)) .or. .not. all(ieee_is_finite(reverse%e))) err = err+1
   if (maxval(abs(forward%e-reverse%e)) <= 100.0_RKIND*epsilon(1.0_RKIND)) err = err+1
end function test_conformal_sgbc_layer_orientation

integer function test_conformal_sgbc_geometry_winding() bind(C) result(err)
   use conformal_m, only: buildSGBCSurfaceMedia
   use NFDETypes_m, only: ConformalPECElements_t, ConformalMedia_t
   use conformal_types_m, only: coord_t, triangle_t
   implicit none
   type(ConformalPECElements_t), dimension(:), pointer :: elements
   type(ConformalMedia_t), dimension(:), allocatable :: media
   type(coord_t) :: c1, c2, c3

   err = 0
   c1 = coord_t(position=[0.75, 0.0, 0.0], id=1)
   c2 = coord_t(position=[0.75, 0.0, 1.0], id=2)
   c3 = coord_t(position=[0.75, 1.0, 0.0], id=3)
   allocate(elements(1))
   allocate(elements(1)%triangles(1), elements(1)%intervals(0))
   elements(1)%triangles(1) = triangle_t(vertices=[c1, c2, c3])

   media = buildSGBCSurfaceMedia(elements)
   call check_faces(media(1), -1, err)

   elements(1)%triangles(1) = triangle_t(vertices=[c1, c3, c2])
   media = buildSGBCSurfaceMedia(elements)
   call check_faces(media(1), 1, err)
contains
   subroutine check_faces(item, expected_sign, err)
      type(ConformalMedia_t), intent(in) :: item
      integer, intent(in) :: expected_sign
      integer, intent(inout) :: err
      integer :: group, face
      if (.not. item%is_sgbc) err = err+1
      do group = 1, item%n_faces_media
         do face = 1, item%face_media(group)%size
            if (.not. item%face_media(group)%faces(face)%is_two_sided) err = err+1
            if (item%face_media(group)%faces(face)%surface_normal_sign /= expected_sign) err = err+1
         end do
      end do
   end subroutine check_faces
end function test_conformal_sgbc_geometry_winding

integer function test_conformal_sgbc_rejects_unsplit_geometry() bind(C) result(err)
   use conformal_m, only: buildSGBCSurfaceMedia
   use NFDETypes_m, only: ConformalPECElements_t, ConformalMedia_t
   use conformal_types_m, only: coord_t, triangle_t
   use Report_m, only: isFatalError, resetFatalError
   implicit none
   type(ConformalPECElements_t), dimension(:), pointer :: elements
   type(ConformalMedia_t), dimension(:), allocatable :: media
   type(coord_t) :: c1, c2, c3

   err = 0
   call resetFatalError()
   c1 = coord_t(position=[0.0, 0.0, 0.0], id=1)
   c2 = coord_t(position=[0.0, 1.0, 0.0], id=2)
   c3 = coord_t(position=[0.0, 0.0, 1.0], id=3)
   allocate(elements(1))
   allocate(elements(1)%triangles(1), elements(1)%intervals(0))
   elements(1)%triangles(1) = triangle_t(vertices=[c1, c2, c3])
   media = buildSGBCSurfaceMedia(elements)
   if (.not. isFatalError()) err = err+1
   call resetFatalError()
end function test_conformal_sgbc_rejects_unsplit_geometry
