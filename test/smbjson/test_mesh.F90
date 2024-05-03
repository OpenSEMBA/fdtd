integer function test_mesh_add_get() bind(C) result(error_cnt)
   use mesh_mod
   use smbjson_testingTools

   implicit none

   type(mesh_t) :: mesh
   logical :: found

   error_cnt = 0

   call mesh%addCoordinate(10, coordinate_t([0.0, 0.0, 0.0]))
   call mesh%addCoordinate(11, coordinate_t([1.0, 0.0, 0.0]))
   call mesh%addCoordinate(12, coordinate_t([2.0, 0.0, 0.0]))
   call mesh%addCoordinate(13, coordinate_t([3.0, 0.0, 0.0]))

   block
      type(node_t) :: expected, obtained      
      expected%coordIds = [10]
      call mesh%addElement(1, expected)
      obtained = mesh%getNode(1, found)
      if (.not. found) error_cnt = error_cnt + 1
      if (any(obtained%coordIds /= expected%coordIds)) error_cnt = error_cnt + 1
   end block

   block 
      type(node_t) :: obtained
      obtained = mesh%getNode(102, found)
      if (found) error_cnt = error_cnt + 1
   end block

   block
      type(polyline_t) :: expected, obtained
      expected%coordIds = [11, 12, 13]
      call mesh%addElement(2, expected)
      obtained = mesh%getPolyline(2, found)
      if (.not. found) error_cnt = error_cnt + 1
      if (any(obtained%coordIds /= expected%coordIds)) error_cnt = error_cnt + 1
   end block

end function

integer function test_mesh_add_get_long_list() bind(C) result(error_cnt)
   use mesh_mod
   use smbjson_testingTools

   implicit none

   type(mesh_t) :: mesh
   logical :: found

   error_cnt = 0
   call mesh%addCoordinate( 1 , coordinate_t([  1, 9, 1]))
   call mesh%addCoordinate( 2 , coordinate_t([ 10, 9, 1]))
   call mesh%addCoordinate( 5 , coordinate_t([ 18, 9, 1])) 
   call mesh%addCoordinate( 6 , coordinate_t([ 10, 2, 1]))    
   call mesh%addCoordinate( 11, coordinate_t([  1, 9, 1]))
   call mesh%addCoordinate( 15, coordinate_t([ 10, 9, 1]))
   call mesh%addCoordinate( 23, coordinate_t([ 18, 9, 1])) 
   call mesh%addCoordinate( 24, coordinate_t([ 10, 2, 1])) 
   call mesh%addCoordinate( 33, coordinate_t([  1, 9, 1]))
   call mesh%addCoordinate( 34, coordinate_t([  1, 9, 1]))
   call mesh%addCoordinate( 35, coordinate_t([  1, 9, 1]))
   call mesh%addCoordinate( 36, coordinate_t([  1, 9, 1]))
   call mesh%addCoordinate( 37, coordinate_t([  1, 9, 1]))
   call mesh%addCoordinate( 38, coordinate_t([  1, 9, 1]))
   call mesh%addCoordinate( 39, coordinate_t([  1, 9, 1]))
   call mesh%addCoordinate( 40, coordinate_t([  1, 9, 1]))
   call mesh%addCoordinate( 41, coordinate_t([ 10, 9, 1]))
   call mesh%addCoordinate( 42, coordinate_t([ 10, 9, 1]))
   call mesh%addCoordinate( 43, coordinate_t([ 10, 9, 1]))
   call mesh%addCoordinate( 44, coordinate_t([ 10, 9, 1]))
   call mesh%addCoordinate( 45, coordinate_t([ 10, 9, 1]))
   call mesh%addCoordinate( 46, coordinate_t([ 10, 9, 1]))
   call mesh%addCoordinate( 47, coordinate_t([ 10, 9, 1]))
   call mesh%addCoordinate( 48, coordinate_t([ 10, 9, 1]))
   call mesh%addCoordinate( 51, coordinate_t([ 18, 9, 1]))
   call mesh%addCoordinate( 52, coordinate_t([ 18, 9, 1]))
   call mesh%addCoordinate( 59, coordinate_t([ 10, 2, 1])) 
   call mesh%addCoordinate( 60, coordinate_t([ 10, 2, 1])) 
   call mesh%addCoordinate( 61, coordinate_t([ 10, 2, 1])) 
   call mesh%addCoordinate( 62, coordinate_t([ 10, 2, 1])) 
   call mesh%addCoordinate( 63, coordinate_t([ 10, 2, 1]))
   call mesh%addCoordinate( 64, coordinate_t([ 10, 2, 1]))
   call mesh%addCoordinate( 65, coordinate_t([ 10, 2, 1]))
   call mesh%addCoordinate( 66, coordinate_t([ 10, 2, 1]))


   block
      type(coordinate_t) :: expected, obtained      
      expected%position = [10,2,1]
      obtained = mesh%getCoordinate(65, found)
      if (.not. found) error_cnt = error_cnt + 1
      if (any(obtained%position /= expected%position)) error_cnt = error_cnt + 1
   end block

   block
      type(coordinate_t) :: expected, obtained      
      expected%position = [10,2,1]
      obtained = mesh%getCoordinate(66, found)
      if (.not. found) error_cnt = error_cnt + 1
      if (any(obtained%position /= expected%position)) error_cnt = error_cnt + 1
   end block

   call mesh%printCoordHashInfo()

end function

integer function test_mesh_node_to_pixel() bind(C) result(err)
   use mesh_mod
   use cells_mod
   use smbjson_testingTools

   implicit none

   type(mesh_t) :: mesh

   err = 0

   call mesh%addCoordinate(1, coordinate_t([0.0, 0.0, 0.0]))

   block
      ! Valid pixel conversion.
      type(pixel_t) :: pix
      pix = mesh%convertNodeToPixel(node_t([1]))
      if (.not. pixel_t(cell=[0, 0, 0], tag=1) == pix) err = err + 1
   end block
end function

integer function test_mesh_polyline_to_linel() bind(C) result(err)
   use mesh_mod
   use cells_mod
   use smbjson_testingTools

   implicit none

   type(mesh_t) :: mesh

   err = 0

   call mesh%addCoordinate(1, coordinate_t([0.0, 0.0, 0.0]))
   call mesh%addCoordinate(2, coordinate_t([3.0, 0.0, 0.0]))
   call mesh%addCoordinate(3, coordinate_t([3.0, 1.0, 0.0]))
   call mesh%addCoordinate(4, coordinate_t([3.0, 0.0, 0.0]))
   call mesh%addCoordinate(5, coordinate_t([3.5, 0.0, 0.0]))

   block
      !   Polyline 1 (valid linel conversion)
      !             3
      !             ^
      !             | 4
      !    1  - - > 2
      !       1 2 3
      type(linel_t), dimension(:), allocatable :: ls
      type(polyline_t) :: pl

      pl = polyline_t([1, 2, 3])
      ls = mesh%convertPolylineToLinels(pl)
      if (.not. allocated(ls)) err = err + 1
      if (size(ls) /= 4) err = err + 1

      if (.not. linel_t(cell=[0, 0, 0], orientation=DIR_X, tag=1) == ls(1)) err = err + 1
      if (.not. linel_t(cell=[1, 0, 0], orientation=DIR_X         ) == ls(2)) err = err + 1
      if (.not. linel_t(cell=[2, 0, 0], orientation=DIR_X         ) == ls(3)) err = err + 1
      if (.not. linel_t(cell=[3, 0, 0], orientation=DIR_Y, tag=3) == ls(4)) err = err + 1

      if (mesh%arePolylineSegmentsStructured(pl) .neqv. .true.) err = err + 1
   end block

   block
      !   Polyline 2 (valid linel conversion, reverse orientation)
      !             3
      !             | 1
      !    1  < - - 2
      !       4 3 2
      type(linel_t), dimension(:), allocatable :: ls
      type(polyline_t) :: pl

      pl = polyline_t([3, 2, 1])
      ls = mesh%convertPolylineToLinels(pl)

      if (.not. allocated(ls)) err = err + 1
      if (size(ls) /= 4) err = err + 1

      if (.not. linel_t(cell=[3, 0, 0], orientation=-DIR_Y, tag=3) == ls(1)) err = err + 1
      if (.not. linel_t(cell=[2, 0, 0], orientation=-DIR_X, tag=2) == ls(2)) err = err + 1
      if (.not. linel_t(cell=[1, 0, 0], orientation=-DIR_X       ) == ls(3)) err = err + 1
      if (.not. linel_t(cell=[0, 0, 0], orientation=-DIR_X, tag=1) == ls(4)) err = err + 1

      if (mesh%arePolylineSegmentsStructured(pl) .neqv. .true.) err = err + 1
   end block

   block
      !   Polyline 3 (invalid linel conversion, does not follow grid)
      !        3
      !      /
      !    1
      type(linel_t), dimension(:), allocatable :: ls
      type(polyline_t) :: pl
      pl = polyline_t([1, 3])
      ls = mesh%convertPolylineToLinels(pl)
      if (size(ls) /= 0) err = err + 1

      if (mesh%arePolylineSegmentsStructured(pl) .neqv. .false.) err = err + 1
   end block

   block
      !   Polyline 4 (polyline with coord in same place.)
      !    1  - - > 2,4
      type(linel_t), dimension(:), allocatable :: ls
      type(polyline_t) :: pl
      pl = polyline_t([1, 2, 4])
      ls = mesh%convertPolylineToLinels(pl)

      if (.not. allocated(ls)) err = err + 1
      if (size(ls) /= 3) err = err + 1
      if (.not. linel_t(cell=[0, 0, 0], orientation=DIR_X, tag=1) == ls(1)) err = err + 1
      if (.not. linel_t(cell=[1, 0, 0], orientation=DIR_X         ) == ls(2)) err = err + 1
      if (.not. linel_t(cell=[2, 0, 0], orientation=DIR_X, tag=4) == ls(3)) err = err + 1

      if (mesh%arePolylineSegmentsStructured(pl) .neqv. .true.) err = err + 1
   end block

   block
      !   Polyline 5 (polyline with coord in same place.)
      !    1 < - -  2,4
      type(linel_t), dimension(:), allocatable :: ls
      type(polyline_t) :: pl
      pl = polyline_t([4, 2, 1])
      ls = mesh%convertPolylineToLinels(pl)

      if (.not. allocated(ls)) err = err + 1
      if (size(ls) /= 3) err = err + 1
      if (.not. linel_t(cell=[2, 0, 0], orientation=-DIR_X, tag=4) == ls(1)) err = err + 1
      if (.not. linel_t(cell=[1, 0, 0], orientation=-DIR_X       ) == ls(2)) err = err + 1
      if (.not. linel_t(cell=[0, 0, 0], orientation=-DIR_X, tag=1) == ls(3)) err = err + 1

      if (mesh%arePolylineSegmentsStructured(pl) .neqv. .true.) err = err + 1
   end block

   block
      !   Polyline 6 (polyline with fractional coord.)
      !    1 - - -> 5
      if (mesh%arePolylineSegmentsStructured(polyline_t([1, 5])) .eqv. .true.) err = err + 1
   end block

end function

