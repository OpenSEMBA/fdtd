integer function test_parser_ctor() bind(C) result(err)
   use smbjson_m
   use smbjson_testingTools_m

   implicit none

   character(len=*),parameter :: FILENAME = PATH_TO_TEST_DATA//INPUT_EXAMPLES//'planewave.fdtd.json'
   type(parser_t) :: parser

   parser = parser_t(filename)
   if (parser%isInitialized) then
      err = 0
   else
      err = 1
   end if
end function

integer function test_parser_tools_interval_to_coords() result(err)
   use parser_tools_m
   use smbjson_testingTools_m

   implicit none

   type(cell_interval_t) :: interval
   type(coords_t), dimension(:), allocatable :: cs
   err = 0

   ! +Y oriented linel interval.
   interval = cell_interval_t( &
      startNode=cell_t([10, 0, 10]), endNode=cell_t([10, 5, 10]))
   cs = cellIntervalsToCoords([interval])

   call expect_eq_int(err,  1, size(cs))
   call expect_eq_int(err, 10, cs(1)%Xi)
   call expect_eq_int(err, 10, cs(1)%Xe)
   call expect_eq_int(err,  0, cs(1)%Yi)
   call expect_eq_int(err,  4, cs(1)%Ye)
   call expect_eq_int(err, 10, cs(1)%Zi)
   call expect_eq_int(err, 10, cs(1)%Ze)
   call expect_eq_int(err, +IEY, cs(1)%Or)

   ! -Z oriented linel interval.
   interval = cell_interval_t( &
      startNode=cell_t([10, 10, 5]), endNode=cell_t([10, 10, -1]))
   cs = cellIntervalsToCoords([interval])

   call expect_eq_int(err,  1, size(cs))
   call expect_eq_int(err, 10, cs(1)%Xi)
   call expect_eq_int(err, 10, cs(1)%Xe)
   call expect_eq_int(err, 10, cs(1)%Yi)
   call expect_eq_int(err, 10, cs(1)%Ye)
   call expect_eq_int(err,  0, cs(1)%Zi)
   call expect_eq_int(err,  4, cs(1)%Ze)
   call expect_eq_int(err, -IEZ, cs(1)%Or)

   ! +Y oriented surfel interval.
   interval = cell_interval_t( &
      startNode=cell_t([ 9,  2,  9]), endNode=cell_t([11,  2, 11]) )
   cs = cellIntervalsToCoords([interval])
   call expect_eq_int(err,  1, size(cs))
   call expect_eq_int(err,  9, cs(1)%Xi)
   call expect_eq_int(err, 10, cs(1)%Xe)
   call expect_eq_int(err,  2, cs(1)%Yi)
   call expect_eq_int(err,  2, cs(1)%Ye)
   call expect_eq_int(err,  9, cs(1)%Zi)
   call expect_eq_int(err, 10, cs(1)%Ze)
   call expect_eq_int(err, +IEY, cs(1)%Or)


end function

integer function test_parser_read_mesh() bind(C) result(err)

   use smbjson_m
   use smbjson_testingTools_m

   implicit none

   character(len=*),parameter :: FILENAME = PATH_TO_TEST_DATA//INPUT_EXAMPLES//'mtln.fdtd.json'
   type(parser_t) :: parser
   type(mesh_t) :: mesh
   logical :: found
   type(coordinate_t) :: expected, obtained

   err = 0

   parser = parser_t(filename)
   mesh = parser%readMesh()

   ! call mesh%printCoordHashInfo() !! For debugging only
   expected%position = [10,0,1]

   obtained = mesh%getCoordinate(59, found)
   if (.not. found) err = err + 1
   if (any(obtained%position /= expected%position)) err = err + 1

   obtained = mesh%getCoordinate(64, found)
   if (.not. found) err = err + 1
   if (any(obtained%position /= expected%position)) err = err + 1

   obtained = mesh%getCoordinate(61, found)
   if (.not. found) err = err + 1
   if (any(obtained%position /= expected%position)) err = err + 1


end function

integer function test_parser_read_conformal_volume() bind(C) result(err)

   use smbjson_m
   use smbjson_testingTools_m
   
   implicit none

   character(len=*),parameter :: FILENAME = PATH_TO_TEST_DATA//INPUT_EXAMPLES//'conformal.fdtd.json'
   type(parser_t) :: parser
   type(mesh_t) :: mesh
   type(Parseador_t) :: problem
   logical :: found
   type(conformal_region_t), dimension(:), allocatable :: conformal_regions
   type(cell_region_t), dimension(:), allocatable :: cell_regions
   real, dimension(3) :: expected

   err = 0

   parser = parser_t(filename)
   mesh = parser%readMesh()
   conformal_regions = mesh%getConformalRegions([5])
   if (size(conformal_regions) /= 1) err = err + 1
   if (size(conformal_regions(1)%triangles) /= 24) err = err + 1
   expected = [0.0, 0.0, 4.25]
   if (any(conformal_regions(1)%triangles(1)%vertices(1)%position /= expected)) err = err + 1
   cell_regions = mesh%getCellRegions([5])
   if (size(cell_regions(1)%intervals) /= 1) err = err + 1

   problem = parser%readProblemDescription()
   if (.not. associated(problem%conformalRegs%volumes)) err = err + 1
   if (size(problem%conformalRegs%volumes) /= 1) err = err + 1
end function

integer function test_parser_read_conformal_local_coordinates() bind(C) result(err)

   use smbjson_m
   use smbjson_testingTools_m

   implicit none

   character(len=*),parameter :: FILENAME = PATH_TO_TEST_DATA//INPUT_EXAMPLES//'conformal_local_coordinates.fdtd.json'
   type(parser_t) :: parser
   type(mesh_t) :: mesh
   type(conformal_region_t), dimension(:), allocatable :: conformal_regions
   real, dimension(3) :: expected

   err = 0

   parser = parser_t(filename)
   mesh = parser%readMesh()
   conformal_regions = mesh%getConformalRegions([10, 11])
   if (size(conformal_regions) /= 2) err = err + 1
   if (size(conformal_regions(1)%triangles) /= 1) err = err + 1
   if (size(conformal_regions(2)%triangles) /= 1) err = err + 1

   expected = [1.0, 1.0, 1.0]
   if (any(conformal_regions(1)%triangles(1)%vertices(1)%position /= expected)) err = err + 1
   expected = [3.0, 3.0, 2.0]
   if (any(conformal_regions(2)%triangles(1)%vertices(1)%position /= expected)) err = err + 1
end function

integer function test_parser_reject_conformal_missing_coordinates() bind(C) result(err)

   use smbjson_m
   use smbjson_testingTools_m
   use Report_m, only: isFatalError, resetFatalError

   implicit none

   character(len=*),parameter :: FILENAME = PATH_TO_TEST_DATA//INPUT_EXAMPLES//'conformal_missing_coordinates.fdtd.json'
   type(parser_t) :: parser
   type(mesh_t) :: mesh

   err = 0
   call resetFatalError()
   parser = parser_t(filename)
   mesh = parser%readMesh()
   if (.not. isFatalError()) err = err + 1
   call resetFatalError()
end function

integer function test_parser_reject_conformal_unknown_coordinate() bind(C) result(err)

   use smbjson_m
   use smbjson_testingTools_m
   use Report_m, only: isFatalError, resetFatalError

   implicit none

   character(len=*),parameter :: FILENAME = PATH_TO_TEST_DATA//INPUT_EXAMPLES//'conformal_unknown_coordinate.fdtd.json'
   type(parser_t) :: parser
   type(mesh_t) :: mesh

   err = 0
   call resetFatalError()
   parser = parser_t(filename)
   mesh = parser%readMesh()
   if (.not. isFatalError()) err = err + 1
   call resetFatalError()
end function

integer function test_parser_reject_conformal_duplicate_coordinate() bind(C) result(err)

   use smbjson_m
   use smbjson_testingTools_m
   use Report_m, only: isFatalError, resetFatalError

   implicit none

   character(len=*),parameter :: FILENAME = PATH_TO_TEST_DATA//INPUT_EXAMPLES//'conformal_duplicate_coordinates.fdtd.json'
   type(parser_t) :: parser
   type(mesh_t) :: mesh

   err = 0
   call resetFatalError()
   parser = parser_t(filename)
   mesh = parser%readMesh()
   if (.not. isFatalError()) err = err + 1
   call resetFatalError()
end function

integer function test_parser_reject_conformal_nonpec_material() bind(C) result(err)
   use smbjson_m
   use smbjson_testingTools_m
   use Report_m, only: isFatalError, resetFatalError
   implicit none

   type(parser_t) :: parser
   type(Parseador_t) :: problem

   err = 0
   call resetFatalError()
   parser = parser_t(PATH_TO_TEST_DATA//INPUT_EXAMPLES//'conformal_nonpec_material.fdtd.json')
   problem = parser%readProblemDescription()
   if (.not. isFatalError()) err = err + 1
   call resetFatalError()
end function
