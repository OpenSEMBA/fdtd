integer function test_idchildtable_fhash() bind(C) result(error_cnt)
   use fhash, only: fhash_tbl_t, key=>fhash_key
   implicit none
   type(fhash_tbl_t) :: tbl
   integer :: val
   integer, parameter :: expectedValue=10
   error_cnt = 0
   
   call tbl%set(key('my_key_1'), value=expectedValue)
   call tbl%set(key('my_key_2'), value=1.0)
   call tbl%set(key(123456), value='a string value')
   call tbl%set(key([1,2,3,4,5]), value=.false.)

   call tbl%get(key('my_key_1'),val)

   if (val /= expectedValue) then 
      error_cnt = error_cnt+1
      return
   end if
end function


integer function test_idchildtable() bind(C) result(err)
   use idchildtable_mod
   use smbjson_labels_mod
   use parser_tools_mod, only: json_value_ptr
   use smbjson_testingTools
   use json_module

   implicit none

   type(IdChildTable_t) :: tbl

   character(len=*),parameter :: filename = PATH_TO_TEST_DATA//'cases/holland1981.fdtd.json'
   type(json_file) :: jsonfile
   type(json_core) :: core
   type(json_value), pointer :: root => null()
   type(json_value_ptr) :: mat

   err = 0
   call jsonfile%initialize()
   if (jsonfile%failed()) then
      err = err + 1
      return
   end if

   call jsonfile%load(filename = filename)
   if (jsonfile%failed()) then
      err = err + 1
      return
   end if

   call jsonfile%get_core(core)
   call jsonfile%get('.', root)

   tbl = IdChildTable_t(core, root, J_MATERIALS)

   if (tbl%count() /= 2) err = err + 1
   block
      character (len=:), allocatable :: matType
      logical :: found

      mat = tbl%getId(1)
      call core%get(mat%p, J_TYPE, matType, found)
      if (.not. found) err = err + 1
      if (matType /= J_MAT_TYPE_WIRE) err = err + 1

      mat = tbl%getId(2)
      call core%get(mat%p, J_TYPE, matType, found)
      if (.not. found) err = err + 1
      if (matType /= J_MAT_TYPE_TERMINAL) err = err + 1
   end block
end function

