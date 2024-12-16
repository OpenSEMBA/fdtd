module idchildtable_mod

#ifdef CompileWithSMBJSON    
   use json_module
   use smbjson_labels_mod, only: J_ID
   use fhash, only: fhash_tbl_t, key=>fhash_key
   use parser_tools_mod, only: json_value_ptr

   type :: IdChildTable_t
      private
      type(fhash_tbl_t) :: idToChilds
   contains
      procedure :: getId 
      procedure :: totalSize
      procedure :: checkId
   end type

   interface IdChildTable_t
      module procedure ctor
   end interface

contains

   function ctor(core, root, path) result(res)
      type(json_core) :: core
      type(json_value), pointer :: root
      character (len=*), intent(in) :: path
      type(IdChildTable_t) :: res
      type(json_value), pointer :: jentries, jentry
      integer :: id
      integer :: i
      logical :: found
      integer :: numberOfEntries

      call core%get(root, path, jentries, found)
      if (.not. found) return
      numberOfEntries = core%count(jentries)
      call res%idToChilds%allocate(10*numberOfEntries)
      do i = 1, numberOfEntries
         call core%get_child(jentries, i, jentry)
         call core%get(jentry, J_ID, id)
         call res%idToChilds%set(key(id), json_value_ptr(jentry))
      end do
   end function

   function totalSize(this) result(res)
      class(IdChildTable_t) :: this
      integer :: res
      call this%idToChilds%stats(num_items=res)
   end function
   
   function checkId(this, id) result(stat)
      class(IdChildTable_t) :: this
      integer, intent(in) :: id
      integer :: stat
      call this%idToChilds%check_key(key(id), stat)
   end function

   function getId(this, id) result(res)
      class(IdChildTable_t) :: this
      type(json_value_ptr) :: res
      integer, intent(in) :: id
      integer :: mStat
      class(*), allocatable :: d
      
      nullify(res%p)
      call this%idToChilds%check_key(key(id), mStat)
      if (mStat /= 0) then
         return
      end if

      call this%idToChilds%get_raw(key(id), d)
      select type(d)
       type is (json_value_ptr)
         res = d
      end select
   end function

#endif
end module
