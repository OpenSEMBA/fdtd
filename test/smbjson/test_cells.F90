integer function test_cells() bind(C) result(err)
   use cells_mod
   use smbjson_testingTools

   implicit none

   type(cell_interval_t) :: interval
   err = 0

   ! +X oriented linel interval.
   interval = cell_interval_t(ini=cell_t([1, 1, 1]), end=cell_t([3, 1, 1]))
   if (interval%getType()        /= CELL_TYPE_LINEL) err = err + 1
   if (interval%getOrientation() /= + DIR_X) err = err + 1
   if (interval%getSize()        /= 2) err = err + 1

   ! -X oriented linel interval.
   interval = cell_interval_t(ini=cell_t([2, 1, 1]), end=cell_t([0, 1, 1]))
   if (interval%getType()        /= CELL_TYPE_LINEL) err = err + 1
   if (interval%getOrientation() /= - DIR_X) err = err + 1
   if (interval%getSize()        /= 2) err = err + 1

   ! +Y oriented linel interval.
   interval = cell_interval_t(ini=cell_t([11, 16, 11]), end=cell_t([11, 21, 11]))
   if (interval%getType()        /= CELL_TYPE_LINEL) err = err + 1
   if (interval%getOrientation() /= + DIR_Y) err = err + 1
   if (interval%getSize()        /= 5) err = err + 1


   ! 2x2 surfel, with +Z normal.
   interval = cell_interval_t(ini=cell_t([3, 3, 3]), end=cell_t([5, 5, 3]))
   if (interval%getType()        /= CELL_TYPE_SURFEL) err = err + 1
   if (interval%getOrientation() /= + DIR_Z) err = err + 1
   if (interval%getSize()        /= 4) err = err + 1

   ! 3x1 surfel, with -X normal.
   interval = cell_interval_t(ini=cell_t([5, 5, 5]), end=cell_t([5, 4, 2]))
   if (interval%getType()        /= CELL_TYPE_SURFEL) err = err + 1
   if (interval%getOrientation() /= - DIR_X) err = err + 1
   if (interval%getSize()        /= 3) err = err + 1
end function