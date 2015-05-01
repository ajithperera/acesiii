C  Copyright (c) 2003-2010 University of Florida
C
C  This program is free software; you can redistribute it and/or modify
C  it under the terms of the GNU General Public License as published by
C  the Free Software Foundation; either version 2 of the License, or
C  (at your option) any later version.

C  This program is distributed in the hope that it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C  GNU General Public License for more details.

C  The GNU General Public License is included in this distribution
C  in the file COPYRIGHT.
      subroutine Place_constanta(array_table, narray_table, 
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table,
     *                      address_table, op)

      implicit none
      include 'interpreter.h'
      include 'mpif.h'
      include 'proto_events.h'
      include 'trace.h'
      include 'parallel_info.h'
#ifdef ALTIX
      include 'sheap.h'
#endif
      integer narray_table, nindex_table, nsegment_table,
     *        nblock_map_table
      integer op(loptable_entry)
      integer array_table(larray_table_entry, narray_table)
      integer index_table(lindex_table_entry, nindex_table)
      integer segment_table(lsegment_table_entry, nsegment_table)
      integer block_map_table(lblock_map_entry, nblock_map_table)
      integer*8 address_table(narray_table),get_index_from_base
      integer nscalar_table
      double precision scalar_table(nscalar_table)

      integer eps, i, n
      integer ierr
      integer source1,source2 
      integer array_type, ind_1(mx_array_index), n_1(mx_array_index)
      integer sval1(mx_array_index),sval2(mx_array_index)
      integer sind 
      integer atype,btype
      integer*8 indblk, get_block_index
      integer*8 isource1,isource2
      integer nindex_source1,nindex_source2
      integer block, blkndx, seg,stack
      integer master, pst_get_master
      integer status(mpi_status_size)
      integer find_current_block
      integer ind(mx_array_index)

      integer type(mx_array_index)
      logical msg_present

      include 'int_gen_parms.h'
      include 'epsilon.h'

      integer company_comm, pst_get_company_comm

      double precision xxx(1), input_constant 
      integer ibuf(1)
#ifdef ALTIX
      pointer (iptr, ibuf)
      pointer (dptr, xxx)
#else
      common xxx
      equivalence (ibuf, xxx)
#endif

c-------------------------------------------------------------------------
c   Locate array address in array_table.
c-------------------------------------------------------------------------

#ifdef ALTIX
      iptr = ishptr
      dptr = dshptr
#endif

c---------------------------------------------------------------------------
c   Find the indices of the array block.
c---------------------------------------------------------------------------

      source1 = op(c_result_array)
      source2 = op(c_op1_array)

      atype = array_table(c_array_type, source1)
      btype = array_table(c_array_type, source2)

      If (atype .NE. static_array) then
          print *,'Error: First argument must be static '
          print *, " in Place_constanta.F"
         call abort_job()
      endif

      nindex_source1 = array_table(c_nindex, source1)
      nindex_source2 = array_table(c_nindex, source2)

      if (nindex_source1 .ne. 3) then
         print *,"Error: In the routine Place_constant.F"
         print *,'Error: First array dimensions must be 3'
         print *,'First array has ',nindex_source1,' indices.'
         call abort_job()
      endif

c-------------------------------------------------------------------------
c   Get segment ranges of the array source.
c-------------------------------------------------------------------------

      do i = 1, nindex_source1
         ind_1(i) = array_table(c_index_array1+i-1,source1)
         n_1(i)   = index_table(c_index_size, ind_1(i))  
         type(i)  = index_table(c_index_type, ind_1(i))
         seg      = index_table(c_current_seg,ind_1(i))
         call get_index_segment(ind_1(i), seg, segment_table,
     *                          nsegment_table, index_table,
     *                          nindex_table, sval1(i), sval2(i))
      enddo

      isource1 = get_index_from_base(address_table(source1), xxx, 2)

c-------------------------------------------------------------------------
c   Get segment ranges of the array target.
c-------------------------------------------------------------------------

      source2 = op(c_op1_array)
      btype   = array_table(c_array_type, source2)

      if (btype .ne. scalar_value) then
         print *,'Error: scalar not sent in Place_constant'
         print *,(op(i),i=1,loptable_entry)
         call abort_job()
      endif

      sind =  array_table(c_scalar_index, source2)

      if (sind .lt. 1 .or. sind .gt. nscalar_table) then
         print *,'Scalar table index out of range in ',
     *           ' Place_constant.F line',current_line
         print *,'Index for array ',source2,' is ',sind,' should be ',
     *           'between 1 and ',nscalar_table
         call abort_job()
      endif

      input_constant = scalar_table(sind)

      call assign_to_array_elementa(xxx(isource1),input_constant, 
     *                             sval1(1), sval2(1), 
     *                             sval1(2), sval2(2),
     *                             sval1(3), sval2(3))

      return
      end
C-------------------------------------------------------------------------

      subroutine assign_to_array_elementa(xxx,input_constant,
     *                                   pa1,pa2,pb1,pb2,
     *                                   pc1,pc2)
      implicit none
      include 'interpreter.h'
      include 'int_gen_parms.h'
      include 'mpif.h'
      include 'proto_events.h'
      include 'parallel_info.h'
C
      integer pa1,pa2,pb1,pb2,pc1,pc2,a,b,c
      double precision xxx(pa1:pa2,pb1:pb2,pc1:pc2)
      double precision input_constant
C
      if (pb1 .eq. 1 .and. pb2 .eq. 1) then

      do c = pc1, pc2
      do a = pa1, pa2

        if(a .eq. nalpha_occupied+c) then

           xxx(a,1,c) = input_constant

C#ifdef _DEBUG_LVL0
           write(6,*) " Single excitation(A): " , a,pb1,c, xxx(a,1,c)
C#endif 
        endif

      enddo
      enddo

      else
     
      do c = pc1, pc2
      do b = pb1, pb2
      do a = pa1, pa2

        if(a .eq. nalpha_occupied+c .and. b .eq. nalpha_occupied-c) then

           xxx(a,b,c) = input_constant

C#ifdef _DEBUG_LVL0
           write(6,*) " Single excitation(A): ", a, b, c, xxx(a,b,c)
C#endif 
        endif

      enddo
      enddo
      enddo 

      endif 

      return
      end
