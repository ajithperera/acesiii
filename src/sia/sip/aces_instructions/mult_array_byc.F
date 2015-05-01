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
      subroutine Mult_Array_BYC(array_table, narray_table, 
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table,
     *                      address_table, op)

      implicit none
      include 'interpreter.h'
      include 'mpif.h'
      include 'proto_events.h'
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
      integer array_type, inda(2), indb(2)
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

      double precision xxx(1)
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

      If (atype .NE. static_array .OR. 
     *    btype .NE. static_array) then
          print *,'Error: Both arrays must be static in'
          print *, "Mult_Array_BYC"
         call abort_job()
      endif

      nindex_source1 = array_table(c_nindex, source1)
      nindex_source2 = array_table(c_nindex, source2)

      if (nindex_source1 .ne. 2 .OR. nindex_source2 .ne. 1) then
         print *,"Error: In the routine Mult_Array_BYC SIP"
         print *,'Error: First array dimensions must be 2'
         print *,'Error: Second array dimensions must be 1'
         print *,'First array has ',nindex_source1,' indices.'
         print *,'Second array has ',nindex_source2,' indices.'
         call abort_job()
      endif

c-------------------------------------------------------------------------
c   Get segment ranges of the array source.
c-------------------------------------------------------------------------

      do i = 1, nindex_source1
         ind(i)  = array_table(c_index_array1+i-1,source1)
         inda(i) = index_table(c_index_size,ind(i))
      enddo
      isource1 = get_index_from_base(address_table(source1), xxx, 2)

c-------------------------------------------------------------------------
c   Get segment ranges of the array target.
c-------------------------------------------------------------------------

      do i = 1, nindex_source2
         ind(i)  = array_table(c_index_array1+i-1,source2)
         indb(i) = index_table(c_index_size, ind(i))  ! pick up length of index
      enddo
      isource2 = get_index_from_base(address_table(source2), xxx, 2)

c--------------------------------------------------------------------------

      call mult_212(xxx(isource1),xxx(isource2),inda(1),
     *                    inda(2),indb(1))

      return
      end
C-------------------------------------------------------------------------

      subroutine mult_212(xxx,yyy,pa1,pa2,pb1)

      implicit none
      include 'interpreter.h'
      include 'int_gen_parms.h'
      include 'mpif.h'
      include 'proto_events.h'
      include 'parallel_info.h'
C
      integer pa1,pa2,pb1,pb2,atype,btype,i,j,ierr,n,k
      integer company_comm, pst_get_company_comm,master,pst_get_master
      double precision xxx(1:pa1,1:pa2)
      double precision yyy(1:pb1)
      double precision z
C
#ifdef _DEBUG_LVL0
      if (me .eq. master) then
      write(6,"(a)") "Checking the input arrays"
      Write(*,"(a,1x,i4,1x,i4)") "The input dimensions: ", pb1, pa2
      call output(xxx, 1, pa1, 1, pa2, pa1, pa2,1)
      write(*,*)
      write(*,"(a)") "The occupation vector"
      Write(*,"(6(1x,F3.1))") (yyy(i),i=1,pa1)
      write(6,"(a)") "End checking input array"
      endif 
#endif 
   
      Do j = 1, pa2
         Z = yyy(j) 
         Do i = 1, pa1
            xxx(i,j) = xxx(i,j) * Z
         Enddo 

      Enddo
#ifdef _DEBUG_LVL0
      if (me .eq. master) then
      write(*,*)
      write(*,"(a)") "The scalled array"
      call output(xxx, 1, pa1, 1, pa2, pa1, pa2, 1)
      Write(*,*)
      endif
#endif


      return
      end
