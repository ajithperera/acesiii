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
      subroutine Report_maximum(array_table, narray_table, 
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

      integer eps, i, n, index4_max
      integer ierr
      integer source1,source2 
      integer array_type, inda(2)
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
          print *, " in Report_maximum.F"
         call abort_job()
      endif

      nindex_source1 = array_table(c_nindex, source1)
      nindex_source2 = array_table(c_nindex, source2)

      if (nindex_source1 .ne. 1) then
         print *,"Error: In the routine Report_maximum.F"
         print *,'Error: First array dimensions must be 1'
         print *,'First array has ',nindex_source1,' indices.'
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

      call find_maximum(xxx(isource1),index4_max,inda(1),inda(2))

      scalar_table(sind) = index4_max

      return
      end
C-------------------------------------------------------------------------

      subroutine find_maximum(xxx,index4_max,pa1,pa2)
      implicit none
      include 'interpreter.h'
      include 'int_gen_parms.h'
      include 'mpif.h'
      include 'proto_events.h'
      include 'parallel_info.h'
C
      integer pa1,pa2,index4_max,n
      double precision xxx(pa1:pa2),zq(6)
C
      Call Vstat_extension(xxx,zq,n,index4_max)
      xxx(1) = zq(1)

      return
      end

      SUBROUTINE VSTAT_EXTENSION(V,ZQ,LENGTH,INDEX4_MAX)
C
C RETURNS STATISTICAL INFO ABOUT VECTOR V in ZQ
C     ZQ(1)  Largest absolute magnitude
C     ZQ(2)  Smallest absolute magnitude
C     ZQ(3)  Largest value
C     ZQ(4)  Smallest value
C     ZQ(5)  2-norm
C     ZQ(6)  Dynamic range of the vector (abs. min. - abs. max.)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION V(LENGTH),ZQ(6)
      U=0.D0
      CALL ZERO(ZQ,6)
      ZQ(2)=DABS(V(1))
      ZQ(4)=V(1)
      DO 20 I=1,LENGTH
      ZQ(1)=MAX(ZQ(1),DABS(V(I)))
      ZQ(2)=MIN(ZQ(2),DABS(V(I)))
      ZQ(3)=MAX(ZQ(3),V(I))
      ZQ(4)=MIN(ZQ(4),V(I))
20    U=U+V(I)*V(I)
      If (Length .ne. 0) ZQ(5)=DSQRT(U/LENGTH)
      ZQ(6)=ZQ(2)-ZQ(1)
   
      DO I = 1, LENGTH
         IF (DABS(V(I) - ZQ(1)) .EQ. 0.0D0) Then
            index4_max = I
         ENDIF
      ENDDO 
      
      RETURN
      END
