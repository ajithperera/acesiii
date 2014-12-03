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
      subroutine fock_denominator_dynamic(array_table, narray_table, 
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table, 
     *                      address_table, op)
c--------------------------------------------------------------------------
c   Divides each sample of the block V(a,i,b,j) given by the array argument 
c   by the denominator. Different from the Fock denominator in the sense 
C   that Fock eigenvalues do not need to be at fixed initial value.
c       eps = fock(i)+fock(j)-fock(a)-fock(b)
c--------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'trace.h'
      include 'mpif.h'
      include 'epsilon.h'
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
      integer nscalar_table
      double precision scalar_table(nscalar_table)
      integer*8 address_table(narray_table)

      integer i, j, k
      integer source, index, nindex, ierr
      integer source_type
      integer nindex_source
      integer block, blkndx, seg
      integer find_current_block
      integer*8 indblk, get_block_index
      integer*8 isource, itarget
      integer stack
      
      integer comm

      integer sval1(mx_array_index), sval2(mx_array_index)
      integer tval1(mx_array_index), tval2(mx_array_index)
      integer type(mx_array_index)
      integer n(mx_array_index),ind(mx_array_index)
      integer na1, na2, ni1, ni2, total 
      integer*8 addr, get_index_from_base

      double precision x(1)
#ifdef ALTIX
      pointer (dptr, x)
#else
      common x
#endif

      if (dryrun) return

#ifdef ALTIX
      dptr = dshptr
#endif

c---------------------------------------------------------------------------
c   Find the indices of the array block.
c---------------------------------------------------------------------------
       
      source = op(c_result_array)
c---------------------------------------------------------------------------
c   Look up source's address.
c---------------------------------------------------------------------------

      if (array_table(c_array_type,source) .eq. static_array) then
         addr = address_table(source)
         indblk = get_index_from_base(addr, x, 2)
      else
         block  = find_current_block(source, array_table(1,source),
     *                                index_table, nindex_table,
     *                                segment_table, nsegment_table,
     *                                block_map_table, blkndx)

         stack = array_table(c_array_stack,source)
         isource = get_block_index(source, block, stack, 
     *                              blkndx, x, .true.)
      endif 

c   Check the dimensions of arrays.
c--------------------------------------------------------------------------

      source_type = array_table(c_array_type, source)
      nindex_source = array_table(c_nindex, source)

      if (.not.(nindex_source .ne. 4 .or. nindex_source .ne. 2)) then
         print *,'Error: Source array dimension in'
         print *,' fock_denominator_dynamic must be 4.'
         print *,'Source array has ',nindex_source,' indices.'
         call abort_job()
      endif

c-------------------------------------------------------------------------
c   Get segment ranges of the array source.
c-------------------------------------------------------------------------

      do i = 1, nindex_source
         ind(i)  = array_table(c_index_array1+i-1,source)
         n(i)    = index_table(c_index_size, ind(i))  ! pick up length of index
         type(i) = index_table(c_index_type, ind(i))
         seg     = index_table(c_current_seg,ind(i))
         call get_index_segment(ind(i), seg, segment_table,
     *                          nsegment_table, index_table,
     *                          nindex_table, sval1(i), sval2(i))

      enddo
#ifdef _DEBUG_LVL0
      total = 26
      write(*,"(a)") "The Fock Diagonal arrays"
      write(*,"(6(1X,F12.7))") (epsilon(i), i=1, total)
      write(*,"(6(1X,F12.7))") (epsilonb(i), i=1, total)
      write(*,"(6(1X,i4))") (type(i), i=1,4)
#endif 

C-------------------------------------------------------------------------
      if (nindex_source .eq. 4) then

         call apply_fock_denominator4(x(isource),
     *                                sval1(1),sval2(1),type(1),
     *                                sval1(2),sval2(2),type(2),
     *                                sval1(3),sval2(3),type(3),
     *                                sval1(4),sval2(4),type(4),
     *                                epsilon,epsilonb)

      elseif (nindex_source .eq. 2) then

            if (array_table(c_array_type,source) .eq. static_array) then
                na1 = array_table(c_index_range1, source)
                na2 = array_table(c_index_range2,source)
                ni1 = array_table(c_index_range1+1,source)
                ni2 = array_table(c_index_range2+1,source)
            else
                na1 = sval1(1)
                na2 = sval2(1)
                ni1 = sval1(2)
                ni2 = sval2(2)
            endif 

         call apply_fock_denominator2(x(indblk),
     *                                na1,na2,type(1),
     *                                ni1,ni2,type(2),
     *                                epsilon,epsilonb)

      endif 
      
      return
      end
C-------------------------------------------------------------------------

      subroutine apply_fock_denominator4(x,
     *                                   a1,a2,atype,
     *                                   i1,i2,itype,
     *                                   b1,b2,btype,
     *                                   j1,j2,jtype,
     *                                   epsilon,epsilonb)
      implicit none
      include 'interpreter.h'

      integer a1,a2,b1,b2,i1,i2,j1,j2
      integer atype, itype, btype, jtype

      double precision x(a1:a2,i1:i2,b1:b2,j1:j2)
      double precision epsilon(*), epsilonb(*)

      integer a,b,i,j
      double precision eps, epsa, epsb, epsi, epsj
      double precision val

      do j = j1,j2
         if (jtype .eq. mobindex) then
            epsj = epsilonb(j)
         else
            epsj = epsilon(j)
         endif
      do b = b1,b2
         if (btype .eq. mobindex) then
            epsb = epsilonb(b)
         else
            epsb = epsilon(b)
         endif
      do i = i1, i2
         if (itype .eq. mobindex) then
            epsi = epsilonb(i)
         else
            epsi = epsilon(i)
         endif
      do a = a1,a2
         if (atype .eq. mobindex) then
            epsa = epsilonb(a)
         else
            epsa = epsilon(a)
         endif
         val =  x(a,i,b,j)
         eps = epsi + epsj - epsa - epsb
         if (eps .ne. 0.0) x(a,i,b,j) = val/eps
      enddo
      enddo
      enddo
      enddo
  
      return
      end

      subroutine apply_fock_denominator2(x,na1,na2,atype,
     *                       ni1,ni2,itype,epsilon, epsilonb)
      implicit none
      include 'interpreter.h'

      integer atype, itype
      integer na1, na2, ni1, ni2
      double precision x(na1:na2,ni1:ni2)

      double precision epsilon(*), epsilonb(*)

      integer a,i
      double precision eps, epsa, epsi
      double precision val

      if (atype .ne. itype) then
         write(6,*) ' Error in fock_divide2: atype .ne. itype '
      endif

      do i = ni1, ni2
         if (itype .eq. mobindex) then
            epsi = epsilonb(i)
         else
            epsi = epsilon(i)
         endif
      do a = na1,na2
         if (atype .eq. mobindex) then
            epsa = epsilonb(a)
         else
            epsa = epsilon(a)
         endif
         val =  x(a,i)
         eps = epsi - epsa
         if (eps .ne. 0.0) x(a,i) = val/eps
         if (eps .eq. 0.0) x(a,i) = 0.0d0
      enddo
      enddo

      return
      end

