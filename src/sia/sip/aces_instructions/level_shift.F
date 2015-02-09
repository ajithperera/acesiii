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
      subroutine level_shift (array_table, narray_table, 
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table, 
     *                      address_table, op)
c--------------------------------------------------------------------------
C Augment the diagonals a 2D array by a constant.
c A(i,i) = A(i,i) + Shift where i can be any value bounded the dimension
C of the matrix.
c--------------------------------------------------------------------------

      implicit none
      include 'mpif.h'
      include 'interpreter.h'
      include 'int_gen_parms.h'
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
      integer nscalar_table
      double precision scalar_table(nscalar_table)
      double precision shift  
      integer sind, inda(2)
      integer*8 address_table(narray_table),get_index_from_base

      integer i, j, k
      integer source1,source2,atype,btype,nindex_source1,ierr
      integer block, blkndx, seg
      integer find_current_block
      integer*8 indblk,get_block_index 
      integer*8 isource1
      
      integer comm,n1,n2
      integer type(mx_array_index),ind(mx_array_index)

      double precision x(1)
#ifdef ALTIX
      pointer (dptr, x)
#else
      common x
#endif

#ifdef ALTIX
      dptr = dshptr
#endif

c---------------------------------------------------------------------------
c   Find the indices of the array block.
c---------------------------------------------------------------------------

      source1 = op(c_result_array)
      atype   = array_table(c_array_type, source1)

      if (atype .ne. static_array) then
        print *,'Error: The virt-virt Fock matrix array must ',
     *                   'be static in Level_shift SI.'
         call abort_job()
      endif

      nindex_source1 = array_table(c_nindex, source1)

      do i = 1, nindex_source1
         na      = array_table(c_index_array1+i-1,source1)
         inda(i) = index_table(c_index_size, na) 
         type(i) = index_table(c_index_type, na)
      enddo

      isource1 = get_index_from_base(address_table(source1), x, 2)

c---------------------------------------------------------------------------
c   Get the constant shift.
c---------------------------------------------------------------------------
       
      source2 = op(c_op1_array)
      btype   = array_table(c_array_type, source2)

      if (btype .ne. scalar_value) then 
         write(6,*) 'Error: second argument in level_shift SI ', 
     *               'must be a scalar. '
         call abort_job() 
      endif 

      sind  = array_table(c_scalar_index, source2)
      shift = scalar_table(sind) 

      If (nindex_source1 .eq.2) then

         n1 = Inda(1)
         n2 = Inda(2)

         call apply_shift (x(isource1),n1,n2,type(1),
     *                     type(2),shift)    
      else
         print *,'Error in : level_shift.F SIP:' , 
     *                       'must have 2 indices.'
         print *,'array ',source1,' has ',nindex_source1,' indices.'
         call abort_job()
      endif
      
      return
      end

      subroutine apply_shift (Fmat,na1,na2,ptype,qtype,shift)

      implicit none
      include 'interpreter.h'
      include 'int_gen_parms.h'
      Include 'parallel_info.h'

      integer ptype, qtype
      integer i1,i2
      integer na1, na2, Nocca,Noccb,Ntotala,Ntotalb
      double precision Fmat(na1,na2)
      double precision val, shift

      If (ptype .eq. moaindex .and. qtype .eq. moaindex) Then
          Nocca   = nalpha_occupied
          Ntotala = nalpha_occupied + nalpha_virtual 
      
#ifdef _DEBUG_LVL0
      If (me .eq.0) then
      Write(*,*) 
      Write(6,*) "The Alpha block of the Fock matrix"
      Write(6,"(6(1x,F15.8))") ((Fmat(i1,i2),i1=1, ntotala),
     *          i2=1,ntotala)
      Write(*,*)
      Endif 
#endif
         Do i1 = nalpha_occupied+1, Ntotala
            Fmat(i1,i1) = Fmat(i1,i1) + Shift
         Enddo 

      Else if (ptype .eq. mobindex .and. qtype .eq. mobindex) Then
          Noccb   = nbeta_occupied
          Ntotalb = nbeta_occupied + nbeta_virtual 

#ifdef _DEBUG_LVL0
      If (me .eq.0) then
      Write(*,*)
      Write(6,*) "The Beta block of the Fock matrix"
      Write(6,"(6(1x,F15.8))") ((Fmat(i1,i2),i1=1, ntotalb),
     *                           i2=1,ntotalb)
      Write(*,*)
      Endif
#endif
         Do i1 = nbeta_occupied+1, Ntotalb
            Fmat(i1,i1) = Fmat(i1,i1) + Shift
         Enddo 

      Endif

      return
      end
