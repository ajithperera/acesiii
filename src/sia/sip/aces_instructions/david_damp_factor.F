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
      subroutine David_damp_factor(array_table, narray_table, 
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
      include 'trace.h'
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
      integer ierr,sind
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

      double precision xxx(1), damp
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

      atype = array_table(c_array_type, source1)

      If (atype .NE. static_array) then
          print *,'Error: Source Array must be static DAVID_DAMP_FACTOR'
         call abort_job()
      endif

      nindex_source1 = array_table(c_nindex, source1)

      if (nindex_source1 .ne. 1) then
         print *,"Error: In the routine DAVID_DAMP_FACTOR"
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
         print *,'Error: scalar not sent in DAVID_DAMP_FACTOR'
         print *,(op(i),i=1,loptable_entry)
         call abort_job()
      endif

      sind =  array_table(c_scalar_index, source2)
      if (sind .lt. 1 .or. sind .gt. nscalar_table) then
         print *,'Scalar table index out of range in ',
     *           ' DAVID_DAMP_FACTOR line',current_line
         print *,'Index for array ',source2,' is ',sind,' should be ',
     *           'between 1 and ',nscalar_table
         call abort_job()
      endif

      damp = scalar_table(sind) 
      call damp_factor(xxx(isource1), damp, inda(1))
      scalar_table(sind) = damp

      return
      end
C-------------------------------------------------------------------------

      subroutine damp_factor(xxx,damp, pa1)

      implicit none
      include 'interpreter.h'
      include 'int_gen_parms.h'
      include 'mpif.h'
      include 'proto_events.h'
      include 'parallel_info.h'
C
      integer pa1,pa2,pb1,pb2,atype,btype,i,j,ierr,n,k
      integer company_comm, pst_get_company_comm,master,pst_get_master

      double precision xxx(1:pa1)

      double precision DMPMAX, ZERO, PT25, PT5, TWO, FOUR, FAC, PT2
      double precision DE, DEP, DEAVG, DAMP

      DATA DMPMAX /256.0D+00/
      DATA ZERO,PT25,PT5,TWO,FOUR,FAC /0.0D+00,2.5D-01,0.5D+00,
     1 2.0D+00,4.0D+00,1.6D+01/

C#ifdef _DEBUG_LVL0
      if (me .eq. master) then
      write(6,"(a)") "Checking the input arrays"
      Write(*,"(a,1x,i4)") "The input dimensions: ", pa1
      write(6,"(4(1x,F10.5))") (XXX(i),i=1,pa1), DAMP
      endif 
C#endif 

      DE    = XXX(1)
      DEP   = XXX(2)
      DEAVG = XXX(3)

      IF( DE.GT.ZERO) GO TO 400
      IF(DEP.GT.ZERO) GO TO 300
      IF( DE.GT. DEP) GO TO 200
C
C     ----- DE < 0. , DEP < 0. , DE < DEP -----
C
      IF( ABS(DE).LT.TWO*DEAVG) GO TO 110
      DAMP=FAC* DMAX1(DAMP,DEAVG)
      GO TO 600
  110 IF( ABS(DE).GT.PT5*DEAVG) GO TO 600
      DAMP=DAMP/FAC
      GO TO 600
  200 CONTINUE
C
C     ----- DE < 0. , DEP < 0. , DE > DEP -----
C
      IF(DE.GT.PT25*DEP) GO TO 210
      DAMP=(DE/DEP)**2* DMAX1(DAMP,DEAVG)
      GO TO 600
  210 DAMP=DAMP/FAC
      GO TO 600
  300 CONTINUE
C
C     ----- DE < 0. , DEP > 0. -----
C
      DAMP=FOUR* DMAX1(DAMP,DEAVG)
      IF(-DE.GT.DEAVG) DAMP=DAMP*FAC
      IF(-DE+DEP.GE.DEAVG) GO TO 600
      DAMP=DAMP/FAC
      GO TO 600
  400 CONTINUE
      IF(DEP.GT.ZERO) GO TO 500
C
C     ----- DE > 0. , DEP < 0. -----
C
      DAMP=FOUR* DMAX1(DAMP,DEAVG)
      IF(DE.GT.PT5*DEAVG) DAMP=DAMP*FAC
      IF(DE-DEP.GE.PT2*DEAVG) GO TO 600
      DAMP=DAMP/FAC
      GO TO 600
  500 CONTINUE
C
C     ----- DE > 0. , DEP > 0. -----
C
      DAMP=FOUR* DMAX1(DAMP,DEAVG)
      IF(DE.LT.FOUR*DEP) GO TO 510
      DAMP=FAC* DMAX1(DAMP,DEAVG)
      GO TO 600
  510 IF(DE.GT.PT25*DEP) GO TO 520
      DAMP=DAMP/FAC
      GO TO 600
  520 DAMP=(DE/DEP)**2* DMAX1(DAMP,DEAVG)
C
C     ----- DO NOT LET -DAMP- GO OUT OF CONTROL -----
C
  600 CONTINUE
      DAMP=DMIN1(DAMP,DMPMAX)
   
C#ifdef _DEBUG_LVL0
      if (me .eq. master) then
      write(*,"(a,1x, F10.5)") "The damp factor", damp
      endif
C#endif

      return
      end
