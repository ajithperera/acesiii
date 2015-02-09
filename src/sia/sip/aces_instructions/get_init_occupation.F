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
      subroutine get_init_occupation(array_table, narray_table, 
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
      include 'int_gen_parms.h'
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
      integer source1,source2,scalar
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
      atype   = array_table(c_array_type, source1)

      If (atype .NE. static_array) then
          print *,'Error: Arrays must be static in get_init_occupation'
         call abort_job()
      endif

      nindex_source1 = array_table(c_nindex, source1)

      if (nindex_source1 .ne. 1) then
         print *,'Error: Array dimensions must be 1'
         print *,'First array has ',nindex_source1,' indices.'
         call abort_job()
      endif

      source2 = op(c_op1_array)
      btype   = array_table(c_array_type, source2)

      if (btype .ne. scalar_value) then
         print *,'Error: scalar not sent in  get_init_occupation'
         print *,(op(i),i=1,loptable_entry)
         call abort_job()
      endif

      sind =  array_table(c_scalar_index, source2)
      if (sind .lt. 1 .or. sind .gt. nscalar_table) then
         print *,'Scalar table index out of range in ',
     *           ' get_init_occupation line ',current_line
         print *,'Index for array ',source2,' is ',sind,' should be ',
     *           'between 1 and ',nscalar_table
         call abort_job()
      endif
      scalar = INT (scalar_table(sind))

c-------------------------------------------------------------------------
c   Get segment ranges of the array source.
c-------------------------------------------------------------------------

      do i = 1, nindex_source1
         ind(i)  = array_table(c_index_array1+i-1,source1)
         inda(i) = index_table(c_index_size,ind(i))
      enddo
      isource1 = get_index_from_base(address_table(source1), xxx, 2)

      call built_occupation(xxx(isource1),inda(1),scalar)

      return
      end
C-------------------------------------------------------------------------

      subroutine built_occupation(xxx,pa1,uhf)

      implicit none
      include 'interpreter.h'
      include 'int_gen_parms.h'
      include 'mpif.h'
      include 'proto_events.h'
      include 'parallel_info.h'
      include 'epsilon.h'
C
      integer pa1,pa2,pb1,pb2,atype,btype,i,j,ierr,n,k,iunit
      integer itmp1(2),Itmp2(2),uhf
      integer company_comm, pst_get_company_comm,master,pst_get_master
      logical OCCNUMS_FILE
      character*80 BLANK
      double precision xxx(1:pa1)
      double precision DOCC_NUM
C
C Check whether OCCNUM file exist. If it is read the file as is done 
C in the serial code (but no symmetry)
C
      INQUIRE(FILE="OCCNUMS", EXIST=OCCNUMS_FILE)
      IUNIT = 5

      IF (OCCNUMS_FILE) THEN

C#ifdef _DEBUG_LVL0
      if (me .eq. master) then
      write(*,*)
      write(*,"(a)") "The occupation number file is present"
      write(*,"(a)") "The occupations are read from OCCNUMS file"
      write(*,*)
      endif
C#endif 
         OPEN(UNIT=IUNIT, FILE="OCCNUMS", FORM="FORMATTED")
         READ(IUNIT,10,END=19) Itmp1(1)
         READ(IUNIT,10,END=19) Itmp1(2)

         READ(IUNIT, "(80a)") Blank
        
         READ(IUNIT,10,END=19) Itmp2(1)
         READ(IUNIT,10,END=19) Itmp2(2)

         READ(IUNIT, "(80a)") Blank

         IF (UHF .EQ. 0) THEN
            DO I = 1, Itmp2(1)
               READ(IUNIT, "(F4.2)",END=19) DOCC_NUM
               XXX(i) = DOCC_NUM 
            ENDDO 

         ELSE

            DO I = 1, Itmp2(1)
               READ(IUNIT, "(F4.2)",END=19) DOCC_NUM
            ENDDO 

            READ(IUNIT, "(80a)") Blank

            DO I =1, Itmp2(2)
               READ(IUNIT, "(F4.2)",END=19) DOCC_NUM
               XXX(i) = DOCC_NUM
            ENDDO 

         ENDIF 

      ELSE
           
         IF (UHF .EQ. 0) THEN

            DO I = 1, nalpha_occupied
               XXX(i) = 1.0D0
            ENDDO 

         ELSE 
           
            DO I = 1, nbeta_occupied
               XXX(i) = 1.0D0
            ENDDO 

         ENDIF

      ENDIF 

 10   FORMAT(16I5)
 19   CLOSE(IUNIT)

C#ifdef _DEBUG_LVL0
      if (me .eq. master) then
      write(*,*)
      If (UHF .EQ. 0) THEN
      write(*,"(a)") "The Alpha occupation vector from init_occupy.F"
      Else
      write(*,"(a)") "The Beta occupation vector from init_occupy.F"
      Endif 
      write(*,*)
      Write(*,"(6(1x,F3.1))") (xxx(i),i=1,pa1)
      Write(*,*)
      endif
C#endif

      return
      end
