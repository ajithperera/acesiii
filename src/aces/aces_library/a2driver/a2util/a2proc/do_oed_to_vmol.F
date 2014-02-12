      subroutine do_oed_to_vmol(nrows, ncolumns, erd_index, scalars, 
     &                          ovlp_oed,
     &                          ovlp_vmol)
c-----------------------------------------------------------------------------
c   Using the shell angular momentum and basis function information, this
c   subroutine calculates an array which maps the ACES integral order to 
c   the ERD integral order.  This array may be used to re-order a block 
c   of integrals calculated by the ERD package into a format corresponding
c   to the ACES (VMOL-based) integrals.
c
c-----------------------------------------------------------------------------
      implicit none

      integer erd_index(*)
      double precision scalars(*)
      double precision ovlp_vmol(nrows,ncolumns),
     &                 Ovlp_oed (nrows,ncolumns)

      integer i, j, k, ii, jj, kk  
      integer nrows, n, ncolumns 

#ifdef _DEBUG_LVLM
      Write(6,"(6(1x,i5))") (erd_index(i), i=1, nrows)
      Write(6,*)
      Write(6,"(6(1x,F10.5))") (scalars(i), i=1, nrows)
#endif

      do i = 1, nrows 
         do j = 1, ncolumns 

             ii = erd_index(i) 

             Ovlp_vmol(ii,j) = ovlp_oed(i,j)/scalars(ii)

#ifdef _DEBUG_LVLM1
           write(6,"(a,1x,i2,1x,i2,4(1x,F15.7))")
     &         ' ECP :', ii, j, Ovlp_oed(i,j), Ovlp_vmol(i,j),
     &                                scalars(ii)
#endif

      enddo
      enddo

#ifdef _DEBUG_LVLM
      write(6,"(a)") "The reordered and renormalized SCF vectors"
      call output(ovlp_vmol, 1, nrows, 1, ncolumns, nrows, 
     &            ncolumns,1)
#endif
      return
      end
