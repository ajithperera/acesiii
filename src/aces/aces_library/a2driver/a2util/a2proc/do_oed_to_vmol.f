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







      do i = 1, nrows 
         do j = 1, ncolumns 

             ii = erd_index(i) 

             Ovlp_vmol(ii,j) = ovlp_oed(i,j)/scalars(ii)







      enddo
      enddo






      return
      end
