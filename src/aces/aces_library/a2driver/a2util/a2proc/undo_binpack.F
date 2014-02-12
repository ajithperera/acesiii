      Subroutine Undo_binpack(ca, Nshells, Nprims_shell, 
     &                        Orig_nprim_shell, Reorder_Shell, 
     &                        Iscr, Dscr, Nrows, Ncolumns)

      Implicit Double Precision (A-H, O-Z)

      Dimension Ca(Nrows, Ncolumns), Nprims_shell(Nshells), 
     &          Orig_nprim_shell(Nshells)
      Integer   Reorder_Shell(Nshells)
      Dimension Iscr(Nrows), Dscr(Nrows)

#ifdef _DEBUG_LVL0
      Write(6,*) Nrows, Ncolumns,Nshells,(Nprims_shell(i), i=1,Nshells)
      Write(6,*) (reorder_shell(i), i=1,Nshells)
#endif

      if (Nrows .gt. 0) then

c-------------------------------------------------------------------------
c   Descramble the number of functions per shell.
c-------------------------------------------------------------------------

         do i = 1, nshells
            Orig_nprim_shell(Reorder_Shell(i)) = Nprims_shell(i)
         enddo

c--------------------------------------------------------------------------
c   Build a basis function index array from the shell index array.
c--------------------------------------------------------------------------

         istart = 1
         do i = 1, nshells
            ishell = Reorder_Shell(i)

c---------------------------------------------------------------------------
c   Calculate the starting basis function.
c---------------------------------------------------------------------------

            n = 0
            do j = 1, ishell-1
               n = n + Orig_nprim_shell(j)
            enddo

c----------------------------------------------------------------------------
c   Store the next "nfps_orig(ishell)" indices in the index array.
c----------------------------------------------------------------------------

            do j = 1, Orig_nprim_shell(ishell)
               iScr(istart) = n + j
               istart = istart + 1
            enddo
         enddo

         do j = 1, Ncolumns 

c--------------------------------------------------------------------------
c   Save column "j".
c--------------------------------------------------------------------------

            do i = 1, Nrows
               Dscr(i) = ca(i,j)
            enddo

            do i = 1, Nrows
               ca(iScr(i),j) = Dscr(i)
            enddo

         enddo

      Endif 

      return
      end
