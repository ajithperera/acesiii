      Subroutine Set_lmn_delta(LMNVAL,LVAL)

      Implicit none
      integer MAXANG
      parameter (Maxang = 7)
      integer LMNval(3, (MAXANG+1)*(MAXANG+2)/2)
      integer ii,lval,m,left,l,n,i,i1,NXYZ
C
      II = 0
      DO L = LVAL,0,-1
         LEFT = LVAL - L
         DO M = LEFT,0,-1
            N = LEFT - M
            II = II + 1
            LMNVAL(1,II) = L
            LMNVAL(2,II) = M
            LMNVAL(3,II) = N
         ENDDO
      ENDDO
C
#ifdef _DEBUG_LVL0
      do i=1,3
         Write(6,"(6(1x,I3))") (lmnval(i, i1),i1=1,(Lval+1)*(Lval+2)/2)
      enddo
      Write(6,*)
#endif

      Return
      End
