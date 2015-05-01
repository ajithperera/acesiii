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
         SUBROUTINE  ERD__2D_ATOM_DSHIELD_PQ_INTEGRALS
     +
     +                    ( MIJ,MKL,MIJKL,NGQP,NGQEXQ,
     +                      P,Q,PINVHF,QINVHF,PQPINV,
     +                      RTS,B00,B01,B10,
     +                      SHELLP,SHELLQ,k,
     +                      WTS,
     +                      XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ,
     +                      CASE2D_0,CASE2D_1,
     + 
     +                               INT2DX,
     +                               INT2DY,
     +                               INT2DZ,
     + 
     +                               INT2DX1,
     +                               INT2DY1,
     +                               INT2DZ1 )
     +
C------------------------------------------------------------------------
C  OPERATION   : ERD__2D_ATOM__DSHIELD_PQ_INTEGRALS
C  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT
C  MODULE-ID   : ERD
C  SUBROUTINES : none
C  DESCRIPTION : This operation calculates a full table of 2D PQ X,Y,Z
C                atomic integrals using the Rys vertical recurrence
C                scheme VRR explained below. The present atomic case
C                differs from the general nonatomic case in the sense
C                that only B00,B01 and B10 coefficients are required
C                (hence no C00x and D00x with x=X,Y,Z coefficients
C                need to be ever calculated) and only nonzero 2D PQ
C                integrals arise in case the sum of both indices
C                corresponding to P and Q is even.
C
C                The Rys weight is multiplied to the 2DX PQ integral
C                to reduce overall FLOP count. Note, that the Rys weight
C                factor needs to be introduced only two times for
C                the starting atomic 2DX PQ integrals for the recurrence
C                scheme, namely to the (0,0) and (1,1) elements. The
C                weight factor is then automatically propagated
C                through the vertical transfer equations.
C
C                The recurrence scheme VRR is due to Rys, Dupuis and
C                King, J. Comp. Chem. 4, p.154-157 (1983) and details
C                about the scheme can be found in the general nonatomic
C                2D PQ X,Y,Z integral generation routine.
C
C                The 2D PQ integrals are calculated for all roots (info
C                already present in transmitted VRR coefficients!) and
C                for all exponent quadruples simultaneously and placed
C                into a 3-dimensional array.
C
C
C                  Input:
C
C                    SHELLx      =  maximum shell type for electrons
C                                   1 and 2 (x = P,Q)
C                    NGQEXQ      =  product of # of gaussian quadrature
C                                   points times exponent quadruplets
C                    WTS         =  all quadrature weights
C                    B00,B01,B10 =  cartesian coordinate independent
C                                   VRR expansion coefficients
C
C
C                  Output:
C
C                    INT2Dx      =  all 2D PQ atomic integrals for each
C                                   cartesian component (x = X,Y,Z)
C                    INT2Dx1     =  all 2D^(1) PQ atomic integrals for each
C                                   cartesian component (x = X,Y,Z)
C
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         IMPLICIT  NONE

         INTEGER   CASE2D_0,CASE2D_1
         INTEGER   I,K,N
         INTEGER   I1,I2,K1,K2,KP1
         INTEGER   MIJ,MKL,MIJKL,NGQP,NGQEXQ
         INTEGER   SHELLP,SHELLQ
         INTEGER   COMPONENT  
         INTEGER   XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ

         DOUBLE PRECISION  B0,B1
         DOUBLE PRECISION  ONE
         DOUBLE PRECISION  WEIGHT
         DOUBLE PRECISION  XI,XK,XI1,XK1,XKP1

         DOUBLE PRECISION  P       (1:MIJ)
         DOUBLE PRECISION  Q       (1:MKL)
         DOUBLE PRECISION  PINVHF  (1:MIJ)
         DOUBLE PRECISION  QINVHF  (1:MKL)
         DOUBLE PRECISION  PQPINV  (1:MIJKL)

         DOUBLE PRECISION  B00 (1:NGQEXQ)
         DOUBLE PRECISION  B01 (1:NGQEXQ)
         DOUBLE PRECISION  B10 (1:NGQEXQ)
         DOUBLE PRECISION  WTS (1:NGQEXQ)
         DOUBLE PRECISION  RTS (1:NGQEXQ)

         DOUBLE PRECISION  INT2DX (1:NGQEXQ,0:SHELLP,0:SHELLQ+1)
         DOUBLE PRECISION  INT2DY (1:NGQEXQ,0:SHELLP,0:SHELLQ+1)
         DOUBLE PRECISION  INT2DZ (1:NGQEXQ,0:SHELLP,0:SHELLQ+1)

         DOUBLE PRECISION  INT2DX1 (1:NGQEXQ,0:SHELLP,0:SHELLQ)
         DOUBLE PRECISION  INT2DY1 (1:NGQEXQ,0:SHELLP,0:SHELLQ)
         DOUBLE PRECISION  INT2DZ1 (1:NGQEXQ,0:SHELLP,0:SHELLQ)

         PARAMETER  (ONE   =  1.D0)
C
C
C------------------------------------------------------------------------
C
C             ...jump according to the 4 different cases that can arise:
C
C         write(*,*) "XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ:", XX,XY,XZ,YX,YY,YZ,
C     +               ZX, ZY, ZZ

         IF      (XX .EQ. 1) THEN
              COMPONENT = 1
         ELSE IF (XY .EQ. 1) THEN
              COMPONENT = 2
         ELSE IF (XZ .EQ. 1) THEN
              COMPONENT = 3
         ELSE IF (YX .EQ. 1) THEN
              COMPONENT = 4
         ELSE IF (YY .EQ. 1) THEN
              COMPONENT = 5
         ELSE IF (YZ .EQ. 1) THEN
              COMPONENT = 6
         ELSE IF (ZX .EQ. 1) THEN
              COMPONENT = 7
         ELSE IF (ZY .EQ. 1) THEN
              COMPONENT = 8
         ELSE IF (ZZ .EQ. 1) THEN
              COMPONENT = 9
         ENDIF
C
C        Write(*,*) "ERD__2D_ATOM_COEFFICIENTS"
        CALL    ERD__2D_ATOM_COEFFICIENTS 
     +
     +                    ( MIJ,MKL,MIJKL,
     +                      NGQP,NGQEXQ,
     +                      P,
     +                      PINVHF,QINVHF,PQPINV,
     +                      RTS,
     +                      CASE2D_1,
     +
     +                                B00,B01,B10 )
     +
     +
C        Write(*,*) "ERD__2D_ATOM_PQ_QPLUS1_INTEGRALS:CASE2D",CASE2D_1

        CALL   ERD__2D_ATOM_PQ_QPLUS1_INTEGRALS
     +
     +                    ( SHELLP,SHELLQ,
     +                      NGQEXQ,
     +                      WTS,
     +                      B00,B01,B10,
     +                      CASE2D_1,
     +
     +                               INT2DX,
     +                               INT2DY,
     +                               INT2DZ )

C        Write(*,*) "printing zeroth order 2Ds"
C         do n=1, NGQP
C         do i=0,shellp
C         do k=0,shellq+1
C              write(*,"(3(1x,I3),F15.10)") n,i,k, INT2DX(n,i,k)
C              write(*,"(3(1x,I3),F15.10)") n,i,k, INT2DY(n,i,k)
C              write(*,"(3(1x,I3),F15.10)") n,i,k, INT2DZ(n,i,k)
C         enddo
C         enddo
C         enddo

C        Write(*,*) "ERD__2D_ATOM_PQ_QPLUS1_INTEGRALS",COMPONENT 
        GOTO (1,2,3,4,5,6,7,8,9) COMPONENT

    1     CALL    ERD__2D_ATOM_DSHIELDXX_PQ_INTEGRALS
     +
     +                    ( SHELLP,SHELLQ,
     +                      MIJ,MKL,NGQP,NGQEXQ,
     +                      P,Q,RTS,
     +                                INT2DX,
     +                                INT2DY,
     +                                INT2DZ,
     +
     +                                INT2DX1,
     +                                INT2DY1,
     +                                INT2DZ1 )


C         Write(*,*) "printing xx-1"
C         do n=1, NGQEXQ
C         do i=0,shellp
C         do k=0,shellq
C              write(*,"(3(1x,I3),F15.10)") n,i,k, INT2DX1(n,i,k)
C              write(*,"(3(1x,I3),F15.10)") n,i,k, INT2DY1(n,i,k)
C              write(*,"(3(1x,I3),F15.10)") n,i,k, INT2DZ1(n,i,k)
C         enddo
C         enddo
C         enddo

          RETURN 

    2     CALL    ERD__2D_ATOM_DSHIELDXY_PQ_INTEGRALS
     +
     +                    ( SHELLP,SHELLQ,
     +                      MIJ,MKL,NGQP,NGQEXQ,
     +                      P,Q,RTS,
     +                                INT2DX,
     +                                INT2DY,
     +                                INT2DZ,
     +
     +                                INT2DX1,
     +                                INT2DY1,
     +                                INT2DZ1 )


C         Write(*,*) "printing xy-1"
C         do n=1, NGQEXP
C         do i=0,shellp
C         do k=0,shellq
C              write(*,"(3(1x,I3),F15.10)") n,i,k, INT2DX1(n,i,k)
C              write(*,"(3(1x,I3),F15.10)") n,i,k, INT2DY1(n,i,k)
C              write(*,"(3(1x,I3),F15.10)") n,i,k, INT2DZ1(n,i,k)
C         enddo
C         enddo
C         enddo

          RETURN

    3     CALL    ERD__2D_ATOM_DSHIELDXZ_PQ_INTEGRALS
     +
     +                    ( SHELLP,SHELLQ,
     +                      MIJ,MKL,NGQP,NGQEXQ,
     +                      P,Q,RTS,
     +                                INT2DX,
     +                                INT2DY,
     +                                INT2DZ,
     +
     +                                INT2DX1,
     +                                INT2DY1,
     +                                INT2DZ1 )


C         Write(*,*) "printing xz-1"
C         do n=1, NGQEXQ
C         do i=0,shellp
C         do k=0,shellq
C              write(*,"(3(1x,I3),F15.10)") n,i,k, INT2DX1(n,i,k)
C              write(*,"(3(1x,I3),F15.10)") n,i,k, INT2DY1(n,i,k)
C              write(*,"(3(1x,I3),F15.10)") n,i,k, INT2DZ1(n,i,k)
C         enddo
C         enddo
C         enddo

          RETURN



   4      CALL    ERD__2D_ATOM_DSHIELDYX_PQ_INTEGRALS
     +
     +                    ( SHELLP,SHELLQ,
     +                      MIJ,MKL,NGQP,NGQEXQ,
     +                      P,Q,RTS,
     +                                INT2DX,
     +                                INT2DY,
     +                                INT2DZ,
     +
     +                                INT2DX1,
     +                                INT2DY1,
     +                                INT2DZ1 )
C         Write(*,*) "printing yx-1"
C         do n=1, NGQEXQ
C         do i=0,shellp
C         do k=0,shellq
C              write(*,"(3(1x,I3),F15.10)") n,i,k, INT2DX1(n,i,k)
C              write(*,"(3(1x,I3),F15.10)") n,i,k, INT2DY1(n,i,k)
C              write(*,"(3(1x,I3),F15.10)") n,i,k, INT2DZ1(n,i,k)
C         enddo
C         enddo
C         enddo
C
          RETURN 

    5     CALL    ERD__2D_ATOM_DSHIELDYY_PQ_INTEGRALS
     +
     +                    ( SHELLP,SHELLQ,
     +                      MIJ,MKL,NGQP,NGQEXQ,
     +                      P,Q,RTS,
     +                                INT2DX,
     +                                INT2DY,
     +                                INT2DZ,
     +
     +                                INT2DX1,
     +                                INT2DY1,
     +                                INT2DZ1 )
C        Write(*,*) "printing yy-1"
C         do n=1, NGQEXQ
C         do i=0,shellp
C         do k=0,shellq
C              write(*,"(3(1x,I3),F15.10)") n,i,k, INT2DX1(n,i,k)
C              write(*,"(3(1x,I3),F15.10)") n,i,k, INT2DY1(n,i,k)
C              write(*,"(3(1x,I3),F15.10)") n,i,k, INT2DZ1(n,i,k)
C         enddo
C         enddo
C         enddo

         RETURN

    6     CALL    ERD__2D_ATOM_DSHIELDYZ_PQ_INTEGRALS
     +
     +                    ( SHELLP,SHELLQ,
     +                      MIJ,MKL,NGQP,NGQEXQ,
     +                      P,Q,RTS,
     +                                INT2DX,
     +                                INT2DY,
     +                                INT2DZ,
     +
     +                                INT2DX1,
     +                                INT2DY1,
     +                                INT2DZ1 )
C        Write(*,*) "printing yZ-1"
C         do n=1, NGQEXQ
C         do i=0,shellp
C         do k=0,shellq
C              write(*,"(3(1x,I3),F15.10)") n,i,k, INT2DX1(n,i,k)
C              write(*,"(3(1x,I3),F15.10)") n,i,k, INT2DY1(n,i,k)
C              write(*,"(3(1x,I3),F15.10)") n,i,k, INT2DZ1(n,i,k)
C         enddo
C         enddo
C         enddo

         RETURN


    7    CALL    ERD__2D_ATOM_DSHIELDZX_PQ_INTEGRALS
     +
     +                    ( SHELLP,SHELLQ,
     +                      MIJ,MKL,NGQP,NGQEXQ,
     +                      P,Q,RTS,
     +                                INT2DX,
     +                                INT2DY,
     +                                INT2DZ,
     +
     +                                INT2DX1,
     +                                INT2DY1,
     +                                INT2DZ1 )

C         Write(*,*) "printing zx-1"
C         do n=1, NGQEXQ
C         do i=0,shellp
C         do k=0,shellq
C              write(*,"(3(1x,I3),F15.10)") n,i,k, INT2DX1(n,i,k)
C              write(*,"(3(1x,I3),F15.10)") n,i,k, INT2DY1(n,i,k)
C              write(*,"(3(1x,I3),F15.10)") n,i,k, INT2DZ1(n,i,k)
C         enddo
C         enddo
C         enddo

         RETURN


    8     CALL    ERD__2D_ATOM_DSHIELDZY_PQ_INTEGRALS
     + 
     +                    ( SHELLP,SHELLQ,
     +                      MIJ,MKL,NGQP,NGQEXQ,
     +                      P,Q,RTS,
     +                                INT2DX,
     +                                INT2DY,
     +                                INT2DZ,
     +
     +                                INT2DX1,
     +                                INT2DY1,
     +                                INT2DZ1 )
C        Write(*,*) "printing zy-1"
C         do n=1, NGQEXQ
C         do i=0,shellp
C         do k=0,shellq
C              write(*,"(3(1x,I3),F15.10)") n,i,k, INT2DX1(n,i,k)
C              write(*,"(3(1x,I3),F15.10)") n,i,k, INT2DY1(n,i,k)
C              write(*,"(3(1x,I3),F15.10)") n,i,k, INT2DZ1(n,i,k)
C         enddo
C         enddo
C         enddo
C

          RETURN

    9     CALL    ERD__2D_ATOM_DSHIELDZZ_PQ_INTEGRALS
     +
     +                    ( SHELLP,SHELLQ,
     +                      MIJ,MKL,NGQP,NGQEXQ,
     +                      P,Q,RTS,
     +                                INT2DX,
     +                                INT2DY,
     +                                INT2DZ,
     +
     +                                INT2DX1,
     +                                INT2DY1,
     +                                INT2DZ1 )
C       Write(*,*) "printing zz-1"
C        do n=1, NGQEXQ
C        do i=0,shellp
C        do k=0,shellq
C             write(*,"(3(1x,I3),F15.10)") n,i,k, INT2DX1(n,i,k)
C             write(*,"(3(1x,I3),F15.10)") n,i,k, INT2DY1(n,i,k)
C             write(*,"(3(1x,I3),F15.10)") n,i,k, INT2DZ1(n,i,k)
C        enddo
C        enddo
C        enddo

          RETURN
C
C
C
C             ...ready!
C
C
         RETURN
         END
