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
         SUBROUTINE  ERD__2D_DSHIELD_PQ_INTEGRALS
     +
     +                    ( MIJ,MKL,MIJKL,NGQP,NGQEXQ,
     +                      ATOMAB,ATOMCD,P,Q,
     +                      PX,PY,PZ,QX,QY,QZ,
     +                      PAX,PAY,PAZ,QCX,QCY,QCZ,
     +                      PINVHF,QINVHF,PQPINV,
     +                      RTS,
     +                      SHELLP,SHELLQ,
     +                      WTS,
     +                      XA,YA,ZA,XC,YC,ZC,
     +                      ABX,ABY,ABZ,CDX,CDY,CDZ,
     +                      XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ,
     +                      B00,B01,B10,
     +                      C00X,C00Y,C00Z,
     +                      D00X,D00Y,D00Z,
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
C  OPERATION   : ERD__2D_DSHIELD_PQ_INTEGRALS
C  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT
C  MODULE-ID   : ERD
C  SUBROUTINES : none
C  DESCRIPTION : This operation calculates a full table of 2D PQ X,Y,Z
C                integrals using the Rys vertical recurrence scheme
C                VRR explained below.
C
C                The Rys weight is multiplied to the 2DX PQ integral
C                to reduce overall FLOP count. Note, that the Rys weight
C                factor needs to be introduced only three times for the
C                starting 2DX PQ integrals for the recurrence scheme,
C                namely to the (0,0), (1,0) and (0,1) elements. The
C                weight factor is then automatically propagated
C                through the vertical transfer equations (see below).
C                The recurrence scheme VRR is due to Rys, Dupuis and
C                King, J. Comp. Chem. 4, p.154-157 (1983).
C
C
C                   INT2D (0,0) = 1.D0    (* WEIGHT for the 2DX case)
C                   INT2D (1,0) = C00     (* WEIGHT for the 2DX case)
C                   INT2D (0,1) = D00     (* WEIGHT for the 2DX case)
C
C                   For I = 1,...,SHELLP-1
C                       INT2D (I+1,0) = I * B10 * INT2D (I-1,0)
C                                         + C00 * INT2D (I,0)
C                   For K = 1,...,SHELLQ-1
C                       INT2D (0,K+1) = K * B01 * INT2D (0,K-1)
C                                         + D00 * INT2D (0,K)
C                   For I = 1,...,SHELLP
C                       INT2D (I,1)   = I * B00 * INT2D (I-1,0)
C                                         + D00 * INT2D (I,0)
C                   For K = 2,...,SHELLQ
C                       INT2D (1,K)   = K * B00 * INT2D (0,K-1)
C                                         + C00 * INT2D (0,K)
C                   For K = 2,...,SHELLQ
C                   For I = 2,...,SHELLP
C                       INT2D (I,K)   = (I-1) * B10 * INT2D (I-2,K)
C                                         + K * B00 * INT2D (I-1,K-1)
C                                             + C00 * INT2D (I-1,K)
C
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
C                    XX,XY,XZ    =  identify which tensor component is
C                    YX,YY,YZ    =  is computed. The one that is being
C                    ZX,ZY,ZZ    =  computed carries 1 and the others
C                                    are zero

C                    NGQEXQ      =  product of # of gaussian quadrature
C                                   points times exponent quadruplets
C                    WTS         =  all quadrature weights
C                    B00,B01,B10 =  VRR expansion coefficients
C                                   (cartesian coordinate independent)
C                    C00x,D00x   =  cartesian coordinate dependent
C                                   VRR expansion coefficients
C                                   (x = X,Y,Z)
C
C
C                  Output:
C
C                    INT2Dx      =  all 2D PQ integrals for each
C                                   cartesian component (x = X,Y,Z)
C                    INT2Dx1     =  all 2D(^1) PQ integrals for each
C                                   cartesian component (x = X,Y,Z)
C
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C
C
C             ...include files and declare variables.
C
C
         IMPLICIT  NONE

         LOGICAL  ATOMAB,ATOMCD 

         INTEGER   CASE2D_0,CASE2D_1
         INTEGER   I,K,N
         INTEGER   I1,I2,K1,K2
         INTEGER   MIJ,MKL,MIJKL,NGQP
         INTEGER   NGQEXQ
         INTEGER   SHELLP,SHELLQ
         INTEGER   XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ
         INTEGER   COMPONENT 

         DOUBLE PRECISION  B0,B1
         DOUBLE PRECISION  F,F2
         DOUBLE PRECISION  ONE
         DOUBLE PRECISION  WEIHT

         DOUBLE PRECISION  P       (1:MIJ)
         DOUBLE PRECISION  PX      (1:MIJ)
         DOUBLE PRECISION  PY      (1:MIJ)
         DOUBLE PRECISION  PZ      (1:MIJ)
         DOUBLE PRECISION  PAX     (1:MIJ)
         DOUBLE PRECISION  PAY     (1:MIJ)
         DOUBLE PRECISION  PAZ     (1:MIJ)
         DOUBLE PRECISION  PINVHF  (1:MIJ)

         DOUBLE PRECISION  Q       (1:MKL)
         DOUBLE PRECISION  QX      (1:MKL)
         DOUBLE PRECISION  QY      (1:MKL)
         DOUBLE PRECISION  QZ      (1:MKL)
         DOUBLE PRECISION  QCX     (1:MKL)
         DOUBLE PRECISION  QCY     (1:MKL)
         DOUBLE PRECISION  QCZ     (1:MKL)
         DOUBLE PRECISION  QINVHF  (1:MKL)

         DOUBLE PRECISION  B00  (1:NGQEXQ)
         DOUBLE PRECISION  B01  (1:NGQEXQ)
         DOUBLE PRECISION  B10  (1:NGQEXQ)
         DOUBLE PRECISION  C00X (1:NGQEXQ)
         DOUBLE PRECISION  C00Y (1:NGQEXQ)
         DOUBLE PRECISION  C00Z (1:NGQEXQ)
         DOUBLE PRECISION  D00X (1:NGQEXQ)
         DOUBLE PRECISION  D00Y (1:NGQEXQ)
         DOUBLE PRECISION  D00Z (1:NGQEXQ)
         DOUBLE PRECISION  WTS  (1:NGQEXQ)
         DOUBLE PRECISION  RTS  (1:NGQEXQ)

         DOUBLE PRECISION  PQPINV  (1:MIJKL)
   
         DOUBLE PRECISION XA,YA,ZA,XC,YC,ZC
         DOUBLE PRECISION ABX,ABY,ABZ,CDX,CDY,CDZ

         DOUBLE PRECISION  INT2DX (1:NGQEXQ,0:SHELLP,0:SHELLQ+1)
         DOUBLE PRECISION  INT2DY (1:NGQEXQ,0:SHELLP,0:SHELLQ+1)
         DOUBLE PRECISION  INT2DZ (1:NGQEXQ,0:SHELLP,0:SHELLQ+1)

         DOUBLE PRECISION  INT2DX1 (1:NGQEXQ,0:SHELLP,0:SHELLQ)
         DOUBLE PRECISION  INT2DY1 (1:NGQEXQ,0:SHELLP,0:SHELLQ)
         DOUBLE PRECISION  INT2DZ1 (1:NGQEXQ,0:SHELLP,0:SHELLQ)
  
         DOUBLE PRECISION xzx,xzy,xzz,d0x,d0y,d0z

         PARAMETER  (ONE = 1.D0)
C
C-------------------------------------------------------------------
C

C         write(*,*) "XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ:", XX,XY,XZ,YX,YY,YZ,
C     +               ZX, ZY, ZZ

         xzx = 0.0d0
         xzy = 0.0d0
         xzz = 0.0d0

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
     
         CALL    ERD__2D_COEFFICIENTS
     +
     +                    ( MIJ,MKL,MIJKL,
     +                      NGQP,NGQEXQ,
     +                      ATOMAB,ATOMCD,
     +                      P,Q,
     +                      PX,PY,PZ,QX,QY,QZ,
     +                      PAX,PAY,PAZ,QCX,QCY,QCZ,
     +                      PINVHF,QINVHF,PQPINV,
     +                      RTS,
     +                      CASE2D_1,
     +
     +                                B00,B01,B10,
     +                                C00X,C00Y,C00Z,
     +                                D00X,D00Y,D00Z )
     +
     +
         CALL    ERD__2D_PQ_QPLUS1_INTEGRALS
     +
     +                    ( SHELLP,SHELLQ,
     +                      NGQEXQ,
     +                      WTS,
     +                      B00,B01,B10,
     +                      C00X,C00Y,C00Z,
     +                      D00X,D00Y,D00Z,
     +                      CASE2D_1,
     +
     +                                INT2DX,
     +                                INT2DY,
     +                                INT2DZ )
     +
     +
C         Write(*,*) "Zeroth order 2Ds"
C         do n=1, NGQEXQ
C         do i=0,shellp
C         do k=0,shellq+1
C               d0x = d0x + INT2DX(n,i,k) * INT2DX(n,i,k)
C               d0y = d0y + INT2DY(n,i,k) * INT2DY(n,i,k)
C               d0z = d0z + INT2DZ(n,i,k) * INT2DZ(n,i,k)
C
C              write(*,"(3(1x,I3),F15.10)") n,i,k, INT2DX(n,i,k)
C              write(*,"(3(1x,I3),F15.10)") n,i,k, INT2DY(n,i,k)
C              write(*,"(3(1x,I3),F15.10)") n,i,k, INT2DZ(n,i,k)
C         enddo
C         enddo
C         enddo
C         Write(*,"(F15.10)") d0x
C         Write(*,"(F15.10)") d0y
C         Write(*,"(F15.10)") d0z
C         Write(*,*) "ERD__2D_DSHIELD_PQ_INTEGRALS", COMPONENT 

         GOTO (1,2,3,4,5,6,7,8,9) COMPONENT 
C
   1     CALL   ERD__2D_DSHIELDXX_PQ_INTEGRALS
     +
     +                    ( SHELLP,SHELLQ,
     +                      MIJ,MKL,NGQP,NGQEXQ,
     +                      PX,PY,PZ,QX,QY,QZ,
     +                      XC,YC,ZC,
     +                      P,Q,RTS,
     +                      CASE2D_0,
     +                                INT2DX,
     +                                INT2DY,
     +                                INT2DZ,
     +
     +                                INT2DX1,
     +                                INT2DY1,
     +                                INT2DZ1 )
C        Write(*,*) "printing xx-1"
C         do n=1,  NGQEXQ
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

   2     CALL   ERD__2D_DSHIELDXY_PQ_INTEGRALS
     +
     +                    ( SHELLP,SHELLQ,
     +                      MIJ,MKL,NGQP,NGQEXQ,
     +                      PX,PY,PZ,QX,QY,QZ,
     +                      XC,YC,ZC,
     +                      P,Q,RTS,
     +                      CASE2D_0,
     +                                INT2DX,
     +                                INT2DY,
     +                                INT2DZ,
     +
     +                                INT2DX1,
     +                                INT2DY1,
     +                                INT2DZ1 )
C        Write(*,*) "printing xy-1"
C         do n=1,  NGQEXQ
C         do i=0,shellp`
C         do k=0,shellq
C              write(*,"(3(1x,I3),F15.10)") n,i,k, INT2DX1(n,i,k)
C              write(*,"(3(1x,I3),F15.10)") n,i,k, INT2DY1(n,i,k)
C              write(*,"(3(1x,I3),F15.10)") n,i,k, INT2DZ1(n,i,k)
C         enddo
C         enddo
C         enddo
C
         RETURN

   3     CALL   ERD__2D_DSHIELDXZ_PQ_INTEGRALS
     +
     +                    ( SHELLP,SHELLQ,
     +                      MIJ,MKL,NGQP,NGQEXQ,
     +                      PX,PY,PZ,QX,QY,QZ,
     +                      XC,YC,ZC,
     +                      P,Q,RTS,
     +                      CASE2D_0,
     +                                INT2DX,
     +                                INT2DY,
     +                                INT2DZ,
     +
     +                                INT2DX1,
     +                                INT2DY1,
     +                                INT2DZ1 )
C        Write(*,*) "printing xz-1", SHELLP, SHELLQ
C         do n=1,  NGQEXQ
C         do i=0,shellp
C         do k=0,shellq
C               xzx = xzx + INT2DX1(n,i,k) * INT2DX1(n,i,k)
C               xzy = xzy + INT2DY1(n,i,k) * INT2DY1(n,i,k)
C               xzz = xzz + INT2DZ1(n,i,k) * INT2DZ1(n,i,k)
C                
c              write(*,"(3(1x,I3),F15.10)") n,i,k, INT2DX1(n,i,k)
c              write(*,"(3(1x,I3),F15.10)") n,i,k, INT2DY1(n,i,k)
C              write(*,"(3(1x,I3),F15.10)") n,i,k, INT2DZ1(n,i,k)
C         enddo
C         enddo
C         enddo
C         Write(*,"(F15.10)") xzx
C         Write(*,"(F15.10)") xzy
C         Write(*,"(F15.10)") xzz

         RETURN



   4     CALL    ERD__2D_DSHIELDYX_PQ_INTEGRALS
     +
     +                    ( SHELLP,SHELLQ,
     +                      MIJ,MKL,NGQP,NGQEXQ,
     +                      PX,PY,PZ,QX,QY,QZ,
     +                      XC,YC,ZC,
     +                      P,Q,RTS,
     +                      CASE2D_0,
     +                                INT2DX,
     +                                INT2DY,
     +                                INT2DZ,
     +
     +                                INT2DX1,
     +                                INT2DY1,
     +                                INT2DZ1 )
C         Write(*,*) "printing yx-1"
C         do n=1, NGQP
C         do i=0,shellp
C         do k=0,shellq
C               xzx = xzx + INT2DX1(n,i,k) * INT2DX1(n,i,k)
C               xzy = xzy + INT2DY1(n,i,k) * INT2DY1(n,i,k)
C               xzz = xzz + INT2DZ1(n,i,k) * INT2DZ1(n,i,k)
C
C
C              write(*,"(3(1x,I3),F15.10)") n,i,k, INT2DX1(n,i,k)
C              write(*,"(3(1x,I3),F15.10)") n,i,k, INT2DY1(n,i,k)
C              write(*,"(3(1x,I3),F15.10)") n,i,k, INT2DZ1(n,i,k)
C         enddo
C         enddo
C         enddo
C         Write(*,"(F15.10)") xzx
C         Write(*,"(F15.10)") xzy
C         Write(*,"(F15.10)") xzz

          RETURN

   5     CALL    ERD__2D_DSHIELDYY_PQ_INTEGRALS
     +
     +                    ( SHELLP,SHELLQ,
     +                      MIJ,MKL,NGQP,NGQEXQ,
     +                      PX,PY,PZ,QX,QY,QZ,
     +                      XC,YC,ZC,
     +                      P,Q,RTS,
     +                      CASE2D_0,
     +                                INT2DX,
     +                                INT2DY,
     +                                INT2DZ,
     +
     +                                INT2DX1,
     +                                INT2DY1,
     +                                INT2DZ1 )
C        Write(*,*) "printing yy-1"
C         do n=1,  NGQEXQ
C         do i=0,shellp
C         do k=0,shellq
C              write(*,"(3(1x,I3),F15.10)") n,i,k, INT2DX1(n,i,k)
C              write(*,"(3(1x,I3),F15.10)") n,i,k, INT2DY1(n,i,k)
C              write(*,"(3(1x,I3),F15.10)") n,i,k, INT2DZ1(n,i,k)
C         enddo
C         enddo
C         enddo

         RETURN

   6     CALL    ERD__2D_DSHIELDYZ_PQ_INTEGRALS
     +
     +                    ( SHELLP,SHELLQ,
     +                      MIJ,MKL,NGQP,NGQEXQ,
     +                      PX,PY,PZ,QX,QY,QZ,
     +                      XC,YC,ZC,
     +                      P,Q,RTS,
     +                      CASE2D_0,
     +                                INT2DX,
     +                                INT2DY,
     +                                INT2DZ,
     +
     +                                INT2DX1,
     +                                INT2DY1,
     +                                INT2DZ1 )
C        Write(*,*) "printing yz-1"
C         do n=1,  NGQEXQ
C         do i=0,shellp
C         do k=0,shellq
C              write(*,"(3(1x,I3),F15.10)") n,i,k, INT2DX1(n,i,k)
C              write(*,"(3(1x,I3),F15.10)") n,i,k, INT2DY1(n,i,k)
C              write(*,"(3(1x,I3),F15.10)") n,i,k, INT2DZ1(n,i,k)
C         enddo
C         enddo
C         enddo


         RETURN

   7     CALL    ERD__2D_DSHIELDZX_PQ_INTEGRALS
     +
     +                    ( SHELLP,SHELLQ,
     +                      MIJ,MKL,NGQP,NGQEXQ,
     +                      PX,PY,PZ,QX,QY,QZ,
     +                      XC,YC,ZC,
     +                      P,Q,RTS,
     +                      CASE2D_0,
     +                                INT2DX,
     +                                INT2DY,
     +                                INT2DZ,
     +
     +                                INT2DX1,
     +                                INT2DY1,
     +                                INT2DZ1 )
C         Write(*,*) "printing zx-1"
C         do n=1, NGQP
C         do i=0,shellp
C         do k=0,shellq
C              write(*,"(3(1x,I3),F15.10)") n,i,k, INT2DX1(n,i,k)
C              write(*,"(3(1x,I3),F15.10)") n,i,k, INT2DY1(n,i,k)
C              write(*,"(3(1x,I3),F15.10)") n,i,k, INT2DZ1(n,i,k)
C         enddo
C         enddo
C         enddo

        RETURN


   8     CALL    ERD__2D_DSHIELDZY_PQ_INTEGRALS
     +
     +                    ( SHELLP,SHELLQ,
     +                      MIJ,MKL,NGQP,NGQEXQ,
     +                      PX,PY,PZ,QX,QY,QZ,
     +                      XC,YC,ZC,
     +                      P,Q,RTS,
     +                      CASE2D_0,
     +                                INT2DX,
     +                                INT2DY,
     +                                INT2DZ,
     +
     +                                INT2DX1,
     +                                INT2DY1,
     +                                INT2DZ1 )
C        Write(*,*) "printing zy-1"
C         do n=1, NGQP
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

   9     CALL    ERD__2D_DSHIELDZZ_PQ_INTEGRALS
     +
     +                    ( SHELLP,SHELLQ,
     +                      MIJ,MKL,NGQP,NGQEXQ,
     +                      PX,PY,PZ,QX,QY,QZ,
     +                      XC,YC,ZC,
     +                      P,Q,RTS,
     +                      CASE2D_0,
     +                                INT2DX,
     +                                INT2DY,
     +                                INT2DZ,
     +
     +                                INT2DX1,
     +                                INT2DY1,
     +                                INT2DZ1 )

C        Write(*,*) "printing zz-1"
C         do n=1,  NGQEXQ
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

C
C
C             ...ready!
C
C
         RETURN
         END
