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
         SUBROUTINE  ERD__2D_ATOM_DSHIELDXZ_PQ_INTEGRALS
     +
     +                    ( SHELLP,SHELLQ,
     +                      MIJ,MKL,NGQP,NGQEXQ,
     +                      P,Q,RTS,
     +                               INT2DX,
     +                               INT2DY,
     +                               INT2DZ,
     +
     +                               INT2DX1,
     +                               INT2DY1,
     +                               INT2DZ1 )
C------------------------------------------------------------------------
C  OPERATION   : ERD__2D_ATOM_DSHIELDXZ_PQ_INTEGRALS
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

         INTEGER   CASE2D
         INTEGER   I,K,N
         INTEGER   MIJ,MKL,IJ,KL,ROOT
         INTEGER   I1,I2,K1,K2,K3,KP1,KP2,KP3
         INTEGER   NGQP,NGQEXQ
         INTEGER   SHELLP,SHELLQ

         DOUBLE PRECISION  ONE,HALF
         DOUBLE PRECISION  PIJ,QKL
         DOUBLE PRECISION  P(1:MIJ),Q(1:MKL)
         DOUBLE PRECISION  QROOT
         DOUBLE PRECISION  PQ,PQSUM,TWORHO

         DOUBLE PRECISION  RTS  (1:NGQEXQ)

         DOUBLE PRECISION  INT2DX (1:NGQEXQ,0:SHELLP,0:SHELLQ+1)
         DOUBLE PRECISION  INT2DY (1:NGQEXQ,0:SHELLP,0:SHELLQ+1)
         DOUBLE PRECISION  INT2DZ (1:NGQEXQ,0:SHELLP,0:SHELLQ+1)

         DOUBLE PRECISION  INT2DX1 (1:NGQEXQ,0:SHELLP,0:SHELLQ)
         DOUBLE PRECISION  INT2DY1 (1:NGQEXQ,0:SHELLP,0:SHELLQ)
         DOUBLE PRECISION  INT2DZ1 (1:NGQEXQ,0:SHELLP,0:SHELLQ)

         LOGICAL DEBUG

         PARAMETER  (ONE   =  1.D0)
         PARAMETER  (HALF  =  0.5D0)
C
C
C------------------------------------------------------------------------
C
         N = 0

         DO IJ = 1,MIJ
            PIJ  = HALF / P(IJ)
            DO KL = 1,MKL
               QKL  = HALF / Q(KL)

               PQSUM  = PIJ + QKL
               TWORHO = ONE / PQSUM

               DO ROOT = 1, NGQP

                 N = N + 1
                 QROOT = RTS(N) * TWORHO

                 INT2DX1 (N,0,0) = 0.0D0
                 INT2DY1 (N,0,0) = INT2DY (N,0,0) * QROOT 
                 INT2DZ1 (N,0,0) = INT2DZ (N,0,1)

                 DO K = 1, SHELLQ
                    K1 = K - 1
                    K2 = K + 1
                    INT2DX1 (N,0,K) =-INT2DX (N,0,K1) * QKL * K
                    INT2DY1 (N,0,K) = INT2DY (N,0,K) * QROOT
                    INT2DZ1 (N,0,K) = INT2DZ (N,0,K2)
                 ENDDO

                 DO I = 1, SHELLP
                    I1 = I - 1
                    INT2DX1 (N,I,0) = INT2DX (N,I1,0) * PIJ * I
                    INT2DY1 (N,I,0) = INT2DY (N,I,0) * QROOT
                    INT2DZ1 (N,I,0) = INT2DZ (N,I,1)
                 ENDDO

                 IF (SHELLP .NE. 0 .AND. SHELLQ .NE. 0) THEN

                 DO I = 1, SHELLP
                    DO K = 1, SHELLQ
                       I1 = I - 1
                       K1 = K - 1
                       K2 = K + 1

                       INT2DX1 (N,I,K) = INT2DX (N,I1,K) * I * PIJ  -
     +                                   INT2DX (N,I,K1) * K * QKL
                       INT2DY1 (N,I,K) = INT2DY (N,I,K) * QROOT
                       INT2DZ1 (N,I,K) = INT2DZ (N,I,K2)

                     ENDDO
                 ENDDO

                 ENDIF 

               ENDDO
            ENDDO
         ENDDO
C
C
C             ...ready!
C
C
         RETURN
         END
