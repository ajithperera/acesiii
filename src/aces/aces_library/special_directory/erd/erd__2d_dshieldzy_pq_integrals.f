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
         SUBROUTINE  ERD__2D_DSHIELDZY_PQ_INTEGRALS
     +
     +                    ( SHELLP,SHELLQ,
     +                      MIJ,MKL,NGQP,NGQEXQ,
     +                      PX,PY,PZ,QX,QY,QZ,
     +                      XC,YC,ZC,
     +                      P,Q,RTS,
     +                      CASE2D,
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
C  OPERATION   : ERD__2D_DSHIELDZY_PQ_INTEGRALS
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
C                    NGQEXQ      =  product of # of gaussian quadrature
C                                   points times exponent quadruplets
C                    WTS         =  all quadrature weights
C                    B00,B01,B10 =  VRR expansion coefficients
C                                   (cartesian coordinate independent)
C                    C00x,D00x   =  cartesian coordinate dependent
C                                   VRR expansion coefficients
C                                   (x = X,Y,Z)
C                    CASE2D      =  logical flag for simplifications
C                                   in 2D integral evaluation for
C                                   low quantum numbers
C
C
C                  Output:
C
C                    INT2Dx      =  all 2D PQ integrals for each
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

         INTEGER   CASE2D
         INTEGER   I,K,N
         INTEGER   IJ,KL
         INTEGER   I1,I2,K1,K2
         INTEGER   MIJ,MKL,NGQEXQ, NGQP
         INTEGER   SHELLP,SHELLQ,ROOT

         DOUBLE PRECISION  ONE,HALF,TWO
         DOUBLE PRECISION  P(1:MIJ),Q(1:MIJ)
         DOUBLE PRECISION  PX(1:MIJ),PY(1:MIJ),PZ(1:MIJ)
         DOUBLE PRECISION  QX(1:MIJ),QY(1:MIJ),QZ(1:MIJ)
         DOUBLE PRECISION  XC,YC,ZC
         DOUBLE PRECISION  PIJ,QKL
         DOUBLE PRECISION  PIJX,PIJY,PIJZ,QKLX,QKLY,QKLZ
         DOUBLE PRECISION  PQZ
         DOUBLE PRECISION  PQXC,PIJXC,QKLXC
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

         PARAMETER  (ONE  = 1.D0)
         PARAMETER  (HALF = 0.5D0) 
         PARAMETER  (TWO  = 2.0D0) 
C
C
C------------------------------------------------------------------------
C
C
C             ...jump according to the 4 different cases that can arise:
C
C                  P-shell = s- or higher angular momentum
C                  Q-shell = s- or higher angular momentum
C
C                each leading to simplifications in the VRR formulas.
C                The case present has been evaluated outside this
C                routine and is transmitted via argument.
C
C
         N = 0
         DEBUG = .TRUE.
         IF (DEBUG) THEN
         DO IJ = 1,MIJ
            PIJZ = PZ(IJ)
            PIJ  =  HALF / P(IJ)
            DO KL = 1,MKL
               QKLZ = QZ(KL)
               QKL  =  HALF / Q(KL)

               PQZ  =  PIJZ - QKLZ

               PQSUM  = PIJ + QKL
               TWORHO = ONE / PQSUM

               DO ROOT = 1, NGQP
                  N = N + 1
                  QROOT = RTS(N) * TWORHO

                  INT2DX1 (N,0,0) = INT2DX (N,0,0) * QROOT
                  INT2DY1 (N,0,0) = INT2DY (N,0,1) + YC *
     +                              INT2DY (N,0,0)
                  INT2DZ1 (N,0,0) = INT2DZ (N,0,0) * PQZ

                  DO K = 1, SHELLQ
                     K1 = K - 1
                     K2 = K + 1
                     INT2DX1 (N,0,K) = INT2DX (N,0,K)  * QROOT
                     INT2DY1 (N,0,K) = INT2DY (N,0,K2) + YC *
     +                                 INT2DY (N,0,K)
                     INT2DZ1 (N,0,K) = INT2DZ (N,0,K)  * PQZ -
     +                                 INT2DZ (N,0,K1) * QKL * K
                  ENDDO

                  DO I = 1, SHELLP
                     I1 = I - 1
                     INT2DX1 (N,I,0) = INT2DX (N,I,0)  * QROOT
                     INT2DY1 (N,I,0) = INT2DY (N,I,1)  + YC  *
     +                                 INT2DY (N,I,0)
                     INT2DZ1 (N,I,0) = INT2DZ (N,I,0)  * PQZ +
     +                                 INT2DZ (N,I1,0) * PIJ * I
                  ENDDO

                  IF (SHELLP .NE. 0 .AND. SHELLQ .NE. 0) THEN

                  DO I = 1, SHELLP
                     DO K = 1, SHELLQ
                        K1 = K - 1
                        I1 = I - 1
                        K2 = K + 1

                        INT2DX1 (N,I,K) = INT2DX (N,I,K)  * QROOT
                        INT2DY1 (N,I,K) = INT2DY (N,I,K2) + YC *
     +                                    INT2DY (N,I,K)
                        INT2DZ1 (N,I,K) = INT2DZ (N,I,K)  * PQZ      +
     +                                    INT2DZ (N,I1,K) * PIJ * I  -
     +                                    INT2DZ (N,I,K1) * QKL * K
                     ENDDO
                  ENDDO
                  
                  ENDIF

               ENDDO
            ENDDO
         ENDDO

         RETURN
         ENDIF

         GOTO (1,3,3,2,4,4,2,4,4) CASE2D
C
C
C             ...the case P = s-shell and Q = s-shell.
C
C
    1    DO IJ = 1,MIJ
            PIJZ = PZ(IJ)
            PIJ  = P(IJ)
            DO KL = 1,MKL
               QKLZ = QZ(KL)
               QKL  = Q(KL)

               PQZ  =  PIJZ - QKLZ

               PQ     = PIJ * QKL
               PQSUM  = PIJ + QKL
               TWORHO = TWO * PQ / PQSUM

               DO ROOT = 1, NGQP
                  N = N + 1
CSSS                  QROOT = RTS(N) * (1.0D0 - RTS(N)) * TWORHO
                  QROOT = RTS(N) * TWORHO

                  INT2DX1 (N,0,0) = INT2DX (N,0,0) * QROOT 

                  INT2DY1 (N,0,0) = INT2DY (N,0,1) + YC *
     +                              INT2DY (N,0,0)

                  INT2DZ1 (N,0,0) = INT2DZ (N,0,0) * PQZ

               ENDDO
            ENDDO
         ENDDO

         RETURN
C
C             ...the cases P = s-shell and Q >= p-shell.
C                Evaluate I=0 and K=0,1.
C
C
    2    DO IJ = 1,MIJ
            PIJZ = PZ(IJ)
            PIJ  = HALF / P(IJ)
                
            DO KL = 1,MKL
               QKLZ = QZ(KL)
               QKL  = HALF / Q(KL)

               PQZ  =  PIJZ - QKLZ
           
               PQSUM  = PIJ + QKL
               TWORHO = ONE / PQSUM

               DO ROOT = 1, NGQP
                  N = N + 1
CSSS                  QROOT = RTS(N) * (1.0D0 - RTS(N)) * TWORHO
                  QROOT = RTS(N) * TWORHO


                  INT2DX1 (N,0,0) =  INT2DX (N,0,0) * QROOT
                  INT2DX1 (N,0,1) =  INT2DX (N,0,1) * QROOT 

                  INT2DY1 (N,0,0) =  INT2DY (N,0,1) + YC *
     +                               INT2DY (N,0,0)
                  INT2DY1 (N,0,1) =  INT2DY (N,0,2) + YC *
     +                               INT2DY (N,0,1)

                  INT2DZ1 (N,0,0) =  INT2DZ (N,0,0) * PQZ
                  INT2DZ1 (N,0,1) =  INT2DZ (N,0,1) * PQZ -
     +                               INT2DZ (N,0,0) * QKL

                 DO K = 2, SHELLQ
                     K1 = K + 1
                     K2 = K - 1

                     INT2DX1 (N,0,K) = INT2DX (N,0,K)  * QROOT 

                     INT2DY1 (N,0,K) = INT2DY (N,0,K1) + YC  *
     +                                 INT2DY (N,0,K) 

                     INT2DZ1 (N,0,K) = INT2DZ (N,0,K)  * PQZ -
     +                                 INT2DZ (N,0,K2) * QKL * K
                 ENDDO

               ENDDO
            ENDDO
         ENDDO

         RETURN
C
C
C             ...the cases P >= p-shell and Q = s-shell.
C                Evaluate I=0,1 and K=0.
C
C

    3     DO IJ = 1,MIJ
            PIJZ = PZ(IJ)
            PIJ  = HALF / P(IJ)

            DO KL = 1,MKL
               QKLZ = QZ(KL)
               QKL  = HALF / Q(KL)

               PQZ  =  PIJZ - QKLZ

               PQSUM  = PIJ + QKL
               TWORHO = ONE / PQSUM

               DO ROOT = 1, NGQP

                  N = N + 1
CSSS                  QROOT = RTS(N) * (1.0D0 - RTS(N)) * TWORHO
                  QROOT = RTS(N) * TWORHO


                  INT2DX1 (N,0,0) = INT2DX (N,0,0) * QROOT 
                  INT2DX1 (N,1,0) = INT2DX (N,1,0) * QROOT 

                  INT2DY1 (N,0,0) = INT2DY (N,0,1) +
     +                              INT2DY (N,0,0) * YC
                  INT2DY1 (N,1,0) = INT2DY (N,1,1) +
     +                              INT2DY (N,1,0) * YC

                  INT2DZ1 (N,0,0) = INT2DZ (N,0,0) * PQZ
                  INT2DZ1 (N,1,0) = INT2DZ (N,1,0) * PQZ + 
     +                              INT2DZ (N,0,0) * PIJ
                  DO I = 2, SHELLP
                     I1 = I - 1

                     INT2DX1 (N,I,0) = INT2DX (N,I,0)  * QROOT 

                     INT2DY1 (N,I,0) = INT2DY (N,I,1)  + 
     +                                 INT2DY (N,I,0)  * YC

                     INT2DZ1 (N,I,0) = INT2DZ (N,I,0)  * PQZ +
     +                                 INT2DZ (N,I1,0) * PIJ * I
                  ENDDO

               ENDDO
            ENDDO
         ENDDO

         RETURN
C
C             ...the cases P >= p-shell and Q >= p-shell.
C                Evaluate I=0,SHELLP       I=0
C                         K=0        and   K=0,SHELLQ
C

    4    DO IJ = 1,MIJ
            PIJZ = PZ(IJ)
            PIJ  = HALF / P(IJ)

            DO KL = 1,MKL
               QKLZ = QZ(KL)
               QKL  = HALF / Q(KL)

               PQZ =  PIJZ - QKLZ

               PQSUM  = PIJ + QKL
               TWORHO = ONE / PQSUM

               DO ROOT = 1, NGQP

                  N = N + 1
CSSS                  QROOT = RTS(N) * (1.0D0 - RTS(N)) * TWORHO
                  QROOT = RTS(N) * TWORHO

                  INT2DX1 (N,0,0) = INT2DX (N,0,0) * QROOT 
                  INT2DX1 (N,1,0) = INT2DX (N,1,0) * QROOT
                  INT2DX1 (N,0,1) = INT2DX (N,0,1) * QROOT

                  INT2DY1 (N,0,0) = INT2DY (N,0,1) +
     +                              INT2DY (N,0,0) * YC
                  INT2DY1 (N,1,0) = INT2DY (N,1,1) +
     +                              INT2DY (N,1,0) * YC
                  INT2DY1 (N,0,1) = INT2DY (N,0,2) +
     +                              INT2DY (N,0,1) * YC

                  INT2DZ1 (N,0,0) = INT2DZ (N,0,0) * PQZ
                  INT2DZ1 (N,1,0) = INT2DZ (N,1,0) * PQZ +
     +                              INT2DZ (N,0,0) * PIJ
                  INT2DZ1 (N,0,1) = INT2DZ (N,0,1) * PQZ -
     +                              INT2DZ (N,0,0) * QKL

                  DO I = 2, SHELLP
                     I1 = I - 1

                     INT2DX1 (N,I,0) = INT2DX (N,I,0)  * QROOT 

                     INT2DY1 (N,I,0) = INT2DY (N,I,1)  +
     +                                 INT2DY (N,I,0)  * YC

                     INT2DZ1 (N,I,0) = INT2DZ (N,I,0)  * PQZ +
     +                                 INT2DZ (N,I1,0) * PIJ * I
                  ENDDO

                  DO K = 2, SHELLQ
                     K1 = K + 1
                     K2 = K - 1

                     INT2DX1 (N,0,K) = INT2DX (N,0,K)   * QROOT

                     INT2DY1 (N,0,K) = INT2DY (N,0,K1)  +
     +                                 INT2DY (N,0,K)   * YC 

                     INT2DZ1 (N,0,K) = INT2DZ (N,0,K)   * PQZ -
     +                                 INT2DZ (N,0,K2)  * QKL * K

                  ENDDO

C             ...evaluate I=1,SHELLP and K=1,SHELLQ (if any)
C                in most economical way.
C
C
                 IF (SHELLQ.LE.SHELLP) THEN

                    DO K = 1,SHELLQ
                       K1 = K + 1 
                       K2 = K - 1

                       INT2DX1 (N,1,K) = INT2DX (N,1,K)  * QROOT

                       INT2DY1 (N,1,K) = INT2DY (N,1,K1) +
     +                                   INT2DY (N,1,K)  * YC

                       INT2DZ1 (N,1,K) = INT2DZ (N,1,K)  * PQZ +
     +                                   INT2DZ (N,0,K)  * PIJ -
     +                                   INT2DZ (N,1,K2) * QKL * K

                       DO I = 2,SHELLP
                          I1 = I - 1

                          INT2DX1 (N,I,K) = INT2DX (N,I,K)  * QROOT

                          INT2DY1 (N,I,K) = INT2DY (N,I,K1) +
     +                                      INT2DY (N,I,K)  * YC

                          INT2DZ1 (N,I,K) = INT2DZ (N,I,K)  * PQZ +
     +                                      INT2DZ (N,I1,K) * PIJ * I -
     +                                      INT2DZ (N,I,K2) * QKL * K

                       ENDDO
                    ENDDO

                 ELSE

                    DO I = 1,SHELLP
                       I1 = I - 1

                       INT2DX1 (N,I,1) =  INT2DX (N,I,1)  * QROOT 

                       INT2DY1 (N,I,1) =  INT2DY (N,I,2)  +
     +                                    INT2DY (N,I,1)  * YC

                       INT2DZ1 (N,I,1) =  INT2DZ (N,I,1)  * PQZ + 
     +                                    INT2DZ (N,I1,1) * PIJ * I -
     +                                    INT2DZ (N,I,0)  * QKL 
                       DO K = 2,SHELLQ
                          K1 = K + 1
                          K2 = K - 1

                          INT2DX1 (N,I,K) = INT2DX (N,I,K)  * QROOT

                          INT2DY1 (N,I,K) = INT2DY (N,I,K1) +
     +                                      INT2DY (N,I,K)  * YC

                          INT2DZ1 (N,I,K) = INT2DZ (N,I,K)  * PQZ + 
     +                                      INT2DZ (N,I1,K) * PIJ * I -
     +                                      INT2DZ (N,I,K2) * QKL * K 

                       ENDDO
                    ENDDO

                 END IF

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
