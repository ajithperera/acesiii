         SUBROUTINE  OED__delta_SET_IJ_PAIRS
     +
     +                    ( NPGTOA,NPGTOB,NPGTOAB,
     +                      ATOMIC,EQUALAB,
     +                      SWAPRS,
     +                      RNABSQ,
     +                      ALPHAA,ALPHAB,
     +
     +                               EMPTY,
     +                               NIJ,
     +                               PRIMA,PRIMB,
     +                               RHO )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__OVL_SET_IJ_PAIRS
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : none
C  DESCRIPTION : This operation determines the ij primitive exponent
C                pairs which will actually be considered for
C                contraction. It is here where the internal prescreening
C                (if any) of overlap primitives must be applied.
C
C                There are two ways the ij pairs can be ordered:
C                1) i runs fastest or 2) j runs fastest. The ordering
C                to be chosen depends on the order the contraction will
C                be performed. We have the following rule:
C
C                     If # primitives on A > # primitives on B
C                     (i.e. the SWAPRS = .true. case), the contraction
C                     will be performed first on the A primitives
C                     followed by the B primitives in order to save
C                     intermediate storage space for the 1st half
C                     transformation. Hence each B primitive will have
C                     associated a set of A primitives for the 1st
C                     half contraction and thus the order for the
C                     primitive pairs ij must be such that i runs
C                     fastest. If # primitives on A =< # primitives
C                     on B we have to let j run fastest.
C
C                The prescreening will be based on the fact that
C                for normalized cartesian integrals the primitive
C                (s|s) overlap integrals will be largest in magnitude
C                for all exponent combinations. The normalized (s|s)
C                overlap integral value expression is:
C
C                 Sqrt (8) * (ab)**(3/4) * (a+b)**(-3/2) * rho (a,b)
C
C                where a and b are the i-th and j-th gaussian exponents,
C                and rho (a,b) is the exponential prefactor given by the
C                expression:
C
C                                      - (ab/(a+b)) * R(AB)**2
C                       rho (a,b) = exp
C
C                where R(AB)**2 is the square of the distance between
C                the i-th and j-th atomic centers A and B.
C
C                The prescreening will be done in a 'K2' fashion,
C                running over all exponent pairs. An ij exponent pair
C                will be rejected, if the corresponding upper bound
C                of the (s|s) integral is less than a certain tolerance
C                value TOL.
C
C                Since evaluation of the exponential prefactors is
C                rather costly, the ones surviving the prescreening
C                are transmitted back to the calling routine for
C                further use.
C
C
C                  Input:
C
C                    NPGTOx       =  # of primitives per contraction
C                                    for contraction shells x = A,B
C                    NPGTOAB      =  # of primitive pairs for
C                                    contraction shell pair AB
C                    ATOMIC       =  indicates, if atomic integrals are
C                                    evaluated
C                    EQUALAB      =  indicates, if csh A and csh B are
C                                    considered to be equal
C                    SWAPRS       =  is .true. if the contraction order
C                                    of the primitive pairs AB will
C                                    be performed in reverse order
C                                    BA for efficiency reasons
C                    RNABSQ       =  the square distance R(AB)**2
C                                    between centers A and B
C                    ALPHAx       =  primitive exponents for centers
C                                    x=A,B
C
C                  Output:
C
C                    EMPTY        =  is true, if no ij pairs were
C                                    established
C                    NIJ          =  # of ij pairs after screening
C                    PRIMx        =  i,j labels of primitives for x=A,B
C                    RHO          =  NIJ exponential prefactors rho(a,b)
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

         LOGICAL   ATOMIC
         LOGICAL   EMPTY
         LOGICAL   EQUALAB
         LOGICAL   SWAPRS

         INTEGER   I,J
         INTEGER   NIJ
         INTEGER   NPGTOA,NPGTOB,NPGTOAB

         INTEGER   PRIMA (1:NPGTOAB)
         INTEGER   PRIMB (1:NPGTOAB)

         DOUBLE PRECISION  A,B,AB
         DOUBLE PRECISION  P
         DOUBLE PRECISION  PREFACT
         DOUBLE PRECISION  RHOAB
         DOUBLE PRECISION  RNABSQ
         DOUBLE PRECISION  SS
         DOUBLE PRECISION  ZP75,ONE,NEG1P5,TOL

         DOUBLE PRECISION  ALPHAA  (1:NPGTOA)
         DOUBLE PRECISION  ALPHAB  (1:NPGTOB)

         DOUBLE PRECISION  RHO (*)

         PARAMETER  (ZP75    =  0.75D0)
         PARAMETER  (ONE     =  1.D0)
         PARAMETER  (NEG1P5  = -1.5D0)
         PARAMETER  (TOL     =  1.D-14)
C
C
C------------------------------------------------------------------------
C
C
C         There is no screening for delta integral.

C         Write(*,*) "SWAPRS :", SWAPRS
         NIJ = 0
         IF (EQUALAB) THEN
            DO I = 1,NPGTOA
               DO J = 1,I
                    NIJ = NIJ + 1
                    PRIMA (NIJ) = I
                    PRIMB (NIJ) = J
                 END DO
            END DO
         ELSE
                 IF (SWAPRS) THEN
                     IF (ATOMIC) THEN
                         DO J = 1,NPGTOB
                         DO I = 1,NPGTOA
                            NIJ = NIJ + 1
                            PRIMA (NIJ) = I
                            PRIMB (NIJ) = J
                         END DO
                         END DO
                     ELSE
                         DO J = 1,NPGTOB
                            B = ALPHAB (J)
                            DO I = 1,NPGTOA
                               A = ALPHAA (I)
                               NIJ = NIJ + 1
                               PRIMA (NIJ) = I
                               PRIMB (NIJ) = J
                            END DO
                         END DO
                     END IF
                 ELSE
                     IF (ATOMIC) THEN
                         DO I = 1,NPGTOA
                         DO J = 1,NPGTOB
                            NIJ = NIJ + 1
                            PRIMA (NIJ) = I
                            PRIMB (NIJ) = J
                         END DO
                         END DO
                     ELSE
                         DO I = 1,NPGTOA
                            A = ALPHAA (I)
                            DO J = 1,NPGTOB
                               B = ALPHAB (J)
                               NIJ = NIJ + 1
                               PRIMA (NIJ) = I
                               PRIMB (NIJ) = J
                            END DO
                         END DO
                     END IF
                 END IF
         END IF

C
         IF (NIJ.EQ.0) THEN
             EMPTY = .TRUE.
         END IF
C
C             ...ready!
C
C
         RETURN
         END
