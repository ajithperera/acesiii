
C THIS PROJECTS OUT LINEARLY DEPENDENT PARTS OF THE INTERNAL COORDINATES

c INPUT
c integer NATOM

c OUTPUT
c double SCR(*)
c double SYMQ(*)

c RECORDS
c get 'COMPNSYQ'
c get 'COMPSYMQ'
c put 'COMPSYMQ'

      SUBROUTINE SCHMIDT(NATOM,SCR,SYMQ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION SCR(*),SYMQ(*)

      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD

      DATA TOL /1.D-8/

      NSIZE=3*NATOM

      CALL GETREC(20,'JOBARC','COMPNSYQ',1,NSYMOLD)
      CALL GETREC(20,'JOBARC','COMPSYMQ',NSIZE*NSIZE*IINTFP,SCR)

c   o loop over vectors
      DO I=2,NSYMOLD
         IPOS=(I-1)*NSIZE+1
         CALL GSCHMIDT(SCR(IPOS),SCR,NSIZE,I-1,SYMQ,ZJUNK,1.D-8)
      END DO
      CALL DCOPY(NSIZE*NSIZE,SCR,1,SYMQ,1)

c   o renormalize
      IOFF=1
      DO I=1,NSYMOLD
         X=DNRM2(NSIZE,SYMQ(IOFF),1)
         IF (ABS(X).GT.TOL) THEN
            Z=1.d0/X
            CALL DSCAL(NSIZE,Z,SYMQ(IOFF),1)
         END IF
         IOFF=IOFF+NSIZE
      END DO

      CALL PUTREC(20,'JOBARC','COMPSYMQ',NSIZE*NSIZE*IINTFP,SYMQ)

      RETURN
      END

