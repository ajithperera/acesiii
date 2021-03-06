      SUBROUTINE MO2AO2(ZMO,ZAO,EVEC,SCR,NBAS,ISPIN)
C
C TRANSFORMS A TWO-INDEX QUANTITY (I.E. DENSITY) FROM THE MO TO
C  THE AO BASIS.
C
C INPUT:
C       ZMO  - TWO-INDEX QUANTITY IN MO BASIS (LENGTH: NBAS*NBAS)
C       NBAS - NUMBER OF BASIS FUNCTIONS
C       ISPIN- 1 FOR ALPHA, 2 FOR BETA
C
C OUTPUT:
C       ZAO  - TWO-INDEX QUANTITY IN AO BASIS (LENGTH: NBAS*NBAS)
C       
C SCRATCH:
C       EVEC - HOLDS EIGENVECTOR MATRIX (LENGTH: NBAS*NBAS)
C       SCR  - HOLDS INTERMEDIATE QUANTITY (LENGTH: NBAS*NBAS) 
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*8 LABEL
      DIMENSION EVEC(NBAS,NBAS),ZMO(NBAS,NBAS),ZAO(NBAS,NBAS)
      DIMENSION SCR(NBAS,NBAS)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      IF(ISPIN.EQ.1)LABEL='SCFEVCA0'
      IF(ISPIN.EQ.2)LABEL='SCFEVCB0'
      CALL GETREC(20,'JOBARC',LABEL,IINTFP*NBAS*NBAS,EVEC)
      CALL MTRAN2(EVEC,NBAS)
      CALL MXM(ZMO,NBAS,EVEC,NBAS,SCR,NBAS)
      CALL MTRAN2(EVEC,NBAS)
      CALL MXM(EVEC,NBAS,SCR,NBAS,ZAO,NBAS)
      RETURN
      END
