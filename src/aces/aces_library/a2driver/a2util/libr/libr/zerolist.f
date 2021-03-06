      SUBROUTINE ZEROLIST(ICORE,MAXCOR,DOLIST)
C
C THIS ROUTINE ZEROS OUT A DPD LIST
C
CEND
C SG 9/6/1996
C Modified to take care of cases when the list is not totally symmetric
C
      IMPLICIT INTEGER (A-Z)
      DIMENSION ICORE(MAXCOR)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
C
      MXNMDS = 0
      MXDSSZ = 0
      DO 10 IRREP = 1,NIRREP
        NUMDIS = IRPDPD(IRREP,ISYTYP(2,DOLIST))
        MXNMDS = MAX(NUMDIS,MXNMDS)
        DISSIZ = IRPDPD(IRREP,ISYTYP(1,DOLIST))
        MXDSSZ = MAX(DISSIZ,MXDSSZ)
 10   CONTINUE
      ISIZE = MXNMDS*MXDSSZ*IINTFP
      IF (ISIZE .LE. MAXCOR) THEN
        CALL IZERO(ICORE,ISIZE)
        DO 20 IRREP = 1,NIRREP
          NUMDIS = IRPDPD(IRREP,ISYTYP(2,DOLIST))
          CALL PUTLST(ICORE,1,NUMDIS,1,IRREP,DOLIST)
 20     CONTINUE
      ELSE
        CALL IZERO(ICORE,MXDSSZ*IINTFP)
        DO 30 IRREP = 1,NIRREP
          NUMDIS = IRPDPD(IRREP,ISYTYP(2,DOLIST))
          DO 40 IDIS = 1,NUMDIS
            CALL PUTLST(ICORE,IDIS,1,1,IRREP,DOLIST)
 40       CONTINUE
 30     CONTINUE
      ENDIF
C
      RETURN
      END
