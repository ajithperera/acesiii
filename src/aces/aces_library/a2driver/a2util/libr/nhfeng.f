      SUBROUTINE NHFENG(ICORE,MAXCOR,ETOT,NLIST2,NLIST1,NLIST1A,IUHF,
     &                  DOTAU)
C
C  DRIVER FOR THE CALCULATION OF THE CORRELATION ENERGY FOR A GIVEN SET
C  OF AMPLITUDES WHICH IS USED IN ROHF PERTURBATION THEORY CALCULATIONS
C 
C  ARGUMENTS :  ICORE ..... ICORE ARRAY
C               MAXCOR .... DIMENSION OF ICORE
C               NLIST2 .... OFFSET OF T2 LIST ON MOINTS (WITH RESPECT TO
C                            TYPE)
C               NLIST1 .... OFFSET OF T1 LISTS ON MOINTS (WITH RESPECT TO
C                            SPIN TYPE) FOR F-T1 PIECE
C               NLIST1 .... OFFSET OF T1 LISTS ON MOINTS (WITH RESPECT TO
C                            SPIN TYPE) FOR W-T1**2 PIECE
C               TECORR ..... RETURNS THE CORRELATION ENERGY FOR ALL SPIN CASES
C               IUHF .....  IUHF FLAG
C
CEND
      IMPLICIT INTEGER (A-Z)
      LOGICAL TAU,NONHF,DOTAU
      LOGICAL MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,SING1 
      CHARACTER*2 SPCASE(3)
      DOUBLE PRECISION E,ETOT,FACTOR,ECORR(3),ESPIN,ET2,ETOTT2,
     &                 ESING,SDOT
      DIMENSION ICORE(MAXCOR)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /SWITCH/ MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,SING1,
     &                QCISD
      COMMON /NHFREF/NONHF
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
C
      DIMENSION I0T(2),I0F(2)
C
      EQUIVALENCE (IFLAGS(2),METHOD)
C
      DATA SPCASE /'AA','BB','AB'/
C
      MXCOR=MAXCOR
C
C   ALLOCATE MEMORY FOR T1 AMPLITUDES
C
      I0T(1)=MXCOR+1-NT(1)*IINTFP
      MXCOR=MXCOR-NT(1)*IINTFP
      CALL GETLST(ICORE(I0T(1)),1,1,1,1+NLIST1,90)
      IF(IUHF.EQ.0) THEN
       I0T(2)=I0T(1)
      ELSE
       I0T(2)=I0T(1)-NT(2)*IINTFP
       MXCOR=MXCOR-NT(2)*IINTFP
       CALL GETLST(ICORE(I0T(2)),1,1,1,2+NLIST1,90)
      ENDIF
C
C   ALLOCATE MEMORY FOR f(a,I)
C
      I0F(1)=I0T(2)-NT(1)*IINTFP
      MXCOR=MXCOR-NT(1)*IINTFP
      CALL GETLST(ICORE(I0F(1)),1,1,1,3,93)
      IF(IUHF.EQ.0) THEN  
       I0F(2)=I0F(1)
      ELSE
       I0F(2)=I0F(1)-NT(2)*IINTFP
       MXCOR=MXCOR-NT(2)*IINTFP
       CALL GETLST(ICORE(I0F(2)),1,1,1,4,93)
      ENDIF
C
      ET1F=0.0
      ET2W=0.0
      ET12W=0.0
      ETOT=0.D0
      ETOTT2=0.D0
      FACTOR=1.D0
      IF(IUHF.EQ.0)FACTOR=2.D0
      DO 10 ISPIN=1,IUHF+1
       LISTT=43+ISPIN
       ESPIN=0.D0
C
C  COMPUTE THE F-T1 PIECE
C
       ESING=SDOT(NT(ISPIN),ICORE(I0T(ISPIN)),1,ICORE(I0F(ISPIN)),1)
       ESPIN=ESPIN+ESING
       ETOT=ETOT+FACTOR*ESING
       ETOTT2=ETOTT2+FACTOR*ESING
       ET1F=ET1F+ESING
C
C  COMPUTE THE TAU-W PIECE WITH AN APPROPRIATE TAU
C
       DO 100 IRREP=1,NIRREP
        DISSYT=IRPDPD(IRREP,ISYTYP(1,LISTT))
        NUMSYT=IRPDPD(IRREP,ISYTYP(2,LISTT))
        IF(MIN(NUMSYT,DISSYT).NE.0) THEN
        I001=1
        I002=I001+IINTFP*NUMSYT*DISSYT
        I003=I002+IINTFP*NUMSYT*DISSYT
        I004=I003+NUMSYT
        IF(DOTAU)THEN
         CALL GETLST(ICORE(I0T(ISPIN)),1,1,1,ISPIN+NLIST1A,90)
        ELSE
         CALL IZERO(ICORE(I0T(ISPIN)),NT(ISPIN)*IINTFP)
        ENDIF 
        IF(I004.LT.MXCOR) THEN
         CALL TENER(NLIST2,ET2,E,NUMSYT,DISSYT,ICORE(I001),
     &             ICORE(I002),ICORE(I0T(ISPIN)),ICORE(I0T(ISPIN)),
     &             ISPIN,DOTAU,IRREP,POP(1,ISPIN),POP(1,ISPIN),
     &             VRT(1,ISPIN),VRT(1,ISPIN),ICORE(I003))
        ELSE
         CALL INSMEM('CMPENG',I004,MXCOR)
        ENDIF
        ETOTT2=ETOTT2+FACTOR*ET2
        ETOT=ETOT+FACTOR*E
        ESPIN=ESPIN+E
        ENDIF
100    CONTINUE
       ECORR(ISPIN)=ESPIN
10    CONTINUE
C
C  THE TAU(Ij,Ab) <Ij//Ab> CONTRIBUTION TO THE ENERGY (SPIN CASE AB)
C
      ESPIN=0.D0
      DO 200 IRREP=1,NIRREP
       DISSYT=IRPDPD(IRREP,ISYTYP(1,46))
       NUMSYT=IRPDPD(IRREP,ISYTYP(2,46))
       IF(MIN(NUMSYT,DISSYT).NE.0) THEN
        I001=1
        I002=I001+IINTFP*NUMSYT*DISSYT
        I003=I002+IINTFP*NUMSYT*DISSYT
        I004=I003+NUMSYT
        IF(DOTAU)THEN
         CALL GETLST(ICORE(I0T(1)),1,1,1,1+NLIST1A,90)
         CALL GETLST(ICORE(I0T(2)),1,1,1,2+NLIST1A,90)
        ELSE
         CALL IZERO(ICORE(I0T(1)),NT(1))
         CALL IZERO(ICORE(I0T(2)),NT(2))
        ENDIF 
        IF(I004.LT.MXCOR) THEN
         CALL TENER(NLIST2,ET2,E,NUMSYT,DISSYT,ICORE(I001),
     &              ICORE(I002),ICORE(I0T(1)),ICORE(I0T(2)),3,DOTAU,
     &              IRREP,POP(1,1),POP(1,2),VRT(1,1),VRT(1,2),
     &              ICORE(I003))
        ELSE
         CALL INSMEM('CMPENG',I004,MXCOR)
        ENDIF
        ETOTT2=ETOTT2+ET2
        ETOT=ETOT+E
        ESPIN=ESPIN+E
       ENDIF
200   CONTINUE
      ECORR(3)=ESPIN
      WRITE(*,81)ETOT
81    FORMAT(T3,' The total correlation energy is ',F15.12,' a.u.')
      RETURN
      END