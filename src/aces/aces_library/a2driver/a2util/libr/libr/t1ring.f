









































      SUBROUTINE T1RING(ICORE,MAXCOR,IUHF,LAMBDA)
C
C DRIVER FOR W(MBEJ) <- T1 CONTRIBUTIONS.  SELECTS BETWEEN INCORE
C AND OUT OF CORE ALGORITHMS AND CALLS THE APPROPRIATE ROUTINES
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ONEM,ZILCH
      LOGICAL LAMBDA,AOBASIS,AOLOG
      DIMENSION ICORE(MAXCOR)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFEA(2),NFMI(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /AOLOG/ AOLOG
      DATA ONE  /1.0/   
      DATA ONEM /-1.0/
      DATA ZILCH /0.0/
      integer disttype_vo(2)
      data disttype_vo /11,12/
      integer  aces_list_rows, aces_list_cols
      external aces_list_rows, aces_list_cols
C
      AOBASIS=IFLAGS(93).EQ.2.AND.(.NOT.LAMBDA)
C
C IF CALCULATION IN AOBASIS FIRST READ THE F(e;a)
C
      AOLOG=.FALSE.
      IF (AOBASIS) THEN
       AOLOG=.TRUE.
       I000=1
       DO 1111 ISPIN=1,1+IUHF
        CALL UPDMOI(1,NFMI(ISPIN),6+ISPIN,92,0,0)
        CALL GETLST(ICORE(I000),1,1,1,ISPIN,92)
        CALL PUTLST(ICORE(I000),1,1,1,6+ISPIN,92)
1111   CONTINUE
      ENDIF
C
C SPIN CASES AAAA AND BBBB
C
      INEED=-1
      DO 10 ISPIN=1,1+IUHF
       LSTTAR=53+ISPIN
       LSTINT1=6+ISPIN
       LSTINT2=26+ISPIN
       FULTYP1=18+ISPIN
       FULTYP2=20+ISPIN
       T1SIZE=NT(1)+IUHF*NT(2)
       ABFULL=IRPDPD(1,FULTYP1)
       MNFULL=IRPDPD(1,FULTYP2)
       DO 20 IRREP=1,NIRREP
        TARDSZ=aces_list_rows(IRREP,LSTTAR)
        TARDIS=aces_list_cols(IRREP,LSTTAR)
        TARSIZ=TARDSZ*TARDIS
        if (TARSIZ.LT.0) call trap_intovf('T1RING:TARSIZ',1)
        FULDSZ1=IRPDPD(IRREP,FULTYP1)
        FULDSZ2=IRPDPD(IRREP,FULTYP2)
        INTDSZ1=aces_list_rows(IRREP,LSTINT1)
        INTDIS1=aces_list_cols(IRREP,LSTINT1)
        INTDSZ2=aces_list_rows(IRREP,LSTINT2)
        INTDIS2=aces_list_cols(IRREP,LSTINT2)
        I1=TARDIS*TARSIZ+T1SIZE
        if (I1.LT.0) call trap_intovf('T1RING:I1',1)
        I2=3*TARDSZ
C
        I3A=FULDSZ1*INTDIS1+INTDSZ1
        I3B=FULDSZ2*INTDIS2+INTDSZ2
        if (I3A.LT.0) call trap_intovf('T1RING:I3A',1)
        if (I3B.LT.0) call trap_intovf('T1RING:I3B',1)
C
        I4A=2*ABFULL
        I4B=2*MNFULL
C
        INEEDA=I1+MAX(I2,I3A,I4A)
        INEEDB=I1+MAX(I2,I3B,I4B)
        INEED =MAX(INEED,INEEDA,INEEDB)
20     CONTINUE
10    CONTINUE
      IF(INEED.LT.MAXCOR/IINTFP)THEN
       CALL T1RA01(ICORE,MAXCOR,IUHF,LAMBDA)
      ELSE
       CALL T1RA02(ICORE,MAXCOR,IUHF,LAMBDA)
      ENDIF
C
C SPIN CASES ABBA AND BAAB
C
      INEED=-1
      DO 110 ISPIN=1,1+IUHF
       LSTTAR=57-ISPIN
       LSTINT1=31-ISPIN
       LSTINT2=8+ISPIN+(1-IUHF)
       T1SIZE=NT(1)+IUHF*NT(2)
       FULTYP1=18+ISPIN
       FULTYP2=20+ISPIN
       ABFULL=IRPDPD(1,FULTYP1)
       DO 120 IRREP=1,NIRREP
        TARDSZ=aces_list_rows(IRREP,LSTTAR)
        TARDIS=aces_list_cols(IRREP,LSTTAR)
        TARSIZ=TARDSZ*TARDIS
        if (TARSIZ.LT.0) call trap_intovf('T1RING:TARSIZ',2)
        INTDSZ1=aces_list_rows(IRREP,LSTINT1)
        INTDIS1=aces_list_cols(IRREP,LSTINT1)
        INTDSZ2=aces_list_rows(IRREP,LSTINT2)
        INTDIS2=aces_list_cols(IRREP,LSTINT2)
        I1=TARDIS*TARSIZ+T1SIZE
        if (I1.LT.0) call trap_intovf('T1RING:I1',2)
        I2=3*TARDSZ
C
        I3A=INTDSZ1*INTDIS1+3*MAX(INTDSZ1,INTDIS1)
        I3B=INTDSZ2*INTDIS2+3*MAX(INTDSZ2,INTDIS2)
        if (I3A.LT.0) call trap_intovf('T1RING:I3A',2)
        if (I3B.LT.0) call trap_intovf('T1RING:I3B',2)
C
        I4=2*ABFULL
C
        I5=TARDIS*TARDSZ+3*TARDSZ
        if (I5.LT.0) call trap_intovf('T1RING:I5',2)

        INEEDA=I1+MAX(I2,I3A,I4,I5)
        INEEDB=I1+MAX(I2,I3B,I4,I5)
        INEED =MAX(INEED,INEEDA,INEEDB)
120    CONTINUE
110   CONTINUE
      IF(INEED.LT.MAXCOR/IINTFP)THEN
       CALL T1RB01(ICORE,MAXCOR,IUHF,LAMBDA)
      ELSE
       CALL T1RB02(ICORE,MAXCOR,IUHF,LAMBDA)
      ENDIF
C
C STORE THE VALUES OF SUM(m,f) t(m;f) <ma||fe> IF AOBASIS IS TRUE
C DO THIS BY SUBTRACTING THE ACTUAL F(e;a)'s 
C BY THE F(e;a)'s READ IN BEFORE
C
      IF (AOBASIS) THEN
       I000=1
       DO 2222 ISPIN=1,1+IUHF
        I010=I000+NFMI(ISPIN)*IINTFP
        CALL GETLST(ICORE(I000),1,1,1,ISPIN,92)
        CALL GETLST(ICORE(I010),1,1,1,6+ISPIN,92)
        CALL SAXPY(NFMI(ISPIN),ONEM,ICORE(I010),1,ICORE(I000),1)
        CALL PUTLST(ICORE(I000),1,1,1,6+ISPIN,92)
 2222  CONTINUE
      ENDIF
C
C SPIN CASES ABAB AND BABA (FOR UHF ONLY)
C
c      IF(IUHF.EQ.0 )RETURN
      INEED=-1
      DO 210 ISPIN=1,1+IUHF
       LSTTAR=55+ISPIN
       LSTTMP1=37+ISPIN
       LSTTMP2=40-ISPIN
       LSTINT1=28+ISPIN
       LSTINT2=8+ISPIN
       T1SIZE=NT(1)+IUHF*NT(2)
       DO 220 IRREP=1,NIRREP
        TARDSZ=aces_list_rows(IRREP,LSTTAR)
        TARDIS=aces_list_cols(IRREP,LSTTAR)
        TARSIZ=TARDSZ*TARDIS
        if (TARSIZ.LT.0) call trap_intovf('T1RING:TARSIZ',3)
        TMPDSZ1=IRPDPD(IRREP,disttype_vo(3-ISPIN))
        TMPDIS1=IRPDPD(IRREP,disttype_vo(ISPIN))
        TMPDSZ2=IRPDPD(IRREP,disttype_vo(ISPIN))
        TMPDIS2=IRPDPD(IRREP,disttype_vo(3-ISPIN))
        INTDSZ1=IRPDPD(IRREP,13)
        INTDIS1=IRPDPD(IRREP,disttype_vo(3-ISPIN))
        INTDSZ2=IRPDPD(IRREP,14)
        INTDIS2=IRPDPD(IRREP,disttype_vo(ISPIN))
        I1=TARDIS*TARSIZ+T1SIZE
        if (I1.LT.0) call trap_intovf('T1RING:I1',3)
        I2=3*TARDSZ
C
        I3A=INTDSZ1*INTDIS1+3*MAX(INTDSZ1,INTDIS1)
        I3B=INTDSZ2*INTDIS2+3*MAX(INTDSZ2,INTDIS2)
        if (I3A.LT.0) call trap_intovf('T1RING:I3A',3)
        if (I3B.LT.0) call trap_intovf('T1RING:I3B',3)
C
        I4=TARDIS*TARDSZ+3*TARDSZ
        if (I4.LT.0) call trap_intovf('T1RING:I4',3)
C
        I5A=3*MAX(TMPDSZ1,TMPDIS1)
        I5B=3*MAX(TMPDSZ2,TMPDIS2)

        INEEDA=I1+MAX(I2,I3A,I4,I5A)
        INEEDB=I1+MAX(I2,I3B,I4,I5B)
        INEED =MAX(INEED,INEEDA,INEEDB)
220    CONTINUE
210   CONTINUE
      IF(INEED.LT.MAXCOR/IINTFP)THEN
       CALL T1RC01(ICORE,MAXCOR,IUHF,LAMBDA)
      ELSE
       CALL T1RC02(ICORE,MAXCOR,IUHF,LAMBDA)
      ENDIF
      RETURN
      END