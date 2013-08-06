      SUBROUTINE SYMOFF(IOFFSM,IRREPDO)
C
C COMPUTE OFFSETS INTO AN i,A DISTRIBUTION FOR THIS IRREP (RIGHT INDEX - 1)
C                         I,a DISTRIBUTION FOR THIS IRREP (RIGHT INDEX - 2)
C                         b,M DISTRIBUTION FOR THIS IRREP (RIGHT INDEX - 3)
C                         B,m DISTRIBUTION FOR THIS IRREP (RIGHT INDEX - 4)
C
CEND
      IMPLICIT INTEGER (A-Z)
      DIMENSION IOFFSM(8,4)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFEA(2),NFMI(2)
      IOFF1=0
      IOFF2=0
      IOFF3=0
      IOFF4=0
      DO 1001 IRREPA=1,NIRREP
       IRREPM=IRREPA
       IRREPI=DIRPRD(IRREPA,IRREPDO)
       IRREPB=IRREPI
       IOFFSM(IRREPA,1)=IOFF1
       IOFFSM(IRREPA,2)=IOFF2
       IOFFSM(IRREPM,3)=IOFF3
       IOFFSM(IRREPM,4)=IOFF4
       IOFF1=IOFF1+POP(IRREPI,2)*VRT(IRREPA,1)
       IOFF2=IOFF2+POP(IRREPI,1)  *VRT(IRREPA,2)
       IOFF3=IOFF3+VRT(IRREPB,2)*POP(IRREPM,1)
       IOFF4=IOFF4+VRT(IRREPB,1)  *POP(IRREPM,2)
1001  CONTINUE
      RETURN
      END