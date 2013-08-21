      subroutine a3_symadapt_estate_norbs(ener,iocc,orb,dens,nlorb,
     &                                    edens,scr, nao,nbas,maxcor,
     &                                    iuhf,iroot,Nbfirr, Nirrep)
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      double precision ener(nbas),orb(nao,nbas),scr(maxcor)
      double precision dens(nbas,nbas),edens(nbas,nbas),nlorb(nbas,nbas)
      double precision nelec
      character*2 iroot
      character*8 cscfener(2)
      character*8 cscforb(2)
      character*8 cexxcoef(2)
      character*8 cnumdrop(2)
      character*5 sptype(2)
      double precision  iocc(nbas)
      integer  idrppop(8),idrpvrt(8)
      integer  nocc(8,2), Nbfirr(8)

#include "machsp.com"
#include "info.com"
#include "flags.com"
#include "sym.com"
C
      parameter (one=1.0D0)
      parameter (zilch=0.0D0)
      parameter (DENS_THRESH=1.0D-08)
      data cscfener /'SCFEVLA0','SCFEVLB0'/
      data sptype   /'Alpha','Beta '/
      call aces_com_info
      call aces_com_syminf
      call aces_com_sym
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c the mo ordering is important for getting the occupations right
c
c depending on the place in the program (if the user is not using xaces2)
c it may be either scf or correlated
c
c    so            mo                        mo
c ao CMP2ZMAT   so cscforb(ispin)    =    ao AOBASMOS
c
c the so ordering is zmat
c the mo ordering is correlated if calc is greater than scf 
c or if this is a vibrational calculation
c
      i000=1
      i010=i000+nbas*nbas
      i020=i010+nbas*nbas
      i030=i020+nbas*nbas
      iend=i030+2*nbas*nbas

      if(Iend.gt.maxcor)
     &  call insmem('NLORB-F',iend,maxcor)

      do 10 ispin=1,iuhf+1
C
C Watson
C   Grab the AO density and transform to MO basis
C
        call getrec(20,'JOBARC','SCFVECA0',nbas*nbas*iintfp,scr(i000))
        call getrec(20,'JOBARC','AOOVRLAP',nbas*nbas*iintfp,scr(i010))
#ifdef _DEBUG_LVL0
        write(6,*) "The SCF eigenvectors"
        call output(scr(i000), 1, nbas, 1, nbas, nbas, nbas, 1)
        write(6,*) "The overlap"
        call output(scr(i010), 1, nbas, 1, nbas, nbas, nbas, 1)
#endif

        CALL XGEMM('N','N',nbas,nbas,nbas,
     +                     1.0D0,SCR(I010),  nbas,  ! Overlap
     +                           SCR(I000),  nbas,  ! Coeff
     +                     0.0D0,SCR(I020),  nbas ) ! SC

        CALL XGEMM('T','N',nbas,nbas,nbas,
     +                     1.0D0,SCR(I000),  nbas,  ! Ct
     +                           SCR(I020),  nbas,  ! SC
     +                     0.0D0,DENS     ,  nbas ) ! Ct SC

#ifdef _DEBUG_LVL0
        Write(6,"(a)") "CtSC = 1 check"
        CALL OUTPUT(DENS,1,NBAS,1,NBAS,NBAS,NBAS,1) ! = 1
#endif
          call getrec(20,'JOBARC','EXCDEN'//iroot,
     +                    nbas*nbas*iintfp,edens)
        call getrec(20,'JOBARC',"TDENSITY",nbas*nbas*iintfp,dens)
 
        CALL  DAXPY (NBAS*NBAS, -1.0D0, DENS, 1, EDENS, 1)
        CALL  DSCAL (NBAS*NBAS,0.5D0,EDENs,1)
        CALL  DCOPY (NBAS*NBAS,EDENS,1,SCR(I000),1)

#ifdef _DEBUG_LVL0
        write(6,*) 
        Write(6,"(a)") "Excited state density"
        CALL OUTPUT(scr(I000),1,NBAS,1,NBAS,NBAS,NBAS,1) ! = 1
#endif
        CALL XGEMM('N','N',nbas,nbas,nbas,
     +                     1.0D0,SCR(I000),  nbas,
     +                           SCR(I010),  nbas,
     +                     0.0D0,SCR(I020),  nbas )

        CALL XGEMM('N','N',nbas,nbas,nbas,
     +                     1.0D0,SCR(I010),  nbas,
     +                           SCR(I020),  nbas,
     +                     0.0D0,SCR(I000),  nbas )

        CALL  DCOPY (NBAS*NBAS,SCR(I000),1,DENS,1)

        call ao2mo2(scr(i000),dens,nlorb,scr(i030),nbas,nbas,ispin)

        DO imo = 1,NBAS
        DO jmo = 1,NBAS
           DVAL = DABS (DENS (imo,jmo))
           IF (DVAL .LT. DENS_THRESH) DENS (imo,jmo) = 0.0D0
        ENDDO
        ENDDO

#ifdef _DEBUG_LVL0
        WRITE (6,"(a)") ' CSDSC before eig! '
        CALL OUTPUT(DENS,1,NBAS,1,NBAS,NBAS,NBAS,1) ! = 1
#endif
C
C     Diagonalize MO basis density
C
        call zero(nlorb,NBAS*NBAS)
        call eig (dens,nlorb,1,nbas,-1)

#ifdef _DEBUG_LVL0
        WRITE (6,"(a)") ' Nat. Orbs. MO x MO ! '
        CALL OUTPUT(nlorb,1,NBAS,1,NBAS,NBAS,NBAS,1) ! = 1
#endif
        call dcopy(nbas,dens,nbas+1,iocc,1)
        nelec = 0.0D0
        DO imo = 1, nbas
           nelec = nelec + iocc(imo)
        ENDDO

#ifdef _DEBUG_LVL0
        WRITE (6,"(a)") ' Trace of diagonalized MO density '
        WRITE (6,"(a,1x,F10.5)")  'matrix - ',nelec
#endif 

        call getrec(20,'JOBARC','SCFVECA0',nbas*nbas*iintfp,dens)
        CALL XGEMM('N','N',nbas,nbas,nbas,
     +                     1.0D0,DENS,   nbas,
     +                           NLORB,  nbas,
     +                     0.0D0,SCR(I010),    nbas )

#ifdef _DEBUG_LVL0
        WRITE (6,"(a)") ' Nat. Orbs. AO x MO ! '
        CALL OUTPUT(SCR(I010),1,NBAS,1,NBAS,NBAS,NBAS,1) ! = 1
#endif

        call getrec(20,'JOBARC','CMP2ZMAT',nao*nbas*iintfp,scr(i000))

#ifdef _DEBUG_LVLM
        WRITE (6,"(a)") ' CMP2ZMAT ! '
        CALL OUTPUT(scr(i000),1,NAO,1,NBAS,NAO,NBAS,1) ! = 1
#endif
        call xgemm('n','n',nao,nbas,nbas,one,scr(i000),nao,scr(i010),
     &    nbas,zilch,orb,nao)

#ifdef _DEBUG_LVL0
        WRITE (6,"(a)") ' Nat. Orbs. NAOBFNS x MO ! '
        CALL OUTPUT(orb,1,NAO,1,NBAS,NAO,NBAS,1) ! = 1
#endif 
        call getrec(-1,'JOBARC',cscfener,nbas*iintfp,ener)
C
C Generate the occupation numbers for each irrep based on eigen
C values and the number of basis functions per irrep.
C
        Call Occupy(Nirrep, Nbfirr, Nbas, ener, scr, Nocc(1,1),
     &              1)
        If (Iuhf .EQ. 0) Call Icopy(8, Nocc(1, 1), 1, Nocc(1,2), 1)
C
        Call Get_irreps(orb, Ener, Scr, Imemleft*Iintfp, Nbas,
     &                  Nao, Ispin, Nocc, Iuhf)
10    continue

      return
      end




