      subroutine a3_symadapt_estate_norbs(ener,iocc,orb,dens,nlorb,
     &                                    edens,Oed2AScal, Ioed2Aord,
     &                                    tmp1,scr, nao,nbas,maxcor,
     &                                    iuhf,iroot,Nbfirr, Nirrep)
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)

#include "machsp.com"
#include "info.com"
#include "flags.com"
#include "sym.com"
C
#include "mxatms.par"
#include "maxangshell.par"
C
      parameter (one=1.0D0)
      parameter (zilch=0.0D0)
      parameter (DENS_THRESH=1.0D-08)
C
      double precision ener(nbas),orb(nao,nbas),scr(maxcor)
      double precision dens(nbas*nbas),edens(nbas,nbas),nlorb(nbas*nbas)
      double precision nelec
      double precision Oed2AScale(Nao), Ioed2Aorder(Nao), tmp1(Nao,Nao)
C
      Dimension Nprim_shell(Maxangshell*Mxatms)
      Dimension Orig_nprim_shell(Maxangshell*Mxatms)
      Integer   Reorder_Shell(Maxangshell*Mxatms)
C
      character*2 iroot
      character*6 denstype(2)
      character*6 string
      character*8 cscforb(2), Transform(4), Scfvecs(2)
      character*8 cexxcoef(2)
      character*8 cnumdrop(2)
      character*5 sptype(2)

      double precision  iocc(nbas)
      integer  idrppop(8),idrpvrt(8)
      integer  nocc(8,2), Nbfirr(8)
C
      data denstype /'REOMDN','LEOMDN'/
      data sptype   /'Alpha','Beta '/
      data Transform /"OOTRANSA", "VVTRNASA", "OOTRANSB", "VVTRANSB"/
      data Scfvecs /"SCFEVCA0", "SCFEVCB0"/
C
      call aces_com_info
      call aces_com_syminf
      call aces_com_sym
C
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
      Call Getrec(20, "JOBARC", "NSHELLS" , 1, nshells)
      Call Getrec(20, "JOBARC", "NPRMSHEL", nshells, Nprim_shell)
      Call Getrec(20, "JOBARC", "BNPAKORD", nshells, Reorder_Shell)
C
      Call Getrec(20, "JOBARC", "ERD2A2CS", Nbfns*Iintfp, Oed2AScale)
      Call Getrec(20, "JOBARC", "ERDORDER", Nbfns, Ioed2Aorder)
C
      i000=1
      i010=i000+nbas*nbas
      i020=i010+nbas*nbas
      i030=i020+nbas*nbas
      iend=i030+2*nbas*nbas

      if(Iend.gt.maxcor)
     &  call insmem('a3_symadapt_estate_norbs',iend,maxcor)

      do ispin=1,iuhf+1
         Do Idens = 1, 2
C
C Watson
C   Grab the AO density and transform to MO basis
C
CSSS        call getrec(20,'JOBARC','SCFVECA0',nbas*nbas*iintfp,scr(i000))
CSSS        call getrec(20,'JOBARC','AOOVRLAP',nbas*nbas*iintfp,scr(i010))
C#ifdef _DEBUG_LVL0
CSSS        write(6,*) "The SCF eigenvectors"
CSSS        call output(scr(i000), 1, nbas, 1, nbas, nbas, nbas, 1)
CSSS        write(6,*) "The overlap"
CSSS        call output(scr(i010), 1, nbas, 1, nbas, nbas, nbas, 1)
C#endif
CSSS
CSSS        CALL XGEMM('N','N',nbas,nbas,nbas,
CSSS     +                     1.0D0,SCR(I010),  nbas,  ! Overlap
CSSS     +                           SCR(I000),  nbas,  ! Coeff
CSSS     +                     0.0D0,SCR(I020),  nbas ) ! SC
CSSS
CSSS        CALL XGEMM('T','N',nbas,nbas,nbas,
CSSS     +                     1.0D0,SCR(I000),  nbas,  ! Ct
CSSS     +                           SCR(I020),  nbas,  ! SC
CSSS     +                     0.0D0,DENS     ,  nbas ) ! Ct SC
CSSS
C#ifdef _DEBUG_LVL0
CSSS        Write(6,"(a)") "CtSC = 1 check"
CSSS        CALL OUTPUT(DENS,1,NBAS,1,NBAS,NBAS,NBAS,1) ! = 1
C#endif
          String = Denstype(Idens)
          call getrec(20,'JOBARC',String//iroot,
     +                    nbas*nbas*iintfp,edens)

CSSS        call getrec(20,'JOBARC',"TDENSITY",nbas*nbas*iintfp,dens)
CSSS        CALL  DAXPY (NBAS*NBAS, -1.0D0, DENS, 1, EDENS, 1)
CSSS        CALL  DSCAL (NBAS*NBAS,0.5D0,EDENs,1)
CSSS        CALL  DCOPY (NBAS*NBAS,EDENS,1,SCR(I000),1)

#ifdef _DEBUG_LVL0
        write(6,*) 
        Write(6,"(a,a)") "Excited state density: ", String//iroot
        CALL OUTPUT(EDENS,1,NBAS,1,NBAS,NBAS,NBAS,1) 
#endif

CSSS        CALL XGEMM('N','N',nbas,nbas,nbas,
CSSS     +                     1.0D0,SCR(I000),  nbas,
CSSS     +                           SCR(I010),  nbas,
CSSS     +                     0.0D0,SCR(I020),  nbas )
CSSS
CSSS        CALL XGEMM('N','N',nbas,nbas,nbas,
CSSS     +                     1.0D0,SCR(I010),  nbas,
CSSS     +                           SCR(I020),  nbas,
CSSS     +                     0.0D0,SCR(I000),  nbas )
CSSS        CALL  DCOPY (NBAS*NBAS,SCR(I000),1,DENS,1)
CSSS        call ao2mo2(scr(i000),dens,nlorb,scr(i030),nbas,nbas,ispin)
C#ifdef _DEBUG_LVL0
CSSS        WRITE (6,"(a)") ' CSDSC before eig! '
CSSS        CALL OUTPUT(DENS,1,NBAS,1,NBAS,NBAS,NBAS,1) ! = 1
C#endif

C
C Generate the occupation numbers for each irrep based on eigen
C values (SCF) and the number of basis functions per irrep.
C
        Call Occupy(Nirrep, Nbfirr, Nbas, ener, scr, Nocc(1,Ispin),
     &              1)
        If (Iuhf .EQ. 0) Call Icopy(8, Nocc(1, 1), 1, Nocc(1,2), 1)

#ifdef _DEBUG_LVL0
        Write(6,"(a)") "The number of occupied orbital per irrep"
        Write(6,"(6(1x,i3))") (Nocc(i,1),i=1,Nirrep)
        If (Iuhf .EQ.1) Write(6,"(6(1x,i3))") (Nocc(i,2),i=1,Nirrep)
#endif

      Noccs = 0
      Do irrep = 1, Nirrep
         Noccs = Noccs + Nocc(irrep,Ispin)
      Enddo
      Nvrts = Nbas - Noccs
C 
C Set the Norb depending on whether we do occ-occ or vrt-vrt block of
C the transition density. The current setup is to handle occ-occ block
C first (idens=1).
C
      If (Idens .Eq. 1) Then
         Norbs = Noccs
         Ioff  = 0
      Else
         Norbs = Nvrts
         Ioff  = Noccs
      Endif
   
      Do Jc = 1, Norbs
         Do Ir = 1, Norbs
            indx = ir + (jc-1)*Norbs
            Dens(Indx) = Edens(Ioff+ir, Ioff+jc)
         Enddo
      Enddo

#ifdef _DEBUG_LVL0
        If (idens .Eq.1) Then
           WRITE (6,"(a)") ' The OO block of trans. density matrices'
        Else
           WRITE (6,"(a)") ' The VV block of trans. density matrices'
        Endif
        CALL OUTPUT(dens,1,Norbs,1,Norbs,Norbs,Norbs,1) 
#endif
C
C     Symmetrize and Diagonalize MO basis density
C
        call symmet2(Dens, Norbs)
        call eig (dens,Nlorb,1,Norbs,-1)

#ifdef _DEBUG_LVL0
        Write(6,*)
        WRITE (6,"(a)") 'The eigenvalues and natural orbitals'
        CALL OUTPUT(dens,1,Norbs,1,Norbs,Norbs,NOrbs,1) 
        CALL OUTPUT(NlOrb,1,Norbs,1,NOrbs,Norbs,Norbs,1) 
#endif
        call dcopy(Norbs,dens,1,iocc,1)
        nelec = 0.0D0
        DO imo = 1, Norbs
           nelec = nelec + iocc(imo)
        ENDDO

#ifdef _DEBUG_LVL0
        Write(6,*)
        WRITE (6,"(a)") ' Trace of diagonalized MO density matrix: '
        WRITE (6,"(F10.5)")  nelec
        Write(6,"(a)") "Eigenvalues"
        Write(6,"(6(1x,F10.5))") (iocc(imo), imo=1, Norbs)

#endif 
C
C Write the OO and VV transformation matrices to JOBARC 
C
C        Call putrec(20, "JOBARC", TRANSFORM(Idens + (ISpin-1)), 
C       &            Norbs*Norbs*Iintfp, Nlorb)
C
C Do the OO and VV transformation. Retrive the SCF vectors of correct
C spin type during first iteration of the inner loop (OO block)
C
        call getrec(20,'JOBARC', SCFVECS(Ispin), Nbas*Nbas*iintfp,dens)

        CALL XGEMM('N','N',nbas,norbs,norbs,
     +                     1.0D0,DENS,   nbas,
     +                           NLORB,  norbs,
     +                     0.0D0,SCR(I010),    nbas )

CSSS       If (Idens .EQ. 1) Call Dcopy (Nbas*Norbs, SCR(I010), 1, Dens, 1)
       
#ifdef _DEBUG_LVL0
        If (Idens .EQ. 1) Then
           WRITE (6,"(a)") ' The OO rotated MO vectors'
           CALL OUTPUT(SCR(I010),1,NBAS,1,Norbs,NBAS,Norbs,1) 
        Else
           WRITE (6,"(a)") ' The VV rotated MO vectors'
           CALL OUTPUT(SCR(I010),1,NBAS,1,Norbs,NBAS,Norbs,1) 
        Endif
#endif
        Enddo

C Before we proceed further we need to convert the vectors from 
C OED (ACESIII) scaling and order to ACESII scaling and order. 
C
C       Call Do_oed_to_vmol(Norbs, Ioed2Aorder, Oed2AScale, Tmp1, 
C    &                      Scr(I010))
C      If (Iuhf .EQ. 1) then
C         Call Do_oed_to_vmol(Norbs, Ioed2Aorder, Oed2AScale, Tmp1,
C     &                      Scr(I010))
C      Endif

#ifdef _NOSKIP

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
        Call Get_irreps(orb, Ener, Scr, Imemleft*Iintfp, Nbas,
     &                  Nao, Ispin, Nocc, Iuhf)
#endif
      Enddo

      return
      end




