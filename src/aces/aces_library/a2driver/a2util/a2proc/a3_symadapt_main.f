













































































































































































































































































































































































































































































































































      subroutine A3_symadapt_main(Work, Icrsiz)
c
      implicit double precision (a-h,o-z)
c


c machsp.com : begin

c This data is used to measure byte-lengths and integer ratios of variables.

c iintln : the byte-length of a default integer
c ifltln : the byte-length of a double precision float
c iintfp : the number of integers in a double precision float
c ialone : the bitmask used to filter out the lowest fourth bits in an integer
c ibitwd : the number of bits in one-fourth of an integer

      integer         iintln, ifltln, iintfp, ialone, ibitwd
      common /machsp/ iintln, ifltln, iintfp, ialone, ibitwd
      save   /machsp/

c machsp.com : end










c This common block contains the IFLAGS and IFLAGS2 arrays for JODA ROUTINES
c ONLY! The reason is that it contains both arrays back-to-back. If the
c preprocessor define MONSTER_FLAGS is set, then the arrays are compressed
c into one large (currently) 600 element long array; otherwise, they are
c split into IFLAGS(100) and IFLAGS2(500).

c iflags(100)  ASVs reserved for Stanton, Gauss, and Co.
c              (Our code is already irrevocably split, why bother anymore?)
c iflags2(500) ASVs for everyone else

      integer        iflags(100), iflags2(500)
      common /flags/ iflags,      iflags2
      save   /flags/




C MXATMS     : Maximum number of atoms currently allowed
C MAXCNTVS   : Maximum number of connectivites per center
C MAXREDUNCO : Maximum number of redundant coordinates.
C
      INTEGER MXATMS, MAXCNTVS, MAXREDUNCO
      PARAMETER (MXATMS=200, MAXCNTVS = 10, MAXREDUNCO = 3*MXATMS)

c These parameters are gathered from vmol and vdint and are used by ecp
c as well. It just so happens that the vmol parameters do not exist in
c vdint and vice versa. LET'S TRY TO KEEP IT THAT WAY!

c VMOL PARAMETERS ------------------------------------------------------

C     MAXPRIM - Maximum number of primitives for a given shell.
      INTEGER    MAXPRIM
      PARAMETER (MAXPRIM=72)

C     MAXFNC  - Maximum number of contracted functions for a given shell.
C               (vmol/readin requires this to be the same as MAXPRIM)
      INTEGER    MAXFNC
      PARAMETER (MAXFNC=MAXPRIM)

C     NHT     - Maximum angular momentum
      INTEGER    NHT
      PARAMETER (NHT=7)

C     MAXATM  - Maximum number of atoms
      INTEGER    MAXATM
      PARAMETER (MAXATM=100)

C     MXTNPR  - Maximum total number of primitives for all symmetry
C               inequivalent centers.
      INTEGER    MXTNPR
      PARAMETER (MXTNPR=MAXPRIM*MAXPRIM)

C     MXTNCC  - Maximum total number of contraction coefficients for
C               all symmetry inequivalent centers.
      INTEGER    MXTNCC
      PARAMETER (MXTNCC=180000)

C     MXTNSH  - Maximum total number of shells for all symmetry
C               inequivalent centers.
      INTEGER    MXTNSH
      PARAMETER (MXTNSH=200)

C     MXCBF   - Maximum number of Cartesian basis functions for the
C               whole system (NOT the number of contracted functions).
c mxcbf.par : begin

c MXCBF := the maximum number of Cartesian basis functions (limited by vmol)

c This parameter is the same as MAXBASFN. Do NOT change this without changing
c maxbasfn.par as well.

      INTEGER MXCBF
      PARAMETER (MXCBF=1000)
c mxcbf.par : end

c VDINT PARAMETERS -----------------------------------------------------

C     MXPRIM - Maximum number of primitives for all symmetry
C              inequivalent centers.
      INTEGER    MXPRIM
      PARAMETER (MXPRIM=MXTNPR)

C     MXSHEL - Maximum number of shells for all symmetry inequivalent centers.
      INTEGER    MXSHEL
      PARAMETER (MXSHEL=MXTNSH)

C     MXCORB - Maximum number of contracted basis functions.
      INTEGER    MXCORB
      PARAMETER (MXCORB=MXCBF)

C     MXORBT - Length of the upper or lower triangle length of MXCORB.
      INTEGER    MXORBT
      PARAMETER (MXORBT=MXCORB*(MXCORB+1)/2)

C     MXAOVC - Maximum number of subshells per center.
      INTEGER    MXAOVC,    MXAOSQ
      PARAMETER (MXAOVC=32, MXAOSQ=MXAOVC*MXAOVC)

c     MXCONT - ???
      INTEGER    MXCONT
      PARAMETER (MXCONT=MXAOVC)

C
      character*32 szFile
      Character*2 Iroot(10)
      logical bExist, Spherical, Estate_calcs, Dens_diff 
      Character*4 Comp_pgrp, Full_pgrp
      Integer Roots
      Dimension Atommass(Mxatms), Iatmchrg(Mxatms),  
     &          Fucoord(3,Mxatms), Coord(3,Mxatms),
     &          Norbits_fullG(Mxatms), NOrbits_compG(Mxatms),
     &          Nbfns_4irrep(8)
C
      Dimension Work(Icrsiz/iintfp)
C
      Data Ione, Ieight, Iunit /1, 8, 10/
      Data Iroot /"01", "02", "03", "04", "05", "06", "07",
     &            "08", "09", "10"/
C    
      Iuhf = 1
      If (iflags(11).eq.0) iuhf = 0 
C
      Maxcor   = Icrsiz/Iintfp
      Length   = 0
c
c Read the JOBARC file for basic data of the molecule. 
c
      Call Getrec(20, 'JOBARC', 'NREALATM', Ione, Nreal_atoms)
      Call Getrec(20, 'JOBARC', 'NATOMS  ', Ione, Natoms)
      Call Getrec(20, 'JOBARC', 'COORD   ', 3*Natoms*Iintfp, Fucoord)
      Call Getrec(20, 'JOBARC', 'ATOMMASS', Natoms*Iintfp, Atommass)
      Call Getrec(20, 'JOBARC', 'ATOMCHRG', Natoms, Iatmchrg)
      Call Getrec(20, 'JOBARC', 'COMPNIRR', Ione, Nirrep)
      Call Getrec(20, 'JOBARC', 'NBASTOT ', Ione, Nbfns)
      Call Getrec(20, 'JOBARC', 'NAOBASFN', Ione, Naobfns)
      Call Getrec(20, 'JOBARC', 'NUMBASIR', Nirrep, Nbfns_4irrep)
C
      I = 1
      Do Iatm = 1, Natoms
         If (Iatmchrg(Iatm) .NE. 0) Then
            Coord(1, I) = Fucoord(1, Iatm)
            Coord(2, I) = Fucoord(2, Iatm)
            Coord(3, I) = Fucoord(3, Iatm)
            I = I + 1
         Endif
      Enddo
C
      Call Getcrec(20, 'JOBARC', 'COMPPTGP', 4, Comp_pgrp)
      Call Getcrec(20, 'JOBARC', 'FULLPTGP', 4, Full_pgrp)
      
      Spherical = .False. 
      If (iflags(62).eq.1) Spherical = .True.

      Write(6,*)
      Print*, "Variable at symadapt main"
      Write(*, '(a)'), "Nreal_atoms, Natoms, Nirrep, Nbfns,"
      Write(*, '(a,5i4)')"Naobfns: ", Nreal_atoms, 
     &                    Natoms, Nirrep, Nbfns,Naobfns
      Write(6,*) "Spherical :", Spherical 
      Write(6,"(a,2A4)")"Comp_pgrp, Full_pgrp: ", Comp_pgrp, Full_pgrp
      Write(6,*)
C
      IScfvec_a = 1
      IScfvec_b = IScfvec_a  + Naobfns*Naobfns
      ITmp1     = IScfvec_b  + Naobfns*Naobfns
      ITmp2     = Itmp1      + Naobfns*Naobfns
      ITmp2_a   = ITmp2      + Naobfns*Naobfns
      ITmp2_b   = ITmp2_a    + Naobfns*Naobfns
      IScfevl_a = ITmp2_b    + Naobfns*Naobfns
      IScfevl_b = IScfevl_a  + Naobfns
      Ioed2a_sc = IScfevl_b  + Naobfns
      Ioed2a_or = Ioed2a_sc  + Naobfns
      Inext     = Ioed2a_or  + Naobfns

      If (Inext .GE. Maxcor) Call Insmem("A3_symadapt_main", Inext,
     &                                    Maxcor)

      Imemleft = (Maxcor - Inext)

      Call a3_symadapt_scfvecs(Work(IScfvec_a), Work(IScfvec_b), 
     &                         Work(IScfevl_a), Work(IScfevl_b),
     &                         Work(Itmp1), Work(Itmp2),
     &                         Work(Ioed2a_sc), Work(Ioed2a_or),
     &                         Work(Itmp2_a), Work(Itmp2_b),
     &                         Nbfns,
     &                         Naobfns, Nbfns_4irrep, Nirrep, Iuhf,
     &                         Spherical, Work(Inext), Imemleft)
C
C Assign symmetry labels for excited states (excited state natural 
C orbitals. Bulk of the work necessary to do this was done by Tom
C Watson and I am simply moving things around to put all togther
C into a nice organized form. 
C
      Estate_calcs = (Iflags(87) .GT. 1) 

      if (Estate_calcs) Then
         call Getrec(20,'JOBARC','EESYMINF',1, Nroots)

         Iexsvec_a  = 1
         Igdens     = Iexsvec_a  + Naobfns*Naobfns
         Iedens     = Igdens     + Naobfns*Naobfns
         Inorbs     = Iedens     + Naobfns*Naobfns
         Iexseval_a = Inorbs     + Naobfns*Naobfns
         Itmp1      = Iexseval_a + Naobfns
         Iocc       = Itmp1      + Naobfns*Naobfns 
         Ioed2a_sc  = Iocc       + Naobfns
         Ioed2a_or  = Ioed2a_sc  + Naobfns
         Io_mos     = Ioed2a_or  + Naobfns
         It_mos     = Io_mos     + Naobfns
         IEng       = It_mos     + Naobfns
         Inext      = IEng       + Naobfns

         If (Inext .GE. Maxcor) Call Insmem("A3_symadapt_main", Inext,
     &                                     Maxcor)
      
         Imemleft = (Maxcor - Inext)

         Do Roots = 1, Nroots 
            Call a3_symadapt_estate_norbs(Work(Iexseval_a), Work(Iocc),
     &                                    Work(Iexsvec_a), 
     &                                    Work(Igdens), Work(Inorbs), 
     &                                    Work(Iedens), 
     &                                    Work(Ioed2a_sc), 
     &                                    Work(Ioed2a_or),
     &                                    Work(Itmp1),
     &                                    Work(Inext), 
     &                                    Naobfns, Nbfns, Imemleft,
     &                                    Iuhf, Iroot(Roots),
     &                                    Nbfns_4irrep, Nirrep,
     &                                    Work(Io_mos), 
     &                                    Work(It_mos),
     &                                    Work(Ieng))
         Enddo

      Endif

      Return
      End

