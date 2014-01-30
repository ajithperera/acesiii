         SUBROUTINE  OED__GENER_delta_BATCH
     +
     +                    ( IMAX,ZMAX,
     +                      NALPHA,NCOEFF,NCSUM,
     +                      NCGTO1,NCGTO2,
     +                      NPGTO1,NPGTO2,
     +                      SHELL1,SHELL2,
     +                      X1,Y1,Z1,X2,Y2,Z2,
     +                      Xn,Yn,Zn,ALPHA,
     *                      CC,CCBEG,CCEND,
     +                      SPHERIC,
     +                      ICORE,
     +                                NBATCH,
     +                                NFIRST,
     +                                ZCORE )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__GENER_delta_BATCH
C  MODULE      : delta INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : OED__delta_CSGTO
C  DESCRIPTION : Main operation that drives the calculation of a batch
C                of contracted electron overlap integrals.
C
C
C                  Input (x = 1 and 2):
C
C                    IMAX,ZMAX    =  maximum integer,flp memory
C                    NALPHA       =  total # of exponents
C                    NCOEFF       =  total # of contraction coeffs
C                    NCSUM        =  total # of contractions
C                    NCGTOx       =  # of contractions for csh x
C                    NPGTOx       =  # of primitives per contraction
C                                    for csh x
C                    SHELLx       =  the shell type for csh x
C                    Xy,Yy,Zy     =  the x,y,z-coordinates for centers
C                                    y = 1 and 2
C                    XN,YN,ZN     =  the x,y,z-coordinates for all
C                                    nuclear attraction centers
C                    ALPHA        =  primitive exponents for csh
C                                    1 and 2 in that order
C                    CC           =  contraction coefficient for csh
C                                    1 and 2 in that order, for each
C                                    csh individually such that an
C                                    (I,J) element corresponds to the
C                                    I-th primitive and J-th contraction
C                    CC(BEG)END   =  (lowest)highest nonzero primitive
C                                    index for contractions for csh
C                                    1 and 2 in that order. They are
C                                    different from (1)NPGTOx only for
C                                    segmented contractions
C                    SPHERIC      =  is true, if spherical integrals
C                                    are wanted, false if cartesian
C                                    ones are wanted
C                    ICORE        =  integer scratch space
C                    ZCORE (part) =  flp scratch space
C
C
C                  Output:
C
C                    NBATCH       =  # of integrals in batch
C                    NFIRST       =  first address location inside the
C                                    ZCORE array containing the first
C                                    integral
C                    ZCORE        =  full batch of contracted (1|2)
C                                    overlap integrals over cartesian
C                                    or spherical gaussians starting
C                                    at ZCORE (NFIRST)
C
C
C
C                            !!! IMPORTANT !!!
C
C                For performance tuning, please see the include file
C                'oed__tuning.inc'.
C
C
C  Adopted from OED__GENER_OVL_BATC  to obtain delta integrals : Prakash Verma
c  
C------------------------------------------------------------------------
C
C             ...include files and declare variables.
C
C
         IMPLICIT    NONE

         INCLUDE     'oed__tuning.inc'

         LOGICAL     SPHERIC

         INTEGER     IMAX,ZMAX
         INTEGER     NALPHA,NCOEFF,NCSUM
         INTEGER     NBATCH,NFIRST
         INTEGER     NCGTO1,NCGTO2
         INTEGER     NPGTO1,NPGTO2
         INTEGER     SHELL1,SHELL2

         INTEGER     CCBEG (1:NCSUM)
         INTEGER     CCEND (1:NCSUM)
         INTEGER     ICORE (1:IMAX)

         DOUBLE PRECISION  X1,Y1,Z1,X2,Y2,Z2

         DOUBLE PRECISION  ALPHA (1:NALPHA)
         DOUBLE PRECISION  CC    (1:NCOEFF)
         DOUBLE PRECISION  ZCORE (1:ZMAX)
         
         integer lmnvala(3,36), lmnvalb(3,36)
         double precision  Xn, Yn, Zn
C
C
C------------------------------------------------------------------------
C
C
C             ...call csgto routine.
C
C

         CALL  OED__delta_CSGTO
     +
     +              ( IMAX,ZMAX,
     +                NALPHA,NCOEFF,NCSUM,
     +                NCGTO1,NCGTO2,
     +                NPGTO1,NPGTO2,
     +                SHELL1,SHELL2,
     +                X1,Y1,Z1,X2,Y2,Z2,
     +                Xn,Yn,Zn,
     +                ALPHA,
     +                CC,CCBEG,CCEND,
     +                L1CACHE,TILE,NCTROW,
     +                SPHERIC,
     +                ICORE,
     +                          NBATCH,
     +                          NFIRST,
     +                          ZCORE )
     +
     +
C
C
C             ...ready!
C
C
         RETURN
         END
