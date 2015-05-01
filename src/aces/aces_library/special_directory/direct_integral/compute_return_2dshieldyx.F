      subroutine compute_return_2dshieldyx(a1,a2,b1,b2,c1,c2,d1,d2,scr,
     *                 maxblk, iscr, coords,coeffs,alphas, ccbeg, ccend,
     *                 DSHIELD_YX,nsend)

c---------------------------------------------------------------------------
c   The block of integrals (a1:a2,b1:b2,c1:c2,d1:d2) is computed for the
c   following 'types' of integrals based on atomic labels.
c---------------------------------------------------------------------------

      implicit none

      include 'mpif.h'
      include 'int_gen_parms.h'
      include 'machine_types.h'

      integer max_dim_coeff
      parameter (max_dim_coeff = 5000)
      integer a1, a2, b1, b2, c1, c2, d1, d2
      integer aa1,aa2,bb1,bb2,cc1,cc2,dd1,dd2
      integer adim, bdim, cdim, ddim
      integer m1, m2, n1, n2, r1, r2, s1, s2
      integer i, j, n, m, r, s
      integer matm, natm, ratm, satm
      integer a,b,c,d

      integer num_to_do, nsend
      integer nints, maxblk
      integer nalpha_pack, npcoeff_pack
      integer ncsum, next, nfirst
      integer me, ierr
      integer nc1, nc2, nd1, nd2
      integer imin, zmin, iblk, zblk
      Integer XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ
      integer jatom, jx, jcenter
      integer*8 arg64(25)
      integer iscr(*)
      integer ccbeg(*), ccend(*)
      integer ccbeg_pack(max_dim_coeff), ccend_pack(max_dim_coeff)
      integer*8 ccbeg_pack64(max_dim_coeff), ccend_pack64(max_dim_coeff)

      logical skip
      logical mn_symmetry
      logical rs_symmetry
      logical mn_rs_symmetry
      logical*8 l8true, l8spherical
      logical spherical

      double precision x1,y1,z1
      double precision x2,y2,z2
      double precision x3,y3,z3
      double precision x4,y4,z4
      double precision xx_s,yx_s,yy_s,zx_s,zy_s,zz_s
      double precision xy_s,xz_s,yz_s
      double precision xx_l,yx_l,yy_l,zx_l,zy_l,zz_l
      double precision xy_l,xz_l,yz_l
      double precision coords(3,*), coeffs(*), alphas(*)
      double precision DSHIELD_YX(a1:a2,b1:b2,c1:c2,d1:d2)
      double precision scr(*)
      double precision sz_max, delta
      double precision itol, bmax, dtemp, emax
      double precision alpha_pack(max_dim_coeff),
     *                 pcoeff_pack(max_dim_coeff)

      common /Imax_com/sz_max(max_shells,max_shells), delta
      common /d2int_com/jatom, jx, jcenter

      save me,alpha_pack, pcoeff_pack, ccbeg_pack, ccend_pack,
     *     ccbeg_pack64, ccend_pack64

      call mpi_comm_rank(mpi_comm_world, me, ierr)
c      print *,'Task ',me,' computing integrals for ',a1,a2,b1,b2,
c     *     c1,c2,d1,d2
c      call c_flush_stdout()

      adim = a2-a1+1
      bdim = b2-b1+1
      cdim = c2-c1+1
      ddim = d2-d1+1
      l8true = .true.
      spherical = (ispherical .eq. 1)
      l8spherical = spherical

c Set the integral tolerance

      call set_itol(delta,itol)

      nsend = adim*bdim*cdim*ddim
      if (nsend .lt. 0) then
         print *,'ERROR IN INTEGRAL WORKER ',me,' nsend = ',nsend
         print *,'adim,bdim,cdim,ddim = ',adim,bdim,cdim,ddim
         call mpi_abort(mpi_comm_world, ierr)
      endif

      do d = d1,d2
      do c = c1,c2
      do b = b1,b2
      do a = a1,a2
         DSHIELD_YX(a,b,c,d) = 0.d0
      enddo
      enddo
      enddo
      enddo
c-----------------------------------------------------------------------
c   Find the shell blocks for which we shall loop through.
c-----------------------------------------------------------------------

      call lookup_shell(end_nfps, nshells, a1, m1)
      call lookup_shell(end_nfps, nshells, a2, m2)
      call lookup_shell(end_nfps, nshells, b1, n1)
      call lookup_shell(end_nfps, nshells, b2, n2)
      call lookup_shell(end_nfps, nshells, c1, r1)
      call lookup_shell(end_nfps, nshells, c2, r2)
      call lookup_shell(end_nfps, nshells, d1, s1)
      call lookup_shell(end_nfps, nshells, d2, s2)

CSSS         Write(6,*)
CSSS         Write(6,*) "Loop over shell quads."
CSSS         Write(6,"(4(a,1x,I4,1x,I4,1x))") "m1 and m2: ", m1, m2,
CSSS     &                                    "n1 and N2: ", n1, n2,
CSSS     &                                    "r1 and r2: ", r1, r2,
CSSS     &                                    "s1 and s2: ", s1, s2

      do m = m1, m2
         if (m .eq. 1) then
             aa1 = 1
         else
             aa1 = end_nfps(m-1) + 1
         endif
         aa2 = end_nfps(m)

         x1 = coords(1,m)
         y1 = coords(2,m)
         z1 = coords(3,m)

         do n = n1, n2
            if (n .eq. 1) then
               bb1 = 1
            else
               bb1 = end_nfps(n-1) + 1
            endif
            bb2 = end_nfps(n)

            x2 = coords(1,n)
            y2 = coords(2,n)
            z2 = coords(3,n)

            do r = r1, r2
              if (r .eq. 1) then
                  cc1 = 1
              else
                  cc1 = end_nfps(r-1) + 1
              endif
              cc2 = end_nfps(r)

              x3 = coords(1,r)
              y3 = coords(2,r)
              z3 = coords(3,r)

              do s = s1, s2
                 x4 = coords(1,s)
                 y4 = coords(2,s)
                 z4 = coords(3,s)

                 call pack_coeffs(alphas, ixalpha, coeffs, ixpcoef,
     *                          ncfps, npfps, m, n,
     *                          r, s, alpha_pack, nalpha_pack,
     *                          pcoeff_pack, npcoeff_pack,
     *                          ccbeg, ccend, indx_cc,
     *                          ccbeg_pack, ccend_pack)

c---------------------------------------------------------------------------
c   Calling sequence for ERD version 2.
c---------------------------------------------------------------------------
       ncsum = ncfps(m) + ncfps(n) + ncfps(r) + ncfps(s)
C
C      Write(6,*) "@-entry to DSHIELD_BATCH"
C      Write(6,"((a,4(1x,i4)))") "m, n, r, s:", m, n, r, s
C      Write(6,"(a,2(1x,I4))") "The Nalpha_pack and Npcoeff_pack:",
C     &                         nalpha_pack, npcoeff_pack
C      Write(6,"((a,4(1x,i4)))") "Ivanagmom:", Ivangmom(m),
C     &      Ivangmom(n), Ivangmom(r), Ivangmom(s)
C      Write(6,"((a,4(1x,i4)))") "Ncfps:", Ncfps(m),
C     &      Ncfps(n), Ncfps(r), Ncfps(s)
C      Write(6,"((a,4(1x,i4)))") "Npfps:", Npfps(m),
C     &      Npfps(n), Npfps(r), Npfps(s)
C      Write(6,*)
C      Write(6,"(4(1x,F15.7))") (Alpha_pack(i),i=1,Nalpha_pack)
C      Write(6,*)
C      Write(6,"(4(1x,F15.7))")(Pcoeff_pack(i),i=1,npcoeff_pack)
C      Write(6,*)
C      Write(6,"(4(1x,I4))")(CCbeg_pack(i),i=1,Nalpha_pack)
C      Write(6,*)
C      Write(6,"(4(1x,I4))")(CCend_pack(i),i=1,Nalpha_pack)
C      Write(6,*)
C      Write(6,"(a,3(1x,I10))") "The intmax,Zmax,ncsum:",intmax, zmax,
C     &                       ncsum
C      Write(6,*) "The coordiantes"
C      Write(6,*)
C      Write(6,"(3(1x,F10.5))") x1, y1, z1
C      Write(6,"(3(1x,F10.5))") x2, y2, z2
C      Write(6,"(3(1x,F10.5))") x3, y3, z3
c      Write(6,"(3(1x,F10.5))") x4, y4, z4

           XX=0
           XY=0
           XZ=0
           YX=1
           YY=0
           YZ=0
           ZX=0
           ZY=0
           ZZ=0
           call ERD__GENER_DSHIELD_BATCH(intmax, zmax,
     *                nalpha_pack, npcoeff_pack, ncsum,
     *                ncfps(m),ncfps(n), ncfps(r), ncfps(s),
     *                npfps(m),npfps(n), npfps(r), npfps(s),
     *                ivangmom(m), ivangmom(n),
     *                ivangmom(r), ivangmom(s), x1,y1,z1,
     *                x2,y2,z2,x3,y3,z3,x4,y4,z4,
     *                xx,xy,xz,yx,yy,yz,zx,zy,zz,
     *                alpha_pack,
     *                pcoeff_pack, ccbeg_pack, ccend_pack,
     *                spherical, .true., iscr,
     *                nints,
     *                nfirst, scr)

           if (nints .gt. 0) then
              if (s .eq. 1) then
                    dd1 = 1
              else
                   dd1 = end_nfps(s-1) + 1
              endif

              dd2 = end_nfps(s)
              call add_4cent_integrals(DSHIELD_YX,
     *                                 a1,a2,b1,b2,c1,
     *                                 c2,d1,d2,scr(nfirst),
     *                                 aa1,aa2,bb1,bb2,cc1,cc2,
     *                                 dd1,dd2,1.d0)
           endif 

      enddo   ! s
      enddo   ! r
      enddo   ! n
      enddo   ! m

C#ifdef _DEBUG_LVL0
      Write(6,*) "The <munu |(r1-r2)*r2/r12| lam sig> integral matrix"
      yx_s = 0.0d0
      yx_l = 0.0d0

      do d = d1,d2
      do c = c1,c2
      do b = b1,b2
      do a = a1,a2

C      Write(6, "(5(1x,F15.10))") (DSHIELD_YX(a,b,c,d), a=a1, a2)
      yx_s = yx_s + DSHIELD_YX(a,b,c,d) * DSHIELD_YX(a,b,c,d)
      yx_l = yx_l + DSHIELD_YX(a,b,c,d)
C      write(*,*)
      enddo 
      enddo 
      enddo 
      enddo 
      Write(6,"(a,a,F20.10,1x,F20.10)")"Check sum of the YX-DSHIELD",
     &                      " integrals = " ,yX_s,yX_l
C#endif

      return
      end

