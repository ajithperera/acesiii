      subroutine compute_return_angmomy(out_ly,da1,da2,db1,db2,nsend)

      implicit none

      include 'mpif.h'
      include 'int_gen_parms.h'
      include 'machine_types.h'
      include 'dbugcom.h'
      include 'parallel_info.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

c--------------------------------------------------------------------------
      integer nsend
      integer da1, da2, db1, db2
      integer adim, bdim
      double precision out_ly(da1:da2, db1:db2)
c.........................................................................
      double precision dbuf(1), temp
      integer ibuf(1)
#ifdef ALTIX
      pointer (iptr, ibuf)
      pointer (dptr, dbuf)
#else
      common dbuf
      equivalence (ibuf(1), dbuf(1))
#endif

      logical spherical

      double precision x1,y1,z1
      double precision x2,y2,z2
      double precision xc,yc,zc

cPV we need this dipole derivative integrals 
      integer iatom
      integer DERCX,DERCY,DERCZ
      integer DER1X, DER1Y,DER1Z,IXDERC,DER2X,DER2Y,DER2Z
      integer XMOMDER, YMOMDER, ZMOMDER
      integer a,b,m,n,m1,m2,n1,n2
      integer Ea1,Ea2,Eb1,Eb2
c end PV
      integer nalpha,ncoeff,nints,nfirst
      integer ncsum
      integer max_dim_coeff
      parameter (max_dim_coeff = 5000)
      integer ccbeg_pack(max_dim_coeff), ccend_pack(max_dim_coeff)
      double precision alpha_pack(max_dim_coeff), 
     *                 pcoeff_pack(max_dim_coeff)

      spherical = (ispherical .eq. 1)
  
      adim = da2-da1+1
      bdim = db2-db1+1

      nsend = adim*bdim
      if (nsend .lt. 0) then
         print *,'ERROR IN INTEGRAL WORKER ',me,' nsend = ',nsend
         print *,'adim,bdim = ',adim,bdim
         call Errex
      endif

c----------------------------------------------------------------------------
c   Clear the output array.
c----------------------------------------------------------------------------

      do b = db1, db2
      do a = da1, da2
         out_ly(a,b) = 0.d0
      end do
      end do
C
c-----------------------------------------------------------------------
c   Find the shell blocks for which we shall loop through.
c-----------------------------------------------------------------------

       call lookup_shell(end_nfps, nshells, da1, m1)
       call lookup_shell(end_nfps, nshells, da2, m2)
       call lookup_shell(end_nfps, nshells, db1, n1)
       call lookup_shell(end_nfps, nshells, db2, n2)
C
C Set the center ar which the angular momentum integrals are
C evaluated. Most cases it is at the origin (0.0,0.0,0.0),
C but it is not necessary.

      XC = 0.0D0
      YC = 0.0D0
      ZC = 0.0D0
C
C The Y component of the angular momentum; -{<mu|(y-C)|derz(nu)> -
C <mu|(z-C)|derx(nu)>}
C
      do m = m1, m2

         x1 = dbuf(icoord+3*(m-1))
         y1 = dbuf(icoord+3*(m-1)+1)
         z1 = dbuf(icoord+3*(m-1)+2)

         do n = n1, n2

            x2 = dbuf(icoord+3*(n-1))
            y2 = dbuf(icoord+3*(n-1)+1)
            z2 = dbuf(icoord+3*(n-1)+2)

            ncsum = ncfps(m) + ncfps(n) 

            call pack_coeffs_oed(dbuf(ialpha), ixalpha, dbuf(ipcoeff),
     *                          ixpcoef,ncfps, npfps,
     *                          m, n,
     *                          alpha_pack, nalpha,
     *                          pcoeff_pack, ncoeff,
     *                          ibuf(iccbeg),ibuf(iccend), indx_cc,
     *                          ccbeg_pack, ccend_pack,max_dim_coeff)

CSSS      Write(6,*) "@-entry to OED_batch"
CSSS      Write(6,"((a,2(1x,i4)))") "m, n :", m, n

      DER1X = 0
      DER1Y = 0
      DER1Z = 0
      DER2X = 0
      DER2Y = 0
      DER2Z = 1

      XMOMDER = 1
      YMOMDER = 0
      ZMOMDER = 0
  
      call oed__gener_xyz_derv_batch(intmax,zmax,
     *                 nalpha,ncoeff,
     *                 ncsum, ncfps(m), ncfps(n), npfps(m),npfps(n),
     *                 ivangmom(m),ivangmom(n),x1,y1,z1,x2,y2,z2,
     *                 der1x,der1y,der1z,        
     *                 der2x,der2y,der2z,
     *                 xmomder,ymomder,zmomder,
     *                 xc,yc,zc,alpha_pack, 
     *                 pcoeff_pack, ccbeg_pack, ccend_pack,
     *                 spherical, .true., ibuf(ipkgscr),
     +                 nints, 
     +                 nfirst,
     +                 dbuf(dpkgscr))

      if(nints .gt. 0) then
         Ea2 = end_nfps(m)
         if(m .eq. 1) then
            Ea1=1
         else
            Ea1 = end_nfps(m-1)+1
         end if
         Eb2 = end_nfps(n)
         if(n .eq. 1) then
           Eb1 = 1
         else
           Eb1 = end_nfps(n-1)+1
         end if

      call add_integrals2(out_ly,da1,da2,db1,db2,
     *                   dbuf(dpkgscr+nfirst-1),Ea1,Ea2,Eb1,Eb2,-1.0d0)
      endif 

      DER1X = 0
      DER1Y = 0
      DER1Z = 0
      DER2X = 1
      DER2Y = 0
      DER2Z = 0

      XMOMDER = 0
      YMOMDER = 0
      ZMOMDER = 1

      call oed__gener_xyz_derv_batch(intmax,zmax,
     *                 nalpha,ncoeff,
     *                 ncsum, ncfps(m), ncfps(n), npfps(m),npfps(n),
     *                 ivangmom(m),ivangmom(n),x1,y1,z1,x2,y2,z2,
     *                 der1x,der1y,der1z,
     *                 der2x,der2y,der2z,
     *                 xmomder,ymomder,zmomder, 
     *                 xc,yc,zc,alpha_pack, 
     *                 pcoeff_pack, ccbeg_pack, ccend_pack,
     *                 spherical, .true., ibuf(ipkgscr),
     +                 nints, 
     +                 nfirst,
     +                 dbuf(dpkgscr))

      if(nints .gt. 0) then
         Ea2 = end_nfps(m)
         if(m .eq. 1) then
            Ea1=1
         else
            Ea1 = end_nfps(m-1)+1
         end if
         Eb2 = end_nfps(n)
         if(n .eq. 1) then
           Eb1 = 1
         else
           Eb1 = end_nfps(n-1)+1
         end if

      call add_integrals2(out_ly,da1,da2,db1,db2,
     *                   dbuf(dpkgscr+nfirst-1),Ea1,Ea2,Eb1,Eb2,1.0d0)
       endif ! nints

       end do !n
       end do !m

       if (dbg) write(6,*) " Finished Ly angular mom. integral comp."

#ifdef _DEBUG_LVL0
      Write(6,*) "The <mu| L_y |nu > integral matrix"
      do b = db1, db2
         Write(6,*)
         Write(6, "(5(1x,F15.10))") (out_ly(a,b), a=da1, da2)
      enddo
#endif

      return
      end

