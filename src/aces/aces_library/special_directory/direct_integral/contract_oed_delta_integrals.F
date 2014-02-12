C  Copyright (c) 2003-2010 University of Florida
C
C  This program is free software; you can redistribute it and/or modify
C  it under the terms of the GNU General Public License as published by
C  the Free Software Foundation; either version 2 of the License, or
C  (at your option) any later version.

C  This program is distributed in the hope that it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C  GNU General Public License for more details.

C  The GNU General Public License is included in this distribution
C  in the file COPYRIGHT.
      subroutine contract_oed_delta_integrals(density, z,
     *                       da1, da2, db1, db2)
c---------------------------------------------------------------------------
      implicit none
      include 'int_gen_parms.h'
      include 'machine_types.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      integer type
      integer da1, da2, db1, db2
      double precision density(da1:da2, db1:db2)
      double precision z(da1:da2, db1:db2)
      double precision dbuf(1)
      integer ibuf(1)
#ifdef ALTIX
      pointer (iptr, ibuf)
      pointer (dptr, dbuf)
#else
      common dbuf
      equivalence (ibuf, dbuf)
#endif
      integer i, ix, m, n, nfirst, nints, nalpha, ncoeff, iatom, ncsum
      integer component
      integer a, b, a1,a2,b1,b2,lastao, igrad
      integer m1, m2, n1, n2
      integer der1x, der1y, der1z
      integer der2x, der2y, der2z
      integer dercx, dercy, dercz
      integer ixderc
      integer ind(2)
      integer flopcount
      logical spherical
   
      integer max_dim_coeff
      parameter (max_dim_coeff = 5000)
      integer ccbeg_pack(max_dim_coeff), ccend_pack(max_dim_coeff)
      double precision alpha_pack(max_dim_coeff), 
     *                 pcoeff_pack(max_dim_coeff)
      double precision xn(max_centers), yn(max_centers), 
     *                 zn(max_centers)
      double precision x1,y1,z1,x2,y2,z2, y

c      double precision A_iso(max_centers)
      double precision X_n,Y_n,Z_n

      common /flags/ iflags
      integer iflags(100)


#ifdef ALTIX
      iptr = ishptr
      dptr = dshptr
#endif

      spherical = (ispherical .eq. 1)
      lastao    = end_nfps(nshells)

      ind(1) = 1
      ind(2) = 2

      if (intpkg .ne. flocke_package) then
         print *,'Error: contract_oed_delta_integrals requires ERD ',
     *           'integral package.'
        call abort_job()
      endif

c----------------------------------------------------------------------
C This will zero out the A_sio for each block and will work only for
C segment size 1 (ie. small basis set). This initilization is moved
C to beta.F, Ajith Perera, 03/2013.
CSSS       do m=1,max_centers
CSSS           A_iso(m) = 0.d0
CSSS       end do

c-------------------------------------------------------------------------
c   Find the shells we need to loop over for this block of density data.
c--------------------------------------------------------------------------

      call lookup_shell(end_nfps, nshells, da1, m1)
      call lookup_shell(end_nfps, nshells, da2, m2)
      call lookup_shell(end_nfps, nshells, db1, n1)
      call lookup_shell(end_nfps, nshells, db2, n2)
      
c-------------------------------------------------------------------------
c   Build an array of coordinates per atom.
c-------------------------------------------------------------------------

      do m = 1, nshells
         iatom = atom(m)
         xn(iatom) = dbuf(icoord+3*(m-1))
         yn(iatom) = dbuf(icoord+3*(m-1)+1)
         zn(iatom) = dbuf(icoord+3*(m-1)+2)
      enddo

c-------------------------------------------------------------------------
c   Calculate the integrals over the necessary shell blocks.
c-------------------------------------------------------------------------

       do m = m1, m2
       do n = n1, n2

            ncsum = ncfps(m) + ncfps(n) 
            x1 = dbuf(icoord+3*(m-1))
            y1 = dbuf(icoord+3*(m-1)+1)
            z1 = dbuf(icoord+3*(m-1)+2)
 
            x2 = dbuf(icoord+3*(n-1))
            y2 = dbuf(icoord+3*(n-1)+1)
            z2 = dbuf(icoord+3*(n-1)+2)

            call pack_coeffs_oed(dbuf(ialpha), ixalpha, dbuf(ipcoeff), 
     *                       ixpcoef, ncfps, npfps, m, n,
     *                       alpha_pack, nalpha, pcoeff_pack,
     *                       ncoeff, ibuf(iccbeg), ibuf(iccend), 
     *                       indx_cc,
     *                       ccbeg_pack, ccend_pack,
     *                       max_dim_coeff)

         do ixderc = 1, ncenters 
               X_n=Xn(ixderc)
               Y_n=Yn(ixderc)
               Z_n=Zn(ixderc)

               do b = db1, db2
               do a = da1, da2
                  z(a,b) = 0.
               enddo
               enddo

               call oed__gener_delta_batch(intmax, zmax, nalpha, ncoeff,
     *                 ncsum, ncfps(m), ncfps(n), npfps(m),npfps(n),
     *                 ivangmom(m), ivangmom(n), x1,y1,z1,x2,y2,z2,
     *                 X_n,Y_n,Z_n,alpha_pack,
     *                 pcoeff_pack, ccbeg_pack, ccend_pack,spherical,
     *                 ibuf(ipkgscr),
     *                 nints, nfirst,
     *                 dbuf(dpkgscr))


               if (nints .gt. 0) then
                  a2 = end_nfps(m)
                  if (m .eq. 1) then
                     a1 = 1
                  else
                     a1 = end_nfps(m-1)+1
                  endif

                  b2 = end_nfps(n)
                  if (n .eq. 1) then
                     b1 = 1
                  else
                     b1 = end_nfps(n-1)+1
                  endif
CSSS       write(*,*) 'm=',ivangmom(m),'n=',ivangmom(n)
c              write(*,*) 'a1=',a1,a2,b1,b2
                 call move_integrals2_delta(z, da1,da2,db1,db2,
     *                 dbuf(dpkgscr+nfirst-1),a1,a2,b1,b2)
C
C This need to moved to cwork0222. It does exactly what this does 
C in generic form. Ajith Perera, 03/2013
C
                call move_integrals2_del(y,z,density, da1,da2,db1,db2,
     *                 a1,a2,b1,b2)
c--------------------------------------------------------------------------

c                 call cwork0222(y,
c     *                          density,da1,da2,db1,db2,
c     *                          da1,da2,db1,db2, ind,
c     *                          z,da1,da2,db1,db2,
c     *                          da1,da2,db1,db2,ind(1), ind(2),
c     *                          flopcount)
c---------------------------------------------------------------------------
c   Sum into the gradient.
c---------------------------------------------------------------------------

                A_iso(ixderc) = A_iso(ixderc) + y

               endif ! nints 

            enddo ! ixderc  


         enddo   ! n shells
         enddo   ! m shells

CSSS      do ixderc = 1, ncenters
CSSS        write(*,*) 'Intermediate A-ten=',ixderc,A_iso(ixderc)
CSSS      end do
      return
      end
