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
      subroutine print_esr_tensors(x, nindex, type, bval,
     *                             eval, bdim, edim)
c---------------------------------------------------------------------------
c   Prints the value of a scalar variable.  The scalar to be printed
c   is defined in the c_result_array field of the op argument.
c----------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'mpif.h'
      include 'trace.h'
      include 'parallel_info.h'
      include 'dbugcom.h'

      double precision x, scalar
      integer nindex, type(*), bval(*), eval(*)
      integer bdim(*), edim(*)
     
      scalar = x

      If (me .eq. 0) Then
          If (scalar .eq. 100.0) then
          call doprint_esr_hpcc_tensor(scalar)
          Endif
      Endif 
 
      return 
      end 

      subroutine doprint_esr_hpcc_tensor(scalar) 

c---------------------------------------------------------------------------

       implicit none
       include 'interpreter.h'
       include 'int_gen_parms.h'
       include 'parallel_info.h'
       include 'atomdata.h'  ! PV for Hyperfine constant
cPV......................................................
       double precision gfactor,abun
       integer niso, isotope, i, Scalar
       double precision g_e,au2ev,ev2Mhz,Beta_e, Beta_N,con,fac,pi
       double precision Eight_pi_by_3,C_light,M_pe,const_Aiso
cPv.......................................................

       Open(Unit=66,File="summary.out")

       write(66,*)
       write(66,*) '-------------------------------------------------'
       Write(66,*) 'Final Aiso in a.u. : Prakash'
        do i =1, ncenters
      write(66,'(A,1x,I3,1x,A,1x,F5.1,2x,A,1x,F24.16)')'atom:',i,'Atomic
     * no:',charge(i),'Aiso = ',A_iso(i)
          end do
       write(66,*) '-------------------------------------------------'

       write(66,*)
cPv.........................................................................

       write(66,*)
       write(66,*)'Final Aiso in MHz  : Prakash'

       au2ev = 27.2113961d+00      ! Conversion from au to eV
       ev2Mhz= 2.41798836d+08      ! Conversion from eV to MHz

       g_e      =  2.002319304386d+00  ! Electron g-factor
       C_light  = 137.0359895d+00      ! Speed of light in au
       M_pe     = 1836.152701d+00      ! Proton-electron mass ratio
       Beta_e   = 1.0d0/(2.0d0*C_light)
       Beta_N   = Beta_e/M_pe          ! Beta_N =1/2*c*Mpe
       pi       = acos(-1.0d0)
       Eight_pi_by_3 = (8.0d0*pi/3.0d0)

       fac = Eight_pi_by_3 * g_e * Beta_e * Beta_N
       con  = fac * au2ev * ev2Mhz

       write(66,*)'Aiso = 8*pi/3  * g_e * g_N * Beta_e * Beta_N *CalcTe'
       Write(66,'(A,1x,F25.10)') 'Conversion from au to eV : ', au2ev
       Write(66,'(A,1x,F25.10)') 'Conversion from eV to Mhz: ', ev2Mhz
       Write(66,'(A,1x,F25.10)') 'g_e                      : ', g_e
       Write(66,'(A,1x,F25.10)') 'Beta_e                   : ', Beta_e
       Write(66,'(A,1x,F25.10)') 'Beta_N                   : ', Beta_N
       Write(66,*) 
       do i = 1, ncenters
        write(66,'(A,1x,I3,2x,A,1x,F10.4,1x,F10.4,1x,F10.4)')
     * 'Atom', i,'Coord =',acenter(i,1),acenter(i,2),acenter(i,3)
        call atom_gfac(charge(i),gfactor,isotope,abun,niso)
        write(66,'(A,2x,F5.1,5x,A,1x,I3)') 'Atomic no:',
     *  charge(i),':  No. of isotopes:',niso
        write(66,'(A,2x,I5,5x,A,1x,F10.4)')'Isotope  :',isotope
     * , ':  Abundance       :',abun

        write(66,'(A,1x,F10.4)') 'Nuclear g-factor  = ',gfactor
        const_Aiso = con * gfactor
        write(66,'(A,1x,F10.4)') 'Conversion factor = ', const_Aiso
        Write(66,'(A,1x,F10.4)') 'Aiso in MHz       = ',const_Aiso*
     *                                                 A_iso(i)
        write(66,*)
      end do
cPV........................................................................
      return
      end

