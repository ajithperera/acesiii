      subroutine print_cints(cints, NXYZET, NCGTOAB)

      Implicit double precision (A-H, O-Z)

      Dimension Cints(NXYZET, NCGTOAB)
     
      Call output(Cints, 1, NXYZET, 1, NCGTOAB, NXYZET, NCGTOAB, 1)

      Return 
      End
       
    
      
