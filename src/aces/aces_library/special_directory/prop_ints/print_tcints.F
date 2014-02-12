      subroutine print_tcints(cints, NCGTOAB, NXYZET)

      Implicit double precision (A-H, O-Z)

      Dimension Cints(NCGTOAB, NXYZET)

      Call output(Cints, 1, NCGTOAB, 1, NXYZET, NCGTOAB, NXYZET, 1)
    
      Return 
      End
       
    
      
