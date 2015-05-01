       Subroutine debug(PRIMA, PRIMB, PRIMC, PRIMD, N, M)

       Implicit none

       INTEGER     PRIMA (1:N)
       INTEGER     PRIMB (1:N)
       INTEGER     PRIMC (1:M)
       INTEGER     PRIMD (1:M)

       INTEGER N, M, ij


       Write(6,*) (PRIMA(ij), ij=1,n)
       Write(6,*) (PRIMB(ij), ij=1,n)
       Write(6,*) (PRIMC(ij), ij=1,m)
       Write(6,*) (PRIMD(ij), ij=1,m)

       Return
       end
   
