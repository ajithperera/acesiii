
c The Broyden-Fletcher-Goldfarb-Shanno update.

      SUBROUTINE BFGS_UPDATE(VC,VP,H,SCRATCH,STEP,TBT,NXM6)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      DIMENSION VC(NXM6),H(NXM6*NXM6),SCRATCH(NXM6*NXM6),
     &          TBT(3*NXM6*NXM6),STEP(NXM6),VP(NXM6)

      N2=1
      N3=1+NXM6
      N4=1+2*NXM6
C
      NX2=1+NXM6*NXM6
C
C SCRATCH(N2) = DG
C
      CALL VADD(SCRATCH(N2),VC,VP,NXM6,-1.D0)
C
C TBT(1) =  [1/DG{dot}DX][DGxDG]
C
      CALL MODMATMUL(TBT(1),SCRATCH(N2),SCRATCH(N2),NXM6,1,NXM6,
     $               NXM6,1,NXM6)
      X0=1.D0/xdot(NXM6,SCRATCH(N2),1,STEP,1)
      CALL xscal(NXM6*NXM6,X0,TBT(1),1)
C
C SCRATCH(N3) = HDX, SCRATCH(N4) = DXH, TBT(NX2)= HDX{x}DXH
C
      CALL MODMATMUL(SCRATCH(N3),H,STEP,NXM6,NXM6,1,NXM6,NXM6,1)
      CALL MODMATMUL(SCRATCH(N4),STEP,H,1,NXM6,NXM6,1,NXM6,NXM6)
      CALL MODMATMUL(TBT(NX2),SCRATCH(N3),SCRATCH(N4),NXM6,1,NXM6,
     $               NXM6,1,NXM6)
C
C Z0 = 1/DXHDX, TBT(NX2) = [1/DXHDX][HDX{x}DXH]
C
      Z0=1.D0/xdot(NXM6,STEP,1,SCRATCH(N3),1)
      CALL xscal(NXM6*NXM6,Z0,TBT(NX2),1)
C
C H(updated) = H(old) + [1/DG{dot}DX][DGxDG] - [1/DXHDX][HDX{x}DXH]
C
      CALL VADD(TBT(1),TBT(1),TBT(NX2),NXM6*NXM6,-1.D0)
      CALL VADD(H,H,TBT(1),NXM6*NXM6,1.D0)
C
C DON'T WRITE ENTIRE HESSIAN UNLESS SPECIFICALLY REQUESTED.
C FOR NOW, JUST USE AN IN-CODE PARAMETER.
C




      RETURN
      END
