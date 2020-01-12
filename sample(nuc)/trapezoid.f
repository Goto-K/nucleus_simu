      IMPLICIT REAL*8(A-H,O-Z)
CCC
      PI=DACOS(-1D0)
CCC
      XMAX=PI
      NXMAX=1000
      DX=PI/NXMAX
CCC
      SUM=0D0
      DO 10 NX=2,NXMAX
         X=DX*(NX-1)
         FX=DSIN(X)
         SUM=SUM+FX*DX
 10   CONTINUE
      SUM=SUM+(DSIN(0D0)+DSIN(PI))*DX/2D0
CCC
      WRITE(*,100) SUM
 100  FORMAT('INTEGRAL OF SIN(X) [0:PI]:',F10.5)

      STOP
      END
