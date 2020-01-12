      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NDIM=3000)
      DIMENSION Y(NDIM),A(NDIM)
      PARAMETER (HC0=197.3269601D0,AC0=931.494013D0,EC0=137.03599976D0)
      COMMON /PHYSCONST/HC,AC,EC    
c$$$CCC
c$$$      HC=HC0
c$$$      AC=AC0
c$$$      EC=EC0
c$$$CCC
c$$$      AMASS1=1.0086649D0
c$$$      AMASS2=1.0078250D0
c$$$      DREM=AMASS1*AMASS2/(AMASS1+AMASS2)*AC
CCC
CCC     Y"(X) = A(X)Y(X)
CCC
CCC     A(X) = -1 
CCC     - Y(X) = SIN(X) OR COS(X)
CCC     - Y(0) = 0 -> Y(X) = SIN(X) 
CCC     - Y(1) = 1 -> Y(X) = COS(X) 
CCC 
      PI = DACOS(-1D0)
      RAD = 180D0/PI
      DX = 10D0
      DXR= DX/RAD 
      XMAX  = 3600D0
      NXMAX = XMAX/DX+1
      IF(NXMAX.GT.NDIM) THEN
      WRITE(*,*) "DIMENSION ERROR"
      STOP
      END IF
      Y(1)=0D0
      Y(2)=DSIN(DX/RAD)
c$$$      Y(1)=1D0
c$$$      Y(2)=DCOS(DX/RAD)
      A(1)=-1D0
      A(2)=-1D0      
      DO NX = 3,NXMAX
         X  = (NX-1)*DX 
         A(NX) = -1D0
         F1 = (2D0+5D0/6D0*DXR**2*A(NX-1))*Y(NX-1)
         F2 = (1D0-1D0/12D0*DXR**2*A(NX-2))*Y(NX-2)
         F3 = 1D0-1D0/12D0*DXR**2*A(NX)
         Y(NX) = (F1-F2)/F3
         WRITE(*,1000) X, Y(NX), DSIN(X/RAD)
      END DO
 1000 FORMAT(F10.5,2(1X,1PE15.8))
      STOP
      END
