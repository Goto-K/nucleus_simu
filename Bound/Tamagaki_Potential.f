CCC   *
CCC   *     VARIATIONAL METHOD WITH THE GAUSSIAN BASIS FUNCTIONS
CCC   *
      IMPLICIT NONE
      INTEGER   I,K,LI,NR
      REAL*8    NDIM,NMAX
      REAL  ::  VC(3) = (/-5.,-230.,2000./)
      REAL  ::  IC(3) = (/2.5,0.942,0.477/)
      REAL  ::  VT(3) = (/-7.5,-67.5,67.5/)
      REAL  ::  IT(3) = (/2.5,1.2,0.477/)
      REAL*8    HC,HC0,AC,AC0,EC,MASS1,MASS2,LS,PI,R0,RMAX
      REAL*8    DX,GR,DFC1,DFC3,SF0,DR,R,A_0,DREM,NRMAX
      REAL*8    WVC,WVT
      PARAMETER (HC0=197.3269601D0,AC=931.494013D0,EC=137.03599976D0)
      !HC=hc,DREM=muc^2
CCC
      PI=ACOS(-1D0)
      HC=HC0*10
CCC   MASS1, MASS2
      MASS1=0.511006D+6
      MASS2=1.0086649D0*AC*1D+6
CCC   REDEUCED MASS
      DREM=MASS1*MASS2/(MASS1+MASS2)
CCC   ANGULAR MOMENTUM 
      LS=0
CCC   RANGE PARAMETER
      R0=0.1D0
      RMAX=20D0
      NMAX=20
      DX=(RMAX/R0)**(1D0/(NMAX-1D0))
CCC     
CCC   -- (2*LS+1)!!, (2*LS+3)!! --
CCC
      DFC1=1D0
      DFC3=1D0
      DO LI=1,2*LS+1,2
      DFC1=LI*DFC1
      END DO
      DFC3=DFC1*(2*LS+3)
CCC     
CCC   -- LS! --
CCC
      SF0=1D0
      DO LI=1,LS
      SF0=LI*SF0
      END DO
CCC   --------------------------------
CCC   --------------------------------   

      DR=0.005D0
      RMAX=10D0
      NRMAX=RMAX/DR+1
      A_0=HC/DREM*EC       !ボーア半径

      DO NR=1,NRMAX
      R=DR*(NR-1)

      WVC=0
      WVT=0

      DO K=1,3
      WVC=WVC+VC(K)*EXP(-(R/(IC(K))**2))
      WVT=WVT+VT(K)*EXP(-(R/(IT(K))**2))
      END DO

      WRITE(*,*) R,WVC,WVT
      END DO

      STOP
      END
