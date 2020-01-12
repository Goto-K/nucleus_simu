      IMPLICIT REAL*8(A-H,O-Z)
      parameter (NRDIM=1000)
      DIMENSION WF(NRDIM,3,3)
      PARAMETER (HC0=197.3269601D0,AC0=931.494013D0,EC0=137.03599976D0)
CCC
      HC=HC0
      AC=AC0
      EC=EC0
CCC
      AMASS1=1.0086649D0
      AMASS2=1.0078250D0
      DREM=AMASS1*AMASS2/(AMASS1+AMASS2)*AC
CCC
      DR=0.1D0
      NRMAX = 301
      ANUE = DSQRT(DREM/HC**2)
      OPEN(10,FILE="wf-ho-ex.d")
      DO N=0,2
         FACT=FACT*N
         IF(N.EQ.0) FACT=1D0
      DO LS=0,2
      GAMMA1 = GAMMAF(N+LS+1,2)
      DO NR=1,NRMAX
      R=(NR-1)*DR
      Z=ANUE**2*R**2
      CALL LAGUERRE(N,LS,Z,FLAC)
      C0 = DSQRT(2d0*ANUE**(2*LS+3)*FACT)/DSQRT(GAMMA1)
      WF(NR,LS+1,N+1) = C0*R**(LS+1)*DEXP(-ANUE**2*R**2/2D0)*FLAC
      END DO
      END DO
      END DO
CCC
      WRITE(*,110) ((LS,N,LS=0,2),N=0,2)
      DO NR=1,NRMAX
      R=(NR-1)*DR
      WRITE(*,100) R,((WF(NR,LS+1,N+1),LS=0,2),N=0,2)
      END DO
 100  FORMAT(F8.3,4X,10(1x,1PE12.5))
 110  FORMAT('#',3X,'R',3X,'(LS,N)=',1X,9('(',I2,',',I2,')',6x))
CCC      
      STOP
      END

      SUBROUTINE LAGUERRE(N,LS,Z,FLAG)
      IMPLICIT REAL*8(A-H,O-Z)
CCC
CCC   かなり適当です. N, Zが大きいときはダメかも.
CCC   恐らく EXP(-0.5*Z)を付けたもので定義した方が良い.
CCC
      MAXDF = 2*N+2*LS+1
      DFACT1=1D0
      DO I=1,MAXDF,2
      DFACT1 = DFACT1*I
      END DO
      FAC1 = DFACT1/2D0**N
CCC      
      FLAG=0D0
      DO KAP=0,N
      MAXDF2 = 2*KAP+2*LS+1
      DFACT2=1D0
      DO I=1,MAXDF2,2
      DFACT2=DFACT2*I
      END DO
CCC
      FACTK=1D0
      DO I=1,KAP
      FACTK=FACTK*I
      END DO
CCC
      FACTKN=1D0
      DO I=1,N-KAP
      FACTKN=FACTKN*I
      END DO
CCC
      FACZ = (-2D0*Z)**KAP
      FLAG = FLAG + FAC1/DFACT2/FACTK/FACTKN*FACZ
      END DO

      RETURN
      END

CCC
      FUNCTION GAMMAF(N,N0)
      IMPLICIT REAL*8(A-H,O-Z)
CCC   N IS A POSITIVE INTEGER
CCC 
CCC   -> IF N0 = 1, GAMMA(N) = (N-1)!
CCC   -> IF N0 = 2, GAMMA(N + 1/2) = (2N)!/(4**N*N!)*DSQRT(PI)
CCC 
      PI=DACOS(-1D0)
CCC   IF N0 = 1, GAMMA(N) = (N-1)!
      IF(N0.EQ.1) THEN
      FACT = 1D0      
      IF(N.EQ.1) THEN
      FACT=1D0
      ELSE
      DO I=1,N-1
      FACT=FACT*I
      END DO
      ENDIF
CCC   IF N0 = 2, GAMMA(N + 1/2) = (2N)!/(4**N*N!)*DSQRT(PI)
      ELSE IF(N0.EQ.2) THEN
      FACT1=1D0
      FACT2=1D0
      N2=N*2
c
      DO I=1,N
      FACT1=FACT1*I
      END DO
C
      DO I=1,N2
      FACT2=FACT2*I
      END DO
c
      FACT=FACT2*DSQRT(PI)/FACT1/4D0**N
      ENDIF
CCC   ---------------------------------------
      GAMMAF=FACT
      RETURN
      END
