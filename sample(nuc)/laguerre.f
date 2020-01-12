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
