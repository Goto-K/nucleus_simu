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
