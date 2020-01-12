CCC
CCC
CCC
      SUBROUTINE LNEQ(ND,N,A,C)
      IMPLICIT REAL*8(A-H,O-Z)
C      COMPLEX*16 A,C,AM,BM
      DIMENSION A(ND,ND),C(ND),NW(ND)
      EPS = 1D-75
      UZ=0D0
      NW(1)=1
      IF(N.LT.1) RETURN
      IF(N.EQ.1) THEN
         C(1)=C(1)/A(1,1)
         RETURN
      END IF
      DO 1 I=1,N
         NW(I)=I
 1    CONTINUE
      N1=N-1
      DO 2 I=1,N1
         KM=NW(I)
         K=I
         AM=A(I,I)
         I1=I+1
         
         DO 3 J=I1,N
            IF(ABS(A(I,J)).GT.ABS(AM)) THEN
               K=J
               AM=A(I,J)
            ENDIF
 3       CONTINUE
         IF(K.NE.I) THEN
            NW(I)=NW(K)
            NW(K)=KM
            DO 5 J=1,N
               BM=A(J,I)
               A(J,I)=A(J,K)
               A(J,K)=BM
 5          CONTINUE
         END IF
         IF(ABS(AM).LE.ABS(C(I))*EPS) GO TO 24
         DO 8 J=I,N
            A(I,J)=A(I,J)/AM
 8       CONTINUE
         C(I)=C(I)/AM
         DO 7 J=I1,N
            C(J)=C(J)-C(I)*A(J,I)
 7       CONTINUE
         DO 9 K=I1,N
            DO 9 J=I1,N
               A(J,K)=A(J,K)-A(I,K)*A(J,I)
 9          CONTINUE
            DO 11 J=I1,N
               A(J,I)=UZ
 11         CONTINUE
 2       CONTINUE
         IF(ABS(A(N,N)).LE.ABS(C(N))*EPS) GO TO 24
         C(N)=C(N)/A(N,N)
         DO 30 I=1,N1
            I2=N-I
            AM=UZ
            DO 31 J=1,I
               N2=N-J+1
               AM = AM + A(I2,N2)*C(N2)
 31         CONTINUE
            C(I2)=C(I2)-AM
 30      CONTINUE
         DO 19 I=1,N
            K=NW(I)
            A(1,K)=C(I)
 19      CONTINUE
         DO 20 I=1,N
            C(I)=A(1,I)
 20      CONTINUE
         NW(1)=N
 18      CONTINUE
         RETURN
 24      NW(1)=-1
         RETURN
         END
