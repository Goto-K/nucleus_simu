CCC      
CCC   --------------------------------------------
CCC   --------------------------------------------
CCC   ハウスホルダー + QR法による対角化
CCC
      subroutine hqr(ND,N,G,X,E)

      implicit none
      integer   i,i1,j,k,k1,l,nd,n,n1,m
      complex*16,dimension(nd,nd)  :: 
     &          G,X
      complex*16,dimension(0:nd)  :: 
     &          WA,WB,WC
      complex*16,dimension(nd)  :: 
     &          E
      real*8    sa,sb,eps,an,a,t
      complex*16   iv,ac,uz,ui,uh,ta,am,sh,p,
     &             pc,pd,ss,cc,sq,xa,xb,sn,sum
      complex*16,parameter :: ci=cmplx(0d0,1d0)

CCC
      IV=1
CCC
      AC=1.0D-14
      UZ=0.0D0
      UI=1.0D0
      UH=0.5D0
      WA(0)=UZ
      IF(N.LE.1) THEN
         X(1,1)=UI
         E(1)=G(1,1)
         RETURN
      ENDIF
      IF(IV.NE.0) THEN
         DO 1 I=1,N
         DO 1 J=1,N
           X(J,I)=UZ
 1      CONTINUE
        DO 2 I=1,N
           X(I,I)=UI
 2      CONTINUE
      ENDIF
      IF(N.EQ.2) GO TO 5
      DO 3 K=1,N-2
         E(K)=G(K,K)
         K1=K+1
         SA=UZ
         DO 4 I=K1,N
            SA=SA+G(I,K)**2
 4       CONTINUE
         SB=SQRT(SA)
         A=G(K1,K)
         IF(A.GE.real(UZ)) SB=-SB
         WA(K)=SB
         IF(SA.LE.CDABS(AC*G(K,K))) GO TO 3
         SN=UI/(SA-A*SB)
         G(K1,K)=A-SB
         DO 6 I=K1,N
            SUM=UZ
            DO 7 J=K1,I
               SUM = SUM+G(I,J)*G(J,K)
 7          CONTINUE
            IF(I.EQ.N) GO TO 6
            I1=I+1
            DO 8 J=I1,N
               SUM = SUM+G(J,I)*G(J,K)
 8          CONTINUE
 6          WB(I)=SUM*SN
         SUM=UZ
         DO 9 I=K1,N
            SUM=SUM+G(I,K)*WB(I)
 9       CONTINUE
         T=UH*SN*SUM
         DO 10 I=K1,N
            WB(I)=WB(I)-T*G(I,K)
 10      CONTINUE
         DO 11 I=K1,N
         DO 11 J=K1,I
            G(I,J)=G(I,J)
     1           -G(I,K)*WB(J)-WB(I)*G(J,K)
 11      CONTINUE
CCC   固有ベクトルの処理
         IF(IV.NE.0)THEN
            DO 12 I=2,N
               SUM =UZ
               DO 13 J=K1,N
                  SUM=SUM+X(I,J)*G(J,K)
 13            CONTINUE
               WB(I)=SN*SUM
 12         CONTINUE
            DO 14 J=K1,N
            DO 14 I=2,N
               X(I,J)=X(I,J)-WB(I)*G(J,K)
 14         CONTINUE
         END IF
 3    CONTINUE
 5    CONTINUE
C
      N1=N-1
      E(N1)=G(N1,N1)
      E(N)=G(N,N)
      WA(N1)=G(N,N1)
CCC   三重対角行列の対角化
      AN=UZ
      TA=UZ
      WA(N)=UZ
      DO 17 K=1,N
         T=CDABS(E(K))+TA
         TA=CDABS(WA(K))
         T=T+TA
         AN=MAX(AN,T)
 17   CONTINUE
      EPS=AN*AC
      AM=UZ
      M=N
 15   IF(M.EQ.0) GO TO 21
      I=M-1
      K=I
      IF(CDABS(WA(K)).GT.EPS) GO TO 16
      AM=UZ
      M=K
      GO TO 15
 16   I=I-1
      IF(CDABS(WA(I)).LE.EPS.OR.I.EQ.0) GO TO 18
      K=I
      GO TO 16
CCC   原点シフト SH
 18   SH=UZ
      IF(CDABS(E(M)-AM).LT.real(UH*CDABS(E(M))).OR.
     1     M.EQ.K+1) SH=E(M)+UH*WA(M-1)
      AM=E(M)
      WB(K-1)=UI
      WC(K-1)=UZ
      P=E(K)-SH
      DO 19 I=K,M
         I1=I-1
         IF(I.GE.K+1) P=(E(I)-SH)*WB(I1)
     1        -WA(I1)*WC(I1)*WB(I-2)
         PD=PC
         PC=P*WB(I1)
         CC=WA(I)
         SQ=SS
         SS=DSQRT(real(P**2+CC**2))
         WB(I)=P/SS
         WC(I)=CC/SS
         IF(I.GE.K+1) E(I-1)=E(I)+PD-PC
         IF(I.GE.K+2) WA(I-2)=WC(I-2)*SQ
 19   CONTINUE
      E(M)=PC+SH
      WA(M-1)=WC(M-1)*P
CCC   固有ベクトルの処理
      IF(IV.NE.0) THEN
         DO 20 I=K,M-1
         DO 20 J=1,N
            XA=X(J,I)
            XB=X(J,I+1)
            X(J,I)=XA*WB(I)+XB*WC(I)
            X(J,I+1)=-XA*WC(I)+XB*WB(I)
 20      CONTINUE
      END IF
      GO TO 15
 21   RETURN
      END

