CCC   *
CCC   *     VARIATIONAL METHOD WITH THE GAUSSIAN BASIS FUNCTIONS
CCC   *
      IMPLICIT NONE
      INTEGER   I,J,K,L,LI,NR,L1,L2
      REAL*8    NDIM,NMAX
      PARAMETER (NDIM=100,LMAX=2)
      REAL      VEC(NDIM,NDIM),ENE(NDIM),VECS(NDIM,NDIM)
      REAL      VEC2(NDIM,NDIM),VEC3(NDIM,NDIM),ENE2(NDIM)
      REAL      VEC4(NDIM,NDIM),VEC5(NDIM,NDIM),VEC6(NDIM,NDIM)
      REAL      VEC7(NDIM,NDIM),VEC8(NDIM,NDIM),ENE3(NDIM)
      REAL      AN(NDIM,NDIM,LMAX,LMAX),AN0(NDIM,NDIM,LMAX,LMAX)
      REAL      H(NDIM,NDIM,LMAX,LMAX),T(NDIM,NDIM,LMAX,LMAX)
      REAL      VS(NDIM,NDIM,LMAX,LMAX),VCX(NDIM,NDIM,LMAX,LMAX)
      REAL      VTX(NDIM,NDIM,LMAX,LMAX)
      REAL      ALP(NDIM),CN(NDIM)
      REAL*8    HC,HC0,AC,AC0,EC,MASS1,MASS2,LS,PI,R0,RMAX
      REAL*8    DX,GR,DFC1,DFC3,SF0,DR,R,A_0,DREM,NRMAX

      REAL  ::  VC(3) = (/-5.,-230.,2000./)
      REAL  ::  IC(3) = (/2.5,0.942,0.477/)
      REAL  ::  VT(3) = (/-7.5,-67.5,67.5/)
      REAL  ::  IT(3) = (/2.5,1.2,0.477/)
c      DIMENSION TEMP1(NDIM,NDIM),TEMP2(NDIM,NDIM)
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
      DO 20 LI=1,2*LS+1,2
      DFC1=LI*DFC1
 20   CONTINUE
      DFC3=DFC1*(2*LS+3)
CCC     
CCC   -- LS! --
CCC
      SF0=1D0
      DO LI=1,LS
      SF0=LI*SF0
      END DO
CCC
CCC   --- CALCULATION OF MATRIX ELEMENTS ---
CCC
      DO 110 I=1,NMAX
C
      GR(I)=R0*DX**(I-1)
      ALP(I)=1D0/GR**2
      CN(I)=DSQRT(2D0**(LS+2)/DFC1
     &     *DSQRT(2D0*ALP(I))**(2*LS+3)/DSQRT(PI))
      !fai_i=Nnr^(l+1)exp[-a_i*r^2]
      DO 120 J=1,I

CCC   --------- NORM MATRIX ELEMENT --------------

      DO L1=0,LMAX
      DO L2=0,LMAX

      IF(L1=0&&L2=0)THEN
      AN(I,J,L1,L2)
     &  =DSQRT(PI)/DSQRT(GR(I)+GR(J))/4/(GR(I)+GR(J))
      ELSE IF(L1=2&&L2=2)THEN
      AN(I,J,L1,L2)
     &  =DSQRT(PI)/DSQRT(GR(I)+GR(J))*15/16/(GR(I)+GR(J))
      ELSE
      AN(I,J,L1,L2)=0
      END IF

      AN(J,I,L1,L2)=AN(I,J,L1,L2)
CCC
      AN0(I,J,L1,L2)=AN(I,J,L1,L2)
      AN0(J,I,L1,L2)=AN(I,J,L1,L2)

CCC   ---------------T(I,J,L1,L2)-----------------------

      IF(L1=0&&L2=0)THEN
      T(I,J,L1,L2)
     & =-HC**2/DREM*DSQRT(PI)/DSQRT(GR(I)+GR(J))/4/(GR(I)+GR(J))
     & *(-6)*GR(I)*GR(J)/(GR(I)+GR(J))
      ELSE IF(L1=2&&L2=2)THEN
      T(I,J,L1,L2)
     & =-HC**2/DREM*DSQRT(PI)/DSQRT(GR(I)+GR(J))*15/16/(GR(I)+GR(J))
     & *(-14)*GR(I)*GR(J)/(GR(I)+GR(J))
      ELSE
      T(I,J,L1,L2)=0
      END IF

      T(J,I,L1,L2)=T(I,J,L1,L2)

CCC   ---------------VCX(I,J,L1,L2)-----------------
CCC   ---------------VCX=ΣVC*exp[-IC*r^2]-----------

      VCX(I,J,L1,L2)=0

      DO K=1,3

      IF(L1=0&&L2=0)THEN
      VCX(I,J,L1,L2)=VCX(I,J,L1,L2)
     & +VC(K)*DSQRT(PI)/DSQRT(GR(I)+GR(J)+IC(K))/4/(GR(I)+GR(J)+IC(K))
      ELSE IF(L1=2.and.L2=2)THEN
      VCX(I,J,L1,L2)=VCX(I,J,L1,L2)
     & +VC(K)*DSQRT(PI)/DSQRT(GR(I)+GR(J)+IC(K))*15/16/(GR(I)+GR(J)+IC(K))
      ELSE IF(L1==0&&L2==2.or.L1=&&L2=2)THEN
      VCX(I,J,L1,L2)=VCX(I,J,L1,L2)
     & +VC(K)*DSQRT(PI)/DSQRT(GR(I)+GR(J)+IC(K))*3/8/(GR(I)+GR(J)+IC(K))
      ELSE
      VCX(I,J,L1,L2)=0 
      END IF

      END DO

      VCX(J,I,L1,L2)=VCX(I,J,L1,L2)

CCC   ---------------VTX(I,J,L1,L2)-----------------
CCC   ---------------VTX=ΣVT*exp[-IT*r^2]-----------

      VTX(I,J,L1,L2)=0

      DO K=1,3

      IF(L1=0&&L2=0)THEN
      VTX(I,J,L1,L2)=VTX(I,J,L1,L2)
     & +VT(K)*DSQRT(PI)/DSQRT(GR(I)+GR(J)+IT(K))/4/(GR(I)+GR(J)+IT(K))
      ELSE IF(L1=2&&L2=2)THEN
      VTX(I,J,L1,L2)=VTX(I,J,L1,L2)
     & +VT(K)*DSQRT(PI)/DSQRT(GR(I)+GR(J)+IT(K))*15/16/(GR(I)+GR(J)+IT(K))
      ELSE IF(L1=0&&L2=2||L1=&&L2=2)THEN
      VTX(I,J,L1,L2)=VTX(I,J,L1,L2)
     & +VT(K)*DSQRT(PI)/DSQRT(GR(I)+GR(J)+IT(K))*3/8/(GR(I)+GR(J)+IT(K))
      ELSE
      VTX(I,J,L1,L2)=0 
      END IF

      END DO

      VTX(J,I,L1,L2)=VTX(I,J,L1,L2)

CCC   ---------------VLSX(I,J,L1,L2)-----------------

      VLSX(I,J,L1,L2)=0

CCC   ---------------VS(I,J,L1,L2)------------------

      IF(L1=0&&L2=0)THEN
      VS(I,J,L1,L2)
     & =VCX(I,J,L1,L2)
      ELSE IF(L1=2&&L2=2)THEN
      VS(I,J,L1,L2)
     & =VCX(I,J,L1,L2)-2*VTX(I,J,L1,L2)-3*VLSX(I,J,L1,L2)
      ELSE IF(L1=0&&L2=2||L1=&&L2=2)THEN
      VS(I,J,L1,L2)
     & =SQRT(8)*VTX(I,J,L1,L2)
      ELSE
      VS(I,J,L1,L2)=0 
      END IF

      VS(J,I,L1,L2)=VS(I,J,L1,L2)

      H(I,J,L1,L2)=T(I,J,L1,L2)+VX(I,J,L1,L2)+VS(I,J,L1,L2)

CCC   ----------------END DO----------------------
      END DO
      END DO
 120  CONTINUE
 110  CONTINUE
CCC   ------------N * VEC = ENE * VEC-------------

      CALL HQR(NDIM,NMAX,AN,VEC,ENE)
CCC
CCC   ----------------ΣENE = sum-------------------
CCC
      sum = 0
 1100 FORMAT('a:Eigen value  =',I3,F10.5)
 1200 FORMAT('x:Eigen vector = (',F10.5,9(1X,F10.5),')')
      DO I=1,NMAX
c      WRITE(*,1100) I,ENE(I)
c      WRITE(*,1200) (VEC(J,I),J=1,NMAX)
c      WRITE(*,*)
      sum = sum + ENE(I)
      END DO
CCC
c      WRITE(*,*) sum

CCC
CCC   ------UNITARY MATRIX(VEC * VEC' = I)--------
CCC
      DO I=1,NMAX
      DO J=1,NMAX
      VECS(I,J) = 0
      DO K=1,NMAX
      VECS(I,J) = VECS(I,J) + VEC(I,K) * VEC(J,K)
      END DO
      END DO
      END DO

 1300 FORMAT(F10.5,9(1X,F10.5))
      DO I=1,NMAX
c      WRITE(*,1300) (VECS(J,I),J=1,NMAX)
      END DO

CCC   -------------VEC2 = VEC / √ENE-------------
CCC   -----------基底ベクトルの規格化----------------

      DO I=1,NMAX
      DO J=1,NMAX
      VEC2(I,J) = VEC(I,J)/DSQRT(ENE(J))
      END DO
      END DO

CCC   --------VEC2' * N * VEC2 = VEC3 = I--------
CCC   -------------------------------------------
      
      DO I=1,NMAX
      DO J=1,NMAX
      VEC3(I,J) = 0d0
      DO K=1,NMAX
      DO L=1,NMAX
      VEC3(I,J) = VEC3(I,J)+VEC2(K,I)*AN0(K,L)*VEC2(L,J)
      END DO
      END DO
      END DO
      END DO

      DO I=1,NMAX
c      WRITE(*,1300) (VEC3(J,I),J=1,NMAX)
      END DO

CCC   --------VEC2' * H * VEC2 = VEC4------------

      DO I=1,NMAX
      DO J=1,NMAX
      VEC4(I,J) = 0d0
      DO K=1,NMAX
      DO L=1,NMAX  
      VEC4(I,J) = VEC4(I,J)+VEC2(K,I)*H(K,L)*VEC2(L,J)
      END DO
      END DO
      END DO
      END DO

CCC   ---------VEC4 * VEC5 = ENE2 * VEC5---------

      CALL HQR(NDIM,NMAX,VEC4,VEC5,ENE2)
      DO I=1,NMAX
c      WRITE(*,*) I,ENE2(I)
      END DO

      DO I=1,NMAX
      DO K=1,NMAX
      VEC6(I,K) = 0d0
      DO J=1,NMAX 
      VEC6(I,K) = VEC6(I,K)+VEC5(J,K)*VEC2(I,J)
      END DO
      END DO
      END DO

CCC   --------------基底状態--------------------

      DR=0.1D0
      RMAX=10D0
      NRMAX=RMAX/DR+1
      K=5
      A_0=HC/DREM*EC       !ボーア半径

      DO NR=1,NRMAX
      R=DR*(NR-1)
      WF=0D0
      DO I=1,NMAX
      WF=WF+VEC6(I,K)*CN(I)*R**(LS+1)*DEXP(-ALP(I)*R**2)
      END DO
      WRITE(10,*) R,WF,2D0*R*DEXP(-R/A_0)/DSQRT(A_0)**3
      END DO

      DO K=1,NMAX
      I=0
      DO J=1,NMAX
      IF(ENE2(K).GE.ENE2(J)) I=I+1
      END DO
      ENE3(I)=ENE2(K)
      DO J=1,NMAX
      VEC7(J,I)=VEC6(J,K)
      END DO
      END DO

      DO I=1,10
c      WRITE(*,*) I,ENE3(I)
c      WRITE(*,*) (VEC7(J,I),J=1,NMAX)
c      WRITE(*,*)
      WRITE(*,*) I,ENE3(I),(2*(I-1)+LS+3D0/2)
      END DO

      STOP
      END

CCC      
CCC   --------------------------------------------
CCC   --------------------------------------------
CCC   ハウスホルダー + QR法による対角化
CCC
      SUBROUTINE HQR(ND,N,G,X,E)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION G(ND,ND),X(ND,ND),E(ND)
      DIMENSION WA(0:ND),WB(0:ND),WC(0:ND)
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
         SB=DSQRT(SA)
         A=G(K1,K)
         IF(A.GE.UZ) SB=-SB
         WA(K)=SB
         IF(SA.LE.DABS(AC*G(K,K))) GO TO 3
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
         T=DABS(E(K))+TA
         TA=DABS(WA(K))
         T=T+TA
         AN=MAX(AN,T)
 17   CONTINUE
      EPS=AN*AC
      AM=UZ
      M=N
 15   IF(M.EQ.0) GO TO 21
      I=M-1
      K=I
      IF(DABS(WA(K)).GT.EPS) GO TO 16
      AM=UZ
      M=K
      GO TO 15
 16   I=I-1
      IF(DABS(WA(I)).LE.EPS.OR.I.EQ.0) GO TO 18
      K=I
      GO TO 16
CCC   原点シフト SH
 18   SH=UZ
      IF(DABS(E(M)-AM).LT.UH*DABS(E(M)).OR.
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
         SS=DSQRT(P**2+CC**2)
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

