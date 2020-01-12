CCC   *
CCC   *     VARIATIONAL METHOD WITH THE GAUSSIAN BASIS FUNCTIONS
CCC   *
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (MDIM=100)
      DIMENSION VEC(MDIM,MDIM),ENE(MDIM)
      DIMENSION XN(MDIM),CN(MDIM)
      PARAMETER (HC0=197.3269601D0,AC0=931.494013D0,EC0=137.03599976D0)
      COMMON /PHYSCONST/HC,AC,EC    
CCC
      HC=HC0
      AC=AC0
      EC=EC0
CCC
      AMASS1=1.0086649D0
      AMASS2=1.0078250D0
      DREM=AMASS1*AMASS2/(AMASS1+AMASS2)*AC
CCC   ANGULAR MOMENTUM AND POTENTIAL PARAMETER
      LS=0
      VG0=-72.15D0
      RG0=1.484D0
CCC   RANGE PARAMETER
      R0=0.1D0
      RMAX=50D0
      NMAX=30
      DX=(RMAX/R0)**(1D0/(NMAX-1D0))
CCC      
      CALL GAUSD(R0,DX,DREM,LS,MDIM,NMAX,VG0,RG0,
     1     ENE,VEC,CN,XN)
CCC
      DO I=1,NMAX
      WRITE(*,'(I3,1X,1PE15.8)') I,ENE(I)
      END DO
CCC
CCC
      STOP
      END


CCC   ---------------------------------------------------
CCC   *                                                 *
CCC   ***** DIAGONALIZATION USING NORMAL GAUSS BASE *****
CCC   *                                                 *
CCC   ---------------------------------------------------
CCC
CCC   2011.03 RENEWAL
CCC
      SUBROUTINE GAUSD(XSTART,DX,DREM,LS,MDIM,MMAX,VG0,RG0,
     1     W,AZ,CN,XN)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION ALP(MDIM)
      DIMENSION A(MDIM,MDIM),B(MDIM,MDIM)
      DIMENSION AZ(MDIM,MDIM),W(MDIM)
      DIMENSION XN(MDIM),CN(MDIM)
      DIMENSION TEMP1(MDIM,MDIM),TEMP2(MDIM,MDIM),WORK(MDIM-1)
      COMMON /PHYSCONST/HC,AC,EC    
C
      PI=DACOS(-1D0)
CCC     
CCC   -- (2*LS+1)!!, (2*LS+3)!! --
CCC
      DKL1=1D0
      DKL2=1D0
      DO 20 LI=1,2*LS+1,2
      DKL1=LI*DKL1
 20   CONTINUE
      DKL2=DKL1*(2*LS+3)
CCC
CCC   -- (LS+1)! --
CCC
      QLS1=1D0
      DO 30 LI=1,LS+1
      QLS1=QLS1*LI
 30   CONTINUE
CCC
CCC   --- CALCULATION OF MATRIX ELEMENTS ---
CCC
      DO 110 I=1,MMAX
C
      XN(I)=XSTART*DX**(I-1)
      ALP(I)=1D0/XN(I)**2
      CN(I)=DSQRT(2D0**(LS+2)/DKL1
     &     *DSQRT(2D0*ALP(I))**(2*LS+3)/DSQRT(PI))
      DO 120 J=1,I
CCC
CCC   --------- CALCULATE NORM MATRIX ELEMENT --------------
CCC
      B(I,J)=CN(I)*CN(J)*DKL1/2D0**(LS+2)
     &     *DSQRT(PI)/DSQRT(ALP(I)+ALP(J))**(2*LS+3)
      B(J,I)=B(I,J)
CCC
CCC   --------- CALCULATE HAMILTONIAN MATRIX ELEMENT --------------
CCC
CCC   --- KINETIC TERM ---
      HT=CN(I)*CN(J)*(DKL2/2**(LS+3)*HC**2/DREM*2D0
     &     *ALP(I)*ALP(J)
     &     *DSQRT(PI)/DSQRT((ALP(I)+ALP(J)))**(2*LS+5))
C
CCC   ***** <<< CENTRAL TERM >>> *****
C
      HVC=0D0
      HVC = CN(I)*CN(J)*VG0*DKL1/2D0**(LS+2)*DSQRT(PI)
     &        /DSQRT((ALP(I)+ALP(J)+1D0/RG0**2))**(2*LS+3)
CCC   --- TOTAL HAMILTONIAN MATRIX ---
      CH1 = HT + HVC
CCC
      A(I,J)=CH1
      A(J,I)=A(I,J)
CCC
 120  CONTINUE
 110  CONTINUE
CCC
CCC   --- PJACOB (SOBROUTINE OF DIAGONALIZATION) ---
CCC
C     *      [HIJ][CJ] = \LAMBDA[NIJ][CJ]
C     *     INPUT
C     *     A = Hij
C     *     B = Nij
C     *
C     *     OUTPUT
C     *      \LAMBDA -> W
C     *      CJ -> AZ
C     *
      CALL PJACOB(A,B,W,AZ,MMAX,MDIM,1.0D0,  TEMP1,TEMP2,MDIM)
CCC
CCC   --- FOR LAPACK
CCC
C$$$      ITYPE=1
C$$$      LDA=MDIM
C$$$      LDB=MDIM
C$$$      LWORK=3D0*MMAX-1D0
C$$$      CALL DSYGV(ITYPE,'V','L',MMAX,A,LDA,B,LDB,W,WORK,
C$$$     &     LWORK,INFO)
C$$$
C$$$      IF(INFO.NE.0) STOP 9
c$$$      DO I=1,MMAX
c$$$      DO J=1,MMAX
c$$$         AZ(I,J)=A(I,J)
c$$$      END DO
c$$$      END DO
C
      RETURN
      END


C     ***************************************************************
CCC
CCC        ***** PJACOB *****
CCC
      SUBROUTINE PJACOB(H,FN,E,VEC,N,NDIM,FUGO ,  W,TEMP,KDIM)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION H(NDIM,NDIM),FN(NDIM,NDIM),E(NDIM),VEC(KDIM,KDIM)
      DIMENSION W(KDIM,KDIM),TEMP(KDIM,KDIM),SN(200),MAE(200)
      DIMENSION ANORM(200)
CCC ===========================================
      IF(N.GT.1) GO TO 70
      E(1)=H(1,1)/FN(1,1)
      VEC(1,1)=1.0/DSQRT(FN(1,1)*FUGO)
      RETURN
   70 CONTINUE
      DO 1 I=1,N
      DO 1 K=I,N
      TEMP(I,K)=FN(I,K)
    1 TEMP(K,I)=FN(I,K)
      CALL JACOBD(TEMP,KDIM,N,1.0D-14,VEC,ILL)
      IF(ILL.EQ.0) GO TO 50
      WRITE(6,100)ILL
  100 FORMAT(1H0,5X,5H*****,5X,4HILL=,I10)
      STOP
   50 CONTINUE
C
      DO 2 I=1,N
      ANORM(I)=TEMP(I,I)
    2 SN(I)=DSQRT(TEMP(I,I)*FUGO)
      DO 3 I=1,N
      DO 3 K=1,N
    3 W(I,K)=VEC(I,K)/SN(K)
      DO 4 I=1,N
      DO 4 K=I,N
      S=0.0
      DO 5 J=1,N
      DO 5 M=1,N
    5 S=S + W(J,I)*H(J,M)*W(M,K)
      TEMP(I,K)=S
    4 TEMP(K,I)=S
      CALL JACOBD(TEMP,KDIM,N,1.0D-14,VEC,ILL)
      IF(ILL.EQ.0) GO TO 60
      WRITE(6,100)ILL
      STOP
   60 CONTINUE
C
      DO 6 I=1,N
      E(I)=TEMP(I,I)*FUGO
    6 MAE(I)=I
      DO 10 I=1,N
      DO 10 K=1,N
      S=0.0
      DO 15 J=1,N
   15 S=S + W(I,J)*VEC(J,K)
   10 TEMP(I,K)=S
C
      DO 20 I=1,N
      DO 25 K=I,N
      IF( E(K).GE.E(I) ) GO TO 25
      EI=E(I)
      E(I)=E(K)
      E(K)=EI
      MAEI=MAE(I)
      MAE(I)=MAE(K)
      MAE(K)=MAEI
   25 CONTINUE
   20 CONTINUE
C
      DO 30 I=1,N
      W(I,1)=ANORM(I)
      DO 30 K=1,N
   30 VEC(I,K)=TEMP(I,MAE(K))
      RETURN
      END
C
C***   JACOBD   ***
C
      SUBROUTINE JACOBD(A,NDIM,N,RHO,S,ILL)
        IMPLICIT REAL*8(A-H,O-Z)
      DOUBLE PRECISION A,S,RHO,DEL,B,AN,ANORM,FNORM,THRESH,C,D,G,G1,
     1     ABSD,COV,SIV,APP,AQQ,U,V,ETA,ZERO
      DIMENSION A(NDIM,NDIM),S(NDIM,NDIM)
      DIMENSION DEL(200),B(200)
      ZERO=1.0D-30
      N1=N-1
      DO 9000 I=1,N1
      I1=I+1
      DO 9000 J=I1,N
      IF(DABS(A(I,J)-A(J,I)).GT.ZERO) GO TO 100
 9000 CONTINUE
      IF(N.GE.2.AND.NDIM.GE.N.AND.RHO.GT.0.0.AND.N.LE.150) GO TO 2
  100 ILL=30000
      RETURN
    2 IND=0
      ILL=0
      AN=N
C     COMPUTE INITIAL NORM(SQ. ROOT OF SUM OF SQUARES OF ALL
C     OFF-DIAGONAL ELEMENTS).
    3 ANORM=0.
      NM1=N-1
      DO 300 I=1,NM1
      IP1=I+1
      DO 300 J=IP1,N
  300 ANORM=ANORM+A(I,J)**2
      ANORM=DSQRT(2.0D0*ANORM)
C     COMPUTE FINAL NORM.
    4 FNORM=RHO*ANORM/AN
C     INITIALIZE MATRIX S(INITIALLY,  S  IS TO BE IDENTITY MATRIX).
      DO 400 I=1,N
      DO 400 J=1,N
  400 S(I,J)=0.
      DO 401 I=1,N
  401 S(I,I)=1.0D0
C     SET UP COUNTER TO BOUND NUMBER OF TRANSFORMATIONS PERFORMED.
C     THIS IS TO SERVE AS A SAFETY MEASURE TO AVOID INFINITE LOOPING.
      ISWEEP=0
C     SET INITIAL THRESHOLD.
    5 THRESH=ANORM/AN
C     SAVE DIAGONAL ELEMENTS. CORRECTION TO THE DIAGONAL ELEMENTS IS
C     DONE COLLECTIVELY AT THE END OF EACH ONE SWEEP.  DEL IS THE NAME
C     OF ARRAY TO STORE CORRECTIONS.
    6 DO 600 I=1,N
      B(I)=A(I,I)
  600 DEL(I)=0.
C     START OF SWEEP.
    7 IQ=2
    8 IP=1
    9 IF(DABS(A(IP,IQ))-THRESH) 10,10,19
   10 IF(IP-IQ+1)  11,12,12
   11 IP=IP+1
      GO TO 9
   12 IF(IQ-N) 13,14,14
   13 IQ=IQ+1
      GO TO 8
   14 IF(IND) 15,17,15
   15 DO 1500 I=1,N
 1500 A(I,I)=B(I)+DEL(I)
C     STOP IF 100 SWEEPS HAVE ALREADY BEEN MADE. HERE, THE NUMBER 100
C     FOR THE MAX. NUMBER OF SWEEPS WAS RATHER ARBITRARILY CHOSEN, BUT
C     IT WOULD NORMALLY BE LARGE ENOUGH.
      IF(ISWEEP-100) 16,1501,1501
 1501 ILL=1
      RETURN
   16 ISWEEP=ISWEEP+1
      IND=0
      GO TO 6
   17 IF(THRESH-FNORM) 18,18,1700
 1700 THRESH=THRESH/AN
      GO TO 7
   18 RETURN
C     DO TRANSFORMATION.
   19 IND=1
      C=2.0D0*A(IP,IQ)
      D=A(IQ,IQ)-A(IP,IP)
      G=C*C+D*D
      G1=DSQRT(G)
      ABSD=DABS(D)
      COV=DSQRT((1.0D0+ABSD/G1)/2.0D0)
      SIV=C/(2.0D0*COV*G1)
      IF(D) 1900,21,21
 1900 SIV=-SIV
   21 APP=A(IP,IP)
      AQQ=A(IQ,IQ)
   22 DO 2200 I=1,N
      U=A(IP,I)
      V=A(IQ,I)
      A(IP,I)=COV*U-SIV*V
      A(IQ,I)=COV*V+SIV*U
      U=S(I,IP)
      V=S(I,IQ)
      S(I,IP)=COV*U-SIV*V
 2200 S(I,IQ)=COV*V+SIV*U
   23 A(IP,IQ)=0.
      A(IQ,IP)=0.
C     COMPUTE CORRECTIONS TO A(IP,IP) AND A(IQ,IQ) ASSOCIATED WITH
C     THIS TRANSFORMATION.
   24 ETA=C*C/(2.0D0*(G1+ABSD))
      IF(D) 2400,25,25
 2400 ETA=-ETA
C     COMPUTE CUMULATIVE CORRECTIONS.
   25 DEL(IP)=DEL(IP)-ETA
      DEL(IQ)=DEL(IQ)+ETA
C     COMPUTE DIAGONAL ELEMENTS A(IP,IP) ADD A(IQ,IQ).
   26 A(IP,IP)=APP-ETA
      A(IQ,IQ)=AQQ+ETA
C     UPDATE IP-TH AND IQ-TH COLUMNS OF MATRIX A.
   27 DO 2700 I=1,N
      A(I,IP)=A(IP,I)
 2700 A(I,IQ)=A(IQ,I)
      GO TO 10
      END
