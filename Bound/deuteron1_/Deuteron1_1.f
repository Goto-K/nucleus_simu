CCC   *
CCC   *     VARIATIONAL METHOD WITH THE GAUSSIAN BASIS FUNCTIONS
CCC   *
      program classical_mechanics
      implicit none
      integer   i,j,k,l,l1,l2
      real*8    ndim,nmax
      parameter (ndim=100,lmax=3)
      real*8    N(ndim,ndim,lmax,lmax),N0(ndim,ndim,lmax,lmax)
     &          T(ndim,ndim,lmax,lmax),V(ndim,ndim,lmax,lmax)
     &          VC(ndim,ndim,lmax,lmax),VT(ndim,ndim,lmax,lmax)
     &          VLS(ndim,ndim,lmax,lmax),H(ndim,ndim,lmax,lmax)
     &          b(ndim)
      real*8    hc0,ac,ec,pi,hc,massp,massn,rem,ls,r0,rmax,dx
      parameter (hc0=197.326901d0,ac=931.494013d0,ec=137.03599976d0)
      real*8  ::  CVC(3) = (/-5.,-230.,2000./)
      real*8  ::  ic(3) = (/2.5,0.942,0.477/)
      real*8  ::  CVT(3) = (/-7.5,-67.5,67.5/)
      real*8  ::  it(3) = (/2.5,1.2,0.477/)


      pi = cos(-1.0d0)
      hc = hc0 * 10
CCC   proton,neutron,redeuced mass      
      massp = 1.0086649d0 * ac * 1d+6
      massn = 1.00137841898 * massp
      rem = massp*massn / (massp+massn)
CCC   angular momentum
      ls = 0
CCC   range parameter
      r0 = 0.1d0
      rmax = 20d0
      nmax = 20
      dx = (rmax/r0)**(1.0d0/(nmax-1.0d0))

CCC
CCC  calculation of matrix elements
CCC

      do i=1,nmax    

      b(i) = r0 * dx**(i-1)
      
CCC   norm matrix

      do j=1,i

      do l1=1,lmax
      do l2=1,lmax

      if(l1==1.and.l2==1)
      N(i,j,l1,l2)
     & = sqrt(pi)/sqrt(b(i)+b(j))/4/(b(i)+b(j))
      else if(l1==3.and.l2==3)
      N(i,j,l1,l2)
     & = sqrt(pi)/sqrt(b(i)+b(j))*15/16/(b(i)+b(j))**3
      else
      N(i,j,l1,l2)=0
      end if

      N(j,i,l1,l2)=N(i,j,l1,l2)
      N0(i,j,l1,l2)=N(i,j,l1,l2)
      N0(j,i,l1,l2)=N(j,i,l1,l2)

CCC   kinetic energy matrix
      
      if(l1==1.and.l2==1)
      T(i,j,l1,l2)
     & = -hc**2/rem * sqrt(pi)/sqrt(b(i)+b(j))/4/(b(i)+b(j))
     & * (-6)*b(i)*b(j)/(b(i)+b(j))
      else if(l1==3.and.l2==3)
      T(i,j,l1,l2)
     & = sqrt(pi)/sqrt(b(i)+b(j))*15/16/(b(i)+b(j))
     & * (-14)*b(i)*b(j)/(b(i)+b(j))
      else
      T(i,j,l1,l2)=0
      end if

      T(j,i,l1,l2)=T(i,j,l1,l2)

CCC   potential energy matrix
CCC   V = VC + VT*S12 + VLS*l(*)s
CCC
CCC   VC(i,j,l1,l2)
      
      do k=1,3
      
      if(l1==1.and.l2==1)
      VC(i,j,l1,l2)+VC(i,j,l1,l2)
     & = CVC(k) * sqrt(pi)/sqrt(b(i)+b(j)+ic(k))/4/(b(i)+b(j)+ic(k))
      else if(l1==3.and.l2==3)
      VC(i,j,l1,l2)+VC(i,j,l1,l2)
     & = CVC(k) * sqrt(pi)/sqrt(b(i)+b(j)+ic(k))*15/16/(b(i)+b(j)+ic(k))
      else if(l1==1.and.l2==3)
      VC(i,j,l1,l2)+VC(i,j,l1,l2)
     & = CVC(k) * sqrt(pi)/sqrt(b(i)+b(j)+ic(k))*3/8/(b(i)+b(j)+ic(k))
      else if(l1==3.and.l2==1)
      VC(i,j,l1,l2)+VC(i,j,l1,l2)
     & = CVC(k) * sqrt(pi)/sqrt(b(i)+b(j)+ic(k))*3/8/(b(i)+b(j)+ic(k))
      else
      VC(i,j,l1,l2)=0
      end if

      end do

      VC(j,i,l1,l2)=VC(i,j,l1,l2)

CCC   VT(i,j,l1,l2)
      
      do k=1,3
      
      if(l1==1.and.l2==1)
      VT(i,j,l1,l2)+VT(i,j,l1,l2)
     & = CVT(k) * sqrt(pi)/sqrt(b(i)+b(j)+it(k))/4/(b(i)+b(j)+it(k))
      else if(l1==3.and.l2==3)
      VT(i,j,l1,l2)+VT(i,j,l1,l2)
     & = CVT(k) * sqrt(pi)/sqrt(b(i)+b(j)+it(k))*15/16/(b(i)+b(j)+it(k))
      else if(l1==1.and.l2==3)
      VT(i,j,l1,l2)+VT(i,j,l1,l2)
     & = CVT(k) * sqrt(pi)/sqrt(b(i)+b(j)+it(k))*3/8/(b(i)+b(j)+it(k))
      else if(l1==3.and.l2==1)
      VT(i,j,l1,l2)+VT(i,j,l1,l2)
     & = CVT(k) * sqrt(pi)/sqrt(b(i)+b(j)+it(k))*3/8/(b(i)+b(j)+it(k))
      else
      VT(i,j,l1,l2)=0
      end if

      end do

      VT(j,i,l1,l2)=VT(i,j,l1,l2)

CCC   VLS(i,j,l1,l2)

      VLS(i,j,l1,l2)=0
      VLS(j,i,l1,l2)=0

CCC
CCC   V(i,j,l1,l2)

      if(l1==1.and.l2==1)
      V(i,j,l1,l2)
     & = VC(i,j,l1,l2)
      else if(l1==3.and.l2==3)
      V(i,j,l1,l2)
     & = VC(i,j,l1,l2)-2*VT(i,j,l1,l2)-3*VLS(i,j,l1,l2)
      else if(l1==1.and.l2==3)
      V(i,j,l1,l2)
     & = sqrt(8)*VT(i,j,l1,l2)
      else if(l1==3.and.l2==1)
      V(i,j,l1,l2)
     & = sqrt(8)*VT(i,j,l1,l2)
      else
      VT(i,j,l1,l2)=0
      end if

      V(j,i,l1,l2)=V(i,j,l1,l2
      
CCC   hamiltonian matrix

      H(i,j,l1,l2) = T(i,j,l1,l2) + V(i,j,l1,l2)
      H(j,i,l1,l2) = H(i,j,l1,l2)

      end do
      end do
      end do
      end do

CCC
CCC   HQR
CCC

      call HQR(dim,nmax,N,VEC1,ENE1)

      stop
      end

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

