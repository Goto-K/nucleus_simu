CCC   *
CCC   *     VARIATIONAL METHOD WITH THE GAUSSIAN BASIS FUNCTIONS
CCC   *
      program positronium
      implicit none
      integer   i,j,k,l,li
      integer   ndim,nmax,ls
      parameter (ndim=100)
      real*8,dimension(ndim,ndim)  ::  
     &          N,N0,T,V,VC,VT,VLS,H
      real*8,dimension(ndim,ndim)  :: 
     &          evec1,evec2,vec0,vec1,vec2,vec3,vec4,
     &          hvec,vecs
      real*8,dimension(ndim)  :: 
     &          CN,b,eig1,ene1,ene2,es1,es2,ed1,ed2
      real*8    hc0,ac,ec,pi,hc,masse,massp,rem,
     &          r0,rmax,dx,r,dr,dfc1,dfc3,sf0,wf
      integer   nr,nrmax
      parameter (hc0=197.326901d0,ac=931.494013d0,ec=137.03599976d0)
CCC   hc0[MeV*fm],ac[Mev],ec=hc/e^2


      pi = acos(-1.0d0)
      hc = hc0*10  ![eV*Å]
CCC   proton,neutron,redeuced mass      
      masse = 0.511006d0 * 1d+6
      massp = 0.511006d0 * 1d+6 ![eV*Å]
      rem = masse*massp / (masse+massp)
C      write(*,*) 'rem=',rem
CCC   angular momentum
      ls = 0
CCC   range parameter[Å]
      r0 = 0.1d0
      rmax = 20d0
      nmax = 50
      dx = (rmax/r0)**(1.0d0/(nmax-1.0d0))

CCC

C      do ls=0,1

CCC
CCC   (2*LS+1)!!,(2*LS+3)!!,LS!
CCC

      dfc1 = 1d0
      dfc3 = 1d0
      do li=1,2*ls+1,2
      dfc1=li*dfc1
      end do
      dfc3 = dfc1*(2*ls+3)

      sf0 = 1d0
      do li=1,ls
      sf0 = li*sf0
      end do

CCC
CCC  calculation of matrix elements
CCC

      do i=1,nmax    

      b(i) = 1.0d0 / (r0 * dx**(i-1))**2

      CN(i)=sqrt(2d0**(ls+2)/dfc1
     &           * sqrt(((2d0*b(i))**(2*ls+3))/pi))
      
CCC   norm matrix

      do j=1,i

      N(i,j)
     & = CN(i)*CN(j)/2d0**(ls+2)
     & * sqrt((pi)/((b(i)+b(j))**(2*ls+3)))

      N(j,i)=N(i,j)
      N0(i,j)=N(i,j)
      N0(j,i)=N0(i,j)

CCC   kinetic energy matrix
      
      T(i,j)
     & = CN(i)*CN(j)*(hc**2)/rem*dfc3/2d0**(ls+2)
     & * b(i)*b(j) * sqrt((pi)/((b(i)+b(j))**(2*ls+5)))

      T(j,i)=T(i,j)

CCC   potential energy matrix
CCC
CCC   coulomb
CCC
CCC   V(i,j)

      V(i,j)
     & = -CN(i)*CN(j)*hc/ec*sf0
     & / (2d0*(b(i)+b(j))**(ls+1))

      V(j,i)=V(i,j)
      
CCC   hamiltonian matrix

      H(i,j) = T(i,j) + V(i,j)
      H(j,i) = H(i,j)

      end do
      end do

CCC
CCC   HQR(NORM)
CCC

      call hqr(ndim,nmax,N,evec1,eig1)

CCC   |f> = Σ_a evec/√eig |fai>

      do i=1,nmax
       do j=1,nmax
        vecs(i,j)=0
        do k=1,nmax
         vecs(i,j) = vecs(i,j) + evec1(i,k) * evec1(j,k)
        end do
       end do
      end do

      do i=1,nmax
C      write(*,'(F10.5,9(1X,F10.5))') (vecs(j,i),j=1,nmax)
      end do
      
C      write(*,*)
      
CCC   standardization

      do i=1,nmax
       do j=1,nmax
        vec1(i,j) = evec1(i,j) / sqrt(eig1(j))
       end do
      end do

      do i=1,nmax
       do j=1,nmax
        vec2(i,j) = 0d0
        do k=1,nmax
         do l=1,nmax
         vec2(i,j) = vec2(i,j) + vec1(k,i)*N0(k,l)*vec1(l,j)
         end do
        end do
       end do
      end do

      do i=1,nmax
C      write(*,'(F10.5,9(1X,F10.5))') (vec2(i,j),j=1,nmax)
      end do

C      write(*,*)

CCC   new hamiltonian

      do i=1,nmax
       do j=1,nmax
        hvec(i,j) = 0d0
        do k=1,nmax
         do l=1,nmax
         hvec(i,j) = hvec(i,j) + vec1(k,i)*H(k,l)*vec1(l,j)
         end do
        end do
       end do
      end do

      do i=1,nmax
C      write(*,'(F10.5,9(1X,F10.5))') (hvec(i,j),j=1,nmax)
      end do
C      write(*,*)

CCC
CCC   HQR(HAMILTONIAN)
CCC

      call hqr(ndim,nmax,hvec,evec2,ene1)

      do i=1,nmax
C      write(*,*) i,ene1(i)
      end do

CCC   new base vactor

      do i=1,nmax
       do k=1,nmax
        vec3(i,k)=0d0
        do j=1,nmax
         vec3(i,k) = vec3(i,k) + evec2(j,k) * vec1(i,j)
        end do
       end do
      end do

CCC   sort

      do k=1,nmax
       i=0
       do j=1,nmax
        if(ene1(k).ge.ene1(j)) i=i+1
       end do
       ene2(i)=ene1(k)
       do j=1,nmax
        vec4(j,i)=vec3(j,k)
       end do
      end do

      write(*,*) 'L=',ls
      do i=1,10
C      write(*,*) i,ene2(i)
      end do
C      write(*,*)

CCC   ground state

      dr = 0.005d0
      rmax = 10d0
      nrmax = rmax/dr + 1
      k=1

CCC   u(r)/r = wf
CCC   u(r) = Σ Cn*r^(2l+1)*exp[-bn*r^2]

      do nr=1,nrmax
      r=dr*(nr-1)
      wf=0d0
      do i=1,nmax
       wf = wf + vec4(i,k)*CN(i)*r**(ls)*exp(-b(i)*r**2)
      end do

      write(*,*) r,wf
      end do

CCC
C      end do   ! ls
CCC

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

