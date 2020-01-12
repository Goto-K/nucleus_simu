CCC   *
CCC   *     VARIATIONAL METHOD WITH THE GAUSSIAN BASIS FUNCTIONS
CCC   *
      program deuteron1_2
      implicit none
      integer   i,j,k,l,l1,l2
      integer   ndim,nmax,lmax
      parameter (ndim=1000,lmax=3)
      real*8,dimension(ndim,ndim,lmax,lmax)  ::  
     &          N,N0,T,V,VC,VT,VLS,H
      real*8,dimension(ndim,ndim)  :: 
     &          n_b1,n_b2,evec1,evec2,vec0,vec1,vec2,vec3,vec4,
     &          h_b,hvec_b,vecs,vs1,vs2,vd1,vd2
      real*8,dimension(ndim)  :: 
     &          CN0,CN2,b,eig1,ene1,ene2,es1,es2,ed1,ed2,wf0,wf2
      real*8    hc0,ac,ec,pi,hc,massp,massn,rem,ls,
     &          r0,rmax,dx,r,dr
      integer   nr,nrmax
      parameter (hc0=197.326901d0,ac=931.494013d0,ec=137.03599976d0)
CCC   hc0[MeV*fm],ac[Mev],ec=hc/e^2
CCC   Potential Models of Nuclear Forces at Small Distances,R.Tamagaki,1968
      real*8  ::  CVC(3) = (/-5d0,-230d0,2000d0/)
      real*8  ::  ic(3) = (/2.5,0.942,0.447/)
      real*8  ::  CVT(3) = (/-7.5,-67.5,67.5/)
      real*8  ::  it(3) = (/2.5,1.2,0.477/)
CCC   CVC,CVT[MeV],ic,it[fm]


      pi = acos(-1.0d0)
C      hc = hc0*10  ![eV*Å]
CCC   proton,neutron,redeuced mass      
      massp = 1.00727646688d0 * ac 
      massn = 1.00866491578d0 * ac
      rem = massp*massn / (massp+massn)
CCC   angular momentum
C      ls = 0
CCC   range parameter
      r0 = 0.1d0
      rmax = 20d0
      nmax = 40
      dx = (rmax/r0)**(1.0d0/(nmax-1.0d0))

CCC
CCC  calculation of matrix elements
CCC

      do i=1,nmax    

      b(i) = 1.0d0 / (r0 * dx**(i-1))**2

      CN0(i)=sqrt(4d0*sqrt(((2d0*b(i))**3)/pi))
      CN2(i)=sqrt(16d0/15d0*sqrt(((2d0*b(i))**7)/pi))
      
CCC   |Fai> = U(r)/r|S> + W(r)/r|D>
CCC   X(r) = Σ_L Σ_n C(L)_n * r^(L+1) * exp(-b_n*r^2) |L>
CCC   |fai> = r^(L+1) * exp(-b_n*r^2) |L>
CCC   N = <fai|fai>       

CCC   norm matrix

      do j=1,i

      do l1=1,lmax
      do l2=1,lmax

      if(l1==1.and.l2==1)then
      N(i,j,l1,l2)
     & = CN0(i)*CN0(j)/4d0*sqrt((pi)/((b(i)+b(j))**3))
      n_b1(i,j)=N(i,j,1,1)
      n_b1(j,i)=n_b1(i,j)
      n_b2(i,j)=n_b1(i,j)
      n_b2(j,i)=n_b2(i,j)
      else if(l1==3.and.l2==3)then
      N(i,j,l1,l2)
     & = CN2(i)*CN2(j)*15d0/16d0*sqrt((pi)/((b(i)+b(j))**7))
      n_b1(nmax+i,nmax+j)=N(i,j,3,3)
      n_b1(nmax+j,nmax+i)=n_b1(nmax+i,nmax+j)
      n_b2(nmax+i,nmax+j)=n_b1(nmax+i,nmax+j)
      n_b2(nmax+j,nmax+i)=n_b2(nmax+i,nmax+j)
      else
      N(i,j,l1,l2)=0
      end if

CCC   kinetic energy matrix
      
      if(l1==1.and.l2==1)then
      T(i,j,l1,l2)
     & = CN0(i)*CN0(j)*(hc0**2)/(2*rem)
     & * 3d0/2d0*b(i)*b(j) * sqrt((pi)/((b(i)+b(j))**5))
      else if(l1==3.and.l2==3)then
      T(i,j,l1,l2)
     & = CN2(i)*CN2(j)*(hc0**2)/(2*rem)
     & * 105d0/8d0*b(i)*b(j) * sqrt((pi)/((b(i)+b(j))**9))
      else
      T(i,j,l1,l2)=0
      end if

      T(j,i,l1,l2)=T(i,j,l1,l2)

CCC   potential energy matrix
CCC
CCC   V = VC + VT*S12 + VLS*l(*)s
CCC
CCC   VC(i,j,l1,l2)
      
      do k=1,3
      
      if(l1==1.and.l2==1)then
      VC(i,j,l1,l2)=VC(i,j,l1,l2)
     & + CN0(i)*CN0(j)*CVC(k)
     & * (sqrt(pi)/(sqrt(b(i)+b(j)+1d0/ic(k)**2))**3)/4d0
      else if(l1==3.and.l2==3)then
      VC(i,j,l1,l2)=VC(i,j,l1,l2)
     & + CN2(i)*CN2(j)*CVC(k)
     & * (sqrt(pi)/(sqrt(b(i)+b(j)+1d0/ic(k)**2))**7)*15d0/16d0
      else if(l1==1.and.l2==3)then
      VC(i,j,l1,l2)=VC(i,j,l1,l2)
     & + CN0(i)*CN2(j)*CVC(k)
     & * (sqrt(pi)/(sqrt(b(i)+b(j)+1d0/ic(k)**2))**5)*3d0/8d0
      else if(l1==3.and.l2==1)then
      VC(i,j,l1,l2)=VC(i,j,l1,l2)
     & + CN2(i)*CN0(j)*CVC(k)
     & * (sqrt(pi)/(sqrt(b(i)+b(j)+1d0/ic(k)**2))**5)*3d0/8d0
      else
      VC(i,j,l1,l2)=0
      end if

      end do

      VC(j,i,l1,l2)=VC(i,j,l1,l2)

CCC   VT(i,j,l1,l2)
      
      do k=1,3
      
      if(l1==1.and.l2==1)then
      VT(i,j,l1,l2)=VT(i,j,l1,l2)
     & + CN0(i)*CN0(j)*CVT(k)
     & * (sqrt(pi)/(sqrt(b(i)+b(j)+1d0/it(k)**2))**3)/4d0
      else if(l1==3.and.l2==3)then
      VT(i,j,l1,l2)=VT(i,j,l1,l2)
     & + CN2(i)*CN2(j)*CVT(k)
     & * (sqrt(pi)/(sqrt(b(i)+b(j)+1d0/it(k)**2))**7)*15d0/16d0
      else if(l1==1.and.l2==3)then
      VT(i,j,l1,l2)=VT(i,j,l1,l2)
     & + CN0(i)*CN2(j)*CVT(k)
     & * (sqrt(pi)/(sqrt(b(i)+b(j)+1d0/it(k)**2))**5)*3d0/8d0
      else if(l1==3.and.l2==1)then
      VT(i,j,l1,l2)=VT(i,j,l1,l2)
     & + CN2(i)*CN0(j)*CVT(k)
     & * (sqrt(pi)/(sqrt(b(i)+b(j)+1d0/it(k)**2))**5)*3d0/8d0
      else
      VT(i,j,l1,l2)=0
      end if

      end do

      VT(j,i,l1,l2)=VT(i,j,l1,l2)

CCC   VLS(i,j,l1,l2)

      VLS(i,j,l1,l2)=0
      VLS(j,i,l1,l2)=0

      end do
      end do
      end do
      end do

CCC

      do i=1,nmax 
      do j=1,i
      do l1=1,lmax
      do l2=1,lmax

CCC
CCC   V(i,j,l1,l2)


      if(l1==1.and.l2==1)then
      V(i,j,l1,l2)
     & = VC(i,j,l1,l2)
      else if(l1==3.and.l2==3)then
      V(i,j,l1,l2)
     & = VC(i,j,l1,l2)-2*VT(i,j,l1,l2)-3*VLS(i,j,l1,l2)
      else if(l1==1.and.l2==3)then
      V(i,j,l1,l2)
     & = sqrt(8d0)*VT(i,j,l1,l2)
      else if(l1==3.and.l2==1)then
      V(i,j,l1,l2)
     & = sqrt(8d0)*VT(i,j,l1,l2)
      else
      V(i,j,l1,l2)=0
      end if

      V(j,i,l1,l2)=V(i,j,l1,l2)
      
CCC   hamiltonian matrix

      H(i,j,l1,l2) = T(i,j,l1,l2) + V(i,j,l1,l2)
      H(j,i,l1,l2) = H(i,j,l1,l2)

      h_b(i,j)=H(i,j,1,1)
      h_b(j,i)=h_b(i,j)
      h_b(i,nmax+j)=H(i,j,1,3)
      h_b(nmax+j,i)=h_b(i,nmax+j)
      h_b(nmax+i,j)=H(i,j,3,1)
      h_b(j,nmax+i)=h_b(nmax+i,j)
      h_b(nmax+i,nmax+j)=H(i,j,3,3)
      h_b(nmax+j,nmax+i)=h_b(nmax+i,nmax+j)

      end do
      end do
      end do
      end do

CCC
CCC   HQR(NORM)
CCC

      call hqr(ndim,2*nmax,n_b1,evec1,eig1)

CCC   |f> = Σ_a evec/√eig |fai>

      do i=1,2*nmax
       do j=1,2*nmax
        vecs(i,j)=0
        do k=1,2*nmax
         vecs(i,j) = vecs(i,j) + evec1(i,k) * evec1(j,k)
        end do
       end do
      end do

      do i=1,2*nmax
C      write(*,'(F10.5,9(1X,F10.5))') (vecs(j,i),j=1,2*nmax)
      end do
      
      write(*,*)
      
CCC   standardization

      do i=1,2*nmax
       do j=1,2*nmax
        vec1(i,j) = evec1(i,j) / sqrt(eig1(j))
       end do
      end do

      do i=1,2*nmax
       do j=1,2*nmax
        vec2(i,j) = 0d0
        do k=1,2*nmax
         do l=1,2*nmax
         vec2(i,j) = vec2(i,j) + vec1(k,i)*n_b2(k,l)*vec1(l,j)
         end do
        end do
       end do
      end do

      do i=1,2*nmax
C      write(*,'(F10.5,9(1X,F10.5))') (vec2(i,j),j=1,2*nmax)
      end do

      write(*,*)

CCC   new hamiltonian

      do i=1,2*nmax
       do j=1,2*nmax
        hvec_b(i,j) = 0d0
        do k=1,2*nmax
         do l=1,2*nmax
         hvec_b(i,j) = hvec_b(i,j) + vec1(k,i)*h_b(k,l)*vec1(l,j)
         end do
        end do
       end do
      end do

      do i=1,2*nmax
C      write(*,'(F10.5,9(1X,F10.5))') (hvec_b(i,j),j=1,2*nmax)
      end do
C      write(*,*)

CCC
CCC   HQR(HAMILTONIAN)
CCC

      call hqr(ndim,2*nmax,hvec_b,evec2,ene1)

      do i=1,2*nmax
C      write(*,*) i,ene1(i)
      end do

CCC   new base vactor

      do i=1,2*nmax
       do k=1,2*nmax
        vec3(i,k)=0d0
        do j=1,2*nmax
         vec3(i,k) = vec3(i,k) + evec2(j,k) * vec1(i,j)
        end do
       end do
      end do

CCC   sort

      do k=1,2*nmax
       i=0
       do j=1,2*nmax
        if(ene1(k).ge.ene1(j)) i=i+1
       end do
       ene2(i)=ene1(k)
       do j=1,2*nmax
        vec4(j,i)=vec3(j,k)
       end do
      end do

      do i=1,10
C      write(*,*) i,ene2(i)
      end do

CCC   ground state

      dr = 0.005d0
      rmax = 10d0
      nrmax = rmax/dr + 1
C      k=1

      do nr=1,nrmax
      r=dr*(nr-1)
      wf0=0d0
      wf2=0d0
      do k=1,10
      do i=1,nmax
       wf0(k)=wf0(k)+vec4(i,k)*CN0(i)*exp(-b(i)*r**2)
       wf2(k)=wf2(k)+vec4(nmax+i,k)*CN2(i)*r**2*exp(-b(i)*r**2)
      end do
      end do

      write(*,*) r,wf0(1),wf2(1)
      end do


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

