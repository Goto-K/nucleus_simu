CCC   *
CCC   *     VARIATIONAL METHOD WITH THE GAUSSIAN BASIS FUNCTIONS
CCC   *
      program nnp1_1
      implicit none
      integer   i,j,k,l,l1,l2,a1,a2
      integer   ndim,lmax,amax
      parameter (ndim=100,lmax=3,amax=5)
C      real*8,dimension(ndim,ndim,ndim,ndim)  ::  
C                N,H
      real*8,dimension(ndim,ndim)  :: 
                N1,H
      real*8,dimension(ndim)  :: 

      real*8    hc0,ac,ec,pi,hc,massp,massn,rem,ls,
     &          r0,rmax,dx,dy,r,dr,lx1,lx3,lxy1,lxy3,ly1,ly3,lyx1,lyx3
      integer   nr,nrmax
      parameter (hc0=197.326901d0,ac=931.494013d0,ec=137.03599976d0)
CCC   hc0[MeV*fm],ac[Mev],ec=hc/e^2

CCC   α parameter
      real    ::  lx(amax)=(/0,0,2,0,2/)
      real    ::  Ly(amax)=(/0,0,0,2,2/)
      real    ::  nxmax(amax)=(/15,15,15,15,15/)
      real    ::  Nymax(amax)=(/15,15,15,15,15/)
      real*8  ::  x0(amax)=(/0.05,0.05,0.1,0.1,0.1/)
      real*8  ::  xmax(amax)=(/15.0,15.0,15.0,15.0,15.0/)
      real*8  ::  y0(amax)=(/0.3,0.3,0.3,0.3,0.3/)
      real*8  ::  ymax(amax)=(/9.0,9.0,9.0,9.0,9.0/)

CCC   Potential Models of Nuclear Forces at Small Distances,R.Tamagaki,1968
      real*8  ::  CVC(3) = (/-5d0,-230d0,2000d0/)
      real*8  ::  ic(3) = (/2.5,0.942,0.447/)
      real*8  ::  CVT(3) = (/-7.5,-67.5,67.5/)
      real*8  ::  it(3) = (/2.5,1.2,0.477/)
CCC   CVC,CVT[MeV],ic,it[fm]


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC   Three-Body

CCC   1p,2n,3n/r1=[2n~3n]
CCC   three-body hamiltonian
CCC   H=T+V(r1)+V(r2)+V(r3)
CCC   wave function  α:channel
CCC   Ψ=Φ(r1,R1)+Φ(r2,R2)+Φ(r3,R3)
CCC   Φ(ri,Ri)=Σ_α u_α(ri,Ri) y_α(jk,i)
CCC   y_α(jk,i)=[[Y_l(ri)(*)Y_L(Ri)]_Λ (*) [x_s(jk)(*)x_1/2(i)]_Σ]_{JK}
CCC             (*)[η_t(jk)(*)η_1/2(i)]_TTz
CCC   u_α(ri,Ri)=Σ_{n,N} A_{αnN} φ_nl(ri) Ψ_NL(Ri)

CCC   Ψ=Σ_{α,n,N} A_{α,n,N} [Φ_{α,n,N}(r1,R1)+Φ_{α,n,N}(r2,R2)+Φ_{α,n,N}(r3,R3)]
CCC   Φ_{α,n,N}(ri,Ri)=φ_nl(ri) Ψ_NL(Ri) y_α(jk,i)
CCC   < Φ_{α’,n,N}(r1,R1) | H-E | Ψ >
CCC   Σ_{α,n,N} (H_{α’,α}-EN_{α’,α}) A_α = 0
CCC    H_{α’,α}=< Φ_{α’,n,N}(1) | H | Φ_{α,n,N}(1)+Φ_{α,n,N}(2)+Φ_{α,n,N}(3) >
CCC           (i)=(ri,Ri)
CCC    H_{α’,α}=H_{α,α’} / ↑ not orthogonal to one another unless the quantum numbers
CCC                          (Λ,Σ) are different between them

CCC   y_α(jk,i)=[[Y_l(ri)(*)Y_L(Ri)]_Λ (*) [x_s(jk)(*)x_1/2(i)]_Σ]_{JK}
CCC             (*)[η_t(jk)(*)η_1/2(i)]_TTz
CCC   "pauli principle" l+s+t=odd , t=1 {nn,np,pp},t=0 {np}

CCC   α l L Λ s Σ   t (particle)
CCC   1 0 0 0 0 1/2 1 (nn,np)
CCC   2 0 0 0 1 1/2 0 (np)
CCC   3 2 0 2 1 3/2 0 (np)
CCC   4 0 2 2 1 3/2 0 (np)
CCC   5 2 2 0 1 1/2 0 (np)

CCC   AV14 potential
CCC   v14 = Σ [v_π + v_I + v_s]*P
CCC   P=1,τj*τk,σj*σk,(σj*σk)(τj*τk),Sjk,Sjk(τj*τk),L*S,(L*S)(τj*τk),O(L^2)

CCC   (τj*τk) |y_α(jk,i)> = 1/2(t(t+1)-3/2)
CCC   (σj*σk) |y_α(jk,i)> = 1/2(s(s+1)-3/2)
CCC   Sjk
CCC   L(*)S
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC   Harmonic Oscillator

CCC   r1,R1
CCC   r2=-r1/2-R1,R2=3r1/4-R1/2
CCC   r3=-r1/2+R1,R2=-3r1/4-R1/2
CCC   <1|V1|1>,<1|V1|2>,<1|V1|3>,<1|V2|1>,<1|V2|2>,<1|V2|3>,
CCC   <1|V3|1>,<1|V3|2>,<1|V3|3>
CCC   <1|V1|1>=r1^(la+lb+4)*exp[-(bi+bj)*r1^2]*R1^(La+Lb+2)*exp[-(Bi+Bj)*R1^2]
CCC
CCC   <1|V1|2>=(1/4)^qb*(9/16)^Qb*r1^(la+2qb+2Qb+4)*exp[-(bi+(1/4)bj+(9/16)Bj)*r1^2]
CCC           *(1/4)^Qb*R1^(La+2qb+2Qb+2)*exp[-(bj+Bi+(1/4)Bj)*R1^2]
CCC           =<1|V1|3>=<1|V2|2>=<1|V3|3>
CCC
CCC   <1|V2|1>=(1/4)*r1^(la+lb+4)*exp[-(bi+bj)*r1^2]
CCC           *R1^(La+Lb+4)*exp[-(Bi+Bj)*R1^2]
CCC           =<1|V3|1>
CCC
CCC   <1|V2|3>=(1/4)*(1/4)^qb*(9/16)^Qb*r1^(la+2qb+2Qb+4)
CCC                                    *exp[-(bi+(1/4)bj+(9/16)Bj)*r1^2]
CCC           *(1/4)^Qb*R1^(La+2lb+2Qb+4)*exp[-(bj+Bi+(1/4)Bj)*R1^2]
CCC           =<1|V3|2>
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC



      do a1=1,amax
      do a2=1,amax

      pi = acos(-1.0d0)
C      hc = hc0*10  ![eV*Å]
CCC   proton,neutron,redeuced mass      
C      massp = 1.00727646688d0 * ac 
C      massn = 1.00866491578d0 * ac
C      rem = massp*massn / (massp+massn)
      rem = ac**2 / (2*ac)
CCC   angular momentum
C      ls = 0
CCC   range parameter
      dx1 = (xmax(a1)/x0(a1))**(1.0d0/(nxmax(a1)-1.0d0))
      dx2 = (xmax(a2)/x0(a2))**(1.0d0/(nxmax(a2)-1.0d0))
      dy1 = (ymax(a1)/y0(a1))**(1.0d0/(Nymax(a1)-1.0d0))
      dy2 = (ymax(a2)/y0(a2))**(1.0d0/(Nymax(a2)-1.0d0))

CCC
CCC   (la+lb+1)!!,(la+lb+3)!!,(La+Lb+1)!!,(La+Lb+3)!!,
CCC   (la+lb+Lb+1)!!,(la+lb+Lb+3)!!,(La+lb+Lb+1)!!,(La+lb+Lb+3)!!
CCC

      lx1 = 1d0
      do li=1,lx(a1)+lx(a2)+1,2
      lx1 = li * lx1
      end do
      lx3 = lx1 * (lx(a1)+lx(a2)+3)

      ly1 = 1d0
      do li=1,ly(a1)+ly(a2)+1,2
      ly1 = li * ly1
      end do
      ly3 = ly1 * (ly(a1)+ly(a2)+3)

      lxy1 = 1d0
      do li=1,lx(a1)+lx(a2)+ly(a2)+1,2
      lxy1 = li * lxy1
      end do
      lxy3 = lxy1 * (lx(a1)+lx(a2)+ly(a2)+3)

      lyx1 = 1d0
      do li=1,ly(a1)+lx(a2)+ly(a2)+1,2
      lyx1 = li * lxy1
      end do
      lyx3 = lyx1 * (ly(a1)+lx(a2)+ly(a2)+3)


CCC
CCC  calculation of matrix elements
CCC

      do i=1,nxmax(a1)
      do k=1,nymax(a1)

      bx1(i) = 1.0d0 / (r0(a1) * dx1**(i-1))**2
      bx2(i) = 1.0d0 / (r0(a2) * dx2**(i-1))**2
      by1(i) = 1.0d0 / (r0(a1) * dy1**(k-1))**2
      by2(i) = 1.0d0 / (r0(a2) * dy2**(k-1))**2

CCC   norm matrix

      do j=1,i
      do l=1,k

      N(i,j,k,l)=r1^(la+lb+4)*exp[-(bi+bj)*r1^2]*R1^(La+Lb+2)*exp[-(Bi+Bj)*R1^2]
      N(i,j,k,l)=N(j,i,l,k)

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

