CCC   *
CCC   *     VARIATIONAL METHOD WITH THE GAUSSIAN BASIS FUNCTIONS
CCC   *
      program positronium
      implicit none
      integer   i,j,k,l,li
      integer   ndim,nmax,ls
      parameter (ndim=100)
      real*8,dimension(ndim,ndim,2)  ::  
     &          VSS,H
      real*8,dimension(ndim,ndim)  ::  
     &          N,N0,T,V,VC
      real*8,dimension(ndim,ndim)  :: 
     &          evec1,evec2,evec3,
     &          vec0,vec1,vec2,vec3,vec4,vec5,vec6,
     &          hvec0,hvec1,vecs
      real*8,dimension(ndim)  :: 
     &          CN,b,eig1,ene1,ene2,ene3,ene4
      real*8    hc0,ac,ec,pi,hc,ge,mu_b,masse,massp,rem,sigma,
     &          r0,rmax,dx,r,dr,dfc1,dfc3,sf0,wf0,wf1
      integer   nr,nrmax
      parameter (hc0=197.326901d0,ac=931.494013d0,ec=137.03599976d0)
      parameter (ge=2.002319d0,sigma=1d-12)
CCC   hc0[MeV*fm],ac[Mev],ec=hc/e^2


      pi = acos(-1.0d0)
      hc = hc0*10  ![eV*Å]
CCC   proton,neutron,redeuced mass      
      masse = 0.511006d0 * 1d+6
      massp = 1.00727646688d0*ac * 1d+6 ![eV*Å]
      rem = masse*massp / (masse+massp)
C      write(*,*) 'rem=',rem
CCC   bohr magneton
      mu_b = sqrt(hc/ec)*hc / (2*masse)
CCC   angular momentum
      ls = 0
CCC   range parameter[Å]
      r0 = 0.1d0
      rmax = 20d0
      nmax = 40
      dx = (rmax/r0)**(1.0d0/(nmax-1.0d0))

CCC

C      do ls=0,2

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
CCC   VC(i,j)

      VC(i,j)
     & = -CN(i)*CN(j)*hc/ec*sf0
     & / (2d0*(b(i)+b(j))**(ls+1))

      VC(j,i)=VC(i,j)
     
CCC   spin-spin coupling
CCC  
CCC   VSS=-2/3*(μe(*)μp)
CCC   μe{->}=-ge*mu_b*S1{->}/h   psは逆符号
CCC   μe(*)μp=-(ge*mu_b/h)^2*S1(*)S2
CCC   S12 (↑↑/↑↓/↓↑/↓↓)
CCC   S=1 |↑↑> |↓↓> |↑↓>+|↓↑>  S=0 |↑↓>-|↓↑>
CCC   S1(*)S2=(S^2-S1^2-S2^2)/2 := (0-3/4-3/4)/2,or,(2-3/4-3/4)/2
CCC   δ(x) = lim_{sigma->0} 1/(√2pi*σ) * exp(-x^2/(2*σ^2))
CCC   VSS(i,j,s)

C      VSS(i,j,1)=
C     &   CN(i)*CN(j)/(sqrt(2d0*pi)*sigma)*(ge*mu_b/hc)**2
C     & * (sqrt(pi)/(sqrt(b(i)+b(j)+1d0/(sqrt(2d0)*sigma)))**3)*2d0/12d0
C     & * 3d0/4d0

C      VSS(i,j,2)=
C     & - CN(i)*CN(j)/(sqrt(2d0*pi)*sigma)*(ge*mu_b/hc)**2
C     & * (sqrt(pi)/(sqrt(b(i)+b(j)+1d0/(sqrt(2d0)*sigma)))**3)*2d0/12d0
C     & * 1d0/4d0

      VSS(i,j,1)=
     &   CN(i)*CN(j)/(sqrt(2d0*pi)*sigma)*hc**2*4d0/6d0*mu_b**2*(-3)
     & * (sqrt(pi)/(sqrt(b(i)+b(j)+1d0/(sqrt(2d0)*sigma)))**3)/4d0
     & * masse/massp

      VSS(i,j,2)=
     &   CN(i)*CN(j)/(sqrt(2d0*pi)*sigma)*hc**2*4d0/6d0*mu_b**2*(1)
     & * (sqrt(pi)/(sqrt(b(i)+b(j)+1d0/(sqrt(2d0)*sigma)))**3)/4d0
     & * masse/massp

CCC   hamiltonian matrix

      H(i,j,1) = T(i,j) + VC(i,j) + VSS(i,j,1)
      H(i,j,2) = T(i,j) + VC(i,j) + VSS(i,j,2)
      H(j,i,1) = H(i,j,1)
      H(j,i,2) = H(i,j,2)

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
      write(*,'(F10.5,9(1X,F10.5))') (n0(i,j),j=1,nmax)
      end do

C      write(*,*)

CCC   new hamiltonian

      do i=1,nmax
       do j=1,nmax
        hvec0(i,j) = 0d0
        hvec1(i,j) = 0d0
        do k=1,nmax
         do l=1,nmax
         hvec0(i,j) = hvec0(i,j) + vec1(k,i)*H(k,l,1)*vec1(l,j)
         hvec1(i,j) = hvec1(i,j) + vec1(k,i)*H(k,l,2)*vec1(l,j)
         end do
        end do
       end do
      end do

      do i=1,nmax
C      write(*,'(F10.5,9(1X,F10.5))') (hvec0(i,j),j=1,nmax)
      end do
C      write(*,*)

      do i=1,nmax
C      write(*,'(F10.5,9(1X,F10.5))') (hvec1(i,j),j=1,nmax)
      end do
C      write(*,*)

CCC
CCC   HQR(HAMILTONIAN)
CCC

      call hqr(ndim,nmax,hvec0,evec2,ene1)
      call hqr(ndim,nmax,hvec1,evec3,ene2)

      do i=1,nmax
C      write(*,*) i,ene1(i)
      end do

      do i=1,nmax
C      write(*,*) i,ene2(i)
      end do

CCC   new base vactor

      do i=1,nmax
       do k=1,nmax
        vec3(i,k)=0d0
        vec4(i,k)=0d0
        do j=1,nmax
         vec3(i,k) = vec3(i,k) + evec2(j,k) * vec1(i,j)
         vec4(i,k) = vec4(i,k) + evec3(j,k) * vec4(i,j)
        end do
       end do
      end do

CCC   sort

      do k=1,nmax
       i=0
       do j=1,nmax
        if(ene1(k).ge.ene1(j)) i=i+1
       end do
       ene3(i)=ene1(k)
       do j=1,nmax
        vec5(j,i)=vec3(j,k)
       end do
      end do

      do k=1,nmax
       i=0
       do j=1,nmax
        if(ene2(k).ge.ene2(j)) i=i+1
       end do
       ene4(i)=ene2(k)
       do j=1,nmax
        vec6(j,i)=vec4(j,k)
       end do
      end do

      write(*,*) 'L=',ls
      do i=1,10
      write(*,*) i,ene3(i),ene4(i)
      end do
      write(*,*)

CCC   ground state

      dr = 0.005d0
      rmax = 10d0
      nrmax = rmax/dr + 1
      k=1

      do nr=1,nrmax
      r=dr*(nr-1)
      wf0=0d0
      wf1=0d0
      do i=1,nmax
       wf0 = wf0 + vec4(i,k)*CN(i)*r**(ls)*exp(-b(i)*r**2)
       wf1 = wf1 + vec6(i,k)*CN(i)*r**(ls)*exp(-b(i)*r**2)
      end do

C      write(*,*) r,wf0,wf1
      end do

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

