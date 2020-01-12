CCC   *
CCC   *     VARIATIONAL METHOD WITH THE GAUSSIAN BASIS FUNCTIONS
CCC   *
      program three_body_ho
C      Use nucleus_simu/lapack/aux_routines/
C     &    lapack_example_aux,Only:nagf_blas_damax_val,
C     &                            nagf_file_print_matrix_real_gen
C      Use Desktop/fortran/nucleus_simu/lapack/aux_routines/
C     &    lapack_interfaces,Only:ddisna,dsyev
C      Use Desktop/fortran/nucleus_simu/lapack/aux_routines/
C     &    lapack_precision,Only:dp
      implicit none
      integer   i,j,k,l,l1,l2,a1,a2
      integer   ndim,lmax,amax,nmax,omax,lwork,info
      parameter (ndim=1000,lmax=3,amax=5,
     &           nmax=6,omax=6,lwork=2*ndim-1)
      real*8,dimension(ndim,ndim)  :: 
     &          N1,N2,T,V,H,t1,hvec,evec1,evec2,
     &          cn1,cn2,
     &          vec0,vec1,vec2,vec3,vec4,vecs,tvec1
      real*8,dimension(ndim)  :: 
     &          bx,By,eig1,ene1,ene2
      real*8,dimension(3*ndim-2)  ::  rwork1,rwork2
      real*8    hc0,ac,ec,pi,hc,massp,massn,rem1,rem2,ls,
     &          r0,rmax,r,dr,x0,xmax,x,dx,y0,ymax,y,dy,wf
      complex*16,parameter :: ci=cmplx(0d0,1d0)
      integer   nx,nxmax,ny,nymax
      parameter (hc=197.326901d0,ac=931.494013d0,ec=137.03599976d0)
CCC   hc[MeV*fm],ac[Mev],ec=hc/e^2


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
CCC           =<1|V1|3>
CCC
CCC   <1|V2|1>=(1/4)*r1^(la+lb+4)*exp[-(bi+bj)*r1^2]
CCC           *R1^(La+Lb+4)*exp[-(Bi+Bj)*R1^2]
CCC           =<1|V3|1>
CCC
CCC   <1|V2|3>=(1/4)*(1/4)^qb*(9/16)^Qb*r1^(la+2qb+2Qb+4)
CCC                                    *exp[-(bi+(1/4)bj+(9/16)Bj)*r1^2]
CCC           *(1/4)^Qb*R1^(La+2lb+2Qb+4)*exp[-(bj+Bi+(1/4)Bj)*R1^2]
CCC           =<1|V3|2>=<1|V2|2>=<1|V3|3>
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC   l=0,L=0
CCC   <1|V1|1>=r1^(4)*exp[-(bi+bj)*r1^2]*R1^(2)*exp[-(Bo+Bp)*R1^2]
CCC
CCC   <1|V1|2>=r1^(4)*exp[-(bi+(1/4)bj+(9/16)Bp)*r1^2]
CCC           *R1^(2)*exp[-(bj+Bo+(1/4)Bp)*R1^2]
CCC           =<1|V1|3>
CCC
CCC   <1|V2|1>=(1/4)*r1^(4)*exp[-(bi+bj)*r1^2]
CCC           *R1^(4)*exp[-(Bo+Bp)*R1^2]
CCC           =<1|V3|1>
CCC
CCC   <1|V2|3>=r1^(4)*exp[-(bi+(1/4)bj+(9/16)Bp)*r1^2]
CCC           *R1^(4)*exp[-(bj+Bi+(1/4)Bj)*R1^2]
CCC           =<1|V3|2>=<1|V2|2>=<1|V3|3>
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC



      pi = acos(-1.0d0)
      hc0 = hc*10  ![eV*Å]
CCC   proton,neutron,redeuced mass      
C      massp = 1.00727646688d0 * ac 
C      massn = 1.00866491578d0 * ac
C      rem = massn*massn / (massn+massn)
      rem1 = ac**2 / (2d0*ac)
      rem2 = rem1*ac / (rem1+ac)
CCC   angular momentum
C      l,L = 0
CCC   range parameter
      x0 = 0.1d0
      xmax = 20d0
      y0 = 0.1d0
      ymax = 20d0

      dx = (xmax/x0)**(1.0d0/(nmax-1.0d0))
      dy = (ymax/y0)**(1.0d0/(omax-1.0d0))

CCC
CCC  calculation of matrix elements
CCC

      do i=1,nmax
      do k=1,omax

      bx(i) = 1.0d0 / (x0 * dx**(i-1))**2
      By(k) = 1.0d0 / (y0 * dy**(k-1))**2

      cn1(i,k) = sqrt(2d0**2/1d0 * sqrt((((2d0*bx(i))**3)/pi)))
     &         * sqrt(2d0**2/1d0 * sqrt((((2d0*By(k))**3)/pi)))
      cn2(i,k) = sqrt(2d0**2/1d0 
     &         * sqrt((((5d0/4d0*bx(i)+9d0/16d0*By(k))**3)/pi)))
     &         * sqrt(2d0**2/1d0
     &         * sqrt((((bx(i)+5d0/4d0*By(k))**3)/pi)))
   
CCC   norm matrix

      do j=1,i
      
      if (j.eq.i) then
       do l=1,k

      N1(omax*(i-1)+k,omax*(j-1)+l)

C      = r1^(2)*exp[-(bi+bj)*r1^2]*R1^(2)*exp[-(Bk+Bl)*R1^2]
C      + 2 * r1^(2)*exp[-(bi+(1/4)bj+(9/16)Bl)*r1^2]*R1^(2)*exp[-(bj+Bk+(1/4)Bl)*R1^2]
C    <1|1>+<1|2>+<1|3>

C     & = cn1(i,k)*cn1(j,l) 
     & = 1d0/2d0**2 * sqrt((pi)/((bx(i)+bx(j))**3)) 
     & * 1d0/2d0**2 * sqrt((pi)/((By(k)+By(l))**3))

     & + 2d0
C     & * cn2(i,k)*cn2(j,l)
     & * 1d0/2d0**2 * sqrt((pi)/((bx(i)+bx(j)/4d0+9d0/16d0*By(l))**3)) 
     & * 1d0/2d0**2 * sqrt((pi)/((By(k)+bx(j)+By(l)/4d0)**3))

      N2(omax*(i-1)+k,omax*(j-1)+l)=N1(omax*(i-1)+k,omax*(j-1)+l)

CCC   kinetic energy matrix
      
      T(omax*(i-1)+k,omax*(j-1)+l)

C    <1|T(r1)|1>+<1|T(r1)|2>+<1|T(r1)|3>

     & = hc0**2 / (2d0*rem1)
C     & = cn1(i,k)*cn1(j,l)*hc0**2 / (2d0*rem1)
     & * 3d0/2d0*bx(i)*bx(j) * sqrt((pi)/(bx(i)+bx(j))**5)
     & * 1d0/2d0**2 * sqrt((pi)/((By(k)+By(l))**3))
     & + 2d0 
     & * hc0**2 / (2d0*rem1)
C     & * cn2(i,k)*cn2(j,l)*hc0**2 / (2d0*rem1)
     & * 3d0/2d0*bx(i)*(bx(j)/4d0+9d0/16d0*By(l))
     & * sqrt((pi)/(bx(i)+bx(j)/4d0+9d0/16d0*By(l))**5)
     & * 1d0/2d0**2 * sqrt((pi)/((bx(j)+By(k)+By(l)/4d0)**3))

C    <1|T(R1)|1>+<1|T(R1)|2>+<1|T(R1)|3>

     & + hc0**2 / (2d0*rem2)
C     & + cn1(i,k)*cn1(j,l)*hc0**2 / (2d0*rem2)
     & * 3d0/2d0*By(k)*By(l) * sqrt((pi)/(By(k)+By(l))**5)
     & * 1d0/2d0**2 * sqrt((pi)/((bx(i)+bx(j))**3)) 
     & + 2d0
     & * hc0**2 / (2d0*rem2)
C     & * cn2(i,k)*cn2(j,l)*hc0**2 / (2d0*rem2)
     & * 3d0/2d0*By(k)*(bx(j)+By(l)/4d0)
     & * sqrt((pi)/(bx(j)+By(k)+By(l)/4d0)**5)
     & * 1d0/2d0**2 * sqrt((pi)/((bx(i)+bx(j)/4d0+9d0/16d0*By(l))**3)) 

CCC   potential energy matrix
CCC
CCC   V = V(r1) + V(r2) + V(r3)
CCC   omega=1/hc0
CCC   <1|V1|1> + 2*<1|V1|2> + 2*<1|V2|1> + 4*<1|V2|3>
      
      V(omax*(i-1)+k,omax*(j-1)+l)

     & = rem1 / (2*hc0**2)
C     & = cn1(i,k)*cn1(j,l)*rem1 / (2*hc0**2)
     & * 3d0/8d0 * sqrt((pi)/(bx(i)+bx(j))**5)
     & * 1d0/2d0**2 * sqrt((pi)/((By(k)+By(l))**3))

     & + 2d0
     & * rem1 / (2*hc0**2)
C     & * cn2(i,k)*cn2(j,l)*rem1 / (2*hc0**2)
     & * 3d0/8d0 * sqrt((pi)/(bx(i)+bx(j)/4d0+9d0/16d0*By(l))**5)
     & * 1d0/2d0**2 * sqrt((pi)/((bx(j)+By(k)+By(l)/4d0)**3))

     & + 2d0
C     & * rem1 / (2*hc0**2)
     & * cn1(i,k)*cn1(j,l)*rem1 / (2*hc0**2)
     & * 1d0/4d0 * 3d0/8d0 * sqrt((pi)/(bx(i)+bx(j))**5)
     & * 3d0/8d0 * sqrt((pi)/(By(k)+By(l))**5)

     & + 4d0
     & * rem1 / (2*hc0**2)
C     & * cn2(i,k)*cn2(j,l)*rem1 / (2*hc0**2)
     & * 1d0/4d0 * 3d0/8d0
     & * sqrt((pi)/(bx(i)+bx(j)/4d0+9d0/16d0*By(l))**5)
     & * 3d0/8d0 * sqrt((pi)/(bx(j)+By(k)+By(l)/4d0)**5)

CCC   hamiltonian matrix

      H(omax*(i-1)+k,omax*(j-1)+l)
     & = T(omax*(i-1)+k,omax*(j-1)+l)+V(omax*(i-1)+k,omax*(j-1)+l)

       end do

      else
       do l=1,omax

      N1(omax*(i-1)+k,omax*(j-1)+l)

C      = r1^(2)*exp[-(bi+bj)*r1^2]*R1^(2)*exp[-(Bk+Bl)*R1^2]
C      + 2 * r1^(2)*exp[-(bi+(1/4)bj+(9/16)Bl)*r1^2]*R1^(2)*exp[-(bj+Bk+(1/4)Bl)*R1^2]
C    <1|1>+<1|2>+<1|3>

     & = 1d0/2d0**2 * sqrt((pi)/((bx(i)+bx(j))**3)) 
C     & = cn1(i,k)*cn1(j,l)/2d0**2 * sqrt((pi)/((bx(i)+bx(j))**3)) 
     & * 1d0/2d0**2 * sqrt((pi)/((By(k)+By(l))**3))

     & + 2d0
C     & * cn2(i,k)*cn2(j,l)
     & * 1d0/2d0**2 * sqrt((pi)/((bx(i)+bx(j)/4d0+9d0/16d0*By(l))**3)) 
     & * 1d0/2d0**2 * sqrt((pi)/((By(k)+bx(j)+By(l)/4d0)**3))

      n2(omax*(i-1)+k,omax*(j-1)+l)=n1(omax*(i-1)+k,omax*(j-1)+l)

CCC   kinetic energy matrix
      
      T(omax*(i-1)+k,omax*(j-1)+l)

C    <1|T(r1)|1>+<1|T(r1)|2>+<1|T(r1)|3>

     & = hc0**2 / (2d0*rem1)
C     & = cn1(i,k)*cn1(j,l)*hc0**2 / (2d0*rem1)
     & * 3d0/2d0*bx(i)*bx(j) * sqrt((pi)/(bx(i)+bx(j))**5)
     & * 1d0/2d0**2 * sqrt((pi)/((By(k)+By(l))**3))
     & + 2d0 
     & * hc0**2 / (2d0*rem1)
C     & * cn2(i,k)*cn2(j,l)*hc0**2 / (2d0*rem1)
     & * 3d0/2d0*bx(i)*(bx(j)/4d0+9d0/16d0*By(l))
     & * sqrt((pi)/(bx(i)+bx(j)/4d0+9d0/16d0*By(l))**5)
     & * 1d0/2d0**2 * sqrt((pi)/((bx(j)+By(k)+By(l)/4d0)**3))

C    <1|T(R1)|1>+<1|T(R1)|2>+<1|T(R1)|3>

     & + hc0**2 / (2d0*rem2)
C     & + cn1(i,k)*cn1(j,l)*hc0**2 / (2d0*rem2)
     & * 3d0/2d0*By(k)*By(l) * sqrt((pi)/(By(k)+By(l))**5)
     & * 1d0/2d0**2 * sqrt((pi)/((bx(i)+bx(j))**3)) 
     & + 2d0
     & * hc0**2 / (2d0*rem2)
C     & * cn2(i,k)*cn2(j,l)*hc0**2 / (2d0*rem2)
     & * 3d0/2d0*By(k)*(bx(j)+By(l)/4d0)
     & * sqrt((pi)/(bx(j)+By(k)+By(l)/4d0)**5)
     & * 1d0/2d0**2 * sqrt((pi)/((bx(i)+bx(j)/4d0+9d0/16d0*By(l))**3)) 

CCC   potential energy matrix
CCC
CCC   V = V(r1) + V(r2) + V(r3)
CCC   omega=1/hc0
CCC   <1|V1|1> + 2*<1|V1|2> + 2*<1|V2|1> + 4*<1|V2|3>
      
      V(omax*(i-1)+k,omax*(j-1)+l)

     & = rem1 / (2*hc0**2)
C     & = cn1(i,k)*cn1(j,l)*rem1 / (2*hc0**2)
     & * 3d0/8d0 * sqrt((pi)/(bx(i)+bx(j))**5)
     & * 1d0/2d0**2 * sqrt((pi)/((By(k)+By(l))**3))

     & + 2d0
     & * rem1 / (2*hc0**2)
C     & * cn2(i,k)*cn2(j,l)*rem1 / (2*hc0**2)
     & * 3d0/8d0 * sqrt((pi)/(bx(i)+bx(j)/4d0+9d0/16d0*By(l))**5)
     & * 1d0/2d0**2 * sqrt((pi)/((bx(j)+By(k)+By(l)/4d0)**3))

     & + 2d0
     & * rem1 / (2*hc0**2)
C     & * cn1(i,k)*cn1(j,l)*rem1 / (2*hc0**2)
     & * 1d0/4d0 * 3d0/8d0 * sqrt((pi)/(bx(i)+bx(j))**5)
     & * 3d0/8d0 * sqrt((pi)/(By(k)+By(l))**5)

     & + 4d0
     & * rem1 / (2*hc0**2)
C     & * cn2(i,k)*cn2(j,l)*rem1 / (2*hc0**2)
     & * 1d0/4d0 * 3d0/8d0
     & * sqrt((pi)/(bx(i)+bx(j)/4d0+9d0/16d0*By(l))**5)
     & * 3d0/8d0 * sqrt((pi)/(bx(j)+By(k)+By(l)/4d0)**5)

CCC   hamiltonian matrix

      H(omax*(i-1)+k,omax*(j-1)+l)
     & = T(omax*(i-1)+k,omax*(j-1)+l)+V(omax*(i-1)+k,omax*(j-1)+l)

      end do

      end if

      end do
      end do
      end do

CCC
CCC
     
      do i=1,nmax*omax
      do j=1,i

      N1(j,i) = N1(i,j)
      N2(j,i) = N2(i,j)

      H(j,i) = H(i,j)

      end do
      end do

CCC
CCC
CCC   DSYEV(NORM)
CCC

      call dsyev('V','U',nmax*omax,N1,ndim,eig1,rwork1,lwork,info)
      write(*,*) info

CCC   |f> = Σ_a evec/√eig |fai>

      do i=1,nmax*omax
       do j=1,nmax*omax
        vecs(i,j)=0
        do k=1,nmax*omax
         vecs(i,j) = vecs(i,j) + N1(i,k) * N1(j,k)
        end do
       end do
      end do

      do i=1,nmax*omax
C      write(*,'(F10.5,9(1X,F10.5))') (vecs(j,i),j=1,nmax*omax)
      end do
      
CCC   standardization

      do i=1,nmax*omax
       do j=1,nmax*omax
        vec1(i,j) = N1(i,j) / sqrt(eig1(j))
       end do
      end do

      do i=1,nmax*omax
       do j=1,nmax*omax
        vec2(i,j) = 0d0
        do k=1,nmax*omax
         do l=1,nmax*omax
         vec2(i,j) = vec2(i,j) + vec1(k,i)*N2(k,l)*vec1(l,j)
         end do
        end do
       end do
      end do

      do i=1,nmax*omax
C      write(*,'(F10.5,9(1X,F10.5))') (N2(i,j),j=1,nmax*omax)
C      write(*,*) (N2(i,j),j=1,nmax*omax)
C      write(*,'(2(F10.5))') (vec1(i,j),j=1,nmax*omax)
C      write(*,*) eig1(i)
      end do

C      write(*,*)

CCC   new hamiltonian

      do i=1,nmax*omax
       do j=1,nmax*omax
        hvec(i,j) = 0d0
        do k=1,nmax*omax
         do l=1,nmax*omax
         hvec(i,j) = hvec(i,j) + vec1(k,i)*H(k,l)*vec1(l,j)
         tvec1(i,j) = tvec1(i,j)
     &              + vec1(k,i)*H(k,l)*vec1(l,j)
         end do
        end do
       end do
      end do

      do i=1,nmax*omax
C      write(*,'(F10.5,9(1X,F10.5))') (hvec(i,j),j=1,nmax*omax)
      end do
      do i=1,nmax*omax
C      write(*,'(F10.5,9(1X,F10.5))') (tvec1(i,j),j=1,nmax*omax)
C      write(*,'(2(F10.5)') (tvec1(i,j),j=1,nmax*omax)
      end do
C      write(*,*)

CCC
CCC   DSYEV(HAMILTONIAN)
CCC

      call dsyev('V','U',nmax*omax,hvec,ndim,ene1,rwork2,lwork,info)

      do i=1,nmax**2
C      write(*,*) i,ene1(i)
      end do

CCC   new base vactor

      do i=1,nmax*omax
       do k=1,nmax*omax
        vec3(i,k)=0d0
        do j=1,nmax*omax
         vec3(i,k) = vec3(i,k) + hvec(j,k) * vec1(i,j)
        end do
       end do
      end do

CCC   sort

      do k=1,nmax*omax
       i=0
       do j=1,nmax*omax
        if(ene1(k).ge.ene1(j)) i=i+1
       end do
       ene2(i)=ene1(k)
       do j=1,nmax*omax
        vec4(j,i)=vec3(j,k)
       end do
      end do

      do i=1,5
      write(*,*) i,ene2(i)
      end do

CCC   ground state

      dx = 0.005d0
      xmax = 7d0
      nxmax = xmax/dx + 1
      dy = 0.005d0
      ymax = 7d0
      nymax = ymax/dy + 1
      k=1

      do nx=1,nxmax
      do ny=1,nymax
      x=dx*(nx-1)
      y=dy*(ny-1)
      wf=0d0
      do i=1,nmax
      do j=1,omax
       wf=wf + vec4((i-1)*omax+j,k)*cn1(i,j)
     &       * exp(-bx(i)*x**2)*exp(-By(j)*y**2)
      end do
      end do

C      write(*,*) x,y,wf
      end do
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

