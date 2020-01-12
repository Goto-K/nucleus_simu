CCC
CCC   PLANE WAVE CALCULATION FOR NEUTRON ELASTIC SCATTERING
CCC   BY NUCLEON-NUCLEON DELTA INTERACTION
CCC
      program main
      implicit none
CCC
      integer   i,j,k,l
      integer   ndim,nmax
      parameter (ndim=1000)
      real*8,dimension(ndim,ndim)  :: 
     &          XX
      real*8,dimension(ndim)  :: 
     &          XXX
      real*8    hc0,ac,ec,pi,piad,hc,at,elab,r0,v0,
     &          nthmax,thmax,thmin,dth,tht,rad,fk,a,a2,sig,
     &          mp,rho0,c0,f,q,rabs,rmax,dr,Rp,rho,ad
      parameter (hc0=197.326901d0,ac=931.494013d0,ec=137.03599976d0,
     &           thmax=180.0d0,thmin=0.0d0,dth=1.0d0)
CCC   function
      real*8    fitg
CCCC   hc0[MeV*fm],ac[Mev],ec=hc/e^2,t[kinetic energy],rho0[1/fm^3]

CCC
CCC 
      pi=acos(-1.0d0)
      piad=pi/180.0d0

      rmax=50.0d0
      dr=0.1d0

      elab=65.0d0 ![MeV]
      v0=2500.0d0

CCC   208Pb
      at=208.0d0
      Rp=6.620d0 ![fm]
      rabs=0.0d0
      ad=0.546d0 ![fm]


      nmax =  (thmax-thmin)/dth + 1.01
      mp=1.0d0*ac
      fk = sqrt(2.0d0*mp*elab)/hc0
      rho0=3.0d0/(4.0d0*pi*Rp**3.0d0)*at
CCC   rho0 = m/(4/3pir^3)
      c0=mp**2/(2.0d0*pi*hc0**2)**2*v0**2*10.0d0
CCC   sig = mp/(2pih^2) * v0 * |rho(q)|^2
CCC       = c0 * f^2

      do i=1,nmax
      tht = (i-1) * dth + thmin
      rad = tht * piad
      q=2.0d0*fk*sin(rad/2.0d0)

      sig=fitg(rabs,rmax,dr,q,Rp,rho0,ad)**2*c0    ! IN UNIT OF MB

      write(*,*) tht,sig

      end do

      stop
      end
C
C
C
C=======================================================================
      FUNCTION FITG(Rabs,RMAX,DR,q,Rp,rho0,ad)
      IMPLICIT REAL*8 (A-H,O-Z)
C=======================================================================
C
CCC---------------------------------------------------------------------
C
      FITG=0.0D0
      NRMAX=(RMAX-Rabs)/DR+1.01
      DO IR=1,NRMAX
       R=(IR-1)*DR+Rabs
       IF(Q.LT.1.0D-10) THEN
         FITG=FITG+(rho0/(1+exp((r-Rp)*ad)))*r
        ELSE
         FITG=FITG+(rho0/(1+exp((r-Rp)*ad)))*SIN(Q*R)/Q
       END IF
      END DO
C
      IF(Q.LT.1.0D-10) THEN
        FITG=FITG-(rho0/(1+exp((rabs-Rp)*ad)))*Rabs*0.5D0
     1           -(rho0/(1+exp((rmax-Rp)*ad)))*RMAX*0.5D0
       ELSE
        FITG=FITG-(rho0/(1+exp((rabs-Rp)*ad)))*SIN(Q*Rabs)/Q*0.5D0
     1           -(rho0/(1+exp((rmax-Rp)*ad)))*SIN(Q*RMAX)/Q*0.5D0
      END IF
C
      FITG=FITG*DR
C
      RETURN
      END

