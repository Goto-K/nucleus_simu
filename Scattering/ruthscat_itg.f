CCC
CCC   RUTHERFORD SCATTERING
CCC   CALCULATED BY NUMERICAL INTEGRATION
CCC   WITH SCREENED COULOMB POTENTIAL
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
      real*8    hc0,ac,ec,pi,piad,hc,z1,z2,am1,t,
     &          nthmax,thmax,thmin,dth,tht,rad,fk,a,a2,sig,
     &          rscr,rmax,dr,q
      parameter (hc0=197.326901d0,ac=931.494013d0,ec=137.03599976d0,
     &           thmax=180.0d0,thmin=0.0d0,dth=1.0d0)
CCC   function
      real*8    fitg
CCC   hc0[MeV*fm],ac[Mev],ec=hc/e^2,t[kinetic energy],rscr[screening radius]

CCC
CCC 
      pi=acos(-1.0d0)
      piad=pi/180.0d0

      z1=4.0d0
      z2=78.0d0
      t=27.7d0 ![MeV]
      am1=4.0d0
      rscr=5.0d0 ![fm]
      rmax=50.0d0
      dr=0.1d0

      nmax =  (thmax-thmin)/dth + 1.01
      fk = sqrt(2.0d0*am1*ac*t)/hc0
      a = z1*z2*hc0/ec*am1*ac*2.0d0/hc0**2
      a2 = a**2

      do i=1,nmax
      tht = (i-1) * dth + thmin
      rad = tht * piad
      q=2.0d0*fk*sin(rad/2.0d0)
C
      sig = a2 * fitg(0.0d0,rmax,dr,q,rscr)**2 * 10.0d0

      write(*,*) tht,sig

      end do

      stop
      end
C
C
C
C=======================================================================
      FUNCTION FITG(RMIN,RMAX,DR,Q,RSCR)
      IMPLICIT REAL*8 (A-H,O-Z)
C=======================================================================
C
CCC---------------------------------------------------------------------
C
      FITG=0.0D0
      NRMAX=(RMAX-RMIN)/DR+1.01
      DO IR=1,NRMAX
       R=(IR-1)*DR+RMIN
       IF(Q.LT.1.0D-10) THEN
         FITG=FITG+EXP(-R/RSCR)*R
        ELSE
         FITG=FITG+EXP(-R/RSCR)*SIN(Q*R)/Q
       END IF
      END DO
C
      IF(Q.LT.1.0D-10) THEN
        FITG=FITG-EXP(-RMIN/RSCR)*RMIN*0.5D0
     1           -EXP(-RMAX/RSCR)*RMAX*0.5D0
       ELSE
        FITG=FITG-EXP(-RMIN/RSCR)*SIN(Q*RMIN)/Q*0.5D0
     1           -EXP(-RMAX/RSCR)*SIN(Q*RMAX)/Q*0.5D0
      END IF
C
      FITG=FITG*DR
C
      RETURN
      END
