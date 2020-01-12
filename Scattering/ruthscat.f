CCC
CCC   RUTHERFORD SCATTERING
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
      real*8    hc0,ac,ec,pi,piad,hc,z1,z2,t,
     &          nthmax,thmax,thmin,dth,tht,rad,a,sig
      parameter (hc0=197.326901d0,ac=931.494013d0,ec=137.03599976d0,
     &           thmax=180.0d0,thmin=0.0d0,dth=1.0d0)
CCC   hc0[MeV*fm],ac[Mev],ec=hc/e^2,t[kinetic energy]

CCC
CCC 
      pi=acos(-1.0d0)
      piad=pi/180.0d0

      z1=4
      z2=78
      t=27.7 ![MeV]

      nmax =  (thmax-thmin)/dth + 1.01

      a = z1*z2*hc0/ec/t

      do i=1,nmax
      tht = (i-1) * dth + thmin
      rad = tht * piad

       if (rad.lt.1.0d-10) rad=1.0d-10

      sig = a**2 / (16.0d0*sin(rad/2.0d0)**4.0d0) * 10.0d0

      write(*,*) tht,sig

      end do

      stop
      end

