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
     &          nthmax,thmax,thmin,dth,tht,rad,k,a,a2,sig,
     &          mp,rho0,c0,f,q
      parameter (hc0=197.326901d0,ac=931.494013d0,ec=137.03599976d0,
     &           thmax=180.0d0,thmin=0.0d0,dth=1.0d0)
CCC   function
      real*8    fitg
CCC   hc0[MeV*fm],ac[Mev],ec=hc/e^2,t[kinetic energy],rho0[1/fm^3]

CCC
CCC 
      pi=acos(-1.0d0)
      piad=pi/180.0d0

      at=12.0d0 ![MeV]
      elab=65.0d0 ![MeV]
      r0=3.8d0
      v0=268.0d0

      nmax =  (thmax-thmin)/dth + 1.01
      mp=1.0d0*ac
      k = sqrt(2.0d0*mp*elab)/hc0
      rho0=3.0d0/(4.0d0*pi*r0**3.0d0)*at
CCC   rho0 = m/(4/3pir^3)
      c0=mp**2/(2.0d0*pi*hc0**2)**2*v0**2*10.0d0
CCC   sig = mp/(2pih^2) * v0 * |rho(q)|^2
CCC       = c0 * f^2

      do i=1,nmax
      tht = (i-1) * dth + thmin
      rad = tht * piad
      q=2.0d0*k*sin(rad/2.0d0)

       if(q.le.1.0d-7) then
         f=at
       else
         f=rho0*4.0d0*pi*(-r0*q*cos(r0*q)+sin(r0*q))/q**3
       end if

      sig=f**2*c0    ! IN UNIT OF MB

      write(*,*) tht,sig

      end do

      stop
      end

