*
* Array sizes
*
      integer matom,mgp,mion,mtar
      parameter(matom=20000,mgp=1e6,mion=100,mtar=100)
*
* I/O units
*
      integer luin,luout
      parameter(luin=10,luout=11)
*
* Constants
*
      double precision ang,eps,pi,qe
      parameter(ang=1.0d-10,eps=8.854187817e-12,pi=3.1415926536d0)
      parameter(qe=1.6021766208d-19)
*
* Grid related and other parameters
*
      integer np
      parameter(np=16)
      double precision big,rho,sigma,small
      parameter(big=1.0d20,rho=8.0d0,sigma=2.0d0,small=1.0d-1)
*
      double precision rho2,sigma2
      parameter(rho2=rho**2,sigma2=sigma**2)
