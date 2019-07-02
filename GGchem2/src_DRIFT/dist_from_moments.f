*-----------------------------------------------------------------------
      SUBROUTINE DIST_FROM_MOMENTS1(rho,L1,L2,L3,L4,a1,a2,N1,N2,L0)
*-----------------------------------------------------------------------
***  calculates the free coefficients a1,a2,N1,N2 of a double        ***
***  delta-peaked size distribution function of dust particles       ***
***                                                                  ***
***  f(a) = N1*Dirac(a-a1) + N2*Dirac(a-a2)                          ***
***                                                                  ***
***  INPUT:  rhoLj = dust moments \int f(V) V^{j/3} dV [cm^{j-3}]    ***
***  OUTPUT: a1,a2 = peak centers  [cm]                              ***
***          N1,N2 = dust particle densities [cm^-3]                 ***
***          L0    = estimated L0 from fit                           ***
*-----------------------------------------------------------------------
      implicit none
      real*8,intent(IN)  :: rho,L1,L2,L3,L4
      real*8,intent(OUT) :: a1,a2,N1,N2,L0
      real*8,parameter :: pi=3.141592653589793d0
      real*8,save  :: const(4)
      real*8  :: K0,K1,K2,K3,K4,den,pp,qq,dum,VIETA1,VIETA2,testa1
      integer :: j
      logical :: ok,IS_NAN
      logical,save :: firstCall=.true.
*-----------------------------------------------------------------------
*  Die Formelfunktion zur Loesung quadratischer Gleichungen mit Vieta
      VIETA1(pp,qq) = qq/(-pp/2.d0-DSQRT(pp**2/4.d0-qq))
      VIETA2(pp,qq) = qq/(-pp/2.d0+DSQRT(pp**2/4.d0-qq))
*-----------------------------------------------------------------------

      if (firstCall) then
        do j=1,4
          const(j) = (3.d0/(4.d0*pi))**(DBLE(j)/3.d0)
        enddo  
        firstCall=.false.
      endif
      
      K1 = const(1)*rho*L1
      K2 = const(2)*rho*L2
      K3 = const(3)*rho*L3
      K4 = const(4)*rho*L4
      
      den = K2**2-K1*K3

      ok  = (den.ne.0.d0)
      if (ok) then
        pp = (K1*K4 - K2*K3) / den
        qq = (K3**2 - K4*K2) / den
        if (pp.gt.0.d0) then
          a2 = VIETA1(pp,qq)
        else  
          a2 = VIETA2(pp,qq)
        endif  
        a1 = (a2*K2-K3) / (a2*K1-K2)
        N1 = (a2*K1-K2)**3 / (a2*K2-K3) / (K3-2.d0*a2*K2+a2**2*K1)
        N2 = (K3*K1-K2**2) / (K3-2.d0*a2*K2+a2**2*K1) / a2
        ok = ((a1.gt.0.d0).and.(a2.gt.0.d0).and.
     &        (N1.gt.0.d0).and.(N2.gt.0.d0).and.(pp**2/4.d0.gt.qq))
        ok = ok.and.(.not.IS_NAN(N1)).and.(.not.IS_NAN(N2))
        ok = ok.and.(N1-1.ne.N1).and.(N2-1.ne.N2)
      endif  
      if (.not.ok) then
        a1 = K2/K1
        N1 = K1**2/K2
        a2 = a1
        N2 = 0.d0
      endif  
      if (N2.gt.N1) then
        dum = a1
        a1  = a2
        a2  = dum
        dum = N1
        N1  = N2
        N2  = dum
      endif    

      K0 = N1+N2
      L0 = K0/rho

c     --- verify ---      
c     write(*,*) a1,a2, N1,N2
c     write(*,*) K1, a1*N1    + a2   *N2
c     write(*,*) K2, a1**2*N1 + a2**2*N2
c     write(*,*) K3, a1**3*N1 + a2**3*N2
c     write(*,*) K4, a1**4*N1 + a2**4*N2

      RETURN
      end


*-----------------------------------------------------------------------
      SUBROUTINE DIST_FROM_MOMENTS2(rho,L1,L2,L3,L4,A,B,C,L0)
*-----------------------------------------------------------------------
***  calculates the coefficients A,B,C of the                        ***
***  size distribution function of dust particles according to       ***
***                                                                  ***
***  f(a) = a^B * exp( A - C*a )   [cm^-4]                           ***
***                                                                  ***
***  INPUT:  rhoLj = dust moments \int f(V) V^{j/3} dV [cm^{j-3}]    ***
***  OUTPUT: A,B,C = fit coefficients                                ***
***          L0    = estimated L0 from fit                           ***
*-----------------------------------------------------------------------
      implicit none
      real*8,intent(IN)  :: rho,L1,L2,L3,L4
      real*8,intent(OUT) :: A,B,C,L0
      real*8,parameter :: pi=3.141592653589793d0
      real*8,save  :: const(4)
      real*8  :: K0,K1,K2,K3,K4,GAMMALN
      integer :: j
      logical :: ok
      logical,save :: firstCall=.true.

      if (firstCall) then
        do j=1,4
          const(j) = (3.d0/(4.d0*pi))**(DBLE(j)/3.d0)
        enddo  
        firstCall=.false.
      endif

      K1 = const(1)*rho*L1
      K2 = const(2)*rho*L2
      K3 = const(3)*rho*L3
      K4 = const(4)*rho*L4

      B = (2.d0*K1*K3-3.d0*K2**2) / (K2**2-K1*K3)
      ok = (B>0.d0)
      if (.not.ok) then
        B = 1.d+10
      endif
      B = DMIN1(B,1.d+10)
      C = (2.d0+B)*K1/K2
      A = LOG(K1)+(2.d0+B)*LOG(C)-GAMMALN(2.d0+B)

      K0 = exp(A+GAMMALN(B+1.d0)-(B+1.d0)*LOG(C))

c     --- verify ---      
c     write(*,*) A,B,C
c     write(*,*) rho*L0,K0
c     write(*,*) K1, exp(A+GAMMALN(B+2.d0)-(B+2.d0)*LOG(C)) 
c     write(*,*) K2, exp(A+GAMMALN(B+3.d0)-(B+3.d0)*LOG(C)) 
c     write(*,*) K3, exp(A+GAMMALN(B+4.d0)-(B+4.d0)*LOG(C)) 
c     write(*,*) K4, exp(A+GAMMALN(B+5.d0)-(B+5.d0)*LOG(C)) 

      L0 = K0/rho

      RETURN
      end


*-----------------------------------------------------------------------
      SUBROUTINE DIST_FROM_MOMENTS3(rho,L1,L2,L3,L4,N,a1,sig,L0,sol)
*-----------------------------------------------------------------------
***  calculates the coefficients A,a1,N of the                       ***
***  size distribution function of dust particles according to       ***
***                                                                  ***
***  f(a) = N / (SQRT(pi)*sig) * exp(-((a-a1)/sig)^2)   [cm^-4]      ***
***                                                                  ***
***  INPUT:  rhoLj = dust moments \int f(V) V^{j/3} dV [cm^{j-3}]    ***
***  OUTPUT: N     = total particle density [cm^-3]                  ***
***          a1    = mean particle size [cm]                         ***
***          sig   = standard deviation [cm]                         ***
***          L0    = estimated L0 from fit                           ***
*-----------------------------------------------------------------------
      implicit none
      real*8,intent(IN)  :: rho,L1,L2,L3,L4
      real*8,intent(OUT) :: N,a1,sig
      real*8,intent(INOUT) :: L0
      integer,intent(OUT) :: sol
      real*8,parameter :: pi=3.141592653589793d0
      real*8,save  :: const(4)
      real*8  :: K0,K1,K2,K3,K4,C,den,pp,qq,VIETA1,VIETA2
      real*8  :: a2,N1,N2,C1,C2,qual1,qual2
      integer :: j
      logical :: ok
      logical,save :: firstCall=.true.
*-----------------------------------------------------------------------
*  Die Formelfunktion zur Loesung quadratischer Gleichungen mit Vieta
      VIETA1(pp,qq) = qq/(-pp/2.d0-DSQRT(pp**2/4.d0-qq))
      VIETA2(pp,qq) = qq/(-pp/2.d0+DSQRT(pp**2/4.d0-qq))
*-----------------------------------------------------------------------

      if (firstCall) then
        do j=1,4
          const(j) = (3.d0/(4.d0*pi))**(DBLE(j)/3.d0)
        enddo  
        firstCall=.false.
      endif

      K1 = const(1)*rho*L1
      K2 = const(2)*rho*L2
      K3 = const(3)*rho*L3
      K4 = const(4)*rho*L4

      den= K3
      pp = -3.d0*K1*K2/den
      qq = 2.d0*K1**3/den
      N  = VIETA2(pp,qq)
      a1 = K1/N
      C  = DSQRT(N/2.d0/(K2-N*a1**2))
      ok = (N.gt.0.d0).and.(pp**2/4.d0.gt.qq).and.(K2.gt.N*a1**2)
      if (.not.ok) then
        N  = VIETA1(pp,qq)
c       a1 = K1/N
c       C  = DSQRT(N/2.d0/(K2-N*a1**2))
c       ok = (N.gt.0.d0).and.(pp**2/4.d0.gt.qq).and.(K2.gt.N*a1**2)
        if (.not.ok) then
          a1 = const(1)*L1/L0
          N  = K1/a1
          C  = DSQRT(N/2.d0/(K2-N*a1**2))
          if (K2.le.N*a1**2) C=1.d+99
          sol = 3
c         write(*,*) "no proper solution"          
        else  
          sol = 2
c         write(*,*) "solution 2"          
        endif  
      else  
        sol = 1  
c       write(*,*) "solution 1"
      endif  
      sig= 1.0/C
      K0 = N

c     --- verify ---      
c     write(*,*) N,a1,sig
c     write(*,*) rho*L0,K0
c     write(*,*) K1, N*a1
c     write(*,*) K2, N*a1**2 + N/(2.d0*C**2)
c     write(*,*) K3, N*a1**3 + 3.d0*N*a1/(2.d0*C**2) 
c     write(*,*) K4, N*a1**4 + 3.d0*N/(4.d0*C**4) + 3.d0*N*a1**2/C**2 

      L0 = K0/rho

      RETURN
      end


*-----------------------------------------------------------------------
      REAL*8 FUNCTION GAMMALN(x)
*-----------------------------------------------------------------------
***   nat.log. of the GAMMA-Function, source: Num.Recipies           ***
***   GAMMA(x+1) = x*GAMMA(x), GAMMA(n+1) = n!                       ***
*-----------------------------------------------------------------------
      integer j
      real*8 x,y,ser,tmp
      real*8,save :: cof(6),stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     & 24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     & -.5395239384953d-5,2.5066282746310005d0/
      y   = x
      tmp = x+5.5d0
      tmp = (x+0.5d0)*DLOG(tmp)-tmp
      ser = 1.000000000190015d0
      do j=1,6
        y = y+1.d0
        ser = ser+cof(j)/y
      enddo 
      GAMMALN = tmp + DLOG(stp*ser/x)
      RETURN
      end
