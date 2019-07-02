**********************************************************************
      SUBROUTINE OUTPUT(NN,z,delz,dz,y,limex_nstep)
**********************************************************************
*****                                                            *****
*****  macht die File-Ausgaben                                   *****
*****                                                            *****
**********************************************************************
      use drift_data,ONLY: NNUC,NEPS,NMOM,NDUST,NSPECIES,no_converge,
     &                     imp_nuklea,imp_growth,imp_drift,imp_misch,
     &                     eps,elnr,LL,T,rho,nH,p,pel,g,wmix,mu,rhod,
     &                     pi,bk,amu,bar,Vl,sizedist,chi,chinet,bmix,
     &                     dust_nel,dust_el,dust_vol,dust_nu,elcode,
     &                     eps0,nsp,Sat,Jst,Nst,Jev,Rnull,Lref,
     &                     NMOLE,cmol,nmol,nat,pred,dust_nam,elnam,
     &                     PHOENIX_MODE
      implicit none
      integer,intent(IN):: NN,limex_nstep
      integer ::  H,He,C,N,O,Ne,Na,Mg,Al,Si,S,K,Ca,Cr,Mn,Fe,Ni
      real*8,intent(IN) :: z,delz,dz,y(NN)
      real*8 :: f(NMOM+NDUST)
      real*8 :: fcond(NEPS),effSat(NDUST)
      real*8 :: nmono,cc,aquer,vdrquer,tau
      real*8 :: L0test,L0,L1,L2,L3,L4
      real*8 :: N1,N2,a1,a2,a3,sigma,AAA,BBB,CCC
      real*8 :: nTiO2,Jstern,Nstern,SS,Jsum,dep
      real*8 :: lamb,mukin,Kn,Re,Hp,dVdt,tau_sink,tau_grow,tau_nuc
      integer:: i,j,el,el2,sol,stindex
      integer,save :: Ti,TiO2
      logical :: changed
      logical,save :: firstCall=.true.
*     -----------------------------------------------------------------
      data H/1/, He/2/, C/6/, N/7/, O/8/, Ne/10/, Na/11/, Mg/12/,Al/13/
      data Si/14/, S/16/, K/19/, Ca/20/, Cr/24/, Mn/25/, Fe/26/
      data Ni/28/, Ti/22/
*     -----------------------------------------------------------------

      if (firstCall) then
        Ti = 22
        TiO2 = stindex(CMOL,NMOLE,'TIO2     ')
        firstCall=.false.
      endif    

      call FF(NN,z,y,f,-1)
      if (no_converge) then
        write(*,*) "*** no convergence in STATIC_EPS called from OUTPUT"  
        RETURN
      endif

      if (.not.PHOENIX_MODE) then
        open(90,file='restart.dat',status='unknown')
        rewind(90)
        write(90,*) z,dz,y(:)
        write(90,*) eps(elnr(:))
        write(90,*) pred(:)
        close(90)
      endif  

      !------------------------------------
      ! ***  effektive Uebersaettigung  ***
      !------------------------------------
c     nTiO2 = nmol(TiO2)
c     CALL KLASS_NUK_TIO2(T,nTiO2,Jstern,Nstern,SS)
      call EFF_SUPERSAT(nH,T,bmix,effSat)
      L0 = LL(0)
      L1 = LL(1)
      L2 = LL(2)
      L3 = LL(3)
      L4 = LL(4)
      if ((L0.lt.0.d0).or.(L1.lt.0.d0).or.
     &    (L2.lt.0.d0).or.(L3.lt.0.d0)) then
        write(*,*) 'negative Momente',LL(0:3)
        stop
      endif

      !-----------------------------------------------------
      ! ***  mittlere Teilchengroesse und -driftgeschw.  ***
      !-----------------------------------------------------
      cc = DSQRT(2.d0*bk*T/(mu*amu))
      if (L0.gt.0.d0) then
        aquer = (3.d0/(4.d0*pi))**(1.d0/3.d0) * L1/L0
      else
        aquer = (3.d0*Vl/(4.d0*pi))**(1.d0/3.d0)
      endif
      vdrquer = DSQRT(pi)/2.d0 * g*rhod/(rho*cc) * aquer
      
      !---------------------
      ! ***  Zeitskalen  ***
      !---------------------
      lamb  = 1.856d-9/rho                     ! mean free path
      mukin = 5.877d-6 * DSQRT(T)              ! kinenatische Viskositaet
      Kn    = lamb/(2.d0*aquer)                ! Knudsen number
      Re    = 2.d0*aquer*rho*vdrquer/mukin     ! Reynolds number
      Hp    = bk*T/(mu*amu*g)                  ! scale height
      dVdt  = 4.d0*pi*aquer**2*DABS(chinet)/(36.d0*pi)**(1.d0/3.d0)
      tau_sink = Hp/vdrquer                    ! sinking time scale
      tau_grow = 4.d0*pi/3.d0*aquer**3/dVdt    ! growth time scale
      Jsum = 0.0
      do  i=1,NNUC
       Jsum = Jsum + Jst(i)    
      enddo
      tau_nuc  = rho*L0/Jsum                   ! nucleation time scale

      !-----------------------------
      ! ***  Groessenverteilung  ***
      !-----------------------------
      if ((L1.gt.0.d0).and.(L2.gt.0.d0).and.(L3.gt.0.d0)) then
        if (sizedist.eq.1) then
          call DIST_FROM_MOMENTS1(rho,L1,L2,L3,L4,a1,a2,N1,N2,L0test)
        else if (sizedist.eq.2) then
          call DIST_FROM_MOMENTS2(rho,L1,L2,L3,L4,AAA,BBB,CCC,L0test)
        else if (sizedist.eq.3) then
          L0test = L0
          call DIST_FROM_MOMENTS3(rho,L1,L2,L3,L4,N1,a1,sigma,L0test,
     &                            sol)
        endif
      else
        a1 = 0.0
        a2 = 0.0
        N1 = 0.0
        N2 = 0.0
        AAA = 0.0
        BBB = 0.0
        CCC = 0.0
        sigma = 0.0
      endif  

      !------------------------------------------
      ! ***  Kondensationsgrade der Elemente  ***
      !------------------------------------------
      do i=1,NEPS
        fcond(i) = 0.d0
      enddo
      do i=1,NDUST
        nmono = bmix(i)*rho*L3/dust_vol(i)
        do j=1,dust_nel(i)
          el  = dust_el(i,j)
          el2 = elcode(el)
          fcond(el2) = fcond(el2) + nmono*dust_nu(i,j)/(nH*eps0(el))
        enddo
      enddo
c     tau = y(NN-1)
      tau = 0.0

      write(91,1000) z, T, nH, rho, p, wmix, tau, 
     &               (nsp(i),i=1,NSPECIES), pel
      write(92,1000) z, p, 
     &               (DMIN1(1.d+30,effSat(i)),i=1,NDUST), 
     &               Nst, DMAX1(1.d-99,Jst), Jev
      write(93,1000) z,p,L0,L1,L2,L3,aquer,vdrquer,
     &  (DMAX1(1.d-99, bmix(i)),i=1,NDUST),
     &                  (chi(i),i=1,NDUST),
     &  (DMAX1(1.d-99,fcond(i)),i=1,NEPS),
     &            (eps(elnr(i)),i=1,NEPS),
     &   DMAX1(1.d-99,rhod*L3), rhod, chinet
      write(94,1000) z, p, imp_nuklea, imp_growth, imp_drift, imp_misch,
     &               rho,T,aquer*1.d+4,vdrquer*1.d-2,Kn,
     &               DMIN1(1.d+30,tau_nuc),tau_grow,tau_sink,1.d0/wmix
      write(95,2000) Rnull-z, T, rho, p, L0, L1, L2, L3,
     &               ( bmix(i),i=1,NDUST),
     &               ( eps(elnr(i)),i=1,NEPS)
      write(100,2000) Rnull-z, T, rho, p, L0, L1, L2, L3,
     &               ( bmix(i),i=1,NDUST),
     &               ( eps(elnr(i)),i=1,NEPS)
      if (sizedist.eq.1) then
        if (a2.gt.a1) then
          a3=a1
          a1=a2
          a2=a3
        endif
        write(96,1000) p/bar, pel/bar, T, a1, a2, N1, N2
      else if (sizedist.eq.2) then
        write(96,1000) p/bar, pel/bar, T, AAA, BBB, CCC
      else if (sizedist.eq.3) then
        write(96,1002) p/bar, pel/bar, T, N1, a1, sigma, sol
      endif  
!     write(*,1010) y
      write(*,1001) p/bar, limex_nstep, delz, T, 
     &              Jst, chinet, aquer, L3/Lref(3)
      write(97,1000) z, p, T, nH, nmol(1:100)
      write(98,1000) z, p, T, nH, nmol(101:NMOLE)
      write(99,1000) z, p, T, nH, nat(He),nat(H),nat(C),nat(N),nat(O),
     & nat(Si),nat(Mg),nat(Al),nat(Fe),nat(S),nat(Na),nat(K),nat(Ti),
     & nat(Ca)

#ifndef IFORT
!      call FLUSH(91)
!      call FLUSH(92)
!      call FLUSH(93)
!      call FLUSH(94)
!      call FLUSH(95)
!      call FLUSH(96)
!      call FLUSH(97)
!      call FLUSH(98)
!      call FLUSH(99)
       call FLUSH(100)
#endif

      !---------------------------------------
      ! ***  Entsteifung der Verdampfung?  ***
      !---------------------------------------
      call DESTIFF(nH,T,bmix)

      RETURN
 1000 format(999(1pE12.4))
 1001 format(1pE10.3,i4,99(1pE10.3))
 1002 format(5(1pE13.5),i3)
 1010 format(99(1pE9.2))
 2000 format(99(1pE20.12))
      end


