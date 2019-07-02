*************************************************************************
      PROGRAM static_weather

c  21.01.2008: FeH in smchem  ChH 
c  23.01.2009: CaSiO3 in supersat.f ChH 
c  23.11.2009: Na[s]  in supersat.f ChH 
*************************************************************************
      use drift_data,ONLY: NMOM,NDUST,NEPS,Teff,logg,eps0,beta,
     &                     muH,Rnull,Rlay,zlay,play,Nlayers,
     &                     amu,bk,pi,T,nH,rho,g,mu,wmix,LL,chinet,chi,
     &                     PHOENIX_MODE,firstStaticEpsCall,restart,
     &                     epsmerk,elfunk_calls,
     &                     sizedist,welche,abschluss,no_converge,
     &                     Nl,Vl,bar,cref,aref,Vref,Lref,pref,
     &                     dust_nam,dust_vol,pred
      implicit none
      integer,parameter:: NN=NMOM+NDUST
      real*8,parameter :: tol= 3.e-4
      real*8 :: V0,muref,epsref,Hp,znull,zende,pnull
      real*8 :: z1,z2,dz,dzmin,emerk(NEPS),z1merk,z2merk,dzmerk
      logical:: ex,abbruch,ok
      integer:: i,n,Al2O3,IRON,limex_nstep,startlay
      character(len=1) :: char1
      COMMON /LIMEX_NSTEP/ limex_nstep

* ----------------------------------------
* ***  variables for the LIMEX Solver  ***
* ----------------------------------------
      external                     :: FCN, JACOBI
      real*8, DIMENSION(NN)        :: y,ys,ymerk
      real*8, DIMENSION(NN), save  :: rtol,atol
      integer, DIMENSION(30), save :: iopt
      real*8, DIMENSION(5), save   :: ropt
      integer, DIMENSION(NN), save :: ipos
      integer, DIMENSION(3)        :: ifail
      integer :: nstep,nfcall1,nfcall2,njac
*
* ------------------------------
* ***  Parameter (vor INIT)  ***
* ------------------------------
      beta = 2.2d0          ! Overshoot parameter d log wmix / d log p 
      abschluss = 8         ! Wahl der Abschlussbedingung, 1-8
      sizedist  = 1         ! Wahl der analytischen Darstellung von f(a)
                            ! 1 -  f(a) = N1*Dirac(a-a1) + N2*Dirac(a-a2)
                            ! 2 -  f(a) = a^B * exp( A - C*a )   [cm^-4]
                            ! 3 -  f(a) = N / (SQRT(pi)*sig) * exp(-((a-a1)/sig)^2)

* -----------------------------------------------------------------------
* welche = 1 -- PHOENIX Structure (function of r: von innen nach aussen!) 
*        = 2 -- AMES    Structure (function of z: von aussen nach innen!)
*        = 3 -- test cases 2006 workshop
*        = 4 -- MARCS
* -----------------------------------------------------------------------
      welche = 4
      if (welche.eq.1) then 
        write(*,*) 'goto PHOENIX'
        call INIT_PHOENIX(beta) 
      else if (welche.eq.2) then 
        write(*,*)  'goto AMES'
        call INIT_AMES(beta)    
      else if (welche.eq.3) then 
        write(*,*) 'goto Test Cases'
        call INIT_TestCases(beta)   
      else if (welche.eq.4) then
        write(*,*) 'goto MARCS'
        call INIT_MARCS(beta)         
      endif

      call INIT(eps0)
      call INIT_DUSTCHEM
      Al2O3 = 0
      IRON  = 0
      do i=1,NDUST
        if (dust_nam(i).eq."Al2O3[s]") Al2O3=i 
        if (dust_nam(i).eq."Fe[s]") IRON=i 
        !write(*,*) dust_nam(i),Al2O3
      enddo
      if ((Al2O3.eq.0).or.(IRON.eq.0)) then
        write(*,*) "*** dust species Al2O3 or Fe not found."
        stop
      endif  

* ------------------------------
* ***  Parameter (nach INIT)  ***
* ------------------------------
      V0 = dust_vol(1)      ! TiO2-Monomervolumen
      Nl = 1000 
      Vl = DBLE(Nl)*V0      ! untere Integrationsgrenze fuer Momente       
*
* --------------------------
* ***  Referenzgroessen  ***
* --------------------------
      muref   = 2.35*amu                  ! H2-reich
      pref    = 1.d+6                     ! 1 bar
      cref    = DSQRT(2.d0*bk*Teff/muref) ! bei Teff
      epsref  = 0.5*eps0(12)              ! Mg-Haeufigkeit/2
      Vref    = 7.d-23                    ! MonomerVolumen von Mg2SiO4
      Lref(3) = epsref*Vref/muH           ! typisches Staubvolumen/g
      aref    = 1.d-4                     ! 1mic typischer Staubkornradius
      Vref    = 4.d0*pi/3.d0*aref**3       
      do i=0,NMOM
        Lref(i) = Lref(3) * Vref**(DBLE(i-3)/3.d0)
      enddo
*
* --------------------
* ***  Startwerte  ***
* --------------------
      Hp = bk*Teff/(10.d0**logg*muref)
      if (welche.eq.1) then          ! new PHOENIX  -- function of R
        startlay = 2  
        Rnull = Rlay(startlay)    
        znull = 0.d0
        zende = Rnull-Rlay(Nlayers)
      else  if (welche.eq.2) then    ! 2000 AMES    -- function of z
        startlay = 2
        Rnull = zlay(startlay)
        znull = 0.d0
        zende = zlay(Nlayers)
      else  if (welche.eq.3) then    ! TestCases    -- function of z
        startlay = 2
        Rnull = zlay(startlay)
        znull = 0.d0
        zende = zlay(Nlayers)
      else  if (welche.eq.4) then    ! MARCS - function of R
        startlay = 2  
        Rnull = Rlay(startlay)    
        znull = 0.d0
        zende = Rnull-Rlay(Nlayers)      
      endif 
      dzmin = (zende-znull)*1.d-9    ! minimum stepsize
      dz    = dzmin                  ! start stepsize
      call THERMO(znull,pnull,T,nH,rho,g,mu,wmix,welche)
      do i=1,NMOM
        y(i) = 0.d0                  ! staubfrei
      enddo
      do i=1,NDUST
        y(i+NMOM) = 0.d0
      enddo  
c     y(NN-1) = 0.d0                 ! optische Tiefe des Staubes
c     y(NN)   = pnull/pref
      write(*,*) "intergration starts at p=",play(startlay)/bar,"bar"
      write(*,*)

* -------------------
* ***  restart ?  ***
* -------------------
      PHOENIX_MODE = .false.
      firstStaticEpsCall = .true.
      restart=.false.  
      if (.not.PHOENIX_MODE) then
        inquire(file='restart.dat', exist=ex)
        if (ex) then
          open(90,file='restart.dat',status='old')
          read(90,*) znull,dz,y(:)
          read(90,*) epsmerk(:)
          read(90,*) pred(:)
          close(90)
          write(*,*) "RESTART from restart.dat ..." 
          write(*,*) "znull=",znull
          write(*,*)
          firstStaticEpsCall = .false.
          restart=.true.  
        endif  
      endif  

      call KOPF

* -----------------------------------------
* ***  parameters for the LIMEX-solver  ***
* -----------------------------------------
      iopt(1)  = 0       ! how much output? (0=no, 1=standard, 2=more)
      iopt(2)  = 0       ! unit number for output (0:default=6)
      iopt(3)  = 0       ! solution output? (0=no)
      iopt(4)  = 0       ! unit number for solution output (0:default=6)
      iopt(5)  = 1       ! nonsigular matrix BB (0=singular, 1=nonsingular)
      iopt(6)  = 0       ! determination of initial values for FF,BB (1=yes)
      iopt(7)  = 0       ! analytic Jacobian? (0=numerical, 1=analytic)
      iopt(8)  = NN      ! Lower bandwidth of the Jacobian (Ndim=full)
      iopt(9)  = NN      ! Upper bandwidth of the Jacobian (Ndim=full)
      iopt(10) = 0       ! reuse of the Jacobian? (1=yes) 
      iopt(11) = 1       ! Switch for error toleranz (0=scalars, 1=vectors)
      iopt(12) = 0       ! Return after one integration time step? (0=no)
      iopt(13) = 0       ! Dense output option (0=no dense output)
      iopt(14) = 0       ! The number of equidistant dense output points
      iopt(15) = 0       ! unit for dense output
      iopt(16) = 0       ! type of call (0=first call, 1=continuation call)
      iopt(17) = 0       ! bevahior at t1 (0=exactly up to t1)
      iopt(18) = 0       ! Generation of of a PostScript plot of the Jacobian
      ropt(1)  = Hp/200. ! maximum allowed stepsize 
      ropt(2)  = 0.0     ! maximal distance of dense outputs (if iopt(13)=3)
      ropt(3)  = 0.0     ! upper limit for t (if iopt(17)=1)
      do i=1,NN
        ipos(i) = 1       ! prevent YY(i)<0 if ipos(i)=1
        rtol(i) = tol     ! relative tolerance
        atol(i) = 1.d-25  ! absolute tolerance
      enddo
      no_converge = .false.
c     ys(:) = 0.0                ! don't provide ys

      call FF(NN,znull,y,ys,-2)  ! provide ys
      if (no_converge) then
        write(*,*) "*** no congergence in STATIC_EPS at initial height." 
        write(*,*) "***   ===> abort program." 
        stop
      endif

* ----------------------------------------
* ***  the main loop with LIMEX calls  ***
* ----------------------------------------
      z1 = znull
      do 
        z2 = z1 + DMIN1(2.d-3*(zende-znull), 10.d0*dz)
        if (y(NMOM+Al2O3)+y(NMOM+IRON).gt.0.7*y(4)) then
          z2 = z1 + DMIN1(5.d-4*(zende-znull), 10.d0*dz)
        endif  
        ymerk(:) = y(:)
        emerk(:) = epsmerk(:)
        z1merk   = z1
        z2merk   = z2
        dzmerk   = dz
        ropt(1)  = z2-z1                       ! maximum allowed stepsize 
        limex_nstep = 0
        do i=1,NN
          if (y(i).gt.0.d0) then
            atol(i)=tol*y(i)  
          endif  
          if ((y(i).lt.1.d-99).and.(y(i)+ys(i)*dz.gt.0.d0)) then
            atol(i)=tol*(y(i)+ys(i)*dz)
          endif  
c         if (i.gt.NMOM) atol(i)=atol(4)
c         write(*,*) y(i),atol(i)
        enddo
        no_converge = .false.

        !iopt(3) = 1
        !iopt(1) = 2
        !write(*,*) z1,z2,dz

        call LIMEX ( NN, FCN, JACOBI, z1, z2, y, ys, 
     &               rtol, atol, dz, iopt, ropt, ipos, ifail )
        !write(*,*) z1,z2,dz,ifail(1),no_converge
        !if (T.gt.1130.0) read(*,'(a1)') char1

        ok = .false.
        n  = ifail(1)
        if ((.not.no_converge).and. 
     &      ((n.eq.-46).or.(n.eq.-48).or.(n.eq.-50))) then
          !-------------------------------------------
          !***  Loesung ist ok, nur nicht fertig:  ***
          !***  z1 = last successful point         ***
          !***  reduziere interne Schrittweite     *** 
          !-------------------------------------------
          ifail(1) = 0             
          dz = DMIN1(0.2*dzmerk,0.1*(z1-z1merk))          
          if (dz<dzmin) ifail(1)=n
        endif  
        if ((ifail(1).eq.0).and.(.not.no_converge)) then
          !----------------------------------------
          !***  der Normalfall: Loesung ist ok  ***
          !----------------------------------------
          call OUTPUT(NN,z1,ropt(1),dz,y,limex_nstep)
          ok = (.not.no_converge)
          if (ok) then      
            if ((T.gt.3000.d0).and.(LL(3).eq.0.d0)) exit
            dz = DMIN1(2.0*dzmerk,dz)            ! Schrittweite vergroessern
            iopt(16) = 0                         ! new initial call
            iopt(10) = 1                         ! do reuse the Jacobian
          endif
        endif
        if (.not.ok) then  
          !----------------------------------------
          !***  Loesung ist nicht ok:           ***
          !***  reduziere interne Schrittweite  ***
          !----------------------------------------
          write(*,*) 'Limex integration failed ',ifail(1),no_converge
          z1 = z1merk
          dz = 0.2*dzmerk
          y(:)  = ymerk(:)
          ys(:) = 0.0
          epsmerk(:) = emerk(:)
          iopt(16) = 0                         ! new initial call
          iopt(10) = 0                         ! no reuse of the Jacobian
          no_converge = .false.
          call FF(NN,z1,y,ys,-2) 
          if (no_converge) then
            write(*,*) "*** no congergence in STATIC_EPS "
            write(*,*) "*** after unsuccessful integration." 
            write(*,*) "===> abort program." 
            stop
          endif
        endif

        !--------------------------------------------------
        !***  regulaeres Ende, wenn Schrittweite klein  ***
        !***  und alle chi_s kleiner als Null           ***
        !--------------------------------------------------
        abbruch = (dz.lt.dzmin)
        do i=1,NDUST
          abbruch = abbruch.and.(chi(i).lt.0.d0)
        enddo
        if (abbruch) then
          write(*,*)
          write(*,*) "regular end of integration dz=",dz 
          write(*,'(a6,99(a10))') " kind=",('  '//dust_nam(:))
          write(*,'(a6,99(1pE10.2))') "  chi=",chi
          write(*,*)
          exit
        endif  
        abbruch = (dz.lt.dzmin).and.(ifail(1).eq.0)
     &            .and.(limex_nstep.eq.0)
        if (abbruch) then
          write(*,*)
          write(*,*) "*** step size too small dz=",dz 
          write(*,'(a6,99(a10))') " kind=",('  '//dust_nam(:))
          write(*,'(a6,99(1pE10.2))') "  chi=",chi
          write(*,*) "*** IS THIS OK? ***"
          write(*,*)
          exit
        endif  
c       if ((dz.le.dzmin).and.(.not.(Jst1.gt.0.d0.and.dz.gt.0.d0))) then
c         write(*,*) "*** step size too small."  
c         write(*,'(99(a10))') " *** kind=",('  '//dust_nam(:))
c         write(*,'(a10,99(1pE10.2))') " ***  chi=",chi
c         !exit
c       endif
c       if (elfunk_calls.gt.500000) exit  
      enddo  

      close(91)
      close(92)
      close(93)
      close(94)
      close(95)
      close(96)
*
      write(*,*) '... finished.'
      write(*,*) 'elfunk_calls= ',elfunk_calls
      write(*,*) 'znull, zende= ',znull,z1
      STOP             
      end
