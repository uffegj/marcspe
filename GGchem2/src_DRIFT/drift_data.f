      module drift_data

      implicit none
      !--------------------------
      !***  nature constants  ***
      !--------------------------
      real*8,parameter :: bk=1.380622D-16, rgas=8.31434D+7 
      real*8,parameter :: hplanck=6.626196D-27, me=9.109558D-28 
      real*8,parameter :: grav=6.6732D-08, elad=4.80325D-10
      real*8,parameter :: sigg=5.66961D-5, cl=2.997925D+10
      real*8,parameter :: pi=3.141592653589793D0
      !---------------------
      !***  other units  ***
      !---------------------
      real*8,parameter :: amu=1.660531D-24, Msun=1.989D+33
      real*8,parameter :: Lsun=3.826D+33, Rsun=6.9599D+10
      real*8,parameter :: day=8.64D+4, yr=3.15576D+7
      real*8,parameter :: pc=3.085D+18, bar=1.D+6, joule=1.D+7
      real*8,parameter :: km=1.D+5, eV=1.602192D-12, AA=1.D-8

      !--------------------
      !***  dimensions  ***
      !--------------------
!      integer,parameter :: NELEM=36, NMOLE=167, NMOM=4, NEPS=8  ! dust12
      integer,parameter :: NELEM=36, NMOLE=167, NMOM=4, NEPS=6  ! dust7
      integer,parameter :: NPOINTS=500
!      integer,parameter :: NSPECIES=31, NDUST=12, NNUC=1, NREAC=60  ! dust12
      integer,parameter :: NSPECIES=25, NDUST=7, NNUC=1, NREAC=32  ! dust7
      integer,parameter :: NN=NMOM+NDUST
      integer,parameter :: maxElementCount=200, maxLayers=512

      !------------------------------------
      ! ***  constant global variables  ***
      !------------------------------------
      character(len=2),save  :: elnam(NELEM) 
      character(len=10),save :: spnam(NSPECIES),dust_nam(NDUST)
      character(len=10),save :: cmol(NMOLE)
      character(len=8),save  :: nuc_nam(NNUC)
      integer,save :: welche,abschluss,sizedist,Nl,Nlayers
      integer,save :: elnr(NEPS),elcode(NELEM),spnr(NSPECIES)
      integer,save :: dust_nel(NDUST),dust_nu(NDUST,5),dust_el(NDUST,5)
      integer,save :: nuc_nel(NNUC),nuc_nu(NNUC,5),nuc_el(NNUC,5)
      integer,save :: neduct(NREAC),nprod(NREAC)
      integer,save :: reac_nu(NREAC,5,2),reac_sp(NREAC,5,2)
      logical,save :: keysp(NSPECIES)
      real*8,save :: eps0(NELEM),logg,Teff,mixLength,beta,Rnull,Vl
      real*8,save :: Rlay(maxLayers),Tlay(maxLayers),play(maxLayers)
      real*8,save :: rholay(maxLayers),glay(maxLayers),mulay(maxLayers)
      real*8,save :: wmixlay(maxLayers),zlay(maxLayers),pconv
      real*8,save :: cref,aref,Vref,Lref(0:NMOM),pref
      real*8,save :: spmass(NSPECIES),dust_rho(NDUST),dust_vol(NDUST)
      real*8,save :: mass(NELEM),muH
      real*8,save :: pred(NDUST)

      !--------------------------------------------------
      ! ***  physical conditions at a certain height  ***
      !--------------------------------------------------
      real*8,save :: eps(NELEM),T,rho,nH,p,pel,wmix,mu,g,rhoL2
      real*8,save :: nel,nat(NELEM),nion(NELEM),nmol(NMOLE)
      real*8,save :: Sat(NDUST),Nst,Jst(NNUC),Jev(NNUC)
      real*8,save :: chinet,chi(NDUST),nsp(NSPECIES)
      real*8,save :: LL(0:NMOM),L4(1:NDUST),bmix(NDUST),rhod
      real*8,save :: imp_nuklea,imp_growth,imp_drift,imp_misch

      !----------------------------
      ! ***  control variables  ***
      !----------------------------      
      logical,save :: restart,PHOENIX_MODE,firstStaticEpsCall
      real*8,save  :: epsmerk(NEPS)
      logical,save :: no_converge,J_is_zero
      integer*8,save :: elfunk_calls

      end module
