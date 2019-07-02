!-----------------------------------------------------------------------
! Intiates ggchem input
! ERC 2018
!-----------------------------------------------------------------------
      program init_ggchem_ERC

      implicit real*8 (a-h,o-z)
      character :: s

      real :: tt


      tt=2700
      ptot=516
      open(unit=70,file='marcs2ggchem.in', status='replace')
      write(70, '(19a)') '# selected elements'
      write(70, '(97a)') 'H He C N O Na Mg Si Fe Al Ca Ti S Cl K Li el'
      write(70, '(38a)') '# name of files with molecular kp-data'
      write(70,'(a54)') 'dispol_BarklemCollet.dat             
     &   ! dispol_file '
      write(70,'(a53)')
     & 'dispol_StockKitzmann_withoutTsuji.dat   ! dispol_file2'
      write(70,'(a54)') 'dispol_WoitkeRefit.dat         
     &         ! dispol_file3'
      write(70,'(a64)') '# abundance options 1=EarthCrust, 2=Ocean,
     & 3=Solar, 4=Meteorites'
      write(70,'(a34)') '3                     ! abund_pick'
      write(70,'(a27)') '# equilibrium condensation?'
      write(70,'(a36)') '.false.               ! model_eqcond'
      write(70,'(a15)') '# model options'
      write(70,'(a42)') '1                     ! model_dim  (0,1,2)'
      write(70,'(a36)') '.true.                ! model_pconst'
!      open(unit=10,file='s',status='old')

c Der er noget galt med read delen herunder!! 
      s='2700.'
      read(s,'(a)') TEFF
!      write(70,"(F4.1,a26)") tt, '                ! Tmax [K]'
!      write(70,"(F4.1,a48)") tt-1,'                ! Tmin [K]     
      write(70,*) tt, '                ! Tmax [K]'
      write(70,*) tt-1,'                ! Tmin [K] 

     & (if model_dim>0)'
      write(70,*) ptot,'                   ! pmax [bar]   
     & (if pconst=.true.)'
      write(70,*) ptot,'                   ! pmin [bar]'
      write(70,'(a31)') '100                   ! Npoints'
      write(70,'(a33)') '5                     ! NewBackIt'
      write(70,'(a29)') '1000.0                ! Tfast'
      close(70)
      print *,'Calling ggchem'
      end

!-----------------------------------------------------------------------

