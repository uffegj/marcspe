**********************************************************************
      SUBROUTINE INIT_TestCases(beta)
**********************************************************************
      use drift_data,ONLY: NELEM,maxElementCount,maxLayers,
     &                     bk,amu,bar,Nl,abschluss,sizedist,
     &                     eps0,logg,Teff,mixLength,Rnull,
     &                     Rlay,Tlay,play,rholay,glay,mulay,
     &                     wmixlay,zlay,pconv,Nlayers
      implicit none
      real*8,intent(IN) :: beta
      real*8  epsPhoenix(maxElementCount), vconvlay(maxLayers)
      integer Z(maxElementCount)
      integer elementCount,i,j
      real*8  dum,zz,p,T,nH,rho,g,mu,wmix,err,pnorm,Tnorm
      real*8  epsH,Hp,grad,zmax
      logical flag_conv(maxLayers),conv
      integer H,He,C,N,O,Ne,Na,Mg,Al,Si,S,K,Ca,Cr,Mn,Fe,Ni,Ti
*     -----------------------------------------------------------------
      data H/1/, He/2/, C/6/, N/7/, O/8/, Ne/10/, Na/11/, Mg/12/,Al/13/
      data Si/14/, S/16/, K/19/, Ca/20/, Cr/24/, Mn/25/, Fe/26/
      data Ni/28/, Ti/22/
*     -----------------------------------------------------------------

      write(*,*) 
      write(*,*) "reading test case  structure ..."
      write(*,*) "============================="

      open(42, file='TestCase.dat', status='old')

      read(42,*) 
      read(42,*) 
      read(42,*) Teff, logg, mixLength   
      read(42,*) 
      read(42,*) 
      read(42,*) 
      read(42,*) 
      read(42,*) 
      read(42,*) 
      read(42,*) Nlayers
      do i=1, Nlayers
        read (42,*) Tlay(i), play(i) , dum, vconvlay(i), dum,
     &              zlay(i), rholay(i), mulay(i)
        glay(i) = 10.d0**logg
      enddo
      if (ABS(Teff-600.0)<1.0) then
        zlay = zlay*1.e+5   ! Korrektur z[km]->z[cm]
        zmax = zlay(1)
        zlay = zmax-zlay    ! umdrehen
      endif  
      print*, "Nlayers= ", Nlayers
      close(42)
*
*     -----------------------------------
*     ***  calculate wmix from vconv  ***
*     -----------------------------------
      conv = .true.
      do i=Nlayers,1,-1
        if (conv.and.vconvlay(i).gt.0.d0) then
          Hp = bk*Tlay(i)/(glay(i)*mulay(i)*amu)
          wmixlay(i) = vconvlay(i)/(mixLength*Hp)  
        else 
          if (conv) then
            pconv = play(i)
            write(*,'(" convection zone ends at p[bar]=",0pF6.2,
     &                " (T=",0pF8.2,"K)")') play(i)/bar,Tlay(i)
          endif  
          wmixlay(i) = 1.d-99
          conv = .false.
        endif    
c       write(*,*)  play(i)/bar,conv,wmixlay(i)
      enddo
*     ------------------------------------------
*     ***  limit log(wmix)-gradient to beta  ***     
*     ------------------------------------------
      do i=Nlayers-1,1,-1
        grad = DLOG(wmixlay(i)/wmixlay(i+1))/DLOG(play(i)/play(i+1))
        if (grad.gt.beta) then
          wmixlay(i) = wmixlay(i+1)*DEXP(beta*DLOG(play(i)/play(i+1)))
        endif    
        grad = DLOG(wmixlay(i)/wmixlay(i+1))/DLOG(play(i)/play(i+1))
        write(*,1100) play(i)/bar,Tlay(i),1.d0/wmixlay(i),grad
      enddo
*
*     -----------------------
*     ***  init epsilons  ***
*     -----------------------
      eps0(:)  = 1.d-99
      eps0(H)  = 10.d0**(12.00-12.0)
      eps0(He) = 10.d0**(10.99-12.0) 
      eps0(C)  = 10.d0**( 8.55-12.0)
      eps0(N)  = 10.d0**( 7.97-12.0)
      eps0(O)  = 10.d0**( 8.87-12.0)
      eps0(Na) = 10.d0**( 6.33-12.0)
      eps0(Mg) = 10.d0**( 7.58-12.0)
      eps0(Al) = 10.d0**( 6.47-12.0)
      eps0(Si) = 10.d0**( 7.55-12.0)
      eps0(S)  = 10.d0**( 7.21-12.0)
      eps0(K)  = 10.d0**( 5.12-12.0)
      eps0(Ca) = 10.d0**( 6.36-12.0)
      eps0(Ti) = 10.d0**( 5.02-12.0)
      eps0(Fe) = 10.d0**( 7.50-12.0)
*
*     ------------------------------
*     ***  check THERMO-routine  ***
*     ------------------------------
      err = 0.d0
      Rnull = 0.d0
      do i=3,Nlayers
        zz = 0.1*zlay(i)+0.9*zlay(i-1)  
        call THERMO(zz,p,T,nH,rho,g,mu,wmix,3)
c       write(*,*) play(i-1),p,play(i)
c       write(*,*) (DLOG(p)-DLOG(play(i-1)))/(rr-Rlay(i-1)),
c    &             (DLOG(play(i))-DLOG(p))/(Rlay(i)-rr)
        pnorm = DEXP(0.1d0*DLOG(play(i))+0.9d0*DLOG(play(i-1)))
        Tnorm = 0.1d0*Tlay(i)+0.9d0*Tlay(i-1)
        err = DMAX1(err,DABS(pnorm/p-1.d0),DABS(Tnorm/T-1.d0)) 
      enddo    
      write(*,*) "maximum interpolation error=",err
      write(*,*) "finished reading TestCase structure."
      if (err.gt.1.d-10) stop "*** too large interpolation errors!"
*
      RETURN
    3 FORMAT(f10.3,3(2X, 1E11.3),2x,D11.3,2x,D12.5,D12.5,D12.5)
  100 format (f12.3, 1x, f12.3, 1x, f12.3)
  200 format (i5)
  210 format (8(i5, 1x))
  300 format (8(e15.8, 1x)) 
 1000 format(' Teff=',0pF8.3,' logg=',0pF5.2,' mixLengthPara=',0pF5.2)
 1100 format (8(1pe15.8, 1x)) 
      end
