**********************************************************************
      SUBROUTINE INIT_PHOENIX(beta)
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
      real*8  dum,rr,p,T,nH,rho,g,mu,wmix,err,pnorm,Tnorm
      real*8  epsH,Hp,grad
      logical flag_conv(maxLayers),conv
      integer H,He,C,N,O,Ne,Na,Mg,Al,Si,S,K,Ca,Cr,Mn,Fe,Ni,Ti
*     -----------------------------------------------------------------
      data H/1/, He/2/, C/6/, N/7/, O/8/, Ne/10/, Na/11/, Mg/12/,Al/13/
      data Si/14/, S/16/, K/19/, Ca/20/, Cr/24/, Mn/25/, Fe/26/
      data Ni/28/, Ti/22/
*     -----------------------------------------------------------------


      write(*,*) 
      write(*,*) "reading PHOENIX structure ..."
      write(*,*) "============================="

      open(42, file='2Drift.data', status='old')
c     write(*,*) "geoeffnet"
      do i=1,5 
         read(42,*) 
c        write(*,*) i, "leere Zeile"
      enddo
c     write(*,*) "weiter"      
      read(42,100) Teff, logg, mixLength   
      write(*,1000) Teff, logg, mixLength
      read(42,*)
      read(42,*)
      read(42,200) Nlayers
      write(*,*) "Nlayers", Nlayers
      read(42,*)
      read(42,*)
      read(42,200) elementCount
c     write(*,*) "elementCount", element)*kT)Count
      read(42,*)
      read(42,*)
      read(42,210) (Z(i), i=1,  elementCount)
c     write(*,*) (Z(i), i=1,  elementCount)
      read(42,*)
      read(42,*)
      read(42,300) (Rlay(i),i=1,Nlayers)
c     write(*,*) "nach Rlay"
      read(42,*)
      read(42,*)
      read(42,300) (Tlay(i),i=1,Nlayers)
c     write(*,*) "nach Tlay"
      read(42,*)
      read(42,*)
      read(42,300) (play(i),i=1,Nlayers)
c     write(*,*) "nach play"
      read(42,*)
      read(42,*)
      read(42,300) (rholay(i),i=1,Nlayers)
c     write(*,*) "nach rholay"
      read(42,*)
      read(42,*)
      read(42,300) (glay(i),i=1,Nlayers)
c     write(*,*) "nach glay"
      read(42,*)
      read(42,*)
      read(42,300) (vconvlay(i),i=1,Nlayers)
c     write(*,*) "nach vconvlay"
      read(42,*)
      read(42,*)
      read(42,'(8(l2, 1x))') (flag_conv(i),i=1,Nlayers)
      read(42,*)
      read(42,*)
      read(42,300) (mulay(i),i=1,Nlayers)
      read(42,*)
      read(42,*)
c      read(42,300) (epsPhoenix(i), i=1,elementCount)
c      epsH = epsPhoenix(1)
c      do i=1,elementCount
c        if ((Z(i).ge.1).and.(Z(i).le.NELEM)) then
c          eps0(Z(i)) = epsPhoenix(i)/epsH
cc          eps0(Z(4)) = epsPhoenix(6)*1.5            ! 6 = O, 4 = C
c        endif  
c      enddo
      close(42)


*
*     -----------------------
*     ***  init epsilons  ***
*     -----------------------
      eps0(:)  = 1.d-99
      eps0(H)  = 10.d0**(12.00-12.0)   ! 12.0
      eps0(He) = 10.d0**(10.93-12.0)   ! 10.99
      eps0(C)  = 10.d0**( 8.38-12.0)   ! 8.55
      eps0(N)  = 10.d0**( 7.77-12.0)   ! 7.97
!      eps0(O)  = 10.d0**( 8.65-12.0)   ! 8.87
      eps0(O)  = 4.5693E-04   ! 10.d0**( 8.62-12.0)   ! 8.66
      eps0(Na) = 10.d0**( 6.15-12.0)   ! 6.33
!      eps0(Mg) = 10.d0**( 7.52-12.0)   ! 7.58
      eps0(Mg) = 3.3884E-05   !10.d0**( 7.53-12.0)
!      eps0(Al) = 10.d0**( 6.35-12.0)
      eps0(Al) = 2.3442E-06 !10.d0**( 6.37-12.0)   ! 6.47
!      eps0(Si) = 10.d0**( 7.50-12.0)   ! 7.55
      eps0(Si) = 3.2359E-05   !10.d0**( 7.51-12.0)
      eps0(S)  = 10.d0**( 7.13-12.0)   ! 7.21
      eps0(K)  = 10.d0**( 5.06-12.0)   ! 5.12
      eps0(Ca) = 10.d0**( 6.29-12.0)   ! 6.36
      eps0(Ti) = 10.d0**( 4.88-12.0)   ! 5.02
!      eps0(Fe) = 10.d0**( 7.44-12.0)   ! 7.50
      eps0(Fe) = 2.8184E-05 !10.d0**( 7.45-12.0)
*
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
c       write(*,*) play(i)/bar,conv,wmixlay(i)
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
c       write(*,1100) play(i)/bar,Tlay(i),1.d0/wmixlay(i),grad
      enddo
*
*     ------------------------------
*     ***  check THERMO-routine  ***
*     ------------------------------
      err = 0.d0
      Rnull = Rlay(2)   
      do i=3,Nlayers
        rr = 0.1d0*Rlay(i)+0.9d0*Rlay(i-1)  
        call THERMO(Rnull-rr,p,T,nH,rho,g,mu,wmix,1)
c       write(*,*) play(i-1),p,play(i)
c       write(*,*) (DLOG(p)-DLOG(play(i-1)))/(rr-Rlay(i-1)),
c    &             (DLOG(play(i))-DLOG(p))/(Rlay(i)-rr)
        pnorm = DEXP(0.1d0*DLOG(play(i))+0.9d0*DLOG(play(i-1)))
        Tnorm = 0.1d0*Tlay(i)+0.9d0*Tlay(i-1)
        err = DMAX1(err,DABS(pnorm/p-1.d0),DABS(Tnorm/T-1.d0))          
      enddo    
      write(*,*) "maximum interpolation error=",err
      write(*,*) "finished reading PHOENIX structure."
      if (err.gt.1.d-10) stop "*** too large interpolation errors!"
*
      RETURN
  100 format (f12.3, 1x, f12.3, 1x, f12.3)
  200 format (i5)
  210 format (8(i5, 1x))
  300 format (8(e15.8, 1x)) 
 1000 format(' Teff=',0pF8.3,' logg=',0pF5.2,' mixLengthPara=',0pF5.2)
 1100 format (8(1pe15.8, 1x)) 
      end
