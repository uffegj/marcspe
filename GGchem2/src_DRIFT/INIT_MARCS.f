**********************************************************************
      SUBROUTINE INIT_MARCS(beta)
**********************************************************************
      use drift_data,ONLY: NELEM,maxElementCount,maxLayers,
     &                     bk,amu,bar,Nl,abschluss,sizedist,
     &                     eps0,logg,Teff,mixLength,Rnull,
     &                     Rlay,Tlay,play,rholay,glay,mulay,
     &                     wmixlay,zlay,pconv,Nlayers
      implicit none
      real*8,intent(IN) :: beta
      real*8  epsMarcs(maxElementCount), vconvlay(maxLayers)
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
      write(*,*) "reading MARCS structure ..."
      write(*,*) "============================="

      open(42, file='marcs2drift.dat', status='old')
      do i=1,5 
         read(42,*) 
      enddo
      read(42,100) Teff, logg, mixLength   
      write(*,1000) Teff, logg, mixLength
      read(42,*)
      read(42,*)
      read(42,200) Nlayers
      write(*,*) "Nlayers", Nlayers
      read(42,*)
      read(42,*)
      read(42,200) elementCount
      read(42,*)
      read(42,*)
      read(42,210) (Z(i), i=1,  elementCount)
      read(42,*)
      read(42,*)
      do i=1,Nlayers
        read(42,'(5x,6e16.8,a6,e16.8)') Rlay(i),Tlay(i), 
     &    play(i),rholay(i),glay(i),vconvlay(i),flag_conv(i),mulay(i)
      end do
      read(42,*)
      read(42,*)
      read(42,'(8f6.2)') (epsMarcs(i), i=1,elementCount)
      close(42)

      eps0(:)  = 1.d-99
      eps0(H)  = 10.d0**(epsMarcs(1)-12.0)
      eps0(He) = 10.d0**(epsMarcs(2)-12.0)
      eps0(C)  = 10.d0**(epsMarcs(6)-12.0)
      eps0(N)  = 10.d0**(epsMarcs(7)-12.0)
      eps0(O)  = 10.d0**(epsMarcs(8)-12.0)
      eps0(Na) = 10.d0**(epsMarcs(11)-12.0)
      eps0(Mg) = 10.d0**(epsMarcs(12)-12.0)
      eps0(Al) = 10.d0**(epsMarcs(13)-12.0)
      eps0(Si) = 10.d0**(epsMarcs(14)-12.0)
      eps0(S)  = 10.d0**(epsMarcs(16)-12.0)
      eps0(K)  = 10.d0**(epsMarcs(19)-12.0)
      eps0(Ca) = 10.d0**(epsMarcs(20)-12.0)
      eps0(Ti) = 10.d0**(epsMarcs(21)-12.0)
      eps0(Fe) = 10.d0**(epsMarcs(25)-12.0)

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
        pnorm = DEXP(0.1d0*DLOG(play(i))+0.9d0*DLOG(play(i-1)))
        Tnorm = 0.1d0*Tlay(i)+0.9d0*Tlay(i-1)
        err = DMAX1(err,DABS(pnorm/p-1.d0),DABS(Tnorm/T-1.d0))          
      enddo    
      write(*,*) "maximum interpolation error=",err
      write(*,*) "finished reading MARCS structure."
      if (err.gt.1.d-10) stop "*** too large interpolation errors!"

      RETURN
  100 format (3(f12.3, 1x))
  200 format (i5)
  210 format (8(i5, 1x))
  300 format (8(e15.8, 1x)) 
 1000 format(' Teff=',0pF8.3,' logg=',0pF5.2,' mixLengthPara=',0pF5.2)
 1100 format (8(1pe15.8, 1x)) 
      end
