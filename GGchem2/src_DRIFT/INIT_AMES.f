**********************************************************************
      SUBROUTINE INIT_AMES(beta)
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
      real*8  epsH,Hp,grad, C_to_O
      logical flag_conv(maxLayers),conv

      write(*,*) 
      write(*,*) "reading AMES structure ..."
      write(*,*) "============================="

      open(43, file='AMES.data', status='old')
      do i=1,5 
         read(43,*) 
      enddo
      read(43,*) Teff, logg, mixLength   
      write(*,1000) Teff, logg, mixLength
      read(43,*)
      read(43,*)
      read(43,*) Nlayers
      read(43,*)
      read(43,*)
      read(43,*) elementCount
      read(43,*)
      read(43,*)
      read(43,*) (Z(i), i=1, elementCount)
c      write(*,*) Nlayers, elementCount, (Z(i), i=1, elementCount)
      read(43,*)
      read(43,*)
      read(43,*)
      do i=1,Nlayers
        read(43,3) j, dum, Tlay(i), play(i), dum, rholay(i), mulay(i),
     &                Rlay(i), zlay(i), dum
        glay(i) = 10.d0**logg
      enddo
      read(43,*)
      read(43,*)
      do i=1,Nlayers
         read(43,4) j, dum,dum,dum,dum,dum,dum,dum,dum, vconvlay(i),
     &                 dum,dum,dum,dum, dum !flag_conv(i)
      enddo 
      read(43,*)
      read(43,*)
      read(43,*) (epsPhoenix(i), i=1,elementCount)
      epsH = epsPhoenix(1)
      do i=1,elementCount
        if ((Z(i).ge.1).and.(Z(i).le.NELEM)) then
          eps0(Z(i)) = epsPhoenix(i)/epsH
        endif  
      enddo
*     ------------------
*     change C/O ration:
*     ------------------
*     use these line for changing C/O, else leave it as it is:
*     --------------------------------------------------------
c      C_to_O = 1.5
c      eps0(6)= C_to_O * eps0(8)
      write(*,*) "C/O =",   eps0(6)/eps0(8) 
      close(43)
 
*
*     -----------------------------------
*     ***  calculate wmix from vconv  ***
*     -----------------------------------
      conv = .true.
      do i=Nlayers,1,-1
        if (conv.and.vconvlay(i).gt.0.d0) then
          Hp = bk*Tlay(i)/(glay(i)*mulay(i)*amu)
c          write(*,*) "vconvlay(i), Hp", vconvlay(i), Hp
          wmixlay(i) = vconvlay(i)/(mixLength*Hp)  
        else 
          if (conv) then
            pconv = play(i)
            write(*,*) "convection zone ends at p[bar]=",play(i)/bar
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
c       write(*,*) play(i)/bar,wmixlay(i),grad
        if (grad.gt.beta) then
          wmixlay(i) = wmixlay(i+1)*DEXP(beta*DLOG(play(i)/play(i+1)))
        endif    
      enddo    
*
*     ------------------------------
*     ***  check THERMO-routine  ***
*     ------------------------------
      err = 0.d0
      Rnull = 0.0
      do i=3,Nlayers
        zz = 0.1*zlay(i)+0.9*zlay(i-1)  
        call THERMO(zz,p,T,nH,rho,g,mu,wmix,2)
c       write(*,*) play(i-1),p,play(i)
c       write(*,*) (DLOG(p)-DLOG(play(i-1)))/(rr-Rlay(i-1)),
c    &             (DLOG(play(i))-DLOG(p))/(Rlay(i)-rr)
        pnorm = DEXP(0.1*DLOG(play(i))+0.9*DLOG(play(i-1)))
        Tnorm = 0.1*Tlay(i)+0.9*Tlay(i-1)
        err = DMAX1(err,DABS(pnorm/p-1.d0),DABS(Tnorm/T-1.d0))          
      enddo    
      write(*,*) "maximum interpolation error=",err
      write(*,*) "finished reading AMES structure."
      if (err.gt.1.d-10) stop "*** too large interpolatiojn errors!"
*
      RETURN
 1000 format(' Teff=',0pF8.2,' logg=',0pF5.2,' mixLengthPara=',0pF5.2)
    3 FORMAT(2X, I2, 2X, 1D10.3, 2x, f10.4, 3(2X, 1D10.3), 2x, f6.3,
     &3(2X, 1PE10.3))
    4 FORMAT(3X, I2, 13(1X, 1D8.3), 1x, I2)
      end
