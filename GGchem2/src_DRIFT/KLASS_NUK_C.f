!*********************************************************************
      SUBROUTINE KLASS_NUK_C(T,nC1,nC2,nC2H,nC2H2,nC3,Jst,Nst,SS)
!*********************************************************************
!*****                                                           *****
!*****  berechnet die Kohlenstoff-Keimbildungsrate               *****
!*****  mit Hilfe von klassischer Nukleationstheorie.            *****
!*****                                                           *****
!*****  EINGABE:   T = Gastemperatur [K]                         *****
!*****           nC1 = Teilchendichte des Monomers C1 [cm^-3]    *****
!*****           nC2 = Teilchendichte des Dimers   C2 [cm^-3]    *****
!*****           nC3 = Teilchendichte des Trimers  C3 [cm^-3]    *****
!*****                                                           *****
!*****  AUSGABE: Jst = Keimbildungsrate [cm^-3 s^-1]             *****
!*****           Nst = Groesse des kritischen Clusters           *****
!*****            SS = Uebersaettigungsverhaeltnis               *****
!*****                                                           *****
!*********************************************************************
       implicit none
       real*8 :: T,nC1,nC2,nC2H,nC2H2,nC3
       real*8 :: Jst,Nst
       real*8 :: pi,bk,amu,f0,molg,sigma,Nl,alf1,alf2,alf3,a0,Tdispol
       real*8,parameter   :: mH=1.008, mC=12.01
       real*8,parameter   :: mC1=mC, mC2=2.0*mC, mC3=3.0*mC
       real*8,parameter   :: mC2H=2.0*mC+mH, mC2H2=2.0*mC+2.0*mH
c      !--------------------------------------------------------------
c      !***  Materialkonstanten von Graphit:                       ***
c      !***  sigma = Oberflaechenspannung in [erg/cm2]             ***
c      !***  Nl = Clustergroesse, bei der sich Theta halbiert      ***
c      !***  a0 = Monomerradius in cm                              ***
c      !***  f0 = 4*pi*a0^2 Monomeroberflaeche in [cm2]            ***
c      !***  alf1,alf2,alf3 = Sticking-coefficients fuer C1,C2,C3  ***
c      !***    ((((( nach GAIL:   sigma=1400   Nl=5 )))))          ***
c      !-------------------------------------------------------------- 
c      !real,parameter :: sigma=2500.e0, Nl=30.e0

       data pi/3.14159265358979D+0/, bk/1.38066D-16/, amu/1.66055D-24/
       data sigma/1400.d0/,Nl/5.d0/,a0/1.28e-8/
       data f0/2.07e-15/,molg/79.898d0/
       data alf1/0.37e0/, alf2/0.34e0/, alf3/0.08e0/
       data Tdispol/300.d0/
       save


       real*8 :: Tvar,SS,th,ex,psat,slog,xx
       real*8 :: thetun,thetaN,x0,x1,x2,x3,dgdn,fst
       real*8 :: zeldof,vth,beta,betaN,Ast

c      !---------------------------------------------------------
c      !***  Dampfdruck von C1 (Gas) ueber Graphit (fest)     ***
c      !***  pvap:  eigener Fit an Janaf(1985)                ***
c      !---------------------------------------------------------
       Tvar = MAX(Tdispol,T)

c Gail Dampfdruck
c       xx   = 5040.0/Tvar
c       ex   = 3.24595E+01 - 1.68624E+01*xx    - 5.17205E-02*xx**2 
c     &                         + 3.99686E-03*xx**3 - 1.00638E-04*xx**4
c       psat = EXP(ex)
c ???????????????
       ex   = +1.01428e+06/Tvar -7.23043e+05 +1.63039e+02*Tvar 
     &        -1.75890e-03*Tvar**2 +9.97416e-08*Tvar**3
       psat = 1.e+6*EXP(ex/(8.31441e0*Tvar))
c
       SS   = nC1*bk*Tvar/psat

       slog = LOG(SS)
       if (SS.le.1.e0) then
         Jst = 0.e+0
         Nst = 0.e+0
         goto 500
       endif

      !---------------------------------------------------------------
      !***  Groesse des krit. Clusters nach dem Troepfchenmodel    ***
      !---------------------------------------------------------------
       thetun = f0*sigma/bk
       x0     = 2.e0*thetun/(3.e0*Tvar*slog)
       x1     = Nl**(1.e0/3.e0)
       x2     = x1/x0
       x3     = 0.5e0*(1.e0+SQRT(1.e0+2.e0*x2)) - x2
       Nst    = 1.e0 + (x0*x3)**3 
       if (Nst.le.1.e0) Nst=1.000000001e0

      !------------------------------------------------
      !***  Teilchenhaeufigkeit des krit. Clusters  ***
      !------------------------------------------------
       x0     = x0*slog
       x2     = (Nst-1.e0)**(1.e0/3.e0) + x1
       x3     = 1.e0/(Nst-1.e0)**(2.e0/3.e0)
       dgdn   = x0*x3* ( 1.e0/x2**2 + x1/x2**3 ) / 3.e0
       zeldof = SQRT(dgdn/(2.e0*pi))
       thetaN = thetun/(1.e0+(Nl/(Nst-1.e0))**(1.e0/3.e0))
       x1     = (Nst-1.e0)*slog - (thetaN/Tvar)*(Nst-1.e0)**(2.e0/3.e0)
       if (x1.lt.-300.0) then
         Jst = 0.e+0
         Nst = 0.e+0
         goto 500
       endif
       fst = nc1*EXP(x1)

       !----------------------------------
       !***  Wachstumsgeschwindigkeit  ***
       !----------------------------------
        vth  = SQRT(bk*T/(2.e0*pi*amu))
c CHILD_NUCL
c        beta = vth*( 1.0*alf1*nC1/SQRT(mC1)              
c     &       + 2.0*alf2*SQRT(2.0)*2.0              
c     &       * ( nC2  /SQRT(mC2)                   
c     &       + nC2H /SQRT(mC2H)                  
c     &       + nC2H2/SQRT(mC2H2) ) )
c ?????????
        beta = vth*( 1.0*alf1*nC1  /SQRT(mC1)            
     &       +2.0*alf2*nC2  /SQRT(mC2)            
     &       +2.0*alf2*nC2H /SQRT(mC2H)           
     &       +2.0*alf2*nC2H2/SQRT(mC2H2)         
     &       +3.0*alf3*nC3  /SQRT(mC3) )

        Ast = f0*Nst**(2.e0/3.e0)

       !--------------------------
       !***  Keimbildungsrate  ***
       !--------------------------
       Jst = fst*beta*Ast*zeldof

 500    continue
       RETURN
       end
