*********************************************************************
      SUBROUTINE KLASS_NUK_TIO2(T,nTiO2,Jstern,Nstern,SS)
*********************************************************************
*****                                                           *****
*****  berechnet die TiO2-Keimbildungsrate                      *****
*****  mit Hilfe von klassischer Nukleationstheorie.            *****
*****                                                           *****
*****  INPUT:   T = Gastemperatur [K]                           *****
*****                                                           *****
*****  OUTPUT: Jstern = Keimbildungsrate [cm^-3 s^-1]           *****
*****                   nucleation rate
*****          Nstern = Groesse des kritischen Clusters         *****
*****                   critical cluster size 
*****              SS = Uebersaettigungsverhaeltnis             *****
*****                   supersaturation ratio 
*****                                                           *****
*********************************************************************
      implicit none
      real*8   T,nTiO2
      real*8   Jstern,Nstern,SS
      real*8   pi,bk,amu,f0,molg,sigma,Nf,alfa
      real*8   ex,psat,slog
      real*8   thetun,thetaN,x0,x1,x2,x3,dgdn,nst
      real*8   zeldof,vth,beta,fNst
      data pi/3.14159265358979D+0/, bk/1.38066D-16/, amu/1.66055D-24/
*
*     --------------------------------------------------------------
*     ***  Materialkonstanten von TiO2(fest):                    ***
*     ***  sigma = Oberflaechenspannung in erg/cm^2              ***
*                = 620   erg/cm^2 Jeong
*                = 480.6 erg/cm^2 Graham Lee
*     ***  Nf = Clustergroesse, bei der sich Theta halbiert      ***
*     ***  f0 = 4*pi*a0^2 Monomeroberflaeche in cm^2             ***
*     ***  alfa = Sticking-coefficients                          ***
*     --------------------------------------------------------------
!      data sigma/480.d0/,Nf/0.d0/
      data sigma/620.d0/,Nf/0.d0/
      data f0/4.808362d-15/,molg/79.898d0/,alfa/1.d0/
      save
*
*     ---------------------------------------------------
*     ***  Dampfdruck von TiO2 (Gas) ueber TiO2(fest) ***
*     ***  nach JANAF-Tafeln, elektr. Version 1985    ***
*     ***  Fit zwischen 500K und 2500K  (Woitke)      ***
*     ---------------------------------------------------
      psat = DEXP(35.8027-74734.7/T)
      SS   = nTiO2*bk*T/psat
      slog = DLOG(SS)
      if (SS.le.1.d0) then
        Jstern = 0.d+0
        Nstern = 9.d+99
        goto 500
      end if  
*
*     ---------------------------------------------------------------
*     ***  Groesse des krit. Clusters nach dem Troepfchenmodel    ***
*     ***  Size of critical cluster according to droplet model    ***
*     ---------------------------------------------------------------
      thetun = f0*sigma/bk
      x0     = 2.d0*thetun/(3.d0*T*slog)
      x1     = Nf**(1.d0/3.d0)
      x2     = x1/x0
      x3     = 0.5d0*(1.d0+DSQRT(1.d0+2.d0*x2)) - x2
      Nstern = 1.d0 + (x0*x3)**3 
      if (Nstern.le.1.d0) Nstern=1.000000001D0
*
*     ------------------------------------------------
*     ***  Teilchenhaeufigkeit des krit. Clusters  ***
*     ***  Number density of critical cluster      ***
*     ------------------------------------------------
      x0     = x0*slog
      x2     = (Nstern-1.d0)**(1.d0/3.d0) + x1
      x3     = 1.d0/(Nstern-1.d0)**(2.d0/3.d0)
      dgdn   = x0*x3* ( 1.d0/x2**2 + x1/x2**3 ) / 3.d0
      zeldof = DSQRT(dgdn/(2.d0*pi))
      thetaN = thetun/(1.d0+(Nf/(Nstern-1.d0))**(1.d0/3.d0))
      x1     = (Nstern-1.d0)*slog - (thetaN/T)
     &         *(Nstern-1.d0)**(2.d0/3.d0)
      nst    = nTiO2*DEXP(x1)
*
*     ----------------------------------
*     ***  Wachstumsgeschwindigkeit  ***
*     ***  Growth velocity           ***
*     ----------------------------------
      vth  = DSQRT(bk*T/(2.d0*pi*molg*amu))
      beta = vth*alfa*nTiO2
      fNst = f0*Nstern**(2.d0/3.d0)
*
*     --------------------------
*     ***  Keimbildungsrate  ***
*     ***  Nucleation rate   ***
*     --------------------------
      Jstern = nst*beta*fNst*zeldof

c!!!!!!!!!!!!!!!!!!!!!!!!
c     test response of model to J*
c!!!!!!!!!!!!!!!!!!!!!!!!

c      Jstern = Jstern * 1.d-04

*
 500  continue
      RETURN
      end
