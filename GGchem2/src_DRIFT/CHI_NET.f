*********************************************************************
      REAL*8 FUNCTION CHI_NET(T,Sat,L4,chi,bmix)
*********************************************************************
*****                                                           *****
*****  berechnet die Wachstumegeschwindigkeit durch             *****
*****  SiO- und TiO2-Addition                                   *****
*****                                                           *****
*****  EINGABE:   T = Gastemperatur [K]                         *****
*****           Sat = Uebersaettigungsverhaeltnisse             *****
*****         L4(:) = Vs/Vtot*L4                                *****
*****                                                           *****
*****  AUSGABE: chinet = Wachstumsgeschw [cm s^-1]              *****
*****       bmix(:) = Vs/Vtot                                   *****
*****                                                           *****
*********************************************************************
      use drift_data,ONLY: pi,bk,NSPECIES,NREAC,NDUST,spnr,nat,nmol,
     &                     dust_rho,dust_vol,spmass,pred,
     &                     neduct,reac_nu,reac_sp,keysp
      implicit none      
      real*8,intent(IN) :: T,Sat(NDUST),L4(NDUST)
      real*8,intent(OUT):: chi(NDUST),bmix(NDUST)
      integer :: r,i,j,sp,dustnr
      real*8 :: nsp(NSPECIES),vrel,Sr,sum,chinet,delV,dchi
      real*8 :: eq_fak,rate,rat2,stoi,mr,msum,rsum
      real*8,save :: alpha,Afak 
      logical,save :: firstCall=.true.
      
      if (firstCall) then
        alpha = 1.d0
        Afak  = (36.d0*pi)**(1.d0/3.d0)
        firstCall=.false.
      endif
*
*     --------------------------------------
*     ***  calculate particle densities  ***
*     --------------------------------------
      do i=1,NSPECIES
        j = spnr(i)
        if (j.gt.1000) then
          nsp(i) = nat(j-1000)
        else
          nsp(i) = nmol(j)
        endif
      enddo
*
*     ------------------------
*     ***  calculate bmix  ***
*     ------------------------
      sum  = 0.d0
      do i=1,NDUST
        chi(i) = 0.d0
        sum = sum + L4(i)
      enddo
      do i=1,NDUST
        if (sum.gt.0.d0) then
          bmix(i) = L4(i)/sum
        else
          bmix(i) = 0.d0
          if (i.eq.1) bmix(i)=1.0     ! => rhod=rhod(TiO2) if no dust
        endif
      enddo
 
*     --------------------------------
*     ***  growth and evaporation  ***
*     --------------------------------
      chinet = 0.d0
      do r=1,NREAC
        dustnr = reac_sp(r,1,2)
        rate   = 9.d+99
        eq_fak = 9.d+99
        msum   = 0.d0
        rsum   = 0.d0
        do j=1,neduct(r)
          sp = reac_sp(r,j,1)                       ! index for growth species
          if (keysp(sp)) then
            vrel = DSQRT(bk*T/(2.d0*pi*spmass(sp))) ! rel. thermal velocity
	    stoi = 1.0/DBLE(reac_nu(r,j,1))         ! stoichiom. coeff. of sp
	    rat2 = alpha*nsp(sp)*vrel*stoi          ! surface reac/cm^2/s
	    if (rat2.lt.rate) then
	      rate = rat2                           ! minimum is key 
              mr   = stoi*reac_nu(r,1,2)            ! power for Sr
            endif  
            rsum = rsum + 1.d0/rat2**2
            msum = msum + (stoi*reac_nu(r,1,2))/rat2**2
          endif
        enddo
        mr = msum/rsum
        !Sr = Sat(dustnr)
        Sr = Sat(dustnr)**mr                        ! reaction supersat.ratio
        !eq_fak = 1.d0 - bmix(dustnr)/Sr            ! equilibrium factor 
        !eq_fak = DMAX1(-3.0,eq_fak)                ! de-stiffen evaporation
        eq_fak = 1.d0 - (bmix(dustnr)/Sr)**pred(dustnr)     
        delV   = dust_vol(dustnr)*reac_nu(r,1,2)    ! volume increment by reac
        dchi   = Afak * rate * eq_fak * delV
        chinet = chinet + dchi
        chi(dustnr) = chi(dustnr) + dchi
      enddo
*
      CHI_NET = chinet
      RETURN 
      END 
