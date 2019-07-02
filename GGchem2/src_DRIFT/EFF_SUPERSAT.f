*********************************************************************
      SUBROUTINE EFF_SUPERSAT(nH,T,bmix,effSat)
*********************************************************************
*****                                                           *****
*****  EINGABE:  nH = total hydrogen density [cm^-3]            *****
*****             T = Gastemperatur [K]                         *****
*****       bmix(:) = Vs/Vtot                                   *****
*****                                                           *****
*****  AUSGABE: effSat = effektive Uebersaettigungsverhaelt.    *****
*****                                                           *****
*********************************************************************
      use drift_data,ONLY: pi,bk,NSPECIES,NREAC,NDUST,spnr,nat,nmol,
     &                     dust_rho,dust_vol,spmass,
     &                     neduct,reac_nu,reac_sp,keysp,pred
      implicit none      
      real*8,intent(IN) :: nH,T,bmix(NDUST)
      real*8,intent(OUT):: effSat(NDUST)
      integer :: r,i,j,sp,dustnr
      real*8 :: Sat(NDUST),sum1(NDUST),sum2(NDUST)
      real*8 :: nsp(NSPECIES),vrel,Sr,sum,delV
      real*8 :: eq_fak,rate,rat2,stoi,mr,msum,rsum,dum
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
      call GGCHEM(nH,T,dum)
      do i=1,NSPECIES
        j = spnr(i)
        if (j.gt.1000) then
          nsp(i) = nat(j-1000)
        else
          nsp(i) = nmol(j)
        endif
      enddo
*
*     ----------------------------------------------
*     ***  calculate pure supersaturation ratio  ***
*     ----------------------------------------------
      call SUPERSAT(nH,T,Sat)

*     --------------------------------
*     ***  growth and evaporation  ***
*     --------------------------------
      sum1 = 0.d0
      sum2 = 0.d0
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
        Sr = Sat(dustnr)**mr                        ! reaction supersat.ratio
        !eq_fak = 1.d0 - bmix(dustnr)/Sr            ! equilibrium factor 
        !eq_fak = DMAX1(-100.0,eq_fak)              ! de-stiffen evaporation
        eq_fak = 1.d0 - (bmix(dustnr)/Sr)**pred(dustnr) 
        delV   = dust_vol(dustnr)*reac_nu(r,1,2)    ! volume increment by reac
        sum1(dustnr) = sum1(dustnr) + delV * rate * eq_fak
        sum2(dustnr) = sum2(dustnr) + delV * rate
      enddo
*
*     ----------------------------------------------------
*     ***  calculate effective supersaturation ratios  ***
*     ----------------------------------------------------
      do i=1,NDUST
        if (sum2(i)-sum1(i).ne.0.d0) then 
c         effSat(i) = 1.d0/(1.d0-sum1(i)/sum2(i))
          effSat(i) = sum2(i)/(sum2(i)-sum1(i))
        else
          effSat(i) = Sat(i)
        endif  
      enddo

      RETURN 
      END 
