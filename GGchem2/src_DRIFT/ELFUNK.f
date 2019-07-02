**********************************************************************
      SUBROUTINE ELFUNK(epsi,jj,deps,FUNC)
**********************************************************************
      use drift_data,ONLY: NELEM,NEPS,NDUST,NNUC,NREAC,NSPECIES,
     &                     pi,bk,spmass,elfunk_calls,
     &                     eps,eps0,elnr,spnr,nat,nmol,
     &                     dust_nel,dust_el,dust_nu,dust_rho,dust_vol,
     &                     elcode,nuc_nel,nuc_el,nuc_nu,Nl,
     &                     neduct,nprod,reac_nu,reac_sp,keysp,
     &                     T,nH,rhoL2,wmix,L4,pred
      implicit none
      integer,intent(IN):: jj
      real*8,intent(IN) :: epsi(NEPS),deps
      real*8,intent(OUT):: FUNC(NEPS)
      integer :: i,j,r,sp,el,dustnr
      real*8 :: vrel,rate,Sr,sum,eq_fak,rat2,stoi,mr
      real*8 :: Sat(NDUST), bmix(NDUST), chi(NDUST)
      real*8 :: Jst(NNUC), Nst(NNUC)
      real*8 :: nsp(NSPECIES),norm(NEPS)
      real*8 :: chinet,dchi,delV,msum,rsum,dum
      real*8,save :: alpha,Afak
      LOGICAL,save :: firstCall=.true.
      COMMON /OUTHELP/Sat,bmix,chi,chinet,Jst

      if (firstCall) then
        alpha = 1.d0
        Afak  = (36.d0*pi)**(1.d0/3.d0)
        elfunk_calls = 0
        firstCall=.false.
      endif
      elfunk_calls = elfunk_calls + 1
 
      do i=1,NELEM
        eps(i) = eps0(i)
      enddo
      do i=1,NEPS
        eps(elnr(i)) = epsi(i)
      enddo
      eps(elnr(jj)) = eps(elnr(jj)) + deps
*
*     ----------------
*     ***  mixing  ***
*     ----------------
      do i=1,NEPS
        el = elnr(i)
        FUNC(i) = nH*(eps0(el)-eps(el))*wmix
        norm(i) = nH*eps0(el)*wmix
c       write(*,*) "Misch",i,elnam(el),FUNC(i)
      enddo

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

*     ------------------------
*     ***  calculate bmix  ***
*     ------------------------
      sum = 0.d0
      do i=1,NDUST
        sum = sum + L4(i)
      enddo
      do i=1,NDUST
        if (sum.gt.0.d0) then
          bmix(i) = L4(i)/sum
        else
          bmix(i) = 0.d0
        endif
      enddo
*
*     -------------------------------
*     ***  consumption by growth  ***
*     -------------------------------
      chinet = 0.d0
      chi = 0.d0
      CALL SUPERSAT(nH,T,Sat)
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
	    rat2 = alpha*nsp(sp)*vrel*stoi          ! surface reac./cm^2/s
	    if (rat2.lt.rate) then
	      rate = rat2                           ! minimum is key 
              mr   = stoi*reac_nu(r,1,2)            ! power for Sr
            endif  
            rsum = rsum + 1.d0/rat2**2
            msum = msum + (stoi*reac_nu(r,1,2))/rat2**2
          endif
        enddo
c       write(*,*) dustnr,msum/rsum,mr
        mr = msum/rsum
        !Sr = Sat(dustnr)    
        Sr = Sat(dustnr)**mr                        ! reaction supersat.ratio
        !eq_fak = 1.d0 - bmix(dustnr)/Sr            ! equilibrium factor 
        !eq_fak = DMAX1(-100.0,eq_fak)              ! de-stiffen evaporation
        eq_fak = 1.d0 - (bmix(dustnr)/Sr)**pred(dustnr) 
        rate = rate * reac_nu(r,1,2)
        do j=1,dust_nel(dustnr)
          el = elcode(dust_el(dustnr,j))
          stoi = dust_nu(dustnr,j)
          FUNC(el) = FUNC(el) - Afak*rhoL2*rate*eq_fak*stoi
c         norm(el) = norm(el) + Afak*rhoL2*rate*stoi
c         write(*,*) "Wachs",r,el,elnam(elnr(el)),
c    &               -Afak*rhoL2*rate*eq_fak*stoi
          delV   = dust_vol(dustnr)                 ! volume increment by reac
          dchi   = Afak * rate * eq_fak * delV
          chinet = chinet + dchi
          chi(dustnr) = chi(dustnr) + dchi
        enddo
      enddo

      

*     -----------------------------------
*     ***  consumption by nucleation  ***
*     -----------------------------------
      CALL NUCLEATION(T,Jst,Nst)
c     if (chinet.lt.0.d0) Jst=0.d0
      do i=1,NNUC
        do j=1,nuc_nel(i)
          el = elcode(nuc_el(i,j))
          FUNC(el) = FUNC(el) - nuc_nu(i,j)*DBLE(nl)*Jst(i)
c         norm(el) = norm(el) + nuc_nu(i,j)*DBLE(nl)*Jst(i)
c         write(*,*) "Nucle",el,elnam(elnr(el)),
c    &               -nuc_nu(i,j)*DBLE(nl)*Jst(i)
        enddo
      enddo

*     ---------------------------------------
*     ***  make deviations dimensionless  ***
*     ---------------------------------------
      do i=1,NEPS
        FUNC(i) = FUNC(i)/norm(i)
      enddo

c     write(*,'(99(1pE11.4))') FUNC
*
      RETURN 
      end
