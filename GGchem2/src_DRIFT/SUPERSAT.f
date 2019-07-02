*********************************************************************
      SUBROUTINE SUPERSAT(nH,T,Sat)
*********************************************************************
      use drift_data,ONLY: NDUST,NMOLE,cmol,nat,nmol,
     &                     dust_nam,dust_nel,dust_el,dust_nu
      implicit none
      real*8,intent(IN) :: nH,T
      real*8,intent(OUT):: Sat(NDUST)
      real*8,parameter :: bk=1.380662d-16,atm=1.013d+6,rgas=8.31434d+0 
      real*8,parameter :: bar=1.d+6,cal=4.184d0 
      real*8,parameter :: mmHg=1.333d+3 ! mmHg --> dyn/cm2
      real*8  :: TT,kT,dG,lbruch,lresult,pst,psat,term,ex
      integer :: i,j,stindex,el,Na
      integer,save :: TiO2,SiO
      logical,save :: firstCall=.true.
*
      data Na/11/ 
*
      if (firstCall) then
        TiO2 = stindex(cmol,NMOLE,'TIO2     ')
        SiO  = stindex(cmol,NMOLE,'SIO      ')
        firstCall=.false.
      endif    
      TT  = DMAX1(T,300.d0)
      kT  = bk*TT
      do i=1,NDUST
        if      (dust_nam(i).eq.'TiO2[s]   ') then
c         ---------------------
c         ***  eigener Fit  ***
c         ---------------------
c         pst = bar
c         dG = +5.04970d+5/TT
c    &         -1.92151d+6 
c    &         +4.64839d+2*TT 
c    &         -9.05840d-3*TT**2
c    &         +9.69662d-7*TT**3
c         dG = dG/(rgas*TT)
c         bruch = 1.d0
c         do j=1,dust_nel(i)
c           el    = dust_el(i,j)
c           bruch = bruch * (nat(el)*kT/pst)**dust_nu(i,j)
c         enddo
c         result = bruch*EXP(-dG)
c         Sat(i) = REAL(result,kind=8)
c         write(*,*) nH,T,Sat(i)
          psat = DEXP(35.8027d0 - 74734.7d0/TT)
          Sat(i) = nmol(TiO2)*kT/psat
c         write(*,*) nH,T,Sat(i)
c         write(*,*) nmol(TiO2)
        else if (dust_nam(i).eq.'Mg2SiO4[s]') then
c         ------------------------------
c         ***  Sharp & Huebner 1990  ***
c         ------------------------------
          pst = atm
          dG  = 7.52334d+4/TT 
     &         -9.38369d+5 
     &         +2.47581d+2*TT
     &         -3.14980d-3*TT**2
          dG = dG/(rgas/cal*TT)
          lbruch = 0.d0
          do j=1,dust_nel(i)
            el     = dust_el(i,j)
            lbruch = lbruch + DBLE(dust_nu(i,j))*DLOG(nat(el)*kT/pst)
          enddo
          Sat(i) = DEXP(lbruch-dG)
c         write(*,*) nH,T,dG,lbruch,Sat(i)
         else if (dust_nam(i).eq.'MgO[s]') then 
c         ------------------------------
c         ***  Sharp & Huebner 1990  ***
c         ------------------------------
          pst = atm
          dG  = 4.38728d+3/TT
     &         -2.38741d+5
     &         +6.86582d+1*TT 
     &         -1.19852d-3*TT**2
     &         +5.72304d-8*TT**3
          dG = dG/(rgas/cal*TT)
          lbruch = 0.d0
          do j=1,dust_nel(i)
            el     = dust_el(i,j)
            lbruch = lbruch + DBLE(dust_nu(i,j))*DLOG(nat(el)*kT/pst)
          enddo
          Sat(i) = DEXP(lbruch-dG)
        else if (dust_nam(i).eq.'MgSiO3[s]') then
c         ------------------------------
c         ***  Sharp & Huebner 1990  ***
c         ------------------------------
          pst = atm
          dG  = 8.74400E+3/TT 
     &         -6.92565E+5
     &         +1.77877E+2*TT
     &         -1.41412E-3*TT**2
          dG = dG/(rgas/cal*TT)
          lbruch = 0.d0
          do j=1,dust_nel(i)
            el    = dust_el(i,j)
            lbruch = lbruch + DBLE(dust_nu(i,j))*DLOG(nat(el)*kT/pst)
          enddo
          Sat(i) = DEXP(lbruch-dG)
        else if (dust_nam(i).eq.'SiO2[s]   ') then
c         ------------------------------
c         ***  Sharp & Huebner 1990  ***
c         ------------------------------
          pst = atm
          dG =  -4.44364d+5 
     &          +1.08531d+2*TT
     &          -6.59213d-4*TT**2
          dG = dG/(rgas/cal*TT)
          lbruch = 0.d0
          do j=1,dust_nel(i)
            el     = dust_el(i,j)
            lbruch = lbruch + DBLE(dust_nu(i,j))*DLOG(nat(el)*kT/pst)
          enddo
          Sat(i) = DEXP(lbruch-dG)
c         write(*,*) nH,T,dG,bruch,result,Sat(i)
        else if (dust_nam(i).eq.'SiO[s]   ') then
c         -------------------
c         ***  Nuth 200x  ***
c         -------------------
          psat = DEXP(17.56d0 - 40760.0d0/TT) * atm
          Sat(i) = nmol(SiO)*kT/psat
c         write(*,*) nH,T,"SiO",psat,nmol(SiO),Sat(i)
        else if (dust_nam(i).eq.'Fe[s]') then 
c         ---------------------
c         ***  eigener Fit  ***
c         ---------------------
          pst = bar
          dG = 7.37828d+5/TT 
     &        -4.22183d+5 
     &        +1.71919d+2*TT
     &        -1.76037d-2*TT**2 
     &        +2.31459d-6*TT**3
          dG = dG/(rgas*TT)
          el = dust_el(i,1)
          lbruch = DLOG(nat(el)*kT/pst)
          Sat(i) = DEXP(lbruch-dG)
c         write(*,*) nat(el),T,dG,Sat(i)
        else if (dust_nam(i).eq.'Al2O3[s]') then 
c         ------------------------------
c         ***  Sharp & Huebner 1990  ***
c         ------------------------------
          pst = atm
          dG  =-7.32976d+5 
     &         +1.84782d+2*TT 
     &         -2.57313d-3*TT**2
          dG = dG/(rgas/cal*TT)
          lbruch = 0.d0
          do j=1,dust_nel(i)
            el     = dust_el(i,j)
            lbruch = lbruch + DBLE(dust_nu(i,j))*DLOG(nat(el)*kT/pst)
          enddo
          Sat(i) = DEXP(lbruch-dG)
        else if (dust_nam(i).eq.'CaTiO3[s]') then 
c         ------------------------------
c         ***  Sharp & Huebner 1990  ***
c         ------------------------------
          pst = atm
          dG  = 1.19107d+4/TT 
     &         -7.30327d+5 
     &         +1.75930d+2*TT
     &         -2.84630d-3*TT**2 
     &         +1.10392d-7*TT**3
          dG = dG/(rgas/cal*TT)
          lbruch = 0.d0
          do j=1,dust_nel(i)
            el     = dust_el(i,j)
            lbruch = lbruch + DBLE(dust_nu(i,j))*DLOG(nat(el)*kT/pst)
          enddo
          Sat(i) = DEXP(lbruch-dG)
c         write(*,*) nH,T,'CaTiO3',dG,lbruch,Sat(i)
        else if (dust_nam(i).eq.'CaSiO3[s]') then 
c         ------------------------------
c         ***  Sharp & Huebner 1990  ***
c         ------------------------------
          pst = atm
          dG  = 6.37937d+4/TT 
     &         -7.21819d+5 
     &         +1.77647d+2*TT
     &         -2.59254d-3*TT**2 
          dG = dG/(rgas/cal*TT)
          lbruch = 0.d0
          do j=1,dust_nel(i)
            el     = dust_el(i,j)
            lbruch = lbruch + DBLE(dust_nu(i,j))*DLOG(nat(el)*kT/pst)
          enddo
          Sat(i) = DEXP(lbruch-dG)
c         write(*,*) nH,T,'CaSiO3',dG,lbruch,Sat(i)
        else if (dust_nam(i).eq.'Fe2SiO4[s]') then 
c         ------------------------------
c         ***  Sharp & Huebner 1990  ***
c         ------------------------------
          pst = atm
          dG  = 6.84477d+4/TT
     &         -9.02146d+5
     &         +2.51045d+2*TT 
     &         -4.44028d-3*TT**2
          dG = dG/(rgas/cal*TT)
          lbruch = 0.d0
          do j=1,dust_nel(i)
            el     = dust_el(i,j)
            lbruch = lbruch + DBLE(dust_nu(i,j))*DLOG(nat(el)*kT/pst)
          enddo
          Sat(i) = DEXP(lbruch-dG)
c          write(*,*) nH,T,'Fe2SiO4',dG,lbruch,Sat(i)
        else if (dust_nam(i).eq.'FeO[s]') then 
c         ------------------------------
c         ***  Sharp & Huebner 1990  ***
c         ------------------------------
          pst = atm
          dG  =-6.05591d+4/TT
     &         -2.22156d+5
     &         +6.78769d+1*TT 
     &         -1.52287d-3*TT**2
          dG = dG/(rgas/cal*TT)
          lbruch = 0.d0
          do j=1,dust_nel(i)
            el     = dust_el(i,j)
            lbruch = lbruch + DBLE(dust_nu(i,j))*DLOG(nat(el)*kT/pst)
          enddo
          Sat(i) = DEXP(lbruch-dG)
c         write(*,*) nH,T,'FeO',dG,lbruch,Sat(i)
        else if (dust_nam(i).eq.'FeS[s]') then 
c         ------------------------------
c         ***  Sharp & Huebner 1990  ***
c         ------------------------------
          pst = atm
          dG  =-6.83787d+4/TT
     &         -1.8862d+5
     &         +6.7059d+1*TT 
     &         -2.318938d-3*TT**2
          dG = dG/(rgas/cal*TT)
          lbruch = 0.d0
          do j=1,dust_nel(i)
            el     = dust_el(i,j)
            lbruch = lbruch + DBLE(dust_nu(i,j))*DLOG(nat(el)*kT/pst)
          enddo
          Sat(i) = DEXP(lbruch-dG)
c         write(*,*) nH,T,'FeS',dG,lbruch,Sat(i)
        else if (dust_nam(i).eq.'Fe2O3[s]') then 
c         ------------------------------
c         ***  Sharp & Huebner 1990  ***
c         ------------------------------
          pst = atm
          dG  =-9.09265d+1/TT
     &         -5.75261d+5
     &         +1.85908d+2*TT 
     &         -7.14322d-3*TT**2
     &         +8.63059d-7*TT**3
          dG = dG/(rgas/cal*TT)
          lbruch = 0.d0
          do j=1,dust_nel(i)
            el     = dust_el(i,j)
            lbruch = lbruch + DBLE(dust_nu(i,j))*DLOG(nat(el)*kT/pst)
          enddo
          Sat(i) = DEXP(lbruch-dG)
c         write(*,*) nH,T,'Fe2O3',dG,lbruch,Sat(i)
        else if (dust_nam(i).eq.'Na[s]') then 
c         ----------------------------------------------------------
c         ***  Na[s] Martinez, Ferguson, Heist, Nuth III (2005)  ***
c              psat in mm Hg!
c         ----------------------------------------------------------          
          psat = mmHg*DEXP(10.86423d0 - 5619.406d0/TT  + 3.45d-6*TT)
          Sat(i) = nat(Na)*kT/psat
        else if (dust_nam(i).eq.'Na2SiO3[s]') then 
c         ------------------------------
c         ***  Na2SiO3[s] Sharp & Huebner 1990  ***
c         ------------------------------
          pst = atm
          dG  = 4.63483d+4/TT
     &         -7.12298d+5
     &         +2.06928d+2*TT 
     &         -4.88925d-3*TT**2
          dG = dG/(rgas/cal*TT)
          lbruch = 0.d0
          do j=1,dust_nel(i)
            el     = dust_el(i,j)
            lbruch = lbruch + DBLE(dust_nu(i,j))*DLOG(nat(el)*kT/pst)
          enddo
          Sat(i) = DEXP(lbruch-dG)
c         write(*,*) nH,T,'Na2SiO3',dG,lbruch,Sat(i)
        else if (dust_nam(i).eq.'SiC[s]') then
          !----------------------------------------
          !***  SiC[s]: eigener Fit nach JANAF  ***
          !----------------------------------------
          pst = bar
          dG  = 6.73337E+5/TT
     &         -1.24381E+6
     &         +3.21779E+2*TT
     &         -4.54405E-3*TT**2
     &         +2.69711E-7*TT**3
          dG = dG/(rgas*TT)
          lbruch = 0.0
          do j=1,dust_nel(i)
            el     = dust_el(i,j)
            term   = nat(el)*kT/pst
!	    print*, "SiC, el, nat", el, nat(el)
            lbruch = lbruch + DLOG(term)*dust_nu(i,j)
          enddo
          Sat(i) = EXP(lbruch-dG)
!	  write(*,*) "Sat(i)", Sat(i)
         else if (dust_nam(i).eq.'TiC[s]') then
c         ----------------------------------------
c         ***  TiC[s]: eigener Fit nach JANAF  ***
c         ----------------------------------------
          pst = bar
          dG =  +1.11878E+5/TT
     &          -1.37612E+6
     &          +3.20666E+2*TT
     &          -4.63379E-3*TT**2
     &          +1.85306E-7*TT**3
          dG = dG/(rgas*TT)
          lbruch = 0.0
          do j=1,dust_nel(i)
            el     = dust_el(i,j)
            term   = nat(el)*kT/pst
!	    print*, "TiC, el, nat(el)", el, nat(el)
            lbruch = lbruch + DLOG(term)*dust_nu(i,j)
          enddo
          Sat(i) = EXP(lbruch-dG)
!	print*, "Sat(i)", Sat(i)
        else if (dust_nam(i).eq.'C[s]') then
!         -------------------------------
!         *** Woitke via JANAF tables ***
!         -------------------------------
          ex   = 1.01428e+06/TT
     &           - 7.23043e+05
     &           + 1.63039e+02*TT
     &           - 1.75890e-03*TT**2
     &           + 9.97416e-08*TT**3
          psat = 1.e+6*EXP(ex/(rgas*TT))
	  Sat(i) = nat(6)*kT/psat	!C is 6
!	  write(*,*) "C, Sat(C)", nat(6), Sat(i)
        else
           stop 'unbekannte Staubsorte in SUPERSAT.'
        endif
       enddo
      RETURN 
      end 
