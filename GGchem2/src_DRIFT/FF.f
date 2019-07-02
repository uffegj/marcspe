**********************************************************************
      SUBROUTINE FF(NN,z,y,f,iflag)
**********************************************************************
*****                                                            *****
*****  Berechnet rechte Seite der gewoehnlichen Dgl: y' = f(y,z) *****
*****                                                            *****
*****  EINGABE:   NN = Dimension von y() und f()                 *****
*****              z = unabhaengige Variable                     *****
*****            y() = Funktionsvektor bei z                     *****
*****          iflag = type of call (1=ordinary, 2=Jacobian,     *****
*****                  -1=from OUTPUT -2=from static_weather)    *****
*****                                                            *****
*****  AUSGABE:  f() = dy()/dz                                   *****
*****                                                            *****
**********************************************************************
      use drift_data,ONLY: NMOM,NDUST,NEPS,NSPECIES,NNUC,
     &                     welche,abschluss,no_converge,
     &                     pi,bk,amu,cref,Lref,Vl,dust_rho,
     &                     elnr,spnr,nat,nmol,nsp,
     &                     T,rho,nH,p,pel,Sat,Nst,Jst,rhod,eps,
     &                     wmix,mu,g,chinet,chi,bmix,LL,L4,rhoL2,
     &                     imp_nuklea,imp_growth,imp_drift,imp_misch,
     &                     nuc_nam
      implicit none
      integer,intent(IN):: NN,iflag
      real*8,intent(IN) :: z,y(NN)
      real*8,intent(OUT):: f(NN)
      real*8 :: epsi(NEPS)
      real*8 :: cc,sum,xi,xi1
      real*8 :: ABSCHLUSS_BED,CHI_NET,Jsum
      real*8,save :: xi_const
      logical:: IS_NAN
      integer:: i,j
      integer,save :: TiO2,Carb
      logical,save :: firstCall=.true.

      if (firstCall) then
        xi_const  = DSQRT(pi)/2.d0*(3.d0/(4.d0*pi))**(1.d0/3.d0)
        firstCall = .false.
        do i=1,NNUC
          if (nuc_nam(i).eq.'TiO2   ') then 
            TiO2=i 
          else if (nuc_nam(i).eq.'Carb   ') then 
            Carb=i 
          else 
            write(*,*) "NO NUKLEATION SPECIES FOUND in FF.f"
          endif
        enddo
      endif   

      if (no_converge) then
        f(:) = 0.d0
        RETURN
      endif

*     --------------------------------------------------------
*     ***  Dichte rho,n<H>,T,mu,wmix als Funktion von z    ***
*     ***  Staubmomente aus Variablen berechnen            ***
*     ***  Staubmaterialdichte rhod berechnen              ***
*     --------------------------------------------------------
c     p = y(NN)*pref
      call THERMO(z,p,T,nH,rho,g,mu,wmix,welche)
      cc = DSQRT(2.d0*bk*T/(mu*amu))
      do i=1,NMOM
        LL(i) = cc/cref * y(i)*Lref(i)
      enddo
      sum  = 0.d0
      rhod = 0.d0
      do i=1,NDUST
        L4(i) = cc/cref * y(NMOM+i)*Lref(4)
        sum = sum + L4(i)
      enddo
      do i=1,NDUST
        if (sum.gt.0.d0) then
          bmix(i) = L4(i)/sum
        else
          bmix(i) = 0.d0
          if (i.eq.1) bmix(i)=1.0     ! => rhod=rhod(TiO2) if no dust
        endif
        rhod = rhod + bmix(i)*dust_rho(i)
      enddo
      
      LL(:) = LL(:)/rhod
      L4(:) = L4(:)/rhod

      LL(0) = ABSCHLUSS_BED(abschluss,Vl,NMOM,LL)

      rhoL2 = rho*LL(2)

 
*     ----------------------------------------
*     ***  statische Elementhaeufigkeiten  ***
*     ----------------------------------------
      call STATIC_EPS(epsi,iflag)
      if (no_converge) then
        f(:) = 0.d0
        RETURN
      endif    
      do i=1,NEPS
        eps(elnr(i)) = epsi(i)
      enddo

*     -----------------------------------
*     ***  GG-Chemie und Keimbildung  ***
*     -----------------------------------
      call GGCHEM(nH,T,pel)
      do i=1,NSPECIES
        j = spnr(i)
        if (j.gt.1000) then
          nsp(i) = nat(j-1000)
        else
          nsp(i) = nmol(j)
        endif
      enddo 
      CALL NUCLEATION(T,Jst,Nst)
      Jsum= 0.0
      do i=1, NNUC
        Jsum = Jsum + Jst(i)
      enddo

*     ----------------------------
*     ***  Wachstum und Drift  ***
*     ----------------------------
      CALL SUPERSAT(nH,T,Sat)
      chinet = CHI_NET(T,Sat,L4,chi,bmix)
c     if (chinet.lt.0.d0) Jst=0.d0
      xi  = xi_const*g
      xi1 = 1.d0/xi
      f(1) = Jsum * xi1
      do i=1,NMOM-1
        f(i+1) = ( Vl**(DBLE(i)/3.d0)*Jsum
     &           + DBLE(i)/3.d0*chinet*rho*LL(i-1) ) * xi1
      enddo
      do i=1,NDUST
        f(NMOM+i) = chi(i)*rho*LL(2) * xi1
      enddo

!   Achtung: hier muss die Nucleationsspecies als          !
!   erstes in der Liste der Solids in  DustChem.dat stehen ! 
      f(NMOM+TiO2) = f(NMOM+TiO2) + Jst(TiO2)*Vl * xi1  ! TiO2
c     f(NMOM+Carb) = f(NMOM+Carb) + Jst(Carb)*Vl * xi1  ! Carbon

*     ------------------
*     ***  Mischung  ***
*     ------------------
      do i=1,NMOM
        f(i) = f(i) - rho*LL(i-1)*wmix * xi1
      enddo

      do i=1,NDUST
        f(NMOM+i) = f(NMOM+i) - rho*LL(3)*wmix*bmix(i) * xi1
      enddo


c      imp_nuklea = Vl*Jst * xi1
      imp_nuklea = Vl*Jsum * xi1
      imp_growth = chinet*rho*LL(2) * xi1
      imp_drift  = -f(4)
      imp_misch  = -rho*LL(3)*wmix * xi1
 
*     -----------------------------------------
*     ***  optische Tiefe des TiO2-Staubes  ***
*     -----------------------------------------
c     f(NN-1) = 3.d0/4.d0 * QPlanck_TiO2(T) * rho*LL(3)  ! Planck-mean
c     f(NN-1) = 3.d0/4.d0 * QRoss_TiO2(T) * rho*LL(3)    ! Rosseland-mean
c     f(NN-1) = 3.d0/4.d0 * 5.6d0 * rho*LL(3)            ! typisch 
 
*     ----------------------------
*     ***  hydrostatisches GG  ***
*     ----------------------------
c     f(NN) = g*rho
 
*     ------------------------------
*     ***  dimensionslos machen  ***
*     ------------------------------
      do i=1,NMOM
        f(i) = f(i)*cref/Lref(i)
      enddo
      do i=1,NDUST
        f(NMOM+i) = f(NMOM+i)*cref/Lref(4)
      enddo
c     f(NN) = f(NN)/pref
*
      do i=1,NN
        if (IS_NAN(f(i))) then
          write(*,*) '*** NaN in SUBROUTINE FF'
          write(*,*) 'rho and Temp:'
          write(*,*) rho,T
          write(*,*) 'gas eps:'
          write(*,*) epsi
          write(*,*) 'Teilchendichten: '
          write(*,*) (nsp(j),j=1,NSPECIES)  
          write(*,*) 'Staubmomente: '
          write(*,*) (LL(j),j=0,4)
          write(*,*) 'Jst(TiO2), Jst(Carb) and chinet: '
          write(*,*) Jst(TiO2),Jst(Carb),chinet
          write(*,*) 'yy = '
          write(*,*) (y(j),j=1,NN)  
          write(*,*) 'FF = '
          write(*,*) (f(j),j=1,NN)  
          stop
        endif
      enddo
    
      RETURN 
      END 
