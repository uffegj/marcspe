***********************************************************************
      SUBROUTINE GGCHEM (nHges,Tg, pel)
***********************************************************************
*****                                                             *****
*****  Ruft lediglich SMCHEM auf (mit kompatiber Datenstruktur)   *****
*****                                                             *****
***********************************************************************
      use drift_data,ONLY: NELEM,NMOLE,nel,nat,nion,nmol,eps
      implicit  none
      real*8,intent(IN) :: nHges,Tg
      real*8,intent(out):: pel
      integer   nml,nelm
      parameter (nelm=17, nml=167)
      integer   H,He,C,N,O,Ne,Na,Mg,Al,Si,S,K,Ca,Cr,Mn,Fe,Ni,Ti,Li,Fl,Cl
      real*8    epsi(nelm),anmono(nelm)
      integer   HII,CII,NII,OII,NaII,MgII,AlII,KII,TiII,SII,SiII,FeII
      integer   CaII
      common/ionnumm/HII,CII,NII,OII,NaII,MgII,AlII,KII,TiII,SII,SiII,
     &               FeII,CaII
*     -----------------------------------------------------------------
      data H/1/, He/2/, Li/3/, C/6/, N/7/, O/8/, Fl/9/, Ne/10/, Na/11/ 
      data Mg/12/, Al/13/, Si/14/, S/16/, Cl/17/, K/19/, Ca/20/
      data Cr/24/, Mn/25/, Fe/26/, Ni/28/, Ti/22/
*     -----------------------------------------------------------------
*
      epsi( 1) = eps(He)
      epsi( 2) = 0.d0
      epsi( 3) = 1.d0
      epsi( 4) = eps( C)
      epsi( 5) = eps( N)
      epsi( 6) = eps( O)
      epsi( 7) = eps(Si)
      epsi( 8) = eps(Mg)
      epsi( 9) = eps(Al)
      epsi(10) = eps(Fe)
      epsi(11) = eps( S)
      epsi(12) = eps(Na)
      epsi(13) = eps( K)
      epsi(14) = eps(Ti)
      epsi(15) = eps(Ca)
      epsi(16) = eps(Li)
      epsi(17) = eps(Cl)
*      
      call SMCHEM ( nHges, Tg, epsi, anmono, nmol, pel )
*
      nat(He)  = anmono( 1)
      nel      = anmono( 2)
      nat( H)  = anmono( 3)
      nat( C)  = anmono( 4)
      nat( N)  = anmono( 5)
      nat( O)  = anmono( 6)
      nat(Si)  = anmono( 7)
      nat(Mg)  = anmono( 8)
      nat(Al)  = anmono( 9)
      nat(Fe)  = anmono(10)
      nat( S)  = anmono(11)
      nat(Na)  = anmono(12)
      nat( K)  = anmono(13)
      nat(Ti)  = anmono(14)
      nat(Ca)  = anmono(15)
      nat(Li)  = anmono(16)
      nat(Cl)  = anmono(17)
      nion( H) = nmol(HII)
      nion( C) = nmol(CII)
      nion( N) = nmol(NII)
      nion( O) = nmol(OII)
      nion(Na) = nmol(NaII)
      nion(Mg) = nmol(MgII)
      nion(Al) = nmol(AlII)
      nion( K) = nmol(KII)
      nion(Ti) = nmol(TiII)
      nion(S ) = nmol(SII)
      nion(Si) = nmol(SiII)
      nion(Fe) = nmol(FeII)
      nion(Ca) = nmol(CaII)

c      print*, 'in GGchem: nat(6)', nat(6)
c      write(*,*) nHges,Tg,eps(C)/eps(O)
*             
      RETURN
      end 
