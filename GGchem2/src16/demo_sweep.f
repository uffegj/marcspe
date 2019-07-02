***********************************************************************
      SUBROUTINE DEMO_SWEEP
***********************************************************************
      use PARAMETERS,ONLY: Tmin,Tmax,pmin,pmax,nHmin,nHmax,useDatabase,
     >                     model_eqcond,model_pconst,Npoints,
     >                     remove_condensates
      use CHEMISTRY,ONLY: NELM,NMOLE,elnum,cmol,catm,el,charge
      use DUST_DATA,ONLY: NELEM,NDUST,elnam,eps0,bk,bar,muH,amu,
     >                    dust_nel,dust_el,dust_nu,dust_nam,dust_mass,
     >                    dust_Vol,mass,mel
      use EXCHANGE,ONLY: nel,nat,nion,nmol,mmol,H,C,N,O,W
      implicit none
      integer,parameter :: qp = selected_real_kind ( 33, 4931 )
      real :: p,Tg,rhog,rhod,dustV,nHges,nges,mges,kT,pges
      real :: nTEA,pTEA,mu,muold,Jstar,Nstar,
     >        ntotat1,ntotat2,ntotmol,ntot,ppel
      real(kind=qp) :: eps(NELEM),eps00(NELEM),ppat(NELM),ppmol(NMOLE)
      real(kind=qp) :: Sat(NDUST),eldust(NDUST),out(NDUST)
      real(kind=qp) :: fac,e_reservoir(NELEM),d_reservoir(NDUST)
      integer :: i,j,jj,k,l,NOUT,iW,stindex
      character(len=5000) :: species,NISTspecies,elnames
      character(len=20) :: frmt,name,short_name(NDUST),test1,test2
      character(len=4) :: sup
      character(len=2) :: test3
      character(len=1) :: char
      integer :: verbose=0
      logical :: isOK,hasW,same,TEAinterface=.false.

      real ::
     > ppN2,ppCH,ppCO,ppCN,ppC2,ppNO,ppC2H2,ppHCN,ppC2H,ppC3,ppCS,ppH2O,
     > ppOH,ppTiO,ppSiO,ppCH4,ppNH,ppSiH,ppFeH,ppVO,ppZrO,ppMgH,ppNH3,
     > ppCO2,ppTiH,ppCaH,ppCrH,ppLiH,ppH,ppO2,ppHm,ppH2,ppH2p,ppHS,
     > ppC3H,ppSiC,ppSiC2,ppNS,ppSiN,ppSO,ppS2,ppSiS,ppLaO,ppCH2,
     > ppCH3,ppSi2C,ppSiO2,ppH2S,ppCaOH,ppCHNO,ppSiF2


      integer ::
     > iN2,iCH,iCO,iCN,iC2,iNO,iC2H2,iHCN,iC2H,iC3,iCS,iH2O,
     > iOH,iTiO,iSiO,iCH4,iNH,iSiH,iFeH,iVO,iZrO,iMgH,iNH3,iCO2,
     > iTiH,iCaH,iCrH,iLiH,iH,iO2,iHm,iH2,iH2p,iHS,iC3H,iSiC,iSiC2,iNS,
     > iSiN,iSO,iS2,iSiS,iLaO,iCH2,iCH3,iSi2C,iSiO2,iH2S,iCaOH,
     > iCHNO,iSiF2


      common /partialpressure/
     > ppN2,ppCH,ppCO,ppCN,ppC2,ppNO,ppC2H2,ppHCN,ppC2H,ppC3,ppCS,ppH2O,
     > ppOH,ppTiO,ppSiO,ppCH4,ppNH,ppSiH,ppFeH,ppVO,ppZrO,ppMgH,ppNH3,
     > ppCO2,ppTiH,ppCaH,ppCrH,ppLiH,ppH,ppO2,ppHm,ppH2,ppH2p,ppHS,
     > ppC3H,ppSiC,ppSiC2,ppNS,ppSiN,ppSO,ppS2,ppSiS,ppLaO,ppCH2,
     > ppCH3,ppSi2C,ppSiO2,ppH2S,ppCaOH,ppCHNO,ppSiF2


      !----------------------------
      ! ***  open output files  ***
      !----------------------------
      do i=1,NDUST
        name = dust_nam(i) 
        j=index(name,"[s]")
        short_name(i) = name
        if (j>0) short_name(i)=name(1:j-1)
      enddo
      eps  = eps0
      NOUT = NELM
      if (charge) NOUT=NOUT-1
      open(unit=70,file='Static_Conc.dat',status='replace')
      write(70,1000) 'H',eps( H), 'C',eps( C),
     &               'N',eps( N), 'O',eps( O)
      write(70,*) NOUT,NMOLE,NDUST,Npoints
      write(70,2000) 'Tg','nHges','pges','el',
     &               (trim(elnam(elnum(j))),j=1,el-1),
     &               (trim(elnam(elnum(j))),j=el+1,NELM),
     &               (trim(cmol(i)),i=1,NMOLE),
     &               ('S'//trim(short_name(i)),i=1,NDUST),
     &               ('n'//trim(short_name(i)),i=1,NDUST),
     &               ('eps'//trim(elnam(elnum(j))),j=1,el-1),
     &               ('eps'//trim(elnam(elnum(j))),j=el+1,NELM),
     &               'dust/gas','dustVol/H','Jstar(W)','Nstar(W)'

      if (TEAinterface) then
      !--- TEA automated choice from dispol_large.dat ---
      species = "H_g He_ref C_g N_g O_g Si_g S_g Na_g "
     &        //"Ca_g Cl_g Ti_g K_g Al_g Mg_g Fe_g Li_g "
     &        //"F_g P_g Cr_g Mn_g Ni_g W_g"
      elnames = " H He C N O Na Mg Si Fe Al Ca Ti S Cl K Li "
     &        //" F P V Cr Mn Ni Zr W " 
      do i=1,NMOLE
        name = trim(cmol(i))
        isOK = .true.
        if (index(name,"+")>0) isOK=.false.
        if (index(name,"-")>0) isOK=.false.
        !print*,trim(name),isOK
        jj   = len(trim(name))
        j    = 1
        do while(j<=jj) 
          test1 = name(j:j)
          test2 = name(j:j+1)
          l = iachar(name(j:j))
          if (l<iachar("0").or.l>iachar("9")) then
            test3 = ''
            if (trim(test2)=='LI') test3='Li'
            if (trim(test2)=='MG') test3='Mg'
            if (trim(test2)=='AL') test3='Al'
            if (trim(test2)=='SI') test3='Si'
            if (trim(test2)=='CL') test3='Cl'
            if (trim(test2)=='CA') test3='Ca'
            if (trim(test2)=='TI') test3='Ti'
            if (trim(test2)=='FE') test3='Fe'
            if (trim(test2)=='HE') test3='He'
            if (trim(test2)=='CR') test3='Cr'
            if (trim(test2)=='MN') test3='Mn'
            if (trim(test2)=='NI') test3='Ni'
            if (trim(test2)=='ZR') test3='Zr'
            if (test3.ne.'') then
              name(j:j+1) = test3(1:2) 
              test1 = test3
              j = j+1
            endif
            if (index(elnames," "//trim(test1)//" ")<=0) isOK=.false.
            !print*,trim(name),j,trim(test1),isOK
          endif  
          j = j+1
        enddo  
        if (trim(name)=="C3H") isOK=.false.  
        if (isOK) then
          sup = "_g"
          if (trim(name)=="HO2")   name="HOO" 
          if (trim(name)=="C2N2")  name="(CN)2" 
          if (trim(name)=="CHNO")  name="HNCO" 
          if (trim(name)=="HNO3")  name="HONO2"
          if (trim(name)=="N2H2")  name="HNNH"
          if (trim(name)=="CNO")   name="NC101"
          if (trim(name)=="H2SO4") name="O2S(OH)2"
          if (trim(name)=="SN")    name="NS"
          if (trim(name)=="S2O")   name="SSO"
          if (trim(name)=="H2") sup="_ref"
          if (trim(name)=="O2") sup="_ref"
          if (trim(name)=="N2") sup="_ref"
          species=trim(species)//" "//trim(name)//trim(sup)
        endif  
      enddo

      !--- TEA explicit choice ---  
      species = "H_g He_ref H2_ref "                                       ! H and He
     &  //"C_g N_g O_g Na_g Mg_g Si_g Al_g Ca_g Ti_g S_g Fe_g Cl_g "       ! atoms
     &  //"K_g Li_g F_g P_g V_g Cr_g Mn_g Ni_g Zr_g W_g "
     &  //"C2_g C2H2_g C2H4_g C2H_g C3_g C4_g C5_g CH2_g CH3_g "           ! C
     &  //"CH4_g CH_g "
     &  //"N2_ref C2N_g C4N2_g (CN)2_g CN_g CNN_g HCN_g HNNH_g "           ! N
     &  //"N2H4_g N3_g NC101_g NCN_g NH2_g NH3_g NH_g "
     &  //"O2_ref C2H4O_g C2O_g C3O2_g CO2_g CO_g H2CO_g H2O_g HCO_g "     ! O
     &  //"HNCO_g HNO_g HONO2_g HOO_g HOOH_g N2O3_g N2O4_g N2O5_g "
     &  //"N2O_g NO2_g NO3_g NO_g O3_g OH_g "
     &  //"Na2_g (NaCN)2_g NaCN_g NaH_g NaO_g (NaOH)2_g NaOH_g "           ! Na
     &  //"Mg2_g MgH_g MgN_g MgO_g Mg(OH)2_g MgOH_g "                      ! Mg
     &  //"Si2C_g Si2_g Si2N_g Si3_g SiC2_g SiC_g "                        ! Si (SiH2,SiH3 missing)
     &  //"SiH4_g SiH_g SiN_g SiO2_g SiO_g "                               ! (Si(CH3)4_g not working)
     &  //"Fe(CO)5_g FeO_g Fe(OH)2_g "                                     ! Fe (FeH missing)
     &  //"Al2_g Al2O_g AlC_g AlH_g AlN_g (AlO)2_g AlO2_g AlO_g "          ! Al
     &  //"AlOH_g OAlH_g OAlOH_g "
     &  //"Ca2_g CaO_g Ca(OH)2_g CaOH_g "                                  ! Ca (CaH missing)
     &  //"TiO2_g TiO_g "                                                  ! Ti (TiS,TiH,TiC,TiN,TiC2 missing)
     &  //"AlS_g CaS_g COS_g CS2_g CS_g FeS_g H2S_g HS_g MgS_g "           ! S
     &  //"Na2SO4_g NS_g O2S(OH)2_g "
     &  //"S2_g S3_g S4_g S7_g S8_g "                                      ! (S5,S6 not working)
     &  //"SiS_g SO2_g SO3_g SO_g SSO_g "
     &  //"AlCl2_g (AlCl3)2_g AlCl3_g AlCl_g C2Cl2_g C2Cl4_g C2Cl6_g "     ! Cl
     &  //"C2HCl_g CaCl2_g CaCl_g CCl2_g CCl3_g CCl4_g CCl_g CH2Cl2_g "
     &  //"CH3Cl_g CHCl3_g CHCl_g Cl2O_g ClCN_g ClO2_g ClO_g "
     &  //"ClSSCl_g COCl2_g COCl_g (FeCl2)2_g FeCl2_g (FeCl3)2_g "
     &  //"FeCl3_g FeCl_g HCl_g HOCl_g (MgCl2)2_g MgCl2_g MgCl_g "
     &  //"(NaCl)2_g NaCl_g NO2Cl_g OAlCl_g ONCl_g OTiCl_g S2Cl_g "
     &  //"SCl2_g SCl_g SiCl2_g SiCl3_g SiCl4_g SiCl_g "                   ! (SiCH3Cl3_g not working)
     &  //"SiH2Cl2_g SiH3Cl_g SiHCl3_g SO2Cl2_g TiCl2_g TiCl3_g "
     &  //"TiCl4_g TiCl_g TiOCl2_g "
     &  //"K2_g K2SO4_g (KCl)2_g KCl_g (KCN)2_g KCN_g KH_g KO_g "          ! K
     &  //"(KOH)2_g KOH_g "
     &  //"Li2_g Li2O_g Li2SO4_g (LiCl)2_g (LiCl)3_g LiCl_g LiH_g "        ! Li
     &  //"LiN_g (LiO)2_g LiOCl_g LiO_g (LiOH)2_g LiOH_g LiONa_g "
     &  //"LiON_g "
     &  //"AlCl2F_g AlClF2_g AlClF_g AlF2_g (AlF3)2_g AlF3_g AlF_g "       ! F
     &  //"C2F2_g C2F4_g C2F6_g C2HF_g CaF2_g CaF_g CCl2F2_g CCl3F_g "
     &  //"CClF3_g CF2_g CF3CN_g CF3_g CF3OF_g CF3SF5_g CF4_g CF_g "
     &  //"CH2ClF_g CH2F2_g CH3F_g CHCl2F_g CHClF2_g CHF3_g CHF_g "
     &  //"ClF3_g ClF5_g ClF_g ClO3F_g COClF_g COF2_g COF_g FCN_g "
     &  //"FeF2_g FeF3_g FeF_g FONO2_g FSSF_g H3F3_g H4F4_g "
     &  //"H5F5_g H6F6_g H7F7_g HCOF_g (HF)2_g HF_g HOF_g HSO3F_g "
     &  //"(KF)2_g KF_g Li2ClF_g LiAlF4_g (LiF)2_g (LiF)3_g LiF_g "
     &  //"LiOF_g MgClF_g (MgF2)2_g MgF2_g MgF_g N2F4_g NaAlF4_g "
     &  //"(NaF)2_g NaF_g NF2_g NF3_g NF3O_g NF_g NO2F_g O2F_g "
     &  //"OAlF2_g OAlF_g OF2_g OF_g ONF_g OSF2_g OSiF2_g OTiF2_g "
     &  //"OTiF_g S2F10_g SClF5_g SF2_g SF3_g SF4_g SF5_g SF6_g SF_g "
     &  //"SiCH3F3_g SiCl3F_g SiClF3_g SiF2_g SiF3_g SiF4_g SiF_g "
     &  //"SiH2F2_g SiH3F_g SiHF3_g SO2ClF_g SO2F2_g SSF2_g TiF2_g "
     &  //"TiF3_g TiF4_g TiF_g "
     &  //"CHP_g CP_g OPCl3_g P2_g (P2O3)2_g (P2O5)2_g P4_g P4S3_g "       ! P
     &  //"PCl3_g PCl5_g PCl_g PF2_g PF3_g PF5_g PF_g PH2_g PH3_g "
     &  //"PH_g PN_g PO2_g POCl2F_g POClF2_g POF3_g PO_g PSF3_g PSF_g "
     &  //"PS_g SPCl3_g "
     &  //"VN_g VO2_g VO_g "                                               ! V
     &  //"CrN_g CrO2_g CrO3_g CrO_g "                                     ! Cr (CrH missing)
     &  //""                                                               ! Mn (all missing)
     &  //"NiCl2_g NiCl_g NiS_g "                                          ! Ni (NiH,NiO,NiF missing)
     &  //""                                                               ! (Ni(CO)4_g not working)
     &  //"ZrCl2_g ZrCl3_g ZrCl4_g ZrCl_g ZrF2_g ZrF3_g ZrF4_g ZrF_g "     ! Zr
     &  //"ZrH_g ZrN_g ZrO2_g ZrO_g "
     &  //"O2W(OH)2_g OWCl4_g W3O8_g WCl2_g WCl4_g (WCl5)2_g WCl5_g "      ! W
     &  //"WCl6_g WCl_g WF4O_g WF6_g WF_g WO2Cl2_g WO2_g (WO3)2_g "
     &  //"(WO3)3_g (WO3)4_g WO3_g WO_g "

      !--- all NIST species ---
      NISTspecies = "H2_ref H_g He_ref "                                   ! H and He
     &  //"C2_g C2H2_g C2H4_g C2H_g C3_g C4_g C5_g C_g CH2_g CH3_g "       ! C
     &  //"CH4_g CH_g "
     &  //"N2_ref C2N_g C4N2_g (CN)2_g CN_g CNN_g HCN_g HNNH_g N2H4_g "    ! N
     &  //"N3_g NC101_g NCN_g N_g NH2_g NH3_g NH_g "
     &  //"O2_ref C2H4O_g C2O_g C3O2_g CO2_g CO_g H2CO_g H2O_g HCO_g "     ! O
     &  //"HNCO_g HNO_g HONO2_g HOO_g HOOH_g N2O3_g N2O4_g N2O5_g "
     &  //"N2O_g NO2_g NO3_g NO_g O3_g O_g OH_g "
     &  //"Na2_g (NaCN)2_g NaCN_g Na_g NaH_g NaO_g (NaOH)2_g NaOH_g "      ! Na
     &  //"Mg2_g Mg_g MgH_g MgN_g MgO_g Mg(OH)2_g MgOH_g "                 ! Mg
     &  //"Si2C_g Si2_g Si2N_g Si3_g SiC2_g SiC_g Si(CH3)4_g Si_g "        ! Si
     &  //"SiH4_g SiH_g SiN_g SiO2_g SiO_g "    
     &  //"Fe(CO)5_g Fe_g FeO_g Fe(OH)2_g "                                ! Fe
     &  //"Al2_g Al2O_g AlC_g Al_g AlH_g AlN_g (AlO)2_g AlO2_g AlO_g "     ! Al
     &  //"AlOH_g OAlH_g OAlOH_g "
     &  //"Ca2_g Ca_g CaO_g Ca(OH)2_g CaOH_g "                             ! Ca
     &  //"Ti_g TiO2_g TiO_g "                                             ! Ti
     &  //"AlS_g CaS_g COS_g CS2_g CS_g FeS_g H2S_g HS_g MgS_g "           ! S
     &  //"Na2SO4_g NS_g O2S(OH)2_g S2_g S3_g S4_g S5_g S6_g S7_g S8_g "
     &  //"S_g SiS_g SO2_g SO3_g SO_g SSO_g "
     &  //"AlCl2_g (AlCl3)2_g AlCl3_g AlCl_g C2Cl2_g C2Cl4_g C2Cl6_g "     ! Cl
     &  //"C2HCl_g CaCl2_g CaCl_g CCl2_g CCl3_g CCl4_g CCl_g CH2Cl2_g "
     &  //"CH3Cl_g CHCl3_g CHCl_g Cl2O_g ClCN_g Cl_g ClO2_g ClO_g "
     &  //"ClSSCl_g COCl2_g COCl_g (FeCl2)2_g FeCl2_g (FeCl3)2_g "
     &  //"FeCl3_g FeCl_g HCl_g HOCl_g (MgCl2)2_g MgCl2_g MgCl_g "
     &  //"(NaCl)2_g NaCl_g NO2Cl_g OAlCl_g ONCl_g OTiCl_g S2Cl_g "
     &  //"SCl2_g SCl_g SiCH3Cl3_g SiCl2_g SiCl3_g SiCl4_g SiCl_g "
     &  //"SiH2Cl2_g SiH3Cl_g SiHCl3_g SO2Cl2_g TiCl2_g TiCl3_g "
     &  //"TiCl4_g TiCl_g TiOCl2_g "
     &  //"K2_g K2SO4_g (KCl)2_g KCl_g (KCN)2_g KCN_g K_g KH_g KO_g "      ! K
     &  //"(KOH)2_g KOH_g "
     &  //"Li2_g Li2O_g Li2SO4_g (LiCl)2_g (LiCl)3_g LiCl_g Li_g LiH_g "   ! Li
     &  //"LiN_g (LiO)2_g LiOCl_g LiO_g (LiOH)2_g LiOH_g LiONa_g "
     &  //"LiON_g "
     &  //"AlCl2F_g AlClF2_g AlClF_g AlF2_g (AlF3)2_g AlF3_g AlF_g "       ! F
     &  //"C2F2_g C2F4_g C2F6_g C2HF_g CaF2_g CaF_g CCl2F2_g CCl3F_g "
     &  //"CClF3_g CF2_g CF3CN_g CF3_g CF3OF_g CF3SF5_g CF4_g CF_g "
     &  //"CH2ClF_g CH2F2_g CH3F_g CHCl2F_g CHClF2_g CHF3_g CHF_g "
     &  //"ClF3_g ClF5_g ClF_g ClO3F_g COClF_g COF2_g COF_g FCN_g "
     &  //"FeF2_g FeF3_g FeF_g F_g FONO2_g FSSF_g H3F3_g H4F4_g "
     &  //"H5F5_g H6F6_g H7F7_g HCOF_g (HF)2_g HF_g HOF_g HSO3F_g "
     &  //"(KF)2_g KF_g Li2ClF_g LiAlF4_g (LiF)2_g (LiF)3_g LiF_g "
     &  //"LiOF_g MgClF_g (MgF2)2_g MgF2_g MgF_g N2F4_g NaAlF4_g "
     &  //"(NaF)2_g NaF_g NF2_g NF3_g NF3O_g NF_g NO2F_g O2F_g "
     &  //"OAlF2_g OAlF_g OF2_g OF_g ONF_g OSF2_g OSiF2_g OTiF2_g "
     &  //"OTiF_g S2F10_g SClF5_g SF2_g SF3_g SF4_g SF5_g SF6_g SF_g "
     &  //"SiCH3F3_g SiCl3F_g SiClF3_g SiF2_g SiF3_g SiF4_g SiF_g "
     &  //"SiH2F2_g SiH3F_g SiHF3_g SO2ClF_g SO2F2_g SSF2_g TiF2_g "
     &  //"TiF3_g TiF4_g TiF_g "
     &  //"CHP_g CP_g OPCl3_g P2_g (P2O3)2_g (P2O5)2_g P4_g P4S3_g "       ! P
     &  //"PCl3_g PCl5_g PCl_g PF2_g PF3_g PF5_g PF_g P_g PH2_g PH3_g "
     &  //"PH_g PN_g PO2_g POCl2F_g POClF2_g POF3_g PO_g PSF3_g PSF_g "
     &  //"PS_g SPCl3_g "
     &  //"V_g VN_g VO2_g VO_g "                                           ! V
     &  //"Cr_g CrN_g CrO2_g CrO3_g CrO_g "                                ! Cr
     &  //"Mn_g "                                                          ! Mn
     &  //"NiCl2_g NiCl_g Ni(CO)4_g Ni_g NiS_g "                           ! Ni 
     &  //"ZrCl2_g ZrCl3_g ZrCl4_g ZrCl_g ZrF2_g ZrF3_g ZrF4_g ZrF_g "     ! Zr
     &  //"Zr_g ZrH_g ZrN_g ZrO2_g ZrO_g "
     &  //"O2W(OH)2_g OWCl4_g W3O8_g WCl2_g WCl4_g (WCl5)2_g WCl5_g "      ! W
     &  //"WCl6_g WCl_g WF4O_g WF6_g WF_g W_g WO2Cl2_g WO2_g (WO3)2_g "
     &  //"(WO3)3_g (WO3)4_g WO3_g WO_g "

      !species = NISTspecies    ! want all NIST species?

      write(frmt,'("(A",I4.4,")")') len(trim(species))
      open(unit=71,file='ggchem.atm',status='replace')
      write(71,'(A8)') "#SPECIES"
      write(71,frmt) trim(species) 
      write(71,*)
      write(71,'(A8)') "#TEADATA"
      write(71,'(A10,A8,99(A18))') "#Pressure ","Temp",
     &         (trim(elnam(elnum(j))),j=1,el-1),
     &         (trim(elnam(elnum(j))),j=el+1,NELM)
      endif

      iW = stindex(dust_nam,NDUST,'W[s]')
      hasW = (iW>0)
      e_reservoir = 0.Q0
      d_reservoir = 0.Q0
      eps00 = eps0

      !----------------------------------------
      ! ***  run chemistry on linear track  ***
      !----------------------------------------
      mu = muH
      do i=1,1!Npoints
        fac = REAL(i-1)/REAL(Npoints-1) 
        Tg  = EXP(LOG(Tmax)+fac*LOG(Tmin/Tmax))
        same = (Tmin==Tmax)
        if (model_pconst) then
          p = EXP(LOG(pmax)+fac*LOG(pmin/pmax))
          same = same.and.(pmin==pmax)
        else  
          nHges = EXP(LOG(nHmax)+fac*LOG(nHmin/nHmax))
          same = same.and.(nHmin==nHmax)
        endif  
        if (same) then
          eps0(C) = eps0(O) * (0.3+1.1*fac) 
          eps(C)  = eps0(C)
          eps(O)  = eps0(O)
          print*,"C/O=",eps(C)/eps(O)
        endif   
        eldust = 0.Q0

        !--- run chemistry (+phase equilibrium)    ---
        !--- iterate to achieve requested pressure ---
        do 
          if (model_pconst) nHges = p*mu/(bk*Tg)/muH
          if (model_eqcond) then
            call EQUIL_COND(nHges,Tg,eps,Sat,eldust,verbose)
          endif  
          call GGCHEM(nHges,Tg,eps,.false.,0)
          kT = bk*Tg
          nges = nel
          mges = nel*mel
          do j=1,NELEM
            nges = nges + nat(j)
            mges = mges + nat(j)*mass(j)
          enddo
          do j=1,NMOLE
            nges = nges + nmol(j)
            mges = mges + nmol(j)*mmol(j)
          enddo
          pges = nges*kT
          muold = mu
          mu = nHges/pges*(bk*Tg)*muH
          if (.not.model_pconst) exit
!          print '("mu=",2(1pE12.5))',muold/amu,mu/amu
          if (ABS(mu/muold-1.0)<1.E-5) exit
        enddo  

        !--- remove all condensates and put them in reservoir? ---
        if (remove_condensates) then
          fac = 1.Q+0
          do j=1,NDUST
            d_reservoir(j) = d_reservoir(j) + fac*eldust(j)
            do jj=1,dust_nel(j)
              k = dust_el(j,jj)
              e_reservoir(k) = e_reservoir(k) 
     &                       + fac*dust_nu(j,jj)*eldust(j)
            enddo
          enddo  
          !do j=1,NELM
          !  if (j==el) cycle 
          !  k = elnum(j)
          !  print'(A3,2(1pE18.10))',elnam(k),eps(k)/eps00(k),
     &    !                  (eps(k)+e_reservoir(k))/eps00(k)
          !enddo
          eps0(:) = eps(:) + (1.Q0-fac)*e_reservoir(:)
          eldust(:) = d_reservoir(:)
        endif  

        !--- compute supersat ratios and nucleation rates ---
        call SUPERSAT(Tg,nat,nmol,Sat)
        if (hasW) then
          call NUCLEATION('W',Tg,dust_vol(iW),nat(W),
     &                    Sat(iW),Jstar,Nstar)
        else
          Jstar = 0
          Nstar = 9.e+99
        endif  

        !--- compute dust/gas density ratio ---
        rhog  = mges
        rhod  = 0.0
        dustV = 0.0
        do jj=1,NDUST
          rhod  = rhod  + nHges*eldust(jj)*dust_mass(jj)
          dustV = dustV + eldust(jj)*dust_Vol(jj)
          out(jj) = LOG10(MIN(1.Q+300,MAX(1.Q-300,Sat(jj))))
          if (ABS(Sat(jj)-1.Q0)<1.E-10) out(jj)=0.Q0
        enddo  

!        print'(i4," Tg[K] =",0pF8.2,"  n<H>[cm-3] =",1pE10.3)',
!     >        i,Tg,nHges

!        write(*,1010) ' Tg=',Tg,' n<H>=',nHges,
!     &                ' p=',pges/bar,' mu=',mu/amu,
!     &                ' dust/gas=',rhod/rhog
!        print*
        write(70,2010) Tg,nHges,pges,
     &       LOG10(MAX(1.Q-300, nel)),
     &      (LOG10(MAX(1.Q-300, nat(elnum(jj)))),jj=1,el-1),
     &      (LOG10(MAX(1.Q-300, nat(elnum(jj)))),jj=el+1,NELM),
     &      (LOG10(MAX(1.Q-300, nmol(jj))),jj=1,NMOLE),
     &      (out(jj),jj=1,NDUST),
     &      (LOG10(MAX(1.Q-300, eldust(jj))),jj=1,NDUST),
     &      (LOG10(eps(elnum(jj))),jj=1,el-1),
     &      (LOG10(eps(elnum(jj))),jj=el+1,NELM),
     &       LOG10(MAX(1.Q-300, rhod/rhog)),
     &       LOG10(MAX(1.Q-300, dustV)),
     &       LOG10(MAX(1.Q-300, Jstar)), 
     &       MIN(999999.99999,Nstar)










!------------ START OF MODIFICATION TO DEMO_SWEEP-----------!

        ntotat1=0
        ntotat2=0
        ntotmol=0
        do k=1,el-1
          ntotat1=ntotat1+nat(elnum(k))
        enddo
        do k=el+1,NELM
          ntotat2=ntotat2+nat(elnum(k))
        enddo
        do k=1,NMOLE
          ntotmol=ntotmol+nmol(k)
        enddo

        ntot=nel+ntotat1+ntotat2+ntotmol

        ppel=(nel/ntot)*pges
        ppat=(nat/(ntot))*pges
        ppmol=(nmol/ntot)*pges

!        Finding indices for molecules 

         do k=1,NMOLE
           if (trim(cmol(k))=="N2")    iN2 = k
           if (trim(cmol(k))=="CH")    iCH = k
           if (trim(cmol(k))=="CO")    iCO = k
           if (trim(cmol(k))=="CN")    iCN = k
           if (trim(cmol(k))=="C2")    iC2 = k
           if (trim(cmol(k))=="NO")    iNO = k
           if (trim(cmol(k))=="C2H2")  iC2H2 = k
           if (trim(cmol(k))=="HCN")   iHCN = k
           if (trim(cmol(k))=="C2H")   iC2H = k
           if (trim(cmol(k))=="C3")    iC3 = k
           if (trim(cmol(k))=="CS")    iCS = k
           if (trim(cmol(k))=="H2O")   iH2O = k
           if (trim(cmol(k))=="OH")    iOH = k
           if (trim(cmol(k))=="TIO")  iTiO = k
           if (trim(cmol(k))=="SIO")   iSiO= k
           if (trim(cmol(k))=="CH4")   iCH4 = k
           if (trim(cmol(k))=="NH")    iNH = k
           if (trim(cmol(k))=="SIH")   iSiH = k
           if (trim(cmol(k))=="FEH")   iFeH = k
           if (trim(cmol(k))=="VO")    iVO = k
           if (trim(cmol(k))=="ZRO")   iZrO = k
           if (trim(cmol(k))=="MGH")   iMgH = k
           if (trim(cmol(k))=="NH3")   iNH3 = k
           if (trim(cmol(k))=="CO2")   iCO2 = k
           if (trim(cmol(k))=="TIH")   iTiH = k
           if (trim(cmol(k))=="CAH")   iCaH = k
           if (trim(cmol(k))=="CRH")   iCrH = k
           if (trim(cmol(k))=="LIH")   iLiH = k
           if (trim(elnam(k))=="H")      iH = k
           if (trim(cmol(k))=="O2")     iO2 = k
           if (trim(cmol(k))=="H-")     iHm = k
           if (trim(cmol(k))=="H2")     iH2 = k
           if (trim(cmol(k))=="H2+")   iH2p = k
           if (trim(cmol(k))=="HS")     iHS = k
           if (trim(cmol(k))=="C3H")   iC3H = k
           if (trim(cmol(k))=="SIC")   iSiC = k
           if (trim(cmol(k))=="SIC2") iSiC2 = k
           if (trim(cmol(k))=="NS")     iNS = k
           if (trim(cmol(k))=="SIN")   iSiN = k
           if (trim(cmol(k))=="SO")     iSO = k
           if (trim(cmol(k))=="S2")     iS2 = k
           if (trim(cmol(k))=="SIS")   iSiS = k
           if (trim(cmol(k))=="LAO")   iLaO = k
           if (trim(cmol(k))=="CH2")   iCH2 = k
           if (trim(cmol(k))=="CH3")   iCH3 = k
           if (trim(cmol(k))=="SI2C") iSi2C = k
           if (trim(cmol(k))=="SIO2") iSiO2 = k
           if (trim(cmol(k))=="H2S")   iH2S = k
           if (trim(cmol(k))=="CAOH") iCaOH = k
           if (trim(cmol(k))=="CHNO") iCHNO = k
           if (trim(cmol(k))=="SIF2") iSiF2 = k
         enddo



!     Writing partial pressures for molecules and saving in common block
!    /partialpressures/

        ppN2  = ppmol(iN2)
        ppCH  = ppmol(iCH)
        ppCO  = ppmol(iCO)
        ppCN  = ppmol(iCN)
        ppC2  = ppmol(iC2)
        ppNO  = ppmol(iNO)
        ppC2H2 = ppmol(iC2H2)
        ppHCN = ppmol(iHCN)
        ppC2H = ppmol(iC2H)
        ppC3  = ppmol(iC3)
        ppCS  = ppmol(iCS)
        ppH2O = ppmol(iH2O)
        ppOH  = ppmol(iOH)
        ppTiO = ppmol(iTiO)
        ppSiO = ppmol(iSiO)
        ppCH4 = ppmol(iCH4)
        ppNH  = ppmol(iNH)
        ppSiH = ppmol(iSiH)
        ppFeH = ppmol(iFeH)
        ppVO  = ppmol(iVO)
        ppZrO = ppmol(iZrO)
        ppMgH = ppmol(iMgH)
        ppNH3 = ppmol(iNH3)
        ppCO2 = ppmol(iCO2)
        ppTiH = ppmol(iTiH)
        ppCaH = ppmol(iCaH)
        ppCrH = ppmol(iCrH)
        ppLiH = ppmol(iLiH)
        ppH   = ppat(iH)
        ppO2  = ppmol(iO2)
        ppHm  = ppmol(iHm)
        ppH2  = ppmol(iH2)
        ppH2p = ppmol(iH2p)
        ppHS  = ppmol(iHS)
        ppC3H  = ppmol(iC3H)
        ppSiC  = ppmol(iSiC)
        ppSiC2  = ppmol(iSiC2)
        ppNS  = ppmol(iNS)
        ppSiN  = ppmol(iSiN)
        ppSO  = ppmol(iSO)
        ppS2  = ppmol(iS2)
        ppSiS  = ppmol(iSiS)
        ppCaH  = ppmol(iCaH)
        ppLaO  = ppmol(iLaO)
        ppCH2  = ppmol(iCH2)
        ppCH3  = ppmol(iCH3)
        ppSi2C  = ppmol(iSi2C)
        ppSiO2  = ppmol(iSiO2)
        ppH2S  = ppmol(iH2S)
        ppCaOH  = ppmol(iCaOH)
        ppCHNO  = ppmol(iCHNO)
        ppSiF2  = ppmol(iSiF2)



      open(unit=778,file='PartialPressures.dat',status='replace')
      write(778,*)ppN2,ppCH,ppCO,ppCN,ppC2,ppNO,ppC2H2,ppHCN,ppC2H,ppC3,
     & ppCS,ppH2O,ppOH,ppTiO,ppSiO,ppCH4,ppNH,ppSiH,ppFeH,ppVO,
     & ppZrO,ppMgH,ppNH3,ppCO2,ppTiH,ppCaH,ppCrH,ppLiH,ppH,ppO2,
     & ppHm,ppH2,ppH2p,ppHS,ppC3H,ppSiC,ppSiC2,ppNS,ppSiN,ppSO,
     & ppS2,ppSiS,ppLaO,ppCH2,ppCH3,ppSi2C,ppSiO2,ppH2S,ppCaOH,
     & ppCHNO,ppSiF2,
!      close(778)

     & iN2,iCH,iCO,iCN,iC2,iNO,iC2H2,iHCN,iC2H,iC3,iCS,iH2O,
     & iOH,iTiO,iSiO,iCH4,iNH,iSiH,iFeH,iVO,iZrO,iMgH,iNH3,
     & iCO2,iTiH,iCaH,iCrH,iLiH,iH,iO2,iHm,iH2,iH2p,iHS,
     & iC3H,iSiC,iSiC2,iNS,iSiN,iSO,iS2,iSiS,iLaO,iCH2,
     & iCH3,iSi2C,iSiO2,iH2S,iCaOH,iCHNO,iSiF2
!     & 'N2 ','CH ','CO ','CN ','C2 ','NO ','C2H2 ','HCN ','C2H ','C3 ',
!     & 'CS ','H2O ','OH ','TiO ','SiO ','CH4 ','NH ','SiH ','FeH ','VO ',
!     & 'ZrO ','MgH ','NH3 ','CO2 ','TiH ','CaH ','CrH ','LiH ','H ',
!     & 'O2 ','Hm ','H2 ','H2p ','HS ','C3H ','SiC ','SiC2 ','NS ','SiN ',
!     & 'SO ','S2 ','SiS ','LaO ','CH2 ','CH3 ','Si2C ','SiO2 ','H2S ',
!     & 'CaOH ','CHNO ','SiF2'        
       close(778)


       open(unit=709, file='runetest.dat',status='replace')
       write(709,*)
     & iFeH,cmol
       close(709)


       open(unit=990,file='GGchem_ppel',status='replace')
       write(990,*) ppel
       close(990)


        if (TEAinterface) then
          nTEA = 0.0                           ! no electrons
          do j=1,NELEM
            nTEA = nTEA + nat(j)               ! atoms
          enddo
          do j=1,NMOLE
            if (index(cmol(j),"+")>0) cycle    ! no ions
            if (index(cmol(j),"-")>0) cycle    ! no cations
            nTEA = nTEA + nmol(j)              ! molecules
          enddo
          pTEA = nTEA*bk*Tg
          write(71,'(1pE10.4,0pF8.2,99(1pE18.11))') pTEA/bar,Tg,
     &        (eps(elnum(jj)),jj=1,el-1),
     &        (eps(elnum(jj)),jj=el+1,NELM)
        endif

        if (verbose>0) read(*,'(a1)') char

      enddo

      close(70)
      close(71)

      !write(*,*)
      if (TEAinterface) write(*,frmt) trim(species)

 1000 format(4(' eps(',a2,') = ',1pD8.2))
 1010 format(A4,0pF8.2,3(a6,1pE9.2),1(a11,1pE9.2))
 2000 format(9999(1x,A19))
 2010 format(0pF20.6,2(1pE20.6),9999(0pF20.7))








      end


   
