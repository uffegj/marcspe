**********************************************************************
      SUBROUTINE KOPF
**********************************************************************
*****                                                            *****
*****  offnet die Files                                          *****
*****                                                            *****
**********************************************************************
      use drift_data,ONLY: NELEM,NDUST,NEPS,NSPECIES,NNUC,
     &                     sizedist,abschluss,beta,restart,
     &                     bar,Teff,logg,mixLength,pconv,Nl,Vl,
     &                     elnam,spnam,dust_nam,nuc_nam,
     &                     eps0,elnr,NMOLE,cmol
      implicit none
      integer :: i,Nelement
      if (restart) then
        open(unit=91, file='out3_thermo.dat', position='APPEND')
        open(unit=92, file='out3_nuclea.dat', position='APPEND')
        open(unit=93, file='out3_dust.dat'  , position='APPEND')
        open(unit=94, file='out3_imp.dat'   , position='APPEND')
        open(unit=95, file='2Phoenix.data'  , position='APPEND')
        open(unit=96, file='out3_dist.dat'  , position='APPEND')
        open(unit=97, file='out3_chem1.dat' , position='APPEND')
        open(unit=98, file='out3_chem2.dat' , position='APPEND')
        open(unit=99, file='out3_chem3.dat' , position='APPEND')
        open(unit=100,file='drift2marcs.dat', position='APPEND')
      else    
        open(unit=91, file='out3_thermo.dat', status='replace')
        open(unit=92, file='out3_nuclea.dat', status='replace')
        open(unit=93, file='out3_dust.dat'  , status='replace')
        open(unit=94, file='out3_imp.dat'   , status='replace')
        open(unit=95, file='2Phoenix.data'  , status='replace')
        open(unit=96, file='out3_dist.dat'  , status='replace')
        open(unit=97, file='out3_chem1.dat' , status='replace')
        open(unit=98, file='out3_chem2.dat' , status='replace')
        open(unit=99, file='out3_chem3.dat' , status='replace')
        open(unit=100,file='drift2marcs.dat', status='replace')
        write(91,1010) Teff,logg,mixLength,abschluss,beta,pconv/bar
        write(92,1010) Teff,logg,mixLength,abschluss,beta,pconv/bar
        write(93,1010) Teff,logg,mixLength,abschluss,beta,pconv/bar
        write(94,1010) Teff,logg,mixLength,abschluss,beta,pconv/bar
        write(96,1010) Teff,logg,mixLength,abschluss,beta,pconv/bar
        write(97,1010) Teff,logg,mixLength,abschluss,beta,pconv/bar
        write(98,1010) Teff,logg,mixLength,abschluss,beta,pconv/bar
        write(91,*) NEPS,NDUST
        write(92,*) NEPS,NDUST
        write(93,*) NEPS,NDUST
        write(94,*) NEPS,NDUST
        write(96,*) NEPS,NDUST
        write(97,*) NEPS,NMOLE
        write(98,*) NEPS,NMOLE
        write(91,*)
        write(92,*)
        write(93,*)
        write(94,*)
        write(96,*)
        write(97,*)
        write(98,*)
        write(91,1000) 'z','T','n<H>','rho','p', 'wmix','tau ', 
     &   ('n_'//spnam(i),i=1,NSPECIES), 'pel'
        write(92,1000) 'z','p',('Seff_'//dust_nam(i),i=1,NDUST), 
     &   'N*_TiO2','J*_TiO2'
        write(93,1000) 'z','p','L0','L1','L2','L3','<a>','<vdr>', 
     &   ('b_'//dust_nam(i),i=1,NDUST), 
     &   ('c_'//dust_nam(i),i=1,NDUST), 
     &   ('fcond_'//elnam(elnr(i)),i=1,NEPS), 
     &   ('eps_'//elnam(elnr(i)),i=1,NEPS),
     &   'rhod/rho', 'rhod', 'chinet'
        write(94,1000) 'z', 'p', 'nuklea', 'growth', 'drift', 'misch',
     &                 'rho','T','aquer[mic]','vdr[m/s]','Kn',
     &                 'tau_nucl','tau_grow','tau_sink','tau_mix'
        if (sizedist.eq.1) then
          write(96,1000) 'p', 'pel', 'T', 'a1', 'a2', 'N1', 'N2'
        else if (sizedist.eq.2) then
          write(96,1000) 'p', 'pel', 'T', 'A', 'B', 'C'
        else if (sizedist.eq.3) then
          write(96,1000) 'p', 'pel', 'T', 'N', 'a1', 'sigma', 'sol-Nr.'
        endif  
        write(97,1000) 'z','p','T','n<H>',cmol(1:100)
        write(98,1000) 'z','p','T','n<H>',cmol(101:NMOLE)
        write(99,1000) 'z','p','T','n<H>','He','H','C','N','O','Si',
     & 'Mg','Al','Fe','S','Na','K','Ti','Ca'
        write(*,1001) 'p', 'dz', 'T', 'J*', 'chinet', '<a>', 'fcond'
*       ---------------------------
*       ***  PHOENIX-Interface  ***
*       ---------------------------
 !       write(95,*)
 !       write(95,*) 'Stellar Parameter:'
 !       write(95,*) '=================='
 !       write(95,1011) Teff,logg,mixLength,beta
 !       Nelement = 0
 !       do i=1,NELEM
 !         if (elnam(i).ne.'  ') Nelement=Nelement+1
 !       enddo 
 !       write(95,'(i4," deep element abundances")') Nelement
 !       do i=1,NELEM
 !         if (elnam(i).ne.'  ') write(95,*) i,elnam(i),eps0(i)
 !       enddo  
 !       write(95,*)
 !       write(95,*) 'Dust Parameter:'
 !       write(95,*) '==============='
 !       write(95,2010) Nl,Vl,abschluss
 !       write(95,'(i4," nucleation species")') NNUC
 !       do i=1,NNUC
 !         write(95,*) nuc_nam(i)
 !       enddo   
 !       write(95,'(i4," dust species")') NDUST
 !       do i=1,NDUST
 !         write(95,*) dust_nam(i)
 !       enddo
 !       write(95,'(i4," affected elements")') NEPS
 !       do i=1,NEPS
 !         write(95,*) elnr(i),elnam(elnr(i))
 !       enddo
 !       write(95,*)
 !       write(95,2000) 'R','T','rho','p','L0','L1','L2','L3',
 !    &   ('Vol%_'//dust_nam(i),i=1,NDUST), 
 !    &   ('eps_'//elnam(elnr(i)),i=1,NEPS)   
     
*       ---------------------------
*       ***  MARCS-Interface  ***
*       ---------------------------
        write(100,*)
        write(100,*) 'Stellar Parameters:'
        write(100,*) '=================='
        write(100,1011) Teff,logg,mixLength,beta
        Nelement = 0
        do i=1,NELEM
          if (elnam(i).ne.'  ') Nelement=Nelement+1
        enddo 
        write(100,'(i4," deep element abundances")') Nelement
        do i=1,NELEM
          if (elnam(i).ne.'  ') write(100,*) i,elnam(i),eps0(i)
        enddo  
        write(100,*)
        write(100,*) 'Dust Parameters:'
        write(100,*) '==============='
        write(100,2010) Nl,Vl,abschluss
        write(100,'(i4," nucleation species")') NNUC
        do i=1,NNUC
          write(100,*) nuc_nam(i)
        enddo   
        write(100,'(i4," dust species")') NDUST
        do i=1,NDUST
          write(100,*) dust_nam(i)
        enddo
        write(100,'(i4," affected elements")') NEPS
        do i=1,NEPS
          write(100,*) elnr(i),elnam(elnr(i))
        enddo
        write(100,*)
        write(100,2000) 'R','T','rho','p','L0','L1','L2','L3',
     &   ('Vol%_'//dust_nam(i),i=1,NDUST), 
     &   ('eps_'//elnam(elnr(i)),i=1,NEPS)   
      endif  

      RETURN
 1000 format(999(A12))
 1001 format(99(A11))
 1010 format(' Teff = ',0pF8.2,'  logg = ',0pF5.2,
     &       '  mixLengthPara = ',0pF5.2,'  Abschluss = ',i2,
     &       '  beta = ',0pF5.2,'  pconv[bar] = ',0pF5.2)
 1011 format(' Teff = ',0pF8.2,'  logg = ',0pF5.2,
     &       '  mixLengthPara = ',0pF5.2,'  OvershotPara = ',0pF5.2)
 2000 format(99(A20))
 2010 format(' Nl = ',i5,'  Vl = ',1pE15.6,'  Abschluss = ',i2)
      end





