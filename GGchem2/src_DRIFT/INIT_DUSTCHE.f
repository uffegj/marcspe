**********************************************************************
      SUBROUTINE INIT_DUSTCHEM
**********************************************************************
      use drift_data,ONLY: NEPS,NELEM,NSPECIES,NDUST,NNUC,NREAC,NMOLE,
     &                     amu,dust_nel, dust_nu, dust_el,
     &                     elnr, elcode, spnr, 
     &                     nuc_nel, nuc_nu, nuc_el,
     &                     neduct, nprod, reac_nu, reac_sp,
     &                     eps, nel, nat, nion, nmol, muH, mass,
     &                     spmass, keysp, dust_rho, dust_vol,
     &                     cmol,elnam,spnam,dust_nam,nuc_nam,pred
      implicit none
      integer :: i,j,k,kk,nn,mm,stindex,sp,el,pel
      integer :: elem_cons(NELEM),spstoich(NSPECIES,10)
      integer :: spzahl(NSPECIES),spelem(NSPECIES,10)
      real*8  :: dmass
      character:: zeile*200
      character(len=10) :: name10
      character(len=2)  :: name2
      logical is_atom,found

      call GGCHEM(1.d+15, 2000.d0,pel)              ! damit cmol vorliegt
 
      write(*,*) 
      write(*,*) "reading DustChem12 data ..."
c      write(*,*) "reading DustChem_Crich data ..."
      write(*,*) "========================="

c      open(12, file='DustChem6.dat', status='old')
      open(12, file='data/DRIFT/DustChem7.dat', status='old')
c      open(12, file='DustChem12.dat', status='old')
c      open(12, file='DustChem15.dat', status='old')
c      open(12, file='DustChem_Crich.dat', status='old')
 
      write(*,*) '--- involved elements ---'
      read(12,1000) zeile
      read(12,1000) zeile
      read(12,*) i
      if (NEPS.ne.i) stop 'Parameter NEPS inconsistent.'
      do i=1,NEPS
         read(12,1010) name2
         found = .false. 
         do j=1,NELEM
            if (name2.eq.elnam(j)) then
               found = .true. 
               elnr(i) = j
               elcode(j) = i
            endif
         enddo
         if (.not.found) stop 'Element not found.'
         write(*,*) elcode(elnr(i)),' ',name2,elnr(i)
      enddo
 
      write(*,*) '--- involved gas species ---'
      read(12,1000) zeile
      read(12,1000) zeile
      read(12,1000) zeile
      read(12,*) i
      if (NSPECIES.lt.i) stop 'Parameter NSPECIES too small.'
      do i=1,NSPECIES
         read(12,1020) is_atom,keysp(i),mm,name10
         spnam(i)  = name10
         spnr(i)   = 0
         if (is_atom) then
            spzahl(i) = 1
            name2 = name10(1:2)
            found = .false.
            do j=1,NELEM
               if (name2.eq.elnam(j)) then
                  spelem(i,1) = j
                  spstoich(i,1) = 1
                  spnr(i) = j + 1000
                  spmass(i) = mass(j)
                  found = .true.
               endif
            enddo
            if (.not.found) stop 'Element of atomic species not found.'
         else
            spnr(i) = stindex(cmol,NMOLE,name10)
            spmass(i) = 0.d0
            spzahl(i) = mm
            do j=1,mm
               read(12,1030) nn,name2
               found = .false. 
               do kk=1,NELEM
                  if (elnam(kk).eq.name2) then
                     spelem(i,j)   = kk
                     spstoich(i,j) = nn 
                     spmass(i) = spmass(i) + DBLE(nn)*mass(kk)
                     found = .true.
                  endif   
               enddo
               if (.not.found) stop 'Element in molecule not found.'
            enddo
         endif
         if (spnr(i).eq.0) stop 'Molecular species not found.'
         write(*,1050) spnam(i),spnr(i),spmass(i)/amu
      enddo
 
      write(*,*) '--- dust species ---'
      read(12,1000) zeile
      read(12,1000) zeile
      read(12,1000) zeile
      read(12,*) i
      if (NDUST.ne.i) stop 'Parameter NDUST incorrect.'
      do i=1,NDUST
         read(12,*) dust_nam(i)
         read(12,*) dust_rho(i)
         read(12,*) dust_nel(i)
         dmass = 0.d0
         do j=1,dust_nel(i)
            read(12,1030) dust_nu(i,j),name2
            found = .false. 
            do kk=1,NELEM
               if (elnam(kk).eq.name2) then
                  dust_el(i,j) = kk
                  dmass = dmass + dust_nu(i,j)*mass(kk)
                  found = .true.
               endif
            enddo
            if (.not.found) stop 'Element in dust species not found'
         enddo
         dust_vol(i) = dmass/dust_rho(i)
         write(*,1060) dust_nam(i), dust_rho(i), dust_vol(i), 
     &      (dust_nu(i,j),elcode(dust_el(i,j)),j=1,dust_nel(i))
      enddo
      if (dust_nam(1).ne.'TiO2[s]   ') then
         stop 'TiO2[s] must be first dust species.'
      endif   

      write(*,*) '--- nucleation species ---'
      read(12,1000) zeile
      read(12,1000) zeile
      read(12,1000) zeile
      read(12,*) i
      if (NNUC.ne.i) stop 'Parameter NNUC incorrect.'
      do i=1,NNUC
         read(12,*) nuc_nam(i)
         read(12,*) nuc_nel(i)
         do j=1,nuc_nel(i)
            read(12,1030) nuc_nu(i,j),name2
            found = .false.
            do kk=1,NELEM
               if (elnam(kk).eq.name2) then
                  nuc_el(i,j) = kk
                  found = .true.
               endif
            enddo
            if (.not.found) stop 'Element for nuc. species not found'
         enddo
         write(*,1070) nuc_nam(i), 
     &      (nuc_nu(i,j),elcode(nuc_el(i,j)),j=1,nuc_nel(i))
      enddo
 
      write(*,*) '--- growth reactions ---'
      read(12,1000) zeile
      read(12,1000) zeile
      read(12,1000) zeile
      read(12,*) i
      if (NREAC.ne.i) stop 'Parameter NREAC incorrect.'
      do i=1,NREAC
         do j=1,NELEM
           elem_cons(j) = 0
         enddo  
         read(12,*) neduct(i), nprod(i)
         do j=1,neduct(i)
            read(12,1040) reac_nu(i,j,1),name10
            found = .false.
            do kk=1,NSPECIES
               if (name10.eq.spnam(kk)) then
                  reac_sp(i,j,1) = kk
                  found = .true.
               endif   
            enddo
            if (.not.found) then
              write(*,*) name10 
              stop 'educt in reaction not found'
            endif  
            sp = reac_sp(i,j,1)
            do k=1,spzahl(sp)
              el = spelem(sp,k)  
              elem_cons(el) = elem_cons(el) 
     &                      + reac_nu(i,j,1)*spstoich(sp,k)
            enddo  
         enddo
         read(12,1040) reac_nu(i,1,2),name10
         found = .false.
         do kk=1,NDUST
            if (name10.eq.dust_nam(kk)) then 
              reac_sp(i,1,2) = kk
              sp = reac_sp(i,1,2)
              found = .true.
              do k=1,dust_nel(sp)
                el = dust_el(sp,k)
                elem_cons(el) = elem_cons(el) 
     &                        - reac_nu(i,1,2)*dust_nu(sp,k)
              enddo  
            endif  
         enddo
         if (.not.found) then
            write(*,*) name10
            stop 'dust species in reaction not found'
         endif   
         do j=2,nprod(i)
            read(12,1040) reac_nu(i,j,2),name10
            found = .false.
            do kk=1,NSPECIES
               if (name10.eq.spnam(kk)) then
                 reac_sp(i,j,2) = kk
                 sp = reac_sp(i,j,2)
                 found = .true.
                 do k=1,spzahl(sp)
                   el = spelem(sp,k)  
                   elem_cons(el) = elem_cons(el) 
     &                           - reac_nu(i,j,2)*spstoich(sp,k)
                 enddo  
               endif  
            enddo
            if (.not.found) then
               write(*,*) name10
               stop 'product in reaction not found'
            endif   
         enddo
         if (neduct(i).eq.1) then
           write(*,2011) (reac_nu(i,j,1),spnam(reac_sp(i,j,1)),
     &     j=1,neduct(i)),reac_nu(i,1,2), dust_nam(reac_sp(i,1,2)), 
     &     (reac_nu(i,j,2),spnam(reac_sp(i,j,2)),j=2,nprod(i))
         elseif (neduct(i).eq.2) then
           write(*,2021) (reac_nu(i,j,1),spnam(reac_sp(i,j,1)),
     &     j=1,neduct(i)),reac_nu(i,1,2), dust_nam(reac_sp(i,1,2)), 
     &     (reac_nu(i,j,2),spnam(reac_sp(i,j,2)),j=2,nprod(i))
         elseif (neduct(i).eq.3) then
           write(*,2031) (reac_nu(i,j,1),spnam(reac_sp(i,j,1)),
     &     j=1,neduct(i)),reac_nu(i,1,2), dust_nam(reac_sp(i,1,2)), 
     &     (reac_nu(i,j,2),spnam(reac_sp(i,j,2)),j=2,nprod(i))
         endif
         do j=1,NELEM
            if (elem_cons(j).ne.0) then
               write(*,*) elnam(j),elem_cons(j) 
               stop 'element conservation in reaction violated.'
            endif  
         enddo    
      enddo
      close(12)
      write(*,*)
*
      !-------------------------------------------------------------
      ! ***  init Potenzen zur "Entsteifung der Uebersaettigung" ***
      !-------------------------------------------------------------
      do j=1,NDUST
        pred(j)=1.0
      enddo  
*
      RETURN 
 1000 format(a200)
 1010 format(a2)
 1020 format(2(l1,1x),i1,1x,a10)
 1030 format(i2,1x,a2)
 1040 format(i2,1x,a10)
 1050 format(1x,a10,i4,' mass=',0pf7.3," amu")
 1060 format(1x,a10," rhod=",0pf6.3," V0=",1pe11.3,2x,99(i1,"x",i2,1x))
 1070 format(1x,a10,99(i1,"x",i2,1x))
 2011 format(1(I2,1x,a8),22x,'->',I2,1x,a10,99(I2,1x,a8))
 2021 format(2(I2,1x,a8),11x,'->',I2,1x,a10,99(I2,1x,a8))
 2031 format(3(I2,1x,a8)    ,'->',I2,1x,a10,99(I2,1x,a8))
      end 
