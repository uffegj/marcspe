**********************************************************************
      SUBROUTINE DESTIFF(nH,T,bmix)
**********************************************************************
*****                                                            *****
*****  Entsteifung der Verdampfung                               *****
*****                                                            *****
**********************************************************************
      use drift_data,ONLY: NNUC,NEPS,NDUST,
     &                     eps,elnr,dust_nam,elnam,pred
      implicit none
      real*8,intent(IN) :: nH,T,bmix(NDUST)
      real*8,dimension(NDUST) :: hSat,hbmix,hchi,hchi1,hchi2
      real*8 :: effSat(NDUST),epsi(NDUST),dchideps(NDUST)
      real*8 :: FF(NEPS),hJst(NNUC),hchinet,dep
      real*8,parameter :: DCHIDEPS_CRIT=0.5,SEFF_CRIT=1.001
      integer :: i,j
      logical :: changed
      COMMON /OUTHELP/hSat,hbmix,hchi,hchinet,hJst

      call EFF_SUPERSAT(nH,T,bmix,effSat)
      do j=1,NEPS
        epsi(j)=eps(elnr(j)) 
      enddo
      do
        changed=.false.
        dchideps=0.0
        do j=1,NEPS
          dep = DMAX1(1.d-22,epsi(j)*1.e-8)
          if (epsi(j)-20.0*dep.lt.0.d0) dep=0.05*epsi(j)
          CALL ELFUNK(epsi,j,+dep,FF)
          hchi1=hchi
          CALL ELFUNK(epsi,j,-dep,FF)
          hchi2=hchi
          do i=1,NDUST
            dchideps(i)=MAX(dchideps(i),(hchi1(i)-hchi2(i))/dep)
          enddo   
        enddo
        do i=1,NDUST
          if ((effSat(i).lt.SEFF_CRIT).and.
     &        (dchideps(i).gt.DCHIDEPS_CRIT)) then 
            changed = .true.
            pred(i) = 0.9*pred(i)           
          endif  
        enddo
        if (.not.changed) exit
        write(*,'(5x,99(1x,a7))') dust_nam(1:NDUST)  
        write(*,'(a5,99(1pE8.1))') "bmix",(bmix(i),i=1,NDUST)
        write(*,'(a5,99(1pE8.1))') "Seff",(effSat(i),i=1,NDUST)
        write(*,'(a5,99(1pE8.1))') "pred",pred(1:NDUST)
        write(*,'(5x,99(1x,a7))') dust_nam(1:NDUST)  
        do j=1,NEPS
          dep = DMAX1(1.d-22,epsi(j)*1.e-8)
          if (epsi(j)-20.0*dep.lt.0.d0) dep=0.05*epsi(j)
          CALL ELFUNK(epsi,j,+dep,FF)
          hchi1=hchi
          CALL ELFUNK(epsi,j,-dep,FF)
          hchi2=hchi
          write(*,'(a5,99(0pF8.4))') elnam(elnr(j)),
     &             ((hchi2(i)-hchi1(i))/dep,i=1,NDUST)
        enddo
      enddo  
      END
