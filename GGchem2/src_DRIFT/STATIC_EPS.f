**********************************************************************
      SUBROUTINE STATIC_EPS(epsi,iflag)
**********************************************************************
      use drift_data,ONLY: NEPS,NDUST,T,nH,pel,rhoL2,wmix,L4,
     &                     eps0,epsmerk,elnr,elnam,nat,nion,nmol,
     &                     NNUC,nuc_nel,nuc_el,elcode,dust_nam,
     &                     chinet,chi,bmix,Sat,
     &                     Jst,Nst,firstStaticEpsCall,
     &                     J_is_zero,no_converge
      implicit none
      real*8,intent(OUT) :: epsi(NEPS)
      integer,intent(IN) :: iflag
      real*8 :: dep, fak, del, old, new, emin, emax
      real*8 :: FF0(NEPS), FF1(NEPS), FF2(NEPS), FF(NEPS), xnew(NEPS)
      real*8 :: deps(NEPS), DF(NEPS,NEPS), DF0(NEPS,NEPS)
      real*8 :: delta,CHI_NET,Jeff,f1,f2,qual
      real*8,parameter :: DELMAX=1.d-13, EVAR=1.d-8
      real*8 :: work(NEPS*(NEPS+1))
      character(len=1) :: char1
      logical :: info,IS_NAN,changed
      logical,save :: nucleates(NEPS)
      integer :: ind,indx(NEPS)
      integer :: it,itmax,i,j,el

      if (firstStaticEpsCall) then
        do i=1,NEPS
          epsi(i) = 0.999d0*eps0(elnr(i))
        enddo
        nucleates = .false.
c       -----------------------------------------------------------------
c       *** reduce nucleating elements to get suitable initial values ***
c       -----------------------------------------------------------------
        do i=1,NNUC
          do j=1,nuc_nel(i)
            el = elcode(nuc_el(i,j))
            if (elnam(elnr(el)).ne.'O') nucleates(el)=.true.
          enddo
        enddo
        do
          CALL ELFUNK(epsi,1,0.d0,FF)
c         write(*,*) "eps=",epsi
c         write(*,*) " FF=",FF
          changed = .false.
          do i=1,NEPS
            if (nucleates(i).and.(FF(i).lt.-0.9d0)) then
              epsi(i) = epsi(i)/10.d0
              changed = .true.
            endif
            epsmerk(i) = epsi(i)
          enddo 
          if (.not.changed) exit 
        enddo
        itmax = 10000
        firstStaticEpsCall = .false.
      else
        do i=1,NEPS
          epsi(i) = epsmerk(i)
        enddo
        itmax = 30
      endif

      if (no_converge) return

      it   = 0
      info = .false.
      J_is_zero = .false.
      Jeff = 0.0
      CALL GGCHEM(nH,T,pel)
      CALL SUPERSAT(nH,T,Sat)
      CALL NUCLEATION(T,Jst,Nst)
      chinet = CHI_NET(T,Sat,L4,chi,bmix)
      do i=1, NNUC
       Jeff = Jeff + Jst(i)
      enddo
      if (Jeff.lt.1.d-70) J_is_zero = .true.
c     if (chinet.lt.0.d0) J_is_zero = .true.
      if (info) write(*,*) "iflag=",iflag

      do
        !if(it.gt.10) info=.true.

c       ---------------------------------------
c       ***  determine Newton-Raphson step  ***
c       ---------------------------------------
        CALL ELFUNK(epsi,1,0.d0,FF)
        if (info) write(*,*) FF
        do j=1,NEPS
          dep = DMAX1(1.d-22,epsi(j)*EVAR)
          if (epsi(j)-20.0*dep.lt.0.d0) dep=0.05*epsi(j)
          CALL ELFUNK(epsi,j,+dep,FF1)
          CALL ELFUNK(epsi,j,-dep,FF2)
          do i=1,NEPS
c           if(info) write(*,*) j,i,FF1(i),FF2(i) 
            DF(i,j) = (FF2(i)-FF1(i))/(2.d0*dep)
          enddo
        enddo
        
        FF0 = FF
        DF0 = DF
        CALL SGEIR(DF,NEPS,NEPS,FF,1,ind,work,indx)
c A, LDA, N, V, ITASK, IND, WORK, IWORK

        deps = FF
        if (ind.eq.-4) then
          FF = FF0
          DF = DF0
          CALL GAUSS(NEPS,DF,deps,FF)
        endif  


        do i=1,NEPS
          if (IS_NAN(deps(i)).or.
     &        (DABS(deps(i)).gt.1.d+29*epsi(i))) then
            stop
            write(*,*) '*** WARNING NaN in STATIC_EPS:',
     &                 i,epsi(i),deps(i),FF0(i)/DF0(i,i)
            if (DF0(i,i).ne.0.d0) then
              deps(i) = FF0(i)/DF0(i,i)
            endif  
            if (IS_NAN(deps(i)).or.
     &        (DABS(deps(i)).gt.1.d+29*epsi(i))) then
              deps(i) = 0.999*(epsmerk(i)-epsi(i))
            endif  
c           write(*,*) "==> deps=",deps(i)
c           info = .true.
          endif
        enddo

        !---------------------------------------------
        ! ***  check solution of linear eq. system ***
        !---------------------------------------------
c       do i=1,NEPS
c         f1 = 0.0
c         do j=1,NEPS
c           f1 = f1 + DF0(i,j)*deps(j)
c         enddo
c         f2 = FF0(i)
c         qual = ABS(f1-f2)/ABS(f2)
c         write(*,*) "qual",i,qual
c       enddo  
 
c       -----------------------------------
c       ***  disable too large changes  ***
c       -----------------------------------
        fak = 1.0
        do i=1,NEPS
c         if (info.and.(iflag.ne.2)) then 
          if (info) then 
            write(*,*) elnam(elnr(i)), epsi(i), deps(i)
          endif
          del  = deps(i)
          old  = epsi(i)
          new  = epsi(i)+del
          emax = 3.0*old 
          emin = 0.3*old
          if (new.gt.emax) fak=DMIN1(fak,(emax-old)/del)
          if (new.lt.emin) fak=DMIN1(fak,(emin-old)/del)
        enddo
        if (info) write(*,*) 'fak=',fak
        deps = fak*deps

c       --------------------------------------------------------
c       ***  pullback NR step to achieve global convergence  ***
c       --------------------------------------------------------
        call PULLBACK(NEPS,epsi,deps,FF0,xnew,info)
        delta = 0.d0
        do i=1,NEPS
          delta = DMAX1(delta,DABS(xnew(i)/epsi(i)-1.d0))
          epsi(i) = xnew(i)
        enddo  
c       if (info.and.(iflag.ne.2)) then
        if (info) then
          write(*,*) it, delta
          read(*,'(a1)') char1
        endif
        it = it + 1
        if (it.ge.itmax) exit
        if ((fak.gt.0.999).and.(delta.lt.DELMAX)) exit

      enddo

      if (info) write(*,*) it,' Iterationen'
c     if (info.and.(iflag.ne.2)) read(*,'(a1)') char1
      if (info) read(*,'(a1)') char1
      if ((it.lt.itmax).and.(delta.lt.DELMAX)) then
        if (iflag.ne.2) then        ! call from Limex to determine Jacobian
          do i=1,NEPS
            epsmerk(i) = epsi(i)
c           write(*,1010) elnam(elnr(i)),epsi(i)/eps0(elnr(i))
          enddo
        endif
      else  
        write(*,*) '*** keine Konvergenz in STATIC_EPS!'
        write(*,*) '*** badness = ',delta
        no_converge = .true.
        do i=1,NEPS
          epsi(i) = epsmerk(i) 
        enddo
        call DESTIFF(nH,T,bmix) 
      endif  
*
      RETURN 
 1010 format(a5,3(1pE12.5))
      END 

!========================================================================
      subroutine PULLBACK(N,xx,dx,Fold,xnew,info)
!========================================================================
!     pull back Newton-Raphson solution, until quality improves
!     Qual(F(xx+fac*dx)) -> Min, fac<=1
!========================================================================
      use drift_data,ONLY: NDUST,NNUC
      implicit none
      integer,intent(in) :: N
      real,intent(in) :: xx(N),Fold(N)  
      real,intent(in) :: dx(N)
      real,intent(out) :: xnew(N)
      logical,intent(IN) :: info
      real :: Fnew(N),fac,qold,qnew
      real*8 :: Sat(NDUST), bmix(NDUST), chi(NDUST)
      real*8 :: Jst(NNUC), chinet
      integer,parameter :: itmax=20
      integer :: i,it
      COMMON /OUTHELP/Sat,bmix,chi,chinet,Jst

      qold=0.0
      do i=1,N
        qold=qold+Fold(i)**2
      enddo  
      xnew=xx+dx
      fac=1.0
      do it=1,itmax
        xnew = xx + fac*dx
        call FUNCV(.false.,N,xnew,Fnew)
        qnew = 0.0
        do i=1,N
          qnew = qnew + Fnew(i)**2
        enddo  
        if(info) then
          !write(*,'(99(1pE11.4))') Jst
          !write(*,'(99(1pE11.4))') bmix
          !write(*,'(99(1pE11.4))') chi
          !write(*,'(99(1pE11.4))') Sat
          !write(*,'(99(1pE11.4))') Fnew
        endif  
        if (info) write(*,*) it,qold,qnew
        if (qnew.le.qold) exit
        fac=fac*0.6
      enddo
c     if (it.ge.itmax) xnew=xx  ! diese Zeile bewirkt, dass in
                                ! STATIC_EPS die NR-Iteration abgebrochen
                                ! und xnew als Loesung akzeptiert wird
                                ! ==> sehr gefaehrlich
      return
      end      
