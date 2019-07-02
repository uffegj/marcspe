*********************************************************************
      REAL*8 FUNCTION JEVAP(chi,Vlow,L0,L1,L2,L3)
*********************************************************************
***  Approximate evaporation rate Jev at Vlow [cm^-3/s].          ***
***  The idea is to choose Jev such that all (some) dust moments  ***
***  become zero after the same evaporation time t, if chi stays  ***
***  constant.                                                    ***
*********************************************************************
      implicit none
      real*8,intent(IN) :: chi,Vlow,L0,L1,L2,L3
      real*8,save :: a
      real*8 :: t,dt,func,dfunc,c0,c1,c2,c3,A0,A1,A2
      real*8 :: VIETA1,VIETA2,ppp,qqq
      integer :: it
      logical,save :: firstCall=.true.
*-----------------------------------------------------------------------
*  Die Formelfunktion zur Loesung quadratische Gleichungen mit Vieta
      VIETA1(ppp,qqq) = qqq/(-ppp/2.d0-DSQRT(ppp**2/4.d0-qqq))
      VIETA2(ppp,qqq) = qqq/(-ppp/2.d0+DSQRT(ppp**2/4.d0-qqq))
*-----------------------------------------------------------------------
      JEVAP = 0.d0
      return 
      
      if (firstCall) then
        a = Vlow**(1.d0/3.d0)  
        firstCall = .false.
      endif  

      if (chi.gt.0.d0) RETURN
      if (L0.le.0.d0) RETURN
      if ((L1.le.a*L0).or.(L2.le.a*L1).or.(L3.le.a*L2)) then
         !write(*,*) "unphysical moments in JEVAP" 
         !write(*,*) L0,L1,L2,L3
         !stop
         RETURN
      endif

      !--- L0(t)=0 and L1(t)=0 ---
      t = 6.d0*(a-L1/L0)/chi  
      JEVAP = -L0/t
      !write(*,*) t,JEVAP
      !write(*,*) L0,L0+JEVAP*t
      !write(*,*) L1,L1+JEVAP*a*t+1.0/3.0*chi*L0*t+1.0/6.0*chi*JEVAP*t**2
      RETURN

      !--- L1(t)=0 and L2(t)=0 ---
      c0  = 27.d0*L2 - 27.d0*a**2*L0 
      c1  = 18.d0*L1 -  9.d0*a*L0
      c2  =  2.d0*L0 
      ppp = c1/c0
      qqq = c2/c0
      t   = VIETA2(ppp,qqq)
      t   = VIETA1(ppp,qqq)
      t   = t/chi  
      JEVAP = -L0/t
      write(*,*) c0,c1,c2
      write(*,*) ppp,qqq,ppp**2/4.0-qqq
      write(*,*) t,JEVAP
      write(*,*) L0,L0 + JEVAP*t
      write(*,*) L1,L1 + JEVAP*a*t + chi*L0*t/3.0 + chi*JEVAP*t**2/6.0
      write(*,*) L2,L2 + JEVAP*a**2*t + chi**2*L0*t**2/9.0 
     &                 + chi*a*JEVAP*t**2/3.0 + chi**2*JEVAP*t**3/27.0
     &                 + chi*L1*t*2.0/3.0
      stop
      RETURN

      !--- L0(t)=0 and L3(t)=0 ---
      c0 = 36.d0*L3 - 36.d0*a**3*L0 
      c1 = 36.d0*L2 - 18.d0*a**2*L0
      c2 = 12.d0*L1 -  4.d0*a*L0 
      c3 = L0
      t  = 6.d0*(a-L1/L0)  
      t  = t*3.d0
      do it=1,100          
        func  = c0 + c1*t + c2*t**2 + c3*t**3
        dfunc = c1 + c2*2.d0*t + c3*3.d0*t**2
        dt = -func/dfunc
        !write(*,*) it,t,dt,func
        t  = t+dt
        if (ABS(dt).lt.1.e-13*ABS(t)) exit
      enddo  
      t = t/chi  
      JEVAP = -L0/t

      !--- L0(t)=0 and L3(t)=0 ---
      !c0 = -540.0*a**3*L0 + 1620.0*a**2*L1 - 1620.0*a*L2 + 540.0*L3 
      !c1 = +270.0*a**2*L0 + 270.0*L2 - 540.0*a*L1
      !c2 = +36.0*L1 - 36.0*a*L0
      !c3 = L0
      !A0 =  -9 * (-12*chi*t*a*L0 + chi**2*t**2*L0 + 12*chi*t*L1 &
      !   &        +30*L2 + 30*a**2*L0 -  60*a*L1) / chi**3 / t**3
      !A1 =  36 * (-16*chi*t*a*L0 + chi**2*t**2*L0 + 16*chi*t*L1 &
      !   &            +45*L2 + 45*a**2*L0 -  90*a*L1) / chi**3 / t**4
      !A2 = -30 * (-18*chi*t*a*L0 + chi**2*t**2*L0 + 18*chi*t*L1 &
      !   &          +54*L2 + 54*a**2*L0 - 108*a*L1) / chi**3 / t**5
      !JEVAP = -A0

      if ((t.lt.0.0).or.(JEVAP.gt.0.0)) then
        write(*,*) "cj=",c0,c1,c2,c3
        !write(*,*) "Aj=",A0,A1,A2 
        write(*,*) "it,t,Jev=",it,t,JEVAP
      endif    

      RETURN
      end 
