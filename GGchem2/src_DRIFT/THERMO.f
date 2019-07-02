**********************************************************************
      SUBROUTINE THERMO(z,p,T,nH,rho,g,mu,wmix,welche)
**********************************************************************
      use drift_data,ONLY: Rnull,muH,Nlayers,zlay,
     &                     Rlay,Tlay,play,rholay,glay,mulay,wmixlay
      implicit none
      real*8,intent(IN)  :: z
      integer,intent(IN) :: welche
      real*8,intent(OUT) :: p,T,nH,rho,g,mu,wmix
      real*8 :: r,rfak1,rfak2
      integer,save :: i=2 

      if ((welche.eq.1) .or. (welche.eq.4)) then               ! PHOENIX/MARCS
        r = Rnull-z
        do while ((r.lt.Rlay(i+1)).and.(i.lt.Nlayers-1))
          i = i+1
        enddo   
        do while ((r.ge.Rlay(i)).and.(i.gt.2))
          i = i-1
        enddo  
        rfak1 = (r-Rlay(i+1))/(Rlay(i)-Rlay(i+1))
c       write(*,*) i,Rnull-Rlay(i),z,Rnull-Rlay(i+1),rfak1
      else if ((welche.eq.2).or.(welche.eq.3)) then        ! AMES/TestCases
        r = Rnull+z                 
        do while ((r.gt.zlay(i+1)).and.(i.lt.Nlayers-1))   
          i = i+1
        enddo  
        do while ((r.le.zlay(i)).and.(i.gt.2))
          i = i-1
        enddo  
        rfak1 = (r-zlay(i+1))/(zlay(i)-zlay(i+1))
c       write(*,*) "in THERMO i,zlay(i),r,zlay(i+1)=", 
c    &             i,zlay(i),r,zlay(i+1),rfak1
      endif     
      rfak2 = 1.d0-rfak1 
*
*     -------------------------------
*     ***  linear interpolations  ***
*     -------------------------------
      T  =  Tlay(i)*rfak1 +  Tlay(i+1)*rfak2
      mu = mulay(i)*rfak1 + mulay(i+1)*rfak2
      g  =  glay(i)*rfak1 +  glay(i+1)*rfak2
      
*     ----------------------------
*     ***  log interpolations  ***
*     ----------------------------
      p    = DEXP(DLOG(   play(i))*rfak1 + DLOG(   play(i+1))*rfak2) 
      rho  = DEXP(DLOG( rholay(i))*rfak1 + DLOG( rholay(i+1))*rfak2) 
      wmix = DEXP(DLOG(wmixlay(i))*rfak1 + DLOG(wmixlay(i+1))*rfak2)      
      nH   = rho/max(1.d-99,muH)     ! sonst gibt es ein NaN waehrend INIT   
*
      RETURN
      end
