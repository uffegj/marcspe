*
*
***********************************************************************
      SUBROUTINE FDJAC(info,n,x,fvec,np,df) 
***********************************************************************
      LOGICAL info
      INTEGER n,np,NMAX
      REAL*8 df(np,np),fvec(n),x(n),EPS 
      PARAMETER (NMAX=30,EPS=1.D-5) 
      INTEGER i,j 
      REAL*8 h,temp,f1(NP),f2(NP) 

c     write(*,*) "FDJAC:",x
      do j=1,n 
        temp=x(j) 
c       h=EPS
        h=DMAX1(1.d-18,temp*EPS)
        if (temp-h.lt.0.d0) h=0.5*temp
        if (temp.eq.0.d0) stop 'PANIC: x=0 in FDJAC'
        x(j)=temp+h
        call FUNCV(.false.,n,x,f2) 
        x(j)=temp-h
        call FUNCV(.false.,n,x,f1) 
        x(j)=temp 
        do i=1,n 
          df(i,j)=(f2(i)-f1(i))/(2.D0*h) 
          if (info) write(*,*) i,j,df(i,j)
        enddo 
      enddo 
      RETURN 
      END 
