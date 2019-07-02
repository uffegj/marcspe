************************************************************************
      SUBROUTINE FUNCV(info,N,epsi,fvec)
************************************************************************
      use drift_data,ONLY: NEPS
      implicit none
      logical,intent(IN):: info
      integer,intent(IN):: N
      real*8,intent(IN) :: epsi(N)
      real*8,intent(OUT):: fvec(N)
      CALL ELFUNK(epsi,1,0.d0,fvec)
      if (info) then
        write(*,*)  
        write(*,*) epsi(1:NEPS)
        write(*,*) fvec(1:NEPS)
      endif  
      RETURN
      end
