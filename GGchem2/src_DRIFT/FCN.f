
!=========================================================================
      subroutine FCN ( N, Nz, t, y, FUNC, BB, ir, ic, iflag )
!=========================================================================
!*****                                                               *****  
!*****  wrapper function for LIMEX-calls to evaluate FUNC,BB in      *****
!*****                 BB(t,y) * y' = FUNC(t,y)                      *****
!*****                                                               *****
!*************************************************************************
      implicit none
      integer, intent(IN)  :: N
      integer, intent(OUT) :: Nz
      real*8 , intent(IN)  :: t
      real*8 , dimension(n), intent(IN)  :: y
      real*8 , dimension(n), intent(OUT) :: FUNC
      real*8 , dimension(*), intent(OUT) :: BB
      integer, dimension(*), intent(OUT) :: ir,ic
      integer, intent(INOUT)             :: iflag
      integer :: i
      Nz = N
      do i=1,N
         ir(i)=i     ! non-zero entries of BB in coordinate format
         ic(i)=i
         BB(i)=1.0   ! BB is and remains unity
      enddo

      call FF(N,t,y,FUNC,iflag)

      iflag = 0
      return
      end 
