
 
!=========================================================================
      subroutine JACOBI(N, t, Y, FF, JAC, ldjac, ml, mu, forb, jacinfo)
!=========================================================================
!*****                                                               *****  
!*****  wrapper function for LIMEX-calls to calculate the Jacobian   *****
!*****                                                               *****
!*************************************************************************
      implicit none
      integer N, ldjac, ml, mu, forb, jacinfo
      real*8 t, Y(n), FF(n), JAC(ldjac,*)
      write(*,*) "*** jacobi has been called!"
      stop
      return
      end 
