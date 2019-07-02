*********************************************************************
      SUBROUTINE NUCLEATION(T,Jst,Nstern)
*********************************************************************
      use drift_data,ONLY: NNUC,NMOLE,nmol,cmol,J_is_zero,
     &                     nuc_nam,nat
      implicit none
      real*8,intent(IN) :: T
      real*8,intent(OUT):: Nstern,Jst(NNUC)
      real*8 :: nTiO2,SS,Jstern
      real*8 :: nC1,nC2,nC2H,nC2H2,nC3
      INTEGER :: stindex,i
      INTEGER,save :: TiO2,C,C2,C3,C2H,C2H2
      logical,save :: firstCall=.true.
*
      if (firstCall) then
        TiO2 = STINDEX(CMOL,NMOLE,'TIO2     ')
        C    = 6
        C2   = STINDEX(CMOL,NMOLE,'C2       ') 
        C2H  = STINDEX(CMOL,NMOLE,'C2H      ')
        C3   = STINDEX(CMOL,NMOLE,'C3       ')
        C2H2 = STINDEX(CMOL,NMOLE,'C2H2     ')
        firstCall=.false.
      endif    
*
      do i=1,NNUC
        if      (nuc_nam(i).eq.'TiO2   ') then
        nTiO2 = nmol(TiO2)
        CALL KLASS_NUK_TIO2(T,nTiO2,Jstern,Nstern,SS)
        if (J_is_zero) Jstern=0.d0
        Jst(i) = Jstern
      else if (nuc_nam(i).eq.'Carb') then
       nC1  = nat(C)
       nC2  = nmol(C2)
       nC3  = nmol(C3)
       nC2H = nmol(C2H)
c       print*, nC1,nC2,nC3,nC2H
       CALL KLASS_NUK_C(T,nC1,nC2,nC2H,nC2H2,nC3,Jstern,Nstern,SS)
c       print*, Jstern, Nstern, SS
       if (J_is_zero) Jstern=0.d0
        Jst(i) = Jstern
      endif
      enddo

      RETURN 
      END 
