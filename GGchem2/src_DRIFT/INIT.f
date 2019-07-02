**********************************************************************
      SUBROUTINE INIT(eps0)
**********************************************************************
*****                                                            *****
*****   Initialisiert Elementhaeufigkeiten                       *****
*****   - Lambert + Rao (1994):                                  *****
*****     (JAA 15, 47, solare Werte); He/H=0.1; Ne,Mg,K,Cr,Ti    *****
*****     und Mn nach nach Allen (1973)                          *****
*****   - Anders + Grevesse (1989):                              *****
*****     Geochimica et Cosmochemica Acta Vol 53, pp 197--214    *****
*****     ("Photosphere")                                        *****
*****   - wie in MARCS-Code                                      *****
*****   - wie in Tsuji-Chemie                                    *****
*****                                                            *****
**********************************************************************
      use drift_data,ONLY: NELEM,eps,mass,muH,elnam,amu
      implicit  none
      real*8,intent(IN) :: eps0(NELEM)
      integer i,H,He,C,N,O,Ne,Na,Mg,Al,Si,S,K,Ca,Cr,Mn,Fe,Ni,Ti
*     -----------------------------------------------------------------
      data H/1/, He/2/, C/6/, N/7/, O/8/, Ne/10/, Na/11/, Mg/12/,Al/13/
      data Si/14/, S/16/, K/19/, Ca/20/, Cr/24/, Mn/25/, Fe/26/
      data Ni/28/, Ti/22/
*     -----------------------------------------------------------------
      write(*,*) 
      write(*,*) "elemental abundances ..."
      write(*,*) "========================"
      do i=1,NELEM
         elnam(i) = '  '
         eps(i)  = -99.d0
         mass(i) = 0.d0
      enddo
      elnam(1)  = 'H'
      elnam(6)  = 'C '
      elnam(7)  = 'N '
      elnam(8)  = 'O '
      elnam(11) = 'Na'
      elnam(12) = 'Mg'
      elnam(13) = 'Al'
      elnam(14) = 'Si'
      elnam(16) = 'S '
      elnam(19) = 'K '
      elnam(20) = 'Ca'
      elnam(22) = 'Ti'
      elnam(26) = 'Fe'
       
*     -------------------------
*     Anders + Grevesse (1989):
*     -------------------------
      eps(H)  = 12.00 D0
      eps(He) = 10.99 D0
      eps(C)  =  8.56 D0
      eps(N)  =  8.05 D0
      eps(O)  =  8.93 D0
      eps(Ne) =  8.09 D0
      eps(Na) =  6.33 D0
      eps(Mg) =  7.58 D0
      eps(Al) =  6.47 D0
      eps(Si) =  7.55 D0
      eps(S)  =  7.21 D0
      eps(K)  =  5.12 D0
      eps(Ca) =  6.36 D0
      eps(Cr) =  5.67 D0
      eps(Mn) =  5.39 D0
      eps(Fe) =  7.67 D0
      eps(Ni) =  6.25 D0
      eps(Ti) =  4.99 D0

c      ------
c      MARCS:
c      ------
c      eps(H)  = 12.00 D0
c      eps(He) = 10.99 D0
c      eps(C)  =  8.55 D0
c      eps(N)  =  8.05 D0
c      eps(O)  =  8.93 D0
c      eps(Ne) =  8.09 D0
c      eps(Na) =  6.33 D0
c      eps(Mg) =  7.58 D0
c      eps(Al) =  6.47 D0
c      eps(Si) =  7.55 D0
c      eps(S)  =  7.21 D0
c      eps(K)  =  5.12 D0
c      eps(Ca) =  6.36 D0
c      eps(Cr) =  5.67 D0
c      eps(Mn) =  5.40 D0
c      eps(Fe) =  7.51 D0
c      eps(Ni) =  6.25 D0
c      eps(Ti) =  4.99 D0

c      ------
c      Tsuji:
c      ------
c      eps(H)  = 12.00 D0
c      eps(He) = 11.00 D0
c      eps(C)  =  8.56 D0
c      eps(N)  =  8.05 D0
c      eps(O)  =  8.93 D0
c      eps(Ne) =  8.09 D0
c      eps(Na) =  6.33 D0
c      eps(Mg) =  7.58 D0
c      eps(Al) =  6.47 D0
c      eps(Si) =  7.55 D0
c      eps(S)  =  7.21 D0
c      eps(K)  =  5.12 D0
c      eps(Ca) =  6.36 D0
c      eps(Cr) =  5.67 D0
c      eps(Mn) =  5.40 D0
c      eps(Fe) =  7.51 D0
c      eps(Ni) =  6.25 D0
c      eps(Ti) =  4.99 D0

*     --------------------
*     ***  Atommassen  ***
*     --------------------
      mass(H)  = 1.008  D0 * amu
      mass(He) = 4.0026 D0 * amu
      mass(C)  = 12.011 D0 * amu
      mass(N)  = 14.007 D0 * amu
      mass(O)  = 15.999 D0 * amu
      mass(Ne) = 20.18  D0 * amu            
      mass(Na) = 22.990 D0 * amu
      mass(Mg) = 24.312 D0 * amu
      mass(Al) = 26.98  D0 * amu
      mass(Si) = 28.086 D0 * amu
      mass(S)  = 32.064 D0 * amu
      mass(K)  = 39.10  D0 * amu  
      mass(Ca) = 40.08  D0 * amu
      mass(Cr) = 52.00  D0 * amu
      mass(Mn) = 54.94  D0 * amu
      mass(Fe) = 55.847 D0 * amu
      mass(Ni) = 58.71  D0 * amu
      mass(Ti) = 47.90  D0 * amu

      muH = 0.d0
      do i=1,NELEM
        eps(i) = 10.d0 ** (eps(i)-12.d0)
        if (elnam(i).ne.'  ') then
          write(*,1000) elnam(i),eps(i),eps0(i),
     &                  DABS(eps(i)-eps0(i))/DMIN1(eps(i),eps0(i))
        endif  
        eps(i) = eps0(i)
        muH = muH + mass(i)*eps(i)
      enddo
      write(*,*) 'rho = n<H> *',muH/amu
*
      RETURN
 1000 format(1x,a2,1x,2(1pE13.6),0pF7.4)
      end
