***********************************************************************
      REAL*8 FUNCTION DISPOL (T,arr)
***********************************************************************
*****                                                             *****
*****   Berechnet das Dissoziationspolynom                        *****
*****                                                             *****
*****   EINGABE:   T = 5040 / Gastemperatur[K]                    *****
*****            arr = Feld mit den Koeffizienten des Polynoms    *****
*****                                                             *****
***********************************************************************
      implicit none
      real*8   T,f,arr(0:4)
      f = arr(0) + T*(arr(1) + T*(arr(2) + T*(arr(3) + T*arr(4) ) ) )
      DISPOL = DEXP(f)
      end
