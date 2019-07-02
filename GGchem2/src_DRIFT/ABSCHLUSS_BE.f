**********************************************************************
      REAL*8 FUNCTION ABSCHLUSS_BED(abschluss,Vl,NMOM,LL)
**********************************************************************
      implicit none
      integer abschluss,NMOM
      real*8 LL(0:NMOM)
      real*8 Vl,L0,L1,L2,L3,L4,L5
      real*8 ppp,qqq,Lsol1,Lsol2
      real*8 alpha,beta,aaa,bbb,ccc,eee,ddd,del,delc
      real*8 y0,y1,y2,y3,y4,abest,qbest
      real*8 al
      real*8 VIETA1,VIETA2
      real*8 c0,c1,c2,c3,func,df,dL0,delta,fak,cmerk
      real*8 a1,a2,N1,N2
      integer it
      CHARACTER*1 char
      data L0/0.d0/, cmerk/1.d0/
      save L0,cmerk
*-----------------------------------------------------------------------
*  Die Formelfunktion zur Loesung quadratische Gleichungen mit Vieta
      VIETA1(ppp,qqq) = qqq/(-ppp/2.d0-DSQRT(ppp**2/4.d0-qqq))
      VIETA2(ppp,qqq) = qqq/(-ppp/2.d0+DSQRT(ppp**2/4.d0-qqq))
*-----------------------------------------------------------------------
      L1 = LL(1)
      L2 = LL(2)
      L3 = LL(3)
      L4 = LL(4)
c     L5 = LL(5)

      if (abschluss.eq.1) then
*       --------------------------------------------------------------
*       ***  Abschlussbedingung  2 L1^3 + L0(L0 L3 - 3 L1 L2) = 0  *** 
*       --------------------------------------------------------------
        if (L3.eq.0.d0) then
          if (L2.ne.0.d0) L0 = 2.d0/3.d0*L1**2/L2 
        else
          ppp = -3.d0*L1*L2/L3
          qqq = 2.d0*L1**3/L3
          if ((ppp.eq.0.d0).and.(qqq.eq.0.d0)) then
            L0 = 0.d0
          else if (ppp**2/4.d0-qqq.lt.0.d0) then 
            write(*,*) 'nur komplexe Loesungen der Abschlussbed.'
            write(*,*) L1,L2,L3,ppp,qqq
            stop
          else
            Lsol1 = VIETA1(ppp,qqq)
            Lsol2 = VIETA2(ppp,qqq)
            L0 = Lsol2
          endif
        endif
        if ((L0.lt.0.d0).or.(L0.ne.L0)) then
          write(*,*) 'keine positive reelle Loesung der Abschlussbed.'
          write(*,*) L1,L2,L3,ppp,qqq
          stop
        endif
      else if (abschluss.eq.2) then
        al  = Vl**(1.d0/3.d0)     ! nicht al, sondern nur Abkuerzung
*       -------------------------------------------------------
*       ***  L0 ( 4 al^3 L0^2 - 12 al^2 L0 L1 + 18 al L1^2  ***
*       ***     - 6 al L0 L2 + 9 L1 L2 - L0 L3 ) = 12 L1^3  ***
*       -------------------------------------------------------
        c3  =   4.d0*al**3
        c2  = -12.d0*al**2*L1 -6.d0*al*L2 - L3
        c1  =  18.d0*al*L1**2 +9.d0*L1*L2
        c0  = -12.d0*L1**3
        fak = 0.5
        it  = 0
 100    continue
          func  = c0 + c1*L0 +      c2*L0**2 +      c3*L0**3
          df    =      c1    + 2.d0*c2*L0    + 3.d0*c3*L0**2
          if (df.eq.0.d0) goto 110
          dL0   = -func/df
          delta = DABS(dL0/L0)
c         write(*,1000) L0,func,df,dL0,delta
          if (L0.eq.0.d0) then
            L0 = L0+dL0
          else
            L0 = DMIN1(DMAX1(L0+dL0,L0*fak),L0/fak)
          endif
          it = it + 1
          if (it.gt.200) then
            write(*,*) 'keine Konvergenz bei Abschlussbedingung'
            stop
          endif
        if (delta.gt.1.d-10) goto 100
 110    continue
c       write(*,*) it,'Iterationen',L0
      else if (abschluss.eq.3) then
*       --------------------------------------------------
*       ***  y(x) = alpha*x^beta  mit  y(x)=Lx/L(x-1)  ***
*       --------------------------------------------------
        L0 = 0.d0
        if ((L1.gt.0.d0).and.(L2.gt.0.d0).and.(L3.gt.0.d0)) then
          beta  = DLOG(L3*L1/L2**2)/DLOG(3.d0/2.d0)
          alpha = L2/L1/(2.d0)**beta
          L0    = L1/alpha
        endif
      else if (abschluss.eq.4) then
*       -----------------------------------------------
*       ***  y(x) = a*(x+c)^b  mit  y(x)=L(x+1)/Lx  ***
*       -----------------------------------------------
        L0 = 0.d0
        if ((L1.ne.0.d0).and.(L2.ne.0.d0).and.(L3.ne.0.d0)) then
          y1  = DLOG(L2/L1)
          y2  = DLOG(L3/L2)
          y3  = DLOG(L4/L3)
          ccc = 1.d0
          fak = 1.2d0
          it  = 0
 200      continue
            it = it + 1
            if ((it.gt.200)) then
              write(*,*) 'keine Konvergenz bei Abschlussbedingung'
              ccc = 1.d0
              goto 201
            endif
            func = (y1-y2)*DLOG(1.d0/4.d0+ccc) 
     &           + (y3-y1)*DLOG(1.d0/3.d0+ccc)
     &           + (y2-y3)*DLOG(1.d0/2.d0+ccc)
            df   = (y1-y2)/(1.d0/4.d0+ccc)
     &           + (y3-y1)/(1.d0/3.d0+ccc)
     &           + (y2-y3)/(1.d0/2.d0+ccc)
            if (df.eq.0.d0) then
              delc = ccc*1.01d0
            else
              delc = -func/df
            endif
            delta = DABS(delc/ccc)
            write(*,1000) ccc,func,df,delc,delta
            ccc = DMIN1(DMAX1(ccc+delc,ccc/fak),ccc*(fak+1.d-2))
          if (delta.gt.1.d-10) goto 200
 201      bbb = (y2-y1)/(DLOG(1.d0/3.d0+ccc)-DLOG(1.d0/2.d0+ccc))
          aaa = (L2/L1)/(1.d0/2.d0+ccc)**bbb
          alpha = aaa*(1.d0+ccc)**bbb
          if (alpha.ne.alpha) then
            write(*,*) 'geht nicht'
            ccc = 1.d0
            goto 201
          endif
          L0 = L1/alpha
          write(*,*) alpha,L2/L1,L3/L2,L4/L3
          write(*,*) DLOG(aaa)+bbb*DLOG(1.d0/2.d0+ccc),y1
          write(*,*) DLOG(aaa)+bbb*DLOG(1.d0/3.d0+ccc),y2
          write(*,*) DLOG(aaa)+bbb*DLOG(1.d0/4.d0+ccc),y3
          read(*,'(a1)') char
          cmerk = ccc
        endif
      else if (abschluss.eq.5) then
*       ----------------------------------------------------------
*       ***  y(x) = a+b/(1+exp(x-c))  mit  y(x)=ln(L(x+1)/Lx)  ***
*       ----------------------------------------------------------
        L0 = 0.d0
        if ((L1.ne.0.d0).and.(L2.ne.0.d0).and.(L3.ne.0.d0)) then
          y1  = DLOG(L2/L1)
          y2  = DLOG(L3/L2)
          y3  = DLOG(L4/L3)
          eee = DEXP(1.d0)
          aaa = (y1*(y3+eee*y3-y2)-eee*y2*y3) / (eee*y1-y2-eee*y2+y3)
          bbb = (1.d0-eee)*(y1-aaa)*(y2-aaa) / (y1-aaa-eee*(y2-aaa))
          ccc = 1.d0-DLOG((bbb-y1+aaa)/(y1-aaa))
          alpha = DEXP(aaa+bbb/(1.d0+DEXP(-ccc)))
          if ((alpha.gt.1.d-99).and.(alpha.lt.1.d+99)) then
c           write(*,*) aaa,bbb,ccc
c           write(*,*) aaa+bbb/(1.d0+DEXP(1.d0-ccc)),y1
c           write(*,*) aaa+bbb/(1.d0+DEXP(2.d0-ccc)),y2
c           write(*,*) aaa+bbb/(1.d0+DEXP(3.d0-ccc)),y3
c           read(*,'(a1)') char
          else
c           write(*,*) 'nehme Abschluss 3'
            beta  = DLOG(L3*L1/L2**2)/DLOG(3.d0/2.d0)
            alpha = L2/L1/(2.d0)**beta
          endif
c         write(*,*) alpha,L2/L1,L3/L2,L4/L3
          L0 = L1/alpha
        endif
      else if (abschluss.eq.6) then
*       ----------------------------------------------------------
*       ***  y(x) = b+c/(1+exp(2+ax))  mit  y(x)=ln(L(x+1)/Lx)  ***
*       ----------------------------------------------------------
        L0 = 0.d0
        if ((L1.ne.0.d0).and.(L2.ne.0.d0).and.(L3.ne.0.d0)) then
          y1 = DLOG(L2/L1)
          y3 = DLOG(L4/L3)
          y4 = DLOG(L5/L4)
          if ((cmerk.eq.0.d0).or.
     &        (DMAX1(y1,y3,y4).gt.DMIN1(y1,y3,y4)+0.01d0)) then
            cmerk = 0.d0
            aaa = 1.d0
            fak = 1.2d0
            it  = 0
            abest = aaa
            qbest = 1.d+99
 300        continue
              it = it + 1
              if ((it.gt.200)) then
                write(*,*) y1,y3,y4
                write(*,*) 'keine Konvergenz bei Abschlussbedingung 6'
                write(*,*) 'badness = ',qbest
                aaa = abest
                goto 301
              endif
              func = y4 - ( (1.d0+DEXP(aaa)+DEXP(2.d0*aaa))
     &             * (1.d0+DEXP(2.d0+3.d0*aaa)) * y3
     &             - DEXP(2.d0*aaa)*(1.d0+DEXP(2.d0+aaa)) * y1 )       
     &             / ((1.d0+DEXP(aaa))*(1.d0+DEXP(2.d0+4.d0*aaa)))
              df   = (y3-y1)*DEXP(2.d0*aaa) 
     &             * (-2.d0-2.d0*DEXP(2.d0+2.d0*aaa)
     &                +2.d0*DEXP(2.d0+4.d0*aaa)+2.d0*DEXP(4.d0+6.d0*aaa)
     &                +DEXP(2.d0+5.d0*aaa)*(3.d0+DEXP(2.d0))
     &                -DEXP(aaa)*(1.d0+3.d0*DEXP(2.d0)))
     &             / ((1.d0+DEXP(aaa))*(1.d0+DEXP(2.d0+4.d0*aaa)))**2 
              if (DABS(func).lt.qbest) then
                qbest = DABS(func)
                abest = aaa
              endif
              ddd   = -func/df
              delta = DABS(ddd/aaa)
              write(*,1000) aaa,func,df,ddd,delta
              aaa = DMAX1(DMIN1(aaa+ddd,aaa*(1.d-2+fak)),aaa/fak)
            if (delta.gt.1.d-10) goto 300
 301        bbb = (y1 + DEXP(-2.d0-aaa)*(y1-y3) - DEXP(2.d0*aaa)*y3)
     &           /(1.d0-DEXP(2.d0*aaa))
            ccc = (y1-y3)/(1.d0/(1.d0+DEXP(2.d0+aaa))
     &                    -1.d0/(1.d0+DEXP(2.d0+3.d0*aaa)))
            y0 = bbb+ccc/(1.d0+DEXP(2.d0))
            y2 = bbb+ccc/(1.d0+DEXP(2.d0+aaa*2.d0))
            L0 = L1/DEXP(y0)
            write(*,1000) L1/L0,L2/L1,L3/L2,L4/L3,L5/L4
            write(*,*) aaa,bbb,ccc
            write(*,*) 'Y1:',y1,bbb+ccc/(1.d0+DEXP(2.d0+aaa*1.d0))
            write(*,*) 'Y3:',y3,bbb+ccc/(1.d0+DEXP(2.d0+aaa*3.d0))
            write(*,*) 'Y4:',y4,bbb+ccc/(1.d0+DEXP(2.d0+aaa*4.d0))
            write(*,*) 'Quatlitaet= ',DABS(1.d0-DLOG(L3/L2)/y2)
            read(*,'(a1)') char
          else
c           write(*,*) 'nehme Abschluss 3'
            beta  = DLOG(L3*L1/L2**2)/DLOG(3.d0/2.d0)
            alpha = L2/L1/(2.d0)**beta
            L0 = L1/alpha
          endif
        endif
      else if (abschluss.eq.7) then
*       ----------------------------------------------------------------------
*       ***  y(x) = ((a+bx)+(x/d)^2*(a+cx))/(1+(x/d)^2)  mit  y(x)=ln(Lx)  ***
*       ----------------------------------------------------------------------
        L0 = 0.d0
        if ((L1.ne.0.d0).and.(L2.ne.0.d0).and.(L3.ne.0.d0)) then
          y1 = DLOG(L1)
          y2 = DLOG(L2)
          y3 = DLOG(L3)
          y4 = DLOG(L4)
          if ((cmerk.eq.0.d0).or.(DABS(y4-y3-y2+y1).gt.0.01d0)) then
            cmerk = 0.d0
            ddd = 2.d0
            fak = 1.2d0
            it  = 0
 400        continue
              it = it + 1
              if ((it.gt.200)) then
                write(*,*) y1,y2,y3,y4
                write(*,*) 'keine Konvergenz bei Abschlussbedingung 7'
                stop
              endif
              func = (-24.d0*y1+144.d0*y2-216.d0*y3 +96.d0*y4)
     &             + ( 26.d0*y1-228.d0*y2+378.d0*y3-176.d0*y4)*ddd
     &             + (-15.d0*y1 -60.d0*y2+165.d0*y3 -90.d0*y4)*ddd**2
     &             + ( 25.d0*y1 -45.d0*y2 +15.d0*y3  +5.d0*y4)*ddd**3
     &             + (  9.d0*y1 -24.d0*y2 +21.d0*y3  -6.d0*y4)*ddd**4
     &             + ( -1.d0*y1  +3.d0*y2  -3.d0*y3  +1.d0*y4)*ddd**5
              df   = ( 26.d0*y1-228.d0*y2+378.d0*y3-176.d0*y4)
     &        + (-15.d0*y1 -60.d0*y2+165.d0*y3 -90.d0*y4)*2.d0*ddd
     &        + ( 25.d0*y1 -45.d0*y2 +15.d0*y3  +5.d0*y4)*3.d0*ddd**2
     &        + (  9.d0*y1 -24.d0*y2 +21.d0*y3  -6.d0*y4)*4.d0*ddd**3
     &        + ( -1.d0*y1  +3.d0*y2  -3.d0*y3  +1.d0*y4)*5.d0*ddd**4
              del   = -func/df
              delta = DABS(del/ddd)
c             write(*,1000) ddd,func,df,del,delta
              ddd = DMAX1(DMIN1(ddd+del,ddd*(1.d-2+fak)),ddd/fak)
            if (delta.gt.1.d-10) goto 400
            ccc = -((-1.d0+ddd)*(6.d0+ddd)*(1.d0+ddd**2)*y1
     &          -2.d0*(4.d0+ddd**2)*(-3.d0+ddd*(4.d0+ddd))*y2
     &          +(9.d0+ddd**2)*(-2.d0+ddd*(3.d0+ddd))*y3)
     &          /(2.d0*(6.d0+ddd*(-11.d0+(-6.d0+ddd)*ddd)))
            bbb = (ccc*(-4.d0+ddd**2*(-7.d0+3.d0*ddd))
     &          -(1.d0+ddd**2)*(4.d0+ddd**2)*(y1-y2))
     &          /(ddd**2*(-2.d0+ddd*(3.d0+ddd)))
            aaa = (-(-2.d0+ddd)*(3.d0*ccc-3.d0*ccc*ddd+y1+ddd**2*y1)
     &          +(-1.d0+ddd)*(4.d0+ddd**2)*y2)
     &          /(-2.d0+ddd*(3.d0+ddd))
            L0 = DEXP(aaa-bbb*ddd)
c           write(*,*) aaa,bbb,ccc,ddd
c           write(*,*) 'Y1:',y1,((aaa+bbb*(1.d0-ddd))+(1.d0/ddd)**2
c    &                      *(aaa+ccc*(1.d0-ddd)))/(1.d0+(1.d0/ddd)**2)
c           write(*,*) 'Y2:',y2,((aaa+bbb*(2.d0-ddd))+(2.d0/ddd)**2
c    &                      *(aaa+ccc*(2.d0-ddd)))/(1.d0+(2.d0/ddd)**2)
c           write(*,*) 'Y3:',y3,((aaa+bbb*(3.d0-ddd))+(3.d0/ddd)**2
c    &                      *(aaa+ccc*(3.d0-ddd)))/(1.d0+(3.d0/ddd)**2)
c           write(*,*) 'Y4:',y4,((aaa+bbb*(4.d0-ddd))+(4.d0/ddd)**2
c    &                      *(aaa+ccc*(4.d0-ddd)))/(1.d0+(4.d0/ddd)**2)
c           read(*,'(a1)') char
          else
c           write(*,*) 'nehme Abschluss 3'
            beta  = DLOG(L3*L1/L2**2)/DLOG(3.d0/2.d0)
            alpha = L2/L1/(2.d0)**beta
            L0 = L1/alpha
          endif
        endif
      else if (abschluss.eq.8) then
*       -------------------------------------------------
*       ***  f(a) = N1*Dirac(a-a1) + N2*Dirac(a-a2)   ***
*       -------------------------------------------------
        L0 = 0.d0
        if ((L1.gt.0.d0).and.(L2.gt.0.d0).and.(L3.gt.0.d0)) then
          call DIST_FROM_MOMENTS1(1.d0,L1,L2,L3,L4,a1,a2,N1,N2,L0)  
        endif
      else
        write(*,*) 'Abschlussbedingung = ',abschluss,' => hae?'
        stop
      endif

      ABSCHLUSS_BED = L0
      RETURN 
 1000 format(99(1pE12.5))
      end
