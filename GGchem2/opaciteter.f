I begyndelsen af hovedprogrammet:


      call oldsta
      call mainb

      io=-1
      call initjn(io)
      call modjon(jontyp,io)
      call initab(io)
      





      if(idust .eq. 1) then
        if(newmod .ne. 0) stop 'Error: Dust only works for NEWMOD=0'
        if(jontyp .ne. 3) stop 'Error: Dust only works for JONTYP=3'
        if(jump .ne. 2 .and. 4) stop 'Error: Dust only works for JUMP=2'
        if(isph .eq. 1) stop 'Error: Dust only works for ISPH=0'
        
        print *, 'Reading DRIFT file...'
        call drift2marcs
        print *, 'Done.'
        print *
        
        call dust_eps
        call dust_opac
      end if

      molh = 0                           ! molh=1 => only h,h2,h2+ in molec. eq.      

      metals = 0
      do im = 1,nosmol
        write(6,'(a18,i3,a45)') ' im,molname(im) = ', im, molname(im)
        if(molname(im).eq. 'ATOM') metals = im
        if(molname(im).eq.'CH4 ' .and. jump.eq.0) then 
          stop ' P(CH4) not comp. by old Marcs eq.; set JUMP>0'
        end if
      end do
      write(6,*) ' IM for metals = ', metals
      if(metals.ge.1) call osinit        ! atomic lines included in computation

      read(5,outlist)
      pfd=itmax.lt.0
      if(pfd) itmax=-itmax
      pfe=itmax.le.nprint      

      if(newmod.eq.1) call startm
      if(newmod.eq.2) print*,' newmo=2; no iteration; trans old mod'
      if(newmod.eq.3) call scale(22)
      if(newmod.eq.4) then
        call resume(22,1)
        if (isph.eq.1) then 
          call tryck_sph
        else 
          call tryck
        end if
      end if
      if(newmod.eq.5) call coscal(22)
      if(newmod.eq.6) call oldarc(22)
      if(newmod.eq.7) read(22)
      if(newmod.eq.8) then
        read(22)
        call scale(22)
      end if
      if(newmod.eq.9) open(66,file='model.lim',form='unformatted')
      if(newmod.ne.1) then
        if(isph.eq.1) then 
          call presnt_sph
        else 
          call presnt
        end if
      end if

      print *, 'Initiation finished.'
      
! Iterate model
      lun=22
      itstop=.false.
      quit=.false.
      onemor=.false.

      do it=1,itmax 
        write(6,*) ' **** NOW ITERATION ',it,' IS STARTING ****'
        write(66,*)' **** NOW ITERATION ',it,' IS STARTING ****'
        write(111,*) 'Iteration' , it
        print *
        write(*,'(a23,i2,a8)')' Iteration #           ',it,' started'
        call gettime(1)
        
        pf=it.gt.itmax-nprint
        if(quit) pf=.true.

        if (isph.eq.1) then 
          call solve_sph(1)
        else 
          call solve(1)
        end if
        write(6,*) ' **** NOW AFTER CALL TO SOLVE **** '
        write(6,*) 'pe(1)= ',ppe(1)
        
        if(newmod.eq.2) go to 102
        if(newmod.eq.9) goto 101
        
        if(isph.eq.1) then 
          call matrix_sph
        else 
          call matrix
        end if
        
        write(6,*) 'pe(1)= ',ppe(1)
        write(6,*) ' **** NOW AFTER CALL TO MATRIX ***** '
      
        call newsta
        
        if(itstop .and. onemor) quit=.true.
        if(itstop) onemor=.true.
        if(quit) exit
      end do


! Save model
!      call savemodel
      call marcs2drift
      
      if(noarch.ge.2) go to 101
      call modjon(3,1)
      if(noarch.eq.1) goto 101
      call newsta
102   continue
      call archiv(22)
      write(6,*) ' **** NOW calling listmo  **** '
      call listmo(1,22,isph)

! End program
101   continue
      write(7,*) ' *** NORMAL END ***'
      call gettime(1)
      print *
      print *, 'Successful termination of execution!'
      print *

      stop
      end
      
C
      SUBROUTINE ABSKO(NEWT,NT,TSKAL,PESKAL,ISETA,J,ABSK,SPRID)
      implicit real*8 (a-h,o-z)
C
C        THE ROUTINE ADMINISTERS THE COMPUTATION OF ABSORPTION
C        COEFFICIENTS. IT CALLS THE ROUTINES, GIVING THE PROPER THERMO-
C        DYNAMIC INFORMATION ( J O N ) , THE DETAILS OF THE ABSORPTION
C        MECHANISMS ( D E T A B S ) AND THE FACTORS FOR THE INTERPOLATION
C        IN T  ( T A B S ) . IT CHOOSES (IF NECESSARY READS) THE RIGHT SET
C        OF ABSORPTION-COEFFICIENT DATA (ISETA), STATEMENT NO. 5 AND MAKES
C        THE INTERPOLATION IN T, STATEMENTS NOS. 10-18, AND THE SUMMATION
C        OF A ROSSELAND MEAN, IF INDICATED BY J = 0, STATEMENTS NOS. 25-28
C
C        NEWT SHOULD BE GT 1 THE FIRST TIME THIS ROUTINE IS USED,
C                       EQ 1 WHEN A NEW SET OF T-PE IS USED,
C                       EQ 0 OTHERWISE.
C
C        NT IS THE NUMBER OF T-PE POINTS. THE TEMPERATURES T SHOULD BE EX-
C        PRESSED IN KELVIN, THE ELECTRON PRESSURES PE IN DYNES PER CM2
C        ISETA IS THE WAVELENGTH-SET NUMBER, J THE WAVELENGTH NUMBER IN THAT
C        SET. J EQUAL TO ZERO INDICATES THAT A ROSSELAND MEAN IS WANTED.
C        THIS MEAN IS COMPUTED USING THE WAVELENGTH POINTS OF THE ACTUAL
C        SET (ISETA) AND THE QUADRATURE WEIGHTS GIVEN IN ROSW.
C        IN ABSK AND SPRID THE ABSORPTION AND SCATTERING COEFFICIENTS PER GRAM
C        OF STELLAR MATTER ARE STORED.
C
C        DIMENSIONS NECESSARY
C        AB(NKOMP),ABSK(1),FAKT(NKOMP),FAKTP(IFADIM),NTPO(NTO),PE(NT),PESKAL(1),
C        ROSW(MAX(NL)),SPRID(1),SUMW(NT),T(NT),TSKAL(1),XLA(MAX(NL)),
C        XLA3(MAX(NL))
C        THE DIMENSIONS ARE LOWER LIMITS.
C        DIMENSIONS OF ARRAYS IN COMMONS /CA2/,/CA3/ AND /CFIL/ ARE COMMENTED
C        ON IN SUBROUTINE INABS, THOSE OF ARRAYS IN COMMON /CA4/ IN SUBROUTINE
C        TABS.
C        NKOMP IS THE NUMBER OF 'COMPONENTS'
C        NL(I) IS THE NUMBER OF WAVELENGTHS IN WAVELENGTH SET I
C        NT IS THE NUMBER OF T-PE POINTS SUPPLIED IN TSKAL AND PESKAL
C        NTO IS THE NUMBER OF POINTS IN THOSE SCALES FOR WHICH A DETAILED
C              PRINT-OUT IS WANTED.
C
      include 'parameter.inc'
C
C      PARAMETER (KFADIM=4000,IFADIM=1000)
      DIMENSION TSKAL(NDP),PESKAL(NDP),ABSK(NDP),SPRID(NDP)
      DIMENSION FAKTP(ifadim)
      DIMENSION SUMW(NDP)
      COMMON/UTPUT/IREAD,IWRIT
      COMMON/CA2/ABKOF(4000),KOMPLA(600),KOMPR,KOMPS,NKOMP
      COMMON/CA3/ILOGTA(30),NULL
      COMMON/CA4/AFAK(KFADIM),NOFAK(IFADIM),NPLATS(IFADIM)
      COMMON/CA5/AB(30),FAKT(30),PE(NDP),T(NDP),XLA(20),XLA3(20),RO,
     &           SUMABS,SUMSCA,VIKTR,ISET,NLB
      COMMON/CFIL/IRESET(10),ISLASK,IREAT
      COMMON/COUTR/NTO,NTPO(10)
      COMMON/CROS/ROSW(20)
      COMMON /CARC3/ F1P,F3P,F4P,F5P,HNIC,PRESMO(33)
      COMMON /CARC4/ PROV(30),NPROVA,NPROVS,NPROV
      COMMON /TIO/PTIO(NDP),ROsav(NDP),POXG1(NDP)
      COMMON/CI4/ TMOLIM,IELEM(16),ION(16,5),MOLH,JUMP
      COMMON/CMOL1/EH,FE,FH,FHE,FC,FCE,FN,FNE,FO,FOE,FK,FKE,FS,FSE
     &             ,FT,FTE
      COMMON /DENSTY/ ROTEST(NDP),PRH2O(NDP)
      common /fullequilibrium/ partryck(ndp,maxmol),
     &  xmettryck(ndp,maxmet),xiontryck(ndp,maxmet)
      COMMON /CMOLRAT/ FOLD(NDP,8),MOLOLD,KL

      INTEGER MOLH, JUMP
C
C
      ISET=ISETA
      IF(NEWT.GT.1)ISETP=-1
      IF(NEWT.EQ.0)GO TO 5
C
C        FACTORS ONLY DEPENDENT ON T-PE
C
C      write(6,1311) nt,tskal
C1311  format(' before tabs: nt,tskal:',i3,7f8.1,5(/9f8.1))
      CALL TABS(NT,TSKAL)
C      write(6,1312) nt,tskal
C1312  format(' after  tabs: nt,tskal:',i3,7f8.1,5(/9f8.1))
      IFAK=1
      KFAK=1
      JP=0
      KP=1
C
C        LOOP OVER THE T-PE POINTS ('THE FIRST NTP-LOOP')
      DO4 NTP=1,NT
      T(NTP)=TSKAL(NTP)
      PE(NTP)=PESKAL(NTP)
C        IS PRINT-OUT WANTED FOR T-PE POINT NO. NTP
      IOUTR=0
      IF(KP.GT.NTO)GO TO 3
      IF(NTP.EQ.NTPO(KP))GO TO 1
      GO TO 3
    1 IOUTR=1
      KP=KP+1
    3 CONTINUE
C
      if(j.le.0) then
        molhs=molh
        molh =0
      endif
C      write(6,*) ' calling jon, j,kl,ntp = ',j,kl,ntp
C       write(6,1313) ntp,t(ntp),pe(ntp)
C1313  format( ' before jon: k,t(k),pe(k): ',i3,f8.1,1p3e11.3)

      if(pe(ntp).le.1.e-33) then
       write(6,*) ' ***** now we are in trouble:'                   
       write(6,*) ' maxmol,maxmet = ',maxmol,maxmet
       write(6,*) ' ntp,nt,j = ',ntp,nt,j
       write(6,1320) tskal
       write(6,1321) peskal
1320  format(' tskal:',7f8.1,5(/9f8.1))
1321  format(' peskal:',1p7e10.2,6(/8e10.2))
       write(6,1322) t
       write(6,1323) pe
1322  format(' t:',7f8.1,5(/9f8.1))
1323  format(' pe:',1p7e10.2,6(/8e10.2))
      end if

      CALL JON(T(NTP),PE(NTP),1,PG,RO,DUM,IOUTR)

C1314  format( ' after jon: k,t(k),pe(k),pg,ro: ',i3,f8.1,1p3e11.3)
C      write(6,1314) ntp,t(ntp),pe(ntp),pg,ro

      if(j.le.0) then
        molh=molhs
        rosav(ntp)=ro
        poxg1(ntp)=pe(ntp)*foe

**********18.12.94 

        IF (JUMP.GE.1) THEN
         prh2o(ntp)=partryck(ntp,4)
         ELSE
         prh2o(ntp)=presmo(4)
        ENDIF
        
         if (ntp.le.0) then
         write(6,*) 'ntp (must be NE 0 for dimension)',ntp
         write(6,*) 'prh2o(ntp)      : ', prh2o(ntp)
         write(6,*) 'partryck(ntp,4) : ',partryck(ntp,4)
         write(6,*) 'should be equal'
         end if
*   well, it is / 11.1.95

*
* it works but it starts with ntp = 0 --> is that right ???
* yes it is, because ntp=0 means optical depth=0 and that means you start on
* the surface of the star
*
***********
      endif

      CALL DETABS(J,0,NTP,IOUTR)
C
C        WE STORE THE FAKT ARRAY, MADE IN JON-DETABS IN LONGER ARRAYS NAMELY
C                  IN AFAK FOR TEMPERATURE-INDEPENDENT COMPONENTS
C                  IN FAKTP FOR TEMPERATURE-DEPENDENT ONES
      DO2 KOMP=1,KOMPR
      AFAK(KFAK)=FAKT(KOMP)
    2 KFAK=KFAK+1
      DO4 KOMP=KOMPS,NKOMP
      FAKTP(IFAK)=FAKT(KOMP)
      KFAK=KFAK+NOFAK(IFAK)
    4 IFAK=IFAK+1
C        END OF 'THE FIRST NTP-LOOP'
C
C        READING  OF A NEW WAVELENGTH SET IF INDICATED BY ISET
    5 IF(ISET.EQ.ISETP)GO TO 6
      IREADP=IRESET(ISET)
   51 READ(IREADP,END=52)ISETP,NLB,XLA,XLA3,NABKOF,ABKOF,NKOMPL,KOMPLA
      GO TO 5
   52 REWIND IREADP
      GO TO 51
C        ROSSELAND MEAN OR NOT
    6 IF(J.GT.0)GO TO 9
    7 J1=1
      J2=NLB
      DO8 NTP=1,NT
      SUMW(NTP)=0.
    8 ABSK(NTP)=0.
      GO TO 10
    9 J1=J
      J2=J
C
C        INTERPOLATION IN T
C        LOOP OVER ALL THE WAVELENGTHS IN CASE OF ROSSELAND MEAN. THIS
C        LOOP ENDS IN STATEMENT NO. 26
   10 CONTINUE
      DO26 JP=J1,J2
      KFAK=1
      IFAK=1
      KP=1
C
C        LOOP OVER THE T-PE POINTS ('THE SECOND NTP-LOOP')
      DO26 NTP=1,NT
C
C        IS PRINT-OUT WANTED FOR T-PE POINT NO. NTP
      IOUTR=0
      IF(KP.GT.NTO)GO TO 93
      IF(NTP.EQ.NTPO(KP))GO TO 92
      GO TO 93
   92 IOUTR=1
      KP=KP+1
      IF(KP.EQ.2)IOUTR=2
   93 CONTINUE
      IU=JP
C
C        COMPONENTS WITH ABSORPTION COEFFICIENTS INDEPENDENT OF THE
C        TEMPERATURE
C
      DO14 KOMP=1,KOMPR
      IF(KOMPLA(IU).LE.0)GO TO 12
C        THE VECTOR KOMPLA IS DETERMINED IN SUBROUTINE INABS.
C        KOMPLA GREATER THAN ZERO GIVES THE INDEX IN ABKOF, WHERE THE TABLE FOR
C        THIS COMPONENT AND WAVELENGTH BEGINS.
C        KOMPLA LESS THAN OR EQUAL TO ZERO INDICATES THAT THE ACTUAL ABSORPTION
C        COEFFICIENT FOR THIS COMPONENT AND WAVELENGTH IS ZERO, AS FOUND IN SUB-
C        ROUTINE INABS.
C
   11 INDEX=KOMPLA(IU)
      AB(KOMP)=AFAK(KFAK)*ABKOF(INDEX)
      GO TO 13
   12 AB(KOMP)=0.
   13 KFAK=KFAK+1
   14 IU=IU+NLB
C
C        COMPONENTS WITH T-DEPENDENT ABSORPTION COEFFICIENTS
      DO19 KOMP=KOMPS,NKOMP
      NOP=NOFAK(IFAK)
      IF(NOP.EQ.0)GO TO 17
      IF(KOMPLA(IU).LE.0)GO TO 17
   15 INDEX=NPLATS(IFAK)-1+KOMPLA(IU)
C        THE VECTOR NPLATS IS DETERMINED BY SUBROUTINE TABS. IT GIVES THE ARRAY
C        INDEX OF THE TEMPERATURE AT WHICH THE INTERPOLATION IN ABKOF
C        BEGINS. NOFAK, GIVING INFORMATION ON THE T-INTERPOLATION AND
C        POSSIBLY INDICATING THAT AB=0 (NOFAK=0) IS ALSO DETERMINED BY TABS
C
C        INTERPOLATION
      DELSUM=0.
      DO16 NP=1,NOP
      DELSUM=DELSUM+AFAK(KFAK)*ABKOF(INDEX)
      KFAK=KFAK+1
   16 INDEX=INDEX+1
C
C        HAS THE INTERPOLATION BEEN MADE ON THE LOGARITHM
      IF(ILOGTA(KOMP).GT.0)DELSUM=EXP(DELSUM)
C        MULTIPLICATION BY FACTOR FROM JON-DETABS
      DELSUM=DELSUM*FAKTP(IFAK)
      IF(DELSUM.GE.0)GO TO 162
C
C        A NEGATIVE INTERPOLATION RESULT
  161 IF(NULL.GT.0)WRITE(IWRIT,200)KOMP,DELSUM,JP,ISET,T(NTP)
  200 FORMAT(4H AB(,I4,11H) NEGATIVE=,E12.4,5X,17HFOR WAVELENGTH NO,I5,
     *5X,6HSET NO,I5,5X,2HT=,F10.4,'  AND THEREFORE PUT =0 ***ABSKO***')
      AB(KOMP)=0.
      GO TO 18
  162 AB(KOMP)=DELSUM
      GO TO 18
   17 AB(KOMP)=0.
      KFAK=KFAK+NOP
   18 IU=IU+NLB
   19 IFAK=IFAK+1
C
C        WE MULTIPLY BY WAVELENGTH-DEPENDENT  FACTORS AND ADD UP. THIS IS
C        DONE IN DETABS.
      CALL DETABS(J,JP,NTP,IOUTR)
C
      IF(J.LE.0)GO TO 25
   24 ABSK(NTP)=SUMABS
      SPRID(NTP)=SUMSCA
      GO TO 26
C
C        SUMMATION TO GET A ROSSELAND MEAN
   25 CONTINUE
C use only central wavelengthinterval for Rosseland to avoid too
C large "abitrariness" on the value due to very small kap_nu in the
C wings of the planck function.
      if (xla(jp).le.5000. .or. xla(jp).ge.1.e5) go to 26
      IF(J.EQ.0) ABSK(NTP)=ABSK(NTP)+ROSW(JP)*VIKTR/(SUMABS+SUMSCA)
      IF(J.LT.0) ABSK(NTP)=ABSK(NTP)+ROSW(JP)*VIKTR/SUMABS
      SUMW(NTP)=SUMW(NTP)+ROSW(JP)*VIKTR
   26 CONTINUE
C
C        END OF 'THE SECOND NTP-LOOP'
C
      IF(J.GT.0)GO TO 29
      DO28 NTP=1,NT
      SPRID(NTP)=0.
   28 ABSK(NTP)=SUMW(NTP)/ABSK(NTP)
C
   29 CONTINUE
      RETURN
      END

C*
C*NEW PDS MEMBER FOLLOWS
C*
      SUBROUTINE ACCEL(CORR,N)
      implicit real*8 (a-h,o-z)
C
      DIMENSION CORR(40), CORRO(40), KONV(40), CORRUT(40)
      DATA NCALL/0/
C
      IF(NCALL.EQ.0) THEN
        DO 10 I=1,N
          KONV(I) = 0
          CORRO(I) = CORR(I)
   10   CONTINUE
        NCALL = NCALL + 1
        RETURN
      ENDIF
C
      MAXCOR=0
      DO 40 I=1,N
C
        IF(CORR(I)*CORRO(I).GT.0.) THEN
          KONV(I)=KONV(I) + 1
        ELSE
          KONV(I)=0
        ENDIF
C
        MAXCOR=MAX0(MAXCOR,KONV(I))
        CORRO(I)   = CORR(I)
   40 CONTINUE
C
      IF(MAXCOR.LT.2) RETURN
C
      IF(KONV(1).EQ.5) THEN
        CORRUT(1) = 5.*CORR(1)
      ELSE IF(KONV(2).EQ.5) THEN
        CORRUT(1) = 2.*CORR(1)
      ENDIF
      DO 50 I=2,N-1
        IF(KONV(I).EQ.5) THEN
          CORRUT(I) = 5.*CORR(I)
        ELSE IF(KONV(I+1).EQ.5 .OR. KONV(I-1).EQ.5) THEN
          CORRUT(I) = 2.*CORR(I)
        ENDIF
   50 CONTINUE
      IF(KONV(N).EQ.5) THEN
        CORRUT(N) = 5.*CORR(N)
      ELSE IF(KONV(N-1).EQ.5) THEN
        CORRUT(N) = 2.*CORR(N)
      ENDIF
C
      DO 60 I=1,N
        CORR(I) = CORRUT(I)
        KONV(I) = 0
   60 CONTINUE
C
      WRITE(7,65) (CORR(I),I=1,N)
   65 FORMAT(1X,10F7.1)
C
      RETURN
      END
C

****************************************************************

****************************************************************


C
      SUBROUTINE DETABS(J,JP,NTP,IOUTR)
      implicit real*8 (a-h,o-z)
C
C        THIS ROUTINE GIVES THE DETAILS OF THE ABSORPTION MECHANISMS.
C        CHANGES IN THE ABSORPTION-COEFFICIENT PROGRAM ARE EXPECTED TO
C        BE CONFINED TO THE TABLES AND TO THIS ROUTINE.
C        DETABS HAS TWO PURPOSES
C        1. JP=0   DETERMINATION OF WAVELENGTH-INDEPENDENT FACTORS (DEP. ON
C                  T, PE AND THE COMPONENT) STORED IN FAKT.
C        2. JP= THE ACTUAL WAVELENGTH NUMBER.
C                  MULTIPLICATION OF AB, COMPUTED IN SUBROUTINE ABSKO,
C                  BY WAVELENGTH-DEPENDENT FACTORS. SUMMATION OF THE TOTAL
C                  ABSORPTION AND SCATTERING COEFFICIENTS ( SUMABS AND
C                  SUMSCA ).
C
C        N O T E .  BEFORE A CALL ON DETABS FOR PURPOSE 1, SUBROUTINE
C        JON MUST HAVE BEEN CALLED.
C
C        IF J IS LESS THAN OR EQUAL TO ZERO, THE WEIGHT FOR A ROSSELAND MEAN
C        WILL BE COMPUTED AND STORED IN VIKTR (THE WEIGHT BEING 1/VIKTR).
C        NTP IS THE ARRAY INDEX OF THE T-PE POINT.
C        IF IOUTR IS GREATER THAN ZERO AT A CALL WITH JP GREATER THAN ZERO
C        (PART TWO OF THE ROUTINE), DETAILS OF THE ABSORPTION COEFFICIENTS
C        ARE PRINTED. IF IOUTR IS GREATER THAN ONE, A TABLE HEADING IS ALSO
C        PRINTED.
C
C
C        CONTENTS OF COMMON/CI5/, COMMUNICATING PHYSICAL INFORMATION FROM
C        SUBROUTINE JON.
C             ABUND  ABUNDANCES
C             ANJON  FRACTIONS OF IONIZATION
C             H      QUANTUM NUMBER OF THE HIGHEST EXISTING HYDROGENIC LEVEL
C             PART   PARTITION FUNCTIONS
C             DXI    DECREASE OF IONIZATION ENERGY OF HYDROGEN IN ELECTRON-VOLTS
C             F1     N(HI)/N(H)
C             F2     N(HII)/N(H)
C             F3     N(H-)/N(H)
C             F4     N(H2+)/N(H)
C             F5     N(H2)/N(H)
C             XKHM   'DISSOCIATION CONSTANT' OF H-
C             XMH    MASS OF THE HYDROGEN ATOM IN GRAMS
C             XMY    GRAMS OF STELLAR MATTER/GRAMS OF HYDROGEN
C
C        DIMENSIONS NECESSARY
C        ABUND(NEL),ANJON(NEL,MAX(NJ)),ELS(NT),H(5),HREST(NT),PART(NEL,MAX(NJ)),
C        PROV(NKOMP)
C        THE DIMENSIONS ARE LOWER LIMITS. DIMENSIONS IN COMMON /CA5/ ARE
C        COMMENTED ON IN SUBROUTINE ABSKO.
C        NEL IS THE NUMBER OF CHEMICAL ELEMENTS INITIATED IN SUBROUTINE INJON
C        NJ(I) IS THE NUMBER OF STAGES OF IONIZATION, INCLUDING THE NEUTRAL
C             STAGE, FOR ELEMENT I
C        NKOMP IS THE NUMBER OF COMPONENTS, NOT INCLUDING THOSE ADDED BY
C             ANALYTICAL EXPRESSIONS AFTER STATEMENT NO. 13.
C        NT   IS THE NUMBER OF TEMPERATURES-ELECTRON PRESSURES GIVEN AT THE
C             CALL OF SUBROUTINE ABSKO.
C
C
      include 'parameter.inc'
C
      DIMENSION ELS(NDP),HREST(NDP)
      DIMENSION FAKRAY(NDP)
      DIMENSION PHTVA(NDP),PHEL(NDP),H2RAY(NDP)
      COMMON /CLIN/lin_cia
      COMMON/CI5/abmarcs(17,ndp),ANJON(17,5),H(5),PART(17,5),
     *DXI,F1,F2,F3,F4,F5,XKHM,XMH,XMY(ndp)
      COMMON/CA2/RCA2DUM(4000),ICA2DUM(602),NKOMP
      COMMON/CA5/AB(30),FAKT(30),PE(NDP),T(NDP),XLA(20),XLA3(20),RO,
     *SUMABS,SUMSCA,VIKTR,ISET,NLB
      COMMON/UTPUT/IREAD,IWRIT
      COMMON /CARC4/ PROV(30),NPROVA,NPROVS,NPROV
      CHARACTER*8 SOURCE,ABNAME
      COMMON /CHAR/ ABNAME(30),SOURCE(30)
      COMMON /CMOLRAT/ FOLD(NDP,8),MOLOLD,KL
      COMMON /CCIATEST/ CIATEST(44,NDP)
      LOGICAL FIRST
      CHARACTER*8 NHMIN,NH2PR,NHEPR,NELS,NHRAY,NH2RAY
      DATA FIRST/.TRUE./
      DATA NH2PR/'H2PR'/,NHEPR/'HEPR'/,NELS/'ELSC'/,NHRAY/'H-SC'/
     & ,NH2RAY/'H2SC'/,NHMIN/'H-'/
C
C        SAVE ABSORPTION COMPONENT NAMES THE FIRST TIME DETABS IS CALLED
C
      IF (.NOT.FIRST) GO TO 1
      DO 2 KOMP=19,NKOMP
      ABNAME(1)=NHMIN
    2 ABNAME(KOMP-16)=ABNAME(KOMP)
      NPROVA=NKOMP-16+2
      ABNAME(NPROVA-1)=NH2PR
      ABNAME(NPROVA)=NHEPR
      NPROVS=3
      NPROV=NPROVA+NPROVS
      ABNAME(NPROV-2)=NELS
      ABNAME(NPROV-1)=NHRAY
      ABNAME(NPROV)=NH2RAY
      write(6,*) ' nkomp,nprova,nprovs = ',nkomp,nprova,nprovs
      write(6,*) ' the names of the abs and scat components are:'
      write(6,1233) (abname(npr),npr=1,nprov)
1233  format(3(/10a8))
      FIRST=.FALSE.
    1 CONTINUE
C
      TETA=5040./T(NTP)
      IF(JP.GT.0)GO TO 7
C      write(50,*) ' JP  XLA  Omega H2_CIA  He_CIA'
C
C        1. COMPUTATION OF WAVELENGTH-INDEPENDENT QUANTITIES
C
      HN=1./(XMH*XMY(ntp))
      HNH=F1*HN
C        H-
      FAKT(1)=PE(NTP)*HNH*1.E-17/XKHM
      FAKT(18)=PE(NTP)*HNH*2.E-26/PART(1,1)
C        HI
      TETA31=31.30364*TETA
      XFAKH=2.0898E-26/PART(1,1)*HNH
      NNIV=15
      XNIV=15.
      IF(H(1).LT.XNIV)NNIV=INT(H(1))
      DO3 M=1,NNIV
      XM2=M*M
      XM3=XM2*DFLOAT(M)
    3 FAKT(M+1)=XFAKH*EXP(-TETA31*(1.-1./XM2))/XM3
      FAKT(NNIV+1)=FAKT(NNIV+1)*MIN(H(1)-NNIV,1.0D+0)
      IF(NNIV.GE.15)GO TO 6
    4 N1=NNIV+1
      DO5 M=N1,15
    5 FAKT(M+1)=0.
C
C        FREE-FREE HI ABSORPTION
    6 UMC=2.3026*DXI*TETA
      EXPJ=XFAKH*EXP(-TETA31+UMC)/(2.*TETA31)
      ADDF=EXP(TETA31/((DFLOAT(NNIV)+0.5)**2)-UMC)-1.
      IF(H(1).LT.XNIV+0.5)ADDF=0.
      FAKT(17)=EXPJ
      HREST(NTP)=EXPJ*ADDF
C
C        H+H
      FAKT(19)=(HNH*1.E-25)*(HNH*1.E-25)*RO
C        H2+
      FAKT(20)=(HNH*1.E-20)**2*RO*ANJON(1,2)/ANJON(1,1)
C        H2-
      FAKT(21)=PE(NTP)*F5*HN
C        C I
      FAKT(22)=ANJON(3,1)*abmarcs(3,ntp)*HN*9./PART(3,1)
C        MG I
      FAKT(23)=ANJON(8,1)*abmarcs(8,ntp)*HN/PART(8,1)
C        AL I
      FAKT(24)=ANJON(9,1)*abmarcs(9,ntp)*HN*6./PART(9,1)
C        SI I
      FAKT(25)=ANJON(10,1)*abmarcs(10,ntp)*HN*9./PART(10,1)
C        HE I
      FAKT(26)=ANJON(2,1)*abmarcs(2,ntp)*HN/PART(2,1)
C        HE-
      FAKT(27)=PE(NTP)*ANJON(2,1)*abmarcs(2,ntp)*HN
C        ELECTRON SCATTERING
      ELS(NTP)=4.8206E-9*PE(NTP)/(T(NTP)*RO)
      PH2=F5*HN*RO
      PH2=PH2*1.38E-16*0.987E-6*273.
      PHTVA(NTP)=PH2*PH2/RO
      PHEL(NTP)=(abmarcs(2,ntp)/abmarcs(1,ntp))*RO*HN/1.008
      PHEL(NTP)=PHEL(NTP)*1.38E-16*0.987E-6*273.
      PHEL(NTP)=PHEL(NTP)*PH2/RO
      CIATEST(1,ntp) = anjon(2,1)
      CIATEST(2,ntp) = ro
      CIATEST(3,ntp) = hn
      CIATEST(4,ntp) = part(2,1)
      CIATEST(5,ntp) = ph2
      CIATEST(6,ntp) = phtva(ntp)
      CIATEST(7,ntp) = abmarcs(1,ntp)
      CIATEST(8,ntp) = abmarcs(2,ntp)
      PHELX=(abmarcs(2,ntp)/abmarcs(1,ntp))*RO*HN/1.008
      CIATEST(9,ntp) = phelx
      PH2X=F5*HN*RO
      CIATEST(10,ntp) = ph2x
      CIATEST(11,ntp) = phel(ntp)

C        RAYLEIGH SCATTERING
      FAKRAY(NTP)=HNH*2./PART(1,1)
      H2RAY(NTP)=F5*HN
      RETURN
C        N O T E . APART FROM VECTORS HREST AND ELS, NONE OF THE
C        TEMPERATURE- OR PRESSURE-DEPENDENT VARIABLES DEFINED ABOVE CAN
C        GENERALLY BE USED AT THE NEXT VISIT BELOW.
C        ANY SET OF FACTORS WHICH IS WANTED SHOULD BE STORED IN AN ARRAY WITH
C        DIMENSION = NT, LIKE HREST AND ELS, OR IN FAKT, WHERE THE DATA FOR
C        FURTHER USE ARE STORED IN SUBR. ABSKO.
C
C        2. WAVELENGTH-DEPENDENT FACTORS. SUMMATION.
C        CORRECTION FOR STIMULATED EMISSION
    7 EXPA=EXP(-28556.*TETA/XLA(JP))
   11 STIM=1.-EXPA
C
C        ABSORPTION
      SUMABS=0.
C        H I
      DO12 KOMP=2,17
      SUMABS=SUMABS+AB(KOMP)
   12 CONTINUE
      SUMABS=(SUMABS+HREST(NTP))*XLA3(JP)
      PROV(2)=SUMABS
C        H-
      HMIN=AB(1)+AB(18)/STIM
      SUMABS=SUMABS+HMIN
      PROV(1)=HMIN
C        H+H, H2+, HE I, C I, MG I, AL I, SI I
      DO13 KOMP=19,NKOMP
      SUMABS=SUMABS+AB(KOMP)
      PROV(KOMP-16)=AB(KOMP)
   13 CONTINUE
C
C        HERE FURTHER ABSORPTION MECHANISMS, GIVEN IN TABLES OR BY
C        ANALYTICAL EXPRESSIONS, MAY BE ADDED.
      OMEGA=1./XLA(JP)*1.E+8
      CALL H2OPAC(OMEGA,T(NTP),PROPAC)
      H2PRES=PROPAC*PHTVA(NTP)
      PROV(NPROVA-1)=H2PRES
      CALL HEOPAC(OMEGA,T(NTP),PROPAC)
      HEPRES=PROPAC*PHEL(NTP)
      PROV(NPROVA)=HEPRES
c	temporary!!!!!!!!!!!!!!!!!!!
C         if(ntp.eq.30)	write(50,2221) jp, xla(jp), omega, 
C    1                          h2pres, hepres
2221	format(i5,2f11.3,2e12.4)
      H2PRES=H2PRES/STIM
      PROV(NPROVA-1)=H2PRES
      HEPRES=HEPRES/STIM
      PROV(NPROVA)=HEPRES
C      if(lin_cia.eq.1 .and.(xla(jp).ge.1250..and.xla(jp).le.250000.))
C      if(lin_cia.eq.1 .and.(xla(jp).ge.5000..and.xla(jp).le.125000.))
      if(lin_cia.eq.1)
     *     SUMABS=SUMABS+H2PRES+HEPRES
      SUMABS=SUMABS*STIM
C
C        SCATTERING
      XRAY=MAX(XLA(JP),1026.0D+0)
      XRAY2=1./(XRAY*XRAY)
      RAYH=XRAY2*XRAY2*(5.799E-13+XRAY2*(1.422E-6+XRAY2*2.784))*
     *FAKRAY(NTP)
      RAYH2=XRAY2*XRAY2*(8.14E-13+XRAY2*(1.28E-6+XRAY2*1.61))*H2RAY(NTP)
      SUMSCA=ELS(NTP)+RAYH+RAYH2
      PROV(NPROV-2)=ELS(NTP)
      PROV(NPROV-1)=RAYH
      PROV(NPROV)=RAYH2
C
      IF(J.GT.0)GO TO 15
C
C        WEIGHT FOR A ROSSELAND MEAN
   14 VIKTR=EXPA/(STIM*STIM*(XLA3(JP)*1E-3)**2)
   15 CONTINUE
C
      IF(IOUTR-1)23,21,20
C
C        **** PRINT-OUT ****
   20 WRITE(IWRIT,200)XLA(JP),(ABNAME(KP),KP=1,NPROV)
  200 FORMAT(' WAVEL.=',F7.0,'    ABS       SCAT  ',6A6,/10A6)
   21 DO22 KP=1,NPROVA
   22 PROV(KP)=PROV(KP)/SUMABS*STIM
      DO24 KP=1,NPROVS
   24 PROV(NPROVA+KP)=PROV(NPROVA+KP)/SUMSCA
      if(lin_cia.ne.1) 
     *   WRITE(IWRIT,*) ' (Linskys CIA is not included in SUMABS)'
      WRITE(IWRIT,201)T(NTP),SUMABS,SUMSCA,(PROV(KP),KP=1,NPROV)
  201 FORMAT(' T=',F7.1,1X,1p2E10.3,0p6F7.4,//10F7.4)
   23 CONTINUE
      RETURN
      END


********************************************************************

*********************************************************************

C
      SUBROUTINE INITAB (IOUTS)
      implicit real*8 (a-h,o-z)
      include 'parameter.inc'
C
      CHARACTER MOLNAME*4,OSFIL*60,SAMPLING*3,INWNFIL*60
      NAMELIST /INPUTOS/ WNB,WNSTEP,WNEND
     *         ,INTVOS,nchrom,NOSMOL,MOLNAME,OSFIL,masabs
     *         ,losresl,osresl,wnos_first,wnos_last,kos_step
     *         ,LISTWN,INWNFIL,NEWC3
     *         ,newosatom,newosatomlist
      COMMON/COSWR/osresl,losresl,listwn
      COMMON/CNEWC3 /NEWC3
      COMMON /COSLIST/ WNB(25),WNSTEP(25),WNEND,INTVOS
      COMMON/COS/WNOS(NWL),CONOS(NDP,NWL),WLOS(NWL),WLSTEP(NWL)
     *    ,KOS_STEP,NWTOT,NOSMOL,NEWOSATOM,NEWOSATOMLIST
     *    ,nchrom,OSFIL(30),MOLNAME(30),SAMPLING
      COMMON/CA1/DUMA(120),IDUM(240),NEXTT,NUTZT
      COMMON/CFIL/IRESET(10),ISLASK,IREAT
      COMMON/CXLSET/XL(20,10),NSET,NL(10)
      COMMON /CROS/WROS(20)
      COMMON/UTPUT/IREAD,IWRIT
      COMMON/COUTR/NTO,NTPO(10)
      COMMON/CVAAGL/XLA(500),W(500),NLB
      COMMON/CLINE4/ILINE
      common /cmasabs/ masabs(3)
C
C LOGICAL UNITS
      ISLASK=11
      IREAT=9
      ISAVE=IWRIT
      DO 1 I=1,10
 1    IRESET(I)=10
C
C ZEROSET WEIGHTS
      DO 4 I=1,500
 4    W(I)=0.
      DO 5 I=1,20
 5    WROS(I)=0.
C
C WAVELENGTH SETS
C ROSSELAND
      CALL VAAGL(NLBRO,XL,WROS)
      NL(1)=NLBRO
      CALL VAAGL(NLB,XLA,W)
      LMAX=15
      IFIRST=2
      ILAST=10
      CALL SETDIS(NLB,XLA,LMAX,IFIRST,ILAST)
C 8000. $NGSTR@M STANDARD
      NL(1)=NL(1)+1
      NL1=NL(1)
      XL(NL1,1)=8000.
      IF(IOUTS.LT.0) IWRIT=4
      CALL INABS(IOUTS)
      IF(IOUTS.GE.-1)IWRIT=ISAVE
      NL(1)=NL(1)-1
C
C CONTROL PARAMETTERS
      NTO=0
      NEXTT=0
      NUTZT=0
C
C  Calculate the wave numbers, WN, for OS opacity and the total number, 
C  NWTOT, of OS frequency points used. This wavenumber scale should be 
C  identical to the one used when the OS-tables were created.
C
      READ(5,INPUTOS)
C     WRITE(6,INPUTOS)
      IF(NOSMOL.GT.MAXOSMOL) STOP ' Increase dimension for MAXOSMOL'
C

C  Calculate (or read, if listwn=1) the wave numbers, WN, for OS opacity
C  If losresl = 1, the os wavenumbers are computed based on a specified
C  and the total number, NWTOT, of OS frequency points for the OS-tables.
C  spectral resolution, osresl. Else it is computed in fixed steps inside a
C  number of prespecified intervals.
C
      IF (LISTWN .eq. 1) go to 225   ! read an existing OS - wn list


C  compute an OS - wn list:

      if (losresl.eq.0) go to 228     ! compute list with fixed steps

C  use OS list with fixed resolution, osresl, through spectrum:
      wnos(1) = wnos_first
      step = 1.d0 + 1.d0/osresl
      do 240 k = 1,nwl-1
      wnkj = wnos(k)
      do 242 kj = 1,kos_step
242   wnkj = wnkj * step
      nwtot = k + 1
      wnos(nwtot) = wnkj
      if (wnos(nwtot) .gt. wnos_last) go to 241
240   continue
C we come here only if dimension for the OS is too small:
      wnos1 = wnos_first
      step = 1.d0 + 1.d0/osresl
      do 740 k = 1,100000
      wnkj = wnos1
      do 742 kj = 1,kos_step
742   wnkj = wnkj * step
      nwtot = k + 1
      wnos1 = wnkj
      if (wnos1 .gt. wnos_last) go to 741
740   continue
741   continue
      print*,' given,needed dimensions -nwl,nwtot='
     &         ,nwl,nwtot
      stop ' error: increase dimension nwl in parameter.inc for wnos '
241   continue

      write(6,*) ' We did a resolution based set of os wavenumbers'
      write(6,245) nwtot,osresl/dfloat(kos_step)
     & ,kos_step*step,wnos(1),wnos(nwtot)
     & ,1.e4/wnos(nwtot),1.e4/wnos(1)
245   format(' Total ',i6,' OS wavenumbers for Marcs radiative transf.'
     & ,/' Approximate resolution is',f8.0,' (corresponding to ',
     & 'step factor',f8.5,')'
     & ,/' OS interval:',f7.1,'-',f9.1,' cm^-1 (=',f5.2,'-',f5.1,'mu)')


      go to 226

228   continue
      L=0
      WNOS(1)=WNB(1)
      WNB(INTVOS+1) = WNEND
      DO 100 I = 1,INTVOS
      WNLAST = WNB(I+1)-WNSTEP(I)
102   L=L+1
      L1=L+1
      WNOS(L1)=WNOS(L)+WNSTEP(I)
      IF ( WNOS(L1).LE.WNLAST ) GO TO 102
100   CONTINUE
      NWTOT = L1

      go to 226

225   continue

C     read an existing OS - wn list:

      open (unit=54,file=inwnfil,status='old')
      do 200 l=1,100000
      lbp = 10912-l+1
      read(54,*,end=201) wnos(lbp),idum
      nwx = l
200   continue
201   continue
      nwtot = nwx
      im = 0
      do 210 i=1,nwx
      if (wnos(i).lt.wnb(1)) then
       im = i
       go to 210
      end if
        iu = i-im
        wnos(iu) = wnos(i)
        nwtot = iu
      if (wnos(iu).gt.wnend) then
       nwtot = iu-1
       go to 211
      end if
210   continue
211   continue
      close(unit=54)


226   continue

      DO 103 I=1,NWTOT
      L=NWTOT-I+1
103   WLOS(L)=1.D8/WNOS(I)
      DO 104 L=2,NWTOT-1
      LP=L+1
      LM=L-1
104   WLSTEP(L) = ( WLOS(LP)-WLOS(LM) ) / 2.
      WLSTEP(1) = WLOS(2)-WLOS(1)
      WLSTEP(NWTOT) = WLOS(NWTOT)-WLOS(NWTOT-1)


      IF (NWTOT.GT.NWL)  STOP ' DIMENSION NWL TOO SMALL'
      WRITE(6,12) NWTOT,WNOS(1),WNOS(NWTOT)
12    FORMAT (' Marcs OS wavenumbers established:',/
     &    ' There will be ',I7,' frequency points from wnos(1)=',
     &    f9.2,' to wnos(nwtot)=',f9.2,' in radiative transfer')

C
C The wavelengths for calculation of the line-opacity must be inside
C the vavelengths where the continuum-opacity is calculated.
C The first and last line-opacity wavelength are WLOS(1) and 
C WLOS(NWTOT), respectively.
C The first and last continuum-opacity wavelength are XL(1,2) and 
C XL(NL(NLB),NLB), respectively.
C
      IF (WLOS(1) .LE. XL(1,2)) THEN
           PRINT*, ' Error-message from INITAB.FOR:'
           PRINT*, ' WLOS(1) = ',WLOS(1)
           PRINT*, ' XL(1,2) = ',XL(1,2)
           STOP ' WLOS(1) < first cont. point'
      END IF
      IF (WLOS(NWTOT) .GT. XL(NL(NSET),NSET)) THEN
           PRINT*, ' Error-message from INITAB.FOR:'
           PRINT*, 'NLB = ',NLB
           PRINT*, 'NSET = ',NSET
           PRINT*, 'NL(NSET) = ',NL(NSET)
           PRINT*, 'XL(NL(NSET),NSET) = ',XL(NL(NSET),NSET)
           STOP ' WLOS(last) > last cont. point'
      END IF
C
      RETURN
      END
C
      SUBROUTINE INITJN(IOUTS)
      implicit real*8 (a-h,o-z)
      include 'parameter.inc'
C
      COMMON /UTPUT/IREAD,IWRIT
      common /cdustdata/ dabstable(max_lay,nwl),dscatable(max_lay,nwl),
     *    eps(max_eps,max_lay),temp(max_lay),pgas(max_lay),n_lay,idust
C
C INITJN INITIATES THE JON BLOCK FROM LUN 9
C
      ISAVE=IREAD
      IREAD=9              !unit=9 is jonabs.dat
*      CALL DUMIN
      CALL INJON(IOUTS)
      IREAD=ISAVE
      RETURN
      END
C
      SUBROUTINE INJON(IOUTS)
      implicit real*8 (a-h,o-z)
      include 'parameter.inc'
C
C
C        THIS ROUTINE READS DATA NECESSARY FOR THE COMPUTATION OF IONIZATION
C        EQUILIBRIA ETC. (IN SUBR. JON).
C        1. NEL= THE NUMBER OF CHEMICAL ELEMENTS CONSIDERED.
C           A=   THE RATIO OF THE NUMBER OF HYDROGEN NUCLEI TO THE NUMBER OF
C                NUCLEI OF METALLIC ELEMENTS.
C           NMET=THE NUMBER OF METALLIC ELEMENTS IN THE LIST OF CHEMICAL
C                ELEMENTS CONSIDERED. THE LAST NMET ELEMENTS OF THE LIST ARE
C                CONSIDERED TO BE METALLIC, FOR THE CALCULATION OF THE
C                QUANTITY A (DEFINED ABOVE).
C        2. IEL  IS THE ARRAY WHICH WILL CONTAIN THE SYMBOLS FOR THE CHEMICAL
C                ELEMENTS CONSIDERED.
C           ABUND IS THE ARRAY WHICH WILL CONTAIN THE PREVAILING ABUNDANCES
C                THE CHEMICAL ELEMENTS CONSIDERED AT INPUT. THESE ABUNDANCES
C                ARE EXPRESSED ON A LOGARITHMIC SCALE (BASE 10) AND NEED NOT BE
C                NORMALIZED. THE ABUNDANCES ARE MODIFIED IN THIS SUBROUTINE
C                SO THAT THE RIGHT VALUE OF A (DEFINED ABOVE) IS OBTAINED.
C        3. AI   IS THE ARRAY WHICH WILL CONTAIN THE ATOMIC WEIGHTS OF THE
C                ELEMENTS CONSIDERED.
C        4. DATA FOR THE COMPUTATION OF THE PARTITION FUNCTIONS IS READ NEXT
C           NJ(I)= THE NUMBER OF STAGES OF IONIZATION CONSIDERED FOR ELEMENT I.
C         FOR EACH STAGE OF IONIZATION JA THE FOLLOWING QUANTITIES ARE READ
C           G0(JA)=THE STATISTICAL WEIGHT OF THE GROUND LEVEL,
C           NK(JA)=THE NUMBER OF ELECTRON CONFIGURATIONS CONSIDERED.
C         FOR EACH ELECTRON CONFIGURATION JB THE FOLLOWING QUANTITIES ARE READ
C           XION(JB)=THE IONIZATION ENERGY IN ELECTRON VOLTS,
C           G2(JB)=THE STATISTICAL WEIGHT (2L+1)*(2J+1)
C           XL(JB)=THE LOWEST QUANTUM NUMBER OF THE ASYMPTOTIC (HYDROGENIC) PART
C                OF THE PARTITION FUNCTION,
C           NL(JB)=THE NUMBER OF TERMS IN THE (APPROXIMATE) EXPRESSION FOR THE
C                'MIDDLE PART' OF THE PARTITION FUNCTION ('QPRIME').
C           ALFA IS AN ARRAY WHICH WILL CONTAIN THE 'STATISTICAL WEIGHTS' OF
C                THE (APPROXIMATE) EXPRESSIONS FOR THE 'MIDDLE PARTS' OF THE
C                PARTITION FUNCTIONS.
C           GAMMA IS AN ARRAY CONTAINING THE CORRESPONDING 'EXCITATION
C                POTENTIALS' (EXPRESSED IN ELECTRON VOLTS).
C         FOR THE METHOD USED SEE TRAVING ET AL., ABH. HAMB. VIII, I (1966).
C        5. ELEMENTS AND STAGES OF IONIZATION THAT SHOULD BE DISREGARDED ARE
C           INDICATED BY IELEM(I)=0 FOR ELEMENT I AND BY ION(I,J)=0 FOR
C           IONIZATION STAGE J. OTHERWISE INDICATORS SHOULD BE = 1.
C        6. NQFIX IS THE NUMBER OF PARTITION FUNCTIONS THAT SHOULD BE CONSTANT.
C                THE VALUES ARE READ INTO THE VECTOR PARCO AND AN INDICATION IS
C                MADE IN IQFIX.  IQFIX(I,J)=0 MEANS THAT THE PARTITION FUNCTION
C                FOR ELEMENT I, STAGE OF IONIZATION J, IS CONSIDERED TO BE
C                CONSTANT.
C           NQTEMP IS THE NUMBER OF PARTITION FUNCTIONS  THAT SHOULD BE
C                PRESSURE-INDEPENDENT AND INTERPOLATED IN T. VALUES OF FOUR
C                TEMPERATURES (TPARF, THE SAME FOR ALL ELEMENTS) AND
C                CORRESPONDING PARTITION FUNCTIONS (PARF) ARE READ. IQFIX(I,J)=1
C                MEANS THAT A PRESSURE-INDEPENDENT PARTITION FUNCTION FOR INTER-
C                POLATION IN T IS GIVEN.
C        7. IFISH IS A PARAMETER FOR THE CHOICE OF THE ASYMPTOTIC PARTITION 
C                FUNCTION. IFISH=0 MEANS THAT THE ASYMPTOTIC PART WILL BE EVALU-
C                ATED FOLLOWING BASCHEK ET AL., ABH. HAMB. VIII,26 (1966). IFISH
C                =1 MEANS THAT IT WILL BE EVALUATED FOLLOWING FISCHEL AN SPARKS
C                ASTROPHYS. J. 164, 356 (1971).
C        8. TMOLIM IS THE HIGHER TEMPERATURE LIMIT BEYOND WHICH MOLECULES WILL
C                NOT BE CONSIDERED
C
C        MOREOVER SOME INITIATING WORK IS DONE FOR SUBR. JON. UNLOGARITHMIC
C        ABUNDANCES ARE NORMALIZED ON HYDROGEN, XMY AND SUMH (DEFINED BELOW)
C        ARE COMPUTED AND SOME FURTHER QUANTITIES ARE EVALUATED AT THE END.
C        A DETAILED PRINTOUT IS GIVEN IF IOUTS IS EQUAL TO ONE. AFTER INJON
C        HAS BEEN CALLED ONCE, A NEW DETAILED PRINTOUT IS OBTAINED IF
C        INJON IS CALLED WITH IOUTS GREATER THAN ONE.
C
C        DIMENSIONS NECESSARY
C        ABUND(NEL),AI(NEL),ALFA(LMAX),ANJON(NEL,MAX(NJ)),FL2(5),F1Q(3),
C        GAMMA(LMAX),G0(JMAX),G2(KMAX),H(5),IEL(NEL),IELEM(NEL),
C        ION(NEL,MAX(NJ)),IQFIX(NEL,MAX(NJ)),JAMEM(NEL,MAX(NJ)),JBBEG(JMAX)
C        JCBEG(JMAX),NJ(NEL),NK(JMAX),NL(KMAX),PARCO(JMAX),PARF(4*JMAX),
C        PARPP(4),PARPT(4),PARQ(4*JMAX),PART(NEL,MAX(NJ)),SHXIJ(5),TPARF
C        XION(KMAX),XIONG(NEL,MAX(NJ)),XL(KMAX)
C        THE DIMENSIONS ARE LOWER LIMITS
C        JMAX IS THE TOTAL NUMBER OF STAGES OF IONIZATION, INCLUDING NEU
C             ATOMS.
C        KMAX IS THE TOTAL NUMBER OF ELECTRON CONFIGURATIONS.
C        LMAX IS THE TOTAL NUMBER OF TERMS IN THE (APPROXIMATE) EXPRESSI
C             FOR THE MIDDLE PART OF THE PARTITION FUNCTIONS ('QPRIME'),
C             ACCORDING TO TRAVING ET AL., CITED ABOVE.
C        NEL  IS THE NUMBER OF CHEMICAL ELEMENTS.
C        NJ(I) IS THE NUMBER OF STAGES OF IONIZATION, INCLUDING THE NEUT
C             STAGE, FOR ELEMENT I.
C
      DIMENSION AI(16),F1Q(3),F2Q(2),PARF(180),PARPP(4),PARPT(4)
      DIMENSION JAMEM(16,5)
      COMMON/CI1/FL2(5),PARCO(45),PARQ(180),SHXIJ(5),TPARF(4),
     *XIONG(16,5),EEV,ENAMN(ndp),SUMH(ndp),XKBOL,NJ(16),IEL(16),
     *SUMM(ndp),NEL
      COMMON/CI9/AI
      COMMON/CI3/ALFA(300),GAMMA(300),G0(45),G2(80),XION(80),XL(80),
     *JBBEG(45),JCBEG(45),NK(45),NL(80),IFISH
      COMMON/CI4/ TMOLIM,IELEM(16),ION(16,5),MOLH,JUMP
      COMMON/CI5/abmarcs(17,ndp),ANJON(17,5),H(5),PART(17,5),
     *DXI,F1,F2,F3,F4,F5,XKHM,XMH,XMY(ndp)
      COMMON/CI6/TP,IQFIX(16,5),NQTEMP
      COMMON/UTPUT/IREAD,IWRIT
      common /tsuji/ nattsuji,nmotsuji,parptsuji(500),abtsuji(17,ndp)
      INTEGER MOLH, JUMP
      character sunz*1
      common/cabinit/abinit(natms),kelem(natms),nelem
      common /statec/ppr(ndp),ppt(ndp),pp(ndp),gg(ndp),zz(ndp),dd(ndp),
     *  vv(ndp),ffc(ndp),ppe(ndp),tt(ndp),tauln(ndp),ro(ndp),ntau,iter
      dimension mx_elm(17),abundatms_inp(natms),sum(ndp),fakt(ndp)
      common /cdustdata/ dabstable(max_lay,nwl),dscatable(max_lay,nwl),
     *    eps(max_eps,max_lay),temp(max_lay),pgas(max_lay),n_lay,idust
      namelist /abundances/sunz,zscale,abundatms_inp
      data mx_elm /1, 2,6,7,8,10,11,12,13,14,16,19,20,23,25,27,21/
C                  H He C N O Ne Na Mg Al Si  S  K Ca Cr Fe Ni Ti

C
      ISAVE=IWRIT
      IWRIT=6
C
        IF(IOUTS.GT.1) GO TO 25
C
C        READING OF THE ABUNDANCES AND THEI ASSOCIATED QUANTITIES
C        **** 1 ****
C
C     READ(5,100) NEL,A,NMET        !->A and NMET seems to be always 0!
       NEL = 16
       A = 0.
       NMET = 0
C
C        **** 2 ****
C
C
      READ(IREAD,110)(IEL(I),I=1,NEL)    !the 16 element names H - Ni
C     READ(5,101) (ABUND(I),I=1,17)      !element #17 is Ti
C The abundances are now read from elabund.dat i gem_init, and can be
C adjusted collectively or individually in input ($abundances). The 
C 17 abundances usually read here from input file are brought into injon
C and other routines from gem_init by common CI5.
      open(unit=2,file='data/elabund.dat',status='old',readonly)
      do i=1,natms
        read(2,*,iostat=io) kelem(i), abinit(i)
        if(io .ne. 0) exit
      end do
      nelem = i-1
      write(6,131) abinit(1),abinit(2),(abinit(i),i=6,8),abinit(25)
131   format('Abundances of H,He,C,N,O,Fe:',/8f7.2)
      read(5,abundances)
      if(sunz.ne.'y' .and. sunz.ne.'Y') then
        print *, 'Error: cannot scale abundances in this version'
        stop
      end if
      write(6,*) ' abund(1-17):'
      write(6,101) (abinit(i),i=1,17)
      
      do i=1,17
        abmarcs(i,1:ntau) = abinit(mx_elm(i))
        abtsuji(i,1:ntau) = abmarcs(i,1:ntau)
      end do
      
!      abmarcs(3,30:ntau) = 8.43
!      abmarcs(4,30:ntau) = 7.83
!      abmarcs(5,30:ntau) = 8.69
!      abtsuji = abmarcs
C
C        **** 3 ****
C
      READ(IREAD,102)(AI(I),I=1,NEL)      !the atomic weight 1.0,4.0,...,58.7
      NU=NEL
      SUM(1:ntau)=0.
      SUMM(1:ntau)=0.
      FAKT(1:ntau)=1.
C
C        THE ABUNDANCES ARE CONVERTED FROM A LOGARITHMIC SCALE TO A DIRECT
C        SCALE, AND ARE THEN NORMALIZED ON HYDROGEN. XMY=GRAMS OF STELLAR MATTER
C        /GRAMS OF HYDROGEN. SUMH=NUMBER OF OTHER NUCLEI/NUMBER OF HYDROGEN
C        NUCLEI.
C        SUMM=NUMBER OF NUCLEI OTHER THAN H, C, N, O / NUMBER OF HYDROGEN
C
      if(idust .eq. 0) then
      if(nmet.le.0)go to 22      !->which seems to be always the case
      nu=nel-nmet+1
      do1 i=nu,nel
      abmarcs(i,1:ntau)=10.**abmarcs(i,1:ntau)
    1 sum(1:ntau)=abmarcs(i,1:ntau)+sum(1:ntau)
      abmarcs(17,1:ntau)=10.**abmarcs(17,1:ntau)

      fakt(1:ntau)=sum(1:ntau)*a/10.**abmarcs(1,1:ntau)
      nu=nu-1
   22 do2 i=1,nu
      abmarcs(i,1:ntau)=10.**abmarcs(i,1:ntau)*fakt(1:ntau)  !->fakt==1., so this is just sum abund
    2 sum(1:ntau)=sum(1:ntau)+abmarcs(i,1:ntau)
      abmarcs(17,1:ntau)=10.**abmarcs(17,1:ntau)
      xmy(1:ntau)=0.
      aha=abmarcs(1,1)
      do3 i=1,nel
      abmarcs(i,1:ntau)=abmarcs(i,1:ntau)/aha
      summ(1:ntau)=summ(1:ntau)+abmarcs(i,1:ntau)
    3 xmy(1:ntau)=xmy(1:ntau)+abmarcs(i,1:ntau)*ai(i)   !AI(I)=atomic weight, so XMY=#AU/H_nuclei
      abmarcs(17,1:ntau)=abmarcs(17,1:ntau)/aha
      xmy(1:ntau)=xmy(1:ntau)/ai(1)            !AI(1) = 1.008, so on AU scale.
      sumh(1:ntau)=sum(1:ntau)/aha-1.
      summ(1:ntau)=summ(1:ntau)-abmarcs(1,1:ntau)-abmarcs(3,1:ntau)-
     &  abmarcs(4,1:ntau)-abmarcs(5,1:ntau)
      end if
C
C        **** 4 ****
C
C        READING OF DATA FOR THE PARTITION FUNCTIONS.
C        FOR THE SYMBOLS, SEE ABOVE.
C      
      READ(IREAD,103)(NJ(I),I=1,NEL)
      JA=1
      JB=1
      JC1=1
      DO11 I=1,NEL
      NJP=NJ(I)
      DO11 J=1,NJP
      JAMEM(I,J)=JA
      JBBEG(JA)=JB
      JCBEG(JA)=JC1
C        JBBEG AND JCBEG ARE INDICATORS USED BY FUNCTION QTRAV
C
      READ(IREAD,104)G0(JA),NK(JA)
      NKP=NK(JA)
      IQFIX(I,J)=2
C        IQFIX(I,J)=2 MEANS THAT A 'FULL' PARTITION FUNCTION SHOULD BE
C        COMPUTED. THIS MAY BE CHANGED UNDER **** 7 ****.
C
      JA=JA+1
      DO11 K=1,NKP
      READ(IREAD,105)XION(JB),G2(JB),XL(JB),NL(JB)
      IF(K.GT.1)GO TO 9
      XIONG(I,J)=XION(JB)
C        XIONG IS THE IONIZATION ENERGY IN ELECTRON VOLTS FOR THE GROUND STATE,
C        USED IN THE COMPUTATION OF IONIZATION EQUILIBRIA IN SUBROUTINE JON.
C
    9 CONTINUE
      JC2=NL(JB)+JC1-1
      JBM=JB
      JB=JB+1
      IF(NL(JBM).LE.0)GO TO 10
      READ(IREAD,106)(GAMMA(L),ALFA(L),L=JC1,JC2)
   10 JC1=JC2+1
   11 CONTINUE
C
C        **** 5 ****
C
C        READING OF THE INDICATORS OF THE ELEMENTS AND THE STAGES OF IONIZATION
C        TO BE DISREGARDED.
      DO12 I=1,NEL
      NJP=NJ(I)
      READ(IREAD,107)IELEM(I),(ION(I,J),J=1,NJP)
   12 CONTINUE
C
C        **** 6 ****
C
C        SPECIFICATION OF THOSE PARTITION FUNCTIONS GIVEN AS CONSTANTS.
C        INDICATION IN IQFIX.
      READ(IREAD,103)NQFIX
      IF(NQFIX.LE.0)GO TO 15
   13 DO14 I=1,NQFIX
      READ(IREAD,109)I1,J1,PARCOP
      JA=JAMEM(I1,J1)
      PARCO(JA)=PARCOP
   14 IQFIX(I1,J1)=0
   15 CONTINUE

C
C        SPECIFICATION OF THOSE PARTITION FUNCTIONS TO BE INTERPOLATED IN T.
C        INDICATION IN IQFIX.
      READ(IREAD,103)NQTEMP
      IF(NQTEMP.EQ.0)GO TO 20
      READ(IREAD,101)TPARF
      DO17 I=1,NQTEMP
      READ(IREAD,109)I1,J1,(PARPP(K),K=1,4)
      IQFIX(I1,J1)=1
C
C        PREPARATION FOR INTERPOLATION OF PARTITION FUNCTIONS IN T (CONCLUDED
C        IN SUBROUTINE JON).
      DO16 K=1,3
   16 F1Q(K)=(PARPP(K+1)-PARPP(K))/(TPARF(K+1)-TPARF(K))
      DO161 K=1,2
  161 F2Q(K)=(F1Q(K+1)-F1Q(K))/(TPARF(K+2)-TPARF(K))
      F3Q=(F2Q(2)-F2Q(1))/(TPARF(4)-TPARF(1))
      PARPT(1)=PARPP(1)
      PARPT(2)=F1Q(1)
      PARPT(3)=F2Q(1)
      PARPT(4)=F3Q
      JA=JAMEM(I1,J1)
      DO17 K=1,4
      JK=(JA-1)*4+K
      PARQ(JK)=PARPT(K)
   17 PARF(JK)=PARPP(K)
C        PARQ IS IN COMMON/CI1/ AND IS USED IN SUBROUTINE JON. PARF IS JUST
C        USED BELOW.
C
   20 CONTINUE
C
C        **** 7, 8 ****
C
C        THE PARAMETERS IFISH AND TMOLIM. INITIATING WORK FOR SUBROUTINE JON.
C        WHEN MOLH IS GREATER THAN ZERO THE MOLECULAR FORMATION WILL BE
C        IN SUBR. MOLEQ (ONLY H2 AND H2+), ELSE MORE COMPLETE MOLECULAR
C        FORMATION WILL BE EVALUATED IN SUBR. MOL.
C        --> not H too, because it is an atom !!
C
      READ(IREAD,100)IFISH

C TMOLIM, MOLH are given in JONABS.DAT
C here is set : IREAD=9 --> = JONABS.DAT

      READ(IREAD,4528) TMOLIM,MOLH
 4528 FORMAT(F10.0,I5)
C      write(6,*) 'MOLH after read from JONABS.DAT = ', MOLH
C      write(6,*) 'TMOLIM after read  from JONABS.DAT = ',TMOLIM
C      write(6,*) 'JUMP after read  from JONABS.DAT = ',jump
      DO21 J=1,5
      FLJ=J
      FL2(J)=31.321*FLJ*FLJ
   21 SHXIJ(J)=SQRT(13.595*FLJ)
C
C        EEV=THE ELECTRON VOLT (EXPRESSED IN TERMS OF ERGS)
C        XMH=THE MASS OF THE HYDROGEN ATOM (EXPRESSED IN GRAMS)
C        XKBOL=BOLTZMANN'S CONSTANT (EXPRESSED IN ERGS PER KELVIN
      EEV=1.602095E-12
      XMH=1.67339E-24
      XKBOL=1.38053E-16
      if(idust.eq.0) enamn(1:ntau)=eev/(xmh*xmy(1:ntau))
      TP=0.
C        TP IS THE TEMPERATURE AT THE 'PRECEDING' CALL OF JON.
C
C        **** PRINT-OUT ****
C
* here: iwrit = 6 = mxms7.out
C      WRITE(IWRIT,201)
C      WRITE(IWRIT,202)
C      WRITE(IWRIT,203)
      DO33 I=1,NEL
      NJP=NJ(I)
C      WRITE(IWRIT,204) IEL(I),ABUND(I),IELEM(I),
C     &                 (ION(I,J),IQFIX(I,J),J=1,NJP)
   33 CONTINUE
      WRITE(IWRIT,207)
      WRITE(IWRIT,208)
      IF(IFISH.EQ.1)WRITE(IWRIT,211)
      IF(IFISH.EQ.0)WRITE(IWRIT,210)
      IF(NQTEMP.GT.0.OR.NQFIX.GT.0)WRITE(IWRIT,214)
      IF(NQTEMP.GT.0)WRITE(IWRIT,209)TPARF
      JA=1
      DO32 I=1,NEL
      NJP=NJ(I)
      DO32 J=1,NJP
      JP=J-1
C      IF(IQFIX(I,J).EQ.0)WRITE(IWRIT,205)IEL(I),JP,PARCO(JA)
      JK1=(JA-1)*4+1
      JK2=(JA-1)*4+4
C      IF(IQFIX(I,J).EQ.1)WRITE(IWRIT,206)IEL(I),JP,(PARF(JK),JK=JK1,JK2)
   32 JA=JA+1
      IF(NQTEMP.GT.0)WRITE(IWRIT,215)
      WRITE(IWRIT,212)TMOLIM
      IF(MOLH.LE.0) WRITE(IWRIT,216)
      IF(MOLH.GT.0) WRITE(IWRIT,217)
      WRITE(IWRIT,213)XMY(1),SUMH(1)


      IF(IOUTS.LE.0)GO TO 40
   25 CONTINUE

   40 CONTINUE
C
      IWRIT=ISAVE
C
      RETURN
C
  100 FORMAT(I10,F10.4,I10)
  101 FORMAT(6F10.4)
  102 FORMAT(6F10.4)
  103 FORMAT(12I5)
  104 FORMAT(F5.0,I5)
  105 FORMAT(F6.3,F4.0,F5.1,I5)
  106 FORMAT(4(F10.3,F10.4))
  107 FORMAT(I10,5I5)
  108 FORMAT(2F10.4)
  109 FORMAT(2I5,4F10.4)
  110 FORMAT(16A3)
  201 FORMAT(1H1,'D A T A  F R O M  S U B R O U T I N E  I N J O N')
  202 FORMAT(1H0,30X,1HI,14X,2HII,13X,3HIII,12X,2HIV)
  203 FORMAT(1H ,'    ABUNDANCE   IELEM      ION   PF       ION   PF  ',
     *'      ION   PF       ION   PF')
  204 FORMAT(1H ,A2,E12.4,I5,5X,5(I5,I5,5X))
  205 FORMAT(1H ,A2,', STAGE OF IONIZATION=',I2,
     *' PARTITION FUNCTION (CONSTANT)=',F10.3)
  206 FORMAT(1H ,A2,', STAGE OF IONIZATION=',I2,' PART.FUNC. (T-DEP.) ='
     *,4F10.3)
  207 FORMAT(1H0,'IELEM AND ION = 1 OR 0 MEANS ELEMENT AND IONIZATION',
     *' STAGE SHOULD BE CONSIDERED OR DISREGARDED RESP.')
  208 FORMAT(1H0,'PF=2 FULL PART. FUNC., =1 PART. FUNC. TO BE',
     *' INTERPOLATED IN T, =0 CONSTANT PART. FUNC.')
  209 FORMAT(1H ,'T-DEPENDENT PARTITION FUNCTIONS GIVEN FOR T =  ',4F10.
     *0)
  210 FORMAT(1H,
     *'ASYMPTOTIC PARTS OF PART. FUNC. FOLLOWING BASCHEK ET AL. 1966')
  211 FORMAT(1H,'ASYMPTOTIC PARTS OF PART. FUNC. FOLLOWING FISCHEL AND',
     *' SPARKS 1971')
  212 FORMAT(1H0,'MOLECULES CONSIDERED BELOW T=',F7.0,' KELVIN')
  213 FORMAT(1H0,'XMY=GRAMS STELLAR MATTER/GRAMS OF HYDROGEN=',F7.4,5X,
     *'SUMH=NUMBER OF OTHER ATOMS/NUMBER OF H=',F8.5)
  214 FORMAT(1H0,'PARTITION FUNCTIONS SUPPLIED BY THE USER')
  215 FORMAT(1H ,'IF T OUTSIDE RANGE FOR INTERPOLATIONS DETAILED PART.',
     *' FUNCTIONS ARE COMPUTED')
  216 FORMAT(1H0,'MOLECULES CONSIDERED:  H2, H2+, H2O, OH, CH, CO, CN,',
     *' C2, N2, O2, NO, NH')
  217 FORMAT(1H0,'MOLECULES CONSIDERED:  H2, H2+')
      END
C
      SUBROUTINE INP3(X,Y,XINT,YINT)
      implicit real*8 (a-h,o-z)
C
C        NEWTONINTERPOLATION, TREPUNKTS
C        OBS *** INGEN SPAERR MOT EXTRAPOLATION *****
C
      DIMENSION X(3),Y(3),F1(2)
      DO1 K=1,2
    1 F1(K)=(Y(K+1)-Y(K))/(X(K+1)-X(K))
      F2=(F1(2)-F1(1))/(X(3)-X(1))
      YINT=Y(1)+(XINT-X(1))*F1(1)+(XINT-X(1))*(XINT-X(2))*F2
      RETURN
      END
C
      SUBROUTINE JON(T,PE,IEPRO,PG,RO,E,IOUTR)
      implicit real*8 (a-h,o-z)
C
C
C        THIS ROUTINE COMPUTES IONIZATION EQUILIBRIA FOR A GIVEN TEMPERATURE 
C        (T, EXPRESSED IN KELVIN) AND A GIVEN ELECTRON PRESSURE (PE, IN
C        DYNES PER CM2). THE FRACTIONS OF IONIZATION ARE PUT IN THE ANJON VECTOR
C        AND THE PARTITION FUNCTIONS ARE PUT IN PART. IF IEPRO IS GREATER THAN
C        ZERO, THE GAS PRESSURE (PG,IN DYNES PER CM2), DENSITY (RO, IN GRAMS
C        PER CM3) AND INNER ENERGY (E, IN ERGS PER GRAM) ARE ALSO EVALUATED.
C        N O T E . RADIATION PRESSURE IS NOT INCLUDED IN E.
C
C        THE ENERGIES OF IONIZATION ARE REDUCED BY DXI, FOLLOWING BASCHEK ET 
C        AL., ABH. HAMB. VIII, 26 EQ. (10). THESE REDUCTIONS ARE ALSO MADE IN
C        THE COMPUTATION OF E.
C        THE ENERGY OF DISSOCIATION FOR H- HAS BEEN REDUCED BY 2*DXI, FOLLOWING
C        TARAFDAR AND VARDYA, THIRD HARV. SMITHS. CONF., PAGE 143. THE FORMATION
C        OF MOLECULES IS CONSIDERED FOR T LESS THAN TMOLIM.
C
C        IF IOUTR IS GREATER THAN ZERO, A DETAILED PRINT-OUT WILL BE GIVEN.
C
C
C        THE FUNCTION  QTRAV AND SUBROUTINE MOLEQ ARE CALLED.
C        THEY CALL QAS AND MOLFYS RESPECTIVELY.
C
C        DIMENSIONS NECESSARY
C        A(5),DQ(4),F(MAX(NJ)),PFAK(MAX(NJ)),RFAK(JMAX)
C        DIMENSIONS OF ARRAYS IN COMMONS /CI1/,/CI4/,/CI5/ AND /CI6/ ARE
C        COMMENTED ON IN SUBROUTINE INJON.
C        JMAX IS THE TOTAL NUMBER OF STAGES OF IONIZATION, INCLUDING NEUTRAL
C             ATOMS.
C        NJ(I) IS THE NUMBER OF STAGES OF IONIZATION, INCLUDING THE NEUTRAL
C             STAGE, FOR ELEMENT I.
C
C
      include 'parameter.inc'
C
      DIMENSION DQ(4),F(5),PFAK(5),RFAK(45)
      COMMON/CI1/FL2(5),PARCO(45),PARQ(180),SHXIJ(5),TPARF(4),
     *XIONG(16,5),EEV,ENAMN(ndp),SUMH(ndp),XKBOL,NJ(16),IEL(16),
     *SUMM(ndp),NEL
      COMMON/CI4/ TMOLIM,IELEM(16),ION(16,5),MOLH,JUMP
      COMMON/CI5/abmarcs(17,ndp),ANJON(17,5),H(5),PART(17,5),
     *DXI,F1,F2,F3,F4,F5,XKHM,XMH,XMY(ndp)
      COMMON/CI6/TP,IQFIX(16,5),NQTEMP
      COMMON/CI7/A(5),PFISH,ITP
      COMMON/UTPUT/IREAD,IWRIT
      COMMON/RABELL/XXRHO(NDP),XYRHO
      COMMON/CI8/YYPG,YYRHO,YYE
      COMMON/CMOL1/EH,FE,FH,FHE,FC,FCE,FN,FNE,FO,FOE,FK,FKE,FS,FSE
     *             ,FT,FTE
      COMMON/CMOL2/PK(33),NMOL
      COMMON/CARC3/F1P,F3P,F4P,F5P,HNIC,PRESMO(33)
      common /fullequilibrium/ partryck(ndp,maxmol),
     &  xmettryck(ndp,maxmet),xiontryck(ndp,maxmet)
      COMMON/CPHYDRO/PHYDRO
      COMMON /CMOLRAT/ FOLD(NDP,8),MOLOLD,KL
      COMMON /PJONINF/ P_MOL(NDP), P_NEU_HCNO(NDP), P_ION_HCNO(NDP),
     & P_NEU_HE(NDP),P_ION_HE(NDP), P_NON_HHECNO(NDP), PG_JON(NDP), 
     & HN_JON(NDP), RO_JON(NDP), P6_JON(NDP)

      INTEGER MOLH, JUMP
C      REAL PHYDRO
C        STATEMENT FUNCTION FOR 10.**
      EXP10(X)=EXP(2.302585*X)
C
      ITP=1
C
C        IS T=THE TEMPERATURE OF THE PRECEDING CALL
      IF(ABS((T-TP)/T).LT.1.E-8)GO TO 53
   51 ITP=0
C
C        SOME QUANTITIES, ONLY DEPENDENT ON T
      TETA=5040./T
      TETA25=1.202E9/(TETA*TETA*SQRT(TETA))
      DO52 J=1,5
   52 A(J)=FL2(J)*TETA
C        A=ALFA(BASCHEK ET AL., CITED ABOVE)
C
      IF(NQTEMP.EQ.0)GO TO 53
C
C        PREPARATION FOR INTERPOLATION OF PARTITION FUNCTIONS IN T
      DQ(1)=1.
      DQ(2)=T-TPARF(1)
      DQ(3)=DQ(2)*(T-TPARF(2))
      DQ(4)=DQ(3)*(T-TPARF(3))
C
C        SOME QUANTITIES ALSO DEPENDENT ON PE
C        THE PFAK FACTORS ARE USED IN THE SAHA EQUATION. H(J) IS THE
C        QUANTUM NUMBER OF THE CUT OF THE PARTITION FUNCTIONS (ACCORDING
C        TO BASCHEK ET AL., CITED ABOVE) FOR J-1 TIMES IONIZED ATOMS. H IS
C        USED IN QAS.
C
C        XNEL= THE ELECTRON (NUMBER) DENSITY (PER CM3)
C        PFISH= P(FISCHEL AND SPARKS, ASTROPHYS. J. 164, 359 (1971)) IS USED IN
C        FUNCTION QAS.
C
   53 DXI=4.98E-4*TETA*SQRT(PE)
      DUM=TETA25/PE
      DIM=EXP10(DXI*TETA)
      PFAK(1)=DIM*DUM
      SQDXI=1./SQRT(DXI)
      H(1)=SHXIJ(1)*SQDXI
      DO54 J=2,5
      PFAK(J)=PFAK(J-1)*DIM
   54 H(J)=SHXIJ(J)*SQDXI
      XNEL=PE/(XKBOL*T)
      PFISH=4.2E3/XNEL**0.166666667
C
C        PARTITION FUNCTIONS AND IONIZATION EQUILIBRIA
C
      XNECNO=0.
      XNENH=0.
      EJON=0.
      JA=1
C
C        BEGINNING OF LOOP OVER ELEMENTS ('THE I-LOOP').
      DO24 I=1,NEL
      NJP=NJ(I)
C
C        SHOULD ELEMENT NO. I BE CONSIDERED
      IF(IELEM(I).GT.0)GO TO 9
      DO 55 J=1,NJP
      ANJON(I,J)=0.
      PART(I,J)=0.
   55 CONTINUE
      GO TO 23
C
C        BEGINNING OF LOOP OVER STAGES OF IONIZATION ('THE J-LOOP')
    9 DO19 J=1,NJP
      JM1=J-1
C
C        SHOULD STAGE OF IONIZATION NO. J BE CONSIDERED
      IF(ION(I,J).GT.0)GO TO 10
      ANJON(I,J)=0.
      PART(I,J)=0.
      GO TO 18
C
C        WHICH KIND OF PARTITION FUNCTION SHOULD BE COMPUTED
C
   10 IF(IQFIX(I,J)-1)14,11,13
   11 IF(T.LT.TPARF(1).OR.T.GT.TPARF(4))GO TO 13
      PARTP=PART(I,J)
      IF(ITP.GT.0)GO TO 15
C
C        PARTITION FUNCTIONS TO BE INTERPOLATED IN T
      JPARF=(JA-1)*4+1
      PARTP=0.
      DO12 IP=1,4
      PARTP=PARTP+PARQ(JPARF)*DQ(IP)
   12 JPARF=JPARF+1
      GO TO 15
C
C        PARTITION FUNCTIONS FOLLOWING TRAVING ET AL., ABH. HAMB. VIII,1 (1966)
   13 PARTP=QTRAV(TETA,H(J),J,JA)
      GO TO 15
C
C        THE PARTITION FUNCTION IS CONSTANT
   14 PARTP=PARCO(JA)
   15 PART(I,J)=PARTP
   
C
C        IONIZATION EQUILIBRIA AND TOTAL NUMBER OF ELECTRONS
C
      IF(J.LE.1)GO TO 19
      IF(ITP.GT.0)GO TO 17
      RFAK(JA)=EXP10(-XIONG(I,JM1)*TETA)
   17 F(JM1)=PFAK(JM1)*RFAK(JA)*PARTP/PART(I,J-1)
      GO TO 19
   18 IF(J.GT.1)F(JM1)=0.
   19 JA=JA+1
C        END OF 'THE J-LOOP'
C
      FIL=1.
      DO20 J=2,NJP
      LL=NJP-J+1
   20 FIL=1.+F(LL)*FIL
      ANJON(I,1)=1./FIL
      XNEN=0.
      DO21 J=2,NJP
      JM1=J-1
      ANJON(I,J)=ANJON(I,JM1)*F(JM1)
      IF(I.LE.1)GO TO 24
      FLJM1=JM1
   21 XNEN=ANJON(I,J)*FLJM1+XNEN
      IF(I.GT.2.AND.I.LT.6) XNECNO=XNECNO+XNEN*abmarcs(I,kl)
      XNENH=XNEN*abmarcs(I,kl)+XNENH
C        XNENH=NUMBER OF ELECTRONS FROM ELEMENTS OTHER THAN HYDROGEN (Q IN
C        MIHALAS, METH. COMP. PHYS. 7, 1 (1967), EQ. (35))
C        XNECNO=NUMBER OF ELECTRONS FROM ELEMENTS OTHER THAN H, C, N, O
CUGJ9902:XNECNO=NUMBER OF ELECTRONS FROM C, N, O (I>2.and.I<6) !!!!
C
C
C        COMPUTATION OF THE ENERGY OF IONIZATION (EJON). HYDROGEN IS NOT
C        INCLUDED.
C
      XERG=0.
C        XERG= THE ENERGY OF IONIZATION PER ATOM (IN ELECTRON VOLTS)
C
      DO22 J=2,NJP
      JM1=J-1
      FLJM1=JM1
   22 XERG=ANJON(I,J)*(XIONG(I,JM1)-DXI*FLJM1)+XERG
      EJON=XERG*abmarcs(I,kl)+EJON
      GO TO 24
   23 JA=JA+NJP
   24 CONTINUE
C        END OF 'THE I-LOOP'
C
C
      XNECNO=XNENH-XNECNO
      TP=T
      IF(IEPRO.LE.0)GO TO 71
C
C        COMP. OF PRESSURE, DENSITY AND INNER ENERGY
C
      XIH=XIONG(1,1)-DXI
      XIHM=0.747-2.*DXI
C        XIH AND XIHM ARE THE ENERGIES OF IONIZATION FOR H AND H- RESPECTIVELY
C        (IN ELECTRON VOLTS).
C
      XKHM=TETA25*2.*EXP10(-TETA*XIHM)
C        XKHM = THE 'DISSOCIATION CONSTANT' FOR H-.
C
      HJONH=ANJON(1,2)/ANJON(1,1)
C

*********
C 15.12.94 Ch.Helling
C
C JUMP is set in JONABS.DAT JUMP = 1 => MOLH = 1 :
C                                       call only MOLEQ
C                                       = only H2, H2+ as molecules cosidered
C                                                                            
C( gives the possibility to use only TSUJI-code instaed of MOL for the molecular
C  equilibrium )
C
C                           Jump = 0 : call MOL
C                                      = all molecules cosidered in old MARCS-routine
C                                        MOL
C                                        ( there are problems to calculate lower 
C                                          temperatures )
C
C 25/1/01, UGJ:
C JUMP = 2: JFF routine corresponding to Tsuji's method
C JUMP = 3: JFF Gibs minimalisation method (892 molecules, ions and atoms).
C JUMP = 4: Added GGchem code by ERC
*********


      IF(T.GT.TMOLIM)GO TO 42
C     IF(JUMP.eq.1) MOLH=1
* the former step is necessary because often is MOLH=0 on other places !
      IF(MOLH.LE.0) GO TO 45


*********
C
C => if MOLH < or = 0 than is MOLEQ allways " overjumped " ( !! )
C                     and a more complete moleculare formation will
C                     be evaluated in MOL
C => if MOLH > 0  then only MOLEQ is supposed to work and only the 
C           (= 1) moleculare formation of H2 and H2+ is considered 
C                 ( => total hydrogen pressure [H, H2, H2+] 
C                   => possibility to use only Tsujis eqilibrium routines EQMOL,
C                      called in TEST_TSUJI )
*********
C
C        FORMATION OF MOLECULES. ONLY H2 AND H2+
   41 CALL MOLEQ(T,PE,HJONH,XIH,XKHM,XIHM,XNENH,F1,F2,F3,F4,F5,FE,FSUM,
     *   EH)

*
* F1 = N(HI)/N(H)   ==> F1P = P(HI)
* F2 = N(HII)/N(H)  ==> F2P = P(HII)
* F3 = N(H-)/N(H)   ==> F3P = P(H-)
* F4 = N(H2+)/N(H)  ==> F4P = P(H2+)
* F5 = N(H2)/N(H)   ==> F5P = P(H2)
*
*      print*, 'I am after MOLEQ'
      FEPE=PE/FE
      F1P=F1*FEPE
      F3P=F3*FEPE
      F4P=F4*FEPE
      F5P=F5*FEPE
      PHYDRO=FSUM*PE/FE
      GO TO 43
C        FORMATION OF MOLECULES COMPOSED OF H,C,N,O
   45 IF(ANJON(3,1).LE.0..OR.ANJON(4,1).LE.0..OR.ANJON(5,1).LE.0.)
     * GOTO 41
*      print*,' now I am in MOL ,  jumped over MOLEQ'
      HJONC=ANJON(3,2)/ANJON(3,1)
      HJONN=ANJON(4,2)/ANJON(4,1)
      HJONO=ANJON(5,2)/ANJON(5,1)
      ABUC=abmarcs(3,kl)/abmarcs(1,kl)
      ABUN=abmarcs(4,kl)/abmarcs(1,kl)
      ABUO=abmarcs(5,kl)/abmarcs(1,kl)
C      write(7,*) ' i jon'
C      write(7,*) ' T,PE,HJONH,HJONC,HJONN,HJONO,ABUC,ABUO,ABUN,XIH'
C     *,' XKHM,XIHM,XNECNO,F1,F2,F3,F4,F5 = '
C      write(7,*) T,PE,HJONH,HJONC,HJONN,HJONO,ABUC,ABUO,ABUN,XIH,XKHM
C     *,XIHM,XNECNO,F1,F2,F3,F4,F5
C      write(7,*) ' now call mol'
      CALL MOL(T,PE,HJONH,HJONC,HJONN,HJONO,ABUC,ABUO,ABUN,XIH,XKHM,XIHM
     *,XNECNO,F1,F2,F3,F4,F5)
      SUMPMO=0.
      PRESMO(1)=FHE*PK(1)
      PRESMO(2)=FHE*FHE*PK(2)
      PRESMO(3)=FHE*FHE*HJONH*PK(3)
      PRESMO(4)=FHE*FHE*FOE*PK(4)
      PRESMO(5)=FHE*FOE*PK(5)
      PRESMO(6)=FHE*FCE*PK(6)
      PRESMO(7)=FCE*FOE*PK(7)
      PRESMO(8)=FCE*FNE*PK(8)
      PRESMO(9)=FCE*FCE*PK(9)
      PRESMO(10)=FNE*FNE*PK(10)
      PRESMO(11)=FOE*FOE*PK(11)
      PRESMO(12)=FNE*FOE*PK(12)
      PRESMO(13)=FNE*FHE*PK(13)
      PRESMO(14)=FCE*FCE*FHE*FHE*PK(14)
      PRESMO(15)=FHE*FCE*FNE*PK(15)
      PRESMO(16)=FCE*FCE*FHE*PK(16)
      PRESMO(17)=0.0
      PRESMO(18)=FHE*FSE*PK(18)
      PRESMO(19)=FKE*FHE*PK(19)
      PRESMO(20)=FCE*FCE*FCE*FHE*PK(20)
      PRESMO(21)=FCE*FCE*FCE*PK(21)
      PRESMO(22)=FCE*FSE*PK(22)
      PRESMO(23)=FKE*FCE*PK(23)
      PRESMO(24)=FKE*FCE*FCE*PK(24)
      PRESMO(25)=FNE*FSE*PK(25)
      PRESMO(26)=FKE*FNE*PK(26)
      PRESMO(27)=FKE*FOE*PK(27)
      PRESMO(28)=FSE*FOE*PK(28)
      PRESMO(29)=FSE*FSE*PK(29)
      PRESMO(30)=FKE*FSE*PK(30)
      PRESMO(31)=FTE*FOE*PK(31)
      PRESMO(32)=FTE*FOE*FOE*PK(32)
      PRESMO(33)=FTE*FCE*FCE*PK(33)
      DO 30 I=1,NMOL
      PRESMO(I)=PRESMO(I)*PE
   30 SUMPMO=SUMPMO+PRESMO(I)
C SUMPMO = sum partial pressures [dyn/cm2] of all molecules
C SUMPA = sum partial pressures of neutral atoms H,C,N,O not in molecules
C SUMPI = sum partial pressure from ionized H,C,N,O
C SUMM=SUMM-ABUND(1)-ABUND(3)-ABUND(4)-ABUND(5) in INJON is number
C of nuclei other than H,C,N,O per H atom (abund are normalized to abund(1) ).
C PE/FE = ro*k_boltz*T/(XMY*XMH) = conversion from #atoms/#H to P[dyn/cm3]
C PE*SUMM/FE = sum partial pressures from all nuclei not HCNO.
      SUMPA=PE*(FHE+FCE+FNE+FOE)
      SUMPI=PE*(FHE*HJONH+FCE*HJONC+FNE*HJONN+FOE*HJONO)
      HNIC=PE*FHE
      HPNIC=HNIC*HJONH
      PG=PE+SUMPMO+SUMPA+SUMPI+PE*SUMM(kl)/FE
      P_MOL(KL) = SUMPMO
      P_NEU_HCNO(KL) = SUMPA
      P_ION_HCNO(KL) = SUMPI
      P_NEU_HE(KL) = 
     &   PE*abmarcs(1,kl)*abmarcs(2,kl)/FE*ANJON(2,1)
      P_ION_HE(KL) = 
     &   PE*abmarcs(1,kl)*abmarcs(2,kl)/FE*ANJON(2,2)
      PHE = P_ION_HE(KL) + P_NEU_HE(KL)
      P_NON_HHECNO(KL) = PE*SUMM(kl)/FE-PHE
      PG_JON(KL) = PG
      HN_JON(KL) = 1./(XMH*XMY(kl))
      P6_JON(KL) = HNIC+0.42*PHE+0.85*PRESMO(2)
      GOTO 46
C
C        NO MOLECULES
   42 F2=ANJON(1,2)
      FE=XNENH+F2
      F1=ANJON(1,1)
      F3=0.
      F4=0.
      F5=0.
      FSUM=1.
      EH=-XIH*F1


**********
C
C ==>
C
C PG = PE + PE*FSUM/FE + PE*SUMH/FE   ( = line 43 )
C    = PE + PHYDRO + PE*SUMH/FE
C
C SUMH = number of other nuclei / number of hydrogen nuclei
C   FE = PE/PH   with PH = NH*kT ( NH number of hydrogen nuclei per cm3 )
C   => FE = number of other nuclei * kT
C         = ficticious pressure of the other nuclei, except hydrogen and e-
C
C   => PE*SUMH/FE = gaspressure contributed by all the atomic species
C
C   PHYGRO = P(HI) + P(HII) + P(H-) + P(H2+) + P(H2) 
C          = P(H)  + P(H+)  + P(H-) + P(H2+) + P(H2)
C 
**********

 43   CONTINUE 

*      print*,' PE           : ',PE
*      print*,'XKBOLZ, T : ', XKBOL,T

*       PHEI=(10**(abund(1) + abund(2)))*XKBOL*T
*       print*,' PE*SUMH/FE = PHEI ? ', PE*SUMH/FE, log10(PHEI)  
*   no it isn't !!
* PHEI partial pressure og He

*      print*,'successful jump around MOL, one line after 43 -> PG= ',PG
46    continue
C
C     write(6,*) ' t,pe,yypg,yyrho,yye,eh,ejon,enamn,p6_jon(kl),phe = '
C     write(6,1342) t,pe,YYPG,YYRHO,YYE,EH,EJON,ENAMN,p6_jon(kl),phe

      RO=PE*XMY(kl)*(XMH/XKBOL)/(FE*T)

      IF (JUMP.GE.1 .and. MOLH.eq.1) THEN

       PG=PE*(1.+(FSUM+SUMH(kl))/FE)
           XNHE = abmarcs(2,kl) / (XMH*XMY(kl))
           PHE = ro * 1.38053e-16 * T * XNHE
       P6_JON(KL) = xmettryck(kl,1)+0.42*PHE+0.85*partryck(kl,2)
C      write(6,1351) kl,p6_jon(kl),xmettryck(kl,1),PHE,partryck(kl,2)
C1351  format('p6 etx in jon,jump1:',i3,1p5e11.3)
      P_NEU_HE(KL) = PE*abmarcs(1,kl)*abmarcs(2,kl)/FE*ANJON(2,1)/
     &  (ANJON(2,1)+ANJON(2,2))
      P_ION_HE(KL) = PE*abmarcs(1,kl)*abmarcs(2,kl)/FE*ANJON(2,2)/
     &  (ANJON(2,1)+ANJON(2,2))
      PHE = P_ION_HE(KL) + P_NEU_HE(KL)
      PG_JON(KL) = PG
      HN_JON(KL) = 1./(XMH*XMY(kl))
C     write(6,*) 
C    & kl,p6_jon(kl),p_neu_he(kl),p_ion_he(kl),pg_jon(kl),hn_jon(kl)
 
      ENDIF


      XYRHO=RO
      E=1.5*PG/RO+(EH+EJON)*ENAMN(kl)
      YYPG=PG
      YYRHO=RO
      YYE=E
      RO_JON(KL) = RO
C partial pressure of He is put into presmo(17)
        XNHE = abmarcs(2,kl) / (XMH*XMY(kl)) 
        PRESMO(17) = 1.38053d-16 * T * XNHE * RO

C      write(6,*) ' again: t,pe,yypg,yyrho,yye,eh,etc:'
C      write(6,1342) t,pe,YYPG,YYRHO,YYE,EH,EJON,ENAMN,p6_jon(kl),phe
C1342  format(f8.0,1p9e12.5)

      IF(IOUTR.LE.0)GO TO 71
C
C        **** PRINT-OUT ****
C
* iwrit = 6 = mxms7.out
      WRITE(IWRIT,204)T,PE,PG,RO,E
      WRITE(IWRIT,201)
      WRITE(IWRIT,202)
      DO93 I=1,NEL
      NJP=NJ(I)
      WRITE(IWRIT,203)IEL(I),abmarcs(I,1),(ANJON(I,J),J=1,NJP)
      WRITE(IWRIT,207)(PART(I,J),J=1,NJP)
   93 CONTINUE
      IF(T.GT.TMOLIM)GO TO 44
      IF(MOLH.LE.0)GOTO 47
      WRITE(IWRIT,205)F1P,F3P,F5P,F4P
      GO TO 71
   47 CONTINUE

***************18.12.94 Ch.H

      IF (JUMP.EQ.0) THEN
       WRITE(IWRIT,208)HNIC,(PRESMO(I),I=1,13)
       ELSE
       WRITE(IWRIT,209) xmettryck(kl,1),(partryck(kl,I),I=1,2),
     &  (PARTRYCK(kl,I),I=4,13)
      ENDIF

************************

      GOTO 71
   44 WRITE(IWRIT,206)
C
   71 CONTINUE
c
c      print*,'end of jon. fo,fe,foe ',fo,fe,foe
c
      RETURN
  201 FORMAT(1H0,'ELEMENT  ABUNDANCE  IONIZATION FRACTIONS',17X,
     *'PARTITION FUNCTIONS')
  202 FORMAT(1H ,23X,1HI,7X,2HII,6X,3HIII,5X,2HIV,12X,1HI,9X,2HII,8X,
     *3HIII,7X,2HIV)
  203 FORMAT(6H      ,A2,E12.4,4F8.4)
  204 FORMAT(3H0T=,F7.1,5X,3HPE=,E12.4,5X,3HPG=,E12.4,5X,3HRO=,E12.4,
     *5X,2HE=,E12.4)
  205 FORMAT(1H0,'PARTIAL PRESSURES'/4X,'H',8X,'H-',7X,'H2',7X,'H2+'/1X,
     *4(1PE9.2))
  206 FORMAT(1H0,'NO MOLECULES CONSIDERED IN MARCS-ROUTINES -- ', 
     *' Tsuji routine for molecules used')
  207 FORMAT(1H+,56X,4E10.3) 
  208 FORMAT(1H0,'PARTIAL PRESSURES'/4X,'H',8X,'H-',7X,'H2',7X,'H2+',6X,
     *'H2O',6X,'OH',7X,'CH',7X,'CO',7X,'CN',7X,'C2',7X,'N2',7X,'O2',7X,
     *'NO',7X,'NH'/1X,14(1PE9.2))
 209  FORMAT(1HO,'PARTIAL PRESSURES'/4X,'H',8X,'H-',7X,'H2',6X,
     *'H2O',6X,'OH',7X,'CH',7X,'CO',7X,'CN',7X,'C2',7X,'N2',7X,'O2',7X,
     *'NO',7X,'NH'/1X,13(1PE9.2))
      END
C

**********************************************************************

**********************************************************************

C
      SUBROUTINE OPAC(J,X,S)
      implicit real*8 (a-h,o-z)
C
C        THIS ROUTINE ADMINISTERS COMPUTATION OF OPACITIES.
C        J IS THE WAVELENGTH NUMBER (IN THE XL ARRAY), X AND S ARE THE
C        NORMALIZED ABSORPTION AND SCATTERING COEFFICIENTS, RESPECTIVELY
C
      include 'parameter.inc'
C
      DIMENSION PRXC(NDP),PRXO(NDP),PRXW(NDP),PRXT(NDP),PRSC(NDP)
      DIMENSION ABSK(NDP),SPRID(NDP),ABSK1(NDP),SPRID1(NDP)
      DIMENSION X(NDP),S(NDP)
      DIMENSION V(NDP),CON(NDP)
      EQUIVALENCE (V(1),CON(1))
      COMMON /TAUC/TAU(NDP),DUMT(NDP),JJTAU
      COMMON /STATEC/PPR(NDP),PPT(NDP),PP(NDP),GG(NDP),ZZ(NDP),DD(NDP),
     &      VV(NDP),FFC(NDP),PE(NDP),T(NDP),TAULN(NDP),RO(NDP),NTAU,ITER
      COMMON /ROSSC/XKAPR(NDP),CROSS(NDP)
      COMMON/CXLSET/XL(20,10),NSET,NL(10)
      COMMON /CVAAGL/XLB(500),W(500),NLB
      COMMON/CLINE1/XLINLO,XLINUP,TSKAL(30),
     &             PESKAL(30),IPEBEG(30),IPEEND(30),LINUN,NTSKAL,NPSKAL
      COMMON/CLINE2/FACTOR(NDP,2,2),IT(NDP),IPE1(NDP),IPE2(NDP)
       COMMON/CLINE4/ILINE
      COMMON /CXMAX/XMAX
      COMMON /ODFCD/  CONODF(NDP,8)
      COMMON /CPF/PF,PFE,PFD,FIXROS,ITSTOP
C
      COMMON /TIO/ PTIO(NDP),ROSAV(NDP),POXG1(NDP)
      COMMON /CMOL1/DMUDMU(9),FOE,XMUDMUD(6)
      COMMON /CI4/dumdum,IDUMDUM(96),MOLH,JUMP
      COMMON /CARC3/ F1P,F3P,F4P,F5P,HNIC,PRESMO(33)
      COMMON /DENSTY/ BPZ(NDP),PRH2O(NDP)
C      COMMON/COPPR/oppr(15,3,120,3),jvxmax,itxmax  !15mol,10dpt,100wn
      INTEGER MOLH, JUMP
C      
! Dust
      common /cdustdata/ dabstable(max_lay,nwl),dscatable(max_lay,nwl),
     *    eps(max_eps,max_lay),temp(max_lay),pgas(max_lay),n_lay,idust      
      common /cdustopac/ dust_abs(ndp,nwl), dust_sca(ndp,nwl)
C
      LOGICAL PF,PFE,PFD,FIXROS,ITSTOP,FIRST
C                                                              
      DIMENSION PELOG(NDP)
c
      CHARACTER MOLNAME*4,OSFIL*60,SAMPLING*3
      COMMON/COS/WNOS(NWL),CONOS(NDP,NWL),WLOS(NWL),WLSTEP(NWL)
     *    ,KOS_STEP,NWTOT,NOSMOL,NEWOSATOM,NEWOSATOMLIST
     *    ,nchrom,OSFIL(30),MOLNAME(30),SAMPLING
C      COMMON/COPPRR/ xconop(120,10),xlineop(120,10)    !100wn,10dpt
C      COMMON/CONLIN/rconop(nwl,ndp),rlineop(nwl,ndp)
      COMMON/COPsum/ SSUM(NDP),XSUM(NDP),CONSUM(NDP)
c
       DATA JJ/4/
       DATA FIRST/.TRUE./
C      EXP10(A)=EXP(2.302585*A)
C
      IF (J.EQ.1) THEN
         IF (WLOS(1).LT.XL(1,2)) THEN
           PRINT*,' WLOS(1),XL(1,2),XL(2,2) = ',WLOS(1),XL(1,2),XL(2,2)
           STOP ' WLOS(1) (==1.E8/WNEND) < first continuums point'
         END IF
         IF (NOSMOL.GT.0) THEN
            write(66,*) 'OPAC J=1 before OSTABLOOK; ntau=',ntau
            write(66,*) 't{1-ntau}, pe{1.ntau}, pg{1-ntau}:'
            write(66,1266) (t(imo),imo=1,ntau)
            write(66,1267) (pe(imo),imo=1,ntau)
            write(66,1267) ((pp(imo)-ppr(imo)-ppt(imo)),imo=1,ntau)
1266        format(10f8.1)
1267        format(1p7e11.4)

            CALL OSTABLOOK

            write(66,*) ' now in OPAC for J=1 after OSTABLOOK'
         END IF
      END IF
c
C ROSSELAND OPACITY OR NOT
       JTAU=JJTAU
C
C MURIEL EST GENTILLE
C
C***********************************************************************
C SPECIAL
      IF (J.EQ.1) THEN
      DO 1234 K=1,JTAU
      PRXC(K)=0.0
      PRXO(K)=0.0
      PRXT(K)=0.0
      PRXW(K)=0.0
      PRSC(K)=0.0
 1234 CONTINUE
      ENDIF
C END OF SPECIAL
C
      IF(J.EQ.1) JJ=4
      IF(J.NE.1)GO TO 9
      ISWITCH=1
CUGJ      IF(ILINE.LE.0) GO TO 1
C      DO 2 K=1,JTAU
C2     PELOG(K)=log10(PE(K))
C      CALL DUBINT(NTSKAL,TSKAL,NPSKAL,PESKAL,IPEBEG,IPEEND,JTAU,T,
C     &            PELOG,FACTOR,IT,IPE1,IPE2)
C1     CONTINUE
C
C        COMPUTATION OF CONTINUOUS ABSORPTION COEFFICIENTS AND INTERPOLATION
      NEWT=1
      JMEM=1
      IMEM=2
C   ******** HERE WE ASSUME THAT THE FIRST SET IS USED FOR ROSSELAND MEAN
      JMEM1=2
      IMEM1=2
      CALL ABSKO(NEWT,JTAU,T,PE,IMEM,JMEM,ABSK,SPRID)
      NEWT=0
      CALL ABSKO(NEWT,JTAU,T,PE,IMEM1,JMEM1,ABSK1,SPRID1)
c
    9 IF(WLOS(J).LE.XL(JMEM1,IMEM1))GO TO 11
c
C    9 IF(XLB(J).LE.XL(JMEM1,IMEM1))GO TO 11
C      IF(ILINE.GT.0.AND.XLB(J).GE.XLINLO.AND.XLB(J).LE.XLINUP)
C    *       READ(LINUN)XXXX,YYYY
C      IF(ILINE.GT.0.AND.XLB(J).GE.XLINLO.AND.XLB(J).LE.XLINUP)
C    *       WRITE(7,7651) XXXX,YYYY
C        NEW COMPUTATION OF CONTINOUOS ABSORPTION COEFFICIENT
      JMEM=JMEM1
      IMEM=IMEM1
      JMEM1=JMEM1+1
      IF (IMEM.GT.10 .OR. JMEM1.GE.20) 
     *                STOP ' WNOS(1) < last continuums point '
      DO8 K=1,JTAU
      ABSK(K)=ABSK1(K)
    8 SPRID(K)=SPRID1(K)
      IF(JMEM1.LE.NL(IMEM1))GO TO 10
      JMEM1=1
      IMEM1=IMEM1+1
   10 CALL ABSKO(NEWT,JTAU,T,PE,IMEM1,JMEM1,ABSK1,SPRID1)
      GO TO 9
C        INTERPOLATION
11    CONTINUE
C      IF ( (J.EQ.1 .OR. J.EQ.NWTOT) .AND. FIRST) THEN
C       WRITE(6,*) 'J,WLOS(J), XL(J,I), XL(J1,I1) = '
C       WRITE(6,*) J,WLOS(J),XL(JMEM,IMEM),XL(JMEM1,IMEM1)
C      END IF
      DO12 K=1,JTAU
      DIFXL=(WLOS(J)-XL(JMEM,IMEM))/(XL(JMEM1,IMEM1)-XL(JMEM,IMEM))
C      DIFXL=(XLB(J)-XL(JMEM,IMEM))/(XL(JMEM1,IMEM1)-XL(JMEM,IMEM))
      X(K)=(ABSK(K)+(ABSK1(K)-ABSK(K))*DIFXL)/XKAPR(K)
      PRXC(K)=X(K)
      S(K)=(SPRID(K)+(SPRID1(K)-SPRID(K))*DIFXL)/XKAPR(K)
      PRSC(K)=S(K)
   12 CONTINUE  
         
C
C        COMPUTATION OF LINE-ABSORPTION COEFFICIENTS
C
      IF (J.EQ.1) THEN
      DO 19 K=1,JTAU
          CONSUM(K) = 0.
          XSUM(K) = 0.
          SSUM(K) = 0.
19    CONTINUE
      END IF
C
      DO 468 K=1,JTAU
          CONSUM(K) = CONSUM(K) + CONOS(K,J)/XKAPR(K)
          XSUM(K) = XSUM(K) + X(K)
          SSUM(K) = SSUM(K) + S(K)

C        if(k.eq.7) then
C          itx=1
C        else if (k.eq.27) then
C          itx = 2
C        else if (k.eq.47) then
C          itx = 3
C        else
C          itx=4
C        end if

C         jvw = nwtot-j+1
C         if(jvw/100*100.eq.jvw   .and. itx.le.3) then
C           jvx = max(1,jvw/100)
C
C           xconop(jvx,itx) = x(k)*XKAPR(K)
C           xlineop(jvx,itx) = conos(k,j)
C          end if

C          rconop(j,k) = xkapr(k)*x(k)
C          rlineop(j,k) = conos(k,j)

          X(K)=X(K)+CONOS(K,J)/XKAPR(K)
468   CONTINUE
C
C

! Computation of dust absorption & scattering
      if(idust .eq. 1) then
        do k=1,jtau
          x(k) = x(k)+dust_abs(k,j)/xkapr(k)
          s(k) = s(k)+dust_sca(k,j)/xkapr(k)
        end do
      end if


   25 CONTINUE
      DO 17 K=1,JTAU
17    X(K)=X(K)*XMAX/(X(K)+XMAX)
50    FORMAT(' X',6(' ***',1PE12.5))
C
C
      FIRST = .FALSE.

      RETURN
      END
C
C
      SUBROUTINE PEMAKE(T,PE,PG,PEX)
      implicit real*8 (a-h,o-z)
C
C 'PEMAKE-R' IS COMPATIBLE WITH 'PEMAKE' AND SOMEWHAT FASTER. A MODIFIED
C REGULA-FALSI PROCEDURE IS USED ON THE LOG-LOG PG-PE RELATION.
C 76.03.08  *NORD*
C
      DATA IT,N,EPS/0,20,1.E-3/
C
C START

      A=log(PE)
      PEX=PE
      CALL JON(T,PE,1,FA,RO,E,0)
C******WRITE(7,40) T,PG,PE,FA
40    FORMAT(' T,PG,PE,PGP=',4E11.4)
      IT=IT+1
      FA=log(FA/PG)
      IF(ABS(FA).LT.EPS) GOTO 101
      B=A-0.69*FA
      PEX=EXP(B)
C ONE PEMAKE ITERATION, CF. PEMAKE
      CALL JON(T,PEX,1,FB,RO,E,0)
C******WRITE(7,40) T,PG,PE,FB
      IT=IT+1
      FB=log(FB/PG)
      IF(ABS(FB).LT.EPS) GOTO 101
      X=B
C
C LOOP
      DO 100 I=1,N
      XOLD=X
C
C INTERPOLATE TO FIND NEW X
      X=A-(B-A)/(FB-FA)*FA
      PEX=EXP(X)
      IF(ABS(X-XOLD).LT.EPS) GOTO 101
      CALL JON(T,PEX,1,FX,RO,E,0)
C******WRITE(7,40) T,PG,PEX,FX
      IT=IT+1
      FX=log(FX/PG)
C
C CHECK IF A OR B CLOSEST TO X
      IF(ABS(A-X).LT.ABS(B-X)) GOTO 102
      A=X
      FA=FX
      GOTO 100
102   B=X
      FB=FX
C
C END OF LOOP
100   CONTINUE
      WRITE(7,51) N,T,PE,PG,A,B,FA,FB,EPS
51    FORMAT('0***PEMAKE, MAX ITER.: N,T,PE,PG,A,B,FA,FB,EPS=',
     * /,1X,I2,8E11.4)
      RETURN
C
C NORMAL END
101   CONTINUE
      RETURN
C
C COUNT ENTRY
      ENTRY PECNT
      WRITE(7,52) IT
52    FORMAT('0TOTAL NUMBER OF CALLS TO JON FROM PEMAKE-R =',I5)
      RETURN

      END
C
********************************************************

**********************************************************


      SUBROUTINE OSTABLOOK

      implicit real*8 (a-h,o-z)
C
       include 'parameter.inc'
C
C OS-TABle-LOOK-up is called 3 time for each itteration (from OPAC with
C J=1, called from SOLVE with T,Pe, T+dT,Pe, T,Pe+dPe). 
C It looks up the OS-tables for all molecules and calculates the
C total line opacity at each depth layer. It returns the array
C CONOS(NDP,NWL) which is the total molecular opacity at each wavenumber
C XL(NWL) and each depth of temperature T(NDP).
C The OS-tables can be generated by use of OS.F, and all tables
C are assumed to have been generated for the same wavenumbers (but
C not necessarily for the same temperatures). OSTABLOOK checks that the
C wavenumbers are identical from one table to another. If this is not the
C case, a warning is written in the output, and nearest wavenumber is used.
C It is in principle ok to use wavenumbers that differs a bit from one
C species to another, because the OS is a statistical approach, but care 
C has to be taken if the wavenumber lists are too different. It would be 
C statistically wrong to interpolate in the OS tables to the adopted OS
C wavenumbers for the specific iteration.
C This subroutine is written by UGJ in mid-1990'ties and has been 
C modified several times later.
C
C        Namelist INPUTOS contains the variables
C        NOSMOL = number of molecules considered in the opacity
C        MOLNAME(30) = the names of the molecules considered
C
C        The header of each OS table is assumed to contain a namelist,
C        INPUTOSMOL, with the followig variables:
C        NWNOS = number of OS wavenumber values for given molecule.
C        NTMOL, TMOL(mtemp) = number of temperatures and the actual 
C        temperatures for which kap is tabulated for that molecule.
C        MOLID which identifies the molecule, and ensures that the
C        right molecules are combined with the right partial pressures.
C        VKMS = the microturbulence value (only in OSs after July 1997)
C        RATIS = 13C/12C (only in OSs produced after July 1997; after
C        Aug 1998 substituted with kiso,reliso,vkms -- only for identification)
C        l_per_stellar (=0 or 1) identifies whether the absorption coefficient 
C        is in units of cm^2/mole (absorption coefficient) or cm^2/g* ("opacity").
C        Usually l_per_stellar should be 0 for molecules but 1 for atoms.
C
      DIMENSION OPJV(MTEMP),AKAPMOL(NDP)
     *      ,TMOL(MTEMP),OPLN(MTEMP),FXLN(MTEMP)
     *      ,reliso(5),OPT(MTEMP),DADT(MTEMP),WT(3,MTEMP)
      dimension xx(ndp), yy (ndp)
C      dimension fi(ndp), fxi(ndp), fyi(ndp), fxyi(ndp) !fi==akapmol
      dimension fxi(ndp), fyi(ndp), fxyi(ndp)
      CHARACTER MOLNAME*4,OSFIL*60,MOLID*4,SAMPLING*3
      LOGICAL FIRST
C
      NAMELIST /INPUTOSMOL/ MOLID, KTEMP, TMOL, NWNOS
     &   ,VKMS, KISO, RELISO, RATIS, JDERIV, L_PER_STELLAR, lchrom
C
      COMMON /COPINF/ SUMOP(maxosmol,NDP),SUMKAP(maxosmol,NDP)
      common/eostab/ xmin, ymin, dx, dy, f(mtemp,mpe)
     *       , fx(mtemp,mpe), fy(mtemp,mpe), fxy(mtemp,mpe), nx, ny
      COMMON /STATEC/PPR(NDP),PPT(NDP),PP(NDP),GG(NDP),ZZ(NDP),DD(NDP),
     & VV(NDP),FFC(NDP),PE(NDP),T(NDP),TAULN(NDP),RO(NDP),NTAU,ITER
      common /ckdtpe/dpex,kdtpe
c      COMMON /STATEC/NTAU,PPR(NDP),PPT(NDP),PP(NDP),GG(NDP),ZZ(NDP)
c     *,DD(NDP),VV(NDP),FFC(NDP),PE(NDP),T(NDP),TAULN(NDP),GEFF(NDP)
c     *,RRO(NDP),RO(NDP),PGZ(NDP)
      COMMON /ROSSC/XKAPR(NDP),CROSS(NDP)
      COMMON/CXLSET/XL(20,10),NSET,NL(10)
      COMMON /CVAAGL/XLB(500),W(500),NLB
      COMMON/COS/WNOS(NWL),CONOS(NDP,NWL),WLOS(NWL),WLSTEP(NWL)
     *    ,KOS_STEP,NWTOT,NOSMOL,NEWOSATOM,NEWOSATOMLIST
     *    ,nchrom,OSFIL(30),MOLNAME(30),SAMPLING
C      COMMON/COPPR/oppr(15,3,120,3),jvxmax,itxmax !15mol,10dpt,100wn
      COMMON/CARC3/F1P,F3P,F4P,F5P,HNIC,PRESMO(33)
      COMMON/CI4/ TMOLIM,IELEM(16),ION(16,5),MOLH,JUMP
      common/ci5/abmarcs(17,ndp),anjon(17,5),h(5),part(17,5),
     *dxi,f1,f2,f3,f4,f5,xkhm,xmh,xmy(ndp)
      COMMON /CMOLRAT/ FOLD(NDP,8),MOLOLD,KL
      COMMON /CPOLY/FACPLY,MOLTSUJI
      COMMON /TIO/ PTIOOLD(NDP),ROOLD(NDP),POXG1OLD(NDP)
      COMMON /KETIO/ PKETIO(NDP),ROKE(NDP),POXG1(NDP)
      common /fullequilibrium/ partryck(ndp,maxmol),
     &  xmettryck(ndp,maxmet),xiontryck(ndp,maxmet),partp(ndp,0:maxmol),
     &  partpp(ndp,0:maxmol)
      common /cmtest/pgm1(ndp),pgm2(ndp),pem1(ndp),pem2(ndp)
     *     ,tm1(ndp),tm2(ndp),pgos(ndp)
      common /cmasabs/ masabs(3)
      common /cosexp/ lops,nops
      common /cindiam/tdiam1,tdiam2,fdiam1,fdiam2,
     *    tc2h21,tc2h22,fc2h21,fc2h22
      common /catoms_head/ p6(9),t_at(17)
      common/catoms/opjva(nwl,9,3:17),fxln_at(nwl,9,3:17)
C opjva is ln(opac_atoms[cm2/g-*]); fxln is d(opjva)/dln(p6)
      dimension oplna(9),fxlna(9),oplna_t(15),dadta(15),wta(3,15)
     & ,oplna_tau(ndp),tat(15),p6ln(9)
      COMMON /CCIATEST/ CIATEST(44,NDP)
      COMMON /PJONINF/ P_MOL(NDP), P_NEU_HCNO(NDP), P_ION_HCNO(NDP),
     & P_NEU_HE(NDP),P_ION_HE(NDP), P_NON_HHECNO(NDP), PG_JON(NDP), 
     & HN_JON(NDP), RO_JON(NDP), P6_JON(NDP)
      COMMON /CKMOL/KMOL(MAXOSMOL)

C     parameter(nspec=892)
      common /cgem/pres_gem(ndp,nspec)
      common /cgemnames/natms_gem,nions_gem,nspec_gem,name_gem(nspec)
C atms,ions,spec ~ highest index of neutral atoms, ions, species total
      character name_gem*8
      common /cabink/abink(ndp,nspec)
      dimension ptot(ndp),pp_sum(ndp),pg(ndp)
      dimension trpe(ndp), trphe(ndp)
      common /cprespp/prespp(ndp,nspec)           !gem computed pp
      common /ctotabk/totabk(ndp,natms)
      dimension pe_gem(ndp),ptot1(ndp),dptot(ndp)
     &  ,pe1(ndp),dptot2(ndp),dpe2(ndp),dpe(ndp)

      INTEGER MOLH, JUMP
      DATA FIRST/.TRUE./

!     Partial pressures from GG-chem subroutine DEMO_SWEEP (ERC JUMP = 4)
      common /partialpressure/
     > ppN2,ppCH,ppCO,ppCN,ppC2,ppNO,ppC2H2,ppHCN,ppC2H,ppC3,ppCS,ppH2O,
     > ppOH,ppTiO,ppSiO,ppCH4,ppNH,ppSiH,ppFeH,ppVO,ppZrO,ppMgH,ppNH3,
     > ppCO2,ppTiH,ppCaH,ppCrH,ppLiH,ppH,ppO2,ppHm,ppH2,ppH2p,ppHS,
     > ppC3H,ppSiC,ppSiC2,ppNS,ppSiN,ppSO,ppS2,ppSiS,ppLaO,ppCH2,
     > ppCH3,ppSi2C,ppSiO2,ppH2S,ppCaOH,ppCHNO,ppSiF2

      real ::
     > ppN2,ppCH,ppCO,ppCN,ppC2,ppNO,ppC2H2,ppHCN,ppC2H,ppC3,ppCS,ppH2O,
     > ppOH,ppTiO,ppSiO,ppCH4,ppNH,ppSiH,ppFeH,ppVO,ppZrO,ppMgH,ppNH3,
     > ppCO2,ppTiH,ppCaH,ppCrH,ppLiH,ppH,ppO2,ppHm,ppH2,ppH2p,ppHS,
     > ppC3H,ppSiC,ppSiC2,ppNS,ppSiN,ppSO,ppS2,ppSiS,ppLaO,ppCH2,
     > ppCH3,ppSi2C,ppSiO2,ppH2S,ppCaOH,ppCHNO,ppSiF2


C
      if (first) then
      TOSREAD = 0.
      TPART = 0.
C      if (nosmol.gt.15) stop ' increase dimension for oppr '
C      if (nwtot/1000.gt.120) stop ' increase dimension for oppr '
C      if (ndp/10+1.gt.10) stop ' increase dimension for oppr '
      end if

      WRITE(6,140)
140   FORMAT(' time when starting OSTABLOOK:')
      call timex1
C      WRITE(56,140)
C      tos1 = etime(t0)
C
      do 276 ip6 = 1,9
276   p6ln(ip6) = 2.3025851*dfloat(ip6)     !ln(p6[dyn/cm2])
      DO 59 IT=1,NTAU
      DO 59 JV=1,NWTOT
         CONOS(IT,JV) = 0.
59    CONTINUE
C
C CALCULATE THE PARTIAL PRESSURE
      TBPART=SECOND()
C      TBPART=0.
      MOLM=MOLH
      MOLH=0
      div = 1./(xmy(1)*xmh)
C
C       write(6,*) ' ((abund(ik),anjon(ik,1-3)),ik=1,17) = '
C       WRITE(6,*) ' K,T(K),al4,al7,al31,ax1,ax2,ax3,ax4,ax5'
C       write(6,*) ' ((log10( max(presmo(ik),1.d-99) )),ik=1,11)'
C       write(6,*) ' ((log10( max(presmo(ik),1.d-99) )),ik=12,22)'
C       write(6,*) ' ((log10( max(presmo(ik),1.d-99) )),ik=23,33)'
C
C        if(first) then
C        write(40,*)' k,t(k),pgp,pe(k),pp(k),ppr(k),ppt(k),pgos(k)'
C        end if

        DO 2111 K=1,NTAU

!        CALL JON(T(K),PE(K),1,PGP,RO(K),EP,0)
!        ptot(k) = pe(k)+pgp
CV20    ptot(k)=pe(k)+pp(k)-ppr(k)-ppt(k)

        ptot(k)=pp(k)-ppr(k)-ppt(k)
2111    continue


       IF (JUMP.eq.3) THEN
       do 720 i=1,ntau
       pg(i)=pp(i)-ppr(i)-ppt(i)
       phei = (pg(i)-pe(i)) * totabk(i,3)/(totabk(i,2)+totabk(i,3))
       trpe(i) = 0.1d0 * pe(i)
       ptot(i) = pg(i)
       ptot(i) = 0.1d0 * ptot(i)     !1N/m^2 = 10 dynes/cm^2
       trphe(i) = (0.1d0 * phei) / ptot(i)
720   CONTINUE
782    format(i2,f8.2,1pe12.3,12x,e12.3,12x,e12.3)

C add here pieces from kemi.f:

       kexit = 0

2106  CONTINUE
C  Nulstilling:
      write(6,*) 'kdtpe = ',kdtpe

      do 3212 jk=1,ndp
      do 3212 jsp=1,nspec
3212  abink(jk,jsp) = 0.




      call gemcmsetup                            ! Make comp. matrix
      call gibbsread                             ! Read polynomials


C In Marcs opac (and hence ostablook) is called from solve with
C TT,PPE for dktpe=0
C TT+TT*dtx,PPE for dktpe=1
C TT,PPE+PPE*dpex for dktpe=2 (dpex=0.001)
      if (kdtpe.eq.0) then
        do 1520 i=1,ntau
        ptot1(i) = ptot(i)
        pg(i) = 10.*ptot1(i)
1520    continue
      else if (kdtpe.eq.1) then
        do 1526 i=1,ntau
        ptot(i) = ptot1(i)
        pg(i) = 10.*ptot1(i)
1526    continue
      else if (kdtpe.eq.2) then
        dptotx = 0.5*dpex
C       dptotx = 0.0005
C       dptotx = 0.001
        do 1521 i=1,ntau
        if(kexit.eq.0) dptot2(i) = dptotx  
        dptot(i) = dptot2(i) * ptot1(i)
        ptot(i) = ptot1(i) + dptot(i)
        pg(i) = 10.*ptot(i)

C       dpe2(i) = dpe(i)
C                         at every iterat dptot2 is the %added ptot1 which gives ptot(i)
C                         and ptot(i) gives rise to pe_gem(i) via call tstgem
C                         The aim is now to iterate so %added to pe1 is dpex (=0.001)
C                         at first iteration (kexit=0) dptot2(i)==dptotx
C                         dpe2(i) is saved from previous dpe(i) for comparison
1521    continue
      end if
C after this introduce: kdtpe=2 => loop over tstgem to find dPe/dPg


      call tstgem(t,ptot,ntau)


       DO 3151 kd=1,NTAU
       pe_gem(kd) = pg(kd) * abink(kd,1)
       pe(kd) = pe_gem(kd)
       p_particles = pg(kd)                      !pg(input) in dyn/cm^2
       pp_sum(kd) = 0.
       DO 3153 km=1,nspec
       PRESPP(kd,km) = p_particles * abink(kd,km)  !(pg)*rel.pp = pp in dyn/cm^2
       pp_sum(kd) = prespp(kd,km) + pp_sum(kd)
3153   continue
3151   continue
       write(6,2125) t(1),pe(1),pe_gem(1),pg(1),pp_sum(1)
       write(6,2126) t(ntau),pe(ntau),pe_gem(ntau),pg(ntau),pp_sum(ntau)
2126   format(32x,f8.2,1p2e14.5,2e11.3)
2125   format('T,pe,pe_gem,pg,pg_gem after gem:',f8.2,1p2e14.5,2e11.3)

C in Marcs kdtpe=0 => we use plain T,Pe as input  = first time around
      IF (kdtpe.eq.0) THEN
        do 1530 i=1,ntau
        pe1(i) = pe_gem(i)        !the Pe computed from Ptot1 in dyn/cm2
1530    continue
C So, first time (kdtpe=0) we have saved Pg,Pe in ptot1(i),pe1(i)

      ELSE IF (kdtpe.ge.2) THEN

         write(6,*)'kdtpe,kexit=',kdtpe,kexit
         kex = 0
         devmax = 1.
         devmin = 1.
         do 1533 i=1,ntau
C        if(i.le.2) write(6,2131) i,t(i),dpex,dpe(i),dptot2(i),dptotx
         dpe(i) = ( pe_gem(i)-pe1(i) )/pe1(i)
         if(dpex/dpe(i).gt.devmax) then
           devmax = dpex/dpe(i)
           imax = i
         end if
         if(dpex/dpe(i).lt.devmin) then
           devmin = dpex/dpe(i)
           imin = i
         end if
C        devmax = max(devmax,dpex/dpe(i))
C        devmin = min(devmin,dpex/dpe(i))
         if(kex.eq.0 .and.
     &    (dpex/dpe(i).gt.1.05 .or. dpex/dpe(i).lt.0.95)) then
          kex = 1
          kexit=kexit+1
         end if
1533     continue

         write(6,*)'dpex/dpe(i) min,max = ',devmin,devmax,'  --for:'
2131     format(i2,f7.0,1p6e14.5)
         write(6,2230)
     &   imin,t(imin),dpex,dpe(imin),dpex/dpe(imin),dptot2(imin)
         write(6,2230)
     &   imax,t(imax),dpex,dpe(imax),dpex/dpe(imax),dptot2(imax)
2230     format('i,T,dpex,dpe,dpex/dpe,dptot2:',i3,f7.0,1pe11.3,6e13.4)

         if(kex.eq.1 .and.kexit.le.10) then    !if in one or more depth point dPe is 5% wrong
          do 1543 i=1,ntau
          dx = dpex/dpe(i)
          if(kexit.le.1) go to 1546
C                              ensure that the counter-correction is not too big
          if(dpex/dpe2(i).lt.1.0) then
             if(dpex/dpe(i).gt.1.0) then
               dx = min(dx,0.5*(1.+dpe2(i)/dpex))
             else if(dpex/dpe(i).lt.1.0) then
               dx = max(dx,0.5*(1.+dpex/dpe2(i)))
             end if
          else if(dpex/dpe2(i).gt.1.0) then
             if(dpex/dpe(i).lt.1.0) then
               dx = max(dx,0.5*(1.+dpe2(i)/dpex))
             else if(dpex/dpe(i).gt.1.0) then
               dx = min(dx,0.5*(1.+dpex/dpe2(i)))
             end if
          end if
          
1546      continue
          dptot2(i) = dx * dptot2(i)
          dpe2(i) = dpe(i)
C         write(6,2131) i,t(i),dpex/dpe(i),dptot2(i)/dptotx
1543      continue
          go to 2106      !one more iteration
         end if
         go to 2104
      END IF

2104  CONTINUE

       END IF         !end computation of Gibbs min. for jump=3


       DO 11 K=1,NTAU
        KL=K         
C KL connects FOLD in CMOLRAT with depth in call to JON...
C PRESMO is in units of dyn/cm2; partp is in units of mol/g(stellar)
C Partryck is in units of dyn/cm2
        if (JUMP.le.3) then
          CALL JON(T(K),PE(K),1,PGP,RO(K),EP,0)

C         ptot(k)=pp(k)-ppr(k)-ppt(k)
C         ptot(k)=0.1*ptot(k)
        
C
C unit = 40 = tsuji.out
* 

!         write(6,*) 'Post JON values'
!         write(6,*) 'pe(1)= ',pe(1),' pp(1)= ',pp(1) 
!         write(6,*) 'ptot(1)= ',ptot(1)

          pgos(k) = pp(k) - ppr(k) - ppt(k)
!         write(6,*) 'pgos(',k,')= ',pgos(k)
C         if(first)
C      &  write(40,1218) k,t(k),pgp,pe(k),pp(k),ppr(k),ppt(k),pgos(k)
1218      format(i3,f8.0,1p6e11.3)
          if (JUMP.eq.1) then
            call test_Tsuji(t(k),pgp)
          else if (JUMP.eq.2) then
            call test_Tsuji_jf(t(k),pgp)
          end if

C ------------ USING GGCHEM TO COMPUTE PARTIALPRESSURES ------------  
          else if (JUMP.eq.4) then
            CALL JON(T(K),PE(K),1,PGP,RO(K),EP,0)
            call init_ggchem_ERC(t(k),ptot(k))
          end if
        
C

        ROKE(K) = RO(K)      !for transfer in KE's TiO common-block
C        epso = fold(k,4)*div*ro(k)
C        poxg1(k) = epso*1.3804e-16*t(k)   !- P(OI) dyn/cm2
C        poxg1(k) = fold(k,4)*8.24673e7*ro(k)*t(k)/xmy
        ph = 8.24673e7*ro(k)*t(k)/xmy(k)
        pox = abmarcs(5,k)/(1+anjon(5,2)/anjon(5,1))*ph
        pomol = presmo(4)+presmo(5)+presmo(7)+2.*presmo(11)+
     *          presmo(12)+presmo(27)+presmo(28)  !H2O,OH,CO,O2,NO,SiO,SO
        poxg1(k) = pox - pomol
        TRIX=T(K)*RO(K)
C
C Molecular equilibrium from Tsuji's routines (JUMP.eq.1)
C or from the old MARCS equilibrium routines ?
C ( --> also moltsuji called )
C UNITS:
C PRESMO and PARTRYCK are in units of dyn/cm2; 
C 1/R = 1.2027e-8 [mol K / erg], where R is the molar gas constant.
C 1/TRIX * 1/R = 1/(T*RO*R) is in units of [cm3 g*-1 K-1 mol K erg-1]
C PARTP is therefore in units of dyn cm-2 cm3 g*-1 mol erg-1
C But dyn = erg cm-1, so PARTP becomes in units [mol/g*], g*=g_stellar_material. 
C When PARTP is multiplied with the molecular absorption coefficient in
C units cm2/mol, the opacity comes out in units cm2/g*.
C
C For CIA the absorption coefficient is listed in units of cm-1/amagat^2,
C because it depends on the number of H2 (or H2-He) pairs, which in
C itself is pressure dependent.
C We therefore divide presmo^2 [dyn/cm2]^2 with 3.7095e3^2*T^2*ro
C in order to get partp in units cm^3*amagat^2/g*  in this case.
C XMH is g_Hydrogen per Hydrogen atom, XMY is #AU/H = g* per g_Hydrogen; 
C i.e., XMH XMY is g*/hydrogen_atom, and 1/(XMH XMY) is #H/g*
C ABUND(2) is #He_atoms/H_nuclei (normalized in INJON) so
C XNHE = ABUND(2) / (XMH XMY) is therefore #He_atoms/g*
C In the ideal gas equation P = (n/V)kT, where n/V is number
C of particles pr cm^3, k=1.38053e-16[erg/K], T[K], so P is in
C units of [cm-3*erg]=[dyn/cm^2]. 
C In our case we have "n/V" in #atoms/g*. Multiplying with ro[g*/cm3]
C brings XNHE*ro into units [#He atoms/cm^3].
C I.e., PRESMO(17)*ro = 1.38053e-16 * T(K) * XNHE * ro is in units 
C dyn/cm2 (and not presmo(17) alone).

        if (JUMP.EQ.4) then
       
C --------- PARTIAL PRESSURES FROM GGCHEM SUBROUTINE DEMO_SWEEP ------
        PARTP(K,0)  = ppH            / TRIX*1.20274D-8    ! H
        PARTP(K,1)  = ppHm           / TRIX*1.20274D-8    ! H-
        PARTP(K,3)  = ppH2p          / TRIX*1.20274D-8    ! H2+
        PARTP(K,6)  = ppCH           / TRIX*1.20274D-8    ! CH
        PARTP(K,7)  = ppCO           / TRIX*1.20274D-8    ! CO
        PARTP(K,8)  = ppCN           / TRIX*1.20274D-8    ! CN
        PARTP(K,9)  = ppC2           / TRIX*1.20274D-8    ! C2
        PARTP(K,11) = ppO2           / TRIX*1.20274D-8    ! O2
        partp(k,12) = ppNO           / trix*1.20274d-8    ! NO
        PARTP(K,14) = ppC2H2         / TRIX*1.20274D-8    ! C2H2
        PARTP(K,15) = ppHCN          / TRIX*1.20274D-8    ! HCN
        PARTP(K,16) = ppC2H          / TRIX*1.20274D-8    ! C2H
        PARTP(K,21) = ppC3           / TRIX*1.20274D-8    ! C3
        PARTP(K,22) = ppCS           / TRIX*1.20274D-8    ! CS
        PARTP(K,4)  = ppH2O          / TRIX*1.20274D-8    ! H2O
        PARTP(K,5)  = ppOH           / TRIX*1.20274D-8    ! OH 
        PARTP(K,31) = ppTiO          / TRIX*1.20274D-8    ! TiO
        PARTP(K,27) = ppSiO          / TRIX*1.20274D-8    ! SiO
        PARTP(K,39) = ppCH4          / TRIX*1.20274D-8    ! CH4

        PARTP(K,2)=(ppH2/3.7095D3/T(K) )**2 /RO(K) ! H2-H2 CIA??
           XNHE = abmarcs(2,k) / (XMH*XMY(k))
           PARTP(K,17) = 1.38053e-16 * T(K) * XNHE
        PARTP(K,17)=
     *     PARTP(K,2)*PARTP(K,17)/( 3.7095D3 * T(K) )**2  ! H2-HE CIA
           PARTP(K,17) = PARTP(K,17) * RO(K)
           PHE = PARTP(K,17)
           P6_JON(K) = xmettryck(k,1)+0.42*PHE+0.85*ppH2

        PARTP(K,13) = ppNH          / TRIX*1.20274D-8    ! NH
        PARTP(K,19) = ppSiH         / TRIX*1.20274D-8    ! SiH
        PARTP(K,216)= ppFeH         / TRIX*1.20274D-8    ! FeH
        PARTP(K,32) = ppVO          / TRIX*1.20274D-8    ! VO
        PARTP(K,33) = ppZrO         / TRIX*1.20274D-8    ! ZrO
        PARTP(K,34) = ppMgH         / TRIX*1.20274D-8    ! MgH
        PARTP(K,36) = ppCaH         / TRIX*1.20274D-8    ! CaH 

        PARTP(K,43) = ppNH3         / TRIX*1.20274D-8    ! NH3
        PARTP(K,46) = ppCO2         / TRIX*1.20274D-8    ! CO2
        partp(k,214)= ppTiH         / trix*1.20274d-8    ! TiH
        partp(k,215)= ppCaH         / trix*1.20274d-8    ! CaH
        partp(k,217)= ppCrH         / trix*1.20274d-8    ! CrH
        partp(k,71) = ppLiH         / trix*1.20274d-8    ! LiH
        partp(k,10) = ppH2          / trix*1.20274d-8    ! H2
        PARTP(K,18) = ppHS          / TRIX*1.20274D-8    ! HS
        PARTP(K,20) = ppC3H         / TRIX*1.20274D-8    ! C3H
        PARTP(K,23) = ppSiC         / TRIX*1.20274D-8    ! SiC
        PARTP(K,24) = ppSiC2        / TRIX*1.20274D-8    ! SiC2
        PARTP(K,25) = ppNS          / TRIX*1.20274D-8    ! NS
        PARTP(K,26) = ppSiN         / TRIX*1.20274D-8    ! SiN
        PARTP(K,28) = ppSO          / TRIX*1.20274D-8    ! SO
        PARTP(K,29) = ppS2          / TRIX*1.20274D-8    ! S2
        PARTP(K,30) = ppSiS         / TRIX*1.20274D-8    ! SiS
        PARTP(K,37) = ppLaO         / TRIX*1.20274D-8    ! La2
        PARTP(K,40) = ppCH2         / TRIX*1.20274D-8    ! CH2
        PARTP(K,41) = ppCH3         / TRIX*1.20274D-8    ! CH3
        PARTP(K,57) = ppSi2C        / TRIX*1.20274D-8    ! Si2C
        PARTP(K,58) = ppSiO2        / TRIX*1.20274D-8    ! SiO2
        PARTP(K,59) = ppH2S         / TRIX*1.20274D-8    ! H2S
        PARTP(K,67) = ppCaOH        / TRIX*1.20274D-8    ! CaOH
        PARTP(K,107)= ppCHNO        / TRIX*1.20274D-8    ! CHNO
        PARTP(K,135)= ppSiF2        / TRIX*1.20274D-8    ! SiF2
        PARTP(K,218)= ppN2          / TRIX*1.20274D-8    ! N2 


        CIATEST(12,k) = presmo(2)
        CIATEST(13,k) = t(k)
        CIATEST(14,k) = xnhe
        CIATEST(15,k) = presmo(17)
        CIATEST(16,k) = ro(k)
        CIATEST(17,k) = partp(k,2)
        CIATEST(18,k) = partp(k,17)
        CIATEST(19,k) = abmarcs(2,k)
        CIATEST(20,k) = xmh
        CIATEST(21,k) = xmy(k)
        CIATEST(22,k) = trix



!Defining Partial pressures in a unit suitable to print in listmo:

        PARTPP(K,0)  = ppH                                 ! H
        PARTPP(K,1)  = ppHm                                ! H-
        PARTPP(K,3)  = ppH2p                               ! H2+
        PARTPP(K,6)  = ppCH                                ! CH
        PARTPP(K,7)  = ppCO                                ! CO
        PARTPP(K,8)  = ppCN                                ! CN
        PARTPP(K,9)  = ppC2                                ! C2
        PARTPP(K,11) = ppO2                                ! O2
        PARTPP(k,12) = ppNO                                ! NO
        PARTPP(K,14) = ppC2H2                              ! C2H2
        PARTPP(K,15) = ppHCN                               ! HCN
        PARTPP(K,16) = ppC2H                               ! C2H
        PARTPP(K,21) = ppC3                                ! C3
        PARTPP(K,22) = ppCS                                ! CS
        PARTPP(K,4)  = ppH2O                               ! H2O
        PARTPP(K,5)  = ppOH                                ! OH 
        PARTPP(K,31) = ppTiO                               ! TiO
        PARTPP(K,27) = ppSiO                               ! SiO
        PARTPP(K,39) = ppCH4                               ! CH4

        PARTPP(K,2)=(ppH2/3.7095D3/T(K) )**2 /RO(K) ! H2-H2 CIA??
           XNHE = abmarcs(2,k) / (XMH*XMY(k))
           PARTPP(K,17) = 1.38053e-16 * T(K) * XNHE
        PARTPP(K,17)=
     *     PARTPP(K,2)*PARTPP(K,17)/( 3.7095D3 * T(K) )**2  ! H2-HE CIA
           PARTPP(K,17) = PARTPP(K,17) * RO(K)
           PHE = PARTPP(K,17)
           P6_JON(K) = xmettryck(k,1)+0.42*PHE+0.85*ppH2

        PARTPP(K,13) = ppNH                              ! NH
        PARTPP(K,19) = ppSiH                             ! SiH
        PARTPP(K,216)= ppFeH                             ! FeH
        PARTPP(K,32) = ppVO                              ! VO
        PARTPP(K,33) = ppZrO                             ! ZrO
        PARTPP(K,34) = ppMgH                             ! MgH
        PARTPP(K,36) = ppCaH                             ! CaH 

        PARTPP(K,43) = ppNH3                             ! NH3
        PARTPP(K,46) = ppCO2                             ! CO2
        PARTPP(k,214)= ppTiH                             ! TiH
        PARTPP(k,215)= ppCaH                             ! CaH
        PARTPP(k,217)= ppCrH                             ! CrH
        PARTPP(k,71) = ppLiH                             ! LiH
        PARTPP(k,10) = ppH2                              ! H2
        PARTPP(K,18) = ppHS                              ! HS
        PARTPP(K,20) = ppC3H                             ! C3H
        PARTPP(K,23) = ppSiC                             ! SiC
        PARTPP(K,24) = ppSiC2                            ! SiC2
        PARTPP(K,25) = ppNS                              ! NS
        PARTPP(K,26) = ppSiN                             ! SiN
        PARTPP(K,28) = ppSO                              ! SO
        PARTPP(K,29) = ppS2                              ! S2
        PARTPP(K,30) = ppSiS                             ! SiS
        PARTPP(K,37) = ppLaO                             ! La2
        PARTPP(K,40) = ppCH2                             ! CH2
        PARTPP(K,41) = ppCH3                             ! CH3
        PARTPP(K,57) = ppSi2C                            ! Si2C
        PARTPP(K,58) = ppSiO2                            ! SiO2
        PARTPP(K,59) = ppH2S                             ! H2S
        PARTPP(K,67) = ppCaOH                            ! CaOH
        PARTPP(K,107)= ppCHNO                            ! CHNO
        PARTPP(K,135)= ppSiF2                            ! SiF2
        PARTPP(K,218)= ppN2                              ! N2 


        CIATEST(12,k) = presmo(2)
        CIATEST(13,k) = t(k)
        CIATEST(14,k) = xnhe
        CIATEST(15,k) = presmo(17)
        CIATEST(16,k) = ro(k)
        CIATEST(17,k) = PARTPP(k,2)
        CIATEST(18,k) = PARTPP(k,17)
        CIATEST(19,k) = abmarcs(2,k)
        CIATEST(20,k) = xmh
        CIATEST(21,k) = xmy(k)
        CIATEST(22,k) = trix



! ------- END OF PP definition for jump = 4 ---------

        go to 2101
        
        end if


        if (JUMP.EQ.3) then               !i.e., from Gibbs minimalisation

        PARTP(K,6) = prespp(K,153) /TRIX*1.20274D-8    ! CH
        PARTP(K,7) = prespp(K,152) /TRIX*1.20274D-8    ! CO
        PARTP(K,8) = prespp(K,154) /TRIX*1.20274D-8    ! CN
        PARTP(K,9) = prespp(K,132) /TRIX*1.20274D-8    ! C2
        PARTP(K,14)= prespp(K,575) /TRIX*1.20274D-8    ! C2H2
        PARTP(K,15)= prespp(K,372) /TRIX*1.20274D-8    ! HCN  (actually CHN)
        PARTP(K,16)= prespp(K,367) /TRIX*1.20274D-8    ! C2H
        PARTP(K,21)= prespp(K,364) /TRIX*1.20274D-8    ! C3
        PARTP(K,22)= prespp(K,160) /TRIX*1.20274D-8    ! CS
        PARTP(K,4) = prespp(K,382) /TRIX*1.20274D-8    ! H2O
        PARTP(K,5) = prespp(K,155) /TRIX*1.20274D-8    ! OH 
        PARTP(K,31)= prespp(K,201) /TRIX*1.20274D-8    ! TiO
        PARTP(K,27)= prespp(K,196) /TRIX*1.20274D-8    ! SiO
        PARTP(K,39)= prespp(K,580) /TRIX*1.20274D-8    ! CH4

C       if(k.eq.1.or.k.eq.27.or.k.eq.47) write(6,784) k
C     &     ,prespp(k,152), prespp(k,153), prespp(k,154)
C     &     ,prespp(k,132), prespp(k,575), prespp(k,372)
C     &     ,prespp(k,367), prespp(k,364), prespp(k,160)
C     &     ,prespp(k,382), prespp(k,155), prespp(k,201)
C     &     ,prespp(k,196), prespp(k,580), prespp(k,128)
C784    format(i3,1p7e9.2,/5x,8e9.2)

C        if(k.eq.1) write(6,2129) abink(k,152),abink(k,372),abink(k,382)
C     &        ,abink(k,201)
C2129    format(' abink{CO,HCN,H2O,TiO} in OSTABLOOK:',1p4e12.3)

        PARTP(K,2)=(PARTRYCK(K,2)/3.7095D3/T(K) )**2 /RO(K) ! H2-H2 CIA
           XNHE = abmarcs(2,k) / (XMH*XMY(k))
           PARTRYCK(K,17) = 1.38053e-16 * T(K) * XNHE
        PARTP(K,17)=
     *     PARTRYCK(K,2)*PARTRYCK(K,17)/( 3.7095D3 * T(K) )**2  ! H2-HE CIA
           PARTRYCK(K,17) = PARTRYCK(K,17) * RO(K)
           PHE = PARTRYCK(K,17)
           P6_JON(K) = xmettryck(k,1)+0.42*PHE+0.85*partryck(k,2)

        go to 2101

        end if    !Gibbs

        if (JUMP.GE.1) then               !i.e., from Tsuji
        PARTP(K,6) =PARTRYCK(K,6)/TRIX*1.20274D-8     ! CH
        PARTP(K,7) =PARTRYCK(K,7)/TRIX*1.20274D-8     ! CO
        PARTP(K,8) =PARTRYCK(K,8)/TRIX*1.20274D-8     ! CN
        PARTP(K,9) =PARTRYCK(K,9)/TRIX*1.20274D-8     ! C2
        partp(k,12)=partryck(k,12)/trix*1.20274d-8    ! NO
        PARTP(K,14)=PARTRYCK(K,14)/TRIX*1.20274D-8    ! C2H2
        PARTP(K,15)=PARTRYCK(K,15)/TRIX*1.20274D-8    ! HCN
        PARTP(K,16)=PARTRYCK(K,16)/TRIX*1.20274D-8    ! C2H
        PARTP(K,21)=PARTRYCK(K,21)/TRIX*1.20274D-8    ! C3
        PARTP(K,22)=PARTRYCK(K,22)/TRIX*1.20274D-8    ! CS
        PARTP(K,4) =PARTRYCK(K,4)/TRIX*1.20274D-8     ! H2O
        PARTP(K,5) =PARTRYCK(K,5)/TRIX*1.20274D-8     ! OH 
        PARTP(K,31)=PARTRYCK(K,31)/TRIX*1.20274D-8    ! TiO
        PARTP(K,27)=PARTRYCK(K,27)/TRIX*1.20274D-8    ! SiO
        PARTP(K,39)=PARTRYCK(K,39)/TRIX*1.20274D-8    ! CH4
        PARTP(K,2)=(PARTRYCK(K,2)/3.7095D3/T(K) )**2 /RO(K) ! H2-H2 CIA
           XNHE = abmarcs(2,k) / (XMH*XMY(k))
           PARTRYCK(K,17) = 1.38053e-16 * T(K) * XNHE
        PARTP(K,17)=
     *     PARTRYCK(K,2)*PARTRYCK(K,17)/( 3.7095D3 * T(K) )**2  ! H2-HE CIA
           PARTRYCK(K,17) = PARTRYCK(K,17) * RO(K)
           PHE = PARTRYCK(K,17)
           P6_JON(K) = xmettryck(k,1)+0.42*PHE+0.85*partryck(k,2)
        PARTP(K,13)=PARTRYCK(K,13)/TRIX*1.20274D-8    ! NH
        PARTP(K,19)=PARTRYCK(K,19)/TRIX*1.20274D-8    ! SiH
        PARTP(K,216)=PARTRYCK(K,216)/TRIX*1.20274D-8  ! FeH
        PARTP(K,32)=PARTRYCK(K,32)/TRIX*1.20274D-8    ! VO
        PARTP(K,33)=PARTRYCK(K,33)/TRIX*1.20274D-8    ! ZrO
        PARTP(K,34)=PARTRYCK(K,34)/TRIX*1.20274D-8    ! MgH
        PARTP(K,36)=PARTRYCK(K,36)/TRIX*1.20274D-8
        PARTP(K,43)=PARTRYCK(K,43)/TRIX*1.20274D-8    ! NH3
        PARTP(K,46)=PARTRYCK(K,46)/TRIX*1.20274D-8    ! CO2
        partp(k,214)=partryck(k,214)/trix*1.20274d-8  ! TiH
        partp(k,215)=partryck(k,215)/trix*1.20274d-8  ! CaH
        partp(k,217)=partryck(k,217)/trix*1.20274d-8  ! CrH
        partp(k,71) =partryck(k,71)/trix*1.20274d-8   ! LiH
        partp(k,10) = partryck(k,2)/trix*1.20274d-8   ! H2

        else

        PARTP(K,6)=PRESMO(6)/TRIX*1.20274D-8      ! CH
        PARTP(K,7)=PRESMO(7)/TRIX*1.20274D-8      ! CO
        PARTP(K,8)=PRESMO(8)/TRIX*1.20274D-8      ! CN
        PARTP(K,9)=PRESMO(9)/TRIX*1.20274D-8      ! C2
        PARTP(K,14)=PRESMO(14)/TRIX*1.20274D-8    ! C2H2
        PARTP(K,15)=PRESMO(15)/TRIX*1.20274D-8    ! HCN
        PARTP(K,16)=PRESMO(16)/TRIX*1.20274D-8    ! C2H
        PARTP(K,21)=PRESMO(21)/TRIX*1.20274D-8    ! C3
        PARTP(K,22)=PRESMO(22)/TRIX*1.20274D-8    ! CS
        PARTP(K,4)=PRESMO(4)/TRIX*1.20274D-8      ! H2O
        PARTP(K,5)=PRESMO(5)/TRIX*1.20274D-8      ! OH
        PARTP(K,31)=PRESMO(31)/TRIX*1.20274D-8    ! TiO
        PARTP(K,27)=PRESMO(27)/TRIX*1.20274D-8    ! SiO
        PARTP(K,2)=(PRESMO(2)/3.7095D3/T(K) )**2 /RO(K) ! H2-H2 CIA
        XNHE = abmarcs(2,k) / (XMH*XMY(k)) 
        PRESMO(17) = 1.38053e-16 * T(K) * XNHE
        PARTP(K,17)=PRESMO(2)*PRESMO(17)/( 3.7095D3 * T(K) )**2  ! H2-HE CIA
        PRESMO(17) = PRESMO(17) * RO(K)
        PARTP(K,13)=PRESMO(13)/TRIX*1.20274D-8    ! NH
        PARTP(K,19)=PRESMO(19)/TRIX*1.20274D-8    ! SiH
C       PARTP(K,x)=PRESM(K,x)/TRIX*1.20274D-8      ! FeH

        CIATEST(12,k) = presmo(2)
        CIATEST(13,k) = t(k)
        CIATEST(14,k) = xnhe
        CIATEST(15,k) = presmo(17)
        CIATEST(16,k) = ro(k)
        CIATEST(17,k) = partp(k,2)
        CIATEST(18,k) = partp(k,17)
        CIATEST(19,k) = abmarcs(2,k)
        CIATEST(20,k) = xmh
        CIATEST(21,k) = xmy(k)
        CIATEST(22,k) = trix

        end if


2101  CONTINUE

      open(unit=81,file='MOL2.dat',status='replace')
      write(81,'(4A8)') 'k    ', 'CO    ', 'H2O    ','TiO    '
      DO 10 NM=1,k
        write(81,'(I3,3f8.3)') NM,LOG10(partp(NM,7)),LOG10(partp(NM,4))
     &  , LOG10(partp(NM,31))
10    CONTINUE
      close(81)

C temporerely the OSs for CO, CN, and C2 exists in cm^2/g as well
C as in cm^2/mol. Above, when partp was computed from presmo, it was 
C assumed that the OS opacity listed was in units of cm^2/mol for 
C all molecules.
C OSs (ODFs) for CO, CN, C2 was originally introduced by Querci, and used 
C in the original version of the MARCS code in units of cm^2/g.
C If old OSs (i.e., masabs=1) are used convert them here from cm^2/g to 
C cm^2/mol in a 
C consistent manner with the OSGEP computation to avoid 12C/13C errors.
C In future only cm^2/mol should
C be used, because it doesn't introduce an error in the conversion
C when 12C/13C is varied as do in principle the cm^2/g OSs.

       if (masabs(1).eq.1)  PARTP(K,7)=PARTP(K,7)*28.1856  ! CO
       if (masabs(2).eq.1)  PARTP(K,8)=PARTP(K,8)*26.1902  ! CN
       if (masabs(3).eq.1)  PARTP(K,9)=PARTP(K,9)*24.1116  ! C2
         
   11 CONTINUE
      if(first)write(6,*)' after call to JON loop (11) in OSTABLOOK'
C
C Now FO = FOLD(k,4) = P(OI)/P(H) is known for all depths
C We compute the partial pressure PKETIO(k) in dyn/cm2 from
C KE's subroutine, which assumes that the formation of TiO
C do not change the oxygen available for other molecules.
C
C     CALL TIOEQ(NTAU,T,PE)
C
      MOLH=MOLM
      TEPART=SECOND()
C      TEPART=0.
      TPART = TPART + TEPART-TBPART
C
C CORRECT FOR THE PARTIAL PRESSURE
C
      TBREAD=SECOND()
C      TBREAD=0.
C

C Now consider each molecule in turn in loop 100:


       write(56,*)
       write(56,1502)iter,ntau
1502   format('From ostablook for iteration',i3
     &   ,' in depth 1,15,30,ntau (=',i2,'):')
       write(56,1503) t(1),t(15),t(30),t(ntau)
       pg1 = pp(1)-ppr(1)-ppt(1)
       pg2 = pp(15)-ppr(15)-ppt(15)
       pg3 = pp(30)-ppr(30)-ppt(30)
       pg4 = pp(ntau)-ppr(ntau)-ppt(ntau)
C      write(56,1504) pg(1),pg(15),pg(30),pg(ntau)
       write(56,1504) pg1,pg2,pg3,pg4
       write(56,1505) pe(1),pe(15),pe(30),pe(ntau)
1503   format('T(k):',4f11.3)
1504   format('Pg(k):',1p4e13.4)
1505   format('Pe(k):',1p4e13.4)


      DO 100 NM=1,NOSMOL
C
      DO 102 IT=1,NTAU
      SUMOP(NM,IT) = 0.
      SUMKAP(NM,IT) = 0.
102   CONTINUE

      if (molname(nm).eq.'ATOM') then
         molid = 'ATOM'
         IDMOL = 1
         go to 215
      end if
      OPEN(UNIT=3,FILE=OSFIL(NM),STATUS='OLD',readonly)
      READ(3,INPUTOSMOL)
      IF (MOLID.EQ.'AJ1 ') THEN
         MOLID = 'C2  '
      ELSE IF (MOLID.EQ.'AJ2 ') THEN
         MOLID = 'CN  '
      ELSE IF (MOLID.EQ.'AJ3 ') THEN
         MOLID = 'HCN '
      ELSE IF (MOLID.EQ.'AJ4 ') THEN
         MOLID = 'C2H2'
      ELSE IF (MOLID.EQ.'AJ5 ') THEN
         MOLID = 'H2O '
      ELSE IF (MOLID.EQ.'AJ6 ') THEN
         MOLID = 'TiO '
      ELSE IF (MOLID.EQ.'AJ7 ') THEN
         MOLID = 'CH  '
      ELSE IF (MOLID.EQ.'H2OS' .or. MOLID.EQ.'swmx') THEN
         MOLID = 'H2O '
      END IF
C     IF (FIRST) WRITE(6,INPUTOSMOL)
      IF (MOLNAME(NM).NE.MOLID) THEN
         print 570, NM,MOLNAME(NM),MOLID
570      FORMAT(' MOLECULE IN INPUT,FILE = ',I3,2X,A4,2X,A4)
         STOP ' MAYBE THE WRONG OPACITY FILE HAS BEEN OPENED ? '
      END IF
      NWMOL=0
      ISKIP=0
      IF (MOLID.EQ.'CH  ') THEN
          IDMOL = 6
      ELSE IF (MOLID.EQ.'CO  ') THEN
          IDMOL = 7
      ELSE IF (MOLID.EQ.'CN  ') THEN
          IDMOL = 8
      ELSE IF (MOLID.EQ.'C2  ') THEN
          IDMOL = 9
      ELSE IF (MOLID.EQ.'C2H2') THEN
          IDMOL = 14
      ELSE IF (MOLID.EQ.'DIAM') THEN
          IDMOL = 14                  !diamond grains
      ELSE IF (MOLID.EQ.'AMPC') THEN
          IDMOL = 14                  !amorphous carbon grains
      ELSE IF (MOLID.EQ.'HCN ') THEN
          IDMOL = 15
      ELSE IF (MOLID.EQ.'C2H ') THEN
          IDMOL = 16
      ELSE IF (MOLID.EQ.'C3  ') THEN
          IDMOL = 21
      ELSE IF (MOLID.EQ.'CS  ') THEN
          IDMOL = 22
      ELSE IF (MOLID.EQ.'CH4 ') THEN
          IDMOL = 39
      ELSE IF (MOLID.EQ.'H2O ') THEN
          IDMOL = 4
      ELSE IF (MOLID.EQ.'OH  ') THEN
          IDMOL = 5
      ELSE IF ((MOLID.EQ.'TIO ').OR.(MOLID.EQ.'TiO ')) THEN
          IDMOL = 31
      ELSE IF ((MOLID.EQ.'SIO ').OR.(MOLID.EQ.'SiO ')) THEN
          IDMOL = 27
      ELSE IF (MOLID.EQ.'H2H2') THEN
          IDMOL = 2
      ELSE IF (MOLID.EQ.'H2HE'.OR.(MOLID.EQ.'H2He')) THEN
          IDMOL = 17
      else if (molid.eq.'NO  ') then
          idmold = 12
      ELSE IF (MOLID.EQ.'NH  ') THEN
          IDMOL = 13
      ELSE IF (MOLID.EQ.'SIH '.OR.(MOLID.EQ.'SiH ')) THEN
          IDMOL = 19
      ELSE IF (MOLID.EQ.'VO  ') THEN
          IDMOL = 32
      ELSE IF (MOLID.EQ.'ZRO '.OR.(MOLID.EQ.'ZrO ')) THEN
          IDMOL = 33
      ELSE IF (MOLID.EQ.'MGH '.OR.(MOLID.EQ.'MgH ')) THEN
          IDMOL = 34
      ELSE IF (MOLID.EQ.'CAH '.OR.(MOLID.EQ.'CaH ')) THEN
          IDMOL = 215
      else if (molid.eq.'TIH '.or.(molid.eq.'TiH ')) then
          idmol = 214
      else if (molid.eq.'FEH ' .or. (molid.eq.'FeH ')) then
          idmol = 216
      else if (molid.eq.'CRH ' .or. (molid.eq.'CrH ')) then
          idmol = 217
      else if (molid.eq.'LIH ' .or. (molid.eq.'LiH ')) then
          idmol = 71
      else if (molid.eq.'H2  ') then
          idmol = 10
      ELSE IF (MOLID.EQ.'NH3 ') THEN
          IDMOL = 43
      ELSE IF (MOLID.EQ.'CO2 ') THEN
          IDMOL = 46
      ELSE IF (MOLID.EQ.'at1d') THEN    !atomic abs.coef. as function of temp.only
          IDMOL = 1
      ELSE
          STOP ' ********** molecule not found ********* '
      END IF

215   continue

C connect the opacity index with the partial pressure index
      KMOL(NM) = IDMOL
      if(first) write(6,620) nm,molid,kmol(nm)
620   format (' nm=',i3,' for molecule: ',a4,' with presmo index',i3)
C     
C in case of metal lines, use mihala's two dimensional interpolation 
C routine to find kappa(metals) for given T,Pe
C      IF (MOLID.EQ.'ATOM') THEN                   
C ~~~~~~~~  atoms
C
C
CC              opjvit(ip6) = 2.3025851*opl(ip6,itemp) !log_10(listed) -> ln(for intp)
CC280           opjva(jv,ip6,itemp) = opjvit(ip6)
CCC  The derivative, DADP(ip6)=d(ln(opac))/d(p6ln), at each p6 and temp.
CC              CALL TB04AT(np6,np6,p6ln,opjvit,DADP,WP)
CC              do 282 ip6=1,np6
CC282           fxln(jv,ip6,itemp) = dadp(ip6)   !d(ln(opac))/d(ln(p6)) for t=temp
CC
C
C      xmin = LOG (TMOL(1))
C      nx = KTEMP
C      ymin = 1.
C      ny = 1
C      dx = LOG(TMOL(2))-LOG(TMOL(1))
C      dy = 1.
C      DO 106 IT=1,NTAU
C      xx(it) = LOG( T(IT) )
C      yy(it) = 1.
C106   CONTINUE   
C      END IF
C
c     write(6,*) ' nwtot = ',nwtot
      jvs = nwtot
      kon0wr = 0
      kread = 0
      kosr = 1


C we have now identified the molecule and starts reading the 
C abs.koef. data in order to compute the total opacity in 
C NWTOT OS opacity frequency points

C
C If we compute a photospheric-chromospheric spectrum (from an
C illuminated atmosphere), find the temperature minimum
          itmax2= 0
          itmin = 1
          if(nchrom.eq.1) then
             do 1612 it = 2,ntau
             if(t(it) .lt. t(it-1)) itmin = it
1612         continue
C            tmin = t(itmin)
C if there is a secondary maximum, i.e. a cooler region above the chromosphere,
C then for now ignore that in the spectrum (this is for strongly illuminated
C atmospheres). The maximum temperature (i.e. the top of the chromosphere)
C is at itmax2 and we compute the chromospheric atomic opacity above itmax2 
C with the corresponding gas pressure from the crhomospheric at the same 
C temperature below itmax2.
             do 1622 it = 2,itmin
             if(t(it) .gt. t(it-1)) itmax2 = it
1622         continue
          end if
C


      DO 110 JV=1,NWTOT

      JVN = MIN(NWTOT,JV+1)
      JVPR = MAX(1,JV-1)
      DW = ( WNOS(JVN)-WNOS(JVPR) )/2.

      IF (MOLID.EQ.'ATOM') go to 315
      JVLAST = JV-1
      LEPS = 0
      IF (ISKIP.EQ.1) GO TO 111
116   CONTINUE
      do 1161 kos = 1,kosr
      READ(3,*,END=199) WN,(OPJV(IT),IT=1,KTEMP)
      if (jderiv.eq.1) READ(3,*) (FX(IT,1),IT=1,KTEMP)
1161  continue
      kread = kread + 1
111   CONTINUE

C Generally the abs.coef. data will have computed to be used
C for the OS frequencies we have also chosen for the model
C atmosphere computations. If, however,
C some abs.coef. data are NOT computed for the OS frequencies
C we use for this model atmosphere computation, then make sure
C for each OS frequency to take a nearby abs.coef. OPJV set.
C The following lines of computing should ensure this when there 
C are too many as well as too few abs.coef. values compared to
C the number of OS frequency values.

      ISKIP=0
      EPSILON = (WNOS(JVN)-WNOS(JV))/100.
      IF (kread.eq.1) then
         if(WN .GT. WNOS(JV)+EPSILON ) THEN
          ISKIP=1
          GO TO 110 !i.e. OS freq.JV < first wn
         end if
      END IF
      IF (WN .LT. WNOS(1)-EPSILON) THEN
         JVFIRST = JV + 1
         GO TO 116 !i.e. wn < first OS freq.
      END IF
      IF (WN .LT. WNOS(JV)-EPSILON) THEN
         LEPS = LEPS + 1
         GO TO 116 !i.e. wn < present OS freq.
      END IF 
      IF (WN .GT. WNOS(JV)+EPSILON .and. LEPS.LE.1) THEN
          ISKIP=1
      END IF
C
C ..now we have found the first OSWN on mol.list correspond to Marcs-OS-wn,
C hereafter use each kos_step of the os-wn for the radiative transfer:
      kosr = kos_step

         IF (jderiv .EQ. 1) go to 315  !0~derivatives are not in OS file

C     READ(3,*,END=199) WN,(OPJV(IT),IT=1,KTEMP)
C     if (jderiv.eq.1) READ(3,*) (FX(IT,1),IT=1,KTEMP)


         DO 312 IT=1,KTEMP
         IF (OPJV(IT).EQ.0) THEN
           DO 314 ITN=1,KTEMP
           FX(ITN,1)=0.
314        CONTINUE
           GO TO 315
         END IF
312      OPT(IT)=LOG( OPJV(IT) )

C  Calculate the derivative, FX/OPJV = DADT(IT)=d(ln(opt)/d(tmol), 
C  at each temperature.

         CALL TB04AT(MTEMP,KTEMP,TMOL,OPT,DADT,WT)

         DO 313 IT=1,KTEMP
313      FX(IT,1) = DADT(IT)*OPJV(IT)    !fx = d(opt)/d(tmol)

315      CONTINUE

C
      NWMOL=NWMOL+1
      SUMDWN=SUMDWN+ABS(WN-WNOS(JV))
C
      IF (MOLID.EQ.'ATOM') THEN                   
C ~~~~~~~~  atoms
C      DO 112 IT=1,KTEMP
C112   f (it,1) = LOG ( OPJV(IT) )
C      call interp( ntau, xx, yy, AKAPMOL, fxi, fyi, fxyi )
C
C ..begin  do loop 120 for evaluating kappa in each depth layer
      II=-1
C     common/catoms/opjva(nwl,9,3:17),fxln_at(nwl,9,3:17)
      MTA = 15
      KTA = 15
      MPA = 9
      KPA = 9

      DO 1210 IDP=1,NTAU
      TINT=T(IDP)
C     PINT = PGOS(IDP)
      PINT = P6_JON(IDP)

      if (pint.le.0.0) then
      write(6,*) ' P6 = 0 ????? '
       WRITE(6,*) ' I,P_MOL(I), P_NEU_HCNO(I), P_ION_HCNO(I)'
     & ,' P_NEU_HE(I), P_ION_HE(I), P_NON_HHECNO(I), P6_JON(I)'
     & ,' PG_JON(I), HN_JON(I), RO_JON(I)'
      DO 2064 I=1,NTAU
       WRITE(6,2062) I,P_MOL(I), P_NEU_HCNO(I), P_ION_HCNO(I)
     & ,P_NEU_HE(I), P_ION_HE(I), P_NON_HHECNO(I), P6_JON(I)
     & ,PG_JON(I), HN_JON(I), RO_JON(I)
2064  continue
 2062 FORMAT(i5,1P15e11.3)
      end if

      pintln = log(pint)
C      if(first) write(6,1217) idp,tint,pint,jv,wnos(jv)
1217  format(' Atoms: i,t,p=',i3,f9.0,1pe12.3/,0p,' jv,wnos:',i3,f10.2)
      do 1211 it = 1,15
      tat(it) = t_at(it+2)
      do 1212 ip = 1,9
      oplna(ip) = opjva(jv,ip,it+2)
      fxlna(ip) = fxln_at(jv,ip,it+2)
1212  continue
C      if(first)write(6,1219) tat(it),(oplna(ip),ip=1,9)
1219  format(' t,op(t,p6),p6=1,..1.e9:',f9.0,1p9e10.3)
      oplna_t(it) = TG01B(II,MPA,KPA,p6ln,OPLNA,FXLNA,pintln)
C      if(first)write(6,1220) it,pint,oplna_t(it)
1220  format(' it,pg_int,op(pg,t)=',i3,1pe11.3,e12.3)
1211  continue
      CALL TB04AT(MTA,KTA,tat,oplna_t,dadta,wta)
      oplna_tau(idp) = TG01B(II,MTA,KTA,TAT,OPLNA_T,DADTA,TINT)
      akapmol(idp) = exp ( oplna_tau(idp) )
C      if(first)write(6,1225) idp,tint,pint,akapmol(idp)
1225  format(i3,f9.0,1p2e10.3)

      SUMOP(NM,IDP) = SUMOP(NM,IDP) + akapmol(idp)*DW

C absko in cm/g*

1210  CONTINUE

C after loop 1210 akapmol(idp), idp=1,ntau is atomic abs.coef. in cm^2/g*
C      if (first)write(6,1216) (sumop(idp),idp=1,ntau,10)
C1216  format(' atomic sumop(tau,10) = ',1p8e12.3)

      ELSE                                  
C ~~~~~~~~ molecules (or one-dimensional atomic abs.coef. cm^2/g*(t))
C
      NULOP = 0
      II=-1
      DO 123 IT=1,KTEMP
      IF (OPJV(IT).EQ.0.) THEN
       OPJV(IT) = 1.E-25
       NULOP = NULOP + 1
      END IF
123   CONTINUE
C
      IF(NULOP.EQ.KTEMP) GO TO 110  ! no opacity for this OS waven.
C
      ZSUM = 0.
      DO 122 IT=1,KTEMP
C     IF (FACPLY.NE.1.) OPJV(IT) = FACPLY*OPJV(IT)
      OPLN(IT)=LOG(OPJV(IT))
      FXLN(IT)=FX(IT,1)/OPJV(IT)
      ZSUM = ZSUM + FXLN(IT)
122   continue
      t1 = tmol(1)
      t2 = tmol(2)
      op1 = opjv(1)
      op2 = opjv(2)
      iop2 = 2
C

C ..begin  do loop 120 for evaluating kappa in each depth layer
C for the 1-dimensional atomic OS we need to split it in two if there
C is a chromosphere caused by an illuminated atmosphere, and interpolate 
C in the chromospheric OS for temperatures above the temperature 
C minimum, and in the photospheric OS for layers below the temeperature 
C minimum.

      DO 120 IT=1,NTAU
      TINT=T(IT)

C if opac is == 0 at some temperatures spline may give
C negative opacities at some depth layers. To avoid this the
C derivatives are set to zero at all temperatures, and we use 
C linear interpolation in opact here:
      IF (ZSUM.EQ.0.) THEN
       do 1221 iop = iop2,ktemp
       if (tmol(iop).gt.tint) then
         t1 = tmol(iop-1)
         t2 = tmol(iop)
         op1 = opjv(iop-1)
         op2 = opjv(iop)
         iop2 = iop
         go to 1222
       end if
1221   continue
1222   continue
      AKAPMOL(IT) = (tint-t1)/(t2-t1)*(op2-op1) + op1
C .. but normally the derivative for spline interpol in temp is given:
       ELSE
      AKAPMOL(IT) = TG01B(II,MTEMP,KTEMP,TMOL,OPLN,FXLN,TINT)
      AKAPMOL(IT) = EXP ( AKAPMOL(IT) )
       END If

C Her er (27/8-96) insat en lille rutine der for T < Tdiam indsaetter
C absorptions koefficienten for X gange diamanterne og Y gange C2H2.

      IF (MOLID.EQ.'DIAM' .or. MOLID.EQ.'AMPH') THEN
C                              for t < tdiam1 fakdiam = fdiam1 (from input)
      if (t(it).le.tdiam1) then
        fakdiam = fdiam1
      else if (t(it).ge.tdiam2) then
C                              for t > tdiam2 fakdiam = fdiam2 (from input)
        fakdiam = fdiam2
      else
C                               i.e., lin.intp. for tdiam1 < t(it) < tdiam2
        fakdiam = (t(it)-tdiam1)/(tdiam2-tdiam1)*(fdiam2-fdiam1)+fdiam1
      end if
           AKAPMOL(IT) = fakdiam * AKAPMOL(IT)
C       if (jv/100*100.eq.jv)
C    &  write(6,*) ' diam: jv,it,t,fakt=',jv,it,tint,fakdiam
      END IF

      IF (MOLID.EQ.'C2H2' .and. 
     &           (fc2h21.ne.1..or.fc2h22.ne.1.)) THEN
C                             for t < tc2h21 fakc2h2 = fc2h21 (from input)
      if (t(it).le.tc2h21) then  
        fakc2h2 = fc2h21
      else if (t(it).ge.tc2h22) then
C                              for t > tc2h22 fakc2h2 = fc2h22 (from input)
        fakc2h2 = fc2h22         
      else
C                               i.e., lin.intp. for tc2h21 < t(it) < tc2h22
        fakc2h2 = (t(it)-tc2h21)/(tc2h22-tc2h21)*(fc2h22-fc2h21)+fc2h21
      end if
            AKAPMOL(IT) = fakc2h2 * AKAPMOL(IT)
C       if (jv/100*100.eq.jv)
C    &  write(6,*) ' c2h2: jv,it,t,fakt=',jv,it,tint,fakc2h2
      END IF

      SUMOP(NM,IT) = SUMOP(NM,IT) + AKAPMOL(IT)*DW

C akapmol(k) in cm^2/mol; sumop in cm/mol integrated over all OS range
120   CONTINUE

C ..end do loop 120 for evaluating kappa in each depth layer

1202  format(i3,i3,3f6.0,1p3e9.2)
C
      END IF
C ~~~~~~~~ end if for metals/molecules
C
C NB: we shift to index for wavelength compared to the listing
C     in the OS-tables   
C this is because the planck function and the calls to opac etc
C are for (increasing) wavelengths. Also the indexing of WLOS is
C reversed compared to the indexing of WNOS.

      JVL=NWTOT-JV+1
C
C experiment with shifting the opacities lops os-steps compared to the
C real computation (jvl = jvl + lops),
C and using only each nops step of the opacity:     
      jvl = jvl + lops
      if (jvl.gt.nwtot) go to 110  
C
      DO 108 IT=1,NTAU
        if (molid.eq.'ATOM' .or. molid.eq.'at1d') then
           f1 = -1.
        else
        F1=LOG10(max(PARTP(IT,IDMOL),1.d-99))
        end if
        F2=LOG10(max(AKAPMOL(IT),1.d-99))
      IF (F1+F2.LE.-33) THEN
         kon0wr = kon0wr + 1
         if (first) then
          if (kon0wr.le.1) then
           WRITE(6,6633) idmol,it,t(it),partp(it,idmol),akapmol(it)
          else if (it/10*10.eq.it .and.jvl/500*500.eq.jvl) then
           WRITE(6,6634) kon0wr,idmol,it,partp(it,idmol),akapmol(it)
          end if
         end if
6633     format(' CONS=0; idmol,it,t=',2i4,f7.0
     #          ,' because partp,akapmol=',1p2e10.2)
6634     format(' CONS=/0; con0wr,idmol,it=',i7,2i4
     #          ,'; partp,akapmol=',1p2e12.3)
      ELSE IF (F1+F2.GE.33) THEN
         PRINT*,' OVERFLOW FOR JV,IT,IDMOL = ',JV,IT,IDMOL
         STOP ' OVERFLOW IN OSTABLOOK '
      ELSE
         if (jvl/nops*nops.eq.jvl)jvs = jvl
        if (molid.eq.'ATOM' .or. molid.eq.'at1d') then
         if (molid.eq.'at1d'.and.nchrom.eq.1 ) then
           if (lchrom.eq.1.and.it.ge.itmin) akapmol(it) = 0.0
           if (lchrom.eq.0.and.it.lt.itmin) akapmol(it) = 0.0
         end if
         CONOS(IT,JVL)=AKAPMOL(IT)+CONOS(IT,JVS)
         sumkap(nm,it) = sumkap(nm,it) + AKAPMOL(IT) *dw
        else
         CONOS(IT,JVL)=PARTP(IT,IDMOL)*AKAPMOL(IT)+CONOS(IT,JVS)
         sumkap(nm,it) = 
     &         sumkap(nm,it) + AKAPMOL(IT)*partp(it,idmol)* dw
        end if
C         CONOS(IT,JVL)=PARTP(IT,IDMOL)*AKAPMOL(IT)+CONOS(IT,JVL)
      END IF

C conos is in cm^2/g*, so sumkap(imol) is in units cm/g* 

C      itxmax = 3
C       if (it.eq.7) then 
C           itx = 1
C       else if (it.eq.27) then
C           itx = 2
C       else if (it.eq.40) then
C           itx = 3
C       else 
C           itx = 4
C       end if
C
C      if (jv/1000*1000.eq.jv .and. (itx.le.3) ) then
C         jvx = max(1,jv/1000)
C         oppr(nm,itx,jvx,1) = PARTP(IT,IDMOL)*AKAPMOL(IT)
C         oppr(nm,itx,jvx,2) = PARTP(IT,IDMOL)
C         oppr(nm,itx,jvx,3) = AKAPMOL(IT)
C         jvxmax = jvx
C      end if

108   CONTINUE
C
110   CONTINUE                    ! next wavenumber, wnos(jv)

199   CONTINUE                    ! no more data for this molecule

      IF (first) THEN
      if (kon0wr.gt.0)
     # WRITE(6,6634) kon0wr,idmol,it,partp(it,idmol),akapmol(it)
      write(6,1404) NWMOL,NWTOT
1404  format(' We assigned opacity data for',
     &   i6,' OS frequency points, out of',i6,' OS freq. in total')
      if (nwmol.ne. jvlast-jvfirst+1)
     & write(6,*) ' jvlast,jvfirst =', jvlast,jvfirst

      SUMDWN=SUMDWN/NWMOL
      IF (SUMDWN.NE.0. .and. molid.ne.'ATOM') 
     &    WRITE(6,1401) NWMOL,MOLID,SUMDWN
1401  FORMAT(' SSD (sum wn-wnmol/',I6,') for ',A4,' was:',1PE12.4)

        write(6,1402) MOLID
        write(6,1403) (SUMOP(NM,K),K=1,NTAU)
        write(6,14021) MOLID
        write(6,1403) (SUMKAP(NM,K),K=1,NTAU)
      if (molid.eq.'ATOM') then
         write(6,*) ' For atoms: (pe(i),pgos(i),p6_jon(i),i=1,ntau):'
         write(6,1216) (pe(idp),pgos(idp),p6_jon(idp),idp=1,ntau)
         WRITE(6,*) ' I,P_MOL(I), P_NEU_HCNO(I), P_ION_HCNO(I)'
     &   ,' P_NEU_HE(I), P_ION_HE(I), P_NON_HHECNO(I), P6_JON(I)'
     &   ,' PG_JON(I), HN_JON(I), RO_JON(I)'
        DO 3064 I=1,NTAU
         WRITE(6,3062) I,P_MOL(I), P_NEU_HCNO(I), P_ION_HCNO(I)
     &   ,P_NEU_HE(I), P_ION_HE(I), P_NON_HHECNO(I), P6_JON(I)
     &   ,PG_JON(I), HN_JON(I), RO_JON(I)
3064    continue
3062    FORMAT(i5,1P15e11.3)
      end if
      END IF                 !end of write out for first call to ostablook
1216  format(1p3e12.3)

C for control of the run and effect of adding dPe and dT, write in unit=56:
       kmw = kmol(nm)
C      write(56,1506) molid,partp(1,idmol),partp(15,idmol)
       write(56,1506) molname(nm),partp(1,idmol),partp(15,idmol)
     &                ,partp(30,idmol),partp(ntau,idmol)
       write(56,1507) sumop(nm,1),sumop(nm,15)
     &                ,sumop(nm,30),sumop(nm,ntau)
       write(56,1508) sumkap(nm,1),sumkap(nm,15)
     &                ,sumkap(nm,30),sumkap(nm,ntau)
1506   format('part.pres.[mol/g*] of ',a4,':',1p4e12.4)
1507   format('abs.coef.[cm/mol]: ',8x,1p4e12.4)
1508   format('opacity[cm/g*]: ',11x,1p4e12.4)


      CLOSE(3)
100   CONTINUE                                    
C next molecule   (/or all atoms)

       do 1509 nmw=1,nosmol
       kmw = kmol(nmw)
       sumpp1 = sumpp1 + partp(1,kmw)
       sumpp2 = sumpp2 + partp(15,kmw)
       sumpp3 = sumpp3 + partp(30,kmw)
       sumpp4 = sumpp4 + partp(ntau,kmw)
1509   continue
       write(56,1510) nosmol,sumpp1,sumpp2,sumpp3,sumpp4
       write(56,1511) sumpp1/pg1,sumpp2/pg2,sumpp3/pg3,sumpp4/pg4
1510   format('sum part.pres. of ',i2,' OS molec:',1p4e10.3)
1511   format('ditto relative Pg:',1p4e10.3)
C      tos2 = etime(t0)
C      tos = tos2  - tos1
C      WRITE(56,1540) tos
C1540   FORMAT(' time processing opacities:',f8.1,'sec.')

      TEREAD=SECOND()
      TOSREAD = TOSREAD + TEREAD-TBREAD
C
1402   FORMAT(' INTEGRATED ABS.COF. IN CM/MOL FOR ',A4)
14021  FORMAT(' INTEGRATED OPACITY IN CM/g* FOR ',A4)
1403   FORMAT(1P8E10.3)
C
      FIRST = .FALSE.
C      CALL TIME1
      WRITE(6,141)
141   FORMAT(' total and accumulated time processing OS tables:')
      call timex1

      return
      END
C
C
