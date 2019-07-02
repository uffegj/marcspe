!         version 2018      by ERC: coupling to GGCHEM to model down to 100 K
!         version 2015      by DJ: coupling to DRIFT to model clouds
!         version 00-03-22  by UGJ includes JF's equilibrium and C2H opac.
!         version 97-03-14  by Uffe Graae Jorgensen: double precision
!         version 94-12-15  by Christiane Helling: Tsuji's molecular eq.
!  VMS    version 93-01-01  by Uffe Graae Jorgensen (Opacity Sampling, Sphericity)
!         version 90-03-07  by Aake Nordlund (stripped all line code, ODF)
!  Cyber  version 1975      by Gustafsson et al. (ODF)
!
! Iteration control parameters:
! NEWMOD can take the following values:
!     000=> start from a stored model
!     001=> start a new model
!     002=> no iterations; transfer in old model
!     003=> continue with new teff,logg,a/h
!     004=> continue with new tau-scale
!     005=> combine old model with temp from other file 
!     006=> start from a binary stored model
!     009=> one iteration to compute the limb functions, no correction on
!           the model, be careful to use the same parameters in the input file.
! If NEWMOD=1, one reads the last 3 entries of the input file, namely: 
!     taucnv=> approx. tau where convection starts
!     dtblnk=> approx. backwarming and surface cooling
!     taubln=> turnover between backw. and surf. cool.
!
! JONTYP can be:
!     001=> constant partition functions
!     002=> full pf, Bashek et al.
!     003=> full pf, Fischel et al., all possible ionisation states
!  
! JUMP can be:
!     1-3 => subroutine MOL doesn't work => no presmo
!     1 => Tsuji-routines, partryck for molecules, xmettryck for atoms
!     2 => JFF routines after same principle as Tsuji jump=1 routines
!     3 => JFF Gibbs minimising routines => GEM-package
!     4 => ERC added GGCHEM code: Equilibrium chemistry down to 100 K
!-----------------------------------------------------------------------
      program scmarcs
      
      implicit real*16 (a-h,o-z)
      include 'parameter.inc'
      integer molh, jump
      character molname*4,osfil*60,sampling*3,pp_sph*3
      logical pf,pfe,pfd,fixros,itstop,quit,onemor
      common/cos/wnos(nwl),conos(ndp,nwl),wlos(nwl),wlstep(nwl),
     *  kos_step,nwtot,nosmol,newosatom,newosatomlist,
     *  nchrom,osfil(30),molname(30),sampling
      common /statec/ppr(ndp),ppt(ndp),pp(ndp),gg(ndp),zz(ndp),dd(ndp),
     *  vv(ndp),ffc(ndp),ppe(ndp),tt(ndp),tauln(ndp),ro(ndp),ntau,iter
      common /carc3/ f1p,f3p,f4p,f5p,hnic,presmo(33)
      common /cit/it,itmax
      common /ci4/dum,idum(96),molh,jump 
      common /cpf/pf,pfe,pfd,fixros,itstop
      common /carc1/istral,idrab1,idrab2,idrab3,idrab4,idrab5,
     *  idrab6,iarch
      common /newmo/newmod
      common /cmolrat/ fold(ndp,8),molold,kl
      common /fullequilibrium/ partryck(ndp,maxmol),
     *  xmettryck(ndp,maxmet),xiontryck(ndp,maxmet),
     *  partp(ndp,0:maxmol)
      common /cisph/isph
      common /clist/nlte
      common/ci5/abmarcs(17,ndp),anjon(17,5),h(5),part(17,5),
     *  dxi,f1,f2,f3,f4,f5,xkhm,xmh,xmy(ndp)
      namelist /outlist/ masslinf
      common /coutlist/pplist(843)
      common /clevprint/ prj2(ndp),masslinf
      common/cabinit/abinit(natms),kelem(natms),nelem
      common /cdustdata/ dabstable(max_lay,nwl),dscatable(max_lay,nwl),
     *    eps(max_eps,max_lay),temp(max_lay),pgas(max_lay),n_lay,idust
      common /cdustopac/ dust_abs(ndp,nwl), dust_sca(ndp,nwl)


! Initiation
      call gettime(0)
      
      open(unit=5,file='mxms7.input',status='old',readonly)
      open(unit=6,file='mxms7.output',status='unknown')
      open(unit=7,file='mxmodel.dat',status='unknown')
      open(unit=9,file='data/jonabs.dat',status='old',readonly)
      open(unit=16,file='arcivaaa.dat',status='old',readonly)
      open(unit=111,file='partf.dat')

      read(5,'(6(7x,i3,5x))') itmax,nprint,newmod,noarch,jontyp,idust
      write(6,'(5(7x,i3,5x))') itmax,nprint,newmod,noarch,jontyp,idust
      
      read(5,'(7x,i3,12x,a3,2(12x,i3),12x,a3)') jump,sampling,molold,
     *  nlte,pp_sph
      write(6,'(7x,i3,12x,a3,2(12x,i3),12x,a3)') jump,sampling,molold,
     *  nlte,pp_sph

      if(pp_sph.eq.'sph' .or. pp_sph.eq.'SPH') isph = 1

      write(7,'(a8,i3,a11,a3,a8,i2)') ' itemax=',itmax,' sampling= ',
     *  sampling,' molold=',molold

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
      implicit real*16 (a-h,o-z)
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
     &  xmettryck(ndp,maxmet),xiontryck(ndp,maxmet),
     *  partp(ndp,0:maxmol)
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
      implicit real*16 (a-h,o-z)
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
      SUBROUTINE AINV3(A,M)
      implicit real*16 (a-h,o-z)
C
C***** THIS SUBROUTINE EVALUATES THE INVERSE OF A
C***** SQUARE M*M MATRIX A
C***** THE DIMENSIONS OF THE ARRAYS MAY BE INTEGER VARIABLES WHEN USED**
C***** ON THE IBM 7094,BUT THEY MUST BE INTEGER CONSTANTS ON THE IBM 113
C
C
C      implicit real*16 (A-H,O-Z)
C
      DIMENSION A(8,8),C(8),IND(8)
C
CC
C
  100 AMAX=0.0
      DO 2 I=1,M
      IND(I)=I
      IF(ABS(A(I,1))-AMAX)2,2,3
    3 AMAX=ABS(A(I,1))
      I4=I
    2 CONTINUE
      MM=M-1
      DO 111 J=1,MM
      IF(I4-J)6,6,4
    4 ISTO=IND(J)
      IND(J)=IND(I4)
      IND(I4)=ISTO
      DO 5 K=1,M
      STO=A(I4,K)
      A(I4,K)=A(J,K)
      A(J,K)=STO
    5 CONTINUE
    6 AMAX=0.0
      J1=J+1
      DO 11 I=J1,M
      A(I,J)=A(I,J)/A(J,J)
      DO 10 K=J1,M
      A(I,K)=A(I,K)-A(I,J)*A(J,K)
      IF (K-J1)14,14,10
   14 IF(ABS(A(I,K))-AMAX)10,10,17
   17 AMAX=ABS(A(I,K))
      I4=I
   10 CONTINUE
   11 CONTINUE
  111 CONTINUE
      DO 140 I1=1,MM
      I=M+1-I1
      I2=I-1
      DO 41 J1=1,I2
      J=I2+1-J1
      J2=J+1
      W1=-A(I,J)
      IF(I2-J2)141,43,43
   43 DO 42 K=J2,I2
      W1=W1-A(K,J)*C(K)
   42 CONTINUE
  141 C(J)=W1
   41 CONTINUE
      DO 40 K=1,I2
      A(I,K)=C(K)
   40 CONTINUE
  140 CONTINUE
      DO 150 I1=1,M
      I=M+1-I1
      I2=I+1
      W=A(I,I)
      DO 56 J=1,M
      IF (I-J)52,53,54
   52 W1=0.0
      GO TO 55
   53 W1=1.0
      GO TO 55
   54 W1=A(I,J)
   55 IF(I1-1)156,156,57
   57 DO 58 K=I2,M
      W1=W1-A(I,K)*A(K,J)
   58 CONTINUE
  156 C(J)=W1
   56 CONTINUE
      DO 50 J=1,M
      A(I,J)=C(J)/W
   50 CONTINUE
  150 CONTINUE
C        DENNA RUTIN LOESER TAANSPORTEKVATIONEN MED FEAUTRIERS METOD,
      DO 60 I=1,M
   63 IF(IND(I)-I)61,60,61
   61 J=IND(I)
      DO 62 K=1,M
      STO=A(K,I)
      A(K,I)=A(K,J)
      A(K,J)=STO
   62 CONTINUE
      ISTO=IND(J)
      IND(J)=J
      IND(I)=ISTO
      GO TO 63
   60 CONTINUE
      RETURN
      END
*
      subroutine algebn(n)
      implicit real*16 (a-h,o-z)
*
*  algebn does the algebraic manipulations when solving the radiation
*  transoport+flux equation system.  The routine is called once per
*  wavelength point.  The results are subtracted from the matrices O
*  and Q, which describe the dependence of the flux and radiation
*  pressure equations on the mean intensities.  The matrices are
*  arranged as follows:
*
*  a2(1)  a3(1)                               b2(1)  b3(1)
*  a1(2)  a2(2)  a3(2)                        b1(2)  b2(2)  b3(2)
*         a1(3)  a2(3)  a3(3)                        b1(3)  b2(3)  b3(3)
*                  .................                         ............
*                                a1(n)  a2(n)
*
*
*  d2(1)
*  d1(2)  d2(2)
*         d1(3)  d2(3)
*                ..................
*                                d1(n)  d2(n)
*
*  c(1)
*         c(2)
*                c(3)
*                 .................
*                                c(n)
*  i.e.,
*
*                   A*dJ + B *dT + E *dPe = 0
*  sum over wavel:  C*dJ + O1*dT + O2*dPe = 0
*  sum over wavel:  D*dJ + Q1*dT + Q2*dPe = 0
*
*  where A and B are tridiagonal, D is bi-diagonal, C is diagonal,
*  and the O and Q matrices are full.
*
*  n is the rank of the matrices.  Operation count 17 n**2 flops.
*
*  Modifierad nov-72/aake.
*  Subtracts d*ainverse*b from q and c*ainverse*b from o, where c is a
*  diagonal matrix and o is a full matrix.
*
*  Modified oct-78/aake.
*  b is now tridiagonal.
*
*  Modified june-79/aake.
*  Following stein (priv. comm. -79), a2 contains the row sum of
*  elements instead of the diagonal eelement.  This eliminates the
*  need for double precision.
*
*  Modified may-89/aake.
*  Bug in the inversion of A fixed.  Syntax cleaned up to show the
*  function more clearly (fortran 77 syntax).
*
      include 'parameter.inc'
c
      dimension a1(ndp),a2(ndp),a3(ndp),d1(ndp),d2(ndp),c(ndp),h(ndp),
     * b1(ndp),b2(ndp),b3(ndp),g1(ndp),q1(ndp,ndp),o1(ndp,ndp),
     * e1(ndp),e2(ndp),e3(ndp),g2(ndp),q2(ndp,ndp),o2(ndp,ndp)
      common/space1/a1,a2,a3,d1,d2,b1,b2,b3,c,q1,o1,
     * e1,e2,e3,q2,o2
*
*  a1(k) stores the fraction of row k-1 to subtract from row k
*  a2(k) stores the normalization factor for row k
*  a3(k) stores the fraction of (normalized) row k+1 to subtract
*        from (normalized) row k
*
      do 1 k=1,n-1
*
*  scaling factor to apply to row k
*
        f=1./(a2(k)-a3(k))
*
*  fraction of row k to add to row k+1
*
        a1(k+1)=a1(k+1)*f
*
*  the new rowsum in row k+1
*
        a2(k+1)=a2(k+1)-a1(k+1)*a2(k)
*
*  scale super-diagonal, save scale factor
*
        a3(k)=a3(k)*f
        a2(k)=f
    1 continue
*
*  last scaling factor
*
      a2(n)=1./a2(n)
*
      do 7 l=1,n
*
*  We first initialize the l'th column vector in Ainv*B.  At this point
*  we also normalize by multiplying each row with a2(k).
*
      do 3 k=1,l-2
        g1(k)=0.
        g2(k)=0.
    3 continue        
*
      if (l.eq.1) then
        g1(l  )  =b2(l  )
        g2(l  )  =e2(l  )
      else
        g1(l-1)=b3(l-1)
        g2(l-1)=e3(l-1)
        g1(l  )  =b2(l  )-a1(l  )*g1(l-1)
        g2(l  )  =e2(l  )-a1(l  )*g2(l-1)
      endif
      if (l.lt.n) then
        g1(l+1)=b1(l+1)-a1(l+1)*g1(l  )
        g2(l+1)=e1(l+1)-a1(l+1)*g2(l  )
      endif
*
*  Eliminate the sub-diagonal in A by adding a fraction of row k-1 to
*  row k.  [2 flops average].
*
      do 4 k=l+2,n
        g1(k)=-a1(k)*g1(k-1)
        g2(k)=-a1(k)*g2(k-1)
    4 continue
*
*  Eliminate the super-diagonal by adding a fraction of (normailzed)
*  row k+1 to (normalized) row k.  [3 flops average].
*
      g1(n)=a2(n)*g1(n)
      g2(n)=a2(n)*g2(n)
      do 6 k=n-1,1,-1
        g1(k)=a2(k)*g1(k)-a3(k)*g1(k+1)
        g2(k)=a2(k)*g2(k)-a3(k)*g2(k+1)
    6 continue
*
*  Subtract the results from the Q and O matrices.
*  [12 flops]
*
      o1(1,l)=o1(1,l)-c(1)*g1(1)
      o2(1,l)=o2(1,l)-c(1)*g2(1)
      q1(1,l)=q1(1,l)-d2(1)*g1(1)
      q2(1,l)=q2(1,l)-d2(1)*g2(1)
      do 7 k=2,n
        o1(k,l)=o1(k,l)-c(k)*g1(k)
        o2(k,l)=o2(k,l)-c(k)*g2(k)
        q1(k,l)=q1(k,l)-d1(k)*g1(k-1)-d2(k)*g1(k)
        q2(k,l)=q2(k,l)-d1(k)*g2(k-1)-d2(k)*g2(k)
    7 continue
*
      return
      end
C
      SUBROUTINE ARCHIV(LUN)
      implicit real*16 (a-h,o-z)
C
C        THIS ROUTINE STORES ALL INTERESTING INFORMATION ON A MODEL ON
C        FORTRAN FILE IARCH. MOREOVER, IT PUNCES CARDS FOR BELL'S USE.
C
      include 'parameter.inc'
C
      CHARACTER*1  DAY(10)
      CHARACTER*10 DAG,KLOCK
      DIMENSION ABSKA(NDP),SPRIDA(NDP),ABSKTR(20),SPRTR(20)
      DIMENSION PRESMP(33)     !345??? was 33 ??? (UGJ/Sep.10.98)
      EQUIVALENCE (DAG,DAY(3))
C
C        COMMONS SHARED BY MAIN PROGRAM
C        THESE SHOULD BE MODIFIED BY NORDLUND FOR HIS PURPOSES.
      COMMON /TAUC/TAU(NDP),DTAULN(NDP),JTAU
      COMMON/STATEC/PRAD(NDP),PTURB(NDP),P(NDP),GG(NDP),ZZ(NDP),DD(NDP),
     &              VV(NDP),FLUXC(NDP),PE(NDP),T(NDP),
     &              TAUDUM(NDP),RO(NDP),NTAU,ITMAX
      COMMON /ROSSC/XKAPR(NDP),CROSS(NDP)
      COMMON /CTEFF/TEFF,FLUX
      COMMON/CG/G,KONSG
      COMMON /MIXC/PALFA,PBETA,PNY,PY
      COMMON /CSTYR/MIHAL,NOCONV
      COMMON/CLINE4/ILINE
      COMMON/CARC1/ISTRAL,IDRAB1,IDRAB2,IDRAB3,IDRAB4,IDRAB5,IDRAB6,
     &             IARCH
      COMMON/CFIL/IRESET(10),IDUM3,IDUM4
      COMMON /CVAAGL/XLB(500),W(500),NLB
      COMMON/CXLSET/XL(20,10),IDUM6,NL(10)
      COMMON /CSPHER/DIFLOG,RADIUS,RR(NDP),NCORE /CTAUM/TAUM
C
C        COMMON SHARED BY SOLVE
      COMMON /CARC2/TKORRM(NDP),FCCORR(NDP),FLUXME(NWL),TAU5(NDP),INORD
C
C        COMMONS SHARED BY JON
      COMMON/CI5/abmarcs(17,ndp),ANJON(17,5),H(5),PART(17,5),
     * DXI,F1,F2,F3,F4,F5,XKHM,XMH,XMY(ndp)
      COMMON/CI1/FL2(5),PARCO(45),PARQ(180),SHXIJ(5),TPARF(4),
     *XIONG(16,5),EEV,ENAMN(ndp),SUMH(ndp),XKBOL,
     *SUMM(ndp),NJ(16),IEL(16),NEL
      COMMON/CI4/ TMOLIM,IELEM(16),ION(16,5),MOLH,JUMP
      COMMON/CMOL2/PPKDUM(33),NMOL
      COMMON/CARC3/F1P,F3P,F4P,F5P,HNIC,PRESMO(33)
      COMMON /CMOLRAT/ FOLD(NDP,8),MOLOLD,KL
C
C        COMMON SHARED BY DETABS
      COMMON/CARC4/PROV(30),NPROVA,NPROVS,NPROV
      COMMON /CHAR/ ABNAME(30),SOURCE(30)
      common /cmtest/pgm1(ndp),pgm2(ndp),pem1(ndp),pem2(ndp)
     *     ,tm1(ndp),tm2(ndp),pgos(ndp)
      common /fullequilibrium/ partryck(ndp,maxmol),
     &  xmettryck(ndp,maxmet),xiontryck(ndp,maxmet),
     *  partp(ndp,0:maxmol)

      CHARACTER*8 SOURCE,ABNAME
      INTEGER MOLH, JUMP
C
C
      IARCH=LUN
      ISTAN2=1
      JSTAN2=NL(1)+1
      DO 60 K=1,JTAU
60    FCCORR(K)=3.14159*FCCORR(K)
C
C        WHICH PROGRAM WAS USED AND WHEN
C
C NOT SUPPORTED ON APOLLO      CALL ISODAT(DAY,KLOCK)
      WRITE(IARCH)INORD,DAG,KLOCK
C
C        STORE MODEL PARAMETERS
      WRITE(IARCH)TEFF,FLUX,G,PALFA,PNY,PY,PBETA,ILINE,ISTRAL,MIHAL,
     &            IDRAB1,IDRAB2,IDRAB3,IDRAB4,IDRAB5,IDRAB6,
     &            ITMAX,NEL,(abmarcs(I,1),I=1,NEL)
      WRITE(IARCH)JTAU,NCORE,DIFLOG,TAUM,RADIUS,(RR(K),K=1,JTAU)
C
C        STORE LAST TEMPERATURE-CORRECTION ARRAY
      WRITE(IARCH)JTAU,(TKORRM(K),K=1,JTAU),(FCCORR(K),K=1,JTAU)
C
C        COMPUTE AND STORE THERMODYNAMIC QUANTITIES AND DEPTHSCALES.
C        PUNCH BELL'S CARDS.
C
      do 2001 k=1,ndp
      K1=MIN0(K+1,ndp)
      PGm1(k)=P(K)-PRAD(K)-0.5*(PTURB(K)+PTURB(K1))
      pem1(k)=pe(k)
      tm1(k)=t(k)
 2001 continue
C
      DO 20 K=1,JTAU
      KL=K
C
C        JUST PRELIMINARY WE CHANGE TO FULL MOLECULAR EQUILIBRIUM. THIS
C        COMPUTED FOR DEMONSTRATION - IT IS NOT CONSISTENT WITH THE PRES
C        DENSITIES ETC GIVEN IN THE LISTING (WHERE ONLY H2 IS CONSIDERED
C        PORTANT). THE TAU(STANDARD) VALUES ARE NOT QUITE CONSISTENT, HO
C        THIS IS CERTAINLY QUITE UNIMPORTANT FOR MOST MODELS.
      MOLHO=MOLH
      MOLH=0
      CALL ABSKO(1,1,T(K),PE(K),ISTAN2,JSTAN2,ABSKA(1),SPRIDA(1))

****************18.12.04 Ch.H
*
* if JUMP >= 1 routine MOL dosn't work => no presmo
*    JUMP=1: Tsuji-routins are working instaed of MOL => partryck for molecules
*                                                xmettryck fot atoms
*    JUMP=2: JFF routines after same principle as Tsuji jump=1 routines
*    JUMP=3: JFF Gibs minimerings routines => call gem_package
*    JUMP=4: Added GGchem code by ERC
* ==>> more explanations in routine JON !!!
*
******

C This following piece doesn't really make sense if JUMP .ne.0 !!!
C -- but on the other hand probably it doesn't make harm, but be careful with the 
C  meaning of presmp  -- study later!!! (UGJ 25/1/01):

      IF (JUMP.GE.1) THEN 
 100  DO 101 I=1,33          !345
C  also in Tsujis routine is the number of molecules called NMOL
C  for no confusing I didn't use this name here
       PRESMP(I)=MAX(PARTRYCK(KL,I),1.D-99)
 101  CONTINUE
      ELSE 
C  NMOL=33 is number of molecules considered in the old MARCS chem.equilibrium
      DO 10 I=1,NMOL
       PRESMP(I)=MAX(PRESMO(I),1.D-99)
10    CONTINUE
      ENDIF
     
      HNIP=HNIC
C        CHANGE BACK AGAIN
CUGJ981018      MOLH=MOLHO
      CALL TERMO(T(K),PE(K),PRAD(K),PTOT,RRO,CP,CV,AGRAD,Q,U2)
      FORE=(ABSKA(1)+SPRIDA(1))/XKAPR(K)
      FURE=1./(XKAPR(K)*RRO)
      IF(K.GT.1)GO TO 11
      TAUS=FORE*TAU(1)
      Z=0.
      GO TO 12
   11 TAUS=TAUS+(TAU(K)-TAU(K-1))*(FORE+FOREM)*0.5
      Z=Z+(TAU(K)-TAU(K-1))*(FURE+FUREM)*0.5
   12 FOREM=FORE
      FUREM=FURE
C
      K1=MIN0(K+1,JTAU)
      PG=P(K)-PRAD(K)-0.5*(PTURB(K)+PTURB(K1))
      EMU=(1.38*RRO*T(K))/(1.67E-8*PG)
      if (U2.lt.0.) then
         write(7,*) ' in Archiv, U2= ',U2
         U2 = -U2
      end if
      U=SQRT(U2)
C
C INTERPOLATE THE CONVECTIVE FLUX AND VELOCITY TO THE DEPTH POINTS
C LOGARITHMICALLY
      IF(K.GT.1) GO TO 13
C K=1
      FCONV=0.
      V=0.
      GO TO 15
13    IF(K.EQ.JTAU) GO TO 14
C K"1,<JTAU
      YA=(TAU(K)-TAU(K-1))/(TAU(K+1)-TAU(K-1))
      YB=1.-YA
      FCONV=YA*FLUXC(K+1)+YB*FLUXC(K)
      IF(FLUXC(K).GT.0..AND.FLUXC(K+1).GT.0.) FCONV=
     &   EXP(YA*log(FLUXC(K+1))+YB*log(FLUXC(K)))
      V=YA*VV(K+1)+YB*VV(K)
      IF(VV(K).GT.0..AND.VV(K+1).GT.0.) V=
     &   EXP(YA*log(VV(K+1))+YB*log(VV(K)))
      GO TO 15
14    CONTINUE
C K=JTAU
      YA=(2.*TAU(K)-TAU(K-1)-TAU(K-2))/(TAU(K)-TAU(K-2))
      YB=1.-YA
      FCONV=YA*FLUXC(K)+YB*FLUXC(K-1)
      IF(FCONV.GT.FLUX) FCONV=FLUX
      V=YA*VV(K)+YB*VV(K-1)
15    CONTINUE
      ANCONV=FCONV/FLUX
C
C partial pressure of He is put into presmo(17)
        XNHE = abmarcs(2,k) / (XMH*XMY(k)) 
        PRESMO(17) = 1.38053d-16 * T(K) * XNHE * RRO
        PRESMP(17)=MAX(PRESMO(17),1.D-99)
C      write(7,7412) k,abund(2),xmh,xmy,xnhe,rro,t(k),presmp(17)
C7412  format(i3,1p7e11.4)

      WRITE(IARCH)K,TAU(K),TAUS,Z,T(K),PE(K),PG,PRAD(K),PTURB(K),
     &            XKAPR(K),RRO,EMU,CP,CV,AGRAD,Q,U,V,ANCONV,HNIP,
     &            NMOL,(PRESMP(I),I=1,NMOL)
   20 CONTINUE
C
      do 2002 k=1,ndp
      K1=MIN0(K+1,ndp)
      PGm2(k)=P(K)-PRAD(K)-0.5*(PTURB(K)+PTURB(K1))
      pem2(k)=pe(k)
      tm2(k)=t(k)
2002  continue
C
C        STORE TYPICAL IONIZATION EQUILIBRIA AND ABSORPTION COEFFICIENTS
C
      NLP=NL(1)+1
      WRITE(IARCH)(NJ(I),I=1,NEL),NLP,(XL(J,1),J=1,NLP),
     &   NPROV,NPROVA,NPROVS,(ABNAME(KP),SOURCE(KP),KP=1,NPROV)
      IF (NLP.GE.21) WRITE(7,*) ' *** ERROR: NLP>20 ****'
C
CUGJ STORE ABSORBTIONKOEFFICIENT, SCATTERING AND LAMBDA FOR SPECTRUM 
CUGJ CALCULATIONS
C
      DO 25 K=1,JTAU
      DO 27 J=1,NLP
      CALL ABSKO(1,1,T(K),PE(K),1,J,ABSKA(1),SPRIDA(1))
      ABSKTR(J)=ABSKA(1)
27    SPRTR(J)=SPRIDA(1)
      WRITE(IARCH) (ABSKTR(J),J=1,NLP)
      WRITE(IARCH) (SPRTR(J),J=1,NLP)
25    CONTINUE
C
      DO40 K=1,JTAU
      KL=K
      TAUK=log10(TAU(K))+10.01
      KTAU=TAUK
C     IF(ABS(TAUK-KTAU).GT.0.02) GO TO 40
      DO32 J=1,NLP
      CALL ABSKO(1,1,T(K),PE(K),1,J,ABSKA(1),SPRIDA(1))
      IF(J.GT.1)GO TO 31
      DO30 I=1,NEL
      NJP=NJ(I)
      WRITE(IARCH)K,TAU(K),T(K),PE(K),IEL(I),abmarcs(I,1),
     &    (ANJON(I,JJ),JJ=1,NJP),(PART(I,JJ),JJ=1,NJP)
   30 CONTINUE
   31 WRITE(IARCH)K,TAU(K),(PROV(KP),KP=1,NPROV),ABSKA(1),SPRIDA(1)
   32 CONTINUE
   40 CONTINUE
C
C        STORE FLUXES
C
      WRITE(IARCH)NLB,(XLB(J),FLUXME(J),J=1,NLB),(W(J),J=1,NLB)
C
      OPEN(UNIT=33,FILE='FLUX.DAT',STATUS='unknown')
      DO 324 J=1,NLB
      FLUXME(J)=3.14159*FLUXME(J)
      WRITE(33,333) XLB(J),FLUXME(J)
  324 CONTINUE
      CLOSE(33)
333   FORMAT(1P2E15.5)
C
      WRITE(6,50) LUN
50    FORMAT(' ARCHIV RECORDS WRITTEN ON LUN',I3)
C
C
      RETURN
      END
C
      FUNCTION BPL(T,X)
      implicit real*16 (a-h,o-z)
C
C T*X must be greater than 1.7e6 to give BPL > 1.e-37
C BPL therefore limited to be > 1.e-30. UGJ 900510
C (generally changed to > 1.e-20 due to conv. problems. UGJ 961230)
C
      COMMON /BPLC/EX,X5
      DATA CP/1.191E27/,C2/1.438E8/

      X5=((X**2)**2)*(X/CP)
      EX=EXP(-C2/(T*X))
      BPL=EX/((1.-EX)*X5)
      BPL = MAX(1.0D-99,BPL)
      RETURN
C
      ENTRY DIVBP(T,X)
      X6=X5*X
      TEX=T*(1.-EX)
      BPL=C2*(EX/TEX)/(TEX*X6)
      BPL = MAX(1.0D-99,BPL)
      RETURN
      END
C     MARK 4 RELEASE NAG COPYRIGHT 1974
C     MARK 4.5 REVISED
C     THIS ROUTINE ATTEMPTS TO SOLVE A REAL POLYNOMIAL EQUATION
C     HAVING N COEFFICIENTS (DEGREE  EQUALS  N-1) USING THE SEARCH
C     ALGORITHM PROPOSED IN GRANT AND HITCHINS (1971) TO
C     LIMITING MACHINE PRECISION.  ON ENTRY THE COEFFICIENTS
C     OF THE POLYNOMIAL ARE HELD IN THE ARRAY A(N), WITH A(0)
C     HOLDING THE COEFFICIENT OF THE HIGHEST POWER.  ON NORMAL
C     ENTRY THE PARAMETER IFAIL HAS VALUE 0 (HARD FAIL) OR 1
C     (SOFT FAIL) AND WILL BE ZERO ON SUCCESFUL EXIT WITH
C     THE CALCULATED ESTIMATES OF THE ROOTS HELD AS
C     REZ(K)+IIMZ(K), K EQUALS 1(1)N-1, IN APPROXIMATE DECREASING
C     ORDER OF MODULUS.  THE VALUE OF TOL IS OBTAINED BY
C     CALLING THE NAG ROUTINE X02AAF.
C     ABNORMAL EXITS WILL BE INDICATED BY IFAIL HAVING
C     VALUE 1 OR 2.  THE FORMER IMPLIES THAT EITHER A(1) EQUALS 0
C     OR N.LT.2 OR N.GT.100.  FOR IFAIL  EQUALS  2, A POSSIBLE
C     SADDLE
C     POINT HAS BEEN DETECTED.  THE NUMBER OF COEFFICIENTS
C     OF THE REDUCED POLYNOMIAL IS STORED IN N AND ITS
C     COEFFICIENTS ARE STORED IN A(1) TO A(N), THE ROOTS
C     THUS FAR BEING STORED IN THE ARRAYS REZ AND IMZ
C     STARTING WITH REZ(N)+IIMZ(N).  AN IMMEDIATE RE-ENTRY
C     IS POSSIBLE WITH IFAIL UNCHANGED AND WITH A NEW
C     STARTING POINT FOR THE SEARCH HELD IN REZ(1)+IIMZ(1).
C     REF - J.I.M.A., VOL.8., PP122-129 (1971).
CCCC
C
C
      SUBROUTINE C02AEF(A, N, REZ, IMZ, TOL, IFAIL)
      implicit real*16 (a-h,o-z)
C
C      implicit real*16 (A-H,O-Z)
C
      INTEGER IFAIL, IND, N, I, K, II, I2, JTEMP
      REAL*16  IMZ,J,JX,NFUN
      DIMENSION A(N),B(100),C(100),REZ(N),IMZ(N)
      LOGICAL SAT,FLAG
      COMMON /AC02AE/ X, Y, R, RX, J, JX, SAT
C     ETEXT/DTEXT
C     DATA ONE/1.0/,A1P5/1.5/,ZERO/0.0/,P4Z1/1.0E-5/
      DATA ONE /1.0/, A1P5 /1.5/, ZERO /0.0/, P4Z1 /1.0E-5/
C     DATA TWO/2.0/,P5/0.5/,P2Z1/1.0E-3/,P1/0.1/
      DATA TWO /2.0/, P5 /0.5/, P2Z1 /1.0E-3/, P1 /0.1/
C     DATA P3Z2/2.0E-4/,FOUR/4.0/
C     DATA P3Z2/2.0E-4/,FOUR/4.0/
      DATA P3Z2 /2.0E-4/, FOUR /4.0/
      XXX = X02AAF(XXX)
C     THE ABOVE TEST WAS ADDED AT 4.5 TO PREVENT TOL BEING TOO
C     SMALL

      IF (TOL.LT.XXX) TOL = XXX
      FAC = ONE
      FLAG = IFAIL.EQ.2
      IF (FLAG) IFAIL = 1
      IND = 0
      TOL2 = TOL**A1P5
      IF (A(1).NE.ZERO .AND. N.GE.2 .AND. N.LE.100) GO TO 20
      IND = 0
      GO TO 720
   20 IF (A(N).NE.0.0) GO TO 40
      REZ(N-1) = ZERO
      IMZ(N-1) = ZERO
      N = N - 1
      GO TO 20
   40 SCALE = ZERO
C     FUNCTION/DFUNCTION
      DO 60 I=1,N
       IF (ABS(A(I)).GE.P4Z1) SCALE = SCALE + LOG(ABS(A(I)))
C     FUNCTION/DFUNCTION
   60 CONTINUE
C      K = IDINT(SCALE/(DBLE(N)*LOG(TWO))+P5)
      K = INT(SCALE/(DBLE(N)*LOG(TWO))+P5)
      SCALE = TWO**(-K)
      DO 80 I=1,N
       A(I) = A(I)*SCALE
       B(I) = A(I)
C     TEST FOR LOW ORDER POLYNOMIAL FOR EXPLICIT SOLUTION
   80 CONTINUE
      IF (N.GT.3) GO TO 100
      GO TO (720, 560, 580), N
  100 DO 160 I=2,N
       II = N - I + 2
       DO 120 K=2,II
       I2 = II - K + 1
       C(K-1) = B(II)*B(K) - B(1)*B(I2)
  120  CONTINUE
       IF (C(II-1).LT.-TOL) GO TO 200
       T = ONE
       IF (C(II-1).GE.ONE) T = ONE/C(II-1)
       JTEMP = II - 1
       DO 140 K=1,JTEMP
       B(K) = C(K)*T
  140  CONTINUE
  160 CONTINUE
      FAC = FAC*TWO
      SCALE = ONE
      I = N
  180 I = I - 1
      IF (I.LT.1) GO TO 100
      SCALE = SCALE*TWO
      A(I) = A(I)*SCALE
      B(I) = A(I)
      GO TO 180
  200 IF (.NOT.FLAG) GO TO 220
      X = REZ(1)
      Y = IMZ(1) + TOL
      FLAG = .FALSE.
      GO TO 240
  220 X = P2Z1
      Y = P1
  240 CALL C02AEZ(A, N, TOL)
      FUN = R*R + J*J
  260 G = RX*RX + JX*JX
      IF (G.GE.FUN*TOL2) GO TO 300
      IND = 0
      SCALE = ONE
      I = N
  280 I = I - 1
      IF (I.LT.1) GO TO 720
      SCALE = SCALE*FAC
      A(I) = A(I)/SCALE
      GO TO 280
  300 S1 = -(R*RX+J*JX)/G
      S2 = (R*JX-J*RX)/G
C     FUNCTION/DFUNCTION
      SIG = P3Z2
      S = SQRT(S1*S1+S2*S2)
      IF (S.LE.ONE) GO TO 320
      S1 = S1/S
      S2 = S2/S
C     VALID DIRECTION OF SEARCH HAS BEEN DETERMINED, NOW
C     PROCEED TO DETERMINE SUITABLE STEP
      SIG = SIG/S
  320 X = X + S1
      Y = Y + S2
  340 CALL C02AEZ(A, N, TOL)
      IF (SAT) GO TO 380
      NFUN = R*R + J*J
      IF (FUN-NFUN.GE.SIG*FUN) GO TO 360
      S1 = P5*S1
      S2 = P5*S2
      S = P5*S
      SIG = P5*SIG
      X = X - S1
      Y = Y - S2
      GO TO 340
  360 FUN = NFUN
      GO TO 260
  380 FUN = ONE/TOL2
      K = 0
C     FUNCTION/DFUNCTION
      IMZ(N-1) = Y*FAC
C     CHECK POSSIBILITY OF REAL ROOT
      IF (ABS(Y).GT.P1) GO TO 420
      S1 = Y
      Y = ZERO
      CALL C02AEZ(A, N, TOL)
      Y = S1
C     REAL ROOT ACCEPTED AND BOTH BACKWARD AND FORWARD DEFLATIONS
C     ARE PERFORMED WITH LINEAR FACTOR
      IF (.NOT.SAT) GO TO 420
      REZ(N-1) = X*FAC
      IMZ(N-1) = ZERO
      N = N - 1
      B(1) = A(1)
      C(N) = -A(N+1)/X
      DO 400 I=2,N
       B(I) = A(I) + X*B(I-1)
       II = N - I + 1
       C(II) = (C(II+1)-A(II+1))/X
  400 CONTINUE
C     COMPLEX ROOT ACCEPTED AND BOTH BACKWARD AND FORWARD
C     DEFLATIONS ARE PERFORMED WITH QUADRATIC FACTOR
      GO TO 460
  420 REZ(N-1) = X*FAC
      REZ(N-2) = X*FAC
      IMZ(N-2) = -IMZ(N-1)
      N = N - 2
      R = TWO*X
      J = -(X*X+Y*Y)
      B(1) = A(1)
      B(2) = A(2) + R*B(1)
      C(N) = -A(N+2)/J
      C(N-1) = -(A(N+1)+R*C(N))/J
      IF (N.EQ.2) GO TO 460
      DO 440 I=3,N
       B(I) = A(I) + R*B(I-1) + J*B(I-2)
       II = N - I + 1
       C(II) = -(A(II+2)-C(II+2)+R*C(II+1))/J
C     MATCHING POINT FOR COMPOSITE DEFLATION
  440 CONTINUE
C     FUNCTION/DFUNCTION
  460 DO 480 I=1,N
       NFUN = ABS(B(I)) + ABS(C(I))
C     FUNCTION/DFUNCTION
       IF (NFUN.LE.TOL) GO TO 480
       NFUN = ABS(B(I)-C(I))/NFUN
       IF (NFUN.GE.FUN) GO TO 480
       FUN = NFUN
       K = I
  480 CONTINUE
      IF (K.EQ.1) GO TO 520
      JTEMP = K - 1
      DO 500 I=1,JTEMP
       A(I) = B(I)
  500 CONTINUE
  520 A(K) = P5*(B(K)+C(K))
      IF (K.EQ.N) GO TO 40
      JTEMP = K + 1
      DO 540 I=JTEMP,N
       A(I) = C(I)
  540 CONTINUE
      GO TO 40
  560 REZ(1) = -A(2)/A(1)*FAC
      IMZ(1) = ZERO
      GO TO 700
  580 R = A(2)*A(2) - FOUR*A(1)*A(3)
      IF (R.GT.ZERO) GO TO 600
      REZ(2) = -P5*A(2)/A(1)*FAC
C     FUNCTION/DFUNCTION
      REZ(1) = REZ(2)
      IMZ(2) = P5*SQRT(-R)/A(1)*FAC
      IMZ(1) = -IMZ(2)
      GO TO 700
  600 IMZ(1) = ZERO
      IMZ(2) = ZERO
C     FUNCTION/DFUNCTION
      IF (A(2)) 620, 640, 660
  620 REZ(1) = P5*(-A(2)+SQRT(R))/A(1)*FAC
      GO TO 680
  640 REZ(1) = -P5*A(2)/A(1)*FAC
C     FUNCTION/DFUNCTION
      GO TO 680
  660 REZ(1) = P5*(-A(2)-SQRT(R))/A(1)*FAC
  680 REZ(2) = A(3)/(REZ(1)*A(1))*FAC*FAC
  700 N = 1
  720 IFAIL = IND
      RETURN
      END
C
C     MARK 4 RELEASE NAG COPYRIGHT 1974.
C     MARK 4.5 REVISED
C     EVALUATES R,RX,J,JX AT THE POINT X+IY AND APPLIES THE ADAMS
C     TEST.
C     THE BOOLEAN VARIABLE SAT IS GIVEN THE VALUE TRUE IF THE TEST
C     IS
C     SATISFIED.
CCCC
C
C
      SUBROUTINE C02AEZ(A, N, TOL)
      implicit real*16 (a-h,o-z)
C
C      implicit real*16 (A-H,O-Z)
C
      INTEGER N, K
      REAL*16  J,JX
      DIMENSION A(N)
      LOGICAL SAT
C
C     ETEXT/DTEXT
      COMMON /AC02AE/ X, Y, R, RX, J, JX, SAT
      DATA TWO/2.0/,ZERO/0.0/,P8/0.8/,TEN /10.0/,A8/8.0/
C

      P = -TWO*X
C     FUNCTION/DFUNCTION
      Q = X*X + Y*Y
      T = SQRT(Q)
      A2 = ZERO
      B2 = ZERO
      B1 = A(1)
C     FUNCTION/DFUNCTION
      A1 = A(1)
      C = ABS(A1)*P8
      N = N - 2
      DO 20 K=2,N
       A3 = A2
       A2 = A1
C     FUNCTION/DFUNCTION
       A1 = A(K) - P*A2 - Q*A3
       C = T*C + ABS(A1)
       B3 = B2
       B2 = B1
       B1 = A1 - P*B2 - Q*B3
   20 CONTINUE
      N = N + 2
      A3 = A2
      A2 = A1
      A1 = A(N-1) - P*A2 - Q*A3
      R = A(N) + X*A1 - Q*A2
      J = A1*Y
      RX = A1 - TWO*B2*Y*Y
C     FUNCTION/DFUNCTION
      JX = TWO*Y*(B1-X*B2)
C     FUNCTION/DFUNCTION
      C = T*(T*C+ABS(A1)) + ABS(R)
      SAT = (SQRT(R*R+J*J)).LT.((TEN*C-A8*(ABS(R)+ABS(A1)*T)+TWO* ABS(X*
     *A1))*TOL)
      RETURN
      END
C
      SUBROUTINE CLOCK
      implicit real*16 (a-h,o-z)
C
C  TIME SINCE LAST CALL,ACCUMULATED EXECUTION TIME
C  AND TIME REMAINING OF REQUESTED CPU TIME (for Cyber machines ?)
C
C
C     DATA IS/0/
CC
CC
C     IF (IS.NE.1) THEN
C       CALL MSLEFT(MS)
C       IS=1
C       MPP=MS
C     END IF
CC
C     CALL MSLEFT(MP)
C     SEC=(MPP-MP)/1000.
C     ACC=(MS-MP)/1000.0
C     REM=MP/1000.0
C     MPP=MP
CCC
      SEC=0.
      ACC=0.
        REM=0.
CCC
C
C      WRITE(6,50) SEC,ACC,REM
C 50   FORMAT(' TIME ',F8.3,4X,'ACCUMULATED ',F8.3,4X,
C     &           'REMAINING ',F8.3)
C
      RETURN  
      END
C
C
C         
      subroutine cclock
      implicit real*16 (a-h,o-z)
      call clock
      END
C
      SUBROUTINE COSCAL(LUN)
      implicit real*16 (a-h,o-z)
C
      include 'parameter.inc'
C
      COMMON /STATEC/DUM1(9*NDP),TT(NDP),DUM2(NDP),RO(NDP),NTAU,ITER
      COMMON /CISPH/ISPH
      DIMENSION TP(NDP)
C
C READ OLD MODEL

      CALL OLDSTA
C
C READ SCALE MODEL
      LUN1=LUN+1
      READ(LUN1) (DUMS,K=1,360),(TP(K),K=1,NTAU)
C
C PRINT
      WRITE(7,50) LUN,LUN1
50    FORMAT('0OLD MODEL FROM LUN',I3,' COMBINED WITH TEMPERATURE '
     &,'FROM LUN',I3,' FOR TAU < 1.0')
C
C ASSEMBLE
      K1=24
      DO 100 K=1,K1
100   TT(K)=TP(K)
      K2=K1+1
      TTK2=TT(K2)
      DO 101 K=K2,NTAU
101   TT(K)=TT(K)+TP(K2)-TTK2
C
C INTEGRATE PRESSURE EQUATION
        if (isph.eq.1) then 
          CALL TRYCK_sph
        else 
          CALL TRYCK
        end if
      ITER=0
      RETURN
      END
C
      SUBROUTINE DETABS(J,JP,NTP,IOUTR)
      implicit real*16 (a-h,o-z)
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
C
      SUBROUTINE DUBINT(NXSKAL,XSKAL,NYSKAL,YSKAL,IYBEG,IYEND,N,X,Y,
     &                                            FACTOR,IX,IY1,IY2)
      implicit real*16 (a-h,o-z)
C
C        THIS ROUTINE COMPUTES INTERPOLATION FACTORS FOR INTERPOLATION I
C        NOT NECESSARILY EQUIDISTANT, TWO DIMENSIONAL TABLE. SCALES
C           XSKAL (NXSKAL POINTS)
C           YSKAL(NYSKAL POINTS)
C        THE TABLE IS ONLY DEFINED FOR Y-VALUES YSKAL(L), WHERE L IS
C        WITHIN THE INTERVAL IYBEG(NX) TO IYEND(NX) FOR A GIVEN XSKAL(NX
C        ARGUMENTS ARE GIVEN IN X AND Y (N POINTS).
C        RESULTING FACTOR FOR POINT K IS PUT IN FACTOR(K,1-2,1-2)
C        STARTING POINTS AT INTERPOLATION IN IX(K), REFERRING TO THE XSC
C        IN IY1(K) AND IY2(K), REFERRING TO THE RESTRICTED Y SCALES DEFI
C        XSCALE(IX(K)) AND XSCALE(IX(K)+1), RESPECTIVELY.
C        INTERPOLATIONS AND EXTRAPOLATIONS ARE  L I N E A R .
C
      include 'parameter.inc'
C
      DIMENSION XSKAL(30),YSKAL(30),X(NDP),Y(NDP),FACTOR(NDP,2,2),
     &          IYBEG(30),IYEND(30),IY1(NDP),IY2(NDP),IX(NDP)
C
      DO5 K=1,N
      DO1 J=2,NXSKAL
      JMEM1=J
      IF(X(K).LT.XSKAL(J))GO TO 2
    1 CONTINUE
    2 IX(K)=JMEM1-1
      JM1=JMEM1-1
      DO3 J=2,NYSKAL
      JMEM2=J
      IF(Y(K).LT.YSKAL(J))GO TO 4
    3 CONTINUE
    4 IY=JMEM2-1
      IY1(K)=IY+1-IYBEG(JM1)
      IY1(K)=MIN0(IY1(K),IYEND(JM1)-IYBEG(JM1))
      IY1(K)=MAX0(IY1(K),1)
      IY2(K)=IY+1-IYBEG(JMEM1)
      IY2(K)=MIN0(IY2(K),IYEND(JMEM1)-IYBEG(JMEM1))
      IY2(K)=MAX0(IY2(K),1)
      I1=IY1(K)+IYBEG(JM1)-1
      I2=IY2(K)+IYBEG(JMEM1)-1
      DX=(X(K)-XSKAL(JMEM1-1))/(XSKAL(JMEM1)-XSKAL(JMEM1-1))
      DY1=(Y(K)-YSKAL(I1))/(YSKAL(I1+1)-YSKAL(I1))
      DY2=(Y(K)-YSKAL(I2))/(YSKAL(I2+1)-YSKAL(I2))
      FACTOR(K,1,1)=(1.-DX-DY1+DX*DY1)
      FACTOR(K,2,1)=(1.-DY2)*DX
      FACTOR(K,1,2)=(1.-DX)*DY1
    5 FACTOR(K,2,2)=DX*DY2
      RETURN
      END
C
      SUBROUTINE DUMIN
      implicit real*16 (a-h,o-z)
C
C 'DUMIN' READS THE FIRST SETS OF CARDS IN THE JONABS-DATA. THESE DATA
C ARE USED IN ATMOS BUT NOT IN MARCS.  *NORD*
C
      COMMON /UTPUT/IREAD,IWRIT
C
      READ(IREAD,50)I,J
      N=4+I/8+J/16+2*(I/16)
      DO 100 I=1,N
100   READ(IREAD,51)A
      DO 101 J=1,3
      READ(IREAD,50)N
      DO 101 I=1,N
101   READ(IREAD,51)A
      RETURN
50    FORMAT(2I5)
51    FORMAT(A4)
      END
C
      FUNCTION FOUR(Y,X,K,N)
      implicit real*16 (a-h,o-z)
C
      include 'parameter.inc'
C
      DIMENSION X(NDP),Y(NDP)
C
C FOURPOINT LAGRANGE INTERPOLATION TO FIND Y(K-.5). Y AND X OF LENGTH N.
C
C START ADRESS AND NONCENTERING
      IF(K.EQ.2) GOTO 3
      KK=MIN0(MAX0(K-3,0),N-4)
      II=0
      IF(K.LE.2) II=-1
      IF(K.EQ.N) II=1
      XX=.5*( X(KK+II+2)+X(KK+II+3))
C
C I LOOP
      FOUR=0.
      DO 1 I=1,4
      PROD=Y(KK+I)
C
C J LOOP
      DO 2 J=1,4
      IF(J.EQ.I) GO TO 2
      PROD=PROD*(XX-X(KK+J))/(X(KK+I)-X(KK+J))
2     CONTINUE
C
1     FOUR=FOUR+PROD
      RETURN
C
C LINEAR INTERPOLATION AT K=2. 780605.
3     FOUR=.5*(Y(1)+Y(2))
      RETURN
      END
C
      SUBROUTINE GAUSI(K,A,B,AI,XMYI)
      implicit real*16 (a-h,o-z)
C
C        RUTINEN GER VIKTER OCH INTEGRATIONSPUNKTER FOER GAUSSINTEGRATIO
C        MELLAN A OCH B - B AER OEVRE GRAENS , A NEDRE. KAELLA FOER DATA
C        LOWAN, DAVIDS, LEVENSON,  BULL AMER MATH SOC  48 SID 739  (1942
C        AI=VIKTER, XMYI=INTEGRATIONSPUNKTER.
C        INTEGRATIONSORDNING K.  K VAELJES MELLAN 2 OCH 10.
C
      DIMENSION AI(K),XMYI(K),AP(29),XMYP(29),INDOV(9)
      DOUBLE PRECISION AP,XMYP
C               10 DATAKORT FOER AP, 9 FOER XMYP OCH 1 FOER INDOV
      DATA AP/1.0,0.55555555555555,.88888888888888,.347854845137
     *,0.65214515486254,0.23692688505618,0.47862867049936,
     * 0.56888888888888,0.17132449237917,0.36076157304813,
     * 0.46791393457269,0.12948496616887,0.27970539148927,
     * 0.38183005050511,0.41795918367346,0.10122853629037,
     * 0.22238103445337,0.31370664587788,0.36268378337836,
     * 0.08127438836157,0.18064816069485,0.26061069640293,
     * 0.31234707704000,0.33023935500126,0.06667134430868,
     * 0.14945134915058,0.21908636251598,0.26926671930999,
     * 0.29552422471475/
      DATA XMYP/
     *0.57735026918962,.77459666924148,.0,0.86113631159405,
     *0.33998104358485,.90617984593866,.53846931010568,.0,
     *0.93246951420315,.66120938646626,.23861918608319,
     *0.94910791234275,.74153118559939,.40584515137739,.0,
     *0.96028985649753,.79666647741362,.52553240991632,
     *0.18343464249565,.96816023950762,.83603110732663,
     *0.61337143270059,.32425342340380,.0,0.97390652851717,
     *0.86506336668898,.67940956829902,.43339539412924,
     *0.14887433898163/
      DATA INDOV/1,3,5,8,11,15,19,24,29/
      IF(K.EQ.1)GO TO 7
      KUD=0
      FLK=DFLOAT(K)/2.
      K2=K/2
      FK=DFLOAT(K2)
      IF(ABS(FLK-FK)-1.E-7)2,1,1
    1 K2=K2+1
      KUD=1
    2 IOEV=INDOV(K-1)
      INED=IOEV-K2
      DO3 I=1,K2
      IP=INED+I
      XMYI(I)=-XMYP(IP)*(B-A)*0.5+(B+A)*0.5
    3 AI(I)=(B-A)*0.5*AP(IP)
      K2=K2+1
      DO4 I=K2,K
      IP=IOEV+K2-I
      IF(KUD)6,6,5
    5 IP=IP-1
    6 CONTINUE
      XMYI(I)= XMYP(IP)*(B-A)*0.5+(B+A)*0.5
    4 AI(I)=(B-A)*0.5*AP(IP)
      RETURN
    7 XMYI(1)=(B+A)*0.5
      AI(1)=B-A
      RETURN
      END
C
      SUBROUTINE H2OPAC(OMEGA,T,PROPAC)
      implicit real*16 (a-h,o-z)
C
      A=7.02391+1.3380*log10(T)
      A=10.**A
      A=1./A
      B=91.67+0.1033*T
      C=(15.57906-2.06158*log10(T)-0.477352*(log10(T))**2)/1.E7
      D=2.31317+3.8856E-4*T
      D=10**D
      OMEGAC=274.3+.2762*T
C     WRITE(7,4) A,B,C,D,OMEGAC
    4 FORMAT(3X,5E15.5)
      IF (OMEGA.GE.OMEGAC)GOTO1
      OMEGAT=A*OMEGA**2*EXP(-OMEGA/B)
      GOTO2
    1 OMEGAT=C*EXP(-OMEGA/D)
    2 CONTINUE
C     WRITE(7,3) T,OMEGA,OMEGAT
    3 FORMAT(3X,F10.0,2E15.5)
    5 CONTINUE
      E=4.2432E-6-2.8854E-7*log(T)
      F=1.2171E+5+258.28*T
      G=2.5830E-4-4.3429E-8*T
      H=1.1332E-2-1.1943E-3*log(T)
      OMEGAP=-2973.3+600.73*log(T)
      OMEGT2=1.5*OMEGAP
C     WRITE(7,11) E,F,G,H,OMEGAP,OMEGT2
   11 FORMAT(3X,6E15.5)
      IF(OMEGA.GT.OMEGT2)GOTO13
      OMEGAR=E*EXP(-(OMEGA-OMEGAP)**2/F)
      GOTO14
   13 OMEGAR=G*EXP(-H*OMEGA)
   14 CONTINUE
C     WRITE(7,15) T,OMEGA,OMEGAR
   15 FORMAT(3X,F10.0,2E15.5)
   12 CONTINUE
      GAM=6.0273E-10+2.2905E-13*T+4.0848E-17*T**2
      W1=363.96+1.3530*T-3.5807E-4*T**2+3.3618E-8*T**3
      XJ=161.45/T-2.6996-1.9537E-4*T
      XJ=10.**XJ
      W2=-108626./T+697.59+0.14353*T
      XT=28.765/T-9.0461+1.1552E-4*T
      XT=10.**XT
      XG=1.4860+0.44462*log10(T)
      XG=10.**XG
      XNYC=-5.1972+2.1*log10(T)
      XNYC=-(10.**XNYC-4172.)
      W232=1.5*W2
C     WRITE(7,21) GAM,W21,XJ,W2,XT,XG,XNYC,W232
   21  FORMAT(3X,8E15.5)
      XNY=OMEGA
      IF(XNY.GE.XNYC)GOTO23
      OMEGAK=GAM*W1**2*EXP(XJ*(XNY-XNYC))/((XNY-XNYC)**2+W1**2)
      GOTO24
   23 IF(XNY.GT.(XNYC+W232))GOTO25
      OMEGAK=GAM*W2**2/((XNY-XNYC)**2+W2**2)
      GOTO24
   25 OMEGAK=XT*EXP(-(XNY-XNYC)/XG)
   24 OMEGAK=OMEGAK*OMEGA
C     WRITE(7,26) T,XNY,OMEGAK
   26 FORMAT(3X,2F10.0,E15.5)
   22 CONTINUE
      PROPAC=OMEGAR+OMEGAT+OMEGAK
      RETURN
      END
C
      SUBROUTINE HEOPAC(OMEGA,T,PROPAC)
      implicit real*16 (a-h,o-z)
C
      A=7.02391+1.3380*log10(T)
      A=10.**A
      A=1./A
      B=91.67+0.1033*T
      C=(15.57906-2.06158*log10(T)-0.477352*(log10(T))**2)/1.E7
      D=2.31317+3.8856E-4*T
      D=10**D
      OMEGAC=274.3+.2762*T
      IF (OMEGA.GE.OMEGAC)GOTO1
      OMEGAK=A*OMEGA**2*EXP(-OMEGA/B)
      GOTO2
    1 OMEGAK=C*EXP(-OMEGA/D)
    2 CONTINUE
      OMEGAK=1.78*OMEGAK
      E=4.2432E-6-2.8854E-7*log(T)
      F=1.2171E+5+258.28*T
      G=2.5830E-4-4.3429E-8*T
      H=1.1332E-2-1.1943E-3*log(T)
      OMEGAP=-2973.3+600.73*log(T)
      OMEGAT=1.5*OMEGAP
      IF(OMEGA.GT.OMEGAT)GOTO13
      OMEGAX=E*EXP(-(OMEGA-OMEGAP)**2/F)
      GOTO12
   13 OMEGAX=G*EXP(-H*OMEGA)
   12 CONTINUE
      OMEGAX=.10*OMEGAX
      DEL2=-4.033E+4+263.93*T
      DEL=SQRT(DEL2)
      APRIM=-7.7245+0.4246*log10(T)
      APRIM=10.**APRIM
      AA=1./(1.125E+9+1.5866E+4*T+24.267*T**2)
      BB=1.2044+.4956*log10(T)
      BB=10.**BB
      OMEGAZ=4161.1
      IF(OMEGA.GE.OMEGAZ)GOTO22
      OMEGAY=(APRIM*DEL*OMEGA*EXP((OMEGA-OMEGAZ)/(.6952*T)))/((OMEGA-
     & OMEGAZ)**2 +DEL2)
      GOTO23
   22 IF(OMEGA.GT.(OMEGAZ+1.5*DEL))GOTO24
      OMEGAY=(APRIM*DEL*OMEGA)/((OMEGA-OMEGAZ)**2+DEL2)
      GOTO23
   24 OMEGAY=AA*OMEGA*EXP(-(OMEGA-OMEGAZ)/BB)
   23 CONTINUE
      PROPAC=OMEGAK+OMEGAX+OMEGAY
      RETURN
      END
C
      SUBROUTINE INABS(IOUTS)
      implicit real*16 (a-h,o-z)
C
C        THIS ROUTINE  READS ABSORPTION COEFFICIENT TABLES AND INTER/EXTRA-
C        POLATES THEM TO OUR WAVELENGTHS GIVEN IN XL. THE INTERPOLATION IS
C        PERFORMED SEPARATELY FOR EACH WAVELENGTH SET.
C
C        NKOMP IS THE NUMBER OF COMPONENTS IN THE FULL TABLE.
C        NEXTL SHOULD BE GREATER THAN ZERO IF A PRINT-OUT IS WANTED ON EXTRA-
C              POLATION IN WAVELENGTH,
C        NUTZL IF PRINT-OUT IS WANTED WHEN WE PUT THE COEFFICIENT =0 OUTSIDE THE
C              WAVELENGTH REGION OF THE TABLES.
C        NEXTT AND NUTZL ARE THE CORRESPONDING QUANTITIES ON INTERPOLATION IN
C              T, MADE IN SUBROUTINE TABS.
C        NULL  SHOULD BE GREATER THAN ZERO IF A PRINT-OUT IS WANTED (FROM SUB-
C              ROUTINE ABSKO) WHEN A COEFFICIENT IS FOUND TO BE LESS THAN ZERO
C              ON INTERPOLATION IN T AND THEREFORE PUT EQUAL TO ZERO.
C
C        FOR EACH COMPONENT THE FOLLOWING PARAMETERS MUST BE SPECIFIED
C        ABNAME IS THE NAME OF, OR A SYMBOL FOR, THE ABSORPTION MECHANISM.
C        SOURCE INDICATES THE SOURCE OR REFERENCE OF THE DATA
C
C        1. PARAMETERS FOR THE WAVELENGTH INTERPOLATION.
C          ILOGL SHOULD BE GREATER THAN ZERO IF INTERPOLATION IN WAVELENGTH IS
C                TO BE PERFORMED ON THE LOGARITHMIC ABSORPTION COEFFICIENTS
C                (WITH SUBSEQUENT EXPONENTIATION OF THE RESULTS - HERE IF ILOGT
C                IS EQUAL TO ZERO OR IN SUBROUTINE ABSKO IF ILOGT IS GREATER
C                THAN ZERO). OTHERWISE INTERPOLATION IN WAVELENGTH IS MADE
C                DIRECTLY ON THE ABSORPTION COEFFICIENTS THEMSELVES.
C          KVADL SHOULD BE GREATER THAN ZERO IF QUADRATIC INTERPOLATION IN
C                WAVELENGTH IS WANTED. OTHERWISE INTERPOLATION WILL BE LINEAR
C          MINEX SHOULD BE GT 0 IF LINEAR EXTRAPOLATION (INSTEAD OF PUTTING THE
C                COEFFICIENT = 0) IS WANTED TOWARDS SHORTER WAVELENGTHS.
C          MAXEX, CORRESPONDING TOWARDS LONGER WAVELENGTHS.
C          NLATB IS THE NUMBER OF WAVELENGTH POINTS OF THE ABSORPTION COEFFI-
C                CIENT TABLE TO BE READ.
C          XLATB ARE THOSE WAVELENGTHS. THEY SHOULD BE GIVEN IN INCREASING ORDER
C
C        2. PARAMETERS FOR THE TEMPERATURE INTERPOLATION.
C          ILOGT, KVADT, MINET, MAXET AND NTETB ARE THE T-INTERPOLATION
C                ANlogUES TO ILOGL-NLATB.
C          ITETA IS PUT GREATER THAN ZERO WHEN TETA VALUES (TETA=5040./T) ARE
C                GIVEN IN XTET INSTEAD OF TEMPERATURES.
C          XTET ARE THE TEMPERATURE (TETA) VALUES OF THE ABSORPTION
C                COEFFICIENT TABLE TO BE READ. THE XTET VALUES SHOULD BE GIVEN
C                IN INCREASING ORDER AND EQUIDISTANTLY, HOWEVER (IELMAX-1)
C                CHANGES OF THE INTERVAL ARE ALLOWED. THE PROGRAM CHECKS  THAT
C                THIS NUMBER IS NOT EXCEEDED.
C        XKAP IS THE ABSORPTION COEFFICIENT TABLE FOR THE ACTUAL COMPONENT. THE
C                WAVELENGTHS INCREASES MORE RAPIDLY THAN T (TETA).
C
C        THE TABLES FOR  T E M P E R A T U R E - I N D E P E N D E N T
C        C O M P O N E N T S  S H O U L D  B E  P U T  F I R S T .
C           THE RESULTING TABLE IS PUT IN ABKOF. HERE T (TETA) INCREASES MORE
C        RAPIDLY THAN XLA, WHICH INCREASES MORE RAPIDLY THAN KOMP. IF THE RESULT
C        OF THE INTERPOLATION IS ZERO FOR A CERTAIN XLA(J) AND KOMP, THIS IS NOT
C        PUT IN ABKOF. INSTEAD A NOTE IS MADE IN KOMPLA (KOMPLA(NLB*(KOMP-
C        1)+J) IS PUT EQUAL TO ZERO). OTHERWISE THE KOMPLA VALUE TELLS WHERE IN
C        ABKOF THE TABLE FOR THE COMPONENT KOMP AND THE WAVELENGTH J BEGINS.
C
C        A DETAILED PRINT-OUT IS GIVEN IF IOUTS IS GREATER THAN ZERO.
C
C
C        DIMENSIONS NECESSARY
C        ABKOF(NABDIM),ABNAME(NKOMP),DELT(NKOMP,IELMAX),IDEL(NKOMP),
C        IDISKV(MAX(NLATB)),ILOGTA(NKOMP),IRESET(NSET),ISVIT(NKOMP),ITETA(NKOMP)
C        KOMPLA(MAX(NL)*NKOMP),KVADT(NKOMP),MAXET(NKOMP),MINET(NKOMP),
C        NL(NSET),NTAET(NKOMP),NTM(NKOMP,IELMAX),SOURCE(NKOMP),
C        TBOLT(NKOMP,IELMAX),XKAP(MAX(NLATB),MAX(NTETB)),XL(MAX(NL),NSET)
C        XLA(MAX(NL)),XLA3(MAX(NL)),XLATB(MAX(NLATB)),XTET(MAX(NTETB)),
C        XTETP(MAX(NTETB))
C
C        THE DIMENSIONS ARE LOWER LIMITS
C        IELMAX IS THE MAXIMUM NUMBER OF DIFFERENT T INTERVALS (GIVEN BELOW) IN
C               ANY ABSORPTION COEFFICIENT TABLE.
C        NABDIM IS THE DIMENSION OF THE ABKOF ARRAY (GIVEN BELOW).
C        NKOMP IS THE NUMBER OF 'COMPONENTS', I.E. EQUAL TO THE NUMBER OF
C               DIFFERENT ABSORPTION COEFFICIENT TABLES TO BE READ.
C        NL(I)  IS THE NUMBER OF WAVELENGTHS IN THE WAVELENGTH SET I.
C        NLATB(KOMP) IS THE NUMBER OF WAVELENGTH POINTS IN THE TABLE TO BE READ
C               FOR THE COMPONENT KOMP.
C        NSET   IS THE NUMBER OF WAVELENGTH SETS.
C        NTETB  IS THE NUMBER OF TEMPERATURE POINTS IN THE TABLE FOR THE COM-
C               PONENT BEING CONSIDERED.
C
C
      DIMENSION IDISKV(40),XLATB(40),XTET(30),NTAET(30),XKAP(40,30),
     *XLA3(20),XLA(20),XTETP(30)
      COMMON/CARC4/PROV(30),NDUM(3)
      COMMON /CHAR/ ABNAME(30),SOURCE(30)
      COMMON/UTPUT/IREAD,IWRIT
      COMMON/CA1/DELT(30,2),TBOT(30,2),IDEL(30),ISVIT(30),ITETA(30),
     *KVADT(30),MAXET(30),MINET(30),NTM(30,2),NEXTT,NUTZT
      COMMON/CA2/ABKOF(4000),KOMPLA(600),KOMPR,KOMPS,NKOMP
      COMMON/CA3/ILOGTA(30),NULL
      COMMON/CFIL/IRESET(10),ISLASK,IREAT
      COMMON/CXLSET/XL(20,10),NSET,NL(10)
      CHARACTER*8 ABNAME,SOURCE
C
C        IELMAX IS THE MAXIMUM NUMBER OF DIFFERENT T INTERVALS IN THE XKAP-
C        TABLE. THE DIMENSIONS OF TBOT, DELT AND NTM ARE AFFECTED BY THIS NUMBER
      IELMAX=2
C        THE DIMENSION OF THE ABKOF ARRAY
      NABDIM=4000
      DO705 L=1,30
  705 XTETP(L)=0.
      IF(IOUTS.GT.0)WRITE(IWRIT,229)
C
      READ(IREAT,101)NKOMP,NEXTL,NUTZL,NEXTT,NUTZT,NULL
C
c           print *,'inabs,    NKOMP=',nkomp
C
      KOMPR=0
      REWIND ISLASK
C
C        LOOP OVER COMPONENTS STARTS (THE 'FIRST KOMP-LOOP')
      DO720 KOMP=1,NKOMP
      READ(IREAT,105)ABNAME(KOMP),SOURCE(KOMP)
      READ(IREAT,102)ILOGL,KVADL,MINEX,MAXEX,NLATB
      READ(IREAT,103)(XLATB(J),J=1,NLATB)
C
C        WE FIND THE DISCONTINUITIES IN WAVELENGTH
C        A DISCONTINUITY IN A TABLE IS DEFINED BY TWO WAVELENGTH POINTS
C        WITHIN LESS THAN TWO ANGSTROEMS.
      IDISK=0
      IDISKV(1)=0
      DO700 J=2,NLATB
      IDISKV(J)=0
      IF((XLATB(J)-XLATB(J-1)).GE.2.)GO TO 700
      IDISKV(J-1)=1
      IDISKV(J)=1
      IDISK=1
  700 CONTINUE
C
C        CONTINUE READING
      READ(IREAT,102)ILOGT,KVADT(KOMP),MINET(KOMP),MAXET(KOMP),NTETB,
     * ITETA(KOMP)
C
c          print*, 'inabs    NTETB= ',ntetb
C
      ILOGTA(KOMP)=ILOGT
      IF(NTETB.GT.1)GO TO 702
  701 KOMPR=KOMPR+1
      GO TO 703
  702 READ(IREAT,103)(XTET(L),L=1,NTETB)
C
C        FINALLY THE ABSORPTION COEFFICIENT TABLE IS READ
  703 DO 704 K=1,NTETB
  704 READ(IREAT,104)(XKAP(JJ,K),JJ=1,NLATB)
C
C        WE TAKE THE LOGARITHMS BEFORE THE WAVELENGTH INTERPOLATION
C        IF ILOGL IS GREATER THAN ZERO.
      IF(ILOGL.LT.1)GO TO 712
  710 DO 711 K=1,NTETB
      DO 711 JJ=1,NLATB
      IF(XKAP(JJ,K).GT.0.)GO TO 711
C
C        A COEFFICIENT FOR WHICH THE LOGARITHM SHOULD BE TAKEN IS ZERO
      WRITE(IWRIT,207)JJ,K,XKAP(JJ,K),KOMP
      XKAP(JJ,K)=1.E-20
  711 XKAP(JJ,K)=log(XKAP(JJ,K))
  712 CONTINUE
C
C        PREPARATION OF THE T-INTERPOLATION IN SUBROUTINE TABS
C
C        WE FIND OUT WHETHER ISVIT(KOMP) CAN BE CHOSEN GREATER THAN ZERO. THIS
C        IS THE CASE IF THE T SCALE AND MAXET, MINET AND KVADT ARE IDENTICAL
C        WITH THOSE OF THE PREVIOUS COMPONENT. IF ISVIT IS GREATER THAN ZERO
C        THE TIME SPENT IN SUBR. TABS WILL BE DECREASED.
      ISVIT(KOMP)=0
      IF(NTETB.LE.1)GO TO 719
      DO 721 L=1,NTETB
      IF(XTET(L).NE.XTETP(L))GO TO 722
  721 CONTINUE
      IF(NTETB.NE.NTETBP)GO TO 722
      IF(MAXET(KOMP).NE.MAXETP) GO TO 722
      IF(MINET(KOMP).NE.MINETP) GO TO 722
      IF(KVADT(KOMP).NE.KVADTP)GO TO 722
      ISVIT(KOMP)=1
  722 CONTINUE
C
C        WE REMEMBER TEMPERATURES ETC. FOR NEXT COMPONENT
      DO723 L=1,NTETB
  723 XTETP(L)=XTET(L)
      NTETBP=NTETB
      MAXETP=MAXET(KOMP)
      MINETP=MINET(KOMP)
      KVADTP=KVADT(KOMP)
C
C        WE FIND THE INTERVALS IN THE T (TETA) SCALE
      TBOT(KOMP,1)=XTET(1)
      DELT(KOMP,1)=XTET(2)-XTET(1)
      NTM(KOMP,1)=1
      IDEL(KOMP)=1
      IF(NTETB.EQ.2)GO TO 719
C
      J=1
      LF=1
      DO714 L=3,NTETB
      DIFF=XTET(L)-XTET(L-1)
      IF(ABS(1.-DIFF/DELT(KOMP,J)).LT.1.E-4)GO TO 714
      J=J+1
      IF(J.GT.IELMAX)GO TO 715
      TBOT(KOMP,J)=XTET(L-1)
      DELT(KOMP,J)=DIFF
      NTM(KOMP,J-1)=LF
      LF=0
  714 LF=LF+1
      NTM(KOMP,J)=LF
      IDEL(KOMP)=J
      GO TO 719
C        TOO MANY DIFFERENT INTERVALS IN THE T-TABLE FOR THIS COMPONENT
  715 WRITE(IWRIT,203)KOMP,IELMAX
      WRITE(IWRIT,206)(XTET(L),L=1,NTETB)
      STOP 'INABS 1'
C
  719 NTAET(KOMP)=NTETB
C        ALL DATA NECESSARY BELOW FOR THIS COMPONENT ARE STORED ON UNIT
C        ISLASK
      WRITE(ISLASK)KVADL,MINEX,MAXEX,NLATB,ILOGL,IDISK,(IDISKV(J),J=1,
     *NLATB),(XLATB(J),J=1,NLATB),NTETB,ILOGT,(XTET(L),L=1,NTETB)
     *,((XKAP(JJ,K),JJ=1,NLATB),K=1,NTETB)
C
C        **** PRINT-OUT ****
      IF(IOUTS.LE.0)GO TO 7
      WRITE(IWRIT,211)KOMP,ABNAME(KOMP),SOURCE(KOMP),ABNAME(KOMP)
      WRITE(IWRIT,212)
      WRITE(IWRIT,213)XLATB(1),XLATB(NLATB)
      IF(MINEX.EQ.0)WRITE(IWRIT,214)
      IF(MINEX.GT.0)WRITE(IWRIT,215)
      IF(KVADL.EQ.0)WRITE(IWRIT,216)
      IF(KVADL.GT.0)WRITE(IWRIT,217)
      IF(ILOGL.EQ.0)WRITE(IWRIT,218)
      IF(ILOGL.GT.0)WRITE(IWRIT,219)
      IF(MAXEX.EQ.0)WRITE(IWRIT,220)
      IF(MAXEX.GT.0)WRITE(IWRIT,221)
      IF(IDISK.GT.0)WRITE(IWRIT,222)
      IF(NTETB-1)8,8,9
    8 WRITE(IWRIT,230)
      GO TO 7
    9 CONTINUE
      WRITE(IWRIT,223)
      WRITE(IWRIT,213)XTET(1),XTET(NTETB)
      IF(MINET(KOMP).EQ.0)WRITE(IWRIT,214)
      IF(MINET(KOMP).GT.0)WRITE(IWRIT,215)
      IF(KVADT(KOMP).EQ.0)WRITE(IWRIT,216)
      IF(KVADT(KOMP).GT.0)WRITE(IWRIT,217)
      IF(ILOGTA(KOMP).EQ.0)WRITE(IWRIT,218)
      IF(ILOGTA(KOMP).GT.0)WRITE(IWRIT,219)
      IF(MAXET(KOMP).EQ.0)WRITE(IWRIT,220)
      IF(MAXET(KOMP).GT.0)WRITE(IWRIT,221)
      IF(ISVIT(KOMP).GT.0)WRITE(IWRIT,224)
      WRITE(IWRIT,231)
    7 CONTINUE
  720 CONTINUE
C        END OF 'THE FIRST KOMP-LOOP'
C
      KOMPS=KOMPR+1
C
C
C        WE BUILD THE ABKOF ARRAY. INTERPOLATION IN WAVELENGTH.
C
C        LOOP OVER WAVELENGTH SETS ('THE ISET-LOOP')
      DO 70 ISET=1,NSET
      REWIND ISLASK
      NLB=NL(ISET)
      DO1 J=1,NLB
      XLA(J)=XL(J,ISET)
    1 XLA3(J)=XLA(J)**3
      INDEX=1
C
C        LOOP OVER COMPONENTS STARTS ('THE SECOND KOMP-LOOP')
      DO60 KOMP=1,NKOMP
      READ(ISLASK)KVADL,MINEX,MAXEX,NLATB,ILOGL,IDISK,(IDISKV(J),J=1,
     *NLATB),(XLATB(J),J=1,NLATB),NTETB,ILOGT,(XTET(L),L=1,NTETB)
     *,((XKAP(JJ,K),JJ=1,NLATB),K=1,NTETB)
      JI=1
      LAMBI=1
C
C        LOOP OVER WAVELENGTHS ('THE J-LOOP') STARTS
      DO60 J=1,NLB
C        SEARCHING IN WAVELENGTH
      IU=NLB*(KOMP-1)+J
      KOMPLA(IU)=INDEX
      DO24 JJ=1,NLATB
      IHELP=JJ
      IF(XLA(J)-XLATB(JJ))25,24,24
   24 LAMBI=JJ
   25 CONTINUE
      IF(IHELP-1)45,45,26
   26 IF(KVADL)33,33,27
   33 IF(NLATB-LAMBI-1)41,31,31
   27 IF(NLATB-LAMBI-1)41,28,29
C
C        QUADRATIC INTERPOLATION
   28 LAMBI=LAMBI-1
   29 CONTINUE
C        ARE DISCONTINUITIES PRESENT
      IF(IDISK.LE.0)GO TO 299
      IF(IDISKV(LAMBI+1).LE.0)GO TO 299
      IF(XLA(J).GT.XLATB(LAMBI+1))GO TO 292
  291 IF(IDISKV(LAMBI).GT.0)GO TO 31
      IF(LAMBI.EQ.1)GO TO 31
      LAMBI=LAMBI-1
      GO TO 299
  292 LAMBI=LAMBI+1
      IF(IDISKV(LAMBI+1).GT.0)GO TO 31
      IF(LAMBI+1.EQ.NLATB)GO TO 31
  299 CONTINUE
      DXX1=XLA(J)-XLATB(LAMBI)
      DXX2=XLA(J)-XLATB(LAMBI+1)
      DXX3=XLA(J)-XLATB(LAMBI+2)
      DX21=XLATB(LAMBI+1)-XLATB(LAMBI)
      DX32=XLATB(LAMBI+2)-XLATB(LAMBI+1)
      DX31=XLATB(LAMBI+2)-XLATB(LAMBI)
      A1=DXX2*DXX3/(DX21*DX31)
      A2=DXX1*DXX3/(DX21*DX32)
      A3=DXX1*DXX2/(DX31*DX32)
C
      DO30 K=1,NTETB
      ABKOF(INDEX)=A1*XKAP(LAMBI,K)-A2*XKAP(LAMBI+1,K)+A3*
     &XKAP(LAMBI+2,K)
   30 INDEX=INDEX+1
      GO TO 59
C
C        LINEAR INTER- AND EXTRAPOLATION
   31 A2=(XLA(J)-XLATB(LAMBI))/(XLATB(LAMBI+1)-XLATB(LAMBI))
      A1=1.-A2
      DO32 K=1,NTETB
      ABKOF(INDEX)=A1*XKAP(LAMBI,K)+A2*XKAP(LAMBI+1,K)
   32 INDEX=INDEX+1
      GO TO 59
C
C        TOO GREAT A WAVELENGTH - OUTSIDE THE TABLE
   41 IF(MAXEX)50,50,42
   42 LAMBI=LAMBI-1
      IF(NEXTL.GT.0)WRITE(IWRIT,201) KOMP,XLA(J)
      GO TO 31
C
C        TOO SMALL A WAVELENGTH - OUTSIDE THE TABLE
   45 IF(MINEX)50,50,46
   46 IF(NEXTL.GT.0)WRITE(IWRIT,201)KOMP,XLA(J)
      GO TO 31
C
C        ABS. COEFF. IS PUT = ZERO
   50 KOMPLA(IU)=0
      IF(NUTZL.GT.0)WRITE(IWRIT,202)KOMP,XLA(J)
      GO TO 60
C
   59 IF(ILOGL.LT.1)GO TO 592
      IF(ILOGT.GT.0)GO TO 60
C
C        LOGARITHMIC INTERPOLATION ONLY IN WAVELENGTH
      LIP=INDEX-NTETB
      LAP=INDEX-1
      DO 591 LL=LIP,LAP
  591 ABKOF(LL)=EXP(ABKOF(LL))
C
  592 CONTINUE
C
      IF(ILOGT.LE.0)GO TO 60
C        WE TAKE THE LOGARITHM BEFORE THE T INTERPOLATION IF ILOGT GT 0
      LIP=INDEX-NTETB
      LAP=INDEX-1
      DO593 LL=LIP,LAP
      IF(ABKOF(LL).GT.0.)GO TO 593
C
C        IMPOSSIBLE TO TAKE THE LOGARITHM OF A NEGATIVE COEFFICIENT
      LUS=LL-LIP+1
      WRITE(IWRIT,208)LL,ABKOF(LL),KOMP,J,ISET,LUS
      ABKOF(LL)=1.E-20
  593 ABKOF(LL)=log(ABKOF(LL))
   60 CONTINUE
C        END OF 'THE J-LOOP'
C        END OF 'THE SECOND KOMP-LOOP'
C
C        WRITE THE DATA OF THE SET ISET ON UNIT IRESET(ISET)
      NABKOF=INDEX-1
      NKOMPL=IU
      IREADP=IRESET(ISET)
      WRITE(IREADP)ISET,NLB,XLA,XLA3,NABKOF,ABKOF,NKOMPL,KOMPLA
C
      END FILE IREADP
      BACKSPACE IREADP
C
C        CHECK DIMENSION OF ABKOF
      IF(IOUTS.GT.0) WRITE(IWRIT,204)NABKOF,ISET
      IF(NABKOF.LE.NABDIM)GO TO 70
C        TOO SMALL DIMENSION FOR ABKOF
      WRITE(IWRIT,205)NABDIM
      STOP 'INABS 2'
   70 CONTINUE
C
C      WRITE(6,*) 'THE 10X20 WAVELENGTHS XL(20,10) (from INABS)'
C      DO 1002 IS=1,10
C1002  WRITE(6,1001) (XL(J,IS),J=1,20)
C1001  FORMAT (8F10.2)
C
C        END OF 'THE ISET-LOOP'
C
C
      DO 71 ISET=1,NSET
      IREADP=IRESET(ISET)
      REWIND IREADP
71    CONTINUE
      IF(IOUTS.LE.0) GOTO 74
C
C        **** PRINT-OUT ****
C        ON WAVELENGTH SETS AND FILES
      WRITE(IWRIT,225)IREAT,ISLASK
      WRITE(IWRIT,226)
      DO73 M=1,NSET
      NP=NL(M)
      WRITE(IWRIT,227)M,IRESET(M)
      WRITE(IWRIT,228)(XL(J,M),J=1,NP)
   73 CONTINUE
      WRITE(IWRIT,232)
   74 CONTINUE
  101 FORMAT(8X,I2,5(9X,I1))
  102 FORMAT(4(9X,I1),8X,I2,9X,I1)
  103 FORMAT(6F10.0)
  104 FORMAT(6E10.3)
  105 FORMAT(2A8)
  201 FORMAT(' EXTRAPOLATION FOR COMPONENT',I5,5X,'AND WAVELENGTH=',F10.
     *3,5X,'***INABS***')
  202 FORMAT(
     *' ABS.COEFF. PUT=0 AT WAVELENGTH-INTER/EXTRAPOLATION FOR COMP. ',
     *I5,5X,'AND WAVELENGTH=',F10.3,5X,'***INABS***')
  203 FORMAT(
     *' TOO MANY DIFFERENT INTERVALS IN THE T-(TETA-)TABLE FOR COMP. ',
     *I5,5X,'MAX IS',I5,5X,'***INABS***')
  204 FORMAT(' NECESSARY DIMENSION FOR ABKOF=',I5,5X,'IN SET',I5,5X,
     *'***INABS***')
  205 FORMAT(' DIMENSION ALLOWED =',I5,5X,'TOO SMALL     ***INABS***')
  206 FORMAT(6H XTET=,10E12.4)
  207 FORMAT(' XKAP(',I2,',',I2,')=',E12.5,
     *' PUT = 1.E-77 BEFORE LOG HAS BEEN TAKEN.    ***INABS***')
  208 FORMAT(' ABKOF(',I4,')=',E12.5,' FOR COMPONENT ',I2,' WAVELENGTH',
     *I3,' SET ',I2,' XTET NR ',I2,' ABKOF PUT=1.E-77   ***INABS***')
  211 FORMAT('0************  COMPONENT NO',I3,', ',A8,' SOURCE  ',A8,
     *'****************',5X,A8)
  212 FORMAT('0     I N T E R P O L A T I O N  I N  W A V E L E N G T H'
     *)
  213 FORMAT(1H ,15X,F10.2,25X,F10.2)
  214 FORMAT(1H+,'   KAPPA=0   - ')
  215 FORMAT(1H+,' LIN. EXTRAP.- ')
  216 FORMAT(1H+,25X,' -LIN. INT. ')
  217 FORMAT(1H+,25X,' -QUAD. INT.')
  218 FORMAT(1H+,37X,'  KAPPA    - ')
  219 FORMAT(1H+,37X,' LOG(KAPPA)- ')
  220 FORMAT(1H+,60X,' -   KAPPA=0')
  221 FORMAT(1H+,60X,' -LIN. EXTRAP.')
  222 FORMAT(1H0,' DISCONTINUITIES PRESENT.')
  223 FORMAT('0     I N T E R P O L A T I O N  I N  T  ( T E T A )')
  224 FORMAT('0 T SCALE ETC. IDENTICAL WITH PRECEEDING COMPONENT')
  225 FORMAT(1H0,'F I L E S  U S E D  B Y  T H E  A B S - B L O C K'//
     * '  INITIAL FILE ',I3,', PRELIMINARY FILE',I3)
  226 FORMAT(1H0,'SET           WAVELENGTHS',81X,'FILE')
  227 FORMAT (1H ,I2,105X,I2)
  228 FORMAT(1H ,5X,10F10.2)
  229 FORMAT(1H1,'D A T A  F R O M  S U B R O U T I N E  I N A B S')
  230 FORMAT(1H0,'     N O  T-  ( T E T A - ) D E P E N D E N C E')
  231 FORMAT(1H0)
  232 FORMAT(1H1)
      RETURN
      END
C
      SUBROUTINE INITAB (IOUTS)
      implicit real*16 (a-h,o-z)
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
      implicit real*16 (a-h,o-z)
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
      implicit real*16 (a-h,o-z)
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
     *XIONG(16,5),EEV,ENAMN(ndp),SUMH(ndp),XKBOL,
     *SUMM(ndp),NJ(16),IEL(16),NEL
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
      implicit real*16 (a-h,o-z)
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
      implicit real*16 (a-h,o-z)
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
     *XIONG(16,5),EEV,ENAMN(ndp),SUMH(ndp),XKBOL,
     *SUMM(ndp),NJ(16),IEL(16),NEL
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
     &  xmettryck(ndp,maxmet),xiontryck(ndp,maxmet),
     *  partp(ndp,0:maxmol)
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
      SUBROUTINE KAP5(T,PE,ABSK)
      implicit real*16 (a-h,o-z)
C
      COMMON /CXLSET/XL(20,10),NSET,NL(10)
C
C COMPUTE KAPPA(5000.). 73.10.17 *NORD*.
      CALL ABSKO(1,1,T,PE,1,NL(1)+1,ABSK,SPRID)
      RETURN
      END
C
       FUNCTION LENSTR(STRING)
      implicit real*16 (a-h,o-z)
*
* Returns the length of a string not counting trailing blanks
*
       CHARACTER*(*)  STRING
*
       DO 10 I = LEN(STRING), 1, -1
         IF(STRING(I:I) .NE. ' ') THEN
           LENSTR = I
           RETURN
         ENDIF
   10  CONTINUE
       LENSTR = 0
       RETURN
       END
         
        

C
      SUBROUTINE LISTMO(MO,IARCH,ISPH)
      implicit real*16 (a-h,o-z)
C        THIS ROUTINE PRINTS A NUMBER (MO) OF MODELS, WHICH ARE STORED ON
C        FORTRAN UNIT IARCH.
C
      include 'parameter.inc'
C
      DIMENSION TKORRM(NDP),FCORR(NDP),TAU(NDP),TAUS(NDP),
     *PE(NDP),PG(NDP),PRAD(NDP),PTURB(NDP),XKAPR(NDP),RO(NDP),
     *CP(NDP),CV(NDP),AGRAD(NDP),Q(NDP),U(NDP),V(NDP),ANCONV(NDP),
     *PRESMO(33,NDP),FCONV(NDP),RR(NDP),Z(NDP),EMU(NDP),HNIC(NDP)
     *,NJ(16),XLR(20),IEL(16),PROV(20,20),
     *ABSKA(20),SPRIDA(20),XLB(500),PEP(16)
     *,ABNAME(30),SOURCE(30),ABSKTR(NDP),SPRTR(NDP)
C     *,ABNAME(30),SOURCE(30),ABSKTR(NDP),SPRTR(NDP),DUMMY(11*NDP+2)
      DIMENSION W(500),UW(12),BW(21),VW(25),SUMPMOL(NDP)
      dimension xdp(8),kdp(8)
C      REAL*8 ROSSO,PTAUO
      COMMON /CVAAGL/XLB,W,NLAM
      COMMON /CARC2/T(NDP),FC(NDP),FLUXME(NWL),TAU5(NDP),INORD
      CHARACTER*10 DAG,KLOCK
      CHARACTER*8 ABNAME,SOURCE
      DIMENSION WAVFLX(10)
      COMMON /COPINF/ SUMOP(maxosmol,NDP),SUMKAP(maxosmol,NDP)
      COMMON /UTPUT/IREAD,IWRIT
C      COMMON /TIO/PTIO(NDP),DUMDUM(2)
      COMMON /KETIO/ptio(ndp),roke(ndp),poxg1(ndp)
      COMMON/CI5/abmarcs(17,ndp),ANJON(17,5),H(5),PART(17,5),
     * DXI,F1,F2,F3,F4,F5,XKHM,XMH,XMY(ndp)
      COMMON /CLEVETAT/GEFF(NDP),PPRG(NDP),AMLOSS
      COMMON /CLEVPRINT/ PRJ2(NDP),masslinf
      CHARACTER MOLNAME*4,OSFIL*60,SAMPLING*3
      COMMON /COSLIST/ WNB(25),WNSTEP(25),WNEND,INTVOS
      COMMON/COS/WNOS(NWL),CONOS(NDP,NWL),WLOS(NWL),WLSTEP(NWL)
     *    ,KOS_STEP,NWTOT,NOSMOL,NEWOSATOM,NEWOSATOMLIST
     *    ,nchrom,OSFIL(30),MOLNAME(30),SAMPLING
      COMMON/COSWR/osresl,losresl,listwn
C      COMMON/COPPR/oppr(15,3,120,3),jvxmax,itxmax  !15mol,10dpt,100wn
C      COMMON/COPPRR/xconop(120,10),xlineop(120,10)    !100wn,10dpt
      COMMON /Cspec/spec(nwl,3),ispec
      COMMON /CLIST/NLTE
      COMMON /MASSE/RELM
      COMMON /CLIN/lin_cia
      COMMON /CNEWC3 /NEWC3
      COMMON /CG/GRAV,KONSG
      COMMON /CSTYR/MIHAL,NOCONV /CXMAX/XMAX /CTAUM/TAUM
      COMMON /MIXC/PALFA,PBETA,PNY,PY /CVFIX/VFIX                          
      COMMON /CPOLY/FACPLY,MOLTSUJI
      COMMON /CROSSOS/ ROSSO(NDP),PTAUO(NDP)
      COMMON /CINDIAM/TDIAM1,TDIAM2,FDIAM1,FDIAM2,
     *    TC2H21,TC2H22,FC2H21,FC2H22
      DATA UW/0.145,0.436,0.910,1.385,1.843,2.126,2.305,2.241,1.270,
     *0.360,0.128,0.028/,BW/0.003,0.026,0.179,0.612,1.903,2.615,2.912,
     *3.005,2.990,2.876,2.681,2.388,2.058,1.725,1.416,1.135,0.840,0.568,
     *0.318,0.126,0.019/,VW/0.006,0.077,0.434,1.455,2.207,2.703,2.872,
     *2.738,2.505,2.219,1.890,1.567,1.233,0.918,0.680,0.474,0.312,0.200,
     *0.132,0.096,0.069,0.053,0.037,0.022,0.012/
      COMMON /FULLEQUILIBRIUM/ PARTRYCK(NDP,MAXMOL),
     &  XMETTRYCK(NDP,MAXMET),XIONTRYCK(NDP,MAXMET),
     * PARTP(NDP,0:MAXMOL),
     & PARTPP(NDP,0:MAXMOL)
      common /tsuji/ nattsuji,nmotsuji,parptsuji(500),abtsuji(17,ndp)
      COMMON /CMTEST/PGM1(NDP),PGM2(NDP),PEM1(NDP),PEM2(NDP)
     *     ,TM1(NDP),TM2(NDP),PGOS(NDP)
      COMMON/COPsum/ SSUM(NDP),XSUM(NDP),CONSUM(NDP)
      COMMON/CI4/ TMOLIM, IELEM(16), ION(16,5), MOLH, JUMP
      COMMON /ROSSC/CXKAPR(NDP),CROSS(NDP)
      COMMON /COSEXP/ LOPS,NOPS
      DATA A,B/.34785485,.65214515/
      COMMON /CCIATEST/ CIATEST(44,NDP)
      COMMON /PJONINF/ P_MOL(NDP), P_NEU_HCNO(NDP), P_ION_HCNO(NDP),
     & P_NEU_HE(NDP),P_ION_HE(NDP), P_NON_HHECNO(NDP), PG_JON(NDP), 
     & HN_JON(NDP), RO_JON(NDP), P6_JON(NDP)
      COMMON /CKMOL/KMOL(MAXOSMOL)   !connects OS-molecule with presmo-index
      data xdp / -5., -4., -3., -2., -1., 0., 1., 2./
      character*5 name_mol, name_listmo
      COMMON /CMOLNAME/NAME_MOL(maxmol),NAME_LISTMO(maxmol)
      COMMON /STATEC/PPR(NDP),PPT(NDP),PP(NDP),GG(NDP),ZZ(NDP),DD(NDP),
     &VV(NDP),FFC(NDP),PPE(NDP),TT(NDP),TAULN(NDP),RO_ST(NDP),NTAU,ITER

      common /cgem/pres_gem(ndp,nspec)
      common /cgemnames/natms_gem,nions_gem,nspec_gem,name_gem(nspec)
C atms,ions,spec ~ highest index of neutral atoms, ions, species total
      character name_gem*8
      common /cabink/abink(ndp,nspec)
      common /cprespp/prespp(ndp,nspec)           !gem computed pp
      dimension phe(ndp), trpe(ndp), trphe(ndp)
      dimension printpp(ndp,nspec), printname(nspec)
      character printname*8
      dimension listpp(nspec)
      dimension ptot(ndp),pp_sum(ndp)
      common /ctotabk/totabk(ndp,natms)
      dimension pe_gem(ndp)
*
      FLUMAG(I)=-2.5*log10(FLUXME(I))-STMAGN
*
C_ursa      listpp(1:3)=1               ! e-,H,He
C_ursa      listpp(4:6)=0
C_ursa      listpp(7:9)=1               ! C,N,O
C_ursa      listpp(10:11)=0
C_ursa      listpp(12:15)=1             ! Na,Mg,Al,Si
C_ursa      listpp(16:19)=0
C_ursa      listpp(20:22)=1             ! K,Ca,Ti
C_ursa      listpp(23:25)=0
C_ursa      listpp(26)=1                ! Fe
C_ursa      listpp(27:49)=0
C_ursa      listpp(50:52)=1             ! H+,H-,He+
C_ursa      listpp(53:66)=0
C_ursa      listpp(67)=1                ! Na+
C_ursa      listpp(68)=0
C_ursa      listpp(69:70)=1             ! Mg+,Al+
C_ursa      listpp(71:80)=0
C_ursa      listpp(81)=1                ! K+
C_ursa      listpp(82)=0
C_ursa      listpp(83:84)=1             ! Ca+,Ti+
C_ursa      listpp(85:87)=0
C_ursa      listpp(88)=1                ! Cr+
C_ursa      listpp(89:90)=0
C_ursa      listpp(91)=1                ! Fe+
C_ursa      listpp(92:94)=0
C_ursa      listpp(95)=1                ! Ni+
C_ursa      listpp(96:127)=0
C_ursa      listpp(128)=1               ! H2
C_ursa      listpp(129:131)=0
C_ursa      listpp(132:134)=1           ! C2,N2,O2
C_ursa      listpp(135:151)=0
C_ursa      listpp(152:155)=1           ! CO,CH,CN,OH
C_ursa      listpp(156:159)=0
C_ursa      listpp(160)=1               ! CS
C_ursa      listpp(161:164)=0
C_ursa      listpp(165)=1               ! HS
C_ursa      listpp(166:169)=0
C_ursa      listpp(170)=1               ! MgH
C_ursa      listpp(171)=0
C_ursa      listpp(172)=1               ! SiH
C_ursa      listpp(173:178)=0
C_ursa      listpp(179)=1               ! SiC
C_ursa      listpp(180:195)=0
C_ursa      listpp(196)=1               ! SiO
C_ursa      listpp(197:200)=0
C_ursa      listpp(201)=1               ! TiO
C_ursa      listpp(202:326)=0
C_ursa      listpp(327:329)=1           ! H2+,H2-,C2-
C_ursa      listpp(330:340)=0
C_ursa      listpp(341)=1               ! SiH+
C_ursa      listpp(342:344)=0
C_ursa      listpp(345)=1               ! NaO-
C_ursa      listpp(346:363)=0
C_ursa      listpp(364)=1               ! C3
C_ursa      listpp(365:366)=0
C_ursa      listpp(367)=1               ! C2H
C_ursa      listpp(368:370)=0
C_ursa      listpp(371:373)=1           ! CH2,HCN,CHO
C_ursa      listpp(374:376)=0
C_ursa      listpp(377:378)=1           ! CO2,SiC2
C_ursa      listpp(379:381)=0
C_ursa      listpp(382)=1               ! H2O
C_ursa      listpp(383:574)=0
C_ursa      listpp(575)=1               ! C2H2
C_ursa      listpp(576:578)=0
C_ursa      listpp(579:580)=1           ! CH3,CH4
C_ursa      listpp(581)=0
C_ursa      listpp(582)=1               ! NH3
C_ursa      listpp(583:845)=0
C_ursa      listpp(846)=1               ! SiH4
C_ursa      listpp(847:892)=0

      if (ispec.eq.1) then
       open (unit=29,file='spectrum.dat',status='unknown')
       do 290 i=1,nwtot
       write(29,295) (spec(i,k),k=1,3)
290    continue
       close(29)
295    format(f10.2,1p2e12.4)
      end if

      CALL ROSSOS

      IREAD=5
      IWRIT=7
C
      REWIND IARCH
      DO 1 IMO=1,MO
C      READ(IARCH) DUMMY
      READ(IARCH) INORD,DAG,KLOCK
      READ(IARCH) TEFF,FLUX,G,PALFA,PNY,PY,PBETA,ILINE,ISTRAL,MIHAL,
     &            IDRAB1,IDRAB2,IDRAB3,IDRAB4,IDRAB5,IDRAB6,
     &            ITMAX,NEL,(abmarcs(I,1),I=1,NEL)
      WRITE(7,219)
      GLOG=log10(G)
      WRITE(7,300) TEFF,GLOG,IDRAB1,IDRAB2,IDRAB3,IDRAB4,IDRAB5,IDRAB6
      WRITE(7,283)
      FNORD=0.1*INORD
      WRITE(7,282) FNORD,DAG,KLOCK
      WRITE(7,284)
      WRITE(7,200)
C        CONVERT TO 'PHYSICAL FLUX'
      FLUX=3.14159*FLUX
      corat = abmarcs(3,1)/abmarcs(5,1)
      if (corat.gt.99. .or. corat.lt.0.1) then
      WRITE(7,2017) TEFF,FLUX,G,GLOG,
     &  abmarcs(5,1)/8.51138e-4,abmarcs(3,1)/abmarcs(5,1),PALFA,PNY,PY
      else
      WRITE(7,201) TEFF,FLUX,G,GLOG,
     &  abmarcs(5,1)/8.51138e-4,abmarcs(3,1)/abmarcs(5,1),PALFA,PNY,PY
      end if
      WRITE(7,2011) NOCONV 
      IF(PBETA.LE.0.) WRITE(7,256)
      IF(PBETA.GT.0.1) WRITE(7,257) PBETA
      IF(ISTRAL.LT.1) WRITE(7,250) ISTRAL
      IF(ISTRAL.GE.1) WRITE(7,251) ISTRAL
      IF(ILINE.LT.1) WRITE(7,252) ILINE
      IF(ILINE.GE.1) WRITE(7,253) ILINE
      WRITE(7,2531) RELM
      IF (isph.ne.1) THEN 
      WRITE(7,2539)
      ENDIF
      WRITE(7,2532) NOSMOL
2532  format(' Following ',i3,' molecules are included in opacity:')
      WRITE(7,2533) (MOLNAME(I),I=1,NOSMOL)
2533  FORMAT(18(2X,A4))
      write(7,*) ' from the following inputfiles:'
      do 2166 i=1,nosmol
2166  write(7,2015) OSFIL(i)
      if (fdiam1.ne.0..or.fdiam2.ne.0..or.fc2h21
     &    .ne.1..or.fc2h22.ne.1)
     &write(7,2541) tdiam1,tdiam2,100.*fdiam1,100.*fdiam2,
     *    tc2h21,tc2h22,100.*fc2h21,100.*fc2h22

      IF (losresl.eq.1) THEN
        write(7,*) ' We did a resolution based set of os wavenumbers'
        write(6,2551) nwtot,osresl/dfloat(kos_step)
     &  ,wnos(1),wnos(nwtot),1.e4/wnos(nwtot),1.e4/wnos(1)
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
2551   format(' Total ',i6,' OS wavenumbers for Marcs radiative transf.'
     &  ,/' Approximate resolution was',f8.0
     &  ,/' OS interval:',f7.1,'-',f9.1,' cm^-1 (=',f5.2,'-',f5.1,'mu)')
      ELSE IF (listwn) THEN
        write(7,*) ' the OS wavenumbers were from an input list'
        write(6,2552) nwtot
     &  ,wnos(1),wnos(nwtot),1.e4/wnos(nwtot),1.e4/wnos(1)
2552   format(' Total ',i6,' OS wavenumbers for Marcs radiative transf.'
     &  ,/' OS interval:',f7.1,'-',f9.1,' cm^-1 (=',f5.2,'-',f5.1,'mu)')
      ELSE
        WRITE(7,2534) INTVOS
        WRITE(7,2535) WNB(1),WNEND
        WRITE(7,2536) 1.e8/WNEND,1.e4/WNB(1)
        WRITE(7,2537)
        write(7,2538) (WNB(I),I=1,INTVOS)
        write(7,2538) (WNSTEP(I),I=1,INTVOS)
      END IF

      if (lin_cia.eq.1) write(7,2529)
      WRITE(7,2012) XMAX, TAUM, facply, moltsuji
      if (jump.eq.2.or.4) then
       WRITE(7,2018)
      else if (jump.eq.3) then
       WRITE(7,2019)
      else if (jump.eq.1) then
       WRITE(7,2013)
      if (newc3.eq.1) write(7,*) ' For C3, Irwins Kp value was used'
      else if (jump.eq.0) then
       WRITE(7,2014)
       write(7,*) 
     *  ' Irwins Kp value was not used (only impl.for Tsuji eq.)'
      end if
      if (lops.ne.0  .or.  nops.ne.1)
     * WRITE(7,2016) lops,nops
2016  format (' the os-absorption was shifted',i3,' os-steps,',
     *    ' and only each',i4,' os-value was used (rest==0)')
      WRITE(7,254) MIHAL,IDRAB1,IDRAB2,IDRAB3,IDRAB4,IDRAB5,IDRAB6
      WRITE(7,202)
      DO 2 I=1,NEL
    2 abmarcs(I,1:ntau)=log10(abmarcs(I,1:ntau))+12.
      abmarcs(17,1:ntau)=log10(abmarcs(17,1:ntau))+12.
      WRITE(7,203) (abmarcs(I,1),I=1,NEL+1)
      WRITE(7,255) ITMAX
      READ(IARCH)JTAU,NCORE,DIFLOG,TAUM,RADIUS,(RR(K),K=1,JTAU)
      READ(IARCH)JTAU,(TKORRM(I),I=1,JTAU),(FCORR(K),K=1,JTAU)
      NTPO=0
      DO 3 K=1,JTAU
      READ(IARCH) KR,TAU(K),TAUS(K),Z(K),T(K),PE(K),PG(K),PRAD(K),
     &            PTURB(K),XKAPR(K),RO(K),EMU(K),CP(K),CV(K),
     &            AGRAD(K),Q(K),U(K),V(K),ANCONV(K),HNIC(K),NMOL,
     &            (PRESMO(J,K),J=1,NMOL)
      TAUK=log10(TAU(K))+10.01
      KTAU=TAUK
      IF(ABS(TAUK-KTAU).GT.0.02) GO TO 31
      IF(KTAU.EQ.10) K0=K
      NTPO=NTPO+1
   31 CONTINUE
    3 CONTINUE
      WRITE(7,204)
      WRITE(7,300) TEFF,GLOG,IDRAB1,IDRAB2,IDRAB3,IDRAB4,IDRAB5,IDRAB6
      WRITE(7,205)
      DO 4 I=1,JTAU
      FCONV(I)=ANCONV(I)*FLUX
      WRITE(7,206) I,TAU(I),T(I),TKORRM(I),FCONV(I),FCORR(I),I
    4 CONTINUE
C
      IF (NLTE.EQ.0) GO TO 4000
C*
C* 90-05-13 START OF MODIFICATIONS (MATS CARLSSON)
C* PRINT MULTI ATMOSPHERIC FILE: ATMOS.MULTI
C*
      OPEN(33,FILE='ATMOS.MULTI',STATUS='NEW',CARRIAGE CONTROL='LIST')
      WRITE(33,400) TEFF
  400 FORMAT(' MARCS MODEL ATMOSPHERE, TEFF=',F10.2/' TAU(5000) SCALE')
      WRITE(33,410) G
  410 FORMAT('*'/'* LG G'/F6.2)
      WRITE(33,420) JTAU
  420 FORMAT('*'/'* NDEP'/I3)
      WRITE(33,430)
  430 FORMAT('*'/'*LG TAU(5000)    TEMPERATURE        NE         V',
     * '              VTURB')
      DO 450 I=1,JTAU
        IF(TAUS(I).GT.0.0) THEN
          TAULG=LOG10(TAUS(I))
        ELSE
          TAULG=2.*LOG10(TAUS(I+1))-LOG10(TAUS(I+2))
        ENDIF
        WRITE(33,440) TAULG,T(I),PE(I)/T(I)/1.380662E-16,0.,2.
  440   FORMAT(1P,5E14.6)
  450 CONTINUE
      CLOSE(33)
C*
C* 90-05-13 END OF MODIFICATIONS
C*
4000  CONTINUE

      WRITE(7,207)
      WRITE(7,300) TEFF,GLOG,IDRAB1,IDRAB2,IDRAB3,IDRAB4,IDRAB5,IDRAB6
      WRITE(7,208)
      Z0=Z(1)
      DO 5 I=1,JTAU
      Z(I)=Z(I)-Z0
      WRITE(7,209) I,TAU(I),TAUS(I),Z(I),T(I),PE(I),PG(I),PRAD(I),
     &             PTURB(I),XKAPR(I),I
      IF (T(I).GT.TEFF) then
         iint = i
         if (t(i)-teff .gt. teff-t(i-1)) iint = i-1
      END IF
    5 CONTINUE

C      masslinf = 1                    !now (sept.2006) in namelist outlist 
      if (masslinf.eq.0) go to 4001
      WRITE(7,2071)
      WRITE(7,300) TEFF,GLOG,IDRAB1,IDRAB2,IDRAB3,IDRAB4,IDRAB5,IDRAB6
      WRITE(7,2081)
      DO 51 I=1,JTAU
      IF (I.LE.JTAU-1) THEN
      ROMIDT = (RO(I)+RO(I+1) )/2.
      DPGDZ = ( PG(I+1)-PG(I) ) / ( Z(I)-Z(I+1) ) / ROMIDT
      DPRDZ = ( PRAD(I+1)-PRAD(I) ) / ( Z(I)-Z(I+1) ) / ROMIDT
      END IF
      RO(I) = MAX (RO(I),1.D-99)
      RI = SQRT(RELM/G) * 1.152E13 + Z(JTAU) - Z(I)
      RILOG = LOG10(RI)
      VLOG = log10(AMLOSS) + 24.700 - 2.*RILOG - LOG10(RO(I))
      VLOSS = 10.**VLOG *1.e-5    !velocity in km/s
      ZINT = Z(IINT) - Z(I)
      WRITE(7,2091) I,TAU(I),LOG10(TAU(I)),RO(I),ZINT,T(I),PPRG(I),
     *         DPRDZ,VLOSS,ROSSO(I),PTAUO(I),cross(I)
   51 CONTINUE
4001  CONTINUE

      WRITE(7,210)
      WRITE(7,300) TEFF,GLOG,IDRAB1,IDRAB2,IDRAB3,IDRAB4,IDRAB5,IDRAB6
      WRITE(7,211)
      DO 6 I=1,JTAU
      WRITE(7,212) I,TAU(I),RO(I),EMU(I),CP(I),CV(I),AGRAD(I),Q(I),U(I),
     &             V(I),ANCONV(I),I
    6 CONTINUE
C

      DO 772 I=1,JTAU
        HNIC(I) = max(1.d-99,HNIC(I))
        HNIC(I)=LOG10(HNIC(I))
           SUMPMOL(I) = 0.
        DO 771 J=1,33
           PRESMO(J,I)=max(1.d-99,PRESMO(J,I))
           if(j.eq.1.or.j.eq.17.or.j.eq.32.or.j.eq.33) then
                 PRESMO(J,I)=LOG10(PRESMO(J,I))
                 go to 771   !(P_He is in presmo(17), 1 is H-, 32,33 extr)
           end if
           SUMPMOL(I) = SUMPMOL(I) + PRESMO(J,I)
           PRESMO(J,I)=LOG10(PRESMO(J,I))
771      CONTINUE
772   CONTINUE
C Molecular partial pressures from the Marcs equilibrium
C UGJ 14/11/98: Test with partial pressures from the two
C routines show them almost identical for giants (for very cool
C giants and dwarfs they are probably a bit different from 
C one another). You can write both sets by deleting the first
C "if(jump.ne.1)" loop. but be aware that the spectrum program
C may have difficult finding P(ZrO), P6 and other part.pressures then.
C     
      IF (JUMP.GT.0) GO TO 2101


      WRITE(7,213)
      WRITE(7,300) TEFF,GLOG,IDRAB1,IDRAB2,IDRAB3,IDRAB4,IDRAB5,IDRAB6

      WRITE(7,214)

      DO 7 I=1,JTAU
        WRITE(7,215) I,HNIC(I),(PRESMO(JJ,I),JJ=1,13),PRESMO(31,I),I
7     CONTINUE


      WRITE(7,213)
      WRITE(7,300) TEFF,GLOG,IDRAB1,IDRAB2,IDRAB3,IDRAB4,IDRAB5,IDRAB6

      WRITE(7,216)
      DO 8 I=1,JTAU
      WRITE(7,217) I,(PRESMO(J,I),J=14,16),(PRESMO(J,I),J=18,30),I
8     CONTINUE


      GO TO 2109

2101  CONTINUE     ! go here if jump>0 (i.e. not old marcs chem. equilibrium.

      IF (JUMP.EQ.3) GO TO 2102

C here: JF's/Tsuji's molecular partial pressures from JANAF polynomial fits:

C for Tsuji's eq names are not read in (MOL(J) as real*8...)

      IF (jump.eq.1) THEN

       open(unit=42,file='data/name_mol.dat'
     &    ,status='old')
       do 442 i=1,207
       read(42,443,end=449) j,name_mol(i)
442    continue
449    continue
443    format(i5,2x,a5)
       close(42)

       do 10 j=4,16
        name_listmo(j) = name_mol(j-1)
  10    continue
       name_listmo(1) = name_mol(1)
       name_listmo(2) = name_mol(2)
       name_listmo(3) = 'H2+   '
       name_listmo(17) = '     '
       name_listmo(20) = 'C3H  '
       name_listmo(36) = 'CaH  '
       name_listmo(37) = 'LaO  '
       name_listmo(69) = 'TiS  '

      do 720 j=18,19
        name_listmo(j) = name_mol(j-2)
720   continue
      do 721 j=21,35
        name_listmo(j) = name_mol(j-3)
721   continue
      do 722 j=38,68
        name_listmo(j) = name_mol(j-5)
722   continue
      do 724 j=70,213
        name_listmo(j) = name_mol(j-6)
724   continue

      END IF

      IF (jump.eq.4) THEN

       open(unit=42,file='data/name_mol.dat'
     &    ,status='old')
       do 725 i=1,207
       read(42,726,end=727) j,name_mol(i)
725    continue
727    continue
726    format(i5,2x,a5)
       close(42)

       do 15 j=4,16
        name_listmo(j) = name_mol(j-1)
  15    continue
       name_listmo(1) = name_mol(1)
       name_listmo(2) = name_mol(2)
       name_listmo(3) = 'H2+   '
       name_listmo(17) = '     '
       name_listmo(20) = 'C3H  '
       name_listmo(36) = 'CaH  '
       name_listmo(37) = 'LaO  '
       name_listmo(69) = 'TiS  '

      do 728 j=18,19
        name_listmo(j) = name_mol(j-2)
728   continue
      do 729 j=21,35
        name_listmo(j) = name_mol(j-3)
729   continue
      do 730 j=38,68
        name_listmo(j) = name_mol(j-5)
730   continue
      do 731 j=70,213
        name_listmo(j) = name_mol(j-6)
731   continue

      END IF

      IF (JUMP.LE.3) then
! H, H-, H2, H2+, H2O, OH, CH, CO, CN, C2, N2, O2, NO, NH, TiO
! 0   1   2   3     4   5   6   7   8   9  10  11  12  13   31
        WRITE(7,2131)
        WRITE(7,300) TEFF,GLOG,IDRAB1,IDRAB2,IDRAB3,IDRAB4,IDRAB5,IDRAB6
        write(7,2145)(name_listmo(i),i=1,13),name_listmo(31)
2145    FORMAT(' K',5X,'P(H)',14(3x,a5))
        do i=1,jtau
!                partryck(I,3)=1.0E-20      !pressure of H2+ not calculated
                write(7,215) i,(log10(xmettryck(i,1))),
     &         (log10(partryck(i,jj)),jj=1,13),(log10(partryck(i,31))),i
        end do

! C2H2, HCN, C2H, HS, SiH, C3H, C3, CS, SiC, SiC2, NS, SiN, SiO, SO, S2, SiS
!  14    15  16  18 19  20  21 22 23  24   25  26  27 28 29  30
        write(7,2131)
        write(7,300) teff,glog,idrab1,idrab2,idrab3,idrab4,idrab5,idrab6
        write(7,2146)(name_listmo(i),i=14,16),(name_listmo(i),i=18,30)
2146    format(' K ',16(3x,a5))
        do 811 i=1,jtau
                 write(7,217) i,(log10(partryck(i,j)),j=14,16),
     *          (log10(partryck(i,j)),j=18,30),i
811     continue

! VO, ZrO, MgH, CaH, LaO, CH4, CH2, CH3, NH3, CO2, Si2C, SiO2, H2S, CaOH, CHNO, SiF2
C 32  33  34  36  37 39  40  41  43  46  57   58   59   67  107 135
        WRITE(7,2131)
        WRITE(7,300) TEFF,GLOG,IDRAB1,IDRAB2,IDRAB3,IDRAB4,IDRAB5,IDRAB6
        WRITE(7,2147)(name_listmo(i),i=32,34)
     &          ,(name_listmo(i),i=36,37)
     &          ,(name_listmo(i),i=39,41)
     &          ,name_listmo(43)
     &          ,name_listmo(46)
     &          ,(name_listmo(i),i=57,59)
     &          ,name_listmo(67)
     &          ,name_listmo(107)
     &          ,name_listmo(135)
2147    FORMAT(' K ',16(3x,a5))
        DO 812 I=1,JTAU
                 WRITE(7,217) I,(log10(partryck(I,J)),J=32,34),
     *          (log10(partryck(I,J)),J=36,37),
     *          (log10(partryck(I,J)),J=39,41),
     *          (log10(partryck(I,J)),J=43,43),
     *          (log10(partryck(I,J)),J=46,46),
     *          (log10(partryck(I,J)),J=57,59),
     *          (log10(partryck(I,J)),J=67,67),
     *          (log10(partryck(I,J)),J=107,107),
     *          (log10(partryck(I,J)),J=135,135),I
812     CONTINUE

! TiH, CaH, FeH, CrH, LiH
! 214, 215, 216, 217, 71
        WRITE(7,2131)
        WRITE(7,300) TEFF,GLOG,IDRAB1,IDRAB2,IDRAB3,IDRAB4,IDRAB5,IDRAB6
        WRITE(7,2148)(name_listmo(i),i=214,217),name_listmo(71)
        DO I=1,JTAU
                WRITE(7,'(I3,5F8.3,2X,I3)')I,(log10(partryck(I,J)),
     *           J=214,217),log10(partryck(I,71)),I
        end do

C  AlO SiH4 SO2 YO  C5H C4H  C4  C5 C6H C4H6 C6H2 C6H4 C10H7 HC5N HC7N HC9N
C  136 144  162 258 302 306 315 316 317 318  319  320  321   322   323  324
!      WRITE(7,2131)
!      WRITE(7,300) TEFF,GLOG,IDRAB1,IDRAB2,IDRAB3,IDRAB4,IDRAB5,IDRAB6
!      write(7,2171)
!      DO 813 I=1,JTAU
!         WRITE(7,217) I,(log10(partryck(I,J)),J=136,136),
!     *   (log10(partryck(I,J)),J=144,144),
!     *   (log10(partryck(I,J)),J=162,162),
!     *   (log10(partryck(I,J)),J=258,258),
!     *   (log10(partryck(I,J)),J=302,302),
!     *   (log10(partryck(I,J)),J=306,306),     
!     *   (log10(partryck(I,J)),J=315,324),I
!813   CONTINUE

C HC11N C4H4S C4H4O C5H5N C6H4 C6H5O C6H6O C10H8 C10H16 C14H10a C14H10p C18H12t
C  325   326  327   328   329   330   331   332    333    334     335    336
C C18H12b C22H14c C24H12 C106H56
C  337     338     339    340
!      WRITE(7,2131)
!      WRITE(7,300) TEFF,GLOG,IDRAB1,IDRAB2,IDRAB3,IDRAB4,IDRAB5,IDRAB6
!      WRITE(7,2181)
!      DO 814 I=1,JTAU
!      WRITE(7,217) I,(log10(partryck(I,J)),J=325,340),I
!      presmo(2,i) = log10(max(1.d-99,partryck(i,2)))
!      do 1912 n= 1,nosmol
!      if (molname(n).eq.'ATOM'
!     &       .or.molname(n).eq.'H2H2'
!     &       .or.molname(n).eq.'H2HE'
!     &       .or.molname(n).eq.'H2He'
!     &                          ) go to 1912
!      presmo(kmol(n),i) =  log10(max(1.d-99,partryck(i,kmol(n))))
!1912  continue
!      HNIC(i) = log10(max(1.d-99,xmettryck(I,1)))
!814   CONTINUE
      END IF




!------------ WRITING PP FOR JUMP = 4 -----------


      IF (JUMP.EQ.4) then
        do j=1,jtau
                partpp(j,2)=partpp(j,10)
                partpp(j,10)=partpp(j,218)
        end do
! H, H-, H2, H2+, H2O, OH, CH, CO, CN, C2, N2, O2, NO, NH, TiO
! 0   1   2   3     4   5   6   7   8   9  10  11  12  13   31
        WRITE(7,2132)
        WRITE(7,300) TEFF,GLOG,IDRAB1,IDRAB2,IDRAB3,IDRAB4,IDRAB5,IDRAB6
        write(7,2142) (name_listmo(i),i=1,13),name_listmo(31)
2142    FORMAT(' K',5X,'P(H)',14(3x,a5))
        do i=1,jtau
                 write(7,215) i,(log10(partpp(i,0))),
     &         (log10(partpp(i,jj)),jj=1,13),(log10(partpp(i,31))),i
        end do

! C2H2, HCN, C2H, HS, SiH, C3H, C3, CS, SiC, SiC2, NS, SiN, SiO, SO, S2,
! SiS
!  14    15  16  18 19  20  21 22 23  24   25  26  27 28 29  30
        write(7,2132)
        write(7,300) teff,glog,idrab1,idrab2,idrab3,idrab4,idrab5,idrab6
        write(7,2143)(name_listmo(i),i=14,16),(name_listmo(i),i=18,30)
2143    format(' K ',16(3x,a5))
        do 813 i=1,jtau
                 write(7,217) i,(log10(partpp(i,j)),j=14,16),
     *           (log10(partpp(i,j)),j=18,30),i
813     continue

! VO, ZrO, MgH, CaH, LaO, CH4, CH2, CH3, NH3, CO2, Si2C, SiO2, H2S,
! CaOH, CHNO, SiF2
C 32  33  34  36  37 39  40  41  43  46  57   58   59   67  107 135
        WRITE(7,2132)
        WRITE(7,300) TEFF,GLOG,IDRAB1,IDRAB2,IDRAB3,IDRAB4,IDRAB5,IDRAB6
        WRITE(7,2144)(name_listmo(i),i=32,34)
     &          ,(name_listmo(i),i=36,37)
     &          ,(name_listmo(i),i=39,41)
     &          ,name_listmo(43)
     &          ,name_listmo(46)
     &          ,(name_listmo(i),i=57,59)
     &          ,name_listmo(67)
     &          ,name_listmo(107)
     &          ,name_listmo(135)
2144    FORMAT(' K ',16(3x,a5))
        DO 814 I=1,JTAU
                 WRITE(7,217) I,(log10(partpp(I,J)),J=32,34),
     *          (log10(partpp(I,J)),J=36,37),
     *          (log10(partpp(I,J)),J=39,41),
     *          (log10(partpp(I,J)),J=43,43),
     *          (log10(partpp(I,J)),J=46,46),
     *          (log10(partpp(I,J)),J=57,59),
     *          (log10(partpp(I,J)),J=67,67),
     *          (log10(partpp(I,J)),J=107,107),
     *          (log10(partpp(I,J)),J=135,135),I
814     CONTINUE

! TiH, CaH, FeH, CrH, LiH
! 214, 215, 216, 217, 71
        name_listmo(214) = 'TiH  '
        name_listmo(215) = 'CaH  '
        name_listmo(216) = 'FeH  '
        name_listmo(217) = 'CrH  '
        WRITE(7,2132)
        WRITE(7,300) TEFF,GLOG,IDRAB1,IDRAB2,IDRAB3,IDRAB4,IDRAB5,IDRAB6
        WRITE(7,2148) name_listmo(214),name_listmo(215),
     &          name_listmo(216),name_listmo(217),name_listmo(71)
        DO I=1,JTAU
                WRITE(7,'(I3,5F8.3,2X,I3)')I,(log10(partpp(I,J)),
     *           J=214,217),log10(partpp(I,71)),I
        end do
2148    FORMAT(' K ',5(3x,a5))
      END IF
      GO TO 2109

!------------ END OF  WRITING PP FOR JUMP = 4 -----------



2102  CONTINUE     ! go here if jump>2 (i.e. not old marcs chem.eq. and not JANAF fits).
      IF (JUMP.GT.4) GO TO 2103   !here one can put other chem.eq. possibilities
!       ERC: was GT3 before, but changed to GT4 to have Pgas printed in
!       SUB routine init_ggchem
C output from GEM:

       do 120 i=1,jtau
       pg(i)=pp(i)-ppr(i)-ppt(i)
CV20   phe(i) = pg(i) * totabk(i,3)/(totabk(i,2)+totabk(i,3))
       phe(i) = (pg(i)-pe(i)) * totabk(i,3)/(totabk(i,2)+totabk(i,3))
       trpe(i) = 0.1d0 * pe(i)
CV20   ptot(i) = pg(i) + pe(i)
       ptot(i) = pg(i)
       ptot(i) = 0.1d0 * ptot(i)     !1N/m^2 = 10 dynes/cm^2
       trphe(i) = (0.1d0 * phe(i)) / ptot(i)
120   CONTINUE
782    format(i3,f8.1,1p4e12.3)

      call tstgem(t,ptot,jtau)
C

       DO 3151 kd=1,jTAU
CV20   pe_gem(kd) = pg(kd) * abink(kd,1) / (1.d0 - abink(kd,1))
       pe_gem(kd) = pg(kd) * abink(kd,1)
CV20   p_particles = pg(kd) + pe_gem(kd)         !pg(input)+pe(computed) in dyn/cm^2
       p_particles = pg(kd)                      !pg(input)+pe(computed) in dyn/cm^2
       pp_sum(kd) = 0.
       DO 3153 km=1,nspec
       PRESPP(kd,km) = p_particles * abink(kd,km)  !(pg+pe)*rel.pp = pp in dyn/cm^2
       if(km.gt.1) pp_sum(kd) = prespp(kd,km) + pp_sum(kd)
3153   continue
3151   continue

        write(6,*) ' pg/pp_sum {1-jtau} in listmo after GEM call:'
        write(6,2126) ((pg(kd)/pp_sum(kd)),kd=1,jtau)
        write(6,*) ' pe/pe_gem {1-jtau}:'
        write(6,2126) ((pe(kd)/pe_gem(kd)),kd=1,jtau)
2126    format(1p8e10.3)

       do 2209 i=1,jtau
CV20   ptot(i)=pe(i)+pg(i)
       ptot(i)=pg(i)
2209   continue

      k = 0
      do 2201 j=1,nspec
      if (listpp(j).eq.1) then
       k = k + 1
       printname(k) = name_gem(j)
       do 2202 i=1,jtau
C2202    printpp(i,k)=abink(i,j)*ptot(i)
2202    printpp(i,k)=prespp(i,j)
      end if
2201  continue

      kpp = k

      write(7,2210) (printname(j),j=1,12) !kpp)
      do 2203 i=1,jtau
       write(7,2211) (printpp(i,j),j=1,12) !kpp)
 2203 continue

      write(7,2210) (printname(j),j=13,24)
      do 2204 i=1,jtau
       write(7,2211) (printpp(i,j),j=13,24)
 2204 continue

      write(7,2210) (printname(j),j=25,36)
      do 2205 i=1,jtau
       write(7,2211) (printpp(i,j),j=25,36)
 2205 continue

      write(7,2210) (printname(j),j=37,48)
      do 2206 i=1,jtau
       write(7,2211) (printpp(i,j),j=37,48)
 2206 continue

      write(7,2210) (printname(j),j=49,kpp)
      do 2207 i=1,jtau
       write(7,2211) (printpp(i,j),j=49,kpp)
 2207 continue

2210  format(/4x,12(a8,2x))
!2211  format(12f8.3)
2211  format(1p12e10.3)

      GO TO 2109
2103  CONTINUE
2109  CONTINUE     ! go here when finished the partial pressure writing



1905  FORMAT(I3,18(2X,A4))
      WRITE(7,1904) NOSMOL
1904  FORMAT
     *(I3,' MOLECULES(/+ATOMS) HAVE BEEN CONSIDERED IN THE OPACITY')

!      WRITE(7,*) ' '
!      WRITE(7,1910)
!1910  format(/' P A R T I A L   P R E S S U R E S ',
!     * '  O F   O P A C I T Y  SP ETC')
C      WRITE(7,300) TEFF,GLOG,IDRAB1,IDRAB2,IDRAB3,IDRAB4,IDRAB5,IDRAB6

C      write(7,1919)
C1919  format(' K log10:  P6  P6_jon  P6_TS  PG  P_sum_mol  SUM_P_OS'
C     &  ,' P_Ts_add  P_neu15  P_ion15  2*P_HI  2*P_He  P_HeI  P_He+')
C      WRITE(7,*) ' '
C     jumptemp = 1
C     if(jumptemp.eq.1) then
C       write(7,*) ' temporarily we dont write these quantities'
C       go to 1922
C     end if

!      DO 1920 k=1,JTAU
!      posmol = 0.
!      do 1911 n= 1,nosmol
!        print *, n, kmol(n)
!      if (molname(n).eq.'ATOM'
!     &       .or.molname(n).eq.'H2H2'
!     &       .or.molname(n).eq.'H2HE'
!     &       .or.molname(n).eq.'H2He'
!     &                          ) go to 1911
!      posmol = 10.**max(-20.D0,presmo(kmol(n),k)) + posmol
!1911  continue
!         addtsmol = 0.
!      DO 1918 n=32,342
!1918     addtsmol = addtsmol + partryck(k,n)
!C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
!        P6 =  10.**max(-20.0D0,HNIC(K))
!     &      + 0.42 * 10.**max(-20.0D0,PRESMO(17,K))
!     &      + 0.85 * 10.**max(-20.0D0,PRESMO(2,K))
!        P6TS =  xmettryck(k,1) + 0.42*xmettryck(k,2) 
!     &      + 0.85*10.**max(-20.0D0,PRESMO(2,K))
!        Pneu_15 = P_NEU_HCNO(K) + P_NEU_HE(K)
!        Pion_15 = P_ION_HCNO(K) + P_ION_HE(K)
!
!C        WRITE(7,1921) K
!C     &   ,log10(max(1.d-99,P6))
!C     &   ,log10(max(1.d-99,P6_JON(K)))
!C     &   ,log10(max(1.d-99,P6TS))
!C     &   ,log10(max(1.d-99,PG(K)))
!C     &   ,log10(max(1.d-99,SUMPMOL(K)))
!C     &   ,log10(max(1.d-99,POSMOL))
!C     &   ,log10( max(1.d-99,addtsmol) )
!C     &   ,log10(max(1.d-99,Pneu_15))
!C     &   ,log10(max(1.d-99,Pion_15))
!C     &   ,log10(max(1.d-99,xmettryck(K,1))),HNIC(K)
!C     &   ,log10(max(1.d-99,xmettryck(K,2))),PRESMO(17,K)
!C     &   ,log10(max(1.d-99,P_NEU_HE(K)))
!C     &   ,log10(max(1.d-99,P_ION_HE(K)))
!1920  CONTINUE
!1921  FORMAT(i3,15f8.3)
!1922  CONTINUE

      WRITE(7,1900)
1900  format(' A B S O R P T I O N   C O E F F I C I E N T S  [CM/MOL]')
      WRITE(7,300) TEFF,GLOG,IDRAB1,IDRAB2,IDRAB3,IDRAB4,IDRAB5,IDRAB6

      WRITE(7,1906) (MOLNAME(I),I=1,NOSMOL)
1906  FORMAT(' K',25(2X,A4))
      DO 2074 k=1,JTAU
      do 2076 i=1,nosmol
2076  SUMOP(I,K)= max(1.d-99, SUMOP(I,K) )
      write(7,2073) k,( log10(SUMOP(I,K)),I=1,NOSMOL )
2074  continue

      WRITE(7,1901)
1901  format(' I N T E G R A T E D   O P A C I T Y   [CM/G*]')
      WRITE(7,300) TEFF,GLOG,IDRAB1,IDRAB2,IDRAB3,IDRAB4,IDRAB5,IDRAB6

      WRITE(7,1906) (MOLNAME(I),I=1,NOSMOL)
      DO 2072 k=1,JTAU
      do 2077 i=1,nosmol
2077  SUMKAP(I,K)= max(1.d-99, SUMKAP(I,K) )
      write(7,2073) k,( log10(SUMKAP(I,K)),I=1,NOSMOL )
2072  continue
2073  format(i3,25f6.2)
C End of molecular partial pressures

      WRITE(7,260)

      READ(IARCH)(NJ(I),I=1,NEL),NLP,(XLR(I),I=1,NLP)
     & ,NPROV,NPROVA,NPROVS,(ABNAME(KP),SOURCE(KP),KP=1,NPROV)
C
C
CUGJ STORE ABSORBTIONCOEFFICIENT, SCATTERING AND LAMBDA FOR SPECTRUM 
CUGJ CALCULATIONS
C
      WRITE(7,*) ' '
      WRITE(7,*) ' '
      WRITE(7,*) 'LAMBDA FOR ABS-. AND SCAT.COEF. = '
      WRITE(7,*) NLP
      WRITE(7,530) (XLR(KLAM),KLAM=1,NLP)
530   FORMAT(10F12.3)
531   FORMAT(I4,E12.3)
532   FORMAT(1P10E12.5)
      WRITE(7,*) ' '
      WRITE(7,*) 'ABSORBTION COEF. AND SCATTERING COEF. I EACH TAU:'
      WRITE(7,*) ' '
      DO 81 KT=1,JTAU
      READ(IARCH) (ABSKTR(KLAM),KLAM=1,NLP)
      READ(IARCH) (SPRTR(KLAM),KLAM=1,NLP)
      WRITE(7,531) KT,TAU(KT)
      WRITE(7,532) (ABSKTR(KLAM),KLAM=1,NLP)
      WRITE(7,532) (SPRTR(KLAM),KLAM=1,NLP)
81    CONTINUE
      WRITE(7,*) ' '
C
C first, identify the depth points where logtau_ross is -4, -3, ..., 2
C for printing of opacities
       kdp(1) = 1
       kn = 2
      DO 233 K=1,JTAU
       taul = log10(tau(k))
      if (taul.gt.xdp(kn)) then
       k1 = max(k-1,1)
       kdp(kn) = k
       if ( abs(taul-xdp(kn)) .gt. abs( xdp(kn)-log10(tau(k1)) )  )
     *     kdp(kn) = k1
       if (kn.eq.8) go to 234
       kn = kn + 1
      end if
233   continue
234   continue
       kdp(8) = jtau

      HCK=143922240.
C     DO 22 KTAU=1,NTPO
      DO 22 KTAU=1,JTAU
      DO 20 IE=1,NEL
      NJP=NJ(IE)
      READ(IARCH) KR,TAUI,TI,PEI,IEL(IE),abmarcs(IE,1),
     &            (ANJON(IE,JJ),JJ=1,NJP),(PART(IE,JJ),JJ=1,NJP)
   20 CONTINUE
      DO 21 KLAM=1,NLP
      READ(IARCH) KR,TAUIL,(PROV(J,KLAM),J=1,NPROV),
     &            ABSKA(KLAM),SPRIDA(KLAM)
   21 CONTINUE

C write only for selected layers:
      do 236 kd=1,8
      if (ktau.eq.kdp(kd)) go to 238
236   continue
      go to 22
238   continue

      IF(KTAU.LE.1) WRITE(7,218)
      IF(KTAU.GT.1) WRITE(7,219)
      WRITE(7,300) TEFF,GLOG,IDRAB1,IDRAB2,IDRAB3,IDRAB4,IDRAB5,IDRAB6
      WRITE(7,220) TAUI
      WRITE(7,221) TI,PEI,kr,log10(tau(kr))
      WRITE(7,222)
      PESUM=0.
      DO 32 IE=1,NEL
      PEP(IE)=0.
      NJP=NJ(IE)
      IF(NJP.LT.2) GOTO 32
      DO 33 JJ=2,NJP
   33 PEP(IE)=PEP(IE)+abmarcs(IE,1)*ANJON(IE,JJ)*(JJ-1)
   32 PESUM=PESUM+PEP(IE)
      DO 19 IE=1,NEL
      NJP=NJ(IE)
      PEP(IE)=PEP(IE)/PESUM
      WRITE(7,223)IEL(IE),abmarcs(IE,1),PEP(IE),(ANJON(IE,JJ),JJ=1,NJP)
C     WRITE(7,224) (PART(IE,JJ),JJ=1,NJP)
C
C
   19 CONTINUE
      WRITE(7,2250)
      WRITE(7,225) (ABNAME(KP),KP=1,NPROV)
C      iwrop = 76 + ktau
C      WRITE(iwrop,2252) KR,log10(TAUI)
C      WRITE(iwrop,2253) 
C      WRITE(iwrop,2251) (ABNAME(KP),KP=1,NPROV-2),ABNAME(NPROV)
      HCKT=HCK/TI
      DO 18 KLAM=1,NLP
      ABKLA=ABSKA(KLAM)
      STIM=1.-EXP(-HCKT/XLR(KLAM))
C      if (lin_cia.ne.1) 
C     *   abkla = abkla+(prov(nprova,klam)+prov(nprova-1,klam))*stim
C      WRITE(iwrop,2261)XLR(KLAM)/1.d4
C     *    ,log10( max(1.0d-99,ABSKA(KLAM)) )
C     *    ,log10( max(1.0d-99,SPRIDA(KLAM)) )
C     *    ,(log10( max(1.0d-99,PROV(J,KLAM)) ),J=1,NPROVA-2)
C     *    ,(log10( max(1.0d-99,PROV(J,KLAM)) ),J=NPROVA+1,NPROVA+NPROVS)
      DO 180 J=1,NPROVA-2
  180 PROV(J,KLAM)=PROV(J,KLAM)/ABKLA*STIM
      if (lin_cia.ne.1) 
     *   abkla = abkla+(prov(nprova,klam)+prov(nprova-1,klam))*stim
      DO 1801 J=NPROVA-1,NPROVA
 1801 PROV(J,KLAM)=PROV(J,KLAM)/ABKLA*STIM
      DO 181 J=1,NPROVS
  181 PROV(NPROVA+J,KLAM)=PROV(NPROVA+J,KLAM)/SPRIDA(KLAM)
      WRITE(7,226) XLR(KLAM),ABSKA(KLAM),SPRIDA(KLAM),
     &             (PROV(J,KLAM),J=1,NPROV)
   18 CONTINUE
   22 CONTINUE

C      write(7,288)
C      WRITE(7,300) TEFF,GLOG,IDRAB1,IDRAB2,IDRAB3,IDRAB4,IDRAB5,IDRAB6
C      write(7,287) (molname(nm),nm=1,nosmol)
CC      write(7,*) ' followed by partp for the same molecules:'
CC      write(7,*) ' .. akapmol for the same molecules:'
C      do 280 itx = 1,itxmax
CC      it = max(1,(itx-1)*10)
C        if(itx.eq.1) then
C          it=7
C        else if (itx.eq.2) then
C          it = 27
C        else if (itx.eq.3) then
C          it = 40
C        end if
C      iwrop = 76 + itx
C      write(iwrop,2871) nosmol,(molname(nm),nm=1,nosmol)
C      do 280 jvx = jvxmax,1,-1
C      jv = jvx * 1000
C      write(iwrop,2863) 1.d4/wnos(jv)
C     *     ,log10( max(1.0d-99,xconop(jvx,itx)) )
C     *     ,log10( max(1.0d-99,xlineop(jvx,itx)) )
C     *     ,(log10( max(1.0d-99,oppr(nm,itx,jvx,1)) ),nm=1,nosmol)
C     *     ,(log10( max(1.0d-99,oppr(nm,itx,jvx,1)) ),nm=1,5)
C280   continue

C        WRITE(7,55)
C        AWTOT=DFLOAT(NWTOT)
C        DO 5518 K=1,JTAU
C        WRITE(6,54) K,SSUM(K)/AWTOT,XSUM(K)/AWTOT,CONSUM(K)/AWTOT
C5518      CONTINUE
C54    FORMAT (I4,1P3E14.4)
C55    FORMAT (///' average continiums and line absorption coefficient,'
C     * ,/' scatter<s(k)>, cont.abs.<x(k)> and line.abs.<conos(k)>'
C     * ,' over the whole spectrum :')


C - skip output of colours:
      NOCOL=1
      IF(NOCOL.EQ.1) GO TO 1
C
C
      READ(IARCH) NLB,(XLB(J),FLUXME(J),J=1,NLB),(W(J),J=1,NLB)
C        CONVERT TO 'PHYSICAL' FLUXES
      DO 24 J=1,NLB
   24 FLUXME(J)=3.14159*FLUXME(J)
      DO 25 J=1,NLB
      JNORM=J
      IF(XLB(J).GT.5000.) GOTO 26
   25 CONTINUE
   26 STMAGN=-2.5*log10(FLUXME(JNORM))
C
C        WRITE(7,FLUXES)
C
      WRITE(7,521)
      N=0
      DO 52 I=1,NLB,4
      WAVFLX(N+1)=(XLB(I)+XLB(I+1)+XLB(I+2)+XLB(I+3))/4.
      WAVFLX(N+2)=100.*(A*(FLUXME(I)+FLUXME(I+3))+B*(FLUXME(I+1)
     1+FLUXME(I+2)))
      N=N+2
      IF(N.LT.10) GO TO 52
      WRITE(7,520) WAVFLX
      N=0
   52 CONTINUE
      IF(N.NE.0) WRITE(7,520)(WAVFLX(I),I=1,N)
  520 FORMAT(1X,5(F12.0,E13.5))
  521 FORMAT('1  BELL''S FLUX              (CENTER WAVE & FLUX)'//)
C
C
      WRITE(7,219)
      WRITE(7,300) TEFF,GLOG,IDRAB1,IDRAB2,IDRAB3,IDRAB4,IDRAB5,IDRAB6
      WRITE(7,227)
      LLB=NLB
   60 IF(LLB.LT.NLB) WRITE(7,219)
      IF(LLB.LT.NLB) WRITE(7,300) TEFF,GLOG,IDRAB1,IDRAB2,IDRAB3,
     &                            IDRAB4,IDRAB5,IDRAB6
      IF(LLB.LT.NLB) WRITE(7,227)
      WRITE(7,228)
      IRAD = 0
      IF(LLB.LT.NLB) GOTO 65
      I=0
      I2=50
      I3=100
      I4=150
      GOTO 650
   65 CONTINUE
      I = I + 149
      I2 = I2 + 149
      I3 = I3 + 149
      I4 = I4 + 149
  650 CONTINUE
      ILI=2
      IF(LLB.GE.51) ILI=3
      IF(LLB.GE.101)ILI=4
      IF(LLB.GE.151)ILI=5
   66 I4=I4+1
   67 I3=I3+1
   68 I2=I2+1
   69 I=I+1
      IRAD = IRAD + 1
      IF(IRAD.GT.50) GOTO 70
      GOTO(70,61,62,63,64),ILI
   64 WRITE(7,229) XLB(I),FLUXME(I),FLUMAG(I),XLB(I2),FLUXME(I2),
     &             FLUMAG(I2),XLB(I3),FLUXME(I3),FLUMAG(I3),
     &             XLB(I4),FLUXME(I4),FLUMAG(I4)
      IF(I4.GE.NLB) ILI=4
      GOTO 66
   63 WRITE(7,229) XLB(I),FLUXME(I),FLUMAG(I),XLB(I2),FLUXME(I2),
     &             FLUMAG(I2),XLB(I3),FLUXME(I3),FLUMAG(I3)
      IF(I3.GE.NLB) ILI=3
      GOTO 67
   62 WRITE(7,229) XLB(I),FLUXME(I),FLUMAG(I),XLB(I2),FLUXME(I2),
     &             FLUMAG(I2)
      IF(I2.GE.NLB) ILI=2
      GOTO 68
   61 WRITE(7,229) XLB(I),FLUXME(I),FLUMAG(I)
      IF(I.GE.NLB) GOTO 70
      GOTO 69
   70 LLB = LLB - 200
      IF(LLB.GT.0) GOTO 60
C
C        COMPUTE U, B, V, R, I AND COLOURS
      SUMCOL=0.
      SUMNOR=0.
      WRITE(7,270)
      DO 40 I=1,NLB
      XXLB=XLB(I)
      IF(XXLB.LT.3000..OR.XXLB.GT.4200.) GOTO 40
      INDEX=INT(XXLB/100.)-29
      VIKT=W(I)*UW(INDEX)
      SUMNOR=SUMNOR+VIKT
      SUMCOL=SUMCOL+VIKT*FLUXME(I)
   40 CONTINUE
      UFLUX=SUMCOL/SUMNOR
      SUMCOL=0.
      SUMNOR=0.
      DO 41 I=1,NLB
      XXLB=XLB(I)
      IF(XXLB.LT.3500..OR.XXLB.GT.5600.) GOTO 41
      INDEX=INT(XXLB/100.)-34
      VIKT=W(I)*BW(INDEX)
      SUMNOR=SUMNOR+VIKT
      SUMCOL=SUMCOL+VIKT*FLUXME(I)
   41 CONTINUE
      BFLUX=SUMCOL/SUMNOR
      SUMCOL=0.
      SUMNOR=0.
      DO 42 I=1,NLB
      XXLB=XLB(I)
      IF(XXLB.LT.4700..OR.XXLB.GT.7200.) GOTO 42
      INDEX=INT(XXLB/100.)-46
      VIKT=W(I)*VW(INDEX)
      SUMNOR=SUMNOR+VIKT
      SUMCOL=SUMCOL+VIKT*FLUXME(I)
   42 CONTINUE
      VFLUX=SUMCOL/SUMNOR
      DO 43 I=1,NLB
      JNORM=I
      IF(XLB(I).GE.7000.) GOTO 44
   43 CONTINUE
   44 RFLUX=FLUXME(JNORM)
      XINTCALL = 9000.0D+0
      CALL TINT(NLB,XLB,FLUXME,XINTCALL,XIFLUX)
      UMAG=-2.5*log10(UFLUX)-STMAGN
      BMAG=-2.5*log10(BFLUX)-STMAGN
      VMAG=-2.5*log10(VFLUX)-STMAGN
      RMAG=-2.5*log10(RFLUX)-STMAGN
      XIMAG=-2.5*log10(XIFLUX)-STMAGN
      UB=UMAG-BMAG
      BV=BMAG-VMAG
      UV=UMAG-VMAG
      RI=RMAG-XIMAG
      VR=VMAG-RMAG
      VI=VMAG-XIMAG
      WRITE(7,300) TEFF,GLOG,IDRAB1,IDRAB2,IDRAB3,IDRAB4,IDRAB5,IDRAB6
      WRITE(7,2698)
      WRITE(7,271) UB,BV,UV
      WRITE(7,272) XLB(JNORM),RI,VR,VI
      WRITE(7,273) UMAG,BMAG,VMAG,RMAG,XIMAG
    1 CONTINUE
  200 FORMAT(//' M O D E L  P A R A M E T E R S  E T C .'//)
  282 FORMAT(' * THE FOLLOWING MODEL WAS COMPUTED BY  S C M A R C S  ',
     +F5.1,2X,A9,2X,A9,' *')
  283 FORMAT(1X,81('*'))
  284 FORMAT(1X,81('*')////)
  201 FORMAT(' EFFECTIVE TEMPERATURE  ',F12.0,'  KELVIN'/' TOTAL ',
     +'FLUX (= SIGMA*TEFF**4)',1PE11.3,'  ERGS/S/CM**2'/' ACCELERATION',
     +' OF GRAVITY',0PF12.1,'  CM/S**2 (i.e., log(g) = ',f5.2,')'
     +/'Z/Zo (i.e. Oxygen/Oxygen_o) =',1PE8.1,'  C/O =',0PF6.2
     +//' CONVECTION PARAMETERS'/
     +' PALFA (L/HP)=',F5.2,',  PNY (NY)=',F5.2,',  PY (Y)=',F6.3)
 2017 FORMAT(' EFFECTIVE TEMPERATURE  ',F12.0,'  KELVIN'/' TOTAL ',
     +'FLUX (= SIGMA*TEFF**4)',1PE11.3,'  ERGS/S/CM**2'/' ACCELERATION',
     +' OF GRAVITY',0PF12.5,'  CM/S**2 (i.e., log(g) = ',f5.2,')'
     +/'Z/Zo (i.e. Oxygen/Oxygen_o) =',1PE8.1,'  C/O =',E8.1
     +//' CONVECTION PARAMETERS'/
     +' PALFA (L/HP)=',F5.2,',  PNY (NY)=',F5.2,',  PY (Y)=',F6.3)
 2011 FORMAT(' Convection was, however, excluded in the uppermost'
     +    ,I4,' layers')
 2012 FORMAT(' XMAX, TAUM, facply, moltsuji = ',1pe8.1,0p2f8.1,i3)
 2013 FORMAT(' Molecular equilibrium was treated by Tsujis routine')
 2014 FORMAT(' Molecular equilibrium was treated by the Marcs routin')
 2015 FORMAT(A60)
 2018 FORMAT(' Molecular equilibrium was treated by JF fit to JANAF 99')
 2019 FORMAT(' Molecular equilibrium was by Gibbs Energy Minimisation')
  256 FORMAT(' TURBULENCE PRESSURE IS NEGLECTED (PBETA<=0.)'/)
  257 FORMAT(' TURBULENCE PRESSURE IS INCLUDED AND PBETA=',F5.2/)
  250 FORMAT(' CONVECTION HAS BEEN INCLUDED IN THIS MODEL; ISTRAL=',I3)
  251 FORMAT(' CONVECTION HAS  N O T  BEEN INCLUDED IN THIS MODEL; '
     *  ,'ISTRAL=',I3)
  252 FORMAT(' LINE BLANKETING HAS  N O T  BEEN INCLUDED IN THIS MODEL;'
     * ,' ILINE=',I3)
 2529 FORMAT(' CIA described by Linskys formulas was included')
  253 FORMAT(' LINE BLANKETING HAS BEEN INCLUDED IN THIS MODEL;'
     * ,' ILINE=',I3)
 2531 FORMAT(' THE MASS WAS ASSUMED TO BE RELM =',F5.1)
 2539 FORMAT(' -- but this computation is plane parallel ')
 2541 FORMAT(' For',f5.0'<T<',f5.0,'K',f6.1,
     * '-',f6.1,' % of C2H2 was assumed to be grains, and'/,
     *  ' for',f5.0'<T<',f5.0,'K',f6.1,
     * '-',f6.1,' % was assumed to be molcular C2H2')
 2534 FORMAT(' THE OS SAMPLING HAVE BEEN PERFORMED IN',I3,' INTERVALS')
 2535 FORMAT(' FROM',F6.1,' CM-1 TO',F9.1,' CM-1')
 2536 FORMAT(' (..which is identical to from',F9.1,' Å to',F6.1,' mu)')
 2537 FORMAT(' THE CHOSEN INTERVAL-BEGINNING AND STEPLENGTH IN CM-1:')
 2538 FORMAT(8F9.1)
  254 FORMAT(' THE STROEMGREN EQUATION HAS BEEN USED FOR THE UPPERMOST',
     *I4,' POINTS'/' IDENTIFICATION    ',6A4)
  255 FORMAT(////' THE FOLLOWING MODEL WAS OBTAINED AFTER',I3,
     *' ITERATION(S)')
  202 FORMAT(/////'  LOG. ABUNDANCES USED IN MODEL CALCULATIONS'/2X,'H'
     *,5X,'HE',4X,'C',5X,'N',5X,'O',5X,'NE',4X,'NA',4X,'MG',4X,'AL',4X,
     *'SI',4X,'S',5X,'K',5X,'CA',4X,'CR',4X,'FE',4X,'NI',4X,'TI')
  203 FORMAT (20F6.2)
  204 FORMAT(//' C O R R E C T I O N S  I N  T H E  L A S T ',
     *' I T E R A T I O N')
  205 FORMAT(2X,' K',4X,'TAUROSS',6X,'T',7X,'DELTA(T)',5X,'FCONV',
     *4X,'DELTA(FCONV)',2X,'K')
  206 FORMAT(I4,1PE12.3,0P2F11.2,1PE12.3,E12.2,I4)
  207 FORMAT(//' M O D E L  A T M O S P H E R E   (CGS UNITS)')
  208 FORMAT(3X,' K',4X,'TAUROSS',3X,'TAU(5000)',2X,'GEOM. DEPTH',5X,
     *'T',8X,'PE',10X,'PG',9X,'PRAD',8X,'PTURB',5X,'KAPPAROSS',5X,'K')
  209 FORMAT(I5,1P2E12.4,E13.4,0PF9.2,1P5E12.4,I4)
 2061 FORMAT(//' S U M   O F   P A R T I A L    P R E S S U R E S')
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789
 2063 FORMAT('  K  P_MOL, P_NEU_HCNO, P_ION_HCNO, P_NEU_HE'
     &,' P_ION_HE, P_NON_HHECNO, P6_ION PG_JON, P6 PHe')
 2071 FORMAT(//' M A S S  L O S S AND S P H E R E AND  L E V I T A T',
     *' I O N')
 2081 FORMAT(2X,' K',2X,'TauRoss',2X,'log(tau)',2X,'density,g/cm3',2X,
     *'Z,cm',2x,'T',7x,'arad',5X,'dpr/dz',4X,'V-loss,km/s',3x,'rossos'
     * ,4x,'p-rossos',6x,'cross')
 2091 FORMAT(I5,1PE11.2,0PF6.2,1P2E11.2,0PF7.0,1P2E11.2,0PF8.4,
     * 1P3E12.3)
  210 FORMAT(//' T H E R M O D Y N A M I C A L  Q U A N T I T I E S ',
     *' AND  C O N V E C T I O N  (CGS UNITS)')
  211 FORMAT(2X,' K',3X,'TAUROSS',3X,'DENSITY',5X,'MU',9X,'CP',10X,
     *'CV',7X,'ADGRAD',9X,'Q',7X,'SOUND VEL.',2X,'CONV. VEL.',1X,
     *'FCONV/F',3X,'K')
  212 FORMAT(I4,1P2E11.3,0PF8.3,6(1PE12.3),0PF9.5,I3)
  213 FORMAT(//' L O G A R I T H M I C  M O L E C U L A R ',
     *' P A R T I A L  P R E S S U R E S   (CGS UNITS)')
 2131 FORMAT(//' L O G A R I T H M I C  M O L E C U L A R ',
     *' P A R T I A L  P R E S S U R E S   (Tsuji)')
 2132 FORMAT(//' L O G A R I T H M I C  M O L E C U L A R ',
     *' P A R T I A L  P R E S S U R E S   (GGCHEM)')
  214 FORMAT(' K',3X,'P(H)',2X,'P(H-)',2X,'P(H2)',2X,'P(H2+)',1X,
     *'P(H2O)',1X,'P(OH)',2X,'P(CH)',2X,'P(CO)',2X,'P(CN)',2X,'P(C2)',2X
     *,'P(N2)',2X,'P(O2)',2X,'P(NO)',2X,'P(NH)',1X,'P(TIO)',4X,'K')
  215 FORMAT(I3,15f8.3,2X,I3)
 2141 FORMAT(' K',2X,'TiO',5X,'TiO2',4X,'TiC2',4X,'TiO-KE',5X,
     *    'P(He)',3x,'P(HI)',3x,'P(Hx2)',3X,'P6',6X,'PG',5X,' K')
 2151 FORMAT(I3,9F8.3,2X,I3)
  216 FORMAT(' K',4X,'C2H2',4X,'HCN',5X,'C2H',5X,'HS',6X,'SIH',5X
     * ,'C3H',5X,'C3',6X,'CS',6X,'SIC',5X,'SIC2',4X,'NS',6X,'SIN',5X
     * ,'SIO',5X,'SO',6X,'S2',6X,'SIS',6X,'K')

 2161 FORMAT(' K',3x,' VO ',3x,'ZrO',4x,'CaH',4x,'LaO',5x,'CH4',4x
     * ,'CH2',4x,'CH3',4x,'NH3',4x,'CO2',4x,'Si2C',3x,'SiO2',3x
     * ,'H2S',4x,'CS2',4x,'CaOH',3x,'NO2',4x,'AlH',5x,'K ')

  217 FORMAT(I3,16F8.3,2X,I3)
 2171 FORMAT(' K',3x'AlO ',3x,'SiH4',3x,'SO2',4x,'YO',5x,'C5H',4x
     * ,'C4H',4x,'C4',5x,'C5 ',4x,'C6H ',3x,'C4H6',3x,'C6H2',3x,'C6H4'
     * ,3x,'C10H7',2x,'HC5N',3x,'HC7N',3x,'HC9N',3x,'K')

 2181 FORMAT(' K',3x,'HC11N',2x,'C4H4S',2x,'C4H4O',2x,'C5H5N',2x
     * ,'C6H4',3x,'C6H5O',2x,'C6H6O',2x,'C10H8',2x,'C10H16',1x,
     * 'C14H10a','C14H10p','C18H12t','C18H12b','C22H14c','C24H12'
     * ,1x,'C106H56','K')

  218 FORMAT(//' I O N I Z A T I O N  C O N D I T I O N S  AND ',
     *' A B S O R P T I O N  M E C H A N I S M S')
  219 FORMAT(1H1)
  220 FORMAT(/' ****************'/' * TAU=',F8.4,' *'/
     *' ****************')
  221 FORMAT(/' T=',F7.0,'  PE=',1PE9.2,' lag:',I3,' lg10_tau(kr)='
     #    ,0pf7.2)
  222 FORMAT(' ELEMENT  ABUNDANCE'
     *                      ,' ELCONT IONIZATION FRACTIONS',17X,
     *'PARTITION FUNCTIONS'/30X,'I',7X,'II',6X,'III',5X,'IV',12X,
     *'I',9X,'II',8X,'III',7X,'IV')
  223 FORMAT(6X,A2,1PE12.3,0PF6.3,4(F8.4))
  224 FORMAT(1H+,62X,4(1PE10.2))
 2250 FORMAT(////16X,'F R A C T I O N S  O F  C O N T I N U O U S  A B',
     *' S O R P T I O N  A N D  S C A T T E R I N G')
  225 FORMAT(' WAVELENGTH   ABS       SCAT    ',16(1X,A5))
2251  FORMAT('  W-mu log->ABS SCAT',3x,16A7)
2252  FORMAT(' k-depth-point in atmosphere:',i3,' ~log(tau_ross)=',f8.2)
2253  FORMAT(' Units: log_10 of values of opacity in cm^2/g*' )
  226 FORMAT(F11.0,1p2E10.3,0p16F6.3)
2261  FORMAT(F6.2,18f7.2)
2871  format(' wn_mu cont.op line-op. log_partp*akapmol for the',
     *        i4,' molecules:',15(2x,a4,1x))

288   format(///' O P A C I T Y   OF   O S - M O L E C U L E S',
     *  '  in log-cm^2/g-*')
287   format('  log(tau)  cont.op  line-op. partp*akapmol for:'
     *        ,15(2x,a4,1x))
286   format(i3,f7.2,15f7.2)
2861   format(i3,2f7.2,1p11e10.3,/10e10.3)
2862   format(1p13e10.3)
2863   format(f7.2,16f7.2)
285   format(' os-wn',i5,':',f8.1,' (=',f8.4,' mu)')
C285   format(' os-wn',i5,':',f8.1)
  227 FORMAT(//' F L U X E S (PHYSICAL FLUXES IN ERGS/S/CM**2/ANGSTROM)'
     */)
  228 FORMAT(4(6X,'LAMBDA',5X,'FLUX',7X,'MAGN'))
  229 FORMAT(4(F12.0,1PE11.2,0PF9.3))
  260 FORMAT(////' MOLECULES OTHER THAN H2 AND H2+ HAVE NOT BEEN',
     *' CONSIDERED IN THE TEMPERATURE - ELECTRON PRESSURE - PRESSURE',
     *' BALANCE')
 2698 FORMAT(1H0/' (UBV TRANSMISSION FUNCTIONS AFTER MATTHEWS +',
     *' SANDAGE.  AIR MASS = 1)')
  270 FORMAT(' G I A N T  L I N E  C O L O U R S')
  271 FORMAT(8X,'U - B =',F8.3//9X,'B - V =',F8.3//9X,'U - V =',F8.3
     *////)
  272 FORMAT(' TENTATIVE CONTINUUM COLOURS ',
     *' (R AT',F6.0,'A AND I AT 9000A)'//9X,
     *'R - I =',F8.3,//9X,'V - R =',F8.3//9X,'V - I =',F8.3)
  273 FORMAT(////' U =',F10.3,'  B =',F10.3,'  V =',F10.3,'  R =',F10.3,
     *'  I =',F10.3)
  300 FORMAT(/' TEFF=',F6.0,' LOG G=',F5.1,2X,6A4/)
      RETURN
         END
C
      SUBROUTINE MAINB
      implicit real*16 (a-h,o-z)
C
      include 'parameter.inc'
C
      COMMON /UTPUT/IREAD,IWRIT
      COMMON /CG/GRAV,KONSG /CTEFF/TEFF,FLUX
      COMMON /CSTYR/MIHAL,NOCONV /CXMAX/XMAX /CTAUM/TAUM
      COMMON /MIXC/PALFA,PBETA,PNY,PY /CVFIX/VFIX                          
      COMMON /CANGLE/XMY(6),XMY2(6),H(6),MMY
      COMMON /CARC1/ISTRAL,IDRAB1,IDRAB2,IDRAB3,IDRAB4,IDRAB5,IDRAB6,
     &              IARCH
      COMMON /CLINE4/ILINE
      COMMON /NATURE/BOLTZK,CLIGHT,ECHARG,HPLNCK,PI,PI4C,RYDBRG,
     *STEFAN
      COMMON /CSPHER/TAURAT,RADIUS,RR(NDP),NCORE
      COMMON /CPOLY/FACPLY,MOLTSUJI
      COMMON /CMETBL/METBL
      COMMON /MASSE/RELM
      COMMON /CLEVETAT/GEFF(NDP),PPRG(NDP),AMLOSS
      COMMON /CLIN/lin_cia
      COMMON /CORRECT/KORT,KPP,TCONV
      COMMON /TAUC/TAU(NDP),DTAULN(NDP),JTAU
      COMMON /STATEC/DUM1(10*NDP),TAULN(NDP),RO(NDP),NTAU,ITER
      COMMON /Cspec/spec(nwl,3),ispec
      common /cosexp/ lops,nops
      common /cindiam/tdiam1,tdiam2,fdiam1,fdiam2,
     *    tc2h21,tc2h22,fc2h21,fc2h22
      DATA TSUN,GSUN,RSUN/5800.,4.44,7E10/

C      print*,' in MAINB, vor schreiben nach mxms7.output'
C
C WRITE OUT THE PARAMETERS CHOSEN
      WRITE(6,*) ' THE PARAMETERS CHOSEN IN parameter.inc:'
      WRITE(6,*) ' NDP,NRAYS,mkomp,mkompr,ifadim,NWL,mtemp,mpe='
      WRITE(6,666) NDP,NRAYS,mkomp,mkompr,ifadim,NWL,mtemp,mpe
666   FORMAT(8I7)
      WRITE(6,*) ' NTAU,JTAU at beginning of MAINB =',NTAU,JTAU
C
C INITIATIONS
      HPLNCK=6.62554E-27
      BOLTZK=1.38054E-16
      CLIGHT=2.997925E10
      ECHARG=4.80298E-10
      RYDBRG=1.097373E5
      STEFAN=5.675E-5
      CLIGHT=2.99793E10
      PI=3.14159265
      PI4C=PI*4./CLIGHT
C
C LOGICAL UNITS
      IREAD=5
      IWRIT=7
C
C TEMPERATURE, GRAVITATION
      READ(5,62) TEFF,G,IDRAB1,IDRAB2,IDRAB3,IDRAB4,IDRAB5,IDRAB6,RELM
      write(6,621)TEFF,G,IDRAB1,IDRAB2,IDRAB3,IDRAB4,IDRAB5,IDRAB6,RELM
      
      FLUX=STEFAN*TEFF**4/PI
      GRAV=10.**G
      RELRAD=0.5*(log10(RELM)+GSUN-G)
      RELLUM=4.*log10(TEFF/TSUN)+2.*RELRAD
      RELATM=-3.+log10(TSUN/TEFF/RELM)+0.5*RELLUM
      BOLMAG=4.72-2.5*RELLUM
      RADIUS=RSUN*10.**RELRAD
      DO 100 K=1,NDP
100   RR(K)=RADIUS
C
C KONV, MIHAL
CUGJ      ILINE=0
      ISTRAL=0
      READ(5,631) NOCONV,XMAX,TAUM,FACPLY,METBL
      write(6,632) NOCONV,XMAX,TAUM,FACPLY,METBL
      IF(METBL.GT.1) METBL=1
      IF(TAUM.EQ.0.) TAUM=50.
      IF(XMAX.EQ.0.) XMAX=1.E10
      READ(5,63) ILINE,AMLOSS,MOLTSUJI,LIN_CIA,ISPEC
      write(6,633) ILINE,AMLOSS,MOLTSUJI,LIN_CIA,ISPEC
      READ(5,634) tdiam1,tdiam2,fdiam1,fdiam2
      READ(5,634) tc2h21,tc2h22,fc2h21,fc2h22
      if (fdiam1.ne.0..or.fdiam2.ne.0..or.fc2h21
     &    .ne.1..or.fc2h22.ne.1) then
      write(7,*)' tdiam1,tdiam2,fdiam1,fdiam2 ='
      write(7,*) tdiam1,tdiam2,fdiam1,fdiam2
      write(7,*)' tc2h21,tc2h22,fc2h21,fc2h22 ='
      write(7,*) tc2h21,tc2h22,fc2h21,fc2h22
      end if
      READ(5,51) MIHAL,KONSG,KORT,KPP,TCONV
      write(6,51) MIHAL,KONSG,KORT,KPP,TCONV
      read(5,635) lops,nops
      write(6,*) ' lops = ',lops,'   nops = ',nops
C
C CONVECTION PARAMETERS
      READ(5,50) PALFA,PBETA,PNY,PY,VFIX
      write(6,50) PALFA,PBETA,PNY,PY,VFIX
C
C MY POINTS  (warning: there are other variables with the same names xmy,h...)
C                      
C
      READ(5,51) MMY,NCORE,KDIFF
      write(6,51) MMY,NCORE,KDIFF
      TAURAT=10.**(0.1*KDIFF-0.01)
C
      ACALL = 0.0D+0
      BCALL = 1.0D+0
      CALL GAUSI(MMY,ACALL,BCALL,H,XMY)
      DO 110 I=1,MMY
110   XMY2(I)=XMY(I)**2
C
C TAU SCALE
      CALL TAUSCA
C
C NTAU is number of depth points in input-model. JTAU is number of depth
C points demanded in the input-file (input tauscale)
C
      WRITE(6,*) ' NTAU,JTAU in MAINB = ',NTAU,JTAU
      IF (NTAU.NE.JTAU) CALL SCALEMOD
      WRITE(6,*) 'NTAU,JTAU at end of MAINB =',NTAU,JTAU
      DO 200 K=1,NTAU
      TAULN(K)=log(TAU(K))
200   CONTINUE
C
      WRITE(6,60)
      WRITE(6,61)TEFF,RELLUM,BOLMAG
      WRITE(6,52)G,RELATM,RELRAD
      WRITE(6,65)RELM
      WRITE(6,53)PALFA
      WRITE(6,54)PBETA
      WRITE(6,55)PNY
      WRITE(6,56)PY
      WRITE(6,59)NOCONV
      WRITE(6,57)MIHAL
      WRITE(6,571)KONSG,KORT,KPP
      WRITE(6,58)MMY
      WRITE(6,69)FACPLY
      WRITE(6,66)XMAX
      WRITE(6,67)TAUM
      WRITE(6,68)METBL
      WRITE(6,64)NCORE,KDIFF
      RETURN
50    FORMAT(5(7X,F8.0))
51    FORMAT(4(7X,I3,5X),7X,F8.0)
52    FORMAT(20X,'LOG G  =',F10.2,10X,'LOG (ATM/R) =',F5.2,10X,
     & 'LOG (R/RSUN)=',F5.2)
53    FORMAT(/20X,'PALFA  =',F10.2)
54    FORMAT(20X,'PBETA  =',F10.2)
55    FORMAT(20X,'PNY    =',F10.0)
56    FORMAT(20X,'PY     =',F10.3)
57    FORMAT(20X,'MIHAL  =',I5)
571   FORMAT(20X,'KONSG  =',I5,'  KORT  =',I5,'  KPP  =',I5)
58    FORMAT(20X,'MYPNTS =',I10,'(no use for spherical)')
59    FORMAT(/20X,'NOCONV =',I10)
60    FORMAT('0* MODEL PARAMETERS')
61    FORMAT(/20X,'TEFF   =',F10.0,10X,'LOG (L/LSUN)=',F5.2,10X,
     & 'BOLOM. MAGN.=',F6.2)
62    FORMAT(2(7X,F8.0),7X,6A4,6X,F8.0)
621   FORMAT(2(7X,F8.1),7X,6A4,6X,F8.1)
631   FORMAT(7X,I3,12X,E8.1,7X,f8.0,7x,f8.0,7x,I3)
632   FORMAT(7X,I3,12X,1pE8.1,7X,0pf8.1,7x,f8.1,7x,I3)
63    FORMAT(7X,I3,12X,E8.1,7X,I3,2(12X,I3))
633   FORMAT(7X,I3,12X,1pE8.1,3(7X,I3))
634   FORMAT(2(7X,F8.1),2(7X,f8.4))
635   FORMAT(7X,I3,12x,i3)
64    FORMAT(/20X,'NCORE  =',I10/20X,'KDIFF  =',I10)
65    FORMAT(20X,'M/MSUN =',F10.1)
66    FORMAT(20X,'XMAX   =',1PE10.2)
67    FORMAT(20X,'TAUM   =',F10.2)
68    FORMAT(20X,'METBL  =',I10)
69    FORMAT(20X,'FACPLY =',F10.3)
      END
C
      SUBROUTINE MATINV(A,N)
      implicit real*16 (a-h,o-z)
C
C 'MATINV' IS A STANDARD ROUTINE FOR MATRIX INVERSION (USED IN THE
C MIHALAS CODE). 'MATINV' ASSUMES THAT A(J,J) IS NONZERO AND MAKES
C NO PIVOTING. THIS IS SOMETIMES ADVANTAGEOUS FOR THE NUMERICAL
C ACCURACY (IN THE MARCS CODE FOR INSTANCE).
C
      include 'parameter.inc'
      DIMENSION A(NDP,NDP)
C
      IF(N.EQ.1)GOTO 25
      DO 5 I=2,N
      IM1=I-1
      DO 2 J=1,IM1
      JM1=J-1
      DIV=A(J,J)
      SUM=.0
      IF(JM1.LT.1)GOTO 2
      DO 1 L=1,JM1
    1 SUM=SUM+A(I,L)*A(L,J)
    2 A(I,J)=(A(I,J)-SUM)/DIV
      DO 4 J=I,N
      SUM=.0
      DO 3 L=1,IM1
    3 SUM=SUM+A(I,L)*A(L,J)
      A(I,J)=A(I,J)-SUM
    4 CONTINUE
    5 CONTINUE
      DO 13  II=2,N
      I=N+2-II
      IM1=I-1
      IF(IM1.LT.1)GOTO 13
      DO 12 JJ=1,IM1
      J=I-JJ
      JP1=J+1
      SUM=.0
      IF(JP1.GT.IM1)GOTO 12
      DO 11 K=JP1,IM1
   11 SUM=SUM+A(I,K)*A(K,J)
   12 A(I,J)=-A(I,J)-SUM
   13 CONTINUE
      DO 17 II=1,N
      I=N+1-II
      DIV=A(I,I)
      IP1=I+1
      IF(IP1.GT.N)GOTO 17
      DO 16 JJ=IP1,N
      J=N+IP1-JJ
      SUM=.0
      DO 15 K=IP1,J
   15 SUM=SUM+A(I,K)*A(K,J)
      A(I,J)=-SUM/DIV
   16 CONTINUE
   17 A(I,I)=1./A(I,I)
      DO 24 I=1,N
      DO 23 J=1,N
      K0=MAX0(I,J)
      IF(K0.EQ.J)GOTO 22
      SUM=.0
   20 DO 21 K=K0,N
   21 SUM=SUM+A(I,K)*A(K,J)
      GOTO 23
   22 SUM=A(I,K0)
      IF(K0.EQ.N)GOTO 23
      K0=K0+1
      GOTO 20
   23 A(I,J)=SUM
   24 CONTINUE
      RETURN
   25 A(1,1)=1./A(1,1)
      RETURN
      END
C
      SUBROUTINE MODJON(JONTYP,IOUTS)
      implicit real*16 (a-h,o-z)
C
C 'JONTYP' DETERMINES THE QUALITY OF THE IONIZATION EQUILIBRIUM COMPUTED
C IN JON.
C
C JONTYP.LE.0  A NUMBER OF ELEMENTS AND IONIZATION STAGES ARE DISREGARDE
C        EQ.1  PARTITION FUNCTIONS ARE TAKEN CONSTANT OR TEMPERATURE DEP
C              ONLY, AS SPECIFIED IN INJON-DATA. ELEMENTS AND IONIZATION
C              STAGES AS SPECIFIED IN INJON-DATA.
C       .EQ.2  FULL PARTITION FUNCTIONS, ASYMTOTIC PART ACCORDING TO
C              BASHEK ET AL.
C       .GE.3  FULL PARTITION FUNCTIONS, ASYMTOTIC PART ACCORDING TO
C              FISCHEL ET AL. ALL POSSIBLE ELEMENTS AND IONIZATION STAGE
C              ARE TAKEN INTO ACCOUNT.
C
C EACH INCREASE OF JONTYP CORRESPONDS TO AN ADDING OF THE QUALITY CITED
C ON TOP OF PREVIOS ONES.
C
C IOUTS .EQ.0  NO PRINT OUT OF NEW JON-PARAMETERS
C       .GT.0  PRINTOUT OF NEW JON-PARAMETERS
C
C ENTRY POINT 'MODMOL' ADDED 76.03.23  *NORD*
C
C MOLTYP.LE.0  GIVES FULL HANDLING OF MOLECULES
C        GT.0  GIVES H-MOLECULES ONLY
C

      COMMON /CI3/DUM(885),IDUM(215),IFISH
      COMMON /CI4/TMOLIM,IELEM(16),ION(16,5),MOLH,JUMP
      COMMON /CI6/TP,IQFIX(16,5),NQTEMP
      DIMENSION IELS(16),IONS(16,5),IQFS(16,5)
      DATA JONOLD/1/,IMAX/16/,JMAX/5/
      INTEGER MOLH, JUMP
C
C CHECK IF READY
      IF(JONTYP.EQ.JONOLD) RETURN
      IF(JONTYP.GT.3.AND.JONOLD.EQ.3) RETURN
      IF(JONTYP.LT.0.AND.JONOLD.EQ.0) RETURN
      TP=0.
      IF(JONTYP.LT.JONOLD) GO TO 60
C
C ZERO TO ONE
      IF(JONOLD.GT.0) GO TO 20
      DO 10 I=1,IMAX
      IELEM(I)=IELS(I)
      DO 10 J=1,JMAX
      ION(I,J)=IONS(I,J)
10    CONTINUE
      JONOLD=1
      IF(JONOLD.EQ.JONTYP) GO TO 90
C
C ONE TO TWO
20    IF(JONOLD.GT.1) GO TO 30
      DO 21 I=1,IMAX
      DO 21 J=1,JMAX
      IQFS(I,J)=IQFIX(I,J)
      IQFIX(I,J)=2
21    CONTINUE
      IFISH=0
      JONOLD=2
      IF(JONOLD.EQ.JONTYP) GO TO 90
C
C TWO TO THREE
30    IF(JONOLD.GT.2) GO TO 90
      IFISH=1
      DO 31 I=1,IMAX
      IELS(I)=IELEM(I)
      IELEM(I)=1
      DO 31 J=1,JMAX
      IONS(I,J)=ION(I,J)
      ION(I,J)=1
31    CONTINUE
      JONOLD=3
      GO TO 90
C
C THREE TO TWO
60    IF(JONOLD.LT.3) GO TO 70
      IFISH=0
      DO 61 I=1,IMAX
      IELEM(I)=IELS(I)
      DO 61 J=1,JMAX
      ION(I,J)=IONS(I,J)
61    CONTINUE
      JONOLD=2
      IF(JONOLD.EQ.JONTYP) GO TO 90
C
C TWO TO ONE
70    IF(JONOLD.LT.2) GO TO 80
      DO 71 I=1,IMAX
      DO 71 J=1,JMAX
      IQFIX(I,J)=IQFS(I,J)
71    CONTINUE
      JONOLD=1
      IF(JONOLD.EQ.JONTYP) GO TO 90
C
C ONE TO ZERO
80    IF(JONOLD.LT.1) GO TO 90
      DO 81 I=1,IMAX
      IELS(I)=IELEM(I)
      DO 81 J=1,JMAX
      IONS(I,J)=ION(I,J)
81    CONTINUE
      ION(2,2)=0
      ION(3,3)=0
      IELEM(4)=0
      IELEM(6)=0
      IELEM(7)=0
      ION(8,3)=0
      ION(9,3)=0
      ION(9,4)=0
      ION(10,3)=0
      IELEM(11)=0
      IELEM(12)=0
      IELEM(13)=0
      IELEM(14)=0
      IELEM(16)=0
      ION(15,3)=0
      JONOLD=0
C
C RETURN
90    CONTINUE
      IF(IOUTS.GT.0) CALL INJON(2)
      RETURN
C
C ENTRY MODMOL
      ENTRY MODMOL(MOLTYP,IOUTS)
      IF(MOLTYP.LE.0) MOLH=0
      IF(MOLTYP.GT.0) MOLH=1
      IF(IOUTS.GT.0) CALL INJON(2)
      WRITE(7,101) MOLTYP
101   FORMAT('0MOLTYP=',I1)
      RETURN
      END
C
      SUBROUTINE MOL(T,PE,G2,GC,GN,GO,ABUC,ABUO,ABUN,
     &               XIH,XKHM,XIHM,XNEN,F1,F2,F3,F4,F5)
      implicit real*16 (a-h,o-z)
      include 'parameter.inc'
C
C        THIS ROUTINE COMPUTES DISSOCIATION EQUILIBRIA FOR THE MOLECULES H2,H2+
C        H2O,OH,CH,CO,CN,C2,O2,N2,NH AND NO WITH H,H-,H+,C,C+,O,O+,N,N+ CON-
C        SIDERED, USING A NEWTON-RAPHSON SCHEME. SOME FEATURES COME FROM THE
C        MONSTER AND FROM MIHALAS, METH. COMP. PHYS. 7,1.
C
C        G2=N(HII)/N(HI), GC=N(CII)/N(CI) ETC.
C        ABUC= THE NUMBER OF CARBON NUCLEI PER HYDROGEN NUCLEUS, ABUO AND ABUN
C        ARE THE CORRESPONDING VALUES FOR OXYGEN AND NITROGEN.
C        XIH = THE IONIZATION ENERGY OF HYDROGEN
C        XIHM= THE 'DISSOCIATION ENERGY' OF H-
C        XKHM= THE 'DISSOCIATION CONSTANT' OF H-
C        XNEN= THE NUMBER OF ELECTRONS PER UNIT VOLUME FROM ELEMENTS OTHER THAN
C        HYDROGEN, CARBON, OXYGEN AND NITROGEN.
C
C        THE SUBSCRIPT IN AKD(I),PK(I)  HAS THE FOLLOWING MEANING
C        I=1 H-, 2 H2, 3 H2+, 4 H2O, 5 OH, 6 CH, 7 CO, 8 CN, 9 C2, 10 N2,
C         11 O2, 12 NO, 13 NH, 14 C2H2, 15 HCN, 16 C2H, 17 -, 18 HS
C         19 SIH, 20 C3H, 21 C3, 22 CS, 23 SIC, 24 SIC2, 25 NS
C         26 SIN, 27 SIO, 28 SO, 29 S2, 30 SIS, 31 TiO, 32 TiO2, 33 TiC2
C        WHEREAS AKA, AK1 - AK4  HAVE  I=I-1  (CF. LOOP 100)
C
C        NMOL IS THE NUMBER OF MOLECULES CONSIDERED
C
C        THIS ROUTINE CALLS MOLMAT AND AINV3
C        THE DATA FOR COMPUTING THE DISSOCIATION CONSTANTS ARE FROM TSUJI
C        (ASTRON. ASTROPHYS.  23,411 (1973))
C        WITH CORRECTIONS OF NH(TO D0=3.46 EV),
C        BG 831114.
C
C     841210 KE:  CHANGED DISS.EN. FOR CS  WITH -0.5 EV ( IN AK1 )
C     950625 UGJ: CHANGED DISS.EN. FOR TIO TO 6.89 EV (Costes 95) BY CHANGING
C                 AK1(30) FROM 8.5956 (Tsuji_73_D0=7.26eV) to 8.2256
C     950625 UGJ: CHANGED DISS.EN. FOR CN TO 7.77 EV (Costes 94) BY CHANGING
C                 AK1(6) FROM 8.2793 (Tsuji_73_D0=7.9eV) to 8.1493
C
C        THE ROUTINE GIVES FH, FC, FO, FN, FE. FH=P(HI)/PH, FC=P(CI)/PH ETC.,
C        WHERE PH=NH*KT (NH IS THE NUMBER OF HYDROGEN NUCLEI PER CM3).
C
      DOUBLE PRECISION AKA(32),AK1(32),AK2(32),AK3(32),AK4(32),AKD(33),
     & PK(33),TH,WKH,WKC,WKN,WKO,FH,FC,FN,FO,FE,FS,FK,FT,FF(8),
     & CAM,ROOT,EMAX,DIFF,A(8,8),F(8),D(8),S,R,PELLE(4),PALLE(4),
     & R1,R2,R3,R4,R5,ROT(4),ZIM(4),X02AAF
C     & CAM,ROOT,EMAX,DIFF,A(7,7),F(7),D(7),S,R,PELLE(4),PALLE(4),
C ARGUMENTS TO STATEMENT FUNCTION 880916, apollo
      DOUBLE PRECISION X,DX
      EQUIVALENCE (FF(1),FH),(FF(2),FC),(FF(3),FO),(FF(4),FN),
     & (FF(5),FE),(FF(6),FS),(FF(7),FK),(FF(8),FT)
      dimension qr(22),tqr(22),ptio(ndp),poxg1(ndp)
     &          ,xnti1(ndp),xnti2(ndp)
      DIMENSION B1(5),B2(5),DIS(10)
      COMMON /CMOL1/EH,FFE,FFH,FHE,FFC,FCE,FFN,FNE,FFO,FOE,FFK,FKE,
     &              FFS,FSE,FFT,FTE
      COMMON /CMOL2/PPK(33),NMOL
      COMMON /CI5/ abmarcs(17,ndp),ANJON(17,5),DUMT(99)
      COMMON /CMOLRAT/ FOLD(NDP,8),MOLOLD,KL
      data tqr/1000., 1500., 2000., 2500., 3000., 3500., 4000., 4500.,
     &         5000., 5500., 6000., 6500., 7000., 8000., 9000., 10000.,
     &        12000., 14000., 16000., 18000., 20000., 50000./
      data qr/1.708, 1.911, 2.035, 2.097, 2.112, 2.091, 2.042, 1.972,
     &        1.887, 1.789, 1.684, 1.574, 1.462, 1.241, 1.032, 0.884,
     &        0.541, 0.330, 0.194, 0.111, 0.062, 0.001/
      DATA B1/2.6757,1.4772,0.60602,0.12427,0.00975/,
     *B2/2.9216,2.0036,1.7231,0.82685,0.15253/,
     *DIS/9.50,4.38,3.47,11.11,7.90,6.12,9.76,5.12,6.51,3.21/
      DATA AKA/12.739D0,11.206998D0,25.420D0,12.371D0,12.135D0,13.820D0,
     *12.805D0,12.804D0,13.590D0,13.228D0,12.831D0,12.033D0,
     +38.184D0,25.635D0,25.063D0,
     +0.0D0,12.019D0,11.852D0,40.791D0,25.230D0,13.436D0,12.327D0,
     +25.623D0,12.543D0,12.399D0,13.413D0,12.929D0,12.960D0,
     +13.182D0,13.398D0,27.901D0,27.018D0/
      DATA AK1/5.1172D0,2.7942767D0,10.522D0,5.0578D0,4.0760D0,11.795D0,
     *8.1493D0,6.5178D0,10.585D0,5.5181D0,7.1964D0,4.0935D0,
     +17.365D0,13.833D0,12.291D0,
     +0.0D0,4.2922D0,3.7418D0,21.762D0,14.445D0,8.0574D0,5.0419D0,
     +13.085D0,5.9563D0,5.4876D0,8.8710D0,6.0100D0,5.0952D0,
     +7.1147D0,8.2256D0,14.031D0,13.534D0/
      DATA AK2/1.2572D-1,-7.9196803D-2,1.6939D-1,1.3822D-1,1.2768D-1,
     *1.7217D-1,6.4162D-2,9.7719D-2,2.2067D-1,6.9935D-2,1.7349D-1,
     +1.3629D-1,2.1512D-2,1.3827D-1,-1.9036D-2,
     +0.0000D+0,1.4913D-1,1.5999D-1,9.3377D-1,1.2547D-1,1.8754D-1,
     +1.3941D-1,-5.5227D-2,2.0901D-1,9.5301D-2,1.5042D-1,1.6253D-1,
     +1.8027D-1,1.9300D-1,0.40873D0,0.42156D0,0.45875D0/
      DATA AK3/1.4149D-2,-2.4790744D-2,1.8368D-2,1.6547D-2,1.5473D-2,
     *2.2888D-2,7.3627D-3,1.2739D-2,2.9997D-2,8.1511D-3,2.3065D-2,
     +1.6643D-2,8.8961D-5,1.8122D-2,-4.4498D-3,
     +0.0000D+0,1.8666D-2,2.0629D-2,1.3863D-1,1.7390D-2,2.5507D-2,
     +1.9363D-2,-9.3363D-3,2.8986D-2,1.3369D-2,1.9581D-2,2.1665D-2,
     +2.4324D-2,2.5826D-2,5.7937D-2,6.1271D-2,6.6158D-2/
      DATA AK4/6.3021D-4,0.D0,8.1730D-4,7.7224D-4,7.2661D-4,1.1349D-3,
     *3.4666D-4,6.2603D-4,1.4993D-3,3.7970D-4,1.1380D-3,7.8691D-4,
     +-2.8720D-5,9.1645D-4,-2.3073D-4,
     +0.0000D+0,8.9438D-4,9.9897D-4,7.4549D-3,8.8394D-4,1.2735D-3,
     +9.6202D-4,-4.9876D-4,1.4621D-3,6.9396D-4,9.4828D-4,1.0676D-3,
     +1.2049D-3,1.2648D-3,2.9287D-3,3.1476D-3,3.3834D-3/
C
C STATEMENT FUNCTION FOR CORRECTIONS LIMITED TO FACTOR TWO
      ASUM(X,DX)=DMIN1(2.D0*X,DMAX1(0.5D0*X,X+DX))
C      write(7,*) ' 1,molold = ',molold
C      write(7,*) ' MOLOLD,KL,FOLD(1-NDP,1-8) = '
C      write(7,*) MOLOLD,KL,((FOLD(i,j),i=1,ndp),j=1,8)
C
C COMPUTATION OF DISSOCIATION CONSTANTS K(AB) (AKD) AND PE/K(AB) (PK).
      TMEM=T
C.... IF (T.LT.1200.)T=1200.
C                              THIS MAY BE NECESSARY ON SOME MACHINES
C
C For CH4 it is necessary to do the Pe/K(CH4) in logarithm, because
C K(CH4) can be bigger than 39. For the other molecules it is
C practical too; therefore this change 27-3-89/  UGJ.
C
C
C      TH=5040./T
C      AKD(1)=LOG10(XKHM)
C      PK(1)=PE/XKHM
C      PELOG=LOG(PE)
C      DO 100 J=2,33
C        M=J-1
C        AKDJ=(AKA(M)-(AK1(M)-(AK2(M)-(AK3(M)-AK4(M)*TH)*TH)*
C     &   TH)*TH)
C        PK(J)=PELOG-AKDJ
C100   CONTINUE
C
C
c      if (kl/10*10.eq.kl) then
c      write(6,*) ' mol: kl,T,PE,G2,GC,GN,GO,ABUC,ABUO,ABUN,',
c     &              ' XIH,XKHM,XIHM,XNEN,F1,F2,F3,F4,F5 = '
c      write(6,*) kl,T,PE,G2,GC,GN,GO,ABUC,ABUO,ABUN,
c     &               XIH,XKHM,XIHM,XNEN,F1,F2,F3,F4,F5
c      end if

      TH=5040.D0/T
      AKD(1)=XKHM
      PK(1)=PE/XKHM
      PELOG=LOG10(PE)
      DO 100 J=2,33
        M=J-1
C        AKD(J)=10.D0**(AKA(M)-(AK1(M)-(AK2(M)-(AK3(M)-AK4(M)*TH)*TH)*
C     &   TH)*TH)
C        PK(J)=PE/AKD(J)
        AKD(J)=AKA(M)-(AK1(M)-(AK2(M)-(AK3(M)-AK4(M)*TH)*TH)*TH)*TH
        PK(J)=PELOG-AKD(J)
100   CONTINUE
C      PK(4)=PE*PK(4)
C      PK(14)=PE*PE*PK(14)
C      PK(15)=PE*PK(15)
C      PK(16)=PE*PK(16)
C      PK(20)=PE*PE*PK(20)
C      PK(21)=PE*PK(21)
C      PK(24)=PE*PK(24)
C
      PK(4)=PELOG+PK(4)
      PK(14)=2.0*PELOG+PK(14)
      PK(15)=PELOG+PK(15)
      PK(16)=PELOG+PK(16)
      PK(20)=2.0*PELOG+PK(20)
      PK(21)=PELOG+PK(21)
      PK(24)=PELOG+PK(24)
      PK(32)=PELOG+PK(32)
      PK(33)=PELOG+PK(33)
      DO 101 J=2,33
101   PK(J)=10.D0**PK(J)
C
C      write(7,*) ' now in mol'
C      write(7,*) (pk(i),i=1,33)
C
C      ANJON(17,2) = 0.
C      ANJON(17,1) = 1.
C
        thta=5040./t
*
* find the ratio  u(Ti II)/u(Ti I)  by interpolating in array qr
*
        j=1
        do 662 jj=1,22
          if(t.le.tqr(jj)) goto 663
          j=jj
662     continue
663     uratio=qr(j) + (qr(j+1)-qr(j))/(tqr(j+1)-tqr(j))*(t-tqr(j))
*
* q1 is the ratio  n(Ti II)/n(Ti I)
*
        q1=.6667*t**2.5*uratio*exp(-6.82*11605./t)/pe
*
*
* q2 is the ratio  n(TiO)/n(Ti I)
*
         tiokp = 10.**akd(31)
C        q2=poxg1(i)/tiokp
*
* xnti1 is n(Ti I)    (i.e. per cm**3)
*
C        xnti1(i)=epsti*ro(i)*div/(1.+q1+q2)
C        ptio(i)=1.38e-16*t*xnti1(i)*q2
C        xnti2(i)=xnti1(i)*q1
         ptio(kl) =  0.
         poxg1(kl) = 0.
         xnti1(kl) = 0.
         xnti2(kl) = 0.
C
      ANJON(17,2) = 1./(1.+1./q1)
      ANJON(17,1) = 1./(1.+q1)
C
C        if (kl/10*10.eq.kl .or.kl.eq.1) then
C         write(6,664) kl,t,pe,q1,q2,anjon(17,1),anjon(17,2)
C     &        ,tiokp
C     &        ,log10( max(ptio(kl),1.d-99) )
C     &        ,poxg1(kl),tiokp,xnti1(kl),xnti2(kl)
C664     format(' MOL:',i3,f7.0,1p6e9.2 )
C        end if
C      write(7,*) ' abund(1-17): '
C      write(7,654) (abund(i),i=1,17)
C      write(7,*) ' anjon(1-17,1) -neutral: '
C      write(7,654) (anjon(i,1),i=1,17)
C      write(7,*) ' anjon(1-17,2) -1.ionis: '
C      write(7,654) (anjon(i,2),i=1,17)
C654   format(1p6e10.3)
C
C COMPUTATION OF STARTING VALUES FOR FH,FC ETC.
      XNENSK=XNEN-abmarcs(11,kl)*ANJON(11,2)-abmarcs(10,kl)*ANJON(10,2)
      GS=ANJON(11,2)/ANJON(11,1)
      GK=ANJON(10,2)/ANJON(10,1)
      GT=ANJON(17,2)/ANJON(17,1)
      AS=abmarcs(11,kl)/abmarcs(1,kl)
      AK=abmarcs(10,kl)/abmarcs(1,kl)
      AT=abmarcs(17,kl)/abmarcs(1,kl)
C
C START WITH VALUES OF FH,FC,FN,FO,FS,FK,FE FROM LAST ITTERATION
C      write(7,*) ' molold = ',molold
C      write(7,*) ' we now fix it to molold == 1'
C      molold = 1
C      write(7,*) ' still in mol: t,pe,g2,... ='
C      write(7,*) T,PE,G2,GC,GN,GO,ABUC,ABUO,ABUN,
C     &               XIH,XKHM,XIHM,XNEN,F1,F2,F3,F4,F5
      IF (MOLOLD.EQ.1) THEN
         FH = FOLD(KL,1)
         FC = FOLD(KL,2)
         FN = FOLD(KL,3)
         FO = FOLD(KL,4)
         FS = FOLD(KL,5)
         FK = FOLD(KL,6)
         FE = FOLD(KL,7)
         FT = FOLD(KL,8)
C         write(7,*)' input in mol:kl,fold(kl,1-8)= ',
C     *             kl,(fold(kl,i),i=1,8)
         GO TO 159
      END IF
C
C ...OR ESTIMATE THESE VALUES...
      WKH=1.D0+G2+PK(1)
      WKC=1.D0+GC
      WKO=1.D0+GO
      WKN=1.D0+GN
C
C FH, FROM H=HI+HII+H2+H-
      FE=XNENSK
      CAM=FE*WKH/(2.D0*PK(2))
      ROOT=SQRT(ABS(CAM*CAM+FE/PK(2)))
      FH=-CAM+ROOT
      FE=FH*(G2-PK(1))+XNENSK
      IF (FE.GT.0.D0) GO TO 110
      R=PK(2)/PK(1)/PK(1)
      CAM=(WKH*XNENSK/PK(1)-1.D0-2.D0*XNENSK*R)
      S=R-WKH/PK(1)
      CAM=CAM/2.D0/S
      ROOT=SQRT(CAM*CAM-R*XNENSK*XNENSK/S)
      FE=-CAM-ROOT
      IF (FE.LT.0.) FE=-CAM+ROOT
      FH=(XNENSK-FE)/PK(1)
110   CONTINUE
C
C FN, FROM N=NI+NII+N2
      R1=PK(10)/FE
      R2=WKN
      R3=-abmarcs(4,kl)
      FN=(-R2+SQRT(R2*R2-4.D0*R1*R3))/(2.D0*R1)
C
C FC, FROM C=CI+CII+CO+C2+C2H2+HCN+CN AND O=OI+OII+CO+H2O
      NDEG =3
      IF (abmarcs(5,kl).GE.abmarcs(3,kl)) GO TO 130
      R1=PK(7)/(FE*WKO)
      R2=2.D0*FH*FH*PK(14)/(FE*FE*FE)+2.D0*PK(9)/FE
      R3=WKC+FN*PK(8)/FE+FH*FN*PK(15)/(FE*FE)
      R4=WKO+FH*FH*PK(4)/(FE*FE)
C                                  R4 = WATER CORRECTION, 4-JAN-82/NORDLUND
      R1LOG=LOG10(R1)
      R2LOG=LOG10(R2)
      SUMLOG=R1LOG+R2LOG
      IF (SUMLOG.GE.34.) THEN
         RADD=(34.-SUMLOG)/2.
c         WRITE(6,*) 'R1,R2,RADD =',R1,R2,RADD
         R1=R1*10.**RADD
         R2=R2*10.**RADD
      END IF
      PELLE(1)=R1*R2
      R1LOG=LOG10(R1)
      R2LOG=LOG10(R2)
      R3LOG=LOG10(R3)
      R4LOG=LOG10(R4)
      SUMLOG=R1LOG+R2LOG+R3LOG+R4LOG
      IF (SUMLOG.GE.34.) THEN
         RADD=(34.-SUMLOG)/4.
C         WRITE(6,*) 'R1,R2,R3,R4,RADD =',R1,R2,R3,R4,RADD
         R1=R1*10.**RADD
         R2=R2*10.**RADD
         R3=R3*10.**RADD
         R4=R4*10.**RADD
      END IF
      PELLE(2)=R2*R4+R1*R3
      PELLE(3)=R3*R4+R1*(abmarcs(5,kl)-abmarcs(3,kl))
      PELLE(4)=-abmarcs(3,kl)*R4
C                                  R4 = WATER CORRECTION, 4-JAN-82/NORDLUND
*
* Original PELLE calculations
*
*      PELLE(1)=R1*R2
*      PELLE(2)=R2*R4+R1*R3
*      PELLE(3)=R3*R4+R1*(ABUND(5)-ABUND(3))
*      PELLE(4)=-ABUND(3)*R4
*
* The checks below are intended to stop machines with a small exponent
* range (i.e., VAX, Alliant etc.) from evaluating a result which will
* cause an overflow exception. Since Mol doesn't contribute much to
* the total execution time, the overhead is negligible.
*
C      if (dlog10(dabs(r1))+dlog10(dabs(r2)).gt.36.) then
C        pelle(1)=1.e+36*(r1/dabs(r1))*(r2/dabs(r2))
C      else
C        pelle(1)=r1*r2
C      end if
*
C      if (dlog10(dabs(r2))+dlog10(dabs(r4)).gt.35. .or. 
C     +    dlog10(dabs(r1))+dlog10(dabs(r3)).gt.35.) then
C        pelle(2)=1.e+36
C      else
C        pelle(2)=r2*r4+r1*r3
C      end if
*
C      if (dlog10(dabs(r3))+dlog10(dabs(r4)).gt.35. .or.
C     +    dlog10(dabs(r1))+log10(abs(abund(5)-abund(3))).gt.35.)
C     +    then
C        pelle(3)=1.e+36
C      else
C        pelle(3)=r3*r4+r1*(abund(5)-abund(3))
C      end if
*
C      if (dlog10(dabs(r4))+log10(abs(abund(3))).gt.36.) then
C        pelle(4)=1.e+36
C      else
C        pelle(4)=-abund(3)*r4
C      end if
*C.... PRINT 120,PELLE
120   FORMAT(1X,5D10.3)
      DO 121 IP=1,4
        PALLE(IP)=PELLE(IP)
121   CONTINUE
      NPOL=NDEG+1
      IFAIL=0
      CALL C02AEF(PELLE,NPOL,ROT,ZIM,X02AAF(DUM),IFAIL)
      IFLAG=0
      DO 122 III=1,3
        IF (ROT(III).LE.0.D0.OR.ROT(III).GT.abmarcs(3,kl)) GO TO 122
        FC=ROT(III)
        IFLAG=IFLAG+1
122   CONTINUE
      FO=abmarcs(5,kl)/(R4+FC*R1)
C.... PRINT 120,FC**3*PALLE(1),FC**2*PALLE(2),FC*PALLE(3),PALLE(4)
      IF (IFLAG.EQ.1) GO TO 135
      IF (IFLAG.EQ.0) PRINT 123
      IF (IFLAG.GE.2) PRINT 124
123   FORMAT('0MOL: NO ROOT FOUND FOR CARBON INITIAL RATIO FC')
124   FORMAT('0MOL: SEVERAL ROOTS FOR CARBON INITIAL RATIO FC')
      STOP ' stop at loop 122 in mol: root problem in mol '
C
C FO, FROM O=OI+OII+CO+H2O+SIO, C=CI+CII+CO+CN+HCN, AND SI=SI(I)+SI(II)+SIO
130   CONTINUE
      R1=PK(7)/FE
      R2=PK(27)/FE
      R3=WKC+FN*PK(8)/FE+FH*FN*PK(15)/(FE*FE)
      R4=1.D0+GK
      R5=WKO+FH*FH*PK(4)/(FE*FE)
      R1LOG=LOG10(R1)
      R2LOG=LOG10(R2)
      R3LOG=LOG10(R3)
      R4LOG=LOG10(R4)
      R5LOG=LOG10(R5)
      SUMLOG=R1LOG+R2LOG+R5LOG
      IF (SUMLOG.GE.34.) THEN
         RADD=(34.-SUMLOG)/3.
C         WRITE(6,*) 'R1,R2,R5,RADD =',R1,R2,R5,RADD
         R1=R1*10.**RADD
         R2=R2*10.**RADD
         R5=R5*10.**RADD
      END IF
      PALLE(1)=R1*R2*R5
      SUMLOG=MAX(R1LOG+R4LOG,R2LOG+R3LOG)+R5LOG
      IF (SUMLOG.GE.34.) THEN
         RADD=(34.-SUMLOG)/3.
C         WRITE(6,*) 'R1,R2,R3,R4,R5 =',R1,R2,R3,R4,R5
         R1=R1*10.**RADD
         R4=R4*10.**RADD
         R2=R2*10.**RADD
         R3=R3*10.**RADD
         R5=R5*10.**RADD
      END IF
      PALLE(2)=(R1*R4+R2*R3)*R5+R1*R2*
     *    (abmarcs(3,kl)+abmarcs(10,kl)-abmarcs(5,kl))
      SUMLOG=R3LOG+R4LOG+R5LOG
      IF (SUMLOG.GE.34.) THEN
         RADD=(34.-SUMLOG)/3.
         WRITE(6,*) 'R3,R4,R5 =',R3,R4,R5
         R3=R3*10.**RADD
         R4=R4*10.**RADD
         R5=R5*10.**RADD
      END IF
      PALLE(3)=R5*R4*R3-(R1*R4+R2*R3)*abmarcs(5,kl)+R1*R4*abmarcs(3,kl)
     & +R2*R3*abmarcs(10,kl)
      PALLE(4)=-R3*R4*abmarcs(5,kl)
      DO 131 IP=1,4
        PELLE(IP)=PALLE(IP)
131   CONTINUE
C.... PRINT 120,PELLE
      NPOL=NDEG+1
      IFAIL=0
      CALL C02AEF(PELLE,NPOL,ROT,ZIM,X02AAF(DUM),IFAIL)
      IFLAG=0
      DO 132 III=1,3
        IF (ROT(III).LE.0.D0.OR.ROT(III).GT.abmarcs(5,kl)) GO TO 132
        FO=ROT(III)
        IFLAG=IFLAG+1
132   CONTINUE
      FC=abmarcs(3,kl)/(WKC+R1*FO+FN*PK(8)/FE+FN*FH*PK(15)/(FE*FE))
C.... PRINT 120,FO**3*PALLE(1),FO**2*PALLE(2),FO*PALLE(3),PALLE(4)
      IF (IFLAG.EQ.1) GO TO 135
      IF (IFLAG.EQ.0) PRINT 133
      IF (IFLAG.GE.2) PRINT 134
133   FORMAT('0MOL: NO ROOT FOUND FOR OXYGEN INITIAL RATIO FO')
134   FORMAT('0MOL: SEVERAL ROOTS FOR OXYGEN INITIAL RATIO FO')
      STOP ' stop at loop 132 in mol;  root problems in mol '
135   CONTINUE
C
C FK, FROM SI=SI(I)+SI(II)+SIS+SIH+SIC2+SIO, AND S=S(I)+S(II)+CS+SIS
      R1=1.0D0+GS+FC*PK(22)/FE
      R2=PK(30)/FE
      R3=1.D0+GK+FH*PK(19)/FE+FC*FC*PK(24)/FE
      R4=PK(27)/FE
      R5=1.0D0+GO+FC*PK(7)/FE
      ALFA=abmarcs(10,kl)-abmarcs(11,kl)-abmarcs(5,kl)
      PALLE(1)=R2*R3*R4
      PALLE(2)=R3*(R1*R4+R2*R5)-ALFA*R2*R4
      PALLE(3)=R1*R3*R5-ALFA*(R1*R4+R2*R5)-abmarcs(11,kl)*R1*R4-
     &  abmarcs(5,kl)*R2*R5
      PALLE(4)=-R1*R5*abmarcs(10,kl)
      NDEG=3
      IFAIL=0
      NPOL=NDEG+1
      CALL C02AEF(PALLE,NPOL,ROT,ZIM,X02AAF(DUM),IFAIL)
      IFLAG=0
      DO 140 III=1,3
        IF (ROT(III).LE.0.D0.OR.ROT(III).GT.abmarcs(10,kl)) GO TO 140
        IFLAG=IFLAG+1
        FK=ROT(III)
140   CONTINUE
C
C FS, FROM S=S(I)+S(II)+SIS+CS
      FS=abmarcs(11,kl)/(R1+FK*R2)
      IF (IFLAG.EQ.1) GO TO 152
      IF (IFLAG.LE.0) PRINT 150
      IF (IFLAG.GT.1) PRINT 151
150   FORMAT('0MOL: NO ROOT FOUND FOR SILICON INITIAL RATIO FK')
151   FORMAT('0MOL: SEVERAL ROOTS FOR SILICON INITIAL RATIO FK')
      STOP ' stop at loop 140;  root problems in mol '
152   CONTINUE
C
159   CONTINUE  
C  if use of old FH,FC,.. go directly here
C
C NEWTON-RAPHSON IMPROVEMNT USING ALL RELEVANT MOLECULES, ATOMS AND IONS.
C DIFF GIVES THE APPROX. ACCURACY TO WHICH FH ETC. HAVE TO CONVERGE
C BEFORE THE ITERATIONS ARE STOPPED.
C
      DIFF=1.D-3
       M=8
C       M=7
C      if (kl.eq.10) write(7,160) FH,FC,FN,FO,FS,FK,FT,FE
C160   FORMAT('P(HI)/P(H),P(CI)/P(H)...=FH,C,N,O,S,Si,Ti,E='/,1P8E9.2)
C      write(7,*) ' now in mol'
C      write(7,*) (pk(i),i=1,33)
C      write(7,*) ' FH,FC,FN,FO,FS,FK,FT,FE = ',FH,FC,FN,FO,FS,FK,FT,FE
C      write(7,*) ' f,a = ',f,a
      DO 163 J=1,500
        CALL MOLMAT(PK,G2,GC,GN,GO,GS,GK,GT,ABUC,ABUN,ABUO
     &             ,AS,AK,AT,FH,FC,FN,FO,FS,FK,FT,FE,XNENSK,F,A)
C        CALL MOLMAT(PK,G2,GC,GN,GO,GS,GK,ABUC,ABUN,ABUO
C     &             ,AS,AK,FH,FC,FN,FO,FS,FK,FE,XNENSK,F,A)
        CALL AINV3(A,M)
        EMAX=0.D0
        DO 162 L=1,M
          D(L)=0.D0
          DO 161 LL=1,M
            D(L)=D(L)-A(L,LL)*F(LL)
161       CONTINUE
          FF(L)=ASUM(FF(L),D(L))
          FF(L) = max(FF(L),1.d-99)
          if (D(L).le.1.d-99) go to 162
          EMAX=DMAX1(ABS(D(L)/FF(L)),EMAX)
162     CONTINUE
Ctemp UGJ 020308:   IF (EMAX.LT.DIFF) GO TO 170
163   CONTINUE
        GO TO 170
      write(7,*)' MOL: THE DESIRED ACCURACY WAS NOT ACHIEVED AFTER'
      WRITE(7,*)' 500 ITERATIONS. LAST VALUES AND DIFFERENCES WERE'
      WRITE(7,164) FF
      WRITE(7,164) D
164   FORMAT(1P8D10.3)
      STOP ' stop at molecular equilibrium iteration in mol '
C
C COMPUTATION OF THE INNER ENERGY. DEH2 AND DEH2P ARE THE SUM OF
C DISSOCIATION, ROTATION AND VIBRATION ENERGIES (IN EV PER MOLECULE) FOR
C H2 AND H2+. DIS(I) IS THE DISSOCIATION ENERGY FOR THE MOLECULE (I+3)
C IN THE LIST OF MOLECULES (VALUES ARE FROM TSUJI). FOR THESE MOLECULES
C THE ROTATION AND VIBRATION ENERGIES ARE NEGLECTED.
C
170   CONTINUE
c
C      write(7,168) FH,FC,FN,FO,FS,FK,FT,FE
C168   FORMAT(' After molmath call: FH,C,N,O,S,Si,Ti,E='/,1P8E9.2)
c      PRINT 160,FH,FC,FN,FO,FS,FK,FT,FE,J
c
   98 TETA=TH
      DEH2=(B1(1)-(B1(2)-(B1(3)-(B1(4)-B1(5)*TETA)*TETA)*TETA)*TETA)*
     *8.617E-5*T-4.476
      DEH2P=(B2(1)-(B2(2)-(B2(3)-(B2(4)-B2(5)*TETA)*TETA)*TETA)*TETA)*
     *8.617E-5*T-2.648
      FHE=FH/FE
      FCE=FC/FE
      FOE=FO/FE
      FNE=FN/FE
      FKE=FK/FE
      FSE=FS/FE
      FTE=FT/FE
      EH2=(-2.*XIH+DEH2)*FHE*FH*PK(2)
      EH2P=(DEH2P-XIH)*FHE*FH*G2*PK(3)
      EHM=-(XIHM+XIH)*FH*PK(1)
      EHJ=-XIH*FH
      EH2O=-(2.*XIH+DIS(1))*FHE*FH*FO*PK(4)
      EOH=-(XIH+DIS(2))*FOE*FH*PK(5)
      ECH=-(XIH+DIS(3))*FHE*FC*PK(6)
      ECO=-DIS(4)*FCE*FO*PK(7)
      ECN=-DIS(5)*FCE*FN*PK(8)
      EC2=-DIS(6)*FCE*FC*PK(9)
      EN2=-DIS(7)*FNE*FN*PK(10)
      EO2=-DIS(8)*FOE*FO*PK(11)
      ENO=-DIS(9)*FNE*FO*PK(12)
      ENH=-(DIS(10)+XIH)*FNE*FH*PK(13)
      EH=EH2+EH2P+EHM+EHJ+EH2O+EOH+ECH+ECO+ECN+EC2+EN2+EO2+ENO+ENH
C                            NOTE THAT ENERGIES INCLUDE ONLY MOLECULES 1-13
C
C PICK UP SINGLE PRECISION VALUES
      FFH=FH
      FFC=FC
      FFN=FN
      FFO=FO
      FFK=FK
      FFS=FS
      FFT=FT
      FFE=FE
c      PRINT*,' FFO,FFE,FOE ',FFO,FFE,FOE
      NMOL=33
      DO 180 I=1,NMOL
        PPK(I)=PK(I)
180   CONTINUE
      F1=FH
      F2=G2*FH
      F3=FH*PK(1)
      F4=FH*FHE*G2*PK(3)
      F5=FH*FHE*PK(2)
      T=TMEM
C
      FOLD(KL,1) = FH
      FOLD(KL,2) = FC
      FOLD(KL,3) = FN
      FOLD(KL,4) = FO
      FOLD(KL,5) = FS
      FOLD(KL,6) = FK
      FOLD(KL,7) = FE
      FOLD(KL,8) = FT
C
C      if (kl.eq.1 .or. kl/10*10.eq.kl) then
C        write(6,661) kl,akd(31),pk(31),ft,fe,fte
C      end if
C661   format (' MOL:k,akd,pk,ft,fe,fte:',i3,5e9.2)
C
      RETURN
      END
C
      SUBROUTINE MOLEQ(T,PE,G2,XIH,XKHM,XIHM,XNENH,F1,F2,F3,
     &                    F4,F5,FE,FSUM,EH)
      implicit real*16 (a-h,o-z)
C
C        THIS ROUTINE COMPUTES DISSOCIATION EQUILIBRIA FOR THE MOLECULES
C        H2 AND H2+ WITH H+, H AND H- CONSIDERED. IT MAINLY FOLLOWS MIHA
C        METH. COMP. PHYS. 7, 1 (1967).
C
C        THE INNER ENERGY OF THE HYDROGEN GAS, EH, IS ALSO EVALUATED.
C
C        XIH=THE IONIZATION ENERGY OF HYDROGEN
C        XKHM=THE 'DISSOCIATION CONSTANT' OF H-
C        XIHM=THE 'DISSOCIATION ENERGY' OF H-.
C        XNENH=THE NUMBER OF ELECTRONS PER UNIT VOLUME FROM ELEMENTS OTH
C        HYDROGEN (Q IN MIHALAS'S ARTICLE)
C        G2,F1,F2 ETC. SEE REF.
C
C
      COMMON/UTPUT/IREAD,IWRIT
C
C        CALL MOLFYS FOR PHYSICAL DATA
      CALL MOLFYS(T,XKH2,XKH2P,DEH2,DEH2P)
C
C        CALCULATION OF THE EQUILIBRIUM
      G3=PE/XKHM
      G4=PE/XKH2P
      G5=PE/XKH2
      A=1.+G2+G3
      E=G2*G4/G5
      B=2.*(1.+E)
      C=G5
      D=G2-G3
      C1=C*B*B+A*D*B-E*A*A
      C2=2.*A*E-D*B+A*B*XNENH
      C3=-(E+B*XNENH)
      CAM=C2/(2.*C1)
      ROOT=SQRT(CAM*CAM-C3/C1)
      F1D=-CAM+ROOT
      IF(F1D.GT.1.0)F1D=-CAM-ROOT
      F5D=(1.0-A*F1D)/B
      F4D=E*F5D
      F3D=G3*F1D
      F2D=G2*F1D
      FED=F2D-F3D+F4D+XNENH
      FSUMD=F1D+F2D+F3D+F4D+F5D
      F1=F1D
      F2=F2D
      F3=F3D
      F4=F4D
      F5=F5D
      FE=FED
      FSUM=FSUMD
C
C        CALCULATION OF THE ENERGIES
      EH2=(-2.*XIH+DEH2)*F5
      EH2P=(-XIH+DEH2P)*F4
      EHM=-(XIHM+XIH)*F3
      EHJ=-XIH*F1
      EH=EHJ+EHM+EH2+EH2P
    1 CONTINUE
      RETURN
      END
C
      SUBROUTINE MOLFYS(T,XKH2,XKH2P,DEH2,DEH2P)
      implicit real*16 (a-h,o-z)
C
C        THIS ROUTINE GIVES DISSOCIATION CONSTANTS XKH2 (=N(H I)*N(H I)/(NH2))
C        AND XKH2P (=N(H I)*N(H II)/N(H2+)), EXPRESSED IN NUMBER PER CM3 AND
C        THE SUM OF DISSOCIATION, ROTATION AND VIBRATION ENERGIES, DEH2 AND
C        DEH2P FOR H2 AND H2+, RESPECTIVELY (EXPRESSED IN ERGS PER MOLECULE)
C        THE DATA ARE FROM VARDYA, M.N.R.A.S. 129, 205 (1965) AND EARLIER
C        REFERENCES. THE DISSOCIATION CONSTANT FOR H2 IS FROM TSUJI,
C        ASTRON. ASTROPHYS. 1973.
C
      DIMENSION A1(5),A2(4),B1(5),B2(5),TE(5)
      DATA A1/12.739,-5.1172,1.2572E-1,-1.4149E-2,6.3021E-4/,
     *A2/11.20699 ,-2.794276 ,-0.079196   ,0.024790   /,
     *B1/2.6757,-1.4772,0.60602,-0.12427,0.009750/,
     *B2/2.9216,-2.0036,1.7231,-0.82685,0.15253/

      TEX=5040./T
      TE(1)=1.
      DO1 K=1,4
    1 TE(K+1)=TE(K)*TEX
      XKH2=0.
      XKH2P=0.
      DEH2=0.
      DEH2P=0.
      DO2 K=1,4
      XKH2=A1(K)*TE(K)+XKH2
      XKH2P=A2(K)*TE(K)+XKH2P
      DEH2=B1(K)*TE(K)+DEH2
    2 DEH2P=B2(K)*TE(K)+DEH2P
      XKH2=A1(5)*TE(5)+XKH2
      DEH2=(B1(5)*TE(5)+DEH2)*8.617E-5*T-4.476
      DEH2P=(B2(5)*TE(5)+DEH2P)*8.617E-5*T-2.648
      XKH2=10.**XKH2
      XKH2P=10.**XKH2P
      RETURN
      END
C
      SUBROUTINE MOLMAT(PK,GGH,GGC,GGN,GGO,GGS,GGK,GGT,AAC,AAN,AAO,
     &               AAS,AAK,AAT,FH,FC,FN,FO,FS,FK,FT,FE,XYNEN,F,A)
      implicit real*16 (a-h,o-z)
C
C     DOUBLE PRECISION PK,FH,FC,FN,FO,FS,FK,FT,FE,F,A
C     DOUBLE PRECISION FHE,FCE,FNE,FOE,FSE,FKE,FTE,H,HH,C,CC,CCC,XN,
C    &                 XNN,O,OO,S,SS,XK,T,GH,GC,GN,GO,GS,GK,GT,XNEN,
C    &                 AC,AN,AO,AS,AK,AT
C
      DIMENSION PK(33),F(8),A(8,8)
CCC
C CONVERT TO DOUBLE PRECISION
C
Cugj test:
C      write(6,*) ' (pk(i),i=1,33):'
C      write(6,*) (pk(i),i=1,33)
C      write(6,*)
C      write(6,*) ' FH,FC,FN,FO,FS,FK,FT,FE = ',FH,FC,FN,FO,FS,FK,FT,FE
C      write(6,*) ' f(1),a(1,1) = ',f(1),a(1,1)
C      write(6,*) ' f(8),a(8,8) = ',f(8),a(8,8)
C     GH=DBLE(GGH)
C     GC=DBLE(GGC)
C     GN=DBLE(GGN)
C     GO=DBLE(GGO)
C     GS=DBLE(GGS)
C     GK=DBLE(GGK)
C     GT=DBLE(GGT)
C     XNEN=DBLE(XYNEN)
C     AC=DBLE(AAC)
C     AN=DBLE(AAN)
C     AO=DBLE(AAO)
C     AS=DBLE(AAS)
C     AK=DBLE(AAK)
C     AT=DBLE(AAT)
      GH=GGH
      GC=GGC
      GN=GGN
      GO=GGO
      GS=GGS
      GK=GGK
      GT=GGT
      XNEN=XYNEN
      AC=AAC
      AN=AAN
      AO=AAO
      AS=AAS
      AK=AAK
      AT=AAT
CCC
      FHE=FH/FE
      FCE=FC/FE
      FNE=FN/FE
      FOE=FO/FE
      FSE=FS/FE
      FKE=FK/FE
      FTE=FT/FE
      H=1.0+GH+PK(1)+FCE*PK(6)+FNE*PK(13)+FOE*PK(5)+FSE*PK(18)
     *  +FCE*FNE*PK(15)+FCE*FCE*(PK(16)+FCE*PK(20))+FKE*PK(19)
      HH=2.0*FHE*(PK(2)+GH*PK(3)+FCE*FCE*PK(14)+FOE*PK(4))
      C=1.0+GC+FHE*PK(6)+FNE*PK(8)+FOE*PK(7)+FSE*PK(22)+FKE*PK(23)+FHE*
     *  FNE*PK(15)
      CC=2.0*FCE*(PK(9)+FHE*FHE*PK(14)+FHE*PK(16)+FKE*PK(24)+FTE*PK(33))
      CCC=3.0*FCE*FCE*(PK(21)+FHE*PK(20))
      XN=1.0+GN+FHE*PK(13)+FCE*PK(8)+FOE*PK(12)+FSE*PK(25)+FKE*PK(26)+
     *   FHE*FCE*PK(15)
      XNN=2.0*FNE*PK(10)
      O=1.0+GO+FHE*PK(5)+FCE*PK(7)+FNE*PK(12)+FSE*PK(28)+FKE*PK(27)+
     *  FTE*PK(31)+FHE*FHE*PK(4)
      OO=2.0*FOE*(PK(11)+FTE*PK(32))
      S=1.0+GS+FHE*PK(18)+FCE*PK(22)+FNE*PK(25)+FOE*PK(28)+FKE*PK(30)
      SS=2.0*FSE*PK(29)
      XK=1.0+GK+FHE*PK(19)+FCE*PK(23)+FNE*PK(26)+FOE*PK(27)+FSE*
     *   PK(30)+FCE*FCE*PK(24)
      T=1.0+GT+FOE*PK(31)+FOE*FOE*PK(32)+FCE*FCE*PK(33)
      F(1)=FH*(H+HH)-1.0
      F(2)=FC*(C+CC+CCC)-AC
      F(3)=FO*(O+OO)-AO
      F(4)=FN*(XN+XNN)-AN
      F(5)=FH*(GH-PK(1)+FHE*GH*PK(3))+FC*GC+FN*GN+FO*GO+FS*GS+FK*GK
     *     +XNEN-FE
      F(6)=FS*(S+SS)-AS
      F(7)=FK*XK-AK
      F(8)=FT*T-AT
      A(1,1)=H+2.0*HH
      A(1,2)=FHE*(PK(6)+FNE*PK(15)+2.0*FCE*PK(16)+3.0*FCE*FCE*PK(
     *       20)+4.0*FCE*FHE*PK(14))
      A(1,3)=FHE*(PK(5)+2.0*FHE*PK(4))
      A(1,4)=FHE*(PK(13)+FCE*PK(15))
      A(1,5)=-FHE*(2.0*FHE*(PK(2)+GH*PK(3))+FCE*PK(6)+FNE*PK(13)+
     *       FOE*PK(5)+FSE*PK(18)+FKE*PK(19))-2.0*FHE*(FCE*(FNE*
     *       PK(15)+FCE*PK(16))+2.0*FHE*FOE*PK(4))-3.0*FCE*FCE*
     *       FHE*(FCE*PK(20)+2.0*FHE*PK(14))
      A(1,6)=FHE*PK(18)
      A(1,7)=FHE*PK(19)
      A(1,8)=0.
      A(2,1)=FCE*(PK(6)+FNE*PK(15)+FCE*(4.0*FHE*PK(14)+2.0*PK(16)+
     *       3.0*FCE*PK(20)))
      A(2,2)=C+2.0*CC+3.0*CCC
      A(2,3)=FCE*PK(7)
      A(2,4)=FCE*(PK(8)+FHE*PK(15))
      A(2,5)=-FCE*(2.0*FCE*PK(9)+FHE*PK(6)+FNE*PK(8)+FOE*PK(7)+
     *       FSE*PK(22)+FKE*PK(23))-2.0*FCE*(FHE*FNE*PK(15)+FCE*(3.0
     *       *FCE*PK(21)+2.0*FHE*PK(16)+2.0*FKE*PK(24)))-3.0*FCE
     *       *FCE*FHE*(2.0*FHE*PK(14)+3.0*FCE*PK(20))
      A(2,6)=FCE*PK(22)
      A(2,7)=FCE*(PK(23)+2.0*FCE*PK(24))
      A(2,8)=4.*FCE*PK(33)
      A(3,1)=FOE*(PK(5)+FHE*PK(4)*2.0)
      A(3,2)=FOE*PK(7)
      A(3,3)=O+2.0*OO
      A(3,4)=FOE*PK(12)
      A(3,5)=-FOE*(2.0*FOE*PK(11)+FHE*PK(5)+FCE*PK(7)+FNE*PK(12)+FSE
     *       *PK(28)+FKE*PK(27)+2.0*FHE*FHE*PK(4))
      A(3,6)=FOE*PK(28)
      A(3,7)=FOE*PK(27)
      A(3,8)=FOE*PK(31)+4.*FOE*PK(32)
      A(4,1)=FNE*(PK(13)+FCE*PK(15))
      A(4,2)=FNE*(PK(8)+FHE*PK(15))
      A(4,3)=FNE*PK(12)
      A(4,4)=XN+2.0*XNN
      A(4,5)=-FNE*(2.0*FNE*PK(10)+FHE*PK(13)+FCE*PK(8)+FOE*PK(12)+
     *       FSE*PK(25)+FKE*PK(26)+2.0*FHE*FCE*PK(15))
      A(4,6)=FNE*PK(25)
      A(4,7)=FNE*PK(26)
      A(4,8)=0.
      A(5,1)=GH*(1.0+2.0*FHE*PK(3))-PK(1)
      A(5,2)=GC
      A(5,3)=GO
      A(5,4)=GN
      A(5,5)=-GH*FHE*FHE*PK(3)-1.0
      A(5,6)=GS
      A(5,7)=GK
      A(5,8)=GT
      A(6,1)=FSE*PK(18)
      A(6,2)=FSE*PK(22)
      A(6,3)=FSE*PK(28)
      A(6,4)=FSE*PK(25)
      A(6,5)=-FSE*(2.0*FSE*PK(29)+FHE*PK(18)+FCE*PK(22)+FNE*PK(25)+
     *       FOE*PK(28)+FKE*PK(30))
      A(6,6)=S+2.0*SS
      A(6,7)=FSE*PK(30)
      A(6,8)=0.
      A(7,1)=FKE*PK(19)
      A(7,2)=FKE*(PK(23)+2.0*FCE*PK(24))
      A(7,3)=FKE*PK(27)
      A(7,4)=FKE*PK(26)
      A(7,5)=-FKE*(FHE*PK(19)+FCE*PK(23)+FNE*PK(26)+FOE*PK(27)+FSE*
     *       PK(30)+2.0*FCE*FCE*PK(24))
      A(7,6)=FKE*PK(30)
      A(7,7)=XK
      A(7,8)=0.
      A(8,1)=0.
      A(8,2)=FTE*FCE*PK(33)
      A(8,3)=FTE*(PK(31)+FOE*PK(32))
      A(8,4)=0.
      A(8,5)=-FTE*(FOE*PK(31)+2.*FOE*FOE*PK(32)+2.*FCE*FCE*PK(33))
      A(8,6)=0.
      A(8,7)=0.
      A(8,8)=T
CCC
      RETURN
      END
C
      SUBROUTINE MONTON(XI,N)
      implicit real*16 (a-h,o-z)
C
C      PARAMETER(NDIM=100)
      DIMENSION XI(N)
C
      DO 50 I=2,N
      IF(XI(I).LT.XI(I-1)) GOTO 50
      XI(I)=XI(I-1)-.05
50    CONTINUE
C
C
      RETURN
C
C
      E  N  D
C
      SUBROUTINE MULT(A,B,C,D,N,M)
      implicit real*16 (a-h,o-z)
C
C MULT MULTIPLIES MATRIX C (N*M) WITH SQUARE MATRIX B (N*N) AND PLACES
C THE RESULT IN MATRIX A. MATRIX D (N*M) IS USED FOR SCRATCH.
C FIRST DIMENSION MUST BE NDP.
C TIMING 850 MS FOR 40*40 MATRICES AND FORH.
C
      include 'parameter.inc'
      DIMENSION A(NDP,M),B(NDP,N),C(NDP,M),D(NDP,M)
C
C ZEROSET
      DO 100 J=1,M
      DO 100 I=1,N
100   D(I,J)=0.
C
C MULTIPLY
      DO 200 J=1,M
      DO 200 I=1,N
      DO 200 K=1,N
200   D(I,J)=B(I,K)*C(K,J)+D(I,J)
C
C RESTORE
      DO 300 J=1,M
      DO 300 I=1,N
300   A(I,J)=D(I,J)
C
      RETURN
      END
C
      SUBROUTINE OLDARC(IARCH)
      implicit real*16 (a-h,o-z)
C
C OLDARC RESTARTS FROM OLD ARCHIV FILE ON LUN 'IARCH'
C
      include 'parameter.inc'
C
      COMMON /STATEC/PPR(NDP),PPT(NDP),PP(NDP),GG(NDP),ZZ(NDP),DD(NDP),
     &   VV(NDP),FFC(NDP),PPE(NDP),TT(NDP),TAULN(NDP),RO(NDP),NTAU,ITER
      COMMON /ROSSC/ROSS(NDP),CROSS(NDP)
      COMMON /TAUC/TAU(NDP),DUMTAU(NDP),JTAU
      COMMON /CTEFF/TEFF,FLUX /CG/G,KONSG
      DIMENSION ABUND(20)
C
C READ OLD ARCHIV FILE

      READ(IARCH) DUM
      READ(IARCH) TEF,FLX,GD,PALFA,PNY,PY,PBETA,ILINE,ISTRAL,MIHAL,
     &            IDRAB1,IDRAB2,IDRAB3,IDRAB4,IDRAB5,IDRAB6,
     &            ITER,NEL,(ABUND(I),I=1,NEL)
      ITER=0
      READ(IARCH) DUM
      READ(IARCH) NTAU
      PI=3.14159
      DO 100 K=1,NTAU
      READ(IARCH) A,A,A,ZZ(K),TT(K),PPE(K),PG,PPR(K),PPT(K),ROSS(K)
      PP(K)=PG+PPR(K)+PPT(K)
      GG(K)=0.
      ZZ(K)=0.
      READ(IARCH) RRO,EMU,CP,CV,AGRAD,Q,U,V,ANCONV
      FFC(K)=FLUX*ANCONV
      VV(K)=V
      DD(K)=0.
      IF(VV(K).GT.0.) DD(K)=2.*PI*FFC(K)/(PALFA*RRO*CP*TT(K)*VV(K))
      READ(IARCH) DUM
      READ(IARCH) DUM
100   CONTINUE
C
C INTERPOLATE SURFACE TEMPERATURE
      YB=TAU(1)/TAU(2)
      YA=1.-YB
      TT(1)=YA*TT(1)+YB*TT(2)
C
C PRINT MODEL INFORMATION
      WRITE(7,51) IARCH
51    FORMAT('0MODEL FROM OLD ARCHIV FILE ON LUN',I3)
      GRAV=log10(G)
      WRITE(7,50) TEFF,GRAV,IDRAB1,IDRAB2,IDRAB3,IDRAB4,IDRAB5,IDRAB6
50    FORMAT(' TEFF=',F6.0,'  log10(G)=',F5.2,'  IDENTIFICATION ',6A4)
C
C MOVE INTEGERS TO HALF INTEGERS
      DO 110 K=2,NTAU
      FFC(NTAU+2-K)=0.5*(FFC(NTAU+2-K)+FFC(NTAU+1-K))
      VV(NTAU+2-K)=0.5*(VV(NTAU+2-K)+VV(NTAU+1-K))
      DD(NTAU+2-K)=0.5*(DD(NTAU+2-K)+DD(NTAU+1-K))
      ZZ(NTAU+2-K)=0.5*(ZZ(NTAU+2-K)+ZZ(NTAU+1-K))
110   CONTINUE
C
      REWIND IARCH
      RETURN
      END
C
      SUBROUTINE OLDSTA
      implicit real*16 (a-h,o-z)
      include 'parameter.inc'
C
C 'OLDSTA/NEWSTA' SAVES MODEL DATA ON AN ASCII FILE.
C
c      COMMON /STATEC/PPR(NDP),PPT(NDP),PP(NDP),GG(NDP),ZZ(NDP),DD(NDP),
c     &      VV(NDP),FFC(NDP),PE(NDP),T(NDP),TAULN(NDP),RO(NDP),NTAU,ITER
      COMMON /STATEC/B(11*NDP),RO(NDP),NTAU,ITER
      COMMON /CMOLRAT/ FOLD(NDP,8),MOLOLD,KL
C
C READ OLD STATE FROM unit 16

      READ(16,*) NTAU,ITER
      ITER=0
C      print*,'ntau,iter',ntau,iter
      DO 100 I=1,11
      IST=(I-1)*NDP
      READ(16,*) (B(IK),IK=IST+1,IST+NTAU)
C      WRITE(6,*) (B(IK),IK=IST+1,IST+NTAU)
100   CONTINUE
      WRITE(6,*) 'NTAU IN OLDSTA,MOLOLD = ',NTAU,MOLOLD
      MINTAU=10.*NDP+1
      MAXTAU=10.*NDP+NTAU
      DO 104 IK=MINTAU,MAXTAU
104   B(IK)=2.3025851*B(IK)    
C transform tauscale to ln for use in SOLVE
      IF (MOLOLD.EQ.1) READ(16,*) ((FOLD(K,I),I=1,8),K=1,NTAU)
C      write(7,*) ' oldsta; molold = ',molold
C      write(7,*) ' MOLOLD,KL,FOLD(1-NDP,1-8) = '
C      write(7,*) MOLOLD,KL,((FOLD(i,j),i=1,ndp),j=1,8)
      RETURN
      END
C
C WRITE NEW STATE
      SUBROUTINE NEWSTA
      implicit real*16 (a-h,o-z)
      include 'parameter.inc'
C
C 'OLDSTA/NEWSTA' SAVES MODEL DATA ON AN ASCII FILE.
C
      COMMON /STATEC/B(11*NDP),RO(NDP),NTAU,ITER
      COMMON /CMOLRAT/ FOLD(NDP,8),MOLOLD,KL
      COMMON /TAUC/TAU(NDP),DTAULN(NDP),JTAU
      COMMON /CIT/IT,ITMAX
C
C      REWIND 17

      OPEN(UNIT=17,FILE='arcivaab.dat',STATUS='UNKNOWN')
      WRITE(17,*) NTAU,IT
      MINTAU=10.*NDP+1
      MAXTAU=10.*NDP+NTAU
      DO 104 IK=MINTAU,MAXTAU
      K=IK-MINTAU+1
      B(IK)=log(TAU(K))
      B(IK)=0.4342945*B(IK)    
C transform tauscale to log10 for store
104   CONTINUE
      DO 100 I=1,11
      IST=(I-1)*NDP
      WRITE(17,*) (B(IK),IK=IST+1,IST+NTAU)
100   CONTINUE
      do 120 K=1,NTAU
120   WRITE(17,122) (FOLD(K,I),I=1,8)
122   FORMAT (1P7E11.3)
      CLOSE(17)
      RETURN
      END
C
C
      SUBROUTINE ONFROM(LUN,LR)
      implicit real*16 (a-h,o-z)
C
      include 'parameter.inc'
C
      COMMON /STATEC/A(10*NDP),TAULN(NDP),RO(NDP),NTAU,ITER
      COMMON /SPACE1/APP(13*NDP),TAULNP(NDP),SP1DUM((2*NDP-5)*NDP)
     &     ,space1dum2(3*ndp),space1dum3(2*ndp*ndp)
      DIMENSION AP(10*NDP)
      EQUIVALENCE (AP(1),APP(1))
      DATA LREC/0/
C
C RESUME PERMITS RESTART WITH NEW TAU-SCALE.
C
      REWIND LUN
      LREC=0
1     READ(LUN,END=2,ERR=99) APP,TAULNP,NTAUP
C
C OLD RECORD FORMAT
      GO TO 98
C
C NEW RECORD FORMAT
99    BACKSPACE LUN
      READ(LUN,END=2) AP,TAULNP,NTAUP,ITER
98    CONTINUE
      LREC=LREC+1
      GO TO 1
2     BACKSPACE LUN
C
C PRINT HEADING
      WRITE(7,52) LUN,LREC
52    FORMAT('0SAVED VALUES FROM LUN',I3,', RECORD',I2)
C
201   CONTINUE
      KP=2
      DO 100 K=1,NTAU
C
101   IF(TAULN(K).LE.TAULNP(KP)) GO TO 102
      IF(KP.EQ.NTAUP) GO TO 102
      KP=KP+1
      GO TO 101
C
102   CONTINUE
      P=(TAULN(K)-TAULNP(KP-1))/(TAULNP(KP)-TAULNP(KP-1))
      Q=1.-P
C
      DO 100 I=1,10
      J=(I-1)*NDP
      B=AP(KP+J)
      C=AP(KP+J-1)
      IF(B.LE.0..OR.C.LE.0.) GO TO 103
      A(K+J)=EXP(P*log(B)+Q*log(C))
      GO TO 100
103   A(K+J)=P*B+Q*C
100   CONTINUE
C
      RETURN
C
      ENTRY SAVEON(LUN,LR)
C
      LREC=LREC+1
      WRITE(LUN) A,TAULN,NTAU,ITER
C
C PRINT MESSAGE
      WRITE(7,51) LUN,LREC
51    FORMAT('0SAVED ON LUN',I3,', RECORD',I2)
C
      RETURN
C
C RESUME
      ENTRY RESUME(LUN,LR)
      DO 200 I=1,LR
200   READ(LUN) AP,TAULNP,NTAUP,ITER
      LREC=LR
      WRITE(7,52) LUN,LREC
      GO TO  201
C
      END
C
      SUBROUTINE OPAC(J,X,S)
      implicit real*16 (a-h,o-z)
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
      implicit real*16 (a-h,o-z)
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
      FUNCTION QAS(H,XL,A,Z,PFISH,IFISH)
      implicit real*16 (a-h,o-z)
C
C        THIS ROUTINE COMPUTES THE ASYMPTOTIC PARTS OF THE PARTITION
C        FUNCTIONS FOLLOWING
C           BASCHEK ET AL., ABH. HAMB. VIII, 26 (1966) IF IFISH = 0
C           FISCHEL AND SPARKS, AP. J. 164, 359 (1971) IF IFISH = 1
C           (APPROXIMATING THE ZETA FUNCTIONS BY INTEGRALS).
C
C        XL=QUANTUM NUMBER FOR THE FIRST LEVEL OF THE ASYMPTOTIC PART
C        H=QUANTUM NUMBER OF THE CUT (FOR IFISH=0)
C        A=DZ(FISCHEL AND SPARKS)=ALFA(BASCHEK ET AL.)
C        PFISH=P(FISCHEL AND SPARKS), ONLY NECESSARY IF IFISH = 1
C
C
      COMMON/UTPUT/IREAD,IWRIT
C
C        WHICH TYPE
      IF(IFISH.GT.0)GO TO 1
C
C        BASCHEK ET AL.
      QAS=0.333333*(H*(H+1.)*(H+0.5)-XL*(XL+1.)*(XL+0.5)) +
     *                       A*(H-XL)+0.5*A*A*(H-XL)/(H*XL)
      RETURN
C
C        FISCHEL AND SPARKS
    1 P=PFISH*Z
C
C        FISCHEL AND SPARKS, EQ. (26)
      P2=P*P
      P3=P2*P
      IF(P.LE.XL)GO TO 2
      XLM1=XL-1.
      R2=XLM1*XLM1
      R3=R2*XLM1
      QAS=1.3333333*P3+0.5*P2+0.16666667*P+1.33333333*A*P-0.4*A*A/P-
     *0.33333333*R3-0.5*R2-0.16666667*XLM1-A*XLM1+0.5*A*A/XL
      RETURN
C
C        FISCHEL AND SPARKS, EQ. (27)
    2 AXL2=A/(XL*XL)
      QAS=P3*P/XL*(1.+AXL2*(0.33333333+0.1*AXL2))
      RETURN
      END
C
      FUNCTION QTRAV(TETA,HP,J,JA)
      implicit real*16 (a-h,o-z)
C
C        HERE THE PARTITION FUNCTIONS ACCORDING TO TRAVING ET AL., ABH. HAMB.
C        STERNW. VIII, 1 (1966) ARE COMPUTED. THE SYMBOLS ARE GIVEN
C        IN THE COMMENTS AT THE BEGINNING OF SUBROUTINE INJON.
C        FUNCTION QAS IS CALLED.
C
C        DIMENSIONS NECESSARY
C        A(5),ASDE(KMAX),H(5),QPRIM(KMAX)
C        KMAX IS THE TOTAL NUMBER OF ELECTRON CONFIGURATIONS.
C        DIMENSIONS OF ARRAYS IN COMMON /CI3/ ARE COMMENTED ON IN SUBROUTINE
C        INJON.
C
      DIMENSION ASDE(80),H(5),QPRIM(80)
      COMMON/CI3/ALFA(300),GAMMA(300),G0(45),G2(80),XION(80),XL(80),
     *JBBEG(45),JCBEG(45),NK(45),NL(80),IFISH
      COMMON/CI7/A(5),PFISH,ITP
C
C
C        STATEMENT FUNCTION FOR 10.**
      EXP10(X)=EXP(2.302585*X)
C
      FLJ=J
      JB=JBBEG(JA)
      JC1=JCBEG(JA)
      NKP=NK(JA)
      QSUM=0.
C
C        WE START THE LOOP OVER DIFFERENT ELECTRON CONFIGURATIONS ('THE K-LOOP')
      DO5 K=1,NKP
      JC2=NL(JB)+JC1-1
C
C        IS TETA=PRECEDING TETA
      IF(ITP.GT.0)GO TO 4
      PRA=XION(JB)*TETA
      IF(PRA.LT.12.)GO TO 1
      ASDE(JB)=0.
      GO TO 2
    1 ASDE(JB)=G2(JB)*EXP10(-PRA)
C
    2 QPRIM(JB)=0.
      IF(NL(JB).LE.0)GO TO 4
      DO3 L=JC1,JC2
      PRE=GAMMA(L)*TETA
      IF(PRE.GT.12.)GO TO 3
      QPRIM(JB)=QPRIM(JB)+ALFA(L)*EXP10(-PRE)
    3 CONTINUE
    4 JC1=JC2+1
      QSUM=QPRIM(JB)+ASDE(JB)*QAS(HP,XL(JB),A(J),FLJ,PFISH,IFISH)
     *             +QSUM
    5 JB=JB+1
C        END OF 'THE K-LOOP'
      QTRAV=G0(JA)+QSUM
C
      RETURN
      END
C
      FUNCTION ROSSOP(T,PE)
      implicit real*16 (a-h,o-z)
C
C 'ROSSOP' CALCULATES THE ROSSELAND MEAN OPACIY AS DEFINED ON THE
C FIRST WAVELENGTH SCALE IN 'INITAB' INPUT. THIS VERSION WITH
C SCATTERING.                     (13/02/89)
C
      COMMON /CPF/PF,PFE,PFD,FIXROS,ITSTOP
      LOGICAL PF,PFE,PFD,FIXROS,ITSTOP
      DATA NEWT/2/
C
      CALL ABSKO(NEWT,1,T,PE,1,0,RSP,DUM)
      NEWT=1
      ROSSOP=RSP
      RETURN
C
      END
C
C
      SUBROUTINE ROSSOS
      implicit real*16 (a-h,o-z)
C
      include 'parameter.inc'
C
      DIMENSION X(NDP),S(NDP),SUMW(NDP),sumabs(ndp)
     *   ,sumwy(ndp),sumwyxs(ndp),dummy(ndp)
      CHARACTER MOLNAME*4,OSFIL*60,SAMPLING*3
C      REAL*8 Y,YA,SUMW,ROSSO,PTAUO
      COMMON/COS/WNOS(NWL),CONOS(NDP,NWL),WLOS(NWL),WLSTEP(NWL)
     *    ,KOS_STEP,NWTOT,NOSMOL,NEWOSATOM,NEWOSATOMLIST
     *    ,nchrom,OSFIL(30),MOLNAME(30),SAMPLING
      COMMON /STATEC/PPR(NDP),PPT(NDP),PP(NDP),GG(NDP),ZZ(NDP),DD(NDP),
     & VV(NDP),FFC(NDP),PPE(NDP),TT(NDP),TAULN(NDP),RO(NDP),NTAU,ITER
      COMMON /TAUC/TAU(NDP),DLNTAU(NDP),JTAU 
      COMMON /CG/GRAV,KONSG
      COMMON /CMOLRAT/ FOLD(NDP,8),MOLOLD,KL
      COMMON /CROSSOS/ ROSSO(NDP),PTAUO(NDP)
      COMMON /ROSSC/XKAPR(NDP),CROSS(NDP)
C      COMMON/CONLIN/rconop(nwl,ndp),rlineop(nwl,ndp)

C CALCULATE DETAILED ROSSELAND MEAN
C      KL=1
C      DUMMY(1)=ROSSOP(TT(1),PPE(1))

      write(7,*) ' ******** in rossos: '
      write(7,*) ' ntau,nwtot = ',ntau,nwtot
      DO 116 K=1,NTAU
C        KL=K
        SUMW(K)=0.
        ROSSO(K)=0.
        dummy(k) = 1.
C        sumwy(k) = 0.
C        sumwyxs(k) = 0.
C        sumabs(k) = 0.
C        DUMMY(k)=ROSSOP(TT(K),PPE(K))
116   CONTINUE
C          write(7,*) 
C     *     ' j,k,wlos(j),ya,xkapr(k),xkapr(k)*x(k),xkapr(k)*s(k)',
C     *     ' rosso(k),rconop(j,k),rlineop(j,k),tt(k) '
      DO 117 J=1,NWTOT
        CALL OPAC(J,X,S)
        Y=((WLOS(J)/1.E4)**2)**3
        DO 117 K=1,NTAU
          YA=EXP(-1.438E8/(TT(K)*WLOS(J)))
          YA=YA/(1.-YA)**2/Y
          SUMW(K)=SUMW(K)+WLSTEP(J)*YA
        if (wlos(j).le.5000. .or. wlos(j).ge.1.e5) go to 117
          ROSSO(K)=ROSSO(K)+WLSTEP(J)*YA/(xkapr(k)*(X(K)+S(K)))
C          sumwy(k) = sumwy(k) + WLSTEP(J)*YA
C          sumwyxs(k) = sumwyxs(k) + WLSTEP(J)*YA/(XKAPR(K)*(X(K)+S(K)))
C          sumabs(k)=sumabs(k)+XKAPR(K)*(x(k)+s(k))
C         if (j/50*50.eq.j  .and.  (k.eq.1.or.k.eq.7)) write(7,1171) 
C     *     j,k,wlos(j),ya,xkapr(k),xkapr(k)*x(k),xkapr(k)*s(k)
C     *     ,rosso(k),rconop(j,k),rlineop(j,k),tt(k)
117   CONTINUE
1171  format(i5,i3,1p8e12.3,0pf8.0)
C
C      write(7,*) ' k,sum(x(k)+s(k)),rosso(k),sumw(k),sumwy,sumwyxs = '
      DO 111 K=1,NTAU
C        write(7,*) k,sumabs(k),rosso(k),sumw(k),sumwy(k),sumwyxs(k)
        ROSSO(K)=SUMW(K)/ROSSO(K)
        PTAUO(K)=GRAV*TAU(K)/ROSSO(K)
111   CONTINUE
C
C      write(7,*)' k,pp,ptauo,rosso,tt ='
C      do 112 k=1,ntau
C      write(7,1172)k,pp(k),ptauo(k),rosso(k),tt(k)
C112   continue
1172  format(i5,1p3e12.3,0pf8.0)
C      write(7,*) ' ******** end from rossos '

      RETURN
      END
C
      SUBROUTINE SCALE(LUN)
      implicit real*16 (a-h,o-z)
C
C THIS ROUTINE SCALES THE MODEL PRESENTLY IN THE COMMON /STATEC/ (AND
C IN THE MODEL FILE) TO THE NEW EFFECTIVE TEMPERATURE AND GRAVITY GIVEN
C IN THE ORDINARY COMMON'S /CTEFF/ AND /CG/.
C
      include 'parameter.inc'
C
      COMMON /STATEC/PPR(NDP),PPT(NDP),PP(NDP),GG(NDP),ZZ(NDP),DD(NDP),
     *VV(NDP),FFC(NDP),PPE(NDP),TT(NDP),TAULN(NDP),ROSTAT(NDP),NTAU,ITER
      COMMON /ROSSC/ROSS(NDP),CROSS(NDP)
      COMMON /CTEFF/TEFF,FLUX /CG/GRAV,KONSG
      COMMON /NATURE/BOLTZK,CLIGHT,ECHARG,HPLNCK,PI,PI4C,RYDBRG,
     *STEFAN
      COMMON /MIXC/PALFA,PBETA,PNY,PY    /CSTYR/ MIHAL,NOCONV
      COMMON /TAUC/TAU(NDP),DTAULN(NDP),JTAU
      COMMON /CMOLRAT/ FOLD(NDP,8),MOLOLD,KL
      COMMON /CISPH/ISPH
C
C READ OLD MODEL

      CALL OLDSTA
      ITER=0
C
C TRY TO READ OLD MODEL FILE ARCHIV RECORDS
      READ(LUN,END=100,ERR=100) IDUM
      READ(LUN,END=100,ERR=100) TEFOLD,FDUM,GRVOLD
      GO TO 110
C
C UNSUCCESSFUL, NO OLD ARCHIV DATA
100   CONTINUE
      DO 101 K=1,NTAU
      IF(TAU(K).GE.0.95) GO TO 102
101   CONTINUE
102   TEFOLD=TT(K)/1.08
      GRVOLD=1.
C EMPIRICAL FITTING
C
C PRINT
110   G1=log10(GRAV)
      G2=log10(GRVOLD)
      WRITE(7,50) TEFF,G1,TEFOLD,G2
50    FORMAT('0SCALING TO TEFF=',F6.0,' log10(G)=',F5.2,' FROM TEFF='
     &,F6.0,' log10(G)=',F5.2)
C
C SCALE THE TEMPERATURE
      DO 111 K=1,NTAU
      PPR(K)=PPR(K)*(TEFF/TEFOLD)**4
      FFC(K)=FFC(K)*(TEFF/TEFOLD)**4
111   TT(K)=TT(K)*TEFF/TEFOLD
C
C TRY TO AVOID DIVIDE BY ZERO IN TRYCK
      DO 112 K=1,NTAU
        CROSS(K)=1.
  112 CONTINUE
C
C INTEGRATE PRESSURE EQUATION
        if (isph.eq.1) then 
          CALL TRYCK_sph
        else 
          CALL TRYCK
        end if
C
C SCALE DD AND VV
      DO 120 K=NOCONV,NTAU
      KL=K
      IF(FFC(K).EQ.0.) GO TO 120
      TMEAN=.5*(TT(K)+TT(K-1))
      PEMEAN=.5*(PPE(K)+PPE(K-1))
      PRMEAN=.5*(PPR(K)+PPR(K-1))
      ROSSMN=.5*(ROSS(K)+ROSS(K-1))
      CALL TERMON(TMEAN,PEMEAN,PRMEAN,PG,PGT,PGPE,RO,ROT,ROPE,CP,ADIA,Q)
      HSCALE=(PG+PRMEAN)/GRAV/RO
      OMEGA=PALFA*HSCALE*RO*ROSSMN
      THETA=OMEGA/(1.+PY*OMEGA**2)
      GAMMA=CP*RO/(8.*STEFAN*TMEAN**3*THETA)
      DD(K)=(GRAV*HSCALE*Q/PNY*(PALFA**2*RO*CP*TMEAN/(2.*PI*FFC(K)))**2
     &)**(-.333333)
      VV(K)=PALFA*SQRT(GRAV*HSCALE*Q/PNY*DD(K))
      GG(K)=GAMMA*VV(K)
120   CONTINUE
      RETURN
      END
C 
      SUBROUTINE SETDIS(NLB,XLB,LMAX,IFIRST,ILAST)
      implicit real*16 (a-h,o-z)
C
C        THIS SUBROUTINE DISTRIBUTES THE NLB WAVELENGTHS GIVEN IN XLB IN
C        WAVELENGTH SETS, WITH MAX. LMAX WAVENLENGTHS IN EACH. THE FIRST
C        SET NUMBER IS IFIRST, THE LAST ILAST. IF MORE SETS ARE NECESSAR
C        EXECUTION IS STOPPED WITH A PRINT-OUT.
C
      DIMENSION XLB(500)
      COMMON/CXLSET/XL(20,10),NSET,NL(10)
      COMMON/CLINE1/XLINLO,XLINUP,TSKAL(30),
     &             PESKAL(30),IPEBEG(30),IPEEND(30),LINUN,NTSKAL,NPSKAL
      COMMON/CLINE3/GLAMD(100),JLBDS
      COMMON/UTPUT/IREAD,IWRIT
        COMMON/CLINE4/ILINE
C
      NLP=1
      NSET=IFIRST
      DO2 J=1,NLB
       IF(ILINE.LE.0)GOTO1
      IF(XLB(J).LT.XLINLO.OR.XLB(J).GT.XLINUP)GO TO 1
      DO11 K=1,JLBDS
      KMEM=K
      IF(XLB(J).LE.GLAMD(K))GO TO 12
   11 CONTINUE
   12 IF(J.EQ.NLB)GO TO 1
      IF(XLB(J+1).LE.GLAMD(KMEM))GO TO 2
    1 CONTINUE
      XL(NLP,NSET)=XLB(J)
      IF(J.EQ.NLB)GO TO 3
      NLP=NLP+1
      IF(NLP.LE.LMAX)GO TO 2
      NL(NSET)=NLP-1
      NLP=1
      NSET=NSET+1
      IF(NSET.GT.ILAST)GO TO 4
    2 CONTINUE
C
    4 WRITE(IWRIT,200)
      STOP 'SETDIS 1'
C
    3 NL(NSET)=NLP
      RETURN
C
  200 FORMAT(1H ,'TOO FEW SETS ALLOWED OR TOO MANY WAVELENGTH POINTS',
     &' WANTED IN SETDIS')
      END
C
      SUBROUTINE STARTM
      implicit real*16 (a-h,o-z)
C
C 'STARTM' FINDS A STARTING MODEL, WITH CONVECTIVE FLUX, AND WITH
C APPROXIMATELY THE CORRECT EFFECTIVE TEMPERATURE.
C
C DIMENSIONS
      include 'parameter.inc'
C
      DIMENSION PTAU(NDP)
C
C STATE VARIABLES
      COMMON /STATEC/PPR(NDP),PPT(NDP),PP(NDP),GG(NDP),
     &               ZZ(NDP),DD(NDP),VV(NDP),FFC(NDP),
     &               PPE(NDP),TT(NDP),TAULN(NDP),ROSTAT(NDP),NTAU,ITER
C
C CONNECTIONS VIA COMMON.
C THE COMMENTED COMMONS MUST BE INITIATED OUTSIDE THIS ROUTINE BEFORE IT
C IS CALLED.
C JTAU=NUMBER OF TAUPOINTS, TAU=TAUSCALE.
C MIHAL=LOWER LIMIT OF RADIATIVE EQUILIBRIUM CONDITION, TAUMAX NOT USED.
C PALFA,PBETA,PNY,PNY = MIXING LENGTH THEORY COEFFICIENTS.
C GRAV=SURFACE GRAVITY, TEFF=EFFECTIVE TEMPERATURE, FLUX=STEFAN*TEFF**4/
      COMMON /TAUC/TAU(NDP),DTAULN(NDP),JTAU
      COMMON /CSTYR/MIHAL,NOCONV /DEBUG/KDEBUG
      COMMON /MIXC/PALFA,PBETA,PNY,PY /CVFIX/VFIX
      COMMON /CG/GRAV,KONSG /CTEFF/TEFF,FLUX
      COMMON /NATURE/BOLTZK,CLIGHT,ECHARG,HPLNCK,PI,PI4C,RYDBRG,
     *STEFAN
C OWN COMMONS
      COMMON /ROSSC/ROSS(NDP),CROSS(NDP)
      COMMON /CARC1/ISTRAL,IDRAB1,IDRAB2,IDRAB3,IDRAB4,IDRAB5,
     &              IDRAB6,IARCH
      COMMON /CI8/PGC,RHOC,EC
      COMMON /CMOLRAT/ FOLD(NDP,8),MOLOLD,KL
      DATA IEDIT/8/
C
C STATEMENT FUNCTIONS
      TF(TAUX)=TEFF*EFF*(.75*(TAUCNV+TAUX))**.25-
     &          EFF*DTBLNK*(1.-2.*TAUX/(TAUBLN+TAUX))
C
C TIME

      CALL CLOCK
C
C START UP
      READ(5,66) TAUCNV,DTBLNK,TAUBLN
      EFF=1.09*(0.75*(TAUCNV+1.))**(-0.25)
      TSURF=TF(TAU(1))
      IF(TSURF/TEFF.LT.0.60) EFF=EFF*0.60*TEFF/TSURF
      DT1=1.E10
      TEOLD=-DT1
94    CONTINUE
      WRITE(6,48) IEDIT,ITER,IDRAB1,IDRAB2,IDRAB3,IDRAB4,IDRAB5,IDRAB6
      WRITE(6,61) TAUCNV,DTBLNK,TAUBLN,EFF
      WRITE(6,52)
C
C THE PRESSURE BOUNDARY CONDITION IS PP(1)/PP(2) = PGFACT. TO FIND
C APPROXIMATE STARTING VALUES WE START WITH AN ARBITRARY ELECTRON PRESSU
C AND INTEGRATE THE PRESSURE EQUATION FROM 1 TO 2. THIS GIVES NEW GAS PR
C AND NEW ELECTRON PRESSURES. ITERATION IS PERFORMED UNTIL SUFFICIENT AC
C IS OBTAINED.
      TT(1)=TF(TAU(1))
      TT(2)=TF(TAU(2))
      PPR(1)=4./3.*STEFAN*TT(1)**4/CLIGHT
      PPR(2)=4./3.*STEFAN*TT(2)**4/CLIGHT
      PPE(1)=1.E-3
      PGFACT=(TAU(1)/TAU(2))**0.667
      PPE(2)=PPE(1)/SQRT(PGFACT)
      DPEL=0.
      DFDPEL=3.
      DO 90 K=1,20
      KL=1
      ROSS(1)=ROSSOP(TT(1),PPE(1))
      PG1=PGC
      PG2=PG1/PGFACT
      KL=2
      CALL PEMAKE(TT(2),PPE(2),PG2,PPE(2))
      ROSS(2)=ROSSOP(TT(2),PPE(2))
      DPG=GRAV*DTAULN(2)*.5*(TAU(1)/ROSS(1)+TAU(2)/ROSS(2))
      CALL ZEROF(log((PG2-PG1)/DPG),DPEL,DFDPEL)
      PPE(1)=PPE(1)*EXP(DPEL)
      IF(ABS(DPEL).LT.0.001) GOTO 88
90    CONTINUE
88    CONTINUE
      PP(1)=PG1+PPR(1)
      PP(2)=PG2+PPR(2)
      KL=1
      ROSS(1)=ROSSOP(TT(1),PPE(1))
      KL=2
      ROSS(2)=ROSSOP(TT(2),PPE(2))
      PTAU(1)=GRAV*TAU(1)/ROSS(1)
      PTAU(2)=GRAV*TAU(2)/ROSS(2)
C
C TAU LOOP
      RO=0.
      DO 99 K=1,NTAU
      KL=K
      PPT(K)=0.
      ZZ(K)=0.
      GG(K)=0.
      DD(K)=0.
      VV(K)=0.
      FFC(K)=0.
      IF(K.LE.2) GO TO 99
C
C NEW TEMPERATURE
      TT(K)=TF(TAU(K))
      PPR(K)=1.33*STEFAN*TT(K)**4/CLIGHT
      PPE(K)=PPE(K-1)*TAU(K)/TAU(K-1)
      PPT(K)=PPT(K-1)
C
C ITERATE THREE TIMES
      DO 92 IT=1,3
C
C COMPUTE ROSSELAND MEAN AND NEW PRESSURE
      KL=K
      ROSS(K)=ROSSOP(TT(K),PPE(K))
      PTAU(K)=GRAV*TAU(K)/ROSS(K)
      PP(K)=PP(K-1)+.5*DTAULN(K)*(PTAU(K)+PTAU(K-1))
      DLNP=log(PP(K)/PP(K-1))
C
C NEW ELECTRON PRESSURE
      PPT(K)=MIN(0.5*PP(K),PPT(K))
      PG=PP(K)-PPR(K)-PPT(K)
      PPEK=PPE(K)
      KL=K
      CALL PEMAKE(TT(K),PPEK,PG,PPE(K))
      IF(K.LE.NOCONV) GO TO 97
C
C CONSIDER CONVECTION
      TMEAN=.5*(TT(K)+TT(K-1))
      PEMEAN=.5*(PPE(K)+PPE(K-1))
      PRMEAN=.5*(PPR(K)+PPR(K-1))
      KL=K
      CALL TERMON(TMEAN,PEMEAN,PRMEAN,PG,PGT,PGPE,RO,ROT,ROPE,CP,ADIA,Q)
      IF(IT.EQ.1) GO TO 95
C
C USE OLD GRAD-ADIA AT SECOND IT-LOOP
      IF(FFC(K).EQ.0.) GO TO 97
      GO TO 96
95    CONTINUE
      GRAD=log(TT(K)/TT(K-1))/DLNP
      IF((GRAD.LE.ADIA.AND.VFIX.EQ.0.).OR.K.LE.NOCONV.OR.PALFA.EQ.0.)
     * GOTO 97
C
C IF CONVECTION ACCORDING TO THESE CONDITIONS AND ESTIMATED CONVECTIVE
C FLUX GREATER THAN TOTAL FLUX, ADJUST GRADIENT UNTIL ESTIMATED CONVECTI
C FLUX EQUALS TOTAL FLUX.
      HSCALE=(PG+PRMEAN)/GRAV/RO
      OMEGA=PALFA*HSCALE*RO*ROSS(K-1)
      THETA=OMEGA/(1.+PY*OMEGA**2)
      GAMMA=CP*RO/(8.*STEFAN*TT(K-1)**3*THETA)
      VV(K)=VVMLT(GRAD-ADIA,GRAV*HSCALE*Q*PALFA**2/PNY,GAMMA**2)
      GG(K)=GAMMA*VV(K)
      DD(K)=GG(K)/(1.+GG(K))*(GRAD-ADIA)
      FFC(K)=PALFA*RO*CP*TT(K-1)*VV(K)*DD(K)/2./PI
      FFCM=FLUX*TAU(K)/10.
      IF(FFCM.GT.FLUX) FFCM=FLUX
      IF(FFC(K).LT.FFCM.AND.VFIX.EQ.0.) GO TO 97
      DD(K)=(GRAV*HSCALE*Q/PNY*(PALFA**2*RO*CP*TT(K-1)/(2.*PI*FFCM))**2
     &)**(-.333333)
      FFC(K)=FFCM
      VV(K)=PALFA*SQRT(GRAV*HSCALE*Q/PNY*DD(K))
      GG(K)=GAMMA*VV(K)
C
C ADJUST TT(K) TO GIVE GRAD-ADIA
96    CONTINUE
      GRAD=ADIA+(1.+GG(K))/GG(K)*DD(K)
      TT(K)=TT(K-1)*EXP(GRAD*DLNP)
      EFF=EFF*TT(K)/TF(TAU(K))
      PPR(K)=1.33*STEFAN*TT(K)**4/CLIGHT
      PPT(K)=MIN(0.5*PP(K),PBETA*RO*VV(K)**2)
      PG=PP(K)-PPR(K)-PPT(K)
      PPEK=PPE(K)
      KL=K
      CALL PEMAKE(TT(K),PPEK,PG,PPE(K))
C
97    CONTINUE
      IF(ABS(TAU(K)-1.).GT.0.01) GOTO 93
C
C IMPROVE FIT TO EFFECTIVE TEMPERATURE
      TE=TT(K)-1.07*TEFF-DTBLNK
      IF(ABS(TE).LT.100.) GOTO 93
      DT1=DT1*TE/(TEOLD-TE)
      TT(1)=TT(1)+DT1
      EFF=EFF*TT(1)/TF(TAU(1))
      TEOLD=TE
      GO TO 94
93    CONTINUE
C
C END OF IT-LOOP
92    CONTINUE
C
C END OF TAU LOOP, PRINT.
99    WRITE(6,51) K,TAU(K),PPR(K),PPT(K),PP(K),GG(K),DD(K),VV(K),FFC(K),
     *PPE(K),TT(K),K
C
C TIME
      CALL CLOCK
C
      ITER=0
      RETURN
C
C FORMATS
45    FORMAT(' TIME',I6,' MSEC')
48    FORMAT('1',74X,'STARTM(',I1,')',5X,'ITERATION',I3,5X,6A4)
49    FORMAT(13('1234567890'),'123')
51    FORMAT(I8,1P10E12.4,I4)
52    FORMAT(T12,'TAU',T24,'PRAD',T36,'PTURB',T48,'PTOT',T60,'GAMMA',
     *T72,'DELTA',T84,'VCONV',T96,'FCONV',T108,'PE',T120,'TEMP')
61    FORMAT(' STARTING VALUES, TAUCNV,DTBLNK,TAUBLN,EFF=',4E10.3)
66    FORMAT(7(7X,F8.0))
      END
C
      SUBROUTINE TABS(NT,T)
      implicit real*16 (a-h,o-z)
C
C        THIS ROUTINE COMPUTES  FACTORS FOR INTERPOLATION IN T (TETA IF
C        ITETA(KOMP) IS GREATER THAN ZERO) IN THE ABKOF TABLE, INITIATED BY
C        SUBROUTINE INABS. CONCERNING THE OTHER CONTROL INTEGERS, SEE INABS.
C        THE RESULTING FACTORS ARE PUT IN AFAK. THE NUMBER OF FACTORS FOR
C        THE COMPONENT KOMP AT TEMPERATURE T(NTP) IS GIVEN IN
C        NOFAK((NKOMP-KOMPR)*(NTP-1)+KOMP-KOMPR). HERE KOMPR IS THE NUMBER
C        OF COMPONENTS WITH T-INDEP. COEFFICIENTS. NOFAK=0 MEANS THAT THE
C        ABSORPTION COEFFICIENT SHOULD BE =0. NPLATS (INDEX AS FOR NOFAK)
C        GIVES THE ARRAY INDEX OF THE TEMPERATURE POINT AT WHICH THE
C        INTERPOLATION IN ABKOF SHOULD START.
C
C        NT=NUMBER OF TEMPERATURES
C        T= ARRAY OF TEMPERATURES
C
C        DIMENSIONS NECESSARY
C        AFAK(KFADIM),NOFAK(IFADIM),NPLATS(IFADIM),T(1)
C        THE DIMENSIONS ARE LOWER LIMITS. DIMENSIONS OF ARRAYS IN COMMON /CA1/
C        AND /CA2/ ARE COMMENTED ON IN SUBROUTINE INABS.
C        IFADIM SHOULD BE AT LEAST =(NKOMP-KOMPR)*NT, WHERE NKOMP IS THE NUMBER
C               OF COMPONENTS, KOMPR THE NUMBER OF TEMPERATURE-INDEPENDENT
C               COMPONENTS AND NT THE NUMBER OF TEMPERATURE POINTS (IN THE PARA-
C               METER LIST).
C        KFADIM SHOULD BE AT LEAST =KOMPR*NT+(NKOMP-KOMPR)*NT*NUM, WHERE NUM IS
C               BETWEEN 2 AND 3 AND DEPENDENT ON THE TYPE OF TEMPERATURE
C               INTERPOLATION USED.
C
C
      include 'parameter.inc'
C
C      PARAMETER (IFADIM=1000,KFADIM=4000)
      DIMENSION T(NT)
      COMMON/UTPUT/IREAD,IWRIT
      COMMON/CA1/DELT(30,2),TBOT(30,2),IDEL(30),ISVIT(30),ITETA(30),
     *KVADT(30),MAXET(30),MINET(30),NTM(30,2),NEXTT,NUTZT
      COMMON/CA2/ABKOF(4000),KOMPLA(600),KOMPR,KOMPS,NKOMP
      COMMON/CA4/AFAK(KFADIM),NOFAK(IFADIM),NPLATS(IFADIM)
C
      IFAK=1
      KFAK=1
      NSVIT=1
C        THIS IS JUST A DUMMY STATEMENT TO GIVE NSVIT A FORMAL VALUE
C
      DO81 NTP=1,NT
      TP=T(NTP)
      KFAK=KFAK+KOMPR
      DO81 KOMP=KOMPS,NKOMP
      IF(ISVIT(KOMP).GT.0)GO TO (51,61,70),NSVIT
      IF(ITETA(KOMP).LE.0)GO TO 2
    1 TS=5040./T(NTP)
      GO TO 3
    2 TS=T(NTP)
C
C        SEARCHING
    3 IF((TS-TBOT(KOMP,1)).GE.0.)GO TO 10
      IF(MINET(KOMP).LE.0)GO TO 70
C
C        EXTRAPOLATION DOWNWARDS
      IF(NEXTT.GT.0)WRITE(IWRIT,200)TS,KOMP
      INTA=1
      AP=(TS-TBOT(KOMP,1))/DELT(KOMP,1)
      IP=0
      GO TO 60
C
C        SEARCHING CONTINUES
   10 INTAP=1
      IDP=IDEL(KOMP)
      DO11 I=1,IDP
      AP=(TS-TBOT(KOMP,I))/DELT(KOMP,I)
      IP=INT(AP)
      INTA=IP+INTAP
      INAP=NTM(KOMP,I)-1+INTAP
      IF(INTA.LE.INAP) GO TO 20
   11 INTAP=INAP+1
      IF(MAXET(KOMP).LE.0)GO TO 70
C
C        EXTRAPOLATION DOWNWARDS
      IF(NEXTT.GT.0)WRITE(IWRIT,200)TS,KOMP
      INTA=INAP
      IP=NTM(KOMP,IDP)-1
      GO TO 60
C
   20 IF(KVADT(KOMP).LE.0)GO TO 60
C
C        QUADRATIC INTERPOLATION
   21 IF(INTA.LT.INAP)GO TO 50
      INTA=INTA-1
      IP=IP-1
C
   50 DXX1=AP-DFLOAT(IP)
      DXX2=DXX1-1.
      DXX3=DXX1-2.
      A1=DXX2*DXX3*0.5
      A2=-DXX1*DXX3
      A3=DXX1*DXX2*0.5
   51 AFAK(KFAK)=A1
      AFAK(KFAK+1)=A2
      AFAK(KFAK+2)=A3
      NPLATS(IFAK)=INTA
      NOFAK(IFAK)=3
      IFAK=IFAK+1
      KFAK=KFAK+3
      NSVIT=1
      GO TO 80
C
C        LINEAR INTER/EXTRAPOLATION
   60 A2=AP-DFLOAT(IP)
      A1=1.-A2
   61 AFAK(KFAK)=A1
      AFAK(KFAK+1)=A2
      NPLATS(IFAK)=INTA
      NOFAK(IFAK)=2
      IFAK=IFAK+1
      KFAK=KFAK+2
      NSVIT=2
      GO TO 80
C
C        OUTSIDE TABLE. ABS.COEFF. SHOULD BE = 0
   70 IF(NUTZT.GT.0)WRITE(IWRIT,201)TS,KOMP
      NOFAK(IFAK)=0
      IFAK=IFAK+1
      NSVIT=3
C
   80 CONTINUE
      IF(KFAK.GT.KFADIM)GO TO 90
      IF(IFAK.GT.IFADIM+1)GO TO 91
   81 CONTINUE
C
      GO TO 92
   90 WRITE(IWRIT,202)KFAK,KFADIM,NT
      STOP 'TABS 1'
   91 WRITE(IWRIT,203)IFAK,IFADIM,NT
      STOP 'TABS 2'
   92 CONTINUE
C
  200 FORMAT(33H EXTRAPOLATION IN TABS, T (TETA)=,E12.5,5X,
     *12HCOMPONENT NO,I5)
  201 FORMAT(24H ZERO IN TABS, T (TETA)=,E12.5,5X,12HCOMPONENT NO,I5)
  202 FORMAT(6H KFAK=,I5,5X,11H GT KFADIM=,I5,5X,12HIN TABS, NT=,I5)
  203 FORMAT(6H IFAK=,I5,5X,11H GT IFADIM=,I5,5X,12HIN TABS, NT=,I5)
      RETURN
      END
C
      SUBROUTINE TAET(T,PE,PG,RO,E)
      implicit real*16 (a-h,o-z)
C
C TAET SIMULATION

      CALL JON(T,PE,1,PG,RO,E,0)

      RETURN
      END
C
      SUBROUTINE TAUSCA
      implicit real*16 (a-h,o-z)
C
C 'TAUSCA' INITIATES A TAU SCALE FROM INPUT LOGTAU AND LOGTAU-DIFFERENCE
C *NORD*
C
      include 'parameter.inc'
C
      DIMENSION TAULNX(NDP)
      COMMON /TAUC/TAU(NDP),DTAULN(NDP),JTAU
      COMMON /STATEC/DUM1(10*NDP),TAULN(NDP),RO(NDP),NTAU,ITER

      K=1
C
C READ LOGTAU AND LOGTAU-DIFFERENCE
      READ(5,50) T1,D1
C
1     CONTINUE
      READ(5,50) T2,D2
      TLIM=T2-.5*D1
      T=T1
2     CONTINUE
C
C NEW TAU POINT
      TAU(K)=10.**T
      K=K+1
      T=T+D1
      IF(T.LT.TLIM) GO TO 2
C
C NEW LOGTAU-DIFFERENCE
      T1=T2
      D1=D2
      IF(D1.GT.0.) GO TO 1
C
C END OF TAU SCALE
      TAU(K)=10.**T
      JTAU=K
CUGJ      NTAU=JTAU  commented 12.6.90 to allow NTAU(input-model) .ne. JTAU
C
      DO 3 K=1,JTAU
3     TAULNX(K)=log(TAU(K))
      DO 4 K=2,JTAU
4     DTAULN(K)=TAULNX(K)-TAULNX(K-1)
C
      RETURN
50    FORMAT(5(7X,F8.0))
      END
C
      SUBROUTINE TERMO(T,PE,PRAD,P,RO,CP,CV,TGRAD,Q,U2)
      implicit real*16 (a-h,o-z)
C
C        RUTINEN BERAEKNAR  OVANSTAAENDE STORHETER UTGAAENDE
C        FRAAN T, PE OCH PRAD (STRAALNINGSTRYCKET). INGEN AV STOR-
C        HETERNA GES ELLER ERHAALLES LOGARITMERAD. METODEN MED DIF-
C        FERENSFORMLER AER VARDAYAS
C        VI ANVAENDER RUTINEN TAET.
C        *** OBSERVERA. OCKSAA TAET MAASTE ARBETA MED OCH GE
C        I C K E  L O G A R I T M E R A D E  S T O R H E T E R . *******
C
C VERSION OF 73.02.05. DOUBLE PRECISION ADDED. *NORD*
C
C 13-OCT-1998 14:19:53.82:
C vi fick negativa Cv i flera modeller. Det beror
C(troligen) paa att naer du beraeknar de partiella derivatorna du
C behoever foer att beraekna Cv saa tar man ett alltfoer kort steg
C ( 0.001 av T resp P ). Vi oekade till 0.005 eller 0.01 och daa gick
C det mycket baettre! (Antagligen fluktuerar part. deriv. starkt eller
C saa aer det helt enkelt numeriskt brus vid laaga temperaturer).
C Kjell

C
      DIMENSION PGH(4),ROH(4),EH(4),HH(4),PH(4),EP(4),PRAH(4)
C
      DEREP=0.001
      DERET=0.001
      DELPE=DEREP*PE
      DELT=DERET*T
      PINV=1./(2.*DELPE)
      TINV=1./(2.*DELT)
C
      CALL JON(T,PE,1,PG,RO,EPP,0)
      P=PRAD+PG
C
      PEP=PE-DELPE
      CALL JON(T,PEP,1,PGH(1),ROH(1),EP(1),0)
      PRAH(1)=PRAD
      PEP=PE+DELPE
      CALL JON(T,PEP,1,PGH(2),ROH(2),EP(2),0)
      PRAH(2)=PRAD
      TP=T-DELT
      CALL JON(TP,PE,1,PGH(3),ROH(3),EP(3),0)
      PRAH(3)=PRAD*(1.-4.*DERET)
      TP=T+DELT
      CALL JON(TP,PE,1,PGH(4),ROH(4),EP(4),0)
      PRAH(4)=PRAD*(1.+4.*DERET)
C
      DO1 I=1,4
      EH(I)=3.*PRAH(I)/ROH(I)+EP(I)
      PH(I)=PRAH(I)+PGH(I)
    1 HH(I)=EH(I)+PH(I)/ROH(I)
C
      DET=(EH(4)-EH(3))*TINV
      DEP=(EH(2)-EH(1))*PINV
      DROT=(ROH(4)-ROH(3))*TINV
      DROP=(ROH(2)-ROH(1))*PINV
      DHT=(HH(4)-HH(3))*TINV
      DHP=(HH(2)-HH(1))*PINV
      DPT=(PH(4)-PH(3))*TINV
      DPP=(PH(2)-PH(1))*PINV
      DPGT=(PGH(4)-PGH(3))*TINV
      DPGP=(PGH(2)-PGH(1))*PINV
C
      CV=DET-DEP*DROT/DROP
      CP=DHT-DHP*DPT/DPP
      HJALP=DROT-DROP*DPT/DPP
      TGRAD=-P*HJALP/(CP*RO*RO)
      Q=-T/RO*(DROT-DROP*DPGT/DPGP)
      U2=CP*DPP/(CV*DROP)
      if (U2.lt.0.) then
          write(7,*) ' U2 becomes negative in TERMO '
          write(7,*) ' T,RO,PE,PRAD,P,CV,CP,Q,DPP,DROP,TGRAD,U2:'
          write(7,777) T,RO,PE,PRAD,P,CV,CP,Q,DPP,DROP,TGRAD,U2
777       format(1p6e12.3)
      end if
C
      RETURN
      END
C
      SUBROUTINE TERMON(T,PE,PRAD,PG,DPGT,DPGP,RO,DROT,DROP,CP,TGRAD,Q)
      implicit real*16 (a-h,o-z)
C
C        RUTINEN BERAEKNAR  OVANSTAAENDE STORHETER UTGAAENDE
C        FRAAN T, PE OCH PRAD (STRAALNINGSTRYCKET). INGEN AV STOR-
C        HETERNA GES ELLER ERHAALLES LOGARITMERAD. METODEN MED DIF-
C        FERENSFORMLER AER VARDAYAS
C        VI ANVAENDER RUTINEN TAET.
C        *** OBSERVERA. OCKSAA TAET MAASTE ARBETA MED OCH GE
C        I C K E  L O G A R I T M E R A D E  S T O R H E T E R . *******
C
C VERSION OF 73.02.05. DOUBLE PRECISION ADDED. *NORD*
C
C
      DIMENSION PGH(4),ROH(4),EH(4),HH(4),PH(4),EP(4),PRAH(4)
C
      DEREP=0.01
      DERET=0.001
      DELPE=DEREP*PE
      DELT=DERET*T
      PINV=1./(2.*DELPE)
      TINV=1./(2.*DELT)
C
      CALL JON(T,PE,1,PG,RO,EPP,0)
      P=PRAD+PG
C
      PEP=PE-DELPE
      CALL JON(T,PEP,1,PGH(1),ROH(1),EP(1),0)
      PRAH(1)=PRAD
      PEP=PE+DELPE
      CALL JON(T,PEP,1,PGH(2),ROH(2),EP(2),0)
      PRAH(2)=PRAD
      TP=T-DELT
      CALL JON(TP,PE,1,PGH(3),ROH(3),EP(3),0)
      PRAH(3)=PRAD*(1.-4.*DERET)
      TP=T+DELT
      CALL JON(TP,PE,1,PGH(4),ROH(4),EP(4),0)
      PRAH(4)=PRAD*(1.+4.*DERET)
C
      DO1 I=1,4
      EH(I)=3.*PRAH(I)/ROH(I)+EP(I)
      PH(I)=PRAH(I)+PGH(I)
    1 HH(I)=EH(I)+PH(I)/ROH(I)
C
      DET=(EH(4)-EH(3))*TINV
      DEP=(EH(2)-EH(1))*PINV
      DROT=(ROH(4)-ROH(3))*TINV
      DROP=(ROH(2)-ROH(1))*PINV
      DHT=(HH(4)-HH(3))*TINV
      DHP=(HH(2)-HH(1))*PINV
      DPT=(PH(4)-PH(3))*TINV
      DPP=(PH(2)-PH(1))*PINV
      DPGT=(PGH(4)-PGH(3))*TINV
      DPGP=(PGH(2)-PGH(1))*PINV
C
      CV=DET-DEP*DROT/DROP
      CP=DHT-DHP*DPT/DPP
      HJALP=DROT-DROP*DPT/DPP
      TGRAD=-P*HJALP/(CP*RO*RO)
      Q=-T/RO*(DROT-DROP*DPGT/DPGP)
      U2=CP*DPP/(CV*DROP)
C
      RETURN
      END
C
      SUBROUTINE TINT(N,X,Y,XINT,YINT)
      implicit real*16 (a-h,o-z)
C
C        DENNA ENDIMENSIONELLA INTERPOLATIONSRUTIN ARBETAR MED SUCCESIV
C        HALVERING. I PRINCIP BAADE INTER- OCH EXTRAPOLERAR DEN VILLIGT.
C        VARNING SKRIVS DOCK UT VID EXTRAPOLATION. INTERPOLATIONERNA SKE
C        MED TREPUNKTSFORMEL, SUBR.  I N P 3  .
C        ***** OBSERVERA. TABELLERNA SKALL VARA V A E X A N D E ******
C
      DIMENSION X(N),Y(N),ARG(3),FUNK(3)
      COMMON/UTPUT/IREAD,IWRIT
C
      NOEV=N
      NED=1
    1 NP=(NED+NOEV)/2
      IF(XINT-X(NP))2,2,3
    2 NOEV=NP
      GO TO 4
    3 NED=NP
    4 IF(NOEV-NED-2)5,5,1
    5 IF(NOEV-2)7,7,6
    6 J=NOEV-3
      GO TO 8
    7 J=NOEV-2
    8 DO9 K=1,3
      JP=K+J
      ARG(K)=X(JP)
    9 FUNK(K)=Y(JP)
      IF(ARG(1)-XINT)11,11,12
   11 IF(ARG(3)-XINT)12,13,13
   12 WRITE(IWRIT,200)XINT,ARG
  200 FORMAT(37H VARNING. EXTRAPOLATION I TINT. XINT=,E12.5,3X,4HARG=,
     *3E12.5)
   13 CALL INP3(ARG,FUNK,XINT,YINT)
      RETURN
      END
C*
C*NEW PDS MEMBER FOLLOWS
C*
      SUBROUTINE SOLVE(NEW)
      implicit real*16 (a-h,o-z)
C
C SOLVE PERFORMS ONE NEWTON-RAPSON ITERATION ON THE MODELATMOSPHERE PROB
C INCLUDING LOCAL CONVECTION.THE STATE OF THE ATMOSPHERE IS DESCRIBED BY
C NUMBER OF VARIABLES SUCH AS TEMPERATURE,ELECTRON PRESSURE,TOTAL PRESSU
C CONVECTIVE FLUX ETC..TO EACH VARIABLE CORRESPONDS A CERTAIN CONDITIONA
C EQUATION WICH DETERMINES THAT VARIABLE,ASSUMING THE OTHER VARIABLES BE
C KNOWN.
C
C NAMING CONVENTION.THE VARIABLES HAVE NAMES WITH A DOUBLE OCCURANCE OF
C FIRST LETTER,CORRECTIONS (TO BEE COMPUTED IN THIS ITERATION) HAVE SING
C OCCURANCE OF FIRST LETTER.RIGHTHANDSIDENAMES BEGIN WITH R.
C
C VARIABLES ARE CENTERED ON INTEGER AND HALFINTEGER TAU-POINTS AS INDICA
C BY 'I' OR 'H' ON THE FOLLOWING COMMENT CARDS.
C
C VARIABLE CORRECTION                                           HALF-INT
C PPR      PR         RADIATION PRESSURE                             I
C PPT      PT         TURBULENT PRESSURE                        H
C PP       P          TOTAL PRESSURE                                 I
C GG       G          CONVECTIVE EFFICIENCY,GAMMA               H
C ZZ       Z          GEOMETRIC HEIGTH                          H
C DD       D          GRADIENT DIFFERENCE,DELTA-DELTAPRIME      H
C VV       V          CONVECTIVE VELOCITY                       H
C FFC      FC         CONVECTIVE FLUX                           H
C PPE      PE         ELECTRON PRESSURE                              I
C TT       T          TEMPERATURE                                    I
C XJ                  MEAN INTENSITY                                 I
C
C NAMES OF THE PARTIAL DERIVATIVES ARE FORMED WITH A FIRST PART FROM THE
C EQUATION TO WICH IT BELONGS,AND A SECOND PART WICH IS THE NAME OF THE
C VARIABLE WITH RESPECT TO WICH THE DERIVATIVE IS TAKEN.(DOUBLE OCCURANC
C IF NECESSARY TO AVOID CONFUSION).
C
C PRPR,PRT
C PTPT,PTV,PTPE,PTT
C PPP,PPPE,PPTT
C GGG,GV,GPE,GT
C ZZZ,ZPE,ZT
C DDD,DP,DG,DPE,DT
C VVV,VZ,VD,VPE,VT
C FCFC,FCP,FCD,FCV,FCPE,FCT
C PEPE,PEPR,PEPT,PEP,PET
C TXJ,TFC,TTT
C
C STATE VARIABLES
      include 'parameter.inc'
      COMMON /STATEC/PPR(NDP),PPT(NDP),PP(NDP),GG(NDP),ZZ(NDP),DD(NDP),
     *VV(NDP),FFC(NDP),PPE(NDP),TT(NDP),TAULN(NDP),stro(ndp),NTAU,ITER
      common /ckdtpe/dpex,kdtpe
C
C DIMENSIONS
      DIMENSION PTAU(NDP),ROSSP(NDP),SUMW(NDP),ROSST(NDP),ROSSPE(NDP)
     *,XL(500),W(500)
     *,XJ1(NDP),XJ2(NDP),XJ3(NDP),XJT1(NDP),XJT2(NDP),XJT3(NDP)
     *,XJPE1(NDP),XJPE2(NDP),XJPE3(NDP)
     *,PR(NDP),PRT(NDP,NDP),PRJ(NDP)
     *,PT(NDP),PTV(NDP),PTPE(NDP),PTT(NDP)
     *,P(NDP),PPPE(2*NDP),PPTT(2*NDP)
     *,GV(NDP),GPE(NDP),GT(NDP)
     *,DP(2*NDP),DG(NDP),DPE(NDP,NDP),DT(NDP,NDP),DV(NDP)
      DIMENSION D(NDP),DTS(2*NDP),DPS(2*NDP),DPES(2*NDP)
     *,V(NDP),VD(NDP),VPE(NDP,NDP),VT(NDP,NDP)
     *,FC(NDP),FCD(NDP),FCV(NDP),FCPE(NDP,NDP),FCT(NDP,NDP)
     *,PE(NDP),PEPE(NDP,NDP),PET(NDP,NDP)
     *,T(NDP),TTT(NDP,NDP),TPE(NDP,NDP),TJ1(NDP),TJ2(NDP)
     *,PRPE(NDP,NDP)
     *,SCRATC(NDP,NDP),DBPL(NDP)
     *,XT(NDP),ST(NDP),DLNX(NDP),XLOG(NDP)
     *,XPE(NDP),SPE(NDP)
     *,RPR(NDP),RP(NDP),RD(NDP),RV(NDP),RFC(NDP),RPE(NDP),RT(NDP)
      LOGICAL NEWV
      real*16 a,b,c,aa,bb,cc,aaa,ccc
      character*24 idmodl
C
C CONNECTIONS VIA COMMON.
C THE COMMENTED COMMONS MUST BE INITIATED OUTSIDE THIS ROUTINE BEFORE IT
C IS CALLED.
C JTAU=NUMBER OF TAUPOINTS, TAU=TAUSCALE.
C NLAM=NUMBER OF LAMBDAPOINTS, XL=LAMBDAPOINTS, W=INTEGRATIONWEIGHTS.
C MIHAL=LOWER LIMIT OF RADIATIVE EQUILIBRIUM CONDITION, TAUMAX NOT USED.
C PALFA,PBETA,PNY,PNY = MIXING LENGTH THEORY COEFFICIENTS.
C GRAV=SURFACE GRAVITY, TEFF=EFFECTIVE TEMPERATURE, FLUX=STEFAN*TEFF**4/
      CHARACTER MOLNAME*4,OSFIL*60,SAMPLING*3
      COMMON/COS/WNOS(NWL),CONOS(NDP,NWL),WLOS(NWL),WLSTEP(NWL)
     *    ,KOS_STEP,NWTOT,NOSMOL,NEWOSATOM,NEWOSATOMLIST
     *    ,nchrom,OSFIL(30),MOLNAME(30),SAMPLING
      COMMON /CLEVETAT/GEFF(NDP),PPRG(NDP),AMLOSS
      COMMON /CMOLRAT/ FOLD(NDP,8),MOLOLD,KL
      common /CPRINT/NPRINT
      COMMON /TAUC/TAU(NDP),DTAULN(NDP),JTAU
      COMMON /CVAAGL/XL,W,NLAM
      COMMON /CSTYR/MIHAL,NOCONV /DEBUG/KDEBUG
      COMMON /MIXC/PALFA,PBETA,PNY,PY /CVFIX/VFIX
      COMMON /CG/GRAV,KONSG /CTEFF/TEFF,FLUX
      COMMON /NATURE/BOLTZK,CLIGHT,ECHARG,HPLNCK,PI,PI4C,RYDBRG,
     * STEFAN
      COMMON /CPF/PF,PFE,PFD,FIXROS,ITSTOP
      LOGICAL PF,PFE,PFD,FIXROS,ITSTOP
CUGJ FFR in excess     COMMON /CSPHER/DIFLOG,RADIUS,RR(NDP),NCORE,FFR(NDP)
      dimension ffr(ndp)
      COMMON /CSPHER/DIFLOG,RADIUS,RR(NDP),NCORE
C OWN COMMONS
      COMMON /CTRAN/X(NDP),S(NDP),BPLAN(NDP),XJ(NDP),HFLUX(NDP),XK(NDP)
     &  ,dumtran(4*ndp),idumtran(3)
      COMMON /CANGLE/XMU(6),XMU2(6),H(6),MMU_PP
      COMMON /CSURF/HSURF,Y1(NRAYS)
      COMMON /ROSSC/ROSS(NDP),CROSS(NDP) /RHOC/RHO(NDP)
      COMMON /CARC1/ISTRAL,IDRAB1,IDRAB2,IDRAB3,IDRAB4,IDRAB5,IDRAB6,
     &              IARCH
      COMMON /CARC2/T,FC,FLUXME(NWL),TAU5(NDP),INORD
      COMMON /Cspec/spec(nwl,3),ispec
      COMMON /CI8/PGC,ROC,EC
      COMMON /NEWMO/NEWMOD
      COMMON /MASSE/RELM
      COMMON /CORRECT/KORT,KPP,TCONV
      COMMON /CIT/IT,ITMAX
C
C SPACE ALLOCATION
      COMMON /SPACE1/XJ1,XJ2,XJ3,TJ1,TJ2,XJT1,XJT2,XJT3,PRJ,TTT,PRT
     &  ,XJPE1,XJPE2,XJPE3,TPE,PRPE
      COMMON /SPACE2_PP/FCT,FCPE,PET,PEPE
! Dust
      common /cdustdata/ dabstable(max_lay,nwl),dscatable(max_lay,nwl),
     *    eps(max_eps,max_lay),temp(max_lay),pgas(max_lay),n_lay,idust      
      common /cdustopac/ dust_abs(ndp,nwl), dust_sca(ndp,nwl)

C
C DATA
      DATA IVERS,IEDIT/21,1/
C
C IN THIS SECTION THE MEAN INTENSITY IS ELIMINATED IN THE TRANSPORT EQUA
C LEAVING THE EXPLICIT TEMPERATURE DEPENDANCE OF FLUX AND RADIATION PRES
C IN THE MATRICES TTT AND PRT.
C
      IF (mmu_pp.GT.nrays) STOP  
     & ' solve: increase nrays in parameter.inc'
      IF (mmu_pp.GT.6) STOP ' solve: increase dimension to mmu_pp'
      ITER=ITER+1
      INORD=IEDIT+10*IVERS
      FNORD=.1*INORD
C
C ZEROSET
      DO 110 I=1,NTAU
      PPRG(I) = 0.
      RT(I)=0.
      RPR(I)=0.
      ROSSP(I)=0.
      SUMW(I)=0.
      DO 110 J=1,NTAU
      TTT(I,J)=0.
      TPE(I,J)=0.
      PRT(I,J)=0.
      PRPE(I,J)=0.
110   CONTINUE
      kdtpe = 0
C
C CALCULATE DETAILED ROSSELAND MEAN
      REWIND 11
      KL=1
      DUMMY=ROSSOP(TT(1),PPE(1))
      PGA=PGC
      DO 116 K=1,NTAU
        KL=K
        DUMMY=ROSSOP(TT(K),PPE(K))
        ROSS(K)=1.0
        RHO(K)=ROC
        SUMW(K)=0.
        ROSSP(K)=0.
116   CONTINUE
C       write(7,*) 'j,wlos(j),wlstep(j),ya,y,x(k),s(k) (k=10) in SOLVE='
      DO 117 J=1,NWTOT
        CALL OPAC(J,X,S)
        WRITE(11) X,S
        Y=((WLOS(J)/1.E4)**2)**3
        DO 117 K=1,NTAU
          YA=EXP(-1.438E8/(TT(K)*WLOS(J)))
          YA=YA/(1.-YA)**2/Y
          SUMW(K)=SUMW(K)+WLSTEP(J)*YA
        if (wlos(j).le.5000. .or. wlos(j).ge.1.e5) go to 117
          ROSSP(K)=ROSSP(K)+WLSTEP(J)*YA/(ROSS(K)*(X(K)+S(K)))
C         if (j/500*500.eq.j  .and.  k.eq.10)
C     *         write(7,1171) j,wlos(j),wlstep(j),ya,y,x(k),s(k)
117   CONTINUE
1171  format(i5,1p8e12.3)
      REWIND 11
      
C
C TEMPERATURE AND ELECTRON PRESSURE PERTURBATIONS.
C KEEP THEM SMALL, TO STAY ON THE LINEAR PART.
      DTX=0.001
      DPEX=0.001
      DO 111 K=1,NTAU
        KL=K
        T(K)=TT(K)*DTX
        PE(K)=PPE(K)*DPEX
        TT(K)=TT(K)+T(K)
        ROSSP(K)=SUMW(K)/ROSSP(K)
        SUMW(K)=0.
        ROSST(K)=0.
111   CONTINUE
      kdtpe = 1
C
C FIRST WAVELENGTH LOOP, TO CALCULATE XT,ST AND SAVE.
      REWIND 12
      DO 112 J=1,NWTOT
        CALL OPAC(J,X,S)
        WRITE(12) X,S
        Y=((WLOS(J)/1.E4)**2)**3
        DO 112 K=1,NTAU
          YA=EXP(-1.438E8/(TT(K)*WLOS(J)))
          YA=YA/(1.-YA)**2/Y
          SUMW(K)=SUMW(K)+WLSTEP(J)*YA
        if (wlos(j).le.5000. .or. wlos(j).ge.1.e5) go to 112
          ROSST(K)=ROSST(K)+WLSTEP(J)*YA/(ROSS(K)*(X(K)+S(K)))
112   CONTINUE
      DO 113 K=1,NTAU
        ROSST(K)=SUMW(K)/ROSST(K)
        SUMW(K)=0.
        ROSSPE(K)=0.
        TT(K)=TT(K)-T(K)
        PPE(K)=PPE(K)+PE(K)
113   CONTINUE
        kdtpe = 2            !information to tstgem about computing dpg/dpe
      REWIND 12
      CALL CLOCK

C
C SECOND WAVELENGTH LOOP, TO CALCULATE XPE,SPE AND SAVE.
      REWIND 14
        KL=1
      DUMMY=ROSSOP(TT(1),PPE(1))
      PGPE=PGC
      DO 114 J=1,NWTOT
        CALL OPAC(J,X,S)
        WRITE(14) X,S
        Y=((WLOS(J)/1.E4)**2)**3
        DO 114 K=1,NTAU
          YA=EXP(-1.438E8/(TT(K)*WLOS(J)))
          YA=YA/(1.-YA)**2/Y
          SUMW(K)=SUMW(K)+WLSTEP(J)*YA
        if (wlos(j).le.5000. .or. wlos(j).ge.1.e5) go to 114
          ROSSPE(K)=ROSSPE(K)+WLSTEP(J)*YA/(ROSS(K)*(X(K)+S(K)))
114   CONTINUE
      kdtpe = 3
      REWIND 14
      CALL CLOCK
C
C FROM THIS POINT ON, ROSS() HOLDS THE TRUE ROSSELAND MEAN.  CROSS HOLDS
C THE RATIO OF THE TRUE TO APPROXIMATE MEANS, WHICH ARE NEEDED IN TRYCK.
C      write (6,*) 'rosseland values'
C      write(7,*)' k,tt(k),ppe(k),ross(k),cross(k) in solve ='
      DO 123 K=1,NTAU
        KL=K
        ROSSPE(K)=SUMW(K)/ROSSPE(K)
        PPE(K)=PPE(K)-PE(K)
        CROSS(K)=ROSSP(K)/ROSSOP(TT(K),PPE(K))
        ROSS(K)=ROSSP(K)
        PTAU(K)=GRAV*TAU(K)/ROSS(K)
C        write (6,'(1x,i3,4(1pe12.3))')
C     &     k,ross(k),cross(k),rosst(k),rosspe(k)
C      write(7,1171)k,tt(k),ppe(k),ross(k),cross(k)
123   CONTINUE
      CALL CLOCK
      

C
C RIGHT HAND SIDE IN PRESSURE EQUATION
        KL=1
      CALL TAET(TT(1),PPE(1),PG,RO,DUM)
      DLNP=1./(1.+(ROSSPE(1)-ROSS(1))/ROSS(1)*PGA/(PGPE-PGA))
      RP(1)=GRAV*TAU(1)/(ROSS(1)*DLNP)+PPR(1)-PP(1)
C SIMPSONS RULE
      DO 101 K=2,NTAU
      F0=PTAU(K-1)
      F1=FOUR(PTAU,TAULN,K,NTAU)
      F2=PTAU(K)
      RP(K)=(F0+4.*F1+F2)*DTAULN(K)/6.-(PP(K)-PP(K-1))
101   CONTINUE
C
C TIME
      CALL CLOCK
      MSA=0
C     CALL MSLEFT(MSA)
C
C WAVELENGTH LOOP
      IF(PF) WRITE(6,48) FNORD,ITER,idmodl
      IF(PF) WRITE(6,59)
      FTOT=0.
      DO 150 J=1,NWTOT
      DO 130 K=1,NTAU
      BPLAN(K)=BPL(TT(K),WLOS(J))
130   DBPL(K)=DIVBP(TT(K),WLOS(J))
C
C CALCULATE OPACITY DERIVATIVES AT CONSTANT GAS PRESSURE
      READ (11) X,S
      READ (12) XT,ST
      READ (14) XPE,SPE
      DO 131 K=1,NTAU
        X(K)=X(K)/ROSS(K)
        S(K)=S(K)/ROSS(K)
        XT(K)=XT(K)/ROSST(K)
        ST(K)=ST(K)/ROSST(K)
        XPE(K)=XPE(K)/ROSSPE(K)
        SPE(K)=SPE(K)/ROSSPE(K)
C        if (j.eq.100.and.k.eq.1) write (7,'(1x,7(1pe10.2))')
C     &   x(k),t(k),xt(k),log(xt(k)/x(k)),log(xt(k)/x(k))*tt(k)/t(k)
        XLOG(K)=log10(X(K))
        XT(K)=log(XT(K)/X(K))*X(K)/T(K)
        ST(K)=log(ST(K)/S(K))*S(K)/T(K)
        XPE(K)=log(XPE(K)/X(K))*X(K)/PE(K)
        SPE(K)=log(SPE(K)/S(K))*S(K)/PE(K)
        DLNX(K)=XT(K)*(TT(K)/2.3)/X(K)
C        if (j.eq.100) write (7,'(1x,7(1pe10.2))')
C     &    x(k),xt(k)/x(k)*tt(k),xpe(k)/x(k)*ppe(k)
131   CONTINUE
C
C TIME
      MS=MSA
C     CALL MSLEFT(MSA)
      MSOPAC=MS-MSA
C
C SOLVE TRANSPORTEQUATION WITH OLD STRATIFICATION.
      CALL TRANEQ
c ??
      DO 132 K=1,NTAU
      IF(XT(K)*(XJ(K)-BPLAN(K)).GT.X(K)*DBPL(K))
     & XT(K)=X(K)*DBPL(K)/(XJ(K)-BPLAN(K))
132   CONTINUE
      MS=MSA
C     CALL MSLEFT(MSA)
      MSTRAN=MS-MSA
C
C FLUX TO PRINT
      HSURF=MAX(HSURF,1.D-99)
      HFLUX1=4.*PI*HSURF
      HFLUX2=4.*PI*HFLUX(NTAU-1)
      FLUXME(J)=HFLUX1/PI
      GFLUX1=4.*HSURF*WLOS(J)**2/CLIGHT
      GFLUX2=4.*HFLUX(NTAU-1)*WLOS(J)**2/CLIGHT
      FFLUX1=-2.5*log10(MAX(1.D-99,GFLUX1))
      FFLUX2=-2.5*log10(MAX(1.D-99,GFLUX2))
      spec(j,1) = wlos(j)
      spec(j,2) = hsurf
      spec(j,3) = fluxme(j)
C
C INITIATE MATRICES, TAU LOOP.
      DO 140 K=1,NTAU
C EQUATION OF RADIATIVE TRANSPORT
      IF(K.GT.1) GO TO 142
C UPPER BOUNDARY
      FKB=XK(1)/XJ(1)
      FKC=XK(2)/XJ(2)
      DC=(X(1)+X(2)+S(1)+S(2))*(TAU(2)-TAU(1))
      SQ3=SQRT(3.)
      SOURCE=(X(1)*BPLAN(1)+S(1)*XJ(1))/(X(1)+S(1))
      YA=0.
      YB=0.
      YC=0.
      DO 141 I=1,MMU_PP
      TAU1=TAU(1)*(X(1)+S(1))/XMU(I)
      EXP1=EXP(-TAU1)
      YA=YA+H(I)*XMU(I)*(1.-EXP1)
      YC=YC+H(I)*XMU(I)*EXP1*TAU1*(XPE(1)+SPE(1))/(X(1)+S(1))
141   YB=YB+H(I)*XMU(I)*EXP1*TAU1*(XT(1)+ST(1))/(X(1)+S(1))
      XG=HFLUX(1)+YA*SOURCE
      Y=DC**2*X(1)/(8.*(X(1)+S(1)))
      XJ1(K)=0.
      FH=HFLUX(1)/XJ(1)
      FG=XG/XJ(1)
      XJ2(K)=(FKC-FKB)-(.5*DC*(FG-YA*S(1)/(X(1)+S(1)))+Y)
      XJ3(K)=FKC
      XSC=X(2)+S(2)+X(1)+S(1)
      XJT2(1)=-(HFLUX(1)*.5*DC+Y*(XJ(K)-BPLAN(K))*2.)*(XT(1)+ST(1))/XSC
     * -Y*(XJ(1)-BPLAN(1))*(S(1)*XT(1)-X(1)*ST(1))/X(1)/(X(1)+S(1))
     * +(YA*0.5*DC*X(1)/(X(1)+S(1))+Y)*DBPL(1)
     * +YB*0.5*DC*SOURCE
     * -YA*.5*DC*(XJ(1)-BPLAN(1))*(S(1)*XT(1)-X(1)*ST(1))/(X(1)+S(1))**2
      XJT3(1)=-(HFLUX(1)*.5*DC+Y*(XJ(1)-BPLAN(1))*2.)/XSC*(XT(2)+ST(2))
c ??
      XJPE2(1)=-(HFLUX(1)*.5*DC+Y*(XJ(K)-BPLAN(K))*2.)*(XPE(1)+SPE(1))/
     * XSC
     * -Y*(XJ(1)-BPLAN(1))*(S(1)*XPE(1)-X(1)*SPE(1))/X(1)/(X(1)+S(1))
     * +YC*0.5*DC*SOURCE
     * -YA*.5*DC*(XJ(1)-BPLAN(1))*(S(1)*XPE(1)-X(1)*SPE(1))/(X(1)+S(1))
     * **2
      XJPE3(1)=-(HFLUX(1)*.5*DC+Y*(XJ(1)-BPLAN(1))*2.)/XSC*(XPE(2)+SPE(2
     *  ))
      GO TO 144
142   DA=DC
      FKA=FKB
      FKB=FKC
      IF(K.EQ.NTAU) GO TO 143
C INNER PARTS
      FKC=XK(K+1)/XJ(K+1)
      DC=(X(K+1)+S(K+1)+X(K)+S(K))*(TAU(K+1)-TAU(K))
      DB=DA+DC
C NOTE THAT DA AND DC ARE 2* AND DB 4* THE NORMAL DTAU'S
      A=8./(DA*DB)
      C=8./(DC*DB)
      B=-A-C
      AA=A*XK(K-1)
      BB=B*XK(K)
      CC=C*XK(K+1)
      AAA=-((AA+BB*A/(A+C))*(1.+DA/DB)+(CC+BB*C/(A+C))*DA/DB)
      CCC=-((CC+BB*C/(A+C))*(1.+DC/DB)+(AA+BB*A/(A+C))*DC/DB)
      XSA=X(K-1)+S(K-1)+X(K)+S(K)
      XSC=X(K+1)+S(K+1)+X(K)+S(K)
      ABK=X(K)/(X(K)+S(K))
      XJT1(K)=AAA/XSA*(XT(K-1)+ST(K-1))
      XJT2(K)=(AAA/XSA+CCC/XSC)*(XT(K)+ST(K))+ABK*DBPL(K)
     & +(XT(K)*S(K)-ST(K)*X(K))/(X(K)+S(K))**2*(BPLAN(K)-XJ(K))
      XJT3(K)=CCC/XSC*(XT(K+1)+ST(K+1))
c ??
      XJPE1(K)=AAA/XSA*(XPE(K-1)+SPE(K-1))
      XJPE2(K)=(AAA/XSA+CCC/XSC)*(XPE(K)+SPE(K))
     & +(XPE(K)*S(K)-SPE(K)*X(K))/(X(K)+S(K))**2*(BPLAN(K)-XJ(K))
      XJPE3(K)=CCC/XSC*(XPE(K+1)+SPE(K+1))
      XJ1(K)=FKA*8./DA/DB
      XJ2(K)=(FKA-FKB)*8./(DA*DB)+(FKC-FKB)*8./(DC*DB)-ABK
      XJ3(K)=FKC*8./DB/DC
      GO TO 144
143   CONTINUE
C LOWER BOUNDARY
      XJ1(K)=0.
      XJ2(K)=-1.
      XJ3(K)=0.
      XJT1(K)=0.
      XJT2(K)=DBPL(K)
      XJPE1(K)=0.
      XJPE2(K)=0.
      XJPE3(K)=0.
144   CONTINUE
C
C TEMPERATURE EQUATION
      IF(K.GT.MIHAL) GO TO 145
C RADIATIVE EQUILIBRIUM
      Y=-WLSTEP(J)*X(K)
      IF(K.GT.2) Y=Y*DB/(X(K)+S(K))
      RT(K)=RT(K)-Y*(XJ(K)-BPLAN(K))
      TJ2(K)=Y
      TJ1(K)=0.
C...      TTT(K,K)=TTT(K,K)+MAX(0.0D+0,-Y*DBPL(K)+Y*(XJ(K)-BPLAN(K))
C...     & *XT(K)/X(K))
      TTT(K,K)=TTT(K,K)-Y*DBPL(K)+Y*(XJ(K)-BPLAN(K))*XT(K)/X(K)
      TPE(K,K)=TPE(K,K)+Y*(XJ(K)-BPLAN(K))*XPE(K)/X(K)
      GO TO 146
145   CONTINUE
C FLUX CONSTANCY
      Y=8.*WLSTEP(J)/DA
      RT(K)=RT(K)-Y*(XK(K)-XK(K-1))
      TJ1(K)=-FKA*Y
      TJ2(K)=FKB*Y
      TTT(K,K-1)=TTT(K,K-1)-Y*(XK(K)-XK(K-1))*(XT(K-1)+ST(K-1))/XSA
      TTT(K,K)=TTT(K,K)-Y*(XK(K)-XK(K-1))*(XT(K)+ST(K))/XSA
      TPE(K,K-1)=TPE(K,K-1)-Y*(XK(K)-XK(K-1))*(XPE(K-1)+SPE(K-1))/XSA
      TPE(K,K)=TPE(K,K)-Y*(XK(K)-XK(K-1))*(XPE(K)+SPE(K))/XSA
146   CONTINUE
C
C EQUATION OF RADIATION PRESSURE
      Y=PI4C*WLSTEP(J)
      RPR(K)=RPR(K)+Y*XK(K)
      PRJ(K)=-Y*FKB
      PPRG(K) = PPRG(K) + Y*HFLUX(K)*(X(K)+S(K))*ROSS(K)
C
C END OF TAU LOOP
140   CONTINUE
C
C ELIMINATE THIS WAVELENGTH
      IF(NEW.EQ.1) CALL ALGEBN(NTAU)
C
C TIME
      MS=MSA
C     CALL MSLEFT(MSA)
      MS=MS-MSA
C
C END OF WAVELENGTH LOOP
      HW1=HFLUX1*WLSTEP(J)
      FTOT=FTOT+HW1
      HW2=HFLUX2*WLSTEP(J)
      WAVEN=1.E4/WLOS(J)
      HSURF=MAX(1.D-99,HSURF)
      TRAD1=1.438E8/(WLOS(J)*
     #   log(1.0D+0+1.191D7*0.25D+0/HSURF*WAVEN**5))
      X01=log10(X(01))
      X25=log10(X(25))
      S01=log10(S(01))
      S25=log10(S(25))
C      IF(PF) WRITE(7,58) J,WLOS(J),WLSTEP(J),HFLUX1,HW1,WAVEN,GFLUX1
C     * ,FFLUX1
C     * ,TRAD1,X01,X25,S01,S25,HW2,DLNX(01),DLNX(25),MSOPAC,MSTRAN,MS
C      IF(PFD) WRITE(7,30) (XLOG(K),K=1,26)
C      IF(PFD) WRITE(7,30) (DLNX(K),K=1,26)
30    FORMAT(1X,26F5.2)
150   CONTINUE
      if (newmod .eq. 2) go to 900
      TEFFP=TEFF*(FTOT/FLUX/PI)**.25
      WRITE(7,65) FTOT,TEFFP
      DO 154 K=1,NTAU
      Y=0.
      DO 155 L=1,NTAU
155   Y=Y+TTT(K,L)
      IF(PFD) WRITE(7,151) K,Y,(TTT(K,L),L=1,NTAU)
      IF(K.LE.MIHAL.AND.Y.LE.0.)TTT(K,K)=TTT(K,K)-Y
  154 CONTINUE
151   FORMAT(' TTT='/(I3,E11.4,4(/10E11.4)))
156   CONTINUE
C
C TIME
      CALL CLOCK
C
C PRINT PRESSURE EQUATION
      IF(PF) WRITE(7,48) FNORD,ITER,idmodl
      IF(PF) WRITE(7,62)
      IF(PF) WRITE(7,63)
      DO 161 K=1,NTAU
        KL=K
      ROSST(K)=(ROSST(K)-ROSS(K))/T(K)
      ROSSPE(K)=(ROSSPE(K)-ROSS(K))/PE(K)
C
C TAU SCALES
      IF(K.GT.1) GO TO 162
      CALL KAP5(TT(1),PPE(1),ABSK)
      TAU5(1)=TAU(1)*ABSK/ROSS(1)
      TAUP=TAU(1)/CROSS(1)
      YC=ABSK/ROSS(K)
      YD=1./CROSS(K)
      GO TO 163
162   CONTINUE
      CALL KAP5(TT(K),PPE(K),ABSK)
      YA=ABSK/ROSS(K)
      YB=1./CROSS(K)
      TAU5(K)=TAU5(K-1)+.5*(YA+YC)*(TAU(K)-TAU(K-1))
      TAUP=TAUP+.5*(YB+YD)*(TAU(K)-TAU(K-1))
      YC=YA
      YD=YB
163   CONTINUE
      ROSPE=ROSSPE(K)*PPE(K)/ROSS(K)
      ROST=ROSST(K)*TT(K)/ROSS(K)
      IF(PF) WRITE(7,51) K,TAU(K),PTAU(K),ROSS(K),ROSPE,ROST
     * ,YC,YD,TAU5(K),TAUP,PT(K),K
161   CONTINUE
900   CONTINUE
      CALL CLOCK
      RETURN
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C+DEF,D=MATRIX
C
C IN THIS SECTION WE COMPUTE MATRIX ELEMENTS FOR THE REST OF THE PROBLEM
C
C
      ENTRY MATRIX
C TIME
C     CALL MSLEFT(MSA)
C
C ZEROSET
      IF(PF) WRITE(7,48) FNORD,ITER,idmodl
      IF(PF) WRITE(7,57)
      IF(PF) WRITE(7,55)
      GRAD=0.
      DO 201 I=1,NTAU
      DO 201 J=1,NTAU
      VT(I,J)=0.
      VPE(I,J)=0.
      FCT(I,J)=0.
      FCPE(I,J)=0.
      DT(I,J)=0.
      DPE(I,J)=0.
      PET(I,J)=0.
      PEPE(I,J)=0.
201   CONTINUE
C
C VFIX OPTION
      IF(VFIX.EQ.0.) GOTO 230
      DO 231 KK=1,NTAU
      K=(1+NTAU)-KK
      IF(K.LE.NOCONV) GOTO 230
231   VV(K)=VFIX*1.E5
230   CONTINUE
C
C TAU LOOP
      DO 200 K=1,NTAU
        KL=K
C
C TERMODYNAMICAL QUANTITIES WITH PARTIAL DERIVATIVES
      K1=MAX0(1,K-1)
      TMEAN=.5*(TT(K)+TT(K1))
      PEMEAN=.5*(PPE(K)+PPE(K1))
      PMEAN=.5*(PP(K)+PP(K1))
      PRMEAN=.5*(PPR(K)+PPR(K1))
      ROSSMN=.5*(ROSS(K)+ROSS(K1))
      IF(K.GT.NOCONV) GO TO 213
C NO CONVECTION
      CALL TERMON(TT(K),PPE(K),PPR(K),PG,PGT,PGPE,RO,ROT,ROPE,CP,ADIA,Q)
      RO1=RO
      PG1=PG
      PG1T=0.
      PG1PE=0.
      CPT=0.
      CPPE=0.
      ADIAT=0.
      ADIAPE=0.
      QT=0.
      QPE=0.
      GO TO 212
C CONVECTION
213   CONTINUE
      DERET=.01
      DEREP=.15
      TDELT=TMEAN*DERET
      PEDELT=PEMEAN*DEREP
      T1=TMEAN+TDELT
      PE1=PEMEAN+PEDELT
      CALL TERMON(T1,PEMEAN,PRMEAN,Y,YA,YB,YC,YD,YE,CP1,ADIA1,Q1)
      CALL TERMON(TMEAN,PE1,PRMEAN,Y,YA,YB,YC,YD,YE,CP2,ADIA2,Q2)
      CALL TERMON(TMEAN,PEMEAN,PRMEAN,PG1,PG1T,PG1PE,RO,ROT,ROPE,CP,
     &ADIA,Q)
      CALL TERMON(TT(K),PPE(K),PPR(K),PG,PGT,PGPE,RO1,Y,YA,YB,YC,YD)
      CPPE=(CP2-CP)/PEDELT
      ADIAPE=(ADIA2-ADIA)/PEDELT
      QPE=(Q2-Q)/PEDELT
      CPT=(CP1-CP)/TDELT
      ADIAT=(ADIA1-ADIA)/TDELT
      QT=(Q1-Q)/TDELT
212   CONTINUE
      RHO(K)=RO1
C
C DEPTH SCALE
      IF(K.EQ.1) ZZ(K)=0.
      IF(K.GT.1) ZZ(K)=ZZ(K-1)-.5*DTAULN(K)*
     & (TAU(K)/ROSS(K)/RHO(K)+TAU(K-1)/ROSS(K-1)/RHO(K-1))
      IF(TAULN(K).LT.0.0) KK0=K
C
C RADIATION PRESSURE
      RPR(K)=RPR(K)-PPR(K)
CC    PRJ IS NOW FREE
CC    PRT IS ALREADY INITIATED
CC    PRPR IS UNITY
C
C TOTAL PRESSURE
      IF(K.GT.1) GO TO 202
      Y=GRAV*TAU(1)/(ROSS(1)*DLNP)
      PPTT(1)=Y*ROSST(1)/ROSS(1)
      PPPE(1)=Y*ROSSPE(1)/ROSS(1)
      GO TO 203
202   Y=GRAV*DTAULN(K)*.5
      YY=Y*TAU(K)/ROSS(K)**2
      Y=Y*TAU(K-1)/ROSS(K-1)**2
      PPPE(K)=YY*ROSSPE(K)
      PPTT(K)=YY*ROSST(K)
      PPPE(K+NTAU-1)=Y*ROSSPE(K-1)
      PPTT(K+NTAU-1)=Y*ROSST(K-1)
203   CONTINUE
C
C CONVECTION EFFICIENCY GAMMA
      HSCALE=(PG1+PRMEAN)/GRAV/RO
      OMEGA=PALFA*HSCALE*RO*ROSSMN
      IF(PALFA.EQ.0.) OMEGA=HSCALE*RO*ROSSMN
      Y=PY*OMEGA**2
      YY=(Y-1.)/(Y+1.)
      THETA=OMEGA/(1.+Y)
      GAMMA=-CP*RO/(8.*STEFAN*TMEAN**3*THETA)
CC    GGG IS UNITY
      GV(K)=GAMMA
      IF(PBETA.GT.0.) VV(K)=MIN(VV(K),SQRT(0.5*PP(K)/PBETA/RO))
      GG(K)=-GAMMA*VV(K)
      ROSPEM=.5*(ROSSPE(K)+ROSSPE(K1))
      ROSSTM=.5*(ROSST(K)+ROSST(K1))
      GPE(K)=-GG(K)*(CPPE/CP+ROPE/RO+YY*(ROSPEM/ROSSMN+PG1PE/PG1))
      GT(K)=-GG(K)*(CPT/CP+ROT/RO-3./TMEAN+YY*(ROSSTM/ROSSMN+PG1T/PG1))
CC    RG IS ZERO
C
C GRADIENT DIFFERENCE
      IF(K.LE.NOCONV) GO TO 206
      DELP=PP(K)-PP(K-1)
      DELT=TT(K)-TT(K-1)
      PM=PP(K)+PP(K-1)
      TM=TT(K)+TT(K-1)
      Y=1.+GG(K)
      YY=-GG(K)/Y
      GRAD=log(TT(K)/TT(K-1))/log(PP(K)/PP(K-1))
      NEW V=DD(K).GT.0..AND.VV(K).EQ.0..AND.PALFA.GT.0..AND.K.GT.2
      IF(.NOT.NEWV) GO TO 263
      VV(K)=SQRT(GRAV*HSCALE*Q*PALFA**2*DD(K)/PNY)
      IF(PF) WRITE(7,64) K,VV(K)
      GO TO 203
263   CONTINUE
      NEW V=GRAD.GE.ADIA.AND.VV(K).EQ.0..AND.PALFA.GT.0..AND.K.GT.2
      IF(.NOT.NEWV) GO TO 204
      VV(K)=VVMLT(GRAD-ADIA,GRAV*HSCALE*Q*PALFA**2/PNY,GAMMA**2)
      IF(PF) WRITE(7,64) K,VV(K)
      GO TO 203
204   CONTINUE
      YYY=GRAD-ADIA
C DDD IS UNITY
C ******* NEXT STATEMENT FIXES T80G4M0 BUT NOT TESTED FOR ALL MODELS
      IF(DD(K).LE.0..AND.VV(K).GT.0.) DD(K)=-YY*YYY
      RD(K)=-YY*YYY-DD(K)
      DG(K)=-YYY/Y**2
      DT(K,K)=YY*(GRAD*(1./DELT-1./TM)-.5*ADIAT)
      DT(K,K-1)=YY*(GRAD*(-1./DELT-1./TM)-.5*ADIAT)
      DP(K)=YY*GRAD*(1./PM-1./DELP)
      DP(K+NTAU-1)=YY*GRAD*(1./PM+1./DELP)
      DPE(K,K)=-.5*YY*ADIAPE
      DPE(K,K-1)=-.5*YY*ADIAPE
      GO TO 205
206   DD(K)=0.
      RD(K)=0.
      DG(K)=0.
      DP(K)=0.
      DP(K+NTAU-1)=0.
205   CONTINUE
C
C VFIX OPTION
      IF(VFIX.EQ.0.) GOTO 280
      Y=0.
      IF(K.GT.NOCONV) Y=-VV(K)
      GOTO 207
280   CONTINUE
C
C CONVECTIVE VELOCITY
CC    VVV IS UNITY
      Y=0.
      IF(DD(K).LE.0.) GOTO 207
      Y=-SQRT(GRAV*HSCALE*Q*PALFA**2*DD(K)/PNY)
      VD(K)=Y*.5/DD(K)
      IF(-Y.GT.2.*VV(K)) VD(K)=VD(K)*2.
      VT(K,K)=.25*Y*(QT/Q+PGT/PG-ROT/RO)
      VT(K,K-1)=VT(K,K)
      VPE(K,K)=.25*Y*(QPE/Q+PGPE/PG-ROPE/RO)
      VPE(K,K-1)=VPE(K,K)
      GO TO 208
207   VD(K)=0.
208   CONTINUE
      RV(K)=-Y-VV(K)
C
C TURBULENT PRESSURE.
CC    PTPT IS UNITY
      Y=-PBETA*VV(K)**2
      PPT(K)=-RO*Y
      PTT(K)=Y*ROT
      PTPE(K)=Y*ROPE
      PTV(K)=-PBETA*2.*VV(K)*RO
C
C CONVECTIVE FLUX
      Y=-CP*RO*PALFA*TMEAN/2./PI
CC    FCFC IS UNITY
      YY=Y*VV(K)*DD(K)
      RFC(K)=-YY-FFC(K)
      FCD(K)=Y*VV(K)
      FCV(K)=Y*DD(K)
      IF(K.LE.NOCONV) GO TO 217
      FCT(K,K)=.5*YY*(CPT/CP+ROT/RO+1./TMEAN)
      FCT(K,K-1)=FCT(K,K)
      FCPE(K,K)=.5*YY*(CPPE/CP+ROPE/RO)
      FCPE(K,K-1)=FCPE(K,K)
217   CONTINUE
C
C ELECTRON PRESSURE
209   RPE(K)=PP(K)-PG-PPR(K)
      PET(K,K)=PGT
      PEPE(K,K)=PGPE
CC    PEPR=PEPT=1.  PEP=-1.
210   CONTINUE
C
C TEMPERATURE
CC    TTT IS ALREADY INITIATED
      IF(K.GT.MIHAL) GO TO 261
C STR"MGREN CONDITION
      RT(K)=RT(K)+(FFC(K+1)-FFC(K))
      GO TO 262
261   CONTINUE
C FLUXCONSTANCY
      RT(K)=RT(K)+FLUX-FFC(K)
262   CONTINUE
C
C END OF TAU LOOP
      PGT=PGT*TT(K)/PG
      PGPE=PGPE*PPE(K)/PG
      IF(PF) WRITE(7,51) K,TAU(K),HSCALE,ADIA,GRAD,CP,Q,PG,RO,PGPE
     * ,PGT,K
200   CONTINUE
C
C SUBTRACT CENTERED TURBULENT PRESSURE
      DO 216 K=2,NTAU
      K1=MIN0(K+1,NTAU)
216   RPE(K)=RPE(K)-.5*(PPT(K)+PPT(K1))
C
C SUBTRACT ZZ(TAU=1) FROM ZZ
      KK0 = MAX(1,KK0)
      ZZ0=ZZ(KK0)
      DO 283 K=1,NTAU
283   ZZ(K)=ZZ(K)-ZZ0
C
C TIME
      CALL CLOCK
C
C PRINT
      IF(PFE.OR.PFD) WRITE(7,48) FNORD,ITER,idmodl
      IF(PFE.OR.PFD) WRITE(7,54)
      IF(PFE.OR.PFD) WRITE(7,52)
      IF(PFE.OR.PFD)WRITE(7,51)(I,TAU(I),RPR(I),PPT(I),RP(I),GG(I),RD(I)
     * ,RV(I),RFC(I),RPE(I),RT(I),I,I=1,NTAU)
C
C SAVE DD-MATRICES
      DO 270 K=1,NTAU
      DPS(K)=DP(K)
      DPES(K)=DPE(K,K)
      DTS(K)=DT(K,K)
      D(K)=RD(K)
      IF(K.EQ.1) GO TO 270
      DPS(K+NTAU-1)=DP(K+NTAU-1)
      DPES(K+NTAU-1)=DPE(K,K-1)
      DTS(K+NTAU-1)=DT(K,K-1)
270   CONTINUE
C
C+DEF,D=ELIMIN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C GAUSS ELIMINATION TO UPPER TRIANGULAR FORM.
C
C TIME
      CALL CLOCK
C
C RADIATION PRESSURE
      RP(1)=RP(1)+RPR(1)
      DO 320 I=2,NTAU
      RPE(I)=RPE(I)-RPR(I)
      DO 301 J=1,NTAU
      PET(I,J)=PET(I,J)-PRT(I,J)
301   CONTINUE
C
C TURBULENT PRESSURE
      PET(I,I)=PET(I,I)-PTT(I)
      PEPE(I,I)=PEPE(I,I)-PTPE(I)
310   CONTINUE
C
C CONVECTIVE EFFICIENCY
      DT(I,I)=DT(I,I)-.5*DG(I)*GT(I)
      DT(I,I-1)=DT(I,I-1)-.5*DG(I)*GT(I)
      DPE(I,I)=DPE(I,I)-.5*DG(I)*GPE(I)
      DPE(I,I-1)=DPE(I,I-1)-.5*DG(I)*GPE(I)
      DV(I)=-DG(I)*GV(I)
C
C TOTAL PRESSURE
      DPI=-DP(I+NTAU-1)
      DPE(I,I-1)=DPE(I,I-1)-DPI*PPPE(I+NTAU-1)
      DT(I,I-1)=DT(I,I-1)-DPI*PPTT(I+NTAU-1)
      DPE(I,I)=DPE(I,I)-DPI*PPPE(I)
      DT(I,I)=DT(I,I)-DPI*PPTT(I)
      RD(I)=RD(I)-DPI*RP(I)
      DP(I)=DP(I)-DPI
320   CONTINUE
C
      RPI=0.
      DV(1)=0.
      DO 321 I=1,NTAU
C WE HAVE THE MATRICES       PPP       PPPE  PPTT
C AND                        PEPP      PEPE  PET
C WE WANT TO SUBTRACT FROM PEPE AND PET THE PRODUCTS OF PEPP*PPP-INVERS
C PPPE AND PET RESPECTIVELY. PPP IS BIDIAGONAL WITH UNITY ON THE DIAGONA
C MINUS UNITY ON THE SUBDIAGONAL. ITS INVERS IS A MATRIX WITH UNITY EVER
C UNDER AND ON THE DIAGONAL. PEPP IS MINUS UNITY. THUS THE PRODUCT OF PE
C INVERS WITH PPPE IS MATRIX OF THE FOLLOWING TYPE. IN EVERY COLUMN EACH
C ELEMENT IS THE SUM OF ALL ELEMENTS ABOVE THAT POINT (AND INCLUDING) IN
C PPPE MATRIX. SIMILARILY FOR THE PPTT-MATRIX.
      RPI=RPI+RP(I)
      RPE(I)=RPE(I)+RPI
      RD(I)=RD(I)-DP(I)*RPI
      Y=PPPE(I)
      YY=PPTT(I)
      YY=YY+PRT(1,I)
      JMIN=MAX0(I,1)
      DO 321 J=JMIN,NTAU
      IF(J.NE.I+1) GO TO 322
      Y=Y+PPPE(I+NTAU)
      YY=YY+PPTT(I+NTAU)
322   CONTINUE
      PEPE(J,I)=PEPE(J,I)+Y
      PET(J,I)=PET(J,I)+YY
      DPE(J,I)=DPE(J,I)-DP(J)*Y
      DT(J,I)=DT(J,I)-DP(J)*YY
321   CONTINUE
C
C GRADIENT DIFFERENCE
      DO 350 I=2,NTAU
      DVI=DV(I)
      VDI=VD(I)
      FCDI=FCD(I)
      RFC(I)=RFC(I)-FCDI*RD(I)
      RV(I)=RV(I)-VDI*RD(I)
      VVVI=MAX(0.5d+0,1.0D+0-VDI*DVI)
      FCV(I)=FCV(I)-FCDI*DVI
      DO 340 J=1,NTAU
      VT(I,J)=VT(I,J)-VDI*DT(I,J)
      VPE(I,J)=VPE(I,J)-VDI*DPE(I,J)
      FCT(I,J)=FCT(I,J)-FCDI*DT(I,J)
      FCPE(I,J)=FCPE(I,J)-FCDI*DPE(I,J)
340   CONTINUE
C
C TURBULENT VELOCITY
      PEVI=-PTV(I)
      RV(I)=RV(I)/VVVI
      FCVI=FCV(I)
      RFC(I)=RFC(I)-FCVI*RV(I)
      RPE(I)=RPE(I)-PEVI*RV(I)
      DO 350 J=1,NTAU
      VT(I,J)=VT(I,J)/VVVI
      VPE(I,J)=VPE(I,J)/VVVI
      FCT(I,J)=FCT(I,J)-FCVI*VT(I,J)
      FCPE(I,J)=FCPE(I,J)-FCVI*VPE(I,J)
      PET(I,J)=PET(I,J)-PEVI*VT(I,J)
      PEPE(I,J)=PEPE(I,J)-PEVI*VPE(I,J)
350   CONTINUE
C
C TIME
      CALL CLOCK
C
C CONVECTIVE FLUX
      DO 360 I=1,NTAU
      RT(I)=RT(I)-RFC(I)
      DO 361 J=1,NTAU
      TTT(I,J)=TTT(I,J)-FCT(I,J)
      TPE(I,J)=TPE(I,J)-FCPE(I,J)
361   CONTINUE
      IF(I.GT.MIHAL) GO TO 360
      RT(I)=RT(I)+RFC(I+1)
      DO 362 J=1,NTAU
      TTT(I,J)=TTT(I,J)+FCT(I+1,J)
      TPE(I,J)=TPE(I,J)+FCPE(I+1,J)
362   CONTINUE
360   CONTINUE
C
C ELECTRON PRESSURE
      CALL MATINV(PEPE,NTAU)
      CALL MULT(PET,PEPE,PET,SCRATC,NTAU,NTAU)
      CALL MULT(RPE,PEPE,RPE,SCRATC,NTAU,1)
      DO 374 I=1,NTAU
      DO 374 J=1,NTAU
      SUMA=0.
      DO 375 L=1,NTAU
      SUMA=SUMA+TPE(I,L)*PET(L,J)
375   CONTINUE
      TTT(I,J)=TTT(I,J)-SUMA
      RT(I)=RT(I)-TPE(I,J)*RPE(J)
374   CONTINUE
C
C TIME
      CALL CLOCK
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C BACKSUBSTITUTION IN GAUSS ELIMINATION SCHEME.
C
C INITIATE
      CALL MATINV(TTT,NTAU)
      DO 400 I=1,NTAU
      PE(I)=RPE(I)
      FC(I)=RFC(I)
      V(I)=RV(I)
      PT(I)=0.
      PR(I)=RPR(I)
      P(I)=RP(I)
C
C SOLVE FOR TEMPERATURE CORRECTION
      T(I)=0.
      DO 400 J=1,NTAU
      T(I)=T(I)+TTT(I,J)*RT(J)
400   CONTINUE
C      WRITE(2,68) ITER
C      WRITE(2,67) (T(I),I=1,NTAU)
      WRITE(6,68) IT
      WRITE(6,*) ' TEMPERATURE CORRECTION WANTED BY MARCS: '
      WRITE(6,67) (T(I),I=1,NTAU)
C---

C If the maximum temperature correction the program suggest to any model layer
C is smaller than TCONV (from input), then we will finish with computing one 
C more iteration where we will write the output etc. 
C UGJ 25/11/2004:
C In order to save time in interpolating in two-dimensional (i.e. atomic) 
C absorption coefficient files and in order not to compute a new atomic
C absorption coefficient file for each iteration, we start with a one-dimesional
C atomic file for a model with near-by parameters, and then recompute a new 
C (possible temporary) 1D atomic file for the almost converged model and use that 
C in the final iteration (if this cause a total sum of TCORMX of the following 
C iterations to be larger than 20.*TCONV then we repeat and make one more 
C atomic file for the more coverged model). In this way the opacity can be 
C calculated in the same way for the molecular (Doppler broadening only, and
C hence 1D absorption coefficients, dependent on temperature only) and the 
C atomic (T,Pg dependent Vogit profiles) absorption coefficient.

      TCORMX=ABS(T(1))
      DO 405 I=2,NTAU
 405  TCORMX=MAX(TCORMX,ABS(T(I)))
         if(newosatom.eq.2) then 
                     tcormxend = tcormxend + tcormx
         end if
      IF(TCORMX.LE. 5.*TCONV )  THEN 
         if(newosatom.eq.1) then 
C                    call osatom
                     newosatom = 2
                     tcormxend = 0.
         end if
      END IF
      IF(TCORMX.LE. TCONV )  THEN 
         ITSTOP=.TRUE.
C        if(newosatom.eq.2 .and. tcormxend.gt.20.0*tconv) call osatom
      END IF
      PRINT406, TCORMX,IT
406   FORMAT(' Max corr. to T wanted was',F6.1,' K for iteration',I2)
C
C
C     CALL DISPLA(' TCORMX ',TCORMX)
C---
C
C CHECK T CORR
      DO 401 I=1,NTAU
      T(I)=T(I)/SQRT(1.+25.*(T(I)/TT(I))**2)
401   CONTINUE
      TCORMX=ABS(T(1))
      DO 407 I=2,NTAU
 407  TCORMX=MAX(TCORMX,ABS(T(I)))
      PRINT4061, TCORMX
4061  FORMAT(' Max corr. to T wanted for kort=1 was',F6.1)
C
      IF (KORT.GE.4.AND.IT.GE.4  .OR. KORT.GE.5) THEN
      PPK = DFLOAT(KPP)
      DO 4011 I=1,NTAU
      if (abs(T(I)).GT.abs(PPK))
     - T(I)=PPK*t(i)/abs(t(i))
4011  CONTINUE
C
      TCORMX=ABS(T(1))
      DO 4071 I=2,NTAU
4071  TCORMX=MAX(TCORMX,ABS(T(I)))
      PRINT408, TCORMX,KORT
408   FORMAT(' Maximum correction applied was',F6.1,
     *     ' for applied kort =',I2)
      END IF
      
!      T(1:ntau) = 0.75*T(1:ntau)
      WRITE(6,*) ' TEMPERATURE CORRECTION APPLIED WAS: '
      WRITE(6,67) (T(I),I=1,NTAU)
C
C SUBTRACT TEMPERATURE
      DO 410 I=1,NTAU
      PT(I)=PT(I)-PTT(I)*T(I)
      P(I)=P(I)-PPTT(I)*T(I)
      IF(I.GT.1) P(I)=P(I)-PPTT(I+NTAU-1)*T(I-1)
      DO 410 J=1,NTAU
      V(I)=V(I)-VT(I,J)*T(J)
      PE(I)=PE(I)-PET(I,J)*T(J)
      FC(I)=FC(I)-FCT(I,J)*T(J)
      PR(I)=PR(I)-PRT(I,J)*T(J)
410   CONTINUE
C
C SUBTRACT ELECTRON PRESSURE
      DO 420 I=1,NTAU
      PT(I)=PT(I)-PTPE(I)*PE(I)
      P(I)=P(I)-PPPE(I)*PE(I)
      IF(I.GT.1) P(I)=P(I)-PPPE(I+NTAU-1)*PE(I-1)
      DO 420 J=1,NTAU
      V(I)=V(I)-VPE(I,J)*PE(J)
      FC(I)=FC(I)-FCPE(I,J)*PE(J)
420   CONTINUE
C
C TOTAL PRESSURE
      P(1)=P(1)+PR(1)-RPR(1)
      DO 425 I=2,NTAU
      P(I)=P(I)+P(I-1)
425   CONTINUE
C
C SUBTRACT VELOCITY
      DO 430 I=1,NTAU
      PT(I)=PT(I)-PTV(I)*V(I)
430   CONTINUE
C
C PRINT
      IF(PFE.OR.PFD) WRITE(7,48) FNORD,ITER,idmodl
      IF(PFE.OR.PFD) WRITE(7,56)
      IF(PFE.OR.PFD) WRITE(7,52)
      DO 440 I=1,NTAU
C
C SOLVE FOR CONVECTIVE EFFICIENCY AND GRADIENT DIFFERENCE
      I1=MAX0(I-1,1)
      GI=-.5*(GT(I)*(T(I)+T(I1))+GPE(I)*(PE(I)+PE(I1)))-GV(I)*V(I)
      GG(I)=GG(I)+GI
      D(I)=D(I)-DG(I)*GI-DPS(I)*P(I)-DPES(I)*PE(I)-DTS(I)*T(I)
      IF(I.GT.1) D(I)=D(I)-DPS(I+NTAU-1)*P(I-1)-DPES(I+NTAU-1)
     &*PE(I-1)-DTS(I+NTAU-1)*T(I-1)
      IF(PFE.OR.PFD)WRITE(7,51)I,TAU(I),PR(I),PT(I),P(I),GI,D(I),V(I)
     * ,FC(I),PE(I),T(I),I
440   CONTINUE
C
C TIME
      CALL CLOCK
C
C APPLY CORRECTIONS
      IPRESS=0
      DO 450 I=1,NTAU
      TT(I)=TT(I)+T(I)
      VV(I)=MAX(VV(I)+V(I),0.0D+0)
      FFC(I)=FFC(I)+FC(I)
      PPT(I)=MAX(PPT(I)+PT(I),0.0D+0)
      PPT(I)=MIN(PPT(I),0.5*PP(I))
      PPR(I)=PPR(I)+PR(I)
      DD(I)=DD(I)+D(I)
C
C IF TOO VIOLENT CHANGES TO PPE OR PP, SET IPRESS FOR AN EXTRA
C PRESSURE INTEGRATION AFTER CORRECTIONS HAVE BEEN APPLIED.
      IF(ABS(P(I)/PP(I)).LT.0.5.AND.ABS(PE(I)/PPE(I)).LT.0.5) GOTO 451
      IPRESS=1
      GOTO 450
C
451   PPE(I)=PPE(I)+PE(I)
      PP(I)=PP(I)+P(I)
450   CONTINUE
      IF(IPRESS.EQ.1) CALL TRYCK

C
C PRINT PRESENT STATE OF THE ATMOSPHERE
      ENTRY PRESNT
      IF(.NOT.PFE) RETURN
      INORD=IEDIT+10*IVERS
      FNORD=.1*INORD
      WRITE(7,48) FNORD,ITER,idmodl
      WRITE(7,53)
      WRITE(7,52)
      I=0
      Y=0.
      TT0=(TAU(2)*TT(1)-TAU(1)*TT(2))/(TAU(2)-TAU(1))
      WRITE(7,51) I,Y,PPR(1),Y,PPR(1),Y,Y,Y,Y,Y,TT0,I
      WRITE(7,51) (I,TAU(I),PPR(I),PPT(I),PP(I),GG(I),DD(I)
     * ,VV(I),FFC(I),PPE(I),TT(I),I,I=1,NTAU)
C
C TIME
      CALL CLOCK
      RETURN
C
C FORMATS
45    FORMAT(' TIME',I6,' MSEC')
48    FORMAT('1MARCS',F5.1,5X,'SOLVE(61) 12-NOV-80',5X,
     & '                   ',5X,'ITERATION',I3,5X,A24/)
49    FORMAT(13('1234567890'),'123')
50    FORMAT(4(/2X,1P10E13.5))
51    FORMAT(I3,1P10E12.4,I4)
52    FORMAT(7X,'TAU',9X,'PRAD',8X,'PTURB',7X,'PTOT',8X,'GAMMA',
     & 7X,'DELTA',7X,'VCONV',7X,'FCONV',7X,'PE',10X,'TEMP')
53    FORMAT(' STATE OF MODEL ATMOSPHERE')
54    FORMAT(' *=R.H. SIDES     *',23X,'*',23X,'*',
     ,11X,'*',11X,'*',11X,'*',11X,'*')
55    FORMAT(6X,'TAU',9X,'HSCALE',6X,'ADIA',8X,'GRAD',8X,'CP',10X,
     *'Q',11X,'PG',10X,'RO',10X,'PGPE',8X,'PGT')
56    FORMAT(' CORRECTIONS')
57    FORMAT(' THERMODYNAMICALS')
58    FORMAT(1X,I3,F9.1,F8.1,1P2E10.3,0PF6.3,1PE10.3,0PF8.3,
     * F7.0,4F5.1,1PE10.3,2(0PF6.1),2I5,I4)
59    FORMAT('   J',4X,'WAVEL',2X,'WEIGHT',2X,'F(WAVEL)',3X,
     * 'CONTRIB',1X,'WAVEN',2X,'F(WAVEN)',4X,'MAGN',2X,'TRAD',3X,
     * 'X01',2X,'X25',2X,'S01',2X,'S25',3X,'CONTRIB',2X,
     * 'DX01  DX25 OPAC TRAN ALGB')
60    FORMAT(' SAVED VALUES FROM LOGICAL UNIT',I3)
61    FORMAT(' STARTING VALUES, QTEMP=',F5.2)
62    FORMAT(' PRESSURE EQUATION')
63    FORMAT(6X,'TAU',9X,'PTAU',8X,'ROSS',8X,'ROSSPE',6X,'ROSST',7X,
     &'KAP5',8X,'ROSSP',7X,'TAU5',8X,'TAUP')
64    FORMAT(' K=',I2,4X,'NEW V =',1PE10.3)
65    FORMAT(' TOTAL FLUX=',1PE11.4,' ERGS/CM**2/S  TEFF=',0PF6.0,' K')
66    FORMAT(7X,F8.0)
67    FORMAT(1X,10F7.1)
68    FORMAT(' ITERATION',I3)
      END
C*
C*NEW PDS MEMBER FOLLOWS
C*
      SUBROUTINE TRANEQ
      implicit real*16 (a-h,o-z)
C
C TRANEQ SOLVES THE TRANSFER EQUATION INCLUDING CONTINUUM SCATTERING.
C FEATURES:
C
C 1. CANNONS PERTURBATION TECHNIQUE IS USED ON THE ANGULAR QUADRATURE.
C    THE BASIC IDEA IN THIS TECHNIQUE IS TO REPLACE THE INVERSION OF
C    A COMPLICATED (MMU_PP ORDER) OPERATOR WITH THE INVERSION OF A SIMPLE
C    OPERATOR (ONE POINT=EDDINGTON APPROXIMATION), PLUS ITERATION ON
C    THE ERROR.
C 2. AITKEN EXTRAPOLATION ACCELLERATES THE CONVERGENCE.
C 3. A TRICK DUE TO ROBERT STEIN (PRIV. COMM., 1979) IS USED TO
C    ELIMINATE THE NEED FOR DOUBLE PRECISION STORAGE OF THE MATRIX
C    ELEMENTS. THE IDEA IS TO STORE THE (SMALL) SUM OF THE THREE
C    MATRIX ELEMENTS ON A ROW, INSTEAD OF THE (LARGE) DIAGONAL ELE-
C    MENT.
C 4. THE SOLUTION IS A CUBIC SPLINE, RATHER THAN A PIECE-WISE
C    QUADRATIC FUNCTION. THIS IS ACCOMPLISHED WITH THE CORRECTION
C    TERMS AD AND BD IN SUBROUTINE TRANFR.
C 5. THE SCATTERING IS TREATED AS DIPOLE SCATTERING INSTEAD OF THE NORMA
C    USED ISOTROPIC APPROXIMATION. THIS CAN BE DONE VERY SIMPLY IN THE
C    ITERATING CANNON SCHEME.
C 6. A BOUNDARY CONDITION WHICH INCLUDES AN ESTIMATED INFALLING
C    RADIATION MAKES THE SOLUTION GOOD ALSO FOR VALUES OF X+S
C    LARGE COMPARED WITH 1./TAU(1). A LOGARITHMIC TAU-SCALE
C    SHOULD BE USED.
C
C THIS VERSION OF TRANEQ IS COMPATIBLE WITH PREVIOUS TRANEQ'S.
C 79.06.21 *NORD*
C
      include 'parameter.inc'
      COMMON /TAUC/TAU(NDP),DTAULN(NDP),JTAU
      COMMON /CTRAN/X(NDP),S(NDP),BPLAN(NDP),XJ(NDP),XH(NDP),XK(NDP)
     &  ,dumtran(4*ndp),idumtran(3)
      COMMON /SPACE2_PP/SOURCE(NDP),ERROR(NDP),DUM(3*NDP),P(NDP)
     & ,SP1(NDP,6),SP2(NDP,6),SP3(NDP,6),AN(NDP),AD(NDP),BD(NDP)
     & ,FACT(NDP),DSO(NDP),SP2DUM((4*NDP-29)*NDP)
      DIMENSION A(7)
C
C INITIATE
      DO 100 K=1,JTAU
      FACT(K)=1.
      DSO(K)=0.
      XJ(K)=0.
      XK(K)=0.
      ERROR(K)=BPLAN(K)*X(K)/(X(K)+S(K))
100   SOURCE(K)=0.
C
C CALCULATE THE MATRIX ELEMENTS
      CALL TRANFR
      CALL TRANSC
C
C ITERATION LOOP
      ITMAX=7
      DO 110 IT=1,ITMAX
110   A(IT)=0.
      DO 140 IT=1,ITMAX
      ITM=IT
C
C SOLVE THE CONTINUUM SCATTERING PROBLEM IN THE EDDINGTON APPROXIMATION
      CALL SCATTR
      DO 120 K=1,JTAU
      XJ(K)=XJ(K)+P(K)
      XK(K)=XK(K)+.333333*P(K)
C
C AITKEN EXTRAPOLATION USED FOR CONVERGENCE ACCELLERATION
      DS=ERROR(K)+P(K)*S(K)/(X(K)+S(K))
      IF(DSO(K).NE.0.) 
     #   FACT(K)=MIN(1.25D+0,MAX(0.8D+0,FACT(K)-DS/DSO(K)))
      DS=DS/FACT(K)
      IF(IT.GE.2) DSO(K)=DS
120   SOURCE(K)=SOURCE(K)+DS
C
C SOLVE THE TRANSFER EQUATION WITH GIVEN SOURCE FUNCTION
      CALL FORMAL
C
C CHECK ERROR IN SOURCE FUNCTION
      DO 130 K=1,JTAU
130   A(IT)=MAX(A(IT),ABS(ERROR(K)/SOURCE(K)))
C
C END OF ITERATION LOOP
      IF(A(IT).LT.0.0001) GO TO 141
140   CONTINUE
      WRITE(6,50) (A(IT),IT=1,ITM)
50    FORMAT(' MAXFEL =',12F10.7)
141   CONTINUE
C
      RETURN
      END
C*
C*NEW PDS MEMBER FOLLOWS
C*
      SUBROUTINE TRANFR
      implicit real*16 (a-h,o-z)
C
C FORMAL SOLVES THE TRANSFER EQUATION WITH GIVEN SOURCE FUNCTION 'SOURCE
C 'ERROR' IS THE RESULTING ERROR IN THE DEFINITION OF THE CONTINUUM
C SCATTERING SOURCE FUNCTION. TRANSFR CALCULATES THE MATRIX ELEMENTS
C OF THE PROBLEM. FLUX AND INTENSITIES AT TAU=0 ARE RETURNED IN /CSURF/.
C 79.06.21 *NORD*
C
      include 'parameter.inc'
      COMMON /CANGLE/XMU(6),XMU2(6),H(6),MMU_PP
      COMMON /CTRAN/X(NDP),S(NDP),BPLAN(NDP),XJ(NDP),XH(NDP),XK(NDP)
     &  ,dumtran(4*ndp),idumtran(3)
      COMMON /TAUC/TAU(NDP),DTAULN(NDP),JTAU
      COMMON /SPACE2_PP/SOURCE(NDP),ERROR(NDP),DUM(3*NDP),P(NDP)
     *,SP1(NDP,6),SP2(NDP,6),SP3(NDP,6),AN(NDP),AD(NDP),BD(NDP)
     *,FACT(NDP),DSO(NDP),C(6),T(6),EX(6),SP2DUM((4*NDP-29)*NDP-18)
      COMMON /CSURF/HSURF,Y1(NRAYS)
C
C MU LOOP
      JTAU1=JTAU-1
      JTAU2=JTAU-2
      IF (mmu_pp.GT.nrays) STOP
     &    ' tranfr:increase nrays in parameter.inc'
      IF (mmu_pp.GT.6) STOP ' tranfr: increase dimensions incl mmu_pp'
      DO 110 I=1,MMU_PP
C
C K=1
      DTAUB=.5*(X(1)+S(1)+X(2)+S(2))*(TAU(2)-TAU(1))/XMU(I)
      A=1./DTAUB
      B=A**2
      SP2(1,I)=1.+2.*A
      SP3(1,I)=-2.*B
C LET P BE THE EVEN PART OF THE INTENSITY, THEN
C
C         P(2)= P(1) + D*P'(1) + .5*D2*P''(1)
C OR      P(2)= P(1) + D*(P(1)-I(1,-MU)) + .5*D2*(P(1)-S(1)) .
C WHERE   I(1,-MU) = S(1)*(1.-EXP(-T))
C
C THE DIFFERENCE AS COMPARED TO THE USUAL SECOND ORDER BOUNDARY CONDITIO
C IS THE ADDITIONAL TERM   I(1,-MU)=S(1)*(1.-EXP(-T)). THUS THE COEFFICI
C FOR S(1) IN THE FIRST EQUATION SHOULD BE CHANGED AS FOLLOWS
C         S(1)=S(1)*(1.+C*(1.-EXP(-T))
C WHERE   C=2./D
C *NORD* 751009
      C(I)=2.*A
      T(I)=TAU(1)*(X(1)+S(1))/XMU(I)
      IF(T(I).LE.0.1) THEN
        EX(I)=T(I)*(1.-.5*T(I)*(1.-.3333*T(I)))
      ELSE IF(T(I).LE.675.0) THEN
        EX(I)=1.0-EXP(-T(I))
      ELSE IF(T(I).GT.675.0) THEN
        EX(I)=1.0
      END IF
C
C K=2,JTAU-1
      DO 100 K=2,JTAU1
      DTAUA=DTAUB
      DTAUB=.5*(X(K)+S(K)+X(K+1)+S(K+1))*(TAU(K+1)-TAU(K))/XMU(I)
      DTAUC=.5*(DTAUA+DTAUB)
      AD(K)=.166667*DTAUA/DTAUC
      BD(K)=.166667*DTAUB/DTAUC
      SP1(K,I)=-1./(DTAUA*DTAUC)+AD(K)
      SP2(K,I)=1.
100   SP3(K,I)=-1./(DTAUB*DTAUC)+BD(K)
C
C K=JTAU
      SP2(JTAU,I)=1.
C
C END OF MU LOOP
110   CONTINUE
C
C ELIMINATE SUBDIAGONAL, SAVE FACTORS IN SP1
      DO 121 I=1,MMU_PP
      DO 120 K=1,JTAU2
      SP1(K,I)=-SP1(K+1,I)/(SP2(K,I)-SP3(K,I))
      SP2(K+1,I)=SP2(K+1,I)+SP1(K,I)*SP2(K,I)
120   SP2(K,I)=SP2(K,I)-SP3(K,I)
121   SP2(JTAU-1,I)=SP2(JTAU-1,I)-SP3(JTAU-1,I)
C
      RETURN
C
      ENTRY FORMAL
C
C ZEROSET
      DO 130 K=1,JTAU
      AN(K)=(3.*XK(K)-XJ(K))/8.*S(K)/(X(K)+S(K))
      XK(K)=0.
130   XJ(K)=0.
C
C MU LOOP
      XH(1)=0.
      HSURF=0.
      DO 170 I=1,MMU_PP
C
C INITIATE APPROXIMATIVE SOURCE FUNCTION
      P(1)=SOURCE(1)+AN(1)*(3.*XMU2(I)-1.)
C NOTE THE ANISOTROPIC SCATTERING CORRECTION
      S0=P(1)
      P(1)=P(1)*(1.+C(I)*EX(I))
      DO 140 K=2,JTAU1
140   P(K)=(1.-AD(K)-BD(K))*(SOURCE(K)+AN(K)*(3.*XMU2(I)-1.))
     & +AD(K)*(SOURCE(K-1)+AN(K-1)*(3.*XMU2(I)-1.))
     & +BD(K)*(SOURCE(K+1)+AN(K+1)*(3.*XMU2(I)-1.))
      P(JTAU)=SOURCE(JTAU)
C
C ACCUMULATE RIGHT HAND SIDE
      DO 150 K=1,JTAU2
150   P(K+1)=P(K+1)+SP1(K,I)*P(K)
C
C BACKSUBSTITUTE
      DO 160 K=1,JTAU1
      P(JTAU-K)=(P(JTAU-K)-SP3(JTAU-K,I)*P(JTAU-K+1))/SP2(JTAU-K,I)
      XK(JTAU-K)=XK(JTAU-K)+H(I)*P(JTAU-K)*XMU2(I)
160   XJ(JTAU-K)=XJ(JTAU-K)+H(I)*P(JTAU-K)
C
C END OF MU LOOP
      XK(JTAU)=XK(JTAU)+H(I)*P(JTAU)*XMU2(I)
      R1=P(1)-S0*EX(I)
      XH(1)=XH(1)+H(I)*XMU(I)*R1
      P0=P(1)*(1.-EX(I))+.5*S0*EX(I)**2
      HSURF=HSURF+H(I)*XMU(I)*P0
      Y1(I)=2.*P0
C HSURF AND Y1(6) ARE THE FLUX AND INTENSITIES AT THE SURFACE
170   CONTINUE
      XJ(JTAU)=P(JTAU)
C
C 'XJ' IS THE NEW MEAN INTENSITY
      DO 180 K=1,JTAU
180   ERROR(K)=(X(K)*BPLAN(K)+S(K)*XJ(K))/(X(K)+S(K))-SOURCE(K)
C
C FLUX AND SECOND MOMENT
      DO 190 K=2,JTAU
190   XH(K)=2.*(XK(K)-XK(K-1))/(X(K)+S(K)+X(K-1)+S(K-1))/
     /(TAU(K)-TAU(K-1))
C
      RETURN
      END
C*
C*NEW PDS MEMBER FOLLOWS
C*
      SUBROUTINE TRANSC
      implicit real*16 (a-h,o-z)
C
C SCATTR SOLVES THE TRANSFER EQUATION INCLUDING CONTINUUM SCATTERING
C IN THE EDDINGTON APPROXIMATION, I.E., USING ONLY ONE MU POINT.
C 'ERROR' IS THE INHOMOGENEOUS TERM OF THE EQUATION, AND 'P' CONTAINS
C THE ESTIMATED MEAN INTENSITY ON EXIT. TRANSC CALCULATES THE MATRIX
C ELEMENTS FOR SCATTR.
C 79.06.21 *NORD*
C
      include 'parameter.inc'
      COMMON /CTRAN/X(NDP),S(NDP),BPLAN(NDP),XJ(NDP),HFLUX(NDP),XK(NDP)
     &  ,dumtran(4*ndp),idumtran(3)
      COMMON /TAUC/TAU(NDP),DTAULN(NDP),JTAU
      COMMON /SPACE2_PP/SOURCE(NDP),ERROR(NDP),SP1(NDP),SP2(NDP),
     +            SP3(NDP),P(NDP),SP2DUM((4*NDP-6)*NDP)
      DATA XMU,XMU2/0.5773503,0.3333333/
C
C K=1
      DTAUB=.5*(X(1)+S(1)+X(2)+S(2))*(TAU(2)-TAU(1))/XMU
      A=1./DTAUB
      B=A**2
      SP2(1)=2.*A+X(1)/(S(1)+X(1))
      SP3(1)=-2.*B
      C=2.*A
      T=TAU(1)*(X(1)+S(1))/XMU
      EX=T*(1.-.5*T*(1.-.33333*T))
      IF(T.GT.0.1) EX=1.-EXP(-T)
      SP2(1)=SP2(1)-C*S(1)/(X(1)+S(1))*EX
C
C K=2,JTAU-1
      JTAU1=JTAU-1
      DO 100 K=2,JTAU1
      DTAUA=DTAUB
      DTAUB=.5*(X(K)+S(K)+X(K+1)+S(K+1))*(TAU(K+1)-TAU(K))/XMU
      DTAUC=.5*(DTAUA+DTAUB)
      A=1./(DTAUA*DTAUC)
      B=1./(DTAUB*DTAUC)
      SP1(K)=-A
      SP2(K)=X(K)/(S(K)+X(K))
      SP3(K)=-B
100   CONTINUE
C
C K=JTAU
      SP2(JTAU)=X(JTAU)/(X(JTAU)+S(JTAU))
C
C ELIMINATE SUBDIAGONAL
      JTAU2=JTAU-2
      DO 110K=1,JTAU2
      SP1(K)=-SP1(K+1)/(SP2(K)-SP3(K))
      SP2(K+1)=SP2(K+1)+SP1(K)*SP2(K)
110   SP2(K)=SP2(K)-SP3(K)
      SP2(JTAU-1)=SP2(JTAU-1)-SP3(JTAU-1)
      RETURN
C
      ENTRY SCATTR
C
C INITIATE INHOMOGENOUS TERMS
      DO 120 K=1,JTAU
120   P(K)=ERROR(K)
      DSDT=0.
C PRELIM
      P(1)=P(1)*(1.+C*EX)-DSDT*C*(EX-T*(1.-EX))
C
C ACCUMULATE INHOMOGENOUS TERMS
      DO 130 K=1,JTAU2
130   P(K+1)=P(K+1)+SP1(K)*P(K)
C

C BACKSUBSTITUTE
      P(JTAU)=P(JTAU)/SP2(JTAU)
      DO 140 K=1,JTAU1
140   P(JTAU-K)=(P(JTAU-K)-SP3(JTAU-K)*P(JTAU-K+1))/SP2(JTAU-K)
C
      RETURN
      END
C*
C*NEW PDS MEMBER FOLLOWS
C*
      SUBROUTINE TRYCK
      implicit real*16 (a-h,o-z)
C
C TRYCK IS A FAST PRESSURE INTEGRATION ROUITINE. IT IS FAST BECAUSE OF
C TWO REASONS: 1) IT INTEGRATES THE DIFFFERENTIAL EQUATION FOR LN(P)
C AS A FUNCTION OF LN(TAU). 2) IT ITERATES DIRECTLY ON THE ELECTRON
C PRESSURE, KEEPING THE NUMBER OF CALLS TO ABSKO TO A MINIMUM.
C ASSUMING A POWER LAW BEHAVIOUR OF PP,TT,PPE,ETC.: PP=C*TAU**DP,ETC., O
C CAN SHOW THAT DP=(1.+DT*(ROSSPE*PGT/PGPE-ROSST)/(1.+ROSSPE/PGPE).
C THE ANSATZ FOR PP IMPLIES PP(1)=TAU(1)*GRAV/(ROSS(1)*DP), WHICH
C SERVES AS A BOUNDARY CONDITION.
C 790516 *NORD*
C
C TRYCK/KOL IS A VERSION THAT USES CROSS(K) AS A FACTOR ON ROSSOP.
C 801105 *NORD*
C
      include 'parameter.inc'
      COMMON /STATEC/PPR(NDP),PPT(NDP),PP(NDP),GG(NDP),ZZ(NDP),DD(NDP),
     & VV(NDP),FFC(NDP),PPE(NDP),TT(NDP),TAULN(NDP),RO(NDP),NTAU,ITER
      COMMON /TAUC/TAU(NDP),DLNTAU(NDP),JTAU 
      COMMON /CG/GRAV,KONSG
      COMMON /ROSSC/ROSS(NDP),CROSS(NDP)
      COMMON /CROSSOS/ ROSSO(NDP),PTAUO(NDP)
      COMMON /CI8/PGC,RHOC,EC
      COMMON /CSPHER/TAURAT,RADIUS,RR(NDP),NCORE
      COMMON /CMOLRAT/ FOLD(NDP,8),MOLOLD,KL
      DATA EPS,RELT,RELPE,PEDEF/1.E-3,1.E-3,1.E-3,1./
CUGJ98DATA EPS,RELT,RELPE,PEDEF/1.E-3,1.E-3,1.E-3,1./

C
C START
C     CALL MSLEFT(MSA)
      MSA=0
C      WRITE(7,101)
C101   FORMAT('1PRESSURE INTEGRATION'/'  K',6X,'TAU',10X,'TT',
C     & 9X,'PPE',8X,'PTOT',8X,'ROSS',10X,'DP',9X,'NABSKO')
      DT=0.
C USE 'DLNT/DLNTAU'=DT=0. TO BE COMPATIBLE WITH SOLVE. OTHERWISE
C DT=(TT(2)/TT(1)-1.)/DLNTAU(2)
      NABSKO=0
      KK=1
      IF(PPE(1).LE.0.) PPE(1)=PEDEF

C
C ITERATE ON BOUNDARY CONDITION, USING PARTIAL DERIVATIVES
100   CONTINUE
      KL=1
C      call rossos
      ROSS(1)=CROSS(1)*ROSSOP(TT(1),PPE(1))
C      ross(1)=rosso(1)
      PP(1)=PGC+PPT(1)+PPR(1)
      PG=PGC
C      tt(1) = tt(1)*(1.+relt)
C      call rossos
C      tt(1) = tt(1)/(1.+relt)
      ROSST=CROSS(1)*ROSSOP(TT(1)*(1.+RELT),PPE(1))
C      rosst = rosso(1)
      PGT=PGC
C      ppe(1) = ppe(1)*(1.+relpe)
C      call rossos
C      ppe(1) = ppe(1)/(1.+relpe)
      ROSSPE=CROSS(1)*ROSSOP(TT(1),PPE(1)*(1.+RELPE))
C      rosspe = rosso(1)
      PGPE=PGC
C      write(7,*) ' ross(1),rosst,rosspe,pg,pgt,pgpe = '
C      write(7,*) ross(1),rosst,rosspe,pg,pgt,pgpe

      PGT=(PGT/PG-1.)/RELT
      PGPE=(PGPE/PG-1.)/RELPE
      ROSST=(ROSST/ROSS(1)-1.)/RELT
      ROSSPE=(ROSSPE/ROSS(1)-1.)/RELPE
      NABSKO=NABSKO+3
      DP=(1.+DT*(ROSSPE*PGT/PGPE-ROSST))/(1.+ROSSPE/PGPE)
      DP=MAX(DP,0.1D+0)
      testz = GRAV*TAU(1)/(PG*ROSS(1)*DP)
C      if (testz.le.0.d0 .or. pgpe*pgt.le.0.) then
C        write(6,*) ' in tryck with problems:'
C        write(6,*) 'pg,ross(1),dp,pp(1),tt(1),ppe(1),'
C    *  ,'testz,1./(pgpe+rosspe),cross(1),pgpe,pgt,rosspe,rosst'
C23456789 12345678901234567890123456789012345678901234567890123456789 12
C        write(6,1212) pg,ross(1),dp,pp(1),tt(1),ppe(1)
C    *  ,testz,1./(pgpe+rosspe),cross(1),pgpe,pgt,rosspe,rosst
1212     format(1p13e11.3)
C      end if
      DLNPE=log(GRAV*TAU(1)/(PG*ROSS(1)*DP))/(PGPE+ROSSPE)
C      write(7,*) ' k-tau=1: dp,dlnpe = ',dp,dlnpe
      PPE(1)=PPE(1)*EXP(DLNPE)
      IF(ABS(DLNPE).GT.EPS) GOTO 100
C
C END BOUNDARY CONDITION
C      call rossos
      ROSS(1)=CROSS(1)*ROSSOP(TT(1),PPE(1))
C      write(7,*) 'final ross(1),rosso(1),cross(1) ='
C      write(7,*) ross(1),rosso(1),cross(1)
C      ross(1)=rosso(1)
      NABSKO=NABSKO+1
      PP(1)=PGC+PPT(1)+PPR(1)
C2021  format(' final PP(1),PGC,PPT(1),PPR(1):',1p4e12.3)
C      write (7,2021) PP(1),PGC,PPT(1),PPR(1)
C      WRITE(7,102) KK,TAU(1),TT(1),PPE(1),PP(1),ROSS(1),DP,NABSKO
C102   FORMAT(I3,6E12.5,I12)
C
C TAU LOOP
      DPE=(DP-DT*PGT)/PGPE
      DEDLNP=-(PGPE*PG/PP(1)+.5*DLNTAU(2)*GRAV*TAU(1)/(PP(1)*ROSS(1))*
     & (PGPE*PG/PP(1)+ROSSPE))
C      write(7,*) ' before ntau loop 110: dpe,dedlnp =',dpe,dedlnp
C      write(7,2022)
C2022  format ('k nabsko ppe(k) error dlnpe pgc',
C     *   ' pp(k) ross(k) rosso(k) cross(k)')

      DO 110 K=2,NTAU
      PPE(K)=PPE(K-1)*EXP(DPE*DLNTAU(K))
      NABSKO=0
C
C ITERATION LOOP
      DLNPE=0.
111   CONTINUE
      KL=K
C      call rossos
      ROSS(K)=CROSS(K)*ROSSOP(TT(K),PPE(K))
!      ross(k) = max(ross(k),1e-5)
C      rossim = ross(k)
C      ross(k)=rosso(k)
      PP(K)=PGC+PPT(K)+PPR(K)
      NABSKO=NABSKO+1
      ERROR=(.5*DLNTAU(K)*GRAV*(TAU(K-1)/(PP(K-1)*ROSS(K-1))+
     & TAU(K)/(PP(K)*ROSS(K)))-log(PP(K)/PP(K-1)))
      CALL ZEROF(ERROR,DLNPE,DEDLNP)
      PPE(K)=PPE(K)*EXP(DLNPE)
C      write(7,2023) k,nabsko,ppe(k),error,dlnpe,pgc,pp(k)
C     *    ,rossim,rosso(k),cross(k)
C2023  format(2i4,1p7e12.3,0pf8.3)
      IF(ABS(DLNPE).GT.EPS) GOTO 111
C
C END TAU LOOP
C      call rossos
      ROSS(K)=CROSS(K)*ROSSOP(TT(K),PPE(K))
C      ross(k)=rosso(k)
      NABSKO=NABSKO+1
      PP(K)=PGC+PPT(K)+PPR(K)
      DP=GRAV*TAU(K)/(PGC*ROSS(K))
      DPE=log(PPE(K)/PPE(K-1))/DLNTAU(K)
C      WRITE(7,102) K,TAU(K),TT(K),PPE(K),PP(K),ROSS(K),DP,NABSKO
110   CONTINUE

C      write(7,*) ' k,tt,ppe,pp,ross in tryck after iterations:'
C      DO 211 K=1,NTAU
C211   write(7,1171)k,tt(k),ppe(k),pp(k),ross(k)
C
C END
C     CALL MSLEFT(MSB)
      MSB=0
      MSB=MSA-MSB
C      WRITE(7,120) MSB
C120   FORMAT('0TRYCK TIME=',I5,' MS')
      RETURN
      END
C
      FUNCTION TRQUAD(N,X,F,W)
      implicit real*16 (a-h,o-z)
C
      DIMENSION X(N),F(N),W(2*N)
C was : dim x(1) etc...
C
C TRAPEZOIDAL QUADRATURE PLUS NEXT ORDER CORRECTION FOR NON-
C -EQUIDISTANT GRID.
      N1=N-1
      Q=0.
      DO 100 K=2,N
      W(K)=X(K)-X(K-1)
      W(N+K)=(F(K)-F(K-1))/W(K)
100   Q=Q+W(K)*(F(K-1)+F(K))
      Q=Q*6.
      DO 101 K=2,N1
101   Q=Q+(W(K+1)-W(K))*(W(K)*W(N+K+1)+W(K+1)*W(N+K))
      W1=((W(2)+0.5*W(3))*W(N+2)-0.5*W(2)*W(N+3))*2.0/(W(2)+W(3))
      WN=((W(N)+0.5*W(N1))*W(N+N)-0.5*W(N)*W(N+N1))*2.0/(W(N)+W(N1))
      Q=0.083333333*(Q+W(2)**2*W1-W(N)**2*WN)
      TRQUAD=Q
      RETURN
      END
C
C
C
C_ugj950523:  here the sphereical part begins @@@@@
C
C
      SUBROUTINE SOLVE_sph(NEW)
      implicit real*16 (a-h,o-z)
C
C SOLVE PERFORMS ONE NEWTON-RAPSON ITERATION ON THE MODELATMOSPHERE PROBLEM
C INCLUDING LOCAL CONVECTION.THE STATE OF THE ATMOSPHERE IS DESCRIBED BY A
C NUMBER OF VARIABLES SUCH AS TEMPERATURE,ELECTRON PRESSURE,TOTAL PRESSURE,
C CONVECTIVE FLUX ETC..TO EACH VARIABLE CORRESPONDS A CERTAIN CONDITIONAL
C EQUATION WICH DETERMINES THAT VARIABLE,ASSUMING THE OTHER VARIABLES BEEING
C KNOWN.
C
C NAMING CONVENTION.THE VARIABLES HAVE NAMES WITH A DOUBLE OCCURANCE OF THE
C FIRST LETTER,CORRECTIONS (TO BEE COMPUTED IN THIS ITERATION) HAVE SINGLE
C OCCURANCE OF FIRST LETTER.RIGHTHANDSIDENAMES BEGIN WITH R.
C
C VARIABLES ARE CENTERED ON INTEGER AND HALFINTEGER TAU-POINTS AS INDICATED
C BY 'I' OR 'H' ON THE FOLLOWING COMMENT CARDS.
C
C VARIABLE CORRECTION                                           HALF-INTEGER
C PPR      PR         RADIATION PRESSURE                             I
C PPT      PT         TURBULENT PRESSURE                        H
C PP       P          TOTAL PRESSURE                                 I
C GG       G          CONVECTIVE EFFICIENCY,GAMMA               H
C ZZ       Z          GEOMETRIC HEIGTH                          H
C DD       D          GRADIENT DIFFERENCE,DELTA-DELTAPRIME      H
C VV       V          CONVECTIVE VELOCITY                       H
C FFC      FC         CONVECTIVE FLUX                           H
C PPE      PE         ELECTRON PRESSURE                              I
C TT       T          TEMPERATURE                                    I
C XJ                  MEAN INTENSITY                                 I
C
C NAMES OF THE PARTIAL DERIVATIVES ARE FORMED WITH A FIRST PART FROM THE
C EQUATION TO WICH IT BELONGS,AND A SECOND PART WICH IS THE NAME OF THE
C VARIABLE WITH RESPECT TO WICH THE DERIVATIVE IS TAKEN.(DOUBLE OCCURANCE
C IF NECESSARY TO AVOID CONFUSION).
C
C PRPR,PRT
C PTPT,PTV,PTPE,PTT
C PPP,PPPE,PPTT
C GGG,GV,GPE,GT
C ZZZ,ZPE,ZT
C DDD,DP,DG,DPE,DT
C VVV,VZ,VD,VPE,VT
C FCFC,FCP,FCD,FCV,FCPE,FCT
C PEPE,PEPR,PEPT,PEP,PET
C TXJ,TFC,TTT
C
C  PARAMETRIZED VERSION (in NDP and NRAYS).     B.PLEZ 20-NOV-88. 
C
      include 'parameter.inc'
C
C STATE VARIABLES
      COMMON /STATEC/PPR(NDP),PPT(NDP),PP(NDP),GG(NDP),ZZ(NDP),DD(NDP),
     &VV(NDP),FFC(NDP),PPE(NDP),TT(NDP),TAULN(NDP),ROSTAT(NDP),NTAU,ITER
      common /ckdtpe/dpex,kdtpe
C
C DIMENSIONS
      DIMENSION PTAU(NDP),ROSSP(NDP),SUMW(NDP),ROSST(NDP),ROSSPE(NDP)
     *,XL(500),W(500)
     *,XJ1(NDP),XJ2(NDP),XJ3(NDP),XJT1(NDP),XJT2(NDP),XJT3(NDP)
     *,XJPE1(NDP),XJPE2(NDP),XJPE3(NDP)
     *,PR(NDP),PRT(NDP,NDP),PRJ(NDP)
     *,PRPE(NDP,NDP)
     *,PT(NDP),PTV(NDP),PTPE(NDP),PTT(NDP)
     *,P(NDP),PPPE(2*NDP),PPTT(2*NDP)
     *,GV(NDP),GPE(NDP),GT(NDP)
     *,DP(2*NDP),DG(NDP),DPE(NDP,NDP),DT(NDP,NDP),DV(NDP)
      DIMENSION D(NDP),DTS(2*NDP),DPS(2*NDP),DPES(2*NDP)
     *,V(NDP),VD(NDP),VPE(NDP,NDP),VT(NDP,NDP)
     *,FC(NDP),FCD(NDP),FCV(NDP),FCPE(NDP,NDP),FCT(NDP,NDP)
     *,PE(NDP),PEPE(NDP,NDP),PET(NDP,NDP)
     *,T(NDP),TTT(NDP,NDP),TPE(NDP,NDP),TJ1(NDP),TJ2(NDP)
     *,SCRATC(NDP,NDP)
     *,X(NDP),S(NDP),BPLAN(NDP),DBPL(NDP),XJ(NDP),HFLUX(NDP),XK(NDP)
     *,XT(NDP),ST(NDP),DLNX(NDP),XLOG(NDP)
     *,XPE(NDP),SPE(NDP)
     *,RPR(NDP),RP(NDP),RD(NDP),RV(NDP),RFC(NDP),RPE(NDP),RT(NDP)
     *,TLAST(NDP),TVD(NDP),IA(5)
      LOGICAL NEWV
C
C CONNECTIONS VIA COMMON.
C THE COMMENTED COMMONS MUST BE INITIATED OUTSIDE THIS ROUTINE BEFORE IT
C IS CALLED.
C JTAU=NUMBER OF TAUPOINTS, TAU=TAUSCALE.
C NLAM=NUMBER OF LAMBDAPOINTS, XL=LAMBDAPOINTS, W=INTEGRATIONWEIGHTS.
C MIHAL=LOWER LIMIT OF RADIATIVE EQUILIBRIUM CONDITION, TAUMAX NOT USED.
C PALFA,PBETA,PNY,PNY = MIXING LENGTH THEORY COEFFICIENTS.
C GRAV=SURFACE GRAVITY, TEFF=EFFECTIVE TEMPERATURE, FLUX=STEFAN*TEFF**4/PI
      CHARACTER MOLNAME*4,OSFIL*60,SAMPLING*3
      COMMON/COS/WNOS(NWL),CONOS(NDP,NWL),WLOS(NWL),WLSTEP(NWL)
     *    ,KOS_STEP,NWTOT,NOSMOL,NEWOSATOM,NEWOSATOMLIST
     *    ,nchrom,OSFIL(30),MOLNAME(30),SAMPLING
      COMMON /CLEVETAT/GEFF(NDP),PPRG(NDP),AMLOSS
      COMMON /CLEVPRINT/ PRJ2(NDP),masslinf
      COMMON /CMOLRAT/ FOLD(NDP,8),MOLOLD,KL
      common /CPRINT/NPRINT
      COMMON /TAUC/TAU(NDP),DTAULN(NDP),JTAU
      COMMON /CVAAGL/XL,W,NLAM
      COMMON /CSTYR/MIHAL,NOCONV /DEBUG/KDEBUG
      COMMON /MIXC/PALFA,PBETA,PNY,PY /CVFIX/VFIX
      COMMON /CG/GRAV,KONSG /CTEFF/TEFF,FLUX
      COMMON /NATURE/BOLTZK,CLIGHT,ECHARG,HPLNCK,PI,PI4C,RYDBRG,
     * STEFAN
      COMMON /CPF/PF,PFE,PFD,FIXROS,ITSTOP
      LOGICAL PF,PFE,PFD,FIXROS,ITSTOP
CUGJ FFR in excess     COMMON /CSPHER/NCORE,DIFLOG,RADIUS,RR(NDP),FFR(NDP)
      dimension ffr(ndp)
      COMMON /CSPHER/DIFLOG,RADIUS,RR(NDP),NCORE
C OWN COMMONS
      COMMON /CTRAN/X,S,BPLAN,XJ,HFLUX,XK,FJ(NDP),SOURCE(NDP),TAUS(NDP)
     & ,DTAUS(NDP),JTAU0,JTAU1,ISCAT
      COMMON /CSURF/HSURF,Y1(NRAYS)
      COMMON /ROSSC/ROSS(NDP),CROSS(NDP) /RHOC/RHO(NDP)
      COMMON /CARC1/ISTRAL,IDRAB1,IDRAB2,IDRAB3,IDRAB4,IDRAB5,IDRAB6,
     &              IARCH
      COMMON /CARC2/T,FC,FLUXME(NWL),TAU5(NDP),INORD
      COMMON /Cspec/spec(nwl,3),ispec
      COMMON /CI8/PGC,ROC,EC
      COMMON /NEWMO/NEWMOD
      COMMON /MASSE/RELM
      COMMON /CORRECT/KORT,KPP,TCONV
      COMMON /CIT/IT,ITMAX
C
C SPACE ALLOCATION
      COMMON /SPACE1/XJ1,XJ2,XJ3,TJ1,TJ2,XJT1,XJT2,XJT3,PRJ,TTT,PRT
     &  ,XJPE1,XJPE2,XJPE3,TPE,PRPE
      COMMON /SPACE2/ SPACEDUM1(NDP*7+NDP*NRAYS*5+NRAYS*2),
     &       SPACEDUM2(NDP*2),PFEAU(NRAYS,NDP),XMU(NRAYS,NDP),
     &       MMU(NDP),KSPACE2DUM(NRAYS+1)
C WAS:
C      COMMON /SPACE2/FCT,FCPE,PET,PEPE,
C    &   SPACEDUM(NDP*10+NDP*NRAYS*7+NRAYS*3+1-NDP*NDP*4)  
C --->WRONG, IT'S + OR - THAT.
C      DIMENSION XMU(NRAYS,NDP),PFEAU(NRAYS,NDP),MMU(NDP),
C     &          SP2(NDP*10+NDP*NRAYS*7+NRAYS*3+1)
C      EQUIVALENCE (FCT,SP2), (MMU,SP2(NDP*7+NDP*NRAYS*5+NRAYS*3+2)),
C     &            (PFEAU,SP2(NDP*10+NDP*NRAYS*5+NRAYS*3+2)),
C     &            (XMU(1,1),SP2(NDP*10+NDP*NRAYS*6+NRAYS*3+2))
C WAS (WITH NRAYS=21 AND NDP=40):
C SP2(5505)
C      EQUIVALENCE (FCT,SP2),(MMU,SP2(4545)),(PFEAU,SP2(4665)),
C     & (XMU(1,1),SP2(5505))
C
C WAS:   XMU(1)    DT(1)   SCRATC(1)     DPE(1)   TPE(1) @@@
      EQUIVALENCE (DT(1,1),SCRATC(1,1))
      EQUIVALENCE (XT(1),D(1)),(ST(1),V(1)),(DLNX(1),FC(1))
      DATA IVERS,IEDIT,ITER/21,3,0/
      DATA NFIRST/0/
C
C IN THIS SECTION THE MEAN INTENSITY IS ELIMINATED IN THE TRANSPORT EQUATIONS
C LEAVING THE EXPLICIT TEMPERATURE DEPENDANCE OF FLUX AND RADIATION PRESSURE
C IN THE MATRICES TTT AND PRT.
C
      ITER=ITER+1
      INORD=IEDIT+10*IVERS
      FNORD=.1*INORD
C
      IF (NEWMOD.EQ.9) THEN
        WRITE(66) RADIUS,TEFF,GRAV,RELM
        WRITE(66) NLAM,(XL(NLL),W(NLL),NLL=1,NLAM)
      ENDIF
C
C ZEROSET
      DO 110 I=1,NTAU
      PPRG(I)=0.
      RT(I)=0.
      RPR(I)=0.
      FFR(I)=0.
      ROSSP(I)=0.
      SUMW(I)=0.
      DO 110 J=1,NTAU
        TTT(I,J)=0.
        TPE(I,J)=0.
        PRT(I,J)=0.
        PRPE(I,J)=0.
110   CONTINUE
      kdtpe = 0
C
C CALCULATE DETAILED ROSSELAND MEAN
      REWIND 11
      KL=1
      DUMMY=ROSSOP(TT(1),PPE(1))
      PGA=PGC
      DO 116 K=1,NTAU
        KL=K
	DUMMY=ROSSOP(TT(K),PPE(K))
	ROSS(K)=1.0
        RHO(K)=ROC
	SUMW(K)=0.
	ROSSP(K)=0.
116   CONTINUE
C      DO 117 J=1,NLAM
      DO 117 J=1,NWTOT
        CALL OPAC(J,X,S)
        WRITE(11) X,S
C        Y=((XL(J)/1.E4)**2)**3
        Y=((WLOS(J)/1.E4)**2)**3
        DO 117 K=1,NTAU
C          YA=EXP(-1.438E8/(TT(K)*XL(J)))
          YA=EXP(-1.438E8/(TT(K)*WLOS(J)))
          YA=YA/(1.-YA)**2/Y
C          SUMW(K)=SUMW(K)+W(J)*YA
C          ROSSP(K)=ROSSP(K)+W(J)*YA/(ROSS(K)*(X(K)+S(K)))
          SUMW(K)=SUMW(K)+WLSTEP(J)*YA
        if (wlos(j).le.5000. .or. wlos(j).ge.1.e5) go to 117
          ROSSP(K)=ROSSP(K)+WLSTEP(J)*YA/(ROSS(K)*(X(K)+S(K)))
117   CONTINUE
      REWIND 11
C
C TEMPERATURE AND ELECTRON PRESSURE PERTURBATIONS.
C KEEP THEM SMALL, TO STAY ON THE LINEAR PART.
      DTX=0.001
      DPEX=0.001
      DO 111 K=1,NTAU
        T(K)=TT(K)*DTX
        PE(K)=PPE(K)*DPEX
        TT(K)=TT(K)+T(K)
	ROSSP(K)=SUMW(K)/ROSSP(K)
	SUMW(K)=0.
	ROSST(K)=0.
111   CONTINUE
      kdtpe = 1
C
C FIRST WAVELENGTH LOOP, TO CALCULATE XT,ST AND SAVE.
      REWIND 12
C      DO 112 J=1,NLAM
      DO 112 J=1,NWTOT
        CALL OPAC(J,X,S)
        WRITE(12) X,S
C        Y=((XL(J)/1.E4)**2)**3
        Y=((WLOS(J)/1.E4)**2)**3
        DO 112 K=1,NTAU
C          YA=EXP(-1.438E8/(TT(K)*XL(J)))
          YA=EXP(-1.438E8/(TT(K)*WLOS(J)))
          YA=YA/(1.-YA)**2/Y
C          SUMW(K)=SUMW(K)+W(J)*YA
C          ROSST(K)=ROSST(K)+W(J)*YA/(ROSS(K)*(X(K)+S(K)))
          SUMW(K)=SUMW(K)+WLSTEP(J)*YA
        if (wlos(j).le.5000. .or. wlos(j).ge.1.e5) go to 112
          ROSST(K)=ROSST(K)+WLSTEP(J)*YA/(ROSS(K)*(X(K)+S(K)))
112   CONTINUE
      DO 113 K=1,NTAU
	ROSST(K)=SUMW(K)/ROSST(K)
	SUMW(K)=0.
	ROSSPE(K)=0.
        TT(K)=TT(K)-T(K)
        PPE(K)=PPE(K)+PE(K)
113   CONTINUE
        kdtpe = 2            !information to tstgem about computing dpg/dpe
      REWIND 12
      CALL CLOCK
C
C SECOND WAVELENGTH LOOP, TO CALCULATE XPE,SPE AND SAVE.
      REWIND 14
      KL=1
      DUMMY=ROSSOP(TT(1),PPE(1))
      PGPE=PGC
C      DO 114 J=1,NLAM
      DO 114 J=1,NWTOT
        CALL OPAC(J,X,S)
        WRITE(14) X,S
C        Y=((XL(J)/1.E4)**2)**3
        Y=((WLOS(J)/1.E4)**2)**3
        DO 114 K=1,NTAU
C          YA=EXP(-1.438E8/(TT(K)*XL(J)))
          YA=EXP(-1.438E8/(TT(K)*WLOS(J)))
          YA=YA/(1.-YA)**2/Y
C          SUMW(K)=SUMW(K)+W(J)*YA
C          ROSSPE(K)=ROSSPE(K)+W(J)*YA/(ROSS(K)*(X(K)+S(K)))
          SUMW(K)=SUMW(K)+WLSTEP(J)*YA
        if (wlos(j).le.5000. .or. wlos(j).ge.1.e5) go to 114
          ROSSPE(K)=ROSSPE(K)+WLSTEP(J)*YA/(ROSS(K)*(X(K)+S(K)))
114   CONTINUE
      kdtpe = 3
      REWIND 14
      CALL CLOCK
C
C FROM THIS POINT ON, ROSS() HOLDS THE TRUE ROSSELAND MEAN.  CROSS HOLDS
C THE RATIO OF THE TRUE TO APPROXIMATE MEANS, WHICH ARE NEEDED IN TRYCK.
      IF(NPRINT.GE.2) write (6,*) 'rosseland values'
      DO 123 K=1,NTAU
        KL=K
	ROSSPE(K)=SUMW(K)/ROSSPE(K)
        PPE(K)=PPE(K)-PE(K)
	CROSS(K)=ROSSP(K)/ROSSOP(TT(K),PPE(K))
        ROSS(K)=ROSSP(K)
        PTAU(K)=GRAV*TAU(K)/ROSS(K)
	IF(NPRINT.GE.2) write (6,'(1x,i3,4(1pe12.3))')
     &     k,ross(k),cross(k),rosst(k),rosspe(k)
123   CONTINUE
      CALL CLOCK
C
C RIGHT HAND SIDE IN PRESSURE EQUATION
      KL=1
      CALL TAET(TT(1),PPE(1),PG,RO,DUM)
      DLNP=1./(1.+(ROSSPE(1)-ROSS(1))/ROSS(1)*PGA/(PGPE-PGA))
      GRVR=GRAV*(RADIUS/RR(1))**2
      IF (KONSG.EQ.1) GRVR=GRAV   !test for effect of varying g
      RP(1)=GRVR*TAU(1)/(ROSS(1)*DLNP)+PPR(1)-PP(1)
C SIMPSONS RULE
      DO 101 K=2,NTAU
      F0=PTAU(K-1)
      F1=FOUR(PTAU,TAULN,K,NTAU)
      F2=PTAU(K)
      RP(K)=(F0+4.*F1+F2)*DTAULN(K)/6.-(PP(K)-PP(K-1))
101   CONTINUE
C
C CALCULATE RADII
      RR(1)=0.
      DO 102 K=2,NTAU
      IF (TAU(K).LT.0.67) K0=K
102   RR(K)=RR(K-1)-0.5*DTAULN(K)
     & *(TAU(K)/(ROSS(K)*RHO(K))+TAU(K-1)/(ROSS(K-1)*RHO(K-1)))
C...      MIHAL=K0
C...      WRITE(7,104) MIHAL
C...104   FORMAT('0MIHAL HAS BEEN SET EQUAL TO',I3)
      Y=RADIUS-RR(K0)
      DO 103 K=1,NTAU
      RR(K)=RR(K)+Y
      IF (KONSG.EQ.1) RR(K)=RADIUS      !study effect of const R and g
103   PTAU(K)=PTAU(K)*(RADIUS/RR(K))**2
C
      WRITE(6,*) 'NTAU IN SOLVE_sph = ',NTAU
      WRITE(6,*) 'K,ROSS(K),RHO(K),DTAULN(K),RR(K)'
      do 104 k=1,ntau
      WRITE(6,*) K,ROSS(K),RHO(K),DTAULN(K),RR(K)
104   CONTINUE
C
C TIME
      CALL CLOCK
      MSA=0
C
C WAVELENGTH LOOP
      IF(NPRINT.GE.2 .AND. PF) 
     &     WRITE(6,48) FNORD,ITER,IDRAB1,IDRAB2,IDRAB3,IDRAB4,
     &                   IDRAB5,IDRAB6
      IF(NPRINT.GE.2 .AND. PF) WRITE(6,59)
      FTOT=0.
C      DO 150 J=1,NLAM
      DO 150 J=1,NWTOT
      DO 130 K=1,NTAU
C      BPLAN(K)=BPL(TT(K),XL(J))
C130   DBPL(K)=DIVBP(TT(K),XL(J))
      BPLAN(K)=BPL(TT(K),WLOS(J))
130   DBPL(K)=DIVBP(TT(K),WLOS(J))
      IF (J.EQ.1 .OR. J.EQ.NWTOT) THEN
         WRITE(6,*) ' J,WLOS(J),BPLAN(K): ',J,WLOS(J)
         WRITE(6,*) (BPLAN(K),K=1,40)
      END IF
C
C CALCULATE OPACITY DERIVATIVES AT CONSTANT GAS PRESSURE
      READ (11) X,S
      READ (12) XT,ST
      READ (14) XPE,SPE
      DO 131 K=1,NTAU
        X(K)=X(K)/ROSS(K)
        S(K)=S(K)/ROSS(K)
        XT(K)=XT(K)/ROSST(K)
        ST(K)=ST(K)/ROSST(K)
        XPE(K)=XPE(K)/ROSSPE(K)
        SPE(K)=SPE(K)/ROSSPE(K)
C	if (j.eq.100.and.k.eq.1) write (7,'(1x,7(1pe10.2))')
C     &   x(k),t(k),xt(k),log(xt(k)/x(k)),log(xt(k)/x(k))*tt(k)/t(k)
        XLOG(K)=log10(X(K))
        XT(K)=log(XT(K)/X(K))*X(K)/T(K)
        ST(K)=log(ST(K)/S(K))*S(K)/T(K)
        XPE(K)=log(XPE(K)/X(K))*X(K)/PE(K)
        SPE(K)=log(SPE(K)/S(K))*S(K)/PE(K)
        DLNX(K)=XT(K)*(TT(K)/2.3)/X(K)
C	if (j.eq.100) write (7,'(1x,7(1pe10.2))')
C     &    x(k),xt(k)/x(k)*tt(k),xpe(k)/x(k)*ppe(k)
131   CONTINUE
C TIME
      MS=MSA
      MSA=0
      MSOPAC=MS-MSA
C
C SOLVE TRANSPORTEQUATION WITH OLD STRATIFICATION.
      CALL TRANEQ_sph
C
C ??
      DO 132 K=1,NTAU
      IF(XT(K).LT.0.0.AND.XT(K)*(XJ(K)-BPLAN(K)).GT.X(K)*DBPL(K))
     & XT(K)=X(K)*DBPL(K)/(XJ(K)-BPLAN(K))
132   CONTINUE
      MS=MSA
      MSA=0
      MSTRAN=MS-MSA
C
C FLUX TO PRINT
C      HFLUX1=4.*PI*HSURF*(RR(1)/RADIUS)**2
      HFLUX1=4.*PI*HSURF
      HFLUX1=MAX(1.0D-99,HFLUX1)
      HFLUX2=4.*PI*HFLUX(NTAU)*(RR(NTAU)/RADIUS)**2
      FLUXME(J)=HFLUX1/PI
C      GFLUX1=HFLUX1/PI*XL(J)**2/CLIGHT
C      GFLUX2=HFLUX2/PI*XL(J)**2/CLIGHT
      GFLUX1=HFLUX1/PI*WLOS(J)**2/CLIGHT
      GFLUX2=HFLUX2/PI*WLOS(J)**2/CLIGHT
      FFLUX1=-2.5*log10(MAX(GFLUX1,1.0D-99))
      FFLUX2=-2.5*log10(MAX(GFLUX2,1.0D-99))
      spec(j,1) = wlos(j)
      spec(j,2) = hsurf
      spec(j,3) = fluxme(j)
C
C OUTPUT IN CASE OF NEWMOD=9 : STORE LIMB FUNCTION
      IF (NEWMOD.EQ.9) THEN
        WRITE(66) MMU(1),(XMU(NLX,1),Y1(NLX),NLX=1,MMU(1))
        WRITE(66) RR(1)
        GOTO 150
      ENDIF
C
C SUM UP RADIATIVE FLUXES
      DO 133 K=1,NTAU
C      FFR(K)=FFR(K)+W(J)*HFLUX(K)
      FFR(K)=FFR(K)+WLSTEP(J)*HFLUX(K)
      if(NPRINT.GE.2 .and. hflux(k).le.0.0) then
C         write (7,*) 'J,K,W(J),HFLUX(K),FFR(K),X(K),S(K)',
C     &    J,K,W(J),HFLUX(K),FFR(K),X(K),S(K)
         write (7,*) 'J,K,WLSTEP(J),HFLUX(K),FFR(K),X(K),S(K)',
     &    J,K,WLSTEP(J),HFLUX(K),FFR(K),X(K),S(K)
      endif 
133   CONTINUE
C
C UPPER BOUNDARY
      DO 140 K=1,NTAU
      IF (K.GT.1) GO TO 143
      PB=XJ(1)/FJ(1)
      PC=XJ(2)/FJ(2)
      IF (TAUS(1).LT.0.1) GO TO 141
      EX=EXP(-TAUS(1))
      EX1=1.-EX
      GO TO 142
141   EX1=TAUS(1)*(1.-0.5*TAUS(1)*(1.-0.333333*TAUS(1)))
      EX=1.-EX1
142   YA=DTAUS(2)*(EX1+0.5*DTAUS(2))/(X(1)+S(1))
      YB=((1.+DTAUS(2))*(PB-SOURCE(1))+SOURCE(1)*EX)*DTAUS(2)
     & /(X(1)+S(1)+X(2)+S(2))
      XJ2(1)=DTAUS(2)+0.5*DTAUS(2)**2-YA*S(1)*FJ(1)
      XJ3(1)=-1.
      XJT2(1)=-YA*X(1)*DBPL(1)+YB*(XT(1)+ST(1))
     & -YA*(BPLAN(1)-XJ(1))/(X(1)+S(1))*(S(1)*XT(1)-X(1)*ST(1))
     & -DTAUS(2)*SOURCE(1)*EX*TAUS(1)/(X(1)+S(1))*(XT(1)+ST(1))
      XJT3(1)=YB*(XT(2)+ST(2))
      XJPE2(1)=YB*(XPE(1)+SPE(1))
     & -YA*(BPLAN(1)-XJ(1))/(X(1)+S(1))*(S(1)*XPE(1)-X(1)*SPE(1))
     & -DTAUS(2)*SOURCE(1)*EX*TAUS(1)/(X(1)+S(1))*(XPE(1)+SPE(1))
      XJPE3(1)=YB*(XPE(2)+SPE(2))
      GO TO 170
C
C INTERNAL POINTS
143   PA=PB
      PB=PC
      IF (K.EQ.JTAU) GO TO 144
      PC=XJ(K+1)/FJ(K+1)
      DTAUSK=0.5*(DTAUS(K+1)+DTAUS(K))
      XJ1(K)=1./DTAUS(K)
      XJ3(K)=1./DTAUS(K+1)
      XJ2(K)=-DTAUSK*(1.-FJ(K)*S(K)/(X(K)+S(K)))
      YA=((PB-PA)/DTAUS(K)-(PB-SOURCE(K))*0.5*DTAUS(K))
     & /(X(K)+S(K)+X(K-1)+S(K-1))
      YB=(-(PC-PB)/DTAUS(K+1)-(PB-SOURCE(K))*0.5*DTAUS(K+1))
     & /(X(K)+S(K)+X(K+1)+S(K+1))
      XJT1(K)=YA*(XT(K-1)+ST(K-1))
      XJT3(K)=YB*(XT(K+1)+ST(K+1))
      XJT2(K)=DTAUSK*X(K)/(X(K)+S(K))*DBPL(K)+(YA+YB)*(XT(K)+ST(K))
     & +DTAUSK*(BPLAN(K)-XJ(K))/(X(K)+S(K))**2*(S(K)*XT(K)-X(K)*ST(K))
      XJPE1(K)=YA*(XPE(K-1)+SPE(K-1))
      XJPE3(K)=YB*(XPE(K+1)+SPE(K+1))
      XJPE2(K)=(YA+YB)*(XPE(K)+SPE(K))
     & +DTAUSK*(BPLAN(K)-XJ(K))/(X(K)+S(K))**2*(S(K)*XPE(K)-X(K)*SPE(K))
      GO TO 170
C
C OPTICALLY THICK POINTS
144   XJ1(K)=0.
      XJ2(K)=-1.
      XJ3(K)=0.
      XJT1(K)=0.
      XJT2(K)=DBPL(K)
      XJT3(K)=0.
      XJPE1(K)=0.
      XJPE2(K)=0.
      XJPE3(K)=0.
C
C TEMPERATURE EQUATION
170   IF (K.GT.MIHAL) GO TO 171
C      Y=W(J)*X(K)
      Y=WLSTEP(J)*X(K)
      IF (K.GT.2) Y=Y*((TAU(K+1)-TAU(K))*(X(K)+S(K)+X(K+1)+S(K+1))
     & +(TAU(K)-TAU(K-1))*(X(K)+S(K)+X(K-1)+S(K-1)))/(X(K)+S(K))
      RT(K)=RT(K)+Y*(XJ(K)-BPLAN(K))
      TJ2(K)=-Y*FJ(K)
      TJ1(K)=0.
c ??
***      TTT(K,K)=TTT(K,K)+AMAX1(0.,Y*DBPL(K)+Y*(BPLAN(K)-XJ(K))
***     & *XT(K)/X(K))
      TTT(K,K)=TTT(K,K)+Y*DBPL(K)+Y*(BPLAN(K)-XJ(K))*XT(K)/X(K)
      TPE(K,K)=TPE(K,K)+Y*(BPLAN(K)-XJ(K))*XPE(K)/X(K)
      GO TO 172
C171   Y=4.*W(J)
171   Y=4.*WLSTEP(J)
      RT(K)=RT(K)-Y*HFLUX(K)
      XHK=0.557*(PB-PA)/DTAUS(K)
      FH=HFLUX(K)/XHK
C
C DEBUG
      IF (FH.GT.0.7.AND.FH.LT.1.4) GO TO 174
      IF(NPRINT.GE.2)
     &WRITE(6,173) JTAU0,JTAU1,K,FH,HFLUX(K),XHK,DTAUS(K),BPLAN(K-1)
     & ,BPLAN(K),XJ(K-1),XJ(K),XK(K-1),XK(K),PA,PB
173   FORMAT('0J0,J1,K,FH,HF,XH,DTAU,B,J,K,P=',3I3,1P,4E10.3/10E10.3)
      IF (K.GT.JTAU1) GO TO 174
      NMU=MMU(K-1)
      IF(NPRINT.GE.2) WRITE(6,175) (XMU(I,K-1),I=1,NMU)
      IF(NPRINT.GE.2) WRITE(6,175) (PFEAU(I,K-1),I=1,NMU)
      NMU=MMU(K)
      IF(NPRINT.GE.2) WRITE(6,175) (XMU(I,K),I=1,NMU)
      IF(NPRINT.GE.2) WRITE(6,175) (PFEAU(I,K),I=1,NMU)
175   FORMAT(10E11.4)
174   CONTINUE
C
      TJ1(K)=-Y*FH*0.557/DTAUS(K)
      TJ2(K)=Y*FH*0.557/DTAUS(K)
      TTT(K,K-1)=TTT(K,K-1)-Y*HFLUX(K)*(XT(K-1)+ST(K-1))
     & /(X(K)+S(K)+X(K-1)+S(K-1))
      TTT(K,K)=TTT(K,K)-Y*HFLUX(K)*(XT(K)+ST(K))
     & /(X(K)+S(K)+X(K-1)+S(K-1))
      TPE(K,K-1)=TPE(K,K-1)-Y*HFLUX(K)*(XPE(K-1)+SPE(K-1))
     & /(X(K)+S(K)+X(K-1)+S(K-1))
      TPE(K,K)=TPE(K,K)-Y*HFLUX(K)*(XPE(K)+SPE(K))
     & /(X(K)+S(K)+X(K-1)+S(K-1))
172   CONTINUE
C
C EQUATION OF RADIATIVE PRESSURE
C      Y=PI4C*W(J)
      Y=PI4C*WLSTEP(J)
      RPR(K)=RPR(K)+Y*XK(K)
      PRJ(K)=-Y*XK(K)*FJ(K)/XJ(K)
      PRJ2(K)=PRJ(K)
CUGJ 900503: The gradient of the radiative pressure
      PPRG(K) = PI4C*WLSTEP(J)*HFLUX(K)*(X(K)+S(K))*ROSS(K) + PPRG(K)
C
C END OF TAU LOOP
140   CONTINUE
C
C ELIMINATE THIS WAVELENGTH
      IF(NEW.EQ.1) CALL ALGEBN(NTAU)
C
C TIME
      MS=MSA
      MSA=0
      MS=MS-MSA
C
C END OF WAVELENGTH LOOP
C      HW1=HFLUX1*W(J)
      HW1=HFLUX1*WLSTEP(J)
      FTOT=FTOT+HW1
      HW2=HFLUX2*WLSTEP(J)
C      HW2=HFLUX2*W(J)
C      WAVEN=1.E4/XL(J)
C      TRAD1=1.438E8/(XL(J)*log(1.+1.191E7*PI/HFLUX1*WAVEN**5))
      WAVEN=1.E4/WLOS(J)
      TRAD1=1.438E8/(WLOS(J)*
     #  log(1.0D+0+1.191D7*PI/HFLUX1*WAVEN**5))
      X01=log10(X(01))
      X25=log10(X(25))
      S01=log10(S(01))
      S25=log10(S(25))
C      IF(PF) WRITE(7,58) J,XL(J),W(J),HFLUX1,HW1,WAVEN,GFLUX1,FFLUX1
      IF(PF.and.NPRINT.GE.2) 
     * WRITE(6,58) J,WLOS(J),WLSTEP(J),HFLUX1,HW1,WAVEN,GFLUX1,FFLUX1
     * ,TRAD1,X01,X25,S01,S25,HW2,DLNX(01),DLNX(25)
     * ,JTAU0,JTAU1,XMU(ISCAT,JTAU0)
C...      IF(PFD) WRITE(7,30) (XLOG(K),K=1,26)
C...      IF(PFD) WRITE(7,30) (DLNX(K),K=1,26)
C...30    FORMAT(1X,26F5.2)
150   CONTINUE
C
C OUTPUT IN CASE NEWMOD=9 : STORE FLUXES
      IF (NEWMOD.EQ.9) THEN
        WRITE(66) (FLUXME(NLL),NLL=1,NWTOT)
        GOTO 900
      ENDIF
      TEFFP=TEFF*(FTOT/(FLUX*PI))**.25
      Y=FTOT/PI
      IF(NPRINT.GE.2) WRITE(6,65) FTOT,Y,TEFFP
      DO 154 K=1,NTAU
      Y=0.
      DO 155 L=1,NTAU
155   Y=Y+TTT(K,L)
      IF(PFD.AND.PFE .and. NPRINT.GE.2 ) 
     * WRITE(6,151) K,Y,TTT(K,K),(TTT(K,L),L=1,NTAU)
151   FORMAT(' TTT='/(I3,2(1PE10.3),20(/8(1PE10.3))))
c ??
      IF (K.LE.MIHAL.AND.Y.LT.0.0) TTT(K,K)=TTT(K,K)-Y
154   CONTINUE
      DO 153 K=1,NTAU
153   FFR(K)=FFR(K)*4./FLUX*(RR(K)/RADIUS)**2
      IF(NPRINT.GE.2) WRITE(6,152) (FFR(K),K=1,NTAU)
152   FORMAT('0FFR='/(1P,8E10.3))
C
C TIME
      CALL CLOCK
C
C PRINT PRESSURE EQUATION
      IF(PF .AND. NPRINT.GE.2) 
     &    WRITE(6,48) FNORD,ITER,IDRAB1,IDRAB2,IDRAB3,IDRAB4,
     &                   IDRAB5,IDRAB6
      IF(PF .AND. NPRINT.GE.2) WRITE(6,62)
      IF(PF .AND. NPRINT.GE.2) WRITE(6,63)
      DO 161 K=1,NTAU
      KL=K
      ROSST(K)=(ROSST(K)-ROSS(K))/T(K)
      ROSSPE(K)=(ROSSPE(K)-ROSS(K))/PE(K)
C
C TAU SCALES
      IF(K.GT.1) GO TO 162
      CALL KAP5(TT(1),PPE(1),ABSK)
      TAU5(1)=TAU(1)*ABSK/ROSS(1)
      TAUP=TAU(1)*CROSS(1)
      YC=ABSK/ROSS(1)
      YD=CROSS(1)
      GO TO 163
162   CONTINUE
      CALL KAP5(TT(K),PPE(K),ABSK)
      YA=ABSK/ROSS(K)
      YB=CROSS(K)
      TAU5(K)=TAU5(K-1)+.5*(YA+YC)*(TAU(K)-TAU(K-1))
      TAUP=TAUP+.5*(YB+YD)*(TAU(K)-TAU(K-1))
      YC=YA
      YD=YB
163   CONTINUE
      ROSPE=ROSSPE(K)*PPE(K)/ROSS(K)
      ROST=ROSST(K)*TT(K)/ROSS(K)
      IF(PF.and.NPRINT.GE.2 ) 
     *  WRITE(6,51) K,TAU(K),PTAU(K),ROSS(K),ROSPE,ROST
     * ,YC,YD,TAU5(K),TAUP,RR(K),K
161   CONTINUE
      CALL CLOCK
900   CONTINUE
      RETURN
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C+DEF,D=MATRIX
C
C IN THIS SECTION WE COMPUTE MATRIX ELEMENTS FOR THE REST OF THE PROBLEM
C
C
      ENTRY MATRIX_sph
C TIME
      MSA=0
C
C ZEROSET
      IF(PF .AND. NPRINT.GE.2) 
     &    WRITE(6,48) FNORD,ITER,IDRAB1,IDRAB2,IDRAB3,IDRAB4,
     &                   IDRAB5,IDRAB6
      IF(PF .AND. NPRINT.GE.2) WRITE(6,57)
      IF(PF .AND. NPRINT.GE.2) WRITE(6,55)
      GRAD=0.
      DO 201 I=1,NTAU
      DO 201 J=1,NTAU
      VT(I,J)=0.
      VPE(I,J)=0.
      FCT(I,J)=0.
      FCPE(I,J)=0.
      DT(I,J)=0.
      DPE(I,J)=0.
      PET(I,J)=0.
      PEPE(I,J)=0.
201   CONTINUE
C
C VFIX OPTION
      IF(VFIX.EQ.0.) GOTO 230
      DO 231 KK=2,NTAU
      K=(1+NTAU)-KK
      IF(K.LE.NOCONV) GOTO 230
231   VV(K)=MAX(VV(K+1)*EXP(-DTAULN(K+1)/VFIX),VV(K))
230   CONTINUE
C
C TAU LOOP
      DO 200 K=1,NTAU
      KL=K
C
C TERMODYNAMICAL QUANTITIES WITH PARTIAL DERIVATIVES
      K1=MAX0(1,K-1)
      TMEAN=.5*(TT(K)+TT(K1))
      PEMEAN=.5*(PPE(K)+PPE(K1))
      PMEAN=.5*(PP(K)+PP(K1))
      PRMEAN=.5*(PPR(K)+PPR(K1))
      ROSSMN=.5*(ROSS(K)+ROSS(K1))
      IF(K.GT.NOCONV) GO TO 213
C NO CONVECTION
      CALL TERMON(TT(K),PPE(K),PPR(K),PG,PGT,PGPE,RO,ROT,ROPE,CP,ADIA,Q)
      RO1=RO
      PG1=PG
      PG1T=0.
      PG1PE=0.
      CPT=0.
      CPPE=0.
      ADIAT=0.
      ADIAPE=0.
      QT=0.
      QPE=0.
      GO TO 212
C CONVECTION
213   CONTINUE
      DERET=.01
      DEREP=.15
      TDELT=TMEAN*DERET
      PEDELT=PEMEAN*DEREP
      T1=TMEAN+TDELT
      PE1=PEMEAN+PEDELT
      CALL TERMON(T1,PEMEAN,PRMEAN,Y,YA,YB,YC,YD,YE,CP1,ADIA1,Q1)
      CALL TERMON(TMEAN,PE1,PRMEAN,Y,YA,YB,YC,YD,YE,CP2,ADIA2,Q2)
      CALL TERMON(TMEAN,PEMEAN,PRMEAN,PG1,PG1T,PG1PE,RO,ROT,ROPE,CP,
     &ADIA,Q)
      CALL TERMON(TT(K),PPE(K),PPR(K),PG,PGT,PGPE,RO1,Y,YA,YB,YC,YD)
      CPPE=(CP2-CP)/PEDELT
      ADIAPE=(ADIA2-ADIA)/PEDELT
      QPE=(Q2-Q)/PEDELT
      CPT=(CP1-CP)/TDELT
      ADIAT=(ADIA1-ADIA)/TDELT
      QT=(Q1-Q)/TDELT
212   CONTINUE
      RHO(K)=RO1
C
C DEPTH SCALE
      IF(K.EQ.1) ZZ(K)=0.
      IF(K.GT.1) ZZ(K)=ZZ(K-1)-.5*DTAULN(K)*
     & (TAU(K)/ROSS(K)/RHO(K)+TAU(K-1)/ROSS(K-1)/RHO(K-1))
      IF(TAULN(K).LT.0.0) KK0=K
C
C RADIATION PRESSURE
      RPR(K)=RPR(K)-PPR(K)
CC    PRJ IS NOW FREE
CC    PRT IS ALREADY INITIATED
CC    PRPR IS UNITY
C
C TOTAL PRESSURE
      GRVR=GRAV*(RADIUS/RR(K))**2
      GEFF(K)=GRVR
      IF (KONSG.EQ.1) GRVR=GRAV   !test for effect of varying g
      IF(K.GT.1) GO TO 202
      Y=GRVR*TAU(1)/(ROSS(1)*DLNP)
      PPTT(1)=Y*ROSST(1)/ROSS(1)
      PPPE(1)=Y*ROSSPE(1)/ROSS(1)
      GO TO 203
202   Y=GRVR*DTAULN(K)*.5
      YY=Y*TAU(K)/ROSS(K)**2
      Y=Y*TAU(K-1)/ROSS(K-1)**2
      PPPE(K)=YY*ROSSPE(K)
      PPTT(K)=YY*ROSST(K)
      PPPE(K+NTAU-1)=Y*ROSSPE(K-1)
      PPTT(K+NTAU-1)=Y*ROSST(K-1)
203   CONTINUE
C
C CONVECTION EFFICIENCY GAMMA
      HSCALE=(PG1+PRMEAN)/GRVR/RO
      OMEGA=PALFA*HSCALE*RO*ROSSMN
      IF(PALFA.EQ.0.) OMEGA=HSCALE*RO*ROSSMN
      Y=PY*OMEGA**2
      YY=(Y-1.)/(Y+1.)
      THETA=OMEGA/(1.+Y)
      GAMMA=-CP*RO/(8.*STEFAN*TMEAN**3*THETA)
CC    GGG IS UNITY
      GV(K)=GAMMA
      IF(PBETA.GT.0.) VV(K)=MIN(VV(K),SQRT(0.5*PP(K)/PBETA/RO))
      GG(K)=-GAMMA*VV(K)
      ROSPEM=.5*(ROSSPE(K)+ROSSPE(K1))
      ROSSTM=.5*(ROSST(K)+ROSST(K1))
      GPE(K)=-GG(K)*(CPPE/CP+ROPE/RO+YY*(ROSPEM/ROSSMN+PG1PE/PG1))
      GT(K)=-GG(K)*(CPT/CP+ROT/RO-3./TMEAN+YY*(ROSSTM/ROSSMN+PG1T/PG1))
CC    RG IS ZERO
C
C GRADIENT DIFFERENCE
      IF(K.LE.NOCONV) GO TO 206
      DELP=PP(K)-PP(K-1)
      DELT=TT(K)-TT(K-1)
      PM=PP(K)+PP(K-1)
      TM=TT(K)+TT(K-1)
      Y=1.+GG(K)
      YY=-GG(K)/Y
      GRAD=log(TT(K)/TT(K-1))/log(PP(K)/PP(K-1))
      NEWV=DD(K).GT.0..AND.VV(K).EQ.0..AND.PALFA.GT.0..AND.K.GT.2
      IF(.NOT.NEWV) GO TO 263
      VV(K)=SQRT(GRVR*HSCALE*Q*PALFA**2*DD(K)/PNY)
      IF(PF .AND. NPRINT.GE.2) WRITE(6,64) K,VV(K)
      GO TO 203
263   CONTINUE
      NEWV=GRAD.GE.ADIA.AND.VV(K).EQ.0..AND.PALFA.GT.0..AND.K.GT.2
      IF(.NOT.NEWV) GO TO 204
      VV(K)=VVMLT(GRAD-ADIA,GRVR*HSCALE*Q*PALFA**2/PNY,GAMMA**2)
      IF(PF .AND. NPRINT.GE.2) WRITE(6,64) K,VV(K)
      GO TO 203
204   CONTINUE
      YYY=GRAD-ADIA
C DDD IS UNITY
C ******* NEXT STATEMENT FIXES T80G4M0 BUT NOT TESTED FOR ALL MODELS
      IF(DD(K).EQ.0..AND.VV(K).GT.0.) DD(K)=-YY*YYY
      RD(K)=-YY*YYY-DD(K)
      DG(K)=-YYY/Y**2
      DT(K,K)=YY*(GRAD*(1./DELT-1./TM)-.5*ADIAT)
      DT(K,K-1)=YY*(GRAD*(-1./DELT-1./TM)-.5*ADIAT)
      DP(K)=YY*GRAD*(1./PM-1./DELP)
      DP(K+NTAU-1)=YY*GRAD*(1./PM+1./DELP)
      DPE(K,K)=-.5*YY*ADIAPE
      DPE(K,K-1)=-.5*YY*ADIAPE
      GO TO 205
206   DD(K)=0.
      RD(K)=0.
      DG(K)=0.
      DP(K)=0.
      DP(K+NTAU-1)=0.
205   CONTINUE
C
C VFIX OPTION
      IF(VFIX.EQ.0.) GOTO 280
      Y=0.
      IF(K.GT.NOCONV) Y=-VV(K)
      GOTO 207
280   CONTINUE
C
C CONVECTIVE VELOCITY
CC    VVV IS UNITY
      Y=0.
      IF(DD(K).LE.0.) GOTO 207
      Y=-SQRT(GRVR*HSCALE*Q*PALFA**2*DD(K)/PNY)
      VD(K)=Y*.5/DD(K)
      IF(-Y.GT.2.*VV(K)) VD(K)=VD(K)*2.
      VT(K,K)=.25*Y*(QT/Q+PGT/PG-ROT/RO)
      VT(K,K-1)=VT(K,K)
      VPE(K,K)=.25*Y*(QPE/Q+PGPE/PG-ROPE/RO)
      VPE(K,K-1)=VPE(K,K)
      GO TO 208
207   VD(K)=0.
208   CONTINUE
      RV(K)=-Y-VV(K)
C
C TURBULENT PRESSURE.
CC    PTPT IS UNITY
      Y=-PBETA*VV(K)**2
      PPT(K)=-RO*Y
      PTT(K)=Y*ROT
      PTPE(K)=Y*ROPE
      PTV(K)=-PBETA*2.*VV(K)*RO
C
C CONVECTIVE FLUX
      Y=-CP*RO*PALFA*TMEAN/2./PI
CC    FCFC IS UNITY
      YY=Y*VV(K)*DD(K)
      RFC(K)=-YY-FFC(K)
      FCD(K)=Y*VV(K)
      FCV(K)=Y*DD(K)
      IF(K.LE.NOCONV) GO TO 217
      FCT(K,K)=.5*YY*(CPT/CP+ROT/RO+1./TMEAN)
      FCT(K,K-1)=FCT(K,K)
      FCPE(K,K)=.5*YY*(CPPE/CP+ROPE/RO)
      FCPE(K,K-1)=FCPE(K,K)
217   CONTINUE
C
C ELECTRON PRESSURE
209   RPE(K)=PP(K)-PG-PPR(K)
      PET(K,K)=PGT
      PEPE(K,K)=PGPE
CC    PEPR=PEPT=1.  PEP=-1.
210   CONTINUE
C
C TEMPERATURE
CC    TTT IS ALREADY INITIATED
      IF(K.GT.MIHAL) GO TO 261
C STR"MGREN CONDITION
      RT(K)=RT(K)+(FFC(K+1)-FFC(K))
      GO TO 262
261   CONTINUE
C FLUXCONSTANCY
      RT(K)=RT(K)+FLUX*(RADIUS/RR(K))**2-FFC(K)
262   CONTINUE
C
C END OF TAU LOOP
      PGT=PGT*TT(K)/PG
      PGPE=PGPE*PPE(K)/PG
      IF(PF .AND. NPRINT.GE.2) 
     *   WRITE(6,51) K,TAU(K),HSCALE,ADIA,GRAD,CP,Q,PG,RO,PGPE
     * ,PGT,K
200   CONTINUE
C
C SUBTRACT CENTERED TURBULENT PRESSURE
      DO 216 K=2,NTAU
      K1=MIN0(K+1,NTAU)
216   RPE(K)=RPE(K)-.5*(PPT(K)+PPT(K1))
C
C SUBTRACT ZZ(TAU=1) FROM ZZ
      KK0 = MAX(1,KK0)
      ZZ0=ZZ(KK0)
      DO 283 K=1,NTAU
283   ZZ(K)=ZZ(K)-ZZ0
C
C TIME
      CALL CLOCK
C
C PRINT
      IF (NPRINT.GE.2) THEN
      IF(PF.OR.PFD) WRITE(6,48) FNORD,ITER,IDRAB1,IDRAB2,IDRAB3,IDRAB4,
     &                          IDRAB5,IDRAB6
      IF(PF.OR.PFD) WRITE(6,54)
      IF(PF.OR.PFD) WRITE(6,52)
      IF(PF.OR.PFD) 
     &  WRITE(6,51) (I,TAU(I),RPR(I),PPT(I),RP(I),GG(I),RD(I)
     & ,RV(I),RFC(I),RPE(I),RT(I),I,I=1,NTAU)
      END IF
C
C SAVE DD-MATRICES
      DO 270 K=1,NTAU
      DPS(K)=DP(K)
      DPES(K)=DPE(K,K)
      DTS(K)=DT(K,K)
      D(K)=RD(K)
      IF(K.EQ.1) GO TO 270
      DPS(K+NTAU-1)=DP(K+NTAU-1)
      DPES(K+NTAU-1)=DPE(K,K-1)
      DTS(K+NTAU-1)=DT(K,K-1)
270   CONTINUE
C
C+DEF,D=ELIMIN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C GAUSS ELIMINATION TO UPPER TRIANGULAR FORM.
C
C TIME
      CALL CLOCK
C
C RADIATION PRESSURE
      RP(1)=RP(1)+RPR(1)
      DO 320 I=2,NTAU
      RPE(I)=RPE(I)-RPR(I)
      DO 301 J=1,NTAU
      PET(I,J)=PET(I,J)-PRT(I,J)
301   CONTINUE
C
C TURBULENT PRESSURE
      PET(I,I)=PET(I,I)-PTT(I)
      PEPE(I,I)=PEPE(I,I)-PTPE(I)
310   CONTINUE
C
C CONVECTIVE EFFICIENCY
      DT(I,I)=DT(I,I)-.5*DG(I)*GT(I)
      DT(I,I-1)=DT(I,I-1)-.5*DG(I)*GT(I)
      DPE(I,I)=DPE(I,I)-.5*DG(I)*GPE(I)
      DPE(I,I-1)=DPE(I,I-1)-.5*DG(I)*GPE(I)
      DV(I)=-DG(I)*GV(I)
C
C TOTAL PRESSURE
      DPI=-DP(I+NTAU-1)
      DPE(I,I-1)=DPE(I,I-1)-DPI*PPPE(I+NTAU-1)
      DT(I,I-1)=DT(I,I-1)-DPI*PPTT(I+NTAU-1)
      DPE(I,I)=DPE(I,I)-DPI*PPPE(I)
      DT(I,I)=DT(I,I)-DPI*PPTT(I)
      RD(I)=RD(I)-DPI*RP(I)
      DP(I)=DP(I)-DPI
320   CONTINUE
C
      RPI=0.
      DV(1)=0.
      DO 321 I=1,NTAU
C WE HAVE THE MATRICES       PPP       PPPE  PPTT
C AND                        PEPP      PEPE  PET
C WE WANT TO SUBTRACT FROM PEPE AND PET THE PRODUCTS OF PEPP*PPP-INVERS WITH
C PPPE AND PET RESPECTIVELY. PPP IS BIDIAGONAL WITH UNITY ON THE DIAGONAL AND
C MINUS UNITY ON THE SUBDIAGONAL. ITS INVERS IS A MATRIX WITH UNITY EVERYWHERE
C UNDER AND ON THE DIAGONAL. PEPP IS MINUS UNITY. THUS THE PRODUCT OF PEPP*PPP-
C INVERS WITH PPPE IS MATRIX OF THE FOLLOWING TYPE. IN EVERY COLUMN EACH
C ELEMENT IS THE SUM OF ALL ELEMENTS ABOVE THAT POINT (AND INCLUDING) IN THE
C PPPE MATRIX. SIMILARILY FOR THE PPTT-MATRIX.
      RPI=RPI+RP(I)
      RPE(I)=RPE(I)+RPI
      RD(I)=RD(I)-DP(I)*RPI
      Y=PPPE(I)
      YY=PPTT(I)
      YY=YY+PRT(1,I)
      JMIN=MAX0(I,1)
      DO 321 J=JMIN,NTAU
      IF(J.NE.I+1) GO TO 322
      Y=Y+PPPE(I+NTAU)
      YY=YY+PPTT(I+NTAU)
322   CONTINUE
      PEPE(J,I)=PEPE(J,I)+Y
      PET(J,I)=PET(J,I)+YY
      DPE(J,I)=DPE(J,I)-DP(J)*Y
      DT(J,I)=DT(J,I)-DP(J)*YY
321   CONTINUE
C
C GRADIENT DIFFERENCE
      DO 350 I=2,NTAU
      DVI=DV(I)
      VDI=VD(I)
      FCDI=FCD(I)
      RFC(I)=RFC(I)-FCDI*RD(I)
      RV(I)=RV(I)-VDI*RD(I)
      VVVI=MAX(0.5D+0,1.0D+0-VDI*DVI)
      FCV(I)=FCV(I)-FCDI*DVI
      DO 340 J=1,NTAU
      VT(I,J)=VT(I,J)-VDI*DT(I,J)
      VPE(I,J)=VPE(I,J)-VDI*DPE(I,J)
      FCT(I,J)=FCT(I,J)-FCDI*DT(I,J)
      FCPE(I,J)=FCPE(I,J)-FCDI*DPE(I,J)
340   CONTINUE
C
C TURBULENT VELOCITY
      PEVI=-PTV(I)
      RV(I)=RV(I)/VVVI
      FCVI=FCV(I)
      RFC(I)=RFC(I)-FCVI*RV(I)
      RPE(I)=RPE(I)-PEVI*RV(I)
      DO 350 J=1,NTAU
      VT(I,J)=VT(I,J)/VVVI
      VPE(I,J)=VPE(I,J)/VVVI
      FCT(I,J)=FCT(I,J)-FCVI*VT(I,J)
      FCPE(I,J)=FCPE(I,J)-FCVI*VPE(I,J)
      PET(I,J)=PET(I,J)-PEVI*VT(I,J)
      PEPE(I,J)=PEPE(I,J)-PEVI*VPE(I,J)
350   CONTINUE
C
C TIME
      CALL CLOCK
C
C CONVECTIVE FLUX
      DO 360 I=1,NTAU
      RT(I)=RT(I)-RFC(I)
      DO 361 J=1,NTAU
      TTT(I,J)=TTT(I,J)-FCT(I,J)
      TPE(I,J)=TPE(I,J)-FCPE(I,J)
361   CONTINUE
      IF(I.GT.MIHAL) GO TO 360
      RT(I)=RT(I)+RFC(I+1)
      DO 362 J=1,NTAU
      TTT(I,J)=TTT(I,J)+FCT(I+1,J)
      TPE(I,J)=TPE(I,J)+FCPE(I+1,J)
362   CONTINUE
360   CONTINUE
C
C ELECTRON PRESSURE
      CALL MATINV(PEPE,NTAU)
      CALL MULT(PET,PEPE,PET,SCRATC,NTAU,NTAU)
      CALL MULT(RPE,PEPE,RPE,SCRATC,NTAU,1)
      DO 374 I=1,NTAU
      DO 374 J=1,NTAU
      SUMA=0.
      DO 375 L=1,NTAU
      SUMA=SUMA+TPE(I,L)*PET(L,J)
375   CONTINUE
      TTT(I,J)=TTT(I,J)-SUMA
      RT(I)=RT(I)-TPE(I,J)*RPE(J)
374   CONTINUE
C
C TIME
      CALL CLOCK
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C BACKSUBSTITUTION IN GAUSS ELIMINATION SCHEME.
C
C INITIATE
      CALL MATINV(TTT,NTAU)
      DO 400 I=1,NTAU
      PE(I)=RPE(I)
      FC(I)=RFC(I)
      V(I)=RV(I)
      PT(I)=0.
      PR(I)=RPR(I)
      P(I)=RP(I)
C
C SOLVE FOR TEMPERATURE CORRECTION
      T(I)=0.
      DO 400 J=1,NTAU
      T(I)=T(I)+TTT(I,J)*RT(J)
400   CONTINUE
      WRITE(6,68) ITER
      WRITE(6,67) (T(I),I=1,NTAU)
      WRITE(7,68) ITER
      WRITE(7,67) (T(I),I=1,NTAU)
C---
      TCORMX=ABS(T(1))
      DO 405 I=2,NTAU
 405  TCORMX=MAX(TCORMX,ABS(T(I)))
      IF(TCORMX.LE. tconv)  ITSTOP=.TRUE.
      PRINT406, TCORMX,ITER
406   FORMAT(' Max corr. to T wanted was',F6.1,' K for iteration',I2)
C
C CHECK T CORR
C
      IF(KORT.LE.0) GO TO 412 !no limiting in temperature correction
C
C LIMITING T CORRECTION TO ONE TENTH OF THE LOCAL TEMPERATURE
C
      DO 401 I=1,NTAU
       T(I)=T(I)/SQRT(1.+100.*(T(I)/TT(I))**2)
 401  CONTINUE
C
      IF(KORT.LE.1) GO TO 412 !no temperature inversion fix in corr
C
C Make sure correction doesn't impose a temperature inversion 
C (UGJ/900614)
C
      DO 402 I=2,NTAU
       IF (TT(I)+T(I).LE.TT(I-1)+T(I-1)) THEN
         IF (I.EQ.2) THEN
           T(I)=TT(I-1)+T(I-1)-TT(I)+2.
          ELSE
           T(I)=1.1*(TT(I-1)+T(I-1))-0.1*(TT(I-2)+T(I-2))-TT(I)
         END IF
       END IF
 402  CONTINUE
C
      IF(KORT.LE.2) GO TO 412 !no .le.-limitits in temperature correction
C
C Avoid growing oscillations in temperature correction
C (UGJ/900614)
C
      WRITE(7,*) 'ITER,NFIRST in SOLVE_sph = ',ITER,NFIRST
      IF (NFIRST.EQ.0) THEN
         DO 403 I=1,NTAU
403      TLAST(I)=T(I)
         NFIRST=1
      END IF
      DO 404 I=1,NTAU
      TAVDIF=0.
      IA(1)=MAX(1,I-2)
      IA(2)=MAX(2,I-1)
      IA(3)=MAX(3,I)
      IA(4)=MIN(MAX(4,I+1),NTAU)
      IA(5)=MIN(MAX(5,I+2),NTAU)
      DO 409 KIA=1,5
409   TAVDIF=TAVDIF+ABS( T(IA(KIA)) )
      TVD(I)=TAVDIF*0.2
      IF( T(I)/TLAST(I).LE.-1.0 .AND. ABS(T(I)).GT.TAVDIF ) THEN
        T(I)=-TLAST(I)
      END IF
404   CONTINUE
C
C
      IF(KORT.LE.3) GO TO 412 ! damping in temperature correction
C
C Damp oscillations in temperature correction
C (UGJ/900719)
C
      DO 407 I=1,NTAU
      IF ( T(I)/TLAST(I).LE.-0.6  .AND. ABS(T(I)).GT.TVD(I) ) 
     -     T(I)=-0.3*TLAST(I)
407   CONTINUE
C
C
      IF(KORT.LE.4) GO TO 412 ! damping in temperature correction
C
C Limit temperature correction to KPP
C (AB/950519)
C
      PPK = DFLOAT(KPP)
      DO 413 I=1,NTAU
      IF (ABS(T(I)).GT.ABS(PPK)) 
     -     T(I)=PPK*T(I)/ABS(T(I))
413   CONTINUE
C 
C
412   continue
C
C
      TCORMX=ABS(T(1))
      DO 414 I=2,NTAU
414   TCORMX=MAX(TCORMX,ABS(T(I)))
      PRINT415, TCORMX,ITER
415   FORMAT(' Max corr. applied to T was',F6.1,' K for iteration',I2)
C
C
      DO 408 I=1,NTAU
408   TLAST(I)=T(I)
C
      WRITE(6,69) ITER
      WRITE(6,*) ' applied corrections: '
      WRITE(6,67) (T(I),I=1,NTAU)
      WRITE(7,*) ' applied corrections: '
      WRITE(7,67) (T(I),I=1,NTAU)
C
C SUBTRACT TEMPERATURE
      DO 410 I=1,NTAU
      PT(I)=PT(I)-PTT(I)*T(I)
      P(I)=P(I)-PPTT(I)*T(I)
      IF(I.GT.1) P(I)=P(I)-PPTT(I+NTAU-1)*T(I-1)
      DO 410 J=1,NTAU
      V(I)=V(I)-VT(I,J)*T(J)
      PE(I)=PE(I)-PET(I,J)*T(J)
      FC(I)=FC(I)-FCT(I,J)*T(J)
      PR(I)=PR(I)-PRT(I,J)*T(J)
410   CONTINUE
C
C SUBTRACT ELECTRON PRESSURE
      DO 420 I=1,NTAU
      PT(I)=PT(I)-PTPE(I)*PE(I)
      P(I)=P(I)-PPPE(I)*PE(I)
      IF(I.GT.1) P(I)=P(I)-PPPE(I+NTAU-1)*PE(I-1)
      DO 420 J=1,NTAU
      V(I)=V(I)-VPE(I,J)*PE(J)
      FC(I)=FC(I)-FCPE(I,J)*PE(J)
420   CONTINUE
C
C TOTAL PRESSURE
      P(1)=P(1)+PR(1)-RPR(1)
      DO 425 I=2,NTAU
      P(I)=P(I)+P(I-1)
425   CONTINUE
C
C SUBTRACT VELOCITY
      DO 430 I=1,NTAU
      PT(I)=PT(I)-PTV(I)*V(I)
430   CONTINUE
C
C PRINT
      IF (NPRINT.GE.2) THEN
      IF(PFE.OR.PFD) WRITE(6,48) FNORD,ITER,IDRAB1,IDRAB2,IDRAB3,
     &                           IDRAB4,IDRAB5,IDRAB6
      IF(PFE.OR.PFD) WRITE(6,56)
      IF(PFE.OR.PFD) WRITE(6,52)
      END IF
      DO 440 I=1,NTAU
C
C SOLVE FOR CONVECTIVE EFFICIENCY AND GRADIENT DIFFERENCE
      I1=MAX0(I-1,1)
      GI=-.5*(GT(I)*(T(I)+T(I1))+GPE(I)*(PE(I)+PE(I1)))-GV(I)*V(I)
      GG(I)=GG(I)+GI
      D(I)=D(I)-DG(I)*GI-DPS(I)*P(I)-DPES(I)*PE(I)-DTS(I)*T(I)
      IF(I.GT.1) D(I)=D(I)-DPS(I+NTAU-1)*P(I-1)-DPES(I+NTAU-1)
     &*PE(I-1)-DTS(I+NTAU-1)*T(I-1)
      IF(PFE.OR.PFD.AND.NPRINT.GE.2) 
     *  WRITE(6,51) I,TAU(I),PR(I),PT(I),P(I),GI,D(I),V(I)
     * ,FC(I),PE(I),T(I),I
440   CONTINUE
C
C TIME
      CALL CLOCK
C
C APPLY CORRECTIONS
      IPRESS=0
      DO 450 I=1,NTAU
      TT(I)=TT(I)+T(I)
      VV(I)=MAX(VV(I)+V(I),0.0D+0)
      FFC(I)=FFC(I)+FC(I)
      PPT(I)=MAX(PPT(I)+PT(I),0.0D+0)
      PPT(I)=MIN(PPT(I),0.5*PP(I))
      PPR(I)=PPR(I)+PR(I)
      DD(I)=DD(I)+D(I)
C
C IF TOO VIOLENT CHANGES TO PPE OR PP, SET IPRESS FOR AN EXTRA
C PRESSURE INTEGRATION AFTER CORRECTIONS HAVE BEEN APPLIED.
      IF(ABS(P(I)/PP(I)).LT.0.5.AND.ABS(PE(I)/PPE(I)).LT.0.5) GOTO 451
      IPRESS=1
      GOTO 450
C
451   PPE(I)=PPE(I)+PE(I)
      PP(I)=PP(I)+P(I)
450   CONTINUE
C
C WE MUST NOT USE TRYCK TOWARDS THE END OF ITERATIONS, BECAUSE ROSS(K)
C IS DEFINED BY OPAC CALLS, RATHER THAN ROSSOP CALLS.
      IF(IPRESS.EQ.1) CALL TRYCK_sph
C
C PRINT PRESENT STATE OF THE ATMOSPHERE
*************************************************************************
      ENTRY PRESNT_sph
*
      IF(.NOT.PFE) RETURN
      INORD=IEDIT+10*IVERS
      FNORD=.1*INORD
      IF (NPRINT.GE.2) THEN
      WRITE(6,48) FNORD,ITER,IDRAB1,IDRAB2,IDRAB3,IDRAB4,
     &            IDRAB5,IDRAB6
      WRITE(6,53)
      WRITE(6,52)
      END IF
      I=0
      Y=0.
      TT0=(TAU(2)*TT(1)-TAU(1)*TT(2))/(TAU(2)-TAU(1))
      IF (NPRINT.GE.2) THEN
      WRITE(6,51) I,Y,PPR(1),Y,PPR(1),Y,Y,Y,Y,Y,TT0,I
      WRITE(6,51) (I,TAU(I),PPR(I),PPT(I),PP(I),GG(I),DD(I)
     * ,VV(I),FFC(I),PPE(I),TT(I),I,I=1,NTAU)
      END IF
C
C TIME
      CALL CLOCK
      RETURN
C
C FORMATS
45    FORMAT(' TIME',I6,' MSEC')
48    FORMAT('1SCMARCS',F5.1,5X,'SOLVE/SPH(53)  7-NOV-80',5X,
     & '.....................................ITERATION',I3,
     & 5X,6A4/)
49    FORMAT(13('1234567890'),'123')
50    FORMAT(4(/2X,1P10E13.5))
51    FORMAT(I3,1P10E12.4,I4)
52    FORMAT(T7,'TAU',T19,'PRAD',T31,'PTURB',T43,'PTOT',T55,'GAMMA',
     *T67,'DELTA',T79,'VCONV',T91,'FCONV',T103,'PE',T115,'TEMP')
53    FORMAT(' STATE OF MODEL ATMOSPHERE')
54    FORMAT(' *=R.H. SIDES     *',23X,'*',23X,'*',
     ,11X,'*',11X,'*',11X,'*',11X,'*')
55    FORMAT(6X,'TAU',9X,'HSCALE',6X,'ADIA',8X,'GRAD',8X,'CP',10X,
     *'Q',11X,'PG',10X,'RO',10X,'PGPE',8X,'PGT')
56    FORMAT(' CORRECTIONS')
57    FORMAT(' THERMODYNAMICALS')
58    FORMAT(I4,F9.1,F8.1,1P2E10.3,0PF6.3,1PE10.3,0PF8.3,
     * F7.0,4F5.1,1PE10.3,2(0PF6.1),2I4,F6.3)
59    FORMAT('   J',4X,'WAVEL',2X,'WEIGHT',2X,'F(WAVEL)',3X,
     * 'CONTRIB',1X,'WAVEN',2X,'F(WAVEN)',4X,'MAGN',2X,'TRAD',3X,
     * 'X01',2X,'X25',2X,'S01',2X,'S25',3X,'CONTRIB',2X,
     * 'DX01  DX25 OPAC TRAN ALGB')
60    FORMAT('0SAVED VALUES FROM LOGICAL UNIT',I3)
61    FORMAT(' STARTING VALUES, QTEMP=',F5.2)
62    FORMAT(' PRESSURE EQUATION')
63    FORMAT(6X,'TAU',9X,'PTAU',8X,'ROSS',8X,'ROSSPE',6X,'ROSST',7X,
     &'KAP5',8X,'ROSSP',7X,'TAU5',8X,'TAUP',8X,'RADIUS')
64    FORMAT(' K=',I2,4X,'NEW V =',1PE10.3)
65    FORMAT('0TOTAL FLUX=',1PE11.4,' ERGS/CM**2/S    FLUX/PI=',1PE11.4,
     & ' ERGS/CM**2/S    TEFF=',0PF6.0,' K')
67    FORMAT(1X,10F7.1)
68    FORMAT(' ITERATION',I3,'     COMPUTED CORRECTION')
69    FORMAT(' ITERATION',I3,'     APPLIED  CORRECTION')
      END
C
      SUBROUTINE TRANEQ_sph
      implicit real*16 (a-h,o-z)
C
C TRANEQ SOLVES THE TRANSFER EQUATION INCLUDING CONTINUUM SCATTERING.
C FEATURES:
C
C 1. CANNONS PERTURBATION TECHNIQUE IS USED ON THE ANGULAR QUADRATURE.
C    THE BASIC IDEA IN THIS TECHNIQUE IS TO REPLACE THE INVERSION OF
C    A COMPLICATED (MMU ORDER) OPERATOR WITH THE INVERSION OF A SIMPLE
C    OPERATOR (ONE POINT=EDDINGTON APPROXIMATION), PLUS ITERATION ON
C    THE ERROR.
C 2. A TRICK DUE TO ROBERT STEIN (PRIV. COMM., 1979) IS USED TO
C    ELIMINATE THE NEED FOR DOUBLE PRECISION STORAGE OF THE MATRIX
C    ELEMENTS. THE IDEA IS TO STORE THE (SMALL) SUM OF THE THREE
C    MATRIX ELEMENTS ON A ROW, INSTEAD OF THE (LARGE) DIAGONAL ELE-
C    MENT.
C 3. THE SOLUTION IS A CUBIC SPLINE, RATHER THAN A PIECE-WISE
C    QUADRATIC FUNCTION. THIS IS ACCOMPLISHED WITH THE CORRECTION
C    TERMS AD AND BD IN SUBROUTINE TRANFR.
C 4. A BOUNDARY CONDITION WHICH INCLUDES AN ESTIMATED INFALLING
C    RADIATION MAKES THE SOLUTION GOOD ALSO FOR VALUES OF X+S
C    LARGE COMPARED WITH 1./TAU(1). A LOGARITHMIC TAU-SCALE
C    SHOULD BE USED.
C
C THIS VERSION OF TRANEQ IS COMPATIBLE WITH PREVIOUS TRANEQ'S.
C 79.06.21 *NORD*
C
      include 'parameter.inc'
C
      PARAMETER (ITMAX=12)
      COMMON /CTRAN/X(NDP),S(NDP),BPLAN(NDP),XJ(NDP),XH(NDP),XK(NDP)
     & ,FJ(NDP),SOURCE(NDP),TAUS(NDP),DTAUS(NDP),JTAU0,JTAU1,ISCAT
      COMMON /TAUC/TAU(NDP),DTAULN(NDP),JTAU   
      COMMON /ROSSC/ROSS(NDP),CDUMM(NDP) 
      COMMON /RHOC/RHO(NDP)
      COMMON /SPACE2/ERROR(NDP),FACT(NDP),DSO(NDP),
     &  P(NDP),DUM(NDP,3),
     &  SP1(NDP,NRAYS),SP2(NDP,NRAYS),SP3(NDP,NRAYS),AD(NDP,NRAYS),
     &  BD(NDP,NRAYS),EX(NRAYS),
     &  PIMPAC(NRAYS),
     &  TAUT(NDP),DTAUT(NDP),
     &  PFEAU(NRAYS,NDP),XMU(NRAYS,NDP),MMU(NDP),KIMPAC(NRAYS),NIMPAC
      COMMON /CSPHER/DIFLOG,RADIUS,RR(NDP),NCORE
      COMMON /TRDBUG/IDEBUG
      LOGICAL DEBUG
      DIMENSION A(ITMAX)
      DATA DEBUG/.FALSE./
C
C INITIATE, XJ IS SET TO THE DIFFUSION LIMIT VALUE
109   IDEBUG=0
      DO 100 K=1,JTAU
      IF (K.GT.1) GO TO 101
      DTAUB=(TAU(2)-TAU(1))*0.5*(X(2)+S(2)+X(1)+S(1))
      DBPLB=(BPLAN(2)-BPLAN(1))/DTAUB
      D2BPL=0.
      GO TO 102
101   IF (K.EQ.JTAU) GO TO 102
      DTAUA=DTAUB
      DTAUB=(TAU(K+1)-TAU(K))*0.5*(X(K)+S(K)+X(K+1)+S(K+1))
      DBPLA=DBPLB
      DBPLB=(BPLAN(K+1)-BPLAN(K))/DTAUB
      DTAUC=0.5*(DTAUA+DTAUB)
      D2BPL=(DBPLB-DBPLA)/DTAUC
102   XH(K)=D2BPL
      XJ(K)=BPLAN(K)+0.333333*(X(K)+S(K))/X(K)*D2BPL
      XK(K)=0.333333*BPLAN(K)+(0.2+0.111111*S(K)/X(K))*D2BPL
      FJ(K)=1.
100   SOURCE(K)=BPLAN(K)
CUGJ      IF (DEBUG) write(6,103) X,S,BPLAN,XJ,XH,XK
CUGJ103   FORMAT(' X,S,B,XJ,D2B,XK='/(/4(1X,1P,10E12.4/)))
C
C CALCULATE THE MATRIX ELEMENTS
      CALL TRRAYS_sph
      CALL TRANFR_sph
      CALL FORMAL_sph
      NIMP1=NIMPAC+1
C_temp      IF (DEBUG) PRINT 132,XJ,SOURCE,ERROR,FJ
C_temp     & ,((PFEAU(I,K),K=1,NDP),I=1,NIMP1)
      IF (IDEBUG.GT.1) GO TO 150
C
C ITERATION LOOP
      DO 110 IT=1,ITMAX
110   A(IT)=0.
      DO 140 IT=1,ITMAX
      ITM=IT
C
C SOLVE THE CONTINUUM SCATTERING PROBLEM IN THE EDDINGTON APPROXIMATION
      CALL TRANSC_sph
      CALL SCATTR_sph
C_temp      IF (DEBUG) PRINT 122,EX(ISCAT),DUM,P,DTAUS
C_temp 122   FORMAT(' EX,SP1,SP2,SP3,P,DTAUS=',E10.4/(/4(1X,1P,10E12.4/)))
C
C CORRECTION TO THE SOURCE FUNCTION
      DO 120 K=1,JTAU1
      P(K)=ERROR(K)+P(K)*FJ(K)*S(K)/(X(K)+S(K))
      A(IT)=MAX(A(IT),ABS(P(K)/SOURCE(K)))
120   CONTINUE
C
C CHECK ERROR IN SOURCE FUNCTION
      IF (A(IT).LT.0.001) GO TO 141
      DO 130 K=1,JTAU1
130   SOURCE(K)=SOURCE(K)+P(K)
C
C SOLVE THE TRANSFER EQUATION WITH GIVEN SOURCE FUNCTION
      CALL FORMAL_sph
      NTAU=KIMPAC(ISCAT)
C
C NOTE THAT FJ() SHOULD ONLY BE PICKED UP ABOVE JTAU0.  THE ISCAT
C BECOMES TO INCLINED BELOW JTAU0.
      DO 131 K=1,JTAU0
131   FJ(K)=XJ(K)/PFEAU(ISCAT,K)
C_temp      IF (DEBUG) write(6,132) XJ,SOURCE,ERROR,FJ
C_temp     & ,((PFEAU(I,K),K=1,NDP),I=1,NIMP1)
132   FORMAT(' XJ,SO,ERR,FJ,PF='/(/4(1X,1P,10E12.4/)))
      IF (IDEBUG.GT.1) GO TO 150
C
C END OF ITERATION LOOP
140   CONTINUE
C
C NOT CONVERGED
      IDEBUG=1
C      WRITE (13) JTAU,TAU,X,S,BPLAN,RADIUS,RR,RHO,ROSS
      write(6,*)' not converged in traneq_sph after ',itm,' iterations'
      WRITE(6,142) (A(IT),IT=1,ITM)
142   FORMAT(' MAXFEL =',1P11E11.2)
C
C CONVERGED, IF IN FIRST ITERATION, HAVE TO CALCULATE FJ().
141   IF (ITM.GT.1) GO TO 143
      NTAU=KIMPAC(ISCAT)
      DO 144 K=1,NTAU
144   FJ(K)=XJ(K)/PFEAU(ISCAT,K)
143   CONTINUE
C
C CALCULATE MOMENTS, AND CHECK DEBUG CONTROL
      CALL TRMOM_sph
      IF (DEBUG.AND.IDEBUG.GT.1) STOP ' stop in traneq_sph at 150 '
150   continue
      IF (DEBUG.AND.IDEBUG.EQ.1) IDEBUG=0
      DEBUG=IDEBUG.GT.1
      IF (DEBUG) GO TO 109
C
      RETURN
      END
C
      SUBROUTINE TRANFR_sph
      implicit real*16 (a-h,o-z)
C
C FORMAL SOLVES THE TRANSFER EQUATION WITH GIVEN SOURCE FUNCTION 'SOURCE'.
C 'ERROR' IS THE RESULTING ERROR IN THE DEFINITION OF THE CONTINUUM
C SCATTERING SOURCE FUNCTION. TRANSFR CALCULATES THE MATRIX ELEMENTS
C OF THE PROBLEM. INTENSITIES AT TAU=0 ARE RETURNED IN /CSURF/.
C 80.08.05 *NORD*
C
      include 'parameter.inc'
C
      COMMON /CTRAN/X(NDP),S(NDP),BPLAN(NDP),XJ(NDP),XH(NDP),XK(NDP),
     &  FJ(NDP),SOURCE(NDP),TAUS(NDP),DTAUS(NDP),JTAU0,JTAU1,ISCAT
      COMMON /TAUC/TAU(NDP),DTAULN(NDP),JTAU
      COMMON /SPACE2/ERROR(NDP),FACT(NDP),DSO(NDP),
     &  P(NDP),DUM(NDP,3),
     &  SP1(NDP,NRAYS),SP2(NDP,NRAYS),SP3(NDP,NRAYS),AD(NDP,NRAYS),
     &  BD(NDP,NRAYS),EX(NRAYS),
     &  PIMPAC(NRAYS),
     &  TAUT(NDP),DTAUT(NDP),
     &  PFEAU(NRAYS,NDP),XMU(NRAYS,NDP),MMU(NDP),KIMPAC(NRAYS),NIMPAC
      COMMON /CSPHER/DIFLOG,RADIUS,RR(NDP),NCORE
      COMMON /CSURF/HSURF,YSURF(NRAYS)
      COMMON /ROSSC/ROSS(NDP),CDUMM(NDP) /RHOC/RHO(NDP)
      COMMON /TRDBUG/IDEBUG
C
      DIMENSION FUN(NRAYS),DER(NRAYS*2),DMU(NRAYS+1)
      DIMENSION TAULOG(NDP)
C
C MU LOOP
      DO 131 I=1,NIMPAC
      NTAU=KIMPAC(I)
      NTAU1=NTAU-1
C
C CALCULATE DTAUT ALONG THE RAY
      ZOLD=0.0
      DO 100 K=1,NTAU
      Z=SQRT(RR(NTAU-K+1)**2-PIMPAC(I)**2)
      XMU(I,NTAU-K+1)=-Z/RR(NTAU-K+1)
      IF (K.EQ.1) GO TO 100
      DZ=Z-ZOLD
      DZDR=DZ/(RR(NTAU-K+1)-RR(NTAU-K+2))
      DTAUT(NTAU-K+2)=DZDR*0.5*(X(NTAU-K+1)+S(NTAU-K+1)
     & +X(NTAU-K+2)+S(NTAU-K+2))*(TAU(NTAU-K+2)-TAU(NTAU-K+1))
100   ZOLD=Z
      TAUT(1)=DZDR*(X(1)+S(1))*TAU(1)
      DO 101 K=2,NTAU
101   TAUT(K)=TAUT(K-1)+DTAUT(K)
C
C SAVE THE TAU SCALE FOR THE RADIAL RAY (I=1).  THIS IS USED IN
C THE EXTRA/INTERPOLATION OF THE PFEAU FOR MU=0., FURTHER DOWN.
C
      IF (I.EQ.1) THEN
        DO 102 K=1,NTAU
102       TAULOG(K)=log(TAUT(K))
      ENDIF
C
C K=1
      A=1./DTAUT(2)
      B=A**2
      SP2(1,I)=1.+2.*A
      SP3(1,I)=-2.*B
      EX(I)=TAUT(1)*(1.-0.5*TAUT(1)*(1.-0.333333*TAUT(1)))
      IF (TAUT(1).GT.0.1) EX(I)=1.-EXP(-TAUT(1))
      SP2(1,I)=SP2(1,I)/(1.+2.*A*EX(I))
      SP3(1,I)=SP3(1,I)/(1.+2.*A*EX(I))
C
C K=2,NTAU-1
      DO 110 K=2,NTAU1
      DTAUC=0.5*(DTAUT(K)+DTAUT(K+1))
C...      AD(K,I)=0.1666667*DTAUT(K)/DTAUC
C...      BD(K,I)=0.1666667*DTAUT(K+1)/DTAUC
      AD(K,I)=0.
      BD(K,I)=0.
      SP1(K,I)=-1./(DTAUT(K)*DTAUC)+AD(K,I)
      SP2(K,I)=1.
110   SP3(K,I)=-1./(DTAUT(K+1)*DTAUC)+BD(K,I)
C
C K=NTAU
      AD(NTAU,I)=0.0
      BD(NTAU,I)=0.0
      SP1(NTAU,I)=0.0
      SP2(NTAU,I)=1.0
      SP3(NTAU,I)=0.0
      IF (I.LE.NCORE) GO TO 120
      AD(NTAU,I)=0.3333333
      SP1(NTAU,I)=0.3333333-2./DTAUT(NTAU)**2
120   CONTINUE
C
C ELIMINATE SUBDIAGONAL
      DO 130 K=1,NTAU1
      SP1(K,I)=-SP1(K+1,I)/(SP2(K,I)-SP3(K,I))
      SP2(K+1,I)=SP2(K+1,I)+SP1(K,I)*SP2(K,I)
130   SP2(K,I)=SP2(K,I)-SP3(K,I)
131   CONTINUE
C
C FIND A GOOD RAY FOR TRANSC.  THE XMU VALUES ARE INCREASING FROM
C -1.0 TOWARDS 0.0, AND WE THUS ALWAYS FIND A VALUE FOR ISCAT.
C
C      IMAX=MIN0(MMU(JTAU0),NCORE)
C
C NCORE MUST BE SMALLER THAN MMU(K) FOR ALL K, SO EFFECTIVELY THE
C MIN0(MMU(),NCORE) IS EQUAL TO NCORE.  SHOULD ONE REQUIRE ISCAT TO
C BE A CORE RAY?  PROBABLY NOT, SINCE WITH A LARGE TAUM (ALLOWED),
C THE ISCAT RAY IS FORCED TOWARDS SMALLER IMPACT PARAMETERS, WHICH
C DEGRADES PERFORMANCE.
C
      IMAX=MMU(JTAU0)
      TMP1=1.
      DO 132 I=1,IMAX
C
C      IF (XMU(I,JTAU0).LT.-0.577) ISCAT=I
C
C  Better to look for the ray which is closest to the Eddington angle
C
	TMP2=ABS(XMU(I,JTAU0)+0.577)
	IF (TMP2.LT.TMP1) THEN
	  TMP1=TMP2
	  ISCAT=I
	ENDIF
132   CONTINUE
      RETURN
C---------------------------------------------------------------------
C
      ENTRY FORMAL_sph
C
C MU LOOP
      DO 170 I=1,NIMPAC
      NTAU=KIMPAC(I)
      NTAU1=NTAU-1
C
C INITIATE
      P(1)=SOURCE(1)
      DO 140 K=2,NTAU
140   P(K)=(1.-AD(K,I)-BD(K,I))*SOURCE(K)+AD(K,I)*SOURCE(K-1)+
     & BD(K,I)*SOURCE(K+1)
      IF(I.LE.NCORE) P(JTAU1)=SOURCE(JTAU1)+XMU(I,JTAU1)**2*XH(JTAU1)
C
C ACCUMULATE RIGHT HAND SIDE
      DO 150 K=1,NTAU1
150   P(K+1)=P(K+1)+SP1(K,I)*P(K)
C
C BACKSUBSTITUTE
      PFEAU(I,NTAU)=P(NTAU)/SP2(NTAU,I)
      DO 160 K=1,NTAU1
      PFEAU(I,NTAU-K)=(P(NTAU-K)-
     & SP3(NTAU-K,I)*PFEAU(I,NTAU-K+1))/SP2(NTAU-K,I)
      IF (PFEAU(I,NTAU-K).LE.0.0) GO TO 230
160   CONTINUE
C
C END MU LOOP
      YSURF(I)=2.*(1.-EX(I))*PFEAU(I,1)+EX(I)**2*SOURCE(1)
      FUN(I)=-XMU(I,1)*(PFEAU(I,1)-SOURCE(1)*EX(I))
      IF (YSURF(I).LE.0.0) GO TO 231
170   CONTINUE
C
C  INTERPOLATE TO PFEAU AT MU=0 FOR THOSE K THAT HAVE NO RAY
*  THIS IS THE ORIGINAL CODE, WHICH EXTRAPOLATES IN MU, FOR CONSTANT K.
*  THIS WAS FOUND TO PRODUCE BAD RESULTS FOR THE FLUX, WHICH IS A
*  DERIVATIVE OF XK, WHICH IS THE SECOND MOMENT OF THE FUNCTION
*  BEING EXTRAPOLATED TO ZERO.
*      DO 181 K=1,JTAU1
*      II=MMU(K)
*      IF (KIMPAC(II).EQ.K) GO TO 181
*      PX=-XMU(II-2,K)/(XMU(II-1,K)-XMU(II-2,K))
*      QX=1.-PX
*      PFEAU(II,K)=EXP(log(PFEAU(II-2,K))*QX+log(PFEAU(II-1,K))*PX)
*181   CONTINUE
*
*  INTERPOLATE/EXTRAPOLATE IN K FOR PFEAU AT MU=0
*  THIS SHOWS WHAT IS INTENDED TO BE DONE, IF USING AN EXTERNAL
*  INTERPOLATION ROUTINE.
*       DO 181 I=1,NIMPAC
*         K=KIMPAC(I)
*         F(I)=PFEAU(I,K)
*         X(I)=log(TAU(K))
*181   CONTINUE
*      DO 182 K=1,JTAU1
*        XX(K)=log(TAU(K))
*182   CONTINUE
*      CALL INTERP(NIMPAC,X,F,JTAU1,XX,FF)
*      DO 183 K=1,JTAU1
*        I=MMU(K)
*        PFEAU(I,K)=FF(K)
*183   CONITNUE
C
C  DO IT INLINE:  EXTRAPOLATE/INTERPOLATE THE PFEAU IN log(TAU), FOR MU=0.
C   START WITH 3 OUTERMOST POINTS (YOU PROBABLY NOTICED THAT KIMPAC(NIMPAC)=4)
C   JUMP OVER THE RAYS (WHERE PFEAU(MU=0.) IS ALREADY CALCULATED), AND
C   INTERPOLATE THE REST OF THE TIME.
C
      IF (NCORE.EQ.NIMPAC) THEN
C
C  TEMPORARY SECURITY TO AVOID DIVIDE BY ZERO ERRORS IN COMPUTATION OF PX,
C  WHEN THE OPACITY IS SO LARGE THAT THERE ARE ONLY NCORE RAYS.
C  IN THAT CASE, ONE USES THE OLD INTERPOLATION ROUTINE
C
      write(6,*) ' NCORE = NIMPAC. INTERPOLATION IN MU USED'
      DO 1810 K=1,JTAU1
        II=MMU(K)
        IF (KIMPAC(II).EQ.K) GO TO 1810
        PX=-XMU(II-2,K)/(XMU(II-1,K)-XMU(II-2,K))
        QX=1.-PX
        PFEAU(II,K)=EXP(log(PFEAU(II-2,K))*QX+log(PFEAU(II-1,K))*PX)
      if(pfeau(ii,k).lt.0) 
     &       write(6,*) 'k',k,ii,pfeau(ii,k),' <0 pfeau'
1810   CONTINUE

      ELSE
C
C HERE IS THE NORMAL ONE.  EXTRAPOLATE FOR K<4, THEN INTERPOLATE
C
      I=NIMPAC
      DO 181 K=1,JTAU1
        IF (K.NE.KIMPAC(I)) THEN
          PX=(TAULOG(K)-TAULOG(KIMPAC(I)))/
     &       (TAULOG(KIMPAC(I-1))-TAULOG(KIMPAC(I)))
          QX=1.-PX
CUGJ          PFEAU(MMU(K),K)=EXP( log( PFEAU(MMU(KIMPAC(I)),KIMPAC(I)) )
CUGJ     &       *QX + log( PFEAU(MMU(KIMPAC(I-1)),KIMPAC(I-1)) )*PX  )
          IF ( (PFEAU(MMU(KIMPAC(I)),KIMPAC(I)).LE.1.E-20) .OR.
     &    (PFEAU(MMU(KIMPAC(I-1)),KIMPAC(I-1)).LE.1.E-20) ) THEN
CUGJ             WRITE(6,1812) I,K,MMU(K),KIMPAC(I-1),KIMPAC(I)
CUGJ     &                     ,PFEAU(MMU(K),K)
CUGJ     &                     ,PFEAU(MMU(KIMPAC(I)),KIMPAC(I))
CUGJ     &                     ,PFEAU(MMU(KIMPAC(I-1)),KIMPAC(I-1)),PX
CUGJ     &                ,TAULOG(KIMPAC(I-1)),TAULOG(KIMPAC(I)),TAULOG(K)
             PFEAU(MMU(KIMPAC(I)),KIMPAC(I))=
     &                MAX(1.0D-99,PFEAU(MMU(KIMPAC(I)),KIMPAC(I)))
             PFEAU(MMU(KIMPAC(I-1)),KIMPAC(I-1))=
     &                MAX(1.0D-99,PFEAU(MMU(KIMPAC(I-1)),KIMPAC(I-1)))
          END IF
          PFEAU(MMU(K),K)=log( PFEAU(MMU(KIMPAC(I)),KIMPAC(I)) )
     &       *QX + log( PFEAU(MMU(KIMPAC(I-1)),KIMPAC(I-1)) )*PX
          IF(PFEAU(MMU(K),K).GE.85.) THEN
CUGJ             WRITE(6,1812) I,K,MMU(K),KIMPAC(I-1),KIMPAC(I)
CUGJ     &                     ,PFEAU(MMU(K),K)
CUGJ     &                     ,PFEAU(MMU(KIMPAC(I)),KIMPAC(I))
CUGJ     &                     ,PFEAU(MMU(KIMPAC(I-1)),KIMPAC(I-1)),PX
CUGJ     &                ,TAULOG(KIMPAC(I-1)),TAULOG(KIMPAC(I)),TAULOG(K)
1812         FORMAT(5I3,1P7E9.2)
             PFEAU(MMU(K),K)=85.
          END IF
          PFEAU(MMU(K),K)=EXP( PFEAU(MMU(K),K) )
        ENDIF
        IF (K.EQ.(KIMPAC(I-1)-1)) THEN
          I=I-1
        ENDIF
181   CONTINUE

      ENDIF

C
C CALCULATE MEAN INTENSITY
      DO 190 K=1,JTAU1
      XJ(K)=0.
      NMU=MMU(K)
      DO 191 I=2,NMU
      DMU(I)=XMU(I,K)-XMU(I-1,K)
      DER(I)=(PFEAU(I,K)-PFEAU(I-1,K))/DMU(I)
191   XJ(K)=XJ(K)+DMU(I)*(PFEAU(I,K)+PFEAU(I-1,K))
      XJ(K)=XJ(K)*6.
      NMU1=NMU-1
      DO 192 I=2,NMU1
192   XJ(K)=XJ(K)+(DMU(I+1)-DMU(I))*(DMU(I)*DER(I+1)+DMU(I+1)*DER(I))
      XJ(K)=XJ(K)+DMU(2)**2*DER(2)-DMU(NMU)**2*DER(NMU)
      XJ(K)=XJ(K)*0.083333333
      ERROR(K)=(XJ(K)*S(K)+BPLAN(K)*X(K))/(X(K)+S(K))-SOURCE(K)
190   CONTINUE
      RETURN
C--------------------------------------------------------------------
C
      ENTRY TRMOM_sph
C
C FLUX AT TAU(1)
      XH(1)=TRQUAD_sph(MMU(1),XMU,FUN,DER)
C
C  CALCULATE SECOND MOMENT XK
C
      DO 201 K=1,JTAU1
      NMU=MMU(K)
      DO 200 I=1,NMU
200   FUN(I)=PFEAU(I,K)*XMU(I,K)**2
201   XK(K)=TRQUAD_sph(NMU,XMU(1,K),FUN,DER)
C
C CALCULATE FIRST MOMENT, XH, FROM MOMENT RELATION
      DO 211 K=JTAU1+1,JTAU
211   XH(K)=( XK(K)-XK(K-1)+
     &       (XJ(K)+XJ(K-1)-3.*(XK(K)+XK(K-1)))*
     &       (RR(K-1)-RR(K))/(RR(K)+RR(K-1))    )*2.0/
     &          ( (TAU(K)-TAU(K-1)) * (X(K)+S(K)+X(K-1)+S(K-1)) )
C
C CALCULATE FIRST MOMENT BY USING R=DP/DTAU.  THIS IS MORE ACCURATE
C IN THE OPTICALLY THIN PARTS.
      ZOLD=0.0
      DO 212 K=2,JTAU1
        NMU=MMU(K)
        DO 213 I=1,NMU-1
          DZDR=(XMU(I,K)*RR(K)-XMU(I,K-1)*RR(K-1))/(RR(K-1)-RR(K))
          DTAU=DZDR*0.5*(X(K-1)+S(K-1)+X(K)+S(K))*(TAU(K)-TAU(K-1))
          DMU(I)=
     &      -(XMU(I,K)*RR(K)+XMU(I,K-1)*RR(K-1))/(RR(K)+RR(K-1))
          FUN(I)=DMU(I)*(PFEAU(I,K)-PFEAU(I,K-1))/DTAU
213     CONTINUE
        FUN(NMU)=0.
        DMU(NMU)=0.
        XH(K)=-TRQUAD_sph(NMU,DMU,FUN,DER)
212   CONTINUE
C
C SURFACE FLUX
      NMU=MMU(1)
      PX=-XMU(NMU-2,1)/(XMU(NMU-1,1)-XMU(NMU-2,1))
      QX=1.-PX
      YSURF(NMU)=EXP(log(YSURF(NMU-2))*QX+log(YSURF(NMU-1))*PX)
      DO 220 I=1,NMU
220   FUN(I)=-XMU(I,1)*YSURF(I)
      HSURF=0.5*TRQUAD_sph(NMU,XMU,FUN,DER)
      RETURN
C------------------------------------------------------------------
C
C EMERGENCY EXIT
230   KK=NTAU-K
      IDEBUG=2
C      PRINT 232,I,KK,JTAU0,JTAU1,NTAU
C232   FORMAT('0NON-POSITIVE RESULT AT I,K,J0,J1,N=',5I3)
      GO TO 233
231   KK=0
      IDEBUG=3
C      PRINT 232,I,KK,JTAU0,JTAU1,NTAU
C233   PRINT 237,NCORE,ISCAT,DIFLOG,RADIUS,EX(I)
233   continue
237   FORMAT(' NCORE,ISCAT,DIFLOG,RADIUS,EX(I)=',2I3,3G12.3)
C      WRITE (13,*) JTAU,TAU,X,S,BPLAN,RADIUS,RR,RHO,ROSS
C      PRINT 236,TAU,X,S,BPLAN,RR,RHO,ROSS,SOURCE
C     & ,(YSURF(I),(PFEAU(I,K),K=1,39),I=1,NMU)
236   FORMAT('0TAU=',4(/10E12.4)/'0X=',4(/10E12.4)/'0S=',4(/10E12.4)
     & /'0BPLAN=',4(/10E12.4)/'0RR=',4(/10E12.4)/'0RHO=',4(/10E12.4)
     & /'0ROSS=',4(/10E12.4)/'0SOURCE=',4(/10E12.4)
     & /'0YSURF,PFEAU='/(10E12.4))
      write(6,238)
c      PRINT 238
238   FORMAT(' DEBUG INFORMATION WRITTEN ON UNIT 13')
      RETURN
      END
C
      SUBROUTINE TRANSC_sph
      implicit real*16 (a-h,o-z)
C
C SCATTR SOLVES THE TRANSFER EQUATION INCLUDING CONTINUUM SCATTERING
C IN THE EDDINGTON APPROXIMATION, I.E., USING ONLY ONE RAY.
C 'ERROR' IS THE INHOMOGENEOUS TERM OF THE EQUATION, AND 'P' GIVES THE
C ESTIMATED MEAN INTENSITY CORRECTION (FJ*P). TRANSC CALCULATES THE MATRIX
C ELEMENTS FOR SCATTR.
C 79.06.21 *NORD*
C
      include 'parameter.inc'
C
      COMMON /CTRAN/X(NDP),S(NDP),BPLAN(NDP),XJ(NDP),XH(NDP),XK(NDP)
     & ,FJ(NDP),SOURCE(NDP),TAUS(NDP),DTAUS(NDP),JTAU0,JTAU1,ISCAT
      COMMON /TAUC/TAU(NDP),DTAULN(NDP),JTAU
      COMMON /SPACE2/ERROR(NDP),FACT(NDP),DSO(NDP),
     &  P(NDP),SP1(NDP),SP2(NDP),SP3(NDP),
     &  DUM(NDP,NRAYS,3),AD(NDP,NRAYS),BD(NDP,NRAYS),EX(NRAYS),
     &  PIMPAC(NRAYS),
     &  TAUT(NDP),DTAUT(NDP),
     &  PFEAU(NRAYS,NDP),XMU(NRAYS,NDP),MMU(NDP),KIMPAC(NRAYS),NIMPAC
      COMMON /CSPHER/DIFLOG,RADIUS,RR(NDP),NCORE
C
C MU LOOP
      NTAU=KIMPAC(ISCAT)
      NTAU1=NTAU-1
C
C CALCULATE TAUS ALONG THE RAY, SPIRAL AT EDDINGTON ANGLE AT DEPTH.
      Z=SQRT(RR(JTAU0)**2-PIMPAC(ISCAT)**2)
      ZOLD=Z
      DO 101 K=2,JTAU
      IF (JTAU-K+2.GT.JTAU0) GO TO 103
      Z=SQRT(RR(JTAU-K+1)**2-PIMPAC(ISCAT)**2)
      DZ=Z-ZOLD
      DZDR=DZ/(RR(JTAU-K+1)-RR(JTAU-K+2))
      GO TO 104
103   DZDR=1.732
104   DTAUS(JTAU-K+2)=DZDR*0.5*(X(JTAU-K+1)+S(JTAU-K+1)
     & +X(JTAU-K+2)+S(JTAU-K+2))*(TAU(JTAU-K+2)-TAU(JTAU-K+1))
101   ZOLD=Z
C
      TAUS(1)=DZDR*(X(1)+S(1))*TAU(1)
C
C K=1
      A=1./DTAUS(2)
      B=A**2
      SP2(1)=1.-FJ(1)*S(1)/(X(1)+S(1))+2.*A
      SP3(1)=-2.*B
      T=TAUS(1)
      EX(ISCAT)=TAUS(1)*(1.-0.5*TAUS(1)*(1.-0.333333*TAUS(1)))
      IF (TAUS(1).GT.0.1) EX(ISCAT)=1.-EXP(-TAUS(1))
      SP2(1)=SP2(1)-2.*A*EX(ISCAT)*FJ(1)*S(1)/(X(1)+S(1))
      SP2(1)=SP2(1)/(1.+2.*A*EX(ISCAT))
      SP3(1)=SP3(1)/(1.+2.*A*EX(ISCAT))
C
C K=2,NTAU-1
      DO 100 K=2,NTAU1
      DTAUC=0.5*(DTAUS(K)+DTAUS(K+1))
      SP1(K)=-1./(DTAUS(K)*DTAUC)
      SP2(K)=1.-FJ(K)*S(K)/(X(K)+S(K))
100   SP3(K)=-1./(DTAUS(K+1)*DTAUC)
C
C K=NTAU
      SP1(NTAU)=0.0
      SP2(NTAU)=X(NTAU)/(X(NTAU)+S(NTAU))
      SP3(NTAU)=0.0
C
C ELIMINATE SUBDIAGONAL
      DO 120 K=1,NTAU1
      SP1(K)=-SP1(K+1)/(SP2(K)-SP3(K))
      SP2(K+1)=SP2(K+1)+SP1(K)*SP2(K)
120   SP2(K)=SP2(K)-SP3(K)
121   CONTINUE
      RETURN
C--------------------------------------------------------------------
C
      ENTRY SCATTR_sph
C
C MU LOOP
      NTAU=KIMPAC(ISCAT)
      NTAU1=NTAU-1
C
C ACCUMULATE RIGHT HAND SIDE
      P(1)=ERROR(1)
      DO 150 K=1,NTAU1
150   P(K+1)=ERROR(K+1)+SP1(K)*P(K)
C
C BACKSUBSTITUTE
      P(NTAU)=P(NTAU)/SP2(NTAU)
      DO 160 K=1,NTAU1
160   P(NTAU-K)=(P(NTAU-K)-SP3(NTAU-K)*P(NTAU-K+1))/SP2(NTAU-K)
C
      RETURN
      END
C
      SUBROUTINE TRYCK_sph
      implicit real*16 (a-h,o-z)
C
C TRYCK IS A FAST PRESSURE INTEGRATION ROUITINE. IT IS FAST BECAUSE OF
C TWO REASONS: 1) IT INTEGRATES THE DIFFFERENTIAL EQUATION FOR LN(P)
C AS A FUNCTION OF LN(TAU). 2) IT ITERATES DIRECTLY ON THE ELECTRON
C PRESSURE, KEEPING THE NUMBER OF CALLS TO ABSKO TO A MINIMUM.
C ASSUMING A POWER LAW BEHAVIOUR OF PP,TT,PPE,ETC.: PP=C*TAU**DP,ETC., ONE
C CAN SHOW THAT DP=(1.+DT*(ROSSPE*PGT/PGPE-ROSST)/(1.+ROSSPE/PGPE).
C THE ANSATZ FOR PP IMPLIES PP(1)=TAU(1)*GRVR/(ROSS(1)*DP), WHICH
C SERVES AS A BOUNDARY CONDITION.
C 790516 *NORD*
C
      include 'parameter.inc'
C
      COMMON /STATEC/PPR(NDP),PPT(NDP),PP(NDP),GG(NDP),ZZ(NDP),DD(NDP),
     & VV(NDP),FFC(NDP),PPE(NDP),TT(NDP),TAULN(NDP),RO(NDP),NTAU,ITER
      COMMON /TAUC/TAU(NDP),DLNTAU(NDP),JTAU 
      COMMON /CG/GRAV,KONSG
      COMMON /ROSSC/ROSS(NDP),CROSS(NDP)
      COMMON /CI8/PGC,RHOC,EC
      COMMON /CSPHER/TAURAT,RADIUS,RR(NDP),NCORE
      COMMON /CMOLRAT/ FOLD(NDP,8),MOLOLD,KL
      DATA EPS,RELT,RELPE,PEDEF/1.E-3,1.E-3,1.E-3,1./
C
C START
      MSA=0
C      WRITE(7,101)
C101   FORMAT('1PRESSURE INTEGRATION'/'  K',6X,'TAU',10X,'TT',
C     & 9X,'PPE',8X,'PTOT',8X,'ROSS',10X,'DP',9X,'NABSKO')
      DT=0.
C USE 'DLNT/DLNTAU'=DT=0. TO BE COMPATIBLE WITH SOLVE. OTHERWISE
C DT=(TT(2)/TT(1)-1.)/DLNTAU(2)
      NABSKO=0
      KK=1
      IF(PPE(1).LE.0.) PPE(1)=PEDEF
C
C ITERATE ON BOUNDARY CONDITION, USING PARTIAL DERIVATIVES
      GRVR=GRAV*(RADIUS/RR(1))**2
C test for effect of constant gravity...
      IF (KONSG.EQ.1) GRVR=GRAV   
100   CONTINUE
      KL=1
      ROSS(1)=CROSS(1)*ROSSOP(TT(1),PPE(1))
      PP(1)=PGC+PPT(1)+PPR(1)
      PG=PGC
      ROSST=CROSS(1)*ROSSOP(TT(1)*(1.+RELT),PPE(1))
      PGT=PGC
      ROSSPE=CROSS(1)*ROSSOP(TT(1),PPE(1)*(1.+RELPE))
      PGPE=PGC
      PGT=(PGT/PG-1.)/RELT
      PGPE=(PGPE/PG-1.)/RELPE
      ROSST=(ROSST/ROSS(1)-1.)/RELT
      ROSSPE=(ROSSPE/ROSS(1)-1.)/RELPE
      NABSKO=NABSKO+3
      DP=(1.+DT*(ROSSPE*PGT/PGPE-ROSST))/(1.+ROSSPE/PGPE)
      DP=MAX(DP,0.1D+0)
      DLNPE=log(GRVR*TAU(1)/(PG*ROSS(1)*DP))/(PGPE+ROSSPE)
      PPE(1)=PPE(1)*EXP(DLNPE)
      IF(ABS(DLNPE).GT.EPS) GOTO 100
C
C END BOUNDARY CONDITION
      ROSS(1)=CROSS(1)*ROSSOP(TT(1),PPE(1))
      NABSKO=NABSKO+1
      PP(1)=PGC+PPT(1)+PPR(1)
C      WRITE(7,102) KK,TAU(1),TT(1),PPE(1),PP(1),ROSS(1),DP,NABSKO
C102   FORMAT(I3,6E12.5,I12)
C
C TAU LOOP
      DPE=(DP-DT*PGT)/PGPE
      DEDLNP=-(PGPE*PG/PP(1)+.5*DLNTAU(2)*GRVR*TAU(1)/(PP(1)*ROSS(1))*
     & (PGPE*PG/PP(1)+ROSSPE))
      DO 110 K=2,NTAU
      GRVR=GRAV*(RADIUS/RR(K))**2
C test for effect of constant gravity....
      IF (KONSG.EQ.1) GRVR=GRAV   
      PPE(K)=PPE(K-1)*EXP(DPE*DLNTAU(K))
      NABSKO=0
C
C ITERATION LOOP
      DLNPE=0.
111   CONTINUE
      KL=K
      ROSS(K)=CROSS(K)*ROSSOP(TT(K),PPE(K))
      PP(K)=PGC+PPT(K)+PPR(K)
      NABSKO=NABSKO+1
      ERROR=(.5*DLNTAU(K)*GRVR*(TAU(K-1)/(PP(K-1)*ROSS(K-1))+
     & TAU(K)/(PP(K)*ROSS(K)))-log(PP(K)/PP(K-1)))
C       print*,'error,dedlnp ', error,dedlnp
      CALL ZEROF(ERROR,DLNPE,DEDLNP)
      PPE(K)=PPE(K)*EXP(DLNPE)
      IF(ABS(DLNPE).GT.EPS) GOTO 111
C
C END TAU LOOP
      ROSS(K)=CROSS(K)*ROSSOP(TT(K),PPE(K))
      NABSKO=NABSKO+1
      PP(K)=PGC+PPT(K)+PPR(K)
      DP=GRVR*TAU(K)/(PGC*ROSS(K))
      DPE=log(PPE(K)/PPE(K-1))/DLNTAU(K)
C      WRITE(7,102) K,TAU(K),TT(K),PPE(K),PP(K),ROSS(K),DP,NABSKO
110   CONTINUE
C
C END
      MSB=0
      MSB=MSA-MSB
C     PRINT 120,MSB
120   FORMAT(' TRYCK_sph TIME=',I5,' MS')
      RETURN
      END
C
C
      FUNCTION TRQUAD_sph(N,X,F,W)
      implicit real*16 (a-h,o-z)
C
      DIMENSION X(N),F(N),W(2*N)
C was : dim x(1) etc...
C
C TRAPEZOIDAL QUADRATURE PLUS NEXT ORDER CORRECTION FOR NON-
C -EQUIDISTANT GRID.
      N1=N-1
      Q=0.
      DO 100 K=2,N
      W(K)=X(K)-X(K-1)
      W(N+K)=(F(K)-F(K-1))/W(K)
100   Q=Q+W(K)*(F(K-1)+F(K))
      Q=Q*6.
      DO 101 K=2,N1
101   Q=Q+(W(K+1)-W(K))*(W(K)*W(N+K+1)+W(K+1)*W(N+K))
      W1=((W(2)+0.5*W(3))*W(N+2)-0.5*W(2)*W(N+3))*2.0/(W(2)+W(3))
      WN=((W(N)+0.5*W(N1))*W(N+N)-0.5*W(N)*W(N+N1))*2.0/(W(N)+W(N1))
      Q=0.083333333*(Q+W(2)**2*W1-W(N)**2*WN)
      TRQUAD_sph=Q
      RETURN
      END
C
C
      SUBROUTINE TRRAYS_sph
      implicit real*16 (a-h,o-z)
C
      include 'parameter.inc'
C
      COMMON /CTRAN/X(NDP),S(NDP),BPLAN(NDP),XJ(NDP),XH(NDP),XK(NDP)
     & ,FJ(NDP),SOURCE(NDP),TAUS(NDP),DTAUS(NDP),JTAU0,JTAU1,ISCAT
      COMMON /TAUC/TAU(NDP),DTAULN(NDP),JTAU
      COMMON /SPACE2/ERROR(NDP),FACT(NDP),DSO(NDP),
     &  P(NDP),DUM(NDP,3),
     &  SP1(NDP,NRAYS),SP2(NDP,NRAYS),SP3(NDP,NRAYS),AD(NDP,NRAYS),
     &  BD(NDP,NRAYS),EX(NRAYS),
     &  PIMPAC(NRAYS),
     &  TAUT(NDP),DTAUT(NDP),
     &  PFEAU(NRAYS,NDP),XMU(NRAYS,NDP),MMU(NDP),KIMPAC(NRAYS),NIMPAC
      COMMON /CSPHER/TLIM,RADIUS,RR(NDP),NCORE
      COMMON /CSTYR/MIHAL,NOCONV
      COMMON /CTAUM/TAUM
*
*  Distribute rays, based on two depth indices, JTAU0 and JTAU1.
*
*  JTAU0 is the largest depth index for which TAU < SQRT((X+S)/X),
*  and represents the "surface", where radiation is released.
*  Above JTAU0, linearized perturbations in the radiation field are
*  represented by a single ray, which is chosen in TRANFR as the ray
*  which is closest to the Eddington angle 1./sqrt(3.) at JTAU0.
*
*  JTAU1 is the largest depth index for which TAU < TAUM*SQRT((X+S)/X),
*  and represents the border to the "core", where diffusion is a good
*  approximation.  The detailed formal solution of the radiative transfer
*  on a set of parallel rays is only performed above JTAU1.  TAUM should
*  be at least 50 to 100.
*
      TAUT(1)=TAU(1)*(X(1)+S(1))
      JTAU0=1
      JTAU1=1
      ICASE=1
      DO 100 K=2,JTAU
      IF (X(K)*BPLAN(K).LE.1.E-20) write(6,*)' X,BPL= ',X(K),BPLAN(K)
      TAUT(K)=TAUT(K-1)+0.5*(X(K)+S(K)+X(K-1)+S(K-1))*(TAU(K)-TAU(K-1))
      GO TO (101,102,103,104,105),ICASE
101   IF (TAUT(K).LT.1.0) GO TO 105
      ICASE=2
102   IF (TAUT(K).LT.SQRT((X(K)+S(K))/X(K))) GO TO 105
      ICASE=3
103   IF (TAUT(K).LT.TAUM) GO TO 105
      ICASE=4
104   IF (TAUT(K).LT.TAUM*SQRT((X(K)+S(K))/X(K))) GO TO 105
      ICASE=5
105   IF (ABS(XH(K))/BPLAN(K).GT.0.10) ICASE=MIN0(ICASE,4)
      GO TO (106,106,107,107,100),ICASE
106   JTAU0=K
107   JTAU1=K
100   CONTINUE
      JTAU1=MAX0(JTAU1,3)
*
*  Distribute rays in the core.  This is done in such a way that the
*  mu values at depth JTAU0 are equidistant.  The rationale for this
*  is to get a good representation of the the "core" part of the
*  radiation as a function of mu, in the optically thin regions.
*
      RRK=RR(JTAU0)
      NCORE1=NCORE-1
      IF (NCORE1.EQ.0) PRINT*,' NCORE1= ',NCORE1
      DR=(RR(JTAU0)-SQRT(RR(JTAU0)**2-RR(JTAU1)**2))/NCORE1
      DO 110 I=1,NCORE1
      PIMPAC(I)=SQRT(RR(JTAU0)**2-RRK**2)
      RRK=RRK-DR
110   KIMPAC(I)=JTAU1
C
C RAYS IN ATMOSPHERE
      I=NCORE
      KI=JTAU1
120   KIMPAC(I)=KI
      PIMPAC(I)=RR(KI)
      KIP=KI
121   KI=KI-1
      IF (TAU(KI).LE.1.E-20) write(6,*)' KI,TAU(KI) = ',KI,TAU(KI)
      IF (TAU(KIP)/TAU(KI).LT.TLIM.AND.KI.GT.1) GO TO 121
      I=I+1
      IF (KI.GE.3) GO TO 120
      KIMPAC(I)=0
      PIMPAC(I)=RR(1)
      NIMPAC=I-1
      IF (NIMPAC.LT.NRAYS) GO TO 131
      PRINT 122,NIMPAC,NRAYS,NCORE,TLIM
122   FORMAT(' ** NMBR OF RAYS (NIMPAC) TOO LARGE =',I3
     &  ,' parameter NRAYS= ',I3,'  NCORE,TLIM =',
     & I3,F6.3)
      STOP ' stop in trrays_sph at 122 '
C
C FIND THE NUMBER OF MU-PNTS FOR EACH K, PLUS ONE EXTRA FOR MU=0.0
131   II=NIMPAC+1
      DO 130 K=1,JTAU1
      MMU(K)=II
      XMU(II,K)=0.0
      IF (K+1.EQ.KIMPAC(II-1)) II=II-1
130   CONTINUE
      RETURN
      END
C
C
C_ugj950523:  here the sphereical part ends@@@@@
C
      SUBROUTINE VAAGL(NLB,XL,W)
      implicit real*16 (a-h,o-z)
C
C        DENNA RUTIN BERAEKNAR VAAGLAENGDSPUNKTER OCH MOTSVARANDE
C        VIKTER.  GLAMD(I),I=1,JLBDS  AER DISKONTINUITETER ELLER DEL-
C        NINGSPUNKTER I VAAGLAENGDSLED, INKLUSIVE VAAGLAENGDSSKALANS
C        AENDPUNKTER.  MLD(I), I=1,JLBDS-1   AER DET OENSKADE ANTALET
C        VAAGLAENGDSPUNKTER I RESP. INTERVALL.   DESSA STORHETER
C        INLAESES I RUTINEN.
C        VI ANVAENDER SUBR.  G A U S I .
C        ****   OBSERVERA. JLBDS FAAR HOEGST VARA 15 MED HAER BRUKADE DI
C        OCH FORMATSATSEN 100 ******
C
      DIMENSION XL(500),W(500),GLAMD(100),MLD(100),XLP(10),WP(10)
      COMMON/UTPUT/IREAD,IWRIT
      COMMON/CLINE3/GLAMD,JLBDS
      COMMON/LDOPAC/ ALES,BLES
C
      READ(IREAD,100)JLBDS,ALESX,BLESX
      IF(JLBDS.GT.1) GOTO 31
C ALLOW ONE POINT STANDARD OPACITY
      READ(IREAD,102)XL(1)
      W(1)=1.
      NLB=1
      RETURN
31    CONTINUE
      IF(JLBDS.GT.10)  GO TO 77
C SET UV OPACITY CONSTANTS
      ALES=ALESX
      BLES=BLESX
   77 CONTINUE
      JP=JLBDS-1
      DO3 K=1,JP
    3 READ(IREAD,102)GLAMD(K),MLD(K)
      READ(IREAD,102)GLAMD(JLBDS)
      I=0
      DO2 K=1,JP
      IF(MLD(K).GT.0)GO TO 21
      JIP=2
      XLP(1)=GLAMD(K)
      WP(1)=(GLAMD(K+1)-GLAMD(K))*0.5
      XLP(2)=GLAMD(K+1)
      WP(2)=WP(1)
      MLD(K)=1
      IF(K.EQ.JP)GO TO 22
      IF(MLD(K+1).EQ.0)GO TO 23
   22 MLD(K)=2
      GO TO 23
   21 CONTINUE
      CALL GAUSI(MLD(K),GLAMD(K),GLAMD(K+1),WP,XLP)
      JIP=MLD(K)
   23 CONTINUE
      DO1 J=1,JIP
      IP=I+J
      XL(IP)=XLP(J)
      W(IP)=WP(J)+W(IP)
    1 CONTINUE
    2 I=I+MLD(K)
      NLB=I
C      WRITE(6,*) ' # cont. points, NLB, selected in VAAGL is: ',NLB
C      WRITE(6,*) ' I,WN(I),XL(I),W(I) = '
C      WRITE(6,*) ' (i.e. the continuums wavenumbers, -lengts selected '
C     *,'in VAAGL.for, and their acumulated gausian integration weigts)'
C      DO 502 I=1,NLB
C      WRITE(6,*) I,1.E8/XL(I),XL(I),W(I)
C502   CONTINUE
  100 FORMAT(I5,5X,2F10.0)
  102 FORMAT(F10.5,I5)
      RETURN
      END
C
      FUNCTION VVMLT(A,B,C)
      implicit real*16 (a-h,o-z)
C
C COMPUTE SQRT(B*(A+.5/B/C-SQRT(.5/B/C*(2.*A+.5/B/C)))) WITH SPECIAL
C CONSIDERATION ON SMALL C-VALUES. 73.12.02  *NORD*
C
      D=.5/B/C
      IF(D.GT.20.*A) GO TO 1
      E=A+D-SQRT(D*(2.*A+D))
      GO TO 2
1     E=.5*MAX(A,0.0D+0)**2/D
2     VVMLT=SQRT(B*E)
      RETURN
      END
C
      DOUBLE PRECISION   FUNCTION X02AAF(X)
      implicit real*16 (a-h,o-z)
C
C  RETURNS THE SMALLEST VALUE X02AAF SUCH THAT 1+X02AAF>1
C
C CYBER ?
C     X02AAF=7.2E-15
C IBM DOUBLE PREC.
C      X02AAF = 2.2205D-16
C VAX DOUBLE PREC.
C     X02AAF=1.E-16
C
C APOLLO/DOMAIN WORKSTATION DN 3000/4000  (TRY CALCULATING IT . . . )
C
      X02AAF = 1.0D-15
C
      RETURN
          END
C
      SUBROUTINE X1MAKE(N,XXI,W,XI)
      implicit real*16 (a-h,o-z)
C
C      implicit real*16 (A-H,O-Z)
      PARAMETER(NDIM=100)
      DIMENSION XI(NDIM),W(NDIM),XXI(NDIM)
C
      I=N
      DO 5 J=1,N-1
        XI(I)=XXI(I-1)
        I=I-1
    5 CONTINUE
C
      F=.5
      XI(1)=XI(2)+F*(XI(2)-XI(3))/(W(3)/W(2)-1.)
C
C
      RETURN
C
C
C
      E  N  D
C
      SUBROUTINE XIINIT(WVAL,K,XIVAL,N,XI,W)
      implicit real*16 (a-h,o-z)
C
C     implicit real*16 (A-H,O-Z)
      PARAMETER (NDIM=100)
      DIMENSION XI(NDIM),W(NDIM)
      DIMENSION A(NDIM),B(NDIM),C(NDIM),ITYPE(NDIM)
C
      D1=(XI(2)-XI(1))/(W(2)-W(1))
      D2=(XI(3)-XI(2))/(W(3)-W(2))
      YPI=D1*D1/D2
      DO 105 I=2,N
        HI=W(I)-W(I-1)
        DI=(XI(I)-XI(I-1))/HI
        IF(YPI/DI.GT.1.) THEN
          ITYPE(I)=0
          C(I)=(1./DI-1./YPI)/DI/HI
          BB=1./YPI-2.*XI(I-1)*C(I)
          A(I)=W(I-1)-XI(I-1)*(BB+C(I)*XI(I-1))
          B(I)=-.5*BB/C(I)
          YPI=1./(BB+2.*C(I)*XI(I))
        ELSE
          ITYPE(I)=1
          C(I)=(DI-YPI)/HI
          B(I)=YPI-2.*W(I-1)*C(I)
          A(I)=XI(I-1)-W(I-1)*(B(I)+W(I-1)*C(I))
          YPI=B(I)+2.*C(I)*W(I)
        END IF
  105 CONTINUE
C
C
      RETURN
C
      ENTRY XIMAKE(WVAL,K,XIVAL,N,XI,W)
C
      IF(WVAL.LE.W(K).OR.K.EQ.N) GOTO 12
C
      DO 10 I=1,N-K-1
        K=K+1
        IF(WVAL.LE.W(K)) GO TO 12
   10 CONTINUE
C
      K=N

   12 IF(ITYPE(K).EQ.1) THEN
      XIVAL=A(K)+WVAL*(B(K)+C(K)*WVAL)
      ELSE
        YYY=B(K)*B(K) - (A(K)-WVAL)/C(K)
        IF(YYY.LT.0.) THEN
C        PRINT 1234,K,B(K),A(K),WVAL,C(K),YYY
C      DO 2727 II7=1,N
C      PRINT 1235,II7,W(II7),XI(II7),A(II7),B(II7),C(II7)
C2727  CONTINUE
C
1234  FORMAT('XIMAKE: K,B(K),A(K),WVAL,C(K),YYY '/I5,1P5E13.5)
1235  FORMAT(I5,F8.5,F10.5,1P3E15.5)
         STOP 'XIMAKE-KRASCH'
      ENDIF
        XIVAL=B(K)-SQRT(B(K)*B(K)-(A(K)-WVAL)/C(K))
      END IF
C
C
      RETURN
C
C
      E  N  D
C
      SUBROUTINE XMETAL(ID,N,W,XI)
      implicit real*16 (a-h,o-z)
C
C      IMPLICIT REAL*8(A-H,O-Z)
C
      include 'parameter.inc'
C
      DIMENSION W(6),XI(6)
      COMMON /ODFAD/ V(NDP,4)
C
      DO 5 I=1,N-2
        XI(I)=V(ID,I)
    5 CONTINUE
C
      FN=.5
      XI(N-1)=XI(N-2)+FN*((XI(N-3)-XI(N-2))*(1.-W(N-2))/
     &                          (W(N-3)-W(N-2)))
C
      RETURN
C
C
C
      E  N  D
C
      SUBROUTINE ZEROF(F,DX,DFDX)
      implicit real*16 (a-h,o-z)
C
C FIND DX=-F/DFDX, TO MAKE F ZERO. IF DX=0 AT ENTRY, THEN DFDX IS A
C START APPROXIMATION, AND IT IS THE FIRST CALL. OTHERWISE USE OLD
C INFO.  780926/NORDLUND.
C
      IF(DX.NE.0.) DFDX=(F-FOLD)/DXOLD
      DX=-F/DFDX
C
C TRY TO AVOID OVERFLOW WHEN A JUMP IN TEMPERATURE IS PRESENT 
      IF (DX.GT.10.) DX=10.
C
      FOLD=F
      DXOLD=DX
      RETURN
      END



      SUBROUTINE OSTABLOOK

      implicit real*16 (a-h,o-z)
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
     &  xmettryck(ndp,maxmet),xiontryck(ndp,maxmet),
     *  partp(ndp,0:maxmol),
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

      real*16 ::
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
276   p6ln(ip6) = 2.3025851*float(ip6)     !ln(p6[dyn/cm2])
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
        F1=LOG10(max(PARTP(IT,IDMOL),1.q-99))
        end if
        F2=LOG10(max(AKAPMOL(IT),1.q-99))
      IF (F1+F2.LE.-150) THEN
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
      ELSE IF (F1+F2.GE.150) THEN
         PRINT*,' OVERFLOW FOR JV,IT,IDMOL = ',JV,IT,IDMOL
         PRINT*,' F1,F2 = ',F1, F2
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
      FUNCTION TG01B(II,ND,N,X,F,D,XX)
      implicit real*16 (a-h,o-z)
c     implicit real*16 (A-H),(O-Z)
      DIMENSION X(ND),F(ND),D(ND)
      COMMON /TG01BA/I1,IN,KK
      DATA I1,IN/1,1/
C ROUTINE TO CALCULATE VALUE FXX OF SPLINE IN POINT XX WHEN N KNOTS
C (XI,FI) WITH DERIVATIVE DI ARE GIVEN.
C II<0 => SEARCH THE WHOLE RANGE. II >= 0 => FUNCTION HAS PREVIOUSLY 
C BEEN ENNTERED WITH A SMALLER VALUE OF XX.
C COMMON VALUES I1 AND IN CONTROLS WHAT TO DO IF XX IS OUTSIDE X INTERVAL.
C
C II NEGATIVE, RESET
      IF(II.LT.0) KK=2  
C   
C CHECK IF OUTSIDE  
      IF(XX.LT.X(1)) GOTO 110   
      IF(XX.GT.X(N)) GOTO 120   
      DO 100 K=KK,N 
      IF(XX.LT.X(K)) GOTO 101   
100   CONTINUE  
      KK=N  
      GOTO 102  
101   KK=K  
C   
C CALCULATE FUNCTION
102   DX=X(KK)-X(KK-1)  
      DF=F(KK)-F(KK-1)  
      P=(XX-X(KK-1))/DX 
      Q=1.-P
      TG01B=Q*F(KK-1)+P*F(KK)+P*Q*  
     & (Q*(D(KK-1)*DX-DF)-P*(D(KK)*DX-DF))
1021  FORMAT(I3,1P5E13.5)
      RETURN
C   
C BEFORE X(1)   
110   TG01B=0.  
      IF(I1.LE.0) RETURN
      TG01B=F(1)
      IF(I1.EQ.1) RETURN
      TG01B=TG01B+(XX-X(1))*D(1)
      IF(I1.EQ.2) RETURN
      DX=X(2)-X(1)  
      D2=2.*(3.*(F(2)-F(1))/DX**2-(2.*D(1)+D(2))/DX)
      TG01B=TG01B+.5*(XX-X(1))**2*D2
      IF(I1.EQ.3) RETURN
      D36=(D(1)+D(2)-2.*(F(2)-F(1))/DX)/DX**2   
      TG01B=TG01B+(XX-X(1))*(XX-X(1))**2*D36
      RETURN
C   
C AFTER X(N)
120   TG01B=0.  
      IF(IN.LE.0) RETURN
      TG01B=F(N)
      IF(IN.EQ.1) RETURN
      TG01B=TG01B+(XX-X(N))*D(N)
      IF(IN.EQ.2) RETURN
      DX=X(N)-X(N-1)
      D2=2.*(-3.*(F(N)-F(N-1))/DX**2+(2.*D(N)+D(N-1))/DX)   
      TG01B=TG01B+.5*(XX-X(N))**2*D2
      IF(IN.EQ.3) RETURN
      D36=(D(N)+D(N-1)-2.*(F(N)-F(N-1))/DX)/DX**2   
      TG01B=TG01B+(XX-X(N))*(XX-X(N))**2*D36
      END
C
C
C
      SUBROUTINE DAYTIM(ADATE,ATIME)
      implicit real*16 (a-h,o-z)
C   This routine returns date and time in 9 character format.
C   The routine is highly machine dependent !
C   This version is for VAX.
C
      CHARACTER*9 ADATE,ATIME
      ADATE=' ' 
      ATIME=' ' 
C_ursa      CALL DATE_AND_TIME(ADATE)  
C_ursa      CALL TIME(ATIME)
C      ADATE(5:5)=CHAR(ICHAR(ADATE(5:5))+32)
C      ADATE(6:6)=CHAR(ICHAR(ADATE(6:6))+32)
      RETURN
      END
C
C
      function second()
      implicit real*16 (a-h,o-z)
c  Added by Bjorn S. Nilsson, nbi, on 19-Feb-1991;  for DEC-stations
c      real times(2)
c      second = etime(times)
      second=0
      return
      end
C
C
C      FUNCTION SECOND(TIME0)
C  This routine returns seconds used in floating format.
C  The routine is highly machine dependent !
C  This version is for VAX.
C
C      DATA ICALL/0/ 
C      IF(ICALL.EQ.0) CALL LIB$INIT_TIMER()  
C      ICALL=1
C      CALL LIB$STAT_TIMER(2,LCSEC)  
C      SECOND=0.01*DFLOAT(LCSEC)  
C      RETURN
C      END
C
      SUBROUTINE TIMEX
      implicit real*4 (a-h,o-z)
C   This routine prints total accumulated time and time spent since
C  last call.
C
      CHARACTER*20 FORM/'(A,F05.2,A,F05.2,A)'/
      REAL(KIND=4), SAVE, DIMENSION(2):: time_last=(/0.,0./), time_0
      REAL(KIND=4), SAVE, DIMENSION(2):: time_now=(/0.,0./)
      ENTRY TIMEX0
C      TIME_LAST=SECOND(DUMTIM)  
      call etime(time_last)
      time_0 = time_last
C      TIME_LAST=SECOND()
      RETURN
      ENTRY TIMEX1
C      TIME_NOW=SECOND(DUMTIM)
C      TIME_NOW=SECOND()
      call etime(time_now)
      L1=5
      IF(TIME_NOW(1).GT.99.) L1=6
      IF(TIME_NOW(1).GT.999.) L1=7
      IF(TIME_NOW(1).GT.9999.) L1=8
      IF(TIME_NOW(1).GT.99999.) L1=9
      IF(TIME_NOW(1).GT.999999.) L1=10
      WRITE(FORM(5:6),'(I2.2)') L1
      DELTA_TIME=TIME_NOW(1)-TIME_LAST(1)
      tot_time = time_now(1) - time_0(1)
      L1=5
      IF(DELTA_TIME.GT.99.) L1=6
      IF(DELTA_TIME.GT.999.) L1=7
      IF(DELTA_TIME.GT.9999.) L1=8
      IF(DELTA_TIME.GT.99999.) L1=9
      IF(DELTA_TIME.GT.999999.) L1=10
      WRITE(FORM(13:14),'(I2.2)') L1
      write(6,FORM) ' Total time spent: ',tot_time,' sec. Added time: ',
     1 DELTA_TIME,' sec.'
      TIME_LAST=TIME_NOW
      RETURN
      END
C
C
      SUBROUTINE TIMEF
      implicit real*4 (a-h,o-z)
C   This routine prints total accumulated time and time spent since
C  last call.
C
      REAL(KIND=4), SAVE, DIMENSION(2):: time_last=(/0.,0./), time_0
      REAL(KIND=4), SAVE, DIMENSION(2):: time_now=(/0.,0./)
      CHARACTER*20 FORM/'(A,F05.2,A,F05.2,A)'/
      ENTRY TIME0
C      TIME_LAST=SECOND(DUMTIM)  
      call etime(time_last)
      time_0 = time_last
C      TIME_LAST=SECOND()
      RETURN
      ENTRY TIME1
C      TIME_NOW=SECOND(DUMTIM)
C      TIME_NOW=SECOND()
      call etime(time_now)
      L1=5
      IF(TIME_NOW(1).GT.99.) L1=6
      IF(TIME_NOW(1).GT.999.) L1=7
      IF(TIME_NOW(1).GT.9999.) L1=8
      IF(TIME_NOW(1).GT.99999.) L1=9
      IF(TIME_NOW(1).GT.999999.) L1=10
      WRITE(FORM(5:6),'(I2.2)') L1
      DELTA_TIME=TIME_NOW(1)-TIME_LAST(1)
      tot_time = time_now(1) - time_0(1)
      L1=5
      IF(DELTA_TIME.GT.99.) L1=6
      IF(DELTA_TIME.GT.999.) L1=7
      IF(DELTA_TIME.GT.9999.) L1=8
      IF(DELTA_TIME.GT.99999.) L1=9
      IF(DELTA_TIME.GT.999999.) L1=10
      WRITE(FORM(13:14),'(I2.2)') L1
      PRINT FORM, ' Total time spent: ',tot_time,' sec. Added time: ',
     1 DELTA_TIME,' sec.'
      TIME_LAST=TIME_NOW
      RETURN
      END
      
!-----------------------------------------------------------------------
! GETTIME: Prints total accumulated time and time spent since last call
! Juncher 2015
!-----------------------------------------------------------------------
      subroutine gettime(io)
      
      implicit real*4 (a-h,o-z)
      include 'parameter.inc'
      integer :: tot_hour, tot_min, tot_sec
      real, save, dimension(2):: time_0, time_last, time_now
      
      if(io .eq. 0) then
        call etime(time_last)
        time_0 = time_last
      else
        call etime(time_now)

        tot_hour = int((time_now(1)-time_0(1))/3600.)
        tot_min  = int((time_now(1)-time_0(1)-tot_hour*3600.)/60.)
        tot_sec  = int(time_now(1)-time_0(1)-tot_hour*3600.-tot_min*60.)
        
        time_last = time_now
        
        if(io .eq. 1) then
          write(*,'(a19,i2,2(a1,i2.2))')' Total time spent: ',
     *      tot_hour, ':', tot_min, ':', tot_sec
        else
          write(6,'(a19,i2,2(a1,i2.2))')' Total time spent: ',
     *      tot_hour, ':', tot_min, ':', tot_sec
        end if
      end if

      return
      end

C
C

      subroutine interp( ni, xx, yy, fi, fxi, fyi, fxyi )
      implicit real*16 (a-h,o-z)
C
       include 'parameter.inc'
c***********************************************************************
c     interpolate pressure or energy density using monotonized         *
c     bicubic hermite interpolation                                    *
C     The routine assumes that there exist a table with values, f(x,y),*
C     of the (natural) logarithm of the function F(x,y) in which we    *
C     we wish to interpolate.                                          *
C     *** do the routine need input values of fx, fy, fxy at the grid **
C     *** points ????????????    ***************************************
C     On calling the routine assumes that the common block /eostab/    *
C     contains the minimum values xmin,ymin of x and y and the         *
C     increments, dx,dy, and the total number of x and y values, nx,ny *
C     The function can now return the values of F at values (xi,yi) in *
C     ni points simoultanously. ni < ndp.                             *
C     The dimension, mtemp,mpe, of the table must be set in the parameter  *
C     statement. This is a pure dimension problem and it must just     *
C     assure that the dimension of the table covers the real table     *
C     The value of ndp, on the other hand, must be greater than       *
C     max(nx,ny) to avoid out-of-range problems.                       *
c***********************************************************************
C ifix changet to idint because our VAX machine doesn't reconise idfix
C UGJ/900417
c
C      parameter ( ndp = 100, mtemp = 6, mpe = 1 )
c

      common/eostab/ xmin, ymin, dx, dy, f (mtemp,mpe)
     *     , fx (mtemp,mpe), fy(mtemp,mpe), fxy(mtemp,mpe), nx, ny
c
      dimension x (mtemp), y (mpe)
c
      dimension xx(ni), yy (ni)
      dimension fi(ni), fxi(ni), fyi(ni), fxyi(ni)
c
      dimension ii (ndp), jj  (ndp)
      dimension xi (ndp), yi  (ndp)
c
      dimension f00(ndp), fx00(ndp), fy00(ndp), fxy00(ndp),
     .          f01(ndp), fx01(ndp), fy01(ndp), fxy01(ndp),
     .          f10(ndp), fx10(ndp), fy10(ndp), fxy10(ndp),
     .          f11(ndp), fx11(ndp), fy11(ndp), fxy11(ndp)
c
      dimension h1 (ndp), h2  (ndp), h3  (ndp), h4   (ndp),
     .          hx1(ndp), hx2 (ndp), hx3 (ndp), hx4  (ndp)
      dimension g1 (ndp), g2  (ndp), g3  (ndp), g4   (ndp),
     .          gy1(ndp), gy2 (ndp), gy3 (ndp), gy4  (ndp)
c
c-----------------------------------------------------------------------
c     construct coordinate mesh
c-----------------------------------------------------------------------
c
      do 1 i = 1, nx
      x(i) = xmin + dfloat( i - 1 ) * dx
    1 continue
c
      do 2 i = 1, ny
      y(i) = ymin + dfloat( i - 1 ) * dy
    2 continue
c
c-----------------------------------------------------------------------
c     evaluate interpolant 
c-----------------------------------------------------------------------
c
      do 3 i = 1, ni
c
c     forbid points off edges of table
      xx(i) = max( x(1), min( xx(i), x(nx) ) )
      yy(i) = max( y(1), min( yy(i), y(ny) ) )
c
c     compute indices of lower left corner of the cell
      ii(i) = min( int( (xx(i) - x(1))/dx + 1.01 ), nx - 1 )
C  at VAX:  idint if *8
      jj(i) = min( int( (yy(i) - y(1))/dy + 1.01 ), ny - 1 )
c
c     compute coordinates ( xi, yi ) relative to lower left corner of cell
c     in units of cell dimensions 
      xi(i) = ( xx(i) - x( ii(i) ) ) / dx
      yi(i) = ( yy(i) - y( jj(i) ) ) / dy
c
c     evaluate basis functions at (xi, yi)
      h 2(i) = - ( 2.0*xi(i) - 3.0 ) * xi(i)**2
      g 2(i) = - ( 2.0*yi(i) - 3.0 ) * yi(i)**2
      h 1(i) = 1.0 - h2(i)
      g 1(i) = 1.0 - g2(i)
      h 3(i) = dx * ( ( xi(i) - 2.0 ) * xi(i) + 1.0 ) * xi(i)
      g 3(i) = dy * ( ( yi(i) - 2.0 ) * yi(i) + 1.0 ) * yi(i)
      h 4(i) = dx * ( xi(i) - 1.0 ) * xi(i)**2
      g 4(i) = dy * ( yi(i) - 1.0 ) * yi(i)**2
      hx2(i) = - 6.0 * ( xi(i) - 1.0 ) * xi(i) / dx
      gy2(i) = - 6.0 * ( yi(i) - 1.0 ) * yi(i) / dy
      hx1(i) = - hx2(i)
      gy1(i) = - gy2(i)
      hx3(i) = ( 3.0*xi(i) - 4.0 ) * xi(i) + 1.0
      gy3(i) = ( 3.0*yi(i) - 4.0 ) * yi(i) + 1.0
      hx4(i) = ( 3.0*xi(i) - 2.0 ) * xi(i)
      gy4(i) = ( 3.0*yi(i) - 2.0 ) * yi(i)
c
c     assemble f, fx, fy, fxy at corners of each cell
      f  00(i) = f  ( ii(i)    , jj(i)     )
      fx 00(i) = fx ( ii(i)    , jj(i)     )
      fy 00(i) = fy ( ii(i)    , jj(i)     )
      fxy00(i) = fxy( ii(i)    , jj(i)     )
      f  01(i) = f  ( ii(i)    , jj(i) + 1 )
      fx 01(i) = fx ( ii(i)    , jj(i) + 1 )
      fy 01(i) = fy ( ii(i)    , jj(i) + 1 )
      fxy01(i) = fxy( ii(i)    , jj(i) + 1 )
      f  10(i) = f  ( ii(i) + 1, jj(i)     )
      fx 10(i) = fx ( ii(i) + 1, jj(i)     )
      fy 10(i) = fy ( ii(i) + 1, jj(i)     )
      fxy10(i) = fxy( ii(i) + 1, jj(i)     )
      f  11(i) = f  ( ii(i) + 1, jj(i) + 1 )
      fx 11(i) = fx ( ii(i) + 1, jj(i) + 1 )
      fy 11(i) = fy ( ii(i) + 1, jj(i) + 1 )
      fxy11(i) = fxy( ii(i) + 1, jj(i) + 1 )
    3 continue
c
      do 4 i = 1, ni
c
c     compute f, fx, fy, and fxy at (xi, yi)
      fi  (i) = g 1(i) * ( h 1(i)*f  00(i) + h 2(i)*f  10(i)
     .                   + h 3(i)*fx 00(i) + h 4(i)*fx 10(i) )
     .        + g 2(i) * ( h 1(i)*f  01(i) + h 2(i)*f  11(i)
     .                   + h 3(i)*fx 01(i) + h 4(i)*fx 11(i) )
     .        + g 3(i) * ( h 1(i)*fy 00(i) + h 2(i)*fy 10(i)
     .                   + h 3(i)*fxy00(i) + h 4(i)*fxy10(i) )
     .        + g 4(i) * ( h 1(i)*fy 01(i) + h 2(i)*fy 11(i)
     .                   + h 3(i)*fxy01(i) + h 4(i)*fxy11(i) )
      fxi (i) = g 1(i) * ( hx1(i)*f  00(i) + hx2(i)*f  10(i)
     .                   + hx3(i)*fx 00(i) + hx4(i)*fx 10(i) )
     .        + g 2(i) * ( hx1(i)*f  01(i) + hx2(i)*f  11(i)
     .                   + hx3(i)*fx 01(i) + hx4(i)*fx 11(i) )
     .        + g 3(i) * ( hx1(i)*fy 00(i) + hx2(i)*fy 10(i)
     .                   + hx3(i)*fxy00(i) + hx4(i)*fxy10(i) )
     .        + g 4(i) * ( hx1(i)*fy 01(i) + hx2(i)*fy 11(i)
     .                   + hx3(i)*fxy01(i) + hx4(i)*fxy11(i) )
      fyi (i) = gy1(i) * ( h 1(i)*f  00(i) + h 2(i)*f  10(i)
     .                   + h 3(i)*fx 00(i) + h 4(i)*fx 10(i) )
     .        + gy2(i) * ( h 1(i)*f  01(i) + h 2(i)*f  11(i)
     .                   + h 3(i)*fx 01(i) + h 4(i)*fx 11(i) )
     .        + gy3(i) * ( h 1(i)*fy 00(i) + h 2(i)*fy 10(i)
     .                   + h 3(i)*fxy00(i) + h 4(i)*fxy10(i) )
     .        + gy4(i) * ( h 1(i)*fy 01(i) + h 2(i)*fy 11(i)
     .                   + h 3(i)*fxy01(i) + h 4(i)*fxy11(i) )
      fxyi(i) = gy1(i) * ( hx1(i)*f  00(i) + hx2(i)*f  10(i)
     .                   + hx3(i)*fx 00(i) + hx4(i)*fx 10(i) )
     .        + gy2(i) * ( hx1(i)*f  01(i) + hx2(i)*f  11(i)
     .                   + hx3(i)*fx 01(i) + hx4(i)*fx 11(i) )
     .        + gy3(i) * ( hx1(i)*fy 00(i) + hx2(i)*fy 10(i)
     .                   + hx3(i)*fxy00(i) + hx4(i)*fxy10(i) )
     .        + gy4(i) * ( hx1(i)*fy 01(i) + hx2(i)*fy 11(i)
     .                   + hx3(i)*fxy01(i) + hx4(i)*fxy11(i) )
c
c     convert ln f to f
      fi  (i) = exp( fi(i) )
c
    4 continue
c
      return
      end

C
      SUBROUTINE SCALEMOD
      implicit real*16 (a-h,o-z)
C
C This routine scales all parameters in one model to a model with another
C number of depth points. The (new) tau-scale from the input file is in
C common TAUC, and the variables to be scaled are those in common STATEC
C
      include 'parameter.inc'
C
      DIMENSION W(3,NDP),SX(NDP),SXDER(NDP)
C
      COMMON /STATEC/PR(NDP),PT(NDP),P(NDP),G(NDP),Z(NDP),D(NDP),
     &V(NDP),FC(NDP),PE(NDP),T(NDP),TAULN(NDP),ROSTAT(NDP),NTAU,ITER
      COMMON /TAUC/TAU(NDP),DTAULN(NDP),JTAU
      COMMON /TG01BA/I1,IN,KK
      COMMON /CMOLRAT/ FOLD(NDP,8),MOLOLD,KL
C
C      DATA I1,IN/2,2/ !0~F==0. outside interpolation interval
C                      1~F==value in first (last) point -||-
C                      2~F from linear extrapolation -||-
C ---------------- ***************** ----------------------------  
C CALCULATE log-log SPLINE DERIVATIVES of variables in STATEC
C
C SxDER IS THE DERIVATIVE IN ALL NTAU NODES WHEN
C TAULN AND SxK ARE GIVEN 
C Calculate interpolated values at depths of new tauscale. Save
C in STATEC and TAUC

      WRITE(6,*) 'NTAU in SCALEMOD = ', NTAU
      if (ntau.gt.ndp) WRITE(6,*) 'NDP = ', NDP
      if (ntau.gt.ndp) STOP  ' increase dimension ndp '
      I1 = 2
      IN = 2
      DO 303 K=1,NTAU
303   write(6,301) K,EXP(TAULN(K)),P(K),G(K),Z(K),PE(K),T(K)
          DO 1 K=1,NTAU
1         Sx(K) = LOG(P(K))
          CALL TB04A(NTAU,tauln,SX,SXDER,W)
          DO 11 K=1,JTAU
          X = LOG(TAU(K))
11        P(K)=EXP(TG01B(-1,NDP,NTAU,TAULN,SX,SXDER,X))
C
C          DO 2 K=1,NTAU
C2         Sx(K) = LOG(G(K))
C          CALL TB04A(NTAU,tauln,SX,SXDER,W)
C          DO 12 K=1,JTAU
C          X = LOG(TAU(K))
C12        G(K)=EXP(TG01B(-1,NDP,NTAU,TAULN,SX,SXDER,X))
C
C          DO 3 K=1,NTAU
C3         Sx(K) = LOG(Z(K))
C          CALL TB04A(NTAU,tauln,SX,SXDER,W)
C          DO 13 K=1,JTAU
C          X = LOG(TAU(K))
C13        Z(K)=EXP(TG01B(-1,NDP,NTAU,TAULN,SX,SXDER,X))
C
          DO 4 K=1,NTAU
4         Sx(K) = LOG(PE(K))
          CALL TB04A(NTAU,tauln,SX,SXDER,W)
          DO 14 K=1,JTAU
          X = LOG(TAU(K))
14        PE(K)=EXP(TG01B(-1,NDP,NTAU,TAULN,SX,SXDER,X))
C
          DO 5 K=1,NTAU
5         Sx(K) = LOG(T(K))
          CALL TB04A(NTAU,tauln,SX,SXDER,W)
          DO 15 K=1,JTAU
          X = LOG(TAU(K))
15        T(K)=EXP(TG01B(-1,NDP,NTAU,TAULN,SX,SXDER,X))
C
C
      IF (MOLOLD.EQ.1) THEN
         DO 26 M=1,8
          DO 6 K=1,NTAU
6         Sx(K) = LOG(FOLD(K,M))
          CALL TB04A(NTAU,tauln,SX,SXDER,W)
          DO 16 K=1,JTAU
          X = LOG(TAU(K))
16        FOLD(K,M)=EXP(TG01B(-1,NDP,NTAU,TAULN,SX,SXDER,X))
26       CONTINUE
      END IF
C
      NTAU=JTAU
      DO 20 K=1,NTAU
20    G(K)=G(1)
C20    TAULN(K)=LOG(TAU(K))
C
      WRITE(6,*) ' NTAU in SCALEMOD = ', NTAU
      WRITE(6,*) ' K    TAU    P-tot    G      Z      Pe      T'
      DO 300 K=1,NTAU
300   write(6,301) K,TAU(K),P(K),G(K),Z(K),PE(K),T(K)
301   FORMAT(I3,1P5E11.3,0PF8.0)
C
      RETURN
      END
c
c
      SUBROUTINE TB04A(N,X,F,D,W)
      implicit real*16 (a-h,o-z)
C
      DIMENSION X(N),F(N),D(N),W(3,N)   
C   
C THIS VERSION OF TB04A IS COMPATIBLE WITH RKU*HARWELL.TB04A, BUT   
C FOR UNKNOWN REASONS 30 % FASTER.  770815
C INPUT ARE X(I), FUNCTION VALUE IN N KNOTS. THEN THE DERIVATIVE IS 
C CALCULATED IN THE KNOTS. W=0 AT SUCCESFULL RETURN, OTHERWISE W=1.
C X(I) SHOULD BE IN STRICTH INCREASING ORDER, X1<X2<...<XN.
C         
C FIRST POINT   
      CXB=1./(X(2)-X(1))
      CXC=1./(X(3)-X(2))
      DFB=F(2)-F(1) 
      DFC=F(3)-F(2) 
      W(1,1)=CXB*CXB
      W(3,1)=-CXC*CXC   
      W(2,1)=W(1,1)+W(3,1)  
      D(1)=2.*(DFB*CXB*CXB*CXB-DFC*CXC*CXC*CXC) 
C   
C INTERIOR POINTS   
      N1=N-1
      DO 100 K=2,N1 
      CXA=CXB   
      CXB=1./(X(K+1)-X(K))  
      DFA=DFB   
      DFB=F(K+1)-F(K)   
      W(1,K)=CXA
      W(3,K)=CXB
      W(2,K)=2.*(CXA+CXB)   
      D(K)=3.*(DFB*CXB*CXB+DFA*CXA*CXA) 
100   CONTINUE  
C   
C LAST POINT
      W(1,N)=CXA*CXA
      W(3,N)=-CXB*CXB   
      W(2,N)=W(1,N)+W(3,N)  
      D(N)=2.*(DFA*CXA*CXA*CXA-DFB*CXB*CXB*CXB) 
C   
C ELIMINATE AT FIRST POINT  
      C=-W(3,1)/W(3,2)  
      W(1,1)=W(1,1)+C*W(1,2)
      W(2,1)=W(2,1)+C*W(2,2)
      D(1)=D(1)+C*D(2)  
      W(3,1)=W(2,1) 
      W(2,1)=W(1,1) 
C   
C ELIMINATE AT LAST POINT   
      C=-W(1,N)/W(1,N-1)
      W(2,N)=W(2,N)+C*W(2,N-1)  
      W(3,N)=W(3,N)+C*W(3,N-1)  
      D(N)=D(N)+C*D(N-1)
      W(1,N)=W(2,N) 
      W(2,N)=W(3,N) 
C   
C ELIMINATE SUBDIAGONAL 
      DO 110 K=2,N  
      C=-W(1,K)/W(2,K-1)
      W(2,K)=W(2,K)+C*W(3,K-1)  
      D(K)=D(K)+C*D(K-1)
110   CONTINUE  
C   
C BACKSUBSTITUTE
      D(N)=D(N)/W(2,N)  
      DO 120 KK=2,N 
      K=(N+1)-KK
      D(K)=(D(K)-W(3,K)*D(K+1))/W(2,K)  
120   CONTINUE  
C   
      RETURN
      END

      subroutine tioeq(NTAU,T,PE)
      implicit real*16 (a-h,o-z)
*
* This subroutine calculates the balance between Ti I, Ti II and TiO. It
* assumes that all titanium exists in one of these forms.
* Method:  Calculate  q1 = n(Ti II)/n(Ti I) and q2 = n(TiO)/n(Ti I).
*          Then  n(Ti I)*(1 + q1 + q2) = n(Ti(tot))
*          n(Ti(tot)) = n(Ti)/n(H) * n(H)/density * density
*                     = abund(Ti)  * 1/(xmy*xmh)  * density
*          Finally, e.g.  p(TiO)=kT*n(Ti I)*q2 etc.
* 
* Output:  ptio (array of partial pressures (cgs units) for TiO at all depths)    
*
* PS: dissoc. Konstant aer updaterat.
*****
      include 'parameter.inc'
      dimension qr(22),tqr(22),T(NDP),PE(NDP)
*
      common/ci5/abmarcs(17,ndp),anjon(17,5),h(5),part(17,5),
     *dxi,f1,f2,f3,f4,f5,xkhm,xmh,xmy(ndp)
      common /ketio/ptio(ndp),ro(ndp),poxg1(ndp)
C      real    absc,abti,abv,abmn,abco                                   /auxabund/
      common /auxabund/ absc,abti,abv,abmn,abco !                        /auxabund/
      common /ti1ti2/ xnti1(ndp),xnti2(ndp)
*
      data tqr/1000., 1500., 2000., 2500., 3000., 3500., 4000., 4500.,
     &         5000., 5500., 6000., 6500., 7000., 8000., 9000., 10000.,
     &        12000., 14000., 16000., 18000., 20000., 50000./
      data qr/1.708, 1.911, 2.035, 2.097, 2.112, 2.091, 2.042, 1.972,
     &        1.887, 1.789, 1.684, 1.574, 1.462, 1.241, 1.032, 0.884,
     &        0.541, 0.330, 0.194, 0.111, 0.062, 0.001/
*
* init (abundance of Ti)
*
      abti = 4.99
*
* epsti is the ratio  n(Ti)/n(H)
*
ccc      epsti=10.**(abti-12)
      epsti=abmarcs(17,1)
      write(6,668) abmarcs(17,1),abmarcs(1,1),abti,10.**(abti-12)
668   format(' Tioqe: abund(17),(1),abti,10.**(abti-12):',1p4e10.3)
*
* div is the ratio  n(H)/rho   (cf. subr. injon)
*
      div = 1./(xmy(1)*xmh)
*
      write(6,*)' TIOEQ:i,t,q1,q2,log(ptio),poxg1'
     &     ,',tiokp,xnti1,xnti2 '

      do 10 i=1,ntau      
        thta=5040./t(i)
*
* find the ratio  u(Ti II)/u(Ti I)  by interpolating in array qr
*
        j=1
        do 2 jj=1,22
          if(t(i).le.tqr(jj)) goto 3
          j=jj
    2   continue
    3   uratio=qr(j) + (qr(j+1)-qr(j))/(tqr(j+1)-tqr(j))*(t(i)-tqr(j))
*
* q1 is the ratio  n(Ti II)/n(Ti I)
*
        q1=.6667*t(i)**2.5*uratio*exp(-6.82*11605./t(i))/pe(i)
*
* tiokp is the equilibrium constant for TiO
* the constant -8.1756 correspond to a dissociation energy of 6.87eV
*   ref: Huber&Herzberg 1979)
*
        tiokp=10.**( 13.398 + thta*(-8.1756 + thta*(.40873 +
     &             thta*(-.057937 + thta*.0029287)))  )
*
* q2 is the ratio  n(TiO)/n(Ti I)
*
        q2=poxg1(i)/tiokp
*
* xnti1 is n(Ti I)    (i.e. per cm**3)
*
        xnti1(i)=epsti*ro(i)*div/(1.+q1+q2)
        ptio(i)=1.38e-16*t(i)*xnti1(i)*q2
        xnti2(i)=xnti1(i)*q1
C
        if (i/10*10.eq.i .or.i.eq.1) then
         write(6,662) i,t(i),q1,q2
     &        ,log10( max(ptio(i),1.d-99) )
     &        ,poxg1(i),tiokp,xnti1(i),xnti2(i)
662     format( i3,f7.0,1p2e9.2,0pf7.2,1p4e9.2 )
C          write(6,*) ' i',i,' t',t(i),' pe',pe(i),
C     &            ' pg',pg,' q1',q1,' q2',q2,' ptio',
C     &             ptio(i),' poxg1(i)',poxg1(i)
C          write(6,*) ' uratio',uratio
        end if
   10 continue
*
*
      return
      end



      FUNCTION TG01BT(II,ND,N,X,F,D,XX) 
      implicit real*16 (a-h,o-z)
      DIMENSION X(ND),F(ND),D(ND)
      COMMON /TG01BA/I1,IN,KK
      I1 = 1     ! we fix extrapolated values to be == end value
      IN = 1
C ROUTINE TO CALCULATE VALUE FXX OF SPLINE IN POINT XX WHEN N KNOTS
C (XI,FI) WITH DERIVATIVE DI ARE GIVEN.
C II<0 => SEARCH THE WHOLE RANGE. II >= 0 => FUNCTION HAS PREVIOUSLY 
C BEEN ENTERED WITH A SMALLER VALUE OF XX.
C COMMON VALUES I1 AND IN CONTROLS WHAT TO DO IF XX IS OUTSIDE X INTERVAL.
C
C II NEGATIVE, RESET
      IF(II.LT.0) KK=2  
C   
C CHECK IF OUTSIDE  
      IF(XX.LT.X(1)) GOTO 110   
      IF(XX.GT.X(N)) GOTO 120   
      DO 100 K=KK,N 
      IF(XX.LT.X(K)) GOTO 101   
100   CONTINUE  
      KK=N  
      GOTO 102  
101   KK=K  
C   
C CALCULATE FUNCTION
102   DX=X(KK)-X(KK-1)  
      DF=F(KK)-F(KK-1)  
      P=(XX-X(KK-1))/DX 
      Q=1.-P
      TG01BT=Q*F(KK-1)+P*F(KK)+P*Q*  
     & (Q*(D(KK-1)*DX-DF)-P*(D(KK)*DX-DF))  
      RETURN
C   
C BEFORE X(1)   
110   TG01BT=0.  
      IF(I1.LE.0) RETURN
      TG01BT=F(1)
      IF(I1.EQ.1) RETURN
      TG01BT=TG01BT+(XX-X(1))*D(1)
      IF(I1.EQ.2) RETURN
      DX=X(2)-X(1)  
      D2=2.*(3.*(F(2)-F(1))/DX**2-(2.*D(1)+D(2))/DX)
      TG01BT=TG01BT+.5*(XX-X(1))**2*D2
      IF(I1.EQ.3) RETURN
      D36=(D(1)+D(2)-2.*(F(2)-F(1))/DX)/DX**2   
      TG01BT=TG01BT+(XX-X(1))*(XX-X(1))**2*D36
      RETURN
C   
C AFTER X(N)
120   TG01BT=0.  
      IF(IN.LE.0) RETURN
      TG01BT=F(N)
      IF(IN.EQ.1) RETURN
      TG01BT=TG01BT+(XX-X(N))*D(N)
      IF(IN.EQ.2) RETURN
      DX=X(N)-X(N-1)
      D2=2.*(-3.*(F(N)-F(N-1))/DX**2+(2.*D(N)+D(N-1))/DX)   
      TG01BT=TG01BT+.5*(XX-X(N))**2*D2
      IF(IN.EQ.3) RETURN
      D36=(D(N)+D(N-1)-2.*(F(N)-F(N-1))/DX)/DX**2   
      TG01BT=TG01BT+(XX-X(N))*(XX-X(N))**2*D36
      END
C
C
C
      SUBROUTINE TB04AT(ND,N,X,F,D,W)
      implicit real*16 (a-h,o-z)
C
      DIMENSION X(ND),F(ND),D(ND),W(3,ND)   
C   
C THIS VERSION OF TB04A IS identical to TB04A, except that it can be
C called with arrays which might be dimensioned larger than the part
C of it which is used for the interpolation.
C INPUT ARE X(I), FUNCTION VALUE IN N KNOTS. THEN THE DERIVATIVE IS 
C CALCULATED IN THE KNOTS. W=0 AT SUCCESFULL RETURN, OTHERWISE W=1.
C X(I) SHOULD BE IN STRICTH INCREASING ORDER, X1<X2<...<XN.
C In the call from SUBROUTINE PROFILE, X(1), X(2),...,X(NTAU) are the 
C temperature values at the NTAU optical depth values, and the derivative
C dPg/dT (=D(N))is calculated in the NTAU points, for later use to 
C calculate Pg in the points where the temperature of the
C absorption coefficient is known.  
C         
C FIRST POINT   
      CXB=1./(X(2)-X(1))
      CXC=1./(X(3)-X(2))
      DFB=F(2)-F(1) 
      DFC=F(3)-F(2) 
      W(1,1)=CXB*CXB
      W(3,1)=-CXC*CXC   
      W(2,1)=W(1,1)+W(3,1)  
      D(1)=2.*(DFB*CXB*CXB*CXB-DFC*CXC*CXC*CXC) 
C   
C INTERIOR POINTS   
      N1=N-1
      DO 100 K=2,N1 
      CXA=CXB   
      CXB=1./(X(K+1)-X(K))  
      DFA=DFB   
      DFB=F(K+1)-F(K)   
      W(1,K)=CXA
      W(3,K)=CXB
      W(2,K)=2.*(CXA+CXB)   
      D(K)=3.*(DFB*CXB*CXB+DFA*CXA*CXA) 
100   CONTINUE  
C   
C LAST POINT
      W(1,N)=CXA*CXA
      W(3,N)=-CXB*CXB   
      W(2,N)=W(1,N)+W(3,N)  
      D(N)=2.*(DFA*CXA*CXA*CXA-DFB*CXB*CXB*CXB) 
C   
C ELIMINATE AT FIRST POINT  
      C=-W(3,1)/W(3,2)  
      W(1,1)=W(1,1)+C*W(1,2)
      W(2,1)=W(2,1)+C*W(2,2)
      D(1)=D(1)+C*D(2)  
      W(3,1)=W(2,1) 
      W(2,1)=W(1,1) 
C   
C ELIMINATE AT LAST POINT   
      C=-W(1,N)/W(1,N-1)
      W(2,N)=W(2,N)+C*W(2,N-1)  
      W(3,N)=W(3,N)+C*W(3,N-1)  
      D(N)=D(N)+C*D(N-1)
      W(1,N)=W(2,N) 
      W(2,N)=W(3,N) 
C   
C ELIMINATE SUBDIAGONAL 
      DO 110 K=2,N  
      C=-W(1,K)/W(2,K-1)
      W(2,K)=W(2,K)+C*W(3,K-1)  
      D(K)=D(K)+C*D(K-1)
110   CONTINUE  
C   
C BACKSUBSTITUTE
      D(N)=D(N)/W(2,N)  
      DO 120 KK=2,N 
      K=(N+1)-KK
      D(K)=(D(K)-W(3,K)*D(K+1))/W(2,K)  
120   CONTINUE  
C   
      RETURN
      END
C

      subroutine osinit
C     program intp_atom
      implicit real*16 (a-h,o-z)
      include 'parameter.inc'
C
      CHARACTER MOLNAME*4,OSFIL*60,SAMPLING*3
      COMMON/COS/WNOS(NWL),CONOS(NDP,NWL),WLOS(NWL),WLSTEP(NWL)
     *    ,KOS_STEP,NWTOT,NOSMOL,NEWOSATOM,NEWOSATOMLIST
     *    ,nchrom,OSFIL(30),MOLNAME(30),SAMPLING
      common/catoms/opjva(nwl,9,3:17),fxln(nwl,9,3:17)
C 
C  program to read the atomic opacities of decreasing wavenumbers
C  constructed by BE, and create an opacity file
C  at the frequency points of interest at increasing wavenumbers.
      
      dimension opjvit(nwl)
     &   ,opl(9,17),opl_pr(9,17),p6ln(9),dadp(9),wp(3,9)
     &   ,sumatdir(9,17)
      common /catoms_head/ p6(9),t_at(17)
c     character osfilatom*39
      data np6,nftemp /9,17/
c     data osfilatom/'/ste1/uffegj/atoms/metals_sun_ascii.x03'/

      write(6,*) ' We did a resolution based set of os wavenumbers'
C      write(6,245) nwtot,osresl/dfloat(kos_step)
C     & ,kos_step*step,wnos(1),wnos(nwtot)
C     & ,1.e4/wnos(nwtot),1.e4/wnos(1)
C245   format(' Total ',i6,' OS wavenumbers for Marcs radiative transf.'
C     & ,/' Approximate resolution is',f8.0,' (corresponding to ',
C     & 'step factor',f8.5,')'
C     & ,/' OS interval:',f7.1,'-',f9.1,' cm^-1 (=',f5.2,'-',f5.1,'mu)')

      DO 2163 I=1,NOSMOL
      IF (MOLNAME(i).eq.'ATOM') THEN
        open(41,file=osfil(i),form='formatted',recl=1246,
     &     status='unknown',readonly)
        write(6,2164) osfil(i)
2164    format(' We opened file',a45,' for reading ATOMic data')
        go to 2165
      END IF
2163  continue
2165  continue

       call atoms_head
       natom = 0
       nerr = 0
       jv = nwtot
       nwatom_os = 0

       read(41,1060)
     &   waven,((opl(ip6,itemp),ip6=1,np6),itemp=1,nftemp)
1060   format(d23.16,153f8.3)
       do 276 ip6 = 1,np6
       p6ln(ip6) = 2.3025851*dfloat(ip6)     !ln(p6[dyn/cm2])
       do 276 it = 1,nftemp
       sumatdir(ip6,it) = 0.                 !init sum op[cm/g*]
276    continue
       write(6,1230)waven,1.e4/waven
1230   format(' first wavenumber in atomic OS:',f10.1
     &           ,'(=',f9.5,' mu)')
       write(6,*) ' nwtot = ',nwtot
       write(6,*) ' wnos(1),wnos(nwtot)=', wnos(1),wnos(nwtot)

       write(6,*) ' In atom_init: jv,waven,wnos(jv),waven_pr'
     &       ,' ((opjva(jv,ip6,itemp),ip6=1,np6),itemp=3,nftemp)'

       do 1065 ia = 1,160000
         waven_pr = waven
         do 270 itemp=1,nftemp
         do 270 ip6=1,np6
270      opl_pr(ip6,itemp)= opl(ip6,itemp)
       read(41,1060,end=1069,err=279)
     &   waven,((opl(ip6,itemp),ip6=1,np6),itemp=1,nftemp)
       dw = waven_pr - waven
       do 2761 ip6 = 1,np6
       do 2761 it = 1,nftemp
       sumatdir(ip6,it) = sumatdir(ip6,it) + dw * 10.**opl(ip6,it) !sum op[cm/g*]
2761   continue
       if (waven.lt.wnos(1)   .or.  waven.gt.wnos(nwtot)) go to 278
       natom = natom + 1
       if (natom.eq.1) waven_1 = waven
       waven_2 = waven

290    continue
       if (waven.le.wnos(jv) .and. waven_pr.gt.wnos(jv)) then
          if (wnos(jv)-waven .le. waven_pr-wnos(jv)) then
              do 281 itemp = 3,nftemp
              do 280 ip6=1,np6
              opjvit(ip6) = 2.3025851*opl(ip6,itemp) !log_10(listed) -> ln(for intp)
              opjva(jv,ip6,itemp) = opjvit(ip6)
280           continue
C  The derivative, DADP(ip6)=d(ln(opac))/d(p6ln), at each p6 and temp.
              CALL TB04AT(np6,np6,p6ln,opjvit,DADP,WP)
              do 282 ip6=1,np6
282           fxln(jv,ip6,itemp) = dadp(ip6)   !d(ln(opac))/d(ln(p6)) for t=temp
281           continue
          else
              do 285 itemp = 3,nftemp
              do 284 ip6=1,np6
              opjvit(ip6) = 2.3025851*opl_pr(ip6,itemp)
284           opjva(jv,ip6,itemp) = opjvit(ip6)
              CALL TB04AT(np6,np6,p6ln,opjvit,DADP,WP)
              do 286 ip6=1,np6
286           fxln(jv,ip6,itemp) = dadp(ip6)
285           continue
          end if
          nwatom_os = nwatom_os + 1
          if(nwatom_os.eq.1 .or. nwatom_os/4000*4000.eq.nwatom_os) 
     &       write(6,1226) jv,waven,wnos(jv),waven_pr
     &       ,((opjva(jv,ip6,itemp),ip6=1,np6,2),itemp=3,nftemp,13)
1226      format(i5,3f10.2/,1p,15(5e10.3/))
          jv = jv -1
          go to 290
       end if

       go to 278
  279  continue
C      write(6,*)' Error in read: ia+1,natom,waven= ',ia+1,natom,waven
C      write(6,1060) waven,((opl(ip6,itemp),ip6=1,np6),
C    &                                      itemp=1,nftemp)
       nerr = nerr + 1
  278  continue
1065   continue
1069   continue

        do 451 it=1,17
        do 453 ip6=1,9
453     sumatdir(ip6,it)= max(1.d-99, sumatdir(ip6,it))
        write(6,452) t_at(it),((log10(sumatdir(ip6,it))-5.d0),ip6=1,9)
451     continue
452     format(f8.1,9f8.3)
       close(unit=11)
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
       write(6,*)' we have now identified opac(p6,temp) in ',nwatom_os
     &           ,' wavenumbers for metalopacities'
       write(6,*) ' between wmin and wmax:',waven_1,waven_2,' cm-1'
       write(6,*) ' there were ',nwtot,' OS frequencies in the interval'
       write(6,*)' the derivatives dadp_it = dln(opac)/dln(p6) are also'
     &           ,' determined in all 9 p6_points (for all wn and temp)'
       write(6,*) ' there were ',natom,' atomic OSlines in the interval'
C       open(unit=3,file='atoms.tmp',status='unknown')
C       do 1070 ia = natom,1,-1
C       write(3,1072) (opjva(ia,k),k=0,6)
C1070   continue
C1072   format(f10.3,1p6e11.4)


       return
       end
   

      subroutine atoms_head
* reads integrated metal opacity file in ascii format
      implicit real*16 (a-h,o-z)
C     implicit none
*
      integer np6,ntemp,maxpnt,nel
      parameter (np6 = 9)
      parameter (ntemp = 17)
      parameter (maxpnt = 153910)
      parameter (nel = 92)
*
      integer ispec(184),itemp,ip6,nwave,i,n,nwav
      real*16 waven
C      real abund(nel)
C      real opl(np6,ntemp)
C      real xite
      real*16 abund(nel)
      real*16 opl(np6,ntemp)
      real*16 xite
      character*60 filename
      character*80 string
*
c     data filename /'/ste1/uffegj/atoms/metals_sun_ascii.x03'/
      common /catoms_head/ p6(9),temp(17)
c     print *,'name of the input file ?'
c     read(*,'(a)') filename(1)
c     open(11,file=filename(1),form='unformatted',status='old')
c     print *,'name of the output file ?'
c     read(*,'(a)') filename(2)
c     open(41,file=filename,form='formatted',recl=1246,
c    &     status='unknown',readonly)
      open(13,file='metals.output',status='unknown')
*
      read(41,1080) string
      write(13,1080) string
* 1) info record
 1080 format(a80)
      write(6,1080)string
      read(41,1080) string
      write(13,1080) string
* 2) info record
      write(6,1080)string
      read(41,1000) ispec
      write(13,1000) ispec
* 3) ispec:  integers identifying all species included:
*            ex: 9201 = U II, 600 = C I    (max number=2*92)
 1000 format(184i5)
      write(6,1000)ispec
      read(41,1010) xite
      write(13,1010) xite
* 4) xite: [km/s]  microturbulence velocity (for line broadening)
 1010 format(f5.2)
      write(6,1010)xite
      read(41,1020) ip6     !1020 or *
      if(ip6.ne.np6) then
        write(6,*) ip6,np6,' Error: ip6 ne np6'
        stop 'Data file not as expected'
      endif
      write(13,1020) np6
* 5) np6: number of P6 values (must always=9)
 1020 format(i3)
      write(6,1020)np6
      read(41,1030) (p6(i),i=1,np6)
      write(13,1030) (p6(i),i=1,np6)
* 6) P6: [dyn] Damping pressures for "van der Waals" line broadening
*        to be computed as P(HI) + 0.42*P(HeI) + 0.85*P(H2)
 1030 format(1p,9e10.2,0p)
      write(6,1030) (p6(i),i=1,np6)
      read(41,1020) itemp
      if(itemp.ne.ntemp) then
        write(6,*) itemp,ntemp,' Error: itemp ne ntemp'
        stop 'Data file not as expected'
      endif
      write(13,1020) ntemp
* 7) ntemp: number of T values (must always=17)
      write(6,1020)ntemp
      read(41,1070) (temp(i),i=1,ntemp)
      write(13,1070) (temp(i),i=1,ntemp)
* 8) temp: [K] temperature values for excitation, ionization, chemical
*       equilibrium, line broadening (LTE)
 1070 format(17f8.0)
      write(6,1070) (temp(i),i=1,ntemp)
      read(41,1040) nwave
      write(13,1040) nwave
* 9) nwave will probably be 153910, the expected number of wavelength points
 1040 format(i10)
      write(6,1040)nwave
      read(41,1080)
      write(13,*) 'dummy record'
* 10) dummy record: unly used in single species files
      write(6,*)'dummy record'
      read(41,1080)
      write(13,*) 'dummy record'
* 11) dummy record: unly used in single species files
      write(6,*)'dummy record'
      read(41,1050,err=3879) abund
      write(13,1050) abund
* 12) abund: logarithmic number abundances used in the compilation of the file
*            on a scale where H = 12.00, 92 values, -99. means not used/known
 1050 format(92f7.2)
      write(6,1050)abund
      goto 3880
 3879 continue
      write(6,*)'abundances not found, probably single element file?'
 3880 continue
*
* 12 header lines ready, write 153910 opacity data recods:
*
      nwav=0
      write(6,*) ' now starting big read; maxpnt,np6,ntemp=',
     &     maxpnt,np6,ntemp
      write(13,*) ' now starting big read; maxpnt,np6,ntemp=',
     &     maxpnt,np6,ntemp
C     do n=1,maxpnt+1
C     read(41,1060,end=98,err=99) 
C    &   waven,((opl(ip6,itemp),ip6=1,np6),itemp=1,ntemp)
* waven: [cm-1] vacuum wave number
* opl: log(10) opacity, opacity in cm2/g stellar matter
*      opl = -30.0 or -40.0 means no significant metal line opacity found
* acurracy: Line data: b-b transition line data of neutral and singly-ionized
*           atoms were all adopted from VALD, January 1998. (cf A&AS 112, 525)
*           For 36 rare species no line data is there or the solar abundance=0
*           For the species where continuous opacities are considered,
*           the line opacity data [cm2/g] is complete to better than 1% of
*           the continuous opacity (at the relevant T and wavel). Based on the
*           experience from these species similar completness levels were
*           adopted for all other species.
C       if (n/10000*10000.eq.n .or.n.eq.1 .or.n.eq.maxpnt)
C    &      write(13,1060) waven,((opl(ip6,itemp),ip6=1,np6),
C    &                                      itemp=1,ntemp)
C1060 format(d22.16,153f8.3)
C       nwav=nwav+1
C     enddo
C     stop 'Error: more than maxpoint data records ???'
C  99 continue
C     write(6,*) ' Error in read at: nwav,waven = ',nwav,waven
C     write(13,*) ' Error in read at: nwav,waven = ',nwav,waven
C     write(13,1060) waven,((opl(ip6,itemp),ip6=1,np6),
C    &                                      itemp=1,ntemp)
C  98 continue
C     write(6,*) ' we read ',nwav,' wavenumber points'
C     write(13,*) ' we read ',nwav,' wavenumber points'
      return
      end


C
C------------------------------------------------------------
C                                                           I
      SUBROUTINE TPGREAD 
C                                                           I
C reads the input T-Pg structure if not a Marcs model       I
      implicit real*16 (a-h,o-z)
      include 'parameter.inc'
      CHARACTER ATMOS*60
      COMMON /CMODEL/TEFF,G,ABUND(16)  
     & ,TAUMOD(NDP),T(NDP),PE(NDP),PG(NDP),PRAD(NDP) 
     & ,PTURB(NDP),XKAPR(NDP),RO(NDP),EMU(NDP),PRESMP(NDP,99)
     & ,XL(20),ABSKA(20,NDP),SPRIDA(20,NDP),NTAU,NLP,ATMOS
      dimension phe(ndp),trpe(ndp),ptot(ndp),trphe(ndp)
      common /ctotabk/totabk(ndp,natms)
C
      open(unit=22,file=atmos,status='old')

       ntau = 0
       do 110 k=1,ndp   
       read(22,*,end=999) t(k), pg(k), pe(k)
       ntau = ntau+1
110    continue
999    continue

       do 120 i=1,ntau
       phe(i) = pg(i) * totabk(i,3)/(totabk(i,2)+totabk(i,3))
       trpe(i) = 0.1d0 * pe(i)
CV20   ptot(i) = pg(i) + pe(i)
       ptot(i) = pg(i)
       ptot(i) = 0.1d0 * ptot(i)     !1N/m^2 = 10 dynes/cm^2
       trphe(i) = (0.1d0 * phe(i)) / ptot(i)
C      if(i.eq.1.or.i.eq.27.or.i.eq.47)
C    &          write(6,782) i,t(i),pg(i),pe(i),trphe(i)
120   CONTINUE
782    format(i3,f8.1,1p4e12.3)

      call tstgem(t,ptot,ntau)
C
      return
      end

C                                                           I
C------------------------------------------------------------
C
C
      SUBROUTINE MODEL
C   
C---------------------------------------------------------- 
C                                                         I 
C  SUBROUTINE TO READ IN MODEL DATA FROM ARCHIV FILE.     I 
C                                                         I 
C  TEFF  : EFFECTIVE TEMPERATURE                          I 
C  G     : ACCELERATION OF GRAVITY                        I 
C  ABUND : ABUNDANCE OF H,HE,C,N,O,...... (NOT LOG() !)   I
C  NTAU  : NUMBER OF DEPTH POINTS                         I 
C  TAU   : STANDARD OPTICAL DEPTH                         I 
C  T       TEMPERATURES                                   I 
C  PE    : ELECTRON PRESSURES                             I 
C  PG    : GAS PRESSURES                                  I 
C  PREAD : RADIATION PRESSURES                            I 
C  PTURB : TURBULENT PRESSURES                            I 
C  XKPAPR: STANDARD ABSORPTION COEFFICIENTS               I 
C  RO    : DENSITIES                                      I 
C  EMU   : MEAN MOLECULAR WEIGHT (NUCLEONS PER PARTICLE)  I 
C  PRESMP: MOLECULAR PARTIAL PRESSURES (C2H2=15, HCN=16)  I
C        : (IN DYN/CM2) - NOT LOG() !                     I 
C  XL    : WAVELENGTHS FOR ABSKA, SPRIDA (ANGSTROEM)      I 
C  ABSKA : CONTINUOUS ABSORPTION COEFFICIENTS             I 
C  SPRIDA: CONTINUOUS SCATTERING COEFFICIENTS             I 
C                                                         I 
C---------------------------------------------------------- 
C
      implicit real*16 (a-h,o-z)
      include 'parameter.inc'
      dimension SRCH(2)
      COMMON /CMODEL/TEFF,G,ABUND(16)  
     & ,TAUMOD(NDP),T(NDP),PE(NDP),PG(NDP),PRAD(NDP) 
     & ,PTURB(NDP),XKAPR(NDP),RO(NDP),EMU(NDP),PRESMP(NDP,99)
     & ,XL(20),ABSKA(20,NDP),SPRIDA(20,NDP),NTAU,NLP,ATMOS
      COMMON /CMODEL2/P_MOL(NDP), P_NEU_HCNO(NDP), P_ION_HCNO(NDP)
     & ,P_NEU_HE(NDP), P_ION_HE(NDP), P_NON_HHECNO(NDP), P6_JON(NDP)
     & ,PG_JON(NDP), P6(NDP), PHE(NDP)
      COMMON /CTRAN/TAU(NDP),X(NDP),BPLAN(NDP),HSURF,Y1(6),JTAU
      CHARACTER SRCH*15,IDEN*15,ID2*115,ATMOS*60
      DATA SRCH /'O R R E C T I O','M O D E L  P A '/
C common for storage of input partial pressures
      character*8 namesave
      common/cppinp/ pressave(ndp,99),nsave,namesave(99)

C extra for inclusion of GEM routines:
      common/cmxread/pel_ion(ndp,16,4), elabund(ndp,16), 
     &                pprel(ndp,nspec)        !elabund=#atoms_i/atom_in_gas
      common /cgem/pres_gem(ndp,nspec)    !stored input model pressures in gem_indexes
      character name_gem*8
      common /cgemnames/natms_gem,nions_gem,nspec_gem,name_gem(nspec)
C atms,ions,spec ~ highest index of neutral atoms, ions, species total
      character*8 name_listmo(16),name_ex(250)
     &           ,name_comp1,name_comp2,name_comp3
      dimension ptot(ndp),trphe(ndp),trpe(ndp)
      dimension klistmo(99)
      common /cmarcsgem/ kmp(0:99),kidn(99),niden
      common /cabink/abink(ndp,nspec)
      real*16 junk1,junk2,junk3,junk4,junk5! Dummies

       DO 2430 I=1,NTAU
       DO 2430 J=1,NSPEC
2430     pprel(i,j) = 1.d-40
        DO 5113 I=1,NDP
        do 5131 j=1,nspec
        pres_gem(i,j) = -45.d0
5131    CONTINUE
        DO 5113 J=1,16
        DO 5114 JJ=1,4
        pel_ion(i,j,jj) = 1.0d-30
5114    CONTINUE
        DO 5113 K=1,99
        PRESMP(I,K) = -99.d0
5113    CONTINUE

       do 470 i=1,16
470    name_listmo(i) = '        '

C
      OPEN (UNIT=22,FILE=ATMOS,STATUS='OLD',readonly)
C
C
10      FORMAT(/' ATMOSPHERE = ',A60,/' (Teff,log(g),C/O) ='
     &         ,F6.0,F5.1,F8.2)
11      FORMAT(A15,A115)
12      FORMAT(29X,F7.0)
13      FORMAT(29X,E11.3)
14      FORMAT(24X,F14.5)
15      FORMAT(16F5.2)
811     FORMAT (A15)
C
C

C
C First count the number of depth points in this atmosphere:
        DO 830 LINE=1,10000
C
        READ(22,811,END=899) IDEN
        IPOS = 0
C        IPOS = INDEX(IDEN,SRCH(1))  !C O R R E C T I O N S  I N  L A S T
        IPOS = INDEX(IDEN,'R R E C')
        if (IDEN.NE.' C O R R E C T ' .and. 
     *     IDEN.NE.'C O R R E C T I') go to 830
C       write(6,811) iden
C       write(6,*) ' ipos = ',ipos
C        IF (IPOS.EQ.0) GO TO 830
        NLIN=0
831     continue
        NLIN = NLIN+1
        READ(22,811,END=899) IDEN
        IPOS = 0
        IPOS = INDEX(IDEN,'TAUROS')
        IF (IPOS.EQ.0 .and. NLIN.LE.99) GO TO 831
        DO 832 I=1,99                   ! 99=dimension of depth variables
        READ(22,811,END=899) IDEN
        IPOS = 0
        IPOS = INDEX(IDEN,'M O D E L  A T')
        IF (IPOS.NE.0) GO TO 833
        IF(IDEN(1:6).EQ.'      ') GO TO 833
        NTAU = I
832     CONTINUE
830     CONTINUE
833     CONTINUE
899     CONTINUE
C
        REWIND(22)
C
        IFLAG = 0
        LOGKILL = 0
        DO 4100 MT=1,10000
        READ(22,11,END=4199) IDEN,ID2
C
        IF (IDEN.EQ.'0M O D E L  P A' .OR. IDEN.EQ.' M O D E L  P A') 
     *    THEN
        IFLAG=IFLAG+1            !IFLAG=1
1005    CONTINUE
        READ(22,12) TEFF
        IF (TEFF.EQ.0.0) GO TO 1005
        READ(22,13) FLUX
        READ(22,14) G
        END IF
C
        IF(IDEN.EQ.'  LOG. ABUNDANC') THEN
        IFLAG=IFLAG+1            !IFLAG=2
        READ(22,11) IDEN,ID2
        READ(22,*) (ABUND(I),I=1,16)
        WRITE(6,10) ATMOS,TEFF,LOG10(G),10.**(ABUND(3)-ABUND(5))
        WRITE(6,41053) (ABUND(IE),IE=1,5),ABUND(15)
        AHA=ABUND(1)
        SUMABUND = 0.D0
        DO 4105 I=1,16
        ABUND(I)=10.**(ABUND(I)-AHA)
        SUMABUND = SUMABUND + ABUND(I)
4105    CONTINUE
        DO 4205 k=1,NDP
        DO 4205 I=1,16
        elabund(k,i) = ABUND(I)/SUMABUND
4205    CONTINUE
C       WRITE(6,41051) (ABUND(IE),IE=1,16)
        END IF
C41051   FORMAT(1P8E10.3)
41053   FORMAT('H,He,C,N,O,Fe:'8f8.3)
C
C
        IF ( IDEN.NE.'1M O D E L  A T' .AND. IDEN.NE.' M O D E L  A T') 
     *    GO TO 4101
        IFLAG=IFLAG+1            !IFLAG=3
1006    READ(22,11,END=4199) IDEN,ID2
        LOOP=LOOP+1
        IF (LOOP.LE.99.AND.IDEN(1:5).NE.'    K'.AND.IDEN(1:5).NE.'   K '
     &       .AND.IDEN(1:5).NE.'  K  '.AND.IDEN(1:5).NE.' K   '
     &       .AND.IDEN(1:5).NE.'0   K'.AND.IDEN(1:5).NE.'0  K '
     &   .AND.IDEN(1:5).NE.'0 K  '.AND.IDEN(1:5).NE.'0K   ') GO TO 1006
        LOOP=0
C        READ(22,11,END=4199) IDEN,ID2
C        IF (IDEN(1:1).EQ.' ') READ(22,11,END=4199) IDEN,ID2
C       WRITE(6,*) ' IFLAG3 = ',IFLAG
        DO 1410 I=1,NTAU
        READ(22,*) K1,TAUMOD(I),TAUS,Z,T(I),PE(I)
     &  ,PG(I),PRAD(I),PTURB(I),XKAPR(I)
        TAU(I)=TAUMOD(I)
C        IF(I.EQ.1.OR.I.EQ.NTAU) WRITE(6,*)
C     &  K1,TAU(I),TAUS,Z,T(I),PE(I),PG(I),XKAPR(I)
1410    CONTINUE
4101    CONTINUE
        JTAU=NTAU
        IF ( IDEN.NE.'1T H E R M O D ' .AND. IDEN.NE.' T H E R M O D ') 
     *    GO TO 4102
1007    READ(22,11,END=4199) IDEN,ID2
        LOOP=LOOP+1
        IF (LOOP.LE.99.AND.IDEN(1:5).NE.'    K'.AND.IDEN(1:5).NE.'   K '
     &       .AND.IDEN(1:5).NE.'  K  '.AND.IDEN(1:5).NE.' K   '
     &       .AND.IDEN(1:5).NE.'0   K'.AND.IDEN(1:5).NE.'0  K '
     &   .AND.IDEN(1:5).NE.'0 K  '.AND.IDEN(1:5).NE.'0K   ') GO TO 1007
        LOOP=0
        IFLAG=IFLAG+1            !IFLAG=4
C       WRITE(6,*) ' IFLAG4 = ',IFLAG
C        READ(22,11,END=4199) IDEN,ID2
C        READ(22,11,END=4199) IDEN,ID2
C        IF (IDEN(1:1).EQ.' ') READ(22,11,END=4199) IDEN,ID2
        DO 4112 I=1,NTAU
        READ(22,*) IK1,TR,RO(I),EMU(I),CP,CV,AGRAD,Q,U,V,FCF,IK2
C        IF(I.EQ.1.OR.I.EQ.NTAU) WRITE(6,*)
C     &  IK1,TR,RO(I),EMU(I),CP,CV,AGRAD
4112    CONTINUE
4102    CONTINUE

        IPOS = INDEX(IDEN,'L O G A R I T')
        IF ((IPOS.EQ.0).OR.(LOGKILL.EQ.1)) GO TO 4107
        IFLAG=IFLAG+1  !IFLAG=5
C       WRITE(6,*) ' IFLAG5 = ',IFLAG
        LOGKILL = 1    !To prevent another search in next loop (AB 1995-05)

        nmol = 0
        nfail = 0
        nsave = 0
2135    CONTINUE        !read (next) block of partial pressures
        ILOOP = 0
5831    READ(22,811,END=999) IDEN
        IF(IDEN.EQ.' A B S O R P T ' .or.
     &               IDEN.EQ.' P A R T I A L') GO TO 4106  !assumed end of PP blocks
        ILOOP = ILOOP + 1
        IPOS = INDEX(IDEN,' K ')
        IF (ILOOP.GE.100) STOP ' ***Error: I found no K in PP'
        IF (IPOS.EQ.0) GO TO 5831
C when here, I identified the line with names of molecules for this block
           backspace(22)
           read(22,2143) (name_listmo(i),i=1,16)
2143  FORMAT(5x,16a8)
           innam = 0
           do 2161 i=1,16
2161       if (name_listmo(i).ne.'        ') innam = innam+1
C innam is number of molecule names in this read block
           niden = nmol

           do 2121 j=1,innam
            kn = 0
            name_comp2 = '        '
            name_comp3 = name_listmo(j)
           do 2127 kc = 1,8
           if(name_comp3(kc:kc).ne.' ') then
            kn = kn+1
            name_comp2(kn:kn) = name_comp3(kc:kc)
           end if
2127       continue
           if(name_comp2.eq.'P(H)    ') name_comp2='H       '
           if(name_comp2.eq.'K       ') name_comp2='XXXXXXXX'
C name_comp2 is now the read name_listmo(j), but starting in position 1
           nsave = nsave+1
           namesave(nsave) = name_comp2

           idsucces = 0
           ipos = 0
           do 2122 i=1,nspec
           name_comp1 = name_gem(i)
C just for the comparison with GEM-names:
C (nothing is changed in any name arrays here)
             if(name_comp1.eq.'HN      ') name_comp1 = 'NH      '
             if(name_comp1.eq.'HNa     ') name_comp1 = 'NaH     '
             if(name_comp1.eq.'HMg     ') name_comp1 = 'MgH     '
             if(name_comp1.eq.'HSi     ') name_comp1 = 'SiH     '
             if(name_comp1.eq.'CSi     ') name_comp1 = 'SiC     '
             if(name_comp1.eq.'ClH     ') name_comp1 = 'HCl     '
             if(name_comp1.eq.'CHN     ') name_comp1 = 'HCN     '
             if(name_comp1.eq.'CSi2    ') name_comp1 = 'Si2C    '
           ipos = index( name_comp2,name_comp1 )
           if(ipos.ne.0) then        !we identified the listmo molecule as a GEM name
              niden = niden + 1
              if(niden.gt.99) stop ' increase dimension invl niden !'
              klistmo(niden) = j
              kidn(niden) = i
              idsucces = 1
C             write(6,2128) name_comp1,name_comp2
C             write(6,2129) niden,i,j,kidn(niden),name_gem(i)
2128          format(' name{gem,listmo}=',2(2x,a8))
2129          format(' niden,i,j,kidn(niden), name_gem(i):',4i4,2x,a8)
           end if
2122    CONTINUE
           if(idsucces.eq.0) then 
                 nfail = nfail + 1
                 name_ex(nfail) = name_comp2
           end if
2121    CONTINUE

C At this point namesave(1-nsave) contain the nsave read pp names from input model.
C The corresponding partial pressures (pp) are being saved in pressave(1-nsave).
C name_gem(kidn{1-niden}) contain those of the input model names which has a gem-name,
C while name_ex(1-nfail) contain the input model names not identified in GEM.
C In contrast to pressave, pres_gem(i,kidn(j),j=1,niden) contain only the subset of
C input model pressures of molecules which are in common between GEM and the input model.
C Note that pres_gem is not the GEM-computed pressures, but input pressures.

        DO 2130 I=1,NTAU
        READ(22,*) K1,(PRESMP(I,K),K=1,innam)
        do 2140 k=1,innam
        js=nsave-innam+k
2140    PRESSAVE(I,js) = 10.d0**PRESMP(I,K)
C       if(i.eq.1) then
C         write(6,2142) (namesave(nsave-innam+k),k=1,innam)
C         write(6,2146) 
C    &  (dlog10( max(1.d-40,pressave(i,nsave-innam+k))),k=1,innam)
CC    &  (pressave(i,nsave-innam+k),k=1,innam)
C       end if
2130    CONTINUE
2142    format(2x,8a8/,8a8)
2146    format(8f8.2/,8f8.2)

         if (niden-nmol.eq.0) then
C         write(6,*) 
C    &    ' we found no gem_molecules in this read block of PP'
          go to 2135
         end if
        DO 2131 I=1,NTAU
        do 2131 j=nmol+1,niden
        pres_gem(i,kidn(j)) = presmp(i,klistmo(j))
2131    CONTINUE

2134    FORMAT(i3,16f8.3)
C         write(6,2143)(name_gem(kidn(j)),j=nmol+1,niden)
          DO 2133 k=1,NTAU,20
C         write(6,2134)k,(pres_gem(k,kidn(j)),j=nmol+1,niden)
2133    CONTINUE

        NMOL = niden

        go to 2135

4106    CONTINUE

        WRITE(6,5110) NMOL,NTAU
5110    FORMAT ( ' We identified the following',i3,
     &     ' molecular part.pres.in',i3,' depth layers:')
        write(6,2144)(name_gem(kidn(j)),j=1,niden)
2144    format(15a5,/15a5,/15a5,/15a5)
        write(6,2141) nfail,(name_ex(nf),nf=1,nfail)
2141    format(' We didnt identify in GEM the following',i3,
     &    ' molecules:'/,10a8/,10a8/,10a8/,10a8/,10a8/,10a8)
        DO 2139 I=1,NTAU
        do 2139 j=1,niden
        presmp(i,j) = pres_gem(i,kidn(j))
2139    CONTINUE

C Now, put into the indexes kmp({1-niden}) the GEM index of the molecules which
C we have absorption coefficients for:
        do 2138 j=1,niden
        if(name_gem(kidn(j)).eq.'CO      ') kmp(0) = kidn(j)
        if(name_gem(kidn(j)).eq.'CH      ') kmp(1) = kidn(j)
        if(name_gem(kidn(j)).eq.'C2      ') kmp(2) = kidn(j)
        if(name_gem(kidn(j)).eq.'SiO     ') kmp(3) = kidn(j)
        if(name_gem(kidn(j)).eq.'CN      ') kmp(4) = kidn(j)
        if(name_gem(kidn(j)).eq.'TiO     ') kmp(5) = kidn(j)
        if(name_gem(kidn(j)).eq.'H2O     ') kmp(6) = kidn(j)
        if(name_gem(kidn(j)).eq.'C2H2    ') kmp(7) = kidn(j)
        if(name_gem(kidn(j)).eq.'HCN     ') kmp(8) = kidn(j)
        if(name_gem(kidn(j)).eq.'C3      ') kmp(9) = kidn(j)
        if(name_gem(kidn(j)).eq.'H2      ') kmp(10) = kidn(j)
        if(name_gem(kidn(j)).eq.'He      ') kmp(11) = kidn(j)
        if(name_gem(kidn(j)).eq.'C2H2    ') kmp(12) = kidn(j)
        if(name_gem(kidn(j)).eq.'CH4     ') kmp(13) = kidn(j)
        if(name_gem(kidn(j)).eq.'CS      ') kmp(14) = kidn(j)
        if(name_gem(kidn(j)).eq.'C2H     ') kmp(15) = kidn(j)
        if(name_gem(kidn(j)).eq.'OH      ') kmp(16) = kidn(j)
2138    continue
        kmp(11) = 3     !He is not in the listmo
        kmp(8) = 372   !HCN as CHN is not in the listmo

C        write(6,*) 
C     &  ' The opacity bearing species were identified as:'
C        do 2136 j=0,16
C        write(6,2137) j,kmp(j),name_gem(kmp(j))
C2137    format(i3,i4,2x,a4)
2136    continue


4107    CONTINUE


        IF(IDEN.NE.' LAMBDA FOR ABS') GO TO 4103
        IFLAG=IFLAG+1            !IFLAG=6
C       WRITE(6,*) ' IFLAG6 = ',IFLAG
        READ(22,*) NLP                      !usually NLP=20 in input, #20 is 8000.AA
C       WRITE(6,811) IDEN
C       WRITE(6,*) ' NLP = ',NLP
        NLP=NLP-1
C       WRITE(6,*) ' NLP = ',NLP
        READ(22,*) (XL(J),J=1,NLP),XLEXTRA   !cont.wavelength scale is 19 freq.
C        WRITE(6,*) (XL(J),J=1,NLP),XLEXTRA
        READ(22,11) IDEN,ID2
        READ(22,11) IDEN,ID2
        DO 4202 I=1,NTAU
        READ(22,*) ITAU,TAUX
C        WRITE(6,*) ITAU,TAUX
        READ (22,*) (ABSKA(K,I),K=1,NLP),ABSEXTRA
        READ (22,*) (SPRIDA(K,I),K=1,NLP),SPREXTRA
4202    CONTINUE
4103    CONTINUE
C
4100     CONTINUE
4199     CONTINUE
C
      CLOSE (22)
      CLOSE (2)


C at this point pres_gem(i,j) are 10^-49.0 if not identified in the input model,
C or the log10(input partial pressure / (pg+pe)) if identified in the input model.
C  pprel(i,j) = 10.**presrel is then the input relative partial pressure in 
C model layer i for all nspec species, j=1,nspec.

      DO 2410 I=1,NTAU
        do 166 j=1,nspec
CV20     presrel= pres_gem(i,j) - dlog10( max(1.d-40,(pg(i)+pe(i))) )
         presrel= pres_gem(i,j) - log10( max(1.q-40,pg(i)) )
C        if(i.eq.1 .or. i.eq.27 .or. i.eq.47) then
C        if (presrel.ge.-9.0d0)
C    & write(6,266) name_gem(j),j,presrel,10.**presrel,pres_gem(i,j)
C        end if
         pprel(i,j) = 10.**presrel
166     continue
266      format(a8,i4,f7.2,f8.5,f8.3)

       phe(i) = log10
     &          (pg(i) * elabund(1,2)/(elabund(1,1)+elabund(1,2)))
       trpe(i) = 0.1d0 * pe(i)
CV20   ptot(i) = pg(i) + pe(i)
       ptot(i) = pg(i)
       ptot(i) = 0.1d0 * ptot(i)     !1N/m^2 = 10 dynes/cm^2
       trphe(i) = (0.1d0 * 10.**phe(i)) / ptot(i)
C      if(i.eq.1.or.i.eq.27.or.i.eq.47)
C    &          write(6,782) i,t(i),pg(i),pe(i),trphe(i)
2410  CONTINUE
782    format(i3,f8.1,1p4e12.3)


       call tstgem(t,ptot,ntau)

C
      IF (IFLAG.LT.5) THEN
        write(6,*)
     &  ' ERROR: problems with reading of model-input. IFLAG=',iflag
        STOP
      ELSE
        write(6,*) 
     &  ' successful reading of input model atmosphere; iflag=',iflag
      END IF
      GO TO 998
999     CONTINUE
      WRITE(6,*) ' MODEL STOP: HOW COULD THIS BE ???????'
998     CONTINUE

C
      RETURN
C  
 
      END
C
C --------------------------------------------------------------
       SUBROUTINE OSOPACITY                                    !
C                                                              |
C --------------------------------------------------------------
C                                                              |
C  'OSOPACITY' reads the OS-files (absorption coefficient in   |
C  a number of wavelength points) for all molecules considered |
C  in the spectrum calculation. 'OSOPACITY' is called from     |
C  MAIN only once; it reads OPAC for all molecules and for all |
C  OS wavenumber points and store them in an array.            |
C  It also establish a wavenumber array, WNOS(1-NWNOS), which  |
C  is stored in common-block CWNOS, and which contains the only|
C  wavenumbers in which a spectrum can be computed from now on.|
C -------------------------------------------------------------|
C
C
      implicit real*16 (a-h,o-z)
      include 'parameter.inc'
C
      CHARACTER ATMOS*60
C
      COMMON /CMODEL/TEFF,G,ABUND(16)
     & ,TAUMOD(NDP),TATMOS(NDP),PE(NDP),PG(NDP),PRAD(NDP)
     & ,PTURB(NDP),XKAPR(NDP),RO(NDP),EMU(NDP),PRESMP(NDP,99)
     & ,XL(20),ABSKA(20,NDP),SPRIDA(20,NDP),NTAU,NLP,ATMOS
 
      common /cabink/abink(ndp,nspec)
      common /cgeminp/tcall(ndp),pcall(ndp),phecall(ndp),pecall(ndp) 
      character name_gem*8
      common /cgemnames/natms_gem,nions_gem,nspec_gem,name_gem(nspec)

      common /cprespp/prespp(ndp,nspec) 


       DO 3151 kd=1,NTAU
       ptrans = 10.d0*pcall(kd)                             !pg+pe in dyn/cm^2
       DO 3153 km=1,nspec
       PRESPP(kd,km) = ptrans * abink(kd,km)      !*relative part.pres. = pp in dyn/cm^2
3153   continue
3151   continue

C
        RETURN
        END
C
C
Ckemi      PROGRAM KEMI
Ckemi      implicit real*16 (a-h,o-z)
CkemiC
Ckemi       include 'parameter.inc'
CkemiC
CkemiC     parameter(nspec=892)
Ckemi      common /cgem/pres_gem(ndp,nspec)
Ckemi      common /cgemnames/natms_gem,nions_gem,nspec_gem,name_gem(nspec)
CkemiC atms,ions,spec ~ highest index of neutral atoms, ions, species total
Ckemi      character name_gem*8,adum*45
Ckemi      common /cabink/abink(ndp,nspec)
Ckemi      dimension ptot(ndp),pp_sum(ndp)
CkemiC     common /cpartryck/partryck(ndp,nspec)           !gem computed pp
Ckemi      common /cpartp/partp(ndp,nspec)           !gem computed pp
Ckemi      common /cprespp/prespp(ndp,nspec)           !gem computed pp
Ckemi      common /ctotabk/totabk(ndp,natms)
Ckemi      dimension presmo(99)
Ckemi      COMMON /CMODEL/TEFF,G,ABUND(17)  
Ckemi     & ,TAUMOD(NDP),T(NDP),PE(NDP),PG(NDP),PRAD(NDP) 
Ckemi     & ,PTURB(NDP),XKAPR(NDP),RO(NDP),EMU(NDP),PRESMP(NDP,99)
Ckemi     & ,XL(20),ABSKA(20,NDP),SPRIDA(20,NDP),NTAU,NLP,ATMOS
Ckemi      common /ckdtpe/dpex,kdtpe
Ckemi
Ckemi      common /fullequilibrium/ partryck(ndp,maxmol),
Ckemi     &  xmettryck(ndp,maxmet),xiontryck(ndp,maxmet)
Ckemi
Ckemi      character namemol*4
Ckemi      dimension namemol(18),pplist(3,ndp,18)
Ckemi      dimension pe_gem(ndp),ptot1(ndp),dptot(ndp),pe1(ndp)
Ckemi     &   ,dptot2(ndp),dpe2(ndp),dpe(ndp)

Ckemi       open(unit=5,file='kemi.input',status='old')
C

      subroutine test_Tsuji(tt,pg)
      implicit real*16 (a-h,o-z)
* NOT THE ORIGINAL VERSION> CHANGED TO STUDY SOME SPECIES
* test of tsuji partial pressure

* version     6.9.94   Ch.Helling

      include 'parameter.inc'
      character*128 molec,met
      common/files/molec,met
      common /fullequilibrium/ partryck(ndp,maxmol),
     &  xmettryck(ndp,maxmet),xiontryck(ndp,maxmet),
     *  partp(ndp,0:maxmol)
      common /tsuji/ nattsuji,nmotsuji,parptsuji(500),abtsuji(17,ndp)
      common/cmolrat/fold(ndp,8),molold,kl

* source file tsuji_big.data

*       met = 'tsuji01.dat'
*       molec = 'tsuji02.dat'
      met ='data/tsuji.atoms'
      molec ='data/tsuji_big.data'

        pein=1
        temp=tt
        ttemp=5040/temp
        pgas=pg

        call eqmol(temp,pein,pgas,pe)

        call takemolec(kl)

        pein=pe
        
C
      return
      end

      subroutine partf(jatom,ion,temp,ntemp,u,ndim)
      implicit real*16 (a-h,o-z)

*
* yields partition functions with polynomial data from
* ref. Irwin, A.W., 1981, ApJ Suppl. 45, 621.
* ln u(temp)=sum(a(i)*(ln(temp))**i) 0<=a<=5
* Input:
* jatom = element number in periodic table
*       NOT implemented for molecules; see ref. in ref.
* ion   = 1 for neutral, 2 for once ionized and 3 for twice ionized
* temp  = array with ntemp values of the temperature
* ntemp = number of temperatures for which partf. is calculated
* Output:
* u     = partf. (linear scale) for iat,ion at the ntemp temperatures
*
C      implicit none
*
      integer ntemp,i,j,k,jatom,ion,iread
      real*16 temp(ndim),u(ndim),sp,spec
      real*16 a(0:5,1:3,1:92),ulog,t,aa(0:5)
*
      save iread,a
*
      data iread /0/
*
      if(ntemp.gt.ndim) then 
          print*,'ntemp,ndim = ',ntemp,ndim
          stop 'ntemp > ndim in partf()'
      end if

      if(iread.ne.1) then
* read data if first call:

        open(67,file= 'data/irwin.dat'
     *    ,status='old',readonly)

        read(67,*)
        read(67,*)
        do 20 j=1,92
          do 10 i=1,3
            if(j.eq.1.and.i.eq.3) goto 10
            sp=qfloat(j)+qfloat(i-1)/100.
            read(67,*) spec,aa
            if(sp.ne.spec) then
              print *,'sp.ne.spec:',sp,spec
              stop 'Unexpected data read'
            else
              do 1 k=0,5
                a(k,i,j)=aa(k)
    1         continue
            endif
   10     continue
   20   continue
        close(67)
        iread=1
      endif


      do 30 i=1,ntemp
        if(temp(i).lt.600.) then
          stop 'partf; temp<600 K'
        else if(temp(i).gt.16000.) then
          print*,i,temp(i)
          stop 'partf; temp>16000 K'
        endif
C        t=log(dble(temp(i)))
        t=log(temp(i))
        ulog=   a(0,ion,jatom)+
     &       t*(a(1,ion,jatom)+
     &       t*(a(2,ion,jatom)+
     &       t*(a(3,ion,jatom)+
     &       t*(a(4,ion,jatom)+
     &       t*(a(5,ion,jatom))))))
        u(i)=exp(ulog)
   30 continue

      return
      end


      subroutine takemolec(kk)
      implicit real*16 (a-h,o-z)
*
* this routine is to be used after a call to jon,
*  -if eqmol has been called also-, in order to get the pressures
* one needs placed in the common 'fullequilibrium'.
* This is the routine to change if one wants more or/and other
* molecular pressure to be kept.
* The values in the common fullequilibrium can then be used for
* computation of opacities, printing, etc.
* kk is the depth point to be adressed.
* 020890 BPlez
*
* 13.09.94 Ch.Helling : change to be able to calculate all possible
*                       molecules in tsuji_big.data
*


      include 'parameter.inc'
      common /tsuji/ nattsuji,nmotsuji,parptsuji(500),abtsuji(17,ndp)
      common /fullequilibrium/ partryck(ndp,maxmol),
     &  xmettryck(ndp,maxmet),xiontryck(ndp,maxmet),
     *  partp(ndp,0:maxmol)
      COMMON/CARC3/F1P,F3P,F4P,F5P,HNIC,PRESMO(33)
* store the metals and ions
* they are not yet indexed like in jon.

      if(nattsuji.gt.maxmet) stop 'takemolec: maxmet too small'
      do 30 j=1,nattsuji
        xmettryck(kk,j)=parptsuji(j)
        xiontryck(kk,j)=parptsuji(j+nattsuji)
30    continue
*
* the indexes from 1 to 33 for partryck correspond to
* the ones in presmo (common carc3).

      partryck(kk,1)=parptsuji(2*nattsuji+1)
      partryck(kk,2)=parptsuji(2*nattsuji+2)
* H2+ (not self-consistent)
      partryck(kk,3)=presmo(3)
      do 10 j=4,16
        partryck(kk,j)=parptsuji(2*nattsuji+j-1)
10    continue
      partryck(kk,17)=1.e-20
      if (maxmol.lt.342) stop 'takemolec: maxmol too small !'
  
* number of molecules in molecBP+CB: 340 / 2.11.94

      do 20 j=18,342
       partryck(kk,j)=parptsuji(2*nattsuji+j-2)
       partryck(kk,j)=max(partryck(kk,j),1.0d-99)
20    continue

*****
      return
      end


*
      subroutine eqmol(tt,pein,pgas,pe)
      implicit real*16 (a-h,o-z)
      include 'parameter.inc'
*     
c  ce programme calcule l'equilibre d'atomes et molecules dont la lise
c   est contenue dans des fichiers METAL.DAT et MOLE.DAT .
c                                                       b. plez   7/7/87
C                                                       modif Pe- 6/11/87
C                                                       modif    22/11/88
* for a change of delta(D0) of the dissociation energy (in eV) of a molecule,
* change the c(2) coefficient in log(Kp) expansion by -delta(D0)
*
* VERSION FOR INCLUSION IN THE MARCS CODE. Uses only one depth point at a time
*                    BPz 010890
* test version for determining minimum input file for molecules
* BPz 161190
* 
* Ch.Helling 070994

        character*128 MET,MOLEC
        common/files/molec,met
        real*16       Z(51)  
ccc     real*4       EPRESLOG(50),PRESLOG(50,350)
        REAL*16 IP,KP,KPLOG,IPI,NNH,NHE
        doubleprecision   MOL
        COMMON/COMFH1/C(600,5),NELEM(5,600),NATOM(5,600),
     2   MMAX(600),PPMOL(600),APMLOG(600),MOL(600),
     3   IP(100),CCOMP(100),UIIDUI(100),P(100),FP(100),KP(100),NELEMX
     4   (50),EPS,SWITER,NMETAL,NMOL,NIMAX
        COMMON/COM6/TETAEF,GLOG,ASASOL,NHE,TIT(9)
        COMMON/VAL/PPG(600,50)
        DIMENSION YA(525),YB(525),YC(525),YD(525),ELEMNT(100),
     2 CCLOG(100),G0(100),G1(100),NATOMM(5),NELEMM(5),XP(50,100),
     3 PPA(51),PC13(51),PPH(51),PB(51),PO(51),PTI(51),ELEM(99)

* commons from injon
      COMMON/CI5/abmarcs(17,ndp),ANJON(17,5),Hz(5),PARTz(17,5),
     & DXIz,F1z,F2z,F3z,F4z,F5z,XKHMz,XMHz,XMYz(ndp)
C      real absc,abti,abv,abmn,abco
C     common /auxabund/ absc,abti,abv,abmn,abco
* common for partial pressures
      common /tsuji/ nattsuji,nmotsuji,parptsuji(500),abtsuji(17,ndp)
C nmotsuji is the number of molecules in Tsujis equilibrium routine (==nmol here)
      logical first,test
      DATA ELEMNT(99),ELEMNT(100)/'E-','H*'/
      data first/.true./

C      print 129, first
C129   format(' in eqmol, first=',l10)

        IND=1
        ECONST=4.342945E-1
        AVO=0.602217E+24
        SPA=0.196E-01
        GRA=0.275423E+05
        AHE=0.100E+00
C
c----lecture de quelques renseignements dans FOR009.DAT------------------------
C 
ccc     READ (9,5000) NMETAL,NIMAX,EPS,SWITER,IT
ccc        nmetal=38
        nimax=3000
        eps=0.001
        switer=1
        it=1
***
*  ! Nbre d'atomes,d'iterations,rcitere de convergence,?,Nbre couches atmospheriques
        IF(NIMAX.GT.0) GOTO 999
        NIMAX=50
        EPS=0.001
c
999     continue
cccc    WRITE(6,6102)
C
C----lecture du fichier contenant les atomes----------------------------------- 
C
        
      if (first) then
*
C     met ='/home/ast/uffegj/tsuji.atoms_AndersGrev'
C     molec ='/home/ast/uffegj/jens/comfits.data'

        OPEN(UNIT=26,FILE=MET,STATUS='OLD',readonly)
        rewind(26)
        read (26,*) nmetal
        write(6,*)
     &'Number of elements read into eqmol from tsuji.atoms: '
     &,nmetal
        nattsuji=nmetal
        DO 1002 I=1,NMETAL
        READ (26,5001) ELEMXI,NELEMI,IPI,IG0,IG1,CCLOGI 
*! nom(code),Nbre d'e-,pot d'ionisation,?,?,log de l'abondancce

        NELEMX(I)=NELEMI
        ELEM(I)=ELEMXI
        ELEMNT(NELEMI)=ELEMXI
        IP(NELEMI)=IPI
        UIIDUI(NELEMI)=DFLOAT(IG1)*0.661/DFLOAT(IG0)
        G0(NELEMI)=IG0
        G1(NELEMI)=IG1
        CCLOG(NELEMI)=CCLOGI
        CCOMP(NELEMI)=EXP(CCLOGI/ECONST)
ckeep for debug
ccc     WRITE(6,6103) ELEMXI,NELEMI,IPI,IG0,IG1,CCLOGI
c 
1002   CONTINUE

      close(26)


C
C Compare the elemental abundances read in from the MARCS input
C file to those read in from tsuji01.dat
C
  
* elements are in Tsuji called with their electron nummber

C     write(6,*)'abtsuji(i),cclog(i)+12:'
C     write(6,*) 'H  : abinp(1)= ',abtsuji(1),
C    & ' AndGrv(1)= ',(cclog(1)+12)
C     write(6,*) 'He : abinp(2)= ',abtsuji(2),
C    & ' AndGrv(2)= ',(cclog(2)+12)
C     write(6,*) 'C  : abinp(3)= ',abtsuji(3),
C    & ' AndGrv(6)= ',(cclog(6)+12)
C     write(6,*) 'N  : abinp(4)= ',abtsuji(4),
C    & ' AndGrv(7)= ',(cclog(7)+12)
C     write(6,*) 'O  : abinp(5)= ',abtsuji(5),
C    & ' AndGrv(8)= ',(cclog(8)+12)
C     write(6,*) 'Na : abinp(7)= ',abtsuji(7),
C    & ' AndGrv(11)= ',(cclog(11)+12)
C     write(6,*) 'Mg : abinp(8)= ',abtsuji(8),
C    & ' AndGrv(12)= ',(cclog(12)+12)
C     write(6,*) 'Al : abinp(9)= ',abtsuji(9),
C    & ' AndGrv(13)= ',(cclog(13)+12)
C     write(6,*) 'Si : abinp(10)= ',abtsuji(10),
C    & ' AndGrv(14)= ',(cclog(14)+12)
C     write(6,*) 'S  : abinp(11)= ',abtsuji(11),
C    & ' AndGrv(16)= ',(cclog(16)+12)
C     write(6,*) 'K  : abinp(12)= ',abtsuji(12),
C    & ' AndGrv(19)= ',(cclog(19)+12)
C     write(6,*) 'Ca : abinp(13)= ',abtsuji(13),
C    & ' AndGrv(20)= ',(cclog(20)+12)
C     write(6,*) 'Cr : abinp(14)= ',abtsuji(14),
C    & ' AndGrv(24)= ',(cclog(24)+12)
C     write(6,*) 'Fe : abinp(15)= ',abtsuji(15),
C    & ' AndGrv(26)= ',(cclog(26)+12)
C     write(6,*) 'Ni : abinp(16)= ',abtsuji(16),
C    & ' AndGrv(28)= ',(cclog(28)+12)
C     write(6,*) 'Ti : abinp(17)= ',abtsuji(17),
C    & ' AndGrv(22)= ',(cclog(22)+12)
    
C     rmag = 10.**(abtsuji(4)-cclog(7)-12.)
C    &  + 10.**(abtsuji(5)-cclog(8)-12.)
C    &  + 10.**(abtsuji(7)-cclog(11)-12.)
      rmag = 10.**(abtsuji(7,1)-cclog(11)-12.)
     &  + 10.**(abtsuji(8,1)-cclog(12)-12.)
     &  + 10.**(abtsuji(9,1)-cclog(13)-12.)
     &  + 10.**(abtsuji(10,1)-cclog(14)-12.)
     &  + 10.**(abtsuji(11,1)-cclog(16)-12.)
     &  + 10.**(abtsuji(12,1)-cclog(19)-12.)
     &  + 10.**(abtsuji(13,1)-cclog(20)-12.)
     &  + 10.**(abtsuji(14,1)-cclog(24)-12.)
     &  + 10.**(abtsuji(15,1)-cclog(26)-12.)
     &  + 10.**(abtsuji(16,1)-cclog(28)-12.)
     &  + 10.**(abtsuji(17,1)-cclog(22)-12.)
      rmaglog = log10(rmag/11.)
      write(6,*)' log average ratio of Marcs-input div. tsuji(~solar) ='
     &    ,rmaglog

C The elemental abundances read in from tsuji01.dat
C are set to those read in from the input file (AB 1995-05)

       do 1102 jab = 3,98
       cclog(jab) = cclog(jab) + rmaglog
1102   continue
C       
       cclog(1) = abtsuji(1,1) - 12.
       cclog(2) = abtsuji(2,1) - 12.
       cclog(6) = abtsuji(3,1) - 12.
       cclog(7) = abtsuji(4,1) - 12.
       cclog(8) = abtsuji(5,1) - 12.
       cclog(11)= abtsuji(7,1) - 12.
       cclog(12)= abtsuji(8,1) - 12.
       cclog(13)= abtsuji(9,1) - 12.
       cclog(14)= abtsuji(10,1)- 12.
       cclog(16)= abtsuji(11,1)- 12.
       cclog(19)= abtsuji(12,1)- 12.
       cclog(20)= abtsuji(13,1)- 12.
       cclog(24)= abtsuji(14,1)- 12.
       cclog(26)= abtsuji(15,1)- 12.
       cclog(28)= abtsuji(16,1)- 12.
       cclog(22)= abtsuji(17,1)- 12.
*
C      write(6,*) ' control for right abundances:'
C      write(6,*) 'H  : abinp(1)= ',abtsuji(1),
C     & ' AndGrv(1)= ',(cclog(1)+12)
C      write(6,*) 'He : abinp(2)= ',abtsuji(2),
C     & ' AndGrv(2)= ',(cclog(2)+12)
C      write(6,*) 'C  : abinp(3)= ',abtsuji(3),
C     & ' AndGrv(6)= ',(cclog(6)+12)
C      write(6,*) 'N  : abinp(4)= ',abtsuji(4),
C     & ' AndGrv(7)= ',(cclog(7)+12)
C      write(6,*) 'O  : abinp(5)= ',abtsuji(5),
C     & ' AndGrv(8)= ',(cclog(8)+12)
C      write(6,*) 'Na : abinp(7)= ',abtsuji(7),
C     & ' AndGrv(11)= ',(cclog(11)+12)
C      write(6,*) 'Mg : abinp(8)= ',abtsuji(8),
C     & ' AndGrv(12)= ',(cclog(12)+12)
C      write(6,*) 'Al : abinp(9)= ',abtsuji(9),
C     & ' AndGrv(13)= ',(cclog(13)+12)
C      write(6,*) 'Si : abinp(10)= ',abtsuji(10),
C     & ' AndGrv(14)= ',(cclog(14)+12)
C      write(6,*) 'S  : abinp(11)= ',abtsuji(11),
C     & ' AndGrv(16)= ',(cclog(16)+12)
C      write(6,*) 'K  : abinp(12)= ',abtsuji(12),
C     & ' AndGrv(19)= ',(cclog(19)+12)
C      write(6,*) 'Ca : abinp(13)= ',abtsuji(13),
C     & ' AndGrv(20)= ',(cclog(20)+12)
C      write(6,*) 'Cr : abinp(14)= ',abtsuji(14),
C     & ' AndGrv(24)= ',(cclog(24)+12)
C      write(6,*) 'Fe : abinp(15)= ',abtsuji(15),
C     & ' AndGrv(26)= ',(cclog(26)+12)
C      write(6,*) 'Ni : abinp(16)= ',abtsuji(16),
C     & ' AndGrv(28)= ',(cclog(28)+12)
C      write(6,*) 'Ti : abinp(17)= ',abtsuji(17),
C     & ' AndGrv(22)= ',(cclog(22)+12)
C      write(6,*) ' carbon/oxygen ratio = ',10.**(cclog(6)-cclog(8))
C
* normalization to H abundance added 11/03/93
        ccomp(1)=exp(cclog(1)/econst)
        do 9876 i=2,98
         ccomp(i)=exp(cclog(i)/econst)
         ccomp(i)=ccomp(i)/ccomp(1)
9876    continue
        ccomp(1)=ccomp(1)/ccomp(1)

C        write(6,*) 'elemental abundances in eqmol (jump=1):'
C        write(6,*)' elem(j),j,log10(ccomp(j))+12.,cclog(j)+12.,ccomp(j)'
        do 9878 jc=1,nmetal
        njc = NELEMX(jc)
C        write(6,9879) elem(jc),njc,log10(ccomp(njc))+12.,
C     *      cclog(njc)+12.,ccomp(njc)
9878    continue
9879    format(1x,a4,i3,2f12.3,1pe12.3)


c
c----lecture des molecules------------------------------------------------------
c
        J=0
        
        OPEN(UNIT=26,FILE=MOLEC,STATUS='OLD',readonly)
1010    J=J+1
        READ (26,5011) MOL(J),(C(J,K),K=1,5),MMAX(J),
     &    (NELEMM(M),NATOMM(M),M=1,4)

cC10H7   1.97800E+02-7.45200E+01-4.26000E+00 6.03100E-01-2.96000E-022 610  17  00 00
cmol(j)='C10H7   ',5 Pk coeff., mmax=#diff_at, 6=C,10,1=H,7~C10H7, 2next_poss_at=0
        
        MMAXJ=MMAX(J)
        IF(MMAXJ.EQ.0) GOTO 1014
C last line is empty, so also MMAXJ=0 when read is formatted
        DO 1012 M=1,MMAXJ
        NELEM(M,J)=NELEMM(M)
        NATOM(M,J)=NATOMM(M)
1012    CONTINUE
1110    GOTO 1010
C
1014    NMOL=J-1
        nmotsuji=nmol
        close(26)
        DO 1400 I=1,NMETAL
        NELEMI=NELEMX(I)
        P(NELEMI)=1.0E-20
1400    CONTINUE
        first=.false.
      endif

        p(99)= pein
        pesave=p(99)

c
c here we have only one depth point
        DO 1020 ITO=1,IT
        p(99)=pesave
* tt(ndp) model temperatures to which equilibrium has to be computed
* ppe(ndp) same for e- pressure
* pgg(ndp) same for gas pressure
* nh(ndp)  no need
        THETA=5040./tt
        TEM=tt
1024    PGLOG=log10(Pgas)
        PG=Pgas
        NNH=PG*AVO/(GRA*(1.+4.*AHE+SPA))
C       NNH=NH(ITO)
C       TO(ITO)=10.**T5L(ITO)
c
        CALL DIE(TEM,PG)
c


        PE=P(99)
        PELOG=log10(PE)
C
cccc    WRITE(6,6300)
cccc    WRITE(6,6091) PGLOG,PELOG,PE,THETA,TEM,Z(ITO)
cccc    WRITE(6,6301)
c
c----les atomes-----------------------------------------------------------------
c
        DO 1303 I=1,NMETAL
        NELEMI=NELEMX(I)
        FPLOG=log10(FP(NELEMI))
        XP(ITO,I)=P(NELEMI)+1.0E-20
        parptsuji(i)=xp(ito,i)
        PLOG=log10(XP(ITO,I))
        PDFPL=PLOG-FPLOG
c
cccc    WRITE (6,6302) ELEMNT(NELEMI),NELEMI,CCLOG(NELEMI),FP(NELEMI),
cccc     1  FPLOG,PLOG,PDFPL
cccc      PRESLOG(ITO,I)=PLOG
        IF(MOD(I,5))1303,1304,1303
1304    CONTINUE
1303    CONTINUE
c
cccc    WRITE(6,6307)
cccc    WRITE(6,5030) DD
cccc1231    WRITE (6,6091) PGLOG,PELOG,PE,THETA,TEM,Z(ITO)
c
cccc      EPRESLOG(ITO)=PELOG
c
cccc      WRITE(6,6992)
c
        IRL=120
c
c----les ions positifs---------------------------------------------------------
c
        DO 1184 I=1,NMETAL
        NELEMI=NELEMX(I)
        YA(I)=ELEMNT(NELEMI)
        PLOG= log10(P(NELEMI)+1.0D-99)
        KPLOG=log10(KP(NELEMI)+1.0D-99)
        YD(I)=KPLOG
        PIONL=PLOG+KPLOG-PELOG
        YB(I)=PIONL
cccc      PRESLOG(ITO,I+NMETAL)=YB(I)
        parptsuji(i+nmetal)=10**yb(i)
        XLOG=PIONL-PGLOG
        YC(I)=XLOG
C
        IF(I.NE.NMETAL) GOTO 1450
        IQ=I/120
        IR = NMETAL-IQ*120
        IRL =IR/3
        GOTO 1460
1450    IF(MOD(I,120)) 1184,1460,1184
C

1460    NBL=0
        DO 1470 K1=1,120,3
        NBL=NBL+1
        K2=K1+1
        K3=K1+2
        IF(NBL.EQ.IRL+1) GOTO 1480
c
cccc1475    WRITE(6,6495) YA(K1),YB(K1), YC(K1), YD(K1),
cccc     2             YA(K2),YB(K2), YC(K2), YD(K2),
cccc     3             YA(K3),YB(K3), YC(K3), YD(K3)
        IF ( MOD(NBL,5)) 1470,1500,1470
1500     CONTINUE
1470     CONTINUE
         GOTO 1184
1480    IRR = IR-IRL*3
        IF (IRR.EQ.0) GOTO 1184
        GOTO (1482,1484), IRR
c
1482    continue
cccc        WRITE(6,6495) YA(K1),YB(K1),YC(K1),YD(K1)
c
        GOTO 1184
c
1484    continue
cccc        WRITE(6,6495) YA(K1),YB(K1),YC(K1),YD(K1),
cccc     2             YA(K2),YB(K2),YC(K2),YD(K2)
1184    CONTINUE
C
        IRL =120
        KD=-119
c
c----les molecules-------------------------------------------------------------
c
	DO 1084 J=1,NMOL
ccc attention la double precision pour mol!!!	YA(J)=MOL(J)
	JCOUNT=JCOUNT+1
	PMOLL=log10(PPMOL(J)+1.0D-99)
	YB(J)=PMOLL
cccc	  PRESLOG(ITO,J+2*NMETAL)=YB(J)
        parptsuji(j+2*nmetal)=10**yb(j)
**************************************
        test=.false.
        if (test) then

        do jjj=1,mmax(j)
          if (nelem(jjj,j).eq.0) then
          print 1111,nelem(jjj,j),
     &         mol(j),parptsuji(j+2*nmetal)*natom(jjj,j)/
     &        (pgas*10**cclog(nelem(jjj,j))),
     &           log10(parptsuji(j+2*nmetal))
1111      format(i3,2x,a8,e10.2,x,f6.2)
         endif
        enddo

        endif
**************************************
        XLOG=PMOLL-PGLOG
        YC(J)=XLOG
        PPG(J,ITO)=XLOG
        YD(J)=APMLOG(J)
C
        IF(J.NE.NMOL) GOTO 2450
        IQ=J/120
        IR=NMOL-IQ*120
        IRL=IR/3
        GOTO 2460
2450    IF(MOD(J,120)) 2184,2460,2184
2460    NBL=0
c
cccc    WRITE(6,6092)
c
        KD=KD+120
        KF=KD+119
        DO 2470 K1=KD,KF,3
        NBL=NBL+1
        K2=K1+1
        K3=K1+2
        IF(NBL.EQ.IRL+1) GOTO 2480
c
cccc2475    WRITE(6,6496) YA(K1), YB(K1),YC(K1), YD(K1),
c
cccc     2             YA(K2),YB(K2),YC(K2),YD(K2),
cccc     3             YA(K3),YB(K3),YC(K3),YD(K3)
        IF(MOD(NBL,5)) 2470,2500,2470
2500    CONTINUE
2470    CONTINUE
        GOTO 2184
2480    IRR=IR-IRL*3
        IF(IRR.EQ.0) GOTO 2184
        GOTO (2482,2484),IRR
c
2482    continue
cccc        WRITE(6,6496) YA(K1),YB(K1),YC(K1),YD(K1)
c
        GOTO 2184
c
2484    continue
cccc        WRITE(6,6496) YA(K1),YB(K1),YC(K1),YD(K1),
cccc     2	       YA(K2),YB(K2),YC(K2),YD(K2)
c
2184     CONTINUE
1084    CONTINUE
1020    CONTINUE
c          DO I=1,4
c
cccc    WRITE(6,51) (XP(ITX,I),ITX=1,IT)
c
c          END DO
C   LES VALEURS SUIVANTES VIENNENT DE BB
        DO ITX=1,IT
         PPH(ITX)=XP(ITX,1)
          PPA(ITX)=XP(ITX,3)
           PB(ITX)=XP(ITX,4)
            PC13(ITX)=XP(ITX,5)
             PO(ITX)=XP(ITX,6)
              PTI(ITX)=XP(ITX,16)
        END DO
C
C-------ecriture des fichiers de resultats--------------------------------
C
C---1)  fichier pression des elements dans les couches atmospheriques-----
c
ccc		write(40)   (Z(I),I=1,IT)
ccc		write(40)   (5040./TETA(I),I=1,IT)
ccc		write(40)   (PGG(I),I=1,IT)
ccc		write(40)   (EPRESLOG(I),I=1,IT)
ccc	do I=1,2*NMETAL+NMOL
ccc		write(40) (PRESLOG(K,I),K=1,IT)
ccc	enddo
ccc	close(40)   
C
C---2)  fichier de reference: No/element----------------------------------
C
ccc	write(30,5031)  COM
ccc	write(30,5031)  NOM
ccc	write(30,500)   NMETAL,NMOL,IT
ccc		do I=1,NMETAL,7
ccc		  write(30,600) (J+4,ELEM(J),J=I,MIN(I+6,NMETAL))
ccc		enddo
ccc		do I=NMETAL+1,2*NMETAL,7
ccc		  write(30,620) (J+4,ELEM(J-NMETAL),J=I,
ccc     &             MIN(I+6,2*NMETAL))
ccc		enddo
ccc		do I=2*NMETAL+1,2*NMETAL+NMOL,7
ccc		  write(30,650) (J+4,MOL(J-2*NMETAL),J=I,MIN(I+6,
ccc     &             2*NMETAL+NMOL))
ccc		enddo
ccc	close(30)
C
C------------formats--------------------------------------------------------
c
  51    FORMAT(7E11.4)
 100    FORMAT (1X,'ENTRER UN TITRE ---> ',$)
 200    FORMAT (1X,'QUEL MODELE UTILISONS NOUS? ---> ',$)
 300    FORMAT (1X,'NOM DU FICHIER DE SORTIE ---> ',$)
 400    FORMAT (1X,'COMMENTAIRE? ---> ',$)
 430    FORMAT (1X,'FICHIER ATOMIQUE UTILISE? ---> ',$)
 460    FORMAT (1X,'FICHIER MOLECULAIRE UTILISE? ---> ',$)	
 501    FORMAT (I4)
 500    FORMAT (3(1X,I5))
 600    FORMAT (7(1X,I4,1X,A4,1X))
 620    FORMAT (7(1X,I4,1X,A4,'+',1X))
 650    FORMAT (7(1X,I4,1X,A8,1X))
5000    FORMAT (2I5,2F10.5,I10)
5001    FORMAT (A4,I6,F10.3,2I5,F10.5)
5011    FORMAT (A8,E12.5,4E12.5,I1,(I2,I3),3(I2,I2))
5021    FORMAT (F11.3,E12.5,E13.6)
5030    FORMAT (A)
5031    FORMAT (1X,A)
6031    FORMAT(1H1,20A4/)
6091    FORMAT (/,10X,'LOG PG=',F8.4,20X,'LOG PE=',F8.4,20X,'PE=',E13.6
     &  /,10X,'THETA =',F8.4,20X,'TEMP. =',F8.0,20X,'PROF. =',E14.6/)
6102    FORMAT (1H0, ' ELEMENT  ATOMIC NUMBER       I.P.        G(0)   '
     &,' G(1)    LOGN/NH')
6103    FORMAT(1H,5X,A4,8X,I5,3X,F10.3,5X,2I5,3X,F10.5)
6300    FORMAT(1H0,' EQUILIBRIUM PARTIAL PRESSURES OF THE GASEOUS',
     &    ' ATOMS'  ///)
6301    FORMAT(1H0,'ELEMENT',3X,'LOG N(E)/N(H)',4X,'P(E)',6X,'LOG P(E)'
     &,2X,'LOG P(A)',2X,'LOG P(A)/P(E)'/)
6302    FORMAT(1H,1X,A4,1X,I2,4X,F10.3,1X,E11.4,2F10.3,4X,F10.3)
6307    FORMAT(1H0,/' P(E) **** FICTITIOUS PRESSURE OF THE NUCLEUS',
     &' OF THE ELEMENT',/' P(A) ****PARTIAL PRESSURE OF THE',
     &' MONATOMIC GAS OF THE ELEMENT')
6092    FORMAT (1H1,3(5X,'MOLECULE   LOG P   LOG P/PG  LOG KP')//)
6495    FORMAT ( 3(8X,A4,'+',3F9.3))
6496    FORMAT ( 3(9X,A4,3F9.3))
6992    FORMAT (///,3(9X,'ION     LOG P   LOGP/PG  LOG KP ')//)
7000    FORMAT (21(A9,I3))
C
1100    return
        END



C DIE SUBROUTINE DE EQMOL

C
        SUBROUTINE DIE(TEM,PG)
        implicit real*16 (a-h,o-z)


        REAL*16 IP,KP
        double precision      MOL
        COMMON/CNEWC3 /NEWC3
        COMMON/COMFH1/C(600,5),NELEM(5,600),NATOM(5,600),
     2  MMAX(600),PPMOL(600),APMLOG(600),MOL(600),
     3 IP(100),CCOMP(100),UIIDUI(100),P(100),FP(100),KP(100),NELEMX(50),
     4  EPS,SWITER,NMETAL,NMOL,NIMAX
        DIMENSION FX(100),DFX(100),Z(100),PREV(100),WA(50)
        LOGICAL  CWRITE
        DATA  CWRITE /.TRUE./
        ECONST=4.342945E-1
        EPSDIE=5.0E-3
        T=5040.0/TEM
        PGLOG=log10(PG)
C
C    HEH=HELIUM/HYDROGENE RATIO BY NUMBER
        HEH=CCOMP(2)/CCOMP(1)
C
C    EVALUATION OF LOG KP(MOL)
        DO 1025 J=1,NMOL
        APLOGJ=C(J,5)
        DO 1026 K=1,4
        KM5=5-K
        APLOGJ=APLOGJ*T + C(J,KM5)
1026    CONTINUE
        APMLOG(J)=APLOGJ

1025    CONTINUE
        DHH=(((0.1196952E-02*T-0.2125713E-01)*T+0.1545253E+00)*T
     1  -0.5161452E+01)*T+0.1277356E+02
        DHH=EXP(DHH/ECONST)
C
C  EVALUATION OF THE IONIZATION CONSTANTS
        TEM25=TEM**2*SQRT(TEM)
        DO 1060 I=1,NMETAL
        NELEMI = NELEMX(I)
*
* calculation of the partition functions following Irwin (1981)
C in calls to partf() the dimension of tem and g0 (or g1) must be identical
        ndim = 1
        call partf(nelemi,1,tem,1,g0,ndim)
        call partf(nelemi,2,tem,1,g1,ndim)
        uiidui(nelemi)=g1/g0*0.6665
        
*
*****
        KP(NELEMI) =UIIDUI(NELEMI)*TEM25*EXP(-IP(NELEMI)*T/ECONST)
1060    CONTINUE
        HKP=KP(1)
        IF(T-0.6) 1084,1072,1072
C
C   PRELIMINARY VALUE OF PH AT HIGH TEMPERATURES
1084    PPH=SQRT(HKP*(PG/(1.0+HEH)+HKP))-HKP
        PH=PPH**2/HKP
        GOTO 1102
C
C   PRELIMINARY VALUE OF PH AT LOW TEMPERATURES
1072    IF(PG/DHH-0.1) 1073,1073,1074
1073    PH=PG/(1.0+HEH)
        GOTO 1102
1074    PH=0.5 * (SQRT(DHH*(DHH+4.0 *PG/(1.0+HEH)))-DHH)
C
C  EVALUATION OF THE FICTITIOUS PRESSURES OF HYDROGENE
C     PG=PH+PHH+2.0*PPH+HEH*(PH+2.0*PHH+PPH)
1102    U=(1.0+2.0*HEH)/DHH
        Q=1.0+HEH
        R=(2.0+HEH)*SQRT(HKP)
        S=-1.0*PG
        X=SQRT(PH)
        ITERAT=0
1103    F=((U*X**2+Q)*X+R)*X+S
        DF=2.0*(2.0*U*X**2+Q)*X+R
        XR=X-F/DF
        IF(ABS((X-XR)/XR)-EPSDIE) 1105,1105,1106
1106    ITERAT=ITERAT+1
***        print*,'iteration hydrogen ( iterat )  ',iterat
        IF(ITERAT-50) 1104,1104,1107
1107    continue
***        print*,'nach continue ( if )'
        WRITE(6,6108) TEM,PG,X,XR,PH
6108    FORMAT(1H1, ' NOT CONVERGE IN DIE '/// 'TEM=',F9.2,5X,'PG=',
     &  E12.5,5X,'X1=',E12.5,5X,'X2=',E12.5,5X,'PH=',E12.5/////)
        GOTO 1105
1104    X=XR
***        print*,'x=xr  ',x
        GOTO 1103
1105    PH=XR**2
        PHH=PH**2/DHH
        PPH=SQRT(HKP*PH)
        FPH=PH+2.0*PHH+PPH
cccc    WRITE(6,6109) TEM,T,PG,FPH,PH
6109    FORMAT(///,5H TEM=,F10.2,10X,7H THETA=,F8.4,10X,'PG=',E13.5,10X,
     2  'FPH=',E12.3,10X,'PH=',E12.3///)
C   P(100)=PH+
        P(100)=PPH
C
C   EVALUATION OF THE FICTITIOUS PRESSURE OF EACH ELEMENT
        DO 1070 I=1,NMETAL
        NELEMI=NELEMX(I)
        FP(NELEMI)=CCOMP(NELEMI)*FPH
cccc    print*,'elem,fp,fph ', nelemi,fp(nelemi),fph
1070    CONTINUE
C
C CHECK OF INITIALIZATION
        PE=P(99)
        IF(PH-P(1)) 1402,1402,1401
1401    DO 1403 I=1,NMETAL
        NELEMI=NELEMX(I)
        P(NELEMI)=FP(NELEMI)*EXP(-5.0*T/ECONST)
        if (P(nelemi).lt.1.e-20) then
         p(nelemi)=1.e-20
CCC      print*,'element  ',nelemi,' set to tractable pressure'
        endif
1403    CONTINUE
        P(1)=PH
C
C    RUSSELL EQUATIONS
1402    NITER = 0
1040    DO 1030 I=1,NMETAL
        NELEMI=NELEMX(I)
        FX(NELEMI)=-FP(NELEMI)+P(NELEMI)*(1.0+KP(NELEMI)/PE)
        DFX(NELEMI)=1.0+KP(NELEMI)/PE
1030    CONTINUE
* ??? ok???
        SPnION=0.0
        DO 1041 J=1,NMOL
        MMAXJ=MMAX(J)
        PMOLJL=-APMLOG(J)
        DO 1042 M=1,MMAXJ
        NELEMJ=NELEM(M,J)
        NATOMJ=NATOM(M,J)
        PMOLJL=PMOLJL+DFLOAT(NATOMJ)*log10(P(NELEMJ))
1042    CONTINUE

        IF(PMOLJL-(PGLOG+1.0)) 1046,1046,1047
1047    DO 1048 M=1,MMAXJ
        NELEMJ=NELEM(M,J)
        NATOMJ=NATOM(M,J)
        P(NELEMJ)=1.0E-2*P(NELEMJ)
        PMOLJL=PMOLJL+DFLOAT(NATOMJ)*(-2.0)
1048    CONTINUE
1046    PMOLJ=EXP(PMOLJL/ECONST)
        DO 1044 M=1,MMAXJ
        NELEMJ=NELEM(M,J)
        NATOMJ=NATOM(M,J)
        ATOMJ=DFLOAT(NATOMJ)
        IF(NELEMJ.EQ.99) SPNION=SPNION+PMOLJ
        DO 1043 I=1,NMETAL
        NELEMI=NELEMX(I)
        IF(NELEMJ.EQ.NELEMI) GO TO 1045
        GOTO 1043
1045    FX(NELEMI)=FX(NELEMI)+ATOMJ*PMOLJ
        DFX(NELEMI)=DFX(NELEMI)+ATOMJ**2*PMOLJ/P(NELEMI)
1043    CONTINUE
1044    CONTINUE
        PPMOL(J)=PMOLJ
1041    CONTINUE
C
C   SOLUTION OF THE RUSSELL EQUATIONS BY NEWTON-RAPHSON METHOD
        DO 2001 I=1,NMETAL
        NELEMI=NELEMX(I)
        WA(I)=log10(P(NELEMI)+1.0D-99)
2001    CONTINUE
        IMAXP1=NMETAL+1
        WA(IMAXP1)=log10(PE+1.0D-99)                  
        DELTA = 0.0 
        DO 1050 I=1,NMETAL
        NELEMI=NELEMX(I)
        PREV(NELEMI)=P(NELEMI)-FX(NELEMI)/DFX(NELEMI)
        PREV(NELEMI)=ABS(PREV(NELEMI))
        IF(PREV(NELEMI).LT.1.0E-20) PREV(NELEMI)=1.0E-20
        Z(NELEMI)=PREV(NELEMI)/P(NELEMI)
        DELTA=DELTA+ABS(Z(NELEMI)-1.0)
        IF(SWITER)  2500,2500,2501
2501    P(NELEMI)=(PREV(NELEMI)+P(NELEMI))*0.5
        GOTO 1050
2500    P(NELEMI)=PREV(NELEMI)
1050    CONTINUE
C
C   IONIZATION EQUILIBRIUM
        PEREV =0.0
        DO 1061 I=1,NMETAL
        NELEMI = NELEMX(I)
***        print*,'kp',kp(nelemi),'nel',nelemi,'P',p(nelemi),'perev',perev
        PEREV=PEREV+KP(NELEMI)*P(NELEMI)
1061    CONTINUE
        PEREV=SQRT(PEREV/(1.0+SPNION/PE))
        DELTA=DELTA+ABS((PE-PEREV)/PE)
        PE=(PEREV+PE)*0.5
ccc        print*,'iter,pe ',niter,pe
        P(99)=PE
        IF(DELTA-EPS) 1051,1051,1052
1052    NITER=NITER+1
***        print*,'iteration ioneqil.( niter ) ',niter
***        print*,'nimax ',nimax 
        
        IF(NITER-NIMAX) 1040,1040,1054
1054    WRITE(6,6055) NIMAX
        
        print*,'error ',delta,' convergence criterium ',eps
6055    FORMAT(1H0,'*DOES NOT CONVERGE AFTER ',I4,' ITERATIONS')
C ATOMC NUMBER 99= ELECTRON
1051    RETURN
        END


      subroutine test_Tsuji_jf(tt,pg)
      implicit real*16 (a-h,o-z)
* NOT THE ORIGINAL VERSION> CHANGED TO STUDY SOME SPECIES
* test of tsuji partial pressure

* version     6.9.94   Ch.Helling

      include 'parameter.inc'
      character*128 molec,met
      common/files/molec,met
      common /fullequilibrium/ partryck(ndp,maxmol),
     &  xmettryck(ndp,maxmet),xiontryck(ndp,maxmet),
     *  partp(ndp,0:maxmol)
      common /tsuji/ nattsuji,nmotsuji,parptsuji(500),abtsuji(17,ndp)
      common/cmolrat/fold(ndp,8),molold,kl
      common/statec/dum1(10*ndp),tauln(ndp),ro(ndp),ntau,iter


	common/jff/atms(300)

!        write(6,*)  ' input abundances again after transfer'
!	do 798 i=1,17
!	write(6,799) abtsuji(i)
!798     continue
!799     format(f12.3)

* source file tsuji_big.data

      met ='data/tsuji.atoms'
      molec ='data/comfits.data'


        pein=1
        temp=tt
        ttemp=5040./temp
        pgas=pg
!Trying to find out what temperature tt is:

      !open(unit=200, file='checktt.txt',status='old',position='append')
      !write(200,*) 'tt=', tt
      !close(200)

      !open(unit=201,file='checkTEFF.txt',status='old',position='append')
      !write(201,*) 'TEFF', TEFF
      !close(201)

      !open(unit=202,file='checkP.txt',status='old',position='append')
      !write(202,*) 'P=', pgas
      !close(202)

        call eqmol_jf(temp,pein,pgas,pe)

        call takemolec_jf(kl)

        pein=pe


100       format(f9.2,f6.3,f10.3,E10.3,3(1x,f7.2),7x,14(1x,f7.2))
200       format(A9,A6,A10,A10,18(1x,A7))


       return
      end


      subroutine takemolec_jf(kk)
      implicit real*16 (a-h,o-z)
      
*
* this routine is to be used after a call to jon,
*  -if eqmol has been called also-, in order to get the pressures
* one needs placed in the common 'fullequilibrium'.
* This is the routine to change if one wants more or/and other
* molecular pressure to be kept.
* The values in the common fullequilibrium can then be used for
* computation of opacities, printing, etc.
* kk is the depth point to be adressed.
* 020890 BPlez
*
* 13.09.94 Ch.Helling : change to be able to calculate all possible
*                       molecules in tsuji_big.data
*


      include 'parameter.inc'
      logical first
      character*5 name_mol, name_listmo
      common /fullequilibrium/ partryck(ndp,maxmol),
     &  xmettryck(ndp,maxmet),xiontryck(ndp,maxmet),
     *  partp(ndp,0:maxmol)
      common /tsuji/ nattsuji,nmotsuji,parptsuji(500),abtsuji(17,ndp)
      COMMON/CARC3/F1P,F3P,F4P,F5P,HNIC,PRESMO(33)
      common /cmolname/name_mol(maxmol),name_listmo(maxmol)
      data presmo(3)/0.d0/
      data first/.true./

	common/jff/atms(300)
      common/cmolrat/fold(ndp,8),molold,kl

* store the metals and ions
* they are not yet indexed like in jon.

*      print*,'nu in TAKEMOLEC' 
       if (maxmol.lt.342) stop 'takemolec: maxmol too small !'
      if (first) then
       open(unit=42,file='data/name_mol.dat'
     &    ,status='old')
       i=1
       do 442
       read(42,443,end=449) j,name_mol(i)
       i = i + 1
442    continue
449    continue
C       do 4431 j=1,181,15
C       write(6,4432) ((i+j),i=0,14)
C4431   write(6,443) (name_mol(i+j),i=0,14)
C       write(6,443) (name_mol(i),i=196,207)
443    format(2x,15a5)
C4432   format(15i5)
       close(42)
       first = .false.
      end if

      if(nattsuji.gt.maxmet) stop 'takemolec: maxmet too small'
      do 30 j=1,nattsuji
        xmettryck(kk,j)=parptsuji(j)
C       write(6,*) 'metal tryk ',j,log10(xmettryck(kk,j))
        xiontryck(kk,j)=parptsuji(j+nattsuji)
C       write(6,*) 'ion tryk ',j,log10(xiontryck(kk,j))
30    continue
*
*      print*,'nu after looking for right atoms'
*
* the indexes from 1 to 33 for partryck correspond to
* the ones in presmo (common carc3).

C     do 888 j=1,207
C        partryck(kk,j)=parptsuji(2*nattsuji+j)
!	 write(6,*) 'molekyl tryk ',log10(partryck(kk,j))
C888  continue

C       write(6,*) 'k,j,name,partryck(kk,j)='

       partryck(kk,1)=parptsuji(2*nattsuji+1)
       partryck(kk,2)=parptsuji(2*nattsuji+2)
       name_listmo(1) = name_mol(1)
       name_listmo(2) = name_mol(2)
        j=1
C        write(6,887) kk,j,name_listmo(1),partryck(kk,1)
        j=2
C        write(6,887) kk,j,name_listmo(2),partryck(kk,2)
* H2+ (not self-consistent)
!      partryck(kk,3)=presmo(3)
       partryck(kk,3)=1.D-99
       do 10 j=4,16
  	partryck(kk,j)=parptsuji(2*nattsuji+j-1)
        name_listmo(j) = name_mol(j-1)
C        write(6,887) kk,j,name_listmo(j),partryck(kk,j)
  10    continue
       partryck(kk,17)=1.D-99
       partryck(kk,20)=1.D-99   !C3H not in new list
       partryck(kk,36)=1.D-99   !CaH not in new list
       partryck(kk,37)=1.D-99   !LaO not in new list
       partryck(kk,69)=1.D-99   !TiS not in new list
       name_listmo(3) = 'H2+   '
       name_listmo(17) = '     '
       name_listmo(20) = 'C3H  '
       name_listmo(36) = '     '
       name_listmo(37) = 'LaO  '
       name_listmo(69) = 'TiS  '

      do 20 j=18,19
        partryck(kk,j)=parptsuji(2*nattsuji+j-2)
        partryck(kk,j)=max(partryck(kk,j),1.0d-99)
        name_listmo(j) = name_mol(j-2)
C        write(6,887) kk,j,name_listmo(j),partryck(kk,j)
20    continue
887    format(2i3,2x,a5,2x,1pe12.3)

      do 21 j=21,35
        partryck(kk,j)=parptsuji(2*nattsuji+j-3)
        partryck(kk,j)=max(partryck(kk,j),1.0d-99)
        name_listmo(j) = name_mol(j-3)
C        write(6,887) kk,j,name_listmo(j),partryck(kk,j)
21    continue

      do 22 j=38,68
        partryck(kk,j)=parptsuji(2*nattsuji+j-5)
        partryck(kk,j)=max(partryck(kk,j),1.0d-99)
        name_listmo(j) = name_mol(j-5)
22    continue

      do 23 j=70,nmotsuji+6
        partryck(kk,j)=parptsuji(2*nattsuji+j-6)
        partryck(kk,j)=max(partryck(kk,j),1.0d-99)
23    continue

      do 24 j=70,nmotsuji+6
        name_listmo(j) = name_mol(j-6)
24    continue


*****

      return
      end


*
      subroutine eqmol_jf(tt,pein,pgas,pe)
      implicit real*16 (a-h,o-z)
      include 'parameter.inc'
*     
c  ce programme calcule l'equilibre d'atomes et molecules dont la lise
c   est contenue dans des fichiers METAL.DAT et MOLE.DAT .
c                                                       b. plez   7/7/87
C                                                       modif Pe- 6/11/87
C                                                       modif    22/11/88
* for a change of delta(D0) of the dissociation energy (in eV) of a molecule,
* change the c(2) coefficient in log(Kp) expansion by -delta(D0)
*
* VERSION FOR INCLUSION IN THE MARCS CODE. Uses only one depth point at a time
*                    BPz 010890
* test version for determining minimum input file for molecules
* BPz 161190
* 
* Ch.Helling 070994

        character*128 MET,MOLEC
	character*8 name                           ! molecule name
	real*16 coeff0,coeff1,coeff2,coeff3,coeff4  ! logK fit coefficients
	integer nats,typ1,typ2,typ3,typ4           ! # atoms in all, types
	integer n1,n2,n3,n4		           ! # atoms of types 1-4
	integer natt                               ! # different atoms
        common/files/molec,met
        real*16       Z(51)  
ccc     real*4       EPRESLOG(50),PRESLOG(50,350)
        REAL*16 IP,KP,KPLOG,IPI,NNH,NHE
        doubleprecision   MOL
        COMMON/COMFH1/C(600,5),NELEM(5,600),NATOM(5,600),
     2   MMAX(600),PPMOL(600),APMLOG(600),MOL(600),
     3   IP(100),CCOMP(100,ndp),UIIDUI(100),P(100),FP(100),KP(100),
     4   NELEMX(50),EPS,SWITER,NMETAL,NMOL,NIMAX
        COMMON/COM6/TETAEF,GLOG,ASASOL,NHE,TIT(9)
        COMMON/VAL/PPG(600,50)
        DIMENSION YA(525),YB(525),YC(525),YD(525),ELEMNT(100),
     2 CCLOG(100,ndp),G0(100),G1(100),NATOMM(5),NELEMM(5),XP(50,100),
     3 PPA(51),PC13(51),PPH(51),PB(51),PO(51),PTI(51),ELEM(99)

* commons from injon
      COMMON/CI5/abmarcs(17,ndp),ANJON(17,5),Hz(5),PARTz(17,5),
     & DXIz,F1z,F2z,F3z,F4z,F5z,XKHMz,XMHz,XMYz(ndp)
C      real absc,abti,abv,abmn,abco
      common /auxabund/ absc,abti,abv,abmn,abco
* common for partial pressures
      common /tsuji/ nattsuji,nmotsuji,parptsuji(500),abtsuji(17,ndp)
      common/cmolrat/fold(ndp,8),molold,kl
      common/statec/dum1(10*ndp),tauln(ndp),ro(ndp),ntau,iter
      logical first,test
      DATA ELEMNT(99),ELEMNT(100)/'E-','H*'/
      data first/.true./

	common/jff/atms(300)

        IND=1
        ECONST=4.342945E-1
        AVO=0.602217E+24
        SPA=0.196E-01
        GRA=0.275423E+05
        AHE=0.100E+00


C
c----lecture de quelques renseignements dans FOR009.DAT------------------------
C 
ccc     READ (9,5000) NMETAL,NIMAX,EPS,SWITER,IT
ccc        nmetal=38
        nimax=3000
        eps=0.001
        switer=1
        it=1

! Iteration parameters


***
*  ! Nbre d'atomes,d'iterations,rcitere de convergence,?,Nbre couches atmospheriques
        IF(NIMAX.GT.0) GOTO 999
        NIMAX=50
        EPS=0.005
c
999     continue
cccc    WRITE(6,6102)
C
C----lecture du fichier contenant les atomes----------------------------------- 
C
        
      if (first) then
*
C     met ='/home/ast/uffegj/tsuji.atoms_AndersGrev'
C     molec ='/home/ast/uffegj/jens/comfits.data'

        OPEN(UNIT=26,FILE=MET,STATUS='OLD',readonly)
        rewind(26)
        read (26,*) nmetal
        write(6,*)
     &'Number of elements read into eqmol from tsuji.atoms: '
     &,nmetal
        nattsuji=nmetal
        DO 1002 I=1,NMETAL
        READ (26,5001) ELEMXI,NELEMI,IPI,IG0,IG1,CCLOGI 
*! nom(code),Nbre d'e-,pot d'ionisation,?,?,log de l'abondancce

        NELEMX(I)=NELEMI        ! Atomic number
        ELEM(I)=ELEMXI          ! Symbol
        ELEMNT(NELEMI)=ELEMXI   ! Symbol
        IP(NELEMI)=IPI          ! Ionization potential
        UIIDUI(NELEMI)=DFLOAT(IG1)*0.661/DFLOAT(IG0)
        G0(NELEMI)=IG0          ! ?
        G1(NELEMI)=IG1          ! ?
        CCLOG(NELEMI,1:ntau)=CCLOGI    ! Logarithm of abundance, H=0
!        CCOMP(NELEMI)=EXP(CCLOGI/ECONST)
ckeep for debug
C           WRITE(6,6103) ELEMXI,NELEMI,IPI,IG0,IG1,CCLOGI
c 
1002   CONTINUE

      close(26)


C
C Compare the elemental abundances read in from the MARCS input
C file to those read in from tsuji01.dat
C
  
* elements are in Tsuji called with their electron nummber

C     write(6,*) 'H  : abinp(1)= ',abtsuji(1),
C    & ' AndGrv(1)= ',(cclog(1)+12)
C     write(6,*) 'He : abinp(2)= ',abtsuji(2),
C    & ' AndGrv(2)= ',(cclog(2)+12)
C     write(6,*) 'C  : abinp(3)= ',abtsuji(3),
C    & ' AndGrv(6)= ',(cclog(6)+12)
C     write(6,*) 'N  : abinp(4)= ',abtsuji(4),
C    & ' AndGrv(7)= ',(cclog(7)+12)
C     write(6,*) 'O  : abinp(5)= ',abtsuji(5),
C    & ' AndGrv(8)= ',(cclog(8)+12)
C     write(6,*) 'Na : abinp(7)= ',abtsuji(7),
C    & ' AndGrv(11)= ',(cclog(11)+12)
C     write(6,*) 'Mg : abinp(8)= ',abtsuji(8),
C    & ' AndGrv(12)= ',(cclog(12)+12)
C     write(6,*) 'Al : abinp(9)= ',abtsuji(9),
C    & ' AndGrv(13)= ',(cclog(13)+12)
C     write(6,*) 'Si : abinp(10)= ',abtsuji(10),
C    & ' AndGrv(14)= ',(cclog(14)+12)
C     write(6,*) 'S  : abinp(11)= ',abtsuji(11),
C    & ' AndGrv(16)= ',(cclog(16)+12)
C     write(6,*) 'K  : abinp(12)= ',abtsuji(12),
C    & ' AndGrv(19)= ',(cclog(19)+12)
C     write(6,*) 'Ca : abinp(13)= ',abtsuji(13),
C    & ' AndGrv(20)= ',(cclog(20)+12)
C     write(6,*) 'Cr : abinp(14)= ',abtsuji(14),
C    & ' AndGrv(24)= ',(cclog(24)+12)
C     write(6,*) 'Fe : abinp(15)= ',abtsuji(15),
C    & ' AndGrv(26)= ',(cclog(26)+12)
C     write(6,*) 'Ni : abinp(16)= ',abtsuji(16),
C    & ' AndGrv(28)= ',(cclog(28)+12)
C     write(6,*) 'Ti : abinp(17)= ',abtsuji(17),
C    & ' AndGrv(22)= ',(cclog(22)+12)
    
!      rmag = 10.**(abtsuji(4,kl)-cclog(7)-12.)
!     &  + 10.**(abtsuji(5,kl)-cclog(8)-12.)
!     &  + 10.**(abtsuji(7,kl)-cclog(11)-12.)
!     &  + 10.**(abtsuji(8,kl)-cclog(12)-12.)
!     &  + 10.**(abtsuji(9,kl)-cclog(13)-12.)
!     &  + 10.**(abtsuji(10,kl)-cclog(14)-12.)
!     &  + 10.**(abtsuji(11,kl)-cclog(16)-12.)
!     &  + 10.**(abtsuji(12,kl)-cclog(19)-12.)
!     &  + 10.**(abtsuji(13,kl)-cclog(20)-12.)
!     &  + 10.**(abtsuji(14,kl)-cclog(24)-12.)
!     &  + 10.**(abtsuji(15,kl)-cclog(26)-12.)
!     &  + 10.**(abtsuji(16,kl)-cclog(28)-12.)
!     &  + 10.**(abtsuji(17,kl)-cclog(22)-12.)
!      rmaglog = log10(rmag/13.)
!      write(6,*)' log average ratio of Marcs-input div. tsuji(~solar) ='
!     &    ,rmaglog

C The elemental abundances read in from tsuji01.dat
C are set to those read in from the input file (AB 1995-05)

!       do 1102 jab = 3,98
!       cclog(jab) = cclog(jab) + rmaglog
!1102   continue
C       
       cclog(1,1:ntau) = abtsuji(1,1:ntau) - 12.
       cclog(2,1:ntau) = abtsuji(2,1:ntau) - 12.
       cclog(6,1:ntau) = abtsuji(3,1:ntau) - 12.
       cclog(7,1:ntau) = abtsuji(4,1:ntau) - 12.
       cclog(8,1:ntau) = abtsuji(5,1:ntau) - 12.
       cclog(11,1:ntau)= abtsuji(7,1:ntau) - 12.
       cclog(12,1:ntau)= abtsuji(8,1:ntau) - 12.
       cclog(13,1:ntau)= abtsuji(9,1:ntau) - 12.
       cclog(14,1:ntau)= abtsuji(10,1:ntau)- 12.
       cclog(16,1:ntau)= abtsuji(11,1:ntau)- 12.
       cclog(19,1:ntau)= abtsuji(12,1:ntau)- 12.
       cclog(20,1:ntau)= abtsuji(13,1:ntau)- 12.
       cclog(24,1:ntau)= abtsuji(14,1:ntau)- 12.
       cclog(26,1:ntau)= abtsuji(15,1:ntau)- 12.
       cclog(28,1:ntau)= abtsuji(16,1:ntau)- 12.
       cclog(22,1:ntau)= abtsuji(17,1:ntau)- 12.
*
C      write(6,*) ' control for right abundances:'
C      write(6,*) 'H  : abinp(1)= ',abtsuji(1),
C     & ' AndGrv(1)= ',(cclog(1)+12)
C      write(6,*) 'He : abinp(2)= ',abtsuji(2),
C     & ' AndGrv(2)= ',(cclog(2)+12)
C      write(6,*) 'C  : abinp(3)= ',abtsuji(3),
C     & ' AndGrv(6)= ',(cclog(6)+12)
C      write(6,*) 'N  : abinp(4)= ',abtsuji(4),
C     & ' AndGrv(7)= ',(cclog(7)+12)
C      write(6,*) 'O  : abinp(5)= ',abtsuji(5),
C     & ' AndGrv(8)= ',(cclog(8)+12)
C      write(6,*) 'Na : abinp(7)= ',abtsuji(7),
C     & ' AndGrv(11)= ',(cclog(11)+12)
C      write(6,*) 'Mg : abinp(8)= ',abtsuji(8),
C     & ' AndGrv(12)= ',(cclog(12)+12)
C      write(6,*) 'Al : abinp(9)= ',abtsuji(9),
C     & ' AndGrv(13)= ',(cclog(13)+12)
C      write(6,*) 'Si : abinp(10)= ',abtsuji(10),
C     & ' AndGrv(14)= ',(cclog(14)+12)
C      write(6,*) 'S  : abinp(11)= ',abtsuji(11),
C     & ' AndGrv(16)= ',(cclog(16)+12)
C      write(6,*) 'K  : abinp(12)= ',abtsuji(12),
C     & ' AndGrv(19)= ',(cclog(19)+12)
C      write(6,*) 'Ca : abinp(13)= ',abtsuji(13),
C     & ' AndGrv(20)= ',(cclog(20)+12)
C      write(6,*) 'Cr : abinp(14)= ',abtsuji(14),
C     & ' AndGrv(24)= ',(cclog(24)+12)
C      write(6,*) 'Fe : abinp(15)= ',abtsuji(15),
C     & ' AndGrv(26)= ',(cclog(26)+12)
C      write(6,*) 'Ni : abinp(16)= ',abtsuji(16),
C     & ' AndGrv(28)= ',(cclog(28)+12)
C      write(6,*) 'Ti : abinp(17)= ',abtsuji(17),
C     & ' AndGrv(22)= ',(cclog(22)+12)
C       write(6,*) ' carbon/oxygen ratio = ',10.**(cclog(6)-cclog(8))

* normalization to H abundance added 11/03/93
        ccomp(1,1:ntau)=exp(cclog(1,1:ntau)/econst)
        do 9876 i=2,98
         ccomp(i,1:ntau)=exp(cclog(i,1:ntau)/econst)
         ccomp(i,1:ntau)=ccomp(i,1:ntau)/ccomp(1,1:ntau)
9876    continue
        ccomp(1,1:ntau)=ccomp(1,1:ntau)/ccomp(1,1)

C        write(6,*)' abundances computed in eqmol_jf (jump=2):'
C        write(6,*)' elem(j),j,log10(ccomp(j))+12.,cclog(j)+12.,ccomp(j)'
        do 9878 jc=1,nmetal
        njc = NELEMX(jc)
C        write(6,9879) elem(jc),njc,log10(ccomp(njc))+12.,
C     *      cclog(njc)+12.,ccomp(njc)
9878    continue
9879    format(1x,a4,i3,2f12.3,1pe12.3)


c
c----lecture des molecules-----------------------------------------------------

! >>>
        J=0

        OPEN(UNIT=26,FILE=MOLEC,STATUS='OLD',readonly)
        open(unit=27,file='data/tsuji_comp.data'
     &      ,status='old',readonly)
	rewind(26)
	rewind(27)

1010    J=J+1
!        READ (26,5011) MOL(J),(C(J,K),K=1,5),MMAX(J),
!     &    (NELEMM(M),NATOMM(M),M=1,4)
*  !  nom(code),coeff. polynom.,?,4(?,?)
 
	read(26,5009,end=1015) mol(j),coeff0,coeff1,coeff2,coeff3,coeff4,nats
cwrite(6,5008) mol(j),coeff0,coeff1,coeff2,coeff3,coeff4,nats
1015	read(27,5010,end=1014) name,natt,typ1,n1,typ2,n2,typ3,n3,typ4,n4
c   	write(6,5010) name,natt,typ1,n1,typ2,n2,typ3,n3,typ4,n4
5008    format(2x,a8,2x,f10.3,4f15.4,i6)
!5009    format(a8,f8.3,4f13.4,i6)
5009    format(a8,5e15.5,i6)
5010    format(a8,i1,8i3)


Cmol(j)=name
        C(J,1)=coeff0
        C(J,2)=coeff1
        C(J,3)=coeff2
        C(J,4)=coeff3
        C(J,5)=coeff4
        mmax(j)=natt
	atms(j)=nats
        nelemm(1)=typ1
        nelemm(2)=typ2
        nelemm(3)=typ3
        nelemm(4)=typ4
        natomm(1)=n1
        natomm(2)=n2
        natomm(3)=n3
        natomm(4)=n4	
        
        MMAXJ=MMAX(J)
        IF(MMAXJ.EQ.0) GOTO 1014
        DO 1012 M=1,MMAXJ
        NELEM(M,J)=NELEMM(M)
        NATOM(M,J)=NATOMM(M)
1012    CONTINUE
1110    GOTO 1010
C
1014    NMOL=J-1
        nmotsuji=nmol
        close(26)
	close(27)
        DO 1400 I=1,NMETAL
        NELEMI=NELEMX(I)
        P(NELEMI)=1.0D-99
1400    CONTINUE
        first=.false.
      endif

        p(99)= pein
        pesave=p(99)

c
c here we have only one depth point
        DO 1020 ITO=1,IT
        p(99)=pesave
* tt(ndp) model temperatures to which equilibrium has to be computed
* ppe(ndp) same for e- pressure
* pgg(ndp) same for gas pressure
* nh(ndp)  no need
        THETA=5040./tt
        TEM=tt
1024    PGLOG=log10(Pgas)
        PG=Pgas
        NNH=PG*AVO/(GRA*(1.+4.*AHE+SPA))
C       NNH=NH(ITO)
C       TO(ITO)=10.**T5L(ITO)

c
        CALL DIE_jf(TEM,PG)
c


        PE=P(99)
        PELOG=log10(PE)
C
!        WRITE(6,6300)
!	WRITE(6,6091) PGLOG,PELOG,PE,THETA,TEM,Z(ITO)
!        WRITE(6,6301)
c
c----les atomes----------------------------------------------------------------
c
        DO 1303 I=1,NMETAL
        NELEMI=NELEMX(I)
        FPLOG=log10(FP(NELEMI))
        XP(ITO,I)=P(NELEMI)+1.0D-99
        parptsuji(i)=xp(ito,i)
        PLOG=log10(XP(ITO,I))
        PDFPL=PLOG-FPLOG
c
c       WRITE (6,6302) ELEMNT(NELEMI),NELEMI,CCLOG(NELEMI),FP(NELEMI),
c    1  FPLOG,PLOG,PDFPL
c        PRESLOG(ITO,I)=PLOG
        IF(MOD(I,5))1303,1304,1303
1304    CONTINUE
1303    CONTINUE
c
!        WRITE(6,6307)
!        WRITE(6,5030) DD
c1231    WRITE (6,6091) PGLOG,PELOG,PE,THETA,TEM,Z(ITO)
c
cccc      EPRESLOG(ITO)=PELOG
c
cccc      WRITE(6,6992)
c
        IRL=120
c
c----les ions positifs---------------------------------------------------------
c
        DO 1184 I=1,NMETAL
        NELEMI=NELEMX(I)
        YA(I)=ELEMNT(NELEMI)
        PLOG= log10(P(NELEMI)+1.0D-99)
        KPLOG=log10(KP(NELEMI)+1.0D-99)
        YD(I)=KPLOG
        PIONL=PLOG+KPLOG-PELOG
        YB(I)=PIONL
cccc      PRESLOG(ITO,I+NMETAL)=YB(I)
        parptsuji(i+nmetal)=10**yb(i)
!	write(6,*) (i+nmetal),log10(parptsuji(i+nmetal))
        XLOG=PIONL-PGLOG
        YC(I)=XLOG
C
        IF(I.NE.NMETAL) GOTO 1450
        IQ=I/120
        IR = NMETAL-IQ*120
        IRL =IR/3
        GOTO 1460
1450    IF(MOD(I,120)) 1184,1460,1184
C

1460    NBL=0
        DO 1470 K1=1,120,3
        NBL=NBL+1
        K2=K1+1
        K3=K1+2
        IF(NBL.EQ.IRL+1) GOTO 1480
c
cccc1475    WRITE(6,6495) YA(K1),YB(K1), YC(K1), YD(K1),
cccc     2             YA(K2),YB(K2), YC(K2), YD(K2),
cccc     3             YA(K3),YB(K3), YC(K3), YD(K3)
        IF ( MOD(NBL,5)) 1470,1500,1470
1500     CONTINUE
1470     CONTINUE
         GOTO 1184
1480    IRR = IR-IRL*3
        IF (IRR.EQ.0) GOTO 1184
        GOTO (1482,1484), IRR
c
1482    continue
cccc        WRITE(6,6495) YA(K1),YB(K1),YC(K1),YD(K1)
c
        GOTO 1184
c
1484    continue
cccc        WRITE(6,6495) YA(K1),YB(K1),YC(K1),YD(K1),
cccc     2             YA(K2),YB(K2),YC(K2),YD(K2)
1184    CONTINUE
C
        IRL =120
        KD=-119
c
c----les molecules-------------------------------------------------------------
c
	DO 1084 J=1,NMOL
ccc attention la double precision pour mol!!!	YA(J)=MOL(J)
	JCOUNT=JCOUNT+1
	PMOLL=log10(PPMOL(J)+1.0D-99)
	YB(J)=PMOLL
cccc	  PRESLOG(ITO,J+2*NMETAL)=YB(J)
        parptsuji(j+2*nmetal)=10**yb(j)
**************************************
C        test=.false.
C        if (test) then
C
C        do jjj=1,mmax(j)
C          if (nelem(jjj,j).eq.0) then
C          print 1111,nelem(jjj,j),
C     &         mol(j),parptsuji(j+2*nmetal)*natom(jjj,j)/
C     &        (pgas*10**cclog(nelem(jjj,j))),
C     &           log10(parptsuji(j+2*nmetal))
C1111      format(i3,2x,a8,e10.2,x,f6.2)
C         endif
C        enddo
C
C        endif
**************************************
        XLOG=PMOLL-PGLOG
        YC(J)=XLOG
        PPG(J,ITO)=XLOG
        YD(J)=APMLOG(J)
C
        IF(J.NE.NMOL) GOTO 2450
        IQ=J/120
        IR=NMOL-IQ*120
        IRL=IR/3
        GOTO 2460
2450    IF(MOD(J,120)) 2184,2460,2184
2460    NBL=0
c
cccc    WRITE(6,6092)
c
        KD=KD+120
        KF=KD+119
        DO 2470 K1=KD,KF,3
        NBL=NBL+1
        K2=K1+1
        K3=K1+2
        IF(NBL.EQ.IRL+1) GOTO 2480
c
cccc2475    WRITE(6,6496) YA(K1), YB(K1),YC(K1), YD(K1),
c
cccc     2             YA(K2),YB(K2),YC(K2),YD(K2),
cccc     3             YA(K3),YB(K3),YC(K3),YD(K3)
        IF(MOD(NBL,5)) 2470,2500,2470
2500    CONTINUE
2470    CONTINUE
        GOTO 2184
2480    IRR=IR-IRL*3
        IF(IRR.EQ.0) GOTO 2184
        GOTO (2482,2484),IRR
c
2482    continue
cccc        WRITE(6,6496) YA(K1),YB(K1),YC(K1),YD(K1)
c
        GOTO 2184
c
2484    continue
cccc        WRITE(6,6496) YA(K1),YB(K1),YC(K1),YD(K1),
cccc     2	       YA(K2),YB(K2),YC(K2),YD(K2)
c
2184     CONTINUE
1084    CONTINUE
1020    CONTINUE
c          DO I=1,4
c
cccc    WRITE(6,51) (XP(ITX,I),ITX=1,IT)
c
c          END DO
C   LES VALEURS SUIVANTES VIENNENT DE BB
        DO ITX=1,IT
         PPH(ITX)=XP(ITX,1)
          PPA(ITX)=XP(ITX,3)
           PB(ITX)=XP(ITX,4)
            PC13(ITX)=XP(ITX,5)
             PO(ITX)=XP(ITX,6)
              PTI(ITX)=XP(ITX,16)
        END DO
C
C-------ecriture des fichiers de resultats--------------------------------
C
C---1)  fichier pression des elements dans les couches atmospheriques-----
c
ccc		write(40)   (Z(I),I=1,IT)
ccc		write(40)   (5040./TETA(I),I=1,IT)
ccc		write(40)   (PGG(I),I=1,IT)
ccc		write(40)   (EPRESLOG(I),I=1,IT)
ccc	do I=1,2*NMETAL+NMOL
ccc		write(40) (PRESLOG(K,I),K=1,IT)
ccc	enddo
ccc	close(40)   
C
C---2)  fichier de reference: No/element----------------------------------
C
ccc	write(30,5031)  COM
ccc	write(30,5031)  NOM
ccc	write(30,500)   NMETAL,NMOL,IT
ccc		do I=1,NMETAL,7
ccc		  write(30,600) (J+4,ELEM(J),J=I,MIN(I+6,NMETAL))
ccc		enddo
ccc		do I=NMETAL+1,2*NMETAL,7
ccc		  write(30,620) (J+4,ELEM(J-NMETAL),J=I,
ccc     &             MIN(I+6,2*NMETAL))
ccc		enddo
ccc		do I=2*NMETAL+1,2*NMETAL+NMOL,7
ccc		  write(30,650) (J+4,MOL(J-2*NMETAL),J=I,MIN(I+6,
ccc     &             2*NMETAL+NMOL))
ccc		enddo
ccc	close(30)
C
C------------formats--------------------------------------------------------
c
  51    FORMAT(7E11.4)
 100    FORMAT (1X,'ENTRER UN TITRE ---> ',$)
 200    FORMAT (1X,'QUEL MODELE UTILISONS NOUS? ---> ',$)
 300    FORMAT (1X,'NOM DU FICHIER DE SORTIE ---> ',$)
 400    FORMAT (1X,'COMMENTAIRE? ---> ',$)
 430    FORMAT (1X,'FICHIER ATOMIQUE UTILISE? ---> ',$)
 460    FORMAT (1X,'FICHIER MOLECULAIRE UTILISE? ---> ',$)	
 501    FORMAT (I4)
 500    FORMAT (3(1X,I5))
 600    FORMAT (7(1X,I4,1X,A4,1X))
 620    FORMAT (7(1X,I4,1X,A4,'+',1X))
 650    FORMAT (7(1X,I4,1X,A8,1X))
5000    FORMAT (2I5,2F10.5,I10)
5001    FORMAT (A4,I6,F10.3,2I5,F10.5)
5011    FORMAT (A8,E12.5,4E12.5,I1,(I2,I3),3(I2,I2))
5021    FORMAT (F10.3,E12.5,E13.6)
5030    FORMAT (A)
5031    FORMAT (1X,A)
6031    FORMAT(1H1,20A4/)
6091    FORMAT (/,10X,'LOG PG=',F8.4,20X,'LOG PE=',F8.4,20X,'PE=',E13.6
     &  /,10X,'THETA =',F8.4,20X,'TEMP. =',F8.0,20X,'PROF. =',E14.6/)
6102    FORMAT (1H0, ' ELEMENT  ATOMIC NUMBER       I.P.        G(0)   '
     &,' G(1)    LOGN/NH')
6103    FORMAT(1H,5X,A4,8X,I5,3X,F10.3,5X,2I5,3X,F10.5)
6300    FORMAT(1H0,' EQUILIBRIUM PARTIAL PRESSURES OF THE GASEOUS',
     &    ' ATOMS'  ///)
6301    FORMAT(1H0,'ELEMENT',3X,'LOG N(E)/N(H)',4X,'P(E)',6X,'LOG P(E)'
     &,2X,'LOG P(A)',2X,'LOG P(A)/P(E)'/)
6302    FORMAT(1H,1X,A4,1X,I2,4X,F10.3,1X,E11.4,2F10.3,4X,F10.3)
6307    FORMAT(1H0,/' P(E) **** FICTITIOUS PRESSURE OF THE NUCLEUS',
     &' OF THE ELEMENT',/' P(A) ****PARTIAL PRESSURE OF THE',
     &' MONATOMIC GAS OF THE ELEMENT')
6092    FORMAT (1H1,3(5X,'MOLECULE   LOG P   LOG P/PG  LOG KP')//)
6495    FORMAT ( 3(8X,A4,'+',3F9.3))
6496    FORMAT ( 3(9X,A4,3F9.3))
6992    FORMAT (///,3(9X,'ION     LOG P   LOGP/PG  LOG KP ')//)
7000    FORMAT (21(A9,I3))
C

C       write(6,*) ' there were nmetal, nmol, 2*nmet+nmol:'
C    &     ,nmetal,nmol,2*nmetal+nmol

1100    return
        END



C DIE SUBROUTINE DE EQMOL

C
        SUBROUTINE DIE_jf(TEM,PG)
        implicit real*16 (a-h,o-z)
        include 'parameter.inc'

        REAL*16 IP,KP
        double precision      MOL
        COMMON/COMFH1/C(600,5),NELEM(5,600),NATOM(5,600),
     2  MMAX(600),PPMOL(600),APMLOG(600),MOL(600),
     3  IP(100),CCOMP(100,ndp),UIIDUI(100),P(100),FP(100),KP(100),
     4  NELEMX(50),EPS,SWITER,NMETAL,NMOL,NIMAX
        DIMENSION FX(100),DFX(100),Z(100),PREV(100),WA(50)
        DIMENSION prev_patom(100),prev_PPMOL(600)
        common/cmolrat/fold(ndp,8),molold,kl
	common/jff/atms(300)

        ECONST=4.342945E-1
        EPSDIE=5.0E-3

        T=5040.0/TEM

        PGLOG=log10(PG)
C
C    HEH=HELIUM/HYDROGENE RATIO BY NUMBER
        HEH=CCOMP(2,kl)/CCOMP(1,kl)

C
!+++

C    EVALUATION OF LOG KP(MOL)
!	write(6,*) 'Temperature:', tem
!	write(6,*) 'Gibbs Energy:'
        DO 1025 J=1,NMOL
        APLOGJ=C(J,5)
        DO 1026 K=1,4
        KM5=5-K
        APLOGJ=APLOGJ*(TEM/1000.) + C(J,KM5)      ! dG from fit
1026    CONTINUE
!	write(6,*) j,(1000.*aplogj)
	APLOGJ=1000.*APLOGJ/(8.31451*TEM)         ! B=log10(EXP(dG/RT))
	APLOGJ=0.43429448*APLOGJ
!	write(6,*) 'atms(j)=',atms(j)
	APLOGJ=APLOGJ+(6.*(atms(j)-1.))            ! -6(n-1) pressure
	APMLOG(J)=APLOGJ
cwrite(6,*) 'apmlog(',j,')', apmlog(j), ' from comfits.data'

************************************************************************
***        if (j.eq.29) then 
***        blabla= 11.4047+log10(t)* (-1.1484+log10(t)*(0.6478+
***     &          log10(t)*(-0.6737))) -t*6.92  +1.
**** +1. is conversion SI to cgs for pressure
***        print *,' tsuji:', aplogj,' Sauval et Tatum:', blabla
***ccc        apmlog(j)=blabla
***      endif
************************************************************************

*       

1025    CONTINUE
        DHH=(((0.1196952E-02*T-0.2125713E-01)*T+0.1545253E+00)*T
     1  -0.5161452E+01)*T+0.1277356E+02
        DHH=EXP(DHH/ECONST)

C
C  EVALUATION OF THE IONIZATION CONSTANTS
        TEM25=TEM**2*SQRT(TEM)
        DO 1060 I=1,NMETAL
        NELEMI = NELEMX(I)

*
* calculation of the partition functions following Irwin (1981)
C in calls to partf() the dimension of tem and g0 (or g1) must be identical
        ndim = 1
        call partf(nelemi,1,tem,1,g0,ndim)
        call partf(nelemi,2,tem,1,g1,ndim)
        uiidui(nelemi)=g1/g0*0.6665
        
*
*****
        KP(NELEMI) =UIIDUI(NELEMI)*TEM25*EXP(-IP(NELEMI)*T/ECONST)
1060    CONTINUE
        HKP=KP(1)
        IF(T-0.6) 1084,1072,1072
C
C   PRELIMINARY VALUE OF PH AT HIGH TEMPERATURES
1084    PPH=SQRT(HKP*(PG/(1.0+HEH)+HKP))-HKP
        PH=PPH**2/HKP
        GOTO 1102
C
C   PRELIMINARY VALUE OF PH AT LOW TEMPERATURES
1072    IF(PG/DHH-0.1) 1073,1073,1074
1073    PH=PG/(1.0+HEH)
        GOTO 1102
1074    PH=0.5 * (SQRT(DHH*(DHH+4.0 *PG/(1.0+HEH)))-DHH)

C
C  EVALUATION OF THE FICTITIOUS PRESSURES OF HYDROGENE
C     PG=PH+PHH+2.0*PPH+HEH*(PH+2.0*PHH+PPH)
1102    U=(1.0+2.0*HEH)/DHH
        Q=1.0+HEH
        R=(2.0+HEH)*SQRT(HKP)
        S=-1.0*PG
        X=SQRT(PH)
        ITERAT=0
1103    F=((U*X**2+Q)*X+R)*X+S
        DF=2.0*(2.0*U*X**2+Q)*X+R
        XR=X-F/DF
        IF(ABS((X-XR)/XR)-EPSDIE) 1105,1105,1106
1106    ITERAT=ITERAT+1
***        print*,'iteration hydrogen ( iterat )  ',iterat
        IF(ITERAT-50) 1104,1104,1107
1107    continue
***        print*,'nach continue ( if )'
!        WRITE(6,6108) TEM,PG,X,XR,PH
6108    FORMAT(1H1, ' NOT CONVERGE IN DIE '/// 'TEM=',F9.2,5X,'PG=',
     &  E12.5,5X,'X1=',E12.5,5X,'X2=',E12.5,5X,'PH=',E12.5/////)
        GOTO 1105
1104    X=XR
***        print*,'x=xr  ',x
        GOTO 1103
1105    PH=XR**2
        PHH=PH**2/DHH
        PPH=SQRT(HKP*PH)
        FPH=PH+2.0*PHH+PPH
!        WRITE(6,6109) TEM,T,PG,FPH,PH
6109    FORMAT(///,5H TEM=,F10.2,10X,7H THETA=,F8.4,10X,'PG=',E13.5,10X,
     2  'FPH=',E12.3,10X,'PH=',E12.3///)
C   P(100)=PH+
        P(100)=PPH
C
C   EVALUATION OF THE FICTITIOUS PRESSURE OF EACH ELEMENT
        DO 1070 I=1,NMETAL
        NELEMI=NELEMX(I)
        FP(NELEMI)=CCOMP(NELEMI,kl)*FPH
!        print*,'elem,fp,fph ', nelemi,fp(nelemi),fph
1070    CONTINUE
C
C CHECK OF INITIALIZATION
        PE=P(99)
        IF(PH-P(1)) 1402,1402,1401
1401    DO 1403 I=1,NMETAL
        NELEMI=NELEMX(I)
        P(NELEMI)=FP(NELEMI)*EXP(-5.0*T/ECONST)
        if (P(nelemi).lt.1.D-99) then
         p(nelemi)=1.D-99
CCC      print*,'element  ',nelemi,' set to tractable pressure'
        endif
1403    CONTINUE
        P(1)=PH


        sumatoms = 0.
        summol =  0.
        DO 10504 I=1,NMETAL
        NELEMI=NELEMX(I)
        prev_patom(nelemi) = 0.
10504   CONTINUE
        DO 10505 J=1,NMOL
        prev_PPMOL(J) = 0.
10505   CONTINUE
C        write(6,*) ' it mol atm pe maxm maxa m+at/Pg:'
C
C    RUSSELL EQUATIONS
1402    NITER = 0
1040    DO 1030 I=1,NMETAL
        NELEMI=NELEMX(I)
        FX(NELEMI)=-FP(NELEMI)+P(NELEMI)*(1.0+KP(NELEMI)/PE)
        DFX(NELEMI)=1.0+KP(NELEMI)/PE
1030    CONTINUE
* ??? ok???
        SPnION=0.0
        prev_sumatoms = sumatoms
        prev_summol = summol
        DO 10502 I=1,NMETAL
        NELEMI=NELEMX(I)
        prev_patom(nelemi) = P(NELEMI)
10502   CONTINUE
        DO 10503 J=1,NMOL
        prev_PPMOL(J) = PPMOL(J)
10503   CONTINUE
        prev_pe = pe
C
C   IONIZATION EQUILIBRIUM
        summol = 0.
        amax_mol_delta = 0.
        bmax_mol_delta = 0.
        DO 1041 J=1,NMOL
        MMAXJ=MMAX(J)
        PMOLJL=-APMLOG(J)
        DO 1042 M=1,MMAXJ
        NELEMJ=NELEM(M,J)
        NATOMJ=NATOM(M,J)
        PMOLJL=PMOLJL+DFLOAT(NATOMJ)*log10(P(NELEMJ))
1042    CONTINUE

        IF(PMOLJL-(PGLOG+1.0)) 1046,1046,1047
1047    DO 1048 M=1,MMAXJ
        NELEMJ=NELEM(M,J)
        NATOMJ=NATOM(M,J)
        P(NELEMJ)=1.0E-2*P(NELEMJ)
        PMOLJL=PMOLJL+DFLOAT(NATOMJ)*(-2.0)
1048    CONTINUE
1046    PMOLJ=EXP(PMOLJL/ECONST)
        DO 1044 M=1,MMAXJ
        NELEMJ=NELEM(M,J)
        NATOMJ=NATOM(M,J)
        ATOMJ=DFLOAT(NATOMJ)
        IF(NELEMJ.EQ.99) SPNION=SPNION+PMOLJ
        DO 1043 I=1,NMETAL
        NELEMI=NELEMX(I)
        IF(NELEMJ.EQ.NELEMI) GO TO 1045
        GOTO 1043
1045    FX(NELEMI)=FX(NELEMI)+ATOMJ*PMOLJ
        DFX(NELEMI)=DFX(NELEMI)+ATOMJ**2*PMOLJ/P(NELEMI)
1043    CONTINUE
1044    CONTINUE
        PPMOL(J)=PMOLJ
        summol = summol + ppmol(j)
        amol_delta = abs( 1.-prev_ppmol(j)/ppmol(j) )
        amax_mol_delta = max ( amax_mol_delta, amol_delta )
        if(ppmol(j)/pg.gt.1.e-9) 
     &       bmax_mol_delta = max ( bmax_mol_delta, amol_delta )
1041    CONTINUE
C
C   SOLUTION OF THE RUSSELL EQUATIONS BY NEWTON-RAPHSON METHOD
        DO 2001 I=1,NMETAL
        NELEMI=NELEMX(I)
        WA(I)=log10(P(NELEMI)+1.0D-99)
2001    CONTINUE
        IMAXP1=NMETAL+1
        WA(IMAXP1)=log10(PE+1.0D-99)                  
        DELTA = 0.0 
        DO 1050 I=1,NMETAL
        NELEMI=NELEMX(I)
        PREV(NELEMI)=P(NELEMI)-FX(NELEMI)/DFX(NELEMI)
        PREV(NELEMI)=ABS(PREV(NELEMI))
        IF(PREV(NELEMI).LT.1.0D-99) PREV(NELEMI)=1.0D-99
        Z(NELEMI)=PREV(NELEMI)/P(NELEMI)
        DELTA=DELTA+ABS(Z(NELEMI)-1.0)
        IF(SWITER)  2500,2500,2501
2501    P(NELEMI)=(PREV(NELEMI)+P(NELEMI))*0.5
        GOTO 1050
2500    P(NELEMI)=PREV(NELEMI)
1050    CONTINUE

        sumatoms = 0.
        amax_atom_delta = 0.
        bmax_atom_delta = 0.
        DO 10501 I=1,NMETAL
        NELEMI=NELEMX(I)
        sumatoms = sumatoms + P(NELEMI)
        atom_delta = abs( 1.-prev_patom(nelemi)/P(NELEMI) )
        amax_atom_delta = max( atom_delta, amax_atom_delta)
        if(p(nelemi)/pg.gt.1.e-9) 
     &       bmax_atom_delta = max( atom_delta, bmax_atom_delta)
10501   CONTINUE
C
C   IONIZATION EQUILIBRIUM
        PEREV =0.0
        DO 1061 I=1,NMETAL
        NELEMI = NELEMX(I)
***        print*,'kp',kp(nelemi),'nel',nelemi,'P',p(nelemi),'perev',perev
        PEREV=PEREV+KP(NELEMI)*P(NELEMI)
1061    CONTINUE
        PEREV=SQRT(PEREV/(1.0+SPNION/PE))
        DELTA=DELTA+ABS((PE-PEREV)/PE)
        PE=(PEREV+PE)*0.5
ccc        print*,'iter,pe ',niter,pe
        P(99)=PE

C        if( abs(1.-(summol+sumatoms)/pg).le.0.002  .and.
C     &     abs(bmax_mol_delta)+abs(bmax_atom_delta).le.0.002) then
C        write(6,82) niter,summol,sumatoms,pe, bmax_mol_delta,
C     &    bmax_atom_delta,(summol+sumatoms)/pg
C82      format(i4,1p6e9.2)
C82      format(' it mol atm pe maxm maxa m+at/Pg',i4,1p6e9.2)
C        end if

        IF(DELTA-EPS) 1051,1051,1052
1052    NITER=NITER+1
***        print*,'iteration ioneqil.( niter ) ',niter
***        print*,'nimax ',nimax 
        
        IF(NITER-NIMAX) 1040,1040,1054
1054    WRITE(6,6055) NIMAX
        
        print*,'error ',delta,' convergence criterium ',eps
6055    FORMAT(1H0,'*DOES NOT CONVERGE AFTER ',I4,' ITERATIONS')
C ATOMC NUMBER 99= ELECTRON


1051    CONTINUE

C        write(6,80) tem,pg
C80      format(' temp, gas-pressure in input to die:',f8.1,1pe13.5)
C        write(6,*) ' sum-p-molecules, sum-p-atoms, pe, sum-these-2/pg:'
C        write(6,81) summol,sumatoms,pe,(summol+sumatoms)/pg
C81      format(1p4e13.5)
C        write(6,83) 
C     &           abs((sumatoms-prev_sumatoms)/prev_sumatoms),
C     &           abs((summol-prev_summol)/prev_summol),
C     &           abs((pe-prev_pe)/prev_pe)
C83      format(' d_sumatoms/prev,d_summol/prev,d_pe/prev:',1p3e12.4)
        amax_atom_delta = 0.
        bmax_atom_delta = 0.
        cmax_atom_delta = 0.
        dmax_atom_delta = 0.
        DO 87 I=1,NMETAL
        NELEMI=NELEMX(I)
        atom_delta = abs( 1.-prev_patom(nelemi)/P(NELEMI) )
        amax_atom_delta = max( atom_delta, amax_atom_delta)
        if(p(nelemi)/pg.gt.1.e-7) 
     &       bmax_atom_delta = max( atom_delta, bmax_atom_delta)
        if(p(nelemi)/pg.gt.1.e-9) 
     &       cmax_atom_delta = max( atom_delta, cmax_atom_delta)
        if(p(nelemi)/pg.gt.1.e-11) 
     &       dmax_atom_delta = max( atom_delta, dmax_atom_delta)
87      CONTINUE
        amax_mol_delta = 0.
        bmax_mol_delta = 0.
        cmax_mol_delta = 0.
        dmax_mol_delta = 0.
        DO 88 J=1,NMOL
        amol_delta = abs( 1.-prev_ppmol(j)/ppmol(j) )
        amax_mol_delta = max ( amax_mol_delta, amol_delta )
        if(ppmol(j)/pg.gt.1.e-7) 
     &       bmax_mol_delta = max ( bmax_mol_delta, amol_delta )
        if(ppmol(j)/pg.gt.1.e-9) 
     &       cmax_mol_delta = max ( cmax_mol_delta, amol_delta )
        if(ppmol(j)/pg.gt.1.e-11) 
     &       dmax_mol_delta = max ( dmax_mol_delta, amol_delta )
88      CONTINUE
C        write(6,86)amax_mol_delta, amax_atom_delta
C86      format(' max_mol_delta, max_atom_delta: ',1p2e12.4)
C        write(6,90)bmax_mol_delta, bmax_atom_delta
C90      format(' ditto for pp contributing > 1.e-7 pg: ',1p2e12.4)
C        write(6,91)cmax_mol_delta, cmax_atom_delta
C91      format(' ditto for pp contributing > 1.e-9 pg: ',1p2e12.4)
C        write(6,92)dmax_mol_delta, dmax_atom_delta
C92      format(' ditto for pp contributing > 1.e-11 pg: ',1p2e12.4)

        RETURN
        END

C
C------------------------------------------------------------
C                                                           I
      SUBROUTINE GEM_INIT
C                                                           I
C reads the names and indiexes etc to be used in Gibbs      I
C minimum energy routines                                   I
C------------------------------------------------------------
      implicit real*16 (a-h,o-z)
      include 'parameter.inc'

      dimension abund(natms),isubst(natms)
      COMMON/CI5/ABUND_marcs(17),ANJON(17,5),H(5),PART(17,5),DXI,
     &   F1,F2,F3,F4,F5,XKHM,XMH,XMY
C  the 17 marcs_abundances: H  HE C  N  O  NE NA MG AL SI S  K  CA CR FE NI Ti
C  and their corresponding GEM numbers, mx_elmi, of the 48 elements in the input 
C  array abund(i), or the stored (normalized) totab(i), or identical 
C  totabk(i),i=2,natms; natms=49
      dimension mx_elm(17)
      data mx_elm /2,3,7,8,9,11,12,13,14,15,17,20,21,24,26,28,22/
C                 H HE C N O NE NA MG AL SI S  K  CA CR FE NI Ti
      COMMON/cabtsuji/abtsuji(17)

      common /cabink/abink(ndp,nspec)
      common /cabundinp/abundatms_inp(natms)
      common /ctotabk/totabk(ndp,natms)
      character name_gem*8
      common /cgemnames/natms_gem,nions_gem,nspec_gem,name_gem(nspec)
      character sunz*1
      NAMELIST /ABUNDANCES/sunz,zscale,abundatms_inp
C atms,ions,spec ~ highest index of neutral atoms, ions, species total
       natms_gem = natms
       nions_gem = 127
       nspec_gem = nspec
        open(unit=2,file='data/gfits.data'
     &   ,status='old',readonly)   ! extract molecular etc names from here
        do 2125 i=1,nspec
            read(2,766) name_gem(i)
            read(2,*) ajunk1,ajunk2,ajunk3,ajunk4,ajunk5   ! Fit data
2125    continue
766     format(a8)
        close(unit=2)

       
        open(unit=2,file='data/elabund.dat'
     &   ,status='old',readonly)   ! read the elemental abundances
        do 125 i=2,natms
            read(2,*,end=126) kelem, abund(i)
125     continue
126	continue
        write(6,131) abund(2),abund(3),(abund(i),i=7,9),abund(26)
131     format('gem_init solar inp.abund(i),i=H,He,C,N,O,Fe:',/8f7.2)

C if sunz = yes, we use solar abundances from unit 2
        do 300 i=1,natms
300     abundatms_inp(i) =  999.9

        READ(5,ABUNDANCES)
        if(sunz.eq.'y'.or.sunz.eq.'Y') then
          write(6,*) 
     &    'We adopted solar abundances from unit=2=elabund.dat'
          go to 309
        end if

        if(zscale.ne.1.) then
          do 310 i=4,natms
          abund(i) = log10(zscale) + abund(i)    !do not scale 1,2,3 ~ e-,H,He
310       continue
        end if

          nchange_ab = 0
        do 320 i=1,natms
        if(abundatms_inp(i).ne.999.9) then
          abund(i) = abundatms_inp(i)
          nchange_ab = nchange_ab + 1
          isubst(nchange_ab) = i
        end if
320     continue
        write(6,322) nchange_ab
322     format('We adopted the following',i3
     &        ,' individual abundances from input:')
C_ursa        write(6,324)((isubst(i),name_gem(isubst(i))
C_ursa     &                      ,abund(isubst(i))),i=1,nchange_ab)
324     format(8(i4,'~',a2,': ',f6.3))
        if(natms-nchange_ab-3.gt.0) then
            if(zscale.eq.1) write(6,323) natms-nchange_ab-3
            if(zscale.ne.1) write(6,325) natms-nchange_ab-3,zscale
        end if
323     format('The other',i3,' abundances were solar.')
325     format('The other',i3
     &    ,' abundances were scaled solar, with zscale =',f8.4)

309     continue

        write(6,132) abund(2),abund(3),(abund(i),i=7,9),abund(26)
132     format('Adopted abundances of H,He,C,N,O,Fe:',8f7.2)

C   the 17 marcs_abundances: H  HE C  N  O  NE NA MG AL SI S  K  CA CR FE NI Ti
        do 160 i=1,17
        abund_marcs(i) = abund(mx_elm(i))
        abtsuji(i) = abund_marcs(i)
160     continue


        aha=abund(2)
        sumabund = 0.d0
        do 4105 i=2,natms
        abund(i)=10.**(abund(i)-aha)
        sumabund = sumabund + abund(i)
4105    continue
        do 4205 k=1,ndp
        do 4205 i=2,natms
        totabk(k,i) = abund(i)/sumabund
4205    continue
        close(unit=2)
C        write(6,130) totabk(1,2), totabk(1,3), totabk(1,7), totabk(1,8),
C     &      totabk(1,9), totabk(1,26)
C130     format('Rel.abund (atoms-i/atoms-all) totabk(1,i) from gem_init'
C     &    ,' for i = H,He,C,N,O,Fe:',/1p6e9.2)
C        write(6,231) (name_gem(i),i=1,natms)
C        write(6,232) (abund(i),i=1,natms)
C231     format(10a8)
C232     format(10f8.5)



      RETURN
      END

C
C                                                           I
C
C
      subroutine tstgem (tcal,pcal,ntau)

!-----------------------------------------------------------------
!
!     This program is designed to test the GEM subroutines.
!
!     27/9-2000 JFF
!
!-----------------------------------------------------------------

      implicit real*16 (a-h,o-z)
      include 'parameter.inc'

      integer i,j,un1,un2,ie,k,keq,ipos,kdp,ntau
      real*16 tiny,small,t,p
      parameter(un1=46,un2=47)            ! Unit numbers for I/O
      parameter(tiny=1.0d-40,small=1.0d-12)
      real*16 totab(natms)                 ! Bulk elemental abundances
      real*16 totabk                       ! ditto, but possibly depth variable
      common /ctotabk/totabk(ndp,natms)
      real*16 abin(nspec),about(nspec)     ! Input and output solutions
      real*16 gibbs(nspec,5)               ! Gibbs' coefficient array
      real*16 compmtx(natms,nspec)         ! Composition matrix
      logical first,firstit,mol_eq_scr            ! Is this the initial run?
      character*8 name(nspec)
      real*16 junk1,junk2,junk3,junk4,junk5! Dummies
      real*16 initsum,finalsum             ! Total fictitious and partial
                                          ! pressures
      common /cabink/abink(ndp,nspec)
      dimension tcal(ndp),pcal(ndp)
      common /cgeminp/tcall(ndp),pcall(ndp)
      common /cmarcsgem/ kmp(0:99),kgem(99),niden
      common /cgem/pres_gem(ndp,nspec) 
      common /cgemnames/natms_gem,nions_gem,nspec_gem,name_gem(nspec)
C atms,ions,spec ~ highest index of neutral atoms, ions, species total
      character name_gem*8
      common/gemcom/gibbs,compmtx
      common/cmxread/pel_ion(ndp,16,4), elabund(ndp,16), 
     &                pprel(ndp,nspec)        !elabund=#atoms_i/atom_in_gas
      data firstit, mol_eq_scr /.true., .true./
C     mol_eq_scr = .true.  => compute mol.equilibrium from scratch
      common /cnrit/ nrit(ndp)

      do 445 k=1,ntau
      tcall(k) = tcal(k)
      pcall(k) = pcal(k)
445   continue


      first=.true.
      initsum=0.0d0
      do 10 i=1,natms
         initsum=initsum+totab(i)              ! Initial abundance sum
 10   continue
      open(unit=un1,file='tstgem.ud',status='unknown')
      open(unit=un2,file='data/gfits.data',
     &                                    status='old',readonly)
         do 121 i=1,nspec
            read(un2,166) name(i)
            read(un2,*) junk1,junk2,junk3,junk4,junk5   ! Fit data
 121      continue
166      format(a8)
      do 13 j=1,ntau
         totab(1) = tiny
         do 131 i=2,natms
131      totab(i) = totabk(j,i)
         finalsum=0.0d0
C it practically doesn't matter whether the following lines are included 
C or not (with abink(j-1,i)), since abin(i) is almost the same as abink(j-1,i) 
C (while abink(j,i) at this place is often zero).
C???         if(firstit .eq. .false.) then
C???         do 3111 i=1,nspec
C???3111     abin(i) = abink(j,i)
C???         end if

C         if(j.eq.1) write(6,*)'abin,abink(j,i),abink(j-1,i) in tstgem:'
C         if(j.eq.2 .or.j.eq.ntau) 
C     &                                      write(6,*)'j=',j
C      do 3111 ik=1,nspec
C         if(j.eq.2 .or.j.eq.ntau) then
C         if(abin(ik).ge.1.e-5 .or. abink(j,ik).ge.1.e-5
C     &   .or. abink(j-1,ik).ge.1.e-5
C     &   .or. ik.le.3 .or.ik.eq.128 .or.ik.eq.133 .or.ik.eq.152)
C     & write(6,3113) ik,name_gem(ik),abin(ik),abink(j,ik),abink(j-1,ik)
C         end if
C3111    continue
C3113    format(i3,2x,a8,2x,1p3e12.3)




         t=tcall(j)                          ! Temperature scale
         p=pcall(j)                          ! Pressure    scale
         call gem(t,p,totab,abin,about,first,j)  ! Actual call to GEM
         do 11 i=1,nspec
            finalsum=finalsum+about(i)         ! Final abundance sum 
 11      continue
         do 12 i=1,nspec
            abin(i)=about(i)                            ! Recycle data
            junk1=about(i)/finalsum                     ! Normalize
            abink(j,i) = about(i)/finalsum 
 12      continue
 13   continue
      close(un1)
      close(un2)
      firstit=.false.

      write(6,*)'Number of iterations in gem for depth 1-ntau:'
      write(6,6112) (nrit(k),k=1,ntau)
6112  format(20i4)

      return
      end



!-----------------------------------------------------------------
!
! Inside MARCS main program do the following: Read Gibbs' energy 
! fits (GIBBSREAD), composition matrix (GEMCMSETUP) and abundances 
! (GEMABSET) into right variables.
! After completion/convergence convert abundances back into MARCS
! format. Set NSPEC and NATMS to their respective values {(246,34) 
! for the common subsample and (49,892) for the full JANAF set}.
!              JFF 18/4-2000
!              JFF 19/9-2000 correction of NSPEC and NATMS values.
! Note, NSPEC is including NATMS, the number of non-elemental 
! species is therefore: NSPEC - NATMS = 212 or 843.
!              JFF 21/9-2000 
!
! Note, the common subsample uses another file-set than the full
! sample. These are not yet converted to GEM-readable format.
!              JFF 27/9-2000
!
! The common sample are now converted to GEM-format, in the files
! TSUCOMP and TSUDATA.
!              JFF 6/11-2000
!
! NSPEC=255, NATMS=48 with TSUCOMP and TSUDATA
!
!-----------------------------------------------------------------

      subroutine gem(temp,pres,totabund,abund1,abund2,first,ktau)

!-----------------------------------------------------------------
!
!     This is the control subroutine, the one called from the
!     main program.
!
!     The algorithm is taken from:
!      White, Johnson & Dantzig; 
!      "Chemical Equilibrium in Complex Mixtures";
!      Journal of Chemical Physics, vol. 28, #5, May 1958, p. 751
!
!     Henceforth referred to as "WJD"
!
!   IMPORTANT!
!     The subroutines 'GEMCMSETUP' and 'GIBBSREAD' should be moved
!     outside this subroutine 'GEM' and be called once from the
!     main program if more than one equilibrium composition is 
!     needed.
!
!     22/09-2000: Jens Falkesgaard
!     10/01-2000: Charge conservation-part of MASSCONS updated
!
!-----------------------------------------------------------------

      implicit real*16 (a-h,o-z)
      include 'parameter.inc'

!--external variables

      real*16 temp,pres                              ! T and P in WJD
      real*16 totabund(natms)                        ! b in WJD
      real*16 abund1(nspec),abund2(nspec)            ! y and x in WJD
      real*16 gibbs(nspec,5)                         ! Gibbs' coeffs.
      real*16 compmtx(natms,nspec)                   ! a in WJD
      logical first

!--internal variables

      integer i,j,count,conprbl,un                  ! iteration params.
      parameter(un=40)                              ! unit number
      logical sofar,newstep                         ! to step or not...
      real*16 absum1,absum2                          ! y-bar and x-bar
      real*16 gef(nspec)                             ! c in WJD
      real*16 pge(nspec)                             ! partial Gibbs' energies
      real*16 eq18a(natms+1,natms+1)                 ! r in WJD (A-matrix)
      real*16 eq18b(natms+1)                         ! B-vector in WJD
      real*16 indx(natms+1)                          ! vector for solvelin
      real*16 delta(nspec)                           ! relative changes
      real*16 oldnorms(6)                            ! delta vector norms
      real*16 norm,avgnorm,lambda

      common/gemcom/gibbs,compmtx
      character name_gem*8
      common /cgemnames/natms_gem,nions_gem,nspec_gem,name_gem(nspec)
      common /cnrit/ nrit(ndp)
      common /cabink/abink(ndp,nspec)

!      print*,'In GEM'
C     call gemcmsetup                            ! Make comp. matrix
C     call gibbsread                             ! Read polynomials

       if(first) call gemabset(temp,absum1,totabund,abund1,first,ktau) 
CUGJ      call gemabset(temp,absum1,totabund,abund1,first,ktau) ! Sets up abundances


      count=1

      do 10 i=1,6
         oldnorms(i)=1.0d0                       ! Initialize norms
 10   continue

 11   continue


      if(count.eq.3000.and.ktau.eq.1) then
           write(6,*) ' count,ktau=',count,ktau
           sum = 0.0
        do 4112 i=1,nspec
           sum = sum + abund1(i)
           if(abund1(i).ge.1.e-5 .or. abund2(i).ge.1.e-5
     &       .or. i.le.3 .or.i.eq.128 .or.i.eq.133 .or.i.eq.152)
     &       write(6,4113) i,name_gem(i),abund1(i),abund2(i)
4113    format(i3,2x,a8,2x,1p2e12.3)
4112    continue
           print*,'sum of abundances (abund1):',sum
      else if(count.ge.3001) then                     ! Prevent infinite loop
         print *,'Convergence failure!'
           write(6,*) ' count,ktau=',count,ktau
           sum = 0.0
        do 12 i=1,nspec
           sum = sum + abund2(i)
           if(abund1(i).ge.1.e-5 .or. abund2(i).ge.1.e-5 
     &       .or. i.le.3 .or.i.eq.128 .or.i.eq.133 .or.i.eq.152)
     &       write(6,4113) i,name_gem(i),abund1(i),abund2(i)
C          print*, i,abund2(i)                  ! Prints last composition
12      continue
           print*,'sum of abundances (abund2):',sum
         stop
      endif


      call gemgibbs(temp,pres,gef,abund1,absum1,pge) ! Computes partial Gibbs' 
                                                     ! energies
      call gemmatrx(abund1,totabund,pge,eq18a,eq18b) ! Set up system matrix

      j=natms+1                                  ! natms is common; j avoids
                                                 ! case-specific code in the
                                                 ! general subrt. SOLVELIN.

      call solvelin(eq18a,j,j,indx,eq18b)            ! Solve the system
      call gemstep(pge,gef,eq18b,absum1,absum2,abund1,abund2) ! Use solution
      call masscons(abund2,totabund)             ! Enforce mass conservation
      sofar=.true.
      newstep=.true.
      conprbl=0                             ! # species that didn't converge
      norm=0.0d0
      avgnorm=0.0d0
      do 13 i=1,6
         avgnorm=avgnorm+oldnorms(i)        ! average of old delta norms
 13   continue
      avgnorm=avgnorm/6.0d0
      do 14 i=1,nspec
         delta(i)=abund2(i)-abund1(i)
         norm=norm+(delta(i)*delta(i))      ! norm-square of delta
         delta(i)=abs(delta(i)/abund1(i))
         if((delta(i).gt.1.5d-2)) then      ! If rel. change too big (*)
            sofar=.false.
            conprbl=conprbl+1               ! This species didn't converge
!            write(un,*) 'Convergence problems for species: ',i,delta(i)
         endif
 14   continue
C     print*,'Convergence problems for ',conprbl,' species.'
      norm=sqrt(norm)
!      print*, '|delta|= ',norm,', avgnorm= ',avgnorm

      if(norm.gt.avgnorm.and.norm.ne.0.0d0) then ! |delta| was too large
         lambda=avgnorm/norm                      
         do 15 i=1,nspec
            delta(i)=(abund2(i)-abund1(i))*8.5d-2*lambda ! reduce delta
            abund2(i)=abund1(i)+delta(i)
 15      continue
      endif
      do 16 i=1,5
      oldnorms(i)=oldnorms(i+1)                  ! Update 'old' norms
 16   continue
      oldnorms(6)=norm
      if(.not.sofar) then 
         newstep=.true.                          ! (*) do another loop.
       else 
         newstep=.false.                         ! If acceptable stop (+)
      endif
      if(newstep) then
         absum1=0.0d0
         do 17 i=1,nspec
            abund1(i)=abund2(i)                  ! Move values
            absum1=absum1+abund2(i)
 17      continue
         count=count+1                           ! Update counter
         goto 11
      endif                                      ! (+) ...here.
!      close(un)
C     print*,'Iteration steps completed: ',count
!      print*,'Mass checksum  = ',msum
!      print*,'Charge checksum= ',qsum

      nrit(ktau) = count

      return
      end

************

      subroutine gemcmsetup

!-----------------------------------------------------------------
!
!     This subroutine reads and sets up the composition matrix
!     for GEM. 
!
!     To be used inside main program, MARCS or elsewhere, before
!     call to GEMMATRX. There is only need for one run of this
!     subroutine, provided the bulk composition does not change.
!
!     21/9-2000 : Jens Falkesgaard
!     12/10-2000: Tested and found OK, JFF
!
!-----------------------------------------------------------------

C     implicit none
      implicit real*16 (a-h,o-z)
      include 'parameter.inc'

!--external variables
     
      real*16 gibbs(nspec,5)                         ! not used here
      real*16 compmtx(natms,nspec)                   ! a in WJD

!--internal variables

      integer i,j,un,tmpint
      parameter(un=41)                              ! unit number
      real*16 tmpreal                                ! dummy variable
      character*16 name
      integer listno,nsp                            ! # and # of elements

      common/gemcom/gibbs,compmtx

C     open(unit=un,file='/home/ast/falkesgd/src/comp',
      open(unit=un,file='data/comp',
     &                                   status='old',readonly)
      do 12 i=1,nspec
         read(un,166) name                            ! name, discard
         read(un,*) listno                          ! number, discard
         if((listno+1).ne.(i)) print *,'Consistency error!'
         read(un,*) nsp                             ! read # species
         do 11 j=1,nsp                              ! for each
            read(un,*) tmpint, tmpreal              ! read data
            compmtx(tmpint+1,i)=tmpreal             ! & assign
 11         continue
 12      continue
166   format(a16)
      close(un)
      return
      end

***********************

      subroutine gibbsread

!-----------------------------------------------------------------
!
!     This subroutine reads Gibbs' energy fit coefficients from
!     the file 'gfits.data'.
!
!     21/9-2000 : Jens Falkesgaard
!     12/10-2000: Tested OK, JFF
!
!-----------------------------------------------------------------

!--external variables

      implicit real*16 (a-h,o-z)
      include 'parameter.inc'
      real*16 gibbs(nspec,5)                   ! Gibbs' coefficients
      real*16 compmtx(natms,nspec)             ! composition matrix

!--internal variables

      integer i,un
      parameter(un=42)                        ! unit number
      character*8 name
      real*16 coef0,coef1,coef2,coef3,coef4    ! Gibbs' coefficients, tmp

      common/gemcom/gibbs,compmtx

      open(unit=un,file='data/gfits.data',
     &                                   status='old',readonly)
      do 11 i=1,nspec
         read(un,166) name
         read(un,*) coef0,coef1,coef2,coef3,coef4
         gibbs(i,1)=coef0
         gibbs(i,2)=coef1                     ! Read and assign for each
         gibbs(i,3)=coef2                     ! species.
         gibbs(i,4)=coef3
         gibbs(i,5)=coef4
 11   continue
166   format(a8)
      close(un)
      return
      end

***********************

      subroutine gemabset(temp,absum1,totabund,abund1,first,ktau)

!-----------------------------------------------------------------
!
!     This subroutine sets the initial abundances.
!     It is only called during the first run.
!     Additionally it computes y-bar from WJD.
!
!     25/9-2000 : Jens Falkesgaard,
!                  explicit values changed to double precision
!
!     12/10-2000: GEMABSET tests nOK! abund1(CO) is not set!!
!                 Fixed, JFF
!     13/10-2000: Re-fixed CO error, JFF
!     30/10-2000: Non-atomic part removed.
!
!-----------------------------------------------------------------

      implicit real*16 (a-h,o-z)
      include 'parameter.inc'

!--external variables

      real*16 temp                                   ! T in WJD
      real*16 absum1                                 ! y-bar in WJD
      real*16 gibbs(nspec,5)                         ! Gibbs' coeffs.
      real*16 totabund(natms)                        ! b in WJD
      real*16 abund1(nspec)                          ! y in WJD
      real*16 compmtx(natms,nspec)                   ! Composition matrix
      logical first                                 ! First time?

!--internal variables

      common/gemcom/gibbs,compmtx
      common /cabink/abink(ndp,nspec)
      character name_gem*8
      common /cgemnames/natms_gem,nions_gem,nspec_gem,name_gem(nspec)


      do 11 i=1,natms
         if(ktau.eq.1) then
         abund1(i)=totabund(i)        ! Set elemental abundances
         else
         abund1(i)=abink(ktau-1,i)
         end if
 11   continue
      do 12 i=natms+1,nspec
         if(ktau.eq.1) then
         abund1(i)=1.0d-40            ! Set molecular and ionic abs.
         else
         abund1(i)=abink(ktau-1,i)
         end if
 12   continue
      absum1=0.0d0
      do 13 i=1,nspec
         absum1=absum1+abund1(i)      ! find y-bar
 13   continue
      first=.false.


      return
      end

***********************

      subroutine gemgibbs(temp,pres,gef,abund1,absum1,pge)

!-----------------------------------------------------------------
!
!     This subroutine computes the partial Gibbs' energies f_i(Y)
!     used in WJD, equations 2 and 15.
!
!     25/9-2000 : Jens Falkesgaard, 
!                  explicit values changed to double precision
!
!-----------------------------------------------------------------

      implicit real*16 (a-h,o-z)
      include 'parameter.inc'

!--external variables

      real*16 absum1
      real*16 temp,pres                              ! T and P in WJD
      real*16 gef(nspec)                             ! c in WJD
      real*16 pge(nspec)                             ! partial Gibbs' energies
      real*16 abund1(nspec)                          ! y in WJD
      real*16 gibbs(nspec,5)                         ! Gibbs' coeffs.
      real*16 compmtx(natms,nspec)                   ! composition matrix

!--internal variables

      integer i,j,un
      parameter(un=43)                              ! unit number
      real*16 tt
      real*16 c(nspec)

      common/gemcom/gibbs,compmtx
      
!      print*,'In GEMGIBBS'
!      open(unit=un,file='gemgibbs.ud',status='unknown')
      tt=temp/1.0d3
      do 11 i=1,nspec
         c(i)=0.0d0                                 ! set to zero
 11   continue   
      do 13 i=1,nspec
         do 12 j=5,1,-1
            c(i)=(c(i)*tt)+gibbs(i,j)               ! using JFF Gibbs' fit
 12      continue
         c(i)=1.0d3*c(i)
         c(i)=c(i)/(8.3145d0*temp)                  ! eq. 2 WJD
         c(i)=c(i)+log(pres/1.0d5)                  !input in SI (N/m^2) -> atm
C        c(i)=c(i)+log(pres/1.0d6)                  !input in cgs -> atmospheres
         gef(i)=c(i)
         pge(i)=0.0d0                               ! Partial Gibbs' energy
 13   continue
   
      do 14 i=1,nspec                               ! eq. 15 WJD
         abund1(i) = max(abund1(i),1.0d-40)
C        if(abund1(i).lt.0.0d0) abund1(i)=1.0d-40
         pge(i)=abund1(i)*(c(i)+log(abs(abund1(i))/absum1))
!         write(6,*) 'pge,abund1(',i,')= ',pge(i),abund1(i)
!         print*, 'pge(',i,')= ',pge(i)
 14   continue
!      close(un)
      return
      end

***********************

      subroutine gemmatrx(abund1,totabund,pge,eq18a,eq18b)

!-----------------------------------------------------------------
!
!     This subroutine sets up the matrix equation AX=B (equations
!     17 and 18 in WJD) for use in SOLVELIN.
!
!     2000: Jens Falkesgaard
!
!-----------------------------------------------------------------

      implicit real*16 (a-h,o-z)
      include 'parameter.inc'

!--external variables

      real*16 abund1(nspec)                          ! y in WJD
      real*16 totabund(natms)                        ! b in WJD
      real*16 pge(nspec)                             ! partial Gibbs' energies
      real*16 eq18a(natms+1,natms+1)                 ! r in WJD (A-matrix)
      real*16 eq18b(natms+1)                         ! B-vector in WJD
      real*16 gibbs(nspec,5)                         ! Gibbs' coefficients
      real*16 compmtx(natms,nspec)                   ! a in WJD

!--internal variables

      common/gemcom/gibbs,compmtx

      do 12 j=1,natms+1                             ! (natms+1)*(natms+1)
         do 11 k=1,natms+1                          ! matrix...
            eq18a(j,k)=0.0d0                        ! initialize to zero
 11      continue
 12   continue   

      do 13 i=1,natms+1                             ! (natms+1) vector...
         eq18b(i)=0.0d0                             ! initialize to zero
 13   continue

      do 16 j=1,natms
         do 15 k=j,natms
            do 14 i=1,nspec
               eq18a(j,k)=eq18a(j,k)+compmtx(j,i)*compmtx(k,i)*abund1(i)
 14         continue                      ! implementation of eq. 17 WJD
            eq18a(k,j)=eq18a(j,k)         ! upper/lower symmetry
 15      continue
 16   continue

      do 17 i=1,natms
         eq18a(i,natms+1)=totabund(i)     ! setting right and lower edge,
         eq18a(natms+1,i)=eq18a(i,natms+1)! lower right element is zero
 17   continue

      do 19 j=1,natms
         do 18 i=1,nspec                  ! setting B-vector elements
            eq18b(j)=eq18b(j)+compmtx(j,i)*pge(i)
 18      continue
 19   continue
      do 20 i=1,nspec                     ! total Gibbs' energy
         eq18b(natms+1)=eq18b(natms+1)+pge(i)
 20   continue
!      print*,'G= ',eq18b(natms+1)
      return
      end

***********************

      subroutine solvelin(a,n,np,indx,b)

!-----------------------------------------------------------------
!
!     This subroutine solves the matrix equation AX=B for X and
!     returns the solution in place of the original B column.
!
!     The routine is based upon 'LUDCMP' and LUBKSB' from 
!     Numerical Recipies in Fortran, 2nd. Ed. 
!     Since this equation is only to be solved once for each 
!     iteration these two routines have been fused into one.
!
!     ref. Numerical Recipies pp. 36-40
!
!     2000: Jens Falkesgaard
!
!-----------------------------------------------------------------

      implicit real*16 (a-h,o-z)
      include 'parameter.inc'

!--internal variables

      integer n,np,nmax
      parameter (nmax=100)
      integer i,ii,imax,j,k,ll
      real*16 d,tiny,aamax,dum,sum
      real*16 a(np,np),vv(nmax),b(n)
      real*16 indx(n)

      tiny=1.0d-40
      d=1.0d0             ! LU-decomposition routine from Numerical Recipies
      do 12 i=1,n         ! modified 7/2-2000, JFF; the two routine fused
         aamax=0.0d0      ! 25/9-2000, JFF; tiny changed to double precision
         do 11 j=1,n
            if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
 11      continue   
         if (aamax.eq.0.0d0) pause 'Singular matrix!'
         vv(i)=1./aamax
 12   continue
      do 19 j=1,n     
         do 14 i=1,j-1
            sum=a(i,j)
            do 13 k=1,i-1
               sum=sum-a(i,k)*a(k,j)
 13         continue
            a(i,j)=sum
 14   continue
      aamax=0.0d0
      do 16 i=j,n
         sum=a(i,j)
         do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
 15      continue
         a(i,j)=sum
         dum=vv(i)*abs(sum)
         if (dum.ge.aamax) then
            imax=i
            aamax=dum
         endif
 16   continue
      if(j.ne.imax) then
         do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
 17      continue
         d=-d
         vv(imax)=vv(j)
      endif
      indx(j)=imax
      if (a(j,j).eq.0.0d0) a(j,j)=tiny
      if (j.ne.n) then
         dum=1./a(j,j)
         do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
 18      continue
      endif
 19   continue

      ii=0               ! Original beginning of 'LUBKSB'
      do 21 i=1,n
         ll=indx(i)
         sum=b(ll)
         b(ll)=b(i)
         if (ii.ne.0) then
            do 20 j=ii,i-1
               sum=sum-a(i,j)*b(j)
 20         continue
         else if (sum.ne.0.0d0) then
            ii=i
         endif
         b(i)=sum
 21   continue
      do 23 i=n,1,-1
         sum=b(i)
         do 22 j=i+1,n
            sum=sum-a(i,j)*b(j)
 22      continue
         b(i)=sum/a(i,i)
 23   continue
      return
      end

***********************

      subroutine gemstep(pge,gef,eq18b,absum1,absum2,abund1,abund2)

!-----------------------------------------------------------------
!
!     This subroutine computes the revised set of abundance 
!     values according to eqs. 14 and 20 in WJD.
!
!     Furthermore it sets abundances that are too low to trace
!     level: 1.0d-40.
!
!     2000: Jens Falkesgaard
!
!-----------------------------------------------------------------

      implicit real*16 (a-h,o-z)
      include 'parameter.inc'

!--external variables

      real*16 absum1,absum2                    ! y-bar and x-bar in WJD
      real*16 abund1(nspec),abund2(nspec)      ! y(i) and x(i) in WJD
      real*16 gef(nspec)                       ! actual Gibbs' energies
      real*16 pge(nspec)                       ! partial Gibbs' energies
      real*16 eq18b(natms+1)
      real*16 compmtx(natms,nspec)             ! composition matrix
      real*16 gibbs(nspec,5)                   ! Gibbs' coefficients

!--internal variables

      integer i,j,un
      parameter (un=44)                       ! unit number
      real*16 eq14sum
      real*16 delta(nspec)

      common/gemcom/gibbs,compmtx

!      open(unit=un,file='gemstep.ud',status='unknown')
      absum2=(eq18b(natms+1)+1)*absum1             ! eq. 19 in WJD
      do 12 i=1,nspec
         eq14sum=0.0d0
         do 11 j=1,natms
            eq14sum=eq14sum+eq18b(j)*compmtx(j,i)
 11      continue
         abund2(i)=-pge(i)+(abund1(i)/absum1)*absum2+(eq14sum*abund1(i))
CUGJ: I set abund2>1.e-40 because abund1 is already so here (???)
CUGJ     if(abund2(i).lt.1.0d-40) abund2(i)=1.0d-40 
         delta(i)=abund2(i)-abund1(i)              ! how large change?
!         write(un,*) 'delta(',i,')= ',delta(i)
 12   continue                                     ! eq. 14 in WJD
      do 13 i=1,nspec
         if(abund2(i).lt.1.0d-40) abund2(i)=1.0d-40 
                                                   ! Underflow protection
 13   continue
!      close(un)
      return
      end

***********************

      subroutine masscons(abund2,totabund)

!-----------------------------------------------------------------
!
!     This subroutine enforces mass and charge conservation 
!     onto the solution suggested by GEMSTEP.
!
!     16/10-2000: Jens Falkesgaard
!     26/10-2000: Charge conservation incorporated
!     10/01-2001: Updated charge conservation
!
!-----------------------------------------------------------------

      implicit real*16 (a-h,o-z)
      include 'parameter.inc'

!--external variables

      real*16 abund2(nspec)                    ! solution suggested by GEMSTEP
      real*16 totabund(natms)                  ! bulk elemental abundances
      real*16 compmtx(natms,nspec)             ! composition matrix
      real*16 gibbs(nspec,5)                   ! Gibbs' coefficients

!--internal variables

      integer i,j,un
      parameter(un=45)                        ! unit number
      real*16 sum(natms)
      real*16 mcdelta(natms)                   ! total mole change
      real*16 qsum,negions,posions,ratio       ! charge checksums

      common/gemcom/gibbs,compmtx

C      print*, 'In MASSCONS'

!      open(unit=un,file='masscons.ud',status='unknown')
      do 11 i=1,natms                            ! check for each element
         sum(i)=0.0d0
         do 10 j=1,nspec                         ! sum up over all species
            sum(i)=sum(i)+abund2(j)*compmtx(i,j) ! remember stoichiometry
 10      continue
C         print*, 'i,totabund,sum:',i,sum(i),totabund(i)
         mcdelta(i)=(sum(i)-totabund(i))/totabund(i)
C         print*, 'mcdelta(',i,')= ',mcdelta(i)
 11   continue

      do 13 i=2,natms                            ! don't mass check electrons
         if(abs(mcdelta(i)).gt.5.0d-2) then 
            ratio=sum(i)/totabund(i)
            do 12 j=2,nspec                      ! if unbalanced...
               if(compmtx(i,j).ne.0.0d0) then    ! ...then adjust
                  if(ratio.gt.1.0d0) abund2(j)=abund2(j)/ratio   
                  if(ratio.lt.1.0d0) abund2(j)=abund2(j)/ratio
               endif
               if(abs(abund2(j)).lt.1.0d-40) abund2(j)=1.0d-40
 12         continue                             ! underflow protection
         endif
 13   continue

      qsum=1.0d-50
      negions=1.0d-50
      posions=1.0d-50

      do 14 i=1,nspec                     ! check for all species
         if(compmtx(1,i).eq.1.0d0) then
            negions=negions+abund2(i)
         endif
         if(compmtx(1,i).eq.-1.0d0) then
            posions=posions+abund2(i)
         endif
         qsum=qsum+abund2(i)*compmtx(1,i) ! allways electron contribution
 14   continue

      ratio=posions/negions

      if(abs(qsum).ge.1.0d-30) then            ! if nOK then adjust
         if(ratio.lt.1.0d0) then                ! if positive imbalance:
            do 15 i=1,nspec
               if(compmtx(1,i).eq.1.0d0) abund2(i)=abund2(i)*ratio
 15         continue                            ! adjust
         endif
         if(ratio.gt.1.0d0) then                ! in negative imbalance:
            do 16 i=1,nspec
               if(compmtx(1,i).eq.-1.0d0) abund2(i)=abund2(i)/ratio
 16         continue                            ! ajust
         endif
      endif
!      close(un)
!      print*, 'Exiting MASSCONS'
      return
      end

!-----------------------------------------------------------------
!
!     End of package
!
!-----------------------------------------------------------------

      SUBROUTINE atomw
C
C Mean atomic weights are stored into the vector watom, AB2001
C Species 1:92 are included, as in irwin.dat
C Deuterium is added as number 93, because it is treated as a seperate
C atom in gfits.data of molecules. Actually W(H)=1.007825, while Earth
C mixture of H and D has W(H+D)=1.0079 as given for watom(1).
C

      implicit real*16 (a-h,o-z)
      character atomname*2
      common /cwatom/watom(200)
      common /catomname/atomname(200)


      DO i=1,200
         watom(i) = 0.0
      ENDDO

      watom(1)  =   1.0079 !H  Hydrogen
      watom(2)  =   4.0026 !He Helium
      watom(3)  =   6.941  !Li Lithium
      watom(4)  =   9.0121 !Be Beryllium
      watom(5)  =  10.81   !B  Boron
      watom(6)  =  12.011  !C  Carbon
      watom(7)  =  14.0067 !N  Nitrogen
      watom(8)  =  15.9994 !O  Oxygen
      watom(9)  =  18.9984 !F  Fluorine
      watom(10) =  20.179  !Ne Neon
      watom(11) =  22.9897 !Na Sodium
      watom(12) =  24.305  !Mg Magnesium
      watom(13) =  26.9814 !Al Aluminum
      watom(14) =  28.0855 !Si Silicon
      watom(15) =  30.9737 !P  Phosphorus
      watom(16) =  32.06   !S  Sulfur
      watom(17) =  35.453  !Cl Chlorine
      watom(18) =  39.948  !Ar Argon
      watom(19) =  39.0983 !K  Potassium
      watom(20) =  40.08   !Ca Calcium
      watom(21) =  44.9559 !Sc Scandium
      watom(22) =  47.88   !Ti Titanium
      watom(23) =  50.9415 !V  Vanadium
      watom(24) =  51.996  !Cr Chromium
      watom(25) =  54.9380 !Mn Manganese
      watom(26) =  55.847  !Fe Iron
      watom(27) =  58.9332 !Co Cobalt
      watom(28) =  58.96   !Ni Nickel
      watom(29) =  63.546  !Cu Copper
      watom(30) =  65.38   !Zn Zinc
      watom(31) =  69.72   !Ga Gallium
      watom(32) =  72.59   !Ge Germanium
      watom(33) =  74.9216 !As Arsenic
      watom(34) =  78.96   !Se Selenium
      watom(35) =  79.904  !Br Bromine
      watom(36) =  83.80   !Kr Krypton
      watom(37) =  85.4678 !Rb Rubidium
      watom(38) =  87.62   !Sr Strontium
      watom(39) =  88.9059 !Y  Yttrium
      watom(40) =  91.22   !Zr Zirconium
      watom(41) =  92.9064 !Nb Niobium
      watom(42) =  95.94   !Mo Molybdenum
      watom(43) =  97.907  !Tc Technetium
      watom(44) = 101.07   !Ru Ruthenium
      watom(45) = 102.9055 !Rh Rhodium
      watom(46) = 106.42   !Pd Palladium
      watom(47) = 107.868  !Ag Silver
      watom(48) = 112.41   !Cd Cadmium
      watom(49) = 114.82   !In Indium
      watom(50) = 118.69   !Sn Tin
      watom(51) = 121.75   !Sb Antimony
      watom(52) = 127.60   !Te Tellurium
      watom(53) = 126.9045 !I  Iodine
      watom(54) = 131.29   !Xe Xenon
      watom(55) = 132.9054 !Cs Cesium
      watom(56) = 137.33   !Ba Barium
      watom(57) = 138.9055 !La Lanthanum
      watom(58) = 140.12   !Ce Cerium
      watom(59) = 140.9077 !Pr Praseodymium
      watom(60) = 144.24   !Nd Neodymium
      watom(61) = 144.913  !Pm Promethium
      watom(62) = 150.36   !Sm Samarium
      watom(63) = 151.96   !Eu Europium
      watom(64) = 157.25   !Gd Gadolinium
      watom(65) = 158.9254 !Tb Terbium
      watom(66) = 162.50   !Dy Dysprosium
      watom(67) = 164.9304 !Ho Holmium
      watom(68) = 167.26   !Er Erbium
      watom(69) = 168.9342 !Tm Thulium
      watom(70) = 173.04   !Yb Ytterbium
      watom(71) = 174.967  !Lu Lutetium
      watom(72) = 178.49   !Hf Hafnium
      watom(73) = 180.9479 !Ta Tantalum
      watom(74) = 183.85   !W  Tungsten
      watom(75) = 186.207  !Re Rhenium
      watom(76) = 190.2    !Os Osmium
      watom(77) = 192.22   !Ir Iridium
      watom(78) = 195.08   !Pt Platinum
      watom(79) = 196.9665 !Au Gold
      watom(80) = 200.59   !Hg Mercury
      watom(81) = 204.383  !Tl Thallium
      watom(82) = 207.2    !Pb Lead
      watom(83) = 208.9804 !Bi Bismuth
      watom(84) = 208.982  !Po Polonium
      watom(85) = 209.987  !At Astatine
      watom(86) = 222.018  !Rn Radon
      watom(87) = 223.020  !Fr Francium
      watom(88) = 226.0254 !Ra Radium
      watom(89) = 227.0278 !Ac Actinium
      watom(90) = 232.0381 !Th Thorium
      watom(91) = 231.0359 !Pa Protactinium
      watom(92) = 238.051  !U  Uranium
      watom(93) = 2.01410  !D  Deuterium

      RETURN
      END

C---------------------------------------------------------

      SUBROUTINE atomnam
C
C 1 or 2 character names of atoms.
C Species 1:92 are included, as in irwin.dat
C


      implicit REAL*16 (a-h,o-z)
      character atomname*2
      common /cwatom/watom(200)
      common /catomname/atomname(200)

      DO i=1,200
         atomname(i) = '  '
      ENDDO

      atomname(1)  = 'H ' ! Hydrogen
      atomname(2)  = 'He' ! Helium
      atomname(3)  = 'Li' ! Lithium
      atomname(4)  = 'Be' ! Beryllium
      atomname(5)  = 'B ' ! Boron
      atomname(6)  = 'C ' ! Carbon
      atomname(7)  = 'N ' ! Nitrogen
      atomname(8)  = 'O ' ! Oxygen
      atomname(9)  = 'F ' ! Fluorine
      atomname(10) = 'Ne' ! Neon
      atomname(11) = 'Na' ! Sodium
      atomname(12) = 'Mg' ! Magnesium
      atomname(13) = 'Al' ! Aluminum
      atomname(14) = 'Si' ! Silicon
      atomname(15) = 'P ' ! Phosphorus
      atomname(16) = 'S ' ! Sulfur
      atomname(17) = 'Cl' ! Chlorine
      atomname(18) = 'Ar' ! Argon
      atomname(19) = 'K ' ! Potassium
      atomname(20) = 'Ca' ! Calcium
      atomname(21) = 'Sc' ! Scandium
      atomname(22) = 'Ti' ! Titanium
      atomname(23) = 'V ' ! Vanadium
      atomname(24) = 'Cr' ! Chromium
      atomname(25) = 'Mn' ! Manganese
      atomname(26) = 'Fe' ! Iron
      atomname(27) = 'Co' ! Cobalt
      atomname(28) = 'Ni' ! Nickel
      atomname(29) = 'Cu' ! Copper
      atomname(30) = 'Zn' ! Zinc
      atomname(31) = 'Ga' ! Gallium
      atomname(32) = 'Ge' ! Germanium
      atomname(33) = 'As' ! Arsenic
      atomname(34) = 'Se' ! Selenium
      atomname(35) = 'Br' ! Bromine
      atomname(36) = 'Kr' ! Krypton
      atomname(37) = 'Rb' ! Rubidium
      atomname(38) = 'Sr' ! Strontium
      atomname(39) = 'Y ' ! Yttrium
      atomname(40) = 'Zr' ! Zirconium
      atomname(41) = 'Nb' ! Niobium
      atomname(42) = 'Mo' ! Molybdenum
      atomname(43) = 'Tc' ! Technetium
      atomname(44) = 'Ru' ! Ruthenium
      atomname(45) = 'Rh' ! Rhodium
      atomname(46) = 'Pd' ! Palladium
      atomname(47) = 'Ag' ! Silver
      atomname(48) = 'Cd' ! Cadmium
      atomname(49) = 'In' ! Indium
      atomname(50) = 'Sn' ! Tin
      atomname(51) = 'Sb' ! Antimony
      atomname(52) = 'Te' ! Tellurium
      atomname(53) = 'I ' ! Iodine
      atomname(54) = 'Xe' ! Xenon
      atomname(55) = 'Cs' ! Cesium
      atomname(56) = 'Ba' ! Barium
      atomname(57) = 'La' ! Lanthanum
      atomname(58) = 'Ce' ! Cerium
      atomname(59) = 'Pr' ! Praseodymium
      atomname(60) = 'Nd' ! Neodymium
      atomname(61) = 'Pm' ! Promethium
      atomname(62) = 'Sm' ! Samarium
      atomname(63) = 'Eu' ! Europium
      atomname(64) = 'Gd' ! Gadolinium
      atomname(65) = 'Tb' ! Terbium
      atomname(66) = 'Dy' ! Dysprosium
      atomname(67) = 'Ho' ! Holmium
      atomname(68) = 'Er' ! Erbium
      atomname(69) = 'Tm' ! Thulium
      atomname(70) = 'Yb' ! Ytterbium
      atomname(71) = 'Lu' ! Lutetium
      atomname(72) = 'Hf' ! Hafnium
      atomname(73) = 'Ta' ! Tantalum
      atomname(74) = 'W ' ! Tungsten
      atomname(75) = 'Re' ! Rhenium
      atomname(76) = 'Os' ! Osmium
      atomname(77) = 'Ir' ! Iridium
      atomname(78) = 'Pt' ! Platinum
      atomname(79) = 'Au' ! Gold
      atomname(80) = 'Hg' ! Mercury
      atomname(81) = 'Tl' ! Thallium
      atomname(82) = 'Pb' ! Lead
      atomname(83) = 'Bi' ! Bismuth
      atomname(84) = 'Po' ! Polonium
      atomname(85) = 'At' ! Astatine
      atomname(86) = 'Rn' ! Radon
      atomname(87) = 'Fr' ! Francium
      atomname(88) = 'Ra' ! Radium
      atomname(89) = 'Ac' ! Actinium
      atomname(90) = 'Th' ! Thorium
      atomname(91) = 'Pa' ! Protactinium
      atomname(92) = 'U ' ! Uranium
      atomname(93) = 'D ' ! Deuterium

      RETURN
      END

!-----------------------------------------------------------------------
! Reads in the dust data from the DRIFT output file and makes a table
! of dust opacities
! Juncher 2015
!-----------------------------------------------------------------------
      subroutine drift2marcs

      implicit real*16 (a-h,o-z)
      include 'parameter.inc'
      real*16    :: L0(max_lay), L1(max_lay), object
      real*16, parameter :: pi = 3.14159265
      complex*16 :: M_inc(max_inc,nwl),M_eff0,M_eff(nwl,max_lay)
      logical   :: first
      dimension :: elnr(max_eps), rho(max_lay),V_inc(max_inc,max_lay),
     *             a(max_lay)
      dimension :: wn(2000),p(2),step(2),var(2)
      common/cos/wnos(nwl),conos(ndp,nwl),wlos(nwl),wlstep(nwl),
     *    kos_step,nwtot
      common /cdustdata/ dabstable(max_lay,nwl),dscatable(max_lay,nwl),
     *    eps(max_eps,max_lay),temp(max_lay),pgas(max_lay),n_lay,idust

! Read in dust data from DRIFT
      open(110,file='drift2marcs.dat',status='old',readonly)
      do i=1,24
        read(110,*)
      end do
      read(110,'(i4)') n_inc
      if(n_inc .gt. max_inc) then
        print *, 'Error: increase max_inc.'
        stop
      end if
      do i=1,n_inc
        read(110,*)
      end do
      read(110,'(i4)') n_eps
      if(n_eps .gt. max_eps) then
        print *, 'Error: increase max_eps.'
        stop
      end if
      do i=1,n_eps
        read(110,*) elnr(i)
      end do
      read(110,*)
      read(110,*)
      i = 1
      do
        read(110,'(20x,5e20.12,40x,99e20.12)',iostat=io) 
     *    temp(i),rho(i),pgas(i),L0(i),L1(i),V_inc(1:n_inc,i),
     *    eps(1:n_eps,i)
        if(io .lt. 0) exit
        if(L0(i) .le. 0.) then    
          a(i) = 0.
        else
          a(i) = (3./4./pi)**(1./3.)*L1(i)/L0(i) 
        end if
        i = i + 1
      end do
      close(110)
      n_lay = i-1
      if(n_lay .gt. max_lay) then
        print *, 'Warning: only used the first 1000 layers from DRIFT.'
      end if
      
! Read in optical constants
      call optical_data(n_inc,M_inc)

! Make opacity table
      do i=1,nwtot
        do j=1,n_lay
          if(i .eq. 1) then
            M_eff0 = (0.0,0.0)

            do k=1,n_inc
              M_eff0 = M_eff0 + V_inc(k,j)*M_inc(k,i)
            end do
          else
            M_eff0 = M_eff(i-1,j)
          end if

          call NR(M_inc(:,i),V_inc(:,j),M_eff0,M_eff(i,j),n_inc)

          x = 2.*3.141593*a(j)*1e8/wlos(i)      ! a: cm -> 
          if(x .eq. 0.) then
            q_abs = 0.
            q_sca = 0.
          else
            call mie(M_eff(i,j),x,q_abs,q_sca)
          end if
          dabstable(j,i) = L0(j)*pi*a(j)**2*q_abs
          dscatable(j,i) = L0(j)*pi*a(j)**2*q_sca
          
        end do
      end do
            
      return
      end
      
!-----------------------------------------------------------------------
! Reads in optical data for solids
! Juncher 2015
!-----------------------------------------------------------------------
      subroutine optical_data(n_inc,M_inc)
      
      implicit real*16 (a-h,o-z)
      include 'parameter.inc'
      integer :: i, j, k, io, nlines, counter, n_inc
      real*16  :: a, b, nkdata(3,3000)
      real*16  :: ndata(max_inc,nwl), kdata(max_inc,nwl)
      complex*16 :: M_inc(max_inc,nwl)
      character :: filename(12)*50
      common/cos/wnos(nwl),conos(ndp,nwl),wlos(nwl),wlstep(nwl),
     *    kos_step,nwtot
      
      filename(1) = 'data/DRIFT/nk_TiO2.txt'	   ! TiO2
      filename(2) = 'data/DRIFT/nk_Mg2SiO4.txt'	   ! Mg2SiO4
      filename(3) = 'data/DRIFT/nk_SiO2.txt'	   ! SiO2
      filename(4) = 'data/DRIFT/nk_Fe.txt'	   ! Fe
      filename(5) = 'data/DRIFT/nk_Al2O3.txt'	   ! Al2O3
      filename(6) = 'data/DRIFT/nk_MgO.txt'	   ! MgO     (non-conducting)
      filename(7) = 'data/DRIFT/nk_MgSiO3.txt'     ! MgSiO3  (non-conducting)
      
      do i=1,n_inc
! n values
        open(unit=1,file=filename(i))
        read(1,*) nlines
        counter = 1
        do j=1,nlines
          read(1,*) nkdata(1:3,counter)
          if(nkdata(2,counter) .gt. 0.) then
            counter = counter + 1
          end if
        end do
        close(1)
        nlines = counter - 1
        nkdata(1,:) = 10000.*nkdata(1,:)     ! micron > 
        
        do j=1,nwtot
          if(wlos(j) .lt. nkdata(1,1)) then
            ndata(i,j) = nkdata(2,1)
            cycle
          end if
          if(wlos(j) .gt. nkdata(1,nlines)) then
            if(i .le. 5) then
              print *, 'Missing n-data for inclusion #', i
              stop
            else
              ndata(i,j) = nkdata(2,nlines)
            end if
            cycle
          end if
          do k=1,nlines-1
            if(wlos(j).ge.nkdata(1,k).and.wlos(j).le.nkdata(1,k+1))then
              a=(nkdata(2,k+1)-nkdata(2,k))/(nkdata(1,k+1)-nkdata(1,k))
              b=nkdata(2,k) - a*nkdata(1,k)
              ndata(i,j) = a*wlos(j)+b
              exit
            end if
          end do
        end do
        
! k values
        open(unit=1,file=filename(i))
        read(1,*) nlines
        counter = 1
        do j=1,nlines
          read(1,*) nkdata(1:3,counter)
          if(nkdata(3,counter) .gt. 0.) then
            counter = counter + 1
          end if
        end do
        close(1)
        nlines = counter - 1
        nkdata(1,:) = 10000.*nkdata(1,:)     ! micron > 
        do j=1,nwtot
          if(wlos(j) .lt. nkdata(1,1)) then
            kdata(i,j) = nkdata(3,1)
            cycle
          end if
          if(wlos(j) .gt. nkdata(1,nlines)) then
            if(i .le. 5) then
              print *, 'Missing k-data for inclusion #', n_inc
            else
              kdata(i,j) = nkdata(3,nlines)*wlos(i)/nkdata(1,nlines)
            end if
            cycle
          end if
          do k=1,nlines-1
            if(wlos(j).ge.nkdata(1,k).and.wlos(j).le.nkdata(1,k+1))then
              a=(nkdata(3,k+1)-nkdata(3,k))/(nkdata(1,k+1)-nkdata(1,k))
              b=nkdata(3,k) - a*nkdata(1,k)
              kdata(i,j) = a*wlos(j)+b
              exit
            end if
          end do
        end do
      end do
      
      do i=1,nwtot
        M_inc(1:n_inc,i) = dcmplx(ndata(1:n_inc,i),kdata(1:n_inc,i))
      end do
      
      end



      
!-----------------------------------------------------------------------
! Minimizes a multidimensional function using the Newton-Raphson method
! Juncher 2015
!-----------------------------------------------------------------------   
      subroutine NR(M_inc,V_inc,M_eff0,M_eff,n_inc)
      
      implicit real*16 (a-h,o-z)
      include 'parameter.inc'
      integer,parameter :: itmax = 30
      complex*16 :: M_inc(n_inc), M_eff0, M_eff
      real*16 :: V_inc(n_inc),F(2),F1(2),F2(2),J(2,2)
      real*16 :: corr(2),xold(2),xnew(2)
      
      do i=1,itmax
        call bruggeman(M_eff0,M_inc,V_inc,n_inc,F)

        if(F(1)**2+F(2)**2 .lt. 1e-13) exit
        call jacobian(M_eff0,M_inc,V_inc,J,n_inc)
        call gauss(2,J,corr,F)
        corr = -corr
        
        xold(1) = real(M_eff0)
        xold(2) = imag(M_eff0)
        xnew = xold + corr
        
        if(xnew(1).gt.0 .and. xnew(2).gt.0) then
          M_eff0 = cmplx(xnew(1),xnew(2))
        else
          print *, 'xnew is unphysical: ', xnew(1),xnew(2)
          stop
        end if
        
        if(i .eq. itmax) then
          print *, 'Reached max number of iterations.'
          print *, '|F| = ', F(1)**2+F(2)**2
        end if

      end do
      
      M_eff = M_eff0
      
      return
      end


!-----------------------------------------------------------------------
! Finds the Jacobian of F(M_eff0)
! Juncher 2015
!-----------------------------------------------------------------------
      subroutine jacobian(M_eff0,M_inc,V_inc,J,n_inc)

      implicit none
      integer :: n_inc
      real*16 :: dx, dy, V_inc(n_inc), J(2,2), F1(2), F2(2)
      complex*16 :: M_eff0, M_inc(n_inc)

      dx = 1.e-5*real(M_eff0)
      call bruggeman(M_eff0+cmplx(dx,0.0),M_inc,V_inc,n_inc,F1)
      call bruggeman(M_eff0-cmplx(dx,0.0),M_inc,V_inc,n_inc,F2)
      J(1,1) = (F1(1)-F2(1))/(2.0*dx)
      J(2,1) = (F1(2)-F2(2))/(2.0*dx)

      dy = 1.e-5*imag(M_eff0)
      call bruggeman(M_eff0+cmplx(0.,dy),M_inc,V_inc,n_inc,F1)
      call bruggeman(M_eff0-cmplx(0.,dy),M_inc,V_inc,n_inc,F2)
      J(1,2) = (F1(1)-F2(1))/(2.0*dy)
      J(2,2) = (F1(2)-F2(2))/(2.0*dy)

      return
      end
      
!-----------------------------------------------------------------------
! Bruggeman's equation for F
! Juncher 2015
!-----------------------------------------------------------------------
      subroutine bruggeman(M_eff0,M_inc,V_inc,n_inc,F)

      implicit none
      integer :: i, n_inc
      real*16  :: V_inc(n_inc), F(2)
      complex*16 :: M_eff0, M_inc(n_inc), Fcmplx
      
      Fcmplx = cmplx(0.0,0.0)
      do i=1,n_inc
        Fcmplx = Fcmplx + V_inc(i)*(M_inc(i)**2-M_eff0**2)/
     *                    (M_inc(i)**2+2.*M_eff0**2)
      end do

      F(1) = real(Fcmplx)
      F(2) = imag(Fcmplx)

      return
      end

!-----------------------------------------------------------------------
! Solves a matrix equation of the form: A*x = b using Gaussian 
! elimination and back-substitution
! Juncher 2015
!-----------------------------------------------------------------------
      subroutine gauss(N,a,x,b)

      implicit none
      integer :: N,i,j,k,kmax,piv(N),istat
      real*16  :: a(N,N),x(N),b(N),c,amax

      do i=1,N-1
        kmax = i
        amax = ABS(a(i,i))
        do k=i+1,N
          if(ABS(a(k,i)) .gt. amax) then
            kmax = k
            amax = ABS(a(k,i))
          endif
        end do
  
        if(kmax.ne.i) then
          do j=1,N
            c         = a(i,j)
            a(i,j)    = a(kmax,j)
            a(kmax,j) = c 
          end do
          c       = b(i)
          b(i)    = b(kmax)
          b(kmax) = c 
        endif

        do k=i+1,N
          c = a(k,i) / a(i,i)
          a(k,i) = 0.e0
          do j=i+1,N
            a(k,j) = a(k,j) - c * a(i,j)
          end do
          b(k) = b(k) - c * b(i)
        end do
      end do

      do i=N,1,-1
        c = 0.e0
        if (i.lt.N) then
          do j=i+1,N
            c = c + a(i,j) * x(j)
          end do
        end if
        x(i) = (b(i) - c) / a(i,i)
      end do

      return
      end 

      
      
!-----------------------------------------------------------------------
! Calculates the calculates efficiency factors for absorption and
! scattering (based on the Bohren-Huffman Mie subroutine)
! Juncher 2015
!-----------------------------------------------------------------------
      subroutine mie(M_eff,x,q_abs,q_sca)
      
      implicit real*16 (a-h,o-z)
      complex*16 refrel,M_eff,d(2000000),xi,xi1,an,bn
      
      refrel = M_eff/1.   ! ref. index of sphere / ref. index of medium
      q_sca = 0.
      q_ext = 0.
      q_sca1 = -1.
      q_ext1 = -1.
      
      npoints = max(x+4.*x**(1./3.)+2.,abs(x*refrel))
      if(npoints .gt. 2000000) then
        print *, 'increase max dimension of d'
        stop
      end if
      d(npoints+15) = (0.,0.)
      do i=npoints+14,1,-1
        d(i) = (i+1)/(x*refrel) - 1./(d(i+1)+(i+1)/(x*refrel))
      end do

      psi0 = cos(x)
      psi1 = sin(x)
      chi0 = -sin(x)
      chi1 = cos(x)
      xi1 = dcmplx(psi1,-chi1)
  
      do i=1,npoints
        psi = (2.0*i-1.)*psi1/x-psi0
        chi = (2.0*i-1.)*chi1/x-chi0
        xi = dcmplx(psi,-chi)
          
        an = ((d(i)/refrel+i/x)*psi-psi1)/((d(i)/refrel+i/x)*xi-xi1)
        bn = ((d(i)*refrel+i/x)*psi-psi1)/((d(i)*refrel+i/x)*xi-xi1)

        q_sca = q_sca + (2.*i+1)*(abs(an)**2+abs(bn)**2)
        q_ext = q_ext + (2.*i+1)*real(an+bn)
            
        if(abs(q_sca-q_sca1).lt.1e-8.and.abs(q_ext-q_ext1).lt.1e-8) then
          exit
        end if
    
        psi0 = psi1
        psi1 = psi
        chi0 = chi1
        chi1 = chi
        xi1  = dcmplx(psi1,-chi1)
    
        q_sca1 = q_sca
        q_ext1 = q_ext
      end do
  
      q_sca = (2./x**2)*q_sca
      q_ext = (2./x**2)*q_ext
      if(imag(refrel) .eq. 0.) then
        q_abs = 0.
      else
        q_abs = q_ext - q_sca
      end if
 
      return
      end

!-----------------------------------------------------------------------
! Calculates the dust opacity for each wavelength at each depth layer 
! by interpolating the dust opacity table.
! Juncher 2015
!-----------------------------------------------------------------------
      subroutine dust_opac

      implicit real*16 (a-h,o-z)
      include 'parameter.inc'
      
      common /statec/ppr(ndp),ppt(ndp),pp(ndp),gg(ndp),zz(ndp),dd(ndp),
     * vv(ndp),ffc(ndp),ppe(ndp),tt(ndp),tauln(ndp),ro(ndp),ntau,iter
      common/ci1/fl2(5),parco(45),parq(180),shxij(5),tparf(4),
     *  xiong(16,5),eev,enamn(ndp),sumh(ndp),xkbol,
     *  summ(ndp),nj(16),iel(16),nel
      common /cdustdata/ dabstable(max_lay,nwl),dscatable(max_lay,nwl),
     *    eps(max_eps,max_lay),temp(max_lay),pgas(max_lay),n_lay,idust
      common /cdustopac/ dust_abs(ndp,nwl), dust_sca(ndp,nwl)
      common/cos/wnos(nwl),conos(ndp,nwl),wlos(nwl),wlstep(nwl),
     *    kos_step,nwtot
      dimension pg(ndp),diff(ndp,max_lay),imin(ndp),itemp(12)
      common /cmolrat/ fold(ndp,8),molold,kl
      
      dust_abs(1:ntau,1:nwtot) = 0.
      dust_sca(1:ntau,1:nwtot) = 0.
      
      do i=1,ntau
        kl = i
        call jon(tt(i),ppe(i),1,pgx,rox,dumx,0)
        pg(i) = pgx
      end do
      
      do i=1,ntau  
        do j=1,n_lay
          diff(i,j) = abs((tt(i)-temp(j))/tt(i)) + 
     *                abs((pg(i)-pgas(j))/pg(i))
        end do
        imin(i) = minloc(diff(i,1:n_lay),1)
        
        if(tt(i) .le. temp(n_lay)) then
          dust_abs(i,1:nwtot) = dabstable(imin(i),1:nwtot)
          dust_sca(i,1:nwtot) = dscatable(imin(i),1:nwtot)
        end if
      end do

! Dust opacity file for SYNTOS
      step = (log10(tt(ntau))-log10(tt(1)))/11
      itemp(1) = 1
      do i=2,11
         xtemp = 10**(log10(tt(1)) + (i-1)*step)
        do j=1,ntau-1
          if(tt(j).le.xtemp .and. tt(j+1).gt.xtemp) then
            itemp(i) = j
            exit
          end if
        end do
      end do
      itemp(12) = ntau
      
      open(unit=8,file='osdust_abs.dat')
      write(8,*) " $INPUTOSMOL"
      write(8,*) " MOLID   = 'dabs',"
      write(8,*) " KTEMP   = ", ntau, ","
      write(8,*) " TMOL    = ", tt(1:ntau), ","
      write(8,*) " NWNOS   = ", nwtot, ","
      write(8,*) " VKMS    = 3.0"
      write(8,*) " $END"
      do i=1,nwtot
       write(8,'(f10.3,1p100e11.3)')wnos(i),(dust_abs(1:ntau,i))
      end do
      close(8)

      open(unit=8,file='osdust_sca.dat')
      write(8,*) " $INPUTOSMOL"
      write(8,*) " MOLID   = 'dsca',"
      write(8,*) " KTEMP   = ", ntau, ","
      write(8,*) " TMOL    = ", tt(1:ntau), ","
      write(8,*) " NWNOS   = ", nwtot, ","
      write(8,*) " VKMS    = 3.0,"
      write(8,*) " L_PER_STELLAR   = 1"
      write(8,*) " $END"
      do i=1,nwtot
       write(8,'(f10.3,1p100e11.3)')wnos(i),(dust_sca(1:ntau,i))
      end do
      close(8)

      return
      end

!-----------------------------------------------------------------------
! Updates the initial element abundances (read from elabund.dat) with
! the dust depleted abundances
! From DRIFT: Mg Si Ti O Fe Al
! 17 MARCS: H HE C N O NE NA MG AL SI S K CA CR FE NI Ti
! Juncher 2015
!-----------------------------------------------------------------------
      subroutine dust_eps
      
      implicit real*16 (a-h,o-z)
      include 'parameter.inc'

      common/statec/ppr(ndp),ppt(ndp),pp(ndp),gg(ndp),zz(ndp),dd(ndp),
     * vv(ndp),ffc(ndp),ppe(ndp),tt(ndp),tauln(ndp),ro(ndp),ntau,iter
      common /cdustdata/ dabstable(max_lay,nwl),dscatable(max_lay,nwl),
     *    eps(max_eps,max_lay),temp(max_lay),pgas(max_lay),n_lay,idust
      common/ci5/abmarcs(17,ndp),anjon(17,5),h(5),part(17,5),dxi,
     *           f1,f2,f3,f4,f5,xkhm,xmh,xmy(ndp)
      common/ci1/fl2(5),parco(45),parq(180),shxij(5),tparf(4),
     *  xiong(16,5),eev,enamn(ndp),sumh(ndp),xkbol,
     *  summ(ndp),nj(16),iel(16),nel
      common/ci9/ai(16)
      common /tsuji/ nattsuji,nmotsuji,parptsuji(500),abtsuji(17,ndp)
      dimension pg(ndp),diff(ndp,max_lay),imin(ndp),sum(ndp)
      common /cmolrat/ fold(ndp,8),molold,kl
      data eev/1.602095e-12/,xmh/1.67339e-24/      

      do i=1,ntau
        do j=1,n_lay-1
          if(tt(i).lt.temp(1)) then
            abmarcs(8,i)  = log10(eps(1,1))+12.       ! 12 Mg
            abmarcs(10,i) = log10(eps(2,1))+12.       ! 14 Si
            abmarcs(17,i) = log10(eps(3,1))+12.       ! 22 Ti
            abmarcs(5,i)  = log10(eps(4,1))+12.       !  8  O
            abmarcs(15,i) = log10(eps(5,1))+12.       ! 26 Fe
            abmarcs(9,i)  = log10(eps(6,1))+12.       ! 13 Al
            exit
          end if
             
          if(tt(i).ge.temp(j) .and. tt(i).le.temp(j+1)) then
            abmarcs(8,i)  = log10(eps(1,j))+12.       ! 12 Mg
            abmarcs(10,i) = log10(eps(2,j))+12.       ! 14 Si
            abmarcs(17,i) = log10(eps(3,j))+12.       ! 22 Ti
            abmarcs(5,i)  = log10(eps(4,j))+12.       !  8  O
            abmarcs(15,i) = log10(eps(5,j))+12.       ! 26 Fe
            abmarcs(9,i)  = log10(eps(6,j))+12.       ! 13 Al
            exit
          end if
        end do
      end do

      abtsuji = abmarcs
      
      write(6,'(a25)') 'Dust depleted abundances:'
      write(6,'(7a6)') 'k', 'Mg','Si','Ti','O','Fe','Al'
      do i=1,ntau,10
        write(6,'(i6,6f6.2)') i, abmarcs(8,i),abmarcs(10,i),
     *    abmarcs(17,i),abmarcs(5,i),abmarcs(15,i),abmarcs(9,i)
      end do
      
      sum(1:ntau) = 0.
      xmy(1:ntau) = 0.
      
      do i=1,nel
        abmarcs(i,1:ntau)=10.**abmarcs(i,1:ntau)
        sum(1:ntau)=SUM(1:ntau)+abmarcs(I,1:ntau)
      end do
      abmarcs(17,1:ntau)=10.**abmarcs(17,1:ntau)
          
      xmy(1:ntau)=0.
      aha=abmarcs(1,1)
      do i=1,nel
        abmarcs(i,1:ntau)=abmarcs(i,1:ntau)/aha
        summ(1:ntau)=summ(1:ntau)+abmarcs(i,1:ntau)
        xmy(1:ntau)=xmy(1:ntau)+abmarcs(i,1:ntau)*ai(i)
      end do
      abmarcs(17,1:ntau)=abmarcs(17,1:ntau)/aha
      xmy(1:ntau)=xmy(1:ntau)/ai(1)
      sumh(1:ntau)=sum(1:ntau)/aha-1.
      summ(1:ntau)=summ(1:ntau)-abmarcs(1,1:ntau)-abmarcs(3,1:ntau)-
     *  abmarcs(4,1:ntau)-abmarcs(5,1:ntau)
     
      enamn(1:ntau) = eev/(xmh*xmy(1:ntau))
            
      return
      end

!-----------------------------------------------------------------------
! Writes a MARCS output file to be read by DRIFT
! Juncher 2015
!-----------------------------------------------------------------------
      subroutine marcs2drift
      
      implicit real*16 (a-h,o-z)
      include 'parameter.inc'
      character flag(ndp)
      dimension pg(ndp),surfgrav(ndp),v(ndp),emu(ndp),rad(ndp)
      dimension flip_rad(ndp), abundances(ndp,100), pg2(ndp)
      common /tauc/tau(ndp),dtauln(ndp),jtau
      common /cg/grav,konsg /cteff/teff,flux
      common /masse/relm
      common /mixc/palfa,pbeta,pny,py /cvfix/vfix
      common /statec/ppr(ndp),ppt(ndp),pp(ndp),gg(ndp),zz(ndp),dd(ndp),
     & vv(ndp),ffc(ndp),ppe(ndp),tt(ndp),tauln(ndp),ro(ndp),ntau,iter
      common /cstyr/mihal,noconv
      common /rossc/xkapr(ndp),cross(ndp)
      common /cabinit/abinit(natms),kelem(natms),nelem
      common /tsuji/ nattsuji,nmotsuji,parptsuji(500),abtsuji(17,ndp)
      common /cdustdata/ dabstable(max_lay,nwl),dscatable(max_lay,nwl),
     *    eps(max_eps,max_lay),temp(max_lay),pgas(max_lay),n_lay,idust
      common /cmolrat/ fold(ndp,8),molold,kl

      sun_rad = 6.96342e10   ! cm
      
! Radius
      relr = sqrt(relm/grav*10**(4.44))

      flip_rad(1:ntau) = 0.
      do k=1,ntau
        kl = k
        furem = fure
        call termo(tt(k),ppe(k),ppr(k),ptot,rro,cp,cv,agrad,q,u2)
        fure = 1./(xkapr(k)*rro)
        if(k .eq. 1) cycle
        flip_rad(k) = flip_rad(k-1) + (tau(k)-tau(k-1))*(fure+furem)*0.5
      end do
      rad(1:ntau) = 0.
      do k=1,ntau
        rad(k) = flip_rad(ntau-k+1)+relr*sun_rad
      end do
      
! Gas pressure
      do k=1,ntau
        kl = k
        call jon(tt(k),ppe(k),1,pgx,rox,dumx,0)
        pg(k) = PGx
      end do

! Surface gravity
      do k=1,ntau
        surfgrav(k) = (relr*sun_rad/rad(k))*grav
      end do 

! Convective velocity
      do k=1,ntau
        if(k .gt. 1) go to 13
        v(k) = 0.
        go to 15
        
13      if(k .eq. ntau) go to 14
        ya=(tau(k)-tau(k-1))/(tau(k+1)-tau(k-1))
        yb=1.-ya
        v(k)=ya*vv(k+1)+yb*vv(k)
        if(vv(k).gt.0..and.vv(k+1).gt.0.) v(k)=
     &    exp(ya*log(vv(k+1))+yb*log(vv(k)))
        go to 15
        
14      continue
        ya=(2.*tau(k)-tau(k-1)-tau(k-2))/(tau(k)-tau(k-2))
        yb=1.-ya
        v(k)=ya*vv(k)+yb*vv(k-1)
15      continue
      end do

! Convection flag
      do k=1,ntau
        if(v(k) .eq. 0.) then
          flag(k) = 'F'
        else
          flag(k) = 'T'
        end if
      end do   
      
! Mean molecular mass
      do k=1,ntau
        kl = k
        call termo(tt(k),ppe(k),ppr(k),ptot,rro,cp,cv,agrad,q,u2)
        emu(k) = (1.38*rro*tt(k))/(1.67e-8*pg(k))
      end do  
      
! Save to file
      open(unit=33,file='marcs2drift.dat')
      write(33,'(a2)') ' !'
      write(33,'(a43)')' ! MARCS output file to be read in by DRIFT'
      if(idust .eq. 0) then
        write(33,'(a2)') ' ! Dust not included'
      else
        write(33,'(a2)') ' ! Dust included'
      end if
      write(33,'(a2)') ' !'
      write(33,'(a32,a19)') ' ! Model parameters: Teff, logg,',
     *  ' mixing, overshoot:'
      write(33,'(f12.3,f13.3,f13.3,f13.3)') teff, log10(grav), 
     *  palfa, 2.200
      write(33,'(a2)') ' !'
      write(33,'(a31)') ' ! Number of atmosphere layers:'
      write(33,'(i5)') ntau
      write(33,'(a2)') ' !'
      write(33,'(a42)') ' ! Number of elements in abundances table:'
      write(33,'(i5)') nelem
      write(33,'(a2)') ' !'
      write(33,'(a32)') ' ! Z of the considered elements:'
      do i=1,nelem,8
        if(i .gt. nelem-8) then
          write(33,'(8(i5.2,1x))') (kelem(j), j=i,nelem)
        else
          write(33,'(8(i5,1x))') (kelem(j), j=i,i+7)
        end if
      end do      
      write(33,'(a2)') ' !'
      write(33,'(a5,6a16,a6,a16)') ' !  #', 'Rad [cm]', 
     *  'Temp [K]', 'Pgas [dyn cm-2]', 'Ro [g cm-3]', 'g [cm s-2]',
     *  'v_conv [cm s-1]', 'Flag', 'mu [amu]'
      do k=1,ntau
        write(33,'(i5,6e16.8,a6,e16.8)') k, rad(k), tt(k), pg(k),
     *  ro(k), surfgrav(k), v(k), flag(k), emu(k)
      end do
      write(33,'(a2)') ' !'
      write(33,'(a41)') ' ! Initial Element abundances for each Z:'
      do i=1,nelem-1,8
        if(i .gt. nelem-9) then
          write(33,'(8f6.2)') (abinit(j), j=i,nelem-1)
        else
          write(33,'(8f6.2)') (abinit(j), j=i,i+7)
        end if
      end do
      close(33)
      
      return 
      end

    
!-----------------------------------------------------------------------
! Saves model data to savemodel.dat
! Juncher 2015
!-----------------------------------------------------------------------
      subroutine savemodel
      
      implicit real*16 (a-h,o-z)
      include 'parameter.inc'
      
      character molname*4,osfil*60
      dimension rad(ndp), pg(ndp), v(ndp)
      common /cg/grav,konsg /cteff/teff,flux
      common /mixc/palfa,pbeta,pny,py
      common /cmolrat/ fold(ndp,8),molold,kl
      common/cxlset/xl(20,10),idum6,nl(10)
      common /cos/wnos(nwl),conos(ndp,nwl),wlos(nwl),wlstep(nwl)
     *    ,kos_step,nwtot,nosmol,newosatom,newosatomlist
     *    ,nchrom,osfil(30),molname(30),sampling
      common /cdustdata/ dabstable(max_lay,nwl),dscatable(max_lay,nwl),
     *    eps(max_eps,max_lay),temp(max_lay),pgas(max_lay),n_lay,idust
      common /tauc/tau(ndp),dtauln(ndp),jtau
      common /rossc/xkapr(ndp),cross(ndp)
      common /statec/ppr(ndp),ppt(ndp),pp(ndp),gg(ndp),zz(ndp),dd(ndp),
     & vv(ndp),ffc(ndp),ppe(ndp),tt(ndp),tauln(ndp),ro(ndp),ntau,iter
      common /cabinit/abinit(ndp,natms)
      common /tsuji/ nattsuji,nmotsuji,parptsuji(500),abtsuji(17,ndp)


! Absorption coefficients
      istan2 = 1
      jstan2 = nl(1)+1
      do k=1,ntau
        kl = k
        call absko(1,1,tt(k),ppe(k),istan2,jstan2,abska,sprida)
      end do
      
! Radius
      rad(1:ntau) = 0.
      do k=1,ntau
        kl = k
        furem = fure      
        call termo(tt(k),ppe(k),ppr(k),ptot,rro,cp,cv,agrad,q,u2)
        fure = 1./(xkapr(k)*rro)
        if(k .eq. 1) cycle
        rad(k) = rad(k-1) + (tau(k)-tau(k-1))*(fure+furem)*0.5
      end do
      
! Gas pressure
      do k=1,ntau
        kl = k
        call jon(tt(k),ppe(k),1,pgx,rox,dumx,0)
        pg(k) = PGx
      end do
      
! Convective velocity
      do k=1,ntau
        if(k .gt. 1) go to 13
        v(k) = 0.
        go to 15
        
13      if(k .eq. ntau) go to 14
        ya=(tau(k)-tau(k-1))/(tau(k+1)-tau(k-1))
        yb=1.-ya
        v(k)=ya*vv(k+1)+yb*vv(k)
        if(vv(k).gt.0..and.vv(k+1).gt.0.) v(k)=
     &    exp(ya*log(vv(k+1))+yb*log(vv(k)))
        go to 15
        
14      continue
        ya=(2.*tau(k)-tau(k-1)-tau(k-2))/(tau(k)-tau(k-2))
        yb=1.-ya
        v(k)=ya*vv(k)+yb*vv(k-1)
15      continue
      end do

! Write to file
      open(unit=18,file='savemodel.dat')
      
      write(18,'(a42)') '!-----------------------------------------'
      write(18,'(a38)') '! Model atmosphere computed with MARCS'
      write(18,'(a42)') '!-----------------------------------------'
      write(18,*)
      write(18,'(a18)') '! Model parameters'
      write(18,'(a4,2x,a1,2x,f9.1)')  'Teff', '=', teff
      write(18,'(a4,2x,a1,2x,f9.1)')  'logg', '=', log10(grav)
      write(18,'(a4,2x,a1,2x,f9.1)')  'Z/Z0', '=', 1.0
      write(18,'(a4,2x,a1,2x,e10.3)') 'Flux', '=', flux
      write(18,'(a4,2x,a1,2x,f9.1)')  'alfa', '=', palfa
      write(18,'(a4,2x,a1,2x,i9)')    'Ntau', '=', ntau
      write(18,'(a4,2x,a1,2x,i9)')    'Nelm', '=', natms
      write(18,'(a4,2x,a1,2x,i9)')    'Nmol', '=', nosmol
      if(idust .eq. 1) then
        write(18,'(a4,2x,a1,2x,i9)')   'Dust', '=', 7
      else
        write(18,'(a4,2x,a1,2x,i9)')   'Dust', '=', 0
      end if
      write(18,*)
      write(18,'(a31)') '! Molecules included in opacity'
      write(18,'(*(a5))') molname(1:nosmol)
      write(18,*)
      write(18,'(a30)') '! Model atmosphere (CGS units)'
      write(18,'(5x,2a12,a9,6a12)') 'TauRoss', 'Geom. depth',
     *  'T','Pe', 'Pg', 'Prad', 'Rho', 'KappaRoss', 'Conv. vel.'
      do k=1,ntau
        write(18,'(i5,1p2e12.4,0pf9.2,1p6e12.4)') k, tau(k), rad(k), 
     *     tt(k), ppe(k), pg(k), ppr(k), ro(k), xkapr(k), v(k)
      end do
      write(18,*)
      write(18,'(28a)') '! Initial element abundances'
      write(18,'(a6,16(1x,a5))') '!    H', 'He', 'C', 'N', 'O', 'Ne',
     *   'Na','Mg', 'Al', 'Si', 'S', 'K', 'Ca', 'Ti', 'Cr', 'Fe', 'Ni'
      write(18,'(17f6.2)') abinit(1,1:2), abinit(1,6:8), 
     *  abinit(1,10:14), abinit(1,16), abinit(1,19:20), 
     *  abinit(1,21), abinit(1,23), abinit(1,25), abinit(1,27)
      write(18,*)
      if(idust .eq. 1) then
        write(18,'(34a)') '! Final layer dependent abundances'
        write(18,'(7a6)') '!     ', 'O', 'Mg', 'Al', 'Si', 'Ti', 'Fe'
        do k=1,ntau,10
          write(18,'(i6,6f6.2)') k, log10(abtsuji(k,8))+12., 
     *      log10(abtsuji(k,11:13))+12., log10(abtsuji(k,20))+12.,
     *      log10(abtsuji(k,24))+12.
        end do
      end if
      close(18)

      return
      end

!-----------------------------------------------------------------------
! Intiates ggchem input
! ERC 2018i
!-----------------------------------------------------------------------
      subroutine init_ggchem_ERC(tt,ptot)
     
      implicit none
      real*16,intent(in) :: tt,ptot
      real*16 :: bar

      print*,' now in init_ggchem_ERC(tt,ptot)'
      bar=1.Q+6
      open(unit=70,file='marcs2ggchem.in', status='replace')
      write(70, '(19a)') '# selected elements'
      write(70, '(99a)') 'H He C N O Na Mg Si F Fe Al Ca Cr Ti S Cl K Li
     & V Zr el'
      write(70, '(38a)') '# name of files with molecular kp-data'
      write(70,'(a54)') 'dispol_BarklemCollet.dat             
     &   ! dispol_file '
      write(70,'(a54)')
     & 'dispol_StockKitzmann_withoutTsuji.dat   ! dispol_file2'
      write(70,'(a54)') 'dispol_WoitkeRefit.dat         
     &         ! dispol_file3'
      write(70,'(a64)') '# abundance options 1=EarthCrust, 2=Ocean,
     & 3=Solar, 4=Meteorites'
      write(70,'(a34)') '3                     ! abund_pick'
      write(70,'(a27)') '# equilibrium condensation?'
      write(70,'(a36)') '.false.               ! model_eqcond'
      write(70,'(a15)') '# model options'
      write(70,'(a42)') '1                     ! model_dim  (0,1,2)'
      write(70,'(a36)') '.true.                ! model_pconst'
!      open(unit=10,file='s',status='old')

!      write(70,"(F4.1,a26)") tt, '                ! Tmax [K]'
!      write(70,"(F4.1,a48)") tt-1,'                ! Tmin [K]     
      write(70,*) tt, '                ! Tmax [K]'
      write(70,*) tt-1,'                ! Tmin [K] 

     & (if model_dim>0)'
      write(70,*) ptot/bar,'                   ! pmax [bar]   
     & (if pconst=.true.)'
      write(70,*) ptot/bar,'                   ! pmin [bar]'
      write(70,'(a31)') '100                   ! Npoints'
      write(70,'(a33)') '5                     ! NewBackIt'
      write(70,'(a29)') '1000.0                ! Tfast'
      close(70)
!      print *,'Calling ggchem'
!      end

!-----------------------------------------------------------------------


      call ggchem_ERC
     

      return
      end
      


!-----------------------------------------------------------------------
! ggchem equilibrium computation
! ERC 2018
!-----------------------------------------------------------------------
      subroutine ggchem_ERC
!      
      implicit real*16 (a-h,o-z)

!     Partial pressures from GG-chem subroutine DEMO_SWEEP (ERC J=4)
      common /partialpressure/
     > ppN2,ppCH,ppCO,ppCN,ppC2,ppNO,ppC2H2,ppHCN,ppC2H,ppC3,ppCS,ppH2O,
     > ppOH,ppTiO,ppSiO,ppCH4,ppNH,ppSiH,ppFeH,ppVO,ppZrO,ppMgH,ppNH3,
     > ppCO2,ppTiH,ppCaH,ppCrH,ppLiH,ppH,ppO2,ppHm,ppH2,ppH2p,ppHS,
     > ppC3H,ppSiC,ppSiC2,ppNS,ppSiN,ppSO,ppS2,ppSiS,ppLaO,ppCH2,
     > ppCH3,ppSi2C,ppSiO2,ppH2S,ppCaOH,ppCHNO,ppSiF2

      real*16 ::
     > ppN2,ppCH,ppCO,ppCN,ppC2,ppNO,ppC2H2,ppHCN,ppC2H,ppC3,ppCS,ppH2O,
     > ppOH,ppTiO,ppSiO,ppCH4,ppNH,ppSiH,ppFeH,ppVO,ppZrO,ppMgH,ppNH3,
     > ppCO2,ppTiH,ppCaH,ppCrH,ppLiH,ppH,ppO2,ppHm,ppH2,ppH2p,ppHS,
     > ppC3H,ppSiC,ppSiC2,ppNS,ppSiN,ppSO,ppS2,ppSiS,ppLaO,ppCH2,
     > ppCH3,ppSi2C,ppSiO2,ppH2S,ppCaOH,ppCHNO,ppSiF2

      include 'parameter.inc'

      !print *, 'Calling GGchem'
      call system('./GGchem/ggchem marcs2ggchem.in')


      open(unit=777,file='PartialPressures.dat')
      read(777,*) ppN2,ppCH,ppCO,ppCN,ppC2,ppNO,ppC2H2,ppHCN,ppC2H,ppC3,
     & ppCS,ppH2O,ppOH,ppTiO,ppSiO,ppCH4,ppNH,ppSiH,ppFeH,ppVO,
     & ppZrO,ppMgH,ppNH3,ppCO2,ppTiH,ppCaH,ppCrH,ppLiH,ppH,ppO2,
     & ppHm,ppH2,ppH2p,ppHS,ppC3H,ppSiC,ppSiC2,ppNS,ppSiN,ppSO,
     & ppS2,ppSiS,ppLaO,ppCH2,ppCH3,ppSi2C,ppSiO2,ppH2S,ppCaOH,
     & ppCHNO,ppSiF2
      close(777)





      open(unit=778,file='test.dat',status='replace')

      write(778,*)ppN2,ppCH,ppCO,ppCN,ppC2,ppNO,ppC2H2,ppHCN,ppC2H,ppC3,
     & ppCS,ppH2O,ppOH,ppTiO,ppSiO,ppCH4,ppNH,ppSiH,ppFeH,ppVO,
     & ppZrO,ppMgH,ppNH3,ppCO2,ppTiH,ppCaH,ppCrH,ppLiH,ppH,ppO2,
     & ppHm,ppH2,ppH2p,ppHS,ppC3H,ppSiC,ppSiC2,ppNS,ppSiN,ppSO,
     & ppS2,ppSiS,ppLaO,ppCH2,ppCH3,ppSi2C,ppSiO2,ppH2S,ppCaOH,
     & ppCHNO,ppSiF2
      close(778)



      return
      end



