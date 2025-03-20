C     ERZEUGUNG DES FILES "ODFOUT" FUER PC !!!!
C  -------------------
      subroutine cko2odf (ckostr, outstr, ODF_new)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      character*255 outstr, ckostr
      double precision ODF_new(19,19,19)
C  -------------------
C
C     Calculates ODF from texture components (the C-coefficients) read
C     from file 'ckostr'.  The result is written to file 'outstr' and
C     also returned in 'ODF_new'.
C
C
C     Subroutinen-Version
C     (c) Okt. 2001, O. Engler
C
C     EINGABE  =  "ckostr*60"
C     AUSGAGE  =  "outstr*60"
C
C  --------------------
c      use msflib
      
      LOGICAL LCKN,LOGA,LOGB
      CHARACTER*1 FORM
      CHARACTER*4 YPRO
      CHARACTER*5 ODF1,ODF2,ODF3
      CHARACTER*7 XLEV
c      character*255 outstr, ckostr, act_drive, lib_drive
      character*255 lib_drive
      CHARACTER*10 ACT1,ACT2,ACT3,ACT11
C  --------------------
      COMMON /SET/ INP,LIB1,LIB2,LIB3,LIB4,LIB5,IOUT
      COMMON /PAR/ INU,NPHS,NAS,NAB,NPH,NPS,KSL
      COMMON /PAR2/ FORM(28)
      COMMON /LLL/ LMIN,LMAX,LD,LMAY
      COMMON /FUN/ TLMN(16983)
      COMMON /FUL/ TSKL(16983)
      COMMON /FUK/ XK(34,3,8),YK(34,18,4),ALF(27),BET(27)
      COMMON /FUK2/ IHKL(27)
      COMMON /EXC/ CLMN(629)
      COMMON /MOC/ CMOD(629),FM
      COMMON /ACT/ ACT1(8),ACT2(8),ACT3(8),ACT11(8)
      COMMON /MMM/ XFI(54,3),XYZ(55,4),Q(27),ZN(27),
     $             WS2(27),WS3(27),WS4(27),WS5(27),WS0,YFI(27,3),H(27)
      COMMON /ODF1/ ODF1,ODF2,ODF3,YPRO,XLEV
      COMMON /ODF2/ YLEV(10), ODF(19,19,19)
      COMMON /LIN/ LCKN,LOGA,LOGB
      common /drives/ lib_drive


C  --------------------
      monitor=1
      INP  = 26
      LIB1 = 28
      LIB2 = 29
      LIB3 = 36
      LIB4 = 37
      LIB5 = 38
      IOUT = 27
C  -----
c      i4=getdrivedirqq (act_drive)                  ! store actual directory 
c      i4=changedirqq (lib_drive)                        ! change to library-directory 

      OPEN (UNIT=inp,FILE=ckostr)

      OPEN (UNIT=lib1,FILE='libmln.dat',ERR=9010, STATUS='OLD')
      OPEN (UNIT=lib2,FILE='libxk.dat', ERR=9010, STATUS='OLD')
      OPEN (UNIT=lib3,FILE='liba1.dat', ERR=9010, STATUS='OLD')
      OPEN (UNIT=lib4,FILE='libbu.dat', ERR=9010, STATUS='OLD')
      OPEN (UNIT=lib5,FILE='liba23.dat',ERR=9010, STATUS='OLD')
      OPEN (UNIT=iout,FILE=outstr)
c      i4=changedirqq (act_drive)                        ! return to actual directory 
C  -----
      LD   = 1
      LMIN = 1
C  -----
      CALL ATK
      IF ( LCKN ) CALL CKN
      CALL ODFS
C   -----
      CLOSE (inp)
      CLOSE (iout)
      CLOSE (lib1)
      CLOSE (lib2)
      CLOSE (lib3)
      CLOSE (lib4)
      CLOSE (lib5)
C  --------------------
      do 8999 i=1,19
        do 8999 j=1,19
          do 8999 k=1,19
              ODF_new(i,j,k)=odf(i,j,k)
8999      continue
      return

c 9010 message = messageboxqq 
c     #('at least one of the libraries LIBMLN,LIBXK,LIBA1,LIBBU,LIBA23 is
c     #& missing'C, 'Error'C, mb$iconquestion)
 9010 continue
      print *, 'at least one of the libraries LIBMLN,LIBXK,LIBA1,',
     $     'LIBBU,LIBA23 is missing'
      STOP 
      END
C
C ******************************************************************************
C
      SUBROUTINE ATK
C  --------------------
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL LCKN,LOGA,LOGB
      CHARACTER*1 FORM
c      CHARACTER*4 YPRO,XJAHR
      CHARACTER*4 YPRO
      CHARACTER*6 XLACT,XWS0
      CHARACTER*5 XDAT,XTIM,ODF1,ODF2,ODF3
c      CHARACTER*8 XMOD,XMINI,DATUM,ZEIT,XALL,XNPH
      CHARACTER*8 XMOD,XMINI,DATUM,ZEIT,XNPH
      CHARACTER*10 ACT1,ACT2,ACT3,ACT11
      CHARACTER XNAB*12,XNPS*16,XKSL*14,XLEV*7
C  --------------------
c      DIMENSION JDT(3),MO(34)
      DIMENSION MO(34)
      DIMENSION FIX(162),FIY(81)
      DIMENSION XC(18,3)
      DIMENSION IDATTIM(8)
C  --------------------
      COMMON /SET/ INP,LIB1,LIB2,LIB3,LIB4,LIB5,IOUT
      COMMON /PAR/ INU,NPHS,NAS,NAB,NPH,NPS,KSL
      COMMON /PAR2/ FORM(28)
      COMMON /LLL/ LMIN,LMAX,LD,LMAY
      COMMON /FUN/ TLMN(16983)
      COMMON /FUL/ TSKL(16983)
      COMMON /FUK/ XK(34,3,8),YK(34,18,4),ALF(27),BET(27)
      COMMON /FUK2/ IHKL(27)
      COMMON /ACT/ ACT1(8),ACT2(8),ACT3(8),ACT11(8)
      COMMON /MMM/ XFI(54,3),XYZ(55,4),Q(27),ZN(27),
     $             WS2(27),WS3(27),WS4(27),WS5(27),WS0,YFI(27,3),H(27)
      COMMON /ODF1/ ODF1,ODF2,ODF3,YPRO,XLEV
      COMMON /ODF2/ YLEV(10), ODF(19,19,19)
      COMMON /LIN/ LCKN,LOGA,LOGB
      COMMON /EXC/ CLMN(629)
      COMMON /OUT1/ XMOD,XMINI,XDAT,DATUM,XTIM,ZEIT,XLACT
      COMMON /OUT2/ NMOD
C  --------------------
      DATA MO /0,0,0,1,0,1,0,1,1,1,0,2,1,1,1,2,1,2,1,2,2,2,1,3,2,2,2,3,
     $         2,3,2,3,3,3/
C  --------------------
      LCKN  = .TRUE.
      PI = DACOS(-1.0D0)
      XB = PI / 180.D0
C  -----
      READ(INP,7000) ACT1
      READ(INP,7011) NMOD,XMOD,XMINI,XDAT,DATUM,XTIM,ZEIT,XLACT,LMAX
      IF ( XDAT.EQ.'     ' ) XDAT = 'DATE '
      IF ( XTIM.EQ.'     ' ) XTIM = 'TIME '
      IF ( NMOD.EQ.0 ) LCKN = .FALSE.
C  -----
c     CALL DATIM(IDATTIM)
C     WRITE(XJAHR,9100) IDATTIM(8)
c     IDATTIM(8) = IDATTIM(8) - 1900
c      call idate (idattim(7), idattim(6), idattim(8))
      call date_and_time(VALUES=idattim)
      idattim(1) = mod(idattim(1), 100)
      WRITE(DATUM,9200) IDATTIM(3),IDATTIM(2),IDATTIM(1)
      WRITE(ZEIT,9210) IDATTIM(5),IDATTIM(6),IDATTIM(7)
c      call time (zeit)
c 9100 FORMAT(A4)
 9200 FORMAT(I2,'.',I2,'.',I2)
 9210 FORMAT(I2,':',I2,':',I2)
C  -----
      IF ( LCKN ) THEN
        READ(INP,7001) NAB,XNAB,NPH,XNPH,NPS,XNPS,KSL,XKSL,XWS0,WS0
        INU = NAB + NPH + NPS
        IF ( INU.GT.27 ) STOP ' TOO MANY COMPONENTS'
        READ(INP,7000) ACT2
        NAP = NAB + NPH
        IK = 1
        DO 11 IT=1,INU
          IT2 = 2 * IT
          IF ( IT.LE.NAP ) READ(INP,7002) FORM(IT),(XFI(IK,K),K=1,3),
     $                     Q(IT),(XYZ(IT2,KK),KK=1,4),
     $                     (XYZ(IT2+1,LL),LL=1,4),
     $                     ZN(IK),WS2(IT),WS3(IT),WS4(IT)
          IF ( IT.LE.NAB.OR.IT.GT.NAP ) GOTO 10
            IP = IT - NAB
            READ(INP,7015) IHKL(IP),ALF(IP),BET(IP)
 10       CONTINUE
          IF ( IT.LE.NAP ) GOTO 20
            READ(INP,7002) FORM(IT),(XFI(IK+1,K),K=1,3),Q(IT),
     $                     (XYZ(IT2,KK),KK=1,4),(XYZ(IT2+1,LL),LL=1,4),
     $                     ZN(IK+1),WS2(IT),WS3(IT),WS4(IT)
            IP = IT - NAB
            READ(INP,7015) IHKL(IP),ALF(IP),BET(IP)
            READ(INP,7016) (XFI(IK,K),K=1,3),ZN(IK)
            IK = IK + 1
 20       CONTINUE
          IK = IK + 1
 11     CONTINUE
        READ(INP,7003) (XYZ(1,KK),KK=1,4)
        DO 12 IT=1,KSL
          READ(INP,7014)(YFI(IT,KK),KK=1,3),XWS5,WS5(IT)
 12     CONTINUE
        READ(INP,7000) ACT3
      END IF
C  -----
      READ(INP,7004) ODF1,ODF2,ODF3,YPRO,XLEV,(YLEV(K),K=1,10)
      IF ( (ODF1.EQ.'FUNE ').AND.(ODF2.EQ.'     ').AND.
     $       (ODF3.EQ.'     ') ) LCKN = .FALSE.
C
      IF ( .NOT.LCKN ) THEN
        ODF2 = '     '
        ODF3 = '     '
      END IF
C
      IF ( (ODF1.EQ.'     ').AND.(ODF2.EQ.'     ').AND.
     $       (ODF3.EQ.'     ') ) STOP ' OUTPUT DEFINITION FEHLT !!!'
C
      IF ((ODF1.NE.'FUNE '.AND.ODF1.NE.'REDU '.AND.ODF1.NE.'MODE '.AND.
     $     ODF1.NE.'TRUE '.AND.ODF1.NE.'DIFF '.AND.ODF1.NE.'GHOST'.AND.
     $     ODF1.NE.'TRANS'.AND.ODF1.NE.'GAUSS').AND.
     $   (ODF2.NE.'FUNE '.AND.ODF2.NE.'REDU '.AND.ODF2.NE.'MODE '.AND.
     $   ODF2.NE.'TRUE '.AND.ODF2.NE.'DIFF '.AND.ODF2.NE.'GHOST' ).AND.
     $   (ODF3.NE.'FUNE '.AND.ODF3.NE.'REDU '.AND.ODF3.NE.'MODE '.AND.
     $   ODF3.NE.'TRUE '.AND.ODF3.NE.'DIFF '.AND.ODF3.NE.'GHOST' ))THEN
C       WRITE(monitor,9000)
        ODF1 = 'FUNE '
        ODF2 = '     '
        ODF3 = '     '
      END IF
C
      READ(INP,7005) ACT11,LMAY
      IF ( LMAX.GT.LMAY ) LMAX = LMAY
C  -----
      IF ( LCKN ) THEN
        DO 21 KK=1,4
          XYZ(1,KK) = XYZ(1,KK) * 0.01D0
 21     CONTINUE
        DO 22 IT=1,INU
          IT2 = 2 * IT
          DO 22 KK=1,4
            XYZ(IT2,KK) = XYZ(IT2,KK) * 0.01D0
 22     CONTINUE
C  -----
        NAS = NAB + NPH + 2 * NPS
C  -----
        J = 0
        DO 5 IT=1,NAS
          DO 5 IX=1,3
            J = J + 1
            FIX(J) = XFI(IT,IX) * XB
 5      CONTINUE
C  -----
        J = 0
        DO 6 IT=1,KSL
          DO 6 IX=1,3
            J = J + 1
            FIY(J) = YFI(IT,IX) * XB
 6      CONTINUE
        CALL FUNT(FIX,FIY)
C
        NPHS = NPH + NPS
        IF ( NPHS.EQ.0 )  GOTO 3
          CALL FUNK                                                    !  PHASE
 3      CONTINUE
      END IF
C  -----
        LOGA = XMINI.EQ.'ALLC    '
        LOGB = XMINI.NE.'ALLC    '
        IF ( LOGA ) LD = 1
        IF ( LOGB ) LD = 2
C
        J = 0
        LLD = 1
        IF ( LOGB ) LLD = 2
        DO 1 L=4,LMAX,LLD
          LTL = MO(L)
          IF ( LTL.EQ.0 ) GOTO 1
          LP1 = L / 2 + 1
          DO 2 M1=1,LTL
            READ(INP,1002) (XC(N,M1),N=1,LP1)
 2        CONTINUE
          LPL = MO(L+1)
          IF ( L.EQ.LMAX ) GOTO 44
          IF ( LOGA.OR.LPL.EQ.0 ) GOTO 44
          L1 = L + 1
 44       CONTINUE
          DO 40 N=1,LP1
            DO 40 M1=1,LTL
              J = J + 1
              CLMN(J) = XC(N,M1)
 40       CONTINUE
          IF ( LPL.NE.0.AND.LLD.EQ.2 ) J = J + LPL * LP1
 1      CONTINUE
C  --------------------
 1002 FORMAT(8F9.4)
 7000 FORMAT(8A10)
 7011 FORMAT(I2,A8,9X,A8,12X,A5,A8,4X,A5,A8,1X,A6,I2)
 7001 FORMAT(I2,A12,1X,I2,A8,1X,I2,A16,1X,I2,A14,3X,A6,F4.0)
 7002 FORMAT(A1,3F5.1,F6.2,1X,F6.2,3F4.1,F5.1,3F4.1,F4.0,1X,3F5.2)
 7003 FORMAT(1X,15X,6X,1X,F6.2,3F4.1)
 7004 FORMAT(3A5,10X,A4,2X,A7,10F4.0)
 7005 FORMAT(7A10,A6,I2)
 7014 FORMAT(1X,3F5.1,48X,A5,F7.2)
 7015 FORMAT(I5,13X,F5.1,6X,F5.1)
 7016 FORMAT(1X,3F5.1,42X,F4.0)
c 8000 FORMAT(I10,A10,5F10.3)
c 9000 FORMAT(' falshe Ouputkontrol Definition !!!',/,
c     $       ' drei Angaben sind mglich von den volgenden ',
c     $       '(A5 Format):',/,
c     $       ' FUNE          Fune      -ODF',/,
c     $       ' REDU          Red.Modell-ODF',/,
c     $       ' MODE          Modell    -ODF',/,
c     $       ' TRUE          True      -ODF',/,
c     $       ' DIFF          Diff.     -ODF',/,
c     $       ' GHOST         Ghost     -ODF' )

C  --------------------
      END
C
C ******************************************************************************
      SUBROUTINE FUNT(XFI,SFI)
C  --------------------
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL KLOP
      LOGICAL LOGD,LOGP
      LOGICAL LOGMAX
      LOGICAL LOL,LOS
C  --------------------
      DIMENSION EOS(36)
      DIMENSION MO(35)
      DIMENSION XFI(162)
      DIMENSION SFI(81)
      DIMENSION T(20,18,3)
      DIMENSION BU(34,3,9)
      DIMENSION D(9,18,35),AA(40),A(40)
C  --------------------
      COMMON /SET/ INP,LIB1,LIB2,LIB3,LIB4,LIB5,IOUT
      COMMON /PAR/ INU,NPHS,NAS,NAB,NPH,NPS,KSL
      COMMON /LLL/ LMIN,LMAX,LD,LMAY
      COMMON /FUN/ TLMN(16983)
      COMMON /FUL/ TSKL(16983)
C  --------------------
      DATA MO/0,0,0,1,0,1,0,1,1,1,0,2,1,1,1,2,1,2,1,2,2,2,1
     $       ,3,2,2,2,3,2,3,2,3,3,3,2/
C  --------------------
      REWIND LIB4
      REWIND LIB5
      READ(LIB4,*) BU
      LOGD = .FALSE.
      LOGP = .FALSE.
      I = 1
      K = 1
      DO 4 L=1,LMAX,1
        KLOP = .TRUE.
        LOL = (L / 2 * 2).EQ.L
        LM2 = L / 4 + 1
        LP1 = L / 2 + 1
        LTL = MO(L)
        LOGMAX = L.NE.LMAX
        IF ( LTL.EQ.0 ) GOTO 44
        LMA2 = L + 1
        DO 5 M=1,LM2
          DO 5 N=1,LP1
            READ(LIB5,*) A
            DO 5 IS=1,LMA2
 5      D(M,N,IS) = A(IS)
        INX = NAS
        GOTO 200
 100    CONTINUE
        INX = KSL
 200    CONTINUE
        JJ = 1
        DO 6 J=1,INX
          IF ( KLOP ) GOTO 202
            X3 = 2.D0 * SFI(JJ)
            Z3 = SFI(JJ+1)
            Y3 = 4.D0 * SFI(JJ+2)
            GOTO 101
 202      CONTINUE
            X3 = 2.D0 * XFI(JJ)
            Z3 = XFI(JJ+1)
            Y3 = 4.D0 * XFI(JJ+2)
 101      CONTINUE
          JJ = JJ + 3
          DO 61 IS=1,LMA2
            EOS(IS) = DCOS(DFLOAT(IS - 1) * Z3)
 61       CONTINUE
          DO 7 N=1,LP1
            XNOR = 3.5449077D0
            IF ( N.EQ.1 ) XNOR = 2.50662827D0
            YXN = DFLOAT(N - 1) * X3
            SX = DSIN(YXN)
            CX = DCOS(YXN)
            DO 77 M1=1,LTL
              Q = .0D0
              Z = .0D0
              DO 71 M=1,LM2
                YXM = DFLOAT(M - 1) * Y3
                SY = DSIN(YXM)
                CY = DCOS(YXM)
                BBB = BU(L,M1,M)
                X = .0D0
                Y = .0D0
                DO 8 IS=1,LMA2,2
                  LOS = IS.EQ.LMA2
                  IS1 = IS + 1
                  X = X + D(M,N,IS) * EOS(IS)
                  IF ( LOS.AND.LOL ) GOTO 8
                  Y = Y + D(M,N,IS1) * EOS(IS1)
 8              CONTINUE
                Q = Q + X * CY * BBB
                Z = Z + Y * SY * BBB
 71           CONTINUE
              T(J,N,M1) = XNOR * (CX * Q - SX * Z)
 77         CONTINUE
 7        CONTINUE
 6      CONTINUE
        DO 11 N=1,LP1
          DO 11 M1=1,LTL
            DO 11 J=1,INX
              IF ( KLOP )      TLMN(I) = T(J,N,M1)
              IF ( .NOT.KLOP ) TSKL(K) = T(J,N,M1)
              IF ( KLOP )      I = I + 1
              IF ( .NOT.KLOP ) K = K + 1
 11     CONTINUE
        IF ( .NOT.KLOP ) GOTO 300
          KLOP = .FALSE.
          GOTO 100
 300    CONTINUE
        IF ( LOGD.AND.LOGMAX ) GOTO 55
          GOTO 4
 44     CONTINUE
        DO 45 M=1,LM2
          DO 45 N=1,LP1
 45     READ(LIB5,*) AA
        IF ( LOGD.AND.LOGMAX ) GOTO 55
        GOTO 4
 55     CONTINUE
        LN = L + 1
        LN2 = LN / 2 + 1
        LN4 = LN / 4 + 1
        DO 66  M=1,LN4
          DO 66  N=1,LN2
 66     READ(LIB5,*) AA
 4    CONTINUE
      REWIND LIB4
      REWIND LIB5
C  --------------------
      RETURN
      END
C
C ******************************************************************************
      SUBROUTINE FUNK
C  --------------------
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C  --------------------
      DIMENSION A(40),AA(40)
C  --------------------
      COMMON /SET/ INP,LIB1,LIB2,LIB3,LIB4,LIB5,IOUT
      COMMON /PAR/ INU,NPHS,NAS,NAB,NPH,NPS,KSL
      COMMON /LLL/ LMIN,LMAX,LD,LMAY
      COMMON /FUK/ XK(34,3,8),YK(34,18,4),ALF(27),BET(27)
      COMMON /FUK2/ IHKL(27)
C  --------------------
      DATA RSPI /3.5449077018D0/
      DATA XB /0.017453293D0/
C  --------------------
      DO 2 I=1,NPHS
        DO 2 L=1,34
          DO 2 N=1,18
 2    YK(L,N,I) = 0.0D0
C  -----
       READ(LIB2,*) XK
       DO 33 K=1,5
         READ(LIB3,*) AA
 33    CONTINUE
C  -----
      XNN = 1.0D0 / RSPI * DSQRT(2.0D0)
C  -----
      DO 10 L=4,LMAX,2
        LP1 = L / 2 + 1
        DO 11 N=1,LP1
          READ(LIB3,*) A
          DO 101 I=1,NPHS
            XALF = ALF(I) * XB
            XBET = BET(I) * XB
            X = 0.0D0
            DO 20 IS=1,LP1
              X = X + A(IS) * DCOS(2.0D0 * DFLOAT(IS - 1) * XALF)
 20         CONTINUE
            EPS = 1.0D0
            IF ( N.NE.1 ) EPS = DSQRT(2.0D0)
            YK(L,N,I) = X * EPS * DCOS(2.0D0 * DFLOAT(N - 1) * XBET)*XNN
 101      CONTINUE
 11     CONTINUE
        DO 12 N=1,LP1
          READ(LIB3,*) AA
 12     CONTINUE
 10   CONTINUE
      REWIND LIB2
      REWIND LIB3
C  --------------------
      RETURN
      END
C
C ******************************************************************************
      SUBROUTINE CKN
C  --------------------
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL LOPH
      LOGICAL LOGA,LOGB,LCKN
      LOGICAL LOWP,LOWN
      LOGICAL LOGP
      LOGICAL LOGS,LOGN
      CHARACTER*8 WPER,W234,WSET
      CHARACTER*1 GAUS,BESL,XLOR,PHAS,PAS1,XFUN,FORM
C  --------------------
      DIMENSION MODHKL(8)
      DIMENSION XMLN(18,18)
      DIMENSION MO(34)
      DIMENSION ERCO(34),SRCO(34)                                      ! PHASE
      DIMENSION NRI(27)
      DIMENSION S0(27),S2(27),S3(27),S4(27),F2(27)
      DIMENSION S5(27)
      DIMENSION XLIT(34,27),YLIT(34,27)                                ! PAS1
C  --------------------
      COMMON /SET/ INP,LIB1,LIB2,LIB3,LIB4,LIB5,IOUT
      COMMON /PAR/ INU,NPHS,NAS,NAB,NPH,NPS,KSL                        ! PHASE
      COMMON /PAR2/ FORM(28)
      COMMON /LLL/ LMIN,LMAX,LD,LMAY
      COMMON /FUN/ TLMN(16983)
      COMMON /FUL/ TSKL(16983)
      COMMON /FUK/ XK(34,3,8),YK(34,18,4),ALF(27),BET(27)
      COMMON /FUK2/ IHKL(27)
      COMMON /MOC/ CMOD(629),FM
      COMMON /EXC/ CLMN(629)
      COMMON /MMM/ XFI(54,3),XYZ(55,4),Q(27),ZN(27),
     $             WS2(27),WS3(27),WS4(27),WS5(27),WS0,YFI(27,3),H(27)
      COMMON /LIN/ LCKN,LOGA,LOGB
C  --------------------
      DATA XB   /0.017453293D0 /
      DATA RSPI /3.5449077018D0/
      DATA MO /0,0,0,1,0,1,0,1,1,1,0,2,1,1,1,2,1,2,1,2,2,2,1,3,2,2,2,3,
     $         2,3,2,3,3,3/
      DATA MODHKL /100,110,111,102,112,122,103,113/
      DATA WPER /'WPER    '/
      DATA W234 /'W234    '/
      DATA WSET /'W234    '/
      DATA GAUS /'G'/
      DATA BESL /'B'/
      DATA XLOR /'L'/
      DATA PHAS /'F'/
      DATA PAS1 /'I'/
C  --------------------
        LOWN = WSET.EQ.W234
        LOWP = WSET.EQ.WPER
C
        NPAR = 2 * INU + 1
        LMP1 = LMAX / 2 + 1
        IF ( LOGA ) LD = 1
        IF ( LOGB ) LD = 2
        LLD = 1
        IF ( LOGB ) LLD = 2
C  -----
        DO 6 L=1,LMAX
          ER = 1.0D0 / (2.0D0 * DFLOAT(L) + 1.0D0)
          SR = RSPI * DSQRT(ER)                    
          ERCO(L) = ER
          SRCO(L) = SR                             
 6      CONTINUE
C  -----
        DO 14 IT=1,INU
          Q(IT) = 1.0D0
 14     CONTINUE
        DO 16 IK=1,KSL
          H(IK) = 1.0D0
 16     CONTINUE
        I = 0
        K = 0
        KS = 0
        NABP = NAB + NPH
        NABP1 = NABP + 1
        DO 18 L=4,LMAX,1
          LTL = MO(L)
          IF ( LTL.EQ.0 ) GOTO 18
            LP1 = L / 2 + 1
            LOGP = (L / 2 * 2).EQ.L
            LOGS = LOGA.OR.LOGP
            DO 24 N=1,LP1
              DO 24 M1=1,LTL
                I = I + 1
                IF ( LOGS ) CI = CLMN(I)
                DO 25 IT=1,NABP
                  K = K + 1
                  IF ( LOGS ) Q(IT) = Q(IT) + CI * TLMN(K)
 25             CONTINUE
                IF ( INU.EQ.NABP ) GOTO 27
                  DO 26 IT=NABP1,INU
                    K = K + 2
                    IF ( LOGS ) Q(IT) = Q(IT) + CI * TLMN(K)
 26               CONTINUE
 27             CONTINUE
                DO 55 IK=1,KSL
                  KS = KS + 1
                  IF ( LOGS ) H(IK) = H(IK) + CI * TSKL(KS)
 55             CONTINUE
 24         CONTINUE
 18     CONTINUE
C  -----
        DO 7 K1=1,NPHS                             
          NHKL = IHKL(K1)                          
          DO 8 K2=1,8                              
            IF ( MODHKL(K2).EQ.NHKL ) GOTO 9       
 8        CONTINUE                                 
 9        NRI(K1) = K2                             
 7      CONTINUE                                   
C  -----
        INU2 = 2 * INU
        INU2 = 2 * INU + 1
C  -----
        CON1 = 0.125D0 * RSPI                      
        CON2 = 4.0D0 * DSQRT(2.0D0) / RSPI         
C  -----
        READ(LIB1,*) XMLN                          
        REWIND LIB1                                
c 100  CONTINUE
C  -----
      FM = XYZ(1,1)
      DO 5 IT=2,INU2,2
        FM = FM + XYZ(IT,1)
 5    CONTINUE
C  -----
      XLN = DLOG(2.0D0)
      YLN = DSQRT(XLN)
      ZR = 2.D0 * YLN * XB
      XR = -XB * XB / 4.0D0
C  -----
      DO 10 IT=1,INU
        XFUN = FORM(IT)
        IT1 = 2 * IT
        IT2 = 2 * IT + 1
        XI1 = XYZ(IT1,1)
        XI2 = XYZ(IT2,1)
        IF ( XFUN.EQ.BESL ) GOTO 400
        IF ( XFUN.EQ.XLOR ) GOTO 500
        IF ( XFUN.EQ.PHAS ) GOTO 600               
C  -----
        Z = XR * XI2 * XI2
        SN = RSPI * XI1 / (XI2 * XB * ZN(IT) * (1.0D0 - DEXP(Z)))
        SNG = SN                                      
        IF ( XFUN.EQ.GAUS ) S0(IT) = SN               
        DO 220 L=4,LMAX,1
          XL = DFLOAT(L)
          XL1 = XL + 1.0D0
          XLL = XL * XL
          XL11 = XL1 * XL1
          XLIT(L,IT) = XI1 * (DEXP(Z * XLL) - DEXP(Z * XL11)) /
     $                 (1.0D0 - DEXP(Z))
 220    CONTINUE
        IF ( XFUN.EQ.PAS1 ) GOTO 600                  
        GOTO 404
C  -----
 400    CONTINUE
        STOP ' STANDARD GAUSS FORM IS NOT ALLOWED !!!'
 500    CONTINUE
        B  = ZR * XI2
        CB = DCOS(.25D0 * B)
        C  = CB * CB
        C2 = C * C
        C3 = C2 * C
        U = DSQRT(19.D0 * C2 - 34.D0 * C + 19.D0)
        TAU = (2.D0 * U * DCOS(DACOS((-82.D0 * C3 + 240.D0 *
     $         C2 - 246.D0 * C + 80.D0) / U**3) / 3.D0)
     $         + 5.D0 - 4.D0 * C) / 3.D0
        TE = DSQRT(TAU) - DSQRT(TAU - 1.D0)
        TE2 = TE * TE
        SN = XI1 / ZN(IT) * ((1.D0 + TE2)**2 + 4.D0 * TE2) /
     $       (1.D0 - TE2)**3
        S0(IT) = SN
        DO 202 L=4,LMAX,1
          L2 = 2 * L
          XL21 = DFLOAT(L2) + 1.D0
          XLIT(L,IT) = XI1 * XL21 * TE**L2
 202    CONTINUE
        GOTO 404                                   
C  -----
 600    CONTINUE                                   
        Z = XR * XI2 * XI2                         
        XI22 = XI2 / 2.0D0 * XB                    
        XIER = ERF(XI22)                           
        ITT = IT                                                        !!!!
        IF ( XFUN.EQ.PAS1 ) ITT = ITT + 1                               !!!!
        SN = XI1 / (ZN(ITT) * CON1 * XI2 * XB * DEXP(Z) * XIER)         !!!!
        SNP = SN                                      
        S0(IT) = SN                                
        XRR = CON2 * XI1 / XIER                    
        DO 204 L=4,LMAX,2                          
          LN = L / 2 + 1                           
          LNMAX = LMAX / 2 + 1                     
          Y = 0.0D0                                
          DO 206 NN=LN,LNMAX                       
            EPS = 1.0D0                            
            IF ( NN.NE.1 ) EPS = 2.0D0             
            N = 2 * (NN - 1)                       
            Y = Y + XMLN(LN,NN) * DEXP(DFLOAT(N * N - 1) * Z) * EPS     ! PHASE
 206      CONTINUE                                 
          YLIT(L,IT) = Y * XRR * SRCO(L)              
 204    CONTINUE                                   
C  -----
        IF ( XFUN.NE.PAS1 ) GOTO 404                  
        R = -SNP / SNG                                
        RR = 1.0D0 / (1.0D0 + R)                      
        RRR = R * RR                                  
C  -----
 404    CONTINUE
C
        S2(IT) = Q(IT) - FM
        S3(IT) = FM - XI1 - XYZ(1,1)
        S4(IT) = SN - XI1
 10   CONTINUE
C  -----
      DO 20 IK=1,KSL
        S5(IK) = H(IK) - FM
        H(IK) = FM
 20   CONTINUE
C  -----
      F0 = 1.0D0 - FM
      F1 = F0 * F0 * WS0
      F = 0.0D0
      I = 0
      J = 0
      K = 0
      KS = 0
      DO 15 L=4,LMAX,1
        LTL = MO(L)                                
        IF ( LTL.EQ.0 ) GOTO 15
        LOGP = (L / 2 * 2).EQ.L
        LOGS = LOGA.OR.LOGP
        LOGN = .NOT.LOGS
        LP1 = L / 2 + 1
        YY = 0.0D0
        DO 45 N=1,LP1                              
          DO 45 M1=1,LTL
            I = I + 1
            KP = 1                                 
            Y = 0.0D0
            DO 60 IT=1,INU
              XFUN = FORM(IT)                      
              LOPH = XFUN.EQ.PHAS                  
              IF ( LOPH ) GOTO 62                  
                J = J + 1                          
                T = TLMN(J)                        
                XLTL = XLIT(L,IT) * T              
                IF ( XFUN.EQ.PAS1 ) GOTO 62           
                F2(IT) = XLTL * T                  
                GOTO 66                            
 62           CONTINUE                             
              NRH = NRI(KP)                        
              J = J + 1                            
              YLTL = 0.0D0                            
              T1 = XK(L,M1,NRH)                    
              T2 = YK(L,N,KP)                                                                                    ! 2
              T = TLMN(J)                          
              IF ( LOGP ) YLTL = YLIT(L,IT) * T1 * T2 
              F2(IT) = YLTL * T                       
              KP = KP + 1                          
 66           CONTINUE                             
              IF ( XFUN.EQ.PAS1 ) Y = Y + RR * YLTL + RRR * XLTL        
              IF ( XFUN.EQ.PHAS ) Y = Y + YLTL     
              IF ( XFUN.NE.PHAS.AND.XFUN.NE.PAS1 ) Y = Y + XLTL
 60         CONTINUE
            DO 75 IT=1,INU
              XFUN = FORM(IT)                         
              K = K + 1
              IF ( XFUN.EQ.PAS1 ) K = K + 1           
              FK = Y * TLMN(K)
              FI = F2(IT)
              IF ( LOGS ) S2(IT) = S2(IT) - FK
              S3(IT) = S3(IT) + (FK - FI)
              S4(IT) = S4(IT) - FI
 75         CONTINUE
C  -----
            DO 80 IK=1,KSL
              KS = KS + 1
              FS = Y * TSKL(KS)
              IF ( LOGS ) S5(IK) = S5(IK) - FS
              IF ( LOGS ) H(IK) = H(IK) + FS
 80         CONTINUE
            DQ = 0.0D0
            IF ( LOGN ) GOTO 77
              DC = CLMN(I) - Y
              DQ = DC * DC
 77         CONTINUE
            YY = YY + DQ
            CMOD(I) = Y
 45     CONTINUE
        F = F + YY * ERCO(L)
 15   CONTINUE
C  --------------------
c 1002 FORMAT(8F9.4)
C  --------------------
      RETURN
      END
C
C ******************************************************************************
      SUBROUTINE ODFS
C  --------------------
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL LOGA,LOGB,LCKN
      LOGICAL LOG1,LOG2,LOG3,LOG4,LOG5,LOG6,LOGTR
      CHARACTER*4  YPRO
      CHARACTER*5  FUNY,REDY,MODY,TRUY,DIFY,GHOY,ODF1,ODF2,ODF3
      CHARACTER*5  XTRA,XGAU
      CHARACTER*5  XDAT,XTIM
      CHARACTER*6  XLACT
      CHARACTER*7  XLEV
      CHARACTER*8  XMOD,XMINI,DATUM,ZEIT
      CHARACTER*10 XAFU,XAMO,XAGH,XADI,XATR,XARE,XAME
      CHARACTER*10 ACT1,ACT2,ACT3,ACT11
C  --------------------
      DIMENSION MO(34)
C  --------------------
      COMMON /SET/ INP,LIB1,LIB2,LIB3,LIB4,LIB5,IOUT
      COMMON /LLL/ LMIN,LMAX,LD,LMAY
      COMMON /XAM/ XAME
      COMMON /ACT/ ACT1(8),ACT2(8),ACT3(8),ACT11(8)
      COMMON /EXC/ CLMN(629)
      COMMON /MOC/ CMOD(629),FM
      COMMON /ODF1/ ODF1,ODF2,ODF3,YPRO,XLEV
      COMMON /ODF2/ YLEV(10), ODF(19,19,19)
      COMMON /EVD/ C(34,3,18),NDF
      COMMON /LIN/ LCKN,LOGA,LOGB
      COMMON /DDF/ DF
      COMMON /OUT1/ XMOD,XMINI,XDAT,DATUM,XTIM,ZEIT,XLACT
      COMMON /OUT2/ NMOD
      COMMON /LOGT/ LOGTR
C  --------------------
      DATA MO/0,0,0,1,0,1,0,1,1,1,0,2,1,1,1,2,1,2,1,2,2,2,1,3,2,2,2,3,2,
     $        3,2,3,3,3/
      DATA FUNY /'FUNE '/, REDY /'REDU '/, MODY /'MODEL'/
      DATA TRUY /'TRUE '/, DIFY /'DIFF '/, GHOY /'GHOST'/
      DATA XTRA /'TRANS'/, XGAU /'GAUSS'/
      DATA XAFU /'FUNE   ODF'/
      DATA XAMO /'MODELL ODF'/
      DATA XAGH /'GHOST  ODF'/
      DATA XADI /'DIFF.  ODF'/
      DATA XATR /'TRUE   ODF'/
      DATA XARE /'RED.MODELL'/
C  --------------------
      NDF = 1
      LOG1 = .TRUE.
      LOG2 = .TRUE.
      LOG3 = .TRUE.
      LOG4 = .TRUE.
      LOG5 = .TRUE.
      LOG6 = .TRUE.
      LOGTR = .TRUE.
C  -----
C       EXPERIMENTAL   ODF
C       C-EXP   PACK
C
      XAME = XAFU
      LMIN = 2
      IF ( LOGA ) LMIN = 1
      LD = 2
      IF ( LOGA ) LD = 1
      J = 0
      DO 5 L=4,LMAX,LD
        LTL = MO(L)
        IF ( LTL.EQ.0 ) GOTO 5
        LP1 = L / 2 + 1
        DO 50 N=1,LP1
          DO 50 M1=1,LTL
            J = J + 1
            C(L,M1,N) = CLMN(J)
 50     CONTINUE
        IF ( LOGA ) GOTO 55
          J = J + MO(L+1) * LP1
 55     CONTINUE
 5    CONTINUE
C
C             C-EXP   PRINT
C
      IF ( (ODF1.EQ.FUNY.OR.ODF2.EQ.FUNY.OR.ODF3.EQ.FUNY.OR.
     $      ODF1.EQ.XTRA.OR.ODF1.EQ.XGAU).AND.LOG1 ) THEN
        DFM = 1.0D0
        WRITE(IOUT,7000) ACT1
        WRITE(IOUT,7011) NMOD,XMOD,XMINI,XDAT,DATUM,XTIM,ZEIT,XLACT,LMAX
        IF ( ODF1.EQ.FUNY.OR.ODF2.EQ.FUNY.OR.ODF3.EQ.FUNY ) THEN
          WRITE(IOUT,1001) YPRO,YLEV
          IF ( .NOT.LOGA ) WRITE(IOUT,1002) ACT1(1),ACT1(2),DFM,LMAY
          IF ( LOGA ) WRITE(IOUT,1022) ACT1(1),ACT1(2),DFM,LMAY
        ELSE IF ( ODF1.EQ.XTRA ) THEN
          WRITE(IOUT,1101) YPRO,YLEV
          IF ( .NOT.LOGA ) WRITE(IOUT,1102) ACT1(1),ACT1(2),DFM,LMAY
          IF ( LOGA ) WRITE(IOUT,1112) ACT1(1),ACT1(2),DFM,LMAY
        ELSE IF ( ODF1.EQ.XGAU ) THEN
          WRITE(IOUT,1201) YPRO,YLEV
          WRITE(IOUT,1202) ACT1(1),ACT1(2),DFM,LMAY
        END IF
        CALL CCCC(DFM,YPRO)
        LOG1 = .FALSE.
      END IF
C  -----
C       COMPLETE   MODEL   ODF
C       C-MOD   PACK
C
      I = 0
      DO 1 L=4,LMAX,1
        LTL = MO(L)
        IF ( LTL.EQ.0 ) GOTO 1
        LP1 = L / 2 + 1
        DO 10 N=1,LP1
          DO 10 M1=1,LTL
            I = I + 1
            C(L,M1,N) = CMOD(I)
 10     CONTINUE
 1    CONTINUE
C
C       C-MOD   PRINT
C
      IF ( (ODF1.EQ.MODY.OR.ODF2.EQ.MODY.OR.ODF3.EQ.MODY)
     $     .AND.LOG2 ) THEN
        DFM = FM
        XAME = XAMO
        LD = 1
        LMIN = 1
        WRITE(IOUT,7000) ACT1
        WRITE(IOUT,7011) NMOD,XMOD,XMINI,XDAT,DATUM,XTIM,ZEIT,XLACT,LMAX
        WRITE(IOUT,1003) YPRO,YLEV
        WRITE(IOUT,1004) ACT1(1),ACT1(2),DFM,LMAY
        CALL CCCC(DFM,YPRO)
        LOG2 = .FALSE.
      END IF
C  -----
C       REDUCED   MODEL   ODF
C
      IF ( (ODF1.EQ.REDY.OR.ODF2.EQ.REDY.OR.ODF3.EQ.REDY)
     $     .AND.LOG3 ) THEN
        DFM = FM
        XAME = XARE
        LD = 2
        LMIN = 2
        WRITE(IOUT,7000) ACT1
        WRITE(IOUT,7011) NMOD,XMOD,XMINI,XDAT,DATUM,XTIM,ZEIT,XLACT,LMAX
        WRITE(IOUT,1005) YPRO,YLEV
        WRITE(IOUT,1006) ACT1(1),ACT1(2),DFM,LMAY
        CALL CCCC(DFM,YPRO)
        LOG3 = .FALSE.
      END IF
C  -----
C       TRUE   ODF
C
      IF ( (ODF1.EQ.TRUY.OR.ODF2.EQ.TRUY.OR.ODF3.EQ.TRUY.OR.
     $     ODF1.EQ.DIFY.OR.ODF2.EQ.DIFY.OR.ODF3.EQ.DIFY).AND.
     $     LOG4 ) THEN
        XAME = XATR
        LD = 1
        LMIN = 1
        J = 0
        DO 2 L=4,LMAX,2
          LTL = MO(L)
          IF ( LTL.EQ.0 ) GOTO 2
          LP1 = L / 2 + 1
          DO 20 N=1,LP1
            DO 20 M1=1,LTL
              J = J + 1
              C(L,M1,N) = CLMN(J)
 20       CONTINUE
          J = J + MO(L+1) * LP1
 2      CONTINUE
C
C       C-MIX   PRINT
C
        LOGTR = .FALSE.
        DFM = 1.0D0
        IF ( ODF1.EQ.TRUY.OR.ODF2.EQ.TRUY.OR.ODF3.EQ.TRUY ) THEN
          LOGTR = .TRUE.
          WRITE(IOUT,7000) ACT1
          WRITE(IOUT,7011) NMOD, XMOD, XMINI, XDAT, DATUM, XTIM, ZEIT, 
     &                                    XLACT, LMAX
          WRITE(IOUT,1007) YPRO,YLEV
          WRITE(IOUT,1008) ACT1(1),ACT1(2),DFM,LMAY
        END IF
        CALL CCCC(DFM,YPRO)
        LOG4 = .FALSE.
        LOGTR = .TRUE.
      END IF
C  -----
C       DIFFERENCES   ODF   /EXP-MOD/
C
      IF ( (ODF1.EQ.DIFY.OR.ODF2.EQ.DIFY.OR.ODF3.EQ.DIFY)
     $     .AND.LOG5 ) THEN
        XAME = XADI
        LMIN = 2
        IF ( LOGA ) LMIN = 1
        LD = 2
        IF ( LOGA ) LD = 1
        J = 0
        DO 4 L=4,LMAX,LD
          LTL = MO(L)
          IF ( LTL.EQ.0 ) GOTO 4
          LP1 = L / 2 + 1
          DO 40 N=1,LP1
            DO 40 M1=1,LTL
              J = J + 1
              C(L,M1,N) = CLMN(J) - CMOD(J)
 40       CONTINUE
          IF ( LOGA ) GOTO 44
            J = J + MO(L+1) * LP1
 44       CONTINUE
 4      CONTINUE
C
C       C-DIFF  PRINT
C
        DFM = 1.0D0 - FM
        WRITE(IOUT,7000) ACT1
        WRITE(IOUT,7011) NMOD,XMOD,XMINI,XDAT,DATUM,XTIM,ZEIT,XLACT,LMAX
        WRITE(IOUT,1009) YPRO,YLEV
        WRITE(IOUT,1010) ACT1(1),ACT1(2),DFM,LMAY
        CALL CCCC(DFM,YPRO)
        LOG5 = .FALSE.
      END IF
C  -----
C       GHOST   ODF
C
      IF ( (ODF1.EQ.GHOY.OR.ODF2.EQ.GHOY.OR.ODF3.EQ.GHOY)
     $     .AND.LOG6 ) THEN
        XAME = XAGH
        V = 1.0
        K = 0
        DO 3 L=4,LMAX,1
          LTL = MO(L)
          IF ( LTL.EQ.0 ) GOTO 33
            LP1 = L / 2 + 1
            DO 30 N=1,LP1
              DO 30 M1=1,LTL
                K = K + 1
                C(L,M1,N) = 0.5D0 * V * CMOD(K)
 30         CONTINUE
 33       CONTINUE
          V = -V
 3      CONTINUE
C
C       C-GHOST   PRINT
C
        DFM = 0.5D0
        WRITE(IOUT,7000) ACT1
        WRITE(IOUT,7011) NMOD,XMOD,XMINI,XDAT,DATUM,XTIM,ZEIT,XLACT,LMAX
        WRITE(IOUT,1011) YPRO,YLEV
        WRITE(IOUT,1012) ACT1(1),ACT1(2),DFM,LMAY
        CALL CCCC(DFM,YPRO)
        LOG6 = .FALSE.
      END IF
c 100  CONTINUE
      IF ( .NOT.LOG4.AND..NOT.LOG5 ) WRITE(IOUT,2000) DF
C  --------------------
 1001 FORMAT('FUNE         ODF',9X,A4,2X,'LEVELS:',10F4.0)
 1002 FORMAT(2A10,6X,'EVEN(EXP) C-COEFFICIENTS',2X,'C1=',F6.4,9X,
     $       'LMAX= ',I2)
 1022 FORMAT(2A10,6X,'EVEN+ODD  C-COEFFICIENTS',2X,'C1=',F6.4,9X,
     $       'LMAX= ',I2)
 1101 FORMAT('TRANS        ODF',9X,A4,2X,'LEVELS:',10F4.0)
 1102 FORMAT(2A10,6X,'EVEN(TRANS)   C-COEFFIC.',2X,'C1=',F6.4,9X,
     $       'LMAX= ',I2)
 1112 FORMAT(2A10,6X,'EVEN+ODD(TRANS) C-COEFF.',2X,'C1=',F6.4,9X,
     $       'LMAX= ',I2)
 1201 FORMAT('GAUSS        ODF',9X,A4,2X,'LEVELS:',10F4.0)
 1202 FORMAT(2A10,6X,'EVEN+ODD(GAUSS) C-COEFF.',2X,'C1=',F6.4,9X,
     $       'LMAX= ',I2)
 1003 FORMAT('COMPL. MODEL ODF',9X,A4,2X,'LEVELS:',10F4.0)
 1004 FORMAT(2A10,6X,'EVEN+ODD(MOD) C-COEFFIC.',2X,'C1=',F6.4,9X,
     $       'LMAX= ',I2)
 1005 FORMAT('RED. MODEL   ODF',9X,A4,2X,'LEVELS:',10F4.0)
 1006 FORMAT(2A10,6X,'EVEN(MOD) C-COEFFICIENTS',2X,'C1=',F6.4,9X,
     $       'LMAX= ',I2)
 1007 FORMAT('TRUE         ODF',9X,A4,2X,'LEVELS:',10F4.0)
 1008 FORMAT(2A10,6X,'EVEN(EXP)+ODD(MOD)   CKO',2X,'C1=',F6.4,9X,
     $       'LMAX= ',I2)
 1009 FORMAT('DIFF         ODF',9X,A4,2X,'LEVELS:',10F4.0)
 1010 FORMAT(2A10,6X,'EVEN(EXP)-EVEN(MOD)  CKO',2X,'C1=',F6.4,9X,
     $       'LMAX= ',I2)
 1011 FORMAT('GHOST        ODF',9X,A4,2X,'LEVELS:',10F4.0)
 1012 FORMAT(2A10,6X,'ODD(MOD)  C-COEFFICIENTS',2X,'C1=',F6.4,9X,
     $       'LMAX= ',I2)
 2000 FORMAT(//,10X,'AF = ',F5.1,' %')
 7000 FORMAT(8A10)
 7011 FORMAT(I2,A8,9X,A8,12X,A5,A8,4X,A5,A8,1X,A6,I2)
C  --------------------
      END
C
C ******************************************************************************
      SUBROUTINE CCCC(Z0,YPRO)
C  --------------------
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL LOGD,LOGP,LOGTR
      LOGICAL LOGMAX
      CHARACTER*4 XPF1,XPF2,YPRO
      CHARACTER*10 ACT1,ACT2,ACT3,ACT11,XAME
C  --------------------
      DIMENSION E(18,36),G(18),H(18),SLN(631)
      DIMENSION XP(19,19),MO(35),W(20),AA(40),A(40)
      DIMENSION D(9,18,35)
      DIMENSION BU(34,3,9)
C  --------------------
      COMMON /SET/ INP,LIB1,LIB2,LIB3,LIB4,LIB5,IOUT
      COMMON /ACT/ ACT1(8),ACT2(8),ACT3(8),ACT11(8)
      COMMON /XAM/ XAME
      COMMON /LLL/ LMIN,LMAX,LD,LMAY
      COMMON /EVD/ C(34,3,18),NDF
      COMMON /ODF2/ YLEV(10), ODF(19,19,19)
      COMMON /DDF/ DF
      COMMON /LOGT/ LOGTR
C  --------------------
      DATA MO/0,0,0,1,0,1,0,1,1,1,0,2,1,1,1,2,1,2,1,2,2,2,1,3,2,2,2,3,2
     $        ,3,2,3,3,3,2/
      DATA XPF1 /'PHI1'/
      DATA XPF2 /'PHI2'/
      DATA B /0.087266465D0/
C  --------------------
      XLF1 = 00.0D0
      XLF2 = 90.0D0
      LF1 = INT((XLF1 + .0001D0) / 5.D0) + 1
      LF2 = INT((XLF2 + .0001D0) / 5.D0) + 1
      IF ( NDF.EQ.1 ) THEN
        DO 3 J=1,631
          SLN(J) = DSIN(DFLOAT(J - 1) * B)
 3      CONTINUE
      END IF
      LMT0 = LMAX / 4 + 1
      LMT1 = LMAX / 2 + 1
      LMT2 = LMAX + 1
      DO 1 MM=1,LMT0
        DO 1 N=1,LMT1
          DO 1 IS=1,LMT2
 1    D(MM,N,IS) = 0.0D0
      DO 20 IFI=2,20
        XX = DFLOAT(IFI - 2) * 5.D0
        W(IFI) = XX
 20   CONTINUE
C  -----
      REWIND LIB4
      REWIND LIB5
      READ(LIB4,*) BU
      LOGD = LD.EQ.2
      LOGP = LMIN.EQ.2
      IF ( LOGP ) READ(LIB5,*) AA
      DO 7 L=LMIN,LMAX,LD
        LM4 = L / 4 + 1
        LP1 = L / 2 + 1
        LTL = MO(L)
        LOGMAX = L.NE.LMAX
        IF ( LTL.EQ.0 ) GOTO 77
        LT2 = L + 1
        DO 8 MM=1,LM4
          DO 8 N=1,LP1
            X = 0.D0
            DO 9 M1=1,LTL
              X = X + C(L,M1,N) * BU(L,M1,MM)
 9          CONTINUE
            IF ( N.NE.1 ) X = X * 1.4142135D0
            READ (LIB5,*) A
            DO 10 IS=1,LT2
              D(MM,N,IS) = D(MM,N,IS) + X * A(IS)
 10         CONTINUE
 8      CONTINUE
        IF ( LOGD.AND.LOGMAX ) GOTO 55
        GOTO 7
 77     CONTINUE
        DO 88 MM=1,LM4
          DO 88  N=1,LP1
 88     READ(LIB5,*) AA
        IF ( LOGD.AND.LOGMAX ) GOTO 55
        GOTO 7
 55     CONTINUE
        LN = L + 1
        LN2 = LN / 2 + 1
        LN4 = LN / 4 + 1
        DO 66 MM=1,LN4
          DO 66  N=1,LN2
          READ(LIB5,*) AA
 66     CONTINUE
 7    CONTINUE
      REWIND LIB4
      REWIND LIB5
C  -----
      IF ( YPRO.EQ.XPF2 ) THEN
        LD2 = LMAX / 4 + 1
        IFIP = 2
        LMA1 = LMAX / 2 + 1
        IPS = 4
        JPS = 2
      ELSE
        LD2 = LMAX / 2 + 1
        IFIP = 1
        LMA1 = LMAX / 4 + 1
        IPS = 2
        JPS = 4
       END IF
      LMM =  2 * LMAX
      FMAX = 0.D0
      DO 100 IF2=LF1,LF2
        FM   = 0.D0
        W(1) = 5.D0 * DFLOAT(IF2 - 1)
        DO 11 N=1,LMA1
          IV = 1
          DO 12 IS=1,LMT2
            X = 0.D0
            I = 9 + 9 * IV
            DO 13 M=1,LD2
              IF ( YPRO.EQ.XPF2 ) THEN
                K = M
                J = N
              ELSE
                K = N
                J = M
              END IF
              LH = IPS * (M - 1) * (IF2 - 1) + 1 + I
              X = X + D(K,J,IS) * SLN(LH)
 13         CONTINUE
            E(N,IS) = X
            IV = -IV
 12       CONTINUE
 11     CONTINUE
        DO 14 IF=1,19
          DO 15 N=1,LMA1
            X = 0.D0
            Y = E(N,1)
            DO 16 IS=2,LMT2,2
              LH = (IS - 1) * (IF - 1) + 19
              X = X + E(N,IS) * SLN(LH)
              IF ( IS.EQ.LMT2 ) GOTO 16
              LH1 = IS * (IF - 1) + 19
              Y = Y + E(N,IS+1) * SLN(LH1)
 16         CONTINUE
            G(N) = Y
            H(N) = X
 15       CONTINUE
          DO 17 IF1=1,19
            Z = 0.0D0
            DO 18 N=1,LMA1
              LH = JPS * (N - 1) * (IF1 - 1) + 1
              Z = Z + G(N) * SLN(LH+18) - H(N) * SLN(LH)
 18         CONTINUE
            XP(IF,IF1) = Z0 + 2.50662827D0 * Z
            IF ( DABS(XP(IF,IF1)).GT.DABS(FM) ) FM = XP(IF,IF1)
            IF ( DABS(FM).GT.DABS(FMAX) ) FMAX = FM
 17       CONTINUE
 14     CONTINUE
C  -----
        IF ( LOGTR ) THEN
          WRITE(IOUT,1004) W(1),YPRO,FMAX,FM,YPRO,W(1)
          WRITE(IOUT,1005) (W(IF1),IF1=2,20)
          WRITE(IOUT,1006) (W(IF+1),(XP(IF,IF1),IF1=1,19),
     $                      W(IF+1),IF=1,19)
c          WRITE(13,123) ((XP(IF,IF1),IF1=1,19),IF=1,19)
c123            FORMAT(19F7.2)
            do if=1,19
              do if1=1,19
                  ODF(if1,if,if2)=XP(IF,IF1)
              enddo
          enddo
        END IF
C  -----
C
 100  CONTINUE
      IF ( XAME.EQ.'TRUE   ODF' ) FMEXP = FMAX
      IF ( XAME.EQ.'DIFF.  ODF'.AND.FMEXP.GT.0.D0 ) THEN
        DF = 100.D0 * FMAX / FMEXP
      END IF
C  --------------------
c1004 FORMAT(/,F5.0,'  DEGREE  -  ',A4,'-SECTION',45X,'( FMAX ABS =',
c    $       F7.2,4X,'FMAX-I =',F6.2,' )',4X,A4,' =',F4.0)
 1004 FORMAT(/,F5.0,'  DEGREE  -  ',A4,'-SECTION',t77,'( FMAX ABS=',
     $       F7.2,4X,'FMAX-I =',F7.2,' )',t121,A4,'=',F4.0)
c1005 FORMAT('------- ',19F6.2,' -------',/)
 1005 FORMAT('-------',19F6.2,' -------',/)
 1006 FORMAT(F5.0,2X,19F6.2,F7.0)
C  --------------------
      RETURN
      END
C
cC ******************************************************************************
c      DOUBLE PRECISION FUNCTION ERF(XX)
cC  --------------------
c      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
cC  --------------------
c      DIMENSION P1(4),Q1(4),P2(6),Q2(6),P3(4),Q3(4)
cC  --------------------
c      DATA (P1(I),I=1,4) /242.6679552305318D0,21.97926161829415D0,
c     1         6.9963834881619136D0,-3.560984370181539D-02/,
c     2     (Q1(I),I=1,4) /215.0588758698612D0,91.16490540451490D0,
c     3         15.082797630407739D0,1.0D0/,
c     4     (P2(I),I=1,6) /22.898992851659D0,26.094746956075D0,
c     5         14.571898596926D0,4.26577201070898D0,
c     6         .56437160686381D0,-6.0858151959688D-06/,
c     7     (Q2(I),I=1,6) /22.898985749891D0,51.933570687552D0,
c     8         50.273202863803D0,26.2788795758761D0,
c     9         7.5688482293618D0,1.0D0/,
c     1     (P3(I),I=1,4) /-1.21308276389978D-2,-.1199039552681460D0,
c     2         -.2439110294988626D0,-3.24319519277746D-2/,
c     3     (Q3(I),I=1,4) /4.30026643452770D-02,.489552441961437D0,
c     4         1.437712279371218D0,1.0D0/
c      DATA SQPI /.564189583547756D0/
cC  --------------------
c7      X = DABS(XX)
c      IF ( X.GT.6.0D0 ) GOTO 320
c        X2 = X * X
c      IF ( X.GT.4.0D0 ) GOTO 300
c      IF ( X.GT..46875D0 ) GOTO 200
c        A = X * (P1(1) + X2 * (P1(2) + X2 * (P1(3) + X2 * P1(4))))
c        A = A / (Q1(1) + X2 * (Q1(2) + X2 * (Q1(3) + X2 * Q1(4))))
c        IF ( XX.LT.0.D0 ) A = -A
c        ERF = A
c      GOTO 400
c 200    A = DEXP(-X2) * (P2(1) + X * (P2(2) + X * (P2(3) + X * (P2(4) +
c     $      X * (P2(5) + X * P2(6))))))
c        A = A / (Q2(1) + X * (Q2(2) + X * (Q2(3) + X * (Q2(4) + X *
c     $      (Q2(5) + X * Q2(6))))))
c        ERF = DSIGN((1.0D0 - A),XX)
c      GOTO 400
c 300  XI2 = 1.D0 / X2
c        R = XI2 * (P3(1) + XI2 * (P3(2) + XI2 * (P3(3) + XI2 * P3(4))))
c     $      / (Q3(1) + XI2 * (Q3(2) + XI2 * (Q3(3) + XI2 * Q3(4))))
c        A = DEXP(-X2) * (SQPI + R) / X
c        ERF = DSIGN((1.0D0 - A),XX)
c      GOTO 400
c 320  CONTINUE
c        ERF = XX / X
c 400  CONTINUE
cC  --------------------
c      RETURN
c      END
c
c
C ******************************************************************************
c      subroutine datim ( idat )
c      integer idat
c      idat = 0
c      return
c      end
c
C ******************************************************************************
