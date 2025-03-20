C     *************************************************************************
      subroutine GAUSS_ODF (inp, oristr, ckostr)
      IMPLICIT DOUBLE PRECISION ( A - H, O - Z )
      integer inp
      character*255 oristr, ckostr
C  -------------------
C
C     Calculates texture components (the C-coefficients) from Euler
C     angles.
C
C     Arguments
C     ---------
C       inp    (in)   Unit number of 'oristr'.
C       oristr (in)   Name of orientation .ori file with a short header 
C                     followed by the Euler angles.
C       ckostr (in)   Name of output file to which the C-coefficients are
C                     written.
C
C
C     Format of ori-file
C     ==================
C     <label>
C     PHI2 <lmax>
C     <NEuler>, <weight>, <psi0>
C     <angles>
C     ...
C     
C     Description
C     -----------
C     <label>    A label for the data (only the 10 first characters are read)
C     <lmax>     Order of expansion
C     <NEuler>   Number of Euler angles
C     <weight>   Weight, seems not to be used
C     <psi0>     Bin size bin in degrees
C     <angles>   List of Euler angles, weight and optionally phi0. If 'psi0' 
C                on the line above is zero, its actual value follows the 
C                weight.  Otherwise only 4 numbers are required per line.
C                Hence, each of the 'NEuler' lines should be of the form:
C                    phi1  Phi  phi2  weight  [phi0]
C
C  -------------------
C     BERECHNET GAUSS-MODELLE AUS TEXTURKOMPONENTEN
C     DURCH REIHENENTWICKLUNG MIT GERADEN UND UNGERADEN
C     REIHENENTWICKLUNGSKOEFFIZIENTEN
C     BIBLIOTHEKEN : LIB1="LIBBU"   LIB2="LIBA23"
C     1.ZEILE  UEBERSCHRIFT  ACT(8) FORMAT 8A10
C     2.ZEILE  SCHNITT,LMAX  FORMAT A4,1X,I2  (BLANK MOEGLICH)
C                     (DEFAULT SCHNITT='PHI2', LMAX=22)
C     3.ZEILE  ANZAHL KOMPONENTEN (MAX=240000), MI, PSI0   FREIES FORMAT  , WG
C                     (DEFAULT PSI0=5.; 0.1 = 0.) (KEIN BLANK MOEGLICH)
C     4.ZEILE  GAUSSKOMPONENTEN  (F1,F,F2,MI,PSI0)       FREIES FORMAT
C              MENGENANTEILE IN PROZENT (%) ANGEBEN
C  -------------------
C
C     EINGABE  =  "oristr*60"
C     AUSGAGE  =  "ckostr*60"
C
C  -------------------
C
C     Subroutinen-Version
C     (c) Okt. 2001, O. Engler
C
C  -------------------
c      use msflib

C  -------------------
      CHARACTER*4  XPRO
      CHARACTER*5  XDAT,XTIM
      CHARACTER*6  XLACT,XLMAX
      CHARACTER*7  XLEV
      CHARACTER*8  XMOD,XMINI,DATUM,ZEIT
      CHARACTER*10 label, name
      CHARACTER*16 XNAME
c      character*255 oristr, ckostr, act_drive, lib_drive
      character*255 lib_drive
C  -------------------
      COMMON /SET/ iout, LIB1,LIB2
      COMMON /XFF/ XFI (3,900000), WG (900000), PSI(900000)
      COMMON /LLL/ LMIN,LMAX,INU,LD
      COMMON /LL2/ XPRO
      COMMON /EVD/ C (35,3,18)
      COMMON /TAB/ BU (34,3,9)
      common /drives/ lib_drive
C  -------------------
      DIMENSION YLEV(10)
c      DIMENSION IDATTIM(8)
C  -------------------
      DATA XDAT  /'DATE '/
      DATA XTIM  /'TIME '/
      DATA XLACT /'LACT= '/
      DATA XLMAX /'LMAX= '/
      DATA XLEV  /'LEVELS:'/
      DATA XMINI /'ALLC    '/
      DATA XMOD  /'.MODEL  '/
      DATA XNAME /'GAUSS MODEL  ODF'/
      DATA YLEV  /2.D0,4.D0,7.D0,12.D0,20.D0,30.D0,40.D0,
     $            50.D0,60.D0,70.D0/

C  -----
c      in   = 26
      in   = inp
      LIB1 = 27
      LIB2 = 28
      iout = 29
C  -------------------
c      i4=getdrivedirqq (act_drive)                  ! store actual directory 
c      i4=changedirqq (lib_drive)                        ! change to library-directory 

      open (unit=iout, file=ckostr)
      OPEN ( UNIT=LIB1, ERR=9010, status='old', FILE='libbu.dat')
      OPEN ( UNIT=LIB2, ERR=9010, status='old', FILE='liba23.dat')
c      i4=changedirqq (act_drive)                        ! return to actual directory 
C  -----
      read(in,*) name
      read(in,2) lmax
2      format(4x,i9)
      read(in,*) inu, weight, psi0
      if (inu.gt.900000) goto 9001
      write(12,*) inu,' orientations'

      weight=0.d0
      if (psi0.lt.1d-6) then
        write(12,3) lmax
3        format('using actual psi0',' lmax=',i2)
        do i=1,inu
          read(in,*) (XFI(j,i), j=1,3), WG(i), PSI(i)
            weight=weight+wg(i)
        enddo
       else
        write(12,1) psi0, lmax
1        format ('using psi0=',f4.1,'° lmax=',i2)
        do i=1,inu
          read(in,*) (XFI(j,i), j=1,3), WG(i)
c
c            read(in,*) xfi(3,i), xfi(2,i), xfi(1,i)            ! KUL
c
c            wg(i)=1./inu
            psi(i)=psi0
            weight=weight+wg(i)
        enddo
      endif
      close (in)

      if (weight.lt.1.d-9) then            ! correct weigths with 1/inu
        do i=1,inu
          WG(i) = 1.d0/inu
        enddo
      endif

C  -----
      LMIN = 1
      LD = 1
C  -----
      DATUM = '        '
      ZEIT = '        '
c      LMAX = 22
      XPRO = 'PHI2'
C  -----
c      call idate (idattim(7), idattim(6), idattim(8))
c      WRITE(DATUM,9200) IDATTIM(6),IDATTIM(7),idattim(8)

c     WRITE(ZEIT,9210) IDATTIM(5),IDATTIM(4),IDATTIM(3)
c      call time (zeit)
c 9200 FORMAT(I2,'.',I2,'.',i2)
c 9210 FORMAT(I2,'.',I2,'.',I2)
C  -----
      NMOD = 0
      write(label,1001) name
      WRITE(iout,1000) label
      WRITE(iout,7011) NMOD,XMOD,XMINI,XDAT,DATUM,XTIM,ZEIT,XLACT,LMAX
      WRITE(iout,7004) XNAME,XPRO,XLEV,(YLEV(K),K=1,10)
      WRITE(iout,7005) label,XLMAX,LMAX
C  -----
      CALL IDAL
C  -------------------
c 1001 format(a<len_trim(name)>)
 1001 format(a10)
 1000 FORMAT(A10)
c 1010 FORMAT(A4,1X,I2)
 7004 FORMAT(A16,9X,A4,2X,A7,10F4.0)
 7005 FORMAT(A10,10X,'EVEN+ODD(GAUSS) C-COEFFICIENTS  C1=1.0000',
     $       9X,A6,I2)
 7011 FORMAT(I2,A8,9X,A8,12X,A5,A8,4X,A5,A8,1X,A6,I2)
c 7777 FORMAT(' PROBE : ',3A10)
C  -------------------
      CLOSE (iout)
      CLOSE (LIB1)
      CLOSE (LIB2)

C  -------------------
C     get rid of warning about unused dummy argument 'oristr'
      name(1:1) = oristr(1:1)
      return

c9001  message = messageboxqq 
c     &      ('too many orientations (max. 500000)'C, 'Error'C,
c     &      mb$iconquestion)
 9001 continue
      print *, 'too many orientations (max. 900000)'
      stop
c9010  message = messageboxqq 
c     &      ('Libraries LIBBU and/or LIBA23 missing'C, 'Error'C,
c     &      mb$iconquestion)
 9010 continue
      print *, 'Libraries LIBBU and/or LIBA23 missing'
      stop
      END
C
C  *************************************************************************
      SUBROUTINE OMPI(IPS)
C  -------------------
      IMPLICIT DOUBLE PRECISION ( A - H, O - Z )
C  -------------------
      CHARACTER *4  XPRO
C  -------------------
      COMMON /SET/ iout ,LIB1,LIB2
      COMMON /LLL/ LMIN,LMAX,INU,LD
      COMMON /LL2/ XPRO
      COMMON /EVD/ C (35,3,18)
C  -------------------
      DIMENSION MO (35)
C  -------------------
      DATA MO /0,0,0,1,0,1,0,1,1,1,0,2,1,1,1,2,1,2,1,2,2,2,1
     $        ,3,2,2,2,3,2,3,2,3,3,3,2/
C  -------------------
      V = -1.D0
      DO 1 L=LMIN,LMAX,1
         LTL = MO (L)
         IF ( LTL.EQ.0 ) GOTO 11
         LP1 = L / 2 + 1
         X2L = 2.D0 * DFLOAT (L)
         X2L1 = X2L + 1.D0
         DO 10 MI=1,LTL
         DO 10 NI=1,LP1
            IF ( IPS.EQ.0 ) C (L,MI,NI) = C(L,MI,NI)*5.D-1* V
            IF ( IPS.EQ.1 ) C (L,MI,NI) = C(L,MI,NI)*5.D-1* V/X2L1
            IF ( IPS.EQ.2 ) C (L,MI,NI) = C(L,MI,NI)*5.D-1* V*X2L/X2L1
   10    CONTINUE
   11    V = - V
    1 CONTINUE
C  -------------------
      RETURN
      END

C  *************************************************************************
      SUBROUTINE IDAL
C  -------------------
      IMPLICIT DOUBLE PRECISION ( A - H, O - Z )
C  -------------------
      LOGICAL LOGP
      LOGICAL LOGD
      LOGICAL LOGMAX
      LOGICAL LOL,LOS
      CHARACTER *4  XPRO
C  -------------------
      COMMON /SET/ iout, LIB1,LIB2
      COMMON /XFF/ XFI (3,900000), WG (900000), PSI(900000)
      COMMON /LLL/ LMIN,LMAX,INU,LD
      COMMON /LL2/ XPRO
      COMMON /EVD/ C (35,3,18)
      COMMON /TAB/ BU (34,3,9)
      COMMON D (9,18,36)
C  -------------------
      DIMENSION EOS(36)
      DIMENSION A(40)
      DIMENSION MO(35)
      DIMENSION XC(18,3)
C  -------------------
      DATA MO /0,0,0,1,0,1,0,1,1,1,0,2,1,1,1,2,1,2,1,2,2,2,1
     $        ,3,2,2,2,3,2,3,2,3,3,3,2/
      DATA XB /0.017453293D0/
C  -------------------
      XR = - XB * XB / 4.D0
      WR = 0.D0
C  -----
      DO 100 K=1,INU
         WR = WR + WG (K)
 100  CONTINUE

      XNN = 1.D0 / WR
      LMP1 = LMAX / 2 + 1

      REWIND LIB1
      REWIND LIB2

      READ ( LIB1, * )  BU
      LOGD = LD.EQ.2
      LOGP = LMIN.EQ.2
      IF ( LOGP )  READ ( LIB2, * )  A

C-----------  HAUPTSCHLEIFE   4   -----------------
      DO 4 L=LMIN,LMAX,LD
         iproz=99*(l-lmin)/(lmax-lmin)*(l-lmin)/(lmax-lmin)
c         call progress(iproz,25,550)
         LOL = L / 2 * 2.EQ.L
         LM2 = L / 4 + 1
         LP1 = L / 2 + 1
         LTL = MO (L)
         LOGMAX = L.NE.LMAX
         IF ( LTL.EQ.0 )  GOTO 44
         L2 = L
         LMA2 = L + 1

         DO 5 M=1,LM2
         DO 5 N=1,LP1
            READ  ( LIB2, * )  A
         DO 5 IS=1,LMA2
            D ( M, N, IS) = A (IS)
    5    CONTINUE

         DO 12 M1=1,LTL
         DO 12 N=1,LP1
            XC ( N, M1) = 0.D0
   12    CONTINUE

c         JJ = 1

C-----------  SCHLEIFE   6    --------------------
         DO 6 J=1,INU
            XIFI = WG(J)
            X3   = 2.D0 * XFI ( 1,J ) * XB
            Z3   = XFI ( 2,J ) * XB
            Y3   = 4.D0 * XFI ( 3,J ) * XB
            XLL  = DFLOAT ( 2 * L + 1 ) * XNN
            IF ( PSI(J).EQ.0.D0 )  GOTO 300
            XL   = DFLOAT (L)
            XL2  = XL + 1.D0
            XL21 = XL2 * XL2
            XL4  = XL * XL
            Z    = XR * PSI(J) * PSI(J)
            XXX1 = Z * XL4
            XXX2 = Z * XL21
            XXX3 = Z
            IF ( XXX1.LT.-179.D0 ) XXX1 = -179.D0
            IF ( XXX2.LT.-179.D0 ) XXX2 = -179.D0
            IF ( XXX3.LT.-179.D0 ) XXX3 = -179.D0
            XLEX = (DEXP (XXX1) - DEXP (XXX2) ) / ( 1.D0 - DEXP(XXX3))
            XLL  = XNN * XLEX

  300       CONTINUE
c           JJ = JJ + 3

            DO 61 IS=1,LMA2
               EOS(IS) = DCOS ( DFLOAT (IS-1) * Z3 )
   61       CONTINUE

C-------------   SCHLEIFE  7    -------------------
            DO 7 N=1,LP1
               YXN = DFLOAT ( N - 1 ) * X3
               SX  = DSIN (YXN)
               CX  = DCOS (YXN)
               DO 77 M1=1,LTL
                  Q = 0.D0
                  Z = 0.D0
                  DO 71 M=1,LM2
                    YXM = DFLOAT ( M - 1 ) * Y3
                    SY  = DSIN (YXM)
                    CY  = DCOS (YXM)
                    BBB = BU ( L, M1, M)
                    X = 0.D0
                    Y = 0.D0
                    DO 8 IS=1,LMA2,2
                        LOS = IS.EQ.LMA2
                        IS1 = IS + 1
                        X = X + D ( M, N, IS) * EOS(IS)
                        IF ( LOS.AND.LOL )  GOTO 8
                        Y = Y + D ( M, N, IS1) * EOS(IS1)
    8               CONTINUE
                    Q = Q + X * CY * BBB
                    Z = Z + Y * SY * BBB
   71             CONTINUE
                  XXXX = XIFI * XLL
                  YYYY = CX * Q - SX * Z
                  IF ( XXXX.NE.0.D0 ) THEN
                    IF ( (DLOG10(DABS(XXXX)) +
     $                    DLOG10(DABS(YYYY))).LT.-78.D0 ) XXXX = 0.D0
                  END IF
                  XC (N,M1) = XC (N,M1) + XXXX * YYYY
   77          CONTINUE

    7 CONTINUE
C-------------  SCHLEIFE  7 ENDE  -------------------
    6 CONTINUE
C-------------  SCHLEIFE  6 ENDE  -------------------

         DO 11 N=1,LP1
            XNOR = 3.5449077D0
            IF ( N.EQ.1 )  XNOR = 2.50662827D0
            DO 11 M1=1,LTL
               XC ( N, M1 ) = XC ( N, M1) * XNOR
   11    CONTINUE

         DO 78 M1=1,LTL
         DO 78 N=1,LP1
            C ( L, M1, N) = XC ( N, M1)
   78    CONTINUE

         DO 14 M1=1,LTL
            WRITE ( iout, 5000)  ( XC(N,M1), N=1,LP1)
   14    CONTINUE

         IF ( LOGD.AND.LOGMAX )  GOTO 55
         GOTO 4
   44    CONTINUE

         DO 45 M=1,LM2
         DO 45 N=1,LP1
            READ ( LIB2, * )  A
   45    CONTINUE

         IF ( LOGD.AND.LOGMAX )  GOTO 55
         GOTO 4

   55    LN  = L + 1
         LN2 = LN / 2 + 1
         LN4 = LN / 4 + 1

         DO 66 MM=1,LN4
         DO 66  N=1,LN2
            READ ( LIB2, * )  A
   66    CONTINUE

    4 CONTINUE
C-----------  HAUPTSCHLEIFE  ENDE  ---------------
      REWIND LIB1
      REWIND LIB2
C  -------------------
 5000 FORMAT(8F9.4)
C  -------------------
      RETURN
      END
C
C ******************************************************************
c      subroutine new_progress (ix,iy)
c
c      use msflib
c
c      i4=setcolorRGB(#00303030)            ! dark grey
c      i4=rectangle($GBorder,ix-2,iy-2,ix+218+2,iy+8+2)
c      i4=rectangle($GBorder,ix-1,iy-1,ix+218+1,iy+8+1)
c      i4=setcolorRGB(#00e0e0e0)            ! light grey
c      i4=rectangle($GFillInterior,ix,iy,ix+218,iy+8)
c      return
c      end
c
C  *************************************************************************
c      subroutine progress(iprozent,ix,iy)
c
c      use msflib
c      type (xycoord) xys
c
c      imax= int(20.*iprozent/100.+0.5)
c      i4=setcolorRGB(#900000)                ! dark blue
c      do i=1,imax
c        iix=ix+(i-1)*11
c        i4=rectangle($GFillInterior,iix,iy,iix+9,iy+9)
c      enddo
c
c      return
c      end
c
C ******************************************************************
