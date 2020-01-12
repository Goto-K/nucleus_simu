C*******************************************************
C*****                                             *****
C*****                MAIN PROGRAM                 *****
C*****                                             *****
C*****     RUTHERFORD SCATTERING CROSS SECTION     *****
C*****                                             *****
C*******************************************************
C
C=======================================================================
      PROGRAM MAIN
      IMPLICIT REAL*8(A-H,O-Z)
C=======================================================================
C
      PARAMETER (HC0=197.3269601D0,AC0=931.494013D0,EC0=137.03599976D0)
CCC---------------------------------------------------------------------
      HC=HC0
      AC=AC0
      EC=EC0
C
      PI=ACOS(-1.0D0)
      PIAD=PI/180.0D0
C
C
C-------------------
C---  FILE OPEN  ---
C-------------------
C
      CALL FOPEN
C
C
C-------------------------
C---  PARAMETER INPUT  ---
C-------------------------
C
      CALL INPUT(Z1,Z2,T,KIBO,THMIN,THMAX,DTH)
C
C
C---------------------
C---  CALCULATION  ---
C---------------------
C
      NTHMAX=(THMAX-THMIN)/DTH+1.01
C
      WRITE(KIBO,600)
C170605
  600 FORMAT(2X,'theta',4X,'cs(mb/sr)')
C170605

      A=Z1*Z2*HC/EC/T
      DO ITH=1,NTHMAX
C170514
C       THT=(ITH-1)*DTH
       THT=(ITH-1)*DTH+THMIN
C170514
       RAD=THT*PIAD
       IF(RAD.LT.1.0D-10) RAD=1.0D-10
C
       SIG=A**2/(16.0D0*SIN(RAD/2.0D0)**4.0D0)*10.0D0
C
       WRITE(10,601) THT,SIG     ! in unit of mb
  601  FORMAT(F8.4,E13.5)
C
      END DO
C
      STOP 0
      END
C
C
C
C=======================================================================
      SUBROUTINE INPUT(Z1,Z2,T,KIBO,THMIN,THMAX,DTH)
      IMPLICIT REAL*8 (A-H,O-Z)
C=======================================================================
C
      CHARACTER*50 COMMENT
CCC---------------------------------------------------------------------
C
      KIBAN=5
C
      READ(KIBAN,500) COMMENT
  500 FORMAT(A50)
      READ(KIBAN,501) Z1,Z2,T
  501 FORMAT(5F10.0)
C
      THMIN=0.0D0
      THMAX=180.0D0
      DTH=1.0D0
C
      READ(KIBAN,*)
      READ(KIBAN,501) FKIBO,THMIN,THMAX,DTH
C
      KIBO=FKIBO
C
C
C------------------------------
C---  OUTPUT OF THE INPUTS  ---
C------------------------------
C
      WRITE(6,600)
  600 FORMAT(5X,'***************************************************'
     1         ,'**********'
     2      /5X,'*****                                              '
     3         ,'     *****'
     4      /5X,'*****             PROGRAM RUTHSCAT ver.1.0         '
     5         ,'     *****'
     6      /5X,'*****                                              '
     7         ,'     *****'
     8      /5X,'***************************************************'
     9         ,'**********')
C
      WRITE(6,601) COMMENT
  601 FORMAT(/5X,'USER''S COMMENT:'/10X,A50)
C
      WRITE(6,602) Z1,Z2,T,THMIN,THMAX,DTH
  602 FORMAT(/3X,'-- GENERAL INPUTS --'
     1       /8X,'Z1    Z2    T        :',3F10.3
C170605
     2       /8X,'THMIN THMAX DTH      :',3F10.3)
C170605
C
C
C-------------------------
C---  PARAMETER CHECK  ---
C-------------------------
C
      ISTOP=0
      IF(Z1.LE.0.0D0) THEN
        WRITE(6,901)
        ISTOP=ISTOP+1
       ELSE IF(THMIN.LT.0.0D0) THEN
        WRITE(6,902)
        ISTOP=ISTOP+1
       ELSE IF(THMAX.LE.0.0D0) THEN
        WRITE(6,903)
        ISTOP=ISTOP+1
       ELSE IF(DTH.LE.0.0D0) THEN
        WRITE(6,904)
        ISTOP=ISTOP+1
       ELSE IF(THMIN.GE.THMAX) THEN
        WRITE(6,905)
        ISTOP=ISTOP+1
      END IF
C
      IF(ISTOP.NE.0) STOP 1
  901 FORMAT(/3X,'-- ERROR IN INPUT --:',/8X,'Z1 <= 0')
  902 FORMAT(/3X,'-- ERROR IN INPUT --:',/8X,'THMIN < 0')
  903 FORMAT(/3X,'-- ERROR IN INPUT --:',/8X,'THMAX <= 0')
  904 FORMAT(/3X,'-- ERROR IN INPUT --:',/8X,'DTH <= 0')
  905 FORMAT(/3X,'-- ERROR IN INPUT --:',/8X,'THMIN >= THMAX')
C
      RETURN
C
      END
C
C
C
C-----------------------------------------------------------------------
C---  Subroutines below are made by Y. Iseri (Chiba-Keizai College)  ---
C-----------------------------------------------------------------------
CCC
CCC        ***** FOPEN *****
CCC
      SUBROUTINE FOPEN
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER FNAME*50,STA*8,COMM*60,OFF*1
C ==========================================
C --- UNIT-5 ---
C     WRITE(*,100)
C 100 FORMAT(//1x,'>>> Input File-Name of Unit 5 ==> ',$)
C     READ(*,'(A)') FNAME

C     OPEN(5,FILE=FNAME,STATUS='OLD')

C --- OTHER UNITS ---
      READ(5,'(A)') COMM
      WRITE(*,200) COMM
  200 FORMAT(/1x,'---- (Comment in File) ----'
     1       /1x,' << ',A,' >>'
     2      //1x,'---- (  Open  Files  ) ----'
     3      //1x,'  Unit  Status      File-Name')

    1 READ(5,300) OFF,IUNIT,STA,FNAME
  300 FORMAT(A1,I3,1X,A8,2X,A50)
       IF(OFF.NE.' ') GO TO 1
       IF(IUNIT.GE.100 .OR. IUNIT.LT.0) GO TO 900
       WRITE(*,310) IUNIT,STA,FNAME
  310  FORMAT(1x,3X,I2,3X,A,3X,A)
       OPEN(IUNIT,FILE=FNAME,STATUS=STA)
      GO TO 1

C --- ONE MORE LINE FOR DISCRIMINATION ---
  900 READ(5,'(A)') COMM
      RETURN
      END
