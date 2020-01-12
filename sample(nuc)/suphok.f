CCC
CCC         ***** SUPHOK *****
CCC    THIS PROGRAM IS O.K. ONLY WHEN JISU=3.    ------ CARE ----
CCC
      FUNCTION SUPHOK(XIN,X,Y,NRMAX,NRDIM)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(NRDIM),Y(NRDIM)
CCC ===================================================
      JISU=3
      IZISU=JISU+1
C
      IF(NRMAX.GE.IZISU) GO TO 4000
      WRITE(6,400) NRMAX,JISU
  400 FORMAT(1H0,' NRMAX,JISU=',2I5,10X,'STOP IN SUPHOK')
      STOP
 4000 CONTINUE
C
      DO 1 I=1,NRMAX-1
      FUGO=(XIN-X(I+1))*(XIN-X(I))
      IF(FUGO.GT.0.0) GO TO 1
      JMIN=I-1
      IF(I.EQ.1) JMIN=1
      IF(I.EQ.NRMAX-1) JMIN=NRMAX-3
      GO TO 2
    1 CONTINUE
    2 CONTINUE
C
      X0=X(JMIN)
      X1=X(JMIN+1)
      X2=X(JMIN+2)
      X3=X(JMIN+3)
C
      Y0=Y(JMIN)
      Y1=Y(JMIN+1)
      Y2=Y(JMIN+2)
      Y3=Y(JMIN+3)
C
      XIN0=XIN-X0
      XIN1=XIN-X1
      XIN2=XIN-X2
      XIN3=XIN-X3
C
      X01=X0-X1
      X02=X0-X2
      X03=X0-X3
C
      X10=X1-X0
      X12=X1-X2
      X13=X1-X3
C
      X20=X2-X0
      X21=X2-X1
      X23=X2-X3
C
      X30=X3-X0
      X31=X3-X1
      X32=X3-X2
C
      XX0= XIN1/X01 *XIN2/X02 *XIN3/X03
      XX1= XIN0/X10 *XIN2/X12 *XIN3/X13
      XX2= XIN0/X20 *XIN1/X21 *XIN3/X23
      XX3= XIN0/X30 *XIN1/X31 *XIN2/X32
C
      SUPHOK= XX0*Y0 +XX1*Y1 +XX2*Y2 +XX3*Y3
C
      RETURN
      END
