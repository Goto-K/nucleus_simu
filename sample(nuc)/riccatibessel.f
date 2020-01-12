C=======================================================================
      SUBROUTINE RICCATIBESSEL(L,Z,FJ,FN,FJP,FNP)
      IMPLICIT REAL*8 (A-H,O-Z)
C=======================================================================
C
C-----L  : orbital angular momentum (0-3)
C-----Z  : K * R
C-----
C-----FJ : z*j_L(z)
C-----FN : z*n_L(z)
C-----FJP: d[z*j_L(z)]/dz (= 1/K * d[z*j_L(z)]/dr)
C-----FNP: d[z*n_L(z)]/dz (= 1/K * d[n*j_L(z)]/dr)
C-----
C-----So, if you need d[z*j_L(z)]/dr, put FJP*K in your code.
C-----You can use the following expression.
C-----
C-----          TEXT                          CODE
C-----
C-----   \bar{h}^{(+)}(KR)       -->     DCMPLX(-FN, FJ)
C-----   \bar{h}^{(-)}(KR)       -->     DCMPLX(-FN,-FJ)
C-----
C-----   [\bar{h}^{(+)}(KR)]'    -->     DCMPLX(-FNP, FJP)*FKAY
C-----   [\bar{h}^{(-)}(KR)]'    -->     DCMPLX(-FNP,-FJP)*FKAY
C-----
C
      COMMON /MACHINE/EPSLIM,EXPMAX
CCC---------------------------------------------------------------------
C
      IF(Z.LE.EPSLIM) Z=EPSLIM
      IF(L.EQ.0) THEN
        FJ=SIN(Z)
        FN=-COS(Z)
        FJP=COS(Z)
        FNP=SIN(Z)
       ELSE IF(L.EQ.1) THEN
        FJ=(SIN(Z)-Z*COS(Z))/Z
        FN=-(COS(Z)+Z*SIN(Z))/Z
        FJP=(Z*COS(Z)-SIN(Z))/Z**2+SIN(Z)
        FNP=(COS(Z)-Z**2*COS(Z)+Z*SIN(Z))/Z**2
       ELSE IF(L.EQ.2) THEN
        FJ=((3.0D0-Z**2)*SIN(Z)-3.0D0*Z*COS(Z))/Z**2
        FN=-((3.0D0-Z**2)*COS(Z)+3.0D0*Z*SIN(Z))/Z**2
        FJP=(-Z*(-6.0D0+Z**2)*COS(Z)+3.0D0*(-2.0D0+Z**2)*SIN(Z))/Z**3
        FNP=(-Z*(-6.0D0+Z**2)*SIN(Z)-3.0D0*(-2.0D0+Z**2)*COS(Z))/Z**3
       ELSE IF(L.EQ.3) THEN
        FJ=((15.0D0-6.0D0*Z**2)*SIN(Z)-Z*(15.0D0-Z**2)*COS(Z))/Z**3
        FN=-((15.0D0-6.0D0*Z**2)*COS(Z)+Z*(15.0D0-Z**2)*SIN(Z))/Z**3
        FJP=-(3.0D0*Z*(-15.0D0+2.0D0*Z**2)*COS(Z)
     1        +(45.0D0-21.0D0*Z**2+Z**4)*SIN(Z))/Z**4
        FNP=((45.0D0-21.0D0*Z**2+Z**4)*COS(Z)
     1        +3.0D0*Z*(15.0D0-2.0D0*Z**2)*SIN(Z))/Z**4
       ELSE
        STOP 1
      END IF
C
      RETURN
      END
