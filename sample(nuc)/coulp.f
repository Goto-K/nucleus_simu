      SUBROUTINE COULP(RHO, ETA, L, FF, FFP, GG, GGP, CPH, NW)
C----- Coulomb wave functions for E > 0
C-----   G --> COS(RHO-ETA*Log(2*rho)-L*Pi/2+C@H(L+1))
C-----   F --> SIN(RHO-ETA*Log(2*rho)-L*Pi/2+C@H(L+1))
C-----
C----- RHO: K*R (>=0)
C----- ETA: Sommerfeld parameter
C-----
C----- With a given L (in the argument of COULP), this subroutine
C----- calculates the following values
C-----
C----- G(RHO,ETA;L) ---> GG(L+1)
C----- F(RHO,ETA;L) ---> FF(L+1)
C----- dG(RHO,ETA;L)/dRHO ---> GGP(L+1)
C----- dF(RHO,ETA;L)/dRHO ---> FFP(L+1)
C----- Coulomb phase shift ---> CPH(L+1)
C-----
C----- from 1 to L+1 for the argument of each array. For example,
C----- if you put L=2 in COULP, you will have GG(1), GG(2), and
C----- GG(3), for the irregular Coulomb wave function, corresponding
C----- to L=0, 1, and 2, respectively.
C-----
C-----NOTICE1:
C----- If calculation fails, or a result is larger than 10^70,
C----- you will have NW(L+1)=1. Do not forget to check NW after
C----- you call COULP.
C-----
C-----NOTICE2:
C----- The maximum value of L in the argument of COULP is 10.
C-----
      IMPLICIT REAL*8(A-H, O-Z)
      DIMENSION FF(12),FFP(12),GG(12),GGP(12),CPH(12),NW(12)
      DIMENSION G(18),G1(18),F(18),F1(18)
      SAVE AG, AF,AG1,AF1, G, G1
      DATA AG,AF, AG1, AF1,G, G1
     C/ 0.1223404015123203D+01,0.7063326374590501D+00
     C, 0.7078817733790302D+00,0.4086957324148128D+00
     C,1.0000000000000000D0,0.0495957016858696D0,-0.0088888888888889D0
     C,0.0024551991808936D0,-0.0009108958061339D0,0.0008453619999332D0
     C,-0.0004096926350538D0,0.0007116506205636D0,-0.0004147346525854D0
     C,0.0010949251404694D0,-0.0007378802572622D0,0.0026667876712056D0
     C,-0.0020270233593062D0,0.0094246529246488D0,-0.0079377765835908D0
     C,0.0455657466669352D0,-0.0419619350508973D0,0.2887892996933970D0
     C,-1.0000000000000000D0,0.1728260369359927D0,-0.0003174603174603D0
     C,0.0035812148490633D0,-0.0003117824679729D0,0.0009073966425577D0
     C,-0.0002128570748571D0,0.0006215584099796D0,-0.0002579671336123D0
     C,0.0008175755607792D0,-0.0005060050773457D0,0.0017543295499921D0
     C,-0.0014772097102188D0,0.0055744199870004D0,-0.0060311585140914D0
     C,0.0245908268807787D0,-0.0328733641949943D0,0.1437985328996377D0/
      DO 500 I=1,18
      F(I)=-DABS(G(I))
  500 F1(I)=DABS(G1(I))
      F(1)=1.D0
      L1=L+1
      DO 810 I=1,L1
  810 NW(I)=0
C
      IF(RHO.LE.0.D0) GO TO 801
      ETA2=ETA+ETA
      E2=ETA*ETA
      RHO2=RHO+RHO
      R2=RHO*RHO
      ETARHO=ETA*RHO
      ER2 =ETARHO+ETARHO
C*** COULOMB PHASE SHIFT
      L0=L
      IF(L.LE.10) L0=10
      EL=L0+1
      A=DATAN(ETA/EL)
      B=DSQRT(EL*EL+E2)
      B2=B*B
      B3=B*B2
      B5=B3*B2
      S=(EL-0.5D0)*A+ETA*(DLOG(B)-1.D0)-DSIN(A)/(12.D0*B)
     C+DSIN(3.D0*A)/(360.D0*B3)-DSIN(5.D0*A)/(1260.D0*B5)
      IF(L.GE.50) GO TO 101
      B7=B5*B2
      S=S+DSIN(7.D0*A)/(1680.D0*B7)
      IF(L.GE.40) GO TO 101
      B9=B7*B2
      S=S-DSIN(9.D0*A)/(1188.D0*B9)
      IF(L.GE.20) GO TO 101
      B11=B9*B2
      B13=B11*B2
      S=S+DSIN(11.D0*A)*691.D0/(360360.D0*B11)-DSIN(13.D0*A)/(156.D0*B13
     C)
      IF(L0.EQ.L) GO TO 101
      K=L0
  100 S=S-DATAN(ETA/K)
      K=K-1
      IF(K.EQ.L) GO TO 101
      GO TO 100
  101 K=L+1
  102 CPH(K)=S
      K=K-1
      IF(K.EQ.0) GO TO 103
      S=S-DATAN(ETA/K)
      GO TO 102
  103 CONTINUE
C*** COULOMB WAVE FUNCTION
      RH=RHO
      RH2=RHO2
      D=4.D0*(R2-ER2)+1.D0
      LC=(-1.D0+DSQRT(DABS(D)))/2.D0
      IF(D.LT.1.D0) LC=-1
      RC=E2/10.D0+1.5D0*ABS(ETA)+16.D0
      IF(ETA.EQ.0.D0) GO TO 4
      IF(ETA.LT.6.D0) GO TO 10
      IF(RHO.GE.RC) GO TO 4
      IF(RHO.LT.ETA2) GO TO 5
      IF((RHO-ETA2)/(RC-RHO).LE.1.D0) GO TO 5
   20 RH=RC
      RH2=RC+RC
      GO TO 4
   10 IF(RHO.GE.RC) GO TO 4
      GO TO 20
C*** ASYMPTOTIC EXPANSION WITH LARGE RHO
    4 CONTINUE
      TH=RH-ETA*DLOG(RH2)+CPH(1)
      SS1=1.D0
      SL1=0.D0
      TS1=0.D0
      TL1=1.D0-ETA/RH
      WSS=SS1
      WTS=TS1
      WSL=SL1
      WTL=TL1
      IF(ETA.EQ.0.D0) GO TO 401
      N=0
  400 C=N+1
      D=(C+C)*RH
      A=(N+C)*ETA/D
      B=(E2-N*C)/D
      SS=A*SS1-B*TS1
      TS=A*TS1+B*SS1
      SL=A*SL1-B*TL1-SS/RH
      TL=A*TL1+B*SL1-TS/RH
      SS1=SS
      TS1=TS
      SL1=SL
      TL1=TL
      WSS=WSS+SS
      WTS=WTS+TS
      WSL=WSL+SL
      WTL=WTL+TL
      N=N+1
      IF(DABS(SS/WSS).LE.1.D-15) GO TO 401
      IF(N.GE.200) GO TO 401
      GO TO 400
  401 CT=DCOS(TH)
      ST=DSIN(TH)
      F0=WTS*CT+WSS*ST
      G0=WSS*CT-WTS*ST
      F01=WTL*CT+WSL*ST
      G01=WSL*CT-WTL*ST
      IF(RHO.GE.RC) GO TO 7
      IF(ETA.EQ.0.D0) GO TO 7
      GO TO 6
C*** EXPANSION WITH 2*ETA=RHO
    5 CONTINUE
      RH=ETA2
      E23=ETA**(2.D0/3.D0)
      E43=E23*E23
      G0=0.D0
      F0=0.D0
      G01=0.D0
      F01=0.D0
      IF(ETA.GE.10.D0) GO TO 501
      G0=(((G(18)/E43+G(17))/E23+G(16))/E43+G(15))/E23
      F0=(((F(18)/E43+F(17))/E23+F(16))/E43+F(15))/E23
      G01=(((G1(18)/E23+G1(17))/E43+G1(16))/E23+G1(15))/E43
      F01=(((F1(18)/E23+F1(17))/E43+F1(16))/E23+F1(15))/E43
  501 IF(ETA.GE.20.D0) GO TO 502
      G0=((((G0+G(14))/E43+G(13))/E23+G(12))/E43+G(11))/E23
      F0=((((F0+F(14))/E43+F(13))/E23+F(12))/E43+F(11))/E23
      G01=((((G01+G1(14))/E23+G1(13))/E43+G1(12))/E23+G1(11))/E43
      F01=((((F01+F1(14))/E23+F1(13))/E43+F1(12))/E23+F1(11))/E43
  502 IF(ETA.GE.60.D0) GO TO 503
      G0=(((G0+G(10))/E43+G(9))/E23+G(8))/E43
      F0=(((F0+F(10))/E43+F(9))/E23+F(8))/E43
      G01=(((G01+G1(10))/E23+G1(9))/E43+G1(8))/E23
      F01=(((F01+F1(10))/E23+F1(9))/E43+F1(8))/E23
  503 IF(ETA.GE.300.D0) GO TO 504
      G0=((G0+G(7))/E23+G(6))/E43
      F0=((F0+F(7))/E23+F(6))/E43
      G01=((G01+G1(7))/E43+G1(6))/E23
      F01=((F01+F1(7))/E43+F1(6))/E23
  504 G0=((((G0+G(5))/E23+G(4))/E43+G(3))/E23+G(2))/E43+G(1)
      F0=((((F0+F(5))/E23+F(4))/E43+F(3))/E23+F(2))/E43+F(1)
      G01=((((G01+G1(5))/E43+G1(4))/E23+G1(3))/E43+G1(2))/E23+G1(1)
      F01=((((F01+F1(5))/E43+F1(4))/E23+F1(3))/E43+F1(2))/E23+F1(1)
      E6=DSQRT(DSQRT(E23))
      G0=AG*E6*G0
      F0=AF*E6*F0
      G01=AG1/E6*G01
      F01=AF1/E6*F01
      IF(DABS(1.D0-ETA2/RHO).LE.1.D-14) GO TO 7
C*** SOLVING THE DIFFERENTIAL EQUATION
    6 CONTINUE
      MG=1
      IF((LC.LT.0).OR.(L.GT.LC)) MG=0
  601 CONTINUE
      EA=6.D0
      EB=0.75D0
      A=1.D0-ETA2/RH
      DT=DMIN1(RH*EB,DABS(RHO-RH),(RH*EA)**(1.D0/3.D0))
      IF(DABS(A).GE.1.D-10)
     C DT=DMIN1(DT,DABS(EA/A)**(1.D0/3.D0),DABS(EA*RH/A)**(1.D0/3.D0))
      IF(RHO.LT.RH) DT=-DT
      DT2=DT*DT
      A=(1.D0-ETA2/RH)*DT2
      B=DT/RH
      C=B*DT2
      V2=G0
      V1=G01*DT
      V=-A/2.D0*G0
      G0=V2+V1+V
      G01=V1+V+V
      I1=1
      N=1
      IF(MG.EQ.0) GO TO 610
      U2=F0
      U1=F01*DT
      U=-A/2.D0*F0
      F0=U2+U1+U
      F01=U1+U+U
  610 D=N+2
      V0=V
      V=-(B*N*V+(A*V1+C*V2)/(N+1))/D
      V2=V1
      V1=V0
      G0=G0+V
      W=D*V
      G01=G01+W
      S=DMAX1(1.D0,DABS(G0))
      IF(DABS(G01).GE.1.D+70) GO TO 805
      IF(MG.EQ.0) GO TO 612
      U0=U
      U=-(B*N*U+(A*U1+C*U2)/(N+1))/D
      U2=U1
      U1=U0
      F0=F0+U
      F01=F01+D*U
  612 IF(DABS(W)/S.LE.1.D-15) GO TO 611
      N=N+1
      IF(N.GE.1000) GO TO 805
      GO TO 610
  611 G01=G01/DT
      IF(MG.EQ.1) F01=F01/DT
      RH=RH+DT
      IF(DABS(RHO-RH).LE.1.D-14) GO TO 7
      GO TO 601
C*** ZENKA-SHIKI
    7 CONTINUE
      I1=0
      EW=1.D-60
      IF((LC. LT.-1).OR.(L.GT.LC)) GO TO 710
      GG(1)=G0
      GGP(1)=G01
      FF(1)=F0
      FFP(1)=F01
      K=1
  701 IF(K.GT.L) RETURN
      K1=K+1
      S=ETA/K
      A=K/RHO+S
      B=DSQRT(1.D0+S*S)
      GG(K1)=(A*GG(K)-GGP(K))/B
      GGP(K1)=B*GG(K)-A*GG(K1)
      FF(K1)=(A*FF(K)-FFP(K))/B
      FFP(K1)=B*FF(K)-A*FF(K1)
      K=K+1
      GO TO 701
  710 GG(1)=G0
      GGP(1)=G01
      K=1
      K1=1
  711 IF(K.GT.L+1) GO TO 720
      S=ETA/K
      A=K/RHO+S
      EP3=DMIN1(1.D+70,1.D+75/(A+1.D0))
      IF(DABS(GGP(K1)).GE.EP3) GO TO 715
      B=DSQRT(1.D0+S*S)
      K1=K+1
      GG(K1)=(A*GG(K)-GGP(K))/B
      GGP(K1)=B*GG(K)-A*GG(K1)
      K=K+1
      GO TO 711
  715 I1=K1
      EW=1.D-70
  720 U=GG(K1)*EW
      UP=(GG(K1)+GGP(K1))*EW
      V=U
      VP=GGP(K1)*EW
      C=V
      D=VP
  721 S=ETA/K
      A=K/RHO+S
      B=DSQRT(1.D0+S*S)
      U1=U
      U=(A*U-UP)/B
      UP=B*U1-A*U
      V1=V
      V=(A*V-VP)/B
      VP=B*V1-A*V
      IF(DABS(V/C).GE.1.D+9) GO TO 730
      K=K+1
      IF(K.GE.L+1000) GO TO 801
      GO TO 721
  730 B=U/V
      FFP(K1)=(1.D0+(1.D0-B)*D/C)/C*EW
      FF(K1)=(1.D0-B)/C*EW
      K=K1
  740 K1=K-1
      IF(K.LE.1) GO TO 745
      S=ETA/K1
      A=K1/RHO+S
      B=DSQRT(1.D0+S*S)
      FF(K1)=(A*FF(K)+FFP(K))/B
      FFP(K1)=A*FF(K1)-B*FF(K)
      K=K-1
      GO TO 740
  745 IF(I1.EQ.0) RETURN
      GO TO 805
C*** SPECIAL CASES OR UNSOLVABLE
  800 DO 802 I=1,L1
      FF(I)=0.D0
      FFP(I)=0.D0
      GG(I)=0.D0
  802 GGP(I)=0.D0
      RETURN
  801 I1=1
  805 DO 803 I=I1,L1
      NW(I)=1
      FF(I)=0.D0
      FFP(I)=0.D0
      GG(I)=1.D+70
  803 GGP(I)=1.D+70
      IF(ETA+DABS(RHO).NE.0.D0) RETURN
      GG(1)=1.D0
      FFP(1)=1.D0
      NW(1)=0
      RETURN
      END
