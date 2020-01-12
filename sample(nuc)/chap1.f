      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(10)
      PARAMETER (AC0=931.494013D0)
      CHARACTER*3 ANCL(10)
CCC
      PI=DACOS(-1D0)
      HC=197.3269601D0
CCC
      AP   = 1.00782503207D0
      AN   = 1.00866491574D0
      AD   = 2.01410177785D0
      A4He = 4.00260325415D0
      A11C = 11.011433613D0
      A11B = 11.009305406D0
      A12C = 12.0D0
CCC
      A(1) = 1.00782503207D0
      A(2) = 1.00866491574D0
      A(3) = 2.01410177785D0
      A(4) = 4.00260325415D0
      A(5) = 12.0D0
CCC
      ANCL(1)='P'
      ANCL(2)='N'
      ANCL(3)='D'
      ANCL(4)='4He'
      ANCL(5)='12C'
CCC
      WRITE(*,*) 'Answer 1.1'
      WRITE(*,*) 
      DO I=1,5
         WRITE(*,"(2X,A3,A7,1X,1F15.8,A6)") 
     1        ANCL(I),' =',A(I)*AC0,' [MeV]' 
      END DO
CCC
      WRITE(*,*) 
      WRITE(*,*) 'Answer 1.2'
      WRITE(*,*) 
      WRITE(*,1100) '\rho =',1d0/(4D0*PI*1.1**3/3D0),' [fm^-3]'
CCC
      WRITE(*,*) 
      WRITE(*,*) 'Answer 1.3'
      WRITE(*,*) 
      WRITE(*,1000) "E(n+p)=",(AN+AP-AD)*AC0,' [MeV]'
CCC
      WRITE(*,*) 
      WRITE(*,*) 'Answer 1.4'
      WRITE(*,*) 
      WRITE(*,1000) "E(4Hex3)= ", (A4He*3D0-A12C)*AC0,' [MeV]'
CCC
      WRITE(*,1000) "E(11C+n)= ", (A11C+AN-A12C)*AC0,' [MeV]'
CCC
C      WRITE(*,"(A10,1F15.8)") "B_10B+p= ", (A11B+AP-A12C)*AC0
CCC
      WRITE(*,*) 
      WRITE(*,*) 'Answer 1.5'
      WRITE(*,*) 
      WRITE(*,1000) "E/A= ", (AN*6D0+AP*6D0-A12C)*AC0/12D0,' [MeV]'
CCC
      WRITE(*,*) 
      WRITE(*,*) 'Answer 1.7'
      WRITE(*,*) 
CCC
      FKF=(3D0*PI**2/2D0*0.17D0)**(1./3.)
C      WRITE(*,*) FKF
      WRITE(*,1000) 'E_F=',FKF**2*HC**2/2D0/A(1)/AC0,' [MeV]'
      WRITE(*,1000) 'P_F=',FKF*HC,' [MeV/c]'
      WRITE(*,1000) 'k_F=',FKF,' [1/fm]'
CCC
      WRITE(*,*) 
      WRITE(*,*) 'Answer 1.8'
      WRITE(*,*) 
CCC
      PIMASS=0.5D0*HC
      WRITE(*,1000) 'm_pi=',PIMASS,' [MeV]'
CCC
 1000 FORMAT(2X,A10,1F10.5,A8)
 1100 FORMAT(2X,A10,1F10.5,A8)
      STOP
      END
