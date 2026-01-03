	PROGRAM CYINMOM
      IMPLICIT NONE
	INTEGER MAX
	PARAMETER(MAX = 2500)
	REAL X(MAX),Y(MAX),D,LAMBDA,AN,ROH,PHI,A,B,ET,PHI0
	REAL PI,U0,E0,C,FR,BETA,ETA,W,ER,SIGMA,U,BR,BJ1
	COMPLEX G(2*MAX,2*MAX+1),JC,H,EN(2*MAX),ERC,UKP,EPH,C3,SUM
	INTEGER M,N,I,J,IX,JY,NCELL
C**************************************************************************
	print*,'                                                   '
      print*,'     ษออออออออออออออออออออออออออออออออออออออออออป  '
      print*,'     บ            IN THE NAME OF GOD            บ  '
      print*,'     บ             M.S.Abrishamian              บ  '
      print*,'     บ   K.N. Toosi University of Technology    บ  '
      print*,'     บ               Tehran - Iran              บ  '
      print*,'     บ             P.O.Box 16315-1355           บ  '
      print*,'     บ      Email: msabrish@eetd.kntu.ac.ir     บ  '
      print*,'     ศออออออออออออออออออออออออออออออออออออออออออผ  '
	print*,'                                                   '
C**************************************************************************
C---->cd11.for = TE &  DIELECTRIC & EXTERNAL & NUMERIC (MOM)
C           TE-WAVE SCATTERING BY DIELECTRIC CYLINDER OF ARBITRARY 
C           CROSS-SECTION SHAPE 
C           BY: JACK H. RICHMOND, IEEE TRANSACTION ON ANTENNA AND 
C               PROPAGATION, VOL. AP-14,NO.4 JULY, 1966
C----------------------<< CONSTANTS >>-------------------------------------
	JC = (0.D0,1.D0)               !  SQRT(-1)
	PI = 4.0 * DATAN(1.D0)         !  PI = 3.1415
	C = 3.E8                       !  SPEED OF LIGHT
	E0 = 1.E-9/(36.*PI)            !  E0=PERMITTIVITY OF FREE SPACE
	U0 = 4.*PI*1.E-7               !  E0 =PERMEABLITY OF FEEE SPACE
	FR = 1.E9                      !  FREQUENCY
	W = 2*PI*FR					   !  ANGULAR  FREQUENCY
	LAMBDA = C/FR                  !  WAVELENGTH
	BETA = 2.D0 * PI / LAMBDA      !  BETA = WAVE NUMBER
	ETA = SQRT(U0/E0)              !  FREE SPACE IMPEDANCE = 377[OHMS]
	SIGMA = 0.1D0                 !  SIGMA = CONDUCTIVITY
	ER = 5.D0					   !  RLARTIVE PERMITTIVITY 
	ERC = ER - JC*SIGMA/(W*E0)	   !  COMPLEX PERMITTIVITY 
	 
	A = LAMBDA/4			       !  LONG DIAMETER
	B = LAMBDA/8				   !  SHORT DIAMETER

	AN = 0.02*LAMBDA/CSQRT(ERC)   !  RADIUS OF EACH CELL
	U = BETA*AN					   !  ELETRICAL CONSTANT
	D = AN*SQRT(PI)				   !  SIDE OF EACH SQUARE CELL
	PHI0 = PI/4
	CALL ELIPTIC(A,B,D,X,Y,NCELL)	   !  MAKE SQUARE CELL INSIDE OF CYLINDER 
	PRINT*,' NUMBER OF CELLS =',NCELL		 
C----------------------<< FIELD GENERATION IN EACH CELL>>---------------
	DO I = 1 , NCELL
	  C3 = EXP(JC*BETA*(X(I)*COS(PHI0)+Y(I)*SIN(PHI0))) 
	 G(I,2*NCELL+1) = ETA*C3*SIN(PHI0) 
	 G(I+NCELL,2*NCELL+1)= -ETA*C3*COS(PHI0)	 ! INCIDENT WAVE
	ENDDO
      DO I = 1 , 2
       DO J = 1 , 2
	  IX = (I-1)*NCELL
	  DO M = 1 , NCELL	 
	   IX = IX + 1
	   JY = (J-1) * NCELL
	   DO N = 1 , NCELL	    
	    JY = JY + 1
	    ROH = SQRT( (X(M)-X(N))**2 + (Y(M)-Y(N))**2 )
	    BR = BETA*ROH
	    IF(I.EQ.1.AND.J.EQ.1)THEN
	     IF(M.EQ.N)THEN
	      G(IX,JY) = 1. + (ERC-1.)*( 0.25*JC*PI*U*H(1,U) +1.)
	     ELSE
	      BJ1 = REAL(H(1,U))
	      UKP = JC*PI*AN*BJ1*(ERC-1.)/(2*ROH**3)
	      G(IX,JY) = UKP*( BR*(Y(M)-Y(N))**2*H(0,BR) +    	 ! AMN 
     &       	        ((X(M)-X(N))**2-(Y(M)-Y(N))**2)*H(1,BR) )
	     ENDIF
	    ELSEIF(I.EQ.2.AND.J.EQ.2)THEN
	     IF(M.EQ.N)THEN
	      G(IX,JY) = 1. + (ERC-1.)*( 0.25*JC*PI*U*H(1,U) +1.)
	     ELSE
	      BJ1 = REAL(H(1,U))
	      UKP = JC*PI*AN*BJ1*(ERC-1.)/(2*ROH**3)
	      G(IX,JY) = UKP*( BR*(X(M)-X(N))**2*H(0,BR) +    	 ! DMN 
     &       	         ( (Y(M)-Y(N))**2-(X(M)-X(N))**2 )*H(1,BR) )
	     ENDIF
	    ELSE
	     IF(M.EQ.N)THEN	      
	      G(IX,JY) = (0.0,0.0)
	     ELSE
	      BJ1 = REAL(H(1,U))
	      UKP = JC*PI*AN*BJ1*(ERC-1.)/(2*ROH**3)
	      G(IX,JY) = UKP*(X(M)-X(N))*(Y(M)-Y(N))*(2*H(1,BR)- ! BMN & CMN
     &                 BR*H(0,BR))
	     ENDIF
	    ENDIF	    
	   ENDDO
	  ENDDO
       ENDDO
      ENDDO
C-----------------------------------------------------------------------
      PRINT*,' MATRIX DIMENSION IS=',IX,'      X',JY
	PRINT*,' MATRIX COMPUTATION IS IN PROGRESS'
c	CALL GAUSSC(2*MAX,2*MAX,2*MAX+1,2*NCELL,2*NCELL+1,G,EN)
       CALL CHOLESKY(G,2*MAX,2*MAX,2*NCELL,2*NCELL+1,EN)
C-----------------------------------------------------------------------
	OPEN(1,FILE='sx.p')
	OPEN(2,FILE='sy.p')
	DO I = 1 , NCELL
	 ET = SQRT( CABS(EN(I))**2 + CABS(EN(I+NCELL))**2 )
	 IF(Y(I).EQ.0.0)THEN
	  WRITE(1,41) X(I)*100 , ET
	 ENDIF
	 IF(X(I).EQ.0.0)THEN
	  WRITE(2,41) Y(I)*100 , ET
	  WRITE(*,41) Y(I)*100 , ET	 
	 ENDIF
	ENDDO
	CLOSE(1)
	CLOSE(2)
C-----------------------------------------------------------------------
	OPEN(1,FILE='s2.p')

	BJ1 = REAL(H(1,U))
	DO J = 0 , 360 , 5
	 PHI = PI*J/180
	 SUM = (0.0,0.0)
	 DO I = 1 , NCELL
	  C3 = CEXP(JC*BETA*( X(I)*COS(PHI)+Y(I)*SIN(PHI)) )

	  EPH = EN(I)*SIN(PHI)- EN(I+NCELL)*COS(PHI)

	  SUM = SUM + (ERC-1.)*EPH*U*BJ1*C3
	 ENDDO
	 SUM = SUM/ETA
	 WRITE(*,41)REAL(J),0.5*PI*CABS(SUM)**2
	 WRITE(1,41)REAL(J),0.5*PI*CABS(SUM)**2
	ENDDO
	CLOSE(1)
 41   FORMAT(3X,F10.4,3X,F10.4)
      STOP
	END
C***********************************************************************
	SUBROUTINE ELIPTIC(A,B,D,X,Y,L)	
	IMPLICIT NONE
	REAL A,B,X(*),Y(*),RV,D,XA,YB
	INTEGER MA,MB,L,IX,IY
	MA = A/D + 3
	MB = B/D + 3
	DO IX = -MA , MA
	 DO IY =  -MB , MB
	  XA = IX*D
	  YB = IY*D
	  RV = XA*XA/A/A + YB*YB/B/B
	  IF(RV.LE.1)THEN
	   L = L + 1
	   X(L) = IX*D
	   Y(L) = IY*D
	  ENDIF
	 ENDDO
	ENDDO
	RETURN
	END
C***********************************************************************
	COMPLEX FUNCTION H(N,X)
	IMPLICIT  NONE
	REAL X,T,Y,SIG,BN,BK,BKP1,BKP2,TINY,JN,J1,J0
	REAL Y0,YN,YK,YKM1,YKM2,C(5),PI,GAMMA
	INTEGER M,MM,IU,N,K
	DATA C/17.1032, 0.2639, 0.6487, -0.0018, 0.6457/
	PI = ACOS(-1.D0)
	IU = ABS(X)
	M = INT( C(1)+C(2)*IU+C(3)*IU+C(4)*IU+C(5)*MAX(N,IU) )
	M = 4*INT(M/4)
	MM = M
	TINY = 3.D-35
	K = MM
	BKP1 = 0.0
	BK = TINY*M
	T = 2.D0*BK
	Y = 2.D0*BK/K
	SIG = 1.D0
	DO K = MM-1,0,-1
	 BKP2 = BKP1
	 BKP1 = BK
	 BK = 2*(K+1)*BKP1/X - BKP2 + TINY*M
	 IF(K.EQ.N) BN = BK
	 IF(K.GT.0.AND.MOD(K,2).EQ.0)THEN
	  T = T + 2*BK
	  SIG = -SIG
	  Y = Y + SIG*2.D0*BK/K
	 ENDIF
	 IF(K.EQ.0)T = T + BK
	ENDDO
	GAMMA = 0.57721566490153286D0
	Y = 2.D0/PI*( LOG(X/2.) + GAMMA )*BK - 4.D0/PI*Y
	JN = BN/T
	J0 = BK/T
	J1 = BKP1/T
	Y0 = Y/T
	IF(N.GT.0)THEN
	 YKM1 = Y0
	 YK = (J1*Y0 - 2.D0/(PI*X))/J0
	 DO K = 2 , N
	  YKM2 = YKM1
	  YKM1 = YK
	  YK = 2.D0*(K-1)*YKM1/X - YKM2
	 ENDDO
	 YN = YK
	ELSE
	 YN = Y0
	ENDIF
	H = CMPLX(JN ,-YN)
	RETURN
	END
C***********************************************************************
	SUBROUTINE CHOLESKY(A,MAX,MAX1,N,M,X)
	COMPLEX A(MAX,MAX1),X(MAX),SUM
	DO J = 2 , M
	 A(1,J) = A(1,J)/A(1,1)
	ENDDO
	DO I = 2 , N
	 DO II = I , N
	  SUM = (0.0,0.0)
	  DO K = 1 , I-1
	   SUM = SUM + A(II,K)*A(K,I)
	  ENDDO
	  A(II,I) = A(II,I) - SUM
	 ENDDO
	 DO J = I+1, M
	  SUM = (0.0,0.0)
	  DO K = 1 , I-1
	   SUM = SUM + A(I,K)*A(K,J)
	  ENDDO
	  A(I,J) = (A(I,J) - SUM)/A(I,I)
	 ENDDO
	ENDDO
	X(N) = A(N,M)
	DO NN = 1 , N-1
	 SUM = (0.0,0.0)
	 I = N - NN
	 DO J = I+1 , N
	  SUM = SUM + A(I,J)*X(J)
	 ENDDO
	 X(I) = A(I,M) - SUM
	ENDDO
	RETURN
	END
C***********************************************************************
	SUBROUTINE GAUSSC(IR,IC,IC1,N,M,A,X)
	IMPLICIT NONE
	INTEGER IR,IC,IC1,L,M,N,KP1,IP1,J,K,I,JU,NV
	COMPLEX TV,TU,SUM,A(IR,IC1),X(IC)
	REAL BIG,AB
	L = N-1
	DO K = 1 , L
	 JU = K
	 BIG = CABS(A(K,K))
	 KP1 = K + 1
	 DO I = KP1 , N
	  AB = CABS(A(I,K))
	  IF(ABS(BIG-AB).LT.0.D0)THEN
	   BIG = AB
	   JU = I
	  ENDIF
	 ENDDO
	 IF(JU-K.EQ.0)THEN
	  DO I = KP1 , N
	   TU = A(I,K)/A(K,K)
	   DO J = KP1 , M
	    A(I,J)= A(I,J) - TU*A(K,J)
	   ENDDO
	  ENDDO
	  DO I = KP1 , N
	   A(I,K) = (0.0,0.0)
	  ENDDO
	 ELSE
	  DO J = K , M
	   TV = A(JU,J)
	   A(JU,J) = A(K,J)
	   A(K,J) = TV
	  ENDDO
	 ENDIF
	ENDDO
	X(N) = A(N,M)/A(N,N)
	DO  NV = 1 , L
	 SUM = (0.0,0.0)
	 I = N - NV
	 IP1 = I + 1
	 DO J = IP1 , N
	  SUM = SUM + A(I,J)*X(J)
	 ENDDO
	 X(I) = (A(I,M)-SUM)/A(I,I)
	ENDDO
	RETURN
	END
C***********************************************************************
c         close 
c         clear all
c         load S.p;
c         load c.p
c         plot(S(:,1),S(:,2));grid;hold;
c         plot(c(:,1),c(:,2),'*');
c         xlabel(' angle "Phi" ')
c         ylabel('RCS/Lambda')
c         title('RCS of a dielectric cylinder: TE Mode ')
c         legend('Series', 'MoM ')
c           