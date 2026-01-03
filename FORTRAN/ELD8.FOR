C              .__________________________________________.
C              |           IN THE NAME OF GOD             |
C              |             M.S.Abrishamian              |
C              |   K.N. Toosi University of Technology    |
C              |               Tehran - Iran              |
C              |             P.O.Box 16315-1355           |
C              |      Email: msabrish@eetd.kntu.ac.ir     |
C              |__________________________________________|
c         **********************************************************
c         *   INTERNAL EM FIELD OF DIELECTRIC ELLIPTIC CYLINDER    *
c         *   METHOD = MOM                                         *
c         *   POLARIZATION = TM                                    *
c         *                   2000                                 *
c         **********************************************************
C---->eld8.for = TM &  DIELECTRIC & INTERNAL & NUMERIC (MOM)
	PROGRAM CYINMOM
      IMPLICIT NONE
	INTEGER MAX
	PARAMETER(MAX = 5500)
	REAL X(MAX),Y(MAX),D,LAMBDA,FR,AN,A,B,PHI0
	REAL BETA,PI,ER,SIGMA,W,E0,U,ET,PHI,BR,RMN
	COMPLEX G(MAX,MAX+1),JC,H,EN(MAX),SUM,C1,C2,C3,ERC
	INTEGER NCELL,M,N,IP,I
C----------------------<< CONSTANTS >>----------------------------------
	JC = (0.0,1.0)						 ! SQRT(-1.0)
	PI = 4.0*ATAN(1.0)					 ! PI=3.1415
	E0 = 8.85D-12						 ! FREE SPACE PERMITTIVITY
	FR = 1.D9							 ! FREQUENCY [Hz]
	W = 2*PI*FR							 ! ANGULAR FREQUENCY
	LAMBDA = 3.E8/FR					 ! WAVELENGTH
	BETA = 2.0*PI/LAMBDA				 ! PHASE CONSTANT
	ER = 5.0						     ! RELATIVE PERMITTIVITY
	SIGMA = 0.10D0  					 ! CONDUCTIVITY OF DIELECTRIC
	ERC = ER - JC*SIGMA/(W*E0)			 ! COMPLEX RELATIVE PERMITTIVITY	
	A = LAMBDA/2	                     ! LONG DIAMETER
	B = A/2					     ! SHORT DIAMETER
	PHI0 = 45.				        	 ! DETERMINE ANGLE OF INCIDENT
	WRITE(*,'(3X,F8.2)') A/LAMBDA
	WRITE(*,'(3X,F8.2)') B/LAMBDA
	PHI0 = PHI0*PI/180.0
C-----------------------------------------------------------------------
	NCELL = 0
	AN = 0.025*LAMBDA/CSQRT(ERC)
	D = AN*SQRT(PI)
	CALL ELLIPTIC(A,B,D,X,Y,NCELL)
	PRINT*,' NCELL=',NCELL
c	PAUSE
C-----------------------------------------------------------------------
	U = BETA*AN
	C1 = 1. + 0.5*(ERC-1.)*(JC*PI*U*H(1,U)+2.)
	C2 = JC*PI*U*0.5*(ERC-1)*REAL(H(1,U))
	DO M = 1 , NCELL
	 G(M,NCELL+1)=CEXP(JC*BETA*(X(M)*COS(PHI0)+Y(M)*SIN(PHI0)))
	 DO N = 1 , NCELL
	  IF(M.EQ.N)THEN
	   G(N,N) = C1
	  ELSE
	   RMN = SQRT( (X(M)-X(N))**2 + (Y(M)-Y(N))**2 )
	   BR = BETA*RMN
	   G(M,N) = C2*H(0,BR)
	  ENDIF
	 ENDDO
	ENDDO
C-----------------------------------------------------------------------
	PRINT*,' MATRIX COMPUTATION IS IN PROGRESS'
       CALL CHOLESKY(G,MAX,MAX+1,NCELL,NCELL+1,EN)
	PRINT*,' MATRIX COMPUTATION IS COMPLETED'
C-----------------------------------------------------------------------
	OPEN(1,FILE='s2.p')
	DO IP = 0 , 360  , 5
	 PHI = PI*IP/180
	 SUM = (0.0,0.0)
	 DO I = 1 , NCELL
	  C3 = EXP(-JC*BETA*(X(I)*COS(PHI)+Y(I)*SIN(PHI)))
	  SUM = SUM + (ERC-1)*EN(I)*U*REAL(H(1,U))*C3
	 ENDDO
	 WRITE(*,41)REAL(IP),PI*CABS(SUM)**2/2
	 WRITE(1,41)REAL(IP),PI*CABS(SUM)**2/2
	ENDDO
	CLOSE(1)
C-----------------------------------------------------------------------
	OPEN(2,FILE='sx.p')
	OPEN(3,FILE='sy.p')
	DO I = 1 , NCELL
	 ET = CABS(EN(I))
	 IF(Y(I).EQ.0.0)THEN
	  WRITE(2,41) X(I)*100 , ET
	  PRINT41, X(I)*100 , ET
	 ENDIF
	 IF(X(I).EQ.0.0)THEN
	  WRITE(3,41) Y(I)*100 , ET
	  WRITE(*,41) Y(I)*100 , ET
	 ENDIF
	ENDDO
	CLOSE(2)
	CLOSE(3)
 41   FORMAT(1X,2F15.4)
	STOP
	END
C***********************************************************************
	SUBROUTINE ELLIPTIC(A,B,D,X,Y,L)
	IMPLICIT NONE
	REAL A,B,X(*),Y(*),U,V,D
	INTEGER MAX1,MAX2,L,IX,IY
	MAX1 = A/D+1
	MAX2 = B/D+1
	DO IX = -MAX1 , MAX1
	 DO IY =  -MAX2 , MAX2
	  U = D*IX
	  V = D*IY	
	  IF( (U/A)**2+(V/B)**2.LE.1.0)THEN
	   L = L + 1
	   X(L) = U
	   Y(L) = V
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
