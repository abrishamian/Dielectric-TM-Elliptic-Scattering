C              |           IN THE NAME OF GOD             |
C              |             M.S.Abrishamian              |
C              |   K.N. Toosi University of Technology    |
C              |               Tehran - Iran              |
C              |             P.O.Box 16315-1355           |
C              |      Email: msabrish@eetd.kntu.ac.ir     |
C              |__________________________________________|
c      ********************************************************************
c      *    SCATTERING AND PENETERATION OF EM WAVE BY A LOSSY DIELECTRIC 	*
C      *    ELLIPTIC CYLINDER                                             *
c      *    METHOD = ANALYTIC      (SERIES SOLUTION)                      *
c      *    POLARIZATION = TM                                             *
c      *                                                                  *
c      ********************************************************************
	PROGRAM ELLIPTIC_CYLINDER
	IMPLICIT NONE
	INTEGER IMAX
	PARAMETER(IMAX=50)
	COMPLEX*16 JC,BETA,Q2,Q1,CSI,G(2*IMAX+2,2*IMAX+3)
	COMPLEX*16 XA(2*IMAX+2),XB(2*IMAX+2),EZ,EZA,EZB 
	COMPLEX*16 F1R,D1R,F2R,D2R,CJ,CSF,CSD,RES
	REAL*8 PI,E0,FR,LAMBDA,BETA0,A,B,F,W,ER,SIGMA
	REAL*8 U1,PHI0,U,V,X1,Y1,R,PHI
	INTEGER M,N,I,J,IX,JY,KF,KC,MAX,ID,IU
c---------------------------------------------------------------------------------
	JC = (0.D0,1.D0)               !  SQRT(-1)
	PI = 4.0 * DATAN(1.D0)         !  PI = 3.1415
	FR = 1.D9                      !  FR = FREQUENCY [Hz]
	W = 2*PI*FR                    !  W = ANGULAR FREQUENCY
	E0 = 8.854223E-12              !  E0=PERMITTIVITY OF FREE SPACE
	LAMBDA = 3.D8/FR               !  LAMBDA = WAVELENGTH [m]
	BETA0 = 2.D0*PI/LAMBDA         !  BETA = WAVE NUMBER  [Rad/m]
C--------------------------<< ELLIPTIC CYLLINDER PARAMETERS >>---------------------
      A = LAMBDA/2.                 !  2*A = LONG DIAMETER
	B = LAMBDA/4.                 !  2*B = SHORT DIAMETER
	F = DSQRT(A*A-B*B)            !  2*F = CONFOCAL DISTANCE
	ER = 5.0D0                    !  DIELECTRIC CONSTANT OF ELLIPTIC CYLINDER
	SIGMA = 0.1D0			      !  CONDUCTIVITY OF ELLIPTIC CYLINDER
	BETA = BETA0*CDSQRT(ER-JC*SIGMA/(W*E0))
	U1 = 0.5*DLOG((A+B)/(A-B))	   ! U1=ELLIPTIC BOUNDARY PARAMETER
	Q2 = F*F*BETA0*BETA0/4		    
	Q1 = F*F*BETA*BETA/4
	WRITE(*,47)Q2,Q1
  47  FORMAT(3X,'Q2=',2F10.5,3X,'Q1=',2F15.5)
	PHI0 = 45.0
	MAX = 15 + INT(SQRT(ABS(Q1)))
	PRINT*,MAX
C--------------------<< COMPUTATION OF A^s & A^t >>-----------------
	KF = 1          !  FOR Cem,Mcm
      DO I = 1 , 2
       DO J = 1 , 2
	  IX = (I-1)*(MAX+1)
	  DO N = 0 , MAX 
	   IX = IX + 1
	   JY = (J-1) *(MAX+1)
	   DO M = 0 , MAX
	    JY = JY + 1
	    IF(I.EQ.1.AND.J.EQ.1)THEN					 !  A11
	     KC = 1					   ! FIRST KIND
	     CALL CMTURAD(MAX,KF,KC,M,Q1,U1,F1R,D1R,F2R,D2R)  
	     CALL IMNCS(MAX,KF,M,N,Q2,Q1,CSI)		 ! COMPUTES I_mn^cc(Q2,Q1)
	     G(IX,JY) = F1R*CSI
		ELSEIF(I.EQ.1.AND.J.EQ.2)THEN				 ! A12
		 IF(M.EQ.N)THEN
	      KC = 3     ! BOTH FIRST AND SECOND KIND
	      CALL CMTURAD(MAX,KF,KC,M,Q2,U1,F1R,D1R,F2R,D2R)  
	      G(IX,JY) = -(F1R-JC*F2R)	
	     ELSE
            G(IX,JY) = (0.0D0,0.0D0)
	     ENDIF
          ELSEIF(I.EQ.2.AND.J.EQ.1)THEN			     !  A21
	     KC = 1
	     CALL CMTURAD(MAX,KF,KC,M,Q1,U1,F1R,D1R,F2R,D2R)			
	     CALL IMNCS(MAX,KF,M,N,Q2,Q1,CSI)	
	     G(IX,JY) = D1R*CSI
          ELSEIF(I.EQ.2.AND.J.EQ.2)THEN	    		  ! A22
		 IF(M.EQ.N)THEN	    
	      KC = 3
	      CALL CMTURAD(MAX,KF,KC,M,Q2,U1,F1R,D1R,F2R,D2R)   
	      G(IX,JY) = -(D1R-JC*D2R)		 	   	      
	     ELSE
            G(IX,JY) = (0.0D0,0.0D0)
	     ENDIF
	    ENDIF
	   ENDDO
	  ENDDO
	 ENDDO
	ENDDO 
C--------------------------<< FILLING THE INCIDENT COEFF >>--------------------
      DO I = 1 , 2
	 IX = (I-1)*(MAX+1)
	 DO M = 0 , MAX
	  IF(M.EQ.0)THEN
	   CJ = 1.D0
	  ELSE
	   CJ = (JC)**M
	  ENDIF
	  IX = IX + 1 
	  KC = 1
	  CALL CMTURAD(MAX,KF,KC,M,Q2,U1,F1R,D1R,F2R,D2R)       
	  CALL CMTUANG(MAX,KF,M,Q2,PHI0,CSF,CSD)	      		 	   	      	    
	  IF(I.EQ.1)THEN
	   G(IX,2*MAX+3) = 2*CJ*CSF*F1R
	  ELSE
	   G(IX,2*MAX+3) = 2*CJ*CSF*D1R
	  ENDIF
	 ENDDO
	ENDDO
C-----------------------------------------------------------------------
	PRINT*,' MATRIX INVERSION IS IN PROGRESS'
      CALL CHOLESKY(G,2*IMAX+2,2*IMAX+3,2*MAX+2,2*MAX+3,XA)
	PRINT*,' XA IS COMPUTED'
C----------------------<< INITIALIZATION >>------------------------------
	DO I = 1 , 2*IMAX+2
	 DO J = 1 , 2*IMAX+3
	  G(I,J) = (0.0D0,0.0D0)
	 ENDDO
	ENDDO
C-----------------------------------------------------------------------
C--------------------<< COMPUTATION OF B^s & B^t >>-----------------
	KF = 2                         !  FOR Sem,MSm
      DO I = 1 , 2
       DO J = 1 , 2
	  IX = (I-1)*(MAX+1)
	  DO N = 1 , MAX+1 
	   IX = IX + 1
	   JY = (J-1) *(MAX+1)
	   DO M = 1 , MAX+1
	    JY = JY + 1
	    IF(I.EQ.1.AND.J.EQ.1)THEN
	     KC = 1
	     CALL CMTURAD(MAX,KF,KC,M,Q1,U1,F1R,D1R,F2R,D2R)    
	     CALL IMNCS(MAX,KF,M,N,Q2,Q1,CSI)   ! COMPUTES I_mn^ss(Q2,Q1)
	     G(IX,JY) = F1R*CSI
		ELSEIF(I.EQ.1.AND.J.EQ.2)THEN
		 IF(M.EQ.N)THEN
	      KC = 3
	      CALL CMTURAD(MAX,KF,KC,M,Q2,U1,F1R,D1R,F2R,D2R)    
	      G(IX,JY) = -(F1R-JC*F2R)		 	   
	     ELSE
            G(IX,JY) = (0.0D0,0.0D0)
	     ENDIF
          ELSEIF(I.EQ.2.AND.J.EQ.1)THEN			    
	     KC = 1
	     CALL CMTURAD(MAX,KF,KC,M,Q1,U1,F1R,D1R,F2R,D2R)    
	     CALL IMNCS(MAX,KF,M,N,Q2,Q1,CSI)
	     G(IX,JY) = D1R*CSI
          ELSEIF(I.EQ.2.AND.J.EQ.2)THEN	    
		 IF(M.EQ.N)THEN	    
	      KC = 3
	      CALL CMTURAD(MAX,KF,KC,M,Q2,U1,F1R,D1R,F2R,D2R)   
	      G(IX,JY) = -(D1R-JC*D2R)		 	   	      
	     ELSE
            G(IX,JY) = (0.0D0,0.0D0)
	     ENDIF
	    ENDIF
	   ENDDO
	  ENDDO
	 ENDDO
	ENDDO 
C--------------------------<< FILLING THE INCIDENT COEFF >>--------------------
      DO I = 1 , 2
	 IX = (I-1)*(MAX+1)
	 DO M = 1 , MAX+1
	  CJ = (JC)**M
	  IX = IX + 1 
	  KC = 1				  
	  CALL CMTURAD(MAX,KF,KC,M,Q2,U1,F1R,D1R,F2R,D2R)
	  CALL CMTUANG(MAX,KF,M,Q2,PHI0,CSF,CSD)	      		 	   	      	    
	  IF(I.EQ.1)THEN
	   G(IX,2*MAX+3) = 2*CJ*CSF*F1R
	  ELSE
	   G(IX,2*MAX+3) = 2*CJ*CSF*D1R
	  ENDIF
	 ENDDO
	ENDDO
C-----------------------------------------------------------------------
	PRINT*,' MATRIX INVERSION IS IN PROGRESS'
      CALL CHOLESKY(G,2*IMAX+2,2*IMAX+3,2*MAX+2,2*MAX+3,XB)
	PRINT*,' XB IS COMPUTED'
C-----------------------------------------------------------------------
	OPEN(UNIT= 1,FILE='cx.p')   ! OPEN CX.P FILE
	OPEN(UNIT= 2,FILE='cy.p')
	DO ID = 1 , 2			   ! AXIS INDEX
	 IF(ID.EQ.1)I = 1000*A	   ! ALONG X AXIS
	 IF(ID.EQ.2)I = 1000*B	   ! ALONG Y AXIS
	 DO IU = -I , I	  , 4
	  IF(ID.EQ.1)THEN
	   X1 = IU*0.001
	   Y1 = 0.0D0
	   R = X1
	  ELSE
	   Y1 = IU*0.001D0
	   X1 = 0.0D0
	   R = Y1
	  ENDIF
	  CALL TREC2EL(X1,Y1,F,U,V)  ! V=DEGREE
	  KC = 1   ! FIRST KIND
	  KF = 1   ! Mcm(u,Q1)	
	  EZA = (0.0D0,0.0D0)
	  DO M = 0 , MAX
	   CALL CMTUANG(MAX,KF,M,Q1,V,CSF,CSD)	      		! ce_m(u,Q1)	   	      	    
         CALL CMTURAD(MAX,KF,KC,M,Q1,U,F1R,D1R,F2R,D2R)   ! Mcm(u,Q1)
	   EZA = EZA + XA(M+1)*F1R*CSF
	  ENDDO
	  KF = 2   ! Msm(u,Q1)
	  EZB = (0.0D0,0.0D0)
	  DO M = 1 , MAX+1	     
	   CALL CMTUANG(MAX,KF,M,Q1,V,CSF,CSD)	          ! se_m(v,Q1)     		 	   	      	    
         CALL CMTURAD(MAX,KF,KC,M,Q1,U,F1R,D1R,F2R,D2R) ! Msm(u,Q1)(1)
	   EZB = EZB + XB(M)*F1R*CSF
	  ENDDO
	  EZ = EZA + EZB
	  IF(ID.EQ.1)THEN
	   WRITE(1,31) 100*R,CDABS(EZ)
	   WRITE(*,31) 100*R,CDABS(EZ)
	  ELSE
	   WRITE(2,31) 100*R,CDABS(EZ)
	   WRITE(*,31) 100*R,CDABS(EZ)
	  ENDIF
	 ENDDO
	ENDDO
C----------------<< RCS OF ELLIPTIC DIELECTRIC  CYLINDER  >>-------------------      
      OPEN(3,FILE='S1.P')
	R = 5*LAMBDA 
	DO I = 1 ,  360 , 2
	 PHI = REAL(I)  
	 CALL POL2REC(R,PHI,X1,Y1)
	 CALL TREC2EL(X1,Y1,F,U,V)	 	 
	 KF = 1   ! cm(v,Q2)	
	 EZA = (0.0D0,0.0D0)
	 DO M = 0 , MAX+1
	  CALL CMTUANG(MAX,KF,M,Q2,V,CSF,CSD)	      		 	   	      	    
	  RES = CDEXP(-JC*M*PI/2)        
	  EZA = EZA + XA(MAX+2+M)*RES*CSF
	 ENDDO
	 KF = 2   ! sm(v,Q2)
	 EZB = (0.0D0,0.0D0)
	 DO M = 1 , MAX+1	     
	  CALL CMTUANG(MAX,KF,M,Q2,V,CSF,CSD)	      		 	   	      	    
     	  RES = CDEXP(-JC*M*PI/2)
	  EZB = EZB + XB(MAX+1+M)*RES*CSF
	 ENDDO
	 EZ = EZA + EZB	
	 WRITE(*,31) PHI , 4*CDABS(EZ)**2/BETA0/LAMBDA
	 WRITE(3,31) PHI , 4*CDABS(EZ)**2/BETA0/LAMBDA
	ENDDO
 31   FORMAT(3X,4F10.4)
	STOP
	END
C*******************************************************************************
C     THIS SUBROUTINE COMPUTES THE MULTIPLICATION OF TWO COEFFICENTS     
	SUBROUTINE IMNCS(MAX,KF,M,N,Q1,Q2,CSI)
	IMPLICIT NONE
	INTEGER KD,M,N,KF,K,MAX
	COMPLEX*16 A,FC1(251),FC2(251),Q1,Q2,CSI,LAMBDA(MAX)
	REAL*8	DELTA
C                KD --- Case code
C                       KD=1 for cem(x,q)  ( m = 0,2,4,...)
C                       KD=2 for cem(x,q)  ( m = 1,3,5,...)
C                       KD=3 for sem(x,q)  ( m = 1,3,5,...)
C                       KD=4 for sem(x,q)  ( m = 2,4,6,...)
C                A  --- Characteristic value of Mathieu
C                       functions for given m and q
C       Output:  FC(k) --- Expansion coefficients of Mathieu

	IF (KF.EQ.1.AND.M.EQ.2*INT(M/2)) KD=1
	IF (KF.EQ.1.AND.M.NE.2*INT(M/2)) KD=2
	IF (KF.EQ.2.AND.M.NE.2*INT(M/2)) KD=3
	IF (KF.EQ.2.AND.M.EQ.2*INT(M/2)) KD=4

      CALL CV(MAX,KD,Q2,LAMBDA)  	! A = computing the characteristic values
	IF(KD.EQ.4)THEN
	 A = LAMBDA(INT(M/2))
	ELSE
	 A = LAMBDA(INT(M/2)+1)
	ENDIF	
	CALL FCOEF(KD,M,Q2,A,FC1)	! FC = computing the expansion coefficients
	IF (KF.EQ.1.AND.N.EQ.2*INT(N/2)) KD=1
	IF (KF.EQ.1.AND.N.NE.2*INT(N/2)) KD=2
	IF (KF.EQ.2.AND.N.NE.2*INT(N/2)) KD=3
	IF (KF.EQ.2.AND.N.EQ.2*INT(N/2)) KD=4

      CALL CV(MAX,KD,Q1,LAMBDA)  	! A = computing the characteristic values
	  IF(KD.EQ.4)THEN
	   A = LAMBDA(INT(N/2))
	  ELSE
	   A = LAMBDA(INT(N/2)+1)
	  ENDIF		
	CALL FCOEF(KD,N,Q1,A,FC2)

	CSI = (0.0D0,0.0D0)
	DO K = 1 , 2*MAX
	 IF ((MOD(M,2).EQ.0).AND.(MOD(N,2).EQ.0)) THEN
	  IF(KF.EQ.1)THEN
         CSI = CSI + (1.D0+DELTA(1,K))*FC1(K)*FC2(K)
	  ELSE
          CSI = CSI + FC1(K)*FC2(K)
	  ENDIF
	 ELSEIF ((MOD(M,2).EQ.1).AND.(MOD(N,2).EQ.1)) THEN
	   CSI = CSI + FC1(K)*FC2(K)
	 ELSE
	  CSI = (0.0D0,0.0D0)
	 ENDIF
	ENDDO
	RETURN
	END
C*******************************************************************************
      SUBROUTINE POL2REC(R,PHI,X,Y)
	REAL*8 R,PHI,X,Y
	X = R*DCOSD(PHI)
	Y = R*DSIND(PHI)
	RETURN
	END
C*******************************************************************************
      SUBROUTINE TEL2REC(U,V,F,X,Y)
	REAL*8 U,V,F,X,Y
	X = F*DCOSH(U)*DCOSD(V)
	Y = F*DSINH(U)*DSIND(V)
	RETURN
	END
C*******************************************************************************
      SUBROUTINE TREC2EL(X,Y,F,U,V)
	USE MSIMSLMD
c	use numerical_libraries

	REAL*8 U,V,F,X,Y
      U = 1.0D0/(2.0D0*F)*(DSQRT(Y*Y+(X+F)**2)+
     &   				 DSQRT(Y*Y+(X-F)**2) )
	 if(u.lt.1.0) u=1.0d0
	U = DACOSH(U)
	IF(Y.GE.0.0D0)THEN	 
	 V  = 1.0D0/(2.0D0*F)*(DSQRT(Y*Y+(X+F)**2)-
     &				  DSQRT(Y*Y+(X-F)**2) )
       V = DACOSD(V)
	ELSE	 
	 V  = 1.0D0/(2.0D0*F)*(DSQRT(Y*Y+(X+F)**2)-
     &				  DSQRT(Y*Y+(X-F)**2) )
       V = 360.D0 - DACOSD(V)
	ENDIF
	RETURN
	END
C*******************************************************************************
	REAL*8 FUNCTION DELTA(M,N)
	INTEGER M,N
	IF (M.EQ.N) THEN
	 DELTA = 1.0D0
	ELSE
	 DELTA = 0.0D0
	END IF
	RETURN
	END
C*******************************************************************************
        SUBROUTINE CMTUANG(MAX,KF,M,Q,X,CSF,CSD)
C
C       ===============================================================
C       Purpose: Compute Mathieu functions cem(x,q) and sem(x,q)
C                and their derivatives ( q ò 0 )
C       Input :  KF  --- Function code
C                        KF=1 for computing cem(x,q) and cem'(x,q)
C                        KF=2 for computing sem(x,q) and sem'(x,q)
C                m   --- Order of Mathieu functions
C                q   --- Parameter of Mathieu functions
C                x   --- Argument of Mathieu functions (in degrees)
C       Output:  CSF --- cem(x,q) or sem(x,q)
C                CSD --- cem'x,q) or sem'x,q)
C       Routines called:
C            (1) CVA2 for computing the characteristic values
C            (2) FCOEF for computing the expansion coefficients
C       ===============================================================
C
        IMPLICIT NONE
	  INTEGER KF,M,KD,KM,IC,K,MAX
	  COMPLEX*16 FG,Q,QM,A,CSF,CSD,LAMBDA(MAX)
        DIMENSION FG(251)
	  REAL*8 EPS,RD,XR,X
 	  EPS=1.0D-14
	  IF (KF.EQ.1.AND.M.EQ.2*INT(M/2)) KD=1
	  IF (KF.EQ.1.AND.M.NE.2*INT(M/2)) KD=2
	  IF (KF.EQ.2.AND.M.NE.2*INT(M/2)) KD=3
	  IF (KF.EQ.2.AND.M.EQ.2*INT(M/2)) KD=4
	  CALL CV(MAX,KD,Q,LAMBDA)
	  IF(KD.EQ.4)THEN
	   A = LAMBDA(INT(M/2))
	  ELSE
	   A = LAMBDA(INT(M/2)+1)
	  ENDIF	

	  IF (CDABS(Q).LE.1.0D0) THEN
	   QM=7.5+56.1*CDSQRT(Q)-134.7*Q+90.7*CDSQRT(Q)*Q
	  ELSE
	   QM=17.0+3.1*CDSQRT(Q)-.126*Q+.0037*CDSQRT(Q)*Q
	  ENDIF	   
        KM=INT(QM+0.5*M)
        CALL FCOEF(KD,M,Q,A,FG)
        IC=INT(M/2)+1
        RD=1.74532925199433D-2
        XR=X*RD
        CSF = (0.0D0,0.0D0)
        DO 10 K=1,KM
           IF (KD.EQ.1) THEN
              CSF=CSF+FG(K)*DCOS((2*K-2)*XR)
           ELSE IF (KD.EQ.2) THEN
              CSF=CSF+FG(K)*DCOS((2*K-1)*XR)
           ELSE IF (KD.EQ.3) THEN
              CSF=CSF+FG(K)*DSIN((2*K-1)*XR)
           ELSE IF (KD.EQ.4) THEN
              CSF=CSF+FG(K)*DSIN(2*K*XR)
           ENDIF
           IF (K.GE.IC.AND.CDABS(FG(K)).LT.CDABS(CSF)*EPS) GO TO 15
10         CONTINUE
15      CSD = (0.0D0,0.0D0)
        DO 20 K=1,KM
           IF (KD.EQ.1) THEN
              CSD=CSD-(2*K-2)**2*FG(K)*DCOS((2*K-2)*XR)
           ELSE IF (KD.EQ.2) THEN
              CSD=CSD-(2*K-1)**2*FG(K)*DCOS((2*K-1)*XR)
           ELSE IF (KD.EQ.3) THEN
              CSD=CSD-(2*K-1)**2*FG(K)*DSIN((2*K-1)*XR)
           ELSE IF (KD.EQ.4) THEN
              CSD=CSD-4.0D0*K*K*FG(K)*DSIN(2*K*XR)
           ENDIF
           IF (K.GE.IC.AND.CDABS(FG(K)).LT.CDABS(CSD)*EPS) GO TO 25
20         CONTINUE
25      RETURN
        END
C******************************************************************************
	SUBROUTINE CMTURAD(MAX,KF,KC,M,Q,X,F1R,D1R,F2R,D2R)
C
C       ==============================================================
C       Purpose: Compute modified Mathieu functions of the first and
C                second kinds, Mcm(1)(2)(x,q) and Msm(1)(2)(x,q),
C                and their derivatives
C       Input:   KF --- Function code
C                       KF=1 for computing Mcm(x,q)
C                       KF=2 for computing Msm(x,q)
C                KC --- Function Code
C                       KC=1 for computing the first kind
C                       KC=2 for computing the second kind
C                            or Msm(2)(x,q) and Msm(2)'(x,q)
C                       KC=3 for computing both the first
C                            and second kinds
C                m  --- Order of Mathieu functions
C                q  --- Parameter of Mathieu functions ( q ò 0 )
C                x  --- Argument of Mathieu functions
C       Output:  F1R --- Mcm(1)(x,q) or Msm(1)(x,q)
C                D1R --- Derivative of Mcm(1)(x,q) or Msm(1)(x,q)
C                F2R --- Mcm(2)(x,q) or Msm(2)(x,q)
C                D2R --- Derivative of Mcm(2)(x,q) or Msm(2)(x,q)
C       Routines called:
C            (1) CVA2 for computing the characteristic values
C            (2) FCOEF for computing expansion coefficients
C            (3) JYNB for computing Jn(x), Yn(x) and their
C                derivatives
C       ==============================================================
C
	IMPLICIT NONE
	INTEGER MAX
	COMPLEX*16 FG,BJ1,DJ1,BJ2,DJ2,BY1,DY1,BY2,DY2,Q
	COMPLEX*16 C1,C2,U1,U2,QM,A,F1R,W1,D1R,W2,F2R,D2R
	COMPLEX*16 LAMBDA(MAX)
	REAL*8 EPS,X
	DIMENSION FG(251),BJ1(0:251),DJ1(0:251),BJ2(0:251),DJ2(0:251),
     &            BY1(0:251),DY1(0:251),BY2(0:251),DY2(0:251)
	INTEGER KM,IC,M,KD,NM,KF,KC,K
	EPS=1.0D-14
	IF (KF.EQ.1.AND.M.EQ.2*INT(M/2)) KD=1
	IF (KF.EQ.1.AND.M.NE.2*INT(M/2)) KD=2
	IF (KF.EQ.2.AND.M.NE.2*INT(M/2)) KD=3
	IF (KF.EQ.2.AND.M.EQ.2*INT(M/2)) KD=4
	CALL CV(MAX,KD,Q,LAMBDA)
	IF(KD.EQ.4)THEN
	 A = LAMBDA(INT(M/2))
	ELSE
	 A = LAMBDA(INT(M/2)+1)
	ENDIF	 
	IF (CDABS(Q).LE.1.0D0) THEN
	   QM=7.5+56.1*CDSQRT(Q)-134.7*Q+90.7*CDSQRT(Q)*Q
	ELSE
	   QM=17.0+3.1*CDSQRT(Q)-.126*Q+.0037*CDSQRT(Q)*Q
	ENDIF
	KM=INT(QM+0.5*M)              
	CALL FCOEF(KD,M,Q,A,FG)
	IC=INT(M/2)+1
	IF (KD.EQ.4) IC=M/2
	C1 = DEXP(-X)
	C2 = DEXP(X)
	U1 = CDSQRT(Q)*C1
	U2 = CDSQRT(Q)*C2
	CALL CJYNA(KM,U1,NM,BJ1,DJ1,BY1,DY1)
	CALL CJYNA(KM,U2,NM,BJ2,DJ2,BY2,DY2)
	IF (KC.EQ.2) GO TO 50
	F1R = (0.0D0,0.0D0)
	DO 30 K=1,KM
	   IF (KD.EQ.1) THEN
	      F1R=F1R+(-1)**(IC+K)*FG(K)*BJ1(K-1)*BJ2(K-1)
	   ELSE IF (KD.EQ.2.OR.KD.EQ.3) THEN
	      F1R=F1R+(-1)**(IC+K)*FG(K)*(BJ1(K-1)*BJ2(K)
     &            +(-1)**KD*BJ1(K)*BJ2(K-1))
	   ELSE
	      F1R=F1R+(-1)**(IC+K)*FG(K)*(BJ1(K-1)*BJ2(K+1)
     &            -BJ1(K+1)*BJ2(K-1))
	   ENDIF
	   IF (K.GE.5.AND.CDABS(F1R-W1).LT.CDABS(F1R)*EPS) GO TO 35
30         W1=F1R
35      F1R=F1R/FG(1)
	D1R=(0.0D0,0.0D0)
	DO 40 K=1,KM
	   IF (KD.EQ.1) THEN
	      D1R=D1R+(-1)**(IC+K)*FG(K)*(C2*BJ1(K-1)*DJ2(K-1)
     &            -C1*DJ1(K-1)*BJ2(K-1))
	   ELSE IF (KD.EQ.2.OR.KD.EQ.3) THEN
	      D1R=D1R+(-1)**(IC+K)*FG(K)*(C2*(BJ1(K-1)*DJ2(K)
     &            +(-1)**KD*BJ1(K)*DJ2(K-1))-C1*(DJ1(K-1)*BJ2(K)
     &            +(-1)**KD*DJ1(K)*BJ2(K-1)))
	   ELSE
	      D1R=D1R+(-1)**(IC+K)*FG(K)*(C2*(BJ1(K-1)*DJ2(K+1)
     &            -BJ1(K+1)*DJ2(K-1))-C1*(DJ1(K-1)*BJ2(K+1)
     &            -DJ1(K+1)*BJ2(K-1)))
	   ENDIF
	   IF (K.GE.5.AND.CDABS(D1R-W2).LT.CDABS(D1R)*EPS) GO TO 45
40         W2=D1R
45      D1R=D1R*CDSQRT(Q)/FG(1)
	IF (KC.EQ.1) RETURN
50      F2R=(0.0D0,0.0D0)
	DO 55 K=1,KM
	   IF (KD.EQ.1) THEN
	      F2R=F2R+(-1)**(IC+K)*FG(K)*BJ1(K-1)*BY2(K-1)
	   ELSE IF (KD.EQ.2.OR.KD.EQ.3) THEN
	      F2R=F2R+(-1)**(IC+K)*FG(K)*(BJ1(K-1)*BY2(K)
     &            +(-1)**KD*BJ1(K)*BY2(K-1))
	   ELSE
	      F2R=F2R+(-1)**(IC+K)*FG(K)*(BJ1(K-1)*BY2(K+1)
     &            -BJ1(K+1)*BY2(K-1))
	   ENDIF
	   IF (K.GE.5.AND.CDABS(F2R-W1).LT.CDABS(F2R)*EPS) GO TO 60
55         W1=F2R
60      F2R=F2R/FG(1)
	D2R=(0.0D0,0.0D0)
	DO 65 K=1,KM
	   IF (KD.EQ.1) THEN
	      D2R=D2R+(-1)**(IC+K)*FG(K)*(C2*BJ1(K-1)*DY2(K-1)
     &            -C1*DJ1(K-1)*BY2(K-1))
	   ELSE IF (KD.EQ.2.OR.KD.EQ.3) THEN
	      D2R=D2R+(-1)**(IC+K)*FG(K)*(C2*(BJ1(K-1)*DY2(K)
     &            +(-1)**KD*BJ1(K)*DY2(K-1))-C1*(DJ1(K-1)*BY2(K)
     &            +(-1)**KD*DJ1(K)*BY2(K-1)))
	   ELSE
	      D2R=D2R+(-1)**(IC+K)*FG(K)*(C2*(BJ1(K-1)*DY2(K+1)
     &            -BJ1(K+1)*DY2(K-1))-C1*(DJ1(K-1)*BY2(K+1)
     &            -DJ1(K+1)*BY2(K-1)))
	   ENDIF
	   IF (K.GE.5.AND.CDABS(D2R-W2).LT.CDABS(D2R)*EPS) GO TO 70
65         W2=D2R
70         D2R=D2R*CDSQRT(Q)/FG(1)
	RETURN
	END
C*******************************************************************************
        SUBROUTINE CJYNA(N,Z,NM,CBJ,CDJ,CBY,CDY)
C
C       =======================================================
C       Purpose: Compute Bessel functions Jn(z), Yn(z) and
C                their derivatives for a complex argument
C       Input :  z --- Complex argument of Jn(z) and Yn(z)
C                n --- Order of Jn(z) and Yn(z)
C       Output:  CBJ(n) --- Jn(z)
C                CDJ(n) --- Jn'(z)
C                CBY(n) --- Yn(z)
C                CDY(n) --- Yn'(z)
C                NM --- Highest order computed
C       Rouitines called:
C            (1) CJY01 to calculate J0(z), J1(z), Y0(z), Y1(z)
C            (2) MSTA1 and MSTA2 to calculate the starting 
C                point for backward recurrence
C       =======================================================
C
        IMPLICIT DOUBLE PRECISION (A,B,E,P,R,W,Y)
        IMPLICIT COMPLEX*16 (C,Z)
        DIMENSION CBJ(0:N),CDJ(0:N),CBY(0:N),CDY(0:N)
        PI=3.141592653589793D0
        A0=CDABS(Z)
        NM=N
        IF (A0.LT.1.0D-100) THEN
           DO 5 K=0,N
              CBJ(K)=(0.0D0,0.0D0)
              CDJ(K)=(0.0D0,0.0D0)
              CBY(K)=-(1.0D+300,0.0D0)
5             CDY(K)=(1.0D+300,0.0D0)
           CBJ(0)=(1.0D0,0.0D0)
           CDJ(1)=(0.5D0,0.0D0)
           RETURN
        ENDIF
        CALL CJY01(Z,CBJ0,CDJ0,CBJ1,CDJ1,CBY0,CDY0,CBY1,CDY1)
        CBJ(0)=CBJ0
        CBJ(1)=CBJ1
        CBY(0)=CBY0
        CBY(1)=CBY1
        CDJ(0)=CDJ0
        CDJ(1)=CDJ1
        CDY(0)=CDY0
        CDY(1)=CDY1
        IF (N.LE.1) RETURN
        IF (N.LT.INT(0.25*A0)) THEN
           CJ0=CBJ0
           CJ1=CBJ1
           DO 70 K=2,N
              CJK=2.0D0*(K-1.0D0)/Z*CJ1-CJ0
              CBJ(K)=CJK
              CJ0=CJ1
70            CJ1=CJK
        ELSE
           M=MSTA1(A0,200)
           IF (M.LT.N) THEN
              NM=M
           ELSE
              M=MSTA2(A0,N,15)
           ENDIF
           CF2=(0.0D0,0.0D0)
           CF1=(1.0D-100,0.0D0)
           DO 75 K=M,0,-1
              CF=2.0D0*(K+1.0D0)/Z*CF1-CF2
              IF (K.LE.NM) CBJ(K)=CF
              CF2=CF1
75            CF1=CF
           IF (CDABS(CBJ0).GT.CDABS(CBJ1)) THEN
              CS=CBJ0/CF
           ELSE
              CS=CBJ1/CF2
           ENDIF
           DO 80 K=0,NM
80            CBJ(K)=CS*CBJ(K)
        ENDIF
        DO 85 K=2,NM
85         CDJ(K)=CBJ(K-1)-K/Z*CBJ(K)
        YA0=CDABS(CBY0)
        LB=0
        CG0=CBY0
        CG1=CBY1
        DO 90 K=2,NM
           CYK=2.0D0*(K-1.0D0)/Z*CG1-CG0
           IF (CDABS(CYK).GT.1.0D+290) GO TO 90            
           YAK=CDABS(CYK)
           YA1=CDABS(CG0)
           IF (YAK.LT.YA0.AND.YAK.LT.YA1) LB=K
           CBY(K)=CYK
           CG0=CG1
           CG1=CYK
90      CONTINUE
        IF (LB.LE.4.OR.DIMAG(Z).EQ.0.0D0) GO TO 125
95      IF (LB.EQ.LB0) GO TO 125
        CH2=(1.0D0,0.0D0)
        CH1=(0.0D0,0.0D0)
        LB0=LB
        DO 100 K=LB,1,-1
           CH0=2.0D0*K/Z*CH1-CH2
           CH2=CH1
100        CH1=CH0
        CP12=CH0
        CP22=CH2
        CH2=(0.0D0,0.0D0)
        CH1=(1.0D0,0.0D0)
        DO 105 K=LB,1,-1
           CH0=2.0D0*K/Z*CH1-CH2
           CH2=CH1
105        CH1=CH0
        CP11=CH0
        CP21=CH2
        IF (LB.EQ.NM) CBJ(LB+1)=2.0D0*LB/Z*CBJ(LB)-CBJ(LB-1)
        IF (CDABS(CBJ(0)).GT.CDABS(CBJ(1))) THEN
           CBY(LB+1)=(CBJ(LB+1)*CBY0-2.0D0*CP11/(PI*Z))/CBJ(0)
           CBY(LB)=(CBJ(LB)*CBY0+2.0D0*CP12/(PI*Z))/CBJ(0)
        ELSE
           CBY(LB+1)=(CBJ(LB+1)*CBY1-2.0D0*CP21/(PI*Z))/CBJ(1)
           CBY(LB)=(CBJ(LB)*CBY1+2.0D0*CP22/(PI*Z))/CBJ(1)
        ENDIF
        CYL2=CBY(LB+1)
        CYL1=CBY(LB)
        DO 110 K=LB-1,0,-1
           CYLK=2.0D0*(K+1.0D0)/Z*CYL1-CYL2
           CBY(K)=CYLK
           CYL2=CYL1
110        CYL1=CYLK
        CYL1=CBY(LB)
        CYL2=CBY(LB+1)
        DO 115 K=LB+1,NM-1
           CYLK=2.0D0*K/Z*CYL2-CYL1
           CBY(K+1)=CYLK
           CYL1=CYL2
115        CYL2=CYLK
        DO 120 K=2,NM
           WA=CDABS(CBY(K))
           IF (WA.LT.CDABS(CBY(K-1))) LB=K
120     CONTINUE
        GO TO 95
125     CONTINUE
        DO 130 K=2,NM
130        CDY(K)=CBY(K-1)-K/Z*CBY(K)
        RETURN
        END
C*************************************************************************
        SUBROUTINE CJY01(Z,CBJ0,CDJ0,CBJ1,CDJ1,CBY0,CDY0,CBY1,CDY1)
C
C       ===========================================================
C       Purpose: Compute complex Bessel functions J0(z), J1(z)
C                Y0(z), Y1(z), and their derivatives
C       Input :  z --- Complex argument
C       Output:  CBJ0 --- J0(z)
C                CDJ0 --- J0'(z)
C                CBJ1 --- J1(z)
C                CDJ1 --- J1'(z)
C                CBY0 --- Y0(z)
C                CDY0 --- Y0'(z)
C                CBY1 --- Y1(z)
C                CDY1 --- Y1'(z)
C       ===========================================================
C
        IMPLICIT DOUBLE PRECISION (A,B,E,P,R,W)
        IMPLICIT COMPLEX*16 (C,Z)
        DIMENSION A(12),B(12),A1(12),B1(12)
        PI=3.141592653589793D0
        EL=0.5772156649015329D0
        RP2=2.0D0/PI
        CI=(0.0D0,1.0D0)
        A0=CDABS(Z)
        Z2=Z*Z
        Z1=Z
        IF (A0.EQ.0.0D0) THEN
           CBJ0=(1.0D0,0.0D0)
           CBJ1=(0.0D0,0.0D0)
           CDJ0=(0.0D0,0.0D0)
           CDJ1=(0.5D0,0.0D0)
           CBY0=-(1.0D300,0.0D0)
           CBY1=-(1.0D300,0.0D0)
           CDY0=(1.0D300,0.0D0)
           CDY1=(1.0D300,0.0D0)
           RETURN
        ENDIF
        IF (REAL(Z).LT.0.0) Z1=-Z
        IF (A0.LE.12.0) THEN
           CBJ0=(1.0D0,0.0D0)
           CR=(1.0D0,0.0D0)
           DO 10 K=1,40
              CR=-0.25D0*CR*Z2/(K*K)
              CBJ0=CBJ0+CR
              IF (CDABS(CR/CBJ0).LT.1.0D-15) GO TO 15
10         CONTINUE
15         CBJ1=(1.0D0,0.0D0)
           CR=(1.0D0,0.0D0)
           DO 20 K=1,40
              CR=-0.25D0*CR*Z2/(K*(K+1.0D0))
              CBJ1=CBJ1+CR
              IF (CDABS(CR/CBJ1).LT.1.0D-15) GO TO 25
20         CONTINUE
25         CBJ1=0.5D0*Z1*CBJ1
           W0=0.0D0
           CR=(1.0D0,0.0D0)
           CS=(0.0D0,0.0D0)
           DO 30 K=1,40
              W0=W0+1.0D0/K
              CR=-0.25D0*CR/(K*K)*Z2
              CP=CR*W0
              CS=CS+CP
              IF (CDABS(CP/CS).LT.1.0D-15) GO TO 35
30         CONTINUE
35         CBY0=RP2*(CDLOG(Z1/2.0D0)+EL)*CBJ0-RP2*CS
           W1=0.0D0
           CR=(1.0D0,0.0D0)
           CS=(1.0D0,0.0D0)
           DO 40 K=1,40
              W1=W1+1.0D0/K
              CR=-0.25D0*CR/(K*(K+1))*Z2
              CP=CR*(2.0D0*W1+1.0D0/(K+1.0D0))
              CS=CS+CP
              IF (CDABS(CP/CS).LT.1.0D-15) GO TO 45
40         CONTINUE
45         CBY1=RP2*((CDLOG(Z1/2.0D0)+EL)*CBJ1-1.0D0/Z1-.25D0*Z1*CS)
        ELSE
           DATA A/-.703125D-01,.112152099609375D+00,
     &            -.5725014209747314D+00,.6074042001273483D+01,
     &            -.1100171402692467D+03,.3038090510922384D+04,
     &            -.1188384262567832D+06,.6252951493434797D+07,
     &            -.4259392165047669D+09,.3646840080706556D+11,
     &            -.3833534661393944D+13,.4854014686852901D+15/
           DATA B/ .732421875D-01,-.2271080017089844D+00,
     &             .1727727502584457D+01,-.2438052969955606D+02,
     &             .5513358961220206D+03,-.1825775547429318D+05,
     &             .8328593040162893D+06,-.5006958953198893D+08,
     &             .3836255180230433D+10,-.3649010818849833D+12,
     &             .4218971570284096D+14,-.5827244631566907D+16/
           DATA A1/.1171875D+00,-.144195556640625D+00,
     &             .6765925884246826D+00,-.6883914268109947D+01,
     &             .1215978918765359D+03,-.3302272294480852D+04,
     &             .1276412726461746D+06,-.6656367718817688D+07,
     &             .4502786003050393D+09,-.3833857520742790D+11,
     &             .4011838599133198D+13,-.5060568503314727D+15/
           DATA B1/-.1025390625D+00,.2775764465332031D+00,
     &             -.1993531733751297D+01,.2724882731126854D+02,
     &             -.6038440767050702D+03,.1971837591223663D+05,
     &             -.8902978767070678D+06,.5310411010968522D+08,
     &             -.4043620325107754D+10,.3827011346598605D+12,
     &             -.4406481417852278D+14,.6065091351222699D+16/
           K0=12
           IF (A0.GE.35.0) K0=10
           IF (A0.GE.50.0) K0=8
           CT1=Z1-0.25D0*PI
           CP0=(1.0D0,0.0D0)
           DO 50 K=1,K0
50            CP0=CP0+A(K)*Z1**(-2*K)
           CQ0=-0.125D0/Z1
           DO 55 K=1,K0
55            CQ0=CQ0+B(K)*Z1**(-2*K-1)
           CU=CDSQRT(RP2/Z1)
           CBJ0=CU*(CP0*CDCOS(CT1)-CQ0*CDSIN(CT1))
           CBY0=CU*(CP0*CDSIN(CT1)+CQ0*CDCOS(CT1))
           CT2=Z1-0.75D0*PI
           CP1=(1.0D0,0.0D0)
           DO 60 K=1,K0
60            CP1=CP1+A1(K)*Z1**(-2*K)
           CQ1=0.375D0/Z1
           DO 65 K=1,K0
65            CQ1=CQ1+B1(K)*Z1**(-2*K-1)
           CBJ1=CU*(CP1*CDCOS(CT2)-CQ1*CDSIN(CT2))
           CBY1=CU*(CP1*CDSIN(CT2)+CQ1*CDCOS(CT2))
        ENDIF
        IF (REAL(Z).LT.0.0) THEN
           IF (DIMAG(Z).LT.0.0) CBY0=CBY0-2.0D0*CI*CBJ0
           IF (DIMAG(Z).GT.0.0) CBY0=CBY0+2.0D0*CI*CBJ0
           IF (DIMAG(Z).LT.0.0) CBY1=-(CBY1-2.0D0*CI*CBJ1)
           IF (DIMAG(Z).GT.0.0) CBY1=-(CBY1+2.0D0*CI*CBJ1)
           CBJ1=-CBJ1
        ENDIF
        CDJ0=-CBJ1
        CDJ1=CBJ0-1.0D0/Z*CBJ1
        CDY0=-CBY1
        CDY1=CBY0-1.0D0/Z*CBY1
        RETURN
        END
C*************************************************************************
        SUBROUTINE CJYNB(N,Z,NM,CBJ,CDJ,CBY,CDY)
C				   
C       =======================================================
C       Purpose: Compute Bessel functions Jn(z), Yn(z) and
C                their derivatives for a complex argument
C       Input :  z --- Complex argument of Jn(z) and Yn(z)
C                n --- Order of Jn(z) and Yn(z)
C       Output:  CBJ(n) --- Jn(z)
C                CDJ(n) --- Jn'(z)
C                CBY(n) --- Yn(z)
C                CDY(n) --- Yn'(z)
C                NM --- Highest order computed
C       Routines called:
C                MSTA1 and MSTA2 to calculate the starting
C                point for backward recurrence
C       =======================================================
C
        IMPLICIT DOUBLE PRECISION (A,B,D-H,O-Y)
        IMPLICIT COMPLEX*16 (C,Z)
        DIMENSION CBJ(0:N),CDJ(0:N),CBY(0:N),CDY(0:N),
     &            A(4),B(4),A1(4),B1(4)
        EL=0.5772156649015329D0
        PI=3.141592653589793D0
        R2P=.63661977236758D0
        Y0=DABS(DIMAG(Z))
        A0=CDABS(Z)
        NM=N
        IF (A0.LT.1.0D-100) THEN
           DO 10 K=0,N
              CBJ(K)=(0.0D0,0.0D0)
              CDJ(K)=(0.0D0,0.0D0)
              CBY(K)=-(1.0D+300,0.0D0)
10            CDY(K)=(1.0D+300,0.0D0)
           CBJ(0)=(1.0D0,0.0D0)
           CDJ(1)=(0.5D0,0.0D0)
           RETURN
        ENDIF
        IF (A0.LE.300.D0.OR.N.GT.INT(0.25*A0)) THEN
           IF (N.EQ.0) NM=1
           M=MSTA1(A0,200)
           IF (M.LT.NM) THEN
              NM=M
           ELSE
              M=MSTA2(A0,NM,15)
           ENDIF
           CBS=(0.0D0,0.0D0)
           CSU=(0.0D0,0.0D0)
           CSV=(0.0D0,0.0D0)
           CF2=(0.0D0,0.0D0)
           CF1=(1.0D-100,0.0D0)
           DO 15 K=M,0,-1
              CF=2.0D0*(K+1.0D0)/Z*CF1-CF2
              IF (K.LE.NM) CBJ(K)=CF
              IF (K.EQ.2*INT(K/2).AND.K.NE.0) THEN
                 IF (Y0.LE.1.0D0) THEN
                    CBS=CBS+2.0D0*CF
                 ELSE
                    CBS=CBS+(-1)**(K/2)*2.0D0*CF
                 ENDIF
                 CSU=CSU+(-1)**(K/2)*CF/K
              ELSE IF (K.GT.1) THEN
                 CSV=CSV+(-1)**(K/2)*K/(K*K-1.0D0)*CF
              ENDIF
              CF2=CF1
15            CF1=CF
           IF (Y0.LE.1.0D0) THEN
              CS0=CBS+CF
           ELSE
              CS0=(CBS+CF)/CDCOS(Z)
           ENDIF
           DO 20 K=0,NM
20            CBJ(K)=CBJ(K)/CS0
           CE=CDLOG(Z/2.0D0)+EL
           CBY(0)=R2P*(CE*CBJ(0)-4.0D0*CSU/CS0)
           CBY(1)=R2P*(-CBJ(0)/Z+(CE-1.0D0)*CBJ(1)-4.0D0*CSV/CS0)
        ELSE
           DATA A/-.7031250000000000D-01,.1121520996093750D+00,
     &            -.5725014209747314D+00,.6074042001273483D+01/
           DATA B/ .7324218750000000D-01,-.2271080017089844D+00,
     &             .1727727502584457D+01,-.2438052969955606D+02/
           DATA A1/.1171875000000000D+00,-.1441955566406250D+00,
     &             .6765925884246826D+00,-.6883914268109947D+01/
           DATA B1/-.1025390625000000D+00,.2775764465332031D+00,
     &             -.1993531733751297D+01,.2724882731126854D+02/
           CT1=Z-0.25D0*PI
           CP0=(1.0D0,0.0D0)
           DO 25 K=1,4
25            CP0=CP0+A(K)*Z**(-2*K)
           CQ0=-0.125D0/Z
           DO 30 K=1,4
30            CQ0=CQ0+B(K)*Z**(-2*K-1)
           CU=CDSQRT(R2P/Z)
           CBJ0=CU*(CP0*CDCOS(CT1)-CQ0*CDSIN(CT1))
           CBY0=CU*(CP0*CDSIN(CT1)+CQ0*CDCOS(CT1))
           CBJ(0)=CBJ0
           CBY(0)=CBY0
           CT2=Z-0.75D0*PI
           CP1=(1.0D0,0.0D0)
           DO 35 K=1,4
35            CP1=CP1+A1(K)*Z**(-2*K)
           CQ1=0.375D0/Z
           DO 40 K=1,4
40            CQ1=CQ1+B1(K)*Z**(-2*K-1)
           CBJ1=CU*(CP1*CDCOS(CT2)-CQ1*CDSIN(CT2))
           CBY1=CU*(CP1*CDSIN(CT2)+CQ1*CDCOS(CT2))
           CBJ(1)=CBJ1
           CBY(1)=CBY1
           DO 45 K=2,NM
              CBJK=2.0D0*(K-1.0D0)/Z*CBJ1-CBJ0
              CBJ(K)=CBJK
              CBJ0=CBJ1
45            CBJ1=CBJK
        ENDIF
        CDJ(0)=-CBJ(1)
        DO 50 K=1,NM
50         CDJ(K)=CBJ(K-1)-K/Z*CBJ(K)
        IF (CDABS(CBJ(0)).GT.1.0D0) THEN
           CBY(1)=(CBJ(1)*CBY(0)-2.0D0/(PI*Z))/CBJ(0)
        ENDIF
        DO 55 K=2,NM
           IF (CDABS(CBJ(K-1)).GE.CDABS(CBJ(K-2))) THEN
              CYY=(CBJ(K)*CBY(K-1)-2.0D0/(PI*Z))/CBJ(K-1)
           ELSE
              CYY=(CBJ(K)*CBY(K-2)-4.0D0*(K-1.0D0)/(PI*Z*Z))/CBJ(K-2)
           ENDIF
           CBY(K)=CYY
55      CONTINUE
        CDY(0)=-CBY(1)
        DO 60 K=1,NM
60         CDY(K)=CBY(K-1)-K/Z*CBY(K)
        RETURN
        END
C*******************************************************************************
	INTEGER FUNCTION MSTA1(X,MP)
C
C       ===================================================
C       Purpose: Determine the starting point for backward  
C                recurrence such that the magnitude of    
C                Jn(x) at that point is about 10^(-MP)
C       Input :  x     --- Argument of Jn(x)
C                MP    --- Value of magnitude
C       Output:  MSTA1 --- Starting point   
C       ===================================================
C
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	A0=DABS(X)
	N0=INT(1.1*A0)+1
	F0=ENVJ(N0,A0)-MP
	N1=N0+5
	F1=ENVJ(N1,A0)-MP
	DO 10 IT=1,20             
	   NN=N1-(N1-N0)/(1.0D0-F0/F1)                  
	   F=ENVJ(NN,A0)-MP
	   IF(ABS(NN-N1).LT.1) GO TO 20
	   N0=N1
	   F0=F1
	   N1=NN
 10        F1=F
 20     MSTA1=NN
	RETURN
	END
C*******************************************************************************
	INTEGER FUNCTION MSTA2(X,N,MP)
C
C       ===================================================
C       Purpose: Determine the starting point for backward
C                recurrence such that all Jn(x) has MP
C                significant digits
C       Input :  x  --- Argument of Jn(x)
C                n  --- Order of Jn(x)
C                MP --- Significant digit
C       Output:  MSTA2 --- Starting point
C       ===================================================
C
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	A0=DABS(X)
	HMP=0.5D0*MP
	EJN=ENVJ(N,A0)
	IF (EJN.LE.HMP) THEN
	   OBJ=MP
	   N0=INT(1.1*A0)
	ELSE
	   OBJ=HMP+EJN
	   N0=N
	ENDIF
	F0=ENVJ(N0,A0)-OBJ
	N1=N0+5
	F1=ENVJ(N1,A0)-OBJ
	DO 10 IT=1,20
	   NN=N1-(N1-N0)/(1.0D0-F0/F1)
	   F=ENVJ(NN,A0)-OBJ
	   IF (ABS(NN-N1).LT.1) GO TO 20
	   N0=N1
	   F0=F1
	   N1=NN
10         F1=F
20      MSTA2=NN+10
	RETURN
	END
C*******************************************************************************
	REAL*8 FUNCTION ENVJ(N,X)
	DOUBLE PRECISION X
	ENVJ=0.5D0*DLOG10(6.28D0*N)-N*DLOG10(1.36D0*X/N)
	RETURN
	END
C*******************************************************************************
      SUBROUTINE SORT(EVAL,N)
	IMPLICIT NONE
	INTEGER N,I,J
	COMPLEX*16 EVAL(N),C
	DO I = 1 , N
	 DO J = 1 , N
	  IF( REAL(EVAL(I)) .LT. REAL(EVAL(J)) ) THEN
	   C = EVAL(I)
	   EVAL(I) = EVAL(J)
	   EVAL(J) = C
	  ENDIF
	 ENDDO
	ENDDO	
	RETURN
      END
C*******************************************************************************	 
        SUBROUTINE FCOEF(KD,M,Q,A,FC)
C
C       =====================================================
C       Purpose: Compute expansion coefficients for Mathieu
C                functions and modified Mathieu functions
C       Input :  m  --- Order of Mathieu functions
C                q  --- Parameter of Mathieu functions
C                KD --- Case code
C                       KD=1 for cem(x,q)  ( m = 0,2,4,...)
C                       KD=2 for cem(x,q)  ( m = 1,3,5,...)
C                       KD=3 for sem(x,q)  ( m = 1,3,5,...)
C                       KD=4 for sem(x,q)  ( m = 2,4,6,...)
C                A  --- Characteristic value of Mathieu
C                       functions for given m and q
C       Output:  FC(k) --- Expansion coefficients of Mathieu
C                       functions ( k= 1,2,...,KM )
C                       FC(1),FC(2),FC(3),... correspond to
C                       A0,A2,A4,... for KD=1 case, A1,A3,
C                       A5,... for KD=2 case, B1,B3,B5,...
C                       for KD=3 case and B2,B4,B6,... for
C                       KD=4 case
C       =====================================================
C
        IMPLICIT NONE
	  COMPLEX*16 FC,Q,QM,S,SS,S0,U,V,F,F1,F2,F3,A,SP
        DIMENSION FC(251)
	  INTEGER KM,K,KD,M,KB,I,J
	  DO I = 1 , 251
	   FC(I) = (0.0D0,0.0D0)
	  ENDDO
        IF (CDABS(Q).LE.1.0D0) THEN
           QM=7.5+56.1*CDSQRT(Q)-134.7*Q+90.7*CDSQRT(Q)*Q
        ELSE
           QM=17.0+3.1*CDSQRT(Q)-.126*Q+.0037*CDSQRT(Q)*Q
        ENDIF
        KM=INT(QM+0.5*M)
        IF (CDABS(Q).EQ.0.0D0) THEN
           DO 10 K=1,KM
10            FC(K)= (0.0D0,0.0D0)
           IF (KD.EQ.1) THEN
              FC((M+2)/2)=(1.0D0,0.0D0)
              IF (M.EQ.0) FC(1)= DCMPLX(1.0D0/DSQRT(2.0D0),
     &          1.0D0/DSQRT(2.0D0))
           ELSE IF (KD.EQ.4) THEN
              FC(M/2)= (1.0D0,0.0D0)
           ELSE
              FC((M+1)/2)= (1.0D0,0.0D0)
           ENDIF
           RETURN
        ENDIF
        KB=0
        S = (0.0D0,0.0D0)
        F = (1.0D-100,1.0D-100)
        U = (0.0D0,0.0D0)
        FC(KM)= (0.0D0,0.0D0)
        IF (KD.EQ.1) THEN
           DO 25 K=KM,3,-1
              V=U
              U=F
              F=(A-4.0D0*K*K)*U/Q-V
              IF (CDABS(F).LT.CDABS(FC(K+1))) THEN
                 KB=K
                 FC(1)=(1.0D-100,1.0D-100)
                 SP=(0.0D0,0.0D0)
                 F3=FC(K+1)
                 FC(2)=A/Q*FC(1)
                 FC(3)=(A-4.0D0)*FC(2)/Q-2.0D0*FC(1)
                 U=FC(2)
                 F1=FC(3)
                 DO 15 I=3,KB
                    V=U
                    U=F1
                    F1=(A-4.0D0*(I-1.0D0)**2)*U/Q-V
                    FC(I+1)=F1
                    IF (I.EQ.KB) F2=F1
                    IF (I.NE.KB) SP=SP+F1*F1
15               CONTINUE
                 SP=SP+2.0D0*FC(1)**2+FC(2)**2+FC(3)**2
                 SS=S+SP*(F3/F2)**2
                 S0=CDSQRT(1.0D0/SS)
                 DO 20 J=1,KM
                    IF (J.LE.KB+1) THEN
                       FC(J)=S0*FC(J)*F3/F2
                    ELSE
                       FC(J)=S0*FC(J)
                    ENDIF
20               CONTINUE
                 GO TO 85
              ELSE
                 FC(K)=F
                 S=S+F*F
              ENDIF
25         CONTINUE
           FC(2)=Q*FC(3)/(A-4.0D0-2.0D0*Q*Q/A)
           FC(1)=Q/A*FC(2)
           S=S+2.0D0*FC(1)**2+FC(2)**2
           S0=CDSQRT(1.0D0/S)
           DO 30 K=1,KM
30            FC(K)=S0*FC(K)
        ELSE IF (KD.EQ.2.OR.KD.EQ.3) THEN
           DO 35 K=KM,3,-1
              V=U
              U=F
              F=(A-(2.0D0*K-1)**2)*U/Q-V
              IF (CDABS(F).GE.CDABS(FC(K))) THEN
                 FC(K-1)=F
                 S=S+F*F
              ELSE
                 KB=K
                 F3=FC(K)
                 GO TO 45
              ENDIF
35         CONTINUE
           FC(1)=Q/(A-1.0D0-(-1)**KD*Q)*FC(2)
           S=S+FC(1)*FC(1)
           S0=CDSQRT(1.0D0/S)
           DO 40 K=1,KM
40            FC(K)=S0*FC(K)
           GO TO 85
45         FC(1)=1.0D-100
           FC(2)=(A-1.0D0-(-1)**KD*Q)/Q*FC(1)
           SP=0.0D0
           U=FC(1)
           F1=FC(2)
           DO 50 I=2,KB-1
              V=U
              U=F1
              F1=(A-(2.0D0*I-1.0D0)**2)*U/Q-V
              IF (I.NE.KB-1) THEN
                 FC(I+1)=F1
                 SP=SP+F1*F1
              ELSE
                 F2=F1
              ENDIF
50         CONTINUE
           SP=SP+FC(1)**2+FC(2)**2
           SS=S+SP*(F3/F2)**2
           S0=1.0D0/CDSQRT(SS)
           DO 55 J=1,KM
              IF (J.LT.KB) FC(J)=S0*FC(J)*F3/F2
              IF (J.GE.KB) FC(J)=S0*FC(J)
55         CONTINUE
        ELSE IF (KD.EQ.4) THEN
           DO 60 K=KM,3,-1
              V=U
              U=F
              F=(A-4.0D0*K*K)*U/Q-V
              IF (CDABS(F).GE.CDABS(FC(K))) THEN
                 FC(K-1)=F
                 S=S+F*F
              ELSE
                 KB=K
                 F3=FC(K)
                 GO TO 70
              ENDIF
60         CONTINUE
           FC(1)=Q/(A-4.0D0)*FC(2)
           S=S+FC(1)*FC(1)
           S0=CDSQRT(1.0D0/S)
           DO 65 K=1,KM
65            FC(K)=S0*FC(K)
           GO TO 85
70         FC(1)= (1.0D-100,1.0D-100)
           FC(2)=(A-4.0D0)/Q*FC(1)
           SP=0.0D0
           U=FC(1)
           F1=FC(2)
           DO 75 I=2,KB-1
              V=U
              U=F1
              F1=(A-4.0D0*I*I)*U/Q-V
              IF (I.NE.KB-1) THEN
                 FC(I+1)=F1
                 SP=SP+F1*F1
              ELSE
                 F2=F1
              ENDIF
75         CONTINUE
           SP=SP+FC(1)**2+FC(2)**2
           SS=S+SP*(F3/F2)**2
           S0=1.0D0/CDSQRT(SS)
           DO 80 J=1,KM
              IF (J.LT.KB) FC(J)=S0*FC(J)*F3/F2
              IF (J.GE.KB) FC(J)=S0*FC(J)
80         CONTINUE
        ENDIF
85      IF (REAL(FC(1)).LT.0.0D0) THEN
           DO 90 J=1,KM
90            FC(J)=-FC(J)
        ENDIF
        RETURN
        END
C*******************************************************************************
      SUBROUTINE CV(MAX,KD,Q,EVAL) 
	USE MSIMSLMD
c	use numerical_libraries

	IMPLICIT NONE
	INTEGER KD,MAX
	COMPLEX*16 Q,EVAL(MAX),A(MAX,MAX),EVEC(MAX,MAX)
	INTEGER K,I,J
	IF(KD.EQ.1)THEN
	 K = 0
	 DO I = 1 , MAX
	  DO J = 1 , MAX
	   IF(I.EQ.J)THEN
	    A(I,J)= (2*K)**2
	    K = K + 1
	   ELSEIF(I.EQ.J+1.OR.I+1.EQ.J)THEN
	    A(I,J) = Q
	   ELSE
	    A(I,J)= (0.0D0,0.0D0)
         ENDIF
	  ENDDO
	 ENDDO
	 A(2,1) = 2.0D0*Q
	ELSEIF(KD.EQ.2)THEN
	 K = 0
	 DO I = 1 , MAX
	  DO J = 1 , MAX
	   IF(I.EQ.J)THEN
	    A(I,J)= (2*K+1)**2
	    K = K + 1
	   ELSEIF(I.EQ.J+1.OR.I+1.EQ.J)THEN
	    A(I,J) = Q
	   ELSE
	    A(I,J)= (0.0D0,0.0D0)
         ENDIF
	  ENDDO
	 ENDDO
	 A(1,1) = 1.0D0 + Q
	ELSEIF(KD.EQ.3)THEN
	 K = 0
	 DO I = 1 , MAX
	  DO J = 1 , MAX
	   IF(I.EQ.J)THEN
	    A(I,J)= (2*K+1)**2
	    K = K + 1
	   ELSEIF(I.EQ.J+1.OR.I+1.EQ.J)THEN
	    A(I,J) = Q
	   ELSE
	    A(I,J)= (0.0D0,0.0D0)
         ENDIF
	  ENDDO
	 ENDDO
	 A(1,1) = 1.0D0-Q	
      ELSEIF(KD.EQ.4)THEN
	 K = 0
	 DO I = 1 , MAX
	  DO J = 1 , MAX
	   IF(I.EQ.J)THEN
	    A(I,J)= (2*K+2)**2
	    K = K + 1
	   ELSEIF(I.EQ.J+1.OR.I+1.EQ.J)THEN
	    A(I,J) = Q
	   ELSE
	    A(I,J)= (0.0D0,0.0D0)
         ENDIF
	  ENDDO
	 ENDDO 
      ENDIF                            
      CALL DEVCCG(MAX, A, MAX, EVAL, EVEC, MAX)
	CALL SORT(EVAL,MAX)
	RETURN
      END
C****************************************************************************	  	  
	SUBROUTINE CHOLESKY(A,MAX,MAX1,N,M,X)
	IMPLICIT NONE
	INTEGER MAX,MAX1,M,N,I,II,J,K,NN
	COMPLEX*16 A(MAX,MAX1),X(MAX),SUM
	DO J = 2 , M
	 A(1,J) = A(1,J)/A(1,1)
	ENDDO
	DO I = 2 , N
	 DO II = I , N
	  SUM = (0.0D0,0.0D0)
	  DO K = 1 , I-1
	   SUM = SUM + A(II,K)*A(K,I)
	  ENDDO
	  A(II,I) = A(II,I) - SUM
	 ENDDO
	 DO J = I+1, M
	  SUM = (0.0D0,0.0D0)
	  DO K = 1 , I-1
	   SUM = SUM + A(I,K)*A(K,J)
	  ENDDO
	  A(I,J) = (A(I,J) - SUM)/A(I,I)
	 ENDDO
	ENDDO
	X(N) = A(N,M)
	DO NN = 1 , N-1
	 SUM = (0.0D0,0.0D0)
	 I = N - NN
	 DO J = I+1 , N
	  SUM = SUM + A(I,J)*X(J)
	 ENDDO
	 X(I) = A(I,M) - SUM
	ENDDO
	RETURN
	END
C******************************************************************************
