C
C PURPOSE: MAIN PROGRAM TO COMPUTE EMISSION COEFFS OF M1 CORONAL LINES
C
C INPUTS:
C
C OUTPUTS:
C
C COMMON:
C
C COMMENTS: AUGUST 9, 2000 P. JUDGE
C 
      PROGRAM PGJTEST
C
C  MAIN PROGRAM
C
C  CALCULATES THE EMISSION COEFFICIENTS OF FORBIDDEN CORONAL
C  LINES, INTEGRATED ALONG THE LOS, FOR A SET OF LOS'S IN AN 
C  ASSIGNED GRID OF CARTESIAN COORDINATES.
C
      INCLUDE 'PREC'
      INCLUDE 'PARAM'
      INCLUDE 'GRDPAR'
      INCLUDE 'CCONST'
      INCLUDE 'CBCNST'
      INCLUDE 'CATOM'
      INCLUDE 'CATMOS'
      INCLUDE 'CATMO2'  
      INCLUDE 'CSLINE'
      INCLUDE 'CORON'
      INCLUDE 'CINPUT'
      INCLUDE 'CINTS'
      INCLUDE 'CGRID'
      INCLUDE 'CLU'
      DIMENSION EDENS(MY)
      DIMENSION ETEMP(MZ)
      DIMENSION AA(MJTOT,MJTOT),BB(MJTOT)
      DIMENSION T0(0:3,0:2)
      LOGICAL DISK,IGNORE
      CALL CPTIME('TOTAL ',0,0,5)
C
      CALL START
C
C  LOOP OVER TEMPERATURE AND DENSITY
C
C      
      T1=5.5
      T2=6.5
      E1=6.
      E2=16.
C       RAD    RADIUS OF POINT RELATIVE TO SUN CENTER, UNITS OF SOLAR RADII
C       THET   AZIMUTHAL ANGLE, UNITS RADIANS, IN S FRAME
C       PH     POLAR ANGLE, UNITS RADIANS, IN S FRAME
       RAD=1.5
       THET=0.
       PH=0.5
       B=10     
       T=0.1    
       P=0.5    
       VEL=0.   
       TURBV=10.
C       write(LOUT,*) 'PGJTEST'
       alpha=.1
       bfield=10.
       thetab=0.1
       phib=0.5
       theta=0.3
       pgamma=0.2
       h=.5
       dop=1.e15
C     
C     LARMOR FREQUENCY 
C     
      ALARMOR=(1.0E-8*(0.25E0/PI)*(EE/EM)*BFIELD)/QNORM
C     
      VTURB=1.
      VEL=0.
C HYDROGEN POPULATIONS: ZERO IF NEUTRAL, TOTH FOR PROTONS
      DO I=1,5
         HD(I)=0.
      ENDDO
C
C  adopt geometric grid as the temperature abnd density grid       
C
       NGX=1
       NGY=10
       NGZ=21
       NTEMP=NGY
       NDENS=NGZ
       DO J=1,NTEMP
          ETEMP(J)=T1+ FLOAT((J-1))/FLOAT(NTEMP-1)* (T2-T1)
          ETEMP(J)=10.**ETEMP(J)
          TEMP=ETEMP(J)
          EQI=EQION(TEMP)
C          write(LOUT,*) 'TEMP', J, ALOG10(ETEMP(J)),EQI
          DO K=1,NDENS
             EDENS(K)=E1+ FLOAT((K-1))/FLOAT(NDENS-1)* (E2-E1)
             EDENS(K)=10.**EDENS(K)
             ED=EDENS(K)
             HDENS=0.8*ED
             TOTH=0.8*ED
             TOTN=10.**(ABND-12.)*TOTH*EQI
             HD(6)=TOTH
             DO KR=1,NLINE
                DO M=0,4
                   DO NY=1,NQ(KR)
                      EMISS(KR,M,NY)=ZERO
                   END DO
                END DO
             END DO
             CALL PROFIL
C     
C     LTE POPULATIONS AND COLLISIONAL RATES TO BE STORED
C     
             IF(ICOLL .NE. 0) THEN 
                CALL LTEPOP
                CALL COLCAL
             ENDIF
C     
C     BUILD AND SOLVE STATISTICAL EQUILIBRIUM EQUATIONS
C     
             CALL SE0_BUILD(NDIM)
             CALL SOLVE(NDIM,AA,BB)
             CALL SE0_COPY(NDIM,BB)
C     
C     SOLVE FOR EMERGENT STOKES PROFILES
C     
             CALL T0TENS(T0)
C     
             
             WTX=1.e6
             ED1=ED/1.e8
             DO KR=1,NLINE
                IF ((ALAMB(KR).GE.WLMIN).AND.(ALAMB(KR).LE.WLMAX)) THEN
                   CALL EMISSION(KR,T0)
                   DO NY=1,NQ(KR)
                      DO MS=0,4
                         EMISS(KR,MS,NY) = EMISS(KR,MS,NY)/ED1/ED1
                      END DO
                   END DO
                END IF
             END DO
             CALL TRAP(WTX)
             JJ=J
             KK=K
             CALL OUTS(JJ,KK)
          END DO
       END DO
       CALL CPTIME('TOTAL ',0,1,5)
       CALL STOP('NORMAL END')
       END
      

