C
C PURPOSE: MAIN PROGRAM TO COMPUTE EMISSION COEFFS OF M1 CORONAL LINES
C          A SPATIAL MESH OF (NX . NZ) IS COMPUTED 
C          THEN AT EACH (X,Z) THE SE SYSTEM IS SOLVED AND THE EMISSION
C          COEFFICIENTS STORED FOR (NPHI . NTHETA) DIRECTIONS OF THE MAGNETIC FIELD.
C          X(I) = 1.0 + I * (4.-1.) / NX
C          Z(J) = -3.0 + J * 6. / NK
C     
C     THE ORDER OF LOOPS IS
C     NED
C        Y
C           X      
C              PHI
C                  THETA
C                      LOG B
C                          LINE
C
C INPUTS:
C
C OUTPUTS:
C
C COMMON:
C
C COMMENTS: JAN 9, 2020 P. JUDGE
C 
      PROGRAM DB
C
C  MAIN PROGRAM
C
C  CALCULATES THE EMISSION COEFFICIENTS OF FORBIDDEN CORONAL
C  LINES, INTEGRATED ALONG THE LOS, FOR A SET OF LOS'S IN AN 
C  ASSIGNED GRID OF CARTESIAN COORDINATES.
C
      INCLUDE 'PREC'
      INCLUDE 'PARAM'
      INCLUDE 'CGRID'
      INCLUDE 'CCONST'
      INCLUDE 'CBCNST'
      INCLUDE 'CORON'
      INCLUDE 'CSLINE'
      INCLUDE 'CATMOS'
      INCLUDE 'CATMO2'
      INCLUDE 'CATOM'
      INCLUDE 'CLU'
      INCLUDE 'CINPUT'
      CHARACTER*6 FILEO
      CHARACTER*2 ZZ
      ZZ='DB'
      CALL CPTIME('DB ',0,0,5)
C
      CALL START

      CALL CSTRIP(LPIN,'DB.INPUT')
      READ(LDUMS,*) NGY, NED, NGX, NBPHI, NBTH
      READ(LDUMS,*) EMIN, EMAX,GYMIN, GYMAX, GXMIN, GXMAX,
     *   BPHIMIN, BPHIMAX,BTHMIN, BTHMAX
      CALL CLOSE(LDUMS)
C     
C
      CALL OPEN(LDBSC,'db.hdr',1,'NEW')
C      
C  GRIDS
C     
c     PAR --- updated ngx and ngy to ngx-1 and ngy-1 to keep intervals 
C           symetric with respect to gxmax and gxmin 
C           IFs corrects for cases where only 1 position 
C           in height of depth is requested
      IF (NGX.GT.1) THEN
            GXSTEP=(GXMAX-GXMIN)/(NGX-1)
      ELSE  
            GXSTEP=0
      ENDIF
      IF (NGY.GT.1) THEN
            GYSTEP=(GYMAX-GYMIN)/(NGY-1)
      ELSE
            GYSTEP=0
      ENDIF
C
C ANGULAR GRID FOR B AT EACH "VOXEL"
C
      BPHISTEP=(BPHIMAX-BPHIMIN)/NBPHI
      BTHSTEP=(BTHMAX-BTHMIN)/NBTH


C FIXED AT 1G

      BFIELD=1.
C
C  WRITE THE GRID INFORMATION TO FORMATTED FILE
C
      WRITE(LDBSC,*) NED,NGX,NBPHI,NBTH
      WRITE(LDBSC,313) EMIN,EMAX
      WRITE(LDBSC,313) GXMIN,GXMAX
      WRITE(LDBSC,313) BPHIMIN,BPHIMAX
      WRITE(LDBSC,313) BTHMIN, BTHMAX
 313  FORMAT(2(1X,F10.6))
 314  FORMAT(1X,F10.2)
C 315  FORMAT(1x,A)
C
C  logarithmic multipliers      
C
      XEMULT=10.** ((EMAX-EMIN)/(NED-1))
C
C
      KOUNT=0
      DO KR=1,NLINE
         IF ((ALAMB(KR) .GE. WLMIN).AND.(ALAMB(KR) .LE. WLMAX)) THEN
            KOUNT=KOUNT+1
         ENDIF
      ENDDO
      WRITE(LDBSC,*) KOUNT
      DO KR=1,NLINE
         IF ((ALAMB(KR) .GE. WLMIN).AND.(ALAMB(KR) .LE. WLMAX)) THEN
            WRITE(LDBSC,314)CONVL(ALAMB(KR))
         ENDIF
      ENDDO
      WRITE(LDBSC,*) 0
      CALL CLOSE(LDBSC)
C
C  PLANE Z=0
C
      GZ=0.
C     
C     SOME THERMAL PARAMETERS
C     WRITE (ans, '(a, i3.3)') chr1, val1
C
      TEMP=2.E6
      HSCALE= 50.E8/RSUNCM
      count=0.
 973  format(A,I4.4)
      DO J=1,NGY
         GY=FLOAT(J-1)*GYSTEP+GYMIN
         KK = INT((GY-1.0)*1000.)
         write(FILEO,973) ZZ, KK
         CALL OPEN(LDB,FILEO // '.DAT',-1,'NEW')
         WRITE(*,*)'DB DOING PROJECTED RADIUS Y=  ',J
     *        ,'/',NGY,'...', FILEO
         
         XED=10.**EMIN/XEMULT
         DO IED=1,NED
C
C logarithmic      multiplier      
C
            XED= XEMULT*XED
            WRITE(*,*)'    DB DOING ELECTRON DENSITY=',IED
     *           ,'/',NED,' MULTIPLIER = ',XED
            DO K=1,NGX
               GX=FLOAT(K-1)*GXSTEP+GXMIN
               RB=SQRT(GX*GX+GY*GY+GZ*GZ)
               H=RB-1.0
               WRITE(*,*)'    DB DOING X=',GX
C     
C     ALLEN 1973 PAGE 177 FIRST 2 TERMS PLUS HSE COMPT,
C     IGNORING 3RD NEAR SUN COMPONENT
C     
               BAUMBACH = 1.E8*(0.036/RB**1.5 + 1.55/RB**6.)
               ED0=3.e8*EXP(- H/HSCALE) + BAUMBACH
               ED= ED0 * XED
               TOTH=.8*ED
               PTURBV=30./2./SQRT(2.)
               VTURB=PTURBV/QNORM
               PVEL=0.
               VEL=PVEL/QNORM
               TOTN=10.**(ABND-12.)*TOTH*EQION(TEMP)
               DO I=1,5
                  HD(I)=0.
               ENDDO
               HD(6)=TOTH
C
C SET UP RADIAL DEPENDENCE OF B
C
               B0=0.5/RB/RB + 0.004/(RB-.9970)**3
C     
C     COMPUTE NECESSARY ANGLES TAKEN FROM CORONA.F:
C     THAT DEPEND ONLY ON X,Y,Z  (NO B INFORMATION)
C     
               GX2=GX*GX
               GY2=GY*GY
               GZ2=GZ*GZ
C     
               ALPHA=-ATAN2(GX,SQRT(GY2+GZ2))
               SALPHA=SIN(ALPHA)
               CALPHA=COS(ALPHA)
C     
               BETA=ATAN2(GY,GZ)
               SBETA=SIN(BETA)
               CBETA=COS(BETA)
C     
C     FOR COMMON BLOCK CORONA, T0TENS ANGLES, GEOMETRY NOT B FIELD
C     
C     THETA = 0.5PI + ALPHA, DEFINING REFERENCE DIRECTION FOR LINEAR POLARIZATION
C
               THETA=0.5*PI+ALPHA
               PGAMMA=BETA
C     
C     LOOP OVER BPHI AND BTH
C     
               DO L=1,NBPHI
                  BPHI=BPHIMIN+(L-1)*BPHISTEP
                  DO M=1,NBTH
                     BTH=BTHMIN+(M-1)*BTHSTEP
C     
C     HERE ARE UNIT MAGNETIC FIELD COMPONENTS IN THE CLE X,Y,Z FRAME
C     
                     XB=SIN(BTH)*COS(BPHI)
                     YB=SIN(BTH)*SIN(BPHI)
                     ZB=COS(BTH)
C     
C     rotation by beta around x-axis
C     in z=0 plane beta = pi/2.
C
                     YYB= CBETA*YB -SBETA * ZB
                     YZB=SBETA*YB+CBETA*ZB
C     
                     XXB= CALPHA*XB+SALPHA*YZB
                     ZZB=-SALPHA*XB+CALPHA*YZB
C     
C     COMPUTE POLAR ANGLES OF B IN THE "S'-FRAME"
C     (VAN VLECK ANGLE = 54.7356103173)
C     THETAB AND PHIB RETURNED ARE:
C     These are varthetab and varphib in eq 42 44 of CJ99 & figure 5
C     (angles between LVS and B and planes containing B LVS and k LVS
C     
                     THETAB=ATAN2(SQRT(XXB*XXB+YYB*YYB),ZZB)	
                     PHIB=ATAN2(YYB,XXB)
                     CALL DBE
                     COUNT=COUNT+1.
C                    ALARMOR=(1.0E-8*(0.25E0/PI)*(EE/EM)*BFIELD)
C     *                    /QNORM
                  ENDDO
               ENDDO
            END DO
         END DO
         CALL CLOSE(LDB)
      END DO
      write(*,*)COUNT, " calculations done"
      CALL CPTIME('DB ',0,1,5)
      CALL STOP('NORMAL END')
      END
      

