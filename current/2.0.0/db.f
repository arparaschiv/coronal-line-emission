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
C                       LINE
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
      PARAMETER(MED=100)
      DIMENSION XED(MED)
      CHARACTER*6 FILEO
      CHARACTER*2 ZZ
      ZZ='DB'
      CALL CPTIME('DB ',0,0,5)
C
      CALL START

      CALL CSTRIP(LPIN,'DB.INPUT')
      READ(LDUMS,*) NGY, NED, NGX, NBPHI, NBTH
      READ(LDUMS,*) (XED(I),I=1,NED)
      READ(LDUMS,*) GYMIN, GYMAX, GXMIN, GXMAX, BPHIMIN, BPHIMAX,
     *              BTHMIN, BTHMAX
      READ(LDUMS,*) BFIELD
      CALL CLOSE(LDUMS)
C     
C
      CALL OPEN(LDBSC,'db.hdr',1,'NEW')
C      
C  GRIDS
C      
      GXSTEP=(GXMAX-GXMIN)/NGX
      GYSTEP=(GYMAX-GYMIN)/NGY
C
C ANGULAR GRID FOR B AT EACH "VOXEL"
C
      BTHSTEP=(BTHMAX-BTHMIN)/NBTH
      BPHISTEP=(BPHIMAX-BPHIMIN)/NBPHI
      GZ=0.
C
C  WRITE THE GRID INFORMATION TO FORMATTED FILE
C
      WRITE(LDBSC,*) NED,NGX,NBPHI,NBTH
      WRITE(LDBSC,*) (XED(I),I=1,NED)
      WRITE(LDBSC,*) GXMIN,GXMAX
      WRITE(LDBSC,*) BPHIMIN,BPHIMAX
      WRITE(LDBSC,*) BTHMIN, BTHMAX
      WRITE(LDBSC,*) BFIELD
      WRITE(LDBSC,*) NLINE
      CALL CLOSE(LDBSC)
      GZ=0.
C     
C     LARMOR FREQUENCY 
C     
      ALARMOR=(1.0E-8*(0.25E0/PI)*(EE/EM)*BFIELD)/QNORM
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
         WRITE(*,*)'DB DOING PROJECTED RADIUS Y=',J
     *        ,'/',NGY,'...', FILEO
         DO IED=1,NED
            DO K=1,NGX
               GX=FLOAT(K-1)*GXSTEP+GXMIN
               RB=SQRT(GX*GX+GY*GY+GZ*GZ)
               H=RB-1.0
C     
C     ALLEN 1973 PAGE 177 FIRST 2 TERMS PLUS HSE COMPT,
C     IGNORING 3RD NEAR SUN COMPONENT
C     
               BAUMBACH = 1.E8*(0.036/RB**1.5 + 1.55/RB**6.)
               ED=3.e8*EXP(- H/HSCALE) + BAUMBACH
               ED= ED * XED(IED)
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
               THETA=0.5*PI+ALPHA
               PGAMMA=BETA
C     
C     
C     LOOP OVER BPHI AND BTH, THESE ARE GIVEN IN THE CLE FRAME
C     
               DO L=1,NBPHI
                  BPHI=BPHIMIN+(L-1)*BPHISTEP
                  DO M=1,NBTH
                     BTH=BTHMIN+(M-1)*BTHSTEP
C     
C     HERE ARE MAGNETIC FIELD COMPONENTS IN THE CLE X,Y,Z FRAME
C     
                     XB=SIN(BTH)*COS(BPHI)
                     YB=SIN(BTH)*SIN(BPHI)
                     ZB=COS(BTH)
C     
C     rotation by beta      around x-axis
C     
C     YZB=SBETA*YB+CBETA*ZB
                     YYB= CBETA*YB -SBETA * ZB
C     
                     XXB= CALPHA*XB+SALPHA*(SBETA*YB+CBETA*ZB)
                     ZZB=-SALPHA*XB+CALPHA*(SBETA*YB+CBETA*ZB)
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
                     
                     count=count+1.
                  ENDDO
               ENDDO
            END DO
         END DO
         CALL CLOSE(LDB)
      END DO
      write(*,*)count, " calculations done"
      CALL CPTIME('DB ',0,1,5)
      CALL STOP('NORMAL END')
      END
      

