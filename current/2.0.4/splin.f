C
C ****************************************************************
C
        FUNCTION SPLIN(T, X, Y, N, NNN)
C
C  SPLINE INTERPOLATION OF T IN A GIVEN TABLE OF POINTS(X(I),Y(I)).
C
C Y(I) IS A FUNCTION OF X(I) - X(N)>X(N-1)>....X(1) -OR- X(1)> X(2)
C > ....> X(N) - WHEN X IS NORMALISED X(MIN) = 0, XMAX=1,
C NORM(ABS(X(I+1)-X(I))) < 1.0E-8,(I=1,N-1). THE ARRAY X IS DIVIDED
C INTO OVERLAPPING INTERVALS. EACH INTERVAL CONSISTS OF NNN POINTS,
C THE OVERLAP IS INT(LOG(NNN)+3) POINTS. ( 3<=NNN<=100,N).
C EACH TIME AFTER A CALL WHEN T CHANGES OF INTERVAL OR AFTER A CALL
C WITH ANOTHER SET OF POINTS, A NEW CUBIC SPLINE OF NNN POINTS HAS
C TO BE COMPUTED, THEREFORE THE COMPUTING TIME DEPENDS STRONGLY OF
C THE SUCCESSION OF CALLS.
C FOR EXAMPLE: SUCCESSIVE CALLS, WHERE T IS RANDOMLY CHOSEN WILL COST
C A LOT OF COMPUTING TIME - (T OFTEN CHANGES OF INTERVAL) --->
C TAKE NNN THEN AS LARGE AS POSSIBLE. SUCCESSIVE CALLS WHERE T CHANGES
C SLIGHTLY COST LESS COMPUTING TIME, SO IT WILL BE CLEAR THAT YOUR
C CHOICE OF NNN DEPENDS HEAVILY UPON ITS WAY OF USEAGE.
C
C COMPUTING TIME:   FOR THE COMPUTATION OF SPLINE COEFFICIENTS
C        THIS ROUTINE NEEDS ABOUT 21*N MULTIPLICATIONS OR DIVISIONS
C        FOR THE INTERPOLATION BY MEANS OF THE SPLINE COEFFICIENTS
C        THE ROUTINE NEEDS 9 MULTIPLICATIONS OR DIVISIONS.
C
C EXTRAPOLATION: IF T IS NOT IN [X(1),X(N)] T IS LINEAR EXTRAPOLATED
C
C A PART OF THE ALGORITM FOR THE COMPUTATION OF A NATURAL SPLINE
C IS TAKEN FROM T.N.E. GREVILLE, THE LINEAR SYSTEM IS SOLVED BY
C THE METHOD OF SUCCESSIVE OVERRELAXATION. (ERROR = 1.0E-6)
C
C REF: T.N.E. GREVILLE, MATHEMATICS RESEARCH CENTER, U.S. ARMY,
C        UNIVERSITY OF WISCONSIN. MATHEMATICAL METHODS.
C
C THIS ROUTINE WAS IMPLEMENTED BY E.B.J. VAN DER ZALM, STERRENWACHT
C UTRECHT (MAY 7TH 1981). ADAPTED FOR VAX\11 BY PAUL KUIN AT OXFORD.
C WRITTEN IN REAL*4 BY P.G. JUDGE, OXFORD
C MODIFIED SO THAT IF THE ARGUMENT IS OUTSIDE THE RANGE OF X-VALUES
C THEN THE FIRST (OR LAST) Y-VALUE IS TAKEN
C
C INPUT:         T        ARGUMENT   (REAL)
C                X        ARRAY OF ARGUMENTS        (REAL)
C                Y        ARRAY OF FUNCTION VALUES (REAL)
C                N        LENGTH OF THE ARRAYS X,Y
C                        IF N = 2 ---> LINEAR INTERPOLATION
C                NNN        NUMBER OF POINTS OF THE CUBIC SPLINE
C                        IF NNN > N   ---> NNN = N
C                        IF NNN > 100 ---> NNN = 100
C                        IF NNN = 2 LINEAR INTERPOLATION ADOPTED.
C
C                        IF NNN IS ZERO NNN=7  (ALTERED BY P JUDGE)
C                       (USED TO BE IF NNN OMITTED )
C OUTPUT:        SPLIN INTERPOLATED VALUE AT T
C
      INCLUDE 'PREC'
        DIMENSION X(N), Y(N), S2(100), S3(100), DELY(100), H(100),
     * H2(100), DELSQY(100), C(100), B(100), XX(100), YY(100)
        LOGICAL SEARCH, FOUND
        CHARACTER*4 INIT
        DATA INIT/'NOTO'/
        SAVE
C
C  NN = NUMBER OF POINTS OF THE SPLINE INTERVAL
C
        NN = MIN(7,N)
        IF (NNN .EQ. 0) GO TO 10
C PGJ        IF (%LOC(NNN).EQ.0) GO TO 10
        NN = NNN
        NN = MIN(N,NNN,100)
C
C           IF T IN [X(K)-ERROR,X(K)+ERROR] T = X(K)
C
 10        ERROR = 1.0E-12
        IF (T.NE.0.) ERROR = ABS(1.0E-12*T)
        FOUND = SEARCH(X,T,1,N,K,ERROR)
        IF (.NOT.FOUND) GOTO 20
        SPLIN = Y(K)
        RETURN
 20        IF (N.NE.2.AND.NN.NE.2) GOTO 30
C
C           IF N=2 OR NNN=2 LINEAR INTERPOLATION IS ADOPTED
C
        SPLIN = (Y(K+1)-Y(K))*(T-X(K))/(X(K+1)-X(K)) + Y(K)
        RETURN
 30        IF (N.NE.1) GOTO 40
        SPLIN = Y(1)
        RETURN
C
C           INTERVAL OF THE SPLINE IS COMPUTED
C
 40        IOVER = LOG10(FLOAT(NN))*3
        INT = (K-IOVER)/(NN-2*IOVER) + 1
        I1 = (INT-1)*(NN-2*IOVER)
        I2 = I1 + NN -1
        IF (I1.GE.1) GOTO 50
        I2 = NN
        I1 = 1
 50        IF (I2.LE.N) GOTO 60
        I2 = N
        I1 = N - NN + 1
 60         ISHIFT = I1 - 1
        N1 = NN - 1
C
C           INITIALIZATION CHECK
C
        IF (INIT.NE.'OKAY') GOTO 80
C
C           INTERVAL CHECK
C
        DO 70 I=1,NN
        IF (XX(I).NE.X(I+ISHIFT)) GOTO 90
        IF (YY(I).NE.Y(I+ISHIFT)) GOTO 90
 70   CONTINUE
      GOTO 200
 80   INIT = 'OKAY'
C     
C     START OF THE COMPUTATION OF SPLINE COEFFICIENTS
C     
 90   RMIN = 1.0E37
      RMAX = -1.0E37
      DO 100 I=1,NN
        IF (Y(I+ISHIFT).LT.RMIN)  RMIN = Y(I+ISHIFT)
        IF (Y(I+ISHIFT).GT.RMAX)  RMAX = Y(I+ISHIFT)
 100  CONTINUE
C
C           COMPUTATION OF THE NORM FACTORS XNORM AND YNORM
C
        XNORM = 1/(X(NN+ISHIFT)-X(1+ISHIFT))
        YNORM = RMAX - RMIN
        IF (YNORM.EQ.0) YNORM = RMAX
        IF (YNORM.EQ.0) YNORM = 1
        YNORM = 1 / YNORM
        DO 110 I=1,N1
           H(I) = (X(I+1+ISHIFT)-X(I+ISHIFT))*XNORM
           DELY(I) = ((Y(I+1+ISHIFT)-Y(I+ISHIFT))/H(I))*YNORM
 110        CONTINUE
        DO 120 I=2,N1
           H2(I) = H(I-1) + H(I)
           B(I) = 0.5*H(I-1)/H2(I)
           DELSQY(I) = (DELY(I)-DELY(I-1))/H2(I)
           S2(I) = 2.*DELSQY(I)
           C(I) = 3.*DELSQY(I)
 120        CONTINUE
        S2(1) = 0
        S2(N1+1) = 0
C
C        SOLUTION OF THE LINEAR SYSTEM OF THE SPLINE COEFFICIENTS
C        BY SUCCESSIVE OVERRELAXATION.
C        CONSTANTS: EPS = ERROR CRITERION IN THE ITERATIVE SOLUTION.
C                   OMEGA = RELAXATION COEFFICIENT.
C
        EPSLN = 1.0E-8
        OMEGA = 1.0717968
 130    ETA = 0.
        DO 160 I=2,N1
           W = (C(I)-B(I)*S2(I-1)-(0.5-B(I))*S2(I+1)-S2(I))*OMEGA
           IF ((ABS(W)-ETA) .GT. 0.) ETA = ABS(W)
           S2(I) = S2(I) + W
 160    CONTINUE
        IF ((ETA-EPSLN) .GE. 0.) GOTO 130
        DO 180 I=1,N1
           S3(I) = (S2(I+1)-S2(I))/H(I)
 180    CONTINUE
C     
C     X AND Y STORED IN XX AND YY FOR INTERVAL CHECK
C     
        DO 190 I=1,NN
           XX(I) = X(I+ISHIFT)
           YY(I) = Y(I+ISHIFT)
 190    CONTINUE
 200    I = K - I1 + 1
        HT1 = (T-X(I+ISHIFT))*XNORM
        HT2 = (T-X(I+ISHIFT+1) )*XNORM
        IF ((T-X(1))*XNORM.GT.0) GOTO 210
C     
C     EXTRAPOLATION   T < X(1)
C     
        SPLIN = Y(1)
        RETURN
 210    IF ((T-X(N))*XNORM.GT.0) GOTO 220
C     
C     INTERPOLATION BY MEANS OF THE SPLINE COEFFICIENTS
C     
        SS2 = S2(I) + HT1*S3(I)
        PROD = HT1*HT2
        DELSQS = (S2(I)+S2(I+1)+SS2)/6.
        SPLIN = Y(I+ISHIFT) + (HT1*DELY(I)+PROD*DELSQS)/YNORM
        RETURN
C
C        EXTRAPOLATION   T > X(N)
C
 220        SPLIN =  Y(N)
        RETURN
        END
C
C **********************************************************************
C
      LOGICAL FUNCTION SEARCH(ARR,X,IA,IB,K,ERR)
C
C * SEARCHES THE POINT X WITHIN AN ERROR BOUND -ERR- IN ARRAY ARR      *
C * IF THAT POINT IS NOT FOUND, IT GIVES THE                           *
C * INTERVAL IN THE ARRAY, WHERE THAT POINT FITS.                      *
C * THE INTERVAL K IF X IN [ARR(K),ARR(K+1)), INTERVAL 1 IF X IN       *
C * (-INF,ARR(2)), INTERVAL N-1 IF X IN [ARR(N-1),INF)-(ARR INCREASING)*
C * INTERVAL 1 IF X IN (INF,ARR(2)),INTERVAL N-1 IF X IN [ARR(N-1),    *
C * -INF)-(ARR DECREASING)                                             *
C * IN THE WORST CASE A SEARCH COSTS 2*LOG2(N) CYCLES.                 *
C * THIS ROUTINE IS VERY ECONOMIC, WHEN IN A SEQUENCE OF SEARCHES,     *
C * THE DIFFERENCES BETWEEN THE POINTS WHICH ARE TO BE SEARCHED ARE    *
C * SMALL.                                                             *
C * FOR EXAMPLE: A SEQUENCE OF SEARCHES WHICH GOES ALONG THE ARRAY,    *
C * OR A SEQUENCE OF SEARCHES WHICH REMAINS IN A SMALL AREA OF         *
C * THE ARRAY.                                                         *
C *                                                                    *
C * INPUT:     - ARR         ARRAY WITH A NON DECREASING OR A NON      *
C *                          INCREASING FUCNTION.                      *
C *            - IA,IB       LOWER AND UPPER LIMIT OF THE ARRAY WITHIN *
C *                          IS SEARCHED.                              *
C *            - X           VALUE WHICH YOU WANT TO SEARCH            *
C *            - K           LAST INDEX WHICH WAS FOUND                *
C *                          SEARCH                                    *
C *            - ERR         ERROR BOUND                               *
C *                                                                    *
C * OUTPUT:    - K           INDEX OF THE POINT OR INTERVAL OF THE     *
C *                          ARRAY WHICH IS FOUND.                     *
C *            - SEARCH      TRUE - X IS FOUND WITHIN THE ERRORBOUND   *
C *                                  -ERR-                             *
C *                          FALSE - X IS NOT FOUND.                   *
C *                                                                    *
C * SEND IN BY E.V.D.ZALM, UTRECHT, STERREWACHT, MARCH 24TH 1981       *
C **********************************************************************
      INCLUDE 'PREC'
      DIMENSION ARR(IB)
C
      ONE=1.0
      SEARCH = .TRUE.
      IF (ARR(IB)-ARR(IA).LT.0.0) GOTO 100
      IF (X.LT.ARR(IA).OR.X.GT.ARR(IB)) GOTO 25
      IF (K.LT.IA.OR.K.GT.IB-1) GOTO 25
      L = K
      IP = SIGN(ONE,X-ARR(K))
      IF (IP.LT.0.) GOTO 20
10    L = MIN(K+IP,IB)
      IF (ARR(L).GE.X) GOTO 30
      IP = IP * 2
      K = L            
      GOTO 10
20    K = MAX(L+IP,IA)
      IF (ARR(K).LE.X) GOTO 30
      IP = IP * 2
      L = K
      GOTO 20
25    K = IA
      L = IB
30    IF ((L-K).LE.1) GOTO 50
      I = INT((K+L)/2.)
      IF (ARR(I).GE.X) GOTO 40
      K = I
      GOTO 30
40    L = I
      GOTO 30
100   IF (X.LT.ARR(IB).OR.X.GT.ARR(IA)) GOTO 250
      IF (K.LT.IA.OR.K.GT.IB-1) GOTO 250
      L = K
      IP = SIGN(ONE,ARR(K)-X)
      IF (IP.LT.0.) GOTO 200
110   L = MIN(K+IP,IB)
      IF (ARR(L).LE.X) GOTO 300
      IP = IP * 2
      K = L
      GOTO 110
200   K = MAX(L+IP,IA)
      IF (ARR(K).GE.X) GOTO 300
      IP = IP * 2
      L = K
      GOTO 200
250   K = IA
      L = IB
300   IF ((L-K).LE.1) GOTO 50
      I = INT((K+L)/2.)
      IF (ARR(I).LE.X) GOTO 400
      K = I
      GOTO 300
400   L = I
      GOTO 300
50    IF (ABS(X-ARR(K)).LE.ERR) RETURN
      IF (ABS(X-ARR(L)).GT.ERR) GOTO 60
      K = L
      RETURN
60    SEARCH = .FALSE.
      RETURN
      END
