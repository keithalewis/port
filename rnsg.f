      SUBROUTINE   RNSG(A, ALF, C, DA, IN, IV, L, L1, LA, LIV, LV,
     1                  N, NDA, P, V, Y)
C
C  ***  ITERATION DRIVER FOR SEPARABLE NONLINEAR LEAST SQUARES.
C
C  ***  PARAMETER DECLARATIONS  ***
C
      INTEGER L, L1, LA, LIV, LV, N, NDA, P
      INTEGER IN(2,NDA), IV(LIV)
C     DIMENSION UIPARM(*)
      REAL A(LA,L1), ALF(P), C(L), DA(LA,NDA), V(LV), Y(N)
C
C  ***  PURPOSE  ***
C
C GIVEN A SET OF N OBSERVATIONS Y(1)....Y(N) OF A DEPENDENT VARIABLE
C T(1)...T(N),   RNSG ATTEMPTS TO COMPUTE A LEAST SQUARES FIT
C TO A FUNCTION  ETA  (THE MODEL) WHICH IS A LINEAR COMBINATION
C
C                  L
C ETA(C,ALF,T) =  SUM C * PHI(ALF,T) +PHI   (ALF,T)
C                 J=1  J     J           L+1
C
C OF NONLINEAR FUNCTIONS PHI(J) DEPENDENT ON T AND ALF(1),...,ALF(P)
C (.E.G. A SUM OF EXPONENTIALS OR GAUSSIANS).  THAT IS, IT DETERMINES
C NONLINEAR PARAMETERS ALF WHICH MINIMIZE
C
C                   2    N                      2
C     NORM(RESIDUAL)  = SUM  (Y - ETA(C,ALF,T )).
C                       I=1    I             I
C
C THE (L+1)ST TERM IS OPTIONAL.
C
C
C  ***  PARAMETERS  ***
C
C      A (IN)  MATRIX PHI(ALF,T) OF THE MODEL.
C    ALF (I/O) NONLINEAR PARAMETERS.
C                 INPUT = INITIAL GUESS,
C                 OUTPUT = BEST ESTIMATE FOUND.
C      C (OUT) LINEAR PARAMETERS (ESTIMATED).
C     DA (IN)  DERIVATIVES OF COLUMNS OF A WITH RESPECT TO COMPONENTS
C                 OF ALF, AS SPECIFIED BY THE IN ARRAY...
C     IN (IN)  WHEN   RNSG IS CALLED WITH IV(1) = 2 OR -2, THEN FOR
C                 I = 1(1)NDA, COLUMN I OF DA IS THE PARTIAL
C                 DERIVATIVE WITH RESPECT TO ALF(IN(1,I)) OF COLUMN
C                 IN(2,I) OF A, UNLESS IV(1,I) IS NOT POSITIVE (IN
C                 WHICH CASE COLUMN I OF DA IS IGNORED.  IV(1) = -2
C                 MEANS THERE ARE MORE COLUMNS OF DA TO COME AND
C                   RNSG SHOULD RETURN FOR THEM.
C     IV (I/O) INTEGER PARAMETER AND SCRATCH VECTOR.    RNSG RETURNS
C                 WITH IV(1) = 1 WHEN IT WANTS A TO BE EVALUATED AT
C                 ALF AND WITH IV(1) = 2 WHEN IT WANTS DA TO BE
C                 EVALUATED AT ALF.  WHEN CALLED WITH IV(1) = -2
C                 (AFTER A RETURN WITH IV(1) = 2),   RNSG RETURNS
C                 WITH IV(1) = -2 TO GET MORE COLUMNS OF DA.
C      L (IN)  NUMBER OF LINEAR PARAMETERS TO BE ESTIMATED.
C     L1 (IN)  L+1 IF PHI(L+1) IS IN THE MODEL, L IF NOT.
C     LA (IN)  LEAD DIMENSION OF A.  MUST BE AT LEAST N.
C    LIV (IN)  LENGTH OF IV.  MUST BE AT LEAST 110 + L + P.
C     LV (IN)  LENGTH OF V.  MUST BE AT LEAST
C                 105 + 2*N + JLEN + L*(L+3)/2 + P*(2*P + 17),
C                 WHERE  JLEN = (L+P)*(N+L+P+1),  UNLESS NEITHER A
C                 COVARIANCE MATRIX NOR REGRESSION DIAGNOSTICS ARE
C                 REQUESTED, IN WHICH CASE  JLEN = N*P.
C      N (IN)  NUMBER OF OBSERVATIONS.
C    NDA (IN)  NUMBER OF COLUMNS IN DA AND IN.
C      P (IN)  NUMBER OF NONLINEAR PARAMETERS TO BE ESTIMATED.
C      V (I/O) FLOATING-POINT PARAMETER AND SCRATCH VECTOR.
C              IF A COVARIANCE ESTIMATE IS REQUESTED, IT IS FOR
C              (ALF,C) -- NONLINEAR PARAMETERS ORDERED FIRST,
C              FOLLOWED BY LINEAR PARAMETERS.
C      Y (IN)  RIGHT-HAND SIDE VECTOR.
C
C
C  ***  EXTERNAL SUBROUTINES  ***
C
      REAL  D7TPR,  L7SVX,  L7SVN,  R7MDC
      EXTERNAL  C7VFN, IVSET,  D7TPR, ITSUM,  L7ITV, L7SRT,  L7SVX,
     1          L7SVN,  N2CVP,  N2LRD,  N2RDP,   RN2G,  Q7APL, Q7RAD,
     2         Q7RFH,  R7MDC,  S7CPR, V2AXY, V7CPY, V7PRM,  V7SCL,
     3          V7SCP
C
C  C7VFN... FINISHES COVARIANCE COMPUTATION.
C IVSET.... SUPPLIES DEFAULT PARAMETER VALUES.
C  D7TPR... RETURNS INNER PRODUCT OF TWO VECTORS.
C ITSUM.... PRINTS ITERATION SUMMARY, INITIAL AND FINAL ALF.
C  L7ITV... APPLIES INVERSE-TRANSPOSE OF COMPACT LOWER TRIANG. MATRIX.
C L7SRT.... COMPUTES (PARTIAL) CHOLESKY FACTORIZATION.
C  L7SVX... ESTIMATES LARGEST SING. VALUE OF LOWER TRIANG. MATRIX.
C  L7SVN... ESTIMATES SMALLEST SING. VALUE OF LOWER TRIANG. MATRIX.
C  N2CVP... PRINTS COVARIANCE MATRIX.
C  N2LRD... COMPUTES COVARIANCE AND REGRESSION DIAGNOSTICS.
C  N2RDP... PRINTS REGRESSION DIAGNOSTICS.
C   RN2G... UNDERLYING NONLINEAR LEAST-SQUARES SOLVER.
C  Q7APL... APPLIES HOUSEHOLDER TRANSFORMS STORED BY Q7RFH.
C Q7RFH.... COMPUTES QR FACT. VIA HOUSEHOLDER TRANSFORMS WITH PIVOTING.
C Q7RAD.... QR FACT., NO PIVOTING.
C  R7MDC... RETURNS MACHINE-DEP. CONSTANTS.
C  S7CPR... PRINTS LINEAR PARAMETERS AT SOLUTION.
C V2AXY.... ADDS MULTIPLE OF ONE VECTOR TO ANOTHER.
C V7CPY.... COPIES ONE VECTOR TO ANOTHER.
C V7PRM.... PERMUTES A VECTOR.
C  V7SCL... SCALES AND COPIES ONE VECTOR TO ANOTHER.
C  V7SCP... SETS ALL COMPONENTS OF A VECTOR TO A SCALAR.
C
C  ***  LOCAL VARIABLES  ***
C
      LOGICAL NOCOV
      INTEGER AR1, CSAVE1, D1, DR1, DR1L, DRI, DRI1, FDH0, HSAVE, I, I1,
     1        IPIV1, IER, IV1, J1, JLEN, K, LH, LI, LL1O2, MD, N1, N2,
     2        NML, NRAN, PP, PP1, R1, R1L, RD1, TEMP1
      REAL SINGTL, T
      REAL MACHEP, NEGONE, SNGFAC, ZERO
C
C  ***  SUBSCRIPTS FOR IV AND V  ***
C
      INTEGER AR, CNVCOD, COVMAT, COVREQ, CSAVE, CVRQSV, D, FDH, H,
     1        IERS, IPIVS, IV1SAV, IVNEED, J, LMAT, MODE, NEXTIV, NEXTV,
     2        NFCALL, NFCOV, NFGCAL, NGCALL, NGCOV, PERM, R, RCOND,
     3        RDREQ, RDRQSV, REGD, REGD0, RESTOR, TOOBIG, VNEED
C
C  ***  IV SUBSCRIPT VALUES  ***
C
C/6
C     DATA AR/110/, CNVCOD/55/, COVMAT/26/, COVREQ/15/, CSAVE/105/,
C    1     CVRQSV/106/, D/27/, FDH/74/, H/56/, IERS/108/, IPIVS/109/,
C    2     IV1SAV/104/, IVNEED/3/, J/70/, LMAT/42/, MODE/35/,
C    3     NEXTIV/46/, NEXTV/47/, NFCALL/6/, NFCOV/52/, NFGCAL/7/,
C    4     NGCALL/30/, NGCOV/53/, PERM/58/, R/61/, RCOND/53/, RDREQ/57/,
C    5     RDRQSV/107/, REGD/67/, REGD0/82/, RESTOR/9/, TOOBIG/2/,
C    6     VNEED/4/
C/7
      PARAMETER (AR=110, CNVCOD=55, COVMAT=26, COVREQ=15, CSAVE=105,
     1           CVRQSV=106, D=27, FDH=74, H=56, IERS=108, IPIVS=109,
     2           IV1SAV=104, IVNEED=3, J=70, LMAT=42, MODE=35,
     3           NEXTIV=46, NEXTV=47, NFCALL=6, NFCOV=52, NFGCAL=7,
     4           NGCALL=30, NGCOV=53, PERM=58, R=61, RCOND=53, RDREQ=57,
     5           RDRQSV=107, REGD=67, REGD0=82, RESTOR=9, TOOBIG=2,
     6           VNEED=4)
C/
      DATA MACHEP/-1.E+0/, NEGONE/-1.E+0/, SNGFAC/1.E+2/, ZERO/0.E+0/
C
C++++++++++++++++++++++++++++++++  BODY  ++++++++++++++++++++++++++++++
C
C
      IF (IV(1) .EQ. 0) CALL IVSET(1, IV, LIV, LV, V)
      N1 = 1
      NML = N
      IV1 = IV(1)
      IF (IV1 .LE. 2) GO TO 20
C
C  ***  CHECK INPUT INTEGERS  ***
C
      IF (P .LE. 0) GO TO 370
      IF (L .LT. 0) GO TO 370
      IF (N .LE. L) GO TO 370
      IF (LA .LT. N) GO TO 370
      IF (IV1 .LT. 12) GO TO 20
      IF (IV1 .EQ. 14) GO TO 20
      IF (IV1 .EQ. 12) IV(1) = 13
C
C  ***  FRESH START -- COMPUTE STORAGE REQUIREMENTS  ***
C
      IF (IV(1) .GT. 16) GO TO 370
      LL1O2 = L*(L+1)/2
      JLEN = N*P
      I = L + P
      IF (IV(RDREQ) .GT. 0 .AND. IV(COVREQ) .NE. 0) JLEN = I*(N + I + 1)
      IF (IV(1) .NE. 13) GO TO 10
         IV(IVNEED) = IV(IVNEED) + L
         IV(VNEED) = IV(VNEED) + P + 2*N + JLEN + LL1O2 + L
 10   IF (IV(PERM) .LE. AR) IV(PERM) = AR + 1
      CALL   RN2G(V, V, IV, LIV, LV, N, N, N1, NML, P, V, V, V, ALF)
      IF (IV(1) .NE. 14) GO TO 999
C
C  ***  STORAGE ALLOCATION  ***
C
      IV(IPIVS) = IV(NEXTIV)
      IV(NEXTIV) = IV(NEXTIV) + L
      IV(D) = IV(NEXTV)
      IV(REGD0) = IV(D) + P
      IV(AR) = IV(REGD0) + N
      IV(CSAVE) = IV(AR) + LL1O2
      IV(J) = IV(CSAVE) + L
      IV(R) = IV(J) + JLEN
      IV(NEXTV) = IV(R) + N
      IV(IERS) = 0
      IF (IV1 .EQ. 13) GO TO 999
C
C  ***  SET POINTERS INTO IV AND V  ***
C
 20   AR1 = IV(AR)
      D1 = IV(D)
      DR1 = IV(J)
      DR1L = DR1 + L
      R1 = IV(R)
      R1L = R1 + L
      RD1 = IV(REGD0)
      CSAVE1 = IV(CSAVE)
      NML = N - L
      IF (IV1 .LE. 2) GO TO 50
C
C  ***  IF F.D. HESSIAN WILL BE NEEDED (FOR COVARIANCE OR REG.
C  ***  DIAGNOSTICS), HAVE   RN2G COMPUTE ONLY THE PART CORRESP.
C  ***  TO ALF WITH C FIXED...
C
      IF (L .LE. 0) GO TO 30
      IV(CVRQSV) = IV(COVREQ)
      IF (IABS(IV(COVREQ)) .GE. 3) IV(COVREQ) = 0
      IV(RDRQSV) = IV(RDREQ)
      IF (IV(RDREQ) .GT. 0) IV(RDREQ) = -1
C
 30   N2 = NML
      CALL   RN2G(V(D1), V(DR1L), IV, LIV, LV, NML, N, N1, N2, P,
     1            V(R1L), V(RD1), V, ALF)
      IF (IABS(IV(RESTOR)-2) .EQ. 1 .AND. L .GT. 0)
     1        CALL V7CPY(L, C, V(CSAVE1))
      IV1 = IV(1)
      IF (IV1-2) 40, 150, 230
C
C  ***  NEW FUNCTION VALUE (RESIDUAL) NEEDED  ***
C
 40   IV(IV1SAV) = IV(1)
      IV(1) = IABS(IV1)
      IF (IV(RESTOR) .EQ. 2 .AND. L .GT. 0) CALL V7CPY(L, V(CSAVE1), C)
      GO TO 999
C
C  ***  COMPUTE NEW RESIDUAL OR GRADIENT  ***
C
 50   IV(1) = IV(IV1SAV)
      MD = IV(MODE)
      IF (MD .LE. 0) GO TO 60
         NML = N
         DR1L = DR1
         R1L = R1
 60   IF (IV(TOOBIG) .NE. 0) GO TO 30
      IF (IABS(IV1) .EQ. 2) GO TO 170
C
C  ***  COMPUTE NEW RESIDUAL  ***
C
      IF (L1 .LE. L) CALL V7CPY(N, V(R1), Y)
      IF (L1 .GT. L) CALL V2AXY(N, V(R1), NEGONE, A(1,L1), Y)
      IF (MD .GT. 0) GO TO 120
      IER = 0
      IF (L .LE. 0) GO TO 110
      LL1O2 = L * (L + 1) / 2
      IPIV1 = IV(IPIVS)
      CALL Q7RFH(IER, IV(IPIV1), N, LA, 0, L, A, V(AR1), LL1O2, C)
C
C *** DETERMINE NUMERICAL RANK OF A ***
C
      IF (MACHEP .LE. ZERO) MACHEP =  R7MDC(3)
      SINGTL = SNGFAC * FLOAT(MAX0(L,N)) * MACHEP
      K = L
      IF (IER .NE. 0) K = IER - 1
 70   IF (K .LE. 0) GO TO 90
         T =  L7SVX(K, V(AR1), C, C)
         IF (T .GT. ZERO) T =  L7SVN(K, V(AR1), C, C) / T
         IF (T .GT. SINGTL) GO TO 80
         K = K - 1
         GO TO 70
C
C *** RECORD RANK IN IV(IERS)... IV(IERS) = 0 MEANS FULL RANK,
C *** IV(IERS) .GT. 0 MEANS RANK IV(IERS) - 1.
C
 80   IF (K .GE. L) GO TO 100
 90      IER = K + 1
         CALL  V7SCP(L-K, C(K+1), ZERO)
 100  IV(IERS) = IER
      IF (K .LE. 0) GO TO 110
C
C *** APPLY HOUSEHOLDER TRANSFORMATONS TO RESIDUALS...
C
      CALL  Q7APL(LA, N, K, A, V(R1), IER)
C
C *** COMPUTING C NOW MAY SAVE A FUNCTION EVALUATION AT
C *** THE LAST ITERATION.
C
      CALL  L7ITV(K, C, V(AR1), V(R1))
      CALL V7PRM(L, IV(IPIV1), C)
C
 110  IF(IV(1) .LT. 2) GO TO 220
      GO TO 999
C
C
C  ***  RESIDUAL COMPUTATION FOR F.D. HESSIAN  ***
C
 120  IF (L .LE. 0) GO TO 140
      DO 130 I = 1, L
 130     CALL V2AXY(N, V(R1), -C(I), A(1,I), V(R1))
 140  IF (IV(1) .GT. 0) GO TO 30
         IV(1) = 2
         GO TO 160
C
C  ***  NEW GRADIENT (JACOBIAN) NEEDED  ***
C
 150  IV(IV1SAV) = IV1
      IF (IV(NFGCAL) .NE. IV(NFCALL)) IV(1) = 1
 160  CALL  V7SCP(N*P, V(DR1), ZERO)
      GO TO 999
C
C  ***  COMPUTE NEW JACOBIAN  ***
C
 170  NOCOV = MD .LE. P .OR. IABS(IV(COVREQ)) .GE. 3
      FDH0 = DR1 + N*(P+L)
      IF (NDA .LE. 0) GO TO 370
      DO 180 I = 1, NDA
         I1 = IN(1,I) - 1
         IF (I1 .LT. 0) GO TO 180
         J1 = IN(2,I)
         K = DR1 + I1*N
         T = NEGONE
         IF (J1 .LE. L) T = -C(J1)
         CALL V2AXY(N, V(K), T, DA(1,I), V(K))
         IF (NOCOV) GO TO 180
         IF (J1 .GT. L) GO TO 180
C        ***  ADD IN (L,P) PORTION OF SECOND-ORDER PART OF HESSIAN
C        ***  FOR COVARIANCE OR REG. DIAG. COMPUTATIONS...
         J1 = J1 + P
         K = FDH0 + J1*(J1-1)/2 + I1
         V(K) = V(K) -  D7TPR(N, V(R1), DA(1,I))
 180     CONTINUE
      IF (IV1 .EQ. 2) GO TO 190
         IV(1) = IV1
         GO TO 999
 190  IF (L .LE. 0) GO TO 30
      IF (MD .GT. P) GO TO 240
      IF (MD .GT. 0) GO TO 30
      K = DR1
      IER = IV(IERS)
      NRAN = L
      IF (IER .GT. 0) NRAN = IER - 1
      IF (NRAN .LE. 0) GO TO 210
      DO 200 I = 1, P
         CALL  Q7APL(LA, N, NRAN, A, V(K), IER)
         K = K + N
 200     CONTINUE
 210  CALL V7CPY(L, V(CSAVE1), C)
 220  IF (IER .EQ. 0) GO TO 30
C
C     *** ADJUST SUBSCRIPTS DESCRIBING R AND DR...
C
         NRAN = IER - 1
         DR1L = DR1 + NRAN
         NML = N - NRAN
         R1L = R1 + NRAN
         GO TO 30
C
C  ***  CONVERGENCE OR LIMIT REACHED  ***
C
 230  IF (L .LE. 0) GO TO 350
      IV(COVREQ) = IV(CVRQSV)
      IV(RDREQ) = IV(RDRQSV)
      IF (IV(1) .GT. 6) GO TO 360
      IF (MOD(IV(RDREQ),4) .EQ. 0) GO TO 360
      IF (IV(FDH) .LE. 0 .AND. IABS(IV(COVREQ)) .LT. 3) GO TO 360
      IF (IV(REGD) .GT. 0) GO TO 360
      IF (IV(COVMAT) .GT. 0) GO TO 360
C
C  *** PREPARE TO FINISH COMPUTING COVARIANCE MATRIX AND REG. DIAG. ***
C
      PP = L + P
      I = 0
      IF (MOD(IV(RDREQ),4) .GE. 2) I = 1
      IF (MOD(IV(RDREQ),2) .EQ. 1 .AND. IABS(IV(COVREQ)) .EQ. 1) I = I+2
      IV(MODE) = PP + I
      I = DR1 + N*PP
      K = P * (P + 1) / 2
      I1 = IV(LMAT)
      CALL V7CPY(K, V(I), V(I1))
      I = I + K
      CALL  V7SCP(PP*(PP+1)/2 - K, V(I), ZERO)
      IV(NFCOV) = IV(NFCOV) + 1
      IV(NFCALL) = IV(NFCALL) + 1
      IV(NFGCAL) = IV(NFCALL)
      IV(CNVCOD) = IV(1)
      IV(IV1SAV) = -1
      IV(1) = 1
      IV(NGCALL) = IV(NGCALL) + 1
      IV(NGCOV) = IV(NGCOV) + 1
      GO TO 999
C
C  ***  FINISH COVARIANCE COMPUTATION  ***
C
 240  I = DR1 + N*P
      DO 250 I1 = 1, L
         CALL  V7SCL(N, V(I), NEGONE, A(1,I1))
         I = I + N
 250     CONTINUE
      PP = L + P
      HSAVE = IV(H)
      K = DR1 + N*PP
      LH = PP * (PP + 1) / 2
      IF (IABS(IV(COVREQ)) .LT. 3) GO TO 270
      I = IV(MODE) - 4
      IF (I .GE. PP) GO TO 260
      CALL  V7SCP(LH, V(K), ZERO)
      CALL Q7RAD(N, N, PP, V, .FALSE., V(K), V(DR1), V)
      IV(MODE) = I + 8
      IV(1) = 2
      IV(NGCALL) = IV(NGCALL) + 1
      IV(NGCOV) = IV(NGCOV) + 1
      GO TO 160
C
 260  IV(MODE) = I
      GO TO 300
C
 270  PP1 = P + 1
      DRI = DR1 + N*P
      LI = K + P*PP1/2
      DO 290 I = PP1, PP
         DRI1 = DR1
         DO 280 I1 = 1, I
            V(LI) = V(LI) +  D7TPR(N, V(DRI), V(DRI1))
            LI = LI + 1
            DRI1 = DRI1 + N
 280        CONTINUE
         DRI = DRI + N
 290     CONTINUE
      CALL L7SRT(PP1, PP, V(K), V(K), I)
      IF (I .NE. 0) GO TO 310
 300  TEMP1 = K + LH
      T =  L7SVN(PP, V(K), V(TEMP1), V(TEMP1))
      IF (T .LE. ZERO) GO TO 310
      T = T /  L7SVX(PP, V(K), V(TEMP1), V(TEMP1))
      V(RCOND) = T
      IF (T .GT.  R7MDC(4)) GO TO 320
 310     IV(REGD) = -1
         IV(COVMAT) = -1
         IV(FDH) = -1
         GO TO 340
 320  IV(H) = TEMP1
      IV(FDH) = IABS(HSAVE)
      IF (IV(MODE) - PP .LT. 2) GO TO 330
         I = IV(H)
         CALL  V7SCP(LH, V(I), ZERO)
 330  CALL  N2LRD(V(DR1), IV, V(K), LH, LIV, LV, N, N, PP, V(R1),
     1            V(RD1), V)
 340  CALL  C7VFN(IV, V(K), LH, LIV, LV, N, PP, V)
      IV(H) = HSAVE
C
 350  IF (IV(REGD) .EQ. 1) IV(REGD) = RD1
 360  IF (IV(1) .LE. 11) CALL  S7CPR(C, IV, L, LIV)
      IF (IV(1) .GT. 6) GO TO 999
         CALL  N2CVP(IV, LIV, LV, P+L, V)
         CALL  N2RDP(IV, LIV, LV, N, V(RD1), V)
         GO TO 999
C
 370  IV(1) = 66
      CALL ITSUM(V, V, IV, LIV, LV, P, V, ALF)
C
 999  RETURN
C
C  ***  LAST CARD OF   RNSG FOLLOWS  ***
      END
