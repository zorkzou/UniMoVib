      SUBROUTINE DPOTF2( UPLO, N, A, LDA, INFO )
c
c  -- LAPACK computational routine (version 3.4.2) --
c  -- LAPACK is a software package provided by Univ. of Tennessee,    --
c  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
c     September 2012
c
c     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
c     ..
c     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
c     ..
c
c  =====================================================================
c
c     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
c     ..
c     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            J
      DOUBLE PRECISION   AJJ
c     ..
c     .. External Functions ..
      LOGICAL            LSAME, DISNAN
      DOUBLE PRECISION   DDOT
      EXTERNAL           LSAME, DDOT, DISNAN
c     ..
c     .. External Subroutines ..
      EXTERNAL           DGEMV, DSCAL, XERBLA
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          MAX, SQRT
c     ..
c     .. Executable Statements ..
c
c     Test the input parameters.
c
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DPOTF2', -INFO )
         RETURN
      END IF
c
c     Quick return if possible
c
      IF( N.EQ.0 )
     $   RETURN
c
      IF( UPPER ) THEN
c
c        Compute the Cholesky factorization A = U**T *U.
c
         DO 10 J = 1, N
c
c           Compute U(J,J) and test for non-positive-definiteness.
c
            AJJ = A( J, J ) - DDOT( J-1, A( 1, J ), 1, A( 1, J ), 1 )
            IF( AJJ.LE.ZERO.OR.DISNAN( AJJ ) ) THEN
               A( J, J ) = AJJ
               GO TO 30
            END IF
            AJJ = SQRT( AJJ )
            A( J, J ) = AJJ
c
c           Compute elements J+1:N of row J.
c
            IF( J.LT.N ) THEN
               CALL DGEMV( 'Transpose', J-1, N-J, -ONE, A( 1, J+1 ),
     $                     LDA, A( 1, J ), 1, ONE, A( J, J+1 ), LDA )
               CALL DSCAL( N-J, ONE / AJJ, A( J, J+1 ), LDA )
            END IF
   10    CONTINUE
      ELSE
c
c        Compute the Cholesky factorization A = L*L**T.
c
         DO 20 J = 1, N
c
c           Compute L(J,J) and test for non-positive-definiteness.
c
            AJJ = A( J, J ) - DDOT( J-1, A( J, 1 ), LDA, A( J, 1 ),
     $            LDA )
            IF( AJJ.LE.ZERO.OR.DISNAN( AJJ ) ) THEN
               A( J, J ) = AJJ
               GO TO 30
            END IF
            AJJ = SQRT( AJJ )
            A( J, J ) = AJJ
c
c           Compute elements J+1:N of column J.
c
            IF( J.LT.N ) THEN
               CALL DGEMV( 'No transpose', N-J, J-1, -ONE, A( J+1, 1 ),
     $                     LDA, A( J, 1 ), LDA, ONE, A( J+1, J ), 1 )
               CALL DSCAL( N-J, ONE / AJJ, A( J+1, J ), 1 )
            END IF
   20    CONTINUE
      END IF
      GO TO 40
c
   30 CONTINUE
      INFO = J
c
   40 CONTINUE
      RETURN
c
c     End of DPOTF2
c
      END

      SUBROUTINE DPOTRF( UPLO, N, A, LDA, INFO )
c
c  -- LAPACK computational routine (version 3.4.0) --
c  -- LAPACK is a software package provided by Univ. of Tennessee,    --
c  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
c     November 2011
c
c     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
c     ..
c     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
c     ..
c
c  =====================================================================
c
c     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
c     ..
c     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            J, JB, NB
c     ..
c     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
c     ..
c     .. External Subroutines ..
      EXTERNAL           DGEMM, DPOTF2, DSYRK, DTRSM, XERBLA
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
c     ..
c     .. Executable Statements ..
c
c     Test the input parameters.
c
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DPOTRF', -INFO )
         RETURN
      END IF
c
c     Quick return if possible
c
      IF( N.EQ.0 )
     $   RETURN
c
c     Determine the block size for this environment.
c
      NB = ILAENV( 1, 'DPOTRF', UPLO, N, -1, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.N ) THEN
c
c        Use unblocked code.
c
         CALL DPOTF2( UPLO, N, A, LDA, INFO )
      ELSE
c
c        Use blocked code.
c
         IF( UPPER ) THEN
c
c           Compute the Cholesky factorization A = U**T*U.
c
            DO 10 J = 1, N, NB
c
c              Update and factorize the current diagonal block and test
c              for non-positive-definiteness.
c
               JB = MIN( NB, N-J+1 )
               CALL DSYRK( 'Upper', 'Transpose', JB, J-1, -ONE,
     $                     A( 1, J ), LDA, ONE, A( J, J ), LDA )
               CALL DPOTF2( 'Upper', JB, A( J, J ), LDA, INFO )
               IF( INFO.NE.0 )
     $            GO TO 30
               IF( J+JB.LE.N ) THEN
c
c                 Compute the current block row.
c
                  CALL DGEMM( 'Transpose', 'No transpose', JB, N-J-JB+1,
     $                        J-1, -ONE, A( 1, J ), LDA, A( 1, J+JB ),
     $                        LDA, ONE, A( J, J+JB ), LDA )
                  CALL DTRSM( 'Left', 'Upper', 'Transpose', 'Non-unit',
     $                        JB, N-J-JB+1, ONE, A( J, J ), LDA,
     $                        A( J, J+JB ), LDA )
               END IF
   10       CONTINUE
c
         ELSE
c
c           Compute the Cholesky factorization A = L*L**T.
c
            DO 20 J = 1, N, NB
c
c              Update and factorize the current diagonal block and test
c              for non-positive-definiteness.
c
               JB = MIN( NB, N-J+1 )
               CALL DSYRK( 'Lower', 'No transpose', JB, J-1, -ONE,
     $                     A( J, 1 ), LDA, ONE, A( J, J ), LDA )
               CALL DPOTF2( 'Lower', JB, A( J, J ), LDA, INFO )
               IF( INFO.NE.0 )
     $            GO TO 30
               IF( J+JB.LE.N ) THEN
c
c                 Compute the current block column.
c
                  CALL DGEMM( 'No transpose', 'Transpose', N-J-JB+1, JB,
     $                        J-1, -ONE, A( J+JB, 1 ), LDA, A( J, 1 ),
     $                        LDA, ONE, A( J+JB, J ), LDA )
                  CALL DTRSM( 'Right', 'Lower', 'Transpose', 'Non-unit',
     $                        N-J-JB+1, JB, ONE, A( J, J ), LDA,
     $                        A( J+JB, J ), LDA )
               END IF
   20       CONTINUE
         END IF
      END IF
      GO TO 40
c
   30 CONTINUE
      INFO = INFO + J - 1
c
   40 CONTINUE
      RETURN
c
c     End of DPOTRF
c
      END

      SUBROUTINE DSYGS2( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )
c
c  -- LAPACK computational routine (version 3.4.2) --
c  -- LAPACK is a software package provided by Univ. of Tennessee,    --
c  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
c     September 2012
c
c     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, ITYPE, LDA, LDB, N
c     ..
c     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
c     ..
c
c  =====================================================================
c
c     .. Parameters ..
      DOUBLE PRECISION   ONE, HALF
      PARAMETER          ( ONE = 1.0D0, HALF = 0.5D0 )
c     ..
c     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            K
      DOUBLE PRECISION   AKK, BKK, CT
c     ..
c     .. External Subroutines ..
      EXTERNAL           DAXPY, DSCAL, DSYR2, DTRMV, DTRSV, XERBLA
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          MAX
c     ..
c     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
c     ..
c     .. Executable Statements ..
c
c     Test the input parameters.
c
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( ITYPE.LT.1 .OR. ITYPE.GT.3 ) THEN
         INFO = -1
      ELSE IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSYGS2', -INFO )
         RETURN
      END IF
c
      IF( ITYPE.EQ.1 ) THEN
         IF( UPPER ) THEN
c
c           Compute inv(U**T)*A*inv(U)
c
            DO 10 K = 1, N
c
c              Update the upper triangle of A(k:n,k:n)
c
               AKK = A( K, K )
               BKK = B( K, K )
               AKK = AKK / BKK**2
               A( K, K ) = AKK
               IF( K.LT.N ) THEN
                  CALL DSCAL( N-K, ONE / BKK, A( K, K+1 ), LDA )
                  CT = -HALF*AKK
                  CALL DAXPY( N-K, CT, B( K, K+1 ), LDB, A( K, K+1 ),
     $                        LDA )
                  CALL DSYR2( UPLO, N-K, -ONE, A( K, K+1 ), LDA,
     $                        B( K, K+1 ), LDB, A( K+1, K+1 ), LDA )
                  CALL DAXPY( N-K, CT, B( K, K+1 ), LDB, A( K, K+1 ),
     $                        LDA )
                  CALL DTRSV( UPLO, 'Transpose', 'Non-unit', N-K,
     $                        B( K+1, K+1 ), LDB, A( K, K+1 ), LDA )
               END IF
   10       CONTINUE
         ELSE
c
c           Compute inv(L)*A*inv(L**T)
c
            DO 20 K = 1, N
c
c              Update the lower triangle of A(k:n,k:n)
c
               AKK = A( K, K )
               BKK = B( K, K )
               AKK = AKK / BKK**2
               A( K, K ) = AKK
               IF( K.LT.N ) THEN
                  CALL DSCAL( N-K, ONE / BKK, A( K+1, K ), 1 )
                  CT = -HALF*AKK
                  CALL DAXPY( N-K, CT, B( K+1, K ), 1, A( K+1, K ), 1 )
                  CALL DSYR2( UPLO, N-K, -ONE, A( K+1, K ), 1,
     $                        B( K+1, K ), 1, A( K+1, K+1 ), LDA )
                  CALL DAXPY( N-K, CT, B( K+1, K ), 1, A( K+1, K ), 1 )
                  CALL DTRSV( UPLO, 'No transpose', 'Non-unit', N-K,
     $                        B( K+1, K+1 ), LDB, A( K+1, K ), 1 )
               END IF
   20       CONTINUE
         END IF
      ELSE
         IF( UPPER ) THEN
c
c           Compute U*A*U**T
c
            DO 30 K = 1, N
c
c              Update the upper triangle of A(1:k,1:k)
c
               AKK = A( K, K )
               BKK = B( K, K )
               CALL DTRMV( UPLO, 'No transpose', 'Non-unit', K-1, B,
     $                     LDB, A( 1, K ), 1 )
               CT = HALF*AKK
               CALL DAXPY( K-1, CT, B( 1, K ), 1, A( 1, K ), 1 )
               CALL DSYR2( UPLO, K-1, ONE, A( 1, K ), 1, B( 1, K ), 1,
     $                     A, LDA )
               CALL DAXPY( K-1, CT, B( 1, K ), 1, A( 1, K ), 1 )
               CALL DSCAL( K-1, BKK, A( 1, K ), 1 )
               A( K, K ) = AKK*BKK**2
   30       CONTINUE
         ELSE
c
c           Compute L**T *A*L
c
            DO 40 K = 1, N
c
c              Update the lower triangle of A(1:k,1:k)
c
               AKK = A( K, K )
               BKK = B( K, K )
               CALL DTRMV( UPLO, 'Transpose', 'Non-unit', K-1, B, LDB,
     $                     A( K, 1 ), LDA )
               CT = HALF*AKK
               CALL DAXPY( K-1, CT, B( K, 1 ), LDB, A( K, 1 ), LDA )
               CALL DSYR2( UPLO, K-1, ONE, A( K, 1 ), LDA, B( K, 1 ),
     $                     LDB, A, LDA )
               CALL DAXPY( K-1, CT, B( K, 1 ), LDB, A( K, 1 ), LDA )
               CALL DSCAL( K-1, BKK, A( K, 1 ), LDA )
               A( K, K ) = AKK*BKK**2
   40       CONTINUE
         END IF
      END IF
      RETURN
c
c     End of DSYGS2
c
      END

      SUBROUTINE DSYGST( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )
c
c  -- LAPACK computational routine (version 3.4.0) --
c  -- LAPACK is a software package provided by Univ. of Tennessee,    --
c  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
c     November 2011
c
c     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, ITYPE, LDA, LDB, N
c     ..
c     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
c     ..
c
c  =====================================================================
c
c     .. Parameters ..
      DOUBLE PRECISION   ONE, HALF
      PARAMETER          ( ONE = 1.0D0, HALF = 0.5D0 )
c     ..
c     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            K, KB, NB
c     ..
c     .. External Subroutines ..
      EXTERNAL           DSYGS2, DSYMM, DSYR2K, DTRMM, DTRSM, XERBLA
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
c     ..
c     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
c     ..
c     .. Executable Statements ..
c
c     Test the input parameters.
c
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( ITYPE.LT.1 .OR. ITYPE.GT.3 ) THEN
         INFO = -1
      ELSE IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSYGST', -INFO )
         RETURN
      END IF
c
c     Quick return if possible
c
      IF( N.EQ.0 )
     $   RETURN
c
c     Determine the block size for this environment.
c
      NB = ILAENV( 1, 'DSYGST', UPLO, N, -1, -1, -1 )
c
      IF( NB.LE.1 .OR. NB.GE.N ) THEN
c
c        Use unblocked code
c
         CALL DSYGS2( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )
      ELSE
c
c        Use blocked code
c
         IF( ITYPE.EQ.1 ) THEN
            IF( UPPER ) THEN
c
c              Compute inv(U**T)*A*inv(U)
c
               DO 10 K = 1, N, NB
                  KB = MIN( N-K+1, NB )
c
c                 Update the upper triangle of A(k:n,k:n)
c
                  CALL DSYGS2( ITYPE, UPLO, KB, A( K, K ), LDA,
     $                         B( K, K ), LDB, INFO )
                  IF( K+KB.LE.N ) THEN
                     CALL DTRSM( 'Left', UPLO, 'Transpose', 'Non-unit',
     $                           KB, N-K-KB+1, ONE, B( K, K ), LDB,
     $                           A( K, K+KB ), LDA )
                     CALL DSYMM( 'Left', UPLO, KB, N-K-KB+1, -HALF,
     $                           A( K, K ), LDA, B( K, K+KB ), LDB, ONE,
     $                           A( K, K+KB ), LDA )
                     CALL DSYR2K( UPLO, 'Transpose', N-K-KB+1, KB, -ONE,
     $                            A( K, K+KB ), LDA, B( K, K+KB ), LDB,
     $                            ONE, A( K+KB, K+KB ), LDA )
                     CALL DSYMM( 'Left', UPLO, KB, N-K-KB+1, -HALF,
     $                           A( K, K ), LDA, B( K, K+KB ), LDB, ONE,
     $                           A( K, K+KB ), LDA )
                     CALL DTRSM( 'Right', UPLO, 'No transpose',
     $                           'Non-unit', KB, N-K-KB+1, ONE,
     $                           B( K+KB, K+KB ), LDB, A( K, K+KB ),
     $                           LDA )
                  END IF
   10          CONTINUE
            ELSE
c
c              Compute inv(L)*A*inv(L**T)
c
               DO 20 K = 1, N, NB
                  KB = MIN( N-K+1, NB )
c
c                 Update the lower triangle of A(k:n,k:n)
c
                  CALL DSYGS2( ITYPE, UPLO, KB, A( K, K ), LDA,
     $                         B( K, K ), LDB, INFO )
                  IF( K+KB.LE.N ) THEN
                     CALL DTRSM( 'Right', UPLO, 'Transpose', 'Non-unit',
     $                           N-K-KB+1, KB, ONE, B( K, K ), LDB,
     $                           A( K+KB, K ), LDA )
                     CALL DSYMM( 'Right', UPLO, N-K-KB+1, KB, -HALF,
     $                           A( K, K ), LDA, B( K+KB, K ), LDB, ONE,
     $                           A( K+KB, K ), LDA )
                     CALL DSYR2K( UPLO, 'No transpose', N-K-KB+1, KB,
     $                            -ONE, A( K+KB, K ), LDA, B( K+KB, K ),
     $                            LDB, ONE, A( K+KB, K+KB ), LDA )
                     CALL DSYMM( 'Right', UPLO, N-K-KB+1, KB, -HALF,
     $                           A( K, K ), LDA, B( K+KB, K ), LDB, ONE,
     $                           A( K+KB, K ), LDA )
                     CALL DTRSM( 'Left', UPLO, 'No transpose',
     $                           'Non-unit', N-K-KB+1, KB, ONE,
     $                           B( K+KB, K+KB ), LDB, A( K+KB, K ),
     $                           LDA )
                  END IF
   20          CONTINUE
            END IF
         ELSE
            IF( UPPER ) THEN
c
c              Compute U*A*U**T
c
               DO 30 K = 1, N, NB
                  KB = MIN( N-K+1, NB )
c
c                 Update the upper triangle of A(1:k+kb-1,1:k+kb-1)
c
                  CALL DTRMM( 'Left', UPLO, 'No transpose', 'Non-unit',
     $                        K-1, KB, ONE, B, LDB, A( 1, K ), LDA )
                  CALL DSYMM( 'Right', UPLO, K-1, KB, HALF, A( K, K ),
     $                        LDA, B( 1, K ), LDB, ONE, A( 1, K ), LDA )
                  CALL DSYR2K( UPLO, 'No transpose', K-1, KB, ONE,
     $                         A( 1, K ), LDA, B( 1, K ), LDB, ONE, A,
     $                         LDA )
                  CALL DSYMM( 'Right', UPLO, K-1, KB, HALF, A( K, K ),
     $                        LDA, B( 1, K ), LDB, ONE, A( 1, K ), LDA )
                  CALL DTRMM( 'Right', UPLO, 'Transpose', 'Non-unit',
     $                        K-1, KB, ONE, B( K, K ), LDB, A( 1, K ),
     $                        LDA )
                  CALL DSYGS2( ITYPE, UPLO, KB, A( K, K ), LDA,
     $                         B( K, K ), LDB, INFO )
   30          CONTINUE
            ELSE
c
c              Compute L**T*A*L
c
               DO 40 K = 1, N, NB
                  KB = MIN( N-K+1, NB )
c
c                 Update the lower triangle of A(1:k+kb-1,1:k+kb-1)
c
                  CALL DTRMM( 'Right', UPLO, 'No transpose', 'Non-unit',
     $                        KB, K-1, ONE, B, LDB, A( K, 1 ), LDA )
                  CALL DSYMM( 'Left', UPLO, KB, K-1, HALF, A( K, K ),
     $                        LDA, B( K, 1 ), LDB, ONE, A( K, 1 ), LDA )
                  CALL DSYR2K( UPLO, 'Transpose', K-1, KB, ONE,
     $                         A( K, 1 ), LDA, B( K, 1 ), LDB, ONE, A,
     $                         LDA )
                  CALL DSYMM( 'Left', UPLO, KB, K-1, HALF, A( K, K ),
     $                        LDA, B( K, 1 ), LDB, ONE, A( K, 1 ), LDA )
                  CALL DTRMM( 'Left', UPLO, 'Transpose', 'Non-unit', KB,
     $                        K-1, ONE, B( K, K ), LDB, A( K, 1 ), LDA )
                  CALL DSYGS2( ITYPE, UPLO, KB, A( K, K ), LDA,
     $                         B( K, K ), LDB, INFO )
   40          CONTINUE
            END IF
         END IF
      END IF
      RETURN
c
c     End of DSYGST
c
      END

      SUBROUTINE DSYGV( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, WORK,
     $                  LWORK, INFO )
c
c  -- LAPACK driver routine (version 3.4.0) --
c  -- LAPACK is a software package provided by Univ. of Tennessee,    --
c  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
c     November 2011
c
c     .. Scalar Arguments ..
      CHARACTER          JOBZ, UPLO
      INTEGER            INFO, ITYPE, LDA, LDB, LWORK, N
c     ..
c     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), W( * ), WORK( * )
c     ..
c
c  =====================================================================
c
c     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
c     ..
c     .. Local Scalars ..
      LOGICAL            LQUERY, UPPER, WANTZ
      CHARACTER          TRANS
      INTEGER            LWKMIN, LWKOPT, NB, NEIG
c     ..
c     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
c     ..
c     .. External Subroutines ..
      EXTERNAL           DPOTRF, DSYEV, DSYGST, DTRMM, DTRSM, XERBLA
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          MAX
c     ..
c     .. Executable Statements ..
c
c     Test the input parameters.
c
      WANTZ = LSAME( JOBZ, 'V' )
      UPPER = LSAME( UPLO, 'U' )
      LQUERY = ( LWORK.EQ.-1 )
c
      INFO = 0
      IF( ITYPE.LT.1 .OR. ITYPE.GT.3 ) THEN
         INFO = -1
      ELSE IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
         INFO = -2
      ELSE IF( .NOT.( UPPER .OR. LSAME( UPLO, 'L' ) ) ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
c
      IF( INFO.EQ.0 ) THEN
         LWKMIN = MAX( 1, 3*N - 1 )
         NB = ILAENV( 1, 'DSYTRD', UPLO, N, -1, -1, -1 )
         LWKOPT = MAX( LWKMIN, ( NB + 2 )*N )
         WORK( 1 ) = LWKOPT
c
         IF( LWORK.LT.LWKMIN .AND. .NOT.LQUERY ) THEN
            INFO = -11
         END IF
      END IF
c
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSYGV ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
c
c     Quick return if possible
c
      IF( N.EQ.0 )
     $   RETURN
c
c     Form a Cholesky factorization of B.
c
      CALL DPOTRF( UPLO, N, B, LDB, INFO )
      IF( INFO.NE.0 ) THEN
         INFO = N + INFO
         RETURN
      END IF
c
c     Transform problem to standard eigenvalue problem and solve.
c
      CALL DSYGST( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )
      CALL DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
c
      IF( WANTZ ) THEN
c
c        Backtransform eigenvectors to the original problem.
c
         NEIG = N
         IF( INFO.GT.0 )
     $      NEIG = INFO - 1
         IF( ITYPE.EQ.1 .OR. ITYPE.EQ.2 ) THEN
c
c           For A*x=(lambda)*B*x and A*B*x=(lambda)*x;
c           backtransform eigenvectors: x = inv(L)**T*y or inv(U)*y
c
            IF( UPPER ) THEN
               TRANS = 'N'
            ELSE
               TRANS = 'T'
            END IF
c
            CALL DTRSM( 'Left', UPLO, TRANS, 'Non-unit', N, NEIG, ONE,
     $                  B, LDB, A, LDA )
c
         ELSE IF( ITYPE.EQ.3 ) THEN
c
c           For B*A*x=(lambda)*x;
c           backtransform eigenvectors: x = L*y or U**T*y
c
            IF( UPPER ) THEN
               TRANS = 'T'
            ELSE
               TRANS = 'N'
            END IF
c
            CALL DTRMM( 'Left', UPLO, TRANS, 'Non-unit', N, NEIG, ONE,
     $                  B, LDB, A, LDA )
         END IF
      END IF
c
      WORK( 1 ) = LWKOPT
      RETURN
c
c     End of DSYGV
c
      END

      SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
c
c  -- LAPACK driver routine (version 3.4.0) --
c  -- LAPACK is a software package provided by Univ. of Tennessee,    --
c  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
c     November 2011
c
c     .. Scalar Arguments ..
      CHARACTER          JOBZ, UPLO
      INTEGER            INFO, LDA, LWORK, N
c     ..
c     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), W( * ), WORK( * )
c     ..
c
c  =====================================================================
c
c     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
c     ..
c     .. Local Scalars ..
      LOGICAL            LOWER, LQUERY, WANTZ
      INTEGER            IINFO, IMAX, INDE, INDTAU, INDWRK, ISCALE,
     $                   LLWORK, LWKOPT, NB
      DOUBLE PRECISION   ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA,
     $                   SMLNUM
c     ..
c     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      DOUBLE PRECISION   DLAMCH, DLANSY
      EXTERNAL           LSAME, ILAENV, DLAMCH, DLANSY
c     ..
c     .. External Subroutines ..
      EXTERNAL           DLASCL, DORGTR, DSCAL, DSTEQR, DSTERF, DSYTRD,
     $                   XERBLA
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          MAX, SQRT
c     ..
c     .. Executable Statements ..
c
c     Test the input parameters.
c
      WANTZ = LSAME( JOBZ, 'V' )
      LOWER = LSAME( UPLO, 'L' )
      LQUERY = ( LWORK.EQ.-1 )
c
      INFO = 0
      IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( LOWER .OR. LSAME( UPLO, 'U' ) ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      END IF
c
      IF( INFO.EQ.0 ) THEN
         NB = ILAENV( 1, 'DSYTRD', UPLO, N, -1, -1, -1 )
         LWKOPT = MAX( 1, ( NB+2 )*N )
         WORK( 1 ) = LWKOPT
c
         IF( LWORK.LT.MAX( 1, 3*N-1 ) .AND. .NOT.LQUERY )
     $      INFO = -8
      END IF
c
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSYEV ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
c
c     Quick return if possible
c
      IF( N.EQ.0 ) THEN
         RETURN
      END IF
c
      IF( N.EQ.1 ) THEN
         W( 1 ) = A( 1, 1 )
         WORK( 1 ) = 2
         IF( WANTZ )
     $      A( 1, 1 ) = ONE
         RETURN
      END IF
c
c     Get machine constants.
c
      SAFMIN = DLAMCH( 'Safe minimum' )
      EPS = DLAMCH( 'Precision' )
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN = SQRT( SMLNUM )
      RMAX = SQRT( BIGNUM )
c
c     Scale matrix to allowable range, if necessary.
c
      ANRM = DLANSY( 'M', UPLO, N, A, LDA, WORK )
      ISCALE = 0
      IF( ANRM.GT.ZERO .AND. ANRM.LT.RMIN ) THEN
         ISCALE = 1
         SIGMA = RMIN / ANRM
      ELSE IF( ANRM.GT.RMAX ) THEN
         ISCALE = 1
         SIGMA = RMAX / ANRM
      END IF
      IF( ISCALE.EQ.1 )
     $   CALL DLASCL( UPLO, 0, 0, ONE, SIGMA, N, N, A, LDA, INFO )
c
c     Call DSYTRD to reduce symmetric matrix to tridiagonal form.
c
      INDE = 1
      INDTAU = INDE + N
      INDWRK = INDTAU + N
      LLWORK = LWORK - INDWRK + 1
      CALL DSYTRD( UPLO, N, A, LDA, W, WORK( INDE ), WORK( INDTAU ),
     $             WORK( INDWRK ), LLWORK, IINFO )
c
c     For eigenvalues only, call DSTERF.  For eigenvectors, first call
c     DORGTR to generate the orthogonal matrix, then call DSTEQR.
c
      IF( .NOT.WANTZ ) THEN
         CALL DSTERF( N, W, WORK( INDE ), INFO )
      ELSE
         CALL DORGTR( UPLO, N, A, LDA, WORK( INDTAU ), WORK( INDWRK ),
     $                LLWORK, IINFO )
         CALL DSTEQR( JOBZ, N, W, WORK( INDE ), A, LDA, WORK( INDTAU ),
     $                INFO )
      END IF
c
c     If matrix was scaled, then rescale eigenvalues appropriately.
c
      IF( ISCALE.EQ.1 ) THEN
         IF( INFO.EQ.0 ) THEN
            IMAX = N
         ELSE
            IMAX = INFO - 1
         END IF
         CALL DSCAL( IMAX, ONE / SIGMA, W, 1 )
      END IF
c
c     Set WORK(1) to optimal workspace size.
c
      WORK( 1 ) = LWKOPT
c
      RETURN
c
c     End of DSYEV
c
      END

      SUBROUTINE DSYTRD( UPLO, N, A, LDA, D, E, TAU, WORK, LWORK, INFO )
c
c  -- LAPACK computational routine (version 3.4.0) --
c  -- LAPACK is a software package provided by Univ. of Tennessee,    --
c  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
c     November 2011
c
c     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, LWORK, N
c     ..
c     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), D( * ), E( * ), TAU( * ),
     $                   WORK( * )
c     ..
c
c  =====================================================================
c
c     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
c     ..
c     .. Local Scalars ..
      LOGICAL            LQUERY, UPPER
      INTEGER            I, IINFO, IWS, J, KK, LDWORK, LWKOPT, NB,
     $                   NBMIN, NX
c     ..
c     .. External Subroutines ..
      EXTERNAL           DLATRD, DSYR2K, DSYTD2, XERBLA
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          MAX
c     ..
c     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
c     ..
c     .. Executable Statements ..
c
c     Test the input parameters
c
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      LQUERY = ( LWORK.EQ.-1 )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LWORK.LT.1 .AND. .NOT.LQUERY ) THEN
         INFO = -9
      END IF
c
      IF( INFO.EQ.0 ) THEN
c
c        Determine the block size.
c
         NB = ILAENV( 1, 'DSYTRD', UPLO, N, -1, -1, -1 )
         LWKOPT = N*NB
         WORK( 1 ) = LWKOPT
      END IF
c
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSYTRD', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
c
c     Quick return if possible
c
      IF( N.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
c
      NX = N
      IWS = 1
      IF( NB.GT.1 .AND. NB.LT.N ) THEN
c
c        Determine when to cross over from blocked to unblocked code
c        (last block is always handled by unblocked code).
c
         NX = MAX( NB, ILAENV( 3, 'DSYTRD', UPLO, N, -1, -1, -1 ) )
         IF( NX.LT.N ) THEN
c
c           Determine if workspace is large enough for blocked code.
c
            LDWORK = N
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
c
c              Not enough workspace to use optimal NB:  determine the
c              minimum value of NB, and reduce NB or force use of
c              unblocked code by setting NX = N.
c
               NB = MAX( LWORK / LDWORK, 1 )
               NBMIN = ILAENV( 2, 'DSYTRD', UPLO, N, -1, -1, -1 )
               IF( NB.LT.NBMIN )
     $            NX = N
            END IF
         ELSE
            NX = N
         END IF
      ELSE
         NB = 1
      END IF
c
      IF( UPPER ) THEN
c
c        Reduce the upper triangle of A.
c        Columns 1:kk are handled by the unblocked method.
c
         KK = N - ( ( N-NX+NB-1 ) / NB )*NB
         DO 20 I = N - NB + 1, KK + 1, -NB
c
c           Reduce columns i:i+nb-1 to tridiagonal form and form the
c           matrix W which is needed to update the unreduced part of
c           the matrix
c
            CALL DLATRD( UPLO, I+NB-1, NB, A, LDA, E, TAU, WORK,
     $                   LDWORK )
c
c           Update the unreduced submatrix A(1:i-1,1:i-1), using an
c           update of the form:  A := A - V*W**T - W*V**T
c
            CALL DSYR2K( UPLO, 'No transpose', I-1, NB, -ONE, A( 1, I ),
     $                   LDA, WORK, LDWORK, ONE, A, LDA )
c
c           Copy superdiagonal elements back into A, and diagonal
c           elements into D
c
            DO 10 J = I, I + NB - 1
               A( J-1, J ) = E( J-1 )
               D( J ) = A( J, J )
   10       CONTINUE
   20    CONTINUE
c
c        Use unblocked code to reduce the last or only block
c
         CALL DSYTD2( UPLO, KK, A, LDA, D, E, TAU, IINFO )
      ELSE
c
c        Reduce the lower triangle of A
c
         DO 40 I = 1, N - NX, NB
c
c           Reduce columns i:i+nb-1 to tridiagonal form and form the
c           matrix W which is needed to update the unreduced part of
c           the matrix
c
            CALL DLATRD( UPLO, N-I+1, NB, A( I, I ), LDA, E( I ),
     $                   TAU( I ), WORK, LDWORK )
c
c           Update the unreduced submatrix A(i+ib:n,i+ib:n), using
c           an update of the form:  A := A - V*W**T - W*V**T
c
            CALL DSYR2K( UPLO, 'No transpose', N-I-NB+1, NB, -ONE,
     $                   A( I+NB, I ), LDA, WORK( NB+1 ), LDWORK, ONE,
     $                   A( I+NB, I+NB ), LDA )
c
c           Copy subdiagonal elements back into A, and diagonal
c           elements into D
c
            DO 30 J = I, I + NB - 1
               A( J+1, J ) = E( J )
               D( J ) = A( J, J )
   30       CONTINUE
   40    CONTINUE
c
c        Use unblocked code to reduce the last or only block
c
         CALL DSYTD2( UPLO, N-I+1, A( I, I ), LDA, D( I ), E( I ),
     $                TAU( I ), IINFO )
      END IF
c
      WORK( 1 ) = LWKOPT
      RETURN
c
c     End of DSYTRD
c
      END

      SUBROUTINE DSYTD2( UPLO, N, A, LDA, D, E, TAU, INFO )
c
c  -- LAPACK computational routine (version 3.4.2) --
c  -- LAPACK is a software package provided by Univ. of Tennessee,    --
c  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
c     September 2012
c
c     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
c     ..
c     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), D( * ), E( * ), TAU( * )
c     ..
c
c  =====================================================================
c
c     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO, HALF
      PARAMETER          ( ONE = 1.0D0, ZERO = 0.0D0,
     $                   HALF = 1.0D0 / 2.0D0 )
c     ..
c     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I
      DOUBLE PRECISION   ALPHA, TAUI
c     ..
c     .. External Subroutines ..
      EXTERNAL           DAXPY, DLARFG, DSYMV, DSYR2, XERBLA
c     ..
c     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DDOT
      EXTERNAL           LSAME, DDOT
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
c     ..
c     .. Executable Statements ..
c
c     Test the input parameters
c
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSYTD2', -INFO )
         RETURN
      END IF
c
c     Quick return if possible
c
      IF( N.LE.0 )
     $   RETURN
c
      IF( UPPER ) THEN
c
c        Reduce the upper triangle of A
c
         DO 10 I = N - 1, 1, -1
c
c           Generate elementary reflector H(i) = I - tau * v * v**T
c           to annihilate A(1:i-1,i+1)
c
            CALL DLARFG( I, A( I, I+1 ), A( 1, I+1 ), 1, TAUI )
            E( I ) = A( I, I+1 )
c
            IF( TAUI.NE.ZERO ) THEN
c
c              Apply H(i) from both sides to A(1:i,1:i)
c
               A( I, I+1 ) = ONE
c
c              Compute  x := tau * A * v  storing x in TAU(1:i)
c
               CALL DSYMV( UPLO, I, TAUI, A, LDA, A( 1, I+1 ), 1, ZERO,
     $                     TAU, 1 )
c
c              Compute  w := x - 1/2 * tau * (x**T * v) * v
c
               ALPHA = -HALF*TAUI*DDOT( I, TAU, 1, A( 1, I+1 ), 1 )
               CALL DAXPY( I, ALPHA, A( 1, I+1 ), 1, TAU, 1 )
c
c              Apply the transformation as a rank-2 update:
c                 A := A - v * w**T - w * v**T
c
               CALL DSYR2( UPLO, I, -ONE, A( 1, I+1 ), 1, TAU, 1, A,
     $                     LDA )
c
               A( I, I+1 ) = E( I )
            END IF
            D( I+1 ) = A( I+1, I+1 )
            TAU( I ) = TAUI
   10    CONTINUE
         D( 1 ) = A( 1, 1 )
      ELSE
c
c        Reduce the lower triangle of A
c
         DO 20 I = 1, N - 1
c
c           Generate elementary reflector H(i) = I - tau * v * v**T
c           to annihilate A(i+2:n,i)
c
            CALL DLARFG( N-I, A( I+1, I ), A( MIN( I+2, N ), I ), 1,
     $                   TAUI )
            E( I ) = A( I+1, I )
c
            IF( TAUI.NE.ZERO ) THEN
c
c              Apply H(i) from both sides to A(i+1:n,i+1:n)
c
               A( I+1, I ) = ONE
c
c              Compute  x := tau * A * v  storing y in TAU(i:n-1)
c
               CALL DSYMV( UPLO, N-I, TAUI, A( I+1, I+1 ), LDA,
     $                     A( I+1, I ), 1, ZERO, TAU( I ), 1 )
c
c              Compute  w := x - 1/2 * tau * (x**T * v) * v
c
               ALPHA = -HALF*TAUI*DDOT( N-I, TAU( I ), 1, A( I+1, I ),
     $                 1 )
               CALL DAXPY( N-I, ALPHA, A( I+1, I ), 1, TAU( I ), 1 )
c
c              Apply the transformation as a rank-2 update:
c                 A := A - v * w**T - w * v**T
c
               CALL DSYR2( UPLO, N-I, -ONE, A( I+1, I ), 1, TAU( I ), 1,
     $                     A( I+1, I+1 ), LDA )
c
               A( I+1, I ) = E( I )
            END IF
            D( I ) = A( I, I )
            TAU( I ) = TAUI
   20    CONTINUE
         D( N ) = A( N, N )
      END IF
c
      RETURN
c
c     End of DSYTD2
c
      END

      SUBROUTINE DLASCL( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )
c
c  -- LAPACK auxiliary routine (version 3.4.2) --
c  -- LAPACK is a software package provided by Univ. of Tennessee,    --
c  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
c     September 2012
c
c     .. Scalar Arguments ..
      CHARACTER          TYPE
      INTEGER            INFO, KL, KU, LDA, M, N
      DOUBLE PRECISION   CFROM, CTO
c     ..
c     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
c     ..
c
c  =====================================================================
c
c     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
c     ..
c     .. Local Scalars ..
      LOGICAL            DONE
      INTEGER            I, ITYPE, J, K1, K2, K3, K4
      DOUBLE PRECISION   BIGNUM, CFROM1, CFROMC, CTO1, CTOC, MUL, SMLNUM
c     ..
c     .. External Functions ..
      LOGICAL            LSAME, DISNAN
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           LSAME, DLAMCH, DISNAN
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
c     ..
c     .. External Subroutines ..
      EXTERNAL           XERBLA
c     ..
c     .. Executable Statements ..
c
c     Test the input arguments
c
      INFO = 0
c
      IF( LSAME( TYPE, 'G' ) ) THEN
         ITYPE = 0
      ELSE IF( LSAME( TYPE, 'L' ) ) THEN
         ITYPE = 1
      ELSE IF( LSAME( TYPE, 'U' ) ) THEN
         ITYPE = 2
      ELSE IF( LSAME( TYPE, 'H' ) ) THEN
         ITYPE = 3
      ELSE IF( LSAME( TYPE, 'B' ) ) THEN
         ITYPE = 4
      ELSE IF( LSAME( TYPE, 'Q' ) ) THEN
         ITYPE = 5
      ELSE IF( LSAME( TYPE, 'Z' ) ) THEN
         ITYPE = 6
      ELSE
         ITYPE = -1
      END IF
c
      IF( ITYPE.EQ.-1 ) THEN
         INFO = -1
      ELSE IF( CFROM.EQ.ZERO .OR. DISNAN(CFROM) ) THEN
         INFO = -4
      ELSE IF( DISNAN(CTO) ) THEN
         INFO = -5
      ELSE IF( M.LT.0 ) THEN
         INFO = -6
      ELSE IF( N.LT.0 .OR. ( ITYPE.EQ.4 .AND. N.NE.M ) .OR.
     $         ( ITYPE.EQ.5 .AND. N.NE.M ) ) THEN
         INFO = -7
      ELSE IF( ITYPE.LE.3 .AND. LDA.LT.MAX( 1, M ) ) THEN
         INFO = -9
      ELSE IF( ITYPE.GE.4 ) THEN
         IF( KL.LT.0 .OR. KL.GT.MAX( M-1, 0 ) ) THEN
            INFO = -2
         ELSE IF( KU.LT.0 .OR. KU.GT.MAX( N-1, 0 ) .OR.
     $            ( ( ITYPE.EQ.4 .OR. ITYPE.EQ.5 ) .AND. KL.NE.KU ) )
     $             THEN
            INFO = -3
         ELSE IF( ( ITYPE.EQ.4 .AND. LDA.LT.KL+1 ) .OR.
     $            ( ITYPE.EQ.5 .AND. LDA.LT.KU+1 ) .OR.
     $            ( ITYPE.EQ.6 .AND. LDA.LT.2*KL+KU+1 ) ) THEN
            INFO = -9
         END IF
      END IF
c
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLASCL', -INFO )
         RETURN
      END IF
c
c     Quick return if possible
c
      IF( N.EQ.0 .OR. M.EQ.0 )
     $   RETURN
c
c     Get machine parameters
c
      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
c
      CFROMC = CFROM
      CTOC = CTO
c
   10 CONTINUE
      CFROM1 = CFROMC*SMLNUM
      IF( CFROM1.EQ.CFROMC ) THEN
!        CFROMC is an inf.  Multiply by a correctly signed zero for
!        finite CTOC, or a NaN if CTOC is infinite.
         MUL = CTOC / CFROMC
         DONE = .TRUE.
         CTO1 = CTOC
      ELSE
         CTO1 = CTOC / BIGNUM
         IF( CTO1.EQ.CTOC ) THEN
!           CTOC is either 0 or an inf.  In both cases, CTOC itself
!           serves as the correct multiplication factor.
            MUL = CTOC
            DONE = .TRUE.
            CFROMC = ONE
         ELSE IF( ABS( CFROM1 ).GT.ABS( CTOC ) .AND. CTOC.NE.ZERO ) THEN
            MUL = SMLNUM
            DONE = .FALSE.
            CFROMC = CFROM1
         ELSE IF( ABS( CTO1 ).GT.ABS( CFROMC ) ) THEN
            MUL = BIGNUM
            DONE = .FALSE.
            CTOC = CTO1
         ELSE
            MUL = CTOC / CFROMC
            DONE = .TRUE.
         END IF
      END IF
c
      IF( ITYPE.EQ.0 ) THEN
c
c        Full matrix
c
         DO 30 J = 1, N
            DO 20 I = 1, M
               A( I, J ) = A( I, J )*MUL
   20       CONTINUE
   30    CONTINUE
c
      ELSE IF( ITYPE.EQ.1 ) THEN
c
c        Lower triangular matrix
c
         DO 50 J = 1, N
            DO 40 I = J, M
               A( I, J ) = A( I, J )*MUL
   40       CONTINUE
   50    CONTINUE
c
      ELSE IF( ITYPE.EQ.2 ) THEN
c
c        Upper triangular matrix
c
         DO 70 J = 1, N
            DO 60 I = 1, MIN( J, M )
               A( I, J ) = A( I, J )*MUL
   60       CONTINUE
   70    CONTINUE
c
      ELSE IF( ITYPE.EQ.3 ) THEN
c
c        Upper Hessenberg matrix
c
         DO 90 J = 1, N
            DO 80 I = 1, MIN( J+1, M )
               A( I, J ) = A( I, J )*MUL
   80       CONTINUE
   90    CONTINUE
c
      ELSE IF( ITYPE.EQ.4 ) THEN
c
c        Lower half of a symmetric band matrix
c
         K3 = KL + 1
         K4 = N + 1
         DO 110 J = 1, N
            DO 100 I = 1, MIN( K3, K4-J )
               A( I, J ) = A( I, J )*MUL
  100       CONTINUE
  110    CONTINUE
c
      ELSE IF( ITYPE.EQ.5 ) THEN
c
c        Upper half of a symmetric band matrix
c
         K1 = KU + 2
         K3 = KU + 1
         DO 130 J = 1, N
            DO 120 I = MAX( K1-J, 1 ), K3
               A( I, J ) = A( I, J )*MUL
  120       CONTINUE
  130    CONTINUE
c
      ELSE IF( ITYPE.EQ.6 ) THEN
c
c        Band matrix
c
         K1 = KL + KU + 2
         K2 = KL + 1
         K3 = 2*KL + KU + 1
         K4 = KL + KU + 1 + M
         DO 150 J = 1, N
            DO 140 I = MAX( K1-J, K2 ), MIN( K3, K4-J )
               A( I, J ) = A( I, J )*MUL
  140       CONTINUE
  150    CONTINUE
c
      END IF
c
      IF( .NOT.DONE )
     $   GO TO 10
c
      RETURN
c
c     End of DLASCL
c
      END

      DOUBLE PRECISION FUNCTION DLANSY( NORM, UPLO, N, A, LDA, WORK )
c
c  -- LAPACK auxiliary routine (version 3.4.2) --
c  -- LAPACK is a software package provided by Univ. of Tennessee,    --
c  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
c     September 2012
c
c     .. Scalar Arguments ..
      CHARACTER          NORM, UPLO
      INTEGER            LDA, N
c     ..
c     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), WORK( * )
c     ..
c
c =====================================================================
c
c     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
c     ..
c     .. Local Scalars ..
      INTEGER            I, J
      DOUBLE PRECISION   ABSA, SCALE, SUM, VALUE
c     ..
c     .. External Subroutines ..
      EXTERNAL           DLASSQ
c     ..
c     .. External Functions ..
      LOGICAL            LSAME, DISNAN
      EXTERNAL           LSAME, DISNAN
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          ABS, SQRT
c     ..
c     .. Executable Statements ..
c
      IF( N.EQ.0 ) THEN
         VALUE = ZERO
      ELSE IF( LSAME( NORM, 'M' ) ) THEN
c
c        Find max(abs(A(i,j))).
c
         VALUE = ZERO
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO 20 J = 1, N
               DO 10 I = 1, J
                  SUM = ABS( A( I, J ) )
                  IF( VALUE .LT. SUM .OR. DISNAN( SUM ) ) VALUE = SUM
   10          CONTINUE
   20       CONTINUE
         ELSE
            DO 40 J = 1, N
               DO 30 I = J, N
                  SUM = ABS( A( I, J ) )
                  IF( VALUE .LT. SUM .OR. DISNAN( SUM ) ) VALUE = SUM
   30          CONTINUE
   40       CONTINUE
         END IF
      ELSE IF( ( LSAME( NORM, 'I' ) ) .OR. ( LSAME( NORM, 'O' ) ) .OR.
     $         ( NORM.EQ.'1' ) ) THEN
c
c        Find normI(A) ( = norm1(A), since A is symmetric).
c
         VALUE = ZERO
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO 60 J = 1, N
               SUM = ZERO
               DO 50 I = 1, J - 1
                  ABSA = ABS( A( I, J ) )
                  SUM = SUM + ABSA
                  WORK( I ) = WORK( I ) + ABSA
   50          CONTINUE
               WORK( J ) = SUM + ABS( A( J, J ) )
   60       CONTINUE
            DO 70 I = 1, N
               SUM = WORK( I )
               IF( VALUE .LT. SUM .OR. DISNAN( SUM ) ) VALUE = SUM
   70       CONTINUE
         ELSE
            DO 80 I = 1, N
               WORK( I ) = ZERO
   80       CONTINUE
            DO 100 J = 1, N
               SUM = WORK( J ) + ABS( A( J, J ) )
               DO 90 I = J + 1, N
                  ABSA = ABS( A( I, J ) )
                  SUM = SUM + ABSA
                  WORK( I ) = WORK( I ) + ABSA
   90          CONTINUE
               IF( VALUE .LT. SUM .OR. DISNAN( SUM ) ) VALUE = SUM
  100       CONTINUE
         END IF
      ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
c
c        Find normF(A).
c
         SCALE = ZERO
         SUM = ONE
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO 110 J = 2, N
               CALL DLASSQ( J-1, A( 1, J ), 1, SCALE, SUM )
  110       CONTINUE
         ELSE
            DO 120 J = 1, N - 1
               CALL DLASSQ( N-J, A( J+1, J ), 1, SCALE, SUM )
  120       CONTINUE
         END IF
         SUM = 2*SUM
         CALL DLASSQ( N, A, LDA+1, SCALE, SUM )
         VALUE = SCALE*SQRT( SUM )
      END IF
c
      DLANSY = VALUE
      RETURN
c
c     End of DLANSY
c
      END

      DOUBLE PRECISION FUNCTION DLANST( NORM, N, D, E )
c
c  -- LAPACK auxiliary routine (version 3.4.2) --
c  -- LAPACK is a software package provided by Univ. of Tennessee,    --
c  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
c     September 2012
c
c     .. Scalar Arguments ..
      CHARACTER          NORM
      INTEGER            N
c     ..
c     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), E( * )
c     ..
c
c  =====================================================================
c
c     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
c     ..
c     .. Local Scalars ..
      INTEGER            I
      DOUBLE PRECISION   ANORM, SCALE, SUM
c     ..
c     .. External Functions ..
      LOGICAL            LSAME, DISNAN
      EXTERNAL           LSAME, DISNAN
c     ..
c     .. External Subroutines ..
      EXTERNAL           DLASSQ
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          ABS, SQRT
c     ..
c     .. Executable Statements ..
c
      IF( N.LE.0 ) THEN
         ANORM = ZERO
      ELSE IF( LSAME( NORM, 'M' ) ) THEN
c
c        Find max(abs(A(i,j))).
c
         ANORM = ABS( D( N ) )
         DO 10 I = 1, N - 1
            SUM = ABS( D( I ) )
            IF( ANORM .LT. SUM .OR. DISNAN( SUM ) ) ANORM = SUM
            SUM = ABS( E( I ) )
            IF( ANORM .LT. SUM .OR. DISNAN( SUM ) ) ANORM = SUM
   10    CONTINUE
      ELSE IF( LSAME( NORM, 'O' ) .OR. NORM.EQ.'1' .OR.
     $         LSAME( NORM, 'I' ) ) THEN
c
c        Find norm1(A).
c
         IF( N.EQ.1 ) THEN
            ANORM = ABS( D( 1 ) )
         ELSE
            ANORM = ABS( D( 1 ) )+ABS( E( 1 ) )
            SUM = ABS( E( N-1 ) )+ABS( D( N ) )
            IF( ANORM .LT. SUM .OR. DISNAN( SUM ) ) ANORM = SUM
            DO 20 I = 2, N - 1
               SUM = ABS( D( I ) )+ABS( E( I ) )+ABS( E( I-1 ) )
               IF( ANORM .LT. SUM .OR. DISNAN( SUM ) ) ANORM = SUM
   20       CONTINUE
         END IF
      ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
c
c        Find normF(A).
c
         SCALE = ZERO
         SUM = ONE
         IF( N.GT.1 ) THEN
            CALL DLASSQ( N-1, E, 1, SCALE, SUM )
            SUM = 2*SUM
         END IF
         CALL DLASSQ( N, D, 1, SCALE, SUM )
         ANORM = SCALE*SQRT( SUM )
      END IF
c
      DLANST = ANORM
      RETURN
c
c     End of DLANST
c
      END

      SUBROUTINE DLATRD( UPLO, N, NB, A, LDA, E, TAU, W, LDW )
c
c  -- LAPACK auxiliary routine (version 3.4.2) --
c  -- LAPACK is a software package provided by Univ. of Tennessee,    --
c  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
c     September 2012
c
c     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, LDW, N, NB
c     ..
c     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), E( * ), TAU( * ), W( LDW, * )
c     ..
c
c  =====================================================================
c
c     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, HALF
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, HALF = 0.5D+0 )
c     ..
c     .. Local Scalars ..
      INTEGER            I, IW
      DOUBLE PRECISION   ALPHA
c     ..
c     .. External Subroutines ..
      EXTERNAL           DAXPY, DGEMV, DLARFG, DSCAL, DSYMV
c     ..
c     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DDOT
      EXTERNAL           LSAME, DDOT
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          MIN
c     ..
c     .. Executable Statements ..
c
c     Quick return if possible
c
      IF( N.LE.0 )
     $   RETURN
c
      IF( LSAME( UPLO, 'U' ) ) THEN
c
c        Reduce last NB columns of upper triangle
c
         DO 10 I = N, N - NB + 1, -1
            IW = I - N + NB
            IF( I.LT.N ) THEN
c
c              Update A(1:i,i)
c
               CALL DGEMV( 'No transpose', I, N-I, -ONE, A( 1, I+1 ),
     $                     LDA, W( I, IW+1 ), LDW, ONE, A( 1, I ), 1 )
               CALL DGEMV( 'No transpose', I, N-I, -ONE, W( 1, IW+1 ),
     $                     LDW, A( I, I+1 ), LDA, ONE, A( 1, I ), 1 )
            END IF
            IF( I.GT.1 ) THEN
c
c              Generate elementary reflector H(i) to annihilate
c              A(1:i-2,i)
c
               CALL DLARFG( I-1, A( I-1, I ), A( 1, I ), 1, TAU( I-1 ) )
               E( I-1 ) = A( I-1, I )
               A( I-1, I ) = ONE
c
c              Compute W(1:i-1,i)
c
               CALL DSYMV( 'Upper', I-1, ONE, A, LDA, A( 1, I ), 1,
     $                     ZERO, W( 1, IW ), 1 )
               IF( I.LT.N ) THEN
                  CALL DGEMV( 'Transpose', I-1, N-I, ONE, W( 1, IW+1 ),
     $                        LDW, A( 1, I ), 1, ZERO, W( I+1, IW ), 1 )
                  CALL DGEMV( 'No transpose', I-1, N-I, -ONE,
     $                        A( 1, I+1 ), LDA, W( I+1, IW ), 1, ONE,
     $                        W( 1, IW ), 1 )
                  CALL DGEMV( 'Transpose', I-1, N-I, ONE, A( 1, I+1 ),
     $                        LDA, A( 1, I ), 1, ZERO, W( I+1, IW ), 1 )
                  CALL DGEMV( 'No transpose', I-1, N-I, -ONE,
     $                        W( 1, IW+1 ), LDW, W( I+1, IW ), 1, ONE,
     $                        W( 1, IW ), 1 )
               END IF
               CALL DSCAL( I-1, TAU( I-1 ), W( 1, IW ), 1 )
               ALPHA = -HALF*TAU( I-1 )*DDOT( I-1, W( 1, IW ), 1,
     $                 A( 1, I ), 1 )
               CALL DAXPY( I-1, ALPHA, A( 1, I ), 1, W( 1, IW ), 1 )
            END IF
c
   10    CONTINUE
      ELSE
c
c        Reduce first NB columns of lower triangle
c
         DO 20 I = 1, NB
c
c           Update A(i:n,i)
c
            CALL DGEMV( 'No transpose', N-I+1, I-1, -ONE, A( I, 1 ),
     $                  LDA, W( I, 1 ), LDW, ONE, A( I, I ), 1 )
            CALL DGEMV( 'No transpose', N-I+1, I-1, -ONE, W( I, 1 ),
     $                  LDW, A( I, 1 ), LDA, ONE, A( I, I ), 1 )
            IF( I.LT.N ) THEN
c
c              Generate elementary reflector H(i) to annihilate
c              A(i+2:n,i)
c
               CALL DLARFG( N-I, A( I+1, I ), A( MIN( I+2, N ), I ), 1,
     $                      TAU( I ) )
               E( I ) = A( I+1, I )
               A( I+1, I ) = ONE
c
c              Compute W(i+1:n,i)
c
               CALL DSYMV( 'Lower', N-I, ONE, A( I+1, I+1 ), LDA,
     $                     A( I+1, I ), 1, ZERO, W( I+1, I ), 1 )
               CALL DGEMV( 'Transpose', N-I, I-1, ONE, W( I+1, 1 ), LDW,
     $                     A( I+1, I ), 1, ZERO, W( 1, I ), 1 )
               CALL DGEMV( 'No transpose', N-I, I-1, -ONE, A( I+1, 1 ),
     $                     LDA, W( 1, I ), 1, ONE, W( I+1, I ), 1 )
               CALL DGEMV( 'Transpose', N-I, I-1, ONE, A( I+1, 1 ), LDA,
     $                     A( I+1, I ), 1, ZERO, W( 1, I ), 1 )
               CALL DGEMV( 'No transpose', N-I, I-1, -ONE, W( I+1, 1 ),
     $                     LDW, W( 1, I ), 1, ONE, W( I+1, I ), 1 )
               CALL DSCAL( N-I, TAU( I ), W( I+1, I ), 1 )
               ALPHA = -HALF*TAU( I )*DDOT( N-I, W( I+1, I ), 1,
     $                 A( I+1, I ), 1 )
               CALL DAXPY( N-I, ALPHA, A( I+1, I ), 1, W( I+1, I ), 1 )
            END IF
c
   20    CONTINUE
      END IF
c
      RETURN
c
c     End of DLATRD
c
      END

      SUBROUTINE DLARFT( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT )
c
c  -- LAPACK auxiliary routine (version 3.4.2) --
c  -- LAPACK is a software package provided by Univ. of Tennessee,    --
c  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
c     September 2012
c
c     .. Scalar Arguments ..
      CHARACTER          DIRECT, STOREV
      INTEGER            K, LDT, LDV, N
c     ..
c     .. Array Arguments ..
      DOUBLE PRECISION   T( LDT, * ), TAU( * ), V( LDV, * )
c     ..
c
c  =====================================================================
c
c     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
c     ..
c     .. Local Scalars ..
      INTEGER            I, J, PREVLASTV, LASTV
c     ..
c     .. External Subroutines ..
      EXTERNAL           DGEMV, DTRMV
c     ..
c     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
c     ..
c     .. Executable Statements ..
c
c     Quick return if possible
c
      IF( N.EQ.0 )
     $   RETURN
c
      IF( LSAME( DIRECT, 'F' ) ) THEN
         PREVLASTV = N
         DO I = 1, K
            PREVLASTV = MAX( I, PREVLASTV )
            IF( TAU( I ).EQ.ZERO ) THEN
c
c              H(i)  =  I
c
               DO J = 1, I
                  T( J, I ) = ZERO
               END DO
            ELSE
c
c              general case
c
               IF( LSAME( STOREV, 'C' ) ) THEN
c                 Skip any trailing zeros.
                  DO LASTV = N, I+1, -1
                     IF( V( LASTV, I ).NE.ZERO ) EXIT
                  END DO
                  DO J = 1, I-1
                     T( J, I ) = -TAU( I ) * V( I , J )
                  END DO   
                  J = MIN( LASTV, PREVLASTV )
c
c                 T(1:i-1,i) := - tau(i) * V(i:j,1:i-1)**T * V(i:j,i)
c
                  CALL DGEMV( 'Transpose', J-I, I-1, -TAU( I ), 
     $                        V( I+1, 1 ), LDV, V( I+1, I ), 1, ONE, 
     $                        T( 1, I ), 1 )
               ELSE
c                 Skip any trailing zeros.
                  DO LASTV = N, I+1, -1
                     IF( V( I, LASTV ).NE.ZERO ) EXIT
                  END DO
                  DO J = 1, I-1
                     T( J, I ) = -TAU( I ) * V( J , I )
                  END DO   
                  J = MIN( LASTV, PREVLASTV )
c
c                 T(1:i-1,i) := - tau(i) * V(1:i-1,i:j) * V(i,i:j)**T
c
                  CALL DGEMV( 'No transpose', I-1, J-I, -TAU( I ),
     $                        V( 1, I+1 ), LDV, V( I, I+1 ), LDV, ONE,
     $                        T( 1, I ), 1 )
               END IF
c
c              T(1:i-1,i) := T(1:i-1,1:i-1) * T(1:i-1,i)
c
               CALL DTRMV( 'Upper', 'No transpose', 'Non-unit', I-1, T,
     $                     LDT, T( 1, I ), 1 )
               T( I, I ) = TAU( I )
               IF( I.GT.1 ) THEN
                  PREVLASTV = MAX( PREVLASTV, LASTV )
               ELSE
                  PREVLASTV = LASTV
               END IF
            END IF
         END DO
      ELSE
         PREVLASTV = 1
         DO I = K, 1, -1
            IF( TAU( I ).EQ.ZERO ) THEN
c
c              H(i)  =  I
c
               DO J = I, K
                  T( J, I ) = ZERO
               END DO
            ELSE
c
c              general case
c
               IF( I.LT.K ) THEN
                  IF( LSAME( STOREV, 'C' ) ) THEN
c                    Skip any leading zeros.
                     DO LASTV = 1, I-1
                        IF( V( LASTV, I ).NE.ZERO ) EXIT
                     END DO
                     DO J = I+1, K
                        T( J, I ) = -TAU( I ) * V( N-K+I , J )
                     END DO   
                     J = MAX( LASTV, PREVLASTV )
c
c                    T(i+1:k,i) = -tau(i) * V(j:n-k+i,i+1:k)**T * V(j:n-k+i,i)
c
                     CALL DGEMV( 'Transpose', N-K+I-J, K-I, -TAU( I ),
     $                           V( J, I+1 ), LDV, V( J, I ), 1, ONE,
     $                           T( I+1, I ), 1 )
                  ELSE
c                    Skip any leading zeros.
                     DO LASTV = 1, I-1
                        IF( V( I, LASTV ).NE.ZERO ) EXIT
                     END DO
                     DO J = I+1, K
                        T( J, I ) = -TAU( I ) * V( J, N-K+I )
                     END DO   
                     J = MAX( LASTV, PREVLASTV )
c
c                    T(i+1:k,i) = -tau(i) * V(i+1:k,j:n-k+i) * V(i,j:n-k+i)**T
c
                     CALL DGEMV( 'No transpose', K-I, N-K+I-J,
     $                    -TAU( I ), V( I+1, J ), LDV, V( I, J ), LDV,
     $                    ONE, T( I+1, I ), 1 )
                  END IF
c
c                 T(i+1:k,i) := T(i+1:k,i+1:k) * T(i+1:k,i)
c
                  CALL DTRMV( 'Lower', 'No transpose', 'Non-unit', K-I,
     $                        T( I+1, I+1 ), LDT, T( I+1, I ), 1 )
                  IF( I.GT.1 ) THEN
                     PREVLASTV = MIN( PREVLASTV, LASTV )
                  ELSE
                     PREVLASTV = LASTV
                  END IF
               END IF
               T( I, I ) = TAU( I )
            END IF
         END DO
      END IF
      RETURN
c
c     End of DLARFT
c
      END

      SUBROUTINE DLARFG( N, ALPHA, X, INCX, TAU )
c
c  -- LAPACK auxiliary routine (version 3.4.2) --
c  -- LAPACK is a software package provided by Univ. of Tennessee,    --
c  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
c     September 2012
c
c     .. Scalar Arguments ..
      INTEGER            INCX, N
      DOUBLE PRECISION   ALPHA, TAU
c     ..
c     .. Array Arguments ..
      DOUBLE PRECISION   X( * )
c     ..
c
c  =====================================================================
c
c     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
c     ..
c     .. Local Scalars ..
      INTEGER            J, KNT
      DOUBLE PRECISION   BETA, RSAFMN, SAFMIN, XNORM
c     ..
c     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLAPY2, DNRM2
      EXTERNAL           DLAMCH, DLAPY2, DNRM2
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          ABS, SIGN
c     ..
c     .. External Subroutines ..
      EXTERNAL           DSCAL
c     ..
c     .. Executable Statements ..
c
      IF( N.LE.1 ) THEN
         TAU = ZERO
         RETURN
      END IF
c
      XNORM = DNRM2( N-1, X, INCX )
c
      IF( XNORM.EQ.ZERO ) THEN
c
c        H  =  I
c
         TAU = ZERO
      ELSE
c
c        general case
c
         BETA = -SIGN( DLAPY2( ALPHA, XNORM ), ALPHA )
         SAFMIN = DLAMCH( 'S' ) / DLAMCH( 'E' )
         KNT = 0
         IF( ABS( BETA ).LT.SAFMIN ) THEN
c
c           XNORM, BETA may be inaccurate; scale X and recompute them
c
            RSAFMN = ONE / SAFMIN
   10       CONTINUE
            KNT = KNT + 1
            CALL DSCAL( N-1, RSAFMN, X, INCX )
            BETA = BETA*RSAFMN
            ALPHA = ALPHA*RSAFMN
            IF( ABS( BETA ).LT.SAFMIN )
     $         GO TO 10
c
c           New BETA is at most 1, at least SAFMIN
c
            XNORM = DNRM2( N-1, X, INCX )
            BETA = -SIGN( DLAPY2( ALPHA, XNORM ), ALPHA )
         END IF
         TAU = ( BETA-ALPHA ) / BETA
         CALL DSCAL( N-1, ONE / ( ALPHA-BETA ), X, INCX )
c
c        If ALPHA is subnormal, it may lose relative accuracy
c
         DO 20 J = 1, KNT
            BETA = BETA*SAFMIN
 20      CONTINUE
         ALPHA = BETA
      END IF
c
      RETURN
c
c     End of DLARFG
c
      END

      SUBROUTINE DLASSQ( N, X, INCX, SCALE, SUMSQ )
c
c  -- LAPACK auxiliary routine (version 3.4.2) --
c  -- LAPACK is a software package provided by Univ. of Tennessee,    --
c  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
c     September 2012
c
c     .. Scalar Arguments ..
      INTEGER            INCX, N
      DOUBLE PRECISION   SCALE, SUMSQ
c     ..
c     .. Array Arguments ..
      DOUBLE PRECISION   X( * )
c     ..
c
c =====================================================================
c
c     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
c     ..
c     .. Local Scalars ..
      INTEGER            IX
      DOUBLE PRECISION   ABSXI
c     ..
c     .. External Functions ..
      LOGICAL            DISNAN
      EXTERNAL           DISNAN
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          ABS
c     ..
c     .. Executable Statements ..
c
      IF( N.GT.0 ) THEN
         DO 10 IX = 1, 1 + ( N-1 )*INCX, INCX
            ABSXI = ABS( X( IX ) )
            IF( ABSXI.GT.ZERO.OR.DISNAN( ABSXI ) ) THEN
               IF( SCALE.LT.ABSXI ) THEN
                  SUMSQ = 1 + SUMSQ*( SCALE / ABSXI )**2
                  SCALE = ABSXI
               ELSE
                  SUMSQ = SUMSQ + ( ABSXI / SCALE )**2
               END IF
            END IF
   10    CONTINUE
      END IF
      RETURN
c
c     End of DLASSQ
c
      END

      SUBROUTINE DLASRT( ID, N, D, INFO )
c
c  -- LAPACK computational routine (version 3.4.2) --
c  -- LAPACK is a software package provided by Univ. of Tennessee,    --
c  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
c     September 2012
c
c     .. Scalar Arguments ..
      CHARACTER          ID
      INTEGER            INFO, N
c     ..
c     .. Array Arguments ..
      DOUBLE PRECISION   D( * )
c     ..
c
c  =====================================================================
c
c     .. Parameters ..
      INTEGER            SELECT
      PARAMETER          ( SELECT = 20 )
c     ..
c     .. Local Scalars ..
      INTEGER            DIR, ENDD, I, J, START, STKPNT
      DOUBLE PRECISION   D1, D2, D3, DMNMX, TMP
c     ..
c     .. Local Arrays ..
      INTEGER            STACK( 2, 32 )
c     ..
c     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
c     ..
c     .. External Subroutines ..
      EXTERNAL           XERBLA
c     ..
c     .. Executable Statements ..
c
c     Test the input paramters.
c
      INFO = 0
      DIR = -1
      IF( LSAME( ID, 'D' ) ) THEN
         DIR = 0
      ELSE IF( LSAME( ID, 'I' ) ) THEN
         DIR = 1
      END IF
      IF( DIR.EQ.-1 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLASRT', -INFO )
         RETURN
      END IF
c
c     Quick return if possible
c
      IF( N.LE.1 )
     $   RETURN
c
      STKPNT = 1
      STACK( 1, 1 ) = 1
      STACK( 2, 1 ) = N
   10 CONTINUE
      START = STACK( 1, STKPNT )
      ENDD = STACK( 2, STKPNT )
      STKPNT = STKPNT - 1
      IF( ENDD-START.LE.SELECT .AND. ENDD-START.GT.0 ) THEN
c
c        Do Insertion sort on D( START:ENDD )
c
         IF( DIR.EQ.0 ) THEN
c
c           Sort into decreasing order
c
            DO 30 I = START + 1, ENDD
               DO 20 J = I, START + 1, -1
                  IF( D( J ).GT.D( J-1 ) ) THEN
                     DMNMX = D( J )
                     D( J ) = D( J-1 )
                     D( J-1 ) = DMNMX
                  ELSE
                     GO TO 30
                  END IF
   20          CONTINUE
   30       CONTINUE
c
         ELSE
c
c           Sort into increasing order
c
            DO 50 I = START + 1, ENDD
               DO 40 J = I, START + 1, -1
                  IF( D( J ).LT.D( J-1 ) ) THEN
                     DMNMX = D( J )
                     D( J ) = D( J-1 )
                     D( J-1 ) = DMNMX
                  ELSE
                     GO TO 50
                  END IF
   40          CONTINUE
   50       CONTINUE
c
         END IF
c
      ELSE IF( ENDD-START.GT.SELECT ) THEN
c
c        Partition D( START:ENDD ) and stack parts, largest one first
c
c        Choose partition entry as median of 3
c
         D1 = D( START )
         D2 = D( ENDD )
         I = ( START+ENDD ) / 2
         D3 = D( I )
         IF( D1.LT.D2 ) THEN
            IF( D3.LT.D1 ) THEN
               DMNMX = D1
            ELSE IF( D3.LT.D2 ) THEN
               DMNMX = D3
            ELSE
               DMNMX = D2
            END IF
         ELSE
            IF( D3.LT.D2 ) THEN
               DMNMX = D2
            ELSE IF( D3.LT.D1 ) THEN
               DMNMX = D3
            ELSE
               DMNMX = D1
            END IF
         END IF
c
         IF( DIR.EQ.0 ) THEN
c
c           Sort into decreasing order
c
            I = START - 1
            J = ENDD + 1
   60       CONTINUE
   70       CONTINUE
            J = J - 1
            IF( D( J ).LT.DMNMX )
     $         GO TO 70
   80       CONTINUE
            I = I + 1
            IF( D( I ).GT.DMNMX )
     $         GO TO 80
            IF( I.LT.J ) THEN
               TMP = D( I )
               D( I ) = D( J )
               D( J ) = TMP
               GO TO 60
            END IF
            IF( J-START.GT.ENDD-J-1 ) THEN
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
            ELSE
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
            END IF
         ELSE
c
c           Sort into increasing order
c
            I = START - 1
            J = ENDD + 1
   90       CONTINUE
  100       CONTINUE
            J = J - 1
            IF( D( J ).GT.DMNMX )
     $         GO TO 100
  110       CONTINUE
            I = I + 1
            IF( D( I ).LT.DMNMX )
     $         GO TO 110
            IF( I.LT.J ) THEN
               TMP = D( I )
               D( I ) = D( J )
               D( J ) = TMP
               GO TO 90
            END IF
            IF( J-START.GT.ENDD-J-1 ) THEN
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
            ELSE
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
            END IF
         END IF
      END IF
      IF( STKPNT.GT.0 )
     $   GO TO 10
      RETURN
c
c     End of DLASRT
c
      END

      SUBROUTINE DLASET( UPLO, M, N, ALPHA, BETA, A, LDA )
c
c  -- LAPACK auxiliary routine (version 3.4.2) --
c  -- LAPACK is a software package provided by Univ. of Tennessee,    --
c  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
c     September 2012
c
c     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, M, N
      DOUBLE PRECISION   ALPHA, BETA
c     ..
c     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
c     ..
c
c =====================================================================
c
c     .. Local Scalars ..
      INTEGER            I, J
c     ..
c     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          MIN
c     ..
c     .. Executable Statements ..
c
      IF( LSAME( UPLO, 'U' ) ) THEN
c
c        Set the strictly upper triangular or trapezoidal part of the
c        array to ALPHA.
c
         DO 20 J = 2, N
            DO 10 I = 1, MIN( J-1, M )
               A( I, J ) = ALPHA
   10       CONTINUE
   20    CONTINUE
c
      ELSE IF( LSAME( UPLO, 'L' ) ) THEN
c
c        Set the strictly lower triangular or trapezoidal part of the
c        array to ALPHA.
c
         DO 40 J = 1, MIN( M, N )
            DO 30 I = J + 1, M
               A( I, J ) = ALPHA
   30       CONTINUE
   40    CONTINUE
c
      ELSE
c
c        Set the leading m-by-n submatrix to ALPHA.
c
         DO 60 J = 1, N
            DO 50 I = 1, M
               A( I, J ) = ALPHA
   50       CONTINUE
   60    CONTINUE
      END IF
c
c     Set the first min(M,N) diagonal elements to BETA.
c
      DO 70 I = 1, MIN( M, N )
         A( I, I ) = BETA
   70 CONTINUE
c
      RETURN
c
c     End of DLASET
c
      END

      SUBROUTINE DLARFB( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV,
     $                   T, LDT, C, LDC, WORK, LDWORK )
c
c  -- LAPACK auxiliary routine (version 3.5.0) --
c  -- LAPACK is a software package provided by Univ. of Tennessee,    --
c  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
c     June 2013
c
c     .. Scalar Arguments ..
      CHARACTER          DIRECT, SIDE, STOREV, TRANS
      INTEGER            K, LDC, LDT, LDV, LDWORK, M, N
c     ..
c     .. Array Arguments ..
      DOUBLE PRECISION   C( LDC, * ), T( LDT, * ), V( LDV, * ),
     $                   WORK( LDWORK, * )
c     ..
c
c  =====================================================================
c
c     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
c     ..
c     .. Local Scalars ..
      CHARACTER          TRANST
      INTEGER            I, J
c     ..
c     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
c     ..
c     .. External Subroutines ..
      EXTERNAL           DCOPY, DGEMM, DTRMM
c     ..
c     .. Executable Statements ..
c
c     Quick return if possible
c
      IF( M.LE.0 .OR. N.LE.0 )
     $   RETURN
c
      IF( LSAME( TRANS, 'N' ) ) THEN
         TRANST = 'T'
      ELSE
         TRANST = 'N'
      END IF
c
      IF( LSAME( STOREV, 'C' ) ) THEN
c
         IF( LSAME( DIRECT, 'F' ) ) THEN
c
c           Let  V =  ( V1 )    (first K rows)
c                     ( V2 )
c           where  V1  is unit lower triangular.
c
            IF( LSAME( SIDE, 'L' ) ) THEN
c
c              Form  H * C  or  H**T * C  where  C = ( C1 )
c                                                    ( C2 )
c
c              W := C**T * V  =  (C1**T * V1 + C2**T * V2)  (stored in WORK)
c
c              W := C1**T
c
               DO 10 J = 1, K
                  CALL DCOPY( N, C( J, 1 ), LDC, WORK( 1, J ), 1 )
   10          CONTINUE
c
c              W := W * V1
c
               CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', N,
     $                     K, ONE, V, LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
c
c                 W := W + C2**T * V2
c
                  CALL DGEMM( 'Transpose', 'No transpose', N, K, M-K,
     $                        ONE, C( K+1, 1 ), LDC, V( K+1, 1 ), LDV,
     $                        ONE, WORK, LDWORK )
               END IF
c
c              W := W * T**T  or  W * T
c
               CALL DTRMM( 'Right', 'Upper', TRANST, 'Non-unit', N, K,
     $                     ONE, T, LDT, WORK, LDWORK )
c
c              C := C - V * W**T
c
               IF( M.GT.K ) THEN
c
c                 C2 := C2 - V2 * W**T
c
                  CALL DGEMM( 'No transpose', 'Transpose', M-K, N, K,
     $                        -ONE, V( K+1, 1 ), LDV, WORK, LDWORK, ONE,
     $                        C( K+1, 1 ), LDC )
               END IF
c
c              W := W * V1**T
c
               CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', N, K,
     $                     ONE, V, LDV, WORK, LDWORK )
c
c              C1 := C1 - W**T
c
               DO 30 J = 1, K
                  DO 20 I = 1, N
                     C( J, I ) = C( J, I ) - WORK( I, J )
   20             CONTINUE
   30          CONTINUE
c
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
c
c              Form  C * H  or  C * H**T  where  C = ( C1  C2 )
c
c              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
c
c              W := C1
c
               DO 40 J = 1, K
                  CALL DCOPY( M, C( 1, J ), 1, WORK( 1, J ), 1 )
   40          CONTINUE
c
c              W := W * V1
c
               CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', M,
     $                     K, ONE, V, LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
c
c                 W := W + C2 * V2
c
                  CALL DGEMM( 'No transpose', 'No transpose', M, K, N-K,
     $                        ONE, C( 1, K+1 ), LDC, V( K+1, 1 ), LDV,
     $                        ONE, WORK, LDWORK )
               END IF
c
c              W := W * T  or  W * T**T
c
               CALL DTRMM( 'Right', 'Upper', TRANS, 'Non-unit', M, K,
     $                     ONE, T, LDT, WORK, LDWORK )
c
c              C := C - W * V**T
c
               IF( N.GT.K ) THEN
c
c                 C2 := C2 - W * V2**T
c
                  CALL DGEMM( 'No transpose', 'Transpose', M, N-K, K,
     $                        -ONE, WORK, LDWORK, V( K+1, 1 ), LDV, ONE,
     $                        C( 1, K+1 ), LDC )
               END IF
c
c              W := W * V1**T
c
               CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', M, K,
     $                     ONE, V, LDV, WORK, LDWORK )
c
c              C1 := C1 - W
c
               DO 60 J = 1, K
                  DO 50 I = 1, M
                     C( I, J ) = C( I, J ) - WORK( I, J )
   50             CONTINUE
   60          CONTINUE
            END IF
c
         ELSE
c
c           Let  V =  ( V1 )
c                     ( V2 )    (last K rows)
c           where  V2  is unit upper triangular.
c
            IF( LSAME( SIDE, 'L' ) ) THEN
c
c              Form  H * C  or  H**T * C  where  C = ( C1 )
c                                                    ( C2 )
c
c              W := C**T * V  =  (C1**T * V1 + C2**T * V2)  (stored in WORK)
c
c              W := C2**T
c
               DO 70 J = 1, K
                  CALL DCOPY( N, C( M-K+J, 1 ), LDC, WORK( 1, J ), 1 )
   70          CONTINUE
c
c              W := W * V2
c
               CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', N,
     $                     K, ONE, V( M-K+1, 1 ), LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
c
c                 W := W + C1**T * V1
c
                  CALL DGEMM( 'Transpose', 'No transpose', N, K, M-K,
     $                        ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
c
c              W := W * T**T  or  W * T
c
               CALL DTRMM( 'Right', 'Lower', TRANST, 'Non-unit', N, K,
     $                     ONE, T, LDT, WORK, LDWORK )
c
c              C := C - V * W**T
c
               IF( M.GT.K ) THEN
c
c                 C1 := C1 - V1 * W**T
c
                  CALL DGEMM( 'No transpose', 'Transpose', M-K, N, K,
     $                        -ONE, V, LDV, WORK, LDWORK, ONE, C, LDC )
               END IF
c
c              W := W * V2**T
c
               CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', N, K,
     $                     ONE, V( M-K+1, 1 ), LDV, WORK, LDWORK )
c
c              C2 := C2 - W**T
c
               DO 90 J = 1, K
                  DO 80 I = 1, N
                     C( M-K+J, I ) = C( M-K+J, I ) - WORK( I, J )
   80             CONTINUE
   90          CONTINUE
c
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
c
c              Form  C * H  or  C * H**T  where  C = ( C1  C2 )
c
c              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
c
c              W := C2
c
               DO 100 J = 1, K
                  CALL DCOPY( M, C( 1, N-K+J ), 1, WORK( 1, J ), 1 )
  100          CONTINUE
c
c              W := W * V2
c
               CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', M,
     $                     K, ONE, V( N-K+1, 1 ), LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
c
c                 W := W + C1 * V1
c
                  CALL DGEMM( 'No transpose', 'No transpose', M, K, N-K,
     $                        ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
c
c              W := W * T  or  W * T**T
c
               CALL DTRMM( 'Right', 'Lower', TRANS, 'Non-unit', M, K,
     $                     ONE, T, LDT, WORK, LDWORK )
c
c              C := C - W * V**T
c
               IF( N.GT.K ) THEN
c
c                 C1 := C1 - W * V1**T
c
                  CALL DGEMM( 'No transpose', 'Transpose', M, N-K, K,
     $                        -ONE, WORK, LDWORK, V, LDV, ONE, C, LDC )
               END IF
c
c              W := W * V2**T
c
               CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', M, K,
     $                     ONE, V( N-K+1, 1 ), LDV, WORK, LDWORK )
c
c              C2 := C2 - W
c
               DO 120 J = 1, K
                  DO 110 I = 1, M
                     C( I, N-K+J ) = C( I, N-K+J ) - WORK( I, J )
  110             CONTINUE
  120          CONTINUE
            END IF
         END IF
c
      ELSE IF( LSAME( STOREV, 'R' ) ) THEN
c
         IF( LSAME( DIRECT, 'F' ) ) THEN
c
c           Let  V =  ( V1  V2 )    (V1: first K columns)
c           where  V1  is unit upper triangular.
c
            IF( LSAME( SIDE, 'L' ) ) THEN
c
c              Form  H * C  or  H**T * C  where  C = ( C1 )
c                                                    ( C2 )
c
c              W := C**T * V**T  =  (C1**T * V1**T + C2**T * V2**T) (stored in WORK)
c
c              W := C1**T
c
               DO 130 J = 1, K
                  CALL DCOPY( N, C( J, 1 ), LDC, WORK( 1, J ), 1 )
  130          CONTINUE
c
c              W := W * V1**T
c
               CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', N, K,
     $                     ONE, V, LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
c
c                 W := W + C2**T * V2**T
c
                  CALL DGEMM( 'Transpose', 'Transpose', N, K, M-K, ONE,
     $                        C( K+1, 1 ), LDC, V( 1, K+1 ), LDV, ONE,
     $                        WORK, LDWORK )
               END IF
c
c              W := W * T**T  or  W * T
c
               CALL DTRMM( 'Right', 'Upper', TRANST, 'Non-unit', N, K,
     $                     ONE, T, LDT, WORK, LDWORK )
c
c              C := C - V**T * W**T
c
               IF( M.GT.K ) THEN
c
c                 C2 := C2 - V2**T * W**T
c
                  CALL DGEMM( 'Transpose', 'Transpose', M-K, N, K, -ONE,
     $                        V( 1, K+1 ), LDV, WORK, LDWORK, ONE,
     $                        C( K+1, 1 ), LDC )
               END IF
c
c              W := W * V1
c
               CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', N,
     $                     K, ONE, V, LDV, WORK, LDWORK )
c
c              C1 := C1 - W**T
c
               DO 150 J = 1, K
                  DO 140 I = 1, N
                     C( J, I ) = C( J, I ) - WORK( I, J )
  140             CONTINUE
  150          CONTINUE
c
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
c
c              Form  C * H  or  C * H**T  where  C = ( C1  C2 )
c
c              W := C * V**T  =  (C1*V1**T + C2*V2**T)  (stored in WORK)
c
c              W := C1
c
               DO 160 J = 1, K
                  CALL DCOPY( M, C( 1, J ), 1, WORK( 1, J ), 1 )
  160          CONTINUE
c
c              W := W * V1**T
c
               CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', M, K,
     $                     ONE, V, LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
c
c                 W := W + C2 * V2**T
c
                  CALL DGEMM( 'No transpose', 'Transpose', M, K, N-K,
     $                        ONE, C( 1, K+1 ), LDC, V( 1, K+1 ), LDV,
     $                        ONE, WORK, LDWORK )
               END IF
c
c              W := W * T  or  W * T**T
c
               CALL DTRMM( 'Right', 'Upper', TRANS, 'Non-unit', M, K,
     $                     ONE, T, LDT, WORK, LDWORK )
c
c              C := C - W * V
c
               IF( N.GT.K ) THEN
c
c                 C2 := C2 - W * V2
c
                  CALL DGEMM( 'No transpose', 'No transpose', M, N-K, K,
     $                        -ONE, WORK, LDWORK, V( 1, K+1 ), LDV, ONE,
     $                        C( 1, K+1 ), LDC )
               END IF
c
c              W := W * V1
c
               CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', M,
     $                     K, ONE, V, LDV, WORK, LDWORK )
c
c              C1 := C1 - W
c
               DO 180 J = 1, K
                  DO 170 I = 1, M
                     C( I, J ) = C( I, J ) - WORK( I, J )
  170             CONTINUE
  180          CONTINUE
c
            END IF
c
         ELSE
c
c           Let  V =  ( V1  V2 )    (V2: last K columns)
c           where  V2  is unit lower triangular.
c
            IF( LSAME( SIDE, 'L' ) ) THEN
c
c              Form  H * C  or  H**T * C  where  C = ( C1 )
c                                                    ( C2 )
c
c              W := C**T * V**T  =  (C1**T * V1**T + C2**T * V2**T) (stored in WORK)
c
c              W := C2**T
c
               DO 190 J = 1, K
                  CALL DCOPY( N, C( M-K+J, 1 ), LDC, WORK( 1, J ), 1 )
  190          CONTINUE
c
c              W := W * V2**T
c
               CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', N, K,
     $                     ONE, V( 1, M-K+1 ), LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
c
c                 W := W + C1**T * V1**T
c
                  CALL DGEMM( 'Transpose', 'Transpose', N, K, M-K, ONE,
     $                        C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
c
c              W := W * T**T  or  W * T
c
               CALL DTRMM( 'Right', 'Lower', TRANST, 'Non-unit', N, K,
     $                     ONE, T, LDT, WORK, LDWORK )
c
c              C := C - V**T * W**T
c
               IF( M.GT.K ) THEN
c
c                 C1 := C1 - V1**T * W**T
c
                  CALL DGEMM( 'Transpose', 'Transpose', M-K, N, K, -ONE,
     $                        V, LDV, WORK, LDWORK, ONE, C, LDC )
               END IF
c
c              W := W * V2
c
               CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', N,
     $                     K, ONE, V( 1, M-K+1 ), LDV, WORK, LDWORK )
c
c              C2 := C2 - W**T
c
               DO 210 J = 1, K
                  DO 200 I = 1, N
                     C( M-K+J, I ) = C( M-K+J, I ) - WORK( I, J )
  200             CONTINUE
  210          CONTINUE
c
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
c
c              Form  C * H  or  C * H'  where  C = ( C1  C2 )
c
c              W := C * V**T  =  (C1*V1**T + C2*V2**T)  (stored in WORK)
c
c              W := C2
c
               DO 220 J = 1, K
                  CALL DCOPY( M, C( 1, N-K+J ), 1, WORK( 1, J ), 1 )
  220          CONTINUE
c
c              W := W * V2**T
c
               CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', M, K,
     $                     ONE, V( 1, N-K+1 ), LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
c
c                 W := W + C1 * V1**T
c
                  CALL DGEMM( 'No transpose', 'Transpose', M, K, N-K,
     $                        ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
c
c              W := W * T  or  W * T**T
c
               CALL DTRMM( 'Right', 'Lower', TRANS, 'Non-unit', M, K,
     $                     ONE, T, LDT, WORK, LDWORK )
c
c              C := C - W * V
c
               IF( N.GT.K ) THEN
c
c                 C1 := C1 - W * V1
c
                  CALL DGEMM( 'No transpose', 'No transpose', M, N-K, K,
     $                        -ONE, WORK, LDWORK, V, LDV, ONE, C, LDC )
               END IF
c
c              W := W * V2
c
               CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', M,
     $                     K, ONE, V( 1, N-K+1 ), LDV, WORK, LDWORK )
c
c              C1 := C1 - W
c
               DO 240 J = 1, K
                  DO 230 I = 1, M
                     C( I, N-K+J ) = C( I, N-K+J ) - WORK( I, J )
  230             CONTINUE
  240          CONTINUE
c
            END IF
c
         END IF
      END IF
c
      RETURN
c
c     End of DLARFB
c
      END

      SUBROUTINE DLARF( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
c
c  -- LAPACK auxiliary routine (version 3.4.2) --
c  -- LAPACK is a software package provided by Univ. of Tennessee,    --
c  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
c     September 2012
c
c     .. Scalar Arguments ..
      CHARACTER          SIDE
      INTEGER            INCV, LDC, M, N
      DOUBLE PRECISION   TAU
c     ..
c     .. Array Arguments ..
      DOUBLE PRECISION   C( LDC, * ), V( * ), WORK( * )
c     ..
c
c  =====================================================================
c
c     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
c     ..
c     .. Local Scalars ..
      LOGICAL            APPLYLEFT
      INTEGER            I, LASTV, LASTC
c     ..
c     .. External Subroutines ..
      EXTERNAL           DGEMV, DGER
c     ..
c     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILADLR, ILADLC
      EXTERNAL           LSAME, ILADLR, ILADLC
c     ..
c     .. Executable Statements ..
c
      APPLYLEFT = LSAME( SIDE, 'L' )
      LASTV = 0
      LASTC = 0
      IF( TAU.NE.ZERO ) THEN
!     Set up variables for scanning V.  LASTV begins pointing to the end
!     of V.
         IF( APPLYLEFT ) THEN
            LASTV = M
         ELSE
            LASTV = N
         END IF
         IF( INCV.GT.0 ) THEN
            I = 1 + (LASTV-1) * INCV
         ELSE
            I = 1
         END IF
!     Look for the last non-zero row in V.
         DO WHILE( LASTV.GT.0 .AND. V( I ).EQ.ZERO )
            LASTV = LASTV - 1
            I = I - INCV
         END DO
         IF( APPLYLEFT ) THEN
!     Scan for the last non-zero column in C(1:lastv,:).
            LASTC = ILADLC(LASTV, N, C, LDC)
         ELSE
!     Scan for the last non-zero row in C(:,1:lastv).
            LASTC = ILADLR(M, LASTV, C, LDC)
         END IF
      END IF
!     Note that lastc.eq.0 renders the BLAS operations null; no special
!     case is needed at this level.
      IF( APPLYLEFT ) THEN
c
c        Form  H * C
c
         IF( LASTV.GT.0 ) THEN
c
c           w(1:lastc,1) := C(1:lastv,1:lastc)**T * v(1:lastv,1)
c
            CALL DGEMV( 'Transpose', LASTV, LASTC, ONE, C, LDC, V, INCV,
     $           ZERO, WORK, 1 )
c
c           C(1:lastv,1:lastc) := C(...) - v(1:lastv,1) * w(1:lastc,1)**T
c
            CALL DGER( LASTV, LASTC, -TAU, V, INCV, WORK, 1, C, LDC )
         END IF
      ELSE
c
c        Form  C * H
c
         IF( LASTV.GT.0 ) THEN
c
c           w(1:lastc,1) := C(1:lastc,1:lastv) * v(1:lastv,1)
c
            CALL DGEMV( 'No transpose', LASTC, LASTV, ONE, C, LDC,
     $           V, INCV, ZERO, WORK, 1 )
c
c           C(1:lastc,1:lastv) := C(...) - w(1:lastc,1) * v(1:lastv,1)**T
c
            CALL DGER( LASTC, LASTV, -TAU, WORK, 1, V, INCV, C, LDC )
         END IF
      END IF
      RETURN
c
c     End of DLARF
c
      END

      SUBROUTINE DLARTG( F, G, CS, SN, R )
c
c  -- LAPACK auxiliary routine (version 3.4.2) --
c  -- LAPACK is a software package provided by Univ. of Tennessee,    --
c  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
c     September 2012
c
c     .. Scalar Arguments ..
      DOUBLE PRECISION   CS, F, G, R, SN
c     ..
c
c  =====================================================================
c
c     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D0 )
c     ..
c     .. Local Scalars ..
c     LOGICAL            FIRST
      INTEGER            COUNT, I
      DOUBLE PRECISION   EPS, F1, G1, SAFMIN, SAFMN2, SAFMX2, SCALE
c     ..
c     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          ABS, INT, LOG, MAX, SQRT
c     ..
c     .. Save statement ..
c     SAVE               FIRST, SAFMX2, SAFMIN, SAFMN2
c     ..
c     .. Data statements ..
c     DATA               FIRST / .TRUE. /
c     ..
c     .. Executable Statements ..
c
c     IF( FIRST ) THEN
         SAFMIN = DLAMCH( 'S' )
         EPS = DLAMCH( 'E' )
         SAFMN2 = DLAMCH( 'B' )**INT( LOG( SAFMIN / EPS ) /
     $            LOG( DLAMCH( 'B' ) ) / TWO )
         SAFMX2 = ONE / SAFMN2
c        FIRST = .FALSE.
c     END IF
      IF( G.EQ.ZERO ) THEN
         CS = ONE
         SN = ZERO
         R = F
      ELSE IF( F.EQ.ZERO ) THEN
         CS = ZERO
         SN = ONE
         R = G
      ELSE
         F1 = F
         G1 = G
         SCALE = MAX( ABS( F1 ), ABS( G1 ) )
         IF( SCALE.GE.SAFMX2 ) THEN
            COUNT = 0
   10       CONTINUE
            COUNT = COUNT + 1
            F1 = F1*SAFMN2
            G1 = G1*SAFMN2
            SCALE = MAX( ABS( F1 ), ABS( G1 ) )
            IF( SCALE.GE.SAFMX2 )
     $         GO TO 10
            R = SQRT( F1**2+G1**2 )
            CS = F1 / R
            SN = G1 / R
            DO 20 I = 1, COUNT
               R = R*SAFMX2
   20       CONTINUE
         ELSE IF( SCALE.LE.SAFMN2 ) THEN
            COUNT = 0
   30       CONTINUE
            COUNT = COUNT + 1
            F1 = F1*SAFMX2
            G1 = G1*SAFMX2
            SCALE = MAX( ABS( F1 ), ABS( G1 ) )
            IF( SCALE.LE.SAFMN2 )
     $         GO TO 30
            R = SQRT( F1**2+G1**2 )
            CS = F1 / R
            SN = G1 / R
            DO 40 I = 1, COUNT
               R = R*SAFMN2
   40       CONTINUE
         ELSE
            R = SQRT( F1**2+G1**2 )
            CS = F1 / R
            SN = G1 / R
         END IF
         IF( ABS( F ).GT.ABS( G ) .AND. CS.LT.ZERO ) THEN
            CS = -CS
            SN = -SN
            R = -R
         END IF
      END IF
      RETURN
c
c     End of DLARTG
c
      END

      LOGICAL FUNCTION DISNAN( DIN )
c
c  -- LAPACK auxiliary routine (version 3.4.2) --
c  -- LAPACK is a software package provided by Univ. of Tennessee,    --
c  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
c     September 2012
c
c     .. Scalar Arguments ..
      DOUBLE PRECISION   DIN
c     ..
c
c  =====================================================================
c
c  .. External Functions ..
      LOGICAL DLAISNAN
      EXTERNAL DLAISNAN
c  ..
c  .. Executable Statements ..
      DISNAN = DLAISNAN(DIN,DIN)
      RETURN
      END

      INTEGER FUNCTION ILADLC( M, N, A, LDA )
c
c  -- LAPACK auxiliary routine (version 3.4.2) --
c  -- LAPACK is a software package provided by Univ. of Tennessee,    --
c  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
c     September 2012
c
c     .. Scalar Arguments ..
      INTEGER            M, N, LDA
c     ..
c     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
c     ..
c
c  =====================================================================
c
c     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER ( ZERO = 0.0D+0 )
c     ..
c     .. Local Scalars ..
      INTEGER I
c     ..
c     .. Executable Statements ..
c
c     Quick test for the common case where one corner is non-zero.
      IF( N.EQ.0 ) THEN
         ILADLC = N
      ELSE IF( A(1, N).NE.ZERO .OR. A(M, N).NE.ZERO ) THEN
         ILADLC = N
      ELSE
c     Now scan each column from the end, returning with the first non-zero.
         DO ILADLC = N, 1, -1
            DO I = 1, M
               IF( A(I, ILADLC).NE.ZERO ) RETURN
            END DO
         END DO
      END IF
      RETURN
      END

      INTEGER FUNCTION ILADLR( M, N, A, LDA )
c
c  -- LAPACK auxiliary routine (version 3.4.2) --
c  -- LAPACK is a software package provided by Univ. of Tennessee,    --
c  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
c     September 2012
c
c     .. Scalar Arguments ..
      INTEGER            M, N, LDA
c     ..
c     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
c     ..
c
c  =====================================================================
c
c     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER ( ZERO = 0.0D+0 )
c     ..
c     .. Local Scalars ..
      INTEGER I, J
c     ..
c     .. Executable Statements ..
c
c     Quick test for the common case where one corner is non-zero.
      IF( M.EQ.0 ) THEN
         ILADLR = M
      ELSE IF( A(M, 1).NE.ZERO .OR. A(M, N).NE.ZERO ) THEN
         ILADLR = M
      ELSE
c     Scan up each column tracking the last zero row seen.
         ILADLR = 0
         DO J = 1, N
            I=M
            DO WHILE((A(MAX(I,1),J).EQ.ZERO).AND.(I.GE.1))
               I=I-1
            ENDDO
            ILADLR = MAX( ILADLR, I )
         END DO
      END IF
      RETURN
      END

      LOGICAL FUNCTION DLAISNAN( DIN1, DIN2 )
c
c  -- LAPACK auxiliary routine (version 3.4.2) --
c  -- LAPACK is a software package provided by Univ. of Tennessee,    --
c  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
c     September 2012
c
c     .. Scalar Arguments ..
      DOUBLE PRECISION   DIN1, DIN2
c     ..
c
c  =====================================================================
c
c  .. Executable Statements ..
      DLAISNAN = (DIN1.NE.DIN2)
      RETURN
      END

      DOUBLE PRECISION FUNCTION DLAMCH( CMACH )
c
c  -- LAPACK auxiliary routine (version 3.4.0) --
c  -- LAPACK is a software package provided by Univ. of Tennessee,    --
c  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
c     November 2011
c
c     .. Scalar Arguments ..
      CHARACTER          CMACH
c     ..
c
c     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B
c     ..
c
c =====================================================================
c
c     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
c     ..
c     .. Local Scalars ..
      DOUBLE PRECISION   RND, EPS, SFMIN, SMALL, RMACH
c     ..
c     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          DIGITS, EPSILON, HUGE, MAXEXPONENT,
     $                   MINEXPONENT, RADIX, TINY
c     ..
c     .. Executable Statements ..
c
c
c     Assume rounding, not chopping. Always.
c
      RND = ONE
c
      IF( ONE.EQ.RND ) THEN
         EPS = EPSILON(ZERO) * 0.5
      ELSE
         EPS = EPSILON(ZERO)
      END IF
c
      IF( LSAME( CMACH, 'E' ) ) THEN
         RMACH = EPS
      ELSE IF( LSAME( CMACH, 'S' ) ) THEN
         SFMIN = TINY(ZERO)
         SMALL = ONE / HUGE(ZERO)
         IF( SMALL.GE.SFMIN ) THEN
c
c           Use SMALL plus a bit, to avoid the possibility of rounding
c           causing overflow when computing  1/sfmin.
c
            SFMIN = SMALL*( ONE+EPS )
         END IF
         RMACH = SFMIN
      ELSE IF( LSAME( CMACH, 'B' ) ) THEN
         RMACH = RADIX(ZERO)
      ELSE IF( LSAME( CMACH, 'P' ) ) THEN
         RMACH = EPS * RADIX(ZERO)
      ELSE IF( LSAME( CMACH, 'N' ) ) THEN
         RMACH = DIGITS(ZERO)
      ELSE IF( LSAME( CMACH, 'R' ) ) THEN
         RMACH = RND
      ELSE IF( LSAME( CMACH, 'M' ) ) THEN
         RMACH = MINEXPONENT(ZERO)
      ELSE IF( LSAME( CMACH, 'U' ) ) THEN
         RMACH = tiny(zero)
      ELSE IF( LSAME( CMACH, 'L' ) ) THEN
         RMACH = MAXEXPONENT(ZERO)
      ELSE IF( LSAME( CMACH, 'O' ) ) THEN
         RMACH = HUGE(ZERO)
      ELSE
         RMACH = ZERO
      END IF
c
      DLAMCH = RMACH
      RETURN
c
c     End of DLAMCH
c
      END

      INTEGER FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
c
c  -- LAPACK auxiliary routine (version 3.4.0) --
c  -- LAPACK is a software package provided by Univ. of Tennessee,    --
c  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
c     November 2011
c
c     .. Scalar Arguments ..
      CHARACTER*( * )    NAME, OPTS
      INTEGER            ISPEC, N1, N2, N3, N4
c     ..
c
c  =====================================================================
c
c     .. Local Scalars ..
      INTEGER            I, IC, IZ, NB, NBMIN, NX
      LOGICAL            CNAME, SNAME
      CHARACTER          C1*1, C2*2, C4*2, C3*3, SUBNAM*6
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          CHAR, ICHAR, INT, MIN, REAL
c     ..
c     .. External Functions ..
      INTEGER            IEEECK, IPARMQ
      EXTERNAL           IEEECK, IPARMQ
c     ..
c     .. Executable Statements ..
c
      GO TO ( 10, 10, 10, 80, 90, 100, 110, 120,
     $        130, 140, 150, 160, 160, 160, 160, 160 )ISPEC
c
c     Invalid value for ISPEC
c
      ILAENV = -1
      RETURN
c
   10 CONTINUE
c
c     Convert NAME to upper case if the first character is lower case.
c
      ILAENV = 1
      SUBNAM = NAME
      IC = ICHAR( SUBNAM( 1: 1 ) )
      IZ = ICHAR( 'Z' )
      IF( IZ.EQ.90 .OR. IZ.EQ.122 ) THEN
c
c        ASCII character set
c
         IF( IC.GE.97 .AND. IC.LE.122 ) THEN
            SUBNAM( 1: 1 ) = CHAR( IC-32 )
            DO 20 I = 2, 6
               IC = ICHAR( SUBNAM( I: I ) )
               IF( IC.GE.97 .AND. IC.LE.122 )
     $            SUBNAM( I: I ) = CHAR( IC-32 )
   20       CONTINUE
         END IF
c
      ELSE IF( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN
c
c        EBCDIC character set
c
         IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $       ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $       ( IC.GE.162 .AND. IC.LE.169 ) ) THEN
            SUBNAM( 1: 1 ) = CHAR( IC+64 )
            DO 30 I = 2, 6
               IC = ICHAR( SUBNAM( I: I ) )
               IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $             ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $             ( IC.GE.162 .AND. IC.LE.169 ) )SUBNAM( I:
     $             I ) = CHAR( IC+64 )
   30       CONTINUE
         END IF
c
      ELSE IF( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN
c
c        Prime machines:  ASCII+128
c
         IF( IC.GE.225 .AND. IC.LE.250 ) THEN
            SUBNAM( 1: 1 ) = CHAR( IC-32 )
            DO 40 I = 2, 6
               IC = ICHAR( SUBNAM( I: I ) )
               IF( IC.GE.225 .AND. IC.LE.250 )
     $            SUBNAM( I: I ) = CHAR( IC-32 )
   40       CONTINUE
         END IF
      END IF
c
      C1 = SUBNAM( 1: 1 )
      SNAME = C1.EQ.'S' .OR. C1.EQ.'D'
      CNAME = C1.EQ.'C' .OR. C1.EQ.'Z'
      IF( .NOT.( CNAME .OR. SNAME ) )
     $   RETURN
      C2 = SUBNAM( 2: 3 )
      C3 = SUBNAM( 4: 6 )
      C4 = C3( 2: 3 )
c
      GO TO ( 50, 60, 70 )ISPEC
c
   50 CONTINUE
c
c     ISPEC = 1:  block size
c
c     In these examples, separate code is provided for setting NB for
c     real and complex.  We assume that NB will take the same value in
c     single or double precision.
c
      NB = 1
c
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.
     $            C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'PO' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NB = 32
         ELSE IF( SNAME .AND. C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            NB = 64
         ELSE IF( C3.EQ.'TRD' ) THEN
            NB = 32
         ELSE IF( C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NB = 32
            END IF
         ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NB = 32
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NB = 32
            END IF
         ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NB = 32
            END IF
         END IF
      ELSE IF( C2.EQ.'GB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'PB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'TR' ) THEN
         IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'LA' ) THEN
         IF( C3.EQ.'UUM' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'ST' ) THEN
         IF( C3.EQ.'EBZ' ) THEN
            NB = 1
         END IF
      END IF
      ILAENV = NB
      RETURN
c
   60 CONTINUE
c
c     ISPEC = 2:  minimum block size
c
      NBMIN = 2
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. C3.EQ.
     $       'QLF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 8
            ELSE
               NBMIN = 8
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NBMIN = 2
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NBMIN = 2
            END IF
         END IF
      END IF
      ILAENV = NBMIN
      RETURN
c
   70 CONTINUE
c
c     ISPEC = 3:  crossover point
c
      NX = 0
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. C3.EQ.
     $       'QLF' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NX = 32
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NX = 32
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NX = 128
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NX = 128
            END IF
         END IF
      END IF
      ILAENV = NX
      RETURN
c
   80 CONTINUE
c
c     ISPEC = 4:  number of shifts (used by xHSEQR)
c
      ILAENV = 6
      RETURN
c
   90 CONTINUE
c
c     ISPEC = 5:  minimum column dimension (not used)
c
      ILAENV = 2
      RETURN
c
  100 CONTINUE
c
c     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
c
      ILAENV = INT( REAL( MIN( N1, N2 ) )*1.6E0 )
      RETURN
c
  110 CONTINUE
c
c     ISPEC = 7:  number of processors (not used)
c
      ILAENV = 1
      RETURN
c
  120 CONTINUE
c
c     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
c
      ILAENV = 50
      RETURN
c
  130 CONTINUE
c
c     ISPEC = 9:  maximum size of the subproblems at the bottom of the
c                 computation tree in the divide-and-conquer algorithm
c                 (used by xGELSD and xGESDD)
c
      ILAENV = 25
      RETURN
c
  140 CONTINUE
c
c     ISPEC = 10: ieee NaN arithmetic can be trusted not to trap
c
c     ILAENV = 0
      ILAENV = 1
      IF( ILAENV.EQ.1 ) THEN
         ILAENV = IEEECK( 1, 0.0, 1.0 )
      END IF
      RETURN
c
  150 CONTINUE
c
c     ISPEC = 11: infinity arithmetic can be trusted not to trap
c
c     ILAENV = 0
      ILAENV = 1
      IF( ILAENV.EQ.1 ) THEN
         ILAENV = IEEECK( 0, 0.0, 1.0 )
      END IF
      RETURN
c
  160 CONTINUE
c
c     12 <= ISPEC <= 16: xHSEQR or one of its subroutines. 
c
      ILAENV = IPARMQ( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
      RETURN
c
c     End of ILAENV
c
      END

      INTEGER FUNCTION IPARMQ( ISPEC, NAME, OPTS, N, ILO, IHI, LWORK )
c
c  -- LAPACK auxiliary routine (version 3.4.0) --
c  -- LAPACK is a software package provided by Univ. of Tennessee,    --
c  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
c     November 2011
c
c     .. Scalar Arguments ..
      INTEGER            IHI, ILO, ISPEC, LWORK, N
      CHARACTER          NAME*( * ), OPTS*( * )
c
c  ================================================================
c     .. Parameters ..
      INTEGER            INMIN, INWIN, INIBL, ISHFTS, IACC22
      PARAMETER          ( INMIN = 12, INWIN = 13, INIBL = 14,
     $                   ISHFTS = 15, IACC22 = 16 )
      INTEGER            NMIN, K22MIN, KACMIN, NIBBLE, KNWSWP
      PARAMETER          ( NMIN = 75, K22MIN = 14, KACMIN = 14,
     $                   NIBBLE = 14, KNWSWP = 500 )
      REAL               TWO
      PARAMETER          ( TWO = 2.0 )
c     ..
c     .. Local Scalars ..
      INTEGER            NH, NS
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          LOG, MAX, MOD, NINT, REAL
c     ..
c     .. Executable Statements ..
      IF( ( ISPEC.EQ.ISHFTS ) .OR. ( ISPEC.EQ.INWIN ) .OR.
     $    ( ISPEC.EQ.IACC22 ) ) THEN
c
c        ==== Set the number simultaneous shifts ====
c
         NH = IHI - ILO + 1
         NS = 2
         IF( NH.GE.30 )
     $      NS = 4
         IF( NH.GE.60 )
     $      NS = 10
         IF( NH.GE.150 )
     $      NS = MAX( 10, NH / NINT( LOG( REAL( NH ) ) / LOG( TWO ) ) )
         IF( NH.GE.590 )
     $      NS = 64
         IF( NH.GE.3000 )
     $      NS = 128
         IF( NH.GE.6000 )
     $      NS = 256
         NS = MAX( 2, NS-MOD( NS, 2 ) )
      END IF
c
      IF( ISPEC.EQ.INMIN ) THEN
c
c
c        ===== Matrices of order smaller than NMIN get sent
c        .     to xLAHQR, the classic double shift algorithm.
c        .     This must be at least 11. ====
c
         IPARMQ = NMIN
c
      ELSE IF( ISPEC.EQ.INIBL ) THEN
c
c        ==== INIBL: skip a multi-shift qr iteration and
c        .    whenever aggressive early deflation finds
c        .    at least (NIBBLE*(window size)/100) deflations. ====
c
         IPARMQ = NIBBLE
c
      ELSE IF( ISPEC.EQ.ISHFTS ) THEN
c
c        ==== NSHFTS: The number of simultaneous shifts =====
c
         IPARMQ = NS
c
      ELSE IF( ISPEC.EQ.INWIN ) THEN
c
c        ==== NW: deflation window size.  ====
c
         IF( NH.LE.KNWSWP ) THEN
            IPARMQ = NS
         ELSE
            IPARMQ = 3*NS / 2
         END IF
c
      ELSE IF( ISPEC.EQ.IACC22 ) THEN
c
c        ==== IACC22: Whether to accumulate reflections
c        .     before updating the far-from-diagonal elements
c        .     and whether to use 2-by-2 block structure while
c        .     doing it.  A small amount of work could be saved
c        .     by making this choice dependent also upon the
c        .     NH=IHI-ILO+1.
c
         IPARMQ = 0
         IF( NS.GE.KACMIN )
     $      IPARMQ = 1
         IF( NS.GE.K22MIN )
     $      IPARMQ = 2
c
      ELSE
c        ===== invalid value of ispec =====
         IPARMQ = -1
c
      END IF
c
c     ==== End of IPARMQ ====
c
      END

      INTEGER          FUNCTION IEEECK( ISPEC, ZERO, ONE )
c
c  -- LAPACK auxiliary routine (version 3.4.0) --
c  -- LAPACK is a software package provided by Univ. of Tennessee,    --
c  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
c     November 2011
c
c     .. Scalar Arguments ..
      INTEGER            ISPEC
      REAL               ONE, ZERO
c     ..
c
c  =====================================================================
c
c     .. Local Scalars ..
      REAL               NAN1, NAN2, NAN3, NAN4, NAN5, NAN6, NEGINF,
     $                   NEGZRO, NEWZRO, POSINF
c     ..
c     .. Executable Statements ..
      IEEECK = 1
c
      POSINF = ONE / ZERO
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
c
      NEGINF = -ONE / ZERO
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
c
      NEGZRO = ONE / ( NEGINF+ONE )
      IF( NEGZRO.NE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
c
      NEGINF = ONE / NEGZRO
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
c
      NEWZRO = NEGZRO + ZERO
      IF( NEWZRO.NE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
c
      POSINF = ONE / NEWZRO
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
c
      NEGINF = NEGINF*POSINF
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
c
      POSINF = POSINF*POSINF
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
c
c
c
c
c     Return if we were only asked to check infinity arithmetic
c
      IF( ISPEC.EQ.0 )
     $   RETURN
c
      NAN1 = POSINF + NEGINF
c
      NAN2 = POSINF / NEGINF
c
      NAN3 = POSINF / POSINF
c
      NAN4 = POSINF*ZERO
c
      NAN5 = NEGINF*NEGZRO
c
      NAN6 = NAN5*ZERO
c
      IF( NAN1.EQ.NAN1 ) THEN
         IEEECK = 0
         RETURN
      END IF
c
      IF( NAN2.EQ.NAN2 ) THEN
         IEEECK = 0
         RETURN
      END IF
c
      IF( NAN3.EQ.NAN3 ) THEN
         IEEECK = 0
         RETURN
      END IF
c
      IF( NAN4.EQ.NAN4 ) THEN
         IEEECK = 0
         RETURN
      END IF
c
      IF( NAN5.EQ.NAN5 ) THEN
         IEEECK = 0
         RETURN
      END IF
c
      IF( NAN6.EQ.NAN6 ) THEN
         IEEECK = 0
         RETURN
      END IF
c
      RETURN
      END

      SUBROUTINE DLASR( SIDE, PIVOT, DIRECT, M, N, C, S, A, LDA )
c
c  -- LAPACK auxiliary routine (version 3.4.2) --
c  -- LAPACK is a software package provided by Univ. of Tennessee,    --
c  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
c     September 2012
c
c     .. Scalar Arguments ..
      CHARACTER          DIRECT, PIVOT, SIDE
      INTEGER            LDA, M, N
c     ..
c     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( * ), S( * )
c     ..
c
c  =====================================================================
c
c     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
c     ..
c     .. Local Scalars ..
      INTEGER            I, INFO, J
      DOUBLE PRECISION   CTEMP, STEMP, TEMP
c     ..
c     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
c     ..
c     .. External Subroutines ..
      EXTERNAL           XERBLA
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          MAX
c     ..
c     .. Executable Statements ..
c
c     Test the input parameters
c
      INFO = 0
      IF( .NOT.( LSAME( SIDE, 'L' ) .OR. LSAME( SIDE, 'R' ) ) ) THEN
         INFO = 1
      ELSE IF( .NOT.( LSAME( PIVOT, 'V' ) .OR. LSAME( PIVOT,
     $         'T' ) .OR. LSAME( PIVOT, 'B' ) ) ) THEN
         INFO = 2
      ELSE IF( .NOT.( LSAME( DIRECT, 'F' ) .OR. LSAME( DIRECT, 'B' ) ) )
     $          THEN
         INFO = 3
      ELSE IF( M.LT.0 ) THEN
         INFO = 4
      ELSE IF( N.LT.0 ) THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLASR ', INFO )
         RETURN
      END IF
c
c     Quick return if possible
c
      IF( ( M.EQ.0 ) .OR. ( N.EQ.0 ) )
     $   RETURN
      IF( LSAME( SIDE, 'L' ) ) THEN
c
c        Form  P * A
c
         IF( LSAME( PIVOT, 'V' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 20 J = 1, M - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 10 I = 1, N
                        TEMP = A( J+1, I )
                        A( J+1, I ) = CTEMP*TEMP - STEMP*A( J, I )
                        A( J, I ) = STEMP*TEMP + CTEMP*A( J, I )
   10                CONTINUE
                  END IF
   20          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 40 J = M - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 30 I = 1, N
                        TEMP = A( J+1, I )
                        A( J+1, I ) = CTEMP*TEMP - STEMP*A( J, I )
                        A( J, I ) = STEMP*TEMP + CTEMP*A( J, I )
   30                CONTINUE
                  END IF
   40          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'T' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 60 J = 2, M
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 50 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = CTEMP*TEMP - STEMP*A( 1, I )
                        A( 1, I ) = STEMP*TEMP + CTEMP*A( 1, I )
   50                CONTINUE
                  END IF
   60          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 80 J = M, 2, -1
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 70 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = CTEMP*TEMP - STEMP*A( 1, I )
                        A( 1, I ) = STEMP*TEMP + CTEMP*A( 1, I )
   70                CONTINUE
                  END IF
   80          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'B' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 100 J = 1, M - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 90 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = STEMP*A( M, I ) + CTEMP*TEMP
                        A( M, I ) = CTEMP*A( M, I ) - STEMP*TEMP
   90                CONTINUE
                  END IF
  100          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 120 J = M - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 110 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = STEMP*A( M, I ) + CTEMP*TEMP
                        A( M, I ) = CTEMP*A( M, I ) - STEMP*TEMP
  110                CONTINUE
                  END IF
  120          CONTINUE
            END IF
         END IF
      ELSE IF( LSAME( SIDE, 'R' ) ) THEN
c
c        Form A * P**T
c
         IF( LSAME( PIVOT, 'V' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 140 J = 1, N - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 130 I = 1, M
                        TEMP = A( I, J+1 )
                        A( I, J+1 ) = CTEMP*TEMP - STEMP*A( I, J )
                        A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
  130                CONTINUE
                  END IF
  140          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 160 J = N - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 150 I = 1, M
                        TEMP = A( I, J+1 )
                        A( I, J+1 ) = CTEMP*TEMP - STEMP*A( I, J )
                        A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
  150                CONTINUE
                  END IF
  160          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'T' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 180 J = 2, N
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 170 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP - STEMP*A( I, 1 )
                        A( I, 1 ) = STEMP*TEMP + CTEMP*A( I, 1 )
  170                CONTINUE
                  END IF
  180          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 200 J = N, 2, -1
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 190 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP - STEMP*A( I, 1 )
                        A( I, 1 ) = STEMP*TEMP + CTEMP*A( I, 1 )
  190                CONTINUE
                  END IF
  200          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'B' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 220 J = 1, N - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 210 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = STEMP*A( I, N ) + CTEMP*TEMP
                        A( I, N ) = CTEMP*A( I, N ) - STEMP*TEMP
  210                CONTINUE
                  END IF
  220          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 240 J = N - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 230 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = STEMP*A( I, N ) + CTEMP*TEMP
                        A( I, N ) = CTEMP*A( I, N ) - STEMP*TEMP
  230                CONTINUE
                  END IF
  240          CONTINUE
            END IF
         END IF
      END IF
c
      RETURN
c
c     End of DLASR
c
      END

      SUBROUTINE DORGTR( UPLO, N, A, LDA, TAU, WORK, LWORK, INFO )
c
c  -- LAPACK computational routine (version 3.4.0) --
c  -- LAPACK is a software package provided by Univ. of Tennessee,    --
c  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
c     November 2011
c
c     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, LWORK, N
c     ..
c     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
c     ..
c
c  =====================================================================
c
c     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
c     ..
c     .. Local Scalars ..
      LOGICAL            LQUERY, UPPER
      INTEGER            I, IINFO, J, LWKOPT, NB
c     ..
c     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
c     ..
c     .. External Subroutines ..
      EXTERNAL           DORGQL, DORGQR, XERBLA
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          MAX
c     ..
c     .. Executable Statements ..
c
c     Test the input arguments
c
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LWORK.LT.MAX( 1, N-1 ) .AND. .NOT.LQUERY ) THEN
         INFO = -7
      END IF
c
      IF( INFO.EQ.0 ) THEN
         IF( UPPER ) THEN
            NB = ILAENV( 1, 'DORGQL', ' ', N-1, N-1, N-1, -1 )
         ELSE
            NB = ILAENV( 1, 'DORGQR', ' ', N-1, N-1, N-1, -1 )
         END IF
         LWKOPT = MAX( 1, N-1 )*NB
         WORK( 1 ) = LWKOPT
      END IF
c
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORGTR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
c
c     Quick return if possible
c
      IF( N.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
c
      IF( UPPER ) THEN
c
c        Q was determined by a call to DSYTRD with UPLO = 'U'
c
c        Shift the vectors which define the elementary reflectors one
c        column to the left, and set the last row and column of Q to
c        those of the unit matrix
c
         DO 20 J = 1, N - 1
            DO 10 I = 1, J - 1
               A( I, J ) = A( I, J+1 )
   10       CONTINUE
            A( N, J ) = ZERO
   20    CONTINUE
         DO 30 I = 1, N - 1
            A( I, N ) = ZERO
   30    CONTINUE
         A( N, N ) = ONE
c
c        Generate Q(1:n-1,1:n-1)
c
         CALL DORGQL( N-1, N-1, N-1, A, LDA, TAU, WORK, LWORK, IINFO )
c
      ELSE
c
c        Q was determined by a call to DSYTRD with UPLO = 'L'.
c
c        Shift the vectors which define the elementary reflectors one
c        column to the right, and set the first row and column of Q to
c        those of the unit matrix
c
         DO 50 J = N, 2, -1
            A( 1, J ) = ZERO
            DO 40 I = J + 1, N
               A( I, J ) = A( I, J-1 )
   40       CONTINUE
   50    CONTINUE
         A( 1, 1 ) = ONE
         DO 60 I = 2, N
            A( I, 1 ) = ZERO
   60    CONTINUE
         IF( N.GT.1 ) THEN
c
c           Generate Q(2:n,2:n)
c
            CALL DORGQR( N-1, N-1, N-1, A( 2, 2 ), LDA, TAU, WORK,
     $                   LWORK, IINFO )
         END IF
      END IF
      WORK( 1 ) = LWKOPT
      RETURN
c
c     End of DORGTR
c
      END

      SUBROUTINE DSTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO )
c
c  -- LAPACK computational routine (version 3.4.0) --
c  -- LAPACK is a software package provided by Univ. of Tennessee,    --
c  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
c     November 2011
c
c     .. Scalar Arguments ..
      CHARACTER          COMPZ
      INTEGER            INFO, LDZ, N
c     ..
c     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), E( * ), WORK( * ), Z( LDZ, * )
c     ..
c
c  =====================================================================
c
c     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO, THREE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0,
     $                   THREE = 3.0D0 )
      INTEGER            MAXIT
      PARAMETER          ( MAXIT = 30 )
c     ..
c     .. Local Scalars ..
      INTEGER            I, ICOMPZ, II, ISCALE, J, JTOT, K, L, L1, LEND,
     $                   LENDM1, LENDP1, LENDSV, LM1, LSV, M, MM, MM1,
     $                   NM1, NMAXIT
      DOUBLE PRECISION   ANORM, B, C, EPS, EPS2, F, G, P, R, RT1, RT2,
     $                   S, SAFMAX, SAFMIN, SSFMAX, SSFMIN, TST
c     ..
c     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, DLANST, DLAPY2
      EXTERNAL           LSAME, DLAMCH, DLANST, DLAPY2
c     ..
c     .. External Subroutines ..
      EXTERNAL           DLAE2, DLAEV2, DLARTG, DLASCL, DLASET, DLASR,
     $                   DLASRT, DSWAP, XERBLA
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SIGN, SQRT
c     ..
c     .. Executable Statements ..
c
c     Test the input parameters.
c
      INFO = 0
c
      IF( LSAME( COMPZ, 'N' ) ) THEN
         ICOMPZ = 0
      ELSE IF( LSAME( COMPZ, 'V' ) ) THEN
         ICOMPZ = 1
      ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
         ICOMPZ = 2
      ELSE
         ICOMPZ = -1
      END IF
      IF( ICOMPZ.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( ( LDZ.LT.1 ) .OR. ( ICOMPZ.GT.0 .AND. LDZ.LT.MAX( 1,
     $         N ) ) ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSTEQR', -INFO )
         RETURN
      END IF
c
c     Quick return if possible
c
      IF( N.EQ.0 )
     $   RETURN
c
      IF( N.EQ.1 ) THEN
         IF( ICOMPZ.EQ.2 )
     $      Z( 1, 1 ) = ONE
         RETURN
      END IF
c
c     Determine the unit roundoff and over/underflow thresholds.
c
      EPS = DLAMCH( 'E' )
      EPS2 = EPS**2
      SAFMIN = DLAMCH( 'S' )
      SAFMAX = ONE / SAFMIN
      SSFMAX = SQRT( SAFMAX ) / THREE
      SSFMIN = SQRT( SAFMIN ) / EPS2
c
c     Compute the eigenvalues and eigenvectors of the tridiagonal
c     matrix.
c
      IF( ICOMPZ.EQ.2 )
     $   CALL DLASET( 'Full', N, N, ZERO, ONE, Z, LDZ )
c
      NMAXIT = N*MAXIT
      JTOT = 0
c
c     Determine where the matrix splits and choose QL or QR iteration
c     for each block, according to whether top or bottom diagonal
c     element is smaller.
c
      L1 = 1
      NM1 = N - 1
c
   10 CONTINUE
      IF( L1.GT.N )
     $   GO TO 160
      IF( L1.GT.1 )
     $   E( L1-1 ) = ZERO
      IF( L1.LE.NM1 ) THEN
         DO 20 M = L1, NM1
            TST = ABS( E( M ) )
            IF( TST.EQ.ZERO )
     $         GO TO 30
            IF( TST.LE.( SQRT( ABS( D( M ) ) )*SQRT( ABS( D( M+
     $          1 ) ) ) )*EPS ) THEN
               E( M ) = ZERO
               GO TO 30
            END IF
   20    CONTINUE
      END IF
      M = N
c
   30 CONTINUE
      L = L1
      LSV = L
      LEND = M
      LENDSV = LEND
      L1 = M + 1
      IF( LEND.EQ.L )
     $   GO TO 10
c
c     Scale submatrix in rows and columns L to LEND
c
      ANORM = DLANST( 'M', LEND-L+1, D( L ), E( L ) )
      ISCALE = 0
      IF( ANORM.EQ.ZERO )
     $   GO TO 10
      IF( ANORM.GT.SSFMAX ) THEN
         ISCALE = 1
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L+1, 1, D( L ), N,
     $                INFO )
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L, 1, E( L ), N,
     $                INFO )
      ELSE IF( ANORM.LT.SSFMIN ) THEN
         ISCALE = 2
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L+1, 1, D( L ), N,
     $                INFO )
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L, 1, E( L ), N,
     $                INFO )
      END IF
c
c     Choose between QL and QR iteration
c
      IF( ABS( D( LEND ) ).LT.ABS( D( L ) ) ) THEN
         LEND = LSV
         L = LENDSV
      END IF
c
      IF( LEND.GT.L ) THEN
c
c        QL Iteration
c
c        Look for small subdiagonal element.
c
   40    CONTINUE
         IF( L.NE.LEND ) THEN
            LENDM1 = LEND - 1
            DO 50 M = L, LENDM1
               TST = ABS( E( M ) )**2
               IF( TST.LE.( EPS2*ABS( D( M ) ) )*ABS( D( M+1 ) )+
     $             SAFMIN )GO TO 60
   50       CONTINUE
         END IF
c
         M = LEND
c
   60    CONTINUE
         IF( M.LT.LEND )
     $      E( M ) = ZERO
         P = D( L )
         IF( M.EQ.L )
     $      GO TO 80
c
c        If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
c        to compute its eigensystem.
c
         IF( M.EQ.L+1 ) THEN
            IF( ICOMPZ.GT.0 ) THEN
               CALL DLAEV2( D( L ), E( L ), D( L+1 ), RT1, RT2, C, S )
               WORK( L ) = C
               WORK( N-1+L ) = S
               CALL DLASR( 'R', 'V', 'B', N, 2, WORK( L ),
     $                     WORK( N-1+L ), Z( 1, L ), LDZ )
            ELSE
               CALL DLAE2( D( L ), E( L ), D( L+1 ), RT1, RT2 )
            END IF
            D( L ) = RT1
            D( L+1 ) = RT2
            E( L ) = ZERO
            L = L + 2
            IF( L.LE.LEND )
     $         GO TO 40
            GO TO 140
         END IF
c
         IF( JTOT.EQ.NMAXIT )
     $      GO TO 140
         JTOT = JTOT + 1
c
c        Form shift.
c
         G = ( D( L+1 )-P ) / ( TWO*E( L ) )
         R = DLAPY2( G, ONE )
         G = D( M ) - P + ( E( L ) / ( G+SIGN( R, G ) ) )
c
         S = ONE
         C = ONE
         P = ZERO
c
c        Inner loop
c
         MM1 = M - 1
         DO 70 I = MM1, L, -1
            F = S*E( I )
            B = C*E( I )
            CALL DLARTG( G, F, C, S, R )
            IF( I.NE.M-1 )
     $         E( I+1 ) = R
            G = D( I+1 ) - P
            R = ( D( I )-G )*S + TWO*C*B
            P = S*R
            D( I+1 ) = G + P
            G = C*R - B
c
c           If eigenvectors are desired, then save rotations.
c
            IF( ICOMPZ.GT.0 ) THEN
               WORK( I ) = C
               WORK( N-1+I ) = -S
            END IF
c
   70    CONTINUE
c
c        If eigenvectors are desired, then apply saved rotations.
c
         IF( ICOMPZ.GT.0 ) THEN
            MM = M - L + 1
            CALL DLASR( 'R', 'V', 'B', N, MM, WORK( L ), WORK( N-1+L ),
     $                  Z( 1, L ), LDZ )
         END IF
c
         D( L ) = D( L ) - P
         E( L ) = G
         GO TO 40
c
c        Eigenvalue found.
c
   80    CONTINUE
         D( L ) = P
c
         L = L + 1
         IF( L.LE.LEND )
     $      GO TO 40
         GO TO 140
c
      ELSE
c
c        QR Iteration
c
c        Look for small superdiagonal element.
c
   90    CONTINUE
         IF( L.NE.LEND ) THEN
            LENDP1 = LEND + 1
            DO 100 M = L, LENDP1, -1
               TST = ABS( E( M-1 ) )**2
               IF( TST.LE.( EPS2*ABS( D( M ) ) )*ABS( D( M-1 ) )+
     $             SAFMIN )GO TO 110
  100       CONTINUE
         END IF
c
         M = LEND
c
  110    CONTINUE
         IF( M.GT.LEND )
     $      E( M-1 ) = ZERO
         P = D( L )
         IF( M.EQ.L )
     $      GO TO 130
c
c        If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
c        to compute its eigensystem.
c
         IF( M.EQ.L-1 ) THEN
            IF( ICOMPZ.GT.0 ) THEN
               CALL DLAEV2( D( L-1 ), E( L-1 ), D( L ), RT1, RT2, C, S )
               WORK( M ) = C
               WORK( N-1+M ) = S
               CALL DLASR( 'R', 'V', 'F', N, 2, WORK( M ),
     $                     WORK( N-1+M ), Z( 1, L-1 ), LDZ )
            ELSE
               CALL DLAE2( D( L-1 ), E( L-1 ), D( L ), RT1, RT2 )
            END IF
            D( L-1 ) = RT1
            D( L ) = RT2
            E( L-1 ) = ZERO
            L = L - 2
            IF( L.GE.LEND )
     $         GO TO 90
            GO TO 140
         END IF
c
         IF( JTOT.EQ.NMAXIT )
     $      GO TO 140
         JTOT = JTOT + 1
c
c        Form shift.
c
         G = ( D( L-1 )-P ) / ( TWO*E( L-1 ) )
         R = DLAPY2( G, ONE )
         G = D( M ) - P + ( E( L-1 ) / ( G+SIGN( R, G ) ) )
c
         S = ONE
         C = ONE
         P = ZERO
c
c        Inner loop
c
         LM1 = L - 1
         DO 120 I = M, LM1
            F = S*E( I )
            B = C*E( I )
            CALL DLARTG( G, F, C, S, R )
            IF( I.NE.M )
     $         E( I-1 ) = R
            G = D( I ) - P
            R = ( D( I+1 )-G )*S + TWO*C*B
            P = S*R
            D( I ) = G + P
            G = C*R - B
c
c           If eigenvectors are desired, then save rotations.
c
            IF( ICOMPZ.GT.0 ) THEN
               WORK( I ) = C
               WORK( N-1+I ) = S
            END IF
c
  120    CONTINUE
c
c        If eigenvectors are desired, then apply saved rotations.
c
         IF( ICOMPZ.GT.0 ) THEN
            MM = L - M + 1
            CALL DLASR( 'R', 'V', 'F', N, MM, WORK( M ), WORK( N-1+M ),
     $                  Z( 1, M ), LDZ )
         END IF
c
         D( L ) = D( L ) - P
         E( LM1 ) = G
         GO TO 90
c
c        Eigenvalue found.
c
  130    CONTINUE
         D( L ) = P
c
         L = L - 1
         IF( L.GE.LEND )
     $      GO TO 90
         GO TO 140
c
      END IF
c
c     Undo scaling if necessary
c
  140 CONTINUE
      IF( ISCALE.EQ.1 ) THEN
         CALL DLASCL( 'G', 0, 0, SSFMAX, ANORM, LENDSV-LSV+1, 1,
     $                D( LSV ), N, INFO )
         CALL DLASCL( 'G', 0, 0, SSFMAX, ANORM, LENDSV-LSV, 1, E( LSV ),
     $                N, INFO )
      ELSE IF( ISCALE.EQ.2 ) THEN
         CALL DLASCL( 'G', 0, 0, SSFMIN, ANORM, LENDSV-LSV+1, 1,
     $                D( LSV ), N, INFO )
         CALL DLASCL( 'G', 0, 0, SSFMIN, ANORM, LENDSV-LSV, 1, E( LSV ),
     $                N, INFO )
      END IF
c
c     Check for no convergence to an eigenvalue after a total
c     of N*MAXIT iterations.
c
      IF( JTOT.LT.NMAXIT )
     $   GO TO 10
      DO 150 I = 1, N - 1
         IF( E( I ).NE.ZERO )
     $      INFO = INFO + 1
  150 CONTINUE
      GO TO 190
c
c     Order eigenvalues and eigenvectors.
c
  160 CONTINUE
      IF( ICOMPZ.EQ.0 ) THEN
c
c        Use Quick Sort
c
         CALL DLASRT( 'I', N, D, INFO )
c
      ELSE
c
c        Use Selection Sort to minimize swaps of eigenvectors
c
         DO 180 II = 2, N
            I = II - 1
            K = I
            P = D( I )
            DO 170 J = II, N
               IF( D( J ).LT.P ) THEN
                  K = J
                  P = D( J )
               END IF
  170       CONTINUE
            IF( K.NE.I ) THEN
               D( K ) = D( I )
               D( I ) = P
               CALL DSWAP( N, Z( 1, I ), 1, Z( 1, K ), 1 )
            END IF
  180    CONTINUE
      END IF
c
  190 CONTINUE
      RETURN
c
c     End of DSTEQR
c
      END

      SUBROUTINE DORGQL( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
c
c  -- LAPACK computational routine (version 3.4.0) --
c  -- LAPACK is a software package provided by Univ. of Tennessee,    --
c  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
c     November 2011
c
c     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, LWORK, M, N
c     ..
c     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
c     ..
c
c  =====================================================================
c
c     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
c     ..
c     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWS, J, KK, L, LDWORK, LWKOPT,
     $                   NB, NBMIN, NX
c     ..
c     .. External Subroutines ..
      EXTERNAL           DLARFB, DLARFT, DORG2L, XERBLA
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
c     ..
c     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
c     ..
c     .. Executable Statements ..
c
c     Test the input arguments
c
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      END IF
c
      IF( INFO.EQ.0 ) THEN
         IF( N.EQ.0 ) THEN
            LWKOPT = 1
         ELSE
            NB = ILAENV( 1, 'DORGQL', ' ', M, N, K, -1 )
            LWKOPT = N*NB
         END IF
         WORK( 1 ) = LWKOPT
c
         IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
            INFO = -8
         END IF
      END IF
c
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORGQL', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
c
c     Quick return if possible
c
      IF( N.LE.0 ) THEN
         RETURN
      END IF
c
      NBMIN = 2
      NX = 0
      IWS = N
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
c
c        Determine when to cross over from blocked to unblocked code.
c
         NX = MAX( 0, ILAENV( 3, 'DORGQL', ' ', M, N, K, -1 ) )
         IF( NX.LT.K ) THEN
c
c           Determine if workspace is large enough for blocked code.
c
            LDWORK = N
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
c
c              Not enough workspace to use optimal NB:  reduce NB and
c              determine the minimum value of NB.
c
               NB = LWORK / LDWORK
               NBMIN = MAX( 2, ILAENV( 2, 'DORGQL', ' ', M, N, K, -1 ) )
            END IF
         END IF
      END IF
c
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
c
c        Use blocked code after the first block.
c        The last kk columns are handled by the block method.
c
         KK = MIN( K, ( ( K-NX+NB-1 ) / NB )*NB )
c
c        Set A(m-kk+1:m,1:n-kk) to zero.
c
         DO 20 J = 1, N - KK
            DO 10 I = M - KK + 1, M
               A( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
      ELSE
         KK = 0
      END IF
c
c     Use unblocked code for the first or only block.
c
      CALL DORG2L( M-KK, N-KK, K-KK, A, LDA, TAU, WORK, IINFO )
c
      IF( KK.GT.0 ) THEN
c
c        Use blocked code
c
         DO 50 I = K - KK + 1, K, NB
            IB = MIN( NB, K-I+1 )
            IF( N-K+I.GT.1 ) THEN
c
c              Form the triangular factor of the block reflector
c              H = H(i+ib-1) . . . H(i+1) H(i)
c
               CALL DLARFT( 'Backward', 'Columnwise', M-K+I+IB-1, IB,
     $                      A( 1, N-K+I ), LDA, TAU( I ), WORK, LDWORK )
c
c              Apply H to A(1:m-k+i+ib-1,1:n-k+i-1) from the left
c
               CALL DLARFB( 'Left', 'No transpose', 'Backward',
     $                      'Columnwise', M-K+I+IB-1, N-K+I-1, IB,
     $                      A( 1, N-K+I ), LDA, WORK, LDWORK, A, LDA,
     $                      WORK( IB+1 ), LDWORK )
            END IF
c
c           Apply H to rows 1:m-k+i+ib-1 of current block
c
            CALL DORG2L( M-K+I+IB-1, IB, IB, A( 1, N-K+I ), LDA,
     $                   TAU( I ), WORK, IINFO )
c
c           Set rows m-k+i+ib:m of current block to zero
c
            DO 40 J = N - K + I, N - K + I + IB - 1
               DO 30 L = M - K + I + IB, M
                  A( L, J ) = ZERO
   30          CONTINUE
   40       CONTINUE
   50    CONTINUE
      END IF
c
      WORK( 1 ) = IWS
      RETURN
c
c     End of DORGQL
c
      END

      SUBROUTINE DORGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
c
c  -- LAPACK computational routine (version 3.4.0) --
c  -- LAPACK is a software package provided by Univ. of Tennessee,    --
c  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
c     November 2011
c
c     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, LWORK, M, N
c     ..
c     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
c     ..
c
c  =====================================================================
c
c     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
c     ..
c     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWS, J, KI, KK, L, LDWORK,
     $                   LWKOPT, NB, NBMIN, NX
c     ..
c     .. External Subroutines ..
      EXTERNAL           DLARFB, DLARFT, DORG2R, XERBLA
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
c     ..
c     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
c     ..
c     .. Executable Statements ..
c
c     Test the input arguments
c
      INFO = 0
      NB = ILAENV( 1, 'DORGQR', ' ', M, N, K, -1 )
      LWKOPT = MAX( 1, N )*NB
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORGQR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
c
c     Quick return if possible
c
      IF( N.LE.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
c
      NBMIN = 2
      NX = 0
      IWS = N
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
c
c        Determine when to cross over from blocked to unblocked code.
c
         NX = MAX( 0, ILAENV( 3, 'DORGQR', ' ', M, N, K, -1 ) )
         IF( NX.LT.K ) THEN
c
c           Determine if workspace is large enough for blocked code.
c
            LDWORK = N
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
c
c              Not enough workspace to use optimal NB:  reduce NB and
c              determine the minimum value of NB.
c
               NB = LWORK / LDWORK
               NBMIN = MAX( 2, ILAENV( 2, 'DORGQR', ' ', M, N, K, -1 ) )
            END IF
         END IF
      END IF
c
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
c
c        Use blocked code after the last block.
c        The first kk columns are handled by the block method.
c
         KI = ( ( K-NX-1 ) / NB )*NB
         KK = MIN( K, KI+NB )
c
c        Set A(1:kk,kk+1:n) to zero.
c
         DO 20 J = KK + 1, N
            DO 10 I = 1, KK
               A( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
      ELSE
         KK = 0
      END IF
c
c     Use unblocked code for the last or only block.
c
      IF( KK.LT.N )
     $   CALL DORG2R( M-KK, N-KK, K-KK, A( KK+1, KK+1 ), LDA,
     $                TAU( KK+1 ), WORK, IINFO )
c
      IF( KK.GT.0 ) THEN
c
c        Use blocked code
c
         DO 50 I = KI + 1, 1, -NB
            IB = MIN( NB, K-I+1 )
            IF( I+IB.LE.N ) THEN
c
c              Form the triangular factor of the block reflector
c              H = H(i) H(i+1) . . . H(i+ib-1)
c
               CALL DLARFT( 'Forward', 'Columnwise', M-I+1, IB,
     $                      A( I, I ), LDA, TAU( I ), WORK, LDWORK )
c
c              Apply H to A(i:m,i+ib:n) from the left
c
               CALL DLARFB( 'Left', 'No transpose', 'Forward',
     $                      'Columnwise', M-I+1, N-I-IB+1, IB,
     $                      A( I, I ), LDA, WORK, LDWORK, A( I, I+IB ),
     $                      LDA, WORK( IB+1 ), LDWORK )
            END IF
c
c           Apply H to rows i:m of current block
c
            CALL DORG2R( M-I+1, IB, IB, A( I, I ), LDA, TAU( I ), WORK,
     $                   IINFO )
c
c           Set rows 1:i-1 of current block to zero
c
            DO 40 J = I, I + IB - 1
               DO 30 L = 1, I - 1
                  A( L, J ) = ZERO
   30          CONTINUE
   40       CONTINUE
   50    CONTINUE
      END IF
c
      WORK( 1 ) = IWS
      RETURN
c
c     End of DORGQR
c
      END

      SUBROUTINE DORG2L( M, N, K, A, LDA, TAU, WORK, INFO )
c
c  -- LAPACK computational routine (version 3.4.2) --
c  -- LAPACK is a software package provided by Univ. of Tennessee,    --
c  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
c     September 2012
c
c     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, M, N
c     ..
c     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
c     ..
c
c  =====================================================================
c
c     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
c     ..
c     .. Local Scalars ..
      INTEGER            I, II, J, L
c     ..
c     .. External Subroutines ..
      EXTERNAL           DLARF, DSCAL, XERBLA
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          MAX
c     ..
c     .. Executable Statements ..
c
c     Test the input arguments
c
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORG2L', -INFO )
         RETURN
      END IF
c
c     Quick return if possible
c
      IF( N.LE.0 )
     $   RETURN
c
c     Initialise columns 1:n-k to columns of the unit matrix
c
      DO 20 J = 1, N - K
         DO 10 L = 1, M
            A( L, J ) = ZERO
   10    CONTINUE
         A( M-N+J, J ) = ONE
   20 CONTINUE
c
      DO 40 I = 1, K
         II = N - K + I
c
c        Apply H(i) to A(1:m-k+i,1:n-k+i) from the left
c
         A( M-N+II, II ) = ONE
         CALL DLARF( 'Left', M-N+II, II-1, A( 1, II ), 1, TAU( I ), A,
     $               LDA, WORK )
         CALL DSCAL( M-N+II-1, -TAU( I ), A( 1, II ), 1 )
         A( M-N+II, II ) = ONE - TAU( I )
c
c        Set A(m-k+i+1:m,n-k+i) to zero
c
         DO 30 L = M - N + II + 1, M
            A( L, II ) = ZERO
   30    CONTINUE
   40 CONTINUE
      RETURN
c
c     End of DORG2L
c
      END

      DOUBLE PRECISION FUNCTION DLAPY2( X, Y )
c
c  -- LAPACK auxiliary routine (version 3.4.2) --
c  -- LAPACK is a software package provided by Univ. of Tennessee,    --
c  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
c     September 2012
c
c     .. Scalar Arguments ..
      DOUBLE PRECISION   X, Y
c     ..
c
c  =====================================================================
c
c     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
c     ..
c     .. Local Scalars ..
      DOUBLE PRECISION   W, XABS, YABS, Z
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
c     ..
c     .. Executable Statements ..
c
      XABS = ABS( X )
      YABS = ABS( Y )
      W = MAX( XABS, YABS )
      Z = MIN( XABS, YABS )
      IF( Z.EQ.ZERO ) THEN
         DLAPY2 = W
      ELSE
         DLAPY2 = W*SQRT( ONE+( Z / W )**2 )
      END IF
      RETURN
c
c     End of DLAPY2
c
      END

      SUBROUTINE DLAE2( A, B, C, RT1, RT2 )
c
c  -- LAPACK auxiliary routine (version 3.4.2) --
c  -- LAPACK is a software package provided by Univ. of Tennessee,    --
c  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
c     September 2012
c
c     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B, C, RT1, RT2
c     ..
c
c =====================================================================
c
c     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D0 )
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   HALF
      PARAMETER          ( HALF = 0.5D0 )
c     ..
c     .. Local Scalars ..
      DOUBLE PRECISION   AB, ACMN, ACMX, ADF, DF, RT, SM, TB
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          ABS, SQRT
c     ..
c     .. Executable Statements ..
c
c     Compute the eigenvalues
c
      SM = A + C
      DF = A - C
      ADF = ABS( DF )
      TB = B + B
      AB = ABS( TB )
      IF( ABS( A ).GT.ABS( C ) ) THEN
         ACMX = A
         ACMN = C
      ELSE
         ACMX = C
         ACMN = A
      END IF
      IF( ADF.GT.AB ) THEN
         RT = ADF*SQRT( ONE+( AB / ADF )**2 )
      ELSE IF( ADF.LT.AB ) THEN
         RT = AB*SQRT( ONE+( ADF / AB )**2 )
      ELSE
c
c        Includes case AB=ADF=0
c
         RT = AB*SQRT( TWO )
      END IF
      IF( SM.LT.ZERO ) THEN
         RT1 = HALF*( SM-RT )
c
c        Order of execution important.
c        To get fully accurate smaller eigenvalue,
c        next line needs to be executed in higher precision.
c
         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      ELSE IF( SM.GT.ZERO ) THEN
         RT1 = HALF*( SM+RT )
c
c        Order of execution important.
c        To get fully accurate smaller eigenvalue,
c        next line needs to be executed in higher precision.
c
         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      ELSE
c
c        Includes case RT1 = RT2 = 0
c
         RT1 = HALF*RT
         RT2 = -HALF*RT
      END IF
      RETURN
c
c     End of DLAE2
c
      END

      SUBROUTINE DLAEV2( A, B, C, RT1, RT2, CS1, SN1 )
c
c  -- LAPACK auxiliary routine (version 3.4.2) --
c  -- LAPACK is a software package provided by Univ. of Tennessee,    --
c  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
c     September 2012
c
c     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B, C, CS1, RT1, RT2, SN1
c     ..
c
c =====================================================================
c
c     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D0 )
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   HALF
      PARAMETER          ( HALF = 0.5D0 )
c     ..
c     .. Local Scalars ..
      INTEGER            SGN1, SGN2
      DOUBLE PRECISION   AB, ACMN, ACMX, ACS, ADF, CS, CT, DF, RT, SM,
     $                   TB, TN
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          ABS, SQRT
c     ..
c     .. Executable Statements ..
c
c     Compute the eigenvalues
c
      SM = A + C
      DF = A - C
      ADF = ABS( DF )
      TB = B + B
      AB = ABS( TB )
      IF( ABS( A ).GT.ABS( C ) ) THEN
         ACMX = A
         ACMN = C
      ELSE
         ACMX = C
         ACMN = A
      END IF
      IF( ADF.GT.AB ) THEN
         RT = ADF*SQRT( ONE+( AB / ADF )**2 )
      ELSE IF( ADF.LT.AB ) THEN
         RT = AB*SQRT( ONE+( ADF / AB )**2 )
      ELSE
c
c        Includes case AB=ADF=0
c
         RT = AB*SQRT( TWO )
      END IF
      IF( SM.LT.ZERO ) THEN
         RT1 = HALF*( SM-RT )
         SGN1 = -1
c
c        Order of execution important.
c        To get fully accurate smaller eigenvalue,
c        next line needs to be executed in higher precision.
c
         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      ELSE IF( SM.GT.ZERO ) THEN
         RT1 = HALF*( SM+RT )
         SGN1 = 1
c
c        Order of execution important.
c        To get fully accurate smaller eigenvalue,
c        next line needs to be executed in higher precision.
c
         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      ELSE
c
c        Includes case RT1 = RT2 = 0
c
         RT1 = HALF*RT
         RT2 = -HALF*RT
         SGN1 = 1
      END IF
c
c     Compute the eigenvector
c
      IF( DF.GE.ZERO ) THEN
         CS = DF + RT
         SGN2 = 1
      ELSE
         CS = DF - RT
         SGN2 = -1
      END IF
      ACS = ABS( CS )
      IF( ACS.GT.AB ) THEN
         CT = -TB / CS
         SN1 = ONE / SQRT( ONE+CT*CT )
         CS1 = CT*SN1
      ELSE
         IF( AB.EQ.ZERO ) THEN
            CS1 = ONE
            SN1 = ZERO
         ELSE
            TN = -CS / TB
            CS1 = ONE / SQRT( ONE+TN*TN )
            SN1 = TN*CS1
         END IF
      END IF
      IF( SGN1.EQ.SGN2 ) THEN
         TN = CS1
         CS1 = -SN1
         SN1 = TN
      END IF
      RETURN
c
c     End of DLAEV2
c
      END

      SUBROUTINE DSTERF( N, D, E, INFO )
c
c  -- LAPACK computational routine (version 3.4.0) --
c  -- LAPACK is a software package provided by Univ. of Tennessee,    --
c  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
c     November 2011
c
c     .. Scalar Arguments ..
      INTEGER            INFO, N
c     ..
c     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), E( * )
c     ..
c
c  =====================================================================
c
c     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO, THREE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0,
     $                   THREE = 3.0D0 )
      INTEGER            MAXIT
      PARAMETER          ( MAXIT = 30 )
c     ..
c     .. Local Scalars ..
      INTEGER            I, ISCALE, JTOT, L, L1, LEND, LENDSV, LSV, M,
     $                   NMAXIT
      DOUBLE PRECISION   ALPHA, ANORM, BB, C, EPS, EPS2, GAMMA, OLDC,
     $                   OLDGAM, P, R, RT1, RT2, RTE, S, SAFMAX, SAFMIN,
     $                   SIGMA, SSFMAX, SSFMIN, RMAX
c     ..
c     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLANST, DLAPY2
      EXTERNAL           DLAMCH, DLANST, DLAPY2
c     ..
c     .. External Subroutines ..
      EXTERNAL           DLAE2, DLASCL, DLASRT, XERBLA
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          ABS, SIGN, SQRT
c     ..
c     .. Executable Statements ..
c
c     Test the input parameters.
c
      INFO = 0
c
c     Quick return if possible
c
      IF( N.LT.0 ) THEN
         INFO = -1
         CALL XERBLA( 'DSTERF', -INFO )
         RETURN
      END IF
      IF( N.LE.1 )
     $   RETURN
c
c     Determine the unit roundoff for this environment.
c
      EPS = DLAMCH( 'E' )
      EPS2 = EPS**2
      SAFMIN = DLAMCH( 'S' )
      SAFMAX = ONE / SAFMIN
      SSFMAX = SQRT( SAFMAX ) / THREE
      SSFMIN = SQRT( SAFMIN ) / EPS2
      RMAX = DLAMCH( 'O' )
c
c     Compute the eigenvalues of the tridiagonal matrix.
c
      NMAXIT = N*MAXIT
      SIGMA = ZERO
      JTOT = 0
c
c     Determine where the matrix splits and choose QL or QR iteration
c     for each block, according to whether top or bottom diagonal
c     element is smaller.
c
      L1 = 1
c
   10 CONTINUE
      IF( L1.GT.N )
     $   GO TO 170
      IF( L1.GT.1 )
     $   E( L1-1 ) = ZERO
      DO 20 M = L1, N - 1
         IF( ABS( E( M ) ).LE.( SQRT( ABS( D( M ) ) )*SQRT( ABS( D( M+
     $       1 ) ) ) )*EPS ) THEN
            E( M ) = ZERO
            GO TO 30
         END IF
   20 CONTINUE
      M = N
c
   30 CONTINUE
      L = L1
      LSV = L
      LEND = M
      LENDSV = LEND
      L1 = M + 1
      IF( LEND.EQ.L )
     $   GO TO 10
c
c     Scale submatrix in rows and columns L to LEND
c
      ANORM = DLANST( 'M', LEND-L+1, D( L ), E( L ) )
      ISCALE = 0
      IF( ANORM.EQ.ZERO )
     $   GO TO 10      
      IF( (ANORM.GT.SSFMAX) ) THEN
         ISCALE = 1
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L+1, 1, D( L ), N,
     $                INFO )
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L, 1, E( L ), N,
     $                INFO )
      ELSE IF( ANORM.LT.SSFMIN ) THEN
         ISCALE = 2
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L+1, 1, D( L ), N,
     $                INFO )
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L, 1, E( L ), N,
     $                INFO )
      END IF
c
      DO 40 I = L, LEND - 1
         E( I ) = E( I )**2
   40 CONTINUE
c
c     Choose between QL and QR iteration
c
      IF( ABS( D( LEND ) ).LT.ABS( D( L ) ) ) THEN
         LEND = LSV
         L = LENDSV
      END IF
c
      IF( LEND.GE.L ) THEN
c
c        QL Iteration
c
c        Look for small subdiagonal element.
c
   50    CONTINUE
         IF( L.NE.LEND ) THEN
            DO 60 M = L, LEND - 1
               IF( ABS( E( M ) ).LE.EPS2*ABS( D( M )*D( M+1 ) ) )
     $            GO TO 70
   60       CONTINUE
         END IF
         M = LEND
c
   70    CONTINUE
         IF( M.LT.LEND )
     $      E( M ) = ZERO
         P = D( L )
         IF( M.EQ.L )
     $      GO TO 90
c
c        If remaining matrix is 2 by 2, use DLAE2 to compute its
c        eigenvalues.
c
         IF( M.EQ.L+1 ) THEN
            RTE = SQRT( E( L ) )
            CALL DLAE2( D( L ), RTE, D( L+1 ), RT1, RT2 )
            D( L ) = RT1
            D( L+1 ) = RT2
            E( L ) = ZERO
            L = L + 2
            IF( L.LE.LEND )
     $         GO TO 50
            GO TO 150
         END IF
c
         IF( JTOT.EQ.NMAXIT )
     $      GO TO 150
         JTOT = JTOT + 1
c
c        Form shift.
c
         RTE = SQRT( E( L ) )
         SIGMA = ( D( L+1 )-P ) / ( TWO*RTE )
         R = DLAPY2( SIGMA, ONE )
         SIGMA = P - ( RTE / ( SIGMA+SIGN( R, SIGMA ) ) )
c
         C = ONE
         S = ZERO
         GAMMA = D( M ) - SIGMA
         P = GAMMA*GAMMA
c
c        Inner loop
c
         DO 80 I = M - 1, L, -1
            BB = E( I )
            R = P + BB
            IF( I.NE.M-1 )
     $         E( I+1 ) = S*R
            OLDC = C
            C = P / R
            S = BB / R
            OLDGAM = GAMMA
            ALPHA = D( I )
            GAMMA = C*( ALPHA-SIGMA ) - S*OLDGAM
            D( I+1 ) = OLDGAM + ( ALPHA-GAMMA )
            IF( C.NE.ZERO ) THEN
               P = ( GAMMA*GAMMA ) / C
            ELSE
               P = OLDC*BB
            END IF
   80    CONTINUE
c
         E( L ) = S*P
         D( L ) = SIGMA + GAMMA
         GO TO 50
c
c        Eigenvalue found.
c
   90    CONTINUE
         D( L ) = P
c
         L = L + 1
         IF( L.LE.LEND )
     $      GO TO 50
         GO TO 150
c
      ELSE
c
c        QR Iteration
c
c        Look for small superdiagonal element.
c
  100    CONTINUE
         DO 110 M = L, LEND + 1, -1
            IF( ABS( E( M-1 ) ).LE.EPS2*ABS( D( M )*D( M-1 ) ) )
     $         GO TO 120
  110    CONTINUE
         M = LEND
c
  120    CONTINUE
         IF( M.GT.LEND )
     $      E( M-1 ) = ZERO
         P = D( L )
         IF( M.EQ.L )
     $      GO TO 140
c
c        If remaining matrix is 2 by 2, use DLAE2 to compute its
c        eigenvalues.
c
         IF( M.EQ.L-1 ) THEN
            RTE = SQRT( E( L-1 ) )
            CALL DLAE2( D( L ), RTE, D( L-1 ), RT1, RT2 )
            D( L ) = RT1
            D( L-1 ) = RT2
            E( L-1 ) = ZERO
            L = L - 2
            IF( L.GE.LEND )
     $         GO TO 100
            GO TO 150
         END IF
c
         IF( JTOT.EQ.NMAXIT )
     $      GO TO 150
         JTOT = JTOT + 1
c
c        Form shift.
c
         RTE = SQRT( E( L-1 ) )
         SIGMA = ( D( L-1 )-P ) / ( TWO*RTE )
         R = DLAPY2( SIGMA, ONE )
         SIGMA = P - ( RTE / ( SIGMA+SIGN( R, SIGMA ) ) )
c
         C = ONE
         S = ZERO
         GAMMA = D( M ) - SIGMA
         P = GAMMA*GAMMA
c
c        Inner loop
c
         DO 130 I = M, L - 1
            BB = E( I )
            R = P + BB
            IF( I.NE.M )
     $         E( I-1 ) = S*R
            OLDC = C
            C = P / R
            S = BB / R
            OLDGAM = GAMMA
            ALPHA = D( I+1 )
            GAMMA = C*( ALPHA-SIGMA ) - S*OLDGAM
            D( I ) = OLDGAM + ( ALPHA-GAMMA )
            IF( C.NE.ZERO ) THEN
               P = ( GAMMA*GAMMA ) / C
            ELSE
               P = OLDC*BB
            END IF
  130    CONTINUE
c
         E( L-1 ) = S*P
         D( L ) = SIGMA + GAMMA
         GO TO 100
c
c        Eigenvalue found.
c
  140    CONTINUE
         D( L ) = P
c
         L = L - 1
         IF( L.GE.LEND )
     $      GO TO 100
         GO TO 150
c
      END IF
c
c     Undo scaling if necessary
c
  150 CONTINUE
      IF( ISCALE.EQ.1 )
     $   CALL DLASCL( 'G', 0, 0, SSFMAX, ANORM, LENDSV-LSV+1, 1,
     $                D( LSV ), N, INFO )
      IF( ISCALE.EQ.2 )
     $   CALL DLASCL( 'G', 0, 0, SSFMIN, ANORM, LENDSV-LSV+1, 1,
     $                D( LSV ), N, INFO )
c
c     Check for no convergence to an eigenvalue after a total
c     of N*MAXIT iterations.
c
      IF( JTOT.LT.NMAXIT )
     $   GO TO 10
      DO 160 I = 1, N - 1
         IF( E( I ).NE.ZERO )
     $      INFO = INFO + 1
  160 CONTINUE
      GO TO 180
c
c     Sort eigenvalues in increasing order.
c
  170 CONTINUE
      CALL DLASRT( 'I', N, D, INFO )
c
  180 CONTINUE
      RETURN
c
c     End of DSTERF
c
      END

      SUBROUTINE DORG2R( M, N, K, A, LDA, TAU, WORK, INFO )
c
c  -- LAPACK computational routine (version 3.4.2) --
c  -- LAPACK is a software package provided by Univ. of Tennessee,    --
c  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
c     September 2012
c
c     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, M, N
c     ..
c     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
c     ..
c
c  =====================================================================
c
c     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
c     ..
c     .. Local Scalars ..
      INTEGER            I, J, L
c     ..
c     .. External Subroutines ..
      EXTERNAL           DLARF, DSCAL, XERBLA
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          MAX
c     ..
c     .. Executable Statements ..
c
c     Test the input arguments
c
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORG2R', -INFO )
         RETURN
      END IF
c
c     Quick return if possible
c
      IF( N.LE.0 )
     $   RETURN
c
c     Initialise columns k+1:n to columns of the unit matrix
c
      DO 20 J = K + 1, N
         DO 10 L = 1, M
            A( L, J ) = ZERO
   10    CONTINUE
         A( J, J ) = ONE
   20 CONTINUE
c
      DO 40 I = K, 1, -1
c
c        Apply H(i) to A(i:m,i:n) from the left
c
         IF( I.LT.N ) THEN
            A( I, I ) = ONE
            CALL DLARF( 'Left', M-I+1, N-I, A( I, I ), 1, TAU( I ),
     $                  A( I, I+1 ), LDA, WORK )
         END IF
         IF( I.LT.M )
     $      CALL DSCAL( M-I, -TAU( I ), A( I+1, I ), 1 )
         A( I, I ) = ONE - TAU( I )
c
c        Set A(1:i-1,i) to zero
c
         DO 30 L = 1, I - 1
            A( L, I ) = ZERO
   30    CONTINUE
   40 CONTINUE
      RETURN
c
c     End of DORG2R
c
      END

c-----------------------------------------------------------------------
c---  END
c-----------------------------------------------------------------------