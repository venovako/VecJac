      PROGRAM EHYPOT
      USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_FMA, IEEE_MAX_NUM
      IMPLICIT NONE
      REAL ZERO, ONE
      PARAMETER (ZERO = 0.0, ONE = 1.0)
      COMPLEX C
      REAL CR, CI, A, H
      INTRINSIC ABS, HYPOT, MAX, MIN, SQRT
      WRITE (*,'(A)',ADVANCE='NO') 'C%RE='
      READ (*,*) CR
      WRITE (*,'(A)',ADVANCE='NO') 'C%IM='
      READ (*,*) CI
      C = CMPLX(CR, CI)
      A = ABS(C)
      WRITE (*,1) '          ABS(C)=', A
      H = HYPOT(CR, CI)
      WRITE (*,1) 'HYPOT(C%RE,C%IM)=', H
      CR = ABS(CR)
      CI = ABS(CI)
      A = MAX(CR, CI)
      H = MIN(CR, CI)
      H = IEEE_MAX_NUM(H / A, ZERO)
      H = SQRT(IEEE_FMA(H, H, ONE)) * A
      WRITE (*,1) 'NAIVE(C%RE,C%IM)=', H
 1    FORMAT(A,ES15.9)
      END
