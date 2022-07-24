      PROGRAM ELAEV2
      IMPLICIT NONE
      REAL ZERO, ONE
      PARAMETER (ZERO = 0.0, ONE = 1.0)
      REAL A_R, B_R, B_I, C_R, RT1, RT2, CS1
      COMPLEX A, B, C, SN1, D(2,2)
      EXTERNAL CLAEV2
      INTRINSIC ABS, AIMAG, CMPLX, CONJG, REAL
      WRITE (*,'(A)',ADVANCE='NO') '   A='
      READ (*,*) A_R
      WRITE (*,'(A)',ADVANCE='NO') 'B%RE='
      READ (*,*) B_R
      WRITE (*,'(A)',ADVANCE='NO') 'B%IM='
      READ (*,*) B_I
      WRITE (*,'(A)',ADVANCE='NO') '   C='
      READ (*,*) C_R
      A = CMPLX(A_R, ZERO)
      B = CMPLX(B_R, B_I)
      C = CMPLX(C_R, ZERO)
      CALL CLAEV2(A, B, C, RT1, RT2, CS1, SN1)
      WRITE (*,1) 'RT1=', RT1
      WRITE (*,1) 'RT2=', RT2
      WRITE (*,1) 'CS1=', CS1
      WRITE (*,2) 'SN1=', SN1
      A_R = CS1 * CS1
      B_R = REAL(SN1)
      B_I = AIMAG(SN1)
      C_R = B_R * B_R + B_I * B_I
      D(1,1) = RT1 * A_R + RT2 * C_R
      D(2,1) = SN1 * (CS1 * (RT1 - RT2))
      D(1,2) = CONJG(D(2,1))
      D(2,2) = RT1 * C_R + RT2 * A_R
      WRITE (*,2) 'D(1,1)=', D(1,1)
      WRITE (*,2) 'D(2,1)=', D(2,1)
      WRITE (*,2) 'D(1,2)=', D(1,2)
      WRITE (*,2) 'D(2,2)=', D(2,2)
      WRITE (*,1) 'ORT(U)=', ((C_R + A_R) - ONE)
 1    FORMAT(A,ES16.9)
 2    FORMAT(A,'(',ES16.9,',',ES16.9,')')
      END
