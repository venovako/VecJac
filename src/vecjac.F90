#ifndef VSL
#define VSL 16
#endif
#ifndef VDL
#define VDL 8
#endif
MODULE VECJAC
  USE, INTRINSIC :: ISO_C_BINDING
  IMPLICIT NONE

  ! character arrays must be nul-terminated

  INTERFACE
     FUNCTION OPEN_RO(BN, EX) BIND(C,NAME='open_ro_')
       USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_CHAR, C_INT
       IMPLICIT NONE
       CHARACTER(KIND=C_CHAR), INTENT(IN), TARGET :: BN(*), EX(*)
       INTEGER(KIND=C_INT) :: OPEN_RO
     END FUNCTION OPEN_RO
  END INTERFACE

  INTERFACE
     FUNCTION OPEN_WO(BN, EX) BIND(C,NAME='open_wo_')
       USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_CHAR, C_INT
       IMPLICIT NONE
       CHARACTER(KIND=C_CHAR), INTENT(IN), TARGET :: BN(*), EX(*)
       INTEGER(KIND=C_INT) :: OPEN_WO
     END FUNCTION OPEN_WO
  END INTERFACE

  INTERFACE
     FUNCTION RESIZEF(FD, SZ) BIND(C,NAME='resizef_')
       USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT, C_SIZE_T
       IMPLICIT NONE
       INTEGER(KIND=C_INT), INTENT(IN) :: FD
       INTEGER(KIND=C_SIZE_T), INTENT(IN) :: SZ
       INTEGER(KIND=C_INT) :: RESIZEF
     END FUNCTION RESIZEF
  END INTERFACE

  INTERFACE
     FUNCTION SREAD2(M, N, A, LDA, FD) BIND(C,NAME='sread2_')
       USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: M, N, LDA
       REAL, INTENT(OUT) :: A(LDA,N)
       INTEGER(KIND=C_INT), INTENT(IN) :: FD
       INTEGER(KIND=C_INT) :: SREAD2
     END FUNCTION SREAD2
  END INTERFACE

  INTERFACE
     FUNCTION DREAD2(M, N, A, LDA, FD) BIND(C,NAME='dread2_')
       USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: M, N, LDA
       DOUBLE PRECISION, INTENT(OUT) :: A(LDA,N)
       INTEGER(KIND=C_INT), INTENT(IN) :: FD
       INTEGER(KIND=C_INT) :: DREAD2
     END FUNCTION DREAD2
  END INTERFACE

  INTERFACE
     FUNCTION CREAD2(M, N, A, LDA, FD) BIND(C,NAME='cread2_')
       USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: M, N, LDA
       COMPLEX, INTENT(OUT) :: A(LDA,N)
       INTEGER(KIND=C_INT), INTENT(IN) :: FD
       INTEGER(KIND=C_INT) :: CREAD2
     END FUNCTION CREAD2
  END INTERFACE

  INTERFACE
     FUNCTION ZREAD2(M, N, A, LDA, FD) BIND(C,NAME='zread2_')
       USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: M, N, LDA
       DOUBLE COMPLEX, INTENT(OUT) :: A(LDA,N)
       INTEGER(KIND=C_INT), INTENT(IN) :: FD
       INTEGER(KIND=C_INT) :: ZREAD2
     END FUNCTION ZREAD2
  END INTERFACE

  INTERFACE
     FUNCTION WWRITE1(M, W, FD) BIND(C,NAME='wwrite1_')
       USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT
       USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: REAL128
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: M
       REAL(KIND=REAL128), INTENT(IN) :: W(M)
       INTEGER(KIND=C_INT), INTENT(IN) :: FD
       INTEGER(KIND=C_INT) :: WWRITE1
     END FUNCTION WWRITE1
  END INTERFACE

  INTERFACE
     FUNCTION SWRITE2(M, N, A, LDA, FD) BIND(C,NAME='swrite2_')
       USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: M, N, LDA
       REAL, INTENT(IN) :: A(LDA,N)
       INTEGER(KIND=C_INT), INTENT(IN) :: FD
       INTEGER(KIND=C_INT) :: SWRITE2
     END FUNCTION SWRITE2
  END INTERFACE

  INTERFACE
     FUNCTION DWRITE2(M, N, A, LDA, FD) BIND(C,NAME='dwrite2_')
       USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: M, N, LDA
       DOUBLE PRECISION, INTENT(IN) :: A(LDA,N)
       INTEGER(KIND=C_INT), INTENT(IN) :: FD
       INTEGER(KIND=C_INT) :: DWRITE2
     END FUNCTION DWRITE2
  END INTERFACE

  INTERFACE
     FUNCTION CWRITE2(M, N, A, LDA, FD) BIND(C,NAME='cwrite2_')
       USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: M, N, LDA
       COMPLEX, INTENT(IN) :: A(LDA,N)
       INTEGER(KIND=C_INT), INTENT(IN) :: FD
       INTEGER(KIND=C_INT) :: CWRITE2
     END FUNCTION CWRITE2
  END INTERFACE

  INTERFACE
     FUNCTION ZWRITE2(M, N, A, LDA, FD) BIND(C,NAME='zwrite2_')
       USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: M, N, LDA
       DOUBLE COMPLEX, INTENT(IN) :: A(LDA,N)
       INTEGER(KIND=C_INT), INTENT(IN) :: FD
       INTEGER(KIND=C_INT) :: ZWRITE2
     END FUNCTION ZWRITE2
  END INTERFACE

  INTERFACE
     INTEGER FUNCTION SALLOC2(M, N, A, LDA) BIND(C,NAME='salloc2_')
       USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: M, N
       INTEGER, INTENT(OUT) :: LDA
       TYPE(C_PTR), INTENT(OUT) :: A
     END FUNCTION SALLOC2
  END INTERFACE

  INTERFACE
     INTEGER FUNCTION DALLOC2(M, N, A, LDA) BIND(C,NAME='dalloc2_')
       USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: M, N
       INTEGER, INTENT(OUT) :: LDA
       TYPE(C_PTR), INTENT(OUT) :: A
     END FUNCTION DALLOC2
  END INTERFACE

  INTERFACE
     INTEGER FUNCTION CALLOC2(M, N, A, LDA, AR, LDAR, AI, LDAI) BIND(C,NAME='calloc2_')
       USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: M, N
       INTEGER, INTENT(OUT) :: LDA, LDAR, LDAI
       TYPE(C_PTR), INTENT(OUT) :: A, AR, AI
     END FUNCTION CALLOC2
  END INTERFACE

  INTERFACE
     INTEGER FUNCTION ZALLOC2(M, N, A, LDA, AR, LDAR, AI, LDAI) BIND(C,NAME='zalloc2_')
       USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: M, N
       INTEGER, INTENT(OUT) :: LDA, LDAR, LDAI
       TYPE(C_PTR), INTENT(OUT) :: A, AR, AI
     END FUNCTION ZALLOC2
  END INTERFACE

  INTERFACE
     SUBROUTINE CZFREE(A) BIND(C,NAME='czfree_')
       USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR
       IMPLICIT NONE
       TYPE(C_PTR), INTENT(INOUT) :: A
     END SUBROUTINE CZFREE
  END INTERFACE

  ! ASSUMES (MOD(M,VSL) .EQ. 0)
  INTERFACE
     REAL FUNCTION SDPSCL(M, X, Y, E, F) BIND(C,NAME='sdpscl_')
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: M
       REAL, INTENT(IN) :: X(MAX(M,VSL)), Y(MAX(M,VSL)), E(2), F(2)
     END FUNCTION SDPSCL
  END INTERFACE

  ! ASSUMES (MOD(M,VDL) .EQ. 0)
  INTERFACE
     DOUBLE PRECISION FUNCTION DDPSCL(M, X, Y, E, F) BIND(C,NAME='ddpscl_')
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: M
       DOUBLE PRECISION, INTENT(IN) :: X(MAX(M,VDL)), Y(MAX(M,VDL)), E(2), F(2)
     END FUNCTION DDPSCL
  END INTERFACE

  ! ASSUMES (MOD(M,VSL) .EQ. 0)
  INTERFACE
     COMPLEX FUNCTION CDPSCL(M, XR, XI, YR, YI, E, F) BIND(C,NAME='cdpscl_')
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: M
       REAL, INTENT(IN) :: XR(MAX(M,VSL)), XI(MAX(M,VSL)), YR(MAX(M,VSL)), YI(MAX(M,VSL)), E(2), F(2)
     END FUNCTION CDPSCL
  END INTERFACE

  ! ASSUMES (MOD(M,VDL) .EQ. 0)
  INTERFACE
     DOUBLE COMPLEX FUNCTION ZDPSCL(M, XR, XI, YR, YI, E, F) BIND(C,NAME='zdpscl_')
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: M
       DOUBLE PRECISION, INTENT(IN) :: XR(MAX(M,VDL)), XI(MAX(M,VDL)), YR(MAX(M,VDL)), YI(MAX(M,VDL)), E(2), F(2)
     END FUNCTION ZDPSCL
  END INTERFACE

  ! ASSUMES (MOD(M,VSL) .EQ. 0)
  INTERFACE
     REAL FUNCTION SGSSCL(M, T, X, Y, E, F) BIND(C,NAME='sgsscl_')
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: M
       REAL, INTENT(IN) :: T, E(2), F(2)
       REAL, INTENT(INOUT) :: X(MAX(ABS(M),VSL)), Y(MAX(ABS(M),VSL))
     END FUNCTION SGSSCL
  END INTERFACE

  ! ASSUMES (MOD(M,VDL) .EQ. 0)
  INTERFACE
     DOUBLE PRECISION FUNCTION DGSSCL(M, T, X, Y, E, F) BIND(C,NAME='dgsscl_')
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: M
       DOUBLE PRECISION, INTENT(IN) :: T, E(2), F(2)
       DOUBLE PRECISION, INTENT(INOUT) :: X(MAX(ABS(M),VSL)), Y(MAX(ABS(M),VSL))
     END FUNCTION DGSSCL
  END INTERFACE

  ! ASSUMES (MOD(M,VSL) .EQ. 0)
  INTERFACE
     REAL FUNCTION CGSSCL(M, TR, TI, XR, XI, YR, YI, E, F) BIND(C,NAME='cgsscl_')
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: M
       REAL, INTENT(IN) :: TR, TI, E(2), F(2)
       REAL, INTENT(INOUT) :: XR(MAX(ABS(M),VSL)), XI(MAX(ABS(M),VSL)), YR(MAX(ABS(M),VSL)), YI(MAX(ABS(M),VSL))
     END FUNCTION CGSSCL
  END INTERFACE

  ! ASSUMES (MOD(M,VDL) .EQ. 0)
  INTERFACE
     DOUBLE PRECISION FUNCTION ZGSSCL(M, TR, TI, XR, XI, YR, YI, E, F) BIND(C,NAME='zgsscl_')
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: M
       DOUBLE PRECISION, INTENT(IN) :: TR, TI, E(2), F(2)
       DOUBLE PRECISION, INTENT(INOUT) :: XR(MAX(ABS(M),VSL)), XI(MAX(ABS(M),VSL)), YR(MAX(ABS(M),VSL)), YI(MAX(ABS(M),VSL))
     END FUNCTION ZGSSCL
  END INTERFACE

  ! ASSUMES (MOD(N,VSL) .EQ. 0)
  INTERFACE
     INTEGER FUNCTION SBJAC2(N, A11, A22, A21, C, AT, L1, L2, P) BIND(C,NAME='sbjac2_')
       USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: N
       REAL, INTENT(IN) :: A11(MAX(N,VSL)), A22(MAX(N,VSL)), A21(MAX(N,VSL))
       REAL, INTENT(OUT) :: C(MAX(N,VSL)), AT(MAX(N,VSL)), L1(MAX(N,VSL)), L2(MAX(N,VSL))
       INTEGER(KIND=C_INT), INTENT(OUT) :: P(MAX(N,VSL)/VSL)
     END FUNCTION SBJAC2
  END INTERFACE

  ! ASSUMES (MOD(N,VDL) .EQ. 0)
  INTERFACE
     INTEGER FUNCTION DBJAC2(N, A11, A22, A21, C, AT, L1, L2, P) BIND(C,NAME='dbjac2_')
       USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: N
       DOUBLE PRECISION, INTENT(IN) :: A11(MAX(N,VDL)), A22(MAX(N,VDL)), A21(MAX(N,VDL))
       DOUBLE PRECISION, INTENT(OUT) :: C(MAX(N,VDL)), AT(MAX(N,VDL)), L1(MAX(N,VDL)), L2(MAX(N,VDL))
       INTEGER(KIND=C_INT), INTENT(OUT) :: P(MAX(N,VDL)/VDL)
     END FUNCTION DBJAC2
  END INTERFACE

  ! ASSUMES MOD(N,VSL) .EQ. 0
  INTERFACE
     INTEGER FUNCTION CBJAC2(N, A11, A22, A21R, A21I, C, CAT, SAT, L1, L2, P) BIND(C,NAME='cbjac2_')
       USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: N
       REAL, INTENT(IN) :: A11(MAX(N,VSL)), A22(MAX(N,VSL)), A21R(MAX(N,VSL)), A21I(MAX(N,VSL))
       REAL, INTENT(OUT) :: C(MAX(N,VSL)), CAT(MAX(N,VSL)), SAT(MAX(N,VSL)), L1(MAX(N,VSL)), L2(MAX(N,VSL))
       INTEGER(KIND=C_INT), INTENT(OUT) :: P(MAX(N,VSL)/VSL)
     END FUNCTION CBJAC2
  END INTERFACE

  ! ASSUMES MOD(N,VDL) .EQ. 0
  INTERFACE
     INTEGER FUNCTION ZBJAC2(N, A11, A22, A21R, A21I, C, CAT, SAT, L1, L2, P) BIND(C,NAME='zbjac2_')
       USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: N
       DOUBLE PRECISION, INTENT(IN) :: A11(MAX(N,VDL)), A22(MAX(N,VDL)), A21R(MAX(N,VDL)), A21I(MAX(N,VDL))
       DOUBLE PRECISION, INTENT(OUT) :: C(MAX(N,VDL)), CAT(MAX(N,VDL)), SAT(MAX(N,VDL)), L1(MAX(N,VDL)), L2(MAX(N,VDL))
       INTEGER(KIND=C_INT), INTENT(OUT) :: P(MAX(N,VDL)/VDL)
     END FUNCTION ZBJAC2
  END INTERFACE

  ! ASSUMES MOD(N,VSL) .EQ. 0
  INTERFACE
     REAL FUNCTION SJROT(N, X, Y, C, AT) BIND(C,NAME='sjrot_')
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: N
       REAL, INTENT(IN) :: C, AT
       REAL, INTENT(INOUT) :: X(MAX(ABS(N),VSL)), Y(MAX(ABS(N),VSL))
     END FUNCTION SJROT
  END INTERFACE

  ! ASSUMES MOD(N,VDL) .EQ. 0
  INTERFACE
     DOUBLE PRECISION FUNCTION DJROT(N, X, Y, C, AT) BIND(C,NAME='djrot_')
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: N
       DOUBLE PRECISION, INTENT(IN) :: C, AT
       DOUBLE PRECISION, INTENT(INOUT) :: X(MAX(ABS(N),VDL)), Y(MAX(ABS(N),VDL))
     END FUNCTION DJROT
  END INTERFACE

  ! ASSUMES MOD(N,VSL) .EQ. 0
  INTERFACE
     REAL FUNCTION CJROT(N, XR, XI, YR, YI, C, CAT, SAT) BIND(C,NAME='cjrot_')
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: N
       REAL, INTENT(IN) :: C, CAT, SAT
       REAL, INTENT(INOUT) :: XR(MAX(ABS(N),VSL)), XI(MAX(ABS(N),VSL)), YR(MAX(ABS(N),VSL)), YI(MAX(ABS(N),VSL))
     END FUNCTION CJROT
  END INTERFACE

  ! ASSUMES MOD(N,VDL) .EQ. 0
  INTERFACE
     DOUBLE PRECISION FUNCTION ZJROT(N, XR, XI, YR, YI, C, CAT, SAT) BIND(C,NAME='zjrot_')
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: N
       DOUBLE PRECISION, INTENT(IN) :: C, CAT, SAT
       DOUBLE PRECISION, INTENT(INOUT) :: XR(MAX(ABS(N),VDL)), XI(MAX(ABS(N),VDL)), YR(MAX(ABS(N),VDL)), YI(MAX(ABS(N),VDL))
     END FUNCTION ZJROT
  END INTERFACE

  ! ASSUMES MOD(N,VSL) .EQ. 0
  INTERFACE
     REAL FUNCTION SJROTF(N, X, Y, C, AT) BIND(C,NAME='sjrotf_')
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: N
       REAL, INTENT(IN) :: C, AT
       REAL, INTENT(INOUT) :: X(MAX(ABS(N),VSL)), Y(MAX(ABS(N),VSL))
     END FUNCTION SJROTF
  END INTERFACE

  ! ASSUMES MOD(N,VDL) .EQ. 0
  INTERFACE
     DOUBLE PRECISION FUNCTION DJROTF(N, X, Y, C, AT) BIND(C,NAME='djrotf_')
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: N
       DOUBLE PRECISION, INTENT(IN) :: C, AT
       DOUBLE PRECISION, INTENT(INOUT) :: X(MAX(ABS(N),VDL)), Y(MAX(ABS(N),VDL))
     END FUNCTION DJROTF
  END INTERFACE

  ! ASSUMES MOD(N,VSL) .EQ. 0
  INTERFACE
     REAL FUNCTION CJROTF(N, XR, XI, YR, YI, C, CAT, SAT) BIND(C,NAME='cjrotf_')
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: N
       REAL, INTENT(IN) :: C, CAT, SAT
       REAL, INTENT(INOUT) :: XR(MAX(ABS(N),VSL)), XI(MAX(ABS(N),VSL)), YR(MAX(ABS(N),VSL)), YI(MAX(ABS(N),VSL))
     END FUNCTION CJROTF
  END INTERFACE

  ! ASSUMES MOD(N,VDL) .EQ. 0
  INTERFACE
     DOUBLE PRECISION FUNCTION ZJROTF(N, XR, XI, YR, YI, C, CAT, SAT) BIND(C,NAME='zjrotf_')
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: N
       DOUBLE PRECISION, INTENT(IN) :: C, CAT, SAT
       DOUBLE PRECISION, INTENT(INOUT) :: XR(MAX(N,VDL)), XI(MAX(N,VDL)), YR(MAX(N,VDL)), YI(MAX(N,VDL))
     END FUNCTION ZJROTF
  END INTERFACE

  INTERFACE
     INTEGER FUNCTION SSWP(N, X, Y) BIND(C,NAME='sswp_')
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: N
       REAL, INTENT(INOUT) :: X(MAX(N,VSL)), Y(MAX(N,VSL))
     END FUNCTION SSWP
  END INTERFACE

  INTERFACE
     INTEGER FUNCTION DSWP(N, X, Y) BIND(C,NAME='dswp_')
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: N
       DOUBLE PRECISION, INTENT(INOUT) :: X(MAX(N,VDL)), Y(MAX(N,VDL))
     END FUNCTION DSWP
  END INTERFACE

  INTERFACE
     REAL FUNCTION SNORM2(M, X, E0, F0, E1, F1) BIND(C,NAME='snrom2_')
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: M
       REAL, INTENT(IN) :: X(MAX(M,0))
       REAL, INTENT(OUT) :: E0, F0, E1, F1
     END FUNCTION SNORM2
  END INTERFACE

  INTERFACE
     DOUBLE PRECISION FUNCTION DNORM2(M, X, E0, F0, E1, F1) BIND(C,NAME='dnrom2_')
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: M
       DOUBLE PRECISION, INTENT(IN) :: X(MAX(M,0))
       DOUBLE PRECISION, INTENT(OUT) :: E0, F0, E1, F1
     END FUNCTION DNORM2
  END INTERFACE

  INTERFACE
     REAL FUNCTION CNORM2(M, ZR, ZI, E0, F0, E1, F1) BIND(C,NAME='cnrom2_')
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: M
       REAL, INTENT(IN) :: ZR(MAX(M,0)), ZI(MAX(M,0))
       REAL, INTENT(OUT) :: E0, F0, E1, F1
     END FUNCTION CNORM2
  END INTERFACE

  INTERFACE
     DOUBLE PRECISION FUNCTION ZNORM2(M, ZR, ZI, E0, F0, E1, F1) BIND(C,NAME='znrom2_')
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: M
       DOUBLE PRECISION, INTENT(IN) :: ZR(MAX(M,0)), ZI(MAX(M,0))
       DOUBLE PRECISION, INTENT(OUT) :: E0, F0, E1, F1
     END FUNCTION ZNORM2
  END INTERFACE

#ifdef USE_SLEEF
  INTERFACE
     DOUBLE PRECISION FUNCTION DNORMF(M, X, E0, F0, E1, F1) BIND(C,NAME='dnormf_')
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: M
       DOUBLE PRECISION, INTENT(IN) :: X(MAX(M,0))
       DOUBLE PRECISION, INTENT(OUT) :: E0, F0, E1, F1
     END FUNCTION DNORMF
  END INTERFACE

  INTERFACE
     DOUBLE PRECISION FUNCTION DNORME(M, X, E0, F0, E1, F1) BIND(C,NAME='dnorme_')
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: M
       DOUBLE PRECISION, INTENT(IN) :: X(MAX(M,VDL))
       DOUBLE PRECISION, INTENT(OUT) :: E0, F0, E1, F1
     END FUNCTION DNORME
  END INTERFACE

  INTERFACE
     DOUBLE PRECISION FUNCTION DNORMS(M, X, E0, F0, E1, F1) BIND(C,NAME='dnorms_')
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: M
       DOUBLE PRECISION, INTENT(IN) :: X(MAX(M,VDL))
       DOUBLE PRECISION, INTENT(OUT) :: E0, F0, E1, F1
     END FUNCTION DNORMS
  END INTERFACE
#endif

  INTERFACE
     REAL FUNCTION SNORMX(M, N, A, LDA) BIND(C,NAME='snormx_')
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: M, N, LDA
       REAL, INTENT(IN) :: A(MAX(LDA,VSL),MAX(N,0))
     END FUNCTION SNORMX
  END INTERFACE

  INTERFACE
     DOUBLE PRECISION FUNCTION DNORMX(M, N, A, LDA) BIND(C,NAME='dnormx_')
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: M, N, LDA
       DOUBLE PRECISION, INTENT(IN) :: A(MAX(LDA,VDL),MAX(N,0))
     END FUNCTION DNORMX
  END INTERFACE

  INTERFACE
     REAL FUNCTION CNORMX(M, N, AR, LDAR, AI, LDAI) BIND(C,NAME='cnormx_')
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: M, N, LDAR, LDAI
       REAL, INTENT(IN) :: AR(MAX(LDAR,VSL),MAX(N,0)), AI(MAX(LDAI,VSL),MAX(N,0))
     END FUNCTION CNORMX
  END INTERFACE

  INTERFACE
     DOUBLE PRECISION FUNCTION ZNORMX(M, N, AR, LDAR, AI, LDAI) BIND(C,NAME='znormx_')
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: M, N, LDAR, LDAI
       DOUBLE PRECISION, INTENT(IN) :: AR(MAX(LDAR,VDL),MAX(N,0)), AI(MAX(LDAI,VDL),MAX(N,0))
     END FUNCTION ZNORMX
  END INTERFACE

  INTERFACE
     INTEGER FUNCTION SSCALE(M, N, A, LDA, E) BIND(C,NAME='sscale_')
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: M, N, LDA, E
       REAL, INTENT(INOUT) :: A(MAX(LDA,VSL),MAX(N,0))
     END FUNCTION SSCALE
  END INTERFACE

  INTERFACE
     INTEGER FUNCTION DSCALE(M, N, A, LDA, E) BIND(C,NAME='dscale_')
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: M, N, LDA, E
       DOUBLE PRECISION, INTENT(INOUT) :: A(MAX(LDA,VDL),MAX(N,0))
     END FUNCTION DSCALE
  END INTERFACE

  INTERFACE
     INTEGER FUNCTION CSCALE(M, N, AR, LDAR, AI, LDAI, E) BIND(C,NAME='cscale_')
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: M, N, LDAR, LDAI, E
       REAL, INTENT(INOUT) :: AR(MAX(LDAR,VSL),MAX(N,0)), AI(MAX(LDAI,VSL),MAX(N,0))
     END FUNCTION CSCALE
  END INTERFACE

  INTERFACE
     INTEGER FUNCTION ZSCALE(M, N, AR, LDAR, AI, LDAI, E) BIND(C,NAME='zscale_')
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: M, N, LDAR, LDAI, E
       DOUBLE PRECISION, INTENT(INOUT) :: AR(MAX(LDAR,VDL),MAX(N,0)), AI(MAX(LDAI,VDL),MAX(N,0))
     END FUNCTION ZSCALE
  END INTERFACE

  ! !!! FOR xVJSVD ROUTINES: !!!
  ! if JTRACE is defined, WORK should start with a nul-terminated name of the trace file to be written
  ! set IWORK(1) to 0 if V is to be set to I, or to something else if V is preset

  ! ASSUMES (MOD(N,2) .EQ. 0) .AND. (MOD(N/2,VSL) .EQ. 0)
  INTERFACE
     INTEGER FUNCTION SVJSVD(M, N, G, LDG, V, LDV, ES, FS, JS, STP, SWP, WORK, IWORK) BIND(C,NAME='svjsvd_')
       USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: M, N, LDG, LDV, JS(*)
       INTEGER(KIND=C_INT), INTENT(IN) :: STP, SWP
       REAL, INTENT(INOUT) :: G(MAX(LDG,VSL),MAX(N,0))
       REAL, INTENT(OUT) :: V(MAX(LDV,VSL),MAX(N,0))
       REAL, INTENT(OUT) :: ES, FS, WORK(5*MAX(N,0))
       INTEGER(KIND=C_INT), INTENT(INOUT) :: IWORK(MAX(N,0)/VSL)
     END FUNCTION SVJSVD
  END INTERFACE

  ! ASSUMES (MOD(N,2) .EQ. 0) .AND. (MOD(N/2,VDL) .EQ. 0)
  INTERFACE
     INTEGER FUNCTION DVJSVD(M, N, G, LDG, V, LDV, ES, FS, JS, STP, SWP, WORK, IWORK) BIND(C,NAME='dvjsvd_')
       USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: M, N, LDG, LDV, JS(*)
       INTEGER(KIND=C_INT), INTENT(IN) :: STP, SWP
       DOUBLE PRECISION, INTENT(INOUT) :: G(MAX(LDG,VDL),MAX(N,0))
       DOUBLE PRECISION, INTENT(OUT) :: V(MAX(LDV,VDL),MAX(N,0))
       DOUBLE PRECISION, INTENT(OUT) :: ES, FS, WORK(5*MAX(N,0))
       INTEGER(KIND=C_INT), INTENT(INOUT) :: IWORK(MAX(N,0)/VDL)
     END FUNCTION DVJSVD
  END INTERFACE

  ! ASSUMES (MOD(N,2) .EQ. 0) .AND. (MOD(N/2,VSL) .EQ. 0)
  INTERFACE
     INTEGER FUNCTION CVJSVD(M,N, GR,LDGR,GI,LDGI, VR,LDVR,VI,LDVI, ES,FS, JS,STP,SWP, WORK,IWORK) BIND(C,NAME='cvjsvd_')
       USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: M, N, LDGR, LDGI, LDVR, LDVI, JS(*)
       INTEGER(KIND=C_INT), INTENT(IN) :: STP, SWP
       REAL, INTENT(INOUT) :: GR(MAX(LDGR,VSL),MAX(N,0)), GI(MAX(LDGI,VSL),MAX(N,0))
       REAL, INTENT(OUT) :: VR(MAX(LDVR,VSL),MAX(N,0)), VI(MAX(LDVI,VSL),MAX(N,0))
       REAL, INTENT(OUT) :: ES, FS, WORK(7*MAX(N,0))
       INTEGER(KIND=C_INT), INTENT(INOUT) :: IWORK(MAX(N,0)/VSL)
     END FUNCTION CVJSVD
  END INTERFACE

  ! ASSUMES (MOD(N,2) .EQ. 0) .AND. (MOD(N/2,VDL) .EQ. 0)
  INTERFACE
     INTEGER FUNCTION ZVJSVD(M,N, GR,LDGR,GI,LDGI, VR,LDVR,VI,LDVI, ES,FS, JS,STP,SWP, WORK,IWORK) BIND(C,NAME='zvjsvd_')
       USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: M, N, LDGR, LDGI, LDVR, LDVI, JS(*)
       INTEGER(KIND=C_INT), INTENT(IN) :: STP, SWP
       DOUBLE PRECISION, INTENT(INOUT) :: GR(MAX(LDGR,VDL),MAX(N,0)), GI(MAX(LDGI,VDL),MAX(N,0))
       DOUBLE PRECISION, INTENT(OUT) :: VR(MAX(LDVR,VDL),MAX(N,0)), VI(MAX(LDVI,VDL),MAX(N,0))
       DOUBLE PRECISION, INTENT(OUT) :: ES, FS, WORK(7*MAX(N,0))
       INTEGER(KIND=C_INT), INTENT(INOUT) :: IWORK(MAX(N,0)/VDL)
     END FUNCTION ZVJSVD
  END INTERFACE

  INTERFACE
     INTEGER FUNCTION CMERGE(M, N, AR, LDAR, AI, LDAI, A, LDA) BIND(C,NAME='cmerge_')
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: M, N, LDAR, LDAI, LDA
       REAL, INTENT(IN) :: AR(MAX(LDAR,VSL),MAX(N,0)), AI(MAX(LDAI,VSL),MAX(N,0))
       COMPLEX, INTENT(OUT) :: A(MAX(LDA,VSL/2),MAX(N,0))
     END FUNCTION CMERGE
  END INTERFACE

  INTERFACE
     INTEGER FUNCTION ZMERGE(M, N, AR, LDAR, AI, LDAI, A, LDA) BIND(C,NAME='zmerge_')
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: M, N, LDAR, LDAI, LDA
       DOUBLE PRECISION, INTENT(IN) :: AR(MAX(LDAR,VDL),MAX(N,0)), AI(MAX(LDAI,VDL),MAX(N,0))
       DOUBLE COMPLEX, INTENT(OUT) :: A(MAX(LDA,VDL/2),MAX(N,0))
     END FUNCTION ZMERGE
  END INTERFACE

  INTERFACE
     INTEGER FUNCTION CSPLIT(M, N, A, LDA, AR, LDAR, AI, LDAI) BIND(C,NAME='csplit_')
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: M, N, LDA, LDAR, LDAI
       COMPLEX, INTENT(IN) :: A(MAX(LDA,VSL/2),MAX(N,0))
       REAL, INTENT(OUT) :: AR(MAX(LDAR,VSL),MAX(N,0)), AI(MAX(LDAI,VSL),MAX(N,0))
     END FUNCTION CSPLIT
  END INTERFACE

  INTERFACE
     INTEGER FUNCTION ZSPLIT(M, N, A, LDA, AR, LDAR, AI, LDAI) BIND(C,NAME='zsplit_')
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: M, N, LDA, LDAR, LDAI
       DOUBLE COMPLEX, INTENT(IN) :: A(MAX(LDA,VDL/2),MAX(N,0))
       DOUBLE PRECISION, INTENT(OUT) :: AR(MAX(LDAR,VDL),MAX(N,0)), AI(MAX(LDAI,VDL),MAX(N,0))
     END FUNCTION ZSPLIT
  END INTERFACE

CONTAINS
END MODULE VECJAC
