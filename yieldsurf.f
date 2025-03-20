      SUBROUTINE JACOBI(A,N,NP,D,V,NROT)

!  Purpose: Computes all eigenvalues and eigenvectors of a real
!     symmetric matrix A, which is of size N by N, stored in a
!     physical NP by NP array.  On output, elements of A above the
!     diagonal are destroyed.  D returns the eigenvalues of A in
!     its first N elements.  V is a matrix with the same logical and
!     physical dimensions as A whose columns contain, on output, the
!     normalized eigenvectors of A.  NROT returns the number of Jacobi
!     rotations which were required. 
!
!  Source: W. H. Press et al., "Numerical Recipes", 1989, p. 346.
!
!  Modifications:
!
!     1. Double precision version
!
!  Prepared by J. Applequist, 10/23/91

      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NMAX=100)
      DIMENSION A(NP,NP),D(NP),V(NP,NP),B(NMAX),Z(NMAX)

!     Initialize the identity matrix.

      DO 12 IP=1,N
      DO 11 IQ=1,N
      V(IP,IQ)=0.D0
 11   CONTINUE
      V(IP,IP)=1.D0
 12   CONTINUE

!     Initialize B and D to the diagonal of A.

      DO 13 IP=1,N
      B(IP)=A(IP,IP)
      D(IP)=B(IP)
      Z(IP)=0.D0
 13   CONTINUE
      NROT=0
      DO 24 I=1,50
      SM=0.D0

!     Sum off-diagonal elements.

      DO 15 IP=1,N-1
      DO 14 IQ=IP+1,N
      SM=SM+DABS(A(IP,IQ))
 14   CONTINUE
 15   CONTINUE
      IF (SM.EQ.0.D0) RETURN
      IF (I.LT.4) THEN
      TRESH=0.2D0*SM/N**2
      ELSE
      TRESH=0.D0
      ENDIF
      DO 22 IP=1,N-1
      DO 21 IQ=IP+1,N
      G=100.D0*DABS(A(IP,IQ))

!     After four sweeps, skip the rotation if the off-diagonal
!     element is small.

      IF ((I.GT.4).AND.(DABS(D(IP))+G.EQ.DABS(D(IP))) &
       .AND.(DABS(D(IQ))+G.EQ.DABS(D(IQ)))) THEN
      A(IP,IQ)=0.D0
      ELSE IF (DABS(A(IP,IQ)).GT.TRESH) THEN
        H=D(IQ)-D(IP)
        IF (DABS(H)+G.EQ.DABS(H)) THEN
        T=A(IP,IQ)/H
        ELSE
        THETA=0.5D0*H/A(IP,IQ)
        T=1.D0/(DABS(THETA)+DSQRT(1.D0+THETA**2))
        IF (THETA.LT.0.D0) T=-T
        ENDIF
      C=1.D0/DSQRT(1.D0+T**2)
      S=T*C
      TAU=S/(1.D0+C)
      H=T*A(IP,IQ)
      Z(IP)=Z(IP)-H
      Z(IQ)=Z(IQ)+H
      D(IP)=D(IP)-H
      D(IQ)=D(IQ)+H
      A(IP,IQ)=0.D0
      DO 16 J=1,IP-1
      G=A(J,IP)
      H=A(J,IQ)
      A(J,IP)=G-S*(H+G*TAU)
      A(J,IQ)=H+S*(G-H*TAU)
 16   CONTINUE
      DO 17 J=IP+1,IQ-1
      G=A(IP,J)
      H=A(J,IQ)
      A(IP,J)=G-S*(H+G*TAU)
      A(J,IQ)=H+S*(G-H*TAU)
 17   CONTINUE
      DO 18 J=IQ+1,N
      G=A(IP,J)
      H=A(IQ,J)
      A(IP,J)=G-S*(H+G*TAU)
      A(IQ,J)=H+S*(G-H*TAU)
 18   CONTINUE
      DO 19 J=1,N
      G=V(J,IP)
      H=V(J,IQ)
      V(J,IP)=G-S*(H+G*TAU)
      V(J,IQ)=H+S*(G-H*TAU)
 19   CONTINUE
      NROT=NROT+1
      ENDIF
 21   CONTINUE
 22   CONTINUE
      DO 23 IP=1,N
      B(IP)=B(IP)+Z(IP)
      D(IP)=B(IP)
      Z(IP)=0.D0
 23   CONTINUE
 24   CONTINUE
      WRITE (6,600) 
 600  FORMAT(/'50 ITERATIONS OCCURRED IN SUBROUTINE JACOBI.')
      RETURN
      END SUBROUTINE JACOBI
