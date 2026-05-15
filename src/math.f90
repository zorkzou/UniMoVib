!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! distance between points a and b
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function distance(a,b)
 implicit real(kind=8) (a-h,o-z)
 real(kind=8) :: a(3),b(3),c(3)

 call ASub(3,a,b,c)
 distance=sqrt( dotx(3,c,c) )

 return
end function distance

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! bond angle a-b-c (0~pi rad.)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function BAngle(a,b,c)
 implicit real(kind=8) (a-h,o-z)
 parameter(One=1.d0)
 real(kind=8) :: a(3),b(3),c(3),tmp(3,2)

 call ASub(3,a,b,tmp(1,1))
 call ASub(3,c,b,tmp(1,2))
 d12=sqrt( dotx(3,tmp(1,1),tmp(1,1)) )
 d32=sqrt( dotx(3,tmp(1,2),tmp(1,2)) )
 call AScale(3,One/d12,tmp(1,1),tmp(1,1))
 call AScale(3,One/d32,tmp(1,2),tmp(1,2))
 cphi=dotx(3,tmp(1,1),tmp(1,2))
 ! because of numerical error, |cphi| can be larger than 1
 if(abs(cphi) > One) cphi = sign(One,cphi)
 BAngle=acos(cphi)

 return
end function BAngle

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! C(*) = A(*) - B(*)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine ASub(N,A,B,C)
 implicit real(kind=8) (a-h,o-z)
 real(kind=8) :: A(N),B(N),C(N)

 !Do I = 1,N
 !  C(I) = A(I) - B(I)
 !end do
 C = A - B

 return
end subroutine ASub

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! C(*) = A(*) + B(*)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine AAdd(N,A,B,C)
 implicit real(kind=8) (a-h,o-z)
 real(kind=8) :: A(N),B(N),C(N)

 !Do I = 1,N
 !  C(I) = A(I) + B(I)
 !end do
 C = A + B

 return
end subroutine AAdd

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! generate a unit matrix
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine UMAT(N,U)
 implicit real(kind=8) (a-h,o-z)

 parameter(one=1.d0)
 real(kind=8) :: U(N,N)

 call AClear(N*N,U)
 do I = 1,N
   U(i,i)=one
 end do

 return
end subroutine UMAT

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! B(*) = c * A(*)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine AScale(N,c,A,B)
 implicit real(kind=8) (a-h,o-z)
 real(kind=8) :: A(N),B(N)

 !Do I = 1,N
 !  B(I) = c*A(I)
 !end do
 B = c * A

 return
end subroutine AScale

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! symmetric square matrix --> lower triangular matrix
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Sq2Tr(N,S,T)
 implicit real(kind=8) (a-h,o-z)
 real(kind=8) :: T(*),S(N,N)

 ii=0
 Do i=1,N
   Do j=1,i
     ii=ii+1
     T(ii)=(S(j,i)+S(i,j))*0.5d0
   end Do
 end Do

 return
end subroutine Sq2Tr

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! lower triangular matrix --> symmetric square matrix
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine LT2Sqr(N,T,S)
 implicit real(kind=8) (a-h,o-z)
 real(kind=8) :: T(*),S(N,N)

 k=0
 do i=1,N
   do j=1,i-1
     k=k+1
     S(j,i)=T(k)
     S(i,j)=T(k)
   end do
   k=k+1
   S(i,i)=T(k)
 end do

 return
end subroutine LT2Sqr

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! It returns the maximum element of array A. IMax is the index of the largest element.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ArMax(N,IMax,A)
 implicit real(kind=8) (a-h,o-z)
 real(kind=8) :: A(N)
 Save Zero
 Data Zero/0.0d0/

 If(N < 1) then
   ArMax = Zero
   IMax = 0
 else
   IMax = 1
   Do I = 2, N
     If(A(I) > A(IMax)) IMax = I
   end do
   ArMax = A(IMax)
 endIf

 return
end function ArMax

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! It returns the minimum element of array A. IMin is the index of the smallest element.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Function ArMin(N,IMin,A)
 implicit real(kind=8) (a-h,o-z)
 real(kind=8) :: A(N)
 Save Zero
 Data Zero/0.0d0/

 If(N < 1) then
   ArMin = Zero
   IMin = 0
 else
   IMin = 1
   Do I = 2, N
     If(A(I) < A(IMin)) IMin = I
   end do
   ArMin = A(IMin)
 endIf

 return
end Function ArMin

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! vector A dot_product vector B
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dotx(N,A,B)
 implicit real(kind=8) (a-h,o-z)
 real(kind=8) :: A(N),B(N)

 dotx = dot_product(A,B)

 return
end function dotx

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! calculate the 2nd derivatives of the nuclear repulsion energy: ENR = Sum_i,j(Z_i * Z_j / R_i,j)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine DDerNRE(NAtm,Z,xyz,FFx,DAT,D2)
 implicit real(kind=8) (a-h,o-z)
 real(kind=8) :: Z(*),xyz(3,*),FFx(NAtm*3,NAtm*3),DAT(5),D2(6)

 NAtm3 = NAtm*3
 call AClear(NAtm3*NAtm3,FFx)

 do i=1,NAtm
   ix = 3*(i-1)+1
   iy = 3*(i-1)+2
   iz = 3*(i-1)+3
   do j=1,i-1
     jx = 3*(j-1)+1
     jy = 3*(j-1)+2
     jz = 3*(j-1)+3

 !   (1) x_i - x_j
 !   (2) y_i - y_j
 !   (3) z_i - z_j
 !   (4) Z_i * Z_j / r_i,j^3
 !   (5) 3 * Z_i * Z_j / r_i,j^5
     DAT(1) = xyz(1,i) - xyz(1,j)
     DAT(2) = xyz(2,i) - xyz(2,j)
     DAT(3) = xyz(3,i) - xyz(3,j)
     RR = dotx(3,DAT,DAT)
     DAT(4) = Z(i) * Z(j) / (RR * sqrt(RR))
     DAT(5) = DAT(4) * 3.d0 / RR

 !   2nd derivative terms:
 !   dx_i*dx_j,
 !   dx_i*dy_j, dy_i*dy_j,
 !   dx_i*dz_j, dy_i*dz_j, dz_i*dz_j
     D2(1) = -DAT(5)*DAT(1)*DAT(1) + DAT(4)
     D2(2) = -DAT(5)*DAT(1)*DAT(2)
     D2(3) = -DAT(5)*DAT(2)*DAT(2) + DAT(4)
     D2(4) = -DAT(5)*DAT(1)*DAT(3)
     D2(5) = -DAT(5)*DAT(2)*DAT(3)
     D2(6) = -DAT(5)*DAT(3)*DAT(3) + DAT(4)

 !   2nd derivatives
     FFx(ix,jx) = D2(1)
     FFx(ix,jy) = D2(2)
     FFx(ix,jz) = D2(4)
     FFx(iy,jx) = D2(2)
     FFx(iy,jy) = D2(3)
     FFx(iy,jz) = D2(5)
     FFx(iz,jx) = D2(4)
     FFx(iz,jy) = D2(5)
     FFx(iz,jz) = D2(6)

     FFx(jx,ix) = D2(1)
     FFx(jy,ix) = D2(2)
     FFx(jz,ix) = D2(4)
     FFx(jx,iy) = D2(2)
     FFx(jy,iy) = D2(3)
     FFx(jz,iy) = D2(5)
     FFx(jx,iz) = D2(4)
     FFx(jy,iz) = D2(5)
     FFx(jz,iz) = D2(6)

   end do
 end do

 return
end subroutine DDerNRE

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Wrapper routine for DGEMM.
! Note! MAXN should be the same as NX.
! MATRIX MULTIPLICATION PACKAGE FOR SQUARE MATRICES ONLY.
! IQ=
! 1  C=A*B
! 2  C=A(TRANSPOSE)*B
! 3  C=A*B(TRANSPOSE)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine MPACMF(A,B,C,MAXN,NX,IQ)
 implicit real(kind=8) (a-h,o-z)
 real(kind=8) :: A(MAXN,MAXN),B(MAXN,MAXN),C(MAXN,MAXN)

 call AClear(NX*NX,C)                                         ! Zero the product

 if(IQ == 1) then
   call DGEMM('N','N',NX,NX,NX,1.d0,A,NX,B,NX,0.d0,C,NX)      ! C=A*B
 else if(IQ == 2) then
   call DGEMM('T','N',NX,NX,NX,1.d0,A,NX,B,NX,0.d0,C,NX)      ! C=A(TRANSPOSE)*B
 else if(IQ == 3) then
   call DGEMM('N','T',NX,NX,NX,1.d0,A,NX,B,NX,0.d0,C,NX)      ! C=A*B(TRANSPOSE)
 end if

 return
end subroutine MPACMF

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Calculate square-root matrix X=S^1/2 and/or its (general) inverse
! Xi. Here S must be a SYMMETRIC matrix.
! Mode > 0: Calculate X
!      < 0: Calculate Xi
!      = 0: Calculate both
! Scratch Scr(N,M), M = max(2*N,4)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine SqrtMp(Intact,Mode,N,S,X,Xi,E,A,Scr)
 implicit real(kind=8) (a-h,o-z)
 real(kind=8) :: S(N,N),X(N,N),Xi(N,N),E(N,N),A(N),Scr(N,*)
 logical :: Intact

 call ACopy(N*N,S,E)
 LWork=2*N*max(N,2)
 call DSYEV('V','L', N, E, N, A, Scr, LWork, Info)
 if(INFO /= 0) call XError(Intact,"SqrtMp")

 if(Mode == 0)then
   call pDiagSq(N,0,X,A)
   call MPACMF(E,X,Scr,N,N,1)
   call MPACMF(Scr,E,X,N,N,3)
   call pDiagSq(N,1,Xi,A)
   call MPACMF(E,Xi,Scr,N,N,1)
   call MPACMF(Scr,E,Xi,N,N,3)
 else if(Mode > 0)then
   call pDiagSq(N,0,X,A)
   call MPACMF(E,X,Scr,N,N,1)
   call MPACMF(Scr,E,X,N,N,3)
 else
   call pDiagSq(N,1,Xi,A)
   call MPACMF(E,Xi,Scr,N,N,1)
   call MPACMF(Scr,E,Xi,N,N,3)
 end if

 return
end subroutine SqrtMp

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! General inverse enhanced version of DiagSqrt.
! A is a diagonal matrix with diagonal terms
! Indx = 0, sqrt(B)
!     /= 0, sqrt(B)**-1 if B(i) > 0.
! It's assumed that B(i) >= 0. It doesn't work if there are negative elements in B, which is not checked.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine pDiagSq(N,Indx,A,B)
 implicit real(kind=8) (a-h,o-z)
 parameter(eps=1.d-12)
 real(kind=8) :: A(N,N),B(N)

 call AClear(N*N,A)
 if(Indx==0) then
   Do i=1,N
     A(i,i)=sqrt(abs(B(i)))
   endDo
 else
   Do i=1,N
     if(abs(B(i)) > eps) A(i,i)=1.d0/sqrt(abs(B(i)))
   endDo
 endIf

 return
end subroutine pDiagSq

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! B(*) = A(*)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine ACopy(N,A,B)
 implicit real(kind=8) (a-h,o-z)
 real(kind=8) :: A(N),B(N)

 !Do I = 1,N
 !  B(I) = A(I)
 !end do
 B = A

 return
end subroutine ACopy

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Routine to clear N elements in array A.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine AClear(N,A)
 implicit real(kind=8) (a-h,o-z)
 real(kind=8) :: A(N)

 !Do I = 1,N
 !  A(I) = 0.0d0
 !end do
 A = 0.0d0

 return
end subroutine AClear

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Routine to clear N characters
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine CClear(N,CA)
 implicit real(kind=8) (a-h,o-z)
 Character*(*) :: CA

 Do I = 1,N
   CA(I:I) = " "
 end do

 return
end subroutine CClear

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! X is (general) inverse of square matrix S.
! ISYMM = 1: S is symmetric
!         0: S is not symmetric
!        -1: S is symmetric, but S should be symmetrized first
! Scr(N,M): scratch, M = max(2*N,4)
! S and X can be the same matrix.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine GInvM(Intact,ISYMM,N,S,X,E,A,Scr)
 implicit real(kind=8) (a-h,o-z)
 real(kind=8) :: S(N,N),X(N,N),E(N,N),A(N),Scr(N,N,*)
 logical :: Intact

 if(ISYMM == 0)then
 ! E = S * S^T will be diagonalized
   call MPACMF(S,S,E,N,N,3)
 else if(ISYMM > 0)then
 ! E = S will be diagonalized
   call ACopy(N*N,S,E)
 else
 ! E = (S + S^T)/2 will be diagonalized
   call Symtrz(N,S,E)
 end if
 LWork=2*N*max(N,2)
 call DSYEV('V','L', N, E, N, A, Scr, LWork, Info)
 if(INFO /= 0) call XError(Intact,"GInvM")

 call DiagIv(N,Scr(1,1,1),A)
 call MPACMF(E,Scr(1,1,1),Scr(1,1,2),N,N,1)
 if(ISYMM == 0)then
 ! Scr(:,:,1) = Inv(S * S^T)
   call MPACMF(Scr(1,1,2),E,Scr(1,1,1),N,N,3)
 ! X = S^T * Inv(S * S^T)
   call MPACMF(S,Scr(1,1,1),Scr(1,1,2),N,N,2)
   call ACopy(N*N,Scr(1,1,2),X)
 else
   call MPACMF(Scr(1,1,2),E,X,N,N,3)
 end if

 return
end subroutine GInvM

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! A is a diagonal matrix with diagonal terms 1/B(i).
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine DiagIv(N,A,B)
 implicit real(kind=8) (a-h,o-z)
 parameter(eps=1.d-12)
 real(kind=8) :: A(N,N),B(N)

 call AClear(N*N,A)
 Do i=1,N
   if(abs(B(i)) > eps) A(i,i)=1.d0/B(i)
 endDo

 return
end subroutine DiagIv

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! symmetrize A by B = (A + A')/2
!
! A is (nearly) symmetric; A and B can be the same matrix
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Symtrz(N,A,B)
 implicit real(kind=8) (a-h,o-z)
 parameter(half=0.5d0)
 real(kind=8) :: A(N,N),B(N,N)

 do i=1,N
   do j=1,i-1
     B(j,i)=(A(j,i)+A(i,j))*half
     B(i,j)=B(j,i)
   end do
   B(i,i)=A(i,i)
 end do

 return
end subroutine Symtrz

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! save square S(3,3,N) to L.T. A(6,N)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine S9to6(N,S,A)
 implicit real(kind=8) (a-h,o-z)
 real(kind=8) :: A(6,N),S(9,N)

 do i=1,N
   call Sq2Tr(3,S(1,i),A(1,i))
 end do

 return
end subroutine S9to6

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! c = a - b, where a<0 or b<0 is an imaginary value
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine cplxsub(a,b,c)
 implicit real(kind=8) (a-h,o-z)
 parameter(Zero=0.0d0)
 real(kind=8) :: c(2)

 c = Zero
 if(a >= Zero) then
   c(1) = a
 else
   c(2) = a
 end if

 if(b >= Zero) then
   c(1) = c(1) - b
 else
   c(2) = c(2) - b
 end if

 return
end subroutine cplxsub

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Solve A * L = L * e, where A is symmetric, and L will be saved in A. The size of W must be 3N or larger.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine DiagS1(Intact,N,A,E,W)
 implicit real(kind=8) (a-h,o-z)
 real(kind=8) :: A(N,N),E(N),W(*)
 logical :: Intact

 LWORK=3*N
 Call DSYEV('V','L',N,A,N,E,W,LWORK,INFO)
 if(INFO /= 0) call XError(Intact,"DiagS1")

 return
end subroutine DiagS1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! B = A + A^T. A is symmetric; A and B can be the same matrix.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine TrAdd(N,A,B)
 Implicit Real*8(A-H,O-Z)
 real(kind=8) :: A(N,N),B(N,N)

 do i=1,N
   do j=1,i-1
     B(j,i)=A(j,i)+A(i,j)
     B(i,j)=B(j,i)
   end do
   B(i,i)=A(i,i)+A(i,i)
 end do

 return
end subroutine TrAdd

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! AT(N,M) = A(M,N)^T
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Transp(M,N,A,AT)
 implicit real(kind=8) (a-h,o-z)
 real(kind=8) :: A(M,N),AT(N,M)

 do i=1,N
   do j=1,M
     AT(i,j) = A(j,i)
   end do
 end do

 return
end subroutine Transp

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! C(I,J) = Sum(K) A(I,K) * B(K,J); C = A * B
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine MMpyMF(L,M,N,A,B,C)
 implicit real(kind=8) (a-h,o-z)
 real(kind=8) :: A(L,M), B(M,N), C(L,N)
 Save Zero,One
 Data Zero/0.0d0/,One/1.0d0/

 call AClear(L*N,C)
 call DGEMM('N','N',L,N,M,One,A,L,B,M,Zero,C,L)

 return
end subroutine MMpyMF

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! B = A - I; A and B can be the same matrix.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine MSubI(N,A,B)
 implicit real(kind=8) (a-h,o-z)
 real(kind=8) A(N,N),B(N,N)

 Do I = 1, N
   Do J = 1, N
     B(J,I) = A(J,I)
     if(I == J) B(J,I) = B(J,I) - 1.d0
   end do
 end do

 return
end subroutine MSubI

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! remove numerical noise
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RmNoise(N,tol,A)
 implicit real(kind=8) (a-h,o-z)
 real(kind=8) :: A(N)

 Do I = 1,N
   if(abs(A(I)) < tol) A(I)=0.0d0
 end do

 return
end subroutine RmNoise

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! B(*) = B(*) + c * A(*)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine AccAB(N,c,A,B)
 implicit real(kind=8) (a-h,o-z)
 real(kind=8) :: A(*),B(*)

 Do I = 1,N
   B(I) = B(I) + c*A(I)
 end do

 return
end subroutine AccAB

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! sum shell; if A is not a one-dimensional array, the intrinsic function Sum doesn't work.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ASum(A,N)
 implicit real(kind=8) (a-h,o-z)
 real(kind=8) :: A(N)

 ASum = Sum(A)

 return
end function ASum

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! multiplication; If IfRm0=.True., get rid of tiny values.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function AMultip(A,N,IfRm0,eps)
 implicit real(kind=8) (a-h,o-z)
 real(kind=8) :: A(N)
 logical :: IfRm0

 AMultip = 1.0d0
 do i=1,N
   if(IfRm0 .and. abs(A(i)) < eps) cycle
   AMultip = AMultip*A(i)
 end do

 return
end function AMultip

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! MULT33 performs a eulerian rotation on a symmetry operation.
!
! The symmetry operation stored in ELEM is subjected to the operation stored in FMAT: ELEM = FMAT * ELEM * FMAT^T
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine mult33(fmat,elem,SC1)
 implicit real(kind=8) (a-h,o-z)
 real(kind=8) :: fmat(3,3),elem(3,3),SC1(3,3)

 call AClear(9,SC1)
 do i = 1, 3
   do j = 1, 3
     do k = 1, 3
       do l = 1, 3
         SC1(i,j) = SC1(i,j) + fmat(i,l)*fmat(j,k)*elem(l,k)
       end do
     end do
   end do
 end do
 call Acopy(9,SC1,elem)

 return
end subroutine mult33

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! renormalization of an array
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine renorm(N,A)
 implicit real(kind=8) (a-h,o-z)
 parameter(one=1.0d0)
 real(kind=8) :: A(N)

 v = one/sqrt( dotx(N,A,A) )
 call ascale(N,v,A,A)

 return
end subroutine renorm

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! cross product c = a x b
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine CrossX(a,b,c)
 implicit real(kind=8) (a-h,o-z)
 real(kind=8) :: a(3),b(3),c(3)

 c(1) = a(2)*b(3) - a(3)*b(2)
 c(2) = a(3)*b(1) - a(1)*b(3)
 c(3) = a(1)*b(2) - a(2)*b(1)

 return
end subroutine CrossX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! It performs a 3x3 rotation operation:
! mode >=0: XYZ = XYZ * ROT
!      < 0: XYZ = XYZ * ROT^-1 = XYZ * ROT^T
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine rotopr(Natom,mode,ROT,XYZ,SCR)
 implicit real(kind=8) (a-h,o-z)
 real(kind=8) :: XYZ(3,Natom),ROT(3,3),SCR(3)

 if(mode >= 0) then
   do i=1,Natom
     call acopy(3,XYZ(1,i),SCR)
     call aclear(3,XYZ(1,i))
     do j=1,3
       do k=1,3
         XYZ(j,i)=XYZ(j,i)+ROT(k,j)*SCR(k)
       end do
     end do
   end do
 else
   do i=1,Natom
     call acopy(3,XYZ(1,i),SCR)
     call aclear(3,XYZ(1,i))
     do j=1,3
       do k=1,3
         XYZ(j,i)=XYZ(j,i)+ROT(j,k)*SCR(k)
       end do
     end do
   end do
 end if

 return
end subroutine rotopr

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Move the Cartesian coordinates to the geometric center.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine Shift(n,o,xyz)
 implicit real(kind=8) (a-h,o-z)
 dimension :: xyz(3,*), o(3)

 o = 0.0d0
 do i=1,n
   o = o + xyz(:,i)
 end do
 o = o /dble(n)
 do i=1,n
   xyz(:,i) = xyz(:,i) - o
 end do

 return
end Subroutine Shift

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! construct Kearsley's Q-matrix. See the last equation in
! S.K.Kearsley, On the orthogonal transformation used for structural comparisons, Acta Cryst. A45, 208 (1989)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine conqmt(iprint,iout,n,xp,xm,q)
 implicit real(kind=8) (a-h,o-z)
 dimension :: xp(3,n), xm(3,n), q(4,4)

 q = 0.0d0

 do i = 1, n
   q(1,1) = q(1,1) + xm(1,i)*xm(1,i) + xm(2,i)*xm(2,i) + xm(3,i)*xm(3,i)
   q(2,2) = q(2,2) + xp(2,i)*xp(2,i) + xp(3,i)*xp(3,i) + xm(1,i)*xm(1,i)
   q(3,3) = q(3,3) + xp(1,i)*xp(1,i) + xp(3,i)*xp(3,i) + xm(2,i)*xm(2,i)
   q(4,4) = q(4,4) + xp(1,i)*xp(1,i) + xp(2,i)*xp(2,i) + xm(3,i)*xm(3,i)

   q(2,1) = q(2,1) + xp(2,i)*xm(3,i) - xm(2,i)*xp(3,i)
   q(3,1) = q(3,1) + xm(1,i)*xp(3,i) - xp(1,i)*xm(3,i)
   q(4,1) = q(4,1) + xp(1,i)*xm(2,i) - xm(1,i)*xp(2,i)
   q(3,2) = q(3,2) + xm(1,i)*xm(2,i) - xp(1,i)*xp(2,i)
   q(4,2) = q(4,2) + xm(1,i)*xm(3,i) - xp(1,i)*xp(3,i)
   q(4,3) = q(4,3) + xm(2,i)*xm(3,i) - xp(2,i)*xp(3,i)

   q(1,2) = q(2,1)
   q(1,3) = q(3,1)
   q(1,4) = q(4,1)
   q(2,3) = q(3,2)
   q(2,4) = q(4,2)
   q(3,4) = q(4,3)
 end do

 if(iprint > 0) then
   write(iout,"(/,' Q-matrix:')")
   write(iout,"(4x,4f16.6)") q
 end if

 return
end subroutine conqmt

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Constructing the best fit rotation matrix. See the 2nd equation in
! S.K.Kearsley, On the orthogonal transformation used for structural comparisons, Acta Cryst. A45, 208 (1989)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine conrot(iprint,iout,q,r)
 implicit real(kind=8) (a-h,o-z)
 dimension         :: q(4), r(3,3)

 r(1,1) = q(1)*q(1) + q(2)*q(2) - q(3)*q(3) - q(4)*q(4)
 r(2,1) = 2 * (q(2)*q(3) + q(1)*q(4))
 r(3,1) = 2 * (q(2)*q(4) - q(1)*q(3))
 r(1,2) = 2 * (q(2)*q(3) - q(1)*q(4))
 r(2,2) = q(1)*q(1) - q(2)*q(2) + q(3)*q(3) - q(4)*q(4)
 r(3,2) = 2 * (q(3)*q(4) + q(1)*q(2))
 r(1,3) = 2 * (q(2)*q(4) + q(1)*q(3))
 r(2,3) = 2 * (q(3)*q(4) - q(1)*q(2))
 r(3,3) = q(1)*q(1) - q(2)*q(2) - q(3)*q(3) + q(4)*q(4)

 if(iprint > 0) then
   write(iout,"(/,' Rotation matrix for the best fit:')")
   do i = 1, 3
     write(iout,"(4x,3f16.6)") r(i,:)
   end do
 end if

 return
end subroutine conrot

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! It performs a 3x3 rotation operation of Cartesian vectors:
! mode >=0: XYZ' = XYZ * ROT
!      < 0: XYZ' = XYZ * ROT^-1 = XYZ * ROT^T
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine rotvec(NV,mode,ROT,XYZ,SCR)
 implicit double precision (a-h,o-z)
 dimension XYZ(3,NV),ROT(3,3),SCR(3)

 if(mode >= 0) then
   do i=1,NV
     SCR = XYZ(:,i)
     XYZ(:,i) = 0.0d0
     do j=1,3
       do k=1,3
         XYZ(j,i)=XYZ(j,i)+ROT(k,j)*SCR(k)
       end do
     end do
   end do
 else
   do i=1,NV
     SCR = XYZ(:,i)
     XYZ(:,i) = 0.0d0
     do j=1,3
       do k=1,3
         XYZ(j,i)=XYZ(j,i)+ROT(j,k)*SCR(k)
       end do
     end do
   end do
 end if

 return
end subroutine rotvec

!--- END

