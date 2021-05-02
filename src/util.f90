!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! obtain atomic indices
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine obt_idx(iout,flags,subsystem_idx,NAtm)
 implicit real(kind=8) (a-h,o-z)
 integer :: subsystem_idx(*),flags(*)

 J = 1
 Do I=1,NAtm
   if(flags(I) == 1)then
     subsystem_idx(J) = I
     J = J + 1
   end if
 End Do

 write(iout, '(10I5)'),(subsystem_idx(i),i=1,J-1)

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Take data of subsystem
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine PartMX(NAtm,AMass,XYZ,ZA,subsystem_idx,NAtm_sub,AMass_sub,XYZ_sub,ZA_sub)
 implicit real(kind=8) (a-h,o-z)
 real(kind=8) :: AMass(*),AMass_sub(NAtm_sub),XYZ(3,*),XYZ_sub(3,NAtm_sub),ZA(*),ZA_sub(NAtm_sub)
 integer :: subsystem_idx(*)

 Do I=1,NAtm_sub
    idx_p = subsystem_idx(I)
    AMass_sub(I) = AMass(idx_p)
    ZA_sub(I) = ZA(idx_p)
    Do J=1,3
       XYZ_sub(J,I) = XYZ(J,idx_p)
    End Do
 End Do

 !PRINT '(10F7.2)',(AMass_sub(i),i=1,3)
 !PRINT '(10F7.2)',(XYZ_sub(i,2),i=1,3)
 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Mode = 0 : returns nuclear charge zchar for an element symbol "el"; iza is not used.
!     /= 0 : returns element symbol "el" for nuclear charge iza; zchar is not used.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine ElemZA(Mode,el,iza,zchar)
 implicit real(kind=8) (a-h,o-z)
 parameter (maxza=120)
 character*3 :: el,atomlib(maxza)
 data (atomlib(i),i=1,maxza) / &
  'H  ','HE ','LI ','BE ','B  ','C  ','N  ','O  ','F  ','NE ',   'NA ','MG ','AL ','SI ','P  ','S  ','CL ','AR ','K  ','CA ', &
  'SC ','TI ','V  ','CR ','MN ','FE ','CO ','NI ','CU ','ZN ',   'GA ','GE ','AS ','SE ','BR ','KR ','RB ','SR ','Y  ','ZR ', &
  'NB ','MO ','TC ','RU ','RH ','PD ','AG ','CD ','IN ','SN ',   'SB ','TE ','I  ','XE ','CS ','BA ','LA ','CE ','PR ','ND ', &
  'PM ','SM ','EU ','GD ','TB ','DY ','HO ','ER ','TM ','YB ',   'LU ','HF ','TA ','W  ','RE ','OS ','IR ','PT ','AU ','HG ', &
  'TL ','PB ','BI ','PO ','AT ','RN ','FR ','RA ','AC ','TH ',   'PA ','U  ','NP ','PU ','AM ','CM ','BK ','CF ','ES ','FM ', &
  'MD ','NO ','LR ','RF ','DB ','SG ','BH ','HS ','MT ','DS ',   'RG ','CN ','NH ','FL ','MC ','LV ','TS ','OG ','UUE','UBN'/
 save atomlib

 if (Mode == 0) then

   call charl2u(el)
   zchar = 0.d0
   do i=1,maxza
     if(index(el,atomlib(i)) /= 0)then
       zchar = dble(i)
       exit
     end if
   end do

 else

   el = "???"
   if(iza > 0 .and. iza <= maxza) el = adjustl(atomlib(iza))
   call charu2l(el(2:3))

 end if

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Mode = 0 : mass of the most abundant isotope or the longest lived isotope
!     /= 0 : averaged mass of isotopes
! negative mass if IZ is out of range.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function EleMas(Mode,IZ)
 implicit real(kind=8) (a-h,o-z)
 parameter (maxza=120)
 real(kind=8) :: amas1(maxza),amas2(maxza)
 ! the most abundant or the longest lived isotopic masses from https://en.wikipedia.org/wiki/Extended_periodic_table
 data (amas1(i),i=1,maxza) /  &
    1.007825037d0,   4.002603250d0,   7.016004500d0,   9.012182500d0,  11.009305300d0,  12.000000000d0,&
   14.003074008d0,  15.994914640d0,  18.998403250d0,  19.992439100d0,  22.989769700d0,  23.985045000d0,&
   26.981541300d0,  27.976928400d0,  30.973763400d0,  31.972071800d0,  34.968852729d0,  39.962383100d0,&
   38.963707900d0,  39.962590700d0,  44.955913600d0,  47.947946700d0,  50.943962500d0,  51.940509700d0,&
   54.938046300d0,  55.934939300d0,  58.933197800d0,  57.935347100d0,  62.929599200d0,  63.929145400d0,&
   68.925580900d0,  73.921178800d0,  74.921595500d0,  79.916520500d0,  78.918336100d0,  83.911506400d0,&
   84.911700000d0,  87.905600000d0,  88.905400000d0,  89.904300000d0,  92.906000000d0,  97.905500000d0,&
   98.906300000d0, 101.903700000d0, 102.904800000d0, 105.903200000d0, 106.905090000d0, 113.903600000d0,&
  114.904100000d0, 117.901800000d0, 120.903800000d0, 129.906700000d0, 126.900400000d0, 131.904200000d0,&
  132.905429000d0, 137.905000000d0, 138.906100000d0, 139.905300000d0, 140.907400000d0, 141.907500000d0,&
  144.912700000d0, 151.919500000d0, 152.920900000d0, 157.924100000d0, 158.925000000d0, 163.928800000d0,&
  164.930300000d0, 165.930400000d0, 168.934400000d0, 173.939000000d0, 174.940900000d0, 179.946800000d0,&
  180.948000000d0, 183.951000000d0, 186.956000000d0, 189.958600000d0, 192.963300000d0, 194.964800000d0,&
  196.966600000d0, 201.970600000d0, 204.974500000d0, 207.976600000d0, 208.980400000d0, 208.982500000d0,&
  210.987500000d0, 222.017500000d0, 223.019800000d0, 226.025400000d0, 227.027800000d0, 232.038200000d0,&
  231.035900000d0, 238.050800000d0, 237.048000000d0, 244.064204000d0, 243.061381100d0, 247.070354000d0,&
  247.070307000d0, 251.079587000d0, 252.082980000d0, 257.095106000d0, 258.098431000d0, 259.101030000d0,&
  266.119830000d0, 267.121790000d0, 268.125670000d0, 269.128630000d0, 270.133360000d0, 269.133750000d0,&
  278.156310000d0, 281.164510000d0, 282.169120000d0, 285.177120000d0, 286.182210000d0, 289.190420000d0,&
  290.195980000d0, 293.204490000d0, 294.210460000d0, 295.216240000d0, 295.000000000d0, 299.000000000d0/
 ! averaged isotopic masses
 data (amas2(i),i=1,maxza) /  &
    1.007940000d0,   4.002602000d0,   6.941000000d0,   9.012183100d0,  10.811000000d0,  12.010700000d0,&
   14.006700000d0,  15.999400000d0,  18.998403163d0,  20.179700000d0,  22.989769280d0,  24.305000000d0,&
   26.981538500d0,  28.085500000d0,  30.973761998d0,  32.065000000d0,  35.453000000d0,  39.948000000d0,&
   39.098300000d0,  40.078000000d0,  44.955908000d0,  47.867000000d0,  50.941500000d0,  51.996100000d0,&
   54.938044000d0,  55.845000000d0,  58.933194000d0,  58.693400000d0,  63.546000000d0,  65.380000000d0,&
   69.723000000d0,  72.640000000d0,  74.921595000d0,  78.971000000d0,  79.904000000d0,  83.798000000d0,&
   85.467800000d0,  87.620000000d0,  88.905840000d0,  91.224000000d0,  92.906370000d0,  95.950000000d0,&
   98.907200000d0, 101.070000000d0, 102.905500000d0, 106.420000000d0, 107.868200000d0, 112.414000000d0,&
  114.818000000d0, 118.710000000d0, 121.760000000d0, 127.600000000d0, 126.904470000d0, 131.293000000d0,&
  132.905451960d0, 137.327000000d0, 138.905470000d0, 140.116000000d0, 140.907660000d0, 144.242000000d0,&
  144.900000000d0, 150.360000000d0, 151.964000000d0, 157.250000000d0, 158.925350000d0, 162.500000000d0,&
  164.930330000d0, 167.259000000d0, 168.934220000d0, 173.054000000d0, 174.966800000d0, 178.490000000d0,&
  180.947880000d0, 183.840000000d0, 186.207000000d0, 190.230000000d0, 192.217000000d0, 195.084000000d0,&
  196.966569000d0, 200.590000000d0, 204.383300000d0, 207.200000000d0, 208.980400000d0, 208.982400000d0,&
  209.987100000d0, 222.017600000d0, 223.019700000d0, 226.024500000d0, 227.027700000d0, 232.037700000d0,&
  231.035880000d0, 238.028910000d0, 237.048200000d0, 239.064200000d0, 243.061400000d0, 247.070400000d0,&
  247.070300000d0, 251.079600000d0, 252.083000000d0, 257.059100000d0, 258.098400000d0, 259.101000000d0,&
  262.109700000d0, 261.108800000d0, 262.114100000d0, 266.121900000d0, 264.120100000d0, 265.000000000d0,&
  268.138800000d0, 269.000000000d0, 272.000000000d0, 277.000000000d0, 280.000000000d0, 280.000000000d0,&
  280.000000000d0, 280.000000000d0, 280.000000000d0, 280.000000000d0, 280.000000000d0, 280.000000000d0/
 save amas1,amas2

 if(IZ < 1 .or. IZ > maxza) then

   EleMas = -1.d0

 else if (Mode == 0) then

   EleMas = amas1(IZ)

 else

   EleMas = amas2(IZ)

 end if

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! take atomic masses from library
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine MasLib(Mode,NAtm,AMass,ZA)
 implicit real(kind=8) (a-h,o-z)
 real(kind=8) :: AMass(*),ZA(*)

 do i=1,NAtm
   iza = nint(ZA(i))
   AMass(i) = EleMas(Mode,iza)
 end do

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! take the APT element(s) in the out of plane (ip=1) as the atomic IR charge
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function AIRCrg(ip,AtmAPT)
 implicit real(kind=8) (a-h,o-z)
 real(kind=8) :: AtmAPT(3,3),ip(3)

 AIRCrg = 0.0d0
 do i=1,3
   if(ip(i) == 1) AIRCrg = AIRCrg + AtmAPT(i,i)
 end do
 AIRCrg = AIRCrg / dble(sum(ip))

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! check whether the molecule is planar (or linear as a special case).
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine planar(NAtm,XYZ,ip)
 implicit real(kind=8) (a-h,o-z)
 parameter(tol=1.0d-6)
 real(kind=8) :: XYZ(3,*),ip(3)

 ip = 0

 do ix = 1, 3
   x = 0.0d0
   do i=1,NAtm
     x = x + abs(XYZ(ix,i))
   end do
   if(x <= tol) ip(ix) = 1
 end do

 return
end

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
end

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
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! replace numbers by space
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine rmnumb(N,cha)
 implicit real(kind=8) (a-h,o-z)
 character*(*) :: cha

 do i=1,N
   if((ichar(cha(i:i)) >= 48).and.(ichar(cha(i:i)) <= 57)) cha(i:i)=' '
 end do

 return
end

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
end

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
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! CHA --> cha
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine charu2l(cha)
 implicit real(kind=8) (a-h,o-z)
 character*(*) :: cha
 character*1  :: U2L

 do i=1,len_trim(cha)
   cha(i:i)=U2L(cha(i:i))
 end do

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! cha --> CHA
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine charl2u(cha)
 implicit real(kind=8) (a-h,o-z)
 character*(*) :: cha
 character*1  :: L2U

 do i=1,len_trim(cha)
   cha(i:i)=L2U(cha(i:i))
 end do

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! l --> L
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function L2U(letter)
 implicit real(kind=8) (a-h,o-z)
 character*1 :: letter,L2U

 if( ichar(letter) >= 97 .and. ichar(letter) <= 122 )then
   L2U=char(ichar(letter)-32)
 else
   L2U=letter
 endif

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! L --> l
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function U2L(letter)
 implicit real(kind=8) (a-h,o-z)
 character*1 :: letter,U2L

 if((ichar(letter) >= 65).and.(ichar(letter) <= 90))then
   U2L=char(ichar(letter)+32)
 else
   U2L=letter
 endif

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! position of the first non-space character in a string
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function nonspace(string)
 implicit real(kind=8) (a-h,o-z)
 character*(*) :: string

 length=LEN_TRIM(string)
 if(length <= 1) then
   i=length
 else
   do i=1,length
     if(string(i:i) /= ' ') exit
   end do
 endif

 nonspace=i

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! read an error message and stop
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine XError(Intact,inf)
 implicit real(kind=8) (a-h,o-z)
 logical :: Intact
 character*(*) :: inf

 write(*,"(/,' *** Error! ',a)")trim(inf)

 call estop(Intact)

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! read an <ENTER> and stop
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine estop(Intact)
 implicit real(kind=8) (a-h,o-z)
 logical :: Intact

 if(Intact) then
   write(*,"(/,1x,70('='),/,' Press <ENTER> to exit',/)")
   read(*,*)
 end if

 stop

 return
end

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
end

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
end

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
end

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
end

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
end

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
end

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
end

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
end

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

 call AClear(NX*NX,C)                                          ! Zero the product

 goto (10,20,30), IQ
 10 call DGEMM('N','N',NX,NX,NX,1.d0,A,NX,B,NX,0.d0,C,NX)      ! C=A*B
 return

 20 call DGEMM('T','N',NX,NX,NX,1.d0,A,NX,B,NX,0.d0,C,NX)      ! C=A(TRANSPOSE)*B
 return

 30 call DGEMM('N','T',NX,NX,NX,1.d0,A,NX,B,NX,0.d0,C,NX)      ! C=A*B(TRANSPOSE)
 return
end

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
end

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
end

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
end

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
end

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
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! open a data file
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine DFOpen(iout,ifile,Intact,ifstop,ftype,fname)
 implicit real(kind=8) (a-h,o-z)
 character*100 :: fname
 character*4 :: ftype
 logical :: Intact, ifstop

 ! check fname
 if(LEN_TRIM(fname) == 0)then
   if(ifstop)then
     write(*,"(/,1x,a4)")ftype
     call XError(Intact,"The above data file is not specified!")
   else
     write(iout,"(/,1x,a4,' is not specified!')")ftype
   end if
   return
 end if

 istr=nonspace(fname)
 iend=LEN_TRIM(fname)
 open(ifile,file=fname(istr:iend),status='old',err=100)
 write(iout,"(1x,a4,' file:',7x,a)")ftype,fname(istr:iend)
 if(Intact) write(*,"(1x,a4,' file:',7x,a)")ftype,fname(istr:iend)
 return

 100 continue
 if(ifstop)then
   write(*,"(/,1x,a4,' = ',a)")ftype,fname(istr:iend)
   call XError(Intact,"The above data file does not exist!")
 else
   write(iout,"(/,1x,a4,' is not specified!')")ftype
 end if

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Calculate force constant matrix using vibrational frequencies (Freq in a.u.) and normal modes by
!
!   FFX * LL = M * LL * EE
!   FFX = M * LL * EE * LL^T * (LL * LL^T)^-1
!   FFX = 0.5 * (FFX + FFX^T)
!
! LL and EE: full dimensional (3N x 3N) L and E (i.e. frequency square).
!
! AL saves LL, which should be mass-unweighted. The size of AL is NAtm3*NAtm3, that is, the rot. + trans. modes should also be
! included, but their element values can be all zero. If NVib < NAtm3, rot. + trans. modes should be after vib. modes.
!
! The size of WORK is max(2*NAtm3,4)*NAtm3.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Frq2FFX(NAtm3,NVib,AMass,Freq,AL,FFX,Scr,WORK,EIG)
 implicit real(kind=8) (a-h,o-z)
 real(kind=8) :: AMass(*),AL(NAtm3*NAtm3),Freq(NVib),FFX(*),Scr(NAtm3,NAtm3),WORK(*),EIG(NAtm3)

 ! LL^T * (LL * LL^T)^-1 --> FFX
 call GInvM(.false.,0,NAtm3,AL,FFX,Scr,EIG,WORK)

 ! M * LL * EE --> Scr
 call ACopy(NAtm3*NVib,AL,Scr)
 ! If NVib = NAtm3, the calculation is in full space. But if NVib < NAtm3, Scr(:,NVib+1:) are zero since rot. + trans. modes
 ! have been projected out (i.e. their eigenvalues are zero)
 call AClear((NAtm3-NVib)*NAtm3, Scr(1,NVib+1))
 do i=1,NVib
   ee=sign(Freq(i)*Freq(i), Freq(i))
   do j=1,NAtm3
     ja=(j-1)/3+1
     Scr(j,i)=Scr(j,i)*ee*AMass(ja)
   end do
 end do

 ! Scr * FFX --> WORK
 call MPACMF(Scr,FFX,WORK,NAtm3,NAtm3,1)

 ! symmetrize FFx
 call Symtrz(NAtm3,WORK,FFX)

 return
end

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
end

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
end

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
end

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
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! generate m.w. normal modes of translations and rotations
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine TRVec(Intact,NAtm,NTR,Imiss,AMass,XYZCM,AL,ROT,ALtmp)
 implicit real(kind=8) (a-h,o-z)
 parameter(tol=1.d-10)
 real(kind=8) :: AMass(*),XYZCM(3,*),AL(3,NAtm,*),ROT(3,3),ALtmp(3,NAtm,6)
 logical :: Intact

 NAtm3=3*NAtm
 call AClear(NAtm3*6,ALtmp)

 do i=1,NAtm
 ! rotate
   RX=dotx(3,XYZCM(1,i),ROT(1,1))
   RY=dotx(3,XYZCM(1,i),ROT(1,2))
   RZ=dotx(3,XYZCM(1,i),ROT(1,3))

   sqrtm=sqrt(AMass(i))
   ALtmp(1,I,1)=sqrtm
   ALtmp(2,I,2)=sqrtm
   ALtmp(3,I,3)=sqrtm
   ALtmp(1,I,4)=(RY*ROT(1,3)-RZ*ROT(1,2))*sqrtm
   ALtmp(2,I,4)=(RY*ROT(2,3)-RZ*ROT(2,2))*sqrtm
   ALtmp(3,I,4)=(RY*ROT(3,3)-RZ*ROT(3,2))*sqrtm
   ALtmp(1,I,5)=(RZ*ROT(1,1)-RX*ROT(1,3))*sqrtm
   ALtmp(2,I,5)=(RZ*ROT(2,1)-RX*ROT(2,3))*sqrtm
   ALtmp(3,I,5)=(RZ*ROT(3,1)-RX*ROT(3,3))*sqrtm
   ALtmp(1,I,6)=(RX*ROT(1,2)-RY*ROT(1,1))*sqrtm
   ALtmp(2,I,6)=(RX*ROT(2,2)-RY*ROT(2,1))*sqrtm
   ALtmp(3,I,6)=(RX*ROT(3,2)-RY*ROT(3,1))*sqrtm
 end do

 ! renormalization
 NTR=0
 Imiss = 0
 do i=1,6
   X = dotx(NAtm3,ALtmp(1,1,i),ALtmp(1,1,i))
   if(X > tol)then
     NTR=NTR+1
     X=1.d0/sqrt(X)
     call AScale(NAtm3,X,ALtmp(1,1,i),AL(1,1,NTR))
   else
     Imiss = i
   end if
 end do

 if(NTR /= 5 .and. NTR /= 6)then
   write(*,"(/,' NTR = ',i2)")NTR
   call XError(Intact,"Wrong NTR in TRVec!")
 end if
 if(NTR == 5 .and. Imiss < 4)then
   write(*,"(/,' Imiss = ',i2)")Imiss
   call XError(Intact,"Wrong Imiss in TRVec!")
 end if

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! principal moment of inertia
!
! The size of WORK should be 9 or larger.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine MIner(Intact,NAtm,AMass,XYZCM,RI,E,WORK)
 implicit real(kind=8) (a-h,o-z)
 real(kind=8) :: AMass(*),XYZCM(3,NAtm),RI(3,3),E(3),WORK(9)
 logical :: Intact

 call AClear(9,RI)
 do i=1,NAtm
   W=AMass(i)
   X=XYZCM(1,i)
   Y=XYZCM(2,i)
   Z=XYZCM(3,i)
   RI(1,1) = RI(1,1) + W * (Y*Y + Z*Z)
   RI(2,2) = RI(2,2) + W * (X*X + Z*Z)
   RI(3,3) = RI(3,3) + W * (X*X + Y*Y)
   RI(1,2) = RI(1,2) - W * X * Y
   RI(1,3) = RI(1,3) - W * X * Z
   RI(2,3) = RI(2,3) - W * Y * Z
 end do
 RI(2,1) = RI(1,2)
 RI(3,1) = RI(1,3)
 RI(3,2) = RI(2,3)

 call DiagS1(Intact,3,RI,E,WORK)

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Calculate the principal moment of inertia (m.w. is not shown)
! | yy+zz -xy   -xz   |                  | xx xy xz |
! | -xy   xx+zz -yz   | = (xx+yy+zz)*I - | xy yy yz |
! | -xz   -yz   xx+yy |                  | xz yz zz |
! ie, M(1) = (xx+yy+zz)*I - M(2)
! Different from subroutine MIner, MIner2 calculates M(2) which has some advantages. The number of zero eigenvalues (e1, e2, and
! e3) can be
!   3:  atom symmetry
!   2:  linear symmetry
!   1:  planar symmetry
!   0:  other symmetry
!       e1 = e2 = e3: cubic symmetry
!       e1 = e2 or e2 = e3: symmetry with 2-fold degeneracy
!
! The size of WORK should be 9 or larger.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine MIner2(Intact,NAtm,AMass,XYZCM,RI,E,WORK)
 implicit real(kind=8) (a-h,o-z)
 real(kind=8) :: AMass(NAtm),XYZCM(3,NAtm),RI(3,3),E(3),WORK(9)
 logical :: Intact

 call AClear(9,RI)
 do i=1,NAtm
   W=AMass(i)
   X=XYZCM(1,i)
   Y=XYZCM(2,i)
   Z=XYZCM(3,i)
   RI(1,1) = RI(1,1) + W * X * X
   RI(2,2) = RI(2,2) + W * Y * Y
   RI(3,3) = RI(3,3) + W * Z * Z
   RI(1,2) = RI(1,2) + W * X * Y
   RI(1,3) = RI(1,3) + W * X * Z
   RI(2,3) = RI(2,3) + W * Y * Z
 end do
 RI(2,1) = RI(1,2)
 RI(3,1) = RI(1,3)
 RI(3,2) = RI(2,3)

 call DiagS1(Intact,3,RI,E,WORK)

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! c = a - b, where a<0 or b<0 is a imaginary value
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
end

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
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! center of mass
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine MassCent(NAtm,AMass,XYZ,XYZCM,CM)
 implicit real(kind=8) (a-h,o-z)
 real(kind=8) :: AMass(*),XYZ(3,*),XYZCM(3,*),CM(3)

 call AClear(3,CM)
 Weight=0.d0
 do i=1,NAtm
   Weight=Weight+AMass(i)
   CM(1)=CM(1)+XYZ(1,i)*AMass(i)
   CM(2)=CM(2)+XYZ(2,i)*AMass(i)
   CM(3)=CM(3)+XYZ(3,i)*AMass(i)
 end do
 Weight=1.d0/Weight
 call AScale(3,Weight,CM,CM)

 do i=1,NAtm
   XYZCM(1,i)=XYZ(1,i)-CM(1)
   XYZCM(2,i)=XYZ(2,i)-CM(2)
   XYZCM(3,i)=XYZ(3,i)-CM(3)
 end do

 return
end

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
end

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
end

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
end

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
end

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
end

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
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! read an integral from a word
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function IfrmCha(N,Word)
 implicit real(kind=8) (a-h,o-z)
 character*(N) :: Word,CTmp

 CTmp = Word
 do i=1,N
   j=ichar(CTmp(i:i))
   if((j < 48) .or. (j > 57)) CTmp(i:i)=" "
 end do
 read(CTmp,*,Err=100,End=100)IfrmCha
 return

 100   IfrmCha = 0

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! sum shell; if A is not a one-dimensional array, the intrinsic function Sum doesn't work.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ASum(A,N)
 implicit real(kind=8) (a-h,o-z)
 real(kind=8) :: A(N)

 ASum = Sum(A,N)

 return
end

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
end

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
end

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
end

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
end

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
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! work space of N integers
! Since -cpp cannot be recognized by very old gfortran compilers, integer(kind=8) is always assumed.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function intwsp(N)
 implicit real(kind=8) (a-h,o-z)

! compile with "-cpp -D_I8_" (ifort -i8) or "-cpp"
!#ifdef _I8_
! intwsp = N
!#else
! intwsp = (N + 1) / 2
!#endif
 intwsp = N

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! work space of N characters
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ichawsp(N)
 implicit real(kind=8) (a-h,o-z)

 ichawsp = (N - 1) / 8 + 1

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! counting irreps
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine CountIrrep(Intact,iport,NAtm3,NVib,NClass,Irreps,ModMap)
 implicit real(kind=8) (a-h,o-z)
 logical :: Intact
 character*4 :: Irreps(NAtm3)
 dimension :: ModMap(NVib)

 rewind(iport)
 do i = 1, NAtm3
   read(iport,"(a4)") Irreps(i)
   if(i <= NVib .and. Irreps(i) == "****") call XError(Intact,"Unknown irreps of model system.")
 end do

 NClass = 0
 ModMap = 0
 do i = 1, NVib
   if(ModMap(i) == 0) then
     NClass = NClass + 1
     ModMap(i) = NClass
     ! save irreps of each class
     Irreps(NClass) = Irreps(i)
     do j = i+1, NVib
       if(Irreps(i) == Irreps(j)) ModMap(j) = NClass
     end do
   end if
 end do

 return
end

!--- END
