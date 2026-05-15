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
end subroutine obt_idx

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
end subroutine PartMX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Mode = 0 : returns nuclear charge iza for an element symbol "el".
!     /= 0 : returns element symbol "el" for nuclear charge iza.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine ElemZA(Mode,el,iza)
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
   iza = 0
   do i=1,maxza
     if(index(el,atomlib(i)) /= 0)then
       iza = i
       exit
     end if
   end do

 else

   el = "???"
   if(iza > 0 .and. iza <= maxza) el = adjustl(atomlib(iza))
   call charu2l(el(2:3))

 end if

 return
end subroutine ElemZA

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
   84.911789738d0,  87.905612257d0,  88.905848300d0,  89.904704400d0,  92.906378100d0,  97.905404820d0,&
   98.906254700d0, 101.904349300d0, 102.905504000d0, 105.903486000d0, 106.905097000d0, 113.903358500d0,&
  114.903878000d0, 119.902194700d0, 120.903815700d0, 129.906224400d0, 126.904473000d0, 131.904153500d0,&
  132.905451933d0, 137.905247200d0, 138.906353300d0, 139.905438700d0, 140.907652800d0, 141.907723300d0,&
  144.912749000d0, 151.919732400d0, 152.921230300d0, 157.924103900d0, 158.925346800d0, 163.929174800d0,&
  164.930322100d0, 165.930293100d0, 168.934213300d0, 173.938862100d0, 174.940771800d0, 179.946550000d0,&
  180.947995800d0, 183.950931200d0, 186.955753100d0, 191.961480700d0, 192.962926400d0, 194.964791100d0,&
  196.966568700d0, 201.970643000d0, 204.974427500d0, 207.976652100d0, 208.980398700d0, 208.982430400d0,&
  209.987148000d0, 222.017577700d0, 223.019735900d0, 226.025409800d0, 227.027752100d0, 232.038055300d0,&
  231.035884000d0, 238.050788200d0, 237.048173400d0, 244.064204000d0, 243.061381100d0, 247.070354000d0,&
  247.070307000d0, 251.079587000d0, 252.082980000d0, 257.095106000d0, 258.098431000d0, 259.101030000d0,&
  266.119830000d0, 267.121790000d0, 268.125670000d0, 269.128630000d0, 270.133360000d0, 277.151900000d0,&
  278.156310000d0, 281.164510000d0, 283.170540000d0, 285.177120000d0, 286.182210000d0, 289.190420000d0,&
  290.195980000d0, 293.204490000d0, 294.210460000d0, 294.213920000d0, 295.000000000d0, 299.000000000d0/
 ! average isotopic masses
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
end function EleMas

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
end subroutine MasLib

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
end function AIRCrg

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
end subroutine planar

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
end subroutine rmnumb

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
end subroutine charu2l

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
end subroutine charl2u

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
end function L2U

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
end function U2L

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
end function nonspace

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
end subroutine XError

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
end subroutine estop

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
end subroutine Frq2FFX

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
end subroutine TRVec

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
end subroutine MIner

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
end subroutine MIner2

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
end subroutine MassCent

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
end function IfrmCha

!xxx!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!xxx!
!xxx! work space of N integers
!xxx! Since -cpp cannot be recognized by very old gfortran compilers, integer(kind=8) is always assumed.
!xxx!
!xxx!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!xxxfunction intwsp(N)
!xxx implicit real(kind=8) (a-h,o-z)
!xxx
!xxx! compile with "-cpp -D_I8_" (ifort -i8) or "-cpp"
!xxx!#ifdef _I8_
!xxx! intwsp = N
!xxx!#else
!xxx! intwsp = (N + 1) / 2
!xxx!#endif
!xxx intwsp = N
!xxx
!xxx return
!xxxend function intwsp

!xxx!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!xxx!
!xxx! work space of N characters
!xxx!
!xxx!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!xxxfunction ichawsp(N)
!xxx implicit real(kind=8) (a-h,o-z)
!xxx
!xxx ichawsp = (N - 1) / 8 + 1
!xxx
!xxx return
!xxxend function ichawsp

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
end subroutine CountIrrep

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! gradient information
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine GradInfo(iout,ifbdfchk,NAtm,NAtm3,NVib,amass,za,xyz,grd,ffx,al,Rslt,dx,scr,Elm)
 implicit real(kind=8) (a-h,o-z)
 parameter(eps=1.0d-10)
 logical :: ifbdfchk
 dimension :: amass(NAtm), za(NAtm), xyz(NAtm3), grd(3,NAtm), ffx(*), al(NAtm3,*), Rslt(NAtm3,*), dx(NAtm3), &
   scr(NAtm3,NAtm3)
 character*3 :: Elm
 allocatable :: avec(:,:), gvib(:), grslt(:)

 write(iout,"(//,1x,30('*'),/,' ***  Gradient Information  ***',/,1x,30('*'))")
 if(NAtm < 2) then
   write(iout,"(' N.A.')")
   return
 end if

 write(iout,"(/,' Cartesian gradients (a.u.)',/,1x,68('-'),/,3x, &
    'No.   Atom    ZA                 X             Y             Z',/,1x,68('-'))")
 do i=1,NAtm
   j = nint(za(i))
   call ElemZA(1,Elm,j)
   write(iout,"(i6,4x,a3,1x,i5,8x,3f14.8)") i,Elm,j,grd(:,i)
 end do
 write(iout,"(1x,68('-'))")

 allocate(avec(NAtm3,NVib), gvib(NVib), grslt(5))

 ! Restoring the purely vibrational eigenvectors (LL) of the vibrational secular equation.
 ! Translations and rotations have been projected out.
 !
 ! F * L = M * L * E  -->  FF * LL = LL * E,
 ! where FF = X * F * X, LL = Y * L, X = M^-1/2, and Y = M^1/2.
 !
 ! For the renormalized vibrational normal modes stored in the first NVib columns of AL,
 ! L = AL * XR
 ! where XR = MR^-1/2 and MR = AL' * M * AL.
 ! So, LL = Y * AL * XR, saved in avec.
 !
 ! Note: AL is affected by isotopes and therefore cannot be used for geometry optimization.
 do i = 1, NVib
   am = 0.0d0
   jm = 0
   do j = 1, NAtm3
     ja = (j-1)/3 + 1
     avec(j,i) = al(j,i) * sqrt(amass(ja) / Rslt(i,2))
     if(am < abs(avec(j,i))) then
       am = abs(avec(j,i))
       jm = j
     end if
   end do
   ! Flip the phase if the element with the largest absolute value (the first one among
   ! equal values) is negative.
   if(avec(jm,i) < 0.0d0) avec(:,i) = -avec(:,i)
 end do

 ! gradients in the vibrational space
 call MMpyMF(1,NAtm3,NVib,grd,avec,gvib)

 ! displacement required to find the minimum on the surface using Newton-Raphson step
 ! dX = F^-1 * g'.
 ! See Eq. (3) in JCP, 111, 10806 (1999), where H --> F and f --> g'.
 ! It is worse than the RFO (rational function optimization) step (see Eq. (9)) but is much simpler.
 !
 ! In the vibrational space, the above equation needs to be modified as follows
 ! dX = FF^-1 * G' = LL' * E^+1 * LL * LL' * g' = LL' * E^-1 * g'.
 ! (it needs to be derived in NAtm3-dimensional space, and then the rotational and translational modes
 ! can be discarded.)
 dx = 0.0d0
 do i = 1, NAtm3
   do j = 1, NVib
     eige = Rslt(j,1) / Rslt(j,2)
     if(abs(eige) < eps) eige = sign(eps,eige)
     dx(i) = dx(i) + avec(i,j) * gvib(j) / eige
   end do
 end do

 !write(iout,"(/,' Displacements (a.u.)',/,1x,68('-'),/,3x, &
 !   'No.   Atom    ZA                 X             Y             Z',/,1x,68('-'))")
 !do i=1,NAtm
 !  j = nint(za(i))
 !  call ElemZA(1,Elm,j)
 !  write(iout,"(i6,4x,a3,1x,i5,8x,3f14.8)") i,Elm,j,dx(3*i-2:3*i)
 !end do
 !write(iout,"(1x,68('-'))")

 ! Rotate XYZ_new (= XYZ + dx) to the orientation that maximizes overlap with XYZ,
 ! then recalculate displacements. This step may reduce the RMS displacement.
 scr(:,1) = xyz
 scr(:,2) = xyz + dx
 call rotmole(iout,NAtm,scr(1,1),scr(1,2),scr(1,3))
 dx = scr(:,2) - scr(:,1)

 !write(iout,"(/,' Rotated displacements (a.u.)',/,1x,68('-'),/,3x, &
 !   'No.   Atom    ZA                 X             Y             Z',/,1x,68('-'))")
 !do i=1,NAtm
 !  j = nint(za(i))
 !  call ElemZA(1,Elm,j)
 !  write(iout,"(i6,4x,a3,1x,i5,8x,3f14.8)") i,Elm,j,dx(3*i-2:3*i)
 !end do
 !write(iout,"(1x,68('-'))")

 grslt = 0.0d0

 do i = 1, NAtm3
   ia = (i-1)/3 + 1
   ix = i - 3*ia +3
   grslt(1) = max(grslt(1),abs(dx(i)))
   grslt(3) = max(grslt(3),abs(grd(ix,ia)))
 end do
 grslt(2) = dotx(NAtm3,dx,dx)
 grslt(2) = sqrt(grslt(2)/dble(NAtm3))
 grslt(4) = dotx(NAtm3,grd,grd)
 grslt(4) = sqrt(grslt(4)/dble(NAtm3))
 ! dE: see Eq. (1) in JCP, 111, 10806 (1999).
 call MMpyMF(1,NAtm3,NAtm3,dx,ffx,scr)
 grslt(5) = abs(0.5d0*dotx(NAtm3,scr,dx) + dotx(NAtm3,grd,dx))

 call prtconv(iout,grslt)

 ! print check data for bdf
 if(ifbdfchk) then
   call bdfchk(iout,.false.)
   write(iout,"('  CHECKDATA:BDFOPT:CONVERGE:',4f16.6)") grslt(1:4)
   write(iout,"('  CHECKDATA:BDFOPT:CONVERGE:',d16.2)") grslt(5)
   call bdfchk(iout,.true.)
 end if

 deallocate(avec, gvib, grslt)

 return
end subroutine GradInfo

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! gradient information
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!subroutine GradInfo_old(iout,ifbdfchk,NAtm,NAtm3,NVib,za,grd,ffx,al,fcon,dx,scr,Elm)
! implicit real(kind=8) (a-h,o-z)
! logical :: ifbdfchk
! dimension :: za(NAtm), grd(3,NAtm), ffx(*), al(NAtm3,*), fcon(NAtm3), dx(3,*), scr(NAtm3,NAtm3)
! character*3 :: Elm
! allocatable :: grslt(:)
!
! write(iout,"(//,1x,30('*'),/,' ***  Gradient Information  ***',/,1x,30('*'))")
! if(NAtm < 2) then
!   write(iout,"(' N.A.')")
!   return
! end if
!
! write(iout,"(/,' Cartesian gradients (a.u.)',/,1x,68('-'),/,3x, &
!    'No.   Atom    ZA                 X             Y             Z',/,1x,68('-'))")
! do i=1,NAtm
!   j = nint(za(i))
!   call ElemZA(1,Elm,j)
!   write(iout,"(i6,4x,a3,1x,i5,8x,3f14.8)") i,Elm,j,grd(:,i)
! end do
! write(iout,"(1x,68('-'))")
!
! ! Newton-Raphson step dX = -F^-1 * g,
! ! see Eq. (3) in JCP, 111, 10806 (1999), where H --> F and f --> -g.
! ! It is worse than the RFO (rational function optimization) step (see Eq. (9)) but is much simpler.
!
! ! F^-1 calculation:
! ! From F L = M L E, we have FF = LL E LL' where FF = X F X, LL = X^-1 L, and X = M^-1/2.
! ! Then F^-1 = L E^-1 L'.
! ! However AL saves renormalized L0 = L mu^-1/2, and E^-1 = mu k^-1, so F^-1 = L0 k^-1 L0'.
! do i = 1, NVib
!   x = sign(max(abs(fcon(i)), 1.0d-8), fcon(i))
!   scr(:,i) = (1.0d0/x) * al(:,i)
! end do
! call DGEMM('N','T',NAtm3,NAtm3,NVib,1.d0,scr,NAtm3,al,NAtm3,0.d0,dx,NAtm3)  ! F^-1 --> dx
! call MMpyMF(NAtm3,NAtm3,1,dx,grd,scr)
! call AScale(NAtm3,-1.0d0,scr,dx)
!
! allocate(grslt(5))
! grslt = 0.0d0
!
! do i = 1, NAtm
!   do j = 1, 3
!     grslt(1) = max(grslt(1),abs(dx(j,i)))
!     grslt(2) = grslt(2) + dx(j,i) * dx(j,i)
!     grslt(3) = max(grslt(3),abs(grd(j,i)))
!     grslt(4) = grslt(4) + grd(j,i) * grd(j,i)
!   end do
! end do
! grslt(2) = sqrt(grslt(2)/dble(NAtm3))
! grslt(4) = sqrt(grslt(4)/dble(NAtm3))
! ! dE: see Eq. (1) in JCP, 111, 10806 (1999), where f = -g.
! call MMpyMF(1,NAtm3,NAtm3,dx,ffx,scr)
! grslt(5) = abs(0.5d0*dotx(NAtm3,scr,dx) + dotx(NAtm3,grd,dx))
!
! call prtconv(iout,grslt)
!
! ! print check data for bdf
! if(ifbdfchk) then
!   call bdfchk(iout,.false.)
!   write(iout,"('  CHECKDATA:BDFOPT:CONVERGE:',4f16.6)") grslt(1:4)
!   write(iout,"('  CHECKDATA:BDFOPT:CONVERGE:',d16.2)") grslt(5)
!   call bdfchk(iout,.true.)
! end if
!
! deallocate(grslt)
!
! return
!end subroutine GradInfo_old

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! rotate the probe coordinates (in xyz1) with the best fit superimposition of the target coordinates (in xyz0).
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine rotmole(iout,natom,xyz0,xyz1,scr)
 implicit real(kind=8) (a-h,o-z)
 dimension :: xyz0(3,natom), xyz1(3,natom), scr(*)
 allocatable   :: rmat(:)

 ! calculate rotation matrix
 allocate(rmat(9))
 call qrotmole(0,iout,natom,xyz0,xyz1,rmat,scr)

 ! rotate the probe coordinates
 call rotvec(natom,-1,rmat,xyz1,scr)

 deallocate(rmat)

 return
end subroutine rotmole

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! This subroutine determines the rotation matrix for the best fit superimposition of two molecular Coordinates.
!
! The basic method was described in S. K. Kearsley, Acta Cryst. A45, 208 (1989), and coded in the PDBSUP program by B.Rupp and
! S.Parkin at LLNL (1996).
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine qrotmole(iprint,iout,natm,cartarget,carprobe,rmat,scr)
 implicit real(kind=8) (a-h,o-z)
 dimension cartarget(3,natm), carprobe(3,natm), rmat(3,3), scr(4,4)
 allocatable :: dxp(:,:), dxm(:,:), qmat(:,:), eige(:)

 allocate(dxp(3,natm), dxm(3,natm), qmat(4,4), eige(4))

 ! translate the molecules to their geometric centers
 call shift(natm,scr,cartarget)
 call shift(natm,scr,carprobe)
 if(iprint > 0) then
   write(iout,"(' Cartesian coordinates of target molecule:')")
   do i = 1, natm
     write(iout,"(3f20.12)") cartarget(:,i)
   end do
   write(iout,"(' Cartesian coordinates of probe molecule:')")
   do i = 1, natm
     write(iout,"(3f20.12)") carprobe(:,i)
   end do
 end if

 ! coordinate differences: plus and minus
 dxp = carprobe + cartarget
 dxm = carprobe - cartarget

 ! construct Kearsley's Q-matrix
 call conqmt(iprint,iout,natm,dxp,dxm,qmat)

 ! diagonalize Q: Q * L = L * A
 call DiagS1(.false.,4,qmat,eige,scr)
 ! Acta Cryst. (1989). A45, 208-210.
 ! See the bottom right corner of page 209 about r.m.s. deviation.
 if(iprint > 0) write(imol,"(/,' RMS deviation: ',f7.4)") sqrt(abs(eige(1))/dble(natm))

 ! construct the best fit rotation matrix using the eigenvectors
 call conrot(iprint,iout,qmat,rmat)

 deallocate(dxp, dxm, qmat, eige)

end subroutine qrotmole

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! print BDF CHECKDATA block.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine bdfchk(iout,ifend)
 implicit real(kind=8) (a-h,o-z)
 logical :: ifend

 if(.not. ifend) then
   write(iout,"(/,1x,39('+'),' DATA CHECK ',40('+'))")
 else
   write(iout,"(1x,37('+'),' END DATA CHECK ',38('+')),/")
 end if

 return
end subroutine bdfchk

!--- END
