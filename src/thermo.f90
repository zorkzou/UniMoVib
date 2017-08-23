!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Thermochemistry calculation. Input:
! NAtom   : #atoms
! NAtm3   : NAtom*3
! NVib    : #vib. modes
! PGNAME  : point group symmetry symbol (without and with masses)
! AMass   : (array) atomic masses (a.m.u)
! XYZ     : (array) Cartesian coordinates (a.u.)
! Freq    : (array) vib. frequencies (a.u.) in Freq(1:NVib,IFrq),
!            where IFrq depends on IExpt
! Sc1,Sc2 : scratch. Size = NAtm3*NAtm3
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Thermochem(iinp,iout,Intact,NAtom,NAtm3,NVib,IFAtom,IExpt,PGNAME,AMass,XYZ,Freq,Sc1,Sc2)
implicit real(kind=8) (a-h,o-z)
parameter(tolr=0.01d0)
parameter( pi       = acos(-1.0d0),                                  &
!          the parameters are taken from Gaussian09 manual
           au2ang   = 0.52917720859d0,                               &   ! a.u. (length) to Angstrom
           clight   = 2.99792458d+08,                                &   ! m/s
           planck   = 6.62606896d-34,                                &   ! Js
           avogadro = 6.02214179d+23,                                &   ! particle/mol
           boltzman = 1.38065040d-23,                                &   ! J/K
           amu2kg   = 1.660538782d-27,                               &   ! amu to kilogram
           Atm2Pa   = 1.01325d+05,                                   &   ! pressure unit: 1 atm = 101325 pa
           Cal2Jou  = 4.184d0,                                       &   ! Cal to Joule
           H2Jou    = 0.435974394d-17,                               &   ! Hartree to Joule
!          derived factors
           b2met    = au2ang*1.0d-10,                                &   ! Bohr to metre
           H2kJ     = H2Jou*avogadro*1.0d-03,                        &   ! Hartree to kJ/mol
           H2kcal   = H2kJ/Cal2Jou,                                  &   ! Hartree to kcal/mol
           Rval     = avogadro*boltzman/(1.0d3*H2kJ),                &   ! R =8.31447247 J/(K*Mol) = 3.16681476d-6 au/K
           au2wn    = sqrt(avogadro*H2Jou/1.0d1)/(2*pi*clight*b2met) )   ! au (freq) to wavenumber (=5140.4871)

real(kind=8) :: AMass(NAtom),XYZ(3,NAtom),Freq(NAtm3,*),Sc1(3,*),Sc2(*)
character*4 :: PGNAME(2),PG
character*1 :: L2U
real(kind=8) :: energy(3),entropy(4)
logical :: Intact,IFAtom,IFLin

! Eel     : Total electronic energy from Q.C. calculation (a.u.)
! temp    : temperature (K)
! press   : pressure (atm)
! scale   : frequency scale factor
! PG =1   : use the point group without mass
!     2   : use the point group with mass (default)
!     xxxx: specify the name of point group, for example, D10h
namelist/Thermo/Eel,NDeg,temp,press,scale,PG

write(iout,"(//,1x,45('*'),/, ' ***   Thermal Contributions to Energies   ***',/, 1x,45('*'))")

Eel=0.d0
NDeg=1
temp=298.15d0
press=1.d0
scale=1.d0
PG="2    "

! to rot temperatures
cf1 = planck * planck / (boltzman*8.d0*Pi*Pi*b2met*b2met*amu2kg)
! rot temperatures to rot constants (GHZ)
cf2 = 1.0d-9 * boltzman / planck
! rot temperatures to rot constants (CM^-1)
cf3 = boltzman / (planck * 1.0d2 * clight)
! temperature to Hartree
cf4 = cf3/au2wn

rewind(iinp)
read(iinp,Thermo,err=2000,end=100)

100   temp=abs(temp)
press=abs(press)
! mass of molecule
VMas = ASum(AMass,NAtom)

if(IFAtom) then

  write(iout,"(/,                              &
  ' Temperature          :',f13.6,' Kelvin',/, &
  ' Pressure             :',f13.6,' Atm',/,    &
  ' Atomic mass          :',f13.6,' AMU',/,    &
  ' Elec. total energy   :',f13.6,' Hartree' )")temp,press,VMas,Eel

else
  IFrq = 3
  scale=abs(scale)
  ! expt. freq.
  if(IExpt .eq. 1) IFrq = 5

  ! PG = 1 or 2 (default)?
  IPG = 0
  read(PG,*,Err=300)IPG
  if(IPG < 1 .or. IPG > 2) IPG = 2
  PG = PGNAME(IPG)
  goto 310
300     IPG = 0    ! PG is symmetry symbol
310     continue

  if(PG == "****")PG = "C1  "
  PG(1:1)=L2U(PG(1:1))
  call charu2l(PG(2:))
  ! delete initial spaces
  PG = adjustl(PG)

! Rot. symmetry number
  call NRotSym(iout,PG,NSigma,IFLin)

  write(iout,"(/, &
  ' Temperature               :',4x,f13.6,4x,'Kelvin',/,   &
  ' Pressure                  :',4x,f13.6,4x,'Atm',/,      &
  ' Molecular mass            :',4x,f13.6,4x,'AMU',/,      &
  ' Elec. total energy        :',4x,f13.6,4x,'Hartree',/,  &
  ' Scale factor of Frequency :',4x,f13.6,/,             &
  ' Rot. symmetry number      :',4x,i6,/,                &
  ' The ',a4,' point group is used to calculate rotational entropy.' )") temp,press,VMas,Eel,scale,NSigma,PG

! calculate principal axes and moments of inertia --> Sc1(:,1:4)
  call RotCons(Intact,NAtom,AMass,XYZ,Sc1,Sc2)
  call RmNoise(12,1.0d-8,Sc1)
  write(iout,"(/, ' Principal axes and moments of inertia in atomic units:',/,37x,'1',19x,'2',19x,'3',/, &
    5x,'Eigenvalues --  ',4x,3f20.6)")Sc1(1:3,1)
  write(iout,"(11x,'X',13x,3f20.6)")(Sc1(1,i),i=2,4)
  write(iout,"(11x,'Y',13x,3f20.6)")(Sc1(2,i),i=2,4)
  write(iout,"(11x,'Z',13x,3f20.6)")(Sc1(3,i),i=2,4)

! calculate rot constants
  do i=1,3
    if(Sc1(i,1) > tolr) then
      Sc1(i,1) = cf1/Sc1(i,1)
    else
      Sc1(i,1) = 0.d0
    end if
  end do
  if(IFLin)then
    write(iout,"(/,' Rotational temperature',2x,f20.6,'    Kelvin')") Sc1(3,1)
    write(iout,"(' Rotational constant',5x,f20.6,'    cm**-1')") Sc1(3,1)*cf3
    write(iout,"(25x,f20.6,'    GHz')")Sc1(3,1)*cf2
  else
    write(iout,"(/,' Rotational temperatures',1x,3f20.6,'    Kelvin')") Sc1(1:3,1)
    write(iout,"(' Rot. constants A, B, C',2x,3f20.6,'    cm**-1')")Sc1(1:3,1)*cf3
    write(iout,"(25x,3f20.6,'    GHz')")Sc1(1:3,1)*cf2
  end if

end if

! translation
energy(1)=1.5d0*Rval*temp
entropy(1)=boltzman*temp
entropy(1)=(entropy(1)/(Atm2Pa*max(press,1.d-5))) * (2.d0*Pi*1.0d-3*VMas*entropy(1)/(avogadro*planck*planck))**1.5d0
entropy(1)=Rval*(log(entropy(1))+5.d0/2.d0)

! rotation
if(IFAtom)then
  energy(2)=0.d0
  entropy(2)=0.d0
else if(IFLin)then
  energy(2)=Rval*temp
  entropy(2)=( temp/max(Sc1(3,1),1.0d-6) )/dble(NSigma)
  entropy(2)=Rval*(1.0d0 + log(entropy(2)) )
else
  energy(2)=Rval*temp*1.5d0
  entropy(2)=Pi*temp*temp*temp/AMultip(Sc1,3,.True.,tolr/cf1)
  entropy(2)=Rval*(1.5d0 + log(sqrt(entropy(2))/dble(NSigma)) )
end if

! vibration
energy(3)=0.d0
entropy(3)=0.d0
do i=1,NVib
  ! at low temperature, no contributions from vibration
  if(Temp <= 10.d0) exit

  scale0=scale
  ! expt. freq should not be scaled
  if(Freq(i,6) > 0.0d0) scale0=1.0d0
  VT = Freq(i,IFrq)*scale0/cf4
  VTT= VT/temp
  ! neglect small freq. because it leads to big errors
  if(VT <= 1.d0) cycle

  energy(3)=energy(3) + Rval * VT / (exp(VTT) - 1.d0)
  entropy(3)=entropy(3) + Rval * ( VTT/(exp(VTT) - 1.d0) - log(1.d0 - exp(-VTT)) )
end do

! electronic contribution
!energy(4) = Eel
entropy(4)=Rval*log(dble(NDeg))

! ZPE
zpe=0.d0
do i=1,NVib
  ! neglect imag. freq.
  if(Freq(i,IFrq) <= 0.d0) cycle

  scale0=scale
  ! expt. freq should not be scaled
  if(Freq(i,6) > 0.0d0) scale0=1.0d0
  zpe = zpe + Freq(i,IFrq)*scale0
end do
zpe=zpe*0.5d0*Rval/cf4

! print results
eu=zpe+ASum(energy,3)

eh=eu+Rval*temp
eg=eh-temp*ASum(entropy,4)
write(iout,"(/,                                                 &
' Thermal correction energies',30x,'Hartree',12x,'kcal/mol',/,  &
' Zero-point Energy                          :',2f20.6,/,       &
' Thermal correction to Energy               :',2f20.6,/,       &
' Thermal correction to Enthalpy             :',2f20.6,/,       & 
' Thermal correction to Gibbs Free Energy    :',2f20.6)") zpe,zpe*H2kcal,eu,eu*H2kcal,eh,eh*H2kcal,eg,eg*H2kcal
if(abs(Eel) > tolr)write(iout,"(/,                        &
' Sum of electronic and zero-point Energies  :',f20.6,/,  &
' Sum of electronic and thermal Energies     :',f20.6,/,  &
' Sum of electronic and thermal Enthalpies   :',f20.6,/,  &
' Sum of electronic and thermal Free Energies:',f20.6)")zpe+Eel,eu+Eel,eh+Eel,eg+Eel

return
2000  call XError(Intact,"Please check your input of $Thermo!")
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! calculate rot constants in cm-1
!
! Eigenvalues and Eigenvectors are saved in Sc1.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RotCons(Intact,NAtom,AMass,XYZ,Sc1,Sc2)
implicit real(kind=8) (a-h,o-z)
real(kind=8) :: AMass(NAtom),XYZ(3,NAtom),Sc1(*),Sc2(*)
logical :: Intact

! center of mass ---> Sc2(1:3); XYZ in CMCS --> Sc2
call MassCent(NAtom,AMass,XYZ,Sc2,Sc1)

! principal moment of inertia;
! Only 3+9+9 elements in of Sc1 will be used.
! If NAtom > 1, the size of Sc1 is large enough
IE = 1          ! Eigenvalues
IR = IE + 3     ! Eigenvectors
IW = IR + 9
call MIner(Intact,NAtom,AMass,Sc2,Sc1(IR),Sc1(IE),Sc1(IW))

return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! calculate rot. symmetry number
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine NRotSym(iout,PG,NSigma,IFLin)
implicit real(kind=8) (a-h,o-z)
character*4 :: PG,LinLib(6)
logical :: IFLin
save LinLib
data LinLib/'Dooh','D*h','Dxh','Coov','C*v','Cxv'/

IFLin = .false.
Do i=1,6
  if(index(PG,LinLib(i)) > 0) then
    IFLin = .True.
    exit
  end if
end do

! Ci, Cs
if(PG(1:2) == "Ci" .or. PG(1:2) == "Cs")then
  NSigma = 1
! Dooh
else if(IFLin .and. PG(1:1) == "D")then
  NSigma = 2
! Coov
else if(IFLin .and. PG(1:1) == "C")then
  NSigma = 1
! T, Td, Th
else if(PG(1:1) == "T")then
  NSigma = 12
! O, Oh
else if(PG(1:1) == "O")then
  NSigma = 24
! I, Ih
else if(PG(1:1) == "I")then
  NSigma = 60
! Sn (n=2m)
else if(PG(1:1) == "S")then
  NSigma = IfrmCha(4,PG)/2
! Cn, Cnv, Cnh
else if(PG(1:1) == "C")then
  NSigma = IfrmCha(4,PG)
! Dn, Dnv, Dnh
else if(PG(1:1) == "D")then
  NSigma = IfrmCha(4,PG)*2
! Unknown
else
  write(iout,"(/, ' Unknown point group ',a4,'!',/,  &
    ' SIGMA = 1 will be assumed, but it may lead to significant errors in entropy',/, ' and free energy.')")PG
  NSigma = 1
end if

if(NSigma < 1)then
  write(iout,"(/, ' Unreasonable NSigma is detected: ',i5,/,  &
    ' SIGMA = 1 will be used instead, but it may lead to significant errors in',/,  &
    ' entropy and free energy.')")NSigma
  NSigma = 1
end if

return
end

!--- END

