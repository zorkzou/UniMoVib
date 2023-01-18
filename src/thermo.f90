!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Thermochemistry calculation. Input:
! ifbdfchk: print check data for bdf
! NAtom   : #atoms
! NAtm3   : NAtom*3
! NVib    : #vib. modes
! PGNAME  : point group symmetry symbol (without and with masses)
! AMass   : (array) atomic masses (a.m.u)
! XYZ     : (array) Cartesian coordinates (a.u.)
! Freq    : (array) vib. frequencies (a.u.) in Freq(1:NVib,IFrq), where IFrq =3 (IExpt /= 1) or 5 (IExpt == 1)
! Sc1,Sc2 : scratch. Size = NAtm3*NAtm3
! ctmp    : scratch for characters.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Thermochem(iinp,iout,Intact,ifbdfchk,NAtom,NAtm3,NVib,IFAtom,IExpt,PGNAME,AMass,XYZ,Freq,Sc1,Sc2,ctmp)
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
character*200 :: ctmp
character*4 :: PGNAME(2),PG
character*1 :: L2U
real(kind=8) :: energy(3),entropy(4)
logical :: Intact,ifbdfchk,IFAtom,IFLin,IFRdT,IfRdP

! Eel     : Total electronic energy from Q.C. calculation (a.u.)
! temp    : temperature (K)
! press   : pressure (atm)
! scale   : frequency scale factor
! PG =1   : use the point group without mass
!     2   : use the point group with mass (default)
!     xxxx: specify the name of point group, for example, D10h
namelist/Thermo/Eel,NDeg,temp,press,scale,sctol,PG
allocatable  :: freqtmp(:)

write(iout,"(//,1x,45('*'),/, ' ***   Thermal Contributions to Energies   ***',/, 1x,45('*'))")

Eel=0.d0
NDeg=1
temp=298.15d0
press=1.d0
scale=1.d0
sctol=0.d0
IFRdT=.false.
IfRdP=.false.
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

100   continue
if(temp < 0.d0)then
  IFRdT=.true.
  temp = 298.15d0
end if
if(press < 0.d0)then
  IfRdP=.true.
  press=1.d0
end if
! mass of molecule
VMas = ASum(AMass,NAtom)

if(IFAtom) then

  write(iout,"(/, &
  ' Atomic mass               :',4x,f13.6,4x,'AMU',/,    &
  ' Electronic total energy   :',4x,f13.6,4x,'Hartree' )")VMas,Eel

else
  scale=abs(scale)

  allocate(freqtmp(NVib))   ! scaled frequencies
  call ScaleF(NVib,NAtm3,IExpt,scale,sctol,Freq,freqtmp)

  ! PG = 1 or 2 (default)?
  IPG = 0
  read(PG,*,Err=300)IPG
  if(IPG < 1 .or. IPG > 2) IPG = 2
  PG = PGNAME(IPG)
  goto 310
  300  IPG = 0    ! PG is symmetry symbol
  310  continue

  if(PG == "****")PG = "C1  "
  PG(1:1)=L2U(PG(1:1))
  call charu2l(PG(2:))
  ! delete initial spaces
  PG = adjustl(PG)

! Rot. symmetry number
  call NRotSym(iout,PG,NSigma,IFLin)

  write(iout,"(/, &
  ' Molecular mass            :',4x,f13.6,4x,'AMU',/,    &
  ' Electronic total energy   :',4x,f13.6,4x,'Hartree',/,&
  ' Scaling factor of Freq.   :',4x,f13.6,/,             &
  ' Tolerance of scaling      :',4x,f13.6,4x,'cm^-1',/,  &
  ' Rotational symmetry number:',4x,i6,/,                &
  ' The ',a4,' point group is used to calculate rotational entropy.' )") VMas,Eel,scale,sctol,NSigma,PG

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
    write(iout,"(' Rotational constant',5x,f20.6,'    cm^-1')") Sc1(3,1)*cf3
    write(iout,"(25x,f20.6,'    GHz')")Sc1(3,1)*cf2
  else
    write(iout,"(/,' Rotational temperatures',1x,3f20.6,'    Kelvin')") Sc1(1:3,1)
    write(iout,"(' Rot. constants A, B, C',2x,3f20.6,'    cm^-1')")Sc1(1:3,1)*cf3
    write(iout,"(25x,3f20.6,'    GHz')")Sc1(1:3,1)*cf2
  end if

end if

! loop of Temperature & Pressure
IRd = 0
do while(.true.)

  if(IRd > 0) then
    if(.NOT. IFRdT .and. .NOT. IFRdP)then
      exit
    else
      read(iinp,"(a200)",end=1000,err=1000)ctmp
      if(len_trim(ctmp) == 0) goto 1000
      if(IFRdT .and. IFRdP)then
        read(ctmp,*,end=1000,err=1000) temp,press
      else if(IFRdT)then
      	read(ctmp,*,end=1000,err=1000) temp
      else if(IFRdP)then
        read(ctmp,*,end=1000,err=1000) press
      end if
    end if
    if(temp < 0.d0) temp = 298.15d0
    if(press < 0.d0) press=1.d0
  end if
  IRd = IRd + 1

  write(iout,"(//,' #',i4,4x,'Temperature = ',f15.5,' Kelvin',9x,'Pressure = ',f15.5,' Atm',/, &
    1x,84('=') )") IRd,temp,press

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

    VT = freqtmp(i)/cf4
    VTT= VT/temp
    ! neglect small freq. because it leads to big errors
    if(VT <= 1.d0) cycle

    energy(3)=energy(3) + Rval * VT / (exp(VTT) - 1.d0)
    entropy(3)=entropy(3) + Rval * ( VTT/(exp(VTT) - 1.d0) - log(1.d0 - exp(-VTT)) )
  end do

  ! electronic contribution
  !energy(4) = Eel
  entropy(4)=Rval*log(dble(NDeg))
  !write(iout,"(4f22.6)")entropy*temp*H2kcal

  ! ZPE
  zpe=0.d0
  do i=1,NVib
    ! neglect imag. freq.
    if(freqtmp(i) <= 0.d0) cycle

    zpe = zpe + freqtmp(i)
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
  write(iout,"(1x,84('=') )")

 ! print check data for bdf
 if(ifbdfchk) then
   call bdfchk(iout,.false.)
   write(iout,"('  CHECKDATA:BDFOPT:THERMO:  ',4f16.2)") zpe*H2kcal,eu*H2kcal,eh*H2kcal,eg*H2kcal
   call bdfchk(iout,.true.)
 end if

end do

1000 continue
if(.not. IFAtom) deallocate(freqtmp)

return
2000  call XError(Intact,"Please check your input of $Thermo!")
end subroutine Thermochem

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

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Put vibrational theoretical/experimental frequencies into freqtmp. The real theoretical ones may be scaled.
!
! Freq(:,3) : theoretical frequencies
! Freq(:,5) : theoretical (and experimentally corrected if Freq(:,6) = 1.0) frequencies if IExpt = 1
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine ScaleF(NVib,NAtm3,IExpt,fscale,sctol,Freq,freqtmp)
 Implicit Real*8(A-H,O-Z)
 dimension Freq(NAtm3,*),freqtmp(*)

 if(IExpt == 1) then
   do ivib=1,NVib
     freqtmp(ivib) = Freq(ivib,5)
     if(Freq(ivib,6) < 0.0d0 .and. freqtmp(ivib) > max(sctol,0.0d0)) freqtmp(ivib) = freqtmp(ivib) * fscale
   end do
 else
   do ivib=1,NVib
     freqtmp(ivib) = Freq(ivib,3)
     if(freqtmp(ivib) > max(sctol,0.0d0)) freqtmp(ivib) = freqtmp(ivib) * fscale
   end do
 end if

 return
end

!--- END

