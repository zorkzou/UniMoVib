!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Read NAtm
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RdNAtm1(idt0,idt1,Intact,IOP,NAtm,ctmp)
implicit real(kind=8) (a-h,o-z)
logical :: Intact
dimension :: IOP(*)
character*100 :: tag, ctmp

NAtm=0

select case(IOP(1))

  case(-1) ! atom
    NAtm=1
    return

  case(-2) ! UniMoVib (alm)
    call RdNAtmALM(idt0,NAtm)

  case(-3) ! xyz
    call RdNAtmXYZ(idt0,NAtm,ctmp)

  case(1)  ! Gaussian
    call RdNAtmGauss(idt0,NAtm,tag,ctmp)

  case(2,3)  ! Gamess, Firefly
    call RdNAtmGAMES(idt1,NAtm,tag,ctmp)

  case(4)  ! ORCA
    call RdNAtmORCA(idt0,NAtm,tag,ctmp)

  case(5)  ! CFour
    call RdNAtmCFour(idt1,NAtm,ctmp)

  case(6)  ! Molpro
    call RdNAtmMolp(idt0,NAtm,tag,ctmp)

  case(7)  ! QChem
    call RdNAtmQChem(idt0,NAtm,tag,ctmp)

  case(8)  ! NWChem
    call RdNAtmNWChem(idt0,NAtm,tag,ctmp)

  case(9)  ! GAMESS-UK
    call RdNAtmGMSUK(idt0,NAtm,tag,ctmp)

  case(10) ! Turbomole
    call RdNAtmTurbom(idt0,NAtm,tag,ctmp)

  case(11) ! deMon
    call RdNAtmDeMon(idt0,NAtm,tag,ctmp)

  case(12) ! PQS
    call RdNAtmPQS(idt0,NAtm,tag,ctmp)

  case(13) ! MOPAC
    call RdNAtmMOPAC(idt0,NAtm,tag,ctmp)

  case(14) ! AMPAC/AMSOL
    call RdNAtmAMPAC(idt0,NAtm,tag,ctmp)

  case(15) ! DALTON
    call RdNAtmDALTON(idt0,NAtm,tag,ctmp)

  case(16) ! FHI-AIMS
    call RdNAtmAIMS(idt0,NAtm,ctmp)

  case(17) ! CP2K
    call RdNAtmCP2K(idt0,NAtm,tag,ctmp)

  case(18) ! ADF
    call RdNAtmADF(idt0,NAtm,tag,ctmp)

  case(19) ! Hyperchem
    call RdNAtmHyper(idt0,NAtm,tag,ctmp)

  case(20) ! Jaguar
    call RdNAtmJag(idt0,NAtm,tag,ctmp)

  case(21) ! MOLDEN
    call RdNAtmMOLDEN(idt0,NAtm,tag,ctmp)

  case(22) ! Crystal
    call RdNAtmCry(idt0,Intact,NAtm,tag,ctmp)

  case(23) ! Spartan
    call RdNAtmSptn(idt0,Intact,NAtm,tag,ctmp)

  case(24) ! PSI
    call RdNAtmPSI(idt0,Intact,NAtm,tag,ctmp)

  case(25) ! DMol3
    call RdNAtmDMol(idt0,Intact,NAtm,tag,ctmp)

  case(26) ! ACES
    call RdNAtmACES(idt0,Intact,NAtm,tag,ctmp)

  case default
    call XError(Intact,"N.Y.I. in RdNAtm1!")

end select

if(Natm < 1) call XError(Intact,"NAtm < 1.")

if(Natm == 1) call XError(Intact,'Specify QCPROG="ATOMCALC" for atomic calculation.')

return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Read Q.C. data
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RdData1(iout,idt0,idt1,idt2,Intact,IOP,Infred,IRaman,NAtm,ctmp,AMass,ZA,XYZ,FFx,APT,DPol,Scr1,Scr2,Scr3,Scr4)
implicit real(kind=8) (a-h,o-z)
logical :: Intact
dimension :: IOP(*)
real(kind=8) :: AMass(*), ZA(*), XYZ(*), FFx(*), APT(*), DPol(*), Scr1(*), Scr2(*), Scr3(*), Scr4(*)
character*100 :: tag,ctmp

NAtm3=3*NAtm
Infred=0
IRaman=0

select case(IOP(1))

  case(-1) ! atom
    return

  case(-2) ! UniMoVib (ALM); the size of Scr4 should be 9*NAtm3 at least
    call RdALMode(idt0,iout,Intact,Infred,IRaman,NAtm,ctmp,AMass,ZA,XYZ,FFx,APT,DPol,Scr4)

  case(-3) ! xyz
    call RdXYZ(idt0,iout,Intact,NAtm,ctmp,ZA,XYZ,FFx,Scr1,Scr2,Scr3,Scr4)
!   the most abundant isotopic masses are used
    call MasLib(0,NAtm,AMass,ZA)

  case(1)  ! Gaussian
    call RdGauss(idt0,iout,tag,ctmp,Intact,IOP(4),Infred,IRaman,NAtm,AMass,ZA,XYZ,FFx,APT,DPol,Scr1)
!   the most abundant isotopic masses are assumed
    if(AMass(1) < 0.d0) call MasLib(0,NAtm,AMass,ZA)

  case(2,3)  ! Gamess,Firefly
    call RdGAMES(idt0,idt1,tag,ctmp,Intact,Infred,NAtm,AMass,ZA,XYZ,FFx,APT)

  case(4)  ! ORCA
    call RdORCA(idt0,ctmp,Intact,Infred,NAtm,AMass,ZA,XYZ,FFx,APT)

  case(5)  ! CFour
    call ChkAFRQ(idt0,tag,ctmp,Intact)
    call RdCFour(idt0,idt1,tag,ctmp,Intact,Infred,NAtm,AMass,ZA,XYZ,FFx,APT,Scr1)

  case(6)  ! Molpro
    call RdMolp(idt0,tag,ctmp,Intact,Infred,NAtm,AMass,ZA,XYZ,FFx,APT)

  case(7)  ! QChem
    call RdQChem(idt0,tag,ctmp,Intact,IOP(4),NAtm,AMass,ZA,XYZ,FFx,APT,Scr1)
!   the most abundant isotopic masses are assumed
    if(AMass(1) < 0.d0) call MasLib(0,NAtm,AMass,ZA)

  case(8)  ! NWChem
    call RdNWChem(idt0,idt1,idt2,tag,ctmp,Intact,Infred,NAtm,AMass,ZA,XYZ,FFx,APT,Scr1)

  case(9)  ! GAMESS-UK
    call RdGMUK(idt0,tag,ctmp,Intact,Infred,NAtm,AMass,ZA,XYZ,FFx,APT)

  case(10) ! Turbomole
    call RdTurbm(idt0,idt1,tag,ctmp,Intact,Infred,NAtm,AMass,ZA,XYZ,FFx,APT)

  case(11) ! deMon
    call RdDeMon(idt0,tag,ctmp,Intact,NAtm,AMass,ZA,XYZ,FFx)

  case(12) ! PQS
    call RdPQS(idt0,idt1,idt2,tag,ctmp,Intact,Infred,NAtm,AMass,ZA,XYZ,FFx,APT,Scr1)

  case(13) ! MOPAC
    call RdMopac(idt0,tag,ctmp,Intact,NAtm,ZA,XYZ,FFx)
!   the most abundant isotopic masses are assumed
!   call MasLib(0,NAtm,AMass,ZA)
!   the averaged isotopic masses (default in MOPAC)
    call MasLib(1,NAtm,AMass,ZA)

  case(14) ! AMPAC/AMSOL
    call RdAMPAC(idt0,tag,ctmp,Intact,NAtm,ZA,XYZ,FFx)
!   the most abundant isotopic masses are assumed
!   call MasLib(0,NAtm,AMass,ZA)
!   the averaged isotopic masses (default in AMPAC)
    call MasLib(1,NAtm,AMass,ZA)

  case(15) ! Dalton
    call RdDalton(idt0,tag,tag(51:),ctmp,Intact,Infred,NAtm,AMass,ZA,XYZ,FFx,APT)

  case(16) ! FHI-AIMS
    call RdAIMS(idt0,idt1,idt2,ctmp,Intact,Infred,NAtm,AMass,ZA,XYZ,FFx,APT)

  case(17) ! CP2K
    call RdCP2K(idt0,tag,ctmp,Intact,NAtm,AMass,ZA,XYZ,FFx)

  case(18) ! ADF
    NWK=2*max(NAtm3,2)*NAtm3
    call RdADF(idt0,tag,ctmp,Intact,Infred,NAtm,AMass,ZA,XYZ,FFx,APT,Scr1,Scr2,Scr3,NWK,Scr4)

  case(19) ! Hyperchem
    call RdHyper(idt0,tag,ctmp,Intact,NAtm,AMass,ZA,XYZ,FFx,Scr1,Scr2,Scr3,Scr4)

  case(20) ! Jaguar
    call RdJaguar(idt0,tag,ctmp,Intact,NAtm,AMass,ZA,XYZ,FFx,Scr1,Scr2,Scr3,Scr4)

  case(21) ! MOLDEN
    call RdMOLDEN(idt0,tag,ctmp,Intact,NAtm,AMass,ZA,XYZ,FFx,Scr1,Scr2,Scr3,Scr4)

  case(22) ! Crystal
    call RdCry(idt0,tag,ctmp,Intact,NAtm,AMass,ZA,XYZ,FFx,Scr1,Scr2,Scr3,Scr4)

  case(23) ! Spartan
    call RdSptn(idt0,tag,ctmp,Intact,Infred,NAtm,ZA,XYZ,FFx,APT,Scr1)
!   the most abundant isotopic masses are assumed
    call MasLib(0,NAtm,AMass,ZA)

  case(24) ! PSI
    call RdPSI(idt0,tag,ctmp,Intact,NAtm,AMass,ZA,XYZ,FFx)

  case(25) ! DMOL3
    call RdDMol(idt0,tag,ctmp,Intact,NAtm,AMass,ZA,XYZ,FFx,Scr1,Scr2,Scr3,Scr4)

  case(26) ! ACES
    call RdACES(idt0,tag,ctmp,Intact,NAtm,AMass,ZA,XYZ,FFx,Scr1,Scr2,Scr3,Scr4)

  case default
    call XError(Intact,"N.Y.I. in RdData1!")

end select

! symmetrize FFx
call Symtrz(NAtm3,FFx,FFx)

return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Read NAtm from UniMoVib (ALM) data file
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RdNAtmALM(ifchk,NAtm)
implicit real(kind=8) (a-h,o-z)
logical :: Intact

i=0

rewind(ifchk)
read(ifchk,"(/)",err=100,end=100)
read(ifchk,*,err=100,end=100) i
100   NAtm=i

return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Read NAtm from XYZ file
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RdNAtmXYZ(ifchk,NAtm,ctmp)
implicit real(kind=8) (a-h,o-z)
character*100 :: ctmp

i=0

rewind(ifchk)
read(ifchk,"(a100)",end=100)ctmp
read(ctmp,*,err=100)i
100   NAtm=i

return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Read NAtm from Gaussian fchk
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RdNAtmGauss(ifchk,NAtm,tag,ctmp)
implicit real(kind=8) (a-h,o-z)
character*100 :: ctmp
character*44 :: tag

tag='Number of atoms                            I'

rewind(ifchk)
10    read(ifchk,"(a100)",end=100)ctmp
if(index(ctmp,tag)==0)goto 10

read(ctmp(45:100),*)NAtm
100   continue

return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Read NAtm from *.out of Gamess or Firefly
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RdNAtmGAMES(iout,NAtm,tag,ctmp)
implicit real(kind=8) (a-h,o-z)
character*(*) :: ctmp
character*56 :: tag

i=0
tag='ATOM      ATOMIC                      COORDINATES (BOHR)'

rewind(iout)
10    read(iout,"(a100)",end=100)ctmp
if(index(ctmp,tag)==0)goto 10

read(iout,*)
do while(.true.)
  read(iout,"(a100)",end=100)ctmp
  if(len_trim(ctmp)==0) exit
  i=i+1
end do
100   NAtm=i

return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Read NAtm from *.hess of ORCA
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RdNAtmORCA(ifchk,NAtm,tag,ctmp)
implicit real(kind=8) (a-h,o-z)
character*(*) :: ctmp
character*6 :: tag

tag='$atoms'

rewind(ifchk)
10    read(ifchk,"(a6)",end=100)ctmp
if(index(ctmp,tag) == 0)goto 10

read(ifchk,*)NAtm
100   continue

return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Read NAtm from GRD of CFour
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RdNAtmCFour(igeom,NAtm,ctmp)
implicit real(kind=8) (a-h,o-z)
character*100 :: ctmp

i=0
rewind(igeom)
read(igeom,"(a100)",err=100,end=100)ctmp
if(len_trim(ctmp) == 0)goto 100

read(ctmp,*,err=100) i

100   NAtm=i

return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Read NAtm from *.out of Molpro
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RdNAtmMolp(ifchk,NAtm,tag,ctmp)
implicit real(kind=8) (a-h,o-z)
character*(*) :: ctmp
character*56 :: tag

i=0
tag=' PROGRAM * HESSIAN'

rewind(ifchk)
10    read(ifchk,"(a18)",end=100)ctmp
if(index(ctmp,tag)==0)goto 10

tag='  Nr  Atom  Charge       X              Y              Z'
20    read(ifchk,"(a56)",end=100)ctmp
if(index(ctmp,tag)==0)goto 20
read(ifchk,*)
do while(.true.)
  read(ifchk,"(a100)",end=100)ctmp
  if(len_trim(ctmp)==0)exit
  i=i+1
end do
100   NAtm=i

return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Read NAtm from QChem fchk
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RdNAtmQChem(ifchk,NAtm,tag,ctmp)
implicit real(kind=8) (a-h,o-z)
character*100 :: ctmp
character*44 :: tag

tag='Number of atoms                            I'

rewind(ifchk)
10    read(ifchk,"(a100)",end=100)ctmp
if(index(ctmp,tag)==0)goto 10

read(ctmp(45:100),*)NAtm
100   continue

return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Read NAtm from NWChem fchk
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RdNAtmNWChem(ifchk,NAtm,tag,ctmp)
implicit real(kind=8) (a-h,o-z)
character*100 :: ctmp
character*69 :: tag

i=0

tag='     atom    #        X              Y              Z            mass'
rewind(ifchk)
10    read(ifchk,"(a100)",end=100)ctmp
if(index(ctmp,tag)==0)goto 10

tag=' --------------------------------------------------------------------'
read(ifchk,*)
do while(.true.)
  read(ifchk,"(a100)",end=100)ctmp
  if(index(ctmp,tag) /= 0)exit
  i=i+1
end do
100   NAtm=i

return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Read NAtm from GAMESS-UK's *.out
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RdNAtmGMSUK(ifchk,NAtm,tag,ctmp)
implicit real(kind=8) (a-h,o-z)
character*100 :: ctmp
character*99 :: tag

i=0

tag='     atom              x                   y                   z             nucleus  atomic weight'
rewind(ifchk)
10    read(ifchk,"(a100)",end=100)ctmp
if(index(ctmp,tag)==0)goto 10

read(ifchk,*,end=100)
do while(.true.)
  read(ifchk,"(a100)",end=100)ctmp
  if(len_trim(ctmp) == 0)exit
  i=i+1
end do
100   NAtm=i

return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Read NAtm from Turbomole's fchk (aoforce.out)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RdNAtmTurbom(ifchk,NAtm,tag,ctmp)
implicit real(kind=8) (a-h,o-z)
character*100 :: ctmp
character*72 :: tag

i=0

tag='                    atomic coordinates            atom    charge  isotop'
rewind(ifchk)
10    read(ifchk,"(a100)",end=100)ctmp
if(index(ctmp,tag)==0)goto 10

do while(.true.)
  read(ifchk,"(a100)",end=100)ctmp
  if(len_trim(ctmp) == 0)exit
  i=i+1
end do
100   NAtm=i

return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Read NAtm from deMon's *.out
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RdNAtmDeMon(ifchk,NAtm,tag,ctmp)
implicit real(kind=8) (a-h,o-z)
character*100 :: ctmp
character*17 :: tag

i=0

tag=' NUMBER OF ATOMS:'
rewind(ifchk)
10    read(ifchk,"(a100)",end=100)ctmp
if(index(ctmp,tag)==0)goto 10

read(ctmp(18:),*)NAtm

100   return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Read NAtm from PQS fchk
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RdNAtmPQS(ifchk,NAtm,tag,ctmp)
implicit real(kind=8) (a-h,o-z)
character*100 :: ctmp
character*12 :: tag

i=0

tag='$coordinates'
rewind(ifchk)
10    read(ifchk,"(a100)",end=100)ctmp
if(index(ctmp,tag)==0)goto 10

tag='$end'
do while(.true.)
  read(ifchk,"(a100)",end=100)ctmp
  if(index(ctmp,tag)/=0)exit
  i=i+1
end do
100   NAtm=i

return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Read NAtm from MOPAC output
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RdNAtmMOPAC(ifchk,NAtm,tag,ctmp)
implicit real(kind=8) (a-h,o-z)
character*100 :: ctmp
character*44 :: tag

i=0

tag='ORIENTATION OF MOLECULE IN FORCE CALCULATION'
rewind(ifchk)
10    read(ifchk,"(a100)",end=100)ctmp
if(index(ctmp,tag) == 0)goto 10

read(ifchk,"(//)",end=100)
do while(.true.)
  read(ifchk,"(a100)",end=100)ctmp
  if(len_trim(ctmp) == 0)exit
  i=i+1
end do
100   NAtm=i

return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Read NAtm from AMPAC/AMSOL output
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RdNAtmAMPAC(ifchk,NAtm,tag,ctmp)
implicit real(kind=8) (a-h,o-z)
character*100 :: ctmp
character*44 :: tag

i=0

tag='ORIENTATION OF MOLECULE IN FORCE CALCULATION'
rewind(ifchk)
10    read(ifchk,"(a100)",end=100)ctmp
if(index(ctmp,tag) == 0)goto 10

read(ifchk,"(//)",end=100)
do while(.true.)
  read(ifchk,"(a100)",end=100)ctmp
  if(len_trim(ctmp) == 0)exit
  i=i+1
end do
100   NAtm=i

return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Read NAtm from DALTON output
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RdNAtmDALTON(ifchk,NAtm,tag,ctmp)
implicit real(kind=8) (a-h,o-z)
character*100 :: ctmp
character*23 :: tag

i=0

tag='Molecular geometry (au)'
rewind(ifchk)
10    read(ifchk,"(a100)",end=100)ctmp
if(index(ctmp,tag) == 0)goto 10

read(ifchk,"(/)",end=100)
do while(.true.)
  read(ifchk,"(a100)",end=100)ctmp
  if(len_trim(ctmp) == 0)exit
  i=i+1
end do
100   NAtm=i

return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Read NAtm from FHI-AIMS output
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RdNAtmAIMS(ifchk,NAtm,ctmp)
implicit real(kind=8) (a-h,o-z)
character*100 :: ctmp

i=0

rewind(ifchk)
do while(.true.)
  read(ifchk,"(a100)",end=100)ctmp
  if(len_trim(ctmp) == 0)exit
  i=i+1
end do
100   NAtm=i

return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Read NAtm from CP2K output
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RdNAtmCP2K(ifchk,NAtm,tag,ctmp)
implicit real(kind=8) (a-h,o-z)
character*100 :: ctmp
character*42 :: tag

i=0

tag='Total number of            - Atomic kinds:'
rewind(ifchk)
10    read(ifchk,"(a100)",end=100)ctmp
if(index(ctmp,tag) == 0)goto 10

tag='                             - Atoms:'
read(ifchk,"(a100)",end=100)ctmp
if(index(ctmp,tag) == 0)goto 100
read(ctmp(38:),*,end=100)i
100   NAtm=i

return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Read NAtm from ADF tape21 or tape13
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RdNAtmADF(ifchk,NAtm,tag,ctmp)
implicit real(kind=8) (a-h,o-z)
character*100 :: ctmp
character*17 :: tag(2)

i=0
! #atoms (including dummy atoms)
tag(1)='Geometry'
tag(2)='nr of atoms'
rewind(ifchk)
do while(.true.)
  read(ifchk,"(a100)",end=100,err=100)ctmp
  if(ctmp(1:8) == tag(1)(1:8))then
    read(ifchk,"(a100)",end=100,err=100)ctmp
    if(ctmp(1:11) == tag(2)(1:11))then
      read(ifchk,"(/,a100)",end=100,err=100)ctmp
      read(ctmp,*,end=100,err=100)i
      goto 100
    end if
  end if
end do
100   NAtm=i

id=0
! #dummy atoms
tag(2)='nr of dummy atoms'
rewind(ifchk)
do while(.true.)
  read(ifchk,"(a100)",end=200,err=200)ctmp
  if(ctmp(1:8) == tag(1)(1:8))then
    read(ifchk,"(a100)",end=200,err=200)ctmp
    if(ctmp(1:17) == tag(2)(1:17))then
      read(ifchk,"(/,a100)",end=200,err=200)ctmp
      read(ctmp,*,end=200,err=200)id
      goto 200
    end if
  end if
end do
200   NAtm=NAtm-id

return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Read NAtm from Hyperchem log file
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RdNAtmHyper(ifchk,NAtm,tag,ctmp)
implicit real(kind=8) (a-h,o-z)
character*100 :: ctmp
character*69 :: tag

i=0

tag='Atom  Z     Charge            Coordinates(Angstrom)              Mass'
rewind(ifchk)
10    read(ifchk,"(a100)",end=100)ctmp
if(index(ctmp,tag) == 0)goto 10

read(ifchk,*,end=100)
do while(.true.)
  read(ifchk,"(a100)",end=100)ctmp
  if(len_trim(ctmp) == 0) exit
  i = i + 1
end do
100   NAtm=i

return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Read NAtm from Jaguar output
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RdNAtmJag(ifchk,NAtm,tag,ctmp)
implicit real(kind=8) (a-h,o-z)
character*100 :: ctmp
character*15 :: tag

i=0

tag='final geometry:'
rewind(ifchk)
do while(.true.)
  read(ifchk,"(a100)",end=20)ctmp
  if(index(ctmp,tag) /= 0) goto 40
end do

20    tag='Input geometry:'
rewind(ifchk)
do while(.true.)
  read(ifchk,"(a100)",end=100)ctmp
  if(index(ctmp,tag) /= 0) goto 40
end do

40    read(ifchk,"(/)",end=100)
do while(.true.)
  read(ifchk,"(a100)",end=100)ctmp
  if(len_trim(ctmp) == 0) exit
  i = i + 1
end do
100   NAtm=i

return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Read NAtm from MOLDEN file
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RdNAtmMOLDEN(ifchk,NAtm,tag,ctmp)
implicit real(kind=8) (a-h,o-z)
character*100 :: ctmp
character*10 :: tag

i=0

tag='[FR-COORD]'
rewind(ifchk)
10    read(ifchk,"(a100)",end=100)ctmp
call charl2u(ctmp)
if(index(ctmp,tag) == 0)goto 10

do while(.true.)
  read(ifchk,"(a100)",end=100)ctmp
  if(len_trim(ctmp) == 0) exit
  if(index(ctmp,"[") /= 0) exit
  if(index(ctmp,"]") /= 0) exit
  i = i + 1
end do
100   NAtm=i

return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Read NAtm from Crystal output
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RdNAtmCry(ifchk,Intact,NAtm,tag,ctmp)
implicit real(kind=8) (a-h,o-z)
character*100 :: ctmp
character*60 :: tag
logical :: Intact

i=0

tag='MOLECULAR CALCULATION'
rewind(ifchk)
do while(.true.)
  read(ifchk,"(a100)",end=1010,err=1010)ctmp
  if(index(ctmp,tag(1:21)) /= 0) exit
end do

tag='FRQFRQFRQFRQFRQFRQFRQFRQFRQFRQFRQFRQFRQFRQFRQFRQFRQFRQFRQFRQ'
do while(.true.)
  read(ifchk,"(a100)",end=1020,err=1020)ctmp
  if(index(ctmp,tag(1:60)) /= 0) exit
end do

tag='- ATOMS IN THE UNIT CELL:'
do while(.true.)
  read(ifchk,"(a100)",end=1030,err=1030)ctmp
  i = index(ctmp,tag(1:25))
  if(i /= 0)then
    read(ctmp(i+25:100),*)NAtm
    exit
  end if
end do

return
1010  call XError(Intact,"This is not a molecular calculation.")
1020  call XError(Intact,"This is not a frequency calculation.")
1030  call XError(Intact,"Failed to read #Atoms.")
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Read NAtm from Spartan output
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RdNAtmSptn(ifchk,Intact,NAtm,tag,ctmp)
implicit real(kind=8) (a-h,o-z)
character*100 :: ctmp
character*7 :: tag
logical :: Intact

i=0

tag='HESSIAN'
rewind(ifchk)
do while(.true.)
  read(ifchk,"(a100)",end=1010,err=1010)ctmp
  Idx = index(ctmp,tag(1:7))
  if(Idx /= 0) then
!  	the section of HESSIAN options (length = 0) should be skipped
    length = len_trim(ctmp(Idx+7:20))
    if(length > 0)then
      read(ifchk,"(a100)",end=1010,err=1010) ctmp
      read(ctmp,*,end=1010,err=1010) NAtm
      NAtm = NAtm/3
      return
    end if
  end if
end do

return
1010  call XError(Intact,"Failed to read #Atoms.")
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Read NAtm from PSI output
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RdNAtmPSI(ifchk,Intact,NAtm,tag,ctmp)
implicit real(kind=8) (a-h,o-z)
character*100 :: ctmp
character*18 :: tag
logical :: Intact

i=0

tag='Number of atoms is'
rewind(ifchk)
do while(.true.)
  read(ifchk,"(a100)",end=100,err=100)ctmp
  Idx = index(ctmp,tag(1:18))
  if(Idx /= 0) then
    read(ctmp(Idx+18:100),*,end=1010,err=1010) X
    NAtm = NINT(X)
    return
  end if
end do

100   tag='Number of atoms:'
rewind(ifchk)
do while(.true.)
  read(ifchk,"(a100)",end=1010,err=1010)ctmp
  Idx = index(ctmp,tag(1:16))
  if(Idx /= 0) then
    read(ctmp(Idx+16:100),*,end=1010,err=1010) NAtm
    return
  end if
end do

return
1010  call XError(Intact,"Failed to read #Atoms.")
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Read NAtm from DMol3 output
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RdNAtmDMol(ifchk,Intact,NAtm,tag,ctmp)
implicit real(kind=8) (a-h,o-z)
character*100 :: ctmp
character*9 :: tag
logical :: Intact

i=0

tag='N_atoms ='
rewind(ifchk)
do while(.true.)
  read(ifchk,"(a100)",end=1010,err=1010)ctmp
  Idx = index(ctmp,tag(1:9))
  if(Idx /= 0) then
    read(ctmp(Idx+9:100),*,end=1010,err=1010) NAtm
    return
  end if
end do

return
1010  call XError(Intact,"Failed to read #Atoms.")
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Read NAtm from ACES output
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RdNAtmACES(ifchk,Intact,NAtm,tag,ctmp)
implicit real(kind=8) (a-h,o-z)
character*100 :: ctmp
character*58 :: tag
logical :: Intact

i=0

tag='Symbol    Number           X              Y              Z'
rewind(ifchk)
do while(.true.)
  read(ifchk,"(a100)",end=1010,err=1010)ctmp
  if(index(ctmp,tag(1:58)) /= 0) exit
end do

tag='----------------------------------------------------------'
read(ifchk,*,end=1020,err=1020)
do while(.true.)
  read(ifchk,"(a100)",end=1020,err=1020)ctmp
  if(index(ctmp,tag(1:47)) /= 0) exit
  i=i+1
end do
NAtm=i

return
1010  call XError(Intact,"Failed to find Cartesian coordinates.")
1020  call XError(Intact,"Failed to read Atoms.")
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Read data from UniMoVib (ALM) data file
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RdALMode(ifchk,iout,Intact,Infred,IRaman,NAtm,ctmp,AMass,ZA,XYZ,FFx,APT,DPol,Scr)
implicit real(kind=8) (a-h,o-z)
parameter(ang2au=1.d0/0.52917720859d0)
real(kind=8) :: AMass(*),ZA(*),XYZ(*),FFx(*),APT(*),DPol(6,*),Scr(*)
character*100 :: ctmp
logical :: Intact

NAtm3 = NAtm * 3
NAtm9 = NAtm * 9
NSS = NAtm3 * NAtm3
NTT = NAtm3 * (NAtm3+1) / 2
Infred = 0
IRaman = 0

rewind(ifchk)

! read mass
read(ifchk,"(///)",err=1010,end=1010)
read(ifchk,*,err=1010,end=1010)(AMass(i),i=1,NAtm)

! read IZ
read(ifchk,*,err=1020,end=1020)
read(ifchk,*,err=1020,end=1020)(ZA(i),i=1,NAtm)

! read Cartesian coordinates in a.u. (XYZ) or in Ang. (XYZANG)
read(ifchk,"(a100)",err=1030,end=1030)ctmp
call charl2u(ctmp)
Scr(1)=1.d0
if(index(ctmp,"XYZANG") /= 0) Scr(1)=ang2au
read(ifchk,*,err=1030,end=1030)(XYZ(i),i=1,NAtm3)
! Ang --> a.u.
call AScale(NAtm3,Scr(1),XYZ,XYZ)

! read square or L.T. FFX in a.u.
read(ifchk,"(a100)",err=1050,end=1050)ctmp
call charl2u(ctmp)
if(index(ctmp,"FFXLT") == 0) then
  read(ifchk,*,err=1050,end=1050)(FFX(i),i=1,NSS)
else
  read(ifchk,*,err=1051,end=1051)(Scr(i),i=1,NTT)
  call LT2Sqr(NAtm3,Scr,FFx)
end if

! read APT in a.u.
read(ifchk,"(a100)",err=1060,end=1060)ctmp
call charl2u(ctmp)
if(index(ctmp,"NOAPT") == 0) then
  Infred = 1
  read(ifchk,*,err=1060,end=1060)(APT(i),i=1,NAtm9)
end if

! read DPR in a.u.
read(ifchk,"(a100)",err=1070,end=1070)ctmp
call charl2u(ctmp)
if(index(ctmp,"NODPR") == 0) then
  IRaman = 1
  if(index(ctmp,"DPRSQ") == 0 .and. index(ctmp,"DPR") /= 0) then
    read(ifchk,*,err=1070,end=1070)((DPol(j,i),j=1,6),i=1,NAtm3)
  else if(index(ctmp,"DPRSQ") /= 0) then
    read(ifchk,*,err=1071,end=1071)(Scr(i),i=1,NAtm3*9)
    call S9to6(NAtm3,Scr,DPol)
  end if
end if

return

1010  call XError(Intact,"Please check AMASS data!")
1020  call XError(Intact,"Please check ZA data!")
1030  call XError(Intact,"Please check XYZ data!")
1050  call XError(Intact,"Please check FFX data!")
1051  call XError(Intact,"Please check FFXLT data!")
1060  call XError(Intact,"Please check APT data!")
1070  call XError(Intact,"Please check DPR data!")
1071  call XError(Intact,"Please check DPRSQ data!")
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Read data from XYZ data file
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RdXYZ(ifchk,iout,Intact,NAtm,Elem,ZA,XYZ,FFx,S1,S2,S3,S4)
implicit real(kind=8) (a-h,o-z)
parameter(ang2au=1.d0/0.52917720859d0)
real(kind=8) :: ZA(*),XYZ(3,*),FFx(*),S1(*),S2(*),S3(*),S4(*)
character*3 :: Elem
logical :: Intact

NAtm3 = NAtm * 3

rewind(ifchk)

read(ifchk,"(/)",err=1010,end=1010)
do i=1,NAtm
  read(ifchk,*,err=1010,end=1010)Elem,XYZ(:,i)
  call ElemZA(0,Elem,Elem,ZA(i))
end do

! Ang to a.u.
call AScale(NAtm3,ang2au,XYZ,XYZ)

! FFX is simulated by CZ * [E"(NRE) * E"(NRE)]^(1/2), which is positive definite.
! CZ = 0.5/Max(ZA)
CZ = 0.5d0 / ArMax(NAtm,i,ZA)
call DDerNRE(NAtm,ZA,XYZ,FFx,S1(1),S1(6))
call MPACMF(FFx,FFx,S1,NAtm3,NAtm3,1)
call SqrtMp(Intact,1,NAtm3,S1,FFx,FFx,S2,S3,S4)
call AScale(NAtm3*NAtm3,CZ,FFx,FFx)

return
1010  call XError(Intact,"Please check XYZ file!")
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Read data from Gaussian fchk
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RdGauss(ifchk,iout,tag,ctmp,Intact,IRdMas,Infred,IRaman,NAtm,AMass,ZA,XYZ,FFx,APT,DPol,Scr)
implicit real(kind=8) (a-h,o-z)
logical :: Intact

real(kind=8) :: AMass(*),ZA(*),XYZ(*),FFx(*),APT(3,*),DPol(6,*),Scr(*)
character*100 :: ctmp
character*49 :: tag

NAtm3 = NAtm * 3
NTT = NAtm3*(NAtm3+1)/2
Infred=0
IRaman=0

!! read nuclear charges; they lead to errors in the case of ECP
!tag='Nuclear charges                            R   N='
!     read Atomic numbers
tag='Atomic numbers                             I   N='
rewind(ifchk)
001   read(ifchk,"(a100)",end=010)ctmp
if(index(ctmp,tag)==0)goto 001
read(ifchk,*)(ZA(i),i=1,NAtm)

! read Cartesian coordinates in standard orientation (a.u.)
tag='Current cartesian coordinates              R   N='
rewind(ifchk)
101   read(ifchk,"(a100)",end=110)ctmp
if(index(ctmp,tag)==0)goto 101
read(ifchk,"(5e16.8)")(XYZ(i),i=1,NAtm3)

! read FFx (a.u.)
tag='Cartesian Force Constants                  R   N='
rewind(ifchk)
201   read(ifchk,"(a100)",end=210)ctmp
if(index(ctmp,tag)==0)goto 201
read(ifchk,"(5e16.8)")(Scr(i),i=1,NTT)
call LT2Sqr(NAtm3,Scr,FFx)

! read APT (a.u.), optional
tag='Dipole Derivatives                         R   N='
rewind(ifchk)
301   read(ifchk,"(a100)",end=350)ctmp
if(index(ctmp,tag)==0) goto 301
read(ifchk,"(5e16.8)")((APT(j,i),j=1,3),i=1,NAtm3)
Infred= 1

! read DPol (a.u.), optional
350   tag='Polarizability Derivatives                 R   N='
rewind(ifchk)
351   read(ifchk,"(a100)",end=400)ctmp
if(index(ctmp,tag)==0) goto 351
read(ifchk,"(5e16.8)")((DPol(j,i),j=1,6),i=1,NAtm3)
IRaman= 1

! read atomic masses (1): masses were saved by freq(SaveNormalModes)
! G09 and higher versions only!
400   if(IRdMas == 2) goto 600
tag='Vib-AtMass                                 R   N='
rewind(ifchk)
401   read(ifchk,"(a100)",end=500)ctmp
if(index(ctmp,tag)==0)goto 401
  read(ifchk,"(5e16.8)")(AMass(i),i=1,NAtm)
goto 600

500   AMass(1)=-1.d0

600   continue

return
010   call XError(Intact,"No nuclear charges found!")
110   call XError(Intact,"No Cartesian coordinates found!")
210   call XError(Intact,"No Cartesian force constants found!")
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Read data from *.dat + *.out of Gamess or PUNCH + *.out of Firefly
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RdGAMES(ifchk,iout,tmp,ctmp,Intact,Infred,NAtm,AMass,ZA,XYZ,FFx,APT)
implicit real(kind=8) (a-h,o-z)
parameter(an2br=0.52917720859d0,au2deb=2.541746d0)
real(kind=8) :: AMass(*),ZA(*),XYZ(3,*),FFx(NAtm*3,*),APT(3,*)
character*100 :: ctmp
character*56 :: tmp
logical :: Intact

NAtm3 = NAtm * 3
Infred= 0

! read nuclear charges and Cartesian coordinates in standard orientation (a.u.)
rewind(iout)
tmp='ATOM      ATOMIC                      COORDINATES (BOHR)'
101   read(iout,"(a100)",end=110,err=120)ctmp
if(index(ctmp,tmp)==0)goto 101

read(iout,*)
do i=1,NAtm
  read(iout,"(a100)",end=120)ctmp
  read(ctmp,*,err=120)tmp,ZA(i),XYZ(1,i),XYZ(2,i),XYZ(3,i)
end do

! read Cartesian Force Constant matrix (a.u.)
rewind(ifchk)
201   read(ifchk,"(a6)",end=210)ctmp
if(index(ctmp,' $HESS')==0)goto 201
read(ifchk,*)
do i=1,NAtm3
  read(ifchk,"(5x,5e15.8)",err=210,end=210)(FFx(j,i),j=1,NAtm3)
end do

! read APT (a.u.)
! rewind(ifchk)
301   read(ifchk,"(a7)",end=401)ctmp
if(index(ctmp,' $DIPDR')==0)goto 301
read(ifchk,"(1x,3e15.8)")((APT(j,i),j=1,3),i=1,NAtm3)
! Deb/Ang --> a.u.
factor=an2br/au2deb
call AScale(NAtm3*3,factor,APT,APT)
Infred= 1

! read atomic masses
401   rewind(ifchk)
read(ifchk,"(a13)",end=410)ctmp
if(index(ctmp,'ATOMIC MASSES')==0)goto 401
read(ifchk,*,err=410,end=410)(AMass(i),i=1,NAtm)

return
110   call XError(Intact,"No Cartesian coordinates found!")
120   call XError(Intact,"Cartesian coordinates are wrong!")
210   call XError(Intact,"No Cartesian force constants found!")
410   call XError(Intact,"No atomic masses found!")
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Read data from ORCA *.hess
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RdORCA(ifchk,ctmp,Intact,Infred,NAtm,AMass,ZA,XYZ,FFx,APT)
implicit real(kind=8) (a-h,o-z)
real(kind=8) :: AMass(*),ZA(*),XYZ(3,*),FFx(NAtm*3,*),APT(3,*)
character*100 :: ctmp
logical :: Intact

NAtm3 = NAtm * 3
Infred= 0

! read nuclear charges, atomic masses, and Cartesian coordinates in standard orientation (a.u.)
rewind(ifchk)
001   read(ifchk,"(a6)",end=010)ctmp
if(index(ctmp,"$atoms")==0)goto 001
read(ifchk,*)
do i=1,NAtm
  read(ifchk,"(1x,a3,f10.4,1x,3f13.6)")ctmp,AMass(i),(XYZ(j,i),j=1,3)
  call ElemZA(0,ctmp,ctmp,ZA(i))
end do

! read FFx (a.u.)
rewind(ifchk)
101   read(ifchk,"(a8)",end=110)ctmp
if(index(ctmp,"$hessian")==0)goto 101
read(ifchk,*)
NBlock=(NAtm3-1)/6+1
do i=1,NBlock
  iv1=(i-1)*6+1
  iv2=min(i*6,NAtm3)
  read(ifchk,*)
  do j=1,NAtm3
    read(ifchk,"(11x,6f11.6)")(FFx(k,j),k=iv1,iv2)
  end do
end do

! read APT (a.u.)
rewind(ifchk)
201   read(ifchk,"(a19)",end=900)ctmp
if(index(ctmp,"$dipole_derivatives")==0)goto 201
read(ifchk,*)
read(ifchk,"(3f13.6)")((APT(j,i),j=1,3),i=1,NAtm3)
Infred= 1

900  continue
return
010  call XError(Intact,"No $atoms found!")
110  call XError(Intact,"No $hessian found!")
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! CFour: check whether this is an analytic frequency calculation
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine ChkAFRQ(ifchk,tag,ctmp,Intact)
implicit real(kind=8) (a-h,o-z)
character*100 :: ctmp
character*32 :: tag
logical :: Intact

tag='       VIBRATION            IVIB'
rewind(ifchk)
001   read(ifchk,"(a100)",end=010)ctmp
if(index(ctmp,tag) == 0)goto 001
if(index(ctmp,"ANALYTIC") == 0) call XError(Intact,"This is not an analytic frequency calculation!")

return
010   call XError(Intact,"Is this a CFour freq. calculation?")
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Read data from CFour *.out
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RdCFour(ifchk,igeom,tag,ctmp,Intact,Infred,NAtm,AMass,ZA,XYZ,FFx,APT,Scr)
implicit real(kind=8) (a-h,o-z)
real(kind=8) :: AMass(*),ZA(*),XYZ(3,*),FFx(NAtm*3,*),APT(3,*),Scr(*)
character*100 :: ctmp
character*77 :: tag,tag1
logical :: Intact,ifc1

NAtm3 = NAtm * 3
ifc1 = .False.
Infred= 0

! read nuclear charges and Cartesian coordinates (a.u.)
rewind(igeom)
read(igeom,*)
do i=1,NAtm
! read(igeom,"(4f20.10)")ZA(i),(XYZ(j,i),j=1,3)
  read(igeom,*)ZA(i),(XYZ(j,i),j=1,3)
end do
call ACopy(NAtm3,XYZ,Scr)

! C1 symmetry?
tag='   The computational point group is '
rewind(ifchk)
001   read(ifchk,"(a100)",end=030)ctmp
if(index(ctmp,tag(1:36))==0)goto 001
if(ctmp(37:39) == "C1 ") ifc1 = .True.

! read unordered Cartesian coordinates (a.u.)
004   tag=' Z-matrix   Atomic            Coordinates (in bohr)'
! rewind(ifchk)
005   read(ifchk,"(a100)",end=010)ctmp
if(index(ctmp,tag(1:51))==0)goto 005
i=0
tag=' --------------------------------------------------'
read(ifchk,"(/)")
do while(.true.)
  read(ifchk,"(a100)",end=020)ctmp
  if(index(ctmp,tag(1:51)) > 0)exit
  read(ctmp,"(12x,i4)")iz
  if(iz <= 0) cycle
  i=i+1
  if(i <= NAtm) read(ctmp,"(20x,3f15.8)")(XYZ(j,i),j=1,3)
end do
if(i > NAtm) goto 020

! read atomic masses
tag='  masses used (in AMU) in vibrational analysis:'
! rewind(ifchk)
101   read(ifchk,"(a100)",end=110)ctmp
if(index(ctmp,tag)==0)goto 101
read(ifchk,*)(AMass(i),i=1,NAtm)

! reorder atomic masses, and recover reordered Cartesian coordinates
call ROMass(Intact,NAtm,Scr,XYZ,AMass,Scr(NAtm3+1))
call ACopy(NAtm3,Scr,XYZ)

! read FFx (a.u.)
tag='                            Molecular hessian'
rewind(ifchk)
401   read(ifchk,"(a100)",end=410)ctmp
if(index(ctmp,tag)==0)goto 401
! if(ifc1) goto 403
tag ='             #1 x        #1 y        #1 z        #1 x        #1 y        #1 z'
tag1='             #1 x        #1 y        #1 z        #2 x        #2 y        #2 z'
! if not C1 symmetry, skip hessian in each symmetry
402   read(ifchk,"(a100)",end=411)ctmp
if(index(ctmp,'Symmetry')/=0)then
  read(ifchk,*)
  read(ifchk,*)
  goto 402
end if
ctmp(12:13)="  "
ctmp(24:25)="  "
ctmp(36:37)="  "
ctmp(48:49)="  "
ctmp(60:61)="  "
ctmp(72:73)="  "
if(index(ctmp,tag)==0 .and. index(ctmp,tag1)==0) goto 402

403   NBlock=(NAtm3-1)/6+1
do i=1,NBlock
  iv1=(i-1)*6+1
  iv2=min(i*6,NAtm3)
  do j=iv1,NAtm3
    if(mod(j,3) == 1)read(ifchk,*)
    read(ifchk,"(7x,6f12.6)",end=420)(FFx(k,j),k=iv1,min(j,iv2))
  end do
  read(ifchk,"(//)",end=420)
end do
do i=1,NAtm3
  do j=1,i-1
    FFx(i,j)=FFx(j,i)
  end do
end do

! read APT (a.u.)
tag='                     Total dipole moment derivatives'
! rewind(ifchk)
501   read(ifchk,"(a100)",end=900)ctmp
if(index(ctmp,tag)==0)goto 501

if(ifc1)then
  tag='                 Ex             Ey             Ez'
else
  tag='                       Ex          Ey          Ez'
end if
502   read(ifchk,"(a100)",end=900)ctmp
if(index(ctmp,tag)==0)goto 502
do i=1,NAtm3
  if(mod(i,3) == 1)read(ifchk,*)
  if(ifc1)then
    read(ifchk,"(21x,3f15.8)")(APT(j,i),j=1,3)
  else
    read(ifchk,"(15x,3f12.6)")(APT(j,i),j=1,3)
  end if
end do
Infred= 1

900  continue

return
010   call XError(Intact,"No Cartesian coordinates found!")
020   call XError(Intact,"Please check Cartesian coordinates!")
030   call XError(Intact,"Unknown group symmetry!")
110   call XError(Intact,"No atomic masses found!")
410   call XError(Intact,"No Cartesian force constants found (1)!")
411   call XError(Intact,"No Cartesian force constants found (2)!")
412   call XError(Intact,"No Cartesian force constants found (3)!")
420   call XError(Intact,"Please check Force Constant matrix!")
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! CFour: reorder atomic masses according to old and new geometries.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine ROMass(Intact,NAtm,Xnew,Xold,AMass,AMtmp)
implicit real(kind=8) (a-h,o-z)
parameter(tol=1.d-6)
real(kind=8) :: Xnew(3,*),Xold(3,*),AMass(*),AMtmp(*)
logical :: Intact

call ACopy(NAtm,AMass,AMtmp)
do i=1,NAtm
  do j=1,NAtm
    x = distance(Xnew(1,i),Xold(1,j))
    if(x < tol) then
      AMass(i) = AMtmp(j)
      exit
    end if
  end do
end do

return
call XError(Intact,"Cartesian coordinates are different!")
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Read data from Molpro *.out
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RdMolp(ifchk,tag,ctmp,Intact,Infred,NAtm,AMass,ZA,XYZ,FFx,APT)
implicit real(kind=8) (a-h,o-z)
parameter(au2deb=2.541746d0,au2ang=0.52917720859d0)
real(kind=8) :: AMass(*),ZA(*),XYZ(3,*),FFx(NAtm*3,*),APT(3,*)
character*100 :: ctmp
character*18 :: hesstag
character*61 :: tag
logical :: Intact

NAtm3 = NAtm * 3
Infred= 0
hesstag=' PROGRAM * HESSIAN'

! read nuclear charges and Cartesian coordinates (a.u.)
! NOTE:
! If ECP basis set is used, effective charges will be printed in the output file.
tag='  Nr  Atom  Charge       X              Y              Z'
rewind(ifchk)
001   read(ifchk,"(a100)",end=010)ctmp
if(index(ctmp,hesstag)==0)goto 001
002   read(ifchk,"(a100)",end=020)ctmp
if(index(ctmp,tag)==0)goto 002
read(ifchk,*)
do i=1,NAtm
  read(ifchk,"(7x,a3,7x,3f15.9)",end=030)ctmp,(XYZ(j,i),j=1,3)
  call rmnumb(3,ctmp)
  call ElemZA(0,ctmp,ctmp,ZA(i))
end do

! read FFx (a.u.)
tag=' Force Constants (Second Derivatives of the Energy) in [a.u.]'
rewind(ifchk)
101   read(ifchk,"(a100)")ctmp
if(index(ctmp,hesstag)==0)goto 101
102   read(ifchk,"(a100)",end=110)ctmp
if(index(ctmp,tag)==0)goto 102
NBlock=(NAtm3-1)/5+1
do i=1,NBlock
  iv1=(i-1)*5+1
  iv2=min(i*5,NAtm3)
  read(ifchk,*)
  do j=iv1,NAtm3
    read(ifchk,"(22x,5f12.7)",end=120)(FFx(k,j),k=iv1,min(j,iv2))
  end do
end do
do i=1,NAtm3
  do j=1,i-1
    FFx(i,j)=FFx(j,i)
  end do
end do

! read APT (a.u.)
tag=' Dipole Moment Derivatives [debye/ang]'
rewind(ifchk)
201   read(ifchk,"(a100)")ctmp
if(index(ctmp,hesstag)==0)goto 201
202   read(ifchk,"(a100)",end=300)ctmp
if(index(ctmp,tag)==0)goto 202
NBlock=(NAtm3-1)/8+1
do i=1,NBlock
  iv1=(i-1)*8+1
  iv2=min(i*8,NAtm3)
  read(ifchk,*)
  do j=1,3
    read(ifchk,"(10x,8f14.7)",end=300)(APT(j,k),k=iv1,iv2)
  end do
end do
! debye/ang --> a.u.
call AScale(NAtm3*3,au2ang/au2deb,APT,APT)
Infred= 1

300  continue
! read atomic masses
tag=' Atomic Masses'
rewind(ifchk)
301   read(ifchk,"(a100)")ctmp
if(index(ctmp,hesstag)==0)goto 301
302   read(ifchk,"(a100)",end=310)ctmp
if(index(ctmp,tag)==0)goto 302
read(ifchk,"(8x,10f12.6)")(AMass(i),i=1,NAtm)

return
010   call XError(Intact,"No 'PROGRAM * HESSIAN' found!")
020   call XError(Intact,"No Cartesian coordinates found!")
030   call XError(Intact,"Please check Cartesian coordinates!")
110   call XError(Intact,"No Cartesian force constants found!")
120   call XError(Intact,"Please check Force Constant matrix!")
210   call XError(Intact,"No dipole derivatives found!")
220   call XError(Intact,"Please check Dipole Moment Derivatives!")
310   call XError(Intact,"No atomic masses found!")
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Read data from QChem fchk
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RdQChem(ifchk,tag,ctmp,Intact,IRdMas,NAtm,AMass,ZA,XYZ,FFx,APT,Scr)
implicit real(kind=8) (a-h,o-z)
real(kind=8) :: AMass(*),ZA(*),XYZ(*),FFx(*),APT(3,*),Scr(*)
character*100 :: ctmp
character*49 :: tag
logical :: Intact

NAtm3 = NAtm * 3
NTT = NAtm3*(NAtm3+1)/2

! read nuclear charges
tag='Nuclear charges                            R   N='
rewind(ifchk)
001   read(ifchk,"(a100)",end=010)ctmp
if(index(ctmp,tag)==0)goto 001
read(ifchk,"(5e16.8)")(ZA(i),i=1,NAtm)

! read Cartesian coordinates in standard orientation (a.u.)
tag='Current cartesian coordinates              R   N='
rewind(ifchk)
101   read(ifchk,"(a100)",end=110)ctmp
if(index(ctmp,tag)==0)goto 101
read(ifchk,"(5e16.8)")(XYZ(i),i=1,NAtm3)

! read FFx (a.u.)
tag='Cartesian Force Constants                  R   N='
rewind(ifchk)
201   read(ifchk,"(a100)",end=210)ctmp
if(index(ctmp,tag)==0)goto 201
read(ifchk,"(5e16.8)")(Scr(i),i=1,NTT)
call LT2Sqr(NAtm3,Scr,FFx)

! read APT (a.u.)
! not available!!!

AMass(1)=-1.d0

return
010   call XError(Intact,"No nuclear charges found!")
110   call XError(Intact,"No Cartesian coordinates found!")
210   call XError(Intact,"No Cartesian force constants found!")
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Read data of NWChem
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RdNWChem(ifchk,ihess,iddip,tag,ctmp,Intact,Infred,NAtm,AMass,ZA,XYZ,FFx,APT,Scr)
implicit real(kind=8) (a-h,o-z)
real(kind=8) :: AMass(*),ZA(*),XYZ(3,*),FFx(*),APT(*),Scr(*)
character*100 :: ctmp
character*69 :: tag
logical :: Intact,ifopen

NAtm3 = NAtm * 3
Infred= 0

! read nuclear charges, Cartesian coordinates (a.u.), and atomic masses
tag='     atom    #        X              Y              Z            mass'
rewind(ifchk)
001   read(ifchk,"(a100)",end=010)ctmp
if(index(ctmp,tag)==0)goto 001
read(ifchk,*)
do i=1,NAtm
  read(ifchk,"(4x,a3,7x,4d15.6)",end=020)ctmp,(XYZ(j,i),j=1,3),AMass(i)
  call ElemZA(0,ctmp,ctmp,ZA(i))
end do

! read FFx (a.u.)
NTT=NAtm3*(NAtm3+1)/2
rewind(ihess)
read(ihess,*,err=110)(Scr(i),i=1,NTT)
call LT2Sqr(NAtm3,Scr,FFx)

! read APT (a.u.)
inquire(unit=iddip,opened=ifopen)
if(ifopen)then
  rewind(iddip)
  read(iddip,*,err=210)(APT(i),i=1,NAtm3*3)
  Infred= 1
end if

return
010   call XError(Intact,"No Cartesian coordinates found!")
020   call XError(Intact,"Please check Cartesian coordinates!")
110   call XError(Intact,"Please check Force Constant matrix!")
210   call XError(Intact,"Please check Dipole Moment Derivatives!")
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Read data from GAMESS-UK *.out
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RdGMUK(ifchk,tag,ctmp,Intact,Infred,NAtm,AMass,ZA,XYZ,FFx,APT)
implicit real(kind=8) (a-h,o-z)
parameter(au2deb=2.541746d0,au2ang=0.52917720859d0,ang2au=1.d0/au2ang)
real(kind=8) :: AMass(*),ZA(*),XYZ(3,*),FFx(NAtm*3,*),APT(3,*)
character*100 :: ctmp
character*99 :: tag
logical :: Intact

NAtm3 = NAtm * 3
Infred= 0

! read nuclear charges, Cartesian coordinates (a.u.), and  atomic masses
tag='     atom              x                   y                   z             nucleus  atomic weight'
rewind(ifchk)
001   read(ifchk,"(a100)",end=010)ctmp
if(index(ctmp,tag(1:99))==0)goto 001
read(ifchk,*)
do i=1,NAtm
  read(ifchk,"(14x,3f20.10,2x,a2,4x,f17.10)",end=020) (XYZ(j,i),j=1,3),ctmp,AMass(i)
  if(ctmp(1:1) == ' ')then
    ctmp(1:1)=ctmp(2:2)
    ctmp(2:2)=' '
  end if
  ctmp(3:3)=' '
  call ElemZA(0,ctmp,ctmp,ZA(i))
end do
! ang --> a.u.
call AScale(NAtm3,ang2au,XYZ,XYZ)

! read FFx (a.u.)
tag='                    * total cartesian 2nd derivative matrix *'
rewind(ifchk)
101   read(ifchk,"(a100)",end=110)ctmp
if(index(ctmp,tag(1:61))==0)goto 101
read(ifchk,*)
NBlock=(NAtm3-1)/9+1
do i=1,NBlock
  iv1=(i-1)*9+1
  iv2=min(i*9,NAtm3)
  do j=1,9
    read(ifchk,*)
  end do
  do j=1,NAtm3
    read(ifchk,"(19x,9f9.5)",end=120)(FFx(k,j),k=iv1,iv2)
  end do
end do

! read APT (a.u.)
tag='                    *     (debye/angstrom)     *'
! rewind(ifchk)
201   read(ifchk,"(a100)",end=300)ctmp
if(index(ctmp,tag)==0)goto 201
do j=1,3
  read(ifchk,*)
end do
NBlock=(NAtm3-1)/8+1
do i=1,NAtm3
  if(mod(i,3) == 1)then
    read(ifchk,*)
    read(ifchk,*)
  end if
  read(ifchk,"(23x,3f16.8)",end=300)(APT(j,i),j=1,3)
end do
! debye/ang --> a.u.
call AScale(NAtm3*3,au2ang/au2deb,APT,APT)
Infred= 1

300  continue

return
010   call XError(Intact,"No Cartesian coordinates found!")
020   call XError(Intact,"Please check Cartesian coordinates!")
110   call XError(Intact,"No Cartesian force constants found!")
120   call XError(Intact,"Please check Force Constant matrix!")
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Read data from Turbomole's aoforce.out and dipgrad
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RdTurbm(ifchk,iddip,tag,ctmp,Intact,Infred,NAtm,AMass,ZA,XYZ,FFx,APT)
implicit real(kind=8) (a-h,o-z)
real(kind=8) :: AMass(*),ZA(*),XYZ(3,*),FFx(NAtm*3,*),APT(*)
character*100 :: ctmp
character*72 :: tag
logical :: Intact,ifopen

NAtm3 = NAtm * 3
Infred= 0

! read nuclear charges and Cartesian coordinates (a.u.)
tag='                    atomic coordinates            atom    charge  isotop'
rewind(ifchk)
001   read(ifchk,"(a100)",end=010)ctmp
if(index(ctmp,tag)==0)goto 001
do i=1,NAtm
  read(ifchk,"(6x,3f14.8,8x,f8.4)",end=020)(XYZ(j,i),j=1,3),ZA(i)
end do

! read atomic masses
tag="            ('*' denotes special isotop !)"
! rewind(ifchk)
101   read(ifchk,"(a100)",end=110)ctmp
if(index(ctmp,tag(1:42))==0)goto 101
do i=1,NAtm
  read(ifchk,"(36x,f9.5)",err=120,end=120)AMass(i)
end do

! read FFx (a.u.)
tag='          CARTESIAN FORCE CONSTANT MATRIX (hartree/bohr**2)'
! rewind(ifchk)
201   read(ifchk,"(a100)",end=210)ctmp
if(index(ctmp,tag(1:59))==0)goto 201
read(ifchk,*)
read(ifchk,*)
NBlock=(NAtm3-1)/6+1
do i=1,NBlock
  iv1=(i-1)*6+1
  iv2=min(i*6,NAtm3)
  read(ifchk,*)
  read(ifchk,*)
  read(ifchk,*)
  do j=iv1,NAtm3
    read(ifchk,"(14x,6f10.7)",err=220,end=220) (FFx(k,j),k=iv1,min(j,iv2))
  end do
end do
do i=1,NAtm3
  do j=1,i-1
    FFx(i,j)=FFx(j,i)
  end do
end do

! read APT (a.u.)
inquire(unit=iddip,opened=ifopen)
if(ifopen)then
  rewind(iddip)
  read(iddip,*)
  read(iddip,*,err=310)(APT(i),i=1,NAtm3*3)
  Infred= 1
else
  call AClear(NAtm3*3,APT)
end if

return
010   call XError(Intact,"No Cartesian coordinates found!")
020   call XError(Intact,"Please check Cartesian coordinates!")
110   call XError(Intact,"No atomic masses found!")
120   call XError(Intact,"Please check atomic masses!")
210   call XError(Intact,"No Cartesian force constants found!")
220   call XError(Intact,"Please check Force Constant matrix!")
310   call XError(Intact,"Please check Dipole Moment Derivatives!")
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Read data from deMon output
! NOTE: FFX is in input orientation, but the Cartesian coordinates in input orientation are printed only at the very beginning and
! the end of OPT calculation.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RdDeMon(ifchk,tag,ctmp,Intact,NAtm,AMass,ZA,XYZ,FFx)
implicit real(kind=8) (a-h,o-z)
parameter(au2ang=0.52917720859d0)
real(kind=8) :: AMass(*),ZA(*),XYZ(3,*),FFx(NAtm*3,*)
character*100 :: ctmp
character*44 :: tag
logical :: Intact

NAtm3 = NAtm * 3

! read Cartesian coordinates in input orientation (a.u.), Z, and mass
!
! search optimized coordinates first in OPT + FREQ calculation
tag=' FINAL INPUT ORIENTATION IN ANGSTROM'
rewind(ifchk)
11    read(ifchk,"(a100)",end=20)ctmp
if(index(ctmp,tag)==0)goto 11
goto 90

! FREQ with Z-MAT
20    tag=' Z-MATRIX INPUT ORIENTATION IN ANGSTROM'
rewind(ifchk)
21    read(ifchk,"(a100)",end=30)ctmp
if(index(ctmp,tag)==0)goto 21
goto 90

! FREQ with CART.
30    tag=' INPUT ORIENTATION IN ANGSTROM'
rewind(ifchk)
31    read(ifchk,"(a100)",end=40)ctmp
if(index(ctmp,tag)==0)goto 31

90    read(ifchk,"(//)")
do i=1,NAtm
! Note: the formats of ver.3.x and ver.4.x are a little different!
  read(ifchk,"(10x,a100)")ctmp
  read(ctmp,*)(XYZ(j,i),j=1,3),IZ,AMass(i)
  ZA(i)=dble(IZ)
end do
! ang --> a.u.
call AScale(NAtm3,1.d0/au2ang,XYZ,XYZ)

! read tag of frequency calculation
tag=' *** MOLECULAR FREQUENCY ANALYSIS ***'
! rewind(ifchk)
101   read(ifchk,"(a100)",end=110)ctmp
if(index(ctmp,tag)==0)goto 101

! read FFx (a.u.) in input orientation
tag=' *** FORCE CONSTANTS [HARTREE/(BOHR**2)] ***'
201   read(ifchk,"(a100)",end=210)ctmp
if(index(ctmp,tag)==0)goto 201
read(ifchk,"(/)")

NBlock=(NAtm3-1)/4+1
do i=1,NBlock
  iv1=(i-1)*4+1
  iv2=min(i*4,NAtm3)
  do j=1,NAtm3
    read(ifchk,"(8x,4f17.10)",end=220)(FFx(k,j),k=iv1,iv2)
  end do
  read(ifchk,"(///)")
end do

return
40    call XError(Intact,"No Cartesian coordinates found!")
110   call XError(Intact,"No tag of frequency calculation found!")
210   call XError(Intact,"No Cartesian force constants found!")
220   call XError(Intact,"Please check Force Constant matrix!")
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Read data of PQS
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RdPQS(ifchk,ihess,iddip,tag,ctmp,Intact,Infred,NAtm,AMass,ZA,XYZ,FFx,APT,Scr)
implicit real(kind=8) (a-h,o-z)
real(kind=8) :: AMass(*),ZA(*),XYZ(3,*),FFx(*),APT(*),Scr(*)
character*100 :: ctmp
character*20 :: tag
logical :: Intact,ifopen

NAtm3 = NAtm * 3
Infred= 0

! read nuclear charges, Cartesian coordinates (a.u.), and atomic masses
tag='$coordinates'
rewind(ifchk)
001   read(ifchk,"(a100)",end=010)ctmp
if(index(ctmp,tag)==0)goto 001
do i=1,NAtm
  read(ifchk,*,end=020)ctmp,(XYZ(j,i),j=1,3),Zval,AMass(i)
  call ElemZA(0,ctmp,ctmp,ZA(i))
end do

! read FFx (a.u.)
NTT=NAtm3*(NAtm3+1)/2
tag='$hessian'
rewind(ihess)
101   read(ihess,"(a100)",end=110)ctmp
if(index(ctmp,tag)==0)goto 101
read(ihess,*,err=110)ctmp
read(ihess,*,err=110)(Scr(i),i=1,NTT)
call LT2Sqr(NAtm3,Scr,FFx)

! read APT (a.u.)
inquire(unit=iddip,opened=ifopen)
if(ifopen)then
  tag='$dipole derivatives'
  rewind(iddip)
  201  read(iddip,"(a100)",end=210)ctmp
  if(index(ctmp,tag)==0)goto 201
  read(iddip,*,err=210)(APT(i),i=1,NAtm3*3)
  Infred= 1
else
  call AClear(NAtm3*3,APT)
end if

return
010   call XError(Intact,"No Cartesian coordinates found!")
020   call XError(Intact,"Please check Cartesian coordinates!")
110   call XError(Intact,"Please check Force Constant matrix!")
210   call XError(Intact,"Please check Dipole Moment Derivatives!")
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Read data from MOPAC *.out
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RdMopac(ifchk,tag,ctmp,Intact,NAtm,ZA,XYZ,FFx)
implicit real(kind=8) (a-h,o-z)
parameter(dy2au=1.d0/15.56893d0,ang2au=1.d0/0.52917720859d0)
real(kind=8) :: ZA(*),XYZ(3,*),FFx(NAtm*3,*)
character*100 :: ctmp
character*44 :: tag
logical :: Intact,IFV70

NAtm3 = NAtm * 3
IFV70 = .false.

! read nuclear charges and Cartesian coordinates (a.u.)
tag='ORIENTATION OF MOLECULE IN FORCE CALCULATION'
rewind(ifchk)
001   read(ifchk,"(a100)",end=010)ctmp
if(index(ctmp,tag)==0)goto 001

read(ifchk,"(//)",end=020)
do i=1,NAtm
  read(ifchk,"(15x,a100)",err=020,end=020)ctmp
! MOPAC 7.0 or higher versions?
  read(ctmp,*,err=005)IZ

! Version = 6.0 or 7.0
  IFV70 = .True.
  read(ctmp,*,err=020,end=020)IZA,(XYZ(j,i),j=1,3)
  ZA(i) = dble(IZA)
  cycle

! Version > 7.0
  005  if(IFV70) goto 020
  call rmnumb(3,ctmp)
  ctmp(1:3) = adjustl(ctmp(1:3))
  call ElemZA(0,ctmp,ctmp,ZA(i))
  read(ctmp(4:100),*,err=020,end=020)(XYZ(j,i),j=1,3)
end do
! Ang --> a.u.
call AScale(NAtm3,ang2au,XYZ,XYZ)

! read FFx (a.u.)
tag='FULL FORCE MATRIX, INVOKED BY "DFORCE"'
! rewind(ifchk)
101   read(ifchk,"(a100)",end=110)ctmp
if(index(ctmp,tag)==0)goto 101

NBlock=(NAtm3-1)/6+1
do i=1,NBlock
  iv1=(i-1)*6+1
  iv2=min(i*6,NAtm3)
! skip 3 or 4 lines
  read(ifchk,"(//)",end=120)
  if(.not. IFV70)read(ifchk,*,end=120)

  do j=iv1,NAtm3
    read(ifchk,"(11x,a100)",end=120)ctmp
    read(ctmp,*,end=120)(FFx(k,j),k=iv1,min(j,iv2))
  end do
end do
do i=1,NAtm3
  do j=1,i-1
    FFx(i,j)=FFx(j,i)
  end do
end do
! mdyn/Ang --> a.u.
if(IFV70) then
  call AScale(NAtm3*NAtm3,dy2au*2.d0,FFx,FFx)
else
  call AScale(NAtm3*NAtm3,dy2au,FFx,FFx)
end if

return
010   call XError(Intact,"No geom. in standard orientation found!")
020   call XError(Intact,"Please check Cartesian coordinates!")
110   call XError(Intact,"No Cartesian force constants found!")
120   call XError(Intact,"Please check Force Constant matrix!")
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Read data from AMPAC(2.x)/AMSOL *.out
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RdAMPAC(ifchk,tag,ctmp,Intact,NAtm,ZA,XYZ,FFx)
implicit real(kind=8) (a-h,o-z)
parameter(dy2au=1.d0/15.56893d0,ang2au=1.d0/0.52917720859d0)
real(kind=8) :: ZA(*),XYZ(3,*),FFx(NAtm*3,*)
character*100 :: ctmp
character*44 :: tag
logical :: Intact

NAtm3 = NAtm * 3

! read nuclear charges and Cartesian coordinates (a.u.)
tag='ORIENTATION OF MOLECULE IN FORCE CALCULATION'
rewind(ifchk)
001   read(ifchk,"(a100)",end=010)ctmp
if(index(ctmp,tag)==0)goto 001

read(ifchk,"(//)",end=020)
do i=1,NAtm
  read(ifchk,"(a100)",err=020,end=020)ctmp
  read(ctmp(7:),*,err=020,end=020)ZA(i),(XYZ(j,i),j=1,3)
end do
! Ang --> a.u.
call AScale(NAtm3,ang2au,XYZ,XYZ)

! read FFx (a.u.)
tag='FULL FORCE MATRIX, INVOKED BY "DFORCE"'
! rewind(ifchk)
101   read(ifchk,"(a100)",end=110)ctmp
if(index(ctmp,tag)==0)goto 101

NBlock=(NAtm3-1)/6+1
do i=1,NBlock
  iv1=(i-1)*6+1
  iv2=min(i*6,NAtm3)
! skip 3 lines
  read(ifchk,"(//)",end=120)

  do j=iv1,NAtm3
    read(ifchk,"(6x,a100)",end=120)ctmp
    read(ctmp,*,end=120)Itmp,(FFx(k,j),k=iv1,min(j,iv2))
  end do
end do
do i=1,NAtm3
  do j=1,i-1
    FFx(i,j)=FFx(j,i)
  end do
end do
! mdyn/Ang --> a.u.
call AScale(NAtm3*NAtm3,dy2au*2.d0,FFx,FFx)

return
010   call XError(Intact,"No geom. in standard orientation found!")
020   call XError(Intact,"Please check Cartesian coordinates!")
110   call XError(Intact,"No Cartesian force constants found!")
120   call XError(Intact,"Please check Force Constant matrix!")
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Read data from Dalton output
!
! 2015.12.23: analytic Hessian and APT of Dalton2015 are supported.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RdDalton(ifchk,tag,ta2,ctmp,Intact,Infred,NAtm,AMass,ZA,XYZ,FFx,APT)
implicit real(kind=8) (a-h,o-z)
real(kind=8) :: AMass(*),ZA(*),XYZ(3,*),FFx(NAtm*3,*),APT(3,*)
character*100 :: ctmp
character*44 :: tag,ta2
logical :: Intact,IfAna,ifgeom

NAtm3 = NAtm * 3
IfAna = .false.
ifgeom= .false.
Infred= 0

rewind(ifchk)

! read the Cartesian coordinates in a.u.
! in the case of opt+freq, the final coordinates should be read.
tag='Molecular geometry (au)'
do while(.true.)
  read(ifchk,"(a100)",end=101)ctmp
  if(index(ctmp,tag) == 0) cycle

  read(ifchk,"(/)")
  do i=1,NAtm
    read(ifchk,"(a100)",end=101)ctmp
    read(ctmp(10:),*)(XYZ(j,i),j=1,3)
  end do
  ifgeom= .true.
end do
101   if(.not. ifgeom) goto 110

rewind(ifchk)

! read FFx (a.u.)
tag='Cartesian Hessian in GSPHES'
ta2='Molecular Hessian (au)'
1001  read(ifchk,"(a100)",end=1010)ctmp
if(index(ctmp,tag)+index(ctmp,ta2) == 0)goto 1001
if(index(ctmp,ta2) /= 0) IfAna = .true.
read(ifchk,*)

NBlock=(NAtm3-1)/6+1
do i=1,NBlock
  iv1=(i-1)*6+1
  iv2=min(i*6,NAtm3)
  read(ifchk,"(//)",end=1020)
  do j=iv1,NAtm3
    read(ifchk,"(9x,6f12.6)",end=1020)(FFx(k,j),k=iv1,min(j,iv2))
    if(mod(j,3) == 0)read(ifchk,*)
  end do
end do
do i=1,NAtm3
  do j=1,i-1
    FFx(i,j)=FFx(j,i)
  end do
end do

! read APT if ana. Hess.
if(IfAna)then
  ta2='Dipole moment gradient (au)'
  1101  read(ifchk,"(a100)",end=1110)ctmp
  if(index(ctmp,ta2)==0)goto 1101
  read(ifchk,"(//)")
  do i=1,NAtm3
    if(mod(i,3) == 1) read(ifchk,*,end=1110)
    read(ifchk,"(a100)",end=1110)ctmp
    read(ctmp(15:),*,err=1110)(APT(j,i),j=1,3)
  end do
  Infred= 1
end if

! read IZ and mass
tag='Isotopic Masses'
1201  read(ifchk,"(a100)",end=1210)ctmp
if(index(ctmp,tag)==0)goto 1201
read(ifchk,"(/)")
do i=1,NAtm
  read(ifchk,"(a100)",end=1210)ctmp
  read(ctmp(38:),*,err=1210)AMass(i)
  ipo=nonspace(ctmp)
  call ElemZA(0,ctmp(ipo:),ctmp(ipo:),ZA(i))
end do

return
0110  call XError(Intact,"No Cartesian coordinates found!")
1010  call XError(Intact,"No Cartesian force constants found!")
1020  call XError(Intact,"Please check Force Constant matrix!")
1110  call XError(Intact,"Please check APT!")
1210  call XError(Intact,"Please check Isotopic Masses!")
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Read data from FHI-AIMS output
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RdAIMS(ifchk,ihess,iddip,ctmp,Intact,Infred,NAtm,AMass,ZA,XYZ,FFx,APT)
implicit real(kind=8) (a-h,o-z)
parameter(ang2au=1.d0/0.52917720859d0,ev2au=1.d0/(ang2au*ang2au*27.2114d0))
real(kind=8) :: AMass(*),ZA(*),XYZ(3,*),FFx(*),APT(*)
character*100 :: ctmp
logical :: Intact,ifopen

NAtm3 = NAtm * 3
Infred= 0

rewind(ifchk)
! read IZ, mass, and Cartesian coordinates in Angstrom
do i=1,NAtm
  read(ifchk,*,err=110)AMass(i),(XYZ(j,i),j=1,3),ctmp
  call ElemZA(0,ctmp,ctmp,ZA(i))
end do
! Ang --> a.u.
call AScale(NAtm3,ang2au,XYZ,XYZ)

rewind(ihess)
! read FFx (eV/Ang**2)
read(ihess,*,err=210)(FFx(i),i=1,NAtm3*NAtm3)
! eV/Ang**2 --> a.u.
call AScale(NAtm3*NAtm3,ev2au,FFx,FFx)

rewind(iddip)
! read APT (au)
inquire(unit=iddip,opened=ifopen)
if(ifopen)then
  read(iddip,*,err=310)(APT(i),i=1,NAtm3*3)
  Infred= 1
else
  call AClear(NAtm3*3,APT)
end if

return
110   call XError(Intact,"Please check Cartesian coordinates!")
210   call XError(Intact,"Please check Force Constant matrix!")
310   call XError(Intact,"Please check APT matrix!")
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Read data from CP2K output
!
! APT is not available at present.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RdCP2K(ifchk,tag,ctmp,Intact,NAtm,AMass,ZA,XYZ,FFx)
implicit real(kind=8) (a-h,o-z)
parameter(ang2au=1.d0/0.52917720859d0,amu2au=1.660538921d-27/9.10938291d-31)
real(kind=8) :: AMass(*),ZA(*),XYZ(3,*),FFx(NAtm*3,*)
character*100 :: ctmp
character*50 :: tag
logical :: Intact

NAtm3 = NAtm * 3

rewind(ifchk)

! read IZ, mass, and Cartesian coordinates in a.u.
tag=' MODULE QUICKSTEP:  ATOMIC COORDINATES IN angstrom'
101   read(ifchk,"(a100)",end=110)ctmp
if(index(ctmp,tag)==0)goto 101
read(ifchk,"(//)")
do i=1,NAtm
  read(ifchk,"(a100)",end=110)ctmp
  read(ctmp(18:),*)ZA(i),(XYZ(j,i),j=1,3),FFx(1,1),AMass(i)
end do
! Ang --> a.u.
call AScale(NAtm3,ang2au,XYZ,XYZ)

! read m.w. FFx (a.u.) (x 1.d6)
tag=' VIB| Hessian in cartesian coordinates'
201   read(ifchk,"(a100)",end=210)ctmp
if(index(ctmp,tag)==0)goto 201

NBlock=(NAtm3-1)/5+1
do i=1,NBlock
  iv1=(i-1)*5+1
  iv2=min(i*5,NAtm3)
  read(ifchk,"(/)",end=220)
  do j=1,NAtm3
    read(ifchk,"(12x,5f13.6)",end=220)(FFx(k,j),k=iv1,iv2)
  end do
end do

! mass unweighting
fc=amu2au*1.d-6
do i=1,NAtm3
  ia = (i-1)/3+1
  do j=1,NAtm3
    ja = (j-1)/3+1
    FFx(j,i)=FFx(j,i)*sqrt(AMass(ia)*AMass(ja))*fc
  end do
end do

return
110   call XError(Intact,"Please check Cartesian coordinates!")
210   call XError(Intact,"No Cartesian force constants found!")
220   call XError(Intact,"Please check Force Constant matrix!")
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Read data from ADF tape21 or tape13.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RdADF(ifchk,tag,ctmp,Intact,Infred,NAtm,AMass,ZA,XYZ,FFx,APT,DAT1,IDAT,SCR,NWork,SCRFFX)
implicit real(kind=8) (a-h,o-z)
! we assume #atomtypes < NAtm3, which should be safe
real(kind=8) :: AMass(*),ZA(*),XYZ(3,*),FFx(*),APT(*),DAT1(NAtm*3,*),SCR(*),SCRFFX(NWork)
dimension :: IDAT(NAtm*3,*)
character*100 :: ctmp
character*28 :: tag(2)
logical :: Intact,anafc

NAtm3 = NAtm * 3
Infred= 0

! #atomtypes (including dummy atoms)
i=0
tag(1)='Geometry'
tag(2)='nr of atomtypes'
rewind(ifchk)
do while(.true.)
  read(ifchk,"(a100)",end=100,err=100)ctmp
  if(ctmp(1:8) == tag(1)(1:8))then
    read(ifchk,"(a100)",end=100,err=100)ctmp
    if(ctmp(1:15) == tag(2)(1:15))then
      read(ifchk,"(/,a100)",end=100,err=100)ctmp
      read(ctmp,*,end=100,err=100)i
      goto 100
    end if
  end if
end do
100   NTyp=i
if(NTyp < 1)call XError(Intact,"nr of atomtypes < 1!")
if(NTyp > NAtm3)call XError(Intact,"nr of atomtypes > NAtm3!")

! read atomic mass --> DAT1(:,1)
tag(2)='mass'
rewind(ifchk)
do while(.true.)
  read(ifchk,"(a100)",end=210,err=210)ctmp
  if(ctmp(1:8) == tag(1)(1:8))then
    read(ifchk,"(a100)",end=210,err=210)ctmp
    if(ctmp(1:4) == tag(2)(1:4))then
      read(ifchk,*,end=210,err=210)
      read(ifchk,*,end=210,err=210)(DAT1(i,1),i=1,NTyp)
      goto 200
    end if
  end if
end do
200   continue

! read IZ --> DAT1(:,2)
tag(2)='atomtype total charge'
rewind(ifchk)
do while(.true.)
  read(ifchk,"(a100)",end=310,err=310)ctmp
  if(ctmp(1:8) == tag(1)(1:8))then
    read(ifchk,"(a100)",end=310,err=310)ctmp
    if(ctmp(1:21) == tag(2)(1:21))then
      read(ifchk,*,end=310,err=310)
      read(ifchk,*,end=310,err=310)(DAT1(i,2),i=1,NTyp)
      goto 300
    end if
  end if
end do
300   continue

! read atom order index --> IDAT(:,1)
tag(2)='atom order index'
rewind(ifchk)
do while(.true.)
  read(ifchk,"(a100)",end=410,err=410)ctmp
  if(ctmp(1:8) == tag(1)(1:8))then
    read(ifchk,"(a100)",end=410,err=410)ctmp
    if(ctmp(1:16) == tag(2)(1:16))then
      read(ifchk,*,end=410,err=410)NVal
!     Note: only the first half part of IDAT(:,1) will be used
      read(ifchk,*,end=410,err=410)(IDAT(i,1),i=1,NVal)
      goto 400
    end if
  end if
end do
400   continue

! read atomtype index --> IDAT(:,2)
tag(2)='fragment and atomtype index'
rewind(ifchk)
do while(.true.)
  read(ifchk,"(a100)",end=460,err=460)ctmp
  if(ctmp(1:8) == tag(1)(1:8))then
    read(ifchk,"(a100)",end=460,err=460)ctmp
    if(ctmp(1:27) == tag(2)(1:27))then
      read(ifchk,*,end=460,err=460)NVal
!     Note: only the second half part of IDAT(:,2) will be used
      read(ifchk,*,end=460,err=460)(IDAT(i,2),i=1,NVal)
      goto 450
    end if
  end if
end do
450   NAtm0 = NVal/2    ! NAtm0 contains dummy atoms

! Read Cartesian coordinates in the Geometry part (in a.u.)
! Because of some bugs, Freq%xyz may be the initial geometry before the optimization, so Geometry%xyz should be used instead.
ia=0
! tag(1)='Freq'
tag(2)='xyz'
rewind(ifchk)
do while(.true.)
  read(ifchk,"(a100)",end=510,err=510)ctmp
! if(ctmp(1:4) == tag(1)(1:4))then
  if(ctmp(1:8) == tag(1)(1:8))then
    read(ifchk,"(a100)",end=510,err=510)ctmp
    if(ctmp(1:3) == tag(2)(1:3) .and. len_trim(ctmp(1:9)) == 3) then   ! skip Geometry % xyz InputOrder
      read(ifchk,*,end=510,err=510)NVal
      if(NAtm0 * 3 /= NVal) call XError(Intact,"NAtm0*3 /= NVal in Freq % xyz!")
      do i=1,NAtm0
        read(ifchk,"(a100)",end=510,err=510)ctmp
        ityp=IDAT(NAtm0+i,2)
!       skip dummy atoms
        if(DAT1(ityp,1) < 0.8d0 .or. DAT1(ityp,2) < 0.8d0) cycle
        ia=ia+1
        read(ctmp,*,end=100,err=100)(XYZ(j,ia),j=1,3)
        AMass(ia)=DAT1(ityp,1)
        ZA(ia)=DAT1(ityp,2)
      end do
      goto 500
    end if
  end if
end do
500   if(ia /= NAtm) call XError(Intact,"Wrong NAtm in Freq % xyz!")

! read analytic FFx (in a.u.)
anafc = .true.
tag(1)='Hessian'
tag(2)='Analytical Hessian'
rewind(ifchk)
do while(.true.)
  read(ifchk,"(a100)",end=650,err=610)ctmp
  if(ctmp(1:7) == tag(1)(1:7))then
    read(ifchk,"(a100)",end=650,err=610)ctmp
    if(ctmp(1:18) == tag(2)(1:18))then
      read(ifchk,*,end=610,err=610)NVal
      if(NAtm3*NAtm3 /= NVal) call XError(Intact,"NAtm3*NAtm3 /= NVal in Hessian % Analytical Hessian!")
      read(ifchk,*,end=610,err=610)(FFX(i),i=1,NVal)
      goto 680
    end if
  end if
end do
650   anafc = .false.

tag(1)='GeoOpt'
tag(2)='Hessian_CART'
rewind(ifchk)
do while(.true.)
  read(ifchk,"(a100)",end=621,err=621)ctmp
  if(ctmp(1:6) == tag(1)(1:6))then
    read(ifchk,"(a100)",end=621,err=621)ctmp
    if(ctmp(1:12) == tag(2)(1:12))then
      read(ifchk,*,end=621,err=621)NVal
      if(NWork < NVal) call XError(Intact,"NWK in sub. RdData1 is too small!")
      read(ifchk,*,end=621,err=621)(SCRFFX(i),i=1,NVal)
      if(NAtm0 == NAtm)then
!       no dummy atom
        if(NAtm3*NAtm3 /= NVal) call XError(Intact,"NAtm3*NAtm3 /= NVal in GeoOpt % Hessian_CART!")
        call ACopy(NVal,SCRFFX,FFX)
      else
!       remove the blocks of dummy atoms
        if(NAtm0*NAtm0*9 /= NVal) call XError(Intact,"NAtm0*NAtm0*9 /= NVal in GeoOpt % Hessian_CART!")
        call RmDummyFFX(NAtm0,NAtm0*3,IDAT(1,1),IDAT(NAtm0+1,2),DAT1(1,1),DAT1(1,2),SCRFFX,FFX)
      end if
      goto 680
    end if
  end if
end do
680   continue

if(anafc) then
! read analytic APT (in a.u.)
  tag(1)='Hessian'
  tag(2)='Analytical Dipole Derivative'
  rewind(ifchk)
  do while(.true.)
    read(ifchk,"(a100)",end=710,err=710)ctmp
    if(ctmp(1:7) == tag(1)(1:7))then
      read(ifchk,"(a100)",end=710,err=710)ctmp
      if(ctmp(1:28) == tag(2)(1:28))then
        read(ifchk,*,end=710,err=710)NVal
        if(NAtm3*3 /= NVal) call XError(Intact,"NAtm3*3 /= NVal in Hessian % Ana. Dipole Der.!")
        read(ifchk,*,end=710,err=710)(APT(i),i=1,NVal)
        Infred= 1
        goto 750
      end if
    end if
  end do
else
! read numerical APT (in a.u.)
! if symmetry is used with equivalent atoms, some APT elements will be zero!!!
  write(*,"(/,' <<< Note >>>',/, &
     ' For numerical Hessian, if ADF uses symmetry and if there are symmetry-',/, &
     ' equivalent atoms, the normal and adiabatic IR intensities cannot be',/, &
     ' calculated correctly. Please check normal IR intensities first!')")

  tag(1)='Freq'
  tag(2)='Dipole derivatives'
  rewind(ifchk)
  do while(.true.)
    read(ifchk,"(a100)",end=720,err=720)ctmp
    if(ctmp(1:4) == tag(1)(1:4))then
      read(ifchk,"(a100)",end=720,err=720)ctmp
      if(ctmp(1:18) == tag(2)(1:18))then
        read(ifchk,*,end=720,err=720)NVal
        if(NAtm3*3 /= NVal) call XError(Intact,"NAtm3*3 /= NVal in Freq % Dipole derivatives!")
        read(ifchk,*,end=720,err=720)(SCR(i),i=1,NVal)
        goto 700
      end if
    end if
  end do
! actually the APT^T was read above
  700  call Transp(NAtm3,3,SCR,APT)
  Infred= 1
end if
750   continue

return
210   call XError(Intact,"Please check Geometry % mass!")
310   call XError(Intact,"Please check Geometry % atomtype total charge!")
410   call XError(Intact,"Please check Geometry % atom order index!")
460   call XError(Intact,"Please check Geometry % atomtype index!")
510   call XError(Intact,"Please check Freq % xyz!")
610   call XError(Intact,"Please check Hessian % Analytical Hessian!")
621   call XError(Intact,"Please check GeoOpt % Hessian_CART!")
710   call XError(Intact,"Please check Hessian % Ana. Dipole Der.!")
720   call XError(Intact,"Please check Freq % Dipole derivatives!")
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Save a new Hessian matrix without the blocks of dummy atom
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RmDummyFFX(NAtm0,NAtm03,Idx1,Idx2,AtomM,AtomZ,FFX0,FFX)
implicit real(kind=8) (a-h,o-z)
dimension :: Idx1(NAtm0),Idx2(NAtm0)
real(kind=8) :: AtomM(*),AtomZ(*),FFX(*),FFX0(NAtm03,NAtm03)

ip = 0
do i=1,NAtm03
  ia = (i-1)/3+1
  it = Idx2(Idx1(ia))
  if(AtomM(it) < 0.8d0 .or. AtomZ(it) < 0.8d0) cycle
  do j=1,NAtm03
    ja = (j-1)/3+1
    jt = Idx2(Idx1(ja))
    if(AtomM(jt) < 0.8d0 .or. AtomZ(jt) < 0.8d0) cycle
    ip = ip + 1
    FFX(ip) = FFX0(j,i)
  end do
end do

return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Read data from Hyperchem log file
!
! APT is not available at present.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RdHyper(ifchk,tag,ctmp,Intact,NAtm,AMass,ZA,XYZ,FFx,SC1,SC2,SC3,WORK)
implicit real(kind=8) (a-h,o-z)
parameter(One=1.d0,ang2au=One/0.52917720859d0,wn2au=One/5140.48714376d0)
real(kind=8) :: AMass(*),ZA(*),XYZ(3,*),FFx(*),SC1(*),SC2(NAtm*3,*),SC3(*),WORK(*)
character*100 :: ctmp
character*69 :: tag
logical :: Intact

NAtm3 = NAtm * 3

rewind(ifchk)

! read IZ, mass, and Cartesian coordinates in Ang.
tag='Atom  Z     Charge            Coordinates(Angstrom)              Mass'
101   read(ifchk,"(a100)",end=110)ctmp
if(index(ctmp,tag)==0)goto 101
read(ifchk,*)
do i=1,NAtm
  read(ifchk,"(a100)",end=110)ctmp
  read(ctmp(4:),*)ZA(i),SC1(1),(XYZ(j,i),j=1,3),AMass(i)
end do
! Ang --> a.u.
call AScale(NAtm3,ang2au,XYZ,XYZ)

! read NAtm3 frequencies and normal modes --> SC1 and SC2
tag='Frequencies of the Normal Modes and Transformation Matrix.'
201   read(ifchk,"(a100)",end=210)ctmp
if(index(ctmp,tag(1:58))==0)goto 201

tag='    mu '
NBlock=(NAtm3-1)/6+1
do i=1,NBlock
  do while(.true.)
    read(ifchk,"(a100)",end=220,err=220)ctmp
    if(index(ctmp,tag(1:7)) /= 0) exit
  end do
  iv1=(i-1)*6+1
  iv2=min(i*6,NAtm3)
  read(ctmp(9:),*,end=220,err=220)(SC1(j),j=iv1,iv2)
  read(ifchk,*,end=220,err=220)
  do j=1,NAtm3
    read(ifchk,"(a100)",end=220,err=220)ctmp
    read(ctmp(9:),*,end=220,err=220)(SC2(j,k),k=iv1,iv2)
  end do
end do
! cm-1 --> a.u.
call AScale(NAtm3,wn2au,SC1,SC1)
! normal modes: mass-unweighting and renormalization
do i=1,NAtm3
  do j=1,NAtm3
    ja=(j-1)/3+1
    SC2(j,i)=SC2(j,i)/sqrt(AMass(ja))
  end do
  X = dotx(NAtm3,SC2(1,i),SC2(1,i))
  X = One/sqrt(X)
  call AScale(NAtm3,X,SC2(1,i),SC2(1,i))
end do

! calculate force constant matrix --> FFX
call Frq2FFX(NAtm3,NAtm3,AMass,SC1,SC2,FFX,SC3,WORK,SC1(NAtm3+1))

return
110   call XError(Intact,"Please check Cartesian coordinates!")
210   call XError(Intact,"No frequency results found!")
220   call XError(Intact,"Please check frequencies and normal modes!")
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Read data from Jaguar output
!
! APT is not available at present.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RdJaguar(ifchk,tag,ctmp,Intact,NAtm,AMass,ZA,XYZ,FFx,SC1,SC2,SC3,WORK)
implicit real(kind=8) (a-h,o-z)
parameter(One=1.d0,ang2au=One/0.52917720859d0,wn2au=One/5140.48714376d0)
real(kind=8) :: AMass(*),ZA(*),XYZ(3,*),FFx(*),SC1(*),SC2(NAtm*3,*),SC3(*),WORK(*)
character*100 :: ctmp
character*53 :: tag
logical :: Intact

NAtm3 = NAtm * 3
NCol = 6
! For old version!!!
! NCol = 7

! read IZ and Cartesian coordinates in Ang.
tag='final geometry:'
rewind(ifchk)
do while(.true.)
  read(ifchk,"(a100)",end=20)ctmp
  if(index(ctmp,tag(1:15)) /= 0) goto 40
end do

20    tag='Input geometry:'
rewind(ifchk)
do while(.true.)
  read(ifchk,"(a100)")ctmp
  if(index(ctmp,tag(1:15)) /= 0) goto 40
end do

40    read(ifchk,"(/)",end=110,err=110)
do i=1,NAtm
  call CClear(15,tag)
  read(ifchk,"(a100)",end=110,err=110)ctmp
  read(ctmp,*,end=110,err=110)tag,(XYZ(j,i),j=1,3)
  call rmnumb(15,tag)
  call CClear(15,ctmp)
  read(tag,*,end=110,err=110)ctmp
  call ElemZA(0,ctmp,ctmp,ZA(i))
end do
! Ang --> a.u.
call AScale(NAtm3,ang2au,XYZ,XYZ)

! read atomic masses
! if masses are not printed, the most abundant isotopic masses are assumed,
! But if this is not true, the re-calculated force constant matrix from normal
! modes (and then the frequencies, normal modes, and so on) will be wrong.
tag='vdw   vdw2    cov     mass  grid  daf  charge   basis'
rewind(ifchk)
AMass(1)=-One
do while(.true.)
  read(ifchk,"(a100)",end=300,err=300)ctmp
  if(index(ctmp,tag(1:53)) /= 0) exit
end do
ibgn=index(ctmp,"   mass")
do i=1,NAtm
  read(ifchk,"(a100)",end=210,err=210)ctmp
  read(ctmp(ibgn:),*,end=210,err=210)AMass(i)
end do

! the most abundant isotopic masses are assumed
300   if(AMass(1) < 0.d0) call MasLib(0,NAtm,AMass,ZA)

! generate rot. + trans. modes
! center of mass ---> WORK; XYZ in CMCS --> Sc1
call MassCent(NAtm,AMass,XYZ,SC1,WORK)
! principal moment of inertia; Eigenvector --> Sc3(1:9)
call MIner(.False.,NAtm,AMass,Sc1,Sc3,Sc3(10),WORK)
! generate m.w. vectors of translations and rotations --> WORK
call TRVec(.False.,NAtm,NTR,Imiss,AMass,Sc1,WORK,Sc3,Sc2)
NVib = NAtm3 - NTR
! save NTR rot. & trans. modes to SC2(:,Nvib+1:)
call ACopy(NAtm3*NTR,WORK,SC2(1,NVib+1))

! read NVib frequencies and normal modes --> SC1 and SC2
call CClear(53,tag)
tag='normal modes in cartesian coordinates:'
rewind(ifchk)
do while(.true.)
  read(ifchk,"(a100)",end=310,err=310)ctmp
  if(index(ctmp,tag(1:38)) /= 0) exit
end do

call CClear(53,tag)
tag='frequencies'
NBlock=(NVib-1)/NCol+1
do i=1,NBlock
  iv1=(i-1)*NCol+1
  iv2=min(i*NCol,NVib)
  do while(.true.)
    read(ifchk,"(a100)",end=320,err=320)ctmp
    ibgn = index(ctmp,tag(1:11))
    if(ibgn /= 0) then
      read(ctmp(ibgn+11:),*,end=320,err=320)(SC1(j),j=iv1,iv2)
    else if(index(ctmp,"X") /= 0) then
      exit
    end if
  end do
  read(ctmp(15:),*,end=320,err=320)(SC2(1,k),k=iv1,iv2)
  do j=2,NAtm3
    read(ifchk,"(a100)",end=320,err=320)ctmp
    read(ctmp(15:),*,end=320,err=320)(SC2(j,k),k=iv1,iv2)
  end do
end do
! cm-1 --> a.u.
call AScale(NVib,wn2au,SC1,SC1)
call AClear(NAtm3-NVib,SC1(NVib+1))

! normal modes: mass-unweighting (R.+T. only) and renormalization
do i=1,NAtm3
  if(i > NVib)then
    do j=1,NAtm3
      ja=(j-1)/3+1
      SC2(j,i)=SC2(j,i)/sqrt(AMass(ja))
    end do
  end if
  X = dotx(NAtm3,SC2(1,i),SC2(1,i))
  X = One/sqrt(X)
  call AScale(NAtm3,X,SC2(1,i),SC2(1,i))
end do

! calculate force constant matrix --> FFX
call Frq2FFX(NAtm3,NAtm3,AMass,SC1,SC2,FFX,SC3,WORK,SC1(NAtm3+1))

return
110   call XError(Intact,"Please check Cartesian coordinates!")
210   call XError(Intact,"Please check atomic masses!")
310   call XError(Intact,"No frequency results found!")
320   call XError(Intact,"Please check frequencies and normal modes, or try NCol = 7!")
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Read data from MOLDEN file. APT is not available.
!
! The normal mode can be mass-unweighted or mass-weighted, which will be detected by the program.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RdMOLDEN(ifchk,tag,ctmp,Intact,NAtm,AMass,ZA,XYZ,FFx,SC1,SC2,SC3,WORK)
implicit real(kind=8) (a-h,o-z)
parameter(One=1.d0,wn2au=One/5140.48714376d0)
real(kind=8) :: AMass(*),ZA(*),XYZ(3,*),FFx(*),SC1(*),SC2(NAtm*3,*),SC3(*),WORK(*)
character*100 :: ctmp
character*15 :: tag
logical :: Intact,full3n,DoUnWt

NAtm3 = NAtm * 3

! read IZ and Cartesian coordinates in a.u.
tag='[FR-COORD]'
rewind(ifchk)
do while(.true.)
  read(ifchk,"(a100)",end=060,err=060)ctmp
  call charl2u(ctmp)
  if(index(ctmp,tag(1:10)) /= 0) exit
end do

do i=1,NAtm
  read(ifchk,"(a100)",end=160,err=160)ctmp
  read(ctmp,*,end=160,err=160)tag,(XYZ(j,i),j=1,3)
  call rmnumb(15,tag)
  call CClear(15,ctmp)
  read(tag,*,end=160,err=160)ctmp
  call ElemZA(0,ctmp,ctmp,ZA(i))
end do

! atomic masses: the most abundant isotopic masses are assumed, and they must have been used in your frequency calculation, otherwise the
! calculated FFX below will be wrong!
call MasLib(0,NAtm,AMass,ZA)

! read NVib or NAtm3 frequencies in cm-1 --> SC3(1:NAtm3)
call CClear(15,tag)
tag='[FREQ]'
rewind(ifchk)
do while(.true.)
  read(ifchk,"(a100)",end=260,err=260)ctmp
  call charl2u(ctmp)
  if(index(ctmp,tag(1:6)) /= 0) exit
end do

NFrq=0
do i=1,NAtm3
  read(ifchk,"(a100)",end=360,err=360)ctmp
  if(len_trim(ctmp) == 0) exit
  if(index(ctmp,"[") /= 0) exit
  if(index(ctmp,"]") /= 0) exit
  NFrq=NFrq+1
  read(ctmp,*,end=360,err=360)SC3(NFrq)
end do
! cm-1 --> a.u.
call AScale(NFrq,wn2au,SC3,SC3)
full3n = .true.
if(NFrq < NAtm3) full3n = .false.

! generate rot. + trans. modes
if(.not. full3n)then
  ! center of mass ---> WORK; XYZ in CMCS --> Sc1
  call MassCent(NAtm,AMass,XYZ,SC1,WORK)
  ! principal moment of inertia; Eigenvector --> Sc3(NAtm3+1:NAtm3+9)
  call MIner(.False.,NAtm,AMass,Sc1,Sc3(NAtm3+1),Sc3(NAtm3+10),WORK)
  ! generate m.w. vectors of translations and rotations --> WORK
  call TRVec(.False.,NAtm,NTR,Imiss,AMass,Sc1,WORK,Sc3(NAtm3+1),Sc2)
  if(NFrq + NTR /= NAtm3) call XError(Intact,"NFrq + NTR /= NAtm3!")
  ! save NTR rot. & trans. modes to SC2(:,NFrq+1:)
  call ACopy(NAtm3*NTR,WORK,SC2(1,NFrq+1))
end if

! read NVib or NAtm3 normal modes --> SC2
! reordered frequencies --> SC1
tag='[FR-NORM-COORD]'
rewind(ifchk)
do while(.true.)
  read(ifchk,"(a100)",end=460,err=460)ctmp
  call charl2u(ctmp)
  if(index(ctmp,tag(1:15)) /= 0) exit
end do

do i=1,NFrq
  read(ifchk,"(a100)",end=560,err=560)ctmp
  read(ctmp,*,end=560,err=560)tag,ifrq
  if(ifrq < 1 .or. ifrq > NFrq)then
    write(*,"('  NFrq = ',i6,'  IFrq = ',i6)")NFrq,ifrq
    call XError(Intact,"IFrq in [FR-NORM-COORD] is out of range!")
  end if
  SC1(i) = SC3(ifrq)
  read(ifchk,*,end=570,err=570)(SC2(j,i),j=1,NAtm3)
end do
! trans. and rot. frequencies are zero
if(.not. full3n) call AClear(NAtm3-NFrq,SC1(NFrq+1))

! normal modes: do mass-unweighting (if necessary) and renormalization of R.+T. and V.
DoUnWt = .false.
! Check whether the normal modes have been mass weighted:
! calculate U = L^T * L; if L is mass weighted, U = I.
call Transp(NAtm3,NFrq,SC2,SC3)
call MMpyMF(NFrq,NAtm3,NFrq,SC3,SC2,WORK)
call MSubI(NFrq,WORK,WORK)
Amax=ArMax(NFrq*NFrq,i,WORK)
Amin=ArMin(NFrq*NFrq,i,WORK)
Amax=max(abs(Amax),abs(Amin))
if(Amax < 1.d-4) DoUnWt = .true.

do i=1,NAtm3
  if(i > NFrq .or. DoUnWt)then
    do j=1,NAtm3
      ja=(j-1)/3+1
      SC2(j,i)=SC2(j,i)/sqrt(AMass(ja))
    end do
  end if
  X = dotx(NAtm3,SC2(1,i),SC2(1,i))
  X = One/sqrt(X)
  call AScale(NAtm3,X,SC2(1,i),SC2(1,i))
end do

! calculate force constant matrix --> FFX
call Frq2FFX(NAtm3,NAtm3,AMass,SC1,SC2,FFX,SC3,WORK,SC1(NAtm3+1))

return
060   call XError(Intact,"[FR-COORD] is not found!")
160   call XError(Intact,"Please check the [FR-COORD] section!")
260   call XError(Intact,"[FREQ] is not found!")
360   call XError(Intact,"Please check the [FREQ] section!")
460   call XError(Intact,"[FR-NORM-COORD] is not found!")
560   call XError(Intact,"Please check the [FR-NORM-COORD] section!")
570   call XError(Intact,"Please check your normal modes!")
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Read data from Crystal output
!
! APT is not available at present.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RdCry(ifchk,tag,ctmp,Intact,NAtm,AMass,ZA,XYZ,FFx,SC1,SC2,SC3,WORK)
implicit real(kind=8) (a-h,o-z)
parameter(One=1.d0,ang2au=One/0.52917720859d0,wn2au=One/5140.48714376d0)
real(kind=8) :: AMass(*),ZA(*),XYZ(3,*),FFx(*),SC1(*),SC2(NAtm*3,*),SC3(*),WORK(*)
character*100 :: ctmp
character*60 :: tag
logical :: Intact

NAtm3 = NAtm * 3
NCol = 6

! read IZ and Cartesian coordinates in Ang.
tag='MOLECULAR CALCULATION'
rewind(ifchk)
do while(.true.)
  read(ifchk,"(a100)",end=1010,err=1010)ctmp
  if(index(ctmp,tag(1:21)) /= 0) exit
end do

tag='FRQFRQFRQFRQFRQFRQFRQFRQFRQFRQFRQFRQFRQFRQFRQFRQFRQFRQFRQFRQ'
do while(.true.)
  read(ifchk,"(a100)",end=1020,err=1020)ctmp
  if(index(ctmp,tag(1:60)) /= 0) exit
end do

tag='- ATOMS IN THE UNIT CELL:'
do while(.true.)
  read(ifchk,"(a100)",end=1030,err=1030)ctmp
  if(index(ctmp,tag(1:25)) /= 0) exit
end do

tag='X(ANGSTROM)         Y(ANGSTROM)         Z(ANGSTROM)'
read(ifchk,"(a100)",end=1030,err=1030)ctmp
if(index(ctmp,tag(1:51)) == 0) goto 1010

read(ifchk,*,end=1040,err=1040)
do i=1,NAtm
  read(ifchk,"(a100)",end=1040,err=1040)ctmp
  read(ctmp(8:100),*,end=1040,err=1040)IZ,tag,(XYZ(j,i),j=1,3)
  ZA(i) = dble(IZ)
end do
! Ang --> a.u.
call AScale(NAtm3,ang2au,XYZ,XYZ)

! read atomic masses
tag='ATOMS ISOTOPIC MASS (AMU) FOR FREQUENCY CALCULATION'
do while(.true.)
  read(ifchk,"(a100)",end=1050,err=1050)ctmp
  if(index(ctmp,tag(1:51)) /= 0) exit
end do
read(ifchk,*,end=1050,err=1050)
read(ifchk,"(4(10x,f10.4))",end=1050,err=1050)(AMass(i),i=1,NAtm)

! read NAtm3 trans.+rot.+vib. frequencies and normal modes --> SC1 and SC2
tag='MODES         EIGV          FREQUENCIES     IRREP  IR'
do while(.true.)
  read(ifchk,"(a100)",end=1060,err=1060)ctmp
  if(index(ctmp,tag(1:53)) /= 0) exit
end do
read(ifchk,*,end=1060,err=1060)
do i=1,NAtm3
  read(ifchk,"(a100)",end=1065,err=1065)ctmp
  ctmp(6:6)=" "
  call CClear(14,ctmp(11:))
  read(ctmp(1:100),*)i1,i2,X
  SC1(i1:i2)=X
  if(i2 >= NAtm3) exit
end do

! cm-1 --> a.u.
call AScale(NAtm3,wn2au,SC1,SC1)

tag='NORMAL MODES NORMALIZED TO CLASSICAL AMPLITUDES'
do while(.true.)
  read(ifchk,"(a100)",end=1070,err=1070)ctmp
  if(index(ctmp,tag(1:47)) /= 0) exit
end do
NBlock=(NAtm3-1)/NCol+1
do i=1,NBlock
  iv1=(i-1)*NCol+1
  iv2=min(i*NCol,NAtm3)
  read(ifchk,"(//)",end=1075,err=1075)
  do j=1,NAtm3
    read(ifchk,"(a100)",end=1075,err=1075)ctmp
    read(ctmp(15:100),*,end=1075,err=1075)(SC2(j,k),k=iv1,iv2)
  end do
end do

! calculate force constant matrix --> FFX
call Frq2FFX(NAtm3,NAtm3,AMass,SC1,SC2,FFX,SC3,WORK,SC1(NAtm3+1))

return
1010  call XError(Intact,"This is not a molecular calculation.")
1020  call XError(Intact,"This is not a frequency calculation.")
1030  call XError(Intact,"Failed to read #Atoms.")
1040  call XError(Intact,"Please check Cartesian coordinates!")
1050  call XError(Intact,"Please check atomic masses!")
1060  call XError(Intact,"No frequencies found!")
1065  call XError(Intact,"Failed to read frequencies.")
1070  call XError(Intact,"No normal modes found!")
1075  call XError(Intact,"Failed to read normal modes.")
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Read data from Spartan *.smol
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RdSptn(ifchk,tag,ctmp,Intact,Infred,NAtm,ZA,XYZ,FFx,APT,SCR)
implicit real(kind=8) (a-h,o-z)
real(kind=8) :: ZA(*),XYZ(3,*),FFx(*),APT(*),SCR(*)
character*100 :: ctmp
character*14 :: tag
logical :: Intact

NAtm3 = NAtm * 3
NTT = NAtm3*(NAtm3+1)/2
Infred= 0

! read FFx (a.u.)
tag='HESSIAN'
rewind(ifchk)
do while(.true.)
  read(ifchk,"(a100)",end=1010,err=1010)ctmp
  Idx = index(ctmp,tag(1:7))
  if(Idx /= 0) then
   	! the section of HESSIAN options (length = 0) should be skipped
    length = len_trim(ctmp(Idx+7:20))
    if(length > 0) goto 100
  end if
end do
100   read(ifchk,*,end=1020,err=1020)
read(ifchk,*,end=1020,err=1020)(Scr(i),i=1,NTT)
call LT2Sqr(NAtm3,Scr,FFx)

! read nuclear charges, and Cartesian coordinates in a.u.
tag='GEOMETRY2'
do while(.true.)
  read(ifchk,"(a100)",end=1030,err=1030)ctmp
  if(index(ctmp,tag(1:9)) /= 0) exit
end do
do i=1,NAtm
  read(ifchk,"(a100)",end=1040,err=1040)ctmp
  read(ctmp,*,end=1040,err=1040)IZ,(XYZ(j,i),j=1,3)
  ZA(i) = dble(IZ)
end do

! read APT (a.u.); optional
tag='VALUE DIPDERIV'
call AClear(NAtm3*3,APT)
do while(.true.)
  read(ifchk,"(a100)",end=300,err=300)ctmp
  if(index(ctmp,tag(1:14)) /= 0) exit
end do
read(ifchk,*,end=300,err=300)(APT(i),i=1,3*NAtm3)
Infred= 1

300   continue

return
1010  call XError(Intact,"Failed to read BEGINARCHIVE.")
1020  call XError(Intact,"Failed to read HESSIAN.")
1030  call XError(Intact,"Failed to read GEOMETRY2.")
1040  call XError(Intact,"Please check Cartesian coordinates!")
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Read data from PSI log file
!
! APT is not available at present.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RdPSI(ifchk,tag,ctmp,Intact,NAtm,AMass,ZA,XYZ,FFx)
implicit real(kind=8) (a-h,o-z)
real(kind=8) :: AMass(*),ZA(*),XYZ(3,*),FFx(NAtm*3,*)
character*100 :: ctmp
character*60 :: tag
character*3 :: Elem
logical :: Intact

NAtm3 = NAtm * 3

rewind(ifchk)

! read Cartesian coordinates in a.u.
tag='## Reference geometry (Symmetry 0) ##'
do while(.true.)
  read(ifchk,"(a100)",end=1010,err=1010)ctmp
  if(index(ctmp,tag(1:37)) /= 0) exit
end do
read(ifchk,"(///)",end=1020,err=1020)
do i=1,NAtm
  read(ifchk,"(a100)",end=1020,err=1020)ctmp
  read(ctmp,*,end=1020,err=1020)IX,(XYZ(j,i),j=1,3)
end do

! read ZA
tag='Center              X                  Y                   Z'
do while(.true.)
  read(ifchk,"(a100)",end=1030,err=1030)ctmp
  if(index(ctmp,tag(1:60)) /= 0) exit
end do
read(ifchk,*,end=1040,err=1040)
do i=1,NAtm
  read(ifchk,"(a100)",end=1040,err=1040)ctmp
  read(ctmp,*,end=1040,err=1040)Elem
  call ElemZA(0,Elem,Elem,ZA(i))
end do

! read FFx (a.u.)
tag='Force Constants in cartesian coordinates.'
do while(.true.)
  read(ifchk,"(a100)",end=1050,err=1050)ctmp
  if(index(ctmp,tag(1:41)) /= 0) exit
end do
NBlock=(NAtm3-1)/5+1
do i=1,NBlock
  iv1=i*5-4
  iv2=min(i*5,NAtm3)
  read(ifchk,"(//)",end=1060,err=1060)
  do j=1,NAtm3
    read(ifchk,*,end=1060,err=1060)IX,(FFx(k,j),k=iv1,iv2)
  end do
end do

! read mass
tag='Nuclear masses:'
do while(.true.)
  read(ifchk,"(a100)",end=1070,err=1070)ctmp
  if(index(ctmp,tag(1:15)) /= 0) exit
end do
read(ifchk,*,end=1080,err=1080)(AMass(i),i=1,NAtm)

return
1010  call XError(Intact,"No Cartesian coordinates found!")
1020  call XError(Intact,"Please check Cartesian coordinates!")
1030  call XError(Intact,"No element symbols found!")
1040  call XError(Intact,"Failed to read element symbols.")
1050  call XError(Intact,"No Hessian found!")
1060  call XError(Intact,"Failed to read Hessian.")
1070  call XError(Intact,"No atomic masses found!")
1080  call XError(Intact,"Failed to read atomic masses.")
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Read data from DMol3 output
!
! APT is not available at present.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RdDMol(ifchk,tag,ctmp,Intact,NAtm,AMass,ZA,XYZ,FFx,SC1,SC2,SC3,WORK)
implicit real(kind=8) (a-h,o-z)
real(kind=8) :: AMass(*),ZA(*),XYZ(3,*),FFx(*),SC1(*),SC2(NAtm*3,*),SC3(*),WORK(*)
character*200 :: ctmp
character*36 :: tag
character*3 :: Elem
logical :: Intact,found

NAtm3 = NAtm * 3
NCol = 9

! read ZA and Cartesian coordinates in au.
! the last geometry before Vibrations Section should be used
tag='++ Entering Vibrations Section ++'
found=.false.
rewind(ifchk)
do while(.true.)
  read(ifchk,"(a100)",end=1010,err=1010)ctmp
  if(index(ctmp,tag(1:33)) /= 0 .and. found) exit
  if(index(ctmp,"$coordinates") /= 0) then
    found=.true.
    do i=1,NAtm
      read(ifchk,*,err=1020,end=1020)Elem,XYZ(:,i)
      call ElemZA(0,Elem,Elem,ZA(i))
    end do
  end if
end do

! read NAtm3 or Nvib frequencies (a.u.) --> SC1
tag='vibrational frequencies, intensities'
do while(.true.)
  read(ifchk,"(a100)",end=1030,err=1030)ctmp
  if(index(ctmp,tag(1:36)) /= 0) exit
end do
read(ifchk,*,end=1040,err=1040)
do i=1,NAtm3
  read(ifchk,"(a100)",end=1040,err=1040)ctmp
  read(ctmp(1:100),*,end=1040,err=1040)Nvib,SC1(i)
  if(Nvib == NAtm3) exit
end do
Nvib = i

! read normal modes --> SC2
tag='Frequencies (cm-1) and normal modes'
do while(.true.)
  read(ifchk,"(a100)",end=1050,err=1050)ctmp
  if(index(ctmp,tag(1:36)) /= 0) exit
end do
NBlock=(Nvib-1)/NCol+1
do i=1,NBlock
  iv1=(i-1)*NCol+1
  iv2=min(i*NCol,Nvib)
  read(ifchk,"(/)",end=1060,err=1060)
  do j=1,NAtm3
    read(ifchk,"(a200)",end=1060,err=1060)ctmp
    read(ctmp(6:200),*,end=1060,err=1060)(SC2(j,k),k=iv1,iv2)
  end do
  read(ifchk,"(/)",end=1060,err=1060)
end do
! trans. and rot. mode elements are zero
if(Nvib < NAtm3) call AClear((NAtm3-Nvib)*NAtm3,SC2(1,Nvib+1))

! read atomic masses
tag='   Zero point vibrational energy: '
do while(.true.)
  read(ifchk,"(a100)",end=1080,err=1080)ctmp
  if(index(ctmp,tag(1:34)) /= 0) exit
end do
tag='Has Mass'
read(ifchk,*,end=1090,err=1090)
do i=1,NAtm
  read(ifchk,"(a100)",end=1090,err=1090)ctmp
  j=index(ctmp,tag(1:8))
  if(j == 0) goto 1090
  read(ctmp(j+8:),*,end=1090,err=1090)AMass(i)
end do

! normal modes: do mass-unweighting and renormalization
do i=1,Nvib
  do j=1,NAtm3
    ja=(j-1)/3+1
    SC2(j,i)=SC2(j,i)/sqrt(AMass(ja))
  end do
  X = dotx(NAtm3,SC2(1,i),SC2(1,i))
  X = 1.d0/sqrt(X)
  call AScale(NAtm3,X,SC2(1,i),SC2(1,i))
end do

! calculate force constant matrix --> FFX
call Frq2FFX(NAtm3,Nvib,AMass,SC1,SC2,FFX,SC3,WORK,SC1(NAtm3+1))

return
1010  call XError(Intact,"No Cartesian coordinates found.")
1020  call XError(Intact,"Failed to read Cartesian coordinates.")
1030  call XError(Intact,"No frequencies found!")
1040  call XError(Intact,"Failed to read frequencies.")
1050  call XError(Intact,"No normal modes found!")
1060  call XError(Intact,"Failed to read normal modes.")
1080  call XError(Intact,"No atomic mass section found.")
1090  call XError(Intact,"Failed to read atomic masses.")
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Read data from ACES output
! APT is not available for ACES3 and numerical freq by ACES2.
! For analytical freq by ACES2, APT is printed, but the atoms are reordered. So APT will not be read.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RdACES(ifchk,tag,ctmp,Intact,NAtm,AMass,ZA,XYZ,FFx,SC1,SC2,SC3,WORK)
implicit real(kind=8) (a-h,o-z)
parameter(One=1.d0,ang2au=One/0.52917720859d0,wn2au=One/5140.48714376d0)
real(kind=8) :: AMass(*),ZA(*),XYZ(3,*),FFx(*),SC1(*),SC2(3,NAtm,*),SC3(*),WORK(*)
character*200 :: ctmp
character*98 :: tag
logical :: Intact

NAtm3 = NAtm * 3
NCol = 3

! read ZA and Cartesian coordinates (in Ang.)
tag='Symbol    Number           X              Y              Z'
rewind(ifchk)
do while(.true.)
  read(ifchk,"(a100)",end=1010,err=1010)ctmp
  if(index(ctmp,tag(1:47)) /= 0) exit
end do
read(ifchk,*,end=1020,err=1020)
do i=1,NAtm
  read(ifchk,"(a100)",end=1020,err=1020)ctmp
  read(ctmp(13:100),*,err=1020,end=1020)IZ,XYZ(:,i)
  ZA(i)=dble(IZ)
end do
! ang --> a.u.
call AScale(NAtm3,ang2au,XYZ,XYZ)
! the most abundant isotopic masses are assumed
call MasLib(0,NAtm,AMass,ZA)

! read Nvib frequencies (cm-1) --> SC1
tag='Normal Coordinate Analysis'
do while(.true.)
  read(ifchk,"(a100)",end=1040,err=1040)ctmp
  if(index(ctmp,tag(1:26)) /= 0) exit
end do
tag='VIBRATION'
read(ifchk,"(5/)",end=1050,err=1050)
Nvib = 0
do i=1,NAtm3
  read(ifchk,"(a100)",end=1050,err=1050)ctmp
  if(index(ctmp,tag(1:9)) /= 0) then
    Nvib = Nvib + 1
    if(ctmp(25:25) == "i") ctmp(25:25)= " "
    if(ctmp(34:34) == "i") ctmp(34:34)= " "
    read(ctmp(12:100),*,end=1050,err=1050)SC1(Nvib)
  end if
end do
! cm-1 --> a.u.
call AScale(NAtm3,wn2au,SC1,SC1)

! read m.w. normal modes --> SC2
! This part doesn't work for the first column of normal modes in CFour output.
tag='                                   Normal Coordinates'
do while(.true.)
  read(ifchk,"(a100)",end=1060,err=1060)ctmp
  if(index(ctmp,tag(1:53)) /= 0) exit
end do
NBlock=(Nvib-1)/NCol+1
tag='VIBRATIONX       Y       Z'
do i=1,NBlock
  iv1=(i-1)*NCol+1
  iv2=min(i*NCol,Nvib)
  do while(.true.)
    read(ifchk,"(a100)",end=1070,err=1070)ctmp
    if(index(ctmp,tag(1:9)) /= 0) exit
  end do
  do j=1,NAtm
    read(ifchk,"(a100)",end=1070,err=1070)ctmp
    if(index(ctmp,tag(10:26)) /= 0) read(ifchk,"(a100)",end=1070,err=1070)ctmp
    read(ctmp(5:100),*,end=1070,err=1070) ((SC2(ix,j,k),ix=1,3),k=iv1,iv2)
  end do
end do
! trans. and rot. mode elements are zero
if(Nvib < NAtm3)call AClear((NAtm3-Nvib)*NAtm3,SC2(1,1,Nvib+1))

! normal modes: do mass-unweighting and renormalization
do i=1,Nvib
  do j=1,NAtm
    do k=1,3
      SC2(k,j,i)=SC2(k,j,i)/sqrt(AMass(j))
    end do
  end do
  X = dotx(NAtm3,SC2(1,1,i),SC2(1,1,i))
  X = 1.d0/sqrt(X)
  call AScale(NAtm3,X,SC2(1,1,i),SC2(1,1,i))
end do

! calculate force constant matrix --> FFX
call Frq2FFX(NAtm3,Nvib,AMass,SC1,SC2,FFX,SC3,WORK,SC1(NAtm3+1))

return
1010  call XError(Intact,"No Cartesian coordinates found.")
1020  call XError(Intact,"Failed to read Cartesian coordinates.")
1040  call XError(Intact,"No frequencies found!")
1050  call XError(Intact,"Failed to read frequencies.")
1060  call XError(Intact,"No normal modes found!")
1070  call XError(Intact,"Failed to read normal modes.")
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Save UniMoVib (ALM) data file (ialm) and/or localmode.dat (part I).
!
! Format of ALM file:
!
! One text line
! NAtm
! AMass
! ZA
! XYZ (in a.u.)
! FFX (square matrix in a.u.)
! APT (in a.u.)
! DPR(6,NAtm3) (in a.u.)
!
! If the frequencies are corrected by the experimental ones (saved in Freq(:,5)), the force constant matrix will be recalculated.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine SavALM(ifalm,ifloc,ialm,iloc,Infred,IRaman,NAtm,NAtm3,NVib,AMass,ZA,XYZ,FFX,APT,DPol,IExpt,AL,Freq,SC1,SC2,WORK,EIG)
implicit real(kind=8) (a-h,o-z)
logical :: ifalm,ifloc
real(kind=8) :: AMass(*),ZA(*),XYZ(*),FFX(*),APT(*),DPol(*),AL(*),Freq(NAtm3,*), SC1(*),SC2(*),WORK(*),EIG(*)

NAtm9 = NAtm * 9
NSS = NAtm3 * NAtm3

if(ifalm) then
  rewind(ialm)

  write(ialm,"(' UniMoVib DATA FILE (THIS TITLE CAN BE MODIFIED)')")

  write(ialm,"('NATM')")
  write(ialm,"(i5)")NAtm

  write(ialm,"('AMASS')")
  write(ialm,"(5d20.10)")(AMass(i),i=1,NAtm)

  write(ialm,"('ZA')")
  write(ialm,"(5d20.10)")(ZA(i),i=1,NAtm)

  write(ialm,"('XYZ')")
  write(ialm,"(5d20.10)")(XYZ(i),i=1,NAtm3)

  write(ialm,"('FFX')")
end if

if(ifloc) then
  rewind(iloc)

  write(iloc,"(' LOCALMODE DATA FILE (THIS TITLE CAN BE MODIFIED)')")

  write(iloc,"(' $CONTRL NAtm=',i5.5,' NVib=',i5.5,' $END')")NAtm,NVib

  write(iloc,"(' $AMASS  $END')")
  write(iloc,"(5d20.10)")(AMass(i),i=1,NAtm)

  write(iloc,"(' $ZA  $END')")
  write(iloc,"(5d20.10)")(ZA(i),i=1,NAtm)

  write(iloc,"(' $XYZ  $END')")
  write(iloc,"(5d20.10)")(XYZ(i),i=1,NAtm3)

  write(iloc,"(' $FFX  $END')")
end if

if(IExpt == 0)then
  if(ifalm) write(ialm,"(5d20.10)")(FFX(i),i=1,NSS)
  if(ifloc) write(iloc,"(5d20.10)")(FFX(i),i=1,NSS)
else
  ! calculate experimentally corrected force constant matrix --> SC1
  call Frq2FFX(NAtm3,NVib,AMass,Freq(1,5),AL,SC1,SC2,WORK,EIG)
  if(ifalm) write(ialm,"(5d20.10)")(SC1(i),i=1,NSS)
  if(ifloc) write(iloc,"(5d20.10)")(SC1(i),i=1,NSS)
end if

if(ifalm) then

  if(Infred == 0) then
    write(ialm,"('NOAPT')")
  else
    write(ialm,"('APT')")
    write(ialm,"(5d20.10)")(APT(i),i=1,NAtm9)
  end if

  if(IRaman == 0) then
    write(ialm,"('NODPR')")
  else
    write(ialm,"('DPR')")
    write(ialm,"(5d20.10)")(DPol(i),i=1,6*NAtm3)
  end if

  write(ialm,"(/)")
end if

if(ifloc) then
! vib + trans & rot
  write(iloc,"(' $NMMODE  $END')")
  write(iloc,"(5d20.10)")(AL(i),i=1,NSS)

  if(Infred == 0) then
    write(iloc,"(' $APT NOAPT=.TRUE. $END')")
  else
    write(iloc,"(' $APT NOAPT=.FALSE. $END')")
    write(iloc,"(5d20.10)")(APT(i),i=1,NAtm9)
  end if

  if(IRaman == 0) then
    write(iloc,"(' $DPR NODPR=.TRUE. $END')")
  else
    write(iloc,"(' $DPR NODPR=.FALSE. $END')")
    write(iloc,"(5d20.10)")(DPol(i),i=1,6*NAtm3)
  end if
end if

return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Save localmode.dat (part II).
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine SavLOC(iloc,irep,NAtm3,IExpt,Rslt,PGNAME,ctmp)
implicit real(kind=8) (a-h,o-z)
real(kind=8) :: Rslt(NAtm3,*)
character*4 :: PGNAME(2),ctmp

Ifrq = 3
if(IExpt == 1) Ifrq = 5

write(iloc,"(' $RSLT  $END')")
do i=1, NAtm3
! k, mr, freq, I.R. Int in a.u.; see subroutine PrtNFq for the conversion factors
  X = Rslt(i,1)
  if(IExpt == 1) X = sign(Rslt(i,Ifrq)*Rslt(i,Ifrq), Rslt(i,Ifrq)) * Rslt(i,2)
  write(iloc,"(5d20.10)") X,Rslt(i,2),Rslt(i,Ifrq),Rslt(i,4)
end do

write(iloc,"(' $SYMM  $END')")
write(iloc,"(2a4)")PGNAME(1),PGNAME(2)
rewind(irep)
do i=1,NAtm3
  read(irep,"(a4)")ctmp
  write(iloc,"(a4)")ctmp
end do

return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Save Molden file (imdn).
!
! The frequencies are saved in Freq(:,3) or Freq(:,5) in a.u.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine SavMDN(imdn,NAtm,NAtm3,NVib,ZA,XYZ,IExpt,AL,Freq)
implicit real(kind=8) (a-h,o-z)
parameter(au2wn=5140.48714376d0,au2ang=0.52917720859d0)
real(kind=8) :: ZA(*),XYZ(3,*),AL(NAtm3,*),Freq(NAtm3,*)
character*3 :: Elm

rewind(imdn)
write(imdn,"('[Molden Format]')")
write(imdn,"('[Atoms] Angs')")
do i=1,NAtm
  iza = nint(ZA(i))
  call ElemZA(1,Elm,iza,Elm)
  write(imdn,"(a3,2x,2i6,2x,3f20.10)") Elm,i,iza,(XYZ(j,i)*au2ang,j=1,3)
end do
write(imdn,"('[GTO]',//)")

write(imdn,"('[FREQ]')")
if(IExpt == 0)then
  j=3
else
  j=5
end if
do i=1,NVib
  write(imdn,"(f14.4)")Freq(i,j)*au2wn
end do
write(imdn,"('[FR-COORD]')")
do i=1,NAtm
  iza = nint(ZA(i))
  call ElemZA(1,Elm,iza,Elm)
  write(imdn,"(a3,2x,3f20.10)")Elm,(XYZ(j,i),j=1,3)
end do
write(imdn,"('[FR-NORM-COORD]')")
do i=1,NVib
  write(imdn,"(2x,'vibration',1x,i8)")i
  write(imdn,"(3f20.10)")(AL(j,i),j=1,NAtm3)
end do

write(imdn,"(/)")

return
end

!--- END

