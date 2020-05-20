subroutine RdGSVA(iinp,iout,NAtm,flags,NAtm_sub)
implicit real(kind=8) (a-h,o-z)
character*1000 ::subsystem!,pair(2)
character*1 ::hyphen
character(20), allocatable :: strarray(:),pair(:)
character(20) :: tmpstr,before,after
integer,allocatable :: subsystem_idx(:)
integer :: flags(*)
namelist/GSVA/subsystem

hyphen = '-'

Do I=1,NAtm
   flags(I) = 0
End Do

rewind(iinp)
read(iinp,GSVA,end=100,err=10)
goto 100
10    call XError(.True.,"$GSVA is wrong!")

100   continue

!debug
!print *,trim(subsystem)

n = count(transfer(subsystem, 'a', len(subsystem)) == ",")
allocate(strarray(n+1))
read(subsystem, *) strarray(1:n+1)
!print *, 'nvalues=', n+1
!print '(a)', strarray(1:n+1)


Do I=1,n+1
   if(index(strarray(I),hyphen) == 0 )then
     read (strarray(I),'(I10)') IAt
     if ((IAt < 1) .or. (IAt > NAtm))then
        call XError(.True.,"Invalid atom label for GSVA!")
     end if  
     !print *,IAt  
     flags(IAt) = 1 
   else
     !print *,strarray(I)
     tmpstr = trim(strarray(I))
     !print *,tmpstr
     !print *,len(tmpstr)
     Ihyph = index(tmpstr,hyphen)
     before = tmpstr(1:Ihyph-1)
     after  = tmpstr(Ihyph+1:len(tmpstr))
     
     read (before,'(I10)') IAt1     
     read (after,'(I10)') IAt2
     !print *,IAt1
     !print *,IAt2

     if ((IAt1 < 1) .or. (IAt1 > NAtm) .or. (IAt1 >= IAt2) )then
        call XError(.True.,"Invalid atom label for GSVA!")
     end if  
     if ((IAt2 < 1) .or. (IAt2 > NAtm))then
        call XError(.True.,"Invalid atom label for GSVA!")
     end if  

     Do J=IAt1,IAt2
        flags(J) = 1
     End Do 

   end if
End Do

NAtm_sub=0
Do I=1,NAtm
  if(flags(I) == 1)then
     NAtm_sub=NAtm_sub+1  
  end if
End Do


write(iout,"(//,' -- Generalized Subsystem Vibrational Analysis (GSVA) --',//,&
' Subsystem has',1x,1I5,2x, 'atoms:')")(NAtm_sub)

!debug
!PRINT '(10I10)',(flags(j),j=1,NAtm)

return 
end 

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Print Head-1
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine head1(iout,ver,dat)
implicit real(kind=8) (a-h,o-z)
character*5   :: ver
character*12  :: dat

write(iout,"(1x,70('='),/, &
' ==       UniMoVib: A unified interface for molecular harmonic       ==',/, &
' ==         vibrational frequency calculations.                      ==',/, &
' ==',66x,'==',/, &
' ==       Ver. ',a5,', ',a12,35x,'==',/, &
' ==       Webpage: https://github.com/zorkzou/UniMoVib               ==',/, &
1x,70('='),/ )") ver,dat

return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Print Head-2
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine head2(iout,ver,dat)
implicit real(kind=8) (a-h,o-z)
character*5   :: ver
character*12  :: dat

write(iout,"(1x,119('='),/, &
'       _/     _/                _/     _/            _/     _/',/, &
'      _/     _/            _/  _/|  _/_/            _/     _/  _/  _/        UniMoVib:          A unified interface for',/, &
'     _/     _/   _/_/_/       _/ |_/ _/   _/_/_/   _/     _/      _/           molecular harmonic vibrational frequency',/, &
'    _/     _/  _/    _/  _/  _/  _/ _/  _/    _/  _/     _/  _/  _/_/_/        calculations.',/, &
'   _/     _/  _/    _/  _/  _/     _/  _/    _/  _/    _/   _/  _/   _/',/, &
'  _/     _/  _/    _/  _/  _/     _/  _/    _/   _/  _/    _/  _/   _/       Ver. ',a5,', ',a12,/, &
'   _/_/_/   _/    _/  _/  _/     _/   _/_/_/     _/_/     _/  _/_/_/         https://github.com/zorkzou/UniMoVib',/, &
1x,119('='),/ )") ver,dat

return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Assign ports of I/O
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine AsgnIO(Intact,ctmp,iinp,iout)
implicit real(kind=8) (a-h,o-z)
parameter(NArg=1)
logical :: Intact
character*100 :: ctmp

! read arguments
do i=1,NArg+1
  call get_command_argument(i,ctmp)
!   call getarg(i,ctmp)
  call charl2u(ctmp)
  istr=nonspace(ctmp)
  iend=len_trim(ctmp)

  if(i == 1)then
    if(iend == 0) then
      Intact = .True.
    else if(ctmp(istr:iend) == '-B') then
      Intact = .False.
    else
      Intact = .False.
      call XError(Intact,"Argument #1 cannot be recognized!")
    end if
  else if(i == 2)then
    if(iend /= 0) then
      Intact = .False.
      call XError(Intact,"Argument #2 cannot be recognized!")
    end if
  end if
end do

if(Intact)then
  iinp=42
  iout=43
else
  iinp=5
  iout=6
end if

return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Read input file name, and define the file name of output
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine OpInp(iinp,iout,Intact,fname,cname)
implicit real(kind=8) (a-h,o-z)
logical :: Intact
character*200 :: fname, cname

100   write(*,"(/,' Type in file name within 100 characters (default: job.inp):',/,' > ',$)")
read(*,"(a100)")fname

if(LEN_TRIM(fname) == 0) fname='job.inp'
istr=nonspace(fname)
iend=LEN_TRIM(fname)
open(iinp,file=fname(istr:iend),status='old',err=110)
goto 200

110   write(*,"(//,' This input file does not exist: ',a)")fname(istr:iend)
write(*,"(' Please try again.',/)")
goto 100

! define the output file name
200   write(*,"(//,' Input file:',6x,a)")fname(istr:iend)
idot=index(fname(istr:iend), '.', .True.)
if(idot == 0) idot = iend + 1
open(iout,file=fname(istr:idot-1)//'.out',status='new',err=210)
write(*,"(' Output file:',5x,a)")fname(istr:idot-1)//'.out'
goto 300

210   write(*,"(//,' The output file exists: ',a,/,' I am not sure whether it is still useful for you.')") &
  fname(istr:idot-1)//'.out'
call XError(Intact,"Please delete or rename it first.")

300   rewind(iout)
cname=fname(istr:idot-1)

return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Read the $Contrl group.
!
! QCProg = IOP(1) =
!       XYZ (*)        -3,       UniMoVib (ALM) -2,       AtomCalc       -1,
!       Gaussian        1,       GAMESS          2,       Firefly         3,       ORCA            4,       CFour           5,
!       Molpro          6,       QChem           7,       NWChem          8,       GAMESSUK        9,       TURBOMOLE      10,
!       deMon          11,       PQS            12,       MOPAC          13,       AMPAC/AMSOL    14,       Dalton         15,
!       FHI-AIMS       16,       CP2k           17,       ADF            18,       HyperChem      19,       Jaguar         20,
!       MOLDEN         21,       Crystal        22,       Spartan        23,       PSI4           24,       DMol3          25,
!       ACES           26
!
! (*) For debugging only.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RdContrl(iinp,iout,iudt,imdn,iloc,igau,Intact,NOp,IOP,qcprog,cname)
implicit real(kind=8) (a-h,o-z)
parameter(NProg=26)
dimension :: IOP(NOp)
logical :: Intact,ifconc,ifexp,ifsave,ifmolden,iflocal,ifrdnm,ifapprx,ifgauts,ifsymtz,ifgsva
namelist/Contrl/qcprog,ifconc,Isotop,ISyTol,ifexp,ifsave,ifmolden,iflocal,ifrdnm,ifapprx,ifgauts,ifsymtz,ifgsva
character*200 :: qcprog, cname
character*9,allocatable  :: DATFMT(:)

IOP=0
! default values
qcprog   = 'GAUSSIAN'
Isotop   = 0
ISyTol   = 0
IApprox  = 0
ifconc   = .false.
ifexp    = .false.
ifsave   = .false.
ifmolden = .false.
iflocal  = .false.
ifrdnm   = .false.
ifapprx  = .false.
ifgauts  = .false.
ifsymtz  = .false.
ifgsva   = .false.

rewind(iinp)
read(iinp,Contrl,end=100,err=10)
goto 100
10    call XError(Intact,"$Contrl is wrong!")

100   continue
!>>>  qcprog --> IOP(1)
if(LEN_TRIM(qcprog) == 0) qcprog = 'GAUSSIAN'
call charl2u(qcprog)

!>>> Atomic thermochemistry calculation
if(index(qcprog,'ATOMCALC') /= 0) then
  IOP(1)=-1
  write(iout,"(/,' Atomic thermochemistry calculation')")
  if(Intact) write(*,"(' Atomic thermochemistry calculation')")
  IOP(4)= 2
! ignore the following options
  return

!>>> UniMoVib (ALM)
else if(index(qcprog,'UNIMOVIB') /= 0 .or. index(qcprog,'ALM') /= 0) then
  IOP(1)=-2
  write(iout,"(/,' Data Format:',5x,'UniMoVib')")
  if(Intact) write(*,"(' Data Format:',5x,'UniMoVib')")

!>>> XYZ
else if(index(qcprog,'XYZ') /= 0) then
  IOP(1)=-3
  write(iout,"(/,' Data Format:',5x,'XYZ')")
  if(Intact) write(*,"(' Data Format:',5x,'XYZ')")

!>>> GAUSSIAN
else if(index(qcprog,'GAUSSIAN') /= 0) then
  IOP(1)=1

!>>> FIREFLY
else if(index(qcprog,'FIREFLY') /= 0 .or. index(qcprog,'PCGAMESS') /= 0 .or. index(qcprog,'PC-GAMESS') /= 0) then
  IOP(1)=3

!>>> GAMESSUK
else if(index(qcprog,'GAMESSUK') /= 0 .or. index(qcprog,'GAMESS-UK') /= 0) then
  IOP(1)=9

!>>> GAMESS
else if(index(qcprog,'GAMESS') /= 0) then
  IOP(1)=2

!>>> ORCA
else if(index(qcprog,'ORCA') /= 0) then
  IOP(1)=4

!>>> CFOUR
else if(index(qcprog,'CFOUR') /= 0) then
  IOP(1)=5

!>>> MOLPRO
else if(index(qcprog,'MOLPRO') /= 0) then
  IOP(1)=6

!>>> QCHEM
else if(index(qcprog,'QCHEM') /= 0 .or. index(qcprog,'Q-CHEM') /= 0) then
  IOP(1)=7

!>>> NWCHEM
else if(index(qcprog,'NWCHEM') /= 0) then
  IOP(1)=8

!>>> TURBOMOLE
else if(index(qcprog,'TURBOMOLE') /= 0) then
  IOP(1)=10

!>>> DEMON2K
else if(index(qcprog,'DEMON') /= 0) then
  IOP(1)=11

!>>> PQS
else if(index(qcprog,'PQS') /= 0) then
  IOP(1)=12

!>>> MOPAC
else if(index(qcprog,'MOPAC') /= 0) then
  IOP(1)=13

!>>> AMPAC/AMSOL
else if(index(qcprog,'AMPAC') /= 0 .or. index(qcprog,'AMSOL') /= 0) then
  IOP(1)=14

!>>> DALTON
else if(index(qcprog,'DALTON') /= 0) then
  IOP(1)=15

!>>> FHI-AIMS
else if(index(qcprog,'AIMS') /= 0) then
  IOP(1)=16

!>>> CP2k
else if(index(qcprog,'CP2K') /= 0) then
  IOP(1)=17

!>>> ADF
else if(index(qcprog,'ADF') /= 0) then
  IOP(1)=18

!>>> Hyperchem
else if(index(qcprog,'HYPERCHEM') /= 0) then
  IOP(1)=19

!>>> JAGUAR
else if(index(qcprog,'JAGUAR') /= 0) then
  IOP(1)=20

!>>> MOLDEN
else if(index(qcprog,'MOLDEN') /= 0) then
  IOP(1)=21

!>>> CRYSTAL
else if(index(qcprog,'CRYSTAL') /= 0) then
  IOP(1)=22

!>>> SPARTAN
else if(index(qcprog,'SPARTAN') /= 0) then
  IOP(1)=23

!>>> PSI4
else if(index(qcprog,'PSI') /= 0) then
  IOP(1)=24

!>>> DMol3
else if(index(qcprog,'DMOL') /= 0) then
  IOP(1)=25

!>>> ACES
else if(index(qcprog,'ACES') /= 0) then
  IOP(1)=26

!>>> Unknown
else
  call XError(Intact,"Unknown quantum chemistry program or data format!")
end if

if (IOP(1) > 0) then
  allocate(DATFMT(NProg))
  DATFMT(1:NProg)=(/"Gaussian ","GAMESS   ","Firefly  ","ORCA     ","CFOUR    ",&
                    "MOLPRO   ","Q-CHEM   ","NWCHEM   ","GAMESS-UK","TURBOMOLE",&
                    "deMon2K  ","PQS      ","MOPAC    ","AMSOL    ","DALTON   ",&
                    "FHI-AIMS ","CP2K     ","ADF      ","HyperChem","Jaguar   ",&
                    "MOLDEN   ","CRYSTAL  ","SPARTAN  ","PSI4     ","Dmol3    ",&
                    "ACES     "/)
  write(iout,"(/,' QC Program:',6x,a)") trim( DATFMT(IOP(1)) )
  if(Intact) write(*,"(' QC Program:',6x,a)") trim( DATFMT(IOP(1)) )
	deallocate(DATFMT)
end if

!>>>  ifconc --> IOP(2)
IOP(2) = 1
if(ifconc) IOP(2) = 0

!>>>  ifexp --> IOP(3) = 0 (ifexp=.false.) or 1 (ifexp=.true.)
if(ifexp) IOP(3) = 1

!>>>  Isotop --> IOP(4) = 0, 1, or 2
if(Isotop == 1 .or. Isotop == 2) IOP(4) = Isotop

!>>>  ISyTol=MN --> IOP(5), where M ~ [1,9] and N ~ [-9,9]
N = MOD(ISyTol,10)
N = SIGN(N,ISyTol)
M = ABS(ISyTol)/10
M = MOD(M,10)
if(M == 0) M = 1
IOP(5) = SIGN(M*10+ABS(N),N)

!>>>  ifsave --> IOP(6)
if(ifsave) IOP(6) = 1
! UniMoVib (ALM) data file
if(IOP(6) == 1)then
  istr=nonspace(cname)
  iend=LEN_TRIM(cname)

  open(iudt,file=cname(istr:iend)//'.umv')
  write(iout,"(' UniMoVib file:',3x,a)") cname(istr:iend)//'.umv'
  if(Intact) write(*,"(' UniMoVib file:',3x,a)") cname(istr:iend)//'.umv'
end if

!>>>  ifmolden --> IOP(7) = 0 (ifmolden=.false.) or 1 (ifmolden=.true.)
if(ifmolden)then
  IOP(7) = 1
  istr=nonspace(cname)
  iend=LEN_TRIM(cname)
  open(imdn,file=cname(istr:iend)//'.molden',status='new',iostat=ioper)
  if(ioper == 0)then
    write(iout,"(' MOLDEN file:',5x,a)") cname(istr:iend)//'.molden'
    if(Intact) write(*,"(' MOLDEN file:',5x,a)") cname(istr:iend)//'.molden'
  else
    open(imdn,file=cname(istr:iend)//'_new.molden',status='new',iostat=ioper)
    if(ioper > 0) then
      write(iout,"(/,' You have to delete or rename this file first:',5x,a)") cname(istr:iend)//'_new.molden'
      call XError(Intact,"New molden file exist!")
    end if
    write(iout,"(' MOLDEN file:',5x,a)") cname(istr:iend)//'_new.molden'
    if(Intact) write(*,"(' MOLDEN file:',5x,a)") cname(istr:iend)//'_new.molden'
  end if
end if

!>>>  iflocal --> IOP(8) = 0 (iflocal=.false.) or 1 (iflocal=.true.)
if(iflocal)then
  IOP(8) = 1
  open(iloc,file='localmode.dat')
  write(iout,"(' LOCALMODE file:',2x,'localmode.dat')")
  if(Intact) write(*,"(' LOCALMODE file:',2x,'localmode.dat')")
end if

!>>>  ifrdnm --> IOP(9) = 0 (ifrdnm=.false.) or 1 (ifrdnm=.true.)
if(ifrdnm) then
! Compatibility: IFExp
  if(IOP(3) /= 0) call XError(Intact,"IFRdNM is incompatible with IFExp=.True.!")
! Compatibility: Isotop
  if(IOP(4) /= 0) call XError(Intact,"IFRdNM is incompatible with Isotop!")
! Compatibility: QCProg
  if(IOP(1) /= 1) call XError(Intact,"At present, IFRdNM supports QCProg=Gaussian only!")
  IOP(9) = 1
end if

!>>>  ifapprx --> IOP(10) = 0 (ifapprx=.false.) or 1 (ifapprx=.true.)
if(ifapprx) then
! Compatibility: ifrdnm
  if(IOP(9) /= 0) call XError(Intact,"IFApprx is incompatible with IFRdNM!")
  IOP(10) = 1
end if

!>>>  ifgauts --> IOP(11) = 0 (ifgauts=.false.) or 1 (ifgauts=.true.)
if(ifgauts) then
  IOP(11) = 1
  open(igau,file='gaussian-ts.gjf')
  write(iout,"(' TS input file:',3x,'gaussian-ts.gjf')")
  if(Intact) write(*,"(' TS input file:',3x,'gaussian-ts.gjf')")
end if

!>>>  ifsymtz --> IOP(12) = 0 (ifsymtz=.false.) or 1 (ifsymtz=.true.)
if(ifsymtz) then
! Compatibility: ifrdnm
  if(IOP(9) /= 0) call XError(Intact,"IFSymtz is incompatible with IFRdNM!")
  IOP(12) = 1
end if

!GSVA
if(ifgsva) IOP(15) = 1

return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Read $QCData group
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RdQCDt(iinp,iout,idt0,idt1,idt2,ibmt,Intact,IOP)
implicit real(kind=8) (a-h,o-z)
dimension :: IOP(*)
logical :: Intact
character*100 :: fchk, hess, ddip, geom, bmat
namelist/QCData/fchk,hess,ddip,geom,bmat

! atomic calculation
if(IOP(1) == -1) return

fchk=' '
hess=' '
ddip=' '
geom=' '
bmat=' '
rewind(iinp)
read(iinp,QCData,end=100,err=10)
goto 100
10    call XError(Intact,"$QCData is wrong!")

100   continue
! default name of $(fchk)
if(LEN_TRIM(fchk) == 0) then
  if(IOP(1) == 3) then
    fchk='PUNCH'
  else if(IOP(1) == 10) then
    fchk='aoforce.out'
  else if(IOP(1) == 11) then
    fchk='deMon.out'
  end if
end if

! open $(fchk)
call DFOpen(iout,idt0,Intact,.True.,'FCHK',fchk)

! Gamess, Firefly: open $(geom)
if(IOP(1) == 2 .or. IOP(1) == 3) then
  if(LEN_TRIM(geom) == 0) geom='gamess.out'
  call DFOpen(iout,idt1,Intact,.True.,'GEOM',geom)

! Cfour: open $(geom); optional
else if(IOP(1) == 5) then
  if(LEN_TRIM(geom) == 0) geom='GRD'
  call DFOpen(iout,idt1,Intact,.False.,'GEOM',geom)

! NWChem, PQS, FHI-AIMS: open $(hess) and optional $(ddip)
else if(IOP(1) == 8 .or. IOP(1) == 12 .or. IOP(1) == 16) then
! open $(hess)
  call DFOpen(iout,idt1,Intact,.True.,'HESS',hess)
! open $(ddip)
  call DFOpen(iout,idt2,Intact,.False.,'DDIP',ddip)

! Turbomole: open optional $(ddip)
else if(IOP(1) == 10) then
  if(LEN_TRIM(ddip) == 0) ddip='dipgrad'
! open $(ddip)
  call DFOpen(iout,idt1,Intact,.False.,'DDIP',ddip)

! IApprox: open $(bmat)
else if(IOP(10) /= 0) then
  call DFOpen(iout,ibmt,Intact,.True.,'BMAT',bmat)
end if

return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Read atomic masses from input
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RdIsot(iinp,iout,Intact,IRdMas,NAtm,ctmp,AMass)
implicit real(kind=8) (a-h,o-z)
real(kind=8) :: AMass(*)
character*100 :: ctmp
namelist/Isomas/null
logical :: Intact

if(IRdMas == 0) return

rewind(iinp)
read(iinp,Isomas,end=1000,err=1000)

if(IRdMas == 1)then
  nmass=0
  do while(.true.)
    read(iinp,"(a100)",end=100,err=1010)ctmp

    if(len_trim(ctmp) == 0) goto 100
    read(ctmp,*,err=1010)iatom,am
    if(iatom < 1 .or. iatom > NAtm)then
      write(iout,"(/,2x,'IAtom = ',i4)")iatom
      goto 1020
    end if
    AMass(iatom) = am
    nmass=nmass+1
    cycle

    100 exit
  end do
  if(nmass < 1 .or. nmass > NAtm)then
    write(iout,"(/,2x,'NMass = ',i4)")nmass
    goto 1030
  end if
else
  read(iinp,*,err=1010,end=1010)(AMass(i),i=1,NAtm)
end if

return
1000  if(NAtm == 1) then
  call XError(Intact,"Atomic mass is not defined correctly!")
else
  call XError(Intact,"The $Isomas group cannot be read!")
end if

1010  call XError(Intact,"Please check your input of atomic masses!")
1020  call XError(Intact,"IAtom in $Isomas is out of range!")
1030  call XError(Intact,"#Mass in $Isomas is out of range!")
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! check data
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine ChkDat(iout,Intact,NAtm,AMass,ZA,XYZ)
implicit real(kind=8) (a-h,o-z)
parameter(tolmass=1.d-2,tolz=0.99d0,tolr=0.5d0)
real(kind=8) :: AMass(*),ZA(*),XYZ(3,*)
logical :: Intact

! check AMass and ZA
do i=1,NAtm
  if(AMass(i) < tolmass)then
    write(iout,"(/,' Mass(',i3,') = ',f9.5)")i,AMass(i)
    call XError(Intact,"Atomic mass is not reasonable!")
  end if
  if(ZA(i) < tolz)then
    write(iout,"(/,' Z(',i3,') = ',f5.1)")i,ZA(i)
    call XError(Intact,"Dummy atom is not allowed!")
  end if
end do

! check XYZ
do i=1,NAtm-1
  do j=i+1,NAtm
    r=distance(XYZ(1,i),XYZ(1,j))
    if(r < tolr)then
      write(iout,"(/,' R(',i3,',',i3,') = ',f6.4)")i,j,r
      call XError(Intact,"Distance is too short!")
    end if
  end do
end do

call RmNoise(NAtm*3,1.0d-8,XYZ)

return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! print geometry and probably atomic IR charges
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine PrtCoord(iout,Infred,NAtm,AMass,ZA,XYZ,APT)
implicit real(kind=8) (a-h,o-z)
parameter(au2ang=0.52917720859d0)
real(kind=8) :: AMass(*),ZA(*),XYZ(3,*),APT(3,3,*), ip(3)
character*3 :: Elm

call planar(NAtm,XYZ,ip)

if(Infred == 0 .or. sum(ip) == 0) then
  write(iout,"(//,' Cartesian coordinates (Angstrom)',/,1x,90('-'),/,3x, &
    'No.   Atom    ZA                 X             Y             Z                Mass',/,1x,90('-'))")

  do i=1,NAtm
    iza = nint(ZA(i))
    call ElemZA(1,Elm,iza,Elm)
    write(iout,"(i6,4x,a3,1x,i5,8x,3f14.8,8x,f14.8)") i,Elm,iza,(XYZ(j,i)*au2ang,j=1,3),AMass(i)
  end do

  write(iout,"(1x,90('-'))")
else
  write(iout,"(//,' Cartesian coordinates (Angstrom)',/,1x,108('-'),/,3x, &
    'No.   Atom    ZA                 X             Y             Z                Mass               IR charge',/,1x,108('-'))")

  do i=1,NAtm
    iza = nint(ZA(i))
    Charge = AIRCrg(ip,APT(1,1,i))
    call ElemZA(1,Elm,iza,Elm)
    write(iout,"(i6,4x,a3,1x,i5,8x,3f14.8,8x,f14.8,9x,f9.4)") i,Elm,iza,(XYZ(j,i)*au2ang,j=1,3),AMass(i),Charge
  end do

  write(iout,"(1x,108('-'),/,' Reference of IR charge:',/, &
    ' A. Milani, M. Tommasini, C. Castiglioni, Theor. Chem. Acc. 131, 1139 (2012).')")
end if

return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! save irreps in ida1 to ida2.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine dumpir(ida1,ida2,c4)
 implicit real(kind=8) (a-h,o-z)
 character*4   :: c4

 rewind(ida1)
 read(ida1,"(2i6)") NClass, Nrt
 do i = 1, NClass
   read(ida1,"(i6,1x,a4)") n, c4
   write(ida2,"(a4)") (adjustl(c4),j=1,n)
 end do

 do i = 1, Nrt
   read(ida1,"(a4)") c4
   write(ida2,"(a4)") adjustl(c4)
 end do

 return
end

!--- END

