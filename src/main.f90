!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! UniMoVib: a unified interface for molecular harmonic vibrational frequency calculations.
!
! Wenli Zou,  Email: qcband@gmail.com
! Institute of Modern Physics, Northwest University, Xi'an, China
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
program UniMoVib
implicit real(kind=8) (a-h,o-z)
parameter (NOption=18, Nrslt=8, iudt=47, imdn=48, iloc=49, idt0=51, idt1=52, idt2=53, irep=54, ibmt=55)
dimension     :: IOP(NOption)
character*5   :: ver
character*12  :: dat
character*200 :: ctmp, cname
logical       :: Intact, ifopen
real(kind=8),allocatable :: AMass(:), ZA(:), XYZ(:), FFx(:), APT(:), DPol(:), AL(:), Rslt(:), Scr1(:), Scr2(:), Scr3(:), &
  Scr4(:), Work(:)

ver="1.2.0"
dat="Jan 28, 2018"

!-----------------------------------------------------------------------
! 1. Assign I/O
!-----------------------------------------------------------------------
call AsgnIO(Intact,ctmp,iinp,iout)

!-----------------------------------------------------------------------
! 2. Read input file name, and define the file name of output
!-----------------------------------------------------------------------
if(Intact) then
  call head1(6,ver,dat)
  call OpInp(iinp,iout,Intact,ctmp,cname)
else
  cname='job'
end if
open(irep,file='SYMFIL1.rep')

call head2(iout,ver,dat)

!-----------------------------------------------------------------------
! 3. Read $contrl
!    IOP(1)   name of quantum chemistry program: 1 (Gaussian),
!             2 (GAMESS), 3 (Firefly), ...; see RdContrl.
!    IOP(2)   0 (IFConc=.True.) or 1 (IFConc=.False.)
!    IOP(3)   0 (IFExp=.False.) or 1 (IFExp=.True.)
!    IOP(4)   0, 1, or 2 (Isotop)
!    IOP(5)   0~ (ISyTol)
!    IOP(6)   0 (IFSAVE=.False.) or 1 (IFSAVE=.True.)
!    IOP(7)   0 (IFMOLDEN=.False.), 1 (IFMOLDEN=.True.)
!    IOP(8)   0 (IFLOCAL=.False.), 1 (IFLOCAL=.True.)
!    IOP(9)   0 (IFRdNM=.False.), 1 (IFRdNM=.True.)
!    IOP(10)  0 (IFApprx=.False.), NPar (IFApprx=.True.)
!-----------------------------------------------------------------------
call RdContrl(iinp,iout,iudt,imdn,iloc,Intact,NOption,IOP,ctmp,cname)

!-----------------------------------------------------------------------
!  4. Read $qcdata
!-----------------------------------------------------------------------
call RdQCDt(iinp,iout,idt0,idt1,idt2,ibmt,Intact,IOP)

!-----------------------------------------------------------------------
! 5. Read NAtm
!-----------------------------------------------------------------------
call RdNAtm1(idt0,idt1,Intact,IOP,NAtm,ctmp)

!-----------------------------------------------------------------------
! 6. Read data: AMass, ZA, XYZ, FFx, APT
!-----------------------------------------------------------------------
NAtm3=3*NAtm
NSS=NAtm3*NAtm3
NWK=2*max(NAtm3,2)*NAtm3
allocate(AMass(NAtm), ZA(NAtm), XYZ(NAtm3), FFx(NSS), APT(NAtm3*3), DPol(NAtm3*6), AL(NSS), stat=ierr)
  if(ierr /= 0) call XError(Intact,"Insufficient Memory (1)!")
allocate(Rslt(NAtm3*Nrslt), Scr1(NSS), Scr2(NSS), stat=ierr)
  if(ierr /= 0) call XError(Intact,"Insufficient Memory (2)!")
if(IOP(9) /= 1) allocate(Scr3(NSS), Scr4(NSS), Work(NWK), stat=ierr)
  if(ierr /= 0) call XError(Intact,"Insufficient Memory (3)!")

APT=0.d0
call RdData1(iout,idt0,idt1,idt2,ibmt,Intact,IOP,Infred,IRaman,NAtm,ctmp,AMass,ZA,XYZ,FFx,APT,DPol,Scr1,Scr2,Scr3,Work)

! read atomic masses from input
call RdIsot(iinp,iout,Intact,IOP(4),NAtm,ctmp,AMass)

if(IOP(1) /= -1) then
! check data
  call ChkDat(iout,Intact,NAtm,AMass,ZA,XYZ)
! print geometry and probably atomic IR charges
  call PrtCoord(iout,Infred,NAtm,AMass,ZA,XYZ,APT)
end if

!-----------------------------------------------------------------------
! 7. Solve Secular equation in Cartesian coordinates
!    Symmetry is also analyzed therein.
!-----------------------------------------------------------------------
call SolvSec(iinp,iout,idt0,irep,iudt,imdn,iloc,Intact,IOP,Infred,IRaman,NAtm,NVib,ctmp,AMass,ZA,XYZ,FFx,APT,DPol,AL,Rslt, &
  Scr1,Scr2,Scr3,Scr4,Work)

!-----------------------------------------------------------------------
! xx. the last step
!-----------------------------------------------------------------------
deallocate(AMass, ZA, XYZ, FFx, APT, DPol, AL, Rslt, Scr1, Scr2)
if(IOP(9) /= 1) deallocate(Scr3, Scr4, Work)

call fdate(ctmp)
write(*,"(//,' UniMoVib job terminated correctly! ',a)")trim(ctmp)
write(*,"(//)")
if(Intact)write(iout,"(//,' UniMoVib job terminated correctly, ',a)") trim(ctmp)

close(iinp)
close(iout)
if(IOP(6) == 1) close(iudt)
if(IOP(7) == 1) close(imdn)
if(IOP(8) == 1) close(iloc)
inquire(unit=idt0,opened=ifopen)
if(ifopen) close(idt0)
inquire(unit=idt1,opened=ifopen)
if(ifopen) close(idt1)
inquire(unit=idt2,opened=ifopen)
if(ifopen) close(idt2)
close(irep,status='delete')

call estop(Intact)

end

!--- END
