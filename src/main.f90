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
parameter (NOption=10, iudt=47, imdn=48, idt0=51, idt1=52, idt2=53, irep=54)
dimension     :: IOP(NOption)
character*5   :: ver
character*12  :: dat
character*200 :: ctmp, cname
logical       :: Intact, ifopen
real(kind=8),allocatable :: AMass(:), ZA(:), XYZ(:), FFx(:), APT(:), AL(:), Scr1(:), Scr2(:), Scr3(:), Scr4(:), Work(:)

ver="1.0.0"
dat="AUG 13, 2017"

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
!-----------------------------------------------------------------------
call RdContrl(iinp,iout,iudt,imdn,Intact,NOption,IOP,ctmp,cname)

!-----------------------------------------------------------------------
!  4. Read $qcdata
!-----------------------------------------------------------------------
call RdQCDt(iinp,iout,idt0,idt1,idt2,Intact,IOP)

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
allocate(AMass(NAtm), ZA(NAtm), XYZ(NAtm3), FFx(NSS), APT(NAtm3*3), AL(NSS))
allocate(Scr1(NSS), Scr2(NSS), Scr3(NSS), Scr4(NSS), Work(NWK))

APT=0.d0
call RdData1(iout,idt0,idt1,idt2,Intact,IOP,NAtm,ctmp,AMass,ZA,XYZ,FFx,APT,Scr1,Scr2,Scr3,Work)

! read atomic masses from input
call RdIsot(iinp,iout,Intact,IOP(4),NAtm,ctmp,AMass)

! check data and print geometry
if(IOP(1) /= -1) call ChkDat(iout,Intact,NAtm,AMass,ZA,XYZ)

!-----------------------------------------------------------------------
! 7. Solve Secular equation in Cartesian coordinates
!    Symmetry is also analyzed therein.
!-----------------------------------------------------------------------
call SolvSec(iinp,iout,irep,iudt,imdn,Intact,IOP,NAtm,NVib,ctmp,AMass,ZA,XYZ,FFx,APT,AL,Scr1,Scr2,Scr3,Scr4,Work,Eig)

!-----------------------------------------------------------------------
! xx. the last step
!-----------------------------------------------------------------------
deallocate(AMass, ZA, XYZ, FFx, APT, AL, Scr1, Scr2, Scr3, Scr4, Work)

call fdate(ctmp)
write(*,"(//,' Job terminated correctly! ',a)")trim(ctmp)
if(Intact)write(iout,"(//,' Job terminated correctly, ',a)") trim(ctmp)

close(iinp)
close(iout)
if(IOP(6) == 1) close(iudt)
if(IOP(8) == 1) close(imdn)
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
