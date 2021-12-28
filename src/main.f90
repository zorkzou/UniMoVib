!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! UniMoVib: a unified interface for molecular harmonic vibrational frequency calculations.
!
! Wenli Zou,  Email: qcband@gmail.com
! Institute of Modern Physics, Northwest University, Xi'an, Shaanxi, China
!
! and
!
! Yunwen Tao, Email: ywtao.smu@gmail.com
! Department of Chemistry, Southern Methodist University, Dallas, TX, USA
!
! The UniMoVib program was originally written by Wenli Zou in FORTRAN 77 during 2014 and 2015 at Southern Methodist University
! (SMU), Dallas, Texas, within the framework of the LocalMode program of the Computational and Theoretical Chemistry Group (CATCO)
! of SMU. This work was supported by the NSF grants CHE 1152357 and CHE 1464906. Guidance from the late Dr. Dieter Cremer is
! acknowledged. After being rewritten in Fortran 90 in the spring of 2017, UniMoVib has been released as a stand-alone program.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
program UniMoVib
  implicit real(kind=8) (a-h,o-z)
  parameter (NOption=18, Nrslt=8, iudt=47, imdn=48, iloc=49, igau=50, idt0=51, idt1=52, idt2=53, irep=54, ireo=55, &
                         ibmt=56, isva=57, ixyz=58, imdf=59, imds=60 )
  dimension     :: IOP(NOption)
  character*5   :: ver
  character*12  :: dat
  character*200 :: ctmp, cname, tag
  logical       :: Intact, ifopen
  allocatable   :: AMass(:), ZA(:), XYZ(:), Grd(:), FFx(:), APT(:), DPol(:), AL(:), Rslt(:)
  allocatable   :: Scr1(:), Scr2(:), Scr3(:), Scr4(:), Work(:)
  logical, allocatable   :: LScr(:)
  integer, allocatable   :: IScr(:)
  character*1, allocatable :: CScr(:)
  integer, allocatable   :: subsystem_idx(:),flags(:)
  allocatable   :: AMass_sub(:), XYZ_sub(:), ZA_sub(:)   

  ver="1.4.4"
  dat="Dec 28, 2021"

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
!    IOP(11)  0 (IFGauTS=.False.), 1 (IFGauTS=.True.)
!    IOP(12)  0 (IFSymtz=.False.), 1 (IFSymtz=.True.)
!    IOP(15)  0 (IFGSVA=.False.), 1 (IFGSVA=.True.)
!    IOP(16)  0 (IFPYVIBMS=.False.), 1 (IFPYVIBMS=.True.)
!-----------------------------------------------------------------------
  call RdContrl(iinp,iout,iudt,imdn,iloc,igau,Intact,NOption,IOP,ctmp,cname)

!-----------------------------------------------------------------------
!  4. Read $qcdata
!-----------------------------------------------------------------------
  call RdQCDt(iinp,iout,idt0,idt1,idt2,ibmt,Intact,IOP)

!-----------------------------------------------------------------------
! 5. Read NAtm
!-----------------------------------------------------------------------
  call RdNAtm1(iinp,idt0,idt1,Intact,IOP,NAtm,tag,ctmp)

!-----------------------------------------------------------------------
! 6. Read data: AMass, ZA, XYZ, FFx, APT, and DPol
!-----------------------------------------------------------------------
  NAtm3=3*NAtm
  NSS=NAtm3*NAtm3
  NWK=2*NSS
  allocate(AMass(NAtm), ZA(NAtm), XYZ(NAtm3), Grd(NAtm3), FFx(NSS), APT(NAtm3*3), DPol(NAtm3*6), AL(NSS), stat=ierr)
    if(ierr /= 0) call XError(Intact,"Insufficient Memory (1)!")
  allocate(Rslt(NAtm3*Nrslt), Scr1(NSS), Scr2(NSS), LScr(NAtm3), IScr(9*NAtm), CScr(12*NAtm), stat=ierr)
    if(ierr /= 0) call XError(Intact,"Insufficient Memory (2)!")
  if(IOP(9) /= 1) allocate(Scr3(NSS), Scr4(NSS), Work(NWK), stat=ierr)
    if(ierr /= 0) call XError(Intact,"Insufficient Memory (3)!")

  APT=0.d0
  call RdData1(iinp,iout,idt0,idt1,idt2,ibmt,Intact,IOP,Infred,IRaman,IGrd,NAtm,tag,ctmp,AMass,ZA,XYZ,Grd,FFx,APT,DPol, &
    Scr1,Scr2,Scr3,Work)

! read atomic masses from input
  call RdIsot(iinp,iout,Intact,IOP(4),NAtm,ctmp,AMass)

  if(IOP(1) /= -1) then
!   check data
    call ChkDat(iout,Intact,NAtm,AMass,ZA,XYZ)
!   print geometry and probably atomic IR charges
    call PrtCoord(iout,Infred,NAtm,AMass,ZA,XYZ,APT,iop(16),ixyz)
  end if

!-----------------------------------------------------------------------
! 7. Solve Secular equation in Cartesian coordinates
!    Symmetry is also analyzed therein.
!-----------------------------------------------------------------------
  call SolvSec(iinp,iout,idt0,irep,ireo,iudt,imdn,iloc,igau,imdf,Intact,IOP,Infred,IRaman,IGrd,NAtm,NVib,ctmp,AMass,ZA,XYZ, &
    Grd,FFx,APT,DPol,AL,Rslt,LScr,IScr,CScr,Scr1,Scr2,Scr3,Scr4,Work)

!-----------------------------------------------------------------------
! 8. GSVA - Generalized Subsystem Vibrational Analysis
!    J. Chem. Theory Comput. 2018, 14(5), 2558-2569
!-----------------------------------------------------------------------
  if(IOP(15) == 1)then
   
   allocate(flags(NAtm))
   call RdGSVA(iinp,iout,NAtm,flags,NAtm_sub)
   !GSVA: partition the data vectors for the subsystem
    
   allocate(subsystem_idx(NAtm_sub),AMass_sub(NAtm_sub),ZA_sub(NAtm_sub),XYZ_sub(3*NAtm_sub))
   call obt_idx(iout,flags,subsystem_idx,NAtm) 
   
   !PRINT '(5I5)',(subsystem_idx(i),i=1,3) 
   call PartMX(NAtm,AMass,XYZ,ZA,subsystem_idx,NAtm_sub,AMass_sub,XYZ_sub,ZA_sub) 
 
   call GSVA_engine(iout,isva,imds,IOP,NAtm,subsystem_idx,NAtm_sub,AMass_sub,XYZ_sub,ZA_sub,FFx)
  end if

!-----------------------------------------------------------------------
! 99. the last step
!-----------------------------------------------------------------------
  deallocate(AMass, ZA, XYZ, Grd, FFx, APT, DPol, AL, Rslt, Scr1, Scr2, LScr, IScr, CScr)
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
  if(IOP(11)== 1) close(igau)
  if(IOP(15)== 1) close(isva)
  if(IOP(16)== 1) close(ixyz)
  if(IOP(16)== 1) close(imdf)
  if(IOP(15)== 1 .and. IOP(16)== 1) close(imds)
  inquire(unit=idt0,opened=ifopen)
  if(ifopen) close(idt0)
  inquire(unit=idt1,opened=ifopen)
  if(ifopen) close(idt1)
  inquire(unit=idt2,opened=ifopen)
  if(ifopen) close(idt2)
  close(irep,status='delete')
  inquire(unit=ireo,opened=ifopen)
  if(ifopen) close(ireo,status='delete')

  call estop(Intact)

end

!--- END

