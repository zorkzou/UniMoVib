!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! UniMoVib: a unified interface for molecular harmonic vibrational frequency calculations.
!
! Wenli Zou,  Email: qcband@gmail.com
! Institute of Modern Physics, Northwest University, Xi¡¯an, China
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
program UniMoVib
implicit real(kind=8) (a-h,o-z)
parameter (NOption=10, NArray=10)
dimension     :: IOP(NOption), IArray(NArray)
character*5   :: ver
character*12  :: dat
character*100 :: ctmp

ver="1.0.0"
dat="MAY 01, 2017"

! 1. Assign I/O
call AsgnIO(Intact,ctmp,iinp,iout)

! 2. Set up IOP and IArray:
!
!    IOP(1)   name of quantum chemistry program: 1 (Gaussian),
!             2 (GAMESS), 3 (Firefly), ...; see RdContrl.
!    IOP(2)   IPrint
!    IOP(3)   0 (IFExp=.False.) or 1 (IFExp=.True.)
!    IOP(4)   0, 1, or -1 (ISymm)
!    IOP(5)   0, 1, or 2 (Isotop)
!    IOP(6)   0~ (ISyTol)
!    IOP(7)   0 (IFSAVE=.False.) or 1 (IFSAVE=.True.)
!    IOP(8)   0, 1, or 2 (MWMode); for MOLDEN (read) only
!    IOP(9)   0 (IFMOLDEN=.False.), 1 (IFMOLDEN=.True.), or -1 (IFMOLDEN=.True.
!             but the file exists; fails)
!    END ---> IOP(NOption)
!
!    IArray:  positions of array: 1 (AMass), 2 (ZA), 3 (XYZ), 4 (FFx;
!             square matrix), 5 (APT), 6 (AL: vib + trans & rot),
!             7 (K_corr)



call head(iout,ver,dat)



end
