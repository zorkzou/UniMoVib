subroutine GSVA_engine(iout,IOP,NAtm,subsystem_idx,NAtm_sub,AMass_sub,XYZ_sub,ZA_sub,FFx)
implicit real(kind=8) (a-h,o-z)
dimension :: IOP(*),FFx(*),AMass_sub(*),XYZ_sub(*),AMass_sub_bak(NAtm_sub),ZA_sub(*)
integer :: subsystem_idx(*)
allocatable :: Scr1(:),Scr2(:),Scr3(:),Scr4(:),AL(:)
allocatable :: FFx_inv(:),Scra(:),Scrb(:),Scrc(:),Scrd(:),prod(:),prod_inv(:)!,Rslt(:)
real(kind=8), dimension (:, :), ALLOCATABLE :: VibSp,VibSp_t,VibSp_full,VibSp_full_t,compl_m1, &
compl_m2,compl_inv,newHes_1,FFx_sub,Rslt

 NAtm3_sub = 3*NAtm_sub
 NSS_sub = NAtm3_sub*NAtm3_sub
 allocate(Scr2(NSS_sub),Scr1(NSS_sub),Scr3(NSS_sub))

!It is required to use identical mass for all atoms
 call acopy(NAtm_sub,AMass_sub,AMass_sub_bak) 
 Do J=1,NAtm_sub 
   AMass_sub(J) = 1.0E0
 End Do

 call MassCent(NAtm_sub,AMass_sub,XYZ_sub,Scr1,Scr2) 
      !PRINT '(10F7.2)',(Scr2(i),i=1,3)

 call MIner(.True.,NAtm_sub,AMass_sub,Scr1,Scr1(NAtm3_sub+1),Scr2,Scr2(4))

 allocate(AL(NSS_sub))
 call TRVec(.True.,NAtm_sub,NTR_sub,Im,AMass_sub,Scr1,AL,Scr1(NAtm3_sub+1),Scr2)
 NVib_sub=NAtm3_sub-NTR_sub
      !PRINT '(10I1)', NVib_sub
 
 call GSorth(.True.,.True.,NAtm3_sub,NTR_sub,AL,Scr3)
      !PRINT '(18F7.2)',(AL(i),i=1+NAtm3_sub*5,NAtm3_sub+NAtm3_sub*5) 
 allocate(VibSp(3*NAtm_sub,NVib_sub))
 call acopy(3*NAtm_sub*NVib_sub,AL,VibSp) 
 allocate(VibSp_t(NVib_sub,3*NAtm_sub))
 call Transp(3*NAtm_sub,NVib_sub,VibSp,VibSp_t) 


!test print for matlab script - to replace b1,b2
 Do J=1,3
    !PRINT '(9E17.8)',(VibSp(i,J),i=1,9)
    !PRINT '(9E17.8)',(VibSp_t(J,i),i=1,9)
 End Do 
 
 allocate(VibSp_full(3*NAtm,NVib_sub))
 call aclear(3*NAtm*NVib_sub,VibSp_full)

 !map from VibSp to VibSp_full
 Do I=1,NAtm_sub
    IPos1 = 3*(subsystem_idx(I)-1)+1
    IPos2 = 3*(subsystem_idx(I)-1)+2
    IPos3 = 3*(subsystem_idx(I)-1)+3
    Do J=1,NVib_sub
       VibSp_full(IPos1,J) = VibSp(3*(I-1)+1,J)
       VibSp_full(IPos2,J) = VibSp(3*(I-1)+2,J)
       VibSp_full(IPos3,J) = VibSp(3*(I-1)+3,J) 
    End Do
 End Do
 allocate(VibSp_full_t(NVib_sub,3*NAtm))
 call Transp(3*NAtm,NVib_sub,VibSp_full,VibSp_full_t)

 Do J=1,3
     !PRINT '(18E17.8)',(VibSp_full(i,J),i=1,18)
     !PRINT '(18E17.8)',(VibSp_full_t(J,i),i=1,18)
 End Do 
 

 !moore-penrose inverse of full hessian
 allocate(FFx_inv(3*NAtm*3*NAtm),Scra(3*NAtm*3*NAtm),Scrb(3*NAtm),Scrc(2*3*NAtm*3*NAtm)) 

 call GInvM(.false.,-1,3*NAtm,FFx,FFx_inv,Scra,Scrb,Scrc)
 !PRINT '(5F18.5)',(FFx_inv(i),i=1,18) ! to check 'pinv(f)' result
 !
 call aclear(3*NAtm*3*NAtm,Scrc)
 call MPACMF(FFx_inv,FFx,Scrc,3*NAtm,3*NAtm,1)
 !PRINT '(5E10.2)',(Scrc(i),i=1,18) ! to check unit matrix 
 
 !obtain FFx_sub
 ! VibSp_full_t (NVib_sub, 3N)  * FFx_inv(3N,3N) -> compl_m1(NVib_sub,3N)
 allocate(compl_m1(NVib_sub,3*NAtm))
 call MMpyMF(NVib_sub,3*NAtm,3*NAtm,VibSp_full_t,FFx_inv,compl_m1)  
 ! compl_m1(NVib_sub,3N) * VibSp_full(3N,NVib_sub) -> compl_m2(NVib_sub,NVib_sub)
 allocate(compl_m2(NVib_sub,NVib_sub))
 call MMpyMF(NVib_sub,3*NAtm,NVib_sub,compl_m1,VibSp_full,compl_m2)
 Do J=1,3
     !PRINT '(3E17.8)',(compl_m2(i,J),i=1,3)
 End Do
 ! compl_m2**-1 -> compl_inv
 allocate(compl_inv(NVib_sub,NVib_sub))
 call GInvM(.false.,-1,NVib_sub,compl_m2,compl_inv,Scra,Scrb,Scrc)  
 Do J=1,3
     !PRINT '(3E17.8)',(compl_inv(i,J),i=1,3)
 End Do
 ! VibSp(3*NAtm_sub,NVib_sub) * compl_inv(NVib_sub,NVib_sub) -> newHes_1(3*NAtm_sub,NVib_sub)
 allocate(newHes_1(3*NAtm_sub,NVib_sub))  
 call MMpyMF(3*NAtm_sub,NVib_sub,NVib_sub,VibSp,compl_inv,newHes_1)
 ! newHes_1(3*NAtm_sub,NVib_sub) * VibSp_t(NVib_sub,3*NAtm_sub) -> FFx_sub(3*NAtm_sub,3*NAtm_sub)
 allocate(FFx_sub(3*NAtm_sub,3*NAtm_sub))
 call MMpyMF(3*NAtm_sub,NVib_sub,3*NAtm_sub,newHes_1,VibSp_t,FFx_sub) 
 Do J=1,3
     !PRINT '(9E17.8)',(FFx_sub(i,J),i=1,9)
 End Do

 !obtain freq and normal mode - vibrational analysis with correct atomic mass 
 call acopy(NAtm_sub,AMass_sub_bak,AMass_sub)
 call MassCent(NAtm_sub,AMass_sub,XYZ_sub,Scr1,Scr2)
 call MIner(.True.,NAtm_sub,AMass_sub,Scr1,Scr1(NAtm3_sub+1),Scr2,Scr2(4))
 call TRVec(.True.,NAtm_sub,NTR_sub,Im,AMass_sub,Scr1,AL,Scr1(NAtm3_sub+1),Scr2)
 NVib_sub=NAtm3_sub-NTR_sub
 call GSorth(.True.,.True.,NAtm3_sub,NTR_sub,AL,Scr3)
 
 call VibSEq(.True.,NAtm_sub,NAtm3_sub,NVib_sub,AMass_sub,FFx_sub,AL,Scr1,Scr2,Scr3)
 call RmNoise(NSS_sub,1.0d-8,AL)
 allocate(Rslt(3*NAtm_sub,8)) ! 8 - from main.f90
 Infred = 0
 IRaman = 0 
 call NormFq(iout,Infred,IRaman,NAtm_sub,NAtm3_sub,NVib_sub,IOP(2),AMass_sub,Scra,FFx_sub,Scra,Scra,AL,Rslt,Scr1,Scr2) 

 call PrtNFq_gsva(iout,NAtm_sub,NAtm3_sub,NVib_sub,ZA_sub,subsystem_idx,AL,Rslt)
 !PRINT '(9E17.8)',(Rslt(i),i=1,3)
end 

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Solve Secular equation in Cartesian coordinates
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine SolvSec(iinp,iout,idt0,irep,ireo,iudt,imdn,iloc,igau,Intact,IOP,Infred,IRaman,IGrd,NAtm,NVib,ctmp,AMass,ZA,XYZ,Grd,  &
  FFx,APT,DPol,AL,Rslt,Scr1,Scr2,Scr3,Scr4,Work)
 implicit real(kind=8) (a-h,o-z)
 dimension :: IOP(*),AMass(*),ZA(*),XYZ(*),Grd(*),FFx(*),APT(*),DPol(*),AL(*),Rslt(*),Scr1(*),Scr2(*),Scr3(*),Scr4(*),Work(*)
 character*4 :: PGNAME(2)
 character*100 :: ctmp
 logical :: Intact,IFAtom

 IFAtom = IOP(1) == -1
 NAtm3=3*NAtm
 NSS=NAtm3*NAtm3

 if(IFAtom) then
   NVib = 0
   goto 1000
 end if

 ! center of mass ---> Scr2(1:3); XYZ in CMCS --> Scr1(1:NAtm3)
 call MassCent(NAtm,AMass,XYZ,Scr1,Scr2)

 ! principal moment of inertia; Eigenvector --> Scr1(NAtm3+1:NAtm3+9)
 call MIner(Intact,NAtm,AMass,Scr1,Scr1(NAtm3+1),Scr2,Scr2(4))

 ! generate m.w. vectors of translations and rotations --> AL
 call TRVec(Intact,NAtm,NTR,Im,AMass,Scr1,AL,Scr1(NAtm3+1),Scr2)
 NVib=NAtm3-NTR
 ! Scr1 and Scr2 will be destroyed

 if (IOP(9) == 1) then
 ! read vibrational normal modes
   call RdNmMod(idt0,Intact,NAtm3,NTR,NVib,IOP,AMass,AL,Scr1)
 else
   ! generate m.w. vectors of vibrations by Gram-Schmidt orthogonalization
   call GSorth(Intact,.True.,NAtm3,NTR,AL,Scr3)
        !PRINT '(18F7.2)',(AL(i),i=1+NAtm3*5,NAtm3+NAtm3*5)
   ! construct secular equation in pure. vib. subspace, do diagonalization, and renormalize the mass-unweighted eigenvectors
   if(IOP(12) == 0) then
     call VibSEq(Intact,NAtm,NAtm3,NVib,AMass,FFx,AL,Scr1,Scr2,Scr3)
   else
     write(iout,"(//,' IFSymtz=T : vibrational normal modes will be symmetrized and reordered.')")
     call VibSEqSymm(iout,ireo,Intact,NAtm,NAtm3,NVib,IOP(5),AMass,XYZ,FFx,AL,Scr1,Scr2,Scr3,Scr4,Work(1),Work(1+NSS),ctmp)
   end if
 end if
 call RmNoise(NSS,1.0d-8,AL)

 ! calculate freq. and IR int. of normal modes. The results are saved in Rslt.
 call NormFq(iout,Infred,IRaman,NAtm,NAtm3,NVib,IOP(2),AMass,ZA,FFx,APT,DPol,AL,Rslt,Scr1,Scr2)

 ! symmetry analysis; irreps of normal modes will be saved in file irep
 if(IOP(12) == 0) then
   call symdrv(iout,irep,NAtm,IOP(5),.false.,.true.,XYZ,ZA,AMass,AL,PGNAME)
 else
   call symdrv(iout,irep,NAtm,IOP(5),.false.,.false.,XYZ,ZA,AMass,AL,PGNAME)
   call dumpir(ireo,irep,ctmp)
 end if

 ! print normal modes results saved in Rslt
 call PrtNFq(iout,irep,Infred,IRaman,NAtm,NAtm3,NVib,(IOP(1)==-3),IOP(2),IOP(10),ZA,AL,Rslt,Scr2)

 ! experimental frequencies
 if(IOP(3) == 1) call RdExFq(iinp,iout,irep,Intact,NAtm3,NVib,Rslt,Scr2,Scr3,Scr4,Work,ctmp)

 ! save plain UniMoVib (ALM) data file and/or localmode.dat (part I)
 if(IOP(6) == 1 .or. IOP(8) == 1) call SavALM((IOP(6) == 1),(IOP(8) == 1),iudt,iloc,Infred,IRaman,IGrd,NAtm,NAtm3,NVib,AMass, &
   ZA,XYZ,Grd,FFx,APT,DPol,IOP(3),AL,Rslt,Scr2,Scr3,WORK,Scr4)

 ! save a molden file
 if(IOP(7) == 1) call SavMDN(imdn,NAtm,NAtm3,NVib,ZA,XYZ,IOP(3),AL,Rslt)

 ! save localmode.dat (part II)
 if(IOP(8) == 1) call SavLOC(iloc,irep,NAtm3,IOP(3),Rslt,PGNAME,ctmp)

 ! save gjf for TS optimization
 if(IOP(11)== 1) call SavGJF(igau,NAtm,NAtm3,NVib,IOP(3),Rslt,XYZ,ZA,AMass,FFx,AL,Scr2,Scr3,WORK,Scr4)

 1000  continue
 ! Thermochemistry calculation. Frequencies are saved in Rslt(1:NVib,3/5) in a.u.
 call Thermochem(iinp,iout,Intact,NAtm,NAtm3,NVib,IFAtom,IOP(3),PGNAME,AMass,XYZ,Rslt,Scr1,Scr2,ctmp)

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Read experimental frequencies from input
!
! Reslt(:,1~4): see NormFq
!
! Reslt(:,5): expt. corrected freq.;  Reslt(:,6): expt. corrected (1.0) or not (-1.0)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RdExFq(iinp,iout,irep,Intact,NAtm3,NVib,Reslt,expf,IRNAME,Ifexpf,dfreq,ctmp)
 implicit real(kind=8) (a-h,o-z)
 parameter(au2wn=5140.48714376d0,epsfrq=5.0d-5,one=1.d0)
 real(kind=8) :: Reslt(NAtm3,*),expf(*),dfreq(*)
 character*100 :: ctmp
 character*4 :: IRNAME(NAtm3)
 namelist/ExpFrq/Mode
 logical :: Intact,OK,Degen,imgfrq,Ifexpf(*)

 write(iout,"(//, 1x,36('*'),/, 1x,'***     Frequency Correction     ***',/, 1x,36('*'))")

 rewind(iinp)
 Mode=0
 read(iinp,ExpFrq,end=1000,err=1000)
 if(Mode /= 0 .and. Mode /= 1) Mode=0

 if(Mode == 0)then
   read(iinp,*,err=1010)(expf(i),i=1,NVib)
   Ifexpf(1:NVib) = .true.
   Reslt(:,6)=one
 else
   ! theor. frequencies are saved in Reslt(1:NVib,3) in a.u.
   do i=1,NVib
     expf(i)=Reslt(i,3)*au2wn
   end do
   nfrq=0
   Ifexpf(1:NVib) = .false.
   do while(.true.)
     read(iinp,"(a100)",end=100,err=1010)ctmp

     if(len_trim(ctmp) == 0) goto 100
     read(ctmp,*,err=1010)ifq,freq
     if(ifq < 1 .or. ifq > NVib)then
       write(iout,"(/,2x,'Ifq = ',i4)")ifq
       goto 1020
     end if
     expf(ifq) = freq
     Ifexpf(ifq) = .true.
     Reslt(ifq,6)=one
     nfrq=nfrq+1
     cycle

     100  exit
   end do
   if(nfrq < 1 .or. nfrq > NVib)then
     write(iout,"(/,2x,'NFreq = ',i4)")nfrq
     goto 1030
   end if
 end if

 ! read irreps
 rewind(irep)
 do i=1,NVib
   read(irep,"(a4)")IRNAME(i)
 end do

 ! frequencies of non-degenerate states cannot be degenerate!
 200   OK = .true.
 do i=1,NVib
   Degen = IRNAME(i)(1:1) /= "A" .and. IRNAME(i)(1:1) /= "B" .and. IRNAME(i)(1:2) /= "Sg" .and. IRNAME(i)(1:2) /= "SG"

   do j=i+1,NVib
     if( ( abs( expf(i) - expf(j) ) < epsfrq) .and. ( ( IRNAME(i) /= IRNAME(j) ) .or. (.not. Degen) ) )then
       expf(j) = expf(j) + epsfrq
       OK = .false.
     end if
   end do
 end do
 if(.not. OK) goto 200

 ! check frequencies
 imgfrq = .false.
 do i=1,NVib
   if(expf(i) < -9.0d3 .or. expf(i) > 9.0d3)then
     write(iout,"(/,2x,'Ifq = ',i4,', Freq = ',f16.4)")i,expf(i)
     goto 1040
   end if
   if(expf(i) < 0.0d0 .or. Reslt(i,3) < 0.0d0) imgfrq = .true.
 end do

 ! print corrected vib. frequencies
 if(imgfrq) then
   write(iout,"(/, 1x,65('-'),/, 1x,'No.      Symm    Expt.Freq    Theo.Freq     Re(Corr)     Im(Corr)',/, 1x,65('-'))")
   do i=1,NVib
     call cplxsub(expf(i),Reslt(i,3)*au2wn,dfreq)
     if(Ifexpf(i)) write(iout,"(1x,i5,4x,a4,4(3x,f10.4))")i,IRNAME(i),expf(i),Reslt(i,3)*au2wn,dfreq(1),dfreq(2)
   end do
   write(iout,"(1x,65('-'),/)")
 else
   write(iout,"(/, 1x,52('-'),/, 1x,'No.      Symm    Expt.Freq    Theo.Freq        Corr.',/, 1x,52('-'))")
   do i=1,NVib
     if(Ifexpf(i)) write(iout,"(1x,i5,4x,a4,3(3x,f10.4))")i,IRNAME(i),expf(i),Reslt(i,3)*au2wn,expf(i)-Reslt(i,3)*au2wn
   end do
   write(iout,"(1x,52('-'),/)")
 end if
 write(iout,"(' <<< NOTE >>>',/,' The corrected frequencies will be used in the following analysis.')")

 do i=1,NVib
   ! cm^-1 --> a.u.
   expf(i)=expf(i)/au2wn
   ! save expt. freq in Reslt(:,5), which will be used in thermochemistry
   ! calculation; the size of Reslt is large enough because NAtm3 >= 6
   Reslt(i,5)=expf(i)
 end do

 return
 1000  call XError(Intact,"The $ExpFrq group cannot be read!")
 1010  call XError(Intact,"Please check your input of exp. freq.s!")
 1020  call XError(Intact,"Ifq in $ExpFrq is out of range!")
 1030  call XError(Intact,"#Freq in $ExpFrq is out of range!")
 1040  call XError(Intact,"Frequency in $ExpFrq is out of range!")
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! calculate normal frequencies, IR intensities, and Raman activities:
!
! Reslt(:,1): k;  Reslt(:,2): mr;  Reslt(:,3): theor. freq;  Reslt(:,4): I.R. Int
!
! Reslt(:,5): expt. corrected freq if Reslt(:,6) = 1.0; see RdExFq.
!
! Reslt(:,7): Raman scattering activity; Reslt(:,8): Depolarization ratio
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine NormFq(iout,Infred,IRaman,NAtm,NAtm3,NVib,IPrint,AMass,ZA,FFx,APT,DPol,AL,Reslt,Scr1,Scr2)
 implicit real(kind=8) (a-h,o-z)
 real(kind=8) :: AMass(*),ZA(*),FFX(NAtm3,NAtm3),APT(3,*),DPol(6,*),AL(NAtm3,*),Reslt(NAtm3,*),Scr1(NAtm3,NAtm3),Scr2(NAtm3,NAtm3)
 parameter(Tol=1.0d-5)

 NMod=NVib
 if(IPrint == 1)NMod=NAtm3

 ! k --> Reslt(:,1)
 do i=1,NMod
   call MMpyMF(1,NAtm3,NAtm3,AL(1,i),FFx,Scr1)
   Reslt(i,1) = dotx(NAtm3,Scr1,AL(1,i))
 end do
 ! mr --> Reslt(:,2)
 call AClear(NAtm3*NAtm3,Scr2)
 do i=1,NAtm3
   im=(i-1)/3+1
   Scr2(i,i)=AMass(im)
 end do
 do i=1,NMod
   call MMpyMF(1,NAtm3,NAtm3,AL(1,i),Scr2,Scr1)
   Reslt(i,2) = dotx(NAtm3,Scr1,AL(1,i))
 end do
 ! omega --> Reslt(:,3)
 do i=1,NMod
   Reslt(i,3) = Reslt(i,1)/Reslt(i,2)
   Reslt(i,3) = sign(sqrt(abs(Reslt(i,3))), Reslt(i,3))
 end do

 ! I.R. Int --> Reslt(:,4)
 if(Infred == 1) then
 ! APT^T * APT --> Scr2
   call Transp(3,NAtm3,APT,Scr1)
   call MMpyMF(NAtm3,3,NAtm3,Scr1,APT,Scr2)
   do i=1,NMod
     call MMpyMF(1,NAtm3,NAtm3,AL(1,i),Scr2,Scr1)
     Reslt(i,4) = dotx(NAtm3,Scr1,AL(1,i)) / Reslt(i,2)
     Reslt(i,4) = abs(Reslt(i,4))
   end do
 end if

 Reslt(:,6)=-1.0d0

 ! Raman --> Reslt(:,7:8)
 if(IRaman == 1) then
   call Transp(6,NAtm3,DPol,Scr1)
   do i=1,NMod
     X = 1.0d0/sqrt(Reslt(i,2))
     do j=1,6
       Scr2(j,1)=dotx(NAtm3,Scr1(1,j),AL(1,i))
     end do
     call AScale(6,X,Scr2(1,1),Scr2(1,1))
 !   [ 1  2  4]
 !   [ 2  3  5]
 !   [ 4  5  6]
 !   alpha^2 --> Scr2(1,2)
     Scr2(1,2) = (Scr2(1,1) + Scr2(3,1) + Scr2(6,1)) / 3.0d0
     Scr2(1,2) = Scr2(1,2) * Scr2(1,2)
 !   beta^2 --> Scr2(4,2)
     Scr2(2,2) = (Scr2(1,1) - Scr2(3,1))**2 + (Scr2(3,1) - Scr2(6,1))**2 + (Scr2(6,1) - Scr2(1,1))**2
     Scr2(3,2) = Scr2(2,1)*Scr2(2,1) + Scr2(4,1)*Scr2(4,1) + Scr2(5,1)*Scr2(5,1)
     Scr2(4,2) = Scr2(2,2) * 0.5d0 + Scr2(3,2) * 3.0d0
 !   Raman scattering activity in a.u.
     Reslt(i,7) = 45.0d0 * Scr2(1,2) + 7.0d0 * Scr2(4,2)
 !   depolarization ratio
     Reslt(i,8) = 45.0d0 * Scr2(1,2) + 4.0d0 * Scr2(4,2)
     if (Reslt(i,8) > Tol) then
       Reslt(i,8) = 3.0d0 * Scr2(4,2) / Reslt(i,8)
     else
       Reslt(i,8) = 0.0d0
     end if
   end do
 end if

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Gram-Schmidt Orthogonalization
! M: the first M vectors in Vec are non-zero ones.
! The size of Scr is N*N if Reorder is .True. or N if .false.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine GSorth(Intact,Reorder,N,M,Vec,Scr)
 implicit real(kind=8) (a-h,o-z)
 parameter(One=1.0D0,Tol=0.1d0)
 real(kind=8) :: Vec(N,N), Scr(*)
 logical :: Intact,Reorder

 NGen = 0
 IGen = 0
 Do Idx = 1, N + M
   ! we will calculate the NGen-th vector
   NGen = NGen + 1

   ! define the NGen-th vector (maybe non-othogonal)
   if(Idx <= M)then
     ! copy the Idx-th non-zero vector to Scr
     Call ACopy(N,Vec(1,Idx),Scr)
   else
     ! take a generator vector from unit matrix
     ! it maybe fails and leads to a new zero vector, so all the column vectors may be tested.
     Call AClear(N,Scr)
     IGen = IGen + 1
     if(IGen > N) call XError(Intact,"IGen > N in GSorth!")
     Scr(IGen) = One
   end if

   !              i-1 [ <u_j, v_i> ]
   ! u_i = v_i - Sum  [------------] u_j
   !              j=1 [ <u_j, u_j> ]
   ! Since u_j has been normalized,
   !
   ! u_i = v_i - Sum_{j=1}^{i-1} {<u_j, v_i> * u_j}
   !
   do Jdx = 1, NGen-1
     X = -dotx(N,Vec(1,Jdx),Scr)
     Call AccAB(N,X,Vec(1,Jdx),Scr)
   end do

   ! normalization
   X = dotx(N,Scr,Scr)
   X = sqrt(X)

   ! new vector?
   if(X > Tol) then
     X = One/X
     Call AScale(N,X,Scr,Vec(1,NGen))
   else
     NGen = NGen - 1
   end if
   if(NGen == N) exit

 end do
 if(NGen /= N) call XError(Intact,"NGen /= N in GSorth!")

 ! reorder
 if(Reorder)then
   Call ACopy(N*N,Vec,Scr)
   Idx=N*M+1
   Jdx=N-M+1
   Call ACopy(N*(N-M),Scr(Idx),Vec)
   Call ACopy(N*M,Scr,Vec(1,Jdx))
 end if

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! sub. VibSEq with symmetry constraint
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine VibSEqSymm(iout,ireo,Intact,NAtm,NAtm3,NVib,ISyTol,AMass,XYZ,FFx,AL,FF0,AL0,IrrMap,Scr1,Scr2,Scr3,ctmp)
 implicit real(kind=8) (a-h,o-z)
 parameter(One=1.d0)
 real(kind=8) :: AMass(*),XYZ(*),FFX(NAtm3,NAtm3),AL(NAtm3,NAtm3),FF0(NAtm3,NAtm3),AL0(NAtm3,NAtm3),Scr1(NAtm3,NAtm3),  &
                 Scr2(NAtm3,NAtm3),Scr3(NAtm3,NAtm3)
 integer      :: IrrMap(NAtm3,3)
 logical      :: Intact
 character*100 :: ctmp
 character*4, allocatable :: Irreps(:)

 !  --------------
 ! | Model system |
 !  --------------

 AL0 = AL
 Nrt = NAtm3 - NVib

 ! simulated FF0 by atomic masses.
 call DDerNRE(NAtm,AMass,XYZ,Scr2,Scr1(1,1),Scr1(1,3))
 call TrAdd(NAtm3,Scr2,FF0)
 ! call MPACMF(Scr2,Scr2,FF0,NAtm3,NAtm3,2)

 ! FF0(m.w.) --> Scr1
 do i=1,NAtm3
   im=(i-1)/3+1
   do j=1,i-1
     jm=(j-1)/3+1
     Scr1(j,i)=FF0(j,i)/sqrt(AMass(im)*AMass(jm))
     Scr1(i,j)=Scr1(j,i)
   end do
   Scr1(i,i)=FF0(i,i)/AMass(im)
 end do
 ! AL0(vib)^T --> Scr2
 call Transp(NAtm3,NVib,AL0,Scr2)
 ! AL0(vib)^T * FF0(m.w.) * AL0(vib) --> Scr1
 call MMpyMF(NVib,NAtm3,NAtm3,Scr2,Scr1,Scr3)
 call MMpyMF(NVib,NAtm3,NVib,Scr3,AL0,Scr1)
 ! diagonalization; eigenvector --> Scr1
 ! The work space is Scr3, which is large enough.
 call DiagS1(Intact,NVib,Scr1,Scr2,Scr3)

 ! transform eigenvectors (saved in Scr1) into Cart. coordinates: AL0(NAtm3,NVib)*Scr1 --> Scr2
 call MMpyMF(NAtm3,NVib,NVib,AL0,Scr1,Scr2)
 Call ACopy(NAtm3*Nrt,AL0(1,NVib+1),Scr2(1,NVib+1))
 ! replace the mass-weighted Gram-Schmidt vectors by the m.w. model vectors
 Call ACopy(NAtm3*NVib,Scr2,AL)

 ! mass-unweighting and renormalization --> AL0
 do i=1,NAtm3
   do j=1,NAtm3
     jm=(j-1)/3+1
     Scr2(j,i) = Scr2(j,i)/sqrt(AMass(jm))
   end do
   X = dotx(NAtm3,Scr2(1,i),Scr2(1,i))
   X = One/sqrt(X)
   call AScale(NAtm3,X,Scr2(1,i),AL0(1,i))
 end do

 ! symmetry analysis; irreps of normal modes will be saved in ireo
 open(ireo,file='SYMFIL2.rep')
 call symdrv(iout,ireo,NAtm,ISyTol,.true.,.true.,XYZ,AMass,AMass,AL0,ctmp)
 ! reorder m.w. vib. modes in AL according to irreps --> AL0
 IrrMap = 0
 allocate(Irreps(NAtm3))
 call CountIrrep(Intact,ireo,NAtm3,NVib,NClass,Irreps,IrrMap)
 !write(iout,"(10i4)")IrrMap(1:NVib,1)
 ivib = 0
 IrrMap(1,2) = 1
 do iclass = 1, NClass
   if(iclass > 1) IrrMap(iclass,2) = 1 + IrrMap(iclass-1,3)
   do i=1,NVib
     if(IrrMap(i,1) == iclass) then
       ivib = ivib + 1
       AL0(:,ivib) = AL(:,i)
     end if
   end do
   IrrMap(iclass,3) = ivib
   ! write(iout,"(i6,4x,a4)") iclass, Irreps(iclass)
 end do
 !write(iout,"(10i4)")IrrMap(1:NClass,2)
 !write(iout,"(10i4)")IrrMap(1:NClass,3)

 !  -------------
 ! | Real system |
 !  -------------

 ! FFX(m.w.) --> FF0
 do i=1,NAtm3
   im=(i-1)/3+1
   do j=1,i-1
     jm=(j-1)/3+1
     FF0(j,i)=FFx(j,i)/sqrt(AMass(im)*AMass(jm))
     FF0(i,j)=FF0(j,i)
   end do
   FF0(i,i)=FFx(i,i)/AMass(im)
 end do

 do iclass = 1, NClass
   ! AL0 --> AL for iclass
   ivib = IrrMap(iclass,2)
   NVib0 = IrrMap(iclass,3) - ivib + 1
   if(NVib0 < 2) then
     Call ACopy(NAtm3,AL0(1,ivib),AL(1,ivib))
   else
     ! AL0(vib)^T --> Scr2
     call Transp(NAtm3,NVib0,AL0(1,ivib),Scr2)
     ! AL0(vib)^T * FFX(m.w.) * AL0(vib) --> Scr1
     call MMpyMF(NVib0,NAtm3,NAtm3,Scr2,FF0,Scr3)
     call MMpyMF(NVib0,NAtm3,NVib0,Scr3,AL0(1,ivib),Scr1)
     ! diagonalization; eigenvector --> Scr1
     ! The work space is Scr3, which is large enough.
     call DiagS1(Intact,NVib0,Scr1,Scr2,Scr3)
     ! transform eigenvectors (saved in Scr1) into Cart. coordinates: AL0(NAtm3,NVib0)*Scr1 --> AL(NAtm3,NVib0)
     call MMpyMF(NAtm3,NVib0,NVib0,AL0(1,ivib),Scr1,AL(1,ivib))
   end if
 end do

 ! reorder vib. normal modes according to the lowest eigenvalues of each class --> AL0
 ! NOTE: symdrv doesn't work well if the degenerate irreps are not continuous after reordering.
 do i=1,NVib
   call MMpyMF(1,NAtm3,NAtm3,AL(1,i),FF0,Scr1)
   Scr3(i,1) = dotx(NAtm3,Scr1,AL(1,i))
 end do
 ! lowest eigenvalues of each class
 do iclass = 1, NClass
   Scr3(iclass,2) = Scr3(IrrMap(iclass,2),1)
 end do
 ! save a data file with reordered irreps
 rewind(ireo)
 write(ireo,"(2i6)") NClass, Nrt
 Emax = ArMax(NClass,jclass,Scr3(1,2)) * 2.0d0
 ivib = 1
 do iclass = 1, NClass
   Emin = ArMin(NClass,jclass,Scr3(1,2))
   Scr3(jclass,2) = Emax
   jvib = IrrMap(jclass,2)
   NVib0 = IrrMap(jclass,3) - jvib + 1
   write(ireo,"(i6,1x,a4)") NVib0, Irreps(jclass)
   Call ACopy(NAtm3*NVib0,AL(1,jvib),AL0(1,ivib))
   ivib = ivib + NVib0
 end do
 Call ACopy(NAtm3*Nrt,AL(1,NVib+1),AL0(1,NVib+1))
 write(ireo,"(a4)") Irreps(NVib+1:NAtm3)

 ! mass-unweighting and renormalization
 do i=1,NAtm3
   do j=1,NAtm3
     jm=(j-1)/3+1
     Scr2(j,i) = AL0(j,i)/sqrt(AMass(jm))
   end do
   X = dotx(NAtm3,Scr2(1,i),Scr2(1,i))
   X = One/sqrt(X)
   call AScale(NAtm3,X,Scr2(1,i),AL(1,i))
 end do

 deallocate(Irreps)

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! construct secular equation in pure. vib. subspace, do diagonalization, and renormalize the mass-unweighted eigenvectors
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine VibSEq(Intact,NAtm,NAtm3,NVib,AMass,FFx,AL,Scr1,Scr2,Scr3)
 implicit real(kind=8) (a-h,o-z)
 parameter(One=1.d0)
 real(kind=8) :: AMass(*),FFX(NAtm3,NAtm3),AL(NAtm3,*),Scr1(NAtm3,NAtm3),Scr2(NAtm3,NAtm3),Scr3(NAtm3,NAtm3)
 logical :: Intact

 ! FFX(m.w.) --> Scr1
 do i=1,NAtm3
   im=(i-1)/3+1
   do j=1,i-1
     jm=(j-1)/3+1
     Scr1(j,i)=FFx(j,i)/sqrt(AMass(im)*AMass(jm))
     Scr1(i,j)=Scr1(j,i)
   end do
   Scr1(i,i)=FFx(i,i)/AMass(im)
 end do
 ! AL(vib)^T --> Scr2
 call Transp(NAtm3,NVib,AL,Scr2)
 ! AL(vib)^T * FFX(m.w.) * AL(vib) --> Scr1
 call MMpyMF(NVib,NAtm3,NAtm3,Scr2,Scr1,Scr3)
 call MMpyMF(NVib,NAtm3,NVib,Scr3,AL,Scr1)
 ! diagonalization; eigenvector --> Scr1
 ! The work space is Scr3, which is large enough.
 call DiagS1(Intact,NVib,Scr1,Scr2,Scr3)

 ! transform eigenvectors (saved in Scr1) into Cart. coordinates: AL(NAtm3,NVib)*Scr1 --> Scr2
 call MMpyMF(NAtm3,NVib,NVib,AL,Scr1,Scr2)
 Call ACopy(NAtm3*(NAtm3-NVib),AL(1,NVib+1),Scr2(1,NVib+1))

 ! mass-unweighting and renormalization
 do i=1,NAtm3
   do j=1,NAtm3
     jm=(j-1)/3+1
     Scr2(j,i) = Scr2(j,i)/sqrt(AMass(jm))
   end do
   X = dotx(NAtm3,Scr2(1,i),Scr2(1,i))
   X = One/sqrt(X)
   call AScale(NAtm3,X,Scr2(1,i),AL(1,i))
 end do

 return
end


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! print results of normal modes for GSVA
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine PrtNFq_gsva(iout,NAtm,NAtm3,NVib,ZA,subsystem_idx,AL,Reslt)
implicit real(kind=8) (a-h,o-z)
parameter(au2wn=5140.48714376d0,au2dy=15.56893d0,cf=31.22307d0,au2ang4=0.52917720859d0**4)
real(kind=8) :: AL(3,NAtm,*),Reslt(NAtm3,*),ZA(*)
integer :: subsystem_idx(*)

write(iout,"(//, 1x,65('*'),/, 1x,'***  Properties of Intrinsic Fragmental Vibrations from GSVA  ***',/, 1x,65('*'))")

write(iout,"(/, ' Results of vibrations:',/,' Normal frequencies (cm^-1), reduced masses (AMU), force constants (mDyn/A)')")

 if(1 == 1)then
   NLine=(NVib-1)/3+1
   do i=1,NLine
     iv1=(i-1)*3+1
     iv2=min(i*3,Nvib)
     ncol=iv2-iv1+1
     write(iout,"(/,18x,3i34)")(j,j=iv1,iv2)
     !write(iout,"('          Irreps',2x,3(30x,a4))")(trim(IRNAME(j)),j=iv1,iv2)
     write(iout,"('     Frequencies',2x,3f34.4)")(Reslt(j,3)*au2wn,j=iv1,iv2)
     write(iout,"('  Reduced masses',2x,3f34.4)")(Reslt(j,2),j=iv1,iv2)
     write(iout,"(' Force constants',2x,3f34.4)")(Reslt(j,1)*au2dy,j=iv1,iv2)
     !if(Infred == 1) then
     !  write(iout,"('  IR intensities',2x,3f34.4)")(Reslt(j,4)*cf*cf,j=iv1,iv2)
     !end if
     !if(IRaman == 1) then
     !  write(iout,"(' Raman sc. activ',2x,3f34.4)")(Reslt(j,7)*au2ang4,j=iv1,iv2)
     !  write(iout,"(' Depolar. ratios',2x,3f34.4)")(Reslt(j,8),j=iv1,iv2)
     !end if
     if(ncol == 1)then
       write(iout,"('        Atom  ZA',2x,  13x,'X         Y         Z') ")
     else if(ncol == 2)then
       write(iout,"('        Atom  ZA',2x,2(13x,'X         Y         Z'))")
     else
       write(iout,"('        Atom  ZA',2x,3(13x,'X         Y         Z'))")
     end if
     do ia=1,NAtm
       write(iout,"(8x,2i4,2x,3(4x,3f10.5))")subsystem_idx(ia),nint(ZA(ia)),((AL(ix,ia,j),ix=1,3),j=iv1,iv2)
     end do
   end do
 end if



!PRINT '(3E17.8)',(Reslt(j,3)*au2wn,j=1,3) 
!PRINT '(3E17.8)',(ZA(j),j=1,3)

return
end


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! print results of normal modes
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine PrtNFq(iout,irep,Infred,IRaman,NAtm,NAtm3,NVib,IfSim,IPrint,IApprox,ZA,AL,Reslt,IRNAME)
 implicit real(kind=8) (a-h,o-z)
 parameter(au2wn=5140.48714376d0,au2dy=15.56893d0,cf=31.22307d0,au2ang4=0.52917720859d0**4)
 real(kind=8) :: ZA(*),AL(3,NAtm,*),Reslt(NAtm3,*)
 character*4 :: IRNAME(NAtm3)
 logical :: IfSim

 ! read irreps
 rewind(irep)
 do i=1,NAtm3
   read(irep,"(a4)")IRNAME(i)
 end do

 write(iout,"(//, 1x,36('*'),/, 1x,'***  Properties of Normal Modes  ***',/, 1x,36('*'))")

 if(IApprox /= 0) then
    write(iout,"(/, 1x,75('!'),/,  &
      ' !! #Internal Coordinates = ',i6,40x,'!!',/,  &
      ' !! NOTE: approximate Hessian is used, so the results are not meaningful. !!',/, 1x,75('!'))") IApprox
 else if(IfSim) then
    write(iout,"(/, 1x,75('!'),/,  &
      ' !!  NOTE: simulated Hessian is used, so the results are not meaningful.  !!',/, 1x,75('!'))")
 end if

 write(iout,"(/, ' Results of vibrations:',/,' Normal frequencies (cm^-1), reduced masses (AMU), force constants (mDyn/A)')",  &
   advance='NO')
 if(Infred == 1) write(iout,"(', IR intensities (km/mol)')",advance='NO')
 if(IRaman == 1) write(iout,"(',',/,' Raman scattering activities (A^4/AMU), depolarization ratios of Raman scattered light')", &
   advance='NO')
 write(iout,*)

 if(IPrint == 1)then
   NLine=(NVib-1)/3+1
   do i=1,NLine
     iv1=(i-1)*3+1
     iv2=min(i*3,Nvib)
     ncol=iv2-iv1+1
     write(iout,"(/,18x,3i34)")(j,j=iv1,iv2)
     write(iout,"('          Irreps',2x,3(30x,a4))")(trim(IRNAME(j)),j=iv1,iv2)
     write(iout,"('     Frequencies',2x,3f34.4)")(Reslt(j,3)*au2wn,j=iv1,iv2)
     write(iout,"('  Reduced masses',2x,3f34.4)")(Reslt(j,2),j=iv1,iv2)
     write(iout,"(' Force constants',2x,3f34.4)")(Reslt(j,1)*au2dy,j=iv1,iv2)
     if(Infred == 1) then
       write(iout,"('  IR intensities',2x,3f34.4)")(Reslt(j,4)*cf*cf,j=iv1,iv2)
     end if
     if(IRaman == 1) then
       write(iout,"(' Raman sc. activ',2x,3f34.4)")(Reslt(j,7)*au2ang4,j=iv1,iv2)
       write(iout,"(' Depolar. ratios',2x,3f34.4)")(Reslt(j,8),j=iv1,iv2)
     end if
     if(ncol == 1)then
       write(iout,"('        Atom  ZA',2x,  13x,'X         Y         Z') ")
     else if(ncol == 2)then
       write(iout,"('        Atom  ZA',2x,2(13x,'X         Y         Z'))")
     else
       write(iout,"('        Atom  ZA',2x,3(13x,'X         Y         Z'))")
     end if
     do ia=1,NAtm
       write(iout,"(8x,2i4,2x,3(4x,3f10.5))")ia,nint(ZA(ia)),((AL(ix,ia,j),ix=1,3),j=iv1,iv2)
     end do
   end do

   write(iout,"(/,' Results of translations and rotations:')")

   NLine=2
   do i=1,NLine
     iv1=Nvib+(i-1)*3+1
     iv2=Nvib+min(i*3,NAtm3-Nvib)
     ncol=iv2-iv1+1
     write(iout,"(/,18x,3i34)")(j-Nvib,j=iv1,iv2)
     write(iout,"('          Irreps',2x,3(30x,a4))")(trim(IRNAME(j)),j=iv1,iv2)
     write(iout,"('     Frequencies',2x,3f34.4)")(Reslt(j,3)*au2wn,j=iv1,iv2)
     write(iout,"('  Reduced masses',2x,3f34.4)")(Reslt(j,2),j=iv1,iv2)
     write(iout,"(' Force constants',2x,3f34.4)")(Reslt(j,1)*au2dy,j=iv1,iv2)
     if(Infred == 1) then
       write(iout,"('  IR intensities',2x,3f34.4)")(Reslt(j,4)*cf*cf,j=iv1,iv2)
     end if
     if(IRaman == 1) then
       write(iout,"(' Raman sc. activ',2x,3f34.4)")(Reslt(j,7)*au2ang4,j=iv1,iv2)
       write(iout,"(' Depolar. ratios',2x,3f34.4)")(Reslt(j,8),j=iv1,iv2)
     end if
     if(ncol == 2)then
         write(iout,"('        Atom  ZA',2x,2(13x,'X         Y         Z'))")
     else
         write(iout,"('        Atom  ZA',2x,3(13x,'X         Y         Z'))")
     end if
     do ia=1,NAtm
       write(iout,"(8x,2i4,2x,3(4x,3f10.5))")ia,nint(ZA(ia)),((AL(ix,ia,j),ix=1,3),j=iv1,iv2)
     end do
   end do

 else
   NLine=(NVib-1)/10+1
   do i=1,NLine
     iv1=(i-1)*10+1
     iv2=min(i*10,Nvib)
     ncol=iv2-iv1+1
     write(iout,"(/,20x,10i10)")(j,j=iv1,iv2)
     write(iout,"('          Irreps',4x,10(6x,a4))")(trim(IRNAME(j)),j=iv1,iv2)
     write(iout,"('     Frequencies',4x,10f10.2)")(Reslt(j,3)*au2wn,j=iv1,iv2)
 !    write(iout,"('  Reduced masses',4x,10f10.4)")(Reslt(j,2),j=iv1,iv2)
     write(iout,"(' Force constants',4x,10f10.4)")(Reslt(j,1)*au2dy,j=iv1,iv2)
     if(Infred == 1) then
       write(iout,"('  IR intensities',4x,10f10.4)")(Reslt(j,4)*cf*cf,j=iv1,iv2)
     end if
     if(IRaman == 1) then
       write(iout,"(' Raman sc. activ',4x,10f10.4)")(Reslt(j,7)*au2ang4,j=iv1,iv2)
       write(iout,"(' Depolar. ratios',4x,10f10.4)")(Reslt(j,8),j=iv1,iv2)
     end if
   end do
 end if

 return
end

!--- END
