!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! driver of symgrp + virrep.
!
! input:
!   iout   = port of output file
!   irep   = port of irrep. file
!   Natom  = #atoms
!   ITol   = (MN) symmetry tolerance M * 10^N * 0.001
!   IFModel= (logical) model system or not
!   IFVMOD = (logical) analyze vib. normal coordinates or not
!   coord  = cartesian coordinates of atoms in Angstrom
!   za     = atomic nuclear charge numbers
!   amass  = atomic masses
!   vmod   = mass-unweighted normal coordinates (vib + tran + rot)
!
! output:
!   PGNAME = mass-independent and dependent point group names
!   The irrep. names of the mass-dependent point group are saved in irep
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine symdrv(iout,irep,Natom,ITol,IFModel,IFVMOD,coord,za,amass,vmod,PGNAME)
 implicit real(kind=8) (a-h,o-z)
 real(kind=8) :: coord(3*Natom),za(Natom),amass(Natom),vmod(*)
 parameter(NOper=20,NGPS=57,NTBS=38)
 character*4 :: PGNAME(2)
 logical :: IFModel,IFVMOD
 real(kind=8),allocatable :: AScr(:)
 character*1,allocatable :: CScr(:)
 integer,allocatable :: IScr(:)

 N = MOD(ITol,10)
 N = SIGN(N,ITol)
 M = ABS(ITol)/10
 tol = dble(M)*10.d0**N   ! it will be multiplied by 0.001 later.

 NATM3 = Natom*3
 NS    = max(9,NATM3)

! IZ    = 1
! IAM   = IZ    + Natom
! IXYZ  = IAM   + Natom
! IIEM  = IXYZ  + NATM3
! IELM  = IIEM  + intwsp(NOper)
! IICY  = IELM  + NOper*9
! ICUB  = IICY  + intwsp(6)
! IROT  = ICUB  + 9
! IEIG  = IROT  + 9
! IEND  = IEIG  + 3

 IZ    = 1
 IAM   = IZ    + Natom
 IXYZ  = IAM   + Natom
 IELM  = IXYZ  + NATM3
 ICUB  = IELM  + NOper*9
 IROT  = ICUB  + 9
 IEIG  = IROT  + 9
 IENA  = IEIG  + 3

 IIEM  = 1
 IICY  = IIEM  + NOper
 IENI  = IICY  + 6

 IENC  = 1

! NOPE  = IEND
! NREP  = NOPE  + intwsp(NGPS)
! NTAB  = NREP  + intwsp(NGPS)
! NALG  = NTAB  + intwsp(NTBS)
! NALO  = NALG  + intwsp(764)
! IARP  = NALO  + intwsp(348)
! IEND  = IARP  + ichawsp(4*406)

 NOPE  = IENI
 NREP  = NOPE  + NGPS
 NTAB  = NREP  + NGPS
 NALG  = NTAB  + NTBS
 NALO  = NALG  + 764
 IENI  = NALO  + 348

 IARP  = IENC
 IENC  = IARP  + 4*406

! for vib. mode analysis

! IJEM  = IEND
! ICMT  = IJEM  + intwsp(NOper*Natom)
! ITCH  = ICMT  + NOper*NATM3
! IGRP  = ITCH  + NOper
! IJX   = IGRP  + NOper*5
! IJY   = IJX   + ichawsp(NOper*4)
! IRRP  = IJY   + intwsp(6)
! IEND  = IRRP  + ichawsp(NATM3*4)

 ICMT  = IENA
 ITCH  = ICMT  + NOper*NATM3
 IGRP  = ITCH  + NOper
 IENA  = IGRP  + NOper*5

 IJEM  = IENI
 IJY   = IJEM  + NOper*Natom
 IENI  = IJY   + 6

 IJX   = IENC
 IRRP  = IJX   + NOper*4
 IENC  = IRRP  + NATM3*4

! scratch

 ITMP  = IENI
 IENI  = ITMP  + Natom

 ISC1  = IENA
 ISC2  = ISC1  + NS
 IENA  = ISC2  + NS

 allocate(AScr(IENA),IScr(IENI),CScr(IENC))

 i1=1
 if(IFModel) i1=2

 do igp=i1,2
   ! 1: Point group for Cartesian coordinates + ZA
   !    this point group symmetry is used for electronic states and orbitals
   ! 2: Point group for Cartesian coordinates + ZA + M
   !    this point group symmetry is used for vib. modes
   AScr = 0.0d0
   IScr = 0
   CScr = ' '

   call ACopy(Natom,za,AScr(IZ))
   if(igp == 1) then
     call ACopy(Natom,AScr(IZ),AScr(IAM))
   else
     call ACopy(Natom,amass,AScr(IAM))
   end if
   call ACopy(NATM3,coord,AScr(IXYZ))
   call symini(Natom,NOper,PGNAME(igp),CScr(IRRP),IScr(IIEM),AScr(IELM),IScr(IJEM),AScr(ICUB))
   call pgsini(IScr(NTAB),IScr(NALG),IScr(NALO),CScr(IARP))
   call symgrp(Natom,NOper,tol,AScr(IZ),AScr(IAM),AScr(IXYZ),IScr(IIEM),AScr(IELM),IScr(IJEM),AScr(ICUB), &
    AScr(IROT),IScr(IICY),AScr(IEIG),IScr(NOPE),IScr(NREP),IScr(NTAB),IScr(NALG),IScr(NALO),CScr(IARP),   &
    nclass,nirred,PGNAME(igp),CScr(IJX),IScr(IJY),AScr(IGRP),IScr(ITMP),AScr(ISC1),AScr(ISC2))

   if(IFModel) then
     if(PGNAME(igp) == "****") then
       call XError(.false.,"The symmetry of model system is not supported.")
     else if(PGNAME(igp) == "C1  ") then
       call XError(.false.,"Model system with C1 symmetry doesn't make sense!")
     end if
   else
     if(igp == 1) then
       write(iout,"(//,1x,'<<< POINT GROUP SYMMETRY >>>')")
       write(iout,"(' Electronic Wavefunctions      :  ',A4)") PGNAME(igp)
     else
       write(iout,"(' Nuclear & Total Wavefunctions :  ',A4)") PGNAME(igp)
     end if
   end if
 end do

 if(.NOT. IFVMOD) return

 ! analyze the irreps of vib. mode for the point group PGNAME(2).
 ! the irrep. symbols will be saved in file irep.
 call virrep(iout,irep,Natom,NOper,vmod,AScr(IZ),AScr(IAM),AScr(IXYZ),IScr(IIEM),AScr(IELM),IScr(IJEM),AScr(ICUB), &
   AScr(IROT),nclass,nirred,CScr(IJX),IScr(IJY),AScr(IGRP),AScr(ICMT),AScr(ITCH),CScr(IRRP),AScr(ISC1),AScr(ISC2))

 deallocate(AScr,IScr,CScr)

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Determine the irrep of a vib. mode for the given point group.
!
! The source codes are taken from MAKOPR and SYMOIR in MOPAC 7.1, which is fully in the public domain.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine virrep(iout,irep,Natom,NOper,vmod,ZA,AMS,XYZ,IELEM,ELEM,JELEM,CUB,ROT,  nclass,nirred,JX,JY,group,carmat,tchar, &
  IRNAME,SC1,SC2)
 implicit real(kind=8) (a-h,o-z)
 real(kind=8) :: vmod(3*Natom,*),ZA(Natom),AMS(Natom),XYZ(3,Natom),ELEM(3,3,*),CUB(3,3),ROT(3,3),group(*),carmat(*),tchar(*), &
   SC1(*),SC2(*)
 dimension :: IELEM(*),JELEM(Natom,*),JY(*)
 character*4 :: IRNAME(Natom*3),JX(*)
 logical :: Intact

 ! make the symmetry-operations for a given Point-Group
 call makopr(iout,Natom,ZA,AMS,XYZ,IELEM,ELEM,JELEM,CUB,ROT,nclass,  JY,  SC1,SC2)

 ! vibrational analyses
 NATM3 = Natom*3
 call symoir(nclass,nirred,Natom,NATM3,NOper,IRNAME,vmod,ELEM,JELEM,ROT,  JX,group,carmat,tchar,  SC1,SC2)

 do i=1,NATM3
   write(irep,"(a4)")adjustl(IRNAME(i))
 end do

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Calculates the Irreducible Representation of an eigenvector.
!
! GROUP(i,j) is the Character of Irreducible Representation i in Class j.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine symoir(nclass,nirred,Natom,NATM3,NOper,IRNAME,vmod,ELEM,JELEM,ROT,  JX,group,carmat,tchar,  SC1,SC2)
 implicit real(kind=8) (a-h,o-z)
 character*4 :: JX(*),IRNAME(NATM3)
 real(kind=8) :: vmod(NATM3,*),ELEM(3,3,*),ROT(3,3),group(NOper,*),carmat(NATM3,*),tchar(*),  SC1(*),SC2(*)
 dimension :: JELEM(Natom,*)

 if(nclass == 1)then
   do i = 1, NATM3
     IRNAME(i) = jx(1)
   end do
   Return
 end if

 ! Calculate all the characters
 do j = 1, nclass
   do i = 1, NATM3
     carmat(i,j) = charvi(Natom,NATM3,vmod(1,i),j,ELEM(1,1,j),JELEM(1,j),ROT,SC1,SC2)
   end do
 end do

 i = 0
 70  continue
 ik = i + 1
 call AClear(nclass,tchar)
 90  continue
 i = i + 1
 if (i > NATM3) return

 do j=1,nclass
   tchar(j)=tchar(j)+carmat(i,j)
 end do
 if (tchar(1) > 5.1D0) go to 70
 do k=1,nirred
   do j=1,nclass
     check=abs(tchar(j)-group(k,j))
     if(check > 0.1d0) goto 140
   end do
   do j=ik,i
     IRNAME(j)=jx(k)
   end do
   goto 70
 140  continue
 end do
 goto 90

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! characters of Irreducible Representation of a vib. mode
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function charvi(Natom,NATM3,vmod,ioper,E,jelem,R,SC1,SC2)
 implicit real(kind=8) (a-h,o-z)
 Dimension :: jelem(Natom)
 real(kind=8) :: vmod(NATM3),E(3,3),R(3,3),SC1(*),SC2(*)
 real(kind=8) :: h(3),p(3)

 charvi = 1.D0
 if (ioper == 1) return

 call AClear(NATM3,SC1)
 call AClear(NATM3,SC2)

 do iatom = 1, Natom
   jatom = jelem(iatom)
   call ACopy(3,vmod(iatom*3-2),p)
   h(1) = r(1,1)*p(1) + r(2,1)*p(2) + r(3,1)*p(3)
   h(2) = r(1,2)*p(1) + r(2,2)*p(2) + r(3,2)*p(3)
   h(3) = r(1,3)*p(1) + r(2,3)*p(2) + r(3,3)*p(3)
   p(1) = e(1,1)*h(1) + e(1,2)*h(2) + e(1,3)*h(3)
   p(2) = e(2,1)*h(1) + e(2,2)*h(2) + e(2,3)*h(3)
   p(3) = e(3,1)*h(1) + e(3,2)*h(2) + e(3,3)*h(3)
   call ACopy(3,h,SC1(iatom*3-2))
   call ACopy(3,p,SC2(jatom*3-2))
 end do

 ! evaluate the term <psi(i)|R(theta)|psi(i)>
 charvi = dotx(NATM3,SC1,SC2)

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! make the symmetry-operations for a given Point-Group
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine makopr(iout,Natom,ZA,AMS,XYZ,IELEM,ELEM,JELEM,CUB,ROT,nclass,  JY,  SC1,SC2)
 implicit real(kind=8) (a-h,o-z)
 parameter(tol = 0.2d0)
 real(kind=8) :: ZA(Natom),AMS(Natom),XYZ(3,Natom),ELEM(3,3,*),CUB(3,3),ROT(3,3),  SC1(*),SC2(*)
 dimension :: IELEM(*),JELEM(Natom,*),JY(*)

 call rotopr(Natom,1,ROT,XYZ,SC1)
 if (nclass < 2) return

 ! construct the Operations corresponding to the Classes, which are stored in ELEM.
 do i = 2, nclass
   call bldsym(jy(i),CUB,ELEM(1,1,i),SC1)
 end do

 ! Use a more tolerant criterion for recognizing operations because the point-group has already been identified.
 do i = 2, nclass
   ! renew jelem(:,2:)
   call chi(Natom,tol,ZA,AMS,XYZ,ELEM(1,1,i),IELEM(i),jelem(1,i),nunmov,SC1,SC2)
   if(IELEM(i) == 1) cycle
   write(iout,"(/,' IClass = ',i4)")i
   call XError(.false.,"The Point Group doesn't work!")
 end do

 call rotopr(Natom,-1,ROT,XYZ,SC1)

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! constructs the point-group symmetry operations
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine bldsym(ioper,CUB,ELEM,SC1)
 implicit real(kind=8) (a-h,o-z)

 parameter(twopi=6.2831853071796d0,one=1.0d0)
 real(kind=8) :: CUB(3,3),ELEM(3,3),SC1(*)
 dimension :: j(3,20)
 save j
 data j(1:3, 1)/ 1,-1,-1/
 data j(1:3, 2)/-1, 1,-1/
 data j(1:3, 3)/-1,-1, 1/
 data j(1:3, 4)/ 1, 1,-1/
 data j(1:3, 5)/ 1,-1, 1/
 data j(1:3, 6)/-1, 1, 1/
 data j(1:3, 7)/-1,-1,-1/
 data j(1:3, 8)/ 3, 0, 1/
 data j(1:3, 9)/ 4, 0, 1/
 data j(1:3,10)/ 5, 0, 1/
 data j(1:3,11)/ 6, 0, 1/
 data j(1:3,12)/ 7, 0, 1/
 data j(1:3,13)/ 8, 0, 1/
 data j(1:3,14)/ 4, 0,-1/
 data j(1:3,15)/ 6, 0,-1/
 data j(1:3,16)/ 8, 0,-1/
 data j(1:3,17)/10, 0,-1/
 data j(1:3,18)/12, 0,-1/
 data j(1:3,19)/ 5, 0,-1/
 data j(1:3,20)/ 0, 0,-1/

 call AClear(9,ELEM)
 do i = 1, 3
   ELEM(i,i) = dble(j(i,ioper))
 end do
 if (ioper /= 20) then
   if (j(1,ioper) >= 2) then
     angle = twopi/dble(j(1,ioper))
     ELEM(1,1) = cos(angle)
     ELEM(2,2) = ELEM(1,1)
     ELEM(2,1) = sin(angle)
     ELEM(1,2) = -ELEM(2,1)
   endif
   if (ioper == 8 .or. ioper == 15) call mult33(CUB,ELEM,SC1)
 else
   ELEM(1,2) = one
   ELEM(2,1) = one
 endif

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Determine the point group symmetry.
!
! The source codes are taken from MOLSYM in MOPAC 7.1, which is fully in the public domain.
!
! Input
!   XYZ    = cartesian coordinates of atoms in Angstrom
!   ZA     = atomic Z
!   AM     = ZA (for symm of molecules) or atomic mass (for symm of molecular vibrations)
!
! Output
!   PGNAME = name of the Group
!   NCLASS = number of Classes in the Group
!   JY     =
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine symgrp(Natom,NOper,tol,ZA,AMS,XYZ,IELEM,ELEM,JELEM,CUB,ROT,ICYC,EIG,  NOPE,NREP,NTAB,NALLG,NALLOP,ALLREP, &
  nclass,nirred,PGNAME,JX,JY,GROUP, Itmp,SC1,SC2)
 implicit real(kind=8) (a-h,o-z)
 Dimension :: IELEM(NOper),JELEM(Natom,*),ICYC(6),NOPE(*),NREP(*),NTAB(*),NALLG(*),NALLOP(*),JY(*),Itmp(*)
 real(kind=8) :: ZA(Natom),AMS(Natom),XYZ(3,Natom),ELEM(3,3,NOper),CUB(3,3),ROT(3,3),EIG(3),GROUP(NOper,5),SC1(*),SC2(*)
 character*4 :: PGNAME,JX(*),ALLREP(*)
 LOGICAL :: ATOM,LINEAR,CUBIC,DEGEN

 pi=acos(-1.0d0)

 ! calculate the principal moment of inertia using    ZA * c1 + AMS * c2,
 ! where c1 = 1/sqrt(3) and c2 = 1 - c1
 SC1(1) = sqrt(3.d0)
 SC1(2) = SC1(1) - 1.d0
 SC1(3) = 1.d0 / SC1(1)
 call AScale(Natom,SC1(2),AMS,SC2)
 call AAdd(Natom,ZA,SC2,SC2)
 call AScale(Natom,SC1(3),SC2,SC2)

 ! center the molecule
 call MassCent(Natom,SC2,XYZ,XYZ,SC1)

 ! principal moment of inertia
 ! ROT: eigenvectors, EIG: eigenvalues
 call MIner2(.false.,Natom,SC2,XYZ,ROT,EIG,SC1)
 SC1(1)=sum(EIG)/3.d0
 ! define the tolerance
 symtol=1.0d-3
 symtol=symtol*tol

 if(EIG(1) > 0.1d0 .and. SC1(1) > 1.0d1)then
   SC1(1) = 1.0d0/SC1(1)
   call AScale(3,SC1(1),EIG,EIG)
 end if

 ! right-hand coordinate system
 call CrossX(ROT(1,1),ROT(1,2),ROT(1,3))

 ! check special groups
 ATOM=(EIG(3) < symtol)
 !xxx LINEAR=(EIG(2) < symtol)
 LINEAR=(EIG(2) < symtol*1.0d-3)
 CUBIC=((EIG(3)-EIG(1)) < symtol) .and. (.NOT. ATOM)
 DEGEN = .false.

 if(ATOM) call XError(.false.,"Atom calculation is not supported.")

 if(LINEAR) then
   call rotopr(Natom,1,ROT,XYZ,SC1)
   IELEM(20)=1
   goto 2000
 end if

 100   if(.NOT. CUBIC)then
   ! for nonlinear and non-cubic molecule, search two-fold degeneracy and put the two-fold degenerate eigenvalues
   ! in the xy plane (high-symmetry axis should be z)
   if((EIG(3)-EIG(2)) < symtol)then
     call acopy(3,ROT(1,1),SC2)
     call acopy(3,ROT(1,2),ROT(1,1))
     call acopy(3,ROT(1,3),ROT(1,2))
     call acopy(3,SC2,ROT(1,3))
     SC2(1)=EIG(1)
     EIG(1)=EIG(2)
     EIG(2)=EIG(3)
     EIG(3)=SC2(1)
   end if
 else
   ! define z-axis of a cubic system
   call cubzaxis(Natom,CUBIC,symtol,ROT,XYZ,ICYC,SC2,Itmp,SC1,SC2(4))
   if(.NOT. CUBIC) goto 100
   IELEM(19)=1
 end if

 ! if EIG(2) < symtol: quasi-linear geometry (Abelian symmetry)
 if((EIG(2)-EIG(1)) < symtol .and. EIG(2) > symtol) DEGEN = .true.
 call rotopr(Natom,1,ROT,XYZ,SC1)

 ! this part is for non-linear and non-abilian symmetry only
 if(DEGEN)then
   ! Check for the existance of a Cn(8~13) or Sn(14~18)
   MaxCn = 2
   do i = 8, 18
     call chi(Natom,symtol,ZA,AMS,XYZ,ELEM(1,1,i),IELEM(i),jelem(1,i),nunmov,SC1,SC2)
     if (IELEM(i) /= 1 .or. i > 13) cycle
     MaxCn = i - 5
   end do

   ! use two adjacent equivalent atoms (i and k, not on the z axis) to define the x axis
   do i = 1, Natom
     if (abs(XYZ(1,i)) <= symtol .and. abs(XYZ(2,i)) <= symtol) cycle
     SC1(1) = ( XYZ(1,i)*XYZ(1,i) + XYZ(2,i)*XYZ(2,i) )
     SC1(3) = 1.0d4
     k = 0
     do j = i + 1, Natom
       if (abs(ZA(i)-ZA(j)) > symtol) cycle
       if (abs(AMS(i)-AMS(j)) > symtol) cycle
       if (abs(abs(XYZ(3,i))-abs(XYZ(3,j))) > symtol) cycle
       SC1(2) = ( XYZ(1,j)*XYZ(1,j) + XYZ(2,j)*XYZ(2,j) )
       if (abs(SC1(1) - SC1(2)) > symtol) cycle
       ! minimal distance between atoms i and j in the xy plane:
       SC1(2) = (XYZ(1,i)-XYZ(1,j))**2 + (XYZ(2,i)-XYZ(2,j))**2
       if (SC1(2) <= SC1(3)) then
         k = j
         SC1(3) = SC1(2)
       end if
     end do
     exit
   end do
   if (k == 0) then
     ! system does not have a Cn axis!  Go back and treat it as an Abelian system
     DEGEN = .false.
     go to 1000
   end if

   ! center between atoms i and k (P; the factor 0.5 is neglected)
   call Aadd(2,XYZ(1,i),XYZ(1,k),SC1)
   ! in the xy plane, calculate the distance OP
   ! the new x axis is defined by OP.
   SC1(3) = sqrt(SC1(1)*SC1(1)+SC1(2)*SC1(2))
   ! rotate atoms
   sina = SC1(2)/SC1(3)
   cosa = SC1(1)/SC1(3)
   call rotmol(Natom,XYZ,sina,cosa,1,2,ROT,SC1)

   ! is there a Sigma(XZ) plane of symmetry?
   call chi(Natom,symtol,ZA,AMS,XYZ,ELEM(1,1,5),IELEM(5),jelem(1,5),nunmov,SC1,SC2)
   if(IELEM(5) == 0)then
     ! is there a C2(X) axis of rotation?
     call chi(Natom,symtol,ZA,AMS,XYZ,ELEM(1,1,1),IELEM(1),jelem(1,1),nunmov,SC1,SC2)
     if(IELEM(1) == 1)then
       ! check for an improper axis of rotation
       theta = pi/dble(2*MaxCn)
       sina = sin(theta)
       cosa = cos(theta)
       icheck = 0
       180  call rotmol(Natom,XYZ,sina,cosa,1,2,ROT,SC1)
       if (icheck > 0)then
         DEGEN = MaxCn /= 2
         go to 1000
       end if
       call chi(Natom,symtol,ZA,AMS,XYZ,ELEM(1,1,5),IELEM(5),jelem(1,5),nunmov,SC1,SC2)
       if (IELEM(5) == 1) go to 1000
       icheck = 1
       sina = -sina
       go to 180
     end if
   end if
 end if

 1000  continue
 if(CUBIC) call orient(Natom,symtol,ZA,AMS,XYZ,ROT,CUB,ELEM,IELEM,JELEM,SC1,SC2)

 ! this part is for abilian symmetries only
 if(.NOT. DEGEN)then
   ! icyc(i) is the number of atoms unmoved by the operation-i
   do i = 1, 6
     call chi(Natom,symtol,ZA,AMS,XYZ,ELEM(1,1,i),IELEM(i),jelem(1,i),nunmov,SC1,SC2)
     ICYC(i) = (1 + nunmov)*IELEM(i)
   end do

   ! Determine the principal axis of Abelian Groups by atom counts
   naxes = IELEM(1) + IELEM(2) + IELEM(3)
   if (naxes < 1) then
     iz = 3
     if (ICYC(5) > ICYC(4)) iz = 2
     if (ICYC(6) > ICYC(7-iz)) iz = 1
   else if (naxes == 1) then
     if (IELEM(1) == 1) then
       iz = 1
     else if (IELEM(2) == 1) then
       iz = 2
     else if (IELEM(3) == 1) then
       iz = 3
     end if
   else
     iz = 1
     if (ICYC(2) > ICYC(1)) iz = 2
     if (ICYC(3) > ICYC(iz)) iz = 3
   end if
   ICYC(7-iz) = -1
   ix = 1
   if (ICYC(5) > ICYC(6)) ix = 2
   if (ICYC(4) > ICYC(7-ix)) ix = 3
   iy = 6 - ix - iz

   ! re-orient the molecule so that the principal axis is z.
   call newaxis(ix,iy,ROT,SC1)
   call rotopr(Natom,-1,ROT,XYZ,SC2)
   call ACopy(9,SC1,ROT)
   call rotopr(Natom, 1,ROT,XYZ,SC2)
 end if

 2000  continue
 ! re-calculate the first 7 Characters.
 do i = 1, 7
   call chi(Natom,symtol,ZA,AMS,XYZ,ELEM(1,1,i),IELEM(i),jelem(1,i),nunmov,SC1,SC2)
 end do

 call rotopr(Natom,-1,ROT,XYZ,SC2)

 SC2(1) = sum(EIG)
 EIG = SC2(1) - EIG

 call cartab(PGNAME,NOper,NOPE,NREP,NTAB,NALLG,NALLOP,ALLREP,IELEM,nclass,nirred,JX,JY,GROUP)

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! In cubic systems the three moments of inertia are identical. Therefore, use the nearest atom to the center to determine the
! principal axis (Z).
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine orient(Natom,symtol,ZA,AMS,XYZ,ROT,CUB,ELEM,IELEM,JELEM,SC1,SC2)
 implicit real(kind=8) (a-h,o-z)
 parameter(sqhlf=sqrt(0.5d0))
 dimension :: IELEM(*),JELEM(Natom,*)
 real(kind=8) :: ZA(Natom),AMS(Natom),XYZ(3,Natom),ROT(3,3),CUB(3,3),ELEM(3,3,*),SC1(*),SC2(*)
 real(kind=8) :: wink(2)

 wink(1)=acos(sqrt(1.d0/3.d0))
 wink(2)=atan(3.d0-sqrt(5.d0))

 if (ielem(8) == 1) then
   ! Check for a S4 or C5 axis
   do i = 1, 2
     jota = 18 - 4*i
     wink2 = wink(i)
     sina = sin(wink2)
     cosa = cos(wink2)
     call rotmol(Natom,XYZ,sina,cosa,1,3,ROT,SC1)
     call chi(Natom,symtol,ZA,AMS,XYZ,ELEM(1,1,jota),IELEM(jota),jelem(1,jota),nunmov,SC1,SC2)
     if (ielem(jota) == 1) exit
     if (i == 1) then
       call chi(Natom,symtol,ZA,AMS,XYZ,ELEM(1,1,3),IELEM(3),jelem(1,3),nunmov,SC1,SC2)
       if (ielem(3) == 1) exit
     end if
     wink2 = -wink2
     sinb = sin(2.D0*wink2)
     cosb = cos(2.D0*wink2)
     call rotmol(Natom,XYZ,sinb,cosb,1,3,ROT,SC1)
     call chi(Natom,symtol,ZA,AMS,XYZ,ELEM(1,1,jota),IELEM(jota),jelem(1,jota),nunmov,SC1,SC2)
     if (ielem(jota) == 1) exit
     if (i == 1) then
       call chi(Natom,symtol,ZA,AMS,XYZ,ELEM(1,1,3),IELEM(3),jelem(1,3),nunmov,SC1,SC2)
       if (ielem(3) == 1) exit
     end if
     call rotmol(Natom,XYZ,sina,cosa,1,3,ROT,SC1)
   end do
   call chi(Natom,symtol,ZA,AMS,XYZ,ELEM(1,1,9),IELEM(9),jelem(1,9),nunmov,SC1,SC2)
   ! Check on all IELEM registers
   if (ielem(10) == 1) call chi(Natom,symtol,ZA,AMS,XYZ,ELEM(1,1,17),IELEM(17),jelem(1,17),nunmov,SC1,SC2)
 else
   ! No C3 axis, therefore not T, Td, Th, O, or Oh.
   wink2 = -wink(1)
   if (ielem(10) == 1) wink2 = -wink(2)
   sina = -sin(wink2)
   cosa = cos(wink2)
   call rotmol(Natom,XYZ,sina,cosa,1,3,ROT,SC1)
   call chi(Natom,symtol,ZA,AMS,XYZ,ELEM(1,1,8),IELEM(8),jelem(1,8),nunmov,SC1,SC2)
   call rotmol(Natom,XYZ,-sina,cosa,1,3,ROT,SC1)
   if (ielem(8) == 0) then
     if (ielem(9) == 0) then
       wink2 = -wink2
     else
       call rotmol(Natom,XYZ,sqhlf,sqhlf,1,2,ROT,SC1)
     end if
   end if
 end if

 j = sum(ielem(1:17))
 if (j == 2 .and. ielem(1)+ielem(8) == 2) return

 CUB(1,1) = cos(wink2)
 CUB(3,3) = CUB(1,1)
 CUB(1,3) = sin(wink2)
 CUB(3,1) = -CUB(1,3)
 call mult33(CUB,ELEM(1,1,8),SC1)
 call mult33(CUB,ELEM(1,1,15),SC1)
 call chi(Natom,symtol,ZA,AMS,XYZ,ELEM(1,1,8),IELEM(8),jelem(1,8),nunmov,SC1,SC2)
 call chi(Natom,symtol,ZA,AMS,XYZ,ELEM(1,1,15),IELEM(15),jelem(1,15),nunmov,SC1,SC2)

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! It constructs Character Tables for Point Groups
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine cartab(PGNAME,NOper,nope,nrep,ntab,nallg,nallop,allrep,ielem,nclass,nirred,jx,jy,Group)
 implicit real(kind=8) (a-h,o-z)
 parameter(NGPS=57,NTBS=38,One=1.0d0,Two=2.0d0,dpi=6.283185307179d0)
 real(kind=8) :: GROUP(NOper,5)
 dimension :: nope(NGPS),nrep(NGPS),ntab(*),nallg(*),nallop(*),ielem(*),jy(*)
 character*4 :: PGNAME,jx(*),allrep(*)
 logical :: SYMDEC

 ! Generate array to hold addresses of groups
 ! NOPE will apply to operations
 ! NREP will apply to group name and irreducible representations
 nope(1) = 1
 nrep(1) = 1
 do i = 2, NGPS
   nope(i) = nope(i-1) + nallop(nope(i-1)) + 4
   nrep(i) = nrep(i-1) + nallop(nope(i-1)+1) + 1
 end do

 ! Generate array to hold addresses of character tables
 j = 1
 do i = 1, NTBS
   k = ntab(i)
   ntab(i) = j
   j = j + k
 end do

 ! Identify the Point Group of the molecule
 do igroup = NGPS, 1, -1
   i = nallop(nope(igroup)+3)
   if ( SYMDEC(i,ielem) ) exit
 end do

 ! IGROUP:  The number of the point-group of the system.
 ! ISTART:  Starting address of the group.
 ! IGP   :  Starting address of the names used in the group.
 ! PGNAME:  The name of the group.
 ! NCLASS:  Number of operations used to represent the group, plus 1.
 ! NIRRED:  Number of Irreducible Representations.
 ! GROUP :  The full point-group character table for group IGROUP.
 istart = nope(igroup)
 igp = nrep(igroup)
 PGNAME = allrep(igp)
 nclass = nallop(istart) + 1
 nirred = nallop(istart+1)

 do i=1,nclass
   group(1,i)=One
 end do
 do i=1,nirred
   jx(i)=allrep(igp+i)
 end do
 jy(1) = 0
 do i=2,nclass
   jy(i)=nallop(istart+2+i)
 end do
 istart=ntab(nallop(istart+2))-1
 do i=2,nirred
    do k=1,nclass
       istart=istart+1
       buff=dble(nallg(istart))
       if(buff >= 10.d0) then
          nzz=nallg(istart)
          nz=nzz/10
          fz=dble(nz)
          fn=dble(nzz-10*nz)
          buff=Two*cos(dpi*fn/fz)
       endif
     group(i,k)=buff
   end do
 end do

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! SYMDEC matches up a set of symmetry operations for the system with a point-group.  The set of symmetry operations are stored
!   in IELEM.  Point-groups are represented by a number, N1, which is expanded into binary.
!
! A system has the symmetry of a specific point-group if every operation of the point-group is present in the system. Extra
!   operations may also be present, but are ignored. This allows a controlled descent in symmetry.
!   The symmetry groups in the calling routine, CARTAB, are stored in order of symmetry - C1 to R3.
!
! SYMDEC returns a .TRUE. if the system matches point-group N1.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
logical function SYMDEC(n1, ielem)
 implicit real(kind=8) (a-h,o-z)
 dimension :: ielem(20)

 in1 = n1
 do i = 1, 20
   ibin = mod(in1,2)
   if (ielem(i) /= 1 .and. ibin == 1) go to 20
   in1 = in1/2
 end do
 symdec = .TRUE.
 Return

 20    symdec = .FALSE.

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! BLOCK DATA FOR ALL THE POINT-GROUPS USED IN MOPAC
!
! Point-Groups are stored in order of increasing symmetry.  Groups supported are, in order
!
! C1, Cs, Ci, C2, D2, C2v, C2h, D2h, C3, C4, S4, D3, C3v, C3h, C5, C6, S6, C7, C8, S8, C4v, D4, D2d, C5v, C6v, D6,
! C4h, D3h, D3d, C7v, C7h, C8v, D8, D4d, D4h, C5h, D5, D5h, D5d, C6h, D6h, D6d, D7, D7h, D7d, C8h, D8h, T, Td, O,
! Th, Oh, I, Ih, Cv, Dh, R3
!
! Character tables are stored separately from the point groups, so that more than one group can use the same character
! table. The totally symmetric representation is assumed to exist, and has been suppressed from the character tables.
!
! Each point group is represented by two arrays, a I<group symbol> and a R<group symbol> array.  The structure of these
! tables is as follows:
!
!      I<group symbol>
!
! (1):    Number of classes of operation used.
! (2):    Number of Irreducible Representations.
! (3):    The number of the associated character table.
! (4):    A "Magic number" identifying the group.
! (5 on): Numbers indicating the operations of the group
!
!      R<group symbol>
!
! (1):    The name of the group
! (2 on): The names of the irreducible representations.
!
! Character tables are named arbitarily.  Each table has an associated entry in the NTAB array.  NTAB(i) holds the
! number of elements in table i.
!
! Storage of point-groups:  Point-Group Name,Number of Classes, Names of Classes,Number of Irred. Reps., Names of Irred. Reps,
!                           Number of data in point-group, Name of Representation
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine pgsini(ntab,nallg,nallop,allrep)
 implicit real(kind=8) (a-h,o-z)
 character*4 :: allrep(*)
 dimension :: ntab(*), nallg(*), nallop(*)

 ! IG11
 ntab(1) = 1
 nallg(1:1) = (/ 1 /)
 ! IG21
 ntab(2) = 2
 nallg(2:3) = (/ 1, -1 /)
 ! IG42
 ntab(3) = 9
 nallg(4:12) = (/  1,  1, -1,  1, -1,  1,  1, -1, -1 /)
 ! IG83
 ntab(4) = 28
 nallg(13:40) = (/  1,  1, -1,  1,  1, -1,  1,  1,  1, -1, -1,  1,  1,  1,  1, -1,  1,  1, -1, -1,  1, -1,  1, -1,  1, -1, -1,  &
                   -1 /)
 ! IG31
 ntab(5) = 2
 nallg(41:42) = (/ 2, -1 /)
 ! IG41
 ntab(6) = 4
 nallg(43:46) = (/ 1, -1, 2, 0 /)
 ! IG51
 ntab(7) = 4
 nallg(47:50) = (/ 2, 51, 2, 52 /)
 ! IG61
 ntab(8) = 6
 nallg(51:56) = (/ 1, -1, 2, 1, 2, -1 /)
 ! IG71
 ntab(9) = 6
 nallg(57:62) = (/ 2, 71, 2, 72, 2, 73 /)
 ! IG81
 ntab(10) = 8
 nallg(63:70) = (/ 1, -1, 2, 81, 2, 0, 2, 83 /)
 ! IG84
 ntab(11) = 12
 nallg(71:82) = (/ 1, 1, -1, 1, -1, 1, 1, -1, -1, 2, 0, 0 /)
 ! IG52
 ntab(12) = 9
 nallg(83:91) = (/ 1, 1, -1, 2, 51, 0, 2, 52, 0 /)
 ! IG62
 ntab(13) = 20
 nallg(92:111) = (/ 1, 1, 1, -1, 1, 1, -1, 1, 1, 1, -1, -1, 2, -1, 0, -2, 2, -1, 0, 2 /)
 ! IG72
 ntab(14) = 12
 nallg(112:123) = (/ 1, 1, -1, 2, 71, 0, 2, 72, 0, 2, 73, 0 /)
 ! IG82
 ntab(15) = 18
 nallg(124:141) = (/ 1, 1, -1, 1, -1, -1, 1, -1, 1, 2, 81, 0, 2, 0, 0, 2, 83, 0 /)
 ! IG4H
 ntab(16) = 36
 nallg(142:177) = (/ 1, 1, 1, -1, 1, 1, -1, 1, 1, 1, -1, -1, 1, -1, 1, 1, 1, -1, 1, -1, 1, -1, -1, 1, 1, -1, -1, -1, 2, 0, 0, &
                    -2, 2, 0,  0, 2 /)
 ! IGC5H
 ntab(17) = 15
 nallg(178:192) = (/ 1, 1, -1, 2, 51, 2, 2, 51, -2, 2, 52, 2, 2, 52, -2 /)
 ! IGD5G
 ntab(18) = 28
 nallg(193:220) = (/ 1, 1, -1, -1, 1, 1, 1, -1, 1, 1, -1, 1, 2, 51, 2, 0, 2, 51, -2, 0, 2, 52, 2, 0, 2, 52, -2, 0 /)
 ! IGC6H
 ntab(19) = 21
 nallg(221:241) = (/ 1, 1, -1, 1, -1, 1, 1, -1, -1, 2, 1, 2, 2, 1, -2, 2, -1, 2, 2, -1, -2 /)
 ! IGD6H
 ntab(20) = 44
 nallg(242:285) = (/ 1, 1, 1, -1, 1, 1, -1, 1, 1, 1, -1, -1, 1, -1, 1, 1, 1, -1, 1, -1, 1, -1, -1, 1, 1, -1, -1, -1, 2, 1, 0, &
                     2, 2, 1,  0,-2, 2, -1, 0, 2, 2, -1,  0,-2 /)
 ! IGD6D
 ntab(21) = 24
 nallg(286:309) = (/ 1, 1, -1, 1, -1, -1, 1, -1, 1, 2, 121, 0, 2, 1, 0, 2, 0, 0, 2, -1, 0, 2, 125, 0 /)
 ! IGC7H
 ntab(22) = 21
 nallg(310:330) = (/ 1, 1, -1, 2, 71, 2, 2, 71, -2, 2, 72, 2, 2, 72, -2, 2, 73, 2, 2, 73, -2 /)
 ! IGD7H
 ntab(23) = 36
 nallg(331:366) = (/ 1, 1, -1, -1, 1, 1, 1, -1, 1, 1, -1, 1, 2, 71, 2, 0, 2, 71, -2, 0, 2, 72, 2, 0, 2, 72, -2, 0, 2, 73, 2, &
                     0, 2, 73, -2, 0 /)
 ! IGD8H
 ntab(24) = 52
 nallg(367:418) = (/ 1, 1, -1, -1, 1, 1, 1, -1, 1, 1, -1, 1, 1, -1, -1, -1, 1, -1, 1, 1, 1, -1, -1, 1, 1, -1, 1, -1, 2, 81,-2, &
                     0, 2, 81,  2, 0, 2, 0,  2, 0, 2,  0,-2, 0,  2, 83, -2, 0,  2,83, 2, 0 /)
 ! IGC8H
 ntab(25) = 27
 nallg(419:445) = (/ 1, 1, -1, 1, -1, 1, 1, -1, -1, 2, 81, 2, 2, 81, -2, 2, 0, 2, 2, 0, -2, 2, 83, 2, 2, 83, -2 /)
 ! IGT
 ntab(26) = 6
 nallg(446:451) = (/ 2, -1, 2, 3, 0, -1 /)
 ! IGTH
 ntab(27) = 20
 nallg(452:471) = (/ 1, 1, -1, 1, 2, -1, 2, 2, 2, -1, -2, 2, 3, 0, 3, -1, 3, 0, -3, -1 /)
 ! IGTD
 ntab(28) = 16
 nallg(472:487) = (/ 1, -1, 1, 1, 2, 0, -1, 2, 3, -1, 0, -1, 3, 1, 0, -1 /)
 ! IGOH
 ntab(29) = 36
 nallg(488:523) = (/ 1, 1, -1, 1, 1, -1, 1, 1, 1, -1, -1, 1, 2, 0, 2, -1, 2, 0, -2, -1, 3, -1, 3, 0, 3, -1, -3, 0, 3, 1, 3, 0, &
                     3, 1, -3, 0 /)
 ! IGI
 ntab(30) = 12
 nallg(524:535) = (/ 3, -1, 101, 3, -1, 103, 4, 0, -1, 5, 1, 0 /)
 ! IGIH
 ntab(31) = 36
 nallg(536:571) = (/ 1, 1, -1, 1, 3, 101, 3, 0, 3, 101, -3, 0, 3, 103, 3, 0, 3, 103, -3, 0, 4, -1, 4, 1, 4, -1, -4, 1, 5, 0, &
                     5,-1,  5, 0,-5, -1 /)
 ! IGCV
 ntab(32) = 20
 nallg(572:591) = (/ 1, 1, 1, -1, 2, 71, 61, 0, 2, 72, 62, 0, 2, 73, 63, 0, 2, 74, 64, 0 /)
 ! IGDH
 ntab(33) = 55
 nallg(592:646) = (/ 1, 1, 1, -1, 1, 1, 1, 1, 1, -1, 1, 1, 1, -1, -1, 2, 71, 61, 2, 0, 2, 71, 61, -2, 0, 2, 72, 62, 2, 0, 2, &
                    72,62,-2,  0, 2,73,63, 2, 0,  2,73,63,-2,  0,  2,74, 64,  2, 0, 2,74, 64, -2,  0 /)
 ! IGC4H
 ntab(34) = 15
 nallg(647:661) = (/ 1, 1, -1, 1, -1, 1, 1, -1, -1, 2, 0, 2, 2, 0, -2 /)
 ! IGD3
 ntab(35) = 6
 nallg(662:667) = (/ 1, 1, -1, 2, -1, 0 /)
 ! IGC3H
 ntab(36) = 9
 nallg(668:676) = (/ 2, -1, 2, 1, 1, -1, 2, -1, -2 /)
 ! IGO
 ntab(37) = 12
 nallg(677:688) = (/ 1, -1, 1, 2, 0, -1, 3, -1, 0, 3, 1, 0 /)
 ! IGR3
 ntab(38) = 60
 nallg(689:764) = (/ 1, 1, 1, -1, 3, 101, 0, 3, 3, 101, 0, -3, 5, 0, -1, 5, 5, 0, -1, -5, 7, 106, 1, 7, 7, 106, 1, -7, 9, -1, &
                     0, 9, 9, -1, 0, -9, 11, 1, -1, 11, 11, 1, -1, -11, 13, 101, 1, 13, 13, 101, 1, -13, 15, 0, 0, 15, 15, 0, &
                     0, -15, 17, 106, -1, 17, 17, 106, -1, -17, 19, -1, 1, 19, 19, -1, 1, -19 /)

 ! C1
 allrep(1:2) = (/ "C1  ", "A   " /)
 nallop(1:4) = (/ 0, 1, 1, 0 /)
 ! CS
 allrep(3:5) = (/ "Cs  ", "A'  ", 'A"  ' /)
 nallop(5:9) = (/ 1, 2, 2, 8, 4 /)
 ! CI
 allrep(6:8) = (/ "Ci  ", "Ag  ", "Au  " /)
 nallop(10:14) = (/ 1, 2, 2, 64, 7 /)
 ! C2
 allrep(9:11) = (/ "C2  ", "A   ", "B   " /)
 nallop(15:19) = (/ 1, 2, 2, 4, 3 /)
 ! D2
 allrep(12:16) = (/ "D2  ", "A   ", "B1  ", "B2  ", "B3  " /)
 nallop(20:25) = (/ 2, 4, 3, 7, 3, 2 /)
 ! C2V
 allrep(17:21) = (/ "C2v ", "A1  ", "A2  ", "B1  ", "B2  " /)
 nallop(26:31) = (/ 2, 4, 3, 52, 3, 5 /)
 ! C2H
 allrep(22:26) = (/ "C2h ", "Ag  ", "Au  ", "Bg  ", "Bu  " /)
 nallop(32:37) = (/ 2, 4, 3, 76, 3, 7 /)
 ! D2H
 allrep(27:35) = (/ "D2h ", "Ag  ", "B1g ", "B2g ", "B3g ", "Au  ", "B1u ", "B2u ", "B3u " /)
 nallop(38:44) = (/ 3, 8, 4, 127, 3, 2, 7 /)
 ! C3
 allrep(36:38) = (/ "C3  ", "A   ", "E   " /)
 nallop(45:49) = (/ 1, 2, 5, 128, 8 /)
 ! C4
 allrep(39:42) = (/ "C4  ", "A   ", "B   ", "E   " /)
 nallop(50:54) = (/ 1, 3, 6, 260, 9 /)
 ! S4
 allrep(43:46) = (/ "S4  ", "A   ", "B   ", "E   " /)
 nallop(55:59) = (/ 1, 3, 6, 8196, 14 /)
 ! D3
 allrep(47:50) = (/ "D3  ", "A1  ", "A2  ", "E   " /)
 nallop(60:65) = (/ 2, 3, 35, 129, 8, 1 /)
 ! C3V
 allrep(51:54) = (/ "C3v ", "A1  ", "A2  ", "E   " /)
 nallop(66:71) = (/ 2, 3, 35, 144, 8, 5 /)
 ! C3H
 allrep(55:59) = (/ "C3h ", "A'  ", "E'  ", 'A"  ', 'E"  ' /)
 nallop(72:77) = (/ 2, 4, 36, 136, 8, 4 /)
 ! C5
 allrep(60:63) = (/ "C5  ", "A   ", "E1  ", "E2  " /)
 nallop(78:82) = (/ 1, 3, 7, 512, 10 /)
 ! C6
 allrep(64:68) = (/ "C6  ", "A   ", "B   ", "E1  ", "E2  " /)
 nallop(83:87) = (/ 1, 4, 8, 1156, 11 /)
 ! S6
 allrep(69:73) = (/ "S6  ", "Ag  ", "Au  ", "Eu  ", "Eg  " /)
 nallop(88:92) = (/ 1, 4, 8, 16576, 15 /)
 ! C7
 allrep(74:78) = (/ "C7  ", "A   ", "E1  ", "E2  ", "E3  " /)
 nallop(93:97) = (/ 1, 4, 9, 2048, 12 /)
 ! C8
 allrep(79:84) = (/ "C8  ", "A   ", "B   ", "E1  ", "E2  ", "E3  " /)
 nallop(98:102) = (/ 1, 5, 10, 4356, 13 /)
 ! S8
 allrep(85:90) = (/ "S8  ", "A   ", "B   ", "E1  ", "E2  ", "E3  " /)
 nallop(103:107) = (/ 1, 5, 10, 33028, 16 /)
 ! C4V
 allrep(91:96) = (/ "C4v ", "A1  ", "A2  ", "B1  ", "B2  ", "E   " /)
 nallop(108:113) = (/ 2, 5, 11, 308, 9, 5 /)
 ! D4
 allrep(97:102) = (/ "D4  ", "A1  ", "A2  ", "B1  ", "B2  ", "E   " /)
 nallop(114:119) = (/ 2, 5, 11, 263, 9, 1 /)
 ! D2D
 allrep(103:108) = (/ "D2d ", "A1  ", "A2  ", "B2  ", "B1  ", "E   " /)
 nallop(120:125) = (/ 2, 5, 11, 8244, 14, 5 /)
 ! C5V
 allrep(109:113) = (/ "C5v ", "A1  ", "A2  ", "E1  ", "E2  " /)
 nallop(126:131) = (/ 2, 4, 12, 528, 10, 5 /)
 ! C6V
 allrep(114:120) = (/ "C6v ", "A1  ", "B1  ", "A2  ", "B2  ", "E1  ", "E2  " /)
 nallop(132:138) = (/ 3, 6, 13, 1204, 8, 5, 3 /)
 ! D6
 allrep(121:127) = (/ "D6  ", "A1  ", "B1  ", "A2  ", "B2  ", "E1  ", "E2  " /)
 nallop(139:145) = (/ 3, 6, 13, 1159, 8, 1, 3 /)
 ! C4H
 allrep(128:134) = (/ "C4h ", "Ag  ", "Au  ", "Bg  ", "Bu  ", "Eg  ", "Eu  " /)
 nallop(146:151) = (/ 2, 6, 34, 8524, 9, 7 /)
 ! D3H
 allrep(135:141) = (/ "D3h ", "A1' ", 'A1" ', "A2' ", 'A2" ', 'E"  ', "E'  " /)
 nallop(152:158) = (/ 3, 6, 13, 153, 8, 1, 4 /)
 ! D3D
 allrep(142:148) = (/ "D3d ", "A1g ", "A1u ", "A2g ", "A2u ", "Eu  ", "Eg  " /)
 nallop(159:165) = (/ 3, 6, 13, 16594, 8, 2, 7 /)
 ! C7V
 allrep(149:154) = (/ "C7v ", "A1  ", "A2  ", "E1  ", "E2  ", "E3  " /)
 nallop(166:171) = (/ 2, 5, 14, 2064, 12, 5 /)
 ! C7H
 allrep(155:163) = (/ "C7h ", "A'  ", 'A"  ', "E1' ", 'E1" ', "E2' ", 'E2" ', "E3' ", 'E3" ' /)
 nallop(172:177) = (/ 2, 8, 22, 2056, 12, 4 /)
 ! C8V
 allrep(164:171) = (/ "C8v ", "A1  ", "A2  ", "B2  ", "B1  ", "E1  ", "E2  ", "E3  " /)
 nallop(178:183) = (/ 2, 7, 15, 4404, 13, 5 /)
 ! D8
 allrep(172:179) = (/ "D8  ", "A1  ", "A2  ", "B2  ", "B1  ", "E1  ", "E2  ", "E3  " /)
 nallop(184:189) = (/ 2, 7, 15, 4359, 13, 1 /)
 ! D4D
 allrep(180:187) = (/ "D4d ", "A1  ", "A2  ", "B1  ", "B2  ", "E1  ", "E2  ", "E3  " /)
 nallop(190:195) = (/ 2, 7, 15, 33076, 16, 5 /)
 ! D4H
 allrep(188:198) = (/ "D4h ", "A1g ", "A1u ", "A2g ", "A2u ", "B1g ", "B1u ", "B2g ", "B2u ", "Eg  ", "Eu  " /)
 nallop(196:202) = (/ 3, 10, 16, 8575, 9, 1, 4 /)
 ! C5H
 allrep(199:205) = (/ "C5h ", "A'  ", 'A"  ', "E1' ", 'E1" ', "E2' ", 'E2" ' /)
 nallop(203:208) = (/ 2, 6, 17, 520, 10, 4 /)
 ! D5
 allrep(206:210) = (/ "D5  ", "A1  ", "A2  ", "E1  ", "E2  " /)
 nallop(209:214) = (/ 2, 4, 12, 513, 10, 1 /)
 ! D5H
 allrep(211:219) = (/ "D5h ", "A1' ", 'A1" ', "A2' ", 'A2" ', "E1' ", 'E1" ', "E2' ", 'E2" ' /)
 nallop(215:221) = (/ 3, 8, 18, 537, 10, 4, 5 /)
 ! D5D
 allrep(220:228) = (/ "D5d ", "A1g ", "A1u ", "A2g ", "A2u ", "E1g ", "E1u ", "E2g ", "E2u " /)
 nallop(222:228) = (/ 3, 8, 18, 66130, 10, 7, 5 /)
 ! C6H
 allrep(229:237) = (/ "C6h ", "Ag  ", "Au  ", "Bg  ", "Bu  ", "E1g ", "E1u ", "E2g ", "E2u " /)
 nallop(229:234) = (/ 2, 8, 19, 17612, 11, 7 /)
 ! D6H
 allrep(238:250) = (/ "D6h ", "A1g ", "A1u ", "A2g ", "A2u ", "B1g ", "B1u ", "B2g ", "B2u ", "E1g ", "E1u ", "E2g ", "E2u " /)
 nallop(235:241) = (/ 3, 12, 20, 17663, 11, 1, 7 /)
 ! D6D
 allrep(251:260) = (/ "D6d ", "A1  ", "A2  ", "B1  ", "B2  ", "E1  ", "E2  ", "E3  ", "E4  ", "E5  " /)
 nallop(242:247) = (/ 2, 9, 21, 140468, 18, 5 /)
 ! D7
 allrep(261:266) = (/ "D7  ", "A1  ", "A2  ", "E1  ", "E2  ", "E3  " /)
 nallop(248:253) = (/ 2, 5, 14, 2049, 12, 1 /)
 ! D7H
 allrep(267:277) = (/ "D7h ", "A1' ", 'A1" ', "A2' ", 'A2" ', "E1' ", 'E1" ', "E2' ", 'E2" ', "E3' ", 'E3" ' /)
 nallop(254:260) = (/ 3, 10, 23, 2073, 12, 4, 5 /)
 ! D7D
 allrep(278:288) = (/ "D7d ", "A1g ", "A1u ", "A2g ", "A2u ", "E1g ", "E1u ", "E2g ", "E2u ", "E3g ", "E3u " /)
 nallop(261:267) = (/ 3, 10, 23, 2130, 12, 7, 5 /)
 ! C8H
 allrep(289:299) = (/ "C8h ", "Ag  ", "Au  ", "Bg  ", "Bu  ", "E1g ", "E1u ", "E2g ", "E2u ", "E3g ", "E3u " /)
 nallop(268:273) = (/ 2, 10, 25, 45388, 13, 7 /)
 ! D8H
 allrep(300:314) = (/ "D8h ", "A1g ", "A1u ", "A2g ", "A2u ", "B2u ", "B2g ", "B1u ", "B1g ", "E1g ", "E1u ", "E2g ", "E2u ", &
                      "E3g ", "E3u " /)
 nallop(274:280) = (/ 3, 14, 24, 45439, 13, 4, 5 /)
 ! T
 allrep(315:318) = (/ "T   ", "A   ", "E   ", "T   " /)
 nallop(281:286) = (/ 2, 3, 26, 262276, 8, 3 /)
 ! TD
 allrep(319:324) = (/ "Td  ", "A1  ", "A2  ", "E   ", "T1  ", "T2  " /)
 nallop(287:293) = (/ 3, 5, 28, 270516, 5, 8, 3 /)
 ! O
 allrep(325:330) = (/ "O   ", "A1  ", "A2  ", "E   ", "T2  ", "T1  " /)
 nallop(294:299) = (/ 2, 5, 37, 262273, 1, 8 /)
 ! TH
 allrep(331:337) = (/ "Th  ", "Ag  ", "Au  ", "Eg  ", "Eu  ", "Tg  ", "Tu  " /)
 nallop(300:306) = (/ 3, 6, 27, 278732, 8, 7, 3 /)
 ! OH
 allrep(338:348) = (/ "Oh  ", "A1g ", "A1u ", "A2g ", "A2u ", "Eg  ", "Eu  ", "T1g ", "T1u ", "T2g ", "T2u " /)
 nallop(307:313) = (/ 3, 10, 29, 287231, 1, 7, 8 /)
 ! I
 allrep(349:354) = (/ "I   ", "A   ", "T1  ", "T2  ", "G   ", "H   " /)
 nallop(314:319) = (/ 2, 5, 30, 262657, 1, 10 /)
 ! IH
 allrep(355:365) = (/ "Ih  ", "Ag  ", "Au  ", "T1g ", "T1u ", "T2g ", "T2u ", "Gg  ",  "Gu  ", "Hg  ", "Hu  " /)
 nallop(320:326) = (/ 3, 10, 31, 344786, 10, 7, 8 /)
 ! CooV
 allrep(366:372) = (/ "Coov", "Sg+ ", "Sg- ", "Pi  ", "Del ", "Phi ", "Gam " /)
 nallop(327:333) = (/ 3, 6, 32, 524340, 12, 11, 5 /)
 ! DooH
 allrep(373:385) = (/ "Dooh", "Sgg+", "Sgu+", "Sgg-", "Sgu-", "Pig ", "Piu ", "Delg", "Delu", "Phig", "Phiu", "Gamg", "Gamu" /)
 nallop(334:341) = (/ 4, 12, 33, 524415, 12, 11, 7, 5 /)
 ! R3
 allrep(386:406) = (/ "R3  ", "S(g)", "S(u)", "P(g)", "P(u)", "D(g)", "D(u)", "F(g)", "F(u)", "G(g)", "G(u)", "H(g)", "H(u)", &
                      "I(g)", "I(u)", "K(g)", "K(u)", "L(g)", "L(u)", "M(g)", "M(u)" /)
 nallop(342:348) = (/ 3, 20, 38, 524992, 10, 8, 7 /)

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! new principal axis
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine newaxis(ix,iy,ROT,RNW)
 implicit real(kind=8) (a-h,o-z)
 real(kind=8) :: ROT(3,3),RNW(3,3)

 call Acopy(3,ROT(:,ix),RNW(:,1))
 call Acopy(3,ROT(:,iy),RNW(:,2))
 RNW(1,3) = ROT(2,ix)*ROT(3,iy) - ROT(3,ix)*ROT(2,iy)
 RNW(2,3) = ROT(3,ix)*ROT(1,iy) - ROT(1,ix)*ROT(3,iy)
 RNW(3,3) = ROT(1,ix)*ROT(2,iy) - ROT(2,ix)*ROT(1,iy)

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! rotate the Cartesian coordinate system I1 and I2, and recalculate rotation matrix and atomic coordinates.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine rotmol(Natom,XYZ,sina,cosa,I1,I2,ROT,SCR)
 implicit real(kind=8) (a-h,o-z)
 real(kind=8) :: XYZ(3,Natom),ROT(3,3),SCR(*)

 ! rotate back
 call rotopr(Natom,-1,ROT,XYZ,SCR)

 ! new ROT
 do I = 1, 3
   SCR(1) = -sina*ROT(I,I1) + cosa*ROT(I,I2)
   ROT(I,I1) = cosa*ROT(I,I1) + sina*ROT(I,I2)
   ROT(I,I2) = SCR(1)
 end do

 ! rotate forward
 call rotopr(Natom, 1,ROT,XYZ,SCR)

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! check whether a symmetry operation is true for the whole system
!
! nunmov: NO. of atoms that are not affected by ELEM
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine chi(Natom,tol,ZA,AMS,XYZ,ELEM,IELEM,jelem,nunmov,SC1,SC2)
 implicit real(kind=8) (a-h,o-z)
 Dimension :: jelem(Natom)
 real(kind=8) :: ZA(Natom),AMS(Natom),XYZ(3,Natom),ELEM(3,3),SC1(*),SC2(*)

 IELEM = 0
 nunmov = 0
 do i = 1, Natom
   call acopy(3,XYZ(1,i),SC1)
   call rotopr(1,-1,ELEM,SC1,SC2)
   ifind=0
   do j = 1, Natom
     if (abs(ZA(i)   -ZA(j))  > tol) cycle
     if (abs(AMS(i)  -AMS(j)) > tol) cycle
     if (abs(XYZ(1,j)-SC1(1)) > tol) cycle
     if (abs(XYZ(2,j)-SC1(2)) > tol) cycle
     if (abs(XYZ(3,j)-SC1(3)) > tol) cycle
     ifind=1
     jelem(i)=j
     if (i == j) nunmov = nunmov + 1
     exit
   end do
   if(ifind == 0) exit
 end do
 if(ifind == 1) IELEM = 1

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! define z-axis of a cubic system
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine cubzaxis(Natom,CUBIC,rtol,ROT,XYZ,ipoly,allr,IdxEq,tmp,SCR)
 implicit real(kind=8) (a-h,o-z)
 real(kind=8) :: XYZ(3,Natom),ROT(3,3),allr(3),tmp(3,3),SCR(*)
 dimension :: ipoly(3),IdxEq(*)
 logical :: CUBIC,OK

 Pi = acos(-1.d0)

 ! for an arbitary non-central atom, calculate the distance
 do i = 1, Natom
   SCR(1) = sqrt( dotx(3,XYZ(1,i),XYZ(1,i)) )
   if(SCR(1) > rtol) exit
 end do
 ! count the number of equivalent atoms
 Neq = 0
 do i = 1, Natom
   SCR(2) = sqrt( dotx(3,XYZ(1,i),XYZ(1,i)) )
   if(SCR(2) < rtol) cycle
   if (SCR(2) > SCR(1) + rtol) cycle
   if(abs(SCR(1)-SCR(2)) < rtol)then
     Neq = Neq + 1
     IdxEq(Neq) = i
   else
     Neq = 0
   end if
   idx0 = i
   SCR(1) = SCR(2)
 end do

 if(Neq == 0) then
   call XError(.false.,"Unreasonable cubic geometry (I).")
 else if(Neq == 1) then
   CUBIC = .False.
   return
 end if

 OK = (Neq == 4 ) .or. (Neq == 6 ) .or. (Neq == 8 ) .or. (Neq == 12) .or. (Neq == 20)

 if(OK) then
   ! Atoms nearest to the center of symmetry outline one of the Platonic
   ! solids (tetrahedron, octahedron, cube, icosahedron, or pentagonal
   ! dodecahedron). Each atom lies on a 3, 4, or 5 fold symmetry axis.

   ! calculate the distance between equivalent atoms
   SCR(1)=1.0d4
   j1 = IdxEq(1)
   do i = 2, Neq
     j2 = IdxEq(i)
     SCR(1) = min(SCR(1), distance(XYZ(1,j1),XYZ(1,j2)) )
   end do

   ! count the number of equidistant nearest neighbors
   Neb = 0
   j1 = IdxEq(1)
   do i = 2, Neq
     j2 = IdxEq(i)
     SCR(2) = distance(XYZ(1,j1),XYZ(1,j2))
     if(abs(SCR(1)-SCR(2)) < rtol) Neb = Neb + 1
   end do
   OK = (Neq ==  4 .and. Neb == 3) .or. (Neq ==  6 .and. Neb == 4) .or. (Neq ==  8 .and. Neb == 3) .or. &
        (Neq == 12 .and. Neb == 5) .or. (Neq == 20 .and. Neb == 3)
   if(OK) then
     ! z axis is defined by the first atom in the equivalent atom set
     call acopy(3,XYZ(1,idx0),ROT(1,3))
     go to 1000
   end if
 end if

 i1=IdxEq(1)
 do i=1,3
   SCR(1)=1.0d4
   ipoly(i)=0
   Loop1: do j1 = 2, Neq
     j = IdxEq(j1)
     do m = 1, i - 1
       if (ipoly(m) /= j) cycle
       cycle  Loop1
     end do
     SCR(2)=distance(XYZ(1,i1),XYZ(1,j))
     if (SCR(1) <= SCR(2)) cycle  Loop1
     SCR(1) = SCR(2)
     ipoly(i) = j
     allr(i) = SCR(2)
   end do Loop1
 end do
 !if(ipoly(1)*ipoly(2)*ipoly(3) == 0) then
 !  CUBIC = .False.
 !  return
 !end if

 ! identify the two atoms in the regular polygon.
 Loop2: do i=1,2
   do j=i+1,3
     if(abs(allr(i)-allr(j)) >= rtol) cycle
     exit Loop2
   end do
 end do Loop2
 if(ipoly(i)*ipoly(j) == 0) then
   CUBIC = .False.
   return
 end if
 do k=1,3
   tmp(k,3)=XYZ(k,ipoly(i))
   tmp(k,2)=XYZ(k,ipoly(j))
   tmp(k,1)=XYZ(k,i1)
 end do
 i3=ipoly(i)
 i2=ipoly(j)
 ipoly(1)=i3
 ipoly(2)=i2

 angle = BAngle(tmp(1,3),tmp(1,1),tmp(1,2))

 if(abs(angle-Pi/3.d0) < rtol)then
   ! triangle
   do i=1,3
     ROT(i,3)=XYZ(i,i1)+XYZ(i,i2)+XYZ(i,i3)
   end do
 else if(abs(angle-Pi*0.5d0) < rtol)then
   ! square
   do i=1,3
     ROT(i,3)=XYZ(i,i2)+XYZ(i,i3)
   end do
 else if(abs(angle-Pi*0.6d0) < rtol)then
   ! pentagon
   do i=1,3
     ROT(i,3)=1.6180341d0*(XYZ(i,i2)+XYZ(i,i3))-XYZ(i,i1)
   end do
 else
   CUBIC = .False.
   return
 endif

 1000  continue
 ! x axis, which makes that two atoms have equal y and z coordinates.
 if(abs(ROT(2,3)) <= abs(ROT(3,3)))then
   ROT(1,1) = ROT(3,3)
   ROT(2,1) = 0.d0
   ROT(3,1) = -ROT(1,3)
 else
   ROT(1,1) = ROT(2,3)
   ROT(2,1) = -ROT(1,3)
   ROT(3,1) = 0.d0
 end if

 ! y axis: y = z X x
 call CrossX(ROT(1,3),ROT(1,1),ROT(1,2))

 call renorm(3,ROT(1,1))
 call renorm(3,ROT(1,2))
 call renorm(3,ROT(1,3))

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Initialization: generate symmetry operations used in this program
!
! 1 C2(X)       6 Sigma(YZ)   11 C6    16 S8
! 2 C2(Y)       7 inversion   12 C7    17 S10
! 3 C2(Z)       8 C3          13 C8    18 S12
! 4 Sigma(XY)   9 C4          14 S4    19 0/1 if cubic
! 5 Sigma(XZ)  10 C5          15 S6    20 0/1 if infinite
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine symini(Natom,NOper,PGNAME,IRNAME,IELEM,ELEM,JELEM,CUB)
 implicit real(kind=8) (a-h,o-z)
 parameter(pi=3.14159265358979d0, dp=pi*2.d0, half=0.5d0, s036=sin(pi/5.d0), c036=cos(pi/5.d0), s045=sin(pi/4.d0), &
   s051=sin(dp/7.d0), c051=cos(dp/7.d0), s060=sin(pi/3.d0), s072=sin(dp/5.d0), c072=cos(dp/5.d0))
 Dimension :: IELEM(NOper),JELEM(Natom,NOper)
 real(kind=8) :: ELEM(3,3,NOper),CUB(3,3)
 character*4 :: PGNAME,IRNAME(*)

 PGNAME="****"
 do i = 1, Natom*3
   IRNAME(i) = PGNAME
 end do
 IELEM = 0
 JELEM = 0
 call UMAT(3,CUB)

 ELEM(:,1, 1)=(/  1.d0,  0.d0,  0.d0/)
 ELEM(:,2, 1)=(/  0.d0, -1.d0,  0.d0/)
 ELEM(:,3, 1)=(/  0.d0,  0.d0, -1.d0/)

 ELEM(:,1, 2)=(/ -1.d0,  0.d0,  0.d0/)
 ELEM(:,2, 2)=(/  0.d0,  1.d0,  0.d0/)
 ELEM(:,3, 2)=(/  0.d0,  0.d0, -1.d0/)

 ELEM(:,1, 3)=(/ -1.d0,  0.d0,  0.d0/)
 ELEM(:,2, 3)=(/  0.d0, -1.d0,  0.d0/)
 ELEM(:,3, 3)=(/  0.d0,  0.d0,  1.d0/)

 ELEM(:,1, 4)=(/  1.d0,  0.d0,  0.d0/)
 ELEM(:,2, 4)=(/  0.d0,  1.d0,  0.d0/)
 ELEM(:,3, 4)=(/  0.d0,  0.d0, -1.d0/)

 ELEM(:,1, 5)=(/  1.d0,  0.d0,  0.d0/)
 ELEM(:,2, 5)=(/  0.d0, -1.d0,  0.d0/)
 ELEM(:,3, 5)=(/  0.d0,  0.d0,  1.d0/)

 ELEM(:,1, 6)=(/ -1.d0,  0.d0,  0.d0/)
 ELEM(:,2, 6)=(/  0.d0,  1.d0,  0.d0/)
 ELEM(:,3, 6)=(/  0.d0,  0.d0,  1.d0/)

 ELEM(:,1, 7)=(/ -1.d0,  0.d0,  0.d0/)
 ELEM(:,2, 7)=(/  0.d0, -1.d0,  0.d0/)
 ELEM(:,3, 7)=(/  0.d0,  0.d0, -1.d0/)

 ELEM(:,1, 8)=(/ -half,  s060,  0.d0/)
 ELEM(:,2, 8)=(/ -s060, -half,  0.d0/)
 ELEM(:,3, 8)=(/  0.d0,  0.d0,  1.d0/)

 ELEM(:,1, 9)=(/  0.d0,  1.d0,  0.d0/)
 ELEM(:,2, 9)=(/ -1.d0,  0.d0,  0.d0/)
 ELEM(:,3, 9)=(/  0.d0,  0.d0,  1.d0/)

 ELEM(:,1,10)=(/  c072,  s072,  0.d0/)
 ELEM(:,2,10)=(/ -s072,  c072,  0.d0/)
 ELEM(:,3,10)=(/  0.d0,  0.d0,  1.d0/)

 ELEM(:,1,11)=(/  half,  s060,  0.d0/)
 ELEM(:,2,11)=(/ -s060,  half,  0.d0/)
 ELEM(:,3,11)=(/  0.d0,  0.d0,  1.d0/)

 ELEM(:,1,12)=(/  c051,  s051,  0.d0/)
 ELEM(:,2,12)=(/ -s051,  c051,  0.d0/)
 ELEM(:,3,12)=(/  0.d0,  0.d0,  1.d0/)

 ELEM(:,1,13)=(/  s045,  s045,  0.d0/)
 ELEM(:,2,13)=(/ -s045,  s045,  0.d0/)
 ELEM(:,3,13)=(/  0.d0,  0.d0,  1.d0/)

 ELEM(:,1,14)=(/  0.d0,  1.d0,  0.d0/)
 ELEM(:,2,14)=(/ -1.d0,  0.d0,  0.d0/)
 ELEM(:,3,14)=(/  0.d0,  0.d0, -1.d0/)

 ELEM(:,1,15)=(/  half,  s060,  0.d0/)
 ELEM(:,2,15)=(/ -s060,  half,  0.d0/)
 ELEM(:,3,15)=(/  0.d0,  0.d0, -1.d0/)

 ELEM(:,1,16)=(/  s045,  s045,  0.d0/)
 ELEM(:,2,16)=(/ -s045,  s045,  0.d0/)
 ELEM(:,3,16)=(/  0.d0,  0.d0, -1.d0/)

 ELEM(:,1,17)=(/  c036,  s036,  0.d0/)
 ELEM(:,2,17)=(/ -s036,  c036,  0.d0/)
 ELEM(:,3,17)=(/  0.d0,  0.d0, -1.d0/)

 ELEM(:,1,18)=(/  s060,  half,  0.d0/)
 ELEM(:,2,18)=(/ -half,  s060,  0.d0/)
 ELEM(:,3,18)=(/  0.d0,  0.d0, -1.d0/)

 call AClear(18,ELEM(1,1,19))

 return
end

!--- END
