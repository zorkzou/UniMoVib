# UniMoVib
A unified interface for molecular harmonic vibrational frequency calculations.

The UniMoVib program was originally written by Wenli Zou in FORTRAN 77 during 2014 and 2015 at Southern Methodist University (SMU), Dallas, Texas, within the framework of the LocalMode (now LModeA) program of the Computational and Theoretical Chemistry Group (CATCO) of SMU. This work was supported by the NSF grants CHE 1152357 and CHE 1464906. Guidance from the late Dr. Dieter Cremer is acknowledged. After being rewritten in Fortran 90 in the spring of 2017, UniMoVib has been released as a stand-alone program.

## Latest Versions
Version 1.4.4 (Dec/28/2021).

1. Improved GCC version 10 (gfortran) compatibility.

Version 1.4.3 (Jul/19/2021).

1. The most abundant or the longest lived isotopic masses of the elements Rb-Np, Hs, Rg, and Og have been updated.

Version 1.4.2 (Jul/03/2021).

1. [xTB](https://github.com/grimme-lab/xtb/) has been supported.

Version 1.4.1 (May/02/2021).

1. [PyVibMS](https://github.com/smutao/PyVibMS) has been supported to visualize vibrational modes (by Y. Tao).
2. The latest version of ifort with MKL may be used in Makefile.

Version 1.4.0 (Mar/22/2021).

1. Raman scattering activities and depolarization ratios may be calculated for the data files saved by Gaussian, Gamess, Firefly, and Orca.
2. The Hessian file saved by the latest version of Orca may be read correctly.

Version 1.3.5 (Nov/20/2020).

1. An ASCII data file of subsystem may be generated by GSVA.

Version 1.3.4 (Oct/24/2020).

1. The format of UniMoVib data file has been updated. See sec. A.1 of the manual.

Version 1.3.3 (Jun/22/2020).

1. Bug fix in symmetry analysis.
2. The manual has been updated.

Version 1.3.2 (May/28/2020).

1. Bug fix in symmetry analysis.

Version 1.3.1 (May/22/2020).

1. Generalized subsystem vibrational analysis (GSVA) by Y. Tao may be performed.
2. Cartesian coordinates may be provided in the input file through `qcprog=xyzinp`.

Version 1.3.0 (Apr/29/2020).

1. Due to Jahn-Teller effects or numerical noise, sometimes the irreps of vibrational normal modes cannot be determined by the program. A new keyword `IFSymtz` has been introduced into the program which may symmetrize the vibrational normal modes.
2. The longest lived isotopic masses have been updated for the elements with Z > 93.

## Features

1. Calculate harmonic vibrational frequencies and (optional) I.R. & Raman intensities from Hessian, coordinates, and other related data generated by quantum chemistry programs or by the user manually. Nearly 30 quantum chemistry programs have been supported.
2. Calculate atomic IR charges of planar and linear molecules. Reference: Theor. Chem. Acc. 131, 1139 (2012).
3. Analyze point group of geometry and irreducible representations of normal modes in full symmetry.
4. Thermochemistry calculation uses the point group in full symmetry, and the results are printed in Gaussian-style.
5. Save data files for animation of normal modes using [Gabedit](http://gabedit.sourceforge.net/), [Molden](https://www.theochem.ru.nl/molden/), or [PyVibMS](https://github.com/smutao/PyVibMS).
6. Set up isotopic masses, temperature, pressure, scale factor and/or experimental frequencies, and so on.
7. Can be used as a third party module for frequency and thermochemistry calculations in a quantum chemistry program, for example, [BDF](http://182.92.69.169:7226/).
8. Interface to [LModeA](https://sites.smu.edu/dedman/catco/) for the local mode analysis (e.g. force constants of chemical bonds, bond angles, and so on).
9. Generalized subsystem vibrational analysis (GSVA). References: J. Chem. Theory Comput. 14, 2558 (2018) and Theor. Chem. Acc. 140,
31 (2021).

## Supported quantum chemistry programs

* [Aces](http://www.qtp.ufl.edu/ACES/)
* [Adf](http://www.scm.com/)
* Ampac 2.x. See [Semichem, Inc.](http://www.semichem.com/)
* [Amsol](http://comp.chem.umn.edu/amsol/)
* [BDF](http://182.92.69.169:7226/)
* [CFour](http://www.cfour.de/)
* [Columbus](http://www.univie.ac.at/columbus/)
* [CP2k](http://www.cp2k.org/)
* [Crystal](http://www.crystal.unito.it/)
* [Dalton](http://daltonprogram.org/)
* [deMon2k](http://www.demon-software.com/public_html/)
* [Dmol3](http://accelrys.com/)
* [Fhi-Aims](https://aimsclub.fhi-berlin.mpg.de/)
* [Firefly](http://classic.chem.msu.su/gran/gamess/)
* [Gabedit](http://gabedit.sourceforge.net/)
* [Gamess](http://www.msg.chem.iastate.edu/gamess/)
* [Gamess-UK](http://www.cfs.dl.ac.uk/)
* [Gaussian](http://www.gaussian.com/)
* [Hyperchem](http://www.hyper.com/)
* [Jaguar](http://www.schrodinger.com/)
* [Molcas](http://www.molcas.org/) and [OpenMolcas](https://gitlab.com/Molcas/OpenMolcas)
* [Molden](https://www.theochem.ru.nl/molden/)
* [Molpro](http://www.molpro.net/)
* [Mopac](http://openmopac.net/). See also Mo-g in [Scigress](http://www.scigress.com/)
* [NWChem](http://www.nwchem-sw.org/index.php/Main_Page)
* [Orca](https://orcaforum.kofo.mpg.de)
* [Pqs](http://www.pqs-chem.com/)
* [Psi4](http://www.psicode.org/)
* [Q-Chem](http://www.q-chem.com/)
* [Spartan](http://www.wavefun.com/)
* [Turbomole](http://www.cosmologic.de/)
* [xTB](https://github.com/grimme-lab/xtb/)
