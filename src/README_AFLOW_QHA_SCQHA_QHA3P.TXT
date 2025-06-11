********************************************************************************
*                                                                              *
*             Aflow STEFANO CURTAROLO - Duke University 2003-2021              *
*                High-Throughput ab-initio Materials Discovery                 *
*                                                                              *
********************************************************************************
LATEST VERSION OF THE FILE:        materials.duke.edu/AFLOW/README_AFLOW_QHA_SCQHA_QHA3P.TXT
********************************************************************************

Quasi-harmonic approximation (QHA), QHA-3phonon (QHA3P), and self consistent quasi harmonic (SCQHA)  
(Gruneisen parameter and Equation of state calculations)
Thermal properties calculated using a QHA, QHA3P, and SCQHA model.
This implementation in AFLOW written by Pinku Nath and Stefano Curtarolo,
and is described in the paper "Computational Materials Science Volume 125, December 2016, Pages 82–91" and arXiv:1807.04669.


    Within QHA following properties can be calculated..
    Gruneisen parameter along high-symmetry path
    Mode Gruneisen in q-mesh
    Average Gruneisen parameter
    Equation of states (Cp, thermal-expansion, Bulk Modulus etc.)

    The above mentioned properties can also be calculated with QHA3P and SCQHA.
    In addition to that SCQHA and QHA3P can also calculate Temperature dependent phonon spectra
    which can not be calculated with QHA method. 

    In this version QHA, QHA3P and SCQHA can be run alone or together.

    This read-me file is divided into following sections:
    SECTION1 ::  All AFLOW options are explained to run QHA, QHA3P and SCQHA
    SECTION2 ::  All output file names
    SECTION3 ::  One example of Si aflow.in to run QHA, QHA3P and SCQHA together
    SECTION4 ::  One example of qsub.sh


********************************* SECTION1 START *****************************************

All the three QHA methods are based on phonon calculations, therefore it is necessary 
to use AFLOW-APL variables while using any of the QHA methods.
Necessary AFLOW-APL variables are shown below:

 
*******************************[APL]******************************************************
[AFLOW_APL]CALC
[AFLOW_APL]ENGINE=DM
[AFLOW_APL]DISMAG=0.015
[AFLOW_APL]SUPERCELL=3x3x3
[AFLOW_APL]ZEROSTATE=y
[AFLOW_APL]POLAR=n 
[AFLOW_APL]DC=y                                //necessary for QHA, QHA3P and SCQHA
[AFLOW_APL]DCINITSG=227                        
[AFLOW_APL]DCUSERPATH=Gamma-X|X-W|W-L|L-Gamma
[AFLOW_APL]DOS=y                               //necessary for QHA, QHA3P and SCQHA
[AFLOW_APL]TP=y                                //Thermal properties for APL as well as QHA, QHA3P and SCQHA
[AFLOW_APL]TPT=0:2000:10                       //temperature range for APL as well as QHA, QHA3P and SCQHA  


NB: In all QHA methods temperature range are given by [AFLOW_APL]TPT option
*******************************[QHA]*****************************************************
1) To perform a AFLOW-QHA isotropic properties calculation the following options can be used:

[AFLOW_QHA]CALC           // QHA option. Default value is OFF         
[AFLOW_QHA]MODE=QHA       // Calculate thermal  properties using QHA,  Default value is QHA
[AFLOW_QHA]EOS=y          // Calculate QHA equation of state.  Default value is OFF

N.B:
 (1)QHA is going to be performed with following default Gruneisen and EOS (equation of state) distortions. 
    (a) GP_DISTORTION=0.03            // [percentage of volume distortion magnitude w.r.t equilibrium volume, default value is 0.03]
    (b) EOS_DISTORTION_RANGE=-3:6:1   // [distortion range to perform phonon and
                                      //  static energies calculations within the range between 3% and 6% with 1% increment.
                                      //  The default value is -3:6:1].

 Follow options in (9) and (12) to change these these default values.

 (2)QHA can also be run with QHA3P and SCQHA options. Check example in SECTION3

*******************************[QHA anisotropy along lattice vector a ]******************

2) To perform a AFLOW-QHA anisotropic  properties along lattice vector a the following options can be used:

[AFLOW_QHA]CALC           // QHA option. Default value is OFF         
[AFLOW_QHA]MODE=QHA_A     // Calculate thermal properties using QHA along lattice vector a. Default value is QHA
[AFLOW_QHA]EOS=y          // Calculate QHA equation of state.  Default value is OFF

N.B:
 (1)QHA is going to be performed with following default Gruneisen and EOS distortions. 
    (a) GP_DISTORTION=0.03            // [percentage of volume distortion magnitude w.r.t equilibrium volume, default value is 0.03]
    (b) EOS_DISTORTION_RANGE=-3:6:1   // [EOS distortion range to perform phonon and static calculations]
 Follow options in (9) and (12) to change these these default values.

 (2) QHA_A is mutually exclusive to QHA option

*******************************[QHA anisotropy along lattice vector b ]******************
3) To perform a AFLOW-QHA anisotropic  properties along lattice vector b the following options can be used:

[AFLOW_QHA]CALC           // QHA option. Default value is OFF         
[AFLOW_QHA]MODE=QHA_B     // Calculate thermal properties using QHA along lattice vector b
[AFLOW_QHA]EOS=y          // Calculate QHA equation of state

N.B:
 (1)QHA_B is going to be performed with following default Gruneisen and EOS distortions. 
    (a) GP_DISTORTION=0.03            // [percentage of volume distortion magnitude w.r.t equilibrium volume, default value is 0.03]
    (b) EOS_DISTORTION_RANGE=-3:6:1   // [EOS distortion range to perform phonon and static calculations]
 Follow options in (9) and (12) to change these these default values.

 (2)QHA_B is mutually exclusive to both QHA and QHA_A options

*******************************[QHA anisotropy along lattice vector c ]******************
4) To perform a AFLOW-QHA anisotropic  properties along lattice vector c the following options can be used:

[AFLOW_QHA]CALC           // QHA option. Default value is OFF         
[AFLOW_QHA]MODE=QHA_C     // Calculate thermal properties using QHA along lattice vector c
[AFLOW_QHA]EOS=y          // Calculate QHA equation of state

N.B:
 (1)QHA_C is going to be performed with following default Gruneisen and EOS distortions. 
    (a) GP_DISTORTION=0.03            // [percentage of volume distortion magnitude w.r.t equilibrium volume, default value is 0.03]
    (b) EOS_DISTORTION_RANGE=-3:6:1   // [EOS distortion range to perform phonon and static calculations]
 Follow options in (9) and (12) to change these these default values.

 (2) QHA_C is mutually exclusive to QHA, QHA_A and QHA_C options

*******************************[QHA3P]**************************************************
5) To perform a AFLOW-QHA3P isotropic properties calculation the following options can be used:
[AFLOW_QHA]CALC           // QHA option. Default value is OFF
[AFLOW_QHA]MODE=QHA3P     // Calculate thermal properties using QHA3P
[AFLOW_QHA]EOS=y          // Calculate QHA3P equation of state

N.B:
 (1)QHA3P is going to be performed with following default Gruneisen and EOS distortions. 
    (a) SCQHA_DISTORTION=3            // [percentage of QHA3P or SCQHA volume distortion magnitude w.r.t equilibrium volume, default value is 3]
    (b) EOS_DISTORTION_RANGE=-3:6:1   // [EOS distortion range to perform phonon and static calculations]
 Follow options in (10) and (12) to change these these default values.

 (2)QHA3P can also be run with SCQHA option. Check example in SECTION3

*******************************[QHA3P anisotropy along lattice vector a ]******************
6) To perform a AFLOW-QHA3P anisotropic  properties along lattice vector a the following options can be used:

[AFLOW_QHA]CALC           // QHA option. Default value is OFF
[AFLOW_QHA]MODE=QHA3P_A   // Calculate thermal properties using QHA3P along lattice vector a
[AFLOW_QHA]EOS=y          // Calculate QHA3P equation of state

N.B:
 (1)QHA3P is going to be performed with following default Gruneisen and EOS distortions. 
    (a) SCQHA_DISTORTION=3            // [percentage of QHA3P or SCQHA volume distortion magnitude w.r.t equilibrium volume, default value is 3]
    (b) EOS_DISTORTION_RANGE=-3:6:1   // [EOS distortion range to perform phonon and static calculations]
 Follow options in (10) and (12) to change these these default values.

 (2)QHA3P can also be run with SCQHA option. Check example in SECTION3

 (3)QHA3P_A is mutually exclusive to QHA3P option

*******************************[QHA3P anisotropy along lattice vector b ]******************
7) To perform a AFLOW-QHA3P anisotropic  properties along lattice vector a the following options can be used:

[AFLOW_QHA]CALC           // QHA option. Default value is OFF
[AFLOW_QHA]MODE=QHA3P_B   // Calculate thermal properties using QHA3P along lattice vector b
[AFLOW_QHA]EOS=y          // Calculate  QHA3P equation of state

N.B:
 (1)QHA3P_B is going to be performed with following default Gruneisen and EOS distortions. 
    (a) SCQHA_DISTORTION=3            // [percentage of QHA3P or SCQHA volume distortion magnitude w.r.t equilibrium volume, default value is 3]
    (b) EOS_DISTORTION_RANGE=-3:6:1   // [EOS distortion range to perform phonon and static calculations]
 Follow options in (10) and (12) to change these these default values.

 (2)QHA3P_B is mutually exclusive to QHA3P and QHA3P_B option

*******************************[SCQHA anisotropy along lattice vector c ]******************
8) To perform a AFLOW-QHA3P anisotropic  properties along lattice vector c the following options can be used:

[AFLOW_QHA]CALC           // QHA option. Default value is OFF
[AFLOW_QHA]MODE=QHA3P_C   // Calculate thermal properties using QHA3P along lattice vector c
[AFLOW_QHA]EOS=y          // Calculate QHA3P equation of state

N.B:
 (1)QHA3P_B is going to be performed with following default Gruneisen and EOS distortions. 
    (a) SCQHA_DISTORTION=3            // [percentage of QHA3P or SCQHA volume distortion magnitude w.r.t equilibrium volume, default value is 3]
    (b) EOS_DISTORTION_RANGE=-3:6:1   // [EOS distortion range to perform phonon and static calculations]
 Follow options in (10) and (12) to change these these default values.

 (2)QHA3P_C is mutually exclusive to QHA3P, QHA3P_B and QHA3P_C options

*******************************[SCQHA]**************************************************
5) To perform a AFLOW-SCQHA isotropic properties calculation the following options can be used:
[AFLOW_QHA]CALC           // QHA option. Default value is OFF
[AFLOW_QHA]MODE=SCQHA     // Calculate thermal properties using SCQHA
[AFLOW_QHA]EOS=y          // Calculate SCQHA equation of state

N.B:
 (1)SCQHA is going to be performed with following default Gruneisen and EOS distortions. 
    (a) SCQHA_DISTORTION=3            // [percentage of QHA3P or SCQHA volume distortion magnitude w.r.t equilibrium volume, default value is 3]
    (b) EOS_DISTORTION_RANGE=-3:6:1   // [EOS distortion range to perform phonon and static calculations]
 Follow options in (10) and (12) to change these these default values.

 (2)SCQHA can also be run with QHA3P option. Check example in SECTION3

*******************************[SCQHA anisotropy along lattice vector a ]******************
6) To perform a AFLOW-SCQHA anisotropic  properties along lattice vector a the following options can be used:

[AFLOW_QHA]CALC           // QHA option. Default value is OFF
[AFLOW_QHA]MODE=SCQHA_A   // Calculate thermal properties using SCQHA along lattice vector a
[AFLOW_QHA]EOS=y          // Calculate SCQHA equation of state

N.B:
 (1)SCQHA_A is going to be performed with following default Gruneisen and EOS distortions. 
    (a) SCQHA_DISTORTION=3            // [percentage of QHA3P or SCQHA volume distortion magnitude w.r.t equilibrium volume, default value is 3]
    (b) EOS_DISTORTION_RANGE=-3:6:1   // [EOS distortion range to perform phonon and static calculations]
 Follow options in (10) and (12) to change these these default values.

 (2)SCQHA_A is mutually exclusive to SCQHA option

*******************************[SCQHA anisotropy along lattice vector b ]******************
7) To perform a AFLOW-SCQHA anisotropic  properties along lattice vector a the following options can be used:

[AFLOW_QHA]CALC           // QHA option. Default value is OFF
[AFLOW_QHA]MODE=SCQHA_B   // Calculate thermal properties using SCQHA along lattice vector b
[AFLOW_QHA]EOS=y          // Calculate SCQHA equation of state

N.B:
 (1)SCQHA is going to be performed with following default Gruneisen and EOS distortions. 
    (a) SCQHA_DISTORTION=3            // [percentage of QHA3P or SCQHA volume distortion magnitude w.r.t equilibrium volume, default value is 3]
    (b) EOS_DISTORTION_RANGE=-3:6:1   // [EOS distortion range to perform phonon and static calculations]
 Follow options in (10) and (12) to change these these default values.

 (2)SCQHA can also be run with QHA3P option. Check example in SECTION3

 (3)SCQHA_A is mutually exclusive to SCQHA and SCQHA_A options

*******************************[SCQHA anisotropy along lattice vector c ]******************
8) To perform a AFLOW-SCQHA anisotropic  properties along lattice vector c the following options can be used:

[AFLOW_QHA]CALC           // QHA option. Default value is OFF
[AFLOW_QHA]MODE=SCQHA_C   // Calculate thermal  properties using SCQHA along lattice vector c. Default value is QHA
[AFLOW_QHA]EOS=y          // Calculate SCQHA equation of state. Default value is OFF

N.B:
 (1)SCQHA is going to be performed with following default Gruneisen and EOS distortions. 
    (a) SCQHA_DISTORTION=3            // [percentage of QHA3P or SCQHA volume distortion magnitude w.r.t equilibrium volume, default value is 3]
    (b) EOS_DISTORTION_RANGE=-3:6:1   // [EOS distortion range to perform phonon and static calculations]
 Follow options in (10) and (12) to change these these default values.

 (2)SCQHA can also be run with QHA3P option. Check example in SECTION3

 (3)SCQHA_C is mutually exclusive to SCQHA, SCQHA_A and SCQHA_C options


*****************************************************************************************
9)Some other common options to QHA.
[AFLOW_QHA]GP_DISTORTION=0.03            // [percentage of volume distortion magnitude w.r.t equilibrium volume, default value is 0.03]

*****************************************************************************************
10)Some other common options to QHA3P and SCQHA.
[AFLOW_QHA]SCQHA_DISTORTION=3            // [percentage of QHA3P or SCQHA volume distortion magnitude w.r.t equilibrium volume, default value is 3]


*****************************************************************************************
11)Some other SCQHA options
[AFLOW_QHA]SCQHA_PDIS_T=50,100,600    // Temperature dependent phonon dispersion spectra calculation at 50, 100 and 600K.
                                      // Default value is OFF.
*****************************************************************************************
12)Some other common options to both QHA, QHA3P and SCQHA.

[AFLOW_QHA]EOS_DISTORTION_RANGE=-3:6:1   // distortion range to perform phonon and
                                         // static energies calculations within the range between 3% and 6% with 1% increment.
                                         // The default value is -3:6:1.

[AFLOW_QHA]EOS_STATIC_KPPRA=8000, This is the number of k-points for static energy calculations 
The default value is 8000.

[AFLOW_QHA]NEDOS=2000, This is the number of electronic density of states for static energy calculations.
The default value is 5000.

[AFLOW_QHA]INCLUDE_ELE=y, This is option is for including temperature dependent electronic
effect. Default value is OFF.

********************************* SECTION2 START ***************************************
Describe all output files 

 aflow.qha.distortions.out => It contains QHA, SCQHA, and QHA3P  distortion information

 aflow.apl_hskpoints.out   => It contains highsymmetry q-point information

 aflow.qha.gpdis.out       =>  It contains QHA Gruneisen dispersion curves along high-symmetry q-points

 aflow.qha.gp.mesh.out     => It contains QHA Gruneisen dispersion curves in  q-points mesh

 aflow.qha.avg_gp.out      => It contains QHA average Gruneisen parameters

 aflow.qha.static_energies.out => It contains static energies

 aflow.qha.FVT.out         => It contains QHA free energy, volume and temperature information

 aflow.qha.err_fit.out     => It contains QHA fitting related information

 aflow.qha.thermo.out      => It contains QHA EOS

 aflow.qha.enthalpy.out    => It contains QHA enthalpies

 aflow.mesh.taylor_cofficients.out => It contains Taylor coefficients in q-mesh

 aflow.qha3P.gp.mesh.out  => It contains QHA3P Gruneisen parameter in q-mesh

 aflow.qha3P.avg_gp.out   => It contains QHA3P average Gruneisen parameters

 aflow.scf.thermo.out     => It contains SCQHA EOS
 
 aflow.scf.iter.out       => It contains SCQHA self-consistent volumes

 aflow.scqha.enthalpy.out => It contains SCQHA enthalpies

 aflow.scf.pressure.out   => It contains SCQHA pressures 

 aflow.qha3P.thermo.out   => It contains QHA3P EOS

 aflow.qha3P.FVT.out      => It contains QHA3P free energy, volume and temperature information

 aflow.qha3P.enthalpy.out => It contains QHA3P enthalpies

 aflow.scqha_pdis_T     => It contains SCQHA PDIS (phonon dispersion) at T K

********************************* SECTION2 END *****************************************


********************************* SECTION3 START ***************************************

[AFLOW] **************************************************************************************************************************
[AFLOW]
[AFLOW]                     .o.        .o88o. oooo
[AFLOW]                    .888.       888 `` `888
[AFLOW]                   .8'888.     o888oo   888   .ooooo.  oooo oooo    ooo
[AFLOW]                  .8' `888.     888     888  d88' `88b  `88. `88.  .8'
[AFLOW]                 .88ooo8888.    888     888  888   888   `88..]88..8'
[AFLOW]                .8'     `888.   888     888  888   888    `888'`888'
[AFLOW]               o88o     o8888o o888o   o888o `Y8bod8P'     `8'  `8'  .in
[AFLOW]
[AFLOW] **************************************************************************************************************************
[AFLOW] * Stefano Curtarolo - (aflow V30735)
[AFLOW] * D. Morgan, W. Setyawan, G. Hart, M. Jahnatek, S. Wang, O. Levy, K. Yang, J. Xue,
[AFLOW] * R. Taylor, C. Calderon, C. Toher, R. Chepulskyy, K. Rasch, M. Buongiorno Nardelli
[AFLOW] **************************************************************************************************************************
[AFLOW] Aflow automatically generated (aflow_avasp.cpp)
[AFLOW] **************************************************************************************************************************
[AFLOW]SYSTEM=Si2
#[AFLOW] single element calculation
[AFLOW] **************************************************************************************************************************
[AFLOW] input file for aflow
[AFLOW_MODE=VASP]
[AFLOW] **************************************************************************************************************************
[AFLOW_MODE_ZIP=bzip2]
[AFLOW_MODE_BINARY=vasp53s]
[AFLOW] **************************************************************************************************************************
[AFLOW] **************************************************************************************************************************
#[AFLOW_MODE_MPI]
[AFLOW_MODE_MPI_MODE]NCPUS=MAX
[AFLOW_MODE_MPI_MODE]COMMAND ="mpirun -np"
[AFLOW_MODE_MPI_MODE]AUTOTUNE
[AFLOW_MODE_MPI_MODE]BINARY="mpivasp53s"
[AFLOW] **************************************************************************************************************************
[AFLOW_SYMMETRY]CALC
#[AFLOW_SYMMETRY]SGROUP_WRITE
#[AFLOW_SYMMETRY]SGROUP_RADIUS=7.77
[AFLOW] **************************************************************************************************************************
#[AFLOW_NEIGHBOURS]CALC
[AFLOW_NEIGHBOURS]RADIUS=7.7
[AFLOW_NEIGHBOURS]DRADIUS=0.1
[AFLOW] **************************************************************************************************************************
[AFLOW_APL]CALC
[AFLOW_APL]ENGINE=DM
[AFLOW_APL]DISMAG=0.015
[AFLOW_APL]SUPERCELL=3x3x3
[AFLOW_APL]ZEROSTATE=y
[AFLOW_APL]POLAR=n
[AFLOW_APL]DC=y
[AFLOW_APL]DCINITSG=227
[AFLOW_APL]DCUSERPATH=Gamma-X|X-W|W-L|L-Gamma
[AFLOW_APL]DOS=y
[AFLOW_APL]TP=y
[AFLOW_APL]TPT=0:300:10
[AFLOW] **************************************************************************************************************************
[AFLOW_QHA]CALC
[AFLOW_QHA]MODE=QHA,QHA3P,SCQHA
[AFLOW_QHA]EOS=y
[AFLOW] **************************************************************************************************************************
[VASP_FORCE_OPTION]CHGCAR=OFF                         // ON | OFF (default ON)
[VASP_FORCE_OPTION]SYM=ON                             // ON | OFF  (default ON)
[VASP_FORCE_OPTION]AUTO_PSEUDOPOTENTIALS=potpaw_PBE   // pot_LDA | pot_GGA | potpaw_LDA | potpaw_GGA | potpaw_PBE
[VASP_FORCE_OPTION]ALGO=NORMAL                        // (NORMAL | VERYFAST | FAST | ALL | DAMPED), PRESERVED (default=NORMAL)
[VASP_FORCE_OPTION]TYPE=DEFAULT                       // METAL | INSULATOR | SEMICONDUCTOR | DEFAULT (default DEFAULT)
[VASP_FORCE_OPTION]CONVERT_UNIT_CELL=SPRIM            // SPRIM, SCONV, NIGGLI, MINK, INCELL, COMPACT, WS, CART, FRAC, PRES
[AFLOW] **************************************************************************************************************************
[VASP_INCAR_MODE_EXPLICIT]START
[AFLOW]SYSTEM=Si2
NELM = 120
NELMIN=2
LPLANE=.TRUE.
LREAL=.FALSE.
LSCALU=.FALSE.
#NBANDS=XX                                            # for hand modification
#IALGO=48                                             # for hand modification
[VASP_INCAR_MODE_EXPLICIT]STOP
[AFLOW] **************************************************************************************************************************
[VASP_KPOINTS_MODE_IMPLICIT]
[VASP_KPOINTS_FILE]KSCHEME=M
[VASP_KPOINTS_FILE]KPPRA=1000
[VASP_KPOINTS_FILE]STATIC_KSCHEME=M
[VASP_KPOINTS_FILE]STATIC_KPPRA=1000
[VASP_KPOINTS_FILE]BANDS_LATTICE=FCC
[VASP_KPOINTS_FILE]BANDS_GRID=20
[AFLOW] **************************************************************************************************************************
[VASP_POTCAR_MODE_IMPLICIT]
[VASP_POTCAR_FILE]Si
[AFLOW] potpaw_PBE: Si
[AFLOW] COMPOSITION_PP=|Si2|
[AFLOW] COMPOSITION=|Si2|
[AFLOW] VOLUME(A^3)=|0|
[AFLOW] MASS(amu)=|28.0855|
[AFLOW] **************************************************************************************************************************
[AFLOW] Aflow automatically generated (aflow_avasp.cpp)
[AFLOW] **************************************************************************************************************************
[AFLOW] aflow/aconvasp/aqe/apennsy in
[AFLOW] ************************************************************************************************************************
[AFLOW] aflow/aconvasp/aqe/apennsy in
[AFLOW] ************************************************************************************************************************
[VASP_POSCAR_MODE_EXPLICIT]START
Si1   [FCC,FCC,cF8] (STD_PRIM doi:10.101
   1.22474500000000
     0.0000000000000000    2.2320315365242354    2.2320315365242354
     2.2320315365242354    0.0000000000000000    2.2320315365242354
     2.2320315365242354    2.2320315365242354    0.0000000000000000
Si
   2
Direct
  0.0000000000000000  0.0000000000000000  0.0000000000000000
  0.7500000000000000  0.7500000000000000  0.7500000000000000
[VASP_POSCAR_MODE_EXPLICIT]STOP
[AFLOW] **************************************************************************************************************************


********************************* SECTION3 END *****************************************



********************************* SECTION4 START ***************************************

#!/bin/bash
#
# These commands set up the Grid Environment for your job:
#PBS -q debug
#PBS -l nodes=1:ppn=64
#PBS -l walltime=2:00:00
#PBS -N testing
#PBS -o out.$PBS_JOBID
#PBS -e err.$PBS_JOBID
#PBS -M xx@duke.edu

AFLOW=“…AFLOW_PATH….”

source /opt/intel/bin/compilervars.sh intel64
export PATH=$PATH:/MAIN/bin/VASP:$AFLOW

cd $PBS_O_WORKDIR

# choose random directory
touch log.$PBS_JOBID
$AFLOW/aflow  --monitor >> log.$PBS_JOBID &
$AFLOW/aflow --force -multi  --machine=qrats  -D ./ >> log.$PBS_JOBID

#rename LOCK
ls -d ARUN.APL_PHONON_*/ > list; echo ./ >> list

for dir in $(<list);do
   mv  $dir/LOCK $dir/LOCK.old
done

#postprocession 

$AFLOW/aflow  --monitor >> log.$PBS_JOBID &
$AFLOW/aflow --force -multi  -rsort --machine=qrats  -D ./ >> log.$PBS_JOBID


********************************* SECTION4 END *****************************************
