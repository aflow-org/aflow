****************************************************************************************************
*                                                                                                  *
*                        aflow - STEFANO CURTAROLO Duke University 2003-2021                       *
*                            High-Throughput ab-initio Computing Project                           *
*                                                                                                  *
****************************************************************************************************
LATEST VERSION OF THE FILE:                           materials.duke.edu/AFLOW/README_AFLOW_POCC.TXT
****************************************************************************************************
One example of the high-throughput first-principle calculations for partial occupied
structure (mostly written by Kesong Yang).

1) Prepare your aflow.in input files:

[AFLOW] ************************************************************************************************************************** 
[AFLOW]SYSTEM=Ga0.25In0.75P
[AFLOW] ************************************************************************************************************************** 
[AFLOW] input file for aflow 
[AFLOW_MODE=VASP] 
[AFLOW] ************************************************************************************************************************** 
[AFLOW_MODE_ZIP=xz] 
[AFLOW_MODE_BINARY=vasp46s] 
[AFLOW] ************************************************************************************************************************** 
[AFLOW] ************************************************************************************************************************** 
#[AFLOW_MODE_MPI]
[AFLOW_MODE_MPI_MODE]NCPUS=MAX 
[AFLOW_MODE_MPI_MODE]COMMAND ="mpirun -np" 
[AFLOW_MODE_MPI_MODE]AUTOTUNE 
[AFLOW_MODE_MPI_MODE]BINARY="mpivasp46s" 
[AFLOW] ************************************************************************************************************************** 
[AFLOW_SYMMETRY]CALC 
#[AFLOW_SYMMETRY]SGROUP_WRITE 
#[AFLOW_SYMMETRY]SGROUP_RADIUS=7.77 
[AFLOW] ************************************************************************************************************************** 
#[AFLOW_NEIGHBOURS]CALC 
[AFLOW_NEIGHBOURS]RADIUS=7.7 
[AFLOW_NEIGHBOURS]DRADIUS=0.1 
[AFLOW] ************************************************************************************************************************** 
#[AFLOW_APL]CALC                                      // README_AFLOW_APL.TXT
[AFLOW_APL]ENGINE=DM                                  // README_AFLOW_APL.TXT
[AFLOW_APL]DMAG=0.015                                 // README_AFLOW_APL.TXT
[AFLOW_APL]MINATOMS=100                               // README_AFLOW_APL.TXT
#[AFLOW_APL]SUPERCELL=3x3x3                           // README_AFLOW_APL.TXT
[AFLOW_APL]DC=y                                       // README_AFLOW_APL.TXT
[AFLOW_APL]DPM=y                                      // README_AFLOW_APL.TXT
[AFLOW_APL]ZEROSTATE=y                                // README_AFLOW_APL.TXT
[AFLOW_APL]DOS=y                                      // README_AFLOW_APL.TXT
[AFLOW_APL]TP=y                                       // README_AFLOW_APL.TXT
[AFLOW_APL]TPT=0:2000:10                              // README_AFLOW_APL.TXT
[AFLOW] ************************************************************************************************************************** 
#[AFLOW_QHA]CALC                                      // README_AFLOW_QHA_SCQHA_QHA3P.TXT
[AFLOW_QHA]MODE=QHA3P                                 // README_AFLOW_QHA_SCQHA_QHA3P.TXT
[AFLOW_QHA]EOS=y                                      // README_AFLOW_QHA_SCQHA_QHA3P.TXT
[AFLOW] ************************************************************************************************************************** 
#[AFLOW_AAPL]CALC                                     // README_AFLOW_APL.TXT
[AFLOW_AAPL]TDMAG=0.015                               // README_AFLOW_APL.TXT
[AFLOW_AAPL]CUT_SHELL=4                               // README_AFLOW_APL.TXT
[AFLOW_AAPL]CUT_RAD=4.5                               // README_AFLOW_APL.TXT
[AFLOW_AAPL]SUMRULE=1E-5                              // README_AFLOW_APL.TXT
[AFLOW_AAPL]BTE=FULL                                  // README_AFLOW_APL.TXT
[AFLOW_AAPL]THERMALGRID=21x21x21                      // README_AFLOW_APL.TXT
[AFLOW_AAPL]ISOTOPE=y                                 // README_AFLOW_APL.TXT
[AFLOW_AAPL]CUMULATIVEK=y                             // README_AFLOW_APL.TXT
[AFLOW_AAPL]BOUNDARY=n                                // README_AFLOW_APL.TXT
[AFLOW_AAPL]NANO_SIZE=100                             // README_AFLOW_APL.TXT
[AFLOW_AAPL]TCT=200:700:20                            // README_AFLOW_APL.TXT
[AFLOW] ************************************************************************************************************************** 
[VASP_RUN]RELAX_STATIC_BANDS=2                        // GENERATE STATIC RELAX=N RELAX_STATIC=N STATIC_BANDS RELAX_STATIC_BANDS=N REPEAT_BANDS 
#[VASP_FORCE_OPTION]NEGLECT_NOMIX 
[VASP_FORCE_OPTION]CHGCAR=ON                         // ON | OFF (default ON)
[VASP_FORCE_OPTION]KEEPK 
#[VASP_FORCE_OPTION]KPOINTS=EVEN                      // EVEN | ODD (default none)
#[VASP_FORCE_OPTION]KPOINTS_KSHIFT_GAMMA=EVEN         // EVEN | ODD (default none)
#[VASP_FORCE_OPTION]KPOINTS_KSCHEME=MONKHORST_PACK    // MONKHORST_PACK | GAMMA (manual)
#[VASP_FORCE_OPTION]KPOINTS_GAMMA 
#[VASP_FORCE_OPTION]KPOINTS_IBZKPT 
[VASP_FORCE_OPTION]SYM=ON                             // ON | OFF  (default ON)
[VASP_FORCE_OPTION]AUTO_PSEUDOPOTENTIALS=potpaw_PBE   // pot_LDA | pot_GGA | potpaw_LDA | potpaw_GGA | potpaw_PBE | potpaw_LDA_KIN | potpaw_PBE_KIN  
[VASP_FORCE_OPTION]NBANDS                             // Estimate Bands (better than VASP)
[VASP_FORCE_OPTION]SPIN=OFF   // (ON | OFF  (default ON)), REMOVE_RELAX_1 | _2
[VASP_FORCE_OPTION]AUTO_MAGMOM=ON                    // ON | OFF (default OFF)
[VASP_FORCE_OPTION]RELAX_MODE=FORCES                  // ENERGY | FORCES | ENERGY_FORCES | FORCES_ENERGY (default ENERGY) 
[VASP_FORCE_OPTION]PREC=HIGH                          // (LOW | MEDIUM | NORMAL | HIGH | ACCURATE), PRESERVED (default=MEDIUM)
[VASP_FORCE_OPTION]ALGO=FAST                          // (NORMAL | VERYFAST | FAST | ALL | DAMPED), PRESERVED (default=NORMAL)
[VASP_FORCE_OPTION]RELAX 
#[VASP_FORCE_OPTION]NOTUNE 
[VASP_FORCE_OPTION]TYPE=INSULATOR                     // METAL | INSULATOR | SEMICONDUCTOR | DEFAULT (default DEFAULT) 
[VASP_FORCE_OPTION]CONVERT_UNIT_CELL=SPRIM       // SPRIM, SCONV, NIGGLI, MINK, INCELL, INCOMPACT, WS, CART, FRAC
#[VASP_FORCE_OPTION]VOLUME+=10.0 
#[VASP_FORCE_OPTION]VOLUME*=1.05 
[AFLOW] ************************************************************************************************************************** 
[AFLOW] ************************************************************************************************************************** 
[VASP_INCAR_MODE_EXPLICIT]START 
SYSTEM=Ga0.25In0.75P
#NELM = 120
NELMIN=2
LPLANE=.TRUE.
LREAL=.FALSE.
LSCALU=.FALSE.
PSTRESS=000                                           # for hand modification
#NBANDS=XX                                            # for hand modification
#IALGO=48                                             # for hand modification
[VASP_INCAR_MODE_EXPLICIT]STOP 
[AFLOW] ************************************************************************************************************************** 
[VASP_KPOINTS_MODE_IMPLICIT] 
[VASP_KPOINTS_FILE]KSCHEME=G 
[VASP_KPOINTS_FILE]KPPRA=1000
[VASP_KPOINTS_FILE]STATIC_KSCHEME=G
[VASP_KPOINTS_FILE]STATIC_KPPRA=3000
[VASP_KPOINTS_FILE]BANDS_LATTICE=AUTO
[VASP_KPOINTS_FILE]BANDS_GRID=20
[AFLOW] ************************************************************************************************************************** 
[AFLOW] ************************************************************************************************************************** 
[VASP_POTCAR_MODE_IMPLICIT] 
[VASP_POTCAR_FILE]Ga
[VASP_POTCAR_FILE]In
[VASP_POTCAR_FILE]P
[AFLOW] ************************************************************************************************************************** 
[AFLOW] ************************************************************************************************************************** 
[POCC_MODE_EXPLICIT]START.POCC_STRUCTURE
Ga0.25In0.75P 
1.224745  0.001
   0.00000000000000   2.24076438492221   2.24076438492221
   2.24076438492221   0.00000000000000   2.24076438492221
   2.24076438492221   2.24076438492221   0.00000000000000
1*0.25 1*0.75 1
Direct(2) P
   0.00000000000000   0.00000000000000   0.00000000000000  Ga    
   0.00000000000000   0.00000000000000   0.00000000000000  In    
   0.75000000000000   0.75000000000000   0.75000000000000  P     
[POCC_MODE_EXPLICIT]STOP.POCC_STRUCTURE
[AFLOW] ************************************************************************************************************************** 
[AFLOW_POCC]CALC



2) Once the appropriate aflow.in file is finished, you have to prepare the input with
aflow --pocc_input
This command will analyze your partially occupied structure from aflow.in file, which is between the lines:
[POCC_MODE_EXPLICIT]START.POCC_STRUCTURE
...
[POCC_MODE_EXPLICIT]STOP.POCC_STRUCTURE
In addition, the command will also generate a series of nonequivalent structures according to your partially occupied structure.
The log will be save in "LOG.POCC" file.

3) let your aflow daemon run your new aflow.in. The code will do the calculations for all the structures.
Once in a while you should check if there is an OUTCAR that has not been touched for a while. This indicates incomplete termination of VASP.
If the process gets stuck somewhere, delete the directory of the OUTCAR which is not been finished,
then delete the LOCK file of aflow, and restart aflow. The completed directories will not be modified.
The running step does not require an internet connection, so it is safe to run the displacements in a supercomputer with many nodes.

A simple trick. If you want to speed up the calculation, after aflow starts and creates its LOCK, you can rename the LOCK as LOCK1.
This  will not interfere with the running process since the code keeps a pointer to the physical file.

4) when all the directories have finished, run
4.1) aflow --pocc_dos=DIRECTORY,TEMPERATURE
Firstly, this command will analyze your aflow.in, scan through all the subdirectories, warn you if there are missing calculations.
If this is a complete calculation, then this command will generate the density of states (DOS) plot for the partially occupied structure at 
a specified temperature.

4.2) aflow --pocc_bandgap=DIRECTORY,TEMPERATURE
This command will produce the band gap of the partially occupied structure at the specified temperature. The band gap estimation comes from the 
density of states (DOS).

4.3) aflow --pocc_mag=DIRECTORY,TEMPERATURE
This command will produce the magnetic moment of the partially occupied structure at the specified temperature.

Good luck.
Kesong
July. 18th, 2013

5) Running APL with POCC

Using the APL module with POCC requires a completed POCC run.

5.1) Copy the aflow.in file for POCC into a new aflow.in file, e.g. aflow_apl.in. Inside this file,
uncomment the line "[AFLOW_APL]CALC"

5.2) Run aflow with the new aflow.in file. All subdirectories should now each contain an aflow_apl.in file.

5.3) Run APL inside the POCC subdirectories. For more information on APL, see the APL README. Running
the APL postprocessing step is not necessary.

5.4) Run aflow in the POCC directory with the aflow_apl.in file.

****************************************************************************************************
*                                                                                                  *
*                        aflow - STEFANO CURTAROLO Duke University 2003-2021                       *
*                            High-Throughput ab-initio Computing Project                           *
*                                                                                                  *
****************************************************************************************************
