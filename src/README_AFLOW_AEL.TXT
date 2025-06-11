********************************************************************************
*                                                                                                                     *
*             Aflow STEFANO CURTAROLO - Duke University 2003-2021           *
*                High-Throughput ab-initio Materials Discovery                          *
*                                                                                                                     *
********************************************************************************
LATEST VERSION OF THE FILE:        materials.duke.edu/AFLOW/README_AFLOW_AEL.TXT
********************************************************************************


Elastic properties calculated directly from stress-strain data.
This implementation in AFLOW written by Cormac Toher and Stefano Curtarolo, 
and is described in the paper Phys. Rev. Materials 1, 015401 (2017),
and is developed from the method described in Scientific Data 2, 150009 (2015).
Please cite the papers mentioned above in addition to the standard AFLOW papers 
if you use results obtained using the AFLOW AEL code.

The current version of AEL can calculate the Poisson ratio, bulk modulus, shear modulus,
compliance tensor and stiffness tensor (elastic constants). The results are written in a
variety of formats in different output files. For further details on the output files
written, see points (9) and (10) below. Other results can be written and plotted by
setting the appropriate values in the aflow.in file, as described in point (10) below.

AFLOW AEL expects a well relaxed structure on input. There is no additional relaxation
of the cell shape and size inside the AFLOW AEL, only of the ion positions. So, you
have to prepare such a structure by a separate AFLOW run using a relaxation mode for
energy and/or forces. 

When AFLOW AEL is run, it will first check for the existence of the required DFT
calculations. If these exist, then it will extract the required information and
proceed to calculate the thermal properties. If these do not exist, then it will
create the required directories and aflow.in files, and then stop. By running AFLOW
with the "--D" and "--multi" flags, AFLOW will then automatically descend into these
subdirectories and run the required electronic structure calculations, after which it
will stop. The user will then have to rename or delete the LOCK file in the main
directory, and the run AFLOW again to calculate the thermal properties. This can all
be performed within a submission script as in the following example:

aflow --multi -D $pathname/ >> log.$PBS_JOBID
mv $pathname/LOCK $pathname/LOCK.old
aflow --multi -D $pathname/ >> log.$PBS_JOBID

This procedure is similar to that required to run AFLOW APL and AFLOW AGL.

1) To perform a AFLOW AEL calculation the following line has to be present in aflow.in:

   [AFLOW_AEL]CALC

   AFLOW checks for the existence of this line in aflow.in. If it does not exist, then
   the job is not an AFLOW AEL calculation and no other AEL settings are read. 

   AEL can also be run from within other parts of AFLOW, such as from AGL, in order to
   obtain properties such as Poisson ratio to be used within these other parts of AFLOW.
   In this case, AEL will run and the other AEL settings in the aflow.in file will be
   read, even if the "[AFLOW_AEL]CALC" line is not present.

2) The AFLOW AEL method requires a set of stress vs. strain data, which are obtained
   from DFT calculations of a set of strained structures. There are two types of strains:
   normal strains and shear strains. For each strain type, there are up to three
   independent strain directions. The user has the possibility to set the number of
   strained structures in each direction, and the size of the strain step, for each
   type of strain. These parameters can be set by including the following lines in
   the aflow.in file:
    
   [AFLOW_AEL]NNORMAL_STRAINS=number
   [AFLOW_AEL]NSHEAR_STRAINS=number
   [AFLOW_AEL]NORMAL_STRAIN_STEP=number
   [AFLOW_AEL]SHEAR_STRAIN_STEP=number

   The default setting is 4 strains in each direction and a step size of 0.005. The
   default is to apply both positive and negative strains, i.e. to both stretch and
   compress the lattice. This can be controlled by including the following line in
   the aflow.in file:  

   [AFLOW_AEL]NEG_STRAINS=option

   A value of "ON" uses both positive and negative strains, while a value of "OFF"
   uses positive strains only.

3) The AFLOW AEL method can use the symmetry of the lattice to reduce the number
   of independent strains required for an elasticity calculation. For example, for
   lattices with cubic symmetry, the three orthogonal directions are equivalent so
   that there is only one independent strain direction for each type of strain. The
   number of strain directions can be set automatically by including the following
   line in the aflow.in file:

   [AFLOW_AEL]STRAIN_SYMMETRY=option

   The default is "OFF" so that all three directions are treated as independent for
   each strain type. The number of independent strain directions can also be set
   directly by the user by including the following line in the aflow.in file:

   [AFLOW_AEL]NIND_STRAINDIRS=number

   The default value for the number of strain directions is three for each type of
   strain.

   The lattice symmetry can also be used to check the elastic constant tensor
   calculated by AFLOW AEL. If the tensor does not have the correct symmetry,
   AFLOW AEL will issue a warning for the user. This can be controlled by including
   the following line in the aflow.in file:  

   [AFLOW_AEL]CHECK_ELASTIC_SYMMETRY=option   
   
   The default is "ON" so that the symmetry check is performed. The tolerance used
   for checking the symmetry can be set by including the following line in the aflow.in
   file:  

   [AFLOW_AEL]SYM_TOL_FRAC=number

   The default value is 0.05, and refers to the fractional difference allowed between
   two values which should be equal.   

4) The AFLOW AEL method allows to choose whether or not to calculate and/or fit the
   stress for the zero-strain system. The stress tensor for the zero-strain system can
   be calculated by including the following line in the aflow.in file:  

   [AFLOW_AEL]ORIGIN_STRAIN_CALC=option

   The default is "OFF" so that the zero-strain system is not calculated. The zero-strain
   stress can be included in the fitting points used to obtain the elastic constants by
   including the following line in the aflow.in file:  

   [AFLOW_AEL]ORIGIN_STRAIN_FIT=option

   The default is "OFF" so that the zero-strain system is not calculated. Note that if
   the zero-strain system is fitted without the stress being explicitly calculated,
   AFLOW AEL will use values of 0.0 for all elements of the stress tensor. This requires
   the system to be very well relaxed (including cell size and shape) prior to the AEL
   run in order to be valid.

   The stress-strain data from a previous AFLOW relaxation can also be fitted by including
   the following line in the aflow.in file:

   [AFLOW_AEL]RELAXED_STRUCT_FIT=option
   
   The default is "OFF" so that the previously relaxed system is not calculated. Note
   that this keyword requires that the results of a full AFLOW relaxation run be present
   in the directory in which is AEL calculation is being run. This keyword is useful if
   such a relaxation already exists, so that an additional point can be added to the
   stress-strain data for the fit without the need for additional DFT calculations.

5) The PREC and ALGO settings for the VASP calculations for the strained systems are set
   to ACCURATE and NORMAL by default. These values are chosen as tests of the AFLOW AEL
   method has shown that these settings produce the most reliable results. To read the
   values for PREC and ALGO from the original aflow.in file, the following line should
   be included in the aflow.in file:

   [AFLOW_AEL]PRECACC_ALGONORM=option
 
   The default is "OFF" so that the PREC and ALGO settings in the initial aflow.in file
   are ignored and the setting of PREC=ACCURATE and ALGO=NORMAL are used.  

6) The stress tensor for a particular deformed system as calculated with VASP can be
   read from either the OUTCAR file or from the vasprun.xml file. The stress tensor
   can be read from the vasprun.xml file by including the following line in the aflow.in
   file:  

   [AFLOW_AEL]VASPRUNXML_STRESS=option

   The default is "OFF" so that the stress tensor is read from the OUTCAR file.

7) The AFLOW AEL method allows the user to control what type of AFLOW runs are used to
   calculate the stress tensor for each strained structure, independently of the [VASP_RUN]
   setting contained in the original aflow.in file. The possible settings are RELAX=2,
   RELAX_STATIC=2, and STATIC, and can be set by including the following lines in the
   aflow.in file respectively:

   [AFLOW_AEL]RELAX=option
   [AFLOW_AEL]RELAX_STATIC=option
   [AFLOW_AEL]STATIC=option
   
   The default is RELAX_STATIC=ON, which runs RELAX_STATIC=2 for each strained structure. 

8) Since multiple calculations are performed for each independent direction, there is a
   certain amount of redundancy towards the failure of some of these calculations, and
   therefore options have been implemented within AEL to allow for the skipping of one
   or more failed run in each independent direction. Failed runs can be detected
   automatically (by the absence of the appropriate LOCK or OUTCAR file in the directory),
   and the option of whether or not to skip these runs can be set by including the
   following line in the aflow.in file:

   [AFLOW_AEL]AUTOSKIP_FAILED_ARUNS=option

   The default is "OFF" so that no runs are skipped automatically; instead detection of
   a failed run will cause AEL to exit with an error unless this particular run was
   explicitly listed to be skipped.

   Specific runs can be listed for skipping by including the following line in the aflow.in
   file:

   [AFLOW_AEL]SKIP_FAILED_ARUNS=directory_name_list

   where "directory_name_list" is a string with the list of directory names containing the
   failed runs to be skipped by the AEL fitting algorithm. If there is more than one failed
   run, the individual directory names should be listed on the same line, separated by commas.

   AEL only permits a certain number of runs for any given independent strain direction to
   be skipped. The maximum number of runs permitted to be skipped in any strain direction
   can be set by including the following line in the aflow.in file:

   [AFLOW_AEL]SKIP_ARUNS_MAX=number

   The default value is 1, so that only one run can be skipped per independent direction.
   Note that there is also a hard limit on the number of runs in a given direction which
   can be skipped which is equal to half of the initial number of runs set up in that
   direction.

9) The main output file for the AEL method is aflow.ael.out. It contains the AEL results
   for the Poisson ratio, bulk modulus (Voigt, Reuss, and Voigt-Reuss-Hill (VRH)), shear
   modulus (Voigt, Reuss, and Voigt-Reuss-Hill (VRH)), speed of sound, AEL Debye temperature,
   and elastic anisotropy. It also contains the elastic constants in the form of the stiffness
   tensor and the compliance tensor.
   Results are also written in JSON format in three files: AEL_elastic_tensor.json and
   AEL_energy_structures.json. The file AEL_elastic_tensor.json contains the elastic
   constants in the form of the individual elements of the stiffness and compliance
   tensors. The file AEL_energy_structures.json contains the energies, stresses and pressures
   for each calculated structure, as well as the structures themselves and the distortions
   applied to them.

10) The AFLOW AEL method allows the user to choose whether or not to automatically write
    the directional elastic moduli (Young's moduli, etc.) to a file. The writing of this
    file can be controlled by including the following entry in the aflow.in file:

    [AFLOW_AEL]WRITE_FULL_RESULTS=option

    The default is "OFF", so that these tensors are not written.

    Finally, the values of the Poisson ratio, bulk modulus, shear modulus and elastic
    anisotropy are written in a script-parsable format in the file aellib.dat. The format
    used is the same as that used for the file aflowlib_entry.out. The elastic constants
    (in the form of the stiffness and compliance tensor) are also written to files
    independently of the value of the WRITE_FULL_RESULTS setting.

11) The elastic properties can be calculated at different pressures using AEL, which is
    useful for determining the pressure range in which a material is mechanically stable.
    This calculation is controlled by including the following line in the aflow.in file:

    [AFLOW_AEL]PRESSURE_CALC=option

    The default is "OFF", so that the elastic properties are only calculated at zero
    pressure.

    If this option is set to "ON", AEL will call AGL to run a set of different volume
    calculations, to determine the pressure as a function of cell volume. This data is
    then used to determine the cell volumes corresponding to the required pressures.
    These pressures are provided by using the AGL pressure stepsize and number of
    pressure steps options:

    [AFLOW_AGL]SPRESSURE=number
    [AFLOW_AGL]NPRESSURE=number

    Specific pressures can also be specified for the AEL pressure calculations,
    independently of the AGL calculations, by including the following line in the
    aflow.in file:

    [AFLOW_AEL]SPRESSURE=number
    [AFLOW_AEL]NPRESSURE=number

    Using this option allows a different set of pressures to be using for calculating
    the elastic constants than for the final AGL postprocessing, within a single
    ingetrated workflow.

    Once the appropriate volumes are determined from the AGL calculations, AEL then
    generates the appropriate normal and shear strained structures, and calculates
    the stress tensors with DFT. The elastic properties at different pressures are
    then written in the file aflow.ael_pressure.out. The output data includes a
    statement of whether or not the material is mechanically stable at each pressure.

    Note that this method will only guarantee that the cell pressure will only have
    the set value prior to any relaxation performed. Relaxing the ion positions, as
    is the default for AEL, could result in a change in the cell pressure. Using the
    STATIC option will prevent any relaxation, and will ensure that the set pressure
    is obtained, although this will result in the "clamped ion" elastic tensor being
    obtained.

    Elastic properties calculated at a finite pressure where the materials is stable
    can be used in other parts of the code by including the following line in the
    aflow.in file:

    [AFLOW_AEL]USE_STABLE_PRESSURE=option

    If this option is set to "ON", then the elastic properties at the lowest
    calculated pressure where the material is elastically stable will be
    set as the main elastic properties for this material: these values
    will be written to the main aflow.ael.out file, this value of the
    Poisson ratio will be used if requested by AGL, etc. The default value
    is "OFF", so that the elastic properties calculated at zero pressure
    are used.

12) The elastic tensor can be symmetrized after calculation. This option can be selected
    by including the following line in the aflow.in file:

    [AFLOW_AEL]SYMMETRIZE=option

    The default is "OFF", so that the elastic tensor is not symmetrized. If the elastic tensor
    is asymmetric, AEL will give a warning, and the user may want to symmetrize the tensor.
    If the option is set to "ON", then the tensor is symmetrized according to the following
    equation:

       sc_{ij} = (c_{ij} + c_{ji}) / 2

    where c_{ij} are the elements of the unsymmetrized tensor, and sc_{ij} are the elements of
    the symmetrized tensor.

13) The symmetry for the VASP calculations for the strained structures can be switched on and off
    by including the following line in the aflow.in file:

    [AFLOW_AEL]VASPSYM=option

    The default is "OFF", so that the VASP symmetry is switched off (sets SYM=0 in the INCAR file).
    This is because the VASP symmetry routines can have problems with the distorted cells used for
    calculating the elastic tensor. Setting this option to "ON" switches on the VASP symmetry. 

14) EXAMPLE:  CaO FCC with default AEL settings (note that this will give the same results if all of the "[AFLOW_AEL]" 
   lines except "[AFLOW_AEL]CALC" are removed or commented out):

[AFLOW] **************************************************
[AFLOW] Test file for AEL calculation of CaO (fcc)
[AFLOW] **************************************************
[AFLOW] ************************************************************************************************************************** 
[AFLOW] input file for aflow 
[AFLOW_MODE=VASP] 
[AFLOW] ************************************************************************************************************************** 
[AFLOW_MODE_ZIP=xz] 
[AFLOW_MODE_BINARY=vasp46s] 
[AFLOW] ************************************************************************************************************************** 
[AFLOW] ************************************************************************************************************************** 
[AFLOW_MODE_MPI]
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
[AFLOW_APL]DISMAG=0.015                               // README_AFLOW_APL.TXT
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
#[VASP_RUN]RELAX_STATIC_BANDS=2                        // GENERATE | STATIC | RELAX=N | RELAX_STATIC=N | STATIC_BANDS | RELAX_STATIC_BANDS=N | REPEAT_BANDS [,DS[,DD[,DSCF]]] 
#[VASP_FORCE_OPTION]NEGLECT_NOMIX 
[VASP_FORCE_OPTION]CHGCAR=OFF                         // ON | OFF (default ON)
#[VASP_FORCE_OPTION]KPOINTS=keyword[,keyword]         // EVEN | ODD | KSHIFT_GAMMA_EVEN | KSHIFT_GAMMA_ODD | KSCHEME_MONKHORST_PACK | KSCHEME_GAMMA | GAMMA | KEEPK | IBZKPT
[VASP_FORCE_OPTION]SYM=ON                             // ON | OFF  (default ON)
[VASP_FORCE_OPTION]AUTO_PSEUDOPOTENTIALS=potpaw_PBE   // pot_LDA | pot_GGA | potpaw_LDA | potpaw_GGA | potpaw_PBE | potpaw_LDA_KIN | potpaw_PBE_KIN  
[VASP_FORCE_OPTION]NBANDS                             // Estimate Bands (better than VASP)
#[VASP_FORCE_OPTION]SPIN=ON,REMOVE_RELAX_1             // (ON | OFF  (default ON)), REMOVE_RELAX_1 | _2
#[VASP_FORCE_OPTION]AUTO_MAGMOM=ON                    // ON | OFF (default OFF)
[VASP_FORCE_OPTION]RELAX_MODE=ENERGY                  // ENERGY | FORCES | ENERGY_FORCES | FORCES_ENERGY (default ENERGY) 
[VASP_FORCE_OPTION]PREC=HIGH                          // (LOW | MEDIUM | NORMAL | HIGH | ACCURATE), PRESERVED (default=MEDIUM)
[VASP_FORCE_OPTION]ALGO=FAST                          // (NORMAL | VERYFAST | FAST | ALL | DAMPED), PRESERVED (default=NORMAL)
#[VASP_FORCE_OPTION]NOTUNE 
[VASP_FORCE_OPTION]TYPE=DEFAULT                       // METAL | INSULATOR | SEMICONDUCTOR | DEFAULT (default DEFAULT) 
[VASP_FORCE_OPTION]CONVERT_UNIT_CELL=SPRIM            // SPRIM, SCONV, NIGGLI, MINK, INCELL, COMPACT, WS, CART, FRAC, PRES 
#[VASP_FORCE_OPTION]VOLUME+=10.0 
#[VASP_FORCE_OPTION]VOLUME*=1.05 
[AFLOW] ************************************************************************************************************************** 
[AFLOW] ************************************************************************************************************************** 
[VASP_INCAR_MODE_EXPLICIT]START 
SYSTEM=Ca1O1_ICSD_26959
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
[VASP_KPOINTS_FILE]KPPRA=8000
[VASP_KPOINTS_FILE]STATIC_KSCHEME=M 
[VASP_KPOINTS_FILE]STATIC_KPPRA=10000
[VASP_KPOINTS_FILE]BANDS_LATTICE=FCC
[VASP_KPOINTS_FILE]BANDS_GRID=20
[AFLOW] ************************************************************************************************************************** 
[AFLOW_AEL]CALC
[AFLOW_AEL]NNORMAL_STRAINS=4
[AFLOW_AEL]NSHEAR_STRAINS=4
[AFLOW_AEL]NORMAL_STRAINSTEP=0.005
[AFLOW_AEL]SHEAR_STRAINSTEP=0.005
[AFLOW_AEL]NEG_STRAINS=ON
[AFLOW_AEL]NIND_STRAINS=3
[AFLOW_AEL]STRAIN_SYMMETRY=OFF
[AFLOW_AEL]CHECK_ELASTIC_SYMMETRY=OFF
[AFLOW_AEL]SYMTOLFRAC=0.05
[AFLOW_AEL]WRITE_FULL_RESULTS=OFF
[AFLOW_AEL]ORIGIN_STRAIN_CALC=OFF
[AFLOW_AEL]RELAXED_STRUCT_FIT=OFF
[AFLOW_AEL]ORIGIN_STRAIN_FIT=OFF
[AFLOW_AEL]VASPRUNXML_STRESS=OFF
[AFLOW_AEL]PRECACC_ALGONORM=ON
[AFLOW_AEL]SYMMETRIZE=OFF
[AFLOW_AEL]VASPSYM=ON
[AFLOW] ************************************************************************************************************************** 
[VASP_POSCAR_MODE_EXPLICIT]START 
Ca1O1 #225 - (Ca1O1_ICSD_26959) - Ca1O1 [Ca1O1] cF8 Fm-3m CaO 225 Ca1O1_ICSD_26959 ICSD_26959  (icsd library) (WICKOFF 225 #225)  [FCC,FCC,cF8] (STD_PRIM doi:10.1016/j.commatsci.2010.05.010)
1.224745
   0.00000000000000   1.95061033183634   1.95061033183634
   1.95061033183634   0.00000000000000   1.95061033183634
   1.95061033183634   1.95061033183634   0.00000000000000
1 1 
Direct(2) [A1B1] 
   0.00000000000000   0.00000000000000   0.00000000000000  Ca    
   0.50000000000000   0.50000000000000   0.50000000000000  O     
[VASP_POSCAR_MODE_EXPLICIT]STOP 
[AFLOW] ************************************************************************************************************************** 
[VASP_POTCAR_MODE_IMPLICIT] 
[VASP_POTCAR_FILE]Ca
[VASP_POTCAR_FILE]O
[AFLOW] potpaw_PBE: Ca_sv O 
[AFLOW] COMPOSITION_PP=|Ca_sv1|O1|
[AFLOW] COMPOSITION=|Ca1|O1|
[AFLOW] VOLUME(A^3)=|0|0|
[AFLOW] MASS(amu)=|40.08|15.9994|
[AFLOW] ************************************************************************************************************************** 
[AFLOW] Aflow automatically generated (aflow_avasp.cpp) 
[AFLOW] ************************************************************************************************************************** 
[AFLOW] aflow/aconvasp/aqe/apennsy v(30793) in ICSD/FCC/Ca1O1_ICSD_26959
[AFLOW] ************************************************************************************************************************** 
[AFLOW] SELF-MODIFICATION 
[AFLOW] Recycling CONTCAR of relax1 
[AFLOW] ************************************************************************************************************************** 
[VASP_POSCAR_MODE_EXPLICIT]START
Ca1O1 #225 - (Ca1O1_ICSD_26959) - Ca1O1 
1.224745
   0.00000000000000   1.97560216103176   1.97560216103176
   1.97560216103176   0.00000000000000   1.97560216103176
   1.97560216103176   1.97560216103176   -0.00000000000000
1 1 
Direct(2) [A1B1] 
   0.00000000000000   -0.00000000000000   -0.00000000000000 
   0.50000000000000   0.50000000000000   0.50000000000000 
[VASP_POSCAR_MODE_EXPLICIT]STOP
[AFLOW] ************************************************************************************************************************** 
[VASP_FORCE_OPTION]SPIN=OFF      // Self Correction
[AFLOW] SELF-MODIFICATION 
[AFLOW] Recycling CONTCAR of relax2 
[AFLOW] ************************************************************************************************************************** 
[VASP_POSCAR_MODE_EXPLICIT]START
Ca1O1 #225 - (Ca1O1_ICSD_26959) - Ca1O1 
1.224745
   0.00000000000000   1.97560216103176   1.97560216103176
   1.97560216103176   0.00000000000000   1.97560216103176
   1.97560216103176   1.97560216103176   0.00000000000000
1 1 
Direct(2) [A1B1] 
   0.00000000000000   0.00000000000000   0.00000000000000 
   0.50000000000000   0.50000000000000   0.50000000000000 
[VASP_POSCAR_MODE_EXPLICIT]STOP
[AFLOW] ************************************************************************************************************************** 



