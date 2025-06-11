****************************************************************************************************
*                                                                                                  *
*                        aflow - STEFANO CURTAROLO Duke University 2003-2021                       *
*                            High-Throughput ab-initio Computing Project                           *
*                                                                                                  *
****************************************************************************************************
LATEST VERSION OF THE FILE:                            materials.duke.edu/AFLOW/README_AFLOW_APL.TXT
****************************************************************************************************

Phonon calculations with deformations written by Marco Esters, Michal Jahnatek, Andriy Smolyanyuk,
Corey Oses, Pinku Nath, Jose J. Plata, Stefano Curtarolo.

For References, please see section 9 at the end of the file.

The current version of the Automatic Phonon Library (APL) can calculate the phonon dispersion
curves, phonon density of states, and thermal properties. All are calculated in the harmonic
approximation. Two methods of calculation can be selected: the direct (supercell) method (DM) and
the linear response method (LR).

The Quasi-Harmonic Approximation (QHA) performs APL calculations at several volumes to calculate a
more accurate temperature-dependent free energy. It can be used to calculate thermomechanical
properties such as thermal expansion, the bulk modulus, and heat capacities.

The Automatic Anharmonic Phonon Library (AAPL) is an expansion of APL and can calculate the lattice
thermal conductivity of a material via the Boltzmann Transport Equation using anharmonic
interatomic force constants. 

The default values of APL, QHAm and AAPL can be changed by editing the appropriate values in the
aflow.rc file.

----------------------------------------------------------------------------------------------------

Table of Contents:

1.    Running APL

2.    Choosing a Phonon Calculation Method
2.1.    Common Parameters Used by Both Methods
2.1.1.    The Precision Parameter
2.1.2.    Supercell Parameters
2.1.3.    k-Point Parameters
2.1.4.    Polar Corrections
2.1.5.    Number of Bands
2.2.    Direct Method Parameters
2.3.    Linear Response Method Parameters

3.    Calculating Phonon Properties
3.1.    Phonon Dispersion Curves
3.2.    Phonon Density of States
3.3.    Thermodynamic Properties
3.4.    Plotting Phonon Properties
3.5.    Atomic Displacements
3.5.1.    Mean Square Displacements
3.5.2.    Visualizing Phonon Displacements

4.    Saving and Loading APL Results (Hibernation)

5.    Inheritance from AFLOW

6.    Example APL Input File

7.    The Quasi-Harmonic Approximation
7.1.    Running QHA Calculations
7.2.    QHA Modes
7.3.    Equation of State Calculations
7.4.    Plotting Thermophysical Properties
7.5.    Example QHA Input File

8.    The Automatic Anharmonic Phonon Library
8.1     Running AAPL Calculations
8.2.    Anharmonic Interatomic Force Constants
8.2.1.    Cluster Parameters
8.2.2.    Sum Rule Parameters
8.2.3.    k-Point Parameters
8.3.    Lattice Thermal Conductivity Calculations
8.4.    Plotting Thermal Conductivity
8.5.    Example AAPL Input File

9. References

----------------------------------------------------------------------------------------------------

1. Running APL

To perform a phonon calculation, the following line has to be present in the aflow.in file:

    [AFLOW_APL]CALC

If it does not exist, the job is not interpreted as a phonon calculation and no other phonon
settings are read.

The calculation of phonon-related properties requires very accurate forces. To achieve this, the
input structure needs to be relaxed with very tight convergence criteria for the total energy and
forces. This can be achieved in one of two ways:

  * (VASP only) By setting the following flags:

      [AFLOW_APL]RELAX=ON
      [AFLOW_APL]PREC=PHONONS

    These are also the default values. Setting RELAX=ON makes APL run two relaxation calculations
    with the precision given in PREC. This method is especially useful when a relaxed structure
    already exists. The input and output files are saved in XXXCAR.relax_aplN with N = 1,2.
    For more information on the PREC option, see Section 2.1.1.

    APL can automatically choose the k-point grid for the relaxation to be commensurate with the
    k-point grid for the supercell. For example, if a 3x2x1 supercell uses a 2x2x3 k-point grid,
    the unit cell will be relaxed using a 6x4x3 k-point grid. This minimizes the risk of ghost
    forces in the supercell calculations. This feature can be turned on by using:

      [AFLOW_APL]RELAX_COMMENSURATE=ON

    The default is ON. If this option is OFF, APL will use the k-point settings for AFLOW RELAX
    runs. See README_AFLOW.TXT for more information.

    NOTE: Relaxations will be skipped when a PHPOSCAR file is present (see below).

  * By relaxing the structure in a normal AFLOW relaxation calculation using the flag:

      [VASP_FORCE_OPTION]PREC=PHONONS

    APL can be run afterwards with the resulting structure. [AFLOW_APL]RELAX should be set to OFF.

APL saves the structure that is used to create the supercell in a VASP-5-formatted POSCAR file
named PHPOSCAR. This structure will be loaded first, no matter what is found in the aflow.in file.

2. Choosing a Phonon Calculation Method

The user has the possibility to choose two methods of calculation using the keyword ENGINE: the
direct method (DM), also known as supercell method, and the linear response method (LR). To 
choose the DM, one has to write:

    [AFLOW_APL]ENGINE=DM

This is also a default setting. For the linear response method, it is:

    [AFLOW_APL]ENGINE=LR

Each method has its own parameters, but there are also some common ones that can be used to tune
the whole calculation.

2.1. Common Parameters Used by Both Methods

2.1.1. The Precision Parameter

The precision for the phonon calculations can be set using the flag:

    [AFLOW_APL]PREC=PHONONS

This parameter is similar [VASP_FORCE_OPTION]PREC and guides the precision of the relaxation
calculation when [AFLOW_APL]RELAX is set to ON as well as all APL and AAPL calculations. The
default is PHONONS. The setting of [VASP_FORCE_OPTION]PREC has no influence on the APL and AAPL
calculation.

APL also introduces the new option "PHONONS" for the PREC parameter, which can be used for both
the [AFLOW_APL]PREC and the [VASP_FORCE_OPTION]PREC flag. It adds the following lines in the
INCAR file:

PREC = ACCURATE
LREAL = .FALSE.
EDIFF = 1E-8
EDIFFG = -0.001
NELMIN = 4
NELM = 120
ADDGRID = .TRUE.
ALGO = NORMAL

NOTE: PREC=PHONONS overrides the [VASP_FORCE_OPTION]ALGO flag.

2.1.2. Supercell Parameters

The most important setting is the supercell dimensions (a supercell is required by both methods).
The supercell can be set by the following commands (in order of precedence):

    [AFLOW_APL]SUPERCELL=numberXnumberXnumber
    [AFLOW_APL]MINATOMS_UNIFORM=number
    [AFLOW_APL]MINATOMS=number
    [AFLOW_APL]MINSHELL=number

Examples:

    * Direct specification of the supercell dimensions, e.g.:
            
        [AFLOW_APL]SUPERCELL=4x4x4

    * Require a minimum number of atoms, e.g.:

        [AFLOW_APL]MINATOMS=100

      The code will homogeneously increase the volume of the supercell until it finds a
      SUPERCELL(i,j,k) containing at least the number of atoms specified.

      If it is desired that the supercell is expanded equally in all directions, i.e. to use
      SUPERCELL(i,i,i), use the flag MINATOMS_UNIFORM as in:

        [AFLOW_APL]MINATOMS_UNIFORM=100

       MINATOMS_UNIFORM takes precedence over MINATOMS.

     * Specifying the minimal number of nearest neighbor shells, e.g.:

        [AFLOW_APL]MINSHELL=6

      The default settings is MINSHELL=6 and APL will construct the appropriate supercell. HOWEVER,
      THIS IS NOT THE RECOMMENDED WAY BECAUSE THE INCLUSION OF 6 SHELLS IS SUFFICIENT ONLY FOR VERY
      SIMPLE MATERIALS. There is unfortunately no rule on how to find the ideal number of shells
      for the calculation.

If no option is specified, AFLOW will use the default value for MINATOMS. For best results, it is
recommended to specify the supercell dimensions directly. We found that 4x4x4 is usually enough for
all fcc systems; 3x3x2 for hcp; bcc cells are complicated, but generally, it should be more that
3x3x3. Typically, the supercell with about 100 atoms should be fine for the majority of crystal
systems. The MINATOMS option might be useful for large unit cells.

2.1.3. k-Point Parameters

The k-point mesh can be set separately for APL calculations. The settings are similar to the
general AFLOW k-point parameters:

    * The number of k-points per reciprocal atoms can be set via:

        [AFLOW_APL]KPPRA=number

      By default, AFLOW takes the number of KPOINTS (KPPRA) from the "[VASP_KPOINTS_FILE]KPPRA"
      entry. The "[AFLOW_APL]KPPRA" flag overrides this for APL calculations. See README_AFLOW.TXT
      for more information on KPPRA.

    * The k-point grid can also be set directly using:

        [AFLOW_APL]KPOINTS_GRID=k1xk2xk3

      This takes precedence over KPPRA. KPPRA is used as default.

    * The k-point scheme can be set using:

        [AFLOW_APL]KSCHEME=scheme

      By default, AFLOW takes the KPOINTS scheme from the "[VASP_KPOINTS_FILE]KSCHEME" entry. The
      "[AFLOW_APL]KSCHEME" overrides this setting. See README_AFLOW.TXT for more information on
      KSCHEME.

    * The k-point shift vector can be set using the following format:

        [AFLOW_APL]KPOINTS_SHIFT=X,Y,Z

      Where X, Y, and Z are the comma-separated components of the shift vector.

    * The k-point parity can be set using:

        [AFLOW_APL]KPOINTS=EVEN | ODD

      By default, AFLOW takes the KPOINTS parity from "[VASP_FORCE_OPTION]KPOINTS" entry, if
      specified. "[AFLOW_APL]KPOINTS" overrides this setting. See README_AFLOW.TXT for more
      information on KPOINTS.

These k-point settings will also be used for the structural relaxations (see Section 1).

2.1.4. Polar Corrections

If the studied system is polar, dipole-dipole interactions need to be taken into account. This
can be specified by the keyword POLAR. This keyword can be set in the aflow.in file as follows:

    [AFLOW_APL]POLAR=ON

This instructs APL to calculate the Born effective charge tensor and the dielectric constant
tensor using density functional perturbation theory (DFPT). The default setting is POLAR=OFF, i.e.
no polar calculations are performed.

DFPT calculations require VASP 5.2 or higher. APL will check if the appropriate binary file is
present upon start-up, and will refuse to continue if the appropriate VASP version is not found.
They also require large stack sizes, so it is critical that the stack size is set to unlimited.
Otherwise, VASP will likely crash at some point during the calculation. To set the stack size to
unlimited, the following line needs to be written either in the .bashrc file or the job script:

    ulimit -s unlimited

The dielectric constant tensor and the Born effective charge tensor is by default calculated using
the tag LEPSILON in the VASP INCAR file, but can be calculated using LCALCEPS as well. To switch to
LCALCEPS, the value DEFAULT_APL_USE_LEPSILON in the aflow.rc file needs to be set to false (0).

2.1.5. Number of Bands

If NBANDS is set in the parent aflow.in file, APL will scale the number of bands automatically in
the supercell file. For example, if NBANDS = 20 and a 2x2x2 supercell is used, NBANDS will be set
to 160 in the supercell calculations.

2.2. Direct Method Parameters

The direct method (DM) requires the calculation of the forces for each unique atom in the primitive
cell along three independent directions. However, some of these directions are related by symmetry
and it is thus unnecessary to perform DFT calculations for all three distortions. For example, if
one distortion direction can by symmetry produce three linearly independent vectors, the atom only
needs to be distorted in one direction. The remaining two can be recovered by the symmetry of the
system as is the case for an fcc lattice. APL will try to determine these distortions to reduce
the number of DFT calculations that need to be performed.

Generally, APL will distort an atom along the [100], [010], [001], [110], [101], [011], and the
[111] direction. It will apply all site symmetry operations, and the direction that produces the
highest number of new independent vectors will by used for the DFT calculations.

To filter out potential anharmonic effects from the forces, APL distorts the atom into the positive
and negative direction of the chosen distortions. The final force acting on each is then calculated
as 0.5 * (f(+) - f(-)).

There are multiple parameters that control how the distortions for the direct method are generated.

    * The distortion magnitude in Angstrom is set via the DMAG entry, e.g.:

        [AFLOW_APL]DMAG=0.015

      The default is 0.015 Angstrom, which is generally accepted to be sufficient for all systems.

    * If the atoms shall only be distorted along the cell axes ([100], [010], [001]), the DXYZONLY
      flag needs to be set to true, i.e.:

        [AFLOW_APL]DXYZONLY=ON

      The default is DXYZONLY=OFF, which will also search along the face and body diagonals.
      DXYZONLY is particularly useful for cells where one primitive cell vector is much longer than
      the others.

    * To prevent the reduction of distortions by factor groups (inequivalent atoms), the DINEQUIV_ONLY entry needs to be set to false:

        [AFLOW_APL]DINEQUIV_ONLY=OFF

      In this case, APL will generate VASP runs for distortions of ALL atoms.
    
    * To prevent the reduction of distortions by site symmetry, the DSYMMETRIZE entry needs to be set to false:

        [AFLOW_APL]DSYMMETRIZE=OFF

      In this case, APL will generate VASP runs for three distortion directions per site.

    * To prevent the automatic generation of negative distortions, the DPM entry needs to be set to:

        [AFLOW_APL]DPM=AUTO

      In this case, APL will determine whether the negative distortions are required or if they
      can be recovered by symmetry. This will be decided on a per-site basis. For DPM=OFF, only
      positive distortions will be created.

If the primitive structure is not well relaxed or cannot be relaxed within a given precision, the
unstrained state has some non-zero forces. In this case, it is possible to subtract these forces
from all distorted states. This can be done using the keyword:

    [AFLOW_APL]ZEROSTATE=ON

If activated, APL will generate a calculation of the undistorted supercell and the obtained forces
will be subtracted from all other strained calculations. The default is ZEROSTATE=OFF.
WARNING: This method is not well tested yet and should be used with caution.
NOTE: Using DPM=ON cancels out zerostate forces, so DPM=ON and ZEROSTATE=ON is redundant. However,
      the zerostate calculation can still be used to test for the presence of ghost forces.

2.3.  Linear Response Method Parameters

NOTE: the linear response method is currently untested!

The use of the linear response method as implemented in VASP requires VASP version 5.2 or higher.
APL will check if the appropriate binary file is present upon start-up, and will refuse to continue
if the appropriate VASP version is not found.

LR calculations require large stack sizes, so it is critical that the stack size is set to
unlimited. Otherwise, VASP will likely crash at some point during the calculation. To set the stack
size to unlimited, the following line needs to be written either in the .bashrc or the job script:

    ulimit -s unlimited

In contrast to the DM method, APL is not responsible for finding the appropriate set of distortion
vectors for LR calculations as it is all done by VASP. Hence, long run DFT calculations are to be
expected. Please note that VASP is using the mode XYZONLY so that more distortions are generated
than in the DM method.

3. Calculating Phonon Properties

APL can calculate phonon dispersion curves, phonon density of states, and various thermodynamic
properties. They are activated by the following keywords:

    [AFLOW_APL]DC=ON
    [AFLOW_APL]DOS=ON
    [AFLOW_APL]TP=ON

The default is 'ON' for all. To switch off these calculations, the parameter needs to be set to
'OFF'. Each calculation mode has its own set of parameters, which are described in the following
sections.

3.1. Phonon Dispersion Curves

The calculated dispersion curves are saved in the file aflow.apl.phonon_dispersion.out, where the
1st column is the ID of points in the sub-path, the 2nd is the position in the whole path, and the
remaining columns are the energies/frequencies of the phonon branches.

The phonon dispersion curve calculations have three adjustable parameters: the dispersion path,
the point density, and the frequency unit.

    * There are two ways to set the path of the phonon dispersion curve: by using the default path
      of the Bravais lattice of the cell as specified in [doi=10.1016/j.commatsci.2010.05.010], or
      by specifying the path manually. The keyword to guide this choice is:

        [AFLOW_APL]DCPATH=LATTICE|MANUAL

      LATTICE (default) chooses the default path, points and labels of the lattice. A different
      path can be chosen using the following format.

             [AFLOW_APL]DCUSERPATH=G-X|X-U|K-G|G-L

      The labels must correspond to the labels for the lattice as set by the AFLOW standard.

      NOTE: DCPATH=LATTICE will project the phonon dispersions to the primitive cell, so the
      labels must correspond to a high-symmetry point in the primitive lattice.
      
      MANUAL lets the user specify which path and points to use. Three parameters need
      to be set to use MANUAL.

      1. The coordinates of each path end point can be set in fractional or Cartesian coordinates
         using the DCINITCOORDSFRAC and DCINITCOORDSCART keywords, respectively:

             [AFLOW_APL]DCINITCOORDSFRAC=0,0,0;0.5,0.5,0.5;...
             [AFLOW_APL]DCINITCOORDSCART=0,0,0;0.5,0.5,0.5;...

         Each coordinate triplet has to be separated by a semicolon. Fractional coordinates take
         priority over Cartesian coordinates.

      2. The labels of each end point need to be defined using DCINITCOORDSLABELS as in:

             [AFLOW_APL]DCINITCOORDSLABELS=G,L,M,X...

         Where the first label corresponds to the first coordinate triplet in DCINITCOORDSFRAC or
         DCINITCOORDSCART. The labels have to be separated by a comma, not a semicolon.

      3. Finally, the path can be specified with the aforementioned labels using DCUSERPATH.

    * The number of points for each sub-path is set with the DCPOINTS keyword, e.g.:

        [AFLOW_APL]DCPOINTS=100

    * For the frequency, a variety of options are available. They can be set as:

        [AFLOW_APL]FREQFORMAT=format

      Where format needs to be replaced with the appropriate format specifier. These are HERTZ,
      THZ, MEV, and RECIPROCAL_CM or CM-1.

      These frequencies do not include the factor 2 pi. To include it, OMEGA needs to be added
      to the format specifier, e.g.:

        [AFLOW_APL]FREQFORMAT=THZ|OMEGA

      A colon (:), a semicolon(;), or a comma (,) may also be used as delimiters.

      At last, it can be specified whether imaginary frequencies should be set to 0.0 or set to a
      negative value. If they are set to negative values, ALLOW_NEGATIVE needs to be added to the
      format specifier, e.g.:

        [AFLOW_APL]FREQFORMAT=THZ|ALLOW_NEGATIVE

      If imaginary frequencies should be set to 0.0, ALLOW_NEGATIVE must not be included, e.g.:

        [AFLOW_APL]FREQFORMAT=THZ

      Sets all imaginary frequencies to 0.0 while all real frequencies are given in THz.
      The ALLOW_NEGATIVE specifier may be combined with OMEGA, e.g. in:

        [AFLOW_APL]FREQFORMAT=THZ|OMEGA|ALLOW_NEGATIVE

      The default specifier is THZ|ALLOW_NEGATIVE.

3.2.  Phonon Density of States

The phonon density of states (pDOS) are saved into the file aflow.apl.phonon_dos.out. Besides the
pDOS, this file also contains the frequencies in THz, the wavenumbers in inverse cm, and the
energies in meV.

There are multiple parameters for the pDOS: the Brillouin zone integration method, the k-point
mesh, the number of energy points, and the smearing width.

    * There are two implementations for the integration of the Brillouin zone: the Linear Tetrahedron
      method (LT) and the Root Sampling method (RS). They can be set as:

        [AFLOW_APL]DOSMETHOD=LT(RS)

      Default is DOSMETHOD=LT.

    * The mesh of q-points used for numerical integration should be set as:

        [AFLOW_APL]DOSMESH=21x21x21

      This is also the default setting. APL will automatically find the irreducible part of the
      of BZ (the list of irreducible q-points with appropriate summation weights). Only for these
      q-points is a dynamical matrix constructed and diagonalized. For the LT method, the list of
      irreducible tetrahedrons with their weights is constructed and used in summation.

    * The number of points in the pDOS can be set via the keyword DOSPOINTS, e.g.:

        [AFLOW_APL]DOSPOINTS=2000

      This is also the default setting.

    * The pDOS can be smeared using a Gaussian function. To specify the width of this function, the
      DOSSMEAR tag can be set, e.g.:

        [AFLOW_APL]DOSSMEAR=0.05

      This is especially important for the RS method. The default is 0.05 for the RS method and
      0.0 for the LT method (i.e. no smearing).

Projected DOS can be calculated using the keyword DOS_PROJECT:

         [AFLOW_APL]DOS_PROJECT=ON

The default is OFF. Projections along specific directions can be added using DOSPROJECTIONS:

        [AFLOW_APL]DOSPROJECTIONS_CART=1,0,0;0,1,0;0,0,1
                                 or
        [AFLOW_APL]DOSPROJECTIONS_FRAC=-1,1,1;1,-1,1;1,1,-1

Where DOSPROJECTIONS_CART and DOSPROJECTIONS_FRAC give the projection directions in Cartesian and
fractional coordinates, respectively. It is not possible to choose both options. If neither option
is populated, atom-projected DOS will be calculated.

Each projection vector triplet needs to be separated by a semicolon, and the coordinates inside the
triplets must be separated by commas. The projection vectors have to be given in Cartesian
coordinates. Per default, both are turned off.

Projected DOS are currently only available for the linear tetrahedron method.

3.3. Thermodynamic Properties

APL can calculate a variety of thermodynamic properties. These are the zero point vibrational
energy (U0), the internal energy (U), the vibrational free energy (Fv), the vibrational
entropy (Sv), and the isochoric specific heat (cv) for a given temperature range. The results are
written into the file aflow.apl.thermodynamic_properties.out.

The temperature range can be set as:

        [AFLOW_APL]TPT=TSTART:TEND:TSTEP

Where TSTART and TEND are the initial and final temperature, respectively, and TSTEP is the
temperature step size. The default setting is 0:2000:10, i.e. it calculates the properties from 0 K
to 2000 K in 10 K steps.

APL will also output mean square displacements for each atom inside that temperature range.

3.4. Plotting Phonon Properties

AFLOW can plot phonon dispersions, phonon DOS, and thermodynamic properties using the command line.
This feature requires gnuplot version 5 and pdflatex.

To plot phonon dispersions, use:

  aflow --plotphdisp|--plotphonondispersion|--pphdis[=directory,[Emin,[Emax]]] [--keep=gpl] [--print=pdf|eps|gif|jpg|png] [--title=] [--unit=THz|Hz|eV|meV|rcm|cm-1] [--outfile=]
    Defaults.
     Directory: default = ./
     Emin: The minimum value of the plot. Default value: minimum frequency
     Emax: The maximum value of the plot. Default value: maximum frequency
    Optional arguments:
     --keep=gpl: keep the gnuplot script
     --print: specify the format of the output image (default: pdf)
     --title: the title of the plot in LaTeX format (default: title in PHEIGENVAL file)
     --unit: the unit of the phonon frequencies/energies (default: THz)
     --outfile: the name of the output image file
    Example:
     aflow --plotphdisp=./ --print=png --unit=cm-1

The files PHEIGENVAL, PHKPOINTS, and PHPOSCAR must reside in the directory.

Phonon DOS can be plotted in a similar fashion:

  aflow --plotphdos[=directory,[Emin,[Emax[,DOSSCALE]]]] [--keep=gpl] [--print=pdf|eps|gif|jpg|png] [--title=] [--unit=THz|Hz|eV|meV|rcm|cm-1] [--outfile=]
    Defaults.
     DOSSCALE: The scaling factor of the DOS plot. Default value: no scaling
    Other options: same as plotphdisp

The file PHDOSCAR must be present in the directory. To plot projected DOS, --projection=atom needs
to be added to the command and PHPOSCAR needs to be in the directory.

A combined dispersion-DOS plot can be created using:

  aflow --plotphdispdos[=directory,[Emin,[Emax[,DOSSCALE]]]] [--keep=gpl] [--print=pdf|eps|gif|jpg|png] [--title=] [--unit=THz|Hz|eV|meV|rcm|cm-1] [--outfile=]

Thermodynamic properties can be plotted using:

  aflow --plotthermo[=directory[,Tmin[,Tmax]]] [--keep=gpl] [--print=pdf|eps|gif|jpg|png] [--title=]
    Defaults.
     Directory: default = ./
     Tmin: The minimum temperature of the plot. Default value: 0 K
     Tmax: The maximum temperature of the plot. Default value: maximum temperature
    Optional arguments:
     --keep=gpl: keep the gnuplot script
     --print: specify the format of the output image (default: pdf)
     --title: the title of the plot in LaTeX format (default: SYSTEM in out file)
     --outfile: the name of the output image file
    Example:
     aflow --plotthermo=./,100,300

This command outputs four separate files for the internal energy, vibrational free energy,
vibrational entropy, and isochoral heat capacity. The file aflow.apl.thermal_properties.out must be
present in the directory.


3.5. Atomic Displacements

3.5.1. Mean Square Displacements

AFLOW can calculate the mean square displacement vectors for each atom in the unit cell by adding
the following to the aflow.in file:

    [AFLOW_APL]DISPLACEMENTS=ON

This is also the default option.

The temperature range is taken from the parameter TPT (see section 3.3.) and the mesh density from
DOS_MESH (see section 3.2.).

3.5.2. Visualizing Phonon Displacements

APL provides an interface to the following output formats that can be used to visualize atomic
displacements:

  * XCRYSDEN
  * V_sim

The commands to create these output files are shown below:

  aflow --visualize_phonons --qpoint=|--qpoints=|--q= [--format=XCRYSDEN|V_SIM] [--amplitude=] [--branches=] [--periods=] [--steps=] [--scell=] -D DIRECTORY
        Outputs a phonon visualization file. Currently supported formats: xcrysden and v_sim.

        Options for all formats:
          --qpoint:    q-points for the phonons to be visualized as a comma-separated list.
                       Multiple points can be specified.

          --format:    The output file format.
                       Default: DEFAULT_APL_ADISP_SCENE_FORMAT in .aflow.rc.

          --amplitude: The amplitude of the displacements.
                       Default: DEFAULT_APL_ADISP_AMPLITUDE in .aflow.rc.


        Options for XCRYSDEN:
          --branches:  The branches of the phonon modes as a comma-separated list (1-based).
                       Default: all.

          --periods:   The number of periods displayed in the output file.
                       Default: DEFAULT_APL_ADISP_NPERIODS in .aflow.rc

          --steps:     The number of steps per period in the output file.
                       Default: DEFAULT_APL_ADISP_NSTEPS in .aflow.rc

          --scell:     The supercell to be displayed in the output file.
                       Default: unit cell.

        Examples:
          aflow --visualize_phonons --qpoint=0.5,0,0 --steps=40 --amplitude=0.3 --scell=3x3x3 -D ./
          aflow --visualize_phonons --q=0,0,0,0.5,0.5,0.5 --branches=4,5,6 --format=XCRYSDEN --scell=4x4x2
          aflow --visualize_phonons --q=0,0,0 --branches=1,2,3 --format=V_SIM --amplitude=0.2

For XCRYSDEN, APL creates an animated structure (axsf) file for each phonon mode. These files
contain animations that show the motions of the atoms for each of the modes. These animations can
be used by XCRYSDEN to create mpeg movies or gif images. The parameter "steps" determines the
number of steps for each period, and "periods" determines how many periods are output into the
file. For gifs, one period is sufficient. The parameter "scell" can be used to show the
displacements in a supercell instead of the unit cell, which is important for modes off the Gamma
point.

V_sim files will be output as aflow.apl.displacements.ascii, which contains all modes for the
desired q-points. The file can be opened using V_sim or ASCII-phonons. These programs are capable
of doing supercells and animations on their own and do not require steps, periods, or supercell
parameters.

3.6. Group Velocities

AFLOW can calculate group velocities by adding the following to the aflow.in file:

    [AFLOW_APL]GROUP_VELOCITY=ON

This is also the default option.

The q-points for the group velocities are taken from the DOS_MESH keyword (see section 3.2.).


4. Saving and Loading APL Results APL Calculations (Hibernation)

In order to avoid recalculating the force constant matrix with every run, APL has implemented a
hibernation system which saves the force constants into the file aflow.apl.harmonicIFCs.xml.
Born effective charges and the dielectric tensor are stored in aflow.apl.polar.xml.

This behavior can be switched on(off) using the HIBERNATE keyword: 

    [AFLOW_APL]HIBERNATE=ON(OFF)

It is turned on by default. If settings in the aflow.in file change the calculation of the force
constants, these files must be deleted for the changes to take effect.


5. Inheritance from AFLOW

The generated ARUN.APL***/aflow.in input files inherit some options from the master aflow.in
If present, the following options are passed to the local aflow.in:

    [VASP_FORCE_OPTION]AUTO_PSEUDOPOTENTIALS=...
    [VASP_FORCE_OPTION]SPIN=...
    [VASP_FORCE_OPTION]TYPE=...
    [VASP_FORCE_OPTION]AUTO_MAGMOM=...
    [VASP_FORCE_OPTION]BADER=...
    [VASP_FORCE_OPTION]ELF=...
    [VASP_FORCE_OPTION]ABMIX=...
    [VASP_FORCE_OPTION]IGNORE_AFIX=...

Other options are added by default:

    [VASP_FORCE_OPTION]WAVECAR=OFF
    [VASP_FORCE_OPTION]CHGCAR=OFF
    [VASP_FORCE_OPTION]ALGO=NORMAL

To avoid NPAR warnings from VASP 5.2+, the following option is added:

    [VASP_FORCE_OPTION]IGNORE_AFIX=NPARC

6.  Example Input File

The following is an example APL input file for PbS using the linear response method. For the
direct method, use [AFLOW_APL]ENGINE=DM and uncomment #[AFLOW_APL]DMAG=0.015.

[AFLOW] **********************************************************************************
[AFLOW]SYSTEM=PbS
[AFLOW] **********************************************************************************
[AFLOW_MODE=VASP]
[AFLOW_MODE_ZIP=xz]
[AFLOW_MODE_BINARY=vasp54s]
[AFLOW] **********************************************************************************
[AFLOW_MODE_MPI]
[AFLOW_MODE_MPI_MODE]NCPUS=MAX
[AFLOW_MODE_MPI_MODE]COMMAND ="mpirun -np"
[AFLOW_MODE_MPI_MODE]AUTOTUNE
[AFLOW_MODE_MPI_MODE]BINARY="mpivasp54s"
[AFLOW] **********************************************************************************
[AFLOW_APL]CALC
[AFLOW_APL]SUPERCELL=4x4x4
[AFLOW_APL]DC=ON
[AFLOW_APL]DCPATH=LATTICE
[AFLOW_APL]DCPOINTS=100
[AFLOW_APL]DINEQUIV_ONLY=ON
[AFLOW_APL]DISPLACEMENTS=ON
[AFLOW_APL]DMAG=0.015
[AFLOW_APL]DOS=ON
[AFLOW_APL]DOSMESH=21x21x21
[AFLOW_APL]DOSMETHOD=LT
[AFLOW_APL]DOSPOINTS=2000
[AFLOW_APL]DOS_PROJECT=OFF
[AFLOW_APL]DPM=ON
[AFLOW_APL]DSYMMETRIZE=ON
[AFLOW_APL]DXYZONLY=OFF
[AFLOW_APL]ENGINE=DM
[AFLOW_APL]FREQFORMAT=THZ|ALLOW_NEGATIVE
[AFLOW_APL]GROUP_VELOCITY=ON
[AFLOW_APL]HIBERNATE=ON
[AFLOW_APL]KPPRA=3000
[AFLOW_APL]KSCHEME=G
[AFLOW_APL]MINATOMS=175
[AFLOW_APL]POLAR=OFF
[AFLOW_APL]PREC=PHONONS
[AFLOW_APL]RELAX=ON
[AFLOW_APL]RELAX_COMMENSURATE=ON
[AFLOW_APL]TP=ON
[AFLOW_APL]TPT=0:2000:10
[AFLOW_APL]ZEROSTATE=OFF
[AFLOW] **********************************************************************************
[VASP_RUN]RELAX=2
#[VASP_FORCE_OPTION]KPOINTS=KEEPK,EVEN,KSHIFT_GAMMA_EVEN,KSCHEME_GAMMA,GAMMA,IBZKPT
[VASP_FORCE_OPTION]SYM=ON
[VASP_FORCE_OPTION]AUTO_PSEUDOPOTENTIALS=potpaw_PBE
[VASP_FORCE_OPTION]NBANDS
[VASP_FORCE_OPTION]SPIN=OFF
[VASP_FORCE_OPTION]RELAX_MODE=FORCES
[VASP_FORCE_OPTION]PREC=ACCURATE
[VASP_FORCE_OPTION]ALGO=FAST
[VASP_FORCE_OPTION]RELAX
[VASP_FORCE_OPTION]TYPE=DEFAULT
[VASP_FORCE_OPTION]CONVERT_UNIT_CELL=SPRIM,MINK
[AFLOW] **********************************************************************************
[VASP_INCAR_MODE_EXPLICIT]START
SYSTEM=PbS
[VASP_INCAR_MODE_EXPLICIT]STOP
[AFLOW] **********************************************************************************
[VASP_KPOINTS_MODE_IMPLICIT]
[VASP_KPOINTS_FILE]KSCHEME=G
[VASP_KPOINTS_FILE]KPPRA=400
[AFLOW] **********************************************************************************
[VASP_POTCAR_MODE_IMPLICIT]
[VASP_POTCAR_FILE]Pb
[VASP_POTCAR_FILE]S
[AFLOW] **********************************************************************************
[VASP_POSCAR_MODE_EXPLICIT]START
PbS   [FCC,FCC,cF8]
1.224745
   0.00000000000000   2.45291416610433   2.45291416610433
   2.45291416610433   0.00000000000000   2.45291416610433
   2.45291416610433   2.45291416610433   0.00000000000000
1 1
Direct(2) [A1B1]
   0.00000000000000   0.00000000000000   0.00000000000000  Pb
   0.50000000000000   0.50000000000000   0.50000000000000  S
[VASP_POSCAR_MODE_EXPLICIT]STOP
[AFLOW] **********************************************************************************

***************************************************************************************************
NB:
Any shared tags between APL/QHA/AAPL can be equivalently specified with [AFLOW_APL], [AFLOW_QHA],
[AFLOW_AAPL] in the aflow.in. This allows users the flexibility to first run an APL calculation and
then proceed with a subsequent analysis (QHA/AAPL) without the need to change the tags.
***************************************************************************************************


7. The Quasi-Harmonic Approximation

7.1. Running QHA Calculations

The Quasi-Harmonic Approximation (QHA) calculates harmonic phonons at different volumes to compute
the volume-dependent free energies. This data is then fitted to an equation of state (EOS) to
obtain a temperature-dependent free-energy curve. This curve is the basis for a variety of
thermomechanical properties such as thermal expansion, the bulk modulus, and heat capacities.

To use QHA, the following needs to be added the aflow.in file:

    [AFLOW_QHA]CALC

This overrides all other APL CALC options. A QHA calculation requires three steps:

    1. In a first AFLOW run, QHA will create subdirectories for an APL run and a STATIC run for each
       selected volume. If the APL option REALX is set to ON, QHA will relax the structure before
       creating the directories.
    2. Run all DFT calculations required for each APL and STATIC directory.
    3. Post-processing: run AFLOW again in the QHA directory to calculate thermomechanical
       properties.

QHA uses many parameters that are also used inside APL. See section 2. for a discussion on APL
options.


7.2. QHA Modes

QHA can either calculate only Grueneisen parameters or the full equation of state. To calculate
Grueneisen parameters, use:

    [AFLOW_QHA]GP_FINITE_DIFF=ON|OFF
    [AFLOW_QHA]GP_DISTORTION=3

GP_DISTORTION is the volume expansion and compression in %. This will create three APL directories,
at 97%, 100%, and 103% volume. No static calculations are needed.

To run equation of state calculation, the following needs to be present in the aflow.in file:

    [AFLOW_QHA]EOS=ON
    [AFLOW_QHA]EOS_DISTORTION_RANGE=VSTART:VEND:VSTEP

Where VSTART and VEND are the initial and final volumes and VSTEP is the step size in %. VSTART must
be negative. These volumes will be used for the static calculations. The number of phonon
calculations depend on the MODE parameter:

    [AFLOW_QHA]MODE=QHA|SCQHA|QHA3P|QHANP

QHA runs standard QHA and uses EOS_DISTORTION_RANGE. SCQHA (self-consistent QHA) and QHA3P (3-phonon
QHA) only need three volumes and use GP_DISTORTION. QHANP (N-th order QHA) calculates the
vibrational free energy based on a Taylor expansion of order N:

    [AFLOW_QHA]TAYLOR_EXPANSION_ORDER=N

It generates (2N + 1) phonon calculations.


7.3. Equation of State Calculations

When EOS=ON, QHA can be used to calculate thermomechanical properties. For this, the obtained free
energies need to be fit to an equation of state:

    [AFLOW_QHA]EOS_MODEL=M|BM2|BM3|BM4|SJ

Here, M is the Murnaghan, BMX are the X-th order Birch-Murnaghan, and SJ the Stabilized Jellium
model. The default is SJ. QHA will use the [AFLOW_APL]TPT option for the temperature range over
which to fit.

Electronic contributions to the free energy can be included via:

    [AFLOW_QHA]INCLUDE_ELEC_CONTRIB=ON|OFF
    [AFLOW_QHA]SOMMERFELD_EXPANSION=ON|OFF

Setting SOMMERFELD_EXPANSION to ON will have QHA use the Sommerfeld model to calculate electronic
free energies. If OFF, it will integrate the electronic structure.

It is possible that at some volumes, the structure becomes dynamically unstable and has imaginary
frequencies. QHA will not process volume sets that contain such structures and the set of volumes
needs to be adjusted. Imaginary frequencies can also be ignored using:

    [AFLOW_QHA]IGNORE_IMAGINARY=ON|OFF

Temperature-dependent phonon dispersion can also be calculated. The temperatures are set using:

    [AFLOW_QHA]T1,T2,T3,...

Where T1, etc. are temperatures in K. This is not supported for QHANP.


7.4. Plotting Thermomechanical Properties

Thermophysical properties can be plotted using:

  aflow --plotthermoqha[=directory[,Tmin[,Tmax]]] --eosmodel=SJ|BM2|BM3|BM4|M [--keep=gpl] [--print=pdf|eps|gif|jpg|png] [--title=]
    Defaults.
     Directory: default = ./
     Tmin: The minimum temperature of the plot. Default value: 0 K
     Tmax: The maximum temperature of the plot. Default value: maximum temperature
    --eosmodel: The equation of state. Default: SJ
    Optional arguments:
     --keep=gpl: keep the gnuplot script
     --print: specify the format of the output image (default: pdf)
     --title: the title of the plot in LaTeX format (default: SYSTEM in out file)
     --outfile: the name of the output image file
    Example:
     aflow --plotthermo=./,100,300

This command outputs separate files for the volume, free energy, bulk modulus, volumetric thermal
expansion coefficient, isochoric and isobaric heat capacity, and average Grueneisen parameter. The
file aflow.qha.thermo.out must be present in the directory and must contain data for the selected
equation of state.

The dispersion of the Grueneisen parameter can be plotted via:

  aflow --plotgrdisp|--plotgruneisendispersion|--plotgrueneisendispersion[=directory,[gmin,[gmax]]] [--keep=gpl] [--print=pdf|eps|gif|jpg|png] [--title=] [--unit=THz|Hz|eV|meV|rcm|cm-1] [--outfile=]
    Defaults.
     Directory: default = ./
     gmin: The minimum value of the plot. Default value: minimum Grueneisen parameter
     gmax: The maximum value of the plot. Default value: maximum Grueneisen parameter
    Optional arguments:
     --keep=gpl: keep the gnuplot script
     --print: specify the format of the output image (default: pdf)
     --title: the title of the plot in LaTeX format (default: title in gp.disp.out file)
     --unit: the unit of the phonon frequencies/energies (default: THz)
     --outfile: the name of the output image file
    Example:
     aflow --plotgrdisp=./ --print=png

The files PHPOSCAR, aflow.qha.kpoints.out, and aflow.qha.qp.disp.out must be present in the
directory.


7.5. Example QHA Input File

[AFLOW] *******************************************************************************************
[AFLOW]SYSTEM=MgO
[AFLOW] *******************************************************************************************
[AFLOW_MODE_ZIP=xz]
[AFLOW_MODE_BINARY=vasp46s]
[AFLOW_MODE_MPI]
[AFLOW_MODE_MPI_MODE]NCPUS=MAX
[AFLOW_MODE_MPI_MODE]COMMAND ="mpirun -np"
[AFLOW_MODE_MPI_MODE]AUTOTUNE
[AFLOW_MODE_MPI_MODE]BINARY="mpivasp46s"
[AFLOW] *******************************************************************************************
[AFLOW_APL]SUPERCELL=4x4x4
[AFLOW_APL]DC=ON
[AFLOW_APL]DCPATH=LATTICE
[AFLOW_APL]DCPOINTS=100
[AFLOW_APL]DINEQUIV_ONLY=ON
[AFLOW_APL]DISPLACEMENTS=ON
[AFLOW_APL]DMAG=0.015
[AFLOW_APL]DOS=ON
[AFLOW_APL]DOSMESH=21x21x21
[AFLOW_APL]DOSMETHOD=LT
[AFLOW_APL]DOSPOINTS=2000
[AFLOW_APL]DOS_PROJECT=OFF
[AFLOW_APL]DPM=ON
[AFLOW_APL]DSYMMETRIZE=ON
[AFLOW_APL]DXYZONLY=OFF
[AFLOW_APL]ENGINE=DM
[AFLOW_APL]FREQFORMAT=THZ|ALLOW_NEGATIVE
[AFLOW_APL]GROUP_VELOCITY=ON
[AFLOW_APL]HIBERNATE=ON
[AFLOW_APL]KPPRA=3000
[AFLOW_APL]KSCHEME=G
[AFLOW_APL]MINATOMS=175
[AFLOW_APL]POLAR=OFF
[AFLOW_APL]PREC=PHONONS
[AFLOW_APL]RELAX=ON
[AFLOW_APL]RELAX_COMMENSURATE=ON
[AFLOW_APL]TP=ON
[AFLOW_APL]TPT=0:2000:10
[AFLOW_APL]ZEROSTATE=OFF
[AFLOW] *******************************************************************************************
[AFLOW_QHA]CALC
[AFLOW_QHA]EOS=ON
[AFLOW_QHA]EOS_DISTORTION_RANGE=-3:9:1
[AFLOW_QHA]EOS_MODEL=SJ
[AFLOW_QHA]GP_DISTORTION=1
[AFLOW_QHA]GP_FINITE_DIFF=OFF
[AFLOW_QHA]IGNORE_IMAGINARY=OFF
[AFLOW_QHA]INCLUDE_ELEC_CONTRIB=OFF
[AFLOW_QHA]MODE=QHA
[AFLOW_QHA]PDIS_T=300
[AFLOW_QHA]SOMMERFELD_EXPANSION=OFF
[AFLOW_QHA]TAYLOR_EXPANSION_ORDER=2
[AFLOW] *******************************************************************************************
[VASP_FORCE_OPTION]SYM=ON
[VASP_FORCE_OPTION]AUTO_PSEUDOPOTENTIALS=potpaw_PBE
[VASP_FORCE_OPTION]NBANDS
[VASP_FORCE_OPTION]SPIN=OFF
[VASP_FORCE_OPTION]TYPE=DEFAULT
[VASP_FORCE_OPTION]PREC=ACCURATE
[VASP_FORCE_OPTION]ALGO=NORMAL
[AFLOW] *******************************************************************************************
[VASP_POTCAR_MODE_IMPLICIT]
[VASP_POTCAR_FILE]Mg
[VASP_POTCAR_FILE]O
[AFLOW] *******************************************************************************************
[VASP_KPOINTS_MODE_IMPLICIT]
[VASP_KPOINTS_FILE]STATIC_KSCHEME=M
[VASP_KPOINTS_FILE]STATIC_KPPRA=10000
[VASP_KPOINTS_FILE]BANDS_LATTICE=AUTO
[VASP_KPOINTS_FILE]BANDS_GRID=20
[AFLOW] *******************************************************************************************
[VASP_POSCAR_MODE_EXPLICIT]START
MgO
1.224745
   0.00000000000000   1.73138099985724   1.73138099985724
   1.73138099985724   0.00000000000000   1.73138099985724
   1.73138099985724   1.73138099985724   0.00000000000000
Mg O
1 1
Direct(2) [A1B1]
   0.00000000000000   0.00000000000000   0.00000000000000  Mg
   0.50000000000000   0.50000000000000   0.50000000000000  O
[VASP_POSCAR_MODE_EXPLICIT]STOP
[AFLOW] *******************************************************************************************


8. The Automatic Anharmonic Phonon Library

8.1. Running AAPL Calculations

The Automatic Anharmonic Phonon Library (AAPL) is the anharmonic extension to APL and calculates
anharmonic interatomic force constants (IFCs). These can then be used to calculate the lattice
thermal conductivity by solving the Boltzmann Transport Equation (BTE). AAPL can use three-phonon
processes to calculate thermal conductivities, Grueneisen parameters, and the scattering phase
space of a material.

To use AAPL, the following needs to be added to the aflow.in file:

    [AFLOW_AAPL]CALC

This overrides all other APL CALC options. An AAPL calculation requires three separate AFLOW runs:

    1. In the first run, the distortions for the anharmonic IFCs (and harmonic IFCs, if not
       already present) are created. If the option REALX is set to ON, APL will relax the structure
       before creating the distortions.

    2. The DFT calculations are performed. However, due to the large number of DFT calculations, it
       is recommended to run each calculation in its own separate AFLOW run.

    3. Post-processing: the harmonic and anharmonic IFCs and the thermal conductivity tensor are
       calculated.

AAPL shares the supercell parameters and the HIBERNATE keyword with APL. AAPL hibernate file store
the clusters and the anharmonic IFCs. IT IS STRONGLY RECOMMENDED TO USE THE HIBERNATE FUNCTION! The
number of calculations for AAPL is very large, so processing the DFT calculations can take a long
time whereas reading the hibernate files takes only a couple of seconds. Moreover, if LR is chosen
for the APL calculations, the DMAG keyword still needs to be specified for AAPL.

8.2. Anharmonic Interatomic Force Constants

Anharmonic IFCs are calculated using the direct method. It is necessary to define a cut-off radius
beyond which the interactions between atoms are considered negligible. Within this shell, clusters
of atoms (three atoms for third order IFCs, four atoms for fourth order IFCs, etc.) are created and
distorted. The number of clusters and distortions is then reduced by symmetry to save computational
resources. The resulting distorted structures are then output into their own aflow.in files inside
a subfolder. These folder all have the name:

    ARUN.AAPL_N_XX_Aat1Dd1-Aat2Dd2-...

Where the N in N is replaced with the phonon process order; at1, at2, etc. are the atomic indices
of the distorted atoms; and d1, d2, etc. are the distortion indices corresponding to the atom. They
range from 0 to 5 and correspond to Cartesian +x, +y, +z, -x, -y, -z. Note that when one distortion
index appears twice on the same atom, the atom is not distorted twice into the same direction.

Once the DFT calculations are finished, the forces are read and the IFCs are determined using a
self-consistent field. 

8.2.1. Cluster Parameters

The supercell and the distortion magnitude are defined the same way as in APL. See Section 2.1.1.
for the supercell parameters and Section 2.2. for the distortion magnitude.

The cut-off radius described above can be defined in two ways:

  * Define a cut-off radius directly:

        [AFLOW_AAPL]CUT_RAD=5.0

  * Define the number of coordination shells that should be included in the interaction volume.
    This can be achieved via the CUT_SHELL parameter.

        [AFLOW_AAPL]CUT_SHELL=6

    This is also the default value. The cut-off radius is then determined automatically.

It is possible to use both parameters. In this case, AAPL determines the cut-off radius based on
the CUT_SHELL parameter and chooses whichever value is greater.

8.2.2. Sum Rule Parameters

The iterative procedure behind the calculation of the anharmonic IFCs has a few adjustable
parameters.

    * The most important parameter is the convergence criterion for the sum rule, which can
      be set using the SUMRULE parameter, e.g.:

          [AFLOW_AAPL]SUMRULE=1E-7

      This is also the default value.

    * The maximum number of iterations can be set using the parameter SUMRULE_MAX_ITER:

          [AFLOW_AAPL]SUMRULE_MAX_ITER=2000

      This is also the default value, which should be more than enough.

    * The last adjustable parameter is the mixing coefficient, which controls how much of the
      previous iteration step is mixed into the new IFC values:

          [AFLOW_AAPL]MIXING_COEFFICIENT=0.0

      It defaults to 0.0, i.e. there is no mixing between two iteration steps. However, if
      convergence cannot be reached because the corrections are too large, this parameter should
      be adjusted.

8.2.3. k-Point Parameters

The KPPRA parameter can be set separately for AAPL calculations. This may be desirable to reduce
the computational cost. The parameter to set is KPPRA_AAPL:

    [AFLOW_AAPL]KPPRA_AAPL=1000

If no entry is found or if the parameter is smaller than 1, the KPPRA value for APL will be used.

8.3. Lattice Thermal Conductivity Calculations

The lattice thermal conductivity tensor is calculated by solving the Boltzmann Transport
Equation (BTE) for a specific temperature range. To set the temperature range, the parameter TCT
needs to be set similar to TPT in APL (see Section 3.3.):

    [AFLOW_AAPL]TCT=50:550:10

The default setting is 50:550:50. The starting temperature cannot be set to 0 K because the thermal
conductivity is infinite. If it is set to 0 K, AAPL will simply skip it.

The BTE can be solved using either the Relaxation Time Approximation (RTA) or an iterative
procedure, which is done by the BTE keyword:

    [AFLOW_AAPL]BTE=RTA|FULL

The default is BTE=FULL.

The BTE gets solved along a q-point grid which can be defined using THERMALGRID:

    [AFLOW_AAPL]THERMALGRID=q1xq2xq3

Where q1, q2, and q3 are the number of q-points along each dimension. The default is 21x21x21.

There are various factors which can influence the phonon scattering rates besides anharmonic phonon
processes. Two mechanisms that are implemented in AAPL: isotope scattering and grain boundary
scattering.

    * Isotope scattering can be switched on by setting the parameter ISOTOPE:

          [AFLOW_AAPL]ISOTOPE=ON

      The default is ISOTOPE=OFF.

    * Grain boundary scattering can be included in two ways:

        1. By introducing a grain boundary scattering term into the scattering rate. This is done
           by setting BOUNDARY to true:

              [AFLOW_AAPL]BOUNDARY=ON

           The default is BOUNDARY=OFF.

        2. By only including phonon modes with a mean free path below a threshold into the thermal
           conductivity calculation. This is done using the keyword CUMULATIVEK:

               [AFLOW_AAPL]CUMULATIVEK=ON

           The default is CUMULATIVEK=OFF.

        If both are set to true, BOUNDARY will be used, i.e. CUMULATIVEK will be set to false.

        Both methods require an additional parameter, NANO_SIZE:

            [AFLOW_AAPL]NANO_SIZE=300

        Which corresponds to the grain size/mean free path threshold in nm. The default is 100 nm.

8.4. Plotting Thermal Conductivity

Thermal conductivity can be plotted using the command line:

  aflow --plottcond[=directory[,Tmin[,Tmax]]] [--keep=gpl] [--print=pdf|eps|gif|jpg|png] [--title=] [--outfile=]

  Defaults.
   Directory: default = ./
   Tmin: The minimum temperature of the plot. Default value: 0 K
   Tmax: The maximum temperature of the plot. Default value: maximum temperature
  Optional arguments:
   --keep=gpl: keep the gnuplot script
   --print: specify the format of the output image (default: pdf)
   --title: the title of the plot in LaTeX format (default: SYSTEM in out file)
   --outfile: the name of the output image file
  Example:
   aflow --plottcond=./,100,300

aflow.aapl.tcond.out must be present in the directory.

This feature requires gnuplot version 5 and pdflatex.

8.5. Example AAPL Input File

The following is an example AAPL input file for NaF. Four-phonon processes are included.

[AFLOW] *******************************************************************************************
[AFLOW]System=NaF
[AFLOW] *******************************************************************************************
[AFLOW_MODE=VASP]
[AFLOW_MODE_ZIP=bzip2]
[AFLOW_MODE_BINARY=vasp46s]
[AFLOW_MODE_MPI]
[AFLOW_MODE_MPI_MODE]NCPUS=MAX
[AFLOW_MODE_MPI_MODE]COMMAND="mpirun -np"
[AFLOW_MODE_MPI_MODE]AUTOTUNE
[AFLOW_MODE_MPI_MODE]BINARY="mpivasp46s"
[AFLOW] *******************************************************************************************
[AFLOW_AAPL]CALC
[AFLOW_AAPL]CUT_SHELL=4
#[AFLOW_AAPL]CUT_RAD=5.0
[AFLOW_AAPL]SUMRULE=1E-7
[AFLOW_AAPL]SUMRULE_MAX_ITER=1E-7
[AFLOW_AAPL]MIXING_COEFFICIENT=1.0
[AFLOW_AAPL]THERMALGRID=26x26x26
[AFLOW_AAPL]ISOTOPE=ON
[AFLOW_AAPL]BTE=FULL
[AFLOW_AAPL]TCT=50:550:10
#[AFLOW_AAPL]CUMULATIVEK=OFF
#[AFLOW_AAPL]BOUNDARY=OFF
#[AFLOW_AAPL]NANO_SIZE=100
[AFLOW] *******************************************************************************************
[AFLOW_APL]ENGINE=DM
[AFLOW_APL]RELAX=ON
[AFLOW_APL]PREC=PHONONS
[AFLOW_APL]DMAG=0.015
[AFLOW_APL]POLAR=ON
[AFLOW_APL]MINATOMS=100
[AFLOW_APL]HIBERNATE=ON
[AFLOW_APL]KPPRA=8000
[AFLOW_APL]KSCHEME=G
[AFLOW_APL]KPOINTS=EVEN
[AFLOW] *******************************************************************************************
[VASP_FORCE_OPTION]SYM=ON
[VASP_FORCE_OPTION]AUTO_PSEUDOPOTENTIALS=potpaw_PBE
[VASP_FORCE_OPTION]NBANDS
[VASP_FORCE_OPTION]SPIN=OFF
[VASP_FORCE_OPTION]TYPE=DEFAULT
[VASP_FORCE_OPTION]PREC=ACCURATE
[VASP_FORCE_OPTION]ALGO=NORMAL
[AFLOW] *******************************************************************************************
[VASP_POTCAR_MODE_IMPLICIT]
[VASP_POTCAR_FILE]F
[VASP_POTCAR_FILE]Na
[AFLOW] *******************************************************************************************
[AFLOW_SYMMETRY]CALC
[VASP_POSCAR_MODE_EXPLICIT]START
NaF
1.224745
   0.00000000000000   1.91118771995223   1.91118771995223
   1.91118771995223   0.00000000000000   1.91118771995223
   1.91118771995223   1.91118771995223   0.00000000000000
F Na 
1 1 
Direct(2) [A1B1] 
   0.00000000000000   0.00000000000000   0.00000000000000  F     
   0.50000000000000   0.50000000000000   0.50000000000000  Na    
[VASP_POSCAR_MODE_EXPLICIT]STOP
[AFLOW] *******************************************************************************************


9. References

APL:

M. Esters, C. Oses, D. Hicks, M. J. Mehl, M. Jahnatek, M. D. Hossain, J.-P. Maria, D. W. Brenner,
C. Toher, and S. Curtarolo, Settling the matter of the role of vibrations in the stability of
high-entropy carbides, Nat. Commun. 12, 5747 (2021). DOI: 10.1038/s41467-021-25979-5.

QHA:

P. Nath, J. J. Plata, D. Usanmaz, R. Al Rahal Al Orabi, M. Fornari, M. Buongiorno Nardelli,
C. Toher, and S. Curtarolo, High-throughput prediction of finite-temperature properties using the
quasi-harmonic approximation, Comput. Mater. Sci. 125, 8291 (2016).
DOI: 10.1016/j.commatsci.2016.07.043

P. Nath, D. Usanmaz, D. Hicks, C. Oses, M. Fornari, M. B. Nardelli, C. Toher, and S. Curtarolo,
AFLOWi-QHA3P: Robust and automated method to compute thermodynamic properties of solids, Phys. Rev.
Materials 3, 073801 (2019). DOI: 10.1103/PhysRevMaterials.3.073801

AAPL:

J. J. Plata, P. Nath, D. Usanmaz, J. Carrete, C. Toher, M. de Jong, M. D. Asta, M. Fornari,
M. Buongiorno Nardelli, and S. Curtarolo, An ecient and accurate framework for calculating lattice
thermal conductivity of solids: AFLOW-AAPL Automatic Anharmonic Phonon Library, npj Comput. Mater.
3, 45 (2017). DOI: 10.1038/s41524-017-0046-7

****************************************************************************************************
*                                                                                                  *
*                        aflow - STEFANO CURTAROLO Duke University 2003-2021                       *
*                            High-Throughput ab-initio Computing Project                           *
*                                                                                                  *
****************************************************************************************************
