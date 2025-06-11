****************************************************************************************************
*                                                                                                  *
*                        aflow - STEFANO CURTAROLO Duke University 2003-2021                       *
*                            High-Throughput ab-initio Computing Project                           *
*                                                                                                  *
****************************************************************************************************
LATEST VERSION OF THE FILE:                       materials.duke.edu/AFLOW/README_AFLOW_SYM.TXT
****************************************************************************************************

AFLOW-SYM README:

README written by: David Hicks (david.hicks@duke.edu)

Citation info: D. Hicks, C. Oses, E. Gossett, G. Gomez, R.H. Taylor
               C. Toher, M.J. Mehl, O. Levy, and S. Curtarolo,
               AFLOW-SYM: platform for the complete, automatic and self-consistent
               symmetry analysis of crystals, Acta Cryst. A74, 184-203 (2018).


****************************************************************************************************
GENERAL OVERVIEW

AFLOW-SYM is a robust symmetry analysis tool catered to - but not limited to -
high-throughput frameworks.  With built-in self-consistency symmetry validation and
an adaptive tolerance scheme, AFLOW-SYM resolves a globally consistent symmetry
profile.

The AFLOW-SYM functions fully characterize the symmetry of a crystal, including
the space group, lattice types, Pearson symbol, Wyckoff positions, point group,
factor group, symmetrically equivalent atoms, etc.  The AFLOW-SYM commands, along
with important usage information, are listed below.

README Sections:
  - LIST OF AFLOW-SYM COMMANDS
  - SUPPORTED GEOMETRY FILE FORMATS
  - SELF-CONSISTENT SYMMETRY
  - TOLERANCE
  - TOLERANCE SPECTRUM
  - MAGNETIC ANALYSIS
  - CELL/SETTING CHOICES
  - OUTPUT FORMAT
  - AFLOW COMMANDS USAGE


****************************************************************************************************
LIST OF AFLOW-SYM COMMANDS:

  aflow --aflow-sym               : Returns symmetry operations of the point groups (lattice, crystal,
                                    reciprocal lattice, and dual of crystal), factor group, space group,
                                    and site symmetry; along with symmetrically equivalent/unique atoms.
  aflow --aflowSG|--space_group   : Returns space group in the ITC convention
  aflow --cif                     : Returns structure in CIF file format
  aflow --edata                   : Returns extended crystallographic symmetry data (edata) of the
                                    crystal, lattice, reciprocal lattice, and superlattice
  aflow --equivalent              : Returns symmetrically equivalent/unique atoms
  aflow --factorgroup             : Returns the coset representatives of the factor group (unit cell
                                    symmetry operations)
  aflow --lattice_data            : Returns the lattice type data (lattice, type/variation, Pearson, etc.)
  aflow --lattice_type            : Returns the lattice type corresponding to the crystal (lattice + atoms) symmetry
  aflow --lattice_lattice_type    : Returns lattice type corresponding to the lattice symmetry
  aflow --pearson_symbol          : Returns Pearson symbol
  aflow --pointgroup              : Returns point group operations of the lattice
  aflow --pointgroup_crystal      : Returns point group operations of the crystal
  aflow --pointgroup_crystal_data : Returns the point group descriptors (e.g., family, class, type, etc.)
  aflow --pointgroupkcrystal      : Returns dual of crystal point group operations
  aflow --pointgroupklattice      : Returns point group operations of the reciprocal lattice
  aflow --pointgroupk_Patterson   : Returns point group operations of the Patterson symmetry (symmorphic/centrosymmetric)
  aflow --reciprocal_lattice_data : Returns the reciprocal lattice type data (lattice, type/variation, Pearson, etc.)
  aflow --sc                      : Returns AFLOW standard conventional cell
  aflow --sgdata                  : Returns extended space group information (number, labels, classes,
                                    settings, Wyckoff positions, etc.)
  aflow --sitepointgroup          : Returns site symmetry operations for atoms in the unit cell
  aflow --sp                      : Returns AFLOW standard primitive cell
  aflow --spacegroup              : Returns space group operations out to a certain radius
  aflow --superlattice_data       : Returns the superlattice type data (lattice, type/variation, Pearson, etc.)
  aflow --wyccar                  : Returns Wyckoff positions in a POSCAR file format
  aflow --wyckoff                 : Returns Wyckoff positions (multiplicity, letter, site symmetry, position)

For details on usage, see the "AFLOW COMMANDS USAGE" section below.


****************************************************************************************************
SUPPORTED GEOMETRY FILE FORMATS:

AFLOW-SYM reads the following geometry file types:

  - Crystallographic Input File (CIF)
  - VASP (POSCAR or CONTCAR)
  - QuantumESPRESSO
  - ABINIT
  - FHI-AIMS
  - ELK

More file types will be added in later versions.


****************************************************************************************************
SELF-CONSISTENT SYMMETRY:

AFLOW-SYM includes self-consistency symmetry checks to ensure descriptions are commensurate
(e.g., valid point group or space group, the cardinality of the factor group is homomorphic to
the point group of the crystal, etc.).

There are some commands that validate the full symmetry profile (i.e., check descriptors
across different analyses). The following commands provide the most extensive self-consistency
checks and are recommended for classifying the symmetry of a material:

  aflow --aflow-sym : Validates symmetry operations between point groups, factor group,
                      space group, site symmetry, and equivalent/unique atoms
  aflow --edata     : Validates symmetry between the space group, lattice types, Pearson
                      symbol, etc.

The remaining commands are performed independently and are guaranteed to be self-consistent;
however, they are not validated against the full symmetry profile.

For example, calculating the space group (--aflowSG) and Pearson symbol (--pearson) separately
does not guarantee the two results will be commensurate with one another.
However, the space group and Pearson will be validated against one another when running the
extended crystallographic data command (--edata).

****************************************************************************************************
TOLERANCE:

AFLOW-SYM offers an adaptive tolerance scheme to combat issues arising from
numerical instabilities.  If any of the symmetry analyses yield incommensurate
symmetry, a radial tolerance scan is performed around the inital input tolerance.
All symmetry properties are then recalculated from the new tolerance.

The initial tolerance can be specified by appending it to the command with an "=".

  Ex: aflow --aflow-sym=<tolerance>

  Preset AFLOW tolerances:
    tight : minimum_interatomic_distance/100.0 (default)
    loose : minimum_interatomic_distance/10.0

    Ex: aflow --aflow-sym=tight  OR  aflow --aflow-sym=loose

  User defined tolerances:
    Any value between 1e-10 and the minimum interatomic distance

    Ex: aflow --aflow-sym=1e-2

    Note, the value must be below the minimum interatomic distance,
    otherwise an error will be thrown.

No tolerance scan option:
  The tolerance scan functionality can be suppresed by appending the option flag
  [--no_scan] to the command.

****************************************************************************************************
TOLERANCE SPECTRUM:

Users can specify a tolerance range (and number of steps) to see how symmetry
descriptors change as a function of tolerance (e.g., Figure 2 of the AFLOW-SYM article).

This can be done with the following --tolerance_spectrum option:

  Ex: aflow --aflowSG --tolerance_spectrum=0.001:0.1:20

where the following three fields are colon-separated:

  field 1 : starting tolerance (must be > 0)
  field 2 : stopping tolerance (must be <= minimum interatomic distance)
  field 3 : number of steps    (must be an integer)

The rountine automatically determines the step size and intermediate tolerances
based on these three inputs.

By default, the tolerance scan will be employed if incommensurate symmetry is found for any
of the tolerance values. To suppress this, use the --no_scan option. This will reveal any
"confusion" regions in the tolerance spectrum, where the symmetry could not be determined.

Relevant for the following commands:

  aflow --aflowSG
  aflow --pearson_symbol

****************************************************************************************************
MAGNETIC ANALYSIS:

The symmetry of magnetic systems is calculated by adding the magnetic moment
information (symmetry of lattice, atoms, and spin). The analysis works for both collinear
and non-collinear spin.

The magnetic symmetry analysis is performed by specifying the magnetic moment(s)
via the optional [--magmom] flag.

Magnetic moment input formats:
   1) Explicitly via comma separated magnetic moments (m1, m2, ..., mn) in the same
      order as the input file
        - Ex: aflow --aflow-sym --magmom=1.0,-1.0,1.0
   2) Read from the VASP INCAR by specifying path to file
        - Ex: aflow --aflow-sym --magmom=path_to_incar
   3) Read from the VASP OUTCAR by specifying path to file
        - Ex: aflow --aflow-sym --magmom=path_to_outcar
   4) Magnetic moment readers for other ab-initio codes will be added in later versions.

   AFLOW autodetects collinear or non-collinear based on the input.

Relevant for the following commands:

  aflow --aflow-sym
  aflow --aflowSG
  aflow --edata
  aflow --equivalent
  aflow --factorgroup
  aflow --lattice_data
  aflow --pointgroup_crystal
  aflow --pointgroup_crystal_data
  aflow --pointgroupkcrystal
  aflow --sgdata
  aflow --sitepointgroup
  aflow --spacegroup
  aflow --wyccar
  aflow --wyckoff

NOTE: This analysis does not account for time reversal symmetry (spin flips), which is
      necessary to find the magnetic space group.  The magnetic symmetry analysis
      here considers the symmetry-breaking feature of the spin-degree of freedom in
      magnetic systems.
      The magnetic space groups will be added in a later version.

****************************************************************************************************
CELL/SETTING CHOICES:

AFLOW-SYM can represent the structure and Wyckoff positions with respect to common cell/setting
choices listed in the ITC with the option flag [--setting=1|2|aflow].
The aflow setting corresponds to the choices listed in the AFLOW Prototype Encyclopedia.

Different cell/setting choice exist for:
  - centrosymmetric space groups with two origin choices (24 space groups, e.g., sg=227):
      --setting=1     : origin centered on high site symmetry (not inversion site)
      --setting=2     : origin centered on inversion site
      --setting=aflow : origin centered on inversion site
  - rhombohedral space groups (e.g., space group 166):
      --setting=1     : rhombohedral setting
      --setting=2     : hexagonal setting
      --setting=aflow : rhombohedral setting
  - monoclinic space groups (e.g., space group 5):
      --setting=1     : unique axis-b
      --setting=2     : unique axis-c
      --setting=aflow : unique axis-b

More information on these setting choices is located in ITC Vol. A on pg. 17 (5th Edition).

Relevant for the following commands (i.e., when the conventional cell/Wyckoff positions is shown):
  aflow --cif
  aflow --edata
  aflow --sgdata
  aflow --wyccar
  aflow --wyckoff


****************************************************************************************************
OUTPUT FORMAT:

AFLOW-SYM prints the results in a standard text output or JavaScript Object Notation (JSON).
The output format is specified by appending the optional flag [--print=<format>].

  Text output (default):
    Ex: aflow --edata --print=txt

  JSON output:
    Ex: aflow --edata --print=json


****************************************************************************************************
AFLOW COMMANDS USAGE:

  aflow --aflow-sym|--AFLOW-SYM|--AFLOWSYM|--aflowSYM|--aflowsym|--full_symmetry|--full_sym|--fullsym[=<tolerance_value>|=tight|=loose] [--magmom=[m1,m2,...|INCAR|OUTCAR]] [--no_scan] [--print=text|json] [--screen_only] < file
        Calculates the full symmetry of a crystal: point group lattice, point group crystal, factor group,
        space group, site point group, inequivalent/equivalent atoms, and the point group of the klattice.
        This routine calculates the full symmetry suite at a certain tolerance. If symmetry rules are broken,
        the code automatically changes tolerance and recalculates from the beginning.

        Options:

          [=<tolerance_value>|=tight|=loose]       : Specifies the symmetry tolerance: <tolerance_value>: user-input, "tight": min_interatomic_dist/100, "loose": min_interatomic_dist/10 (option is appended to command without spaces)
          [--magmom=[m1,m2,...|INCAR|OUTCAR]]      : Perform crystal spin symmetry analysis, where spin degree of freedom breaks symmetry
          [--no_scan]                              : Suppresses tolerance scan.
          [--print=text|json]                      : Outputs the results into human-readable text (default) or JSON, respectively.
          [--screen_only]                          : Prints symmetry elements to string and does not write to a file.

  aflow --aflowSG|--space_group|--sg|--aflowSG_number|--space_group_number|--sg_number|--aflowSG_label|--space_group_label|--sg_label[=<tolerance_value>|=tight|=loose] [--magmom=[m1,m2,...|INCAR|OUTCAR]] [--no_scan] [--tolerance_spectrum=start:stop:nsteps] < file
        Calculates the space group of the crystal in the ITC convention.
        This routine consists of an adaptive tolerance; if any symmetry
        rules are broken, the tolerance is changed and everything is recalculated.
        The space group number and international symbol are printed by default.
        The commands ending in "_number" and "_label" only return the space group number and label, respectively.
        Can be called inside aflow via:
              xstructure.SpaceGroup_ITC()        (default tolerance)
              xstructure.SpaceGroup_ITC(tol)     (user defined tolerance)
        This is the AFLOW counterpart to platonSG and findsymSG.

        Options:

          [=<tolerance_value>|=tight|=loose]       : Specifies the symmetry tolerance: <tolerance_value>: user-input, "tight": min_interatomic_dist/100, "loose": min_interatomic_dist/10 (option is appended to command without spaces)
          [--magmom=[m1,m2,...|INCAR|OUTCAR]]      : Performs crystal spin symmetry analysis, where spin degree of freedom breaks symmetry
          [--no_scan]                              : Suppresses tolerance scan.
          [--tolerance_spectrum=start:stop:nsteps] : Performs symmetry analysis for a range of tolerance values

  aflow --cif[=<tolerance_value>|=tight|=loose] [--no_scan] [--no_symmetry] [--setting=1|2|aflow] < file
        Converts a geometry file to a Crystallographic Information File (CIF).
        By default, a symmetry analysis is performed via AFLOW-SYM; however, this can be suppresed
        with --no_symmetry, printing the positions explicitly (i.e., representing as space group 1).

        Options:

          [=<tolerance_value>|=tight|=loose]       : Specifies the symmetry tolerance: <tolerance_value>: user-input, "tight": min_interatomic_dist/100, "loose": min_interatomic_dist/10 (option is appended to command without spaces)
          [--no_scan]                              : Suppresses tolerance scan.
          [--no_symmetry]                          : Ignores symmetry analysis and returns the CIF with space group 1.
          [--setting=1|2|aflow]                    : Specifies the space group setting. The `aflow' setting follows the choices of the AFLOW Prototype Encyclopedia.
                                                     Different settings exist for the following space groups:
                                                       - centrosymmetric space groups with two origin choices (24 space groups, e.g., sg=227):
                                                           --setting=1     : origin centered on high site symmetry (not inversion site)
                                                           --setting=2     : origin centered on inversion site
                                                           --setting=aflow : origin centered on inversion site
                                                       - rhombohedral space groups (e.g., space group 166):
                                                           --setting=1     : rhombohedral setting
                                                           --setting=2     : hexagonal setting
                                                           --setting=aflow : rhombohedral setting
                                                       - monoclinic space groups (e.g., space group 5):
                                                           --setting=1     : unique axis-b
                                                           --setting=2     : unique axis-c
                                                           --setting=aflow : unique axis-b

  aflow --edata[=<tolerance_value>|=tight|=loose] [--magmom=[m1,m2,...|INCAR|OUTCAR]] [--no_scan] [--print=text|json] [--setting=1|2|aflow] [--suppress_Wyckoff|--suppress_Wyckoff_printing] < file
        Calculates the extended crystallographic data information (lattice, crystal, reciprocal, and superlattice descriptors).
        This routine consists of symmetry consistency checks. If symmetry rules are broken,
        the code automatically changes tolerance and recalculates from the beginning.

        Options:

          [=<tolerance_value>|=tight|=loose]       : Specifies the symmetry tolerance: <tolerance_value>: user-input, "tight": min_interatomic_dist/100, "loose": min_interatomic_dist/10 (option is appended to command without spaces)
          [--magmom=[m1,m2,...|INCAR|OUTCAR]]      : Perform crystal spin symmetry analysis, where spin degree of freedom breaks symmetry
          [--no_scan]                              : Suppresses tolerance scan.
          [--print=text|json]                      : Outputs the results into human-readable text (default) or JSON, respectively.
          [--setting=1|2|aflow]                    : Specifies the space group setting. The `aflow' setting follows the choices of the AFLOW Prototype Encyclopedia.
                                                     Different settings exist for the following space groups:
                                                       - centrosymmetric space groups with two origin choices (24 space groups, e.g., sg=227):
                                                           --setting=1     : origin centered on high site symmetry (not inversion site)
                                                           --setting=2     : origin centered on inversion site
                                                           --setting=aflow : origin centered on inversion site
                                                       - rhombohedral space groups (e.g., space group 166):
                                                           --setting=1     : rhombohedral setting
                                                           --setting=2     : hexagonal setting
                                                           --setting=aflow : rhombohedral setting
                                                       - monoclinic space groups (e.g., space group 5):
                                                           --setting=1     : unique axis-b
                                                           --setting=2     : unique axis-c
                                                           --setting=aflow : unique axis-b
          [--suppress_Wyckoff]                     : Suppresses Wyckoff position output (i.e., space group info only)

  aflow --equivalent | --equiv | --iatoms[=<tolerance_value>|=tight|=loose] [--magmom=[m1,m2,...|INCAR|OUTCAR]] [--no_scan] [--print=text|json] < file
        Calculate the point/factor/space groups to identify equivalent
        and inequivalent atoms. In the output, each structure atom has
        number_label_of_atom[equivalent_to_label]* (* if inequivalent). (SC201107).
        This routine consists of symmetry consistency checks. If symmetry rules are broken,
        the code automatically changes tolerance and recalculates from the beginning.

        Options:

          [=<tolerance_value>|=tight|=loose]       : Specifies the symmetry tolerance: <tolerance_value>: user-input, "tight": min_interatomic_dist/100, "loose": min_interatomic_dist/10 (option is appended to command without spaces)
          [--magmom=[m1,m2,...|INCAR|OUTCAR]]      : Perform crystal spin symmetry analysis, where spin degree of freedom breaks symmetry
          [--no_scan]                              : Suppresses tolerance scan.
          [--print=text|json]                      : Outputs the results into human-readable text (default) or JSON, respectively.

  aflow --factorgroup|--fgroup[=<tolerance_value>|=tight|=loose] [--magmom=[m1,m2,...|INCAR|OUTCAR]] [--no_scan] [--print=text|json] [--screen_only] < file
        Calculates factor group symmetry of the cell {R|t} and writes
        in the aflow.fgroup.out file. See documentation of aflow.
        The point group is required for the factor group, therefore
        the aflow.pgroup.out file will be generated as well.
        This routine consists of symmetry consistency checks. If symmetry rules are broken,
        the code automatically changes tolerance and recalculates from the beginning.

        Options:

          [=<tolerance_value>|=tight|=loose]       : Specifies the symmetry tolerance: <tolerance_value>: user-input, "tight": min_interatomic_dist/100, "loose": min_interatomic_dist/10 (option is appended to command without spaces)
          [--magmom=[m1,m2,...|INCAR|OUTCAR]]      : Perform crystal spin symmetry analysis, where spin degree of freedom breaks symmetry
          [--no_scan]                              : Suppresses tolerance scan.
          [--print=text|json]                      : Outputs the results into human-readable text (default) or JSON, respectively.
          [--screen_only]                          : Prints symmetry elements to string and does not write to a file.
  
  aflow --lattice_data|--data_lattice|--real_lattice_data|--data_real_lattice[=<tolerance_value>|=tight|=loose] [--magmom=[m1,m2,...|INCAR|OUTCAR]] [--no_scan] [--print=text|json] < file
        Returns the lattice type data, including the lattice vectors, lattice parameters (Angstroms/Bohrs and degrees)
        volume, c/a ratio, lattice type, lattice variation, lattice system, and Pearson symbol.

        Options:

          [=<tolerance_value>|=tight|=loose]       : Specifies the symmetry tolerance: <tolerance_value>: user-input, "tight": min_interatomic_dist/100, "loose": min_interatomic_dist/10 (option is appended to command without spaces)
          [--magmom=[m1,m2,...|INCAR|OUTCAR]]      : Perform crystal spin symmetry analysis, where spin degree of freedom breaks symmetry
          [--no_scan]                              : Suppresses tolerance scan.
          [--print=text|json]                      : Outputs the results into human-readable text (default) or JSON, respectively.

  aflow --lattice_type | --lattice | --lattice_crystal[=<tolerance_value>|=tight|=loose] < file
        Returns the lattice type and the conventional lattice type of the CRYSTAL following the tables
        of Setyawayn-Curtarolo [http://dx.doi.org/10.1016/j.commatsci.2010.05.010].

        Options:

          [=<tolerance_value>|=tight|=loose]       : Specifies the symmetry tolerance: <tolerance_value>: user-input, "tight": min_interatomic_dist/100, "loose": min_interatomic_dist/10 (option is appended to command without spaces)

  aflow --lattice_lattice_type | --lattice_lattice[=<tolerance_value>|=tight|=loose] < file
        Returns the lattice type and the conventional lattice type of the LATTICE following the tables
        of Setyawayn-Curtarolo [http://dx.doi.org/10.1016/j.commatsci.2010.05.010].

        Options:

          [=<tolerance_value>|=tight|=loose]       : Specifies the symmetry tolerance: <tolerance_value>: user-input, "tight": min_interatomic_dist/100, "loose": min_interatomic_dist/10 (option is appended to command without spaces)

  aflow --pearson_symbol|--pearson|--Pearson_symbol|--Pearson [--magmom=[m1,m2,...|INCAR|OUTCAR]] [--no_scan] [--tolerance_spectrum=start:stop:nsteps] < file
        Returns the Pearson symbol of the structure.

        Options:
          [=<tolerance_value>|=tight|=loose]       : Specifies the symmetry tolerance: <tolerance_value>: user-input, "tight": min_interatomic_dist/100, "loose": min_interatomic_dist/10 (option is appended to command without spaces)
          [--magmom=[m1,m2,...|INCAR|OUTCAR]]      : Performs crystal spin symmetry analysis, where spin degree of freedom breaks symmetry
          [--no_scan]                              : Suppresses tolerance scan.
          [--tolerance_spectrum=start:stop:nsteps] : Performs symmetry analysis for a range of tolerance values

  aflow --pointgroup | --pgroup[=<tolerance_value>|=tight|=loose] [--no_scan] [--print=text|json] [--screen_only] < file
        Calculates the point group symmetry of the lattice {R} and writes
        it in the aflow.pgroup.out file. See documentation of aflow.
        This routine consists of symmetry consistency checks. If symmetry rules are broken,
        the code automatically changes tolerance and recalculates from the beginning.

        Options:

          [=<tolerance_value>|=tight|=loose]       : Specifies the symmetry tolerance: <tolerance_value>: user-input, "tight": min_interatomic_dist/100, "loose": min_interatomic_dist/10 (option is appended to command without spaces)
          [--no_scan]                              : Suppresses tolerance scan.
          [--print=text|json]                      : Outputs the results into human-readable text (default) or JSON, respectively.
          [--screen_only]                          : Prints symmetry elements to string and does not write to a file.

  aflow --pointgroup_crystal|--pgroup_crystal|--pgroup_xtal[=<tolerance_value>|=tight|=loose] [--magmom=[m1,m2,...|INCAR|OUTCAR]] [--no_scan] [--print=text|json] [--screen_only] < file
        Calculate the point group symmetry of the crystal {R+atoms} and
        writes it the aflow.pgroup_xtal.file.
        This routine consists of symmetry consistency checks. If symmetry rules are broken,
        the code automatically changes tolerance and recalculates from the beginning.

        Options:

          [=<tolerance_value>|=tight|=loose]       : Specifies the symmetry tolerance: <tolerance_value>: user-input, "tight": min_interatomic_dist/100, "loose": min_interatomic_dist/10 (option is appended to command without spaces)
          [--magmom=[m1,m2,...|INCAR|OUTCAR]]      : Perform crystal spin symmetry analysis, where spin degree of freedom breaks symmetry
          [--no_scan]                              : Suppresses tolerance scan.
          [--print=text|json]                      : Outputs the results into human-readable text (default) or JSON, respectively.
          [--screen_only]                          : Prints symmetry elements to string and does not write to a file.
  
  aflow --point_group_crystal_data|--pointgroup_crystal_data|--pgroupxtal_data|--pgroup_xtal_data[=<tolerance_value>|=tight|=loose] [--magmom=[m1,m2,...|INCAR|OUTCAR]] [--no_scan] [--print=text|json] < file
        Returns the crystallographic point group data, including the crystal family, crystal system, crystal class,
        point group symbols (Hermann Mauguin and Schoenflies), point group orbifold, point group type, point group
        order, and point group structure.

        Options:

          [=<tolerance_value>|=tight|=loose]       : Specifies the symmetry tolerance: <tolerance_value>: user-input, "tight": min_interatomic_dist/100, "loose": min_interatomic_dist/10 (option is appended to command without spaces)
          [--magmom=[m1,m2,...|INCAR|OUTCAR]]      : Perform crystal spin symmetry analysis, where spin degree of freedom breaks symmetry
          [--no_scan]                              : Suppresses tolerance scan.
          [--print=text|json]                      : Outputs the results into human-readable text (default) or JSON, respectively.

  aflow --pointgroupkcrystal|--pointgroupk_xtal[=<tolerance_value>|=tight|=loose] [--magmom=[m1,m2,...|INCAR|OUTCAR]] [--no_scan] [--print=text|json] [--screen_only] < file
        Calculates the dual of the crystal point group symmetry operations and writes
        it in the aflow.pgroupk_xtal.out file. See documentation of aflow.
        This routine consists of symmetry consistency checks. If symmetry rules are broken,
        the code automatically changes tolerance and recalculates from the beginning.

        Options:

          [=<tolerance_value>|=tight|=loose]       : Specifies the symmetry tolerance: <tolerance_value>: user-input, "tight": min_interatomic_dist/100, "loose": min_interatomic_dist/10 (option is appended to command without spaces)
          [--magmom=[m1,m2,...|INCAR|OUTCAR]]      : Perform crystal spin symmetry analysis, where spin degree of freedom breaks symmetry
          [--no_scan]                              : Suppresses tolerance scan.
          [--print=text|json]                      : Outputs the results into human-readable text (default) or JSON, respectively.
          [--screen_only]                          : Prints symmetry elements to string and does not write to a file.

  aflow --pointgroupklattice|--pgroupk[=<tolerance_value>|=tight|=loose] [--no_scan] [--print=text|json] [--screen_only] < file
        Calculates the point group symmetry of the klattice {K} and writes
        it in the aflow.pgroupk.out file. See documentation of aflow.
        This routine consists of symmetry consistency checks. If symmetry rules are broken,
        the code automatically changes tolerance and recalculates from the beginning.

        Options:

          [=<tolerance_value>|=tight|=loose]       : Specifies the symmetry tolerance: <tolerance_value>: user-input, "tight": min_interatomic_dist/100, "loose": min_interatomic_dist/10 (option is appended to command without spaces)
          [--no_scan]                              : Suppresses tolerance scan.
          [--print=text|json]                      : Outputs the results into human-readable text (default) or JSON, respectively.
          [--screen_only]                          : Prints symmetry elements to string and does not write to a file.

  aflow --pointgroupk_Patterson|--pgroupk_Patterson[=<tolerance_value>|=tight|=loose] [--magmom=[m1,m2,...|INCAR|OUTCAR] [--no_scan] [--print=text|json] [--screen_only] < file
        Calculates the Patterson point group symmetry operations and writes
        it in the aflow.pgroupk_Patterson.out file. See documentation of aflow.
        This routine consists of symmetry consistency checks. If symmetry rules are broken,
        the code automatically changes tolerance and recalculates from the beginning.

        Options:

          [=<tolerance_value>|=tight|=loose]       : Specifies the symmetry tolerance: <tolerance_value>: user-input, "tight": min_interatomic_dist/100, "loose": min_interatomic_dist/10 (option is appended to command without spaces)
          [--magmom=[m1,m2,...|INCAR|OUTCAR]]      : Perform crystal spin symmetry analysis, where spin degree of freedom breaks symmetry
          [--no_scan]                              : Suppresses tolerance scan.
          [--print=text|json]                      : Outputs the results into human-readable text (default) or JSON, respectively.
          [--screen_only]                          : Prints symmetry elements to string and does not write to a file.
  
   aflow --reciprocal_lattice_data|--reciprocallattice_data|--klattice_data|--data_reciprocal_lattice[=<tolerance_value>|=tight|=loose] [--no_scan] [--print=text|json] < file
        Returns the reciprocal lattice type data, including the lattice vectors, lattice parameters (Angstroms/Bohrs and degrees)
        volume, c/a ratio, lattice type, lattice variation, lattice system, and Pearson symbol.

        Options:

          [=<tolerance_value>|=tight|=loose]       : Specifies the symmetry tolerance: <tolerance_value>: user-input, "tight": min_interatomic_dist/100, "loose": min_interatomic_dist/10 (option is appended to command without spaces)
          [--no_scan]                              : Suppresses tolerance scan.
          [--print=text|json]                      : Outputs the results into human-readable text (default) or JSON, respectively.

  aflow --sc|--std_conv|--standard_conventional|--sconv < file
   	    Output POSCAR in a standard conventional lattice
        DOI: 10.1016/j.commatsci.2010.05.010

  aflow --sgdata|--space_group_data[=<tolerance_value>|=tight|=loose] [--magmom=[m1,m2,...|INCAR|OUTCAR]] [--no_scan] [--print=text|json] [--setting=1|2|aflow] [--suppress_Wyckoff|--suppress_Wyckoff_printing] < file
        Outputs the following space group information: space group number, space group labels
        (Hermann Mauguin, Hall, Schoenflies), Laue class, crystal class, ITC setting, ITC origin shift,
        general Wyckoff position, and the WYCCAR (similar to --edata).
        This routine consists of symmetry consistency checks. If symmetry rules are broken,
        the code automatically changes tolerance and recalculates from the beginning.

        Options:

          [=<tolerance_value>|=tight|=loose]       : Specifies the symmetry tolerance: <tolerance_value>: user-input, "tight": min_interatomic_dist/100, "loose": min_interatomic_dist/10 (option is appended to command without spaces)
          [--magmom=[m1,m2,...|INCAR|OUTCAR]]      : Perform crystal spin symmetry analysis, where spin degree of freedom breaks symmetry
          [--no_scan]                              : Suppresses tolerance scan.
          [--print=text|json]                      : Outputs the results into human-readable text (default) or JSON, respectively.
          [--setting=1|2|aflow]                    : Specifies the space group setting. The `aflow' setting follows the choices of the AFLOW Prototype Encyclopedia.
                                                     Different settings exist for the following space groups:
                                                       - centrosymmetric space groups with two origin choices (24 space groups, e.g., sg=227):
                                                           --setting=1     : origin centered on high site symmetry (not inversion site)
                                                           --setting=2     : origin centered on inversion site
                                                           --setting=aflow : origin centered on inversion site
                                                       - rhombohedral space groups (e.g., space group 166):
                                                           --setting=1     : rhombohedral setting
                                                           --setting=2     : hexagonal setting
                                                           --setting=aflow : rhombohedral setting
                                                       - monoclinic space groups (e.g., space group 5):
                                                           --setting=1     : unique axis-b
                                                           --setting=2     : unique axis-c
                                                           --setting=aflow : unique axis-b
          [--suppress_Wyckoff]                     : Suppresses Wyckoff position output (i.e., space group info only)

  aflow --sitepointgroup|--agroup[=<tolerance_value>|=tight|=loose] [--magmom=[m1,m2,...|INCAR|OUTCAR]] [--no_scan] [--print=text|json] [--screen_only] < file
        Calculates the site point group symmetry for every atom in the
        unit cell and writes it in the aflow.agroup.out file.
        This routine consists of symmetry consistency checks. If symmetry rules are broken,
        the code automatically changes tolerance and recalculates from the beginning.

        Options:

          [=<tolerance_value>|=tight|=loose]       : Specifies the symmetry tolerance: <tolerance_value>: user-input, "tight": min_interatomic_dist/100, "loose": min_interatomic_dist/10 (option is appended to command without spaces)
          [--magmom=[m1,m2,...|INCAR|OUTCAR]]      : Perform crystal spin symmetry analysis, where spin degree of freedom breaks symmetry
          [--no_scan]                              : Suppresses tolerance scan.
          [--print=text|json]                      : Outputs the results into human-readable text (default) or JSON, respectively.
          [--screen_only]                          : Prints symmetry elements to string and does not write to a file.

  aflow --sp|--std_prim|--standard_primitive|--sprim < file
   	    Output POSCAR in a standard primitive lattice
        DOI: 10.1016/j.commatsci.2010.05.010

  aflow --spacegroup[=<tolerance_value>|=tight|=loose] [--magmom=[m1,m2,...|INCAR|OUTCAR]] [--no_scan] [--print=text|json] [--radius=<radius>] [--screen_only] < file
        Calculates space group symmetry of the cell {R|t+T} with
        translations as big as |T| and writes it in the aflow.sgroup.out
        file. See documentation of aflow.
        Be careful because the size of the space group increases
        as the radius^3 times the size of the factor group.
        The point and factor groups are required for the space group,
        therefore the aflow.pgroup.out and aflow.fgroup.out files will
        be generated as well.
        This routine consists of symmetry consistency checks. If symmetry rules are broken,
        the code automatically changes tolerance and recalculates from the beginning.

        Options:

          [=<tolerance_value>|=tight|=loose]       : Specifies the symmetry tolerance: <tolerance_value>: user-input, "tight": min_interatomic_dist/100, "loose": min_interatomic_dist/10 (option is appended to command without spaces)
          [--magmom=[m1,m2,...|INCAR|OUTCAR]]      : Perform crystal spin symmetry analysis, where spin degree of freedom breaks symmetry
          [--no_scan]                              : Suppresses tolerance scan.
          [--print=text|json]                      : Outputs the results into human-readable text (default) or JSON, respectively.
          [--radius=<radius>]                      : Specifies radius (Angstroms) to print symmetry operations out to a given distance from the origin
          [--screen_only]                          : Prints symmetry elements to string and does not write to a file.
  
  aflow --superlattice_data|--data_superlattice[=<tolerance_value>|=tight|=loose] [--no_scan] [--print=text|json] < file
        Returns the superlattice type data, including the lattice vectors, lattice parameters (Angstroms/Bohrs and degrees)
        volume, c/a ratio, lattice type, lattice variation, lattice system, and Pearson symbol.

        Options:

          [=<tolerance_value>|=tight|=loose]       : Specifies the symmetry tolerance: <tolerance_value>: user-input, "tight": min_interatomic_dist/100, "loose": min_interatomic_dist/10 (option is appended to command without spaces)
          [--no_scan]                              : Suppresses tolerance scan.
          [--print=text|json]                      : Outputs the results into human-readable text (default) or JSON, respectively.

  aflow --wyccar[=<tolerance_value>|=tight|=loose] [--magmom=[m1,m2,...|INCAR|OUTCAR]] [--no_scan] [--setting=1|2|aflow] < file
        Prints the wyccar file (the poscar file with wyckoff positions) for the input POSCAR
        The wyccar format presents the lattice in standard conventional form and the wyckoff
        positions with their associated labels and site symmetry. If variability in the
        Wyckoff positions exists, the algorithm finds the Wykoff positions with the "smallest"
        lettering scheme ("a" being the smallest).

        Options:

          [=<tolerance_value>|=tight|=loose]       : Specifies the symmetry tolerance: <tolerance_value>: user-input, "tight": min_interatomic_dist/100, "loose": min_interatomic_dist/10 (option is appended to command without spaces)
          [--magmom=[m1,m2,...|INCAR|OUTCAR]]      : Perform crystal spin symmetry analysis, where spin degree of freedom breaks symmetry
          [--no_scan]                              : Suppresses tolerance scan.
          [--setting=1|2|aflow]                    : Specifies the space group setting. The `aflow' setting follows the choices of the AFLOW Prototype Encyclopedia.
                                                     Different settings exist for the following space groups:
                                                       - centrosymmetric space groups with two origin choices (24 space groups, e.g., sg=227):
                                                           --setting=1     : origin centered on high site symmetry (not inversion site)
                                                           --setting=2     : origin centered on inversion site
                                                           --setting=aflow : origin centered on inversion site
                                                       - rhombohedral space groups (e.g., space group 166):
                                                           --setting=1     : rhombohedral setting
                                                           --setting=2     : hexagonal setting
                                                           --setting=aflow : rhombohedral setting
                                                       - monoclinic space groups (e.g., space group 5):
                                                           --setting=1     : unique axis-b
                                                           --setting=2     : unique axis-c
                                                           --setting=aflow : unique axis-b

  aflow --Wyckoff|--Wyckoff_positions|--wyckoff|--wyckoff_positions[=<tolerance_value>|=tight|=loose] [--letters_only] [--magmom=[m1,m2,...|INCAR|OUTCAR]] [--multiplicities_only] [--no_scan] [--setting=1|2|aflow] [--site_symmetries_only] < file
        Prints the Wyckoff positions (multiplicity, letter, site symmetry and representative
        coordinate) based on the input geometry. If variability in the Wyckoff positions exists,
        the algorithm finds the Wykoff positions with the "smallest" lettering scheme ("a" being the smallest).

        Options:

          [=<tolerance_value>|=tight|=loose]       : Specifies the symmetry tolerance: <tolerance_value>: user-input, "tight": min_interatomic_dist/100, "loose": min_interatomic_dist/10 (option is appended to command without spaces)
          [--letters_only]                         : Prints the Wyckoff letters only (colons separate elements, commas separate Wyckoff positions for a given element)
          [--magmom=[m1,m2,...|INCAR|OUTCAR]]      : Perform crystal spin symmetry analysis, where spin degree of freedom breaks symmetry
          [--multiplicities_only]                  : Prints the Wyckoff multiplicities only (colons separate elements, commas separate Wyckoff positions for a given element)
          [--no_scan]                              : Suppresses tolerance scan.
          [--setting=1|2|aflow]                    : Specifies the space group setting. The `aflow' setting follows the choices of the AFLOW Prototype Encyclopedia.
                                                     Different settings exist for the following space groups:
                                                       - centrosymmetric space groups with two origin choices (24 space groups, e.g., sg=227):
                                                           --setting=1     : origin centered on high site symmetry (not inversion site)
                                                           --setting=2     : origin centered on inversion site
                                                           --setting=aflow : origin centered on inversion site
                                                       - rhombohedral space groups (e.g., space group 166):
                                                           --setting=1     : rhombohedral setting
                                                           --setting=2     : hexagonal setting
                                                           --setting=aflow : rhombohedral setting
                                                       - monoclinic space groups (e.g., space group 5):
                                                           --setting=1     : unique axis-b
                                                           --setting=2     : unique axis-c
                                                           --setting=aflow : unique axis-b
          [--site_symmetries_only]                 : Prints the Wyckoff site symmetries only (colons separate elements, commas separate Wyckoff positions for a given element)


For additional information contact: David Hicks (david.hicks@duke.edu)

****************************************************************************************************
*                                                                                                  *
*                        aflow - STEFANO CURTAROLO Duke University 2003-2021                       *
*                            High-Throughput ab-initio Computing Project                           *
*                                                                                                  *
****************************************************************************************************
