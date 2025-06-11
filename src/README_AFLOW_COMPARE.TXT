********************************************************************************
*                                                                              *
*             Aflow STEFANO CURTAROLO - Duke University 2003-2021              *
*                High-Throughput ab-initio Materials Discovery                 *
*                                                                              *
********************************************************************************
LATEST VERSION OF THE FILE:    materials.duke.edu/AFLOW/README_AFLOW_COMPARE.TXT
********************************************************************************

AFLOW-XtalFinder (COMPARE STRUCTURE) README

README written by: David Hicks (david.hicks@duke.edu)

Citation info: D. Hicks, C. Toher, D.C. Ford, F. Rose,
               C. De Santo O. Levy, M.J. Mehl, and S. Curtarolo,
               AFLOW-XtalFinder: a reliable choice to identify crystalline
               prototypes, npj Comput. Mater. 7, 30 (2020).
               [doi=10.1038/s41524-020-00483-4]

********************************************************************************
GENERAL OVERVIEW

AFLOW-XtalFinder identifies and classifies crystallographic prototypes, regardless of their representation.
The code maps structures to their ideal prototype designation, determining their internal degrees of
freedom consistent with the International Tables for Crystallography (ITC).

To ensure uniqueness, structures are compared based on
symmetry (isopointal analysis),
atom environment (isoconfigurational snapshots), and
geometric structure (complete isoconfigurational analysis).
Structural similarity is quantified via a crystal misfit criteria.

The AFLOW-XtalFinder functionalities are described below.

----------------------------------------------------------------------
PROTOTYPING TOOLS
Structures are automatically cast into the AFLOW standard representation
consistent with the AFLOW Prototype Encyclopedia.
The routine determines the internal degrees of freedom for the prototype,
e.g., esseneite:

AFLOW label   : ABC6D2_mC40_15_e_e_3f_f
params        : a,b/a,c/a,beta,y1,y2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6
params values : 4.3814,0.901127,0.548523,105.81,0.6918,0.0942,0.1112,0.9123,0.1341,0.3635,0.7411,0.3201,0.3532,0.9802,0.9876,0.2883,0.9057,0.2236

----------------------------------------------------------------------
QUANTITATIVE SIMILARITY METRIC
The code returns a similarity metric (i.e. misfit value) between the two structures.
Particular ranges of the misfit are adopted from Burzlaff and Malinovsky's work, characterizing structural similarity:
  - misfit <= epsilon_match                  : structures are similar/compatible
  - epsilon_match < misfit <= epsilon_family : structures are within the same family
  - misfit > epsilon_family                  : structures are not compatible/no match found
By default:
  epsilon_match = 0.1 and
  epsilon_family = 0.2
In general, structures that match within 0.1 will have enthalpies within 5 meV/atom.
Structures that match beyond 0.2 have poor structure-property correlation.
See XtalFinder article for details.

----------------------------------------------------------------------
SUPER-TYPE COMPARISON MODES
The module offers four comparison modes:

  1) material-type   : maps same element types onto one another, revealing duplicate compounds
  2) structure-type  : maps structures regardless of the element types, identifying crystallographic prototypes
  3) decoration-type : creates and compares all possible atom decorations for a given structure, determining unique and equivalent decorations
  4) magnetic-type   : maps compounds by element types and magnetic moments, distinguishing distinct spin configurations

----------------------------------------------------------------------
MULTIPLE COMPARISONS
AFLOW-XtalFinder can compare multiple compounds in a single command.
Structures are automatically grouped into comparison sets via stoichiometry,
species (if applicable), symmetry (space group number and Wyckoff positions), and atomic environment.
Compounds are compared until each structure is matched or all possible comparisons have been exhausted.

See AFLOW COMMANDS USAGE for more details.

----------------------------------------------------------------------
INPUT FILE TYPES
AFLOW-XtalFinder automatically detects the input geometry file type.
The supported formats are:
  - VASP (POSCAR or CONTCAR)
  - QuantumESPRESSO
  - FHI-AIMS
  - ABINIT
  - CIF
  - ELK

----------------------------------------------------------------------
COMPARE TO ESTABLISHED STRUCTURES
Users can compare input structures to the AFLOW.org repository and the AFLOW Prototype Encyclopedia.
The functionality is useful for identifying equivalent compounds within repositories and relating their properties.

----------------------------------------------------------------------
MULTITHREADED CAPABILITIES
Many of the AFLOW-XtalFinder functions are parallelized: symmetry analysis, environment analysis, and comparisons.
Use the --np=<number> flag to specify the number of threads for multithreaded comparisons.

****************************************************************************************************
LIST OF AFLOW-XtalFinder COMMANDS:

  aflow --compare2database          : Returns AFLOW database entries that match with an input structure
  aflow --compare2prototypes        : Returns AFLOW prototypes that match with an input structure
  aflow --compare_database_entries  : Compares entries in the AFLOW.org repository
  aflow --compare_materials         : Compares materials (mapping alike atomic species) and
                                      returns their level of similarity
  aflow --compare_structures        : Compares structures (ignore species in mapping) and returns their level
                                      of similarity
  aflow --get_isopointal_prototypes : Returns all isopointal prototypes (i.e., similar symmetry) to an input structure
  aflow --prototype                 : Returns AFLOW label/parameters/parameter values of the structure
  aflow --unique_atom_decorations   : Returns the unique atom decorations for the structure

For details on usage, see the "AFLOW COMMANDS USAGE" section below.

********************************************************************************
NOTE ON STRUCTURE SIZE VS SPEED

The number of atoms in the system generally increases the computational expense
of the comparisons.  Namely, it is dependent on the least frequently occurring
atom (LFA) in the crystal. If LFA>20 or there are many LFA types, then the comparison
time will increase. By default, the algorithm terminates as soon as a match is found
(unless you request the optimal match via the --optimize_match option).
However, this does not improve the speed of the comparison when a match cannot be found,
since we must explore the entire search space.

METHODS FOR ENHANCING COMPUTATION SPEED:
  1) Increase the number of processors (--np=xx)
  2) Convert structures to a primitive and/or more orthogonal lattice representations
     via the following options (conversions can be combined):
       --primitivize : convert the crystal to a primitive representation
       --minkowski   : convert the crystal to the Minkowski representation
       --niggli      : convert the crystal to the Niggli representation

********************************************************************************
OVERVIEW OF COMPARISON ALGORITHM

Steps:
    1) Input two structures: structure1 and structure2
    2) Scale volumes of two structures to be commensurate
    3) Determine the least frequently occurring atom (LFA)
    4) Shift LFA for both structures to their respective origins;
       these will be used to indicate our lattice
    5) Create a supercell grid of structure2
    6) Search for all possible lattices in the
       supercell and see if we can match it with structure1
       using Burzlaff's criteria. (This part is parallel if more than one
       processor is specified.)
    7) Once a commensurate lattice is found (i.e., similar lattice vector lengths
       and volume to structure1), transform structure2 to its new
       representation.
    8) Find the best match (smallest misfit (mis)):
       If mis <= epsilon_match:
         Structures are similar. Print out the new representation of
         structure2 and the figure of misfit
       else if epsilon_match < mis <= epsilon_family:
         Structures are in the same family (possible symmetric
         group-subgroup relation). Print out the new representation
         of structure2 and the figure of misfit
       else mis>0.2:
         Structures are not the same.

********************************************************************************
OVERVIEW OF MULTIPLE COMPARISON SCHEME

Steps:
    1) Read in all files in file list, directory, or structures concatenated into a file.
    2) Obtain stoichiometry/species/space group.
    3) Compare structures that have the same stoichiometry, species (if applicable), space group, and
       LFA environment by picking a representative structure and comparing it with the other structures.
    4) If any of the structures did not match the representative, create another group of
       comparisons for the ones that did not match. Continue until all structures are
       matched or until all comparisons are exhausted.
    5) Print out results in both a JSON file (for easy integration into a Python workflow) and text file.

********************************************************************************
AFLOW COMMANDS USAGE:

The AFLOW commands are listed below. Command-specific options are listed immediately below the command.
General options - applicable to all commands - is listed at the bottom of this section.

    aflow --compare2database [GENERAL_COMPARISON_OPTIONS] [COMPARE2DATABASE_OPTIONS] < file
        Compares a structure (file) to entries in the AFLOW database, returning similar
        compounds and quantifying their levels of similarity (misfit values).
        Material properties can be extracted from the database and printed (via AFLUX), highlighting structure-property relationships.
        Performs material-type comparisons or structure-type comparisons (by adding the --structure_comparison option).

        [COMPARE2DATABASE_OPTIONS] : Options specific to this command:

            [--catalog=lib1|lib2|lib3|lib4|lib6|lib7|icsd]   : Restrict the entries to compare with to a specific catalog/library
                                                               (e.g., `lib1', `lib2', `lib3', `icsd', etc.).
            [--properties=enthalpy_atom,natoms,...]          : Specify the properties via their API keyword to print the corresponding
                                                               values with the comparison results.
            [--relaxation_step=original|relax1|most_relaxed] : Compare geometries from a particular step in the DFT relaxation.
                                                               Allows numeric or string input (e.g., 0->original, 1->relax1, 2->most_relaxed).


    aflow --compare_database_entries [GENERAL_COMPARISON_OPTIONS] [COMPARE_DATABASE_ENTRIES_OPTIONS]
        Compares structures in the AFLOW.org repository and groups entries into structurally similar materials.
        Entries can be queried by arity, alloy systems, stoichiometry, and space group.
        Material properties can be extracted from the database and printed (via AFLUX), highlighting structure-property relationships.
        Performs material-type comparisons or structure-type comparisons (by adding the --structure_comparison option).

        [COMPARE_DATABASE_ENTRIES_OPTIONS] : Options specific to this command:

            [--alloy=AgAlMn...]                              : Specify the alloy system(s) to extract from the database and compare.
            [--nspecies=3]                                   : Specify the arity group (number of species) to extract from the database and compare.
            [--catalog=lib1|lib2|lib3|lib4|lib6|lib7|icsd]   : Restrict the entries to compare with to a specific catalog/library
                                                               (e.g., `lib1', `lib2', `lib3', `icsd', etc.).
            [--properties=enthalpy_atom,natoms,...]          : Specify the properties via their API keyword to print the corresponding
                                                               values with the comparison results.
            [--relaxation_step=original|relax1|most_relaxed] : Compare geometries from a particular step in the DFT relaxation.
                                                               Allows numeric or string input (e.g., 0->original, 1->relax1, 2->most_relaxed).
            [--space_group=225,186,227,...]                  : Specify the space group(s) of the structures to extract and compare.
                                                               More than one space group can be specified.
            [--stoichiometry=1,2,3,...]                      : Specify the stoichiometric ratio of the structures to extract and compare.


    aflow --compare2prototypes [GENERAL_COMPARISON_OPTIONS] [COMPARE2PROTOTYPES_OPTIONS] < file
        Compares a structure (file) to those in the AFLOW Prototype Encyclopedia, returning the similar
        prototypes and quantifying their levels of similarity (misfit values).

        [COMPARE2PROTOTYPES_OPTIONS] : Options specific to this command:

            [--catalog=aflow|htqc|all] : Restrict the entries to compare with to a specific catalog/library (e.g., `aflow' or `htqc').


    aflow --compare_materials=<f1>,<f2>,...| -D <path> | -F=<filename> [GENERAL_COMPARISON_OPTIONS]
        Compares compounds comprised of the same atomic species and with
        commensurate stoichiometric ratios, i.e., material-type comparison, and returns their level of similarity (misfit value).
        This method identifies unique and duplicate materials.
        There are three input types:
            aflow --compare_materials=<f1>,<f2>,...  : append geometry files (<f1>,<f2>,...) to compare,
            aflow --compare_materials -D <path>      : specify path to directory (<path>) containing geometry files to compare, and
            aflow --compare_materials -F=<filename>  : specify file (<filename>) containing compounds between delimiters
                                                       [VASP_POSCAR_MODE_EXPLICIT]START and VASP_POSCAR_MODE_EXPLICIT]STOP.
                                                       Additional delimiters will be included in later versions.


    aflow --compare_structures=<f1>,<f2>,...| -D <path> | -F=<filename> [GENERAL_COMPARISON_OPTIONS]
        Compares compounds with commensurate stoichiometric ratios and no
        requirement of the atomic species, i.e., structure-type comparison, and returns their level of similarity (misfit value).
        This method identifies unique and duplicate prototypes.
        There are three input types:
            aflow --compare_structures=<f1>,<f2>,...  : append geometry files (<f1>,<f2>,...) to compare,
            aflow --compare_structures -D <path>      : specify path to directory (<path>) containing geometry files to compare, and
            aflow --compare_structures -F=<filename>  : specify file (<filename>) containing compounds between delimiters
                                                       [VASP_POSCAR_MODE_EXPLICIT]START and VASP_POSCAR_MODE_EXPLICIT]STOP.
                                                       Additional delimiters will be included in later versions.


    aflow --prototype [PROTOTYPE_OPTIONS] < file
        Converts a structure (file) into its corresponding standard AFLOW prototype label.
        The parameter variables (degrees of freedom) and parameter values are also listed.
        Information about the label and parameters are described in the AFLOW Prototype Encyclopedia:
        Part 1 (doi=10.1016/j.commatsci.2017.01.017),
        Part 2 (doi=10.1016/j.commatsci.2018.10.043),
        or README_AFLOW_ANRL.TXT.

        [PROTOTYPE_OPTIONS] : Options specific to this command:

            [--setting=1|2|aflow]    : Specify the space group setting for the conventional cell/Wyckoff positions.
                                       The `aflow' setting follows the choices of the AFLOW Prototype Encyclopedia:
                                         - axis-b for monoclinic space groups,
                                         - rhombohedral setting for rhombohedral space groups, and
                                         - origin centered on the inversion for centrosymmetric space groups
                                       (default: aflow).
            [--print_element_names]  : Print the element names from the input structure
            [--print_atomic_numbers] : Print the atomic number for the elements in the input structure
            [--print=TEXT|JSON]      : Print the output into human-readable text or JSON, respectively.


    aflow --get_isopointal_prototype [ISOPOINTAL_OPTIONS] < file
        Returns the prototype labels that are isopointal (i.e. similar symmetry) as the input structure.
        In this case, the same symmetry means that the structures have commensurate:
          1) space groups (equivalent or enantiomorphic pairs) and
          2) Wyckoff positions (same multiplicity and site symmetry)
        Can be used to quickly determine if structures have similar AFLOW prototype designations.

        [ISOPOINTAL_OPTIONS] : Options specific to this command:

            [--catalog=aflow|htqc|all] : Restrict the entries to compare with to a specific catalog/library (e.g., `aflow' or `htqc').
            [--print=TEXT|JSON]      : Print the output into human-readable text or JSON, respectively.


    aflow --unique_atom_decorations [GENERAL_COMPARISON_OPTIONS] < file
        Determines the unique and duplicate atom decorations for a given structure, grouping
        equivalent decorations (if applicable).


    [GENERAL_COMPARISON_OPTIONS] : Generic options for all commands (unless indicated otherwise):

        [--misfit_match=0.1]                    : Specify the misfit threshold for matching structures (default: 0.1).
        [--misfit_family=0.2]                   : Specify the misfit threshold for structures within the same family (default: 0.2).
        [--np=16|--num_proc=16]                 : Allocate the number of processors/threads for the task.
        [--optimize_match]                      : Explore all lattice and origin choices to find the best matching representation, i.e., minimizes misfit value.
        [--no_scale_volume]                     : Suppresses volume rescaling during structure matching; identifies differences due to volume expansion or compression of a structure.
        [--ignore_symmetry]                     : Neglects symmetry analysis (both space group and Wyckoff positions) for grouping comparisons.
        [--ignore_Wyckoff]                      : Neglects Wyckoff symmetry (site symmetry) for filtering comparisons, but considers the space group number.
        [--ignore_environment]                  : Neglects LFA environment analysis for filtering comparisons.
        [--keep_unmatched]                      : Retains misfit information of unmatched structures (i.e., misfit>0.1).
        [--match_to_aflow_prototypes]           : Identifies matching AFLOW prototypes to the representative structure.  The option does not apply to
                                                  --unique_atom_decorations or --compare2prototypes (redundant).
        [--magmom=<m1,m2,...|INCAR|OUTCAR>:...] : Specifies the magnetic moment for each structure (collinear or non-collinear) delimited by colons,
                                                  signaling a magnetic-type comparison. The option does not apply to --compare_structures since the atom type is neglected.
                                                  AFLOW-XtalFinder supports three input formats for the magnetic moment:
                                                   1) explicit declaration via comma-separated string m_{1},m_{2},...m_{n} (m_{1,x},m_{1,y},m_{1,z},m_{2,x},...m_{n,z} for non-collinear)
                                                   2) read from a VASP INCAR, or
                                                   3) read from a VASP OUTCAR.
                                                  Additional magnetic moment readers for other ab initio codes will be available in future versions.
        [--add_aflow_prototype_designation]     : Casts prototype structure into the AFLOW standard designation.  The option does not apply to commands
                                                  --unique_atom_decorations or --prototype! (redundant).
        [--remove_duplicate_compounds]          : For structure-type comparisons, duplicate compounds are identified first
                                                  (via a material-type comparison without volume scaling), then remaining unique compounds are compared, removing duplicate bias.
        [--ICSD]                                : Only if ICSD structures: ensure the representative structure is the entry with the smallest ICSD number
        [--print_mapping|--print_map|--print]   : For comparing two structures, additional comparison information is printed, including atom mappings, distances between matched atoms,
                                                  and the transformed structures in the closest matching representation.
        [--print=TEXT|JSON]                     : For comparing multiple structures, the results printed to the relevant files are formatted into human-readable
                                                  text or JSON, respectively.  By default, AFLOW-XtalFinder writes the output in both formats.
        [--quiet|--q]                           : Suppresses the log information for the comparisons.
        [--screen_only]                         : Prints the comparison results to the screen and does not write to any files.
        [--primitivize]                         : Primitivize structures (enhance comparison speed).
        [--minkowski]                           : Perform Minkowski lattice reduction on the structures (enhance comparison speed).
        [--niggli]                              : Perform Niggli lattice reduction on the structures (enhance comparison speed).


For additional information contact: David Hicks (david.hicks@duke.edu)
// ***************************************************************************

********************************************************************************
*                                                                              *
*             Aflow STEFANO CURTAROLO - Duke University 2003-2021              *
*                High-Throughput ab-initio Materials Discovery                 *
*                                                                              *
********************************************************************************
