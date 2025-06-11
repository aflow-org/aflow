****************************************************************************************************
*                                                                                                  *
*                        aflow - STEFANO CURTAROLO Duke University 2003-2021                       *
*                            High-Throughput ab-initio Computing Project                           *
*                                                                                                  *
****************************************************************************************************
LATEST VERSION OF THE FILE:                           materials.duke.edu/AFLOW/README_AFLOW_ANRL.TXT
****************************************************************************************************

AFLOW/ANRL Prototypes

AFLOW implementation: Stefano Curtarolo and David Hicks

Citation Information:
    - M. J. Mehl, D. Hicks, C. Toher, O. Levy, R. M. Hanson, G. L. W. Hart, and S. Curtarolo,
      The AFLOW Library of Crystallographic Prototypes: Part 1, Comp. Mat. Sci. 136, S1–S828 (2017).
    - D. Hicks, M. J. Mehl, E. Gossett, C. Toher, O. Levy, R. M. Hanson, G. L. W. Hart, and S. Curtarolo,
      The AFLOW Library of Crystallographic Prototypes: Part 2, Comp. Mat. Sci. 161 Supplement, S1-S1011 (2019).
    - D. Hicks, M. J. Mehl, M. Esters, C. Oses, O. Levy, G. L. W. Hart, C. Toher, and S. Curtarolo,
      The AFLOW Library of Crystallographic Prototypes: Part 3, in press, Comp. Mat. Sci. (2021).

********************************************************************************
GENERAL OVERVIEW

The AFLOW/ANRL Library consists of 1,100 (part 1: 288, part 2: 302, part3: 510) crystallographic prototypes from the papers listed above.
These prototypes have been programmed into AFLOW. The list of prototypes and how to call them can be
found in aflow_anrl.cpp. The prototypes are generally listed in the order they appear in the
paper (low space group number to high space group number).

Note: This implementation of the geometric prototypes allows for the creation of different
structure prototypes depending on the parameters chosen. Therefore, these prototypes are quite general.

----------------------------------------------------------------------
[New feature as of V3.2.5]
Any possible prototype can be created with the new symbolic prototype generator as long as it follows
the label and internal degrees of freedom scheme.
Prototypes no longer need to be tabulated (hard-coded) in the AFLOW software, and are now created on-the-fly.

****************************************************************************************************
PROTOTYPE COMMAND SYNTAX:

The general format for the command is as follows (optional arguments are enclosed in square brackets "[]"):

  aflow --proto=label:A:B:C... --params=a,b/a,c/a,alpha,beta,gamma,x1,y1,z1,x2,y2,z2,... [--hex] [--vasp|--aims|--qe|--abinit|--elk|--wyccar] [--add_equations|--equations_only]

Example:

  aflow --proto=ABC4_oP12_16_ag_cd_2u --params=5.61,1.01069518717,1.61319073084,0.2,0.26,0.125,0.74,0.8,0.63

****************************************************************************************************
PROTOTYPE ALIASES:

Particular structure prototypes from the AFLOW Library can be generated with a short-hand alias,
rather than explicitly providing the parameters with --params (new feature as of V3.1.220).

For prototypes with no degrees of freedom (other than the lattice parameter), the alias is the
label itself, since there is no other possible structure that is generated from this label.
  Example: quaternary Heusler

    aflow --proto=ABCD_cF16_216_c_d_b_a

  The Wyckoff positions do not contain any variables.

For prototypes with degrees of freedom, the alias is the label and a number.
  Example: A_hR2_166_c (two degrees of freedom) has three structures in the library

    aflow --proto=A_hR2_166_c-001 :    rhombohedral graphite
    aflow --proto=A_hR2_166_c-002 :    alpha-As
    aflow --proto=A_hR2_166_c-003 :    beta-O

----------------------------------------------------------------------
AUTOMATIC VOLUME SCALING
The alias scheme uses the atomic scaled volume - as opposed to the ANRL lattice parameter -
for more "realistic" volumes when decorating with different atomic species.
To force the ANRL defined lattice parameter, add the following flag to the
command: --use_anrl_lattice_parameter.

By default, AFLOW will check if the input parameters match with any of the predefined aliases/enumerations.
If there is a match, then the label is appended with the appropriate enumeration suffix and will
be printed in the aflow.in and/or geometry file title.

********************************************************************************
DESCRIPTION COMMAND SYNTAX:

  aflow --proto=label:A:B:C... --params=a,b/a,c/a,alpha,beta,gamma,x1,y1,z1,x2,y2,z2,... [--hex] [--vasp|--aims|--qe|--abinit|--elk|--wyccar] [--add_equations|--equations_only]

A description of each input is described below:

  label:      The label is a string that provides a great deal of information regarding the
              geometric prototype. The information is separated by underscores into different
              fields.
                - Field 1: Stoichiometry
                - Field 2: Pearson symbol
                - Field 3: Space group number
                - Field 4: Wyckoff letters (in alphabetic order) associated with species A
                - Field 5: Wyckoff letters (in alphabetic order) associated with species B
                - Field 6: Wyckoff letters (in alphabetic order) associated with species C
                -  ...
                -  ...
                -  ...

              Note: for space group 47, we denote the 27th Wyckoff position, "alpha", as a capital "A"
              e.g., A_oP8_47_A

  A:B:C...:   The prototype can be decorated with different atoms by appending the elements
              to the label separated by colons ":".

  parameters: The parameters are necessary to create the structure prototype.
              The parameter information is separated by commas. The parameters are
              type double, which specify the degrees of freedom in the prototype.

                - a,b/a,c/a :        Magnitudes of lattice vectors. Note, there may not be a need
                                     to specify all magnitudes, since higher symmetry prototypes
                                     have lattice vector restrictions:
		                     TRI/MCL/MCLC/ORC/ORCC/ORCF/ORCI: a,b/a,c/a
                                     TET/BCT/HEX/RHL: a,c/a
		                     CUB/FCC/BCC: a

                - alpha,beta,gamma : Angles between lattice vectors. Note, there may not be a need
                                     to specify all angles, since higher symmetry prototypes
                                     have lattice vector restrictions:
                                     TRI: alpha,beta,gamma
                                     MCL/MCLC: beta
	                             Rest: None

                - x1,y1,z1,x2,... :  Coordinates of Wyckoff positions, which have a degree of
                                     variability. Caution: Order these parameters with care.
                                     The order of the Wyckoff positions is found in the paper.

  [--hex]:    For rhombohedral systems, there are two choices for representing the protoypes:
              rhombohedral (RHL) or hexagonal (HEX).
	              - RHL cell: a=b=c and alpha=beta=gamma (better for ab-initio, smaller than HEX)
	              - HEX cell: a=b and gamma=120 (this cell is 3 times larger than the RHL)
              The default representation is RHL. To show the prototype in the HEX representation,
              we add the --hex flag.
              NOTE: Only applicable for rhombohedral systems.

  [--vasp|--aims|--qe|--abinit|--elk|--wyccar]: The structure can be represented in a variety of ab-initio formats,
                                                e.g., VASP, FHI-AIMS, Quantum Espresso, ABINIT, ELK, and WYCKAR.
                                                The POSCAR (VASP) format is the default.
	
  [--add_equations|--equations_only]:     The symbolic math representations (i.e., equations) of the prototype can also be printed.

                                            - [--add_equations]  : appends the symbolic math representation to the output

                                            - [--equations_only] : only prints the symbolic math represention; does not require
                                                                   parameters to be specified (i.e., can neglect --params=...)

                                          Currently, the symbolic math representation is only supported for the VASP and AIMS formats.
                                          Other formats will be implemented in the future.


For additional information contact: David Hicks (david.hicks@duke.edu)

AFLOW/ANRL Prototypes
****************************************************************************************************
