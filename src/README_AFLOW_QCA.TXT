****************************************************************************************************
*                                                                                                  *
*                        aflow - STEFANO CURTAROLO Duke University 2003-2022                       *
*                            High-Throughput ab-initio Computing Project                           *
*                                                                                                  *
****************************************************************************************************
LATEST VERSION OF THE FILE:                       materials.duke.edu/AFLOW/README_AFLOW_QCA.TXT
****************************************************************************************************

AFLOW-QCA README:

README written by: Simon Divilov (simon.divilov@duke.edu)

Citation info: Lederer, Toher, Vecchio, Curtarolo, Acta Materialia 159, 364 (2018)

****************************************************************************************************
GENERAL OVERVIEW

AFLOW-QCA is an automated tool to calculate phase equilibria properties of alloys for a given 
parent lattice. Alloy data is downloaded from the AFLOW database and cluster expansion is 
performed to calculate the energies of additional derivative structures using the Alloy Theoretic 
Automated Toolkit [J. Phase Equilib. 23, 348 (2002)].

Using the generalized quasi-chemical approximation, the equi-concentration transition temperature 
and the binodal curve are calculated, based on the LTVC method [Acta Mater. 159, 364 (2018)].

****************************************************************************************************
LIST OF AFLOW-QCA COMMANDS:

  aflow --quasi_chem_approx|--qca --plattice=|--plat=fcc --elements=|--elem=Au,Pt[,Zn] [qca_options] [--directory=[DIRECTORY]]
        Loads the alloy data from the AFLOW database, calculates the equi-concentration transition
        temperature, and returns the output as a TXT (default, see --format_data).
        --quasi_chem_approx|--qca : Necessary argument, enters the mode to use the QCA.
        --plattice|--plat : Necessary argument, specifies the parent lattice of the alloy. Argument is restricted to
            only FCC (face centered cubic), BCC (body centered cubic) and HCP (hexagonal closed pack for only binaries).
        --elements|--elem : Necessary argument, specifies the elements present in the alloy.
        [--directory=[DIRECTORY]] : Optional argument, specifies the directory where to run the
            calculation.
            Default: "./"

        GENERAL OPTIONS:
        [--usage] : Returns the usage commands and options.
            Default: OFF
        [--screen_only] : Outputs directly to screen (stdout) instead of writing to a file.
            Default: OFF
        [--image_only|--image] : Reads only the output data and plot the results. Plotting
            is only implemented for binary alloys.
            Default: OFF
        [--aflowlib_directory=|--aflowlib_dir=] : Specifies a directory (including
            its subdirectories) from where to read additional aflowlib files.
            Default: ""
        [--print=|--p=|--output=|--o=txt|json|pdf|eps|png] : Selects the output format.
            Default: txt

        BINODAL OPTIONS:
        [--binodal] : Calculates the binodal curve for the given concentration and temperature
            profiles.
            Default: OFF
        [--use_sg] : Uses the space groups for structure comparison.
            Default: OFF
        [--aflow_max_num_atoms=4] : Specifies the maximum number of atoms in the unit cell of
            the AFLOW calculations.
            Default: 4
        [--max_num_atoms=|--mna=8] : Specifies the maximum number of atoms in the alloy.
            Default: 8
        [--cv_cutoff=|--cv_cut=0.05] : Specifies the cross validation (CV) score cut-off for
            the cluster expansion. The calculation will terminate prematurely if the CV
            score is below the cut-off.
            Default: 0.05
            Units: eV
        [--conc_curve_range=|--conc_curve=0,1,1,0] : Specifies the concentration curve
            to evaluate the equilibria properties. The format of this argument is:
            <c_1i>,<c_1f>,<c_2i>,<c_2f>,...,<c_[K]i>,<c_[K]f> where K is number
            of elements, <c_[k]i> is the initial concentration and <c_[k]f> is the 
            final concentration. The concentration is defined on [0,1].
            Default: NONE
        [--conc_npts=20] : Specifies the number of points used to determine the sub-division of 
            the concentration profile.
            Default: 20
        [--temp_range=|--temp=300,5000] : Specifies the temperature range to evaluate the
            equilibria properties. The format of this argument is: <T_i>,<T_f>
            where <T_i> and <T_f> is the initial and final temperature, respectively.
            Default: 300,2000
            Units: K
        [--temp_npts=100] : Specifies the number of points used to determine the sub-division of
            the temperature profile.
            Default: 150

For additional information contact: Simon Divilov (simon.divilov@duke.edu)

****************************************************************************************************
*                                                                                                  *
*                        aflow - STEFANO CURTAROLO Duke University 2003-2022                       *
*                            High-Throughput ab-initio Computing Project                           *
*                                                                                                  *
****************************************************************************************************
