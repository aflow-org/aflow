********************************************************************************
*                                                                              *
*             Aflow STEFANO CURTAROLO - Duke University 2003-2021              *
*                High-Throughput ab-initio Materials Discovery                 *
*                                                                              *
********************************************************************************
LATEST VERSION OF THE FILE:    materials.duke.edu/AFLOW/README_AFLOW_GFA.TXT
********************************************************************************

AFLOW-GFA README

This README was written by Denise Ford and the code was written by Denise Ford
and Eric Perim.

Citation info: - D.C. Ford, D. Hicks, C. Oses, C. Toher, and S. Curtarolo, 
               	 Metallic glasses for biodegradable implants, submitted (2019);
		             arXiv: 1902.00485.
               - E. Perim, D. Lee, Y. Liu, C. Toher, P. Gong, Y. Li, W. N.
	       	       Simmons, O. Levy, J.J. Vlassak, J. Schroers, and S. Curtarolo,
	       	       Spectral descriptors for bulk metallic glasses based on the
	       	       thermodynamics of competing crystalline phases, Nat. Commun. 7
	       	       (2016) 12315.
		 
********************************************************************************

OVERVIEW

This algorithm computes the spectrum of glass-forming ability (GFA) for an alloy
system based on the enthalpy and structural similarity of crystalline phases. 

Structures and formation enthalpies are loaded from the AFLOW repository.

Local atomic environments are calculated for each unique atom in each structure,
and the similarity between two structures is characterized by the number of
atomic environments that they have in common.  Local variations in the composition
of an alloy system are captured by taking linear combinations of pairs (binary
systems) or triplets (ternary systems) of structures at stoichiometries other than
the global stoichiometry, while requiring that the local stoichiometries balance
to the global stoichiometry.  Weights are assigned to each contribution from a
Gaussian distribution based on the distance between the local and global
stoichiometries.

The structures and enthalpies of metastable states (anything with a formation
enthalpy above the convex hull) are compared to those of the ground state
(defined by the convex hull).  The structures of metastable states are also
compared to each other. A Boltzmann factor is used to determine the enthalpy
proximity between a metastable state and the ground state.

Finally, the GFA is calculated over a grid of stoichiometries as the multiplication
of the functions describing the structural similarity and enthalpy similarity 
normalized by the sum of the weights.

AFLOW_CHULL is used to calculate the convex hull and AFLOW_COMPARE is used to
remove duplicate structures from the analysis.

********************************************************************************

USAGE:
  aflow --gfa --alloy=XX --ae_file=YY --cutoff_energy=ZZ > output.out

    XX is the case sensitive alloy system (e.g. CaCu).

    YY is the name of a file containing the atomic environments for the
	  system (e.g. All_atomic_environments_read.dat).  The name must not be
	  All_atomic_environments.dat.  The code will calculate the atomic
	  environments if no file is specified.

    ZZ is the formation enthalpy per atom in eV to use as a cutoff for
    including structures in the analysis.  Set this based on the
    expected glass transition temperature of the system.  The default is
    0.05 eV ~ 580 K.

  It is suggested to direct the standard out to a file (e.g. output.out)
	because it contains useful information about the calculation that
	you may want to keep.  
           
********************************************************************************

OUTPUT:
  standard out - contains information about the calculation, such as which
    atomic environments were read/calculated and the GFA calculation at each
	  stoichiometry (weights, number of combinations (pseudostructures), etc.)
	   
  GFA_entries.dat - contains the stoichiometries and formation enthalpies
	  of the structures used in the GFA calculation.

  GFA_xx.dat - contains the GFA at each stoichiometry on the grid for the
	  alloy system.

  All_atomic_environments.dat - contains the atomic environments for each
	  species in each structure used in the GFA calculation. The format is:

    AUID = aaaa   atom = species
	  # atoms  # environments
	  # types of vertices in each environment (one entry for each environment)
	  # occurrences of each environment (one entry for each environment)
	  # each vertex   # triangles   # squares (# types vertices x 3 entries)

	  For a description of atomic environments, see P. Villars, Factors
	  governing crystal structures, in: J.H. Westbrook, R.L. Fleisher (Eds.),
	  Crystal Structures of Intermetallic Compounds, Wiley, New York, 2000,
	  pp. 1-49.
	   
********************************************************************************
