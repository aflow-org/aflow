********************************************************************************
*                                                                              *
*             Aflow STEFANO CURTAROLO - Duke University 2003-2021              *
*                High-Throughput ab-initio Materials Discovery                 *
*                                                                              *
********************************************************************************
LATEST VERSION OF THE FILE:        materials.duke.edu/AFLOW/README_AFLOW_CCE.TXT
********************************************************************************

Written by Rico Friedrich, Corey Oses, and Marco Esters, 2018-2021

USER INSTRUCTIONS:

(i) GENERAL INFORMATION:
Implementation to obtain corrected DFT formation enthalpies based on the coordination corrected 
enthalpies (CCE) methodology described in:
Friedrich et al., Coordination corrected ab initio formation enthalpies, npj Comput. Mater. 5, 59 (2019);
https://doi.org/10.1038/s41524-019-0192-1
and
Friedrich et al., Automated coordination corrected enthalpies with AFLOW-CCE, Phys. Rev. Mater. 5, 043803 (2021)
https://doi.org/10.1103/PhysRevMaterials.5.043803
Please cite these articles when using this method and/or this implementation.
The corrections depend on the number of cation-anion bonds and on the cation oxidation state.
The current implementation requires a structure file as input and will return the CCE corrections for
the DFT formation enthalpies for PBE, LDA and SCAN if the oxidation numbers can be determined.
The oxidation numbers can be provided by the user or can be automatically determined from the
structure using Allen electronegativities.
In addition, the user can provide precalculated DFT formation enthalpies for PBE, LDA and SCAN and the
implementation will directly calculate and return the CCE formation enthalpies at both 298.15 and 0 K.
If DFT formation enthalpies are provided as additional input, the functionals corresponding to those
enthalpies are required. They need to be given in the same sequence as the precalculated DFT formation
enthalpies they correspond to.
However, when you specify only one formation enthalpy, no functional is required as it is assumed to
be PBE (default).
Details for the options available are explained in the following.

(ii) AVAILABLE OPTIONS:
--cce                            Prints user instructions.

--cce=STRUCTURE_FILE_PATH        Provide the path to the structure file. It can be in any structure
                                 format that AFLOW supports, e.g. VASP POSCAR, QE, AIMS, ABINIT, ELK, and CIF.
                                 For VASP, a VASP5 POSCAR is required or if a VASP4 POSCAR is used, the species
                                 must be written on the right side next to the coordinates for each atom
                                 just as for the EXAMPLE INPUT STRUCTURE FOR ROCKSALT MgO below.

--enthalpies_formation_dft=|--dfte=
                                 Provide a comma separated list of precalculated DFT formation enthalpies,
                                 they are assumed to be: (i) negative for compounds lower in enthalpy
                                 than the elements, (ii) in eV/cell. Currently, corrections are available
                                 for PBE, LDA and SCAN.

--functionals=|--func=|--functional=
                                 Provide a comma separated list of functionals for which corrections
                                 should be returned. If used together with --enthalpies_formation_dft,
                                 the functionals must be in the same sequence as the DFT formation
                                 enthalpies they correspond to. Available functionals are:
                                 (i) PBE, (ii) LDA or (iii) SCAN. Default: PBE (if only one DFT formation
                                 enthalpy is provided).

--oxidation_numbers=|--ox_nums=|--oxidation_number=
                                 Provide as a comma separated list the oxidation numbers. It is
                                 assumed that: (i) one is provided for each atom of the structure and
                                 (ii) they are in the same sequence as the corresponding atoms in the
                                 provided structure file.

--get_cce_correction|--get_cce_cor|--poscar2cce < STRUCTURE_FILE_PATH
                                 Determines the CCE corrections for the structure in STRUCTURE_FILE_PATH.
                                 It can be in any structure format that AFLOW supports, e.g. VASP POSCAR,
                                 QE, AIMS, ABINIT, ELK, and CIF. For VASP, a VASP5 POSCAR is required or if a
                                 VASP4 POSCAR is used, the species must be written on the right side next to
                                 the coordinates for each atom just as for the EXAMPLE INPUT STRUCTURE FOR
                                 ROCKSALT MgO below.

--get_oxidation_numbers|--get_ox_nums|--poscar2ox_nums < STRUCTURE_FILE_PATH
                                 Determines the oxidation numbers for the structure in STRUCTURE_FILE_PATH.
                                 It can be in any structure format that AFLOW supports, e.g. VASP POSCAR,
                                 QE, AIMS, ABINIT, ELK, and CIF. For VASP, a VASP5 POSCAR is required or if a
                                 VASP4 POSCAR is used, the species must be written on the right side next to
                                 the coordinates for each atom just as for the EXAMPLE INPUT STRUCTURE FOR
                                 ROCKSALT MgO below.

--get_cation_coordination_numbers|--get_cation_coord_nums|--poscar2cation_coord_nums < STRUCTURE_FILE_PATH
                                 Determines the number of anion neighbors for each cation for the structure
                                 in STRUCTURE_FILE_PATH. It can be in any structure format that AFLOW
                                 supports, e.g. VASP POSCAR, QE, AIMS, ABINIT, ELK, and CIF. For VASP, a
                                 VASP5 POSCAR is required or if a VASP4 POSCAR is used, the species must
                                 be written on the right side next to the coordinates for each atom just
                                 as for the EXAMPLE INPUT STRUCTURE FOR ROCKSALT MgO below.

--print=                         Obtain output in standard format (--print=out), json (--print=json), or as a python script (--print=python).
                                 Default: out.

(iii) EXAMPLE INPUT STRUCTURE FOR ROCKSALT MgO:
Mg1O1   [FCC,FCC,cF8] (STD_PRIM doi:10.1  [FCC,FCC,cF8] (STD_PRIM doi:10.1016/j.commatsci.2010.05.010)
1.224745
   0.00000000000000   1.73568248770103   1.73568248770103
   1.73568248770103   0.00000000000000   1.73568248770103
   1.73568248770103   1.73568248770103   0.00000000000000
Mg O
1 1
Direct(2) [A1B1]
   0.00000000000000   0.00000000000000   0.00000000000000  Mg
   0.50000000000000   0.50000000000000   0.50000000000000  O

(iv) EXAMPLE COMMANDS:
Assuming that AFLOW is in your PATH and you saved the above example structure file for MgO 
in the current directory as POSCAR, the following commands can be executed:

aflow --cce=POSCAR --enthalpies_formation_dft=-5.434,-6.220,-6.249 --functionals=PBE,LDA,SCAN
This will give you the CCE corrections and CCE formation enthalpies for PBE, LDA, and SCAN for MgO.

aflow --cce=POSCAR --enthalpies_formation_dft=-6.220 --functionals=LDA
This gives you only the CCE corrections and CCE formation enthalpies for LDA.

aflow --cce=POSCAR --enthalpies_formation_dft=-5.434
This gives you the CCE corrections and CCE formation enthalpies for PBE with a warning that 
PBE is assumed as functional.

aflow --cce=POSCAR
This gives you the CCE corrections for PBE, LDA, and SCAN and a rough guess of the formation 
enthalpy based on experimental formation enthalpies per bond.

aflow --cce=POSCAR --oxidation_numbers=2,-2
Oxidation numbers for each atom can also be provided as input.

(v) AVAILABLE CORRECTIONS:

Currently, corrections are available for:
OXIDES
Ag +1 from Ag2O           
Al +3 from Al2O3          
As +5 from As2O5          
B  +3 from B2O3            
Ba +2 from BaO            
Be +2 from BeO            
Bi +3 from Bi2O3          
Ca +2 from CaO            
Cd +2 from CdO            
Co +2 from CoO            
Cr +3 from Cr2O3          
Cr +6 from CrO3           
Cs +1 from Cs2O           
Cu +1 from Cu2O
Cu +2 from CuO
Fe +2 from FeO
Fe +3 from Fe2O3          
Ga +3 from Ga2O3          
Ge +4 from GeO2           
Hf +4 from HfO2
Hg +2 from HgO            
In +3 from In2O3          
Ir +4 from IrO2
K  +1 from K2O             
Li +1 from Li2O           
Mg +2 from MgO            
Mn +2 from MnO            
Mn +4 from MnO2           
Mo +4 from MoO2           
Mo +6 from MoO3           
Na +1 from Na2O
Nb +2 from NbO            
Ni +2 from NiO            
Os +4 from OsO2           
Os +8 from OsO4           
Pb +2 from PbO            
Pb +4 from PbO2           
Pd +2 from PdO            
Rb +1 from Rb2O           
Re +4 from ReO2
Re +6 from ReO3
Rh +3 from Rh2O3          
Ru +4 from RuO2           
Sb +3 from Sb2O3          
Sb +5 from Sb2O5
Sc +3 from Sc2O3          
Se +4 from SeO2           
Si +4 from SiO2(al-quartz)
Sn +2 from SnO
Sn +4 from SnO2
Sr +2 from SrO            
Te +4 from TeO2
Ti +2 from TiO            
Ti +3 from Ti2O3          
Ti +4 from TiO2(rutile)   
Tl +1 from Tl2O           
Tl +3 from Tl2O3          
V  +2 from VO              
V  +3 from V2O3            
V  +4 from VO2             
V  +5 from V2O5            
W  +4 from WO2             
W  +6 from WO3             
Y  +3 from Y2O3            
Zn +2 from ZnO            
Zr +4 from ZrO2

for per- and superoxides
O2 -2 from Li2O2
O2 -1 from KO2  

NITRIDES
Al +3 from AlN     
B  +3 from BN
Ca +2 from Ca3N2
Li +1 from Li3N    
Zn +2 from Zn3N2   
Be +2 from Be3N2   
Cr +3 from CrN     
Ga +3 from GaN
Hf +3 from HfN     
In +3 from InN
La +3 from LaN
Mg +2 from Mg3N2   
Nb +3 from NbN     
Sc +3 from ScN     
Si +4 from Si3N4   
Ta +3 from TaN     
Ti +3 from TiN     
V  +3 from VN       
Y  +3 from YN       
Zr +3 from ZrN
