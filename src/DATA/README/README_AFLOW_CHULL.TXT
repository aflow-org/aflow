****************************************************************************************************
*                                                                                                  *
*                        aflow - STEFANO CURTAROLO Duke University 2003-2021                       *
*                            High-Throughput ab-initio Computing Project                           *
*                                                                                                  *
****************************************************************************************************
LATEST VERSION OF THE FILE:                       materials.duke.edu/AFLOW/README_AFLOW_CHULL.TXT
****************************************************************************************************

AFLOW-CHULL README:

README written by: Corey Oses (corey.oses@duke.edu)

Citation info: C. Oses, E. Gossett, D. Hicks, F. Rose, E. Perim, I. Takeuchi, S. Sanvito, 
               O. Levy, C. Toher, and S. Curtarolo,
               AFLOW-CHULL: Cloud-based platform for autonomous phase stability analysis,
               J. Chem. Inf. Model. 58(12), 2477-2490 (2018)


****************************************************************************************************
GENERAL OVERVIEW

AFLOW-CHULL is an automated thermodynamic stability analysis tool. Compound data is downloaded from the
AFLOW REST-API/AFLUX Search-API and the minimum energy surface is calculated using a modified
variant (dimensionally-iterative) of the qhull convex hull algorithm:
ACM Trans. on Mathematical Software, 22(4):469-483, Dec 1996.
Various thermodynamic properties are extracted, including off-hull (full decomposition reaction, 
distances to the hull) and on-hull (equilibrium phases, equivalent ground state structures) features.
Determination of the stability criterion is also available, which has served as a powerful
screening tool for the discovery of new permanent magnets: Science Adv. 3(4), e1602241 (2017).

Other important analyses running under the hood include:
 - structure comparison analysis via Burzlaff criteria (--readme=compare)
 - outlier detection and removal (interquartile range).

Various output formats are provided, including plain text, json, PDF, and Jupyter.

****************************************************************************************************
LIST OF AFLOW-CHULL COMMANDS:
  
  aflow --convex_hull=|--chull --alloy=MnPdPt[,AlCuZn,...] [--np=1] [chull_options] [--destination=DIRECTORY]
      Queries the AFLOW API for relevant entries (see --load_library), calculates
      the convex hull, and returns the information as a PDF (default, see --output).
      --chull : Necessary flag for entering mode for calculating convex hull.
      --alloy : Necessary argument, specifies system. This code is not dimension specific, 
          i.e., you can calculate the convex hull for any n-ary system.
          There are two input modes:  raw (comma-separated) and combinatorial (colon- and comma-separated).
          Raw input:  --alloy=MnPdPt,AlCuZn.
          Combinatorial input:  --alloy=Ag,Au:Mn.  This is interpreted as --alloy=AgMn,AuMn.
      [--np=1] : Number of processors for multiple alloy inputs.
      [--destination=DIRECTORY] : Optional argument, specify the directory for the output.  
          Default is "./".
      [--usage] : Returns usage commands and options.
      [--print=|--p=|--output=|--o=latex|pdf|png|json|text|jupyter|jupyter2|jupyter3] : Select the output format. Latex/PDF are the same (.pdf).
          JSON and plain text have the following extensions:  .json and .txt.
          jupyter|jupyter2|jupyter3 creates a jupter notebook json file that plots a convex hull for the specified alloy (jupyter3 is default).
          Default is Latex/PDF.
      [--screen_only] : Output is direct to screen (stdout) instead of writing to a file. 
          All logging output is suppressed.
      [--keep=tex|--keep_tex|--keeptex|--tex] : Latex/PDF output mode only. Will keep latex .tex and put it in
          --destination=DIRECTORY.
      [--keep=log|--keep_log|--keeplog|--log] : Prints a log file of relevant output.
      [--keep=tex,log] : Keep both.
      
      LOADING OPTIONS:
      [--load_library=|--loadlibrary=|--ll=icsd|lib1|lib2|lib3] : Specify libraries
          from which to load.  Default:  icsd, lib2, and lib3.
      [--load_API|--load_api|--loadapi|--lapi|--api] : Force loading entries from the API (default unless
          on nietzsche, aflowlib, or habana).
      [--load_entries_entry_output|--loadentriesentryoutput|--leo] : Get full output for
          all entries loaded from the AFLOW API.
      [--neglect=|--ban=aflow:bb0d45ab555bc208,aflow:fb9eaa58604ce774] :
          Ban specific entries from the convex hull calculation, done by AUID.
      [--see_neglect|--seeneglect|--sn] : Show why entries were neglected.
      [--remove_extreme_points=|--removeextremepoints=|--remove_extrema=|--removeextrema=|--rep=-1000] :
          Exclude points based on Hf/Ts (floor/ceiling). Units are meV/atom / K.
      [--entropic_temperature|--entropictemperature|--entroptemp] : Calculate the Ts
          convex hull (upper-half). Default is Hf hull (lower-half).  
      [--include_paw_gga|--paw_gga] : Include calculations run with PAW_GGA. Easily coupled
          with default PAW_PBE (same level of theory, differently parametrized).

      ANALYSIS OPTIONS:
      [--distance_to_hull=|--distancetohull=|--distance2hull=|--dist2hull=|--d2h=aflow:bb0d45ab555bc208,aflow:fb9eaa58604ce774] :
          Calculates a structure's distance below/above (Hf/Ts) the hull.
      [--stability_criterion=|--stabilitycriterion=|--stable_criterion=|--scriterion=|--sc=aflow:bb0d45ab555bc208,aflow:fb9eaa58604ce774] :
          Calculates the stability criterion of the ground state structure. Will return
          a warning if structure is not a ground state one. It removes the point from the 
          hull, calculates the new hull, and calculates the distance of this point from
          below/above (Hf/Ts) hull.
      [--n+1_enthalpy_gain=|--n+1_energy_gain=|--n+1enthalpygain=|--n+1energygain=|--n+1egain=|--n1egain=|--n+1_enthalpygain=|--n+1+energygain=|--n+1_egain=|--nplus1=aflow:bb0d45ab555bc208,aflow:fb9eaa58604ce774] :
          Calculates the N+1 enthalpy gain of the ground state structure. Will return
          a warning if structure is not a ground state one. It removes all points having the 
          same dimensionality from the hull, calculates the new hull, and calculates the 
          distance of this point from below/above (Hf/Ts) hull.
      [--hull_formation_enthalpy=|--hull_enthalpy=|--hull_energy=0.25,0.25] : Returns the value of the convex hull surface at the specified coordinate/concentration.
          For stoichiometric hulls, provide the reduced composition form, i.e., for the compound Mn2PdPt, use
          --hull_formation_enthalpy=0.5,0.25, where the composition of the last component is implicitly 1-sum(0.5+0.25).
          Providing an additional component (e.g., --dist2hull=0.5,0.25,0.1) results in a rigid shift in the energy axis and,
          thus, a constant added to the final distance calculated (0.1 is in units of eV for formation energy hulls).
      [--skip_structure_comparison|--skipstructruecomparison|--skipstructcomp|--ssc] : Structure comparison analysis for determination of
          equivalent ground state structures is skipped (can be time consuming).
      [--skip_stability_criterion_analysis|--skipstabilitycriterionanalysis|--skip_scriterion|--skipscriterion|--sscriterion] : Analysis
          of stability criterion for ground state structures is skipped (can be time consuming).
      [--skip_n_plus_1_enthalpy_gain_analysis|--skip_n_plus_1_energy_gain_analysis|--skipnplus1enthalpygainanalysis|--skipnplus1energygainanalysis|--skip_nplus1|--skipnplus1|--snp1|--snpo] : Analysis
          of N+1 enthalpy gain for ground state structures is skipped (can be time consuming).
      [--include_skewed_hulls|--include_skewed|--ish] : Process hull despite skewed endpoints (abs(ground state endpoint energy) < 15 meV).
      [--include_unreliable_hulls|--include_unreliable|--iuh] : Process hull despite poor statistics (binary hull count < 200).
      [--include_outliers|--io] : Include any detected outliers (output/warning is still shown).
      [--strict_outlier_analysis|--soa] : Terminate hull analysis if binary hull statistics not satisfied.
      [--force] : Force output despite warnings.

      LATEX/PDF/PNG OPTIONS:
      [--image_only|--imageonly|--image] : Latex/PDF output mode only. Similar to --document_only, 
          but the image dimensions are not necessarily for a standard page. Preferred option for 
          importing into papers/presentations.
      [--no_document|--nodocument|--no_doc|--nodoc|--full_page_image|--fullpageimage] :
          Latex/PDF output mode only. Generates a PDF of just the convex hull illustration (2/3D).
          Dimensions: 8.5x11 inches (standard page), landscape.
      [--document_only|--documentonly|--doc_only|--doconly|--doc] : Latex/PDF output mode only.
          Exclude convex hull illustration. This is the default for quaternary systems and above.
      [--latex_output|--latexoutput] : Latex/PDF output mode only. See full latex output.
          Good for troubleshooting.
      [--latex_interactive|--latexinteractive] : Latex/PDF output mode only. Allows you
          to interact with latex as it compiles. Good for troubleshooting.
      [--light_contrast|--lightcontrast|--lc] : Latex/PDF output mode only. Modifies the convex 
          hull illustration color scheme to be lighter.
      [--large_font|--largefont|--large|--lf] : Latex/PDF output mode only. Prints a "larger" font size 
          for convex hull illustration.
      [--png_resolution=|--pngresolution=|--pngr=300] : Set PNG resolution (-density in convert).
      [--plot_iso_max_latent_heat|--plot_isomax|--iso_max|--isomax] : Plot iso-max-latent-heat lines for ground state structures.
          See dx.doi.org/10.1063/1.4902865.

For additional information contact: Corey Oses (corey.oses@duke.edu)

****************************************************************************************************
*                                                                                                  *
*                        aflow - STEFANO CURTAROLO Duke University 2003-2021                       *
*                            High-Throughput ab-initio Computing Project                           *
*                                                                                                  *
****************************************************************************************************
