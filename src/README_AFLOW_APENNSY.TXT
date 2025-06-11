****************************************************************************************************
*                                                                                                  *
*                        aflow - STEFANO CURTAROLO Duke University 2003-2021                       *
*                            High-Throughput ab-initio Computing Project                           *
*                                                                                                  *
****************************************************************************************************
LATEST VERSION OF THE FILE:                        materials.duke.edu/AFLOW/README_AFLOW_APENNSY.TXT
****************************************************************************************************

 APENNSY MODE (SC 2009-2014: the pennsy project)
  aflow --apennsy ... then --help OR
        Load library options are:
        --lib2    (loads LIB2 of binary alloys)	
        --alloy Element1Element1...  
        --list (gives a list of available alloys)
          One of the two options MUST be specified to enter in apennsy mode.
 MODIFIERS 
       --server=aflow.org
        Connects to aflow.org and download the list of servers containing
        thermodynamic informations. 
        Uses the aflowlib REST-API to download the data (default if apennsy
        does not find the LIB* directories in the server).
       --neglect=1,2,3,4,5,662.AB,
        Neglects structure /1/ /2/ ... /662.AB/ ... from the phase diagram
       --noclean
        In this case you need to specify the correct alloy with pseudopotential, 
        e.g. --alloy AgMo_pv, and you will receive the file with the right name.
       Output options are:
       --hull  (prepare gnuplot/matlab code for the convex hulls)
             it will run as "gnuplot code" or  "matlab -nodesktop -r file ..."
       --shull (prepare matlab code for the small convex hulls)
             it will run as "gnuplot code" or "matlab -nodesktop -r file ..."
       --matlab  (use matlab to plot)
       --gnuplot (use gnuplot to plot DEFAULT)
       --energy (prepare the latex code for the pdf analisys)
       --nolatex (remove latex stuff from output)
       --print=html   (makes html code when appropriate)
       --print=hyperlinks   (add hyperlinks to the LATEX/HTML code)
       --reference (|ref) gives the reference paper for that system
       --update (loads the alloy, makes the picture and makes the pdf)  
       --keep="filextensions_separated_by_commas"
              "tex"   keeps *.tex during the process (if appropriate)
              "eps"   keeps *.eps during the process (if appropriate)
              "dvi"   keeps *.dvi during the process (if appropriate)
              "toc"   keeps *.toc during the process (if appropriate)
              "GPL or gnuplot" keeps gnuplot code during the process  (if appropriate)
              "MAT or matlab"  keeps matlab code during the process  (if appropriate)
              "jpg"   keeps *.jpg during the process (if appropriate)
              "png"   keeps *.png during the process (if appropriate)
              "gif"   keeps *.gif during the process (if appropriate)

       --cite       add appropriate \\cite{} after the alloy title
       --snapshot 
         wrap up containing --update, --keep=tex, --cite, --print=hyperlinks
       --data  (prints raw data in self-explanatory form for other codes)
       --web (prints input file for the web awrapper page)
            --all --fcc --bcc --hcp
       --uncle (prints input file for uncle cluster expansion program)
             --fcc --bcc --hcp
             --Htot | --enthalpy_total
               (prints the total enthalpy of the unit cell)
             --Hat | --enthalpy_atom
               (prints the enthalpy per atom in the unit cell)
               (if no enthalpy mode is chosen, then Htot is the default).
             --Hfat | --enthalpy_formation_atom
               (prints the enthalpy per atom in the unit cell)
               (if no enthalpy mode is chosen, then Htot is the default).
      --mix | --miscibility
              (makes the aflow_nomix.cpp file for automatic miscibility
              determination, cutoff is MISCIBILIT_SYSTEM_CUTOFF
              and it is in aflow.h)

      --information 
           Prints information of each calculation in seconds 
           time(secs) cores(int) time*cores(secs) mem(MB)
 
      --experiments
           Prints the prototype of the aflow_mix_experiments.cpp file
      --miedema
           Prints the Miedema predictions
      --humerothery
           Prints the Hume-Rothery predictions
      --statistics
           Prints statistics
      --order
           Order project (values in meV)
      --rules
           Rules project for LIBRARYU
      --protocheck
	   Check the relaxed POSCARs and fix name in phase diagram
      --oss=cout 
           Prints in "cout" mode.

     Output pictures format are (through matlab code):
     --print=eps   (not necessary, it will always write eps puctures)    // XHOST.vflag
     --print=jpg   (make the matlab code to create the jpg pic)          // XHOST.vflag
     --print=pdf   (make the matlab code to create the pdf pic)          // XHOST.vflag
     --print=gif   (make the matlab code to create the gif pic)          // XHOST.vflag
     --print=png   (make the matlab code to create the png pic)          // XHOST.vflag
     --quiet (writes less on the picture, usefull for websites)

EXAMPLE
If you have gnuplot and wget available, this command should produce a pdf:
aflow --alloy AgCd --hull --energy

****************************************************************************************************
*                                                                                                  *
*                        aflow - STEFANO CURTAROLO Duke University 2003-2021                       *
*                            High-Throughput ab-initio Computing Project                           *
*                                                                                                  *
****************************************************************************************************
