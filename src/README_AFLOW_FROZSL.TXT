****************************************************************************************************
*                                                                                                  *
*                        aflow - STEFANO CURTAROLO Duke University 2003-2021                       *
*                            High-Throughput ab-initio Computing Project                           *
*                                                                                                  *
****************************************************************************************************
LATEST VERSION OF THE FILE:                         materials.duke.edu/AFLOW/README_AFLOW_FROZSL.TXT
****************************************************************************************************
Simple calculation with the Frozen Phonon code called FROZSL by Harold Stokes and
the high-throughpout wrapper functions of aflow.

1) Create an aflow.in of the type:

[AFLOW] **************************************************
[AFLOW] SYSTEM = AgZr
[AFLOW] **************************************************
[AFLOW_MODE=VASP]
[VASP] ***************************************************
[AFLOW_MODE_ZIP=xz]
[AFLOW_MODE_BINARY=vasp46s]
[AFLOW] **************************************************
[AFLOW] **************************************************
[AFLOW] **************************************************
#[VASP_RUN]GENERATE]
[VASP_RUN]STATIC
[VASP_FORCE_OPTION]KPOINTS=KEEPK
[VASP_FORCE_OPTION]WAVECAR=OFF
[VASP_FORCE_OPTION]CHGCAR=OFF
[VASP_FORCE_OPTION]SYM=ON
#[VASP_FORCE_OPTION]SYM=OFF
[VASP_FORCE_OPTION]AUTO_PSEUDOPOTENTIALS=potpaw_PBE
[VASP_FORCE_OPTION]NBANDS
#[VASP_FORCE_OPTION]SPIN=ON
[VASP_FORCE_OPTION]SPIN=OFF
#[VASP_FORCE_OPTION]LSCOUPLING=ON
#[VASP_FORCE_OPTION]AUTO_MAGMOM=ON
[VASP_FORCE_OPTION]PREC=ACCURATE
#[VASP_FORCE_OPTION]PREC=HIGH
#[VASP_FORCE_OPTION]PREC=MEDIUM
#[VASP_FORCE_OPTION]PREC=NORMAL
[VASP_FORCE_OPTION]ALGO=FAST
#[VASP_FORCE_OPTION]ALGO=VERYFAST
#[VASP_FORCE_OPTION]ALGO=NORMAL
[VASP_FORCE_OPTION]RELAX
#[VASP_FORCE_OPTION]NOTUNE
[VASP_FORCE_OPTION]TYPE=METAL
#[VASP_FORCE_OPTION]TYPE=INSULATOR
#[VASP_FORCE_OPTION]TYPE=SEMICONDUCTOR
[AFLOW] **************************************************
[VASP_INCAR_MODE_EXPLICIT]START
SYSTEM = AgZr.B11
#NELM = 120      # for bands
#NELMIN=2        # for bands
#LPLANE=.TRUE.   # for bands
#LREAL=.FALSE.   # for bands
#LSCALU=.FALSE.  # for bands
#PSTRESS=000     # for hand modification
#NBANDS=XX       # for hand modification
#IALGO=48        # for hand modification
[VASP_INCAR_MODE_EXPLICIT]STOP
[AFLOW] **************************************************  
[VASP_KPOINTS_MODE_IMPLICIT]
[VASP_KPOINTS_FILE]KSCHEME=G              // put G if cF/hP
[VASP_KPOINTS_FILE]KPPRA=8192             // KKPRA for the relax calcs
[VASP_KPOINTS_FILE]STATIC_KSCHEME=G       // put G if cF/hP
[VASP_KPOINTS_FILE]STATIC_KPPRA=8192      // KPPRA for static calcs
[AFLOW] **************************************************  
[AFLOW] **************************************************
[VASP_POTCAR_MODE_IMPLICIT]
[VASP_POTCAR_FILE]Ag
[VASP_POTCAR_FILE]Zr
[AFLOW] **************************************************
[AFLOW] no need of POSCARs as it is calculated automatically
[VASP_POSCAR_MODE_EXPLICIT]START
[VASP_POSCAR_MODE_EXPLICIT]STOP
[AFLOW] **************************************************
[AFLOW] put here the input file that you would give to the web
[AFLOW] you can use aflow --edata and --platonSG to get the Bohrs and the SG#
[FROZSL_MODE_EXPLICIT]START.FROZSL_STRUCTURE
B11 (gamma-CuTi) Phonons in AgZr on the delta line
   129                     ! P4/nmm - D_{4h}7 #129 (c2)
5.98773760769 5.98773760769 11.29335160415 90.0 90.0 90.0  ! in Bohr
   1                ! # Displacements
   0.1000           ! Displacement
   1                ! Polynomial order
   2
   2
Ag   107.87    c  0.25 0.25 0.615323036
Zr   91.224    c  0.25 0.25 0.142462423
  1                 ! Number of modes
DT 1/4              ! Put your K point phonons
  0
[FROZSL_MODE_EXPLICIT]STOP.FROZSL_STRUCTURE
[AFLOW] **************************************************
[AFLOW] if you want the Born effective charges and the dielectric corrections,
[AFLOW] put them here in FROZSL format
#[FROZSL_MODE_EXPLICIT]START.FROZSL_DIELECTRIC
#[FROZSL_MODE_EXPLICIT]STOP.FROZSL_DIELECTRIC
[AFLOW] **************************************************
[AFLOW_FROZSL]CALC


2) Once the appropriate aflow.in file is setup, you run aflow.
This will analyze your aflow.in and add frozsl output within lines
[AFLOW_FROZSL]CALC.START
...
[AFLOW_FROZSL]CALC.STOP
VASP poscars within
[VASP_POSCAR_MODE_EXPLICIT]START.mode*_disp*
....
[VASP_POSCAR_MODE_EXPLICIT]STOP.mode*_disp*
and produce the various subdirectories with the ARUN.mode*_disp*/aflow.in.
You can look all the process inside aflow.in below the line
[AFLOW_FROZSL]CALC

If something goes wrong, or if you you change your mind about phonons paths and/or parameters,
you have to delete ALL the text below [AFLOW_FROZSL]CALC (but keep this line !!),
fix your aflow.in and repeat step 2.

3) let your aflow daemon run your new aflow.ins representing displacements associated to the
irreducible representations one by one in sub-directories labeled as ARUN.mode*_disp*.
Once in a while you should check if there is an OUTCAR that has not been touched for a while. This indicates incomplete termination of VASP.
If the process gets stuck somewhere, delete the directory of the OUTCAR which is not been finished,
then delete the LOCK file of aflow, and restart aflow. The completed directories will not be modified.

4) when all the directories have finished, run
aflow --frozsl_output
This command will analyze your aflow.in, scan through all the subdirectories, warn you if there are missing calculations.
The eigen-energies and eigen-modes of the phonons will be printed on the screen and saved in the eigen.out file.

Good luck.

Stefano & Kesong




********************************************************************************


Written by Dr. Kesong Yang (kesong.yang@gmail.com) & Prof. Stefano Curtarolo (stefano@duke.edu).
09-12-2011

This file describes how to make FINDSYM/FROZSL work within gfortran compiler.

----------------------------------------FINDSYM----------------------------------------
Change all the function "dcosd" and "dsind" to "dcos" and "dsin", respectively, in the files
lattparam_to_cart.f, and wavelength_to_rgb.f, and multiply it with "2.0*pi/360.0D0".
Examples:
basis(1,2)=abc(2)*cos(abc(6))
basis(2,2)=abc(2)*sin(abc(6))

is changed to:
basis(1,2)=abc(2)*dcos(abc(6)*2.0*pi/360.0D0)
basis(2,2)=abc(2)*dsin(abc(6)*2.0*pi/360.0D0)
----------------------------------------FINDSYM----------------------------------------

----------------------------------------FORZSL----------------------------------------
1. Change all the function "dacosd" as "dacos", and multiply it with "360.0D0/2.0/pi" in the file "frozsl_init.f".
Make sure to include the constant value of pi by adding one statement "include 'dpi_real.f'" after
the "implicit none " line in the beginning of the file.

Examples:
alpha2=dacosd(dot(2,3)/b2/c2)
beta2=dacosd(dot(3,1)/c2/a2)
gamma2=dacosd(dot(1,2)/a2/b2)

are changed to:
alpha2=dacos(dot(2,3)/b2/c2)*360.0D0/2.0/pi
beta2=dacos(dot(3,1)/c2/a2)*360.0D0/2.0/pi
gamma2=dacos(dot(1,2)/a2/b2)*360.0D0/2.0/pi

2. Comment the assignment statement "dwyxyz2=0" at 1133 line of the "frozsl_init.f" file!
This is not necessary and not supported in gfortran!

3. Change the statement "CALL CARPOV(U,W,Y,KKX,NV,N,maxeig)" to "CALL CARPOV(U,W,Y,dble(KKX),NV,N,maxeig)" at line 203 of the file "hfroz.f".

Finally, note that you need to add -llapack package in the linux OS.

----------------------------------------FORZSL----------------------------------------

****************************************************************************************************
*                                                                                                  *
*                        aflow - STEFANO CURTAROLO Duke University 2003-2021                       *
*                            High-Throughput ab-initio Computing Project                           *
*                                                                                                  *
****************************************************************************************************
