# ***************************************************************************
# *                                                                         *
# *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
# *                                                                         *
# ***************************************************************************
#
#  Copyright 2003-2019 - Stefano Curtarolo - AFLOW.ORG consortium
#
#  This file is part of AFLOW software.
#
#  AFLOW is free software: you can redistribute it and/or modify it under
#  the terms of the GNU General Public License as published by the Free
#  Software Foundation, either version 3 of the License, or (at your option)
#  any later version.
#
#  This program is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# ***************************************************************************

VERSNUMBER=3.2.14

# MAKEFILE FOR AFLOW3
# FOR FROZSL -lmkl_lapack95_ilp64 -lmkl_core -lmkl_intel_lp64 -lmkl_sequential

#system variables
# https://stackoverflow.com/questions/714100/os-detecting-makefile
UNAME=$(shell uname)

# for aflow
CPP=g++ -std=c++17  -DAFLOW_MULTITHREADS_ENABLE #HE20220124 ensure c++ standard is consistent
#-O2 -frepo
# for sqllite
CC=gcc
ARCH=
ifeq ($(UNAME),Darwin)
CPP=clang++ -std=c++17 -DAFLOW_MULTITHREADS_ENABLE #HE20220124 ensure c++ standard is consistent
CC=clang
ARCH=-D_MACOSX_
endif

#CPP=g++-fsf-7 
#CPP=g++-fsf-4.9 -std=gnu++17
# -std=gnu++17
# -O3
#CPP=g++ -std=c++0x-4.7 -std=c++0x
# -std=gnu++0x
#CCFLAGS= -DCYGWIN -Wall -W -fno-inline
#ME20200520 - need to add assembler flags for Cygwin or aurostd.o won't compile
#[CO20200520 - OBSOLETE]WUNINITIALIAZED=-Wmaybe-uninitialized #-Wuninitialized on mac
W_AMBIG_OBJ=
ifeq ($(OS), Windows_NT)
W_AMBIG_OBJ=-Wa,-mbig-obj
endif
CCFLAGS=-Wall -W -fno-inline $(W_AMBIG_OBJ)
#ME20200514 - Allow inlining for threaded files (big performance penalty otherwise)
CCFLAGS_MT=-Wall -W -DAFLOW_MULTITHREADS_ENABLE
#[CO20191217 - put in OPTS_MT]CCFLAGS_MT=$(CCFLAGS) -O2
#CCFLAGS_GSA=$(CCFLAGS) -O1	#CO OBSOLETE 170601
INCLUDE = -I ./
#INCLUDE = -I ./ -I /opt/intel/mkl100/include

# -fno-inline -fno-implicit-inline-templates-frepo
OPTS=
#-O2
#-O1
#-O3
# -march=opteron
OPTS_MT=-O2
TODAY=-DTODAY=\"$(TIME)\"
VERS=-DAFLOW_VERSION=\"$(VERSNUMBER)\" $(TODAY)

# -DHOSTNAME=\"`hostname`\" -DYYY=\"2003-2012\"

#BIN_MARYLOU=/fslhome/fslcollab8/group/bin
BIN_MARYLOU=/fslhome/fslcollab8/bin
#BIN_RANGER=/share/home/00470/tg457283/bin
BIN=~/bin
BIN_KRAKEN=/nics/a/proj/aflow/bin
BIN_OHAD=/home/aflow/bin

EXTRA_FREE=../AFLOW3_FREE/EXTRA
ICSD_LIB=$(EXTRA_FREE)/ICSD
NIST_LIB=/common/NIST/2014
AFLOW_LIB=/common/AFLOW/LIBS
#NIST_LIB=/common/NIST/2010
PROJECT_DIRECTORY_GNDSTATE=/common/LIB2

# -static
OPTSA=-O3
# -static
DAY=$(shell date +%Y%m%d)
#TIME="$(shell date +%Y-%m-%d %H:%M:%S)
TIME="$(shell date +%Y-%m-%d)"
LIBS=-lpthread
# -static
#CO START 170614 easy integration
LIBS_AFLOW=-ldl $(LIBS)
FLAGS_STATIC=-static
LIBS_AFLOW_STATIC=$(LIBS_AFLOW) $(FLAGS_STATIC)
FLAGS_STATIC_PRE=
FLAGS_STATIC_POST=
#CO END 170614 easy integration
# -lglut -lGLU -lGL -D_HOST_=\"`hostname`\"
#AFLOW_VERSION := $(shell cat aflow.cpp | /home/auro/bin/getvalue AFLOW_VERSION)-$(TIME)
#AFLOW_VERSIONM := $(shell cat Makefile | /home/auro/bin/getvalue MAKEFILE)-$(TIME)
#DEBUG=-DXSTR_DEBUG
DEBUG=
BACKUP=/common/AFLOW/SRC/

# ME190823 - SQLite options
# See https://www.sqlite.org/compile.html for compile-time options
# The link recommends -DSQLITE_THREADSAFE=0, but that would disable
# all threading capabilities
SQLITE_OPTS=-DSQLITE_THREADSAFE=1

# Library MKL dynamic link
#LIB_MKL= -L /opt/intel/mkl100/lib/32 -lmkl_lapack -lmkl_core -lmkl_intel -lmkl_sequential
# Library MKL static link
#LIB_MKL= -L /opt/intel/mkl100/lib/32/ /opt/intel/mkl100/lib/32/libmkl_lapack.a \
# /opt/intel/mkl100/lib/32/libmkl_core.a \
# /opt/intel/mkl100/lib/32/libguide.a

#CO20200508 START - FIGURE OUT HOW TO GET RID OF THESE
#AUROSTD_CPPS = ./AUROSTD/aurostd.cpp ./AUROSTD/aurostd_main.cpp ./AUROSTD/aurostd_argv.cpp ./AUROSTD/aurostd_xoption.cpp ./AUROSTD/aurostd_boot.cpp ./AUROSTD/aurostd_crc64.cpp ./AUROSTD/aurostd_xscalar.cpp ./AUROSTD/aurostd_xcomplex.cpp ./AUROSTD/aurostd_xvector.cpp ./AUROSTD/aurostd_xmatrix.cpp ./AUROSTD/aurostd_xtensor.cpp ./AUROSTD/aurostd_xcombos.cpp ./AUROSTD/aurostd_xerror.cpp
#AUROSTD_HPPS = ./AUROSTD/aurostd.h ./AUROSTD/aurostd_argv.h ./AUROSTD/aurostd_xcomplex.h ./AUROSTD/aurostd_xmatrix.h ./AUROSTD/aurostd_xoption.h ./AUROSTD/aurostd_xrandom.h ./AUROSTD/aurostd_xscalar.h ./AUROSTD/aurostd_xtensor.h ./AUROSTD/aurostd_xvector.h ./AUROSTD/aurostd_xcombos.h AUROSTD/aurostd_xerror.h
#AFLOW_HPPS = aflow.h aflow_pflow.h aflow_agl_debye.h aflow_apennsy.h aflow_chull.h aflow_pocc.h aflow_pocc_old.h aflow_cce.h aflow_xgndstate.h aflow_xphases.h aflow_xvaspin.h APL/aflow_apl.h APL/aflow_qha_operations.h
#AFLOWLIB_HPP = aflowlib.h
#CO20200508 STOP - FIGURE OUT HOW TO GET RID OF THESE

#CO20170614
#ME20191125 -- Updated SQLite to version 3.30.1 (2019-10-11, https://www.sqlite.org/2019/sqlite-amalgamation-3300100.zip)
SQLITE_C = ./SQLITE/sqlite3.c
SQLITE_OBJ=$(SQLITE_C:.c=.o)

#README_SRC= aflow_readme.cpp

#MAKEFILES=Makefile Makefile.aflow	#CO20200508 - special makefile variable, DO NOT USE
MFS=Makefile Makefile.aflow

#for incrementing within a recipe
COUNTER_FILE=.aflow_counter

#CO20200508 - all MUST come first
#CO20200508 - .PHONY rules do NOT get 'Makefile' as dependence
.PHONY: all
all: aflow aflow_data 
	$(MAKE) check
include Makefile.aflow
aflow: $(SQLITE_OBJ) $(AFLOW_DEPS) $(MFS)
#	$(MAKE) -C ../AFLOW2_FREE/EXTRA/APL2AGR/ && cd ../../
#	cp -f aflow_apl*.o aflow_aapl*.o aflow_qha*.o ./APL/
#	cp -f aflow_anrl*.o ./ANRL/
#	cp sqlite3.o ./SQLITE/
#	$(CPP) $(LIBS_AFLOW) -o aflow *.o $(LIB_MKL) $(LIBS_AFLOW)
#	$(CPP) $(LIBS_AFLOW) -o aflow `ls *.o | grep -v aflow_data` $(LIB_MKL) $(LIBS_AFLOW)
	@echo "linking aflow*.o"
	@$(CPP) $(LIBS_AFLOW) -o aflow $(SQLITE_OBJ) $(AFLOW_OBJS) $(LIB_MKL) $(LIBS_AFLOW)
	ln -sf aflow aflowd
	ln -sf aflow apennsy
.PHONY: static_no_pdf
static_no_pdf: aflow aflow_data
#	[CO20200513 - KEEP rm command in case we do 'make -j' first]
	rm -rf $(SQLITE_OBJ) aflow aflow_data
	$(MAKE) LIBS_AFLOW="$(LIBS_AFLOW_STATIC)" FLAGS_STATIC_POST="$(FLAGS_STATIC)"
#	$(CC) -DSQLITE_ENABLE_JSON1 -c $(SQLITE_C) $(LIBS) $(FLAGS_STATIC)	#CO20190401 - $(LIBS) -lpthreads not needed, issues with clang++ (mac)
#	$(MAKE) $(SQLITE_OBJ) FLAGS_STATIC_POST="$(FLAGS_STATIC)"
#	rm -f aflow
#	$(MAKE) aflow LIBS_AFLOW="$(LIBS_AFLOW_STATIC)"
#	rm -f aflow_data
#	$(MAKE) aflow_data FLAGS_STATIC_POST="$(FLAGS_STATIC)"
.PHONY: static
static: static_no_pdf
#	tar Jcfvp /www/AFLOW/CURRENT/aflow_$(VERSNUMBER)_amd64_debian_static.tar.xz aflow aflow_data	
#ME190313 - Change compiler to clang because g++ fails to compile on some Macs

#HE20210629 build with the oldest supported c++ standard
.PHONY: min_cpp
min_cpp: CPP+= -std=c++98
min_cpp: aflow aflow_data

#HE20210318 - add debug build
#For debugging build aflow with
#  make debug
#then start it using GNU Debugger
#  gdb aflow
#see also https://www.gnu.org/software/gdb/documentation/
.PHONY: debug
debug: CPP += -ggdb
debug: aflow aflow_data

#HE20220221 - add coverage build
#
.PHONY: coverage
coverage: CPP += --coverage
coverage: aflow aflow_data

.PHONY: mac
mac:
#	[CO20200508 - OBSOLETE]$(MAKE) CC=clang CPP=clang++ aflow_data
	$(MAKE) CC=clang CPP="clang++ -std=c++17" ARCH=-D_MACOSX_ -j 4
#	[CO20200508 - OBSOLETE]ln -sf aflow aflowd
#	[CO20200508 - OBSOLETE]ln -sf aflow apennsy
#	[CO20200508 - OBSOLETE]$(MAKE) check
#.PHONY: mac_static
#mac_static: mac
#	tar Jcfvp aflow_$(VERSNUMBER)_maxOSX_10.X.X.tar.xz aflow aflow_data
.PHONY: production
production: 
	$(MAKE) OPTS="$(OPTS) -D__XOPTIMIZE" OPTS_MT="$(OPTS_MT) -D__XOPTIMIZE"


###################################################################################################
#CO201706147
$(SQLITE_OBJ): $(SQLITE_C) Makefile
	$(CC) $(SQLITE_OPTS) -c $< -o $@ $(FLAGS_STATIC_POST)
#	$(CC) -DSQLITE_ENABLE_JSON1 -DSQLITE_OMIT_LOAD_EXTENSION -o $@ -c $< $(LIBS)	#-DSQLITE_OMIT_LOAD_EXTENSION doesn't really work amalgamation
#	$(CC) -DSQLITE_ENABLE_JSON1 -o $@ -c $< $(LIBS)
#	$(CC) -DSQLITE_ENABLE_JSON1 -c $< $(LIBS)	#CO20190401 - $(LIBS) -lpthreads not needed, issues with clang++ (mac)
#	$(CC) $(SQLITE_OPTS) -c $<
#	cp sqlite3.o ./SQLITE/
#CO20200508 - catch all rule
#[CO20200508 - do not include AFLOW*.h or AUROSTD*.cpp/AUROSTD*.h].cpp.o: .cpp .h $(MFS) $(AFLOW_HPPS) $(AUROSTD_HPPS) $(AUROSTD_CPPS)
#[CO20200508 - addressing: warning: ignoring prerequisites on suffix rule definition].cpp.o: .cpp .h Makefile
%.o: %.cpp %.h Makefile
	$(CPP) $(VERS) -D_AFLOW_FILE_NAME_=\""$<"\" $(INCLUDE) $(CCFLAGS) $(OPTS) $(ARCH) $< -c -o $@
%.o: %.cpp Makefile
	$(CPP) $(VERS) -D_AFLOW_FILE_NAME_=\""$<"\" $(INCLUDE) $(CCFLAGS) $(OPTS) $(ARCH) $< -c -o $@
#	touch $@
###################################################################################################
# commands to cut/paste to generate stuff SC20200601
pseudo: 
	rm -f aflow_xpseudopotentials_data.cpp
	touch aflow_xpseudopotentials.cpp aflow_xpseudopotentials_data.cpp
	make -j
	rm -f aflow_xpseudopotentials_data.cpp
	./aflow --scrub=POTCAR --FILE `find /common/VASP/pot_LDA/ -name POTCAR |sort ` >> aflow_xpseudopotentials_data.cpp
	./aflow --scrub=POTCAR --FILE `find /common/VASP/pot_GGA/ -name POTCAR |sort ` >> aflow_xpseudopotentials_data.cpp
	./aflow --scrub=POTCAR --FILE `find /common/VASP/potpaw_LDA/ -name POTCAR |sort ` >> aflow_xpseudopotentials_data.cpp
	./aflow --scrub=POTCAR --FILE `find /common/VASP/potpaw_GGA/ -name POTCAR |sort ` >> aflow_xpseudopotentials_data.cpp
	./aflow --scrub=POTCAR --FILE `find /common/VASP/potpaw_PBE/ -name POTCAR |sort ` >> aflow_xpseudopotentials_data.cpp
	./aflow --scrub=POTCAR --FILE `find /common/VASP/potpaw_LDA.54/ -name POTCAR |sort ` >> aflow_xpseudopotentials_data.cpp
	./aflow --scrub=POTCAR --FILE `find /common/VASP/potpaw_PBE.54/ -name POTCAR |sort ` >> aflow_xpseudopotentials_data.cpp
	touch aflow_xpseudopotentials.cpp aflow_xpseudopotentials_data.cpp
	# subst "/common/VASP/" "" aflow_xpseudopotentials_data.cpp
	cat aflow_xpseudopotentials_data.cpp | grep -c TITEL
	cat aflow_xpseudopotentials_data.cpp | grep structure | grep -v N/A | wc
	make -j

###################################################################################################
.PHONY: extra
extra:	findsym frozsl apl2agr platon subst
apl2agr: Makefile
	@echo "MAKING APL2AGR"
	if [ -f $(EXTRA_FREE)/APL2AGR/Makefile ]; then $(MAKE) -C $(EXTRA_FREE)/APL2AGR realclean; fi
	if [ -f $(EXTRA_FREE)/APL2AGR/Makefile ]; then $(MAKE) -C $(EXTRA_FREE)/APL2AGR; fi
	if [ -f $(EXTRA_FREE)/APL2AGR/apl2agr ] && [ -d ../AFLOW3_NONFREE/BIN ]; then cp $(EXTRA_FREE)/APL2AGR/apl2agr ../AFLOW3_NONFREE/BIN; fi
subst: Makefile
	@echo "MAKING SUBST"
	if [ -f $(EXTRA_FREE)/SUBST/Makefile ]; then $(MAKE) -C $(EXTRA_FREE)/SUBST clean; fi
	if [ -f $(EXTRA_FREE)/SUBST/Makefile ]; then $(MAKE) -C $(EXTRA_FREE)/SUBST; fi
	if [ -f $(EXTRA_FREE)/SUBST/subst ] && [ -d ../AFLOW3_NONFREE/BIN ]; then cp $(EXTRA_FREE)/SUBST/subst ../AFLOW3_NONFREE/BIN; fi
findsym: Makefile
	@echo "MAKING FINDSYM"
	if [ -f ../AFLOW3_NONFREE/EXTRA/FINDSYM/Makefile ]; then $(MAKE) -C ../AFLOW3_NONFREE/EXTRA/FINDSYM clean; fi
	if [ -f ../AFLOW3_NONFREE/EXTRA/FINDSYM/Makefile ]; then $(MAKE) -C ../AFLOW3_NONFREE/EXTRA/FINDSYM; fi
	if [ -f ../AFLOW3_NONFREE/EXTRA/FINDSYM/findsym ]; then cp ../AFLOW3_NONFREE/EXTRA/FINDSYM/findsym ../AFLOW3_NONFREE/BIN; fi
frozsl: Makefile
	@echo "MAKING FROZSL"
	if [ -f ../AFLOW3_NONFREE/EXTRA/FROZSL/Makefile ]; then $(MAKE) -C ../AFLOW3_NONFREE/EXTRA/FROZSL clean; fi
	if [ -f ../AFLOW3_NONFREE/EXTRA/FROZSL/Makefile ]; then $(MAKE) -C ../AFLOW3_NONFREE/EXTRA/FROZSL; fi
	if [ -f ../AFLOW3_NONFREE/EXTRA/FROZSL/frozsl ]; then cp ../AFLOW3_NONFREE/EXTRA/FROZSL/frozsl ../AFLOW3_NONFREE/BIN; fi
platon: Makefile
	@echo "MAKING PLATON"
	if [ -f ../AFLOW3_NONFREE/EXTRA/PLATON/51108/Makefile ]; then $(MAKE) -C ../AFLOW3_NONFREE/EXTRA/PLATON/51108 clean; fi
	if [ -f ../AFLOW3_NONFREE/EXTRA/PLATON/51108/Makefile ]; then $(MAKE) -C ../AFLOW3_NONFREE/EXTRA/PLATON/51108; fi
	if [ -f ../AFLOW3_NONFREE/EXTRA/PLATON/51108/platon ]; then cp ../AFLOW3_NONFREE/EXTRA/PLATON/51108/platon ../AFLOW3_NONFREE/BIN; fi
enum: Makefile
	@echo "MAKING ENUM"
	if [ -f $(EXTRA_FREE)/GUS/Makefile ]; then $(MAKE) -C $(EXTRA_FREE)/GUS; fi
.PHONY: extra_clean
extra_clean:
	if [ -f ../AFLOW3_NONFREE/EXTRA/FINDSYM/Makefile ]; then $(MAKE) -C ../AFLOW3_NONFREE/EXTRA/FINDSYM clean; fi
	if [ -f ../AFLOW3_NONFREE/EXTRA/FROZSL/Makefile ]; then $(MAKE) -C ../AFLOW3_NONFREE/EXTRA/FROZSL clean; fi
	if [ -f ../AFLOW3_NONFREE/EXTRA/PLATON/51108/Makefile ]; then $(MAKE) -C ../AFLOW3_NONFREE/EXTRA/PLATON/51108 clean; fi
	if [ -f $(EXTRA_FREE)/GUS/Makefile ]; then $(MAKE) -C $(EXTRA_FREE)/GUS realclean; fi
	if [ -f $(EXTRA_FREE)/APL2AGR/Makefile ]; then $(MAKE) -C $(EXTRA_FREE)/APL2AGR realclean; fi
###################################################################################################

###################################################################################################
aflow_library_icsd.dat.xz: $(ICSD_LIB)/README_LIBRARY_ICSD1.TXT $(ICSD_LIB)/README_LIBRARY_ICSD2.TXT $(ICSD_LIB)/README_LIBRARY_ICSD3.TXT $(ICSD_LIB)/README_LIBRARY_ICSD4.TXT $(ICSD_LIB)/README_LIBRARY_ICSD5.TXT $(ICSD_LIB)/README_LIBRARY_ICSD6.TXT $(ICSD_LIB)/README_LIBRARY_ICSD7.TXT $(ICSD_LIB)/README_LIBRARY_ICSD8.TXT $(ICSD_LIB)/README_LIBRARY_ICSD9.TXT Makefile
	$(MAKE) aflow -j 16
	rm -rf aflow_library_icsd.dat
	echo "[AFLOW] aflow_library_icsd.dat " >> aflow_library_icsd.dat
	echo "[README_LIBRARY_ICSD1.TXT]START " >> aflow_library_icsd.dat
	cat $(ICSD_LIB)/README_LIBRARY_ICSD1.TXT >> aflow_library_icsd.dat
	echo "[README_LIBRARY_ICSD1.TXT]STOP " >> aflow_library_icsd.dat
	echo "[README_LIBRARY_ICSD2.TXT]START " >> aflow_library_icsd.dat
	cat $(ICSD_LIB)/README_LIBRARY_ICSD2.TXT >> aflow_library_icsd.dat
	echo "[README_LIBRARY_ICSD2.TXT]STOP " >> aflow_library_icsd.dat
	echo "[README_LIBRARY_ICSD3.TXT]START " >> aflow_library_icsd.dat
	cat $(ICSD_LIB)/README_LIBRARY_ICSD3.TXT >> aflow_library_icsd.dat
	echo "[README_LIBRARY_ICSD3.TXT]STOP " >> aflow_library_icsd.dat
	echo "[README_LIBRARY_ICSD4.TXT]START " >> aflow_library_icsd.dat
	cat $(ICSD_LIB)/README_LIBRARY_ICSD4.TXT >> aflow_library_icsd.dat
	echo "[README_LIBRARY_ICSD4.TXT]STOP " >> aflow_library_icsd.dat
	echo "[README_LIBRARY_ICSD5.TXT]START " >> aflow_library_icsd.dat
	cat $(ICSD_LIB)/README_LIBRARY_ICSD5.TXT >> aflow_library_icsd.dat
	echo "[README_LIBRARY_ICSD5.TXT]STOP " >> aflow_library_icsd.dat
	echo "[README_LIBRARY_ICSD6.TXT]START " >> aflow_library_icsd.dat
	cat $(ICSD_LIB)/README_LIBRARY_ICSD6.TXT >> aflow_library_icsd.dat
	echo "[README_LIBRARY_ICSD6.TXT]STOP " >> aflow_library_icsd.dat
	echo "[README_LIBRARY_ICSD7.TXT]START " >> aflow_library_icsd.dat
	cat $(ICSD_LIB)/README_LIBRARY_ICSD7.TXT >> aflow_library_icsd.dat
	echo "[README_LIBRARY_ICSD7.TXT]STOP " >> aflow_library_icsd.dat
	echo "[README_LIBRARY_ICSD8.TXT]START " >> aflow_library_icsd.dat
	cat $(ICSD_LIB)/README_LIBRARY_ICSD8.TXT >> aflow_library_icsd.dat
	echo "[README_LIBRARY_ICSD8.TXT]STOP " >> aflow_library_icsd.dat
	echo "[README_LIBRARY_ICSD9.TXT]START " >> aflow_library_icsd.dat
	cat $(ICSD_LIB)/README_LIBRARY_ICSD9.TXT >> aflow_library_icsd.dat
	echo "[README_LIBRARY_ICSD9.TXT]STOP " >> aflow_library_icsd.dat
	cat aflow_library_icsd.dat |sed "s/WYCK_ICSD/WYC_ICSD/g" |sed "s/  / /g"|sed "s/  / /g"|sed "s/  / /g"|sed "s/  / /g"|sed "s/  / /g" | grep -v "//" > aflow_library_icsd.dat2
	mv aflow_library_icsd.dat2 aflow_library_icsd.dat
	xz -9fv -T8 aflow_library_icsd.dat
	cp -f aflow_library_icsd.dat.xz ../AFLOW3_NONFREE/LIBS
	cp -f aflow_library_icsd.dat.xz /common/AFLOW/LIBS
	rm -f aflow_data_*
	$(MAKE) aflow_data
#	rm -f $(EXTRA_FREE)/NIST/ICSD_List_2014.txt*
	./aflow --protos_icsd > $(EXTRA_FREE)/NIST/ICSD_List_2014.txt
	rm -f $(EXTRA_FREE)/NIST/ICSD_List_2014.txt.xz
	xz -9fv -T8 $(EXTRA_FREE)/NIST/ICSD_List_2014.txt
$(ICSD_LIB)/README_LIBRARY_ICSD1.TXT: Makefile
	rm -rf $@
	echo "// ***************************************************************************">>$@
	xzcat $(NIST_LIB)/README_LIBRARY_ICSD1.TXT.xz >>$@
	echo "// ***************************************************************************">>$@
	echo " STRUCTURE ___PROTO_END___">>$@
	echo " END">>$@
	echo "// ***************************************************************************">>$@
$(ICSD_LIB)/README_LIBRARY_ICSD2.TXT: Makefile
	rm -rf $@
	echo "// ***************************************************************************">>$@
	xzcat $(NIST_LIB)/README_LIBRARY_ICSD2.TXT.xz >>$@
	echo "// ***************************************************************************">>$@
	echo " STRUCTURE ___PROTO_END___">>$@
	echo " END">>$@
	echo "// ***************************************************************************">>$@
$(ICSD_LIB)/README_LIBRARY_ICSD3.TXT: Makefile
	rm -rf $@
	echo "// ***************************************************************************">>$@
	xzcat $(NIST_LIB)/README_LIBRARY_ICSD3.TXT.xz >>$@
	echo "// ***************************************************************************">>$@
	echo " STRUCTURE ___PROTO_END___">>$@
	echo " END">>$@
	echo "// ***************************************************************************">>$@
$(ICSD_LIB)/README_LIBRARY_ICSD4.TXT: Makefile
	rm -rf $@
	echo "// ***************************************************************************">>$@
	xzcat $(NIST_LIB)/README_LIBRARY_ICSD4.TXT.xz >>$@
	echo "// ***************************************************************************">>$@
	echo " STRUCTURE ___PROTO_END___">>$@
	echo " END">>$@
	echo "// ***************************************************************************">>$@
$(ICSD_LIB)/README_LIBRARY_ICSD5.TXT: Makefile
	rm -rf $@
	echo "// ***************************************************************************">>$@
	xzcat $(NIST_LIB)/README_LIBRARY_ICSD5.TXT.xz >>$@
	echo "// ***************************************************************************">>$@
	echo " STRUCTURE ___PROTO_END___">>$@
	echo " END">>$@
	echo "// ***************************************************************************">>$@
$(ICSD_LIB)/README_LIBRARY_ICSD6.TXT: Makefile
	rm -rf $@
	echo "// ***************************************************************************">>$@
	xzcat $(NIST_LIB)/README_LIBRARY_ICSD6.TXT.xz >>$@
	echo "// ***************************************************************************">>$@
	echo " STRUCTURE ___PROTO_END___">>$@
	echo " END">>$@
	echo "// ***************************************************************************">>$@
$(ICSD_LIB)/README_LIBRARY_ICSD7.TXT: Makefile
	rm -rf $@
	echo "// ***************************************************************************">>$@
	xzcat $(NIST_LIB)/README_LIBRARY_ICSD7.TXT.xz >>$@
	echo "// ***************************************************************************">>$@
	echo " STRUCTURE ___PROTO_END___">>$@
	echo " END">>$@
	echo "// ***************************************************************************">>$@
$(ICSD_LIB)/README_LIBRARY_ICSD8.TXT: Makefile
	rm -rf $@
	echo "// ***************************************************************************">>$@
	xzcat $(NIST_LIB)/README_LIBRARY_ICSD8.TXT.xz >>$@
	echo "// ***************************************************************************">>$@
	echo " STRUCTURE ___PROTO_END___">>$@
	echo " END">>$@
	echo "// ***************************************************************************">>$@
$(ICSD_LIB)/README_LIBRARY_ICSD9.TXT: Makefile
	rm -rf $@
	echo "// ***************************************************************************">>$@
	xzcat $(NIST_LIB)/README_LIBRARY_ICSD9.TXT.xz >>$@
	echo "// ***************************************************************************">>$@
	echo " STRUCTURE ___PROTO_END___">>$@
	echo " END">>$@
	echo "// ***************************************************************************">>$@
############################################################################################################
aflow_data_calculated.cpp: Makefile
	rm -rf $@
	# vLIBS	
	if [ -f $(EXTRA_FREE)/CALCULATED/vLIBS.dat.xz ]; then echo "// $@ automatic generated from CALCULATED/vLIBS.dat" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/CALCULATED/vLIBS.dat.xz ]; then echo "#include <string>" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/CALCULATED/vLIBS.dat.xz ]; then echo "std::string vAUID,vAURL,vLOOP;" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/CALCULATED/vLIBS.dat.xz ]; then echo "std::string vLIBS=\"\\" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/CALCULATED/vLIBS.dat.xz ]; then xzcat $(EXTRA_FREE)/CALCULATED/vLIBS.dat.xz | grep -v "//" | perl -p -e 's/\n/\\n\\\n/g' >> $@ ; fi
	if [ -f $(EXTRA_FREE)/CALCULATED/vLIBS.dat.xz ]; then echo "\";" >> $@ ; fi
	# done
	-if [ ! -f $@ ]; then rm -rf $@.xz*; wget http://materials.duke.edu/AFLOW/AFLOW3_FREE/$@.xz && xz -dv $@.xz ; fi
	-if [ ! -f $@ ]; then rm -rf $@.xz*; wget http://aflowlib.duke.edu/AFLOW/AFLOW3_FREE/$@.xz && xz -dv $@.xz ; fi
	if [ -f $@ ]; then GUARD_CHECK=$$(grep '_AFLOW_DATA_CALCULATED_CPP_' $@); if [ -z "$$GUARD_CHECK" ]; then mv $@ $@_tmp; echo "#ifndef _AFLOW_DATA_CALCULATED_CPP_" >> $@; echo "#define _AFLOW_DATA_CALCULATED_CPP_" >> $@; cat $@_tmp >> $@ && rm -f $@_tmp; echo "#endif // _AFLOW_DATA_CALCULATED_CPP_" >> $@; fi; fi
	if [ ! -f $@ ]; then exit 1; fi
	touch $@
aflow_data_htqc.cpp: README_LIBRARY_HTQC.TXT README_LIBRARY_HTQC_BORIDES.TXT README_LIBRARY_HTQC_CARBIDES.TXT Makefile
	rm -f $@
	echo "#ifndef _AFLOW_DATA_HTQC_CPP_" >> $@
	echo "#define _AFLOW_DATA_HTQC_CPP_" >> $@
	echo "// $@ automatic generated from README_LIBRARY_HTQC.TXT" >> $@
	echo "#include <string>" >> $@
	echo "std::string Library_HTQC=\"\\" >> $@
	cat README_LIBRARY_HTQC.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> $@
	cat README_LIBRARY_HTQC_BORIDES.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> $@
	cat README_LIBRARY_HTQC_CARBIDES.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> $@
	echo "STRUCTURE ___PROTO_END___" | perl -p -e 's/\n/ \\n\\\n/g' >> $@
	echo " END" | perl -p -e 's/\n/ \\n\\\n/g' >> $@
	echo "\";" >> $@
	echo "#endif // _AFLOW_DATA_HTQC_CPP_" >> $@
#aflow_data_binary.cpp: Makefile
#	rm -rf aflow_data_binary.cpp
#	echo "// aflow_data_binary.cpp automatic generated from PROTO/Binary.exp" >> aflow_data_binary.cpp
#	echo "std::string AFLOW_BinaryRead=\"\\" >> aflow_data_binary.cpp
#	cat $(EXTRA_FREE)/PROTO/Binary.exp | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> aflow_data_binary.cpp
#	echo "\";" >> aflow_data_binary.cpp
#	echo "std::string AFLOW_Binary_Angle_Read=\"\\" >> aflow_data_binary.cpp
#	cat $(EXTRA_FREE)/PROTO/Binary_angle.exp | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> aflow_data_binary.cpp
#	echo "\";" >> aflow_data_binary.cpp
aflow_data_libraries.cpp: Makefile
	rm -rf $@
	echo "#ifndef _AFLOW_DATA_LIBRARIES_CPP_" >> $@
	echo "#define _AFLOW_DATA_LIBRARIES_CPP_" >> $@
	echo "// $@ automatic generated from ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_libX.dat" >> $@
	echo "#include <string>" >> $@
#	echo "std::string aflowlib_icsd=\"\\" >> $@ && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_icsd.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\\Gamma/\\\\Gamma/g' | perl -p -e 's/\\Sigma/\\\\Sigma/g' >> $@ && echo "\";" >> $@
	echo "std::string aflowlib_icsd=\"\";" >> $@
	@echo "LIB0"
#	echo "std::string aflowlib_lib0=\"\\" >> $@ && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib0.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\\Gamma/\\\\Gamma/g' | perl -p -e 's/\\Sigma/\\\\Sigma/g' >> $@ && echo "\";" >> $@	
	echo "std::string aflowlib_lib0=\"\";" >> $@
	@echo "LIB1"
#	echo "std::string aflowlib_lib1=\"\\" >> $@ && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib1.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\\Gamma/\\\\Gamma/g' | perl -p -e 's/\\Sigma/\\\\Sigma/g' >> $@ && echo "\";" >> $@	
	echo "std::string aflowlib_lib1=\"\";" >> $@
	@echo "LIB2"
#	echo "std::string aflowlib_lib2=\"\\" >> $@ && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib2.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\\Gamma/\\\\Gamma/g' | perl -p -e 's/\\Sigma/\\\\Sigma/g' >> $@ && echo "\";" >> $@	
	echo "std::string aflowlib_lib2=\"\";" >> $@
	@echo "LIB3"
#	echo "std::string aflowlib_lib3=\"\\" >> $@ && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib3.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\\Gamma/\\\\Gamma/g' | perl -p -e 's/\\Sigma/\\\\Sigma/g' >> $@ && echo "\";" >> $@
	echo "std::string aflowlib_lib3=\"\";" >> $@
	@echo "LIB4"
#	echo "std::string aflowlib_lib4=\"\\" >> $@ && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib4.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\\Gamma/\\\\Gamma/g' | perl -p -e 's/\\Sigma/\\\\Sigma/g' >> $@ && echo "\";" >> $@
	echo "std::string aflowlib_lib4=\"\";" >> $@
	@echo "LIB5"
#	echo "std::string aflowlib_lib5=\"\\" >> $@ && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib5.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\\Gamma/\\\\Gamma/g' | perl -p -e 's/\\Sigma/\\\\Sigma/g' >> $@ && echo "\";" >> $@
	echo "std::string aflowlib_lib5=\"\";" >> $@
	@echo "LIB6"
#	echo "std::string aflowlib_lib6=\"\\" >> $@ && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib6.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\\Gamma/\\\\Gamma/g' | perl -p -e 's/\\Sigma/\\\\Sigma/g' >> $@ && echo "\";" >> $@
	echo "std::string aflowlib_lib6=\"\";" >> $@
	@echo "LIB7"
#	echo "std::string aflowlib_lib7=\"\\" >> $@ && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib7.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\\Gamma/\\\\Gamma/g' | perl -p -e 's/\\Sigma/\\\\Sigma/g' >> $@ && echo "\";" >> $@
	echo "std::string aflowlib_lib7=\"\";" >> $@
	@echo "LIB8"
#	echo "std::string aflowlib_lib8=\"\\" >> $@ && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib8.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\\Gamma/\\\\Gamma/g' | perl -p -e 's/\\Sigma/\\\\Sigma/g' >> $@ && echo "\";" >> $@
	echo "std::string aflowlib_lib8=\"\";" >> $@
	@echo "LIB9"
#	echo "std::string aflowlib_lib9=\"\\" >> $@ && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib9.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\\Gamma/\\\\Gamma/g' | perl -p -e 's/\\Sigma/\\\\Sigma/g' >> $@ && echo "\";" >> $@
	echo "std::string aflowlib_lib9=\"\";" >> $@
	echo "#endif // _AFLOW_DATA_LIBRARIES_CPP_" >> $@
aflow_data_stokes.cpp: Makefile 
	rm -rf $@
	if [ -f $(EXTRA_FREE)/DATA_FINDSYM/data_space.txt.xz ]; then echo "// $@ automatic generated " >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_FINDSYM/data_space.txt.xz ]; then echo "#include <string>" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_FINDSYM/data_space.txt.xz ]; then echo "std::string FINDSYM_data_space_txt=\"\\" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_FINDSYM/data_space.txt.xz ]; then xzcat $(EXTRA_FREE)/DATA_FINDSYM/data_space.txt.xz | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\0//g' >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_FINDSYM/data_space.txt.xz ]; then echo "\";" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_FINDSYM/data_wyckoff.txt.xz ]; then echo "std::string FINDSYM_data_wyckoff_txt=\"\\" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_FINDSYM/data_wyckoff.txt.xz ]; then xzcat $(EXTRA_FREE)/DATA_FINDSYM/data_wyckoff.txt.xz | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\0//g' >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_FINDSYM/data_wyckoff.txt.xz ]; then echo "\";" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_space.txt.xz ]; then echo "std::string FROZSL_data_space_txt=\"\\" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_space.txt.xz ]; then xzcat $(EXTRA_FREE)/DATA_FROZSL/data_space.txt.xz | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\0//g' >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_space.txt.xz ]; then echo "\";" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_wyckoff.txt.xz ]; then echo "std::string FROZSL_data_wyckoff_txt=\"\\" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_wyckoff.txt.xz ]; then xzcat $(EXTRA_FREE)/DATA_FROZSL/data_wyckoff.txt.xz | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\0//g' >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_wyckoff.txt.xz ]; then echo "\";" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_images.txt.xz ]; then echo "std::string FROZSL_data_images_txt=\"\\" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_images.txt.xz ]; then xzcat $(EXTRA_FREE)/DATA_FROZSL/data_images.txt.xz | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\0//g' >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_images.txt.xz ]; then echo "\";" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_irreps.txt.xz ]; then echo "std::string FROZSL_data_irreps_txt=\"\\" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_irreps.txt.xz ]; then xzcat $(EXTRA_FREE)/DATA_FROZSL/data_irreps.txt.xz | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\0//g' >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_irreps.txt.xz ]; then echo "\";" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_isotropy.txt.xz ]; then echo "std::string FROZSL_data_isotropy_txt=\"\\" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_isotropy.txt.xz ]; then xzcat $(EXTRA_FREE)/DATA_FROZSL/data_isotropy.txt.xz | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\0//g' >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_isotropy.txt.xz ]; then echo "\";" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_little.txt.xz ]; then echo "std::string FROZSL_data_little_txt=\"\\" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_little.txt.xz ]; then xzcat $(EXTRA_FREE)/DATA_FROZSL/data_little.txt.xz | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\0//g' >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_little.txt.xz ]; then echo "\";" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/symmetry2.dat.xz ]; then echo "std::string FROZSL_symmetry2_dat=\"\\" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/symmetry2.dat.xz ]; then xzcat $(EXTRA_FREE)/DATA_FROZSL/symmetry2.dat.xz | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\0//g' >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/symmetry2.dat.xz ]; then echo "\";" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/const.dat.xz ]; then echo "std::string FROZSL_const_dat=\"\\" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/const.dat.xz ]; then xzcat $(EXTRA_FREE)/DATA_FROZSL/const.dat.xz | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\0//g' >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/const.dat.xz ]; then echo "\";" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_PHVASP/phvaspsetup_AFLOW.xz ]; then echo "std::string FROZSL_phvaspsetup_AFLOW=\"\\" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_PHVASP/phvaspsetup_AFLOW.xz ]; then xzcat $(EXTRA_FREE)/DATA_PHVASP/phvaspsetup_AFLOW.xz | perl -p -e 's/\#/_ASCII_23_/g' | perl -p -e 's/\"/_ASCII_22_/g' | perl -p -e 's/\\/_ASCII_5C_/g' | perl -p -e 's/\n/ \\n\\\n/g' >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_PHVASP/phvaspsetup_AFLOW.xz ]; then echo "\";" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_PHVASP/phvaspsetup_POSCAR.xz ]; then echo "std::string FROZSL_phvaspsetup_POSCAR=\"\\" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_PHVASP/phvaspsetup_POSCAR.xz ]; then xzcat $(EXTRA_FREE)/DATA_PHVASP/phvaspsetup_POSCAR.xz | perl -p -e 's/\#/_ASCII_23_/g' | perl -p -e 's/\"/_ASCII_22_/g' | perl -p -e 's/\\/_ASCII_5C_/g' | perl -p -e 's/\n/ \\n\\\n/g' >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_PHVASP/phvaspsetup_POSCAR.xz ]; then echo "\";" >> $@ ; fi
	-if [ ! -f $@ ]; then rm -rf $@.xz*; wget http://materials.duke.edu/AFLOW/AFLOW3_FREE/$@.xz && xz -dv $@.xz ; fi
	-if [ ! -f $@ ]; then rm -rf $@.xz*; wget http://aflowlib.duke.edu/AFLOW/AFLOW3_FREE/$@.xz && xz -dv $@.xz ; fi
	if [ -f $@ ]; then GUARD_CHECK=$$(grep '_AFLOW_DATA_STOKES_CPP_' $@); if [ -z "$$GUARD_CHECK" ]; then mv $@ $@_tmp; echo "#ifndef _AFLOW_DATA_STOKES_CPP_" >> $@; echo "#define _AFLOW_DATA_STOKES_CPP_" >> $@; cat $@_tmp >> $@ && rm -f $@_tmp; echo "#endif // _AFLOW_DATA_STOKES_CPP_" >> $@; fi; fi
	if [ ! -f $@ ]; then exit 1; fi
	touch $@
aflow_data_nist.cpp: Makefile 
	rm -rf $@
	if [ -f $(EXTRA_FREE)/NIST/ElectronStoppingPower.txt.xz ]; then echo "// $@ automatic generated from NIST" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/NIST/ElectronStoppingPower.txt.xz ]; then echo "#include <string>" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/NIST/ElectronStoppingPower.txt.xz ]; then echo "std::string ElectronStoppingPower_txt=\"\\" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/NIST/ElectronStoppingPower.txt.xz ]; then xzcat $(EXTRA_FREE)/NIST/ElectronStoppingPower.txt.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> $@ ; fi
	if [ -f $(EXTRA_FREE)/NIST/ElectronStoppingPower.txt.xz ]; then echo "\";" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/NIST/PhotonCrossSection.txt.xz ]; then echo "std::string PhotonCrossSection_txt=\"\\" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/NIST/PhotonCrossSection.txt.xz ]; then xzcat $(EXTRA_FREE)/NIST/PhotonCrossSection.txt.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> $@ ; fi
	if [ -f $(EXTRA_FREE)/NIST/PhotonCrossSection.txt.xz ]; then echo "\";" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/NIST/PhotonStoppingPower.txt.xz ]; then echo "std::string PhotonStoppingPower_txt=\"\\" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/NIST/PhotonStoppingPower.txt.xz ]; then xzcat $(EXTRA_FREE)/NIST/PhotonStoppingPower.txt.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> $@ ; fi
	if [ -f $(EXTRA_FREE)/NIST/PhotonStoppingPower.txt.xz ]; then echo "\";" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/NIST/ICSD_List_2014.txt.xz ]; then echo "std::string ICSD_List_txt=\"\\" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/NIST/ICSD_List_2014.txt.xz ]; then xzcat $(EXTRA_FREE)/NIST/ICSD_List_2014.txt.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> $@ ; fi
	if [ -f $(EXTRA_FREE)/NIST/ICSD_List_2014.txt.xz ]; then echo "\";" >> $@ ; fi
	-if [ ! -f $@ ]; then rm -rf $@.xz*; wget http://materials.duke.edu/AFLOW/AFLOW3_FREE/$@.xz && xz -dv $@.xz ; fi
	-if [ ! -f $@ ]; then rm -rf $@.xz*; wget http://aflowlib.duke.edu/AFLOW/AFLOW3_FREE/$@.xz && xz -dv $@.xz ; fi
	if [ -f $@ ]; then GUARD_CHECK=$$(grep '_AFLOW_DATA_NIST_CPP_' $@); if [ -z "$$GUARD_CHECK" ]; then mv $@ $@_tmp; echo "#ifndef _AFLOW_DATA_NIST_CPP_" >> $@; echo "#define _AFLOW_DATA_NIST_CPP_" >> $@; cat $@_tmp >> $@ && rm -f $@_tmp; echo "#endif // _AFLOW_DATA_NIST_CPP_" >> $@; fi; fi
	if [ ! -f $@ ]; then exit 1; fi
	touch $@
aflow_data_aflow_potcars.cpp: Makefile 
	rm -rf $@
	echo "#ifndef _AFLOW_DATA_AFLOW_POTCARS_CPP_" >> $@
	echo "#define _AFLOW_DATA_AFLOW_POTCARS_CPP_" >> $@
	echo "// $@ automatic generated from AFLOW_PSEUDOPOTENTIALS.TXT" >> $@
	echo "#include <string>" >> $@
	if [ -f ../AFLOW3_NONFREE/EXTRA/NONE_AFLOW_PSEUDOPOTENTIALS/AFLOW_PSEUDOPOTENTIALS.TXT.xz ]; then echo "bool AFLOW_PSEUDOPOTENTIALS=1;" >> $@ ; else echo "bool AFLOW_PSEUDOPOTENTIALS=0;" >> $@ ; fi
	echo "std::string AFLOW_PSEUDOPOTENTIALS_TXT=\"\\" >> $@
	#if[ -f ../AFLOW3_NONFREE/EXTRA/NONE_AFLOW_PSEUDOPOTENTIALS/AFLOW_PSEUDOPOTENTIALS.TXT.xz ]; then xzcat ../AFLOW3_NONFREE/EXTRA/NONE_AFLOW_PSEUDOPOTENTIALS/AFLOW_PSEUDOPOTENTIALS.TXT.xz | perl -p -e 's/\n/\\n\\\n/g' >> $@ ; fi
	if [ -f ../AFLOW3_NONFREE/EXTRA/NONE_AFLOW_PSEUDOPOTENTIALS/AFLOW_PSEUDOPOTENTIALS.TXT_XZ_BASE64.xz ]; then xzcat ../AFLOW3_NONFREE/EXTRA/NONE_AFLOW_PSEUDOPOTENTIALS/AFLOW_PSEUDOPOTENTIALS.TXT_XZ_BASE64.xz | perl -p -e 's/\n/\\n\\\n/g' >> $@ ; fi
	echo "\";" >> $@
	echo "std::string AFLOW_PSEUDOPOTENTIALS_LIST_TXT=\"\\" >> $@
	#if [ -f ../AFLOW3_NONFREE/EXTRA/NONE_AFLOW_PSEUDOPOTENTIALS/AFLOW_PSEUDOPOTENTIALS.TXT.xz ]; then xzcat ../AFLOW3_NONFREE/EXTRA/NONE_AFLOW_PSEUDOPOTENTIALS/AFLOW_PSEUDOPOTENTIALS.TXT.xz | grep "POTCAR=" | perl -p -e 's/\n/\\n\\\n/g' >> $@ ; fi
	if [ -f ../AFLOW3_NONFREE/EXTRA/NONE_AFLOW_PSEUDOPOTENTIALS/AFLOW_PSEUDOPOTENTIALS_LIST.TXT.xz ]; then xzcat ../AFLOW3_NONFREE/EXTRA/NONE_AFLOW_PSEUDOPOTENTIALS/AFLOW_PSEUDOPOTENTIALS_LIST.TXT.xz | grep "POTCAR=" | perl -p -e 's/\n/\\n\\\n/g' >> $@ ; fi
	echo "\";" >> $@
	echo "#endif // _AFLOW_DATA_AFLOW_POTCARS_CPP_" >> $@
aflow_data_latex.cpp: Makefile
	rm -rf $@
	echo "#ifndef _AFLOW_DATA_LATEX_CPP_" >> $@
	echo "#define _AFLOW_DATA_LATEX_CPP_" >> $@
	echo "// $@ automatic generated from LATEX" >> $@
	echo "#include <string>" >> $@
	echo "std::string f144468a7ccc2d3a72ba44000715efdb=\"\\" >> $@
	if [ -f ../AFLOW3_NONFREE/EXTRA/LATEX/f144468a7ccc2d3a72ba44000715efdb ]; then cat ../AFLOW3_NONFREE/EXTRA/LATEX/f144468a7ccc2d3a72ba44000715efdb* | perl -p -e 's/\\/\\\\/g' | perl -p -e 's/\"/\\\"/g' | perl -p -e 's/\n/ \\n\\\n/g' >> $@ ; fi
	echo "\";" >> $@
	echo "#endif // _AFLOW_DATA_LATEX_CPP_" >> $@
aflow_data_readme.cpp: README_AFLOW_LICENSE_GPL3.TXT README_AFLOW_ACONVASP.TXT README_AFLOW_APL.TXT README_AFLOW_AGL.TXT README_AFLOW_AEL.TXT README_AFLOW_POCC.TXT README_AFLOW_ANRL.TXT README_AFLOW_COMPARE.TXT README_AFLOW_GFA.TXT README_AFLOW.TXT README_PROTO.TXT README_AFLOW_APENNSY.TXT README_AFLOW_FROZSL.TXT README_AFLOW_SCRIPTING.TXT README_AFLOW_SYM.TXT README_AFLOW_CCE.TXT README_AFLOW_CHULL.TXT README_AFLOW_XAFLOW.TXT README_AFLOW_AFLOWRC.TXT README_AFLOW_EXCEPTIONS.TXT README_CONTRIBS.TXT README_AFLOW_VERSIONS_HISTORY.TXT Makefile
	rm -rf $@
	echo "#ifndef _AFLOW_DATA_README_CPP_" >> $@
	echo "#define _AFLOW_DATA_README_CPP_" >> $@
	echo "// $@ automatic generated" >> $@
	echo "#include <string>" >> $@
	echo "std::string README_AFLOW_LICENSE_GPL3_TXT=\"\\" >> $@
	cat README_AFLOW_LICENSE_GPL3.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> $@
	echo "\";" >> $@
	echo "std::string README_AFLOW_TXT=\"\\" >> $@
	cat README_AFLOW.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> $@
	echo "\";" >> $@
	echo "std::string README_AFLOW_VERSIONS_HISTORY_TXT=\"\\" >> $@
	cat README_AFLOW_VERSIONS_HISTORY.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> $@
	echo "\";" >> $@
	echo "std::string README_AFLOW_PFLOW_TXT=\"\\" >> $@
	cat README_AFLOW_ACONVASP.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> $@
	echo "\";" >> $@
	echo "std::string README_AFLOW_APENNSY_TXT=\"\\" >> $@
	cat README_AFLOW_APENNSY.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> $@
	echo "\";" >> $@
	echo "std::string README_AFLOW_SCRIPTING_TXT=\"\\" >> $@
	cat README_AFLOW_SCRIPTING.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> $@
	echo "\";" >> $@
	echo "std::string README_AFLOW_FROZSL_TXT=\"\\" >> $@
	cat README_AFLOW_FROZSL.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> $@
	echo "\";" >> $@
	echo "std::string README_AFLOW_POCC_TXT=\"\\" >> $@
	cat README_AFLOW_POCC.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> $@
	echo "\";" >> $@
	echo "std::string README_AFLOW_APL_TXT=\"\\" >> $@
	cat README_AFLOW_APL.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> $@
	echo "\";" >> $@
	echo "std::string README_AFLOW_AGL_TXT=\"\\" >> $@
	cat README_AFLOW_AGL.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> $@
	echo "\";" >> $@
	echo "std::string README_AFLOW_AEL_TXT=\"\\" >> $@
	cat README_AFLOW_AEL.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> $@
	echo "\";" >> $@
	echo "std::string README_AFLOW_ANRL_TXT=\"\\" >> $@
	cat README_AFLOW_ANRL.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> $@
	echo "\";" >> $@
	echo "std::string README_AFLOW_COMPARE_TXT=\"\\" >> $@
	cat README_AFLOW_COMPARE.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> $@
	echo "\";" >> $@
	echo "std::string README_AFLOW_GFA_TXT=\"\\" >> $@
	cat README_AFLOW_GFA.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> $@
	echo "\";" >> $@
	echo "std::string README_AFLOW_SYM_TXT=\"\\" >> $@
	cat README_AFLOW_SYM.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> $@
	echo "\";" >> $@
	echo "std::string README_AFLOW_CCE_TXT=\"\\" >> $@
	cat README_AFLOW_CCE.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> $@
	echo "\";" >> $@
	echo "std::string README_AFLOW_CHULL_TXT=\"\\" >> $@
	cat README_AFLOW_CHULL.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> $@
	echo "\";" >> $@
	echo "std::string README_AFLOW_EXCEPTIONS_TXT=\"\\" >> $@
	cat README_AFLOW_EXCEPTIONS.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> $@
	echo "\";" >> $@
	echo "std::string README_PROTO_TXT=\"\\" >> $@
	cat README_PROTO.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> $@
	echo "\";" >> $@
	echo "std::string README_AFLOW_XAFLOW_TXT=\"\\" >> $@
	cat README_AFLOW_XAFLOW.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> $@
	echo "\";" >> $@
	echo "std::string README_AFLOW_AFLOWRC_TXT=\"\\" >> $@
	cat README_AFLOW_AFLOWRC.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> $@
	echo "\";" >> $@
	echo "#endif // _AFLOW_DATA_README_CPP_" >> $@
aflow_data_extra.cpp: aflow_library_icsd.dat.xz Makefile
	rm -rf $@
	echo "#ifndef _AFLOW_DATA_EXTRA_CPP_" >> $@
	echo "#define _AFLOW_DATA_EXTRA_CPP_" >> $@
	echo "#include <string>" >> $@
	echo "std::string icsd_1_ary=\"\\" >> $@ && xzcat $(NIST_LIB)/1-ary.icsd.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> $@ && echo "\";" >> $@	
	echo "std::string icsd_2_ary=\"\\" >> $@ && xzcat $(NIST_LIB)/2-ary.icsd.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> $@ && echo "\";" >> $@	
	echo "std::string icsd_3_ary=\"\\" >> $@ && xzcat $(NIST_LIB)/3-ary.icsd.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> $@ && echo "\";" >> $@	
	echo "std::string icsd_4_ary=\"\\" >> $@ && xzcat $(NIST_LIB)/4-ary.icsd.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> $@ && echo "\";" >> $@	
	echo "std::string icsd_5_ary=\"\\" >> $@ && xzcat $(NIST_LIB)/5-ary.icsd.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> $@ && echo "\";" >> $@	
	echo "std::string icsd_6_ary=\"\\" >> $@ && xzcat $(NIST_LIB)/6-ary.icsd.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> $@ && echo "\";" >> $@	
	echo "std::string icsd_7_ary=\"\\" >> $@ && xzcat $(NIST_LIB)/7-ary.icsd.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> $@ && echo "\";" >> $@	
	echo "std::string icsd_8_ary=\"\\" >> $@ && xzcat $(NIST_LIB)/8-ary.icsd.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> $@ && echo "\";" >> $@	
	echo "std::string icsd_9_ary=\"\\" >> $@ && xzcat $(NIST_LIB)/9-ary.icsd.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> $@ && echo "\";" >> $@	
	echo "std::string icsd_10_ary=\"\\" >> $@ && xzcat $(NIST_LIB)/10-ary.icsd.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> $@ && echo "\";" >> $@
	echo "std::string icsd_11_ary=\"\\" >> $@ && xzcat $(NIST_LIB)/11-ary.icsd.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> $@ && echo "\";" >> $@
	echo "std::string icsd_12_ary=\"\\" >> $@ && xzcat $(NIST_LIB)/12-ary.icsd.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> $@ && echo "\";" >> $@
	echo "std::string icsd_13_ary=\"\\" >> $@ && xzcat $(NIST_LIB)/13-ary.icsd.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> $@ && echo "\";" >> $@
	echo "std::string icsd_14_ary=\"\\" >> $@ && xzcat $(NIST_LIB)/14-ary.icsd.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> $@ && echo "\";" >> $@
	echo "std::string icsd_15_ary=\"\\" >> $@ && xzcat $(NIST_LIB)/15-ary.icsd.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> $@ && echo "\";" >> $@
	echo "#endif // _AFLOW_DATA_EXTRA_CPP_" >> $@
#	echo "std::string Library_ICSD=\"\\" >> $@ && xzcat ../AFLOW3_NONFREE/LIBS/$< | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> $@ && echo "\";" >> $@
#	MAKE aflow_data_libraries.cpp, do not make another rule, aflow_data_libraries.cpp is already a rule
	rm -rf aflow_data_libraries.cpp
	echo "#ifndef _AFLOW_DATA_LIBRARIES_CPP_" >> aflow_data_libraries.cpp
	echo "#define _AFLOW_DATA_LIBRARIES_CPP_" >> aflow_data_libraries.cpp
	echo "#include <string>" >> aflow_data_libraries.cpp
	echo "std::string aflowlib_lib0=\"\\" >> aflow_data_libraries.cpp && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib0.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_libraries.cpp && echo "\";" >> aflow_data_libraries.cpp	
	echo "std::string aflowlib_lib1=\"\\" >> aflow_data_libraries.cpp && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib1.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_libraries.cpp && echo "\";" >> aflow_data_libraries.cpp	
	echo "std::string aflowlib_lib2=\"\\" >> aflow_data_libraries.cpp && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib2.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_libraries.cpp && echo "\";" >> aflow_data_libraries.cpp	
	echo "std::string aflowlib_lib3=\"\\" >> aflow_data_libraries.cpp && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib3.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_libraries.cpp && echo "\";" >> aflow_data_libraries.cpp
	echo "std::string aflowlib_lib4=\"\\" >> aflow_data_libraries.cpp && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib4.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_libraries.cpp && echo "\";" >> aflow_data_libraries.cpp
	echo "std::string aflowlib_lib5=\"\\" >> aflow_data_libraries.cpp && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib5.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_libraries.cpp && echo "\";" >> aflow_data_libraries.cpp
	echo "std::string aflowlib_lib6=\"\\" >> aflow_data_libraries.cpp && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib6.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_libraries.cpp && echo "\";" >> aflow_data_libraries.cpp
	echo "std::string aflowlib_lib7=\"\\" >> aflow_data_libraries.cpp && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib7.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_libraries.cpp && echo "\";" >> aflow_data_libraries.cpp
	echo "std::string aflowlib_lib8=\"\\" >> aflow_data_libraries.cpp && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib8.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_libraries.cpp && echo "\";" >> aflow_data_libraries.cpp
	echo "std::string aflowlib_lib9=\"\\" >> aflow_data_libraries.cpp && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib9.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_libraries.cpp && echo "\";" >> aflow_data_libraries.cpp
	echo "std::string aflowlib_icsd=\"\\" >> aflow_data_libraries.cpp && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_icsd.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_libraries.cpp && echo "\";" >> aflow_data_libraries.cpp
	echo "#endif // _AFLOW_DATA_LIBRARIES_CPP_" >> aflow_data_libraries.cpp
#aflow_data: aflow_data_readme.cpp aflow_data_calculated.cpp aflow_data_htqc.cpp aflow_data_libraries.cpp aflow_data_stokes.cpp aflow_data_nist.cpp aflow_data_latex.cpp aflow_data_aflow_potcars.cpp aflow_data.cpp
aflow_data: aflow_data.cpp aflow_data_readme.o aflow_data_calculated.o aflow_data_htqc.o aflow_data_libraries.o aflow_data_stokes.o aflow_data_nist.o aflow_data_latex.o aflow_data_aflow_potcars.o Makefile
	rm -f aflow_data_extra.cpp
	echo "#ifndef _AFLOW_DATA_EXTRA_CPP_" >> aflow_data_extra.cpp
	echo "#define _AFLOW_DATA_EXTRA_CPP_" >> aflow_data_extra.cpp
#	#[CO20200508 - NO LONGER NEEDED, the file will be created now anyway with ifndef guards]echo "" >> aflow_data_extra.cpp
	if [ -f ../AFLOW3_NONFREE/LIBS/aflow_library_icsd.dat.xz ]; then echo "#define defined_Library_ICSD" >> aflow_data_extra.cpp; fi
	if [ -f ../AFLOW3_NONFREE/LIBS/aflow_library_icsd.dat.xz ]; then echo "std::string Library_ICSD=\"\\" >> aflow_data_extra.cpp; fi
	if [ -f ../AFLOW3_NONFREE/LIBS/aflow_library_icsd.dat.xz ]; then xzcat ../AFLOW3_NONFREE/LIBS/aflow_library_icsd.dat.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_extra.cpp; fi
	if [ -f ../AFLOW3_NONFREE/LIBS/aflow_library_icsd.dat.xz ]; then echo "\";" >> aflow_data_extra.cpp; fi
	echo "#endif // _AFLOW_DATA_EXTRA_CPP_" >> aflow_data_extra.cpp
	$(MAKE) aflow_data.o
#	$(CPP) $(VERS) -D_AFLOW_FILE_NAME_=\""aflow_data.cpp"\" $(INCLUDE) $(CCFLAGS) $(OPTS) $(ARCH) aflow_data.cpp -c -o aflow_data.o
	@echo "linking aflow_data*.o"
	@$(CPP) $(VERS) $(INCLUDE) $(CCFLAGS) $(OPTS) -o aflow_data aflow_data.o aflow_data_*.o $(FLAGS_STATIC_POST)
.PHONY: data
data: aflow_data
.PHONY: amir
amir:	aflow_data.cpp aflow_data_readme.cpp aflow_data_calculated.cpp aflow_data_htqc.cpp aflow_data_libraries.cpp aflow_data_stokes.cpp aflow_data_nist.cpp aflow_data_latex.cpp aflow_data_aflow_potcars.cpp aflow_data_extra.cpp
	$(CPP) $(VERS) $(INCLUDE) $(CCFLAGS) $(OPTS) -DICSD -o aflow_data $<
###################################################################################################

###################################################################################################
aflow_matlab_funcs.cpp: Makefile
	rm -rf $@
	if [ -f $(EXTRA_FREE)/MATLAB/param.m ]; then echo "// $@ automatic generated from MATLAB/*" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/MATLAB/param.m ]; then echo "#include <sstream>" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/MATLAB/param.m ]; then echo "// $@ automatic generated from MATLAB/param.m" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/MATLAB/param.m ]; then echo "std::string MATLAB_FUNCS_param(void){" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/MATLAB/param.m ]; then echo "   std::stringstream strstream; strstream << \"\\" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/MATLAB/param.m ]; then cat $(EXTRA_FREE)/MATLAB/param.m | perl -p -e 's/\"/\\\"/g' | perl -p -e 's/\\G/\\\\G/g' | perl -p -e 's/\\S/\\\\S/g' | perl -p -e 's/\n/\" << std::endl;\n   strstream << \"/g' >> $@ ; fi
	if [ -f $(EXTRA_FREE)/MATLAB/param.m ]; then echo "\" << std::endl; return strstream.str();};" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/MATLAB/plotband.m ]; then echo "// $@ automatic generated from MATLAB/plotband.m" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/MATLAB/plotband.m ]; then echo "// std::string MATLAB_FUNCS_plotband(std::string DIRECTORY,std::string OPTION1){" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/MATLAB/plotband.m ]; then echo "std::string MATLAB_FUNCS_plotband(void){" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/MATLAB/plotband.m ]; then echo "   std::stringstream strstream; strstream << \"\\" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/MATLAB/plotband.m ]; then cat $(EXTRA_FREE)/MATLAB/plotband.m | perl -p -e 's/\"/\\\"/g' | perl -p -e 's/\\G/\\\\G/g' | perl -p -e 's/\\S/\\\\S/g' | perl -p -e 's/\n/\" << std::endl;\n   strstream << \"/g' | perl -p -e 's/MATLAB_FUNCS_DIRECTORY/\" << DIRECTORY << \"/g' | perl -p -e 's/DOS_SCALE_MODE/\" << OPTION1 << \"/g' >> $@ ; fi
	if [ -f $(EXTRA_FREE)/MATLAB/plotband.m ]; then echo "\" << std::endl; return strstream.str();};" >> $@ ; fi
	-if [ ! -f $@ ]; then rm -rf $@.xz*; wget http://materials.duke.edu/AFLOW/AFLOW3_FREE/$@.xz && xz -dv $@.xz ; fi
	-if [ ! -f $@ ]; then rm -rf $@.xz*; wget https://aflowlib.duke.edu/AFLOW/AFLOW3_FREE/$@.xz && xz -dv $@.xz ; fi
	if [ -f $@ ]; then GUARD_CHECK=$$(grep '_AFLOW_MATLAB_FUNCS_CPP_' $@); if [ -z "$$GUARD_CHECK" ]; then mv $@ $@_tmp; echo "#ifndef _AFLOW_MATLAB_FUNCS_CPP_" >> $@; echo "#define _AFLOW_MATLAB_FUNCS_CPP_" >> $@; cat $@_tmp >> $@ && rm -f $@_tmp; echo "#endif // _AFLOW_MATLAB_FUNCS_CPP_" >> $@; fi; fi
	if [ ! -f $@ ]; then exit 1; fi
	touch $@
aflow_gnuplot_funcs.cpp: Makefile
	rm -rf $@
	if [ -f $(EXTRA_FREE)/GNUPLOT/plotbz.sh ]; then echo "// $@ automatic generated from GNUPLOT/*" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/GNUPLOT/plotbz.sh ]; then echo "#include <sstream>" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/GNUPLOT/plotbz.sh ]; then echo "// $@ automatic generated from GNUPLOT/plotbz.sh" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/GNUPLOT/plotbz.sh ]; then echo "//std::string GNUPLOT_FUNCS_plotbz(std::string DIRECTORY,std::string OPTION1){" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/GNUPLOT/plotbz.sh ]; then echo "std::string GNUPLOT_FUNCS_plotbz(void){" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/GNUPLOT/plotbz.sh ]; then echo "   std::stringstream strstream; strstream << \"\\" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/GNUPLOT/plotbz.sh ]; then cat $(EXTRA_FREE)/GNUPLOT/plotbz.sh | perl -p -e 's/\\/\\\\/g' | perl -p -e 's/\"/\\\"/g' | perl -p -e 's/\\G/\\\\G/g' | perl -p -e 's/\\S/\\\\S/g' | perl -p -e 's/\n/\" << std::endl;\n   strstream << \"/g' | perl -p -e 's/GNUPLOT_FUNCS_DIRECTORY/\" << DIRECTORY << \"/g' | perl -p -e 's/GNUPLOT_FUNCS_OPTION1/\" << OPTION1 << \"/g' >> $@ ; fi
	if [ -f $(EXTRA_FREE)/GNUPLOT/plotbz.sh ]; then echo "\" << std::endl; return strstream.str();};" >> $@ ; fi
	-if [ ! -f $@ ]; then rm -rf $@.xz*; wget http://materials.duke.edu/AFLOW/AFLOW3_FREE/$@.xz && xz -dv $@.xz ; fi
	-if [ ! -f $@ ]; then rm -rf $@.xz*; wget https://aflowlib.duke.edu/AFLOW/AFLOW3_FREE/$@.xz && xz -dv $@.xz ; fi
	if [ -f $@ ]; then GUARD_CHECK=$$(grep '_AFLOW_GNUPLOT_FUNCS_CPP_' $@); if [ -z "$$GUARD_CHECK" ]; then mv $@ $@_tmp; echo "#ifndef _AFLOW_GNUPLOT_FUNCS_CPP_" >> $@; echo "#define _AFLOW_GNUPLOT_FUNCS_CPP_" >> $@; cat $@_tmp >> $@ && rm -f $@_tmp; echo "#endif // _AFLOW_GNUPLOT_FUNCS_CPP_" >> $@; fi; fi
	if [ ! -f $@ ]; then exit 1; fi
	touch $@
###################################################################################################

###################################################################################################

aflow_mix_nomix_table.cpp: Makefile $(AFLOW_HPPS) $(AUROSTD_HPPS)
	rm -rf nomix.tmp
	find $(PROJECT_DIRECTORY_GNDSTATE)/NOMIX/ -name aflow.in | sed "s/\/common\/$(PROJECT_DIRECTORY_GNDSTATE)\/NOMIX\/LIB2\/LIB\///g" | sed "s/aflow.in//g" >> nomix.tmp
	rm -rf $@
	echo "#ifndef _AFLOW_MIX_NOMIX_TABLE_CPP_" >> $@
	echo "#define _AFLOW_MIX_NOMIX_TABLE_CPP_" >> $@
	echo "// $@ automatic generated " >> $@
	echo "#include <string>" >> $@
	echo "std::string NomixTable=\"\\" >> $@
	cat nomix.tmp | sort | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> $@
	echo "\";" >> $@
	echo "#endif // _AFLOW_MIX_NOMIX_TABLE_CPP_" >> $@
	chmod 755 $@
	rm -rf nomix.tmp
###################################################################################################

###################################################################################################
aflow_phonons.o: aflow_phonons.cpp $(AFLOW_HPPS) $(AUROSTD_HPPS) Makefile
	$(CPP) $(VERS) $(INCLUDE) $(CCFLAGS) $(OPTS) aflow_phonons.cpp -c
aurostd_xtemplates.o: aurostd_xtemplates.cpp $(AFLOW_HPPS) $(AUROSTD_HPPS) Makefile
	$(CPP) $(VERS) $(INCLUDE) $(CCFLAGS) $(OPTS) $< -c
###################################################################################################

###################################################################################################
#HTRESOURCES
aflow_contrib_stefano_htcurriculum.o: aflowlib.h $(AFLOW_HPPS) $(AUROSTD_HPPS) Makefile
	$(CPP) $(VERS) $(INCLUDE) $(CCFLAGS) $(OPTS)  -c
aflow_contrib_stefano_htsecurity.o: aflow_contrib_stefano_htsecurity.cpp aflowlib.h $(AFLOW_HPPS) $(AUROSTD_HPPS) Makefile
	$(CPP) $(VERS) $(INCLUDE) $(CCFLAGS) $(OPTS) $< -c
###################################################################################################

###################################################################################################
.PHONY: install
install: all
	cp -f aflow $(BIN)/aflow.$(VERSNUMBER)
	cp -f aflow_data $(BIN)/
	cp -f $(BIN)/aflow.$(VERSNUMBER) $(BIN)/aflow
	ls -las $(BIN)/aflow*
	$(MAKE) check
.PHONY: install_auro
install_auro: install_user
	./aflow --version | grep AFLOW >> $(BIN)/etc_aflow_status.txt
	chmod 644 $(BIN)/etc_aflow_status.txt
	./aflow --calculated | grep Calcs >> $(BIN)/etc_aflowlib_status.txt
	chmod 644 $(BIN)/etc_aflowlib_status.txt
	./aflow --apool > $(BIN)/etc_aflowlib_apool.txt
	chmod 644 $(BIN)/etc_aflowlib_apool.txt
.PHONY: install_root
install_root: all
	cp -f aflow /usr/local/bin/aflow.$(VERSNUMBER)
	cp -f aflow_data /usr/local/bin/aflow_data
	cp -f /usr/local/bin/aflow.$(VERSNUMBER) /usr/local/bin/aflow
	ls -las /usr/local/bin/aflow*
	$(MAKE) check
.PHONY: niet
niet:	
	dropbox stop
	$(MAKE) clean
	$(MAKE) icsdclean
	$(MAKE) extra_clean
	if [ -f $(EXTRA_FREE)/GUS/Makefile ]; then $(MAKE) -C $(EXTRA_FREE)/GUS realclean; fi
	echo "" > $(ICSD_LIB)/README_LIBRARY_ICSD1.TXT
	echo "" > $(ICSD_LIB)/README_LIBRARY_ICSD2.TXT
	echo "" > $(ICSD_LIB)/README_LIBRARY_ICSD3.TXT
	echo "" > $(ICSD_LIB)/README_LIBRARY_ICSD4.TXT
	echo "" > $(ICSD_LIB)/README_LIBRARY_ICSD5.TXT
	echo "" > $(ICSD_LIB)/README_LIBRARY_ICSD6.TXT
	echo "" > $(ICSD_LIB)/README_LIBRARY_ICSD7.TXT
	echo "" > $(ICSD_LIB)/README_LIBRARY_ICSD8.TXT
	$(MAKE) web
	$(MAKE) backup
	$(MAKE) icsdclean
	$(MAKE) -j 32
	rm -fv aflow_data_calculated.cpp* && $(MAKE) aflow_data_calculated.cpp && cp -f aflow_data_calculated.cpp /www/AFLOW/AFLOW3_FREE/ && xz -9fv -T8 /www/AFLOW/AFLOW3_FREE/aflow_data_calculated.cpp
	rm -fv aflow_data_stokes.cpp* && $(MAKE) aflow_data_stokes.cpp && cp -f aflow_data_stokes.cpp /www/AFLOW/AFLOW3_FREE/ && xz -9fv -T8 /www/AFLOW/AFLOW3_FREE/aflow_data_stokes.cpp
	rm -fv aflow_data_nist.cpp* && $(MAKE) aflow_data_nist.cpp && cp -f aflow_data_nist.cpp /www/AFLOW/AFLOW3_FREE/ && xz -9fv -T8 /www/AFLOW/AFLOW3_FREE/aflow_data_nist.cpp
	rm -fv aflow_matlab_funcs.cpp* && $(MAKE) aflow_matlab_funcs.cpp && cp -f aflow_matlab_funcs.cpp /www/AFLOW/AFLOW3_FREE/ && xz -9fv -T8 /www/AFLOW/AFLOW3_FREE/aflow_matlab_funcs.cpp
	rm -fv aflow_gnuplot_funcs.cpp* && $(MAKE) aflow_gnuplot_funcs.cpp && cp -f aflow_gnuplot_funcs.cpp /www/AFLOW/AFLOW3_FREE/ && xz -9fv -T8 /www/AFLOW/AFLOW3_FREE/aflow_gnuplot_funcs.cpp
	dropbox start
.PHONY: backup
backup: realclean clean icsdclean 
	if [ -f $(EXTRA_FREE)/GUS/Makefile ]; then $(MAKE) -C $(EXTRA_FREE)/GUS realclean ; fi
	mkdir -p aflow.$(VERSNUMBER)
#	mkdir -p aflow.$(VERSNUMBER)/SOURCES
	cp -drvf a*cpp a*.h README* Makefile* APL* ANRL* DATA* AUROSTD* SQLITE* SYMBOLICCPLUSPLUS* aflow.$(VERSNUMBER)
# no EXTRA in SRC	if [ -d $(EXTRA_FREE) ]; then cp -drvf $(EXTRA_FREE) aflow.$(VERSNUMBER) ; fi
	if [ -d SOURCE ]; then cp -drvf SOURCE aflow.$(VERSNUMBER) ; fi
	tar Jcfvp aflow.$(VERSNUMBER).tar.xz aflow.$(VERSNUMBER)
	rm -rfv $(BACKUP)/aflow.$(VERSNUMBER) $(BACKUP)/aflow.$(VERSNUMBER).tar.xz
	if [ -d /www/AFLOW ]; then cp -f aflow.$(VERSNUMBER).tar.xz /www/AFLOW ; fi
	mv aflow.$(VERSNUMBER).tar.xz $(BACKUP)
	rm -rfv aflow.$(VERSNUMBER)
.PHONY: www
www:	web
.PHONY: web
web:	realclean
	rm -rvf aflow3.tar.xz aflow_test.tar.xz aflow3 ~/aflow3.tar.xz ~/aflow_test.tar.xz ~/aflow3
	tar Jcfvp ~/aflow3.tar.xz Makefile* *cpp *.h READ* APL* DATA* ANRL* AUROSTD* SQLITE* SYMBOLICCPLUSPLUS*
# no EXTRA IN SRC	if [ -d $(EXTRA_FREE) ]; then tar Jcfvp ~/aflow3.tar.xz Makefile *cpp *.h *.js *.json *.py *.txt READ* APL ANRL* AUROSTD SQLITE $(EXTRA_FREE) ; fi
# no other SOURCES in SRC	if [ -d SOURCES ]; then tar Jcfvp ~/aflow3.tar.xz Makefile *cpp *.h *.js *.json *.py *.txt READ* APL ANRL* AUROSTD SQLITE $(EXTRA_FREE) SOURCES ; fi
	if [ -d /www/AFLOW ]; then cp README_AFLOW.TXT /www/AFLOW ; fi
	if [ -d /www/AFLOW ]; then cp README_AFLOW_AEL.TXT /www/AFLOW ; fi
	if [ -d /www/AFLOW ]; then cp README_AFLOW_AGL.TXT /www/AFLOW ; fi
	if [ -d /www/AFLOW ]; then cp README_AFLOW_ANRL.TXT /www/AFLOW ; fi
	if [ -d /www/AFLOW ]; then cp README_AFLOW_APENNSY.TXT /www/AFLOW ; fi
	if [ -d /www/AFLOW ]; then cp README_AFLOW_APL.TXT /www/AFLOW ; fi
	if [ -d /www/AFLOW ]; then cp README_AFLOW_QHA_SCQHA_QHA3P.TXT /www/AFLOW ; fi
	if [ -d /www/AFLOW ]; then cp README_AFLOW_ACONVASP.TXT /www/AFLOW ; fi
	if [ -d /www/AFLOW ]; then cp README_AFLOW_COMPARE.TXT /www/AFLOW ; fi
	if [ -d /www/AFLOW ]; then cp README_AFLOW_GFA.TXT /www/AFLOW ; fi
	if [ -d /www/AFLOW ]; then cp README_AFLOW_FROZSL.TXT /www/AFLOW ; fi
	if [ -d /www/AFLOW ]; then cp README_AFLOW_LICENSE_GPL3.TXT /www/AFLOW ; fi
	if [ -d /www/AFLOW ]; then cp README_AFLOW_POCC.TXT /www/AFLOW ; fi
	if [ -d /www/AFLOW ]; then cp README_AFLOW_SCRIPTING.TXT /www/AFLOW ; fi
	if [ -d /www/AFLOW ]; then cp README_AFLOW_SYM.TXT /www/AFLOW ; fi
	if [ -d /www/AFLOW ]; then cp README_AFLOW_CCE.TXT /www/AFLOW ; fi
	if [ -d /www/AFLOW ]; then cp README_AFLOW_CHULL.TXT /www/AFLOW ; fi
	if [ -d /www/AFLOW ]; then cp README_AFLOW_EXCEPTIONS.TXT /www/AFLOW ; fi
	if [ -d /www/AFLOW ]; then cp README_AFLOW_VERSIONS_HISTORY.TXT /www/AFLOW ; fi
	if [ -d /www/AFLOW ]; then cp README_AFLOW_XAFLOW.TXT /www/AFLOW ; fi
	if [ -d /www/AFLOW ]; then cp README_AFLOW_AFLOWRC.TXT /www/AFLOW ; fi
	if [ -d /www/AFLOW ]; then cp README_CONTRIBS.TXT /www/AFLOW ; fi
	if [ -d /www/AFLOW ]; then cp README_PROTO.TXT /www/AFLOW ; fi
	if [ -f ../AFLOW3_NONFREE/EXTRA/FROZSL/Makefile ]; then $(MAKE) -C ../AFLOW3_NONFREE/EXTRA/FROZSL; fi
	if [ -f ../AFLOW3_NONFREE/EXTRA/FROZSL/frozsl ]; then cp ../AFLOW3_NONFREE/EXTRA/FROZSL/frozsl /www/AFLOW ; fi
	if [ -f ../AFLOW3_NONFREE/EXTRA/FROZSL/frozsl_init ]; then cp ../AFLOW3_NONFREE/EXTRA/FROZSL/frozsl_init /www/AFLOW ; fi
	if [ -f ../AFLOW3_NONFREE/EXTRA/FINDSYM/Makefile ]; then $(MAKE) -C ../AFLOW3_NONFREE/EXTRA/FINDSYM; fi
	if [ -f ../AFLOW3_NONFREE/EXTRA/FINDSYM/findsym ]; then cp ../AFLOW3_NONFREE/EXTRA/FINDSYM/findsym /www/AFLOW ; fi
	if [ -d /www/auro ]; then mv -f  ~/aflow3.tar.xz /www/auro ; fi
.PHONY: wget
wget:	realclean
	rm -rf aflow3.tar.xz* aflow_test.tar.xz*
	-if [ ! -f aflow3.tar.xz ]; then wget http://materials.duke.edu/auro/aflow3.tar.xz; fi
	-if [ ! -f aflow3.tar.xz ]; then wget https://aflowlib.duke.edu/AFLOW/AFLOW3_FREE/aflow3.tar.xz; fi
	if [ ! -f aflow3.tar.xz ]; then exit 1; fi
	tar xfvp aflow3.tar.xz
	rm -rf aflow3.tar.xz 
.PHONY: distroget
distroget: wget
	rm -rf TESTS PLATON* FINDSYM* FROZSL*
	cat Makefile | sed "s/Makefile /Makefile /g" | sed "s/\*\.a//g" > M
	mv M Makefile
.PHONY: patch_makefile1
patch_makefile1: Makefile
	cat Makefile | sed "s/ \-Wno\-unused\-but\-set\-variable//g" > M	
	mv M Makefile
.PHONY: patch_makefile2
patch_makefile2: Makefile
	cat Makefile | sed "s/\-O1//g" | sed "s/\-O2//g" | sed "s/\-O3//g" > M	
	mv M Makefile
###################################################################################################

###################################################################################################
.PHONY: clean
clean:
	rm -rf *~ *.o *.o-* *.rpo AUROSTD/*~ AUROSTD/*.o AUROSTD/*.rpo ./APL/*~ ./APL/*o ./APL/*a ./ANRL/*~ ./ANRL/*o ./ANRL/*a ./ANRL*/*~ ./ANRL*/*o ./ANRL*/*a ./SQLITE/*o SYMBOLICCPLUSPLUS/*~ SYMBOLICCPLUSPLUS/*.o aflow aflow_data aflowd aflow2 aflowd2 aflow1 aflowd1 aflow.exe aconvasp apennsy aflow_convasp* *.eps *.pdf *ps *aux *dvi *log *m *tex nohup* core* VASP_TEST* GRND_TEST* xa xb xc xd xe xf aflow.out* core* *~ \#* *.o aflow_pennsy *.eps *.pdf *ps *aux *dvi *log *m *tex *exe aurostd_xtensor_extra.cpp aflow_data_*cpp xaflowauto boo* xlist* _xdat `find . -name *~`
	rm -rfv `find . -name "LOCK*" | grep -v TESTS | grep -v SOURCES` *.out aflow.rasmol.xyz* aflow_potcar.cpp aflow_xpotcar.cpp aflow3.tar.xz ./AFLOWDATA ./LIBRARY* *html *toc *m *.jpg *.gif *.png
	rm -rfv ./APL/*~ ./ANRL/*~ ./ANRL*/*~ ./SQLITE/*~ ./SYMBOLICCPLUSPLUS/*~ ../AFLOW3_NONFREE/BIN/* ./EXTRA/GUS/*~ ./EXTRA/GUS/*.o
.PHONY: realclean
realclean: clean
	$(MAKE) clean_autogen
	rm -rf aurostd.cpp aurostd.h *conflic* 
	if [ -f $(EXTRA_FREE)/GUS/Makefile ]; then $(MAKE) -C $(EXTRA_FREE)/GUS realclean; fi
	if [ -f ../AFLOW3_NONFREE/EXTRA/FINDSYM/Makefile ]; then $(MAKE) -C ../AFLOW3_NONFREE/EXTRA/FINDSYM realclean; fi
	if [ -f ../AFLOW3_NONFREE/EXTRA/FROZSL/Makefile ]; then $(MAKE) -C ../AFLOW3_NONFREE/EXTRA/FROZSL clean; fi
	if [ -f ../AFLOW3_NONFREE/EXTRA/FROZSL/Makefile ]; then rm -fv ../AFLOW3_NONFREE/EXTRA/FROZSL/*conflic*; fi
	if [ -f $(EXTRA_FREE)/APL2AGR/Makefile ]; then $(MAKE) -C $(EXTRA_FREE)/APL2AGR realclean; fi
.PHONY: dataclean
dataclean:
	rm -f aflow_data aflow_data_*
datainstall:
	rm -f aflow_data aflow_data_*
	make data
	cp -f aflow_data $(BIN)/
	ls -las $(BIN)/aflow_data*
redata: datainstall
.PHONY: clean_autogen
clean_autogen:
	rm -rf aflow_data_readme.cpp aflow_readme_aflow.cpp aflow_readme_aconvasp.cpp aflow_readme_apennsy.cpp aflow_readme_proto.cpp aflow_xproto_icsd_lib.cpp aflow_matlab_funcs.cpp aflow_gnuplot_funcs.cpp aflow_data_*cpp *.dat
#CO20230507 - keep python/json/js files here as transitory branches still generate these files
	rm -rf aflow_cce_python.cpp aflow_chull_jupyter_plotter.cpp aflow_chull_python.cpp aflow_sym_python.cpp aflow_xtalfinder_python.cpp aflow_chull_jupyter.cpp aflow_chull_jupyter_requirements.cpp aflowlib_webapp_bands.cpp aflowlib_webapp_entry.cpp
###################################################################################################

###################################################################################################
# icsd
.PHONY: icsdclean
icsdclean: clean realclean
	echo "dummy" > $(ICSD_LIB)//dummy.txt
	rm $(ICSD_LIB)//*
.PHONY: icsd
icsd:	$(ICSD_LIB)/README_LIBRARY_ICSD1.TXT $(ICSD_LIB)/README_LIBRARY_ICSD2.TXT $(ICSD_LIB)/README_LIBRARY_ICSD3.TXT $(ICSD_LIB)/README_LIBRARY_ICSD4.TXT $(ICSD_LIB)/README_LIBRARY_ICSD5.TXT $(ICSD_LIB)/README_LIBRARY_ICSD6.TXT $(ICSD_LIB)/README_LIBRARY_ICSD7.TXT $(ICSD_LIB)/README_LIBRARY_ICSD8.TXT
.PHONY: icsd_test
icsd_test: aflow
	cat $(ICSD_LIB)/README_LIBRARY_ICSD2.TXT | grep PROTOTYPE | sed "s/\[/\n/g" | grep PROTOTYPE | sed "s/PROTOTYPE/\.\/aflow --proto_icsd /g" > xgo
	sh ./xgo > /dev/null
###################################################################################################

###################################################################################################
# PENNSY PENNSY PENNSY PENNSY PENNSY PENNSY PENNSY PENNSY PENNSY PENNSY PENNSY
.PHONY: energy
energy:	all
	./aflow --energy > apennsy.tex
	latex apennsy.tex | grep -v Overfull | grep -v cmtt | grep -v type | grep -v "\[\]" | grep "\["
.PHONY: psenergy
psenergy: all
	./aflow --psenergy > psapennsy.tex
	latex psapennsy.tex
	dvips psapennsy.dvi -o psapennsy.eps
	ps2pdf psapennsy.eps
	cp psapennsy*ps EPS
	cp psapennsy*pdf EPS
.PHONY: MIX
MIX:	aflow
	./aflow --libraryX --mix > aflow_nomix.$(TIME).cpp
	chmod 755 aflow_nomix.$(TIME).cpp
	rm -rf aflow_nomix.cpp
	ln -sf aflow_nomix.$(TIME).cpp aflow_nomix.cpp
	grep -c UNKNOWN aflow_nomix.2*
.PHONY: mix
mix:	MIX
.PHONY: NOMIX
NOMIX:	aflow
	cat aflow_nomix.cpp | grep NOMIX | sed "s/\/\//\n/g" | grep calcs | sed "s/\=/\n/g" | grep calcs | sed "s/ calcs/xxxx/g" | sed "s/ /cd \.\//g" | sed "s/xxxx/ \&\& aflow --multi \&\& cd \.\.\//g" > xnomix
	scp xnomix fslcollab8@mary:LIBS/LIB2/LIB
	rm -rf xnomix
.PHONY: nomix
nomix: NOMIX
.PHONY: EXPERIMENTS
EXPERIMENTS: aflow
	./aflow --libraryX --experiments > aflow_mix_experiments.$(TIME).cpp
.PHONY: experiments
experiments: EXPERIMENTS
.PHONY: TABLE
TABLE:aflow
	./aflow --libraryX --table > aflow_table.$(TIME).tex
.PHONY: table
table: TABLE
.PHONY: xaflowauto
xaflowauto: aflow
	./aflow --libraryX --xaflowauto > xaflowauto
	scp xaflowauto fslcollab8@marylou5.byu.edu:LIBS
.PHONY: boo
boo: aflow
	rm -rf boo
	./aflow --calculated=icsd --remove | sed "s/\.\//\.\/ICSD\//g" >> boo
	./aflow --calculated=lib0 --remove | sed "s/\.\//\.\/LIB0\/LIB\//g" >> boo
	./aflow --calculated=lib1 --remove | sed "s/\.\//\.\/LIB1\/LIB\//g" >> boo
	./aflow --calculated=lib2 --remove | sed "s/\.\//\.\/LIB2\/LIB\//g" >> boo
	./aflow --calculated=lib3 --remove | sed "s/\.\//\.\/LIB3\/LIB\//g" >> boo
	./aflow --calculated=lib4 --remove | sed "s/\.\//\.\/LIB4\/LIB\//g" >> boo
	./aflow --calculated=lib5 --remove | sed "s/\.\//\.\/LIB5\/LIB\//g" >> boo
	./aflow --calculated=lib6 --remove | sed "s/\.\//\.\/LIB6\/LIB\//g" >> boo
	./aflow --calculated=lib7 --remove | sed "s/\.\//\.\/LIB7\/LIB\//g" >> boo
	./aflow --calculated=lib8 --remove | sed "s/\.\//\.\/LIB8\/LIB\//g" >> boo
	./aflow --calculated=lib9 --remove | sed "s/\.\//\.\/LIB9\/LIB\//g" >> boo
#	scp boo aflow@lonsdale.tchpc.tcd.ie:LIB3/
	scp boo fslcollab8@mary:LIBS/
.PHONY: WEB
WEB:	aflow
	rm -rf *eps *pdf apennsy.m *html
	#	./aflow --libraryX --QUIET --print=pdf --hull > apennsy.m
	./aflow --libraryX --QUIET --print=gif --print=jpg --hull --oss=cout > apennsy.m
	./aflow --libraryX --QUIET --print=html > apennsy.html
	/usr/local/bin/matlab -nodesktop -r apennsy
	mv *jpg *gif /www/AFLOW
	rm -rf *eps
.PHONY: AFLOW_OLD
APOOL_OLD: aflow aflow_apennsy_main.cpp .cpp
	rm -rf /www/php/apool.html
	cat /www/php/apool_pre.html >> /www/php/apool.html
#	./aflow --apool_test >> /www/php/apool.html
	./aflow --apool >> /www/php/apool.html
	cat /www/php/apool_post.html >> /www/php/apool.html	
	./aflow --apool | grep available | sed "s/ checked//g" | sed "s/></XXXX\n/g" | grep radio | sed "s/<input type=\"radio\" name=\"alloy\" value=//g" | sed "s/XXXX/,/g" > /home/auro/work/AWRAPPER/awrapper_tables_public.cpp
#
#
#	rm -rf /www/php/apool2.html
#	cat /www/php/apool_pre.html >> /www/php/apool2.html
#	./aflow --apool_private >> /www/php/apool2.html
#	cat /www/php/apool_post.html >> /www/php/apool2.html	
#	#
.PHONY: shull
shull:	all
	./aflow --library2 --shull --oss=cout > apennsys.m
.PHONY: vasp
vasp:	all
	./aflow --VASPIN | sed "s/\t/ /g" | sed "s/  / /g" | sed "s/  / /g" > structures.tex
###################################################################################################

###################################################################################################
#EMAIL VERSION UPDATES
#CO20190131 - START
.PHONY: send_version_updates_email
send_version_updates_email:
	rm -f version_updates_email.txt
	$(MAKE) version_updates_email.txt
	[ -f version_updates_email.txt ] && [ -s version_updates_email.txt ]
	send_aflow_code_update_email $(VERSNUMBER)
	rm -f version_updates_email.txt
#CO20190131 - STOP
#CO20190218 - START
version_updates_email.txt: README_AFLOW_VERSIONS_HISTORY.TXT Makefile
	rm -f $(@:.txt=_tmp.txt)
	@file_content=""; \
	found_content=false; \
	while IFS= read -r line; do \
		if [ "$$found_content" = true ]; then \
			if [ -n "$$(echo $$line | grep "Location" )" ] && [ -z "$$(echo $$line | grep "$(VERSNUMBER)" )" ]; then \
				found_content=false; \
				break; \
			fi; \
			file_content="$${file_content}\n$${line}"; \
		elif [ -n "$$( echo "$$line" | grep "$(VERSNUMBER)" )" ]; then \
			found_content=true; \
			file_content="$$line"; \
		fi; \
	done < README_AFLOW_VERSIONS_HISTORY.TXT; \
	line_count=$$(echo $$file_content | wc -l | sed 's/^ *//'); \
	IFS=; \
	file_content=$$(echo $$file_content | head -n $$(( $$line_count - 1 ))); \
	if [ -z "$$file_content" ]; then echo "ERROR - No relevant version content found. Exiting."; exit 1; fi; \
	echo $$file_content >> $(@:.txt=_tmp.txt);
	rm -f $@;
	echo "Dear AFLOW user," >> $@
	echo "" >> $@
	echo "A new version of the AFLOW C++ code is now available (v$(VERSNUMBER))." >> $@
	echo "Update details:" >> $@
	echo "" >> $@
	cat $(@:.txt=_tmp.txt) >> $@; rm -f $(@:.txt=_tmp.txt)
	echo "" >> $@
	echo "Thanks," >> $@
	echo "AFLOW Dev" >> $@
#CO20190218 - STOP
###################################################################################################

###################################################################################################
#UNIT TESTS
.PHONY: check
check:
#CO20220906 - run one quick instance of aflow to read/write aflow.rc
	./aflow -v
	$(MAKE) check_proto
.PHONY: check_full
check_full:
	$(MAKE) check_aflow_exits
	$(MAKE) makefile
	$(MAKE) check_smith
	$(MAKE) check_proto
	$(MAKE) check_chull
	$(MAKE) check_aflowSG
	$(MAKE) check_AFLOWSYM
	$(MAKE) check_edata
	$(MAKE) check_aflow_xtal_finder
	$(MAKE) check_cce
	$(MAKE) check_pocc
	$(MAKE) unit_test
.PHONY: check_proto
check_proto:
	@echo "Performing 'aflow --proto' test"
	@if [ -z `which md5sum` ]; then echo "Cannot find binary md5sum - check skipped."; else \
		PROTO4SUM=$$(./aflow --proto=4:Bi:Tc | ./aflow --sprim | md5sum | cut -d ' ' -f1) && \
		[ "$${PROTO4SUM#*ffb9c5681045fb80f391e9a931f21d1}" != "$$PROTO4SUM" ] && \
		echo "'aflow --proto' check passed!" && exit 0 || \
		echo "'aflow --proto' test failed. Check compilation of aflow and aflow_data." && exit 1; \
	fi
#ME190313 - only execute check if md5sum is present
#	-if [ -z `which md5sum` ]; then echo "Cannot find binary md5sum - check skipped."; else ./aflow --proto=4:Bi:Tc | ./aflow --sprim | md5sum; echo "#ffb9c5681045fb80f391e9a931f21d1"; fi
#	./aflow --proto=4:Bi:Tc | ./aflow --sprim | md5sum
#	#ffb9c5681045fb80f391e9a931f21d1
#	./aflow --checki | grep FAIL
.PHONY: check_smith
check_smith:
	@if [ ! -f ./aflow ]; then echo "No local aflow binary found (./aflow)" && exit 1; fi
	./aflow --smith_test
.PHONY: check_chull
check_chull:
	@if [ ! -f ./aflow ]; then echo "No local aflow binary found (./aflow)" && exit 1; fi
	@rm -f $(COUNTER_FILE) && echo "STARTING CHULL" >> $(COUNTER_FILE)
	@ALLOY="MnPd" && COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running CHULL test $$COUNTER - text/json output files for $$ALLOY" | tee -a $(COUNTER_FILE) && COMMAND="./aflow --chull --alloy=$$ALLOY --o=j,t" RESPONSE=$$(eval "$$COMMAND") && \
		if [ -z "$${RESPONSE##*REST-API appears to be down*}" ]; then echo "CHULL did not return text/json output files, the REST-API may be down, please check '$$COMMAND'" && exit 0; fi && \
		if [ ! -f "aflow_$${ALLOY}_hull.json" ] || [ ! -f "aflow_$${ALLOY}_hull.txt" ]; then echo "CHULL did not return $$ALLOY text/json output files. Try: '$$COMMAND'" && exit 1; fi && \
		rm -f aflow_$${ALLOY}_hull.json aflow_$${ALLOY}_hull.txt
#travis does not have pdflatex
	@ALLOY="MnPd" && if !([ -n "$$TRAVIS" ] || [ "$$CI" = true ]); then COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running CHULL test $$COUNTER - PDF output for $$ALLOY" | tee -a $(COUNTER_FILE) && COMMAND="./aflow --chull --alloy=$$ALLOY --o=p" RESPONSE=$$(eval "$$COMMAND") && \
		if [ -z "$${RESPONSE##*REST-API appears to be down*}" ]; then echo "CHULL did not return a PDF, the REST-API may be down, please check '$$COMMAND'" && exit 0; fi && \
		if [ ! -f "aflow_$${ALLOY}_hull.pdf" ]; then echo "CHULL did not return a PDF for $$ALLOY. Try: '$$COMMAND'" && exit 1; fi && \
		rm -f aflow_$${ALLOY}_hull.pdf; fi
	@ALLOY="AgAuCd" && COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running CHULL test $$COUNTER - text/json output files for $$ALLOY" | tee -a $(COUNTER_FILE) && COMMAND="./aflow --chull --alloy=$$ALLOY --o=j,t" RESPONSE=$$(eval "$$COMMAND") && \
		if [ -z "$${RESPONSE##*REST-API appears to be down*}" ]; then echo "CHULL did not return text/json output files, the REST-API may be down, please check '$$COMMAND'" && exit 0; fi && \
		if [ ! -f "aflow_$${ALLOY}_hull.json" ] || [ ! -f "aflow_$${ALLOY}_hull.txt" ]; then echo "CHULL did not return $$ALLOY text/json output files. Try: '$$COMMAND'" && exit 1; fi && \
		rm -f aflow_$${ALLOY}_hull.json aflow_$${ALLOY}_hull.txt
#travis does not have pdflatex
	@ALLOY="AgAuCd" && if !([ -n "$$TRAVIS" ] || [ "$$CI" = true ]); then COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running CHULL test $$COUNTER - PDF output for $$ALLOY" | tee -a $(COUNTER_FILE) && COMMAND="./aflow --chull --alloy=$$ALLOY --o=p" RESPONSE=$$(eval "$$COMMAND") && \
		if [ -z "$${RESPONSE##*REST-API appears to be down*}" ]; then echo "CHULL did not return a PDF, the REST-API may be down, please check '$$COMMAND'" && exit 0; fi && \
		if [ ! -f "aflow_$${ALLOY}_hull.pdf" ]; then echo "CHULL did not return a PDF for $$ALLOY. Try: '$$COMMAND'" && exit 1; fi && \
		rm -f aflow_$${ALLOY}_hull.pdf; fi
	@ALLOY="AlCo" && COMPOUND="Al17Co12" && COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running CHULL test $$COUNTER - distance from hull of $$COMPOUND" | tee -a $(COUNTER_FILE) && ANSWER="242." && COMMAND="./aflow --chull --alloy=$$ALLOY --d2h=aflow:ef70a469857d3c47" RESPONSE=$$(eval "$$COMMAND --screen_only --o=t" | cut -d' ' -f2) && \
		if [ -z "$$RESPONSE" ]; then echo "CHULL did not return a distance, the REST-API may be down, please check '$$COMMAND'" && exit 0; fi && \
		if [ -n "$${RESPONSE##$$ANSWER*}" ]; then echo "CHULL did not return a distance of ~$$(echo $$ANSWER + 0 | bc) meV/atom for $$COMPOUND, found $$RESPONSE instead. Try: '$$COMMAND'" && exit 1; fi
	@ALLOY="MnPd" && COMPOUND="Mn2Pd3" && COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running CHULL test $$COUNTER - distance from hull of $$COMPOUND" | tee -a $(COUNTER_FILE) && ANSWER="62." && COMMAND="./aflow --chull --alloy=$$ALLOY --d2h=aflow:9001c322296fc162" RESPONSE=$$(eval "$$COMMAND --screen_only --o=t" | cut -d' ' -f2) && \
		if [ -z "$$RESPONSE" ]; then echo "CHULL did not return a distance, the REST-API may be down, please check '$$COMMAND'" && exit 0; fi && \
		if [ -n "$${RESPONSE##$$ANSWER*}" ]; then echo "CHULL did not return a distance of ~$$(echo $$ANSWER + 0 | bc) meV/atom for $$COMPOUND, found $$RESPONSE instead. Try: '$$COMMAND'" && exit 1; fi
	@ALLOY="TeZr" && COMPOUND="Te2Zr" && COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running CHULL test $$COUNTER - stability criterion of $$COMPOUND" | tee -a $(COUNTER_FILE) && ANSWER="48." && COMMAND="./aflow --chull --alloy=$$ALLOY --sc=aflow:570e80d4daa0a9d4" && RESPONSE=$$(eval "$$COMMAND --screen_only --o=t" | cut -d' ' -f2) && \
		if [ -z "$$RESPONSE" ]; then echo "CHULL did not return a stability criterion, the REST-API may be down, please check '$$COMMAND'" && exit 0; fi && \
		if [ -n "$${RESPONSE##$$ANSWER*}" ]; then echo "CHULL did not return a stability criterion of ~$$(echo $$ANSWER + 0 | bc) meV/atom for $$COMPOUND, found $$RESPONSE instead. Try: '$$COMMAND'" && exit 1; fi
	@ALLOY="MnPd" && COMPOUND="Mn3Pd5" && COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running CHULL test $$COUNTER - stability criterion of $$COMPOUND" | tee -a $(COUNTER_FILE) && ANSWER="13." && COMMAND="./aflow --chull --alloy=$$ALLOY --sc=aflow:4fc97e549ca17d4e" && RESPONSE=$$(eval "$$COMMAND --screen_only --o=t" | cut -d' ' -f2) && \
		if [ -z "$$RESPONSE" ]; then echo "CHULL did not return a stability criterion, the REST-API may be down, please check '$$COMMAND'" && exit 0; fi && \
		if [ -n "$${RESPONSE##$$ANSWER*}" ]; then echo "CHULL did not return a stability criterion of ~$$(echo $$ANSWER + 0 | bc) meV/atom for $$COMPOUND, found $$RESPONSE instead. Try: '$$COMMAND'" && exit 1; fi
	@ALLOY="MnPdPt" && COMPOUND="Mn2PdPt" && COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running CHULL test $$COUNTER - stability criterion of $$COMPOUND" | tee -a $(COUNTER_FILE) && ANSWER="32." && COMMAND="./aflow --chull --alloy=$$ALLOY --sc=aflow:fb9eaa58604ce774" && RESPONSE=$$(eval "$$COMMAND --screen_only --o=t" | cut -d' ' -f2) && \
		if [ -z "$$RESPONSE" ]; then echo "CHULL did not return a stability criterion, the REST-API may be down, please check '$$COMMAND'" && exit 0; fi && \
		if [ -n "$${RESPONSE##$$ANSWER*}" ]; then echo "CHULL did not return a stability criterion of ~$$(echo $$ANSWER + 0 | bc) meV/atom for $$COMPOUND, found $$RESPONSE instead. Try: '$$COMMAND'" && exit 1; fi
	@ALLOY="MnPd" && STOICH="x=0.25" && COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running CHULL test $$COUNTER - $$ALLOY hull energy at $$STOICH" | tee -a $(COUNTER_FILE) && ANSWER="-235." && COMMAND="./aflow --chull --alloy=$$ALLOY --hull_energy=0.25" && RESPONSE=$$(eval "$$COMMAND --screen_only --o=t" | cut -d' ' -f2) && \
		if [ -z "$$RESPONSE" ]; then echo "CHULL did not return a hull energy, the REST-API may be down, please check '$$COMMAND'" && exit 0; fi && \
		if [ -n "$${RESPONSE##$$ANSWER*}" ]; then echo "CHULL did not return a $$ALLOY hull energy of ~$$(echo $$ANSWER + 0 | bc) meV/atom at $$STOICH, found $$RESPONSE instead. Try: '$$COMMAND'" && exit 1; fi
	@ALLOY="CuZr" && STOICH="x=0.5" && COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running CHULL test $$COUNTER - $$ALLOY hull energy at $$STOICH" | tee -a $(COUNTER_FILE) && ANSWER="-158." && COMMAND="./aflow --chull --alloy=$$ALLOY --hull_energy=0.5" && RESPONSE=$$(eval "$$COMMAND --screen_only --o=t" | cut -d' ' -f2) && \
		if [ -z "$$RESPONSE" ]; then echo "CHULL did not return a hull energy, the REST-API may be down, please check '$$COMMAND'" && exit 0; fi && \
		if [ -n "$${RESPONSE##$$ANSWER*}" ]; then echo "CHULL did not return a $$ALLOY hull energy of ~$$(echo $$ANSWER + 0 | bc) meV/atom at $$STOICH, found $$RESPONSE instead. Try: '$$COMMAND'" && exit 1; fi
	@ALLOY="MnPdPt" && COMPOUND="Mn2PdPt" && COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running CHULL test $$COUNTER - N+1 enthalpy gain of $$COMPOUND" | tee -a $(COUNTER_FILE) && ANSWER="32." && COMMAND="./aflow --chull --alloy=$$ALLOY --n+1_enthalpy_gain=aflow:fb9eaa58604ce774" && RESPONSE=$$(eval "$$COMMAND --screen_only --o=t" | cut -d' ' -f2) && \
		if [ -z "$$RESPONSE" ]; then echo "CHULL did not return an enthalpy gain, the REST-API may be down, please check '$$COMMAND'" && exit 0; fi && \
		if [ -n "$${RESPONSE##$$ANSWER*}" ]; then echo "CHULL did not return an enthalpy gain of ~$$(echo $$ANSWER + 0 | bc) meV/atom for $$COMPOUND, found $$RESPONSE instead. Try: '$$COMMAND'" && exit 1; fi
	@echo "CHULL tests passed"
	@rm -f $(COUNTER_FILE)
.PHONY: check_pocc
check_pocc:
	@if [ ! -f ./aflow ]; then echo "No local aflow binary found (./aflow)" && exit 1; fi
	@rm -f $(COUNTER_FILE) && echo "STARTING POCC" >> $(COUNTER_FILE)
#Mg_(x)Zn_(1-x)O and ZnS_(1-x)Se_(x)
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt Mg0.75Zn0.25O" | tee -a $(COUNTER_FILE) && ANSWER=7 && COMMAND="./aflow --proto=201:Mg:O:Zn --pocc_params=P0-1xB_P1-0.75xA-0.25xC | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of wurtzite ZnS0.25Se0.75" | tee -a $(COUNTER_FILE) && ANSWER=7 && COMMAND="./aflow --proto=218:S:Se:Zn --pocc_params=S0-1xC_S1-0.75xB-0.25xA | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
#travis needs aflow-data in the path
	@if [ -n "$$TRAVIS" ] || [ "$$CI" = true ]; then export PATH=$$(pwd):$$PATH; fi && COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - aflow.in generator with different pocc_params" | tee -a $(COUNTER_FILE) && ANSWER=2 && AFLOW_BIN=$$(pwd)/aflow && TEMPDIR=$$(mktemp -d) && AFLOWDATA_DIR="./AFLOWDATA" && COMMAND="$$AFLOW_BIN --aflow_proto=AB_cF8_225_a_b.AB:C:Hf:N:Zr --pocc_params=S0-1xA_S1-0.5xB-0.5xD,S0-0.5xA-0.5xC_S1-0.5xB-0.5xD --generate_aflowin_only" && COUNT=$$(cd $$TEMPDIR && eval "$$COMMAND" >/dev/null 2>&1 && find $$AFLOWDATA_DIR -name aflow.in | wc -l | sed 's/^ *//') && rm -rf $$TEMPDIR && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER new aflow.in's for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@if [ -n "$$TRAVIS" ] || [ "$$CI" = true ]; then export PATH=$$(pwd):$$PATH; fi && COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - bcc/fcc Fe_(x)Cu_(1-x) aflow.in generator with different concentrations" | tee -a $(COUNTER_FILE) && ANSWER=22 && AFLOW_BIN=$$(pwd)/aflow && TEMPDIR=$$(mktemp -d) && AFLOWDATA_DIR="./AFLOWDATA" && COMMAND="$$AFLOW_BIN --aflow_proto=bcc,fcc:Cu:Fe --pocc_params=P0-0.1xA-0.9xB,P0-0.2xA-0.8xB,P0-0.25xA-0.75xB,P0-0.333xA-0.667xB,P0-0.4xA-0.6xB,P0-0.5xA-0.5xB,P0-0.6xA-0.4xB,P0-0.667xA-0.333xB,P0-0.75xA-0.25xB,P0-0.8xA-0.2xB,P0-0.9xA-0.1xB --generate_aflowin_only" && COUNT=$$(cd $$TEMPDIR && eval "$$COMMAND" >/dev/null 2>&1 && find $$AFLOWDATA_DIR -name aflow.in | wc -l | sed 's/^ *//') && rm -rf $$TEMPDIR && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER new aflow.in's for '$$COMMAND', found $$COUNT instead" && exit 1; fi
#bcc Fe_(x)Cu_(1-x)
	@if [ -n "$$TRAVIS" ] || [ "$$CI" = true ]; then export PATH=$$(pwd):$$PATH; fi && COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - bcc Fe0.1Cu0.9 aflow.in generator" | tee -a $(COUNTER_FILE) && ANSWER=18 && AFLOW_BIN=$$(pwd)/aflow && TEMPDIR=$$(mktemp -d) && AFLOWDATA_DIR="./AFLOWDATA/Cu_pvFe_pv/bcc:POCC_P0-0.1xA-0.9xB" && COMMAND="$$AFLOW_BIN --aflow_proto=bcc:Cu:Fe --pocc_params=P0-0.1xA-0.9xB --generate_aflowin_only && cd $$AFLOWDATA_DIR && $$AFLOW_BIN --run --generate_aflowin_only" && COUNT=$$(cd $$TEMPDIR && eval "$$COMMAND" >/dev/null 2>&1 && find ARUN.POCC_* -name aflow.in | wc -l | sed 's/^ *//') && rm -rf $$TEMPDIR && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER new aflow.in's for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@if [ -n "$$TRAVIS" ] || [ "$$CI" = true ]; then export PATH=$$(pwd):$$PATH; fi && COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - bcc Fe0.2Cu0.8 aflow.in generator" | tee -a $(COUNTER_FILE) && ANSWER=5 && AFLOW_BIN=$$(pwd)/aflow && TEMPDIR=$$(mktemp -d) && AFLOWDATA_DIR="./AFLOWDATA/Cu_pvFe_pv/bcc:POCC_P0-0.2xA-0.8xB" && COMMAND="$$AFLOW_BIN --aflow_proto=bcc:Cu:Fe --pocc_params=P0-0.2xA-0.8xB --generate_aflowin_only && cd $$AFLOWDATA_DIR && $$AFLOW_BIN --run --generate_aflowin_only" && COUNT=$$(cd $$TEMPDIR && eval "$$COMMAND" >/dev/null 2>&1 && find ARUN.POCC_* -name aflow.in | wc -l | sed 's/^ *//') && rm -rf $$TEMPDIR && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER new aflow.in's for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@if [ -n "$$TRAVIS" ] || [ "$$CI" = true ]; then export PATH=$$(pwd):$$PATH; fi && COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - bcc Fe0.25Cu0.75 aflow.in generator" | tee -a $(COUNTER_FILE) && ANSWER=7 && AFLOW_BIN=$$(pwd)/aflow && TEMPDIR=$$(mktemp -d) && AFLOWDATA_DIR="./AFLOWDATA/Cu_pvFe_pv/bcc:POCC_P0-0.25xA-0.75xB" && COMMAND="$$AFLOW_BIN --aflow_proto=bcc:Cu:Fe --pocc_params=P0-0.25xA-0.75xB --generate_aflowin_only && cd $$AFLOWDATA_DIR && $$AFLOW_BIN --run --generate_aflowin_only" && COUNT=$$(cd $$TEMPDIR && eval "$$COMMAND" >/dev/null 2>&1 && find ARUN.POCC_* -name aflow.in | wc -l | sed 's/^ *//') && rm -rf $$TEMPDIR && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER new aflow.in's for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@if [ -n "$$TRAVIS" ] || [ "$$CI" = true ]; then export PATH=$$(pwd):$$PATH; fi && COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - bcc Fe0.333Cu0.667 aflow.in generator" | tee -a $(COUNTER_FILE) && ANSWER=3 && AFLOW_BIN=$$(pwd)/aflow && TEMPDIR=$$(mktemp -d) && AFLOWDATA_DIR="./AFLOWDATA/Cu_pvFe_pv/bcc:POCC_P0-0.333xA-0.667xB" && COMMAND="$$AFLOW_BIN --aflow_proto=bcc:Cu:Fe --pocc_params=P0-0.333xA-0.667xB --generate_aflowin_only && cd $$AFLOWDATA_DIR && $$AFLOW_BIN --run --generate_aflowin_only" && COUNT=$$(cd $$TEMPDIR && eval "$$COMMAND" >/dev/null 2>&1 && find ARUN.POCC_* -name aflow.in | wc -l | sed 's/^ *//') && rm -rf $$TEMPDIR && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER new aflow.in's for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@if [ -n "$$TRAVIS" ] || [ "$$CI" = true ]; then export PATH=$$(pwd):$$PATH; fi && COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - bcc Fe0.4Cu0.6 aflow.in generator" | tee -a $(COUNTER_FILE) && ANSWER=9 && AFLOW_BIN=$$(pwd)/aflow && TEMPDIR=$$(mktemp -d) && AFLOWDATA_DIR="./AFLOWDATA/Cu_pvFe_pv/bcc:POCC_P0-0.4xA-0.6xB" && COMMAND="$$AFLOW_BIN --aflow_proto=bcc:Cu:Fe --pocc_params=P0-0.4xA-0.6xB --generate_aflowin_only && cd $$AFLOWDATA_DIR && $$AFLOW_BIN --run --generate_aflowin_only" && COUNT=$$(cd $$TEMPDIR && eval "$$COMMAND" >/dev/null 2>&1 && find ARUN.POCC_* -name aflow.in | wc -l | sed 's/^ *//') && rm -rf $$TEMPDIR && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER new aflow.in's for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@if [ -n "$$TRAVIS" ] || [ "$$CI" = true ]; then export PATH=$$(pwd):$$PATH; fi && COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - bcc Fe0.5Cu0.5 aflow.in generator" | tee -a $(COUNTER_FILE) && ANSWER=2 && AFLOW_BIN=$$(pwd)/aflow && TEMPDIR=$$(mktemp -d) && AFLOWDATA_DIR="./AFLOWDATA/Cu_pvFe_pv/bcc:POCC_P0-0.5xA-0.5xB" && COMMAND="$$AFLOW_BIN --aflow_proto=bcc:Cu:Fe --pocc_params=P0-0.5xA-0.5xB --generate_aflowin_only && cd $$AFLOWDATA_DIR && $$AFLOW_BIN --run --generate_aflowin_only" && COUNT=$$(cd $$TEMPDIR && eval "$$COMMAND" >/dev/null 2>&1 && find ARUN.POCC_* -name aflow.in | wc -l | sed 's/^ *//') && rm -rf $$TEMPDIR && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER new aflow.in's for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@if [ -n "$$TRAVIS" ] || [ "$$CI" = true ]; then export PATH=$$(pwd):$$PATH; fi && COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - bcc Fe0.6Cu0.4 aflow.in generator" | tee -a $(COUNTER_FILE) && ANSWER=9 && AFLOW_BIN=$$(pwd)/aflow && TEMPDIR=$$(mktemp -d) && AFLOWDATA_DIR="./AFLOWDATA/Cu_pvFe_pv/bcc:POCC_P0-0.6xA-0.4xB" && COMMAND="$$AFLOW_BIN --aflow_proto=bcc:Cu:Fe --pocc_params=P0-0.6xA-0.4xB --generate_aflowin_only && cd $$AFLOWDATA_DIR && $$AFLOW_BIN --run --generate_aflowin_only" && COUNT=$$(cd $$TEMPDIR && eval "$$COMMAND" >/dev/null 2>&1 && find ARUN.POCC_* -name aflow.in | wc -l | sed 's/^ *//') && rm -rf $$TEMPDIR && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER new aflow.in's for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@if [ -n "$$TRAVIS" ] || [ "$$CI" = true ]; then export PATH=$$(pwd):$$PATH; fi && COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - bcc Fe0.667Cu0.333 aflow.in generator" | tee -a $(COUNTER_FILE) && ANSWER=3 && AFLOW_BIN=$$(pwd)/aflow && TEMPDIR=$$(mktemp -d) && AFLOWDATA_DIR="./AFLOWDATA/Cu_pvFe_pv/bcc:POCC_P0-0.667xA-0.333xB" && COMMAND="$$AFLOW_BIN --aflow_proto=bcc:Cu:Fe --pocc_params=P0-0.667xA-0.333xB --generate_aflowin_only && cd $$AFLOWDATA_DIR && $$AFLOW_BIN --run --generate_aflowin_only" && COUNT=$$(cd $$TEMPDIR && eval "$$COMMAND" >/dev/null 2>&1 && find ARUN.POCC_* -name aflow.in | wc -l | sed 's/^ *//') && rm -rf $$TEMPDIR && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER new aflow.in's for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@if [ -n "$$TRAVIS" ] || [ "$$CI" = true ]; then export PATH=$$(pwd):$$PATH; fi && COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - bcc Fe0.75Cu0.25 aflow.in generator" | tee -a $(COUNTER_FILE) && ANSWER=7 && AFLOW_BIN=$$(pwd)/aflow && TEMPDIR=$$(mktemp -d) && AFLOWDATA_DIR="./AFLOWDATA/Cu_pvFe_pv/bcc:POCC_P0-0.75xA-0.25xB" && COMMAND="$$AFLOW_BIN --aflow_proto=bcc:Cu:Fe --pocc_params=P0-0.75xA-0.25xB --generate_aflowin_only && cd $$AFLOWDATA_DIR && $$AFLOW_BIN --run --generate_aflowin_only" && COUNT=$$(cd $$TEMPDIR && eval "$$COMMAND" >/dev/null 2>&1 && find ARUN.POCC_* -name aflow.in | wc -l | sed 's/^ *//') && rm -rf $$TEMPDIR && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER new aflow.in's for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@if [ -n "$$TRAVIS" ] || [ "$$CI" = true ]; then export PATH=$$(pwd):$$PATH; fi && COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - bcc Fe0.8Cu0.2 aflow.in generator" | tee -a $(COUNTER_FILE) && ANSWER=5 && AFLOW_BIN=$$(pwd)/aflow && TEMPDIR=$$(mktemp -d) && AFLOWDATA_DIR="./AFLOWDATA/Cu_pvFe_pv/bcc:POCC_P0-0.8xA-0.2xB" && COMMAND="$$AFLOW_BIN --aflow_proto=bcc:Cu:Fe --pocc_params=P0-0.8xA-0.2xB --generate_aflowin_only && cd $$AFLOWDATA_DIR && $$AFLOW_BIN --run --generate_aflowin_only" && COUNT=$$(cd $$TEMPDIR && eval "$$COMMAND" >/dev/null 2>&1 && find ARUN.POCC_* -name aflow.in | wc -l | sed 's/^ *//') && rm -rf $$TEMPDIR && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER new aflow.in's for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@if [ -n "$$TRAVIS" ] || [ "$$CI" = true ]; then export PATH=$$(pwd):$$PATH; fi && COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - bcc Fe0.9Cu0.1 aflow.in generator" | tee -a $(COUNTER_FILE) && ANSWER=18 && AFLOW_BIN=$$(pwd)/aflow && TEMPDIR=$$(mktemp -d) && AFLOWDATA_DIR="./AFLOWDATA/Cu_pvFe_pv/bcc:POCC_P0-0.9xA-0.1xB" && COMMAND="$$AFLOW_BIN --aflow_proto=bcc:Cu:Fe --pocc_params=P0-0.9xA-0.1xB --generate_aflowin_only && cd $$AFLOWDATA_DIR && $$AFLOW_BIN --run --generate_aflowin_only" && COUNT=$$(cd $$TEMPDIR && eval "$$COMMAND" >/dev/null 2>&1 && find ARUN.POCC_* -name aflow.in | wc -l | sed 's/^ *//') && rm -rf $$TEMPDIR && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER new aflow.in's for '$$COMMAND', found $$COUNT instead" && exit 1; fi
#fcc Fe_(x)_Cu_(1-x)
	@if [ -n "$$TRAVIS" ] || [ "$$CI" = true ]; then export PATH=$$(pwd):$$PATH; fi && COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - fcc Fe0.1Cu0.9 aflow.in generator" | tee -a $(COUNTER_FILE) && ANSWER=18 && AFLOW_BIN=$$(pwd)/aflow && TEMPDIR=$$(mktemp -d) && AFLOWDATA_DIR="./AFLOWDATA/Cu_pvFe_pv/fcc:POCC_P0-0.1xA-0.9xB" && COMMAND="$$AFLOW_BIN --aflow_proto=fcc:Cu:Fe --pocc_params=P0-0.1xA-0.9xB --generate_aflowin_only && cd $$AFLOWDATA_DIR && $$AFLOW_BIN --run --generate_aflowin_only" && COUNT=$$(cd $$TEMPDIR && eval "$$COMMAND" >/dev/null 2>&1 && find ARUN.POCC_* -name aflow.in | wc -l | sed 's/^ *//') && rm -rf $$TEMPDIR && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER new aflow.in's for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@if [ -n "$$TRAVIS" ] || [ "$$CI" = true ]; then export PATH=$$(pwd):$$PATH; fi && COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - fcc Fe0.2Cu0.8 aflow.in generator" | tee -a $(COUNTER_FILE) && ANSWER=5 && AFLOW_BIN=$$(pwd)/aflow && TEMPDIR=$$(mktemp -d) && AFLOWDATA_DIR="./AFLOWDATA/Cu_pvFe_pv/fcc:POCC_P0-0.2xA-0.8xB" && COMMAND="$$AFLOW_BIN --aflow_proto=fcc:Cu:Fe --pocc_params=P0-0.2xA-0.8xB --generate_aflowin_only && cd $$AFLOWDATA_DIR && $$AFLOW_BIN --run --generate_aflowin_only" && COUNT=$$(cd $$TEMPDIR && eval "$$COMMAND" >/dev/null 2>&1 && find ARUN.POCC_* -name aflow.in | wc -l | sed 's/^ *//') && rm -rf $$TEMPDIR && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER new aflow.in's for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@if [ -n "$$TRAVIS" ] || [ "$$CI" = true ]; then export PATH=$$(pwd):$$PATH; fi && COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - fcc Fe0.25Cu0.75 aflow.in generator" | tee -a $(COUNTER_FILE) && ANSWER=7 && AFLOW_BIN=$$(pwd)/aflow && TEMPDIR=$$(mktemp -d) && AFLOWDATA_DIR="./AFLOWDATA/Cu_pvFe_pv/fcc:POCC_P0-0.25xA-0.75xB" && COMMAND="$$AFLOW_BIN --aflow_proto=fcc:Cu:Fe --pocc_params=P0-0.25xA-0.75xB --generate_aflowin_only && cd $$AFLOWDATA_DIR && $$AFLOW_BIN --run --generate_aflowin_only" && COUNT=$$(cd $$TEMPDIR && eval "$$COMMAND" >/dev/null 2>&1 && find ARUN.POCC_* -name aflow.in | wc -l | sed 's/^ *//') && rm -rf $$TEMPDIR && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER new aflow.in's for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@if [ -n "$$TRAVIS" ] || [ "$$CI" = true ]; then export PATH=$$(pwd):$$PATH; fi && COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - fcc Fe0.333Cu0.667 aflow.in generator" | tee -a $(COUNTER_FILE) && ANSWER=3 && AFLOW_BIN=$$(pwd)/aflow && TEMPDIR=$$(mktemp -d) && AFLOWDATA_DIR="./AFLOWDATA/Cu_pvFe_pv/fcc:POCC_P0-0.333xA-0.667xB" && COMMAND="$$AFLOW_BIN --aflow_proto=fcc:Cu:Fe --pocc_params=P0-0.333xA-0.667xB --generate_aflowin_only && cd $$AFLOWDATA_DIR && $$AFLOW_BIN --run --generate_aflowin_only" && COUNT=$$(cd $$TEMPDIR && eval "$$COMMAND" >/dev/null 2>&1 && find ARUN.POCC_* -name aflow.in | wc -l | sed 's/^ *//') && rm -rf $$TEMPDIR && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER new aflow.in's for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@if [ -n "$$TRAVIS" ] || [ "$$CI" = true ]; then export PATH=$$(pwd):$$PATH; fi && COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - fcc Fe0.4Cu0.6 aflow.in generator" | tee -a $(COUNTER_FILE) && ANSWER=9 && AFLOW_BIN=$$(pwd)/aflow && TEMPDIR=$$(mktemp -d) && AFLOWDATA_DIR="./AFLOWDATA/Cu_pvFe_pv/fcc:POCC_P0-0.4xA-0.6xB" && COMMAND="$$AFLOW_BIN --aflow_proto=fcc:Cu:Fe --pocc_params=P0-0.4xA-0.6xB --generate_aflowin_only && cd $$AFLOWDATA_DIR && $$AFLOW_BIN --run --generate_aflowin_only" && COUNT=$$(cd $$TEMPDIR && eval "$$COMMAND" >/dev/null 2>&1 && find ARUN.POCC_* -name aflow.in | wc -l | sed 's/^ *//') && rm -rf $$TEMPDIR && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER new aflow.in's for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@if [ -n "$$TRAVIS" ] || [ "$$CI" = true ]; then export PATH=$$(pwd):$$PATH; fi && COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - fcc Fe0.5Cu0.5 aflow.in generator" | tee -a $(COUNTER_FILE) && ANSWER=2 && AFLOW_BIN=$$(pwd)/aflow && TEMPDIR=$$(mktemp -d) && AFLOWDATA_DIR="./AFLOWDATA/Cu_pvFe_pv/fcc:POCC_P0-0.5xA-0.5xB" && COMMAND="$$AFLOW_BIN --aflow_proto=fcc:Cu:Fe --pocc_params=P0-0.5xA-0.5xB --generate_aflowin_only && cd $$AFLOWDATA_DIR && $$AFLOW_BIN --run --generate_aflowin_only" && COUNT=$$(cd $$TEMPDIR && eval "$$COMMAND" >/dev/null 2>&1 && find ARUN.POCC_* -name aflow.in | wc -l | sed 's/^ *//') && rm -rf $$TEMPDIR && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER new aflow.in's for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@if [ -n "$$TRAVIS" ] || [ "$$CI" = true ]; then export PATH=$$(pwd):$$PATH; fi && COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - fcc Fe0.6Cu0.4 aflow.in generator" | tee -a $(COUNTER_FILE) && ANSWER=9 && AFLOW_BIN=$$(pwd)/aflow && TEMPDIR=$$(mktemp -d) && AFLOWDATA_DIR="./AFLOWDATA/Cu_pvFe_pv/fcc:POCC_P0-0.6xA-0.4xB" && COMMAND="$$AFLOW_BIN --aflow_proto=fcc:Cu:Fe --pocc_params=P0-0.6xA-0.4xB --generate_aflowin_only && cd $$AFLOWDATA_DIR && $$AFLOW_BIN --run --generate_aflowin_only" && COUNT=$$(cd $$TEMPDIR && eval "$$COMMAND" >/dev/null 2>&1 && find ARUN.POCC_* -name aflow.in | wc -l | sed 's/^ *//') && rm -rf $$TEMPDIR && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER new aflow.in's for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@if [ -n "$$TRAVIS" ] || [ "$$CI" = true ]; then export PATH=$$(pwd):$$PATH; fi && COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - fcc Fe0.667Cu0.333 aflow.in generator" | tee -a $(COUNTER_FILE) && ANSWER=3 && AFLOW_BIN=$$(pwd)/aflow && TEMPDIR=$$(mktemp -d) && AFLOWDATA_DIR="./AFLOWDATA/Cu_pvFe_pv/fcc:POCC_P0-0.667xA-0.333xB" && COMMAND="$$AFLOW_BIN --aflow_proto=fcc:Cu:Fe --pocc_params=P0-0.667xA-0.333xB --generate_aflowin_only && cd $$AFLOWDATA_DIR && $$AFLOW_BIN --run --generate_aflowin_only" && COUNT=$$(cd $$TEMPDIR && eval "$$COMMAND" >/dev/null 2>&1 && find ARUN.POCC_* -name aflow.in | wc -l | sed 's/^ *//') && rm -rf $$TEMPDIR && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER new aflow.in's for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@if [ -n "$$TRAVIS" ] || [ "$$CI" = true ]; then export PATH=$$(pwd):$$PATH; fi && COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - fcc Fe0.75Cu0.25 aflow.in generator" | tee -a $(COUNTER_FILE) && ANSWER=7 && AFLOW_BIN=$$(pwd)/aflow && TEMPDIR=$$(mktemp -d) && AFLOWDATA_DIR="./AFLOWDATA/Cu_pvFe_pv/fcc:POCC_P0-0.75xA-0.25xB" && COMMAND="$$AFLOW_BIN --aflow_proto=fcc:Cu:Fe --pocc_params=P0-0.75xA-0.25xB --generate_aflowin_only && cd $$AFLOWDATA_DIR && $$AFLOW_BIN --run --generate_aflowin_only" && COUNT=$$(cd $$TEMPDIR && eval "$$COMMAND" >/dev/null 2>&1 && find ARUN.POCC_* -name aflow.in | wc -l | sed 's/^ *//') && rm -rf $$TEMPDIR && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER new aflow.in's for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@if [ -n "$$TRAVIS" ] || [ "$$CI" = true ]; then export PATH=$$(pwd):$$PATH; fi && COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - fcc Fe0.8Cu0.2 aflow.in generator" | tee -a $(COUNTER_FILE) && ANSWER=5 && AFLOW_BIN=$$(pwd)/aflow && TEMPDIR=$$(mktemp -d) && AFLOWDATA_DIR="./AFLOWDATA/Cu_pvFe_pv/fcc:POCC_P0-0.8xA-0.2xB" && COMMAND="$$AFLOW_BIN --aflow_proto=fcc:Cu:Fe --pocc_params=P0-0.8xA-0.2xB --generate_aflowin_only && cd $$AFLOWDATA_DIR && $$AFLOW_BIN --run --generate_aflowin_only" && COUNT=$$(cd $$TEMPDIR && eval "$$COMMAND" >/dev/null 2>&1 && find ARUN.POCC_* -name aflow.in | wc -l | sed 's/^ *//') && rm -rf $$TEMPDIR && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER new aflow.in's for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@if [ -n "$$TRAVIS" ] || [ "$$CI" = true ]; then export PATH=$$(pwd):$$PATH; fi && COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - fcc Fe0.9Cu0.1 aflow.in generator" | tee -a $(COUNTER_FILE) && ANSWER=18 && AFLOW_BIN=$$(pwd)/aflow && TEMPDIR=$$(mktemp -d) && AFLOWDATA_DIR="./AFLOWDATA/Cu_pvFe_pv/fcc:POCC_P0-0.9xA-0.1xB" && COMMAND="$$AFLOW_BIN --aflow_proto=fcc:Cu:Fe --pocc_params=P0-0.9xA-0.1xB --generate_aflowin_only && cd $$AFLOWDATA_DIR && $$AFLOW_BIN --run --generate_aflowin_only" && COUNT=$$(cd $$TEMPDIR && eval "$$COMMAND" >/dev/null 2>&1 && find ARUN.POCC_* -name aflow.in | wc -l | sed 's/^ *//') && rm -rf $$TEMPDIR && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER new aflow.in's for '$$COMMAND', found $$COUNT instead" && exit 1; fi
#testing P vs. S
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - optimized HNF value of cementite C5MoNbTaVW" | tee -a $(COUNTER_FILE) && ANSWER=5 && COMMAND="./aflow --proto=AB3_oP16_62_c_cd-001.AB:C:Mo:Nb:Ta:V:W --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --hnf" && COUNT=$$(eval "$$COMMAND" | grep 'Optimized HNF value = ' | tail -1 | awk -F 'value = ' '{print $$2}' | cut -d ' ' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return an HNF value of $$ANSWER for '$$COMMAND', found $$COUNT instead" && exit 1; fi
#5-metal carbides
	@if [ -n "$$TRAVIS" ] || [ "$$CI" = true ]; then export PATH=$$(pwd):$$PATH; fi && COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - 5-metal rocksalt refractory carbide aflow.in generator" | tee -a $(COUNTER_FILE) && ANSWER=56 && AFLOW_BIN=$$(pwd)/aflow && TEMPDIR=$$(mktemp -d) && AFLOWDATA_DIR="./AFLOWDATA" && COMMAND="$$AFLOW_BIN --aflow_proto=AB_cF8_225_a_b.AB:C:Hf,Mo,Nb,Ta,Ti,V,W,Zr:Hf,Mo,Nb,Ta,Ti,V,W,Zr:Hf,Mo,Nb,Ta,Ti,V,W,Zr:Hf,Mo,Nb,Ta,Ti,V,W,Zr:Hf,Mo,Nb,Ta,Ti,V,W,Zr --pocc_params=P0-1xA_P1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF --generate_aflowin_only" && COUNT=$$(cd $$TEMPDIR && eval "$$COMMAND" >/dev/null 2>&1 && find $$AFLOWDATA_DIR -name aflow.in | wc -l | sed 's/^ *//') && rm -rf $$TEMPDIR && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER new aflow.in's for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - optimized HNF value of rocksalt C2HfNb" | tee -a $(COUNTER_FILE) && ANSWER=2 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Nb --pocc_params=S0-1xA_S1-0.5xB-0.5xC | ./aflow --hnf" && COUNT=$$(eval "$$COMMAND" | grep 'Optimized HNF value = ' | tail -1 | awk -F 'value = ' '{print $$2}' | cut -d ' ' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return an HNF value of $$ANSWER for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C2HfNb" | tee -a $(COUNTER_FILE) && ANSWER=2 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Nb --pocc_params=S0-1xA_S1-0.5xB-0.5xC | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return an HNF value of $$ANSWER for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C5HfMoNbTaTi" | tee -a $(COUNTER_FILE) && ANSWER=49 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Mo:Nb:Ta:Ti --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C5HfMoNbTaV" | tee -a $(COUNTER_FILE) && ANSWER=49 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Mo:Nb:Ta:V --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C5HfMoNbTaW" | tee -a $(COUNTER_FILE) && ANSWER=49 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Mo:Nb:Ta:W --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C5HfMoNbTaZr" | tee -a $(COUNTER_FILE) && ANSWER=49 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Mo:Nb:Ta:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C5HfMoNbTiV" | tee -a $(COUNTER_FILE) && ANSWER=49 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Mo:Nb:Ti:V --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C5HfMoNbTiW" | tee -a $(COUNTER_FILE) && ANSWER=49 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Mo:Nb:Ti:W --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C5HfMoNbTiZr" | tee -a $(COUNTER_FILE) && ANSWER=49 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Mo:Nb:Ti:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C5HfMoNbVW" | tee -a $(COUNTER_FILE) && ANSWER=49 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Mo:Nb:V:W --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C5HfMoNbVZr" | tee -a $(COUNTER_FILE) && ANSWER=49 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Mo:Nb:V:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C5HfMoNbWZr" | tee -a $(COUNTER_FILE) && ANSWER=49 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Mo:Nb:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C5HfMoTaTiV" | tee -a $(COUNTER_FILE) && ANSWER=49 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Mo:Ta:Ti:V --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C5HfMoTaTiW" | tee -a $(COUNTER_FILE) && ANSWER=49 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Mo:Ta:Ti:W --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C5HfMoTaTiZr" | tee -a $(COUNTER_FILE) && ANSWER=49 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Mo:Ta:Ti:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C5HfMoTaVW" | tee -a $(COUNTER_FILE) && ANSWER=49 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Mo:Ta:V:W --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C5HfMoTaVZr" | tee -a $(COUNTER_FILE) && ANSWER=49 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Mo:Ta:V:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C5HfMoTaWZr" | tee -a $(COUNTER_FILE) && ANSWER=49 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Mo:Ta:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C5HfMoTiVW" | tee -a $(COUNTER_FILE) && ANSWER=49 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Mo:Ti:V:W --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C5HfMoTiVZr" | tee -a $(COUNTER_FILE) && ANSWER=49 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Mo:Ti:V:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C5HfMoTiWZr" | tee -a $(COUNTER_FILE) && ANSWER=49 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C5HfMoVWZr" | tee -a $(COUNTER_FILE) && ANSWER=49 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Mo:V:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C5HfNbTaTiV" | tee -a $(COUNTER_FILE) && ANSWER=49 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Nb:Ta:Ti:V --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C5HfNbTaTiW" | tee -a $(COUNTER_FILE) && ANSWER=49 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Nb:Ta:Ti:W --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C5HfNbTaTiZr" | tee -a $(COUNTER_FILE) && ANSWER=49 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Nb:Ta:Ti:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C5HfNbTaVW" | tee -a $(COUNTER_FILE) && ANSWER=49 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Nb:Ta:V:W --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C5HfNbTaVZr" | tee -a $(COUNTER_FILE) && ANSWER=49 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Nb:Ta:V:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C5HfNbTaWZr" | tee -a $(COUNTER_FILE) && ANSWER=49 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Nb:Ta:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C5HfNbTiVW" | tee -a $(COUNTER_FILE) && ANSWER=49 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Nb:Ti:V:W --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C5HfNbTiVZr" | tee -a $(COUNTER_FILE) && ANSWER=49 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Nb:Ti:V:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C5HfNbTiWZr" | tee -a $(COUNTER_FILE) && ANSWER=49 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Nb:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C5HfNbVWZr" | tee -a $(COUNTER_FILE) && ANSWER=49 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Nb:V:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C5HfTaTiVW" | tee -a $(COUNTER_FILE) && ANSWER=49 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Ta:Ti:V:W --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C5HfTaTiVZr" | tee -a $(COUNTER_FILE) && ANSWER=49 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Ta:Ti:V:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C5HfTaTiWZr" | tee -a $(COUNTER_FILE) && ANSWER=49 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Ta:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C5HfTaVWZr" | tee -a $(COUNTER_FILE) && ANSWER=49 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Ta:V:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C5HfTiVWZr" | tee -a $(COUNTER_FILE) && ANSWER=49 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Ti:V:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C5MoNbTaTiV" | tee -a $(COUNTER_FILE) && ANSWER=49 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Mo:Nb:Ta:Ti:V --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C5MoNbTaTiW" | tee -a $(COUNTER_FILE) && ANSWER=49 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Mo:Nb:Ta:Ti:W --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C5MoNbTaTiZr" | tee -a $(COUNTER_FILE) && ANSWER=49 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Mo:Nb:Ta:Ti:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - optimized HNF value of rocksalt C5MoNbTaVW" | tee -a $(COUNTER_FILE) && ANSWER=5 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Mo:Nb:Ta:V:W --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --hnf" && COUNT=$$(eval "$$COMMAND" | grep 'Optimized HNF value = ' | tail -1 | awk -F 'value = ' '{print $$2}' | cut -d ' ' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return an HNF value of $$ANSWER for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C5MoNbTaVW" | tee -a $(COUNTER_FILE) && ANSWER=49 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Mo:Nb:Ta:V:W --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@if [ -n "$$TRAVIS" ] || [ "$$CI" = true ]; then export PATH=$$(pwd):$$PATH; fi && COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - aflow.in generator of rocksalt C5MoNbTaVW" | tee -a $(COUNTER_FILE) && ANSWER=49 && AFLOW_BIN=$$(pwd)/aflow && TEMPDIR=$$(mktemp -d) && AFLOWDATA_DIR="./AFLOWDATA/CMo_pvNb_svTa_pvV_svW_pv:PAW_PBE/AB_cF8_225_a_b.AB:POCC_P0-1xA_P1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF" && COMMAND="$$AFLOW_BIN --aflow_proto=AB_cF8_225_a_b.AB:C:Mo:Nb:Ta:V:W --pocc_params=P0-1xA_P1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF --generate_aflowin_only && cd $$AFLOWDATA_DIR && $$AFLOW_BIN --run --generate_aflowin_only" && COUNT=$$(cd $$TEMPDIR && eval "$$COMMAND" >/dev/null && find ARUN.POCC_* -name aflow.in | wc -l | sed 's/^ *//') && rm -rf $$TEMPDIR && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER new aflow.in's for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C5MoNbTaVZr" | tee -a $(COUNTER_FILE) && ANSWER=49 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Mo:Nb:Ta:V:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C5MoNbTaWZr" | tee -a $(COUNTER_FILE) && ANSWER=49 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Mo:Nb:Ta:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C5MoNbTiVW" | tee -a $(COUNTER_FILE) && ANSWER=49 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Mo:Nb:Ti:V:W --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C5MoNbTiVZr" | tee -a $(COUNTER_FILE) && ANSWER=49 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Mo:Nb:Ti:V:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C5MoNbTiWZr" | tee -a $(COUNTER_FILE) && ANSWER=49 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Mo:Nb:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C5MoNbVWZr" | tee -a $(COUNTER_FILE) && ANSWER=49 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Mo:Nb:V:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C5MoTaTiVW" | tee -a $(COUNTER_FILE) && ANSWER=49 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Mo:Ta:Ti:V:W --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C5MoTaTiVZr" | tee -a $(COUNTER_FILE) && ANSWER=49 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Mo:Ta:Ti:V:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C5MoTaTiWZr" | tee -a $(COUNTER_FILE) && ANSWER=49 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Mo:Ta:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C5MoTaVWZr" | tee -a $(COUNTER_FILE) && ANSWER=49 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Mo:Ta:V:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C5MoTiVWZr" | tee -a $(COUNTER_FILE) && ANSWER=49 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Mo:Ti:V:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C5NbTaTiVW" | tee -a $(COUNTER_FILE) && ANSWER=49 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Nb:Ta:Ti:V:W --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C5NbTaTiVZr" | tee -a $(COUNTER_FILE) && ANSWER=49 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Nb:Ta:Ti:V:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C5NbTaTiWZr" | tee -a $(COUNTER_FILE) && ANSWER=49 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Nb:Ta:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C5NbTaVWZr" | tee -a $(COUNTER_FILE) && ANSWER=49 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Nb:Ta:V:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C5NbTiVWZr" | tee -a $(COUNTER_FILE) && ANSWER=49 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Nb:Ti:V:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C5TaTiVWZr" | tee -a $(COUNTER_FILE) && ANSWER=49 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Ta:Ti:V:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C5CrMoTaVW" | tee -a $(COUNTER_FILE) && ANSWER=49 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ta:V:W --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C5CrMoNbTaW" | tee -a $(COUNTER_FILE) && ANSWER=49 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Nb:Ta:W --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C5CrMoNbVW" | tee -a $(COUNTER_FILE) && ANSWER=49 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Nb:V:W --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C5CrHfTaWZr" | tee -a $(COUNTER_FILE) && ANSWER=49 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Hf:Ta:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C5CrHfMoTiW" | tee -a $(COUNTER_FILE) && ANSWER=49 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Hf:Mo:Ti:W --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C5CrMoTiWZr" | tee -a $(COUNTER_FILE) && ANSWER=49 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
#5-metal borides
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt B10CrHfMoNbTa" | tee -a $(COUNTER_FILE) && ANSWER=84 && COMMAND="./aflow --proto=AB2_hP3_191_a_d-001.BA:B:Cr:Hf:Mo:Nb:Ta --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
#3-metal carbonitrides
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C3N3Hf2Ta2Ti2" | tee -a $(COUNTER_FILE) && ANSWER=1355 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:N:Ta:Ti --pocc_params=S0-0.5xA-0.5xC_S1-0.333xB-0.333xD-0.333xE | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C3N3Ta2Ti2Zr2" | tee -a $(COUNTER_FILE) && ANSWER=1355 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:N:Ta:Ti:Zr --pocc_params=S0-0.5xA-0.5xB_S1-0.333xC-0.333xD-0.333xE | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C3N3Nb2Ta2Ti2" | tee -a $(COUNTER_FILE) && ANSWER=1355 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:N:Nb:Ta:Ti --pocc_params=S0-0.5xA-0.5xB_S1-0.333xC-0.333xD-0.333xE | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt C3N3Hf2Nb2Zr2" | tee -a $(COUNTER_FILE) && ANSWER=1355 && COMMAND="./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:N:Nb:Zr --pocc_params=S0-0.5xA-0.5xC_S1-0.333xB-0.333xD-0.333xE | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
#5-metal perovskite oxides
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt Ba5HfNbSnTiZrO15" | tee -a $(COUNTER_FILE) && ANSWER=49 && COMMAND="./aflow --proto=AB3C_cP5_221_a_c_b.BAC:Ba:Hf:Nb:O:Sn:Ti:Zr --pocc_params=S0-1xD_S1-1xA_S2-0.2xB-0.2xC-0.2xE-0.2xF-0.2xG | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@COUNTER=$$(cat $(COUNTER_FILE) | wc -l | sed 's/^ *//') && echo "Running POCC test $$COUNTER - unique structures count of rocksalt Ba5HfSnTiYZrO15" | tee -a $(COUNTER_FILE) && ANSWER=49 && COMMAND="./aflow --proto=AB3C_cP5_221_a_c_b.BAC:Ba:Hf:O:Sn:Ti:Y:Zr --pocc_params=S0-1xC_S1-1xA_S2-0.2xB-0.2xD-0.2xE-0.2xF-0.2xG | ./aflow --pocc_count_unique" && COUNT=$$(eval "$$COMMAND" | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" != "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for '$$COMMAND', found $$COUNT instead" && exit 1; fi
	@echo "POCC tests passed"
	@rm -f $(COUNTER_FILE)
.PHONY: check_aflowSG
check_aflowSG:
	@if [ ! -f ./aflow ]; then echo "No local aflow binary found (./aflow)" && exit 1; fi
	@echo "Running aflowSG test 1" && ANSWER="Pm-3m #221" && RESULT=$$(./aflow --proto=AB_cP2_221_b_a:Cl:Cs | ./aflow --aflowSG) && \
		if [ "$$RESULT" != "$$ANSWER" ]; then "aflowSG did not return $$ANSWER for './aflow --proto=AB_cP2_221_b_a:Cl:Cs | ./aflow --aflowSG', found $$RESULT instead" && exit 1; fi
	@echo "aflowSG test 1 passed!"
.PHONY: check_AFLOWSYM
check_AFLOWSYM:
	@if [ ! -f ./aflow ]; then echo "No local aflow binary found (./aflow)" && exit 1; fi
	@echo "Running aflowSYM test 1" && ANSWER_PGROUPS=8 && ANSWER_PGROUPKS=8 && ANSWER_FGROUPS=4 && ANSWER_PGROUP_XTALS=4 && ANSWER_PGROUPK_XTALS=4 && ANSWER_AGROUPS=24 && ANSWER_SGROUPS=108 && RESULT=$$(./aflow --proto=ABC2_mP8_10_ac_eh_mn-001:Au:Ag:Te | ./aflow --aflowSYM --screen_only) && COUNT_PGROUPS=$$(echo $$RESULT | awk -F'pgroup ' '{print NF-1}') && COUNT_PGROUPKS=$$(echo $$RESULT | awk -F'pgroupk ' '{print NF-1}') && COUNT_FGROUPS=$$(echo $$RESULT | awk -F'fgroup ' '{print NF-1}') && COUNT_PGROUP_XTALS=$$(echo $$RESULT | awk -F'pgroup_xtal ' '{print NF-1}') && COUNT_PGROUPK_XTALS=$$(echo $$RESULT | awk -F'pgroupk_xtal ' '{print NF-1}') && COUNT_AGROUPS=$$(echo $$RESULT | awk -F'agroup ' '{print NF-1}') && COUNT_SGROUPS=$$(echo $$RESULT | awk -F'sgroup ' '{print NF-1}') && \
		if [ "$$ANSWER_PGROUPS" -ne "$$COUNT_PGROUPS" ] || [ "$$ANSWER_PGROUPKS" -ne "$$COUNT_PGROUPKS" ] || [ "$$ANSWER_FGROUPS" -ne "$$COUNT_FGROUPS" ] || [ "$$ANSWER_PGROUP_XTALS" -ne "$$COUNT_PGROUP_XTALS" ] || [ "$$ANSWER_PGROUPK_XTALS" -ne "$$COUNT_PGROUPK_XTALS" ] || [ "$$ANSWER_AGROUPS" -ne "$$COUNT_AGROUPS" ] || [ "$$ANSWER_SGROUPS" -ne "$$COUNT_SGROUPS" ]; then "aflowSYM did not return one of the correct group cardinalities for './aflow --proto=ABC2_mP8_10_ac_eh_mn-001:Au:Ag:Te | ./aflow --aflowSYM --screen_only'; (ANSWER==FOUND) $$ANSWER_PGROUPS==$$COUNT_PGROUPS, $$ANSWER_PGROUPKS==$$COUNT_PGROUPKS, $$ANSWER_FGROUPS==$$COUNT_FGROUPS, $$ANSWER_PGROUP_XTALS==$$COUNT_PGROUP_XTALS, $$ANSWER_PGROUPK_XTALS==$$COUNT_PGROUPK_XTALS, $$ANSWER_AGROUPS==$$COUNT_AGROUPS, $$ANSWER_SGROUPS==$$COUNT_SGROUPS"  && exit 1; fi
	@echo "aflowSYM test 1 passed!"
.PHONY: check_edata
check_edata:
	@if [ ! -f ./aflow ]; then echo "No local aflow binary found (./aflow)" && exit 1; fi
	@echo "Running edata test 1" && ANSWER=1 RESULT=$$(./aflow --proto=AB2_tI6_139_a_e:Os:Rh --params=2.7095,4.27466,0.662829 | ./aflow --edata | grep "BRAVAIS LATTICE OF THE LATTICE (pgroup)" | wc -l | sed 's/^ *//') && \
		if [ "$$ANSWER" -ne "$$RESULT" ]; then "edata did not finish successfully" && exit 1; fi
	@echo "edata test 1 passed!"
.PHONY: check_aflow_xtal_finder
check_aflow_xtal_finder:
	@if [ ! -f ./aflow ]; then echo "No local aflow binary found (./aflow)" && exit 1; fi
	@echo "Running aflow_xtal_finder - ideal prototyper test 1" && ANSWER="A2B3_hR10_167_c_e" && RESULT=$$(./aflow --proto=A2B3_hR10_167_c_e:Al:O --params=5.22337,2.72958,0.64784,0.9439 | ./aflow --prototype) && AFLOW_LABEL=$$(echo $$RESULT | grep "AFLOW label" | awk '{print $$4}') && \
		if [ $$AFLOW_LABEL != $$ANSWER ]; then "aflow_xtal_finder did not return $$ANSWER for './aflow --proto=A2B3_hR10_167_c_e:Al:O --params=5.22337,2.72958,0.64784,0.9439 | ./aflow --prototype', found $$AFLOW_LABEL instead" && exit 1; fi
	@echo "aflow_xtal_finder - ideal prototyper test 1 passed!"
	@echo "Running aflow_xtal_finder - ideal prototyper test 2" && ANSWER="AB18C8_cF108_225_a_eh_f" && RESULT=$$(./aflow --proto=AB18C8_cF108_225_a_eh_f --params=10.3439,0.325,0.65833,0.16 | ./aflow --prototype) && AFLOW_LABEL=$$(echo $$RESULT | grep "AFLOW label" | awk '{print $$4}') && \
		if [ $$AFLOW_LABEL != $$ANSWER ]; then "aflow_xtal_finder did not return $$ANSWER for './aflow --proto=AB18C8_cF108_225_a_eh_f --params=10.3439,0.325,0.65833,0.16 | ./aflow --prototype', found $$AFLOW_LABEL instead" && exit 1; fi
	@echo "aflow_xtal_finder - ideal prototyper test 2 passed!"
	@echo "Running aflow_xtal_finder - ideal prototyper test 3" && ANSWER="ABC4_mP12_13_e_a_2g" && RESULT=$$(./aflow --proto=ABC4_mP12_13_e_a_2g --params=8.6244,0.500336,1.63361,145.35,0.5182,0.2986,0.0278,0.0003,0.2821,0.4045,0.2366 | ./aflow --prototype) && AFLOW_LABEL=$$(echo $$RESULT | grep "AFLOW label" | awk '{print $$4}') && \
		if [ $$AFLOW_LABEL != $$ANSWER ]; then "aflow_xtal_finder did not return $$ANSWER for './aflow --proto=ABC4_mP12_13_e_a_2g --params=8.6244,0.500336,1.63361,145.35,0.5182,0.2986,0.0278,0.0003,0.2821,0.4045,0.2366 | ./aflow --prototype', found $$AFLOW_LABEL instead" && exit 1; fi
	@echo "aflow_xtal_finder - ideal prototyper test 3 passed!"
	@echo "Running aflow_xtal_finder - compare materials test 1 (different cell sizes)" && ANSWER="MATCH" && RESULT=$$(./aflow --proto=AB_cF8_225_a_b:Cl:Na > AB_cF8_225_a_b.poscar; cat AB_cF8_225_a_b.poscar | ./aflow --supercell=4,1,1 > AB_cF8_225_a_b_411_supercell.poscar; ./aflow --compare_materials=AB_cF8_225_a_b.poscar,AB_cF8_225_a_b_411_supercell.poscar --quiet) && MATCH=$$(echo $$RESULT | awk '{print $$3}') && \
		if [ $$MATCH != $$ANSWER ]; then "aflow_xtal_finder did not return $$ANSWER for './aflow --proto=AB_cF8_225_a_b:Cl:Na > AB_cF8_225_a_b; cat AB_cF8_225_a_b.poscar | ./aflow --supercell=4,1,1 > AB_cF8_225_a_b_411_supercell.poscar; ./aflow --compare_materials=AB_cF8_225_a_b.poscar,AB_cF8_225_a_b_411_supercell.poscar --quiet', found $$MATCH instead" && exit 1; fi
	@echo "aflow_xtal_finder - compare_materials test 1 (different cell sizes) passed!"
	@rm -f AB_cF8_225_a_b.poscar
	@rm -f AB_cF8_225_a_b_411_supercell.poscar
	@echo "Running aflow_xtal_finder - compare materials test 2 (different cell sizes)" && ANSWER="MATCH" && RESULT=$$(./aflow --proto=AB_cF8_225_a_b:Cl:Na > AB_cF8_225_a_b.poscar; cat AB_cF8_225_a_b.poscar | ./aflow --supercell=4,1,1 > AB_cF8_225_a_b_411_supercell.poscar; ./aflow --compare_materials=AB_cF8_225_a_b_411_supercell.poscar,AB_cF8_225_a_b.poscar --quiet) && MATCH=$$(echo $$RESULT | awk '{print $$3}') && \
		if [ $$MATCH != $$ANSWER ]; then "aflow_xtal_finder did not return $$ANSWER for './aflow --proto=AB_cF8_225_a_b:Cl:Na > AB_cF8_225_a_b; cat AB_cF8_225_a_b.poscar | ./aflow --supercell=4,1,1 > AB_cF8_225_a_b_411_supercell.poscar; ./aflow --compare_materials=AB_cF8_225_a_b_411_supercell.poscar,AB_cF8_225_a_b.poscar --quiet', found $$MATCH instead" && exit 1; fi
	@echo "aflow_xtal_finder - compare_materials test 2 (different cell sizes) passed!"
	@rm -f AB_cF8_225_a_b.poscar
	@rm -f AB_cF8_225_a_b_411_supercell.poscar
	@echo "Running aflow_xtal_finder - compare atom decorations test 1" && ANSWER=3 && RESULT=$$(./aflow --proto=ABCD_cF16_216_c_d_b_a | ./aflow --unique_atom_decorations --print_misfit --quiet) && NUM_UNIQUE_DECORATIONS=$$(echo $$RESULT | grep "Unique atom decorations" | cut -d '(' -f2 | cut -d ')' -f1) && \
		if [ "$$NUM_UNIQUE_DECORATIONS" -ne "$$ANSWER" ]; then "aflow_xtal_finder did not return $$ANSWER for './aflow --proto=ABCD_cF16_216_c_d_b_a | ./aflow --unique_atom_decorations --print_misfit --quiet', found $$NUM_UNIQUE_DECORATIONS instead" && exit 1; fi
	@echo "aflow_xtal_finder test - compare atom decorations test 1 passed!"
	@echo "Running aflow_xtal_finder - compare atom decorations test 2 " && ANSWER=3 && RESULT=$$(./aflow --proto=ABC_hP3_156_a_c_b --params=4.346,1.57271,4.44089e-16,0.4227,0.7214 | ./aflow --unique_atom_decorations --print_misfit --optimize --quiet) && NUM_UNIQUE_DECORATIONS=$$(echo $$RESULT | grep "Unique atom decorations" | cut -d '(' -f2 | cut -d ')' -f1) && \
		if [ "$$NUM_UNIQUE_DECORATIONS" -ne "$$ANSWER" ]; then "aflow_xtal_finder did not return $$ANSWER for './aflow --proto=ABC_hP3_156_a_c_b --params=4.346,1.57271,4.44089e-16,0.4227,0.7214 | ./aflow --unique_atom_decorations --print_misfit --optimize --quiet', found $$NUM_UNIQUE_DECORATIONS instead" && exit 1; fi
	@echo "aflow_xtal_finder - compare atom decorations test 1 passed!"
	@echo "Running aflow_xtal_finder - compare2protos test 1" && ANSWER=2 && RESULT=$$(./aflow --proto=A2B3_hR10_167_c_e:Al:O --params=-1,2.72957758313,0.36,0.5561 | ./aflow --compare2protos --print=txt --catalog=anrl --quiet) && NUM_MATCHING_PROTOS=$$(echo $$RESULT | awk -F'duplicate_compounds=' '{print $$2}' | awk '{print $$1}') && \
		if [ "$$NUM_MATCHING_PROTOS" -ne "$$ANSWER" ]; then "aflow_xtal_finder did not return $$ANSWER for './aflow --proto=A2B3_hR10_167_c_e:Al:O --params=-1,2.72957758313,0.36,0.5561 | ./aflow --compare2protos --print=txt --catalog=anrl --quiet', found $$NUM_MATCHING_PROTOS instead" && exit 1; fi
	@echo "aflow_xtal_finder - compare2protos test 1 passed!"
	@echo "Running aflow_xtal_finder - compare magnetic configurations test 1 (antiferro vs antiferro)" && ANSWER="MATCH" && RESULT=$$(./aflow --proto=A_cI2_229_a:Cr | ./aflow --sconv >A_cI2_229_a_sconv.poscar ; ./aflow --compare_materials=A_cI2_229_a_sconv.poscar,A_cI2_229_a_sconv.poscar --magmoms=-1,1:1,-1 --quiet) && MATCH=$$(echo $$RESULT | awk '{print $$3}') && \
		if [ $$MATCH != $$ANSWER ]; then "aflow_xtal_finder did not return $$ANSWER for './aflow --proto=A_cI2_229_a:Cr | ./aflow --sconv >A_cI2_229_a_sconv.poscar ; ./aflow --compare_materials=A_cI2_229_a_sconv.poscar,A_cI2_229_a_sconv.poscar --magmoms=-1,1:1,-1 --quiet', found $$MATCH instead" && exit 1; fi
	@echo "aflow_xtal_finder - compare magnetic configuration test 1 (antiferro vs antiferro) passed!"
	@rm -f A_cI2_229_a_sconv.poscar
	@echo "Running aflow_xtal_finder - compare magnetic configurations test 2 (ferro vs antiferro)" && ANSWER="UNMATCHABLE" && RESULT=$$(./aflow --proto=A_cI2_229_a:Cr | ./aflow --sconv >A_cI2_229_a_sconv.poscar ; ./aflow --compare_materials=A_cI2_229_a_sconv.poscar,A_cI2_229_a_sconv.poscar --magmoms=1,1:1,-1 --quiet) && MATCH=$$(echo $$RESULT | awk '{print $$1}') && \
		if [ $$MATCH != $$ANSWER ]; then "aflow_xtal_finder did not return $$ANSWER for './aflow --proto=A_cI2_229_a:Cr | ./aflow --sconv >A_cI2_229_a_sconv.poscar ; ./aflow --compare_materials=A_cI2_229_a_sconv.poscar,A_cI2_229_a_sconv.poscar --magmoms=1,1:1,-1 --quiet', found $$MATCH instead" && exit 1; fi
	@echo "aflow_xtal_finder - compare magnetic configuration test 2 (ferro vs antiferro) passed!"
	@rm -f A_cI2_229_a_sconv.poscar
	@echo "Running aflow_xtal_finder - compare magnetic configurations test 3 (ferro-up vs ferro-down)" && ANSWER="UNMATCHABLE" && RESULT=$$(./aflow --proto=A_cI2_229_a:Cr | ./aflow --sconv >A_cI2_229_a_sconv.poscar ; ./aflow --compare_materials=A_cI2_229_a_sconv.poscar,A_cI2_229_a_sconv.poscar --magmoms=1,1:-1,-1 --quiet) && MATCH=$$(echo $$RESULT | awk '{print $$1}') && \
		if [ $$MATCH != $$ANSWER ]; then "aflow_xtal_finder did not return $$ANSWER for './aflow --proto=A_cI2_229_a:Cr | ./aflow --sconv >A_cI2_229_a_sconv.poscar ; ./aflow --compare_materials=A_cI2_229_a_sconv.poscar,A_cI2_229_a_sconv.poscar --magmoms=1,1:-1,-1 --quiet', found $$MATCH instead" && exit 1; fi
	@echo "aflow_xtal_finder - compare magnetic configuration test 3 (ferro-up vs ferro-down) passed!"
	@rm -f A_cI2_229_a_sconv.poscar
#RF20200302 - START
.PHONY: check_cce
check_cce:
	@if [ ! -f ./aflow ]; then echo "No local aflow binary found (./aflow)" && exit 1; fi
	@echo "Performing 'aflow --cce' test"
	@echo "Running CCE test 1: MgO"
	@./aflow --proto=AB_cF8_225_a_b:Mg:O --params=4.253 > POSCAR_proto_Mg1O1_ICSD_159378; ANSWER=0.801,0.763,0.015,-0.025,-0.014,-0.054,-6.235,-6.197,-6.235,-6.195,-6.235,-6.195; RESULT=$$(./aflow --cce=POSCAR_proto_Mg1O1_ICSD_159378 --cce_test --dft_formation_energies=-5.434,-6.220,-6.249 --functionals=PBE,LDA,SCAN --oxidation_numbers=2,-2); if [ "$$RESULT" != "$$ANSWER" ]; then echo "CCE test 1 failed!, check CCE on POSCAR_proto_Mg1O1_ICSD_159378" && exit 1; else echo "CCE test 1 passed!"; rm -f POSCAR_proto_Mg1O1_ICSD_159378 ; fi
	@echo "Running CCE test 2: Ti3O5"
	@wget --quiet https://aflowlib.duke.edu/AFLOWDATA/ICSD_WEB/MCLC/O5Ti3_ICSD_75193/CONTCAR.relax.vasp; ANSWER=-3.210,-3.532,3.638,3.338,-3.134,-3.454,-50.754,-50.432,-51.234,-50.934,-50.870,-50.550; mv CONTCAR.relax.vasp CONTCAR.relax.vasp_O5Ti3_ICSD_75193; RESULT=$$(./aflow --cce=CONTCAR.relax.vasp_O5Ti3_ICSD_75193 --cce_test --dft_formation_energies=-53.964,-47.596,-54.004 --functionals=LDA,PBE,SCAN); if [ "$$RESULT" != "$$ANSWER" ]; then echo "CCE test 2 failed!, check CCE on CONTCAR.relax.vasp_O5Ti3_ICSD_75193" && exit 1; else echo "CCE test 2 passed!"; rm -f CONTCAR.relax.vasp_O5Ti3_ICSD_75193 ; fi
	@echo "Running CCE test 3: Mg2SiO4"
	@./aflow --proto=A2B4C_oP28_62_ac_2cd_c:Mg:O:Si --params=10.193,0.586382811734,0.466202295693,0.2774,-0.0085,0.0913,0.7657,0.4474,0.2215,0.094,0.4262,0.1628,0.0331,0.2777 > POSCAR_proto_Mg2O4Si1_ICSD_83793; ANSWER=10.456,9.912,-0.252,-0.824,-0.444,-1.020,-89.916,-89.372,-89.860,-89.288,-90.000,-89.424; RESULT=$$(./aflow --cce=POSCAR_proto_Mg2O4Si1_ICSD_83793 --cce_test --dft_formation_energies=-79.46,-90.112,-90.444 --functionals=PBE,LDA,SCAN --oxidation_numbers=2,2,2,2,2,2,2,2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,4,4,4,4); if [ "$$RESULT" != "$$ANSWER" ]; then echo "CCE test 3 failed!, check CCE on POSCAR_proto_Mg2O4Si1_ICSD_83793" && exit 1; else echo "CCE test 3 passed!"; rm -f POSCAR_proto_Mg2O4Si1_ICSD_83793 ; fi
	@echo "Running CCE test 4: BaTiO3"
	@./aflow --proto=AB3C_cP5_221_a_c_b:Ba:O:Ti --params=4.033 > POSCAR_proto_Ba1O3Ti1_ICSD_187292; ANSWER=2.063,1.981,-0.461,-0.553,-0.692,-0.783,-17.513,-17.431,-17.180,-17.088,-17.010,-16.919; RESULT=$$(./aflow --cce=POSCAR_proto_Ba1O3Ti1_ICSD_187292 --cce_test --dft_formation_energies=-15.450,-17.641,-17.702 --functionals=PBE,LDA,SCAN); if [ "$$RESULT" != "$$ANSWER" ]; then echo "CCE test 4 failed!, check CCE on POSCAR_proto_Ba1O3Ti1_ICSD_187292" && exit 1; else echo "CCE test 4 passed!"; rm -f POSCAR_proto_Ba1O3Ti1_ICSD_187292 ; fi
	@echo "Running CCE test 5: Na2WO4"
	@./aflow --proto=A2BC4_cF56_227_d_a_e:Na:W:O --params=9.0433,0.260612 > POSCAR_proto_Na2O4W1_ICSD_44524; ANSWER=2.014,1.843,-1.765,-1.963,-1.359,-1.561,-32.470,-32.299,-31.827,-31.629,-32.513,-32.311; RESULT=$$(./aflow --cce=POSCAR_proto_Na2O4W1_ICSD_44524 --cce_test --dft_formation_energies=-30.456,-33.592,-33.872 --functionals=PBE,LDA,SCAN --oxidation_numbers=1,1,1,1,-2,-2,-2,-2,-2,-2,-2,-2,6,6); if [ "$$RESULT" != "$$ANSWER" ]; then echo "CCE test 5 failed!, check CCE on POSCAR_proto_Na2O4W1_ICSD_44524" && exit 1; else echo "CCE test 5 passed!"; rm -f POSCAR_proto_Na2O4W1_ICSD_44524 ; fi
	@echo "Running CCE test 6: PbMoO4"
	@wget --quiet https://aflowlib.duke.edu/AFLOWDATA/ICSD_WEB/BCT/Mo1O4Pb1_ICSD_89034/CONTCAR.relax.vasp; ANSWER=-2.922,-3.054,-0.332,-0.430,-4.468,-4.608,-21.832,-21.700,-22.330,-22.232,-21.104,-20.964; mv CONTCAR.relax.vasp CONTCAR.relax.vasp_Mo1O4Pb1_ICSD_89034; RESULT=$$(./aflow --cce=CONTCAR.relax.vasp_Mo1O4Pb1_ICSD_89034 --cce_test --dft_formation_energies=-24.754,-22.662,-25.572 --functionals=SCAN,PBE,LDA); if [ "$$RESULT" != "$$ANSWER" ]; then echo "CCE test 6 failed!, check CCE on CONTCAR.relax.vasp_Mo1O4Pb1_ICSD_89034" && exit 1; else echo "CCE test 6 passed!"; rm -f CONTCAR.relax.vasp_Mo1O4Pb1_ICSD_89034 ; fi
	@echo "Running CCE test 7: Fe3O4"
	@./aflow --proto=A3B4_cF56_227_ad_e:Fe:O --params=8.3886,0.745448 > POSCAR_proto_Fe3O4_ICSD_263007; ANSWER=5.381,5.246,1.832,1.654,-1.084,-1.263,-23.575,-23.440,-23.068,-22.890,-24.320,-24.141; RESULT=$$(./aflow --cce=POSCAR_proto_Fe3O4_ICSD_263007 --cce_test --dft_formation_energies=-18.194,-21.236,-25.404 --functionals=PBE,LDA,SCAN); if [ "$$RESULT" != "$$ANSWER" ]; then echo "CCE test 7 failed!, check CCE on POSCAR_proto_Fe3O4_ICSD_263007" && exit 1; else echo "CCE test 7 passed!"; rm -f POSCAR_proto_Fe3O4_ICSD_263007 ; fi
	@echo "Running CCE test 8: Ca2V2O7"
	@wget --quiet https://aflowlib.duke.edu/AFLOWDATA/ICSD_WEB/TRI/Ca2O7V2_ICSD_421266/CONTCAR.relax.vasp; ANSWER=2.600,2.264,-4.675,-5.073,-4.499,-4.881,-64.422,-64.086,-64.113,-63.715,-65.011,-64.629; mv CONTCAR.relax.vasp CONTCAR.relax.vasp_Ca2O7V2_ICSD_421266; RESULT=$$(./aflow --cce=CONTCAR.relax.vasp_Ca2O7V2_ICSD_421266 --cce_test --dft_formation_energies=-61.822,-68.788,-69.51 --functionals=PBE,LDA,SCAN); if [ "$$RESULT" != "$$ANSWER" ]; then echo "CCE test 8 failed!, check CCE on CONTCAR.relax.vasp_Ca2O7V2_ICSD_421266" && exit 1; else echo "CCE test 8 passed!"; rm -f CONTCAR.relax.vasp_Ca2O7V2_ICSD_421266 ; fi
	@echo "Running CCE test 9: CuFeO2"
	@wget --quiet https://aflowlib.duke.edu/AFLOWDATA/ICSD_WEB/RHL/Cu1Fe1O2_ICSD_246912/CONTCAR.relax.vasp; ANSWER=1.242,1.210,0.146,0.097,-0.266,-0.307,-5.152; mv CONTCAR.relax.vasp CONTCAR.relax.vasp_Cu1Fe1O2_ICSD_246912; RESULT=$$(./aflow --cce=CONTCAR.relax.vasp_Cu1Fe1O2_ICSD_246912 --cce_test); if [ "$$RESULT" != "$$ANSWER" ]; then echo "CCE test 9 failed!, check CCE on CONTCAR.relax.vasp_Cu1Fe1O2_ICSD_246912" && exit 1; else echo "CCE test 9 passed!"; rm -f CONTCAR.relax.vasp_Cu1Fe1O2_ICSD_246912 ; fi
	@echo "Running CCE test 10: MnMoO4"
	@wget --quiet https://aflowlib.duke.edu/AFLOWDATA/ICSD_WEB/MCL/Mn1Mo1O4_ICSD_61078/CONTCAR.relax.vasp; ANSWER=2.452,2.293,-0.852,-1.058,-3.459,-3.669,-24.390,-24.231,-25.458,-25.252,-24.891,-24.681; mv CONTCAR.relax.vasp CONTCAR.relax.vasp_Mn1Mo1O4_ICSD_61078; RESULT=$$(./aflow --cce=CONTCAR.relax.vasp_Mn1Mo1O4_ICSD_61078 --cce_test --dft_formation_energies=-21.938,-26.31,-28.35 --functionals=PBE,LDA,SCAN); if [ "$$RESULT" != "$$ANSWER" ]; then echo "CCE test 10 failed!, check CCE on CONTCAR.relax.vasp_Mn1Mo1O4_ICSD_61078" && exit 1; else echo "CCE test 10 passed!"; rm -f CONTCAR.relax.vasp_Mn1Mo1O4_ICSD_61078 ; fi
	@echo "Running CCE test 11: KO2"
	@wget --quiet https://aflowlib.duke.edu/AFLOWDATA/ICSD_WEB/BCT/K1O2_ICSD_38245/CONTCAR.relax.vasp; ANSWER=0.285,0.259,-0.532,-0.571,-0.084,-0.118,-2.949,-2.923,-2.949,-2.910,-2.949,-2.915; mv CONTCAR.relax.vasp CONTCAR.relax.vasp_K1O2_ICSD_38245; RESULT=$$(./aflow --cce=CONTCAR.relax.vasp_K1O2_ICSD_38245 --cce_test --dft_formation_energies=-2.664,-3.481,-3.033 --functionals=PBE,LDA,SCAN); if [ "$$RESULT" != "$$ANSWER" ]; then echo "CCE test 11 failed!, check CCE on CONTCAR.relax.vasp_K1O2_ICSD_38245" && exit 1; else echo "CCE test 11 passed!"; rm -f CONTCAR.relax.vasp_K1O2_ICSD_38245 ; fi
	@echo "Running CCE test 12: Li2O"
	@./aflow --proto=AB2_cF12_225_a_c:O:Li --params=4.5709 > POSCAR_proto_Li2O1_ICSD_642216; ANSWER=0.613,0.563,-0.123,-0.178,-0.094,-0.149,-6.197,-6.147,-6.197,-6.142,-6.197,-6.142; RESULT=$$(./aflow --cce=POSCAR_proto_Li2O1_ICSD_642216 --cce_test --dft_formation_energies=-5.584,-6.320,-6.291 --functionals=PBE,LDA,SCAN); if [ "$$RESULT" != "$$ANSWER" ]; then echo "CCE test 12 failed!, check CCE on POSCAR_proto_Li2O1_ICSD_642216" && exit 1; else echo "CCE test 12 passed!"; rm -f POSCAR_proto_Li2O1_ICSD_642216 ; fi
	@echo "Running CCE test 13: Li2O2"
	@wget --quiet https://aflowlib.duke.edu/AFLOWDATA/ICSD_WEB/HEX/Li2O2_ICSD_152183/CONTCAR.relax.vasp; ANSWER=1.653,1.518,-0.601,-0.756,0.198,0.048,-13.139,-13.004,-13.141,-12.986,-13.140,-12.990; mv CONTCAR.relax.vasp CONTCAR.relax.vasp_Li2O2_ICSD_152183; RESULT=$$(./aflow --cce=CONTCAR.relax.vasp_Li2O2_ICSD_152183 --cce_test --dft_formation_energies=-11.486,-13.742,-12.942 --functionals=PBE,LDA,SCAN); if [ "$$RESULT" != "$$ANSWER" ]; then echo "CCE test 13 failed!, check CCE on CONTCAR.relax.vasp_Li2O2_ICSD_152183" && exit 1; else echo "CCE test 13 passed!"; rm -f CONTCAR.relax.vasp_Li2O2_ICSD_152183 ; fi
	@echo "Running CCE test 14: CaFe2O4"
	@wget --quiet https://aflowlib.duke.edu/AFLOWDATA/ICSD_WEB/ORC/Ca1Fe2O4_ICSD_159751/CONTCAR.relax.vasp; ANSWER=11.217,10.817,-0.437,-1.000,-3.853,-4.340,-69.230; mv CONTCAR.relax.vasp CONTCAR.relax.vasp_Ca1Fe2O4_ICSD_159751; RESULT=$$(./aflow --cce=CONTCAR.relax.vasp_Ca1Fe2O4_ICSD_159751 --cce_test); if [ "$$RESULT" != "$$ANSWER" ]; then echo "CCE test 14 failed!, check CCE on CONTCAR.relax.vasp_Ca1Fe2O4_ICSD_159751" && exit 1; else echo "CCE test 14 passed!"; rm -f CONTCAR.relax.vasp_Ca1Fe2O4_ICSD_159751 ; fi
	@echo "Running CCE test 15: Ca2Fe2O5"
	@wget --quiet https://aflowlib.duke.edu/AFLOWDATA/ICSD_WEB/ORC/Ca2Fe2O5_ICSD_161509/CONTCAR.relax.vasp; ANSWER=11.602,11.151,-1.072,-1.676,-3.684,-4.214,-81.087; mv CONTCAR.relax.vasp CONTCAR.relax.vasp_Ca2Fe2O5_ICSD_161509; RESULT=$$(./aflow --cce=CONTCAR.relax.vasp_Ca2Fe2O5_ICSD_161509 --cce_test); if [ "$$RESULT" != "$$ANSWER" ]; then echo "CCE test 15 failed!, check CCE on CONTCAR.relax.vasp_Ca2Fe2O5_ICSD_161509" && exit 1; else echo "CCE test 15 passed!"; rm -f CONTCAR.relax.vasp_Ca2Fe2O5_ICSD_161509 ; fi
	@echo "Running CCE test 16: FeTiO3"
	@./aflow --proto=AB3C_hR10_148_c_f_c:Fe:O:Ti --params=5.08264,2.77119,0.648852,0.852345,0.436666,0.768665,0.0556432 > POSCAR_proto_Fe1O3Ti1_ICSD_187688; ANSWER=3.402,3.240,0.414,0.258,-1.258,-1.442,-25.786,-25.624,-25.612,-25.456,-25.808,-25.624; RESULT=$$(./aflow --cce=POSCAR_proto_Fe1O3Ti1_ICSD_187688 --cce_test --dft_formation_energies=-22.384,-25.198,-27.066 --functionals=PBE,LDA,SCAN); if [ "$$RESULT" != "$$ANSWER" ]; then echo "CCE test 16 failed!, check CCE on POSCAR_proto_Fe1O3Ti1_ICSD_187688" && exit 1; else echo "CCE test 16 passed!"; rm -f POSCAR_proto_Fe1O3Ti1_ICSD_187688 ; fi
	@echo "Running CCE test 17: BaO2"
	@wget --quiet https://aflowlib.duke.edu/AFLOWDATA/ICSD_WEB/BCT/Ba1O2_ICSD_24248/CONTCAR.relax.vasp; ANSWER=1.055,1.043,-6.638,-6.626,-6.743; mv CONTCAR.relax.vasp CONTCAR.relax.vasp_Ba1O2_ICSD_24248; RESULT=$$(./aflow --cce=CONTCAR.relax.vasp_Ba1O2_ICSD_24248 --cce_test --dft_formation_energies=-5.583,0 --functionals=PBE+U:ICSD,exp); if [ "$$RESULT" != "$$ANSWER" ]; then echo "CCE test 17 failed!, check CCE on CONTCAR.relax.vasp_Ba1O2_ICSD_24248" && exit 1; else echo "CCE test 17 passed!"; rm -f CONTCAR.relax.vasp_Ba1O2_ICSD_24248 ; fi
	@echo "Running CCE test 18: NaO2"
	@wget --quiet https://aflowlib.duke.edu/AFLOWDATA/ICSD_WEB/CUB/Na1O2_ICSD_87178/CONTCAR.relax.vasp; ANSWER=-0.205,-0.310,-9.363,-9.258,-5.972; mv CONTCAR.relax.vasp CONTCAR.relax.vasp_Na1O2_ICSD_87178; RESULT=$$(./aflow --cce=CONTCAR.relax.vasp_Na1O2_ICSD_87178 --cce_test --dft_formation_energies=-9.568,0 --functionals=PBE+U:ICSD,exp); if [ "$$RESULT" != "$$ANSWER" ]; then echo "CCE test 18 failed!, check CCE on CONTCAR.relax.vasp_Na1O2_ICSD_87178" && exit 1; else echo "CCE test 18 passed!"; rm -f CONTCAR.relax.vasp_Na1O2_ICSD_87178 ; fi
	@echo "Running CCE test 19: HfN"
	@./aflow --proto=AB_cF8_225_a_b:Hf:N --params=4.4938 > POSCAR_proto_Hf1N1_ICSD_53025; ANSWER=0.359,0.338,-0.557,-0.581,-0.086,-0.109,0.356,0.335,-3.872,-3.851,-3.872,-3.848,-3.872,-3.849,-3.872,-3.851,-3.872; RESULT=$$(./aflow --cce=POSCAR_proto_Hf1N1_ICSD_53025 --cce_test --dft_formation_energies=-3.513,-4.429,-3.958,-3.516,0 --functionals=PBE,LDA,SCAN,PBE+U:ICSD,exp); if [ "$$RESULT" != "$$ANSWER" ]; then echo "CCE test 19 failed!, check CCE on POSCAR_proto_Hf1N1_ICSD_53025" && exit 1; else echo "CCE test 19 passed!"; rm -f POSCAR_proto_Hf1N1_ICSD_53025 ; fi
	@echo "Running CCE test 20: LiCaN"
	@wget --quiet https://aflowlib.duke.edu/AFLOWDATA/ICSD_WEB/ORC/Ca1Li1N1_ICSD_107304/CONTCAR.relax.vasp; ANSWER=0.834,0.653,-2.298,-2.474,-1.133,-1.311,0.960,0.792,-8.794,-8.613,-8.814,-8.638,-8.787,-8.609,-8.840,-8.672,-8.518; mv CONTCAR.relax.vasp CONTCAR.relax.vasp_Ca1Li1N1_ICSD_107304; RESULT=$$(./aflow --cce=CONTCAR.relax.vasp_Ca1Li1N1_ICSD_107304 --cce_test --dft_formation_energies=-7.960,-11.112,-9.920,-7.880,0 --functionals=PBE,LDA,SCAN,PBE+U:ICSD,exp); if [ "$$RESULT" != "$$ANSWER" ]; then echo "CCE test 20 failed!, check CCE on CONTCAR.relax.vasp_Ca1Li1N1_ICSD_107304" && exit 1; else echo "CCE test 20 passed!"; rm -f CONTCAR.relax.vasp_Ca1Li1N1_ICSD_107304 ; fi
#RF20200302 - STOP

#[ME20220205 - Obsolete with C++11]
#ifeq ($(UNAME),Linux)
#	@echo "Running aflow_gcc_compatibility test 1: no 'nullptr'" && ANSWER=0 && NUMBER_NULLPTRS=$$(find . -name '*.h' -o -name '*.cpp' | xargs grep -P '^(?=[\s]*+[^//])[^//]*(nullptr)' | wc -l | sed 's/^ *//') && \
#		if [ "$$NUMBER_NULLPTRS" -ne "$$ANSWER" ]; then echo "Found $$NUMBER_NULLPTRS nullptrs in the source code; this is not backwards compatible with older GCC version. Please swap with 'NULL' and recompile." && exit 1; fi
#	@echo "aflow_gcc_compatibility test 1 passed!"
#else
#	@echo "aflow_gcc_compatibility test 1 (no 'nullptr') requires GNU grep, try: find . -name '*.h' -o -name '*.cpp' | xargs grep -P '^(?=[\s]*+[^//])[^//]*(nullptr)'"
#endif
check_aflow_exits:
	@echo "Running aflow_exits test 1: NO 'exit()'" && ANSWER=0 && NUMBER_EXITS=$$(find . -name '*.h' -o -name '*.cpp' | xargs grep 'exit(' | wc -l | sed 's/^ *//') && \
		if [ "$$NUMBER_EXITS" -ne "$$ANSWER" ]; then echo "Found $$NUMBER_EXITS exits in the source code. Please swap with 'aurostd::xerror()' and recompile." && exit 1; fi
unit_test:
	@if [ ! -f ./aflow ]; then echo "No local aflow binary found (./aflow)" && exit 1; fi
	./aflow --unit_test
makefile:
	./aflow --makefile && $(MAKE) -j 4
###################################################################################################
