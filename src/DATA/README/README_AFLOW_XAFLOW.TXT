****************************************************************************************************
*                                                                                                  *
*                        aflow - STEFANO CURTAROLO Duke University 2003-2021                       *
*                            High-Throughput ab-initio Computing Project                           *
*                                                                                                  *
****************************************************************************************************
LATEST VERSION OF THE EXECUTABLE:                 materials.duke.edu/AFLOW/xaflow
LATEST VERSION OF THE README:                     materials.duke.edu/AFLOW/README_XAFLOW.TXT
****************************************************************************************************

XAFLOW README:

README written by: Corey Oses (corey.oses@duke.edu)

****************************************************************************************************
INSTRUCTIONS:

xaflow downloads and installs the newest version of AFLOW automatically.
xaflow needs to be installed first before installing aflow.
N.B. xaflow is a makefile.

Installing xaflow (perform only once):
1. Download xaflow:
      wget materials.duke.edu/AFLOW/xaflow
  - Note that this may require you to first install wget on your machine, 
    if it is not all already available.
  - On mac: brew install wget.
  - On linux: sudo apt-get install wget.
2. Move/copy xaflow somewhere in your path:
      mkdir ~/bin
      export PATH=~/bin:$PATH
      mv xaflow ~/bin/xaflow
  - Add the "export PATH=~/bin:$PATH" line to your ~/.bashrc so 
    it always stays in your path.
3. Ensure the file can be executed: 
      chmod 555 ~/bin/xaflow
4. Check that xaflow is in your path: 
      which xaflow
  - A returned path means xaflow has been successfully installed.

Installing aflow (perform with every new version):
1. Run xaflow:
      xaflow
  - xaflow downloads, unpacks, and makes AFLOW in ~/src/AFLOW.
  - The compilation can be redirected to another path:
      xaflow AWD=/home/src/AFLOW
  - If $(AWD) already exists, it will first move $(AWD) to $(AWD)_OLD.
2. Install the newly compiled executable to your path:
      xaflow install
  - By default, the executable is installed to ~/bin.
  - The installation can be redirected to another path:
      xaflow install ULB=/usr/local/bin
  - As with all makefiles, the AWD and ULB variables can be set simultaneously:
      xaflow install AWD=/home/src/AFLOW ULB=/usr/local/bin
3. Check that aflow is in your path:
      which aflow
  - A returned path means aflow has been successfully installed.
4. Run aflow:
      aflow

****************************************************************************************************

For additional information contact: Stefano Curtarolo (stefano@duke.edu) or 
                                    Corey Oses (corey.oses@duke.edu)

****************************************************************************************************
*                                                                                                  *
*                        aflow - STEFANO CURTAROLO Duke University 2003-2021                       *
*                            High-Throughput ab-initio Computing Project                           *
*                                                                                                  *
****************************************************************************************************
