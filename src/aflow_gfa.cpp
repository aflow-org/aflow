// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *           Aflow ERIC PERIM MARTINS - Duke University 2014-2016          *
// *           Aflow DENISE FORD - Duke University 2016-2019                 *
// *                                                                         *  
// ***************************************************************************
// Original code written by Eric Perim Martins - mostly routines for calculating the
// atomic environments remain.

// Revised and expanded by Denise Ford - includes new DetermineDistinctAEs function,
// partially rewritten GetAtomicEnvironment function, and completely rewritten
// ComputeGFA and CalculateGFA functions.

#include "aflow.h"
#include "aflow_gfa.h"
#include "aflow_chull.h"
#include "aflow_compare_structure.h"

#define _DEBUG_GFA_ false //DX20201005

/********DEFINITIONS**********/

//atomic environments
//[CO20200731 - OBSOLETE]#define _DOUBLE_TOL_ 0.001 //Tolerance for double comparison
#define _MIN_NUM_ATM_ 500.0 //Number of atoms used to define the supercell
#define _MIN_SAMPLING_ 40 //Minimum number of neighbors for acceptable sampling
#define _CON_TOL_ 1.3 //Tolerance when calculating connectivity (multiplying constant < sqrt(2))

//gfa
// [OBSOLETE] #define _ZTOL_ 0.05 //Tolerance for positive formation enthalpies - was originally 0.05
//[CO20190628]#define KbT 0.025 // setting for room temperature - 0.025, could also be human body temperature - 0.0267
//DX20210111 [OBSOLETE - moved to xscalar] #define TEMPERATURE 300.0 //Kelvin //CO20190628 //DX20201005 - needs to be a float/double (not integer)
//DX20210111 [OBSOLETE - moved to xscalar] #define KbT (KBOLTZEV*TEMPERATURE) // setting for room temperature - 0.025, could also be human body temperature - 0.0267 //CO20190628 //DX20201006 - need parentheses around macro multiplication

/********DEFINITIONS**********/



//************** Function for obtaining all distances (beginning)
xmatrix<double> GetAllDistances (bool& samplingOK, xstructure str, xstructure superStr, int cellNumber) {
  //THIS IS FOR CALCULATING THE MAXIMUM GAPS. PART OF THE ATOMIC ENVIRONMENTS COMPUTATION.
  xvector<double> a, b, c, d; //Lattice vectors and distance vector
  double radius;
  //////////////////
  // For k (odd) unit cells in each direction, the central cell equivalent of atom "i" in the unit cell is going to be
  // atom "p*k²+p*k+p+1+(i-1)*k³" = "p(k²+k+1)+1+(i-1)*k³" with p = (int) k/2
  //////////////////
  uint k = (uint) cellNumber; 
  uint p = (uint) cellNumber/2;
  xmatrix<double> distances(1,1,str.atoms.size(),superStr.atoms.size());
  samplingOK=true; 

  //Defines radius to be considered
  for(int i=1;i<=3;i++){            // This radius is important because supercells are not 
    a[i]=str.lattice[1][i];         // spherical, so, without this, gaps around the edge of the 
    b[i]=str.lattice[2][i];         // supercell can be greatly exaggerated
    c[i]=str.lattice[3][i];
  }

  xvector<double> moduli;

  moduli[1]=modulus(a);
  moduli[2]=modulus(b);
  moduli[3]=modulus(c);
  sort(moduli);
  radius = p*moduli[1];            // This definition is eliminating ~50% of the supercell atoms (for small unit cells) 
  int out;

  for (uint i=1;i<=str.atoms.size();i++){
    out=0;
    for (uint j=1;j<=superStr.atoms.size();j++){
      d = (superStr.atoms[(p*(k*k+k+1)+(i-1)*k*k*k)].cpos - superStr.atoms[j-1].cpos); //atoms run from 0 unlike xmatrix and xvector

      if(modulus(d)<=radius)
        distances(i,j)=modulus(d);
      else {
        distances(i,j)=-1.0;  //Outside radius
        out++;}

    }

    if((superStr.atoms.size()-out)<_MIN_SAMPLING_)
      samplingOK = false;
    // cout << "WARNING!!! Sampling may be insufficient!" << endl;

  }
  return distances;
}
//************** Function for obtaining all distances (end)


//************** Function for obtaining the gap values and positions (beginning)
void FindGap(xvector<double>& maxGaps, xvector<double>& gapPositions, xmatrix<double> distances, 
    xstructure str, xstructure superStr, double new_tol){
  /*
     COMPUTES THE MAXIMUM GAP FROM THE DISTANCES. PART OF THE ATOMIC ENVIRONMENT ROUTINE.
     */
  xvector<double> individualDistances(superStr.atoms.size(),1);
  set(gapPositions, 9999.99);

  for(uint i=1;i<=str.atoms.size();i++){ //Must be done for each individual atom

    for(uint j=1;j<=superStr.atoms.size();j++)
      individualDistances(j)=distances(i,j); //Temporary

    sort(individualDistances);

    for(uint j=2;j<=superStr.atoms.size();j++){ 

      if((maxGaps(i)-new_tol) <= (individualDistances(j)-individualDistances(j-1))){
        if((individualDistances(j-1)<gapPositions(i) ||
              maxGaps(i)<(individualDistances(j)-individualDistances(j-1)-new_tol)) && 
            individualDistances(j-1)>_ZERO_TOL_LOOSE_) { //Excludes first gap and negative distances (atoms outside cutoff)
          gapPositions(i)=individualDistances(j-1);
          maxGaps(i)=individualDistances(j)-individualDistances(j-1);
        } 
      }

    }

    if(maxGaps(i)<_ZERO_TOL_LOOSE_) //Just a check
      cout << "WARNING!!! Failed to find gap!" << endl << endl;
  } 
}
//************** Function for obtaining the gap values and positions (end)

//************** Function for obtaining the gap values and positions in case of maximum gap rule failure (beginning)
void FindNewGap(uint i, xvector<double>& maxGaps, xvector<double>& gapPositions, xmatrix<double> distances, 
    xstructure superStr, double new_tol){

  double cutoff = gapPositions(i);
  xvector<double> individualDistances(superStr.atoms.size(),1);

  for(uint j=1;j<=superStr.atoms.size();j++)
    individualDistances(j)=distances(i,j); //Temporary

  sort(individualDistances);

  maxGaps(i)=0;
  gapPositions(i)=9999.666;

  for(uint j=2;j<=superStr.atoms.size();j++){ 

    if(individualDistances(j-1)<(cutoff-new_tol)) //Enforces new AE volume < old AE volume
      if((maxGaps(i)-new_tol) <= (individualDistances(j)-individualDistances(j-1))){
        if((individualDistances(j-1)<gapPositions(i) || 
              maxGaps(i)<(individualDistances(j)-individualDistances(j-1)-new_tol)) && 
            individualDistances(j-1)>_ZERO_TOL_LOOSE_) { //Excludes first gap and negative distances (atoms outside cutoff)
          gapPositions(i)=individualDistances(j-1);
          maxGaps(i)=individualDistances(j)-individualDistances(j-1);
        } 
      }

  }

  if(maxGaps(i)<_ZERO_TOL_LOOSE_){ //Gap rule failure
    for(uint j=2;j<=superStr.atoms.size();j++){ 
      if(individualDistances(j-1)<cutoff){
        maxGaps(i)=individualDistances(j)-individualDistances(j-1);
        gapPositions(i)=individualDistances(j-1);
      }
    }
  }

}
//************** Function for obtaining the gap values and positions in case of maximum gap rule failure (end)


//************** Function for building the AEs (beginning)
void BuildAE(deque<deque<int> > &environment, xvector<double> gapPositions, xstructure str, xstructure superStr,
    int cellNumber, double new_tol){
  xvector<double> d;
  uint k = (uint) cellNumber; 
  uint p = (uint) cellNumber/2;

  for(uint i=1;i<=str.atoms.size();i++){
    for(uint j=1;j<=superStr.atoms.size();j++){
      d = superStr.atoms[p*(k*k+k+1)+(i-1)*k*k*k].cpos - superStr.atoms[j-1].cpos;
      if((gapPositions(i)+new_tol) >= modulus(d) && modulus(d)>_ZERO_TOL_LOOSE_){
        environment[i-1].push_back(j);
      }
    }
  }
}
//************** Function for building the AEs (end)

//************** Function for rebuilding failed AEs (beginning)
void BuildNewAE(int i, deque<deque<int> > &environment, xvector<double> gapPositions, xstructure superStr,
    int cellNumber, double new_tol){
  xvector<double> d;
  uint k = (uint) cellNumber; 
  uint p = (uint) cellNumber/2;

  environment[i-1].clear(); //Erases old AE

  for(uint j=1;j<=superStr.atoms.size();j++){
    d = superStr.atoms[p*(k*k+k+1)+(i-1)*k*k*k].cpos - superStr.atoms[j-1].cpos;
    if((gapPositions(i)+new_tol) >= modulus(d) && modulus(d)>_ZERO_TOL_LOOSE_){
      environment[i-1].push_back(j);
    }
  }
}
//************** Function for rebuilding failed AEs (end)

//////////////Create proper AE container(begin)
void DetermineDistinctAEs(vector<vector<xmatrix<int> > > VatomicEnvironments, vector<vector<vector<int> > > &distinctAE,
    vector<vector<vector<int> > > &AEcoefs, vector<int> &str_remove){

  xmatrix<int> str;
  uint num_species=VatomicEnvironments[0].size();

  for(uint i=0;i<VatomicEnvironments.size();i++){
    for(uint ns=0;ns<num_species;ns++){
      str=VatomicEnvironments[i][ns];
      if(str(1,1)!=0 && str(1,2)==0){
        str_remove[i]=1;
        break;
      }
    }
  }

  vector<int> c_temp;
  for(uint ns=0;ns<num_species;ns++){
    for(uint i=0;i<VatomicEnvironments.size();i++){
      str=VatomicEnvironments[i][ns];
      if(str(1,1)!=0 && str_remove[i]==0){
        vector<int> AE1;
        for(int j=0;j<str(2,1)*3;j++){
          AE1.push_back(str(4,j+1));
        }
        vector<vector<int> > dAEt;
        AEcoefs.push_back(dAEt);
        dAEt.push_back(AE1);
        distinctAE.push_back(dAEt);
        break;
      }
    }
  }

  for(uint ns=0;ns<num_species;ns++){
    for(uint i=0;i<VatomicEnvironments.size();i++){
      if(str_remove[i]==0){
        str=VatomicEnvironments[i][ns];
        c_temp.assign(distinctAE[ns].size(),0);
        if(str(1,1)!=0){
          for(int j=0;j<str(1,2);j++){
            for(uint l=0;l<distinctAE[ns].size();l++){
              int similarVertices=0;
              if(str(2,j+1)==(int)distinctAE[ns][l].size()/3){ //CO20190329
                for(int k=0;k<str(2,j+1);k++){
                  for(uint m=0;m<distinctAE[ns][l].size()/3;m++){
                    if(str((4+j),(3*k+1))==distinctAE[ns][l][3*m] && str((4+j),(3*k+2))==distinctAE[ns][l][3*m+1] && str((4+j),(3*k+3))==distinctAE[ns][l][3*m+2]){
                      similarVertices++;
                    }
                  }
                }
              }
              if(similarVertices==str(2,j+1)){
                c_temp[l]=c_temp[l]+str(3,j+1);
                break;
              }
              else if (l==distinctAE[ns].size()-1){ //new AE
                vector<int> AEtemp;
                for(int n=0;n<str(2,j+1)*3;n++){
                  AEtemp.push_back(str((4+j),(n+1)));
                }
                distinctAE[ns].push_back(AEtemp);
                c_temp.push_back(0);
                for(uint n=0;n<AEcoefs[ns].size();n++){
                  AEcoefs[ns][n].push_back(0);
                }
              }
            }
          }
        }
        AEcoefs[ns].push_back(c_temp);
      }
    }
  }
}

///////////////Create proper AE container(end)

namespace pflow {

  deque<deque<int> > CalcAtomicEnvironment(xstructure& str) {

    CalculateSymmetry(str);
    double new_tol = str.sym_eps;

    //Determines the size of the supercell
    int cellNumber = (int) pow((_MIN_NUM_ATM_/str.atoms.size()),1.0/3.0);
    if((pow((_MIN_NUM_ATM_/str.atoms.size()),1.0/3.0)-(double)cellNumber)>_ZERO_TOL_LOOSE_)
      cellNumber++;
    if(cellNumber<3)
      cellNumber=3;
    if(cellNumber%2==0) //Must always be odd
      cellNumber++;

    //Determines the size of the supercell

    // Creates supercell
    xmatrix<double> supercell(1,1,3,3);
    supercell(1,1)=supercell(2,2)=supercell(3,3)=cellNumber;
    xstructure superStr=GetSuperCell(str, supercell);


    //////////////////
    // For k (odd) unit cells in each direction, the central cell equivalent of atom "i" in the unit cell is going to be
    // atom "p*k²+p*k+p+1+(i-1)*k³" = "p(k²+k+1)+1+(i-1)*k³" with p = (int) k/2
    //////////////////

    uint k = (uint) cellNumber; 
    uint p = (uint) cellNumber/2;

    for(uint i=1; i<=str.atoms.size(); i++) //Just a check
      if(str.atoms[i-1].type!=superStr.atoms[(p*(k*k+k+1)+(i-1)*k*k*k)].type){cout << "WARNING!!! Mapping unsuccessful!" << endl;}

    bool samplingOK;
    xmatrix<double> distances = GetAllDistances(samplingOK, str, superStr, cellNumber);
    while(!samplingOK){
      cellNumber+=2;
      k = (uint) cellNumber; 
      p = (uint) cellNumber/2;
      supercell(1,1)=supercell(2,2)=supercell(3,3)=cellNumber;
      superStr=GetSuperCell(str, supercell);
      cout << "Insufficient sampling." << endl;
      cout << "Now using a " << cellNumber << "x" << cellNumber << "x" << cellNumber << " supercell" << endl;
      distances = GetAllDistances(samplingOK, str, superStr, cellNumber);
    }

    xvector<double> maxGaps(str.atoms.size(),1), gapPositions(str.atoms.size(),1);
    FindGap(maxGaps, gapPositions, distances, str, superStr, new_tol);

    deque<deque<int> > environment, dummyIDD;
    deque<int> dummyID;
    for(uint i=0;i<str.atoms.size();i++)
      environment.push_back(dummyID);
    BuildAE(environment, gapPositions, str, superStr, cellNumber, new_tol);

    uint maxAE=0;
    for(uint i=0; i<str.atoms.size(); i++)
      if(maxAE<environment[i].size())
        maxAE=environment[i].size();

    deque<deque<deque<int> > > aEConnectivity;
    for(uint i=0;i<str.atoms.size();i++){
      aEConnectivity.push_back(dummyIDD);
      for(uint j=0;j<maxAE;j++){	
        aEConnectivity[i].push_back(dummyID);
        for(uint k=0;k<maxAE;k++){
          aEConnectivity[i][j].push_back(0);
        }
      }
    }

    double dist=9999.666;
    xvector<double> atomicDistances(1, maxAE), d, v1, v2, v3, vp;
    deque<deque<int> > underCoord;
    for(uint i=0;i<str.atoms.size();i++)
      underCoord.push_back(dummyID);
    xvector<int> recalcAE(1,str.atoms.size());
    int vecNum=1, sumRecalcAE=1, neighCount=0;
    bool pyramid = false;

    set(recalcAE,1);

    while(sumRecalcAE!=0){

      for(uint i=1;i<=str.atoms.size();i++){//Finds connectivity and raises flag in case of maximum gap rule failure 
        if(recalcAE[i]==1){                                                                                          
          recalcAE[i]=0;                                                                                              

          for(uint j=1;j<=environment[i-1].size();j++){
            for(uint l=1;l<=environment[i-1].size();l++){
              d = superStr.atoms[environment[i-1][j-1]-1].cpos - superStr.atoms[environment[i-1][l-1]-1].cpos;
              if(modulus(d)<dist && j!=l)
                dist=modulus(d);
            }
            for(uint l=1;l<=environment[i-1].size();l++){ //Build connectivity matrix
              d = superStr.atoms[environment[i-1][j-1]-1].cpos - superStr.atoms[environment[i-1][l-1]-1].cpos;
              if(modulus(d)<=(dist*_CON_TOL_) && j!=l)
                aEConnectivity[i-1][j-1][l-1]=1;
            }
            dist=9999.666;
          }

          for(uint j=1;j<=environment[i-1].size();j++){ //Symmetrizes connectivity matrix
            for(uint l=1;l<=environment[i-1].size();l++){
              if(aEConnectivity[i-1][j-1][l-1]!=aEConnectivity[i-1][l-1][j-1]){
                aEConnectivity[i-1][j-1][l-1]=aEConnectivity[i-1][l-1][j-1]=1;
              }
            }
          }

          for(uint j=1;j<=environment[i-1].size();j++){ //Flags undercoordinated atoms (only the ones with less than 3 neighbors)
            neighCount=0;
            for(uint l=1;l<=environment[i-1].size();l++)
              neighCount+=aEConnectivity[i-1][j-1][l-1];
            if(neighCount<3)
              underCoord[i-1].push_back(j);
          }

          for(uint j=1;j<=environment[i-1].size();j++){
            for(uint l=1;l<=environment[i-1].size();l++){ //Checks if there are atoms over polyhedra surfaces
              if(aEConnectivity[i-1][j-1][l-1]==1){//j->l
                if(vecNum==1){//first connection
                  v1 = superStr.atoms[environment[i-1][j-1]-1].cpos - superStr.atoms[environment[i-1][l-1]-1].cpos;
                  vecNum=2;
                }
                else if(vecNum==2){//second connection
                  v2 = superStr.atoms[environment[i-1][j-1]-1].cpos - superStr.atoms[environment[i-1][l-1]-1].cpos;
                  vecNum=3;
                }
                else {//third connection
                  v3 = superStr.atoms[environment[i-1][j-1]-1].cpos - superStr.atoms[environment[i-1][l-1]-1].cpos;
                  vp = vector_product(vector_product(v1,v2),vector_product(v2,v3));
                  if(modulus(vp)<_ZERO_TOL_LOOSE_) //v1,v2 and v2,v3 define the same plane
                    recalcAE[i]=1; //Flag for recalculating AE
                }
              }
            }
            vecNum=1;
          }
        }

      }

      //Finds connectivity and raises flag in case of maximum gap rule failure
      sumRecalcAE=sum(recalcAE);
      if(sumRecalcAE!=0){ //Recalculates failed AEs
        for(uint i=1;i<=str.atoms.size();i++){
          if(recalcAE[i]==1){
            FindNewGap(i, maxGaps, gapPositions, distances, superStr, new_tol);
            BuildNewAE(i, environment, gapPositions, superStr, cellNumber, new_tol);
            for(uint j=1;j<=maxAE;j++){ //Resets connectivity matrix
              for(uint l=1;l<=maxAE;l++){
                aEConnectivity[i-1][j-1][l-1]=0;
              }
            }
          }
        }
      }

    }

    xmatrix<int> polyhedra(1,1,str.atoms.size(),4); //keeps # of: [1]atoms, [2]triangles, [3]rectangles, [4]higher ones (probably useless)
    deque<deque<deque<int> > > vertices;
    for(uint i=0;i<str.atoms.size();i++){
      vertices.push_back(dummyIDD);
      for(uint j=0;j<maxAE;j++){	
        vertices[i].push_back(dummyID);
        for(uint k=0;k<2;k++){
          vertices[i][j].push_back(0);
        }
      }
    }

    for(uint i=1;i<=str.atoms.size();i++)
      for(uint j=1;j<=maxAE;j++){
        vertices[i-1][j-1][0]=0;
        vertices[i-1][j-1][1]=0;
      }

    for(uint i=1;i<=str.atoms.size();i++){ //Categorizes each AE

      polyhedra(i,1)=environment[i-1].size();

      for(uint j=1;j<=environment[i-1].size();j++){
        for(uint l=1;l<=environment[i-1].size();l++){
          if(aEConnectivity[i-1][j-1][l-1]==1){ //j linked to l

            for(uint m=1;m<=environment[i-1].size();m++){
              if(aEConnectivity[i-1][l-1][m-1]==1 && aEConnectivity[i-1][j-1][m-1]==1){//Triangle (j->l->m->j)
                polyhedra(i,2)++;
                vertices[i-1][j-1][0]++;
              }
            }

            for(uint m=1;m<=environment[i-1].size();m++){
              for(uint n=1;n<=environment[i-1].size();n++){
                if(aEConnectivity[i-1][l-1][m-1]==1 && aEConnectivity[i-1][m-1][n-1]==1 
                    && aEConnectivity[i-1][n-1][j-1]==1 && aEConnectivity[i-1][m-1][j-1]!=1 && aEConnectivity[i-1][n-1][l-1]!=1 && 
                    j!=m && l!=n){ //Rectangle (j->l->m->n->j)
                  for(uint o=1;o<=environment[i-1].size();o++){
                    if(aEConnectivity[i-1][l-1][o-1]==1 && aEConnectivity[i-1][m-1][o-1]==1 && 
                        aEConnectivity[i-1][n-1][o-1]==1 && aEConnectivity[i-1][j-1][o-1]==1){ //The tetragon is the base of a pyramid
                      pyramid=true;
                    }
                  }
                  if(!pyramid){
                    polyhedra(i,3)++;
                    vertices[i-1][j-1][1]++;
                  }
                  else {
                    pyramid=false;
                  }
                }
              }
            }
            //No code implemented for higher symmetry polygons yet (probably not necessary)
          }
        }
        //Removes duplicity
        vertices[i-1][j-1][0]=vertices[i-1][j-1][0]/2;
        vertices[i-1][j-1][1]=vertices[i-1][j-1][1]/2;
      }
      //Removes duplicity
      polyhedra(i,2)=polyhedra(i,2)/6; //Factor 2 for clockwise and counter clockwise
      polyhedra(i,3)=polyhedra(i,3)/8; //Factor n for the number of vertices
    }

    //Categorizes each AE

    bool equalAE=false;
    deque<deque<int> >vertType; 
    for(uint i=1;i<=str.atoms.size();i++)
      vertType.push_back(dummyID);

    for(uint i=1;i<=str.atoms.size();i++){ //Identifies and groups vertices
      if(environment[i-1].size()!=0){
        for(uint j=1;j<=(environment[i-1].size()-1);j++){
          for(uint l=1;l<=(environment[i-1].size()-j);l++){
            if(vertices[i-1][j-1][0]==vertices[i-1][j+l-1][0] && vertices[i-1][j-1][1]==vertices[i-1][j+l-1][1]){
              equalAE=true;
            }
          }
          if(!equalAE){
            vertType[i-1].push_back(0);
            vertType[i-1].push_back(vertices[i-1][j-1][0]);
            vertType[i-1].push_back(vertices[i-1][j-1][1]);
          }
          equalAE=false;
        }
      }
      else {
        cout << "no environment to check" << endl;
        continue;
      }

      vertType[i-1].push_back(0); //The last vertex is not run by the above loop
      vertType[i-1].push_back(vertices[i-1][environment[i-1].size()-1][0]);
      vertType[i-1].push_back(vertices[i-1][environment[i-1].size()-1][1]);
      for(uint j=1;j<=environment[i-1].size();j++){ //Counts the #occurrences of each vertex type
        for(uint l=1;l<=(vertType[i-1].size()/3);l++){
          if(vertices[i-1][j-1][0]==vertType[i-1][1+(l-1)*3] && vertices[i-1][j-1][1]==vertType[i-1][2+(l-1)*3]){
            vertType[i-1][(l-1)*3]++;
          }
        }
      }
    }

    return vertType;
  }

  xmatrix<int> GetAtomicEnvironment(xstructure& str, int ns, deque<deque<int> > vertType) {

    bool equalAE=false;
    deque<int> dummyID;
    uint aux=0, distinctAE=0;

    int i_inc=0;
    if(ns==0){i_inc=0;}
    else {
      for(int i=0;i<ns;i++){
        i_inc=i_inc+str.num_each_type.at(i);
      }
    }

    for(int i=(1+i_inc);i<=(str.num_each_type.at(ns)+i_inc-1);i++){
      for(uint j=1;j<=(str.atoms.size()-i);j++){
        if(vertType[i-1].size()==vertType[i+j-1].size()){//AE may be similar
          for(uint l=1;l<=(vertType[i-1].size()/3);l++)
            for(uint m=1;m<=(vertType[i+j-1].size()/3);m++)
              if(vertType[i-1][1+(l-1)*3]==vertType[i+j-1][1+(m-1)*3] && vertType[i-1][(l-1)*3]==vertType[i+j-1][(m-1)*3] &&
                  vertType[i-1][2+(l-1)*3]==vertType[i+j-1][2+(m-1)*3]){
                aux = aux+3;
              }
          if(aux==vertType[i-1].size()){
            equalAE=true;
          }
          aux=0;
        }
      }
      if(!equalAE)
        distinctAE++;
      equalAE=false;
    }
    distinctAE++; //Accounts for the last AE which is not run through by the loop above

    deque<deque<int> > systemAE;
    for(uint i=0;i<=distinctAE;i++)
      systemAE.push_back(dummyID);           //First deque contains global information, the second onwards contain info on each
    systemAE[0].push_back(str.atoms.size()); //specific AE (the number of a type of vertex followed by the numbers of triangles
    //and squares meeting there    
    uint counter = 1;

    for(int i=(1+i_inc);i<=(str.num_each_type.at(ns)+i_inc);i++){ //Gets each AE properties
      for(uint j=1;j<=(str.atoms.size()-i);j++){
        if(vertType[i-1].size()==vertType[i+j-1].size()){//AE may be similar
          for(uint l=1;l<=(vertType[i-1].size()/3);l++){
            for(uint m=1;m<=(vertType[i+j-1].size()/3);m++){
              if(vertType[i-1][1+(l-1)*3]==vertType[i+j-1][1+(m-1)*3] && vertType[i-1][(l-1)*3]==vertType[i+j-1][(m-1)*3] &&
                  vertType[i-1][2+(l-1)*3]==vertType[i+j-1][2+(m-1)*3]){
                aux = aux+3;
              }
            }
          }
          if(aux==vertType[i-1].size()){
            equalAE=true;
          }
          aux=0;
        }
      }
      if(!equalAE || i==str.num_each_type.at(ns)+i_inc){
        for(uint l=1;l<=(vertType[i-1].size()/3);l++){
          systemAE[counter].push_back(vertType[i-1][(l-1)*3]);
          systemAE[counter].push_back(vertType[i-1][1+(l-1)*3]);
          systemAE[counter].push_back(vertType[i-1][2+(l-1)*3]);
        }
        counter++;
      }
      equalAE=false;
    }

    bool undercoord;
    uint AEcount=0;
    deque<deque<int> > AEkeep;
    AEkeep.push_back(systemAE[0]);
    for(uint i=1;i<distinctAE+1;i++){
      undercoord=false;
      for(uint j=0;j<(systemAE[i].size()/3);j++){
        if((systemAE[i][3*j+1]+systemAE[i][3*j+2])<3){
          undercoord=true;
          break;
        }
      }
      if(undercoord==false){
        AEcount++;
        AEkeep.push_back(systemAE[i]);
        //cout << "coordinated AE:" << endl;
        //for(uint j=0;j<systemAE[i].size();j++){
        // cout << systemAE[i][j] << " ";
        //}
        //cout << endl;
      }
      //else {//cout <<"undercoordinated AE:" << endl;
      //for(uint j=0;j<systemAE[i].size();j++){
      //  cout << systemAE[i][j] << " ";
      //}
      //cout << endl;
      //}
    }
    uint col=1;

    AEkeep[0].push_back(AEcount);
    uint maxVert;
    if(AEcount==0){
      cout << "ERROR calculating AE: vertices are undercoordinated." << endl;
      maxVert=1;
    }
    else {
      systemAE[0].push_back(AEcount);       
      maxVert=0;
      for(uint i=1;i<=AEcount;i++)
        if((AEkeep[i].size()/3)>maxVert)
          maxVert=AEkeep[i].size()/3;

      for(uint j=1;j<=AEcount;j++)
        AEkeep[j].push_back(0);

      for(int i=(1+i_inc);i<=(str.num_each_type.at(ns)+i_inc);i++){ //Counts AE occurrences
        for(uint j=1;j<=AEcount;j++){
          if(vertType[i-1].size()==(AEkeep[j].size()-1)){//AE may be similar
            for(uint l=1;l<=(vertType[i-1].size()/3);l++){
              for(uint m=1;m<=(AEkeep[j].size()/3);m++){
                if(vertType[i-1][1+(l-1)*3]==AEkeep[j][1+(m-1)*3] && vertType[i-1][(l-1)*3]==AEkeep[j][(m-1)*3] &&
                    vertType[i-1][2+(l-1)*3]==AEkeep[j][2+(m-1)*3]){
                  aux = aux+3;
                }
              }
            }
            if(aux==vertType[i-1].size()){
              AEkeep[j].back()++;
            }
            aux=0;
          }
        }
      }
    }
    if(AEcount>(3*maxVert)){
      col=AEcount;
    }
    else {
      col=3*maxVert;
    }

    xmatrix<int> environmentSet(1,1,(AEcount+3),col);
    environmentSet(1,1)=str.num_each_type.at(ns);
    environmentSet(1,2)=AEcount;
    for(uint i=1;i<=AEcount;i++){
      environmentSet(2,i)=AEkeep[i].size()/3;
      environmentSet(3,i)=AEkeep[i].back();
      for(uint j=1;j<=(AEkeep[i].size()-1);j++){
        environmentSet(i+3,j)=AEkeep[i][j-1];
      }
    }

    return environmentSet;

  }


  //************** Function for computing GFA (beginning)
  void ComputeGFA(vector<double> Egs, vector<vector<double> > VDecomp_coefs, vector<double>& sampling, vector<double>& gfa,
      vector<vector<uint> > VDecomp_phases, vector<xmatrix<double> > xAEcoefs, uint dimension, uint entries_size,
      vector<uint> distinctAE_size, vector<chull::ChullPoint> vcpt, vector<vector<double> > VStoichE,
      vector<string> species){

    bool LDEBUG=(FALSE || XHOST.DEBUG || _DEBUG_GFA_); //CO20190424
    if(LDEBUG) { //CO20190424
      for(uint i=0;i<Egs.size();i++){cerr << __AFLOW_FUNC__ << " Egs[i=" << i << "]=" << Egs[i] << endl;} //CO20190424
    } //CO20190424

    aurostd::xcombos xc(entries_size, dimension);
    vector<vector<int> > indices;
    while(xc.increment()){
      vector<int> ti = xc.getIndices();
      bool equal_stoich=false;
      for(uint i=0;i<ti.size()-1;i++){
        for(uint j=i+1;j<ti.size();j++){
          for(uint x=0;x<dimension-1;x++){
            if(abs(VStoichE.at(ti[i]).at(x)-VStoichE.at(ti[j]).at(x))>0.000001){
              break;
            }
            else if(x==dimension-2){
              equal_stoich=true;
            }
          }
          if(equal_stoich==true){
            break;
          }
        }
      }
      if(equal_stoich==false){
        indices.push_back(ti);
      }
    }
    cout << "number of possible combinations: " << indices.size() << endl << endl;
    double gfa_max=0; double Tweight=0; 
    vector<double> gstoich_max;

    for(uint X=0;X<vcpt.size();X++){
      vector<double> gstoich; double last=1;
      for(int i=vcpt[X].m_coords.lrows;i<=vcpt[X].m_coords.urows-1;i++){  //CO20190424
        gstoich.push_back(vcpt[X].m_coords[i]);
        last=last-vcpt[X].m_coords[i];
      }
      gstoich.push_back(last);
      cout << "stoichiometry " << X << ": ";
      for(uint i=0;i<gstoich.size();i++){
        cout << gstoich[i] << " ";
      }
      cout << endl << "formation enthalpy of the ground state: " << Egs[X] << " eV" << endl;

      vector<vector<xvector<double> > >Vpseudo_AEcoefs;
      vector<xvector<double> > VpseudoGS_AEcoefs;
      Vpseudo_AEcoefs.resize(dimension);

      int dim_div=dimension;

      for(uint ns=0;ns<dimension;ns++){
        if(gstoich[ns]<0.0000001){dim_div=dim_div-1;}
        xvector<double> ps_temp(distinctAE_size[ns],1);
        VpseudoGS_AEcoefs.push_back(ps_temp);
      }
      xvector<double> pscg_temp(entries_size,1);
      vector<double> weight_avgDP, VEavg;
      for(uint p=0;p<VDecomp_phases[X].size();p++){
        if(VDecomp_phases[X][p]>entries_size-1){cerr << "ERROR in convex hull calculation: decomposition phase is not an entry." << endl;}
        pscg_temp[VDecomp_phases[X][p]+1]=VDecomp_coefs[X][p];
      }
      for(uint ns=0;ns<dimension;ns++){
        xvector<double> psg_temp=aurostd::operator*(pscg_temp,xAEcoefs[ns]);
        double modulus=aurostd::modulus(psg_temp);

        if(modulus==0){
          for(uint a=1;a<distinctAE_size[ns]+1;a++){
            VpseudoGS_AEcoefs[ns][a]=0;
          }
        }
        else {
          for(uint a=1;a<distinctAE_size[ns]+1;a++){
            VpseudoGS_AEcoefs[ns][a]=psg_temp[a]/modulus;
            if(VpseudoGS_AEcoefs[ns][a]>1.0000001){cout << "ERROR: a=" << a << ", VpseudoGS_AEcoefs[ns][a]=" << VpseudoGS_AEcoefs[ns][a] << endl;}
          }
        }
      }
      vector<int> stoich_indices;

      for(uint l=0;l<indices.size();l++){

        bool may_balance=true; vector<vector<int> > side_l;
        for(uint i=0;i<dimension;i++){
          vector<int> side(dimension-1,2);
          for(uint x=0;x<dimension;x++){
            if(gstoich[x]-VStoichE.at(indices[l][i]).at(x)<-0.000001){
              side[x]=0;
            }
            else if(gstoich[x]-VStoichE.at(indices[l][i]).at(x)>0.000001){
              side[x]=1;
            }
          }
          side_l.push_back(side);
        }
        for(uint i=0;i<dimension-1;i++){
          for(uint x=0;x<dimension;x++){
            if(side_l[i][x]!=side_l[i+1][x]){
              break;
            }
            else if(i==dimension-2){
              may_balance=false;
            }
          }
        }

        if(may_balance==true){

          double sigma_r=0.1;

          //CO20190509
          //heuristic, so calculation is feasible for ternaries
          int sigma_weight = 3;
          if(dimension==2){sigma_weight=300;} //include everything

          vector<double> dist;
          bool stop=false;
          for(uint i=0;i<dimension;i++){
            double dt=0;
            for(uint x=0;x<dimension;x++){
              dt=dt+pow(VStoichE.at(indices[l][i]).at(x)-gstoich[x],2);
            }
            dt=sqrt(dt);
            //[CO20190501 - OBSOLETE]if(dt>300*sigma_r)
            if(dt>sigma_weight*sigma_r) //CO20190501
            { //CO20200106 - patching for auto-indenting
              stop=true;
              break;
            }
            else {
              dist.push_back(dt);
            }
          }

          if(stop==false){

            vector<xvector<double> > lhs, rhs; xvector<double> compositionL(dimension); //-1,0);  //CO20190424
            for (uint x=0;x<dimension;x++){
              compositionL[x+compositionL.lrows]=gstoich[x];  //CO20190424
            }
            lhs.push_back(compositionL);
            for(uint i=0;i<dimension;i++){
              xvector<double> compositionR(dimension);  //-1,0);  //CO20190424
              for(uint x=0;x<dimension;x++){
                compositionR[x+compositionR.lrows]=VStoichE.at(indices[l][i]).at(x);  //CO20190424
              }
              rhs.push_back(compositionR);
            }

            xvector<double> stoich_coeffs;
            bool balanced=true;
            try{stoich_coeffs=balanceChemicalEquation(lhs,rhs,true,0.001);}
            catch(aurostd::xerror& vre){
              balanced=false;
              string check=vre.what();
              vector<string> tokens;
              aurostd::string2tokens(check, tokens," ");
              if(tokens[0]!="Found"){
                cout << vre.what() << endl;
              }
            }
            if (balanced==true){

              double weight=0, tdist=0;
              xvector<double> psc_temp(entries_size,1);

              for(uint i=0;i<dimension;i++){
                tdist=tdist+dist[i]*stoich_coeffs[i+2];
                psc_temp[indices[l][i]+1]=stoich_coeffs[i+2];
              }
              weight=weight+exp(-(tdist*tdist)/(2*sigma_r*sigma_r));

              Tweight=Tweight+weight;

              double dotProduct=0;
              for(uint ns=0;ns<dimension;ns++){
                xvector<double> ps_temp(distinctAE_size[ns],1);
                if(gstoich[ns]>0.0000001){
                  for(uint a=1;a<distinctAE_size[ns]+1;a++){
                    for(uint i=0;i<dimension;i++){
                      if(psc_temp[indices[l][i]+1]>0.0000001 && xAEcoefs[ns][indices[l][i]+1][a]>0.0000001){
                        ps_temp[a]=ps_temp[a]+psc_temp[indices[l][i]+1]*xAEcoefs[ns][indices[l][i]+1][a];
                      }
                    }
                  }

                  double modulus=aurostd::modulus(ps_temp);
                  if(modulus<0.0000001){
                    for(uint a=1;a<distinctAE_size[ns]+1;a++){
                      ps_temp[a]=0;
                    }
                  }
                  else {
                    for(uint a=1;a<distinctAE_size[ns]+1;a++){
                      ps_temp[a]=ps_temp[a]/modulus;
                      if(ps_temp[a]>1.0000001){cout << "ERROR: l=" << l << ", ns=" << ns << ", a=" << a << ", ps_temp[a]=" << ps_temp[a] << endl;}
                    }
                  }

                  double dotProduct_temp=aurostd::scalar_product(ps_temp,VpseudoGS_AEcoefs[ns])*gstoich[ns];
                  if(abs(dotProduct_temp)<0.0000001){dotProduct_temp=0;}
                  dotProduct=dotProduct+dotProduct_temp;
                } 
                Vpseudo_AEcoefs[ns].push_back(ps_temp);
              }

              if(dotProduct>1.0000001 || dotProduct<-0.0000001){cout << "ERROR: dotProduct=" << dotProduct << endl;}
              if(abs(dotProduct)<0.0000001){dotProduct=0;}

              double El_ref=0;
              for(uint i=0;i<dimension;i++){
                if(psc_temp[indices[l][i]+1]>0.0000001){
                  El_ref=El_ref+VStoichE.at(indices[l][i]).back()*psc_temp[indices[l][i]+1];
                }
              }
              double gfa_tmp = exp(-abs(El_ref-Egs[X])/kBT_ROOM)*weight*(1.0-dotProduct); //DX20210122 - create variable
              if(LDEBUG){ cerr << __AFLOW_FUNC__ << " exp(-abs(El_ref-Egs[X])/kBT)*weight*(1.0-dotProduct): " << gfa_tmp << endl; }

              gfa[X]+=gfa_tmp;
              if(LDEBUG){ cerr << __AFLOW_FUNC__ << " gfa[X]: " << gfa[X] << endl; }
              sampling[X]=sampling[X]+weight;

              weight_avgDP.push_back(weight);
              if(weight_avgDP.back()<-0.0000001){cout << "ERROR: weight_avgDP.back() = " << weight_avgDP.back() << ", weight = " << weight << endl;}
              VEavg.push_back(El_ref/kBT_ROOM);

            }	  
          }
        }
      }
      if(sampling[X]==0){gfa[X]=0;}
      else {gfa[X]=gfa[X]*100/sampling[X];}

      cout << "GFA before avgDP: " << gfa[X] << ", sum of weights: " << sampling[X] << endl;

      double avgDP=0.0;

      if(gfa[X]>0.0000001){
        double count=0;
        for(uint ns=0;ns<dimension;ns++){
          double count_ns=0, TDP_ns=0;
          if(gstoich[ns]>0.0000001){
            if(Vpseudo_AEcoefs[ns].size()>0){
              vector<xvector<double> > distinctPS;
              vector<double> dPS_coefs, dPS_weights;
              distinctPS.push_back(Vpseudo_AEcoefs[ns][0]); dPS_coefs.push_back(0); dPS_weights.push_back(weight_avgDP[0]);
              cout << "species: " << species[ns] << ", number of pseudostructures: " << Vpseudo_AEcoefs[ns].size() << endl;
              for(uint l=0;l<Vpseudo_AEcoefs[ns].size();l++){
                bool equal_pseudo=false;
                for(uint m=0;m<distinctPS.size();m++){
                  if(abs(weight_avgDP[l]-dPS_weights[m])<0.0000001){
                    for(uint a=1;a<distinctAE_size[ns]+1;a++){
                      if(abs(Vpseudo_AEcoefs[ns][l][a]-distinctPS[m][a])>0.0000001){
                        break;
                      }
                      else if(a==distinctAE_size[ns]){
                        equal_pseudo=true;
                      }
                    }
                  }
                  if(equal_pseudo==true){
                    dPS_coefs[m]=dPS_coefs[m]+1;
                    break;
                  }
                  else if(m==distinctPS.size()-1){
                    distinctPS.push_back(Vpseudo_AEcoefs[ns][l]);
                    dPS_coefs.push_back(0);
                    dPS_weights.push_back(weight_avgDP[l]);
                    if(dPS_weights.back()<-0.0000001){cout << "ERROR: l=" << l << ", weight_avgDP[l] = " << weight_avgDP[l]
                      << ", dPS_weights.back() = " << dPS_weights.back() << endl;}
                  }
                }
              }

              cout << "number of distinct pseudostructures: " << distinctPS.size() << endl;

              for(uint d1=0;d1<distinctPS.size()-1;d1++){
                for(uint d2=d1+1;d2<distinctPS.size();d2++){
                  double DP=0.0, d12=dPS_coefs[d1]*dPS_coefs[d2]*dPS_weights[d1]*dPS_weights[d2];
                  if(dPS_coefs[d1]<0){cout << "ERROR: dPS_coefs[d1] = " << dPS_coefs[d1] << endl;}
                  else if(dPS_coefs[d2]<0){cout << "ERROR: dPS_coefs[d2] = " << dPS_coefs[d2] << endl;}
                  else if(dPS_weights[d1]<0){cout << "ERROR: dPS_weights[d1] = " << dPS_weights[d1] << endl;}
                  else if(dPS_weights[d2]<0){cout  << "ERROR: dPS_weights[d2] = "  << dPS_weights[d2] << endl;}
                  for(uint a=1;a<distinctAE_size[ns]+1;a++){
                    if(distinctPS[d1][a]>0.0000001 && distinctPS[d2][a]>0.0000001){
                      DP=DP+distinctPS[d1][a]*distinctPS[d2][a];
                    }
                  }
                  TDP_ns=TDP_ns+(1-DP)*d12;
                  count_ns=count_ns+d12;
                }
              }

              for(uint d=0;d<distinctPS.size();d++){
                count_ns=count_ns+(dPS_coefs[d]*(dPS_coefs[d]-1)/2)*dPS_weights[d]*dPS_weights[d];
              }
            }
            count=count+count_ns; 
            avgDP=avgDP+TDP_ns*gstoich[ns]/count_ns;
          }
        }

      }

      else {cout << "Skipping avgDP" << endl;}

      gfa[X]=gfa[X]*avgDP*avgDP;

      cout << "avgDP: " << avgDP << ", final GFA=GFA*avgDP^2: " << gfa[X] << endl << endl;

      if(gfa[X]>gfa_max){
        gfa_max=gfa[X];
        gstoich_max=gstoich;
      }

    }
    cout << "maximum GFA: " << gfa_max << ", stoichiometry at maximum: ";
    for(uint ns=0;ns<gstoich_max.size();ns++){
      cout << gstoich_max[ns] << " ";
    }
    cout << endl;
  }

  //************** Function for computing GFA (end)

  void CalculateGFA(aurostd::xoption& vpflow, const string& alloy, const string& AE_file_read, double fe_cut){
    bool LDEBUG=(FALSE || XHOST.DEBUG || _DEBUG_GFA_); //CO20190424

    uint entries_size=0;
    string buffer;
    vector<double> sampling, gfa, tEForm, EForm, Egs;
    vector<string> tAUID, AUID;
    vector<vector<double> > stoichiometries, VDecomp_coefs, tVStoichiometry, VStoichiometry, VStoichE;
    vector<vector<uint> > VDecomp_phases;
    vector<vector<string> > tVSpecies, VSpecies;
    vector<vector<xstructure> > tVStructure, tVStructure2, VStructure;
    vector<vector<xmatrix<int> > > VatomicEnvironments;
    vector<aflowlib::_aflowlib_entry> LIB_entries;
    chull::ConvexHull the_hull;
    vector<chull::ChullPoint> Vpoint, vcpt;
    ostream& oss=cout;
    vector<string> species=aurostd::getElements(alloy,pp_string,true,true,false,oss);  //clean and sort, do not keep_pp //[CO20190712 - OBSOLETE]getAlphabeticVectorString(alloy);
    string input=aurostd::joinWDelimiter(species,""); //getAlphabeticString(alloy); //CO20190712
    vector<chull::CoordGroup> VCoordGroup;

    uint dimension=species.size();

    cout << endl << "Loading entries from LIBX. . ." << endl;

    vpflow.flag("PFLOW::LOAD_ENTRIES_LOAD_LIB2", true);
    if(dimension==3){
      vpflow.flag("PFLOW::LOAD_ENTRIES_LOAD_LIB3", true);
    }
    vpflow.flag("PFLOW::LOAD_ENTRIES_NARIES_MINUS_ONE", true);
    //[CO20190715 - LOAD_ENTRIES_ONLY_ALPHABETICAL -> LOAD_ENTRIES_NON_ALPHABETICAL]vpflow.flag("PFLOW::LOAD_ENTRIES_ONLY_ALPHABETICAL",true);  
    vpflow.flag("PFLOW::LOAD_ENTRIES_LOAD_XSTRUCTURES",true);

    bool quiet=XHOST.QUIET;
    XHOST.QUIET=true;
    loadEntries(vpflow,species,LIB_entries,cout);
    XHOST.QUIET=quiet;

    for(uint i=0;i<LIB_entries.size();i++){

      if(LIB_entries[i].enthalpy_formation_atom!=0){
        tEForm.push_back(LIB_entries[i].enthalpy_formation_atom);
      }
      else if(LIB_entries[i].enthalpy_formation_cell!=0){
        tEForm.push_back(LIB_entries[i].enthalpy_formation_cell/LIB_entries[i].natoms);
      }
      tVStructure.push_back(LIB_entries[i].vstr);
    }

    cout << endl << "Removing duplicates. . ." << endl;

    int entries_less_dup=0;
    for(uint i=0;i<LIB_entries.size();i++){
      bool equal_E=false;
      for(uint j=0;j<LIB_entries.size();j++){
        if(j!=i && LIB_entries[i].vstoichiometry[0]==LIB_entries[j].vstoichiometry[0] && LIB_entries[i].vspecies[0]==LIB_entries[j].vspecies[0]
            && abs(tEForm[j]-tEForm[i])<0.01 && (LIB_entries[i].lattice_system_relax=="" || LIB_entries[j].lattice_system_relax=="" ||
              LIB_entries[i].lattice_system_relax==LIB_entries[j].lattice_system_relax)){
          if(tEForm[j]==tEForm[i] && i<j){
            equal_E=true;
          }
          else if (tEForm[j]<tEForm[i]){
            if(compare::structuresMatch(tVStructure[i].back(),tVStructure[j].back(),true,false,true)==true){
              break;
            }
          }
        }
        if(j==LIB_entries.size()-1 && equal_E==false){
          entries_less_dup++;
          if(tEForm[i] < fe_cut){
            entries_size++;
            tAUID.push_back(LIB_entries[i].auid);
            VStoichiometry.push_back(LIB_entries[i].vstoichiometry);
            VSpecies.push_back(LIB_entries[i].vspecies);
            tVStructure2.push_back(LIB_entries[i].vstr);
            EForm.push_back(tEForm[i]);
          }
        }
      }
    }

    for(uint i=0;i<entries_size;i++){
      if(VSpecies[i].size()<species.size()){
        vector<double> vstoich_temp(species.size(),0);
        for(uint j=0;j<species.size();j++){
          for(uint k=0;k<VSpecies[i].size();k++){
            if(VSpecies[i][k]==species[j]){
              vstoich_temp[j]=VStoichiometry[i][k];
            }
          }
        }
        VStoichiometry[i]=vstoich_temp;
      }
    }
    stoichiometries.push_back(VStoichiometry[0]);
    for(uint i=1;i<entries_size;i++){
      bool GotIt=false;
      for(uint j=0;j<stoichiometries.size();j++){
        for(uint k=0;k<VStoichiometry[i].size();k++){
          if(VStoichiometry[i][k]!=stoichiometries[j][k]){
            break;
          }
          else if(k==VStoichiometry[i].size()-1){
            GotIt=true;
          }
        }
        if(GotIt==true){
          break;
        }
        else if(j==stoichiometries.size()-1){
          stoichiometries.push_back(VStoichiometry[i]);
          break;
        }
      }
    }

    sort(stoichiometries.begin(), stoichiometries.end());

    for (uint i=0;i<entries_size;i++){
      VStoichE.push_back(VStoichiometry[i]); VStoichE[i].push_back(EForm[i]); VStoichE[i].push_back(i);
    }
    sort(VStoichE.begin(), VStoichE.end());

    for(uint i=0;i<VStoichE.size();i++){
      for(uint j=0;j<tAUID.size();j++){
        if(VStoichE[i].back()==j){
          AUID.push_back(tAUID[j]);
          VStructure.push_back(tVStructure2[j]);
          break;
        }
      }
    }

    for(uint i=0;i<VStoichE.size();i++){
      VStoichE[i].pop_back();
    }

    cout << endl << "Getting atomic environments for " << entries_size << " entries. . ." << endl << endl;

    ofstream AEData;
    string AEfilename;

    AEfilename = "All_atomic_environments.dat";
    AEData.open(AEfilename.c_str());

    VatomicEnvironments.resize(entries_size);

    for(uint i=0;i<entries_size;i++){
      xstructure str=VStructure[i].back();

      // Read/write atomic environments file

      if(AE_file_read != "none"){
        uint count = 0, line_count = 0;
        stringstream AE_file_read_sstr;
        aurostd::file2stringstream(AE_file_read,AE_file_read_sstr);
        vector<string> AE_file_read_vstr;
        aurostd::string2vectorstring(AE_file_read_sstr.str(),AE_file_read_vstr);
        vector<vector<string> > Vtokens;
        for(uint ns=0;ns<species.size();ns++){
          Vtokens.clear();
          line_count=0;
          vector<string> tokens;
          for(uint j=0;j<AE_file_read_vstr.size();j++){
            line_count++;
            tokens.clear();
            string AE_line(AE_file_read_vstr[j]);
            stringstream AE_line_sstr(AE_line);
            while(AE_line_sstr >> buffer){tokens.push_back(buffer);}
            if(tokens[0] == "AUID" && tokens[2] == AUID[i] && tokens[5] == species.at(ns)){
              count = 0;
              for(uint k=j+1;k<AE_file_read_vstr.size();k++){
                tokens.clear();
                string AE_line2(AE_file_read_vstr[k]);
                stringstream AE_line2_sstr(AE_line2);
                while(AE_line2_sstr >> buffer){tokens.push_back(buffer);}
                if(tokens[0] == "AUID"){
                  break;
                }
                else {
                  Vtokens.push_back(tokens);
                  count++;
                }
              }
              break;
            }
          }

          if(line_count != AE_file_read_vstr.size()){
            int AEmc = Vtokens[0].size();
            xmatrix<int> ae_tmp(1,1,count,AEmc);
            for(uint k=0;k<count;k++){
              for(uint l=0;l<Vtokens[k].size();l++){
                uint m=k+1, n=l+1;
                ae_tmp[m][n]=aurostd::string2utype<int>(Vtokens[k][l]);
              }
            }
            VatomicEnvironments[i].push_back(ae_tmp);
            if(ns==species.size()-1){
              cout << "Read atomic environment for AUID " << AUID[i] << " from file" << endl;
            }
          }
        }
        if(line_count == AE_file_read_vstr.size()){
          cout << "Couldn't find AUID " << AUID[i] << " in AE file, calculating AE" << endl;
          str.ReScale(1.0);
          deque<deque<int> > vertType=CalcAtomicEnvironment(str);
          int nae=0;
          for(uint ns=0;ns<species.size();ns++){
            if(VStoichE[i][ns]>0.0000001){
              VatomicEnvironments[i].push_back(GetAtomicEnvironment(str, nae, vertType));
              nae++;
            }
            else {
              xmatrix<int> environmentSet(1,1,1,1);
              VatomicEnvironments[i].push_back(environmentSet);
            }
          }
        }
      }
      else {
        cout << "Calculating AE for AUID " << AUID[i] << endl;
        str.ReScale(1.0);
        deque<deque<int> > vertType=CalcAtomicEnvironment(str);
        int nae=0;
        for(uint ns=0;ns<species.size();ns++){
          if(VStoichE[i][ns]>0.0000001){
            VatomicEnvironments[i].push_back(GetAtomicEnvironment(str, nae, vertType));
            nae++;
          }
          else {
            xmatrix<int> environmentSet(1,1,1,1);
            VatomicEnvironments[i].push_back(environmentSet);
          }
        }
      }
      for(uint ns=0;ns<species.size();ns++){
        AEData << "AUID = " << AUID[i] << " atom = " << species.at(ns) << endl << VatomicEnvironments[i][ns] << endl;
      }
    }

    AEData.close();

    vector<int> str_remove (entries_size,0);
    vector<vector<vector<int> > > distinctAE, AEcoefs;

    cout << endl << "Determining distinct atomic environments. . ." << endl << endl;
    DetermineDistinctAEs(VatomicEnvironments, distinctAE, AEcoefs, str_remove);

    for(uint ns=0;ns<distinctAE.size();ns++){cout << "species: " << species[ns] << ", number of distinct atomic environments: " << distinctAE[ns].size() << endl;}

    deque<int> sremove;
    for(uint i=0;i<str_remove.size();i++){
      if(str_remove[i]==1){sremove.push_front(i);}
    }

    vector<uint> distinctAE_size;
    vector<xmatrix<double> > xAEcoefs;
    for(uint ns=0;ns<species.size();ns++){
      distinctAE_size.push_back(distinctAE[ns].size());
      xmatrix<int> xAEcoefsI=aurostd::vectorvector2xmatrix(AEcoefs[ns]);
      xAEcoefs.push_back(aurostd::xdouble(xAEcoefsI));
    }

    for(uint i=0;i<sremove.size();i++){
      AUID.erase(AUID.begin()+sremove[i]);
      VStoichE.erase(VStoichE.begin()+sremove[i]);
    }

    ofstream EntryData;
    string Entryfilename = "GFA_entries.dat";
    EntryData.open(Entryfilename.c_str());

    EntryData << endl << LIB_entries.size() << " entries loaded. " << endl << entries_less_dup << " entries remain after removing duplicates."
      << endl << entries_size << " entries have formation enthalpy < " << fe_cut << " eV." << endl;

    entries_size=VStoichE.size();

    EntryData << entries_size << " entries remain after atomic environments routine." << endl
      << endl << "stoichiometries and formation enthalpies for all structures used in the GFA calculation:" << endl;

    for(uint i=0;i<entries_size;i++){
      for(uint j=0;j<VStoichE[i].size()-1;j++){
        EntryData << VStoichE[i][j] << "  " ;
      }
      EntryData << VStoichE[i].back() << endl;
    }

    EntryData.close();

    if(dimension==2){
      xvector<double> vst(dimension); //1,0); //CO20190424
      for(double i=0;i<101;i++){
        vst[vst.lrows]=i/100; //CO20190424
        vst[vst.lrows+1]=2; //CO20190424
        chull::ChullPoint cpt(vst, cout, true, true, true);
        vcpt.push_back(cpt);
      }
    }
    else if(dimension==3){
      xvector<double> vst(dimension); //2,0); //CO20190424
      for(double i=0;i<21;i++){
        for(double j=0;j<21;j++){
          if(1-i/20-j/20>-0.00001){
            vst[vst.lrows]=i/20;  //CO20190424
            vst[vst.lrows+1]=j/20;  //CO20190424
            vst[vst.lrows+2]=2; //CO20190424
            chull::ChullPoint cpt(vst, cout, true, true, true);
            vcpt.push_back(cpt);
          }
        }
      }
    }
    vector<chull::ChullPoint> tVpoint;
    for(uint i=0;i<VStoichE.size();i++){
      xvector<double> vst(dimension); //-1,0);  //CO20190424
      for(int j=vst.lrows;j<=vst.urows-1;j++){ //CO20190424
        vst[j]=VStoichE[i][j-vst.lrows]; //CO20190424
      }
      vst[vst.urows]=VStoichE[i].back();  //CO20190424
      chull::ChullPoint cpt(vst, cout, true, true, false);
      tVpoint.push_back(cpt);  
    }

    cout << endl << "Calculating convex hull. . ." << endl;

    vpflow.flag("CHULL::INCLUDE_UNRELIABLE_HULLS",true);
    vpflow.flag("FORCE",true);

    XHOST.QUIET=true;
    the_hull = chull::ConvexHull(vpflow, tVpoint, cout, true, true);
    XHOST.QUIET=quiet;

    VCoordGroup = the_hull.m_coord_groups;
    Vpoint = the_hull.m_points;

    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " [1]" << endl;} //CO20190424

    for(uint i=0;i<stoichiometries.size();i++){
      bool equal=false;
      for(uint j=0;j<vcpt.size();j++){
        for(int k=vcpt[j].m_coords.lrows;k<=vcpt[j].m_coords.urows-1;k++){  //CO20190424
          if(vcpt[j].m_coords[k]!=stoichiometries[i][k-vcpt[j].m_coords.lrows]){  //CO20190424
            break;
          }
          else if ((k-vcpt[j].m_coords.lrows)==(int)dimension-2){  //CO20190424
            equal=true;
          }
        }
        if(equal==true){
          break;
        }
        else if(j==vcpt.size()-1){
          xvector<double> vst(dimension); //-1,0);  //CO20190424
          for(int k=vst.lrows;k<=vst.urows-1;k++){  //CO20190424
            vst[k]=stoichiometries[i][k-vst.lrows]; //CO20190424
          }
          vst[vst.urows]=2; //CO20190424
          chull::ChullPoint cpt(vst, cout, true, true, true);
          vcpt.push_back(cpt);
        }
      }
    }

    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " [2]" << endl;} //CO20190424

    sort(vcpt.begin(), vcpt.end());

    VDecomp_coefs.resize(vcpt.size());  VDecomp_phases.resize(vcpt.size());
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " [2.0]" << endl;} //CO20190424
    for(uint i=0;i<vcpt.size();i++){
      bool on_hull=false;
      if(LDEBUG) {cerr << __AFLOW_FUNC__ << " [2.1]" << endl;} //CO20190424
      for(uint j=0;j<VCoordGroup.size();j++){
        for(int k=vcpt[i].m_coords.lrows;k<=vcpt[i].m_coords.urows-1;k++){  //CO20190424
          if(abs(vcpt[i].m_coords[k]-VCoordGroup[j].m_coords[k])>0.00001){
            break;
          }
          else if((k-vcpt[i].m_coords.lrows)==(int)dimension-2 && VCoordGroup[j].m_is_on_hull){  //CO20190424
            on_hull=true;
            VDecomp_coefs[i].push_back(1.0); VDecomp_phases[i].push_back(VCoordGroup[j].m_ref_state);
            Egs.push_back(Vpoint.at(VCoordGroup[j].m_hull_member).getLastCoord());  //CO20190424
            if(LDEBUG) { //CO20190424
              cerr << __AFLOW_FUNC__ << " VCoordGroup[j].m_hull_member=" << VCoordGroup[j].m_hull_member << endl; //CO20190424
              cerr << __AFLOW_FUNC__ << " Egs.back()=" << Egs.back() << endl; //CO20190424
            } //CO20190424
          }
        }
      }
      if(LDEBUG) {cerr << __AFLOW_FUNC__ << " [2.2]" << endl;} //CO20190424
      if(on_hull==false){
        Egs.push_back(vcpt[i].getLastCoord()-the_hull.getDistanceToHull(vcpt[i]));  //CO20190424
        if(LDEBUG) { //CO20190424
          cerr << __AFLOW_FUNC__ << " i=" << i << endl; //CO20190424
          cerr << __AFLOW_FUNC__ << " vcpt[i].m_coords[dimension-1]=" << vcpt[i].getLastCoord() << endl; //CO20190424
          cerr << __AFLOW_FUNC__ << " the_hull.getDistanceToHull(vcpt[i])=" << the_hull.getDistanceToHull(vcpt[i]) << endl; //CO20190424
          cerr << __AFLOW_FUNC__ << " Egs.back()=" << Egs.back() << endl; //CO20190424
        } //CO20190424
        VDecomp_phases[i]=the_hull.getDecompositionPhases(vcpt[i]);
        vector<xvector<double> > lhs, rhs; double others=0; xvector<double> compositionL(VDecomp_phases[i].size()); //-1,0);  //CO20190424

        for(uint j=0;j<VDecomp_phases[i].size();j++){
          if(Vpoint[VDecomp_phases[i][j]].m_is_artificial){
            VDecomp_phases[i][j]=the_hull.m_coord_groups.at(Vpoint.at(VDecomp_phases[i][j]).m_i_coord_group).m_ref_state;
          }
        }
        int k=compositionL.lrows; //CO20190424
        for(int j=vcpt[i].m_coords.lrows;j<=vcpt[i].m_coords.urows-1;j++){  //CO20190424
          if(vcpt[i].m_coords[j]!=0){
            compositionL[k]=vcpt[i].m_coords[j];
            others=others+vcpt[i].m_coords[j];   
            k++;
          }
        }
        if(abs(1-others)>0.00001){compositionL[compositionL.urows]=1-others;} //CO20190424
        lhs.push_back(compositionL);

        for(uint j=0;j<VDecomp_phases[i].size();j++){
          others=0; /*k=compositionR.lrows;*/ xvector<double> compositionR(VDecomp_phases[i].size()); //-1,0); //CO20190424
          k=compositionR.lrows; //CO20190424
          for(int l=vcpt[i].m_coords.lrows;l<=vcpt[i].m_coords.urows-1;l++){  //CO20190424
            if(vcpt[i].m_coords[l]!=0){
              compositionR[k]=Vpoint.at(VDecomp_phases[i][j]).m_coords[l];
              others=others+Vpoint.at(VDecomp_phases[i][j]).m_coords[l];
              k++;
            }
          }
          if(abs(1-others)>0.00001){compositionR[compositionR.urows]=1-others;}  //CO20190424
          rhs.push_back(compositionR);
        }

        xvector<double> decomp_coeffs = balanceChemicalEquation(lhs,rhs,true,0.001);
        for(uint j=1;j<VDecomp_phases[i].size()+1;j++){
          VDecomp_coefs[i].push_back(decomp_coeffs[j+decomp_coeffs.lrows]); //CO20190424
        }
      }
    }

    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " [3]" << endl;} //CO20190424

    cout << endl << "Computing GFA for " << vcpt.size() << " stoichiometries. . ." << endl << endl;

    for(uint i=0;i<vcpt.size();i++){
      gfa.push_back(0);
      sampling.push_back(0);
    }

    ComputeGFA(Egs, VDecomp_coefs, sampling, gfa, VDecomp_phases, xAEcoefs, dimension, entries_size,
        distinctAE_size, vcpt, VStoichE, species);

    ofstream GFAData;
    string GFAfilename;

    GFAfilename = "GFA_"+input+".dat";
    GFAData.open(GFAfilename.c_str());
    for(uint i=0;i<dimension;i++){
      GFAData << species[i] << "  ";
    }
    GFAData << "GFA" << endl;
    for(uint i=0;i<vcpt.size();i++){
      double last=1;
      for(int j=vcpt[i].m_coords.lrows;j<=vcpt[i].m_coords.urows-1;j++){  //CO20190424
        GFAData << vcpt[i].m_coords[j] << "  ";
        last=last-vcpt[i].m_coords[j];
      }
      GFAData << last << "  " << gfa[i] << endl;
    }

    GFAData.close();

    cout << endl << "DONE!" << endl;

  }
}
