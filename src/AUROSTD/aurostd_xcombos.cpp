// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *           Aflow COREY OSES - Duke University 2013-2021                  *
// *                                                                         *
// ***************************************************************************

#ifndef _AUROSTD_XCOMBOS_CPP_
#define _AUROSTD_XCOMBOS_CPP_

#include "aurostd_xcombos.h"

#ifndef _AUROSTD_XERROR_H_
#include "aurostd_xerror.h"
#endif

#define _DEBUG_XCOMBOS_ false

using aurostd::xerror;

namespace aurostd {
  //--------------------------------------------------------------------------------
  // class xcombos
  //
  // A class to create permutations, combinations, and enumerations. xcombos uses
  // different "modes" to dinstinguish between these different types. A mode can
  // be called using the char 'P', 'C', or 'E'.
  //
  // Permutations (mode 'P') are called using a vector that contains the elements
  // to permute.
  //
  // Combinations (mode 'C') of m_choice out of n_choose elements are called using
  // two integers. They may be created with or without repetitions.
  // 
  // Enumerations (mode 'E') are combinations with repetitions that also include all
  // permutations. For example, the combinations with repetitions of the numbers
  // 0 and 1 of length 3 are {0, 0, 0}, {0, 0, 1}, {0, 1, 1}, and {1, 1, 1}.
  // The enumerations with the same parameters are {0, 0, 0}, {0, 0, 1}, {0, 1, 0},
  // {0, 1, 1}, {1, 0, 0}, {1, 0, 1}, {1, 1, 0}, {1, 1, 1}.
  // It is not necessary for all dimension in an enumeration to have the same size.
  // For example, the first dimension may just contain the number 0, which results
  // in the enumerations {0, 0, 0}, {0, 0, 1}, {0, 1, 0}, {0, 1, 1}.
  // If all dimensions are the same size, they are called like combinations but with
  // mode 'E'. If the dimensions have different sizes, they are called with a vector.
  // 
  //--------------------------------------------------------------------------------

  //--------------------------------------------------------------------------------
  // Permutation algorithms
  // Depending on which algorithm is used, the swapping order for permutations can
  // be different, see http://combos.org/perm
  //
  // Two algorithms are currently supported:
  //  1) Shen - lexicographical order (default)
  //            reference: doi:10.1007/BF01940170
  //            first swap: 1234 -> 1243
  //            example use case(s): POSCARs and POCC algorithm
  //  2) Heap - swap left-most first (fastest, minimum number of swaps)
  //            reference: https://en.wikipedia.org/wiki/Heap%27s_algorithm
  //            first swap: 1234 -> 2134
  //            example use case(s): finding representative atom decorations/permutations for prototypes
  //
  // If new algorithms are added, update the algorithms_xcombos enum.

  //------------------------------------------------------------------------------
  // constructor
  //------------------------------------------------------------------------------
  xcombos::xcombos() {free();}
  xcombos::xcombos(const std::vector<int>& vec, bool sort, char mode, algorithm_xcombos algorithm) {reset(vec, sort, mode, algorithm);} //DX20201222 - added algorithm
  xcombos::xcombos(int choice_count,int choose_count, char mode, bool rpt) {reset(choice_count, choose_count, mode, rpt);}
  xcombos::xcombos(const std::vector<int>& vec, char mode) {reset(vec, mode);}
  xcombos::xcombos(const xcombos& b) {copy(b);} // copy PUBLIC
  xcombos::~xcombos() {free();}

  void xcombos::free() {
    m_initialized=FALSE;
    m_input.clear();
    n_choices=0;
    m_choose=0;
    m_mode = '\0'; //ME20180529
    m_algorithm = shen_alg_xcombos; //DX20201222 - use Shen algorithm by default
    m_sets.clear(); //ME20180529
    m_sort=FALSE;
    m_started=FALSE;
    m_exhausted=FALSE;
    m_current.clear();
    m_indices.clear();
    m_p.clear();
    m_x=0;
    m_y=0;
    m_repeat = false; //ME20180509
  }

  void xcombos::copy(const xcombos& b) {  //copy PRIVATE
    m_initialized=b.m_initialized;
    m_input=b.m_input;
    n_choices=b.n_choices;
    m_choose=b.m_choose;
    m_indices = b.m_indices; //ME20180620
    m_mode = b.m_mode; //ME20180529
    m_algorithm = b.m_algorithm; //DX20201222
    m_sets = b.m_sets;
    m_sort=b.m_sort;
    m_started=b.m_started;
    m_exhausted=b.m_exhausted;
    m_current=b.m_current;
    m_p=b.m_p;
    m_x=b.m_x;
    m_y=b.m_y;
    m_repeat = b.m_repeat; //ME20180509
  }

  const xcombos& xcombos::operator=(const xcombos& other) {
    if(this!=&other) {free();copy(other);}
    return *this;
  }

  xcombos& xcombos::operator++() {  //remember, this is PREFIX (faster than POSTFIX)
    if(!m_initialized) {
      throw xerror(__AFLOW_FILE__,"xcombos::operator++()", "Cannot increment uninitialized xcombos class.",  _RUNTIME_INIT_);
    }
    if(m_mode=='P') {incrementPermutation();}
    else {
      if ((m_mode == 'C') && (!m_repeat)) {
        incrementCombinations();
      }
      else {
        incrementEnumerations();
      }
    }
    return *this;
  }

  void xcombos::reset() {  //clear PUBLIC
    if(m_mode == 'P') {
      reset(m_input, m_sort);
    } else {
      if (m_sets.size() > 0) {
        reset(m_sets, m_mode);
      } else {
        reset(n_choices, m_choose, m_mode, m_repeat);
      }
    }
  }

  void xcombos::reset(std::vector<int> vec,bool sort, char mode, algorithm_xcombos algorithm) { //do NOT make input vec a const &, this will screw up reset()
    free();
    m_input=vec;
    m_mode=mode;
    m_algorithm=algorithm; //DX20201222
    if (m_mode == 'p') {m_mode = 'P';}
    if (m_mode != 'P') {
      m_exhausted = TRUE;
      std::cerr << "xcombos::reset: Unrecognized mode " << m_mode << " for permutations." << std::endl;
      return;
    }
    m_sort=sort;
    m_started=FALSE;
    if(!m_input.size()){m_exhausted=TRUE; return;}  //safe
    initialize();
  }

  //ME20180509 - Implemented combinations with repetitions
  void xcombos::reset(int choice_count,int choose_count, char mode, bool rpt) {
    free();
    m_mode = mode;
    if (choice_count<=0 || choose_count<=0 || ((!rpt) && (choice_count<choose_count))) {m_exhausted=TRUE; return;}  //safe
    if (m_mode == 'e') {m_mode = 'E';}
    if (m_mode == 'E') {
      std::vector<int> vec(choose_count, choice_count);
      reset(vec, m_mode);
    } else {
      n_choices=choice_count;
      m_choose=choose_count;
      //[OBSOLETE ME20180705]m_permutation=FALSE;
      m_started=FALSE;
      m_repeat = rpt;
      if (m_mode == 'c') {m_mode = 'C';}
      if (m_mode != 'C') {
        m_exhausted = TRUE;
        std::cerr << "Unrecognized mode " << m_mode << " for combinations." << std::endl;
        return;
      }
      initialize();
    }
  }

  //ME20180529 - Enumerations
  void xcombos::reset(std::vector<int> vec, char mode) {
    free();
    m_mode = mode;
    if (m_mode == 'e') {m_mode = 'E';}
    if (m_mode != 'E') {
      m_exhausted = TRUE;
      std::cerr << "Unrecognized mode " << m_mode << " for enumerations." << std::endl;
      return;
    }
    m_choose = (int) vec.size();
    n_choices = vec[vec.size() - 1];
    m_sets = vec;
    bool all_zero = true;
    for (uint i = 0; i < m_sets.size(); i++) {
      if (m_sets[i] < 0) {
        std::cerr << "xcombos::reset: m_sets cannot contain negative numbers." << std::endl;
        m_exhausted = TRUE;
        return;
      } else if (m_sets[i] > 0) {
        all_zero = false;
      }
    }
    if (all_zero) {
      m_exhausted = TRUE;
      return;
    }
    initialize();
  }

  //ME20180509 - Implemented combinations with repetitions
  void xcombos::initialize() {
    bool LDEBUG=(FALSE || XHOST.DEBUG||_DEBUG_XCOMBOS_);
    if(m_mode == 'P') {
      m_current=m_input;
      if(m_sort) {std::sort(m_current.begin(),m_current.end());}
      //ME20180620
      m_indices.resize(m_input.size());
      for (uint i = 0; i < m_indices.size(); i++) {
        m_indices[i] = i;
      }
      // for Heap's algorithm //DX20201222
      m_p.clear();
      m_p.resize(m_input.size(),0);
      n_choices = (int)m_input.size(); //DX+ME20210111
      m_choose = (int)m_input.size(); //DX+ME20210111
    } else {
      if(LDEBUG) {
        if (m_mode == 'C') {
          std::cerr << __AFLOW_FUNC__ << " combinations: " << n_choices << " choose " << m_choose << ((m_repeat)? " with " : " without ") << "repetitions."; //CO20200404
        } else {
          std::cerr << __AFLOW_FUNC__ << " enumerations of length " << m_choose << std::endl;  //CO20200404
        }
      }
      if ((m_repeat) || (m_mode == 'E')) {
        m_current.resize(m_choose, 0);
      } else {
        m_current.resize(n_choices);
        for (uint i = 0; i < m_current.size(); i++) {
          if((int)i<m_choose) {
            m_current[i] = 1;
          } else {
            m_current[i] = 0;
          }
        }
        initializeCombinationsP();
      }
    }
    m_initialized=TRUE;
  }

  const std::vector<int>& xcombos::getCombo() const {return m_current;} //ME20190703 - use const & (faster)
  int xcombos::getN() const {return n_choices;}
  int xcombos::getM() const {return m_choose;}

  std::vector<int> xcombos::getIndices() const {
    std::vector<int> vout;
    if(m_mode == 'C'){ //only works for combinations, return nothing otherwise
      for(int i=0;i<(int)m_current.size();i++) {
        if(m_current[i]==1) {vout.push_back(i);}
      }
    }
    return vout;
  }

  template<class utype> std::vector<utype> xcombos::applyCombo(const std::vector<utype>& v_items) const {
    std::vector<utype> v_items_new;
    if(!(m_mode=='P' || m_mode=='C')){return v_items_new;}
    //if permutations, then m_current contains indices
    //otherwise, we need to get them
    //ME says combinations with repetitions is same structure as permutations
    std::vector<int> v_indices=m_current;
    if((m_mode=='C') && !(m_repeat)){v_indices.clear();v_indices=getIndices();} // combo indices
    for(uint i=0;i<v_indices.size();i++){
      if(v_indices[i]>=(int)v_items.size()){throw xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index",_INDEX_MISMATCH_);}
      v_items_new.push_back(v_items[v_indices[i]]);
    }
    return v_items_new;
  }

  template<class utype> std::vector<utype> xcombos::applyCombo(const std::vector<std::vector<utype> >& v_items) const {
    std::vector<utype> v_items_new;
    if(m_mode!='E'){return v_items_new;} //only applies to enumerations
    std::vector<int> v_indices=m_current;
    for(uint i=0;i<v_indices.size();i++){
      if(v_indices[i]>=(int)v_items[i].size()){throw xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index",_INDEX_MISMATCH_);}
      v_items_new.push_back(v_items[i][v_indices[i]]);
    }
    return v_items_new;
  }

  bool xcombos::increment() {
    if(!m_initialized) {return FALSE;}     //safety
    if(m_exhausted==TRUE) {return FALSE;}  //safety
    if(!m_started) {m_started=TRUE;return TRUE;} //allows us to use increment() for first grab too!
    ++(*this);
    if(m_exhausted==TRUE) {return FALSE;}
    return TRUE;
  }

  void xcombos::incrementPermutation() {
    bool LDEBUG=(FALSE || XHOST.DEBUG||_DEBUG_XCOMBOS_);
    if(m_exhausted) {return;}
    if(m_algorithm==shen_alg_xcombos){
      //Shen, MK. BIT (1962) 2(228). doi:10.1007/BF01940170
      //note this will generate next permutation in lexicographical order (left to right), e.g., first swap 1234 -> 1243
      int _i=-1;
      int _j=-1;
      for(int i=1;i<m_choose;i++) {if(m_current[i-1]<m_current[i]&&(i>_i)){_i=i;}}
      if(_i==-1) {m_exhausted=TRUE; return;} //stop condition
      for(int j=0;j<m_choose;j++) {if(m_current[_i-1]<m_current[j]&&(j>_j)){_j=j;}}
      if(LDEBUG) {cerr << "xcombos::incrementPermutation(): i=" << _i << "  j=" << _j << " " << endl;}
      std::swap(m_current[_i-1],m_current[_j]);
      std::swap(m_indices[_i-1], m_indices[_j]); //ME20180620
      for(int i=0;i<(m_choose-_i+1)/2;i++) {std::swap(m_current[_i+i],m_current[m_current.size()-i-1]);}
    }
    else if (m_algorithm==heap_alg_xcombos){ //DX20201222 - added this algorithm
      //Heap, B.R. (1963): https://en.wikipedia.org/wiki/Heap%27s_algorithm
      //note this will generate permutation by swapping lowest position index (left-most position), e.g., first swap 1234 -> 2134
      bool found_permutation=false;
      uint count=0, safety = m_choose; //DX20210111 - smarter saftey; each increment should swap no more than m_choose times
      while(!found_permutation&&count<=safety){
        count++;
        if(m_x>=m_choose){ m_exhausted=TRUE; return; } //stop condition
        if(m_p[m_x] < m_x){
          if(m_x%2==0){
            std::swap(m_current[0],m_current[m_x]);
          }
          else{
            std::swap(m_current[m_p[m_x]],m_current[m_x]);
          }
          m_p[m_x]+=1;
          m_x=0;
          found_permutation = true;
          if(LDEBUG) {cerr << "xcombos::incrementPermutation(): [HEAP] m_current: " << aurostd::joinWDelimiter(m_current, ",") << endl;}
        }
        else{
          m_p[m_x]=0;
          m_x+=1;
        }
      }
      if(count>safety){
        throw xerror(__AFLOW_FILE__,"xcombos::incrementPermutation()", "[HEAP] Uncontrolled while loop. Algorithm is not working as expected.",  _RUNTIME_ERROR_);
      }
    }
    else{ //DX2020111
      throw xerror(__AFLOW_FILE__,"xcombos::incrementPermutation()", "Invalide algorithm type, only Shen (shen_alg_xcombos) and Heap's (heap_alg_xcombos) algorithms are available.",  _INPUT_ERROR_);
    }
  }

  //ME20180509 - Implemented combinations with repetitions
  void xcombos::incrementCombinations() {
    if(m_exhausted) {return;}
    //[OBSOLETE ME20180705]if (repeat) { getNextRepeatCombination(); }
    //[OBSOLETE ME20180705]else { //[CO20200106 - close bracket for indenting]}
    setCombinationsIncrementParameters();
    m_current[m_x]=1;
    m_current[m_y]=0;
  }

  void xcombos::incrementEnumerations() {
    if (m_exhausted) {return;}
    getNextEnumeration();
  }

  void xcombos::initializeCombinationsP() {  //combinations only
    bool LDEBUG=(FALSE || XHOST.DEBUG||_DEBUG_XCOMBOS_);
    int i=0;
    m_p.resize(n_choices+2);
    m_p[i++]=-2;
    while(m_choose-i+1>0) {m_p[i]=m_choose-i+1;i++;}
    while(i<n_choices+1) {m_p[i++]=0;}
    m_p[n_choices+1]=n_choices+1;
    if(m_choose==0) {m_p[n_choices]=1;}
    if(LDEBUG) {
      cerr << "xcombos::initializeCombinationsP(): p(0): ";
      for(uint l=0;l<m_p.size();l++) {cerr << m_p[l] << " ";} 
      cerr << endl;
    }
  }

  void xcombos::setCombinationsIncrementParameters() { //combinations only
    if(!m_initialized) {return;}
    if(m_exhausted) {return;}
    bool LDEBUG=(FALSE || XHOST.DEBUG||_DEBUG_XCOMBOS_);
    //Chase, PJ. ACM (1970) 13(6). doi:10.1145/362384.362502
    //implementation modified from twiddle.c by Matthew Belmonte
    //http://www.netlib.no/netlib/toms/382
    //note this will generate next permutation in lexicographical order (left to right)
    //originally, Belmonte programmed right to left
    int i, j, k;
    i=j=k=-1;
    j=n_choices;
    if(LDEBUG) {cerr << "xcombos::setCombinationsIncrementParameters(): j(0)=" << j << endl;}
    while(m_p[j] <= 0) {j--;}
    if(LDEBUG) {cerr << "xcombos::setCombinationsIncrementParameters(): j(1)=" << j << endl;}
    if(m_p[j+1]==0) {
      for(i=j+1; i != n_choices; i++) {m_p[i]=-1;}
      m_p[j]=0;
      m_x=n_choices-1;
      m_p[n_choices]=1;
      m_y=j-1;
      if(LDEBUG) {
        cerr << "xcombos::setCombinationsIncrementParameters(): p(0): ";
        for(uint l=0;l<m_p.size();l++) {cerr << m_p[l] << " ";} 
        cerr << endl;
        cerr << "xcombos::setCombinationsIncrementParameters(): x=" << m_x << "(1)" << endl;
        cerr << "xcombos::setCombinationsIncrementParameters(): y=" << m_y << "(0)" << endl;
      }
    }else {
      if(j<n_choices) {m_p[j+1]=0;}
      j--;
      while(m_p[j]>0) {j--;}
      if(LDEBUG) {cerr << "xcombos::setCombinationsIncrementParameters(): j(2)=" << j << endl;}
      k=j+1;
      i=j;
      if(LDEBUG) {
        cerr << "xcombos::setCombinationsIncrementParameters(): i(0)=" << i << endl;
        cerr << "xcombos::setCombinationsIncrementParameters(): k(0)=" << k << endl;
      }
      while(m_p[i]==0) {m_p[i--]=-1;}
      if(LDEBUG) {cerr << "xcombos::setCombinationsIncrementParameters(): i(1)=" << i << endl;}
      if(m_p[i]==-1) {
        m_p[i]=m_p[k];
        m_x=i-1;
        m_y=k-1;
        m_p[k]=-1;
        if(LDEBUG) {
          cerr << "xcombos::setCombinationsIncrementParameters(): p(1): ";
          for(uint l=0;l<m_p.size();l++) {cerr << m_p[l] << " ";} 
          cerr << endl;
          cerr << "x=" << m_x << "(1)" << endl;
          cerr << "y=" << m_y << "(0)" << endl;
        }
      }else {
        if(LDEBUG) {
          cerr << "xcombos::setCombinationsIncrementParameters(): i=" << i << endl;
          cerr << "xcombos::setCombinationsIncrementParameters(): p(2): ";
          for(uint l=0;l<m_p.size();l++) {cerr << m_p[l] << " ";} 
          cerr << endl;
        }
        if(i==0) {m_exhausted=TRUE; return;}
        //if(i==m_p[n_choices+1]) {m_exhausted=TRUE; cerr << "HERE" << endl;return;}
        else {
          m_p[j]=m_p[i];
          m_p[i]=0;
          m_x=j-1;
          m_y=i-1;
          if(LDEBUG) {
            cerr << "xcombos::setCombinationsIncrementParameters(): p(3): ";
            for(uint l=0;l<m_p.size();l++) {cerr << m_p[l] << " ";} 
            cerr << endl;
            cerr << "x=" << m_x << "(1)" << endl;
            cerr << "y=" << m_y << "(0)" << endl;
          }
        }
      }
    }
    return;
  }

  //ME20180529
  // Creates the next enumeration in lexicographical order
  void xcombos::getNextEnumeration() {
    bool LDEBUG=(FALSE || XHOST.DEBUG || _DEBUG_XCOMBOS_);
    if (m_current[m_choose-1] == n_choices - 1) {
      int p = m_choose - 2;
      if (m_mode == 'E') {
        while ((p >= 0) && (m_current[p] == m_sets[p] - 1)) {
          p--;
        }
      } else {
        while ((p >= 0) && (m_current[p] == n_choices - 1)) {
          p--;
        }
      }
      if (LDEBUG) {cerr << "xcombos::getNextEnumeration: p(0) = " << p << endl;}
      if (p == -1)  {
        m_exhausted = true;
        return;
      } else {
        m_current[p]++;
        for (int i = p + 1; i < m_choose; i++) {
          if (m_mode == 'E') {
            m_current[i] = 0;
          } else {
            m_current[i] = m_current[p];
          }
        }
      }
    } else {
      m_current[m_choose-1]++;
    }
    if (LDEBUG) {
      cerr << "xcombos::getNextEnumeration: m_current(0): ";
      for (uint i = 0; i < m_current.size(); i++) {
        cerr << m_current[i] << " ";
      }
      cerr << endl;
    }
    return; 
  }
} // namespace aurostd

#endif  // _AUROSTD_XCOMBOS_CPP_
// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *           Aflow COREY OSES - Duke University 2013-2021                  *
// *                                                                         *
// ***************************************************************************
