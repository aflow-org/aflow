// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *           Aflow COREY OSES - Duke University 2013-2021                  *
// *                                                                         *
// ***************************************************************************
// Written by Corey Oses
// corey.oses@duke.edu

#ifndef _AUROSTD_XCOMBOS_H_
#define _AUROSTD_XCOMBOS_H_

enum algorithm_xcombos { //DX20210111
  shen_alg_xcombos,                               // lexicographical order (left to right); Shen, MK. BIT (1962) 2(228). doi:10.1007/BF01940170
  heap_alg_xcombos                                // swap left-most positions first (fast); Heap, B.R. (1963): https://en.wikipedia.org/wiki/Heap%27s_algorithm
};

namespace aurostd {
  class xcombos {
    private:
      //NECESSARY PRIVATE CLASS METHODS - START
      void free();
      void copy(const xcombos& b);
      //NECESSARY PRIVATE CLASS METHODS - END

      bool m_initialized;
      vector<int> m_input;
      int n_choices, m_choose;
      //[OBSOLETE ME20180705]bool m_permutation; //otherwise combination
      bool m_sort;

      bool m_started;
      bool m_exhausted; //all possibilities explored
      char m_mode; // C: combinations, E: enumeration, P: permutations,
      algorithm_xcombos m_algorithm; // allow use of different algorithms (Shen, Heap, etc.) //DX20201222

      //for permutations
      vector<int> m_current;

      //for combinations
      vector<int> m_p;
      int m_x,m_y;
      bool m_repeat;  // Set to true if repetitions are allowed

      //for enumerations
      std::vector<int> m_sets; 

      void incrementPermutation();
      void incrementCombinations();
      void incrementEnumerations();
      void initializeCombinationsP();
      void setCombinationsIncrementParameters();
      void getNextEnumeration(); //ME20180529
      void getNextEnumerationEqual();
      void initialize();
    public:
      //NECESSARY PUBLIC CLASS METHODS - START
      //constructors - START
      xcombos();
      xcombos(const vector<int>& vec, bool sort=TRUE, char mode='P', algorithm_xcombos algorithm=shen_alg_xcombos);      //permutations //DX20201222 - added algorithm
      xcombos(int choice_count,int choose_count, char mode='C', bool rpt=FALSE); //combinations, m choose n, m=choice_count, n=choose_count
      xcombos(const vector<int>& vec, char mode); // enumerations
      xcombos(const xcombos& b);
      //constructors - END
      ~xcombos();
      vector<int> m_indices;  // Keeps track of in which position original indices are in current permutation
      const xcombos& operator=(const xcombos& other) ;
      xcombos& operator++();
      void reset(); //reset with same parameters
      void reset(vector<int> vec,bool sort=TRUE, char mode='P', algorithm_xcombos algorithm=shen_alg_xcombos); //reset with permutations //DX20201222 - added algorithm
      void reset(int choice_count,int choose_count, char mode='C', bool rpt=FALSE); //reset with combinations
      void reset(vector<int> vec, char mode); // reset with enumerations
      const vector<int>& getCombo() const; //grab current possibility - ME20190703 use const & (faster)
      int getN() const; //grab n (total count)
      int getM() const; //grab m (choose)
      vector<int> getIndices() const; //get which indicies are 1
      template<class utype> std::vector<utype> applyCombo(const std::vector<utype>& v_items) const;
      template<class utype> std::vector<utype> applyCombo(const std::vector<std::vector<utype> >& v_items) const;
      bool increment(); //nice wrapper for ++, returns m_exhausted
      //NECESSARY PUBLIC CLASS METHODS - STOP
  };
} // namespace aurostd

#endif  // _AUROSTD_XCOMBOS_H_
// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *           Aflow COREY OSES - Duke University 2013-2021                  *
// *                                                                         *
// ***************************************************************************
