// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************

#ifndef _AUROSTD_XCOMBOS_H_
#define _AUROSTD_XCOMBOS_H_

#include <vector>

enum algorithm_xcombos { // DX20210111
  shen_alg_xcombos,                               // lexicographical order (left to right); Shen, MK. BIT (1962) 2(228). doi:10.1007/BF01940170
  heap_alg_xcombos                                // swap left-most positions first (fast); Heap, B.R. (1963): https://en.wikipedia.org/wiki/Heap%27s_algorithm
};

namespace aurostd {
  class xcombos {
    // Written by Corey Oses corey.oses@duke.edu
  private:
      // NECESSARY PRIVATE CLASS METHODS - START
    void free();
    void copy(const xcombos& b);
      // NECESSARY PRIVATE CLASS METHODS - END

    bool m_initialized;
    std::vector<int> m_input;
    int n_choices, m_choose;
    bool m_sort;

    bool m_started;
    bool m_exhausted; // all possibilities explored
    char m_mode; // C: combinations, E: enumeration, P: permutations,
    algorithm_xcombos m_algorithm; // allow use of different algorithms (Shen, Heap, etc.) //DX20201222

      // for permutations
    std::vector<int> m_current;

      // for combinations
    std::vector<int> m_p;
    int m_x, m_y;
    bool m_repeat;  // Set to true if repetitions are allowed

      // for enumerations
    std::vector<int> m_sets;

    void incrementPermutation();
    void incrementCombinations();
    void incrementEnumerations();
    void initializeCombinationsP();
    void setCombinationsIncrementParameters();
    void getNextEnumeration(); // ME20180529
    void getNextEnumerationEqual();
    void initialize();

  public:
      // NECESSARY PUBLIC CLASS METHODS - START
      // constructors - START
    xcombos();
    xcombos(const std::vector<int>& vec, bool sort = true, char mode = 'P', algorithm_xcombos algorithm = shen_alg_xcombos);      // permutations //DX20201222 - added algorithm
    xcombos(int choice_count, int choose_count, char mode = 'C', bool rpt = false); // combinations, m choose n, m=choice_count, n=choose_count
    xcombos(const std::vector<int>& vec, char mode); // enumerations
    xcombos(const xcombos& b);
      // constructors - END
    ~xcombos();
    std::vector<int> m_indices;  // Keeps track of in which position original indices are in current permutation
    const xcombos& operator=(const xcombos& other);
    xcombos& operator++();
    void reset(); // reset with same parameters
    void reset(std::vector<int> vec, bool sort = true, char mode = 'P', algorithm_xcombos algorithm = shen_alg_xcombos); // reset with permutations //DX20201222 - added algorithm
    void reset(int choice_count, int choose_count, char mode = 'C', bool rpt = false); // reset with combinations
    void reset(std::vector<int> vec, char mode); // reset with enumerations
    [[nodiscard]] const std::vector<int>& getCombo() const; // grab current possibility - ME20190703 use const & (faster)
    [[nodiscard]] int getN() const; // grab n (total count)
    [[nodiscard]] int getM() const; // grab m (choose)
    [[nodiscard]] std::vector<int> getIndices() const; // get which indicies are 1
    template <class utype> [[nodiscard]] std::vector<utype> applyCombo(const std::vector<utype>& v_items) const;
    template <class utype> [[nodiscard]] std::vector<utype> applyCombo(const std::vector<std::vector<utype>>& v_items) const;
    bool increment(); // nice wrapper for ++, returns m_exhausted
      // NECESSARY PUBLIC CLASS METHODS - STOP
  };
} // namespace aurostd

#endif  // _AUROSTD_XCOMBOS_H_
// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
