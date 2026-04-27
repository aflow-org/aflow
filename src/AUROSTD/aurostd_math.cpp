// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2025           *
// *                                                                         *
// ***************************************************************************

#include <algorithm>
#include <vector>

namespace aurostd {
    /// @brief Takes a vector of doubles that sums to one and updates small values to zero while maintaining the sum of one
    /// @param v
    /// @param cutoff_ratio optional cutoff for "small" ratios to be removed
    /// @return vector<double>
    /// @authors
    /// @mod {NHA,20251011,created}
  std::vector<double> absorbResidualsToMax(const std::vector<double>& v, const double cutoff_ratio = 1e-8) {
    std::vector<double> new_v;
    double residual_sum = 0; //storage var to sum all small residuals and later add to max value

    for (auto element : v) {
      if (element < cutoff_ratio) {
        residual_sum += element;
        new_v.push_back(0); //set small values to zero
      } else {
        new_v.push_back(element);
      }
    }

      // add residuals to the largest element
    auto maxElement = std::max_element(new_v.begin(), new_v.end());
    if (maxElement != new_v.end()) { // guard against empty vector
      *maxElement += residual_sum;
    }

    return new_v;
  }

}// namespace aurostd

// **************************************************************************
// *                                                                        *
// *             STEFANO CURTAROLO - Duke University 2003-2025              *
// *                                                                        *
// **************************************************************************
