// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
// Written by Nicholas Anderson
// nicholas.anderson@duke.edu

#ifndef AFLOW_NHULL_STATS_H
#define AFLOW_NHULL_STATS_H

#include <map>
#include <string>
#include <utility>
#include <vector>

#include "AUROSTD/aurostd_automatic_template.h"
#include "AUROSTD/aurostd_xserialization.h"
#include "AUROSTD/aurostd_xvector.h"

#include "modules/HULL/aflow_nhull_entry.h"
#include "modules/HULL/aflow_nhull_stats.h"

namespace nhull {
  //***************** STATISTICS FUNCTIONS
  double degenSize(std::vector<std::pair<double, int>>& v);
  double weightedMean(std::vector<std::pair<double, int>>& v);
  double degenVariance(std::vector<std::pair<double, int>>& v);

  // data structure for statistics of 1 compound
  // includes EFA and DEED
  struct compoundStatStore : public JsonSerializable<compoundStatStore> {
    std::map<std::string, int> compound_map;
    std::vector<Entry> members;
    std::vector<std::pair<double, int>> enthalpy_degeneracy_store;
    std::vector<std::pair<double, int>> delta_hull_degeneracy_store;
    aurostd::xvector<double> enthalpies;
    aurostd::xvector<double> delta_hulls;
    int total_members = 0;
    bool compute_statistics = true;
    double enthalpy_mean = 0;
    double enthalpy_median = 0;
    double enthalpy_variance = 0;
    double enthalpy_stdDev = 0;
    double delta_hull_mean = 0;
    double delta_hull_median = 0;
    double delta_hull_variance = 0;
    double delta_hull_stdDev = 0;
    double EFA = 0;
    double DEED = 0;

    compoundStatStore(const std::vector<Entry>& members);

    compoundStatStore() = default;

    // function used to compile statistics
    void compileStats();

    void poccStats();

    // SERIALIZATION MEMBERS Excluding: members, enthalpies, delta_hulls
#define JSON_compoundStatStore_MEMBERS \
  compound_map, total_members, compute_statistics, EFA, DEED, enthalpy_mean, enthalpy_median, enthalpy_variance, enthalpy_stdDev, delta_hull_mean, delta_hull_median, delta_hull_variance, delta_hull_stdDev

  protected:
    [[nodiscard]] std::string getJsonID() const override {
      return "compoundStatStore";
    };

  public:
    ~compoundStatStore() override = default;

    // JSON load/dump
    inline static const std::string JSONID = "compoundStatStore";

    [[nodiscard]] aurostd::JSON::object serialize() const override {
      return aurostd::JSON::object({AST_JSON_GETTER(JSON_compoundStatStore_MEMBERS)});
    }

    compoundStatStore deserialize(const aurostd::JSON::object& jo) override {
      AST_JSON_SETTER(JSON_compoundStatStore_MEMBERS)
      return *this;
    }
  };

  void singlePoccStats(std::vector<Entry>& entries, std::vector<std::pair<std::string, unsigned long long>> pocc_store);

  void statCalcLoop(std::vector<Entry>& entries, std::vector<compoundStatStore>& compounds);

  void compoundStats(std::vector<Entry> entries, const std::string& root_path);
} // namespace nhull

#endif // AFLOW_NHULL_STATS_H
